program fitpot
  use variables
  use parallel
  implicit none

  interface
    subroutine write_vars(cadd)
      character(len=*),intent(in),optional:: cadd
    end subroutine write_vars
    subroutine write_energy_relation(cadd)
      character(len=*),intent(in),optional:: cadd
    end subroutine write_energy_relation
    subroutine write_force_relation(cadd)
      character(len=*),intent(in),optional:: cadd
    end subroutine write_force_relation
  end interface
  
  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world,nnode,ierr)
  call mpi_comm_rank(mpi_comm_world,myid,ierr)
  mpi_world= mpi_comm_world

  if( myid.eq.0 ) then
    call read_input(10,'in.fitpot')
    call write_initial_setting()
  endif
  call sync_input()
  allocate(samples(nsmpl))
  if( nnode.gt.nsmpl ) then
    if( myid.eq.0 ) then
      print *,'[Error] nnode.gt.nsmpl'
      print *,'  nnode,nspml=',nnode,nsmpl
      print *,'you should use less number of node than nsmpl.'
    endif
    call mpi_finalize(ierr)
    stop
  endif

  if( myid.eq.0 ) then
    print *,' nsmpl,nnode=',nsmpl,nnode
  endif
  call get_node2sample()

  call get_dir_list(11)

  call read_samples()
  call read_ref_data()

  call scatter_samples()

  call read_vars()

  select case (trim(cfmethod))
    case ('bfgs','BFGS')
      call bfgs_wrapper()
    case ('check_grad')
      call check_grad()
    case ('test','TEST')
      call test()
    case default
      print *,'unknown fitting_method:',trim(cfmethod)
      stop
  end select

  if( myid.eq.0 ) then
    call write_vars('fin')
    call write_energy_relation('fin')
    call write_force_relation('fin')
    call write_statistics()
  endif
  call mpi_finalize(ierr)

  write(6,'(a,f15.3,a)') ' time function =', timef,' sec'
  write(6,'(a,f15.3,a)') ' time gradient =', timeg,' sec'

end program fitpot
!=======================================================================
subroutine write_initial_setting()
  use variables
  implicit none 
  integer:: i

  write(6,'(a)') '---------------- INITIAL SETTING -----------------'
  write(6,'(2x,a25,2x,i5)') 'num_samples',nsmpl
  write(6,'(2x,a25,2x,i5)') 'num_iteration',nstp
  write(6,'(2x,a25,2x,a)') 'fitting_method',trim(cfmethod)
  write(6,'(2x,a25,2x,a)') 'main_directory',trim(cmaindir)
  write(6,'(2x,a25,2x,a)') 'param_file',trim(cparfile)
  write(6,'(2x,a25,2x,a)') 'run_mode',trim(crunmode)
  write(6,'(2x,a25,2x,es12.3)') 'eps',eps
  write(6,'(2x,a25,2x,es12.3)') 'xtol',xtol
  do i=1,maxnsp
    write(6,'(2x,a25,2x,i2,es15.7)') 'atom_energy',i,eatom(i)
  enddo
  write(6,'(2x,a25,2x,l3)') 'force_match',lfmatch
  write(6,'(2x,a25,2x,l3)') 'penalty',lpena
  write(6,'(2x,a25,2x,es12.3)') 'penalty_weight',pwgt
  write(6,'(2x,a25,2x,a)') 'potential',trim(cpot)
  write(6,'(2x,a25,2x,l3)') 'gradient',lgrad
  write(6,'(2x,a25,2x,l3)') 'grad_scale',lgscale
  write(6,'(2x,a25,2x,es12.3)') 'gscale_factor',gscl
  write(6,'(2x,a25,2x,l3)') 'regularize',lreg
  write(6,'(2x,a25,2x,l3)') 'force_scale',lfscale
  write(6,'(2x,a25,2x,l3)') 'sample_weight',lswgt
  write(6,'(2x,a25,2x,es12.3)') 'sample_weight_beta',swbeta
  write(6,'(a)') '------------------------------------------------'

end subroutine write_initial_setting
!=======================================================================
subroutine get_dir_list(ionum)
  use variables
  use parallel
  implicit none
  integer,intent(in):: ionum
  character(len=128),allocatable:: cdname(:)
  integer:: is

  allocate(cdname(nsmple))
  
  if( myid.eq.0 ) then
    call system('ls '//trim(cmaindir) &
         //' | grep "^[0-9]...." > dir_list.txt')
    open(ionum,file='dir_list.txt',status='old')
    do is=1,nsmpl
      !read(ionum,*,end=999) samples(i)%cdirname
      read(ionum,*,end=999) cdname(is)
    enddo
    close(ionum)
  endif
  
  call mpi_bcast(cdname,128,mpi_character,0,mpi_world,ierr)

  do is=1,nsmpl
    samples(is)%cdirname= cdname(is)
  enddo

!!$  do i=1,nsmpl
!!$    print *,' i,cdirlist(i)=',i,cdirlist(i)
!!$  enddo
  deallocate(cdname)
  return

999 continue
  print *,' Error: num_samples may be wrong.'
  call mpi_finalize(ierr)
  stop

end subroutine get_dir_list
!=======================================================================
subroutine read_samples()
  use variables
  use parallel
  implicit none
  integer:: is
  character*5:: cdir

  if( myid.eq.0 ) then
    do is=1,nsmpl
!    print *,' is=',is,cdirlist(is)
      cdir= samples(is)%cdirname
      call read_pos(12,trim(cmaindir)//'/'//trim(cdir) &
           //'/pos',is)
    enddo
  endif

  !if( myid.eq.0 ) print *,'read_samples done.'

end subroutine read_samples
!=======================================================================
subroutine read_pos(ionum,fname,ismpl)
  use variables
  implicit none 
  integer,intent(in):: ionum,ismpl
  character(len=*),intent(in):: fname

  integer:: i,natm
  real(8):: tmp

  open(ionum,file=trim(fname),status='old')
  read(ionum,*) samples(ismpl)%h0
  read(ionum,*) samples(ismpl)%h(1,1:3)
  read(ionum,*) samples(ismpl)%h(2,1:3)
  read(ionum,*) samples(ismpl)%h(3,1:3)
  read(ionum,*) tmp,tmp,tmp
  read(ionum,*) tmp,tmp,tmp
  read(ionum,*) tmp,tmp,tmp
  read(ionum,*) natm
  samples(ismpl)%natm= natm
  allocate(samples(ismpl)%ra(3,natm),samples(ismpl)%fa(3,natm) &
       ,samples(ismpl)%tag(natm) &
       ,samples(ismpl)%fref(3,natm))
  do i=1,samples(ismpl)%natm
    read(ionum,*) samples(ismpl)%tag(i),samples(ismpl)%ra(1:3,i), &
         tmp,tmp,tmp
!    write(6,'(4es22.14)')  smpl%tag(i),smpl%ra(1:3,i)
  enddo
  close(ionum)
end subroutine read_pos
!=======================================================================
subroutine read_ref_data()
  use variables
  use parallel
  integer:: ismpl,i

  if( myid.eq.0 ) then
    do ismpl=1,nsmpl
      open(13,file=trim(cmaindir)//'/'//samples(ismpl)%cdirname &
           //'/erg.ref',status='old')
      read(13,*) samples(ismpl)%eref
      close(13)

      open(14,file=trim(cmaindir)//'/'//samples(ismpl)%cdirname &
           //'/frc.ref',status='old')
      read(14,*) natm
      if( natm.ne.samples(ismpl)%natm ) then
        print *,'Error: natm in sample is not same as smpl%natm'
        print *,' natm,smpl%natm=',natm,samples(ismpl)%natm
        stop
      endif
      do i=1,natm
        read(14,*) samples(ismpl)%fref(1:3,i)
      enddo
      close(14)
    enddo
  endif

!!$  print *,'read_ref_data done.'

end subroutine read_ref_data
!=======================================================================
subroutine read_vars()
  use variables
  use parallel
  implicit none

  integer:: i

  if( myid.eq.0 ) then
    open(15,file=trim(cmaindir)//'/'//cparfile,status='old')
    read(15,*) nvars, rcut
  endif
  call mpi_bcast(nvars,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(rcut,1,mpi_double_precision,0,mpi_world,ierr)
  allocate(vars(nvars),vranges(2,nvars))
  if( myid.eq.0 ) then
    do i=1,nvars
      read(15,*) vars(i),vranges(1:2,i)
!    print *, vars(i),vranges(1:2,i)
    enddo
  endif

  call mpi_bcast(vars,nvars,mpi_double_precision,0,mpi_world,ierr)

  close(15)
end subroutine read_vars
!=======================================================================
subroutine write_vars(cadd)
  use variables
  implicit none
  character(len=*),intent(in),optional:: cadd
  integer:: i
  character(len=128):: cfname

  cfname= trim(cmaindir)//'/'//trim(cparfile)
  if( present(cadd) ) cfname= trim(cfname)//'.'//trim(cadd)
  

  open(15,file=trim(cfname),status='replace')
  write(15,'(i10,es15.4)') nvars,rcut
  do i=1,nvars
    write(15,'(es23.14e3,2es12.4)') vars(i),vranges(1:2,i)
  enddo
  close(15)
!!$  print *, 'vars written.'

end subroutine write_vars
!=======================================================================
subroutine bfgs_wrapper()
  use variables
  use NN,only:NN_init,NN_get_f,NN_get_g
  implicit none
  integer:: i,m
  real(8):: fval,fmin
  real(8),allocatable:: gval(:),w(:),diag(:)

  !.....for lbfgs library
  integer,parameter:: msave= 7
  integer:: mp,lp,iprint(2),iflag,istp,nwork
  real(8):: bgtol,bstpmin,bstpmax
  logical:: ldiagco
!  external:: lb2
  common /lb3/mp,lp,bgtol,bstpmin,bstpmax
  
  nwork= nvars*(2*msave+1) +2*msave
  allocate(gval(nvars),w(nwork),diag(nvars))

  !.....NN specific code hereafter
  call NN_init()
  call NN_get_f(fval)
  call NN_get_g(gval)
!!$  print *, 'fval=',fval
!!$  print *, 'gval:'
!!$  do i=1,nvars
!!$    print *, 'i,g=',i,gval(i)
!!$  enddo

  !.....bfgs wrapper
  fmin= 1d+20
  m= 5
  ldiagco= .false.
  iprint(1)= 1
  iprint(2)= 0
  istp= 0
  iflag= 0
10 continue
  call NN_get_f(fval)
  call NN_get_g(gval)
  call lbfgs(nvars,m,vars,fval,gval,ldiagco,diag,iprint,eps,xtol &
       ,w,iflag)
!!$  print *, 'istp,fval,iflag=',istp,fval,iflag
!!$  if( fmin.gt.fval ) then
!!$    fmin= fval
!!$    call write_vars()
!!$  endif
  if( iflag.le.0 ) goto 999
  istp= istp +1
  if( istp.gt.nstp ) goto 999
  goto 10

999 continue
  return
end subroutine bfgs_wrapper
!=======================================================================
subroutine check_grad()
  use variables
  use NN,only:NN_init,NN_get_f,NN_get_g
  implicit none
  integer:: iv
  real(8):: f0,ftmp,dv,vmax
  real(8),allocatable:: ganal(:),gnumer(:),vars0(:)

  allocate(gnumer(nvars),ganal(nvars),vars0(nvars))
  call NN_init()
  call NN_get_f(f0)
  call NN_get_g(ganal)

  vars0(1:nvars)= vars(1:nvars)
  vmax= 0d0
  do iv=1,nvars
    vmax= max(vmax,abs(vars0(iv)))
    write(6,'(a,i6,es12.4)') ' iv,vars(iv)=',iv,vars0(iv)
  enddo
  dv= vmax *1d-5
  print *,''
  print *,'deviation [dv] =',dv
  do iv=1,nvars
    vars(1:nvars)= vars0(1:nvars)
!!$    dv= vars(iv) *1d-3
    vars(iv)= vars(iv) +dv
    call NN_get_f(ftmp)
    gnumer(iv)= (ftmp-f0)/dv
  enddo

  write(6,'(a)') '----------------- check_grad ------------------------'
  write(6,'(a)') '     #,    analytical,'// &
       '     numerical,'// &
       '      error [%]'
  do iv=1,nvars
    write(6,'(i6,2es15.4,f15.3)') iv,ganal(iv),gnumer(iv), &
         abs((ganal(iv)-gnumer(iv))/gnumer(iv))*100
  enddo
  write(6,'(a)') '-----------------------------------------------------'
  print *, 'check_grad done.'
end subroutine check_grad
!=======================================================================
subroutine test()
  use variables
  use NN,only:NN_init,NN_get_f,NN_get_g
  implicit none 
  integer:: iv
  real(8):: ft
  real(8),allocatable:: gt(:)

  allocate(gt(nvars))

  call NN_init()
  call NN_get_f(ft)
  call NN_get_g(gt)
  
  print *,'func value=',ft
  print *,'grad values:'
  do iv=1,nvars
    print *,'iv,grad(iv)=',iv,gt(iv)
  enddo
  print *,'test done.'
end subroutine test
!=======================================================================
subroutine write_energy_relation(cadd)
  use variables
  implicit none
  character(len=*),intent(in),optional:: cadd
  character(len=128):: cfname
  
  integer:: ismpl
  type(mdsys):: smpl
  
  cfname='out.erg.smd-vs-dft'
  if( present(cadd) ) cfname= trim(cfname)//'.'//trim(cadd)

  open(90,file=trim(cfname),status='replace')
  do ismpl=1,nsmpl
    smpl= samples(ismpl)
    write(90,'(2es15.7,2x,a)') smpl%eref/smpl%natm &
         ,smpl%epot/smpl%natm,smpl%cdirname
  enddo
  close(90)
  
end subroutine write_energy_relation
!=======================================================================
subroutine write_force_relation(cadd)
  use variables
  implicit none
  character(len=*),intent(in),optional:: cadd
  character(len=128):: cfname

  integer:: ismpl,ia,ixyz
  type(mdsys):: smpl
  
  cfname= 'out.frc.smd-vs-dft'
  if( present(cadd) ) cfname= trim(cfname)//'.'//trim(cadd)

  open(91,file=trim(cfname),status='replace')
  do ismpl=1,nsmpl
    smpl= samples(ismpl)
    do ia=1,smpl%natm
      do ixyz=1,3
        write(91,'(2es15.7,2x,a)') smpl%fref(ixyz,ia)/smpl%natm &
             ,smpl%fa(ixyz,ia)/smpl%natm,smpl%cdirname
      enddo
    enddo
  enddo
  close(91)
  
end subroutine write_force_relation
!=======================================================================
subroutine write_statistics()
  use variables
  implicit none
  integer:: ismpl,ia,l,n,natm
  real(8):: demax,desum,de,rmse,dfmax,dfsum,df
  
  write(6,'(/,a)') '>>>>> statistics:'
  !.....energies
  demax= 0d0
  desum= 0d0
  do ismpl=1,nsmpl
    natm= samples(ismpl)%natm
    de= abs(samples(ismpl)%epot -samples(ismpl)%eref)/natm
    demax= max(demax,de)
    desum=desum +de*de/nsmpl
  enddo
  rmse= sqrt(desum)
  write(6,'(a,f12.3,a)') '  RMSE of energies         =',rmse,' eV/atom'
  write(6,'(a,f12.3,a)') '  Max residual of energies =',demax,' eV/atom'

  !.....forces
  dfmax= 0d0
  dfsum= 0d0
  n= 0
  do ismpl=1,nsmpl
    natm= samples(ismpl)%natm
    do ia=1,natm
      do l=1,3
        df= abs(samples(ismpl)%fa(l,ia)-samples(ismpl)%fref(l,ia))
        dfmax= max(dfmax,df)
        dfsum=dfsum +df*df
        n=n +1
      enddo
    enddo
  enddo
  rmse= sqrt(dfsum/n)
  write(6,'(a,f12.3,a)') '  RMSE of forces           =',rmse,' eV/atom'
  write(6,'(a,f12.3,a)') '  Max residual of forces   =',dfmax,' eV/atom'
  print *,''
end subroutine write_statistics
!=======================================================================
subroutine sync_input()
  use variables
  use parallel
  implicit none
  
  call mpi_bcast(nsmpl,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(nstp,1,mpi_integer,0,mpi_world,ierr)

  call mpi_bcast(cfmethod,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(cmaindir,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(cparfile,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(crunmode,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(cpot,128,mpi_character,0,mpi_world,ierr)

  call mpi_bcast(eps,1,mpi_double_precision,0,mpi_world,ierr)
  call mpi_bcast(xtol,1,mpi_double_precision,0,mpi_world,ierr)
  call mpi_bcast(eatom,maxnsp,mpi_double_precision,0,mpi_world,ierr)
  call mpi_bcast(gscl,1,mpi_double_precision,0,mpi_world,ierr)
  call mpi_bcast(fscl,1,mpi_double_precision,0,mpi_world,ierr)
  call mpi_bcast(swbeta,1,mpi_double_precision,0,mpi_world,ierr)
  call mpi_bcast(pwgt,1,mpi_double_precision,0,mpi_world,ierr)
  
  call mpi_bcast(lfmatch,1,mpi_logical,0,mpi_world,ierr)
  call mpi_bcast(lreg,1,mpi_logical,0,mpi_world,ierr)
  call mpi_bcast(lgrad,1,mpi_logical,0,mpi_world,ierr)
  call mpi_bcast(lgscale,1,mpi_logical,0,mpi_world,ierr)
  call mpi_bcast(lfscale,1,mpi_logical,0,mpi_world,ierr)
  call mpi_bcast(lswgt,1,mpi_logical,0,mpi_world,ierr)
  call mpi_bcast(lpena,1,mpi_logical,0,mpi_world,ierr)
end subroutine sync_input
!=======================================================================
subroutine get_node2sample()
  use variables
  use parallel
  implicit none
  integer:: n,m,ip

  n= nsmpl/nnode
  m= nsmpl -n*nnode
  allocate(nspn(nnode))
  do ip=1,nnode
    nspn(ip)= n
    if( ip.le.m ) nspn(ip)= nspn(ip) +1
  enddo

  !.....compute start and end of sample-id of this process
  isid0= 0
  isid1= 0
  do ip=1,nnode
    isid1= isid1 +nspn(ip)
    if( ip.eq.myid+1 ) then
      isid0= isid1 -nspn(ip) +1
      exit
    endif
  enddo

  if( myid.eq.0 ) print *,'get_node2sample done.'
  return
end subroutine get_node2sample
!=======================================================================
subroutine scatter_samples()
  use variables
  use parallel
  implicit none

  do is=isid0,isid1
    
  enddo
end subroutine scatter_samples
