program fitpot
  use variables
  use parallel
  implicit none
  real(8):: tmp

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
  time0= mpi_wtime()
  call mpi_comm_size(mpi_comm_world,nnode,ierr)
  call mpi_comm_rank(mpi_comm_world,myid,ierr)
  mpi_world= mpi_comm_world
  tcomm= 0d0

  if( myid.eq.0 ) then
    call read_input(10,'in.fitpot')
    call write_initial_setting()
  endif
  call sync_input()

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

  allocate(samples(isid0:isid1))

  call get_dir_list(11)

  call read_samples()
  call read_ref_data()

!!$  call scatter_samples()

  call read_vars()

  select case (trim(cfmethod))
    case ('sd','SD')
      call sd_wrapper()
    case ('cg','CG')
      call cg_wrapper()
    case ('bfgs','BFGS')
      call bfgs_wrapper()
    case ('sequential')
      call sequential_update()
    case ('check_grad')
      call check_grad()
    case ('test','TEST')
      call test()
    case default
      if(myid.eq.0) print *,'unknown fitting_method:',trim(cfmethod)
      stop
  end select

  if( myid.eq.0 ) then
    if( iflag/1000.ne.0 ) then
      print *,'something wrong with 1D line search.'
      print *,'  iflag=',iflag
    else if( iflag/100.ne.0 ) then
      print *,'something wrong with minimization.'
      print *,'  iflag=',iflag
    endif
  endif
  call write_vars('fin')
  call write_energy_relation('fin')
  call write_force_relation('fin')
  call write_statistics()
!!$  write(6,'(a,i4,3f15.3)') ' myid,tfunc,tgrad,tcom=' &
!!$       ,myid,tfunc,tgrad,tcomm
  tmp= tfunc
  call mpi_reduce(tmp,tfunc,1,mpi_double_precision,mpi_max &
       ,0,mpi_world,ierr)
  tmp= tgrad
  call mpi_reduce(tmp,tgrad,1,mpi_double_precision,mpi_max &
       ,0,mpi_world,ierr)
  tmp= tcomm
  call mpi_reduce(tmp,tcomm,1,mpi_double_precision,mpi_max &
       ,0,mpi_world,ierr)
  if( myid.eq.0 ) then
    write(6,'(a,i8)') ' num of func calls=',nfunc
    write(6,'(a,i8)') ' num of grad calls=',ngrad
    write(6,'(a,f15.3,a)') ' time func =', tfunc,' sec'
    write(6,'(a,f15.3,a)') ' time grad =', tgrad,' sec'
    write(6,'(a,f15.3,a)') ' time comm =', tcomm,' sec'
    write(6,'(a,f15.3,a)') ' time      =', mpi_wtime() -time0, ' sec'
  endif
  call mpi_finalize(ierr)

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
  write(6,'(2x,a25,2x,es12.3)') 'xtol',xtol
  write(6,'(2x,a25,2x,es12.3)') 'ftol',ftol
  write(6,'(2x,a25,2x,es12.3)') 'gtol',gtol
  do i=1,maxnsp
    write(6,'(2x,a25,2x,i2,es15.7)') 'atom_energy',i,eatom(i)
  enddo
  write(6,'(2x,a25,2x,l3)') 'force_match',lfmatch
  write(6,'(2x,a25,2x,a)') 'penalty',trim(cpena)
  write(6,'(2x,a25,2x,es12.3)') 'penalty_weight',pwgt
  write(6,'(2x,a25,2x,a)') 'potential',trim(cpot)
  write(6,'(2x,a25,2x,l3)') 'gradient',lgrad
  write(6,'(2x,a25,2x,l3)') 'grad_scale',lgscale
  write(6,'(2x,a25,2x,es12.3)') 'gscale_factor',gscl
  write(6,'(2x,a25,2x,l3)') 'regularize',lreg
  write(6,'(2x,a25,2x,l3)') 'force_scale',lfscale
  write(6,'(2x,a25,2x,l3)') 'sample_weight',lswgt
  write(6,'(2x,a25,2x,es12.3)') 'sample_weight_beta',swbeta
  write(6,'(2x,a25,2x,es12.3)') 'coeff_sequential',seqcoef
  write(6,'(2x,a25,2x,a)') 'line_minimization',trim(clinmin)
  write(6,'(a)') '------------------------------------------------'

end subroutine write_initial_setting
!=======================================================================
subroutine get_dir_list(ionum)
  use variables
  use parallel
  implicit none
  integer,intent(in):: ionum
  integer:: is

  if( .not. allocated(cdirlist)) allocate(cdirlist(nsmpl))
  
  if( myid.eq.0 ) then
    call system('ls '//trim(cmaindir) &
         //' | grep "^[0-9]...." > dir_list.txt')
    open(ionum,file='dir_list.txt',status='old')
    do is=1,nsmpl
      read(ionum,*,end=999) cdirlist(is)
    enddo
    close(ionum)
  endif
  call mpi_barrier(mpi_world,ierr)
  call mpi_bcast(cdirlist,5*nsmpl,mpi_character,0,mpi_world,ierr)
  
  do is=isid0,isid1
    samples(is)%cdirname= cdirlist(is)
!!$    write(6,'(a,2i4,a)') ' myid,ismpl,cdirname=',myid,is,cdirlist(is)
  enddo

  if(myid.eq.0) print*,'get_dir_list done.'
  return

999 continue
  if( myid.eq.0 ) print *,' Error: num_samples may be wrong.'
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
  integer,allocatable,save:: nal(:)

  if( .not. allocated(nalist) ) allocate(nalist(nsmpl))
  allocate(nal(nsmpl))
  nalist(1:nsmpl)= 0
  nal(1:nsmpl)= 0

  do is=isid0,isid1
!    print *,' is=',is,cdirlist(is)
    cdir= samples(is)%cdirname
    call read_pos(12,trim(cmaindir)//'/'//trim(cdir) &
         //'/pos',is)
    nal(is)= samples(is)%natm
  enddo
  call mpi_reduce(nal,nalist,nsmpl,mpi_integer,mpi_sum &
       ,0,mpi_world,ierr)
!!$  if( myid.eq.0 ) then
!!$    do is=1,nsmpl
!!$      print *,' ismpl,natm=',is,nalist(is)
!!$    enddo
!!$  endif
  
  if( myid.eq.0 ) print *,'read_samples done.'
  call mpi_barrier(mpi_world,ierr)
  return
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
  implicit none 
  integer:: ismpl,i,is,jflag,natm

  jflag= 0
  do ismpl=isid0,isid1
    open(13,file=trim(cmaindir)//'/'//samples(ismpl)%cdirname &
         //'/erg.ref',status='old')
    read(13,*) samples(ismpl)%eref
    close(13)
!.....reduce atomic energy from eref
    do i=1,samples(ismpl)%natm
      is= samples(ismpl)%tag(i)
      samples(ismpl)%eref= samples(ismpl)%eref -eatom(is)
    enddo

    open(14,file=trim(cmaindir)//'/'//samples(ismpl)%cdirname &
         //'/frc.ref',status='old')
    read(14,*) natm
    if( natm.ne.samples(ismpl)%natm ) then
      print *,'Error: natm in sample is not same as smpl%natm'
      print *,' myid,ismpl,natm,smpl%natm=',myid,ismpl &
           ,natm,samples(ismpl)%natm
      jflag= jflag +1
    endif
    do i=1,natm
      read(14,*) samples(ismpl)%fref(1:3,i)
    enddo
    close(14)
  enddo

  if( jflag.gt.0 ) then
    call mpi_finalize(ierr)
    stop
  endif

  if(myid.eq.0) print *,'read_ref_data done.'

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
  use parallel
  implicit none
  character(len=*),intent(in),optional:: cadd
  integer:: i
  character(len=128):: cfname

  cfname= trim(cmaindir)//'/'//trim(cparfile)
  if( present(cadd) ) cfname= trim(cfname)//'.'//trim(cadd)

  if( myid.eq.0 ) then
    open(15,file=trim(cfname),status='replace')
    write(15,'(i10,es15.4)') nvars,rcut
    do i=1,nvars
      write(15,'(es23.14e3,2es12.4)') vars(i),vranges(1:2,i)
    enddo
    close(15)
    print *, 'vars written.'
  endif

end subroutine write_vars
!=======================================================================
subroutine bfgs_wrapper()
  use variables
  use NN,only:NN_init,NN_func,NN_grad
  use parallel
  implicit none
  integer:: i,m
  real(8):: fval

  !.....NN specific code hereafter
  call NN_init()
  call bfgs(nvars,vars,fval,xtol,gtol,ftol,nstp &
       ,iprint,iflag,myid,NN_func,NN_grad,cpena,clinmin)

  return
end subroutine bfgs_wrapper
!=======================================================================
subroutine sd_wrapper()
!
!  Steepest descent minimization
!
  use variables
  use NN,only:NN_init,NN_func,NN_grad
  use parallel
  implicit none
  integer:: i,m
  real(8):: fval

  !.....NN specific code hereafter
  call NN_init()
  call steepest_descent(nvars,vars,fval,xtol,gtol,ftol,nstp&
       ,iprint,iflag,myid,NN_func,NN_grad)

  return
end subroutine sd_wrapper
!=======================================================================
subroutine cg_wrapper()
  use variables
  use NN,only:NN_init,NN_func,NN_grad
  use parallel
  implicit none
  integer:: i,m
  real(8):: fval

  !.....NN specific code hereafter
  call NN_init()
  call bfgs(nvars,vars,fval,xtol,gtol,ftol,nstp &
       ,iprint,iflag,myid,NN_func,NN_grad)

  return
end subroutine cg_wrapper
!=======================================================================
subroutine sequential_update()
  use variables
  use NN,only:NN_init,NN_fs,NN_gs,NN_func,NN_grad
  use parallel
  implicit none
  integer,parameter:: nstp_eval= 1
  integer,parameter:: nstp_time= 1
  real(8),allocatable:: gval(:)
  integer:: istp,iv
  real(8):: gnorm,alpha,gmax,vmax,fval,gg
  integer:: ismpl
  common /samplei/ ismpl
  real(8),external:: sprod

  interface
    subroutine write_vars(cadd)
      character(len=*),intent(in),optional:: cadd
    end subroutine write_vars
  end interface

  allocate(gval(nvars))

  if( nnode.ne.1 ) then
    if(myid.eq.0) print *,'[Error] sequential mode only works '// &
         'with single process.'
    call mpi_finalize(ierr)
    stop
  endif

  call NN_init()

  do istp=1,nstp
    if(mod(istp,nstp_eval).eq.0) then
      fval= NN_func(nvars,vars)
      gval= NN_grad(nvars,vars)
      gnorm= 0d0
      do iv=1,nvars
        gnorm= gnorm +gval(iv)*gval(iv)
      enddo
      write(6,'(a,i6,2f20.7,f10.3)') ' istp,f,gnorm,time=',istp,fval &
           ,gnorm ,mpi_wtime()-time0
      call write_vars('tmp')
    else if(mod(istp,nstp_time).eq.0) then
      write(6,'(a,i6,f10.3)') ' istp,time=',istp,mpi_wtime()-time0
    endif
    do ismpl=1,nsmpl
      fval= NN_fs(nvars,vars)
      gval= NN_gs(nvars,vars)
      gnorm= sprod(nvars,gval,gval)
      gval(1:nvars)= -gval(1:nvars)/sqrt(gnorm)
      call quad_interpolate(nvars,vars,gval,fval,xtol,gtol,ftol,alpha &
           ,iprint,iflag,myid,NN_fs)
      if( iflag/100.ne.0 ) then
        iflag= iflag -(iflag/100)*100
        print*,'since quad_interpolate failed, call golden_section.'
        call golden_section(nvars,vars,gval,fval,xtol,gtol,ftol,alpha &
           ,iprint,iflag,myid,NN_fs)
      endif
      vars(1:nvars)=vars(1:nvars) +alpha*gval(1:nvars)
!!$      call NN_get_f(fval)
!!$      call NN_get_g(gval)
!!$      gnorm= 0d0
!!$      do iv=1,nvars
!!$        gnorm= gnorm +gval(iv)*gval(iv)
!!$      enddo
!!$      write(6,'(a,2i6,2es15.7)') ' istp,ismpl,f,gnorm='&
!!$           ,istp,ismpl,fval,gnorm
    enddo
  enddo

end subroutine sequential_update
!=======================================================================
subroutine check_grad()
  use variables
  use NN,only:NN_init,NN_func,NN_grad
  use parallel
  implicit none
  integer:: iv
  real(8):: f0,ftmp,dv,vmax
  real(8),allocatable:: ganal(:),gnumer(:),vars0(:)

  allocate(gnumer(nvars),ganal(nvars),vars0(nvars))
  call NN_init()
  f0= NN_func(nvars,vars)
  ganal= NN_grad(nvars,vars)

  vars0(1:nvars)= vars(1:nvars)
  vmax= 0d0
  do iv=1,nvars
    vmax= max(vmax,abs(vars0(iv)))
    if( myid.eq.0) write(6,'(a,i6,es12.4)') ' iv,vars(iv)=',iv,vars0(iv)
  enddo
  dv= vmax *1d-5
  if( myid.eq.0 ) then
    print *,''
    print *,'deviation [dv] =',dv
  endif
  do iv=1,nvars
    vars(1:nvars)= vars0(1:nvars)
!!$    dv= vars(iv) *1d-3
    vars(iv)= vars(iv) +dv
    ftmp= NN_func(nvars,vars)
    gnumer(iv)= (ftmp-f0)/dv
  enddo

  if( myid.eq.0 ) then
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
  endif
end subroutine check_grad
!=======================================================================
subroutine test()
  use variables
  use NN,only:NN_init,NN_func,NN_grad
  use parallel
  implicit none 
  integer:: iv
  real(8):: ft
  real(8),allocatable:: gt(:)

  allocate(gt(nvars))

  call NN_init()
  ft= NN_func(nvars,vars)
  gt= NN_grad(nvars,vars)

  if( myid.eq.0 ) then
    print *,'func value=',ft
    print *,'grad values:'
    do iv=1,nvars
      print *,'iv,grad(iv)=',iv,gt(iv)
    enddo
    print *,'test done.'
  endif
end subroutine test
!=======================================================================
subroutine write_energy_relation(cadd)
  use variables
  use parallel
  implicit none
  character(len=*),intent(in),optional:: cadd
  character(len=128):: cfname
  
  integer:: ismpl
  
  cfname='out.erg.smd-vs-dft'
  if( present(cadd) ) cfname= trim(cfname)//'.'//trim(cadd)

  if( .not. allocated(erefl) ) allocate(erefl(nsmpl),erefg(nsmpl) &
       ,epotl(nsmpl),epotg(nsmpl))

  erefl(1:nsmpl)= 0d0
  epotl(1:nsmpl)= 0d0
  do ismpl=isid0,isid1
    erefl(ismpl)= samples(ismpl)%eref
    epotl(ismpl)= samples(ismpl)%epot
  enddo
  erefg(1:nsmpl)= 0d0
  epotg(1:nsmpl)= 0d0
  call mpi_reduce(epotl,epotg,nsmpl,mpi_double_precision,mpi_sum &
       ,0,mpi_world,ierr)
  call mpi_reduce(erefl,erefg,nsmpl,mpi_double_precision,mpi_sum &
       ,0,mpi_world,ierr)

  if( myid.eq.0 ) then
    open(90,file=trim(cfname),status='replace')
    do ismpl=1,nsmpl
      write(90,'(2es15.7,2x,a)') erefg(ismpl)/nalist(ismpl) &
           ,epotg(ismpl)/nalist(ismpl),cdirlist(ismpl)
    enddo
    close(90)
  endif
  
end subroutine write_energy_relation
!=======================================================================
subroutine write_force_relation(cadd)
  use variables
  use parallel
  implicit none
  character(len=*),intent(in),optional:: cadd
  character(len=128):: cfname

  integer:: ismpl,ia,ixyz,natm,nmax,nmaxl
  
  cfname= 'out.frc.smd-vs-dft'
  if( present(cadd) ) cfname= trim(cfname)//'.'//trim(cadd)

  nmaxl= 0
  do ismpl=1,nsmpl
    nmaxl= max(nmaxl,nalist(ismpl))
  enddo
  call mpi_allreduce(nmaxl,nmax,1,mpi_integer,mpi_max &
       ,mpi_world,ierr)

  if( .not. allocated(frefl) ) allocate(frefl(3,nmax,nsmpl)&
       ,frefg(3,nmax,nsmpl),fal(3,nmax,nsmpl),fag(3,nmax,nsmpl) )

  frefl(1:3,1:nmax,1:nsmpl)= 0d0
  fal(1:3,1:nmax,1:nsmpl)= 0d0
  do ismpl=isid0,isid1
    natm= samples(ismpl)%natm
    frefl(1:3,1:natm,ismpl)= samples(ismpl)%fref(1:3,1:natm)
    fal(1:3,1:natm,ismpl)= samples(ismpl)%fa(1:3,1:natm)
  enddo
  frefg(1:3,1:nmax,1:nsmpl)= 0d0
  fag(1:3,1:nmax,1:nsmpl)= 0d0
  call mpi_reduce(fal,fag,3*nmax*nsmpl,mpi_double_precision,mpi_sum &
       ,0,mpi_world,ierr)
  call mpi_reduce(frefl,frefg,3*nmax*nsmpl,mpi_double_precision,mpi_sum &
       ,0,mpi_world,ierr)

  if( myid.eq.0 ) then
    open(91,file=trim(cfname),status='replace')
    do ismpl=1,nsmpl
      natm= nalist(ismpl)
      do ia=1,natm
        do ixyz=1,3
          write(91,'(2es15.7,2x,a,i6,i3)') frefg(ixyz,ia,ismpl) &
               ,fag(ixyz,ia,ismpl) &
               ,cdirlist(ismpl),ia,ixyz
        enddo
      enddo
    enddo
    close(91)
  endif
  
end subroutine write_force_relation
!=======================================================================
subroutine write_statistics()
  use variables
  use parallel
  implicit none
  integer:: ismpl,ia,l,n,natm
  real(8):: demax,desum,de,rmse,dfmax,dfsum,df

  if( .not.allocated(erefg) ) then
    if(myid.eq.0) then
      print *,'[Error] write_statistics should be called'// &
           ' after calling write_???_relation.'
    endif
    call mpi_finalize(ierr)
    stop
  endif

  if(myid.eq.0) write(6,'(/,a)') '>>>>> statistics:'
  !.....energies
  demax= 0d0
  desum= 0d0
  do ismpl=1,nsmpl
    natm= nalist(ismpl)
    de= abs(epotg(ismpl) -erefg(ismpl))/natm
    demax= max(demax,de)
    desum=desum +de*de/nsmpl
  enddo
  rmse= sqrt(desum)
  if( myid.eq.0 ) then
    write(6,'(a,f12.3,a)') '  RMSE of energies         =',rmse,' eV/atom'
    write(6,'(a,f12.3,a)') '  Max residual of energies =',demax,' eV/atom'
  endif

  !.....forces
  dfmax= 0d0
  dfsum= 0d0
  n= 0
  do ismpl=1,nsmpl
    natm= nalist(ismpl)
    do ia=1,natm
      do l=1,3
        df= abs(fag(l,ia,ismpl)-frefg(l,ia,ismpl))
        dfmax= max(dfmax,df)
        dfsum=dfsum +df*df
        n=n +1
      enddo
    enddo
  enddo
  rmse= sqrt(dfsum/n)
  if(myid.eq.0) then
    write(6,'(a,f12.3,a)') '  RMSE of forces           =',rmse,' eV/A'
    write(6,'(a,f12.3,a)') '  Max residual of forces   =',dfmax,' eV/A'
    print *,''
  endif
end subroutine write_statistics
!=======================================================================
subroutine sync_input()
  use variables
  use parallel
  implicit none
  
  call mpi_bcast(nsmpl,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(nstp,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(iprint,1,mpi_integer,0,mpi_world,ierr)

  call mpi_bcast(cfmethod,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(cmaindir,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(cparfile,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(crunmode,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(cpot,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(cpena,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(clinmin,128,mpi_character,0,mpi_world,ierr)

  call mpi_bcast(eps,1,mpi_double_precision,0,mpi_world,ierr)
  call mpi_bcast(xtol,1,mpi_double_precision,0,mpi_world,ierr)
  call mpi_bcast(ftol,1,mpi_double_precision,0,mpi_world,ierr)
  call mpi_bcast(gtol,1,mpi_double_precision,0,mpi_world,ierr)
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
end subroutine sync_input
!=======================================================================
subroutine get_node2sample()
  use variables
  use parallel
  implicit none
  integer:: n,m,ip

!.....compute num samples per node (nspn)
  n= nsmpl/nnode
  m= nsmpl -n*nnode
  allocate(nspn(nnode),ispn(nnode))
  do ip=1,nnode
    nspn(ip)= n
    if( ip.le.m ) nspn(ip)= nspn(ip) +1
  enddo
  mynsmpl= nspn(myid-1)

  !.....compute start and end of sample-id of this process
  isid0= 0
  isid1= 0
  do ip=1,nnode
    isid1= isid1 +nspn(ip)
    ispn(ip)= isid1 -nspn(ip) +1
  enddo
  isid0= ispn(myid+1)
  isid1= ispn(myid+1) +nspn(myid+1) -1

!!$  write(6,'(a,20i3)') ' myid,isid0,isid1,nspn(:),ispn(:)=' &
!!$       ,myid,isid0,isid1,nspn(1:nnode),ispn(1:nnode)

  if( myid.eq.0 ) print *,'get_node2sample done.'
  return
end subroutine get_node2sample
!=======================================================================
subroutine scatter_samples()
  use variables
  use parallel
  implicit none
  integer:: i,idst,ismpl,itag,natm

  integer,allocatable,save:: istat(:)

  if( nnode.eq.1 ) return

  if( .not. allocated(istat) ) allocate(istat(mpi_status_size))
  
  if( myid.eq.0 ) then ! send from node-0 (root)
    do ismpl=ispn(2),nsmpl
!.....get node-id to which data will be sent
      do i=1,nnode-1
        if( ismpl.lt.ispn(i+1) ) then
          idst= i-1
          exit
        endif
      enddo
      natm= samples(ismpl)%natm
      call mpi_send(samples(ismpl)%natm,1,mpi_integer &
           ,idst,itag,mpi_world,ierr)
      call mpi_send(samples(ismpl)%cdirname,5,mpi_character &
           ,idst,itag,mpi_world,ierr)
      call mpi_send(samples(ismpl)%h0,1,mpi_double_precision &
           ,idst,itag,mpi_world,ierr)
      call mpi_send(samples(ismpl)%h,9,mpi_double_precision &
           ,idst,itag,mpi_world,ierr)
      call mpi_send(samples(ismpl)%epot,1,mpi_double_precision &
           ,idst,itag,mpi_world,ierr)
      call mpi_send(samples(ismpl)%eref,1,mpi_double_precision &
           ,idst,itag,mpi_world,ierr)
      call mpi_send(samples(ismpl)%tag,natm,mpi_double_precision &
           ,idst,itag,mpi_world,ierr)
      call mpi_send(samples(ismpl)%ra,3*natm,mpi_double_precision &
           ,idst,itag,mpi_world,ierr)
      call mpi_send(samples(ismpl)%fa,3*natm,mpi_double_precision &
           ,idst,itag,mpi_world,ierr)
      call mpi_send(samples(ismpl)%fref,3*natm,mpi_double_precision &
           ,idst,itag,mpi_world,ierr)
    enddo
  else ! receive at nodes except root
    do ismpl=isid0,isid1
      call mpi_recv(samples(ismpl)%natm,1,mpi_integer &
           ,mpi_any_source,itag,mpi_world,istat,ierr)
      natm= samples(ismpl)%natm
      call mpi_recv(samples(ismpl)%cdirname,5,mpi_character &
           ,mpi_any_source,itag,mpi_world,istat,ierr)
      call mpi_recv(samples(ismpl)%h0,1,mpi_double_precision &
           ,mpi_any_source,itag,mpi_world,istat,ierr)
      call mpi_recv(samples(ismpl)%h,9,mpi_double_precision &
           ,mpi_any_source,itag,mpi_world,istat,ierr)
      call mpi_recv(samples(ismpl)%epot,1,mpi_double_precision &
           ,mpi_any_source,itag,mpi_world,istat,ierr)
      call mpi_recv(samples(ismpl)%eref,1,mpi_double_precision &
           ,mpi_any_source,itag,mpi_world,istat,ierr)
      call mpi_recv(samples(ismpl)%tag,natm,mpi_double_precision &
           ,mpi_any_source,itag,mpi_world,istat,ierr)
      call mpi_recv(samples(ismpl)%ra,3*natm,mpi_double_precision &
           ,mpi_any_source,itag,mpi_world,istat,ierr)
      call mpi_recv(samples(ismpl)%fa,3*natm,mpi_double_precision &
           ,mpi_any_source,itag,mpi_world,istat,ierr)
      call mpi_recv(samples(ismpl)%fref,3*natm,mpi_double_precision &
           ,mpi_any_source,itag,mpi_world,istat,ierr)
    enddo
  endif
  if( myid.eq.0 ) then
    print *,'scatter_samples done.'
  endif
end subroutine scatter_samples
!=======================================================================
subroutine gather_samples()
  use variables
  use parallel
  implicit none
  integer:: i,isrc,idst,nmpi,ismpl,itag,natm
  integer,allocatable,save:: istat(:,:),ireq(:)

  if( nnode.eq.1 ) return

  if( .not. allocated(ireq) ) allocate(ireq(nsmpl*2))
  if( .not. allocated(istat) ) allocate(istat(mpi_status_size,nsmpl*2))

  if( myid.eq.0 ) then ! send from node-0 (root)
    nmpi= 0
    do ismpl=ispn(2),nsmpl
!.....get node-id to which data will be sent
      do i=1,nnode-1
        if( ismpl.lt.ispn(i+1) ) then
          isrc= i-1
          exit
        endif
      enddo
      natm= samples(ismpl)%natm
      call mpi_irecv(samples(ismpl)%epot,1,mpi_double_precision &
           ,isrc,itag,mpi_world,ireq,ierr)
      call mpi_irecv(samples(ismpl)%fa,3*natm,mpi_double_precision &
           ,isrc,itag,mpi_world,ireq,ierr)
      nmpi= nmpi +2
    enddo
    call mpi_waitall(nmpi,ireq,istat,ierr)
  else ! receive at nodes except root
    nmpi= 0
    idst= 0
    do ismpl=isid0,isid1
      natm= samples(ismpl)%natm
      call mpi_isend(samples(ismpl)%epot,1,mpi_double_precision &
           ,idst,itag,mpi_world,ireq,ierr)
      call mpi_isend(samples(ismpl)%fa,3*natm,mpi_double_precision &
           ,idst,itag,mpi_world,ireq,ierr)
      nmpi= nmpi +2
    enddo
    call mpi_waitall(nmpi,ireq,istat,ierr)
  endif

end subroutine gather_samples
!=======================================================================
subroutine gather_idata(ibufs,ibufr,ndim)
  use parallel
  implicit none
  integer,intent(in):: ndim,ibufs(ndim)
  integer,intent(inout):: ibufr(ndim)
  integer:: i,isrc,idst,nmpi,ismpl,itag,natm
  integer,allocatable,save:: istat(:,:),ireq(:)

  if( nnode.eq.1 ) then
    ibufr(1:ndim)= ibufs(1:ndim)
    return
  endif

  if( .not. allocated(ireq) ) allocate(ireq(nnode))
  if( .not. allocated(istat) ) allocate(istat(mpi_status_size,nnode))

  if( myid.eq.0 ) then ! send from node-0 (root)
    nmpi= 0
    do isrc=1,nnode-1
      call mpi_irecv(ibufr,ndim,mpi_integer &
           ,isrc,itag,mpi_world,ireq,ierr)
      nmpi= nmpi +1
    enddo
    call mpi_waitall(nmpi,ireq,istat,ierr)
  else ! receive at nodes except root
    nmpi= 0
    idst= 0
    call mpi_isend(ibufs,ndim,mpi_integer &
         ,idst,itag,mpi_world,ireq,ierr)
    nmpi= nmpi +1
    call mpi_waitall(nmpi,ireq,istat,ierr)
  endif
  
end subroutine gather_idata
!=======================================================================
subroutine gather_rdata(bufs,bufr,ndim)
  use parallel
  implicit none
  integer,intent(in):: ndim
  real(8),intent(in):: bufs(ndim)
  real(8),intent(inout):: bufr(ndim)
  integer:: i,isrc,idst,nmpi,ismpl,itag,natm
  integer,allocatable,save:: istat(:,:),ireq(:)

  if( nnode.eq.1 ) then
    bufr(1:ndim)= bufs(1:ndim)
    return
  endif

  if( .not. allocated(ireq) ) allocate(ireq(nnode))
  if( .not. allocated(istat) ) allocate(istat(mpi_status_size,nnode))

  if( myid.eq.0 ) then ! send from node-0 (root)
    nmpi= 0
    do isrc=1,nnode-1
      call mpi_irecv(bufr,ndim,mpi_double_precision &
           ,isrc,itag,mpi_world,ireq,ierr)
      nmpi= nmpi +1
    enddo
    call mpi_waitall(nmpi,ireq,istat,ierr)
  else ! receive at nodes except root
    nmpi= 0
    idst= 0
    call mpi_isend(bufs,ndim,mpi_double_precision &
         ,idst,itag,mpi_world,ireq,ierr)
    nmpi= nmpi +1
    call mpi_waitall(nmpi,ireq,istat,ierr)
  endif
  
end subroutine gather_rdata
