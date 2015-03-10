program fitpot
  use variables
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

  call read_input(10,'in.fitpot')
  allocate(samples(nsmpl))
  call write_initial_setting()
  call get_dir_list(11)

  call read_samples()
  call read_ref_data()

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
      print *,'unknown fitting_method:',trim(cfmethod)
      stop
  end select

  if( iflag/1000.ne.0 ) then
    print *,'something wrong with 1D line search.'
    print *,'  iflag=',iflag
  else if( iflag/100.ne.0 ) then
    print *,'something wrong with minimization.'
    print *,'  iflag=',iflag
  endif
  
  call write_vars('fin')
  call write_energy_relation('fin')
  call write_force_relation('fin')

  call write_statistics()
  
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
  implicit none
  integer,intent(in):: ionum
  integer:: i

  call system('ls '//trim(cmaindir) &
       //' | grep "^[0-9]...." > dir_list.txt')

  open(ionum,file='dir_list.txt',status='old')
  do i=1,nsmpl
    read(ionum,*,end=999) samples(i)%cdirname
  enddo
  close(ionum)

!!$  do i=1,nsmpl
!!$    print *,' i,cdirlist(i)=',i,cdirlist(i)
!!$  enddo
  return

999 continue
  print *,' Error: num_samples may be wrong.'
  stop

end subroutine get_dir_list
!=======================================================================
subroutine read_samples()
  use variables
  implicit none
  integer:: ismpl
  character*5:: cdir

  do ismpl=1,nsmpl
!    print *,' ismpl=',ismpl,cdirlist(ismpl)
    cdir= samples(ismpl)%cdirname
    call read_pos(12,trim(cmaindir)//'/'//trim(cdir) &
         //'/pos',ismpl)
  enddo

  print *,'read_samples done.'

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
  integer:: ismpl,i,is

  do ismpl=1,nsmpl
    open(13,file=trim(cmaindir)//'/'//samples(ismpl)%cdirname &
         //'/erg.ref',status='old')
    read(13,*) samples(ismpl)%eref
    close(13)
!.....reduce atomic energy from eref
    do i=1,samples(ismpl)%natm
      is= samples(ismpl)%tag(i)
      eref= eref -eatom(is)
    enddo

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

  print *,'read_ref_data done.'
  
end subroutine read_ref_data
!=======================================================================
subroutine read_vars()
  use variables
  implicit none

  integer:: i

  open(15,file=trim(cmaindir)//'/'//cparfile,status='old')
  read(15,*) nvars, rcut
  allocate(vars(nvars),vranges(2,nvars))
  do i=1,nvars
    read(15,*) vars(i),vranges(1:2,i)
!    print *, vars(i),vranges(1:2,i)
  enddo

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
  use NN,only:NN_init,NN_func,NN_grad
  implicit none
  integer,parameter:: iprint = 1
  integer:: i,m
  real(8):: fval

  !.....NN specific code hereafter
  call NN_init()
  call bfgs(nvars,vars,fval,xtol,gtol,ftol,nstp &
       ,iprint,iflag,NN_func,NN_grad)

  return
end subroutine bfgs_wrapper
!=======================================================================
subroutine sd_wrapper()
!
!  Steepest descent minimization
!
  use variables
  use NN,only:NN_init,NN_func,NN_grad
  implicit none
  integer,parameter:: iprint = 1
  integer:: i,m
  real(8):: fval

  !.....NN specific code hereafter
  call NN_init()
  call steepest_descent(nvars,vars,fval,xtol,gtol,ftol,nstp&
       ,iprint,iflag,NN_func,NN_grad)

  return
end subroutine sd_wrapper
!=======================================================================
subroutine cg_wrapper()
  use variables
  use NN,only:NN_init,NN_func,NN_grad
  implicit none
  integer,parameter:: iprint = 1
  integer:: i,m
  real(8):: fval

  !.....NN specific code hereafter
  call NN_init()
  call bfgs(nvars,vars,fval,xtol,gtol,ftol,nstp &
       ,iprint,iflag,NN_func,NN_grad)

  return
end subroutine cg_wrapper
!=======================================================================
subroutine sequential_update()
  use variables
  use NN,only:NN_init,NN_get_fs,NN_get_gs,NN_func,NN_grad
  implicit none
  real(8),allocatable:: gval(:)
  integer:: istp,ismpl,iv
  real(8):: gnorm,dsc,gmax,vmax,fval

  allocate(gval(nvars))

  call NN_init()

  do istp=1,nstp
    do ismpl=1,nsmpl
      call NN_get_fs(ismpl,fval)
      call NN_get_gs(ismpl,gval)
      gmax= 0d0
      vmax= 0d0
      do iv=1,nvars
        gmax= max(gmax,abs(gval(iv)))
        vmax= max(vmax,abs(vars(iv)))
      enddo
      dsc= min(vmax/gmax,seqcoef)
      vars(1:nvars)=vars(1:nvars) -dsc*gval(1:nvars)
!!$      call NN_get_f(fval)
!!$      call NN_get_g(gval)
!!$      gnorm= 0d0
!!$      do iv=1,nvars
!!$        gnorm= gnorm +gval(iv)*gval(iv)
!!$      enddo
!!$      write(6,'(a,2i6,2es15.7)') ' istp,ismpl,f,gnorm='&
!!$           ,istp,ismpl,fval,gnorm
    enddo
    fval= NN_func(nvars,vars)
    gval= NN_grad(nvars,vars)
    gnorm= 0d0
    do iv=1,nvars
      gnorm= gnorm +gval(iv)*gval(iv)
    enddo
    write(6,'(a,i6,2es15.7)') ' istp,f,gnorm=',istp,fval,gnorm
  enddo

end subroutine sequential_update
!=======================================================================
subroutine check_grad()
  use variables
  use NN,only:NN_init,NN_func,NN_grad
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
    write(6,'(a,i6,es12.4)') ' iv,vars(iv)=',iv,vars0(iv)
  enddo
  dv= vmax *1d-5
  print *,''
  print *,'deviation [dv] =',dv
  do iv=1,nvars
    vars(1:nvars)= vars0(1:nvars)
!!$    dv= vars(iv) *1d-3
    vars(iv)= vars(iv) +dv
    ftmp= NN_func(nvars,vars)
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
  use NN,only:NN_init,NN_func,NN_grad
  implicit none 
  integer:: iv
  real(8):: ft
  real(8),allocatable:: gt(:)

  allocate(gt(nvars))

  call NN_init()
  ft= NN_func(nvars,vars)
  gt= NN_grad(nvars,vars)
  
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
