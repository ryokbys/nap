program fitpot
!-----------------------------------------------------------------------
!                        Time-stamp: <2015-03-14 11:04:41 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
  use variables
  use parallel
  use minimize
  implicit none
  integer:: ismpl
  real(8):: tmp

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

  !.....divide samples to traning and test sets
  if( myid.eq.0 ) nsmpl_tst= nsmpl*ratio_test
  call mpi_bcast(nsmpl_tst,1,mpi_integer,0,mpi_world,ierr)
  nsmpl_trn= nsmpl -nsmpl_tst
  if( myid.eq.0 ) then
    print *,'nsmpl,train,test=',nsmpl,nsmpl_trn,nsmpl_tst
  endif

  if( nnode.gt.min(nsmpl_trn,nsmpl_tst) ) then
    if( myid.eq.0 ) then
      print *,'[Error] nnode.gt.min(nsmpl_trn,nsmpl_tst)'
      print *,'  nnode,nspml_trn,nsmpl_tst=' &
           ,nnode,nsmpl_trn,nsmpl_tst
      print *,'you should use less number of nodes than nsmpl.'
    endif
    call mpi_finalize(ierr)
    stop
  endif
  call get_node2sample()

  allocate(smpl_trn(isid0_trn:isid1_trn),smpl_tst(isid0_tst:isid1_tst))

  call get_dir_list(11)
!.....store dirname
  do ismpl=isid0_trn,isid1_trn
    smpl_trn(ismpl)%cdirname= cdirlist(ismpl)
  enddo
  do ismpl=isid0_tst,isid1_tst
    smpl_tst(ismpl)%cdirname= cdirlist(ismpl)
  enddo

  call read_samples()
  call read_ref_data()

  call read_vars()
  allocate(gvar(nvars),dvar(nvars))

  select case (trim(cfmethod))
    case ('sd','SD')
      call sd_wrapper()
    case ('cg','CG')
      call cg_wrapper()
    case ('bfgs','BFGS','dfp','DFP')
      call qn_wrapper()
    case ('fs','FS')
      call fs_wrapper()
    case ('sgd','SGD')
      call sgd()
    case ('check_grad')
      call check_grad()
    case ('test','TEST')
      call test()
    case default
      if(myid.eq.0) print *,'unknown fitting_method:',trim(cfmethod)
      call mpi_finalize(ierr)
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
  if(trim(cpena).eq.'lasso' .or. trim(cpena).eq.'glasso' &
       .or.trim(cpena).eq.'ridge') &
       call write_eliminated_vars()
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
  use minimize
  implicit none 
  integer:: i

  write(6,'(a)') '---------------- INITIAL SETTING -----------------'
  write(6,'(2x,a25,2x,i5)') 'num_samples',nsmpl
  write(6,'(2x,a25,2x,i10)') 'num_iteration',niter
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
  write(6,'(2x,a25,2x,es12.3)') 'fscale_factor',fscl
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
    call shuffle_dirlist(nsmpl,cdirlist)
  endif
  call mpi_barrier(mpi_world,ierr)
  call mpi_bcast(cdirlist,5*nsmpl,mpi_character,0,mpi_world,ierr)

  if( myid.eq.0 ) then
    do is=1,nsmpl
      print *,' is,dirlist=',is,cdirlist(is)
    enddo
  endif
  
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
  integer,allocatable:: nal(:)

  if( .not. allocated(nalist) ) allocate(nalist(nsmpl))
  allocate(nal(nsmpl))
  nalist(1:nsmpl)= 0
  nal(1:nsmpl)= 0

!.....training set
  do is=isid0_trn,isid1_trn
!!$    print *,' is=',is,smpl_trn(is)%cdirname
    cdir= smpl_trn(is)%cdirname
    call read_pos(12,trim(cmaindir)//'/'//trim(cdir) &
         //'/pos',is,smpl_trn(is))
    nal(is)= smpl_trn(is)%natm
  enddo
!.....test set
  do is=isid0_tst,isid1_tst
!!$    print *,' is=',is,smpl_tst(is)%cdirname
    cdir= smpl_tst(is)%cdirname
    call read_pos(12,trim(cmaindir)//'/'//trim(cdir) &
         //'/pos',is,smpl_tst(is))
    nal(is)= smpl_tst(is)%natm
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
  deallocate(nal)
  return
end subroutine read_samples
!=======================================================================
subroutine read_pos(ionum,fname,ismpl,smpl)
  use variables
  implicit none 
  integer,intent(in):: ionum,ismpl
  character(len=*),intent(in):: fname
  type(mdsys),intent(inout):: smpl

  integer:: i,natm
  real(8):: tmp

  open(ionum,file=trim(fname),status='old')
  read(ionum,*) smpl%h0
  read(ionum,*) smpl%h(1,1:3)
  read(ionum,*) smpl%h(2,1:3)
  read(ionum,*) smpl%h(3,1:3)
  read(ionum,*) tmp,tmp,tmp
  read(ionum,*) tmp,tmp,tmp
  read(ionum,*) tmp,tmp,tmp
  read(ionum,*) natm
  smpl%natm= natm
  allocate(smpl%ra(3,natm),smpl%fa(3,natm) &
       ,smpl%tag(natm) &
       ,smpl%fref(3,natm))
  do i=1,smpl%natm
    read(ionum,*) smpl%tag(i),smpl%ra(1:3,i), &
         tmp,tmp,tmp
  enddo
  close(ionum)
end subroutine read_pos
!=======================================================================
subroutine read_ref_data()
  use variables
  use parallel
  implicit none 
  integer:: ismpl,i,is,jflag,natm
  character(len=5):: cdir

  jflag= 0
!.....training set
  do ismpl=isid0_trn,isid1_trn
!!$    print *,' ismpl=',ismpl,smpl_trn(ismpl)%cdirname
    cdir=smpl_trn(ismpl)%cdirname
    open(13,file=trim(cmaindir)//'/'//trim(cdir) &
         //'/erg.ref',status='old')
    read(13,*) smpl_trn(ismpl)%eref
    close(13)
!.....reduce atomic energy from eref
    do i=1,smpl_trn(ismpl)%natm
      is= smpl_trn(ismpl)%tag(i)
      smpl_trn(ismpl)%eref= smpl_trn(ismpl)%eref -eatom(is)
    enddo

    open(14,file=trim(cmaindir)//'/'//trim(cdir) &
         //'/frc.ref',status='old')
    read(14,*) natm
    if( natm.ne.smpl_trn(ismpl)%natm ) then
      print *,'Error: natm in sample is not same as smpl%natm'
      print *,' myid,ismpl,natm,smpl%natm=',myid,ismpl &
           ,natm,smpl_trn(ismpl)%natm
      jflag= jflag +1
    endif
    do i=1,natm
      read(14,*) smpl_trn(ismpl)%fref(1:3,i)
    enddo
    close(14)
  enddo
!!$  print *,' reading ref data of training set done.'

!.....test set
  do ismpl=isid0_tst,isid1_tst
!!$    print *,' ismpl=',ismpl,smpl_tst(ismpl)%cdirname
    cdir= smpl_tst(ismpl)%cdirname
    open(13,file=trim(cmaindir)//'/'//trim(cdir) &
         //'/erg.ref',status='old')
    read(13,*) smpl_tst(ismpl)%eref
    close(13)
!.....reduce atomic energy from eref
    do i=1,smpl_tst(ismpl)%natm
      is= smpl_tst(ismpl)%tag(i)
      smpl_tst(ismpl)%eref= smpl_tst(ismpl)%eref -eatom(is)
    enddo

    open(14,file=trim(cmaindir)//'/'//trim(cdir) &
         //'/frc.ref',status='old')
    read(14,*) natm
    if( natm.ne.smpl_tst(ismpl)%natm ) then
      print *,'Error: natm in sample is not same as smpl%natm'
      print *,' myid,ismpl,natm,smpl%natm=',myid,ismpl &
           ,natm,smpl_tst(ismpl)%natm
      jflag= jflag +1
    endif
    do i=1,natm
      read(14,*) smpl_tst(ismpl)%fref(1:3,i)
    enddo
    close(14)
  enddo

!!$  print *,' reading ref data of test set done.'

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
  character(len=*),intent(in):: cadd
  integer:: i
  character(len=128):: cfname

  cfname= trim(cmaindir)//'/'//trim(cparfile)//'.'//trim(cadd)

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
subroutine qn_wrapper()
  use variables
  use NN,only:NN_init,NN_func,NN_grad,NN_restore_standard,NN_analyze
  use parallel
  use minimize
  implicit none
  integer:: i,m
  real(8):: fval

  !.....NN specific code hereafter
  call NN_init()
  do i=1,niter,niter_eval
    call qn(nvars,vars,fval,gvar,dvar,xtol,gtol,ftol,i,i+niter_eval-1 &
         ,iprint,iflag,myid,NN_func,NN_grad,cfmethod)
    call eval_testset(i+niter_eval-1,fval,gvar)
  enddo
  call NN_analyze()
  call NN_restore_standard()

  return
end subroutine qn_wrapper
!=======================================================================
subroutine sd_wrapper()
!
!  Steepest descent minimization
!
  use variables
  use NN,only:NN_init,NN_func,NN_grad
  use parallel
  use minimize
  implicit none
  integer:: i,m
  real(8):: fval

  !.....NN specific code hereafter
  call NN_init()
  call steepest_descent(nvars,vars,fval,gvar,dvar,xtol,gtol,ftol,niter &
       ,iprint,iflag,myid,NN_func,NN_grad)

  return
end subroutine sd_wrapper
!=======================================================================
subroutine cg_wrapper()
  use variables
  use NN,only:NN_init,NN_func,NN_grad
  use parallel
  use minimize
  implicit none
  integer:: i,m
  real(8):: fval

  !.....NN specific code hereafter
  call NN_init()
  call cg(nvars,vars,fval,gvar,dvar,xtol,gtol,ftol,niter &
       ,iprint,iflag,myid,NN_func,NN_grad)

  return
end subroutine cg_wrapper
!=======================================================================
subroutine sgd()
!
! Stochastic gradient decent (SGD)
!
  use variables
  use NN,only:NN_init,NN_fs,NN_gs,NN_func,NN_grad,NN_analyze &
       ,NN_restore_standard
  use parallel
  use minimize
  implicit none
  integer,parameter:: niter_time= 1
  real(8),parameter:: alpha0  = 1d0
  real(8),parameter:: dalpha  = 0.0001d0
  real(8),allocatable:: gval(:),u(:)
  integer:: iter,istp,iv
  real(8):: gnorm,alpha,alpha1,gmax,vmax,fval,gg
  integer:: ismpl
  common /samplei/ ismpl
  real(8),external:: urnd

  allocate(gval(nvars),u(nvars))

!!$  if( nnode.ne.1 ) then
!!$    if(myid.eq.0) print *,'[Error] SGD mode only works '// &
!!$         'with single process.'
!!$    call mpi_finalize(ierr)
!!$    stop
!!$  endif

  alpha1= alpha0

  call NN_init()
  do iter=1,niter
    if(mod(iter,niter_eval).eq.0) then
      fval= NN_func(nvars,vars)
      gval= NN_grad(nvars,vars)
      gnorm= 0d0
      do iv=1,nvars
        gnorm= gnorm +gval(iv)*gval(iv)
      enddo
      if( myid.eq.0 ) then
        write(6,'(a,i6,2es15.7,f10.3)') ' iter,f,gnorm,time=',iter,fval &
             ,gnorm ,mpi_wtime()-time0
        call write_vars('tmp')
      endif
    else if(mod(iter,niter_time).eq.0) then
      if( myid.eq.0 ) then
        write(6,'(a,i6,f10.3)') ' iter,time=',iter,mpi_wtime()-time0
      endif
    endif
    do istp=1,maxmynsmpl_trn
      ismpl= isid0_trn +mynsmpl_trn*urnd()
      fval= NN_fs(nvars,vars)
      gval= NN_gs(nvars,vars)
      u(1:nvars)= -gval(1:nvars)
      alpha= alpha1
      call armijo_search(nvars,vars,u,fval,gval,alpha,iprint &
           ,iflag,myid,NN_fs)
      gnorm= 0d0
      do iv=1,nvars
        gnorm= gnorm +gval(iv)*gval(iv)
      enddo
!!$      if( myid.eq.0 ) then
!!$        write(6,'(a,3es12.4)') 'alpha1,alpha,gnorm=',alpha1,alpha,gnorm
!!$      endif
      vars(1:nvars)=vars(1:nvars) +alpha*u(1:nvars)
    enddo
    alpha1= alpha1*(1d0-dalpha)
  enddo

  fval= NN_func(nvars,vars)
  gval= NN_grad(nvars,vars)
  gnorm= 0d0
  do iv=1,nvars
    gnorm= gnorm +gval(iv)*gval(iv)
  enddo

  call NN_analyze()
  call NN_restore_standard()


end subroutine sgd
!=======================================================================
subroutine fs_wrapper()
  use variables
  use NN,only:NN_init,NN_func,NN_grad,NN_restore_standard,NN_analyze
  use parallel
  use minimize
  implicit none
  integer:: i,m
  real(8):: fval

  !.....NN specific code hereafter
  call NN_init()
  call fs(nvars,vars,fval,gvar,dvar,xtol,gtol,ftol,niter &
       ,iprint,iflag,myid,NN_func,NN_grad)
  call NN_analyze()
  call NN_restore_standard()

  return
end subroutine fs_wrapper
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
  dv= vmax *1d-4
  if( myid.eq.0 ) then
    print *,''
    print *,'deviation [dv] =',dv
  endif
  do iv=1,nvars
    vars(1:nvars)= vars0(1:nvars)
!!$    dv= abs(vars(iv)) *1d-6
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
  use NN,only:NN_init,NN_func,NN_grad,NN_func_tst
  use parallel
  implicit none 
  integer:: iv
  real(8):: ft,ftest
  real(8),allocatable:: gt(:)

  allocate(gt(nvars))

  call NN_init()
  ft= NN_func(nvars,vars)
  ftest= NN_func_tst(nvars,vars)
  gt= NN_grad(nvars,vars)

  if( myid.eq.0 ) then
    print *,'func value     =',ft
    print *,'func_test value=',ftest
    print *,'grad values:'
    do iv=1,nvars
      print *,'iv,grad(iv)=',iv,gt(iv)
    enddo
    print *,'test done.'
  endif

end subroutine test
!=======================================================================
subroutine eval_testset(iter,fv,gv)
  use variables
  use parallel
  use NN, only: NN_func_tst
  implicit none
  integer,intent(in):: iter
  real(8):: fv,gv(nvars)
  
  real(8):: ft
  character(len=5):: cnum

  ft= NN_func_tst(nvars,vars)
  write(cnum,'(i5.5)') iter
  call write_vars(cnum)
  call write_energy_relation(cnum)
  call write_force_relation(cnum)
  call write_statistics()
  
end subroutine eval_testset
!=======================================================================
subroutine write_energy_relation(cadd)
  use variables
  use parallel
  implicit none
  character(len=*),intent(in):: cadd
  character(len=128):: cfname
  
  integer:: ismpl
  
  cfname='out.erg.'//trim(cadd)

  if( .not. allocated(erefl) ) allocate(erefl(nsmpl),erefg(nsmpl) &
       ,epotl(nsmpl),epotg(nsmpl))

  erefl(1:nsmpl)= 0d0
  epotl(1:nsmpl)= 0d0
  do ismpl=isid0_trn,isid1_trn
    erefl(ismpl)= smpl_trn(ismpl)%eref
    epotl(ismpl)= smpl_trn(ismpl)%epot
  enddo
  do ismpl=isid0_tst,isid1_tst
    erefl(ismpl)= smpl_tst(ismpl)%eref
    epotl(ismpl)= smpl_tst(ismpl)%epot
  enddo
  erefg(1:nsmpl)= 0d0
  epotg(1:nsmpl)= 0d0
  call mpi_reduce(epotl,epotg,nsmpl,mpi_double_precision,mpi_sum &
       ,0,mpi_world,ierr)
  call mpi_reduce(erefl,erefg,nsmpl,mpi_double_precision,mpi_sum &
       ,0,mpi_world,ierr)

  if( myid.eq.0 ) then
    open(90,file=trim(cfname),status='replace')
    do ismpl=1,nsmpl_trn
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
  character(len=*),intent(in):: cadd
  character(len=128):: cfname

  integer:: ismpl,ia,ixyz,natm,nmax,nmaxl
  
  cfname= 'out.frc.'//trim(cadd)

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
  do ismpl=isid0_trn,isid1_trn
    natm= smpl_trn(ismpl)%natm
    frefl(1:3,1:natm,ismpl)= smpl_trn(ismpl)%fref(1:3,1:natm)
    fal(1:3,1:natm,ismpl)= smpl_trn(ismpl)%fa(1:3,1:natm)
  enddo
  do ismpl=isid0_tst,isid1_tst
    natm= smpl_tst(ismpl)%natm
    frefl(1:3,1:natm,ismpl)= smpl_tst(ismpl)%fref(1:3,1:natm)
    fal(1:3,1:natm,ismpl)= smpl_tst(ismpl)%fa(1:3,1:natm)
  enddo
  frefg(1:3,1:nmax,1:nsmpl)= 0d0
  fag(1:3,1:nmax,1:nsmpl)= 0d0
  call mpi_reduce(fal,fag,3*nmax*nsmpl,mpi_double_precision,mpi_sum &
       ,0,mpi_world,ierr)
  call mpi_reduce(frefl,frefg,3*nmax*nsmpl,mpi_double_precision,mpi_sum &
       ,0,mpi_world,ierr)

  if( myid.eq.0 ) then
    open(91,file=trim(cfname),status='replace')
    do ismpl=1,nsmpl_trn
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
  real(8):: demax2,desum2,de2,rmse2,dfmax2,dfsum2,df2
  logical,save:: l1st=.true.

  if( .not.allocated(erefg) ) then
    if(myid.eq.0) then
      print *,'[Error] write_statistics should be called'// &
           ' after calling write_???_relation.'
    endif
    call mpi_finalize(ierr)
    stop
  endif

  if(myid.eq.0) write(6,'(a)') '>>>>> statistics:'
  !.....energies of training set
  demax= 0d0
  desum= 0d0
  do ismpl=1,nsmpl_trn
    natm= nalist(ismpl)
    de= abs(epotg(ismpl) -erefg(ismpl))/natm
    demax= max(demax,de)
    desum=desum +de*de/nsmpl_trn
  enddo
  rmse= sqrt(desum)
  !.....energies of test set
  demax2= 0d0
  desum2= 0d0
  do ismpl=nsmpl_trn+1,nsmpl
    natm= nalist(ismpl)
    de2= abs(epotg(ismpl) -erefg(ismpl))/natm
    demax2= max(demax2,de2)
    desum2=desum2 +de2*de2/nsmpl_tst
  enddo
  rmse2= sqrt(desum2)
  if( myid.eq.0 ) then
    write(6,'(a)') ' Training set,'
    write(6,'(a,f12.3,a)') '  RMSE of energies         =',rmse,' eV/atom'
    write(6,'(a,f12.3,a)') '  Max residual of energies =',demax,' eV/atom'
    write(6,'(a)') ' Test set,'
    write(6,'(a,f12.3,a)') '  RMSE of energies         =',rmse2,' eV/atom'
    write(6,'(a,f12.3,a)') '  Max residual of energies =',demax2,' eV/atom'
    write(6,'(a,4f15.7)') '  energy:training(rmse,max),test(rmse,max)=' &
         ,rmse,demax,rmse2,demax2
  endif

  !.....forces of training set
  dfmax= 0d0
  dfsum= 0d0
  n= 0
  do ismpl=1,nsmpl_trn
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
  !.....forces of test set
  dfmax2= 0d0
  dfsum2= 0d0
  n= 0
  do ismpl=nsmpl_trn+1,nsmpl
    natm= nalist(ismpl)
    do ia=1,natm
      do l=1,3
        df2= abs(fag(l,ia,ismpl)-frefg(l,ia,ismpl))
        dfmax2= max(dfmax2,df2)
        dfsum2=dfsum2 +df2*df2
        n=n +1
      enddo
    enddo
  enddo
  rmse2= sqrt(dfsum2/n)
  if(myid.eq.0) then
    write(6,'(a)') ' Training set,'
    write(6,'(a,f12.3,a)') '  RMSE of forces           =',rmse,' eV/A'
    write(6,'(a,f12.3,a)') '  Max residual of forces   =',dfmax,' eV/A'
    write(6,'(a)') ' Test set,'
    write(6,'(a,f12.3,a)') '  RMSE of forces           =',rmse2,' eV/A'
    write(6,'(a,f12.3,a)') '  Max residual of forces   =',dfmax2,' eV/A'
    write(6,'(a,4f15.7)') '  force:training(rmse,max),test(rmse,max)=' &
         ,rmse,dfmax,rmse2,dfmax2
    print *,''
  endif
end subroutine write_statistics
!=======================================================================
subroutine write_eliminated_vars()
  use variables
  use parallel
  implicit none
  integer:: i,i0

  i0= 0
  do i=1,nvars
    if( abs(vars(i)).lt.1d-8 ) then
      i0=i0+1
    endif
  enddo
  if(myid.eq.0) write(6,'(a,i6,a,i6)') ' num of 0-vars = ',i0,'/',nvars
end subroutine write_eliminated_vars
!=======================================================================
subroutine sync_input()
  use variables
  use parallel
  use minimize
  implicit none
  
  call mpi_bcast(nsmpl,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(niter,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(niter_eval,1,mpi_integer,0,mpi_world,ierr)
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

!.....compute num samples for training per node (nspn)
  n= nsmpl_trn/nnode
  m= nsmpl_trn -n*nnode
  allocate(nspn_trn(nnode),ispn_trn(nnode))
  do ip=1,nnode
    nspn_trn(ip)= n
    if( ip.le.m ) nspn_trn(ip)= nspn_trn(ip) +1
  enddo
  mynsmpl_trn= nspn_trn(myid+1)
  call mpi_allreduce(mynsmpl_trn,maxmynsmpl_trn,1,mpi_integer,mpi_max &
       ,mpi_world,ierr)
!!$  print *,' myid,mynsmpl_trn,maxmynsmpl_trn=' &
!!$       ,myid,mynsmpl_trn,maxmynsmpl_trn
  if( myid.eq.0 ) print *,'maxmynsmpl_trn=',maxmynsmpl_trn

  !.....compute start and end of sample-id of this process
  isid0_trn= 0
  isid1_trn= 0
  do ip=1,nnode
    isid1_trn= isid1_trn +nspn_trn(ip)
    ispn_trn(ip)= isid1_trn -nspn_trn(ip) +1
  enddo
  isid0_trn= ispn_trn(myid+1)
  isid1_trn= ispn_trn(myid+1) +nspn_trn(myid+1) -1

!!$  write(6,'(a,20i3)') ' myid,isid0,isid1,nspn(:),ispn(:)=' &
!!$       ,myid,isid0,isid1,nspn(1:nnode),ispn(1:nnode)

!.....compute num samples for test per node (nspn)
  n= nsmpl_tst/nnode
  m= nsmpl_tst -n*nnode
  allocate(nspn_tst(nnode),ispn_tst(nnode))
  do ip=1,nnode
    nspn_tst(ip)= n
    if( ip.le.m ) nspn_tst(ip)= nspn_tst(ip) +1
  enddo
  mynsmpl_tst= nspn_tst(myid+1)
  call mpi_allreduce(mynsmpl_tst,maxmynsmpl_tst,1,mpi_integer,mpi_max &
       ,mpi_world,ierr)
!!$  print *,' myid,mynsmpl_tst,maxmynsmpl_tst=' &
!!$       ,myid,mynsmpl_tst,maxmynsmpl_tst
  if( myid.eq.0 ) print *,'maxmynsmpl_tst=',maxmynsmpl_tst

  !.....compute start and end of sample-id of this process
  isid0_tst= nsmpl_trn
  isid1_tst= nsmpl_trn
  do ip=1,nnode
    isid1_tst= isid1_tst +nspn_tst(ip)
    ispn_tst(ip)= isid1_tst -nspn_tst(ip) +1
  enddo
  isid0_tst= ispn_tst(myid+1)
  isid1_tst= ispn_tst(myid+1) +nspn_tst(myid+1) -1


  if( myid.eq.0 ) print *,'get_node2sample done.'
  return
end subroutine get_node2sample
!=======================================================================
function urnd()
!
!  Uniform random number generator
!      
  implicit none 
  real(8):: urnd
  real(8),save:: dseed= 12345d0
  real(8),save:: d2p31m,d2p31
  data d2p31m/2147483647d0/
  data d2p31 /2147483648d0/

  dseed=dmod(16807d0*dseed,d2p31m)
  urnd=dseed/d2p31
  return
end function urnd
!=======================================================================
subroutine shuffle_dirlist(nsmpl,cdirlist)
!
!  Randomize the order of the cdirlist
!
  implicit none
  integer,intent(in):: nsmpl
  character(len=*):: cdirlist(nsmpl)
  real(8),external:: urnd
  integer:: i,j,k,n
  character(len=5),allocatable:: cdltmp(:)
  
  allocate(cdltmp(nsmpl))
  cdltmp(1:nsmpl)= cdirlist(1:nsmpl)
  n=nsmpl
  do i=1,nsmpl
    j= n*urnd()+1
    cdirlist(i)= cdltmp(j)
    do k=j,n-1
      cdltmp(k)= cdltmp(k+1)
    enddo
    n= n -1
  enddo
  deallocate(cdltmp)
end subroutine shuffle_dirlist

