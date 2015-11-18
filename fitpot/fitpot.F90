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

  if( nnode.gt.nsmpl ) then
    if( myid.eq.0 ) then
      print *,'[Error] nnode.gt.nsmpl'
      print *,'  nnode,nspml=',nnode,nsmpl
      print *,'you should use less number of nodes than nsmpl.'
    endif
    call mpi_finalize(ierr)
    stop
  endif
  call get_node2sample()

  allocate(samples(isid0:isid1))

  call get_dir_list(11)
!.....store dirname
  do ismpl=isid0,isid1
    samples(ismpl)%cdirname= cdirlist(ismpl)
  enddo
  call set_training_test()

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
    case ('gfs')
      call gfs_wrapper()
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
  call write_stats(niter)
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
    write(6,'(a,i10)') ' num of func calls=',nfunc
    write(6,'(a,i10)') ' num of grad calls=',ngrad
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
  write(6,'(2x,a25,2x,es12.3)') 'eps_energy',epse
  write(6,'(2x,a25,2x,es12.3)') 'eps_force',epsf
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

!!$  if( myid.eq.0 ) then
!!$    do is=1,nsmpl
!!$      print *,' is,dirlist=',is,cdirlist(is)
!!$    enddo
!!$  endif
  
  if(myid.eq.0) print*,'get_dir_list done.'
  return

999 continue
  if( myid.eq.0 ) print *,' Error: num_samples may be wrong.'
  call mpi_finalize(ierr)
  stop

end subroutine get_dir_list
!=======================================================================
subroutine set_training_test()
  use variables
  use parallel
  implicit none
  integer:: ismpl,n
  integer,allocatable,dimension(:):: icll

  allocate(icll(nsmpl),iclist(nsmpl))
  icll(1:nsmpl)= 0

  nsmpl_tst= nsmpl*ratio_test
  nsmpl_trn= nsmpl -nsmpl_tst
  call mpi_bcast(nsmpl_trn,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(nsmpl_tst,1,mpi_integer,0,mpi_world,ierr)

  if( myid.lt.mod(nsmpl_tst,nnode) ) then
    myntst= nsmpl_tst/nnode +1
  else
    myntst= nsmpl_tst/nnode
  endif
  myntrn= mynsmpl -myntst
  call mpi_allreduce(myntrn,maxmyntrn,1,mpi_integer,mpi_max &
       ,mpi_world,ierr)

!!$  print *,' myid,myntrn,myntst=',myid,myntrn,myntst
  n=0
  do ismpl=isid0,isid1
    n=n+1
    if( n.le.myntrn ) then
      samples(ismpl)%iclass= 1   ! training data
    else
      samples(ismpl)%iclass= 2   ! test data
    endif
    icll(ismpl)= samples(ismpl)%iclass
  enddo

  iclist(1:nsmpl)= 0
  call mpi_allreduce(icll,iclist,nsmpl,mpi_integer,mpi_max &
       ,mpi_world,ierr)
  
  if( myid.eq.0 ) then
!!$    do ismpl=1,nsmpl
!!$      print *,'ismpl,iclist=',ismpl,iclist(ismpl)
!!$    enddo
    print *,'nsmpl, training, test=',nsmpl,nsmpl_trn,nsmpl_tst
    print *,'set_training_test done.'
  endif
  deallocate(icll)
  return
end subroutine set_training_test
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

  do is=isid0,isid1
    cdir= samples(is)%cdirname
    call read_pos(12,trim(cmaindir)//'/'//trim(cdir) &
         //'/pos',is,samples(is))
    nal(is)= samples(is)%natm
  enddo
  call mpi_reduce(nal,nalist,nsmpl,mpi_integer,mpi_sum &
       ,0,mpi_world,ierr)
  
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
  real(8):: erefminl

  jflag= 0
  erefminl= 0d0
  do ismpl=isid0,isid1
    cdir=samples(ismpl)%cdirname
    open(13,file=trim(cmaindir)//'/'//trim(cdir) &
         //'/erg.ref',status='old')
    read(13,*) samples(ismpl)%eref
    close(13)
!.....reduce atomic energy from eref
    do i=1,samples(ismpl)%natm
      is= samples(ismpl)%tag(i)
      samples(ismpl)%eref= samples(ismpl)%eref -eatom(is)
    enddo
    erefminl= min(erefminl,samples(ismpl)%eref/samples(ismpl)%natm)

    open(14,file=trim(cmaindir)//'/'//trim(cdir) &
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

  erefmin= 0d0
  call mpi_allreduce(erefminl,erefmin,1,mpi_double_precision,mpi_min &
       ,mpi_world,ierr)
  

  if(myid.eq.0) then
    write(6,'(a,es12.4)') ' erefmin = ',erefmin
    print *,'read_ref_data done.'
  endif

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

!!$  cfname= trim(cmaindir)//'/'//trim(cparfile)//'.'//trim(cadd)
  cfname= trim(cparfile)//'.'//trim(cadd)

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
  external:: write_stats

  !.....NN specific code hereafter
  call NN_init()
  call qn(nvars,vars,fval,gvar,dvar,xtol,gtol,ftol,niter &
       ,iprint,iflag,myid,NN_func,NN_grad,cfmethod &
       ,niter_eval,write_stats)
  call NN_analyze("fin")
  if( cpena.eq.'lasso' .or. cpena.eq.'glasso' ) then
    call NN_restore_standard()
  endif

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
  real(8),parameter:: dalpha  = 0.00001d0
!!$  real(8),parameter:: dalpha  = 0.d0
  real(8),allocatable:: gval(:),u(:)
  integer:: iter,istp,iv
  real(8):: gnorm,alpha,alpha1,gmax,vmax,fval,gg
  integer:: ismpl
  common /samplei/ ismpl
  real(8),external:: urnd

  allocate(gval(nvars),u(nvars))

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
      call write_stats(iter)
    else if(mod(iter,niter_time).eq.0) then
      if( myid.eq.0 ) then
        write(6,'(a,i6,f10.3)') ' iter,time=',iter,mpi_wtime()-time0
      endif
    endif
    do istp=1,maxmynsmpl
      ismpl= isid0 +mynsmpl*urnd()
      fval= NN_fs(nvars,vars)
      gval= NN_gs(nvars,vars)
      u(1:nvars)= -gval(1:nvars)
!!$      print '(a,3i5,es12.4)',' istp,ismpl,myid,fval=',istp,ismpl,myid,fval
      alpha= alpha1
      call armijo_search(nvars,vars,u,fval,gval,alpha,iprint &
           ,iflag,myid,NN_fs)
      gnorm= 0d0
      do iv=1,nvars
        gnorm= gnorm +gval(iv)*gval(iv)
      enddo
!!$      if(myid.eq.0) print '(a,i5,3es12.4)',' istp,fval,gnorm,alpha=',istp,fval,gnorm,alpha
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

  call NN_analyze("fin")
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
  call NN_analyze("fin")
  call NN_restore_standard()

  return
end subroutine fs_wrapper
!=======================================================================
subroutine gfs_wrapper()
  use variables
  use NN,only:NN_init,NN_func,NN_grad,NN_restore_standard,NN_analyze
  use parallel
  use minimize
  implicit none
  integer:: i,m
  real(8):: fval
  external:: write_stats,analyze_wrapper

  !.....NN specific code hereafter
  call NN_init()
  call gfs(nvars,vars,fval,gvar,dvar,xtol,gtol,ftol,niter &
       ,iprint,iflag,myid,NN_func,NN_grad,cfmethod,niter_eval &
       ,write_stats,analyze_wrapper)
  call NN_analyze("fin")
  call NN_restore_standard()

  return
end subroutine gfs_wrapper
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

  call write_stats(0)

  if( myid.eq.0 ) then
    print *,'func value     =',ft
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
  character(len=*),intent(in):: cadd
  character(len=128):: cfname
  
  integer:: ismpl,n
  
  cfname='out.erg.'//trim(cadd)

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
    open(90,file=trim(cfname)//'.1',status='replace')
    open(91,file=trim(cfname)//'.2',status='replace')
    do ismpl=1,nsmpl
      erefg(ismpl)= erefg(ismpl)/nalist(ismpl)
      epotg(ismpl)= epotg(ismpl)/nalist(ismpl)
      if( iclist(ismpl).eq.1 ) then
        if( lswgt ) then
          write(90,'(2es15.7,2x,a,es15.7)') erefg(ismpl) &
               ,epotg(ismpl),cdirlist(ismpl) &
               ,exp(-(erefg(ismpl)-erefmin)/abs(erefmin)*swbeta)
        else
          write(90,'(2es15.7,2x,a,es15.7)') erefg(ismpl) &
               ,epotg(ismpl),cdirlist(ismpl) &
               ,abs(erefg(ismpl)-epotg(ismpl))
        endif
      else if( iclist(ismpl).eq.2 ) then
        if( lswgt ) then
          write(91,'(2es15.7,2x,a,es15.7)') erefg(ismpl) &
               ,epotg(ismpl),cdirlist(ismpl) &
               ,exp(-(erefg(ismpl)-erefmin)/abs(erefmin)*swbeta)
        else
          write(91,'(2es15.7,2x,a,es15.7)') erefg(ismpl) &
               ,epotg(ismpl),cdirlist(ismpl) &
               ,abs(erefg(ismpl)-epotg(ismpl))
        endif
!!$        write(91,'(2es15.7,2x,a)') erefg(ismpl)/nalist(ismpl) &
!!$             ,epotg(ismpl)/nalist(ismpl),cdirlist(ismpl)
      endif
    enddo
    close(90)
    close(91)
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
    open(92,file=trim(cfname)//'.1',status='replace')
    open(93,file=trim(cfname)//'.2',status='replace')
    do ismpl=1,nsmpl
      if( iclist(ismpl).eq.1 ) then
        natm= nalist(ismpl)
        do ia=1,natm
          do ixyz=1,3
            write(92,'(2es15.7,2x,a,i6,i3,es15.7)') frefg(ixyz,ia,ismpl) &
                 ,fag(ixyz,ia,ismpl) &
                 ,cdirlist(ismpl),ia,ixyz &
                 ,abs(frefg(ixyz,ia,ismpl)-fag(ixyz,ia,ismpl))
          enddo
        enddo
      else if( iclist(ismpl).eq.2 ) then
        natm= nalist(ismpl)
        do ia=1,natm
          do ixyz=1,3
            write(93,'(2es15.7,2x,a,i6,i3,es15.7)') frefg(ixyz,ia,ismpl) &
                 ,fag(ixyz,ia,ismpl) &
                 ,cdirlist(ismpl),ia,ixyz &
                 ,abs(frefg(ixyz,ia,ismpl)-fag(ixyz,ia,ismpl))
          enddo
        enddo
      endif
    enddo
    close(92)
    close(93)
  endif
  
end subroutine write_force_relation
!=======================================================================
subroutine write_stats(iter)
  use variables
  use parallel
  use NN
  implicit none
  integer,intent(in):: iter
  integer:: ismpl,natm,ntrnl,ntstl,ia,l,ntrn,ntst
  type(mdsys)::smpl
  real(8):: de,df
  real(8):: demaxl_trn,demax_trn,desuml_trn,desum_trn,rmse_trn
  real(8):: demaxl_tst,demax_tst,desuml_tst,desum_tst,rmse_tst
  real(8):: dfmaxl_trn,dfmax_trn,dfsuml_trn,dfsum_trn
  real(8):: dfmaxl_tst,dfmax_tst,dfsuml_tst,dfsum_tst

  demaxl_trn= 0d0
  desuml_trn= 0d0
  demaxl_tst= 0d0
  desuml_tst= 0d0
  do ismpl=isid0,isid1
    smpl= samples(ismpl)
    natm= smpl%natm
    de= abs(smpl%epot -smpl%eref)/natm
    if( smpl%iclass.eq.1 ) then
      demaxl_trn= max(demaxl_trn,de)
      desuml_trn=desuml_trn +de*de
    else if( smpl%iclass.eq.2 ) then
      demaxl_tst= max(demaxl_tst,de)
      desuml_tst=desuml_tst +de*de
    endif
  enddo
  desum_trn= 0d0
  desum_tst= 0d0
  call mpi_reduce(desuml_trn,desum_trn,1 &
       ,mpi_double_precision,mpi_sum,0,mpi_world,ierr)
  call mpi_reduce(desuml_tst,desum_tst,1 &
       ,mpi_double_precision,mpi_sum,0,mpi_world,ierr)
  call mpi_reduce(demaxl_trn,demax_trn,1 &
       ,mpi_double_precision,mpi_max,0,mpi_world,ierr)
  call mpi_reduce(demaxl_tst,demax_tst,1 &
       ,mpi_double_precision,mpi_max,0,mpi_world,ierr)
  rmse_trn= sqrt(desum_trn/nsmpl_trn)
  rmse_tst= sqrt(desum_tst/nsmpl_tst)
  if( myid.eq.0 ) then
    write(6,'(a,i8,f15.2,4f12.7)') '  energy:training(rmse,max)' &
         //',test(rmse,max)=',iter,mpi_wtime()-time0 &
         ,rmse_trn,demax_trn,rmse_tst,demax_tst
  endif

!.....force
  dfmaxl_trn= 0d0
  dfsuml_trn= 0d0
  dfmaxl_tst= 0d0
  dfsuml_tst= 0d0
  ntrnl= 0
  ntstl= 0
  do ismpl=isid0,isid1
    smpl= samples(ismpl)
    natm= smpl%natm
    if( smpl%iclass.eq.1 ) then
      do ia=1,natm
        do l=1,3
          df= abs(smpl%fa(l,ia)-smpl%fref(l,ia))
          dfmaxl_trn= max(dfmaxl_trn,df)
          dfsuml_trn=dfsuml_trn +df*df
          ntrnl=ntrnl +1
        enddo
      enddo
    else if( smpl%iclass.eq.2 ) then
      do ia=1,natm
        do l=1,3
          df= abs(smpl%fa(l,ia)-smpl%fref(l,ia))
          dfmaxl_tst= max(dfmaxl_tst,df)
          dfsuml_tst=dfsuml_tst +df*df
          ntstl=ntstl +1
        enddo
      enddo
    endif
  enddo
  dfsum_trn= 0d0
  dfsum_tst= 0d0
  call mpi_reduce(dfsuml_trn,dfsum_trn,1 &
       ,mpi_double_precision,mpi_sum,0,mpi_world,ierr)
  call mpi_reduce(dfsuml_tst,dfsum_tst,1 &
       ,mpi_double_precision,mpi_sum,0,mpi_world,ierr)
  call mpi_reduce(dfmaxl_trn,dfmax_trn,1 &
       ,mpi_double_precision,mpi_max,0,mpi_world,ierr)
  call mpi_reduce(dfmaxl_tst,dfmax_tst,1 &
       ,mpi_double_precision,mpi_max,0,mpi_world,ierr)
  ntrn= 0
  ntst= 0
  call mpi_reduce(ntrnl,ntrn,1 &
       ,mpi_integer,mpi_sum,0,mpi_world,ierr)
  call mpi_reduce(ntstl,ntst,1 &
       ,mpi_integer,mpi_sum,0,mpi_world,ierr)
  rmse_trn= sqrt(dfsum_trn/ntrn)
  rmse_tst= sqrt(dfsum_tst/ntst)
  if( myid.eq.0 ) then
    write(6,'(a,i8,f15.2,4f12.7)') '  force:training(rmse,max)' &
         //',test(rmse,max)=',iter,mpi_wtime()-time0 &
         ,rmse_trn,dfmax_trn,rmse_tst,dfmax_tst
!    call write_vars('tmp')
  endif

end subroutine write_stats
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
subroutine analyze_wrapper(num)
  use NN
  implicit none 
  integer,intent(in):: num
  character(len=5):: cadd

  write(cadd,'(i5.5)') num
  call NN_analyze(cadd)

end subroutine analyze_wrapper
!=======================================================================
subroutine sync_input()
  use variables
  use parallel
  use minimize
  implicit none
  
  call mpi_bcast(nsmpl,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(niter,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(niter_eval,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(nitergfs,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(iprint,1,mpi_integer,0,mpi_world,ierr)

  call mpi_bcast(cfmethod,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(cmaindir,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(cparfile,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(crunmode,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(cpot,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(cpena,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(clinmin,128,mpi_character,0,mpi_world,ierr)

  call mpi_bcast(epse,1,mpi_double_precision,0,mpi_world,ierr)
  call mpi_bcast(epsf,1,mpi_double_precision,0,mpi_world,ierr)
  call mpi_bcast(xtol,1,mpi_double_precision,0,mpi_world,ierr)
  call mpi_bcast(ftol,1,mpi_double_precision,0,mpi_world,ierr)
  call mpi_bcast(gtol,1,mpi_double_precision,0,mpi_world,ierr)
  call mpi_bcast(eatom,maxnsp,mpi_double_precision,0,mpi_world,ierr)
  call mpi_bcast(gscl,1,mpi_double_precision,0,mpi_world,ierr)
  call mpi_bcast(fscl,1,mpi_double_precision,0,mpi_world,ierr)
  call mpi_bcast(swbeta,1,mpi_double_precision,0,mpi_world,ierr)
  call mpi_bcast(pwgt,1,mpi_double_precision,0,mpi_world,ierr)
  call mpi_bcast(ratio_test,1,mpi_double_precision,0,mpi_world,ierr)
  
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
  n= nsmpl/nnode
  m= nsmpl -n*nnode
  allocate(nspn(nnode),ispn(nnode))
  do ip=1,nnode
    nspn(ip)= n
    if( ip.le.m ) nspn(ip)= nspn(ip) +1
  enddo
  mynsmpl= nspn(myid+1)
  call mpi_allreduce(mynsmpl,maxmynsmpl,1,mpi_integer,mpi_max &
       ,mpi_world,ierr)
  if( myid.eq.0 ) print *,'maxmynsmpl=',maxmynsmpl

  !.....compute start and end of sample-id of this process
  isid0= 0
  isid1= 0
  do ip=1,nnode
    isid1= isid1 +nspn(ip)
    ispn(ip)= isid1 -nspn(ip) +1
  enddo
  isid0= ispn(myid+1)
  isid1= ispn(myid+1) +nspn(myid+1) -1

!!$  print *,'myid,isid0,isid1=',myid,isid0,isid1

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

