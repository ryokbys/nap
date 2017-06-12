program fitpot
!-----------------------------------------------------------------------
!                     Last modified: <2017-06-13 07:54:30 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
  use variables
  use parallel
  use minimize
  implicit none
  integer:: ismpl,ihour,imin,isec
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
    samples(ismpl)%wgt = 1d0
  enddo

  if( nserr.gt.0 ) call set_sample_errors()
  if( test_assigned ) then
    call count_training_test()
  else
    call set_training_test_with_ratio()
    test_assigned= .true.
  endif

  call read_samples()
  call read_ref_data()
  call get_base_energies()
  if( lswgt ) call set_sample_weights()

!.....Subtract energy and forces of other force-fields
!  call subtract_other_FF()
  
  call read_vars()
  allocate(gvar(nvars),dvar(nvars))

  select case (trim(cfmethod))
    case ('sd','SD')
      call sd_wrapper()
    case ('cg','CG')
      call cg_wrapper()
    case ('bfgs','BFGS','dfp','DFP')
      call qn_wrapper()
    case ('lbfgs','LBFGS','L-BFGS')
      call lbfgs_wrapper()
    case ('sa','SA')
      call sa_wrapper()
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
  if( nsmpl.lt.10000 ) then
    call write_force_relation('fin')
  endif
  call write_stats(niter)
  if(trim(cpena).eq.'lasso' .or. trim(cpena).eq.'glasso') &
       call write_eliminated_vars()
!!$  write(6,'(a,i4,3f15.3)') ' myid,tfunc,tgrad,tcom=' &
!!$       ,myid,tfunc,tgrad,tcomm
!!$  tmp= tfunc
!!$  call mpi_reduce(tmp,tfunc,1,mpi_double_precision,mpi_max &
!!$       ,0,mpi_world,ierr)
!!$  tmp= tgrad
!!$  call mpi_reduce(tmp,tgrad,1,mpi_double_precision,mpi_max &
!!$       ,0,mpi_world,ierr)
!!$  tmp= tcomm
!!$  call mpi_reduce(tmp,tcomm,1,mpi_double_precision,mpi_max &
!!$       ,0,mpi_world,ierr)
  if( myid.eq.0 ) then
    write(6,'(a,i10)') ' num of func calls=',nfunc
    write(6,'(a,i10)') ' num of grad calls=',ngrad
    write(6,'(a,f15.3,a)') ' time func =', tfunc,' sec'
    write(6,'(a,f15.3,a)') ' time grad =', tgrad,' sec'
    write(6,'(a,f15.3,a)') ' time comm =', tcomm,' sec'
    tmp = mpi_wtime() -time0
    ihour = int(tmp/3600)
    imin  = int((tmp-ihour*3600)/60)
    isec  = int(tmp -ihour*3600 -imin*60)
    write(6,'(a,f15.3,a,i3,"h",i2.2,"m",i2.2,"s")') &
         ' time      =', tmp, &
         ' sec  = ', ihour,imin,isec
  endif
  call mpi_finalize(ierr)

end program fitpot
!=======================================================================
subroutine write_initial_setting()
  use variables
  use minimize
  use random
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
  write(6,'(a)') ''
  write(6,'(2x,a25,2x,a)') 'penalty',trim(cpena)
  write(6,'(2x,a25,2x,es12.3)') 'penalty_weight',pwgt
  write(6,'(2x,a25,2x,a)') 'potential',trim(cpot)
  write(6,'(2x,a25,2x,l3)') 'gradient',lgrad
  write(6,'(2x,a25,2x,l3)') 'grad_scale',lgscale
  write(6,'(2x,a25,2x,es12.3)') 'gscale_factor',gscl
  write(6,'(2x,a25,2x,a)') 'normalize_input',trim(cnormalize)
  write(6,'(2x,a25,2x,es12.3)') 'freduce_threshold',fred
  write(6,'(2x,a25,2x,l3)') 'sample_weight',lswgt
  write(6,'(2x,a25,2x,es12.3)') 'sample_weight_erg',swerg
  write(6,'(2x,a25,2x,es12.3)') 'coeff_sequential',seqcoef
  write(6,'(2x,a25,2x,a)') 'line_minimization',trim(clinmin)
  write(6,'(a)') ''
  if( trim(cfmethod).eq.'sa' .or. trim(cfmethod).eq.'SA' ) then
    write(6,'(2x,a25,2x,es12.3)') 'sa_temperature',sa_temp0
    write(6,'(2x,a25,2x,es12.3)') 'sa_dxwidth',sa_xw0
    write(6,'(2x,a25,2x,es12.3)') 'random_seed',rseed
  endif
  write(6,'(a)') ''
  write(6,'(2x,a25,2x,i5)') 'individual_weight',nwgtindiv
  do i=1,nwgtindiv
    write(6,'(2x,a25,2x,f6.1)') trim(cwgtindiv(i)),wgtindiv(i)
  enddo
  write(6,'(a)') '------------------------------------------------'

end subroutine write_initial_setting
!=======================================================================
subroutine get_dir_list(ionum)
  use variables
  use parallel
  implicit none
  interface
    subroutine shuffle_dirlist(nsmpl,cdirlist,iclist)
      integer,intent(in):: nsmpl
      character(len=128),intent(inout):: cdirlist(nsmpl)
      integer,optional,intent(inout):: iclist(nsmpl)
    end subroutine shuffle_dirlist
  end interface
  integer,intent(in):: ionum
  integer:: is,ndat
  integer,external:: ndat_in_line
  logical:: lerror = .false.

  if( .not. allocated(cdirlist)) allocate(cdirlist(nsmpl))
  if( .not. allocated(iclist)) allocate(iclist(nsmpl))

  if( myid.eq.0 ) then
    lerror = .true.
    if( len(trim(csmplist)).lt.1 ) then
      print *,'sample list was created by command line...'
      call system('ls '//trim(cmaindir) &
           //' | grep "smpl_" > dir_list.txt')
      open(ionum,file='dir_list.txt',status='old')
    else
      print *,'sample list was given by input.'
      open(ionum,file=trim(csmplist),status='old')
    endif
    ndat = ndat_in_line(ionum,' ')
    if( ndat.eq.1 ) then
      do is=1,nsmpl
        read(ionum,*,end=998) cdirlist(is)
      enddo
      lerror = .false.
      call shuffle_dirlist(nsmpl,cdirlist)
    else if(ndat.eq.2 ) then
      print *,'training and test are determined by input, ',trim(csmplist)
      do is=1,nsmpl
        read(ionum,*,end=998) cdirlist(is),iclist(is)
      enddo
      lerror = .false.
      call shuffle_dirlist(nsmpl,cdirlist,iclist)
    else
      print *,'[Error] ndat should be 1 or 2, ndat = ',ndat
      call mpi_finalize(ierr)
      stop
    endif
    close(ionum)
  endif  ! myid.eq.0
998 continue
  call mpi_bcast(lerror,1,mpi_logical,0,mpi_world,ierr)
  call mpi_barrier(mpi_world,ierr)
  if( lerror ) goto 999
  call mpi_bcast(cdirlist,128*nsmpl,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(ndat,1,mpi_integer,0,mpi_world,ierr)
  if( ndat.eq.2 ) then
    call mpi_bcast(iclist,nsmpl,mpi_integer,0,mpi_world,ierr)
    test_assigned = .true.
  endif

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
subroutine set_training_test_with_ratio()
  use variables
  use parallel
  implicit none
  integer:: ismpl,n
  integer,allocatable,dimension(:):: icll

  allocate(icll(nsmpl))
  if( .not. allocated(iclist) ) allocate(iclist(nsmpl))
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
    print *,'set_training_test_with_ratio done.'
  endif
  deallocate(icll)
  return
end subroutine set_training_test_with_ratio
!=======================================================================
subroutine count_training_test()
!
!  Count number of training set and test set,
!  only if iclist is already set.
!
  use variables
  use parallel
  implicit none
  integer:: i

  myntrn = 0
  myntst = 0
  do i=isid0,isid1
    samples(i)%iclass = iclist(i)
    if( iclist(i).eq.1 ) myntrn= myntrn +1
    if( iclist(i).eq.2 ) myntst= myntst +1
  enddo
  
  if( myid.eq.0 ) then
    nsmpl_trn = 0
    nsmpl_tst = 0
    do i=1,nsmpl
      if( iclist(i).eq.1 ) nsmpl_trn = nsmpl_trn +1
      if( iclist(i).eq.2 ) nsmpl_tst = nsmpl_tst +1
    enddo
    print *,'nsmpl, training, test=',nsmpl,nsmpl_trn,nsmpl_tst
  endif
end subroutine count_training_test
!=======================================================================
subroutine read_samples()
  use variables
  use parallel
  implicit none
  integer:: is
  character*128:: cdir
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
       ,smpl%fref(3,natm), smpl%ifcal(natm),smpl%fabs(natm) &
       ,smpl%va(3,natm),smpl%strsi(3,3,natm) &
       ,smpl%eki(3,3,natm),smpl%epi(natm) &
       ,smpl%chg(natm),smpl%chi(natm))
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
  integer:: ismpl,i,is,jflag,natm,nfrc,nftot,nfrcg,nftotg
  integer:: imax,ifsmpl,nfsmplmax,nfrefdat,ifcal
  character(len=128):: cdir
  real(8):: erefminl,ftmp(3),fmax
  integer,external:: ndat_in_line

  jflag= 0
  erefminl= 0d0
  nftot= 0
  nfrc = 0
  do ismpl=isid0,isid1
    cdir=samples(ismpl)%cdirname
    open(13,file=trim(cmaindir)//'/'//trim(cdir) &
         //'/erg.ref',status='old')
    read(13,*) samples(ismpl)%eref
    close(13)
!.....Subtract (isolated) atomic energy from eref
    samples(ismpl)%naps(1:mspcs) = 0
    do i=1,samples(ismpl)%natm
      is= samples(ismpl)%tag(i)
      samples(ismpl)%naps(is) = samples(ismpl)%naps(is) +1
      samples(ismpl)%eref= samples(ismpl)%eref -eatom(is)
    enddo
    samples(ismpl)%ifcal(1:samples(ismpl)%natm)= 1
!!$    erefminl= min(erefminl,samples(ismpl)%eref/samples(ismpl)%natm)
!    write(6,*) 'ismpl,naps=',ismpl,samples(ismpl)%naps(1:mspcs)

    open(14,file=trim(cmaindir)//'/'//trim(cdir) &
         //'/frc.ref',status='old')
    read(14,*) natm
    if( natm.ne.samples(ismpl)%natm ) then
      print *,'Error: natm in sample is not same as samples(ismpl)%natm'
      print *,' myid,ismpl,natm,samples(ismpl)%natm=',myid,ismpl &
           ,natm,samples(ismpl)%natm
      jflag= jflag +1
    endif
    nfrefdat = ndat_in_line(14,' ')
    do i=1,natm
      nftot= nftot + 1
      if( nfrefdat.eq.3 ) then
        read(14,*) ftmp(1:3)
        ifcal = 1
      else if( nfrefdat.eq.4 ) then
!.....if frc.ref includes ifcal values after each force data,
!.....read 4 values from every line
        read(14,*) ftmp(1:3), ifcal
      endif
      samples(ismpl)%fref(1:3,i)= ftmp(1:3)
      samples(ismpl)%ifcal(i)= ifcal
      samples(ismpl)%fabs(i)= sqrt(ftmp(1)**2 +ftmp(2)**2 +ftmp(3)**2)
    enddo
    close(14)

!.....count nfcal
    samples(ismpl)%nfcal= 0
    do i=1,natm
      if( samples(ismpl)%ifcal(i).eq.1 ) then
        samples(ismpl)%nfcal = samples(ismpl)%nfcal +1
        nfrc = nfrc +1
      endif
    enddo

  enddo

  if( jflag.gt.0 ) then
    call mpi_finalize(ierr)
    stop
  endif

  nfrcg= 0
  nftotg= 0
  call mpi_reduce(nfrc,nfrcg,1,mpi_integer,mpi_sum,0 &
       ,mpi_world,ierr)
  call mpi_reduce(nftot,nftotg,1,mpi_integer,mpi_sum,0 &
       ,mpi_world,ierr)

  if(myid.eq.0) then
!    write(6,'(a,es12.4)') ' erefmin = ',erefmin
    if( lfmatch ) then
      write(6,'(a,i8)') ' number of forces to be used = ',nfrcg
      write(6,'(a,i8)') ' total number of forces      = ',nftotg
    endif
    print *,'read_ref_data done.'
  endif

end subroutine read_ref_data
!=======================================================================
subroutine get_base_energies()
  use variables
  use parallel
  implicit none
  integer:: ismpl,naps(mspcs),natm,ispcs,ielem
  real(8):: ebl(mspcs),erg
  
  ebl(1:mspcs) = 0d0
  do ismpl=isid0,isid1
    naps(1:mspcs) = samples(ismpl)%naps(1:mspcs)
    natm = samples(ismpl)%natm
!!$    write(6,*) ' ismpl,eref =',ismpl,samples(ismpl)%eref
!!$    write(6,*) ' ismpl,naps =',ismpl,naps(1:mspcs)
    ielem = 0
    do ispcs=1,mspcs
      if( naps(ispcs).eq.natm ) then
        ! this system is unary system
        erg = samples(ismpl)%eref/natm
        ielem = ispcs
        exit
      endif
    enddo
    if( ielem.ne.0 ) ebl(ielem) = min(ebl(ielem),erg)
  enddo

  ebase(1:mspcs) = 0d0
  call mpi_allreduce(ebl,ebase,mspcs,mpi_double_precision,mpi_min &
       ,mpi_world,ierr)

  if(myid.eq.0) then
    write(6,'(a)') ' base energies:'
    do ispcs=1,mspcs
      write(6,'(a,i3,es12.4)') '   is, ebase(is) =',ispcs,ebase(ispcs)
    enddo
  endif

end subroutine get_base_energies
!=======================================================================
subroutine read_vars()
  use variables
  use parallel
  use random
  implicit none
  integer:: i
  real(8):: rs0

  if( myid.eq.0 ) then
!!$    open(15,file=trim(cmaindir)//'/'//cparfile,status='old')
    open(15,file=trim(cparfile),status='old')
    read(15,*) nvars, rcut, rc3
  endif
  call mpi_bcast(nvars,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(rcut,1,mpi_double_precision,0,mpi_world,ierr)
  call mpi_bcast(rc3,1,mpi_double_precision,0,mpi_world,ierr)
  allocate(vars(nvars),vranges(2,nvars))
  if( myid.eq.0 ) then
    do i=1,nvars
      read(15,*) vars(i),vranges(1:2,i)
!    print *, vars(i),vranges(1:2,i)
    enddo
    if( cinitv.eq.'gauss' ) then
      rs0 = get_seed()
      call set_seed(vinitrs)
      do i=1,nvars
        vars(i) = vinitsgm*(polarbm()-vinitmu)
      enddo
      call set_seed(rs0)
      write(6,'(a)') ' params are shuffled to give normal distribution'
      write(6,'(a,2es10.2)') '   with mu and sgm =',vinitmu,vinitsgm
    else
      write(6,'(a)') ' params are read from file: '//cparfile
    endif
  endif
  call mpi_bcast(vars,nvars,mpi_double_precision,0,mpi_world,ierr)

  close(15)
end subroutine read_vars
!=======================================================================
subroutine write_vars(cadd)
  use variables
  use parallel
  use NNd, only: NN_standardize, NN_restore_standard
  implicit none
  character(len=*),intent(in):: cadd
  integer:: i
  character(len=128):: cfname

  call NN_restore_standard()

!!$  cfname= trim(cmaindir)//'/'//trim(cparfile)//'.'//trim(cadd)
  cfname= trim(cparfile)//'.'//trim(cadd)

  if( myid.eq.0 ) then
    open(15,file=trim(cfname),status='replace')
    write(15,'(i10,2es15.4)') nvars,rcut,rc3
    do i=1,nvars
      write(15,'(es23.14e3,2es12.4)') vars(i),vranges(1:2,i)
    enddo
    close(15)
!    print *, 'wrote '//trim(cfname)
  endif

  call NN_standardize()

end subroutine write_vars
!=======================================================================
subroutine qn_wrapper()
  use variables
  use NNd,only:NN_init,NN_func,NN_grad,NN_restore_standard,NN_analyze
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
!!$  call NN_restore_standard()

  return
end subroutine qn_wrapper
!=======================================================================
subroutine lbfgs_wrapper()
  use variables
  use NNd,only:NN_init,NN_func,NN_grad,NN_restore_standard,NN_analyze
  use parallel
  use minimize
  implicit none
  integer:: i,m
  real(8):: fval
  external:: write_stats

  !.....NN specific code hereafter
  call NN_init()
  call lbfgs(nvars,vars,fval,gvar,dvar,xtol,gtol,ftol,niter &
       ,iprint,iflag,myid,NN_func,NN_grad,cfmethod &
       ,niter_eval,write_stats)
  call NN_analyze("fin")
!!$  call NN_restore_standard()

  return
end subroutine lbfgs_wrapper
!=======================================================================
subroutine sd_wrapper()
!
!  Steepest descent minimization
!
  use variables
  use NNd,only:NN_init,NN_func,NN_grad
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
  use NNd,only:NN_init,NN_func,NN_grad,NN_restore_standard,NN_analyze
  use parallel
  use minimize
  implicit none
  integer:: i,m
  real(8):: fval
  external:: write_stats

  !.....NN specific code hereafter
  call NN_init()
  call cg(nvars,vars,fval,gvar,dvar,xtol,gtol,ftol,niter &
       ,iprint,iflag,myid,NN_func,NN_grad,cfmethod,niter_eval &
       ,write_stats)
  call NN_analyze("fin")
!!$  call NN_restore_standard()
  
  return
end subroutine cg_wrapper
!=======================================================================
subroutine sa_wrapper()
  use variables
  use NNd,only:NN_init,NN_func,NN_grad,NN_restore_standard,NN_analyze
  use parallel
  use minimize
  implicit none
  integer:: i,m
  real(8):: fval
  external:: write_stats

  !.....NN specific code hereafter
  call NN_init()
  call sa(nvars,vars,fval,xtol,gtol,ftol,niter &
       ,iprint,iflag,myid,NN_func,cfmethod &
       ,niter_eval,write_stats)
  call NN_analyze("fin")
!!$  if( cpena.eq.'lasso' .or. cpena.eq.'glasso' ) then
!!$    call NN_restore_standard()
!!$  endif

  return
end subroutine sa_wrapper
!=======================================================================
subroutine sgd()
!
! Stochastic gradient decent (SGD)
!
  use variables
  use NNd,only:NN_init,NN_fs,NN_gs,NN_func,NN_grad,NN_analyze &
       ,NN_restore_standard
  use parallel
  use minimize
  use random
  implicit none
  integer,parameter:: niter_time= 1
  real(8),parameter:: alpha0  = 1d-2
  real(8),parameter:: dalpha  = 0.001d0
  real(8),parameter:: gmmin   = 0.5d0
  real(8),parameter:: gmmax   = 0.95d0
  real(8),parameter:: dgmm    = 0.01d0
  real(8),parameter:: tiny  = 1d-8
!!$  real(8),parameter:: dalpha  = 0.d0
  real(8),allocatable:: g(:),u(:),gp(:),v(:),g2m(:),v2m(:)
  integer:: iter,istp,iv,i,nsize
  real(8):: gnorm,alpha,alpha1,gmax,vmax,f,gg,fp,gpnorm,gamma,vnorm
  real(8):: ftrn,ftst
  integer:: ismpl

  allocate(g(nvars),u(nvars),gp(nvars),v(nvars)&
       ,g2m(nvars),v2m(nvars))

  if( myid.eq.0 ) then
    write(6,'(a,a)') ' sgd_update: ',csgdupdate
    write(6,'(a,i6)') ' sgd_batch_size = ',nsgdbsize
  endif

  if( nsgdbsize.gt.maxmynsmpl ) then
    if(myid.eq.0) then
      write(6,'(a)') ' Error: nsgdbsize > maxmynsml'
      write(6,'(a,2i6)') '   nsgdbsize,maxmynsmpl = '&
         ,nsgdbsize,maxmynsmpl
      write(6,*) '  sgd_batch_size is too large.'
      write(6,*) '  You should decrease it much less than maxmynsmpl'
      write(6,*) '  Or you should use other minimization method.'
    endif
    call mpi_finalize(ierr)
    stop
  endif

  allocate(ismplsgd(nsgdbsize))

  alpha1= r0sgd
  v(1:nvars)= 0d0
  gamma = gmmin

  call NN_init()
  do iter=1,niter
    if(mod(iter,niter_eval).eq.0) then
      call NN_func(nvars,vars,ftrn,ftst)
      call NN_grad(nvars,vars,g)
      call penalty(cpena,pwgt,nvars,f,g,fp,gp,vars)
      g(1:nvars)= g(1:nvars) +gp(1:nvars)
      gnorm= sqrt(sprod(nvars,g,g))
      if( myid.eq.0 ) then
        write(6,'(a,i6,2es15.7,f10.3)') ' iter,f,gnorm,time=',iter,ftrn &
             ,gnorm ,mpi_wtime()-time0
        call write_vars('tmp')
      endif
      call write_stats(iter)
    else if(mod(iter,niter_time).eq.0) then
      if( myid.eq.0 ) then
        write(6,'(a,i6,f10.3)') ' iter,time=',iter,mpi_wtime()-time0
      endif
    endif

    nsize = 0
    call get_uniq_iarr(mynsmpl,nsgdbsize,ismplsgd)
    do i=1,nsgdbsize
      ismplsgd(i)= ismplsgd(i)+isid0-1
    enddo
    call NN_fs(nvars,vars,ftrn,ftst)
    call NN_gs(nvars,vars,g)
    call penalty(cpena,pwgt,nvars,ftrn,g,fp,gp,vars)
    gnorm= sqrt(sprod(nvars,g,g))
    gpnorm= sqrt(sprod(nvars,gp,gp))
    g(1:nvars)= g(1:nvars) +gp(1:nvars)
    u(1:nvars)= -g(1:nvars)
    if( csgdupdate.eq.'armijo' ) then
      alpha= r0sgd
      call armijo_search(nvars,vars,u,ftrn,ftst,g,alpha,iprint &
           ,iflag,myid,NN_fs)
      vars(1:nvars)=vars(1:nvars) +alpha*u(1:nvars)
      alpha1= alpha1*(1d0-dalpha)
    else if( csgdupdate.eq.'momentum' ) then
      v(1:nvars)= gamma*v(1:nvars) +alpha1*u(1:nvars)
      vars(1:nvars)= vars(1:nvars) +v(1:nvars)
      vnorm = sqrt(sprod(nvars,v,v))
      write(6,'(a,i5,3es12.3)') 'iter,gamma,vnorm=',iter,gamma,vnorm,gnorm
      gamma= min(gamma+dgmm,gmmax)
    else ! default: adadelta
      if(iter.eq.1) then
        alpha = r0sgd
        do i=1,nvars
          g2m(i)= gamma*g2m(i) +(1d0-gamma)*g(i)*g(i)
          v(i)= -alpha/sqrt(g2m(i)+tiny) *g(i)
        enddo
      else
        do i=1,nvars
          v2m(i)= gamma*v2m(i) +(1d0-gamma)*v(i)*v(i)
          g2m(i)= gamma*g2m(i) +(1d0-gamma)*g(i)*g(i)
          v(i)= -sqrt(v2m(i)+tiny)/sqrt(g2m(i)+tiny) *g(i)
        enddo
      endif
      vars(1:nvars)= vars(1:nvars) +v(1:nvars)
!!$      write(6,'(a,i5,10es11.3)') 'iter,v(:)=',iter,v(1:10)
    endif
  enddo

  call NN_func(nvars,vars,ftrn,ftst)
  call NN_grad(nvars,vars,g)
  gnorm= 0d0
  do iv=1,nvars
    gnorm= gnorm +g(iv)*g(iv)
  enddo

  call NN_analyze("fin")
!!$  call NN_restore_standard()

  deallocate(ismplsgd,g,u,gp,v,g2m,v2m)
end subroutine sgd
!=======================================================================
subroutine fs_wrapper()
  use variables
  use NNd,only:NN_init,NN_func,NN_grad,NN_restore_standard,NN_analyze
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
!!$  call NN_restore_standard()

  return
end subroutine fs_wrapper
!=======================================================================
subroutine gfs_wrapper()
  use variables
  use NNd,only:NN_init,NN_func,NN_grad,NN_restore_standard,NN_analyze
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
!!$  call NN_restore_standard()

  return
end subroutine gfs_wrapper
!=======================================================================
subroutine check_grad()
  use variables
  use NNd,only:NN_init,NN_func,NN_grad
  use parallel
  implicit none
  integer:: iv
  real(8):: ftrn0,ftst0,ftmp,dv,vmax,ftst,ftmp1,ftmp2
  real(8),allocatable:: ganal(:),gnumer(:),vars0(:)
  real(8),parameter:: dev  = 1d-5
  real(8),parameter:: tiny = 1d-6

  allocate(gnumer(nvars),ganal(nvars),vars0(nvars))
  call NN_init()
  call NN_func(nvars,vars,ftrn0,ftst0)
  call NN_grad(nvars,vars,ganal)

  vars0(1:nvars)= vars(1:nvars)
  vmax= 0d0
  do iv=1,nvars
    vmax= max(vmax,abs(vars0(iv)))
!!$    if( myid.eq.0) write(6,'(a,i6,es12.4)') ' iv,vars(iv)=',iv,vars0(iv)
  enddo
  dv= vmax *dev
  if( myid.eq.0 ) then
    print *,''
    print *,'deviation [dv] =',dv
  endif
  do iv=1,nvars
    vars(1:nvars)= vars0(1:nvars)
    dv = vars(iv)*dev
    vars(iv)= vars(iv) +dv/2
    call NN_func(nvars,vars,ftmp1,ftst)
    vars(1:nvars)= vars0(1:nvars)
    vars(iv)= vars(iv) -dv/2
    call NN_func(nvars,vars,ftmp2,ftst)
    gnumer(iv)= (ftmp1-ftmp2)/dv
  enddo

  if( myid.eq.0 ) then
    write(6,'(a)') '----------------- check_grad ------------------------'
    write(6,'(a)') '     #,          x,    analytical,'// &
         '     numerical,'// &
         '      error [%]'
    do iv=1,nvars
      write(6,'(i6,es12.4,2es15.4,f15.3)') iv,vars0(iv), &
           ganal(iv) ,gnumer(iv), &
           abs((ganal(iv)-gnumer(iv))/(gnumer(iv)+tiny))*100
    enddo
    write(6,'(a)') '-----------------------------------------------------'
    print *, 'check_grad done.'
  endif
end subroutine check_grad
!=======================================================================
subroutine test()
  use variables
  use NNd,only:NN_init,NN_func,NN_grad
  use parallel
  implicit none 
  integer:: iv
  real(8):: ftrn,ftst
  real(8),allocatable:: g(:)

  allocate(g(nvars))

  call NN_init()
  call NN_func(nvars,vars,ftrn,ftst)
  call NN_grad(nvars,vars,g)

  call write_stats(0)

  if( myid.eq.0 ) then
    print *,'func values (training,test) =',ftrn,ftst
    print *,'grad values (training):'
    do iv=1,nvars
      print *,'iv,grad(iv)=',iv,g(iv)
    enddo
    print *,'test done.'
  endif

  deallocate(g)

end subroutine test
!=======================================================================
subroutine write_energy_relation(cadd)
  use variables
  use parallel
  implicit none
  character(len=*),intent(in):: cadd
  character(len=128):: cfname
  
  integer:: ismpl,n
  logical,save:: l1st = .true.
  
  cfname='out.erg.'//trim(cadd)

  if( .not. allocated(erefl) ) allocate(erefl(nsmpl),erefg(nsmpl) &
       ,epotl(nsmpl),epotg(nsmpl),eerrl(nsmpl),eerrg(nsmpl)&
       ,swgtl(nsmpl),swgtg(nsmpl))

  if( l1st ) then
    erefl(1:nsmpl)= 0d0
    eerrl(1:nsmpl)= 0d0
    swgtl(1:nsmpl)= 0d0
    do ismpl=isid0,isid1
      erefl(ismpl)= samples(ismpl)%eref
      eerrl(ismpl)= samples(ismpl)%eerr
      swgtl(ismpl)= samples(ismpl)%wgt
    enddo
    erefg(1:nsmpl)= 0d0
    eerrg(1:nsmpl)= 0d0
    swgtg(1:nsmpl)= 0d0
    call mpi_reduce(erefl,erefg,nsmpl,mpi_double_precision,mpi_sum &
         ,0,mpi_world,ierr)
    call mpi_reduce(eerrl,eerrg,nsmpl,mpi_double_precision,mpi_sum &
         ,0,mpi_world,ierr)
    call mpi_reduce(swgtl,swgtg,nsmpl,mpi_double_precision,mpi_sum &
         ,0,mpi_world,ierr)
  endif

  epotl(1:nsmpl)= 0d0
  do ismpl=isid0,isid1
    epotl(ismpl)= samples(ismpl)%epot
  enddo
  epotg(1:nsmpl)= 0d0
  call mpi_reduce(epotl,epotg,nsmpl,mpi_double_precision,mpi_sum &
       ,0,mpi_world,ierr)

  if( myid.eq.0 ) then
    open(90,file=trim(cfname)//'.1',status='replace')
    open(91,file=trim(cfname)//'.2',status='replace')
    do ismpl=1,nsmpl
      erefg(ismpl)= erefg(ismpl)/nalist(ismpl)
      epotg(ismpl)= epotg(ismpl)/nalist(ismpl)
      if( iclist(ismpl).eq.1 ) then
        write(90,'(2es15.7,2x,a,3es15.6e3)') erefg(ismpl) &
             ,epotg(ismpl),trim(cdirlist(ismpl)) &
             ,abs(erefg(ismpl)-epotg(ismpl)) &
             ,eerrg(ismpl),swgtg(ismpl)
      else if( iclist(ismpl).eq.2 ) then
        write(91,'(2es15.7,2x,a,3es15.6e3)') erefg(ismpl) &
             ,epotg(ismpl),trim(cdirlist(ismpl)) &
             ,abs(erefg(ismpl)-epotg(ismpl)) &
             ,eerrg(ismpl),swgtg(ismpl)
!!$        write(91,'(2es15.7,2x,a)') erefg(ismpl)/nalist(ismpl) &
!!$             ,epotg(ismpl)/nalist(ismpl),cdirlist(ismpl)
      endif
    enddo
    close(90)
    close(91)
  endif
  l1st = .false.
  
end subroutine write_energy_relation
!=======================================================================
subroutine write_force_relation(cadd)
  use variables
  use parallel
  implicit none
  character(len=*),intent(in):: cadd
  character(len=128):: cfname

  integer:: ismpl,ia,ixyz,natm,nmax,nmaxl
  logical:: l1st = .true.
  
  cfname= 'out.frc.'//trim(cadd)

  nmaxl= 0
  do ismpl=1,nsmpl
    nmaxl= max(nmaxl,nalist(ismpl))
  enddo
  call mpi_allreduce(nmaxl,nmax,1,mpi_integer,mpi_max &
       ,mpi_world,ierr)

  if( .not. allocated(frefl) ) allocate(frefl(3,nmax,nsmpl)&
       ,frefg(3,nmax,nsmpl),fal(3,nmax,nsmpl),fag(3,nmax,nsmpl)&
       ,ferrl(nsmpl),ferrg(nsmpl))

  if( l1st ) then
    frefl(1:3,1:nmax,1:nsmpl)= 0d0
    ferrl(1:nsmpl) = 0d0
    do ismpl=isid0,isid1
      natm= samples(ismpl)%natm
      frefl(1:3,1:natm,ismpl)= samples(ismpl)%fref(1:3,1:natm)
      ferrl(ismpl) = samples(ismpl)%ferr
    enddo
    frefg(1:3,1:nmax,1:nsmpl)= 0d0
    ferrg(1:nsmpl) = 0d0
    call mpi_reduce(frefl,frefg,3*nmax*nsmpl,mpi_double_precision,mpi_sum &
         ,0,mpi_world,ierr)
    call mpi_reduce(ferrl,ferrg,nsmpl,mpi_double_precision,mpi_sum &
         ,0,mpi_world,ierr)
  endif

  fal(1:3,1:nmax,1:nsmpl)= 0d0
  do ismpl=isid0,isid1
    natm= samples(ismpl)%natm
    fal(1:3,1:natm,ismpl)= samples(ismpl)%fa(1:3,1:natm)
  enddo
  fag(1:3,1:nmax,1:nsmpl)= 0d0
  call mpi_reduce(fal,fag,3*nmax*nsmpl,mpi_double_precision,mpi_sum &
       ,0,mpi_world,ierr)

  if( myid.eq.0 ) then
    open(92,file=trim(cfname)//'.1',status='replace')
    open(93,file=trim(cfname)//'.2',status='replace')
    do ismpl=1,nsmpl
      if( iclist(ismpl).eq.1 ) then
        natm= nalist(ismpl)
        do ia=1,natm
          do ixyz=1,3
            write(92,'(2es15.7,2x,a,i6,i3,2es15.6e3)') frefg(ixyz,ia,ismpl) &
                 ,fag(ixyz,ia,ismpl) &
                 ,trim(cdirlist(ismpl)),ia,ixyz &
                 ,abs(frefg(ixyz,ia,ismpl)-fag(ixyz,ia,ismpl))&
                 ,ferrg(ismpl)
          enddo
        enddo
      else if( iclist(ismpl).eq.2 ) then
        natm= nalist(ismpl)
        do ia=1,natm
          do ixyz=1,3
            write(93,'(2es15.7,2x,a,i6,i3,2es15.6e3)') frefg(ixyz,ia,ismpl) &
                 ,fag(ixyz,ia,ismpl) &
                 ,trim(cdirlist(ismpl)),ia,ixyz &
                 ,abs(frefg(ixyz,ia,ismpl)-fag(ixyz,ia,ismpl))&
                 ,ferrg(ismpl)
          enddo
        enddo
      endif
    enddo
    close(92)
    close(93)
  endif
  l1st = .false.
end subroutine write_force_relation
!=======================================================================
subroutine write_stats(iter)
  use variables
  use parallel
  use NNd
  implicit none
  integer,intent(in):: iter
  integer:: ismpl,natm,ntrnl,ntstl,ia,l,ntrn,ntst,nfcal
  type(mdsys)::smpl
  real(8):: de,df
  real(8):: demaxl_trn,demax_trn,desuml_trn,desum_trn,rmse_trn
  real(8):: demaxl_tst,demax_tst,desuml_tst,desum_tst,rmse_tst
  real(8):: dfmaxl_trn,dfmax_trn,dfsuml_trn,dfsum_trn
  real(8):: dfmaxl_tst,dfmax_tst,dfsuml_tst,dfsum_tst
  real(8),save:: rmse_tst_best= 1d+30
  character:: cnum*5
  logical,save:: l1st = .true.

  if( l1st ) then
    if( myid.eq.0 ) then
      write(6,*) '# ENERGY: ITER, TIME, ' &
           //'RMSE(TRAINING), MAX(TRAINING), ' &
           //'RMSE(TEST), RMSE(TEST)'
      write(6,*) '# FORCE:  ITER, TIME, ' &
           //'RMSE(TRAINING), MAX(TRAINING), ' &
           //'RMSE(TEST), RMSE(TEST)'
    endif
  endif

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
  if( nsmpl_tst.ne.0 ) then
    rmse_tst= sqrt(desum_tst/nsmpl_tst)
  else
    rmse_tst= 0.d0
  endif
  if( myid.eq.0 ) then
!!$    write(6,'(a,2i6)') ' nsmpl_trn, nsmpl_tst = ',nsmpl_trn,nsmpl_tst
    write(6,'(a,i8,f15.2,4(1x,f12.7))') ' ENERGY: ' &
         ,iter,mpi_wtime()-time0 &
         ,rmse_trn,demax_trn,rmse_tst,demax_tst
!!$    if( rmse_tst < rmse_tst_best ) then
!!$      rmse_tst_best= rmse_tst
!!$      call write_vars('best')
!!$    endif
  endif
  write(cnum(1:5),'(i5.5)') iter
  call write_vars(cnum)

!.....force
  dfmaxl_trn= 0d0
  dfsuml_trn= 0d0
  dfmaxl_tst= 0d0
  dfsuml_tst= 0d0
  ntrnl= 0
  ntstl= 0
  do ismpl=isid0,isid1
    smpl= samples(ismpl)
    nfcal= smpl%nfcal
    if( nfcal.eq.0 ) cycle
    natm= smpl%natm
    if( smpl%iclass.eq.1 ) then
      do ia=1,natm
        if( smpl%ifcal(ia).eq.0 ) cycle
        do l=1,3
          df= abs(smpl%fa(l,ia)-smpl%fref(l,ia))
          dfmaxl_trn= max(dfmaxl_trn,df)
          dfsuml_trn=dfsuml_trn +df*df
          ntrnl=ntrnl +1
        enddo
      enddo
    else if( smpl%iclass.eq.2 ) then
      do ia=1,natm
        if( smpl%ifcal(ia).eq.0 ) cycle
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
  if( ntst.ne.0 ) then
    rmse_tst= sqrt(dfsum_tst/ntst)
  else
    rmse_tst= 0d0
  endif
  if( myid.eq.0 ) then
    write(6,'(a,i8,f15.2,4(1x,f12.7))') ' FORCE:  ' &
         ,iter,mpi_wtime()-time0 &
         ,rmse_trn,dfmax_trn,rmse_tst,dfmax_tst
!    call write_vars('tmp')
  endif
  
  l1st = .false.
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
  use NNd
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
  use random
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
  call mpi_bcast(cevaltype,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(cpot,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(cpena,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(clinmin,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(cnormalize,128,mpi_character,0,mpi_world,ierr)

  call mpi_bcast(xtol,1,mpi_double_precision,0,mpi_world,ierr)
  call mpi_bcast(ftol,1,mpi_double_precision,0,mpi_world,ierr)
  call mpi_bcast(gtol,1,mpi_double_precision,0,mpi_world,ierr)
  call mpi_bcast(eatom,maxnsp,mpi_double_precision,0,mpi_world,ierr)
  call mpi_bcast(gscl,1,mpi_double_precision,0,mpi_world,ierr)
  call mpi_bcast(fred,1,mpi_double_precision,0,mpi_world,ierr)
  call mpi_bcast(nfpsmpl,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(pwgt,1,mpi_double_precision,0,mpi_world,ierr)
  call mpi_bcast(ratio_test,1,mpi_double_precision,0,mpi_world,ierr)
  call mpi_bcast(rseed,1,mpi_double_precision,0,mpi_world,ierr)
  
  call mpi_bcast(lfmatch,1,mpi_logical,0,mpi_world,ierr)
  call mpi_bcast(lgrad,1,mpi_logical,0,mpi_world,ierr)
  call mpi_bcast(lgscale,1,mpi_logical,0,mpi_world,ierr)

  call mpi_bcast(lswgt,1,mpi_logical,0,mpi_world,ierr)
  call mpi_bcast(swerg,1,mpi_double_precision,0,mpi_world,ierr)

  call mpi_bcast(sa_temp0,1,mpi_double_precision,0,mpi_world,ierr)
  call mpi_bcast(sa_xw0,1,mpi_double_precision,0,mpi_world,ierr)
!.....sgd
  call mpi_bcast(csgdupdate,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(r0sgd,1,mpi_double_precision,0,mpi_world,ierr)
  call mpi_bcast(nsgdbsize,1,mpi_integer,0,mpi_world,ierr)
!.....initialize parameters
  call mpi_bcast(cinitv,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(vinitsgm,1,mpi_double_precision,0,mpi_world,ierr)
  call mpi_bcast(vinitmu,1,mpi_double_precision,0,mpi_world,ierr)
  call mpi_bcast(numtol,1,mpi_integer,0,mpi_world,ierr)
!.....CG
  call mpi_bcast(icgbtype,1,mpi_integer,0,mpi_world,ierr)
!.....L-BFGS
  call mpi_bcast(mstore,1,mpi_integer,0,mpi_world,ierr)
!.....Armijo
  call mpi_bcast(armijo_xi,1,mpi_double_precision,0,mpi_world,ierr)
  call mpi_bcast(armijo_tau,1,mpi_double_precision,0,mpi_world,ierr)
  call mpi_bcast(armijo_maxiter,1,mpi_integer,0,mpi_world,ierr)

  call mpi_bcast(nwgtindiv,1,mpi_integer,0,mpi_world,ierr)
  if( myid.gt.0 ) then
    allocate(cwgtindiv(nwgtindiv),wgtindiv(nwgtindiv))
  endif
  call mpi_bcast(cwgtindiv,128*nwgtindiv,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(wgtindiv,nwgtindiv,mpi_double_precision,0,mpi_world,ierr)

  call mpi_bcast(nserr,1,mpi_integer,0,mpi_world,ierr)
  if( myid.gt.0 ) then
    allocate(cserr(nserr),seerr(nserr),sferr(nserr))
  endif
  call mpi_bcast(cserr,128*nserr,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(seerr,nserr,mpi_double_precision,0,mpi_world,ierr)
  call mpi_bcast(sferr,nserr,mpi_double_precision,0,mpi_world,ierr)
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
subroutine shuffle_dirlist(nsmpl,cdirlist,iclist)
!
!  Randomize the order of the cdirlist
!
  use random
  implicit none
  integer,intent(in):: nsmpl
  character(len=128),intent(inout):: cdirlist(nsmpl)
  integer,optional,intent(inout):: iclist(nsmpl)

  integer:: i,j,k,n
  character(len=128),allocatable:: cdltmp(:)
  integer,allocatable:: icltmp(:)

  allocate(cdltmp(nsmpl))
  if( .not. present(iclist) ) then
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
!!$  do i=1,nsmpl
!!$    print *,i,trim(cdirlist(i))
!!$  enddo
  else
    allocate(icltmp(nsmpl))
    cdltmp(1:nsmpl)= cdirlist(1:nsmpl)
    icltmp(1:nsmpl)= iclist(1:nsmpl)
    n=nsmpl
    do i=1,nsmpl
      j= n*urnd()+1
      cdirlist(i)= cdltmp(j)
      iclist(i)= icltmp(j)
      do k=j,n-1
        cdltmp(k)= cdltmp(k+1)
        icltmp(k)= icltmp(k+1)
      enddo
      n= n -1
    enddo
    deallocate(icltmp)
  endif
  deallocate(cdltmp)
end subroutine shuffle_dirlist
!=======================================================================
subroutine set_individual_weight()
  use variables
  use parallel
  implicit none 
  integer:: ismpl,iwgt,idx

  do ismpl=isid0,isid1
    do iwgt= 1,nwgtindiv
      idx= index(samples(ismpl)%cdirname,trim(cwgtindiv(iwgt)))
      if( idx.ne.0 ) then
        samples(ismpl)%wgt= wgtindiv(iwgt)
      endif
    enddo
  enddo
  if( myid.eq.0 ) then
    write(6,'(a)') ' set_individual_weight done'
  endif
end subroutine set_individual_weight
!=======================================================================
subroutine set_sample_errors()
  use variables
  use parallel
  implicit none 
  integer:: ismpl,iserr,idx

  do ismpl=isid0,isid1
    do iserr= 1,nserr
      idx= index(samples(ismpl)%cdirname,trim(cserr(iserr)))
      if( idx.ne.0 ) then
        samples(ismpl)%eerr= seerr(iserr)
        samples(ismpl)%ferr= sferr(iserr)
      endif
    enddo
  enddo
  if( myid.eq.0 ) then
    write(6,'(a)') ' set_samples_errors done'
  endif
end subroutine set_sample_errors
!=======================================================================
subroutine set_sample_weights()
  use variables
  use parallel
  implicit none
  integer:: naps(mspcs),ismpl,is
  real(8):: erg
  
  do ismpl=isid0,isid1
    naps(1:mspcs) = samples(ismpl)%naps(1:mspcs)
    erg = samples(ismpl)%eref
    do is=1,mspcs
      erg = erg -naps(is)*ebase(is)
    enddo
    erg = erg/samples(ismpl)%natm
    samples(ismpl)%wgt = min(exp(-erg/swerg),1d0)
  enddo
  if( myid.eq.0 ) then
    write(6,'(a)') ' set_sample_weights done'
  endif
end subroutine set_sample_weights
!=======================================================================
subroutine get_uniq_iarr(n,m,iarr)
!
! Create an array with length m which includes a unique set of integers
! randomly chosen from from 1 to n.
!
  use random
  integer,intent(in):: n,m
  integer,intent(out):: iarr(m)
  integer:: jarr(n)
  integer:: i,j,k,l
  
  do i=1,n
    jarr(i) = i
  enddo
  l=n
  do i=1,m
    j = l*urnd()+1
    iarr(i)= jarr(j)
    do k=j,l-1
      jarr(k)= jarr(k+1)
    enddo
    l= l-1
  enddo
  return
end subroutine get_uniq_iarr
!=======================================================================
subroutine subtract_other_FF()
!
!  Subtract energies and forces from other force-fields.
!  This uses force-fields implemented in pmd and the NN potential made
!  should also be used with those force-fields, of course.
!  This routine should be called for each force-field specified, so
!  it could be called several times if several force-fields are taken
!  into account.
!
  use variables
  use parallel
  implicit none

  integer:: ismpl
  type(mdsys):: smpl

  do ismpl=isid0,isid1
    smpl = samples(ismpl)
    call run_pmd(smpl)
  enddo
  
  
end subroutine subtract_other_FF
!=======================================================================
subroutine run_pmd(smpl)
!
!  Run pmd and get energy and forces of the system.
!
  use variables
  use parallel
  implicit none
  type(mdsys),intent(inout):: smpl

  logical,save:: l1st = .true.

  integer:: i,maxstp,nerg,npmd,ifpmd,ifdmp,minstp,n_conv,ifsort, &
       ifcoulomb,nismax,nstps_done,ntdst,nx,ny,nz,numff
  real(8):: am(9),dt,rc,rbuf,dmp,tinit,tfin,ttgt(9),trlx,stgt(3,3),&
       ptgt,srlx,stbeta,strfin,fmv(3,0:9),ptnsr(3,3),epot,ekin,eps_conv
  logical:: ltdst,lstrs,lcellfix(3,3)
  character:: ciofmt*6,ctctl*20,cpctl*20,czload_type*5
  character(len=20):: cffs(1)

  if( l1st ) then
!.....Create MPI COMM for pmd only for the 1st time
    call create_mpi_comm_pmd()
  endif

!.....Every time allocate total arrays according to the num of atoms
!     in the givin sample system.
  maxstp = 0
  nismax = 9
  nerg = 1
  npmd = 1
  am(1:9) = 1d0  ! Since no dynamics, no need of mass
  dt = 5d0
  ciofmt = 'ascii'
  ifpmd = 0
  numff = 1
  cffs(1) = 'NN'
  rc = 5.5d0
  rbuf = 0.2d0
  ifdmp = 0  ! no damping as well
  dmp = 0.99d0
  minstp = 0
  tinit = 0d0
  tfin = 0d0
  ctctl = 'none'
  ttgt(1:9) = 300d0
  trlx = 100d0
  ltdst = .false.
  ntdst = 1
  lstrs = .false.
  cpctl = 'none'
  stgt(1:3,1:3) = 0d0
  ptgt = 0d0
  srlx = 100d0
  stbeta = 1d-1
  strfin = 0d0
  fmv(1:3,0) = (/ 0d0, 0d0, 0d0 /)
  fmv(1:3,1:9) = 1d0
  ptnsr(1:3,1:3) = 0d0
  epot = 0d0
  ekin = 0d0
  n_conv = 1
  czload_type = 'no'
  eps_conv = 1d-3
  ifsort = 1
  iprint = 0
  ifcoulomb = 0
  lcellfix(1:3,1:3) = .false.
  nx = 1
  ny = 1
  nz = 1
  
!.....Run one-shot force calculation to get an energy and forces
  call pmd_core(smpl%h0,smpl%h,smpl%natm,smpl%tag,smpl%ra &
       ,smpl%va,smpl%fa,smpl%strsi,smpl%eki,smpl%epi &
       ,smpl%chg,smpl%chi,maxstp,nerg,npmd &
       ,myid_pmd,mpi_comm_pmd,nnode_pmd,nx,ny,nz &
       ,nismax,am,dt,ciofmt,ifpmd,numff,cffs,rc,rbuf,ifdmp,dmp,minstp &
       ,tinit,tfin,ctctl,ttgt,trlx,ltdst,ntdst,cpctl,stgt,ptgt &
       ,srlx,stbeta,strfin,lstrs,lcellfix &
       ,fmv,ptnsr,epot,ekin,n_conv,ifcoulomb &
       ,czload_type,eps_conv,ifsort,iprint,nstps_done)

!.....Subtract energy and forces from eref and fref, respectively
  smpl%eref = smpl%eref -epot
  do i=1,smpl%natm
    smpl%fref(1:3,i) = smpl%fref(1:3,i) -smpl%fa(1:3,i)
  enddo
  
end subroutine run_pmd
!=======================================================================
subroutine create_mpi_comm_pmd()
!
!  Create MPI COMM for pmd.
!  To create a MPI COMM on each node, first create a MPI GROUP for world,
!  then create a MPI GROUP for this node, and then create the MPI COMM
!  from the GROUP for this node. This is how to create sub communicator
!  in the MPI.
!
  use parallel
  implicit none

  integer:: n = 1
  integer:: iranks(1)
  integer:: mpi_group_world,mpi_group_pmd

  call mpi_comm_group(mpi_world, mpi_group_world,ierr)

  iranks(1) = myid
  call mpi_group_incl(mpi_group_world, 1, iranks, mpi_group_pmd,ierr)
  call mpi_comm_create_group(mpi_world, mpi_group_pmd, 0, mpi_comm_pmd,ierr)

  call mpi_comm_size(mpi_comm_pmd,nnode_pmd,ierr)
  call mpi_comm_rank(mpi_comm_pmd,myid_pmd,ierr)
  
  call mpi_group_free(mpi_group_world,ierr)
  call mpi_group_free(mpi_group_pmd,ierr)
  if( myid.eq.0 ) then
    write(6,'(a)') ''
    write(6,'(a)') 'MPI_COMM_PMD were created for pmd calculations.'
  endif
  
end subroutine create_mpi_comm_pmd
!-----------------------------------------------------------------------
! Local Variables:
! compile-command: "make fitpot"
! End:
