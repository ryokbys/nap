program fitpot
!-----------------------------------------------------------------------
!                     Last modified: <2017-10-26 15:21:35 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
  use variables
  use parallel
  use minimize
  use version
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
    write(6,'(a)') '========================================================================'
    write(6,'(a)') ' FITPOT: A program for FITting interatomic POTential parameters'
    write(6,*) ''
    call write_version()
    call write_authors()
    write(6,'(a)') '========================================================================'
    write(6,*) ''
    call time_stamp(' Job started')
    write(6,*) ''
    write(6,'(a,i6)') ' Number of processes in MPI = ',nnode
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
  call set_max_num_atoms()

  call read_vars()
  allocate(gvar(nvars),dvar(nvars))

!.....Subtract atomic energy
  if( trim(cpot).ne.'vcMorse' ) then
    if( len(trim(crefstrct)).gt.5 ) then
      call subtract_ref_struct_energy()
    else
      call subtract_atomic_energy()
    endif
  endif
  
  if( nswgt.gt.0 ) call set_sample_weights()

!.....Subtract energy and forces of other force-fields
  if( nsubff.gt.0 ) then
    call subtract_FF()
  endif

!.....Set cffs only for pmd calculation
  nff = 1
  allocate(cffs(nff))
  cffs(1) = trim(cpot)
  
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
    case ('ga','GA')
      call ga_wrapper()
    case ('de','DE')
      call de_wrapper()
    case ('pso','PSO')
      call pso_wrapper()
    case ('md','metadynamics')
      call md_wrapper()
    case ('random_search','random')
      call random_search_wrapper()
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
      print *,'Minimization stopped because of convergence or something.'
      print *,'  iflag=',iflag
    else if( iflag/100.ne.0 ) then
      print *,'Minimization stopped because of convergence or something.'
      print *,'  iflag=',iflag
    endif
  endif

!!$  if( myid.eq.0 ) write(6,'(a,100f7.3)')  ' vars beofre stats: ',vars(1:nvars)
  call write_stats(niter)

!!$  call write_energy_relation('subtracted')
!!$  if( nsmpl.lt.nsmpl_outfrc ) then
!!$    call write_force_relation('subtracted')
!!$  endif


!!$!.....Restore subtracted energies and forces to get original reference values
!!$  if( nsubff.gt.0 ) then
!!$    call restore_FF()
!!$  endif
  
  call write_vars('fin')
  call write_energy_relation('fin')
  if( nsmpl.lt.nsmpl_outfrc ) then
    call write_force_relation('fin')
  endif
  call write_stress_relation('fin')

  if(trim(cpena).eq.'lasso' .or. trim(cpena).eq.'glasso') &
       call write_eliminated_vars()
!!$  write(6,'(a,i4,3f15.3)') ' myid,tfunc,tgrad,tcom=' &
!!$       ,myid,tfunc,tgrad,tcomm
!!$  tmp= tfunc
!!$  call mpi_reduce(tmp,tfunc,1,mpi_real8,mpi_max &
!!$       ,0,mpi_world,ierr)
!!$  tmp= tgrad
!!$  call mpi_reduce(tmp,tgrad,1,mpi_real8,mpi_max &
!!$       ,0,mpi_world,ierr)
!!$  tmp= tcomm
!!$  call mpi_reduce(tmp,tcomm,1,mpi_real8,mpi_max &
!!$       ,0,mpi_world,ierr)
  if( myid.eq.0 ) then
    write(6,'(a,i0)') ' Number of func calls = ',nfunc
    write(6,'(a,i0)') ' Number of grad calls = ',ngrad
    write(6,'(a,f15.3,a)') ' Time func = ', tfunc,' sec'
    write(6,'(a,f15.3,a)') ' Time grad = ', tgrad,' sec'
    write(6,'(a,f15.3,a)') ' Time comm = ', tcomm,' sec'
    tmp = mpi_wtime() -time0
    ihour = int(tmp/3600)
    imin  = int((tmp-ihour*3600)/60)
    isec  = int(tmp -ihour*3600 -imin*60)
    write(6,'(a,f15.3,a,i3,"h",i2.2,"m",i2.2,"s")') &
         ' Time      = ', tmp, &
         ' sec  = ', ihour,imin,isec
    write(6,*) ''
    call time_stamp(' Job finished')
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

  write(6,*) ''
  write(6,'(a)') '------------------------------------------------------------------------'
  write(6,'(a)') '                          Initial setting                               '
  write(6,'(a)') '------------------------------------------------------------------------'
  write(6,'(2x,a25,2x,i5)') 'num_samples',nsmpl
  write(6,'(2x,a25,2x,i10)') 'num_iteration',niter
  write(6,'(2x,a25,2x,a)') 'fitting_method',trim(cfmethod)
  write(6,'(2x,a25,2x,a)') 'main_directory',trim(cmaindir)
  write(6,'(2x,a25,2x,a)') 'param_file',trim(cparfile)
!!$  write(6,'(2x,a25,2x,a)') 'run_mode',trim(crunmode)
  write(6,'(2x,a25,2x,es12.3)') 'xtol',xtol
  write(6,'(2x,a25,2x,es12.3)') 'ftol',ftol
  write(6,'(2x,a25,2x,es12.3)') 'gtol',gtol
  if( len(trim(crefstrct)).gt.5 ) then
    write(6,'(2x,a25,2x,a)') 'reference_structure',trim(crefstrct)
  else
    do i=1,maxnsp
      write(6,'(2x,a25,2x,i2,es15.7)') 'atom_energy',i,eatom(i)
    enddo
  endif
  write(6,'(2x,a25,2x,l3)') 'energy_match',lematch
  write(6,'(2x,a25,2x,l3)') 'force_match',lfmatch
  write(6,'(2x,a25,2x,l3)') 'stress_match',lsmatch
  write(6,'(a)') ''
  write(6,'(2x,a25,2x,a)') 'potential',trim(cpot)
  write(6,'(2x,a25,10(2x,a))') 'subtract_force_field',(trim(csubffs(i)),i=1,nsubff)
  write(6,'(a)') ''
  write(6,'(2x,a25,2x,a)') 'penalty',trim(cpena)
  write(6,'(2x,a25,2x,es12.3)') 'penalty_weight',pwgt
  write(6,'(2x,a25,2x,l3)') 'gradient',lgrad
!!$  write(6,'(2x,a25,2x,l3)') 'grad_scale',lgscale
!!$  write(6,'(2x,a25,2x,es12.3)') 'gscale_factor',gscl
  write(6,'(2x,a25,2x,a)') 'normalize_input',trim(cnormalize)
  if( lfmatch ) then
    write(6,'(2x,a25,2x,f0.2)') 'force_limit',force_limit
  endif

  if( nswgt.gt.0 ) then
    write(6,'(2x,a25,2x,i5)') 'sample_weight',nswgt
    do i=1,nswgt
      write(6,'(4x,a23,2x,f8.4,2x,f8.1)') trim(cswgt(i)), swerg0(i), swdenom(i)
    enddo
  endif
  
  write(6,'(2x,a25,2x,es12.3)') 'coeff_sequential',seqcoef
  write(6,'(2x,a25,2x,a)') 'line_minimization',trim(clinmin)
  write(6,'(a)') ''
  write(6,'(2x,a25,2x,i4)') 'sample_error',nserr
  do i=1,nserr
    write(6,'(4x,a23,3(1x,f8.4))') trim(cserr(i)), seerr(i), sferr(i), sserr(i)
  enddo
  write(6,'(a)') ''
  if( trim(cfmethod).eq.'sa' .or. trim(cfmethod).eq.'SA' ) then
    write(6,'(2x,a25,2x,es12.3)') 'sa_temperature',sa_temp0
    write(6,'(2x,a25,2x,es12.3)') 'sa_dxwidth',sa_xw0
    write(6,'(2x,a25,2x,es12.3)') 'random_seed',rseed
  else if( trim(cfmethod).eq.'ga' .or. trim(cfmethod).eq.'GA' ) then
    write(6,'(2x,a25,2x,i3)') 'ga_num_bits',ga_nbits
    write(6,'(2x,a25,2x,i4)') 'ga_num_individuals',ga_nindivs
    write(6,'(2x,a25,2x,i4)') 'ga_num_offsprings',ga_noffsp
    write(6,'(2x,a25,2x,a)') 'ga_fitness',trim(ga_fitness)
    write(6,'(2x,a25,2x,f8.4)') 'ga_mutation_rate',ga_rate_mutate
    write(6,'(2x,a25,2x,es12.3)') 'random_seed',rseed
  else if( trim(cfmethod).eq.'de' .or. trim(cfmethod).eq.'DE' ) then
    write(6,'(2x,a25,2x,i4)') 'de_num_individuals',de_nindivs
    write(6,'(2x,a25,2x,a)') 'de_fitness',trim(de_fitness)
    write(6,'(2x,a25,2x,f8.4)') 'de_fraction',de_frac
    write(6,'(2x,a25,2x,f8.4)') 'de_lambda',de_lambda
    write(6,'(2x,a25,2x,f8.4)') 'de_crossover_rate',de_cross_rate
    write(6,'(2x,a25,2x,f8.4)') 'de_wmin',de_wmin
    write(6,'(2x,a25,2x,f8.4)') 'de_wmax',de_wmax
    write(6,'(2x,a25,2x,es12.3)') 'random_seed',rseed
  else if( trim(cfmethod).eq.'pso' .or. trim(cfmethod).eq.'PSO' ) then
    write(6,'(2x,a25,2x,i4)') 'pso_num_individuals',pso_nindivs
    write(6,'(2x,a25,2x,f8.4)') 'pso_w',pso_w
    write(6,'(2x,a25,2x,f8.4)') 'pso_c1',pso_c1
    write(6,'(2x,a25,2x,f8.4)') 'pso_c2',pso_c2
  endif
!!$  write(6,'(a)') ''
!!$  write(6,'(2x,a25,2x,i5)') 'individual_weight',nwgtindiv
!!$  do i=1,nwgtindiv
!!$    write(6,'(2x,a25,2x,f6.1)') trim(cwgtindiv(i)),wgtindiv(i)
!!$  enddo
  write(6,'(a)') '------------------------------------------------------------------------'

end subroutine write_initial_setting
!=======================================================================
subroutine get_dir_list(ionum)
  use variables
  use parallel
  use fp_common,only: ndat_in_line
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
  logical:: lerror = .false.

  if( .not. allocated(cdirlist)) allocate(cdirlist(nsmpl))
  if( .not. allocated(iclist)) allocate(iclist(nsmpl))

  if( myid.eq.0 ) then
    lerror = .true.
    if( len(trim(csmplist)).lt.1 ) then
      print '(/,a)',' Sample list was created by performing the following'
      print *,'  $ ls '//trim(cmaindir) &
           //' | grep "smpl_" > dir_list.txt'
      call system('ls '//trim(cmaindir) &
           //' | grep "smpl_" > dir_list.txt')
      open(ionum,file='dir_list.txt',status='old')
    else
      print *,'Sample list was given by input.'
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
  
  if(myid.eq.0) print*,'Finished get_dir_list.'
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
  
  if( myid.eq.0 ) write(6,'(/,a)') 'Finished reading samples.'
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
  read(ionum,*) smpl%h(1:3,1)
  read(ionum,*) smpl%h(1:3,2)
  read(ionum,*) smpl%h(1:3,3)
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
       ,smpl%chg(natm),smpl%chi(natm),smpl%fsub(3,natm) &
       ,smpl%symbols(natm),smpl%eatm(natm))
  smpl%chg(1:natm) = 0d0
  smpl%esub= 0d0
  smpl%fsub(1:3,1:natm)= 0d0
  smpl%ssub(1:3,1:3) = 0d0
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
  use fp_common,only: ndat_in_line
  implicit none 
  integer:: ismpl,i,is,jflag,natm,nfrc,nftot,nfrcg,nftotg
  integer:: imax,ifsmpl,nfsmplmax,nfrefdat,ifcal
  character(len=128):: cdir
  real(8):: erefminl,ftmp(3),fmax,ptnsr(3,3)

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
!.....Count numbers of each species
    samples(ismpl)%naps(1:mspcs) = 0
    do i=1,samples(ismpl)%natm
      is= int(samples(ismpl)%tag(i))
      samples(ismpl)%naps(is) = samples(ismpl)%naps(is) +1
!!$      samples(ismpl)%eref= samples(ismpl)%eref  !-eatom(is)
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
!!$    print *,'ismpl,cdirname =',ismpl,trim(cdir)
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
      samples(ismpl)%fabs(i)= sqrt(ftmp(1)**2 +ftmp(2)**2 +ftmp(3)**2)
      if( samples(ismpl)%fabs(i).gt.force_limit ) ifcal = 0
      samples(ismpl)%ifcal(i)= ifcal
!!$      write(6,'(a,2i5,3es12.4)') 'ismpl,i,samples(ismpl)%fref(1:3,i) = ',&
!!$           ismpl,i,samples(ismpl)%fref(1:3,i)
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

!.....Read stress tensor data
    ptnsr(1:3,1:3) = 0d0
    open(15,file=trim(cmaindir)//'/'//trim(cdir) &
         //'/strs.ref',status='old')
    read(15,*) ptnsr(1,1),ptnsr(2,2),ptnsr(3,3), &
         ptnsr(2,3),ptnsr(1,3),ptnsr(1,2)
    close(15)
    ptnsr(2,1) = ptnsr(1,2)
    ptnsr(3,1) = ptnsr(1,3)
    ptnsr(3,2) = ptnsr(2,3)
    samples(ismpl)%sref(1:3,1:3) = ptnsr(1:3,1:3)
    
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
      write(6,'(a,i8)') ' Number of forces to be used = ',nfrcg
      write(6,'(a,i8)') ' Total number of forces      = ',nftotg
    endif
    print *,'Finished read_ref_data.'
  endif

end subroutine read_ref_data
!=======================================================================
subroutine get_base_energies()
!
!  Compute base energies for elements from unary systems.
!
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
  call mpi_allreduce(ebl,ebase,mspcs,mpi_real8,mpi_min &
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
  character(len=128):: fname

  if( trim(cpot).eq.'NN' ) then
    fname = cparfile
  else
    fname = 'in.vars.fitpot'
  endif

  if( myid.eq.0 ) then
    open(15,file=trim(fname),status='old')
    read(15,*) nvars, rcut, rc3
  endif
  call mpi_bcast(nvars,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(rcut,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(rc3,1,mpi_real8,0,mpi_world,ierr)
  allocate(vars(nvars),vranges(2,nvars))
  if( myid.eq.0 ) then
    do i=1,nvars
      read(15,*) vars(i),vranges(1:2,i)
!    print *, vars(i),vranges(1:2,i)
    enddo
    if( cinitv.eq.'gaussian' ) then
      rs0 = get_seed()
      call set_seed(vinitrs)
      do i=1,nvars
        vars(i) = vinitsgm*(polarbm()-vinitmu)
      enddo
      call set_seed(rs0)
      write(6,'(a)') ' Potential parameters are shuffled'&
           //' to give normal distribution'
      write(6,'(a,2es10.2)') '   with mu and sgm =',vinitmu,vinitsgm
    else
      write(6,'(a)') ' Potential parameters are read from file: '//trim(fname)
    endif
    close(15)
  endif
  call mpi_bcast(vars,nvars,mpi_real8,0,mpi_world,ierr)

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

  if( trim(cpot).eq.'NN' .and. .not. &
    (trim(cfmethod).eq.'sa' .or. trim(cfmethod).eq.'SA' .or. &
     trim(cfmethod).eq.'ga' .or. trim(cfmethod).eq.'GA' .or. &
     trim(cfmethod).eq.'de' .or. trim(cfmethod).eq.'DE' .or. &
     trim(cfmethod).eq.'pso' .or. trim(cfmethod).eq.'PSO') ) then
    call NN_restore_standard()
  endif

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

  if( trim(cpot).eq.'NN' .and. .not. &
    (trim(cfmethod).eq.'sa' .or. trim(cfmethod).eq.'SA' .or. &
     trim(cfmethod).eq.'ga' .or. trim(cfmethod).eq.'GA' .or. &
     trim(cfmethod).eq.'de' .or. trim(cfmethod).eq.'DE' .or. &
     trim(cfmethod).eq.'pso' .or. trim(cfmethod).eq.'PSO') ) then
    call NN_standardize()
  endif

end subroutine write_vars
!=======================================================================
subroutine qn_wrapper()
  use variables
  use NNd,only:NN_init,NN_func,NN_grad,NN_restore_standard,NN_analyze
  use parallel
  use minimize
  use fp_common,only: func_w_pmd, grad_w_pmd
  implicit none
  integer:: i,m
  real(8):: fval
  external:: write_stats

  if( trim(cpot).eq.'NN' ) then
!.....NN specific code hereafter
    call NN_init()
    call qn(nvars,vars,fval,gvar,dvar,xtol,gtol,ftol,niter &
         ,iprint,iflag,myid,NN_func,NN_grad,cfmethod &
         ,niter_eval,write_stats)
    call NN_analyze("fin")
    
  else if( trim(cpot).eq.'vcMorse' ) then
    call qn(nvars,vars,fval,gvar,dvar,xtol,gtol,ftol,niter &
         ,iprint,iflag,myid,func_w_pmd,grad_w_pmd,cfmethod &
         ,niter_eval,write_stats)
  endif

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
  use parallel
  use minimize
  use fp_common,only: func_w_pmd, grad_w_pmd
  implicit none
  integer:: i,m
  real(8):: fval
  external:: write_stats

  if( trim(cpot).eq.'vcMorse' .or. trim(cpot).eq.'Morse' .or. &
       trim(cpot).eq.'EAM' .or. trim(cpot).eq.'NN' ) then
    call sa(nvars,vars,fval,vranges,xtol,gtol,ftol,niter &
         ,iprint,iflag,myid,func_w_pmd,cfmethod &
         ,niter_eval,write_stats)
  else
    if(myid.eq.0) print *,'Simulated annealing is not available for '//&
         trim(cpot)
  endif

  return
end subroutine sa_wrapper
!=======================================================================
subroutine md_wrapper()
  use variables
  use parallel
  use minimize
  use fp_common,only: func_w_pmd, grad_w_pmd
  implicit none
  integer:: i,m
  real(8):: fval
  external:: write_stats

  if( trim(cpot).eq.'vcMorse' .or. trim(cpot).eq.'Morse' .or. &
       trim(cpot).eq.'EAM' .or. trim(cpot).eq.'NN' ) then
    call metadynamics(nvars,vars,fval,vranges,xtol,gtol,ftol,niter &
         ,iprint,iflag,myid,func_w_pmd,cfmethod &
         ,niter_eval,write_stats)
  else
    if(myid.eq.0) print *,'Metadynamics is not available for '//&
         trim(cpot)
  endif

  return
end subroutine md_wrapper
!=======================================================================
subroutine ga_wrapper()
  use variables
  use parallel
  use minimize
  use fp_common,only: func_w_pmd, grad_w_pmd
  implicit none
  integer:: i,m
  real(8):: fval
  external:: write_stats,write_energy_relation

  if( trim(cpot).eq.'vcMorse' .or. trim(cpot).eq.'Morse' .or. &
       trim(cpot).eq.'EAM' .or. trim(cpot).eq.'NN' ) then
    call ga(nvars,vars,fval,vranges,xtol,gtol,ftol,niter &
         ,iprint,iflag,myid,func_w_pmd,cfmethod &
         ,niter_eval,write_stats,write_energy_relation)
  else
    if(myid.eq.0) print *,'Genetic Algorithm (GA) is not available for '//&
         trim(cpot)
  endif

  return
end subroutine ga_wrapper
!=======================================================================
subroutine de_wrapper()
  use variables
  use parallel
  use minimize
  use fp_common,only: func_w_pmd, grad_w_pmd
  implicit none
  integer:: i,m
  real(8):: fval
  external:: write_stats, write_energy_relation

  if( trim(cpot).eq.'vcMorse' .or. trim(cpot).eq.'Morse' .or. &
       trim(cpot).eq.'EAM' .or. trim(cpot).eq.'NN' ) then
    call de(nvars,vars,fval,vranges,xtol,gtol,ftol,niter &
         ,iprint,iflag,myid,func_w_pmd,cfmethod &
         ,niter_eval,write_stats, write_energy_relation)
  else
    if(myid.eq.0) print *,'Differential evolution (DE) is not available for '//&
         trim(cpot)
  endif

  return
end subroutine de_wrapper
!=======================================================================
subroutine pso_wrapper()
  use variables
  use parallel
  use minimize
  use fp_common,only: func_w_pmd, grad_w_pmd
  implicit none
  integer:: i,m
  real(8):: fval
  external:: write_stats

  if( trim(cpot).eq.'vcMorse' .or. trim(cpot).eq.'Morse' .or. &
       trim(cpot).eq.'EAM' .or. trim(cpot).eq.'NN' ) then
    call pso(nvars,vars,fval,vranges,xtol,gtol,ftol,niter &
         ,iprint,iflag,myid,func_w_pmd,cfmethod &
         ,niter_eval,write_stats)
  else
    if(myid.eq.0) print *,'Particle Swarm Optimization (PSO) is'//&
         ' not available for '//trim(cpot)
  endif

  return
end subroutine pso_wrapper
!=======================================================================
subroutine random_search_wrapper()
  use variables
  use parallel
  use minimize
  use fp_common,only: func_w_pmd, grad_w_pmd
  implicit none
  integer:: i,m
  real(8):: fval
  external:: write_stats

  if( trim(cpot).eq.'vcMorse' .or. trim(cpot).eq.'Morse' .or. &
       trim(cpot).eq.'EAM' .or. trim(cpot).eq.'NN' ) then
    call random_search(nvars,vars,fval,vranges,xtol,gtol,ftol,niter &
         ,iprint,iflag,myid,func_w_pmd,cfmethod &
         ,niter_eval,write_stats)
  else
    if(myid.eq.0) print *,'Random search is not available for '//&
         trim(cpot)
  endif
  
end subroutine random_search_wrapper
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
  integer:: iter,istp,iv,i,nsize,narmijo
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
           ,iflag,myid,NN_fs,narmijo)
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
  use fp_common,only: func_w_pmd, grad_w_pmd
  implicit none
  integer:: iv
  real(8):: ftrn0,ftst0,ftmp,dv,vmax,ftst,ftmp1,ftmp2
  real(8),allocatable:: ganal(:),gnumer(:),vars0(:)
  real(8),parameter:: dev  = 1d-6
  real(8),parameter:: tiny = 1d-6

  allocate(gnumer(nvars),ganal(nvars),vars0(nvars))

  if( trim(cpot).eq.'NN' ) then
    call NN_init()
    call NN_func(nvars,vars,ftrn0,ftst0)
    call NN_grad(nvars,vars,ganal)
  else if( trim(cpot).eq.'vcMorse' ) then
    call func_w_pmd(nvars,vars,ftrn0,ftst)
    call grad_w_pmd(nvars,vars,ganal)
  endif

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
    dv = max(abs(vars(iv)*dev),dev)
    vars(iv)= vars(iv) +dv/2
    if( trim(cpot).eq.'NN' ) then
      call NN_func(nvars,vars,ftmp1,ftst)
    else if( trim(cpot).eq.'vcMorse' ) then
      call func_w_pmd(nvars,vars,ftmp1,ftst)
    endif
    vars(1:nvars)= vars0(1:nvars)
    vars(iv)= vars(iv) -dv/2
    if( trim(cpot).eq.'NN' ) then
      call NN_func(nvars,vars,ftmp2,ftst)
    else if( trim(cpot).eq.'vcMorse' ) then
      call func_w_pmd(nvars,vars,ftmp2,ftst)
    endif
    gnumer(iv)= (ftmp1-ftmp2)/dv
!!$    write(6,'(a,i5,10es15.7)') 'iv,var1,var2,ftmp1,ftmp2,gnumer = ', &
!!$         iv,vars0(iv)+dev/2,&
!!$         vars0(iv)-dev/2,ftmp1,ftmp2,gnumer(iv)
!!$    print *,''
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
  use fp_common,only: func_w_pmd, grad_w_pmd
  implicit none 
  integer:: iv
  real(8):: ftrn,ftst
  real(8),allocatable:: g(:)

  allocate(g(nvars))

  if( trim(cpot).eq.'NN' ) then
    call NN_init()
    call NN_func(nvars,vars,ftrn,ftst)
    call NN_grad(nvars,vars,g)
  else if( trim(cpot).eq.'vcMorse' .or. trim(cpot).eq.'Morse') then
    call func_w_pmd(nvars,vars,ftrn,ftst)
!!$    call grad_w_pmd(nvars,vars,g)
  endif

!!$  print *,'write_stats, myid=',myid
  call write_stats(0)

  if( myid.eq.0 ) then
    print *,'func values (training,test) =',ftrn,ftst
!!$    print *,'grad values (training):'
!!$    do iv=1,nvars
!!$      print *,'iv,grad(iv)=',iv,g(iv)
!!$    enddo
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
!!$  real(8):: epotsub
  
  cfname='out.erg'

  if( .not. allocated(erefl) ) allocate(erefl(nsmpl),erefg(nsmpl) &
       ,epotl(nsmpl),epotg(nsmpl),eerrl(nsmpl),eerrg(nsmpl)&
       ,swgtl(nsmpl),swgtg(nsmpl),esubl(nsmpl),esubg(nsmpl))

  if( l1st ) then
    erefl(1:nsmpl)= 0d0
    esubl(1:nsmpl)= 0d0
    eerrl(1:nsmpl)= 0d0
    swgtl(1:nsmpl)= 0d0
    do ismpl=isid0,isid1
      erefl(ismpl)= samples(ismpl)%eref
      esubl(ismpl)= samples(ismpl)%esub
      eerrl(ismpl)= samples(ismpl)%eerr
      swgtl(ismpl)= samples(ismpl)%wgt
    enddo
    erefg(1:nsmpl)= 0d0
    esubg(1:nsmpl)= 0d0
    eerrg(1:nsmpl)= 0d0
    swgtg(1:nsmpl)= 0d0
    call mpi_reduce(erefl,erefg,nsmpl,mpi_real8,mpi_sum &
         ,0,mpi_world,ierr)
    call mpi_reduce(esubl,esubg,nsmpl,mpi_real8,mpi_sum &
         ,0,mpi_world,ierr)
    call mpi_reduce(eerrl,eerrg,nsmpl,mpi_real8,mpi_sum &
         ,0,mpi_world,ierr)
    call mpi_reduce(swgtl,swgtg,nsmpl,mpi_real8,mpi_sum &
         ,0,mpi_world,ierr)
    do ismpl=1,nsmpl
      erefg(ismpl)= erefg(ismpl)/nalist(ismpl)
      esubg(ismpl)= esubg(ismpl)/nalist(ismpl)
    enddo
  endif

!!$  if( len(trim(crefstrct)).gt.5 ) then
!!$    if( myid.eq.myidrefsub ) then
!!$      epotsub = samples(isidrefsub)%epot
!!$    endif
!!$    call mpi_bcast(epotsub,1,mpi_integer,myidrefsub,mpi_world,ierr)
!!$  endif

  epotl(1:nsmpl)= 0d0
  do ismpl=isid0,isid1
    epotl(ismpl)= samples(ismpl)%epot +samples(ismpl)%esub
  enddo
  epotg(1:nsmpl)= 0d0
  call mpi_reduce(epotl,epotg,nsmpl,mpi_real8,mpi_sum &
       ,0,mpi_world,ierr)

  if( myid.eq.0 ) then
    open(90,file=trim(cfname)//'.trn.'//trim(cadd),status='replace')
    open(91,file=trim(cfname)//'.tst.'//trim(cadd),status='replace')
    write(90,'(a)') '# eref, epot, cdirname, diff, error, esub, swgt'
    write(91,'(a)') '# eref, epot, cdirname, diff, error, esub, swgt'
    do ismpl=1,nsmpl
      epotg(ismpl)= epotg(ismpl)/nalist(ismpl)
      if( iclist(ismpl).eq.1 ) then
        write(90,'(2es15.7,2x,a,10es13.4e3)') erefg(ismpl) &
             ,epotg(ismpl),trim(cdirlist(ismpl)) &
             ,abs(erefg(ismpl)-epotg(ismpl)) &
             ,eerrg(ismpl),esubg(ismpl),swgtg(ismpl)
      else if( iclist(ismpl).eq.2 ) then
        write(91,'(2es15.7,2x,a,10es12.3e3)') erefg(ismpl) &
             ,epotg(ismpl),trim(cdirlist(ismpl)) &
             ,abs(erefg(ismpl)-epotg(ismpl)) &
             ,eerrg(ismpl),esubg(ismpl),swgtg(ismpl)
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
  
  cfname= 'out.frc'

  nmaxl= 0
  do ismpl=1,nsmpl
    nmaxl= max(nmaxl,nalist(ismpl))
  enddo
  call mpi_allreduce(nmaxl,nmax,1,mpi_integer,mpi_max &
       ,mpi_world,ierr)

  if( .not. allocated(frefl) ) allocate(frefl(3,nmax,nsmpl)&
       ,frefg(3,nmax,nsmpl),fal(3,nmax,nsmpl),fag(3,nmax,nsmpl)&
       ,ferrl(nsmpl),ferrg(nsmpl),fsubl(3,nmax,nsmpl),fsubg(3,nmax,nsmpl))

  if( l1st ) then
    frefl(1:3,1:nmax,1:nsmpl)= 0d0
    fsubl(1:3,1:nmax,1:nsmpl)= 0d0
    ferrl(1:nsmpl) = 0d0
    do ismpl=isid0,isid1
      natm= samples(ismpl)%natm
      frefl(1:3,1:natm,ismpl)= samples(ismpl)%fref(1:3,1:natm)
      fsubl(1:3,1:natm,ismpl)= samples(ismpl)%fsub(1:3,1:natm)
      ferrl(ismpl) = samples(ismpl)%ferr
    enddo
    frefg(1:3,1:nmax,1:nsmpl)= 0d0
    fsubg(1:3,1:nmax,1:nsmpl)= 0d0
    ferrg(1:nsmpl) = 0d0
    call mpi_reduce(frefl,frefg,3*nmax*nsmpl,mpi_real8,mpi_sum &
         ,0,mpi_world,ierr)
    call mpi_reduce(fsubl,fsubg,3*nmax*nsmpl,mpi_real8,mpi_sum &
         ,0,mpi_world,ierr)
    call mpi_reduce(ferrl,ferrg,nsmpl,mpi_real8,mpi_sum &
         ,0,mpi_world,ierr)
  endif

  fal(1:3,1:nmax,1:nsmpl)= 0d0
  do ismpl=isid0,isid1
    natm= samples(ismpl)%natm
    fal(1:3,1:natm,ismpl)= samples(ismpl)%fa(1:3,1:natm) &
         +samples(ismpl)%fsub(1:3,1:natm)
  enddo
  fag(1:3,1:nmax,1:nsmpl)= 0d0
  call mpi_reduce(fal,fag,3*nmax*nsmpl,mpi_real8,mpi_sum &
       ,0,mpi_world,ierr)

  if( myid.eq.0 ) then
    open(92,file=trim(cfname)//'.trn.'//trim(cadd),status='replace')
    open(93,file=trim(cfname)//'.tst.'//trim(cadd),status='replace')
    do ismpl=1,nsmpl
      if( iclist(ismpl).eq.1 ) then
        natm= nalist(ismpl)
        do ia=1,natm
          do ixyz=1,3
            write(92,'(2es15.7,2x,a,i6,i3,3es12.3e3)') frefg(ixyz,ia,ismpl) &
                 ,fag(ixyz,ia,ismpl) &
                 ,trim(cdirlist(ismpl)),ia,ixyz &
                 ,abs(frefg(ixyz,ia,ismpl)-fag(ixyz,ia,ismpl))&
                 ,ferrg(ismpl),fsubg(ixyz,ia,ismpl)
          enddo
        enddo
      else if( iclist(ismpl).eq.2 ) then
        natm= nalist(ismpl)
        do ia=1,natm
          do ixyz=1,3
            write(93,'(2es15.7,2x,a,i6,i3,3es12.3e3)') frefg(ixyz,ia,ismpl) &
                 ,fag(ixyz,ia,ismpl) &
                 ,trim(cdirlist(ismpl)),ia,ixyz &
                 ,abs(frefg(ixyz,ia,ismpl)-fag(ixyz,ia,ismpl))&
                 ,ferrg(ismpl),fsubg(ixyz,ia,ismpl)
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
subroutine write_stress_relation(cadd)
  use variables
  use parallel
  implicit none
  character(len=*),intent(in):: cadd
  character(len=128):: cfname

  integer:: ismpl,ia,ixyz,jxyz,natm,nmax,nmaxl
  logical:: l1st = .true.

  cfname= 'out.strs'

  nmaxl= 0
  do ismpl=1,nsmpl
    nmaxl= max(nmaxl,nalist(ismpl))
  enddo
  call mpi_allreduce(nmaxl,nmax,1,mpi_integer,mpi_max &
       ,mpi_world,ierr)

  if( .not. allocated(srefl) ) allocate(srefl(3,3,nsmpl)&
       ,srefg(3,3,nsmpl),strsl(3,3,nsmpl),strsg(3,3,nsmpl)&
       ,serrl(nsmpl),serrg(nsmpl),ssubl(3,3,nsmpl),ssubg(3,3,nsmpl))

  if( l1st ) then
    srefl(1:3,1:3,1:nsmpl)= 0d0
    ssubl(1:3,1:3,1:nsmpl)= 0d0
    serrl(1:nsmpl) = 0d0
    do ismpl=isid0,isid1
      natm= samples(ismpl)%natm
      srefl(1:3,1:3,ismpl)= samples(ismpl)%sref(1:3,1:3)
      ssubl(1:3,1:3,ismpl)= samples(ismpl)%ssub(1:3,1:3)
      serrl(ismpl) = samples(ismpl)%serr
    enddo
    srefg(1:3,1:3,1:nsmpl)= 0d0
    ssubg(1:3,1:3,1:nsmpl)= 0d0
    serrg(1:nsmpl) = 0d0
    call mpi_reduce(srefl,srefg,3*3*nsmpl,mpi_real8,mpi_sum &
         ,0,mpi_world,ierr)
    call mpi_reduce(ssubl,ssubg,3*3*nsmpl,mpi_real8,mpi_sum &
         ,0,mpi_world,ierr)
    call mpi_reduce(serrl,serrg,nsmpl,mpi_real8,mpi_sum &
         ,0,mpi_world,ierr)
  endif

  strsl(1:3,1:3,1:nsmpl)= 0d0
  do ismpl=isid0,isid1
    strsl(1:3,1:3,ismpl)= samples(ismpl)%strs(1:3,1:3) &
         +samples(ismpl)%ssub(1:3,1:3)
  enddo
  strsg(1:3,1:3,1:nsmpl)= 0d0
  call mpi_reduce(strsl,strsg,3*3*nsmpl,mpi_real8,mpi_sum &
       ,0,mpi_world,ierr)

  if( myid.eq.0 ) then
    open(94,file=trim(cfname)//'.trn.'//trim(cadd),status='replace')
    open(95,file=trim(cfname)//'.tst.'//trim(cadd),status='replace')
    do ismpl=1,nsmpl
      if( iclist(ismpl).eq.1 ) then
        do ixyz=1,3
          do jxyz=ixyz,3
            write(94,'(2es15.7,2x,a,i6,i3,3es12.3e3)') srefg(ixyz,jxyz,ismpl) &
                 ,strsg(ixyz,jxyz,ismpl) &
                 ,trim(cdirlist(ismpl)),ixyz,jxyz &
                 ,abs(srefg(ixyz,jxyz,ismpl)-strsg(ixyz,jxyz,ismpl))&
                 ,serrg(ismpl),ssubg(ixyz,jxyz,ismpl)
          enddo
        enddo
      else if( iclist(ismpl).eq.2 ) then
        do ixyz=1,3
          do jxyz=ixyz,3
            write(95,'(2es15.7,2x,a,i6,i3,3es12.3e3)') srefg(ixyz,jxyz,ismpl) &
                 ,strsg(ixyz,jxyz,ismpl) &
                 ,trim(cdirlist(ismpl)),ixyz,jxyz &
                 ,abs(srefg(ixyz,jxyz,ismpl)-strsg(ixyz,jxyz,ismpl))&
                 ,serrg(ismpl),ssubg(ixyz,jxyz,ismpl)
          enddo
        enddo
      endif
    enddo
    close(94)
    close(95)
  endif
  l1st = .false.
end subroutine write_stress_relation
!=======================================================================
subroutine write_stats(iter)
  use variables
  use parallel
  use NNd
  implicit none
  integer,intent(in):: iter
  integer:: ismpl,natm,ntrnl,ntstl,ia,l,ntrn,ntst,nfcal,ixyz,jxyz
  type(mdsys)::smpl
  real(8):: de,df,epotsub,ds
  real(8):: demaxl_trn,demax_trn,desuml_trn,desum_trn,rmse_trn
  real(8):: demaxl_tst,demax_tst,desuml_tst,desum_tst,rmse_tst
  real(8):: dfmaxl_trn,dfmax_trn,dfsuml_trn,dfsum_trn
  real(8):: dfmaxl_tst,dfmax_tst,dfsuml_tst,dfsum_tst
  real(8):: dsmaxl_trn,dsmax_trn,dssuml_trn,dssum_trn
  real(8):: dsmaxl_tst,dsmax_tst,dssuml_tst,dssum_tst
  real(8),save:: rmse_tst_best= 1d+30
  character(len=128):: cnum
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

!!$!.....Restore subtracted energies and forces to get original reference values
!!$  if( nsubff.gt.0 ) then
!!$    call restore_FF()
!!$  endif

  if( len(trim(crefstrct)).gt.5 ) then
    if( myid.eq.myidrefsub ) then
      epotsub = samples(isidrefsub)%epot +samples(isidrefsub)%esub
      epotsub = epotsub /samples(isidrefsub)%natm
    endif
    call mpi_bcast(epotsub,1,mpi_real8,myidrefsub,mpi_world,ierr)
  endif

  demaxl_trn= 0d0
  desuml_trn= 0d0
  demaxl_tst= 0d0
  desuml_tst= 0d0
  do ismpl=isid0,isid1
    smpl= samples(ismpl)
    natm= smpl%natm
    if( len(trim(crefstrct)).gt.5 ) then
      de= abs(smpl%epot-epotsub*natm+smpl%esub &
           -(smpl%eref-erefsub))/natm
    else
      de= abs(smpl%epot+smpl%esub -smpl%eref)/natm
    endif
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
       ,mpi_real8,mpi_sum,0,mpi_world,ierr)
  call mpi_reduce(desuml_tst,desum_tst,1 &
       ,mpi_real8,mpi_sum,0,mpi_world,ierr)
  call mpi_reduce(demaxl_trn,demax_trn,1 &
       ,mpi_real8,mpi_max,0,mpi_world,ierr)
  call mpi_reduce(demaxl_tst,demax_tst,1 &
       ,mpi_real8,mpi_max,0,mpi_world,ierr)
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
  write(cnum,'(i0)') iter
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
          df= abs(smpl%fa(l,ia)+smpl%fsub(l,ia) -smpl%fref(l,ia))
          dfmaxl_trn= max(dfmaxl_trn,df)
          dfsuml_trn=dfsuml_trn +df*df
!!$          write(6,'(a,3i5,3es12.4)')  'ismpl,ia,l,fa,fref,dfsuml_trn=',&
!!$               ismpl,ia,l,smpl%fa(l,ia),smpl%fref(l,ia),dfsuml_trn
          ntrnl=ntrnl +1
        enddo
      enddo
    else if( smpl%iclass.eq.2 ) then
      do ia=1,natm
        if( smpl%ifcal(ia).eq.0 ) cycle
        do l=1,3
          df= abs(smpl%fa(l,ia)+smpl%fsub(l,ia) -smpl%fref(l,ia))
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
       ,mpi_real8,mpi_sum,0,mpi_world,ierr)
  call mpi_reduce(dfsuml_tst,dfsum_tst,1 &
       ,mpi_real8,mpi_sum,0,mpi_world,ierr)
  call mpi_reduce(dfmaxl_trn,dfmax_trn,1 &
       ,mpi_real8,mpi_max,0,mpi_world,ierr)
  call mpi_reduce(dfmaxl_tst,dfmax_tst,1 &
       ,mpi_real8,mpi_max,0,mpi_world,ierr)
  ntrn= 0
  ntst= 0
  call mpi_reduce(ntrnl,ntrn,1 &
       ,mpi_integer,mpi_sum,0,mpi_world,ierr)
  call mpi_reduce(ntstl,ntst,1 &
       ,mpi_integer,mpi_sum,0,mpi_world,ierr)
  rmse_trn= sqrt(dfsum_trn/ntrn)
!!$  print *,'dfsum_trn,ntrn = ',dfsum_trn,dfsuml_trn,ntrn,rmse_trn
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

!.....stress
  dsmaxl_trn = 0d0
  dssuml_trn = 0d0
  dsmaxl_tst = 0d0
  dssuml_tst = 0d0
  ntrnl = 0
  ntstl = 0
  do ismpl=isid0,isid1
    smpl = samples(ismpl)
    if( smpl%iclass.eq.1 ) then
      do ixyz=1,3
        do jxyz=1,3
          ds = abs(smpl%strs(ixyz,jxyz) +smpl%ssub(ixyz,jxyz) &
               -smpl%sref(ixyz,jxyz))
          dsmaxl_trn = max(dsmaxl_trn,ds)
          dssuml_trn = dssuml_trn +ds
          ntrnl = ntrnl +1
        enddo
      enddo
    else if( smpl%iclass.eq.2 ) then
      do ixyz=1,3
        do jxyz=1,3
          ds = abs(smpl%strs(ixyz,jxyz) +smpl%ssub(ixyz,jxyz) &
               -smpl%sref(ixyz,jxyz))
          dsmaxl_tst = max(dsmaxl_tst,ds)
          dssuml_tst = dssuml_tst +ds
          ntstl = ntstl +1
        enddo
      enddo
    endif
  enddo
  dssum_trn = 0d0
  dssum_tst = 0d0
  call mpi_reduce(dssuml_trn,dssum_trn,1,mpi_real8,mpi_sum, &
       0,mpi_world,ierr)
  call mpi_reduce(dssuml_tst,dssum_tst,1,mpi_real8,mpi_sum, &
       0,mpi_world,ierr)
  call mpi_reduce(dsmaxl_trn,dsmax_trn,1,mpi_real8,mpi_max, &
       0,mpi_world,ierr)
  call mpi_reduce(dsmaxl_tst,dsmax_tst,1,mpi_real8,mpi_max, &
       0,mpi_world,ierr)
  ntrn = 0
  ntst = 0
  call mpi_reduce(ntrnl,ntrn,1 &
       ,mpi_integer,mpi_sum,0,mpi_world,ierr)
  call mpi_reduce(ntstl,ntst,1 &
       ,mpi_integer,mpi_sum,0,mpi_world,ierr)
  rmse_trn= sqrt(dssum_trn/ntrn)
  if( ntst.ne.0 ) then
    rmse_tst = sqrt(dssum_tst/ntst)
  else
    rmse_tst = 0d0
  endif
  if( myid.eq.0 ) then
    write(6,'(a,i8,f15.2,4(1x,f12.7))') ' STRESS: ' &
         ,iter,mpi_wtime()-time0 &
         ,rmse_trn,dsmax_trn,rmse_tst,dsmax_tst
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
  call mpi_bcast(crefstrct,128,mpi_character,0,mpi_world,ierr)

  call mpi_bcast(xtol,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(ftol,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(gtol,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(eatom,maxnsp,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(gscl,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(nfpsmpl,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(pwgt,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(ratio_test,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(rseed,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(force_limit,1,mpi_real8,0,mpi_world,ierr)
  
  call mpi_bcast(lematch,1,mpi_logical,0,mpi_world,ierr)
  call mpi_bcast(lfmatch,1,mpi_logical,0,mpi_world,ierr)
  call mpi_bcast(lsmatch,1,mpi_logical,0,mpi_world,ierr)
  call mpi_bcast(lgrad,1,mpi_logical,0,mpi_world,ierr)
  call mpi_bcast(lgscale,1,mpi_logical,0,mpi_world,ierr)

  call mpi_bcast(fupper_lim,1,mpi_real8,0,mpi_world,ierr)
!.....Simulated annealing
  call mpi_bcast(sa_temp0,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(sa_xw0,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(sa_tau,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(sa_div_best,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(sa_tctrl,128,mpi_character,0,mpi_world,ierr)
!.....Metadynamics
  call mpi_bcast(md_height,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(md_sigma,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(md_ng,1,mpi_integer,0,mpi_world,ierr)
!.....Genetic algorithm
  call mpi_bcast(ga_temp,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(ga_nbits,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(ga_nindivs,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(ga_noffsp,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(ga_rate_mutate,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(ga_fitness,128,mpi_character,0,mpi_world,ierr)
!.....Differential evolution
  call mpi_bcast(de_nindivs,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(de_frac,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(de_lambda,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(de_cross_rate,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(de_wmin,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(de_wmax,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(de_fitness,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(de_algo,128,mpi_character,0,mpi_world,ierr)
!.....Particle swarm optimization
  call mpi_bcast(pso_nindivs,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(pso_w,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(pso_c1,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(pso_c2,1,mpi_real8,0,mpi_world,ierr)
!.....sgd
  call mpi_bcast(csgdupdate,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(r0sgd,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(nsgdbsize,1,mpi_integer,0,mpi_world,ierr)
!.....initialize parameters
  call mpi_bcast(cinitv,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(vinitsgm,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(vinitmu,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(numtol,1,mpi_integer,0,mpi_world,ierr)
!.....CG
  call mpi_bcast(icgbtype,1,mpi_integer,0,mpi_world,ierr)
!.....L-BFGS
  call mpi_bcast(mstore,1,mpi_integer,0,mpi_world,ierr)
!.....Armijo
  call mpi_bcast(armijo_xi,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(armijo_tau,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(armijo_maxiter,1,mpi_integer,0,mpi_world,ierr)

  call mpi_bcast(nserr,1,mpi_integer,0,mpi_world,ierr)
  if( myid.gt.0 ) then
    allocate(cserr(nserr),seerr(nserr),sferr(nserr),sserr(nserr))
  endif
  call mpi_bcast(cserr,128*nserr,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(seerr,nserr,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(sferr,nserr,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(sserr,nserr,mpi_real8,0,mpi_world,ierr)

  call mpi_bcast(nswgt,1,mpi_integer,0,mpi_world,ierr)
  if( myid.gt.0 ) then
    allocate(cswgt(nswgt),swerg0(nswgt),swdenom(nswgt))
  endif
  call mpi_bcast(cswgt,128*nswgt,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(swerg0,nswgt,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(swdenom,nswgt,mpi_real8,0,mpi_world,ierr)

  call mpi_bcast(nsmpl_outfrc,1,mpi_integer,0,mpi_world,ierr)

!.....Force-fields to be subtracted from reference values
  call mpi_bcast(nsubff,1,mpi_integer,0,mpi_world,ierr)
  if( myid.gt.0 ) then
    allocate(csubffs(nsubff))
  endif
!.....TODO: check what happens if numff==0...
  call mpi_bcast(csubffs,20*nsubff,mpi_character,0,mpi_world,ierr)

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
        samples(ismpl)%serr= sferr(iserr)
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
  integer:: ismpl,iswgt,idx,natm
  real(8):: erg
  
  do ismpl=isid0,isid1
    natm = samples(ismpl)%natm
    erg = samples(ismpl)%eref /natm
    samples(ismpl)%wgt = 1d0
    do iswgt=1,nswgt
      idx = index(samples(ismpl)%cdirname,trim(cswgt(iswgt)))
      if( idx.ne.0 ) then
        samples(ismpl)%wgt = exp(-(erg -swerg0(iswgt))/swdenom(iswgt))
      endif
    enddo
!!$    write(6,'(a,i5,a20,3es12.4)') ' ismpl,cdirname,wgt,eref,erg=' &
!!$         ,ismpl,trim(samples(ismpl)%cdirname)&
!!$         ,samples(ismpl)%wgt,samples(ismpl)%eref,erg
  enddo
  if( myid.eq.0 ) then
    write(6,'(a)') ' Set sample weights of Botlzmann factor.'
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
subroutine subtract_FF()
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
  use Coulomb,only: set_paramsdir_Coulomb
  use Morse,only: set_paramsdir_Morse,set_params_vcMorse,set_params_Morse
  implicit none

  integer:: i,ismpl,natm
  logical:: lcalcgrad = .false.
  logical:: luse_Morse = .false.
  logical:: luse_Morse_repul = .false.
  logical:: luse_Coulomb = .false.
  logical,save:: l1st = .true.
  real(8):: epot,strs(3,3)
  real(8),save,allocatable:: frcs(:,:)

  if( l1st ) then
    do i=1,nsubff
      if( index(trim(csubffs(i)),'Morse').ne.0 ) then
        luse_Morse = .true.
      else if( index(trim(csubffs(i)),'Morse_repul').ne.0 ) then
        luse_Morse_repul = .true.
      else if( index(trim(csubffs(i)),'Coulomb').ne.0 .or. &
           index(trim(csubffs(i)),'vcGaussian').ne.0 ) then
        luse_Coulomb = .true.
      endif
    enddo
    if( myid.eq.0 .and. iprint.ne.0 ) then
      print '(/,a)','Force field to be subtracted:'
      do i=1,nsubff
        print *,'  i,FF = ',i,trim(csubffs(i))
      enddo
    endif

    if( .not.allocated(frcs) ) allocate(frcs(3,maxna))

!.....Only at the 1st call, perform pmd to get esubs
    do ismpl=isid0,isid1
      natm = samples(ismpl)%natm
      if( luse_Morse ) then
        call set_paramsdir_Morse(trim(cmaindir)//'/'&
             //trim(samples(ismpl)%cdirname)//'/pmd')
      else if( luse_Morse_repul ) then
        call set_paramsdir_Morse(trim(cmaindir)//'/'&
             //trim(samples(ismpl)%cdirname)//'/pmd')
      else if( luse_Coulomb ) then
        call set_paramsdir_Coulomb(trim(cmaindir)//'/'&
             //trim(samples(ismpl)%cdirname)//'/pmd')
      endif
      call run_pmd(samples(ismpl),lcalcgrad,nvars,gvar,&
           nsubff,csubffs,epot,frcs,strs)
!!$      print *,'myid,ismpl,epot=',myid,ismpl,epot
      samples(ismpl)%esub = epot
      samples(ismpl)%fsub(1:3,1:natm) = frcs(1:3,1:natm)
      samples(ismpl)%ssub(1:3,1:3) = strs(1:3,1:3)
    enddo

!!$    allocate(esubl(nsmpl),esubg(nsmpl))
!!$    esubl(1:nsmpl) = 0d0
!!$    do ismpl=isid0,isid1
!!$      esubl(ismpl) = samples(ismpl)%esub
!!$    enddo
!!$    esubg(1:nsmpl) = 0d0
!!$    call mpi_reduce(esubl,esubg,nsmpl,mpi_real8,mpi_sum,0,mpi_world,ierr)
!!$
!!$    if( myid.eq.0 ) then
!!$      esubg(1:nsmpl)= esubg(1:nsmpl)/nalist(1:nsmpl)
!!$      open(92,file='out.erg.subtract',status='replace')
!!$      do ismpl=1,nsmpl
!!$        write(92,'(i8,es15.7,2x,a)') ismpl,esubg(ismpl),trim(cdirlist(ismpl))
!!$      enddo
!!$      close(92)
!!$    endif
!!$    deallocate(esubl,esubg)
  endif

!!$!.....After the 1st call, subtract esubs calculated at the 1st call
!!$  do ismpl=isid0,isid1
!!$!.....Subtract energy and forces from eref and fref, respectively
!!$    write(6,'(a,i8,1x,a,3es15.7)') 'ismpl,cdirname,eref,epot,esub=',ismpl,&
!!$         trim(samples(ismpl)%cdirname),samples(ismpl)%eref,&
!!$         samples(ismpl)%epot,samples(ismpl)%esub
!!$    samples(ismpl)%eref = samples(ismpl)%eref -samples(ismpl)%esub
!!$    samples(ismpl)%epot = samples(ismpl)%epot -samples(ismpl)%esub
!!$    do i=1,samples(ismpl)%natm
!!$      samples(ismpl)%fref(1:3,i) = samples(ismpl)%fref(1:3,i) &
!!$           -samples(ismpl)%fsub(1:3,i)
!!$      samples(ismpl)%fa(1:3,i) = samples(ismpl)%fa(1:3,i) &
!!$           -samples(ismpl)%fsub(1:3,i)
!!$    enddo
!!$    write(6,'(a,i8,1x,a,2es15.7)') ' ismpl,cdirname,eref=',ismpl, &
!!$         trim(samples(ismpl)%cdirname),samples(ismpl)%eref
!!$!.....TODO: stress should also come here.
!!$
!!$  enddo

  l1st = .false.
  return
end subroutine subtract_FF
!=======================================================================
subroutine restore_FF()
!
!  Restore subtracted energies and forces
!
  use variables
  use parallel
  implicit none

  integer:: i,ismpl

!!$  print *,'restore_FF'
  do ismpl=isid0,isid1
!!$    write(6,*) 'ismpl,eref,epot,esub=',ismpl,samples(ismpl)%eref,&
!!$         samples(ismpl)%epot,samples(ismpl)%esub
    samples(ismpl)%eref = samples(ismpl)%eref +samples(ismpl)%esub
    samples(ismpl)%epot = samples(ismpl)%epot +samples(ismpl)%esub
!!$    write(6,*) 'ismpl,eref,epot,esub=',ismpl,samples(ismpl)%eref,&
!!$         samples(ismpl)%epot,samples(ismpl)%esub
    do i=1,samples(ismpl)%natm
      samples(ismpl)%fref(1:3,i) = samples(ismpl)%fref(1:3,i) &
           +samples(ismpl)%fsub(1:3,i)
      samples(ismpl)%fa(1:3,i) = samples(ismpl)%fa(1:3,i) &
           +samples(ismpl)%fsub(1:3,i)
    enddo
!.....TODO: stress should also come here.
  enddo
  
end subroutine restore_FF
!=======================================================================
subroutine run_pmd(smpl,lcalcgrad,ndimp,pderiv,nff,cffs,epot,frcs, &
     strs)
!
!  Run pmd and get energy and forces of the system.
!  TODO: stress should be returned as well.
!
  use variables,only: rcut,mdsys,maxna,iprint,lsmatch
  use parallel,only: myid_pmd,mpi_comm_pmd,nnode_pmd,myid,mpi_world
  implicit none
  include "../pmd/params_unit.h"
  type(mdsys),intent(inout):: smpl
  integer,intent(in):: ndimp,nff
  real(8),intent(inout):: pderiv(ndimp),epot,frcs(3,maxna)
  real(8),intent(out):: strs(3,3)
  logical,intent(in):: lcalcgrad
  character(len=20),intent(in):: cffs(nff)

  logical,save:: l1st = .true.

  integer:: i,maxstp,nerg,npmd,ifpmd,ifdmp,minstp,n_conv,ifsort, &
       ifcoulomb,nismax,nstps_done,ntdst,nx,ny,nz,iprint_pmd
  real(8):: am(9),dt,rc,rbuf,dmp,tinit,tfin,ttgt(9),trlx,stgt(3,3),&
       ptgt,srlx,stbeta,strfin,fmv(3,0:9),ptnsr(3,3),ekin,eps_conv
  logical:: ltdst,lstrs,lcellfix(3,3),lvc
  character:: ciofmt*6,ctctl*20,cpctl*20,czload_type*5

  if( l1st ) then
!.....Create MPI COMM for pmd only for the 1st time
    call create_mpi_comm_pmd()
    l1st = .false.
  endif

  maxstp = 0
  nismax = 9
  nerg = 1
  npmd = 1
  am(1:9) = 1d0  ! Since no dynamics, no need of mass
  dt = 5d0
  ciofmt = 'ascii'
  ifpmd = 0
  rc = rcut
  rbuf = 0.0d0
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
!!$  lstrs = lsmatch
  cpctl = 'none'
  stgt(1:3,1:3) = 0d0
  ptgt = 0d0
  srlx = 100d0
  stbeta = 1d-1
  strfin = 0d0
  fmv(1:3,0) = (/ 0d0, 0d0, 0d0 /)
  fmv(1:3,1:9) = 1d0
  ptnsr(1:3,1:3) = 0d0
  ekin = 0d0
  n_conv = 1
  czload_type = 'no'
  eps_conv = 1d-3
  ifsort = 1
  ifcoulomb = 0
  lcellfix(1:3,1:3) = .false.
  nx = 1
  ny = 1
  nz = 1
  iprint_pmd = 0
  if( iprint.ge.100 ) iprint_pmd = 10

  lvc = .false.
  do i=1,nff
    if( trim(cffs(i)).eq.'long_Coulomb' .or. &
         trim(cffs(i)).eq.'vcMorse' ) then
      lvc = .true.
    endif
  enddo
  
!.....one_shot force calculation
!!$  print *,'calling one_shot, myid,mpi_world,myid_pmd,mpi_comm_pmd='&
!!$       ,myid,mpi_world,myid_pmd,mpi_comm_pmd
  call one_shot(smpl%h0,smpl%h,smpl%natm,smpl%tag,smpl%ra &
       ,smpl%va,frcs,smpl%strsi,smpl%eki,smpl%epi &
       ,smpl%chg,smpl%chi &
       ,myid_pmd,mpi_comm_pmd,nnode_pmd,nx,ny,nz &
       ,nismax,am,dt,nff,cffs,rc,rbuf,ptnsr,epot,ekin &
       ,ifcoulomb,iprint_pmd,lcalcgrad,ndimp,pderiv,lvc)
  strs(1:3,1:3) = ptnsr(1:3,1:3)*up2gpa*(-1d0)
!!$  print *,'one_shot done, cdirname,epot = ',trim(smpl%cdirname),epot
!!$  print *,'smpl%natm =',smpl%natm
!!$  write(6,'(a,30es12.4)') 'smpl%epi=',(smpl%epi(i),i=1,smpl%natm)

  return
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

!!$  call mpi_comm_group(mpi_world, mpi_group_world,ierr)
  iranks(1) = myid
!!$  call mpi_group_incl(mpi_group_world, 1, iranks, mpi_group_pmd,ierr)
!!$  call mpi_comm_create_group(mpi_world, mpi_group_pmd, 0, mpi_comm_pmd,ierr)
!!$  print *,'myid,iranks(1),mpi_group_pmd,mpi_comm_pmd=',&
!!$       myid,iranks(1),mpi_group_pmd,mpi_comm_pmd

  call mpi_comm_split(mpi_world,myid,myid,mpi_comm_pmd,ierr)

  call mpi_comm_size(mpi_comm_pmd,nnode_pmd,ierr)
  call mpi_comm_rank(mpi_comm_pmd,myid_pmd,ierr)
  call mpi_comm_group(mpi_comm_pmd,mpi_group_pmd,ierr)
!!$  print *,'myid,mpi_world,myid_pmd,mpi_comm_pmd,mpi_group_pmd = ',&
!!$       myid,mpi_world,myid_pmd,mpi_comm_pmd,mpi_group_pmd
  
!!$  call mpi_group_free(mpi_group_world,ierr)
!!$  call mpi_group_free(mpi_group_pmd,ierr)
  if( myid.eq.0 ) then
    write(6,'(a)') ''
    write(6,'(a)') ' MPI_COMM_PMD was created at each node '// &
         'for pmd calculations.'
    write(6,'(a)') ''
  endif
  
end subroutine create_mpi_comm_pmd
!=======================================================================
subroutine subtract_atomic_energy()
  use variables
  use parallel
  implicit none
  integer:: ismpl,is,i
  type(mdsys):: smpl

  do ismpl=isid0,isid1
    smpl = samples(ismpl)
    do i=1,smpl%natm
      is= int(smpl%tag(i))
      samples(ismpl)%eref= samples(ismpl)%eref -eatom(is)
    enddo
  enddo

  if( myid.eq.0 ) then
    write(6,'(a)') ' Subtracted atomic energies from reference energies.'
  endif

end subroutine subtract_atomic_energy
!=======================================================================
subroutine subtract_ref_struct_energy()
  use variables
  use parallel
  implicit none
  integer:: ismpl,myidrefsubl

  do ismpl=isid0,isid1
    if( trim(samples(ismpl)%cdirname).eq.trim(crefstrct) ) then
      myidrefsub = myid
      isidrefsub = ismpl
      erefsub = samples(ismpl)%eref
    endif
  enddo

  myidrefsubl = myidrefsub
  call mpi_allreduce(myidrefsubl,myidrefsub,1,mpi_integer,mpi_max,&
       mpi_world,ierr)
  call mpi_bcast(erefsub,1,mpi_real8,myidrefsub,mpi_world,ierr)

  if( myid.eq.0 .and. iprint.ne.0 ) then
    write(6,'(a,es12.4,a)') ' reference structure energy = ', &
         erefsub,' eV'
  endif
!!$  print *,'myid,myidrefsubl,myidrefsub,erefsub=',myid,&
!!$       myidrefsubl,myidrefsub,erefsub
  
end subroutine subtract_ref_struct_energy
!=======================================================================
subroutine set_max_num_atoms()
  use variables
  use parallel
  integer:: ismpl, maxnal

  maxnal = 0
  do ismpl=isid0,isid1
    na = samples(ismpl)%natm
    maxnal = max(na,maxnal)
  enddo
  maxna = 0
  call mpi_allreduce(maxnal,maxna,1,mpi_integer,mpi_max,mpi_world,ierr)
  if( myid.eq.0 .and. iprint.ne.0 ) then
    print *,'max num of atoms among samples = ',maxna
  endif
  
end subroutine set_max_num_atoms
!-----------------------------------------------------------------------
! Local Variables:
! compile-command: "make fitpot"
! End:
