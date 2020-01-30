program fitpot
!-----------------------------------------------------------------------
!                     Last modified: <2020-01-30 14:54:58 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
  use variables
  use parallel
!!$  use NNd,only:NN_init,NN_func,NN_grad
  use fp_common,only: func_w_pmd, grad_w_pmd, write_dsgnmats &
       ,subtract_FF, restore_FF, normalize
  use minimize
  use version
  use NN2,only: set_iglid_NN2
  use linreg,only: set_iglid_linreg
  use util,only: time_stamp
  implicit none
  integer:: ismpl,ihour,imin,isec
  real(8):: tmp,ftrn0,ftst0

  call mpi_init(ierr)
  time0= mpi_wtime()
  call mpi_comm_size(mpi_comm_world,nnode,ierr)
  call mpi_comm_rank(mpi_comm_world,myid,ierr)
  mpi_world= mpi_comm_world
  tcomm= 0d0

  call init_variables()

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
!.....NN and NN2 are both pointing NN2
    if( trim(cpot).eq.'NN' ) cpot = 'NN2'
!.....Check GDW; GDW works only with ML potentials, which use descriptors
    if( index(cpot,'NN').eq.0  .and. trim(cpot).ne.'linreg' ) then
      if( lgdw ) print *,'Gaussian density weight only works for ML potentials, so unset GDW.'
      lgdw = .false.
    endif
!.....Normalization and GDW cannot be used with SGD
    if( trim(cfmethod).eq.'sgd' .or. trim(cfmethod).eq.'SGD' ) then
      if( lgdw ) print *,'Gaussian density weight is not available with SGD, so unset GDW.'
      lgdw = .false.
      if( lnormalize ) print *,'Normalization is not available with SGD, so unset normalization.'
      lnormalize = .false.
    endif
!.....GDW assumes that G's are normalized
    if( lgdw ) lnormalize = .true.
    call write_initial_setting()
  endif
  call sync_input()

  if( index(cpot,'NN').ne.0 .or. trim(cpot).eq.'linreg' ) call read_params_desc()

  call read_vars()
  allocate(gvar(nvars),dvar(nvars))
  mem = mem +8*size(gvar) +8*size(dvar)

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
  if( trim(cpot).eq.'BVS' ) then
    nff = 2
    allocate(cffs(nff))
    cffs(1) = 'Coulomb'
    cffs(2) = 'Morse'
  else if( trim(cpot).eq.'BVSx' ) then
    nff = 3
    allocate(cffs(nff))
    cffs(1) = 'Coulomb'
    cffs(2) = 'Morse'
    cffs(3) = 'angular'
  else if( trim(cpot).eq.'fpc' ) then
    nff = 2
    allocate(cffs(nff))
    cffs(1) = 'fpc'
    cffs(2) = 'Coulomb'
  else
    nff = 1
    allocate(cffs(nff))
    cffs(1) = trim(cpot)
  endif

  if( (trim(cfmethod).ne.'test' .or. trim(cfmethod).ne.'dsgnmat') .and. &
       (trim(cpot).eq.'linreg' .or. trim(cpot).eq.'NN2') .or. trim(cpot).eq.'DNN' ) then
    lnormalize = .true.
  endif

!.....Initial computations of all samples
  if( trim(cpot).eq.'vcMorse' .or. trim(cpot).eq.'Morse' &
       .or. index(cpot,'BVS').ne.0 .or. trim(cpot).eq.'linreg' &
       .or. trim(cpot).eq.'NN2' .or. trim(cpot).eq.'BMH' &
       .or. trim(cpot).eq.'Abell' .or. trim(cpot).eq.'fpc' &
       .or. trim(cpot).eq.'DNN' ) then
    call func_w_pmd(nvars,vars,ftrn0,ftst0)
!!$    if( lnormalize ) call normalize()
!!$    if( lgdw ) call compute_gdw()
  else
    print *,'ERROR: '//trim(cpot)//' is not available.'
    stop
  endif

!.....NN2 only, and should be called after func_w_pmd
  if( trim(cpot).eq.'NN2' ) then
    call set_iglid_NN2(cpena,cfmethod)
  else if( trim(cpot).eq.'linreg' ) then
    call set_iglid_linreg(cpena,cfmethod)
  endif

  select case (trim(cfmethod))
  case ('dsgnmat')
    if( trim(cpot).ne.'linreg' ) then
      if( myid.eq.0 ) print *,'dsgnmat is only available for linreg' &
           //' and not for '//trim(cpot)
    else if( nnode.ne.1 ) then
      if( myid.eq.0 ) print *,'dsgnmat is not available for parallel run.'
    else
      call write_dsgnmats()
    endif
    goto 99
  case ('sd','SD')
    call sd_wrapper(ftrn0,ftst0)
  case ('cg','CG')
    call cg_wrapper(ftrn0,ftst0)
  case ('bfgs','BFGS','dfp','DFP')
    call qn_wrapper(ftrn0,ftst0)
  case ('sa','SA')
    call sa_wrapper(ftrn0,ftst0)
  case ('ga','GA')
    call ga_wrapper(ftrn0,ftst0)
  case ('de','DE')
    call de_wrapper(ftrn0,ftst0)
  case ('pso','PSO')
    call pso_wrapper(ftrn0,ftst0)
  case ('md','metadynamics')
    call md_wrapper(ftrn0,ftst0)
  case ('random_search','random')
    call random_search_wrapper(ftrn0,ftst0)
  case ('fs','FS')
!!$    call fs_wrapper(ftrn0,ftst0)
    if( myid.eq.0 ) print *,'FS is not available in the current version.'
  case ('gfs')
    call gfs_wrapper(ftrn0,ftst0)
  case ('sgd','SGD')
    call sgd_wrapper(ftrn0,ftst0)
  case ('check_grad')
    call check_grad(ftrn0,ftst0)
  case ('test','TEST')
    call test(ftrn0,ftst0)
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
99 if( myid.eq.0 ) then
    write(6,'(a,2(2x,i0))') ' Number of func and grad calls =',nfunc, ngrad
!!$    write(6,'(a,i0)') ' Number of grad calls = ',ngrad
    write(6,'(a,f13.3,a)') ' Memory/proc = ',dble(mem)/1000/1000,' MB'
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
  write(6,'(2x,a25,2x,i0)') 'num_samples',nsmpl
  write(6,'(2x,a25,2x,i0)') 'num_iteration',niter
  write(6,'(2x,a25,2x,a)') 'fitting_method',trim(cfmethod)
  write(6,'(2x,a25,2x,a)') 'main_directory',trim(cmaindir)
  write(6,'(2x,a25,2x,a)') 'param_file',trim(cparfile)
  write(6,'(2x,a25,2x,a)') 'loss_function',trim(ctype_loss)
!!$  write(6,'(2x,a25,2x,a)') 'run_mode',trim(crunmode)
  write(6,'(2x,a25,2x,es12.3)') 'xtol',xtol
  write(6,'(2x,a25,2x,es12.3)') 'ftol',ftol
  write(6,'(2x,a25,2x,es12.3)') 'gtol',gtol
  write(6,'(2x,a25,2x,i0)') 'numtol',numtol
  write(6,'(2x,a25,10(2x,a3))') 'specorder',(trim(specorder(i)),i=1,nsp)
  if( len(trim(crefstrct)).gt.5 ) then
    write(6,'(2x,a25,2x,a)') 'reference_structure',trim(crefstrct)
  else
    do i=1,nspmax
      if( trim(specorder(i)).ne.'x' ) then
        write(6,'(2x,a25,2x,i2,a4,es15.7)') 'atom_energy',i,specorder(i),eatom(i)
      endif
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
  if( nspcs_neglect.gt.0 ) then
    write(6,'(2x,a25,9(2x,4a))') 'force_neglect_species',(cspcs_neglect(i),i=1,nspcs_neglect)
  endif

  if( nswgt.gt.0 ) then
    write(6,'(2x,a25,2x,i5)') 'sample_weight',nswgt
    do i=1,nswgt
      write(6,'(4x,a23,2x,f8.4,2x,f8.1)') trim(cswgt(i)), swerg0(i), swdenom(i)
    enddo
  endif

  if( lgdw ) then
    write(6,'(2x,a25,2x,l3)') 'gaussian_density_weight',lgdw
    write(6,'(2x,a25,2x,es12.3)') 'GDW_sigma',gdsgm
  endif
  
  write(6,'(2x,a25,2x,es12.3)') 'coeff_sequential',seqcoef
  write(6,'(2x,a25,2x,a)') 'line_minimization',trim(clinmin)
  write(6,'(a)') ''
  write(6,'(2x,a25,2x,i0)') 'sample_error',nserr
  do i=1,nserr
    write(6,'(4x,a23,3(1x,f10.4))') trim(cserr(i)), seerr(i), sferr(i), sserr(i)
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
  use util,only: num_data
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
!!$  integer,external:: num_data
  character(len=128):: ctmp 

  if( .not. allocated(cdirlist)) allocate(cdirlist(nsmpl))
  if( .not. allocated(iclist)) allocate(iclist(nsmpl))

  if( myid.eq.0 ) then
    lerror = .true.
    if( len(trim(csmplist)).lt.1 ) then
      print '(/,a)',' Sample list was created by performing the following command:'
      print *,'  $ ls '//trim(cmaindir) &
           //' | grep "smpl_" > dir_list.txt'
      call system('ls '//trim(cmaindir) &
           //' | grep "smpl_" > dir_list.txt')
      open(ionum,file='dir_list.txt',status='old')
    else
      print *,'Sample list was given by input.'
      open(ionum,file=trim(csmplist),status='old')
    endif
    read(ionum,'(a)',end=998) ctmp
    backspace(ionum)
    ndat = num_data(trim(ctmp),' ')
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
  
  if(myid.eq.0 .and. iprint.gt.1 ) print*,'Finished get_dir_list'
  return

999 continue
  if( myid.eq.0 ) print *,' ERROR: num_samples may be wrong.'
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
    print *,''
    print '(a,3(2x,i0))',' Number of samples (total,training,test) = ', &
         nsmpl,nsmpl_trn,nsmpl_tst
    if( iprint.gt.1 ) print *,'Finished set_training_test_with_ratio'
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
    print '(a,3(2x,i0))',' Number of samples (total,training,test) = ', &
         nsmpl,nsmpl_trn,nsmpl_tst
  endif
end subroutine count_training_test
!=======================================================================
subroutine read_samples()
  use variables
  use parallel
  use pmdio,only: csp2isp
  implicit none

  integer:: is,isp,jsp
  character:: cdir*128, cspmd*3, cspfp*3
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
!.....Specorder in fitpot and that in sample should be the same
    do isp=1,nspmax
      cspmd = samples(is)%specorder(isp)
      cspfp = specorder(isp)
      if( trim(cspfp).ne.trim(cspmd) .and. trim(cspmd).ne.'x' ) then
        print '(a)','ERROR: specorder in the sample is different from that in fitpot.'
        print '(a,2a5)','   failed species in fitpot and the sample: ',trim(cspfp),trim(cspmd)
        print '(a)','   specorder in fitpot and the sample, '//trim(cdir)
        do jsp=1,nspmax
          print '(i5,2a5)',jsp,specorder(jsp), samples(is)%specorder(jsp)
        enddo
        stop
      endif
    enddo
  enddo
  call mpi_reduce(nal,nalist,nsmpl,mpi_integer,mpi_sum &
       ,0,mpi_world,ierr)

  call mpi_barrier(mpi_world,ierr)
  if( myid.eq.0 .and. iprint.gt.1 ) then
    write(6,'(/,a)') ' Finished read_samples'
  endif
  deallocate(nal)
  return
end subroutine read_samples
!=======================================================================
subroutine read_pos(ionum,fname,ismpl,smpl)
  use variables
  use util,only: num_data
  implicit none 
  integer,intent(in):: ionum,ismpl
  character(len=*),intent(in):: fname
  type(mdsys),intent(inout):: smpl

  integer:: i,natm,num
  real(8):: tmp
  character(len=128):: cline
  character(len=10):: c1,copt

  open(ionum,file=trim(fname),status='old')
  do while(.true.)
    read(ionum,'(a)') cline
    if( cline(1:1).eq.'!' .or. cline(1:1).eq.'#' ) then
      if( index(cline,'specorder:').ne.0 ) then
        num = num_data(trim(cline),' ')
        if( num.gt.11 ) stop 'ERROR: number of species exceeds the limit.'
        read(cline,*) c1, copt, smpl%specorder(1:num-2)
!!$        print *,'specorder = ',smpl%specorder(1:num-2)
      endif
    else
      backspace(ionum)
      exit
    endif
  enddo
  read(ionum,*) smpl%h0
  read(ionum,*) smpl%h(1:3,1)
  read(ionum,*) smpl%h(1:3,2)
  read(ionum,*) smpl%h(1:3,3)
  read(ionum,*) tmp,tmp,tmp
  read(ionum,*) tmp,tmp,tmp
  read(ionum,*) tmp,tmp,tmp
  read(ionum,*) natm
  smpl%h(1:3,1:3) = smpl%h(1:3,1:3)*smpl%h0
  smpl%natm= natm
  allocate(smpl%ra(3,natm),smpl%fa(3,natm) &
       ,smpl%tag(natm) &
       ,smpl%fref(3,natm), smpl%ifcal(natm),smpl%fabs(natm) &
       ,smpl%va(3,natm),smpl%strsi(3,3,natm) &
       ,smpl%eki(3,3,natm),smpl%epi(natm) &
       ,smpl%chg(natm),smpl%chi(natm),smpl%tei(natm),smpl%fsub(3,natm) &
       ,smpl%eatm(natm) &
       ,smpl%gwe(nvars),smpl%gwf(nvars,3,natm),smpl%gws(nvars,6))
  mem = mem +8*size(smpl%ra) +8*size(smpl%fa) +8*size(smpl%tag) &
       +8*size(smpl%fref) +4*size(smpl%ifcal) +8*size(smpl%fabs) &
       +8*size(smpl%va) +8*size(smpl%strsi) +8*size(smpl%eki) +8*size(smpl%epi) &
       +8*size(smpl%chg) +8*size(smpl%chi) +8*size(smpl%tei) +8*size(smpl%fsub) &
       +8*size(smpl%eatm) +8*size(smpl%gwe) +8*size(smpl%gwf) +8*size(smpl%gws)
  if( lgdw ) then
    allocate(smpl%gdf(natm),smpl%gdw(natm))
    mem = mem +8*size(smpl%gdf) +8*size(smpl%gdw)
  endif
  smpl%chg(1:natm) = 0d0
  smpl%tei(1:natm) = 0d0
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
  integer:: ismpl,i,is,jflag,natm,nfrc,nftot,nfrcg,nftotg,ispmax,ispmaxl
  integer:: nfrefdat,ifcal
  character(len=128):: cdir
  real(8):: erefminl,ftmp(3),ptnsr(3,3)
  character(len=3):: cspi 

  jflag= 0
  erefminl= 0d0
  nftot= 0
  nfrc = 0
  ispmaxl = 0
  do ismpl=isid0,isid1
    cdir=samples(ismpl)%cdirname
    open(13,file=trim(cmaindir)//'/'//trim(cdir) &
         //'/erg.ref',status='old')
    read(13,*) samples(ismpl)%eref
    close(13)
!.....Count numbers of each species
    samples(ismpl)%naps(1:nspmax) = 0
    ispmax = 0
    do i=1,samples(ismpl)%natm
      is= int(samples(ismpl)%tag(i))
      ispmax = max(ispmax,is)
      samples(ismpl)%naps(is) = samples(ismpl)%naps(is) +1
!!$      samples(ismpl)%eref= samples(ismpl)%eref  !-eatom(is)
    enddo
    ispmaxl = max(ispmaxl,ispmax)
    samples(ismpl)%ispmax = ispmax
    samples(ismpl)%ifcal(1:samples(ismpl)%natm)= 1
!!$    erefminl= min(erefminl,samples(ismpl)%eref/samples(ismpl)%natm)
!    write(6,*) 'ismpl,naps=',ismpl,samples(ismpl)%naps(1:nspmax)

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
      nftot= nftot + 3
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
      is= int(samples(ismpl)%tag(i))
      cspi = samples(ismpl)%specorder(is)
      if( csp_in_neglect(cspi) ) ifcal = 0
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
        nfrc = nfrc +3
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
  maxisp = 0
  call mpi_reduce(nfrc,nfrcg,1,mpi_integer,mpi_sum,0,mpi_world,ierr)
  call mpi_reduce(nftot,nftotg,1,mpi_integer,mpi_sum,0,mpi_world,ierr)
  call mpi_allreduce(ispmaxl,maxisp,1,mpi_integer,mpi_max,0,mpi_world,ierr)

  if( myid.eq.0 ) then
    print *,''
!    write(6,'(a,es12.4)') ' erefmin = ',erefmin
    if( iprint.gt.1 ) print '(a)',' Finished read_ref_data.'
    if( lfmatch ) then
      write(6,'(a,i0)') ' Number of forces to be used = ',nfrcg
      write(6,'(a,i0)') ' Total number of forces      = ',nftotg
    endif
    write(6,'(a,i0)') ' Number of species in all samples = ',maxisp
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
  integer:: ismpl,naps(nspmax),natm,ispcs,ielem
  real(8):: ebl(nspmax),erg
  
  ebl(1:nspmax) = 0d0
  do ismpl=isid0,isid1
    naps(1:nspmax) = samples(ismpl)%naps(1:nspmax)
    natm = samples(ismpl)%natm
!!$    write(6,*) ' ismpl,eref =',ismpl,samples(ismpl)%eref
!!$    write(6,*) ' ismpl,naps =',ismpl,naps(1:nspmax)
    ielem = 0
    do ispcs=1,nspmax
      if( naps(ispcs).eq.natm ) then
        ! this system is unary system
        erg = samples(ismpl)%eref/natm
        ielem = ispcs
        exit
      endif
    enddo
    if( ielem.ne.0 ) ebl(ielem) = min(ebl(ielem),erg)
  enddo

  ebase(1:nspmax) = 0d0
  call mpi_allreduce(ebl,ebase,nspmax,mpi_real8,mpi_min &
       ,mpi_world,ierr)

#ifdef _DEBUG
  if(myid.eq.0) then
    write(6,'(a)') ' Base energies obtained from unary systems:'
    do ispcs=1,nspmax
      write(6,'(a,i3,es12.4)') '   is, ebase(is) =',ispcs,ebase(ispcs)
    enddo
  endif
#endif

end subroutine get_base_energies
!=======================================================================
subroutine qn_wrapper(ftrn0,ftst0)
  use variables
!!$  use NNd,only:NN_init,NN_func,NN_grad,NN_restore_standard,NN_analyze
  use parallel
  use minimize
  use fp_common,only: func_w_pmd, grad_w_pmd
  implicit none
  real(8),intent(in):: ftrn0,ftst0
  real(8):: fval
  external:: write_stats

  if( trim(cpot).eq.'Morse' .or. trim(cpot).eq.'BVS' &
       .or. trim(cpot).eq.'linreg' &
       .or. trim(cpot).eq.'NN2' .or. trim(cpot).eq.'DNN' ) then
    call qn(nvars,vars,fval,gvar,dvar,vranges,xtol,gtol,ftol,niter &
         ,iprint,iflag,myid,func_w_pmd,grad_w_pmd,cfmethod &
         ,niter_eval,write_stats)
  else
    if( myid.eq.0 ) then
      print *,'Warning: BFGS is not available for '&
           //trim(cpot)
    endif
  endif

  return
end subroutine qn_wrapper
!=======================================================================
subroutine sd_wrapper(ftrn0,ftst0)
!
!  Steepest descent minimization
!
  use variables
!!$  use NNd,only:NN_init,NN_func,NN_grad
  use fp_common,only: func_w_pmd, grad_w_pmd
  use parallel
  use minimize
  implicit none
  real(8),intent(in):: ftrn0,ftst0
  real(8):: fval
  external:: write_stats

  call steepest_descent(nvars,vars,fval,gvar,dvar,vranges,xtol,gtol &
       ,ftol,niter,iprint,iflag,myid,func_w_pmd,grad_w_pmd,cfmethod &
       ,niter_eval,write_stats)

  return
end subroutine sd_wrapper
!=======================================================================
subroutine cg_wrapper(ftrn0,ftst0)
  use variables
!!$  use NNd,only:NN_init,NN_func,NN_grad,NN_restore_standard,NN_analyze
  use parallel
  use minimize
  use fp_common,only: func_w_pmd, grad_w_pmd
  implicit none
  real(8),intent(in):: ftrn0,ftst0
  real(8):: fval
  external:: write_stats

  if( trim(cpot).eq.'Morse' .or. trim(cpot).eq.'BVS' &
       .or. trim(cpot).eq.'linreg' &
       .or. trim(cpot).eq.'NN2' .or. trim(cpot).eq.'DNN' ) then
    call cg(nvars,vars,fval,gvar,dvar,vranges,xtol,gtol,ftol,niter &
         ,iprint,iflag,myid,func_w_pmd,grad_w_pmd,cfmethod &
         ,niter_eval,write_stats)
  else
    if( myid.eq.0 ) then
      print *,'Warning: CG is not available for '&
           //trim(cpot)
    endif
  endif
  
  return
end subroutine cg_wrapper
!=======================================================================
subroutine sa_wrapper(ftrn0,ftst0)
  use variables
  use parallel
  use minimize
  use fp_common,only: func_w_pmd, grad_w_pmd
  implicit none
  real(8),intent(in):: ftrn0,ftst0
  real(8):: fval
  external:: write_stats

  if( trim(cpot).eq.'vcMorse' .or. trim(cpot).eq.'Morse' &
       .or. trim(cpot).eq.'EAM' .or. trim(cpot).eq.'NN2' &
       .or. trim(cpot).eq.'DNN' ) then
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
subroutine md_wrapper(ftrn0,ftst0)
  use variables
  use parallel
  use minimize
  use fp_common,only: func_w_pmd, grad_w_pmd
  implicit none
  real(8),intent(in):: ftrn0,ftst0
  real(8):: fval
  external:: write_stats

  if( trim(cpot).eq.'vcMorse' .or. trim(cpot).eq.'Morse' &
       .or. trim(cpot).eq.'EAM' .or. trim(cpot).eq.'NN2' &
       .or. trim(cpot).eq.'DNN' ) then
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
subroutine ga_wrapper(ftrn0,ftst0)
  use variables
  use parallel
  use minimize
  use fp_common,only: func_w_pmd, grad_w_pmd
  implicit none
  real(8),intent(in):: ftrn0,ftst0
  real(8):: fval
  external:: write_stats,write_energy_relation

  if( trim(cpot).eq.'vcMorse' .or. trim(cpot).eq.'Morse' .or. &
       trim(cpot).eq.'EAM' .or. trim(cpot).eq.'NN2' .or. &
       index(cpot,'BVS').ne.0 .or. trim(cpot).eq.'linreg' .or. &
       trim(cpot).eq.'Abell' .or. trim(cpot).eq.'BMH' .or. &
       trim(cpot).eq.'fpc' ) then
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
subroutine de_wrapper(ftrn0,ftst0)
  use variables
  use parallel
  use minimize
  use fp_common,only: func_w_pmd, grad_w_pmd
  implicit none
  real(8),intent(in):: ftrn0,ftst0
  real(8):: fval
  external:: write_stats, write_energy_relation

  if( trim(cpot).eq.'vcMorse' .or. trim(cpot).eq.'Morse' .or. &
       trim(cpot).eq.'EAM' .or. trim(cpot).eq.'NN2' .or. &
       index(cpot,'BVS').ne.0 .or. trim(cpot).eq.'linreg' .or. &
       trim(cpot).eq.'Abell' .or. trim(cpot).eq.'BMH' .or. &
       trim(cpot).eq.'fpc' ) then
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
subroutine pso_wrapper(ftrn0,ftst0)
  use variables
  use parallel
  use minimize
  use fp_common,only: func_w_pmd, grad_w_pmd
  implicit none
  real(8),intent(in):: ftrn0,ftst0
  real(8):: fval
  external:: write_stats

  if( trim(cpot).eq.'vcMorse' .or. trim(cpot).eq.'Morse' .or. &
       trim(cpot).eq.'EAM' .or. trim(cpot).eq.'NN2' .or. &
       index(cpot,'BVS').ne.0 .or. trim(cpot).eq.'linreg' .or. &
       trim(cpot).eq.'BMH' .or. trim(cpot).eq.'Abell' .or. &
       trim(cpot).eq.'fpc'  ) then
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
subroutine random_search_wrapper(ftrn0,ftst0)
  use variables
  use parallel
  use minimize
  use fp_common,only: func_w_pmd, grad_w_pmd
  implicit none
  real(8),intent(in):: ftrn0,ftst0
  real(8):: fval
  external:: write_stats

  if( trim(cpot).eq.'vcMorse' .or. trim(cpot).eq.'Morse' .or. &
       trim(cpot).eq.'EAM' .or. trim(cpot).eq.'NN2' ) then
    call random_search(nvars,vars,fval,vranges,xtol,gtol,ftol,niter &
         ,iprint,iflag,myid,func_w_pmd,cfmethod &
         ,niter_eval,write_stats)
  else
    if(myid.eq.0) print *,'Random search is not available for '//&
         trim(cpot)
  endif
  
end subroutine random_search_wrapper
!=======================================================================
subroutine sgd_wrapper(ftrn0,ftst0)
!
! Wrapper for SGD
!
  use variables
  use minimize
  use parallel
  use fp_common,only: func_w_pmd, grad_w_pmd
  implicit none
  real(8),intent(in):: ftrn0,ftst0
  real(8):: fval
  external:: write_stats

  call sgd(nvars,vars,fval,gvar,dvar,vranges,xtol,gtol,ftol,niter &
       ,iprint,iflag,myid,mynsmpl,isid0,isid1,func_w_pmd,grad_w_pmd,cfmethod &
       ,niter_eval,write_stats)

end subroutine sgd_wrapper
!=======================================================================
subroutine fs_wrapper(ftrn0,ftst0)
  use variables
!!$  use NNd,only:NN_init,NN_func,NN_grad,NN_restore_standard,NN_analyze
  use parallel
  use minimize
  implicit none
  real(8),intent(in):: ftrn0,ftst0

  !.....NN specific code hereafter
!!$  call NN_init()
!!$  call fs(nvars,vars,fval,gvar,dvar,xtol,gtol,ftol,niter &
!!$       ,iprint,iflag,myid,NN_func,NN_grad)
!!$  call NN_analyze("fin")
!!$  call NN_restore_standard()

  return
end subroutine fs_wrapper
!=======================================================================
subroutine gfs_wrapper(ftrn0,ftst0)
  use variables
  use parallel
  use minimize
  use fp_common,only: func_w_pmd, grad_w_pmd
  implicit none
  real(8),intent(in):: ftrn0,ftst0
  real(8):: ftrn,ftst
  external:: write_stats

  if( trim(cpot).eq.'linreg' .or. trim(cpot).eq.'NN2' ) then
    ftrn = ftrn0
    ftst = ftst0
    call gfs(nvars,vars,ftrn,ftst,gvar,dvar,xtol,gtol,ftol,vranges,niter &
         ,iprint,iflag,myid,func_w_pmd,grad_w_pmd,cfmethod &
         ,niter_eval,write_stats)
  else
    if( myid.eq.0 ) then
      print *,'Warning: Group FS is not available for '&
           //trim(cpot)
    endif
  endif

  return
end subroutine gfs_wrapper
!=======================================================================
subroutine check_grad(ftrn0,ftst0)
  use variables
!!$  use NNd,only:NN_init,NN_func,NN_grad
  use parallel
  use fp_common,only: func_w_pmd, grad_w_pmd
  implicit none
  real(8),intent(in):: ftrn0,ftst0
  integer:: iv
  real(8):: dv,vmax,ftst,ftmp1,ftmp2
  real(8),allocatable:: ganal(:),gnumer(:),vars0(:)
  real(8),parameter:: dev  = 1d-5
  real(8),parameter:: tiny = 1d-6

  allocate(gnumer(nvars),ganal(nvars),vars0(nvars))

  call grad_w_pmd(nvars,vars,ganal)

  vars0(1:nvars)= vars(1:nvars)
  vmax= 0d0
  do iv=1,nvars
    vmax= max(vmax,abs(vars0(iv)))
  enddo
  dv= vmax *dev
  if( myid.eq.0 ) then
    print *,''
    print '(a,es12.4)',' Deviation for numerical derivative =',dv
  endif
  if( myid.eq.0 ) then
    write(6,'(a)') '------------------------------ check_grad '&
         //'------------------------------'
    write(6,'(a)') '     #,          x,    analytical,'// &
         '     numerical,'// &
         '      error [%]'
  endif
!.....Loop over variables for numerical derivative
  do iv=1,nvars
    vars(1:nvars)= vars0(1:nvars)
    dv = max(abs(vars(iv)*dev),dev)
    vars(iv)= vars(iv) +dv/2
    call func_w_pmd(nvars,vars,ftmp1,ftst)
    vars(1:nvars)= vars0(1:nvars)
    vars(iv)= vars(iv) -dv/2
    call func_w_pmd(nvars,vars,ftmp2,ftst)
    gnumer(iv)= (ftmp1-ftmp2)/dv
    if( myid.eq.0 ) then
      write(6,'(i6,es12.4,2es15.4,f15.3)') iv,vars0(iv), &
           ganal(iv) ,gnumer(iv), &
           abs((ganal(iv)-gnumer(iv))/(gnumer(iv)+tiny))*100
    endif
  enddo

  deallocate(gnumer,ganal,vars0)
  if( myid.eq.0 ) then
    write(6,'(a)') '----------------------------------------'&
         //'--------------------------------'
    print *, 'Finished check_grad'
  endif
end subroutine check_grad
!=======================================================================
subroutine test(ftrn0,ftst0)
  use variables
!!$  use NNd,only:NN_init,NN_func,NN_grad
  use parallel
  use fp_common,only: func_w_pmd, grad_w_pmd
  implicit none
  real(8),intent(in):: ftrn0,ftst0
  real(8),allocatable:: g(:)

  allocate(g(nvars))

!!$  if( trim(cpot).eq.'NN' ) then
!!$    call NN_init()
!!$    call NN_func(nvars,vars,ftrn,ftst)
!!$    call NN_grad(nvars,vars,g)
  if( trim(cpot).eq.'vcMorse' .or. trim(cpot).eq.'Morse' &
       .or. index(cpot,'BVS').ne.0 .or. trim(cpot).eq.'linreg' &
       .or. index(cpot,'NN').ne.0 ) then
!!$    call func_w_pmd(nvars,vars,ftrn,ftst)
    call grad_w_pmd(nvars,vars,g)
  else
    print *,'ERROR @test: '//trim(cpot)//' is not available for test().'
    stop
  endif

!!$  print *,'write_stats, myid=',myid
!!$  call write_stats(0)

  if( myid.eq.0 ) then
    print '(a,2es15.7)',' Loss func values (training,test) =',ftrn0,ftst0
!!$    print *,'grad values (training):'
!!$    do iv=1,nvars
!!$      print *,'iv,grad(iv)=',iv,g(iv)
!!$    enddo
    print *,'Finished test'
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
  
  integer:: ismpl
  logical,save:: l1st = .true.
!!$  real(8):: epotsub
  
  cfname='out.erg'

  if( .not. allocated(erefl) ) then
    allocate(erefl(nsmpl),erefg(nsmpl) &
       ,epotl(nsmpl),epotg(nsmpl),eerrl(nsmpl),eerrg(nsmpl)&
       ,swgtl(nsmpl),swgtg(nsmpl),esubl(nsmpl),esubg(nsmpl))
    mem = mem +8*size(erefl) +8*size(erefg) +8*size(epotl) +8*size(epotg) &
         +8*size(eerrl) +8*size(eerrg) +8*size(swgtl) +8*size(swgtg) &
         +8*size(esubl) +8*size(esubg)
  endif

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
  logical:: lfcal
  
  cfname= 'out.frc'

  nmaxl= 0
  do ismpl=1,nsmpl
    nmaxl= max(nmaxl,nalist(ismpl))
  enddo
  call mpi_allreduce(nmaxl,nmax,1,mpi_integer,mpi_max &
       ,mpi_world,ierr)

!.....This uses a lot of memory when the number of samples is large.
  if( .not. allocated(frefl) ) then
    allocate(frefl(3,nmax,nsmpl)&
         ,frefg(3,nmax,nsmpl),fal(3,nmax,nsmpl),fag(3,nmax,nsmpl)&
         ,ferrl(nsmpl),ferrg(nsmpl),fsubl(3,nmax,nsmpl) &
         ,fsubg(3,nmax,nsmpl),lfcall(nmax,nsmpl),lfcalg(nmax,nsmpl) &
         ,gdwl(nmax,nsmpl),gdwg(nmax,nsmpl))
    mem = mem +8*size(frefl) +8*size(frefg) +8*size(fal) +8*size(fag) &
         +8*size(ferrl) +8*size(ferrg) +8*size(fsubl) +8*size(fsubg) &
         +4*size(lfcall) +4*size(lfcalg) +8*size(gdwl) +8*size(gdwg)
  endif

  if( l1st ) then
    frefl(1:3,1:nmax,1:nsmpl)= 0d0
    fsubl(1:3,1:nmax,1:nsmpl)= 0d0
    ferrl(1:nsmpl) = 0d0
    lfcall(1:nmax,1:nsmpl) = .true.
    gdwl(:,:) = 0d0
    do ismpl=isid0,isid1
      natm= samples(ismpl)%natm
      frefl(1:3,1:natm,ismpl)= samples(ismpl)%fref(1:3,1:natm)
      fsubl(1:3,1:natm,ismpl)= samples(ismpl)%fsub(1:3,1:natm)
      ferrl(ismpl) = samples(ismpl)%ferr
      do ia=1,natm
        lfcall(ia,ismpl) = samples(ismpl)%ifcal(ia).eq.1
      enddo
      if( lgdw ) gdwl(1:natm,ismpl) = samples(ismpl)%gdw(1:natm)
    enddo
    frefg(1:3,1:nmax,1:nsmpl)= 0d0
    fsubg(1:3,1:nmax,1:nsmpl)= 0d0
    ferrg(1:nsmpl) = 0d0
    lfcalg(1:nmax,1:nsmpl) = .true.
    gdwg(:,:) = 0d0
    call mpi_reduce(frefl,frefg,3*nmax*nsmpl,mpi_real8,mpi_sum &
         ,0,mpi_world,ierr)
    call mpi_reduce(fsubl,fsubg,3*nmax*nsmpl,mpi_real8,mpi_sum &
         ,0,mpi_world,ierr)
    call mpi_reduce(ferrl,ferrg,nsmpl,mpi_real8,mpi_sum &
         ,0,mpi_world,ierr)
    call mpi_reduce(lfcall,lfcalg,nmax*nsmpl,mpi_logical,mpi_land &
         ,0,mpi_world,ierr)
    if( lgdw ) call mpi_reduce(gdwl,gdwg,nmax*nsmpl,mpi_real8,mpi_sum &
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
    write(92,'(a)') '# fref, fpot, cdirname, ia, ixyz, diff, error, fsub, gdw, lfcal'
    write(93,'(a)') '# fref, fpot, cdirname, ia, ixyz, diff, error, fsub, gdw, lfcal'
    do ismpl=1,nsmpl
      if( iclist(ismpl).eq.1 ) then
        natm= nalist(ismpl)
        do ia=1,natm
          lfcal = lfcalg(ia,ismpl)
          do ixyz=1,3
            write(92,'(2es12.4,2x,a,i6,i3,4es11.2e3,l3)') frefg(ixyz,ia,ismpl) &
                 ,fag(ixyz,ia,ismpl) &
                 ,trim(cdirlist(ismpl)),ia,ixyz &
                 ,abs(frefg(ixyz,ia,ismpl)-fag(ixyz,ia,ismpl))&
                 ,ferrg(ismpl),fsubg(ixyz,ia,ismpl),gdwg(ia,ismpl),lfcal
          enddo
        enddo
      else if( iclist(ismpl).eq.2 ) then
        natm= nalist(ismpl)
        do ia=1,natm
          lfcal = lfcalg(ia,ismpl)
          do ixyz=1,3
            write(93,'(2es12.4,2x,a,i6,i3,4es11.2e3,l3)') frefg(ixyz,ia,ismpl) &
                 ,fag(ixyz,ia,ismpl) &
                 ,trim(cdirlist(ismpl)),ia,ixyz &
                 ,abs(frefg(ixyz,ia,ismpl)-fag(ixyz,ia,ismpl))&
                 ,ferrg(ismpl),fsubg(ixyz,ia,ismpl),gdwg(ia,ismpl),lfcal
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

  integer:: ismpl,ixyz,jxyz,natm,nmax,nmaxl
  logical:: l1st = .true.

  cfname= 'out.strs'

  nmaxl= 0
  do ismpl=1,nsmpl
    nmaxl= max(nmaxl,nalist(ismpl))
  enddo
  call mpi_allreduce(nmaxl,nmax,1,mpi_integer,mpi_max &
       ,mpi_world,ierr)

  if( .not. allocated(srefl) ) then
    allocate(srefl(3,3,nsmpl)&
         ,srefg(3,3,nsmpl),strsl(3,3,nsmpl),strsg(3,3,nsmpl)&
         ,serrl(nsmpl),serrg(nsmpl),ssubl(3,3,nsmpl),ssubg(3,3,nsmpl))
    mem = mem +8*size(srefl) +8*size(srefg) +8*size(strsl) +8*size(strsg) &
         +8*size(serrl) +8*size(serrg) +8*size(ssubl) +8*size(ssubg)
  endif

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
            write(94,'(2es15.6e3,2x,a,i6,i3,3es12.3e3)') srefg(ixyz,jxyz,ismpl) &
                 ,strsg(ixyz,jxyz,ismpl) &
                 ,trim(cdirlist(ismpl)),ixyz,jxyz &
                 ,abs(srefg(ixyz,jxyz,ismpl)-strsg(ixyz,jxyz,ismpl))&
                 ,serrg(ismpl),ssubg(ixyz,jxyz,ismpl)
          enddo
        enddo
      else if( iclist(ismpl).eq.2 ) then
        do ixyz=1,3
          do jxyz=ixyz,3
            write(95,'(2es15.6e3,2x,a,i6,i3,3es12.3e3)') srefg(ixyz,jxyz,ismpl) &
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
!!$  use NNd
  implicit none
  integer,intent(in):: iter
  integer:: ismpl,natm,ntrnl,ntstl,ia,l,ntrn,ntst,nfcal,ixyz,jxyz
  type(mdsys)::smpl
  real(8):: de,df,ds
  real(8):: demaxl_trn,demax_trn,desuml_trn,desum_trn,rmse_trn
  real(8):: demaxl_tst,demax_tst,desuml_tst,desum_tst,rmse_tst
  real(8):: dfmaxl_trn,dfmax_trn,dfsuml_trn,dfsum_trn
  real(8):: dfmaxl_tst,dfmax_tst,dfsuml_tst,dfsum_tst
  real(8):: dsmaxl_trn,dsmax_trn,dssuml_trn,dssum_trn
  real(8):: dsmaxl_tst,dsmax_tst,dssuml_tst,dssum_tst
  real(8),save:: etrndnm,etstdnm
  real(8),save:: ftrndnm,ftstdnm
  real(8),save:: strndnm,ststdnm
  real(8):: er2trn,er2tst,fr2trn,fr2tst,sr2trn,sr2tst
  real(8),save:: epotsub
  character(len=128):: cnum
  logical,save:: l1st = .true.

  if( l1st ) then
    if( myid.eq.0 ) then
      write(6,*) '# ENERGY: ITER, TIME, ' &
           //'RMSE(TRAINING), RMSE(TEST), ' &
           //'MAX(TRAINING), MAX(TEST), ' &
           //'R^2(TRAINING), R^2(TEST), ' 
      write(6,*) '# FORCE:  ITER, TIME, ' &
           //'RMSE(TRAINING), RMSE(TEST), ' &
           //'MAX(TRAINING), MAX(TEST), ' &
           //'R^2(TRAINING), R^2(TEST)'
      write(6,*) '# STRESS:  ITER, TIME, ' &
           //'RMSE(TRAINING), RMSE(TEST), ' &
           //'MAX(TRAINING), MAX(TEST), ' &
           //'R^2(TRAINING), R^2(TEST)'
    endif

    call get_r2denom(etrndnm,etstdnm,ftrndnm,ftstdnm,strndnm,ststdnm)
    if( myid.eq.0 .and. iprint.gt.1 ) then
      print '(a,2es12.4)',' Denominator (energy) for trn, tst = ' &
           ,etrndnm, etstdnm
      print '(a,2es12.4)',' Denominator (force)  for trn, tst = ' &
           ,ftrndnm, ftstdnm
      print '(a,2es12.4)',' Denominator (stress) for trn, tst = ' &
           ,strndnm, ststdnm
    endif

  endif
  l1st = .false.

  epotsub = 0d0
  if( len(trim(crefstrct)).gt.5 ) then
    if( myid.eq.myidrefsub ) then
      epotsub = samples(isidrefsub)%epot +samples(isidrefsub)%esub
!!$      epotsub = samples(isidrefsub)%epot
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
    de= abs(smpl%epot-epotsub*natm+smpl%esub &
         -(smpl%eref-erefsub*natm))/natm
!!$    de= abs(smpl%epot-epotsub*natm &
!!$         -(smpl%eref-erefsub*natm-smpl%esub))/natm
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
  er2trn = 0d0
  if( etrndnm.gt.1d-15 ) er2trn = 1d0 -desum_trn/etrndnm
  er2tst = 0d0
  if( etstdnm.gt.1d-15 ) er2tst = 1d0 -desum_tst/etstdnm
  
  if( myid.eq.0 ) then
!!$    write(6,'(a,2i6)') ' nsmpl_trn, nsmpl_tst = ',nsmpl_trn,nsmpl_tst
    write(6,'(a,i8,f15.2,6(1x,f12.7))') ' ENERGY: ' &
         ,iter,mpi_wtime()-time0 &
         ,rmse_trn,rmse_tst &
         ,demax_trn,demax_tst &
         ,er2trn,er2tst
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
!!$          df= abs(smpl%fa(l,ia) -(smpl%fref(l,ia) -smpl%fsub(l,ia)))
          dfmaxl_trn= max(dfmaxl_trn,df)
          dfsuml_trn=dfsuml_trn +df*df
          ntrnl=ntrnl +1
        enddo
      enddo
    else if( smpl%iclass.eq.2 ) then
      do ia=1,natm
        if( smpl%ifcal(ia).eq.0 ) cycle
        do l=1,3
          df= abs(smpl%fa(l,ia)+smpl%fsub(l,ia) -smpl%fref(l,ia))
!!$          df= abs(smpl%fa(l,ia) -(smpl%fref(l,ia)) -smpl%fsub(l,ia))
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
  fr2trn = 0d0
  if( ftrndnm.gt.1d-15 ) fr2trn = 1d0 -dfsum_trn/ftrndnm
  fr2tst = 0d0
  if( ftstdnm.gt.1d-15 ) fr2tst = 1d0 -dfsum_tst/ftstdnm
  if( myid.eq.0 ) then
    write(6,'(a,i8,f15.2,6(1x,f12.7))') ' FORCE:  ' &
         ,iter,mpi_wtime()-time0 &
         ,rmse_trn,rmse_tst &
         ,dfmax_trn,dfmax_tst &
         ,fr2trn,fr2tst
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
          dssuml_trn = dssuml_trn +ds*ds
          ntrnl = ntrnl +1
        enddo
      enddo
    else if( smpl%iclass.eq.2 ) then
      do ixyz=1,3
        do jxyz=1,3
          ds = abs(smpl%strs(ixyz,jxyz) +smpl%ssub(ixyz,jxyz) &
               -smpl%sref(ixyz,jxyz))
          dsmaxl_tst = max(dsmaxl_tst,ds)
          dssuml_tst = dssuml_tst +ds*ds
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
  sr2trn = 0d0
  if( strndnm.gt.1d-15 ) sr2trn = 1d0 -dssum_trn/strndnm
  sr2tst = 0d0
  if( ststdnm.gt.1d-15 ) sr2tst = 1d0 -dssum_tst/ststdnm
  if( myid.eq.0 ) then
    write(6,'(a,i8,f15.2,6(1x,f12.7))') ' STRESS: ' &
         ,iter,mpi_wtime()-time0 &
         ,rmse_trn,rmse_tst &
         ,dsmax_trn,dsmax_tst &
         ,sr2trn,sr2tst
!!$    if( iprint.gt.1 ) then
!!$      print *,' dssum_trn,strndnm,sr2trn=',dssum_trn,strndnm,sr2trn
!!$      print *,' dssum_tst,ststdnm,sr2tst=',dssum_tst,ststdnm,sr2tst
!!$    endif
  endif

  return
end subroutine write_stats
!=======================================================================
subroutine get_r2denom(etrn,etst,ftrn,ftst,strn,stst)
!
!  Compute denominator for R^2 score.
!  It is called only once at the 1st call of the subroutine write_stats.
!  Subtract fixed FF from the reference values (erg,frc,strs),
!  to evaluate R^2 in the range of subtracted data.
!
  use variables
  use parallel
  implicit none
  real(8),intent(out):: etrn,etst,ftrn,ftst,strn,stst
  
  integer:: ismpl,ia,l,ixyz,jxyz,natm,nfcal,ntrn,ntst,ntrnl,ntstl
  type(mdsys)::smpl
  real(8):: eref,esub,tmp
  real(8):: esumltrn,esumltst,esumtrn,esumtst,emtrn,emtst
  real(8):: e2sumltrn,e2sumltst,e2sumtrn,e2sumtst,e2mtrn,e2mtst
  real(8):: fsumltrn,fsumltst,fsumtrn,fsumtst,fmtrn,fmtst
  real(8):: f2sumltrn,f2sumltst,f2sumtrn,f2sumtst,f2mtrn,f2mtst
  real(8):: ssumltrn,ssumltst,ssumtrn,ssumtst,smtrn,smtst
  real(8):: s2sumltrn,s2sumltst,s2sumtrn,s2sumtst,s2mtrn,s2mtst

!.....Energy
  esumltrn= 0d0
  esumltst= 0d0
  e2sumltrn= 0d0
  e2sumltst= 0d0
  do ismpl=isid0,isid1
    smpl= samples(ismpl)
    eref = smpl%eref
    esub = smpl%esub
    natm = smpl%natm
!!$    tmp = (eref-esub)/natm
    tmp = eref/natm
    if( smpl%iclass.eq.1 ) then
      esumltrn = esumltrn +tmp
      e2sumltrn= e2sumltrn +tmp*tmp
    else if( smpl%iclass.eq.2 ) then
      esumltst = esumltst +tmp
      e2sumltst= e2sumltst +tmp*tmp
    endif
  enddo
  esumtrn = 0d0
  e2sumtrn = 0d0
  esumtst = 0d0
  e2sumtst = 0d0
  call mpi_reduce(esumltrn,esumtrn,1,mpi_real8,mpi_sum,0,mpi_world,ierr)
  call mpi_reduce(e2sumltrn,e2sumtrn,1,mpi_real8,mpi_sum,0,mpi_world,ierr)
  call mpi_reduce(esumltst,esumtst,1,mpi_real8,mpi_sum,0,mpi_world,ierr)
  call mpi_reduce(e2sumltst,e2sumtst,1,mpi_real8,mpi_sum,0,mpi_world,ierr)
  emtrn = esumtrn/nsmpl_trn
  e2mtrn = e2sumtrn/nsmpl_trn
  if( nsmpl_tst.ne.0 ) then
    emtst = esumtst/nsmpl_tst
    e2mtst = e2sumtst/nsmpl_tst
  else
    emtst = 0d0
    e2mtst= 0d0
  endif
  if( iprint.gt.1 .and. myid.eq.0 ) then
    print *,'emtrn,emtrn^2,e2mtrn,nsmpl_trn =',emtrn,emtrn**2,e2mtrn,nsmpl_trn
    print *,'emtst,emtst^2,e2mtst,nsmpl_tst =',emtst,emtst**2,e2mtst,nsmpl_tst
  endif
  etrn = (e2mtrn -emtrn**2)*nsmpl_trn
  etst = (e2mtst -emtst**2)*nsmpl_tst

!.....Force
  fsumltrn = 0d0
  f2sumltrn = 0d0
  fsumltst = 0d0
  f2sumltst = 0d0
  ntrnl = 0
  ntstl = 0
  do ismpl=isid0,isid1
    smpl= samples(ismpl)
    nfcal= smpl%nfcal
    if( nfcal.eq.0 ) cycle
    natm = smpl%natm
    if( smpl%iclass.eq.1 ) then
      do ia=1,natm
        if( smpl%ifcal(ia).eq.0 ) cycle
        do l=1,3
!!$          tmp = smpl%fref(l,ia)-smpl%fsub(l,ia)
          tmp = smpl%fref(l,ia)
          fsumltrn = fsumltrn +tmp
          f2sumltrn= f2sumltrn +tmp*tmp
          ntrnl=ntrnl +1
        enddo
      enddo
    else if( smpl%iclass.eq.2 ) then
      do ia=1,natm
        if( smpl%ifcal(ia).eq.0 ) cycle
        do l=1,3
!!$          tmp = smpl%fref(l,ia)-smpl%fsub(l,ia)
          tmp = smpl%fref(l,ia)
          fsumltst = fsumltst +tmp
          f2sumltst= f2sumltst +tmp*tmp
          ntstl=ntstl +1
        enddo
      enddo
    endif
  enddo
  fsumtrn = 0d0
  f2sumtrn = 0d0
  fsumtst = 0d0
  f2sumtst = 0d0
  call mpi_reduce(fsumltrn,fsumtrn,1,mpi_real8,mpi_sum,0,mpi_world,ierr)
  call mpi_reduce(f2sumltrn,f2sumtrn,1,mpi_real8,mpi_sum,0,mpi_world,ierr)
  call mpi_reduce(fsumltst,fsumtst,1,mpi_real8,mpi_sum,0,mpi_world,ierr)
  call mpi_reduce(f2sumltst,f2sumtst,1,mpi_real8,mpi_sum,0,mpi_world,ierr)
  ntrn = 0
  ntst = 0
  call mpi_reduce(ntrnl,ntrn,1 &
       ,mpi_integer,mpi_sum,0,mpi_world,ierr)
  call mpi_reduce(ntstl,ntst,1 &
       ,mpi_integer,mpi_sum,0,mpi_world,ierr)
  fmtrn = fsumtrn/ntrn
  f2mtrn= f2sumtrn/ntrn
  if( ntst.ne.0 ) then
    fmtst = fsumtst/ntst
    f2mtst= f2sumtst/ntst
  else
    fmtst = 0d0
    f2mtst= 0d0
  endif
  if( iprint.gt.1 .and. myid.eq.0 ) then
    print *,'fmtrn,fmtrn^2,f2mtrn,ntrn =',fmtrn,fmtrn**2,f2mtrn,ntrn
    print *,'fmtst,fmtst^2,f2mtst,ntst =',fmtst,fmtst**2,f2mtst,ntst
  endif
  ftrn = (f2mtrn -fmtrn**2) *ntrn
  ftst = (f2mtst -fmtst**2) *ntst
  

  strn = 0d0
  stst = 0d0
!.....Stress
  ssumltrn = 0d0
  s2sumltrn = 0d0
  ssumltst = 0d0
  s2sumltst = 0d0
  ntrnl = 0
  ntstl = 0
  do ismpl=isid0,isid1
    smpl = samples(ismpl)
    if( smpl%iclass.eq.1 ) then
      do ixyz=1,3
        do jxyz=1,3
!!$          tmp = smpl%sref(ixyz,jxyz)-smpl%ssub(ixyz,jxyz)
          tmp = smpl%sref(ixyz,jxyz)
          ssumltrn = ssumltrn +tmp
          s2sumltrn= s2sumltrn +tmp*tmp
          ntrnl = ntrnl +1
        enddo
      enddo
    else if( smpl%iclass.eq.2 ) then
      do ixyz=1,3
        do jxyz=1,3
!!$          tmp = smpl%sref(ixyz,jxyz)-smpl%ssub(ixyz,jxyz)
          tmp = smpl%sref(ixyz,jxyz)
          ssumltst = ssumltst +tmp
          s2sumltst= s2sumltst +tmp*tmp
          ntstl = ntstl +1
        enddo
      enddo
    endif
  enddo
  ssumtrn = 0d0
  s2sumtrn = 0d0
  ssumtst = 0d0
  s2sumtst = 0d0
  call mpi_reduce(ssumltrn,ssumtrn,1,mpi_real8,mpi_sum,0,mpi_world,ierr)
  call mpi_reduce(s2sumltrn,s2sumtrn,1,mpi_real8,mpi_sum,0,mpi_world,ierr)
  call mpi_reduce(ssumltst,ssumtst,1,mpi_real8,mpi_sum,0,mpi_world,ierr)
  call mpi_reduce(s2sumltst,s2sumtst,1,mpi_real8,mpi_sum,0,mpi_world,ierr)
  ntrn = 0
  ntst = 0
  call mpi_reduce(ntrnl,ntrn,1 &
       ,mpi_integer,mpi_sum,0,mpi_world,ierr)
  call mpi_reduce(ntstl,ntst,1 &
       ,mpi_integer,mpi_sum,0,mpi_world,ierr)
  smtrn = ssumtrn/ntrn
  s2mtrn= s2sumtrn/ntrn
  if( ntst.ne.0 ) then
    smtst = ssumtst/ntst
    s2mtst= s2sumtst/ntst
  else
    smtst = 0d0
    s2mtst= 0d0
  endif
  if( iprint.gt.1 .and. myid.eq.0 ) then
    print *,'smtrn,smtrn^2,s2mtrn,ntrn =',smtrn,smtrn**2,s2mtrn,ntrn
    print *,'smtst,smtst^2,s2mtst,ntst =',smtst,smtst**2,s2mtst,ntst
  endif
  strn = (s2mtrn -smtrn**2) *ntrn
  stst = (s2mtst -smtst**2) *ntst
  
  return
end subroutine get_r2denom
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
  use random
  use pmdio,only: nnmax
  implicit none
  
  call mpi_bcast(nsmpl,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(niter,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(niter_eval,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(iprint,1,mpi_integer,0,mpi_world,ierr)

  call mpi_bcast(cfmethod,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(cmaindir,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(cparfile,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(crunmode,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(cevaltype,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(cpot,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(cpena,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(clinmin,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(cfsmode,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(crefstrct,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(ctype_loss,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(cfrc_denom,128,mpi_character,0,mpi_world,ierr)

  call mpi_bcast(xtol,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(ftol,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(gtol,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(eatom,nspmax,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(specorder,3*nspmax,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(gscl,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(nfpsmpl,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(pwgt,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(ratio_test,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(rseed,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(force_limit,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(rc_other,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(rcut,1,mpi_real8,0,mpi_world,ierr)
  
  call mpi_bcast(lematch,1,mpi_logical,0,mpi_world,ierr)
  call mpi_bcast(lfmatch,1,mpi_logical,0,mpi_world,ierr)
  call mpi_bcast(lsmatch,1,mpi_logical,0,mpi_world,ierr)
  call mpi_bcast(lgrad,1,mpi_logical,0,mpi_world,ierr)
  call mpi_bcast(lgscale,1,mpi_logical,0,mpi_world,ierr)
  call mpi_bcast(nspcs_neglect,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(cspcs_neglect,3*nspmax,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(interact,nspmax*nspmax,mpi_logical,0,mpi_world,ierr)
  call mpi_bcast(interact3,nspmax*nspmax*nspmax,mpi_logical,0,mpi_world,ierr)

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
  call mpi_bcast(de_temp,1,mpi_real8,0,mpi_world,ierr)
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
  call mpi_bcast(nsgdbsize,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(sgd_rate0,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(sgd_eps,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(adam_b1,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(adam_b2,1,mpi_real8,0,mpi_world,ierr)
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
!.....gFS
  call mpi_bcast(ninnergfs,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(cread_fsmask,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(cfs_xrefresh,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(maxfsrefresh,1,mpi_integer,0,mpi_world,ierr)

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
!.....GDW
  call mpi_bcast(lgdw,1,mpi_logical,0,mpi_world,ierr)
  call mpi_bcast(gdsgm,1,mpi_real8,0,mpi_world,ierr)
!.....Normalization
  call mpi_bcast(cnormalize,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(lnormalize,1,mpi_logical,0,mpi_world,ierr)
  
  call mpi_bcast(nsmpl_outfrc,1,mpi_integer,0,mpi_world,ierr)

!.....Force-fields to be subtracted from reference values
  call mpi_bcast(nsubff,1,mpi_integer,0,mpi_world,ierr)
  if( myid.gt.0 ) then
    allocate(csubffs(nsubff))
  endif
!.....TODO: check what happens if numff==0...
  call mpi_bcast(csubffs,20*nsubff,mpi_character,0,mpi_world,ierr)

!.....pmdio related
  call mpi_bcast(nnmax,1,mpi_integer,0,mpi_world,ierr)
  
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
  mem = mem +4*size(nspn) +4*size(ispn)
  do ip=1,nnode
    nspn(ip)= n
    if( ip.le.m ) nspn(ip)= nspn(ip) +1
  enddo
  mynsmpl= nspn(myid+1)
  call mpi_allreduce(mynsmpl,maxmynsmpl,1,mpi_integer,mpi_max &
       ,mpi_world,ierr)
  if( myid.eq.0 ) print '(a,i0)',' Max num of samples per process = ',maxmynsmpl

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

  if( myid.eq.0 .and. iprint.gt.1 ) print *,'Finished get_node2sample'
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
        samples(ismpl)%serr= sserr(iserr)
      endif
    enddo
  enddo
  if( myid.eq.0 .and. iprint.gt.1 ) then
    write(6,'(a)') ' Finished set_samples_errors'
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
subroutine subtract_atomic_energy()
  use variables
  use parallel
  use pmdio,only: csp2isp
  implicit none
  integer:: ismpl,is,i,isp
  type(mdsys):: smpl
  character(len=3):: csp 

  do ismpl=isid0,isid1
    smpl = samples(ismpl)
    do i=1,smpl%natm
!.....Convert species-ID in the sample (IS) to that in fitpot (ISP).
      is= int(smpl%tag(i))
      csp = smpl%specorder(is)
      isp = csp2isp(trim(csp),specorder)
!.....EATOM stores atomic energy according to the order of specorder in fitpot.
      samples(ismpl)%eref= samples(ismpl)%eref -eatom(isp)
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

  myidrefsub = -1
  do ismpl=isid0,isid1
    if( trim(samples(ismpl)%cdirname).eq.trim(crefstrct) ) then
      myidrefsub = myid
      isidrefsub = ismpl
      erefsub = samples(ismpl)%eref /samples(ismpl)%natm
    endif
  enddo

  myidrefsubl = myidrefsub
  call mpi_allreduce(myidrefsubl,myidrefsub,1,mpi_integer,mpi_max,&
       mpi_world,ierr)

  if( myidrefsub.lt.0 ) then
    print *,'Error: No reference sample, '//trim(crefstrct)
    call mpi_finalize()
    stop
  endif
  call mpi_bcast(erefsub,1,mpi_real8,myidrefsub,mpi_world,ierr)
  if( myid.eq.0 .and. iprint.ne.0 ) then
    write(6,'(a,es12.4,a)') ' Reference structure energy = ', &
         erefsub,' eV/atom'
  endif
!!$  print *,'myid,myidrefsubl,myidrefsub,erefsub=',myid,&
!!$       myidrefsubl,myidrefsub,erefsub
  
end subroutine subtract_ref_struct_energy
!=======================================================================
subroutine set_max_num_atoms()
  use variables
  use parallel
  integer:: ismpl, maxnal, maxninl

  maxnal = 0
  maxninl = 0
  do ismpl=isid0,isid1
    na = samples(ismpl)%natm
    maxnal = max(na,maxnal)
    maxninl = maxninl +na
  enddo
  maxna = 0
  maxnin = 0
  natot = 0
  call mpi_allreduce(maxnal,maxna,1,mpi_integer,mpi_max,mpi_world,ierr)
  call mpi_allreduce(maxninl,maxnin,1,mpi_integer,mpi_max,mpi_world,ierr)
  call mpi_allreduce(maxninl,natot,1,mpi_integer,mpi_sum,mpi_world,ierr)
  if( myid.eq.0 .and. iprint.ne.0 ) then
    write(6,'(a,i0)') ' Max num of atoms among samples   = ',maxna
    write(6,'(a,i0)') ' Max num of atoms among nodes     = ',maxnin
    write(6,'(a,i0)') ' Total num of atoms among samples = ',natot
  endif
  
end subroutine set_max_num_atoms
!=======================================================================
function string_in_arr(string,narr,array)
  implicit none
  integer,intent(in):: narr
  character(len=*),intent(in):: string
  character(len=*),intent(in):: array(narr)
  logical:: string_in_arr
  integer:: i

  string_in_arr = .false.
  do i=1,narr
    if( trim(string).eq.trim(array(i)) ) then
      string_in_arr = .true.
      return
    endif
  enddo
  return
  
end function string_in_arr
!-----------------------------------------------------------------------
! Local Variables:
! compile-command: "make fitpot"
! End:
