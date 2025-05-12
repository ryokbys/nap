program fitpot
!-----------------------------------------------------------------------
!                     Last modified: <2025-05-10 11:04:27 KOBAYASHI Ryo>
!-----------------------------------------------------------------------
  use variables
  use parallel
!!$  use NNd,only:NN_init,NN_func,NN_grad
  use fp_common,only: func_w_pmd, grad_w_pmd, write_dsgnmats &
       ,subtract_FF, restore_FF, normalize, wrap_ranges, init_fp_common
  use composition
  use minimize
  use version
!!$  use NN2,only: set_iglid_NN2
  use linreg,only: set_iglid_linreg
  use time,only: time_stamp
  use DNN,only: write_tgrads_DNN
  use UF3,only: symmetrize_params_uf3
  use pmdvars,only: specorder_pmd => specorder
  implicit none
  integer:: ismpl,ihour,imin,isec,isp,jsp
  real(8):: tmp,ftrn0,ftst0

  call mpi_init(ierr)
  time0= mpi_wtime()
  call mpi_comm_size(mpi_comm_world,nnode,ierr)
  call mpi_comm_rank(mpi_comm_world,myid,ierr)
  mpi_world= mpi_comm_world
  tcomm= 0d0
  twait= 0d0
  terg = 0d0
  tfrc = 0d0
  tstrs = 0d0

  call init_variables()

  if( myid.eq.0 ) then
    write(6,'(a)') '========================================================================'
    write(6,'(a)') ' FITPOT: A program for FITting interatomic POTential parameters'
    write(6,*) ''
    call write_revision()
    call write_authors()
    write(6,'(a)') '========================================================================'

    write(6,*) ''
    call time_stamp(' Job started')
    write(6,*) ''
    write(6,'(a,i6)') ' Number of processes in MPI = ',nnode
    call read_infitpot(10,'in.fitpot')
!.....NN and NN2 are both pointing NN2
    if( trim(cpot).eq.'NN' .or. trim(cpot).eq.'NN2' ) then
      print *,'ERROR: NN and NN2 potentials are no longer available in fitpot.'
      print *,'       Use DNN instead.'
      stop
    endif
!.....Check GDW; GDW works only with ML potentials, which use descriptors
    if( index(cpot,'NN').eq.0  .and. trim(cpot).ne.'linreg' ) then
      if( lgdw ) print *,'Gaussian density weight only works for ML potentials, so unset GDW.'
      lgdw = .false.
    endif
!!$!.....Normalization and GDW cannot be used with SGD
!!$    if( trim(cfmethod).eq.'sgd' .or. trim(cfmethod).eq.'SGD' ) then
!!$      if( lgdw ) print *,'Gaussian density weight is not available with SGD, so unset GDW.'
!!$      lgdw = .false.
!!$      if( lnormalize ) print *,'Normalization is not available with SGD, so unset normalization.'
!!$      lnormalize = .false.
!!$    endif
!.....GDW assumes that G's are normalized
    if( lgdw ) lnormalize = .true.
    call write_initial_setting()
    call normalize_weights()
  endif
  call sync_input()

!.....Copy specorder in fitpot to that in pmdvars
  specorder_pmd(:) = specorder(:)

!!$!.....read_params_desc in read_params.F90
!!$  if( index(cpot,'NN').ne.0 .or. trim(cpot).eq.'linreg' ) call read_params_desc()

  call read_vars()
  allocate(gvar(nvars),dvar(nvars))
  dmem = dmem +8d0*size(gvar) +8d0*size(dvar)
!!$  if( myid.eq.0 ) print *,'after read_vars, dmem=',dmem

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

  call get_smpl_list(11)
!.....store dirname
  do ismpl=isid0,isid1
    samples(ismpl)%csmplname= csmplist(ismpl)
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
!!$  call read_ref_data()
  if( lwgt_compos ) then
    call read_compos_list()
    call assign_compos_weight()
  endif
  call get_base_energies()
  call set_max_num_atoms()

  if( nswgt.gt.0 ) call set_sample_weights()

!!$!...Probably this is obsolete, since subtracting other forces should be done as preprocess
!!$!.....Subtract energy and forces of other force-fields
!!$  if( nsubff.gt.0 ) then
!!$    call subtract_FF()
!!$  endif

!.....Compute data statistics after reading ref samples and
!     subtracting FFs if needed.
  call get_data_stats()

!!$!.....Set cffs only for pmd calculation
!!$  if( trim(cpot).eq.'BVS' ) then
!!$    nff = 2
!!$    allocate(cffs(nff))
!!$    cffs(1) = 'Coulomb'
!!$    cffs(2) = 'Morse'
!!$  else if( trim(cpot).eq.'BVSx' ) then
!!$    nff = 3
!!$    allocate(cffs(nff))
!!$    cffs(1) = 'Coulomb'
!!$    cffs(2) = 'Morse'
!!$    cffs(3) = 'angular'
!!$  else if( trim(cpot).eq.'fpc' ) then
!!$    nff = 2
!!$    allocate(cffs(nff))
!!$    cffs(1) = 'fpc'
!!$    cffs(2) = 'Coulomb'
!!$  else
!!$    nff = 1
!!$    allocate(cffs(nff))
!!$    cffs(1) = trim(cpot)
!!$  endif
    nff = 1
    allocate(cffs(nff))
    cffs(1) = trim(cpot)

  if( (trim(cfmethod).ne.'test' .or. trim(cfmethod).ne.'dsgnmat') .and. &
       trim(cpot).eq.'linreg' .or. trim(cpot).eq.'dnn' ) then
    lnormalize = .true.
  endif

  call init_fp_common()

!.....Initial computations of all samples
  if( trim(cpot).eq.'linreg' .or. trim(cpot).eq.'dnn' &
       .or. cpot(1:3).eq.'uf3' ) then
!.....Some restriction to parameters in case of UF3 potential.
!.....No need for UF3L.
    if( trim(cpot)=='uf3' ) call symmetrize_params_uf3(nvars,vars)
    call wrap_ranges(nvars,vars,vranges)
    call func_w_pmd(nvars,vars,ftrn0,ftst0)
  else
    print *,'ERROR: '//trim(cpot)//' is not available.'
    stop
  endif

!.....NN2 only, and should be called after func_w_pmd
!!$  if( trim(cpot).eq.'NN2' ) then
!!$    call set_iglid_NN2(cpena,cfmethod)
!!$  else if( trim(cpot).eq.'linreg' ) then
  if( trim(cpot).eq.'linreg' ) then
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
  case ('lbfgs','L-BFGS','LBFGS','bfgs','BFGS')
    call qn_wrapper(ftrn0,ftst0)
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

!.....Compute func value of the best vars
  call func_w_pmd(nvars,vbest,ftrn0,ftst0)
  if( myid.eq.0 ) then
    print '(a,i0,f8.4)', ' Best of test loss (iter,loss): ', &
         ibest,ftst0
  endif
  call write_stats(ibest)

!!$  call write_energy_relation('subtracted')
!!$  if( nsmpl.lt.nsmpl_outfrc ) then
!!$    call write_force_relation('subtracted')
!!$  endif


!!$!.....Restore subtracted energies and forces to get original reference values
!!$  if( nsubff.gt.0 ) then
!!$    call restore_FF()
!!$  endif

  call write_vars(nvars,vbest,vranges,'best')
  call write_vars(nvars,vars,vranges,'fin')
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
    if( dmem > 1d+9 ) then
      write(6,'(a,f19.3,a)') ' Memory/proc = ',dmem/1d+9,' GB'
      write(6,'(a,f17.3,a)') ' Memory(total) = ',nnode*dmem/1d+9,' GB'
    else if( dmem > 1d+6 ) then
      write(6,'(a,f19.1,a)') ' Memory/proc = ',dmem/1d+6,' MB'
      write(6,'(a,f17.1,a)') ' Memory(total) = ',nnode*dmem/1d+6,' MB'
    else
      write(6,'(a,f19.1,a)') ' Memory/proc = ',dmem/1d+3,' kB'
      write(6,'(a,f17.1,a)') ' Memory(total) = ',nnode*dmem/1d+3,' kB'
    endif
    write(6,'(a,f15.3,a)') ' Time func (max) = ', tfunc,' sec'
    write(6,'(a,f15.3,a)') ' Time grad (max) = ', tgrad,' sec'
    write(6,'(a,f15.3,a)') ' Time comm (max) = ', tcomm,' sec'
    write(6,'(a,f15.3,a)') ' Time wait (max) = ', twait,' sec'
    tmp = mpi_wtime() -time0
    ihour = int(tmp/3600)
    imin  = int((tmp-ihour*3600)/60)
    isec  = int(tmp -ihour*3600 -imin*60)
    write(6,'(a,f15.3,a,i3,"h",i2.2,"m",i2.2,"s")') &
         ' Time      = ', tmp, &
         ' sec  = ', ihour,imin,isec
    if( iprint.gt.1 ) then
      print '(a,3(2x,f0.3))',' Time for erg,frc,strs in xxx_w_pmd = ', &
           terg,tfrc,tstrs
      call write_tgrads_DNN(myid)
    endif
    call time_stamp(' Job finished')
  endif
  call mpi_finalize(ierr)

end program fitpot
!=======================================================================
subroutine write_initial_setting()
  use variables
  use minimize
  use random
  use composition
  implicit none 
  integer:: i,isp,jsp

  write(6,*) ''
  write(6,'(a)') '------------------------------------------------------------------------'
  write(6,'(a)') '                          Initial setting                               '
  write(6,'(a)') '------------------------------------------------------------------------'
  write(6,'(2x,a25,2x,i0)') 'num_samples',nsmpl
  write(6,'(2x,a25,2x,i0)') 'num_iteration',niter
  write(6,'(2x,a25,2x,f5.2)') 'test_ratio',ratio_test
  
  print *,''
!!$  write(6,'(2x,a25,2x,a)') 'dataset_directory',trim(cdatasetdir)
  write(6,'(2x,a25,2x,a)') 'sample_file',trim(csmplfile)
  write(6,'(2x,a25,2x,a)') 'param_file',trim(cparfile)
!!$  write(6,'(2x,a25,2x,a)') 'run_mode',trim(crunmode)
  write(6,'(2x,a25,10(2x,a3))') 'specorder',(trim(specorder(i)),i=1,nsp)

  write(6,'(a)') ''
  write(6,'(2x,a25,2x,a)') 'potential',trim(cpot)
  if( nsubff.gt.0 ) then
    write(6,'(2x,a25,10(2x,a))') 'subtract_force_field',(trim(csubffs(i)),i=1,nsubff)
  endif
  print *,''
  write(6,'(2x,a25,2x,a)') 'fitting_method',trim(cfmethod)
  write(6,'(2x,a25,2x,a)') 'loss_function',trim(ctype_loss)
  write(6,'(2x,a25,2x,es12.3)') 'xtol',xtol
  write(6,'(2x,a25,2x,es12.3)') 'ftol',ftol
  write(6,'(2x,a25,2x,es12.3)') 'gtol',gtol
  write(6,'(2x,a25,2x,i0)') 'numtol',numtol
  write(6,'(2x,a25,2x,a)') 'init_params',trim(cinitv)
  write(6,'(2x,a25,2x,f0.2)') 'init_params_rs',vinitrs
  write(6,'(2x,a25,2x,f0.2)') 'init_params_sgm',vinitsgm
  write(6,'(2x,a25,2x,f0.2)') 'init_params_mu',vinitmu
  
  print *,''
  write(6,'(2x,a25,2x,l3)') 'energy_match',lematch
  write(6,'(2x,a25,2x,l3)') 'force_match',lfmatch
  write(6,'(2x,a25,2x,l3)') 'stress_match',lsmatch
  print *,''
  if( lfmatch ) then
    write(6,'(2x,a25,2x,a)') 'force_denom_type',trim(cfrc_denom)
    write(6,'(2x,a25,2x,es12.3)') 'force_limit',force_limit
    if( rate_eval_frc < 1d0 -1d-8 ) write(6,'(2x,a25,2x,f6.3)') &
         'reduce_fmatch',rate_eval_frc
  endif
  if( lsmatch ) then
    write(6,'(2x,a25,2x,a)') 'stress_denom_type',trim(cstrs_denom)
    write(6,'(2x,a25,2x,es12.3)') 'stress_limit',stress_limit
  endif
  write(6,'(2x,a25,2x,es12.3)') 'fval_upper_limit',fupper_lim
  
  write(6,'(a)') ''
  if( trim(cpenalty).ne.'none' ) write(6,'(2x,a25,2x,a)') 'penalty',trim(cpenalty)
  if( trim(cpenalty).eq.'ridge' ) then
    write(6,'(2x,a25,2x,es12.3)') 'penalty_weight',penalty
  else if( trim(cpenalty).eq.'uf3' ) then
    write(6,'(2x,a25,2x,es12.3)') 'pwgt_2b',pwgt2b
    write(6,'(2x,a25,2x,es12.3)') 'pwgt_2b_diff',pwgt2bd
    write(6,'(2x,a25,2x,es12.3)') 'pwgt_2b_short',pwgt2bs
    write(6,'(2x,a25,2x,es12.3)') 'pwgt_3b',pwgt3b
    write(6,'(2x,a25,2x,es12.3)') 'pwgt_3b_diff',pwgt3bd
  endif
!!$  write(6,'(2x,a25,2x,l3)') 'gradient',lgrad
!!$  write(6,'(2x,a25,2x,l3)') 'grad_scale',lgscale
!!$  write(6,'(2x,a25,2x,es12.3)') 'gscale_factor',gscl
  write(6,'(2x,a25,2x,a)') 'normalize_input',trim(cnormalize)

  if( nspcs_neglect.gt.0 ) then
    write(6,'(2x,a25,9(2x,4a))') 'force_neglect_species',(cspcs_neglect(i),i=1,nspcs_neglect)
  endif

  if( lwgt_compos ) then
    print '(2x,a25,2x,l3)','compos_weight',lwgt_compos
    print '(2x,a25,2x,f6.2)','compos_weight_scale',escl_compos
  endif

  if( lgdw ) then
    print *,''
    write(6,'(2x,a25,2x,l3)') 'gaussian_density_weight',lgdw
    write(6,'(2x,a25,2x,es12.3)') 'GDW_sigma',gdsgm
  endif
  if( nswgt > 0 ) then
    write(6,'(2x,a25,2x,i0)') 'sample_weight',nswgt
    do i=1,nswgt
      write(6,'(4x,a23,1x,f10.4)') trim(cswgt(i)), swgt0(i)
    enddo
  endif
  
!!$  write(6,'(2x,a25,2x,es12.3)') 'coeff_sequential',seqcoef
  write(6,'(2x,a25,2x,a)') 'line_minimization',trim(clinmin)
  write(6,'(a)') ''
  if( nserr > 0 ) then
    write(6,'(2x,a25,2x,i0)') 'sample_error',nserr
    do i=1,nserr
      write(6,'(4x,a23,3(1x,f10.4))') trim(cserr(i)), seerr(i), sferr(i), sserr(i)
    enddo
  endif
  if( len(trim(crefstrct)).gt.5 ) then
    print *,''
    write(6,'(2x,a25,2x,a)') 'reference_structure',trim(crefstrct)
  else if( trim(cpot).ne.'dnn' ) then
    do i=1,nspmax
      if( trim(specorder(i)).ne.'x' ) then
        write(6,'(2x,a25,2x,i2,a4,es15.7)') 'atom_energy',i,specorder(i),eatom(i)
      endif
    enddo
  endif
  write(6,'(a)') ''

  if( l_correct_short ) then
    write(6,'(2x,a25,2x,l3)') 'correct_short', l_correct_short
    do isp=1,nsp
      do jsp=isp,nsp
        write(6,'(2x,a25,2x,a,f7.3)') 'short_radii', &
             trim(specorder(isp))//' '//trim(specorder(jsp)), &
             short_radii(isp,jsp)
      enddo
    enddo
    write(6,'(a)') ''
  endif
  
  write(6,'(a)') '------------------------------------------------------------------------'

end subroutine write_initial_setting
!=======================================================================
subroutine get_smpl_list(ionum)
  use variables
  use parallel
  use fp_common,only: ndat_in_line
  use util,only: num_data
  implicit none
  integer,intent(in):: ionum
  integer:: is,ndat,ipos
  logical:: lerror = .false.
!!$  integer,external:: num_data
  character(len=128):: ctmp
  logical:: luse_iclist

  if( .not. allocated(csmplist)) allocate(csmplist(nsmpl))
  if( .not. allocated(iclist)) allocate(iclist(nsmpl))

  if( myid.eq.0 ) then
    lerror = .true.
    if( len(trim(csmplistfile)).lt.1 ) then
      print '(/,a)',' Creating sample list by performing the following command:'
      ctmp = 'ls '//trim(csmplfile)//' > smpl_list.txt'
      print '(a)','   $ '//trim(ctmp)
      call system(trim(ctmp))
!!$      print *,'  $ ls '//trim(cdatasetdir) &
!!$           //' | grep "smpl_" > smpl_list.txt'
!!$      call system('ls '//trim(cdatasetdir) &
!!$           //' | grep "smpl_" > smpl_list.txt')
      open(ionum,file='smpl_list.txt',status='old')
    else
      print *,'Sample list was given by input.'
      open(ionum,file=trim(csmplistfile),status='old')
    endif
    read(ionum,'(a)',end=998) ctmp
    backspace(ionum)
    ndat = num_data(trim(ctmp),' ')
    if( ndat.eq.1 ) then
      do is=1,nsmpl
        read(ionum,'(a)',end=998) ctmp
        csmplist(is) = trim(ctmp)
      enddo
      lerror = .false.
      luse_iclist = .false.
      call shuffle_list(nsmpl,csmplist,iclist,luse_iclist)
    else if(ndat.eq.2 ) then
      print *,'training and test are determined by input, ',trim(csmplistfile)
      do is=1,nsmpl
        read(ionum,'(a)',end=998) ctmp
        ipos = scan(ctmp, ' ')  ! get space position
        csmplist(is) = ctmp(:ipos)  ! get 1st position word as sample file
        backspace(ionum)  ! to read this line again, go back
        read(ionum,*,end=998) ctmp, iclist(is)  ! get iclist(is) only
      enddo
      lerror = .false.
      luse_iclist = .true.
      call shuffle_list(nsmpl,csmplist,iclist,luse_iclist)
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
  call mpi_bcast(csmplist,128*nsmpl,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(ndat,1,mpi_integer,0,mpi_world,ierr)
  if( ndat.eq.2 ) then
    call mpi_bcast(iclist,nsmpl,mpi_integer,0,mpi_world,ierr)
    test_assigned = .true.
  endif

!!$  if( myid.eq.0 ) then
!!$    do is=1,nsmpl
!!$      print *,' is,csmpl=',is,csmplist(is)
!!$    enddo
!!$  endif
  
  if(myid.eq.0 .and. iprint.gt.1 ) print*,'Finished get_smpl_list'
  return

999 continue
  if( myid.eq.0 ) print *,' ERROR@get_smp_list: num_samples may be wrong.'
  call mpi_finalize(ierr)
  stop

end subroutine get_smpl_list
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
!!$    print *,''
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
!
!  Read sample files of the new pmd format, in which erg, frc, strs data
!  are stored not as separate files.
!
  use variables
  use parallel
  implicit none

  integer:: is,isp,jsp,ia,ispmax,ispmaxl,nfrc,nfrcg,nftot,nftotg
  character:: cdir*128, cspmd*3, cspfp*3
  character:: cfname*128
  integer,allocatable:: nal(:)

  if( .not. allocated(nalist) ) allocate(nalist(nsmpl))
  allocate(nal(nsmpl))
  nalist(1:nsmpl)= 0
  nal(1:nsmpl)= 0
  nftot = 0
  nfrc = 0
  ispmaxl = 0
  
  do is=isid0,isid1
    cfname= samples(is)%csmplname
!!$    call read_smpl(12,trim(cdatasetdir)//'/'//trim(cfname),is,samples(is))
    call read_smpl(12,trim(cfname),is,samples(is))
    nal(is)= samples(is)%natm
!.....Num of each species
    samples(is)%naps(1:nspmax) = 0
    ispmax = 0
    do ia=1,samples(is)%natm
      isp = int(samples(is)%tag(ia))
      ispmax = max(ispmax,isp)
      samples(is)%naps(isp) = samples(is)%naps(isp) +1
    enddo
    ispmaxl = max(ispmaxl,ispmax)
    samples(is)%ispmax = ispmax
!.....Num of forces to be calculated
    nfrc = nfrc +samples(is)%nfcal*3
    nftot = nftot +samples(is)%natm*3
!.....Specorder in fitpot and that in sample should be the same
    do isp=1,nspmax
      cspmd = samples(is)%specorder(isp)
      cspfp = specorder(isp)
      if( trim(cspfp).ne.trim(cspmd) .and. trim(cspmd).ne.'x' ) then
        print '(a)','ERROR: specorder in the sample is different from that in in.fitpot.'
        print '(a,2a5)','   failed species in fitpot and the sample: ',trim(cspfp),trim(cspmd)
        print '(a)','   specorder in fitpot and the sample, '//trim(cdir)
        do jsp=1,nspmax
          print '(i5,2a5)',jsp,specorder(jsp), samples(is)%specorder(jsp)
        enddo
        stop
      endif
    enddo
  enddo

  nfrcg = 0
  nftotg = 0
  maxisp = 0
  call mpi_reduce(nfrc,nfrcg,1,mpi_integer,mpi_sum,0,mpi_world,ierr)
  call mpi_reduce(nftot,nftotg,1,mpi_integer,mpi_sum,0,mpi_world,ierr)
  call mpi_allreduce(ispmaxl,maxisp,1,mpi_integer,mpi_max,mpi_world,ierr)

  call mpi_reduce(nal,nalist,nsmpl,mpi_integer,mpi_sum &
       ,0,mpi_world,ierr)

  call mpi_barrier(mpi_world,ierr)
  deallocate(nal)
  
  if( myid.eq.0 ) then
    print *,''
    if( iprint.gt.1 ) print '(a)',' Finished read_smpl'
    if( lfmatch ) then
      write(6,'(a,i0)') ' Number of forces to be used = ',nfrcg
      write(6,'(a,i0)') ' Total number of forces      = ',nftotg
    endif
    write(6,'(a,i0)') ' Number of species in all samples = ',maxisp
  endif
  return
end subroutine read_samples
!=======================================================================
subroutine read_smpl(ionum,fname,ismpl,smpl)
!
!  Read a sample from smpl_XXX file which is a new format since 2024-11-07
!  and contains energy, forces, stress data in it.
!
  use variables
  use util,only: num_data
  use random,only: urnd
  implicit none 
  integer,intent(in):: ionum,ismpl
  character(len=*),intent(in):: fname
  type(mdsys),intent(inout):: smpl

  integer:: i,natm,num,ia,l
  real(8):: tmp,stmp(3,3),ftmp(3)
  character(len=128):: cline
  character(len=10):: c1,copt,ctmp1,ctmp2,ctmp3,ctmp
  logical:: ltmp

  open(ionum,file=trim(fname),status='old')
!.....First, parse options
  do while(.true.)
    read(ionum,'(a)') cline
    if( cline(1:1).eq.'!' .or. cline(1:1).eq.'#' ) then
      if( index(cline,'specorder:').ne.0 ) then
        num = num_data(trim(cline),' ')
        if( num.gt.11 ) stop 'ERROR@read_smpl: number of species exceeds the limit.'
        read(cline,*) c1, copt, smpl%specorder(1:num-2)
!!$        print *,'specorder = ',smpl%specorder(1:num-2)
      else if( (index(cline,'energy:').ne.0 .or. &
           index(cline,'potential_energy:').ne.0 ) .and. &
           .not.index(cline,'kinetic').ne.0 ) then
        read(cline,*) c1, copt, smpl%eref
        smpl%leref_given = .true.
      else if( index(cline,'stress:').ne.0 ) then
        num = num_data(trim(cline),' ')
        if( num.lt.8 ) stop 'ERROR@read_smpl: number of data'&
             //' for stress is not enough, it must be more than 6.'
        read(cline,*) c1, copt, stmp(1,1),stmp(2,2),stmp(3,3), &
             stmp(2,3),stmp(1,3),stmp(1,2)
        stmp(2,1) = stmp(1,2)
        stmp(3,1) = stmp(1,3)
        stmp(3,2) = stmp(2,3)
        smpl%sref(1:3,1:3) = stmp(1:3,1:3)
        smpl%lsref_given = .true.
      else if( index(cline,'auxiliary_data:').ne.0 .or. &
           index(cline,'aux_data:').ne.0 ) then
!.....Assuming that auxiliary_data is "fx fy fz"
        read(cline,*) c1, copt, ctmp1, ctmp2, ctmp3
        if( trim(ctmp1).eq.'fx' .and. trim(ctmp2).eq.'fy' &
             .and. trim(ctmp3).eq.'fz' ) smpl%lfref_given = .true.
      endif
    else
      backspace(ionum)
      exit
    endif
  enddo
!.....Read cell information
  smpl%h(:,:,:) = 0d0
  read(ionum,*) smpl%h0
  read(ionum,*) ((smpl%h(ia,1,l),ia=1,3),l=0,1)
  read(ionum,*) ((smpl%h(ia,2,l),ia=1,3),l=0,1)
  read(ionum,*) ((smpl%h(ia,3,l),ia=1,3),l=0,1)
!!$  read(ionum,*) tmp,tmp,tmp
!!$  read(ionum,*) tmp,tmp,tmp
!!$  read(ionum,*) tmp,tmp,tmp
  smpl%h(1:3,1:3,0) = smpl%h(1:3,1:3,0)*smpl%h0
!.....Read num of atoms
  read(ionum,*) natm
  smpl%natm= natm
!.....Allocate arrays with length of num of atoms.
  allocate(smpl%ra(3,natm),smpl%fa(3,natm) &
       ,smpl%tag(natm) &
       ,smpl%fref(3,natm), smpl%fabs(natm) &
       ,smpl%va(3,natm),smpl%strsi(3,3,natm) &
       ,smpl%eki(3,3,natm),smpl%epi(natm) &
       ,smpl%fsub(3,natm),smpl%eatm(natm), &
       smpl%lfrc_eval(natm))
  dmem = dmem +8d0*size(smpl%ra) +8d0*size(smpl%fa) +8d0*size(smpl%tag) &
       +8d0*size(smpl%fref) +8d0*size(smpl%fabs) &
       +8d0*size(smpl%va) +8d0*size(smpl%strsi) +8d0*size(smpl%eki) +8d0*size(smpl%epi) &
       +8d0*size(smpl%fsub) +8d0*size(smpl%eatm) +4d0*size(smpl%lfrc_eval)
  if( lgdw ) then
    allocate(smpl%gdf(natm),smpl%gdw(natm))
    dmem = dmem +8d0*size(smpl%gdf) +8d0*size(smpl%gdw)
  endif
  smpl%esub= 0d0
  smpl%fsub(1:3,1:natm)= 0d0
  smpl%ssub(1:3,1:3) = 0d0
  if( smpl%lfref_given ) then
!.....NOTE: unit of original forces in smpl files is eV/A/A (scaled by h-mat),
!           thus we need to revert them to eV/A unit.
    do i=1,smpl%natm
      read(ionum,*) smpl%tag(i),smpl%ra(1:3,i), tmp,tmp,tmp,ftmp(1:3)
      smpl%fref(1:3,i) = smpl%h(1:3,1,0)*ftmp(1) &
           +smpl%h(1:3,2,0)*ftmp(2) +smpl%h(1:3,3,0)*ftmp(3)
      smpl%fabs(i) = sqrt(smpl%fref(1,i)**2 +smpl%fref(2,i)**2 &
           +smpl%fref(3,i)**2)
    enddo
!.....Limit number of forces to be evaluated if rate_eval_frc < 1.0
    smpl%lfrc_eval(:) = .true.
    if( rate_eval_frc < 1d0 -1d-8 &
         .and. rate_eval_frc > 0d0 &
         .and. lfmatch ) then
      do i=1,smpl%natm
        if( urnd() > rate_eval_frc ) smpl%lfrc_eval(i) = .false.
      enddo
!.....if all the atoms are excluded, set 1st atom to true
      ltmp = .false.
      do i=1,smpl%natm
        ltmp = smpl%lfrc_eval(i)
      enddo
      if( .not. ltmp ) smpl%lfrc_eval(1)= .true.
    endif
!.....Count num of force to be calculated
    smpl%nfcal = 0
    do i=1,smpl%natm
      if( smpl%lfrc_eval(i) ) smpl%nfcal = smpl%nfcal +1
    enddo
  else  ! if forces are not given (does this happen?)
    do i=1,smpl%natm
      read(ionum,*) smpl%tag(i),smpl%ra(1:3,i), tmp,tmp,tmp
    enddo
  endif
  close(ionum)
  
end subroutine read_smpl
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

#ifdef __DEBUG__
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
  implicit none
  real(8),intent(in):: ftrn0,ftst0

!!$  if( trim(cpot).eq.'Morse' .or. trim(cpot).eq.'BVS' &
!!$       .or. trim(cpot).eq.'linreg' .or. trim(cpot).eq.'dnn' ) then
  if( trim(cpot).eq.'linreg' .or. trim(cpot).eq.'dnn' &
       .or. cpot(1:3).eq.'uf3' ) then
    call qn(nvars,vars,vbest,ibest,fbest,gvar,dvar,vranges,xtol,gtol,ftol,niter &
         ,iprint,iflag,myid,cfmethod &
         ,niter_eval)
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
  use parallel
  use minimize
  implicit none
  real(8),intent(in):: ftrn0,ftst0

  call steepest_descent(nvars,vars,vbest,ibest,fbest,gvar,dvar,vranges,xtol,gtol &
       ,ftol,niter,iprint,iflag,myid,cfmethod &
       ,niter_eval)

  return
end subroutine sd_wrapper
!=======================================================================
subroutine cg_wrapper(ftrn0,ftst0)
  use variables
!!$  use NNd,only:NN_init,NN_func,NN_grad,NN_restore_standard,NN_analyze
  use parallel
  use minimize
  implicit none
  real(8),intent(in):: ftrn0,ftst0

!!$  if( trim(cpot).eq.'Morse' .or. trim(cpot).eq.'BVS' &
!!$       .or. trim(cpot).eq.'linreg' .or. trim(cpot).eq.'dnn' ) then
  if( trim(cpot).eq.'linreg' .or. trim(cpot).eq.'dnn' &
       .or. cpot(1:3).eq.'uf3' ) then
    call cg(nvars,vars,vbest,ibest,fbest,gvar,dvar,vranges,xtol,gtol,ftol,niter &
         ,iprint,iflag,myid,cfmethod &
         ,niter_eval)
  else
    if( myid.eq.0 ) then
      print *,'Warning: CG is not available for '&
           //trim(cpot)
    endif
  endif
  
  return
end subroutine cg_wrapper
!=======================================================================
subroutine sgd_wrapper(ftrn0,ftst0)
!
! Wrapper for SGD
!
  use variables
  use minimize
  use parallel
  implicit none
  real(8),intent(in):: ftrn0,ftst0

  call sgd(nvars,vars,vbest,ibest,fbest,gvar,dvar,vranges,xtol,gtol,ftol,niter, &
       iprint,iflag,myid,mpi_world,mynsmpl,myntrn,isid0,isid1, &
       cfmethod,niter_eval)

end subroutine sgd_wrapper
!=======================================================================
subroutine check_grad(ftrn0,ftst0)
  use variables
!!$  use NNd,only:NN_init,NN_func,NN_grad
  use parallel
  use fp_common,only: func_w_pmd, grad_w_pmd, wrap_ranges
  implicit none
  real(8),intent(in):: ftrn0,ftst0
  integer:: iv
  real(8):: dv,vmax,ftst,ftmp1,ftmp2, absgnum, absgana
  real(8),allocatable:: ganal(:),gnumer(:),vars0(:)
  real(8),parameter:: dev  = 1d-5
  real(8),parameter:: tiny = 1d-8
  real(8):: vtmp1,vtmp2

  allocate(gnumer(nvars),ganal(nvars),vars0(nvars))

  vars0(1:nvars)= vars(1:nvars)
  vmax= 0d0
  do iv=1,nvars
    vmax= max(vmax,abs(vars0(iv)))
  enddo
  dv= vmax *dev

  if( myid.eq.0 ) then
    print *,''
    write(6,'(a)') '------------------------------ check_grad '&
         //'------------------------------'
    print '(a,es12.4)',' Deviation for numerical derivative =',dv
    print *,''
    write(6,'(a)') '     #,          x,    analytical,'// &
         '     numerical,'// &
         '      error [%]'
  endif
  call flush(6)
  
  call wrap_ranges(nvars,vars,vranges)
  call grad_w_pmd(nvars,vars,ganal)

!.....Loop over variables for numerical derivative
  do iv=1,nvars
    if( id_check_grad > 0 .and. id_check_grad <= nvars &
         .and. iv.ne.id_check_grad ) cycle
    vars(1:nvars)= vars0(1:nvars)
    dv = max(abs(vars(iv)*dev),dev)
    vars(iv)= vars(iv) +dv/2
    call wrap_ranges(nvars,vars,vranges)
    call func_w_pmd(nvars,vars,ftmp1,ftst)
    vtmp1 = vars(iv)
    vars(1:nvars)= vars0(1:nvars)
    vars(iv)= vars(iv) -dv/2
    call wrap_ranges(nvars,vars,vranges)
    call func_w_pmd(nvars,vars,ftmp2,ftst)
    gnumer(iv)= (ftmp1-ftmp2)/dv
    vtmp2 = vars(iv)
    if( myid.eq.0 ) then
      absgnum = abs(gnumer(iv))
      absgana = abs(ganal(iv))
      write(6,'(i6,es12.4,2es15.4,f15.3)') iv,vars0(iv), &
           ganal(iv) ,gnumer(iv), &
           abs(ganal(iv)-gnumer(iv)) &
           /(max(max(absgnum,absgana),tiny)) *100
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
  if( trim(cpot).eq.'linreg' .or. index(cpot,'NN').ne.0 ) then
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
    dmem = dmem +8d0*size(erefl) +8d0*size(erefg) +8d0*size(epotl) +8d0*size(epotg) &
         +8d0*size(eerrl) +8d0*size(eerrg) +8d0*size(swgtl) +8d0*size(swgtg) &
         +8d0*size(esubl) +8d0*size(esubg)
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
    write(90,'(a)') '# eref, epot, smplname, diff, error, esub, swgt'
    write(91,'(a)') '# eref, epot, smplname, diff, error, esub, swgt'
    do ismpl=1,nsmpl
      epotg(ismpl)= epotg(ismpl)/nalist(ismpl)
      if( iclist(ismpl).eq.1 ) then
        write(90,'(2es15.7,2x,a,10es13.4e3)') erefg(ismpl) &
             ,epotg(ismpl),trim(csmplist(ismpl)) &
             ,abs(erefg(ismpl)-epotg(ismpl)) &
             ,eerrg(ismpl),esubg(ismpl),swgtg(ismpl)
      else if( iclist(ismpl).eq.2 ) then
        write(91,'(2es15.7,2x,a,10es12.3e3)') erefg(ismpl) &
             ,epotg(ismpl),trim(csmplist(ismpl)) &
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
  real(8):: dmeml
  logical:: l1st = .true.
  
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
         ,fsubg(3,nmax,nsmpl) &
         ,gdwl(nmax,nsmpl),gdwg(nmax,nsmpl))
    dmeml = 8d0*size(frefl) +8d0*size(frefg) +8d0*size(fal) +8d0*size(fag) &
         +8d0*size(ferrl) +8d0*size(ferrg) +8d0*size(fsubl) +8d0*size(fsubg) &
         +8d0*size(gdwl) +8d0*size(gdwg)
    dmem = dmem +dmeml
    if( iprint.gt.1 .and. myid.eq.0 .and. l1st ) then
      print '(a,f0.3,a)',' Memory for write_force_relation = ', &
           dmeml/1000/1000,' MB'
    endif
  endif

  if( l1st ) then
    frefl(1:3,1:nmax,1:nsmpl)= 0d0
    fsubl(1:3,1:nmax,1:nsmpl)= 0d0
    ferrl(1:nsmpl) = 0d0
    gdwl(:,:) = 0d0
    do ismpl=isid0,isid1
      natm= samples(ismpl)%natm
      frefl(1:3,1:natm,ismpl)= samples(ismpl)%fref(1:3,1:natm)
      fsubl(1:3,1:natm,ismpl)= samples(ismpl)%fsub(1:3,1:natm)
      ferrl(ismpl) = samples(ismpl)%ferr
      if( lgdw ) gdwl(1:natm,ismpl) = samples(ismpl)%gdw(1:natm)
    enddo
    frefg(1:3,1:nmax,1:nsmpl)= 0d0
    fsubg(1:3,1:nmax,1:nsmpl)= 0d0
    ferrg(1:nsmpl) = 0d0
    gdwg(:,:) = 0d0
    call mpi_reduce(frefl,frefg,3*nmax*nsmpl,mpi_real8,mpi_sum &
         ,0,mpi_world,ierr)
    call mpi_reduce(fsubl,fsubg,3*nmax*nsmpl,mpi_real8,mpi_sum &
         ,0,mpi_world,ierr)
    call mpi_reduce(ferrl,ferrg,nsmpl,mpi_real8,mpi_sum &
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
    write(92,'(a)') '# 1:fref, 2:fpot, 3:smplname, 4:ia, 5:ixyz, 6:diff,' &
         //' 7:error, 8:fsub, 9:gdw'
    write(93,'(a)') '# 1:fref, 2:fpot, 3:smplname, 4:ia, 5:ixyz, 6:diff,' &
         //' 7:error, 8:fsub, 9:gdw'
    do ismpl=1,nsmpl
      if( iclist(ismpl).eq.1 ) then
        natm= nalist(ismpl)
        do ia=1,natm
          do ixyz=1,3
            write(92,'(2es12.4,2x,a,i6,i3,4es11.2e3,l3)') frefg(ixyz,ia,ismpl) &
                 ,fag(ixyz,ia,ismpl) &
                 ,trim(csmplist(ismpl)),ia,ixyz &
                 ,abs(frefg(ixyz,ia,ismpl)-fag(ixyz,ia,ismpl))&
                 ,ferrg(ismpl),fsubg(ixyz,ia,ismpl),gdwg(ia,ismpl)
          enddo
        enddo
      else if( iclist(ismpl).eq.2 ) then
        natm= nalist(ismpl)
        do ia=1,natm
          do ixyz=1,3
            write(93,'(2es12.4,2x,a,i6,i3,4es11.2e3,l3)') frefg(ixyz,ia,ismpl) &
                 ,fag(ixyz,ia,ismpl) &
                 ,trim(csmplist(ismpl)),ia,ixyz &
                 ,abs(frefg(ixyz,ia,ismpl)-fag(ixyz,ia,ismpl))&
                 ,ferrg(ismpl),fsubg(ixyz,ia,ismpl),gdwg(ia,ismpl)
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
    dmem = dmem +8d0*size(srefl) +8d0*size(srefg) +8d0*size(strsl) +8d0*size(strsg) &
         +8d0*size(serrl) +8d0*size(serrg) +8d0*size(ssubl) +8d0*size(ssubg)
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
                 ,trim(csmplist(ismpl)),ixyz,jxyz &
                 ,abs(srefg(ixyz,jxyz,ismpl)-strsg(ixyz,jxyz,ismpl))&
                 ,serrg(ismpl),ssubg(ixyz,jxyz,ismpl)
          enddo
        enddo
      else if( iclist(ismpl).eq.2 ) then
        do ixyz=1,3
          do jxyz=ixyz,3
            write(95,'(2es15.6e3,2x,a,i6,i3,3es12.3e3)') srefg(ixyz,jxyz,ismpl) &
                 ,strsg(ixyz,jxyz,ismpl) &
                 ,trim(csmplist(ismpl)),ixyz,jxyz &
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
  real(8):: eref1l,eref1,eref2l,eref2,evar,tmp
  real(8):: fref1l,fref1,fref2l,fref2,fvar
  real(8):: sref1l,sref1,sref2l,sref2,svar
  real(8):: demaxl_trn,demax_trn,desuml_trn,desum_trn,rmse_trn
  real(8):: demaxl_tst,demax_tst,desuml_tst,desum_tst,rmse_tst
  real(8):: dfmaxl_trn,dfmax_trn,dfsuml_trn,dfsum_trn
  real(8):: dfmaxl_tst,dfmax_tst,dfsuml_tst,dfsum_tst
  real(8):: dsmaxl_trn,dsmax_trn,dssuml_trn,dssum_trn
  real(8):: dsmaxl_tst,dsmax_tst,dssuml_tst,dssum_tst
  real(8):: er2trn,er2tst,fr2trn,fr2tst,sr2trn,sr2tst
  real(8),save:: epotsub
  character(len=128):: cnum
  logical,save:: l1st = .true.

  if( l1st ) then
    if( myid.eq.0 ) then
      write(6,*) '# ENERGY:  ITER, TIME, ' &
           //'RMSE(TRAINING), RMSE(TEST), ' &
           //'MAX(TRAINING), MAX(TEST), ' &
           //'R^2(TRAINING), R^2(TEST)' 
      write(6,*) '# FORCE:   ITER, TIME, ' &
           //'RMSE(TRAINING), RMSE(TEST), ' &
           //'MAX(TRAINING), MAX(TEST), ' &
           //'R^2(TRAINING), R^2(TEST)'
      write(6,*) '# STRESS:  ITER, TIME, ' &
           //'RMSE(TRAINING), RMSE(TEST), ' &
           //'MAX(TRAINING), MAX(TEST), ' &
           //'R^2(TRAINING), R^2(TEST)'
    endif

!!$    call get_r2denom(etrndnm,etstdnm,ftrndnm,ftstdnm,strndnm,ststdnm)
!!$    if( myid.eq.0 .and. iprint.gt.1 ) then
!!$      print '(a,2es12.4)',' Denominator (energy) for trn, tst = ' &
!!$           ,etrndnm, etstdnm
!!$      print '(a,2es12.4)',' Denominator (force)  for trn, tst = ' &
!!$           ,ftrndnm, ftstdnm
!!$      print '(a,2es12.4)',' Denominator (stress) for trn, tst = ' &
!!$           ,strndnm, ststdnm
!!$    endif

  endif

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
  eref1l = 0d0
  eref2l = 0d0
  do ismpl=isid0,isid1
    smpl= samples(ismpl)
    natm= smpl%natm
    de= abs(smpl%epot-epotsub*natm+smpl%esub &
         -(smpl%eref-erefsub*natm))/natm
    tmp = smpl%eref/natm
    eref1l = eref1l +tmp
    eref2l = eref2l +tmp*tmp
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
  if( l1st ) then
    call mpi_reduce(eref1l,eref1,1 &
         ,mpi_real8,mpi_sum,0,mpi_world,ierr)
    call mpi_reduce(eref2l,eref2,1 &
         ,mpi_real8,mpi_sum,0,mpi_world,ierr)
    evar = eref2/(nsmpl_trn+nsmpl_tst) -(eref1/(nsmpl_trn+nsmpl_tst))**2
  endif
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
!!$  write(cnum,'(i0)') iter
!!$  call write_vars(nvars,vars,vranges,cnum)

!.....force
  dfmaxl_trn= 0d0
  dfsuml_trn= 0d0
  dfmaxl_tst= 0d0
  dfsuml_tst= 0d0
  ntrnl= 0
  ntstl= 0
  fref1l = 0d0
  fref2l = 0d0
  do ismpl=isid0,isid1
    smpl= samples(ismpl)
    nfcal= smpl%nfcal
    if( nfcal.eq.0 ) cycle
    natm= smpl%natm
    if( smpl%iclass.eq.1 ) then
      do ia=1,natm
        if( smpl%lfrc_eval(ia) ) then
          do l=1,3
            df= abs(smpl%fa(l,ia)+smpl%fsub(l,ia) -smpl%fref(l,ia))
            dfmaxl_trn= max(dfmaxl_trn,df)
            dfsuml_trn=dfsuml_trn +df*df
            ntrnl=ntrnl +1
            tmp = smpl%fref(l,ia)
            fref1l = fref1l +tmp
            fref2l = fref2l +tmp*tmp
          enddo
        else
          do l=1,3
            df= abs(smpl%fa(l,ia)+smpl%fsub(l,ia) -smpl%fref(l,ia))
            dfmaxl_tst= max(dfmaxl_tst,df)
            dfsuml_tst=dfsuml_tst +df*df
            ntstl=ntstl +1
            tmp = smpl%fref(l,ia)
            fref1l = fref1l +tmp
            fref2l = fref2l +tmp*tmp
          enddo
        endif
      enddo
    else if( smpl%iclass.eq.2 ) then
      do ia=1,natm
        do l=1,3
          df= abs(smpl%fa(l,ia)+smpl%fsub(l,ia) -smpl%fref(l,ia))
          dfmaxl_tst= max(dfmaxl_tst,df)
          dfsuml_tst=dfsuml_tst +df*df
          ntstl=ntstl +1
          tmp = smpl%fref(l,ia)
          fref1l = fref1l +tmp
          fref2l = fref2l +tmp*tmp
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
  if( l1st ) then
    call mpi_reduce(fref1l,fref1,1 &
         ,mpi_real8,mpi_sum,0,mpi_world,ierr)
    call mpi_reduce(fref2l,fref2,1 &
         ,mpi_real8,mpi_sum,0,mpi_world,ierr)
    fvar = fref2/(ntrn+ntst) -(fref1/(ntrn+ntst))**2 
  endif
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
  endif

!.....stress
  dsmaxl_trn = 0d0
  dssuml_trn = 0d0
  dsmaxl_tst = 0d0
  dssuml_tst = 0d0
  ntrnl = 0
  ntstl = 0
  sref1l = 0d0
  sref2l = 0d0
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
          tmp = smpl%sref(ixyz,jxyz)
          sref1l = sref1l +tmp
          sref2l = sref2l +tmp*tmp
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
          tmp = smpl%sref(ixyz,jxyz)
          sref1l = sref1l +tmp
          sref2l = sref2l +tmp*tmp
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
  if( l1st ) then
    call mpi_reduce(sref1l,sref1,1 &
         ,mpi_real8,mpi_sum,0,mpi_world,ierr)
    call mpi_reduce(sref2l,sref2,1 &
         ,mpi_real8,mpi_sum,0,mpi_world,ierr)
    svar = sref2/(ntrn+ntst) -(sref1/(ntrn+ntst))**2
  endif
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
!!$!.....Finally write out std of reference data
!!$    if( l1st ) print '(a,3es12.4)',' Standard deviation of reference data ' &
!!$         //'(energy, force, stress) = ', &
!!$         sqrt(evar),sqrt(fvar),sqrt(svar)
  endif

  l1st = .false.
  return
end subroutine write_stats
!=======================================================================
subroutine get_data_stats()
!
!  Compute data statistics that are used for denominators in R^2 score.
!  It is called only once after reading ref data.
!  Subtract fixed FF from the reference values (erg,frc,strs),
!  to evaluate R^2 in the range of subtracted data.
!
  use variables
  use parallel
  implicit none
  
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
  evtrn = (e2mtrn -emtrn**2)
  evtst = (e2mtst -emtst**2)
  etrndnm = evtrn *nsmpl_trn
  etstdnm = evtst *nsmpl_tst
  netrn = nsmpl_trn
  netst = nsmpl_tst
  call mpi_bcast(evtrn,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(evtst,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(etrndnm,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(etstdnm,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(netrn,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(netst,1,mpi_integer,0,mpi_world,ierr)

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
        if( smpl%lfrc_eval(ia) ) then
          do l=1,3
            tmp = smpl%fref(l,ia)
            fsumltrn = fsumltrn +tmp
            f2sumltrn= f2sumltrn +tmp*tmp
            ntrnl=ntrnl +1
          enddo
        else
          do l=1,3
            tmp = smpl%fref(l,ia)
            fsumltst = fsumltst +tmp
            f2sumltst= f2sumltst +tmp*tmp
            ntstl=ntstl +1
          enddo
        endif
      enddo
    else if( smpl%iclass.eq.2 ) then
      do ia=1,natm
        do l=1,3
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
  nftrn = 0
  nftst = 0
  call mpi_reduce(ntrnl,nftrn,1 &
       ,mpi_integer,mpi_sum,0,mpi_world,ierr)
  call mpi_reduce(ntstl,nftst,1 &
       ,mpi_integer,mpi_sum,0,mpi_world,ierr)
  fmtrn = fsumtrn/nftrn
  f2mtrn= f2sumtrn/nftrn
  if( nftst.ne.0 ) then
    fmtst = fsumtst/nftst
    f2mtst= f2sumtst/nftst
  else
    fmtst = 0d0
    f2mtst= 0d0
  endif
  if( iprint.gt.1 .and. myid.eq.0 ) then
    print *,'fmtrn,fmtrn^2,f2mtrn,ntrn =',fmtrn,fmtrn**2,f2mtrn,nftrn
    print *,'fmtst,fmtst^2,f2mtst,ntst =',fmtst,fmtst**2,f2mtst,nftst
  endif
  fvtrn = (f2mtrn -fmtrn**2)
  fvtst = (f2mtst -fmtst**2)
  ftrndnm = fvtrn *nftrn
  ftstdnm = fvtst *nftst
  call mpi_bcast(fvtrn,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(fvtst,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(ftrndnm,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(ftstdnm,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(nftrn,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(nftst,1,mpi_integer,0,mpi_world,ierr)

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
        do jxyz=ixyz,3
!!$          tmp = smpl%sref(ixyz,jxyz)-smpl%ssub(ixyz,jxyz)
          tmp = smpl%sref(ixyz,jxyz)
          ssumltrn = ssumltrn +tmp
          s2sumltrn= s2sumltrn +tmp*tmp
          ntrnl = ntrnl +1
        enddo
      enddo
    else if( smpl%iclass.eq.2 ) then
      do ixyz=1,3
        do jxyz=ixyz,3
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
  nstrn = 0
  nstst = 0
  call mpi_reduce(ntrnl,nstrn,1 &
       ,mpi_integer,mpi_sum,0,mpi_world,ierr)
  call mpi_reduce(ntstl,nstst,1 &
       ,mpi_integer,mpi_sum,0,mpi_world,ierr)
  smtrn = ssumtrn/nstrn
  s2mtrn= s2sumtrn/nstrn
  if( nstst.ne.0 ) then
    smtst = ssumtst/nstst
    s2mtst= s2sumtst/nstst
  else
    smtst = 0d0
    s2mtst= 0d0
  endif
  if( iprint.gt.1 .and. myid.eq.0 ) then
    print *,'smtrn,smtrn^2,s2mtrn,ntrn =',smtrn,smtrn**2,s2mtrn,nstrn
    print *,'smtst,smtst^2,s2mtst,ntst =',smtst,smtst**2,s2mtst,nstst
  endif
  svtrn = (s2mtrn -smtrn**2)
  svtst = (s2mtst -smtst**2)
  strndnm = svtrn *nstrn
  ststdnm = svtst *nstst
  call mpi_bcast(svtrn,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(svtst,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(strndnm,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(ststdnm,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(nstrn,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(nstst,1,mpi_integer,0,mpi_world,ierr)
  if( myid.eq.0 .and. iprint.ge.1 ) then
    print '(/a)',' Number of data (total,train,test), standard deviations (train, test):'
    print '(a,3i10,2es12.3)', '   Energy: ',netrn+netst,netrn,netst,&
         sqrt(evtrn),sqrt(evtst)
    print '(a,3i10,2es12.3)', '   Force:  ',nftrn+nftst,nftrn,nftst,&
         sqrt(fvtrn),sqrt(fvtst)
    print '(a,3i10,2es12.3)', '   Stress: ',nstrn+nstst,nstrn,nstst,&
         sqrt(svtrn),sqrt(svtst)
  endif
  
  return
end subroutine get_data_stats
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
  use pmdvars,only: nnmax
  use composition
  implicit none
  
  call mpi_bcast(nsmpl,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(niter,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(niter_eval,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(iprint,1,mpi_integer,0,mpi_world,ierr)

  call mpi_bcast(cfmethod,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(csmplfile,128,mpi_character,0,mpi_world,ierr)
!!$  call mpi_bcast(cdatasetdir,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(cparfile,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(crunmode,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(cevaltype,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(csmplistfile,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(csmplftype,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(cpot,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(cpenalty,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(clinmin,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(cfsmode,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(crefstrct,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(ctype_loss,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(cfrc_denom,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(cstrs_denom,128,mpi_character,0,mpi_world,ierr)

  call mpi_bcast(xtol,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(ftol,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(gtol,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(eatom,nspmax,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(specorder,3*nspmax,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(gscl,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(nfpsmpl,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(penalty,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(pwgt2b,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(pwgt2bd,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(pwgt2bs,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(pwgt3b,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(pwgt3bd,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(ratio_test,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(rseed,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(force_limit,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(stress_limit,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(rc_other,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(rcut,1,mpi_real8,0,mpi_world,ierr)
  
  call mpi_bcast(lematch,1,mpi_logical,0,mpi_world,ierr)
  call mpi_bcast(lfmatch,1,mpi_logical,0,mpi_world,ierr)
  call mpi_bcast(lsmatch,1,mpi_logical,0,mpi_world,ierr)
  call mpi_bcast(nspcs_neglect,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(cspcs_neglect,3*nspmax,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(interact,nspmax*nspmax,mpi_logical,0,mpi_world,ierr)
  call mpi_bcast(interact3,nspmax*nspmax*nspmax,mpi_logical,0,mpi_world,ierr)
  call mpi_bcast(wgte,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(wgtf,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(wgts,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(rate_eval_frc,1,mpi_real8,0,mpi_world,ierr)

  call mpi_bcast(fupper_lim,1,mpi_real8,0,mpi_world,ierr)
!.....sgd
  call mpi_bcast(csgdupdate,128,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(nsgdbsnode,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(sgd_rate_ini,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(sgd_rate_fin,1,mpi_real8,0,mpi_world,ierr)
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
  call mpi_bcast(m_lbfgs,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(fac_inc,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(fac_dec,1,mpi_real8,0,mpi_world,ierr)
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
  if( myid.ne.0 ) then
    allocate(cserr(nserr),seerr(nserr),sferr(nserr),sserr(nserr))
  endif
  call mpi_bcast(cserr,128*nserr,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(seerr,nserr,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(sferr,nserr,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(sserr,nserr,mpi_real8,0,mpi_world,ierr)
  
  call mpi_bcast(nswgt,1,mpi_integer,0,mpi_world,ierr)
  if( myid.ne.0 ) then
    allocate(cswgt(nswgt),swgt0(nswgt))
  endif
  call mpi_bcast(cswgt,128*nswgt,mpi_character,0,mpi_world,ierr)
  call mpi_bcast(swgt0,nswgt,mpi_real8,0,mpi_world,ierr)

  call mpi_bcast(lwgt_compos,1,mpi_logical,0,mpi_world,ierr)
  if( lwgt_compos ) then
    call mpi_bcast(escl_compos,1,mpi_real8,0,mpi_world,ierr)
  endif
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
  if( nsubff.gt.0 ) call mpi_bcast(csubffs,20*nsubff,mpi_character,0,mpi_world,ierr)

!.....Short correction
  call mpi_bcast(l_correct_short,1,mpi_logical,0,mpi_world,ierr)
  call mpi_bcast(short_radii,nspmax*nspmax,mpi_real8,0,mpi_world,ierr)
!!$!.....Repulsion correction
!!$  call mpi_bcast(n_repul_pnts,1,mpi_integer,0,mpi_world,ierr)
!!$  if( n_repul_pnts > 0 ) then
!!$    if( .not.allocated(drepul_tbl) ) &
!!$         allocate(drepul_tbl(n_repul_pnts,nspmax,nspmax))
!!$    call mpi_bcast(valence_chgs,nspmax,mpi_real8,0,mpi_world,ierr)
!!$    call mpi_bcast(core_chgs,nspmax,mpi_real8,0,mpi_world,ierr)
!!$    call mpi_bcast(pwgt_repul,1,mpi_real8,0,mpi_world,ierr)
!!$  endif

!.....pmdio related
  call mpi_bcast(nnmax,1,mpi_integer,0,mpi_world,ierr)

!.....Debugging
  call mpi_bcast(id_check_grad,1,mpi_integer,0,mpi_world,ierr)
  
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
  dmem = dmem +4d0*size(nspn) +4d0*size(ispn)
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
subroutine shuffle_list(nsmpl,clist,iclist,luse_iclist)
!
!  Randomize the order of the clist
!
  use random
  implicit none
  integer,intent(in):: nsmpl
  character(len=128),intent(inout):: clist(nsmpl)
  integer,intent(inout):: iclist(nsmpl)
  logical,intent(in):: luse_iclist 

  integer:: i,j,k,n
  character(len=128),allocatable:: cltmp(:)
  integer,allocatable:: icltmp(:)

  allocate(cltmp(nsmpl))
  if( .not. luse_iclist ) then
    cltmp(1:nsmpl)= clist(1:nsmpl)
    n=nsmpl
    do i=1,nsmpl
      j= n*urnd()+1
      clist(i)= cltmp(j)
      do k=j,n-1
        cltmp(k)= cltmp(k+1)
      enddo
      n= n -1
    enddo
  else
    allocate(icltmp(nsmpl))
    cltmp(1:nsmpl)= clist(1:nsmpl)
    icltmp(1:nsmpl)= iclist(1:nsmpl)
    n=nsmpl
    do i=1,nsmpl
      j= n*urnd()+1
      clist(i)= cltmp(j)
      iclist(i)= icltmp(j)
      do k=j,n-1
        cltmp(k)= cltmp(k+1)
        icltmp(k)= icltmp(k+1)
      enddo
      n= n -1
    enddo
    deallocate(icltmp)
  endif
  deallocate(cltmp)
end subroutine shuffle_list
!=======================================================================
subroutine normalize_weights()
  use variables,only: wgte,wgtf,wgts,lematch,lfmatch,lsmatch
  real(8):: denom

  if( .not. lematch ) wgte = 0d0
  if( .not. lfmatch ) wgtf = 0d0
  if( .not. lsmatch ) wgts = 0d0
  wgte = abs(wgte)
  wgtf = abs(wgtf)
  wgts = abs(wgts)
  denom = wgte + wgtf + wgts
  if( denom < 1d-15 ) stop 'ERROR(normalize_weights): denom < 1d-15'
  wgte = wgte /denom
  wgtf = wgtf /denom
  wgts = wgts /denom
  print '(a,3f7.3)',' Normalized weights = ',wgte,wgtf,wgts
end subroutine normalize_weights
!=======================================================================
subroutine set_sample_errors()
  use variables
  use parallel
  implicit none 
  integer:: ismpl,iserr,idx

  do ismpl=isid0,isid1
    do iserr= 1,nserr
      idx= index(samples(ismpl)%csmplname,trim(cserr(iserr)))
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
      idx = index(samples(ismpl)%csmplname,trim(cswgt(iswgt)))
      if( idx.ne.0 ) then
!!$        samples(ismpl)%wgt = exp(-(erg -swerg0(iswgt))/swdenom(iswgt))
        samples(ismpl)%wgt = swgt0(iswgt)
      endif
    enddo
!!$    write(6,'(a,i5,a20,3es12.4)') ' ismpl,cdirname,wgt,eref,erg=' &
!!$         ,ismpl,trim(samples(ismpl)%cdirname)&
!!$         ,samples(ismpl)%wgt,samples(ismpl)%eref,erg
  enddo
  if( myid.eq.0 ) then
    write(6,'(a)') ' Set sample weights provided in in.fitpot.'
  endif
end subroutine set_sample_weights
!=======================================================================
subroutine subtract_atomic_energy()
  use variables
  use parallel
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
      isp = csp2isp(trim(csp))
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
    if( trim(samples(ismpl)%csmplname).eq.trim(crefstrct) ) then
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
!
!  Set max num of atoms among all samples.
!  Also set max num of atoms whose force are to be used for force-matching.
!
  use variables
  use parallel
  integer:: ismpl, maxnal, maxninl, maxnfl

  maxnal = 0
  maxnfl = 0
  maxninl = 0
  do ismpl=isid0,isid1
    na = samples(ismpl)%natm
    maxnal = max(na,maxnal)
    maxninl = maxninl +na
    maxnfl = max(maxnfl, samples(ismpl)%nfcal)
  enddo
  maxna = 0
  maxnf = 0
  maxnin = 0
  natot = 0
  call mpi_allreduce(maxnal,maxna,1,mpi_integer,mpi_max,mpi_world,ierr)
  call mpi_allreduce(maxnfl,maxnf,1,mpi_integer,mpi_max,mpi_world,ierr)
  call mpi_allreduce(maxninl,maxnin,1,mpi_integer,mpi_max,mpi_world,ierr)
  call mpi_allreduce(maxninl,natot,1,mpi_integer,mpi_sum,mpi_world,ierr)
  if( myid.eq.0 .and. iprint.ne.0 ) then
    write(6,'(a,i0)') ' Max num of atoms in samples   = ',maxna
    write(6,'(a,i0)') ' Max num of atoms in nodes     = ',maxnin
    write(6,'(a,i0)') ' Total num of atoms in dataset = ',natot
    if( lfmatch ) write(6,'(a,i0)') ' Max num of atoms whose forces' &
         //' are used for force-matching = ',maxnf
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
