subroutine read_vars()
  use variables
  use parallel
  use random
  implicit none

!!$  if( index(cpot,'NN').ne.0 ) then
!!$    call read_vars_NN()
!!$  else
!!$    call read_vars_fitpot()
!!$  endif
!.....Always call read_vars_fitpot
  call read_vars_fitpot()

end subroutine read_vars
!=======================================================================
subroutine write_vars(nvars,vars,vranges,cadd)
  use variables,only: cnormalize, lnormalize, cparfile
  use fp_common,only: normalize, restore_normalize
  implicit none
  integer,intent(in):: nvars
  real(8),intent(in):: vars(nvars),vranges(2,nvars)
  character(len=*),intent(in):: cadd
  character(len=128):: cfname

  if( cnormalize(1:4).ne.'none' ) then
!!$    if( trim(cpot).eq.'NN' .and. .not. &
!!$         (trim(cfmethod).eq.'sa' .or. trim(cfmethod).eq.'SA' .or. &
!!$         trim(cfmethod).eq.'ga' .or. trim(cfmethod).eq.'GA' .or. &
!!$         trim(cfmethod).eq.'de' .or. trim(cfmethod).eq.'DE' .or. &
!!$         trim(cfmethod).eq.'pso' .or. trim(cfmethod).eq.'PSO') ) then
!!$      call NN_restore_standard()
!!$    else if( lnormalize ) then
    if( lnormalize ) then
      call restore_normalize()
    endif
  endif

  cfname= trim(cparfile)//'.'//trim(cadd)

!!$  if( index(cpot,'NN').ne.0 ) then
!!$    call write_vars_NN(cfname)
!!$  else
!!$    call write_vars_fitpot(cfname)
!!$  endif
!.....Always call write_vars_fitpot
  call write_vars_fitpot(nvars,vars,vranges,cfname)

  if( cnormalize(1:4).ne.'none' ) then
!!$    if( trim(cpot).eq.'NN' .and. .not. &
!!$         (trim(cfmethod).eq.'sa' .or. trim(cfmethod).eq.'SA' .or. &
!!$         trim(cfmethod).eq.'ga' .or. trim(cfmethod).eq.'GA' .or. &
!!$         trim(cfmethod).eq.'de' .or. trim(cfmethod).eq.'DE' .or. &
!!$         trim(cfmethod).eq.'pso' .or. trim(cfmethod).eq.'PSO') ) then
!!$      call NN_standardize()
!!$    else if( lnormalize ) then
    if( lnormalize ) then
      call normalize()
    endif
  endif

end subroutine write_vars
!=======================================================================
subroutine read_vars_fitpot()
!
!  Read fitpot original param-format file.
!  Linreg potenitla param file is also the same style.
!
  use variables
  use parallel
  use random
  use fp_common,only: wrap_ranges
  implicit none
  integer,parameter:: ionum = 15
  integer:: i
  real(8):: rs0
  character(len=128):: ctmp

  if( myid.eq.0 ) then
    print *,'Read parameters to be optimized from '//trim(cparfile)
    open(ionum,file=trim(cparfile),status='old')
10  read(ionum,'(a)') ctmp
    if( ctmp(1:1).eq.'!' .or. ctmp(1:1).eq.'#' ) goto 10
    read(ctmp,*) nvars, rcut, rc3
    rcut = max(rcut,rc3)
  endif
  call mpi_bcast(nvars,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(rcut,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(rc3,1,mpi_real8,0,mpi_world,ierr)
  allocate(vars(nvars),vranges(2,nvars),vbest(nvars))
  if( myid.eq.0 ) then
    print '(a,i0)',' Number of variables to be optimized = ',nvars
    do i=1,nvars
      read(ionum,*) vars(i),vranges(1:2,i)
    enddo
    if( index(cinitv,'gauss').ne.0 ) then
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
      write(6,'(a)') ' Potential parameters are read from file: '//trim(cparfile)
    endif
    close(ionum)
  endif
  call mpi_bcast(vars,nvars,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(vranges,2*nvars,mpi_real8,0,mpi_world,ierr)
  call wrap_ranges(nvars,vars,vranges)

end subroutine read_vars_fitpot
!=======================================================================
subroutine write_vars_fitpot(nvars,vars,vranges,cfname)
!
!  Write params in original fitpot format.
!
  use variables,only: rcut, rc3
  use parallel,only: myid
  implicit none
  integer,intent(in):: nvars
  real(8),intent(in):: vars(nvars), vranges(2,nvars)
  character(len=128),intent(in):: cfname
  integer,parameter:: ionum = 16
  integer:: i

  if( myid.eq.0 ) then
    open(ionum,file=trim(cfname),status='replace')
    write(ionum,'(i10,2es15.4)') nvars,rcut,rc3
    do i=1,nvars
      write(ionum,'(es24.14e3,2es14.4)') vars(i),vranges(1:2,i)
    enddo
    close(ionum)
!    print *, 'wrote '//trim(cfname)
  endif

end subroutine write_vars_fitpot
!=======================================================================
subroutine read_vars_NN()
!
!  Read in.params.NN2/DNN and convert data to vars array.
!
  use variables
  use parallel
  use random
  implicit none

  integer,parameter:: ionum = 15
  integer:: i
  real(8):: rs0
  character(len=128):: ctmp

  if( myid.eq.0 ) then
    print *,'Read parameters to be optimized from '//trim(cparfile)
    open(ionum,file=trim(cparfile),status='old')
!.....num of symmetry functions, num of node in 1st hidden layer
10  read(ionum,'(a)') ctmp
    if( ctmp(1:1).eq.'!' .or. ctmp(1:1).eq.'#' ) then
      call parse_option_NN(ctmp)
      goto 10
    else
      read(ctmp,*) nn_nl
      backspace(ionum)
    endif
    allocate(nn_nhl(0:nn_nl+1))
    read(ionum,*) nn_nl, (nn_nhl(i),i=0,nn_nl)
    nn_nhl(nn_nl+1)= 1

!.....Bias treatment depends on NN version
    nvars= 0
    if( trim(cpotlow).eq.'nn2' ) then
      do i=1,nn_nl+1
        nvars= nvars +nn_nhl(i-1)*nn_nhl(i)
      enddo
    else if( trim(cpotlow).eq.'dnn' ) then
      do i=1,nn_nl+1
        nvars = nvars +(nn_nhl(i-1)+1)*nn_nhl(i)
      enddo
    endif
!!$!.....
!!$    read(ionum,*) nvars, rcut, rc3
  endif
  call mpi_bcast(nn_nl,1,mpi_integer,0,mpi_world,ierr)
  if( .not. allocated(nn_nhl) ) allocate(nn_nhl(0:nn_nl+1))
  call mpi_bcast(nn_nhl,nn_nl+2,mpi_integer,0,mpi_world,ierr)
  if( nn_nhl(0).ne.nsf_desc ) then
    if( myid.eq.0 ) then
      print *,'ERROR: num of inputs/descriptors are inconsistent '//&
           'between in.params.desc and '//trim(cparfile)
    endif
    stop 1
  endif
  call mpi_bcast(nn_sigtype,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(nn_asig,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(nvars,1,mpi_integer,0,mpi_world,ierr)
  allocate(vars(nvars),vranges(2,nvars))
  if( myid.eq.0 ) then
    print '(a,i0)',' Number of variables to be optimized = ',nvars
    do i=1,nvars
      read(ionum,*) vars(i),vranges(1:2,i)
    enddo
    if( index(cinitv,'gauss').ne.0 ) then
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
      write(6,'(a)') ' Potential parameters are read from file: '//trim(cparfile)
    endif
    close(ionum)
  endif
  call mpi_bcast(vars,nvars,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(vranges,2*nvars,mpi_real8,0,mpi_world,ierr)

end subroutine read_vars_NN
!=======================================================================
subroutine write_vars_NN(cfname)
!
!  Write parameters in NN2/DNN format
!
  use variables
  use parallel
  implicit none
  character(len=*),intent(in):: cfname
  integer:: i
  integer,parameter:: ionum = 16

  if( myid.eq.0 ) then
    open(ionum,file=trim(cfname),status='replace')
    write(ionum,'(a,i2)') '!  sigtype: ',nn_sigtype
    write(ionum,'(a)') '!  '
    write(ionum,'(100i4)') nn_nl, (nn_nhl(i),i=0,nn_nl)
!!$    write(ionum,'(i10,2es15.4)') nvars,rcut,rc3
    do i=1,nvars
      write(ionum,'(es23.14e3,2es12.4)') vars(i),vranges(1:2,i)
    enddo
    close(ionum)
!    print *, 'wrote '//trim(cfname)
  endif

end subroutine write_vars_NN
!=======================================================================
subroutine parse_option_NN(cline)
!
!  Parse options from a comment line.
!  Lines starting from ! or # are treated as comment lines,
!  but options can be given at the comment lines.
!  The option words should be put after these comment characters with
!  one or more spaces between them for example,
!
!  bias:  .true.
!
!  Currently available options are:
!    - "sigtype:" sigmoid type: 1, 2 or 3.
!    - "asig:" coefficient in activation function: [default: 0.01]
!
  use variables
  use parallel
  implicit none
  character(len=*),intent(in):: cline

  character(len=10):: c1,copt
  integer:: iopt
  real(8):: ropt

  ierr = 0
  if( index(cline,'sigtype:').ne.0 ) then
    read(cline,*) c1,copt,iopt
    nn_sigtype = iopt
  else if( index(cline,'asig:').ne.0 ) then
    read(cline,*) c1,copt,ropt
    nn_asig = ropt
  endif

end subroutine parse_option_NN
!=======================================================================
subroutine read_params_desc()
!
!  Read in.params.desc.
!
  use variables
  use parallel
  use util, only: num_data
  implicit none

  integer,parameter:: ionum = 16

  integer:: i,j,k,isp,jsp,ksp,isf,ityp,is1,is2
  logical:: lexist
  real(8):: rc,rcut2,rcut3,wgt
  character(len=128):: ctmp,cfname,cline,cmode
  character(len=3):: ccmb(3),csp

!.....initialize some
  ncnst_type(1) = 2   ! Gaussian
  ncnst_type(2) = 1   ! cosine
  ncnst_type(3) = 1   ! polynomial
  ncnst_type(4) = 2   ! Morse
  ncnst_type(101) = 1 ! angular1 (SW-type, not including fc(rjk))
  ncnst_type(102) = 2 ! angular2 (Behler-type, including fc(rjk))
  ncnst_type(103) = 1 ! cos(cos(thijk))
  ncnst_type(104) = 1 ! sin(cos(thijk))
  ncnst_type(105) = 2 ! exp(-eta*(cos(thijk)-rs)**2)
  ncnst_type(106) = 1 ! angular6 (no exp term in angular1)
  ncomb_type(1:100) = 2    ! pair
  ncomb_type(101:200) = 3  ! triplet

  if( myid.eq.0 ) then
    if( iprint.gt.1 ) print *,'Read in.params.desc...'
!.....read constants at the 1st call
    cfname = 'in.params.desc'
    inquire(file=trim(cfname),exist=lexist)
    if( .not. lexist ) then
      if( myid.eq.0 ) then
        write(6,'(a)') ' ERROR: '//trim(cfname)//' does not exist !!!.'
      endif
      stop
    endif
    open(ionum,file=trim(cfname),status='old')
!.....num of symmetry functions, num of node in 1st hidden layer
10  read(ionum,'(a)') ctmp
    if( ctmp(1:1).eq.'!' .or. ctmp(1:1).eq.'#' ) then
      call parse_option_desc(ctmp,iprint,ierr)
      goto 10
    else
      backspace(ionum)
    endif
!.....Read numbers of species and symmetry functions
    read(ionum,*) nsp_desc, nsf_desc
  endif

!.....Bcast nsp and nsf before allocating arrays
  call mpi_bcast(nsp_desc,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(nsf_desc,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(lcheby,1,mpi_logical,0,mpi_world,ierr)

!.....Allocate arrays of lenths, nsp and/or nsf
  if( .not.allocated(descs) ) then
    allocate(descs(nsf_desc),ilsf2(0:nsf_desc,nspmax,nspmax) &
         ,ilsf3(0:nsf_desc,nspmax,nspmax,nspmax))
  endif
  if( lcheby ) then
    allocate(wgtsp_desc(nspmax))
    do i=1,nsp_desc
      wgtsp_desc(i)= dble(i)
    enddo
  endif
  call mpi_bcast(nsf_desc,nspmax,mpi_real8,0,mpi_world,ierr)

  if( myid.eq.0 ) then
    if( lcheby ) then
!-----------------------------------------------------------------------
!  Input file format for Chebyshev (in.params.desc)
!-----------------------------------------------------------------------
!  ! Chebyshev:   T
!     3   100         ! nsp, nsf
!  2-body   30   5.00   ! 2-body, num of series (nsf2), rcut
!  3-body   20   4.00   ! 3-body, num of series (nsf3), rcut
!  Weight  Artrith    ! type of species-weight
!     La   1.0         ! In case of Artrith, specify (species,weight) pair
!     F   -1.0
!     Ca   2.0
!     ...
!-----------------------------------------------------------------------
!  Thus, users can specify different num of series and rcut for 2- and 3-body.
!  Note that the NSF should be identical to the total number of series.
!  If NSP==1, NSF=NSF2+NSF3, but if NSP>1, NSF=(NSF2+NSF3)*2.
!  Using a factor, NSFF, NSF=(NSF2+NSF3)*NSFF,
!  where NSFF=1 for NSP==1, and NSFF=2 for NSP>1.
!-----------------------------------------------------------------------
      if( iprint.gt.2 ) print *,'reading Chebyshev descriptors...'
      nsff_desc = 1
      nsf2_desc = 0
      nsf3_desc = 0
      if( nsp_desc.gt.1 ) nsff_desc = 2
      cmode = 'none'
      do while(.true.)
        read(ionum,'(a)',end=30) cline
        if( cline(1:1).eq.'!' .or. cline(1:1).eq.'#' ) cycle
        if( index(cline,'Weight').ne.0 .or. &
             index(cline,'weight').ne.0 ) then ! read Weight control
          cmode = 'Weight'  ! currently only Artrith type is available
          if( iprint.gt.0 ) write(6,'(a)') '   species-weight type: '//' Artrith'
        else if( index(cline,'2-body').ne.0 ) then
          cmode = 'none'
          read(cline,*) ctmp, nsf2_desc, rcut2
        else if( index(cline,'3-body').ne.0 ) then
          cmode = 'none'
          read(cline,*) ctmp, nsf3_desc, rcut3
        else
          if( trim(cmode).eq.'Weight' ) then
            read(cline,*,end=30) csp, wgt
            isp = csp2isp(trim(csp))
            if( isp.gt.0 ) then
              wgtsp_desc(isp) = wgt
              if( iprint.gt.0 ) write(6,'(5x,i2,a4,f6.1)') isp, trim(csp), wgt
            endif
          endif
        endif
      enddo
30    continue
!.....Check nsf vs nsf2,nsf3
      if( nsf_desc.ne.nsff_desc*(nsf2_desc+nsf3_desc) ) then
        print *,'ERROR @read_params_desc: nsf != (nsf2+nsf3)*nsff with nsff,nsp=',nsff_desc,nsp_desc
        stop 1
      endif
!.....Set rcut for all isf
      do isf=1,nsf2_desc*nsff_desc
        descs(isf)%rcut = rcut2
        descs(isf)%rcut2 = rcut2**2
      enddo
      do isf=nsf2_desc*nsff_desc+1,nsf_desc
        descs(isf)%rcut = rcut3
        descs(isf)%rcut2 = rcut3**2
      enddo

    else  ! not Chebyshev
      nsf2_desc = 0
      nsf3_desc = 0
      ilsf2(:,:,:) = 0
      ilsf3(:,:,:,:) = 0
      do isf=1,nsf_desc
        read(ionum,*,end=20) ityp,(ccmb(k),k=1,ncomb_type(ityp)) &
             ,rc,(cnst(j),j=1,ncnst_type(ityp))
        descs(isf)%itype = ityp
        isp = csp2isp(trim(ccmb(1)))
        jsp = csp2isp(trim(ccmb(2)))
        descs(isf)%isp = isp
        descs(isf)%jsp = jsp
        descs(isf)%rcut = rc
        descs(isf)%rcut2 = rc*rc
        descs(isf)%nprm = ncnst_type(ityp)
        if( .not.allocated(descs(isf)%prms) ) &
             allocate(descs(isf)%prms(descs(isf)%nprm))
        do j=1,descs(isf)%nprm
          descs(isf)%prms(j) = cnst(j)
        enddo
        if( isp.lt.0 .or. jsp.lt.0 ) cycle
        if( ityp.le.100 ) then  ! 2-body
          nsf2_desc = nsf2_desc + 1
          is1 = min(isp,jsp)
          is2 = max(isp,jsp)
          ilsf2(0,is1,is2) = ilsf2(0,is1,is2) + 1
          ilsf2(ilsf2(0,is1,is2),is1,is2) = isf
        else if( ityp.le.200 ) then  ! 3-body
          nsf3_desc = nsf3_desc + 1
          ksp = csp2isp(trim(ccmb(3)))
          if( ksp.lt.0 ) cycle
          descs(isf)%ksp = ksp
          is1 = min(jsp,ksp)
          is2 = max(jsp,ksp)
          ilsf3(0,isp,is1,is2) = &
               ilsf3(0,isp,is1,is2) +1
          ilsf3(ilsf3(0,isp,is1,is2),isp,is1,is2) = isf
        endif
      enddo  ! isf=1,nsf
20    continue
    endif ! lcheby
    close(ionum)
  endif ! myid.eq.0

!.....Broadcast some
  call bcast_descs()
  call mpi_bcast(nsf2_desc,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(nsf3_desc,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(nsff_desc,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(ilsf2,size(ilsf2),mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(ilsf3,size(ilsf3),mpi_integer,0,mpi_world,ierr)
  if( lcheby ) then
    call mpi_bcast(wgtsp_desc,nspmax,mpi_real8,0,mpi_world,ierr)
  endif

  return
end subroutine read_params_desc
!=======================================================================
subroutine parse_option_desc(cline,iprint,ierr)
!
!  Parse options from a comment line.
!  Lines starting from ! or # are treated as comment lines,
!  and options can be given at the comment lines.
!  The option words should be put after these comment characters with
!  one or more spaces between them for example,
!
!    Chebyshev:  T
!
!  Currently available options are:
!    - "Chebyshev:", toggle switch for Chebyshev polynomial series
!      ex) Chebyshev:  T 
!
  use variables,only: lcheby
  implicit none
  character(len=*),intent(in):: cline
  integer,intent(in):: iprint
  integer,intent(out):: ierr

  character(len=10):: c1,copt

  ierr = 0
  if( index(cline,'Chebyshev:').ne.0 .or. &
       index(cline,'chebyshev:').ne.0 ) then
    read(cline,*) c1, copt, lcheby
    if( iprint.gt.0 ) then
      print '(a)', ' Chebyshev series for descriptors.'
    endif
  endif

end subroutine parse_option_desc
!=======================================================================
subroutine bcast_descs()
  use variables
  use parallel

  integer:: i

  do i=1,nsf_desc
    call mpi_bcast(descs(i)%itype,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(descs(i)%isp,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(descs(i)%jsp,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(descs(i)%ksp,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(descs(i)%rcut,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(descs(i)%rcut2,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(descs(i)%nprm,1,mpi_integer,0,mpi_world,ierr)
    if( .not. lcheby ) then
      if( myid.ne.0 ) then
        if( allocated(descs(i)%prms) ) deallocate(descs(i)%prms)
        allocate(descs(i)%prms(descs(i)%nprm))
      endif
      call mpi_barrier(mpi_world,ierr)
      call mpi_bcast(descs(i)%prms,descs(i)%nprm,mpi_real8,0,mpi_world,ierr)
    endif
  enddo
end subroutine bcast_descs
!=======================================================================
subroutine read_params_ZBL()
!
!  Read in.params.ZBL for fitpot.
!
  use variables
  use parallel
  use force, only: loverlay
  use util, only: num_data
  implicit none

  integer,parameter:: ionum = 17

  integer:: isp,jsp
  real(8):: qnucli,ri,ro
  character(len=128):: cline,cfname,cmode,ctmp
  character(len=5):: cspi,cspj
  real(8),parameter:: qtiny = 1d-10
!!$    integer,external:: num_data

  if( myid.eq.0 ) then
    if( iprint.gt.1 ) print *,'Read in.params.ZBL...'
    
    cmode = ''
    cfname = 'in.params.ZBL'
    open(ionum,file=trim(cfname),status='old')
    zbl_interact(1:nspmax,1:nspmax) = .true.
    zbl_qnucl(1:nspmax) = 0d0
    zbl_rc = 0d0

    if( iprint.ne.0 ) write(6,'(/,a)') ' ZBL parameters:'
    do while(.true.)
      read(ionum,*,end=10) cline
      if( num_data(cline,' ').eq.0 ) cycle
      if( cline(1:1).eq.'#' .or. cline(1:1).eq.'!' ) cycle
!.....Mode detection
      if( trim(cline).eq.'parameters' ) then
        cmode = trim(cline)
        cycle
      else if( trim(cline).eq.'interactions' ) then
        cmode = trim(cline)
        zbl_interact(1:nspmax,1:nspmax) = .false.
        cycle
      endif
!.....Read parameters depending on the mode
      if( trim(cmode).eq.'parameters' ) then
        backspace(ionum)
        read(ionum,*) cspi, qnucli, ri, ro
        isp = csp2isp(trim(cspi))
        if( isp.le.0 ) cycle
        zbl_qnucl(isp) = qnucli
        zbl_ri(isp) = ri
        zbl_ro(isp) = ro
        zbl_rc = max(zbl_rc,ro)
        if( iprint.ne.0 ) then
          write(6,'(a,a3,3(2x,f7.3))') &
               '   csp,qnucl,ri,ro = ',trim(cspi),qnucli,ri,ro
        endif
      else if( trim(cmode).eq.'interactions' ) then
        backspace(ionum)
        read(ionum,*) cspi, cspj
        isp = csp2isp(trim(cspi))
        jsp = csp2isp(trim(cspj))
        if( isp.gt.0 .and. jsp.gt.0 ) then
          zbl_interact(isp,jsp) = .true.
          zbl_interact(jsp,isp) = .true.
        else
          print *,'  interacion read but not used: ',isp,jsp
        endif
      endif
    enddo
10  close(ionum)
    if( iprint.gt.1 ) then
      do isp=1,nspmax
        do jsp=isp,nspmax
          if( zbl_interact(isp,jsp) ) then
            write(6,'(a,2i3,l2)') '   isp,jsp,interact = ',isp,jsp,zbl_interact(isp,jsp)
          endif
        enddo
      enddo
    endif
  endif  ! myid.eq.0

  call mpi_bcast(zbl_rc,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(zbl_qnucl,nspmax,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(zbl_ri,nspmax,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(zbl_ro,nspmax,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(zbl_interact,nspmax*nspmax,mpi_logical,0,mpi_world,ierr)
  return

end subroutine read_params_ZBL
!-----------------------------------------------------------------------
! Local Variables:
! compile-command: "make fitpot"
! End:
