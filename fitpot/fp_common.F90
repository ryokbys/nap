module fp_common
!-----------------------------------------------------------------------
!                     Last modified: <2025-04-14 12:21:15 KOBAYASHI Ryo>
!-----------------------------------------------------------------------
!
! Module that contains common functions/subroutines for fitpot.
!
  use variables,only: cpenalty,penalty,pwgt2b,pwgt2bd,pwgt2bs, &
       pwgt3b,pwgt3bd,repul_radii
  use pmdvars,only: nspmax
  implicit none
  save
  real(8),allocatable:: fdiff(:,:),frcs(:,:),fref(:,:),fsub(:,:),gtrnl(:)
  logical,allocatable:: ldcover(:)  ! whether the parameter is data covered
  real(8),allocatable:: gwe(:),gwf(:,:,:),gws(:,:)
  real(8),allocatable:: gwesub(:)
  real(8):: pdiff(6), ptnsr(3,3)
  real(8):: epotsub

  real(8):: fac_etrn,fac_ftrn,fac_strn,fac_etst,fac_ftst,fac_stst
  logical:: fp_common_initialized= .false.

  integer,parameter:: ivoigt(3,3)= &
       reshape((/ 1, 6, 5, 6, 2, 4, 5, 4, 3 /),shape(ivoigt))

!.....Store loverlay and r_inner/outer arrays
  logical:: overlay

contains
!=======================================================================
  subroutine init_fp_common()
!!$    use variables,only: swgt2trn,swgt2tst,lematch,lfmatch,lsmatch
    use variables,only: lematch,lfmatch,lsmatch,wgte,wgtf,wgts, &
         evtrn,fvtrn,svtrn,netrn,nftrn,nstrn, &
         evtst,fvtst,svtst,netst,nftst,nstst, &
         nspmax,rcut,iprint,nff,cffs
    use parallel,only: myid
    use pmdvars, only: naux,nstp,nx,ny,nz,dt,rbuf,rc1nn,lvc
    use pmdvars,only: iprint_pmd => iprint, rc_pmd => rc
    use element,only: init_element
    use force,only: num_forces, force_list, loverlay

    integer:: nterms,i
    logical:: update_force_list

    fac_etrn = wgte
    fac_ftrn = wgtf
    fac_strn = wgts
    if( netrn > 1 ) fac_etrn = wgte /(evtrn*netrn)
    if( nftrn > 1 ) fac_ftrn = wgtf /(fvtrn*nftrn)
    if( nstrn > 1 ) fac_strn = wgts /(svtrn*nstrn)
    fac_etst = wgte
    fac_ftst = wgtf
    fac_stst = wgts
    if( netst > 1 ) fac_etst = wgte /(evtst*netst)
    if( nftst > 1 ) fac_ftst = wgtf /(fvtst*nftst)
    if( nstst > 1 ) fac_stst = wgts /(svtst*nstst)
    if( myid.eq.0 ) then
      write(6,'(/a)') ' Prefactors for loss function by terms (train,test):'
      write(6,'(a,2es14.3)') '   Energy: ', fac_etrn, fac_etst
      write(6,'(a,2es14.3)') '   Force:  ', fac_ftrn, fac_ftst
      write(6,'(a,2es14.3)') '   Stress: ', fac_strn, fac_stst
    endif

!.....Create MPI COMM for pmd only for the 1st time
    call create_mpi_comm_pmd()
    call init_element()

    nstp = 0
    dt = 1d0
    rbuf = 0.0d0
    rc1nn = 3.0d0
    rc_pmd = rcut

    nx = 1
    ny = 1
    nz = 1
    iprint_pmd = max(0,iprint-10)

    lvc = .false.

!.....Set force_list in the force module
    update_force_list = .false.
    if( nff.ne.num_forces ) then
      update_force_list = .true.
    else
      do i=1,nff
        if( trim(cffs(i)).ne.trim(force_list(i)) ) then
          update_force_list = .true.
          exit
        endif
      enddo
    endif
!.....Update force_list if needed
    if( update_force_list ) then
      num_forces = nff
      do i=1,num_forces
        force_list(i) = trim(cffs(i))
      enddo
    endif
!.....Overlay setting
    if( overlay ) then
!.....Set loverlay variable in force module
      loverlay = overlay
    endif

!.....init_force would perform read_params_XXX for force_XXX.F90 in it.
    call init_force(.true.)
    call set_cauxarr()  ! pmd_core.F90

    fp_common_initialized = .true.

  end subroutine init_fp_common
!=======================================================================
  subroutine calc_swgts(swgts_trn,swgts_tst)
!
!  Calculate prefactor as a sum of sample weights.
!
    use variables,only: samples
    use parallel
    real(8),intent(out):: swgts_trn,swgts_tst

    integer:: ismpl
    real(8):: swgtrn,swgtst
    
!.....set nominator for sample weights
    swgtrn = 0d0
    swgtst = 0d0
    do ismpl=isid0,isid1
      if( samples(ismpl)%iclass.eq.1 ) then
        swgtrn = swgtrn +samples(ismpl)%wgt
      else if(samples(ismpl)%iclass.eq.2 ) then
        swgtst = swgtst +samples(ismpl)%wgt
      endif
    enddo
    swgts_trn = 0d0
    swgts_tst = 0d0
    call mpi_allreduce(swgtrn,swgts_trn,1,mpi_real8,mpi_sum &
         ,mpi_world,ierr)
    call mpi_allreduce(swgtst,swgts_tst,1,mpi_real8,mpi_sum &
         ,mpi_world,ierr)
    
    return
  end subroutine calc_swgts
!=======================================================================
  subroutine wrap_ranges(ndim,x,xranges)
    use variables,only: cpot, nsp, short_radii, l_correct_short
    use uf3,only: uf3_short_correction
    implicit none
    integer,intent(in):: ndim
    real(8),intent(in):: xranges(2,ndim)
    real(8),intent(inout):: x(ndim)

    integer:: i

!.....Short-range correction for uf3
    if( trim(cpot).eq.'uf3' .and. l_correct_short ) then
      if( allocated(ldcover) ) then
        call uf3_short_correction(ndim,x,nsp,short_radii,ldcover)
      endif
    endif

    do i=1,ndim
      if( x(i).lt.xranges(1,i) ) then
        x(i) = xranges(1,i)
      else if( x(i).gt.xranges(2,i) ) then
        x(i) = xranges(2,i)
      endif
    enddo
    return
  end subroutine wrap_ranges
!=======================================================================
  subroutine mask_grad(ndim,vranges,grad)
!
!  Set grad 0.0 if vranges is super narrow so that the var to be fixed.
!
    integer,intent(in):: ndim
    real(8),intent(in):: vranges(2,ndim)
    real(8),intent(inout):: grad(ndim)

    integer:: i
    logical,save,allocatable:: lfix(:)

    if( .not.allocated(lfix) ) then
      allocate(lfix(ndim))
      lfix(:) = .false.
      do i=1,ndim
        if( abs(vranges(1,i)-vranges(2,i)) < 1d-14 ) lfix(i) = .true.
      enddo
    endif

    do i=1,ndim
      if( lfix(i) ) grad(i) = 0d0
    enddo
    return
  end subroutine mask_grad
!=======================================================================
  subroutine func_w_pmd(ndim,x,ftrn,ftst)
!
!  Evaluate loss function value using pmd (actually one_shot routine.)
!
    use variables,only:samples,tfunc, &
         lematch,lfmatch,lsmatch,nfunc,tcomm,twait,mdsys, &
         swgt2trn,swgt2tst,cpot,ismask, &
         nff,cffs,maxna,rcut,force_limit,stress_limit, &
         crefstrct,erefsub,myidrefsub,isidrefsub,iprint, &
         ctype_loss,dmem,cfmethod,cfrc_denom,cstrs_denom, &
         lnormalize,lnormalized,lgdw,lgdwed,terg,tfrc,tstrs, &
         nn_nl, nn_nhl, nn_sigtype, nn_asig, &
         wgte,wgtf,wgts,netrn,nftrn,nstrn,evtrn,fvtrn,svtrn, &
         repul_radii,nsp,specorder
    use parallel
    use descriptor,only: lupdate_gsf,get_descs,get_ints
    use DNN,only: nlayer, nhl, itypesig, asig
    use UF3,only: calc_short_lossfunc
    implicit none
    integer,intent(in):: ndim
    real(8),intent(in):: x(ndim)
    real(8),intent(out):: ftrn,ftst

    integer:: ismpl,natm,ia,ixyz,jxyz,k,nsf,nal,nnl,jfcal,ir
    integer:: isp,jsp
    character(len=3):: csi,csj
    real(8):: ediff,eref,epot,swgt,esub,gsfmem
    real(8):: eerr,ferr,ferri,serr,serri,strs(3,3),absfref,abssref, &
         sref(3,3),ssub(3,3)
    real(8):: ftrnl,ftstl,ftmp,gdw,fpena,fetmp,fftmp,fstmp
    real(8):: tfl,tcl,tfg,tcg,tf0,tc0,tw0,twl,twg,tsmp0
    real(8):: tergl,tfrcl,tstrsl,tmp
!!$    real(8):: fac_e, fac_f, fac_s
    type(mdsys):: smpl
    logical,save:: l1st = .true.
    logical,parameter:: lgrad = .false.
    logical,parameter:: lgrad_done = .false.
    logical:: lfdsgnmat
    character(len=128):: csmplname
    real(8),allocatable,save:: rs(:)
    logical,external:: string_in_arr

    call flush(6)
    nfunc= nfunc +1

    tc0= mpi_wtime()
    call mpi_bcast(x,ndim,mpi_real8,0,mpi_world,ierr)
    tcl= mpi_wtime() -tc0
    tf0= mpi_wtime()

    lfdsgnmat = .false.  ! Initialize lfdsgnmat

    if( l1st ) then
      if( .not.allocated(fdiff) ) then
        allocate(fdiff(3,maxna),frcs(3,maxna),fref(3,maxna),fsub(3,maxna))
        dmem = dmem +8d0*size(fdiff) +8d0*size(frcs)*3
      endif
      if( index(cpot,'NN').ne.0 .or. trim(cpot).eq.'linreg') lupdate_gsf = .true.
      gsfmem = 0d0
    endif

    if( .not. lematch .and. .not.lfmatch .and. .not.lsmatch ) then
      if( myid.eq.0 ) then
        print *,'Nothing to be fitted.'
      endif
      call mpi_finalize(ierr)
      stop
    endif

    if( trim(cpot).eq.'uf3' .and. l1st .and. abs(pwgt2bs)>1d-14 ) then
      repul_radii(:,:) = 1d+10
    endif
    
    do ismpl=isid0,isid1
      tsmp0 = mpi_wtime()
      if( allocated(ismask) ) then
        if( ismask(ismpl).ne.0 ) cycle
      endif
!.....CAUTION: assignment of a structure-type is not pointer copy,
!              rather a copy of its contents.
!              So even if one changes the variable in the copied structure,
!              the variable in the original structure (e.g., samples(ismpl))
!              will not be changed.
      natm= samples(ismpl)%natm
      csmplname = trim(samples(ismpl)%csmplname)
!.....CAUTION 2: argument of a structure-type is a pointer, contrary to the above.
      call pre_pmd(samples(ismpl),ndim,x,nff,cffs,rcut,l1st)

!.....Set lfdsgnmat=.true. to make run_pmd() compute dsgnmat_force related data
      if( trim(cpot).eq.'linreg' .and. &
           l1st .and. lfmatch .and. trim(cfmethod).eq.'dsgnmat' ) then
          lfdsgnmat = .true.
      endif
!!$      print '(a,2i5,f8.4)','func: myid,ismpl,tsmpl 01=',myid,ismpl,mpi_wtime()-tsmp0

      call run_pmd(samples(ismpl),lgrad,lgrad_done,ndim,epot,frcs,strs &
           ,rcut,lfdsgnmat)
!!$      print '(a,2i5,f8.4)','func: myid,ismpl,tsmpl 02=',myid,ismpl,mpi_wtime()-tsmp0
      samples(ismpl)%epot = epot
      samples(ismpl)%fa(1:3,1:natm) = frcs(1:3,1:natm)
      samples(ismpl)%strs(1:3,1:3) = strs(1:3,1:3)
      if( trim(cpot).eq.'linreg' .or. index(cpot,'nn').ne.0 ) then
        if( .not. allocated(samples(ismpl)%gsf) ) then
          call get_ints(nsf,nal,nnl)
          samples(ismpl)%nsf = nsf
          samples(ismpl)%nal = nal
          samples(ismpl)%nnl = nnl
          allocate(samples(ismpl)%gsf(nsf,nal) )
          gsfmem = gsfmem +8d0*size(samples(ismpl)%gsf)
          dmem = dmem +gsfmem
        endif
        call get_descs(samples(ismpl)%nsf,samples(ismpl)%nal, &
             samples(ismpl)%nnl,samples(ismpl)%gsf)
      endif
      if( trim(cpot).eq.'uf3' .and. l1st .and. abs(pwgt2bs)>1d-14 ) then
        call get_shortest_distances(repul_radii)  ! in pmd_core
      endif
!!$      print '(a,2i5,f8.4)','func: myid,ismpl,tsmpl=',myid,ismpl,mpi_wtime()-tsmp0
    enddo  ! ismpl

    if( l1st .and. index(cpot,'nn').ne.0 ) then
      nn_nl = nlayer -1
      if( .not.allocated(nn_nhl) ) allocate(nn_nhl(0:nn_nl))
      nn_nhl(0:nn_nl) = nhl(0:nlayer-1)
      nn_asig = asig
      nn_sigtype = itypesig
    endif

    if( lnormalize .and. .not.lnormalized ) call normalize()
    if( lgdw .and. .not.lgdwed ) call compute_gdw()

    if( len(trim(crefstrct)).gt.5 ) then
      if( myid.eq.myidrefsub ) then
        epotsub = samples(isidrefsub)%epot +samples(isidrefsub)%esub
        epotsub = epotsub /samples(isidrefsub)%natm
      endif
      call mpi_bcast(epotsub,1,mpi_real8,myidrefsub,mpi_world,ierr)
    endif

    ftrnl = 0d0
    ftstl = 0d0
    tergl = 0d0
    tfrcl = 0d0
    tstrsl= 0d0
    do ismpl=isid0,isid1
      if( allocated(ismask) ) then
        if( ismask(ismpl).ne.0 ) cycle
      endif
      smpl = samples(ismpl)
      csmplname= smpl%csmplname
      natm = smpl%natm
      epot = smpl%epot
      ftmp = 0d0
      fetmp = 0d0
      fftmp = 0d0
      fstmp = 0d0
      swgt = smpl%wgt
!.....Energy matching
      if( lematch ) then
        tmp = mpi_wtime()
        eref= smpl%eref
        esub= smpl%esub
!!$        eerr = smpl%eerr
        eerr = 1d0
        if( len(trim(crefstrct)).gt.5 ) then
          ediff= (epot-epotsub*natm+esub &
               -(eref-erefsub*natm))/natm /eerr
        else
          ediff= (epot+esub -eref)/natm /eerr
        endif
        if( trim(ctype_loss).eq.'huber' ) then
          if( abs(ediff).gt.1.d0 ) then
            ediff = 2d0*abs(ediff) -1d0
          else
            ediff = ediff *ediff
          endif
        else ! LS as default
          ediff= ediff*ediff
        endif
        fetmp= fetmp +ediff *swgt
        if( iprint.gt.2 ) then
          write(6,'(a,2i4,1x,a,7es11.3)') ' myid,ismpl,smplname,epot,eref,esub,(epot+esub)/natm= ', &
               myid,ismpl,trim(csmplname),epot,eref,esub,(epot+esub)/natm
        endif
        tergl = tergl +mpi_wtime() -tmp
      endif
!.....Force matching
      if( lfmatch .and. smpl%nfcal.ne.0 ) then
        tmp = mpi_wtime()
        frcs(1:3,1:natm) = smpl%fa(1:3,1:natm)
        fref(1:3,1:natm) = smpl%fref(1:3,1:natm)
        fsub(1:3,1:natm) = smpl%fsub(1:3,1:natm)
!!$        ferr = smpl%ferr
        ferr = 1d0
!!$        ferri = 1d0/ferr
        jfcal = 0
        do ia=1,natm
          if( .not. smpl%lfrc_eval(ia) ) cycle
          jfcal = jfcal +1
          gdw = 1d0
          if( lgdw ) gdw = smpl%gdw(ia)
          absfref = sqrt(fref(1,ia)**2 +fref(2,ia)**2 +fref(3,ia)**2)
          if( cfrc_denom(1:7).eq.'abs2rel' ) then
            if( absfref.lt.force_limit ) then
              ferri = 1d0 /ferr
            else
              ferri = 1d0 /max(absfref,ferr)
            endif
          else if( cfrc_denom(1:3).eq.'rel' ) then
            ferri = 1d0/ max(absfref,ferr)
          else  ! default: abs
            ferri = 1d0/ ferr
          endif
          do ixyz=1,3
            fdiff(ixyz,ia)= (frcs(ixyz,ia)+fsub(ixyz,ia) &
                 -fref(ixyz,ia)) *ferri
            if( trim(ctype_loss).eq.'huber' ) then
              if( abs(fdiff(ixyz,ia)).gt.1d0 ) then
                fdiff(ixyz,ia) = 2d0*abs(fdiff(ixyz,ia)) -1d0
              else
                fdiff(ixyz,ia)= fdiff(ixyz,ia)*fdiff(ixyz,ia)
              endif
            else ! LS as default
              fdiff(ixyz,ia)= fdiff(ixyz,ia)*fdiff(ixyz,ia)
            endif
            fftmp= fftmp +fdiff(ixyz,ia) *swgt *gdw
          enddo
        enddo
        tfrcl = tfrcl +mpi_wtime() -tmp
      endif

!.....Stress matching
      if( lsmatch ) then
        tmp = mpi_wtime()
!.....Compare these ptnsr elements with sref elements
!!$        serr = smpl%serr
        serr = 1d0
        strs(:,:) = smpl%strs(:,:)
        sref(:,:) = smpl%sref(:,:)
        ssub(:,:) = smpl%ssub(:,:)
        pdiff(1:6) = 0d0
        abssref = abs(sref(1,1) + sref(2,2) + sref(3,3))/3
        if( cstrs_denom(1:7).eq.'abs2rel' ) then
          if( abssref.lt.stress_limit ) then
            serri = 1d0/ serr
          else
            serri = 1d0/ max(abssref,serr)
          endif
        else if( cstrs_denom(1:3).eq.'abs' ) then
          serri = 1d0/ max(abssref,serr)
        else  ! default: relative
          serri = 1d0/ serr
        endif
!.....Skip abnormally large stress
!!$        if( stress_limit.lt.0d0 .or. &
!!$             .not.(cstrs_denom(1:7).ne.'abs2rel' .and. abssref.gt.stress_limit) ) then
          do ixyz=1,3
            do jxyz=ixyz,3
              k = ivoigt(ixyz,jxyz)
              pdiff(k)= pdiff(k) +(strs(ixyz,jxyz)  &
                   +ssub(ixyz,jxyz) &
                   -sref(ixyz,jxyz)) *serri
            enddo
          enddo
          if( trim(ctype_loss).eq.'huber' ) then
            do k=1,6
              if( abs(pdiff(k)).gt.1d0 ) then
                pdiff(k) = 2d0*abs(pdiff(k)) -1d0
              else
                pdiff(k)= pdiff(k)*pdiff(k)
              endif
              fstmp= fstmp +pdiff(k) *swgt
            enddo
          else  ! LS as default
            do k=1,6
              pdiff(k)= pdiff(k)*pdiff(k)
              fstmp= fstmp +pdiff(k) *swgt
            enddo
          endif
!!$        endif
        tstrsl = tstrsl +mpi_wtime() -tmp
      endif  ! stress matching

      if( smpl%iclass.eq.1 ) then
        ftrnl = ftrnl +fetmp*fac_etrn +fftmp*fac_ftrn +fstmp*fac_strn
      else if( smpl%iclass.eq.2 ) then
        ftstl = ftstl +fetmp*fac_etst +fftmp*fac_ftst +fstmp*fac_stst
      endif
    enddo  ! ismpl

    terg = terg + tergl
    tfrc = tfrc + tfrcl
    tstrs = tstrs + tstrsl

    tfl = mpi_wtime() -tf0

    tw0 = mpi_wtime()
    call mpi_barrier(mpi_world,ierr)
    twl = mpi_wtime() -tw0

    ftrn= 0d0
    ftst = 0d0
    tc0= mpi_wtime()
    call mpi_allreduce(ftrnl,ftrn,1,mpi_real8,mpi_sum,mpi_world,ierr)
    call mpi_allreduce(ftstl,ftst,1,mpi_real8,mpi_sum,mpi_world,ierr)
    tcl = tcl + (mpi_wtime() -tc0)
!!$    ftrn = ftrn /swgt2trn
!!$    if( swgt2tst.gt.1d-5 ) then
!!$      ftst = ftst /swgt2tst
!!$    endif

!.....Compute repulsion gradient
    if( l1st .and. trim(cpot).eq.'uf3' .and. abs(pwgt2bs)>1d-14 ) then
!.....Merge minimum distances in all nodes
      call mpi_allreduce(mpi_in_place,repul_radii,nspmax*nspmax, &
           mpi_real8,mpi_min,mpi_world,ierr)
!.....Print out repul_radii
      if( myid==0 ) print '(a)',' Shortest distances of pairs in dataset:'
      do isp=1,nsp
!!$        call set_charges(isp, valence_chgs(isp), core_chgs(isp))
        csi = specorder(isp)
        do jsp=isp,nsp
          csj = specorder(jsp)
          if( myid==0 ) then
            print '(a,f10.4)','  '//trim(csi)//'-'//trim(csj)//': ', repul_radii(isp,jsp)
          endif
        enddo
      enddo
!!$      drepul_tbl(:,:,:) = 0d0
!!$      if( .not.allocated(rs) ) allocate(rs(n_repul_pnts))
!!$      do isp=1,nsp
!!$        do jsp=isp,nsp
!!$          do ir=1,n_repul_pnts
!!$!.....rs -- r-points to be evaluated as mid-points of sections
!!$            rs(ir) = repul_radii(isp,jsp)/n_repul_pnts *(dble(ir)-0.5d0)
!!$          enddo
!!$          call calc_repul_grad(isp,jsp,n_repul_pnts,rs, &
!!$               drepul_tbl(:,isp,jsp))
!!$          drepul_tbl(:,jsp,isp) = drepul_tbl(:,isp,jsp)
!!$        enddo
!!$      enddo
    endif

!.....Penalty to function
    call func_penalty(ndim,x,fpena)
    ftrn = ftrn +fpena
!.....Repulsion correction
!!$    if( trim(cpot).eq.'uf3' .and. n_repul_pnts > 0 ) then
!!$      call calc_short_lossfunc(n_repul_pnts,repul_radii,drepul_tbl,frepul)
!!$      ftrn = ftrn +pwgt_repul*frepul
!!$    endif

!.....only the bottle-neck times are taken into account
    tcg = tcl
    tfg = tfl
    twg = twl
    call mpi_reduce(tcl,tcg,1,mpi_real8,mpi_max,0,mpi_world,ierr)
    call mpi_reduce(tfl,tfg,1,mpi_real8,mpi_max,0,mpi_world,ierr)
    call mpi_reduce(twl,twg,1,mpi_real8,mpi_max,0,mpi_world,ierr)
    tcomm= tcomm +tcg
    tfunc= tfunc +tfg
    twait= twait +twg

    if( l1st .and. myid.eq.0 .and. iprint.gt.1 ) then
      print '(a,f0.3,a)',' Memory for gsfs = ',gsfmem/1000/1000,' MB'
    endif
    l1st = .false.
!!$    if( index(cpot,'NN').ne.0 .or. trim(cpot).eq.'linreg' ) lupdate_gsf = .false.

  end subroutine func_w_pmd
!=======================================================================
  subroutine grad_w_pmd(ndim,x,gtrn)
!
!  Evaluate the gradient of loss function value
!  using pmd (actually one_shot routine.)
!
    use variables,only: tgrad,ngrad,tcomm,tgrad,twait, &
         samples,mdsys,swgt2trn,nff,cffs,force_limit,stress_limit, &
         maxna,maxnf,lematch,lfmatch,lsmatch,erefsub,crefstrct, &
         rcut,myidrefsub,isidrefsub,iprint, &
         ctype_loss,cfrc_denom,cstrs_denom,lgdw,dmem,terg,tfrc,tstrs, &
         wgte,wgtf,wgts,netrn,nftrn,nstrn,evtrn,fvtrn,svtrn,cpot, &
         vranges,ismask, repul_radii
    use parallel
!!$    use minimize
    use UF3,only: get_mem_uf3, get_mem_uf3l, dealloc_gwx_uf3, &
         calc_short_lossgrad
    implicit none
    integer,intent(in):: ndim
    real(8),intent(in):: x(ndim)
    real(8),intent(out):: gtrn(ndim)

    integer:: ismpl,natm,k,ia,ixyz,jxyz,iv,iff,i,jfcal,nfcal
    real(8):: tcl,tgl,tcg,tgg,tc0,tg0,tw0,twl,twg,tsmp0, &
         epot,esub,strs(3,3), &
         sref(3,3),ssub(3,3)
    real(8):: ediff,eerr,eref,swgt,ferr,ferri,serr,serri,tmp,gdw
    real(8):: absfref,abssref
    real(8):: tergl, tfrcl, tstrsl, ttmp
    type(mdsys):: smpl
    logical,save:: l1st = .true.
    logical,parameter:: lgrad = .true.
    logical:: lgrad_done = .false.
    logical,parameter:: lfdsgnmat = .false.
    character(len=128):: csmplname
    real(8),allocatable,save:: gpena(:)
    
    logical,external:: string_in_arr

    call flush(6)
    if( .not.allocated(gtrnl) ) then
      allocate(gtrnl(ndim))
      dmem = dmem +8d0*size(gtrnl)
    endif
    if( .not.allocated(gwe) ) then
      allocate(gwe(ndim),gwf(3,ndim,maxnf),gws(6,ndim))
      dmem = dmem +8d0*size(gwe) +8d0*size(gwf) +8d0*size(gws)
!!$      if( myid.eq.0 ) print *,'grad_w_pmd: dmem,size(gwf)=',dmem,size(gwf)
    endif
    if( len(trim(crefstrct)).gt.5 ) then
      if( .not.allocated(gwesub) ) then
        allocate(gwesub(ndim))
        dmem = dmem +8d0*size(gwesub)
      endif
    endif

    ngrad= ngrad +1
    tc0= mpi_wtime()
    call mpi_bcast(x,ndim,mpi_real8,0,mpi_world,ierr)
    tcl= mpi_wtime() -tc0
    tg0= mpi_wtime()

    if( .not. lematch .and. .not.lfmatch .and. .not.lsmatch ) then
      if( myid.eq.0 ) then
        print *,'Nothing to be fitted.'
      endif
      call mpi_finalize(ierr)
      stop
    endif

    if( len(trim(crefstrct)).gt.5 ) then
      if( myid.eq.myidrefsub ) then
        epotsub = samples(isidrefsub)%epot +samples(isidrefsub)%esub
        epotsub = epotsub /samples(isidrefsub)%natm
!!$        gwesub(1:ndim) = samples(isidrefsub)%gwe(1:ndim)
      endif
      call mpi_bcast(epotsub,1,mpi_real8,myidrefsub,mpi_world,ierr)
!!$      call mpi_bcast(gwesub,ndim,mpi_real8,myidrefsub,mpi_world,ierr)
    endif

    gtrnl(1:ndim) = 0d0
    tergl = 0d0
    tfrcl = 0d0
    tstrsl = 0d0
    do ismpl=isid0,isid1
      tsmp0 = mpi_wtime()
      if( allocated(ismask) ) then
        if( ismask(ismpl).ne.0 ) cycle
      endif
      smpl = samples(ismpl)
      natm= smpl%natm
      nfcal = smpl%nfcal
      csmplname = smpl%csmplname
      if( iprint.gt.10 ) print *,'grad_w_pmd: myid,ismpl,iclass,csmplname=', &
           myid,ismpl,smpl%iclass,trim(csmplname)

!.....Since g calc is time consuming,
!.....not calculate g for test set.
      if( smpl%iclass.ne.1 ) cycle

!.....Since the derivative of linear regression does not change
!     even after the variables change, we dont need to recalculate gwx().
!     And thus only copy them calculated at the 1st step,
!     which will reduce a lot of computational cost but require a lot of memory.
      if( (trim(cpot).eq.'uf3' .or. trim(cpot).eq.'linreg') &
           .and. (allocated(samples(ismpl)%gwe) .or. allocated(samples(ismpl)%gwf) &
           .or. allocated(samples(ismpl)%gws) ) ) then
        if( lematch ) gwe(:) = smpl%gwe(:)
        if( lfmatch ) gwf(:,:,1:nfcal) = smpl%gwf(:,:,1:nfcal)
        if( lsmatch ) gws(:,:) = smpl%gws(:,:)
      else
        if( iprint.gt.10 ) print *,'grad_w_pmd: pre_pmd for csmplname: ',trim(csmplname)
        call pre_pmd(samples(ismpl),ndim,x,nff,cffs,rcut,.false.)
        if( iprint.gt.10 ) print *,'grad_w_pmd: run_pmd for csmplname: ',trim(csmplname)
        lgrad_done = smpl%lgrad_done
!.....Note: since lgrad==.true., epot, frcs, strs are not calculated in this run_pmd.
        call run_pmd(samples(ismpl),lgrad,lgrad_done,ndim,epot,frcs,strs,rcut &
             ,lfdsgnmat,gwe,gwf,gws)
        if( (trim(cpot).eq.'uf3' .or. trim(cpot).eq.'linreg') ) then
!!$          allocate(samples(ismpl)%gwe(ndim), samples(ismpl)%gwf(3,ndim,maxnf), &
!!$               samples(ismpl)%gws(6,ndim))
!!$          dmem = dmem +4d0*(size(gwe) +size(gwf) +size(gws))
          if( lematch ) then
            if( .not.allocated(samples(ismpl)%gwe) ) allocate(samples(ismpl)%gwe(ndim))
            samples(ismpl)%gwe(:)= gwe(:)
            dmem = dmem +8d0*size(gwe)
          endif
          if( lfmatch ) then
            if( .not.allocated(samples(ismpl)%gwf) ) allocate(samples(ismpl)%gwf(3,ndim,maxnf))
            samples(ismpl)%gwf(:,:,1:nfcal)= gwf(:,:,1:nfcal)
            dmem = dmem +8d0*size(gwf)
          endif
          if( lsmatch ) then
            if( .not.allocated(samples(ismpl)%gws) ) allocate(samples(ismpl)%gws(6,ndim))
            samples(ismpl)%gws(:,:)= gws(:,:)
            dmem = dmem +8d0*size(gws)
          endif
        endif
      endif
      
!!$    enddo  ! ismpl
!!$
!!$    do ismpl=isid0,isid1
!!$      if( allocated(ismask) ) then
!!$        if( ismask(ismpl).ne.0 ) cycle
!!$      endif
!!$      smpl= samples(ismpl)
!!$!.....Since g calc is time consuming,
!!$!.....not calculate g for test set.
!!$      if( smpl%iclass.ne.1 ) cycle
      
      natm= smpl%natm
      epot= smpl%epot
      swgt= smpl%wgt
!.....Derivative of energy term w.r.t. weights
      if( lematch ) then
        if( iprint > 10 ) print *,'(grad_w_pmd) in lematch'
        ttmp = mpi_wtime()
        eref= smpl%eref
        esub= smpl%esub
!!$        eerr= smpl%eerr
        eerr= smpl%eerr
        eerr = 1d0
        if( len(trim(crefstrct)).gt.5 ) then
          ediff= (epot-epotsub*natm+esub -(eref-erefsub*natm))/natm /eerr
          if( trim(ctype_loss).eq.'LS' ) then
            tmp = 2d0 *ediff
          else  ! Huber
            if( abs(ediff).gt.1d0 ) then
              tmp = 2d0 *sign(1d0,ediff)
            else
              tmp = 2d0 *ediff
            endif
          endif
          gtrnl(1:ndim) = gtrnl(1:ndim) &
               +tmp/natm/eerr *swgt *fac_etrn &
               *(gwe(1:ndim) -gwesub(1:ndim))
!!$               *(smpl%gwe(1:ndim) -gwesub(1:ndim))
        else
          ediff= (epot+esub -eref)/natm /eerr
          if( trim(ctype_loss).eq.'LS' ) then
            tmp = 2d0 *ediff
          else  ! Huber
            if( abs(ediff).gt.1d0 ) then
              tmp = 2d0 *sign(1d0,ediff)
            else
              tmp = 2d0 *ediff
            endif
          endif
          gtrnl(1:ndim) = gtrnl(1:ndim) &
               +tmp*gwe(1:ndim)/natm/eerr *swgt *fac_etrn
!!$               +tmp*smpl%gwe(1:ndim)/natm/eerr *swgt
        endif
        tergl = tergl +mpi_wtime() -ttmp
      endif
!.....Derivative of force term w.r.t. weights
      if( lfmatch ) then
        ttmp = mpi_wtime()
        frcs(1:3,1:natm)= smpl%fa(1:3,1:natm)
        fref(1:3,1:natm)= smpl%fref(1:3,1:natm)
        fsub(1:3,1:natm)= smpl%fsub(1:3,1:natm)
!!$        ferr= smpl%ferr
        ferr= 1d0
        jfcal = 0
        do ia=1,natm
          if( .not. smpl%lfrc_eval(ia) ) cycle
          jfcal = jfcal +1
          gdw = 1d0
          if( lgdw ) gdw = smpl%gdw(ia)
          absfref = sqrt(fref(1,ia)**2 + fref(2,ia)**2 + fref(3,ia)**2)
          if( cfrc_denom(1:7).eq.'abs2rel' ) then
            if( absfref.lt.force_limit ) then
              ferri = 1d0 /ferr
            else
              ferri = 1d0 /max(absfref,ferr)
            endif
          else if( cfrc_denom(1:3).eq.'rel' ) then
            ferri = 1d0/ max(absfref,ferr)
          else  ! default: abs
            ferri = 1d0/ ferr
          endif
          do ixyz=1,3
            fdiff(ixyz,ia)= (frcs(ixyz,ia) + fsub(ixyz,ia) &
                 -fref(ixyz,ia)) *ferri
            if( trim(ctype_loss).eq.'LS' ) then
              tmp = 2d0 *fdiff(ixyz,ia) *ferri
            else  ! Huber
              if( abs(fdiff(ixyz,ia)).gt.1d0 ) then
                tmp = 2d0 *sign(1d0,fdiff(ixyz,ia))
              else
                tmp = 2d0 *fdiff(ixyz,ia) *ferri
              endif
            endif
            gtrnl(1:ndim)= gtrnl(1:ndim) +tmp &
                 *gwf(ixyz,1:ndim,jfcal) *swgt *gdw *fac_ftrn
          enddo  ! ixyz=1,3
        enddo  ! ia=1,natm
        tfrcl = tfrcl +mpi_wtime() -ttmp
      endif
!.....Derivative of stress w.r.t. weights
      if( lsmatch ) then
        ttmp = mpi_wtime()
!!$        serr= smpl%serr
        serr= 1d0
        strs(:,:) = smpl%strs(:,:)
        sref(:,:) = smpl%sref(:,:)
        ssub(:,:) = smpl%ssub(:,:)
        pdiff(1:6) = 0d0
        abssref = abs(sref(1,1) + sref(2,2) + sref(3,3))/3
        if( cstrs_denom(1:7).eq.'abs2rel' ) then
          if( abssref.lt.stress_limit ) then
            serri = 1d0/ serr
          else
            serri = 1d0/ max(abssref,serr)
          endif
        else if( cstrs_denom(1:3).eq.'rel' ) then
          serri = 1d0/ max(abssref,serr)
        else ! default: absolute
          serri = 1d0/ serr
        endif
!.....Skip abnormally large stress?
!!$        if( stress_limit.lt.0d0 .or. &
!!$             .not.(cstrs_denom(1:7).ne.'abs2rel' .and. abssref.gt.stress_limit) ) then
          do ixyz=1,3
            do jxyz=ixyz,3
              k = ivoigt(ixyz,jxyz)
              pdiff(k)= pdiff(k) +(strs(ixyz,jxyz) + ssub(ixyz,jxyz) &
                   - sref(ixyz,jxyz)) *serri
            enddo
          enddo
          do k=1,6
            if( trim(ctype_loss).eq.'LS' ) then
              tmp = 2d0 *pdiff(k) *serri 
            else  ! Huber
              if( abs(pdiff(k)).gt.1d0 ) then
                tmp = 2d0 *sign(1d0,pdiff(k))
              else
                tmp = 2d0 *pdiff(k) *serri 
              endif
            endif
            gtrnl(1:ndim)= gtrnl(1:ndim) +tmp *gws(k,1:ndim) *swgt *fac_strn
          enddo
!!$        endif
        tstrsl = tstrs +mpi_wtime() -ttmp
      endif
!!$      print '(a,3i5,f8.4)','grad: myid,ismpl,nfcal,tsmpl=',myid,ismpl,smpl%nfcal,mpi_wtime()-tsmp0
    enddo  ! ismpl


    terg = terg +tergl
    tfrc = tfrc +tfrcl
    tstrs = tstrs +tstrsl
    tgl= mpi_wtime() -tg0

    tw0 = mpi_wtime()
    call mpi_barrier(mpi_world,ierr)
    twl = mpi_wtime() -tw0

    gtrn(1:ndim) = 0d0
    tc0= mpi_wtime()
!.....TODO: allreduce may be redundant,  only reducing to node-0 is enough
!           if the minimization routine is wrtten so...
    call mpi_allreduce(gtrnl,gtrn,ndim,mpi_real8,mpi_sum,mpi_world,ierr)
    tcl= tcl +mpi_wtime() -tc0

!!$    gtrn(1:ndim)= gtrn(1:ndim) /swgt2trn

!.....Penalty
!!$    print '(a,20f10.4)','  gtrn=',gtrn(4:21)
    if( .not.allocated(gpena) ) allocate(gpena(ndim))
    call grad_penalty(ndim,x,gpena)
    gtrn(:) = gtrn(:) +gpena(:)
!!$    if( trim(cpot).eq.'uf3' .and. n_repul_pnts > 0 ) then
!!$      if( .not.allocated(grepul) ) allocate(grepul(ndim))
!!$      call calc_short_lossgrad(n_repul_pnts,repul_radii,drepul_tbl, &
!!$           ndim,grepul)
!!$      gtrn(:) = gtrn(:) +pwgt_repul*grepul(:)
!!$    endif
    call mask_grad(ndim,vranges,gtrn)

    if( .not.allocated(ldcover) ) then
      allocate(ldcover(ndim))
      dmem = dmem +4d0*size(ldcover)
      ldcover(:) = .false.
      do i=1,ndim
        if( abs(gtrn(i)) > 1d-14 ) ldcover(i) = .true.
      enddo
    endif

!.....only the bottle-neck times are taken into account
    tcg = tcl
    tgg = tgl
    twg = twl
    call mpi_reduce(tcl,tcg,1,mpi_real8,mpi_max,0,mpi_world,ierr)
    call mpi_reduce(tgl,tgg,1,mpi_real8,mpi_max,0,mpi_world,ierr)
    call mpi_reduce(twl,twg,1,mpi_real8,mpi_max,0,mpi_world,ierr)
    tcomm= tcomm +tcg
    tgrad= tgrad +tgg
    twait= twait +twg

    if( l1st ) then
      if( trim(cpot).eq.'uf3' ) then
        dmem = dmem + get_mem_uf3()
      else if( trim(cpot).eq.'uf3l') then
        dmem = dmem + get_mem_uf3l()
      endif
    endif
    
    if( trim(cpot).eq.'uf3') then
      call dealloc_gwx_uf3()
    endif
    l1st = .false.
    if( iprint > 10 ) print *,'(grad_w_pmd) done'
    return
  end subroutine grad_w_pmd
!=======================================================================
  subroutine pre_pmd(smpl,ndim,x,nff,cffs,rc,l1st)
!
!  Preprocesses before running pmd
!
    use variables,only: csmplfile,cpot,nsubff,csubffs,mdsys, &
         maxisp,nn_nl,nn_nhl,nn_sigtype,nn_asig,rc3, &
         interact,interact3,num_interact,iprint, &
         descs,nsf_desc,nsf2_desc,nsf3_desc,nsff_desc,ilsf2,ilsf3, &
         lcheby,cnst,wgtsp_desc,nspmax
    use parallel,only: myid_pmd,mpi_comm_pmd,nnode_pmd
    use pmdvars, only: specorder, am, naux
    use element,only: init_element, atom, get_element
    use DNN,only: set_paramsdir_DNN,set_params_DNN,set_actfunc_DNN
    use linreg,only: set_paramsdir_linreg,set_params_linreg
    use descriptor,only: set_paramsdir_desc,get_descs,get_ints,set_descs &
         ,lupdate_gsf,set_params_desc_new, lfitpot_desc => lfitpot
    use UF3,only: set_params_uf3, print_1b, print_2b, prm2s, n2b, &
         set_params_uf3l
    implicit none
    type(mdsys),intent(inout):: smpl
    integer,intent(in):: ndim, nff
    real(8),intent(in):: x(ndim)
    character(len=20),intent(in):: cffs(nff)
    real(8),intent(in):: rc
    logical,intent(in):: l1st

    integer:: i,is,nsf,nal,nnl,ndimt,ndim0,i2b
    character(len=128):: csmplname,ctype

    logical,external:: string_in_arr
    type(atom):: elem
    logical:: update_force_list
    real(8):: ptnsr(3,3),ekin
    character:: csp*3

!.....Sample dependent information
!.....Since at least one of FF requires mass infomation,
!     set mass info from specorder anyways.
    am(:) = 12d0
    specorder(:) = smpl%specorder
    do is=1,nspmax
      csp = specorder(is)
      if( trim(csp).ne.'x' ) then
        elem = get_element(trim(csp))
        am(is) = elem%mass
      endif
    enddo

    if( l1st .and. .not.allocated(smpl%aux) ) allocate(smpl%aux(naux,smpl%natm))

    do i=1,nff
      if( trim(cffs(i)).eq.'linreg' ) then
!.....Set lfitpot in descriptor module to let it know that it is called from fitpot
        lfitpot_desc = .true.
        call set_params_linreg(ndim,x)
      else if( trim(cffs(i)).eq.'dnn' ) then
!.....Set lfitpot in descriptor module to let it know that it is called from fitpot
        lfitpot_desc = .true.
        call set_params_DNN(ndim,x)
      else if( trim(cffs(i)).eq.'uf3' ) then
        call set_params_uf3(ndim,x)
      else if( trim(cffs(i)).eq.'uf3l' ) then
        if( iprint > 10 ) print *,'pre_pmd: into set_params_uf3l...'
        call set_params_uf3l(ndim,x)
      endif
    
      if( index(cffs(i),'NN').ne.0 .or. trim(cffs(i)).eq.'linreg' ) then
!.....Some potentials use descriptors already computed in the previous steps
!.....and stored in samples (and normalized if specified so.)
        if( .not. lupdate_gsf ) then
          nsf = smpl%nsf
          nal = smpl%nal
          nnl = smpl%nnl
!!$        call set_descs(nsf,nal,nnl,samples(ismpl)%gsf, &
!!$             samples(ismpl)%dgsf,samples(ismpl)%igsf)
          call set_descs(nsf,nal,nnl,smpl%gsf,smpl%dgsf,smpl%igsf)
        endif
      endif
    enddo  ! i=1,nff
    
    return
  end subroutine pre_pmd
!=======================================================================
  subroutine run_pmd(smpl,lgrad,lgrad_done,ndimp,epot,frcs,strs,rc &
       ,lfdsgnmat,gwe,gwf,gws)
!
!  Run pmd and get energy and forces of the system.
!
    use variables,only: mdsys,maxna,iprint,lematch,lfmatch,lsmatch, &
         maxisp,maxnf,nsp
    use parallel,only: myid_pmd,mpi_comm_pmd,nnode_pmd
    use force
    use descriptor,only: get_dsgnmat_force
    use ZBL,only: r_inner,r_outer
    use pmdvars, only: nspmax,naux,nstp,nx,ny,nz,specorder,am,dt,rbuf, &
         rc1nn,lvc
    use pmdvars,only: iprint_pmd => iprint, rc_pmd => rc
    use element
    implicit none
    include "../pmd/params_unit.h"
    type(mdsys),intent(inout):: smpl
    integer,intent(in):: ndimp
    real(8),intent(in):: rc
    real(8),intent(inout):: epot,frcs(3,maxna)
    real(8),intent(out):: strs(3,3)
    logical,intent(in):: lgrad,lgrad_done,lfdsgnmat
    real(8),intent(out),optional:: gwe(ndimp)
    real(8),intent(out),optional:: gwf(3,ndimp,maxnf)
    real(8),intent(out),optional:: gws(6,ndimp)

    logical,save:: l1st = .true.

    integer:: i,is
    real(8):: ptnsr(3,3),ekin
    character:: csp*3 
    type(atom):: elem
    logical:: update_force_list

    logical,external:: string_in_arr

!.....one_shot force calculation
!.....NOTE: unit of forces is eV/Ang. (not scaled by h-mat)
    if( iprint.gt.10 ) print *,'rum_pmd: into oneshot4fp...'
    call oneshot4fp(smpl%h0,smpl%h,smpl%natm,smpl%tag,smpl%ra, &
         smpl%va,frcs,smpl%strsi,smpl%eki,smpl%epi, &
         smpl%aux,ekin,epot,ptnsr,lgrad,lgrad_done,ndimp,maxisp, &
         gwe,gwf,gws,lematch,lfmatch,lsmatch,smpl%nfcal,smpl%lfrc_eval)
!!$!.....Stress definition, negative as compressive, positive as tensile,
!!$!     which is opposite in pmd. So multiply -1 to ptnsr and gws.
!!$!     But this should not be corrected here, rather before converting to sample data.
!!$    strs(1:3,1:3) = -ptnsr(1:3,1:3)
!!$    if( present(gws) ) gws(1:6,1:ndimp) = -gws(1:6,1:ndimp)
    strs(1:3,1:3) = ptnsr(1:3,1:3)
    if( lfdsgnmat ) call get_dsgnmat_force(smpl%dgsfa,mpi_comm_pmd)

    if( lvc ) smpl%charge_set = .true.

    return
  end subroutine run_pmd
!=======================================================================
  subroutine func_penalty(ndim,x,fp)
    use variables,only: cpot
    use UF3,only: calc_penalty_uf3, calc_penalty_uf3l
    integer,intent(in):: ndim
    real(8),intent(in):: x(ndim)
    real(8),intent(out):: fp

    integer:: i

    fp = 0d0
    if( trim(cpenalty).eq.'ridge' ) then
      do i=1,ndim
        fp = fp +x(i)*x(i)
      enddo
      fp = fp *penalty
    else if( trim(cpenalty).eq.'uf3' ) then
      if( trim(cpot).ne.'uf3' ) stop 'potential and penalty is not consistent !'
      call calc_penalty_uf3(ndim,x,pwgt2b,pwgt2bd,pwgt2bs, &
           pwgt3b,pwgt3bd,repul_radii,fp)
    else if( trim(cpenalty).eq.'uf3l' ) then
      if( trim(cpot).ne.'uf3l' ) stop 'potential and penalty is not consistent !'
      call calc_penalty_uf3l(ndim,x,pwgt2b,pwgt2bd,pwgt2bs, &
           pwgt3b,pwgt3bd,repul_radii,fp)
    endif
    return
  end subroutine func_penalty
!=======================================================================
  subroutine grad_penalty(ndim,x,gp)
    use UF3,only: calc_penalty_grad_uf3,calc_penalty_grad_uf3l
    integer,intent(in):: ndim
    real(8),intent(in):: x(ndim)
    real(8),intent(out):: gp(ndim)

    integer:: i

    gp(:) = 0d0
    if( trim(cpenalty).eq.'ridge' ) then
      gp(:) = 2d0*penalty*x(:)
    else if( trim(cpenalty).eq.'uf3' ) then
      call calc_penalty_grad_uf3(ndim,x,pwgt2b,pwgt2bd,pwgt2bs, &
           pwgt3b,pwgt3bd,repul_radii,gp)
    else if( trim(cpenalty).eq.'uf3l' ) then
      call calc_penalty_grad_uf3l(ndim,x,pwgt2b,pwgt2bd,pwgt2bs, &
           pwgt3b,pwgt3bd,repul_radii,gp)
    endif
    return
  end subroutine grad_penalty
!=======================================================================
  subroutine create_mpi_comm_pmd()
!
!  Create MPI COMM for pmd.
!  To create a MPI COMM on each node, first create a MPI GROUP for world,
!  then create a MPI GROUP for this node, and then create the MPI COMM
!  from the GROUP for this node. This is how to create sub communicator
!  in the MPI.
!
    use variables,only: iprint
    use parallel
    use pmdvars,only: myid_md,mpi_md_world,nodes_md
    implicit none

    integer:: iranks(1)
    integer:: mpi_group_pmd

    iranks(1) = myid

    call mpi_comm_split(mpi_world,myid,myid,mpi_comm_pmd,ierr)

    call mpi_comm_size(mpi_comm_pmd,nnode_pmd,ierr)
    call mpi_comm_rank(mpi_comm_pmd,myid_pmd,ierr)
    call mpi_comm_group(mpi_comm_pmd,mpi_group_pmd,ierr)
!.....Store these vars within pmdvars module
    nodes_md = nnode_pmd
    myid_md = myid_pmd
    mpi_md_world = mpi_comm_pmd

    if( myid.eq.0 .and. iprint.gt.0 ) then
      write(6,'(a)') ''
      write(6,'(a)') ' MPI_COMM_PMD was created at each node '// &
           'for pmd calculations.'
    endif
!!$    if( myid.eq.0 ) print *,'MPI_COMM values:'
!!$    print *,'  myid,mpi_world,mpi_comm_pmd=',myid,mpi_world,mpi_comm_pmd

  end subroutine create_mpi_comm_pmd
!=======================================================================
  subroutine subtract_FF_dataset()
!
!  Subtract energies, forces and stresses, obtained from data
!  of the same filename in dataset_subtract.
!  Compared to the original subract_FF() routine, this does not
!  call pmd and compute specific force-field.
!
    use variables,only: samples,dmem,maxna
    use parallel,only: myid,isid0,isid1
    implicit none

    integer:: i,ismpl,natm
    logical,save:: l1st = .true.
    real(8):: epot,strs(3,3)
    real(8),save,allocatable:: frcs(:,:)

    if( l1st ) then

      if( .not.allocated(frcs) ) then
        allocate(frcs(3,maxna))
        dmem = dmem +8d0*size(frcs)
      endif

!.....Only at the 1st call, perform pmd to get (esub,fsub,ssub)
      do ismpl=isid0,isid1
        natm = samples(ismpl)%natm
        samples(ismpl)%esub = epot
        samples(ismpl)%fsub(1:3,1:natm) = frcs(1:3,1:natm)
        samples(ismpl)%ssub(1:3,1:3) = strs(1:3,1:3)
      enddo

    endif


    l1st = .false.
    return
  end subroutine subtract_FF_dataset
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
    use dipole,only: set_paramsdir_dipole
    use Morse,only: set_paramsdir_Morse,set_params_Morse
    use LJ,only: set_paramsdir_LJ
    use ZBL,only: set_params_ZBL
    use Bonny_WRe,only: set_paramsdir_Bonny
    use cspline,only: set_paramsdir_cspline
    use force,only: loverlay
    implicit none

    integer:: i,ismpl,natm
    logical:: lgrad = .false.
    logical:: lgrad_done = .false.
    logical:: luse_Morse = .false.
    logical:: luse_Morse_repul = .false.
    logical:: luse_Coulomb = .false.
    logical:: luse_LJ_repul = .false.
    logical:: luse_ZBL = .false.
!!$    logical:: luse_Bonny_WRe = .false.
!!$    logical:: luse_cspline = .false.
!!$    logical:: luse_dipole = .false.
    logical:: lfdsgnmat = .false.  ! Not to compute dsgnmat for subtracted FFs
    logical,save:: l1st = .true.
    real(8):: epot,strs(3,3)
    real(8),save,allocatable:: frcs(:,:)

    if( l1st ) then
      if( myid.eq.0 .and. iprint.ne.0 ) then
        print '(/,a)',' Force field to be subtracted:'
        do i=1,nsubff
          print *,'  i,FF = ',i,trim(csubffs(i))
        enddo
      endif

      if( .not.allocated(frcs) ) then
        allocate(frcs(3,maxna))
        dmem = dmem +8d0*size(frcs)
      endif

!.....Only at the 1st call, perform pmd to get (esub,fsub,ssub)
      do ismpl=isid0,isid1
        natm = samples(ismpl)%natm
        call pre_pmd(samples(ismpl),nvars,vars,nsubff,csubffs,&
             rc_other,.true.)
        call run_pmd(samples(ismpl),lgrad,lgrad_done,nvars,&
             epot,frcs,strs,rc_other,lfdsgnmat)
!!$      print *,'myid,ismpl,epot=',myid,ismpl,epot
        samples(ismpl)%esub = epot
        samples(ismpl)%fsub(1:3,1:natm) = frcs(1:3,1:natm)
        samples(ismpl)%ssub(1:3,1:3) = strs(1:3,1:3)
!!$        print *,'ismpl,esub=',ismpl,epot
!!$        do i=1,natm
!!$          print *,'ia,fsub=',i,frcs(1:3,i)
!!$        enddo
        if( luse_ZBL ) then
          overlay = loverlay
        endif
      enddo  ! ismpl

    endif


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
  subroutine get_mean_gsf()
!
!  Compute the mean of input symmetric functions.
!  Called from normalize(), since the normalization usually requires mean of G's.
!
    use variables
    use parallel
    implicit none 
    real(8),parameter:: tiny = 1d-15
    integer:: nsuml,nsumg,ismpl,natm,isf,jsf,ia,nsf
    real(8):: gmeanl,tmp,gvarl,dgi,dgj
    real(8),allocatable:: gsfml(:),gsfvl(:),gsfcl(:,:),gsfvsq(:),gsfsl(:)
!    real(8),allocatable:: dgsft(:,:,)

    if( .not. allocated(samples(isid0)%gsf) ) then
      print *,'smpl%gsf is not allocated...'
      stop
    endif

    nsf = samples(isid0)%nsf
    if( .not.allocated(gsfml) ) then
      allocate(gsfml(nsf),gsfvl(nsf),gsfcl(nsf,nsf),&
           gsfvsq(nsf),gsfsl(nsf))
!!$      if( myid.eq.0 ) print *,'get_mean_gsf: dmem=',dmem
    endif
    if( .not. allocated(gsfms) ) then
      allocate(gsfms(nsf),gsfvs(nsf),gsfss(nsf),gsfcorr(nsf,nsf))
      dmem = dmem +8d0*size(gsfms) +8d0*size(gsfvs) +8d0*size(gsfss) &
           +8d0*size(gsfcorr)
!!$      if( myid.eq.0 ) print *,'get_mean_gsf: dmem=',dmem
    endif

!.....compute mean value
    gsfml(1:nsf) = 0d0
    gsfvl(1:nsf) = 0d0
    gsfsl(1:nsf) = 0d0
    gmeanl= 0d0
    gvarl = 0d0
    nsuml= 0
    do ismpl=isid0,isid1
      natm= samples(ismpl)%natm
      nsuml= nsuml +natm
      !.....sum up gsf
      do ia=1,natm
        do isf=1,nsf
          tmp = samples(ismpl)%gsf(isf,ia)
!!$          print *,'myid,ia,isf,tmp=',myid,ia,isf,tmp
          gsfml(isf) = gsfml(isf) +tmp
          gsfsl(isf) = gsfsl(isf) +tmp*tmp
          gmeanl= gmeanl +tmp
          gvarl = gvarl +tmp*tmp
        enddo
      enddo
    enddo
    nsumg= 0
    gsfms(1:nsf) = 0d0
    gsfvs(1:nsf) = 0d0
    gsfss(1:nsf) = 0d0
    call mpi_allreduce(gsfml,gsfms,nsf,mpi_real8 &
         ,mpi_sum,mpi_world,ierr)
    call mpi_allreduce(gsfvl,gsfvs,nsf,mpi_real8 &
         ,mpi_sum,mpi_world,ierr)
    call mpi_allreduce(gsfsl,gsfss,nsf,mpi_real8 &
         ,mpi_sum,mpi_world,ierr)
    call mpi_allreduce(nsuml,nsumg,1,mpi_integer &
         ,mpi_sum,mpi_world,ierr)
    do isf=1,nsf
      gsfms(isf)= gsfms(isf)/nsumg
      gsfss(isf)= gsfss(isf)/nsumg
      gsfvs(isf)= gsfss(isf) -gsfms(isf)**2
      gsfss(isf)= sqrt(gsfss(isf))
    enddo

!.....Correlation coefficients
    gsfcl(:,:) = 0d0
    do ismpl=isid0,isid1
      natm= samples(ismpl)%natm
      do ia=1,natm
        do isf=1,nsf-1
          dgi = samples(ismpl)%gsf(isf,ia) - gsfms(isf)
          do jsf=isf+1,nsf
            dgj = samples(ismpl)%gsf(jsf,ia) - gsfms(jsf)
            gsfcl(jsf,isf) = gsfcl(jsf,isf) +dgi*dgj
          enddo
        enddo
      enddo
    enddo
    gsfcorr(:,:) = 0d0
    call mpi_allreduce(gsfcl,gsfcorr,nsf*nsf,mpi_real8 &
         ,mpi_sum,mpi_world,ierr)
    do isf=1,nsf
      gsfvsq(isf) = sqrt(gsfvs(isf)*nsumg)
      gsfcorr(isf,isf) = 1d0
    enddo
    do isf=1,nsf-1
      do jsf=isf+1,nsf
        gsfcorr(jsf,isf) = gsfcorr(jsf,isf) &
               / gsfvsq(isf) /gsfvsq(jsf)
        gsfcorr(isf,jsf) = gsfcorr(jsf,isf)
      enddo
    enddo
!.....Output correlation
    if( myid.eq.0 .and. iprint.gt.2 ) then
      print *,'Write out.correlation'
      open(20,file='out.correlation',status='replace')
      do isf=1,nsf
        do jsf=1,nsf
          write(20,'(2i6,2es13.4e3)') isf,jsf,gsfcorr(isf,jsf)
        enddo
      enddo
      close(20)
    endif
    
    gsfmean= 0d0
    gsfvar = 0d0
    do isf=1,nsf
      gsfmean= gsfmean +gsfms(isf)
      gsfvar= gsfvar +gsfvs(isf)
    enddo
    gsfmean= gsfmean/nsf
    gsfvar = gsfvar/nsf
!!$    gsfvar = get_variance_input(gsfmean)
    if( myid.eq.0 .and. iprint.gt.0 ) then
      print *,''
      write(6,'(a,es12.3)') ' Mean of input symmetry functions = ',gsfmean
      write(6,'(a,es12.3)') ' Var  of input symmetry functions = ',gsfvar
      if( iprint.gt.1 ) then
        do isf=1,nsf
          write(6,'(a,i5,2es14.4e3)') '   isf,mean,variance= ' &
               ,isf,gsfms(isf),gsfvs(isf)
        enddo
      endif
    endif


    deallocate(gsfml,gsfvl,gsfcl,gsfvsq,gsfsl)
    return
  end subroutine get_mean_gsf
!=======================================================================
  subroutine write_dsgnmats()
!
!  Write design matrices of atoms, energies, forces.
!  Only available for nnode==1.
!
    use variables
    use parallel
    implicit none

    integer:: nsf,natm,nasum,ndat,ismpl,isf,ia,ixyz
    character(len=128):: cnum
    real(8):: gtmp

    nsf = samples(isid0)%nsf

    nasum = 0
    do ismpl=isid0,isid1
      natm= samples(ismpl)%natm
      nasum= nasum +natm
    enddo

    write(6,'(a)') ' Writing out.dsgnmat_atm... '
    write(cnum,'(i0)') nsf
    open(21,file='out.dsgnmat_atm')
    write(21,'(2i10)') nasum,nsf
    do ismpl=isid0,isid1
      natm = samples(ismpl)%natm
      do ia=1,natm
        write(21,'('//trim(cnum)//'es12.3e3)',advance='no') (samples(ismpl)%gsf(isf,ia),isf=1,nsf)
        write(21,'(1x,a)') trim(samples(ismpl)%csmplname)
      enddo
    enddo
    close(21)

!.....Design matrix for lasso ()
    if( lematch ) then
      write(6,'(a)') ' Writing out.dsgnmat_erg... '
      open(22,file='out.dsgnmat_erg')
      open(25,file='out.esubs')
!.....Since now it is only nnode==1, isid1==nsmpl
      write(22,'(a)') '# y_i, (x_{ij},j=1,nsf) of energy matching'
      write(22,'(2i8)') isid1,nsf
      write(25,'(a)') '# esub'
      write(25,'(2i8)') isid1
      do ismpl=isid0,isid1
        write(22,'(es12.3e3)',advance='no') samples(ismpl)%eref &
             -samples(ismpl)%esub
        write(25,'(es12.3e3)') samples(ismpl)%esub
        natm = samples(ismpl)%natm
        do isf=1,nsf
          gtmp = 0d0
          do ia=1,natm
            gtmp = gtmp +samples(ismpl)%gsf(isf,ia)
          enddo
          write(22,'(es12.3e3)',advance='no') gtmp
        enddo
        write(22,*) ''
      enddo
      close(22)
      close(25)
    endif

!.....Design matrix for force-matching
    if( lfmatch ) then
      write(6,'(a)') ' Writing out.dsgnmat_frc... '
      open(23,file='out.dsgnmat_frc')
      write(23,'(a)') '# y_i, (x_{ij},j=1,nsf) of force matching'
      open(26,file='out.fsubs')
      write(26,'(a)') '# fsub'
      ndat = 0
      do ismpl=isid0,isid1
        ndat = ndat +3*samples(ismpl)%natm
      enddo
      write(23,'(2i8)') ndat,nsf
      write(26,'(2i8)') ndat
      do ismpl=isid0,isid1
        natm = samples(ismpl)%natm
        do ia=1,natm
          do ixyz=1,3
            write(23,'(es12.3e3)',advance='no') samples(ismpl)%fref(ixyz,ia) &
                 -samples(ismpl)%fsub(ixyz,ia)
            write(26,'(es12.3e3)') samples(ismpl)%fsub(ixyz,ia)
            do isf=1,nsf
              write(23,'(es12.3e3)',advance='no') -samples(ismpl)%dgsfa(ixyz,isf,ia)
            enddo
            write(23,*) ''
          enddo
        enddo
      enddo
      close(23)
      close(26)
    endif

  end subroutine write_dsgnmats
!=======================================================================
  subroutine normalize()
    use variables,only: lnormalized,cnormalize,iprint
    use parallel
    implicit none

    logical,save:: l1st = .true.

!.....If already done, skip
    if( lnormalized ) return

    if( l1st ) call get_mean_gsf()

    if( cnormalize(1:3).eq.'std' .or. cnormalize(1:8).eq.'standard' .or. &
        cnormalize(1:3).eq.'var' ) then
      if( myid.eq.0 .and. iprint.ne.0 .and. l1st ) &
           print *,'Normalize descriptors wrt standard deviation.'
      call normalize_std()
    else if( cnormalize(1:4).eq.'norm' ) then
      if( myid.eq.0 .and. iprint.ne.0 .and. l1st ) &
           print *,'Normalize descriptors wrt norm.'
      call normalize_norm()
    else if( cnormalize(1:4).eq.'none' ) then
      if( myid.eq.0 .and. iprint.ne.0 .and. l1st ) &
           print *,'No normalization of descriptors.'
    else
      if( myid.eq.0 .and. iprint.ne.0 .and. l1st ) &
           print *,'WARNING: no such normalization, '//trim(cnormalize)
    endif
    l1st = .false.
  end subroutine normalize
!=======================================================================
  subroutine normalize_std()
!
!  Normalize inputs (descriptors) wrt standard deviation.
!
    use variables, only: samples,nvars,vars,vranges&
         ,lnormalized,cpot,gsfvar,gsfvs,sgms,sgmis,sgm_min,iprint &
         ,nn_nhl
    use parallel
    use descriptor,only: set_gscale
    implicit none
    integer:: ismpl,natm,isf,i,iv,ihl0,ihl1
    integer,save:: nsf
    real(8):: sgm,sgmi,sgmax,sgmin
    logical,save:: l1st = .true.

    if( l1st ) then
      nsf = samples(isid0)%nsf
      if( .not. allocated(sgms) ) allocate(sgms(nsf),sgmis(nsf))
      sgmax= 0d0
      do isf=1,nsf
        sgmax= max(gsfvs(isf),sgmax)
      enddo
      sgmax= sqrt(sgmax)
      sgmin= sgmax*sgm_min
      do isf=1,nsf
        sgms(isf) = max(sqrt(gsfvs(isf)),sgmin)
        sgmis(isf)= 1d0/sgms(isf)
      enddo
      sgm = sqrt(gsfvar)
      sgmi= 1d0/sgm
!.....Set gscale in descriptor module
      call set_gscale(nsf,sgmis)
!.....AND scale G's in each sample as well
      do ismpl=isid0,isid1
        natm= samples(ismpl)%natm
        if( .not. allocated(samples(ismpl)%gsf) ) then
          print *,'ERROR: gsf not allocated, myid,ismpl,csmplname='&
               ,myid,ismpl,trim(samples(ismpl)%csmplname)
          stop
        endif
        do isf=1,nsf
          samples(ismpl)%gsf(isf,:)= samples(ismpl)%gsf(isf,:) *sgmis(isf)
!!$          samples(ismpl)%dgsf(:,isf,:,:)= &
!!$               samples(ismpl)%dgsf(:,isf,:,:) *sgmis(isf)
        enddo
      enddo
      if( myid.eq.0 .and. iprint.gt.1 ) then
        do isf=1,nsf
          print '(a,i5,2es12.4)','   isf,sgm,sgmi=',isf,sgms(isf),sgmis(isf)
        enddo
      endif
    endif  ! l1st

    if( trim(cpot).eq.'linreg' ) then
      do i=1,nvars
        vars(i) = vars(i) *sgms(i)
        vranges(1:2,i) = vranges(1:2,i) *sgms(i)
      enddo
    else if( trim(cpot).eq.'dnn' ) then
      iv = 0
      do ihl1=1,nn_nhl(1)
        do ihl0=0,nn_nhl(0)
          iv = iv + 1
          if( ihl0.eq.0 ) cycle  ! Care about bias node
          vars(iv) = vars(iv) *sgms(ihl0)
          vranges(1:2,iv) = vranges(1:2,iv) *sgms(ihl0)
        enddo
      enddo
    endif

    l1st = .false.
    lnormalized = .true.
    
  end subroutine normalize_std
!=======================================================================
  subroutine normalize_norm()
!
!  Normalize inputs (descriptors)
!
    use variables, only: samples,nvars,vars,vranges&
         ,lnormalized,cpot,sgms,sgmis,gsfss,sq_min,iprint &
         ,nn_nhl
    use parallel
    use descriptor,only: set_gscale
    implicit none
    integer:: ismpl,natm,isf,i,iv,ihl0,ihl1
    real(8):: sqmax,sqmin
    integer,save:: nsf
    logical,save:: l1st = .true.

    if( l1st ) then
      nsf = samples(isid0)%nsf
      if( .not. allocated(sgms) ) allocate(sgms(nsf),sgmis(nsf))
      sqmax = 0d0
      do isf=1,nsf
        sqmax = max(gsfss(isf),sqmax)
      enddo
      sqmin = sqmax*sq_min
      do isf=1,nsf
        sgms(isf) = max(gsfss(isf),sqmin)
        sgmis(isf)= 1d0/sgms(isf)
      enddo
!.....standardize G values
      do ismpl=isid0,isid1
        natm= samples(ismpl)%natm
        if( .not. allocated(samples(ismpl)%gsf) ) then
          print *,'ERROR: gsf not allocated, myid,ismpl,csmplname='&
               ,myid,ismpl,trim(samples(ismpl)%csmplname)
          stop
        endif
!.....Set gscale in descriptor module
        call set_gscale(nsf,sgmis)
!.....AND scale G's in each sample as well
        do isf=1,nsf
          samples(ismpl)%gsf(isf,:)= samples(ismpl)%gsf(isf,:) *sgmis(isf)
!!$          samples(ismpl)%dgsf(:,isf,:,:)= &
!!$               samples(ismpl)%dgsf(:,isf,:,:) *sgmis(isf)
        enddo
      enddo
      if( myid.eq.0 .and. iprint.gt.1 ) then
        do isf=1,nsf
          print '(a,i5,2es12.4)','   isf,sgm,sgmi=',isf,sgms(isf),sgmis(isf)
        enddo
      endif
    endif

    if( trim(cpot).eq.'linreg' ) then
      do i=1,nvars
        vars(i) = vars(i) *sgms(i)
        vranges(1:2,i) = vranges(1:2,i) *sgms(i)
      enddo
!!$    else if( trim(cpot).eq.'NN2' ) then
!!$      iv = 0
!!$      do ihl0=1,nn_nhl(0)
!!$        do ihl1=1,nn_nhl(1)  ! NN2 does not use bias...
!!$          iv = iv + 1
!!$          vars(iv) = vars(iv) *sgms(ihl0)
!!$          vranges(1:2,iv) = vranges(1:2,iv) *sgms(ihl0)
!!$        enddo
!!$      enddo
    else if( trim(cpot).eq.'dnn' ) then
      iv = 0
      do ihl1=1,nn_nhl(1)
        do ihl0=0,nn_nhl(0)
          iv = iv + 1
          if( ihl0.eq.0 ) cycle  ! Care about bias node
          vars(iv) = vars(iv) *sgms(ihl0)
          vranges(1:2,iv) = vranges(1:2,iv) *sgms(ihl0)
        enddo
      enddo
    endif

    l1st = .false.
    lnormalized = .true.
    
  end subroutine normalize_norm
!=======================================================================
  subroutine restore_normalize()
!
!  Restore weights by inverse normalization
!
    use variables, only: nvars,vars,vranges&
         ,lnormalized,cnormalize,cpot,sgmis, nn_nhl
    use parallel
    implicit none
    integer:: i,iv,ihl0,ihl1
    real(8):: sgmi

    if( .not. lnormalized ) return

    if( cnormalize(1:3).eq.'std' .or. cnormalize(1:3).eq.'var' .or. &
         cnormalize(1:4).eq.'norm' ) then
      if( trim(cpot).eq.'linreg' ) then
        do i=1,nvars
          vars(i) = vars(i) *sgmis(i)
          vranges(:,i) = vranges(:,i) *sgmis(i)
        enddo
!!$      else if( trim(cpot).eq.'NN2' ) then
!!$        iv = 0
!!$        do ihl0=1,nn_nhl(0)
!!$          sgmi = sgmis(ihl0)
!!$          do ihl1=1,nn_nhl(1)
!!$            iv= iv + 1
!!$            vars(iv)= vars(iv) *sgmi
!!$            vranges(1:2,iv)= vranges(1:2,iv) *sgmi
!!$          enddo
!!$        enddo
      else if( trim(cpot).eq.'dnn' ) then
        iv = 0
        do ihl1=1,nn_nhl(1)
          do ihl0=0,nn_nhl(0)
            iv= iv + 1
            if( ihl0.eq.0 ) cycle ! Care about bias node
            sgmi = sgmis(ihl0)
            vars(iv)= vars(iv) *sgmi
            vranges(1:2,iv)= vranges(1:2,iv) *sgmi
          enddo
        enddo
      endif
    endif
    
    lnormalized = .false.
  end subroutine restore_normalize
!=======================================================================
  subroutine compute_gdw()
!
!  Compute Gaussian density weight (GDW).
!  The Theta function in Ref. [1] almost equals to Theta(x) = x.
!  Ref:
!    [1] W. Jeong, et al., Journal of Physical Chemistry C, 122(39), (2019) 22790–22795
!
    use variables,only: mdsys,samples,gdsgm,lgdwed,maxnin,natot
    use parallel
    implicit none

    type(mdsys):: smpl
    integer:: inode,itag,iprty,nn,ndims,ndimr,inc,ismpl,ia,ja,isf,nsf
    real(8):: dg2,gdwmin,gdwmax
    real(8),allocatable:: xparts(:,:),xpartr(:,:),gdfs(:),gdfr(:)

    if( lgdwed ) return

    nsf = samples(isid0)%nsf
    allocate(xparts(nsf,maxnin),xpartr(nsf,maxnin))
    allocate(gdfs(maxnin),gdfr(maxnin))

!.....Prepare xpartr by copying that of the current node
    inc = 0
    do ismpl=isid0,isid1
      smpl = samples(ismpl)
      do ia=1,smpl%natm
        inc = inc +1
        xpartr(1:nsf,inc) = smpl%gsf(1:nsf,ia)
      enddo
    enddo
    ndimr = inc  ! total num of atoms in the current node
    gdfr(:) = 0d0
    
!.....Repeat sending, recving, and computing the GDF
    inode = mod(myid+1,nnode)
    itag = 10
    iprty = mod(myid,2)
    do nn=1,nnode

!.....Copy recv data to the data to be sent
      gdfs(:) = gdfr(:)
      xparts(1:nsf,1:ndimr) = xpartr(1:nsf,1:ndimr)
      ndims = ndimr

!.....Send and recv data if needed, in case of nn==1, it is not needed
      if( nn.ne.1 ) then
        call mespasi(inode,iprty,ndims,ndimr,1,1,itag,mpi_world)
        call mespasd(inode,iprty,xparts,xpartr,ndims*nsf,ndimr*nsf,itag,mpi_world)
        call mespasd(inode,iprty,gdfs,gdfr,ndims,ndimr,itag,mpi_world)
      endif
      
!.....Compute GDF between given data and data in the current node
      do ismpl=isid0,isid1
        smpl = samples(ismpl)
        do ja=1,ndimr
          do ia=1,smpl%natm
            dg2 = 0d0
            do isf=1,nsf
              dg2 = dg2 +(xpartr(isf,ja) -smpl%gsf(isf,ia))**2
            enddo
            gdfr(ja) = gdfr(ja) +exp(-dg2/nsf/2/gdsgm)
          enddo
        enddo
      enddo

    enddo ! nn=1,nnode

!.....Now gdfr(:) has complete contributions from all the atoms,
!.....so send it back to the original node
    ndims = ndimr
    gdfs(:) = gdfr(:)
    if( nnode.gt.1 ) then
      call mespasi(inode,iprty,ndims,ndimr,1,1,itag,mpi_world)
      call mespasd(inode,iprty,gdfs,gdfr,ndims,ndimr,itag,mpi_world)
    endif
    
!.....Store GDFs and compute GDWs
    inc = 0
    do ismpl=isid0,isid1
      do ia=1,samples(ismpl)%natm
        inc = inc + 1
        samples(ismpl)%gdf(ia) = gdfr(inc)/natot
        samples(ismpl)%gdw(ia) = 1d0/samples(ismpl)%gdf(ia)
      enddo
    enddo

    if( myid.eq.0 ) then
      gdwmin = 1d0 /maxval(gdfr)*natot
      gdwmax = 1d0 /minval(gdfr)*natot
      print '(a,2es12.4)',' Gaussian density weights, min and max = ', &
           gdwmin,gdwmax
    endif
    deallocate(xparts,xpartr,gdfs,gdfr)
    lgdwed = .true.
    
  end subroutine compute_gdw
!=======================================================================
  function ndat_in_line(ionum,delim) result(ndat)
    use util, only: num_data
    implicit none
    integer,intent(in):: ionum
    character(len=1),intent(in):: delim
    integer:: ndat
!!$    integer,external:: num_data
    character:: ctmp*128

    read(ionum,'(a)') ctmp
    ndat = num_data(trim(ctmp),delim)
    backspace(ionum)
    return

  end function ndat_in_line
!=======================================================================
end module fp_common
!-----------------------------------------------------------------------
! Local Variables:
! compile-command: "make fitpot"
! End:
