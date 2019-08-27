module fp_common
!-----------------------------------------------------------------------
!                     Last modified: <2019-08-27 16:44:49 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
!
! Module that contains common functions/subroutines for fitpot.
!
  implicit none
  save
  real(8),allocatable:: fdiff(:,:),frcs(:,:),gtrnl(:)
  real(8),allocatable:: gwe(:),gwf(:,:,:),gws(:,:)
  real(8),allocatable:: gwesub(:)
  real(8):: pdiff(6), ptnsr(3,3)
  real(8):: epotsub

  logical:: fp_common_initialized= .false.

  integer,parameter:: ivoigt(3,3)= &
       reshape((/ 1, 6, 5, 6, 2, 4, 5, 4, 3 /),shape(ivoigt))

!.....Store loverlay and r_inner/outer arrays
  logical:: overlay
  real(8),allocatable:: ri_zbl(:),ro_zbl(:)
  
contains
!=======================================================================
  subroutine init()
    use variables,only: swgt2trn, swgt2tst, samples,lematch,lfmatch,lsmatch
    use parallel

    integer:: ismpl,nterms
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
!!$      print *,'ismpl,wgt,swgtrn=',ismpl,samples(ismpl)%wgt,swgtrn
    enddo
    swgt2trn = 0d0
    swgt2tst = 0d0
    call mpi_allreduce(swgtrn,swgt2trn,1,mpi_real8,mpi_sum &
         ,mpi_world,ierr)
    call mpi_allreduce(swgtst,swgt2tst,1,mpi_real8,mpi_sum &
         ,mpi_world,ierr)
    nterms = 0
    if( lematch ) nterms = nterms + 1
    if( lfmatch ) nterms = nterms + 1
    if( lsmatch ) nterms = nterms + 1
    swgt2trn = swgt2trn*nterms
    swgt2tst = swgt2tst*nterms
    if( myid.eq.0 ) then
      write(6,'(a)') ' Weights to divide loss function:'
      write(6,'(a,f10.1)') '   for training: ',swgt2trn
      write(6,'(a,f10.1)') '   for test:     ',swgt2tst
    endif

    fp_common_initialized = .true.

  end subroutine init
!=======================================================================
  subroutine func_w_pmd(ndim,x,ftrn,ftst)
!
!  Evaluate loss function value using pmd (actually one_shot routine.)
!
    use variables,only:nsmpl,nsmpl_trn,samples,nprcs,tfunc &
         ,lematch,lfmatch,lsmatch,nfunc,tcomm,mdsys,erefmin &
         ,cmaindir,cevaltype,swgt2trn,swgt2tst,cpot &
         ,nff,cffs,nsubff,csubffs,cmaindir,maxna,rcut,rc3 &
         ,crefstrct,erefsub,myidrefsub,isidrefsub,iprint,maxisp &
         ,nn_nl,nn_nhl,nn_sigtype,ctype_loss,interact,memgsf,cfmethod
    use parallel
    use minimize
    use descriptor,only: lupdate_gsf,get_descs,get_ints
    implicit none
    integer,intent(in):: ndim
    real(8),intent(in):: x(ndim)
    real(8),intent(out):: ftrn,ftst

    integer:: ismpl,natm,ia,ixyz,jxyz,idim,k,nsf,nal,nnl,isf,i,ndimt
    real(8):: dn3i,ediff,fscale,eref,epot,swgt,wgtidv,esub
    real(8):: eerr,ferr,ferri,serr,serri,strs(3,3),fref
    real(8):: ftrnl,ftstl,ftmp
    real(8):: edenom,fdenom
    real(8):: tfl,tcl,tfg,tcg,tf0,tc0
    type(mdsys):: smpl
    logical,save:: l1st = .true.
    logical,parameter:: lcalcgrad = .false.
    logical:: lfdsgnmat
    character(len=128):: cdirname,ctype

    logical,external:: string_in_arr

    nfunc= nfunc +1

    tc0= mpi_wtime()
    call mpi_bcast(x,ndim,mpi_real8,0,mpi_world,ierr)
    tcl= mpi_wtime() -tc0
    tf0= mpi_wtime()

    lfdsgnmat = .false.  ! Initialize lfdsgnmat

    if( l1st ) then
      if( .not.fp_common_initialized ) call init()
      if( .not.allocated(fdiff) ) allocate(fdiff(3,maxna),frcs(3,maxna))
      if( trim(cpot).eq.'NN2' .or. trim(cpot).eq.'linreg') lupdate_gsf = .true.
    endif

!!$    print *,'func_w_pmd: 02'
    if( .not. lematch .and. .not.lfmatch .and. .not.lsmatch ) then
      if( myid.eq.0 ) then
        print *,'Nothing to be fitted.'
      endif
      call mpi_finalize(ierr)
      stop
    endif

    print *,'myid,ndim=',myid,ndim
    print *,'myid,x(:)=',myid,x(1:ndim)
!!$    print *,'myid,isid0,isid1=',myid,isid0,isid1
!!$    print *,'func_w_pmd: 03'
    ftrnl = 0d0
    ftstl = 0d0
    do ismpl=isid0,isid1
      natm= samples(ismpl)%natm
      cdirname = trim(samples(ismpl)%cdirname)
      call ready4pmd(ismpl,ndim,x)

!.....Set lfdsgnmat=.true. to make run_pmd() compute dsgnmat_force related data
      if( trim(cpot).eq.'linreg' .and. &
           l1st .and. lfmatch .and. trim(cfmethod).eq.'dsgnmat' ) then
          lfdsgnmat = .true.
      endif
      
      call run_pmd(samples(ismpl),lcalcgrad,ndim,nff,cffs,epot,frcs,strs &
           ,rcut,lfdsgnmat)
      samples(ismpl)%epot = epot
      samples(ismpl)%fa(1:3,1:natm) = frcs(1:3,1:natm)
      samples(ismpl)%strs(1:3,1:3) = strs(1:3,1:3)
      if( trim(cpot).eq.'linreg' .or. trim(cpot).eq.'NN2' ) then
        call get_ints(nsf,nal,nnl)
        samples(ismpl)%nsf = nsf
        samples(ismpl)%nal = nal
        samples(ismpl)%nnl = nnl
        if( .not. allocated(samples(ismpl)%gsf) ) then
          memgsf = memgsf +nsf*nal +3*nsf*(nnl+1)*nal +nsf*(nnl+1) +nal
!!$          if( iprint.gt.1 ) print '(a,2i5,f9.2)',' myid,ismpl,memgsf(MB) = ' &
!!$               ,myid,ismpl,memgsf
          allocate(samples(ismpl)%gsf(nsf,nal) &
               ,samples(ismpl)%dgsf(3,nsf,0:nnl,nal) &
               ,samples(ismpl)%igsf(nsf,0:nnl,nal) )
        endif
!!$        print *,'func: cdirname,nsf,nal,nnl=',trim(cdirname),nsf,nal,nnl
        call get_descs(nsf,nal,nnl,samples(ismpl)%gsf &
             ,samples(ismpl)%dgsf,samples(ismpl)%igsf)
      endif
!!$      if( trim(cpot).eq.'NN2' ) then
!!$        if( .not. allocated(samples(ismpl)%hl1) ) then
!!$          allocate(samples(ismpl)%hl1(nhl(1),nal))
!!$        else if( size(samples(ismpl)%hl1) .ne. nhl(1)*nal ) then
!!$          deallocate(samples(ismpl)%hl1)
!!$          allocate(samples(ismpl)%hl1(nhl(1),nal))
!!$        endif
!!$        call get_NN2_hl1(samples(ismpl)%hl1)
!!$        print *,'hl1 after get_NN2_hl1:'
!!$        do i=1,nhl(1)
!!$          print *,'i,hl1(i,1)=',samples(ismpl)%hl1(i,1)
!!$        enddo
!!$      endif
    enddo

    if( len(trim(crefstrct)).gt.5 ) then
      if( myid.eq.myidrefsub ) then
        epotsub = samples(isidrefsub)%epot +samples(isidrefsub)%esub
        epotsub = epotsub /samples(isidrefsub)%natm
      endif
      call mpi_bcast(epotsub,1,mpi_real8,myidrefsub,mpi_world,ierr)
!!$      print *,'myid,epotsub=', myid,epotsub
    endif

    do ismpl=isid0,isid1
      smpl = samples(ismpl)
      cdirname= smpl%cdirname
      natm = smpl%natm
      epot = smpl%epot
      ftmp = 0d0
      swgt = smpl%wgt
!.....Energy matching
      if( lematch ) then
        eref= smpl%eref
        esub= smpl%esub
        eerr = smpl%eerr
        if( len(trim(crefstrct)).gt.5 ) then
          ediff= (epot-epotsub*natm+esub -(eref-erefsub*natm))/natm /eerr
        else
          ediff= (epot+esub -eref)/natm /eerr
        endif
!!$        print *,'ismpl,epot,esub,eref,ediff,natm,eerr,swgt=' &
!!$             ,ismpl,epot,esub,eref,ediff,natm,eerr,swgt
        if( trim(ctype_loss).eq.'Huber' ) then
          if( abs(ediff).gt.1.d0 ) then
            ediff = 2d0*abs(ediff) -1d0
          else
            ediff = ediff *ediff
          endif
        else ! LS as default
          ediff= ediff*ediff
        endif
        ftmp= ftmp +ediff *swgt
        if( iprint.gt.2 ) then
          write(6,'(a,2i4,1x,a,7es11.3)') ' myid,ismpl,cdirname,epot,eref,esub,(epot+esub)/natm= ', &
               myid,ismpl,trim(cdirname),epot,eref,esub,(epot+esub)/natm
        endif
      endif
!.....Force matching
      if( lfmatch .and. smpl%nfcal.ne.0 ) then
        frcs(1:3,1:natm) = smpl%fa(1:3,1:natm)
        ferr = smpl%ferr
        ferri = 1d0/ferr
        dn3i = 1d0/3/smpl%nfcal
        do ia=1,natm
          if( smpl%ifcal(ia).eq.0 ) cycle
          do ixyz=1,3
            fdiff(ixyz,ia)= (frcs(ixyz,ia)+smpl%fsub(ixyz,ia) &
                 -(smpl%fref(ixyz,ia))) *ferri
!!$            fref = smpl%fref(ixyz,ia)
!!$            fdiff(ixyz,ia)= (frcs(ixyz,ia)+smpl%fsub(ixyz,ia) &
!!$                 -fref) /(abs(fref) +ferr)
            if( trim(ctype_loss).eq.'Huber' ) then
              if( abs(fdiff(ixyz,ia)).gt.1d0 ) then
                fdiff(ixyz,ia) = 2d0*abs(fdiff(ixyz,ia)) -1d0
              else
                fdiff(ixyz,ia)= fdiff(ixyz,ia)*fdiff(ixyz,ia)
              endif
            else ! LS as default
              fdiff(ixyz,ia)= fdiff(ixyz,ia)*fdiff(ixyz,ia)
            endif
            ftmp= ftmp +fdiff(ixyz,ia) *dn3i *swgt
          enddo
        enddo
!!$        write(6,'(a,i6,es12.4)') 'ismpl,ftmp=',ismpl,ftmp
      endif

!.....Stress matching
      if( lsmatch ) then
!.....Compare these ptnsr elements with sref elements
        serr = smpl%serr
        serri = 1d0/serr
        pdiff(1:6) = 0d0
        do ixyz=1,3
          do jxyz=ixyz,3
            k = ivoigt(ixyz,jxyz)
            pdiff(k)= pdiff(k) +(smpl%strs(ixyz,jxyz) +smpl%ssub(ixyz,jxyz) &
                 -smpl%sref(ixyz,jxyz)) *serri
          enddo
        enddo
        if( trim(ctype_loss).eq.'Huber' ) then
          do k=1,6
            if( abs(pdiff(k)).gt.1d0 ) then
              pdiff(k) = 2d0*abs(pdiff(k)) -1d0
            else
              pdiff(k)= pdiff(k)*pdiff(k)
            endif
            ftmp= ftmp +pdiff(k) *swgt /6
          enddo
        else  ! LS as default
          do k=1,6
            pdiff(k)= pdiff(k)*pdiff(k)
            ftmp= ftmp +pdiff(k) *swgt /6
          enddo
        endif
      endif  ! stress matching

      if( smpl%iclass.eq.1 ) then
        ftrnl = ftrnl +ftmp
      else if( smpl%iclass.eq.2 ) then
        ftstl = ftstl +ftmp
      endif
!!$      write(6,'(a,f12.5)') ' ftrnl = ',ftrnl
    enddo  ! ismpl

!!$    call mpi_barrier(mpi_world,ierr)
!    tfunc= tfunc +mpi_wtime() -tf0
    tfl = mpi_wtime() -tf0

!!$    print *,'func_w_pmd: 07'
!!$    print *,'myid,ftrnl=',myid,ftrnl
    tc0= mpi_wtime()
    ftrn= 0d0
    ftst = 0d0
    call mpi_allreduce(ftrnl,ftrn,1,mpi_real8,mpi_sum,mpi_world,ierr)
    call mpi_allreduce(ftstl,ftst,1,mpi_real8,mpi_sum,mpi_world,ierr)
    ftrn = ftrn /swgt2trn
    if( swgt2tst.gt.1d-5 ) then
      ftst = ftst /swgt2tst
    endif
    tcl = tcl + (mpi_wtime() -tc0)

!!$    print *,'myid,ftrn=',myid,ftrn
!.....only the bottle-neck times are taken into account
    call mpi_reduce(tcl,tcg,1,mpi_real8,mpi_max,0,mpi_world,ierr)
    call mpi_reduce(tfl,tfg,1,mpi_real8,mpi_max,0,mpi_world,ierr)
    tcomm= tcomm +tcg
    tfunc= tfunc +tfg
!!$    print *,'myid,tfunc=',myid,tfunc
    l1st = .false.
    if( trim(cpot).eq.'NN2' .or. trim(cpot).eq.'linreg' ) lupdate_gsf = .false.

  end subroutine func_w_pmd
!=======================================================================
  subroutine grad_w_pmd(ndim,x,gtrn)
!
!  Evaluate the gradient of loss function value
!  using pmd (actually one_shot routine.)
!
    use variables,only: nsmpl,nsmpl_trn,tgrad,ngrad,tcomm,tgrad &
         ,samples,mdsys,swgt2trn,swgt2tst,cpot,nff,cffs,nsubff,csubffs &
         ,cmaindir,maxna,lematch,lfmatch,lsmatch,erefsub,crefstrct &
         ,rcut,rc3,myidrefsub,isidrefsub,iprint,maxisp,gscl &
         ,nn_nl,nn_nhl,nn_sigtype,ctype_loss,interact
    use parallel
    use minimize
!!$    use Coulomb,only: set_paramsdir_Coulomb,set_params_Coulomb
!!$    use Morse,only: set_paramsdir_Morse,set_params_vcMorse,set_params_Morse
!!$    use BMH,only: set_paramsdir_BMH,set_params_BMH
!!$    use Abell,only: set_paramsdir_Abell,set_params_Abell
!!$    use fpc,only: set_paramsdir_fpc,set_params_fpc
!!$    use EAM,only: set_paramsdir_EAM,set_params_EAM
!!$    use NN,only: set_paramsdir_NN,set_params_NN
!!$    use NN2,only: set_paramsdir_NN2,set_params_NN2,set_NN2_hl1
!!$    use linreg,only: set_paramsdir_linreg,set_params_linreg
!!$    use descriptor,only: set_paramsdir_desc,set_descs
    implicit none
    integer,intent(in):: ndim
    real(8),intent(in):: x(ndim)
    real(8),intent(out):: gtrn(ndim)

    integer:: ismpl,i,idim,natm,k,ia,ixyz,jxyz,nsf,nal,nnl,ndimt
    real(8):: tcl,tgl,tcg,tgg,tc0,tg0,epot,esub,strs(3,3),dn3i
    real(8):: ediff,eerr,eref,swgt,ferr,ferri,serr,serri,fref,tmp
    type(mdsys):: smpl
    logical,parameter:: lcalcgrad = .true.
    logical,parameter:: lfdsgnmat = .false.
    character(len=128):: cdirname,ctype

    logical,external:: string_in_arr

    if( .not.allocated(gtrnl) ) allocate(gtrnl(ndim))
    if( .not.allocated(gwe) ) allocate(gwe(ndim),gwf(ndim,3,maxna)&
         ,gws(ndim,6))
    if( len(trim(crefstrct)).gt.5 ) then
      if( .not.allocated(gwesub) ) allocate(gwesub(ndim))
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

    do ismpl=isid0,isid1
      smpl = samples(ismpl)
      natm= smpl%natm
      cdirname = smpl%cdirname
      if( iprint.gt.20 ) print *,'myid,ismpl,cdirname=',myid,ismpl,trim(cdirname)
!.....Since g calc is time consuming,
!.....not calculate g for test set.
      if( smpl%iclass.ne.1 ) cycle
      call ready4pmd(ismpl,ndim,x)
      
!!$      if( trim(cpot).eq.'vcMorse' ) then
!!$        call set_paramsdir_Morse(trim(cmaindir)//'/'//trim(cdirname)&
!!$             //'/pmd')
!!$        call set_paramsdir_Coulomb(trim(cmaindir)//'/'//trim(cdirname)&
!!$             //'/pmd')
!!$        call set_params_vcMorse(ndim,x)
!!$      else if( trim(cpot).eq.'Morse' ) then
!!$        call set_paramsdir_Morse(trim(cmaindir)//'/'//trim(cdirname)&
!!$             //'/pmd')
!!$        if( string_in_arr('screened_Coulomb',nsubff,csubffs) .or. &
!!$            string_in_arr('Coulomb',nsubff,csubffs) ) then
!!$          ctype = 'BVS'
!!$        else if( string_in_arr('Ewald_long',nsubff,csubffs) .or.&
!!$             string_in_arr('Ewald',nsubff,csubffs) ) then
!!$          ctype = 'full_Morse'
!!$        else
!!$          ctype = 'full_Morse'
!!$        endif
!!$        call set_params_Morse(ndim,x,ctype,interact)
!!$      else if( trim(cpot).eq.'BMH' ) then
!!$        call set_paramsdir_BMH(trim(cmaindir)//'/'//trim(cdirname)&
!!$             //'/pmd')
!!$        call set_params_BMH(ndim,x,cpot,interact)
!!$      else if( trim(cpot).eq.'Abell' ) then
!!$        call set_paramsdir_Abell(trim(cmaindir)//'/'//trim(cdirname)&
!!$             //'/pmd')
!!$        call set_params_Abell(ndim,x,cpot,interact)
!!$      else if( trim(cpot).eq.'fpc' ) then
!!$        call set_paramsdir_fpc(trim(cmaindir)//'/'//trim(cdirname)&
!!$             //'/pmd')
!!$        call set_paramsdir_Coulomb(trim(cmaindir)//'/'//trim(cdirname)&
!!$             //'/pmd')
!!$        call set_params_Coulomb(1,x(1),cpot, &
!!$             samples(ismpl)%specorder)
!!$        call set_params_fpc(ndim-1,x(2:ndim),cpot,interact)
!!$      else if( trim(cpot).eq.'EAM' ) then
!!$        call set_paramsdir_EAM(trim(cmaindir)//'/'//trim(cdirname)&
!!$             //'/pmd')
!!$        call set_params_EAM(ndim,x)
!!$      else if( trim(cpot).eq.'NN' ) then
!!$        call set_paramsdir_NN(trim(cmaindir)//'/'//trim(cdirname)&
!!$             //'/pmd')
!!$        call set_params_NN(ndim,x,rcut,rc3)
!!$      else if( trim(cpot).eq.'linreg') then
!!$        call set_paramsdir_desc(trim(cmaindir)//'/'//trim(cdirname)&
!!$             //'/pmd')
!!$        call set_paramsdir_linreg(trim(cmaindir)//'/'//trim(cdirname)&
!!$             //'/pmd')
!!$        call set_params_linreg(ndim,x)
!!$        nsf = smpl%nsf
!!$        nal = smpl%nal
!!$        nnl = smpl%nnl
!!$        call set_descs(nsf,nal,nnl,smpl%gsf,smpl%dgsf,smpl%igsf)
!!$      else if( trim(cpot).eq.'NN2') then
!!$        call set_paramsdir_desc(trim(cmaindir)//'/'//trim(cdirname)&
!!$             //'/pmd')
!!$        call set_paramsdir_NN2(trim(cmaindir)//'/'//trim(cdirname)&
!!$             //'/pmd')
!!$        call set_params_NN2(ndim,x,nn_nl,nn_nhl)
!!$        nsf = smpl%nsf
!!$        nal = smpl%nal
!!$        nnl = smpl%nnl
!!$        call set_descs(nsf,nal,nnl,smpl%gsf,smpl%dgsf,smpl%igsf)
!!$      else if( index(cpot,'BVS').ne.0 ) then
!!$        call set_paramsdir_Morse(trim(cmaindir)//'/'//trim(cdirname)&
!!$             //'/pmd')
!!$        call set_paramsdir_Coulomb(trim(cmaindir)//'/'//trim(cdirname)&
!!$             //'/pmd')
!!$        if( trim(cpot).eq.'BVS1' ) then
!!$          call set_params_Coulomb(1,x(1),cpot, &
!!$               samples(ismpl)%specorder)
!!$          call set_params_Morse(ndim-1,x(2:ndim),cpot,interact)
!!$        else if( trim(cpot).eq.'BVS' .or. &
!!$             trim(cpot).eq.'BVS2' .or. trim(cpot).eq.'BVS3' ) then
!!$          ndimt = 1+maxisp
!!$          call set_params_Coulomb(ndimt,x(1),cpot, &
!!$               samples(ismpl)%specorder)
!!$          call set_params_Morse(ndim-ndimt,x(ndimt+1:ndim),cpot,interact)
!!$        endif
!!$      endif
!.....Although epot, frcs, and strs are calculated,
!.....only gs is required.
!!$      print *,'ismpl,cpot,ctype,=',ismpl,trim(cpot) &
!!$           ,trim(ctype)
      call run_pmd(smpl,lcalcgrad,ndim,nff,cffs,epot,frcs,strs,rcut &
           ,lfdsgnmat,gwe,gwf,gws)
      samples(ismpl)%gwe(1:ndim)= gwe(1:ndim)
      samples(ismpl)%gwf(1:ndim,1:3,1:natm)= gwf(1:ndim,1:3,1:natm)
      samples(ismpl)%gws(1:ndim,1:6)= gws(1:ndim,1:6)
    enddo  ! ismpl

    if( len(trim(crefstrct)).gt.5 ) then
      if( myid.eq.myidrefsub ) then
        epotsub = samples(isidrefsub)%epot +samples(isidrefsub)%esub
        epotsub = epotsub /samples(isidrefsub)%natm
        gwesub(1:ndim) = samples(isidrefsub)%gwe(1:ndim)
      endif
      call mpi_bcast(epotsub,1,mpi_real8,myidrefsub,mpi_world,ierr)
      call mpi_bcast(gwesub,ndim,mpi_real8,myidrefsub,mpi_world,ierr)
    endif

    gtrnl(1:ndim) = 0d0
    do ismpl=isid0,isid1
      smpl= samples(ismpl)
!.....Since g calc is time consuming,
!.....not calculate g for test set.
      if( smpl%iclass.ne.1 ) cycle
      natm= smpl%natm
      epot= smpl%epot
      swgt= smpl%wgt
!.....Derivative of energy term w.r.t. weights
      if( lematch ) then
        eref= smpl%eref
        esub= smpl%esub
        eerr= smpl%eerr
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
               +tmp/natm/eerr *swgt &
               *(smpl%gwe(1:ndim) -gwesub(1:ndim))
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
               +tmp*smpl%gwe(1:ndim)/natm/eerr *swgt
        endif
      endif
!.....Derivative of force term w.r.t. weights
      if( lfmatch ) then
        frcs(1:3,1:natm)= smpl%fa(1:3,1:natm)
        ferr= smpl%ferr
        ferri= 1d0/ferr
        dn3i= 1d0/3/smpl%nfcal
        do ia=1,natm
          if( smpl%ifcal(ia).eq.0 ) cycle
          do ixyz=1,3
            fdiff(ixyz,ia)= (frcs(ixyz,ia) +smpl%fsub(ixyz,ia) &
                 -(smpl%fref(ixyz,ia))) *ferri
            if( trim(ctype_loss).eq.'LS' ) then
              tmp = 2d0 *fdiff(ixyz,ia)
            else  ! Huber
              if( abs(fdiff(ixyz,ia)).gt.1d0 ) then
                tmp = 2d0 *sign(1d0,fdiff(ixyz,ia))
              else
                tmp = 2d0 *fdiff(ixyz,ia)
              endif
            endif
            gtrnl(1:ndim)= gtrnl(1:ndim) +tmp &
                 *smpl%gwf(1:ndim,ixyz,ia) *dn3i *swgt *ferri
          enddo
        enddo
      endif
!.....Derivative of stress w.r.t. weights
      if( lsmatch ) then
        serr= smpl%serr
        serri= 1d0/serr
        pdiff(1:6) = 0d0
        do ixyz=1,3
          do jxyz=ixyz,3
            k = ivoigt(ixyz,jxyz)
            pdiff(k) = pdiff(k) +( smpl%strs(ixyz,jxyz) &
                 +smpl%ssub(ixyz,jxyz) &
                 -smpl%sref(ixyz,jxyz) ) *serri
          enddo
        enddo
        do k=1,6
          if( trim(ctype_loss).eq.'LS' ) then
            tmp = 2d0 *pdiff(k)
          else  ! Huber
            if( abs(pdiff(k)).gt.1d0 ) then
              tmp = 2d0 *sign(1d0,pdiff(k))
            else
              tmp = 2d0 *pdiff(k)
            endif
          endif
          gtrnl(1:ndim)= gtrnl(1:ndim) +tmp &
               *smpl%gws(1:ndim,k) *swgt *serri /6
        enddo
      endif
    enddo

    tgl= mpi_wtime() -tg0

    tc0= mpi_wtime()
    gtrn(1:ndim) = 0d0
!.....TODO: allreduce may be redundant,  only reducing to node-0 is enough
!           if the minimization routine is wrtten so...
    call mpi_allreduce(gtrnl,gtrn,ndim,mpi_real8,mpi_sum,mpi_world,ierr)
    tcl= tcl +mpi_wtime() -tc0

    gtrn(1:ndim)= gtrn(1:ndim) /swgt2trn
!!$!.....GSCL==1.0 by default, but users can change
!!$    gtrn(:) = gtrn(:) *gscl

!.....only the bottle-neck times are taken into account
    call mpi_reduce(tcl,tcg,1,mpi_real8,mpi_max,0,mpi_world,ierr)
    call mpi_reduce(tgl,tgg,1,mpi_real8,mpi_max,0,mpi_world,ierr)
    tcomm= tcomm +tcg
    tgrad= tgrad +tgg
    return
  end subroutine grad_w_pmd
!=======================================================================
  subroutine ready4pmd(ismpl,ndim,x)
!
!  Get ready for running pmd.
!
    use variables,only: cmaindir,cpot,nsubff,csubffs,mdsys,samples &
         ,maxisp,nn_nl,nn_nhl,nn_sigtype,ctype_loss,rc3,rcut &
         ,interact,interact3,num_interact
    use parallel
    use Coulomb,only: set_paramsdir_Coulomb, set_params_Coulomb
    use Morse,only: set_paramsdir_Morse,set_params_vcMorse,set_params_Morse
    use BMH,only: set_paramsdir_BMH,set_params_BMH
    use Abell,only: set_paramsdir_Abell,set_params_Abell
    use fpc,only: set_paramsdir_fpc,set_params_fpc
    use angular,only: set_paramsdir_angular,set_params_angular
    use EAM,only: set_paramsdir_EAM,set_params_EAM
    use NN,only: set_paramsdir_NN,set_params_NN
    use NN2,only: set_paramsdir_NN2,set_params_NN2,nl,nhl,nlmax,get_NN2_hl1 &
         ,set_NN2_hl1,set_sigtype_NN2
    use linreg,only: set_paramsdir_linreg,set_params_linreg
    use descriptor,only: set_paramsdir_desc,get_descs,get_ints,set_descs &
         ,lupdate_gsf
    implicit none
    integer,intent(in):: ismpl,ndim
    real(8),intent(in):: x(ndim)

    integer:: nsf,nal,nnl,ndimt,ndim0
    character(len=128):: cdirname,ctype
    type(mdsys):: smpl

    logical,external:: string_in_arr

    smpl = samples(ismpl)
    cdirname = smpl%cdirname
    
    if( trim(cpot).eq.'vcMorse' ) then
      call set_paramsdir_Morse(trim(cmaindir)//'/'//trim(cdirname)&
           //'/pmd')
      call set_paramsdir_Coulomb(trim(cmaindir)//'/'//trim(cdirname)&
           //'/pmd')
      call set_params_vcMorse(ndim,x)
    else if( trim(cpot).eq.'Morse' ) then
      call set_paramsdir_Morse(trim(cmaindir)//'/'//trim(cdirname)&
           //'/pmd')
      if( string_in_arr('screened_Coulomb',nsubff,csubffs) .or. &
           string_in_arr('Coulomb',nsubff,csubffs) ) then
        ctype = 'BVS'
      else if( string_in_arr('Ewald_long',nsubff,csubffs) .or.&
           string_in_arr('Ewald',nsubff,csubffs) ) then
        ctype = 'full_Morse'
      else
        ctype = 'full_Morse'
      endif
      call set_params_Morse(ndim,x,ctype,interact)
    else if( trim(cpot).eq.'BMH' ) then
      call set_paramsdir_BMH(trim(cmaindir)//'/'//trim(cdirname)&
           //'/pmd')
      call set_params_BMH(ndim,x,cpot,interact)
    else if( trim(cpot).eq.'Abell' ) then
      call set_paramsdir_Abell(trim(cmaindir)//'/'//trim(cdirname)&
           //'/pmd')
      call set_params_Abell(ndim,x,cpot,interact)
    else if( trim(cpot).eq.'fpc' ) then
      call set_paramsdir_fpc(trim(cmaindir)//'/'//trim(cdirname)&
           //'/pmd')
      call set_paramsdir_Coulomb(trim(cmaindir)//'/'//trim(cdirname)&
           //'/pmd')
      call set_params_Coulomb(1,x(1),cpot, &
           smpl%specorder)
      call set_params_fpc(ndim-1,x(2:ndim),cpot,interact)
    else if( trim(cpot).eq.'EAM' ) then
      call set_paramsdir_EAM(trim(cmaindir)//'/'//trim(cdirname)&
           //'/pmd')
      call set_params_EAM(ndim,x)
    else if( trim(cpot).eq.'NN' ) then
      call set_paramsdir_NN(trim(cmaindir)//'/'//trim(cdirname)&
           //'/pmd')
      call set_params_NN(ndim,x,rcut,rc3)
    else if( trim(cpot).eq.'linreg' ) then
      call set_paramsdir_desc(trim(cmaindir)//'/'//trim(cdirname)&
           //'/pmd')
      call set_paramsdir_linreg(trim(cmaindir)//'/'//trim(cdirname)&
           //'/pmd')
      call set_params_linreg(ndim,x)
      if( .not. lupdate_gsf ) then
        nsf = smpl%nsf
        nal = smpl%nal
        nnl = smpl%nnl
        call set_descs(nsf,nal,nnl,samples(ismpl)%gsf, &
             samples(ismpl)%dgsf,samples(ismpl)%igsf)
      endif
!!$!.....Set lfdsgnmat=.true. to make run_pmd() compute dsgnmat_force related data
!!$      if( l1st .and. lfmatch .and. trim(cfmethod).eq.'dsgnmat' ) then
!!$        lfdsgnmat = .true.
!!$      endif
    else if( trim(cpot).eq.'NN2' ) then
      call set_paramsdir_desc(trim(cmaindir)//'/'//trim(cdirname)&
           //'/pmd')
      call set_paramsdir_NN2(trim(cmaindir)//'/'//trim(cdirname)&
           //'/pmd')
      call set_params_NN2(ndim,x,nn_nl,nn_nhl)
      call set_sigtype_NN2(nn_sigtype)
      if( .not. lupdate_gsf ) then
        nsf = smpl%nsf
        nal = smpl%nal
        nnl = smpl%nnl
        call set_descs(nsf,nal,nnl,samples(ismpl)%gsf, &
             samples(ismpl)%dgsf,samples(ismpl)%igsf)
      endif
    else if( index(cpot,'BVS').ne.0 ) then
      call set_paramsdir_Morse(trim(cmaindir)//'/'//trim(cdirname)&
           //'/pmd')
      call set_paramsdir_Coulomb(trim(cmaindir)//'/'//trim(cdirname)&
           //'/pmd')
      call set_paramsdir_angular(trim(cmaindir)//'/'//trim(cdirname)&
           //'/pmd')
      if( trim(cpot).eq.'BVS1' ) then
        call set_params_Coulomb(1,x(1),cpot,smpl%specorder)
        call set_params_Morse(ndim-1,x(2:ndim),cpot,interact)
      else if( trim(cpot).eq.'BVS' ) then
        ndim0 = 1
        ndimt = 1+maxisp
        call set_params_Coulomb(ndimt,x(ndim0),cpot,smpl%specorder)
        ndim0 = ndim0 +ndimt
        ndimt = num_interact(2)*3
        call set_params_Morse(ndimt,x(ndim0),cpot,interact)
!.....In case of BVSx, not only Coulomb and Morse but also angular should be added.
      else if( trim(cpot).eq.'BVSx' ) then
        ndim0 = 1
        ndimt = 1+maxisp
        call set_params_Coulomb(ndimt,x(ndim0),cpot,smpl%specorder)
        ndim0 = ndim0 +ndimt
        ndimt = num_interact(2)*3
        call set_params_Morse(ndimt,x(ndim0),cpot,interact)
        ndim0 = ndim0 +ndimt
        ndimt = num_interact(3)*3
!!$        if( myid.eq.0 ) print *,'# of interaction = ',num_interact(2) &
!!$             ,num_interact(3)
        call set_params_angular(ndimt,x(ndim0),'angular1',rc3,interact3)
      else
        print *,'ERROR@ready4pmd: No such BVS FF available in fitpot.'
        stop 1
      endif
    endif

  end subroutine ready4pmd
!=======================================================================
  subroutine run_pmd(smpl,lcalcgrad,ndimp,nff,cffs,epot,frcs,strs,rc &
       ,lfdsgnmat,gwe,gwf,gws)
!
!  Run pmd and get energy and forces of the system.
!
    use variables,only: mdsys,maxna,iprint,lematch,lfmatch,lsmatch&
         ,rc_other,maxisp
    use parallel,only: myid_pmd,mpi_comm_pmd,nnode_pmd,myid,mpi_world
    use force
    use descriptor,only: get_dsgnmat_force
    use ZBL,only: r_inner,r_outer
    use pmdio, only: nspmax,has_specorder,specorder,nnmax
    use element
    implicit none
    include "../pmd/params_unit.h"
    type(mdsys),intent(inout):: smpl
    integer,intent(in):: ndimp,nff
    real(8),intent(in):: rc
    real(8),intent(inout):: epot,frcs(3,maxna)
    real(8),intent(out):: strs(3,3)
    logical,intent(in):: lcalcgrad,lfdsgnmat
    character(len=20),intent(in):: cffs(nff)
    real(8),intent(out),optional:: gwe(ndimp),gwf(ndimp,3,maxna),&
         gws(ndimp,6)

    logical,save:: l1st = .true.

    integer:: i,is,maxstp,nerg,npmd,ifpmd,ifdmp,minstp,n_conv,ifsort, &
         nstps_done,ntdst,nx,ny,nz,iprint_pmd,ifcoulomb
    real(8):: am(nspmax),dt,rbuf,dmp,tinit,tfin,ttgt(nspmax),trlx,stgt(3,3),&
         ptgt,srlx,stbeta,strfin,fmv(3,0:9),ptnsr(3,3),ekin,eps_conv &
         ,rc1nn
    logical:: ltdst,lcellfix(3,3),lvc
    character:: ciofmt*6,ctctl*20,cpctl*20,czload_type*5,boundary*3,csp*3
    type(atom):: elem
    logical:: update_force_list

    logical,external:: string_in_arr

    if( l1st ) then
!.....Create MPI COMM for pmd only for the 1st time
      call create_mpi_comm_pmd()
      call init_element()
      l1st = .false.
    endif
!!$    print *,'cdirname,specorder= ',trim(smpl%cdirname) &
!!$         ,(smpl%specorder(i),i=1,2)

    maxstp = 0
    nerg = 1
    npmd = 1
!.....Since at least one of FF requires mass infomation,
!     set mass info from specorder anyways.
    am(:) = 12d0
    do is=1,nspmax
      csp = smpl%specorder(is)
      if( trim(csp).ne.'x' ) then
        elem = get_element(trim(csp))
        am(is) = elem%mass
      endif
    enddo
    dt = 5d0
    ciofmt = 'ascii'
    ifpmd = 0
    rbuf = 0.0d0
    rc1nn = 3.0d0
    ifdmp = 0  ! no damping as well
    dmp = 0.99d0
    minstp = 0
    tinit = 0d0
    tfin = 0d0
    ctctl = 'none'
    ttgt(1:nspmax) = 300d0
    trlx = 100d0
    ltdst = .false.
    ntdst = 1
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
    boundary = 'ppp'
    eps_conv = 1d-3
    ifsort = 1
    lcellfix(1:3,1:3) = .false.
    nx = 1
    ny = 1
    nz = 1
    iprint_pmd = max(0,iprint-10)

    lvc = .false.
    ifcoulomb = 0
    do i=1,nff
      if( trim(cffs(i)).eq.'Ewald_long' .or. &
           trim(cffs(i)).eq.'vcMorse' ) then
        ifcoulomb = 3
        if( .not. smpl%charge_set ) lvc = .true.
!.....Even if lvc is .false., lvc will be set true at the begining of init_force
!.....in case of vcMorse.
      else if( trim(cffs(i)).eq.'screened_Coulomb' ) then
        ifcoulomb = 1
      else if( trim(cffs(i)).eq.'Ewald' ) then
        ifcoulomb = 2
      endif
    enddo

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
!.....Set r_inner/outer arrays in force_ZBL
      r_inner(:) = ri_zbl(:)
      r_outer(:) = ro_zbl(:)
    endif

!.....one_shot force calculation
    call one_shot(smpl%h0,smpl%h,smpl%natm,smpl%tag,smpl%ra &
         ,smpl%va,frcs,smpl%strsi,smpl%eki,smpl%epi &
         ,smpl%chg,smpl%chi &
         ,myid_pmd,mpi_comm_pmd,nnode_pmd,nx,ny,nz &
         ,smpl%specorder,am,dt,rc,rbuf,rc1nn,ptnsr,epot,ekin &
         ,ifcoulomb,lvc,iprint_pmd,lcalcgrad,ndimp,maxisp &
         ,gwe,gwf,gws &
         ,lematch,lfmatch,lsmatch,boundary)
!.....Stress definition, negative as compressive, positive as tensile
    strs(1:3,1:3) = ptnsr(1:3,1:3) *up2gpa*(-1d0)
    if( present(gws) ) gws(1:ndimp,1:6) = gws(1:ndimp,1:6) *up2gpa*(-1d0)
    if( lfdsgnmat ) call get_dsgnmat_force(smpl%dgsfa,mpi_comm_pmd)
!!$  print *,'one_shot done, cdirname,epot = ',trim(smpl%cdirname),epot
!!$  print *,'smpl%natm =',smpl%natm
!!$  write(6,'(a,30es12.4)') 'smpl%epi=',(smpl%epi(i),i=1,smpl%natm)

    if( lvc ) smpl%charge_set = .true.

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
    use Morse,only: set_paramsdir_Morse,set_params_vcMorse,set_params_Morse
    use LJ,only: set_paramsdir_LJ
    use ZBL,only: set_paramsdir_ZBL, r_inner, r_outer
    use Bonny_WRe,only: set_paramsdir_Bonny
    use cspline,only: set_paramsdir_cspline
    use force,only: loverlay
    implicit none

    integer:: i,ismpl,natm
    logical:: lcalcgrad = .false.
    logical:: luse_Morse = .false.
    logical:: luse_Morse_repul = .false.
    logical:: luse_Coulomb = .false.
    logical:: luse_LJ_repul = .false.
    logical:: luse_ZBL = .false.
    logical:: luse_Bonny_WRe = .false.
    logical:: luse_cspline = .false.
    logical:: luse_dipole = .false.
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

      if( .not.allocated(frcs) ) allocate(frcs(3,maxna))

      do i=1,nsubff
        if( index(trim(csubffs(i)),'Morse').ne.0 ) then
          luse_Morse = .true.
        else if( index(trim(csubffs(i)),'Morse_repul').ne.0 ) then
          luse_Morse_repul = .true.
        else if( index(trim(csubffs(i)),'Ewald').ne.0 .or. &
             index(trim(csubffs(i)),'Coulomb').ne.0 .or. &
             index(trim(csubffs(i)),'vcGaussian').ne.0 ) then
          luse_Coulomb = .true.
        else if( index(trim(csubffs(i)),'LJ_repul').ne.0 ) then
          luse_LJ_repul = .true.
        else if( index(trim(csubffs(i)),'ZBL').ne.0 ) then
          luse_ZBL = .true.
        else if( index(trim(csubffs(i)),'Bonny_WRe').ne.0 ) then
          luse_Bonny_WRe = .true.
        else if( index(trim(csubffs(i)),'cspline').ne.0 ) then
          luse_cspline = .true.
        else if( index(trim(csubffs(i)),'dipole').ne.0 ) then
          luse_dipole = .true.
        endif
      enddo

!.....Only at the 1st call, perform pmd to get esubs
      do ismpl=isid0,isid1
        natm = samples(ismpl)%natm
        if( luse_Morse ) then
          call set_paramsdir_Morse(trim(cmaindir)//'/'&
               //trim(samples(ismpl)%cdirname)//'/pmd')
        endif
        if( luse_Morse_repul ) then
          call set_paramsdir_Morse(trim(cmaindir)//'/'&
               //trim(samples(ismpl)%cdirname)//'/pmd')
        endif
        if( luse_Coulomb ) then
          call set_paramsdir_Coulomb(trim(cmaindir)//'/'&
               //trim(samples(ismpl)%cdirname)//'/pmd')
        endif
        if( luse_dipole ) then
          call set_paramsdir_dipole(trim(cmaindir)//'/'&
               //trim(samples(ismpl)%cdirname)//'/pmd')
        endif
        if( luse_LJ_repul ) then
          call set_paramsdir_LJ(trim(cmaindir)//'/'&
               //trim(samples(ismpl)%cdirname)//'/pmd')
        endif
        if( luse_ZBL ) then
          call set_paramsdir_ZBL(trim(cmaindir)//'/'&
               //trim(samples(ismpl)%cdirname)//'/pmd')
        endif
        if( luse_Bonny_WRe ) then
          call set_paramsdir_Bonny(trim(cmaindir)//'/'&
               //trim(samples(ismpl)%cdirname)//'/pmd')
        endif
        if( luse_cspline ) then
          call set_paramsdir_cspline(trim(cmaindir)//'/'&
               //trim(samples(ismpl)%cdirname)//'/pmd')
        endif
        call run_pmd(samples(ismpl),lcalcgrad,nvars,&
             nsubff,csubffs,epot,frcs,strs,rc_other,lfdsgnmat)
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
          if( .not. allocated(ri_zbl) ) then
            allocate(ri_zbl(size(r_inner)),ro_zbl(size(r_outer)))
            ri_zbl(:) = r_inner(:)
            ro_zbl(:) = r_outer(:)
          endif
        endif
      enddo

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
! Compute the mean of input symmetric functions.
!
    use variables
    use parallel
    implicit none 
    real(8),parameter:: tiny = 1d-15
    integer:: nsuml,nsumg,ismpl,natm,isf,jsf,ia,nsf,ixyz
    real(8):: gmeanl,gmean,tmp,gvarl,dgi,dgj,gtmp
    real(8),allocatable:: gsfml(:),gsfvl(:),gsfcl(:,:),gsfvsq(:),gsfsl(:)
    character(len=128):: cnum
!    real(8),allocatable:: dgsft(:,:,)

    nsf = samples(isid0)%nsf
!!$    print *,'myid,isid0,nsf=',myid,isid0,nsf
    if( .not.allocated(gsfml) ) then
      allocate(gsfml(nsf),gsfvl(nsf),gsfcl(nsf,nsf),&
           gsfvsq(nsf),gsfsl(nsf))
    endif
    if( .not. allocated(gsfms) ) then
      allocate(gsfms(nsf),gsfvs(nsf),gsfss(nsf),gsfcorr(nsf,nsf))
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
      write(6,'(a,es12.3)') ' mean of input symmetry functions = ',gsfmean
      write(6,'(a,es12.3)') ' var  of input symmetry functions = ',gsfvar
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

    write(6,'(a)',advance='no') ' Writing out.dsgnmat_atm... '
    write(cnum,'(i0)') nsf
    open(21,file='out.dsgnmat_atm')
    write(21,'(2i10)') nasum,nsf
    do ismpl=isid0,isid1
      natm = samples(ismpl)%natm
      do ia=1,natm
        write(21,'('//trim(cnum)//'es12.3e3)') (samples(ismpl)%gsf(isf,ia),isf=1,nsf)
      enddo
    enddo
    close(21)
    print *,'done'

!.....Design matrix for lasso ()
    if( lematch ) then
      write(6,'(a)',advance='no') ' Writing out.dsgnmat_erg... '
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
      print *,'done'
    endif

!.....Design matrix for force-matching
    if( lfmatch ) then
      write(6,'(a)',advance='no') ' Writing out.dsgnmat_frc... '
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
      print *,'done'
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
      if( myid.eq.0 .and. iprint.ne.0 ) &
           print *,'No normalization of descriptors.'
    else
      if( myid.eq.0 .and. iprint.ne.0 ) &
           print *,'WARNING: no such normalization, '//trim(cnormalize)
    endif
    l1st = .false.
  end subroutine normalize
!=======================================================================
  subroutine normalize_std()
!
!  Normalize inputs (descriptors) wrt standard deviation.
!
    use variables, only: nsmpl,samples,nvars,nalist,vars,vranges&
         ,lnormalized,cpot,gsfvar,gsfvs,gsfms,sgms,sgmis,sgm_min,iprint
    use parallel
    use NN2, only: nl,nhl,mhl
    implicit none
    integer:: ismpl,ia,natm,isf,i,iv,ihl0,ihl1
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
!.....standardize G values
      do ismpl=isid0,isid1
        natm= samples(ismpl)%natm
        if( .not. allocated(samples(ismpl)%gsf) ) then
          print *,'ERROR: gsf not allocated, myid,ismpl,cdirname='&
               ,myid,ismpl,trim(samples(ismpl)%cdirname)
          stop
        endif
        do isf=1,nsf
          samples(ismpl)%gsf(isf,:)= samples(ismpl)%gsf(isf,:) &
               *sgmis(isf)
          samples(ismpl)%dgsf(:,isf,:,:)= &
               samples(ismpl)%dgsf(:,isf,:,:) *sgmis(isf)
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
    else if( trim(cpot).eq.'NN2' ) then
      iv = 0
      do ihl0=1,nhl(0)
        do ihl1=1,mhl(1)
          iv = iv + 1
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
    use variables, only: nsmpl,samples,nvars,nalist,vars,vranges&
         ,lnormalized,cpot,gsfvar,sgms,sgmis,gsfss,sq_min,iprint
    use parallel
    use NN2, only: nl,nhl,mhl
    implicit none
    integer:: ismpl,ia,natm,isf,i,iv,ihl0,ihl1
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
          print *,'ERROR: gsf not allocated, myid,ismpl,cdirname='&
               ,myid,ismpl,trim(samples(ismpl)%cdirname)
          stop
        endif
        do isf=1,nsf
          samples(ismpl)%gsf(isf,:)= samples(ismpl)%gsf(isf,:) &
               *sgmis(isf)
          samples(ismpl)%dgsf(:,isf,:,:)= &
               samples(ismpl)%dgsf(:,isf,:,:) *sgmis(isf)
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
    else if( trim(cpot).eq.'NN2' ) then
      iv = 0
      do ihl0=1,nhl(0)
        do ihl1=1,mhl(1)
          iv = iv + 1
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
    use variables, only: nsmpl,samples,nvars,nalist,vars,vranges&
         ,lnormalized,cnormalize,cpot,gsfvar,sgms,sgmis
    use parallel
    use NN2,only: nhl,mhl
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
      else if( trim(cpot).eq.'NN2' ) then
        iv = 0
        do ihl0=1,nhl(0)
          sgmi = sgmis(ihl0)
          do ihl1=1,mhl(1)
            iv=iv+1
            vars(iv)= vars(iv) *sgmi
            vranges(1:2,iv)= vranges(1:2,iv) *sgmi
          enddo
        enddo
      endif
    endif
    
    lnormalized = .false.
  end subroutine restore_normalize
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
