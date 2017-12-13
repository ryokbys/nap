module fp_common
!-----------------------------------------------------------------------
!                     Last modified: <2017-12-13 23:32:44 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
!
! Module that contains common functions/subroutines for fitpot.
!
  implicit none
  save
  real(8),allocatable:: fdiff(:,:),frcs(:,:),gtrnl(:)
  real(8):: pdiff(6), ptnsr(3,3)
  real(8):: epotsub
  
  logical:: fp_common_initialized= .false.

  integer,parameter:: ivoigt(3,3)= &
       reshape((/ 1, 6, 5, 6, 2, 4, 5, 4, 3 /),shape(ivoigt))
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
      write(6,'(a)') ' Weights to divide evaluation value:'
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
         ,crefstrct,erefsub,myidrefsub,isidrefsub,iprint
    use parallel
    use minimize
    use Coulomb,only: set_paramsdir_Coulomb
    use Morse,only: set_paramsdir_Morse,set_params_vcMorse,set_params_Morse
    use EAM,only: set_paramsdir_EAM,set_params_EAM
    use NN,only: set_paramsdir_NN,set_params_NN
    implicit none
    integer,intent(in):: ndim
    real(8),intent(in):: x(ndim)
    real(8),intent(out):: ftrn,ftst

    integer:: ismpl,natm,ia,ixyz,jxyz,idim,k
    real(8):: dn3i,ediff,fscale,eref,epot,swgt,wgtidv,esub
    real(8):: eerr,ferr,ferri,serr,serri,strs(3,3)
    real(8):: ftrnl,ftstl,ftmp
    real(8):: edenom,fdenom
    real(8):: tfl,tcl,tfg,tcg,tf0,tc0
    type(mdsys):: smpl
    logical,save:: l1st = .true.
    logical,parameter:: lcalcgrad = .false.
    character(len=128):: cdirname,ctype

    logical,external:: string_in_arr

    nfunc= nfunc +1

!!$    print *,'func_w_pmd: 01'
    tc0= mpi_wtime()
    call mpi_bcast(x,ndim,mpi_real8,0,mpi_world,ierr)
    tcl= mpi_wtime() -tc0
    tf0= mpi_wtime()

    if( l1st ) then
      if( .not.fp_common_initialized ) call init()
      if( .not.allocated(fdiff) ) allocate(fdiff(3,maxna),frcs(3,maxna))
    endif

!!$    print *,'func_w_pmd: 02'
    if( .not. lematch .and. .not.lfmatch .and. .not.lsmatch ) then
      if( myid.eq.0 ) then
        print *,'Nothing to be fitted.'
      endif
      call mpi_finalize(ierr)
      stop
    endif

!!$    print *,'myid,isid0,isid1=',myid,isid0,isid1
!!$    print *,'func_w_pmd: 03'
    ftrnl = 0d0
    ftstl = 0d0
    do ismpl=isid0,isid1
!!$      print *,'func_w_pmd: 04,ismpl,size(samples)=',ismpl,size(samples)
      smpl = samples(ismpl)
      natm= smpl%natm
      cdirname = trim(smpl%cdirname)
      if( trim(cpot).eq.'vcMorse' ) then
        call set_paramsdir_Morse(trim(cmaindir)//'/'//trim(cdirname)&
             //'/pmd')
        call set_paramsdir_Coulomb(trim(cmaindir)//'/'//trim(cdirname)&
             //'/pmd')
        call set_params_vcMorse(ndim,x)
      else if( trim(cpot).eq.'Morse' ) then
        call set_paramsdir_Morse(trim(cmaindir)//'/'//trim(cdirname)&
             //'/pmd')
        if( string_in_arr('screened_Coulomb',nsubff,csubffs) ) then
          ctype = 'BVS'
        else if( string_in_arr('long_Coulomb',nsubff,csubffs) ) then
          ctype = 'full_Morse'
        else
          ctype = 'full_Morse'
        endif
        call set_params_Morse(ndim,x,ctype)
      else if( trim(cpot).eq.'EAM' ) then
        call set_paramsdir_EAM(trim(cmaindir)//'/'//trim(cdirname)&
             //'/pmd')
        call set_params_EAM(ndim,x)
      else if( trim(cpot).eq.'NN' ) then
        call set_paramsdir_NN(trim(cmaindir)//'/'//trim(cdirname)&
             //'/pmd')
        call set_params_NN(ndim,x,rcut,rc3)
      endif
!!$      print *,'func_w_pmd: 04-1, ismpl,lcalcgrad,ndim,nff='&
!!$           ,ismpl,lcalcgrad,ndim,nff
      call run_pmd(smpl,lcalcgrad,ndim,nff,cffs,epot,frcs,strs)
      samples(ismpl)%epot = epot
      samples(ismpl)%fa(1:3,1:natm) = frcs(1:3,1:natm)
      samples(ismpl)%strs(1:3,1:3) = strs(1:3,1:3)
      if( iprint.ge.10 ) then
        write(6,'(a,2i4,1x,a,7es12.4)') ' myid,ismpl,cdirname,epot,strs= ', &
             myid,ismpl,trim(cdirname), &
             epot,strs(1,1),strs(2,2),strs(3,3), &
             strs(2,3),strs(1,3),strs(1,2)
      endif
    enddo

!!$    print *,'func_w_pmd: 05'
    if( len(trim(crefstrct)).gt.5 ) then
      if( myid.eq.myidrefsub ) then
        epotsub = samples(isidrefsub)%epot +samples(isidrefsub)%esub
        epotsub = epotsub /samples(isidrefsub)%natm
      endif
      call mpi_bcast(epotsub,1,mpi_real8,myidrefsub,mpi_world,ierr)
!!$      print *,'myid,epotsub=', myid,epotsub
    endif

    do ismpl=isid0,isid1
!!$      print *,'func_w_pmd: 06,ismpl=',ismpl
      smpl = samples(ismpl)
      natm = smpl%natm
      epot = smpl%epot
      ftmp = 0d0
      swgt = smpl%wgt
      if( lematch ) then
        eref= smpl%eref
        esub= smpl%esub
        eerr = smpl%eerr
!.....Energy matching
        if( len(trim(crefstrct)).gt.5 ) then
          ediff= (epot-epotsub*natm+esub -(eref-erefsub*natm))/natm /eerr
        else
          ediff= (epot+esub -eref)/natm /eerr
        endif
        ediff= ediff*ediff
        ftmp= ftmp +ediff *swgt
!!$        write(6,'(a,2i3,a20,6f12.5)') &
!!$             'myid,ismpl,cdirname,(eref-erefsub),epot,ediff,esub,epotsub,erefsub=',&
!!$             myid,ismpl,trim(samples(ismpl)%cdirname),&
!!$             (eref-erefsub*natm)/natm,epot/natm,&
!!$             abs(epot-epotsub*natm+esub -(eref-erefsub*natm))/natm,esub/natm,&
!!$             epotsub, erefsub
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
            fdiff(ixyz,ia)= fdiff(ixyz,ia)*fdiff(ixyz,ia)
            ftmp= ftmp +fdiff(ixyz,ia) *dn3i *swgt
          enddo
        enddo
!!$        write(6,'(a,i6,es12.4)') 'ismpl,ftmp=',ismpl,ftmp
      endif

!.....Stress matching
      if( lsmatch ) then
!!$!.....Make current ptnsr of the system
!!$        smpl%cptnsr(1:3,1:3)= 0d0
!!$        do ia=1,natm
!!$          do jxyz=1,3
!!$            do ixyz=1,3
!!$              smpl%cptnsr(ixyz,jxyz)=smpl%cptnsr(ixyz,jxyz) &
!!$                   +smpl%strsi(ixyz,jxyz,ia)
!!$            enddo
!!$          enddo
!!$        enddo

!!$!  vol = 1d0
!!$        smpl%cptnsr(1:3,1:3) = smpl%cptnsr(1:3,1:3) !/vol

!.....Compare these ptnsr elements with sref elements
        serr = smpl%serr
        serri = 1d0/serr
        do ixyz=1,3
          do jxyz=ixyz,3
            k = ivoigt(ixyz,jxyz)
            pdiff(k)= (smpl%strs(ixyz,jxyz) +smpl%ssub(ixyz,jxyz) &
                 -smpl%sref(ixyz,jxyz)) *serri
            pdiff(k)= pdiff(k)*pdiff(k)/6
            ftmp= ftmp +pdiff(k) *swgt
          enddo
        enddo
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
         ,rcut,rc3,myidrefsub,isidrefsub,iprint
    use parallel
    use minimize
    use Coulomb,only: set_paramsdir_Coulomb
    use Morse,only: set_paramsdir_Morse,set_params_vcMorse,set_params_Morse
    use EAM,only: set_paramsdir_EAM,set_params_EAM
    use NN,only: set_paramsdir_NN,set_params_NN
    implicit none
    integer,intent(in):: ndim
    real(8),intent(in):: x(ndim)
    real(8),intent(out):: gtrn(ndim)
    
    integer:: ismpl,i,idim,natm,k,ia,ixyz,jxyz
    real(8):: tcl,tgl,tcg,tgg,tc0,tg0,epot,esub,strs(3,3),dn3i
    real(8):: ediff,eerr,eref,swgt,ferr,ferri,serr,serri
    type(mdsys):: smpl
    logical,parameter:: lcalcgrad = .true. 
    character(len=128):: cdirname,ctype

    logical,external:: string_in_arr

!!$    print *,'grad_w_pmd'
    if( .not.allocated(gtrnl) ) allocate(gtrnl(ndim))

    ngrad= ngrad +1
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
!.....Since g calc is time consuming,
!.....not calculate g for test set.
      if( smpl%iclass.ne.1 ) cycle
      if( trim(cpot).eq.'vcMorse' ) then
        call set_paramsdir_Morse(trim(cmaindir)//'/'//trim(cdirname)&
             //'/pmd')
        call set_paramsdir_Coulomb(trim(cmaindir)//'/'//trim(cdirname)&
             //'/pmd')
        call set_params_vcMorse(ndim,x)
      else if( trim(cpot).eq.'Morse' ) then
        call set_paramsdir_Morse(trim(cmaindir)//'/'//trim(cdirname)&
             //'/pmd')
        if( string_in_arr('screened_Coulomb',nsubff,csubffs) ) then
          ctype = 'BVS'
        else if( string_in_arr('long_Coulomb',nsubff,csubffs) ) then
          ctype = 'full_Morse'
        else
          ctype = 'full_Morse'
        endif
        call set_params_Morse(ndim,x,ctype)
      else if( trim(cpot).eq.'EAM' ) then
        call set_paramsdir_EAM(trim(cmaindir)//'/'//trim(cdirname)&
             //'/pmd')
        call set_params_EAM(ndim,x)
      else if( trim(cpot).eq.'NN' ) then
        call set_paramsdir_NN(trim(cmaindir)//'/'//trim(cdirname)&
             //'/pmd')
        call set_params_NN(ndim,x,rcut,rc3)
      endif
!.....Although epot, frcs, and strs are calculated,
!.....only gs is required.
!!$      print *,'ismpl,cpot,ctype,=',ismpl,trim(cpot) &
!!$           ,trim(ctype)
      call run_pmd(smpl,lcalcgrad,ndim,nff,cffs,epot,frcs,strs)
    enddo  ! ismpl

    if( len(trim(crefstrct)).gt.5 ) then
!!$      print *,'Going to get epotsub...'
      if( myid.eq.myidrefsub ) then
        epotsub = samples(isidrefsub)%epot +samples(isidrefsub)%esub
        epotsub = epotsub /samples(isidrefsub)%natm
      endif
      call mpi_bcast(epotsub,1,mpi_real8,myidrefsub,mpi_world,ierr)
    endif

!!$    print *,'Going to compute derivative of loss function'
    gtrnl(1:ndim) = 0d0
    do ismpl=isid0,isid1
!!$      print *,'isid0,isid1,ismpl,size(samples)=',isid0,isid1,ismpl,size(samples)
      smpl= samples(ismpl)
      natm= smpl%natm
      epot= smpl%epot
      swgt= smpl%wgt
!!$      print *,'ismpl,natm,epot=',ismpl,natm,epot
!.....Derivative of energy term w.r.t. weights
      if( lematch ) then
        eref= smpl%eref
        esub= smpl%esub
        eerr= smpl%eerr
        if( len(trim(crefstrct)).gt.5 ) then
          ediff= (epot-epotsub*natm+esub -(eref-erefsub*natm))/natm /eerr
        else
          ediff= (epot+esub -eref)/natm /eerr
        endif
        gtrnl(1:ndim) = gtrnl(1:ndim) &
             +2d0*ediff*smpl%gwe(1:ndim)/natm/eerr *swgt
!!$        print *,'ismpl,ediff,gwe=',ismpl,ediff,smpl%gwe(1:ndim)
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
            gtrnl(1:ndim)= gtrnl(1:ndim) +2d0*fdiff(ixyz,ia) &
                 *smpl%gwf(1:ndim,ixyz,ia) *dn3i *swgt *ferri
          enddo
        enddo
!!$        print *,'ismpl,ediff,gwf=',ismpl,ediff,smpl%gwf(1:ndim,1,1)
      endif
!.....Derivative of stress w.r.t. weights
      if( lsmatch ) then
        serr= smpl%serr
        serri= 1d0/serr
        do ixyz=1,3
          do jxyz=ixyz,3
            k = ivoigt(ixyz,jxyz)
            pdiff(k) = (smpl%strs(ixyz,jxyz) +smpl%ssub(ixyz,jxyz) &
                 -smpl%sref(ixyz,jxyz)) *serri
            gtrnl(1:ndim)= gtrnl(1:ndim) +2d0*pdiff(k) &
                 *smpl%gws(1:ndim,k)/6 *serri
          enddo
        enddo
!!$        print *,'ismpl,ediff,gws=',ismpl,ediff,smpl%gws(1:ndim,1)
      endif
    enddo

    tgl= mpi_wtime() -tg0

    tc0= mpi_wtime()
    gtrn(1:ndim) = 0d0
    call mpi_allreduce(gtrnl,gtrn,ndim,mpi_real8,mpi_sum,mpi_world,ierr)
    tcl= mpi_wtime() -tc0

    gtrn(1:ndim)= gtrn(1:ndim) /swgt2trn

!.....only the bottle-neck times are taken into account
    call mpi_reduce(tcl,tcg,1,mpi_real8,mpi_max,0,mpi_world,ierr)
    call mpi_reduce(tgl,tgg,1,mpi_real8,mpi_max,0,mpi_world,ierr)
    tcomm= tcomm +tcg
    tgrad= tgrad +tgg
    
    return
  end subroutine grad_w_pmd
!=======================================================================
  function ndat_in_line(ionum,delim) result(ndat)
    implicit none
    integer,intent(in):: ionum
    character(len=1),intent(in):: delim
    integer:: ndat
    integer,external:: num_data
    character:: ctmp*128

    read(ionum,'(a)') ctmp
    ndat = num_data(trim(ctmp),delim)
    backspace(ionum)
    return

  end function ndat_in_line
  
end module fp_common
!-----------------------------------------------------------------------
! Local Variables:
! compile-command: "make fitpot"
! End:
