module fp_common
!-----------------------------------------------------------------------
!                     Last modified: <2018-02-14 16:33:48 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
!
! Module that contains common functions/subroutines for fitpot.
!
  implicit none
  save
  real(8),allocatable:: fdiff(:,:),frcs(:,:),gtrnl(:)
  real(8),allocatable:: gwe(:),gwf(:,:,:),gws(:,:)
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
        else if( string_in_arr('Ewald_long',nsubff,csubffs) .or.&
             string_in_arr('Ewald',nsubff,csubffs) ) then
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
        ediff= ediff*ediff
        ftmp= ftmp +ediff *swgt
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
        do k=1,6
          pdiff(k)= pdiff(k)*pdiff(k)
          ftmp= ftmp +pdiff(k) *swgt /6
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
    if( .not.allocated(gwe) ) allocate(gwe(ndim),gwf(ndim,3,maxna)&
         ,gws(ndim,6))

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
        else if( string_in_arr('Ewald_long',nsubff,csubffs) .or.&
             string_in_arr('Ewald',nsubff,csubffs) ) then
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
      call run_pmd(smpl,lcalcgrad,ndim,nff,cffs,epot,frcs,strs&
           ,gwe,gwf,gws)
      samples(ismpl)%gwe(1:ndim)= gwe(1:ndim)
      samples(ismpl)%gwf(1:ndim,1:3,1:natm)= gwf(1:ndim,1:3,1:natm)
      samples(ismpl)%gws(1:ndim,1:6)= gws(1:ndim,1:6)
!!$      print '(a,i3,9es11.3)',' ismpl,gwe=',ismpl,smpl%gwe(1:ndim)
    enddo  ! ismpl

    if( len(trim(crefstrct)).gt.5 ) then
      if( myid.eq.myidrefsub ) then
        epotsub = samples(isidrefsub)%epot +samples(isidrefsub)%esub
        epotsub = epotsub /samples(isidrefsub)%natm
      endif
      call mpi_bcast(epotsub,1,mpi_real8,myidrefsub,mpi_world,ierr)
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
!!$        print *,'ismpl,epot,esub,eref,ediff,natm,eerr,swgt=' &
!!$             ,ismpl,epot,esub,eref,ediff,natm,eerr,swgt
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
          gtrnl(1:ndim)= gtrnl(1:ndim) +2d0*pdiff(k) &
               *smpl%gws(1:ndim,k)/6 *swgt *serri
        enddo
!!$        print *,'ismpl,ediff,gws=',ismpl,ediff,smpl%gws(1:ndim,1)
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

!.....only the bottle-neck times are taken into account
    call mpi_reduce(tcl,tcg,1,mpi_real8,mpi_max,0,mpi_world,ierr)
    call mpi_reduce(tgl,tgg,1,mpi_real8,mpi_max,0,mpi_world,ierr)
    tcomm= tcomm +tcg
    tgrad= tgrad +tgg

    return
  end subroutine grad_w_pmd
!=======================================================================
  subroutine run_pmd(smpl,lcalcgrad,ndimp,nff,cffs,epot,frcs,strs,&
       gwe,gwf,gws)
!
!  Run pmd and get energy and forces of the system.
!
    use variables,only: rcut,mdsys,maxna,iprint,lematch,lfmatch,lsmatch
    use parallel,only: myid_pmd,mpi_comm_pmd,nnode_pmd,myid,mpi_world
    use force
    implicit none
    include "../pmd/params_unit.h"
    type(mdsys),intent(inout):: smpl
    integer,intent(in):: ndimp,nff
    real(8),intent(inout):: epot,frcs(3,maxna)
    real(8),intent(out):: strs(3,3)
    logical,intent(in):: lcalcgrad
    character(len=20),intent(in):: cffs(nff)
    real(8),intent(out),optional:: gwe(ndimp),gwf(ndimp,3,maxna),&
         gws(ndimp,6)

    logical,save:: l1st = .true.

    integer:: i,maxstp,nerg,npmd,ifpmd,ifdmp,minstp,n_conv,ifsort, &
         nismax,nstps_done,ntdst,nx,ny,nz,iprint_pmd,ifcoulomb
    real(8):: am(9),dt,rc,rbuf,dmp,tinit,tfin,ttgt(9),trlx,stgt(3,3),&
         ptgt,srlx,stbeta,strfin,fmv(3,0:9),ptnsr(3,3),ekin,eps_conv
    logical:: ltdst,lcellfix(3,3),lvc
    character:: ciofmt*6,ctctl*20,cpctl*20,czload_type*5,boundary*3
    logical:: update_force_list

    logical,external:: string_in_arr

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
    if( .not.allocated(force_list) ) then
      update_force_list = .true.
    else if( nff.ne.num_forces ) then
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
      if( allocated(force_list) ) deallocate(force_list)
      num_forces = nff
      allocate(force_list(num_forces))
      do i=1,num_forces
        force_list(i) = trim(cffs(i))
      enddo
    endif

!.....one_shot force calculation
    call one_shot(smpl%h0,smpl%h,smpl%natm,smpl%tag,smpl%ra &
         ,smpl%va,frcs,smpl%strsi,smpl%eki,smpl%epi &
         ,smpl%chg,smpl%chi &
         ,myid_pmd,mpi_comm_pmd,nnode_pmd,nx,ny,nz &
         ,nismax,am,dt,rc,rbuf,ptnsr,epot,ekin &
         ,ifcoulomb,lvc,iprint_pmd,lcalcgrad,ndimp &
         ,gwe,gwf,gws &
         ,lematch,lfmatch,lsmatch,boundary)
    strs(1:3,1:3) = ptnsr(1:3,1:3) *up2gpa*(-1d0)
    if( present(gws) ) gws(1:ndimp,1:6) = gws(1:ndimp,1:6) *up2gpa*(-1d0)
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
        call run_pmd(samples(ismpl),lcalcgrad,nvars,&
             nsubff,csubffs,epot,frcs,strs)
!!$      print *,'myid,ismpl,epot=',myid,ismpl,epot
        samples(ismpl)%esub = epot
        samples(ismpl)%fsub(1:3,1:natm) = frcs(1:3,1:natm)
        samples(ismpl)%ssub(1:3,1:3) = strs(1:3,1:3)
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
