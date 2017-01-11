module NN
!-----------------------------------------------------------------------
!                     Last modified: <2017-01-11 15:04:57 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
!.....parameter file name
  save
  character(128),parameter:: cpfname= 'in.params.NN'
  character(128),parameter:: ccfname='in.const.NN'
  character(128),parameter:: cmbfname='in.comb.NN'
!.....NN mode, which determines if it contains bias node
  integer:: mode = -1
  integer,parameter:: maxnl= 2
  integer:: nl,nsp,nsf2,nsf3,ncmb2,ncmb3
!.....number of nodes in each layer
!.....  nhl includes bias nodes whereas mhl does not
  integer:: nhl(0:maxnl+1),mhl(0:maxnl+1)
  integer,allocatable:: nwgt(:)
  real(8),allocatable:: wgt11(:,:),wgt12(:)
  real(8),allocatable:: wgt21(:,:),wgt22(:,:),wgt23(:)
  real(8):: gsfmean,gsfvar

  type smpldata
    real(8),allocatable:: gsf(:,:),dgsf(:,:,:,:)
    real(8),allocatable:: hl1(:,:),hl2(:,:)
    real(8),allocatable:: gsfo(:,:)
  end type smpldata

  type(smpldata),allocatable:: sds(:)

  integer:: maxna
  real(8),allocatable:: fdiff(:,:)
  real(8),allocatable:: gmax(:),gmin(:)

  logical:: lstandard= .false.

  integer:: nterm_trn = 0
  integer:: nterm_tst = 0

contains
!=======================================================================
  subroutine NN_init()
    use variables
    use parallel
    use minimize
    implicit none 
    integer:: itmp,i,nw,natm,ismpl,ihl0,ihl1,itmp2,ndat
    real(8):: swgt,dtmp
    character:: ctmp*128
    integer,external:: ndat_in_line

    tfunc= 0d0
    tgrad= 0d0
    nfunc= 0
    ngrad= 0

    !.....read in.const.NN to get nl,nsp,nhl(:)
    if( myid.eq.0 ) then
      open(20,file=trim(cmaindir)//'/'//trim(ccfname),status='old')
      ndat = ndat_in_line(20,' ')
      if( ndat.eq.4 ) then  ! old in.const.NN file
        ! set mode = 1 and reread 1st line without reading mode
        mode = 1
        read(20,*) nl,nsp,nhl(0:nl)
      else if( ndat.eq.5 ) then
        read(20,*) nl,nsp,nhl(0:nl),mode
      endif
!!$      print *,'nl,nsp,nhl(0:nl),mode = ',nl,nsp,nhl(0:nl),mode
      nsf2= 0
      nsf3= 0
      do while(.true.)
        read(20,*,end=10) itmp !,itmp2,itmp2,dtmp,dtmp
        if( itmp.le.100 ) then
          nsf2=nsf2+1
        else if( itmp.le.200 ) then
          nsf3=nsf3+1
        endif
!!$        print *, nsf2,nsf3,dtmp
      enddo
10    close(20)
      nhl(nl+1)= 1
    endif
    call mpi_bcast(nl,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(nsp,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(nhl,nl+2,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(nsf2,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(nsf3,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(mode,1,mpi_integer,0,mpi_world,ierr)


!!$    ncmb2= nsp +factorial(nsp,2)/factorial((nsp-2),2)/2
!!$    ncmb3= nsp*ncmb2
!.....check number of symmetry functions
    if( nhl(0).ne.nsf2 +nsf3 ) then
      if( myid.eq.0) then
        print *,'[Error] nhl(0).ne.nsf2 +nsf3 '
!!$        print *,' ncmb2,ncmb3,nsf2,nsf3=',ncmb2,ncmb3,nsf2,nsf3
        print *,'   nsf2,nsf3 = ',nsf2,nsf3
        print *,'   nhl(0),nsf2+nsf3=',nhl(0),nsf2+nsf3
      endif
      call mpi_finalize(ierr)
      stop
    endif

!.....correct number of nodes according to mode value
    mhl(0:nl+1) = nhl(0:nl+1)
    if( mode.ge.10 ) then
      nhl(0) = nhl(0) +1
      nhl(1) = nhl(1) +1
      if( nl.eq.2 ) then
        nhl(2) = nhl(2) +1
      endif
    endif
    if( mode.eq.12 ) nhl(0) = nhl(0) +1
    if( myid.eq.0 ) then
      write(6,'(a,5i5)') ' nhl(0:nl+1)=',nhl(0:nl+1)
      write(6,'(a,5i5)') ' mhl(0:nl+1)=',mhl(0:nl+1)
    endif

!.....total number of weights
    allocate(nwgt(nl+1))
    nw= 0
    do i=1,nl+1
      nwgt(i)= nhl(i-1)*mhl(i)
      nw= nw +nwgt(i)
    enddo
    if(myid.eq.0) then
!!$      print *,'nsf2,nsf3,ncmb2,ncmb3=',nsf2,nsf3,ncmb2,ncmb3
      print *,'nsf2,nsf3 = ',nsf2,nsf3
      print *,'nhl(0:nl+1),nw = ',nhl(0:nl+1),nw
    endif

!.....training set
    allocate(sds(isid0:isid1))
    do ismpl=isid0,isid1
      natm= samples(ismpl)%natm
      allocate( sds(ismpl)%gsf(natm,nhl(0)) )
      if( lfmatch ) allocate( sds(ismpl)%dgsf(3,natm,natm,nhl(0)) )
      if( nl.eq.1 ) then
        allocate(sds(ismpl)%hl1(natm,nhl(1)))
      else if( nl.eq.2 ) then
        allocate(sds(ismpl)%hl1(natm,nhl(1)),&
             sds(ismpl)%hl2(natm,nhl(2)))
      endif
    enddo
    call get_bases()

    maxna= 0
    do ismpl=isid0,isid1
      if( maxna.lt.samples(ismpl)%natm )  &
           maxna= samples(ismpl)%natm
    enddo
!!$    print *,' myid, max num of atoms [maxna] =',myid,maxna
    allocate(fdiff(3,maxna))

    gsfmean= get_mean_input()
    gsfvar = get_variance_input(gsfmean)
    if(myid.eq.0) then
      write(6,'(a,es12.3)') ' mean of input     = ',gsfmean
      write(6,'(a,es12.3)') ' variance of input = ',gsfvar
    endif

    if( (cpena.eq.'glasso' .or. cpena.eq.'lasso' .or. &
         cfmethod.eq.'gfs') .and. &
      (cnormalize(1:3).ne.'var' .and. cnormalize(1:3).ne.'max') ) then
      if(myid.eq.0) then
        print *,'Error: no normalization with glasso, lasso, or gfs'
        print *,'   might cause no-good optimization.'
        print *,'   You should use var or max for normalize_input.'
      endif
      call mpi_finalize(ierr)
      stop
    endif
    call NN_standardize()

!.....make groups for group lasso
    if( trim(cpena).eq.'glasso' &
         .or. trim(cfmethod).eq.'gfs') then
      ngl= nhl(0)
      allocate(iglid(nw),glval(0:ngl))
      iglid(1:nw)= 0
      i= 0
      do ihl0=1,nhl(0)
        do ihl1=1,nhl(1)
          i= i +1
          iglid(i)= ihl0
        enddo
      enddo
!.....weights between hidden nodes are not penalized in g-lasso
      do i=nhl(0)*nhl(1)+1,nw
        iglid(i)= -1
      enddo
    endif

!.....set nominator for sample weights
    swgt = 0d0
    do ismpl=isid0,isid1
      swgt = swgt +samples(ismpl)%wgt
    enddo
    swgt2 = 0d0
    call mpi_allreduce(swgt,swgt2,1,mpi_double_precision,mpi_sum &
         ,mpi_world,ierr)
    swgt2 = swgt2*2d0
    if( myid.eq.0 ) write(6,'(a,es12.3)') ' swgt2 = ',swgt2

    if( myid.eq.0 ) print *, 'NN_init done.'
  end subroutine NN_init
!=======================================================================
  function NN_func(ndim,x)
    use variables,only:nsmpl,nsmpl_trn,samples,nprcs,tfunc &
         ,lfmatch,nfunc,tcomm,mdsys,erefmin &
         ,cmaindir,epse,epsf,cevaltype,swgt2
    use parallel
    use minimize
    implicit none
    integer,intent(in):: ndim
    real(8),intent(in):: x(ndim)
    real(8):: NN_func

    integer:: ismpl,natm,ia,ixyz,idim
    real(8):: dn3i,ediff,fscale,eref,swgt,wgtidv
    real(8):: eerr,ferr,ferri
    real(8):: flocal
    real(8):: edenom,fdenom
    real(8):: tfl,tcl,tfg,tcg,tf0,tc0
    type(mdsys):: smpl
    logical:: l1st = .true.

    nfunc= nfunc +1

    tc0= mpi_wtime()
    call mpi_bcast(x,ndim,mpi_double_precision,0,mpi_world,ierr)
    tcl= mpi_wtime() -tc0
    tf0= mpi_wtime()
    call vars2wgts(ndim,x)

    if( l1st ) then
      call count_nterms()
      if( myid.eq.0 ) then
        write(6,*) ' nterm_trn, nterm_tst = ',nterm_trn, nterm_tst
      endif
    endif
    
    do ismpl=isid0,isid1
      if( nl.eq.1 ) then
        call calc_ef1(samples(ismpl),sds(ismpl))
      else if( nl.eq.2 ) then
        call calc_ef2(samples(ismpl),sds(ismpl))
      endif
    enddo

!!$    NN_func= 0d0
    flocal= 0d0
    do ismpl=isid0,isid1
      smpl= samples(ismpl)
      if( smpl%iclass.ne.1 ) cycle
      natm= smpl%natm
      eref= smpl%eref
      eerr = smpl%eerr
      swgt = smpl%wgt
      ediff= (smpl%epot -eref)/natm /eerr
      ediff= ediff*ediff
      flocal= flocal +ediff *swgt
      if( .not. lfmatch ) cycle
      if( smpl%nfcal.eq.0 ) cycle
      ferr = smpl%ferr
      ferri = 1d0/ferr
      dn3i = 1d0/3/smpl%nfcal
      do ia=1,natm
        if( smpl%ifcal(ia).eq.0 ) cycle
        do ixyz=1,3
          fdiff(ixyz,ia)= (smpl%fa(ixyz,ia) &
               -smpl%fref(ixyz,ia)) *ferri
          fdiff(ixyz,ia)= fdiff(ixyz,ia)*fdiff(ixyz,ia)
          flocal= flocal +fdiff(ixyz,ia) *dn3i *swgt
        enddo
      enddo
    enddo

!    tfunc= tfunc +mpi_wtime() -tf0
    tfl = mpi_wtime() -tf0

    tc0= mpi_wtime()
    NN_func= 0d0
    call mpi_allreduce(flocal,NN_func,1,mpi_double_precision &
         ,mpi_sum,mpi_world,ierr)
    tcl = tcl + (mpi_wtime() -tc0)
!    tcomm= tcomm +mpi_wtime() -tc0
!    NN_func= NN_func/nsmpl_trn
    NN_func= NN_func /swgt2

!.....only the bottle-neck times are taken into account
    call mpi_reduce(tcl,tcg,1,mpi_double_precision,mpi_max,0 &
         ,mpi_world,ierr)
    call mpi_reduce(tfl,tfg,1,mpi_double_precision,mpi_max,0 &
         ,mpi_world,ierr)
    tcomm= tcomm +tcg
    tfunc= tfunc +tfg
    l1st = .false.

  end function NN_func
!=======================================================================
  function NN_fs(ndim,x)
    use variables
    use parallel
    use minimize
    implicit none
    integer,intent(in):: ndim
    real(8),intent(in):: x(ndim)
    real(8):: NN_fs
    integer:: natm,ia,ixyz,idim,i
    real(8):: dn3i,ediff,fscale,eref,swgt,flocal,wgtidv
    real(8):: tfl,tcl,tfg,tcg,tf0,tc0
    real(8):: eerr,ferr,ferri
    integer:: ismpl
    type(mdsys):: smpl
    logical,save:: l1st = .true.

    nfunc=nfunc +1

    tc0= mpi_wtime()
    call mpi_bcast(x,ndim,mpi_double_precision,0,mpi_world,ierr)
    tcl= mpi_wtime() -tc0
    tf0= mpi_wtime()
    call vars2wgts(ndim,x)

    if( l1st ) then
      call count_nterms()
      if( myid.eq.0 ) then
        write(6,*) ' nterm_trn, nterm_tst = ',nterm_trn, nterm_tst
      endif
    endif

    do i=1,nsgdbsize
      ismpl= ismplsgd(i)
      if( nl.eq.1 ) then
        call calc_ef1(samples(ismpl),sds(ismpl))
      else if( nl.eq.2 ) then
        call calc_ef2(samples(ismpl),sds(ismpl))
      endif
    enddo

    flocal= 0d0
    do i=1,nsgdbsize
      ismpl= ismplsgd(i)
      smpl= samples(ismpl)
      if( smpl%iclass.ne.1 ) cycle
      natm= smpl%natm
      eref= smpl%eref
      eerr = smpl%eerr
      swgt = smpl%wgt
      ediff= (smpl%epot -eref)/natm /eerr
      ediff= ediff*ediff
      flocal= flocal +ediff*swgt
      if( .not. lfmatch ) cycle
      if( smpl%nfcal.eq.0 ) cycle
      ferr = smpl%ferr
      ferri = 1d0/ferr
      dn3i = 1d0/3/smpl%nfcal
      do ia=1,natm
        if( smpl%ifcal(ia).eq.0 ) cycle
        do ixyz=1,3
          fdiff(1:3,ia)= (smpl%fa(1:3,ia) &
               -smpl%fref(1:3,ia)) *ferri
          fdiff(1:3,ia)= fdiff(1:3,ia)*fdiff(1:3,ia)
          flocal= flocal +fdiff(ixyz,ia)*dn3i *swgt
        enddo
      enddo
    enddo

    tfl = mpi_wtime() -tf0

    tc0= mpi_wtime()
    NN_fs= 0d0
    call mpi_allreduce(flocal,NN_fs,1,mpi_double_precision &
         ,mpi_sum,mpi_world,ierr)
    NN_fs= NN_fs /swgt2
    tcl = tcl + (mpi_wtime() -tc0)

!.....only the bottle-neck times are taken into account
    call mpi_reduce(tcl,tcg,1,mpi_double_precision,mpi_max,0 &
         ,mpi_world,ierr)
    call mpi_reduce(tfl,tfg,1,mpi_double_precision,mpi_max,0 &
         ,mpi_world,ierr)
    tcomm= tcomm +tcg
    tfunc= tfunc +tfg
    l1st = .false.
    return
  end function NN_fs
!=======================================================================
  subroutine calc_ef1(smpl,sds)
    use variables
    use minimize, only: mskgfs
    implicit none
    type(mdsys),intent(inout):: smpl
    type(smpldata),intent(inout):: sds
    integer:: natm,ia,ja,ihl0,ihl1,nfcal
    real(8):: tmp,w1,w2,h1,dh1,ddh,t,dg(3)

    natm= smpl%natm
    sds%hl1(1:natm,1:nhl(1))= 0d0
    smpl%epot =0d0

    if( mode.ge.10 ) then
      sds%hl1(1:natm,nhl(1)) = 1d0
    endif
    
!.....energy
    do ia=1,natm
      if( allocated(mskgfs) ) then
        do ihl1=1,mhl(1)
          tmp= 0d0
          do ihl0=1,nhl(0)
            if( mskgfs(ihl0).ne.0 ) cycle
            tmp= tmp +wgt11(ihl0,ihl1) *sds%gsf(ia,ihl0)
          enddo
          sds%hl1(ia,ihl1)= sigmoid(tmp)
        enddo
      else
        do ihl1=1,mhl(1)
          tmp= 0d0
          do ihl0=1,nhl(0)
            tmp= tmp +wgt11(ihl0,ihl1) *sds%gsf(ia,ihl0)
          enddo
          sds%hl1(ia,ihl1)= sigmoid(tmp)
        enddo
      endif
      do ihl1=1,nhl(1)
        smpl%epot= smpl%epot &
             +wgt12(ihl1)*(sds%hl1(ia,ihl1)-0.5d0)
      enddo
    enddo

    !.....forces
    if( .not.lfmatch ) return
    nfcal= smpl%nfcal
    if( nfcal.eq.0 ) return
    smpl%fa(1:3,1:natm)= 0d0
    if( allocated(mskgfs) ) then
      do ihl1=1,mhl(1)
        w2= wgt12(ihl1)
        do ihl0=1,mhl(0)
          if( mskgfs(ihl0).ne.0 ) cycle
          w1= wgt11(ihl0,ihl1)
          do ja=1,natm
            h1= sds%hl1(ja,ihl1)
            dh1= h1*(1d0-h1)
            t= w1*w2 *dh1
            do ia=1,natm
              dg(1:3)=sds%dgsf(1:3,ia,ja,ihl0)
              smpl%fa(1:3,ia)= smpl%fa(1:3,ia) &
                   -t *dg(1:3)
            enddo
          enddo
        enddo
      enddo
!!$    else if( fred.ge.0d0 .or. &
!!$         (nfpsmpl.gt.0 .and. nfpsmpl.lt.natm) ) then
    else if( nfcal.lt.natm ) then
      do ia=1,natm
        if( smpl%ifcal(ia).eq.0 ) cycle
        do ihl1=1,mhl(1)
          w2= wgt12(ihl1)
          do ihl0=1,mhl(0)
            w1= wgt11(ihl0,ihl1)
            do ja=1,natm
              h1= sds%hl1(ja,ihl1)
              dh1= h1*(1d0-h1)
              t= w1*w2 *dh1
              dg(1:3)=sds%dgsf(1:3,ia,ja,ihl0)
              smpl%fa(1:3,ia)= smpl%fa(1:3,ia) &
                   -t *dg(1:3)
            enddo
          enddo
        enddo
      enddo
    else
      do ihl1=1,mhl(1)
        w2= wgt12(ihl1)
        do ihl0=1,mhl(0)
          w1= wgt11(ihl0,ihl1)
          do ja=1,natm
            h1= sds%hl1(ja,ihl1)
            dh1= h1*(1d0-h1)
            t= w1*w2 *dh1
            do ia=1,natm
              dg(1:3)=sds%dgsf(1:3,ia,ja,ihl0)
              smpl%fa(1:3,ia)= smpl%fa(1:3,ia) &
                   -t *dg(1:3)
            enddo
          enddo
        enddo
      enddo
    endif

  end subroutine calc_ef1
!=======================================================================
  subroutine calc_ef2(smpl,sds)
    use variables
    use minimize, only: mskgfs
    implicit none
    type(mdsys),intent(inout):: smpl
    type(smpldata),intent(inout):: sds
    integer:: natm,ia,ja,ihl0,ihl1,ihl2,nfcal,nl0,nl1,nl2
    real(8):: tmp1,tmp2,w1,w2,w3,h1,h2,dh1,dh2,t

    natm= smpl%natm
    sds%hl1(1:natm,1:nhl(1))= 0d0
    sds%hl2(1:natm,1:nhl(2))= 0d0
    smpl%epot= 0d0

    if( mode.ge.10 ) then
      sds%hl1(1:natm,nhl(1)) = 1d0
      sds%hl2(1:natm,nhl(2)) = 1d0
    endif

!.....energy term
    do ia=1,natm
      do ihl2=1,mhl(2)
        tmp2= 0d0
        if( allocated(mskgfs) ) then
          do ihl1=1,mhl(1)
            tmp1= 0d0
            do ihl0=1,nhl(0)
              if( mskgfs(ihl0).ne.0 ) cycle
              tmp1=tmp1 +wgt21(ihl0,ihl1) *sds%gsf(ia,ihl0)
            enddo
            sds%hl1(ia,ihl1)= sigmoid(tmp1)
          enddo
          do ihl1=1,nhl(1)
            tmp2=tmp2 +wgt22(ihl1,ihl2) *(sds%hl1(ia,ihl1)-0.5d0)
          enddo
        else
          do ihl1=1,mhl(1)
            tmp1= 0d0
            do ihl0=1,nhl(0)
              tmp1=tmp1 +wgt21(ihl0,ihl1) *sds%gsf(ia,ihl0)
            enddo
            sds%hl1(ia,ihl1)= sigmoid(tmp1)
            tmp2=tmp2 +wgt22(ihl1,ihl2) *(sds%hl1(ia,ihl1)-0.5d0)
          enddo
          do ihl1=1,nhl(1)
            tmp2=tmp2 +wgt22(ihl1,ihl2) *(sds%hl1(ia,ihl1)-0.5d0)
          enddo
        endif
        sds%hl2(ia,ihl2)= sigmoid(tmp2)
        smpl%epot= smpl%epot &
             +wgt23(ihl2) *(sds%hl2(ia,ihl2)-0.5d0)
      enddo
      do ihl2=1,nhl(2)
        smpl%epot= smpl%epot &
             +wgt23(ihl2) *(sds%hl2(ia,ihl2)-0.5d0)
      enddo
    enddo
    
    !.....force term
    if( .not.lfmatch ) return
    nfcal= smpl%nfcal
    if( nfcal.eq.0 ) return
    smpl%fa(1:3,1:natm)= 0d0
    if( allocated(mskgfs) ) then
      do ihl2=1,mhl(2)
        w3= wgt23(ihl2)
        do ihl1=1,mhl(1)
          w2= wgt22(ihl1,ihl2)
          do ihl0=1,nhl(0)
            if( mskgfs(ihl0).ne.0 ) cycle
            w1= wgt21(ihl0,ihl1)
            do ja=1,natm
              h1= sds%hl1(ja,ihl1)
              h2= sds%hl2(ja,ihl2)
              dh1= h1*(1d0-h1)
              dh2= h2*(1d0-h2)
              t= w3*dh2 *w2*dh1 *w1
              do ia=1,natm
                smpl%fa(1:3,ia)= smpl%fa(1:3,ia) &
                     -t *sds%dgsf(1:3,ia,ja,ihl0)
              enddo
            enddo
          enddo
        enddo
      enddo
!!$    else if( fred.ge.0d0 .or. &
!!$         (nfpsmpl.gt.0 .and. nfpsmpl.lt.natm) ) then
    else if( nfcal.lt.natm ) then
      do ia=1,natm
        if( smpl%ifcal(ia).eq.0 ) cycle
        do ihl2=1,mhl(2)
          w3= wgt23(ihl2)
          do ihl1=1,mhl(1)
            w2= wgt22(ihl1,ihl2)
            do ihl0=1,nhl(0)
              w1= wgt21(ihl0,ihl1)
              do ja=1,natm
                h1= sds%hl1(ja,ihl1)
                h2= sds%hl2(ja,ihl2)
                dh1= h1*(1d0-h1)
                dh2= h2*(1d0-h2)
                t= w3*dh2 *w2*dh1 *w1
                smpl%fa(1:3,ia)= smpl%fa(1:3,ia) &
                     -t *sds%dgsf(1:3,ia,ja,ihl0)
              enddo
            enddo
          enddo
        enddo
      enddo
    else
      do ihl2=1,mhl(2)
        w3= wgt23(ihl2)
        do ihl1=1,mhl(1)
          w2= wgt22(ihl1,ihl2)
          do ihl0=1,nhl(0)
            w1= wgt21(ihl0,ihl1)
            do ja=1,natm
              h1= sds%hl1(ja,ihl1)
              h2= sds%hl2(ja,ihl2)
              dh1= h1*(1d0-h1)
              dh2= h2*(1d0-h2)
              t= w3*dh2 *w2*dh1 *w1
              do ia=1,natm
                smpl%fa(1:3,ia)= smpl%fa(1:3,ia) &
                     -t *sds%dgsf(1:3,ia,ja,ihl0)
              enddo
            enddo
          enddo
        enddo
      enddo
    endif
  
    
  end subroutine calc_ef2
!=======================================================================
  function NN_grad(ndim,x)
    use variables,only: nsmpl,nsmpl_trn,tgrad,ngrad,tcomm &
         ,samples,mdsys,epse,epsf,swgt2
    use parallel
    use minimize
    implicit none
    integer,intent(in):: ndim
    real(8),intent(in):: x(ndim)
    real(8):: NN_grad(ndim)
    
    integer:: ismpl,i,idim
    real(8),save,allocatable:: gs(:),glocal(:)
    real(8):: gmax,vmax
    real(8):: tcl,tgl,tcg,tgg,tc0,tg0
    type(mdsys):: smpl

    if( .not.allocated(gs) ) allocate(gs(ndim),glocal(ndim))

    ngrad= ngrad +1
    tg0= mpi_wtime()

!!$    NN_grad(1:ndim)= 0d0
    glocal(1:ndim)= 0d0

    do ismpl=isid0,isid1
      if( samples(ismpl)%iclass.ne.1 ) cycle
!!$      write(6,*) ' grad ismpl = ',ismpl
      if( nl.eq.1 ) then
        call grad1(samples(ismpl),sds(ismpl),gs)
      else if( nl.eq.2 ) then
        call grad2(samples(ismpl),sds(ismpl),gs)
      endif
      glocal(1:ndim)= glocal(1:ndim) +gs(1:ndim)
!!$      print *,' ismpl,glocal(16)=',ismpl,glocal(16)
    enddo

    tgl= mpi_wtime() -tg0

    tc0= mpi_wtime()
    NN_grad(1:ndim)= 0d0
    call mpi_allreduce(glocal,NN_grad,ndim,mpi_double_precision &
         ,mpi_sum,mpi_world,ierr)
    tcl= mpi_wtime() -tc0
!    tcomm= tcomm +mpi_wtime() -tc0

    NN_grad(1:ndim)= NN_grad(1:ndim) /swgt2

!.....only the bottle-neck times are taken into account
    call mpi_reduce(tcl,tcg,1,mpi_double_precision,mpi_max,0 &
         ,mpi_world,ierr)
    call mpi_reduce(tgl,tgg,1,mpi_double_precision,mpi_max,0 &
         ,mpi_world,ierr)
    tcomm= tcomm +tcg
    tgrad= tgrad +tgg
    
    return
  end function NN_grad
!=======================================================================
  function NN_gs(ndim,x)
    use variables
    use parallel
    use minimize
    implicit none
    integer,intent(in):: ndim
    real(8),intent(in):: x(ndim)
    real(8):: NN_gs(ndim)
    integer:: i,idim
    real(8),save,allocatable:: gsl(:),gs(:)
    real(8):: gmax,vmax,tg0,tc0
    integer:: ismpl

    if( .not.allocated(gsl) ) allocate(gsl(ndim),gs(ndim))

    ngrad= ngrad +1
    tg0= mpi_wtime()

    gsl(1:ndim)= 0d0

    do i=1,nsgdbsize
      ismpl= ismplsgd(i)
      if( samples(ismpl)%iclass.ne.1 ) cycle
      if( nl.eq.1 ) then
        call grad1(samples(ismpl),sds(ismpl),gs)
      else if( nl.eq.2 ) then
        call grad2(samples(ismpl),sds(ismpl),gs)
      endif
      gsl(1:ndim)= gsl(1:ndim) +gs(1:ndim)
    enddo

    tc0= mpi_wtime()

    NN_gs(1:ndim)= 0d0
    call mpi_allreduce(gsl,NN_gs,ndim,mpi_double_precision &
         ,mpi_sum,mpi_world,ierr)
    tcomm= tcomm +mpi_wtime() -tc0
    NN_gs(1:ndim)= NN_gs(1:ndim) /swgt2

    tgrad= tgrad +mpi_wtime() -tg0
    return
  end function NN_gs
!=======================================================================
  subroutine grad1(smpl,sds,gs)
    use variables
    use minimize, only: mskgfs
    implicit none
    type(mdsys),intent(inout):: smpl
    type(smpldata),intent(inout):: sds
    real(8),intent(inout):: gs(nvars)
    integer:: iv,ivp,nv,ihl1,ia,ja,ihl0,jhl0,natm,ixyz,nfcal
    real(8):: ediff,tmp,h1,w1,w2,dn3i,dh1,ddhg,fscale,eref,swgt,wgtidv
    real(8):: edenom,fdenom
    real(8):: eerr,ferr,ferri
    real(8),save,allocatable:: dgs(:),ab(:),wdg(:,:,:),bms(:,:,:,:)

    if( .not. allocated(dgs) ) then
      allocate(dgs(nvars),ab(3),wdg(3,maxna,maxna) &
           ,bms(3,maxna,maxna,nhl(1)))
    endif

    natm= smpl%natm
    eref= smpl%eref
    eerr= smpl%eerr
    swgt= smpl%wgt
    ediff= (smpl%epot -eref) /natm /eerr
    ediff= 2d0 *ediff /natm /eerr *swgt ! *wgtidv /natm
    gs(1:nvars)= 0d0
    iv= nhl(0)*mhl(1) +nhl(1)
    do ihl1=nhl(1),1,-1
      tmp= 0d0
      do ia=1,natm
        h1= sds%hl1(ia,ihl1)
        tmp= tmp +(h1-0.5d0)
      enddo
      gs(iv)= gs(iv) +ediff*tmp
      iv= iv -1
    enddo
    if( allocated(mskgfs) ) then
      do ihl0=nhl(0),1,-1
        do ihl1=mhl(1),1,-1
          tmp= 0d0
          if( mskgfs(ihl0).ne.0 ) goto 20
          w2= wgt12(ihl1)
          do ia=1,natm
            h1= sds%hl1(ia,ihl1)
            tmp= tmp +w2 *h1*(1d0-h1) *sds%gsf(ia,ihl0)
          enddo
20        continue
          gs(iv)= gs(iv) +ediff*tmp
          iv= iv -1
        enddo
      enddo
    else
      do ihl0=nhl(0),1,-1
        do ihl1=mhl(1),1,-1
          tmp= 0d0
          w2= wgt12(ihl1)
          do ia=1,natm
            h1= sds%hl1(ia,ihl1)
            tmp= tmp +w2 *h1*(1d0-h1) *sds%gsf(ia,ihl0)
          enddo
          gs(iv)= gs(iv) +ediff*tmp
          iv= iv -1
        enddo
      enddo
    endif

    if( .not. lfmatch ) return
    nfcal= smpl%nfcal
    if( nfcal.eq.0 ) return
    dgs(1:nvars)= 0d0
    ferr = smpl%ferr
    ferri= 1d0/ferr
    dn3i= 1d0/(3*nfcal)
    do ia=1,natm
      do ixyz=1,3
        fdiff(ixyz,ia)= (smpl%fa(ixyz,ia) &
             -smpl%fref(ixyz,ia)) !*ferri *ferri *2 *dn3i
      enddo
    enddo
    iv= nhl(0)*mhl(1) +nhl(1)
    if( allocated( mskgfs) ) then
      do ihl1=nhl(1),1,-1
        tmp= 0d0
        if( ihl1.gt.mhl(1) ) goto 30
        do ihl0=1,nhl(0)
          if( mskgfs(ihl0).ne.0 ) cycle
          w1= wgt11(ihl0,ihl1)
          do ja=1,natm
            h1= sds%hl1(ja,ihl1)
            dh1= h1*(1d0-h1)
            do ia=1,natm
              tmp= tmp +w1 *dh1*( &
                   fdiff(1,ia)  *sds%dgsf(1,ia,ja,ihl0) &
                   +fdiff(2,ia) *sds%dgsf(2,ia,ja,ihl0) &
                   +fdiff(3,ia) *sds%dgsf(3,ia,ja,ihl0) &
                   )
            enddo
          enddo
        enddo
30      dgs(iv)= -tmp
        iv= iv -1
      enddo
!!$    else if( fred.ge.0d0 .or. &
!!$         (nfpsmpl.gt.0 .and. nfpsmpl.lt.natm) ) then
    else if( nfcal.lt.natm ) then
      ivp = iv
      do ia=1,natm
        if( smpl%ifcal(ia).eq.0 ) cycle
        iv = ivp
        do ihl1=nhl(1),1,-1
          tmp= 0d0
          if( ihl1.gt.mhl(1) ) goto 40
          do ihl0=1,nhl(0)
            w1= wgt11(ihl0,ihl1)
            do ja=1,natm
              h1= sds%hl1(ja,ihl1)
              dh1= h1*(1d0-h1)
              tmp= tmp +w1 *dh1*( &
                   fdiff(1,ia)  *sds%dgsf(1,ia,ja,ihl0) &
                   +fdiff(2,ia) *sds%dgsf(2,ia,ja,ihl0) &
                   +fdiff(3,ia) *sds%dgsf(3,ia,ja,ihl0) &
                   )
            enddo
          enddo
40        dgs(iv)= dgs(iv) -tmp
          iv= iv -1
        enddo
      enddo
    else
      do ihl1=nhl(1),1,-1
        tmp= 0d0
        if( ihl1.gt.mhl(1) ) goto 50
        do ihl0=1,nhl(0)
          w1= wgt11(ihl0,ihl1)
          do ja=1,natm
            h1= sds%hl1(ja,ihl1)
            dh1= h1*(1d0-h1)
            do ia=1,natm
              tmp= tmp +w1 *dh1*( &
                   fdiff(1,ia)  *sds%dgsf(1,ia,ja,ihl0) &
                   +fdiff(2,ia) *sds%dgsf(2,ia,ja,ihl0) &
                   +fdiff(3,ia) *sds%dgsf(3,ia,ja,ihl0) &
                   )
            enddo
          enddo
        enddo
50      dgs(iv)= -tmp
        iv= iv -1
      enddo
    endif
!.....make bms before computing dgs
    bms(1:3,1:natm,1:natm,1:nhl(1))= 0d0
    if( allocated(mskgfs) ) then
      do ihl1=1,mhl(1)
        do ihl0=1,nhl(0)
          if( mskgfs(ihl0).ne.0 ) cycle
          w1= wgt11(ihl0,ihl1)
          do ja=1,natm
            do ia=1,natm
              bms(1:3,ia,ja,ihl1)= bms(1:3,ia,ja,ihl1) &
                   +w1*sds%dgsf(1:3,ia,ja,ihl0)
            enddo
          enddo
        enddo
      enddo
    else if( fred.ge.0d0 .or. &
         (nfpsmpl.gt.0 .and. nfpsmpl.lt.natm) ) then
      do ia=1,natm
        if( smpl%ifcal(ia).eq.0 ) cycle
        do ihl1=1,mhl(1)
          do ihl0=1,nhl(0)
            w1= wgt11(ihl0,ihl1)
            do ja=1,natm
              bms(1:3,ia,ja,ihl1)= bms(1:3,ia,ja,ihl1) &
                   +w1*sds%dgsf(1:3,ia,ja,ihl0)
            enddo
          enddo
        enddo
      enddo
    else
      do ihl1=1,mhl(1)
        do ihl0=1,nhl(0)
          w1= wgt11(ihl0,ihl1)
          do ja=1,natm
            do ia=1,natm
              bms(1:3,ia,ja,ihl1)= bms(1:3,ia,ja,ihl1) &
                   +w1*sds%dgsf(1:3,ia,ja,ihl0)
            enddo
          enddo
        enddo
      enddo
    endif
!.....then compute dgs wrt w1
    if( allocated(mskgfs) ) then
      do ihl0=nhl(0),1,-1
        do ihl1=mhl(1),1,-1
          tmp= 0d0
          if( mskgfs(ihl0).ne.0 ) goto 10
          w2= wgt12(ihl1)
          do ja=1,natm
            h1= sds%hl1(ja,ihl1)
            dh1= h1*(1d0-h1)
            ddhg= dh1*(1d0-2d0*h1)*sds%gsf(ja,ihl0)
            do ia=1,natm
              ab(1:3)= dh1*sds%dgsf(1:3,ia,ja,ihl0) &
                   +ddhg*bms(1:3,ia,ja,ihl1)
              tmp= tmp +w2 *( &
                   fdiff(1,ia) *ab(1) &
                   +fdiff(2,ia) *ab(2) &
                   +fdiff(3,ia) *ab(3) &
                   )
            enddo
          enddo
10        continue
          dgs(iv)= -tmp
          iv= iv -1
        enddo
      enddo
    else if( fred.ge.0d0 .or. &
         (nfpsmpl.gt.0 .and. nfpsmpl.lt.natm) ) then
      ivp = iv
      do ia=1,natm
        if( smpl%ifcal(ia).eq.0 ) cycle
        iv = ivp
        do ihl0=nhl(0),1,-1
          do ihl1=mhl(1),1,-1
            tmp= 0d0
            w2= wgt12(ihl1)
            do ja=1,natm
              h1= sds%hl1(ja,ihl1)
              dh1= h1*(1d0-h1)
              ddhg= dh1*(1d0-2d0*h1)*sds%gsf(ja,ihl0)
              ab(1:3)= dh1*sds%dgsf(1:3,ia,ja,ihl0) &
                   +ddhg*bms(1:3,ia,ja,ihl1)
              tmp= tmp +w2 *( &
                   fdiff(1,ia) *ab(1) &
                   +fdiff(2,ia) *ab(2) &
                   +fdiff(3,ia) *ab(3) &
                   )
            enddo
            dgs(iv)= dgs(iv) -tmp
            iv= iv -1
          enddo
        enddo
      enddo
    else
      do ihl0=nhl(0),1,-1
        do ihl1=mhl(1),1,-1
          tmp= 0d0
          w2= wgt12(ihl1)
          do ja=1,natm
            h1= sds%hl1(ja,ihl1)
            dh1= h1*(1d0-h1)
            ddhg= dh1*(1d0-2d0*h1)*sds%gsf(ja,ihl0)
            do ia=1,natm
              ab(1:3)= dh1*sds%dgsf(1:3,ia,ja,ihl0) &
                   +ddhg*bms(1:3,ia,ja,ihl1)
              tmp= tmp +w2 *( &
                   fdiff(1,ia) *ab(1) &
                   +fdiff(2,ia) *ab(2) &
                   +fdiff(3,ia) *ab(3) &
                   )
            enddo
          enddo
          dgs(iv)= -tmp
          iv= iv -1
        enddo
      enddo
    endif

    gs(1:nvars)= gs(1:nvars) +dgs(1:nvars)*ferri *ferri *2 *dn3i *swgt
    return
  end subroutine grad1
!=======================================================================
  subroutine grad2(smpl,sds,gs)
    use variables
    use minimize, only: mskgfs
    implicit none
    type(mdsys),intent(inout):: smpl
    type(smpldata),intent(inout):: sds
    real(8),intent(inout):: gs(nvars)
    integer:: iv,ihl0,ihl1,ihl2,ia,ja,natm,nfcal
    real(8):: ediff,tmp,tmp1,tmp2,h1,h2,w1,w2,w3,dn3i,dh1,dh2,t1,t2,t3&
         ,ddh1,ddh2,dh1gsf,fscale,eref,swgt,wgtidv
    real(8):: eerr,ferr,ferri
    real(8),save,allocatable:: dgs(:),w1dg(:,:,:,:),w2sw1dg(:,:,:,:)

    if( .not. allocated(dgs) ) then
      allocate( dgs(nvars),w1dg(3,maxna,maxna,nhl(1)) &
           ,w2sw1dg(3,maxna,maxna,nhl(2)) )
    endif

    natm= smpl%natm
    eref= smpl%eref
    eerr= smpl%eerr
    swgt= smpl%wgt
    ediff= (smpl%epot -eref)*2 /natm/natm /eerr /eerr *swgt
    gs(1:nvars)= 0d0
    iv= nhl(0)*mhl(1) +nhl(1)*mhl(2) +nhl(2)

    do ihl2=nhl(2),1,-1
      tmp= 0d0
      do ia=1,natm
        h2= sds%hl2(ia,ihl2)
        tmp= tmp +(h2-0.5d0)
      enddo
      gs(iv)=gs(iv) +ediff*tmp
      iv=iv -1
    enddo
    do ihl1=nhl(1),1,-1
      do ihl2=mhl(2),1,-1
        tmp= 0d0
        w3= wgt23(ihl2)
        do ia=1,natm
          h2= sds%hl2(ia,ihl2)
          h1= sds%hl1(ia,ihl1)
          tmp= tmp +w3 *h2*(1d0-h2) *(h1-0.5d0)
        enddo
        gs(iv)=gs(iv) +ediff*tmp
        iv=iv -1
      enddo
    enddo
    if( allocated(mskgfs) ) then
      do ihl0=nhl(0),1,-1
        do ihl1=mhl(1),1,-1
          tmp= 0d0
          if( mskgfs(ihl0).ne.0 ) goto 20
          do ia=1,natm
            h1= sds%hl1(ia,ihl1)
            dh1= h1*(1d0-h1)
            dh1gsf= dh1*sds%gsf(ia,ihl0)
            do ihl2=1,nhl(2)
              h2= sds%hl2(ia,ihl2)
              dh2= h2*(1d0-h2)
              w2= wgt22(ihl1,ihl2)
              w3= wgt23(ihl2)
              tmp=tmp +w3*w2 *dh2 *dh1gsf
            enddo
          enddo
20        continue
          gs(iv)=gs(iv) +ediff*tmp
          iv=iv -1
        enddo
      enddo
    else
      do ihl0=nhl(0),1,-1
        do ihl1=mhl(1),1,-1
          tmp= 0d0
          do ia=1,natm
            h1= sds%hl1(ia,ihl1)
            dh1= h1*(1d0-h1)
            dh1gsf= dh1*sds%gsf(ia,ihl0)
            do ihl2=1,nhl(2)
              h2= sds%hl2(ia,ihl2)
              dh2= h2*(1d0-h2)
              w2= wgt22(ihl1,ihl2)
              w3= wgt23(ihl2)
              tmp=tmp +w3*w2 *dh2 *dh1gsf
            enddo
          enddo
          gs(iv)=gs(iv) +ediff*tmp
          iv=iv -1
        enddo
      enddo
    endif

    if( .not. lfmatch ) return
    nfcal= smpl%nfcal
    ferr = smpl%ferr
    ferri= 1d0/ferr
    dgs(1:nvars)= 0d0
    fdiff(1:3,1:natm)= (smpl%fa(1:3,1:natm) &
         -smpl%fref(1:3,1:natm)) *ferri
    dn3i= 1d0/3/natm
    fdiff(1:3,1:natm)= fdiff(1:3,1:natm) *2 *ferri *dn3i *swgt

    iv= nhl(0)*nhl(1) +nhl(1)*nhl(2) +nhl(2)
!.....make w1dg
    w1dg(1:3,1:natm,1:natm,1:nhl(1))= 0d0
    if( allocated(mskgfs) ) then
      do ihl1=1,mhl(1)
        do ihl0=1,nhl(0)
          if( mskgfs(ihl0).ne.0 ) cycle
          w1= wgt21(ihl0,ihl1)
          do ja=1,natm
            do ia=1,natm
              w1dg(1:3,ia,ja,ihl1)= w1dg(1:3,ia,ja,ihl1) &
                   +w1*sds%dgsf(1:3,ia,ja,ihl0)
            enddo
          enddo
        enddo
      enddo
    else
      do ihl1=1,mhl(1)
        do ihl0=1,nhl(0)
          w1= wgt21(ihl0,ihl1)
          do ja=1,natm
            do ia=1,natm
              w1dg(1:3,ia,ja,ihl1)= w1dg(1:3,ia,ja,ihl1) &
                   +w1*sds%dgsf(1:3,ia,ja,ihl0)
            enddo
          enddo
        enddo
      enddo
    endif
!.....make w2sw1dg
    w2sw1dg(1:3,1:natm,1:natm,1:nhl(2))= 0d0
    do ihl2=1,mhl(2)
      do ihl1=1,nhl(1)
        w2= wgt22(ihl1,ihl2)
        do ja=1,natm
          h1= sds%hl1(ja,ihl1)
          dh1= h1*(1d0-h1)
          do ia=1,natm
            w2sw1dg(1:3,ia,ja,ihl2)= w2sw1dg(1:3,ia,ja,ihl2) &
                 +w2*dh1 *w1dg(1:3,ia,ja,ihl1)
          enddo
        enddo
      enddo
    enddo
!.....derivative wrt w3
    do ihl2=nhl(2),1,-1
      tmp= 0d0
      do ja=1,natm
        h2= sds%hl2(ja,ihl2)
        dh2= h2*(1d0-h2)
        do ia=1,natm
          tmp=tmp +dh2 *( &
               fdiff(1,ia)  *w2sw1dg(1,ia,ja,ihl2) &
               +fdiff(2,ia) *w2sw1dg(2,ia,ja,ihl2) &
               +fdiff(3,ia) *w2sw1dg(3,ia,ja,ihl2) &
               )
        enddo
      enddo
      dgs(iv)= -tmp
      iv= iv -1
    enddo
!.....derivative wrt w2
    do ihl1=nhl(1),1,-1
      do ihl2=mhl(2),1,-1
        tmp= 0d0
        w3= wgt23(ihl2)
        do ja=1,natm
          h1= sds%hl1(ja,ihl1)
          dh1= h1*(1d0-h1)
          h2= sds%hl2(ja,ihl2)
          dh2= h2*(1d0-h2)
          ddh2= dh2*(1d0-2d0*h2)
          t1= w3*ddh2*(h1-0.5d0)
          t2= w3*dh1*dh2
          do ia=1,natm
            tmp=tmp +t1*( &
                 fdiff(1,ia) *w2sw1dg(1,ia,ja,ihl2) &
                 +fdiff(2,ia)*w2sw1dg(2,ia,ja,ihl2) &
                 +fdiff(3,ia)*w2sw1dg(3,ia,ja,ihl2) &
                 )
            tmp=tmp +t2*( &
                 fdiff(1,ia) *w1dg(1,ia,ja,ihl1) &
                 +fdiff(2,ia)*w1dg(2,ia,ja,ihl1) &
                 +fdiff(3,ia)*w1dg(3,ia,ja,ihl1) &
                 )
          enddo
        enddo
        dgs(iv)= -tmp
        iv= iv -1
      enddo
    enddo
!.....derivative wrt w1
    if( allocated(mskgfs) ) then
      do ihl0=nhl(0),1,-1
        do ihl1=mhl(1),1,-1
          tmp= 0d0
          if( mskgfs(ihl0).ne.0 ) goto 10
          do ihl2=1,mhl(2)
            w3= wgt23(ihl2)
            w2= wgt22(ihl1,ihl2)
            do ja=1,natm
              h1= sds%hl1(ja,ihl1)
              dh1= h1*(1d0-h1)
              ddh1= dh1*(1d0-2d0*h1)
              h2= sds%hl2(ja,ihl2)
              dh2= h2*(1d0-h2)
              ddh2= dh2*(1d0-2d0*h2)
              t1= w3 *ddh2*w2*dh1*sds%gsf(ja,ihl0)
              t2= w3 *dh2*w2*ddh1*sds%gsf(ja,ihl0)
              t3= w3 *dh2 *w2 *dh1
              do ia=1,natm
                tmp=tmp +t1 *( &
                     fdiff(1,ia)  *w2sw1dg(1,ia,ja,ihl2) &
                     +fdiff(2,ia) *w2sw1dg(2,ia,ja,ihl2) &
                     +fdiff(3,ia) *w2sw1dg(3,ia,ja,ihl2) &
                     )
                tmp=tmp +t2 *( &
                     fdiff(1,ia)  *w1dg(1,ia,ja,ihl1) &
                     +fdiff(2,ia) *w1dg(2,ia,ja,ihl1) &
                     +fdiff(3,ia) *w1dg(3,ia,ja,ihl1) &
                     )
                tmp=tmp +t3 *( &
                     fdiff(1,ia)  *sds%dgsf(1,ia,ja,ihl0) &
                     +fdiff(2,ia) *sds%dgsf(2,ia,ja,ihl0) &
                     +fdiff(3,ia) *sds%dgsf(3,ia,ja,ihl0) &
                     )
              enddo
            enddo
          enddo
10        continue
          dgs(iv)= -tmp
          iv= iv -1
        enddo
      enddo
    else
      do ihl0=nhl(0),1,-1
        do ihl1=mhl(1),1,-1
          tmp= 0d0
          do ihl2=1,mhl(2)
            w3= wgt23(ihl2)
            w2= wgt22(ihl1,ihl2)
            do ja=1,natm
              h1= sds%hl1(ja,ihl1)
              dh1= h1*(1d0-h1)
              ddh1= dh1*(1d0-2d0*h1)
              h2= sds%hl2(ja,ihl2)
              dh2= h2*(1d0-h2)
              ddh2= dh2*(1d0-2d0*h2)
              t1= w3 *ddh2*w2*dh1*sds%gsf(ja,ihl0)
              t2= w3 *dh2*w2*ddh1*sds%gsf(ja,ihl0)
              t3= w3 *dh2 *w2 *dh1
              do ia=1,natm
                tmp=tmp +t1 *( &
                     fdiff(1,ia)  *w2sw1dg(1,ia,ja,ihl2) &
                     +fdiff(2,ia) *w2sw1dg(2,ia,ja,ihl2) &
                     +fdiff(3,ia) *w2sw1dg(3,ia,ja,ihl2) &
                     )
                tmp=tmp +t2 *( &
                     fdiff(1,ia)  *w1dg(1,ia,ja,ihl1) &
                     +fdiff(2,ia) *w1dg(2,ia,ja,ihl1) &
                     +fdiff(3,ia) *w1dg(3,ia,ja,ihl1) &
                     )
                tmp=tmp +t3 *( &
                     fdiff(1,ia)  *sds%dgsf(1,ia,ja,ihl0) &
                     +fdiff(2,ia) *sds%dgsf(2,ia,ja,ihl0) &
                     +fdiff(3,ia) *sds%dgsf(3,ia,ja,ihl0) &
                     )
              enddo
            enddo
          enddo
          dgs(iv)= -tmp
          iv= iv -1
        enddo
      enddo
    endif

    gs(1:nvars)= gs(1:nvars) +dgs(1:nvars)
    return
  end subroutine grad2
!=======================================================================
  subroutine vars2wgts(nvars,vars)
    implicit none
    integer,intent(in):: nvars
    real(8),intent(in):: vars(nvars)
    
    integer:: iv,ihl0,ihl1,ihl2

    if( nl.eq.1 ) then
      if( .not. allocated(wgt11) ) then
        allocate(wgt11(nhl(0),mhl(1)),wgt12(nhl(1)))
      endif
      iv= 0
      do ihl0=1,nhl(0)
        do ihl1=1,mhl(1)
          iv= iv +1
          wgt11(ihl0,ihl1)= vars(iv)
!!$          print *, ' iv,ihl0,ihl1 = ',iv,ihl0,ihl1
        enddo
      enddo
      do ihl1=1,nhl(1)
        iv= iv+1
        wgt12(ihl1)= vars(iv)
!!$        print *, ' iv,ihl1 = ',iv,ihl1
      enddo
    else if( nl.eq.2 ) then
      if( .not. allocated(wgt21) ) then
        allocate(wgt21(nhl(0),mhl(1)),wgt22(nhl(1),mhl(2)) &
             ,wgt23(nhl(2)))
      endif
      iv= 0
      do ihl0=1,nhl(0)
        do ihl1=1,mhl(1)
          iv=iv+1
          wgt21(ihl0,ihl1)= vars(iv)
        enddo
      enddo
      do ihl1=1,nhl(1)
        do ihl2=1,mhl(2)
          iv=iv+1
          wgt22(ihl1,ihl2)= vars(iv)
        enddo
      enddo
      do ihl2=1,nhl(2)
        iv=iv+1
        wgt23(ihl2)= vars(iv)
      enddo
    endif
    
  end subroutine vars2wgts
!=======================================================================
  function factorial(n,m)
!  compute factorial of n, m-times.
    implicit none
    integer,intent(in):: n,m
    real(8):: factorial

    integer:: i

    factorial= 1
    if( n .le. m-1 ) then
      return
    endif
    do i=0,m-1
      factorial= factorial*(n-i)
    enddo
    return
  end function factorial
!=======================================================================
  function sigmoid(x)
    implicit none
    real(8),intent(in):: x
    real(8):: sigmoid

    sigmoid= 1d0/(1d0 +exp(-x))
    return
  end function sigmoid
!=======================================================================
  subroutine get_bases()
    !
    ! Read binary files of out.NN.{gsf,dgsf}, written by pmd.
    !
    use variables
    use parallel
    implicit none

    integer:: itmp,ismpl,natm,ia,ihl0,ja
    character*128:: cdir

    do ismpl=isid0,isid1
      natm= samples(ismpl)%natm
      cdir= samples(ismpl)%cdirname
      !.....gsf
      open(21,file=trim(cmaindir)//'/'//trim(cdir)//'/pmd/out.NN.gsf'&
           ,status='old',form='unformatted')
      read(21) itmp
      do ia=1,natm
        read(21) (sds(ismpl)%gsf(ia,ihl0),ihl0=1,nhl(0))
      enddo
      close(21)
      !.....dgsf
      if( lfmatch ) then
        open(22,file=trim(cmaindir)//'/'//trim(cdir)//'/pmd/out.NN.dgsf'&
             ,status='old',form='unformatted')
        do ia=1,natm
          do ihl0=1,nhl(0)
            read(22) (sds(ismpl)%dgsf(1:3,ja,ia,ihl0),ja=1,natm)
          enddo
        enddo
        close(22)
      endif
    enddo


    if(myid.eq.0) print *, 'get_bases done.'
  end subroutine get_bases
!=======================================================================
  function get_mean_input()
!
! Compute the mean of input symmetric functions.
!
    use variables
    use parallel
    implicit none 
    real(8):: get_mean_input
    integer:: nsuml,nsumg,ismpl,natm,ihl0,ia
    real(8):: gmeanl,gmean

    !.....compute mean value
    gmeanl= 0d0
    nsuml= 0
    do ismpl=isid0,isid1
      natm= samples(ismpl)%natm
      !.....sum up gsf
      do ihl0=1,nhl(0)
        do ia=1,natm
          gmeanl= gmeanl +sds(ismpl)%gsf(ia,ihl0)
          nsuml= nsuml +1
        enddo
      enddo
    enddo
    gmean= 0d0
    nsumg= 0
    call mpi_allreduce(gmeanl,gmean,1,mpi_double_precision &
         ,mpi_sum,mpi_world,ierr)
    call mpi_allreduce(nsuml,nsumg,1,mpi_integer &
         ,mpi_sum,mpi_world,ierr)
    get_mean_input = gmean/nsumg
    return
  end function get_mean_input
!=======================================================================
  function get_variance_input(gmean)
    use variables
    use parallel
    implicit none
    real(8),intent(in):: gmean
    real(8):: get_variance_input
    real(8):: varl,varg
    integer:: nsuml,nsumg
    integer:: ismpl,natm,ihl0,ia

    varl = 0d0
    nsuml= 0
    do ismpl=isid0,isid1
      natm= samples(ismpl)%natm
      !.....sum up gsf
      do ihl0=1,nhl(0)
        do ia=1,natm
          varl = varl &
               +(gmean-sds(ismpl)%gsf(ia,ihl0))&
               *(gmean-sds(ismpl)%gsf(ia,ihl0))
          nsuml= nsuml +1
        enddo
      enddo
    enddo
    varg = 0d0
    nsumg= 0
    call mpi_allreduce(varl,varg,1,mpi_double_precision &
         ,mpi_sum,mpi_world,ierr)
    call mpi_allreduce(nsuml,nsumg,1,mpi_integer &
         ,mpi_sum,mpi_world,ierr)
    varg = varg/(nsumg-1) ! unbiased variance (not sample variance)
    get_variance_input = varg
    return
  end function get_variance_input
!=======================================================================
  subroutine NN_standardize()
    use variables
    use parallel
    implicit none

!.....if the standardize already done, skip
    if( lstandard ) return
    
    if( cnormalize(1:3).eq.'var' ) then
      if(myid.eq.0) print *,'normalize w.r.t. variance'
      call standardize_var()
    else if( cnormalize(1:3).eq.'max' ) then
      if(myid.eq.0) print *,'normalize w.r.t. max not implemented'
      call mpi_finalize(ierr)
      stop
      !call standardize_max()
    endif
  end subroutine NN_standardize
!=======================================================================
  subroutine standardize_max()
!
!  Standardize of inputs is recommended when you use lasso or ridge.
!
    use variables, only: nsmpl,samples,nvars,nalist,vars
    use parallel
    implicit none
    integer:: nsuml,nsumg,ismpl,ia,natm,ihl0,ihl1,iv
    real(8),allocatable:: gmaxl(:),gminl(:)

    allocate(gmax(nhl(0)),gmin(nhl(0)) &
         ,gmaxl(nhl(0)),gminl(nhl(0)))

    gmaxl(1:nhl(0))= 0d0
    gminl(1:nhl(0))= 1d+30
    do ismpl=isid0,isid1
      natm= samples(ismpl)%natm
      !.....sum up gsf
      do ihl0=1,nhl(0)
        do ia=1,natm
          gmaxl(ihl0)= max(gmaxl(ihl0),sds(ismpl)%gsf(ia,ihl0))
          gminl(ihl0)= min(gminl(ihl0),sds(ismpl)%gsf(ia,ihl0))
        enddo
      enddo
    enddo

    gmax(1:nhl(0))= 0d0
    gmin(1:nhl(0))= 1d+30
    call mpi_allreduce(gmaxl,gmax,nhl(0),mpi_double_precision &
         ,mpi_max,mpi_world,ierr)
    call mpi_allreduce(gminl,gmin,nhl(0),mpi_double_precision &
         ,mpi_min,mpi_world,ierr)

!!$    if( myid.eq.0 ) then
!!$      print *,'Max and min of G:'
!!$      do ihl0=1,nhl(0)
!!$        write(6,'(i5,2es12.4)') ihl0,gmax(ihl0),gmin(ihl0)
!!$      enddo
!!$    endif

!.....neglect 0 values
    do ihl0=1,nhl(0)
      if( gmax(ihl0).lt.1d-5 ) gmax(ihl0)=1d0
    enddo

!.....standardize G values
    do ismpl=isid0,isid1
      natm= samples(ismpl)%natm
      allocate(sds(ismpl)%gsfo(natm,nhl(0)))
      do ihl0=1,nhl(0)
        do ia=1,natm
          sds(ismpl)%gsfo(ia,ihl0)= sds(ismpl)%gsf(ia,ihl0)
          sds(ismpl)%gsf(ia,ihl0)= sds(ismpl)%gsf(ia,ihl0) /gmax(ihl0)
        enddo
      enddo
    enddo

    iv=0
    do ihl0=1,nhl(0)
      do ihl1=1,mhl(1)
        iv=iv+1
        vars(iv)= vars(iv)*gmax(ihl0)
      enddo
    enddo

    lstandard= .true.
    deallocate(gmaxl,gminl)
  end subroutine standardize_max
!=======================================================================
  subroutine standardize_norm()
!
!  Standardize of inputs by dividing by L2 norm
!
    use variables, only: nsmpl,samples,nvars,nalist,vars
    use parallel
    implicit none
    integer:: nsuml,nsumg,ismpl,ia,natm,ihl0,ihl1,iv
    real(8),allocatable:: gmaxl(:),gminl(:)

    allocate(gmax(nhl(0)),gmaxl(nhl(0)))

    gmaxl(1:nhl(0))= 0d0
    do ismpl=isid0,isid1
      natm= samples(ismpl)%natm
      !.....sum up gsf
      do ihl0=1,nhl(0)
        do ia=1,natm
          gmaxl(ihl0)= gmaxl(ihl0) +sds(ismpl)%gsf(ia,ihl0)**2
        enddo
      enddo
    enddo

    gmax(1:nhl(0))= 0d0
    call mpi_allreduce(gmaxl,gmax,nhl(0),mpi_double_precision &
         ,mpi_sum,mpi_world,ierr)
    do ihl0=1,nhl(0)
      gmax(ihl0)= sqrt(gmax(ihl0))
      if( gmax(ihl0).lt.1d-8 ) gmax(ihl0)= 1d0
    enddo

!.....standardize G values
    do ismpl=isid0,isid1
      natm= samples(ismpl)%natm
      allocate(sds(ismpl)%gsfo(natm,nhl(0)))
      do ihl0=1,nhl(0)
        do ia=1,natm
          sds(ismpl)%gsfo(ia,ihl0)= sds(ismpl)%gsf(ia,ihl0)
          sds(ismpl)%gsf(ia,ihl0)= sds(ismpl)%gsf(ia,ihl0) /gmax(ihl0)
        enddo
      enddo
    enddo

    iv=0
    do ihl0=1,nhl(0)
      do ihl1=1,mhl(1)
        iv=iv+1
        vars(iv)= vars(iv)*gmax(ihl0)
      enddo
    enddo

    lstandard= .true.
    deallocate(gmaxl)
    if(myid.eq.0) print *,'standardize done.'
  end subroutine standardize_norm
!=======================================================================
  subroutine standardize_var()
!
!  Standardize of inputs is requied when you use lasso or ridge.
!
    use variables, only: nsmpl,samples,nvars,nalist,vars
    use parallel
    implicit none
    integer:: ismpl,ia,natm,ihl0,ihl1,iv
    real(8):: sgm,sgmi
    logical,save:: l1st= .true.

    if( l1st ) then
      do ismpl=isid0,isid1
        natm= samples(ismpl)%natm
        allocate(sds(ismpl)%gsfo(natm,nhl(0)))
      enddo

      sgm = sqrt(gsfvar)
      sgmi= 1d0/sgm
!.....standardize G values
      do ismpl=isid0,isid1
        natm= samples(ismpl)%natm
        do ihl0=1,nhl(0)
          do ia=1,natm
            sds(ismpl)%gsfo(ia,ihl0)= sds(ismpl)%gsf(ia,ihl0)
            sds(ismpl)%gsf(ia,ihl0)= sds(ismpl)%gsf(ia,ihl0) *sgmi
          enddo
        enddo
      enddo
    endif

    iv=0
    do ihl0=1,nhl(0)
      do ihl1=1,mhl(1)
        iv=iv+1
        vars(iv)= vars(iv)*sgm
      enddo
    enddo

    l1st= .false.
    lstandard= .true.
  end subroutine standardize_var
!=======================================================================
  subroutine NN_restore_standard()
!
!  Restore weights by inverse standardization
!
    use variables
    use parallel
    implicit none
    integer:: iv,ihl0,ihl1
    real(8):: sgm

    if( .not. lstandard ) then
!!$      if(myid.eq.0) print *,'NN_restore_standard not needed.'
      return
    endif

    if( cnormalize.eq.'var' ) then
      sgm = sqrt(gsfvar)
      iv= 0
      do ihl0=1,nhl(0)
        do ihl1=1,mhl(1)
          iv=iv+1
          vars(iv)= vars(iv)/sgm
        enddo
      enddo
    else if( cnormalize.eq.'max' ) then
      if(myid.eq.0) then
        print *,'It is not implemented yet.'
        print *,'Something wrong if you come here...'
      endif
      call mpi_finalize(ierr)
      stop
    endif

    !if(myid.eq.0) print *,'NN_restore_standard done.'
    lstandard = .false.
  end subroutine NN_restore_standard
!=======================================================================
  subroutine NN_analyze(cadd)
!
!  Get which input nodes are more/less important.
!
    use variables
    use parallel
    implicit none
    character(len=*),intent(in):: cadd
    character(len=14),parameter:: cfname= 'out.NN_analyze'

    integer,parameter:: ionum=  30
    integer,allocatable:: icmb2(:,:),icmb3(:,:,:),itype(:),nctype(:)&
         ,nstype(:)
    real(8),allocatable:: sumv(:),cnst(:,:),sumvv(:)
    integer:: i,j,k,l,i2,i3,isf,iv,ic,ihl0,ihl1,itmp,icmb(3)

    allocate(sumv(mhl(0)))
    call eval_1st_layer(sumv)


    if( myid.eq.0 ) then
!!$!.....read in.comb.NN file
!!$      allocate(icmb2(nsp,nsp),icmb3(nsp,nsp,nsp))
!!$      open(ionum,file=trim(cmaindir)//'/'//trim(cmbfname),status='old')
!!$      do i2=1,ncmb2
!!$        read(ionum,*) i,j,icmb2(i,j)
!!$        icmb2(j,i)= icmb2(i,j)
!!$      enddo
!!$      do i3=1,ncmb3
!!$        read(ionum,*) i,j,k,icmb3(i,j,k)
!!$        icmb3(i,k,j)= icmb3(i,j,k)
!!$      enddo
!!$      close(ionum)
!.....read in.const.NN file
      allocate(itype(nhl(0)),cnst(2,nhl(0)),nctype(200),nstype(200))
      nctype(1)= 2   ! Gaussian
      nctype(2)= 1   ! cosine
      nctype(3)= 1   ! polynomial
      nctype(4)= 2   ! Morse
      nctype(101)= 1 ! angular
      nstype(1:100)= 2
      nstype(101:200)= 3
      open(ionum,file=trim(cmaindir)//'/'//trim(ccfname),status='old')
      read(ionum,*) itmp
      do isf=1,mhl(0)
        read(ionum,*) itype(isf),(icmb(k),k=1,nstype(itype(isf))) &
             ,(cnst(j,isf),j=1,nctype(itype(isf)))
      enddo
      close(ionum)

      open(ionum+1,file=cfname//"."//trim(cadd),status='replace')
      iv=0
      do isf=1,mhl(0)
        write(ionum+1,'(2i5,f24.14)') isf,itype(isf),sumv(isf)
      enddo
!!$!.....about 2body terms
!!$      do ihl0=1,nsf2*ncmb2
!!$        isf= mod(ihl0-1,nsf2)+1
!!$        ic = (ihl0-1)/nsf2 +1
!!$        do i=1,nsp
!!$          do j=1,nsp
!!$            if(ic.eq.icmb2(i,j)) goto 10
!!$          enddo
!!$        enddo
!!$10      continue
!!$        write(ionum+1,'(f24.14,2x,i1,"-",i1,":",i5,1es12.4,2i8)') &
!!$             sumv(ihl0),i,j,itype(isf),cnst(1,isf),ihl0,ic
!!$      enddo
!!$!.....about 3body terms
!!$      do ihl0=nsf2*ncmb2+1,nsf2*ncmb2+nsf3*ncmb3
!!$        isf= nsf2+mod(ihl0-nsf2*ncmb2-1,nsf3)+1
!!$        ic = (ihl0-nsf2*ncmb2-1)/nsf3 +1
!!$        do i=1,nsp
!!$          do j=1,nsp
!!$            do k=1,nsp
!!$              if(ic.eq.icmb3(i,j,k)) goto 20
!!$            enddo
!!$          enddo
!!$        enddo
!!$20      continue
!!$        write(ionum+1,'(f24.14,2x,i1,"-",i1,"-",i1,":",i5,1es12.4,2i8)') &
!!$             sumv(ihl0),i,j,k,itype(isf),cnst(1,isf),ihl0,ic
!!$      enddo
      close(ionum+1)

!!$      deallocate(icmb2,icmb3)
    endif

    deallocate(sumv)
    call mpi_barrier(mpi_world,ierr)

    return
  end subroutine NN_analyze
!=======================================================================
  subroutine eval_1st_layer(sumv)
!
!  Evaluate contributions from input nodes. Not only the sum of weights,
!  contribution means weights*(gsf value) of every atoms in every samples.
!
    use variables
    use parallel
    implicit none
    real(8),intent(out):: sumv(mhl(0))
    real(8),save,allocatable:: sumvl(:)

    integer:: ismpl,ia,ihl0,ihl1,natm
    integer,save,allocatable:: ncnt(:)

    if( .not.allocated(sumvl) ) then
      allocate(sumvl(mhl(0)),ncnt(mhl(0)))
    endif
    sumvl(1:nhl(0))= 0d0
    ncnt(1:nhl(0))= 0
    if( nl.eq.1 ) then ! 1-layer NN
      do ismpl=isid0,isid1
        natm= samples(ismpl)%natm
        do ia=1,natm
          do ihl0=1,nhl(0)
            do ihl1=1,mhl(1)
              sumvl(ihl0)= sumvl(ihl0) &
                   +abs(wgt11(ihl0,ihl1) *sds(ismpl)%gsf(ia,ihl0))
              ncnt(ihl0)= ncnt(ihl0) +1
            enddo
          enddo
        enddo
      enddo
    else if( nl.eq.2 ) then ! 2-layer NN
      do ismpl=isid0,isid1
        natm= samples(ismpl)%natm
        do ia=1,natm
          do ihl0=1,nhl(0)
            do ihl1=1,mhl(1)
              sumvl(ihl0)= sumvl(ihl0) &
                   +abs(wgt21(ihl0,ihl1) *sds(ismpl)%gsf(ia,ihl0))
              ncnt(ihl0)= ncnt(ihl0) +1
            enddo
          enddo
        enddo
      enddo
    endif
    
    do ihl0=1,mhl(0)
      sumvl(ihl0)= sumvl(ihl0)/ncnt(ihl0)
    enddo

    sumv(1:mhl(0))= 0d0
    call mpi_reduce(sumvl,sumv,mhl(0),mpi_double_precision &
         ,mpi_sum,0,mpi_world,ierr)
    
  end subroutine eval_1st_layer
!=======================================================================
  subroutine count_nterms()
    use variables
    use parallel
    implicit none
    integer:: ismpl,natm,nttrnl,nttstl
    type(mdsys):: smpl

    nttrnl = 0
    nttstl = 0
    do ismpl=isid0,isid1
      smpl = samples(ismpl)
      natm = smpl%natm
      if( smpl%iclass.eq.1 ) then
        nttrnl = nttrnl + 1 +1
      else
        nttstl = nttstl + 1 +1
      endif
    enddo

    nterm_trn = 0
    nterm_tst = 0
    call mpi_allreduce(nttrnl,nterm_trn,1,mpi_integer &
         ,mpi_sum,mpi_world,ierr)
    call mpi_allreduce(nttstl,nterm_tst,1,mpi_integer &
         ,mpi_sum,mpi_world,ierr)
    
  end subroutine count_nterms
!=======================================================================
end module NN
