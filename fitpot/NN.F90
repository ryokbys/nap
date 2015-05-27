module NN
!-----------------------------------------------------------------------
!                        Time-stamp: <2015-03-14 11:04:41 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
!.....parameter file name
  character(128),parameter:: cpfname= 'in.params.NN'
  character(128),parameter:: ccfname='in.const.NN'
  character(128),parameter:: cmbfname='in.comb.NN'
  integer,parameter:: maxnl= 2
  integer,save:: nl,nsp,nsfc,nsf2,nsf3,ncmb2,ncmb3
  integer,save:: nhl(0:maxnl+1)
  integer,save,allocatable:: nwgt(:)
  real(8),save,allocatable:: wgt11(:,:),wgt12(:)
  real(8),save,allocatable:: wgt21(:,:),wgt22(:,:),wgt23(:)

  type smpldata
    real(8),allocatable:: gsf(:,:),dgsf(:,:,:,:)
    real(8),allocatable:: hl1(:,:),hl2(:,:)
    real(8),allocatable:: gsfo(:,:)
  end type smpldata

  type(smpldata),save,allocatable:: sds(:)

  integer,save:: maxna
  real(8),save,allocatable:: fdiff(:,:)
  real(8),save,allocatable:: gmax(:),gmin(:)

  logical,save:: lstandard= .false.

contains
!=======================================================================
  subroutine NN_init()
    use variables
    use parallel
    use minimize
    implicit none 
    integer:: itmp,i,nw,natm,ismpl,ihl0,ihl1

    tfunc= 0d0
    tgrad= 0d0
    nfunc= 0
    ngrad= 0

    !.....read in.const.NN to get nl,nsp,nhl(:)
    if( myid.eq.0 ) then
      open(20,file=trim(cmaindir)//'/'//trim(ccfname),status='old')
      read(20,*) nl,nsp,nhl(0:nl)
      nsf2= 0
      nsf3= 0
      do while(.true.)
        read(20,*,end=10) itmp
        if( itmp.le.100 ) then
          nsf2=nsf2+1
        else if( itmp.le.200 ) then
          nsf3=nsf3+1
        endif
      enddo
10    close(20)
      nsfc= nhl(0)
      nhl(nl+1)= 1
      write(6,'(a,5i5)') 'nhl(0:nl+1)=',nhl(0:nl+1)
    endif
    call mpi_bcast(nl,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(nsp,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(nhl,nl+2,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(nsf2,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(nsf3,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(nsfc,1,mpi_integer,0,mpi_world,ierr)

!.....calc number of weights
    ncmb2= nsp +factorial(nsp,2)/2
    ncmb3= nsp*ncmb2
    nhl(0)= nsf2*ncmb2 +nsf3*ncmb3
    allocate(nwgt(nl+1))
    nw= 0
    do i=1,nl+1
      nwgt(i)= nhl(i-1)*nhl(i)
      nw= nw +nwgt(i)
    enddo
    if(myid.eq.0) then
      print *,'nsf2,nsf3,ncmb2,ncmb3=',nsf2,nsf3,ncmb2,ncmb3
      print *,'nhl(0:nl+1),nw=',nhl(0:nl+1),nw
    endif

!.....training set
    allocate(sds(isid0:isid1))
    do ismpl=isid0,isid1
      natm= samples(ismpl)%natm
      allocate( sds(ismpl)%gsf(natm,nhl(0)) &
           ,sds(ismpl)%dgsf(3,natm,natm,nhl(0)) )
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

    call standardize_var()

!.....make groups for group lasso
    if( trim(cpena).eq.'glasso' ) then
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
      do i=nhl(0)*nhl(1)+1,nw
        iglid(i)= -1
      enddo
    endif

    if( myid.eq.0 ) print *, 'NN_init done.'
  end subroutine NN_init
!=======================================================================
  function NN_func(ndim,x)
    use variables,only:nsmpl,nsmpl_trn,samples,nprcs,tfunc &
         ,lfmatch,lfscale,fscl,nfunc,tcomm,lswgt,swbeta,mdsys
    use parallel
    use minimize
    implicit none
    integer,intent(in):: ndim
    real(8),intent(in):: x(ndim)
    real(8):: NN_func

    integer:: ismpl,natm,ia,ixyz,idim
    real(8):: dn3i,ediff,tf0,tc0,fscale,eref,swgt
    real(8):: flocal
    type(mdsys):: smpl

    nfunc= nfunc +1

    call mpi_bcast(x,ndim,mpi_double_precision,0,mpi_world,ierr)
    tf0= mpi_wtime()
    call vars2wgts(ndim,x)
    
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
      ediff= (smpl%epot -eref)/natm
      ediff= ediff*ediff
      swgt= 1d0
      if( lswgt ) then
        swgt= exp(-eref/natm*swbeta)
      endif
      flocal= flocal +ediff*swgt
      if( .not. lfmatch ) cycle
      fdiff(1:3,1:natm)= (smpl%fa(1:3,1:natm) &
           -smpl%fref(1:3,1:natm))
      dn3i= 1d0 /(3*natm)
      fscale= 1d0
      !.....force-scale makes force contribution same order to energy
      if( lfscale ) fscale= fscl
      fdiff(1:3,1:natm)= fdiff(1:3,1:natm)*fdiff(1:3,1:natm) &
           *dn3i *fscale *swgt
      do ia=1,natm
        do ixyz=1,3
          flocal= flocal +fdiff(ixyz,ia)
        enddo
      enddo
    enddo

    tc0= mpi_wtime()
    NN_func= 0d0
    call mpi_allreduce(flocal,NN_func,1,mpi_double_precision &
         ,mpi_sum,mpi_world,ierr)
    tcomm= tcomm +mpi_wtime() -tc0
    NN_func= NN_func/nsmpl_trn

    tfunc= tfunc +mpi_wtime() -tf0
    return
  end function NN_func
!=======================================================================
!!$  function NN_func_tst(ndim,x)
!!$    use variables, only:nsmpl_tst,nprcs,tfunc,smpl_tst &
!!$         ,lfmatch,lfscale,fscl,nfunc,tcomm,lswgt,swbeta
!!$    use parallel
!!$    use minimize
!!$    implicit none
!!$    integer,intent(in):: ndim
!!$    real(8),intent(in):: x(ndim)
!!$    real(8):: NN_func_tst
!!$
!!$    integer:: ismpl,natm,ia,ixyz,idim
!!$    real(8):: dn3i,ediff,tf0,tc0,fscale,eref,swgt
!!$    real(8):: flocal
!!$
!!$    nfunc= nfunc +1
!!$
!!$    call mpi_bcast(x,ndim,mpi_double_precision,0,mpi_world,ierr)
!!$    tf0= mpi_wtime()
!!$    call vars2wgts(ndim,x)
!!$    
!!$    do ismpl=isid0_tst,isid1_tst
!!$      if( nl.eq.1 ) then
!!$        call calc_ef1(smpl_tst(ismpl),sds_tst(ismpl))
!!$      else if( nl.eq.2 ) then
!!$        call calc_ef2(smpl_tst(ismpl),sds_tst(ismpl))
!!$      endif
!!$    enddo
!!$
!!$    flocal= 0d0
!!$    do ismpl=isid0_tst,isid1_tst
!!$      natm= smpl_tst(ismpl)%natm
!!$      eref= smpl_tst(ismpl)%eref
!!$      ediff= (smpl_tst(ismpl)%epot -eref)/natm
!!$      ediff= ediff*ediff
!!$      swgt= 1d0
!!$      if( lswgt ) then
!!$        swgt= exp(-eref/natm*swbeta)
!!$      endif
!!$      flocal= flocal +ediff*swgt
!!$      if( .not. lfmatch ) cycle
!!$      fdiff(1:3,1:natm)= (smpl_tst(ismpl)%fa(1:3,1:natm) &
!!$           -smpl_tst(ismpl)%fref(1:3,1:natm))
!!$      dn3i= 1d0 /(3*natm)
!!$      fscale= 1d0
!!$      !.....force-scale makes force contribution same order to energy
!!$      if( lfscale ) fscale= fscl
!!$      fdiff(1:3,1:natm)= fdiff(1:3,1:natm)*fdiff(1:3,1:natm) &
!!$           *dn3i *fscale *swgt
!!$      do ia=1,natm
!!$        do ixyz=1,3
!!$          flocal= flocal +fdiff(ixyz,ia)
!!$        enddo
!!$      enddo
!!$    enddo
!!$
!!$    tc0= mpi_wtime()
!!$    NN_func_tst= 0d0
!!$    call mpi_allreduce(flocal,NN_func_tst,1,mpi_double_precision &
!!$         ,mpi_sum,mpi_world,ierr)
!!$    tcomm= tcomm +mpi_wtime() -tc0
!!$    NN_func_tst= NN_func_tst/nsmpl_tst
!!$
!!$    tfunc= tfunc +mpi_wtime() -tf0
!!$    return
!!$  end function NN_func_tst
!!$!=======================================================================
  function NN_fs(ndim,x)
    use variables
    use parallel
    use minimize
    implicit none
    integer,intent(in):: ndim
    real(8),intent(in):: x(ndim)
    real(8):: NN_fs
    integer:: natm,ia,ixyz,idim
    real(8):: dn3i,ediff,tf0,fscale,eref,swgt,flocal,tc0
    integer:: ismpl
    common /samplei/ ismpl
    type(mdsys):: smpl

    nfunc=nfunc +1
    tf0= mpi_wtime()
    call vars2wgts(ndim,x)

    smpl= samples(ismpl)
    if( nl.eq.1 ) then
      call calc_ef1(smpl,sds(ismpl))
    else if( nl.eq.2 ) then
      call calc_ef2(smpl,sds(ismpl))
    endif

    flocal= 0d0
    if( smpl%iclass.ne.1 ) goto 888
    natm= smpl%natm
    eref= smpl%eref
    ediff= (smpl%epot -eref)/natm
    ediff= ediff*ediff
    swgt= 1d0
    if( lswgt ) then
      swgt= exp(-eref/natm*swbeta)
    endif
    flocal= flocal +ediff*swgt
    if( .not. lfmatch ) goto 888
    fdiff(1:3,1:natm)= (smpl%fa(1:3,1:natm) &
         -smpl%fref(1:3,1:natm))
    dn3i= 1d0 /(3*natm)
    fscale= 1d0
!.....force-scale makes force contribution same order to energy
    if( lfscale ) fscale= fscl
    fdiff(1:3,1:natm)= fdiff(1:3,1:natm)*fdiff(1:3,1:natm) &
         *dn3i *fscale *swgt
    do ia=1,natm
      do ixyz=1,3
        flocal= flocal +fdiff(ixyz,ia)
      enddo
    enddo

888 continue

    tc0= mpi_wtime()
    NN_fs= 0d0
!!$    print *,'here01,myid,flocal=',myid,flocal
    call mpi_allreduce(flocal,NN_fs,1,mpi_double_precision &
         ,mpi_sum,mpi_world,ierr)
!!$    print *,'here02,myid,NN_func=',myid,NN_func
    tcomm= tcomm +mpi_wtime() -tc0
    NN_fs= NN_fs/nnode

999 tfunc= tfunc +mpi_wtime() -tf0
    return
  end function NN_fs
!=======================================================================
  subroutine calc_ef1(smpl,sds)
    use variables
    implicit none
    type(mdsys),intent(inout):: smpl
    type(smpldata),intent(inout):: sds
    integer:: natm,ia,ja,ihl0,ihl1
    real(8):: tmp,w1,w2,h1,dh1,ddh,t,dg(3)
    
    natm= smpl%natm
    sds%hl1(1:natm,1:nhl(1))= 0d0
    smpl%epot =0d0
    !.....energy
    do ia=1,natm
      do ihl1=1,nhl(1)
        tmp= 0d0
        do ihl0=1,nhl(0)
          tmp= tmp +wgt11(ihl0,ihl1) *sds%gsf(ia,ihl0)
        enddo
        sds%hl1(ia,ihl1)= sigmoid(tmp)
      enddo
      do ihl1=1,nhl(1)
        smpl%epot= smpl%epot &
             +wgt12(ihl1)*(sds%hl1(ia,ihl1)-0.5d0)
      enddo
    enddo

    !.....forces
    if( .not.lfmatch ) return
    smpl%fa(1:3,1:natm)= 0d0
    do ihl1=1,nhl(1)
      w2= wgt12(ihl1)
      do ihl0=1,nhl(0)
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

  end subroutine calc_ef1
!=======================================================================
  subroutine calc_ef2(smpl,sds)
    use variables
    implicit none
    type(mdsys),intent(inout):: smpl
    type(smpldata),intent(inout):: sds
    integer:: natm,ia,ja,ihl0,ihl1,ihl2
    real(8):: tmp1,tmp2,w1,w2,w3,h1,h2,dh1,dh2,t

    natm= smpl%natm
    sds%hl1(1:natm,1:nhl(1))= 0d0
    sds%hl2(1:natm,1:nhl(2))= 0d0
    smpl%epot= 0d0

    !.....energy term
    do ia=1,natm
      do ihl2=1,nhl(2)
        tmp2= 0d0
        do ihl1=1,nhl(1)
          tmp1= 0d0
          do ihl0=1,nhl(0)
            tmp1=tmp1 +wgt21(ihl0,ihl1) *sds%gsf(ia,ihl0)
          enddo
          sds%hl1(ia,ihl1)= sigmoid(tmp1)
          tmp2=tmp2 +wgt22(ihl1,ihl2) *(sds%hl1(ia,ihl1)-0.5d0)
        enddo
        sds%hl2(ia,ihl2)= sigmoid(tmp2)
        smpl%epot= smpl%epot &
             +wgt23(ihl2) *(sds%hl2(ia,ihl2)-0.5d0)
      enddo
    enddo
    
    !.....force term
    if( .not.lfmatch ) return
    smpl%fa(1:3,1:natm)= 0d0
    do ihl2=1,nhl(2)
      w3= wgt23(ihl2)
      do ihl1=1,nhl(1)
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
    
  end subroutine calc_ef2
!=======================================================================
  function NN_grad(ndim,x)
    use variables,only: nsmpl,nsmpl_trn,tgrad,ngrad,tcomm,samples,mdsys
    use parallel
    use minimize
    implicit none
    integer,intent(in):: ndim
    real(8),intent(in):: x(ndim)
    real(8):: NN_grad(ndim)
    
    integer:: ismpl,i,idim
    real(8),save,allocatable:: gs(:),glocal(:)
    real(8):: gmax,vmax,tc0,tg0
    type(mdsys):: smpl

    if( .not.allocated(gs) ) allocate(gs(ndim),glocal(ndim))

    ngrad= ngrad +1
    tg0= mpi_wtime()

!!$    NN_grad(1:ndim)= 0d0
    glocal(1:ndim)= 0d0

    do ismpl=isid0,isid1
      if( samples(ismpl)%iclass.ne.1 ) cycle
      if( nl.eq.1 ) then
        call grad1(samples(ismpl),sds(ismpl),gs)
      else if( nl.eq.2 ) then
        call grad2(samples(ismpl),sds(ismpl),gs)
      endif
!!$      NN_grad(1:ndim)= NN_grad(1:ndim) +gs(1:ndim)
      glocal(1:ndim)= glocal(1:ndim) +gs(1:ndim)
    enddo

    tc0= mpi_wtime()
    NN_grad(1:ndim)= 0d0
    call mpi_allreduce(glocal,NN_grad,ndim,mpi_double_precision &
         ,mpi_sum,mpi_world,ierr)
    tcomm= tcomm +mpi_wtime() -tc0

    NN_grad(1:ndim)= NN_grad(1:ndim)/nsmpl_trn

!!$    if( lgscale ) then
!!$      gmax= 0d0
!!$      vmax= 0d0
!!$      do i=1,ndim
!!$        vmax= max(vmax,abs(vars(i)))
!!$        gmax= max(gmax,abs(gval(i)))
!!$      enddo
!!$      gval(1:ndim)= gval(1:ndim)/gmax *gscl*vmax
!!$    endif

    tgrad= tgrad +mpi_wtime() -tg0
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
    real(8),save,allocatable:: gsl(:)
    real(8):: gmax,vmax,tg0,tc0
    integer:: ismpl
    common /samplei/ ismpl

    if( .not.allocated(gsl) ) allocate(gsl(nvars))

    ngrad= ngrad +1
    tg0= mpi_wtime()

    gsl(1:nvars)= 0d0
    if( samples(ismpl)%iclass.ne.1 ) goto 888
    if( nl.eq.1 ) then
      call grad1(samples(ismpl),sds(ismpl),gsl)
    else if( nl.eq.2 ) then
      call grad2(samples(ismpl),sds(ismpl),gsl)
    endif

888 continue
    tc0= mpi_wtime()
    NN_gs(1:ndim)= 0d0
    call mpi_allreduce(gsl,NN_gs,ndim,mpi_double_precision &
         ,mpi_sum,mpi_world,ierr)
    tcomm= tcomm +mpi_wtime() -tc0
    NN_gs(1:ndim)= NN_gs(1:ndim)/nnode

    tgrad= tgrad +mpi_wtime() -tg0
    return
  end function NN_gs
!=======================================================================
  subroutine grad1(smpl,sds,gs)
    use variables
    implicit none
    type(mdsys),intent(inout):: smpl
    type(smpldata),intent(inout):: sds
    real(8),intent(inout):: gs(nvars)
    integer:: iv,ihl1,ia,ja,ihl0,jhl0,natm
    real(8):: ediff,tmp,h1,w1,w2,dn3i,dh1,ddhg,fscale,eref,swgt
    real(8),save,allocatable:: dgs(:),ab(:),wdg(:,:,:),bms(:,:,:,:)

    if( .not. allocated(dgs) ) then
      allocate(dgs(nvars),ab(3),wdg(3,maxna,maxna) &
           ,bms(3,maxna,maxna,nhl(1)))
    endif

    natm= smpl%natm
    eref= smpl%eref
    swgt= 1d0
    if( lswgt ) then
      swgt= exp(-eref/natm*swbeta)
    endif
    ediff= (smpl%epot -eref)*2 /natm/natm *swgt
!!$    print *,'ediff=',ediff*natm
    gs(1:nvars)= 0d0
    iv= nhl(0)*nhl(1) +nhl(1)
    do ihl1=nhl(1),1,-1
      tmp= 0d0
      do ia=1,natm
        h1= sds%hl1(ia,ihl1)
        tmp= tmp +(h1-0.5d0)
      enddo
      gs(iv)= gs(iv) +ediff*tmp
      iv= iv -1
    enddo
    do ihl0=nhl(0),1,-1
      do ihl1=nhl(1),1,-1
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

    if( .not. lfmatch ) return
    dgs(1:nvars)= 0d0
    fdiff(1:3,1:natm)= (smpl%fa(1:3,1:natm) &
         -smpl%fref(1:3,1:natm))
    dn3i= 1d0/(3*natm)
    fscale= 1d0
    if( lfscale ) fscale= fscl
    fdiff(1:3,1:natm)= fdiff(1:3,1:natm) *2 *dn3i *fscale *swgt
    iv= nhl(0)*nhl(1) +nhl(1)
    do ihl1=nhl(1),1,-1
      tmp= 0d0
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
      dgs(iv)= -tmp
      iv= iv -1
    enddo
!.....make bms before computing dgs
    bms(1:3,1:natm,1:natm,1:nhl(1))= 0d0
    do ihl1=1,nhl(1)
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
!.....then compute dgs wrt w1
    do ihl0=nhl(0),1,-1
      do ihl1=nhl(1),1,-1
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

    gs(1:nvars)= gs(1:nvars) +dgs(1:nvars)
    return
  end subroutine grad1
!=======================================================================
  subroutine grad2(smpl,sds,gs)
    use variables
    implicit none
    type(mdsys),intent(inout):: smpl
    type(smpldata),intent(inout):: sds
    real(8),intent(inout):: gs(nvars)
    integer:: iv,ihl0,ihl1,ihl2,ia,ja,natm
    real(8):: ediff,tmp,tmp1,tmp2,h1,h2,w1,w2,w3,dn3i,dh1,dh2,t1,t2,t3&
         ,ddh1,ddh2,dh1gsf,fscale,eref,swgt
    real(8),save,allocatable:: dgs(:),w1dg(:,:,:,:),w2sw1dg(:,:,:,:)

    if( .not. allocated(dgs) ) then
      allocate( dgs(nvars),w1dg(3,maxna,maxna,nhl(1)) &
           ,w2sw1dg(3,maxna,maxna,nhl(2)) )
    endif

    natm= smpl%natm
    eref= smpl%eref
    swgt= 1d0
    if( lswgt ) then
      swgt= exp(-eref/natm*swbeta)
    endif
    ediff= (smpl%epot -eref)*2 /natm/natm *swgt
    gs(1:nvars)= 0d0
    iv= nhl(0)*nhl(1) +nhl(1)*nhl(2) +nhl(2)

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
      do ihl2=nhl(2),1,-1
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
    do ihl0=nhl(0),1,-1
      do ihl1=nhl(1),1,-1
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

    if( .not. lfmatch ) return
    dgs(1:nvars)= 0d0
    fdiff(1:3,1:natm)= (smpl%fa(1:3,1:natm) &
         -smpl%fref(1:3,1:natm))
    dn3i= 1d0/(3*natm)
    fscale= 1d0
    if( lfscale ) fscale= fscl
    fdiff(1:3,1:natm)= fdiff(1:3,1:natm) *2 *dn3i *fscale *swgt
    iv= nhl(0)*nhl(1) +nhl(1)*nhl(2) +nhl(2)
!.....make w1dg
    w1dg(1:3,1:natm,1:natm,1:nhl(1))= 0d0
    do ihl1=1,nhl(1)
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
!.....make w2sw1dg
    w2sw1dg(1:3,1:natm,1:natm,1:nhl(2))= 0d0
    do ihl2=1,nhl(2)
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
      do ihl2=nhl(2),1,-1
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
    do ihl0=nhl(0),1,-1
      do ihl1=nhl(1),1,-1
        tmp= 0d0
        do ihl2=1,nhl(2)
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
        allocate(wgt11(nhl(0),nhl(1)),wgt12(nhl(1)))
      endif
      iv= 0
      do ihl0=1,nhl(0)
        do ihl1=1,nhl(1)
          iv= iv +1
          wgt11(ihl0,ihl1)= vars(iv)
        enddo
      enddo
      do ihl1=1,nhl(1)
        iv= iv+1
        wgt12(ihl1)= vars(iv)
      enddo
    else if( nl.eq.2 ) then
      if( .not. allocated(wgt21) ) then
        allocate(wgt21(nhl(0),nhl(1)),wgt22(nhl(1),nhl(2)) &
             ,wgt23(nhl(2)))
      endif
      iv= 0
      do ihl0=1,nhl(0)
        do ihl1=1,nhl(1)
          iv=iv+1
          wgt21(ihl0,ihl1)= vars(iv)
        enddo
      enddo
      do ihl1=1,nhl(1)
        do ihl2=1,nhl(2)
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
    use variables
    use parallel
    implicit none

    integer:: itmp,ismpl,natm,ia,ihl0,ja
    character*5:: cdir

    do ismpl=isid0,isid1
      natm= samples(ismpl)%natm
      cdir= samples(ismpl)%cdirname
      !.....gsf
      open(21,file=trim(cmaindir)//'/'//cdir//'/smd/out.NN.gsf'&
           ,status='old')
      read(21,*) itmp
      do ia=1,natm
        do ihl0=1,nhl(0)
          read(21,*) itmp, itmp, sds(ismpl)%gsf(ia,ihl0)
        enddo
      enddo
      close(21)
      !.....dgsf
      open(22,file=trim(cmaindir)//'/'//cdir//'/smd/out.NN.dgsf'&
           ,status='old')
      do ia=1,natm
        do ihl0=1,nhl(0)
          do ja=1,natm
            read(22,*) itmp,itmp,itmp ,sds(ismpl)%dgsf(1:3,ja,ia,ihl0)
          enddo
        enddo
      enddo
      close(22)
    enddo


    if(myid.eq.0) print *, 'get_bases done.'
  end subroutine get_bases
!=======================================================================
  subroutine standardize_max()
!
!  Standardize of inputs is requied when you use lasso or ridge.
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
      do ihl1=1,nhl(1)
        iv=iv+1
        vars(iv)= vars(iv)*gmax(ihl0)
      enddo
    enddo

    lstandard= .true.
    deallocate(gmaxl,gminl)
    if(myid.eq.0) print *,'standardize done.'
  end subroutine standardize_max
!=======================================================================
  subroutine standardize_var()
!
!  Standardize of inputs is requied when you use lasso or ridge.
!
    use variables, only: nsmpl,samples,nvars,nalist,vars
    use parallel
    implicit none
    integer:: ismpl,ia,natm,ihl0,ihl1,iv
    real(8),allocatable:: gmaxl(:),gminl(:),gmeanl(:),gmean(:)
    integer,allocatable:: nsum(:),nsuml(:)

    allocate(gmax(nhl(0)),gmin(nhl(0)) &
         ,gmaxl(nhl(0)),gminl(nhl(0)) &
         ,gmean(nhl(0)),gmeanl(nhl(0)),nsum(nhl(0)),nsuml(nhl(0)))

    !.....compute mean value
    gmeanl(1:nhl(0))= 0d0
    nsuml(1:nhl(0))= 0
    do ismpl=isid0,isid1
      natm= samples(ismpl)%natm
      !.....sum up gsf
      do ihl0=1,nhl(0)
        do ia=1,natm
          gmeanl(ihl0)= gmeanl(ihl0) +sds(ismpl)%gsf(ia,ihl0)
          nsuml(ihl0)= nsuml(ihl0) +1
        enddo
      enddo
    enddo
    gmean(1:nhl(0))= 0d0
    nsum(1:nhl(0))= 0
    call mpi_allreduce(gmeanl,gmean,nhl(0),mpi_double_precision &
         ,mpi_sum,mpi_world,ierr)
    call mpi_allreduce(nsuml,nsum,nhl(0),mpi_integer &
         ,mpi_sum,mpi_world,ierr)
    gmean(1:nhl(0))= gmean(1:nhl(0))/nsum(1:nhl(0))

    !.....compute variance
    gmaxl(1:nhl(0))= 0d0
    do ismpl=isid0,isid1
      natm= samples(ismpl)%natm
      !.....sum up gsf
      do ihl0=1,nhl(0)
        do ia=1,natm
          gmaxl(ihl0)= gmaxl(ihl0) +(gmean(ihl0)-sds(ismpl)%gsf(ia,ihl0))&
               *(gmean(ihl0)-sds(ismpl)%gsf(ia,ihl0))
        enddo
      enddo
    enddo
    gmax(1:nhl(0))= 0d0
    call mpi_allreduce(gmaxl,gmax,nhl(0),mpi_double_precision &
         ,mpi_sum,mpi_world,ierr)
    !.....get standard deviation
    do ihl0=1,nhl(0)
      gmax(ihl0)= sqrt(gmax(ihl0)/nsum(ihl0))
      if( gmax(ihl0).lt.1d-14 ) gmax(ihl0)=1d0
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
      do ihl1=1,nhl(1)
        iv=iv+1
        vars(iv)= vars(iv)*gmax(ihl0)
      enddo
    enddo

    lstandard= .true.
    deallocate(gmaxl,gminl,gmeanl,gmean,nsuml,nsum)
    if(myid.eq.0) print *,'standardize done.'
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

    if( .not. lstandard ) then
      if(myid.eq.0) print *,'NN_restore_standard not needed.'
      return
    endif

    iv= 0
    do ihl0=1,nhl(0)
      do ihl1=1,nhl(1)
        iv=iv+1
        vars(iv)= vars(iv)/gmax(ihl0)
      enddo
    enddo
    if(myid.eq.0) print *,'NN_restore_standard done.'

  end subroutine NN_restore_standard
!=======================================================================
  subroutine NN_analyze()
!
!  Get which input nodes are more/less important.
!
    use variables
    use parallel
    implicit none
    character(len=14),parameter:: cfname= 'out.NN_analyze'
    character(len=14),parameter:: cfsum = 'out.NN_summary'
    integer,parameter:: ionum=  30
    integer,allocatable:: icmb2(:,:),icmb3(:,:,:),itype(:),nctype(:)
    real(8),allocatable:: sumv(:),cnst(:,:),sumvv(:)
    integer:: i,j,k,l,i2,i3,isf,iv,ic,ihl0,ihl1,itmp

    if( myid.eq.0 ) then
!.....read in.comb.NN file
      allocate(icmb2(nsp,nsp),icmb3(nsp,nsp,nsp))
      open(ionum,file=trim(cmaindir)//'/'//trim(cmbfname),status='old')
      do i2=1,ncmb2
        read(ionum,*) i,j,icmb2(i,j)
        icmb2(j,i)= icmb2(i,j)
      enddo
      do i3=1,ncmb3
        read(ionum,*) i,j,k,icmb3(i,j,k)
        icmb3(i,k,j)= icmb3(i,j,k)
      enddo
      close(ionum)
!.....read in.const.NN file
      allocate(itype(nsfc),cnst(2,nsfc),nctype(200))
      nctype(1)= 2   ! Gaussian
      nctype(2)= 1   ! cosine
      nctype(3)= 1   ! polynomial
      nctype(101)= 1 ! angular
      open(ionum,file=trim(cmaindir)//'/'//trim(ccfname),status='old')
      read(ionum,*) itmp
      do isf=1,nsfc
        read(ionum,*) itype(isf),(cnst(j,isf),j=1,nctype(itype(isf)))
      enddo
      close(ionum)

      open(ionum+1,file=cfname,status='replace')
      iv=0
      allocate(sumv(nhl(0)),sumvv(nsfc))
      do ihl0=1,nhl(0)
        sumv(ihl0)= 0d0
        do ihl1=1,nhl(1)
          iv=iv+1
          sumv(ihl0)=sumv(ihl0) +abs(vars(iv))
        enddo
      enddo
!.....about 2body terms
      do ihl0=1,nsf2*ncmb2
        isf= mod(ihl0-1,nsf2)+1
        ic = (ihl0-1)/nsf2 +1
        do i=1,nsp
          do j=1,nsp
            if(ic.eq.icmb2(i,j)) goto 10
          enddo
        enddo
10      continue
        write(ionum+1,'(f24.14,2x,i1,"-",i1,":",i5,1es12.4,2i8)') &
             sumv(ihl0),i,j,itype(isf),cnst(1,isf),ihl0,ic
      enddo
!.....about 3body terms
      do ihl0=nsf2*ncmb2+1,nsf2*ncmb2+nsf3*ncmb3
        isf= nsf2+mod(ihl0-nsf2*ncmb2-1,nsf3)+1
        ic = (ihl0-nsf2*ncmb2-1)/nsf3 +1
        do i=1,nsp
          do j=1,nsp
            do k=1,nsp
              if(ic.eq.icmb3(i,j,k)) goto 20
            enddo
          enddo
        enddo
20      continue
        write(ionum+1,'(f24.14,2x,i1,"-",i1,"-",i1,":",i5,1es12.4,2i8)') &
             sumv(ihl0),i,j,k,itype(isf),cnst(1,isf),ihl0,ic
      enddo
      close(ionum+1)

!.....summary
      open(ionum+2,file=cfsum,status='replace')
      sumvv(1:nsfc)= 0d0
      do ihl0=1,nsf2*ncmb2
        isf= mod(ihl0-1,nsf2)+1
        ic = (ihl0-1)/nsf2 +1
        sumvv(isf)= sumvv(isf) +sumv(ihl0)
      enddo
      do ihl0=nsf2*ncmb2+1,nsf2*ncmb2+nsf3*ncmb3
        isf= nsf2+mod(ihl0-nsf2*ncmb2-1,nsf3)+1
        ic = (ihl0-nsf2*ncmb2-1)/nsf3 +1
        sumvv(isf)= sumvv(isf) +sumv(ihl0)
      enddo
      do isf=1,nsfc
        if( isf.le.nsf2 ) then
          if( itype(isf).eq.1 ) then
            write(ionum+2,'(f24.14,2x,"2body:",i5,2es12.4)') &
                 sumvv(isf),itype(isf),cnst(1:2,isf)
          else if( itype(isf).eq.2 ) then
            write(ionum+2,'(f24.14,2x,"2body:",i5,1es12.4)') &
                 sumvv(isf),itype(isf),cnst(1,isf)
          endif
        else
          write(ionum+2,'(f24.14,2x,"3body:",i5,1es12.4)') &
               sumvv(isf),itype(isf),cnst(1,isf)
        endif
      enddo
      close(ionum+2)
      deallocate(icmb2,icmb3,sumv,sumvv)
    endif

    call mpi_barrier(mpi_world,ierr)
    return
  end subroutine NN_analyze
!=======================================================================
end module NN
