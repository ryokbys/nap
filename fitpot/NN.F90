module NN
!-----------------------------------------------------------------------
!                        Time-stamp: <2015-02-26 17:13:21 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
!.....parameter file name
  character(128),parameter:: cpfname= 'in.params.NN'
  character(128),parameter:: ccfname='in.const.NN'
  character(128),parameter:: cmbfname='in.comb.NN'
  integer,parameter:: maxnl= 2
  integer:: nl,nsp,nsfc,nsf2,nsf3,ncmb2,ncmb3
  integer:: nhl(0:maxnl+1)
  integer,save,allocatable:: nwgt(:)
  real(8),save,allocatable:: wgt11(:,:),wgt12(:)
  real(8),save,allocatable:: wgt21(:,:),wgt22(:,:),wgt23(:)

  type smpldata
    real(8),allocatable:: gsf(:,:),dgsf(:,:,:,:)
    real(8),allocatable:: hl1(:,:),hl2(:,:)
    real(8),allocatable:: ams1(:,:,:,:),bms1(:,:,:,:)
  end type smpldata

  type(smpldata),save,allocatable:: sds(:)

  integer:: maxna
  real(8),save,allocatable:: fdiff(:,:)

contains
!=======================================================================
  subroutine NN_init()
    use variables
    implicit none 
    integer:: itmp,i,nw,natm,ismpl

    timef= 0.0
    timeg= 0.0

    !.....read in.const.NN to get nl,nsp,nhl(:)
    open(20,file=trim(cmaindir)//'/'//trim(ccfname),status='old')
    read(20,*) nl,nsp,nhl(0:nl)
    nsf2= 0
    nsf3= 0
    do while(.true.)
      read(20,*,end=10) itmp
      if( itmp.eq.1 ) then
        nsf2=nsf2+1
      else if( itmp.eq.2 ) then
        nsf3=nsf3+1
      endif
    enddo
10  close(20)
    nsfc= nhl(0)
    nhl(nl+1)= 1

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
!!$    print *, 'nvars,nw=',nvars,nw

    allocate(sds(nsmpl))
    do ismpl=1,nsmpl
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
    do ismpl=1,nsmpl
      if( maxna.lt.samples(ismpl)%natm )  &
           maxna= samples(ismpl)%natm
    enddo
    allocate(fdiff(3,maxna))

    print *, 'NN_init done.'
  end subroutine NN_init
!=======================================================================
  subroutine NN_get_f(fval)
    use variables
    implicit none
    real(8),intent(out):: fval
    integer:: ismpl,natm,ia,ixyz
    real(4):: t0,t1
    real(8):: dn3i,ediff

    call cpu_time(t0)
    call vars2wgts(nvars,vars)
    
    if( nprcs.eq.1 ) then
      do ismpl=1,nsmpl
        if( nl.eq.1 ) then
          call calc_ef1(ismpl)
        else if( nl.eq.2 ) then
          call calc_ef2(ismpl)
        endif
      enddo
    else
      print *,'[nprcs.ne.1] is not implemented yet in NN_get_f.'
      stop
    endif

    fval= 0d0
    do ismpl=1,nsmpl
      natm= samples(ismpl)%natm
      ediff= (samples(ismpl)%epot -samples(ismpl)%eref)
      ediff= ediff*ediff /natm
      fval= fval +ediff
      if( .not. lfmatch ) cycle
      fdiff(1:3,1:natm)= (samples(ismpl)%fa(1:3,1:natm) &
           -samples(ismpl)%fref(1:3,1:natm))
      dn3i= 1d0 /(3*natm)
      fscl= 1d0
      !.....force-scale makes force contribution same order to energy
      if( lfscale ) fscl= 1d0/(3*natm)
      fdiff(1:3,1:natm)= fdiff(1:3,1:natm)*fdiff(1:3,1:natm) &
           *dn3i *fscl
      do ia=1,natm
        do ixyz=1,3
          fval= fval +fdiff(ixyz,ia)
        enddo
      enddo
    enddo
    
    call cpu_time(t1)
!!$    write(6,'(a,f15.3)')  '>>> time NN_get_f =',t1-t0
    timef= timef +t1-t0
    return
  end subroutine NN_get_f
!=======================================================================
  subroutine calc_ef1(ismpl)
    use variables
    implicit none
    integer,intent(in):: ismpl
    integer:: natm,ia,ja,ihl0,ihl1
    real(8):: tmp,w1,w2,h,dh,ddhg,t
    
    natm= samples(ismpl)%natm
    sds(ismpl)%hl1(1:natm,1:nhl(1))= 0d0
    samples(ismpl)%epot =0d0
    !.....energy
    do ia=1,natm
      do ihl1=1,nhl(1)
        tmp= 0d0
        do ihl0=1,nhl(0)
          tmp= tmp +wgt11(ihl0,ihl1) *sds(ismpl)%gsf(ia,ihl0)
        enddo
        sds(ismpl)%hl1(ia,ihl1)= sigmoid(tmp)
      enddo
      do ihl1=1,nhl(1)
        samples(ismpl)%epot= samples(ismpl)%epot &
             +wgt12(ihl1)*(sds(ismpl)%hl1(ia,ihl1)-0.5d0)
      enddo
    enddo

    !.....forces
    if( .not.lfmatch ) return
    if( .not. allocated(sds(ismpl)%ams1) ) then
      allocate(sds(ismpl)%ams1(3,natm,nhl(1),nhl(0)), &
           sds(ismpl)%bms1(3,natm,nhl(1),nhl(0)))
    endif
    samples(ismpl)%fa(1:3,1:natm)= 0d0
    sds(ismpl)%ams1(1:3,1:natm,1:nhl(1),1:nhl(0))= 0d0
    sds(ismpl)%bms1(1:3,1:natm,1:nhl(1),1:nhl(0))= 0d0
    do ihl1=1,nhl(1)
      w2= wgt12(ihl1)
      do ihl0=1,nhl(0)
        w1= wgt11(ihl0,ihl1)
        do ja=1,natm
          h= sds(ismpl)%hl1(ja,ihl1)
          dh= h*(1d0-h)
          ddhg= w1 *dh*(1d0-2d0*h) *sds(ismpl)%gsf(ja,ihl0)
          t= w1*w2*dh
          do ia=1,natm
            samples(ismpl)%fa(1:3,ia)= samples(ismpl)%fa(1:3,ia) &
                 -t *sds(ismpl)%dgsf(1:3,ia,ja,ihl0)
            sds(ismpl)%ams1(1:3,ia,ihl1,ihl0)= &
                 sds(ismpl)%ams1(1:3,ia,ihl1,ihl0) &
                 +dh *sds(ismpl)%dgsf(1:3,ia,ja,ihl0)
            sds(ismpl)%bms1(1:3,ia,ihl1,ihl0)= &
                 sds(ismpl)%bms1(1:3,ia,ihl1,ihl0) &
                 +ddhg *sds(ismpl)%dgsf(1:3,ia,ja,ihl0)
          enddo
        enddo
      enddo
    enddo
    
  end subroutine calc_ef1
!=======================================================================
  subroutine calc_ef2(ismpl)
    implicit none
    integer,intent(in):: ismpl
    
  end subroutine calc_ef2
!=======================================================================
  subroutine NN_get_g(gval)
    use variables
    implicit none
    real(8),intent(out):: gval(nvars)
    integer:: ismpl,i
    real(4):: t0,t1
    real(8),save,allocatable:: gs(:)
    real(8):: gmax,vmax

    if( .not.allocated(gs) ) allocate(gs(nvars))

    call cpu_time(t0)

    gval(1:nvars)= 0d0

    if( nprcs.eq.1 ) then
      do ismpl=1,nsmpl
        if( nl.eq.1 ) then
          call grad1(ismpl,gs)
        else if( nl.eq.2 ) then
          call grad2(ismpl,gs)
        endif
        gval(1:nvars)= gval(1:nvars) +gs(1:nvars)
      enddo
    else
      print *,'[nprcs.ne.1] is not implemented yet in NN_get_g.'
      stop
    endif

    if( lgscale ) then
      gmax= 0d0
      vmax= 0d0
      do i=1,nvars
        vmax= max(vmax,abs(vars(i)))
        gmax= max(gmax,abs(gval(i)))
      enddo
      gval(1:nvars)= gval(1:nvars)/gmax *gscl*vmax
    endif

    call cpu_time(t1)
!!$    write(6,'(a,f15.3)')  '>>> time NN_get_g =',t1-t0
    timeg= timeg +t1-t0
    return
  end subroutine NN_get_g
!=======================================================================
  subroutine grad1(ismpl,gs)
    use variables
    implicit none
    integer,intent(in):: ismpl
    real(8),intent(inout):: gs(nvars)
    integer:: iv,ihl1,ia,ihl0,natm
    real(8):: ediff,tmp,h1,w1,w2,dgs(nvars),ab(3)

    natm= samples(ismpl)%natm
    ediff= (samples(ismpl)%epot -samples(ismpl)%eref) /natm

    gs(1:nvars)= 0d0
    iv= nhl(0)*nhl(1) +nhl(1)
    do ihl1=nhl(1),1,-1
      tmp= 0d0
      do ia=1,natm
        h1= sds(ismpl)%hl1(ia,ihl1)
        tmp= tmp +(h1-0.5d0)
      enddo
      gs(iv)= gs(iv) +2d0*ediff*tmp
      iv= iv -1
    enddo
    do ihl0=nhl(0),1,-1
      do ihl1=nhl(1),1,-1
        tmp= 0d0
        w2= wgt12(ihl1)
        do ia=1,natm
          h1= sds(ismpl)%hl1(ia,ihl1)
          tmp= tmp +w2 *h1*(1d0-h1) *sds(ismpl)%gsf(ia,ihl0)
        enddo
        gs(iv)= gs(iv) +2d0*ediff*tmp
        iv= iv -1
      enddo
    enddo

    if( .not. lfmatch ) return
    dgs(1:nvars)= 0d0
    fdiff(1:3,1:natm)= (samples(ismpl)%fa(1:3,1:natm) &
         -samples(ismpl)%fref(1:3,1:natm))*2/natm/3
    iv= nhl(0)*nhl(1) +nhl(1)
    do ihl1=nhl(1),1,-1
      tmp= 0d0
      do ihl0=1,nhl(0)
        w1= wgt11(ihl0,ihl1)
        do ia=1,natm
          tmp = tmp +w1 *( &
               fdiff(1,ia) *sds(ismpl)%ams1(1,ia,ihl1,ihl0) &
               +fdiff(2,ia)*sds(ismpl)%ams1(2,ia,ihl1,ihl0) &
               +fdiff(3,ia)*sds(ismpl)%ams1(3,ia,ihl1,ihl0) )
        enddo
      enddo
      dgs(iv)= dgs(iv) -tmp
      iv= iv -1
    enddo
    do ihl0=nhl(0),1,-1
      do ihl1=nhl(1),1,-1
        tmp= 0d0
        w2= wgt12(ihl1)
        do ia=1,natm
          ab(1:3)= sds(ismpl)%ams1(1:3,ia,ihl1,ihl0) &
               +sds(ismpl)%bms1(1:3,ia,ihl1,ihl0)
          tmp= tmp +w2*( &
               fdiff(1,ia) *ab(1) &
               +fdiff(2,ia)*ab(2) &
               +fdiff(3,ia)*ab(3) )
        enddo
        dgs(iv)= dgs(iv) -tmp
        iv= iv -1
      enddo
    enddo

    gs(1:nvars)= gs(1:nvars) +dgs(1:nvars)
    return
  end subroutine grad1
!=======================================================================
  subroutine grad2(ismpl,gs)
    use variables
    implicit none
    integer,intent(in):: ismpl
    real(8),intent(inout):: gs(nvars)
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
    implicit none

    integer:: itmp,ismpl,natm,ia,ihl0,ja
    character*5:: cdir
    
    do ismpl=1,nsmpl
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
    print *, 'get_bases done.'
  end subroutine get_bases
end module NN
