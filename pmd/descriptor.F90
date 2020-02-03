module descriptor
!=======================================================================
! Descriptor module
!=======================================================================
  use pmdio, only: csp2isp, nspmax
  implicit none
  include 'params_unit.h'
  save
!!$ putting mpif.h inclusion here could cause some conflicts
!!$  include "mpif.h"
  character(len=128):: paramsdir = '.'
  
  character(128),parameter:: cpfname = 'in.params.desc'
  integer,parameter:: ionum = 51
  logical:: lprmset_desc = .false.

!!$  real(8),parameter:: pi = 3.14159265358979d0
!!$  integer,parameter:: msp = 9  ! hard-coded max-num of species

  type desc
    integer:: itype
    real(8):: rcut,rcut2
    character(len=3):: cspi,cspj,cspk
    integer:: isp,jsp
    integer:: ksp = -1
    integer:: nprm
    real(8),allocatable:: prms(:)
  end type desc
  type(desc),allocatable:: descs(:)
!.....List of isf's for each pair (ilsf2) and angle (ilsf3)
  integer,allocatable:: ilsf2(:,:,:),ilsf3(:,:,:,:)

!.....Memory used in byte
  integer:: mem
!.....Time consumed
  real(8):: time

  integer:: nsf,nsf2,nsf3,nsff
!!$  integer,allocatable:: itype(:)
!!$  real(8),allocatable:: cnst(:,:),rcs(:),rcs2(:)
!.....symmetry function values and their derivatives
  real(8),allocatable:: gsf(:,:),dgsf(:,:,:,:)
  real(8),allocatable:: gsfi(:),dgsfi(:,:,:)
!.....symmetry function IDs for each pair
  integer(2),allocatable:: igsf(:,:,:),igsfi(:,:)
  logical:: lupdate_gsf = .true.
!.....start and end points of symmetry function IDs for each species pair
  integer,allocatable:: iaddr2(:,:,:),iaddr3(:,:,:,:)
!.....function types and num of constatns for types
  integer,parameter:: max_ncnst = 2
  integer:: ncnst_type(200)
  integer:: ncomb_type(200)
  real(8):: cnst(max_ncnst)

!.....Maximum cutoff radius
  real(8):: rcmax,rcmax2,rc3max

  integer:: nalmax,nnlmax,nal,nnl

!.....Chebyshev
  logical:: lcheby = .false.
  real(8),allocatable:: ts_cheby(:),dts_cheby(:),wgtsp(:)

!.....For group LASSO/FS in minimization
  integer:: ngl
  real(8),allocatable:: glval(:)
  integer,allocatable:: mskgfs(:),msktmp(:),iglid(:)

!.....Whether or not called from fitpot [default: .false.]
  logical:: lfitpot = .false.

contains
  subroutine init_desc()
!
!  Initialize descriptor module
!
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
    ncomb_type(1:100) = 2    ! pair
    ncomb_type(101:200) = 3  ! triplet

    time = 0d0
    mem = 0

  end subroutine init_desc
!=======================================================================
  subroutine make_gsf_arrays(l1st,namax,natm,tag,nnmax,lspr &
       ,myid,mpi_world,iprint)
!
!  Make or update gsf arrays
!
    include "mpif.h"
    logical,intent(in):: l1st
    integer,intent(in):: namax,natm,nnmax,lspr(0:nnmax,namax),iprint
    integer,intent(in):: myid,mpi_world
    real(8),intent(in):: tag(namax)

    integer:: i,nnltmp,ierr,time0
    logical,save:: lrealloc = .false.

    if( .not. lupdate_gsf ) return

    time0 = mpi_wtime()

    if( l1st ) then
!  To reduce the memory usage, compute num of atoms and num of neighbors,
!  and add some margin for those numbers because they can change during
!  the simulation.
      nal = int(natm*1.1)
      nnltmp = 0
      do i=1,natm
        nnltmp = max(nnltmp,lspr(0,i))
      enddo
      nnl = int(nnltmp*1.1)
      if( nal .gt. namax ) then
        write(6,'(a)') ' [Error] nal .gt.namax'
        write(6,'(a,3i10)') '   myid,nal,namax = ',myid,nal,namax
        stop
      endif
      if( nnl.gt.nnmax ) then
        write(6,'(a)') ' [Error] nnl.gt.nnmax'
        write(6,'(a,3i10)') '   myid,nnl,nnmax = ',myid,nnl,nnmax
        stop
      endif
      if( myid.ge.0 ) then
        call mpi_reduce(nal,nalmax,1,mpi_integer,mpi_max,0,mpi_world,ierr)
        call mpi_reduce(nnl,nnlmax,1,mpi_integer,mpi_max,0,mpi_world,ierr)
      else
        nalmax = nal
        nnlmax = nnl
      endif
      if( myid.le.0 .and. iprint.ne.0 ) then
        print *,''
        print *,'make_gsf_arrays @descriptor:'
        write(6,'(a,2i0)') '   Max num of (local atoms *1.1) = ',nalmax
        write(6,'(a,2i0)') '   Max num of (neighbors *1.1)   = ',nnlmax
        write(6,'(a,f10.3,a)') '   gsf size  = ', &
             dble(nsf*nalmax*8)/1000/1000,' MB'
        write(6,'(a,f10.3,a)') '   dgsf size = ', &
             dble(3*nsf*(nnlmax+1)*nalmax*8)/1000/1000,' MB'
        write(6,'(a,f10.3,a)') '   igsf size = ', &
             dble(nsf*(nnlmax+1)*nalmax*2)/1000/1000,' MB'
      endif
      if( allocated(gsf) ) then
        if( size(gsf).lt.nsf*nal .or. size(dgsf).lt.3*nsf*(nnl+1)*nal ) then
          lrealloc = .true.
        else
          lrealloc = .false.
        endif
      else
        lrealloc = .true.
      endif

!.....gsfi,dgsfi,igsfi are independend on number of atoms but on nnlmax
      if( allocated(gsfi) ) deallocate(gsfi,dgsfi,igsfi)
      allocate(gsfi(nsf),dgsfi(3,nsf,0:nnmax),igsfi(nsf,0:nnmax))
      mem = mem +8*size(gsfi) +8*size(dgsfi) +2*size(igsfi)
    endif

!  Since natm and nn can change every step of MD,
!  if natm/nnltmp becomes >nal/nnl, they should be updated and
!  gsf/dgsf as well.
    if( natm.gt.nal ) then
      nal = int(natm*1.1)
      if( nal .gt. namax ) then
        write(6,'(a)') ' [Error] nal .gt.namax'
        write(6,'(a,3i0)') '   myid,nal,namax = ',myid,nal,namax
        stop
      else
        write(6,*) ' Since natm.gt.nal, nal was updated at myid =',myid
      endif
      lrealloc=.true.
    endif

    nnltmp = 0
    do i=1,natm
      nnltmp = max(nnltmp,lspr(0,i))
    enddo
    if( nnltmp.gt.nnl ) then
      nnl = int(nnltmp*1.1)
      if( nnlmax.gt.nnmax ) then
        write(6,'(a)') ' [Error] nnl.gt.nnmax'
        write(6,'(a,3i0)') '   myid,nnl,nnmax = ',myid,nnl,nnmax
        stop
      else
        write(6,*) ' Since nnltmp.gt.nnl, nnl was updated at myid =',myid
      endif
      lrealloc=.true.
    endif

!!$    if( allocated(dgsf).and.lrealloc ) then
!!$      mem = mem -8*size(gsf) -8*size(dgsf) -4*size(igsf)
!!$      deallocate( gsf,dgsf,igsf )
!!$      allocate( gsf(nsf,nal),dgsf(3,nsf,0:nnl,nal) &
!!$           ,igsf(nsf,0:nnl,nal))
!!$      mem = mem +8*size(gsf) +8*size(dgsf) +4*size(igsf)
!!$      lrealloc=.false.
!!$    endif
    if( lrealloc ) then
      if( allocated(gsf) ) then
        mem = mem -8*size(gsf) -8*size(dgsf) -2*size(igsf)
        deallocate( gsf,dgsf,igsf )
      endif
      allocate( gsf(nsf,nal),dgsf(3,nsf,0:nnl,nal) &
           ,igsf(nsf,0:nnl,nal))
      mem = mem +8*size(gsf) +8*size(dgsf) +2*size(igsf)
      lrealloc=.false.
    endif

    time = time +(mpi_wtime() -time0)
    
    return
  end subroutine make_gsf_arrays
!=======================================================================
  subroutine calc_desc(namax,natm,nb,nnmax,h,tag,ra,lspr,rc &
       ,myid,mpi_world,l1st,iprint)

    integer,intent(in):: namax,natm,nb,nnmax,lspr(0:nnmax,namax)
    integer,intent(in):: myid,mpi_world,iprint
    real(8),intent(in):: h(3,3),tag(namax),ra(3,namax),rc
    logical,intent(in):: l1st
    include "mpif.h"

    real(8):: time0

    if( .not.lupdate_gsf ) return
    
    time0 = mpi_wtime()

    if( lcheby ) then
      call calc_desc_cheby(namax,natm,nb,nnmax,h,tag,ra,lspr,rc &
           ,myid,mpi_world,l1st,iprint)
    else ! default
      call calc_desc_default(namax,natm,nb,nnmax,h,tag,ra,lspr,rc &
           ,myid,mpi_world,l1st,iprint)
    endif

    time = time +(mpi_wtime() -time0)
    
    return
  end subroutine calc_desc
!=======================================================================
  subroutine calc_desci(ia,namax,natm,nnmax,h,tag,ra,lspr,rc,iprint)
!
!  Wrapper for descriptor calculation
!
    integer,intent(in):: ia,namax,natm,nnmax,lspr(0:nnmax,namax)
    integer,intent(in):: iprint
    real(8),intent(in):: h(3,3),tag(namax),ra(3,namax),rc

    gsfi(:)= 0d0
    dgsfi(:,:,:)= 0d0
    igsfi(:,:) = 0

!!$    if( iprint.gt.1 .and. ia.eq.1 ) print *,'ia,lupdate_gsf,lfitpot=',ia,lupdate_gsf,lfitpot
    if( lupdate_gsf ) then
      if( lcheby ) then
        call desci_cheby(ia,namax,natm,nnmax,h,tag,ra,lspr,rc,iprint)
      else ! default
        call desci_default(ia,namax,natm,nnmax,h,tag,ra,lspr,rc,iprint)
      endif
      if( lfitpot ) then
        if( .not. allocated(gsf) .or. .not. allocated(dgsf) &
             .or. .not. allocated(igsf) ) then
          print *,'ERROR: either gsf/dgsf/igsf is not allocated, which should not happen'&
               //' when called from fitpot.'
          stop
        endif
        gsf(:,ia) = gsfi(:)
        dgsf(:,:,0:nnl,ia) = dgsfi(:,:,0:nnl)
        igsf(:,0:nnl,ia) = igsfi(:,0:nnl)
      endif
    else ! Not to update gsf by desci_xxx just use gsfs given by fitpot
      if( lfitpot ) then
!!$        if( iprint.gt.1 .and. ia.eq.1 ) print *,'ia,nnl,nal=',ia,nnl,nal
        gsfi(:) = gsf(:,ia)
        dgsfi(:,:,0:nnl) = dgsf(:,:,0:nnl,ia)
        igsfi(:,0:nnl) = igsf(:,0:nnl,ia)
      endif
    endif


  end subroutine calc_desci
!=======================================================================
  subroutine prepare_desci(myid,iprint,rc)
!
!  Preparation for desci calculation
!
    include "mpif.h"
    integer,intent(in):: myid,iprint
    real(8),intent(in):: rc

    integer:: ierr
    
    if( .not. allocated(ts_cheby) ) then
      allocate(ts_cheby(0:max(nsf2,nsf3)),dts_cheby(0:max(nsf2,nsf3)))
      mem = mem +8*size(ts_cheby)*2
    endif

!.....Check the maximumx cutoff and given rc
    if( rc.lt.rcmax ) then
      if( myid.eq.0 ) then
        print *,'ERROR: cutoff radius rc is smaller than rcmax.'
        print *,'  rc,rcmax = ',rc,rcmax
      endif
      call mpi_finalize(ierr)
      stop
    endif
  end subroutine prepare_desci
!=======================================================================
  subroutine calc_desc_default(namax,natm,nb,nnmax,h,tag,ra,lspr,rc &
       ,myid,mpi_world,l1st,iprint)
!
!  Evaluate descriptors (symmetry functions)
!  and their derivatives wrt positions for multi-species system.
!
!  - Cutoff radii are set in each symmetry functions.
!  - If the overlay option is set, use inner and outer cutoff of ZBL potential.
!
!!$    implicit none
    include "mpif.h"
    integer,intent(in):: namax,natm,nb,nnmax,lspr(0:nnmax,namax)
    integer,intent(in):: myid,mpi_world,iprint
    real(8),intent(in):: h(3,3),tag(namax),ra(3,namax),rc
    logical,intent(in):: l1st 

    integer:: isf,ia,jj,ja,kk,ka,is,js,ks,ierr,i,isp,jsp,ksp,ityp,is1,is2 &
         ,ksf,itypp
    real(8):: xi(3),xj(3),xij(3),rij(3),dij,dij2,fcij,eta,rs,texp,driji(3), &
         dfcij,drijj(3),dgdr,xk(3),xik(3),rik(3),dik,fcik,dfcik, &
         driki(3),drikk(3),almbd,spijk,cs,t1,t2,dgdij,dgdik,dgcs, &
         dcsdj(3),dcsdk(3),dcsdi(3),tcos,tpoly,a1,a2,tmorse,dik2,tmp,dtmp, &
         xjk(3),rjk(3),djk,djk2,fcjk,dfcjk,drjkj(3),drjkk(3),dgdjk, &
         ri,ro,xs,z,dz,an,gijk,rcut,rcut2,rcutp,ttmp
    type(desc):: desci
    real(8):: texpij,texpik,eta3,zang,twozeta
!!$    real(8),save:: time2, time3

    if( l1st ) then
!!$      time2 = 0d0
!!$      time3 = 0d0
      if( iprint.gt.1 .and. myid.eq.0 ) then
        print *,'calc_desc_default...'
        print '(a,4es12.4)','  rc,rcmax,rcmax2,rc3max=',rc,rcmax,rcmax2,rc3max
      endif
!.....Check all the rcs and compare them with rc
      if( rc.lt.rcmax ) then
        if( myid.eq.0 ) then
          print *,'ERROR: cutoff radius rc is smaller than rcmax.'
          print *,'  rc,rcmax = ',rc,rcmax
        endif
        call mpi_finalize(ierr)
        stop
      endif
    endif

    gsf(1:nsf,1:nal)= 0d0
    dgsf(1:3,1:nsf,0:nnl,1:nal)= 0d0
    igsf(1:nsf,0:nnl,1:nal) = 0
    itypp = -1
    do ia=1,natm
      xi(1:3)= ra(1:3,ia)
      is= int(tag(ia))
      do jj=1,lspr(0,ia)
        ja= lspr(jj,ia)
        if( ja.eq.ia ) cycle
        xj(1:3)= ra(1:3,ja)
        xij(1:3)= xj(1:3) -xi(1:3)
        rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij2 = rij(1)**2 +rij(2)**2 +rij(3)**2
        if( dij2.ge.rcmax2 ) cycle
        dij = dsqrt(dij2)
        js= int(tag(ja))
        driji(1:3)= -rij(1:3)/dij
        drijj(1:3)= -driji(1:3)
        is1 = min(is,js)
        is2 = max(is,js)
!!$        ttmp = mpi_wtime()
        do ksf=1,ilsf2(0,is1,is2)
          isf = ilsf2(ksf,is1,is2)
          desci = descs(isf)
          ityp = desci%itype
          rcut = desci%rcut
          if( dij.ge.rcut ) cycle
          if( ityp.eq.1 ) then ! Gaussian
            call get_fc_dfc(dij,is,js,rcut,fcij,dfcij)
            eta= desci%prms(1)
            rs = desci%prms(2)
!.....function value
            texp= exp(-eta*(dij-rs)**2)
            dgdr= -2d0*eta*(dij-rs)*texp*fcij +texp*dfcij
            tmp = texp*fcij
            gsf(isf,ia)= gsf(isf,ia) +tmp
!.....derivative
! dgsf(ixyz,isf,jj,ia): derivative of isf-th basis of atom-ia
! by ixyz coordinate of atom-jj. jj=0 means derivative by atom-ia.
            dgsf(1:3,isf,0,ia)= dgsf(1:3,isf,0,ia) +driji(1:3)*dgdr
            dgsf(1:3,isf,jj,ia)= dgsf(1:3,isf,jj,ia) +drijj(1:3)*dgdr
            igsf(isf,0,ia) = 1
            igsf(isf,jj,ia) = 1
          else if( ityp.eq.2 ) then ! cosine
            call get_fc_dfc(dij,is,js,rcut,fcij,dfcij)
            a1 = desci%prms(1)
!.....func value
            tcos= 0.5d0*(1d0+cos(dij*a1))
            dgdr= -0.5d0*a1*sin(dij*a1)*fcij +tcos*dfcij
            tmp = tcos*fcij
            gsf(isf,ia)= gsf(isf,ia) +tmp
!.....derivative
            dgsf(1:3,isf,0,ia)= dgsf(1:3,isf,0,ia) +driji(1:3)*dgdr
            dgsf(1:3,isf,jj,ia)= dgsf(1:3,isf,jj,ia) +drijj(1:3)*dgdr
            igsf(isf,0,ia) = 1
            igsf(isf,jj,ia) = 1
          else if( ityp.eq.3 ) then ! polynomial
            call get_fc_dfc(dij,is,js,rcut,fcij,dfcij)
            a1= desci%prms(1)
!.....func value
            tpoly= 1d0*dij**(-a1)
            dgdr= -a1*dij**(-a1-1d0)*fcij +tpoly*dfcij
            tmp = tpoly*fcij
            gsf(isf,ia)= gsf(isf,ia) +tmp
!.....derivative
            dgsf(1:3,isf,0,ia)= dgsf(1:3,isf,0,ia) +driji(1:3)*dgdr
            dgsf(1:3,isf,jj,ia)= dgsf(1:3,isf,jj,ia) +drijj(1:3)*dgdr
            igsf(isf,0,ia) = 1
            igsf(isf,jj,ia) = 1
          else if( ityp.eq.4 ) then ! Morse-type
            call get_fc_dfc(dij,is,js,rcut,fcij,dfcij)
            a1= desci%prms(1)
            a2= desci%prms(2)
!.....func value
            texp= exp(-a1*(dij-a2))
            tmorse= ((1d0-texp)**2 -1d0)
            dgdr= 2d0*a1*(1d0-texp)*texp*fcij +tmorse*dfcij
            tmp = tmorse*fcij
            gsf(isf,ia)= gsf(isf,ia) +tmp
!.....derivative
            dgsf(1:3,isf,0,ia)= dgsf(1:3,isf,0,ia) +driji(1:3)*dgdr
            dgsf(1:3,isf,jj,ia)= dgsf(1:3,isf,jj,ia) +drijj(1:3)*dgdr
            igsf(isf,0,ia) = 1
            igsf(isf,jj,ia) = 1
          endif
        enddo  ! isf=1,nsf
!!$        time2 = time2 +(mpi_wtime() -ttmp)

!.....3-body forms
        if( dij.gt.rc3max ) cycle
!!$        ttmp = mpi_wtime()
        do kk=1,lspr(0,ia)
          ka= lspr(kk,ia)
          ks= int(tag(ka))
          if( ka.eq.ia .or. ka.le.ja ) cycle
          xk(1:3)= ra(1:3,ka)
          xik(1:3)= xk(1:3) -xi(1:3)
          rik(1:3)= h(1:3,1)*xik(1) +h(1:3,2)*xik(2) +h(1:3,3)*xik(3)
          dik2 = rik(1)*rik(1) +rik(2)*rik(2) +rik(3)*rik(3)
          dik= dsqrt(dik2)
!.....Cosine is common for all the angular SFs
          spijk= rij(1)*rik(1) +rij(2)*rik(2) +rij(3)*rik(3)
          cs= spijk/dij/dik
          dcsdj(1:3)= rik(1:3)/dij/dik -rij(1:3)*cs/dij2
          dcsdk(1:3)= rij(1:3)/dij/dik -rik(1:3)*cs/dik2
          dcsdi(1:3)= -dcsdj(1:3) -dcsdk(1:3)
          gijk = 0d0
          is1 = min(js,ks)
          is2 = max(js,ks)
          do ksf=1,ilsf3(0,is,is1,is2)
            isf = ilsf3(ksf,is,is1,is2)
            desci = descs(isf)
            ityp = desci%itype
            rcut = desci%rcut
            rcut2= desci%rcut2
            if( dij.ge.rcut .or. dik.ge.rcut ) cycle
            if( ityp.eq.101 ) then ! RK's original angular SF
!.....fcij's should be computed after rcs is determined
              call get_fc_dfc(dij,is,js,rcut,fcij,dfcij)
              call get_fc_dfc(dik,is,ks,rcut,fcik,dfcik)
              almbd= desci%prms(1)
              t2= (abs(almbd)+1d0)**2
              driki(1:3)= -rik(1:3)/dik
              drikk(1:3)= -driki(1:3)
!.....function value
              t1= (almbd +cs)**2
              eta3 = 0.5d0 /rcut**2
              texp = exp(-eta3*(dij2+dik2))
              tmp = t1/t2 *texp
              gsf(isf,ia)= gsf(isf,ia) +tmp*fcij*fcik
              gijk = gijk +tmp*fcij*fcik
!.....derivative
              dgdij= dfcij *fcik *tmp &
                   +tmp *(-2d0*eta3*dij) *fcij*fcik 
              dgdik= fcij *dfcik *tmp &
                   +tmp *(-2d0*eta3*dik) *fcij*fcik 
              dgsf(1:3,isf,0,ia)= dgsf(1:3,isf,0,ia) &
                   +dgdij*driji(1:3) +dgdik*driki(1:3)
              dgsf(1:3,isf,jj,ia)= dgsf(1:3,isf,jj,ia) +dgdij*drijj(1:3)
              dgsf(1:3,isf,kk,ia)= dgsf(1:3,isf,kk,ia) +dgdik*drikk(1:3)
              dgcs= 2d0*(almbd+cs)/t2 *fcij*fcik *texp 
              dgsf(1:3,isf,0,ia)= dgsf(1:3,isf,0,ia) +dgcs*dcsdi(1:3)
              dgsf(1:3,isf,jj,ia)= dgsf(1:3,isf,jj,ia) +dgcs*dcsdj(1:3)
              dgsf(1:3,isf,kk,ia)= dgsf(1:3,isf,kk,ia) +dgcs*dcsdk(1:3)
              igsf(isf,0,ia) = 1
              igsf(isf,jj,ia) = 1
              igsf(isf,kk,ia) = 1
            else if( ityp.eq.102 ) then  ! Similar to Behler's angular SF that includes fc(rjk)
!.....djk is required for Behler's angular SF (itype(isf)==102)
              xjk(1:3)= xk(1:3)-xj(1:3)
              rjk(1:3)= h(1:3,1)*xjk(1) +h(1:3,2)*xjk(2) +h(1:3,3)*xjk(3)
              djk2= rjk(1)**2 +rjk(2)**2 +rjk(3)**2
              if( djk2.ge.rcut2 ) cycle
              djk= sqrt(djk2)
!.....fcij's should be computed after rcs is determined
              call get_fc_dfc(dij,is,js,rcut,fcij,dfcij)
              call get_fc_dfc(dik,is,ks,rcut,fcik,dfcik)
              call get_fc_dfc(djk,js,ks,rcut,fcjk,dfcjk)
              almbd= desci%prms(1)
              zang = desci%prms(2)
              driki(1:3)= -rik(1:3)/dik
              drikk(1:3)= -driki(1:3)
              drjkj(1:3)= -rjk(1:3)/djk
              drjkk(1:3)= -drjkj(1:3)
!.....function value
              twozeta = 2d0**(-zang)
              t1= (almbd +cs)**zang *twozeta
              eta3 = 0.5d0 /rcut**2
              texp = exp(-eta3*(dij2+dik2+djk2))
              tmp = t1 *texp
!.....This part is different from itype(isf)==101 by the factor fcjk
              gsf(isf,ia)= gsf(isf,ia) +tmp*fcij*fcik *fcjk
!.....derivative
              dgdij= dfcij *fcik*fcjk *tmp &
                   +tmp *(-2d0*eta3*dij) *fcij*fcik*fcjk
              dgdik= fcij *dfcik*fcjk *tmp &
                   +tmp *(-2d0*eta3*dik) *fcij*fcik*fcjk 
              dgdjk= fcij *fcik *dfcjk *tmp &
                   +tmp *(-2d0*eta3*djk) *fcij*fcik*fcjk 
              dgsf(1:3,isf,0,ia)= dgsf(1:3,isf,0,ia) &
                   +dgdij*driji(1:3) +dgdik*driki(1:3)
              dgsf(1:3,isf,jj,ia)= dgsf(1:3,isf,jj,ia) &
                    +dgdij*drijj(1:3)+dgdjk*drjkj(1:3)
              dgsf(1:3,isf,kk,ia)= dgsf(1:3,isf,kk,ia) &
                    +dgdik*drikk(1:3)+dgdjk*drjkk(1:3)
              dgcs= zang*(almbd +cs)**(zang-1d0) *twozeta *fcij*fcik*fcjk *texp
              dgsf(1:3,isf,0,ia)= dgsf(1:3,isf,0,ia) +dgcs*dcsdi(1:3)
              dgsf(1:3,isf,jj,ia)= dgsf(1:3,isf,jj,ia) +dgcs*dcsdj(1:3)
              dgsf(1:3,isf,kk,ia)= dgsf(1:3,isf,kk,ia) +dgcs*dcsdk(1:3)
              igsf(isf,0,ia) = 1
              igsf(isf,jj,ia) = 1
              igsf(isf,kk,ia) = 1
            else if( ityp.eq.103 ) then ! cos(cos(thijk)*n*pi) w/o fc(rjk)
!.....fcij's should be computed after rcs is determined
              call get_fc_dfc(dij,is,js,rcut,fcij,dfcij)
              call get_fc_dfc(dik,is,ks,rcut,fcik,dfcik)
              an = desci%prms(1)
              driki(1:3)= -rik(1:3)/dik
              drikk(1:3)= -driki(1:3)
!.....Function value
              tmp = cos(cs*an*pi)
              gsf(isf,ia)= gsf(isf,ia) +tmp*fcij*fcik
!.....Derivative
              dgdij= tmp*dfcij*fcik
              dgdik= tmp*fcij*dfcik
              dgsf(1:3,isf,0,ia)= dgsf(1:3,isf,0,ia)&
                   +dgdij*driji(1:3) +dgdik*driki(1:3)
              dgsf(1:3,isf,jj,ia)= dgsf(1:3,isf,jj,ia) +dgdij*drijj(1:3)
              dgsf(1:3,isf,kk,ia)= dgsf(1:3,isf,kk,ia) +dgdik*drikk(1:3)
              dgcs= -an*pi*sin(cs*an*pi) *fcij*fcik
              dgsf(1:3,isf,0,ia)= dgsf(1:3,isf,0,ia) +dgcs*dcsdi(1:3)
              dgsf(1:3,isf,jj,ia)= dgsf(1:3,isf,jj,ia) +dgcs*dcsdj(1:3)
              dgsf(1:3,isf,kk,ia)= dgsf(1:3,isf,kk,ia) +dgcs*dcsdk(1:3)
              igsf(isf,0,ia) = 1
              igsf(isf,jj,ia) = 1
              igsf(isf,kk,ia) = 1
            else if( ityp.eq.104 ) then ! sin(cos(thijk)*n*pi) w/o fc(rjk)
!.....fcij's should be computed after rcs is determined
              call get_fc_dfc(dij,is,js,rcut,fcij,dfcij)
              call get_fc_dfc(dik,is,ks,rcut,fcik,dfcik)
              an = desci%prms(1)
              driki(1:3)= -rik(1:3)/dik
              drikk(1:3)= -driki(1:3)
!.....Function value
              tmp = sin(cs*an*pi)
              gsf(isf,ia)= gsf(isf,ia) +tmp*fcij*fcik
!.....Derivative
              dgdij= tmp*dfcij*fcik
              dgdik= tmp*fcij*dfcik
              dgsf(1:3,isf,0,ia)= dgsf(1:3,isf,0,ia)&
                   +dgdij*driji(1:3) +dgdik*driki(1:3)
              dgsf(1:3,isf,jj,ia)= dgsf(1:3,isf,jj,ia) +dgdij*drijj(1:3)
              dgsf(1:3,isf,kk,ia)= dgsf(1:3,isf,kk,ia) +dgdik*drikk(1:3)
              dgcs= an*pi*cos(cs*an*pi) *fcij*fcik
              dgsf(1:3,isf,0,ia)= dgsf(1:3,isf,0,ia) +dgcs*dcsdi(1:3)
              dgsf(1:3,isf,jj,ia)= dgsf(1:3,isf,jj,ia) +dgcs*dcsdj(1:3)
              dgsf(1:3,isf,kk,ia)= dgsf(1:3,isf,kk,ia) +dgcs*dcsdk(1:3)
              igsf(isf,0,ia) = 1
              igsf(isf,jj,ia) = 1
              igsf(isf,kk,ia) = 1
            else if( ityp.eq.105 ) then ! exp(-eta*(cos(thijk)-c)**2) w/o fc(rjk)
!.....fcij's should be computed after rcs is determined
              call get_fc_dfc(dij,is,js,rcut,fcij,dfcij)
              call get_fc_dfc(dik,is,ks,rcut,fcik,dfcik)
              eta= desci%prms(1)
              rs = desci%prms(2)
              driki(1:3)= -rik(1:3)/dik
              drikk(1:3)= -driki(1:3)
!.....Function value
              tmp = exp(-eta*(cs-rs)**2)
              gsf(isf,ia)= gsf(isf,ia) +tmp*fcij*fcik
!.....Derivative
              dgdij= tmp*dfcij*fcik
              dgdik= tmp*fcij*dfcik
              dgsf(1:3,isf,0,ia)= dgsf(1:3,isf,0,ia)&
                   +dgdij*driji(1:3) +dgdik*driki(1:3)
              dgsf(1:3,isf,jj,ia)= dgsf(1:3,isf,jj,ia) +dgdij*drijj(1:3)
              dgsf(1:3,isf,kk,ia)= dgsf(1:3,isf,kk,ia) +dgdik*drikk(1:3)
              dgcs= -2d0*eta*(cs-rs)*tmp *fcij*fcik
              dgsf(1:3,isf,0,ia)= dgsf(1:3,isf,0,ia) +dgcs*dcsdi(1:3)
              dgsf(1:3,isf,jj,ia)= dgsf(1:3,isf,jj,ia) +dgcs*dcsdj(1:3)
              dgsf(1:3,isf,kk,ia)= dgsf(1:3,isf,kk,ia) +dgcs*dcsdk(1:3)
              igsf(isf,0,ia) = 1
              igsf(isf,jj,ia) = 1
              igsf(isf,kk,ia) = 1
            endif

          enddo ! isf=1,...
        enddo ! kk=1,...
!!$        time3 = time3 +(mpi_wtime() -ttmp)
20      continue
      enddo ! jj=1,...
    enddo ! ia=1,...

    return
  end subroutine calc_desc_default
!=======================================================================
  subroutine desci_default(ia,namax,natm,nnmax,h,tag,ra,lspr,rc,iprint)
!
!  Evaluate descriptors (symmetry functions) of an atom-i only
!  and their derivatives wrt positions for multi-species system.
!
!  - Cutoff radii are set in each symmetry functions.
!  - If the overlay option is set, use inner and outer cutoff of ZBL potential.
!
    integer,intent(in):: ia,namax,natm,nnmax,lspr(0:nnmax,namax)
    integer,intent(in):: iprint
    real(8),intent(in):: h(3,3),tag(namax),ra(3,namax),rc

    integer:: isf,jj,ja,kk,ka,is,js,ks,ierr,i,isp,jsp,ksp,ityp,is1,is2 &
         ,ksf,itypp
    real(8):: xi(3),xj(3),xij(3),rij(3),dij,dij2,fcij,eta,rs,texp,driji(3), &
         dfcij,drijj(3),dgdr,xk(3),xik(3),rik(3),dik,fcik,dfcik, &
         driki(3),drikk(3),almbd,spijk,cs,t1,t2,dgdij,dgdik,dgcs, &
         dcsdj(3),dcsdk(3),dcsdi(3),tcos,tpoly,a1,a2,tmorse,dik2,tmp,dtmp, &
         xjk(3),rjk(3),djk,djk2,fcjk,dfcjk,drjkj(3),drjkk(3),dgdjk, &
         ri,ro,xs,z,dz,an,gijk,rcut,rcut2,rcutp,ttmp
    type(desc):: desci
    real(8):: texpij,texpik,eta3,zang,twozeta
!!$    real(8),save:: time2, time3

    xi(1:3)= ra(1:3,ia)
    is= int(tag(ia))
    do jj=1,lspr(0,ia)
      ja= lspr(jj,ia)
      if( ja.eq.ia ) cycle
      xj(1:3)= ra(1:3,ja)
      xij(1:3)= xj(1:3) -xi(1:3)
      rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
      dij2 = rij(1)**2 +rij(2)**2 +rij(3)**2
      if( dij2.ge.rcmax2 ) cycle
      dij = dsqrt(dij2)
      js= int(tag(ja))
      driji(1:3)= -rij(1:3)/dij
      drijj(1:3)= -driji(1:3)
      is1 = min(is,js)
      is2 = max(is,js)
      do ksf=1,ilsf2(0,is1,is2)
        isf = ilsf2(ksf,is1,is2)
        desci = descs(isf)
        ityp = desci%itype
        rcut = desci%rcut
        if( dij.ge.rcut ) cycle
        if( ityp.eq.1 ) then ! Gaussian
          call get_fc_dfc(dij,is,js,rcut,fcij,dfcij)
          eta= desci%prms(1)
          rs = desci%prms(2)
!.....function value
          texp= exp(-eta*(dij-rs)**2)
          dgdr= -2d0*eta*(dij-rs)*texp*fcij +texp*dfcij
          tmp = texp*fcij
          gsfi(isf)= gsfi(isf) +tmp
!.....derivative
! dgsf(ixyz,isf,jj,ia): derivative of isf-th basis of atom-ia
! by ixyz coordinate of atom-jj. jj=0 means derivative by atom-ia.
          dgsfi(1:3,isf,0)= dgsfi(1:3,isf,0) +driji(1:3)*dgdr
          dgsfi(1:3,isf,jj)= dgsfi(1:3,isf,jj) +drijj(1:3)*dgdr
          igsfi(isf,0) = 1
          igsfi(isf,jj) = 1
        else if( ityp.eq.2 ) then ! cosine
          call get_fc_dfc(dij,is,js,rcut,fcij,dfcij)
          a1 = desci%prms(1)
!.....func value
          tcos= 0.5d0*(1d0+cos(dij*a1))
          dgdr= -0.5d0*a1*sin(dij*a1)*fcij +tcos*dfcij
          tmp = tcos*fcij
          gsfi(isf)= gsfi(isf) +tmp
!.....derivative
          dgsfi(1:3,isf,0)= dgsfi(1:3,isf,0) +driji(1:3)*dgdr
          dgsfi(1:3,isf,jj)= dgsfi(1:3,isf,jj) +drijj(1:3)*dgdr
          igsfi(isf,0) = 1
          igsfi(isf,jj) = 1
        else if( ityp.eq.3 ) then ! polynomial
          call get_fc_dfc(dij,is,js,rcut,fcij,dfcij)
          a1= desci%prms(1)
!.....func value
          tpoly= 1d0*dij**(-a1)
          dgdr= -a1*dij**(-a1-1d0)*fcij +tpoly*dfcij
          tmp = tpoly*fcij
          gsfi(isf)= gsfi(isf) +tmp
!.....derivative
          dgsfi(1:3,isf,0)= dgsfi(1:3,isf,0) +driji(1:3)*dgdr
          dgsfi(1:3,isf,jj)= dgsfi(1:3,isf,jj) +drijj(1:3)*dgdr
          igsfi(isf,0) = 1
          igsfi(isf,jj) = 1
        else if( ityp.eq.4 ) then ! Morse-type
          call get_fc_dfc(dij,is,js,rcut,fcij,dfcij)
          a1= desci%prms(1)
          a2= desci%prms(2)
!.....func value
          texp= exp(-a1*(dij-a2))
          tmorse= ((1d0-texp)**2 -1d0)
          dgdr= 2d0*a1*(1d0-texp)*texp*fcij +tmorse*dfcij
          tmp = tmorse*fcij
          gsfi(isf)= gsfi(isf) +tmp
!.....derivative
          dgsfi(1:3,isf,0)= dgsfi(1:3,isf,0) +driji(1:3)*dgdr
          dgsfi(1:3,isf,jj)= dgsfi(1:3,isf,jj) +drijj(1:3)*dgdr
          igsfi(isf,0) = 1
          igsfi(isf,jj) = 1
        endif
      enddo  ! isf=1,nsf

!.....3-body forms
      if( dij.gt.rc3max ) cycle
      do kk=1,lspr(0,ia)
        ka= lspr(kk,ia)
        ks= int(tag(ka))
        if( ka.eq.ia .or. ka.le.ja ) cycle
        xk(1:3)= ra(1:3,ka)
        xik(1:3)= xk(1:3) -xi(1:3)
        rik(1:3)= h(1:3,1)*xik(1) +h(1:3,2)*xik(2) +h(1:3,3)*xik(3)
        dik2 = rik(1)*rik(1) +rik(2)*rik(2) +rik(3)*rik(3)
        dik= dsqrt(dik2)
!.....Cosine is common for all the angular SFs
        spijk= rij(1)*rik(1) +rij(2)*rik(2) +rij(3)*rik(3)
        cs= spijk/dij/dik
        dcsdj(1:3)= rik(1:3)/dij/dik -rij(1:3)*cs/dij2
        dcsdk(1:3)= rij(1:3)/dij/dik -rik(1:3)*cs/dik2
        dcsdi(1:3)= -dcsdj(1:3) -dcsdk(1:3)
        gijk = 0d0
        is1 = min(js,ks)
        is2 = max(js,ks)
        do ksf=1,ilsf3(0,is,is1,is2)
          isf = ilsf3(ksf,is,is1,is2)
          desci = descs(isf)
          ityp = desci%itype
          rcut = desci%rcut
          rcut2= desci%rcut2
          if( dij.ge.rcut .or. dik.ge.rcut ) cycle
          if( ityp.eq.101 ) then ! RK's original angular SF
!.....fcij's should be computed after rcs is determined
            call get_fc_dfc(dij,is,js,rcut,fcij,dfcij)
            call get_fc_dfc(dik,is,ks,rcut,fcik,dfcik)
            almbd= desci%prms(1)
            t2= (abs(almbd)+1d0)**2
            driki(1:3)= -rik(1:3)/dik
            drikk(1:3)= -driki(1:3)
!.....function value
            t1= (almbd +cs)**2
            eta3 = 0.5d0 /rcut**2
            texp = exp(-eta3*(dij2+dik2))
            tmp = t1/t2 *texp
            gsfi(isf)= gsfi(isf) +tmp*fcij*fcik
            gijk = gijk +tmp*fcij*fcik
!.....derivative
            dgdij= dfcij *fcik *tmp &
                 +tmp *(-2d0*eta3*dij) *fcij*fcik 
            dgdik= fcij *dfcik *tmp &
                 +tmp *(-2d0*eta3*dik) *fcij*fcik 
            dgsfi(1:3,isf,0)= dgsfi(1:3,isf,0) &
                 +dgdij*driji(1:3) +dgdik*driki(1:3)
            dgsfi(1:3,isf,jj)= dgsfi(1:3,isf,jj) +dgdij*drijj(1:3)
            dgsfi(1:3,isf,kk)= dgsfi(1:3,isf,kk) +dgdik*drikk(1:3)
            dgcs= 2d0*(almbd+cs)/t2 *fcij*fcik *texp 
            dgsfi(1:3,isf,0)= dgsfi(1:3,isf,0) +dgcs*dcsdi(1:3)
            dgsfi(1:3,isf,jj)= dgsfi(1:3,isf,jj) +dgcs*dcsdj(1:3)
            dgsfi(1:3,isf,kk)= dgsfi(1:3,isf,kk) +dgcs*dcsdk(1:3)
            igsfi(isf,0) = 1
            igsfi(isf,jj) = 1
            igsfi(isf,kk) = 1
          else if( ityp.eq.102 ) then  ! Similar to Behler's angular SF that includes fc(rjk)
!.....djk is required for Behler's angular SF (itype(isf)==102)
            xjk(1:3)= xk(1:3)-xj(1:3)
            rjk(1:3)= h(1:3,1)*xjk(1) +h(1:3,2)*xjk(2) +h(1:3,3)*xjk(3)
            djk2= rjk(1)**2 +rjk(2)**2 +rjk(3)**2
            if( djk2.ge.rcut2 ) cycle
            djk= sqrt(djk2)
!.....fcij's should be computed after rcs is determined
            call get_fc_dfc(dij,is,js,rcut,fcij,dfcij)
            call get_fc_dfc(dik,is,ks,rcut,fcik,dfcik)
            call get_fc_dfc(djk,js,ks,rcut,fcjk,dfcjk)
            almbd= desci%prms(1)
            zang = desci%prms(2)
            driki(1:3)= -rik(1:3)/dik
            drikk(1:3)= -driki(1:3)
            drjkj(1:3)= -rjk(1:3)/djk
            drjkk(1:3)= -drjkj(1:3)
!.....function value
            twozeta = 2d0**(-zang)
            t1= (almbd +cs)**zang *twozeta
            eta3 = 0.5d0 /rcut**2
            texp = exp(-eta3*(dij2+dik2+djk2))
            tmp = t1 *texp
!.....This part is different from itype(isf)==101 by the factor fcjk
            gsfi(isf)= gsfi(isf) +tmp*fcij*fcik *fcjk
!.....derivative
            dgdij= dfcij *fcik*fcjk *tmp &
                 +tmp *(-2d0*eta3*dij) *fcij*fcik*fcjk
            dgdik= fcij *dfcik*fcjk *tmp &
                 +tmp *(-2d0*eta3*dik) *fcij*fcik*fcjk 
            dgdjk= fcij *fcik *dfcjk *tmp &
                 +tmp *(-2d0*eta3*djk) *fcij*fcik*fcjk 
            dgsfi(1:3,isf,0)= dgsfi(1:3,isf,0) &
                 +dgdij*driji(1:3) +dgdik*driki(1:3)
            dgsfi(1:3,isf,jj)= dgsfi(1:3,isf,jj) &
                 +dgdij*drijj(1:3)+dgdjk*drjkj(1:3)
            dgsfi(1:3,isf,kk)= dgsfi(1:3,isf,kk) &
                 +dgdik*drikk(1:3)+dgdjk*drjkk(1:3)
            dgcs= zang*(almbd +cs)**(zang-1d0) *twozeta *fcij*fcik*fcjk *texp
            dgsfi(1:3,isf,0)= dgsfi(1:3,isf,0) +dgcs*dcsdi(1:3)
            dgsfi(1:3,isf,jj)= dgsfi(1:3,isf,jj) +dgcs*dcsdj(1:3)
            dgsfi(1:3,isf,kk)= dgsfi(1:3,isf,kk) +dgcs*dcsdk(1:3)
            igsfi(isf,0) = 1
            igsfi(isf,jj) = 1
            igsfi(isf,kk) = 1
          else if( ityp.eq.103 ) then ! cos(cos(thijk)*n*pi) w/o fc(rjk)
!.....fcij's should be computed after rcs is determined
            call get_fc_dfc(dij,is,js,rcut,fcij,dfcij)
            call get_fc_dfc(dik,is,ks,rcut,fcik,dfcik)
            an = desci%prms(1)
            driki(1:3)= -rik(1:3)/dik
            drikk(1:3)= -driki(1:3)
!.....Function value
            tmp = cos(cs*an*pi)
            gsfi(isf)= gsfi(isf) +tmp*fcij*fcik
!.....Derivative
            dgdij= tmp*dfcij*fcik
            dgdik= tmp*fcij*dfcik
            dgsfi(1:3,isf,0)= dgsfi(1:3,isf,0)&
                 +dgdij*driji(1:3) +dgdik*driki(1:3)
            dgsfi(1:3,isf,jj)= dgsfi(1:3,isf,jj) +dgdij*drijj(1:3)
            dgsfi(1:3,isf,kk)= dgsfi(1:3,isf,kk) +dgdik*drikk(1:3)
            dgcs= -an*pi*sin(cs*an*pi) *fcij*fcik
            dgsfi(1:3,isf,0)= dgsfi(1:3,isf,0) +dgcs*dcsdi(1:3)
            dgsfi(1:3,isf,jj)= dgsfi(1:3,isf,jj) +dgcs*dcsdj(1:3)
            dgsfi(1:3,isf,kk)= dgsfi(1:3,isf,kk) +dgcs*dcsdk(1:3)
            igsfi(isf,0) = 1
            igsfi(isf,jj) = 1
            igsfi(isf,kk) = 1
          else if( ityp.eq.104 ) then ! sin(cos(thijk)*n*pi) w/o fc(rjk)
!.....fcij's should be computed after rcs is determined
            call get_fc_dfc(dij,is,js,rcut,fcij,dfcij)
            call get_fc_dfc(dik,is,ks,rcut,fcik,dfcik)
            an = desci%prms(1)
            driki(1:3)= -rik(1:3)/dik
            drikk(1:3)= -driki(1:3)
!.....Function value
            tmp = sin(cs*an*pi)
            gsfi(isf)= gsfi(isf) +tmp*fcij*fcik
!.....Derivative
            dgdij= tmp*dfcij*fcik
            dgdik= tmp*fcij*dfcik
            dgsfi(1:3,isf,0)= dgsfi(1:3,isf,0)&
                 +dgdij*driji(1:3) +dgdik*driki(1:3)
            dgsfi(1:3,isf,jj)= dgsfi(1:3,isf,jj) +dgdij*drijj(1:3)
            dgsfi(1:3,isf,kk)= dgsfi(1:3,isf,kk) +dgdik*drikk(1:3)
            dgcs= an*pi*cos(cs*an*pi) *fcij*fcik
            dgsfi(1:3,isf,0)= dgsfi(1:3,isf,0) +dgcs*dcsdi(1:3)
            dgsfi(1:3,isf,jj)= dgsfi(1:3,isf,jj) +dgcs*dcsdj(1:3)
            dgsfi(1:3,isf,kk)= dgsfi(1:3,isf,kk) +dgcs*dcsdk(1:3)
            igsfi(isf,0) = 1
            igsfi(isf,jj) = 1
            igsfi(isf,kk) = 1
          else if( ityp.eq.105 ) then ! exp(-eta*(cos(thijk)-c)**2) w/o fc(rjk)
!.....fcij's should be computed after rcs is determined
            call get_fc_dfc(dij,is,js,rcut,fcij,dfcij)
            call get_fc_dfc(dik,is,ks,rcut,fcik,dfcik)
            eta= desci%prms(1)
            rs = desci%prms(2)
            driki(1:3)= -rik(1:3)/dik
            drikk(1:3)= -driki(1:3)
!.....Function value
            tmp = exp(-eta*(cs-rs)**2)
            gsfi(isf)= gsfi(isf) +tmp*fcij*fcik
!.....Derivative
            dgdij= tmp*dfcij*fcik
            dgdik= tmp*fcij*dfcik
            dgsfi(1:3,isf,0)= dgsfi(1:3,isf,0)&
                 +dgdij*driji(1:3) +dgdik*driki(1:3)
            dgsfi(1:3,isf,jj)= dgsfi(1:3,isf,jj) +dgdij*drijj(1:3)
            dgsfi(1:3,isf,kk)= dgsfi(1:3,isf,kk) +dgdik*drikk(1:3)
            dgcs= -2d0*eta*(cs-rs)*tmp *fcij*fcik
            dgsfi(1:3,isf,0)= dgsfi(1:3,isf,0) +dgcs*dcsdi(1:3)
            dgsfi(1:3,isf,jj)= dgsfi(1:3,isf,jj) +dgcs*dcsdj(1:3)
            dgsfi(1:3,isf,kk)= dgsfi(1:3,isf,kk) +dgcs*dcsdk(1:3)
            igsfi(isf,0) = 1
            igsfi(isf,jj) = 1
            igsfi(isf,kk) = 1
          endif

        enddo ! isf=1,...
      enddo ! kk=1,...
20    continue
    enddo ! jj=1,...

    return
  end subroutine desci_default
!=======================================================================
  subroutine calc_desc_cheby(namax,natm,nb,nnmax,h,tag,ra,lspr,rc &
       ,myid,mpi_world,l1st,iprint)
!-----------------------------------------------------------------------
!  Evaluate descriptors (Chebyshev polynomials)
!  and their derivatives wrt positions.
!-----------------------------------------------------------------------    
!  See Artrith et al., PRB96, 014112 (2017)
!-----------------------------------------------------------------------
    include "mpif.h"
    integer,intent(in):: namax,natm,nb,nnmax,lspr(0:nnmax,namax)
    integer,intent(in):: myid,mpi_world,iprint
    real(8),intent(in):: h(3,3),tag(namax),ra(3,namax),rc
    logical,intent(in):: l1st

    integer:: ia,ja,ka,jj,kk,is,js,ks,n,isf,ierr,isf0
    real(8):: x,xi(3),xj(3),xij(3),rij(3),dij,dij2,xk(3),xik(3) &
         ,rik(3),dik,dik2,driji(3),drijj(3),driki(3),drikk(3) &
         ,fcij,dfcij,fcik,dfcik,spijk,cs,dcsdj(3),dcsdk(3),dcsdi(3) &
         ,dgdcs,dgdij,dgdik,dgdr,wgt,rcut2,rcut3,rcut,wgts(nsf)

    if( l1st ) then
      if( .not. allocated(ts_cheby) ) then
        allocate(ts_cheby(0:max(nsf2,nsf3)),dts_cheby(0:max(nsf2,nsf3)))
        mem = mem +8*size(ts_cheby)*2
      endif

!.....Check the maximumx cutoff and given rc
      if( rc.lt.rcmax ) then
        if( myid.eq.0 ) then
          print *,'ERROR: cutoff radius rc is smaller than rcmax.'
          print *,'  rc,rcmax = ',rc,rcmax
        endif
        call mpi_finalize(ierr)
        stop
      endif
    endif

    gsf(1:nsf,1:nal)= 0d0
    dgsf(1:3,1:nsf,:,1:nal)= 0d0
    igsf(1:nsf,0:nnl,1:nal) = 0

    do ia=1,natm
      xi(1:3)= ra(1:3,ia)
      is= int(tag(ia))
      do jj=1,lspr(0,ia)
        ja= lspr(jj,ia)
        if( ja.eq.ia ) cycle
        xj(1:3)= ra(1:3,ja)
        xij(1:3)= xj(1:3)-xi(1:3)
        rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij2= rij(1)**2 +rij(2)**2 +rij(3)**2
        if( dij2.ge.rcmax2 ) cycle
        dij = sqrt(dij2)
        js= int(tag(ja))
        driji(1:3)= -rij(1:3)/dij
        drijj(1:3)= -driji(1:3)
!.....Rcut for 2-body common for all 2-body
        rcut2 = descs(1)%rcut
        x = 2d0*dij/rcut2 -1d0
        do isf = 1, nsf2*nsff
          wgt = 1d0
          if( mod(isf-1,nsff).eq.1 ) wgt = wgtsp(js)
          wgts(isf) = wgt
        enddo
        call chebyshev(x,nsf2,ts_cheby,dts_cheby)
!.....Since rcut is common in all the 2body terms,
!     this can be called outside the isf-loop.
        call get_fc_dfc(dij,is,js,rcut2,fcij,dfcij)
        do isf=1,nsf2*nsff  ! isf=[1,nsf2] for 2-body
          n = 1 +(isf-1)/nsff
          wgt = wgts(isf)
          gsf(isf,ia) = gsf(isf,ia) +ts_cheby(n)*fcij *wgt
!!$          if( ia.eq.1 ) print *,ia,ja,is,js,isf,wgt,gsf(isf,ia)
!.....Derivative
!     dgsf(ixyz,isf,jj,ia): derivative of isf-th basis of atom-ia
!     by ixyz coordinate of atom-jj.
!     jj=0 means derivative by atom-ia.
          dgdr = (2d0/rcut2*dts_cheby(n)*fcij +ts_cheby(n)*dfcij) *wgt
          dgsf(1:3,isf,0,ia)= dgsf(1:3,isf,0,ia) +driji(1:3)*dgdr
          dgsf(1:3,isf,jj,ia)= dgsf(1:3,isf,jj,ia) +drijj(1:3)*dgdr
          igsf(isf,0,ia) = 1
          igsf(isf,jj,ia)= 1
        enddo

!.....3-body terms
        if( .not. nsf3.gt.0 ) cycle
        isf0 = nsf2*nsff
!.....Skip 3-body if dij > rcij for 3-body which is common for all the 3-body terms
!!$        if( dij2.ge.rcs2(isf0+1) ) cycle
        rcut3 = descs(isf0+1)%rcut2
        rcut = descs(isf0+1)%rcut
        if( dij2.ge.rcut3 ) cycle
!.....Reset fcij and dfcij with rcij for 3-body
        call get_fc_dfc(dij,is,js,rcut,fcij,dfcij)
        do kk=1,lspr(0,ia)
          ka= lspr(kk,ia)
          ks= int(tag(ka))
          if( ka.eq.ia .or. ka.le.ja ) cycle
          xk(1:3)= ra(1:3,ka)
          xik(1:3)= xk(1:3)-xi(1:3)
          rik(1:3)= h(1:3,1)*xik(1) +h(1:3,2)*xik(2) +h(1:3,3)*xik(3)
          dik2= rik(1)**2 +rik(2)**2 +rik(3)**2
          if( dik2.gt.rcut3 ) cycle
          dik= sqrt(dik2)
          spijk= rij(1)*rik(1) +rij(2)*rik(2) +rij(3)*rik(3)
          cs= spijk/dij/dik
          dcsdj(1:3)= rik(1:3)/dij/dik -rij(1:3)*cs/dij**2
          dcsdk(1:3)= rij(1:3)/dij/dik -rik(1:3)*cs/dik**2
          dcsdi(1:3)= -dcsdj(1:3) -dcsdk(1:3)
          driki(1:3)= -rik(1:3)/dik
          drikk(1:3)= -driki(1:3)
          do isf = isf0+1, isf0+nsf3*nsff
            wgt = 1d0
            if( mod(isf-1,nsff).eq.1 ) wgt = wgtsp(js)*wgtsp(ks)            
            wgts(isf) = wgt
          enddo
          call get_fc_dfc(dik,is,ks,rcut,fcik,dfcik)
          call chebyshev(cs,nsf3,ts_cheby,dts_cheby)
          do isf= isf0+1, isf0+nsf3*nsff  ! isf=[isf0+1,isf0+nsf3] for 3-body
            n = 1 +(isf-isf0-1)/nsff
            wgt = wgts(isf)
            gsf(isf,ia) = gsf(isf,ia) +ts_cheby(n)*fcij*fcik *wgt
            dgdij = ts_cheby(n)*dfcij*fcik *wgt
            dgdik = ts_cheby(n)*fcij*dfcik *wgt
            dgdcs = dts_cheby(n)*fcij*fcik *wgt
!!$            if( ia.eq.8 ) then
!!$              print '(2i6,4i4,6es11.3)',ja,ka,js,ks,isf,n,wgt,dgdij,dgdik,dgdcs,ts_cheby(n),dts_cheby(n)
!!$            endif
            dgsf(1:3,isf,jj,ia)= dgsf(1:3,isf,jj,ia) +dgdcs*dcsdj(1:3) &
                 +dgdij*drijj(1:3)
            dgsf(1:3,isf,kk,ia)= dgsf(1:3,isf,kk,ia) +dgdcs*dcsdk(1:3) &
                 +dgdik*drikk(1:3)
            dgsf(1:3,isf,0,ia) = dgsf(1:3,isf, 0,ia) +dgdcs*dcsdi(1:3) &
                 +dgdij*driji(1:3) +dgdik*driki(1:3)
            igsf(isf,0,ia) = 1
            igsf(isf,jj,ia) = 1
            igsf(isf,kk,ia) = 1
          enddo
        enddo ! kk=1,...
      enddo ! jj=1,...
    enddo ! ia=1,natm

    return    
  end subroutine calc_desc_cheby
!=======================================================================
  subroutine desci_cheby(ia,namax,natm,nnmax,h,tag,ra,lspr,rc,iprint)
!-----------------------------------------------------------------------
!  Chebyshev descriptor of an atom-ia
!  See Artrith et al., PRB96, 014112 (2017)
!-----------------------------------------------------------------------
    integer,intent(in):: ia,namax,natm,nnmax,iprint, &
         lspr(0:nnmax,namax)
    real(8),intent(in):: h(3,3),tag(namax),ra(3,namax),rc

    integer:: ja,ka,jj,kk,is,js,ks,n,isf,ierr,isf0
    real(8):: x,xi(3),xj(3),xij(3),rij(3),dij,dij2,xk(3),xik(3) &
         ,rik(3),dik,dik2,driji(3),drijj(3),driki(3),drikk(3) &
         ,fcij,dfcij,fcik,dfcik,spijk,cs,dcsdj(3),dcsdk(3),dcsdi(3) &
         ,dgdcs,dgdij,dgdik,dgdr,wgt,rcut2,rcut3,rcut,wgts(nsf)

    xi(1:3)= ra(1:3,ia)
    is= int(tag(ia))
    do jj=1,lspr(0,ia)
      ja= lspr(jj,ia)
      if( ja.eq.ia ) cycle
      xj(1:3)= ra(1:3,ja)
      xij(1:3)= xj(1:3)-xi(1:3)
      rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
      dij2= rij(1)**2 +rij(2)**2 +rij(3)**2
      if( dij2.ge.rcmax2 ) cycle
      dij = sqrt(dij2)
      js= int(tag(ja))
      driji(1:3)= -rij(1:3)/dij
      drijj(1:3)= -driji(1:3)
!.....Rcut for 2-body common for all 2-body
      rcut2 = descs(1)%rcut
      x = 2d0*dij/rcut2 -1d0
      do isf = 1, nsf2*nsff
        wgt = 1d0
        if( mod(isf-1,nsff).eq.1 ) wgt = wgtsp(js)
        wgts(isf) = wgt
      enddo
      call chebyshev(x,nsf2,ts_cheby,dts_cheby)
!.....Since rcut is common in all the 2body terms,
!     this can be called outside the isf-loop.
      call get_fc_dfc(dij,is,js,rcut2,fcij,dfcij)
      do isf=1,nsf2*nsff  ! isf=[1,nsf2] for 2-body
        n = 1 +(isf-1)/nsff
        wgt = wgts(isf)
        gsfi(isf) = gsfi(isf) +ts_cheby(n)*fcij *wgt
!.....Derivative
!     dgsfi(ixyz,isf,jj): derivative of isf-th basis of atom-ia
!     by ixyz coordinate of atom-jj.
!     jj=0 means derivative by atom-ia.
        dgdr = (2d0/rcut2*dts_cheby(n)*fcij +ts_cheby(n)*dfcij) *wgt
        dgsfi(1:3,isf,0) = dgsfi(1:3,isf,0)  +driji(1:3)*dgdr
        dgsfi(1:3,isf,jj)= dgsfi(1:3,isf,jj) +drijj(1:3)*dgdr
        igsfi(isf,0) = 1
        igsfi(isf,jj)= 1
      enddo

!.....3-body terms
      if( .not. nsf3.gt.0 ) cycle
      isf0 = nsf2*nsff
!.....Skip 3-body if dij > rcij for 3-body which is common for all the 3-body terms
      rcut3 = descs(isf0+1)%rcut2
      rcut = descs(isf0+1)%rcut
      if( dij2.ge.rcut3 ) cycle
!.....Reset fcij and dfcij with rcij for 3-body
      call get_fc_dfc(dij,is,js,rcut,fcij,dfcij)
      do kk=1,lspr(0,ia)
        ka= lspr(kk,ia)
        ks= int(tag(ka))
        if( ka.eq.ia .or. ka.le.ja ) cycle
        xk(1:3)= ra(1:3,ka)
        xik(1:3)= xk(1:3)-xi(1:3)
        rik(1:3)= h(1:3,1)*xik(1) +h(1:3,2)*xik(2) +h(1:3,3)*xik(3)
        dik2= rik(1)**2 +rik(2)**2 +rik(3)**2
        if( dik2.gt.rcut3 ) cycle
        dik= sqrt(dik2)
        spijk= rij(1)*rik(1) +rij(2)*rik(2) +rij(3)*rik(3)
        cs= spijk/dij/dik
        dcsdj(1:3)= rik(1:3)/dij/dik -rij(1:3)*cs/dij**2
        dcsdk(1:3)= rij(1:3)/dij/dik -rik(1:3)*cs/dik**2
        dcsdi(1:3)= -dcsdj(1:3) -dcsdk(1:3)
        driki(1:3)= -rik(1:3)/dik
        drikk(1:3)= -driki(1:3)
        do isf = isf0+1, isf0+nsf3*nsff
          wgt = 1d0
          if( mod(isf-1,nsff).eq.1 ) wgt = wgtsp(js)*wgtsp(ks)            
          wgts(isf) = wgt
        enddo
        call get_fc_dfc(dik,is,ks,rcut,fcik,dfcik)
        call chebyshev(cs,nsf3,ts_cheby,dts_cheby)
        do isf= isf0+1, isf0+nsf3*nsff  ! isf=[isf0+1,isf0+nsf3] for 3-body
          n = 1 +(isf-isf0-1)/nsff
          wgt = wgts(isf)
          gsfi(isf) = gsfi(isf) +ts_cheby(n)*fcij*fcik *wgt
          dgdij = ts_cheby(n)*dfcij*fcik *wgt
          dgdik = ts_cheby(n)*fcij*dfcik *wgt
          dgdcs = dts_cheby(n)*fcij*fcik *wgt
          dgsfi(1:3,isf,jj)= dgsfi(1:3,isf,jj) +dgdcs*dcsdj(1:3) &
               +dgdij*drijj(1:3)
          dgsfi(1:3,isf,kk)= dgsfi(1:3,isf,kk) +dgdcs*dcsdk(1:3) &
               +dgdik*drikk(1:3)
          dgsfi(1:3,isf,0) = dgsfi(1:3,isf, 0) +dgdcs*dcsdi(1:3) &
               +dgdij*driji(1:3) +dgdik*driki(1:3)
          igsfi(isf,0) = 1
          igsfi(isf,jj) = 1
          igsfi(isf,kk) = 1
        enddo
      enddo ! kk=1,...
    enddo ! jj=1,...

    return
  end subroutine desci_cheby
!=======================================================================
  function fc1(r,rin,rout)
    implicit none
    real(8),intent(in):: r,rin,rout
    real(8):: fc1

    if( r.le.rin ) then
      fc1= 1d0
    else if( r.gt.rin .and. r.le.rout ) then
      fc1= 0.5d0 *(cos((r-rin)/(rout-rin)*pi)+1d0)
    else
      fc1= 0d0
    endif
    return
  end function fc1
!=======================================================================
  function dfc1(r,rin,rout)
    implicit none
    real(8),intent(in):: r,rin,rout
    real(8):: dfc1

    if( r.le.rin ) then
      dfc1= 0d0
    else if( r.gt.rin .and. r.le.rout ) then
      dfc1= -0.5d0*pi/(rout-rin) *sin((r-rin)/(rout-rin)*pi)
    else
      dfc1= 0d0
    endif
    return
  end function dfc1
!=======================================================================
  subroutine get_fc_dfc(r,isp,jsp,rcut,fc,dfc)
!
!  Calculate cutoff/switching function depending on r and r_outer
!
    real(8),intent(in):: r
    integer,intent(in):: isp,jsp
    real(8),intent(out):: fc,dfc,rcut

    real(8):: rin,rout
    
    fc= fc1(r,0d0,rcut)
    dfc= dfc1(r,0d0,rcut)
    
    return
  end subroutine get_fc_dfc
!=======================================================================
  subroutine read_params_desc(myid,mpi_world,iprint,specorder)
    use util, only: num_data
    implicit none
    include 'mpif.h'

    integer,intent(in):: myid,mpi_world,iprint
    character(len=3),intent(in):: specorder(nspmax)
!!$    real(8),intent(in):: rcin

    integer:: ierr,i,j,k,nc,ncoeff,nsp,isp,jsp,ksp,isf,ityp &
         ,ihl0,ihl1,ihl2,icmb(3),iap,jap,kap,ndat,is1,is2
    real(8):: rcut2,rcut3,rcut,time0,wgt
    logical:: lexist
    character(len=128):: ctmp,fname,cline,cmode
    character(len=3):: ccmb(3),csp

    if( lprmset_desc ) return

    time0 = mpi_wtime()

    if( myid.eq.0 ) then
      if( iprint.gt.1 ) print *,'read_params_desc...'
!.....read constants at the 1st call
      fname = trim(paramsdir)//'/'//trim(cpfname)
      inquire(file=trim(fname),exist=lexist)
      if( .not. lexist ) then
        if( myid.eq.0 ) then
          write(6,'(a)') ' [Error] '//trim(fname)//' does not exist !!!.'
        endif
        call mpi_finalize(ierr)
        stop
      endif
      open(ionum,file=trim(fname),status='old')
!.....num of symmetry functions, num of node in 1st hidden layer
10    read(ionum,'(a)') ctmp
      if( ctmp(1:1).eq.'!' .or. ctmp(1:1).eq.'#' ) then
        call parse_option(ctmp,iprint,ierr)
        goto 10
      else
        backspace(ionum)
      endif
!.....Read numbers of species and symmetry functions
      read(ionum,*) nsp, nsf
    endif

!.....Bcast nsp and nsf before allocating arrays
    call mpi_bcast(nsp,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(nsf,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(lcheby,1,mpi_logical,0,mpi_world,ierr)
    ngl = nsf

!.....Allocate arrays of lenths, nsp and/or nsf
    if( .not.allocated(descs) ) then
!!$      allocate(itype(nsf),cnst(max_ncnst,nsf),rcs(nsf),rcs2(nsf))
      allocate(descs(nsf),ilsf2(0:nsf,nspmax,nspmax) &
           ,ilsf3(0:nsf,nspmax,nspmax,nspmax))
!.....Also allocate group-LASSO/FS related variables,
!     which are not used in pmd but in fitpot
      allocate(mskgfs(ngl),msktmp(ngl),glval(0:ngl))
      mem = mem +8*ngl +8*ngl +8*(ngl+1)
      mskgfs(1:ngl) = 0d0
    endif
    if( lcheby ) then
      if( .not. allocated(wgtsp) ) allocate(wgtsp(nspmax))
      do i=1,nsp
        wgtsp(i)= dble(i)
      enddo
      call mpi_bcast(wgtsp,nspmax,mpi_real8,0,mpi_world,ierr)
      mem = mem + 8*size(wgtsp)
    endif

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
        nsff = 1
        nsf2 = 0
        nsf3 = 0
        if( nsp.gt.1 ) nsff = 2
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
            read(cline,*) ctmp, nsf2, rcut2
          else if( index(cline,'3-body').ne.0 ) then
            cmode = 'none'
            read(cline,*) ctmp, nsf3, rcut3
          else
            if( trim(cmode).eq.'Weight' ) then
              read(cline,*,end=30) csp, wgt
              isp = csp2isp(trim(csp),specorder)
              if( isp.gt.0 ) then
                wgtsp(isp) = wgt
                if( iprint.gt.0 ) write(6,'(5x,i2,a4,f6.1)') isp, trim(csp), wgt
              endif
!!$              read(ionum,*,end=30) isp, wgtsp(isp)
            endif
          endif
        enddo
30      continue
!.....Check nsf vs nsf2,nsf3
        if( nsf.ne.nsff*(nsf2+nsf3) ) then
          print *,'ERROR@read_params_desc: nsf != (nsf2+nsf3)*nsff with nsff,nsp=',nsff,nsp
          stop 1
        endif
!.....Set rcut for all isf
        do isf=1,nsf2*nsff
          descs(isf)%rcut = rcut2
          descs(isf)%rcut2 = rcut2**2
        enddo
        do isf=nsf2*nsff+1,nsf
          descs(isf)%rcut = rcut3
          descs(isf)%rcut2 = rcut3**2
        enddo

      else  ! not Chebyshev
        nsf2 = 0
        nsf3 = 0
        ilsf2(:,:,:) = 0
        ilsf3(:,:,:,:) = 0
        do isf=1,nsf
          read(ionum,*,end=20) ityp,(ccmb(k),k=1,ncomb_type(ityp)) &
               ,rcut,(cnst(j),j=1,ncnst_type(ityp))
          descs(isf)%itype = ityp
          isp = csp2isp(trim(ccmb(1)),specorder)
          jsp = csp2isp(trim(ccmb(2)),specorder)
          descs(isf)%isp = isp
          descs(isf)%jsp = jsp
          descs(isf)%rcut = rcut
          descs(isf)%rcut2 = rcut*rcut
          descs(isf)%nprm = ncnst_type(ityp)
          if( .not.allocated(descs(isf)%prms) ) &
               allocate(descs(isf)%prms(descs(isf)%nprm))
          do j=1,descs(isf)%nprm
            descs(isf)%prms(j) = cnst(j)
          enddo
          if( isp.lt.0 .or. jsp.lt.0 ) cycle
          if( ityp.le.100 ) then  ! 2-body
            nsf2 = nsf2 + 1
            is1 = min(isp,jsp)
            is2 = max(isp,jsp)
            ilsf2(0,is1,is2) = ilsf2(0,is1,is2) + 1
            ilsf2(ilsf2(0,is1,is2),is1,is2) = isf
          else if( ityp.le.200 ) then  ! 3-body
            nsf3 = nsf3 + 1
            ksp = csp2isp(trim(ccmb(3)),specorder)
            if( ksp.lt.0 ) cycle
            descs(isf)%ksp = ksp
            is1 = min(jsp,ksp)
            is2 = max(jsp,ksp)
            ilsf3(0,isp,is1,is2) = &
                 ilsf3(0,isp,is1,is2) +1
            ilsf3(ilsf3(0,isp,is1,is2),isp,is1,is2) = isf
          endif
        enddo  ! isf=1,nsf
20      continue
      endif ! lcheby
      close(ionum)
    endif ! myid.eq.0

!.....Broadcast cspline data
    call bcast_descs(myid,mpi_world,iprint)
    call mpi_bcast(nsf2,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(nsf3,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(nsff,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(ilsf2,size(ilsf2),mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(ilsf3,size(ilsf3),mpi_integer,0,mpi_world,ierr)

!.....Compute maximum rcut in all descriptors
    rcmax = 0d0
    rc3max = 0d0
    do isf=1,nsf
      rcut = descs(isf)%rcut
      rcmax = max(rcmax,rcut)
!!$      rcs2(i) = rcs(i)**2
      if( descs(isf)%itype.gt.100 ) rc3max = max(rc3max,rcut)
    enddo
    rcmax2 = rcmax**2

    lprmset_desc = .true.

    if( myid.eq.0 .and. iprint.gt.2 ) print *,'read_params_desc done'

    time = time +(mpi_wtime() -time0)

    return
  end subroutine read_params_desc
!=======================================================================
  subroutine parse_option(cline,iprint,ierr)
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
    implicit none
    character(len=*),intent(in):: cline
    integer,intent(in):: iprint
    integer,intent(out):: ierr

    integer:: iopt1,isp
    real(8):: opt1, opt2
    character(len=10):: c1,copt
    logical:: lopt

    ierr = 0
    if( index(cline,'Chebyshev:').ne.0 .or. &
         index(cline,'chebyshev:').ne.0 ) then
      read(cline,*) c1, copt, lopt
      lcheby = lopt
      if( iprint.gt.0 ) then
        print *,''
        print '(a)', ' Chebyshev series for descriptors.'
      endif
    endif
    
  end subroutine parse_option
!=======================================================================
  subroutine bcast_descs(myid,mpi_world,iprint)
    include 'mpif.h'
    integer,intent(in):: myid,mpi_world,iprint

    integer:: i,ierr

    do i=1,nsf
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
  subroutine set_paramsdir_desc(dname)
!
!  Accessor routine to set paramsdir.
!
    implicit none
    character(len=*),intent(in):: dname

    paramsdir = trim(dname)
    return
  end subroutine set_paramsdir_desc
!=======================================================================
  subroutine set_params_desc(descs_in,nsf_in,nsf2_in,nsf3_in, &
       nsff_in, ilsf2_in, ilsf3_in, lcheby_in, cnst_in, wgtsp_in )
!
!  Accessor routine to set desc parameters from outside (fitpot).
!  Curretnly this routine is supposed to be called only on serial run.
!
    integer,intent(in):: nsf_in,nsf2_in,nsf3_in,nsff_in, &
         ilsf2_in(0:nsf_in,nspmax,nspmax), ilsf3_in(0:nsf_in,nspmax,nspmax,nspmax)
    logical,intent(in):: lcheby_in
    real(8),intent(in):: cnst_in(max_ncnst)
    type(desc),intent(in):: descs_in(nsf_in)
    real(8),intent(in),optional:: wgtsp_in(nspmax)

    integer:: i,nsp
    real(8):: rcut

    if( lprmset_desc ) return

    if( lcheby .and. .not. present(wgtsp_in) ) then
      print *,'ERROR: wgtsp_in should be present if lcheby == .true.'
      stop
    endif

!.....Copy and allocate if needed
    nsf = nsf_in
    lcheby = lcheby_in
    if( .not. allocated(descs) ) then
      allocate(descs(nsf),ilsf2(0:nsf,nspmax,nspmax), &
           ilsf3(0:nsf,nspmax,nspmax,nspmax))
      if( lcheby .and. .not. allocated(wgtsp) ) allocate(wgtsp(nspmax))
    endif
    nsf2 = nsf2_in
    nsf3 = nsf3_in
    nsff = nsff_in
    ilsf2(:,:,:) = ilsf2_in(:,:,:)
    ilsf3(:,:,:,:) = ilsf3_in(:,:,:,:)
    cnst(:) = cnst_in(:)
    do i=1,nsf
      descs(i)%itype = descs_in(i)%itype
      descs(i)%isp   = descs_in(i)%isp
      descs(i)%jsp   = descs_in(i)%jsp
      descs(i)%ksp   = descs_in(i)%ksp
      descs(i)%rcut  = descs_in(i)%rcut
      descs(i)%rcut2 = descs_in(i)%rcut2
      descs(i)%nprm  = descs_in(i)%nprm
      if( .not. lcheby ) then
        if( allocated(descs(i)%prms) ) deallocate(descs(i)%prms)
        allocate(descs(i)%prms(descs(i)%nprm))
        descs(i)%prms(:) = descs_in(i)%prms(:)
      endif
    enddo
    if( present(wgtsp_in) ) wgtsp(:) = wgtsp_in(:)

!.....Compute maximum rcut in all descriptors
    rcmax = 0d0
    rc3max = 0d0
    do i=1,nsf
      rcut = descs(i)%rcut
      rcmax = max(rcmax,rcut)
      if( descs(i)%itype.gt.100 ) rc3max = max(rc3max,rcut)
    enddo
    rcmax2 = rcmax**2

    lprmset_desc = .true.
    
    return
  end subroutine set_params_desc
!=======================================================================
  subroutine write_descs(ionum,natm,namax,nnmax,lspr,tag)
!
!   Write out descriptor data (gsf,dgsf,igsf).
!   Buffer atom indices are replaced to resident atom ones.
!
    use util,only: itotOf
    implicit none
    integer,intent(in):: ionum
    integer,intent(in):: natm,namax,nnmax,lspr(0:nnmax,namax)
    real(8),intent(in):: tag(namax)

    integer:: ia,jj,ja,jra,isf,ihl0
    real(8),allocatable:: dgsfo(:,:,:,:)
    integer(2),allocatable:: igsfo(:,:,:)
!!$    integer,external:: itotOf

    open(ionum,file='out.desc.gsf',status='replace',form='unformatted')
    write(ionum) nsf
    do ia=1,natm
      write(ionum) (gsf(ihl0,ia),ihl0=1,nsf)
    enddo
    close(ionum)

    allocate(dgsfo(3,natm,nsf,natm),igsfo(nsf,natm,natm))
!.....reduce d(i)gsf data of buffer atoms to those of resident atoms
    dgsfo(1:3,1:natm,1:nsf,1:natm)= 0d0
    igsfo(1:nsf,1:natm,1:natm)= 0
    do ia=1,natm
      do jj=0,lspr(0,ia)
        if( jj.eq.0 ) then
          ja= ia
        else
          ja= lspr(jj,ia)
        endif
        jra= itotOf(tag(ja))
        do isf=1,nsf
          dgsfo(1:3,jra,isf,ia)= dgsfo(1:3,jra,isf,ia) &
               +dgsf(1:3,isf,jj,ia)
          igsfo(isf,jra,ia) = igsf(isf,jj,ia)
        enddo
      enddo
    enddo
!.....write
    open(ionum+1,file='out.desc.dgsf',status='replace',form='unformatted')
    do ia=1,natm
      do isf=1,nsf
        write(ionum+1) (dgsfo(1:3,jra,isf,ia),jra=1,natm)
      enddo
    enddo
    close(ionum+1)
!.....write igsf to ionum+1
    open(ionum+2,file='out.desc.igsf',status='replace',form='unformatted')
    do ia=1,natm
      do ja=1,natm
        write(ionum+2) (igsfo(isf,ja,ia),isf=1,nsf)
      enddo
    enddo
    close(ionum+2)

    deallocate(dgsfo,igsfo)
  end subroutine write_descs
!=======================================================================
  subroutine get_ints(nsfo,nalo,nnlo)
!
!  Access to some integers.
!
    integer,intent(out):: nsfo,nalo,nnlo

    nsfo = nsf
    nalo = nal
    nnlo = nnl
    return
  end subroutine get_ints
!=======================================================================
  subroutine get_descs(nsfo,nalo,nnlo,gsfo,dgsfo,igsfo)
!
!  Access to descriptors from outside
!
    integer,intent(in):: nsfo,nalo,nnlo
    real(8),intent(out):: gsfo(nsfo,nalo)
    real(8),intent(out),optional:: dgsfo(3,nsfo,0:nnlo,nalo)
    integer(2),intent(out),optional:: igsfo(nsfo,0:nnlo,nalo)

    if( nalo.gt.nal .or. nnlo.gt.nnl ) then
      print *,'ERROR: nalo/nnlo is greater than nal/nnl.'
      print *,'nalo,nal=',nalo,nal
      print *,'nnlo,nnl=',nnlo,nnl
      stop
    endif

    gsfo(:,1:nalo) = gsf(:,1:nalo)
    if( present(dgsfo) ) dgsfo(:,:,0:nnlo,1:nalo) = dgsf(:,:,0:nnlo,1:nalo)
    if( present(igsfo) ) igsfo(:,0:nnlo,1:nalo) = igsf(:,0:nnlo,1:nalo)
    return
  end subroutine get_descs
!=======================================================================
  subroutine set_descs(nsfo,nalo,nnlo,gsfo,dgsfo,igsfo)
!
!  Set descriptors from outside (fitpot).
!
    integer,intent(in):: nsfo,nalo,nnlo
    real(8),intent(in):: gsfo(nsfo,nalo),dgsfo(3,nsfo,0:nnlo,nalo)
    integer(2),intent(in):: igsfo(nsfo,0:nnlo,nalo)

    integer:: isf

    if( nsf.ne.nsfo ) stop 'ERROR @set_descs: nsf.ne.nsfo, which should not happen.'

    if( nal.lt.nalo .or. nnl.lt.nnlo )  then
      if( allocated(gsf) ) deallocate(gsf,dgsf,igsf)
      nsf = nsfo
      nal = nalo
      nnl = nnlo
      allocate( gsf(nsf,nal),dgsf(3,nsf,0:nnl,nal) &
           ,igsf(nsf,0:nnl,nal) )
    endif
    
    gsf(:,1:nalo) = gsfo(:,1:nalo)
    dgsf(:,:,0:nnlo,1:nalo) = dgsfo(:,:,0:nnlo,1:nalo)
    igsf(:,0:nnlo,1:nalo) = igsfo(:,0:nnlo,1:nalo)
    return
  end subroutine set_descs
!=======================================================================
  subroutine get_dsgnmat_force(dgsfa,mpi_world)
!
!  Special routine to compute the contribution of each descriptors
!  to force and return design matrix of force-matching,
!  which is used only in fitpot.
!
    use pmdio,only: namax,nbmax
    use pmdvars,only: natm,nb,lsb,nex,lsrc,myparity,nn &
         ,lspr,tcom
    integer,intent(in):: mpi_world
    real(8),allocatable,intent(out):: dgsfa(:,:,:)

    integer:: ia,jj,ja,isf

    if( allocated(dgsfa) ) then
      if( size(dgsfa).ne.3*nsf*namax ) then
        deallocate(dgsfa)
        allocate(dgsfa(3,nsf,namax))
      endif
    else
      allocate(dgsfa(3,nsf,namax))
    endif

    dgsfa(1:3,1:nsf,1:namax) = 0d0
    do ia=1,natm
      do jj=1,lspr(0,ia)
        ja = lspr(jj,ia)
        do isf=1,nsf
          dgsfa(1:3,isf,ja) = dgsfa(1:3,isf,ja) +dgsf(1:3,isf,jj,ia)
        enddo
      enddo
      do isf=1,nsf
        dgsfa(1:3,isf,ia) = dgsfa(1:3,isf,ia) +dgsf(1:3,isf,0,ia)
      enddo
    enddo
    call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,dgsfa,3*nsf)
    return
  end subroutine get_dsgnmat_force
!=======================================================================
  subroutine chebyshev(x,nmax,ts,dts)
!
!  Compute Chebyshev polynomial series at x Tn(x) up to nmax order.
!  Derivative of Tn(x), dTn(x)/dx, is also computed.
!  Returns ts and dts, which correspond to Tn(x) and dTn(x)
!
    real(8),intent(in):: x
    integer,intent(in):: nmax
    real(8),intent(out):: ts(0:nmax),dts(0:nmax)

    integer:: n

    ts(0) = 1d0
    ts(1) = x
    dts(0)= 0d0
    dts(1)= 1d0
    do n=2,nmax
      ts(n) = 2d0*x*ts(n-1) -ts(n-2)
      dts(n)= 2d0*(ts(n-1) +x*dts(n-1)) -dts(n-2)
    enddo
    return
  end subroutine chebyshev
!=======================================================================
  function mem_descriptor()
    integer:: mem_descriptor

    mem_descriptor = mem
    return
  end function mem_descriptor
!=======================================================================
  function time_descriptor()
    real(8):: time_descriptor

    time_descriptor = time
    return
  end function time_descriptor
!=======================================================================
  subroutine debug_descs()
!
!  Write descs for debugging
!
    integer:: i,isf
    
    do i=1,2
      do isf=1,nsf
        print '(a,2i5,2es12.5)','i,isf,rcut,gsf=' &
             ,i,isf,descs(isf)%rcut,gsf(isf,i)
      enddo
    enddo
  end subroutine debug_descs
end module descriptor
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
