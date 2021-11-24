module linreg
!-----------------------------------------------------------------------
!                     Last modified: <2021-11-24 15:55:33 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
!  Parallel implementation of linear regression potential for pmd
!    - 2014.06.11 by R.K. 1st implementation
!    - 2015.02.03 by R.K. extended to multiple species
!-----------------------------------------------------------------------
  use vector,only: dot
  implicit none
  save
  
  character(len=128):: paramsdir = '.'
!.....parameter file name
  character(128),parameter:: cpfname= 'in.params.linreg'
  character(128),parameter:: ccfname='in.const.linreg'
  character(128),parameter:: cmbfname='in.comb.linreg'

  logical:: lprmset_linreg = .false.

!.....parameters
  real(8),allocatable:: coeff(:)
!.....constants
  integer:: nelem,nexp,nsp
  integer,allocatable:: itype(:),icmb(:,:)
  real(8),allocatable:: cnst(:,:),exps(:)
!.....function types and num of constatns for types
  integer,parameter:: max_ncnst= 3
  integer,parameter:: ncnst_type(1:10)= &
       (/ 3, &  ! Gaussian
          2, &  ! cosine
          1, &  ! angular
          1, &  ! poly-1
          1, &  ! poly-2
          1, &  ! poly-4
          1, &  ! poly-6
          1, &  ! poly-8
          1, &  ! poly-10
          1  &  ! poly-12
          /)
!.....pairs (2) or triplets (3)
  integer,parameter:: max_natm= 3
  integer,parameter:: natm_type(1:10)= &
       (/ 2, &  ! Gaussian
          2, &  ! cosine
          3, &  ! angular
          2, &  ! poly-1
          2, &  ! poly-2
          2, &  ! poly-4
          2, &  ! poly-6
          2, &  ! poly-8
          2, &  ! poly-10
          2  &  ! poly-12
          /)

contains
  subroutine force_linreg_old(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_world,myid,epi,epot,nismax,lstrs,iprint)
    implicit none
    include "mpif.h"
    include "./params_unit.h"
    integer,intent(in):: namax,natm,nnmax,nismax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_world,myid,lspr(0:nnmax,namax),nex(3)
    real(8),intent(in):: ra(3,namax),tag(namax) &
         ,h(3,3),hi(3,3),sv(3,6)
    real(8),intent(inout):: rc
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical:: lstrs

!.....local
    integer:: i,j,k,l,m,n,ixyz,jxyz,is,js,ks,ierr,nbl,ia,ielem &
         ,iwgt
    real(8):: rcin,b_na,at(3),epotl,wgt,aexp,bnai,apot,epott
    real(8),save,allocatable:: fat(:,:),dbna(:,:,:)
!.....1st call
    logical,save:: l1st=.true.

    if( l1st ) then
!.....read in.params.linreg
      call read_params(myid,mpi_world,rcin)
!.....reset rc
      if( myid.eq.0 ) then
        write(6,'(a,f10.5,a,f10.5)') &
             ' Cutoff radius rc may have been changed from '&
             ,rc,' to ',rcin
      endif
      rc= rcin
      allocate(fat(3,namax),dbna(3,nelem,namax))
      l1st= .false.
    endif

#ifdef __FITPOT__
    open(80,file='out.basis.linreg',status='replace')
    write(80,'(3i10)') natm,nelem
    open(81,file='out.dbasis.linreg',status='replace')
    write(81,'(3i10)') natm,nelem
#endif

    epotl= 0d0

    dbna(1:3,1:nelem,1:natm+nb)= 0d0
    do ia=1,natm
      iwgt= 0
#ifdef __3BODY__
      apot= 0d0
#endif
      do ielem=1,nelem
        iwgt= iwgt +1
        wgt= coeff(iwgt)
        aexp= exps(ielem)
        call bfunc(ia,natm,namax,nnmax,ra,lspr,h,tag,dbna,rc &
             ,ielem,aexp,bnai)
#ifdef __FITPOT__
        write(80,'(2i10,f5.1,es23.14e3)') ia,ielem,aexp,bnai
!        write(81,'(2i10,f5.1,3es23.14e3)') ia,ielem,aexp,fat(1:3,ia)
#endif
!        write(6,*) ' ia,ielem,bnai,wgt=',ia,ielem,bnai,wgt
        epotl=epotl +bnai*wgt
        epi(ia)= epi(ia) +bnai*wgt
#ifdef __3BODY__
        if( itype(ielem).eq.3 ) then
          apot= apot +bnai*wgt
        endif
#endif
      enddo
#ifdef __3BODY__
      write(6,'(a,i8,es22.14)') ' 3-body term:',ia,apot
#endif
    enddo

!.....sum up forces
    do ia=1,natm+nb
      do ielem=1,nelem
        wgt= coeff(ielem)
        aa(1:3,ia)= aa(1:3,ia) &
             -dbna(1:3,ielem,ia)*wgt
      enddo
    enddo

#ifdef __FITPOT__
    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,dbna,3*nelem)
    do ia=1,natm
      do ielem=1,nelem
        write(81,'(2i10,3es23.14e3)') ia,ielem,dbna(1:3,ielem,ia)
      enddo
    enddo
    close(80)
    close(81)
#endif

    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,aa,3)

!-----gather epot
    call mpi_allreduce(epotl,epott,1,mpi_double_precision &
         ,mpi_sum,mpi_world,ierr)
    if( iprint.gt.2 ) print *,'linreg epot = ',epott
    epot= epot +epott
    
    return
  end subroutine force_linreg_old
!=======================================================================
  subroutine force_linreg(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rcin,lspr &
       ,mpi_world,myid,epi,epot,nismax,lstrs,iprint,l1st)
    use descriptor,only: gsfi,dgsfi,igsfi,nsf,calc_desci,make_gsf_arrays &
         ,pre_desci
    implicit none
    include "mpif.h"
    include "./params_unit.h"
    integer,intent(in):: namax,natm,nnmax,nismax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_world,myid,lspr(0:nnmax,namax),nex(3)
    real(8),intent(in):: ra(3,namax),tag(namax) &
         ,h(3,3),hi(3,3),sv(3,6)
    real(8),intent(inout):: rcin
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: l1st 
    logical:: lstrs

!.....local
    integer:: i,j,k,l,m,n,ixyz,jxyz,is,js,ks,ierr,nbl,ia,ielem &
         ,iwgt,isf,ja,jj
    real(8):: b_na,at(3),epotl,wgt,aexp,bnai,apot,epott,tmp
    real(8):: xi(3),xj(3),xji(3),sji,rji(3),dji
    real(8),save,allocatable:: aal(:,:),strsl(:,:,:)
    real(8),save:: rcin2 

    call pre_desci(namax,natm,nnmax,lspr,iprint,rcin)
    call make_gsf_arrays(l1st,namax,natm,tag,nnmax,lspr &
         ,myid,mpi_world,iprint)

    if( l1st ) then
      if( allocated(aal) ) deallocate(aal,strsl)
      allocate(aal(3,namax),strsl(3,3,namax))

    endif

    if( size(aal).lt.3*namax ) then
      deallocate(aal,strsl)
      allocate(aal(3,namax),strsl(3,3,namax))
      rcin2 = rcin*rcin
    endif

    aal(1:3,1:namax) = 0d0
    strsl(1:3,1:3,1:namax) = 0d0
    epotl= 0d0

!.....Energy
    do ia=1,natm
      call calc_desci(ia,namax,natm,nnmax,h &
           ,tag,ra,lspr,rcin,iprint)
      do isf=1,nsf
        wgt = coeff(isf)
        tmp = wgt*gsfi(isf)
        epotl = epotl +tmp
        epi(ia) = epi(ia) +tmp
      enddo
!.....Force
      do jj=1,lspr(0,ia)
        ja = lspr(jj,ia)
        do isf=1,nsf
          if( igsfi(isf,jj).eq.0 ) cycle
          wgt = coeff(isf)
          aal(1:3,ja) = aal(1:3,ja) &
               -wgt*dgsfi(1:3,isf,jj)
        enddo
      enddo
!.....Atom-ia
      do isf=1,nsf
        wgt = coeff(isf)
        aal(1:3,ia) = aal(1:3,ia) &
             -wgt*dgsfi(1:3,isf,0)
      enddo
!.....Stress
      if( .not.lstrs ) cycle
      xi(1:3)= ra(1:3,ia)
      is = int(tag(ia))
      do jj=1,lspr(0,ia)
        ja = lspr(jj,ia)
        xj(1:3) = ra(1:3,ja)
        xji(1:3) = xj(1:3)-xi(1:3)
        rji(1:3) = h(1:3,1)*xji(1) +h(1:3,2)*xji(2) +h(1:3,3)*xji(3)
        dji = rji(1)*rji(1) +rji(2)*rji(2) +rji(3)*rji(3)
        if( dji.ge.rcin2 ) cycle
        dji = sqrt(dji)
        if( dji.ge.rcin ) exit
        js = int(tag(ja))
        do isf=1,nsf
          if( igsfi(isf,jj).eq.0 ) cycle
          wgt = coeff(isf)
          do ixyz=1,3
            do jxyz=1,3
              sji = -wgt*dgsfi(jxyz,isf,jj)*rji(ixyz)
              strsl(ixyz,jxyz,ja) = strsl(ixyz,jxyz,ja) +sji
              strsl(ixyz,jxyz,ia) = strsl(ixyz,jxyz,ia) +sji
            enddo
          enddo
        enddo
      enddo
    enddo

    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,aal,3)
    aa(1:3,1:natm) = aa(1:3,1:natm) +aal(1:3,1:natm)
    if( lstrs ) then
      call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
           ,nn,mpi_world,strsl,9)
      strs(1:3,1:3,1:natm) = strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)*0.5d0
    endif

!-----gather epot
    call mpi_allreduce(epotl,epott,1,mpi_real8,mpi_sum,mpi_world,ierr)
    epot= epot +epott
    if( myid.eq.0 .and. iprint.gt.2 ) print *,'epot linreg = ',epott
    
    return
  end subroutine force_linreg
!=======================================================================
  function fc(r,rc)
    implicit none
    real(8),intent(in):: r,rc
    real(8):: fc
    real(8),parameter:: pi= 3.14159265358979d0

    fc= 0.5d0 *(cos(r/rc*pi)+1d0)
    return
  end function fc
!=======================================================================
  function dfc(r,rc)
    implicit none
    real(8),intent(in):: r,rc
    real(8):: dfc
    real(8),parameter:: pi= 3.14159265358979d0

    dfc= -pi/2/rc *sin(r/rc*pi)
    return
  end function dfc
!=======================================================================
  subroutine bfunc(ia,natm,namax,nnmax,ra,lspr,h,tag,dbna,rc &
       ,ielem,aexp,bnai)
!
!  basis function in the linear regression potetnial
!
    implicit none
    integer,intent(in):: ia,natm,namax,nnmax,lspr(0:nnmax,natm) &
         ,ielem
    real(8),intent(in):: ra(3,namax),h(3,3),tag(namax),rc,aexp
    real(8),intent(out):: bnai,dbna(3,nelem,namax)

    integer:: ja,jj,ka,kk,is,js,ks
    real(8):: xi(3),xj(3),xij(3),rij(3),r,dirij(3),djrij(3),tmp &
         ,fcij,xk(3),xik(3),rik(3),rj,rk,fcik,dkrik(3),dirik(3) &
         ,f3,dfcj,dfck,tmp2,cs,acnst,dcosi(3),dcosj(3),dcosk(3)

    bnai= 0d0
    xi(1:3)= ra(1:3,ia)
    is= int(tag(ia))
    if( itype(ielem).eq.3 ) then ! angular (3-body) basis
      if( is.ne.icmb(1,ielem) ) return
      do jj=1,lspr(0,ia)
        ja= lspr(jj,ia)
        if( ja.eq.ia ) cycle
        js= int(tag(ja))
        if( js.ne.icmb(2,ielem) .and. js.ne.icmb(3,ielem) ) cycle
        xj(1:3)= ra(1:3,ja)
        xij(1:3)= xj(1:3)-xi(1:3)
        rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        rj= sqrt(rij(1)**2 +rij(2)**2 +rij(3)**2)
        if( rj.gt.rc ) cycle
        fcij= fc(rj,rc)
        !.....another loop on neighbors for angular basis
        do kk=1,lspr(0,ia)
          ka= lspr(kk,ia)
          if( ka.le.ja ) cycle
          ks= int(tag(ka))
          if( .not.( (js.eq.icmb(2,ielem).and.ks.eq.icmb(3,ielem)) .or. &
               (js.eq.icmb(3,ielem).and.ks.eq.icmb(2,ielem))) ) cycle
          xk(1:3)= ra(1:3,ka)
          xik(1:3)= xk(1:3)-xi(1:3)
          rik(1:3)= h(1:3,1)*xik(1) +h(1:3,2)*xik(2) +h(1:3,3)*xik(3)
          rk= sqrt(rik(1)**2 +rik(2)**2 +rik(3)**2)
          if( rk.gt.rc ) cycle
          fcik= fc(rk,rc)
          bnai= bnai + func3(rij,rj,rik,rk,ielem,is,js,ks) *fcij *fcik 
        enddo
      enddo
      !.....bnai will be used in the following loop.
      !.....Therefore these two loops cannot be merged.
      tmp= aexp *bnai**(aexp-1)
      do jj=1,lspr(0,ia)
        ja= lspr(jj,ia)
        if( ja.eq.ia ) cycle
        js=int(tag(ja))
        if( js.ne.icmb(2,ielem) .and. js.ne.icmb(3,ielem) ) cycle
        xj(1:3)= ra(1:3,ja)
        xij(1:3)= xj(1:3)-xi(1:3)
        rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        rj= sqrt(rij(1)**2 +rij(2)**2 +rij(3)**2)
        if( rj.gt.rc ) cycle
        fcij= fc(rj,rc)
        dfcj= dfc(rj,rc)
        dirij(1:3)= -rij(1:3)/rj
        djrij(1:3)= -dirij(1:3)
        do kk=1,lspr(0,ia)
          ka= lspr(kk,ia)
          if( ka.le.ja ) cycle
          xk(1:3)= ra(1:3,ka)
          ks= int(tag(ka))
          if( .not.( (js.eq.icmb(2,ielem).and.ks.eq.icmb(3,ielem)) .or. &
               (js.eq.icmb(3,ielem).and.ks.eq.icmb(2,ielem))) ) cycle
          xik(1:3)= xk(1:3)-xi(1:3)
          rik(1:3)= h(1:3,1)*xik(1) +h(1:3,2)*xik(2) +h(1:3,3)*xik(3)
          rk= sqrt(rik(1)**2 +rik(2)**2 +rik(3)**2)
          if( rk.gt.rc ) cycle
          fcik= fc(rk,rc)
          dfck= dfc(rk,rc)
          dirik(1:3)= -rik(1:3)/rk
          dkrik(1:3)= -dirik(1:3)
          if( itype(ielem).eq.3 ) then
            acnst= cnst(1,ielem)
            cs= dot(rij,rik)/rj/rk
            f3= func3(rij,rj,rik,rk,ielem,is,js,ks)
            dbna(1:3,ielem,ia)= dbna(1:3,ielem,ia) &
                 +f3*fcik *dfcj*dirij(1:3) *tmp &
                 +f3*fcij *dfck*dirik(1:3) *tmp
            dbna(1:3,ielem,ja)= dbna(1:3,ielem,ja) &
                 +f3*fcik *dfcj*djrij(1:3) *tmp
            dbna(1:3,ielem,ka)= dbna(1:3,ielem,ka) &
                 +f3*fcij *dfck*dkrik(1:3) *tmp
!!$            fat(1:3,ia)= fat(1:3,ia) &
!!$                 -f3*fcik *dfcj*dirij(1:3) *tmp &
!!$                 -f3*fcij *dfck*dirik(1:3) *tmp
!!$            fat(1:3,ja)= fat(1:3,ja) &
!!$                 -f3*fcik *dfcj*djrij(1:3) *tmp
!!$            fat(1:3,ka)= fat(1:3,ka) &
!!$                 -f3*fcij *dfck*dkrik(1:3) *tmp
            dcosj(1:3)= rik(1:3)/rj/rk -rij(1:3)/rj*cs/rj
            dcosk(1:3)= rij(1:3)/rj/rk -rik(1:3)/rk*cs/rk
            dcosi(1:3)= -dcosj(1:3) -dcosk(1:3)
            tmp2= 2d0*(acnst+cs)/(abs(acnst)+1d0)**2
            dbna(1:3,ielem,ia)= dbna(1:3,ielem,ia) &
                 +dcosi(1:3)*tmp2*fcij*fcik *tmp
            dbna(1:3,ielem,ja)= dbna(1:3,ielem,ja) &
                 +dcosj(1:3)*tmp2*fcij*fcik *tmp
            dbna(1:3,ielem,ka)= dbna(1:3,ielem,ka) &
                 +dcosk(1:3)*tmp2*fcij*fcik *tmp
!!$            fat(1:3,ia)= fat(1:3,ia) &
!!$                 -dcosi(1:3)*tmp2*fcij*fcik *tmp
!!$            fat(1:3,ja)= fat(1:3,ja) &
!!$                 -dcosj(1:3)*tmp2*fcij*fcik *tmp
!!$            fat(1:3,ka)= fat(1:3,ka) &
!!$                 -dcosk(1:3)*tmp2*fcij*fcik *tmp
          endif
        enddo
      enddo
      
    else ! 2-body basis
      if( is.ne.icmb(1,ielem) .and. is.ne.icmb(2,ielem) ) return
      do jj=1,lspr(0,ia)
        ja= lspr(jj,ia)
        if( ja.eq.ia ) cycle
        js= int(tag(ja))
        if( .not. ((is.eq.icmb(1,ielem).and.js.eq.icmb(2,ielem)) .or. &
             (is.eq.icmb(2,ielem).and.js.eq.icmb(1,ielem))) ) cycle
        xj(1:3)= ra(1:3,ja)
        xij(1:3)= xj(1:3)-xi(1:3)
        rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        r= sqrt(rij(1)**2 +rij(2)**2 +rij(3)**2)
        if( r.gt.rc ) cycle
        fcij= fc(r,rc)
        bnai= bnai +func2(r,ielem,is,js)*fcij
      enddo
      !.....bnai will be used in the following loop.
      !.....Therefore these two loops cannot be merged.
      do jj=1,lspr(0,ia)
        ja= lspr(jj,ia)
        if( ja.eq.ia ) cycle
        js= int(tag(ja))
        if( .not. ((is.eq.icmb(1,ielem).and.js.eq.icmb(2,ielem)) .or. &
             (is.eq.icmb(2,ielem).and.js.eq.icmb(1,ielem))) ) cycle
        xj(1:3)= ra(1:3,ja)
        xij(1:3)= xj(1:3)-xi(1:3)
        rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        r= sqrt(rij(1)**2 +rij(2)**2 +rij(3)**2)
        if( r.gt.rc ) cycle
        dirij(1:3)= -rij(1:3)/r
        djrij(1:3)= -dirij(1:3)
        tmp= dfunc2(r,ielem,is,js)*fc(r,rc) &
             +func2(r,ielem,is,js)*dfc(r,rc)
        dbna(1:3,ielem,ia)= dbna(1:3,ielem,ia) &
             +dirij(1:3)*tmp *aexp*bnai**(aexp-1)
        dbna(1:3,ielem,ja)= dbna(1:3,ielem,ja) &
             +djrij(1:3)*tmp *aexp*bnai**(aexp-1)
!!$        fat(1:3,ia)= fat(1:3,ia) &
!!$             -dirij(1:3)*tmp *aexp*bnai**(aexp-1)
!!$        fat(1:3,ja)= fat(1:3,ja) &
!!$             -djrij(1:3)*tmp *aexp*bnai**(aexp-1)
      enddo
    endif
    bnai= bnai**aexp

  end subroutine bfunc
!=======================================================================
  function func2(rij,ielem,is,js)
!
!  Calculate selected basis function specified by ielem.
!
!  ***If you add a function, you have to change ncnst_type parameter
!     in this module header.
!
    implicit none
    integer,intent(in):: ielem,is,js
    real(8),intent(in):: rij
    real(8):: func2
    real(8):: a(max_ncnst),r2i,r4i

    func2= 0d0
    if( itype(ielem).eq.1 ) then ! Gaussian-type
      a(1:3)= cnst(1:3,ielem)
      func2= rij**a(1) *exp(-a(2)*(rij-a(3))**2)

    elseif( itype(ielem).eq.2 ) then ! cosine-type
      a(1:2)= cnst(1:2,ielem)
      func2= (1d0+cos(rij*a(1))) *rij**nint(a(2))

    elseif( itype(ielem).eq.4 ) then ! poly-1
      a(1)= cnst(1,ielem)
      func2= a(1) /rij
    elseif( itype(ielem).eq.5 ) then ! poly-2
      a(1)= cnst(1,ielem)
      func2= a(1) /rij**2
    elseif( itype(ielem).eq.6 ) then ! poly-4
      a(1)= cnst(1,ielem)
      r2i= 1d0/rij**2
      func2= a(1) *r2i*r2i
    elseif( itype(ielem).eq.7 ) then ! poly-6
      a(1)= cnst(1,ielem)
      r2i= 1d0/rij**2
      func2= a(1) *r2i*r2i*r2i
    elseif( itype(ielem).eq.8 ) then ! poly-8
      a(1)= cnst(1,ielem)
      r2i= 1d0/rij**2
      r4i= r2i*r2i
      func2= a(1) *r4i*r4i
    elseif( itype(ielem).eq.9 ) then ! poly-10
      a(1)= cnst(1,ielem)
      r2i= 1d0/rij**2
      r4i= r2i*r2i
      func2= a(1) *r4i*r4i*r2i
    elseif( itype(ielem).eq.10 ) then ! poly-12
      a(1)= cnst(1,ielem)
      r2i= 1d0/rij**2
      r4i= r2i*r2i
      func2= a(1) *r4i*r4i*r4i

    endif
    return
  end function func2
!=======================================================================
  function dfunc2(rij,ielem,is,js)
    implicit none
    integer,intent(in):: ielem,is,js
    real(8),intent(in):: rij
    real(8):: dfunc2,tmp,a(max_ncnst),ri,r2i,r4i
    integer:: ia2

    dfunc2= 0d0
    if( itype(ielem).eq.1 ) then ! Gaussian-type
      a(1:3)= cnst(1:3,ielem)
      dfunc2= (a(1)*rij**(a(1)-1d0) -2d0*a(2)*(rij-a(3))*rij**a(1) ) &
           *exp(-a(2)*(rij-a(3))**2)
!!$      dfunc2= (a(1) -2d0*a(2)*(rij-a(3))*rij) *rij**(a(1)-1d0) &
!!$           *exp(-a(2)*(rij-a(3))**2)

    elseif( itype(ielem).eq.2 ) then ! cosine-type
      a(1:2)= cnst(1:2,ielem)
      ia2= nint(a(2))
      dfunc2= -a(1)*sin(rij*a(1)) *rij**ia2 &
           +ia2*rij**(ia2-1) *(1d0+cos(rij*a(1)))

    elseif( itype(ielem).eq.4 ) then ! poly-1
      a(1)= cnst(1,ielem)
      r2i= 1d0/rij/rij
      dfunc2= -a(1) *r2i
    elseif( itype(ielem).eq.5 ) then ! poly-2
      a(1)= cnst(1,ielem)
      r2i= 1d0/rij/rij
      dfunc2= -a(1) *2.0*r2i/rij
    elseif( itype(ielem).eq.6 ) then ! poly-4
      a(1)= cnst(1,ielem)
      r2i= 1d0/rij/rij
      dfunc2= -a(1) *4.0*r2i*r2i/rij
    elseif( itype(ielem).eq.7 ) then ! poly-6
      a(1)= cnst(1,ielem)
      r2i= 1d0/rij/rij
      dfunc2= -a(1) *6.0*r2i*r2i*r2i/rij
    elseif( itype(ielem).eq.8 ) then ! poly-8
      a(1)= cnst(1,ielem)
      r2i= 1d0/rij/rij
      r4i= r2i*r2i
      dfunc2= -a(1) *8.0*r4i*r4i/rij
    elseif( itype(ielem).eq.9 ) then ! poly-10
      a(1)= cnst(1,ielem)
      r2i= 1d0/rij/rij
      r4i= r2i*r2i
      dfunc2= -a(1) *10.0*r4i*r4i*r2i/rij
    elseif( itype(ielem).eq.10 ) then ! poly-12
      a(1)= cnst(1,ielem)
      r2i= 1d0/rij/rij
      r4i= r2i*r2i
      dfunc2= -a(1) *10.0*r4i*r4i*r4i/rij
    endif
    return
  end function dfunc2
!=======================================================================
  function func3(rij,rj,rik,rk,ielem,is,js,ks)
    implicit none
    integer,intent(in):: ielem,is,js,ks
    real(8),intent(in):: rij(3),rj,rik(3),rk
    real(8):: func3,cs,a(max_ncnst)

    func3= 0d0
    if( itype(ielem).eq.3 ) then ! angular
      a(1)= cnst(1,ielem)
      cs= dot(rij,rik)/rj/rk
      func3= (a(1)+cs)**2/(abs(a(1))+1d0)**2
    endif

  end function func3
!=======================================================================
  subroutine read_params(myid,mpi_world,rcin)
    implicit none
    include 'mpif.h'

    integer,intent(in):: myid,mpi_world
    real(8),intent(out):: rcin
    integer:: itmp,ierr,i,j
    logical:: lexist

    
!.....read constants at the 1st call
    inquire(file=trim(ccfname),exist=lexist)
    if( .not. lexist ) then
      if( myid.eq.0 ) then
        write(6,'(a)') ' [Error] '//ccfname//' does not exist !!!.'
        write(6,'(a)') '   The linreg potential needs '//ccfname//'.'
      endif
      call mpi_finalize(ierr)
      stop
    endif
    open(51,file=trim(ccfname),status='old')
    read(51,*) nelem,nexp,nsp
    allocate(itype(nelem),cnst(max_ncnst,nelem),exps(nelem) &
         ,icmb(max_natm,nelem))
    do i=1,nelem
      read(51,*) itype(i),exps(i), &
           (icmb(j,i),j=1,natm_type(itype(i))), &
           (cnst(j,i),j=1,ncnst_type(itype(i)))
    enddo
    close(51)

!.....read parameters at the 1st call
    inquire(file=trim(cpfname),exist=lexist)
    if( .not. lexist ) then
      if( myid.eq.0 ) then
        write(6,'(a)') ' [Error] '//cpfname//' does not exist !!!.'
        write(6,'(a)') '   The linreg potential needs '//cpfname//'.'
      endif
      call mpi_finalize(ierr)
      stop
    endif
    open(50,file=trim(cpfname),status='old')
    read(50,*) itmp,rcin
    if( itmp.ne.nelem ) then
      if(myid.eq.0) then
        write(6,'(a)') ' [Error] nelems in in.const.linreg ' // &
             'and in.paramlinreg are inconsistent.'
      endif
      call mpi_finalize(ierr)
      stop
    endif
    allocate(coeff(nelem))
    do i=1,nelem
      read(50,*) coeff(i)
    enddo
    close(50)


    return
  end subroutine read_params
!=======================================================================
  subroutine read_params_linreg(myid,mpi_world,iprint)
    use descriptor,only: nsf
    implicit none
    include 'mpif.h'

    integer,intent(in):: myid,mpi_world,iprint
    integer:: itmp,ierr,i,j
    logical:: lexist
    character(len=128):: fname

    if( myid.eq.0 ) then
      if( iprint.gt.1 ) print *,'read_params_linreg...'
      fname = trim(paramsdir)//'/'//trim(cpfname)
!.....read parameters at the 1st call
      inquire(file=trim(fname),exist=lexist)
      if( .not. lexist ) then
        write(6,'(a)') ' [Error] '//trim(cpfname)//' does not exist !!!.'
        write(6,'(a)') '   The linreg potential needs '//trim(cpfname)//'.'
!!$        call mpi_finalize(ierr)
        stop
      endif
      open(50,file=trim(cpfname),status='old')
      read(50,*) itmp
      if( itmp.ne.nsf ) then
        write(6,'(a)') ' [Error] nsf in in.params.linreg ' // &
             'and in.params.desc are inconsistent.'
!!$        call mpi_finalize(ierr)
        stop
      endif
    endif

    allocate(coeff(nsf))
    if( myid.eq.0 ) then
      do i=1,nsf
        read(50,*) coeff(i)
      enddo
      close(50)
    endif
    call mpi_bcast(coeff,nsf,mpi_real8,0,mpi_world,ierr)
    

    return
  end subroutine read_params_linreg
!=======================================================================
  subroutine set_paramsdir_linreg(dname)
!
!  Accessor routine to set paramsdir.
!
    implicit none
    character(len=*),intent(in):: dname

    paramsdir = trim(dname)
    return
  end subroutine set_paramsdir_linreg
!=======================================================================
  subroutine set_params_linreg(ndimp,params_in)
!
!  Accesor routine to set linreg parameters from outside.
!  It is supposed to be called from fitpot in a seriral process.
!
    integer,intent(in):: ndimp
    real(8),intent(in):: params_in(ndimp)

    if( .not. allocated(coeff) ) then
      allocate(coeff(ndimp))
    endif

    coeff(1:ndimp) = params_in(1:ndimp)
    lprmset_linreg = .true.
    return
  end subroutine set_params_linreg
!=======================================================================
  subroutine gradw_linreg(namax,natm,tag,ra,nnmax,h,rcin,lspr &
       ,iprint,ndimp,gwe,gwf,gws,lematch,lfmatch,lsmatch,iprm0)
!
!  Derivative of linreg pot w.r.t parameters.
!  - iprm0: The starting point -1 in parameter array for this FF.
!
    use descriptor,only: gsfi,dgsfi,igsfi,nsf,calc_desci,make_gsf_arrays, &
         pre_desci
    use util,only: itotOf
    implicit none
    integer,intent(in):: namax,natm,nnmax,ndimp,iprint,lspr(0:nnmax,namax)&
         ,iprm0
    real(8),intent(in):: tag(namax),ra(3,namax),h(3,3),rcin
    real(8),intent(inout):: gwe(ndimp),gwf(3,ndimp,natm),gws(6,ndimp)
    logical,intent(in):: lematch,lfmatch,lsmatch

    integer:: i,ia,ja,jj,isf,ne,nf,jra
!!$    integer,external:: itotOf
    real(8):: ftmp(3),xi(3),xj(3),xij(3),rij(3)

    call pre_desci(namax,natm,nnmax,lspr,iprint,rcin)
    
    do ia=1,natm
      call calc_desci(ia,namax,natm,nnmax,h,tag,ra,lspr,rcin,iprint)
      if( lematch ) then
        ne = iprm0
        do isf=1,nsf
          ne = ne + 1
          gwe(ne) = gwe(ne) +gsfi(isf)
        enddo ! isf=...
      endif ! lematch
    
      if( lfmatch .or. lsmatch ) then
!!$!.....ja != ia
!!$        do jj=0,lspr(0,ia)
!!$          if( jj.eq.0 ) then
!!$            ja = ia
!!$          else
!!$            ja = lspr(jj,ia)
!!$          endif
!!$          jra = itotOf(tag(ja))
!!$          nf = iprm0
!!$          do isf=1,nsf
!!$            nf = nf + 1
!!$            gwf(1:3,nf,jra) = gwf(1:3,nf,jra) -dgsf(1:3,isf,jj,ia)
!!$          enddo
!!$        enddo
    
        xi(1:3) = ra(1:3,ia)
!.....ja != ia
        do jj=0,lspr(0,ia)
          if( jj.eq.0 ) then
            ja = ia
          else
            ja = lspr(jj,ia)
            xj(1:3) = ra(1:3,ja)
            xij(1:3) = xj(1:3) -xi(1:3)
            rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
          endif
          jra = itotOf(tag(ja))
          if( jj.eq.0 ) then
            nf = iprm0
            do isf=1,nsf
              nf = nf + 1
              ftmp(1:3) = -dgsfi(1:3,isf,jj)
!.....Force
              gwf(1:3,nf,jra) = gwf(1:3,nf,jra) +ftmp(1:3)
!.....No stress contribution for jj==0
            enddo
          else
            nf = iprm0
            do isf=1,nsf
              nf = nf + 1
              ftmp(1:3) = -dgsfi(1:3,isf,jj)
!.....Force
              gwf(1:3,nf,jra) = gwf(1:3,nf,jra) +ftmp(1:3)
!.....Stress
              gws(1,nf) = gws(1,nf) +rij(1)*ftmp(1)
              gws(2,nf) = gws(2,nf) +rij(2)*ftmp(2)
              gws(3,nf) = gws(3,nf) +rij(3)*ftmp(3)
              gws(4,nf) = gws(4,nf) +rij(2)*ftmp(3)
              gws(5,nf) = gws(5,nf) +rij(1)*ftmp(3)
              gws(6,nf) = gws(6,nf) +rij(1)*ftmp(2)
            enddo
          endif
        enddo
      endif ! lfmatch or lsmatch
    enddo ! ia=...
    return
  end subroutine gradw_linreg
!=======================================================================
  subroutine set_iglid_linreg(cpena,cfmethod)
!
!  Initialize some only required for fitpot.
!
    use descriptor,only: ngl,glval,iglid,nsf
    character(len=*),intent(in):: cpena,cfmethod

    integer:: i,isf
    
!.....Make groups for group LASSO/FS
    if( trim(cpena).eq.'glasso' &
         .or. trim(cfmethod).eq.'gfs') then
      if( .not.allocated(iglid) ) allocate(iglid(nsf))
      iglid(1:nsf)= 0
      i= 0
      do isf=1,nsf
        iglid(isf)= isf
      enddo
    endif
    
  end subroutine set_iglid_linreg
end module linreg
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
