module linreg
!-----------------------------------------------------------------------
!                     Last modified: <2017-01-06 09:33:21 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
!  Parallel implementation of linear regression potential for pmd
!    - 2014.06.11 by R.K. 1st implementation
!    - 2015.02.03 by R.K. extended to multiple species
!-----------------------------------------------------------------------
!.....parameter file name
  character(128),parameter:: cpfname= 'in.params.linreg'
  character(128),parameter:: ccfname='in.const.linreg'
  character(128),parameter:: cmbfname='in.comb.linreg'
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
  subroutine force_linreg(namax,natm,tag,ra,nnmax,aa,strs,h,hi,tcom &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_world,myid,epi,epot,nismax,acon,lstrs,iprint)
    implicit none
    include "mpif.h"
    include "./params_unit.h"
    integer,intent(in):: namax,natm,nnmax,nismax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_world,myid,lspr(0:nnmax,namax),nex(3)
    real(8),intent(in):: ra(3,namax),tag(namax),acon(nismax) &
         ,h(3,3),hi(3,3),sv(3,6)
    real(8),intent(inout):: tcom,rc
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical:: lstrs

!.....local
    integer:: i,j,k,l,m,n,ixyz,jxyz,is,js,ks,ierr,nbl,ia,ielem &
         ,iwgt
    real(8):: rcin,b_na,at(3),epotl,wgt,aexp,bnai,apot
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
    epi(1:natm+nb)= 0d0
    strs(1:3,1:3,1:namax)= 0d0
    aa(1:3,1:namax)= 0d0

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
    call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,dbna,3*nelem)
    do ia=1,natm
      do ielem=1,nelem
        write(81,'(2i10,3es23.14e3)') ia,ielem,dbna(1:3,ielem,ia)
      enddo
    enddo
    close(80)
    close(81)
#endif

    call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,aa,3)
!-----reduced force
    do i=1,natm
      at(1:3)= aa(1:3,i)
      aa(1:3,i)= hi(1:3,1)*at(1) +hi(1:3,2)*at(2) +hi(1:3,3)*at(3)
    enddo
!-----multiply 0.5d0*dt**2/am(i)
    do i=1,natm
      is= int(tag(i))
      aa(1:3,i)= acon(is)*aa(1:3,i)
    enddo

!-----gather epot
    epot= 0d0
    call mpi_allreduce(epotl,epot,1,mpi_double_precision &
         ,mpi_sum,mpi_world,ierr)
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
    real(8),external:: sprod

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
            cs= sprod(rij,rik)/rj/rk
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
    real(8),external:: sprod 

    func3= 0d0
    if( itype(ielem).eq.3 ) then ! angular
      a(1)= cnst(1,ielem)
      cs= sprod(rij,rik)/rj/rk
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

end module linreg
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
