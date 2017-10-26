module Bonny_WRe
!-----------------------------------------------------------------------
!                     Last modified: <2017-10-24 14:50:58 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
!  Parallel implementation of EAM poetntial of Bonney et al.
!  See G. Bonny et al., J. Appl. Phys. 121, 165107 (2017).
!-----------------------------------------------------------------------
  implicit none

!.....Max num of species
  integer,parameter:: msp = 9
!.....Number of species
  integer:: nsp = 0

!!$.....Default parameter set for Al
!!$  real(8):: ea_a  = 0.763905d0
!!$  real(8):: ea_b  = 0.075016d0
!!$  real(8):: ea_c  = 0.159472d0
!!$  real(8):: ea_re = 3.389513d-10 /ang
!!$  real(8):: ea_al = 1.755162d+10 *ang
!!$  real(8):: ea_bt = 2.003449d+10 *ang
!!$  real(8):: ea_xi = 0.147699d0
!!$  real(8):: am_al = 26.9815d0

  real(8),allocatable:: ea_a(:),ea_xi(:)
  real(8),allocatable,dimension(:,:):: ea_b,ea_c,ea_re, &
       ea_alp,ea_beta,ea_rc
  logical,allocatable:: interact(:,:)

  logical:: lprmset_EAM = .false.

!.....Parameters from fitpot
  integer:: nprms
  real(8),allocatable:: prms(:)

!.....Types of forms of potential terms
  character(len=128),allocatable:: type_rho(:,:),type_frho(:),type_phi(:,:)
  character(len=128):: default_type_rho = 'exp1'
  character(len=128):: default_type_frho = 'sqrt1'
  character(len=128):: default_type_phi = 'SM'

contains
  subroutine force_Bonny_WRe(namax,natm,tag,ra,nnmax,aa,strs,h,hi,tcom &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint,l1st)
    implicit none
    include "mpif.h"
    include "./params_unit.h"
    integer,intent(in):: namax,natm,nnmax,nismax,lspr(0:nnmax,namax)&
         ,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_md_world,myid_md,nex(3)
    real(8),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),sv(3,6) &
         ,acon(nismax),rc,tag(namax)
    real(8),intent(inout):: tcom
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: l1st
    logical:: lstrs

    integer:: i,j,k,l,m,n,ierr,is,js,ixyz,jxyz
    real(8):: xij(3),rij,rcij,dfi,dfj,drho,drdxi(3),drdxj(3),r,at(3)
    real(8):: x,y,z,xi(3),epotl,epott,tmp,dtmp
    real(8),allocatable,save:: rho(:)
    real(8),allocatable,save:: strsl(:,:,:)

    if( l1st ) then
!.....Check cutoff radius
      if( allocated(rho) ) deallocate(rho)
      allocate(rho(namax+nbmax))
      if( allocated(strsl) ) deallocate(strsl)
      allocate(strsl(3,3,namax))
      do is=1,msp
        do js=is,msp
          if( .not.interact(is,js) ) cycle
          if( rc.lt.ea_rc(is,js) ) then
            if( myid_md.eq.0 ) then
              write(6,*) ' Error: rc is smaller than one of EAM rcs'
            endif
            call mpi_finalize(ierr)
            stop
          endif
        enddo
      enddo
    endif

    if( size(strsl).lt.3*3*namax ) then
      deallocate(strsl)
      allocate(strsl(3,3,namax))
    endif
    if( size(rho).lt.namax+nbmax ) then
      deallocate(rho)
      allocate(rho(namax+nbmax))
    endif
    
    epotl= 0d0
    rho(1:namax)= 0d0
    strsl(1:3,1:3,1:namax) = 0d0

!-----rho(i)
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is = int(tag(i))
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
        js = int(tag(j))
        x= ra(1,j) -xi(1)
        y= ra(2,j) -xi(2)
        z= ra(3,j) -xi(3)
        xij(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
        rij=sqrt(xij(1)*xij(1)+ xij(2)*xij(2) +xij(3)*xij(3))
        rcij= ea_rc(is,js)
        if( rij.gt.rcij ) cycle
!!$        rho(i)= rho(i) +exp(-ea_beta(is,js)*(rij-ea_re(is,js))) &
!!$             *fcut1(rij,rcij)
        
        rho(i) = rho(i) +rhoij(rij,is,js,type_rho(is,js))
      enddo
!!$      rho(i)= dsqrt(rho(i))
    enddo

!-----copy rho of boundary atoms
    call copy_dba_fwd(tcom,namax,natm,nb,nbmax,lsb,nex,&
         lsrc,myparity,nn,sv,mpi_md_world,rho,1)

!-----dE/dr_i
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is = int(tag(i))
!!$      dfi= -0.5d0*ea_a(is)/rho(i)
      dfi = dfrho(is,rho(i),type_frho(is))
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
        if(j.le.i) cycle
        js = int(tag(j))
        x= ra(1,j) -xi(1)
        y= ra(2,j) -xi(2)
        z= ra(3,j) -xi(3)
        xij(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
        rij=sqrt(xij(1)**2+ xij(2)**2 +xij(3)**2)
        rcij = ea_rc(is,js)
        if( rij.gt.rcij ) cycle
        drdxi(1:3)= -xij(1:3)/rij
!!$        r= rij -ea_re(is,js)
!.....2-body term
!!$        tmp= 2d0*ea_b(is,js)*exp(-0.5d0*ea_beta(is,js)*r) &
!!$             -ea_c(is,js)*(1d0+ea_alp(is,js)*r)*exp(-ea_alp(is,js)*r)
!!$        tmp2 = tmp *0.5d0 *fcut1(rij,rcij)
        tmp = 0.5d0 *phi(rij,rcij,is,js,type_phi(is,js))
        epi(i)= epi(i) +tmp
        epi(j)= epi(j) +tmp
        if(j.le.natm) then
          epotl=epotl +tmp +tmp
        else
          epotl=epotl +tmp
        endif
!!$        dphi= -ea_beta(is,js)*ea_b(is,js)*exp(-0.5d0*ea_beta(is,js)*r)*fcut1(rij,rcij)  &
!!$             +ea_c(is,js)*ea_alp(is,js)*ea_alp(is,js)*r &
!!$             *exp(-ea_alp(is,js)*r) *fcut1(rij,rcij) &
!!$             +tmp*dfcut1(rij,rcij)
        dtmp = dphi(rij,rcij,is,js,type_phi(is,js))
        aa(1:3,i)=aa(1:3,i) -dtmp*drdxi(1:3)
        aa(1:3,j)=aa(1:3,j) +dtmp*drdxi(1:3)
!.....Atomic stress for 2-body terms
        if( lstrs ) then
          do ixyz=1,3
            do jxyz=1,3
              strsl(jxyz,ixyz,i)=strsl(jxyz,ixyz,i) &
                   -0.5d0*dtmp*xij(ixyz)*(-drdxi(jxyz))
              strsl(jxyz,ixyz,j)=strsl(jxyz,ixyz,j) &
                   -0.5d0*dtmp*xij(ixyz)*(-drdxi(jxyz))
            enddo
          enddo
        endif
!.....Embedded term
!!$        drhoij= -ea_beta(is,js)*exp(-ea_beta(is,js)*r)*fcut1(rij,rcij) &
!!$             +exp(-ea_beta(is,js)*r)*dfcut1(rij,rcij)
        drho = drhoij(rij,rcij,is,js,type_rho(is,js))
!!$        dfj= -0.5d0*ea_a(js)/rho(j)
        dfj = dfrho(js,rho(j),type_frho(js))
        tmp = (dfi+dfj)*drho
        aa(1:3,i)=aa(1:3,i) -tmp*drdxi(1:3)
        aa(1:3,j)=aa(1:3,j) +tmp*drdxi(1:3)
!.....Atomic stress of many-body contributions
        if( lstrs ) then
          do ixyz=1,3
            do jxyz=1,3
              strsl(jxyz,ixyz,i)=strsl(jxyz,ixyz,i) &
                   -0.5d0*tmp*xij(ixyz)*(-drdxi(jxyz))
              strsl(jxyz,ixyz,j)=strsl(jxyz,ixyz,j) &
                   -0.5d0*tmp*xij(ixyz)*(-drdxi(jxyz))
            enddo
          enddo
        endif
      enddo
      tmp = frho(is,rho(i),type_frho(is))
      epi(i)=epi(i) +tmp
      epotl=epotl +tmp
    enddo

    if( lstrs ) then
!!$      call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
!!$           ,nn,mpi_md_world,strsl,9)
      strs(1:3,1:3,1:natm)= strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)
    endif

!-----gather epot
    call mpi_allreduce(epotl,epott,1,mpi_real8 &
         ,mpi_sum,mpi_md_world,ierr)
    epot= epot +epott

  end subroutine force_Bonny_WRe
!=======================================================================
  function rhoij(rij,js)
!
! Calculate rho of atom j at distance rij.
!
    real(8),intent(in):: rij
    integer,intent(in):: js
    real(8):: rhoij

    rhoij = phi(r)
    if( js.eq.1 ) then  ! Only in case of W, scaling with S
      rhoij = rhoij *bonny_gauge_S
    endif
    return
  end function rhoij
!=======================================================================
  function drhoij(rij,js)
    implicit none
    real(8),intent(in):: rij
    integer,intent(in):: js
    real(8):: drhoij
    drhoij = dphi(r)
    if( js.eq.1 ) then
      drhoij = drhoij *bonny_gauge_S
    endif
    return
  end function drhoij
!=======================================================================
  function phi(r)
    implicit none
    real(8),intent(in):: r
    real(8):: phi
    real(8),external:: hvsd
    phi = bonny_c0 *(bonny_rc -r)**3 *hvsd(bonny_rc -r)
    return
  end function phi
!=======================================================================
  function dphi(r)
    implicit none
    real(8),intent(in):: r
    real(8):: dphi
    real(8),external:: hsvd
    dphi = -3d0 * bonny_c0 *(bonny_rc -r)**2 *hvsd(bonny_rc -r)
  end function dphi
!=======================================================================
  function frho(rho)
    implicit none
    real(8),intent(in):: rho
    real(8):: frho
    
  end function frho
!=======================================================================
  function veq(r)
    implicit none
    real(8),intent(in):: r
    
  end function veq
!=======================================================================
  function xi(x)
    implicit none
    real(8),intent(in):: x
    real(8):: xi

    xi= 0.1818d0*exp(-3.2d0*x) &
         +0.5099d0*exp(-0.9423d0*x) &
         +0.2802d0*exp(-0.4029d0*x) &
         +0.02817d0*exp(-0.2016d0*x)
    return
  end function xi
!=======================================================================
  function dxi(x)
    implicit none
    real(8),intent(in):: x
    real(8):: dxi

    dxi= -0.58176d0*exp(-3.2d0*x) &
         -0.48047877d0*exp(-0.9423d0*x) &
         -0.11289258d0*exp(-0.4029d0*x) &
         -0.005679072d0*exp(-0.2016d0*x)
    return
  end function dxi
!=======================================================================
  function zeta(x)
    implicit none
    real(8),intent(in):: x
    real(8):: zeta

    zeta = (3d0*x**5 -10d0*x**3 +15d0*x +8d0)/16d0
    return
  end function zeta
!=======================================================================
  function dzeta(x)
    implicit none
    real(8),intent(in):: x
    real(8):: dzeta

    dzeta = (15d0*x**4 -30d0*x**2 +15d0)/16d0
    return
  end function dzeta
end module Bonny_WRe
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
