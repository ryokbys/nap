module SM_Al
contains
  subroutine force_SM_Al(namax,natm,tag,ra,nnmax,aa,strs,h,hi,tcom &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,lstrs,iprint)
!-----------------------------------------------------------------------
!  Parallel implementation of EAM force calculation for Al
!    - smoothing is applied to both 2- and many-body terms
!    - rho of boundary atoms are sent to the neighbor nodes
!    - only force on i is calculated, not necessary to send-back
!-----------------------------------------------------------------------
!  See Streitz and Mintmire, PRB 50(16), 11996 (1994).
!  This is implementation of only Al-Al potential in the literature,
!  and non-variable charge implementation.
!-----------------------------------------------------------------------
    implicit none
    include "mpif.h"
    include "./params_unit.h"
    include "params_SM_Al.h"
    integer,intent(in):: namax,natm,nnmax,nismax,lspr(0:nnmax,namax)&
         ,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_md_world,myid_md,nex(3)
    real(8),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),sv(3,6) &
         ,rc,tag(namax)
    real(8),intent(inout):: tcom
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical:: lstrs

    integer:: i,j,k,l,m,n,ierr,is
    real(8):: xij(3),rij,dfi,dfj,drhoij,drdxi(3),drdxj(3),r,dphi,at(3)
    real(8):: x,y,z,xi(3),epotl,epott,tmp
    logical,save:: l1st=.true.
    real(8),allocatable,save:: sqrho(:)
    real(8),save:: exrc,phic,dphic

    if( l1st ) then
      allocate(sqrho(namax+nbmax))
!-----smoothing embeded term
      exrc= exp(-ea_bt*(rc-ea_re))
!-----smoothing 2-body term
      r= rc -ea_re
      phic= 2d0*ea_b*exp(-0.5d0*ea_bt*r) &
           -ea_c*(1d0+ea_al*r)*exp(-ea_al*r)
      dphic= -ea_bt*ea_b*exp(-0.5d0*ea_bt*r)  &
           +ea_c*ea_al*ea_al*r*exp(-ea_al*r)
      l1st=.false.
    endif

    epotl= 0d0
    sqrho(1:natm)= 0d0

!-----rho(i)
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
        x= ra(1,j) -xi(1)
        y= ra(2,j) -xi(2)
        z= ra(3,j) -xi(3)
        xij(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
        rij=sqrt(xij(1)*xij(1)+ xij(2)*xij(2) +xij(3)*xij(3))
        if( rij.gt.rc ) cycle
        sqrho(i)=sqrho(i) +exp(-ea_bt*(rij-ea_re)) &
             -exrc -(rij-rc)*(-ea_bt)*exrc
      enddo
      sqrho(i)= dsqrt(sqrho(i))
    enddo

!-----copy rho of boundary atoms
    call copy_dba_fwd(tcom,namax,natm,nb,nbmax,lsb,nex,&
         lsrc,myparity,nn,sv,mpi_md_world,sqrho,1)
!!$    call copy_rho_ba(tcom,namax,natm,nb,nbmax,lsb,lsrc,myparity,nn,sv &
!!$         ,mpi_md_world,sqrho)

!-----dE/dr_i
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      dfi= -0.5d0*ea_a/sqrho(i)
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
        if(j.le.i) cycle
        x= ra(1,j) -xi(1)
        y= ra(2,j) -xi(2)
        z= ra(3,j) -xi(3)
        xij(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
        rij=sqrt(xij(1)**2+ xij(2)**2 +xij(3)**2)
        if( rij.gt.rc ) cycle
        drdxi(1:3)= -xij(1:3)/rij
        r= rij -ea_re
!---------2-body term
        tmp= 0.5d0 *( 2d0*ea_b*exp(-0.5d0*ea_bt*r) &
             -ea_c*(1d0+ea_al*r)*exp(-ea_al*r) &
             -phic -(rij-rc)*dphic )
        epi(i)= epi(i) +tmp
        epi(j)= epi(j) +tmp
        if(j.le.natm) then
          epotl=epotl +tmp +tmp
        else
          epotl=epotl +tmp
        endif
        dphi= -ea_bt*ea_b*exp(-0.5d0*ea_bt*r)  &
             +ea_c*ea_al*ea_al*r*exp(-ea_al*r) &
             -dphic
        aa(1:3,i)=aa(1:3,i) -dphi*drdxi(1:3)
        aa(1:3,j)=aa(1:3,j) +dphi*drdxi(1:3)
!---------embedded term
        drhoij= -ea_bt*exp(-ea_bt*r) +ea_bt*exrc
        dfj= -0.5d0*ea_a/sqrho(j)
        aa(1:3,i)=aa(1:3,i) -(dfi+dfj)*drhoij*drdxi(1:3)
        aa(1:3,j)=aa(1:3,j) +(dfi+dfj)*drhoij*drdxi(1:3)
      enddo
      epi(i)=epi(i) -ea_a*sqrho(i)
      epotl=epotl -ea_a*sqrho(i)
    enddo

!-----gather epot
    call mpi_allreduce(epotl,epott,1,MPI_DOUBLE_PRECISION &
         ,MPI_SUM,mpi_md_world,ierr)
    epot= epot +epott

!      deallocate(sqrho)
  end subroutine force_SM_Al
end module SM_Al
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
