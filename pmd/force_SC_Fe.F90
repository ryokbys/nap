module SC_Fe
contains
  subroutine force_SC_Fe(namax,natm,tag,ra,nnmax,aa,strs,h,hi,tcom &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
!-----------------------------------------------------------------------
!  Parallel implementation of EAM Sutton-Chen model for Fe (iron).
!    - See PRB 73, 224113 (2006), L.Koci et al.
!    - rho of boundary atoms are sent to the neighbor nodes
!    - only force on i is calculated, not necessary to send-back
!-----------------------------------------------------------------------
    use force, only: copy_dba_fwd
    implicit none
    include "mpif.h"
    include "./params_unit.h"
    include "params_SC_Fe.h"
    integer,intent(in):: namax,natm,nnmax,nismax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_md_world,myid_md,nex(3)
    integer,intent(in):: lspr(0:nnmax,namax)
    real(8),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),sv(3,6) &
         ,acon(nismax),rc,tag(namax)
    real(8),intent(inout):: tcom
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical:: lstrs

    integer:: i,j,k,l,m,n,ierr,is
    real(8):: xij(3),rij,dfi,dfj,drhoij,drdxi(3),drdxj(3),at(3)
    real(8):: x,y,z,xi(3),epotl,epott,phic,dphic,phi,dphi,tmp

    logical,save:: l1st=.true.
    real(8),allocatable,save:: sqrho(:)
    real(8),save:: rhoc,drhoc

    if( l1st ) then
      allocate(sqrho(namax+nbmax))
!.....smoothing embeded term
      rhoc = (sc_a/rc)**sc_m
      drhoc= -sc_m*rhoc/rc
      l1st=.false.
    endif

    epotl= 0d0
    sqrho(1:natm)= 0d0

!c-----make pair list for 2-body term
!      call mk_lspr(namax,natm,nb,nnmax,tag,ra,rc,h,hi
!     &     ,anxi,anyi,anzi,lspr)

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
        sqrho(i)= sqrho(i) +(sc_a/rij)**sc_m &
             -rhoc -(rij-rc)*drhoc
      enddo
      sqrho(i)= dsqrt(sqrho(i))
    enddo

!.....copy rho of boundary atoms
    call copy_dba_fwd(tcom,namax,natm,nb,nbmax,lsb,nex,&
         lsrc,myparity,nn,sv,mpi_md_world,sqrho,1)
!!$    call copy_rho_ba(tcom,namax,natm,nb,nbmax,lsb &
!!$         ,lsrc,myparity,nn,sv,mpi_md_world,sqrho)

!-----smoothing 2-body term
    phic = sc_eps*(sc_a/rc)**sc_n
    dphic= -sc_n*phic/rc

!-----dE/dr_i
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      dfi= -0.5d0*sc_eps*sc_c/sqrho(i)
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
        if(j.le.i) cycle
        x= ra(1,j) -xi(1)
        y= ra(2,j) -xi(2)
        z= ra(3,j) -xi(3)
        xij(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
        rij=sqrt(xij(1)**2+ xij(2)**2 +xij(3)**2)
        drdxi(1:3)= -xij(1:3)/rij
!          drdxj(1:3)=  xij(1:3)/rij
!---------2-body term
        phi= sc_eps*(sc_a/rij)**sc_n
        tmp= 0.5d0 *( phi -phic -(rij-rc)*dphic )
        epi(i)= epi(i) +tmp
        epi(j)= epi(j) +tmp
        if(j.le.natm) then
          epotl=epotl +tmp +tmp
        else
          epotl=epotl +tmp
        endif
        dphi= -sc_n*phi/rij -dphic
        aa(1:3,i)=aa(1:3,i) -dphi*drdxi(1:3)
        aa(1:3,j)=aa(1:3,j) +dphi*drdxi(1:3)
!---------embedded term
        drhoij= -sc_m*(sc_a/rij)**sc_m /rij -drhoc
        dfj= -0.5d0 *sc_eps*sc_c/sqrho(j)
        aa(1:3,i)=aa(1:3,i) -(dfi+dfj)*drhoij*drdxi(1:3)
        aa(1:3,j)=aa(1:3,j) +(dfi+dfj)*drhoij*drdxi(1:3)
      enddo
      epi(i)=epi(i) -sc_eps*sc_c*sqrho(i)
      epotl=epotl -sc_eps*sc_c*sqrho(i)
    enddo

!-----gather epot
    call mpi_allreduce(epotl,epott,1,MPI_DOUBLE_PRECISION &
         ,MPI_SUM,mpi_md_world,ierr)
    epot= epot +epott

  end subroutine force_SC_Fe
end module SC_Fe
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
