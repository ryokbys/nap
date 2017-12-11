module AFS_W
contains
  subroutine force_AFS_W(namax,natm,tag,ra,nnmax,aa,strs,h,hi,tcom &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
!-----------------------------------------------------------------------
!  Parallel implementation of Ackland potential for W based on 
!  Finnis and Sinclair potential.
!    - smoothing is applied to 2-body potential
!    - rho of boundary atoms are sent to the neighbor nodes
!    - only force on i is calculated, not necessary to send-back
!-----------------------------------------------------------------------
!  See
!    - Philo. Mag. A, vol.56(1), 1987, pp.15--30
!    - Philo. Mag. A, vol.50(1), 1984, pp.45--55
!-----------------------------------------------------------------------
    use force, only: copy_dba_fwd, copy_dba_bk
    implicit none
    include "mpif.h"
    include "./params_unit.h"
    include "params_AFS_W.h"
    integer,intent(in):: namax,natm,nnmax,nismax,lspr(0:nnmax,namax)&
         ,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_md_world,myid_md,nex(3)
    real(8),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),sv(3,6) &
         ,acon(nismax),rc,tag(namax)
    real(8),intent(inout):: tcom
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical:: lstrs

    integer:: i,j,k,l,m,n,ierr,is,js,ixyz,jxyz
    real(8):: xij(3),rij,dfi,dfj,drhoij,drdxi(3),drdxj(3),r,at(3)
    real(8):: x,y,z,xi(3),epotl,v2,dv2,dphi,dphj,tmp
    logical,save:: l1st=.true.
    real(8),allocatable,save:: sqrho(:)

    if( l1st ) then
      allocate(sqrho(namax+nbmax))
      l1st=.false.
    endif

    epotl= 0d0
    sqrho(1:natm)= 0d0

!-----rho(i)
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is= int(tag(i))
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
        js= int(tag(j))
        x= ra(1,j) -xi(1)
        y= ra(2,j) -xi(2)
        z= ra(3,j) -xi(3)
        xij(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
        rij=sqrt(xij(1)*xij(1)+ xij(2)*xij(2) +xij(3)*xij(3))
        sqrho(i)= sqrho(i) +phi_IWHe(rij,is,js)
      enddo
      sqrho(i)= dsqrt(sqrho(i))
    enddo

    call copy_dba_fwd(tcom,namax,natm,nb,nbmax,lsb,nex,&
         lsrc,myparity,nn,sv,mpi_md_world,sqrho,1)
!!$    if( myid_md.ge.0 ) then
!!$!-----copy rho of boundary atoms
!!$      call copy_rho_ba(tcom,namax,natm,nb,nbmax,lsb,lsrc,myparity,nn,sv &
!!$           ,mpi_md_world,sqrho)
!!$    else
!!$      call distribute_dba(natm,namax,tag,sqrho,1)
!!$    endif

    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is= int(tag(i))
      dfi= -0.5d0*p_W_A/sqrho(i)
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
        if(j.le.i) cycle
        js= int(tag(j))
        x= ra(1,j) -xi(1)
        y= ra(2,j) -xi(2)
        z= ra(3,j) -xi(3)
        xij(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
        rij=sqrt(xij(1)**2+ xij(2)**2 +xij(3)**2)
        if( rij.gt.rc ) cycle
        drdxi(1:3)= -xij(1:3)/rij
!.....2-body term
        v2= 0.5d0 *v2_IWHe(rij,is,js)
        dv2= dv2_IWHe(rij,is,js)
        epi(i)= epi(i) +v2
        epi(j)= epi(j) +v2
        if(j.le.natm) then
          epotl=epotl +v2 +v2
        else
          epotl=epotl +v2
        endif
        aa(1:3,i)=aa(1:3,i) -dv2*drdxi(1:3)
        aa(1:3,j)=aa(1:3,j) +dv2*drdxi(1:3)
!.....atomic stress for 2-body terms
        do ixyz=1,3
          do jxyz=1,3
            strs(jxyz,ixyz,i)=strs(jxyz,ixyz,i) &
                 -0.5d0*dv2*xij(ixyz)*(-drdxi(jxyz))
            strs(jxyz,ixyz,j)=strs(jxyz,ixyz,j) &
                 -0.5d0*dv2*xij(ixyz)*(-drdxi(jxyz))
          enddo
        enddo
!.....N-body term
        dfj= -0.5d0*p_W_A/sqrho(j)
        dphi= dphi_IWHe(rij,is,js) !/2
        dphj= dphi_IWHe(rij,js,is) !/2
        tmp= (dfi*dphi+dfj*dphj)
        aa(1:3,i)=aa(1:3,i) -tmp*drdxi(1:3)
        aa(1:3,j)=aa(1:3,j) +tmp*drdxi(1:3)
!.....atomic stress of many-body contributions
        do ixyz=1,3
          do jxyz=1,3
            strs(jxyz,ixyz,i)=strs(jxyz,ixyz,i) &
                 -0.5d0*tmp*xij(ixyz)*(-drdxi(jxyz))
            strs(jxyz,ixyz,j)=strs(jxyz,ixyz,j) &
                 -0.5d0*tmp*xij(ixyz)*(-drdxi(jxyz))
          enddo
        enddo
      enddo
      epi(i)=epi(i) -p_W_A*sqrho(i)
      epotl=epotl -p_W_A*sqrho(i)
    enddo

    call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_md_world,strs,9)
!!$    if( myid_md.ge.0 ) then
!!$!-----copy strs of boundary atoms
!!$      call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,lsrc,myparity &
!!$           ,nn,mpi_world,strs,9)
!!$    else
!!$      call reduce_dba_bk(natm,namax,tag,strs,9)
!!$    endif

!-----gather epot
    epot= 0d0
    if( myid_md.ge.0 ) then
      call mpi_allreduce(epotl,epot,1,MPI_DOUBLE_PRECISION &
           ,MPI_SUM,mpi_md_world,ierr)
    else
      epot= epotl
    endif

!      deallocate(sqrho)
  end subroutine force_AFS_W
!=======================================================================
  function v2_IWHe(r,is,js)
!
!  Two-body potential energy
!
    implicit none
    include "./params_unit.h"
    include "params_AFS_W.h"
    real(8),intent(in):: r
    integer,intent(in):: is,js
    real(8):: v2_IWHe,x

    v2_IWHe= 0d0

!.....W-W
    if( r.lt.p_WW_c) then
      v2_IWHe= (r-p_WW_c)*(r-p_WW_c)*(p_WW_c0 +p_WW_c1*r &
           +p_WW_c2*r*r)
    endif
!      if( r.lt.p_WW_b0 ) then
!        v2_IWHe=v2_IWHe
!     &       +p_WW_B*(p_WW_b0-r)**3 *exp(-p_WW_alpha*r)
!      endif

    return
  end function v2_IWHe
!=======================================================================
  function dv2_IWHe(r,is,js)
!
!  Derivative of two-body potential
!
    implicit none
    include "./params_unit.h"
    include "params_AFS_W.h"
    real(8),intent(in):: r
    integer,intent(in):: is,js
    real(8):: dv2_IWHe,ri,x,exar

    dv2_IWHe= 0d0

!.....W-W
    if( r.lt.p_WW_c) then
      dv2_IWHe= 2d0*(r-p_WW_c)*(p_WW_c0 +p_WW_c1*r +p_WW_c2*r*r) &
           +(r-p_WW_c)**2*(p_WW_c1 +2d0*p_WW_c2*r)
    endif
    if( r.lt.p_WW_b0 ) then
      dv2_IWHe=dv2_IWHe &
           +p_WW_B*(p_WW_b0-r)**2 *exp(-p_WW_alpha*r) &
           *(-3d0 -p_WW_alpha*(p_WW_b0-r))
    endif

    return
  end function dv2_IWHe
!=======================================================================
  function phi_IWHe(r,is,js)
!
!  Phi for many-body potential
!
    implicit none
    include "./params_unit.h"
    include "params_AFS_W.h"
    real(8),intent(in):: r
    integer,intent(in):: is,js
    real(8):: phi_IWHe

    phi_IWHe= 0d0

!.....phi_W
    if( r.le.p_W_d ) then
      phi_IWHe= (r-p_W_d)**2
    endif

    return
  end function phi_IWHe
!=======================================================================
  function dphi_IWHe(r,is,js)
    implicit none
    include "./params_unit.h"
    include "params_AFS_W.h"
    real(8),intent(in):: r
    integer,intent(in):: is,js
    real(8):: dphi_IWHe

    dphi_IWHe= 0d0

!.....phi_W
    if( r.le.p_W_d ) then
      dphi_IWHe= 2d0*(r-p_w_d)
    endif

    return
  end function dphi_IWHe
!=======================================================================
end module AFS_W
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
