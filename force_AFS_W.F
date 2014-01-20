      subroutine set_avol()
c
c  Not applicable to parallel calc
c
      use variables
      include './params_au.h'

      avol= h(1,1,0)*h(2,2,0)*h(3,3,0)/natm
      if(myid_md.eq.0) write(6,'(a,es12.4)') ' avol=',avol*bohr2aa**3

      end subroutine set_avol
c=======================================================================
      subroutine get_force(namax,natm,tag,ra,nnmax,aa,strs,h,hi,tcom
     &     ,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
     &     ,mpi_md_world,myid_md,epi,epot,nismax,acon)
c-----------------------------------------------------------------------
c  Parallel implementation of Ackland potential for W based on 
c  Finnis and Sinclair potential.
c    - smoothing is applied to 2-body potential
c    - rho of boundary atoms are sent to the neighbor nodes
c    - only force on i is calculated, not necessary to send-back
c-----------------------------------------------------------------------
c  See
c    - Philo. Mag. A, vol.56(1), 1987, pp.15--30
c    - Philo. Mag. A, vol.50(1), 1984, pp.45--55
c-----------------------------------------------------------------------
      use variables,only: avol
      implicit none
      include "mpif.h"
      include "./params_au.h"
      include "params_AFS_W.h"
      integer,intent(in):: namax,natm,nnmax,nismax
      integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3)
     &     ,nn(6),mpi_md_world,myid_md
      integer,intent(inout):: lspr(0:nnmax,namax)
      real(8),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),sv(3,6)
     &     ,acon(nismax),rc,tag(namax)
      real(8),intent(inout):: tcom
      real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)

      integer:: i,j,k,l,m,n,ierr,is,js,ixyz,jxyz
      real(8):: xij(3),rij,dfi,dfj,drhoij,drdxi(3),drdxj(3),r,at(3)
      real(8):: x,y,z,xi(3),epotl,v2,dv2,dphi,dphj,tmp
      logical,save:: l1st=.true.
      real(8),allocatable,save:: sqrho(:)
      real(8),external:: v2_IWHe,dv2_IWHe,phi_IWHe,dphi_IWHe

      if( l1st ) then
        allocate(sqrho(namax+nbmax))
c.....set atomic volume
        call set_avol()
        l1st=.false.
      endif

      aa(1:3,1:natm)=0d0
      epi(1:natm)= 0d0
      epotl= 0d0
      sqrho(1:natm)= 0d0
      strs(1:3,1:3,1:natm+nb)= 0d0

c-----rho(i)
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

c-----copy rho of boundary atoms
      call copy_rho_ba(tcom,namax,natm,nb,nbmax,lsb,lsrc,myparity,nn,sv
     &     ,mpi_md_world,sqrho)

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
c.....2-body term
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
c.....atomic stress for 2-body terms
          do ixyz=1,3
            do jxyz=1,3
              strs(jxyz,ixyz,i)=strs(jxyz,ixyz,i)
     &             -0.5d0*dv2*xij(ixyz)*(-drdxi(jxyz))
              strs(jxyz,ixyz,j)=strs(jxyz,ixyz,j)
     &             -0.5d0*dv2*xij(ixyz)*(-drdxi(jxyz))
            enddo
          enddo
c.....N-body term
          dfj= -0.5d0*p_W_A/sqrho(j)
          dphi= dphi_IWHe(rij,is,js) !/2
          dphj= dphi_IWHe(rij,js,is) !/2
          tmp= (dfi*dphi+dfj*dphj)
          aa(1:3,i)=aa(1:3,i) -tmp*drdxi(1:3)
          aa(1:3,j)=aa(1:3,j) +tmp*drdxi(1:3)
c.....atomic stress of many-body contributions
          do ixyz=1,3
            do jxyz=1,3
              strs(jxyz,ixyz,i)=strs(jxyz,ixyz,i)
     &             -0.5d0*tmp*xij(ixyz)*(-drdxi(jxyz))
              strs(jxyz,ixyz,j)=strs(jxyz,ixyz,j)
     &             -0.5d0*tmp*xij(ixyz)*(-drdxi(jxyz))
            enddo
          enddo
        enddo
        epi(i)=epi(i) -p_W_A*sqrho(i)
        epotl=epotl -p_W_A*sqrho(i)
      enddo

c-----copy strs of boundary atoms
      call copy_strs_ba(tcom,namax,natm,nb,nbmax,lsb
     &     ,lsrc,myparity,nn,sv,mpi_md_world,strs)
c-----atomic level stress in [Hartree/Bohr^3] assuming 1 Bohr thick
      do i=1,natm
        strs(1:3,1:3,i)= strs(1:3,1:3,i) /avol
c        write(6,'(i5,9es10.2)') i,strs(1:3,1:3,i)
      enddo

c-----reduced force
      do i=1,natm
        at(1:3)= aa(1:3,i)
        aa(1:3,i)= hi(1:3,1)*at(1) +hi(1:3,2)*at(2) +hi(1:3,3)*at(3)
      enddo
c-----multiply 0.5d0*dt**2/am(i)
      do i=1,natm
        is= int(tag(i))
        aa(1:3,i)= acon(is)*aa(1:3,i)
      enddo

c-----gather epot
      epot= 0d0
      call mpi_allreduce(epotl,epot,1,MPI_DOUBLE_PRECISION
     &     ,MPI_SUM,mpi_md_world,ierr)

c      deallocate(sqrho)
      end subroutine get_force
c=======================================================================
      subroutine copy_rho_ba(tcom,namax,natm,nb,nbmax,lsb
     &     ,lsrc,myparity,nn,sv,mpi_md_world,rho)
c-----------------------------------------------------------------------
c     Exchanges boundary-atom data among neighbor nodes
c-----------------------------------------------------------------------
      implicit none
      include "mpif.h"
      integer:: status(MPI_STATUS_SIZE)
c-----in
      integer,intent(in):: namax,natm,nb,nbmax,mpi_md_world
      integer,intent(in):: lsb(0:nbmax,6),lsrc(6),myparity(3),nn(6)
      real(8),intent(in):: sv(3,6)
c-----out
      real(8),intent(inout):: rho(natm+nb),tcom

c-----locals
      integer:: i,j,k,l,m,n,kd,kdd,ku,inode,nsd,nsd3,nrc,nrc3,nbnew,ierr
      real(8):: tcom1,tcom2
      logical,save:: l1st=.true.
      real(8),allocatable,save:: dbuf(:),dbufr(:)

      if( l1st ) then
        allocate(dbuf(nbmax),dbufr(nbmax))
        l1st=.false.
      endif

      nbnew= 0

c-----loop over z, y, & x directions
      do kd=1,3
        tcom1= mpi_wtime()
        do kdd=-1,0
          ku= 2*kd +kdd
          inode= nn(ku)
c---------num. of to-be-sent particles
          nsd= lsb(0,ku)
c---------num. of to-be-recieved particles
          nrc= lsrc(ku)

c---------exchange x
          do i=1,nsd
            j=lsb(i,ku)
            dbuf(i)= rho(j)
          enddo
          call mespasd(inode,myparity(kd),dbuf,dbufr,nsd,nrc,21
     &         ,mpi_md_world)
          do i=1,nrc
            rho(natm+nbnew+i)= dbufr(i)
          enddo

c---------mpi barrier
          call mpi_barrier(mpi_md_world,ierr)
c---------accumulate num. of boundary particles
c          write(6,'(a,2i8)') "nbnew,nrc=",nbnew,nrc
          nbnew=nbnew +nrc
        enddo
        tcom2= mpi_wtime()
        tcom= tcom +tcom2-tcom1
      enddo

      if(nbnew.ne.nb) then
        write(6,'(a,2i8)') "nbnew,(natm+nb)=",nbnew,natm+nb
        stop "error: nbnew.ne.(natm+nb)!!"
      endif
      
      end subroutine copy_rho_ba
c=======================================================================
      subroutine copy_strs_ba(tcom,namax,natm,nb,nbmax,lsb
     &     ,lsrc,myparity,nn,sv,mpi_md_world,strs)
c-----------------------------------------------------------------------
c  Exchanges boundary-atom data among neighbor nodes
c-----------------------------------------------------------------------
      implicit none
      include "mpif.h"
      integer:: status(MPI_STATUS_SIZE)
c-----in
      integer,intent(in):: namax,natm,nb,nbmax,mpi_md_world
      integer,intent(in):: lsb(0:nbmax,6),lsrc(6),myparity(3),nn(6)
      real(8),intent(in):: sv(3,6)
c-----out
      real(8),intent(inout):: strs(9,natm+nb),tcom

c-----locals
      integer:: i,j,k,l,m,n,kd,kdd,ku,inode,nsd,nrc,nbnew,ierr
      real(8):: tcom1,tcom2
      
      logical,save:: l1st=.true.
      real(8),save,allocatable:: dbuf(:,:),dbufr(:,:)

      if( l1st ) then
        allocate(dbuf(9,nbmax),dbufr(9,nbmax))
        l1st=.false.
      endif

      nbnew= 0

c-----loop over z, y, & x directions
      do kd=1,3
        tcom1= mpi_wtime()
        do kdd=-1,0
          ku= 2*kd +kdd
          inode= nn(ku)
c---------num. of to-be-sent particles
          nsd= lsb(0,ku)
c---------num. of to-be-recieved particles
          nrc= lsrc(ku)

c---------exchange strs
          do i=1,nsd
            j=lsb(i,ku)
            dbuf(1:9,i)= strs(1:9,j)
          enddo
          call mespasd(inode,myparity(kd),dbuf,dbufr,9*nsd,9*nrc,21
     &         ,mpi_md_world)
          do i=1,nrc
            strs(1:9,natm+nbnew+i)= dbufr(1:9,i)
          enddo

c---------mpi barrier
          call mpi_barrier(mpi_md_world,ierr)
c---------accumulate num. of boundary particles
c          write(6,'(a,2i8)') "nbnew,nrc=",nbnew,nrc
          nbnew=nbnew +nrc
        enddo
        tcom2= mpi_wtime()
        tcom= tcom +tcom2-tcom1
      enddo

      if(nbnew.ne.nb) then
        write(6,'(a,2i8)') "nbnew,(natm+nb)=",nbnew,natm+nb
        stop "error: nbnew.ne.(natm+nb)!!"
      endif
      
      end subroutine copy_strs_ba
c=======================================================================
      function v2_IWHe(r,is,js)
c
c  Two-body potential energy
c
      implicit none
      include "./params_au.h"
      include "params_AFS_W.h"
      real(8),intent(in):: r
      integer,intent(in):: is,js
      real(8):: v2_IWHe,x

      v2_IWHe= 0d0

c.....W-W
      if( r.lt.p_WW_c) then
        v2_IWHe= (r-p_WW_c)*(r-p_WW_c)*(p_WW_c0 +p_WW_c1*r
     &       +p_WW_c2*r*r)
      endif
c      if( r.lt.p_WW_b0 ) then
c        v2_IWHe=v2_IWHe
c     &       +p_WW_B*(p_WW_b0-r)**3 *exp(-p_WW_alpha*r)
c      endif

      return
      end function v2_IWHe
c=======================================================================
      function dv2_IWHe(r,is,js)
c
c  Derivative of two-body potential
c
      implicit none
      include "./params_au.h"
      include "params_AFS_W.h"
      real(8),intent(in):: r
      integer,intent(in):: is,js
      real(8):: dv2_IWHe,ri,x,exar

      dv2_IWHe= 0d0

c.....W-W
      if( r.lt.p_WW_c) then
        dv2_IWHe= 2d0*(r-p_WW_c)*(p_WW_c0 +p_WW_c1*r +p_WW_c2*r*r)
     &       +(r-p_WW_c)**2*(p_WW_c1 +2d0*p_WW_c2*r)
      endif
      if( r.lt.p_WW_b0 ) then
        dv2_IWHe=dv2_IWHe
     &       +p_WW_B*(p_WW_b0-r)**2 *exp(-p_WW_alpha*r)
     &       *(-3d0 -p_WW_alpha*(p_WW_b0-r))
      endif

      return
      end function dv2_IWHe
c=======================================================================
      function phi_IWHe(r,is,js)
c
c  Phi for many-body potential
c
      implicit none
      include "./params_au.h"
      include "params_AFS_W.h"
      real(8),intent(in):: r
      integer,intent(in):: is,js
      real(8):: phi_IWHe

      phi_IWHe= 0d0

c.....phi_W
      if( r.le.p_W_d ) then
        phi_IWHe= (r-p_W_d)**2
      endif

      return
      end function phi_IWHe
c=======================================================================
      function dphi_IWHe(r,is,js)
      implicit none
      include "./params_au.h"
      include "params_AFS_W.h"
      real(8),intent(in):: r
      integer,intent(in):: is,js
      real(8):: dphi_IWHe

      dphi_IWHe= 0d0

c.....phi_W
      if( r.le.p_W_d ) then
        dphi_IWHe= 2d0*(r-p_w_d)
      endif

      return
      end function dphi_IWHe
c=======================================================================
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make pmd"
c     End:
