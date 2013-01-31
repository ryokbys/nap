      subroutine setup(nismax,am,acon,fack,dt)
      implicit none
      include "params_au.h"
      include "params_Ito_W-He.h"
      integer,intent(in):: nismax
      real(8),intent(in):: dt
      real(8),intent(out):: am(nismax),acon(nismax),fack(nismax)

      integer:: i

c-----atomic masses
      am(1:nismax)= am_W
      am(2)       = am_He

c-----prefactors for normalized accelerration
      do i=1,nismax
        acon(i)= 0.5d0*dt**2/am(i)
      enddo

c-----prefactors for kinetic energy, FACK
      do i=1,nismax
        fack(i)= 0.5d0*am(i)/dt**2
      enddo

      end subroutine setup
c=======================================================================
      subroutine get_force(namax,natm,tag,ra,nnmax,aa,strs,h,hi,tcom
     &     ,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
     &     ,mpi_md_world,myid_md,epi,epot,nismax,acon)
c-----------------------------------------------------------------------
c  Parallel implementation of Ito potential for W and He (IWHe)
c    - smoothing is applied to 2-body potential for W-He and He-He
c    - rho of boundary atoms are sent to the neighbor nodes
c    - only force on i is calculated, not necessary to send-back
c-----------------------------------------------------------------------
c  See RK's molecular dynamics memo
c-----------------------------------------------------------------------
      implicit none
      include "mpif.h"
      include "params_au.h"
      include "params_Ito_W-He.h"
      integer,intent(in):: namax,natm,nnmax,nismax
      integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3)
     &     ,nn(6),mpi_md_world,myid_md
      integer,intent(inout):: lspr(0:nnmax,namax)
      real(8),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),sv(3,6)
     &     ,acon(nismax),rc,tag(namax)
      real(8),intent(inout):: tcom
      real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)

      integer:: i,j,k,l,m,n,ierr,is,js
      real(8):: xij(3),rij,dfi,dfj,drhoij,drdxi(3),drdxj(3),r,dphi,at(3)
      real(8):: x,y,z,xi(3),epotl,v2,dv2
      real(8),save:: v2c_WW,v2c_WHe,v2c_HeHe,dv2c_WW,dv2c_WHe,dv2c_HeHe
      logical,save:: l1st=.true.
      real(8),allocatable,save:: sqrho(:)
      real(8),external:: v2_IWHe,dv2_IWHe,phi_IWHe,dphi_IWHe

      if( l1st ) then
        allocate(sqrho(namax+nbmax))
c.....Smoothing for 2-body terms
        v2c_WW  = v2_IWHe(p_WW_c,1,1)
        v2c_WHe = v2_IWHe(p_WHe_rc,1,2)
        v2c_HeHe= v2_IWHe(p_HeHe_rc,2,2)
        dv2c_WW  = dv2_IWHe(p_WW_c,1,1)
        dv2c_WHe = dv2_IWHe(p_WHe_rc,1,2)
        dv2c_HeHe= dv2_IWHe(p_HeHe_rc,2,2)
        l1st=.false.
      endif

      aa(1:3,1:natm)=0d0
      epi(1:natm)= 0d0
      epotl= 0d0
      sqrho(1:natm)= 0d0
      strs(1:3,1:3,natm+nb)= 0d0

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
        if( is.eq.1 ) sqrho(i)= dsqrt(sqrho(i))
      enddo

c-----copy rho of boundary atoms
      call copy_rho_ba(tcom,namax,natm,nb,nbmax,lsb,lsrc,myparity,nn,sv
     &     ,mpi_md_world,sqrho)

      do i=1,natm
        xi(1:3)= ra(1:3,i)
        is= int(tag(i))
        if( is.eq.1 ) then
          dfi= -0.5d0/sqrho(i)
        elseif( is.eq.2 ) then
          dfi= p_He_A
        endif
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
          if( is.eq.2 .and. js.eq.2 ) then
            v2= v2 -v2c_HeHe -(r-p_HeHe_rc)*dv2c_HeHe
            dv2=dv2 -dv2c_HeHe
          elseif( is.ne.js ) then
            v2= v2 -v2c_WHe -(r-p_WHe_rc)*dv2c_WHe
            dv2=dv2 -dv2c_WHe
          endif
          epi(i)= epi(i) +v2
          epi(j)= epi(j) +v2
          if(j.le.natm) then
            epotl=epotl +v2 +v2
          else
            epotl=epotl +v2
          endif
          aa(1:3,i)=aa(1:3,i) -dv2*drdxi(1:3)
          aa(1:3,j)=aa(1:3,j) +dv2*drdxi(1:3)
c.....N-body term
          dphi= dphi_IWHe(rij,is,js)
          if( js.eq.1 ) then
            dfj= -0.5d0/sqrho(j)
          elseif( js.eq.2 ) then
            dfj= p_He_A
          endif
          aa(1:3,i)=aa(1:3,i) -(dfi+dfj)*dphi*drdxi(1:3)
          aa(1:3,j)=aa(1:3,j) +(dfi+dfj)*dphi*drdxi(1:3)
        enddo
c        write(6,'(a,3es12.4)') 'sqrho(i)=',sqrho(i)
        if( is.eq.1 ) then
          epi(i)=epi(i) -sqrho(i)
          epotl=epotl -sqrho(i)
        elseif( is.eq.2 ) then
          epi(i)=epi(i) +sqrho(i)*p_He_A
          epotl=epotl +sqrho(i)*p_He_A
        endif
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
      function v2_IWHe(r,is,js)
c
c  Two-body potential energy
c
      implicit none
      include "params_au.h"
      include "params_Ito_W-He.h"
      real(8),intent(in):: r
      integer,intent(in):: is,js
      real(8):: v2_IWHe,x
      real(8),external:: fc

      v2_IWHe= 0d0

c.....W-W
      if( is.eq.1 .and. js.eq.1 ) then
        if( r.lt.p_WW_c) then
          v2_IWHe= (r-p_WW_c)*(r-p_WW_c)*(p_WW_c0 +p_WW_c1*r
     &         +p_WW_c2*r*r)
        endif
        if( r.lt.p_WW_b0 ) then
          v2_IWHe=v2_IWHe
     &         +p_WW_B*(p_WW_b0-r)**3 *exp(-p_WW_alpha*r)
        endif

c.....He-He
      elseif( is.eq.2 .and. js.eq.2 ) then
        if( r.lt.p_HeHe_rc ) then
          x= (r-p_HeHe_rs)/p_HeHe_cs
          v2_IWHe= p_HeHe_z1/r *(exp(-p_HeHe_a1*r) +p_HeHe_b0) *fc(x)
        endif

c.....W-He or He-W
      elseif( is.ne.js ) then
        if( r.lt.p_WHe_rc ) then
          x= (r-p_WHe_rs)/p_WHe_cs
          v2_IWHe= p_WHe_z1/r *(exp(-p_WHe_a1*r) +p_WHe_b0) *fc(x)
        endif
      endif

      return
      end function v2_IWHe
c=======================================================================
      function dv2_IWHe(r,is,js)
c
c  Derivative of two-body potential
c
      implicit none
      include "params_au.h"
      include "params_Ito_W-He.h"
      real(8),intent(in):: r
      integer,intent(in):: is,js
      real(8):: dv2_IWHe,ri,x,exar
      real(8),external:: fc,dfc

      dv2_IWHe= 0d0

c.....W-W
      if( is.eq.1 .and. js.eq.1 ) then
        if( r.lt.p_WW_c) then
          dv2_IWHe= 2d0*(r-p_WW_c)*(p_WW_c0 +p_WW_c1*r +p_WW_c2*r*r)
     &         +(r-p_WW_c)**2*(p_WW_c1 +2d0*p_WW_c2*r)
        endif
        if( r.lt.p_WW_b0 ) then
          dv2_IWHe=dv2_IWHe
     &         +p_WW_B*(p_WW_b0-r)**2 *exp(-p_WW_alpha*r)
     &         *(-3d0 -p_WW_alpha*(p_WW_b0-r))
        endif

c.....He-He
      elseif( is.eq.2 .and. js.eq.2 ) then
        if( r.lt.p_HeHe_rc ) then
          x= (r-p_HeHe_rs)/p_HeHe_cs
          ri= 1d0/r
          exar= exp(-p_HeHe_a1*r)
          dv2_IWHe= -p_HeHe_z1*ri*fc(x)
     &         *(ri*(exar+p_HeHe_b0) +p_HeHe_a1*exar)
     &         +p_HeHe_z1*ri*(exar+p_HeHe_b0)*dfc(x)/p_HeHe_cs
        endif

c.....W-He or He-W
      elseif( is.ne.js ) then
        if( r.lt.p_WHe_rc ) then
          x= (r-p_WHe_rs)/p_WHe_cs
          ri= 1d0/r
          exar= exp(-p_WHe_a1*r)
          dv2_IWHe= -p_WHe_z1*ri*fc(x)
     &         *(ri*(exar+p_WHe_b0) +p_WHe_a1*exar)
     &         +p_WHe_z1*ri*(exar+p_WHe_b0)*dfc(x)/p_WHe_cs
        endif
      endif

      return
      end function dv2_IWHe
c=======================================================================
      function phi_IWHe(r,is,js)
c
c  Phi for many-body potential
c
      implicit none
      include "params_au.h"
      include "params_Ito_W-He.h"
      real(8),intent(in):: r
      integer,intent(in):: is,js
      real(8):: phi_IWHe

      phi_IWHe= 0d0

c.....No need to consider, He-He pair
      if( is.eq.2 .and. js.eq.2 ) return

c.....phi_W
      if( js.eq.1 ) then
        if( r.le.p_W_d ) then
          phi_IWHe= p_W_w**2 *(r-p_W_d)**2
     &         *(1d0 +p_W_beta/p_W_d*(r-p_W_d))
        endif
        
c.....phi_He
      elseif( js.eq.2 ) then
        if( r.le.p_He_d ) then
          phi_IWHe= p_He_w**2 *(r-p_He_d)**2
        endif

      endif

      return
      end function phi_IWHe
c=======================================================================
      function dphi_IWHe(r,is,js)
      implicit none
      include "params_au.h"
      include "params_Ito_W-He.h"
      real(8),intent(in):: r
      integer,intent(in):: is,js
      real(8):: dphi_IWHe,rd

      dphi_IWHe= 0d0

c.....No need to consider, He-He pair
      if( is.eq.2 .and. js.eq.2 ) return

c.....phi_W
      if( js.eq.1 ) then
        if( r.le.p_W_d ) then
          rd= (r-p_W_d)
          dphi_IWHe= 2d0*p_W_w**2 *rd 
     &         +3d0 *p_W_w**2 *rd**2 *p_W_beta/p_W_d
        endif
        
c.....phi_He
      elseif( js.eq.2 ) then
        if( r.le.p_He_d ) then
          dphi_IWHe= 2d0*p_He_w**2 *(r-p_He_d)
        endif
      endif

      return
      end function dphi_IWHe
c=======================================================================
      function fc(x)
      implicit none
      real(8),intent(in):: x
      real(8):: fc

      if( x.lt.0d0 ) then
        fc= 1d0
      elseif( x.ge.0d0 .and. x.lt.1d0 ) then
        fc= (-6d0*x*x +15d0*x -10d0)*x*x*x +1d0
      else
        fc= 0d0
      endif
      return
      end function fc
c=======================================================================
      function dfc(x)
      implicit none
      real(8),intent(in):: x
      real(8):: dfc

      if( x.lt.0d0 ) then
        dfc= 0d0
      elseif( x.ge.0d0 .and. x.lt.1d0 ) then
        dfc= -30d0*x*x*(x-1d0)*(x-1d0)
      else
        dfc= 0d0
      endif
      return
      end function dfc
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make pmd"
c     End:
