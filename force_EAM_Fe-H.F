      subroutine setup(nismax,am,acon,fack,dt)
      implicit none
      include "params_au.h"
      include "params_EAM_Fe-H.h"
      integer,intent(in):: nismax
      real(8),intent(in):: dt
      real(8),intent(out):: am(nismax),acon(nismax),fack(nismax)

      integer:: i

c-----atomic masses
      am(1:nismax)= am_fe

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
c  Parallel implementation of EAM Ackland model for Fe (iron).
C    - See Philos. Mag. 83(35) (2003) 3977--3994
c    - rho of boundary atoms are sent to the neighbor nodes
c    - only force on i is calculated, not necessary to send-back
c-----------------------------------------------------------------------
      implicit none
      include "mpif.h"
      include "params_au.h"
      include "params_EAM_Fe-H.h"
      integer,intent(in):: namax,natm,nnmax,nismax
      integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3)
     &     ,nn(6),mpi_md_world,myid_md
      integer,intent(in):: lspr(0:nnmax,namax)
      real(8),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),sv(3,6)
     &     ,acon(nismax),rc,tag(namax)
      real(8),intent(inout):: tcom
      real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)

      integer:: i,j,k,l,m,n,ierr,is,ixyz,jxyz
      real(8):: xij(3),rij,dfi,dfj,drhoij,drdxi(3),drdxj(3),at(3)
      real(8):: x,y,z,xi(3),epotl,tmp
      real(8),external:: fphi,dfphi,fvphi,dfvphi,fpsi,dfpsi,femb,dfemb

      logical,save:: l1st=.true.
      real(8),allocatable,save:: rho(:)
      real(8),save:: rs,avol

      if( l1st ) then
        allocate(rho(namax+nbmax))
        rs=  a_rs /dsqrt(2d0)/z_fe**(1d0/3)
c.....assuming fixed (constant) atomic volume (BCC)
        avol= alcfe**3 /2
        call set_avol()
        l1st=.false.
      endif

      aa(1:3,1:natm)=0d0
      epi(1:natm)= 0d0
      epotl= 0d0
      rho(1:natm)= 0d0
      strs(1:3,1:3,1:natm+nb)= 0d0

c.....rho(i)
      do i=1,natm
        xi(1:3)= ra(1:3,i)
        do k=1,lspr(0,i)
          j=lspr(k,i)
          if(j.eq.0) exit
          x= ra(1,j) -xi(1)
          y= ra(2,j) -xi(2)
          z= ra(3,j) -xi(3)
          xij(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
          rij=dsqrt(xij(1)*xij(1)+ xij(2)*xij(2) +xij(3)*xij(3))
          if( rij.gt.rc_rho ) cycle
          rho(i)=rho(i) +fpsi(rij)
        enddo
      enddo

c.....copy rho of boundary atoms
      call copy_rho_ba(tcom,namax,natm,nb,nbmax,lsb
     &     ,lsrc,myparity,nn,sv,mpi_md_world,rho)

c-----dE/dr_i
      do i=1,natm
        xi(1:3)= ra(1:3,i)
        dfi= dfemb(rho(i))
        do k=1,lspr(0,i)
          j=lspr(k,i)
          if(j.eq.0) exit
          if(j.le.i) cycle
          x= ra(1,j) -xi(1)
          y= ra(2,j) -xi(2)
          z= ra(3,j) -xi(3)
          xij(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
          rij=sqrt(xij(1)**2+ xij(2)**2 +xij(3)**2)
          if( rij.gt.rc_vphi ) cycle
          drdxi(1:3)= -xij(1:3)/rij
c.....2-body term
          tmp= 0.5d0*fvphi(rij,rs)
          epi(i)= epi(i) +tmp
          epi(j)= epi(j) +tmp
          if( j.le.natm ) then
            epotl=epotl +tmp +tmp
          else
            epotl=epotl +tmp
          endif
          tmp= dfvphi(rij,rs)
          aa(1:3,i)=aa(1:3,i) -drdxi(1:3)*tmp
          aa(1:3,j)=aa(1:3,j) +drdxi(1:3)*tmp
c.....atomic stress for 2-body terms
          do ixyz=1,3
            do jxyz=1,3
              strs(jxyz,ixyz,i)=strs(jxyz,ixyz,i)
     &             -0.5d0*tmp*xij(ixyz)*(-drdxi(jxyz))
              strs(jxyz,ixyz,j)=strs(jxyz,ixyz,j)
     &             -0.5d0*tmp*xij(ixyz)*(-drdxi(jxyz))
            enddo
          enddo
c.....Embedded term
          if( rij.gt.rc_rho ) cycle
          dfj= dfemb(rho(j))
          tmp= (dfi+dfj)*dfpsi(rij)
          aa(1:3,i)=aa(1:3,i) -drdxi(1:3)*tmp
          aa(1:3,j)=aa(1:3,j) +drdxi(1:3)*tmp
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
        tmp= femb(rho(i))
        epi(i)=epi(i) +tmp
        epotl=epotl +tmp
      enddo

c-----copy strs of boundary atoms
      call copy_strs_ba(tcom,namax,natm,nb,nbmax,lsb
     &     ,lsrc,myparity,nn,sv,mpi_md_world,strs)
c-----atomic level stress in [Hartree/Bohr^3] assuming 1 Bohr thick
      do i=1,natm
        strs(1:3,1:3,i)= strs(1:3,1:3,i) /avol
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

      end subroutine get_force
c=======================================================================
      function hvsd(x)
c
c  Heaviside's stepwise function
c
      implicit none
      real(8),intent(in):: x
      real(8):: hvsd

      hvsd= 0d0
      if( x.ge.0 ) then
        hvsd= 1d0
        return
      endif
      return 

      end function hvsd
c=======================================================================
      function fphi(x)
c
c  Screening function for r < r1
c
      implicit none
      real(8),intent(in):: x
      real(8):: fphi

      fphi= 0.1818d0*exp(-3.2d0*x)
     &     +0.5099d0*exp(-0.9423d0*x)
     &     +0.2802d0*exp(-0.4029d0*x)
     &     +0.02817d0*exp(-0.2016d0*x)
      return
      end function fphi
c=======================================================================
      function dfphi(x)
c
c  1st derivative of the screening function for r < r1
c
      implicit none 
      real(8),intent(in):: x
      real(8):: dfphi

      dfphi= -0.58176d0*exp(-3.2d0*x)
     &     -0.48047877d0*exp(-0.9423d0*x)
     &     -0.11289258d0*exp(-0.4029d0*x)
     &     -0.005679072d0*exp(-0.2016d0*x)
      return
      end function dfphi
c=======================================================================
      function fvphi(r,rs)
c
c  varphi function for r > r2
c
      implicit none 
      include './params_au.h'
      include './params_EAM_Fe-H.h'
      real(8),intent(in):: r,rs
      real(8):: fvphi
      integer:: i
      real(8),external:: fphi,hvsd

      fvphi=0d0
      if( r.le.r1 ) then
        fvphi=fvphi +z2q2/r*fphi(r/rs)
      else if( r.le.r2 ) then
        fvphi=fvphi +exp(b0 +b1*r +b2*r**2 +b3*r**3) *ev2hrt
      else
        do i=2,14
          fvphi=fvphi +a_vphi(i)*(r_vphi(i)-r)**3
     &         *hvsd(r_vphi(i)-r)
        enddo
      endif
      return
      end function fvphi
c=======================================================================
      function dfvphi(r,rs)
c
c  1st derivative of varphi func for r > r2
c
      implicit none 
      include './params_au.h'
      include './params_EAM_Fe-H.h'
      real(8),intent(in):: r,rs
      real(8):: dfvphi
      integer:: i
      real(8),external:: fphi,dfphi,hvsd

      dfvphi= 0d0
      if( r.le.r1 ) then
        dfvphi=dfvphi -z2q2/r**2*fphi(r/rs) +z2q2/r/rs*dfphi(r/rs)
      else if( r.le.r2 ) then
        dfvphi=dfvphi +(b1 +2d0*b2*r +3d0*b3*r**2)
     &       *exp(b0 +b1*r +b2*r**2 +b3*r**3) *ev2hrt
      else
        do i=2,14
          dfvphi=dfvphi -a_vphi(i)*(r_vphi(i)-r)**2
     &         *hvsd(r_vphi(i)-r)
        enddo
        dfvphi=dfvphi*3d0
      endif
      return
      end function dfvphi
c=======================================================================
      function fpsi(r)
c
c  Cubic spline function for calculating rho
c
      implicit none 
      include './params_au.h'
      include './params_EAM_Fe-H.h'
      real(8),intent(in):: r
      real(8):: fpsi
      integer:: i
      real(8),external:: hvsd

      fpsi=0d0
      do i=1,3
        fpsi=fpsi +a_psi(i)*(r_psi(i)-r)**3
     &       *hvsd(r_psi(i)-r)
      enddo
      return
      end function fpsi
c=======================================================================
      function dfpsi(r)
c
c  1st derivative of the cubic spline func of calculation of rho
c
      implicit none 
      include './params_au.h'
      include './params_EAM_Fe-H.h'
      real(8),intent(in):: r
      real(8):: dfpsi
      integer:: i
      real(8),external:: hvsd

      dfpsi=0d0
      do i=1,3
        dfpsi=dfpsi -a_psi(i)*(r_psi(i)-r)**2
     &       *hvsd(r_psi(i)-r)
      enddo
      dfpsi=dfpsi*3d0
      return
      end function dfpsi
c=======================================================================
      function femb(rho)
c
c  Embedding function
c
      implicit none
      include './params_au.h'
      include './params_EAM_Fe-H.h'
      real(8),intent(in):: rho
      real(8):: femb
      
      femb= (-dsqrt(rho) +a_emb*rho*rho)*ev2hrt
      return
      end function femb
c=======================================================================
      function dfemb(rho)
c
c  1st derivative of the embedding function
c
      implicit none 
      include './params_au.h'
      include './params_EAM_Fe-H.h'
      real(8),intent(in):: rho
      real(8):: dfemb

      dfemb= (-0.5d0/dsqrt(rho)+2d0*a_emb*rho)*ev2hrt
      return
      end function dfemb
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
      subroutine set_avol()
      use variables
      implicit none
      include './params_au.h'
      include './params_EAM_Fe-H.h'

      avol= alcfe**3 /2
      if(myid_md.eq.0) write(6,'(a,es12.4)') ' avol =',avol
      return
      end subroutine set_avol
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make pmd"
c     End:
