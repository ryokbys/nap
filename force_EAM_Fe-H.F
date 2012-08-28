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
      am(2)       = am_h

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
c  Parallel implementation of EAM Ackland model for Fe (iron) and H.
C    - See Philos. Mag. 83(35) (2003) 3977--3994
c    - See also PRB 79, 174101 (2009) for Fe-H
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

      integer:: i,j,k,l,m,n,ierr,is,js,ixyz,jxyz
      real(8):: xij(3),rij,dfi,dfj,drhoij,drdxi(3),drdxj(3),at(3)
      real(8):: x,y,z,xi(3),epotl,tmp,t1,t2,t3,vemb
      real(8),external:: fphi,dfphi,fvphi,dfvphi,fpsi,dfpsi,femb,dfemb
     &     ,rho_hfe,rho_feh,rho_hh,fh,dfh,fvphi_feh,fvphi_hfe
     &     ,dfvphi_feh,fvphi_hh,dfvphi_hh
     &     ,drho_feh,drho_hfe,drho_hh

      logical,save:: l1st=.true.
      real(8),allocatable,save:: rho(:)
      real(8),save:: rs,rs_feh,avol

      if( l1st ) then
        allocate(rho(namax+nbmax))
        rs= a_rs /dsqrt(2d0)/z_fe**(1d0/3)
        rs_feh= a_rs /(z_fe**(2d0/3)+z_h**(2d0/3))
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
        is= int(tag(i))
        xi(1:3)= ra(1:3,i)
        do k=1,lspr(0,i)
          j=lspr(k,i)
          if(j.eq.0) exit
          js= int(tag(j))
          x= ra(1,j) -xi(1)
          y= ra(2,j) -xi(2)
          z= ra(3,j) -xi(3)
          xij(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
          rij=dsqrt(xij(1)*xij(1)+ xij(2)*xij(2) +xij(3)*xij(3))
          if( rij.gt.rc_rho ) cycle ! Fortunately, rc_rho is common
          if( is.eq.1 .and. js.eq.1 ) then ! Fe-Fe
            rho(i)=rho(i) +fpsi(rij)
          else if( is.eq.1 .and. js.eq.2 ) then ! Fe-H
            rho(i)=rho(i) +rho_hfe(rij)
          else if( is.eq.2 .and. js.eq.1 ) then ! H-Fe
            rho(i)=rho(i) +rho_feh(rij)
          else if( is.eq.2 .and. js.eq.2 ) then ! H-H
            rho(i)=rho(i) +rho_hh(rij)
          endif
        enddo
      enddo

c.....copy rho of boundary atoms
      call copy_rho_ba(tcom,namax,natm,nb,nbmax,lsb
     &     ,lsrc,myparity,nn,sv,mpi_md_world,rho)

c-----dE/dr_i
      do i=1,natm
        is= int(tag(i))
        xi(1:3)= ra(1:3,i)
        if( is.eq.1 ) then ! F_{Fe}
          vemb= femb(rho(i))
          dfi= dfemb(rho(i))
        else if( is.eq.2 ) then ! Fe_{H}
          vemb= fh(rho(i))
          dfi= dfh(rho(i))
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
          if( rij.gt.rc_vphi ) cycle
          drdxi(1:3)= -xij(1:3)/rij
c.....2-body term
          if( is.eq.1 .and. js.eq.1 ) then ! Fe-Fe
            t1= 0.5d0 *fvphi(rij,rs)
            t2= dfvphi(rij,rs)
          else if( is.ne.js ) then ! Fe-H
            t1= 0.5d0 *fvphi_feh(rij,rs_feh)
            t2= dfvphi_feh(rij,rs_feh)
          else if( is.eq.2 .and. js.eq.2 ) then ! H-H
            t1= 0.5d0 *fvphi_hh(rij)
            t2= dfvphi_hh(rij)
          endif
          epi(i)= epi(i) +t1
          epi(j)= epi(j) +t1
          if( j.le.natm ) then
            epotl=epotl +t1 +t1
          else
            epotl=epotl +t1
          endif
          aa(1:3,i)=aa(1:3,i) -drdxi(1:3)*t2
          aa(1:3,j)=aa(1:3,j) +drdxi(1:3)*t2
c.....atomic stress for 2-body terms
          do ixyz=1,3
            do jxyz=1,3
              strs(jxyz,ixyz,i)=strs(jxyz,ixyz,i)
     &             -0.5d0*t2*xij(ixyz)*(-drdxi(jxyz))
              strs(jxyz,ixyz,j)=strs(jxyz,ixyz,j)
     &             -0.5d0*t2*xij(ixyz)*(-drdxi(jxyz))
            enddo
          enddo
c.....Embedded term
          if( rij.gt.rc_rho ) cycle
          if( is.eq.1 .and. js.eq.1 ) then ! F_{Fe}, Fe-Fe
            dfj= dfemb(rho(j))
            tmp= (dfi+dfj) *dfpsi(rij)
          else if( is.eq.1 .and. js.eq.2 ) then ! F_{Fe}, Fe-H
            dfj= dfh(rho(j))
            tmp= dfi*drho_hfe(rij) +dfj*drho_feh(rij)
          else if( is.eq.2 .and. js.eq.1 ) then ! F_{H}, H-Fe
            dfj= dfemb(rho(j))
            tmp= dfi*drho_feh(rij) +dfj*drho_hfe(rij)
          else if( is.eq.2 .and. js.eq.2 ) then ! F_{H}, H-H
            dfj= dfh(rho(j))
            tmp= (dfi+dfj) *drho_hh(rij)
          endif
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
        epi(i)=epi(i) +vemb
        epotl=epotl +vemb
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
      subroutine get_force_FeFe(namax,natm,tag,ra,nnmax,aa,strs,h,hi,tcom
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

      end subroutine get_force_FeFe
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
c=======================================================================
      function rho_feh(r)
      implicit none
      include './params_au.h'
      include './params_EAM_Fe-H.h'
      real(8),intent(in):: r
      real(8):: rho_feh

      integer:: i
      real(8),external:: hvsd
      
      rho_feh= 0d0
      do i=1,6
        rho_feh= rho_feh +a_rho_feh(i)*(r_rho_feh(i)-r)**3
     &       *hvsd(r_rho_feh(i)-r)
      enddo
      return
      end function rho_feh
c=======================================================================
      function drho_feh(r)
      implicit none
      include './params_au.h'
      include './params_EAM_Fe-H.h'
      real(8),intent(in):: r
      real(8):: drho_feh
      integer:: i
      real(8),external:: hvsd

      drho_feh=0d0
      do i=1,6
        drho_feh=drho_feh -a_rho_feh(i)*(r_rho_feh(i)-r)**2
     &       *hvsd(r_rho_feh(i)-r)
      enddo
      drho_feh=drho_feh*3d0
      return
      end function drho_feh
c=======================================================================
      function rho_hfe(r)
      implicit none
      include './params_au.h'
      include './params_EAM_Fe-H.h'
      real(8),intent(in):: r
      real(8):: rho_hfe

      integer:: i
      real(8),external:: hvsd
      
      rho_hfe= 0d0
      do i=1,5
        rho_hfe= rho_hfe +a_rho_hfe(i)*(r_rho_hfe(i)-r)**3
     &       *hvsd(r_rho_hfe(i)-r)
      enddo
      return
      end function rho_hfe
c=======================================================================
      function drho_hfe(r)
      implicit none
      include './params_au.h'
      include './params_EAM_Fe-H.h'
      real(8),intent(in):: r
      real(8):: drho_hfe
      integer:: i
      real(8),external:: hvsd

      drho_hfe=0d0
      do i=1,5
        drho_hfe=drho_hfe -a_rho_hfe(i)*(r_rho_hfe(i)-r)**2
     &       *hvsd(r_rho_hfe(i)-r)
      enddo
      drho_hfe=drho_hfe*3d0
      return
      end function drho_hfe
c=======================================================================
      function rho_hh(r)
      implicit none
      include './params_au.h'
      include './params_EAM_Fe-H.h'
      real(8),intent(in):: r
      real(8):: rho_hh

      real(8),external:: fcut

      rho_hh= c_rho_hh *r**2 *exp(-2d0*r) *fcut(r)
      return
      end function rho_hh
c=======================================================================
      function drho_hh(r)
      implicit none
      include './params_au.h'
      include './params_EAM_Fe-H.h'
      real(8),intent(in):: r
      real(8):: drho_hh

      real(8),external:: fcut,dfcut
      real(8):: e1,fc
      
      e1= exp(-2d0*r)
      fc= fcut(r)
      drho_hh= c_rho_hh *(2d0*r*e1*fc -r**2*e1*fc +r**2*e1*dfcut(r))
      return
      end function drho_hh
c=======================================================================
      function fcut(r)
      implicit none
      include './params_au.h'
      include './params_EAM_Fe-H.h'
      real(8),intent(in):: r
      real(8):: fcut

      fcut= exp(aa2bohr/(r-rc_phi_hh))
      return
      end function fcut
c=======================================================================
      function dfcut(r)
      implicit none
      include './params_au.h'
      include './params_EAM_Fe-H.h'
      real(8),intent(in):: r
      real(8):: dfcut
      real(8),external:: fcut

      dfcut= -aa2bohr/(r-rc_phi_hh)**2 *fcut(r)
      return
      end function dfcut
c=======================================================================
      function fh(rho)
      implicit none
      include './params_au.h'
      include './params_EAM_Fe-H.h'
      real(8),intent(in):: rho
      real(8):: fh
      integer:: i
      
      fh= 0d0
      do i=1,6
        fh= fh +a_f(i) *rho**i
      enddo
      return
      end function fh
c=======================================================================
      function dfh(rho)
      implicit none
      include './params_au.h'
      include './params_EAM_Fe-H.h'
      real(8),intent(in):: rho
      real(8):: dfh
      integer:: i

      dfh= 0d0
      do i=1,6
        dfh=dfh +a_f(i) *rho**(i-1) *i
      enddo

      end function dfh
c=======================================================================
      function fvphi_hh(r)
      implicit none
      include './params_au.h'
      include './params_EAM_Fe-H.h'
      real(8),intent(in):: r
      real(8):: fvphi_hh

      real(8):: adag,emol,rho
      real(8),external:: rho_hh,s_hh,fh

      fvphi_hh= 0d0
      if( r.gt.rc_phi_hh ) return
      
      adag= (r-r0_hh)/(r0_hh*almbd_hh)
      emol= -2d0*eb_hh*(1d0+adag)*exp(-adag)
      rho= rho_hh(r)
      fvphi_hh= s_hh(r) *(emol -2d0*fh(rho))

      return
      end function fvphi_hh
c=======================================================================
      function dfvphi_hh(r)
      implicit none
      include './params_au.h'
      include './params_EAM_Fe-H.h'
      real(8),intent(in):: r
      real(8):: dfvphi_hh

      real(8):: adag,emol,rho
      real(8):: dsr,demolr,dfhr,drhor
      real(8),external:: rho_hh,drho_hh,s_hh,fh,dfh

      adag= (r-r0_hh)/(r0_hh*almbd_hh)
      emol= -2d0*eb_hh*(1d0+adag)*exp(-adag)
      rho= rho_hh(r)
      drhor= drho_hh(r)

c      dsr= -0.5d0 *a_tanh_hh *2d0 /cosh(a_tanh_hh*(r-r_tanh_hh))
      dsr= -a_tanh_hh /cosh(a_tanh_hh*(r-r_tanh_hh))
      demolr= 2d0 *eb_hh *adag *exp(-adag) /r0_hh /almbd_hh
      dfhr= drhor *dfh(rho)

      dfvphi_hh= dsr*(emol-2d0*fh(rho)) +s_hh(r)*(demolr-2d0*dfhr)
      return
      end function dfvphi_hh
c=======================================================================
      function s_hh(r)
      implicit none
      include './params_au.h'
      include './params_EAM_Fe-H.h'
      real(8),intent(in):: r
      real(8):: s_hh

      s_hh= 0.5d0 *(1d0 -tanh(a_tanh_hh*(r-r_tanh_hh)))
      return
      end function s_hh
c=======================================================================
      function fvphi_feh(r,rs)
      implicit none 
      include './params_au.h'
      include './params_EAM_Fe-H.h'
      real(8),intent(in):: r,rs
      real(8):: fvphi_feh
      integer:: i
      real(8),external:: fphi,hvsd
      real(8):: rr,r4

      fvphi_feh=0d0
      if( r.le.r1_feh ) then
        fvphi_feh= z2q2_feh/r*fphi(r/rs)
      else if( r.le.r2_feh ) then
        rr=r*r
        r4=rr*rr
        fvphi_feh= exp(b0_feh +b1_feh*r +b2_feh*rr
     &       +b3_feh*rr*r +b4_feh*r4 +b5_feh*r4*r) *ev2hrt
      else
        do i=1,7
          fvphi_feh=fvphi_feh +a_phi_feh(i)*(r_phi_feh(i)-r)**3
     &         *hvsd(r_phi_feh(i)-r)
        enddo
      endif
      return

      end function fvphi_feh
c=======================================================================
      function dfvphi_feh(r,rs)
      implicit none 
      include './params_au.h'
      include './params_EAM_Fe-H.h'
      real(8),intent(in):: r,rs
      real(8):: dfvphi_feh
      integer:: i
      real(8),external:: fphi,dfphi,hvsd
      real(8):: r4,r3,rr

      dfvphi_feh= 0d0
      if( r.le.r1_feh ) then
        dfvphi_feh= -z2q2_feh/r**2*fphi(r/rs)
     &       +z2q2_feh/r/rs*dfphi(r/rs)
      else if( r.le.r2_feh ) then
        rr=r**2
        r3=rr*r
        r4=r3*r
        dfvphi_feh= (b1_feh +2d0*b2_feh*r +3d0*b3_feh*rr
     &       +4d0*b4_feh*r3 +5d0*b5_feh*r4)
     &       *exp(b0_feh +b1_feh*r +b2_feh*rr +b3_feh*r3
     &       +b4_feh*r4 +b5_feh*r4*r) *ev2hrt
      else
        do i=1,7
          dfvphi_feh=dfvphi_feh -a_phi_feh(i)*(r_phi_feh(i)-r)**2
     &         *hvsd(r_phi_feh(i)-r)
        enddo
        dfvphi_feh=dfvphi_feh*3d0
      endif
      return
      end function dfvphi_feh
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make pmd"
c     End:
