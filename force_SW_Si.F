      subroutine setup(nismax,am,acon,fack,dt)
      implicit none
      include "params_au.h"
      include "params_SW_Si.h"
      integer,intent(in):: nismax
      real(8),intent(in):: dt
      real(8),intent(out):: am(nismax),acon(nismax),fack(nismax)

      integer:: i

c-----atomic masses
      am(1:nismax)= am_si

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
      subroutine get_force(namax,natm,tag,ra,nnmax,aa,h,hi,tcom
     &     ,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
     &     ,mpi_world,myid,epi,epot,nismax,acon)
c-----------------------------------------------------------------------
c  Parallel implementation of SW(Si) force calculation for pmd
c    - 2010.03.29 by R.K.
c      Made 1st version.
c    - 2011.04.15 by R.K.
c      Modified to 2 species with differet bond lengthes.
c      Species must be 1 and 2.
c-----------------------------------------------------------------------
      implicit none
      include "mpif.h"
      include "params_au.h"
      include "params_SW_Si.h"
      integer,intent(in):: namax,natm,nnmax,nismax
      integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3)
     &     ,nn(6),mpi_world,myid,lspr(0:nnmax,namax)
      real(8),intent(in):: ra(3,namax),tag(namax),acon(nismax)
     &     ,h(3,3),hi(3,3),sv(3,6),rc
      real(8),intent(inout):: tcom
      real(8),intent(out):: aa(3,namax),epi(namax),epot
c-----local
      integer:: i,j,k,l,m,n,ixyz,is,js,ks,ierr
      real(8):: rij,rik,riji,riki,rij2,rik2,rc2,swli,src2,src
      real(8):: tmp,vexp,df2,drij,csn,tcsn,tcsn2,dhrij,dhrik,dhcsn
      real(8):: drik,dcsn,drijc,drikc,x,y,z
      real(8):: epotl,epotl2,epotl3
      real(8),save,allocatable:: aa2(:,:),aa3(:,:)
      real(8),save,allocatable,dimension(:):: xi,xj,xk,xij,xik,at
c-----For 2-lattice-constant system
      real(8),save,allocatable:: tswrc(:,:),tswt(:,:)
c-----1st call
      logical,save:: l1st=.true.

c-----only at 1st call
      if( l1st ) then
        allocate(aa2(3,namax),aa3(3,namax))
        allocate(xi(3),xj(3),xk(3),xij(3),xik(3),at(3))
        allocate(tswrc(2,2),tswt(2,2))
        tswrc(1,1)= swrc
        tswrc(2,2)= swrc*ratio
        tswrc(1,2)= (tswrc(1,1)+tswrc(2,2))*0.5d0
        tswrc(2,1)= (tswrc(1,1)+tswrc(2,2))*0.5d0
        tswt(1,1)= swt
        tswt(2,2)= swt*ratio
        tswt(1,2)= (tswt(1,1)+tswt(2,2))*0.5d0
        tswt(2,1)= (tswt(1,1)+tswt(2,2))*0.5d0
c-------check rc
        if( int(rc).ne.int(max(tswrc(1,1),tswrc(2,2))*swl) ) then
          write(6,'(1x,a)') "!!! Cutoff radius is not appropriate !!!"
          write(6,'(1x,a,es12.4)') "rc should be"
     &         ,max(tswrc(1,1),tswrc(2,2))*swl*bohr
          call mpi_finalize(ierr)
          stop
        endif
c-------finally set l1st
        l1st=.false.
      endif
      
      epotl= 0d0
      epi(1:namax)= 0d0
      swli= 1d0/swl

c      write(6,'(a)') "2-body term"
c-----2 body term
      epotl2= 0d0
      aa2(1:3,1:namax)=0d0
      do i=1,natm
        xi(1:3)= ra(1:3,i)
        is= int(tag(i))
        do k=1,lspr(0,i)
          j=lspr(k,i)
          if(j.eq.0) exit
          if(j.le.i) cycle
          js= int(tag(j))
          xj= ra(1:3,j)-xi(1:3)
          xij(1:3)= ( h(1:3,1)*xj(1) +h(1:3,2)*xj(2)
     &         +h(1:3,3)*xj(3) )*swli
          rij2= xij(1)*xij(1) +xij(2)*xij(2) +xij(3)*xij(3)
          src= tswrc(is,js)
          src2= src*src
          if( rij2.ge.src2 ) cycle
          rij= dsqrt(rij2)
          riji= 1d0/rij
          drijc= 1d0/(rij-src)
          vexp=exp(swc*drijc)
c---------potential
          tmp= 0.5d0*swe *swa*vexp*(swb*riji**swp -riji**swq)
          epi(i)= epi(i) +tmp
          epotl2= epotl2 +tmp
          if( j.le.natm ) then
            epi(j)= epi(j) +tmp
            epotl2= epotl2 +tmp
          endif
c---------force
          df2= swe *swa*vexp*(-swp*swb*(riji**(swp+1d0))
     &         +swq*(riji**(swq+1d0))
     &         -(swb*(riji**swp) -riji**swq)*swc*drijc*drijc)
          do ixyz=1,3
            drij= -xij(ixyz)*riji
            aa2(ixyz,i)= aa2(ixyz,i) +df2*drij
            aa2(ixyz,j)= aa2(ixyz,j) -df2*drij
          enddo
        enddo
      enddo

c      write(6,'(a)') "3-body term"
c-----3 body term
      epotl3= 0d0
      aa3(1:3,1:namax)=0d0
c-----atom (i)
      do i=1,natm
        xi(1:3)=ra(1:3,i)
        is= int(tag(i))
        do n=1,lspr(0,i)
c---------atom (j)
          j=lspr(n,i)
          if(j.eq.0) exit
          if( j.eq.i ) cycle
          js= int(tag(j))
          xj(1:3)= ra(1:3,j) -xi(1:3)
          xij(1:3)= ( h(1:3,1)*xj(1) +h(1:3,2)*xj(2)
     &         +h(1:3,3)*xj(3) )*swli
          rij2= xij(1)*xij(1) +xij(2)*xij(2) +xij(3)*xij(3)
          src= tswrc(is,js)
          src2= src*src
          if( rij2.ge.src2 ) cycle
          rij= dsqrt(rij2)
          riji= 1d0/rij
          drijc= 1d0/(rij-src)
c---------atom (k)
          do m=1,lspr(0,i)
            k=lspr(m,i)
            if(k.eq.0) exit
            if( k.le.j .or. k.eq.i ) cycle
            ks= int(tag(k))
            xk(1:3)= ra(1:3,k) -xi(1:3)
            xik(1:3)= ( h(1:3,1)*xk(1) +h(1:3,2)*xk(2)
     &           +h(1:3,3)*xk(3) )*swli
            rik2= xik(1)*xik(1)+xik(2)*xik(2)+xik(3)*xik(3)
            src= tswrc(is,ks)
            src2= src*src
            if( rik2.ge.src2 ) cycle
            rik=dsqrt(rik2)
            riki= 1d0/rik
            drikc= 1d0/(rik-src)
c-----------common term
            csn=(xij(1)*xik(1) +xij(2)*xik(2) +xij(3)*xik(3))
     &           * (riji*riki)
            tcsn = csn +1d0/3d0
            tcsn2= tcsn*tcsn
            vexp= dexp(swt*drijc +swt*drikc)
c-----------potential
            tmp= swe *sws *vexp *tcsn2
            epi(i)= epi(i) +tmp
            epotl3= epotl3 +tmp
c-----------force
            dhrij= -sws *swt *vexp *tcsn2 *drijc*drijc
            dhrik= -sws *swt *vexp *tcsn2 *drikc*drikc
            dhcsn= 2d0 *sws *tcsn *vexp
            do l=1,3
              drij= -xij(l)*riji
              drik= -xik(l)*riki
              dcsn= xij(l)*csn*(riji*riji) +xik(l)*csn*(riki*riki)
     &             -(xij(l)+xik(l))*(riji*riki)
              aa3(l,i)=aa3(l,i) +swe*(dhcsn*dcsn +dhrij*drij
     &             +dhrik*drik)
c
              dcsn= -xij(l)*csn*(riji*riji) +xik(l)*(riji*riki)
              aa3(l,j)=aa3(l,j) +swe*(dhcsn*dcsn +dhrij*(-drij))
c
              dcsn= -xik(l)*csn*(riki*riki) +xij(l)*(riji*riki)
              aa3(l,k)=aa3(l,k) +swe*(dhcsn*dcsn +dhrik*(-drik))
            enddo
          enddo
        enddo
      enddo

c-----send back (3-body)forces and potentials on immigrants
      call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,lsrc,myparity
     &     ,nn,mpi_world,aa3,3)
      call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,lsrc,myparity
     &     ,nn,mpi_world,epi,1)

c-----sum
      aa(1:3,1:natm)= -aa2(1:3,1:natm) -aa3(1:3,1:natm)

c-----reduced force
      do i=1,natm
        at(1:3)= aa(1:3,i)
        aa(1:3,i)= ( hi(1:3,1)*at(1) +hi(1:3,2)*at(2)
     &       +hi(1:3,3)*at(3) ) *swli
      enddo
c-----multiply 0.5d0*dt**2/am(i)
      do i=1,natm
        is= int(tag(i))
        aa(1:3,i)= acon(is)*aa(1:3,i)
      enddo

c-----gather epot
      epot= 0d0
      epotl= epotl2 +epotl3
      call mpi_allreduce(epotl,epot,1,MPI_DOUBLE_PRECISION
     &     ,MPI_SUM,mpi_world,ierr)

c      deallocate(aa2,aa3)
c      deallocate(xi,xj,xk,xij,xik,at)
      end subroutine get_force
c=======================================================================
      subroutine get_force_stress(namax,natm,tag,ra,nnmax,aa,h,hi,tcom
     &     ,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
     &     ,mpi_world,myid,epi,epot,nismax,acon,strs)
c-----------------------------------------------------------------------
c  Parallel implementation of SW(Si) force calculation for pmd
c    - 2010.03.29 by R.K.
c      Made 1st version.
c    - 2011.04.15 by R.K.
c      Modified to 2 species with differet bond lengthes.
c      Species must be 1 and 2.
c-----------------------------------------------------------------------
      implicit none
      include "mpif.h"
      include "params_au.h"
      include "params_SW_Si.h"
      integer,intent(in):: namax,natm,nnmax,nismax
      integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3)
     &     ,nn(6),mpi_world,myid,lspr(0:nnmax,namax)
      real(8),intent(in):: ra(3,namax),tag(namax),acon(nismax)
     &     ,h(3,3),hi(3,3),sv(3,6),rc
      real(8),intent(inout):: tcom
      real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
c-----local
      integer:: i,j,k,l,m,n,ixyz,jxyz,is,js,ks,ierr,nbl
      real(8):: rij,rik,riji,riki,rij2,rik2,rc2,src2,src
      real(8):: tmp,tmp1,tmp2,vexp,df2,drij,csn,tcsn,tcsn2,dhrij,dhrik
     &     ,dhcsn,vol,voli,volj,volk
      real(8):: drik,dcsni,dcsnj,dcsnk,drijc,drikc,x,y,z,bl
      real(8):: epotl,epotl2,epotl3
      real(8),save:: swli,a8d3r3
      real(8),save,allocatable:: aa2(:,:),aa3(:,:)
      real(8),save,allocatable,dimension(:):: xi,xj,xk,xij,xik,at,bli
c-----For 2-lattice-constant system
      real(8),save,allocatable:: tswrc(:,:),tswt(:,:)
c-----1st call
      logical,save:: l1st=.true.

c-----only at 1st call
      if( l1st ) then
        allocate(aa2(3,namax),aa3(3,namax))
        allocate(xi(3),xj(3),xk(3),xij(3),xik(3),at(3),bli(namax))
        allocate(tswrc(2,2),tswt(2,2))
        tswrc(1,1)= swrc
        tswrc(2,2)= swrc*ratio
        tswrc(1,2)= (tswrc(1,1)+tswrc(2,2))*0.5d0
        tswrc(2,1)= (tswrc(1,1)+tswrc(2,2))*0.5d0
        tswt(1,1)= swt
        tswt(2,2)= swt*ratio
        tswt(1,2)= (tswt(1,1)+tswt(2,2))*0.5d0
        tswt(2,1)= (tswt(1,1)+tswt(2,2))*0.5d0
c-------check rc
        if( int(rc).ne.int(max(tswrc(1,1),tswrc(2,2))*swl) ) then
          write(6,'(1x,a)') "!!! Cutoff radius is not appropriate !!!"
          write(6,'(1x,a,es12.4)') "rc should be"
     &         ,max(tswrc(1,1),tswrc(2,2))*swl*bohr
          call mpi_finalize(ierr)
          stop
        endif
        swli= 1d0/swl
        a8d3r3= 8d0/(3d0*sqrt(3d0))
c-------finally set l1st
        l1st=.false.
      endif
      
      epotl= 0d0
      epi(1:natm+nb)= 0d0
      strs(1:3,1:3,1:namax)= 0d0

c      write(6,'(a)') "2-body term"
c-----2 body term
      epotl2= 0d0
      aa2(1:3,1:natm+nb)=0d0
      do i=1,natm
        xi(1:3)= ra(1:3,i)
        is= int(tag(i))
        bl= 0d0
        nbl= 0
        do k=1,lspr(0,i)
          j=lspr(k,i)
          if(j.eq.0) exit
          if(j.le.i) cycle
          js= int(tag(j))
          xj= ra(1:3,j)-xi(1:3)
          xij(1:3)= ( h(1:3,1)*xj(1) +h(1:3,2)*xj(2)
     &         +h(1:3,3)*xj(3) )*swli
          rij2= xij(1)*xij(1) +xij(2)*xij(2) +xij(3)*xij(3)
          src= tswrc(is,js)
          src2= src*src
          if( rij2.ge.src2 ) cycle
          rij= dsqrt(rij2)
          nbl=nbl +1
          bl=bl +rij
          riji= 1d0/rij
          drijc= 1d0/(rij-src)
          vexp=exp(swc*drijc)
c---------potential
          tmp= 0.5d0*swe *swa*vexp*(swb*riji**swp -riji**swq)
          epi(i)= epi(i) +tmp
          epotl2= epotl2 +tmp
          if( j.le.natm ) then
            epi(j)= epi(j) +tmp
            epotl2= epotl2 +tmp
          endif
c---------force
          df2= swe *swa*vexp*(-swp*swb*(riji**(swp+1d0))
     &         +swq*(riji**(swq+1d0))
     &         -(swb*(riji**swp) -riji**swq)*swc*drijc*drijc)
          do ixyz=1,3
            drij= -xij(ixyz)*riji
            aa2(ixyz,i)= aa2(ixyz,i) +df2*drij
            aa2(ixyz,j)= aa2(ixyz,j) -df2*drij
          enddo
c-----------Stress
          vol= a8d3r3*rij**3
          voli= 1d0/vol
          if( j.le.natm ) then
            do ixyz=1,3
              tmp= 0.5d0*(-df2*drij)
              do jxyz=1,3
                strs(ixyz,jxyz,i)= strs(ixyz,jxyz,i)
     &               -xij(jxyz)*tmp*voli
                strs(ixyz,jxyz,j)= strs(ixyz,jxyz,j)
     &               -xij(jxyz)*tmp*voli
              enddo
            enddo
          else
            do ixyz=1,3
              tmp= 0.5d0*(-df2*drij)
              do jxyz=1,3
                strs(ixyz,jxyz,i)= strs(ixyz,jxyz,i)
     &               -xij(jxyz)*tmp*voli
c                strs(ixyz,jxyz,j)= strs(ixyz,jxyz,j)
c     &               -xij(jxyz)*tmp
              enddo
            enddo
          endif
        enddo
c        if( nbl.ne.0 ) then
c          bli(i)= bl/nbl
c        else
c          bli(i)= 1d0
c        endif
      enddo

c      write(6,'(a)') "3-body term"
c-----3 body term
      epotl3= 0d0
      aa3(1:3,1:natm+nb)=0d0
c-----atom (i)
      do i=1,natm
        xi(1:3)=ra(1:3,i)
        is= int(tag(i))
        do n=1,lspr(0,i)
c---------atom (j)
          j=lspr(n,i)
          if(j.eq.0) exit
          if( j.eq.i ) cycle
          js= int(tag(j))
          xj(1:3)= ra(1:3,j) -xi(1:3)
          xij(1:3)= ( h(1:3,1)*xj(1) +h(1:3,2)*xj(2)
     &         +h(1:3,3)*xj(3) )*swli
          rij2= xij(1)*xij(1) +xij(2)*xij(2) +xij(3)*xij(3)
          src= tswrc(is,js)
          src2= src*src
          if( rij2.ge.src2 ) cycle
          rij= dsqrt(rij2)
          riji= 1d0/rij
          drijc= 1d0/(rij-src)
          vol= a8d3r3*rij**3
          volj= 1d0/vol
c---------atom (k)
          do m=1,lspr(0,i)
            k=lspr(m,i)
            if(k.eq.0) exit
            if( k.le.j .or. k.eq.i ) cycle
            ks= int(tag(k))
            xk(1:3)= ra(1:3,k) -xi(1:3)
            xik(1:3)= ( h(1:3,1)*xk(1) +h(1:3,2)*xk(2)
     &           +h(1:3,3)*xk(3) )*swli
            rik2= xik(1)*xik(1)+xik(2)*xik(2)+xik(3)*xik(3)
            src= tswrc(is,ks)
            src2= src*src
            if( rik2.ge.src2 ) cycle
            rik=dsqrt(rik2)
            riki= 1d0/rik
            drikc= 1d0/(rik-src)
            vol= a8d3r3*rik**3
            volk= 1d0/vol
c-----------common term
            csn=(xij(1)*xik(1) +xij(2)*xik(2) +xij(3)*xik(3))
     &           * (riji*riki)
            tcsn = csn +1d0/3d0
            tcsn2= tcsn*tcsn
            vexp= dexp(swt*drijc +swt*drikc)
c-----------potential
            tmp= swe *sws *vexp *tcsn2
            epi(i)= epi(i) +tmp
            epotl3= epotl3 +tmp
c-----------force
            dhrij= -sws *swt *vexp *tcsn2 *drijc*drijc
            dhrik= -sws *swt *vexp *tcsn2 *drikc*drikc
            dhcsn= 2d0 *sws *tcsn *vexp
            do l=1,3
              drij= -xij(l)*riji
              drik= -xik(l)*riki
              dcsnj= -xij(l)*csn*(riji*riji) +xik(l)*(riji*riki)
              dcsnk= -xik(l)*csn*(riki*riki) +xij(l)*(riji*riki)
              dcsni= -dcsnj -dcsnk
c              dcsn= xij(l)*csn*(riji*riji) +xik(l)*csn*(riki*riki)
c     &             -(xij(l)+xik(l))*(riji*riki)
              aa3(l,i)=aa3(l,i) +swe*(dhcsn*dcsni +dhrij*drij
     &             +dhrik*drik)
c
              tmp1= swe*(dhcsn*dcsnj +dhrij*(-drij))
              aa3(l,j)=aa3(l,j) +tmp1
c
              tmp2= swe*(dhcsn*dcsnk +dhrik*(-drik))
              aa3(l,k)=aa3(l,k) +tmp2
c-------------Stress
              do jxyz=1,3
                strs(l,jxyz,i)=strs(l,jxyz,i)
     &               -xij(jxyz)*(-tmp1)*volj -xik(jxyz)*(-tmp2)*volk
                strs(l,jxyz,j)=strs(l,jxyz,j)
     &               -xij(jxyz)*(-tmp1)*volj
                strs(l,jxyz,k)=strs(l,jxyz,k)
     &               -xik(jxyz)*(-tmp2)*volk
              enddo
            enddo
          enddo
        enddo
      enddo

c-----send back (3-body)forces, stresses, and potentials on immigrants
      call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,lsrc,myparity
     &     ,nn,mpi_world,aa3,3)
      call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,lsrc,myparity
     &     ,nn,mpi_world,epi,1)
      call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,lsrc,myparity
     &     ,nn,mpi_world,strs,9)

c-----sum
      aa(1:3,1:natm)= -aa2(1:3,1:natm) -aa3(1:3,1:natm)

c-----reduced force
      do i=1,natm
        at(1:3)= aa(1:3,i)
        aa(1:3,i)= ( hi(1:3,1)*at(1) +hi(1:3,2)*at(2)
     &       +hi(1:3,3)*at(3) ) *swli
      enddo
c-----multiply 0.5d0*dt**2/am(i)
      do i=1,natm
        is= int(tag(i))
        aa(1:3,i)= acon(is)*aa(1:3,i)
      enddo

c-----gather epot
      epot= 0d0
      epotl= epotl2 +epotl3
      call mpi_allreduce(epotl,epot,1,MPI_DOUBLE_PRECISION
     &     ,MPI_SUM,mpi_world,ierr)

c      deallocate(aa2,aa3)
c      deallocate(xi,xj,xk,xij,xik,at)
      end subroutine get_force_stress
c=======================================================================
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make pmd"
c     End:
