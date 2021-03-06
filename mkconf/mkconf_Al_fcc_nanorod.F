      program mkconf_Al_nanorod
c-----------------------------------------------------------------------
c  Make atom configuration file of fcc Al nanorod
c-----------------------------------------------------------------------
c  OUTPUT
c  ------
c    - pmd0000
c-----------------------------------------------------------------------
      implicit none
      include "../pmd/params_unit.h"
      include "../pmd/params_EAM_Al.h"
c-----max # of atoms
      integer,parameter::nmax=1000000
c-----# of unit cells
      integer,parameter:: nuc(1:3)= (/ 15, 15, 20 /)
c-----vacuum width in unit of cell
      integer,parameter:: nvac(1:3)= (/ 0, 0, 20 /)
c.....Radius of base of the cylinder
      real(8),parameter:: rad  = 10.d0

      integer:: i,ix,iy,iz,m,inc,is,ifmv
      real(8):: small,al,rad2,rx,ry,r2
      real(8):: ua(3,10)
      real(8):: ra(3,nmax)
      real(8):: h(3,3,0:1)
      real(8):: va(3,nmax)
      real(8):: tag(nmax),epi(nmax),eki(nmax),strs(3,3,nmax)
c-----initial temperature
      real(8),parameter:: tempini= 100d0

c.....In this purpose nuc(1)==nuc(2)
      if( nuc(1).ne.nuc(2) )
     &     stop '[Error] nuc(1).ne.nuc(2) !!!'

      small=1d-7

c-----box size for Al: 4.025 Ang
      al= 4.025d0

c-----simulation box size
      h(1:3,1:3,0:1)= 0d0
      h(1,1,0)= (nuc(1)+nvac(1)) *al
      h(2,2,0)= (nuc(2)+nvac(2)) *al
      h(3,3,0)= (nuc(3)+nvac(3)) *al
c-----unit cell, FCC
c-----Al
      ua(1:3,1)= (/ 0.00d0, 0.00d0, 0.00d0 /)
      ua(1:3,2)= (/ 0.50d0, 0.50d0, 0.00d0 /)
      ua(1:3,3)= (/ 0.50d0, 0.00d0, 0.50d0 /)
      ua(1:3,4)= (/ 0.00d0, 0.50d0, 0.50d0 /)

      rad2= rad*rad
      
      inc=0 
      do ix=0,nuc(1)-1
        do iy=0,nuc(2)-1
          do iz=0,nuc(3)-1
            do m=1,4
              inc=inc+1
              if(inc.gt.nmax)then
                write(*,*)'Error inc>nmax',inc,nmax
                stop
              endif
              ra(1,inc)= (ua(1,m)+dble(ix))/(nuc(1)+nvac(1)) +small
              ra(2,inc)= (ua(2,m)+dble(iy))/(nuc(2)+nvac(2)) +small
              ra(3,inc)= (ua(3,m)+dble(iz)+dble(nvac(3))/2)
     &             /(nuc(3)+nvac(3)) +small
c.....Only inside atoms the radius remain to form cylinder
              rx= (ra(1,inc)-0.5d0)*h(1,1,0)
              ry= (ra(2,inc)-0.5d0)*h(2,2,0)
              r2= rx**2 +ry**2
              if( r2.gt.rad2 ) then
                inc=inc-1
                cycle
              endif
              is= 1
              ifmv= 1
              if( iz.eq.nuc(3)-1 ) then ! top layer
                ifmv= 2
              else if( iz.eq.0 ) then ! bottom layer
                ifmv= 2
              endif
              tag(inc)= 1d0  *is
     &                 +1d-1 *ifmv
     &                 +1d-14*inc
            enddo
          enddo
        enddo
      enddo
      write(6,'(a,i10)') " natm=",inc

      call setv(inc,va,tag,tempini)

c-----scale velocities to reduced unit
      do i=1,inc
        va(1,i)=va(1,i) /h(1,1,0)
        va(2,i)=va(2,i) /h(2,2,0)
        va(3,i)=va(3,i) /h(3,3,0)
      enddo
      
      call write_pmd0_ascii(15,'pmd0000','replace',inc,tag,ra,va
     &     ,h,al,eki,epi,strs)
      
      end program mkconf_Al_nanorod
c=======================================================================
      subroutine setv(natm,va,tag,tempini)
      implicit none
      include "../pmd/params_unit.h"
      include "../pmd/params_EAM_Al.h"
      integer,intent(in):: natm
      real(8),intent(in):: tempini,tag(natm)
      real(8),intent(out):: va(3,natm)

      integer:: i,l,is
      real(8):: dseed,sumvx,sumvy,sumvz,rnd1,rnd2,tmp
      real(8),parameter:: pi = 3.14159265358979d0
      real(8):: facv,am

      facv=dsqrt(2d0*tempini*fkb/am_al)
      am= am_al

c-----velocities in Maxwell-Boltzmann distribution
      dseed=12345
      do i=1,natm
        is= int(tag(i))
        do l=1,3
          call myrnd(rnd1,dseed)
          call myrnd(rnd2,dseed)
          va(l,i)=facv*dsqrt(-dlog(rnd1))*dcos(2d0*pi*rnd2)
        enddo
      enddo
c-----set center of mass motion to zero
      sumvx=0d0
      sumvy=0d0
      sumvz=0d0
      do i=1,natm
        sumvx=sumvx+va(1,i)
        sumvy=sumvy+va(2,i)
        sumvz=sumvz+va(3,i)
      enddo
      do i=1,natm
        va(1,i)=va(1,i)-sumvx/dble(natm)
        va(2,i)=va(2,i)-sumvy/dble(natm)
        va(3,i)=va(3,i)-sumvz/dble(natm)
      enddo

      tmp=0d0
      do i=1,natm
        is= int(tag(i))
        tmp= tmp +0.5d0*am*(va(1,i)**2 +va(2,i)**2 +va(3,i)**2)
      enddo
      write(6,'(a,es12.4)') " ekin=",tmp
      write(6,'(a,es12.4)') " temp.=",tmp*2d0/3d0/fkb/natm

      end subroutine setv
c=======================================================================
      subroutine myrnd(rnd,dseed)
      real*8 rnd,dseed
      real*8 d2p31m,d2p31
      save d2p31m,d2p31
      data d2p31m/2147483647d0/
      data d2p31 /2147483648d0/
      
      dseed=dmod(16807d0*dseed,d2p31m)
      rnd=dseed/d2p31
      return
      end subroutine myrnd
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make 10mkconf"
c     End:
