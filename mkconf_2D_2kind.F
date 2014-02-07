      program mkconf_2D_2kind
c-----------------------------------------------------------------------
c  Create 2D 2-kind system for droplet simulation
c    Write: ini000
c-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include "params_unit.h"
c-----max # of atoms
      integer,parameter::nmax=1000000
c-----# of unit cells
      integer,parameter:: nuc(1:3)= (/ 10, 10, 1 /)
c-----vacuum width in unit of cell
      integer,parameter:: nvac(1:3)= (/ 0, 0, 0 /)
      real(8):: ua(3,10)
      real(8):: tag(nmax),ra(3,nmax),va(3,nmax),eki(nmax),epi(nmax)
     &     ,h(3,3,0:1),strs(3,3,nmax)
      real*8,parameter:: pi=3.141592653589793d0
c-----for Ar
      real*8,parameter:: sgm = 3.41d0

c-----velocity scaling factor at Temp(K)=50
      tempini=50d0
c-----mass of Ar
      am=39.986d0*umass
      facv=dsqrt(2d0*tempini*fkb/am)

      small=1d-7

c-----box size for Ar
      cunit= 2d0**(1d0/6) *sgm *0.996d0

c-----simulation box size
      h(1:3,1:3,0:1)= 0d0
      h(1,1,0)= cunit*(nuc(1)+nvac(1))
      h(2,2,0)= cunit*(nuc(2)+nvac(2)) *sqrt(3d0)
      h(3,3,0)= cunit *5

c-----unit cell, triangular lattice
      ua(1:3,1)= (/ 0.0d0, 0.0d0, 0.5d0 /)
      ua(1:3,2)= (/ 0.5d0, 0.5d0, 0.5d0 /)
      
      inc=0 
      do ix=0,nuc(1)-1
        do iy=0,nuc(2)-1
          do iz=0,nuc(3)-1
            do m=1,2
              inc=inc+1
              if(inc.gt.nmax)then
                write(*,*)'Error inc>nmax',inc,nmax
                stop
              endif
              ra(1,inc)= (ua(1,m)+dble(ix))/(nuc(1)+nvac(1)) +small
              ra(2,inc)= (ua(2,m)+dble(iy))/(nuc(2)+nvac(2)) +small
              ra(3,inc)= (ua(3,m)+dble(iz))/(nuc(3)+nvac(3)) +small
cc-------------shift y by a half of vacuum space
c              ra(2,inc)= ra(2,inc)
c     &             +dble(nvac(2))/(nuc(2)+nvac(2))*0.25d0
              is= 1
c              if( iy.gt.3 ) is= 2 ! top half is species 2
              ifmv= 1
              if( iy.eq.0 ) ifmv= 0 ! fix bottom layer
              tag(inc)= 1d0*is +0.1d0*ifmv +1d-14*inc
            enddo
          enddo
        enddo
      enddo
      write(6,'(a,i10)') " natm=",inc
      write(6,'(a,i10)') " id of inc=",int(mod(tag(inc)*1d14,1d13))
      
c-----velocities in Maxwell-Boltzmann distribution
      dseed=12345
      do i=1,inc
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
      do i=1,inc
        sumvx=sumvx+va(1,i)
        sumvy=sumvy+va(2,i)
        sumvz=sumvz+va(3,i)
      enddo
      do i=1,inc
        va(1,i)=va(1,i)-sumvx/dble(inc)
        va(2,i)=va(2,i)-sumvy/dble(inc)
        va(3,i)=va(3,i)-sumvz/dble(inc)
      enddo
c-----scale velocities to reduced unit
      do i=1,inc
        va(1,i)=va(1,i) /h(1,1,0)
        va(2,i)=va(2,i) /h(2,2,0)
        va(3,i)=va(3,i) /h(3,3,0)
      enddo
      
cc-----output 'ini000' for MD program
c      open(15,file='ini000',form='formatted',status='replace')
c      write(15,'(i10)')inc
c      write(15,'(3es22.14)') (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
c      do i=1,inc
c        write(15,'(7es22.14)') tag(i),ra(1:3,i),va(1:3,i)
c      enddo
c      close(15)
      call write_pmd(15,'pmd00000','replace',inc,tag,ra,va,h
     &     ,eki,epi,strs)

      
c-----output 'akr0000' for Akira visualization
      open(15,file='akr0000',form='formatted',status='replace')
      write(15,'(i10,3i5)') inc,3,0,0
      write(15,'(3es11.3)') ((h(ia,ib,0),ia=1,3),ib=1,3)
      do i=1,inc
        write(15,'(i3,6es11.3)') int(tag(i)),ra(1:3,i),va(1:3,i)
      enddo
      close(15)
      
      end program mkconf_2D_2kind
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
