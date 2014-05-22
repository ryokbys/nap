      program mkconf_Si_bulk
c-----------------------------------------------------------------------
c Make atom configuration file of Si bulk of diamond structure.
c-----------------------------------------------------------------------
c OUTPUT
c ------
c  * pmd00000
c-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include "../src/params_unit.h"
      include "../src/params_SW_Si.h"
c-----max # of atoms
      integer,parameter::nmax=1000000
c-----# of unit cells
      integer,parameter:: nuc(1:3)= (/ 3,3,3 /)
c-----vacuum width in unit of cell
      integer,parameter:: nvac(1:3)= (/ 0, 0, 0 /)
      real(8):: ua(3,10)
      real(8):: ra(3,nmax)
      real(8):: h(3,3,0:1)
      real(8):: va(3,nmax)
      real(8):: tag(nmax),eki(nmax),epi(nmax),strs(3,3,nmax)
      real(8):: si(3)
c.....pi
      real(8),parameter:: pi= 3.1415926535d0
      real(8),parameter:: twopi= 2d0*pi

      small=0d0

c-----box size for Si: 5.4309 Ang
      cunit= 5.427d0
c      cunit= 4.07d0

c-----unit cell, diamond structure
      ua(1:3,1)= (/ 0.00d0, 0.00d0, 0.00d0 /)
      ua(1:3,2)= (/ 0.50d0, 0.50d0, 0.00d0 /)
      ua(1:3,3)= (/ 0.00d0, 0.50d0, 0.50d0 /)
      ua(1:3,4)= (/ 0.50d0, 0.00d0, 0.50d0 /)
      ua(1:3,5)= (/ 0.25d0, 0.25d0, 0.25d0 /)
      ua(1:3,6)= (/ 0.75d0, 0.75d0, 0.25d0 /)
      ua(1:3,7)= (/ 0.25d0, 0.75d0, 0.75d0 /)
      ua(1:3,8)= (/ 0.75d0, 0.25d0, 0.75d0 /)
c-----simulation box size
      h(1:3,1:3,0:1)= 0d0
      h(1,1,0)= cunit*(nuc(1)+nvac(1))
      h(2,2,0)= cunit*(nuc(2)+nvac(2))
      h(3,3,0)= cunit*(nuc(3)+nvac(3))

      q= pi/cunit 
      inc=0 
      do ix=0,nuc(1)-1
        do iy=0,nuc(2)-1
          do iz=0,nuc(3)-1
            do m=1,8
              inc=inc+1
              if(inc.gt.nmax)then
                write(*,*)'Error inc>nmax',inc,nmax
                stop
              endif
              si(1)= (ua(1,m)+dble(ix))/(nuc(1)+nvac(1)) +small
              si(2)= (ua(2,m)+dble(iy))/(nuc(2)+nvac(2)) +small
              si(3)= (ua(3,m)+dble(iz))/(nuc(3)+nvac(3)) +small
              si(1:3)= si(1:3)
     &             +dble(nvac(1:3))/(nuc(1:3)+nvac(1:3))/2
c              si(2)= si(2)+0.5d0/h(1,1,0)*exp(-(si(1)-0.5d0)**2/0.005d0)
cc     &             *sin(q*(si(1)-0.5d0)*h(1,1,0))
c              write(10,'(2es12.4)') si(1)
c     &             ,0.5d0/h(1,1,0) *exp(-(si(1)-0.5d0)**2/0.005d0)
              ra(1:3,inc)= si(1:3)
              is= 1
              ifmv= 1
              tag(inc)= 1d0*is +0.1d0*ifmv +1d-14*inc
            enddo
          enddo
        enddo
      enddo
      write(6,'(a,i10)') " natm=",inc

      va(1:3,1:inc)= 0d0

      call write_pmd0_ascii(15,'pmd00000','replace',inc,tag
     &     ,ra,va,h,cunit,eki,epi,strs)
      
c-----output for Akira
      open(15,file='akr0000',form='formatted',status='replace')
      write(15,'(es15.7)') cunit
      write(15,'(3es11.3)') ((h(ia,ib,l)/cunit,ia=1,3),ib=1,3)
      write(15,'(i10,3i5)')inc,3,0,0
      do i=1,inc
        write(15,'(i3,6es11.3)') int(tag(i)),ra(1:3,i),va(1:3,i)
      enddo
      close(15)
      
      end program mkconf_Si_bulk
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
c     compile-command: "make mkconf_Si_bulk"
c     End:
