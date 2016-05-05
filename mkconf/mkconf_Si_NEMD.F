      program mkconf_Si_NEMD
c-----------------------------------------------------------------------
c Make atom configuration file of crystalline Si system
c to perform the non-equilibrium MD for thermal conductance calculation.
c The system is long in x-direction which is to be the thermal flow
c direction. Vacuum region is added in x-direction.
c-----------------------------------------------------------------------
c OUTPUT
c ------
c  * pmd0000
c-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include "../pmd/params_unit.h"
c-----max # of atoms
      integer,parameter::nmax=1000000
c-----# of unit cells
      integer,parameter:: nuc(1:3)= (/ 100,5,5 /)
c-----vacuum width in unit of cell
      integer,parameter:: nvac(1:3)= (/ 1, 0, 0 /)
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
              ra(1:3,inc)= si(1:3)
              is= 1
              ifmv= 1
c.....change ifmv of atoms at temperature control region
              if( ix.le.3 ) then
                ifmv= 2
              elseif ( ix.ge.nuc(1)-4 ) then
                ifmv= 3
              endif
c.....fix atoms at leftmost and rightmost face
              if( ix.eq.0 .and.
     &             (m.eq.1 .or. m.eq.3 .or. m.eq.5 .or. m.eq.7) ) then
                ifmv= 0
              else if( ix.eq.nuc(1)-1 .and.
     &               (m.eq.2 .or. m.eq.4 .or. m.eq.6 .or. m.eq.8) ) then
                ifmv= 0
              endif
              tag(inc)= 1d0*is +0.1d0*ifmv +1d-14*inc
            enddo
          enddo
        enddo
      enddo
      write(6,'(a,i10)') " natm=",inc

      va(1:3,1:inc)= 0d0

      call write_pmd0_ascii(15,'pmd0000','replace',inc,tag
     &     ,ra,va,h,cunit,eki,epi,strs)
      
c-----output for Akira
      open(15,file='akr0000',form='formatted',status='replace')
      write(15,'(es15.7)') cunit
      write(15,'(3es11.3)') ((h(ia,ib,l)/cunit,ia=1,3),ib=1,3)
      write(15,'(i10,3i5)') inc,3,0,0
      do i=1,inc
        write(15,'(i3,6es11.3)') int(tag(i)),ra(1:3,i),va(1:3,i)
      enddo
      close(15)
      
      end program mkconf_Si_NEMD
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make mkconf_Si_NEMD"
c     End:
