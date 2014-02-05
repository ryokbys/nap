      program mkconf_Si_dislocation
c-----------------------------------------------------------------------
c  Make atomic configuration of diamond structure with:
c    x: [-1, 2,-1]
c    y: [ 1, 1, 1]
c    z: [-1, 0, 1]
c  and create dislocation at the center of the xy plane (z-direction).
c-----------------------------------------------------------------------
c  OUTPUT
c  ------
c    * pm00000-0000
c    * akr0000
c-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include "./params_unit.h"
c-----max # of atoms
      integer,parameter::nmax=1000000
c.....2*pi
      real(8),parameter:: pi   = 3.14159265358979d0
      real(8),parameter:: twopi= 2d0 *pi
c.....(1) shuffle or (2) glide
      integer,parameter:: ifshffl = 1
c-----# of unit cells
      integer,parameter:: nuc(1:3)= (/ 30, 30, 2 /)
      real(8):: ua(3,12)
      real(8):: ra(3,nmax)
      real(8):: h(3,3,0:1)
      real(8):: va(3,nmax)
      real(8):: tag(nmax)
      real(8):: eki(nmax),epi(nmax),strs(3,3,nmax)
      real(8):: xi(3),du(3)

      small=1d-7

c.....Cylinder radius in Ang
      radius= 50d0
c.....anu: Poisson ratio
      anu= 0.25d0

      if( ifshffl.eq.1 ) then ! shuffl
c.....Center position in normalized unit
        xc= (5d0/12 +nuc(1)/2)/nuc(1)
        yc= (11d0/24 +nuc(2)/2)/nuc(2)
      else if( ifshffl.eq.2 ) then ! glide
c.....Center position in normalized unit
        xc= (0.5d0 +nuc(1)/2)/nuc(1)
        yc= (15d0/24 +nuc(2)/2)/nuc(2)
      endif

c-----be: equilibrium bond length: 2.35 Ang
      be= 2.35d0
      write(6,'(a,3es12.4)') ' bond length=',be
c.....al: length of an edge of conventional cubic cell
      al= 4d0*be/sqrt(3d0)
      write(6,'(a,3es12.4)') ' cubic diamond edge length=',al

c-----simulation box size
      h(1:3,1:3,0:1)= 0d0
      h(1,1,0)= nuc(1) *al *sqrt(6d0)/2
      h(2,2,0)= nuc(2) *al *sqrt(3d0)
      h(3,3,0)= nuc(3) *al *sqrt(2d0)/2
      write(6,'(a)') ' h-matrix:'
      write(6,'(2x,3es12.4)') h(1,1:3,0)
      write(6,'(2x,3es12.4)') h(2,1:3,0)
      write(6,'(2x,3es12.4)') h(3,1:3,0)

c-----Unit cell of diamond structure with z-axis towrads (111): type2
      ua(1:3, 1)= (/ 0.00d0, 0.00d0, 0.00d0 /)
      ua(1:3, 2)= (/ 0.00d0, 1d0/4,  0.00d0 /)
      ua(1:3, 3)= (/ 1d0/6,  1d0/3,  1d0/2  /)
      ua(1:3, 4)= (/ 1d0/6,  7d0/12, 1d0/2  /)
      ua(1:3, 5)= (/ 1d0/3,  2d0/3,  0.00d0 /)
      ua(1:3, 6)= (/ 1d0/3, 11d0/12, 0.00d0 /)
      ua(1:3, 7)= (/ 1d0/2,  0.00d0, 1d0/2  /)
      ua(1:3, 8)= (/ 1d0/2,  1d0/4,  1d0/2  /)
      ua(1:3, 9)= (/ 2d0/3,  1d0/3,  0.00d0 /)
      ua(1:3,10)= (/ 2d0/3,  7d0/12, 0.00d0 /)
      ua(1:3,11)= (/ 5d0/6,  2d0/3,  1d0/2  /)
      ua(1:3,12)= (/ 5d0/6, 11d0/12, 1d0/2  /)

c.....Burgers vector
      bv= al *sqrt(2d0)/2
c.....Normalized Burgers vector in x- and z-direction
      bvx= bv *sqrt(3d0)/2
      bvz= bv /2
      bvxn= bvx /h(1,1,0)
      bvzn= bvz /h(3,3,0)
      tan04pi= tan(0.4d0*pi)
      write(6,'(a,es12.4)') ' bv   =',bv
      write(6,'(a,es12.4)') ' bvx  =',bvx
      write(6,'(a,es12.4)') ' bvz  =',bvz
      write(6,'(a,es12.4)') ' bvxn =',bvxn
      write(6,'(a,es12.4)') ' bvzn =',bvzn

      inc=0 
      xmin= 1d10
      xmax= 0d0
      do ix=0,nuc(1)-1
        do iy=0,nuc(2)-1
          do iz=0,nuc(3)-1
            do m=1,12
              xi(1)= (ua(1,m)+dble(ix))/nuc(1) +small
              xi(2)= (ua(2,m)+dble(iy))/nuc(2) +small
              xi(3)= (ua(3,m)+dble(iz))/nuc(3) +small
              x= (xi(1)-xc)*h(1,1,0)
              y= (xi(2)-yc)*h(2,2,0)
              r= sqrt(x**2+y**2)
              if( r.gt.radius ) cycle
c.....Introduce dislocation by subtracting one atomic layer from
c.....lower half of the system.
              if( xi(2).lt.yc .and.
     &             xi(1).gt.xc-bvxn/2 .and. xi(1).lt.xc+bvxn/2 ) cycle
c.....Displace only bottom half atoms to remove the space where atoms
c.....are eliminated. (See, my note on 2014-01-30)
              if( xi(2).lt.yc ) then
                du(1:3)= 0d0
                dx= xi(1)-xc
                if( dx.gt. 0.5d0 ) then dx=dx -1d0
                if( dx.le.-0.5d0 ) then dx=dx +1d0
                if( dx.gt.0d0 ) du(1)= +bvxn*dx -bvxn/2
                if( dx.lt.0d0 ) du(1)= +bvxn*dx +bvxn/2
                dy= (xi(2)-yc) *h(2,2,0)
                du(1)= du(1) *atan(dy/2/bv*(-tan04pi))/(pi/2)
                xi(1)= xi(1) +du(1)
              endif
c.....Displace atoms according to dislocation theory,
c.....only atoms far from the core
              du(1:3)= 0d0
c                du(1)= bvx/twopi *(atan(y/x)
c     &               +x*y/2/(1d0-anu)/(x**2+y**2)
c     &               )
c                du(2)= bvx/twopi
c     &               *((1d0-2d0*anu)/4/(1d0-anu)*log(x**2+y**2)
c     &               +(x**2-y**2)/2/(1-anu)/(x**2+y**2))
              theta= asin(y/(x**2+y**2)) +twopi
              if( x.lt.0d0 ) theta= theta -pi
              du(3)= bvz/twopi *theta
     &             *atan(r/2/bv*(-tan04pi))/(pi/2)
              xi(3)= xi(3) +du(3)/h(3,3,0)
c.....Store xmin and xmax
              xmin= min(xmin,xi(1))
              xmax= max(xmax,xi(1))
c.....Add this atom
              inc=inc+1
              if(inc.gt.nmax)then
                write(*,*)'Error inc>nmax',inc,nmax
                stop
              endif
              call pbc(xi(1))
              call pbc(xi(2))
              call pbc(xi(3))
              ra(1:3,inc)= xi(1:3)
              is= 1
              ifmv= 1
              if( r.gt.radius-2d0*bv ) ifmv= 0
              tag(inc)= 1d0*is +0.1d0*ifmv +1d-14*inc
            enddo
          enddo
        enddo
      enddo
      write(6,'(a,i10)') " natm=",inc

c.....Compress bvx/2 along x
      do i=1,inc
        ra(1,i)=ra(1,i) -bvxn/4*(ra(1,i)-xc)/(xmax-xc)
      enddo

      va(1:3,1:inc)= 0d0

      call write_pmd0_ascii(15,'pmd00000-0000','replace',inc,tag
     &     ,ra,va,h,eki,epi,strs)
      
c-----output 'kvs000' for KVS visualization
      open(15,file='akr0000',form='formatted',status='replace')
      write(15,'(i10)')inc
      write(15,'(3es11.3)') (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
      do i=1,inc
        write(15,'(i3,6es11.3)') int(tag(i)),ra(1:3,i),va(1:3,i)
      enddo
      close(15)
      
      end program mkconf_Si_dislocation
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
c=======================================================================
      subroutine pbc(x)
      implicit none
      real(8),intent(inout):: x

      if( x.lt.0d0 ) x= x+1d0
      if( x.ge.1d0 ) x= x-1d0
      return
      end subroutine pbc
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make 10mkconf"
c     End:
