      program mkconf_2d
c-----------------------------------------------------------------------
c  Create 2D triangle lattice system with an edge dislocation
c  by removing one atomic layer from lower half area.
c-----------------------------------------------------------------------
c  OUTPUT: 
c   * ini000
c-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c-----max # of atoms
      integer,parameter::nmax=10000000
c-----# of unit cells
      integer,parameter:: nuc(1:3)= (/ 1000, 200, 1 /)
c-----vacuum width in unit of cell
      integer,parameter:: nvac(1:3)= (/ 0, 0, 0 /)
      real(8):: ua(3,10)
      real(8):: ra(3,nmax),ra0(3,nmax)
      real(8):: h(3,3,0:1),h0(3,3,0:1)
      real(8):: va(3,nmax),strs(3,3,nmax),eki(nmax),epi(nmax)
      real(8):: tag(nmax)
      real*8,parameter:: pi=3.141592653589793d0
      real*8,parameter:: umass=1822.9d0
      real*8,parameter:: tempau=3.15785d5
      real(8),parameter:: tempev= 11604d0
c-----for Ar
      real*8,parameter:: sgmlj = 3.41d0
c.....Variables
      real(8):: wlyr,fct,yfix

c-----velocity scaling factor at Temp(K)=100
      temp=0d0/tempev
c-----mass of Ar
      am=39.986d0*umass
cc-----mass of Al
c      am=26.9815d0*amu
      facv=dsqrt(2d0*temp/am)

      small=1d-7

c-----box size for Ar
      cunit= 2d0**(1d0/6) *sgmlj *0.996d0
c      cunit= 2d0**(1d0/6) *sgmlj !*1.2d0!*dsqrt(2d0)
cc-----box size for Al
cc      cunit= 4.025d0
c      cunit= 2.8461d0 *0.956d0

c-----simulation box size
      h(1:3,1:3,0:1)= 0d0
      h(1,1,0)= cunit*(nuc(1)+nvac(1))
      h(2,2,0)= sqrt(3d0)* cunit*(nuc(2)+nvac(2))
      h(3,3,0)= 20d0 *(nuc(3)+nvac(3))
      write(6,'(a)') " h-matrix:"
      write(6,'(3es15.7)') h(1:3,1:3,0)
      
c-----unit cell, trianglar lattice
      ua(1:3,1)= (/ 0.0d0, 0.0d0, 0.5d0 /)
      ua(1:3,2)= (/ 0.5d0, 0.5d0, 0.5d0 /)

c.....Values for removing one atomic layer
      wlyr= cunit/h(1,1,0)
      fct = h(1,1,0)/(h(1,1,0)-cunit)

c.....Fixed region top and bottom
      yfix= cunit*sqrt(3d0)/h(2,2,0)

      inc=0 
      do ix=0,nuc(1)-1
        do iy=0,nuc(2)-1
          do iz=0,nuc(3)-1
            do m=1,2
              x= (ua(1,m)+dble(ix))/(nuc(1)+nvac(1)) +small
              y= (ua(2,m)+dble(iy))/(nuc(2)+nvac(2)) +small
              z= (ua(3,m)+dble(iz))/(nuc(3)+nvac(3)) +small

c.....Remove one atomic layer from lower half area and
c.....shift other atoms to compensate space where the layer was.
              if( y.lt.0.5d0 ) then
                if( x.ge.1d0-wlyr ) cycle
                x= x*fct
              endif
              
              inc=inc+1
              if(inc.gt.nmax)then
                write(*,*)'Error inc>nmax',inc,nmax
                stop
              endif

c              ra(1:3,inc)= ra(1:3,inc)
c     &             +dble(nvac(1:3))/(nuc(1:3)+nvac(1:3))/2
              ra(1,inc)= x
              ra(2,inc)= y
              ra(3,inc)= z

              is= 1
              ifmv= 1
cc-------------x right and left
c              if(ra(1,inc).lt.0.015d0) then
c                is= 2
c                ifmv= 2
c              elseif( ra(1,inc).gt.0.985d0) then
c                is= 3
c                ifmv= 3
c              endif
c-------------y top and bottom
              if(ra(2,inc).lt.yfix .or. ra(2,inc).gt.1d0-yfix ) then
                is= 2
                ifmv= 2
              endif
c-------------tag
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
      
c-----add vacuum by scaling
c      h(1,1,0)= h(1,1,0) *1.1d0
      h(2,2,0)= h(2,2,0) *1.1d0
      do i=1,inc
c-------scaling
c        ra(1,i)= ra(1,i) /1.1d0
        ra(2,i)= ra(2,i) /1.1d0
c-------shift mass center
c        ra(1,i)= ra(1,i) +(1d0-1d0/1.1d0)/2
        ra(2,i)= ra(2,i) +(1d0-1d0/1.1d0)/2
      enddo

c-----set reference position
      h0(1:3,1:3,0:1)= h(1:3,1:3,0:1)
      ra0(1:3,1:inc)= ra(1:3,1:inc)

cc-----output 'ini000' for MD program
c      open(15,file='ini000',form='formatted',status='replace')
c      write(15,'(i10)') inc
c      write(15,'(3es22.14)') (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
c      do i=1,inc
c        write(15,'(10es22.14)') tag(i),ra(1:3,i),ra0(1:3,i),va(1:3,i)
c      enddo
c      close(15)
      call write_pmd(15,'pmd00000','replace',inc,tag,ra,va,h
     &     ,eki,epi,strs)
      
c-----output '0000' for Akira visualization
      open(15,file='akr0000',form='formatted',status='replace')
      write(15,'(i10,3i4)') inc,3,0,0
      write(15,'(3es12.4)') ((h(ia,ib,0),ia=1,3),ib=1,3)
      do i=1,inc
        write(15,'(i3,6es12.4)') int(tag(i)),ra(1:3,i),va(1:3,i)
      enddo
      close(15)
      
      end program mkconf_2d
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
