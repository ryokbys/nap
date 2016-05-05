      program mkconf_fcc001
c
c    Create a FCC crystal with (001)-surface
c    Write: ini000
c
      implicit real*8(a-h,o-z)
c-----max # of atoms
      integer,parameter::nmax=10000000
c-----# of unit cells
      integer,parameter:: nuc(1:3)= (/ 100, 100, 4 /)
c-----vacuum thickness in Angstrom
      real(8),parameter:: vacuum = 10d0

      real(8):: ua(3,10)
      real(8):: ra(3,nmax)
      real(8):: h(3,3,0:1),hunit
      real(8):: va(3,nmax)
      integer:: is(nmax),ifmv(nmax)
      real*8,parameter:: pi=3.141592653589793d0
      real*8,parameter:: umass=1822.9d0
      real*8,parameter:: tempau=3.15785d5
cc-----for Ar
c      real*8,parameter:: sgmlj = 3.41d0

c-----velocity scaling factor at Temp(K)=100
      temp=100d0/tempau
cc-----mass of Ar
c      am=39.986d0*umass
c-----mass of Al
      am=26.9815d0*umass
      facv=dsqrt(2d0*temp/am)

      small=1d-7

cc-----box size for Ar
c      cunit= 2d0**(1d0/6) *sgmlj *dsqrt(2d0)
c-----box size for Al
      hunit= 4.025d0

c-----unit cell, FCC
      ua(1:3,1)= (/ 0.0d0, 0.0d0, 0.0d0 /)
      ua(1:3,2)= (/ 0.5d0, 0.5d0, 0.0d0 /)
      ua(1:3,3)= (/ 0.0d0, 0.5d0, 0.5d0 /)
      ua(1:3,4)= (/ 0.5d0, 0.0d0, 0.5d0 /)
      
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
              ra(1,inc)= (ua(1,m)+dble(ix))/nuc(1) +small
              ra(2,inc)= (ua(2,m)+dble(iy))/nuc(2) +small
              ra(3,inc)= (ua(3,m)+dble(iz))/nuc(3) +small
              is(inc)= 1
              ifmv(inc)= 1
c-------------fix bottom atoms
              if(iz.eq.0 .and. m.le.2) ifmv(inc)= 0
            enddo
          enddo
        enddo
      enddo
c-----simulation box size
      h(1:3,1:3,0:1)= 0d0
      h(1,1,0)= hunit*nuc(1)
      h(2,2,0)= hunit*nuc(2)
      h(3,3,0)= hunit*nuc(3)

c-----absolute coordinate only for z
      do i=1,inc
c        ra(1,i)= ra(1,i)*h(1,1,0)
c        ra(2,i)= ra(2,i)*h(2,2,0)
        ra(3,i)= ra(3,i)*h(3,3,0)
      enddo
c-----add water-like molecule only for visualization
      inc= inc+1
      is(inc)= 2
      ra(1,inc)= 0.5d0
      ra(2,inc)= 0.25d0
      ra(3,inc)= h(3,3,0) +5.5d0
      inc= inc+1
      is(inc)= 3
      ra(1,inc)= 0.5d0 +vacuum/10 /h(1,1,0)
      ra(2,inc)= 0.25d0 -vacuum/10 /h(2,2,0)
      ra(3,inc)= h(3,3,0) +4.d0
      inc= inc+1
      is(inc)= 3
      ra(1,inc)= 0.5d0 -vacuum/10 /h(1,1,0)
      ra(2,inc)= 0.25d0 -vacuum/10 /h(2,2,0)
      ra(3,inc)= h(3,3,0) +4.d0
c-----add vacuum
      h(3,3,0)=h(3,3,0) +vacuum
      do i=1,inc
c        ra(1,i)= ra(1,i)/h(1,1,0)
c        ra(2,i)= ra(2,i)/h(2,2,0)
        ra(3,i)= ra(3,i)/h(3,3,0)
      enddo

      
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
      
c-----output 'ini000' for MD program
      open(15,file='pmd0000',form='formatted',status='replace')
      write(15,'(es15.7)') hunit
      write(15,'(3es22.14)') (((h(ia,ib,l)/hunit,ia=1,3),ib=1,3),l=0,1)
      write(15,'(i10)')inc
      do i=1,inc
        write(15,'(2i3,6es22.14)') is(i),ifmv(i),ra(1:3,i),va(1:3,i)
      enddo
      close(15)
      
c-----output 'akr0000' for KVS visualization
      open(15,file='akr0000',form='formatted',status='replace')
      write(15,'(es15.7)') hunit
      write(15,'(3es11.3)') (((h(ia,ib,l)/hunit,ia=1,3),ib=1,3),l=0,1)
      write(15,'(i10,3i5)')inc, 3, 0, 0
      do i=1,inc
        write(15,'(i3,6es11.3)') is(i),ra(1:3,i),va(1:3,i)
      enddo
      close(15)
      
      end program mkconf_fcc001
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
c     compile-command: "gfortran -o mkconf mkconf_fcc.f"
c     End:
