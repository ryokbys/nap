      program mkconf_Si_2lc
c-----------------------------------------------------------------------
c Make atom configuration file of Si bulk of diamond structure.
c 2-lattice-constant artificial Si atoms are included.
c-----------------------------------------------------------------------
c OUTPUT
c ------
c  * ini000
c-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include "./params_au.h"
      include "./params_SW_Si.h"
c-----max # of atoms
      integer,parameter::nmax=1000000
c-----# of unit cells
      integer,parameter:: nuc(1:3)= (/ 5, 5, 6 /)
c-----vacuum width in unit of cell
      integer,parameter:: nvac(1:3)= (/ 0, 0, 0 /)
      real(8):: ua(3,10)
      real(8):: ra(3,nmax)
      real(8):: h(3,3,0:1)
      real(8):: va(3,nmax)
      real(8):: tag(nmax)
c-----initial temperature
      real(8),parameter:: tempini= 300d0

c-----Si mass
      am= am_si
c-----factor to be multiplied to velocity
c      facv=dsqrt(2d0*temp/am)
      facv=dsqrt(2d0*tempini*fkb/am)

      small=1d-7

c-----box size for Si: 5.4309 AA
      cunit= 5.427d-10 /bohr
c      cunit= 4.07d-10 /bohr

c-----unit cell, diamond structure
      ua(1:3,1)= (/ 0.00d0, 0.00d0, 0.00d0 /)
      ua(1:3,2)= (/ 0.50d0, 0.50d0, 0.00d0 /)
      ua(1:3,3)= (/ 0.00d0, 0.50d0, 0.50d0 /)
      ua(1:3,4)= (/ 0.50d0, 0.00d0, 0.50d0 /)
      ua(1:3,5)= (/ 0.25d0, 0.25d0, 0.25d0 /)
      ua(1:3,6)= (/ 0.75d0, 0.75d0, 0.25d0 /)
      ua(1:3,7)= (/ 0.25d0, 0.75d0, 0.75d0 /)
      ua(1:3,8)= (/ 0.75d0, 0.25d0, 0.75d0 /)
      
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
              ra(1,inc)= (ua(1,m)+dble(ix))/(nuc(1)+nvac(1)) +small
              ra(2,inc)= (ua(2,m)+dble(iy))/(nuc(2)+nvac(2)) +small
              ra(3,inc)= (ua(3,m)+dble(iz))/(nuc(3)+nvac(3)) +small
              ra(1:3,inc)= ra(1:3,inc)
     &             +dble(nvac(1:3))/(nuc(1:3)+nvac(1:3))/2
              is= 1
              if( iz.ge.nuc(3)/2 ) then
                is= 2
              endif
              ifmv= 1
              tag(inc)= 1d0*is +0.1d0*ifmv +1d-14*inc
            enddo
          enddo
        enddo
      enddo
      write(6,'(a,i10)') " natm=",inc
c-----simulation box size
      h(1:3,1:3,0:1)= 0d0
      h(1,1,0)= cunit*(nuc(1)+nvac(1))
      h(2,2,0)= cunit*(nuc(2)+nvac(2))
      h(3,3,0)= cunit*(nuc(3)+nvac(3))

      call setv(inc,va,facv)
      tmp=0d0
      do i=1,inc
        tmp= tmp +0.5d0*am*(va(1,i)**2 +va(2,i)**2 +va(3,i)**2)
      enddo
      write(6,'(a,es12.4)') " ekin=",tmp
      write(6,'(a,es12.4)') " temp.=",tmp*2d0/3d0/fkb/inc

c-----scale velocities to reduced unit
      do i=1,inc
        va(1,i)=va(1,i) /h(1,1,0)
        va(2,i)=va(2,i) /h(2,2,0)
        va(3,i)=va(3,i) /h(3,3,0)
      enddo
      
c-----output 'ini000' for MD program
      open(15,file='ini000',form='formatted',status='replace')
      write(15,'(i10)')inc
      write(15,'(3es22.14)') (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
      do i=1,inc
        write(15,'(7es22.14)') tag(i),ra(1:3,i),va(1:3,i)
      enddo
      close(15)
      
c-----output 'kvs000' for KVS visualization
      open(15,file='kvs000',form='formatted',status='replace')
      write(15,'(i10)')inc
      write(15,'(3es11.3)') (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
      do i=1,inc
        write(15,'(i3,6es11.3)') int(tag(i)),ra(1:3,i),va(1:3,i)
      enddo
      close(15)
      
      end program mkconf_Si_2lc
c=======================================================================
      subroutine setv(natm,va,facv)
      implicit none
      integer,intent(in):: natm
      real(8),intent(in):: facv
      real(8),intent(out):: va(3,natm)

      integer:: i,l
      real(8):: dseed,sumvx,sumvy,sumvz,rnd1,rnd2
      real(8),parameter:: pi = 3.14159265358979d0
      
c-----velocities in Maxwell-Boltzmann distribution
      dseed=12345
      do i=1,natm
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
