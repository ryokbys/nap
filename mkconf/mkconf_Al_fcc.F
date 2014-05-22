      program mkconf_Al_fcc
c-----------------------------------------------------------------------
c  Make atom configuration file of fcc Al
c    Output: ini000
c-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include "../src/params_unit.h"
      include "../src/params_EAM_Al.h"
c-----max # of atoms
      integer,parameter::nmax=1000000
c-----# of unit cells
      integer,parameter:: nuc(1:3)= (/ 3,3,3 /)
c-----vacuum width in unit of cell
      integer,parameter:: nvac(1:3)= (/ 0, 0, 0 /)
      real(8):: ua(3,4)
      real(8):: ra(3,nmax)
      real(8):: h(3,3,0:1)
      real(8):: va(3,nmax)
      real(8):: tag(nmax),eki(nmax),epi(nmax),strs(3,3,nmax)
c-----initial temperature
      real(8),parameter:: tempini= 100d0

      small=1d-7

c-----Lattice constant of this potential= 4.0425 Ang
      al= 4.0425d0

c-----unit cell, wurtzite structure
c-----Al
      ua(1:3,1)= (/ 0.00d0, 0.00d0, 0.00d0 /)
      ua(1:3,2)= (/ 0.50d0, 0.50d0, 0.00d0 /)
      ua(1:3,3)= (/ 0.50d0, 0.00d0, 0.50d0 /)
      ua(1:3,4)= (/ 0.00d0, 0.50d0, 0.50d0 /)
      
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
              ra(3,inc)= (ua(3,m)+dble(iz))/(nuc(3)+nvac(3)) +small
              ra(1:3,inc)= ra(1:3,inc)
     &             +dble(nvac(1:3))/(nuc(1:3)+nvac(1:3))/2
              is=1
              ifmv= 1
              tag(inc)= 1d0*is +0.1d0*ifmv +1d-14*inc
            enddo
          enddo
        enddo
      enddo
      write(6,'(a,i10)') " natm=",inc
c-----simulation box size
      h(1:3,1:3,0:1)= 0d0
      h(1,1,0)= (nuc(1)+nvac(1)) *al
      h(2,2,0)= (nuc(2)+nvac(2)) *al
      h(3,3,0)= (nuc(3)+nvac(3)) *al

      call setv(inc,va,tag,tempini)

c-----scale velocities to reduced unit
      do i=1,inc
        va(1,i)=va(1,i) /h(1,1,0)
        va(2,i)=va(2,i) /h(2,2,0)
        va(3,i)=va(3,i) /h(3,3,0)
      enddo
      
      call write_pmd0_ascii(15,'pmd00000','replace',inc,tag,ra,va,h
     &     ,al,eki,epi,strs)
      
c-----output 'akr0000' for KVS visualization
      call write_akr(15,'akr0000',inc,h,cunit,tag,ra,va)
      
      end program mkconf_Al_fcc
c=======================================================================
      subroutine setv(natm,va,tag,tempini)
      implicit none
      include "../src/params_unit.h"
      include "../src/params_EAM_Al.h"
      integer,intent(in):: natm
      real(8),intent(in):: tempini,tag(natm)
      real(8),intent(out):: va(3,natm)

      integer:: i,l,is
      real(8):: dseed,sumvx,sumvy,sumvz,rnd1,rnd2,tmp,facv,am
      real(8),parameter:: pi = 3.14159265358979d0

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
