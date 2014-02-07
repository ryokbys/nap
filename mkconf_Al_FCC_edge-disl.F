      program mkconf_Al_fcc
c-----------------------------------------------------------------------
c  Make atom configuration file of FCC Al in the simulation box
c  of x=(-1,1,0), y=(-1,-1,2), z=(1,1,1).
c-----------------------------------------------------------------------
c  Output: ini000
c-----------------------------------------------------------------------
      implicit real*8(a-h,o-z),integer(i-n)
      include "./params_unit.h"
      include "./params_EAM_Al.h"
c-----max # of atoms
      integer,parameter::nmax=1000000
c-----# of unit cells
      integer,parameter:: nuc(1:3)= (/ 200,3,50 /)
c.....Cutoff radius of the EAM potential
      real(8),parameter:: rc = 7.0d0
      real(8):: ua(3,6)
      real(8):: ra(3,nmax)
      real(8):: h(3,3,0:1),h0(3,3),sorg(3)
      real(8):: va(3,nmax),strs(3,3,nmax),eki(nmax),epi(nmax)
      real(8):: tag(nmax)
c-----initial temperature
      real(8),parameter:: tempini= 100d0

      small=1d-7

c-----Lattice constant of this potential: 4.0170 Ang
      al= 4.017d0
      write(6,'(a,es12.4,a)') " Lattice constant  =",al," [Ang]"
c.....Equilibrium distance between nearest neighbor Al atoms
      daa= al/sqrt(2d0)
      write(6,'(a,es12.4,a)') " Al-Al bond length =",daa," [Ang]"

c.....Unit vectors, h0= (a0,b0,c0) where a0,b0,c0 are column vectors
      h0(1:3,1)= (/ daa, 0d0,           0d0          /)
      h0(1:3,2)= (/ 0d0, daa*sqrt(3d0), 0d0          /)
      h0(1:3,3)= (/ 0d0, 0d0,           daa*sqrt(6d0)/)
c.....Simulation box size
      h(1:3,1:3,0:1)= 0d0
      h(1:3,1,0)= h0(1:3,1) *nuc(1)
      h(1:3,2,0)= h0(1:3,2) *nuc(2)
      h(1:3,3,0)= h0(1:3,3) *nuc(3)
      
c.....Atom positions in the unit cell
      ua(1:3,1)= (/ 0.0d0, 0.0d0, 0.0d0 /)
      ua(1:3,2)= (/ 0.5d0, 0.5d0, 0.0d0 /)
      ua(1:3,3)= (/ 0.5d0, 1d0/6, 1d0/3 /)
      ua(1:3,4)= (/ 0.0d0, 4d0/6, 1d0/3 /)
      ua(1:3,5)= (/ 0.0d0, 2d0/6, 2d0/3 /)
      ua(1:3,6)= (/ 0.5d0, 5d0/6, 2d0/3 /)
      
      tmp1= rc/h(3,3,0)
      tmp2= h0(3,3)/h(3,3,0)
      write(6,'(a,es12.4)') " tmp1 = ",tmp1
      write(6,'(a,es12.4)') " tmp2 = ",tmp2
      inc=0 
      do ix=0,nuc(1)-1
        do iy=0,nuc(2)-1
          do iz=0,nuc(3)-1
            do m=1,6
              inc=inc+1
              if(inc.gt.nmax)then
                write(*,*)'Error inc>nmax',inc,nmax
                stop
              endif
              ra(1,inc)= (ua(1,m)+dble(ix))/nuc(1) +small
              ra(2,inc)= (ua(2,m)+dble(iy))/nuc(2) +small
              ra(3,inc)= (ua(3,m)+dble(iz))/nuc(3) +small
              is=1
              ifmv= 1
              tag(inc)= 1d0*is +0.1d0*ifmv +1d-14*inc
c.....Eliminate top and bottom atoms
              if( ra(3,inc).lt.tmp1 ) then
                inc=inc -1
                cycle
              elseif( 1d0-tmp1.lt.ra(3,inc) ) then
                inc=inc -1
                cycle
              endif
c.....Eliminate atomic plane from half bottom of the system
              if( ra(3,inc).lt.0.5d0 .and. ix.eq.nuc(1)/2 ) then
                inc=inc-1
                cycle
              endif
            enddo
          enddo
        enddo
      enddo
      write(6,'(a,i10)') " natm=",inc
c.....Set ifmv=2 of top and bottom atoms which will be controled
      do i=1,inc
        if( ra(3,i).lt.tmp1+tmp2 ) then
          is=1
          ifmv=2
          tag(i)= 1d0*is +0.1d0*ifmv +1d-14*i
        elseif( 1d0-(tmp1+tmp2).lt.ra(3,i) ) then
          is=1
          ifmv=2
          tag(i)= 1d0*is +0.1d0*ifmv +1d-14*i
        endif
      enddo


      call setv(inc,va,tag,tempini)

c-----scale velocities to reduced unit
      do i=1,inc
        va(1,i)=va(1,i) /h(1,1,0)
        va(2,i)=va(2,i) /h(2,2,0)
        va(3,i)=va(3,i) /h(3,3,0)
      enddo
      
cc-----output 'pmd00000' for MD program
c      open(15,file='pmd00000',form='unformatted',status='replace')
c      write(15,'(i10)') inc
c      write(15,'(3es22.14)') (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
c      do i=1,inc
c        write(15,'(7es22.14)') tag(i),ra(1:3,i),va(1:3,i)
c      enddo
c      close(15)
      call write_pmd(15,'pmd00000','replace',inc,tag,ra,va,h
     &     ,eki,epi,strs)

c-----output 'akr0000' for Akira visualization
      open(15,file='akr0000',form='formatted',status='replace')
      write(15,'(i10,3i5)') inc, 3, 0, 0
      write(15,'(3es11.3)') ((h(ia,ib,0),ia=1,3),ib=1,3)
      do i=1,inc
        write(15,'(i3,6es11.3)') int(tag(i)),ra(1:3,i),va(1:3,i)
      enddo
      close(15)
      
      end program mkconf_Al_fcc
c=======================================================================
      subroutine setv(natm,va,tag,tempini)
      implicit none
      include "./params_unit.h"
      include "./params_EAM_Al.h"
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
