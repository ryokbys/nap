      program mkconf_BCC
c-----------------------------------------------------------------------
c  Create a BCC crystal for edge dislocation calculation.
c-----------------------------------------------------------------------
c  Output
c  ------
c    * pmd000-000
c    * akr000
c-----------------------------------------------------------------------
      implicit real*8(a-h,o-z),integer(i-n)
      include './params_au.h'
      include './params_EAM_Fe-H.h'
c-----max # of atoms
      integer,parameter::nmax=10000000
c-----# of unit cells
c      integer,parameter:: nuc(1:3)= (/ 1, 1, 1 /)
      integer,parameter:: nuc(1:3)= (/ 50, 25, 20 /)
c-----vacuum width in unit of cell
c      integer,parameter:: nvac(1:3)= (/ 5, 5, 5 /)
      integer,parameter:: nvac(1:3)= (/ 0, 1, 0 /)
      integer:: ix,iy,iz,inc,nh
      real(8):: ua(3,12)
      real(8):: tag(nmax),ra(3,nmax),va(3,nmax),eki(nmax),epi(nmax)
     &     ,h(3,3,0:1),strs(3,3,nmax),h0(3,3),s(3),ymax,ymin,dseed,rnd
c.....Parameters
      real(8),parameter:: tempini = 100d0
c.....H density in atomic ppm
      real(8),parameter:: hppm    = 200d0

      small=1d-8

c.....Lattice constant of Fe, see Philos. Mag. 83 (2003) 3977
      cunit= alcfe
      daa = sqrt(3d0)/2 *alcfe
      write(6,'(a,es12.4,a)') " Lattice constant =",cunit," [Bohr]"
      write(6,'(a,es12.4,a)') " Fe-Fe bond length=",daa," [Bohr]"

c.....Unit vectors, h0= (a0,b0,c0) where a0,b0,c0 are column vectors
      h0(1:3,1:3)= 0d0
      h0(1,1)=  1.73205080756888d0 *cunit
      h0(2,2)=  2.44948974278318d0 *cunit
      h0(3,3)=  1.41421356237310d0 *cunit

c.....Atom positions in the unit cell
      ua(1:3,1) =(/ 1.666666d-01, 3.333333d-01, 0.000d+00 /)
      ua(1:3,2) =(/ 0.000000d+00, 0.000000d+00, 0.000d+00 /)
      ua(1:3,3) =(/ 5.000000d-01, 0.000000d+00, 0.000d+00 /)
      ua(1:3,4) =(/ 1.666666d-01, 8.333333d-01, 5.000d-01 /)
      ua(1:3,5) =(/ 0.000000d+00, 5.000000d-01, 5.000d-01 /)
      ua(1:3,6) =(/ 5.000000d-01, 5.000000d-01, 5.000d-01 /)
      ua(1:3,7) =(/ 3.333333d-01, 1.666666d-01, 5.000d-01 /)
      ua(1:3,8) =(/ 8.333333d-01, 1.666666d-01, 5.000d-01 /)
      ua(1:3,9) =(/ 3.333333d-01, 6.666666d-01, 0.000d+00 /)
      ua(1:3,10)=(/ 8.333333d-01, 6.666666d-01, 0.000d+00 /)
      ua(1:3,11)=(/ 6.666666d-01, 3.333333d-01, 0.000d+00 /)
      ua(1:3,12)=(/ 6.666666d-01, 8.333333d-01, 5.000d-01 /)

c.....simulation box size
      h(1:3,1:3,0:1)= 0d0
      h(1,1,0)= h0(1,1)*(nuc(1)+nvac(1))
      h(2,2,0)= h0(2,2)*(nuc(2)+nvac(2))
      h(3,3,0)= h0(3,3)*(nuc(3)+nvac(3))

c.....Width of to-be-removed half plane
      daah= daa/h(1,1,0) +small
      
      nh= 0
      inc=0 
      ymin=1d0
      ymax=0d0
      do ix=0,nuc(1)-1
        do iy=0,nuc(2)-1
          do iz=0,nuc(3)-1
            do m=1,12
              s(1)= (ua(1,m)+dble(ix))/(nuc(1)+nvac(1)) +small
              s(2)= (ua(2,m)+dble(iy))/(nuc(2)+nvac(2)) +small
              s(3)= (ua(3,m)+dble(iz))/(nuc(3)+nvac(3)) +small
              s(1:3)= s(1:3)
     &             +dble(nvac(1:3))/(nuc(1:3)+nvac(1:3))/2
c.....Remove an atomic plane along yz-plane of bottom half
              if( s(2).lt.0.5d0 .and.
     &             s(1).ge.0.5d0-daah/2 .and.
     &             s(1).lt.0.5d0+daah/2 ) then
c.....Add H at the dislocation core
                if( s(2).gt.0.49d0 ) then
                  nh= nh +1
                  inc=inc+1
                  is= 2
                  ifmv= 1
                  ra(1:3,inc)= s(1:3)
                  tag(inc)= 1d0*is +0.1d0*ifmv +1d-14*inc
                endif
                cycle
              endif
cc.....Shift to bury removed area
c              if( s(2).lt.0.5d0 ) then
c                if( s(1).lt.0.5d0 ) then
c                  s(1)=s(1) +daah/2*s(1)/0.5d0
c                else
c                  s(1)=s(1) -daah/2*(1d0-s(1))/0.5d0
c                endif
c              endif
              inc=inc+1
              if(inc.gt.nmax)then
                write(*,*)' [Error] inc>nmax, inc,nmax=',inc,nmax
                stop
              endif
              is= 1
              ifmv= 1
              ra(1:3,inc)= s(1:3)
              tag(inc)= 1d0*is +0.1d0*ifmv +1d-14*inc
c.....Search top and bottom in y-direction
              ymin=min(ymin,s(2))
              ymax=max(ymax,s(2))
            enddo
          enddo
        enddo
      enddo
c.....Set top and bottom atoms to ifmv=2
      do i=1,inc
        if( ra(2,i).lt.ymin+daah ) then
          is=1
          ifmv=2
          tag(i)= 1d0*is +0.1d0*ifmv +1d-14*i
        elseif( ra(2,i).gt.ymax-daah ) then
          is=1
          ifmv=2
          tag(i)= 1d0*is +0.1d0*ifmv +1d-14*i
        endif
      enddo

cc.....Add H atoms to O-sites
c      dseed=1234589d0
c      do ix=0,nuc(1)-1
c        do iy=0,nuc(2)-1
c          do iz=0,nuc(3)-1
c            do m=1,12
c              call myrnd(rnd,dseed)
c              if( rnd.le.hppm/1d+6 ) then
c                nh=nh +1
c                inc=inc +1
c                is= 2
c                ifmv=1
c                tag(inc)= 1d0*is +0.1d0*ifmv +1d-14*inc
c                ra(1,inc)= (0.5d0 +dble(ix))/(nuc(1)+nvac(1)) +small
c                ra(2,inc)= (0.5d0 +dble(iy))/(nuc(2)+nvac(2)) +small
c                ra(3,inc)= (1.0d0/3 +dble(iz))/(nuc(3)+nvac(3)) +small
c              endif
c            enddo
c          enddo
c        enddo
c      enddo

      write(6,'(a,i10)') " num of atom=",inc
      write(6,'(a,i10)') " num of H   =",nh
c      write(6,'(a,i10)') " id of inc=",int(mod(tag(inc)*1d14,1d13))

      call setv(inc,va,tag,tempini)
      
c-----scale velocities to reduced unit
      do i=1,inc
        va(1,i)=va(1,i) /h(1,1,0)
        va(2,i)=va(2,i) /h(2,2,0)
        va(3,i)=va(3,i) /h(3,3,0)
      enddo

      call write_pmd(15,'pmd000-000','replace',inc,tag,ra,va,h
     &     ,eki,epi,strs)
      
c-----output 'akr000' for Akira visualization
      open(15,file='akr000',form='formatted',status='replace')
      write(15,'(i10,3i5)') inc, 3, 0, 0
      write(15,'(3es11.3)') ((h(ia,ib,0),ia=1,3),ib=1,3)
      do i=1,inc
        write(15,'(i3,6es11.3)') int(tag(i)),ra(1:3,i),va(1:3,i)
      enddo
      close(15)
      
      end program mkconf_BCC
c=======================================================================
      subroutine setv(natm,va,tag,tempini)
      implicit none
      include "./params_au.h"
      include "./params_EAM_Fe-H.h"
      integer,intent(in):: natm
      real(8),intent(in):: tempini,tag(natm)
      real(8),intent(out):: va(3,natm)

      integer:: i,l,is
      real(8):: dseed,sumvx,sumvy,sumvz,rnd1,rnd2,tmp,facv(2),am(2)
      real(8),parameter:: pi = 3.14159265358979d0

      facv(1)=dsqrt(2d0*tempini*fkb/am_fe)
      facv(2)=dsqrt(2d0*tempini*fkb/am_h)
      am(1)= am_fe
      am(2)= am_h

c-----velocities in Maxwell-Boltzmann distribution
      dseed=12345
      do i=1,natm
        is= int(tag(i))
        do l=1,3
          call myrnd(rnd1,dseed)
          call myrnd(rnd2,dseed)
          va(l,i)=facv(is)*dsqrt(-dlog(rnd1))*dcos(2d0*pi*rnd2)
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
        tmp= tmp +0.5d0*am(is)*(va(1,i)**2 +va(2,i)**2 +va(3,i)**2)
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
