      program mkconf_H2
      implicit none
      include './params_unit.h'
      include './params_EAM_Fe.h'
      integer:: is,ifmv,natm,i,ia,ib
      real(8):: h(3,3,0:1),ra(3,2),va(3,2),eki(2),epi(2),strs(3,3,2)
     &     ,tag(2),al
      real(8),parameter:: tempini= 100d0

      al= 2.8553d0
      h(1:3,1:3,0)= 0d0
      h(1,1,0)= al*3d0
      h(2,2,0)= al*3d0
      h(3,3,0)= al*3d0

      is= 2
      ifmv= 1
      ra(1:3,1)= 0.5d0
      tag(1)= 1d0*is +0.1d0*ifmv +1d-14*1
      ra(1:3,2)= ra(1:3,1)
      ra(1,2)= ra(1,2) +0.74d-10 /aa /h(1,1,0)
      tag(2)= 1d0*is +0.1d0*ifmv +1d-14*2
      natm= 2

      call setv(natm,va,tag,tempini)
c-----scale velocities to reduced unit
      do i=1,natm
        va(1,i)=va(1,i) /h(1,1,0)
        va(2,i)=va(2,i) /h(2,2,0)
        va(3,i)=va(3,i) /h(3,3,0)
      enddo

      call write_pmd(15,'pmd000-000','replace',natm,tag,ra,va,h
     &     ,eki,epi,strs)
      
c-----output 'akr000' for Akira visualization
      open(15,file='akr000',form='formatted',status='replace')
      write(15,'(i10,3i5)') natm, 3, 0, 0
      write(15,'(3es11.3)') ((h(ia,ib,0),ia=1,3),ib=1,3)
      do i=1,natm
        write(15,'(i3,6es11.3)') int(tag(i)),ra(1:3,i),va(1:3,i)
      enddo
      close(15)

      end program mkconf_H2
c=======================================================================
      subroutine setv(natm,va,tag,tempini)
      implicit none
      include "./params_unit.h"
      include "./params_EAM_Fe-H.h"
      integer,intent(in):: natm
      real(8),intent(in):: tempini,tag(natm)
      real(8),intent(out):: va(3,natm)

      integer:: i,l,is
      real(8):: dseed,sumvx,sumvy,sumvz,rnd1,rnd2,tmp,facv,am
      real(8),parameter:: pi = 3.14159265358979d0

      facv=dsqrt(2d0*tempini*fkb/am_h)
      am= am_h

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
c     compile-command: "ifort mkconf_H2.F util_pmd.F"
c     End:
