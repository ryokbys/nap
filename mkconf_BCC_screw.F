      program mkconf_BCC
c-----------------------------------------------------------------------
c  Create a BCC crystal for screw dislocation calculation.
c    a1=[11-2]
c    a2=[-110]
c    a3=[111]
c  Only one screw dislocation is included. This configuration may be
c  wrong boundary condition rather than dislocation dipole config.
c-----------------------------------------------------------------------
c  Output
c  ------
c    * pmd000-000
c    * akr000
c-----------------------------------------------------------------------
      implicit none
      include './params_au.h'
      include './params_RK_Fe-H.h'
c-----max # of atoms
      integer,parameter::nmax=10000000
c-----vacuum width in unit of cell
c      integer,parameter:: nvac(1:3)= (/ 5, 5, 5 /)
      integer,parameter:: nvac(1:3)= (/ 0, 1, 0 /)
      integer:: ix,iy,iz,nh,i,m,is,ifmv,ia,ib,natm,itmp,nuc(3)
      real(8):: ua(3,12)
      real(8):: tag(nmax),ra(3,nmax),va(3,nmax),eki(nmax),epi(nmax)
     &     ,h(3,3,0:1),strs(3,3,nmax),h0(3,3),s(3),ymax,ymin,dseed,rnd
     &     ,ra0(3,nmax),rcx,rcy,rcu(3),dmax,cunit,daa,dyl
c.....Parameters
      real(8),parameter:: tempini = 1d0
      real(8),parameter:: small   = 1d-5

c.....Lattice constant of Fe, see Philos. Mag. 83 (2003) 3977
      cunit= alcfe
      daa = sqrt(3d0)/2 *alcfe
      write(6,'(a,es12.4,a)') " Lattice constant =",cunit," [Bohr]"
      write(6,'(a,es12.4,a)') " Fe-Fe bond length=",daa," [Bohr]"

c.....System size
 10   write(6,'(a)') ' Input system size nx,ny,nz:'
      read(5,*) nuc(1:3)

cc.....Selection about screw dislocation configuration
c 20   continue
c      write(6,'(a)') '  1: one screw with PBC on x & z'
c      write(6,'(a)') '  2: two screws with PBC on x,y,z'
c      write(6,'(a)') '  3: one screw with linear elastic'
c     &     //' solution BC'
c      read(5,*) itmp
c      if( itmp.ne.1 .and. itmp.ne.2 .and. itmp.ne.3 ) goto 20
c      write(6,'(a,i2)') ' Selected: ',itmp
      itmp= 1

c.....Unit vectors, h0= (a0,b0,c0) where a0,b0,c0 are column vectors
      h0(1:3,1:3)= 0d0
      h0(1,1)=  2.44948974278318d0 *cunit
      h0(2,2)=  1.41421356237310d0 *cunit
      h0(3,3)=  1.73205080756888d0 *cunit

c.....Atom positions in the unit cell
      ua(1:3,1) =(/ 8.33333333333d-01, 5.000d-01, 1.6666666666d-01 /)
      ua(1:3,2) =(/ 3.33333333333d-01, 0.000d+00, 1.6666666666d-01 /)
      ua(1:3,3) =(/ 5.00000000000d-01, 5.000d-01, 0.0000000000d+00 /)
      ua(1:3,4) =(/ 5.00000000000d-01, 5.000d-01, 5.0000000000d-01 /)
      ua(1:3,5) =(/ 6.66666666666d-01, 0.000d+00, 3.3333333333d-01 /)
      ua(1:3,6) =(/ 6.66666666666d-01, 0.000d+00, 8.3333333333d-01 /)
      ua(1:3,7) =(/ 8.33333333333d-01, 5.000d-01, 6.6666666666d-01 /)
      ua(1:3,8) =(/ 0.00000000000d+00, 0.000d+00, 0.0000000000d+00 /)
      ua(1:3,9) =(/ 0.00000000000d+00, 0.000d+00, 5.0000000000d-01 /)
      ua(1:3,10)=(/ 1.66666666666d-01, 5.000d-01, 3.3333333333d-01 /)
      ua(1:3,11)=(/ 1.66666666666d-01, 5.000d-01, 8.3333333333d-01 /)
      ua(1:3,12)=(/ 3.33333333333d-01, 0.000d+00, 6.6666666666d-01 /)

c.....simulation box size
      h(1:3,1:3,0:1)= 0d0
      h(1,1,0)= h0(1,1)*(nuc(1)+nvac(1))
      h(2,2,0)= h0(2,2)*(nuc(2)+nvac(2))
      h(3,3,0)= h0(3,3)*(nuc(3)+nvac(3))

      if( itmp.eq.1 ) then
        call one_screw_pbc(nmax,natm,tag,ra0,ra,ua,12,nuc,nvac,small
     &       ,daa,h)
      else if( itmp.eq.2 ) then
        
      else

      endif
      
      write(6,'(a,i10)') " num of atom=",natm
c      write(6,'(a,i10)') " id of natm=",int(mod(tag(natm)*1d14,1d13))

      call setv(natm,va,tag,tempini)
      
c-----scale velocities to reduced unit
      do i=1,natm
        va(1,i)=va(1,i) /h(1,1,0)
        va(2,i)=va(2,i) /h(2,2,0)
        va(3,i)=va(3,i) /h(3,3,0)
      enddo

c.....output 'akr000' for Akira visualization
      open(15,file='akr000',form='formatted',status='replace')
      write(15,'(i10,3i5)') natm, 3, 0, 0
      write(15,'(3es11.3)') ((h(ia,ib,0),ia=1,3),ib=1,3)
      do i=1,natm
        write(15,'(i3,6es11.3)') int(tag(i)),ra(1:3,i),va(1:3,i)
      enddo
      close(15)

c.....output 'akrref' for Akira visualization
      open(16,file='akrref',form='formatted',status='replace')
      write(16,'(i10,3i5)') natm, 3, 0, 0
      write(16,'(3es11.3)') ((h(ia,ib,0),ia=1,3),ib=1,3)
      do i=1,natm
        write(16,'(i3,6es11.3)') int(tag(i)),ra0(1:3,i)
     &       ,ra(1:3,i)-ra0(1:3,i)
      enddo
      close(16)

c.....output 'pmdref'
      call write_pmd0_ascii(15,'pmdref','replace',natm,tag,ra0,va,h
     &     ,eki,epi,strs)

      do i=1,natm
        if( ra(3,i).ge.1d0 ) ra(3,i)=ra(3,i)-1d0
        if( ra(3,i).lt.0d0 ) ra(3,i)=ra(3,i)+1d0
      enddo

      call write_pmd0_ascii(15,'pmd000-000','replace',natm,tag,ra,va,h
     &     ,eki,epi,strs)
      
      
      end program mkconf_BCC
c=======================================================================
      subroutine setv(natm,va,tag,tempini)
      implicit none
      include "./params_au.h"
      include "./params_RK_Fe-H.h"
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
c=======================================================================
      subroutine one_screw_pbc(nmax,natm,tag,ra0,ra,ua,nua,nuc,nvac
     &     ,small,daa,h)
      implicit none
      include 'params_au.h'
      integer,intent(in):: nmax,nuc(3),nvac(3),nua
      real(8),intent(in):: ua(3,nua),small,daa,h(3,3,0:1)
      integer,intent(out):: natm
      real(8),intent(out):: tag(nmax),ra0(3,nmax),ra(3,nmax)

      integer:: i,m,inc,ix,iy,iz,is,ifmv
      real(8):: dyl,rcu(3),dmax,ymin,ymax,s(3),rcx,rcy,rx,ry,r,theta
      real(8),parameter:: radi  = 20d0 *aa2bohr
      real(8),parameter:: pi    = 3.1415926535d0

c.....y-unit length
      dyl= 1d0/(nuc(2)+nvac(2))
      write(6,'(a,es12.4)') ' dyl=',dyl

c.....Dislocation center in the unit cell
      rcu(1:3)= (/ 1d0/6, 1d0/4, 0d0 /)
c.....max displacement
      dmax= daa/h(3,3,0)
      
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
              inc=inc+1
              if(inc.gt.nmax)then
                write(*,*)' [Error] inc>nmax, inc,nmax=',inc,nmax
                stop
              endif
              is= 1
              ifmv= 1
              ra0(1:3,inc)= s(1:3)
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
        if( ra0(2,i).lt.ymin+dyl ) then
          is=1
          ifmv=2
          tag(i)= 1d0*is +0.1d0*ifmv +1d-14*i
        elseif( ra0(2,i).gt.ymax-dyl ) then
          is=1
          ifmv=2
          tag(i)= 1d0*is +0.1d0*ifmv +1d-14*i
        endif
      enddo

c.....displace to make one screw dislocation at the center
      rcx= (nuc(1)/2 +rcu(1))/nuc(1)
      rcy= (nuc(2)/2-1 +rcu(2))/nuc(2)
      ra(1:3,1:inc)= ra0(1:3,1:inc)
cc.....Fermi function
c      do i=1,inc
c        if( ra(2,i).lt.rcy ) cycle
c        ra(3,i)=ra(3,i) +dmax*1d0/(1d0+exp(-(ra(1,i)-rcx)/0.1d0))
c      enddo
c.....cosine and theta func
      do i=1,inc
        rx= ra(1,i)-rcx
        ry= ra(2,i)-rcy
        if( ry.lt.0d0 ) cycle
        rx= rx*h(1,1,0)
        ry= ry*h(2,2,0)
        r= sqrt(rx*rx+ry*ry)
        if( r.lt.radi ) then
          theta= acos(-rx/r)
          ra(3,i)= ra(3,i) +dmax/pi *theta
        else
          if( rx.le.-radi ) then
            ! do nothing
          elseif( rx.ge.radi) then
            ra(3,i)=ra(3,i) +dmax
          else
            theta= acos(-rx/radi)
            ra(3,i)=ra(3,i) +dmax/pi *theta
          endif
        endif
      enddo

c.....shift up to in y-direction by half of vacuum length
      do i=1,inc
        ra(2,i)=ra(2,i) +dble(nvac(2))/(nuc(2)+nvac(2))/2
        ra0(2,i)=ra(2,i)
      enddo
      
      natm=inc
      return
      end subroutine one_screw_pbc
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make 10mkconf"
c     End:
