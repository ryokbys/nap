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
      include './params_unit.h'
c-----max # of atoms
      integer,parameter::nmax=10000000
c-----vacuum width in unit of cell
c      integer,parameter:: nvac(1:3)= (/ 5, 5, 5 /)
      integer,parameter:: nvac(1:3)= (/ 0, 0, 0 /)
      integer:: ix,iy,iz,nh,i,m,is,ifmv,ia,ib,natm,itmp,nuc(3)
      real(8):: ua(3,12)
      real(8):: tag(nmax),ra(3,nmax),va(3,nmax),eki(nmax),epi(nmax)
     &     ,h(3,3,0:1),strs(3,3,nmax),h0(3,3),s(3),ymax,ymin,dseed,rnd
     &     ,ra0(3,nmax),rcx,rcy,rcu(3),dmax,daa,dyl,hunit
c.....Parameters
      real(8),parameter:: tempini = 1d0
      real(8),parameter:: small   = 1d-8
cc.....lattice constant of bcc-W
c      real(8),parameter:: alcfe   = 3.171d0
c.....lattice constant of bcc-Fe
      real(8),parameter:: alcfe   = 2.835d0

c.....Lattice constant of Fe, see Philos. Mag. 83 (2003) 3977
      hunit= alcfe
      daa = sqrt(3d0)/2 *alcfe
      write(6,'(a,es12.4,a)') " Lattice constant =",hunit," [Ang]"
      write(6,'(a,es12.4,a)') " Fe-Fe bond length=",daa," [Ang]"

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
      itmp= 0

c.....Unit vectors, h0= (a0,b0,c0) where a0,b0,c0 are column vectors
      h0(1:3,1:3)= 0d0
      h0(1,1)=  2.44948974278318d0 *hunit
      h0(2,2)=  1.41421356237310d0 *hunit
      h0(3,3)=  1.73205080756888d0/2 *hunit

c.....Atom positions in the unit cell
      ua(1:3,1) =(/ 0.00000000000d+00, 0.000d+00, 0.0000000000d+00 /)
      ua(1:3,2) =(/ 5.00000000000d-01, 5.000d-01, 0.0000000000d+00 /)
      ua(1:3,3) =(/ 3.33333333333d-01, 0.000d+00, 0.3333333333d+00 /)
      ua(1:3,4) =(/ 8.33333333333d-01, 5.000d-01, 0.3333333333d+00 /)
      ua(1:3,5) =(/ 6.66666666666d-01, 0.000d+00, 0.6666666666d+00 /)
      ua(1:3,6) =(/ 1.66666666666d-01, 5.000d-01, 0.6666666666d+00 /)

c.....simulation box size
      h(1:3,1:3,0:1)= 0d0
      h(1,1,0)= h0(1,1)*(nuc(1)+nvac(1))
      h(2,2,0)= h0(2,2)*(nuc(2)+nvac(2))
      h(3,3,0)= h0(3,3)*(nuc(3)+nvac(3))

      if( itmp.eq.1 ) then
        call one_screw_pbc(nmax,natm,tag,ra0,ra,ua,6,nuc,nvac,small
     &       ,daa,h)
      else if( itmp.eq.2 ) then
        call one_screw_cylinder(nmax,natm,tag,ra0,ra,ua,6
     &       ,nuc,nvac,small,daa,h)
      else
        call bulk(nmax,natm,tag,ra0,ra,ua,6,nuc,nvac,small,daa,h)
      endif
      
      write(6,'(a,i10)') " num of atom=",natm

c.....output 'akr000' for Akira visualization
      open(15,file='akr000',form='formatted',status='replace')
      write(15,'(es15.7)') hunit
      write(15,'(3es11.3)') ((h(ia,ib,0)/hunit,ia=1,3),ib=1,3)
      write(15,'(i10,3i5)') natm, 3, 0, 0
      do i=1,natm
        write(15,'(i3,6es11.3)') int(tag(i)),ra(1:3,i),va(1:3,i)
      enddo
      close(15)

c.....output 'akrref' for Akira visualization
      open(16,file='akrref',form='formatted',status='replace')
      write(16,'(es15.7)') hunit
      write(16,'(3es11.3)') ((h(ia,ib,0)/hunit,ia=1,3),ib=1,3)
      write(16,'(i10,3i5)') natm, 3, 0, 0
      do i=1,natm
        write(16,'(i3,6es11.3)') int(tag(i)),ra0(1:3,i)
     &       ,ra(1:3,i)-ra0(1:3,i)
      enddo
      close(16)

c.....output 'pmdref'
      call write_pmd_ascii(17,'pmdref',natm,h,hunit,tag,ra0,va
     &     ,eki,epi,strs)

      do i=1,natm
        if( ra(3,i).ge.1d0 ) ra(3,i)=ra(3,i)-1d0
        if( ra(3,i).lt.0d0 ) ra(3,i)=ra(3,i)+1d0
      enddo

      call write_pmd_ascii(18,'pmd000-000',natm,h,hunit,tag,ra,va
     &     ,eki,epi,strs)
      
      
      end program mkconf_BCC
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
      include 'params_unit.h'
      integer,intent(in):: nmax,nuc(3),nvac(3),nua
      real(8),intent(in):: ua(3,nua),small,daa,h(3,3,0:1)
      integer,intent(out):: natm
      real(8),intent(out):: tag(nmax),ra0(3,nmax),ra(3,nmax)

      integer:: i,m,inc,ix,iy,iz,is,ifmv
      real(8):: dyl,rcu(3),dmax,ymin,ymax,s(3),rcx,rcy,rx,ry,r,theta,sy
      real(8),parameter:: radi  = 2.5d0
      real(8),parameter:: pi    = 3.1415926535d0

c.....y-unit length
      dyl= 1d0/(nuc(2)+nvac(2))
      write(6,'(a,es12.4)') ' dyl=',dyl

c.....Dislocation center in the unit cell
c.....left
      rcu(1:3)= (/ 2d0/6, 0.1d0/6, 0d0 /)
cc.....right
c      rcu(1:3)= (/ 4d0/6, 1d0/6, 0d0 /)
c.....max displacement
      dmax= daa/h(3,3,0)
      
      inc=0 
      ymin=1d0
      ymax=0d0
      do ix=0,nuc(1)-1
        do iy=0,nuc(2)-1
          do iz=0,nuc(3)-1
            do m=1,nua
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

      ra(1:3,1:inc)= ra0(1:3,1:inc)
c.....displace to make one screw dislocation at the center
      rcx= (nuc(1)/2 +rcu(1))/(nuc(1)+nvac(1))
      rcy= (nuc(2)/2 +rcu(2))/(nuc(2)+nvac(2))
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
c.....PBC
        if( ra(3,i).ge.1d0) ra(3,i)= ra(3,i) -1d0
        if( ra(3,i).lt.0d0) ra(3,i)= ra(3,i) +1d0
      enddo

c.....shift up to in y-direction by half of vacuum length
      sy= dble(nvac(2))/(nuc(2)+nvac(2))/2
      do i=1,inc
        ra(2,i)=ra(2,i) +sy
        ra0(2,i)=ra(2,i)
      enddo

      write(6,'(a,2es12.4)') ' rcx=',rcx
      write(6,'(a,2es12.4)') ' rcy=',rcy+sy
      
      natm=inc
      return
      end subroutine one_screw_pbc
c=======================================================================
      subroutine one_screw_cylinder(nmax,natm,tag,ra0,ra,ua,nua,nuc,nvac
     &     ,small,daa,h)
      implicit none
      include 'params_unit.h'
      integer,intent(in):: nmax,nuc(3),nvac(3),nua
      real(8),intent(in):: ua(3,nua),small,daa,h(3,3,0:1)
      integer,intent(out):: natm
      real(8),intent(out):: tag(nmax),ra0(3,nmax),ra(3,nmax)

      real(8),parameter:: pi  = 3.14159265358979d0

      integer:: ix,iy,iz,m,inc,is,ifmv,i
      real(8):: rvac,rc,rcu(3),rca(3),dmax,ymin,ymax,s(3),dyl
     &     ,rx,ry,r,theta,rcx,rcy,drcx,drcy,x,y,uh(3),rh(3)

      write(6,'(a)') ' h-mat:'
      write(6,'(3es12.4)') h(1,1:3,0)
      write(6,'(3es12.4)') h(2,1:3,0)
      write(6,'(3es12.4)') h(3,1:3,0)

      rvac= max(nvac(1)*h(1,1,0)/(nuc(1)+nvac(1))
     &     ,nvac(2)*h(2,2,0)/(nuc(1)+nvac(1)) )
      write(6,'(a,es12.4)') ' rvac [Ang]=',rvac
      rc= (min(h(1,1,0),h(2,2,0))-rvac)/2
      write(6,'(a,es12.4)') ' rc [Ang]=',rc

c.....y-unit length
      dyl= 1d0/(nuc(2)+nvac(2)) *h(2,2,0)
      write(6,'(a,es12.4)') ' dyl [Ang]=',dyl
c.....max displacement
      dmax= daa/h(3,3,0)

c.....Dislocation center in the unit cell
c.....Left
      rcu(1:3)= (/ 5d0/6, 5d0/6, 0d0 /)
      rcx= (nuc(1)/2-1 +rcu(1))/(nuc(1)+nvac(1))
cc.....Right
c      rcu(1:3)= (/ 1d0/6, 5d0/6, 0d0 /)
c      rcx= (nuc(1)/2 +rcu(1))/(nuc(1)+nvac(1))

c      rcy= (nuc(2)/2 +rcu(2))/(nuc(2)+nvac(2))
      rcy= (nuc(2)/2 -1 +rcu(2))/(nuc(2)+nvac(2))
      write(6,'(a,es12.4)') ' rcx=',rcx
      write(6,'(a,es12.4)') ' rcy=',rcy
      drcx= rcx -0.5d0
      drcy= rcy -0.5d0
c.....max displacement
      dmax= daa/h(3,3,0)

      inc=0 
      ymin=1d0
      ymax=0d0
      do ix=0,nuc(1)-1
        do iy=0,nuc(2)-1
          do iz=0,nuc(3)-1
            do m=1,nua
              s(1)= (ua(1,m)+dble(ix))/(nuc(1)+nvac(1)) +small
              s(2)= (ua(2,m)+dble(iy))/(nuc(2)+nvac(2)) +small
              s(3)= (ua(3,m)+dble(iz))/(nuc(3)+nvac(3)) +small
              x= (s(1)-0.5d0)*h(1,1,0)
              y= (s(2)-0.5d0)*h(2,2,0)
              r= sqrt(x*x +y*y)
              if( r.gt.rc ) cycle
              inc=inc+1
              if(inc.gt.nmax)then
                write(*,*)' [Error] inc>nmax, inc,nmax=',inc,nmax
                stop
              endif
              is= 1
              ifmv= 1
c              if( r.gt.rc-dyl ) ifmv=0
              ra0(1:3,inc)= s(1:3)
              tag(inc)= 1d0*is +0.1d0*ifmv +1d-14*inc
c.....Search top and bottom in y-direction
              ymin=min(ymin,s(2))
              ymax=max(ymax,s(2))
            enddo
          enddo
        enddo
      enddo

c c.....Hydrogen at T-site
c cc.....Hard-core
c c      uh(1:3)= (/ 0d0, 1d0/4, 1d0/2 /)
c c      rh(1)= (nuc(1)/2 +uh(1))/(nuc(1)+nvac(1)) +small
c c.....Easy-core
c       uh(1:3)= (/ 0.0d0, 3d0/4, 1d0/2 /)
c       rh(1)= (nuc(1)/2 +uh(1))/(nuc(1)+nvac(1)) +small
c 
c       rh(2)= (nuc(2)/2 -1 +uh(2))/(nuc(2)+nvac(2)) +small
c       rh(3)= uh(3) +small
c       inc= inc +1
c       is=2
c       ifmv=1
c       ra0(1:3,inc)= rh(1:3)
c       tag(inc)= 1d0*is +0.1d0*ifmv +1d-14*inc
c 
      ra(1:3,1:inc)= ra0(1:3,1:inc)
c.....Displace
      do i=1,inc
        rx= ra(1,i) -rcx
        ry= ra(2,i) -rcy
        r= sqrt(rx*rx +ry*ry)
        theta= acos(rx/r)
        if( asin(ry/r).lt.0 ) theta= 2d0*pi -theta
        ra(3,i)= ra(3,i) +dmax*theta/2/pi
c.....PBC
        if( ra(3,i).ge.1d0) ra(3,i)= ra(3,i) -1d0
        if( ra(3,i).lt.0d0) ra(3,i)= ra(3,i) +1d0
      enddo
      natm= inc
      return
      end subroutine one_screw_cylinder
c=======================================================================
      subroutine bulk(nmax,natm,tag,ra0,ra,ua,nua,nuc,nvac,small,daa,h)
      implicit none
      integer,intent(in):: nmax,nuc(3),nvac(3),nua
      real(8),intent(in):: ua(3,nua),small,daa,h(3,3,0:1)
      integer,intent(out):: natm
      real(8),intent(out):: tag(nmax),ra0(3,nmax),ra(3,nmax)

      integer:: ix,iy,iz,m,is,ifmv
      real(8):: s(3)

      write(6,'(a)') ' h-mat:'
      write(6,'(3es12.4)') h(1,1:3,0)
      write(6,'(3es12.4)') h(2,1:3,0)
      write(6,'(3es12.4)') h(3,1:3,0)

      natm= 0
      do ix=0,nuc(1)-1
        do iy=0,nuc(2)-1
          do iz=0,nuc(3)-1
            do m=1,nua
              s(1)= (ua(1,m)+dble(ix))/(nuc(1)+nvac(1)) +small
              s(2)= (ua(2,m)+dble(iy))/(nuc(2)+nvac(2)) +small
              s(3)= (ua(3,m)+dble(iz))/(nuc(3)+nvac(3)) +small
              natm=natm +1
              is= 1
              ifmv=1
              ra0(1:3,natm)= s(1:3)
              tag(natm)= 1d0*is +0.1d0*ifmv +1d-14*natm
            enddo
          enddo
        enddo
      enddo

      ra(1:3,1:natm)= ra0(1:3,1:natm)
      
      end subroutine bulk
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make 10mkconf"
c     End:
