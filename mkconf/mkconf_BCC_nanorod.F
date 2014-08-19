      program mkconf_bcc_nanorod
c-----------------------------------------------------------------------
c  Make atom configuration file of bcc-Fe nanorod
c-----------------------------------------------------------------------
c  OUTPUT
c  ------
c    - pmd00000
c-----------------------------------------------------------------------
      implicit none
      include "../pmd/params_unit.h"
c-----max # of atoms
      integer,parameter::nmax=1000000
c-----# of unit cells
      integer,parameter:: nuc(1:3)= (/ 30, 30, 50 /)
c-----vacuum width in unit of cell
      integer,parameter:: nvac(1:3)= (/ 0, 0, 20 /)
c.....Radius of base of the cylinder
      real(8),parameter:: rad  = 30.d0
c.....Hydrogen occupancy
      real(8),parameter:: hratio = 0.01d0
c.....bcc-Fe lattice constant
      real(8),parameter:: alcfe= 2.8553d0
      

      integer:: i,ix,iy,iz,m,inc,is,ifmv,in
      real(8):: small,al,rad2,rx,ry,rr,rnd,dseed,x,y,z,rad2plus,hunit
      real(8):: ua(3,10),uh(3,2)
      real(8):: ra(3,nmax)
      real(8):: h(3,3,0:1)
      real(8):: va(3,nmax)
      real(8):: tag(nmax),epi(nmax),eki(nmax),strs(3,3,nmax)

c.....In this purpose nuc(1)==nuc(2)
      if( nuc(1).ne.nuc(2) )
     &     stop '[Error] nuc(1).ne.nuc(2) !!!'

c.....Selection about Hydrogen inclusion
 10   continue
      write(6,'(a)') '  1: without hydrogens'
      write(6,'(a)') '  2: hydrogens inside nanorod'
      write(6,'(a)') '  3: hydrogens at the surface'
      write(6,'(a)') ' Select one of the above 3:'
      read(5,*) in
      if( in.ne.1 .and. in.ne.2 .and. in.ne.3 ) goto 10
      write(6,'(a,i3)') ' Selected:',in

      small=1d-7

c-----box size for Fe: alcfe= 2.8553 A
      hunit= alcfe

c-----simulation box size
      h(1:3,1:3,0:1)= 0d0
      h(1,1,0)= (nuc(1)+nvac(1)) *hunit
      h(2,2,0)= (nuc(2)+nvac(2)) *hunit
      h(3,3,0)= (nuc(3)+nvac(3)) *hunit
c-----unit cell, BCC
      ua(1:3,1)= (/ 0.00d0, 0.00d0, 0.00d0 /)
      ua(1:3,2)= (/ 0.50d0, 0.50d0, 0.50d0 /)
c.....Hydrogen position in bcc, 1:T-site, 2:O-site
      uh(1:3,1)= (/ 0.25d0, 0.50d0, 0.00d0 /)
      uh(1:3,2)= (/ 0.50d0, 0.50d0, 0.00d0 /)

      rad2= rad*rad
      
      inc=0 
      do ix=0,nuc(1)-1
        do iy=0,nuc(2)-1
          do iz=0,nuc(3)-1
            do m=1,2
              x= (ua(1,m)+dble(ix))/(nuc(1)+nvac(1)) +small
              y= (ua(2,m)+dble(iy))/(nuc(2)+nvac(2)) +small
              z= (ua(3,m)+dble(iz)+dble(nvac(3))/2)
     &             /(nuc(3)+nvac(3)) +small
c.....Only atoms inside the radius remain to form cylinder
              rx= (x-0.5d0)*h(1,1,0)
              ry= (y-0.5d0)*h(2,2,0)
              rr= rx**2 +ry**2
              if( rr.gt.rad2 ) cycle
              inc=inc+1
              if(inc.gt.nmax)then
                write(*,*)'Error inc>nmax',inc,nmax
                stop
              endif
              ra(1,inc)= x
              ra(2,inc)= y
              ra(3,inc)= z
              is= 1
              ifmv= 1
              if( iz.eq.nuc(3)-1 ) then ! top layer
                ifmv= 2
              else if( iz.eq.0 ) then ! bottom layer
                ifmv= 2
              endif
              tag(inc)= 1d0  *is
     &                 +1d-1 *ifmv
     &                 +1d-14*inc
            enddo
          enddo
        enddo
      enddo

c.....Hydrogens
      if( in.eq.2 ) then
        dseed= 11627
        do ix=0,nuc(1)-1
          do iy=0,nuc(2)-1
            do iz=0,nuc(3)-1
              call myrnd(rnd,dseed)
              if( rnd.gt.hratio ) cycle
              x= (uh(1,1)+dble(ix))/(nuc(1)+nvac(1)) +small
              y= (uh(2,1)+dble(iy))/(nuc(2)+nvac(2)) +small
              z= (uh(3,1)+dble(iz)+dble(nvac(3))/2)
     &             /(nuc(3)+nvac(3)) +small
c.....Only inside atoms the radius remain to form cylinder
              rx= (x-0.5d0)*h(1,1,0)
              ry= (y-0.5d0)*h(2,2,0)
              rr= rx**2 +ry**2
              if( rr.gt.rad2 ) cycle
              inc=inc+1
              ra(1,inc)= x
              ra(2,inc)= y
              ra(3,inc)= z
              is= 2
              ifmv= 1
              tag(inc)= 1d0  *is
     &                 +1d-1 *ifmv
     &                 +1d-14*inc
            enddo
          enddo
        enddo
      elseif( in.eq.3 ) then
        rad2plus= (rad+2d0)**2
        dseed= 11627
        do ix=0,nuc(1)-1
          do iy=0,nuc(2)-1
            do iz=0,nuc(3)-1
              call myrnd(rnd,dseed)
              if( rnd.gt.hratio ) cycle
              x= (uh(1,1)+dble(ix))/(nuc(1)+nvac(1)) +small
              y= (uh(2,1)+dble(iy))/(nuc(2)+nvac(2)) +small
              z= (uh(3,1)+dble(iz)+dble(nvac(3))/2)
     &             /(nuc(3)+nvac(3)) +small
c.....Only inside atoms the radius remain to form cylinder
              rx= (x-0.5d0)*h(1,1,0)
              ry= (y-0.5d0)*h(2,2,0)
              rr= rx**2 +ry**2
              if( rr.le.rad2 .or. rr.gt.rad2plus) cycle
              inc=inc+1
              ra(1,inc)= x
              ra(2,inc)= y
              ra(3,inc)= z
              is= 2
              ifmv= 1
              tag(inc)= 1d0  *is
     &                 +1d-1 *ifmv
     &                 +1d-14*inc
            enddo
          enddo
        enddo
      endif

      write(6,'(a,i10)') " natm=",inc

      va(1:3,1:inc)= 0d0

      call write_pmd0_ascii(15,'pmd00000','replace',inc,tag,ra,va
     &     ,h,hunit,eki,epi,strs)

c.....output 'akr0000' for Akira visualization
      call write_akr(15,'akr0000',inc,h,hunit,tag,ra,va)
      
      end program mkconf_bcc_nanorod
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
c     compile-command: "make mkconf_BCC_nanorod"
c     End:
