      program mkconf_WHe_compress
c-----------------------------------------------------------------------
c  Create a BCC tungsten(W) crystal, He gas and vacuum region.
c-----------------------------------------------------------------------
c  OUTPUT:
c    - pmd00000
c    - akr0000
c-----------------------------------------------------------------------
      implicit real*8(a-h,o-z),integer(i-n)
      include '../pmd/params_unit.h'
      include './params_Ito_W-He.h'
c-----max # of atoms
      integer,parameter::nmax=10000000
c-----# of unit cells of bcc W crystall
      integer,parameter:: nuc(1:3)= (/ 5, 5, 3 /)
c.....# of unit cells of He gas region per both top and bottom
      integer,parameter:: nhe(1:3)= (/ 0, 0, 1 /)
c.....# of unit cells for vacuume region per both top and bottom
      integer,parameter:: nvac(1:3)= (/ 0, 0, 1 /)
      real(8):: ua(3,10),uhe(3,4)
      real(8):: tag(nmax),ra(3,nmax),va(3,nmax),eki(nmax),epi(nmax)
     &     ,h(3,3,0:1),strs(3,3,nmax)
      logical:: l1st_He=.true.
      real(8),external:: rnd

      small = 1d0/8 /(nuc(1)+nvac(1))
      write(6,'(a,f10.4)') ' small=',small

cc.....Lattice constant of W for 0K
c      cunit= 3.204d0
c.....Lattice constant of W for 1000K
      cunit= 3.226d0

c-----simulation box size
      h(1:3,1:3,0:1)= 0d0
      h(1,1,0)= cunit*(nuc(1)+nhe(1)+nvac(1))
      h(2,2,0)= cunit*(nuc(2)+nhe(2)+nvac(2))
      h(3,3,0)= cunit*(nuc(3)+nhe(3)+nvac(3))

c-----unit cell, BCC
      ua(1:3,1)= (/ 0.0d0, 0.0d0, 0.0d0 /)
      ua(1:3,2)= (/ 0.5d0, 0.5d0, 0.5d0 /)

c.....He atoms, assuming FCC
      uhe(1:3,1)= (/ 0.0d0, 0.0d0, 0.0d0 /)
      uhe(1:3,2)= (/ 0.5d0, 0.5d0, 0.0d0 /)
      uhe(1:3,3)= (/ 0.5d0, 0.0d0, 0.5d0 /)
      uhe(1:3,4)= (/ 0.0d0, 0.5d0, 0.5d0 /)
      
      inc=0 
c      do iz=nvac(3),nucs(3)-nvac(3)-1
      do iz=0,nuc(3)+nhe(3)-1
        do iy=0,nuc(2)-1
          do ix=0,nuc(1)-1
            if( iz.lt.nuc(3) ) then ! W atom
              do m=1,2
cc.....Vacancy?
c              if( iy.eq.nucs(2)/2 .and. ix.eq.nucs(1)/2
c     &             .and. iz.eq.nucs(3)/2 .and. m.eq.1 ) cycle 
                inc=inc+1
                if(inc.gt.nmax)then
                  write(*,*)'Error inc>nmax',inc,nmax
                  stop
                endif
                ra(1,inc)=(ua(1,m)+dble(ix))/(nuc(1)+nhe(1)+nvac(1))
     &               +small
                ra(2,inc)=(ua(2,m)+dble(iy))/(nuc(2)+nhe(2)+nvac(2))
     &               +small
                ra(3,inc)=(ua(3,m)+dble(iz))/(nuc(3)+nhe(3)+nvac(3))
     &               +small
                is= 1
                ifmv= 1
                if( iz.eq.0 ) ifmv= 0
                tag(inc)= 1d0*is +0.1d0*ifmv +1d-14*inc
              enddo
            else ! He atom
              do m=1,4
                inc=inc+1
                if(inc.gt.nmax)then
                  write(*,*)'Error inc>nmax',inc,nmax
                  stop
                endif
                ra(1,inc)=(uhe(1,m)+dble(ix))/(nuc(1)+nhe(1)+nvac(1))
     &               +small
                ra(2,inc)=(uhe(2,m)+dble(iy))/(nuc(2)+nhe(2)+nvac(2))
     &               +small
                ra(3,inc)=(uhe(3,m)+dble(iz))/(nuc(3)+nhe(3)+nvac(3))
     &               +small
                is= 2
c.....Only one W atom is on top of the surface
                if( l1st_He ) then
                  is= 1
                  l1st_He=.false.
                endif
                ifmv= 1
                tag(inc)= 1d0*is +0.1d0*ifmv +1d-14*inc
              enddo
            endif
          enddo
        enddo
      enddo

cc.....Add one atom on top of the surface
c      inc= inc +1
c      ra(1,inc)= rnd()
c      ra(2,inc)= rnd()
c      ra(3,inc)= (0.2d0+nuc(3))/(nuc(3)+nhe(3)+nvac(3)) +small
c      is= 1
c      ifmv= 1
c      tag(inc)= 1d0*is +0.1d0*ifmv +1d-14*inc

      write(6,'(a,i10)') " natm=",inc


      va(1:3,1:inc)= 0d0

      call write_pmd0_ascii(15,'pmd00000','replace',inc,tag,ra,va,h
     &     ,cunit,eki,epi,strs)
      
c-----output 'akr000' for Akira visualization
      open(15,file='akr0000',form='formatted',status='replace')
      write(15,'(es15.7)') cunit
      write(15,'(3es11.3)') ((h(ia,ib,0)/cunit,ia=1,3),ib=1,3)
      write(15,'(i10,3i5)') inc, 4, 0, 0
      do i=1,inc
        write(15,'(i3,7es11.3)') int(tag(i)),ra(1:3,i),va(1:3,i)
     &       ,dble(int((tag(i)-int(tag(i)))*10d0))
      enddo
      close(15)
      
      end program mkconf_WHe_compress
c=======================================================================
      function rnd()
      implicit none
      real(8):: rnd
      real(8),parameter:: d2p31m= 2147483647d0
      real(8),parameter:: d2p31 = 2147483648d0
      real(8),save:: dseed
      logical,save:: l1st=.true.

      if( l1st ) then
        dseed= 1234.0d0
        l1st=.false.
      endif
      
      dseed=dmod(16807d0*dseed,d2p31m)
      rnd=dseed/d2p31
      return
      end function rnd
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make 10mkconf"
c     End:
