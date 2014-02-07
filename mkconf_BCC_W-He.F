      program mkconf_BCC_WHe
c-----------------------------------------------------------------------
c  Create a BCC crystal with W and He atoms
c-----------------------------------------------------------------------
c  OUTPUT:
c    - pmd00000
c    - akr0000
c-----------------------------------------------------------------------
      implicit real*8(a-h,o-z),integer(i-n)
      include './params_unit.h'
      include './params_Ito_W-He.h'
c-----max # of atoms
      integer,parameter::nmax=100000
c-----# of unit cells
c      integer,parameter:: nuc(1:3)= (/ 1, 1, 1 /)
      integer,parameter:: nuc(1:3)= (/ 5,5,5 /)
c-----vacuum width in unit of cell
c      integer,parameter:: nvac(1:3)= (/ 5, 5, 5 /)
      integer,parameter:: nvac(1:3)= (/ 0, 0, 0 /)
      real(8):: ua(3,10)
      real(8):: tag(nmax),ra(3,nmax),va(3,nmax),eki(nmax),epi(nmax)
     &     ,h(3,3,0:1),strs(3,3,nmax),dseed,rnd
c.....vacancy ratio
      real(8),parameter:: vratio  = 0.0001
c.....He ratio
      real(8),parameter:: heratio = 0.050

      small = 1d0/4 /nuc(1)

cc.....Lattice constant of Fe, see Philos. Mag. 83 (2003) 3977
c      cunit= 2.835d0
c.....Lattice constant of W for 1000K
      cunit= 3.202d0

c-----simulation box size
      h(1:3,1:3,0:1)= 0d0
      h(1,1,0)= cunit*(nuc(1)+nvac(1))
      h(2,2,0)= cunit*(nuc(2)+nvac(2))
      h(3,3,0)= cunit*(nuc(3)+nvac(3))

c-----unit cell, BCC
      ua(1:3,1)= (/ 0.0d0, 0.0d0, 0.0d0 /)
      ua(1:3,2)= (/ 0.5d0, 0.5d0, 0.5d0 /)
      
      dseed= 123456d0
      inc=0 
      do ix=0,nuc(1)-1
        do iy=0,nuc(2)-1
          do iz=0,nuc(3)-1
            do m=1,2
c              call myrnd(rnd,dseed)
c              if( rnd.lt.vratio ) cycle
              if( ix.eq.nuc(1)/2 .and. iy.eq.nuc(2)/2
     &             .and. iz.eq.nuc(3)/2 .and. m.eq.2 ) cycle
              inc=inc+1
              if(inc.gt.nmax)then
                write(*,*)'Error inc>nmax',inc,nmax
                stop
              endif
              ra(1,inc)= (ua(1,m)+dble(ix))/(nuc(1)+nvac(1)) +small
              ra(2,inc)= (ua(2,m)+dble(iy))/(nuc(2)+nvac(2)) +small
              ra(3,inc)= (ua(3,m)+dble(iz))/(nuc(3)+nvac(3)) +small
c              ra(1:3,inc)= ra(1:3,inc)
c     &             +dble(nvac(1:3))/(nuc(1:3)+nvac(1:3))/2
              is= 1
c              if( m.eq.2 ) is= 2
              ifmv= 1
              tag(inc)= 1d0*is +0.1d0*ifmv +1d-14*inc
            enddo
            call myrnd(rnd,dseed)
            if( rnd.lt.heratio ) then
c.....Add He atom
              inc= inc+1
              if(inc.gt.nmax)then
                write(*,*)'Error inc>nmax',inc,nmax
                stop
              endif
              is= 2
              ifmv= 1
              tag(inc)= 1d0*is +0.1d0*ifmv +1d-14*inc
cc.....O-site
c              ra(1,inc)= (0.5d0 +dble(ix))/(nuc(1)+nvac(1)) +small
c              ra(2,inc)= (0.5d0 +dble(iy))/(nuc(2)+nvac(2)) +small
c              ra(3,inc)= (0.0d0 +dble(iz))/(nuc(3)+nvac(3)) +small
c.....T-site
              ra(1,inc)= (0.50d0 +dble(ix))/(nuc(1)+nvac(1)) +small
              ra(2,inc)= (0.25d0 +dble(iy))/(nuc(2)+nvac(2)) +small
              ra(3,inc)= (0.00d0 +dble(iz))/(nuc(3)+nvac(3)) +small
            endif
          enddo
        enddo
      enddo

      write(6,'(a,i10)') " natm=",inc

      call write_pmd0_ascii(15,'pmd00000','replace',inc,tag,ra,va,h
     &     ,eki,epi,strs)
      
c-----output 'akr000' for Akira visualization
      open(15,file='akr0000',form='formatted',status='replace')
      write(15,'(i10,3i5)') inc, 3, 0, 0
      write(15,'(3es11.3)') ((h(ia,ib,0),ia=1,3),ib=1,3)
      do i=1,inc
        write(15,'(i3,6es11.3)') int(tag(i)),ra(1:3,i),va(1:3,i)
      enddo
      close(15)
      
      end program mkconf_BCC_WHe
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
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make 10mkconf"
c     End:
