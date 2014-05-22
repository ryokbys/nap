      program mkconf_Al_fcc
c-----------------------------------------------------------------------
c  Make atom configuration file of fcc Al
c    Output: qmcl000
c-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      include "./params_unit.h"
c-----max # of atoms
      integer,parameter::nmax=1000000
c-----# of unit cells
      integer,parameter:: nuc(1:3)= (/ 20,20,1 /)
c-----vacuum width in unit of cell
      integer,parameter:: nvac(1:3)= (/ 0, 0, 0 /)
      real(8):: x,y,z,r
      real(8):: ua(3,4)
      real(8):: ra(3,nmax)
      real(8):: h(3,3,0:1)
      real(8):: va(3,nmax)
      real(8):: tag(nmax),eki(nmax),epi(nmax),strs(3,3,nmax)
c.....Radius
      real(8),parameter:: rc     = 10d0

      small= 0.1d0 /nuc(1)

c-----Lattice constant of this potential= 4.0425 AA
      al= 4.0425d0

c-----unit cell, wurtzite structure
c-----Al
      ua(1:3,1)= (/ 0.00d0, 0.00d0, 0.00d0 /)
      ua(1:3,2)= (/ 0.50d0, 0.50d0, 0.00d0 /)
      ua(1:3,3)= (/ 0.50d0, 0.00d0, 0.50d0 /)
      ua(1:3,4)= (/ 0.00d0, 0.50d0, 0.50d0 /)
      
c-----simulation box size
      h(1:3,1:3,0:1)= 0d0
      h(1,1,0)= (nuc(1)+nvac(1)) *al
      h(2,2,0)= (nuc(2)+nvac(2)) *al
      h(3,3,0)= (nuc(3)+nvac(3)) *al

      inc=0 
      do ix=0,nuc(1)-1
        do iy=0,nuc(2)-1
          do iz=0,nuc(3)-1
            do m=1,4
c              x= (ua(1,m)+dble(ix))/(nuc(1)+nvac(1)) +small
c              y= (ua(2,m)+dble(iy))/(nuc(2)+nvac(2)) +small
c              z= (ua(3,m)+dble(iz))/(nuc(3)+nvac(3)) +small
c              x= (x-0.5d0) *h(1,1,0)
c              y= (y-0.5d0) *h(2,2,0)
c              z= (z-0.5d0) *h(3,3,0)
c              r= sqrt(x**2 +y**2)
c              if( r.gt.rc ) cycle
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
      
      va(1:3,1:inc)= 0d0

      call write_pmd_ascii(15,'qmcl000',inc,h,tag,ra,va
     &     ,eki,epi,strs)
      call write_pmd_ascii(15,'qmclref',inc,h,tag,ra,va
     &     ,eki,epi,strs)
      
c-----output 'akr000' for KVS visualization
      open(15,file='akr000',form='formatted',status='replace')
      write(15,'(i10,3i4)') inc, 3, 0, 0
      write(15,'(3es11.3)') ((h(ia,ib,0),ia=1,3),ib=1,3)
      do i=1,inc
        write(15,'(i3,6es11.3)') int(tag(i)),ra(1:3,i),va(1:3,i)
      enddo
      close(15)
      
      end program mkconf_Al_fcc
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
