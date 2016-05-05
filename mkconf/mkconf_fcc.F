      program mkconf_fcc
c
c    Create a FCC crystal
c    Write: ini000
c
      implicit real*8(a-h,o-z)
      include "../pmd/params_unit.h"
c-----max # of atoms
      integer,parameter::nmax=1000000
c-----# of unit cells
      integer,parameter:: nuc(1:3)= (/ 2,2,2 /)
c-----vacuum width in unit of cell
      integer,parameter:: nvac(1:3)= (/ 0, 0, 0 /)
      real(8):: ua(3,10)
      real(8):: tag(nmax),ra(3,nmax),va(3,nmax),eki(nmax),epi(nmax)
     &      ,h(3,3,0:1),strs(3,3,nmax),hunit
      real*8,parameter:: pi=3.141592653589793d0
c-----for Ar
      real*8,parameter:: sgmlj = 3.41d0

      small=1d-7

c-----box size for Ar
      hunit= 2d0**(1d0/6) *sgmlj *dsqrt(2d0)
cc-----box size for Al
c      cunit= 4.025d-10 /aa

c-----unit cell, FCC
      ua(1:3,1)= (/ 0.0d0, 0.0d0, 0.0d0 /)
      ua(1:3,2)= (/ 0.5d0, 0.5d0, 0.0d0 /)
      ua(1:3,3)= (/ 0.0d0, 0.5d0, 0.5d0 /)
      ua(1:3,4)= (/ 0.5d0, 0.0d0, 0.5d0 /)
      
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
              is= 1
              ifmv= 1
              tag(inc)= 1d0*is +0.1d0*ifmv +1d-14*inc
            enddo
          enddo
        enddo
      enddo
      write(6,'(a,i10)') " natm=",inc
      write(6,'(a,i10)') " id of inc=",int(mod(tag(inc)*1d14,1d13))
c-----simulation box size
      h(1:3,1:3,0:1)= 0d0
      h(1,1,0)= hunit*(nuc(1)+nvac(1))
      h(2,2,0)= hunit*(nuc(2)+nvac(2))
      h(3,3,0)= hunit*(nuc(3)+nvac(3))
      
      va(1:3,1:inc)= 0d0

      call write_pmd0_ascii(15,'pmd0000','replace',inc,tag
     &     ,ra,va,h,hunit,eki,epi,strs)
      
c-----output 'akr0000' for KVS visualization
      open(15,file='akr0000',form='formatted',status='replace')
      write(15,'(es15.7)') hunit
      write(15,'(3es11.3)') ((h(ia,ib,l)/hunit,ia=1,3),ib=1,3)
      write(15,'(i10,3i5)')inc,3,0,0
      do i=1,inc
        write(15,'(i3,6es11.3)') int(tag(i)),ra(1:3,i),va(1:3,i)
      enddo
      close(15)
      
      end program mkconf_fcc
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make mkconf_fcc"
c     End:
