      program mkconf_WHe_compress
c-----------------------------------------------------------------------
c  Create a BCC tungsten(W) crystal, He gas and vacuum region.
c-----------------------------------------------------------------------
c  OUTPUT:
c    - pmd00000
c    - akr0000
c-----------------------------------------------------------------------
      implicit real*8(a-h,o-z),integer(i-n)
      include '../src/params_unit.h'
      include './params_Ito_W-He.h'
c-----max # of atoms
      integer,parameter::nmax=10000000
c-----# of unit cells of bcc W crystall
      integer,parameter:: nuc(1:3)= (/ 5, 5, 3 /)
c.....# of unit cells of He gas region per both top and bottom
      integer,parameter:: nhe(1:3)= (/ 0, 0, 0 /)
c.....# of unit cells for vacuume region per both top and bottom
      integer,parameter:: nvac(1:3)= (/ 0, 0, 2 /)
c.....# of unit cells
      integer,parameter:: nucs(1:3)= nuc(1:3) +2*nhe(1:3) +2*nvac(1:3)
      real(8):: ua(3,10)
      real(8):: tag(nmax),ra(3,nmax),va(3,nmax),eki(nmax),epi(nmax)
     &     ,h(3,3,0:1),strs(3,3,nmax)
c.....Parameters
      real(8),parameter:: tempini = 100d0

      small = 1d0/8 /nucs(1)
      write(6,'(a,f10.4)') ' small=',small

cc.....Lattice constant of W for 0K
c      cunit= 3.204d0
c.....Lattice constant of W for 1000K
      cunit= 3.226d0

c-----simulation box size
      h(1:3,1:3,0:1)= 0d0
      h(1,1,0)= cunit*nucs(1)
      h(2,2,0)= cunit*nucs(2)
      h(3,3,0)= cunit*nucs(3)

c-----unit cell, BCC
      ua(1:3,1)= (/ 0.0d0, 0.0d0, 0.0d0 /)
      ua(1:3,2)= (/ 0.5d0, 0.5d0, 0.5d0 /)
      
      inc=0 
c      do iz=nvac(3),nucs(3)-nvac(3)-1
      do iz=0,nucs(3)-1
        do iy=0,nuc(2)-1
          do ix=0,nuc(1)-1
            do m=1,2
cc.....Vacancy?
c              if( iy.eq.nucs(2)/2 .and. ix.eq.nucs(1)/2
c     &             .and. iz.eq.nucs(3)/2 .and. m.eq.1 ) cycle 
              inc=inc+1
              if(inc.gt.nmax)then
                write(*,*)'Error inc>nmax',inc,nmax
                stop
              endif
              ra(1,inc)= (ua(1,m)+dble(ix))/nucs(1) +small
              ra(2,inc)= (ua(2,m)+dble(iy))/nucs(2) +small
              ra(3,inc)= (ua(3,m)+dble(iz))/nucs(3) +small
c              ra(1:3,inc)= ra(1:3,inc)
c     &             +dble(nvac(1:3))/(nuc(1:3)+nvac(1:3))/2
              is= 1
c              if( iz.le.nvac(3)+nhe(3)-1 ) is= 2
c              if( iz.ge.nvac(3)+nhe(3)+nuc(3) ) is= 2
              ifmv= 1
c              if( iz.eq.nvac(3)+nhe(3) ) ifmv=2
c              if( iz.eq.nvac(3)+nhe(3)+nuc(3)-1 ) ifmv=2
              
              tag(inc)= 1d0*is +0.1d0*ifmv +1d-14*inc
            enddo
          enddo
        enddo
      enddo
      write(6,'(a,i10)') " natm=",inc

      write(6,'(a,es15.7)') ' wbot = ',dble(nvac(3))/nucs(3)
      write(6,'(a,es15.7)') ' wtop = ',1d0-dble(nvac(3))/nucs(3)

      va(1:3,1:inc)= 0d0

      call write_pmd0_ascii(15,'pmd00000','replace',inc,tag,ra,va,h
     &     ,eki,epi,strs)
      
c-----output 'akr000' for Akira visualization
      open(15,file='akr0000',form='formatted',status='replace')
      write(15,'(i10,3i5)') inc, 4, 0, 0
      write(15,'(3es11.3)') ((h(ia,ib,0),ia=1,3),ib=1,3)
      do i=1,inc
        write(15,'(i3,7es11.3)') int(tag(i)),ra(1:3,i),va(1:3,i)
     &       ,dble(int((tag(i)-int(tag(i)))*10d0))
      enddo
      close(15)
      
      end program mkconf_WHe_compress
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make 10mkconf"
c     End:
