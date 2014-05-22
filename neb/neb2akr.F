      program neb2akr
c-----------------------------------------------------------------------
c  Convert neb??? file to akr???
c-----------------------------------------------------------------------
c  Usage:
c    $ ./neb2akr < akr??? > neb???
c-----------------------------------------------------------------------
      implicit none
      integer:: i,natm,ia,ib,l,ix,iy,iz,ncpx,ncpy,ncpz,natmt,m
      real(8):: h(3,3,0:1),ht(3,3,0:1),dcpx,dcpy,dcpz,hunit
      real(8),allocatable:: tag(:),ra(:,:),va(:,:),eki(:),epi(:)
     &     ,strs(:,:,:)
      real(8),allocatable:: rat(:,:),vat(:,:),tagt(:),ekit(:),epit(:)
     &     ,strst(:,:,:)
      real(8),parameter:: bohr =  0.5291772d-10
      real(8),parameter:: rc   =  4.3d-10 /bohr
      
      read(5,*) hunit
      read(5,*) (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
      h(1:3,1:3,0:1)= h(1:3,1:3,0:1)*hunit
      read(5,*) natm
      allocate(tag(natm),ra(3,natm),va(3,natm),eki(natm),epi(natm)
     &     ,strs(3,3,natm))
      do i=1,natm
        read(5,*) tag(i),ra(1:3,i),va(1:3,i),eki(i),epi(i)
     &       ,strs(1:3,1:3,i)
      enddo
      
      ncpx=1
      ncpy=1
      ncpz=1
      if( h(1,1,0).lt.2d0*rc ) then
        ncpx= 1+int(2d0*rc/h(1,1,0))
      elseif( h(2,2,0).lt.2d0*rc ) then
        ncpy= 1+int(2d0*rc/h(2,2,0))
      elseif( h(3,3,0).lt.2d0*rc ) then
        ncpz= 1+int(2d0*rc/h(3,3,0))
      endif
      natmt= natm*ncpx*ncpy*ncpz
c      write(6,'(4i10)') natmt,ncpx,ncpy,ncpz
      ht(1:3,1:3,0:1)= 0d0
      ht(1,1,0)=  h(1,1,0)*ncpx
      ht(2,2,0)=  h(2,2,0)*ncpy
      ht(3,3,0)=  h(3,3,0)*ncpz
      allocate(rat(3,natmt),vat(3,natmt),tagt(natmt)
     &     ,ekit(natmt),epit(natmt),strst(3,3,natmt))

      m=0
      do iz=1,ncpz
        dcpz= dble(iz-1)/ncpz
        do iy=1,ncpy
          dcpy= dble(iy-1)/ncpy
          do ix=1,ncpx
            dcpx= dble(ix-1)/ncpx
            do i=1,natm
              m=m+1
              tagt(m)= tag(i)
              rat(1,m)= ra(1,i)/ncpx +dcpx
              rat(2,m)= ra(2,i)/ncpy +dcpy
              rat(3,m)= ra(3,i)/ncpz +dcpz
              vat(1:3,m)= va(1:3,i)
              ekit(m)= eki(i)
              epit(m)= epi(i)
              strst(1:3,1:3,m)= strs(1:3,1:3,i)
            enddo
          enddo
        enddo
      enddo

      write(6,'(es15.7)') hunit
      write(6,'(3es15.7)') ht(1,1:3,0)/hunit
      write(6,'(3es15.7)') ht(2,1:3,0)/hunit
      write(6,'(3es15.7)') ht(3,1:3,0)/hunit
      write(6,'(i8,3i3)') natmt, 5, 0, 0
      do i=1,natmt
        write(6,'(i3,3es15.7,10es12.4)') int(tagt(i))
     &       ,rat(1:3,i),ekit(i),epit(i)
     &       ,strst(1,1,i),strst(2,2,i),strst(3,3,i)
      enddo

      end program neb2akr
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make neb2akr"
c     End:
