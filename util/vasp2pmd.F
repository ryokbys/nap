      program vasp2pmd
      implicit none
      integer:: i,j,n,natm,ia,ib,l,mis(9)
      logical:: lexist
      character:: ctmp*128
      real(8):: h(3,3,0:1),al
      real(8),allocatable:: tag(:),ra(:,:),va(:,:),eki(:),epi(:)
     &     ,strs(:,:,:)

      integer,parameter:: nis=  2
      real(8),parameter:: aa2bohr = 1d0/0.5291772d0
      real(8),parameter:: bohr2aa = 1d0/aa2bohr

      read(5,*) ctmp
      read(5,*) al
      read(5,*) h(1,1:3,0)
      read(5,*) h(2,1:3,0)
      read(5,*) h(3,1:3,0)
      read(5,*) (mis(n),n=1,nis)
      read(5,*) ctmp
      natm= 0
      do n=1,nis
        natm=natm +mis(n)
      enddo
      allocate(tag(natm),ra(3,natm),va(3,natm),eki(natm),epi(natm)
     &     ,strs(3,3,natm))
      natm= 0
      do n=1,nis
        do i=1,mis(n)
          natm=natm+1
          read(5,*) ra(1:3,natm)
          tag(natm)= n
        enddo
      enddo
      
      h(1:3,1:3,0)= h(1:3,1:3,0)*al
      write(6,'(i6)') natm
      write(6,'(3es22.14)') (((h(ia,ib,l)*aa2bohr,ia=1,3),ib=1,3),l=0,1)
      do i=1,natm
        write(6,'(4es22.14,20es11.3)') tag(i),ra(1:3,i),va(1:3,i)
     &       ,eki(i),epi(i),strs(1:3,1:3,i)
      enddo

      end program vasp2pmd
