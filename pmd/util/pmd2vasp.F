      program pmd2vasp
      implicit none
      integer:: i,natm,lis(9),n,nis,is,ia,ib,l,mis(9)
      logical:: lexist
      real(8):: h(3,3,0:1),hunit
      real(8),allocatable:: tag(:),ra(:,:),va(:,:),eki(:),epi(:)
     &     ,strs(:,:,:)

      read(5,*) hunit
      read(5,*) (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
      h(1:3,1:3,0:1)= h(1:3,1:3,0:1)*hunit
      read(5,*) natm
      allocate(tag(natm),ra(3,natm),va(3,natm),eki(natm),epi(natm)
     &     ,strs(3,3,natm))
      read(5,*) (tag(i),ra(1:3,i),va(1:3,i)
     &     ,eki(i),epi(i),strs(1:3,1:3,i),i=1,natm)

c.....Count num of species
      nis= 0
      mis(1:9)= 0
      do i=1,natm
        lexist= .false.
        is= int(tag(i))
        do n=1,nis
          if( is.eq.lis(n) ) then
            lexist= .true.
            exit
          endif
        enddo
        if(.not.lexist) then
          nis=nis+1
          if( nis.gt.9 ) stop ' [Error] nis.gt.9 !!!'
          lis(nis)= is
          mis(nis)=mis(nis) +1
        else
          mis(n)= mis(n) +1
        endif
      enddo

      write(6,'(a)') ' converted from pmd'
      write(6,'(es11.3)') hunit
      write(6,'(3es22.14)') h(1,1:3,0)/hunit
      write(6,'(3es22.14)') h(2,1:3,0)/hunit
      write(6,'(3es22.14)') h(3,1:3,0)/hunit
      write(6,'(9i4)') (mis(n),n=1,nis)
      write(6,'(a)') ' Direct'
      do n=1,nis
        do i=1,natm
          is= int(tag(i))
          if( is.ne.lis(n) ) cycle
          write(6,'(3es22.14)') ra(1:3,i)
        enddo
      enddo
      
      end program pmd2vasp
