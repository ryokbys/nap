      program poscar2akr
c-----------------------------------------------------------------------
c  Convert POSCAR to Akira format.
c-----------------------------------------------------------------------
c  USAGE:
c    $ poscar2akr < POSCAR > akr000
c-----------------------------------------------------------------------
      implicit none
      include 'params_unit.h'
      integer:: i,natm
      real(8):: al,h(3,3)
      character:: ctmp*128
      real(8),allocatable,dimension(:,:):: ra

      read(5,*) ctmp
      read(5,*) al
      read(5,*) h(1,1:3)
      read(5,*) h(2,1:3)
      read(5,*) h(3,1:3)
      h(1:3,1:3)= h(1:3,1:3) *al
      read(5,*) natm
      read(5,*) ctmp
      allocate(ra(3,natm))
      do i=1,natm
        read(5,*) ra(1:3,i)
      enddo

c.....Write akr000
      write(6,'(i10,3i5)') natm, 1, 0, 0
      write(6,'(3es15.7)') h(1,1:3)
      write(6,'(3es15.7)') h(2,1:3)
      write(6,'(3es15.7)') h(3,1:3)
      do i=1,natm
        write(6,'(i3,3es12.4,9es12.4)') 1,ra(1:3,i),0d0
      enddo

      end program poscar2akr
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make poscar2akr"
c     End:
