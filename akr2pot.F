      program akr2pot
c-----------------------------------------------------------------------
c Reduce atoms from akr file according to their potential energy values.
c   - assuming akr file to be redirect to this program as STDIN
c   - output akr format to STDOUT
c-----------------------------------------------------------------------
      implicit none
c.....Threshold potential energy
      real(8),parameter:: pth = -0.1410d0

      integer:: i,n,itmp,nauxdat,j,ia,ib,l,natm
      real(8):: h(3,3)
      integer,allocatable:: is(:)
      real(8),allocatable:: ra(:,:),auxdat(:,:)

c-----read akr file from STDIN
      read(5,*) natm,nauxdat
      read(5,*) ((h(ia,ib),ia=1,3),ib=1,3)
      allocate(ra(3,natm),auxdat(nauxdat,natm),is(natm))
      n=0
      do i=1,natm
        read(5,*) is(i),ra(1:3,i),auxdat(1:nauxdat,i)
        if( auxdat(2,i).gt.pth ) n=n+1
      enddo

c-----write akr file to STDOUT
      write(6,'(i10,3i4)') n,nauxdat,0,0
      write(6,'(3es12.4)') ((h(ia,ib),ia=1,3),ib=1,3)
      do i=1,natm
        if( auxdat(2,i).gt.pth ) then
          write(6,'(i4,20es11.3)') is(i),ra(1:3,i)
     &         ,auxdat(1:nauxdat,i)
        endif
      enddo

      end program akr2pot
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make akr2pot"
c     End:
