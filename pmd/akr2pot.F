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
      real(8):: h(3,3),hunit
      integer,allocatable:: is(:)
      real(8),allocatable:: ra(:,:),auxdat(:,:)
      integer nargc
      integer,external:: iargc
      character*128:: ciname,coname

c.....Read an argument, ipmd
      nargc=iargc()
      if( nargc.ne.2 ) then
        stop 'Usage: $ ./akr2pot akr0000 pot0000'
      endif
      call getarg(1,ciname)
      call getarg(2,coname)


c-----read akr file from
      open(10,file=trim(ciname),status="old")
      read(10,*) hunit
      read(10,*) ((h(ia,ib)*hunit,ia=1,3),ib=1,3)
      read(10,*) natm,nauxdat
      allocate(ra(3,natm),auxdat(nauxdat,natm),is(natm))
      n=0
      do i=1,natm
        read(10,*) is(i),ra(1:3,i),auxdat(1:nauxdat,i)
        if( auxdat(2,i).gt.pth ) n=n+1
      enddo
      close(10)

c-----write akr file to STDOUT
      open(11,file=trim(coname),status="replace")
      write(11,'(es15.7)') hunit
      write(11,'(3es12.4)') ((h(ia,ib)/hunit,ia=1,3),ib=1,3)
      write(11,'(i10,3i4)') n,nauxdat,0,0
      do i=1,natm
        if( auxdat(2,i).gt.pth ) then
          write(11,'(i4,20es11.3)') is(i),ra(1:3,i)
     &         ,auxdat(1:nauxdat,i)
        endif
      enddo
      close(11)

      end program akr2pot
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make akr2pot"
c     End:
