      program sort_qmcl
c-----------------------------------------------------------------------
c  This program sorts qmcl??? file ascending order of itot
c  in TAG variable.
c-----------------------------------------------------------------------
c  USAGE:
c    $ ./sort qmcl000
c-----------------------------------------------------------------------
      use variables
      implicit none
      integer:: nargc,i
      integer,external:: iargv
      character:: cin*128
      integer,allocatable,dimension(:):: itag
      real(8),allocatable:: opt(:,:)


      nargc= iargc()
      if( nargc.ne.1 ) then
        write(6,'(a)') ' [Error] Num of arguments is wrong !!!'
        write(6,'(a)') '   Usage: $ ./sort qmcl000'
        stop
      endif
      call getarg(1,cin)

      call read_pmd_ascii(10,trim(cin)
     &     ,namax,natm,h,tag,ra,va,eki,epi,strs)

      allocate(opt(18,namax),itag(namax))

      do i=1,natm
        opt(1:3,i)= ra(1:3,i)
        opt(4:6,i)= va(1:3,i)
        opt(7,i)= eki(i)
        opt(8,i)= epi(i)
        opt(9:11,i) = strs(1:3,1,i)
        opt(12:14,i)= strs(1:3,2,i)
        opt(15:17,i)= strs(1:3,3,i)
        opt(18,i)= tag(i)
        itag(i)= itot_tag(tag(i))
      enddo

      call heapsort_itag(natm,namax,itag,18,opt)

      do i=1,natm
        ra(1:3,i)= opt(1:3,i)
        va(1:3,i)= opt(4:6,i)
        eki(i)= opt(7,i)
        epi(i)= opt(8,i)
        strs(1:3,1,i)= opt(9:11,i)
        strs(1:3,2,i)= opt(12:14,i)
        strs(1:3,3,i)= opt(15:17,i)
        tag(i)= opt(18,i)
      enddo

      call write_pmd_ascii(11,trim(cin)
     &     ,natm,h,tag,ra,va,eki,epi,strs)

      end program sort_qmcl
c=======================================================================
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make sort"
c     End:
