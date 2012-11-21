      program sort_pmd
c-----------------------------------------------------------------------
c  This program sorts pmd??? file ascending order of itot
c  in TAG variable.
c-----------------------------------------------------------------------
c  USAGE:
c    $ ./sort pmd000
c-----------------------------------------------------------------------
      use variables
      implicit none
      include 'params_au.h'
      integer:: nargc,i
      integer,external:: iargv,itot_tag
      character:: cin*128
      integer,allocatable,dimension(:):: itag
      real(8),allocatable:: opt(:,:)

      sorg(1:3)= 0d0
      dt= 1d-15/aut

      nargc= iargc()
      if( nargc.ne.1 ) then
        write(6,'(a)') ' [Error] Num of arguments is wrong !!!'
        write(6,'(a)') '   Usage: $ ./sort pmd000'
        stop
      endif
      call getarg(1,cin)

      call read_pmd_ascii(10,trim(cin)
     &     ,namax,natm,h,tag,ra,va,eki,epi,strs)

      allocate(opt(18,namax),itag(namax))

      do i=1,natm
        opt(1:3,i)= ra(1:3,i)
        opt(4:6,i)= va(1:3,i)
        opt(7,i)= eki(1,1,i)
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
        eki(1,1,i)= opt(7,i)
        epi(i)= opt(8,i)
        strs(1:3,1,i)= opt(9:11,i)
        strs(1:3,2,i)= opt(12:14,i)
        strs(1:3,3,i)= opt(15:17,i)
        tag(i)= opt(18,i)
      enddo

      call write_pmd_ascii(11,trim(cin)
     &     ,natm,h,tag,ra,va,eki,epi,strs,sorg,dt)

      end program sort_pmd
c=======================================================================
      function is_tag(tag)
      implicit none
      real(8),intent(in):: tag
      integer:: is_tag
      
      is_tag= int(tag)
      return
      end function is_tag
c=======================================================================
      function ifmv_tag(tag)
      implicit none
      real(8),intent(in):: tag
      integer:: ifmv_tag
      
      ifmv_tag= int((tag -int(tag))*10)

      return
      end function ifmv_tag
c=======================================================================
      function itot_tag(tag)
      implicit none
      real(8),intent(in):: tag
      integer:: itot_tag
      real(8):: tmp
      integer,external:: is_tag,ifmv_tag

      tmp= tag -is_tag(tag) -ifmv_tag(tag)*1d-1
      itot_tag= nint(tmp*1d+14)
      return
      end function itot_tag
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make sortpmd"
c     End:
