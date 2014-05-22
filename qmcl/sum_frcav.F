      program sum_frcav
      use variables
      implicit none
      integer,external:: iargc
      integer:: nargc,i
      character*128:: cfname
      real(8):: fsum

      nargc= iargc()
      if( nargc.ne.1 ) stop ' [main] nargc.ne.1 !!!'
      call getarg(1,cfname)

      call read_frcav(10,trim(cfname))

      fsum=0d0
      do i=1,natm
        fsum=fsum +facc(i)
      enddo

      write(6,'(a,es22.14)') ' fsum,fsum/istp0=',fsum,fsum/istp0

      end program sum_frcav
c=======================================================================
      subroutine read_frcav(ionum,cfname)
      use variables
      implicit none
      integer,intent(in):: ionum
      character(len=*),intent(in):: cfname

      integer:: itmp,i
      real(8):: tmp

      facc(1:natm)= 0d0
      open(ionum,file=trim(cfname),status='old',err=100)
      read(ionum,*) natm,istp0
      do i=1,natm
        read(ionum,*) facc(i),tmp
      enddo
      close(ionum)
      return

 100  write(6,'(a)') ' [read_frcav] '//trim(cfname)//' is not read !'
      return
      end subroutine read_frcav
c=======================================================================
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make sfrcav"
c     End:
