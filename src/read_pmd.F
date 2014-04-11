      program read_pmd
      use variables
      implicit none
      integer:: ia,ib,l,n,i,j

      integer,external:: iargc
      integer:: nargc
      character*128:: cargv
      real*8:: x

      nargc= iargc()
      if( nargc.lt.1 ) stop ' [Error] nargc.lt.1 !'
      call getarg(1,cargv)
      
      open(10,file=trim(cargv),form='unformatted',status='old')
      read(10) n
      read(10) (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
      read(10) (tag(i),ra(1:3,i),va(1:3,i)
     &       ,eki(1,1,i),epi(i),strs(1:3,1:3,i),i=1,n)
      close(10)
      write(6,'(a,i10)') " n=",n
      write(6,'(a)') " h="
      write(6,'(3es15.7)') (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
      do i=1,5
        write(6,'(10es24.14)') tag(i),ra(1:3,i),va(1:3,i)
      enddo

      end program read_pmd
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make rdpmd"
c     End:
