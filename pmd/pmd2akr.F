      program pmd2akr
c-----------------------------------------------------------------------
c Combine pmd files made by parallel program
c   - pmd#####-$$$$  -->  akr$$$$
c     #####: node number
c     $$$$: akira frame number
c USAGE:
c   $ ./pmd2akr #####
c-----------------------------------------------------------------------
      use variables
      implicit none
      integer,parameter:: npmax= 2100000
c.....Num of auxiliary data for coloring atoms in Akira
      integer,parameter:: nauxdat = 5
      integer:: ipmd,n,inode,itmp,i,nnodes,nkvs,ia,ib,l,itag(npmax)
      real(8):: tmp,dat(18,npmax)
      character:: fin*13,fout*7
      integer,external:: itot_tag
      integer nargc
c      integer,external:: iargc
      character*128:: argv

c.....Read an argument, ipmd
c      nargc=iargc()
      nargc= command_argument_count()
      if( nargc.gt.1 ) stop 'Usage: $ ./pmd2akr 10'
      call getarg(1,argv)
      read(argv,*) ipmd

      call read_input(10,'in.pmd')
      nnodes= nx*ny*nz
c      write(6,'(a,4i5)') " nx,ny,nz,nnodes=",nx,ny,nz,nnodes
c      write(6,'(a,i5)') " npmd =",npmd

      if( ipmd.lt.0 .or. ipmd.gt.npmd ) then
        write(6,*) ' Error: ipmd should be [0:',npmd,']'
        stop 
      endif
      
      n= 0
      do inode=0,nnodes-1
        fin="0000/pmd00000"
c-------about node number
        write(fin(9:13),'(i5.5)') inode
c---------about pmd-step number
        write(fin(1:4),'(i4.4)') ipmd
c---------read pmd###-$$$ file
        if( trim(ciofmt).eq.'bin'.or.trim(ciofmt).eq.'binary' ) then
          call read_pmd_bin(20,fin
     &         ,namax,itmp,h,hunit,tag,ra,va,eki,epi,strs)
        else if( trim(ciofmt).eq.'ascii' ) then
          call read_pmd_ascii(20,fin
     &         ,namax,itmp,h,hunit,tag,ra,va,eki,epi,strs)
        endif
        write(6,'(a,i8)') 'itmp=',itmp
        do i=1,itmp
          n=n+1
          dat(1:3,n)= ra(1:3,i)
          dat(4:6,n)= va(1:3,i)
          dat(7,n)= eki(1,1,i)
          dat(8,n)= epi(i)
          dat(9:11,n)= strs(1:3,1,i)
          dat(12:14,n)= strs(1:3,2,i)
          dat(15:17,n)= strs(1:3,3,i)
          dat(18,n)= tag(i)
          itag(n)= itot_tag(tag(i))
c            tagt(n)=tag(i)
c            rat(1:3,n)=ra(1:3,i)
c            vat(1:3,n)=va(1:3,i)
c            op(1,n)=eki(1,1,i)
c            op(2,n)=epi(i)
c            op(3,n)=strs(1,1,i)
c            op(4,n)=strs(2,2,i)
c            op(5,n)=strs(1,2,i)
        enddo
      enddo

cc.....Sort before output
c        call heapsort_itag(n,npmax,itag,18,dat)

c-------output file name
      fout="akr0000"
      write(fout(4:7),'(i4.4)') ipmd
c-------output combined mts-file
      open(30,file=fout,status="replace")
      write(30,'(es15.7)') hunit
      write(30,'(3es15.7)') ((h(ia,ib,0)/hunit,ia=1,3),ib=1,3)
      write(30,'(i10,3i4)') n, nauxdat, 0, 0
      write(30,'(i3,8es22.14)') (int(dat(18,i)),dat(1:3,i)
     &     ,dat(7:8,i),dat(9,i),dat(13,i),dat(17,i),i=1,n)
      close(30)
      write(6,'(a)') " "//fout//" done."

      end program pmd2akr
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
c=======================================================================
c-----for emacs---------------------------------------------------------
c     Local Variables:
c     compile-command: "make pmd2akr"
c     End:
