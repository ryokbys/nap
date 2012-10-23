      program combine_pmd
c-----------------------------------------------------------------------
c Combine pmd files made by parallel program
c   - pmd###-$$$  -->  akr$$$
c     ###: node number
c     $$$: akira frame number
c   - only combine them, not shifting positions
c-----------------------------------------------------------------------
      use variables
      implicit none
      integer,parameter:: npmax= 2000000
c.....Num of auxiliary data for coloring atoms in Akira
      integer,parameter:: nauxdat = 5
      integer:: ipmd,n,inode,itmp,i,nnodes,nkvs,ia,ib,l
      real(8):: tmp,tagt(npmax),rat(3,npmax),vat(3,npmax),op(9,npmax)
      character:: fin*10,fout*6

      write(6,'(a)') " <<< Utility program combine_pmd >>>"

      call read_input(10,'in.pmd')
      nnodes= nx*ny*nz
      write(6,'(a,4i5)') " nx,ny,nz,nnodes=",nx,ny,nz,nnodes
      write(6,'(a,i5)') " npmd =",npmd
      
      do ipmd=0,npmd
        n= 0
        do inode=0,nnodes-1
          fin="pmd000-000"
c---------about node number
          write(fin(4:6),'(i3.3)') inode
c---------about pmd-step number
          write(fin(8:10),'(i3.3)') ipmd
c---------read pmd###-$$$ file
          if( trim(ciofmt).eq.'bin'.or.trim(ciofmt).eq.'binary' ) then
            call read_pmd_bin(20,fin
     &           ,namax,itmp,h,tag,ra,va,eki,epi,strs)
          else if( trim(ciofmt).eq.'ascii' ) then
            call read_pmd_ascii(20,fin
     &           ,namax,itmp,h,tag,ra,va,eki,epi,strs)
          endif
          write(6,'(a,i8)') 'itmp=',itmp
          do i=1,itmp
            n=n+1
            tagt(n)=tag(i)
            rat(1:3,n)=ra(1:3,i)
            vat(1:3,n)=ra(1:3,i)
            op(1,n)=eki(1,1,i)
            op(2,n)=epi(i)
            op(3,n)=strs(1,1,i)
            op(4,n)=strs(2,2,i)
            op(5,n)=strs(1,2,i)
          enddo
        enddo
c-------output file name
        fout="akr000"
        write(fout(4:6),'(i3.3)') ipmd
c-------output combined mts-file
        open(30,file=fout,status="replace")
        write(30,'(i10,3i4)') n, nauxdat, 0, 0
        write(30,'(3es15.7)') ((h(ia,ib,0),ia=1,3),ib=1,3)
        write(30,'(i3,8es12.4)') (int(tagt(i)),rat(1:3,i)
     &       ,op(1:nauxdat,i),i=1,n)
        close(30)
        write(6,'(a)') " "//fout//" done."
      enddo

      end program combine_pmd
c-----for emacs---------------------------------------------------------
c     Local Variables:
c     compile-command: "make 40combine"
c     End:
