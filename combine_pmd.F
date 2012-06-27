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
c.....Num of auxiliary data for coloring atoms in Akira
      integer,parameter:: nauxdat = 5
      integer:: ipmd,n,ncg,inode,itmp,i,nnodes,nkvs,ia,ib,l
      integer:: ifmv(namax)
      real(8):: tmp,tagt(namax),rat(3,namax),vat(3,namax),op(9,namax)
      character:: fin*10,fout*6

      write(6,'(a)') " <<< Utility program combine_pmd >>>"

      call read_input(10,'in.pmd')
      nnodes= nx*ny*nz
      write(6,'(a,4i5)') " nx,ny,nz,nnodes=",nx,ny,nz,nnodes
      write(6,'(a,i5)') " npmd =",npmd
      
      do ipmd=0,npmd
        n= 0
        ncg= 0
        do inode=0,nnodes-1
          fin="pmd000-000"
c---------about node number
          write(fin(4:6),'(i3.3)') inode
c---------about pmd-step number
          write(fin(8:10),'(i3.3)') ipmd
c---------read pmd###-$$$ file
          open(20,file=fin,form='unformatted',status="old")
          read(20) itmp
          ncg= ncg +itmp
          read(20) (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
          read(20) (tag(i),ra(1:3,i),va(1:3,i)
     &         ,eki(i),epi(i),strs(1:3,1:3,i),i=1,itmp)
          do i=1,itmp
            n=n+1
            tagt(n)=tag(i)
            rat(1:3,n)=ra(1:3,i)
            vat(1:3,n)=ra(1:3,i)
            op(1,n)=eki(i)
            op(2,n)=epi(i)
            op(3,n)=strs(1,1,i)
            op(4,n)=strs(2,2,i)
            op(5,n)=strs(1,2,i)
          enddo
          close(20)
        enddo
c-------check n & ncg
        if(n.ne.ncg) stop "n.ne.ncg!"
c-------output file name
        fout="akr000"
        write(fout(4:6),'(i3.3)') ipmd
c-------output combined mts-file
        open(30,file=fout,status="replace")
        write(30,'(i10,3i4)') ncg, nauxdat, 0, 0
        write(30,'(3es15.7)') ((h(ia,ib,0),ia=1,3),ib=1,3)
        write(30,'(i3,8es12.4)') (int(tagt(i)),rat(1:3,i)
     &       ,op(1:nauxdat,i),i=1,ncg)
        close(30)
        write(6,'(a)') " "//fout//" done."
      enddo

      end program combine_pmd
c-----for emacs---------------------------------------------------------
c     Local Variables:
c     compile-command: "make 40combine"
c     End:
