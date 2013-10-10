      program node_conv
c-----------------------------------------------------------------------
c   Node-layout converter program for pmd.
c     Read:  pmd?????-0000 (or some other name)
c     write: pmd?????-0000
c-----------------------------------------------------------------------
      implicit none
      integer,parameter:: namax = 2000000

      integer:: nox,noy,noz,nodeso,nnx,nny,nnz,nodesn,i,j,k,l,m,n
     &     ,nlast,ia,ib,myx,myy,myz
      real(8):: sxog,syog,szog,h(3,3,0:1)
      real(8):: tag(namax),tagt(namax)
      real(8):: ra(3,namax),va(3,namax),eki(namax),epi(namax)
     &     ,strs(3,3,namax)
      real(8):: rat(3,namax),vat(3,namax)
     &     ,ekit(namax),epit(namax),strst(3,3,namax)
      character:: filename*13,key*1

      write(6,'(a)') ' This program pmd?????-0000 files and '
     &     //'split to different node layout.'
      write(6,'(a)') ' Then override pmd?????-0000 files.'
      filename='pmd00000-0000'

c-----read some info from STDIN
      write(6,'(a)') ' Current node-layout: (nx,ny,nz)=?'
      read(5,*)nox,noy,noz
      nodeso=nox*noy*noz
      write(6,'(a)')' New node-layout: (nx,ny,nz)=?'
      read(5,*)nnx,nny,nnz
      nodesn=nnx*nny*nnz

      write(6,'(a,3i3,a,3i3,a)') 'Converts ( ',nox,noy,noz
     &     ,' ) to ( ',nnx,nny,nnz,' )'

      nlast=0
      do i=0,nodeso-1
        myx=i/(noy*noz)
        myy=mod(i/noz,noy)
        myz=mod(i,noz)
        sxog=dble(myx)/dble(nox)
        syog=dble(myy)/dble(noy)
        szog=dble(myz)/dble(noz)
c-------obtain filename
        write(filename(4:8),'(i5.5)') i
cc.....Binary format
c        open(12,file=trim(filename),form='unformatted',status='old')
c        read(12) n
c        if(nlast+n.gt.namax)then
c          write(6,'(a)')' [Error] n.gt.namax !!'
c          write(6,'(a,2i10)') ' nlast+n,namax =',nlast+n,namax
c          stop
c        endif
c        read(12) (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
c        read(12) (tag(j),ra(1:3,j),va(1:3,j)
c     &       ,eki(j),epi(j),strs(1:3,1:3,j),j=nlast+1,nlast+n)
c.....Ascii format
        open(12,file=trim(filename),form='formatted',status='old')
        read(12,*) n
        if(nlast+n.gt.namax)then
          write(6,'(a)')' [Error] n.gt.namax !!'
          write(6,'(a,2i10)') ' nlast+n,namax =',nlast+n,namax
          stop
        endif
        read(12,*) (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
        read(12,*) (tag(j),ra(1:3,j),va(1:3,j)
     &       ,eki(j),epi(j),strs(1:3,1:3,j),j=nlast+1,nlast+n)
c.....Some check
        do j=nlast+1,nlast+n
c          read(12,*) tag(j),ra(1:3,j),va(1:3,j)
          if(  ra(1,j).lt.sxog.or.ra(1,j).ge.sxog+1d0/nox.or.
     &         ra(2,j).lt.syog.or.ra(2,j).ge.syog+1d0/noy.or.
     &         ra(3,j).lt.szog.or.ra(3,j).ge.szog+1d0/noz)then
            write(*,*)'Old node-layout error!'
            write(6,'(i10,3es24.15)') j,ra(1:3,j)
            write(6,'(2es24.15)') sxog, sxog+1d0/nox
            write(6,'(2es24.15)') syog, syog+1d0/noy
            write(6,'(2es24.15)') szog, szog+1d0/noz
            stop
          endif
        enddo
        nlast=nlast+n
        close(12)
      enddo
      write(*,*)' Total number of atoms =',nlast


100   write(*,*)' Overwrite pmd###-000 files? (y/n)'
      read(*,*) key
      if(key.eq.'n')then
        stop
      elseif(key.ne.'y')then
        goto 100
      endif
      
      filename='pmd00000-0000'
      do i=0,nodesn-1
        myx=i/(nny*nnz)
        myy=mod(i/nnz,nny)
        myz=mod(i,nnz)
        sxog=dble(myx)/nnx
        syog=dble(myy)/nny
        szog=dble(myz)/nnz
        
        n=0
        do j=1,nlast
          if(ra(1,j).ge.sxog .and.
     &       ra(1,j).lt.sxog+1d0/(nnx) .and.
     &       ra(2,j).ge.syog .and.
     &       ra(2,j).lt.syog+1d0/(nny) .and.  
     &       ra(3,j).ge.szog .and.
     &       ra(3,j).lt.szog+1d0/(nnz)) then
            n=n+1
            tagt(n)=tag(j)
            rat(1:3,n)=ra(1:3,j)
            vat(1:3,n)=va(1:3,j)
            ekit(n)=eki(j)
            epit(n)=epi(j)
            strs(1:3,1:3,n)=strs(1:3,1:3,j)
          endif
        enddo
        write(6,'(a,i4,i8)') " node-id, num of atoms =",i,n

        write(filename(4:8),'(i5.5)') i
c        call write_pmd0_bin(12,filename,'replace',n,tagt,rat,vat,h
c     &       ,ekit,epit,strst)
        call write_pmd0_ascii(12,filename,'replace',n,tagt,rat,vat,h
     &       ,ekit,epit,strst)
      enddo

      write(6,'(a)') ' You must replace num_nodes_? values in in.pmd'
     &     //' as following...'
      write(6,'(a,i8)') '  num_nodes_x',nnx
      write(6,'(a,i8)') '  num_nodes_y',nny
      write(6,'(a,i8)') '  num_nodes_z',nnz

      end program node_conv
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make 20nconv"
c     End:
