      program modify_pmd_boxsize
      use variables
      implicit none
      integer:: ia,ib,l,n,i,j

      integer,external:: iargc
      integer:: nargc
      character(len=128):: cin,cout,ctmp
      real(8):: x,amat(3,3),d,tmat(3,3)

      nargc= iargc()
      if( nargc.ne.3 ) then
        write(6,'(a)') ' [Error] nargc.lt.1 !'
        write(6,'(a)') ' $ ./boxsize infile 1.01 outfile'
      endif
      call getarg(1,cin)
      call getarg(2,ctmp)
      read(ctmp,*) x
      call getarg(3,cout)
      
c.....Read
c      open(10,file=trim(cin),form='unformatted',status='old')
c      read(10) n
c      read(10) (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
c      read(10) (tag(i),ra(1:3,i),va(1:3,i)
c     &       ,eki(1,1,i),epi(i),strs(1:3,1:3,i),i=1,n)
c      close(10)
      if( trim(ciofmt).eq.'bin' .or. trim(ciofmt).eq.'binary' ) then
        call read_pmd_bin(20,trim(cin)
     &       ,namax,natm,h,tag,ra,va,eki,epi,strs)
      else if( trim(ciofmt).eq.'ascii' ) then
        call read_pmd_ascii(20,trim(cin)
     &     ,namax,natm,h,tag,ra,va,eki,epi,strs)
      endif
      
      write(6,'(a)') ' h-mat before multiplication:'
      write(6,'(3es12.4)') h(1,1:3,0)
      write(6,'(3es12.4)') h(2,1:3,0)
      write(6,'(3es12.4)') h(3,1:3,0)
c.....Uniform deformation
      h(1:3,1:3,0)= h(1:3,1:3,0) *x
c.....Uniaxial deformation
c      h(1,1,0)=h(1,1,0) *x
c      h(2,2,0)=h(2,2,0) *x
cc.....for elastic stiffness strain
c      d= 1d0-x
c      amat(1:3,1:3)= 0d0
c      amat(1:3,1)=(/ 1d0, 0.5d0*d, 0d0 /)
c      amat(1:3,2)=(/ 0.5d0*d, 1d0, 0d0 /)
c      amat(1:3,3)=(/ 0d0, 0d0, 1d0+d**2/(4d0-d**2) /)
c      tmat(1:3,1:3)= matmul(amat(1:3,1:3),h(1:3,1:3,0))
c      h(1:3,1:3,0)= tmat(1:3,1:3)
cc.....for shear strain
c      d= 1d0-x
c      amat(1:3,1:3)= 0d0
c      amat(1:3,1)=(/ 1d0+d, 0d0, 0d0 /)
c      amat(1:3,2)=(/ 0d0, 1d0-d, 0d0 /)
c      amat(1:3,3)=(/ 0d0, 0d0, 1d0+d**2/(1d0-d**2) /)
c      tmat(1:3,1:3)= matmul(amat(1:3,1:3),h(1:3,1:3,0))
c      h(1:3,1:3,0)= tmat(1:3,1:3)
      write(6,'(a,es12.4)') ' x =',x
      write(6,'(a)') ' h-mat after multiplication:'
      write(6,'(3es12.4)') h(1,1:3,0)
      write(6,'(3es12.4)') h(2,1:3,0)
      write(6,'(3es12.4)') h(3,1:3,0)

c.....Write
      if( trim(ciofmt).eq.'bin' .or. trim(ciofmt).eq.'binary' ) then
        call write_pmd_bin(20,trim(cout)
     &       ,natm,h,tag,ra,va,eki,epi,strs,sorg,dt)
      elseif( trim(ciofmt).eq.'ascii' ) then
        call write_pmd_ascii(20,trim(cout)
     &       ,natm,h,tag,ra,va,eki,epi,strs,sorg,dt)
      endif
c      open(11,file=trim(cout),form='unformatted',status='replace')
c      write(11) n
c      write(11) (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
c      write(11) (tag(i),ra(1:3,i),va(1:3,i)
c     &       ,eki(1,1,i),epi(i),strs(1:3,1:3,i),i=1,n)
c      close(11)

      end program modify_pmd_boxsize
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make boxsize"
c     End:
