      program get_core_energy
      implicit none
      character(len=6),parameter:: cref= 'akr000'
      character(len=6),parameter:: cfin= 'akr010'
      real(8),parameter:: bohr= 0.5291772d-10
      real(8),parameter:: rcc = 4.9555d-10 /bohr
      integer:: i,j,na,n,itmp
      real(8):: epot,tmp,h(3,3),x,y
      real(8),allocatable:: ra0(:,:),ra(:,:),epi(:)
      integer,allocatable:: lscl(:)

c      open(10,file=cref,status='old')
c      read(10,*) na
c      read(10,*) h(1:3,1)
c      read(10,*) h(1:3,2)
c      read(10,*) h(1:3,3)
c      allocate(ra0(3,na),ra(3,na),epi(na),lscl(na))
c      do i=1,na
c        read(10,*) itmp,ra0(1:3,i)
c      enddo
c      close(10)
c
c      open(10,file=cfin,status='old')
c      read(10,*) itmp
c      if( itmp.ne.na ) stop ' [Error] itmp.ne.na !!!'
c      read(10,*) h(1:3,1)
c      read(10,*) h(1:3,2)
c      read(10,*) h(1:3,3)
c      do i=1,na
c        read(10,*) itmp,ra(1:3,i),tmp,epi(i)
c      enddo
c      close(10)

      open(10,file='cna010',status='old')
      read(10,*) na
      allocate(ra(3,na),epi(na))
      read(10,*) h(1:3,1)
      read(10,*) h(1:3,2)
      read(10,*) h(1:3,3)
      do i=1,na
        read(10,*) itmp,ra(1:3,i),tmp,epi(i)
      enddo
      close(10)

c      n=0
c      epot=0d0
c      do i=1,na
c        x= (ra0(1,i)-0.5d0)*h(1,1)
c        y= (ra0(2,i)-0.5d0)*h(2,2)
c        if( sqrt(x*x+y*y).lt.rcc ) then
c          n= n+1
c          epot= epot +epi(i)
c          lscl(n)= i
c        endif
c      enddo
c
c      open(10,file='core010',status='replace')
c      write(10,'(4i6)') n, 2, 0, 0
c      write(10,'(3es15.7)') h(1:3,1)
c      write(10,'(3es15.7)') h(1:3,2)
c      write(10,'(3es15.7)') h(1:3,3)
c      do i=1,n
c        j= lscl(i)
c        write(10,'(i4,3es15.7,2es12.4)') 1, ra0(1:3,j), 0d0, epi(j)
c      enddo
c      close(10)

      n=0
      epot= 0d0
      do i=1,na
        y= ra(2,i)
        if( y.gt.0.8d0 .or. y.lt.0.2d0 ) cycle
        n=n+1
        epot=epot +epi(i)
      enddo
      
      write(6,'(a,i6,es22.14)') ' n, epot =',n,epot

      end program get_core_energy
c=======================================================================
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "ifort -o getEcore get_core_energy.F"
c     End:
      
