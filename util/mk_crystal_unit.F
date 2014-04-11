      program mk_crystal_unit
c-----------------------------------------------------------------------
c  Create arbitrary unit cell from Bravais lattice.
c-----------------------------------------------------------------------
      implicit none
      integer,parameter:: nmax= 10000
      integer:: i,j,k,l,m,n,istrct,ixyz,ncell,nua,i1,i2,i3,inc,natm
      character:: clattice*10
      integer:: ib1(3),ib2(3),ib3(3),imax(3),imin(3)
      real(8):: a1(3),a2(3),a3(3),h0(3,3),ua(3,10),al
     &     ,b1(3),b2(3),b3(3),ra(3,nmax),s(3),bb1,bb2,bb3,rb(3,nmax)
     &     ,ab1,ab2,ab3
      logical:: lb1,lb2,lb3

      real(8),parameter:: small= 1d-8

      open(10,file='out.mkxtlunit',status='replace')
      
c.....Choose Bravais lattice
 1    write(6,'(a)') "Choose a Bravais lattice from the list below:"
      write(6,'(a)') " 1: BCC"
      write(6,'(a)') " 2: FCC"
      write(6,'(a)') " 3: HCP"
      write(6,'(a)') " 4: Diamond"
      write(6,'(a)') " 5: Wultzite"
      read(5,*) istrct

c.....Lattice constant
      write(6,'(a)') "Input lattice constant in Angstrom:"
      read(5,*) al

c.....Decide new lattice vectors
 3    write(6,'(a)') "Input lattice vectors (e.g., [-121]= -1 2 1):"
      write(6,'(a)',advance='no') "  b1: "
      read(5,*) ib1(1:3)
      write(6,'(a)',advance='no') "  b2: "
      read(5,*) ib2(1:3)
      write(6,'(a)',advance='no') "  b3: "
      read(5,*) ib3(1:3)
c.....Multiplying factors
 4    write(6,'(a)') "Input coefficient for these vectors:"
      write(6,'(a)',advance='no') "  c1: "
      read(5,*) ab1
      write(6,'(a)',advance='no') "  c2: "
      read(5,*) ab2
      write(6,'(a)',advance='no') "  c3: "
      read(5,*) ab3

      imax(1:3)= 0
      imin(1:3)= 0
      do ixyz=1,3
        imax(ixyz)= max(0,ib1(ixyz))
        imax(ixyz)= max(imax(ixyz),ib2(ixyz))
        imax(ixyz)= max(imax(ixyz),ib3(ixyz))
        imax(ixyz)= max(imax(ixyz),ib1(ixyz)+ib2(ixyz))
        imax(ixyz)= max(imax(ixyz),ib1(ixyz)+ib3(ixyz))
        imax(ixyz)= max(imax(ixyz),ib2(ixyz)+ib3(ixyz))
        imax(ixyz)= max(imax(ixyz),ib1(ixyz)+ib2(ixyz)+ib3(ixyz))
        imin(ixyz)= min(0,ib1(ixyz))
        imin(ixyz)= min(imin(ixyz),ib2(ixyz))
        imin(ixyz)= min(imin(ixyz),ib3(ixyz))
        imin(ixyz)= min(imin(ixyz),ib1(ixyz)+ib2(ixyz))
        imin(ixyz)= min(imin(ixyz),ib1(ixyz)+ib3(ixyz))
        imin(ixyz)= min(imin(ixyz),ib2(ixyz)+ib3(ixyz))
        imin(ixyz)= min(imin(ixyz),ib1(ixyz)+ib2(ixyz)+ib3(ixyz))
      enddo
      write(6,'(a,3i3)') " imax(1:3)=",imax(1:3)
      write(6,'(a,3i3)') " imin(1:3)=",imin(1:3)
      ncell=(imax(1)-imin(1)+1)*(imax(2)-imin(2)+1)*(imax(3)-imin(3)+1)
      write(6,'(a,i4)')  " ncell    =",ncell

c.....Make large enough crystal to be extracted unit cell
      if( istrct.eq.1 ) then    ! BCC
        write(6,'(a)') " Making BCC..."
        a1(1:3) =(/ 1.0d0, 0.0d0, 0.0d0 /)
        a2(1:3) =(/ 0.0d0, 1.0d0, 0.0d0 /)
        a3(1:3) =(/ 0.0d0, 0.0d0, 1.0d0 /)

        nua=2
        ua(1:3,1)=(/ 0.0d0, 0.0d0, 0.0d0 /)
        ua(1:3,2)=(/ 0.5d0, 0.5d0, 0.5d0 /)
      elseif( istrct.eq.2 ) then
        write(6,'(a)') " Making FCC..."
        a1(1:3) =(/ 1.0d0, 0.0d0, 0.0d0 /)
        a2(1:3) =(/ 0.0d0, 1.0d0, 0.0d0 /)
        a3(1:3) =(/ 0.0d0, 0.0d0, 1.0d0 /)

        nua=4
        ua(1:3,1)=(/ 0.0d0, 0.0d0, 0.0d0 /)
        ua(1:3,2)=(/ 0.5d0, 0.5d0, 0.0d0 /)
        ua(1:3,1)=(/ 0.5d0, 0.0d0, 0.5d0 /)
        ua(1:3,2)=(/ 0.0d0, 0.5d0, 0.5d0 /)
      elseif( istrct.eq.4 ) then
        write(6,'(a)') " Making diamond structure..."
        a1(1:3) =(/ 1.0d0, 0.0d0, 0.0d0 /)
        a2(1:3) =(/ 0.0d0, 1.0d0, 0.0d0 /)
        a3(1:3) =(/ 0.0d0, 0.0d0, 1.0d0 /)

        nua=8
        ua(1:3,1)=(/ 0.00d0, 0.00d0, 0.00d0 /)
        ua(1:3,2)=(/ 0.50d0, 0.50d0, 0.00d0 /)
        ua(1:3,3)=(/ 0.50d0, 0.00d0, 0.50d0 /)
        ua(1:3,4)=(/ 0.00d0, 0.50d0, 0.50d0 /)
        ua(1:3,5)=(/ 0.25d0, 0.25d0, 0.25d0 /)
        ua(1:3,6)=(/ 0.75d0, 0.75d0, 0.25d0 /)
        ua(1:3,7)=(/ 0.75d0, 0.25d0, 0.75d0 /)
        ua(1:3,8)=(/ 0.25d0, 0.75d0, 0.75d0 /)
      else
        write(6,'(a)') ' Sorry, the structure is not available...'
        stop
      endif

c.....Get absolute atom positions
      inc=0
      do i3=imin(3),imax(3)
      do i2=imin(2),imax(2)
      do i1=imin(1),imax(1)
        do m=1,nua
          inc= inc +1
          if( inc.gt.nmax ) stop ' [Error] inc.gt.nmax!'
          ra(1,inc)= a1(1)*ua(1,m) +a2(1)*ua(2,m) +a3(1)*ua(3,m)
          ra(2,inc)= a1(2)*ua(1,m) +a2(2)*ua(2,m) +a3(2)*ua(3,m)
          ra(3,inc)= a1(3)*ua(1,m) +a2(3)*ua(2,m) +a3(3)*ua(3,m)
          ra(1,inc)= ra(1,inc) +a1(1)*i1 +a2(1)*i2 +a3(1)*i3
          ra(2,inc)= ra(2,inc) +a1(2)*i1 +a2(2)*i2 +a3(2)*i3
          ra(3,inc)= ra(3,inc) +a1(3)*i1 +a2(3)*i2 +a3(3)*i3
c          write(6,'(a,i5,3es12.4)') " inc,ra=",inc,ra(1:3,inc)
        enddo
      enddo
      enddo
      enddo
      natm= inc

c.....New unit vectors
      b1(1)= (a1(1)*ib1(1) +a2(1)*ib1(2) +a3(1)*ib1(3))*ab1
      b1(2)= (a1(2)*ib1(1) +a2(2)*ib1(2) +a3(2)*ib1(3))*ab1
      b1(3)= (a1(3)*ib1(1) +a2(3)*ib1(2) +a3(3)*ib1(3))*ab1
      b2(1)= (a1(1)*ib2(1) +a2(1)*ib2(2) +a3(1)*ib2(3))*ab2
      b2(2)= (a1(2)*ib2(1) +a2(2)*ib2(2) +a3(2)*ib2(3))*ab2
      b2(3)= (a1(3)*ib2(1) +a2(3)*ib2(2) +a3(3)*ib2(3))*ab2
      b3(1)= (a1(1)*ib3(1) +a2(1)*ib3(2) +a3(1)*ib3(3))*ab3
      b3(2)= (a1(2)*ib3(1) +a2(2)*ib3(2) +a3(2)*ib3(3))*ab3
      b3(3)= (a1(3)*ib3(1) +a2(3)*ib3(2) +a3(3)*ib3(3))*ab3
      bb1= b1(1)**2 +b1(2)**2 +b1(3)**2
      bb2= b2(1)**2 +b2(2)**2 +b2(3)**2
      bb3= b3(1)**2 +b3(2)**2 +b3(3)**2
c      write(6,'(a,es12.4)') " bb1=",bb1
c      write(6,'(a,es12.4)') " bb2=",bb2
c      write(6,'(a,es12.4)') " bb3=",bb3

c.....Check if there is an atom at the point of the new unit vectors
      lb1=.false.
      lb2=.false.
      lb3=.false.
      do i=1,natm
        s(1)= (ra(1,i)*b1(1) +ra(2,i)*b1(2) +ra(3,i)*b1(3))/bb1
        s(2)= (ra(1,i)*b2(1) +ra(2,i)*b2(2) +ra(3,i)*b2(3))/bb2
        s(3)= (ra(1,i)*b3(1) +ra(2,i)*b3(2) +ra(3,i)*b3(3))/bb3
        if( ( nint(s(1)*100).eq.100 .and. nint(s(2)*100).eq.0 .and.
     &        nint(s(3)*100).eq.0 ) ) lb1=.true.
        if( ( nint(s(1)*100).eq.0 .and. nint(s(2)*100).eq.100 .and.
     &        nint(s(3)*100).eq.0 ) ) lb2=.true.
        if( ( nint(s(1)*100).eq.0 .and. nint(s(2)*100).eq.0 .and.
     &        nint(s(3)*100).eq.100 ) ) lb3=.true.
      enddo
      if( .not.lb1 .or. .not.lb2 .or. .not.lb3 ) then
        write(6,'(a)') ' There is no atom at the apex of '//
     &       'the unit cell.'
        write(6,*) '  lb1,lb2,lb3=',lb1,lb2,lb3
        goto 3
      endif

      inc=0
      do i=1,natm
        s(1)= (ra(1,i)*b1(1) +ra(2,i)*b1(2) +ra(3,i)*b1(3))/bb1
        s(2)= (ra(1,i)*b2(1) +ra(2,i)*b2(2) +ra(3,i)*b2(3))/bb2
        s(3)= (ra(1,i)*b3(1) +ra(2,i)*b3(2) +ra(3,i)*b3(3))/bb3
c        write(6,'(a,i5,3es12.4)') " i,s(1:3)=",i,s(1:3)
        if(  s(1).ge.0d0 .and. s(1).lt.1d0 .and.
     &       s(2).ge.0d0 .and. s(2).lt.1d0 .and.
     &       s(3).ge.0d0 .and. s(3).lt.1d0 ) then
          inc=inc +1
          if( inc.gt.nmax ) stop ' [Error] inc.gt.nmax!'
          rb(1:3,inc)= s(1:3)
          write(6,'(a,i5,3es11.3,a,3es11.3)') " inc,ra,rb="
     &         ,inc,ra(1:3,i),' ---> ',s(1:3)
        endif
      enddo

c.....scaling
      b1(1:3)= al*b1(1:3)
      b2(1:3)= al*b2(1:3)
      b3(1:3)= al*b3(1:3)

c.....Output unit cell information
      write(6,'(a)') ''
      write(6,'(a)') " Write following data in your code..."
      write(6,'(a,i5)') " natm=",inc
      write(6,'(a)') " h-mat(b1,b2,b3)="
      write(6,'(3es22.14)') b1(1),b2(1),b3(1)
      write(6,'(3es22.14)') b1(2),b2(2),b3(2)
      write(6,'(3es22.14)') b1(3),b2(3),b3(3)
      write(6,'(a,es22.14)') " |b1|=",sqrt(bb1)*al
      write(6,'(a,es22.14)') " |b2|=",sqrt(bb2)*al
      write(6,'(a,es22.14)') " |b3|=",sqrt(bb3)*al
      do i=1,inc
        write(6,'(i5,3es22.14)') i,rb(1:3,i)
      enddo

c.....Output this configuration in akr format
      open(10,file='akr0000',status='replace')
      write(10,'(i10,3i4)') inc,0,0,0
      write(10,'(3es15.7)') sqrt(bb1)*al, 0d0, 0d0
      write(10,'(3es15.7)') 0d0, sqrt(bb2)*al, 0d0
      write(10,'(3es15.7)') 0d0, 0d0, sqrt(bb3)*al
      do i=1,inc
        write(10,'(4es22.14)') 1d0+i*1d-14,rb(1:3,i)
      enddo
      close(10)

      end program mk_crystal_unit
c=======================================================================
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "ifort -o mkxtlunit mk_crystal_unit.F"
c     End:
