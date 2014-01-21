      program move_bcc_screw
c-----------------------------------------------------------------------
c  Move screw dislocation in bcc crystal to next position.
c  Needs current and reference configurations sorted in the same order.
c-----------------------------------------------------------------------
c  Usage:
c    $ ./sortpmd pmdcur            <--- pmdcur must be sorted!!!
c    $ ./mvscrew pmdcur pmdref pmdnew
c-----------------------------------------------------------------------
      implicit none
      include '../params_au.h'
      integer:: natm,natm0,i,ndf,ia,ib,l
      real(8):: tmp,rcore(2),rcnew(2),rcut,h(3,3,0:1),dcore(2),dr
      real(8),allocatable:: tag(:),tag0(:),ra(:,:),ra0(:,:),va(:,:)
     &     ,eki(:),epi(:),strs(:,:,:),rdf(:,:),udf(:,:)

c.....Getting arguments
      integer,external:: iargc
      integer:: nargc
      character*128:: cfcur,cfref,cfout

c.....Get current and reference pmd files...............................
c.....Get file names
      nargc= iargc()
      if( nargc.ne.3 ) then
        write(6,'(a)') ' [Error] nargc.ne.3 !!!'
        write(6,'(a)') '   Usage: $ mvscrew pmdcur pmdref pmdnew'
        stop
      endif
      call getarg(1,cfcur)
      call getarg(2,cfref)
      call getarg(3,cfout)

c.....Read current atom positions
      open(10,file=trim(cfcur),status='old')
      read(10,*) natm
      allocate(tag(natm),tag0(natm),ra(3,natm),ra0(3,natm),va(3,natm)
     &     ,eki(natm),epi(natm),strs(3,3,natm),rdf(2,natm),udf(3,natm))
      read(10,*) (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
      read(10,*) (tag(i),ra(1:3,i),va(1:3,i)
     &     ,eki(i),epi(i),strs(1:3,1:3,i),i=1,natm)
      close(10)

c.....Read reference atom positions
      open(11,file=trim(cfref),status='old')
      read(11,*) natm0
      if( natm0.ne.natm ) stop ' [Error] natm0.ne.natm !!!'
      read(11,*) tmp,tmp,tmp
      read(11,*) tmp,tmp,tmp
      read(11,*) tmp,tmp,tmp
      read(11,*) tmp,tmp,tmp
      read(11,*) tmp,tmp,tmp
      read(11,*) tmp,tmp,tmp
      read(11,*) (tag0(i),ra0(1:3,i),tmp,tmp,tmp
     &     ,tmp,tmp, tmp,tmp,tmp, tmp,tmp,tmp
     &     ,tmp,tmp,tmp, i=1,natm0)
      close(11)

c.....Write akrcur
      open(80,file='akrcur',status='replace')
      write(80,'(i10,3i4)') natm, 5, 0, 0
      write(80,'(3es22.14)') h(1,1:3,0)
      write(80,'(3es22.14)') h(2,1:3,0)
      write(80,'(3es22.14)') h(3,1:3,0)
      do i=1,natm
        write(80,'(i4,3es15.7,5es12.4)') int(tag(i)),ra(1:3,i)
     &       ,eki(i),epi(i),ra(1:3,i)-ra0(1:3,i)
      enddo
      close(80)


c.....Ask the user and get estimated core postion
 10   write(6,'(a)') 'Please input core position in reduced rep. (x,y):'
      read(5,*) rcore(1:2)
      write(6,'(a,2es12.4)') ' Input core pos    =',rcore(1:2)

c.....Correct core position looking at reference configuration
      call correct_core_pos(rcore,natm,ra0,h)
      write(6,'(a,2es12.4)') ' Corrected core pos=',rcore(1:2)

c.....Get displacement field around current core position...............
c.....Ask the user and get cutoff radius
 20   write(6,'(a)') 'Please input cutoff radius of displacement field'
     &     //' (Ang.):'
      read(5,*) rcut
      rcut=rcut /0.5291772d0

c.....Get displacement field around current core
      call get_disp_field(natm,h,ra,ra0,rcut,rcore,ndf,rdf,udf)

c.....Ask the user and get next core position
 30   write(6,'(a)') 'Please input next core position (x,y):'
      read(5,*) rcnew(1:2)
      dcore(1)= (rcnew(1)-rcore(1)-anint(rcnew(1)-rcore(1)))*h(1,1,0)
      dcore(2)= (rcnew(2)-rcore(2)-anint(rcnew(2)-rcore(2)))*h(2,2,0)
      dr= sqrt(dcore(1)**2+dcore(2)**2)
      if( dr.gt.2.835d-10*sqrt(2d0)/sqrt(3d0) /aa *1.3d0 ) then
        write(6,'(a)') ' New core seems to be too far.'
        write(6,'(a,es12.4,a)') ' dr =',dr *0.5291772d0,' Ang.'
        goto 30
      endif
      call correct_core_pos(rcnew,natm,ra0,h)
      write(6,'(a,2es12.4)') ' Corrected core pos=',rcnew(1:2)

c.....Alter current positions of atoms according to the displacement
c.....field given by the new core position
      call set_disp_field(natm,h,ra,ra0,rcut,rcnew,ndf,rdf,udf)

      do i=1,natm
        if( ra(1,i).lt.0d0 ) ra(1,i)=ra(1,i) +1d0
        if( ra(2,i).lt.0d0 ) ra(2,i)=ra(2,i) +1d0
        if( ra(3,i).lt.0d0 ) ra(3,i)=ra(3,i) +1d0
        if( ra(1,i).ge.1d0 ) ra(1,i)=ra(1,i) -1d0
        if( ra(2,i).ge.1d0 ) ra(2,i)=ra(2,i) -1d0
        if( ra(3,i).ge.1d0 ) ra(3,i)=ra(3,i) -1d0
      enddo

c.....Write new configuration in pmd format
      open(20,file=trim(cfout),status='replace')
      write(20,'(i10)') natm
      write(20,'(3es22.14)') (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
      write(20,'(7es22.14,11es11.3)') (tag(i),ra(1:3,i),va(1:3,i)
     &     ,eki(i),epi(i),strs(1:3,1:3,i),i=1,natm)
      close(20)

c.....Write akrnew
      open(80,file='akrnew',status='replace')
      write(80,'(i10,3i4)') natm, 5, 0, 0
      write(80,'(3es22.14)') h(1,1:3,0)
      write(80,'(3es22.14)') h(2,1:3,0)
      write(80,'(3es22.14)') h(3,1:3,0)
      do i=1,natm
        write(80,'(i4,3es15.7,5es12.4)') int(tag(i)),ra(1:3,i)
     &       ,eki(i),epi(i),ra(1:3,i)-ra0(1:3,i)
      enddo
      close(80)

      end program move_bcc_screw
c=======================================================================
      subroutine correct_core_pos(rcore,natm,ra,h)
c
c  Correct core position looking at the reference configuration.
c
      implicit none
      include '../params_au.h'
      integer,intent(in):: natm
      real(8),intent(in):: ra(3,natm),h(3,3)
      real(8),intent(inout):: rcore(2)

      real(8),parameter:: a  = 2.835d-10 /aa
      real(8),parameter:: a1 = 1.41421356d0 /1.7320508d0 *a
      real(8),parameter:: a2 = 1.7320508d0 *a1
c      real(8),parameter:: rc = (a1+a2)/2
      real(8),parameter:: rc = a1
      real(8),parameter:: rc2= rc*rc

      integer:: i,j,n,nn,lsn(3)
      real(8):: x,y,r2
      logical:: lexist

c.....Look for 3 nearest atoms which construct a triangle
      nn= 0
      lsn(1:3)= 0
      do i=1,natm
        x= (ra(1,i)-rcore(1) -anint(ra(1,i)-rcore(1)))*h(1,1)
        y= (ra(2,i)-rcore(2) -anint(ra(2,i)-rcore(2)))*h(2,2)
        r2= x*x +y*y
        if( r2.gt.rc2 ) cycle
        lexist=.false.
        do n=1,nn
          j=lsn(n)
          if(  nint(ra(1,i)*h(1,1)*100d0).eq.nint(ra(1,j)*h(1,1)*100d0)
     &         .and.
     &         nint(ra(2,i)*h(2,2)*100d0).eq.nint(ra(2,j)*h(2,2)*100d0)
     &         ) then
            lexist=.true.
            exit
          endif
        enddo
        if( lexist ) cycle
        nn=nn+1
        lsn(nn)= i
        if( nn.eq.3 ) exit
      enddo

      if( nn.ne.3 ) then
        write(6,'(a)') ' [Error] nn.ne.3 !!!'
        write(6,'(a,i4)') '  nn=',nn
        do n=1,nn
          write(6,'(a,2i8,2es12.4)') '  n,lsn(n),ra(1:2,j)='
     &         ,n,lsn(n),ra(1:2,lsn(n))
        enddo
        stop
      endif

c.....Get the center of the triangle
      rcore(1:2)= 0d0
      do n=1,nn
        j=lsn(n)
        rcore(1:2)= rcore(1:2) +ra(1:2,j)
      enddo
      rcore(1:2)=rcore(1:2)/3

      return
      end subroutine correct_core_pos
c=======================================================================
      subroutine get_disp_field(natm,h,ra,ra0,rcut,rcore,ndf,rdf,udf)
      implicit none
      integer,intent(in):: natm
      real(8),intent(in):: h(3,3),ra(3,natm),ra0(3,natm),rcore(2),rcut
      integer,intent(out):: ndf
      real(8),intent(out):: rdf(2,natm),udf(3,natm)

      integer:: i,n
      real(8):: x,y,rc2,r2
      logical:: lexist

      rc2= rcut*rcut
      ndf= 0
      do i=1,natm
        x= (ra0(1,i)-rcore(1) -anint(ra0(1,i)-rcore(1)))*h(1,1)
        y= (ra0(2,i)-rcore(2) -anint(ra0(2,i)-rcore(2)))*h(2,2)
        r2= x*x +y*y
        if( r2.gt.rc2 ) cycle
        lexist=.false.
        do n=1,ndf
          if(  nint(x*100d0).eq.nint(rdf(1,n)*100d0)
     &         .and.
     &         nint(y*100d0).eq.nint(rdf(2,n)*100d0) ) then
            lexist=.true.
            exit
          endif
        enddo
        if( lexist ) cycle
        ndf=ndf+1
        rdf(1,ndf)= x
        rdf(2,ndf)= y
        udf(1:3,ndf)= ra(1:3,i)-ra0(1:3,i)-anint(ra(1:3,i)-ra0(1:3,i))
c        write(6,'(a,2i6,5es12.4)') ' i,ndf,x,y,udf(1:3)='
c     &       ,i,ndf,x,y,udf(1:3,ndf)
      enddo
      write(6,'(a,i8)') ' Num of registered atoms in rcut=',ndf

      return
      end subroutine get_disp_field
c=======================================================================
      subroutine set_disp_field(natm,h,ra,ra0,rcut,rcore,ndf,rdf,udf)
      implicit none
      integer,intent(in):: natm,ndf
      real(8),intent(in):: h(3,3),ra0(3,natm),rcore(2),rcut
     &     ,rdf(2,natm),udf(3,natm)
      real(8),intent(out):: ra(3,natm)

      integer:: i,n
      real(8):: x,y,rc2,r2
      logical:: lexist

      rc2= rcut*rcut

      do i=1,natm
        x= (ra0(1,i)-rcore(1) -anint(ra0(1,i)-rcore(1)))*h(1,1)
        y= (ra0(2,i)-rcore(2) -anint(ra0(2,i)-rcore(2)))*h(2,2)
        r2= x*x +y*y
        if( r2.gt.rc2 ) cycle
        lexist=.false.
        do n=1,ndf
          if(  nint(x*100d0).eq.nint(rdf(1,n)*100d0)
     &         .and.
     &         nint(y*100d0).eq.nint(rdf(2,n)*100d0) ) then
            lexist=.true.
            exit
          endif
        enddo
        if( .not. lexist ) then
          write(6,'(a)') ' [Error] There is no correspondance...'
          write(6,'(a,i8,3es12.4)') '  i,x,y,r2=',i,x,y,r2
          stop
        endif
        ra(1:3,i)= ra0(1:3,i) +udf(1:3,n)
      enddo

      return
      end subroutine set_disp_field
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "ifort -o mvscrew move_bcc_screw.F"
c     End:
