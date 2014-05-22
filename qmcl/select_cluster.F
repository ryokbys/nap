      program select_cluster
c
c  Select QM cluster to be used in qmcl program.
c  Write out.idclst and out.idqm files.
c
      use variables
      implicit none
      include 'params_unit.h'
      integer:: i,j,m
      real(8):: xt(3),xr(3),r2,rc2,rcqm2

      call read_input(10,'in.qmcl')
      
      call read_pmd_ascii(20,"qmclref"
     &     ,namax,natm,h,tag,ra,va,eki,epi,strs)
c.....Some conversions
      do i=1,natm
c.....Species of atom-i
        nis= max(nis,int(tag(i)))
        if( nis.gt.nismax ) then
          write(6,'(a)') ' [qmcl] nis.gt.nismax !'
          stop
        endif
      enddo
      write(6,'(a,i8)') ' Number of atoms   =',natm
      write(6,'(a,i8)') ' Number of species =',nis

      if( trim(cqmshape).eq.'circle' ) then
        call mk_clst_circle()
      elseif( cqmshape(1:4).eq.'rect' ) then
        call mk_clst_rect()
      else
        write(6,'(a)')  ' There is no such option '//cqmshape//' !!!'
        stop
      endif

      write(6,'(a,i6)') ' num of cluster atoms      =',nac
      write(6,'(a,i6)') ' num of QM atoms           =',naqm

      open(90,file='out.idclst',status='replace')
      write(90,'(i6)') nac
      do i=1,nac
        write(90,'(i8)') idclst(i)
      enddo
      close(90)

      open(91,file='out.idqm',status='replace')
      write(91,'(i6)') naqm
      do i=1,naqm
        write(91,'(i8)') idqm(i)
      enddo
      close(91)

      end program select_cluster
c=======================================================================
      subroutine mk_clst_circle()
      use variables
      implicit none
      include 'params_unit.h'
      integer:: i,j
      real(8):: rc2,rcqm2,xt(3),xr(3),r2
      
      rc2= rclst**2
      rcqm2= rcqm**2

c.....From the center position and radius, create the cluster
      nac= 0
      naqm= 0
      do i=1,natm
        xt(1:3)= ra(1:3,i)-xc(1:3) -anint(ra(1:3,i)-xc(1:3))
c        xt(1:3)= ra0(1:3,i)-xc(1:3) -anint(ra0(1:3,i)-xc(1:3))
        xr(1:3)= h(1:3,1,0)*xt(1) +h(1:3,2,0)*xt(2) +h(1:3,3,0)*xt(3)
c.....If xc(?) is minus, no length needed in ?-direction
        if( xc(1).lt.0d0 ) xr(1)= 0d0
        if( xc(2).lt.0d0 ) xr(2)= 0d0
        if( xc(3).lt.0d0 ) xr(3)= 0d0
        r2= xr(1)**2 +xr(2)**2 +xr(3)**2
        if( r2.ge.rc2 ) cycle
        nac=nac +1
        idclst(nac)= i
c        write(6,'(a,2i6)') ' nac,idclst=',nac,idclst(nac)
      enddo

      call sort_idclst_by_is()

      do i=1,nac
        j= idclst(i)
        xt(1:3)= ra(1:3,j)-xc(1:3) -anint(ra(1:3,j)-xc(1:3))
        xr(1:3)= h(1:3,1,0)*xt(1) +h(1:3,2,0)*xt(2) +h(1:3,3,0)*xt(3)
c.....If xc(?) is minus, no length needed in ?-direction
        if( xc(1).lt.0d0 ) xr(1)= 0d0
        if( xc(2).lt.0d0 ) xr(2)= 0d0
        if( xc(3).lt.0d0 ) xr(3)= 0d0
        r2= xr(1)**2 +xr(2)**2 +xr(3)**2
        if( r2.ge.rcqm2 ) cycle
        naqm=naqm +1
        idqm(naqm)= i
c        write(6,'(a,2i6)') ' naqm,idqm(naqm)=',naqm,idqm(naqm)
      enddo

      end subroutine mk_clst_circle
c=======================================================================
      subroutine mk_clst_rect()
      use variables
      implicit none
      include 'params_unit.h'
      integer:: i,j
      real(8):: rc2,rcqm2,xt(3),xr(3),r2

      write(6,'(a)') ' creating rectangular QM regions...'
      
c.....Length: AA
      rclx= rclx /2
      rcly= rcly /2
      rclz= rclz /2
      rcqmx= rcqmx /2
      rcqmy= rcqmy /2
      rcqmz= rcqmz /2

      if( rclx.lt.0d0 ) rclx= h(1,1,0)
      if( rcly.lt.0d0 ) rcly= h(2,2,0)
      if( rclz.lt.0d0 ) rclz= h(3,3,0)
      if( rcqmx.lt.0d0 ) rcqmx= h(1,1,0)
      if( rcqmy.lt.0d0 ) rcqmy= h(2,2,0)
      if( rcqmz.lt.0d0 ) rcqmz= h(3,3,0)

c.....From the center position and radius, create the cluster
      nac= 0
      naqm= 0
      do i=1,natm
        xt(1:3)= ra(1:3,i)-xc(1:3) -anint(ra(1:3,i)-xc(1:3))
c        xt(1:3)= ra0(1:3,i)-xc(1:3) -anint(ra0(1:3,i)-xc(1:3))
        xr(1:3)= h(1:3,1,0)*xt(1) +h(1:3,2,0)*xt(2) +h(1:3,3,0)*xt(3)
c.....If xc(?) is minus, no length needed in ?-direction
        if( xc(1).lt.0d0 ) xr(1)= 0d0
        if( xc(2).lt.0d0 ) xr(2)= 0d0
        if( xc(3).lt.0d0 ) xr(3)= 0d0
        if(  abs(xr(1)).gt.rclx .or.
     &       abs(xr(2)).gt.rcly .or.
     &       abs(xr(3)).gt.rclz ) cycle
c        r2= xr(1)**2 +xr(2)**2 +xr(3)**2
c        if( r2.ge.rc2 ) cycle
        nac=nac +1
        idclst(nac)= i
      enddo

      call sort_idclst_by_is()

      do i=1,nac
        j= idclst(i)
        xt(1:3)= ra(1:3,j)-xc(1:3) -anint(ra(1:3,j)-xc(1:3))
        xr(1:3)= h(1:3,1,0)*xt(1) +h(1:3,2,0)*xt(2) +h(1:3,3,0)*xt(3)
c.....If xc(?) is minus, no length needed in ?-direction
        if( xc(1).lt.0d0 ) xr(1)= 0d0
        if( xc(2).lt.0d0 ) xr(2)= 0d0
        if( xc(3).lt.0d0 ) xr(3)= 0d0
        if(  abs(xr(1)).gt.rcqmx .or.
     &       abs(xr(2)).gt.rcqmy .or.
     &       abs(xr(3)).gt.rcqmz ) cycle
        naqm=naqm +1
        idqm(naqm)= i
      enddo

      end subroutine mk_clst_rect
c=======================================================================
      subroutine sort_idclst_by_is()
      use variables
      implicit none
      integer,allocatable:: isc(:)
      integer:: ic,i,j,l,m

      allocate(isc(nac))

      do ic=1,nac
        isc(ic)= int(tag(idclst(ic)))
      enddo

c.....bubble sort because nac is small
      do i=1,nac-1
        do j=2,nac-i+1
          if( isc(j-1).gt.isc(j) ) then
c.....swap j and j-1 elements
            l=isc(j)
            isc(j)= isc(j-1)
            isc(j-1)= l
            m=idclst(j)
            idclst(j)= idclst(j-1)
            idclst(j-1)= m
          endif
        enddo
      enddo

      deallocate(isc)
      end subroutine sort_idclst_by_is
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make selcl"
c     End:
