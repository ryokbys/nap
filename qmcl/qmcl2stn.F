      program qmcl2stn
c-----------------------------------------------------------------------
c  This program converts qmcl file to stn which includes
c  atomic strain in auxiliary data of qmcl format using qmclref file
c  as reference atomic positions.
c  The ordering of atoms in qmclref and qmcl??? files should be the same
c  and each line must correspond one to one.
c-----------------------------------------------------------------------
c  USAGE:
c    $ ./qmcl2stn  qmcl??? stn???
c-----------------------------------------------------------------------
      use variables
      implicit none
      integer:: i,nargc,ia,ib
      integer,external:: iargv
      character:: cin*128,cout*128,ctmp*3
      character*7,parameter:: cref='qmclref' 

      nargc= iargc()
      if( nargc.ne.2 ) then
        write(6,'(a)') ' [Error] Num of arguments too small !!!'
        write(6,'(a)') '   Usage: $ ./qmcl2stn qmcl000 stn000'
        stop
      endif
      call getarg(1,cin)
      call getarg(2,cout)

      if( cout(1:4).eq.'qmcl' ) then
        write(6,'(a)') ' Do you really want to overwrite '//cout//'?'
        read(5,*) ctmp
        if( ctmp(1:1).eq.'y' .or. ctmp(1:1).eq.'Y' ) goto 1
        stop
      endif
 1    continue

      call read_pmd_ascii(10,trim(cref)
     &     ,namax,natm,h0,hunit,tag,ra0,va,eki,epi,strs)
      call read_pmd_ascii(11,trim(cin)
     &     ,namax,natm,h,hunit,tag,ra,va,eki,epi,strs)

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

c.....Calc atomic strain
      call get_atomic_strain()

      open(12,file=trim(cout),status='replace')
      write(12,'(es15.7)') hunit
      write(12,'(3es15.7)') ((h(ia,ib,0)/hunit,ia=1,3),ib=1,3)
      write(12,'(i10,3i4)') natm, 8, 0, 0
      do i=1,natm
        write(12,'(i3,3es12.4,10es11.3)') int(tag(i)),ra(1:3,i)
     &       ,eki(i),epi(i)
     &       ,stn(1,1,i),stn(2,2,i),stn(3,3,i)
     &       ,stn(2,3,i),stn(3,1,i),stn(1,2,i)
      enddo
      close(12)

      end program qmcl2stn
c=======================================================================
      subroutine get_atomic_strain()
c-----------------------------------------------------------------------
c  Calculate strain tensor according to
c    F.Shimizu, et al., Mater. Trans. 48 (2007) pp.2923-2927.
c-----------------------------------------------------------------------
      use variables
      implicit none
      include 'params_unit.h'

      integer:: i,j,l,ixyz,jxyz
      real(8):: x,y,z,rij2
      logical,save:: l1st=.true.
      real(8),save:: rc2
      real(8),allocatable,save:: xi(:),xi0(:),xij(:),xij0(:),ddt(:,:)
     &     ,dd0i(:,:)

      if( l1st ) then
        allocate(xi(3),xi0(3),xij(3),xij0(3),ddt(3,3),dd0i(3,3))
c.....rc2: cutoff for 1st nearest neighbors in FCC (Al)
        rc2= 4.0425d0 *(1d0 +1d0/sqrt(2d0))/2
        rc2= rc2*rc2
cc.....rc2: cutoff for 1st nearest neighbors in BCC-Fe
c        rc2= 2.8553d0 *(1d0 +sqrt(3d0)/2) /2
c        rc2= rc2*rc2
        l1st=.false.
      endif

      call mklspr(sqrt(rc2))

c-----initialize strain tensor
      stn(1:3,1:3,1:natm)= 0d0

      do i=1,natm
        xi(1:3)=ra(1:3,i)
        xi0(1:3)=ra0(1:3,i)
        ddt(1:3,1:3)= 0d0
        dd0i(1:3,1:3)= 0d0
        do l=1,lspr(0,i)
          j=lspr(l,i)
          xij(1)= (ra(1,j)-xi(1)-anint(ra(1,j)-xi(1)))*h(1,1,0)
          xij(2)= (ra(2,j)-xi(2)-anint(ra(2,j)-xi(2)))*h(2,2,0)
          xij(3)= (ra(3,j)-xi(3)-anint(ra(3,j)-xi(3)))*h(3,3,0)
          rij2= xij(1)**2 +xij(2)**2 +xij(3)**2
c---------only 1st nearest neighbors
          if(rij2.gt.rc2) cycle
          xij0(1)= (ra0(1,j)-xi0(1))*h0(1,1,0)
          xij0(2)= (ra0(2,j)-xi0(2))*h0(2,2,0)
          xij0(3)= (ra0(3,j)-xi0(3))*h0(3,3,0)
          do jxyz=1,3
            do ixyz=1,3
              stn(ixyz,jxyz,i)= stn(ixyz,jxyz,i) +xij0(ixyz)*xij(jxyz)
              ddt(ixyz,jxyz)= ddt(ixyz,jxyz) +xij0(ixyz)*xij0(jxyz)
            enddo
          enddo
        enddo
c        if( i.eq.1 ) then
c          write(6,'(a,9es11.3)') ' stn=',stn(1:3,1:3,i)
c          write(6,'(a,9es11.3)') ' ddt=',ddt(1:3,1:3)
c        endif
c-------inverse of ddt ---> dd0i
c        call ludc_inv(2,ddt(1:2,1:2),dd0i(1:2,1:2))
        call ludc_inv(3,ddt,dd0i)
c        if( i.eq.1 ) then
c          write(6,'(a,9es11.3)') ' dd0i=',dd0i(1:3,1:3)
c        endif
c-------J =[dd0i]*dd =[dd0i]*stn
        ddt(1:3,1:3)= stn(1:3,1:3,i)
        do jxyz=1,3
          do ixyz=1,3
            stn(ixyz,jxyz,i)= dd0i(1,ixyz)*ddt(jxyz,1)
     &           +dd0i(2,ixyz)*ddt(jxyz,2) +dd0i(3,ixyz)*ddt(jxyz,3)
          enddo
        enddo
c        if( i.eq.1 ) then
c          write(6,'(a,9es11.3)') ' stn=',stn(1:3,1:3,i)
c        endif
c-------J^{T}J
        ddt(1:3,1:3)= stn(1:3,1:3,i)
        do jxyz=1,3
          do ixyz=1,3
            stn(ixyz,jxyz,i)= ddt(ixyz,1)*ddt(jxyz,1)
     &           +ddt(ixyz,2)*ddt(jxyz,2) +ddt(ixyz,3)*ddt(jxyz,3)
          enddo
        enddo
c-------Lagrange strain tensor
        stn(1,1,i)= stn(1,1,i) -1d0
        stn(2,2,i)= stn(2,2,i) -1d0
        stn(3,3,i)= stn(3,3,i) -1d0
        stn(1:3,1:3,i)= 0.5d0 *stn(1:3,1:3,i)
c        if( i.eq.1 ) then
c          write(6,'(a,9es11.3)') ' stn=',stn(1:3,1:3,i)
c        endif
      enddo

      end subroutine get_atomic_strain
c=======================================================================
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make qmcl2stn"
c     End:
