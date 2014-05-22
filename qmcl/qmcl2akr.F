      program qmcl2akr
c-----------------------------------------------------------------------
c USAGE:
c   $ ./qmcl2akr < qmcl??? > akr???
c-----------------------------------------------------------------------
      implicit none
      integer:: natm,i,ia,ib,l,ncpx,ncpy,ncpz,natmt,m,ix,iy,iz
     &     ,iqm,ic,mis,nac,naqm
      integer,allocatable:: idqm(:),idclst(:)
      real(8):: h(3,3,0:1),ht(3,3,0:1),dcpx,dcpy,dcpz,hunit
      real(8),allocatable:: ra(:,:),tag(:),va(:,:),eki(:),epi(:)
     &     ,strs(:,:,:)
      real(8),allocatable:: rat(:,:),vat(:,:),tagt(:),ekit(:),epit(:)
     &     ,strst(:,:,:)
      real(8),parameter:: rc   =  4.4d0

c.....Read out.idclst and out.idqm
      open(10,file='out.idclst')
      read(10,*) nac
      allocate(idclst(nac))
      do i=1,nac
        read(10,*) idclst(i)
      enddo
      close(10)

      open(11,file='out.idqm')
      read(11,*) naqm
      allocate(idqm(naqm))
      do i=1,naqm
        read(11,*) idqm(i)
      enddo
      close(11)

c.....Read qmcl file from standard input
      read(5,*) hunit
      read(5,*) (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
      h(1:3,1:3,0:1)= h(1:3,1:3,0:1) *hunit
      read(5,*) natm
      allocate(ra(3,natm),tag(natm),va(3,natm),eki(natm),epi(natm)
     &     ,strs(3,3,natm))
      mis= 0
      do i=1,natm
        read(5,*) tag(i),ra(1:3,i),va(1:3,i)
     &     ,eki(i),epi(i),strs(1:3,1:3,i)
        mis= max(mis,int(tag(i)))
      enddo

cc.....Change tag of qm atoms
c      do iqm=1,naqm
c        ic= idqm(iqm)
c        i= idclst(ic)
c        tag(i)= tag(i) +mis
c      enddo

      ncpx=1
      ncpy=1
      ncpz=1
c      if( h(1,1,0).lt.2*rc ) then
c        ncpx= 1+int(2*rc/h(1,1,0))
c      elseif( h(2,2,0).lt.2*rc ) then
c        ncpy= 1+int(2*rc/h(2,2,0))
c      elseif( h(3,3,0).lt.2*rc ) then
c        ncpz= 1+int(2*rc/h(3,3,0))
c      endif
      natmt= natm*ncpx*ncpy*ncpz
c      write(6,'(4i10)') natmt,ncpx,ncpy,ncpz
c      stop
      ht(1:3,1:3,0:1)= 0d0
      ht(1,1,0)=  h(1,1,0)*ncpx
      ht(2,2,0)=  h(2,2,0)*ncpy
      ht(3,3,0)=  h(3,3,0)*ncpz
      allocate(rat(3,natmt),vat(3,natmt),tagt(natmt)
     &     ,ekit(natmt),epit(natmt),strst(3,3,natmt))

      m=0
      do iz=1,ncpz
        dcpz= dble(iz-1)/ncpz
        do iy=1,ncpy
          dcpy= dble(iy-1)/ncpy
          do ix=1,ncpx
            dcpx= dble(ix-1)/ncpx
            do i=1,natm
              m=m+1
              tagt(m)= tag(i)
              rat(1,m)= ra(1,i)/ncpx +dcpx
              rat(2,m)= ra(2,i)/ncpy +dcpy
              rat(3,m)= ra(3,i)/ncpz +dcpz
              vat(1:3,m)= va(1:3,i)
              ekit(m)= eki(i)
              epit(m)= epi(i)
              strst(1:3,1:3,m)= strs(1:3,1:3,i)
            enddo
          enddo
        enddo
      enddo

      write(6,'(es15.7)') hunit
      write(6,'(3es15.7)') ((ht(ia,ib,0)/hunit,ia=1,3),ib=1,3)
      write(6,'(i10,3i4)') natmt, 8, 0, 0
      do i=1,natmt
        write(6,'(i3,3es12.4,10es11.3)') int(tagt(i)),rat(1:3,i)
     &       ,ekit(i),epit(i)
     &       ,strst(1,1,i),strst(2,2,i),strst(3,3,i)
     &       ,strst(2,3,i),strst(3,1,i),strst(1,2,i)
c     &       ,eki(i),epi(i),va(1:3,i)
      enddo

      end program qmcl2akr
c-----for emacs---------------------------------------------------------
c     Local Variables:
c     compile-command: "make qmcl2akr"
c     End:
