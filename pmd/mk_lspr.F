c-----------------------------------------------------------------------
c Routines for making pair-list.
c-----------------------------------------------------------------------
      subroutine mk_lspr_para(namax,natm,nbmax,nb,nnmax,tag,ra,rc
     &     ,h,hi,anxi,anyi,anzi,lspr)
      implicit none
      integer,intent(in):: namax,natm,nbmax,nb,nnmax
      integer,intent(out):: lspr(0:nnmax,natm)
      real(8),intent(in):: ra(3,namax),rc,anxi,anyi,anzi
     &     ,hi(3,3),h(3,3),tag(namax)

      integer:: i,j,k,l,m,n
      integer:: mx,my,mz,kux,kuy,kuz,m1x,m1y,m1z,m1,ic,jc,ierr
      real(8):: xi(3),xij(3),rij(3),rij2
      
      logical,save:: l1st=.true.
      integer,allocatable,save:: lscl(:),lshd(:)
      real(8),save:: rc2,rcx,rcy,rcz,rcxi,rcyi,rczi
      integer,save:: lcx,lcy,lcz,lcx2,lcy2,lcz2,lcyz2,lcxyz2

      if( l1st ) then
        rc2= rc**2
c-----make a linked cell list, LSCL
        lcx=anxi/dsqrt(hi(1,1)**2+hi(1,2)**2+hi(1,3)**2)/rc
        lcy=anyi/dsqrt(hi(2,1)**2+hi(2,2)**2+hi(2,3)**2)/rc
        lcz=anzi/dsqrt(hi(3,1)**2+hi(3,2)**2+hi(3,3)**2)/rc
c.....In case that system is thinner than rc, modify lc?.
c.....but notice this modification does not correct results.
        if( lcx.eq.0 ) lcx=1
        if( lcy.eq.0 ) lcy=1
        if( lcz.eq.0 ) lcz=1
        lcx2= lcx +2
        lcy2= lcy +2
        lcz2= lcz +2
        lcyz2=lcy2*lcz2
        lcxyz2=lcx2*lcyz2
        rcx= anxi/lcx
        rcy= anyi/lcy
        rcz= anzi/lcz
        rcxi=1d0/rcx
        rcyi=1d0/rcy
        rczi=1d0/rcz
c        write(6,'(a,3i8)') ' lcx,lcy,lcz=',lcx,lcy,lcz
c        write(6,'(a,3i8)') ' lcx2,lcy2,lcz2=',lcx2,lcy2,lcz2
c        write(6,'(a,i8)') ' lcxyz2=',lcxyz2
c        write(6,'(a,3es12.4)') ' rcx,rcy,rcz=',rcx,rcy,rcz
c-----allocate LSCL & LSHD after obtaining lcxyz2
        allocate(lscl(namax+nbmax),lshd(lcxyz2))
        l1st=.false.
      endif

c-----reset pair list, LSPR
      lspr(0,1:natm)= 0

c-----reset headers
      lshd(1:lcxyz2)= 0


c-----construct a linked-cell list, LSCL, & a header list, LSHD
      do i=1,natm+nb
c-------assign a vector cell index
        mx=(ra(1,i)+rcx)*rcxi
        my=(ra(2,i)+rcy)*rcyi
        mz=(ra(3,i)+rcz)*rczi
c-------classify residents in inner cells even if they are not
        if(i.le.natm) then
          mx= min(max(mx,1),lcx)
          my= min(max(my,1),lcy)
          mz= min(max(mz,1),lcz)
c-------copied atoms are either in inner or surface cells
        else
          mx= min(max(mx,0),lcx+1)
          my= min(max(my,0),lcy+1)
          mz= min(max(mz,0),lcz+1)
        endif
        m= mx*lcyz2 +my*lcz2 +mz +1
        lscl(i)= lshd(m)
c-------the last one goes to the header
        lshd(m)= i
      enddo
c      write(6,'(a)') ' lscl,lshd done'

c-----make a pair list, LSPR
c-----Scan resident cells
      do 100 mz=1,lcz
      do 100 my=1,lcy
      do 100 mx=1,lcx
        m= mx*lcyz2 +my*lcz2 +mz +1
        if (lshd(m).eq.0) goto 100
        do 102 kuz= -1,1
        do 102 kuy= -1,1
        do 102 kux= -1,1
          m1x= mx +kux
          m1y= my +kuy
          m1z= mz +kuz
          m1=m1x*lcyz2 +m1y*lcz2 +m1z +1
          if (lshd(m1).eq.0) goto 102

          i=lshd(m)
1         continue
          if (natm.lt.i) goto 4

c          ic=is(i)
          ic= int(tag(i))
          xi(1:3)= ra(1:3,i)

          j=lshd(m1)

2         continue
          if (j.eq.i) goto 3
c          jc=is(j)
          jc= int(tag(i))
          xij(1:3)= ra(1:3,j) -xi(1:3)
          rij(1)= h(1,1)*xij(1) +h(1,2)*xij(2) +h(1,3)*xij(3)
          rij(2)= h(2,1)*xij(1) +h(2,2)*xij(2) +h(2,3)*xij(3)
          rij(3)= h(3,1)*xij(1) +h(3,2)*xij(2) +h(3,3)*xij(3)
          rij2= rij(1)**2 +rij(2)**2 +rij(3)**2

c          if (rij2.lt.rcij2(ic,jc)) then
          if(rij2.lt.rc2) then
            lspr(0,i)=lspr(0,i)+1
            if(lspr(0,i).gt.nnmax) then
              write(6,'(a)') " !!!lspr(0,i).gt.nnmax"
              call mpi_finalize(ierr)
              stop
            endif
            lspr(lspr(0,i),i)=j
          endif

c---------Continue until j= 0
3         j=lscl(j)
          if (j.gt.0) goto 2

c---------Continue until i= 0
4         i=lscl(i)
          if (i.gt.0) goto 1

102     continue
100   continue

      end subroutine mk_lspr_para
c=======================================================================
      subroutine mk_lspr_sngl(namax,natm,nnmax,tag,ra,rc,h,hi,lspr)
c
c Make lspr in serial implimentation taking the periodic boundary
c condition into account.
c
      implicit none
      integer,intent(in):: namax,natm,nnmax
      integer,intent(out):: lspr(0:nnmax,natm)
      real(8),intent(in):: ra(3,namax),rc,hi(3,3),h(3,3),tag(namax)

      integer:: i,j,k,l,m,n
      integer:: mx,my,mz,kux,kuy,kuz,m1x,m1y,m1z,m1,ic,jc,ierr
      real(8):: xi(3),xij(3),rij(3),rij2
      
      logical,save:: l1st=.true.
      integer,allocatable,save:: lscl(:),lshd(:)
      real(8),save:: rc2,rcx,rcy,rcz,rcxi,rcyi,rczi
      integer,save:: lcx,lcy,lcz,lcyz,lcxyz

      if( l1st ) then
        rc2= rc**2
c-----make a linked cell list, LSCL
        lcx= 1d0/dsqrt(hi(1,1)**2+hi(1,2)**2+hi(1,3)**2)/rc
        lcy= 1d0/dsqrt(hi(2,1)**2+hi(2,2)**2+hi(2,3)**2)/rc
        lcz= 1d0/dsqrt(hi(3,1)**2+hi(3,2)**2+hi(3,3)**2)/rc
        if( lcx.lt.2 .or. lcy.lt.2 .or. lcz.lt.2 ) then
          write(6,'(a)') ' [Error] mk_lspr_sngl cannot handle'
     &         //' too small system !!!'
          stop
        endif
        lcyz= lcy*lcx
        lcxyz= lcx*lcyz
        rcx= 1d0/lcx
        rcy= 1d0/lcy
        rcz= 1d0/lcz
        rcxi=1d0/rcx
        rcyi=1d0/rcy
        rczi=1d0/rcz
c-----allocate LSCL & LSHD after obtaining lcxyz
        allocate(lscl(namax),lshd(lcxyz))
        l1st=.false.
      endif

c-----reset pair list, LSPR
      lspr(0,1:natm)= 0

c-----reset headers
      lshd(1:lcxyz)= 0


c-----construct a linked-cell list, LSCL, & a header list, LSHD
      do i=1,natm
c-------assign a vector cell index
        mx=(ra(1,i)+rcx)*rcxi
        my=(ra(2,i)+rcy)*rcyi
        mz=(ra(3,i)+rcz)*rczi
        mx= min(max(mx,1),lcx)
        my= min(max(my,1),lcy)
        mz= min(max(mz,1),lcz)
        m= (mx-1)*lcyz +(my-1)*lcz +mz
        lscl(i)= lshd(m)
c-------the last one goes to the header
        lshd(m)= i
      enddo
c      write(6,'(a)') ' lscl,lshd done'

c-----make a pair list, LSPR
c-----Scan resident cells
      do 100 mz=1,lcz
      do 100 my=1,lcy
      do 100 mx=1,lcx
        m= (mx-1)*lcyz +(my-1)*lcz +mz
        if (lshd(m).eq.0) goto 100
        do 102 kuz= -1,1
        do 102 kuy= -1,1
        do 102 kux= -1,1
          m1x= mx +kux
          m1y= my +kuy
          m1z= mz +kuz
          if( m1x.lt.1   ) m1x= m1x +lcx
          if( m1x.gt.lcx ) m1x= m1x -lcx
          if( m1y.lt.1   ) m1y= m1y +lcy
          if( m1y.gt.lcy ) m1y= m1y -lcy
          if( m1z.lt.1   ) m1z= m1z +lcz
          if( m1z.gt.lcz ) m1z= m1z -lcz
          m1=(m1x-1)*lcyz +(m1y-1)*lcz +m1z
          if (lshd(m1).eq.0) goto 102

          i=lshd(m)
1         continue
          if (natm.lt.i) goto 4

          ic= int(tag(i))
          xi(1:3)= ra(1:3,i)

          j=lshd(m1)

2         continue
          if (j.eq.i) goto 3
          jc= int(tag(i))
          xij(1:3)= ra(1:3,j)-xi(1:3) -anint(ra(1:3,j)-xi(1:3))
          rij(1:3)=h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
          rij2= rij(1)**2 +rij(2)**2 +rij(3)**2

          if(rij2.lt.rc2) then
            lspr(0,i)=lspr(0,i)+1
            if(lspr(0,i).gt.nnmax) then
              write(6,'(a)') " [Error] lspr(0,i).gt.nnmax !!!"
              stop
            endif
            lspr(lspr(0,i),i)=j
          endif

c---------Continue until j= 0
3         j=lscl(j)
          if (j.gt.0) goto 2

c---------Continue until i= 0
4         i=lscl(i)
          if (i.gt.0) goto 1

102     continue
100   continue

      end subroutine mk_lspr_sngl
c=======================================================================
      subroutine mk_lspr_brute(namax,natm,nbmax,nb,nnmax,tag,ra,rc
     &     ,h,hi,sgm,lspr)
c
c  Make pair list, lspr, by brute force approach, because the system
c  is supposed to be small. Expand the system to take all the atoms 
c  within given cutoff radius into account.
c
      implicit none
      integer,intent(in):: namax,natm,nnmax,nbmax
      integer,intent(out):: lspr(0:nnmax,natm),nb
      real(8),intent(in):: rc,hi(3,3),h(3,3),sgm(3,3)
      real(8),intent(out):: ra(3,namax),tag(namax)

      integer:: i,j,k,l,m,n,ia,ja,inc,ix,iy,iz
      real(8):: tmp,xi(3),sij(3),xij(3),rij,vol,asgm

      logical,save:: l1st=.true.
      integer,save:: naex,nex(3)
      real(8),save:: rc2

      if( l1st ) then
        rc2= rc**2
        print *, ' rc,rc2=',rc,rc2
        print *, ' h(1:3,1:3): '
        print *, h(1:3,1)
        print *, h(1:3,2)
        print *, h(1:3,3)
        vol= h(1,1)*sgm(1,1) +h(2,1)*sgm(2,1) +h(3,1)*sgm(3,1)
        do i=1,3
          asgm= dsqrt(sgm(1,i)**2 +sgm(2,i)**2 +sgm(3,i)**2)
          nex(i)= int(rc/(vol/asgm))
          print *, ' rc,vol,asgm,nex=',rc,vol,asgm,nex(i)
c          nex(i)= int(1d0/dsqrt(hi(1,i)**2+hi(2,i)**2+hi(3,i)**2)/rc)
c          if( nex(i).lt.1 ) nex(i)= 1
          nex(i)= nex(i) +1
        enddo
        naex= natm *(2*nex(1)+1) *(2*nex(2)+1) *(2*nex(3)+1)
        print *, ' nex(1:3),naex=',nex(1:3),naex
        if( naex.gt.namax ) then
          write(6,'(a)') ' [Error] naex.gt.namax !!!'
          write(6,'(a,4i6)') ' nex(1:3),naex =',nex(1:3),naex
          print *, ' naex =',naex
          print *, ' namax=',namax
          stop
        endif
      endif


c.....Expand system to take into account the image atoms.
      inc= 0
      do iz=-nex(3),nex(3)
      do iy=-nex(2),nex(2)
      do ix=-nex(1),nex(1)
        if( ix.eq.0 .and. iy.eq.0 .and. iz.eq.0 ) cycle
        do ia=1,natm
          inc= inc +1
          ra(1,natm+inc)= ra(1,ia) +ix
          ra(2,natm+inc)= ra(2,ia) +iy
          ra(3,natm+inc)= ra(3,ia) +iz
          tag(natm+inc)= tag(ia)
        enddo
      enddo
      enddo
      enddo
      if( natm+inc.ne.naex ) then
        print *, '[Error] natm+inc .ne. naex !!!'
        print *, '   natm,inc,naex=',natm,inc,naex
        stop 
      endif
      if( inc.gt.nbmax ) then
        print *, '[Error] inc.gt.nbmax !!!'
        print *, '  inc,nbmax =',inc,nbmax
        stop 
      endif
      nb= inc

c.....Search neighbor atoms and make the list
      lspr(0:nnmax,1:natm)= 0
      do ia=1,natm
        xi(1:3)= ra(1:3,ia)
        do ja=1,naex
          if( ja.eq.ia ) cycle
          sij(1:3)= ra(1:3,ja)-xi(1:3)
          xij(1:3)= h(1:3,1)*sij(1) +h(1:3,2)*sij(2) +h(1:3,3)*sij(3)
          rij= xij(1)**2 +xij(2)**2 +xij(3)**2
          if( rij.lt.rc2 ) then
            lspr(0,ia)= lspr(0,ia) +1
            if( lspr(0,ia).gt.nnmax ) then
              write(6,'(a)') ' [Error] lspr(0,ia).gt.nnmax !!!'
              stop
            endif
            lspr(lspr(0,ia),ia)= ja
          endif
        enddo
      enddo

c      if( l1st ) then
c        do ia=1,natm
c          write(6,'(i3,a,50i6)') ia,": ",(lspr(i,ia),i=1,lspr(0,ia))
c        enddo
c      endif

      l1st=.false.
      end subroutine mk_lspr_brute
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make pmd"
c     End:
