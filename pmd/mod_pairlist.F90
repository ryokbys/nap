module pairlist
!-----------------------------------------------------------------------
! Module for pair-list.
!-----------------------------------------------------------------------
  
contains
  subroutine init_pairlist(namax)
    
    
  end subroutine init_pairlist
!=======================================================================
  subroutine mk_lspr_para(namax,natm,nbmax,nb,nnmax,tag,ra,rc,rc1nn &
       ,h,hi,anxi,anyi,anzi,lspr,ls1nn,iprint,l1st)
    implicit none
    integer,intent(in):: namax,natm,nbmax,nb,nnmax,iprint
    integer,intent(out):: lspr(0:nnmax,namax),ls1nn(0:nnmax,namax)
    real(8),intent(in):: ra(3,namax),rc,rc1nn,anxi,anyi,anzi &
         ,hi(3,3),h(3,3),tag(namax)
    logical,intent(in):: l1st

    integer:: i,j,k,l,m,n
    integer:: mx,my,mz,kux,kuy,kuz,m1x,m1y,m1z,m1,ic,jc,ierr
    real(8):: xi(3),xij(3),rij(3),rij2

    integer,allocatable,save:: lscl(:),lshd(:)
    real(8),save:: rc2,rcx,rcy,rcz,rcxi,rcyi,rczi,rc1nn2
    integer,save:: lcx,lcy,lcz,lcx2,lcy2,lcz2,lcyz2,lcxyz2

    if( l1st ) then
      rc2= rc**2
      rc1nn2 = rc1nn**2
!-----make a linked cell list, LSCL
      lcx=anxi/dsqrt(hi(1,1)**2+hi(1,2)**2+hi(1,3)**2)/rc
      lcy=anyi/dsqrt(hi(2,1)**2+hi(2,2)**2+hi(2,3)**2)/rc
      lcz=anzi/dsqrt(hi(3,1)**2+hi(3,2)**2+hi(3,3)**2)/rc
!.....In case that system is thinner than rc, modify lc?.
!.....but notice this modification does not correct results.
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
!        write(6,'(a,3i8)') ' lcx,lcy,lcz=',lcx,lcy,lcz
!        write(6,'(a,3i8)') ' lcx2,lcy2,lcz2=',lcx2,lcy2,lcz2
!        write(6,'(a,i8)') ' lcxyz2=',lcxyz2
!        write(6,'(a,3es12.4)') ' rcx,rcy,rcz=',rcx,rcy,rcz
!-----allocate LSCL & LSHD after obtaining lcxyz2
      if( allocated(lscl) ) deallocate(lscl,lshd)
      allocate(lscl(namax+nbmax),lshd(lcxyz2))
    endif

!-----reset pair list, LSPR
    lspr(0,:)= 0
    ls1nn(0,:)= 0

!-----reset headers
    lshd(1:lcxyz2)= 0


!-----construct a linked-cell list, LSCL, & a header list, LSHD
    do i=1,natm+nb
!-------assign a vector cell index
      mx=(ra(1,i)+rcx)*rcxi
      my=(ra(2,i)+rcy)*rcyi
      mz=(ra(3,i)+rcz)*rczi
!-------classify residents in inner cells even if they are not
      if(i.le.natm) then
        mx= min(max(mx,1),lcx)
        my= min(max(my,1),lcy)
        mz= min(max(mz,1),lcz)
!-------copied atoms are either in inner or surface cells
      else
        mx= min(max(mx,0),lcx+1)
        my= min(max(my,0),lcy+1)
        mz= min(max(mz,0),lcz+1)
      endif
      m= mx*lcyz2 +my*lcz2 +mz +1
      lscl(i)= lshd(m)
!-------the last one goes to the header
      lshd(m)= i
    enddo
!      write(6,'(a)') ' lscl,lshd done'

!-----make a pair list, LSPR
!-----Scan resident cells
!!$  do mz=1,lcz
!!$    do my=1,lcy
!!$      do mx=1,lcx
    do mz=0,lcz+1
      do my=0,lcy+1
        do mx=0,lcx+1
          m= mx*lcyz2 +my*lcz2 +mz +1
          if(lshd(m).eq.0) goto 5
          do kuz= -1,1
            m1z= mz +kuz
            if( m1z.lt.0 .or. m1z.gt.lcz+1 ) cycle
            do kuy= -1,1
              m1y= my +kuy
              if( m1y.lt.0 .or. m1y.gt.lcy+1 ) cycle
              do kux= -1,1
                m1x= mx +kux
                if( m1x.lt.0 .or. m1x.gt.lcx+1 ) cycle
                m1=m1x*lcyz2 +m1y*lcz2 +m1z +1
                if(lshd(m1).eq.0) goto 6

                i=lshd(m)
1               continue
!              if (natm.lt.i) goto 4

                ic= int(tag(i))
                xi(1:3)= ra(1:3,i)

                j=lshd(m1)

2               continue
                if( j.le.i ) goto 3
!          if (j.eq.i) goto 3

                jc= int(tag(j))
                xij(1:3)= ra(1:3,j) -xi(1:3)
                rij(1)= h(1,1)*xij(1) +h(1,2)*xij(2) +h(1,3)*xij(3)
                rij(2)= h(2,1)*xij(1) +h(2,2)*xij(2) +h(2,3)*xij(3)
                rij(3)= h(3,1)*xij(1) +h(3,2)*xij(2) +h(3,3)*xij(3)
                rij2= rij(1)**2 +rij(2)**2 +rij(3)**2

                if( rij2.lt.rc2 ) then
                  lspr(0,i)= lspr(0,i) +1
                  if( lspr(0,i).gt.nnmax ) then
                    write(6,'(a)') " ERROR: lspr(0,i)  > nnmax"
                    write(6,'(a,3i5)') "   nnmax, lspr(0,i) = " &
                         ,nnmax,lspr(0,i)
                    write(6,'(a)') " You should rerun pmd with increased nnmax " &
                         //"with the following in.pmd option,"
                    write(6,'(a,i5)') "   max_num_neighbors   ",nnmax+100
                    call mpi_finalize(ierr)
                    stop
                  endif
                  lspr(lspr(0,i),i)=j
!!$!.....Store i in j's neighbor list
!!$                if( j.le.natm ) then
!!$                  lspr(0,j)= lspr(0,j)+1
!!$                  if( lspr(0,j).gt.nnmax ) then
!!$                    write(6,'(a)') " Error: lspr(0,j) > nnmax"
!!$                    write(6,'(a,3i5)') "  nnmax, lspr(0,j) = " &
!!$                         ,nnmax,lspr(0,j)
!!$                    call mpi_finalize(ierr)
!!$                    stop
!!$                  endif
!!$                  lspr(lspr(0,j),j)=i
!!$                endif
                  lspr(0,j)= lspr(0,j)+1
                  if( lspr(0,j).gt.nnmax ) then
                    write(6,'(a)') " ERROR: lspr(0,j) > nnmax"
                    write(6,'(a,3i5)') "   nnmax, lspr(0,j) = " &
                         ,nnmax,lspr(0,j)
                    write(6,'(a)') " You should rerun pmd with increased nnmax " &
                         //"with the following in.pmd option,"
                    write(6,'(a,i5)') "   max_num_neighbors   ",nnmax+100
                    call mpi_finalize(ierr)
                    stop
                  endif
                  lspr(lspr(0,j),j)=i
                  if( rij2.lt.rc1nn2 ) then
                    ls1nn(0,i)= ls1nn(0,i) +1
                    ls1nn(ls1nn(0,i),i)= j
!.....Regardless of j.ne.natm or not, add a counter data in ls1nn(*,j)
                    ls1nn(0,j)= ls1nn(0,j) +1
                    ls1nn(ls1nn(0,j),j)= i
                  endif
                endif

!---------Continue until j= 0
3               j=lscl(j)
                if (j.gt.0) goto 2

!---------Continue until i= 0
4               i=lscl(i)
                if (i.gt.0) goto 1

6               continue
              enddo
            enddo
          enddo
5         continue
        enddo
      enddo
    enddo

  end subroutine mk_lspr_para
!=======================================================================
  subroutine mk_lspr_sngl(namax,natm,nnmax,tag,ra,rc,rc1nn,h,hi &
       ,lspr,ls1nn,iprint,l1st)
!
! Make lspr in serial implimentation taking the periodic boundary
! condition into account.
!
    implicit none
    integer,intent(in):: namax,natm,nnmax,iprint
    integer,intent(out):: lspr(0:nnmax,namax),ls1nn(0:nnmax,namax)
    real(8),intent(in):: ra(3,namax),rc,rc1nn &
         ,hi(3,3),h(3,3),tag(namax)
    logical,intent(in):: l1st

    integer:: i,j,k,l,m,n
    integer:: mx,my,mz,kux,kuy,kuz,m1x,m1y,m1z,m1,ic,jc,ierr
    real(8):: xi(3),xij(3),rij(3),rij2

    integer,allocatable,save:: lscl(:),lshd(:)
    real(8),save:: rc2,rcx,rcy,rcz,rcxi,rcyi,rczi,rc1nn2
    integer,save:: lcx,lcy,lcz,lcyz,lcxyz

    if( l1st ) then
      rc2= rc**2
      rc1nn2= rc1nn**2
!-----make a linked cell list, LSCL
      lcx= 1d0/dsqrt(hi(1,1)**2+hi(1,2)**2+hi(1,3)**2)/rc
      lcy= 1d0/dsqrt(hi(2,1)**2+hi(2,2)**2+hi(2,3)**2)/rc
      lcz= 1d0/dsqrt(hi(3,1)**2+hi(3,2)**2+hi(3,3)**2)/rc
      if( lcx.lt.2 .or. lcy.lt.2 .or. lcz.lt.2 ) then
        write(6,'(a)') ' [Error] mk_lspr_sngl cannot handle' &
             //' too small system !!!'
        stop
      endif
      lcyz= lcy*lcz
      lcxyz= lcx*lcyz
      rcx= 1d0/lcx
      rcy= 1d0/lcy
      rcz= 1d0/lcz
      rcxi=1d0/rcx
      rcyi=1d0/rcy
      rczi=1d0/rcz
!-----allocate LSCL & LSHD after obtaining lcxyz
      if( allocated(lscl) ) deallocate(lscl,lshd)
      allocate(lscl(namax),lshd(lcxyz))
    endif

!-----reset pair list, LSPR
    lspr(0,1:natm)= 0
    ls1nn(0,1:natm)= 0

!-----reset headers
    lshd(1:lcxyz)= 0


!-----construct a linked-cell list, LSCL, & a header list, LSHD
    do i=1,natm
!-------assign a vector cell index
      mx=(ra(1,i)+rcx)*rcxi
      my=(ra(2,i)+rcy)*rcyi
      mz=(ra(3,i)+rcz)*rczi
      mx= min(max(mx,1),lcx)
      my= min(max(my,1),lcy)
      mz= min(max(mz,1),lcz)
      m= (mx-1)*lcyz +(my-1)*lcz +mz
      lscl(i)= lshd(m)
!-------the last one goes to the header
      lshd(m)= i
    enddo
!      write(6,'(a)') ' lscl,lshd done'

!-----make a pair list, LSPR
!-----Scan resident cells
    do mz=1,lcz
      do my=1,lcy
        do mx=1,lcx
          m= (mx-1)*lcyz +(my-1)*lcz +mz
          if (lshd(m).eq.0) goto 5
          do kuz= -1,1
            do kuy= -1,1
              do kux= -1,1
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
                if (lshd(m1).eq.0) goto 6

                i=lshd(m)
1               continue
                if (natm.lt.i) goto 4

                ic= int(tag(i))
                xi(1:3)= ra(1:3,i)

                j=lshd(m1)

2               continue
!          if (j.eq.i) goto 3
                if( j.le.i ) goto 3
                jc= int(tag(j))
                xij(1:3)= ra(1:3,j)-xi(1:3) -anint(ra(1:3,j)-xi(1:3))
                rij(1:3)=h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
                rij2= rij(1)**2 +rij(2)**2 +rij(3)**2

                if(rij2.lt.rc2) then
                  do l=1,lspr(0,i)
                    if( lspr(l,i).eq.j ) then
                      write(6,'(a)') " [Error] lspr(0,i) already has it !!!"
                      write(6,'(a)') "   Check the simulation cell size."
                      stop
                    endif
                  enddo
                  lspr(0,i)= lspr(0,i) +1
                  if(lspr(0,i).gt.nnmax ) then
                    write(6,'(a)') " [Error] lspr(0,i) > nnmax"
                    print *, '   i,nnmax,lspr(0,i) =' &
                         ,i,nnmax,lspr(0,i)
                    stop
                  endif
                  lspr(lspr(0,i),i)=j
!.....Store i in j's neighbor list
                  lspr(0,j)= lspr(0,j)+1
                  if( lspr(0,j).gt.nnmax ) then
                    write(6,'(a)') " Error: lspr(0,j) > nnmax"
                    write(6,'(a,3i5)') "  nnmax, lspr(0,j) = " &
                         ,nnmax,lspr(0,j)
                    call mpi_finalize(ierr)
                    stop
                  endif
                  lspr(lspr(0,j),j)=i
!.....1st NN neighbors
                  if( rij2.lt.rc1nn2 ) then
                    ls1nn(0,i)= ls1nn(0,i) +1
                    ls1nn(ls1nn(0,i),i)= j
                    ls1nn(0,j)= ls1nn(0,j) +1
                    ls1nn(ls1nn(0,j),j)= i
                  endif
                endif

!---------Continue until j= 0
3               j=lscl(j)
                if (j.gt.0) goto 2

!---------Continue until i= 0
4               i=lscl(i)
                if (i.gt.0) goto 1

6               continue
              enddo
            enddo
          enddo
5         continue
        enddo
      enddo
    enddo

  end subroutine mk_lspr_sngl
!=======================================================================
  subroutine mk_lspr_brute(namax,natm,nbmax,nb,nnmax,tag,ra,rc &
       ,rc1nn,h,hi,sgm,lspr,ls1nn,iprint,l1st)
!
!  Make pair list, lspr, by brute force approach, because the system
!  is supposed to be small. Expand the system to take all the atoms 
!  within given cutoff radius into account.
!
    implicit none
    integer,intent(in):: namax,natm,nnmax,nbmax,iprint
    integer,intent(out):: lspr(0:nnmax,natm),nb,ls1nn(0:nnmax,natm)
    real(8),intent(in):: rc,rc1nn,hi(3,3),h(3,3),sgm(3,3)
    real(8),intent(out):: ra(3,namax),tag(namax)
    logical,intent(in):: l1st

    integer:: i,j,k,l,m,n,ia,ja,inc,ix,iy,iz
    real(8):: tmp,xi(3),sij(3),xij(3),rij,vol,asgm

    integer,save:: naex,nex(3)
    real(8),save:: rc2,rc1nn2

    if( l1st ) then
      rc2= rc**2
      rc1nn2= rc1nn**2
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
!          nex(i)= int(1d0/dsqrt(hi(1,i)**2+hi(2,i)**2+hi(3,i)**2)/rc)
!          if( nex(i).lt.1 ) nex(i)= 1
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


!.....Expand system to take into account the image atoms.
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

!.....Search neighbor atoms and make the list
    lspr(0:nnmax,1:natm)= 0
    ls1nn(0:nnmax,1:natm)= 0
    do ia=1,natm
      xi(1:3)= ra(1:3,ia)
      do ja=1,naex
!          if( ja.eq.ia ) cycle
        if( ja.le.ia ) cycle
        sij(1:3)= ra(1:3,ja)-xi(1:3)
        xij(1:3)= h(1:3,1)*sij(1) +h(1:3,2)*sij(2) +h(1:3,3)*sij(3)
        rij= xij(1)**2 +xij(2)**2 +xij(3)**2
        if( rij.lt.rc2 ) then
          lspr(0,ia)= lspr(0,ia) +1
          if( lspr(0,ia).gt.nnmax ) then
            write(6,'(a)') ' [Error] lspr(0,ia) > nnmax'
            write(6,'(a,3i4)') ' nnmax,lspr(0,ia)=' &
                 ,nnmax,lspr(0,ia)
            stop
          endif
          lspr(lspr(0,ia),ia)= ja
          rij = dsqrt(rij)
!.....Store i in j's neighbor list
          lspr(0,ja)= lspr(0,ja)+1
          if( lspr(0,ja).gt.nnmax ) then
            write(6,'(a)') " Error: lspr(0,ja) > nnmax"
            write(6,'(a,3i5)') "  nnmax, lspr(0,ja) = " &
                 ,nnmax,lspr(0,ja)
            stop
          endif
          lspr(lspr(0,ja),ja)=ia
          if( rij.lt.rc1nn ) then
            ls1nn(0,ia)= ls1nn(0,ia) +1
            ls1nn(ls1nn(0,ia),ia)= ja
          endif
        endif
      enddo
    enddo

!      if( l1st ) then
!        do ia=1,natm
!          write(6,'(i3,a,50i6)') ia,": ",(lspr(i,ia),i=1,lspr(0,ia))
!        enddo
!      endif

  end subroutine mk_lspr_brute
!=======================================================================
!!$  subroutine sort_by_lscl(namax,natm,nbmax,nb,tag,ra,va)
!!$!
!!$! To gather atom data more accessible in memory space,
!!$! make atoms in a cell contiguous in memory.
!!$!
!!$  
!!$    
!!$  end subroutine sort_by_lscl
end module pairlist
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
