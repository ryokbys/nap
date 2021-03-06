module pairlist
!-----------------------------------------------------------------------
! Module for pair-list.
!-----------------------------------------------------------------------
  use memory, only: accum_mem
  implicit none
  save
  
  integer,allocatable:: lscl(:),lshd(:)
  real(8):: rc2,rcx,rcy,rcz,rcxi,rcyi,rczi
  integer:: lcx,lcy,lcz,lcxyz,lcyz,lcx2,lcy2,lcz2,lcyz2,lcxyz2
  real(8),allocatable:: tmparr(:)
  integer:: ndmax

!.....Arrays for Gonnet's algorithm
  real(8),allocatable:: dlist(:)
  integer,allocatable:: ilist(:),ileft(:),ia2ic(:)
  real(8):: dh(3,3),vecc(3,26)
  integer:: idcell(-1:1,-1:1,-1:1)
  
contains
!=======================================================================
  subroutine mk_lscl_para(namax,natm,nbmax,nb,ra,anxi,anyi,anzi &
       ,rc,h,hi,l1st)
!
! Make a linked cell list.
! Codes are slightly different bewteen parallel and single.
!
    integer,intent(in):: namax,natm,nbmax,nb
    real(8),intent(in):: anxi,anyi,anzi,rc,h(3,3),hi(3,3)
    real(8),intent(inout):: ra(3,namax)
    logical,intent(in):: l1st

    integer:: i,mx,my,mz,m

    if( l1st ) then
      rc2= rc**2
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
!-----allocate LSCL & LSHD after obtaining lcxyz2
      if( allocated(lscl) ) then
        call accum_mem('pairlist',-4*size(lscl) -4*size(lshd))
        deallocate(lscl,lshd)
      endif
      allocate(lscl(namax+nbmax),lshd(lcxyz2))
      call accum_mem('pairlist',4*size(lscl)+4*size(lshd))
    endif

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
    
    return
  end subroutine mk_lscl_para
!=======================================================================
  subroutine reorder_arrays(namax,natm,nb,tag,ra,va,aux,naux)
!
!  Reorder arrays (tag,ra,va, and some more) using the cell list
!
    integer,intent(in):: namax,natm,nb,naux
    real(8),intent(inout):: tag(namax),ra(3,namax),va(3,namax)
    real(8),intent(inout):: aux(naux,namax)

!.....Sort arrays
    call sort_by_lscl(namax,natm,nb,1,tag)
    call sort_by_lscl(namax,natm,nb,3,ra)
    call sort_by_lscl(namax,natm,nb,3,va)
    call sort_by_lscl(namax,natm,nb,naux,aux)
!!$    if( luse_charge ) then
!!$      call sort_by_lscl(namax,natm,nb,1,chg)
!!$      call sort_by_lscl(namax,natm,nb,1,chi)
!!$    endif
!!$    if( luse_elec_temp ) then
!!$      call sort_by_lscl(namax,natm,nb,1,tei)
!!$    endif
!!$    if( lclrchg ) then
!!$      call sort_by_lscl(namax,natm,nb,1,clr)
!!$    endif
    
    return
  end subroutine reorder_arrays
!=======================================================================
  subroutine mk_lspr_para(namax,natm,nbmax,nb,nnmax,tag,ra,va &
       ,rc,h,hi,anxi,anyi,anzi,lspr,d2lspr,iprint,l1st)
!
!  Make a pairlist using cell list created before.
!
    implicit none
    integer,intent(in):: namax,natm,nbmax,nb,nnmax,iprint
    integer,intent(out):: lspr(0:nnmax,namax)
    real(8),intent(in):: rc,anxi,anyi,anzi,hi(3,3),h(3,3)
    real(8),intent(inout):: ra(3,namax),tag(namax),va(3,namax)
    real(8),intent(out):: d2lspr(nnmax,namax)
    logical,intent(in):: l1st

    integer:: i,j,k,l,m,n,inc,nni,nnj
    integer:: mx,my,mz,kux,kuy,kuz,m1x,m1y,m1z,m1,ic,jc,ierr,mmax
    real(8):: xi(3),xij(3),rij(3),rij2

    call mk_lscl_para(namax,natm,nbmax,nb,ra,anxi,anyi,anzi,rc &
         ,h,hi,l1st)

    if( l1st ) then
      mmax = 0
      do mz=1,lcz
        do my=1,lcy
          do mx=1,lcx
            m= mx*lcyz2 +my*lcz2 +mz +1
            i = lshd(m)
            inc = 0
            do while(i.gt.0)
              inc = inc + 1
              i = lscl(i)
            enddo
            mmax = max(mmax,inc)
          enddo
        enddo
      enddo
!!$!.....If nnmax.lt.(4*pi/3)*mmax, nnmax would be a bit too small
!!$      if( nnmax.lt.4.2*mmax ) then
!!$        write(6,'(a)') " ================= WARNING ======================="
!!$        write(6,'(a)') "   nnmax is less than 5*(num of atoms in a cell)"
!!$        write(6,'(a,i0)') "   You had better set max_num_meighbors greater than " &
!!$             , int(4.2*mmax)
!!$        write(6,'(a)') " ================================================="
!!$      endif
    endif

!-----reset pair list, LSPR
    lspr(0,:)= 0
    d2lspr(:,:) = 0d0

!-----make a pair list, LSPR
!.....Scan atoms, which would be more efficient with OpenMP than scanning cells
!$omp parallel
!$omp do private(mz,my,mx,m,kuz,kuy,kux,m1z,m1y,m1x,m1,i,j,xi,xij,rij,rij2)
    do i=1,natm
      xi(1:3) = ra(1:3,i)
!-------assign a vector cell index
      mx=(xi(1)+rcx)*rcxi
      my=(xi(2)+rcy)*rcyi
      mz=(xi(3)+rcz)*rczi
      mx= min(max(mx,0),lcx+1)
      my= min(max(my,0),lcy+1)
      mz= min(max(mz,0),lcz+1)
      m= mx*lcyz2 +my*lcz2 +mz +1
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
            if(lshd(m1).eq.0) cycle
            j=lshd(m1)
            do while( j.gt.0 )
              if( j.eq.i ) then
                j = lscl(j)
                cycle
              endif
              xij(1:3)= ra(1:3,j) -xi(1:3)
              rij(1)= h(1,1)*xij(1) +h(1,2)*xij(2) +h(1,3)*xij(3)
              rij(2)= h(2,1)*xij(1) +h(2,2)*xij(2) +h(2,3)*xij(3)
              rij(3)= h(3,1)*xij(1) +h(3,2)*xij(2) +h(3,3)*xij(3)
              rij2= rij(1)**2 +rij(2)**2 +rij(3)**2

              if( rij2.lt.rc2 ) then
                lspr(0,i) = lspr(0,i) +1
                lspr(lspr(0,i),i) = j
                d2lspr(lspr(0,i),i) = rij2
              endif

              j=lscl(j)
            enddo ! while (j.gt.0)
          enddo  ! kux
        enddo  ! kuy
      enddo  ! kux
    enddo  ! i=1,natm
!$omp end do
!$omp end parallel

!.....Scan resident cells, which would be inefficient with OpenMP.
!!$    do mz=0,lcz+1
!!$      do my=0,lcy+1
!!$        do mx=0,lcx+1
!!$          m= mx*lcyz2 +my*lcz2 +mz +1
!!$          if(lshd(m).eq.0) cycle
!!$          do kuz= -1,1
!!$            m1z= mz +kuz
!!$            if( m1z.lt.0 .or. m1z.gt.lcz+1 ) cycle
!!$            do kuy= -1,1
!!$              m1y= my +kuy
!!$              if( m1y.lt.0 .or. m1y.gt.lcy+1 ) cycle
!!$              do kux= -1,1
!!$                m1x= mx +kux
!!$                if( m1x.lt.0 .or. m1x.gt.lcx+1 ) cycle
!!$                m1=m1x*lcyz2 +m1y*lcz2 +m1z +1
!!$                if(lshd(m1).eq.0) cycle
!!$
!!$                i=lshd(m)
!!$                do while( i.gt.0 )
!!$                  if( i.gt.natm ) then
!!$                    i = lscl(i)
!!$                    cycle
!!$                  endif
!!$                  ic= int(tag(i))
!!$                  xi(1:3)= ra(1:3,i)
!!$
!!$                  j=lshd(m1)
!!$                  do while( j.gt.0 )
!!$                    if( j.le.i ) then
!!$                      j = lscl(j)
!!$                      cycle
!!$                    endif
!!$
!!$                    jc= int(tag(j))
!!$                    xij(1:3)= ra(1:3,j) -xi(1:3)
!!$                    rij(1)= h(1,1)*xij(1) +h(1,2)*xij(2) +h(1,3)*xij(3)
!!$                    rij(2)= h(2,1)*xij(1) +h(2,2)*xij(2) +h(2,3)*xij(3)
!!$                    rij(3)= h(3,1)*xij(1) +h(3,2)*xij(2) +h(3,3)*xij(3)
!!$                    rij2= rij(1)**2 +rij(2)**2 +rij(3)**2
!!$
!!$                    if( rij2.lt.rc2 ) then
!!$                      lspr(0,i)= lspr(0,i) +1
!!$                      lspr(0,j)= lspr(0,j) +1
!!$                      lspr(lspr(0,i),i)=j
!!$                      lspr(lspr(0,j),j)=i
!!$                      d2lspr(lspr(0,i),i) = rij2
!!$                      d2lspr(lspr(0,j),j) = rij2
!!$                    endif
!!$
!!$                    j=lscl(j)
!!$                  enddo
!!$
!!$                  i=lscl(i)
!!$                enddo
!!$
!!$              enddo
!!$            enddo
!!$          enddo
!!$        enddo
!!$      enddo
!!$    enddo

!!$!.....Only 1st call
!!$    if( l1st ) then
!!$      mmax = 0
!!$      do i=1,natm
!!$        mmax = max(mmax,lspr(0,i))
!!$      enddo
!!$      print '(a,i0)',' Max num of neighbors at 1st call = ',mmax
!!$    endif

  end subroutine mk_lspr_para
!=======================================================================
  subroutine mk_lspr_gonnet(namax,natm,nbmax,nb,nnmax,tag,ra,va &
       ,rc,h,hi,anxi,anyi,anzi,lspr,d2lspr,iprint,l1st)
!
!  Make a pairlist by Gonnet's algorithm using cell list created before.
!  Ref.
!    1) P. Gonnet, J. Comput. Chem. 28, 570–573 (2007)
!
    implicit none
    integer,intent(in):: namax,natm,nbmax,nb,nnmax,iprint
    integer,intent(out):: lspr(0:nnmax,namax)
    real(8),intent(in):: rc,anxi,anyi,anzi,hi(3,3),h(3,3)
    real(8),intent(inout):: ra(3,namax),tag(namax),va(3,namax)
    real(8),intent(out):: d2lspr(nnmax,namax)
    logical,intent(in):: l1st

    integer:: i,j,k,l,m,n,ii,jj,inc,nleft
    integer:: mx,my,mz,kux,kuy,kuz,m1x,m1y,m1z,m1,ic,jc,ierr,mmax
    real(8):: xi(3),xj(3),xij(3),rij(3),rij2,vc(3)
    real(8):: tmp

    if( .not. allocated(dlist) ) then
      allocate(dlist(namax),ilist(namax),ileft(namax),ia2ic(namax))
      call accum_mem('pairlist',8*size(dlist)+4*size(ilist) &
           +4*size(ileft)+4*size(ia2ic))
    else if( size(dlist).ne.namax ) then
      call accum_mem('pairlist',-8*size(dlist)-4*size(ilist) &
           -4*size(ileft)-4*size(ia2ic))
      deallocate(dlist,ilist,ileft,ia2ic)
      allocate(dlist(namax),ilist(namax),ileft(namax),ia2ic(namax))
      call accum_mem('pairlist',8*size(dlist)+4*size(ilist) &
           +4*size(ileft)+4*size(ia2ic))
    endif

    if( l1st ) then
!.....Make dh matrix
      dh(1:3,1) = h(1:3,1)*anxi
      dh(1:3,2) = h(1:3,2)*anyi
      dh(1:3,3) = h(1:3,3)*anzi
!.....Make vectors of neighboring cells
      idcell(-1:1,-1:1,-1:1) = 0
      inc = 0
      do kuz= -1,1
        do kuy= -1,1
          do kux= -1,1
            if( kuz.eq.0 .and. kuy.eq.0 .and. kux.eq.0 ) cycle
            inc = inc +1
            vecc(1:3,inc) = dh(1:3,1)*kux &
                 +dh(1:3,2)*kuy +dh(1:3,3)*kuz
            rij2 = vecc(1,inc)**2 +vecc(2,inc)**2 +vecc(3,inc)**2
            vecc(1:3,inc) = vecc(1:3,inc)/sqrt(rij2)
            idcell(kux,kuy,kuz) = inc
          enddo
        enddo
      enddo
    endif

    call mk_lscl_para(namax,natm,nbmax,nb,ra,anxi,anyi,anzi,rc &
         ,h,hi,l1st)

!.....Make a converter that converts atom index to cell index
    ia2ic(:) = 0
    do mz=0,lcz+1
      do my=0,lcy+1
        do mx=0,lcx+1
          m= mx*lcyz2 +my*lcz2 +mz +1
          if(lshd(m).eq.0) cycle
          i = lshd(m)
          ia2ic(i) = m
          do while(lscl(i).ne.0)
            i = lscl(i)
            ia2ic(i) = m
          enddo
        enddo
      enddo
    enddo

    dlist(:) = 0d0
    ilist(:) = 0

!-----reset pair list, LSPR
    lspr(0,:)= 0
    d2lspr(:,:) = 0d0

!=====Make a pair list, LSPR, using Gonnet's algorithm=====
!.....Scan resident cells
    do mz=0,lcz+1
      do my=0,lcy+1
        do mx=0,lcx+1
          m= mx*lcyz2 +my*lcz2 +mz +1
          if(lshd(m).eq.0) cycle
          i = lshd(m)
11        continue
          ic= int(tag(i))
          xi(1:3)= ra(1:3,i)

!.....Search for neighbors within the same cell
          j=lshd(m)
12        continue
          if( j.le.i ) goto 13

          jc= int(tag(j))
          xij(1:3)= ra(1:3,j) -xi(1:3)
          rij(1)= h(1,1)*xij(1) +h(1,2)*xij(2) +h(1,3)*xij(3)
          rij(2)= h(2,1)*xij(1) +h(2,2)*xij(2) +h(2,3)*xij(3)
          rij(3)= h(3,1)*xij(1) +h(3,2)*xij(2) +h(3,3)*xij(3)
          rij2= rij(1)**2 +rij(2)**2 +rij(3)**2

          if( rij2.lt.rc2 ) then
            lspr(0,i)= lspr(0,i) +1
            lspr(0,j)= lspr(0,j) +1
            if( lspr(0,i).gt.nnmax .or. lspr(0,j).gt.nnmax ) then
              write(6,'(a)') " ERROR: lspr(0,i or j)  > nnmax"
              write(6,'(a,3i5)') "   nnmax,lspr(0,i),lspr(0,j) = " &
                   ,nnmax,lspr(0,i),lspr(0,j)
              write(6,'(a)') " You should rerun pmd with increased nnmax " &
                   //"with the following in.pmd option,"
              write(6,'(a,i5)') "   max_num_neighbors   ",nnmax+100
              stop
            endif
            lspr(lspr(0,i),i)=j
            lspr(lspr(0,j),j)=i
            d2lspr(lspr(0,i),i) = rij2
            d2lspr(lspr(0,j),j) = rij2
          endif

!.....Continue until j= 0
13        j=lscl(j)
          if (j.gt.0) goto 12

!.....Continue until i= 0
14        i=lscl(i)
          if (i.gt.0) goto 11

!.....Search for neighbors from neighboring cells,
!     where Gonnet's algorithm is used.
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
                if( m1.eq.m ) cycle
                if( lshd(m1).eq.0 ) cycle

!.....Vector towards the neighbor cell
                ic = idcell(kux,kuy,kuz)
                vc(1:3) = vecc(1:3,ic)
!.....Project atoms in cell-m onto the neighboring cell vector
                dlist(:) = 0d0
                ilist(:) = 0
                inc = 0
                i= lshd(m)
                do while( i.ne.0 )
                  xi(1:3) = h(1:3,1)*ra(1,i) &
                     +h(1:3,2)*ra(2,i) +h(1:3,3)*ra(3,i)
                  inc = inc + 1
                  dlist(inc) = vc(1)*xi(1) +vc(2)*xi(2) +vc(3)*xi(3)
                  ilist(inc) = i
                  i = lscl(i)
                enddo
!.....Project atoms in cell-m1 onto the neighboring cell vector 
                j= lshd(m1)
                do while( j.ne.0 )
                  xj(1:3) = h(1:3,1)*ra(1,j) &
                       +h(1:3,2)*ra(2,j) +h(1:3,3)*ra(3,j)
                  inc = inc + 1
                  dlist(inc) = vc(1)*xj(1) +vc(2)*xj(2) +vc(3)*xj(3) -rc
                  ilist(inc) = j
                  j = lscl(j)
                enddo
!.....Sort arrays dlist and ilist according to dlist
                call qsort_list(inc,1,inc,dlist,ilist)
!.....Finally, extract pairs rij<rc to make lspr
                ileft(:) = 0
                nleft = 0
                do ii=1,inc
                  i = ilist(ii)
                  if( ia2ic(i).eq.m1 ) then
                    nleft = nleft +1
                    ileft(nleft) = i
                  else
                    do jj=1,nleft
                      j = ileft(jj)
                      if( j.le.i ) cycle
                      xij(1:3)= ra(1:3,j) -ra(1:3,i)
                      rij(1)= h(1,1)*xij(1) +h(1,2)*xij(2) +h(1,3)*xij(3)
                      rij(2)= h(2,1)*xij(1) +h(2,2)*xij(2) +h(2,3)*xij(3)
                      rij(3)= h(3,1)*xij(1) +h(3,2)*xij(2) +h(3,3)*xij(3)
                      rij2= rij(1)**2 +rij(2)**2 +rij(3)**2
                      if( rij2.lt.rc2 ) then
                        lspr(0,i)= lspr(0,i) +1
                        lspr(0,j)= lspr(0,j) +1
                        lspr(lspr(0,i),i)=j
                        lspr(lspr(0,j),j)=i
                        d2lspr(lspr(0,i),i) = rij2
                        d2lspr(lspr(0,j),j) = rij2
                      endif
                    enddo
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo

!.....Only 1st call
    if( l1st ) then
      mmax = 0
      do i=1,natm
        mmax = max(mmax,lspr(0,i))
      enddo
!!$      print '(a,i0)',' Max num of neighbors at 1st call = ',mmax
    endif

  end subroutine mk_lspr_gonnet
!=======================================================================
  subroutine mk_lspr_sngl(namax,natm,nnmax,tag,ra,rc,h,hi &
       ,lspr,d2lspr,iprint,l1st)
!
! Make lspr in serial implimentation taking the periodic boundary
! condition into account.
!
    implicit none
    integer,intent(in):: namax,natm,nnmax,iprint
    integer,intent(out):: lspr(0:nnmax,namax)
    real(8),intent(in):: ra(3,namax),rc,hi(3,3),h(3,3),tag(namax)
    real(8),intent(out):: d2lspr(nnmax,namax)
    logical,intent(in):: l1st

    integer:: i,j,k,l,m,n
    integer:: mx,my,mz,kux,kuy,kuz,m1x,m1y,m1z,m1,ic,jc,ierr
    real(8):: xi(3),xij(3),rij(3),rij2

!!$    integer,allocatable,save:: lscl(:),lshd(:)
!!$    real(8),save:: rc2,rcx,rcy,rcz,rcxi,rcyi,rczi,rc1nn2
!!$    integer,save:: lcx,lcy,lcz,lcyz,lcxyz

    if( l1st ) then
      rc2= rc**2
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
      if( allocated(lscl) ) then
        call accum_mem('pairlist',-4*size(lscl)-4*size(lshd))
        deallocate(lscl,lshd)
      endif
      allocate(lscl(namax),lshd(lcxyz))
      call accum_mem('pairlist',4*size(lscl)+4*size(lshd))
    endif

!-----reset pair list, LSPR
    lspr(0,1:natm)= 0
    d2lspr(:,:) = 0d0

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
            m1z= mz +kuz
            if( m1z.lt.1   ) m1z= m1z +lcz
            if( m1z.gt.lcz ) m1z= m1z -lcz
            do kuy= -1,1
              m1y= my +kuy
              if( m1y.lt.1   ) m1y= m1y +lcy
              if( m1y.gt.lcy ) m1y= m1y -lcy
              do kux= -1,1
                m1x= mx +kux
                if( m1x.lt.1   ) m1x= m1x +lcx
                if( m1x.gt.lcx ) m1x= m1x -lcx
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
!.....If there is the same index in the neighbor, skip storing this j
                      goto 3
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
                  d2lspr(lspr(0,i),i) = rij2
!.....Store i in j's neighbor list
                  lspr(0,j)= lspr(0,j)+1
                  if( lspr(0,j).gt.nnmax ) then
                    write(6,'(a)') " Error: lspr(0,j) > nnmax"
                    write(6,'(a,3i5)') "  nnmax, lspr(0,j) = " &
                         ,nnmax,lspr(0,j)
                    stop
                  endif
                  lspr(lspr(0,j),j)=i
                  d2lspr(lspr(0,j),j) = rij2
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
       ,h,hi,sgm,lspr,d2lspr,iprint,l1st)
!
!  Make pair list, lspr, by brute force approach, because the system
!  is supposed to be small. Expand the system to take all the atoms 
!  within given cutoff radius into account.
!
    implicit none
    integer,intent(in):: namax,natm,nnmax,nbmax,iprint
    integer,intent(out):: lspr(0:nnmax,natm),nb
    real(8),intent(in):: rc,hi(3,3),h(3,3),sgm(3,3)
    real(8),intent(out):: ra(3,namax),tag(namax),d2lspr(nnmax,namax)
    logical,intent(in):: l1st

    integer:: i,j,k,l,m,n,ia,ja,inc,ix,iy,iz
    real(8):: tmp,xi(3),sij(3),xij(3),rij,vol,asgm,rij2

    integer,save:: naex,nex(3)

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
    d2lspr(:,:) = 0d0
    do ia=1,natm
      xi(1:3)= ra(1:3,ia)
      do ja=1,naex
!          if( ja.eq.ia ) cycle
        if( ja.le.ia ) cycle
        sij(1:3)= ra(1:3,ja)-xi(1:3)
        xij(1:3)= h(1:3,1)*sij(1) +h(1:3,2)*sij(2) +h(1:3,3)*sij(3)
        rij2= xij(1)**2 +xij(2)**2 +xij(3)**2
        if( rij2.lt.rc2 ) then
          lspr(0,ia)= lspr(0,ia) +1
          if( lspr(0,ia).gt.nnmax ) then
            write(6,'(a)') ' [Error] lspr(0,ia) > nnmax'
            write(6,'(a,3i4)') ' nnmax,lspr(0,ia)=' &
                 ,nnmax,lspr(0,ia)
            stop
          endif
          lspr(lspr(0,ia),ia)= ja
          d2lspr(lspr(0,ia),ia) = rij2
!.....Store i in j's neighbor list
          lspr(0,ja)= lspr(0,ja)+1
          if( lspr(0,ja).gt.nnmax ) then
            write(6,'(a)') " Error: lspr(0,ja) > nnmax"
            write(6,'(a,3i5)') "  nnmax, lspr(0,ja) = " &
                 ,nnmax,lspr(0,ja)
            stop
          endif
          lspr(lspr(0,ja),ja)=ia
          d2lspr(lspr(0,ja),ja) = rij2
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
  subroutine update_d2lspr(namax,natm,nnmax,lspr,h,ra,rcut,rbuf,d2lspr)
!
!  Since d2lspr should be updated even if lspr is not changed,
!  update only d2lspr taking the rbuf into account.
!
    integer,intent(in):: namax,natm,nnmax,lspr(0:nnmax,natm)
    real(8),intent(in):: rcut,rbuf,ra(3,namax),h(3,3)
    real(8),intent(inout):: d2lspr(nnmax,natm)

    integer:: i,j,jj
    real(8):: rmin2,xi(3),xij(3),rij(3),dij2
    
    rmin2 = (rcut-rbuf)**2
!$omp parallel
!$omp do private(i,xi,jj,j,xij,rij,dij2)
    do i=1,natm
      xi(1:3) = ra(1:3,i)
      do jj=1,lspr(0,i)
!.....Probably, do not need to consider neighbors close enough compared to rcut
        if( d2lspr(jj,i).lt.rmin2 ) cycle
!.....Only neighbors within a shell around rcut need to be updated 
!.....for the purpose of determining whether or not the neighbors are within cutoff
        j = lspr(jj,i)
        xij(1:3) = ra(1:3,j) -xi(1:3)
        xij(1:3) = xij(1:3) -anint(xij(1:3))
        rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij2 = rij(1)**2 +rij(2)**2 +rij(3)**2
        d2lspr(jj,i) = dij2
      enddo
    enddo
!$omp end do
!$omp end parallel
    return
  end subroutine update_d2lspr
!=======================================================================
  subroutine sort_by_lscl(namax,natm,nb,ndim,arr)
!
! To gather atom data more accessible in memory space,
! make atoms in a cell contiguous in 1D memory.
!
    integer,intent(in):: namax,natm,nb,ndim
    real(8),intent(inout):: arr(ndim,namax)

    integer:: n,nmig,mz,my,mx,i,j,m
    logical:: lmig,lmigx,lmigy,lmigz

    if( .not. allocated(tmparr) ) then
      ndmax = max(3,ndim)
      allocate(tmparr(ndmax*namax))
      call accum_mem('pairlist',8*size(tmparr))
    else if( size(tmparr).lt.ndim*namax ) then
      ndmax = max(ndim,ndmax)
      call accum_mem('pairlist',-8*size(tmparr))
      deallocate(tmparr)
      allocate(tmparr(ndmax*namax))
      call accum_mem('pairlist',8*size(tmparr))
    endif

!  Sort arr taking into account residents only,
!  since reordering migrants after copying migrants causes
!  inconsistency of lsb() array.
!  So this routine should be called before calling bacopy.
    n = 0
    do mz=1,lcz
      do my=1,lcy
        do mx=1,lcx
          m= mx*lcyz2 +my*lcz2 +mz +1
          i = lshd(m)
          if( i.eq.0 ) cycle
          if( i.le.natm ) then
            n = n +1
            tmparr(ndim*(n-1)+1:ndim*n) = arr(1:ndim,i)
          endif
          do while(lscl(i).ne.0)
            i = lscl(i)
            if( i.gt.natm ) cycle
            n = n +1
            tmparr(ndim*(n-1)+1:ndim*n) = arr(1:ndim,i)
          end do
        enddo
      enddo
    enddo
    do i=1,n
      do j=1,ndim
        arr(j,i) = tmparr(ndim*(i-1)+j)
      enddo
    enddo
    
  end subroutine sort_by_lscl
!=======================================================================
  subroutine swap(ndim,i,j,dlist,ilist)
    integer,intent(in):: ndim,i,j
    real(8),intent(inout):: dlist(ndim)
    integer,intent(inout):: ilist(ndim)

    integer:: itmp
    real(8):: tmp

    tmp = dlist(i)
    dlist(i) = dlist(j)
    dlist(j) = tmp

    itmp = ilist(i)
    ilist(i) = ilist(j)
    ilist(j) = itmp
    return
  end subroutine swap
!=======================================================================
  recursive subroutine qsort_list(ndim,il,ir,dlist,ilist)
!
!  Sort dlist and ilist used in Gonnet's algorithm by Quicksort.
!
    integer,intent(in):: ndim,il,ir
    real(8),intent(inout):: dlist(ndim)
    integer,intent(inout):: ilist(ndim)

    integer:: ip,i,j
    real(8):: dip

    if( ir-il.lt.1 ) return
    ip = int((il+ir)/2)
    dip = dlist(ip)
    call swap(ndim,ip,ir,dlist,ilist)
    i = il
    do j=il,ir-1
      if( dlist(j).lt.dip ) then
        call swap(ndim,i,j,dlist,ilist)
        i = i + 1
      endif
    enddo
    call swap(ndim,i,ir,dlist,ilist)
    call qsort_list(ndim,il,i,dlist,ilist)
    call qsort_list(ndim,i+1,ir,dlist,ilist)
    
  end subroutine qsort_list
!=======================================================================
end module pairlist
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd lib"
!     End:
