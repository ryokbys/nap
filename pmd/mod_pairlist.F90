module pairlist
!-----------------------------------------------------------------------
! Module for pair-list.
!-----------------------------------------------------------------------
  use memory, only: accum_mem
  implicit none
  include './const.h'
  save
  private

  public:: mk_lspr_para,mk_lscl_para,sort_arrays, &
       check_lspr, check_lscl
  public:: mk_lspr_gonnet, mk_lspr_sngl, mk_lspr_brute, sort_by_lscl, &
       swap, qsort_list, sort_lspr, set_nnmax
  
  integer,allocatable:: lscl(:),lshd(:)
!$acc declare create(lscl,lshd)
  real(8):: rc2,rcx,rcy,rcz,rcxi,rcyi,rczi
  integer:: lcx,lcy,lcz,lcxyz,lcyz,lcx2,lcy2,lcz2,lcyz2,lcxyz2
  real(8),allocatable:: tmparr(:)
  integer:: ndmax

!.....Cutoff list and
  real(8),allocatable:: rclst(:)
  integer,allocatable:: idxlst_rc(:,:)

!.....Arrays for Gonnet's algorithm
  real(8),allocatable:: dlist(:)
  integer,allocatable:: ilist(:),ileft(:),ia2ic(:)
  real(8):: dh(3,3),vecc(3,26)
  integer:: idcell(-1:1,-1:1,-1:1)
  
contains
!=======================================================================
  subroutine mk_lscl_para()
!
! Make a linked cell list.
! Codes are slightly different bewteen parallel and single.
!
    use pmdvars,only: namax,natm,nbmax,nb,ra,anxi,anyi,anzi,rc,rbuf, &
         h,hi
!!$    integer,intent(in):: namax,natm,nbmax,nb
!!$    real(8),intent(in):: anxi,anyi,anzi,rc,h(3,3),hi(3,3)
!!$    real(8),intent(inout):: ra(3,namax)

    integer:: i,mx,my,mz,m
    real(8):: rcut

    rcut = rc +rbuf

!.....The following code (block) is used to be called once at the 1st call,
!.....but these variables could change when the cell size and shape change.
    rc2= rcut**2
!.....Number of division along each lattice vector
    lcx=anxi/dsqrt(hi(1,1)**2+hi(1,2)**2+hi(1,3)**2)/rcut
    lcy=anyi/dsqrt(hi(2,1)**2+hi(2,2)**2+hi(2,3)**2)/rcut
    lcz=anzi/dsqrt(hi(3,1)**2+hi(3,2)**2+hi(3,3)**2)/rcut
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
    if( .not.allocated(lscl) ) then
      allocate(lscl(namax),lshd(lcxyz2))
      call accum_mem('pairlist',4*size(lscl)+4*size(lshd))
    endif

    if( size(lshd).ne.lcxyz2 ) then
      call accum_mem('pairlist',-4*size(lshd))
      deallocate(lshd)
      allocate(lshd(lcxyz2))
      call accum_mem('pairlist',4*size(lshd))
    endif

    if( size(lscl).ne.namax ) then
      call accum_mem('pairlist',-4*size(lscl))
      deallocate(lscl)
      allocate(lscl(namax))
      call accum_mem('pairlist',4*size(lscl))
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
!$acc update device(lshd,lscl)

    return
  end subroutine mk_lscl_para
!=======================================================================
  subroutine sort_arrays(namax,natm,nb,tag,ra,va,aux,naux)
!
!  Sort arrays (tag,ra,va, and some more) using the cell list
!  NOTE: This is not going to work well for now...,
!    and thus it should not be used...
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
  end subroutine sort_arrays
!=======================================================================
  subroutine mk_lspr_para(l1st)
!
!  Make a pairlist using cell list created in mk_lscl_para.
!  Cutoff radius is already set in the mk_lscl_para, and thus not given
!  as an argument.
!
    use pmdvars,only: namax,natm,nbmax,nb,nnmax,maxnn,tag,ra,va,h,hi,&
         anxi,anyi,anzi,lspr,iprint,rc,rbuf,mpi_md_world,myid_md, &
         ratio_nnmax_update
    implicit none
    include "mpif.h"
    logical,intent(in):: l1st

    integer:: i,j,k,l,m,n,inc,nni,nnj,maxnnl,ierr,ierr_max,ierr_maxg
    integer:: mx,my,mz,kux,kuy,kuz,m1x,m1y,m1z,m1,ic,jc,mmax
    real(8):: xi(3),xij(3),rij(3),rij2
    integer:: nnmax_prev

    call mk_lscl_para()

    call set_nnmax(l1st)

    if( .not.allocated(lspr) ) then
      allocate(lspr(0:nnmax,namax))
      call accum_mem('pairlist',4*size(lspr))
    endif

    maxnn = nnmax
    ierr_maxg = 0

10  continue
    if( ierr_maxg .gt. 0 ) then
      nnmax_prev = nnmax
      nnmax = int(maxnn * ratio_nnmax_update)
      if( myid_md.eq.0 .and. iprint.gt.0 ) then
        write(6,'(a,i4,a,i4)') ' Max num of neighbors is updated from ', &
             nnmax_prev,' to ',nnmax
        call flush(6)
      endif
      ierr_maxg = 0
    endif
    
    if( size(lspr).ne.(nnmax+1)*namax ) then
      call accum_mem('pairlist',-4*size(lspr))
      deallocate(lspr)
      allocate(lspr(0:nnmax,namax))
      call accum_mem('pairlist',4*size(lspr))
    endif

    maxnnl = 0
    ierr_max = 0
!-----make a pair list, LSPR
!.....Scan atoms (not scanning cells)
!$omp parallel
!$omp do private(mz,my,mx,m,kuz,kuy,kux,m1z,m1y,m1x,m1,i,j,xi,xij,rij,rij2) &
!$omp    reduction(max:maxnnl,ierr_max)

!$acc data present(ra,h,lspr,lshd,lscl), &
!$acc      copyin(rcx,rcy,rcz,rcxi,rcyi,rczi,lcx,lcy,lcz, &
!$acc             lcz2,lcyz2,rc2,natm,nnmax), &
!$acc      pcreate(xi,xij,rij)
!$acc kernels
!$acc loop independent private(xi,xij,rij) &
!$acc      reduction(max:maxnnl,ierr_max) gang worker vector
    do i=1,natm
      lspr(:,i) = 0  ! initialize
      xi(1:3) = ra(1:3,i)
!-------assign a vector cell index
      mx=(xi(1)+rcx)*rcxi
      my=(xi(2)+rcy)*rcyi
      mz=(xi(3)+rcz)*rczi
      mx= min(max(mx,0),lcx+1)
      my= min(max(my,0),lcy+1)
      mz= min(max(mz,0),lcz+1)
      m= mx*lcyz2 +my*lcz2 +mz +1
!$acc loop seq
      do kux= -1,1
        m1x= mx +kux
        if( m1x.lt.0 .or. m1x.gt.lcx+1 ) cycle
!$acc loop seq
        do kuy= -1,1
          m1y= my +kuy
          if( m1y.lt.0 .or. m1y.gt.lcy+1 ) cycle
!$acc loop seq
          do kuz= -1,1
            m1z= mz +kuz
            if( m1z.lt.0 .or. m1z.gt.lcz+1 ) cycle
            m1=m1x*lcyz2 +m1y*lcz2 +m1z +1
            if(lshd(m1).eq.0) cycle
            j=lshd(m1)
            do while( j.gt.0 )
              if( j.eq.i ) then
                j = lscl(j)
                cycle
              endif
              xij(1:3)= ra(1:3,j) -xi(1:3)
              rij(1)= h(1,1,0)*xij(1) +h(1,2,0)*xij(2) +h(1,3,0)*xij(3)
              rij(2)= h(2,1,0)*xij(1) +h(2,2,0)*xij(2) +h(2,3,0)*xij(3)
              rij(3)= h(3,1,0)*xij(1) +h(3,2,0)*xij(2) +h(3,3,0)*xij(3)
              rij2= rij(1)**2 +rij(2)**2 +rij(3)**2

              if( rij2.lt.rc2 ) then
                lspr(0,i) = lspr(0,i) +1
                if( lspr(0,i).gt.nnmax ) then
                  ierr_max = 1
                else
                  lspr(lspr(0,i),i) = j
                endif
              endif

              j=lscl(j)
            enddo ! while (j.gt.0)
          enddo  ! kuz
        enddo  ! kuy
      enddo  ! kux
      maxnnl = max(maxnnl,lspr(0,i))
    enddo  ! i=1,natm
!$acc end kernels
!$acc end data
    
!$omp end do
!$omp end parallel

    call mpi_allreduce(ierr_max,ierr_maxg,1,mpi_integer,mpi_max, &
         mpi_md_world,ierr)
    call mpi_allreduce(maxnnl,maxnn,1,mpi_integer,mpi_max, &
         mpi_md_world,ierr)
    if( ierr_maxg.gt.0 ) then
      goto 10
    endif


  end subroutine mk_lspr_para
!=======================================================================
  subroutine mk_lspr_gonnet(namax,natm,nbmax,nb,nnmax,tag,ra,va &
       ,rc,h,hi,anxi,anyi,anzi,lspr,iprint,l1st)
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

    call mk_lscl_para()

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
  subroutine mk_lspr_sngl(namax,natm,nnmax,tag,ra,rc,h,hi &
       ,lspr,iprint,l1st)
!
! Make lspr in serial implimentation taking the periodic boundary
! condition into account.
!
! NOTICE:
! This routine is called from a python script nappy/pmd/pairlist.py,
! so be careful when modfiying the code.
!
    implicit none
    integer,intent(in):: namax,natm,nnmax,iprint
    integer,intent(out):: lspr(0:nnmax,namax)
    real(8),intent(in):: ra(3,namax),rc,hi(3,3),h(3,3),tag(namax)
    logical,intent(in):: l1st

    integer:: i,j,k,l,m,n
    integer:: mx,my,mz,kux,kuy,kuz,m1x,m1y,m1z,m1,ic,jc,ierr
    real(8):: xi(3),xij(3),rij(3),rij2,shft(3)

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

!.....Reset header
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

!.....Reset pair list
    lspr(:,:)= 0
    do i=1,natm
      ierr = 0
      xi(1:3) = ra(1:3,i)
!.....Assign a vector cell index
      mx = (xi(1)+rcx)*rcxi
      my = (xi(2)+rcy)*rcyi
      mz = (xi(3)+rcz)*rczi
      mx= min(max(mx,1),lcx)
      my= min(max(my,1),lcy)
      mz= min(max(mz,1),lcz)
      m= (mx-1)*lcyz +(my-1)*lcz +mz
      do kux=-1,1
        shft(1) = 0d0
        m1x = mx +kux
        if( m1x.lt.1   ) then
          m1x= m1x +lcx
          shft(1) = -1d0
        else if( m1x.gt.lcx ) then
          m1x= m1x -lcx
          shft(1) = +1d0
        endif
        do kuy=-1,1
          shft(2) = 0d0
          m1y= my +kuy
          if( m1y.lt.1   ) then
            m1y= m1y +lcy
            shft(2) = -1d0
          else if( m1y.gt.lcy ) then
            m1y= m1y -lcy
            shft(2) = +1d0
          endif
          do kuz=-1,1
            shft(3) = 0d0
            m1z= mz +kuz
            if( m1z.lt.1   ) then
              m1z= m1z +lcz
              shft(3) = -1d0
            else if( m1z.gt.lcz ) then
              m1z= m1z -lcz
              shft(3) = +1d0
            endif
            m1= (m1x-1)*lcyz +(m1y-1)*lcz +m1z
            if( lshd(m1).eq.0 ) cycle
            j = lshd(m1)
            do while( j.gt.0 .and. ierr.le.0 )
              if( j.eq.i ) then
                j = lscl(j)
                cycle
              endif
              xij(1:3)= ra(1:3,j) -xi(1:3) +shft(1:3)
!!$              xij(1:3)= xij(1:3) -anint(xij(1:3))  !! This is wrong in this case.
              rij(1:3)=h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
              rij2= rij(1)**2 +rij(2)**2 +rij(3)**2
              if( rij2.lt.rc2 ) then
                lspr(0,i) = lspr(0,i) +1
                lspr(lspr(0,i),i) = j
                if( lspr(0,i).eq.nnmax ) ierr = 1
              endif

              j= lscl(j)
            enddo  ! while(j.gt.0)
          enddo  ! kuz
        enddo  ! kuy
      enddo  ! kux
    enddo  ! i=1,natm

  end subroutine mk_lspr_sngl
!=======================================================================
  subroutine mk_lspr_brute(namax,natm,nbmax,nb,nnmax,tag,ra,rc &
       ,h,hi,sgm,lspr,iprint,l1st)
!
!  Make pair list, lspr, by brute force approach, because the system
!  is supposed to be small. Expand the system to take all the atoms 
!  within given cutoff radius into account.
!
    implicit none
    integer,intent(in):: namax,natm,nnmax,nbmax,iprint
    integer,intent(out):: lspr(0:nnmax,natm),nb
    real(8),intent(in):: rc,hi(3,3),h(3,3),sgm(3,3)
    real(8),intent(out):: ra(3,namax),tag(namax)
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
!.....Store i in j's neighbor list
          lspr(0,ja)= lspr(0,ja)+1
          if( lspr(0,ja).gt.nnmax ) then
            write(6,'(a)') " Error: lspr(0,ja) > nnmax"
            write(6,'(a,3i5)') "  nnmax, lspr(0,ja) = " &
                 ,nnmax,lspr(0,ja)
            stop
          endif
          lspr(lspr(0,ja),ja)=ia
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
  subroutine check_lscl(myid,iprint)
    integer,intent(in):: myid,iprint

    integer:: mmax,mx,my,mz,i,inc,m

    if( myid.eq.0 .and. iprint.ge.ipl_info) then
      print *,'check_lscl:'
      print '(a,3i4)','   lcx,y,z = ',lcx,lcy,lcz
      print '(a,5(2x,i0))','   lcx2,y2,z2,yz2,xyz2 = ',lcx2,lcy2,lcz2,lcyz2,lcxyz2
      print '(a,3es12.4)','   rcx,y,z = ',rcx,rcy,rcz
    endif
    
    mmax = 0
    do mz=0,lcz+1
      do my=0,lcy+1
        do mx=0,lcx+1
          m= mx*lcyz2 +my*lcz2 +mz +1
          i = lshd(m)
          inc = 0
          do while(i.gt.0)
            inc = inc + 1
            i = lscl(i)
          enddo
          mmax = max(mmax,inc)
!!$          print '(a,5i6)',' mx,my,mz,m,inc=',m,inc
        enddo
      enddo
    enddo
  end subroutine check_lscl
!=======================================================================
  subroutine check_lspr(namax,natm,nnmax,lspr,iprint,myid,mpi_world)
    integer,intent(in):: namax,natm,nnmax,lspr(0:nnmax,namax),iprint
    integer,intent(in):: myid,mpi_world

    integer:: i,mmax
    real(8):: vratio
    real(8),parameter:: pi = 3.14159265358979d0

    mmax = 0
    do i=1,natm
      mmax = max(mmax,lspr(0,i))
    enddo
!.....vratio = sphere/cube
    vratio = 2d0**2*pi/3d0**4
    vratio = 1.25d0 *vratio
    mmax = int(mmax *vratio)
    if( mmax.gt.nnmax ) then
      print *,'ERROR: Max num of neighbors exceeds the limit:'
      print '(3x,i5,a,i5)', mmax,' > ',nnmax
      print *,'Increase max_num_neighbors in in.pmd...'
      stop
    endif
  end subroutine check_lspr
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
  subroutine sort_lspr(namax,natm,ra,nnmax,h,lspr)
!
!  Sort indices in lspr in ascending order of distance.
!
    integer,intent(in):: namax,natm,nnmax
    real(8),intent(in):: ra(3,namax),h(3,3)
    integer,intent(inout):: lspr(0:nnmax,namax)

    integer:: i,j,jj,nn
    real(8):: xi(3),xij(3),rij(3),dij2
    real(8),allocatable,save:: dists(:)
    integer,allocatable,save:: idxarr(:),itmparr(:)

    if( .not. allocated(dists) ) then
      allocate(dists(nnmax),idxarr(nnmax),itmparr(nnmax))
    else if( size(dists).ne.nnmax ) then
      deallocate(dists,idxarr,itmparr)
      allocate(dists(nnmax),idxarr(nnmax),itmparr(nnmax))
    endif
    
    do i=1,natm
      xi(1:3) = ra(1:3,i)
!.....Create a temporary list of neighbor distances
      dists(:) = 0d0
      nn = lspr(0,i)
      do jj=1,nn
        j = lspr(jj,i)
        xij(1:3) = ra(1:3,j) -xi(1:3)
        rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij2 = rij(1)**2 +rij(2)**2 +rij(3)**2
        dists(jj) = dsqrt(dij2)
      enddo
!.....Sort lspr according to dists
      call arg_heapsort_arr(nn,nnmax,dists,idxarr)
      itmparr(1:nn) = lspr(1:nn,i)
      do jj=1,nn
        lspr(jj,i) = itmparr(idxarr(jj))
      enddo
    enddo
    return
  end subroutine sort_lspr
!=======================================================================
  subroutine set_nnmax(l1st)
!
!  Estimate and set nnmax from lscl. Thus this must be called after 
!  lscl and lshd are computed.
!  The density rho is determined by the max density among lscl. 
!  And the nnmax is determined by 4*pi*rho*rcut**3 *alpha /3,
!  where alpha is a mergin of the estimate, like 1.2.
!
    use pmdvars,only: vol,myid_md,mpi_md_world,nnmax,nxyz,rc,rbuf, &
         iprint,ratio_nnmax_update,ntot
    include "mpif.h"
    logical,intent(in):: l1st
    
    integer:: ierr,ic,i,inc,nmaxl,nmax,nnmax_estimate,nnmax_prev, &
         ncell,ncelltot
    integer:: ix,iy,iz
    real(8):: volc,rho
    real(8),parameter:: pi = 3.14159265358979d0
!!$    logical,save:: l1st = .true. 

    nmaxl = 0
    ncell = 0  ! num. of cells that contain at least one atom
    do ix=1,lcx
      do iy=1,lcy
        do iz=1,lcz
          ic = ix*lcyz2 +iy*lcz2 +iz +1
          i = lshd(ic)
          inc = 0
          do while( i.gt.0 )
            inc = inc +1
            i = lscl(i)
          enddo ! while (i.gt.0)
          nmaxl = max(inc,nmaxl)
          if( inc > 0 ) ncell = ncell +1
        enddo
      enddo
    enddo

    nmax = 0
    call mpi_allreduce(nmaxl,nmax,1,mpi_integer,mpi_max,mpi_md_world, &
         ierr)
    ncelltot = 0
    call mpi_allreduce(ncell,ncelltot,1,mpi_integer,mpi_sum,mpi_md_world, &
         ierr)

    volc = vol/(lcx*lcy*lcz)/nxyz
!!$    rho = dble(nmax)/volc
    rho = dble(ntot)/(volc*ncelltot)
    nnmax_estimate = int(4*pi*(rc+rbuf)**3*rho/3) +1
!!$    print *,'vol,lcx,lcy,lcz=',vol,lcx,lcy,lcz
!!$    print *,'myid,nmaxl,nmax,volc,rho,nnmax_estimate=', &
!!$         myid_md,nmaxl,nmax,volc,rho,nnmax_estimate
!.....If nnmax_estmate < nnmax, do nothing and return
    if( nnmax_estimate.le.nnmax ) return
      
    nnmax_prev = nnmax
    nnmax = int(nnmax_estimate *ratio_nnmax_update)

    if( myid_md.eq.0 .and. iprint.gt.0 ) then
      if( l1st ) then
        write(6,*) ''
        write(6,'(a)') ' Estimation of num of neighbors:'
        write(6,'(a,i5)') '   Max num in link-list cell = ',nmax
        write(6,'(a,f0.1,3x,f6.4)') '   Cell volume and density = ',volc,rho
        write(6,'(a,2(2x,i0))') '   Max num of neighbors, that incl. margin (nnmax) = ', &
             nnmax_estimate, nnmax
        call flush(6)
      else
        write(6,'(a,i4,a,i4)') ' Max num of neighbors is updated from ', &
             nnmax_prev,' to ',nnmax
        call flush(6)
      endif
    endif
!!$    l1st = .false.
    return
  end subroutine set_nnmax
!=======================================================================
end module pairlist
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd lib"
!     End:
