program cluster_analysis
!-----------------------------------------------------------------------
!                     Last-modified: <2022-07-08 16:40:42 KOBAYASHI Ryo>
!-----------------------------------------------------------------------
! Cluster analysis program.
! The cluster analysis is usually performed for large scale systems.
! The computational efficiency is an important issure in this case
! and python program is usually very slow. So here we make a fortran
! program of cluster analysis.
!-----------------------------------------------------------------------
! INPUT:
! ------
!   - in.cluster:    Main input file for `cluster` program
!   - pmdini:        Structure file.
!
! OUTPUT:
! -------
!   - STDOUT
!-----------------------------------------------------------------------
  use pmdvars
  use pmdio,only: read_pmdtot_ascii, get_ntot_ascii
  use pairlist,only: mk_lspr_sngl
  implicit none
  include 'mpif.h'
  integer,parameter:: maxpair = 5
  character(len=20),parameter:: cfinput='in.cluster'

  real(8):: hunit,hmat(3,3,0:1)
  real(8),allocatable:: tagtot(:),rtot(:,:),vtot(:,:),atot(:,:)

  integer:: ia,ic,nc,maxnn,is,js,msp,inc,ict,i,ib,l,n,nacmax
  integer,allocatable:: ictot(:),nacs(:),icouts(:),nhist(:)
  logical:: lpair(nspmax,nspmax),lspc(nspmax),lrecur
  real(8):: rcut,outthd,t0,tmp
  character(len=20):: cnum

  t0 = mpi_wtime()

!.....Read atom configuration
!!$  call read_pmdtot_ascii(10,trim(cpmdini))
  ntot = get_ntot_ascii(20,trim(cpmdini))
  allocate(tagtot(ntot),rtot(3,ntot),vtot(3,ntot))
  call read_pmdtot_ascii(20,trim(cpmdini),ntot,hunit,h,tagtot,rtot,vtot)
  call get_hi(h,hi)

!.....Read input
  rcut = 3.0d0
  outthd = 10
  lrecur = .false.
  call read_in_cluster(11,trim(cfinput),maxpair,rcut,lpair,outthd,lrecur)

!.....Make neighbor list
  allocate(lspr(0:nnmax,ntot),d2lspr(nnmax,ntot))
  call mk_lspr_sngl(ntot,ntot,nnmax,tagtot,rtot,rcut,h,hi, &
       lspr,d2lspr,iprint,.true.)
  maxnn = 0
  msp = 0
  do ia=1,ntot
    maxnn = max(maxnn,lspr(0,ia))
    msp = max(msp,int(tagtot(ia)))
  enddo
  print *,'Max num of neighbors = ',maxnn
  print *,'Max species ID = ',msp

  lspc(:) = .false.
  do is=1,nspmax
    do js=1,nspmax
      if( lpair(is,js) ) then
        lspc(is) = .true.
        exit
      endif
    enddo
  enddo

!.....Perform cluster analysis
  tmp = mpi_wtime()
  allocate(ictot(ntot))
  ictot(:) = 0
  if( lrecur ) then
    print *,'clustering with recursive routine.'
    call clustering2(ntot,tagtot,nnmax,lspr,nspmax,lpair,lspc,ictot,nc)
  else
    print *,'clustering without recursive routine.'
    call clustering1(ntot,tagtot,nnmax,lspr,nspmax,lpair,lspc,ictot,nc)
  endif
  print '(a,f0.3)',' Time for clustering = ',mpi_wtime()-tmp

!.....Comp num of atoms in each cluster
  allocate(nacs(nc))
  nacs(:) = 0
  do ia=1,ntot
    is = int(tagtot(ia))
    if( .not.lspc(is) ) cycle
    ic = ictot(ia)
    nacs(ic) = nacs(ic) + 1
  enddo
!.....Comp histogram vs num of atoms in a cluster
  nacmax = 0
  do ic=1,nc
    nacmax = max(nacmax,nacs(ic))
  enddo
  allocate(nhist(nacmax))
  nhist(:) = 0
  do ic=1,nc
    n = nacs(ic)
    if( n.eq.0 ) cycle
    nhist(n) = nhist(n) +1
  enddo
  open(30,file='out.histogram',status='replace')
  write(30,'(a)') '# num of atoms in cluster, count of clusters'
  do n=1,nacmax
    if( nhist(n).ne.0 ) write(30,'(2i10)') n,nhist(n)
  enddo
  close(30)
  
!.....Out pmdini_ic of cluster size larger than the threshold
  allocate(icouts(nc))
  icouts(:) = -1
  inc = 0
  do ic=1,nc
    if( nacs(ic).gt.outthd ) then
      inc = inc + 1
      icouts(inc) = ic
      print *,'ic,nac = ',ic,nacs(ic)
    endif
  enddo
  do ict=1,inc
    ic = icouts(ict)
    write(cnum,'(i0)') ic
    open(20,file='pmdini_'//trim(cnum),status='replace')
    write(20,'(a)') '!'
    write(20,'(a,9(2x,a))') '!  specorder: ',(trim(specorder(i)),i=1,msp)
    write(20,'(a)') '!'
    write(20,'(es23.14e3)') hunit
    write(20,'(3es23.14e3)') (((h(ia,ib,l)/hunit,ia=1,3) &
         ,ib=1,3),l=0,1)
    write(20,'(i10)') nacs(ic)
    do ia=1,ntot
      if( ictot(ia).eq.ic ) then
        write(20,'(7es23.14e3,11es13.4e3)') tagtot(ia) &
             ,rtot(1:3,ia),vtot(1:3,ia)    ! dt
      endif
    enddo
    close(20)
  enddo

  print '(a,f0.3)',' Time = ',mpi_wtime() -t0

end program cluster_analysis
!=======================================================================
subroutine read_in_cluster(ionum,cfname,maxpair,rcut,lpair,outthd,lrecur)
!-----------------------------------------------------------------------
!  The input format should be like the following.
!-----------------------------------------------------------------------
!  rcut   3.0
!  pairs  Si-Si Si-O
!  out_threshold   10
!  recursive  T
!-----------------------------------------------------------------------
  use pmdvars,only: nspmax
  use pmdio,only: split_pair
  use util, only: num_data,csp2isp
  implicit none 
  integer,intent(in):: ionum,maxpair
  character(len=*),intent(in):: cfname
  real(8),intent(out):: rcut,outthd
  logical,intent(out):: lpair(nspmax,nspmax),lrecur

  integer:: i,nentry,isp1,isp2,npair
  character(len=128):: c1st,cline
  character(len=7):: cpairs(maxpair)
  character(len=3):: csp1,csp2

!.....Set initial lpair as all False
  lpair(:,:) = .false.

  open(ionum,file=trim(cfname),status='old')
  do
    read(ionum,'(a)',end=1) cline
    nentry = num_data(cline,' ')
    if( nentry.eq.0 ) cycle
    if( cline(1:1).eq.'!' .or. cline(1:1).eq.'#' ) cycle
    read(cline,*) c1st
    if( trim(c1st).eq.'rcut' ) then
      if( nentry.ne.2 ) then
        print *,'ERROR: nentry != 2 !!!'
        stop 1
      endif
      read(cline,*) c1st, rcut
    else if( trim(c1st).eq.'pairs' ) then
      npair = nentry -1
      read(cline,*) c1st, (cpairs(i),i=1,npair)
!.....Convert cpairs to lpair
      do i=1,npair
        call split_pair(cpairs(i),csp1,csp2)
        isp1 = csp2isp(csp1)
        isp2 = csp2isp(csp2)
        lpair(isp1,isp2) = .true.
        lpair(isp2,isp1) = .true.  ! symmetrize
      enddo
    else if( trim(c1st).eq.'out_threshold' ) then
      read(cline,*) c1st, outthd
    else if( trim(c1st).eq.'recursive' ) then
      read(cline,*) c1st, lrecur
    endif
  end do
1 close(ionum)

  print *,' Read '//trim(cfname)
  
end subroutine read_in_cluster
!=======================================================================
subroutine get_hi(h,hi)
  implicit none 
  real(8),intent(in):: h(3,3)
  real(8),intent(out):: hi(3,3)

  integer:: i,j,im,ip,jm,jp
  real(8):: vol,sgm(3,3),hit(3,3)

!-----cofactor matrix, SGM
  do j=1,3
    jm=mod(j+1,3)+1
    jp=mod(j,  3)+1
    do i=1,3
      im=mod(i+1,3)+1
      ip=mod(i,  3)+1
      sgm(i,j)=h(ip,jp)*h(im,jm)-h(im,jp)*h(ip,jm)
    enddo
  enddo
!-----MD-box volume
  vol=h(1,1)*sgm(1,1)+h(2,1)*sgm(2,1)+h(3,1)*sgm(3,1)
  do j=1,3
    do i=1,3
      hit(i,j)= sgm(i,j)/vol
    enddo
  enddo
!-----transpose
  do j=1,3
    do i=1,3
      hi(i,j)= hit(j,i)
    enddo
  enddo
  return
end subroutine get_hi
!=======================================================================
subroutine clustering1(ntot,tagtot,nnmax,lspr,nspmax,lpair,lspc,ictot,nc)
!
!  Clustering not using the recursive routine.
!
  implicit none 
  integer,intent(in):: ntot,nnmax,lspr(0:nnmax,ntot),nspmax
  logical,intent(in):: lpair(nspmax,nspmax),lspc(nspmax)
  real(8),intent(in):: tagtot(ntot)
  integer,intent(out):: ictot(ntot),nc

  integer:: iter,icid,num_update,ia,is,icmin0,icmax,icmin,jj,ja,js,jc,ic

  iter = 0
  icid = 0
  do
    num_update = 0
    nc = 0
    do ia=1,ntot
      is = int(tagtot(ia))
      if( .not. lspc(is) ) cycle
      icmin0 = ictot(ia)
      icmax  = ictot(ia)
      icmin = ntot +1
      if( ictot(ia).ne.0 ) icmin = ictot(ia)
      do jj=1,lspr(0,ia)
        ja = lspr(jj,ia)
        js = int(tagtot(ja))
        if( .not. lpair(is,js) ) cycle
        jc = ictot(ja)
        icmin0 = min(icmin0,jc)
        icmax  = max(icmax,jc)
        if( jc.ne.0 ) icmin  = min(icmin,jc)
      enddo
      nc = max(nc,icmax)
!.....No need to update for these atoms
      if( icmax.eq.icmin .and. icmin0.ne.0 ) then
        cycle
!.....Determine minimum cluster-ID among neighbors or new cluster-ID
      else
        if( icmax.eq.0 ) then
          icid = icid + 1
          ic = icid
        else
          ic = icmin
        endif
        num_update = num_update +1
      endif
!.....Set the cluster-ID to all the neighbors
      if( ic.eq.0 ) then
        print *,'ia,ic,icmin0,icmin,icmax,icid=',ia,ic,icmin0,icmin,icmax,icid
        stop 1
      endif
      ictot(ia) = ic
      do jj=1,lspr(0,ia)
        ja = lspr(jj,ia)
        js = int(tagtot(ja))
        if( .not. lpair(is,js) ) cycle
        ictot(ja) = ic
      enddo
    enddo
    iter = iter + 1
    print '(a,2i10)',' iteration, num_update = ',iter,num_update
    if( num_update.eq.0 ) exit

  enddo
  print *,'num of clusters = ',nc
  return
end subroutine clustering1
!=======================================================================
subroutine clustering2(ntot,tagtot,nnmax,lspr,nspmax,lpair,lspc,ictot,nc)
!
!  Clustering using the recursive routine.
!
  implicit none 
  integer,intent(in):: ntot,nnmax,lspr(0:nnmax,ntot),nspmax
  logical,intent(in):: lpair(nspmax,nspmax),lspc(nspmax)
  real(8),intent(in):: tagtot(ntot)
  integer,intent(out):: ictot(ntot),nc

  integer:: is,js,ic,ia

  ic = 0
  do ia=1,ntot
    is = int(tagtot(ia))
    if( .not.lspc(is) ) cycle
    if( ictot(ia).ne.0 ) cycle
    ic = ic +1
!!$    print *,'ic,ia = ',ic,ia
    call neighbor_connection(ntot,ictot,tagtot,nnmax,lspr,nspmax,lpair,ia,is,ic)
  end do

  nc = ic

  return
end subroutine clustering2
!=======================================================================
recursive subroutine neighbor_connection(ntot,ictot,tagtot,nnmax,lspr,nspmax,&
     lpair,ia,is,ic)
  integer,intent(in):: ntot,nnmax,lspr(0:nnmax,ntot),nspmax,ia,is,ic
  logical,intent(in):: lpair(nspmax,nspmax)
  real(8),intent(in):: tagtot(ntot)
  integer,intent(inout):: ictot(ntot)

  ictot(ia) = ic
  
  do jj=1,lspr(0,ia)
    ja = lspr(jj,ia)
    js = int(tagtot(ja))
    if( lpair(is,js) .and. ictot(ja).eq.0 ) then
!!$      print *,'ia,ja,is,js = ',ia,ja,is,js
      call neighbor_connection(ntot,ictot,tagtot,nnmax,lspr,nspmax,lpair,ja,js,ic)
    endif
  enddo
  return
end subroutine neighbor_connection
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make cluster"
!     End:
