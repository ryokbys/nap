module structure
!-----------------------------------------------------------------------
!                     Last modified: <2019-05-09 13:20:46 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
!  Routines of structure analyses.
!-----------------------------------------------------------------------
  implicit none
  save

  integer:: mem = 0
  real(8):: time = 0d0

!.....CNA
  integer,parameter:: idfcc = 1
  integer,parameter:: idhcp = 2
  integer,parameter:: idbcc = 3
  integer,parameter:: idunclass = 0
  integer,parameter:: lmax = 12
  integer,parameter:: mmax = 14
  integer,allocatable:: icommon(:),ibond(:,:),nb(:),idc(:,:,:)
  integer,allocatable:: idcna(:) ! CNA type
  logical:: lcna1st = .true.

!.....adaptive-CNA
  integer,parameter:: nnsortmax = 14
  
contains
!=======================================================================
  subroutine cna_indices(namax,natm,nbndr,nnmax,lspr)
!-----------------------------------------------------------------------
! Computes three CNA indices (l,m,n) and store them in IDC variable.
!-----------------------------------------------------------------------
    implicit none
    integer,intent(in):: namax,natm,nbndr,nnmax,lspr(0:nnmax,namax)

    integer:: i,j,l,m,n,ii,iii,ni,jj,nj,il,jl,n1,n2,iil,nn1,im,iim &
         ,ib1,ib2,iib1,iib2

!.....Initialize three indices yet to be detemined
    do i=1,natm+nbndr
      if( idcna(i).ne.0 ) idc(1:3,1:nnmax,i) = 0
    enddo

!.....Compute and store three indices (LMN)
    do i=1,natm+nbndr

      if( idcna(i).ne.0 ) cycle ! No need to perform CNA for atom-i

!-------for each 1st n.n.
      do ii=1,lspr(0,i)
        j=lspr(ii,i)
!---------j>i only
        if(j.le.i) cycle

!---------count num of common neighbors: L
        l= 0
        icommon(1:lmax)= 0
        do iii=1,lspr(0,i)
          ni=lspr(iii,i)
          if(ni.eq.j) cycle
          do jj=1,lspr(0,j)
            nj=lspr(jj,j)
            if(nj.eq.ni) then
              l=l+1
              if(l.gt.lmax) then
!!$                stop " l.gt.lmax!!!"
                goto 1
              endif
              icommon(l)= ni
              exit
            endif
          enddo
!---------end of counting L
        enddo
1       continue

!---------count num of bonds between common neighbors: M
        m= 0
        ibond(1:2,1:mmax)= 0
!---------for each common neighbor-n1
        do il=1,l
          n1=icommon(il)
!-----------for common neighbor-n2 which must be larger than n1
          do jl=1,l
            n2=icommon(jl)
            if(n2.le.n1) cycle
!-------------scan 1st n.n. of n1
            do iil=1,lspr(0,n1)
              nn1=lspr(iil,n1)
              if(nn1.eq.n2) then
                m=m+1
                if(m.gt.mmax) then
!!$                  print *,'ERROR: m>mmax !!'
!!$                  print *,'  i,jj,j,l,m=',i,jj,j,l,m
!!$                  stop 1
!.....The case (m>mmax) is probably irregular and assigned to
!     unknown structure. So skip looking for further bonds.
                  goto 2
                endif
                ibond(1:2,m)= (/ n1,n2 /)
                exit
              endif
            enddo
          enddo
        enddo
2       continue

!---------count max num of continuous bonds: N
        nb(1:mmax)= 1
!---------this does not distinguish star and chain connections
        do im=1,m-1
          ib1= ibond(1,im)
          ib2= ibond(2,im)
          do iim=im+1,m
            iib1= ibond(1,iim)
            iib2= ibond(2,iim)
!-----if two bonds are connected, increase nb
            if(iib1.eq.ib1 .or. iib2.eq.ib1 &
                 .or. iib1.eq.ib2 .or. iib2.eq.ib2) then
              nb(im)=nb(im) +1
              nb(iim)=nb(iim) +1
            endif
          enddo
        enddo
!---------maximum nb
        n= 0
        do im=1,m
          n= max(nb(im),n)
        enddo

!---------store (LMN) to i
        idc(1:3,ii,i)= (/ l,m,n /)
!---------store (LMN) to j, too
        do jj=1,lspr(0,j)
          if(lspr(jj,j).eq.i) then
            idc(1:3,jj,j)= (/ l,m,n /)
            exit
          endif
        enddo
      enddo  ! do ii=...
    enddo  ! do i=...

  end subroutine cna_indices
!=======================================================================
  subroutine cna(namax,natm,nbndr,nnmax,ls1nn,tag)
!-----------------------------------------------------------------------
! Common Neighbor Analysis using the neighbor list for 1st NN, LS1NN.
! Computes a digit for each atom that represents crystalline structure:
!   0: other
!   1: FCC
!   2: HCP
!   3: BCC
!-----------------------------------------------------------------------
    implicit none
    integer,intent(in):: namax,natm,nbndr,nnmax,ls1nn(0:nnmax,namax)
    real(8),intent(in):: tag(namax)

    integer:: i,j,l,m,n,ii,iii,ni,jj,nj,il,jl,n1,n2,iil,nn1,im,iim &
         ,ib1,ib2,iib1,iib2,n421,n422,n663,n443

!-----for parallel code
    if( lcna1st ) then
      allocate(icommon(lmax),ibond(2,mmax),nb(mmax),idc(3,nnmax,namax)&
           ,idcna(namax))
      mem = mem +4*size(icommon) +4*size(ibond) +4*size(nb) &
           +4*size(idc) +4*size(idcna)
      lcna1st = .false.
    endif

    if( size(idcna).ne.namax ) then
      mem = mem -4*size(idc) -4*size(idcna)
      deallocate(idc,idcna)
      allocate(idc(3,nnmax,namax),idcna(namax))
    endif

    idcna(:) = 0

    call cna_indices(namax,natm,nbndr,nnmax,ls1nn)
    call assign_struct_cna(namax,natm,nbndr,nnmax,ls1nn)

  end subroutine cna
!=======================================================================
  subroutine acna(namax,natm,nbndr,nnmax,lspr,h,ra,tag)
!-----------------------------------------------------------------------
! Adaptive CNA using the full neighbor list, LSPR.
! See, [1] Stukowski, MSMSE 20 (2012) 045021.
! Sort neighbors in the ascending order of their distances.
! Compute a digit for each atom that represents crystalline structure:
!   0: other
!   1: FCC
!   2: HCP
!   3: BCC
!-----------------------------------------------------------------------
    implicit none
    integer,intent(in):: namax,natm,nbndr,nnmax
    real(8),intent(in):: h(3,3),ra(3,namax),tag(namax)
    integer,intent(inout):: lspr(0:nnmax,namax)

    integer:: i,j,l,m,n,ii,iii,ni,jj,kk,nj,il,jl,n1,n2,iil,nn1,im,iim &
         ,ib1,ib2,iib1,iib2,itmp,lm
    real(8):: tmp,xi(3),xj(3),xij(3),rij(3),dij,dsum,dsum1,dsum2
    real(8),allocatable,save:: dists(:,:),rcfccs(:),rcbccs(:)
    integer,allocatable,save:: lsnn(:,:)

!-----for parallel code
    if( lcna1st ) then
      allocate(icommon(lmax),ibond(2,mmax),nb(mmax),idc(3,nnmax,namax)&
           ,idcna(namax),dists(nnmax,namax),lsnn(0:nnsortmax,namax) &
           ,rcfccs(namax),rcbccs(namax))
      mem = mem +4*size(icommon) +4*size(ibond) +4*size(nb) &
           +4*size(idc) +4*size(idcna) +8*size(dists) +4*size(lsnn) &
           +8*size(rcfccs) +8*size(rcbccs)
      lcna1st = .false.
    endif

    if( size(idcna).ne.namax ) then
      mem = mem -4*size(idc) -4*size(idcna) -4*size(lsnn) -8*size(dists) &
           -8*size(rcfccs) -8*size(rcbccs)
      deallocate(idc,idcna)
      allocate(idc(3,nnmax,namax),idcna(namax),lsnn(0:nnsortmax,namax) &
           ,dists(nnmax,namax))
      mem = mem +4*size(idc) +4*size(idcna) +4*size(lsnn) +8*size(dists) &
           +8*size(rcfccs) +8*size(rcbccs)
    endif

    idcna(:) = 0
    dists(:,:) = 0d0

!.....Compute all the distances between neighbors
    do i=1,natm+nbndr
      xi(1:3) = ra(1:3,i)
      do jj=1,lspr(0,i)
        j = lspr(jj,i)
        xj(1:3) = ra(1:3,j)
        xij(1:3) = xj(1:3) -xi(1:3)
        rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij = sqrt(rij(1)**2 +rij(2)**2 +rij(3)**2)
        dists(jj,i) = dij
      enddo
    enddo ! i=...


!.....Compute local cutoff for each atom
    do i=1,natm+nbndr
!.....Sort LSPR in the ascending order of the distance only up to nnsortmax
      do kk=1,min(lspr(0,i),nnsortmax)
        do jj=lspr(0,i),kk+1,-1
!.....Exchange jj-th and (jj+1)-th of lspr and dists if dists(jj)>dists(jj+1)
          if( dists(jj-1,i).gt.dists(jj,i) ) then
            tmp = dists(jj-1,i)
            dists(jj-1,i) = dists(jj,i)
            dists(jj,i) = tmp
            itmp = lspr(jj-1,i)
            lspr(jj-1,i) = lspr(jj,i)
            lspr(jj,i) = itmp
          endif
        enddo ! jj
      enddo ! kk
!!$      if( i.eq.1 .or. i.eq.50 ) then
!!$        do jj=1,nnsortmax
!!$          j = lspr(jj,i)
!!$          print *,'i,jj,j,dist=',i,jj,j,dists(jj,i)
!!$        enddo
!!$      endif
!.....Local cutoff for fcc (use only up to 12 neighbors)
      dsum = 0d0
      do jj=1,12
        dsum = dsum +dists(jj,i)
      enddo
      rcfccs(i) = (1d0+sqrt(2d0))/2 *dsum/12

!.....Local cutoff for bcc (use up to 14 neighbors)
      dsum1 = 0d0
      do jj=1,8
        dsum1 = dsum1 +dists(jj,i)
      enddo
      dsum2 = 0d0
      do jj=9,14
        dsum2 = dsum2 +dists(jj,i)
      enddo
!.....NOTE: rcbcc is different from the definition in Ref [1]
      rcbccs(i) = (1d0+sqrt(2d0))/4 *( 2d0/sqrt(3d0)*dsum1/8 +dsum2/6)
    enddo

!.....Use nearest 12 neighbors to identify FCC and HCP
    lsnn(0,:) = 0
    do i=1,natm+nbndr
      do jj=1,min(nnsortmax,lspr(0,i))
        if( dists(jj,i).gt.rcfccs(i) ) cycle
        lsnn(0,i) = lsnn(0,i) + 1
        lsnn(jj,i) = lspr(jj,i)
      enddo
    enddo
    call cna_indices(namax,natm,nbndr,nnsortmax,lsnn)
    call assign_struct_cna(namax,natm,nbndr,nnsortmax,lsnn)

!.....Use 14 neighbors to identify BCC
    lsnn(0,:) = 0
    do i=1,natm+nbndr
      do jj=1,min(nnsortmax,lspr(0,i))
        if( dists(jj,i).gt.rcbccs(i) ) cycle
        lsnn(0,i) = lsnn(0,i) + 1
        lsnn(jj,i) = lspr(jj,i)
      enddo
    enddo
    call cna_indices(namax,natm,nbndr,nnsortmax,lsnn)
    call assign_struct_cna(namax,natm,nbndr,nnsortmax,lsnn)

  end subroutine acna
!=======================================================================
  subroutine assign_struct_cna(namax,natm,nbndr,nnmax,lspr)
!
!  Assign crystalline structure from CNA indices
!
    integer,intent(in):: namax,natm,nbndr,nnmax,lspr(0:nnmax,namax)

    integer:: i,ii,l,m,n,n421,n422,n663,n443
    
!.....Identify crystalline structure and assign a digit to each atom
    do i=1,natm
      if( idcna(i).ne.0 ) cycle  ! skip already assigned atom
      n421= 0
      n422= 0
      n663= 0
      n443= 0
      do ii=1,lspr(0,i)
        l=idc(1,ii,i)
        m=idc(2,ii,i)
        n=idc(3,ii,i)
        if(l.eq.4 .and. m.eq.2 .and. n.eq.1 ) n421=n421 +1
        if(l.eq.4 .and. m.eq.2 .and. n.eq.2 ) n422=n422 +1
        if(l.eq.6 .and. m.eq.6 .and. n.eq.3 ) n663=n663 +1
        if(l.eq.4 .and. m.eq.4 .and. n.eq.3 ) n443=n443 +1
!!$        print *,'i,ii,l,m,n=',i,ii,l,m,n
      enddo
!.....FCC structure
      if(n421.eq.12 .and. n422.eq.0) then
        idcna(i) = idfcc

!.....HCP structure
      else if(n421.eq.6 .and. n422.eq.6) then
        idcna(i) = idhcp

!.....BCC structure
      else if(n663.eq.8 .and. n443.eq.6) then
        idcna(i) = idbcc
        
!.....Unclassified
      else
        idcna(i) = idunclass
      endif
!!$      print *,'i,idcna=',i,idcna(i)
    enddo
  end subroutine assign_struct_cna
end module structure
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
