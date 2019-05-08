module structure
!-----------------------------------------------------------------------
!                     Last modified: <2019-05-08 15:38:48 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
!  Routines of structure analyses.
!-----------------------------------------------------------------------
  implicit none
  save

  integer:: mem = 0
  real(8):: time = 0d0

!.....CNA
  integer,parameter:: lmax= 12
  integer,parameter:: mmax= 12
  integer,allocatable:: icommon(:),ibond(:,:),nb(:),idc(:,:,:)
  integer,allocatable:: idcna(:) ! CNA type
  logical:: lcna1st = .true.
  
  
contains
  subroutine cna(namax,natm,nbndr,nnmax,lspr,tag)
!-----------------------------------------------------------------------
! Common Neighbor Analysis using the given neighbor list, LSPR.
! Computes a digit for each atom that represents crystalline structure in IDCNA:
!   0: other
!   1: FCC
!   2: HCP
!   3: BCC
!-----------------------------------------------------------------------
    implicit none
    integer,intent(in):: namax,natm,nbndr,nnmax,lspr(0:nnmax,namax)
    real(8),intent(in):: tag(namax)

    integer:: i,j,l,m,n,ii,iii,ni,jj,nj,il,jl,n1,n2,iil,nn1,im,iim &
         ,ib1,ib2,iib1,iib2,n421,n422,n663,n443

!c-----for serial code
!      allocate(icommon(lmax),ibond(2,mmax),nb(mmax) &
!          ,idc(3,nnmax,natm))
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

!-----init three indices
    idc(1:3,1:nnmax,1:natm+nbndr)= 0

!-----for each atom-i, store three indices (LMN)
    do i=1,natm+nbndr

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
              if(l.gt.lmax) stop " l.gt.lmax!!!"
              icommon(l)= ni
              exit
            endif
          enddo
!---------end of counting L
        enddo

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
                if(m.gt.mmax) stop " m.gt.mmax!!"
                ibond(1:2,m)= (/ n1,n2 /)
                exit
              endif
            enddo
          enddo
        enddo

!---------count max num of continuous bonds: N
        nb(1:mmax)= 1
!---------this does not distinguish star and chain connections
        do im=1,m-1
          ib1= ibond(1,im)
          ib2= ibond(2,im)
          do iim=im+1,m
            iib1= ibond(1,iim)
            iib2= ibond(2,iim)
!-------------if two bonds are connected, up nb
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
!-------end of 1st n.n. of i: j
      enddo
!-----end of atom-i
    enddo

!.....Identify crystalline structure and assign a digit to each atom
    do i=1,natm
      n421= 0
      n422= 0
      n663= 0
      n443= 0
      print *,'i,ls1nn=',i,lspr(0:8,i)
      do ii=1,lspr(0,i)
        l=idc(1,ii,i)
        m=idc(2,ii,i)
        n=idc(3,ii,i)
        if(l.eq.4 .and. m.eq.2 .and. n.eq.1 ) n421=n421 +1
        if(l.eq.4 .and. m.eq.2 .and. n.eq.2 ) n422=n422 +1
        if(l.eq.6 .and. m.eq.6 .and. n.eq.3 ) n663=n663 +1
        if(l.eq.4 .and. m.eq.4 .and. n.eq.3 ) n443=n443 +1
        print *,'i,ii,l,m,n=',i,ii,l,m,n
      enddo
      
!.....FCC structure
      if(n421.eq.12 .and. n422.eq.0) then
        idcna(i) = 1

!.....HCP structure
      else if(n421.eq.6 .and. n422.eq.6) then
        idcna(i) = 2

!.....BCC structure
      else if(n663.eq.8 .and. n443.eq.6) then
        idcna(i) = 3
        
!.....Other or unclassified
      else
        idcna(i) = 0
      endif
      print *,'i,idcna=',i,idcna(i)
    enddo

  end subroutine cna
  
end module structure
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
