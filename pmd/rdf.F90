program rdf
!-----------------------------------------------------------------------
!                     Last-modified: <2021-02-27 14:21:08 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
! Compute RDF.
!-----------------------------------------------------------------------
! INPUT:
! ------
!   - pmdini:        Structure file.
!   - in.rdf:        Input file that contains some parameters.
!
! OUTPUT:
! -------
!   - STDOUT should not be redirected to out.rdf
!   - out.rdf:       RDF data
!   - out.sfac:      Structure factor S(q) data obtained from total RDF.
!-----------------------------------------------------------------------
  use pmdio,only: read_pmdtot_ascii, get_ntot_ascii
  use pairlist,only: mk_lspr_sngl
  implicit none
  include 'mpif.h'
  real(8),parameter:: pi = 3.14159265358979d0
  character(len=128),parameter:: cpmdini='pmdini'
  character(len=128),parameter:: cfinput='in.rdf'
  character(len=128),parameter:: cfoutrdf='out.rdf'
  character(len=128),parameter:: cfoutsq='out.sfac'
  integer,parameter:: nspmax = 9

  integer:: ia,ic,nc,maxnn,is,js,msp,inc,ict,i,ib,l,n,nbins
  integer:: nspc(nspmax)
  integer:: ntot
  real(8),allocatable:: tagtot(:),rtot(:,:),vtot(:,:),atot(:,:)
  real(8),allocatable:: stot(:,:,:),epitot(:),ekitot(:,:,:)
  real(8),allocatable:: auxtot(:,:)
  integer,allocatable:: lspr(:,:),ls1nn(:,:)
  real(8),allocatable:: rdfs(:,:,:),denoms(:,:),sqs(:)
  logical:: lpair(nspmax,nspmax),lspc(nspmax)
  real(8):: rcut,hi(3,3),t0,tmp,vol,dr,r
  integer:: jb,jbm
  real(8):: dq,q,qcut,tmp1,tmp2,rm,rho
  character(len=3):: csi,csj
  character(len=20):: cnum

  t0 = mpi_wtime()

!.....Read atom configuration
  ntot = get_ntot_bin(10,trim(cpmdini))
  allocate(tagtot(ntot),rtot(3,ntot),vtot(3,ntot),epitot(ntot) &
       ,ekitot(3,3,ntot),stot(3,3,ntot),atot(3,ntot))
  call read_pmdtot_ascii(10,trim(cpmdini),ntot,hunit,h,tagtot,rtot,vtot)
  call get_hi(h,hi,vol)
  print *,'Num of atoms = ',ntot

!.....Default values
  rcut = 5.0d0
  qcut = 25.0d0
  nbins = 100
  call read_in_rdf(11,trim(cfinput),rcut,qcut,nbins,lpair)
  dr = rcut/nbins
  dq = qcut/nbins
  print *,'nbins = ',nbins
  print *,'rcut  = ',rcut
  print *,'dr    = ',dr
  print *,'qcut  = ',qcut
  print *,'dq    = ',dq

!.....Make neighbor list
  allocate(lspr(0:nnmax,ntot),ls1nn(0:nnmax,ntot))
  call mk_lspr_sngl(ntot,ntot,nnmax,tagtot,rtot,rcut,rcut,h,hi, &
       lspr,ls1nn,iprint,.true.)
  maxnn = 0
  msp = 0
  nspc(:) = 0
  do ia=1,ntot
    is = int(tagtot(ia))
    maxnn = max(maxnn,lspr(0,ia))
    msp = max(msp,is)
    nspc(is) = nspc(is) +1
  enddo
  print *,'Max num of neighbors = ',maxnn
  print *,'Max species ID = ',msp
  do is=1,msp
    print *,'Num of '//trim(specorder(is))//' = ',nspc(is)
  enddo

  lspc(:) = .false.
  do is=1,msp
    do js=1,msp
      if( lpair(is,js) ) then
        lspc(is) = .true.
        exit
      endif
    enddo
  enddo

!.....Compute RDF
  tmp = mpi_wtime()
  allocate(rdfs(0:msp,0:msp,nbins),denoms(0:msp,0:msp))
  call comp_rdf(ntot,tagtot,h,rtot,rcut,nnmax,lspr,msp,nbins,rdfs)
  print '(a,f0.3)',' Time for comp RDF = ',mpi_wtime()-tmp

!.....Aggregate pairwise to get total RDF
!!$  do is=1,msp
!!$    do js=is,msp
!!$      if( is.ne.js) then
!!$        do ib=1,nbins
!!$          rdfs(is,js,ib) = rdfs(is,js,ib) +rdfs(is,js,ib)
!!$        enddo
!!$      endif
!!$    enddo
!!$  enddo

!.....Normalization factors
  denoms(0,0) = 4d0*pi*ntot*(ntot-1)*dr/vol
  do is=1,msp
    denoms(is,is) = 4d0*pi*dr/vol *nspc(is)*(nspc(is)-1)
    do js=is+1,msp
      denoms(is,js) = 4d0*pi*dr/vol *nspc(is)*nspc(js)
    enddo
  enddo

!.....Normalize
  do ib=1,nbins
    r = (ib-0.5d0)*dr
    rdfs(0,0,ib) = rdfs(0,0,ib) /(denoms(0,0)*r*r)
    do is=1,msp
      do js=is,msp
        rdfs(is,js,ib) = rdfs(is,js,ib) /(denoms(is,js)*r*r)
      enddo
    enddo
  enddo


!.....Output RDF
  open(20,file=trim(cfoutrdf),status='replace')
!.....Header
  write(20,'(a)',advance='no') '#  1:distance, 2:total rdf, '
  inc = 2
  do is=1,msp
    csi = specorder(is)
    do js=is,msp
      if( .not.lpair(is,js) ) cycle
      inc = inc +1
      csj = specorder(js)
      write(20,'(i0,a)',advance='no') inc,':'//trim(csi)//'-'//trim(csj)//', '
    enddo
  enddo
  write(20,*) ''
  do ib=1,nbins
    write(20,'(2es14.4)',advance='no') (ib-0.5d0)*dr, rdfs(0,0,ib)
    do is=1,msp
      do js=is,msp
        if( .not.lpair(is,js) ) cycle
        write(20,'(es14.4)',advance='no') rdfs(is,js,ib)
      enddo
    enddo
    write(20,*) ''
  enddo
  close(20)

!.....S(q)
  allocate(sqs(nbins))
  sqs(:) = 0d0
  rho = dble(ntot)/vol
  print *,'ntot,vol,rho=',ntot,vol,rho
  do ib=1,nbins
    q = (ib-0.5d0)*dq
    tmp = 0d0
    do jb=2,nbins
      r = (jb-0.5d0)*dr
      jbm = jb-1
      rm = (jbm-0.5d0)*dr
      tmp1= (rdfs(0,0,jbm)-1d0)*sin(q*rm)/(q*rm)*rm*rm
      tmp2= (rdfs(0,0,jb)-1d0)*sin(q*r)/(q*r)*r*r
      tmp = tmp + 0.5d0*dr *(tmp1+tmp2)
    enddo
    sqs(ib)= 1d0 +4*pi*rho*tmp
  enddo
!.....Output S(Q)
  open(21,file=trim(cfoutsq),status='replace')
  write(21,'(a)') '#  1:wave number, 2:S(Q) '
  do ib=1,nbins
    q = (ib-0.5d0)*dq
    write(21,'(2es14.4)') q, sqs(ib)
  enddo
  close(21)
  print '(a,f0.3)',' Time = ',mpi_wtime() -t0

end program rdf
!=======================================================================
subroutine read_in_rdf(ionum,cfname,rcut,qcut,nbins,lpair)
!-----------------------------------------------------------------------
!  Input format should be like the following:
!-----------------------------------------------------------------------
!  rcut   5.0
!  nbins  200
!  pairs  Si-Si Si-O
!-----------------------------------------------------------------------
!  If pairs is specified, output pairwise RDF in addition to total RDF.
!-----------------------------------------------------------------------
  use pmdvars,only: csp2isp,nspmax
  use util, only: num_data
  implicit none 
  integer,intent(in):: ionum
  character(len=*),intent(in):: cfname
  real(8),intent(inout):: rcut,qcut
  integer,intent(inout):: nbins
  logical,intent(out):: lpair(nspmax,nspmax)

  integer,parameter:: maxpair = 5

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
    else if( trim(c1st).eq.'qcut' ) then
      if( nentry.ne.2 ) then
        print *,'ERROR: nentry != 2 !!!'
        stop 1
      endif
      read(cline,*) c1st, qcut
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
    else if( trim(c1st).eq.'nbins' ) then
      if( nentry.ne.2 ) then
        print *,'ERROR: nentry != 2 !!!'
        stop 1
      endif
      read(cline,*) c1st, nbins
    else if( trim(c1st).eq.'max_num_neighbors' ) then
      if( nentry.ne.2 ) then
        print *,'ERROR: nentry != 2 !!!'
        stop 1
      endif
      read(cline,*) c1st, nnmax
    endif
  end do
1 close(ionum)

  print *,' Read '//trim(cfname)
  
end subroutine read_in_rdf
!=======================================================================
subroutine get_hi(h,hi,vol)
  implicit none 
  real(8),intent(in):: h(3,3)
  real(8),intent(out):: hi(3,3),vol

  integer:: i,j,im,ip,jm,jp
  real(8):: sgm(3,3),hit(3,3)

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
subroutine comp_rdf(ntot,tagtot,h,rtot,rcut,nnmax,lspr,msp,nbins,rdfs)
!
!  Compute RDF.
!
  implicit none
  integer,intent(in):: ntot,nnmax,lspr(0:nnmax,ntot),msp,nbins
  real(8),intent(in):: tagtot(ntot),rtot(3,ntot),rcut,h(3,3)
  real(8),intent(out):: rdfs(0:msp,0:msp,nbins)

  integer:: ia,is,ja,js,jj,ib
  real(8):: dr,rc2,xi(3),xj(3),xij(3),rij(3),dij2,dij

  dr = rcut/nbins
  rc2 = rcut*rcut
  rdfs(:,:,:) = 0d0
  do ia=1,ntot
    is = int(tagtot(ia))
    xi(1:3) = rtot(1:3,ia)
    do jj=1,lspr(0,ia)
      ja = lspr(jj,ia)
      js = int(tagtot(ja))
      xj(1:3)= rtot(1:3,ja)
      xij(1:3)= xj(1:3)-xi(1:3) -anint(xj(1:3)-xi(1:3))
      rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
      dij2 = rij(1)**2 +rij(2)**2 +rij(3)**2
      if( dij2.gt.rc2 ) cycle
      dij = dsqrt(dij2)
      ib = min(int(dij/dr)+1,nbins)
      rdfs(0,0,ib) = rdfs(0,0,ib) + 1d0
      rdfs(is,js,ib) = rdfs(is,js,ib) + 1d0
!!$      if( js.ne.is ) rdfs(js,is,ib) = rdfs(js,is,ib) + 1d0
    enddo
  enddo
  return
end subroutine comp_rdf
!=======================================================================
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make rdf"
!     End:
