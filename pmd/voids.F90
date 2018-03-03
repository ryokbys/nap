!-----------------------------------------------------------------------
!                     Last-modified: <2018-03-03 21:49:52 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
program voids
!-----------------------------------------------------------------------
!  Create volumetric data to detect voids.
!  OUTPUT:
!  - out.voids
!    > The mesh point having an atom within cutoff distance ==> 0
!    > otherwise (void mesh) ==> -1
!-----------------------------------------------------------------------
  use pmdio
  implicit none

  integer:: nargc,ix,iy,iz,nx,ny,nz,nxyz,n,i,j,js,iflag,kux,kuy,kuz
  integer:: lcx,lcy,lcz,lcyz,lcxyz,mx,my,mz,m,m1x,m1y,m1z,m1
  integer,allocatable:: imesh(:,:,:),lscl(:),lshd(:)
  real(8):: vol,sgm(3,3),hi(3,3)
  real(8):: alx,aly,alz,x,y,z,rcx,rcy,rcz,rcxi,rcyi,rczi,pi(3) &
       ,xij(3),rij(3),rij2,rc2,sidelen
  real(8),allocatable:: pmesh(:,:)
  character(len=128):: cusage,cfname,ctmp
!.....Functions
  integer,external:: iargc

  cusage = '  Usage: '//&
       ' $ ./voids pmdfile length cutoff'

!.....Get command-line arguments
  nargc = iargc()
  if( nargc.ne.3 ) then
    print *,'Error: nargc.ne.3 '
    print *, trim(cusage)
    stop
  endif
  call getarg(1,cfname)
  call getarg(2,ctmp)
  read(ctmp,*) sidelen
  call getarg(3,ctmp)
  read(ctmp,*) rc

!.....Read system info
  call read_pmdtot_ascii(10,trim(cfname))

!.....Make volumetric data mesh
  call boxmat(h,hi,vol,sgm)
  call get_lattice_lengths(h,alx,aly,alz)
  nx = int(alx/sidelen) +1
  ny = int(aly/sidelen) +1
  nz = int(alz/sidelen) +1
  nxyz = nx*ny*nz
  print *,'h-mat:'
  print *, h(1,1:3,0)
  print *, h(2,1:3,0)
  print *, h(3,1:3,0)
  print *,'nx,ny,nz,nxyz = ',nx,ny,nz,nxyz
  allocate(imesh(nx,ny,nz),pmesh(3,nxyz))
  imesh(:,:,:) = 0
  n = 0
  do ix=1,nx
    x = (dble(ix)-0.5)/nx
    do iy=1,ny
      y = (dble(iy)-0.5)/ny
      do iz=1,nz
        z = (dble(iz)-0.5)/nz
        n= n+1
        pmesh(1,n) = x
        pmesh(2,n) = y
        pmesh(3,n) = z
      enddo
    enddo
  enddo

!.....Make cell list
  lcx= 1d0/dsqrt(hi(1,1)**2+hi(1,2)**2+hi(1,3)**2)/rc
  lcy= 1d0/dsqrt(hi(2,1)**2+hi(2,2)**2+hi(2,3)**2)/rc
  lcz= 1d0/dsqrt(hi(3,1)**2+hi(3,2)**2+hi(3,3)**2)/rc
  if( lcx.lt.2 .or. lcy.lt.2 .or. lcz.lt.2 ) then
    write(6,'(a)') ' Error: Cannot handle' &
         //' too small system !!!'
    stop
  endif
  lcyz= lcy*lcx
  lcxyz= lcx*lcyz
  print *,'lcx,lcy,lcz,lcxyz = ',lcx,lcy,lcz,lcxyz
  rcx= 1d0/lcx
  rcy= 1d0/lcy
  rcz= 1d0/lcz
  rcxi=1d0/rcx
  rcyi=1d0/rcy
  rczi=1d0/rcz
!-----allocate LSCL & LSHD after obtaining lcxyz
  allocate(lscl(ntot),lshd(lcxyz))
  lshd(1:lcxyz) = 0
  do i=1,ntot
!.....Assign a vector cell index
    mx = (rtot(1,i)+rcx)*rcxi
    my = (rtot(2,i)+rcy)*rcyi
    mz = (rtot(3,i)+rcz)*rczi
    mx = min(max(mx,1),lcx)
    my = min(max(my,1),lcy)
    mz = min(max(mz,1),lcz)
    m = (mx-1)*lcyz +(my-1)*lcz +mz
!!$    print *,'i,rtot,mx,my,mz,m=',i,rtot(1:3,i),mx,my,mz,m
    lscl(i) = lshd(m)
    lshd(m) = i
  enddo

!.....Search distances bewteen cell-mesh points and atoms
  rc2 = rc*rc
  n = 0
  do ix=1,nx
    do iy=1,ny
      do iz=1,nz
        n = n +1
        pi(1:3) = pmesh(1:3,n)
        mx = int(pi(1)*rcxi) +1
        my = int(pi(2)*rcyi) +1
        mz = int(pi(3)*rczi) +1
        m = (mx-1)*lcyz +(my-1)*lcz +mz
        iflag = -1
        do kuz=-1,1
          m1z = mz + kuz
          if( m1z.lt.1 ) m1z= m1z +lcz
          if( m1z.gt.lcz ) m1z= m1z -lcz
          do kuy=-1,1
            m1y = my + kuy
            if( m1y.lt.1 ) m1y= m1y +lcy
            if( m1y.gt.lcy ) m1y= m1y -lcy
            do kux=-1,1
              m1x = mx + kux
              if( m1x.lt.1 ) m1x= m1x +lcx
              if( m1x.gt.lcx ) m1x= m1x -lcx
              m1 = (m1x-1)*lcyz +(m1y-1)*lcz +m1z
              if (lshd(m1).eq.0) cycle

              j = lshd(m1)
2             continue
              xij(1:3) = rtot(1:3,j) -pi(1:3) &
                   -anint(rtot(1:3,j) -pi(1:3))
              rij(1:3) = h(1:3,1,0)*xij(1) &
                   +h(1:3,2,0)*xij(2) &
                   +h(1:3,3,0)*xij(3)
              rij2 = rij(1)**2 +rij(2)**2 +rij(3)**2
              if( rij2.lt.rc2 ) then
                iflag = 0
                exit
              endif
              j = lscl(j)
              if( j.gt.0 ) goto 2
            enddo
            if( iflag.eq.0 ) exit
          enddo
          if( iflag.eq.0 ) exit
        enddo
        imesh(ix,iy,iz) = iflag
      enddo
    enddo
  enddo

  open(10,file='out.voids',status='replace')
  write(10,'(a,3i8)') '# ',nx,ny,nz
  n = 0
  do ix=1,nx
    do iy=1,ny
      do iz=1,nz
        n = n+1
        pi(1:3) = h(1:3,1,0)*pmesh(1,n) &
             +h(1:3,2,0)*pmesh(2,n) &
             +h(1:3,3,0)*pmesh(3,n)
        write(10,'(3es12.4,i4)') pi(1:3),imesh(ix,iy,iz)
      enddo
    enddo
  enddo
  close(10)

  open(11,file='out.ions',status='replace')
  write(11,'(a,i10)') '# ', ntot
  do n=1,ntot
    pi(1:3) = h(1:3,1,0)*rtot(1,n) &
         +h(1:3,2,0)*rtot(2,n) &
         +h(1:3,3,0)*rtot(3,n)
    write(11,'(3es12.4)') pi(1:3)
  enddo
  close(11)

end program voids
!=======================================================================
subroutine get_lattice_lengths(hmat,alx,aly,alz)
!
!  Assuming the cell is orthogonal.
!
  real(8),intent(in):: hmat(3,3)
  real(8),intent(out):: alx,aly,alz

  real(8):: a(3),b(3),c(3)

  a(1:3) = hmat(1:3,1)
  b(1:3) = hmat(1:3,2)
  c(1:3) = hmat(1:3,3)
  alx = sqrt(a(1)**2 +a(2)**2 +a(3)**2)
  aly = sqrt(b(1)**2 +b(2)**2 +b(3)**2)
  alz = sqrt(c(1)**2 +c(2)**2 +c(3)**2)
  return
  
end subroutine get_lattice_lengths
!=======================================================================
subroutine boxmat(h,hi,vol,sgm)
!-----------------------------------------------------------------------
!  setup matrices of MD-box
!    H:   MD-box matrix
!    HI:  inverse MD-box matrix
!    SGM: cofactor matrix
!-----------------------------------------------------------------------
  implicit none
  real(8),intent(in):: h(3,3,0:1)
  real(8),intent(out):: vol,sgm(3,3),hi(3,3)

  real(8):: hit(3,3)
  integer:: i,j,k,im,ip,jm,jp

!-----cofactor matrix, SGM
  do j=1,3
    jm=mod(j+1,3)+1
    jp=mod(j,  3)+1
    do i=1,3
      im=mod(i+1,3)+1
      ip=mod(i,  3)+1
      sgm(i,j)=h(ip,jp,0)*h(im,jm,0)-h(im,jp,0)*h(ip,jm,0)
    enddo
  enddo
!-----MD-box volume
  vol=h(1,1,0)*sgm(1,1)+h(2,1,0)*sgm(2,1)+h(3,1,0)*sgm(3,1)
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
end subroutine boxmat
!=======================================================================

