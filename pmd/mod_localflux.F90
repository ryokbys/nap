module localflux
!
!  Module for evaluation of local flux using the color-charge NEMD.
!
  use pmdvars,only: nspmax
  use pmdmpi,only: nid2xyz
  use memory,only: accum_mem
  implicit none
  include 'mpif.h'
  include "./const.h"
  save

  character(len=128):: paramsdir = '.'
  character(len=128),parameter:: cfoutlf = 'out.lflux'
  integer,parameter:: ionum = 67

  logical:: lflux = .false.  ! Flag to eval local flux.
  logical:: initialized = .false.
  integer:: noutlflux = 1000
  integer:: nskipout
  real(8),allocatable:: fluxg(:,:)  ! Global flux
  real(8),allocatable:: fluxl(:,:),fluxt(:,:)  ! Local flux
  integer:: nlx = 1  ! Num of divisions in a node for local regions
  integer:: nly = 1  ! where flux is evaluated within.
  integer:: nlz = 1
  integer:: nl
  integer:: ng,ngx,ngy,ngz  ! Global num of divisions
  integer,allocatable:: idxl2g(:,:)  ! Conversion of index from local to global
  integer,allocatable:: nargn(:)  ! Num of atoms in local regions
  real(8):: dlx,dly,dlz,dlxi,dlyi,dlzi,evflux(3)

contains
!=======================================================================
  subroutine init_lflux(myid,nx,ny,nz,lclrchg,maxstp,mpi_world,iprint)
!
!  Initialize lflux.
!
    use clrchg,only: clrfield
    integer,intent(in):: myid,mpi_world,nx,ny,nz,iprint,maxstp
    logical,intent(in):: lclrchg
    
    integer:: inc,ix,iy,iz,mx,my,mz,ixyz,nxyz
    integer:: idxlfg,ilfgx,ilfgy,ilfgz,mem
    real(8):: anxi,anyi,anzi,fext
    integer,allocatable:: idxl2gt(:)
    integer:: istat(mpi_status_size),itag,ierr

    if( noutlflux.eq.0 ) stop 'noutlflux must not be 0.'
    nskipout = max(maxstp/noutlflux,1)

!.....Set unit vector of external force on color charge
    fext = sqrt(clrfield(1)**2 +clrfield(2)**2 +clrfield(3)**2)
    evflux(1:3) = clrfield(1:3)/fext

!.....Parallel setting
    anxi= 1d0/nx
    anyi= 1d0/ny
    anzi= 1d0/nz
    nxyz = nx*ny*nz

!.....Number of local-flux regions in a node and global
    nl = nlx*nly*nlz
    ngx = nlx*nx
    ngy = nly*ny
    ngz = nlz*nz
    ng = ngx*ngy*ngz

    dlx = anxi /nlx
    dly = anyi /nly
    dlz = anzi /nlz
    dlxi = 1d0/dlx
    dlyi = 1d0/dly
    dlzi = 1d0/dlz

    allocate(fluxl(3,nl),fluxt(3,nl),idxl2g(nl,0:nxyz-1),nargn(nl))
    allocate(idxl2gt(nl))
    mem = 8*(size(fluxl)+size(fluxt)) +4*(size(idxl2g)+size(nargn))
    call accum_mem('localflux',mem)
    fluxl(:,:) = 0d0
    nargn(:) = 0
    if( myid.eq.0 ) then
      allocate(fluxg(3,ng))
      call accum_mem('localflux',8*size(fluxg))
      fluxg(:,:) = 0d0
    endif

!.....Index conversion from local to global
    call nid2xyz(myid,mx,my,mz)
    inc = 0
    do ix=1,nlx
      ilfgx = mx*nlx +ix
      do iy=1,nly
        ilfgy = my*nly +iy
        do iz=1,nlz
          ilfgz = mz*nlz +iz
          inc = inc +1
          idxlfg = (ilfgx-1)*ngz*ngy +(ilfgy-1)*ngz +ilfgz
          idxl2gt(inc) = idxlfg
!!$          print '(a,6i6)','myid,ilfgx,ilfgy,ilfgz,inc,idxlfg=' &
!!$               ,myid,ilfgx,ilfgy,ilfgz,inc,idxlfg
        enddo
      enddo
    enddo
!.....Gather idxl2gt to make idxl2g in node-0
    if( myid.eq.0 ) then
      idxl2g(1:nl,0) = idxl2gt(1:nl)
      do ixyz=1,nxyz-1
        itag = ixyz
        call mpi_recv(idxl2gt,nl,mpi_integer,ixyz,itag,mpi_world &
             ,istat,ierr)
        idxl2g(1:nl,ixyz) = idxl2gt(1:nl)
      enddo
    else  ! myid.ne.0
      itag = myid
      call mpi_send(idxl2gt,nl,mpi_integer,0,itag,mpi_world,ierr)
    endif

    if( .not.lclrchg ) then
      if( myid.eq.0 ) then
        print *,'ERROR: color charge NEMD is not ON.'
        print *,'       Local flux is only available with color charge NEMD.'
      endif
      stop
    endif

    if( myid.eq.0 .and. iprint.ge.ipl_basic ) then
      print *,''
      print '(a)',' Local flux measuring ON:'
      print '(a,4(1x,i0))','   Number of regions in a node (x,y,z,total) =' &
           ,nlx,nly,nlz,nl
      print '(a,4(1x,i0))','   Number of regions globally (x,y,z,total)  =' &
           ,ngx,ngy,ngz,ng
      print '(a,4f7.3)','   Normalized lengths(x,y,z) of local region = ' &
           ,dlx,dly,dlz
      print '(a,3f7.3)', '   Unit vector along Fext = ',evflux(1:3)
      print '(a)', '   Output file = '//trim(cfoutlf)

      open(ionum,file=trim(cfoutlf),status='replace')
      write(ionum,'(a)') '#'
      write(ionum,'(a,4(1x,i0))') '# Local flux of regions, ngx,ngy,ngz,ng = '&
           ,ngx,ngy,ngz,ng
      write(ionum,'(a)') '# Order of data: idx=(igx-1)*ngy*ngz +(igy-1)*ngz +igz'
      write(ionum,'(a)') '# Entry:  istp, time [fs], (flux(ig),ig=1,ng) [Ang]'
      write(ionum,'(a)') '#'
    endif
    initialized = .true.
    
    deallocate(idxl2gt)
    return
  end subroutine init_lflux
!=======================================================================
  subroutine accum_lflux(namax,natm,hmat,ra,va,clr,istp,dt&
       ,myid,mpi_world,nxyz)
!
!  Accumurate velocities multiplying color charge for evaluation of
!  local flux.
!
    integer,intent(in):: namax,natm,istp,myid,mpi_world,nxyz
    real(8),intent(in):: hmat(3,3),ra(3,namax),va(3,namax),clr(namax),dt

    integer:: i,ilx,ily,ilz,idxl,idxg,ixyz
    real(8):: vi(3),flux
    integer,parameter:: nmpi = 2
    integer:: istat(mpi_status_size),itag,ierr

    nargn(:) = 0
    fluxt(:,:) = 0d0
    do i=1,natm
      if( abs(clr(i)).lt.0.9 ) cycle
!.....Get local region index
      ilx = int((ra(1,i)+dlx)*dlxi)
      ily = int((ra(2,i)+dly)*dlyi)
      ilz = int((ra(3,i)+dlz)*dlzi)
      idxl = (ilx-1)*nly*nlz +(ily-1)*nlz +ilz
      nargn(idxl) = nargn(idxl) +1
!.....Scale velocity
      vi(1:3) = hmat(1:3,1)*va(1,i) +hmat(1:3,2)*va(2,i) &
           +hmat(1:3,3)*va(3,i)
!.....Accumulate velocities in flux of its local region
      fluxt(1:3,idxl) = fluxt(1:3,idxl) +clr(i)*vi(1:3)*dt
    enddo
    do idxl=1,nl
      if( nargn(idxl).eq.0 ) cycle
      fluxl(1:3,idxl) = fluxl(1:3,idxl) +fluxt(1:3,idxl)/nargn(idxl)
    enddo

!.....Output
    if( mod(istp,nskipout).eq.0 ) then
!.....Reduce local fluxes in each node to global local-flux
      if( myid.eq.0 ) then
        do idxl=1,nl
          idxg = idxl2g(idxl,0)
          fluxg(1:3,idxg) = fluxg(1:3,idxg) +fluxl(1:3,idxl)
        enddo
        do ixyz=1,nxyz-1
          itag = ixyz*nmpi -nmpi
          fluxl(:,:) = 0d0
          call mpi_recv(fluxl,3*nl,mpi_real8,ixyz,itag+1,mpi_world &
               ,istat,ierr)
          do idxl=1,nl
            idxg = idxl2g(idxl,ixyz)
            fluxg(1:3,idxg) = fluxg(1:3,idxg) +fluxl(1:3,idxl)
          enddo
        enddo
!.....Write out fluxg only at node-0
        write(ionum,'(i10,es18.10)',advance='no') istp, dt*istp
        do idxg=1,ng
          flux = fluxg(1,idxg)*evflux(1) &
               +fluxg(2,idxg)*evflux(2) +fluxg(3,idxg)*evflux(3)
          write(ionum,'(1x,es11.3)',advance='no') flux
        enddo
        write(ionum,'(a)') ''
        
      else  ! myid.ne.0
        itag = myid*nmpi -nmpi
        call mpi_send(fluxl,3*nl,mpi_real8,0,itag+1,mpi_world,ierr)
      endif
!.....Refresh fluxl
      fluxl(:,:) = 0d0
    endif
    
  end subroutine accum_lflux
!=======================================================================
  subroutine final_lflux(myid)
!
!  Finalize local flux.
!
    integer,intent(in):: myid

    if( myid.eq.0 ) then
      close(ionum)
    endif
  end subroutine final_lflux
end module localflux
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
