module pdens
!
!  Module for evaluation of probability density.
!
  use pmdio,only: csp2isp,nspmax,get_vol
  use pmdmpi,only: nid2xyz
  use element,only: get_cube_info
  implicit none
  include 'mpif.h'
  save

  character(len=128):: paramsdir = '.'
  character(len=10),parameter:: cfprefix = 'out.pdens.'
  character(len=5),parameter:: cfpostfix = '.cube'
  character(len=128):: cfoutpd
  integer,parameter:: ionum = 68

  logical:: lpdens = .false.  ! Flag to eval local flux.
  logical:: initialized = .false.
  character(len=3):: cspc_pdens = 'non'
  integer:: ispc_pdens
  real(8),allocatable:: pdg(:)  ! Global pdens
  real(8),allocatable:: pdl(:),pdt(:)  ! Local pdens
  integer:: npx = 1  ! Num of divisions in a node for local regions
  integer:: npy = 1  ! where flux is evaluated within.
  integer:: npz = 1
  integer:: nl
  integer:: ng,ngx,ngy,ngz  ! Global num of divisions
  integer,allocatable:: idxl2g(:,:)  ! Conversion of index from local to global
  real(8):: dpx,dpy,dpz,dpxi,dpyi,dpzi
  integer:: nacc

contains
!=======================================================================
  subroutine init_pdens(myid,nx,ny,nz,hmat,specorder,mpi_world,iprint)
!
!  Initialize lflux.
!
    integer,intent(in):: myid,mpi_world,nx,ny,nz,iprint
    character(len=3),intent(in):: specorder(nspmax)
    real(8),intent(in):: hmat(3,3)
    
    integer:: inc,ix,iy,iz,mx,my,mz,ixyz,nxyz
    integer:: idxg,ipx,ipy,ipz
    real(8):: anxi,anyi,anzi,fext
    integer,allocatable:: idxl2gt(:)
    integer:: istat(mpi_status_size),itag,ierr

    if( trim(cspc_pdens).eq.'non' ) then
      stop 'ERROR: spcs_pdens must be specified.'
    else
      ispc_pdens = csp2isp(trim(cspc_pdens),specorder)
      cfoutpd = cfprefix//trim(cspc_pdens)//cfpostfix
    endif

!.....Parallel setting
    anxi= 1d0/nx
    anyi= 1d0/ny
    anzi= 1d0/nz
    nxyz = nx*ny*nz

!.....Number of local-flux regions in a node and global
    nl = npx*npy*npz
    ngx = npx*nx
    ngy = npy*ny
    ngz = npz*nz
    ng = ngx*ngy*ngz

    dpx = anxi /npx
    dpy = anyi /npy
    dpz = anzi /npz
    dpxi = 1d0/dpx
    dpyi = 1d0/dpy
    dpzi = 1d0/dpz

    allocate(pdl(nl),pdt(nl),idxl2g(nl,0:nxyz-1))
    allocate(idxl2gt(nl))
    pdl(:) = 0d0

!.....Index conversion from local to global
    call nid2xyz(myid,mx,my,mz)
    inc = 0
    do ix=1,npx
      ipx = mx*npx +ix
      do iy=1,npy
        ipy = my*npy +iy
        do iz=1,npz
          ipz = mz*npz +iz
          inc = inc +1
          idxg = (ipx-1)*ngz*ngy +(ipy-1)*ngz +ipz
          idxl2gt(inc) = idxg
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

    if( myid.eq.0 .and. iprint.gt.0 ) then
      print *,''
      print '(a)',' Probability density measuring ON:'
      print '(a)','   Tracked species = '//trim(cspc_pdens)
      print '(a,4(1x,i0))','   Number of regions in a node (x,y,z,total) =' &
           ,npx,npy,npz,nl
      print '(a,4(1x,i0))','   Number of regions globally (x,y,z,total)  =' &
           ,ngx,ngy,ngz,ng
      print '(a,4f7.4)','   Normalized lengths(x,y,z) of local region = ' &
           ,dpx,dpy,dpz
      print '(a)', '   Output file = '//trim(cfoutpd)

    endif
    nacc = 0
    initialized = .true.
    
    deallocate(idxl2gt)
    return
  end subroutine init_pdens
!=======================================================================
  subroutine accum_pdens(namax,natm,tag,ra)
!
!  Accumurate density of specified species.
!
    integer,intent(in):: namax,natm
    real(8),intent(in):: ra(3,namax),tag(namax)

    integer:: i,ipx,ipy,ipz,idxl,is
    integer,parameter:: nmpi = 2

    do i=1,natm
      is = int(tag(i))
      if( is.ne.ispc_pdens ) cycle
!.....Get local region index
      ipx = int((ra(1,i)+dpx)*dpxi)
      ipy = int((ra(2,i)+dpy)*dpyi)
      ipz = int((ra(3,i)+dpz)*dpzi)
!!$      ipx = min(max(ipx,1),npx)
!!$      ipy = min(max(ipy,1),npy)
!!$      ipz = min(max(ipz,1),npz)
      idxl = (ipx-1)*npy*npz +(ipy-1)*npz +ipz
      pdl(idxl) = pdl(idxl) +1d0
    enddo
    nacc = nacc +1
    
  end subroutine accum_pdens
!=======================================================================
  subroutine final_pdens(myid,mpi_world,nxyz,hmat,natm,tag,ra,specorder)
!
!  Finalize local flux.
!
    include 'params_unit.h'
    integer,intent(in):: myid,mpi_world,nxyz,natm
    real(8),intent(in):: hmat(3,3),ra(3,natm),tag(natm)
    character(len=*),intent(in):: specorder(nspmax) 
    
    integer,parameter:: nmpi = 1
    integer:: idxl,idxg,ixyz,is,num,ia,nums(nspmax)
    integer:: istat(mpi_status_size),itag,ierr
    real(8):: vol,val,ri(3),vals(nspmax)

!.....Reduce local fluxes in each node to global local-flux
    if( myid.eq.0 ) then
      vol = get_vol(hmat)/ng
      allocate(pdg(ng))
      pdg(:) = 0d0
      do idxl=1,nl
        idxg = idxl2g(idxl,0)
        pdg(idxg) = pdg(idxg) +pdl(idxl)
      enddo
      do ixyz=1,nxyz-1
        itag = ixyz*nmpi -nmpi
        pdl(:) = 0d0
        call mpi_recv(pdl,nl,mpi_real8,ixyz,itag+1,mpi_world &
             ,istat,ierr)
        do idxl=1,nl
          idxg = idxl2g(idxl,ixyz)
          pdg(idxg) = pdg(idxg) +pdl(idxl)
        enddo
      enddo
!.....Write out pdens only at node-0
      open(ionum,file=trim(cfoutpd),status='replace')
      write(ionum,'(a)') '# Probability density in Gaussian cube format.'
      write(ionum,'(a)') '# NOTE: length unit in Bohr.'
      write(ionum,'(2x,i0,3(1x,f7.3))') natm, 0d0, 0d0, 0d0
      write(ionum,'(2x,i0,3(1x,es15.7))') ngx,hmat(1:3,1)*ang2bohr/ngx
      write(ionum,'(2x,i0,3(1x,es15.7))') ngy,hmat(1:3,2)*ang2bohr/ngy
      write(ionum,'(2x,i0,3(1x,es15.7))') ngz,hmat(1:3,3)*ang2bohr/ngz
      call get_cube_info(nspmax,specorder,nums,vals)
      do ia=1,natm
        is = int(tag(ia))
        num = nums(is)
        val = vals(is)
        ri(1:3)= hmat(1:3,1)*ra(1,ia) +hmat(1:3,2)*ra(2,ia) +hmat(1:3,3)*ra(3,ia)
        write(ionum,'(i6,f8.3,3(1x,es12.4))') num,val,ri(1:3)*ang2bohr
      enddo
!!$      write(ionum,'(a)') '  1   1.000   0.000  0.000  0.000'
      write(ionum,'(6(2x,es11.3))') (pdg(idxg)/nacc/(vol*ang2bohr**3),idxg=1,ng)
      close(ionum)
      deallocate(pdg)
    else  ! myid.ne.0
      itag = myid*nmpi -nmpi
      call mpi_send(pdl,nl,mpi_real8,0,itag+1,mpi_world,ierr)
    endif
    call mpi_barrier(mpi_world,ierr)
    deallocate(pdl,pdt,idxl2g)
    return
  end subroutine final_pdens
end module pdens
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
