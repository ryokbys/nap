module pdens
!
!  Module for evaluation of probability density.
!
  use pmdmpi
  use mod_precision
  use memory,only: accum_mem
  implicit none
  include "./const.h"
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
  real(rp):: orig_pdens(3),hmat_pdens(3,3)
  real(rp):: sosub(3),shsub(3,3),shsubi(3,3)
  real(rp),allocatable:: pds(:,:,:)
  integer:: npx = 1  ! Num of divisions of subsystem read as ndiv_pdens in read_input
  integer:: npy = 1
  integer:: npz = 1
  integer:: np
  real(rp):: dpx,dpy,dpz,vol,dpxi,dpyi,dpzi
  integer:: nacc

contains
!=======================================================================
  subroutine init_pdens(myid,hmat,mpi_world,iprint)
!
!  Initialize lflux.
!
    use util,only: csp2isp, get_vol
    use vector,only: matinv3,matxvec3
    integer,intent(in):: myid,mpi_world,iprint
    real(rp),intent(in):: hmat(3,3)

    integer:: inc,ix,iy,iz,mx,my,mz,ixyz
    integer:: idxg,ipx,ipy,ipz
    real(rp):: hmati(3,3)
    integer,allocatable:: idxl2gt(:)
    integer:: istat(mpi_status_size),itag,ierr

    if( trim(cspc_pdens).eq.'non' ) then
      stop 'ERROR: spcs_pdens must be specified.'
    else
      ispc_pdens = csp2isp(trim(cspc_pdens))
      cfoutpd = cfprefix//trim(cspc_pdens)//cfpostfix
    endif

!.....Number of local-flux regions in a node and global
    np = npx*npy*npz
    dpx = 1.0_rp /npx
    dpy = 1.0_rp /npy
    dpz = 1.0_rp /npz
    dpxi= 1.0_rp /dpx   ! for efficiency
    dpyi= 1.0_rp /dpy
    dpzi= 1.0_rp /dpz

!.....Reset orig_pdens and hmat_pdens if hmat_pdens is not given
    if( get_vol(hmat_pdens).lt.1e-8_rp ) then
      orig_pdens(:) = 0.0_rp
      hmat_pdens(:,:) = hmat(:,:)
    endif
!.....Sub lattice representation in original hmat
    hmati = matinv3(hmat)
    sosub = matxvec3(hmati,orig_pdens)
    shsub(1:3,1) = matxvec3(hmati,hmat_pdens(1:3,1))
    shsub(1:3,2) = matxvec3(hmati,hmat_pdens(1:3,2))
    shsub(1:3,3) = matxvec3(hmati,hmat_pdens(1:3,3))
    shsubi = matinv3(shsub)

    allocate(pds(npz,npy,npx))
    call accum_mem('pdens',rp*size(pds))
    pds(:,:,:) = 0.0_rp

    if( myid.eq.0 .and. iprint.ge.ipl_basic ) then
      print *,''
      print '(a)',' Probability density measuring ON:'
      print '(a,i0,a)','   Tracked species = ',ispc_pdens,'('//trim(cspc_pdens)//')'
      print '(a,4(1x,i0))','   Number of regions (x,y,z,total) =' &
           ,npx,npy,npz,np
      print '(a,3(1x,f8.2))','   Origin of subsystem: ',orig_pdens
      print '(a)','   Lattice vectors of subsystem:'
      print '(5x,3(1x,f8.2))',hmat_pdens(1:3,1)
      print '(5x,3(1x,f8.2))',hmat_pdens(1:3,2)
      print '(5x,3(1x,f8.2))',hmat_pdens(1:3,3)
      print '(a,3(1x,f8.2))','   Origin of subsystem reduced by hmati: ',sosub(1:3)
      print '(a)','   Lattice vectors of subsystem reduced by hmati:'
      print '(5x,3(1x,f8.2))',shsub(1:3,1)
      print '(5x,3(1x,f8.2))',shsub(1:3,2)
      print '(5x,3(1x,f8.2))',shsub(1:3,3)
      print '(a)','   shsubi:'
      print '(5x,3(1x,f8.2))',shsubi(1:3,1)
      print '(5x,3(1x,f8.2))',shsubi(1:3,2)
      print '(5x,3(1x,f8.2))',shsubi(1:3,3)
      print '(a)', '   Output file = '//trim(cfoutpd)

    endif
    nacc = 0
    initialized = .true.
    return
  end subroutine init_pdens
!=======================================================================
  subroutine accum_pdens(namax,natm,tag_isp,ra,sorg)
!
!  Accumurate density of specified species.
!
    use vector,only: matxvec3
    integer,intent(in):: namax,natm
    integer,intent(in):: tag_isp(namax)
    real(rp),intent(in):: ra(3,namax),sorg(3)

    integer:: i,ipx,ipy,ipz,idx,is
    integer,parameter:: nmpi = 2
    real(rp):: ri(3),sri(3)

    do i=1,natm
      is = tag_isp(i)
      if( is.ne.ispc_pdens ) cycle
!.....Convert from hmat-rep to shsub-rep
      ri(1:3) = ra(1:3,i) +sorg(1:3) -sosub(1:3)
      sri(1:3) = matxvec3(shsubi,ri)
!.....Get subsystem index
      if(  sri(1).lt.0.0_rp .or. sri(1).ge.1.0_rp .or. &
           sri(2).lt.0.0_rp .or. sri(2).ge.1.0_rp .or. &
           sri(3).lt.0.0_rp .or. sri(3).ge.1.0_rp ) cycle
      ipx = int(sri(1)*dpxi) +1
      ipy = int(sri(2)*dpyi) +1
      ipz = int(sri(3)*dpzi) +1
      pds(ipz,ipy,ipx) = pds(ipz,ipy,ipx) +1.0_rp
    enddo
    nacc = nacc +1

  end subroutine accum_pdens
!=======================================================================
  subroutine final_pdens(myid,mpi_world,hmat)
!
!  Finalize prob density
!
    use util,only: get_vol
    include 'params_unit.h'
    integer,intent(in):: myid,mpi_world
    real(rp),intent(in):: hmat(3,3)

    integer,parameter:: nmpi = 1
    integer:: idx,ixyz,is,ia,ix,iy,iz,inc
    integer:: ierr
    real(rp):: vol,dr(3),fac
    real(rp),allocatable:: pdl(:,:,:)

!.....Reduce prob densities in each node to global prob density
    vol = get_vol(hmat)/np
    allocate(pdl(npz,npy,npx))
    call accum_mem('pdens',rp*size(pdl))
    pdl(:,:,:) = 0.0_rp
    call mpi_reduce(pds,pdl,npx*npy*npz,mpi_real_rp,mpi_sum,0,mpi_world,ierr)
!.....Write out pdens only at node-0
    if( myid.eq.0 ) then
      dr(1:3) = hmat_pdens(1,1:3)/npx +hmat_pdens(2,1:3)/npy +hmat_pdens(3,1:3)/npz
      open(ionum,file=trim(cfoutpd),status='replace')
      write(ionum,'(a)') '# Probability density in Gaussian cube format.'
      write(ionum,'(a)') '# NOTE: length unit in Bohr.'
!!$      write(ionum,'(2x,i0,3(1x,es15.7))') 1, (orig_pdens(1:3)+dr(1:3)/2)*ang2bohr
      write(ionum,'(2x,i0,3(1x,es15.7))') 1, orig_pdens(1:3)*ang2bohr
      write(ionum,'(2x,i0,3(1x,es15.7))') npx,hmat_pdens(1:3,1)*ang2bohr/npx
      write(ionum,'(2x,i0,3(1x,es15.7))') npy,hmat_pdens(1:3,2)*ang2bohr/npy
      write(ionum,'(2x,i0,3(1x,es15.7))') npz,hmat_pdens(1:3,3)*ang2bohr/npz
!.....Put a line for dummy atom, which is necessary to be loaded by Ovito.
      write(ionum,'(a)') '  1   1.000   0.000  0.000  0.000'
!.....Volumetric data
      fac = 1.0_rp /nacc /(vol*ang2bohr**3)
      pdl(:,:,:) = pdl(:,:,:) *fac
      inc= 0
      do ix=1,npx
        do iy=1,npy
          do iz=1,npz
            write(ionum,'(2x,es11.3)',advance='no') pdl(iz,iy,ix)
            inc = inc +1
            if( mod(inc,6) .eq. 0 ) then
              write(ionum,*) ''
            endif
          enddo
        enddo
      enddo
!!$      write(ionum,'(6(2x,es11.3))') (pdl(idx)/nacc/(vol*ang2bohr**3),idx=1,np)
      close(ionum)
    endif
    deallocate(pdl)
    call mpi_barrier(mpi_world,ierr)
    return
  end subroutine final_pdens
end module pdens
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd lib"
!     End:
