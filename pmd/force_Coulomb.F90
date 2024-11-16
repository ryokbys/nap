module Coulomb
!-----------------------------------------------------------------------
!                     Last modified: <2024-11-16 22:31:41 KOBAYASHI Ryo>
!-----------------------------------------------------------------------
!  Parallel implementation of Coulomb potential
!
!  For screened Coulomb potential:
!    - Ref. Adams & Rao, Phys. Status Solidi A 208, No.8 (2011)
!    - Cutoff treatment of energy and force shifts.
!  For Ewald Coulomb potential:
!    - ...
!-----------------------------------------------------------------------
  use pmdvars,only: nspmax,nsp
  use util,only: csp2isp, itotOf
  use memory,only: accum_mem
  use vector,only: dot,norm,cross
  implicit none
  include "./const.h"
  save

  character(len=128):: paramsdir = '.'
  character(len=128),parameter:: paramsfname = 'in.params.Coulomb'
  logical:: params_read = .false.
  real(8),parameter:: pi = 3.14159265398979d0

  logical:: lprmset_Coulomb = .false.

  character(len=128):: cchgs,cdist,cterms
!.....Input Keywords: charges, charge_dist, terms
!.....Keywords for charges: fixed, fixed_bvs, variable/qeq
!.....Keywords for charge_dist: point, gaussian
!.....Keywords for terms: full, short, long, direct, direct_cut, screened_cut

  integer,parameter:: ioprms = 20
!.....Coulomb's constant, acc = 1.0/(4*pi*epsilon0) in eV *Ang /e^2
  real(8),parameter:: acc  = 14.3998554737d0
!.....permittivity of vacuum
  real(8),parameter:: eps0 = 0.00552634939836d0  ! e^2 /Ang /eV


!!$  integer,parameter:: nspmax = 9
!!$  integer:: nsp = 1
!.....Flag for existence of the species
  logical:: ispflag(nspmax)
!.....Species charge
  real(8):: schg(nspmax),schg0(nspmax)
!  logical,allocatable:: interact(:,:)
!.....Interactions
  character(len=20):: cinteract = 'full'
  logical:: interact(nspmax,nspmax)
!.....ideal valence charges of species
  real(8):: vid_bvs(nspmax)
!  real(8),allocatable:: vid_bvs(:)
!.....principal quantum numbers of species
  integer:: npq_bvs(nspmax)
!  integer,allocatable:: npq_bvs(:)
!.....covalent radius
  real(8):: rad_bvs(nspmax)
!  real(8),allocatable:: rad_bvs(:)
!.....screening length
  real(8):: rho_scr(nspmax,nspmax)
!  real(8),allocatable:: rho_scr(:,:)
!!$  real(8):: fbvs = 0.74d0 +- 0.04
  real(8):: fbvs = 0.74d0
!.....Dielectric constant
  real(8):: dielec = 1d0
!.....Scaling factor multiplied to E_Coulomb
  real(8):: fscale = 1d0

!.....charge threshold for Coulomb interaction [default: 0.01]
  real(8),parameter:: qthd = 1d-12

!.....Gaussian width of Ewald sum
  real(8):: sgm_ew = 3.5355339d0
  real(8):: sgm(nspmax)
!.....Rho value for screened_cut
!     Default value = 5.0 (Ang), which corresponds to alpha = 0.2 A^{-1}
!     See, C.J. Fennell and J.D. Gezelter, J. Chem. Phys. 124, 234104 (2006).
  real(8):: rho_screened_cut = 5.0d0
  real(8):: vrcs(nspmax,nspmax),dvdrcs(nspmax,nspmax)
  real(8):: rcut = -1d0

!.....Accuracy controlling parameter for Ewald sum
!.....See, http://www.jncasr.ac.in/ccms/sbs2007/lecturenotes/5day10nov/SBS_Ewald.pdf
!.....Exp(-pacc) = 1e-7 when pacc= 18.0
!  real(8),parameter:: pacc   = 18d0
  real(8):: pacc = 9.21034d0  ! exp(-pacc) = 1e-4
!.....real-space cell volume
  real(8):: vol
!.....k-space variables
  real(8):: b1(3),b2(3),b3(3)
  real(8),allocatable:: qcosl(:),qcos(:),qsinl(:),qsin(:),pflr(:,:)
!.....Initial kmax = 20 is hard coded, which has no meaning.
  integer,parameter:: kmaxini = 20
  real(8):: bkmax
  integer:: kmax1,kmax2,kmax3,nk
  logical,allocatable:: lkuse(:,:,:)
!.....kmax threshold
  real(8),parameter:: threshold_kmax = 1d-4

!.....Variable-charge optimization method: damping, cg, matinv
  character(len=20):: chgopt_method = 'damping'
!.....Variable-charge potential variables
  real(8):: vc_chi(nspmax),vc_jii(nspmax),vc_e0(nspmax),vcg_sgm(nspmax) &
       ,qtop(nspmax),qbot(nspmax)
  real(8),parameter:: vcg_lambda = 0.5d0
!.....Convergence criterion for QEq in eV/atom
  real(8):: conv_eps_qeq = 1.0d-8
!.....chgopt-damping related parameters
  character(len=20):: codmp_method = 'damping' ! or FIRE or damping
  integer:: nstp_qeq = 100
  real(8):: dt_codmp = 0.005d0  ! fs
  real(8):: qmass = 0.002d0  ! atomic mass unit
  real(8):: qtot_qeq = 0d0
  integer:: minstp_qeq = 3
  integer:: minstp_conv_qeq = 3
!.....Velocity-damping related parameters
  real(8):: fdamp_codmp = 0.7d0
  real(8):: dfdamp_codmp = 0.9d0
!.....FIRE-related parameters
  real(8):: finc_codmp = 1.1d0
  real(8):: fdec_codmp = 0.5d0
  real(8):: alpha0_codmp = 0.1d0
  real(8):: falpha_codmp = 0.99d0
!.....Gradient descent parameters
  real(8):: dqmax_cogrd = 1d-1
  real(8):: dqeps_cogrd = 1d-4
!.....2nd order potential coeff for bounding q inside [qbot,qtop]
  real(8):: bound_k2 = 1.0d+2
!.....4-th order potential coeff for bounding q inside [qbot,qtop]
  real(8):: bound_k4 = 0d0
!.....Extended lagrangian
  real(8),allocatable:: aauxq(:)
  real(8):: omg2dt2 = 1d0  ! = omg^2*dt^2 = 2 is recommended by Nomura et al.
  real(8):: auxomg2 = 32.d0

  integer,parameter:: ivoigt(3,3)= &
       reshape((/ 1, 6, 5, 6, 2, 4, 5, 4, 3 /),shape(ivoigt))
contains
  subroutine init_coulomb(myid,mpi_md_world,iprint,lvc,specorder)
!
!  Allocate and initialize parameters to be used.
!  This is called when force is 'Coulomb'
!  The *_coulombx methods will replace *_coulomb methods in the future.
!
    include "mpif.h"
    integer,intent(in):: myid,mpi_md_world,iprint
    character(len=3),intent(in):: specorder(nspmax)
    logical,intent(inout):: lvc

    integer:: i,is,ierr,nspl

    rho_scr(1:nspmax,1:nspmax) = 0d0

!!$!.....Get umber of species
!!$    nsp = nspin

    call read_params(myid,mpi_md_world,iprint,specorder)

    if( trim(cchgs).eq.'variable' .or. trim(cchgs).eq.'qeq' ) then
      lvc = .true.
    endif
    
  end subroutine init_coulomb
!=======================================================================
  subroutine init_for_Ewald(h,rc,myid,mpi_world,iprint,l1st)
!
! Initialization for Ewald method, which should be called in both cases of
! variable charge and fixed charge.
!
    real(8),intent(in):: h(3,3),rc
    integer,intent(in):: myid,mpi_world,iprint
    logical,intent(in):: l1st
    
    if( trim(cterms).eq.'full' .or. trim(cterms).eq.'long' ) then
      if( trim(cchgs).eq.'variable' .or. trim(cchgs).eq.'qeq' ) then
        call init_vc_Ewald(h,rc,myid,mpi_world,iprint,l1st)
      else
        call init_fc_Ewald(h,rc,myid,mpi_world,iprint,l1st)
      endif
    endif
    
  end subroutine init_for_Ewald
!=======================================================================
  subroutine init_fc_Ewald(h,rc,myid,mpi_world,iprint,l1st)
!
!  Ewald sum with fixed charge.
!
    implicit none 
    integer,intent(in):: myid,mpi_world,iprint
    real(8),intent(in):: h(3,3),rc
    logical,intent(in):: l1st 

    integer:: i,ik,k1,k2,k3,isp
    real(8):: bk1(3),bk2(3),bk3(3),bk(3),bb2

    sgm_ew = rc/sqrt(2d0*pacc)
    sgm(1:nspmax) = sgm_ew
    bkmax  = 2d0*pacc /rc
    if( myid.eq.0 .and. iprint.ge.ipl_basic .and. l1st ) then
      write(6,'(/,a)') ' Ewald sum parameters:'
      write(6,'(a,f12.4)') '   1/(4*pi*eps0)        = ', acc
      write(6,'(a,f12.4)') '   Accuracy parameter p = ', pacc
      write(6,'(a,f12.4)') '   Gaussian width sgm   = ', sgm_ew
      write(6,'(a,f12.4)') '   real-space cutoff    = ', rc
    endif

    if( .not. (trim(cterms).eq.'full' .or. trim(cterms).eq.'long') ) return
    call get_recip_vectors(h)
    if( myid.eq.0 .and. iprint.ge.ipl_basic .and. l1st ) then
      write(6,'(a,f12.4)') '   k-space cutoff       = ', bkmax
      write(6,'(a)') ' Reciprocal vectors:'
      write(6,'(a,3es12.3)') '   b1 = ',b1(1:3)
      write(6,'(a,3es12.3)') '   b2 = ',b2(1:3)
      write(6,'(a,3es12.3)') '   b3 = ',b3(1:3)
    endif
!.....kmax# is constant during MD run even if h-matrix can change...
    call setup_kspace()
    if( myid.eq.0 .and. iprint.ge.ipl_basic .and. l1st ) then
      write(6,'(a)') ' Number of k-points for Ewald sum:'
      write(6,'(a,i0)') '   kmax1 = ',kmax1
      write(6,'(a,i0)') '   kmax2 = ',kmax2
      write(6,'(a,i0)') '   kmax3 = ',kmax3
      write(6,'(a,i0)') '   total = ',nk
    endif
!!$    print *,'myid,lkuse(:)=',myid,lkuse(:,:,:)
!!$    print *,'myid,b1(1:3)=',myid,b1(1:3)
!!$    print *,'myid,b2(1:3)=',myid,b2(1:3)
!!$    print *,'myid,b3(1:3)=',myid,b3(1:3)
    if( allocated(qcos) ) then
      call accum_mem('force_Coulomb',-8*(size(qcosl)+size(qcos)+size(qsinl) &
           +size(qsin)+size(pflr)))
      deallocate(qcosl,qcos,qsinl,qsin,pflr)
    endif
    allocate(qcosl(nk),qcos(nk),qsinl(nk),qsin(nk),pflr(nk,nspmax))
    call accum_mem('force_Coulomb',8*(size(qcosl)+size(qcos)+size(qsinl) &
         +size(qsin)+size(pflr)))
!.....prefactor for long-range term
    ik = 0
    do k1= -kmax1,kmax1
      bk1(1:3) = k1*b1(1:3)
      do k2= -kmax2,kmax2
        bk2(1:3) = k2*b2(1:3)
        do k3= -kmax3,kmax3
          if( .not. lkuse(k3,k2,k1) ) cycle
          ik= ik+1
          bk3(1:3) = k3*b3(1:3)
          bk(1:3) = bk1(1:3) +bk2(1:3) +bk3(1:3)
          bb2 = norm(bk)
          bb2 = bb2*bb2
          pflr(ik,1:nspmax)= 4d0 *pi /bb2 *exp(-0.5d0 *sgm_ew**2 *bb2)
        enddo
      enddo
    enddo

  end subroutine init_fc_Ewald
!=======================================================================
  subroutine init_vc_Ewald(h,rc,myid,mpi_world,iprint,l1st)
!
!  Since variable-charge potential with Gaussian distribution charges
!  is mostly identical to the long-range part of Ewald summation,
!  the code here is mostly the same as Ewald sum except that vcGaussian
!  requires differnt Gaussian width for different species that should be
!  read from input file in.params.Coulomb.
!
    implicit none
    real(8),intent(in):: h(3,3),rc
    integer,intent(in):: myid,mpi_world,iprint
    logical,intent(in):: l1st

    integer:: i,isp,ik,k1,k2,k3,is
    real(8):: bk1(3),bk2(3),bk3(3),bk(3),bb2,sgm_min,sgm_rcmd

!.....Current implementation does not allow species-dpendent sigma.
    vcg_sgm(1:nspmax) = sgm_ew
    sgm(1:nspmax) = sgm_ew
!.....If long-range term exists, self interaction should be added to Jii.
    if( trim(cterms).eq.'full' .or. trim(cterms).eq.'long' ) then
      do is=1,nspmax
        vc_jii(is) = vc_jii(is) -acc*sqrt(2d0/pi) /vcg_sgm(is)
      enddo
    endif
!.....Detect minimum sigma to determine k-max
    sgm_min = sgm_ew
    bkmax  = sqrt(2d0*pacc) /sgm_min
    if( myid.eq.0 .and. iprint.ne.0 .and. l1st) then
      write(6,'(/,a)') ' Ewald sum parameters:'
      write(6,'(a,f12.4)') '   1/(4*pi*eps0)        = ', acc
      write(6,'(a,f12.4)') '   Accuracy parameter p = ', pacc
      write(6,'(a,f12.4)') '   Gaussian simga       = ', sgm_ew
      sgm_rcmd = rc/sqrt(2d0*pacc)
      write(6,'(a,f12.4)') '   Recommended sigma    = ', sgm_rcmd
      write(6,'(a,f12.4)') '   real-space cutoff    = ', rc
      write(6,'(a,f12.4)') '   k-space cutoff       = ', bkmax
    endif
    call get_recip_vectors(h)
    if( myid.eq.0 .and. iprint.ne.0 .and. l1st ) then
      write(6,'(a)') ' Reciprocal vectors:'
      write(6,'(a,3es12.3)') '   b1 = ',b1(1:3)
      write(6,'(a,3es12.3)') '   b2 = ',b2(1:3)
      write(6,'(a,3es12.3)') '   b3 = ',b3(1:3)
    endif
!.....kmax# is constant during MD run even if h-matrix can change...
    call setup_kspace()
    if( myid.eq.0 .and. iprint.ne.0 .and. l1st ) then
      write(6,'(a)') ' Number of k-points for Ewald sum:'
      write(6,'(a,i8)') '   kmax1 = ',kmax1
      write(6,'(a,i8)') '   kmax2 = ',kmax2
      write(6,'(a,i8)') '   kmax3 = ',kmax3
      write(6,'(a,i8)') '   total = ',nk
      write(6,*) ''
    endif
    if( allocated(qcos) ) then
      call accum_mem('force_Coulomb',-8*(size(qcosl)+size(qcos)+size(qsinl) &
           +size(qsin)+size(pflr)))
      deallocate(qcos,qsin,qcosl,qsinl,pflr)
    endif
    allocate(qcosl(nk),qcos(nk),qsinl(nk),qsin(nk),pflr(nk,nspmax))
    call accum_mem('force_Coulomb',8*(size(qcosl)+size(qcos)+size(qsinl) &
         +size(qsin)+size(pflr)))
!.....prefactor for long-range term
    ik = 0
    do k1= -kmax1,kmax1
      bk1(1:3) = k1*b1(1:3)
      do k2= -kmax2,kmax2
        bk2(1:3) = k2*b2(1:3)
        do k3= -kmax3,kmax3
          if( .not. lkuse(k3,k2,k1) ) cycle
          ik= ik+1
          bk3(1:3) = k3*b3(1:3)
          bk(1:3) = bk1(1:3) +bk2(1:3) +bk3(1:3)
          bb2 = norm(bk)
          bb2 = bb2*bb2
          do isp=1,nsp
            pflr(ik,isp)= 4d0 *pi /bb2 *exp(-0.5d0 *sgm(isp)**2 *bb2)
          enddo
        enddo
      enddo
    enddo

  end subroutine init_vc_Ewald
!=======================================================================
  subroutine read_params(myid,mpi_world,iprint,specorder)
!
!  Read parameter file for any Coulomb potential.
!
    use util, only: num_data
    include "mpif.h"
    integer,intent(in):: myid,mpi_world,iprint
    character(len=3),intent(in):: specorder(nspmax)

    character(len=128):: cmode,cline,ctmp,fname,c1st
    character(len=3):: cname,csp,cspj,cspi
    integer:: i,ierr,jerr,isp,jsp,npq,nentry
    real(8):: chgi,vid,rad,dchi,djii,sgmt,de0,qlow,qup&
         ,vcgjiimin,sgmlim, rin,rout,rhoij, rhoii
    logical:: rhoii_set, rhoij_set, rad_set

    if( myid.eq.0 ) then
!.....Initialization
      interact(1:nspmax,1:nspmax) = .true.
      ispflag(1:nspmax) = .false.
      vc_chi(1:nspmax) = 0d0
      vc_jii(1:nspmax) = 0d0
      vc_e0(1:nspmax) = 0d0
      vcg_sgm(1:nspmax) = 0d0
      qbot(1:nspmax) = 0d0
      qtop(1:nspmax) = 0d0
      rho_scr(:,:) = -1d0
      rad_bvs(:) = -1d0
      cmode = 'none'
!.....File name
      fname = trim(paramsdir)//'/'//trim(paramsfname)
      open(ioprms,file=trim(fname),status='old')
      if( iprint.ge.ipl_basic ) write(6,'(/,a)') ' Coulomb parameters:'
!.....Start reading
      do while(.true.)
        read(ioprms,'(a)',end=10) cline
        nentry = num_data(cline,' ')
        if( nentry.eq.0 ) then
          cmode = 'none'
          cycle
        endif
        if( cline(1:1).eq.'!' .or. cline(1:1).eq.'#' ) cycle
        read(cline,*) c1st
!.....Detect keyword....................................................
!.....But 1st, if a certain mode is already selected...
        if( trim(cmode).eq.'charges' ) then
          backspace(ioprms)
          if( trim(cchgs).eq.'fixed' ) then
!!$            read(ioprms,*) isp, chgi
            read(ioprms,*) csp, chgi
            isp = csp2isp(trim(csp))
            if( isp.gt.0 ) then
              schg0(isp) = chgi
              schg(isp) = schg0(isp)
              ispflag(isp) = .true.
              if( iprint.ge.ipl_basic ) print '(a,a3,i3,f8.4)','   fixed charge: ',trim(csp),isp,chgi
            else
              if( iprint.ge.ipl_info ) then
                print '(a,a3,i3,f8.4)','   fixed charge read but not used: ',trim(csp),isp,chgi
              endif
            end if
          else if( trim(cchgs).eq.'fixed_bvs' ) then
            read(ioprms,*) csp,vid,rad,npq
            isp = csp2isp(trim(csp))
            if( isp.lt.1 .or. isp.gt.nspmax ) then
              if( iprint.ge.ipl_info ) then
                print *,'  fixed_bvs charge read but not used: ',trim(csp)
              endif
              cycle
            endif
            ispflag(isp) = .true.
            vid_bvs(isp) = vid
            rad_bvs(isp) = rad
            npq_bvs(isp) = npq
            if( iprint.ge.ipl_basic ) then
              write(6,'(a,a5,2f7.3,i4)') '   csp,vid,rad,npq =' &
                   ,trim(csp),vid,rad,npq
            endif
          else if( trim(cchgs).eq.'variable' .or. trim(cchgs).eq.'qeq') then
            read(ioprms,*) csp, dchi,djii,de0,qlow,qup
            isp = csp2isp(trim(csp))
            if( isp.lt.1 .or. isp.gt.nspmax ) then
              if( iprint.ge.ipl_info ) then
                print *,'  variable charge read but not used: ',trim(csp)
              endif
              cycle
            endif
            ispflag(isp) = .true.
            vc_chi(isp) = dchi
            vc_jii(isp) = djii
            vc_e0(isp) = de0
            qbot(isp) = qlow
            qtop(isp) = qup
            if( iprint.ge.ipl_basic ) then
              write(6,'(a,a3,4f10.4,2f5.1)') &
                   '   csp,chi,Jii,e0,qbot,qtop = ', &
                   trim(csp),dchi,djii,de0,qlow,qup
            endif
          endif
!!$        else if( trim(cmode).eq.'charge_dist' ) then
!!$        else if( trim(cmode).eq.'terms' ) then
        else if( trim(cmode).eq.'interactions' ) then
          backspace(ioprms)
!!$          read(ioprms,*) isp,jsp
!!$          interact(isp,jsp) = .true.
!!$          interact(jsp,isp) = .true.
          read(ioprms,*) cspi,cspj
          isp = csp2isp(trim(cspi))
          jsp = csp2isp(trim(cspj))
          if( isp.gt.0 .and. isp.le.nspmax .and. jsp.gt.0 .and. jsp.le.nspmax ) then
            interact(isp,jsp) = .true.
            interact(jsp,isp) = .true.
            if( iprint.ge.ipl_basic ) print '(5x,a)',trim(cspi)//'-'//trim(cspj)
          else
            if( iprint.ge.ipl_info ) print *,'  interacion read but not used: ',isp,jsp
          endif

!.....If no cmode is given, detect keyword.
        else if( trim(c1st).eq.'charges' ) then
          cmode = 'charges'
          backspace(ioprms)
          read(ioprms,*) ctmp, cchgs
          if(  trim(cchgs).ne.'fixed' .and. &
               trim(cchgs).ne.'fixed_bvs' .and. &
               trim(cchgs).ne.'variable' .and. &
               trim(cchgs).ne.'qeq' ) then
            print *,'ERROR: charges should have an argument of'//&
                 ' either fixed, fixed_bvs, variable or qeq.'
            stop
          endif
          if( iprint.ge.ipl_info ) print *,trim(ctmp),' '//trim(cchgs)
          cycle
        else if( trim(c1st).eq.'charge_dist' ) then
          cmode = 'charge_dist'
          backspace(ioprms)
          read(ioprms,*) ctmp, cdist
          if(  trim(cdist).ne.'point' .and. &
               trim(cdist).ne.'gaussian' ) then
            print *,'ERROR: charge_dist should have an argument of'//&
                 ' either point or gaussian.'
            stop
          endif
          if( iprint.ge.ipl_info ) print *,trim(ctmp),' '//trim(cdist)
          cycle
        else if( trim(c1st).eq.'terms' ) then
          cmode = 'terms'
          backspace(ioprms)
          read(ioprms,*) ctmp, cterms
          if(  trim(cterms).ne.'full' .and. &
               trim(cterms).ne.'short' .and. &
               trim(cterms).ne.'long' .and. &
               trim(cterms).ne.'direct' .and. &
               trim(cterms).ne.'direct_cut' .and. &
               trim(cterms).ne.'screened_cut' ) then
            print *,'ERROR: terms should have an argument of '//&
                 'either direct_cut, screened_cut, full, short, or long.'
            stop
          endif
          if( iprint.ge.ipl_info ) print *,trim(ctmp),' '//trim(cterms)
          cycle
        else if( trim(c1st).eq.'interactions' ) then
          backspace(ioprms)
          read(ioprms,'(a)') ctmp
          nentry = num_data(ctmp,' ')
          if( nentry.eq.1 ) then
            cmode = 'interactions'
            cinteract = 'given'
            interact(1:nspmax,1:nspmax) = .false.
            if( iprint.ge.ipl_basic ) print *,'  selected interactions:'
          else if( nentry.eq.2 ) then
            backspace(ioprms)
            read(ioprms,*) ctmp, cinteract
            cmode = 'none'
          endif
          cycle
        else if( trim(c1st).eq.'rcut' ) then
          backspace(ioprms)
          read(ioprms,*) ctmp, rcut
          cycle
        else if( trim(c1st).eq.'dielectric' ) then
          backspace(ioprms)
          read(ioprms,*) ctmp, dielec
          cycle
        else if( trim(c1st).eq.'scale_factor' ) then
          backspace(ioprms)
          read(ioprms,*) ctmp, fscale
          cycle
        else if( trim(c1st).eq.'sigma' ) then
          backspace(ioprms)
          read(ioprms,*) ctmp, sgm_ew
          cycle
        else if( trim(c1st).eq.'fbvs' ) then
          backspace(ioprms)
          read(ioprms,*) ctmp, fbvs
          cycle
        else if( trim(c1st).eq.'pacc' ) then
          backspace(ioprms)
          read(ioprms,*) ctmp, pacc
          cycle
        else if( trim(c1st).eq.'rho_screened_cut' ) then
          backspace(ioprms)
          read(ioprms,*) ctmp, rho_screened_cut
          if( iprint.ge.ipl_info ) print *,trim(ctmp),rho_screened_cut
          cycle
        else if( trim(c1st).eq.'rad_screened_cut' ) then
!.....Set rho_screened_cut minus to show rad should be used to determine rho_screened_cut
          rho_screened_cut = -1d0 *abs(rho_screened_cut)
          backspace(ioprms)
          read(ioprms,*) ctmp, csp, rad
          if( iprint.ge.ipl_info ) print '(a,3x,a,1x,f7.4)',trim(ctmp), trim(csp), rad
          isp = csp2isp(trim(csp))
          if( isp.gt.0 ) then
            rad_bvs(isp) = rad
          endif
          cycle
        else if( trim(c1st).eq.'rhoii_screened_cut' ) then
          backspace(ioprms)
          read(ioprms,*) ctmp, cspi, rhoii
          isp = csp2isp(trim(cspi))
          if( isp.le.0 ) cycle
          rho_scr(isp,isp) = rhoii
          if( iprint.ge.ipl_info ) print '(1x,a,3x,a,1x,f7.3)',trim(ctmp),trim(cspi),rhoii
          cycle
        else if( trim(c1st).eq.'rhoij_screened_cut' ) then
          backspace(ioprms)
          read(ioprms,*) ctmp, cspi, cspj, rhoij
          isp = csp2isp(trim(cspi))
          jsp = csp2isp(trim(cspj))
          if( isp.le.0 .or. jsp.le.0 ) cycle
          rho_scr(isp,jsp) = rhoij
          rho_scr(jsp,isp) = rhoij
          if( iprint.ge.ipl_info ) print *,trim(ctmp),trim(cspi) &
               ,trim(cspj),rhoij
          cycle
!.....QEq related parameters hereafter
        else if( trim(c1st).eq.'chgopt_method' ) then
          backspace(ioprms)
          read(ioprms,*) ctmp, chgopt_method
          cycle
        else if( trim(c1st).eq.'codmp_method' ) then
          backspace(ioprms)
          read(ioprms,*) ctmp, codmp_method
          cycle
        else if( trim(c1st).eq.'conv_eps_qeq' ) then
          backspace(ioprms)
          read(ioprms,*) ctmp, conv_eps_qeq
          cycle
        else if( trim(c1st).eq.'nstp_qeq') then
          backspace(ioprms)
          read(ioprms,*) ctmp, nstp_qeq
          cycle
        else if( trim(c1st).eq.'minstp_qeq') then
          backspace(ioprms)
          read(ioprms,*) ctmp, minstp_qeq
          cycle
        else if( trim(c1st).eq.'minstp_conv_qeq') then
          backspace(ioprms)
          read(ioprms,*) ctmp, minstp_conv_qeq
          cycle
        else if( trim(c1st).eq.'dt_codmp') then
          backspace(ioprms)
          read(ioprms,*) ctmp, dt_codmp
          cycle
        else if( trim(c1st).eq.'qmass') then
          backspace(ioprms)
          read(ioprms,*) ctmp, qmass
          cycle
        else if( trim(c1st).eq.'fdamp_codmp') then
          backspace(ioprms)
          read(ioprms,*) ctmp, fdamp_codmp
          cycle
        else if( trim(c1st).eq.'dfdamp_codmp') then
          backspace(ioprms)
          read(ioprms,*) ctmp, dfdamp_codmp
          cycle
        else if( trim(c1st).eq.'finc_codmp') then
          backspace(ioprms)
          read(ioprms,*) ctmp, finc_codmp
          cycle
        else if( trim(c1st).eq.'fdec_codmp') then
          backspace(ioprms)
          read(ioprms,*) ctmp, fdec_codmp
          cycle
        else if( trim(c1st).eq.'alpha0_codmp') then
          backspace(ioprms)
          read(ioprms,*) ctmp, alpha0_codmp
          cycle
        else if( trim(c1st).eq.'falpha_codmp') then
          backspace(ioprms)
          read(ioprms,*) ctmp, falpha_codmp
          cycle
        else if( trim(c1st).eq.'dqmax_cogrd') then
          backspace(ioprms)
          read(ioprms,*) ctmp, dqmax_cogrd
          cycle
        else if( trim(c1st).eq.'dqeps_cogrd') then
          backspace(ioprms)
          read(ioprms,*) ctmp, dqeps_cogrd
          cycle
        else if( trim(c1st).eq.'bound_k4') then
          backspace(ioprms)
          read(ioprms,*) ctmp, bound_k4
          cycle
        else if( trim(c1st).eq.'bound_k2') then
          backspace(ioprms)
          read(ioprms,*) ctmp, bound_k2
          cycle
        else if( trim(c1st).eq.'omg2dt2') then
          backspace(ioprms)
          read(ioprms,*) ctmp, omg2dt2
          cycle
        else
          print *,' No such entry: ',trim(c1st)
          cycle
        endif
      enddo ! while(.true.)

!.....Corrections
      if( trim(cdist).eq.'gaussian' ) then
        if(  trim(cterms).eq.'full' .or. &
             trim(cterms).eq.'short' ) then
          cterms = 'long'
          print *,'  WARNING: terms was corrected to long, because charge_dist is gaussian.'
        endif
      endif

10    close(ioprms)

!.....In case of variable charge, sigma has a lower bound w.r.t. min(Jii)
      if( (trim(cchgs).eq.'variable' .or. trim(cchgs).eq.'qeq') &
           .and. .not. trim(cterms).ne.'screened_cut' &
           .and. .not. trim(cterms).ne.'direct_cut' ) then
        vcgjiimin = 1d+30
        do isp=1,nsp
          vcgjiimin = min(vcgjiimin,vc_jii(isp))
        enddo
        sgmlim = acc*sqrt(2d0/pi)/vcgjiimin
        if( sgm_ew.lt.sgmlim ) then
          if( iprint.ne.0 ) then
            print *,'  WARNING: Since sgm_ew is too small, sgm_ew is replaced by sgmlim'
            print *,'           which is determined by acc*sqrt(2/pi)/Jii.'
            print *,'     original sgm_ew = ', sgm_ew
            print *,'     replaced sgm_ew = ', sgmlim
          endif
          sgm_ew = sgmlim
        endif
      endif

!.....Set screening length in the case of screened_cut
      if( trim(cterms).eq.'screened_cut' ) then
        rhoii_set = .true.
        rhoij_set = .true.
        rad_set = .true.
        do isp=1,nsp
          if( rho_scr(isp,isp).lt.0d0 ) rhoii_set = .false.
          if( rad_bvs(isp).lt.0d0 ) rad_set = .false.
          do jsp=1,nsp
            if( isp.eq.jsp ) cycle
            if( rho_scr(isp,jsp).lt.0d0 ) rhoij_set = .false.
          enddo
        enddo

        if( trim(cchgs).eq.'fixed_bvs' ) then
          if( iprint.ge.ipl_basic ) print '(a)', ' rhoij are set from given rads of fixed_bvs entries.'
          do isp=1,nsp
            do jsp=1,nsp
              rho_scr(isp,jsp) = fbvs*(rad_bvs(isp)+rad_bvs(jsp))
            enddo
          enddo
        else
          if( rad_set ) then
            if( iprint.ge.ipl_basic ) print '(a)', ' rhoij are set from given rad info.'
            do isp=1,nsp
              do jsp=1,nsp
                rho_scr(isp,jsp) = rad_bvs(isp)+rad_bvs(jsp)
              enddo
            enddo
          else if( .not.rhoij_set .and. rhoii_set ) then
            if( iprint.ge.ipl_basic ) print '(a)', ' rhoij are set from given rhoii info as ' &
                 //'rhoij=(rhoi + rhoj)/2.'
            do isp=1,nsp
              do jsp=1,nsp
                if( isp.eq.jsp ) cycle
                rho_scr(isp,jsp) = (rho_scr(isp,isp) +rho_scr(jsp,jsp)) /2
              enddo
            enddo
          else
            if( iprint.ge.ipl_basic ) print '(a,f7.3,a)', ' rhoij are set from a unique rho (=', &
                 rho_screened_cut, ')'
            do isp=1,nsp
              do jsp=1,nsp
                if( rho_scr(isp,jsp).lt.0d0 ) rho_scr(isp,jsp) = rho_screened_cut
              enddo
            enddo
          endif
        endif
      endif

    endif  ! myid.eq.0

!.....Broadcast data just read from the file
    call mpi_bcast(cterms,128,mpi_character,0,mpi_world,ierr)
    call mpi_bcast(cdist,128,mpi_character,0,mpi_world,ierr)
    call mpi_bcast(cchgs,128,mpi_character,0,mpi_world,ierr)

    call mpi_bcast(schg0,nspmax,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(schg,nspmax,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(vc_chi,nsp,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(vc_jii,nsp,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(vc_e0,nsp,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(rcut,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(dielec,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(fscale,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(sgm_ew,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(qbot,nspmax,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(qtop,nspmax,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(vid_bvs,nspmax,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(rad_bvs,nspmax,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(npq_bvs,nspmax,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(cinteract,20,mpi_character,0,mpi_world,ierr)
    call mpi_bcast(interact,nspmax*nspmax,mpi_logical,0,mpi_world,ierr)
    call mpi_bcast(ispflag,nspmax,mpi_logical,0,mpi_world,ierr)

    call mpi_bcast(pacc,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(fbvs,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(rho_screened_cut,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(rho_scr,nspmax*nspmax,mpi_real8,0,mpi_world,ierr)

    call mpi_bcast(chgopt_method,20,mpi_character,0,mpi_world,ierr)
    call mpi_bcast(codmp_method,20,mpi_character,0,mpi_world,ierr)
    call mpi_bcast(conv_eps_qeq,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(nstp_qeq,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(minstp_qeq,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(minstp_conv_qeq,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(dt_codmp,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(qmass,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(qtot_qeq,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(fdamp_codmp,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(dfdamp_codmp,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(finc_codmp,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(fdec_codmp,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(alpha0_codmp,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(falpha_codmp,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(dqmax_cogrd,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(dqeps_cogrd,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(bound_k2,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(bound_k4,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(omg2dt2,1,mpi_real8,0,mpi_world,ierr)

    if( trim(cterms).eq.'screened_cut' .or. trim(cterms).eq.'short' ) then
      if( myid.eq.0 .and. iprint.ge.ipl_basic ) then
        do isp=1,nsp
          cspi = specorder(isp)
          do jsp=isp,nsp
!!$            rho_scr(isp,jsp) = fbvs*(rad_bvs(isp)+rad_bvs(jsp))
!!$            rho_scr(isp,jsp) = 2d0
            if( interact(isp,jsp) ) then
              cspj = specorder(jsp)
              write(6,'(a,2a5,f10.4)') '   cspi,cspj,rho_scr= ',trim(cspi),trim(cspj),rho_scr(isp,jsp)
            endif
          enddo
        enddo
      endif
    endif

    if( (trim(cchgs).eq.'variable' .or. trim(cchgs).eq.'qeq') &
         .and. trim(cinteract).eq.'same_sign' ) then
      if( myid.eq.0 ) print *,'ERROR@read_params: variable/qeq '&
           //'cannot be used with same_sign interactions.'
      call mpi_finalize(ierr)
      stop 1
    endif

    if( myid.eq.0 .and. iprint.gt.ipl_basic ) then
      write(6,'(a)') ' Finished reading '//trim(fname)
    endif
!!$    params_read = .true.

    return
  end subroutine read_params
!=======================================================================
  subroutine force_Coulomb(namax,natm,tag,ra,nnmax,aa,strs &
       ,chg,h,hi,nb,nbmax,lsb,nex,lsrc &
       ,myparity,nn,sv,rc,lspr,sorg &
       ,mpi_md_world,myid,epi,epot,nismax,lstrs,iprint &
       ,l1st,lcell_updated,lvc,specorder)
!
!  Coulomb potential and force computation using Ewald sum method.
!  Currently, the efficiency when parallel computation is not taken into account.
!  So it is not appropriate to use this routine for large system.
!
    implicit none
    include "mpif.h"
    include "./params_unit.h"

    integer,intent(in):: namax,natm,nnmax,nismax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),lspr(0:nnmax,namax),nex(3)
    integer,intent(in):: mpi_md_world,myid
    real(8),intent(in):: ra(3,namax),h(3,3),hi(3,3),rc &
         ,tag(namax),sv(3,6),sorg(3)
    real(8),intent(inout):: chg(namax)
    real(8),intent(inout):: aa(3,namax),epi(namax),strs(3,3,namax),epot
    logical,intent(in):: lstrs,l1st,lcell_updated,lvc
    character(len=3),intent(in):: specorder(nspmax)

    integer:: i,j,ik,is,js,k1,k2,k3,ierr,jj,ixyz,jxyz
    real(8):: elrl,elr,esrl,esr,epotl,epott,qi,qj,tmp,ftmp &
         ,bdotr,terfc,diji,dij,ss2i,sgmsq2,rc2,q2tot,q2loc,bb2 &
         ,e0,q2,sgmi
    real(8),save:: eself,eselfl
    real(8),allocatable,save:: ri(:),bk(:),bk1(:),bk2(:),bk3(:)&
         ,bb(:),dxdi(:),dxdj(:),rij(:),xij(:),xj(:),xi(:)
    real(8),allocatable,save:: strsl(:,:,:),aal(:,:)

    if( l1st ) then
      if( .not.allocated(ri) ) then
        allocate(ri(3),bk(3),bk1(3),bk2(3),bk3(3),bb(3),dxdi(3) &
             ,dxdj(3),rij(3),xij(3),xj(3),xi(3))
      endif
      if( allocated(strsl) ) then
        call accum_mem('force_Coulomb',-8*size(strsl)-8*size(aal))
        deallocate(strsl,aal)
      endif
      allocate(strsl(3,3,namax),aal(3,namax))
      call accum_mem('force_Coulomb',8*size(strsl)+8*size(aal))

      if( trim(cchgs).eq.'fixed_bvs' ) then
        call set_charge_BVS(natm,nb,tag,chg,myid,mpi_md_world,iprint,specorder)
      endif
!.....In case that Coulomb interactions work only between species having
!     charges of the same sign, set interaction array after setting BVS charges.
      if( trim(cinteract).eq.'same_sign' ) then
        interact(:,:) = .false.
        do is=1,nsp
          do js=1,nsp
            if( schg(is)*schg(js).gt.0d0 ) interact(is,js) = .true.
          enddo
        enddo
      endif
    endif

    if( .not.allocated(strsl) ) then
      allocate(strsl(3,3,namax),aal(3,namax))
      call accum_mem('force_Coulomb',8*size(strsl)+8*size(aal))
    else if( size(strsl).lt.3*3*namax ) then
      call accum_mem('force_Coulomb',-8*size(strsl)-8*size(aal))
      deallocate(strsl,aal)
      allocate(strsl(3,3,namax),aal(3,namax))
      call accum_mem('force_Coulomb',8*size(strsl)+8*size(aal))
    endif

    if( trim(cchgs).eq.'fixed' ) then
      do i=1,natm+nb
        is = int(tag(i))
        chg(i) = schg(is)
      enddo
    endif

    strsl(1:3,1:3,1:namax) = 0d0
    aal(:,1:namax) = 0d0
    elrl = 0d0
    esrl = 0d0
    eselfl = 0d0
    elr = 0d0
    esr = 0d0
    eself = 0d0

    if( lvc .or. trim(cterms).eq.'full' .or. trim(cterms).eq.'long' ) then
      call self_term(namax,natm,tag,chg,epi,eselfl,iprint,lvc)
      call mpi_allreduce(eselfl,eself,1,mpi_real8,mpi_sum,mpi_md_world,ierr)
    endif

    if(  trim(cterms).eq.'full' .or. &
         trim(cterms).eq.'short' ) then
      call Ewald_short(namax,natm,tag,ra,nnmax,aal,strsl,chg,h,hi &
           ,lspr,epi,esrl,iprint,lstrs,rc)
      call mpi_allreduce(esrl,esr,1,mpi_real8,mpi_sum,mpi_md_world,ierr)
    endif

    if(  trim(cterms).eq.'full' .or. &
         trim(cterms).eq.'long' ) then
      call Ewald_long(namax,natm,tag,ra,nnmax,aal,strsl,chg,h,hi &
           ,lspr,sorg,epi,elrl,iprint,myid,mpi_md_world,lstrs &
           ,lcell_updated)
      call mpi_allreduce(elrl,elr,1,mpi_real8,mpi_sum,mpi_md_world,ierr)
    endif

    if( trim(cterms).eq.'direct' ) then
      call force_direct(namax,natm,tag,ra,nnmax,aal,strsl,chg,h,hi &
           ,lspr,epi,esrl,iprint,lstrs,rc)
      call mpi_allreduce(esrl,esr,1,mpi_real8,mpi_sum,mpi_md_world,ierr)
    endif
    if( trim(cterms).eq.'direct_cut' ) then
      call force_direct_cut(namax,natm,tag,ra,nnmax,aal,strsl,chg,h,hi &
           ,lspr,epi,esrl,iprint,lstrs,rc,l1st)
      call mpi_allreduce(esrl,esr,1,mpi_real8,mpi_sum,mpi_md_world,ierr)
    endif
    if( trim(cterms).eq.'screened_cut' ) then
      call force_screened_cut(namax,natm,tag,ra,nnmax,aal,strsl,chg,h,hi &
           ,lspr,epi,esrl,iprint,lstrs,rc,l1st)
      call mpi_allreduce(esrl,esr,1,mpi_real8,mpi_sum,mpi_md_world,ierr)
    endif

!.....Scale energies, forces, stresses
    elr = elr *fscale
    esr = esr *fscale
    eself = eself *fscale
    aal(:,1:natm) = aal(:,1:natm) *fscale
    strsl(:,:,1:natm) = strsl(:,:,1:natm) *fscale

!.....Copy local variables to global ones
    aa(1:3,1:natm) = aa(1:3,1:natm) +aal(1:3,1:natm)
    strs(1:3,1:3,1:natm)= strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)

    if( l1st .and. myid.eq.0 .and. iprint.ge.ipl_basic ) then
      if( cterms(1:6).eq.'direct' ) then
        print *,''
        print '(a,f12.4)',' Direct Coulomb energy = ',esr
      else
        print *,''
        print *,'Ewald energy by terms:'
        print '(a,f12.4," eV")','   Self term         = ',eself
        print '(a,f12.4," eV")','   Short-range term  = ',esr
        print '(a,f12.4," eV")','   Long-range term   = ',elr
      endif
    endif

    epott = esr +elr +eself
    epot= epot +epott
    if( myid.eq.0 .and. iprint.ge.ipl_info ) &
         print *,'epot Coulomb,self,short,long = ',epott,eself,esr,elr

  end subroutine force_Coulomb
!=======================================================================
  subroutine force_direct(namax,natm,tag,ra,nnmax,aa,strsl,chg,h,hi &
       ,lspr,epi,esrl,iprint,lstrs,rc)
!
!  Direct Coulomb interaction with cutoff radius without any cutoff treatment.
!
    implicit none
    integer,intent(in):: namax,natm,nnmax,iprint, &
         lspr(0:nnmax,namax)
    real(8),intent(in)::tag(namax),ra(3,namax),chg(namax), &
         h(3,3),hi(3,3),rc
    logical,intent(in):: lstrs
    real(8),intent(inout):: aa(3,namax),strsl(3,3,namax), &
         epi(namax),esrl

    integer:: i,j,jj,is,js,ixyz,jxyz
    real(8):: rc2,ss2i,sqpi,qi,qj,dij,diji,tmp,ftmp
    real(8):: xi(3),xj(3),xij(3),rij(3),dxdi(3),dxdj(3)

!.....Compute direct sum
    rc2 = rc*rc
    esrl = 0d0
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is= int(tag(i))
      qi = chg(i)
      do jj=1,lspr(0,i)
        j = lspr(jj,i)
        if( j.eq.0 ) exit
        if( j.le.i ) cycle
        js = int(tag(j))
        qj = chg(j)
        xj(1:3) = ra(1:3,j)
        xij(1:3)= xj(1:3)-xi(1:3)
        rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij = rij(1)**2 +rij(2)**2 +rij(3)**2
        if( dij.gt.rc2 ) cycle
        dij = sqrt(dij)
        diji = 1d0/dij
        dxdi(1:3)= -rij(1:3)*diji
        dxdj(1:3)=  rij(1:3)*diji
!.....potential
        tmp = 0.5d0 *acc *qi*qj*diji /dielec
        if( j.le.natm ) then
          epi(i)= epi(i) +tmp
          epi(j)= epi(j) +tmp
          esrl = esrl +tmp +tmp
        else
          epi(i)= epi(i) +tmp
          esrl = esrl +tmp
        endif
!.....force
        ftmp = -acc *qi*qj*diji*diji /dielec
        aa(1:3,i)= aa(1:3,i) -dxdi(1:3)*ftmp
        aa(1:3,j)= aa(1:3,j) -dxdj(1:3)*ftmp
!.....stress
        do ixyz=1,3
          do jxyz=1,3
            strsl(jxyz,ixyz,i)= strsl(jxyz,ixyz,i) &
                 -0.5d0 *ftmp*rij(ixyz)*(-dxdi(jxyz))
            strsl(jxyz,ixyz,j)= strsl(jxyz,ixyz,j) &
                 -0.5d0 *ftmp*rij(ixyz)*(-dxdi(jxyz))
          enddo
        enddo
      enddo
    enddo

  end subroutine force_direct
!=======================================================================
  subroutine force_direct_cut(namax,natm,tag,ra,nnmax,aa,strsl,chg,h,hi &
       ,lspr,epi,esrl,iprint,lstrs,rc,l1st)
!
!  Direct Coulomb interaction with cutoff radius.
!  Coulomb potential is shifted by v(rc) and modified by the term (r-rc)*v'(rc).
!  Thus the potential form is:
!    v_cut(r) = v(r) -v(rc) -(r-rc)*v'(rc)
!  where v(r) is the direct Coulomb potential
!    v(r) = k*qi*qj*(1/r)
!
    implicit none
    integer,intent(in):: namax,natm,nnmax,iprint, &
         lspr(0:nnmax,namax)
    real(8),intent(in)::tag(namax),ra(3,namax),chg(namax), &
         h(3,3),hi(3,3),rc
    logical,intent(in):: lstrs,l1st
    real(8),intent(inout):: aa(3,namax),strsl(3,3,namax), &
         epi(namax),esrl

    integer:: i,j,jj,is,js,ixyz,jxyz
    real(8):: rc2,sqpi,qi,qj,dij,diji,tmp,ftmp,terfc&
         ,dvdrc,vrc
    real(8):: xi(3),xj(3),xij(3),rij(3),dxdi(3),dxdj(3)

!.....Compute direct sum
    rc2 = rc*rc
    esrl = 0d0
!$omp parallel
!$omp do private(i,xi,is,qi,jj,j,js,qj,xj,xij,rij,dij,diji,dxdi, &
!$omp     dxdj,vrc,dvdrc,tmp,ftmp,ixyz,jxyz) &
!$omp     reduction(+:esrl)
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is= int(tag(i))
      qi = chg(i)
      do jj=1,lspr(0,i)
        j = lspr(jj,i)
        js = int(tag(j))
        if( .not.interact(is,js) ) cycle
        qj = chg(j)
        xj(1:3) = ra(1:3,j)
        xij(1:3)= xj(1:3)-xi(1:3)
        rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij = rij(1)**2 +rij(2)**2 +rij(3)**2
        if( dij.gt.rc2 ) cycle
        dij = sqrt(dij)
        diji = 1d0/dij
        dxdi(1:3)= -rij(1:3)*diji
        dxdj(1:3)=  rij(1:3)*diji
        vrc = acc*qi*qj/rc
        dvdrc = -acc*qi*qj /rc**2
!.....potential
        tmp = 0.5d0 *(acc *qi*qj*diji -vrc -dvdrc*(dij-rc) )
        tmp = tmp /dielec
        epi(i)= epi(i) +tmp
        esrl= esrl +tmp
!.....force
        ftmp = -acc *qi*qj*diji*diji -dvdrc
        ftmp = ftmp /dielec
        aa(1:3,i)= aa(1:3,i) -dxdi(1:3)*ftmp
!.....stress
        do ixyz=1,3
          do jxyz=1,3
            strsl(jxyz,ixyz,i)= strsl(jxyz,ixyz,i) &
                 -0.5d0 *ftmp*rij(ixyz)*(-dxdi(jxyz))
          enddo
        enddo
      enddo
    enddo
!$omp end do
!$omp end parallel

  end subroutine force_direct_cut
!=======================================================================
  subroutine force_screened_cut(namax,natm,tag,ra,nnmax,aa,strsl &
       ,chg,h,hi,lspr,epi,esrl,iprint,lstrs,rc,l1st)
!
!  Screened Coulomb with cutoff that uses rho_scr.
!  smoothing using vrc and dVdrc where
!    V_smooth(r) = V(r) -V(rc) -(r-rc)*dVdrc
!
    implicit none
    include "mpif.h"
    include "./params_unit.h"
    integer,intent(in):: namax,natm,nnmax,iprint, &
         lspr(0:nnmax,namax)
    real(8),intent(in):: ra(3,namax),h(3,3),hi(3,3),rc,tag(namax)
    real(8),intent(inout):: chg(namax)
    real(8),intent(inout):: strsl(3,3,namax),aa(3,namax)&
         ,epi(namax),esrl
    logical,intent(in):: lstrs,l1st

    integer:: i,j,jj,k,kk,l,m,n,ierr,is,js,ixyz,jxyz,nconnect(4)
    integer:: itot,jtot
    real(8):: xi(3),xj(3),xij(3),rij(3),dij,diji,dedr &
         ,dxdi(3),dxdj(3),x,y,z,epotl,epott,at(3),tmp &
         ,qi,qj,radi,radj,rhoij,terfc,texp &
         ,vrc,dvdrc,terfcc,ri,ro,xs,dz,fc,dfc
    real(8),save:: sqpi,rc2

    if( l1st ) then
      sqpi = 1d0/sqrt(pi)
      rc2 = rc*rc
      do is=1,nspmax
        do js=is,nspmax
          if( .not.interact(is,js) ) cycle
          rhoij = rho_scr(is,js)
          terfcc = erfc(rc/rhoij)
          vrc = acc /rc *terfcc
          dvdrc = -acc /rc &
               *(terfcc/rc +2d0/rhoij *sqpi *exp(-(rc/rhoij)**2))
          vrcs(is,js) = vrc
          vrcs(js,is) = vrc
          dvdrcs(is,js) = dvdrc
          dvdrcs(js,is) = dvdrc
        enddo
      enddo
    endif

    esrl= 0d0
!.....Loop over resident atoms
!$omp parallel
!$omp do private(i,xi,is,qi,jj,j,js,qj,xj,xij,rij,dij,diji,dxdi, &
!$omp     dxdj,rhoij,terfc,vrc,dvdrc,texp,dedr,tmp,ixyz,jxyz) &
!$omp     reduction(+:esrl)
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is= int(tag(i))
      qi= chg(i)
      do jj=1,lspr(0,i)
        j=lspr(jj,i)
!!$        if(j.eq.0) exit
!!$        if(j.le.i) cycle
        js= int(tag(j))
        if( .not.interact(is,js) ) cycle
        qj= chg(j)
        xj(1:3)= ra(1:3,j)
        xij(1:3)= xj(1:3)-xi(1:3)
        rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij= rij(1)**2 +rij(2)**2 +rij(3)**2
        if( dij.gt.rc2 ) cycle
        dij = sqrt(dij)
        diji= 1d0/dij
        dxdi(1:3)= -rij(1:3)*diji
        dxdj(1:3)=  rij(1:3)*diji
        rhoij = rho_scr(is,js)
        terfc = erfc(dij/rhoij)
        vrc = qi*qj*vrcs(is,js)
        dvdrc = qi*qj*dvdrcs(is,js)
        texp = exp(-(dij/rhoij)**2)
        dedr= -acc *qi*qj*diji *(1d0*diji*terfc +2d0/rhoij *sqpi *texp) -dvdrc
        tmp= 0.5d0 *( acc *qi*qj*diji *terfc -vrc -dvdrc*(dij-rc) )
        tmp = tmp /dielec
!.....potential
!!$        if( j.le.natm ) then
!!$          epi(i)= epi(i) +tmp
!!$          epi(j)= epi(j) +tmp
!!$          esrl = esrl +tmp +tmp
!!$        else
!!$          epi(i)= epi(i) +tmp
!!$          esrl = esrl +tmp
!!$        endif
        epi(i)= epi(i) +tmp
        esrl = esrl +tmp
!.....force
        dedr = dedr /dielec
        aa(1:3,i)= aa(1:3,i) -dxdi(1:3)*dedr
!!$        aa(1:3,j)= aa(1:3,j) -dxdj(1:3)*dedr
!.....stress
        do ixyz=1,3
          do jxyz=1,3
            strsl(jxyz,ixyz,i)= strsl(jxyz,ixyz,i) &
                 -0.5d0 *dedr*rij(ixyz)*(-dxdi(jxyz))
!!$              strsl(jxyz,ixyz,j)= strsl(jxyz,ixyz,j) &
!!$                   -0.5d0 *dedr*rij(ixyz)*(-dxdi(jxyz))
          enddo
        enddo
      enddo
    enddo
!$omp end do
!$omp end parallel

  end subroutine force_screened_cut
!=======================================================================
  subroutine Ewald_short(namax,natm,tag,ra,nnmax,aa,strsl,chg,h,hi &
       ,lspr,epi,esrl,iprint,lstrs,rc)
    implicit none
    integer,intent(in):: namax,natm,nnmax,iprint, &
         lspr(0:nnmax,namax)
    real(8),intent(in)::tag(namax),ra(3,namax),chg(namax), &
         h(3,3),hi(3,3),rc
    logical,intent(in):: lstrs
    real(8),intent(inout):: aa(3,namax),strsl(3,3,namax), &
         epi(namax),esrl

    integer:: i,j,jj,is,js,ixyz,jxyz
    real(8):: rc2,sgmsq2,ss2i,sqpi,qi,qj,dij,diji,tmp,ftmp,terfc
    real(8):: xi(3),xj(3),xij(3),rij(3),dxdi(3),dxdj(3)
    real(8),external:: fcut1,dfcut1

!.....Compute direct sum
    rc2 = rc*rc
    sgmsq2 = sqrt(2d0)*sgm_ew
    ss2i = 1d0 /sgmsq2
    sqpi = 1d0 /sqrt(pi)
    esrl = 0d0
!$omp parallel
!$omp do private(i,xi,is,qi,jj,j,js,qj,xj,xij,rij,dij,diji,dxdi, &
!$omp            dxdj,terfc,tmp,ftmp,ixyz,jxyz) &
!$omp     reduction(+:esrl)
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is= int(tag(i))
      qi = chg(i)
      do jj=1,lspr(0,i)
        j = lspr(jj,i)
!!$        if( j.le.i ) cycle
        js = int(tag(j))
        qj = chg(j)
        xj(1:3) = ra(1:3,j)
        xij(1:3)= xj(1:3)-xi(1:3)
        rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij = rij(1)**2 +rij(2)**2 +rij(3)**2
        if( dij.ge.rc2 ) cycle
        dij = sqrt(dij)
        diji = 1d0/dij
        dxdi(1:3)= -rij(1:3)*diji
!!$        dxdj(1:3)=  rij(1:3)*diji
        terfc = erfc(dij*ss2i)
!.....potential
        tmp = 0.5d0 *acc *qi*qj*diji *terfc
        tmp = tmp /dielec
!!$        if( j.le.natm ) then
!!$          epi(i)= epi(i) +tmp
!!$          esrl = esrl +tmp +tmp
!!$!$omp atomic
!!$          epi(j)= epi(j) +tmp
!!$        else
!!$          epi(i)= epi(i) +tmp
!!$          esrl = esrl +tmp
!!$        endif
        epi(i)= epi(i) +tmp
        esrl= esrl +tmp
!.....force
        ftmp = -acc *qj*qi*diji *( diji *terfc &
             +2d0 *sqpi *ss2i *exp(-(dij*ss2i)**2) )
        ftmp = ftmp /dielec
        aa(1:3,i)= aa(1:3,i) -dxdi(1:3)*ftmp
!!$        do ixyz=1,3
!!$!$omp atomic
!!$          aa(ixyz,j)= aa(ixyz,j) +dxdi(ixyz)*ftmp
!!$        enddo
!.....stress
        do ixyz=1,3
          do jxyz=1,3
            strsl(jxyz,ixyz,i)= strsl(jxyz,ixyz,i) &
                 -0.5d0 *ftmp*rij(ixyz)*(-dxdi(jxyz))
!!$!$omp atomic
!!$            strsl(jxyz,ixyz,j)= strsl(jxyz,ixyz,j) &
!!$                 -0.5d0 *ftmp*rij(ixyz)*(-dxdi(jxyz))
          enddo
        enddo
      enddo
    enddo
!$omp end do
!$omp end parallel

  end subroutine Ewald_short
!=======================================================================
  subroutine Ewald_long(namax,natm,tag,ra,nnmax,aa,strsl,chg,h,hi,&
       lspr,sorg,epi,elrl,iprint,myid,mpi_md_world,lstrs,lcell_updated)
!
!  Long-range term of Ewald sum.
!  Not parallelized and not suitable for large system as it is too slow.
!
    implicit none
    include 'mpif.h'
    integer,intent(in):: namax,natm,nnmax,iprint, &
         myid,mpi_md_world,lspr(0:nnmax,namax)
    real(8),intent(in):: tag(namax),ra(3,namax),chg(namax),&
         h(3,3),hi(3,3),sorg(3)
    logical,intent(in):: lstrs,lcell_updated
    real(8),intent(inout):: aa(3,namax),epi(namax),&
         elrl,strsl(3,3,namax)

    integer:: i,is,ik,k1,k2,k3,ierr,ixyz,jxyz,itot
    real(8):: qi,xi(3),ri(3),bk1(3),bk2(3),bk3(3),bb(3),bdotr, &
         tmp,ftmp,cs,sn,texp,bb2,bk
    real(8):: bdk1,bdk2,bdk3,cs1,sn1,cs2,sn2,cs3,sn3,cs10,cs20,cs30, &
         cs1m,cs1mm,sn1m,sn1mm,cs2m,cs2mm,sn2m,sn2mm,cs3m,cs3mm,sn3m,sn3mm
    real(8):: emat(3,3)

!.....Compute reciprocal vectors
    if( lcell_updated ) call get_recip_vectors(h)
!.....Compute structure factor of the local processor
    call calc_qcos_qsin(namax,natm,tag,ra,chg,h,iprint &
         ,myid,mpi_md_world,sorg)

!.....Compute long-range contribution to potential energy
    elrl = 0d0
!.....Long-range contribution to forces and stress
    emat(1:3,1) = (/ 1d0, 0d0, 0d0 /)
    emat(1:3,2) = (/ 0d0, 1d0, 0d0 /)
    emat(1:3,3) = (/ 0d0, 0d0, 1d0 /)
!$omp parallel
!$omp do private(i,xi,is,itot,qi,ri,ik,bdk1,bdk2,bdk3,cs10,cs20,cs30, &
!$omp            cs1m,cs1mm,sn1m,sn1mm,cs2m,cs2mm,sn2m,sn2mm,cs3m,cs3mm, &
!$omp            sn3m,sn3mm,k1,k2,k3,bk1,bk2,bk3,cs1,sn1,cs2,sn2,cs3,sn3, &
!$omp            bb,cs,sn,tmp,ftmp,bk,ixyz,jxyz) &
!$omp     reduction(+:elrl)
    do i=1,natm
      xi(1:3)= ra(1:3,i) +sorg(1:3)
      is= int(tag(i))
      itot = itotOf(tag(i))
      qi = chg(i)
!!$      if( abs(qi).lt.qthd ) cycle
      ri(1:3) = h(1:3,1)*xi(1) +h(1:3,2)*xi(2) +h(1:3,3)*xi(3)
      ik = 0
      bdk1 = dot(b1,ri)
      bdk2 = dot(b2,ri)
      bdk3 = dot(b3,ri)
      cs10 = cos(bdk1)
      cs20 = cos(bdk2)
      cs30 = cos(bdk3)
!!$      if( itot.eq.19 .or. itot.eq.21 ) then
!!$        print '(a,3i4,7es11.3)','myid,i,itot,xi,qi,cs10,cs20,cs30=' &
!!$             ,myid,i,itot,xi(1:3),qi,cs10,cs20,cs30
!!$      endif
      cs1m = 0d0
      cs1mm= 0d0
      sn1m = 0d0
      sn1mm= 0d0
      cs2m = 0d0
      cs2mm= 0d0
      sn2m = 0d0
      sn2mm= 0d0
      cs3m = 0d0
      cs3mm= 0d0
      sn3m = 0d0
      sn3mm= 0d0
      do k1= -kmax1,kmax1
        bk1(1:3) = k1 *b1(1:3)
        if( k1.lt.-kmax1+2 ) then
          cs1 = cos(k1*bdk1)
          sn1 = sin(k1*bdk1)
          if( k1.eq.-kmax1   ) then
            cs1mm= cs1
            sn1mm= sn1
          else if( k1.eq.-kmax1+1 ) then
            cs1m = cs1
            sn1m = sn1
          endif
        else
          cs1 = 2d0*cs1m*cs10 -cs1mm
          sn1 = 2d0*sn1m*cs10 -sn1mm
        endif
        do k2= -kmax2,kmax2
          bk2(1:3) = k2 *b2(1:3)
          if( k2.lt.-kmax2+2 ) then
            cs2 = cos(k2*bdk2)
            sn2 = sin(k2*bdk2)
            if( k2.eq.-kmax2 ) then
              cs2mm = cs2
              sn2mm = sn2
            else if( k2.eq.-kmax2+1 ) then
              cs2m = cs2
              sn2m = sn2
            endif
          else
            cs2 = 2d0*cs2m*cs20 -cs2mm
            sn2 = 2d0*sn2m*cs20 -sn2mm
          endif
          do k3= -kmax3,kmax3
            bk3(1:3) = k3 *b3(1:3)
            if( k3.lt.-kmax3+2 ) then
              cs3 = cos(k3*bdk3)
              sn3 = sin(k3*bdk3)
              if( k3.eq.-kmax3 ) then
                cs3mm = cs3
                sn3mm = sn3
              else if( k3.eq.-kmax3+1 ) then
                cs3m = cs3
                sn3m = sn3
              endif
            else
              cs3 = 2d0*cs3m*cs30 -cs3mm
              sn3 = 2d0*sn3m*cs30 -sn3mm
            endif
            if( .not. lkuse(k3,k2,k1) ) goto 10
            ik= ik +1
            bb(1:3) = bk1(1:3) +bk2(1:3) +bk3(1:3)
!.....Since the computations of cos and sin are rather heavy,
!     compute cos and sin of each direction are computed iterative form
!     and compute cos(bb*ri)=cos(bk1*ri +bk2*ri +bk3*ri),
!     and compute it using cos(bk1*ri)'s and sin(bk1*ri)'s
            cs = cs1*cs2*cs3 -cs1*sn2*sn3 -sn1*cs2*sn3 -sn1*sn2*cs3
            sn = cs1*cs2*sn3 +cs1*sn2*cs3 +sn1*cs2*cs3 -sn1*sn2*sn3
!.....Potential energy per atom
!!$            tmp = 0.5d0 *acc /vol *qi *pflr(ik,is) &
!!$                 *( cs*qcos(ik) +sn*qsin(ik) )
            tmp = 0.5d0 *acc /vol *qi *pflr(ik,is) &
                 *( cs*qcos(ik) +sn*qsin(ik) )
            tmp = tmp /dielec
            epi(i) = epi(i) +tmp
            elrl = elrl +tmp
!.....Forces
!!$            aa(1:3,i)= aa(1:3,i) -acc/vol *qi*bb(1:3) *pflr(ik,is) &
!!$                 *0.5d0*( -sn*qcos(ik) +cs*qsin(ik) )
!!$            aa(1:3,i)= aa(1:3,i) -acc/vol *qi*bb(1:3) *pflr(ik,is) &
!!$                 *( -sn*qcos(ik) +cs*qsin(ik) )
            ftmp = -acc/vol *qi *pflr(ik,is) &
                 *( -sn*qcos(ik) +cs*qsin(ik) ) /dielec
            aa(1:3,i) = aa(1:3,i) +bb(1:3)*ftmp
!!$            if( itot.eq.19 .or. itot.eq.21 ) then
!!$              print '(a,4i4,7es11.3)','myid,i,itot,ik,aa,qcos,qsin,cs,sn=' &
!!$                   ,myid,i,itot,ik,aa(1:3,i) &
!!$                   ,qcos(ik),qsin(ik),cs,sn
!!$            endif
!.....Stress
            bk = norm(bb)
            do ixyz=1,3
              do jxyz=1,3
                strsl(ixyz,jxyz,i) = strsl(ixyz,jxyz,i) +tmp &
                     *( bb(ixyz)*bb(jxyz)/bk*(bk *sgm_ew**2 +2d0/bk) &
                     -emat(ixyz,jxyz))
              enddo
            enddo

10          if( k3.ge.-kmax3+2 ) then
              cs3mm= cs3m
              cs3m = cs3
              sn3mm= sn3m
              sn3m = sn3
            endif
          enddo  ! k3
          if( k2.ge.-kmax2+2 ) then
            cs2mm= cs2m
            cs2m = cs2
            sn2mm= sn2m
            sn2m = sn2
          endif
        enddo  ! k2
        if( k1.ge.-kmax1+2 ) then
          cs1mm= cs1m
          cs1m = cs1
          sn1mm= sn1m
          sn1m = sn1
        endif
      enddo  ! k1
    enddo  ! i
!$omp end do
!$omp end parallel

  end subroutine Ewald_long
!=======================================================================
  subroutine self_term(namax,natm,tag,chg,&
       epi,eselfl,iprint,lvc)
    implicit none
    integer,intent(in):: namax,natm,iprint
    logical,intent(in):: lvc
    real(8),intent(in):: tag(namax),chg(namax)
    real(8),intent(inout):: epi(namax),eselfl

    integer:: i,is
    real(8):: sgmi,qi,q2,tmp

!.....Compute self term.
    if( lvc ) then  ! variable charge
      eselfl = 0d0
      do i=1,natm
        is = int(tag(i))
        qi = chg(i)
        q2 = qi*qi
        tmp = (vc_e0(is) +vc_chi(is)*qi +0.5d0*vc_jii(is)*q2)
        eselfl = eselfl +tmp
        epi(i) = epi(i) +tmp
      enddo
    else if( trim(cterms).eq.'full' .or. &
         trim(cterms).eq.'long') then ! fixed charge
!....If charge per atom is fixed, it is constant, though.
      eselfl = 0d0
      do i=1,natm
        is = int(tag(i))
        q2 = chg(i)*chg(i)
        tmp = -q2 /sgm_ew *acc/sqrt(2d0*pi)
        eselfl = eselfl +tmp
        epi(i) = epi(i) +tmp
      enddo
    endif
    return
  end subroutine self_term
!=======================================================================
  subroutine qforce_short(namax,natm,tag,ra,nnmax,chg,h &
       ,lspr,iprint,rc,fq,esr)
!
!  Compute q-force of short-range term of Ewald sum.
!
    implicit none
    integer,intent(in):: namax,natm,nnmax,iprint,lspr(0:nnmax,namax)
    real(8),intent(in)::tag(namax),ra(3,namax),chg(namax),h(3,3),rc
    real(8),intent(inout):: fq(namax),esr

    integer:: i,j,jj,is,js,ixyz,jxyz
    real(8):: rc2,sgmsq2,ss2i,sqpi,qi,qj,dij,diji,tmp,ftmp,terfc &
         ,sgmi,sgmj,gmmij
    real(8):: xi(3),xj(3),xij(3),rij(3),dxdi(3),dxdj(3)
    real(8),external:: fcut1

!.....Compute direct sum
    rc2 = rc*rc
    sgmsq2 = sqrt(2d0)*sgm_ew
    ss2i = 1d0 /sgmsq2
    sqpi = 1d0 /sqrt(pi)
    esr = 0d0
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is= int(tag(i))
      qi = chg(i)
!!$      if( abs(qi).lt.qthd ) cycle
!!$      sgmi = sgm(is)
      do jj=1,lspr(0,i)
        j = lspr(jj,i)
        if( j.eq.0 ) exit
        if( j.le.i ) cycle
        js = int(tag(j))
        qj = chg(j)
!!$        if( abs(qj).lt.qthd ) cycle
        xj(1:3) = ra(1:3,j)
        xij(1:3)= xj(1:3)-xi(1:3)
        rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij = rij(1)**2 +rij(2)**2 +rij(3)**2
        if( dij.gt.rc2 ) cycle
!!$        sgmj = sgm(js)
!!$        gmmij = 1d0 /sqrt(2d0*(sgmi**2+sgmj**2))
        dij = sqrt(dij)
        diji = 1d0/dij
        dxdi(1:3)= -rij(1:3)*diji
        dxdj(1:3)=  rij(1:3)*diji
        terfc = erfc(dij*ss2i)
!!$        terfc = erfc(gmmij*dij)
!.....potential
        tmp = acc *diji *terfc *fcut1(dij,0d0,rc) /dielec
        if( j.le.natm ) then
          esr = esr +tmp*qi*qj
        else
          esr = esr +0.5d0*tmp*qi*qj
        endif
!.....Force on charge
        fq(i) = fq(i) -tmp*qj
        fq(j) = fq(j) -tmp*qi
      enddo
    enddo

    return
  end subroutine qforce_short
!=======================================================================
  subroutine qforce_direct_cut(namax,natm,tag,ra,nnmax,chg,h &
       ,lspr,iprint,rc,fq,esr,l1st)
!
!  Compute q-force of direct_cut Coulomb potential.
!
    implicit none
    integer,intent(in):: namax,natm,nnmax,iprint, &
         lspr(0:nnmax,namax)
    real(8),intent(in)::tag(namax),ra(3,namax),chg(namax), &
         h(3,3),rc
    real(8),intent(inout):: fq(namax),esr
    logical,intent(in):: l1st 

    integer:: i,j,jj,is,js,ixyz,jxyz
    real(8):: sgmsq2,ss2i,qi,qj,dij,dij2,diji,tmp,ftmp, &
         vrc,dvdrc
    real(8):: xi(3),xj(3),xij(3),rij(3),dxdi(3),dxdj(3)
!!$    real(8),external:: fcut1

    real(8),save:: rc2,sqpi

!.....Compute direct sum
    rc2 = rc*rc
    esr = 0d0
!$omp parallel
!$omp do private(i,xi,xj,xij,rij,is,qi,jj,j,js,qj,dij,diji,vrc,dvdrc,tmp) &
!$omp    reduction(+:esr)
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is= int(tag(i))
      qi = chg(i)
      do jj=1,lspr(0,i)
        j = lspr(jj,i)
        js = int(tag(j))
        if( .not.interact(is,js) ) cycle
        qj = chg(j)
        xj(1:3) = ra(1:3,j)
        xij(1:3)= xj(1:3)-xi(1:3)
        rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij = rij(1)**2 +rij(2)**2 +rij(3)**2
        if( dij.gt.rc2 ) cycle
        dij = sqrt(dij)
        diji = 1d0/dij
        vrc = acc/rc
        dvdrc = -acc/rc2
!.....potential
        tmp = acc*diji -vrc -dvdrc*(dij-rc)
        tmp = tmp /dielec
        esr = esr +0.5d0*tmp*qi*qj
!.....Force on charge
        fq(i) = fq(i) -tmp*qj
      enddo
    enddo
!$omp end do
!$omp end parallel

    return
  end subroutine qforce_direct_cut
!=======================================================================
  subroutine qforce_screened_cut(namax,natm,tag,ra,nnmax,chg,h &
       ,lspr,iprint,rc,fq,esr,l1st)
!
!  Compute q-force of screened_cut Coulomb potential.
!
    implicit none
    integer,intent(in):: namax,natm,nnmax,iprint, &
         lspr(0:nnmax,namax)
    real(8),intent(in)::tag(namax),ra(3,namax),chg(namax), &
         h(3,3),rc
    real(8),intent(inout):: fq(namax),esr
    logical,intent(in):: l1st 

    integer:: i,j,jj,is,js,ixyz,jxyz
    real(8):: sgmsq2,ss2i,qi,qj,dij,dij2,diji,tmp,ftmp,terfc &
         ,sgmi,sgmj,gmmij,rhoij,terfcc,vrc,dvdrc
    real(8):: xi(3),xj(3),xij(3),rij(3),dxdi(3),dxdj(3)
!!$    real(8),external:: fcut1

    real(8),save:: rc2,sqpi

    if( l1st ) then
      rc2 = rc*rc
      sqpi = 1d0 /sqrt(pi)
      do is=1,nspmax
        do js=is,nspmax
          if( .not.interact(is,js) ) cycle
          rhoij = rho_scr(is,js)
          terfcc = erfc(rc/rhoij)
          vrc = acc /rc *terfcc
          dvdrc = -acc /rc &
               *(terfcc/rc +2d0/rhoij *sqpi *exp(-(rc/rhoij)**2))
          vrcs(is,js) = vrc
          vrcs(js,is) = vrc
          dvdrcs(is,js) = dvdrc
          dvdrcs(js,is) = dvdrc
        enddo
      enddo
    endif

!.....Compute direct sum
    sgmsq2 = sqrt(2d0)*sgm_ew
    ss2i = 1d0 /sgmsq2
    esr = 0d0
!$omp parallel
!$omp do private(i,xi,is,qi,jj,j,js,qj,xj,xij,rij,dij,diji,rhoij,terfc,vrc,dvdrc,tmp) &
!$omp    reduction(+:esr)
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is= int(tag(i))
      qi = chg(i)
      do jj=1,lspr(0,i)
        j = lspr(jj,i)
!!$        if( j.eq.0 ) exit
!!$        if( j.le.i ) cycle
        js = int(tag(j))
        if( .not.interact(is,js) ) cycle
        qj = chg(j)
        xj(1:3) = ra(1:3,j)
        xij(1:3)= xj(1:3)-xi(1:3)
        rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij = rij(1)**2 +rij(2)**2 +rij(3)**2
        if( dij.gt.rc2 ) cycle
        dij = sqrt(dij)
        diji = 1d0/dij
        rhoij = rho_scr(is,js)
        terfc = erfc(dij/rhoij)
        vrc = vrcs(is,js)
        dvdrc = dvdrcs(is,js)
!.....potential
        tmp = acc*diji*terfc -vrc -dvdrc*(dij-rc)
        tmp = tmp /dielec
!!$        if( j.le.natm ) then
!!$          esr = esr +tmp*qi*qj
!!$        else
!!$          esr = esr +0.5d0*tmp*qi*qj
!!$        endif
        esr = esr +0.5d0*tmp*qi*qj        
!.....Force on charge
        fq(i) = fq(i) -tmp*qj
!!$        fq(j) = fq(j) -tmp*qi
      enddo
    enddo
!$omp end do
!$omp end parallel

    return
  end subroutine qforce_screened_cut
!=======================================================================
  subroutine qforce_long(namax,natm,tag,ra,chg,h, &
       sorg,mpi_md_world,myid,iprint,fq,elr)
!
!  Derivative of Ewald long-range term w.r.t. charges
!
    implicit none
    include 'mpif.h'
    integer,intent(in):: namax,natm,mpi_md_world,myid,iprint
    real(8),intent(in):: tag(namax),ra(3,namax),chg(namax),h(3,3),sorg(3)
    real(8),intent(inout):: fq(namax),elr

    integer:: i,ik,k1,k2,k3,is,ierr
    real(8):: prefac,bk1(3),bk2(3),bk3(3),bb(3),bb2,xi(3),ri(3), &
         sgmi,sgmi2,qi,bdotr,texp,cs,sn,elrl,eselfl,q2,epotl,tmp

!.....Compute reciprocal vectors
    call get_recip_vectors(h)
    prefac = 1d0 /(2d0*vol*eps0)
!.....Compute structure factor
    call calc_qcos_qsin(namax,natm,tag,ra,chg,h,iprint &
         ,myid,mpi_md_world,sorg)

    ik = 0
    elr = 0d0
    do k1= -kmax1,kmax1
      bk1(1:3) = k1 *b1(1:3)
      do k2= -kmax2,kmax2
        bk2(1:3) = k2 *b2(1:3)
        do k3= -kmax3,kmax3
          if( .not. lkuse(k3,k2,k1) ) cycle
          ik= ik +1
          bk3(1:3) = k3 *b3(1:3)
          bb(1:3) = bk1(1:3) +bk2(1:3) +bk3(1:3)
          bb2 = dot(bb,bb)
          do i=1,natm
            xi(1:3)= ra(1:3,i)
            is= int(tag(i))
            sgmi = sgm_ew
            sgmi2= sgmi*sgmi
            qi = chg(i)
            ri(1:3) = h(1:3,1)*xi(1) +h(1:3,2)*xi(2) +h(1:3,3)*xi(3)
            bdotr = dot(bb,ri)
            texp = exp(-bb2*sgmi2/2)
            cs = cos(bdotr)
            sn = sin(bdotr)
!.....Potential energy per atom
            tmp = prefac/bb2*qi*texp &
                 *(cs*qcos(ik) +sn*qsin(ik))
            tmp = tmp /dielec
            elr = elr +tmp
!.....Force on charge
            fq(i)= fq(i) -2d0 *prefac/bb2 *texp &
                 *(cs*qcos(ik) +sn*qsin(ik)) /dielec
!!$            if( i.eq.1 .and. (ik.gt.5000.and.ik.le.6000) ) then
!!$              print '(a,i6,3i4,10f10.4)','ik,qi,1/bb2,texp,cs,sn,qcos,qsin,fqikkk,fq(i)='&
!!$                   ,ik,k1,k2,k3,qi,1.d0/bb2,texp,cs,sn &
!!$                   ,qcos(ik),qsin(ik) &
!!$                   , -2d0 *prefac/bb2 *texp &
!!$                   *(cs*qcos(ik) +sn*qsin(ik)),fq(i)
!!$            endif
          enddo
        enddo
      enddo
    enddo

!!$    epotl = eselfl +elrl
!!$    call mpi_allreduce(epotl,epot,1,mpi_real8, &
!!$         mpi_sum,mpi_md_world,ierr)
    return
  end subroutine qforce_long
!=======================================================================
  subroutine qforce_self(namax,natm,tag,chg,fq,eself)
!
!  Derivative of self term w.r.t. charges
!
    integer,intent(in):: namax,natm
    real(8),intent(in):: tag(namax),chg(namax)
    real(8),intent(inout):: fq(namax),eself

    integer:: i,is
    real(8):: qi,q2,sgmi,tmp

    eself = 0d0
!$omp parallel
!$omp do private(i,is,qi,q2,sgmi,tmp) &
!$omp    reduction(+:eself)
    do i=1,natm
      is = int(tag(i))
      qi = chg(i)
      q2 = qi*qi
      sgmi = sgm_ew
      tmp = vc_e0(is) +vc_chi(is)*qi +0.5d0*vc_jii(is)*q2
      eself = eself +tmp
      fq(i) = fq(i) -(vc_chi(is) +vc_jii(is)*qi)
    enddo
!$omp end do
!$omp end parallel

  end subroutine qforce_self
!=======================================================================
  subroutine get_qforce(chg,fq,epot,l1st)
!
!  Wrapper routine for calculating forces on charges.
!
    use pmdvars,only: namax,natm,nnmax,tag,ra,h,lspr, &
         rc,sorg,myid_md,mpi_md_world,iprint
    include "mpif.h"
    real(8),intent(in):: chg(namax)
    logical,intent(in):: l1st
    real(8),intent(out):: fq(namax),epot

    integer:: ierr
    real(8):: eclongl,eselfl,ecshortl
    real(8):: eclong,eself,ecshort

    fq(:) = 0d0
    call qforce_self(namax,natm,tag,chg,fq,eselfl)
    call mpi_allreduce(eselfl,eself,1,mpi_real8,mpi_sum,mpi_md_world,ierr)

    if( trim(cterms).eq.'long' ) then
      call qforce_long(namax,natm,tag,ra,chg,h,sorg,mpi_md_world, &
           myid_md,iprint,fq,eclongl)
      call mpi_allreduce(eclongl,eclong,1,mpi_real8,mpi_sum,mpi_md_world,ierr)
    else if( trim(cterms).eq.'full' ) then
      call qforce_short(namax,natm,tag,ra,nnmax,chg,h,lspr,iprint &
           ,rc,fq,ecshortl)
      call qforce_long(namax,natm,tag,ra,chg,h,sorg,mpi_md_world, &
           myid_md,iprint,fq,eclongl)
      call mpi_allreduce(ecshortl,ecshort,1,mpi_real8,mpi_sum,mpi_md_world,ierr)
      call mpi_allreduce(eclongl,eclong,1,mpi_real8,mpi_sum,mpi_md_world,ierr)
    else if( trim(cterms).eq.'direct_cut'  ) then
      call qforce_direct_cut(namax,natm,tag,ra,nnmax,chg,h, &
           lspr,iprint,rc,fq,ecshortl,l1st)
      call mpi_allreduce(ecshortl,ecshort,1,mpi_real8,mpi_sum,mpi_md_world,ierr)
    else if( trim(cterms).eq.'short' .or. trim(cterms).eq.'screened' ) then
      call qforce_short(namax,natm,tag,ra,nnmax,chg,h,lspr,iprint &
           ,rc,fq,ecshortl)
      call mpi_allreduce(ecshortl,ecshort,1,mpi_real8,mpi_sum,mpi_md_world,ierr)
    else if( trim(cterms).eq.'screened_cut'  ) then
      call qforce_screened_cut(namax,natm,tag,ra,nnmax,chg,h, &
           lspr,iprint,rc,fq,ecshortl,l1st)
      call mpi_allreduce(ecshortl,ecshort,1,mpi_real8,mpi_sum,mpi_md_world,ierr)
    endif
    epot = eself +ecshort + eclong
!!$    print '(a,3f12.3)',' get_qforce: eself,eshort,elong=',eself,ecshort,eclong
    return
  end subroutine get_qforce
!=======================================================================
  subroutine set_charge_BVS(natm,nb,tag,chg,myid,mpi_md_world,iprint &
       ,specorder)
!
! Reset actual charges of atoms using effective charge information
! in order to keep charge neutrality in the system.
!
! This would be called only once at the beginning.
!
    include "mpif.h"
    integer,intent(in):: natm,nb,myid,mpi_md_world,iprint
    real(8),intent(in):: tag(natm+nb)
    real(8),intent(out):: chg(natm+nb)
    character(len=3),intent(in):: specorder(nspmax)

    integer,allocatable:: nbvsl(:),nbvs(:)
    integer:: i,is,ierr
    real(8):: sum_anion,sum_cation
    character(len=3):: csp

    allocate(nbvsl(nspmax),nbvs(nspmax))
    nbvsl(1:nspmax) = 0
    nbvs(1:nspmax) = 0
    do i=1,natm
      is = int(tag(i))
      nbvsl(is) = nbvsl(is) +1
    enddo

    call mpi_allreduce(nbvsl,nbvs,nspmax,mpi_integer &
         ,mpi_sum,mpi_md_world,ierr)

    sum_anion = 0d0
    sum_cation = 0d0
    do is=1,nspmax
      if( vid_bvs(is).lt.0d0 ) then  !anion
        sum_anion = sum_anion +vid_bvs(is) *nbvs(is) /sqrt(dble(npq_bvs(is)))
      else if( vid_bvs(is).gt.0d0 ) then !cation
        sum_cation = sum_cation +vid_bvs(is)*nbvs(is)/sqrt(dble(npq_bvs(is)))
      endif
    enddo
    sum_anion = abs(sum_anion)

!.....Compute valence charges of species
    do is=1,nspmax
      if( vid_bvs(is).lt.0d0 ) then  ! anion
        schg(is) = vid_bvs(is)/sqrt(dble(npq_bvs(is))) &
             *sqrt(sum_cation/sum_anion)
      else if( vid_bvs(is).gt.0d0 ) then ! cation
        schg(is)= vid_bvs(is)/sqrt(dble(npq_bvs(is))) &
             *sqrt(sum_anion/sum_cation)
      endif
    enddo

!!$!.....Overwrite charges for debugging...
!!$    schg(1:4) = (/ 1.27443, 0.78466, 1.10968, 3.20338/)

    if( myid.eq.0 .and. iprint.ge.ipl_basic ) then
      print *,''
      print *,'Charges fixed from ideal valences and composition:'
      do is=1,nspmax
        if( vid_bvs(is).eq.0d0 ) cycle
        csp = specorder(is)
        write(6,'(a,a4,2f7.3)') '   csp, V_ideal, V_actual = '&
             ,trim(csp),vid_bvs(is),schg(is)
      enddo
    endif

    do i=1,natm+nb
      is= int(tag(i))
      chg(i) = schg(is)
!!$!.....Negative charge for anion and positive for cation
!!$      if( is .eq. 1 ) then
!!$        chg(i)= -schg(is)
!!$      else
!!$        chg(i)= schg(is)
!!$      endif
    enddo

    deallocate(nbvsl,nbvs)

  end subroutine set_charge_BVS
!=======================================================================
  subroutine get_recip_vectors(h)
    implicit none
    real(8),intent(in):: h(3,3)

    real(8):: a1(3),a2(3),a3(3),a23(3),a12(3),a31(3),pi2

    a1(1:3) = h(1:3,1)
    a2(1:3) = h(1:3,2)
    a3(1:3) = h(1:3,3)
    a23 = cross(a2,a3)
    a31 = cross(a3,a1)
    a12 = cross(a1,a2)
    vol = abs(dot(a1,a23))
    pi2 = 2d0 *pi
    b1(1:3) = pi2 /vol *a23(1:3)
    b2(1:3) = pi2 /vol *a31(1:3)
    b3(1:3) = pi2 /vol *a12(1:3)
    return
  end subroutine get_recip_vectors
!=======================================================================
  subroutine calc_qcos_qsin(namax,natm,tag,ra,chg,h,iprint&
       ,myid,mpi_md_world,sorg)
!
!  Compute qcos and qsin needed for calculation of Ewald long-range term.
!
    implicit none
    include 'mpif.h'
    integer,intent(in):: namax,natm,iprint,myid,mpi_md_world
    real(8),intent(in):: tag(namax),ra(3,namax),chg(namax) &
         ,h(3,3),sorg(3)

    integer:: ik,k1,k2,k3,is,i,ierr
    real(8):: bk1(3),bk2(3),bk3(3),bb(3),xi(3),qi&
         ,ri(3),bdotr,cs,sn
    real(8):: bdk1,bdk2,bdk3,cs1,sn1,cs2,sn2,cs3,sn3,cs10,cs20,cs30, &
         cs1m,cs1mm,sn1m,sn1mm,cs2m,cs2mm,sn2m,sn2mm,cs3m,cs3mm,sn3m,sn3mm
    real(8):: qcmax,qsmax,qclmax,qslmax

!.....Compute structure factor of the local processor
    qcosl(1:nk) = 0d0
    qsinl(1:nk) = 0d0
!!$    ik = 0
!!$    do k1= -kmax1,kmax1
!!$      bk1(1:3) = k1 *b1(1:3)
!!$      do k2= -kmax2,kmax2
!!$        bk2(1:3) = k2 *b2(1:3)
!!$        do k3= -kmax3,kmax3
!!$          if( .not. lkuse(k3,k2,k1) ) cycle
!!$          ik= ik +1
!!$          bk3(1:3) = k3 *b3(1:3)
!!$          bb(1:3) = bk1(1:3) +bk2(1:3) +bk3(1:3)
!!$          do i=1,natm
!!$            xi(1:3)= ra(1:3,i)
!!$            qi = chg(i)
!!$            ri(1:3) = h(1:3,1)*xi(1) +h(1:3,2)*xi(2) +h(1:3,3)*xi(3)
!!$            bdotr = dot(bb,ri)
!!$            qcosl(ik) = qcosl(ik) +qi*cos(bdotr)
!!$            qsinl(ik) = qsinl(ik) +qi*sin(bdotr)
!!$          enddo
!!$        enddo
!!$      enddo
!!$    enddo  ! ia
!!$    print *,'myid,b1=',myid,b1(1:3)
!!$    print *,'myid,b2=',myid,b2(1:3)
!!$    print *,'myid,b3=',myid,b3(1:3)
!!$    print *,'myid,sorg=',myid,sorg(1:3)
    do i=1,natm
      xi(1:3)= ra(1:3,i) +sorg(1:3)
      qi = chg(i)
      ri(1:3) = h(1:3,1)*xi(1) +h(1:3,2)*xi(2) +h(1:3,3)*xi(3)
      ik = 0
      bdk1 = dot(b1,ri)
      bdk2 = dot(b2,ri)
      bdk3 = dot(b3,ri)
      cs10 = cos(bdk1)
      cs20 = cos(bdk2)
      cs30 = cos(bdk3)
      cs1m = 0d0
      cs1mm= 0d0
      sn1m = 0d0
      sn1mm= 0d0
      cs2m = 0d0
      cs2mm= 0d0
      sn2m = 0d0
      sn2mm= 0d0
      cs3m = 0d0
      cs3mm= 0d0
      sn3m = 0d0
      sn3mm= 0d0
      do k1= -kmax1,kmax1
!!$        bk1(1:3) = k1 *b1(1:3)
        if( k1.lt.-kmax1+2 ) then
          cs1 = cos(k1*bdk1)
          sn1 = sin(k1*bdk1)
          if( k1.eq.-kmax1   ) then
            cs1mm= cs1
            sn1mm= sn1
          else if( k1.eq.-kmax1+1 ) then
            cs1m = cs1
            sn1m = sn1
          endif
        else
          cs1 = 2d0*cs1m*cs10 -cs1mm
          sn1 = 2d0*sn1m*cs10 -sn1mm
        endif
        do k2= -kmax2,kmax2
!!$          bk2(1:3) = k2 *b2(1:3)
          if( k2.lt.-kmax2+2 ) then
            cs2 = cos(k2*bdk2)
            sn2 = sin(k2*bdk2)
            if( k2.eq.-kmax2 ) then
              cs2mm = cs2
              sn2mm = sn2
            else if( k2.eq.-kmax2+1 ) then
              cs2m = cs2
              sn2m = sn2
            endif
          else
            cs2 = 2d0*cs2m*cs20 -cs2mm
            sn2 = 2d0*sn2m*cs20 -sn2mm
          endif
          do k3= -kmax3,kmax3
            if( k3.lt.-kmax3+2 ) then
              cs3 = cos(k3*bdk3)
              sn3 = sin(k3*bdk3)
              if( k3.eq.-kmax3 ) then
                cs3mm = cs3
                sn3mm = sn3
              else if( k3.eq.-kmax3+1 ) then
                cs3m = cs3
                sn3m = sn3
              endif
            else
              cs3 = 2d0*cs3m*cs30 -cs3mm
              sn3 = 2d0*sn3m*cs30 -sn3mm
            endif
            if( .not. lkuse(k3,k2,k1) ) goto 10
            ik= ik +1
!!$            bk3(1:3) = k3 *b3(1:3)
!!$            bb(1:3) = bk1(1:3) +bk2(1:3) +bk3(1:3)
!!$            bdotr = dot(bb,ri)
            cs = cs1*cs2*cs3 -cs1*sn2*sn3 -sn1*cs2*sn3 -sn1*sn2*cs3
            sn = cs1*cs2*sn3 +cs1*sn2*cs3 +sn1*cs2*cs3 -sn1*sn2*sn3
            qcosl(ik) = qcosl(ik) +qi*cs
            qsinl(ik) = qsinl(ik) +qi*sn

10          if( k3.ge.-kmax3+2 ) then
              cs3mm= cs3m
              cs3m = cs3
              sn3mm= sn3m
              sn3m = sn3
            endif
          enddo ! k3
          if( k2.ge.-kmax2+2 ) then
            cs2mm= cs2m
            cs2m = cs2
            sn2mm= sn2m
            sn2m = sn2
          endif
        enddo ! k2
        if( k1.ge.-kmax1+2 ) then
          cs1mm= cs1m
          cs1m = cs1
          sn1mm= sn1m
          sn1m = sn1
        endif
      enddo ! k1
    enddo  ! ia
!.....Allreduce qcos and qsin, which could be inefficient and time consuming
    qcos(1:nk) = 0d0
    qsin(1:nk) = 0d0
    call mpi_allreduce(qcosl,qcos,nk,mpi_real8 &
         ,mpi_sum,mpi_md_world,ierr)
    call mpi_allreduce(qsinl,qsin,nk,mpi_real8 &
         ,mpi_sum,mpi_md_world,ierr)

!!$    qcmax = 0d0
!!$    qsmax = 0d0
!!$    qclmax = 0d0
!!$    qslmax = 0d0
!!$    do ik=1,nk
!!$      qcmax = max(qcos(ik),qcmax)
!!$      qsmax = max(qsin(ik),qsmax)
!!$      qclmax = max(qcosl(ik),qclmax)
!!$      qslmax = max(qsinl(ik),qslmax)
!!$    enddo
!!$    print '(a,i2,10es11.3)','myid,qcmax,qclmax,qsmax,qslmax=',myid &
!!$         ,qcmax,qclmax,qsmax,qslmax
    return
  end subroutine calc_qcos_qsin
!=======================================================================
  subroutine setup_kspace()
!
!  Estimate kmax# so that exp(-sgm^2*k^2/2)/k^2 outside kmax# is small enough.
!  And set flags to k-points to be used.
!
    real(8),allocatable:: bbs(:,:,:)
    integer:: k1,k2,k3
    real(8):: bk1(3),bk2(3),bk3(3),bk(3),bb,bb2,bbmax

    allocate(bbs(-kmaxini:kmaxini,-kmaxini:kmaxini,-kmaxini:kmaxini))

    bbmax = 0d0
    do k1= -kmaxini,kmaxini
      bk1(1:3) = k1*b1(1:3)
      do k2= -kmaxini,kmaxini
        bk2(1:3) = k2*b2(1:3)
        do k3= -kmaxini,kmaxini
          if( k1.eq.0 .and. k2.eq.0 .and. k3.eq.0 ) cycle
          bk3(1:3) = k3*b3(1:3)
          bk(1:3) = bk1(1:3) +bk2(1:3) +bk3(1:3)
          bb2 = norm(bk)
          bbs(k3,k2,k1) = bb2
!!$          bb2= bb2*bb2
!!$          bbs(k3,k2,k1) = exp(-0.5d0 *sgm_ew**2 *bb2)/bb2
!!$          bbmax = max(bbmax,bbs(k3,k2,k1))
        enddo
      enddo
    enddo

    kmax1= 1
    kmax2= 1
    kmax3= 1
    do k1= 1,kmaxini
      do k2= 1,kmaxini
        do k3= 1,kmaxini
          if( bbs(k3,k2,k1).le.bkmax ) then
            kmax1= max(kmax1,k1)
            kmax2= max(kmax2,k2)
            kmax3= max(kmax3,k3)
          endif
        enddo
      enddo
    enddo

    if( allocated(lkuse) ) then
      call accum_mem('force_Coulomb',-4*size(lkuse))
      deallocate(lkuse)
    endif
    allocate(lkuse(-kmax3:kmax3,-kmax2:kmax2,-kmax1:kmax1))
    call accum_mem('force_Coulomb',4*size(lkuse))

    lkuse(:,:,:) = .false.
    nk = 0
    do k1= -kmax1,kmax1
      do k2= -kmax2,kmax2
        do k3= -kmax3,kmax3
          if( k1.eq.0 .and. k2.eq.0 .and. k3.eq.0 ) cycle
!!$          print *,'k1,k2,k3,bbs(k3,k2,k1),bkmax =',k1,k2,k3,bbs(k3,k2,k1),bkmax
!!$          if( bbs(k3,k2,k1).gt.bbmax*threshold_kmax ) then
          if( bbs(k3,k2,k1).le.bkmax ) then
            lkuse(k3,k2,k1) = .true.
            nk = nk +1
          endif
        enddo
      enddo
    enddo

    deallocate(bbs)
  end subroutine setup_kspace
!=======================================================================
  subroutine chgopt_matrix(namax,natm,h,ra,tag,chg,nnmax,lspr, &
       rc,sorg,myid,mpi_world,iprint,l1st)
!
!  Charge optimization/equilibration by matrix inversion.
!  THIS CODE IS NOT COMPLETE AND NOT YET AVAILABLE...
!
    integer,intent(in):: namax,natm,nnmax,lspr(0:nnmax,namax)
    integer,intent(in):: myid,mpi_world,iprint
    real(8),intent(in):: h(3,3),ra(3,natm),tag(natm),rc,sorg(3)
    logical,intent(in):: l1st
    real(8),intent(inout):: chg(natm)

    integer:: ierr,i,j,is
    real(8):: epot
    real(8),allocatable,save:: amat(:,:),qvec(:),xvec(:),amati(:,:),fq(:)

    if( l1st ) then
      if( allocated(amat) ) then
        call accum_mem('force_Coulomb',-8*(2*size(amat)+size(qvec) &
             +size(xvec)+size(fq)))
        deallocate(amat,amati,qvec,xvec,fq)
      endif
      allocate(amat(natm+1,natm+1),amati(natm+1,natm+1),qvec(natm+1), &
           xvec(natm+1),fq(namax))
      call accum_mem('force_Coulomb',8*(2*size(amat)+size(qvec)+size(xvec)+size(fq)))
    endif

    amat(1:natm+1,1:natm+1) = 0d0
    qvec(natm+1) = 0d0

!.....Make A-matrix and X-vector for linear equation, A*Q=X
    do i=1,natm
      is = int(tag(i))
!.....Since qforces contain chi(i) that should not be in A-matrix,
!     extract chi from fq
      fq(i) = fq(i) +vc_chi(is)
!.....Instead, put them into X-vector
      xvec(i) = -vc_chi(is)
!.....Q-vector
      qvec(i) = chg(i)
    enddo
!.....MAKING A-MATRIX IN PARALLEL IMPLEMENTATION IS NOT SO EASY.
!.....THIS CODE DOES NOT GO WELL...
!!$    do i=1,natm
!!$      do j=1,natm
!!$        if( chg(i)
!!$      enddo
!!$    enddo
    amat(1:natm,natm+1) = 1d0
    xvec(natm+1) = 0d0  ! charge neutrality, qtot = 0.0

    if( iprint.ge.ipl_debug .and. myid.eq.0 .and. l1st ) then
      print *,'amat:'
      do i=1,natm+1
        write(6,'(10es10.2)') (amat(i,j),j=1,natm+1)
      enddo
      print *,'qvec:'
      write(6,'(10f8.4)') (qvec(j),j=1,natm+1)
      print *,'xvec:'
      write(6,'(10es10.2)') (xvec(j),j=1,natm+1)
    endif

    if( trim(chgopt_method).eq.'cg' ) then
!.....Perform CG optimization to equilibrate the system and get charges
      call cg(natm+1,natm+2,1d-5,amat,xvec,qvec,ierr)
    else if( trim(chgopt_method).eq.'matinv' ) then
!.....Get optimal qvec by matrix inversion
!!$      call choldc_inv(natm+1,amat,amati)
      call ludc_inv(natm+1,amat,amati)
      qvec(:) = 0d0
      do i=1,natm+1
        do j=1,natm+1
          qvec(i) = amati(i,j)*xvec(j)
        enddo
      enddo
    endif

    if( iprint.ge.ipl_debug .and. myid.eq.0 .and. l1st ) then
      print *,'qvec after cg/matinv:'
      write(6,'(10f10.6)') (qvec(j),j=1,natm)
      print *,'A*Q:'
      xvec = matmul(amat,qvec)
      write(6,'(10es10.2)') (xvec(j),j=1,natm+1)
    endif

!.....Restore charge to chg
    do i=1,natm
      chg(i) = qvec(i)
    enddo

    return
  end subroutine chgopt_matrix
!=======================================================================
  subroutine cg(ndim,nstp,eps,amat,b,x,ierr)
    implicit none
    integer,intent(in):: ndim,nstp
    integer,intent(out):: ierr
    real(8),intent(in):: amat(ndim,ndim),b(ndim),eps
    real(8),intent(inout):: x(ndim)

    integer:: istp
    real(8):: xp(ndim),bnrm,rr,rrp,beta,ap(ndim), &
         alpha,r(ndim),xd(ndim),p(ndim),xdnrm

    xp(1:ndim) = x(1:ndim)
    r = matmul(amat,xp)
    r(1:ndim) = b(1:ndim) -r(1:ndim)
    p(1:ndim) = r(1:ndim)
    bnrm = sqrt(dot_product(b,b))
    rr = dot_product(r,r)

    do istp=1,nstp
      ap = matmul(amat,p)
      alpha = rr/dot_product(p,ap)
      x(1:ndim) = xp(1:ndim) +alpha*p(1:ndim)
      r(1:ndim) = r(1:ndim) -alpha*ap(1:ndim)
!.....Check convergence
      xd(1:ndim) = x(1:ndim) -xp(1:ndim)
      xdnrm = sqrt(dot_product(xd,xd))
!      write(6,'(i5,8es10.2)') istp,xdnrm,x(1:5),xdnrm/bnrm
      if( xdnrm/bnrm .lt. eps ) then
        ierr = istp
        return
      endif
      rrp = rr
      rr = dot_product(r,r)
      beta = rr/rrp
      p(1:ndim) = r(1:ndim) +beta*p(1:ndim)
!.....Store previous step values
      xp(1:ndim) = x(1:ndim)
    enddo

!.....In case not converged within nstp, warn something...
    ierr= -1
    return
  end subroutine cg
!=======================================================================
  subroutine set_paramsdir_Coulomb(dname)
!
!  Accessor routine to set paramsdir.
!
    implicit none
    character(len=*),intent(in):: dname

    paramsdir = trim(dname)
    return
  end subroutine set_paramsdir_Coulomb
!=======================================================================
  subroutine set_params_Coulomb(ndimp,prms_in,ctype,specorder,iprint)
!
!  Accessor routine to set Coulomb parameters from outside.
!  This is supposed to be called only on serial run.
!
    integer,intent(in):: ndimp,iprint
    real(8),intent(in):: prms_in(ndimp)
    character(len=*),intent(in):: ctype
    character(len=3),intent(in):: specorder(nspmax)

    integer:: isp,jsp,ns,inc,ipr,myid,mpiw,maxisp

    if( index(ctype,'BVS').ne.0 ) then
!!$      if( ctype(4:4).eq.'2' .or. ctype(4:4).eq.'3' ) then
!.....Not only fbvs, but also rad_bvs are given from fitpot
      ipr = max(0,iprint-10)
      myid = 0
      mpiw = -1
      call read_params(myid,mpiw,ipr,specorder)
      lprmset_Coulomb = .true.

      fbvs = prms_in(1)
      inc = 1
      do isp=1,nspmax
        if( specorder(isp).eq.'x' ) cycle
        inc = inc + 1
        if( inc.gt.ndimp ) then
          print *,'ERROR @set_parmas_Coulomb: inc.gt.ndimp !!!'
          stop
        endif
        rad_bvs(isp) = prms_in(inc)
      enddo

!.....Reset screening length
      do isp=1,nspmax
        if( vid_bvs(isp).eq.0d0 ) cycle
        do jsp=1,nspmax
          if( vid_bvs(jsp).eq.0d0 ) cycle
          rho_scr(isp,jsp) = fbvs*(rad_bvs(isp)+rad_bvs(jsp))
        enddo
      enddo

!!$      else
!.....As of 190819, ctype==BVS means that only fbvs is to be given from fitpot.
!.....Need to read in.params.XX file before going further
!!         ipr = 0
!!         myid = 0
!!         mpiw = -1
!!         call read_paramsx(myid,mpiw,ipr,specorder)
!!         params_read = .true.
!!         lprmset_Coulomb = .true.
!! 
!! !.....set fbvs
!!         fbvs = prms_in(1)
!! !!$      print *,'fbvs=',fbvs
!! !.....Reset screening length
!!         do isp=1,nspmax
!!           if( vid_bvs(isp).eq.0d0 ) cycle
!!           do jsp=1,nspmax
!!             if( vid_bvs(jsp).eq.0d0 ) cycle
!!             rho_scr(isp,jsp) = fbvs*(rad_bvs(isp)+rad_bvs(jsp))
!! !!$          print *,'isp,jsp,rho_scr=',isp,jsp,rho_scr(isp,jsp)
!!           enddo
!!         enddo
!!       endif

    else if( trim(ctype).eq.'fpc' ) then
!.....As of 190818, only one parameter is passed to this routine
!     that scale the original charges per species obtained from in.params.Coulomb
!.....Need to read in.params.Coulomb
      if( .not. params_read ) then
        ipr = 0
        myid = 0
        mpiw = -1
        call read_params(myid,mpiw,ipr,specorder)
        params_read = .true.
      endif
      lprmset_Coulomb = .true.

      do isp=1,nspmax
        if( ispflag(isp) ) then
          schg(isp) = schg0(isp)*prms_in(1)
        endif
      enddo

    endif   ! ctype

  end subroutine set_params_Coulomb
!=======================================================================
  subroutine gradw_Coulomb(namax,natm,nb,tag,ra,chg,nnmax &
       ,h,rc,lspr,epot,iprint,ndimp,gwe,gwf,gws &
       ,lematch,lfmatch,lsmatch,iprm0,myid,mpi_world,specorder)
!
!  Gradient w.r.t. weights.
!  Note: This routine is always called in single run,
!  thus no need of parallel implementation.
!  - iprm0: The starting point -1 in parameter array for this FF.
!
    implicit none
    include "./params_unit.h"
    integer,intent(in):: namax,natm,nb,nnmax,iprint,iprm0 &
         ,myid,mpi_world
    integer,intent(in):: lspr(0:nnmax,namax)
    real(8),intent(in):: ra(3,namax),h(3,3),rc,tag(namax)
    real(8),intent(inout):: epot,chg(namax)
    integer,intent(in):: ndimp
    real(8),intent(inout):: gwe(ndimp),gwf(3,ndimp,natm),gws(6,ndimp)
    logical,intent(in):: lematch,lfmatch,lsmatch
    character(len=3),intent(in):: specorder(nspmax)

    integer:: i,j,k,isp,jsp,ne,nf,ns,maxisp,ixyz,jxyz,jj
    real(8):: rc2,xi(3),xj(3),xij(3),rij(3),dxdi(3),dxdj(3) &
         ,dedr,ddedrho,dedrho,dij,dij2,diji,fac,rhoij,qi,qj,sqpi &
         ,terfc,terfcc,tmp,vrc
    real(8),allocatable:: ge_rho(:),gf_rho(:,:,:),gs_rho(:,:)

    if( .not.allocated(ge_rho) ) then
      allocate(ge_rho(nspmax),gs_rho(nspmax,6),gf_rho(nspmax,3,natm))
      call accum_mem('force_Coulomb',8*(size(ge_rho)+size(gs_rho)+size(gf_rho)))
    endif
    if( size(gf_rho).ne.nspmax*3*natm ) then
      if( allocated(gf_rho) ) then
        call accum_mem('force_Coulomb',-8*size(gf_rho))
        deallocate(gf_rho)
      endif
      allocate(gf_rho(nspmax,3,natm))
      call accum_mem('force_Coulomb',8*size(gf_rho))
    endif

    call set_charge_BVS(natm,nb,tag,chg,myid,mpi_world,iprint,specorder)

!.....Max isp
    maxisp = 0
    do i=1,natm
      maxisp = max(maxisp,int(tag(i)))
    enddo

    rc2 = rc*rc

    ge_rho(1:nspmax) = 0d0
    gf_rho(1:nspmax,1:3,1:natm) = 0d0
    gs_rho(1:nspmax,1:6) = 0d0

!.....Loop over resident atoms
    do i=1,natm
      xi(1:3) = ra(1:3,i)
      isp = int(tag(i))
      if( .not. ispflag(isp) ) cycle
      qi = chg(i)
      do jj=1,lspr(0,i)
        j= lspr(jj,i)
        if( j.eq.0 ) exit
        if( j.le.i ) cycle
        jsp = int(tag(j))
!.....Check if two species interact
        if( .not. interact(isp,jsp) ) cycle
        xj(1:3) = ra(1:3,j)
        xij(1:3) = xj(1:3) -xi(1:3)
        rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij2 = rij(1)**2 +rij(2)**2 +rij(3)**2
        if( dij2.gt.rc2 ) cycle
        dij= sqrt(dij2)
        diji = 1d0/dij
        dxdi(1:3)= -rij(1:3)*diji
        dxdj(1:3)=  rij(1:3)*diji
        rhoij = rho_scr(isp,jsp)
        terfc = erfc(dij/rhoij)
        terfcc = erfc(rc/rhoij)
        qj = chg(j)
        if( j.le.natm ) then
          fac = 1d0
        else
          fac = 0.5
        endif
        if( lematch ) then
          vrc = acc*qi*qj/rc *terfcc
!.....Derivative of potential energy wrt {w}
          dedrho = dvdrho(dij,rhoij,qi,qj) -dvdrho(rc,rhoij,qi,qj) &
               -(dij-rc)*ddvdrho(rc,rhoij,qi,qj)
          ge_rho(isp) = ge_rho(isp) +dedrho*fac
          ge_rho(jsp) = ge_rho(jsp) +dedrho*fac
        endif
!.....Pre-compute some factors required in force and stress derivatives
        if( lfmatch .or. lsmatch ) then
          dedr = -acc *qi*qj*diji *(terfc*diji +2d0/rhoij *sqpi *exp(-(dij/rhoij)**2))
          ddedrho = ddvdrho(dij,rhoij,qi,qj) -ddvdrho(rc,rhoij,qi,qj)
        endif
        if( lfmatch ) then
          gf_rho(isp,1:3,i) = gf_rho(isp,1:3,i) -dxdi(1:3)*ddedrho
          gf_rho(jsp,1:3,i) = gf_rho(jsp,1:3,i) -dxdi(1:3)*ddedrho
          if( j.le.natm ) then
            gf_rho(isp,1:3,j) = gf_rho(isp,1:3,j) -dxdj(1:3)*ddedrho
            gf_rho(jsp,1:3,j) = gf_rho(jsp,1:3,j) -dxdj(1:3)*ddedrho
          endif
        endif
        if( lsmatch ) then
          do ixyz=1,3
            do jxyz=1,3
              k = ivoigt(ixyz,jxyz)
              gs_rho(isp,k) = gs_rho(isp,k) &
                   -fac *rij(ixyz) *(-dxdi(jxyz)) *ddedrho
              gs_rho(jsp,k) = gs_rho(jsp,k) &
                   -fac *rij(ixyz) *(-dxdi(jxyz)) *ddedrho
            enddo
          enddo
        endif
      enddo
    enddo

!.....Tidy up gradient arrays
    if( lematch ) then
      ne = iprm0
      do isp=1,maxisp
        if( .not. ispflag(isp) ) cycle
        ne = ne + 1
        gwe(ne) = gwe(ne) +ge_rho(isp)
      enddo
    endif
    if( lfmatch ) then
      do i=1,natm
        nf = iprm0
        do isp=1,maxisp
          if( .not. ispflag(isp) ) cycle
          nf = nf + 1
          gwf(1:3,nf,i) = gwf(1:3,nf,i) +gf_rho(isp,1:3,i)
        enddo
      enddo
    endif
    if( lsmatch ) then
      do i=1,natm
        ns = iprm0
        do isp=1,maxisp
          if( .not.ispflag(isp) ) cycle
          ns = ns + 1
          gws(1:6,ns) = gws(1:6,ns) +gs_rho(isp,1:6)
        enddo
      enddo
    endif
    return
  end subroutine gradw_Coulomb
!=======================================================================
  function dvdrho(r,rho,qi,qj)
!
!  Derivative screened Coulomb wrt rho.
!
    real(8),intent(in):: r,rho,qi,qj
    real(8):: dvdrho

    dvdrho = acc *qi*qj*2d0 /rho**2 /sqrt(pi) *exp(-(r/rho)**2)
    return
  end function dvdrho
!=======================================================================
  function ddvdrho(r,rho,qi,qj)
!
!  Derivative of dvdr wrt rho.
!  See the memo on 2018-04-11.
!
    real(8),intent(in):: r,rho,qi,qj
    real(8):: ddvdrho

    ddvdrho = -4d0*acc*qi*qj*r *exp(-(r/rho)**2) /sqrt(pi) /rho**4
    return
  end function ddvdrho
!=======================================================================
  subroutine update_vauxq(vauxq)
!
!  Update velocity of auxq
!
    use pmdvars,only: dt,namax,natm
    real(8),intent(inout):: vauxq(namax)

    vauxq(1:natm) = vauxq(1:natm) +0.5d0 *dt *aauxq(1:natm)
    return
  end subroutine update_vauxq
!=======================================================================
  subroutine update_auxq(auxq,vauxq)
!
!  Update auxq
!
    use pmdvars,only: dt,namax,natm
    real(8),intent(in):: vauxq(namax)
    real(8),intent(inout):: auxq(namax)

    auxq(1:natm) = auxq(1:natm) +dt *vauxq(1:natm)
    call impose_qtot(auxq)

    return
  end subroutine update_auxq
!=======================================================================
  subroutine get_aauxq(chg,auxq)
!
!  Compute forces on auxqs.
!  d^2 auxq/dt^2 = omg^2 *(chg -auxq)
!
    use pmdvars,only: dt,namax,natm
    real(8),intent(in):: chg(namax),auxq(namax)

    if( .not.allocated(aauxq) ) then
      allocate(aauxq(namax))
    else if( size(aauxq).ne.namax ) then
      call accum_mem('force_Coulomb',-8*size(aauxq))
      deallocate(aauxq)
      allocate(aauxq(namax))
      call accum_mem('force_Coulomb',8*size(aauxq))
    endif

    auxomg2 = omg2dt2 /dt**2
    aauxq(1:natm) = auxomg2 *(chg(1:natm) -auxq(1:natm))
  end subroutine get_aauxq
!=======================================================================
  subroutine impose_qtot(chg)
!
!  Impose total charge to qtot.
!
    use pmdvars,only: namax,natm,myid_md,mpi_md_world,iprint,ntot
    implicit none
    include 'mpif.h'
    include './const.h'
    real(8),intent(inout):: chg(namax)

    integer:: i,ierr
    real(8):: ql,qg,dq,qdist

    ql = 0d0
    do i=1,natm
      ql = ql +chg(i)
    enddo
    qg = 0d0
    call mpi_allreduce(ql,qg,1,mpi_real8,mpi_sum,mpi_md_world,ierr)
    dq = qg -qtot_qeq
    qdist = dq /ntot
    do i=1,natm
      chg(i) = chg(i) -qdist
    enddo
    return
  end subroutine impose_qtot
!=======================================================================
  subroutine bacopy_chg_fixed(chg)
!-----------------------------------------------------------------------
!  Exchange only chg of boundary-atom
!  This doesnt search using position, just send & recv data of atoms
!    which were listed by 'bacopy'.
!  This is called from dampopt in force_common.F90.
!-----------------------------------------------------------------------
    use pmdvars,only: nx,ny,nz,lsb,lsex,namax,natm,nbmax,nb,nn, &
         myid_md,mpi_md_world,myparity,lsrc,nex,boundary
    use pmdmpi,only: nid2xyz
    implicit none
    include 'mpif.h'
    real(8),intent(inout):: chg(namax)

    integer:: i,j,kd,kdd,ku,ierr,iex,ix,iy,iz
    integer:: inode,nsd,nrc,nbnew
    real(8),save,allocatable:: dbuf(:),dbufr(:)
    logical,save:: l1st=.true.

    if( l1st ) then
      allocate(dbuf(nbmax),dbufr(nbmax))
      l1st=.false.
    endif

    if( size(dbuf).ne.nbmax ) then
      call accum_mem('force_Coulomb',-8*(size(dbuf)+size(dbufr)))
      deallocate(dbuf,dbufr)
      allocate(dbuf(nbmax),dbufr(nbmax))
      call accum_mem('force_Coulomb',8*(size(dbuf)+size(dbufr)))
    endif

    call nid2xyz(myid_md,ix,iy,iz)

    nbnew= 0

!c-----calculate the cut-off lengths
!      do kd=1,3
!        asgm= dsqrt(sgm(1,i)**2 +sgm(2,i)**2 +sgm(3,i)**2)
!        rcv(kd)= rc*asgm/vol
!        nex(kd)= int(rcv(kd)) +1
!      enddo

!-----loop over x, y, & z directions
    do kd=1,3

      if( nex(kd).gt.1 ) then
        do kdd= -1,0
          ku= 2*kd+kdd
          do i=1,lsb(0,ku)
            j= lsb(i,ku)
            iex= lsex(i,ku)
            chg(natm+nbnew+i)= chg(j)
          enddo
          nbnew= nbnew +lsb(0,ku)
        enddo
      else
        do kdd= -1,0
          ku=2*kd+kdd
          inode=nn(ku)
          nsd=lsb(0,ku)
          nrc=lsrc(ku)

!---------Exchange ra and tag
          do i=1,nsd
            j= lsb(i,ku)
            dbuf(i)  = chg(j)
          enddo
          call mespasd(inode,myparity(kd),dbuf,dbufr,nsd,nrc,23 &
               ,mpi_md_world)
          do i=1,nrc
            chg(natm+nbnew+i) = dbufr(i)
          enddo

          call mpi_barrier(mpi_md_world,ierr)
          nbnew=nbnew +nrc
200       continue
        enddo
      endif

100   continue
    enddo

  end subroutine bacopy_chg_fixed
!=======================================================================
  subroutine bacopy_auxq_fixed(auxq,vauxq)
!-----------------------------------------------------------------------
!  Exchange only auxq and vauxq of boundary atoms
!  This doesnt search using position, just send & recv data of atoms
!    which were listed by 'bacopy'.
!-----------------------------------------------------------------------
    use pmdvars,only: nx,ny,nz,lsb,lsex,namax,natm,nbmax,nb,nn, &
         myid_md,mpi_md_world,myparity,lsrc,nex,boundary
    use pmdmpi,only: nid2xyz
    implicit none
    include 'mpif.h'

    real(8),intent(inout):: auxq(namax),vauxq(namax)

    integer:: i,j,kd,kdd,ku,ierr,iex,ix,iy,iz
    integer:: inode,nsd,nrc,nbnew
    real(8),save,allocatable:: dbuf(:,:),dbufr(:,:)
    integer,parameter:: ndim = 2
    logical,save:: l1st=.true.

    if( l1st ) then
      allocate(dbuf(ndim,nbmax),dbufr(ndim,nbmax))
      l1st=.false.
    endif

    if( size(dbuf).ne.nbmax ) then
      call accum_mem('force_Coulomb',-8*(size(dbuf)+size(dbufr)))
      deallocate(dbuf,dbufr)
      allocate(dbuf(ndim,nbmax),dbufr(ndim,nbmax))
      call accum_mem('force_Coulomb',8*(size(dbuf)+size(dbufr)))
    endif
    
    call nid2xyz(myid_md,ix,iy,iz)

    nbnew= 0

!-----loop over x, y, & z directions
    do kd=1,3

      if( nex(kd).gt.1 ) then
        do kdd= -1,0
          ku= 2*kd+kdd
          do i=1,lsb(0,ku)
            j= lsb(i,ku)
            iex= lsex(i,ku)
            auxq(natm+nbnew+i)= auxq(j)
            vauxq(natm+nbnew+i)= vauxq(j)
          enddo
          nbnew= nbnew +lsb(0,ku)
        enddo
      else
        do kdd= -1,0
          ku=2*kd+kdd
          inode=nn(ku)
          nsd=lsb(0,ku)
          nrc=lsrc(ku)

!---------Exchange ra and tag
          do i=1,nsd
            j= lsb(i,ku)
            dbuf(1,i)  = auxq(j)
            dbuf(2,i)  = vauxq(j)
          enddo
          call mespasd(inode,myparity(kd),dbuf,dbufr,nsd*ndim &
               ,nrc*ndim,23,mpi_md_world)
          do i=1,nrc
            auxq(natm+nbnew+i) = dbufr(1,i)
            vauxq(natm+nbnew+i) = dbufr(2,i)
          enddo

          call mpi_barrier(mpi_md_world,ierr)
          nbnew=nbnew +nrc
200       continue
        enddo
      endif

100   continue
    enddo

  end subroutine bacopy_auxq_fixed
!=======================================================================
end module Coulomb
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd lib"
!     End:
