module Coulomb
!-----------------------------------------------------------------------
!                     Last modified: <2018-01-15 18:01:47 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
!  Parallel implementation of Coulomb potential
!  ifcoulomb == 1: screened Coulomb potential
!            == 2: Ewald Coulomb potential
!            == 3: variable-charge with Gaussian distribution
!
!  For screened Coulomb potential:
!    - Ref. Adams & Rao, Phys. Status Solidi A 208, No.8 (2011)
!    - No cutoff treatment
!  For Ewald Coulomb potential:
!    - ...
!-----------------------------------------------------------------------
  implicit none
  save

  character(len=128):: paramsdir = '.'
  character(len=128),parameter:: paramsfname = 'in.params.Coulomb'
  real(8),parameter:: pi = 3.14159265398979d0

  character(len=128):: cchgs,cdist,cterms
!.....Input Keywords: charges, charge_dist, terms
!.....Keywords for charges: fixed, fixed_bvs, variable/qeq
!.....Keywords for charge_dist: point, gaussian
!.....Keywords for terms: full, short/screened, long

  integer,parameter:: ioprms = 20
!.....Coulomb's constant, acc = 1.0/(4*pi*epsilon0)
  real(8),parameter:: acc  = 14.3998554737d0
!.....permittivity of vacuum
  real(8),parameter:: eps0 = 0.00552634939836d0  ! e^2 /Ang /eV


  integer,parameter:: msp = 9
  integer:: nsp
!.....Species charge
  real(8):: schg(msp)
!  logical,allocatable:: interact(:,:)
  logical:: interact(msp,msp)
!.....ideal valence charges of species
  real(8):: vid_bvs(msp)
!  real(8),allocatable:: vid_bvs(:)
!.....principal quantum numbers of species
  integer:: npq_bvs(msp)
!  integer,allocatable:: npq_bvs(:)
!.....covalent radius
  real(8):: rad_bvs(msp)
!  real(8),allocatable:: rad_bvs(:)
!.....screening length
  real(8):: rho_bvs(msp,msp)
!  real(8),allocatable:: rho_bvs(:,:)
  real(8),parameter:: fbvs = 0.74d0

!.....charge threshold for Coulomb interaction [default: 0.01]
  real(8),parameter:: qthd = 1d-12

!.....Gaussian width of Ewald sum [default: 1.0 (Angstrom)]
  real(8):: sgm_ew = 1d0
  real(8):: sgm(msp)
!.....Accuracy controlling parameter for Ewald sum
!.....See, http://www.jncasr.ac.in/ccms/sbs2007/lecturenotes/5day10nov/SBS_Ewald.pdf
!.....Exp(-pacc) = 1e-7 when pacc= 18.0
!  real(8),parameter:: pacc   = 18d0
  real(8),parameter:: pacc = 9.21d0  ! exp(-pacc) = 1e-4
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

!.....Variable-charge potential variables
  real(8):: vcg_chi(msp),vcg_jii(msp),vcg_e0(msp),vcg_sgm(msp) &
       ,qlower(msp),qupper(msp)
  real(8),parameter:: vcg_lambda = 0.5d0

contains
  subroutine initialize_coulomb(natm,nspin,tag,chg,chi, &
       myid,mpi_md_world,ifcoulomb,iprint,h,rc,lvc)
!
!  Allocate and initialize parameters to be used.
!
    include "mpif.h"
    integer,intent(in):: myid,mpi_md_world,natm,nspin &
         ,ifcoulomb,iprint
    real(8),intent(in):: tag(natm),rc,h(3,3),chg(natm)
    real(8),intent(inout):: chi(natm)
    logical,intent(in):: lvc 

    integer:: i,ierr,nspl

!!$    print *,'ifcoulomb @initialize_coulomb = ',ifcoulomb

    
!.....Get umber of species
    nsp = nspin
!!$    nspl = 0
!!$    do i=1,natm
!!$      nspl = max(int(tag(i)),nspl)
!!$    enddo
!!$    call mpi_allreduce(nspl,nsp,1,mpi_integer,mpi_max &
!!$         ,mpi_md_world,ierr)

    if( ifcoulomb.eq.1 ) then  ! screened Coulomb
      call read_params(myid,mpi_md_world,ifcoulomb,iprint,lvc)
    else
      if( lvc ) then
!.....Variable-charge Coulomb with Gaussian distribution charges
!     which ends-up long-range-only Ewald summation
        call init_vc_Ewald(myid,mpi_md_world,ifcoulomb,iprint,h,rc,&
             natm,tag,chi,chg,lvc)
      else
        call init_fc_Ewald(h,rc,myid)
      endif
    endif

  end subroutine initialize_coulomb
!=======================================================================
  subroutine initialize_coulombx(natm,nspin,tag,chg,chi, &
       myid,mpi_md_world,ifcoulomb,iprint,h,rc,lvc)
!
!  Allocate and initialize parameters to be used.
!
    include "mpif.h"
    integer,intent(in):: myid,mpi_md_world,natm,nspin &
         ,ifcoulomb,iprint
    real(8),intent(in):: tag(natm),rc,h(3,3),chg(natm)
    real(8),intent(inout):: chi(natm)
    logical,intent(in):: lvc 

    integer:: i,ierr,nspl

!!$    print *,'ifcoulomb @initialize_coulomb = ',ifcoulomb

    
!.....Get umber of species
    nsp = nspin

    call read_paramsx(myid,mpi_md_world,iprint)

    if( trim(cchgs).eq.'variable' ) then
!.....Variable-charge Coulomb with Gaussian distribution charges
!     which ends-up long-range-only Ewald summation
      call init_vc_Ewald(myid,mpi_md_world,ifcoulomb,iprint,h,rc,&
           natm,tag,chi,chg,lvc)
    else
      call init_fc_Ewald(h,rc,myid)
    endif

  end subroutine initialize_coulombx
!=======================================================================
  subroutine init_fc_Ewald(h,rc,myid)
!
!  Ewald sum with fixed charge.
!
    implicit none 
    integer,intent(in):: myid
    real(8),intent(in):: h(3,3),rc

    integer:: i,ik,k1,k2,k3,isp
    real(8):: bk1(3),bk2(3),bk3(3),bk(3),bb2
    real(8),external:: absv

    sgm_ew = rc/sqrt(2d0*pacc)
    sgm(1:msp) = sgm_ew
    bkmax  = 2d0*pacc /rc
    if( myid.eq.0 ) then
      write(6,'(/,a)') ' Ewald sum parameters:'
      write(6,'(a,f12.4)') '   1/(4*pi*eps0)        = ', acc
      write(6,'(a,f12.4)') '   Accuracy parameter p = ', pacc
      write(6,'(a,f12.4)') '   Gaussian width sgm   = ', sgm_ew
      write(6,'(a,f12.4)') '   real-space cutoff    = ', rc
      write(6,'(a,f12.4)') '   k-space cutoff       = ', bkmax
    endif
    call get_recip_vectors(h)
    if( myid.eq.0 ) then
      write(6,'(a)') ' Reciprocal vectors:'
      write(6,'(a,3es12.3)') '   b1 = ',b1(1:3)
      write(6,'(a,3es12.3)') '   b2 = ',b2(1:3)
      write(6,'(a,3es12.3)') '   b3 = ',b3(1:3)
    endif
!.....kmax# is constant during MD run even if h-matrix can change...
    call setup_kspace()
    if( myid.eq.0 ) then
      write(6,'(a)') ' Number of k-points for Ewald sum:'
      write(6,'(a,i8)') '   kmax1 = ',kmax1
      write(6,'(a,i8)') '   kmax2 = ',kmax2
      write(6,'(a,i8)') '   kmax3 = ',kmax3
      write(6,'(a,i8)') '   total = ',nk
    endif
    if( allocated(qcos) ) deallocate(qcosl,qcos,qsinl,qsin,pflr)
    allocate(qcosl(nk),qcos(nk),qsinl(nk),qsin(nk),pflr(nk,msp))
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
          bb2 = absv(3,bk)
          bb2 = bb2*bb2
          pflr(ik,1:nsp)= 4d0 *pi /bb2 *exp(-0.5d0 *sgm_ew**2 *bb2)
        enddo
      enddo
    enddo

  end subroutine init_fc_Ewald
!=======================================================================
  subroutine init_vc_Ewald(myid,mpi_world,ifcoulomb,iprint,h,rc,&
       natm,tag,chi,chg,lvc)
!
!  Since variable-charge potential with Gaussian distribution charges
!  is mostly identical to the long-range part of Ewald summation,
!  the code here is mostly the same as Ewald sum except that vcGaussian
!  requires differnt Gaussian width for different species that should be
!  read from input file in.params.Coulomb.
!
    implicit none
    integer,intent(in):: myid,mpi_world,ifcoulomb,iprint
    integer,intent(in):: natm
    real(8),intent(in):: h(3,3),rc,tag(natm),chg(natm)
    real(8),intent(inout):: chi(natm)
    logical,intent(in):: lvc

    integer:: i,isp,ik,k1,k2,k3,is
    real(8):: bk1(3),bk2(3),bk3(3),bk(3),bb2,sgm_min
    real(8),external:: absv

    call read_params(myid,mpi_world,ifcoulomb,iprint,lvc)

    do i=1,natm
      is = int(tag(i))
      chi(i) = vcg_chi(is)
    enddo

    if( ifcoulomb.eq.2 ) then
      sgm_ew = rc/sqrt(2d0*pacc)
      sgm(1:msp) = sgm_ew
      bkmax  = 2d0*pacc /rc
      if( myid.eq.0 .and. iprint.ne.0 ) then
        print *,'Sigmas are overwritten by rc/sqrt(2*pacc),'
        print *,'since the full Ewald cannot treat species-dependent sigma.'
      endif
    else
      sgm(1:msp) = vcg_sgm(1:msp)
      sgm_min = 1d+30
      do i=1,nsp
        sgm_min = min(sgm(i),sgm_min)
      enddo
      bkmax  = sqrt(2d0*pacc) /sgm_min
    endif
    if( myid.eq.0 .and. iprint.ne.0 ) then
      write(6,'(/,a)') ' Ewald sum parameters:'
      write(6,'(a,f12.4)') '   1/(4*pi*eps0)        = ', acc
      write(6,'(a,f12.4)') '   Accuracy parameter p = ', pacc
      if( ifcoulomb.eq.2 ) then
        write(6,'(a,f12.4)') '   Gaussian width sgm   = ', sgm_ew
      else
        write(6,'(a)')       '   Gaussian simgas:'
        do i=1,nsp
          write(6,'(a,f12.4)') '     ',sgm(i)
        enddo
        write(6,'(a,f12.4)') '   Minimum sigma        = ', sgm_min
      endif
      write(6,'(a,f12.4)') '   real-space cutoff    = ', rc
      write(6,'(a,f12.4)') '   k-space cutoff       = ', bkmax
    endif
    call get_recip_vectors(h)
    if( myid.eq.0 .and. iprint.ne.0 ) then
      write(6,'(a)') ' Reciprocal vectors:'
      write(6,'(a,3es12.3)') '   b1 = ',b1(1:3)
      write(6,'(a,3es12.3)') '   b2 = ',b2(1:3)
      write(6,'(a,3es12.3)') '   b3 = ',b3(1:3)
    endif
!.....kmax# is constant during MD run even if h-matrix can change...
    call setup_kspace()
    if( myid.eq.0 .and. iprint.ne.0 ) then
      write(6,'(a)') ' Number of k-points for Ewald sum:'
      write(6,'(a,i8)') '   kmax1 = ',kmax1
      write(6,'(a,i8)') '   kmax2 = ',kmax2
      write(6,'(a,i8)') '   kmax3 = ',kmax3
      write(6,'(a,i8)') '   total = ',nk
      write(6,*) ''
    endif
    if( allocated(qcos) ) deallocate(qcos,qsin,qcosl,qsinl,pflr)
    allocate(qcosl(nk),qcos(nk),qsinl(nk),qsin(nk),pflr(nk,msp))
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
          bb2 = absv(3,bk)
          bb2 = bb2*bb2
          do isp=1,nsp
            pflr(ik,isp)= 4d0 *pi /bb2 *exp(-0.5d0 *sgm(isp)**2 *bb2)
          enddo
        enddo
      enddo
    enddo
    
  end subroutine init_vc_Ewald
!=======================================================================
  subroutine read_params(myid,mpi_world,ifcoulomb,iprint,lvc)
!
!  Read parameters
!
    include "mpif.h"
    integer,intent(in):: myid,mpi_world,ifcoulomb,iprint
    character(len=128):: cline,c1st,fname
    character(len=5):: cname
    logical,intent(in):: lvc

    real(8):: vid,rad,dchi,djii,dsgm,de0
    integer:: npq,isp,jsp,ierr,mode

    if( ifcoulomb.eq.1 ) then  ! screened_bvs
      call read_params_sc(myid,mpi_world,ifcoulomb,iprint)
    else if( lvc ) then  ! variable-charge
      call read_params_vc(myid,mpi_world,ifcoulomb,iprint)
    else
      print *,'ERROR: something wrong about force_Coulomb setting.'
      call mpi_finalize(ierr)
      stop
    endif  ! lvc

  end subroutine read_params
!=======================================================================
  subroutine read_params_sc(myid,mpi_world,ifcoulomb,iprint)
!
!  Read params for screened_Coulomb.
!
    implicit none 
    include "mpif.h"
    integer,intent(in):: myid,mpi_world,ifcoulomb,iprint
    
    integer:: mode,ierr,isp,jsp,npq
    real(8):: vid,rad
    character(len=128):: cline,c1st,fname
    character(len=5):: cname
    
!!$    if( allocated(rad_bvs) ) deallocate(rad_bvs,npq_bvs,vid_bvs,rho_bvs)
!!$    allocate(rad_bvs(msp),npq_bvs(msp),vid_bvs(msp) &
!!$         ,rho_bvs(msp,msp))
    if( myid.eq.0 ) then
      fname = trim(paramsdir)//'/'//trim(paramsfname)
      open(ioprms,file=trim(fname),status='old')
      mode = 1
      interact(1:msp,1:msp) = .false.
!.....1st line for check the Coulomb computation type
      read(ioprms,*) c1st
      if( trim(c1st).ne.'screened_bvs' ) then
        write(6,*) 'Error@read_params_sc: something wrong '&
             //'with the Coulomb type: '//trim(c1st)
        stop
      endif
      if( iprint.ne.0 ) then
        write(6,'(/,a)') ' Screened Coulomb parameters:'
      endif
      vid_bvs(1:msp)= 0d0
      do while(.true.)
        read(ioprms,*,end=10) cline
        if( cline(1:1).eq.'#' .or. cline(1:1).eq.'!' ) cycle
        if( trim(cline).eq.'interactions' ) then
          mode = 2
          cycle
        endif
        if( mode.eq.1 ) then
          backspace(ioprms)
          read(ioprms,*) isp, cname, vid, rad, npq
          if( isp.gt.msp ) then
            write(6,*) ' Warning @read_params: since isp is greater than nsp,'&
                 //' skip reading the line.'
            cycle
          endif
          if( iprint.ne.0 ) then
            write(6,'(a,i3,a5,2f7.3,i4)') ' isp,cname,vid,rad,npq =' &
                 ,isp,trim(cname),vid,rad,npq
          endif
          vid_bvs(isp) = vid
          rad_bvs(isp) = rad
          npq_bvs(isp) = npq
        else if( mode.eq.2 ) then
          backspace(ioprms)
          read(ioprms,*) isp, jsp
          interact(isp,jsp) = .true.
          interact(jsp,isp) = interact(isp,jsp)
!!$            if( iprint.gt.0 ) then
!!$              write(6,'(a,2i5,l3)') ' isp,jsp,interact= ',isp,jsp,interact(isp,jsp)
!!$            endif
        endif
      enddo
10    close(ioprms)
      if( iprint.ne.0 ) then
        write(6,'(a)') ' Finished reading '//trim(fname)
        write(6,*) ''
      endif

!.....Set screening length
      do isp=1,nsp
        do jsp=1,nsp
          rho_bvs(isp,jsp) = fbvs*(rad_bvs(isp)+rad_bvs(jsp))
!!$            rho_bvs(isp,jsp) = 2d0
          if( iprint.gt.0 .and. interact(isp,jsp) .and. jsp.ge.isp ) then
            write(6,'(a,2i5,f10.4)') ' isp,jsp,rho_bvs= ',isp,jsp,rho_bvs(isp,jsp)
          endif
        enddo
      enddo
    endif  ! myid

    call mpi_bcast(vid_bvs,msp,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(rad_bvs,msp,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(npq_bvs,msp,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(rho_bvs,msp*msp,mpi_real8 &
         ,0,mpi_world,ierr)
    call mpi_bcast(interact,msp*msp,mpi_logical,0,mpi_world,ierr)
!.....end of screend_bvs
  end subroutine read_params_sc
!=======================================================================
  subroutine read_params_vc(myid,mpi_world,ifcoulomb,iprint)
!
!  Read params for variable-charge.
!
    implicit none 
    include "mpif.h"
    integer,intent(in):: myid,mpi_world,ifcoulomb,iprint

    integer:: isp,ierr
    real(8):: dchi,djii,sgmt,de0,qlow,qup
    character(len=128):: cline,c1st,fname
    character(len=5):: cname

!!$      if( allocated(vcg_chi) ) deallocate(vcg_chi,vcg_jii,sgm,vcg_e0)
!!$      allocate( vcg_chi(msp), vcg_jii(msp), sgm(msp), vcg_e0(msp))
    if( myid.eq.0 ) then
      fname = trim(paramsdir)//'/'//trim(paramsfname)
      open(ioprms,file=trim(fname),status='old')
!.....1st line for check the Coulomb computation type
      read(ioprms,*) c1st
      if( trim(c1st).ne.'variable_charge' ) then
        write(6,'(/,a)') 'Error@read_params: ifcoulomb does not match '&
             //'with Coulomb type: '//trim(c1st)
        stop
      endif
      if( iprint.gt.0 ) write(6,'(/,a)') ' Variable_charge parameters:'
      do while(.true.)
        read(ioprms,*,end=20) cline
        if( cline(1:1).eq.'#' .or. cline(1:1).eq.'!' ) cycle
        backspace(ioprms)
        read(ioprms,*,end=20) isp, cname, dchi,djii,sgmt,de0,qlow,qup
        if( isp.gt.nsp .and. iprint.gt.0 ) then
          write(6,'(a,i2)') ' [WARNING] isp.gt.nsp !!!  isp = ',isp
        else
          vcg_chi(isp) = dchi
          vcg_jii(isp) = djii
          vcg_e0(isp) = de0
          vcg_sgm(isp) = sgmt
          qlower(isp) = qlow
          qupper(isp) = qup
          if( iprint.gt.0 ) then
            write(6,'(a,i3,a3,3f10.4,2f5.1)') '   isp,name,chi,Jii,sgm,e0,qlower,qupper = ', &
                 isp,trim(cname),dchi,djii,vcg_sgm(isp),de0,qlow,qup
          endif
        endif
      enddo  ! do while
20    close(ioprms)
    endif  ! myid
!!$    call mpi_bcast(sgm,nsp,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(vcg_chi,nsp,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(vcg_jii,nsp,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(vcg_e0,nsp,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(vcg_sgm,nsp,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(qlower,nsp,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(qupper,nsp,mpi_real8,0,mpi_world,ierr)

  end subroutine read_params_vc
!=======================================================================
  subroutine read_paramsx(myid,mpi_world,iprint)
!
!  Read parameter file for any Coulomb potential.
!
    include "mpif.h"
    integer,intent(in):: myid,mpi_world,iprint

    integer,external:: num_data
    
    character(len=128):: cmode,cline,ctmp,fname
    character(len=3):: cname
    integer:: i,ierr,jerr,isp,jsp,npq
    real(8):: chgi,vid,rad,dchi,djii,sgmt,de0,qlow,qup

    if( myid.eq.0 ) then
!.....Initialization
      interact(1:msp,1:msp) = .true.
!.....File name
      fname = trim(paramsdir)//'/'//trim(paramsfname)
      open(ioprms,file=trim(fname),status='old')
!.....Start reading
      do while(.true.)
        read(ioprms,*,end=10) cline
        if( num_data(cline,' ').eq.0 ) cycle
        if( cline(1:1).eq.'!' .or. cline(1:1).eq.'#' ) cycle
!.....Detect keyword
        if( trim(cline).eq.'charges' ) then
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
          cycle
        else if( trim(cline).eq.'charge_dist' ) then
          cmode = 'charge_dist'
          backspace(ioprms)
          read(ioprms,*) ctmp, cdist
          if(  trim(cdist).ne.'point' .and. &
               trim(cdist).ne.'gaussian' ) then
            print *,'ERROR: charge_dist should have an argument of'//&
                 ' either point or gaussian.'
            stop
          endif
          cycle
        else if( trim(cline).eq.'terms' ) then
          cmode = 'terms'
          backspace(ioprms)
          read(ioprms,*) ctmp, cterms
          if(  trim(cterms).ne.'full' .and. &
               trim(cterms).ne.'short' .and. &
               trim(cterms).ne.'screened' .and. &
               trim(cterms).ne.'long' ) then
            print *,'ERROR: terms should have an argument of '//&
                 'either full, short/screened or long.'
            stop
          endif
          cycle
        else if( trim(cline).eq.'interactions' ) then
          cmode = 'interactions'
          interact(1:msp,1:msp) = .false.
        endif
!.....Not a keyword, a certain mode should be already selected.
        if( trim(cmode).eq.'charges' ) then
          if( trim(cchgs).eq.'fixed' ) then
            read(ioprms,*) isp, chgi
            schg(isp) = chgi
          else if( trim(cchgs).eq.'fixed_bvs' ) then
            read(ioprms,*) isp,cname,vid,rad,npq
            if( isp.gt.nsp .and. iprint.gt.0 ) then
              print *,'WARNING: isp.gt.nsp !!!  isp = ',isp
            endif
            vid_bvs(isp) = vid
            rad_bvs(isp) = rad
            npq_bvs(isp) = npq
            if( iprint.ne.0 ) then
              write(6,'(a,i3,a5,2f7.3,i4)') ' isp,cname,vid,rad,npq =' &
                   ,isp,trim(cname),vid,rad,npq
            endif
          else if( trim(cchgs).eq.'variable' .or. trim(cchgs).eq.'qeq') then
            read(ioprms,*) isp, cname, dchi,djii,sgmt,de0,qlow,qup
            if( isp.gt.nsp .and. iprint.gt.0 ) then
              print *,'WARNING: isp.gt.nsp !!!  isp = ',isp
            endif
            vcg_chi(isp) = dchi
            vcg_jii(isp) = djii
            vcg_e0(isp) = de0
            vcg_sgm(isp) = sgmt
            qlower(isp) = qlow
            qupper(isp) = qup
            if( iprint.gt.0 ) then
              write(6,'(a,i3,a3,3f10.4,2f5.1)') '   isp,name,chi,Jii,sgm,e0,qlower,qupper = ', &
                   isp,trim(cname),dchi,djii,vcg_sgm(isp),de0,qlow,qup
            endif
          endif
!!$        else if( trim(cmode).eq.'charge_dist' ) then
!!$        else if( trim(cmode).eq.'terms' ) then
        else if( trim(cmode).eq.'interactions' ) then
          read(ioprms,*) isp,jsp
          interact(isp,jsp) = .true.
          interact(jsp,isp) = .true.
        endif
      enddo ! while(.true.)

!.....Corrections
      if( trim(cdist).eq.'gaussian' ) then
        if( trim(cterms).eq.'full' .or. trim(cterms).eq.'screened' .or. &
             trim(cterms).eq.'short' ) then
          cterms = 'long'
          print *,'WARNING: terms was corrected to long, because charge_dist is gaussian.'
        endif
      endif
      
10    close(ioprms)
    endif  ! myid.eq.0

!.....Broadcast data just read from the file
    call mpi_bcast(vcg_chi,nsp,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(vcg_jii,nsp,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(vcg_e0,nsp,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(vcg_sgm,nsp,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(qlower,nsp,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(qupper,nsp,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(vid_bvs,msp,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(rad_bvs,msp,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(npq_bvs,msp,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(rho_bvs,msp*msp,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(interact,msp*msp,mpi_logical,0,mpi_world,ierr)

    if( myid.eq.0 .and. iprint.ne.0 ) then
      write(6,'(a)') ' Finished reading '//trim(fname)
      write(6,*) ''
    endif
    
    return
  end subroutine read_paramsx
!=======================================================================
  subroutine force_screened_Coulomb(namax,natm,tag,ra,nnmax,aa,strs &
       ,chg,h,hi,tcom &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid,epi,epot,nismax,acon,lstrs,iprint,ifcoulomb&
       ,l1st)
    implicit none
    include "mpif.h"
    include "./params_unit.h"
!!$    include "params_BVS_Morse.h"
    integer,intent(in):: namax,natm,nnmax,nismax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),lspr(0:nnmax,namax),nex(3)
    integer,intent(in):: mpi_md_world,myid,ifcoulomb
    real(8),intent(in):: ra(3,namax),h(3,3),hi(3,3),rc &
         ,acon(nismax),tag(namax),sv(3,6)
    real(8),intent(inout):: chg(namax)
    real(8),intent(inout):: tcom
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: lstrs,l1st


    integer:: i,j,k,l,m,n,ierr,is,js,ixyz,jxyz,nconnect(4)
    real(8):: xi(3),xj(3),xij(3),rij(3),dij,diji,dedr &
         ,dxdi(3),dxdj(3),x,y,z,epotl,epott,at(3),tmp &
         ,qi,qj,radi,radj,rhoij,terfc,texp
    real(8),allocatable,save:: strsl(:,:,:)

    if( l1st ) then
      call set_charge_BVS(natm,nb,tag,chg,myid,mpi_md_world,iprint)
      if( allocated(strsl) ) deallocate(strsl)
      allocate(strsl(3,3,namax))
!!$      do i=1,natm
!!$        print *,'i,chg(i)=',i,chg(i)
!!$      enddo
    endif

    if( size(strsl).lt.3*3*namax ) then
      deallocate(strsl)
      allocate(strsl(3,3,namax))
    endif

    epotl= 0d0
    strsl(1:3,1:3,1:namax) = 0d0
!!$    write(6,'(a,30f7.3)') 'chgs =',chg(1:natm)

!.....Loop over resident atoms
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is= int(tag(i))
!!$      nconnect(i) = 0
      qi= chg(i)
      if( abs(qi).lt.qthd ) cycle
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
        if(j.le.i) cycle
        js= int(tag(j))
        if( .not.interact(is,js) ) cycle
        qj= chg(j)
        if( abs(qj).lt.qthd ) cycle
        xj(1:3)= ra(1:3,j)
        xij(1:3)= xj(1:3)-xi(1:3)
        rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij= sqrt(rij(1)**2 +rij(2)**2 +rij(3)**2)
        if( dij.gt.rc ) cycle
!!$        nconnect(i)= nconnect(i) +1
        diji= 1d0/dij
        dxdi(1:3)= -rij(1:3)*diji
        dxdj(1:3)=  rij(1:3)*diji
        rhoij = rho_bvs(is,js)
        terfc = erfc(dij/rhoij)
!.....potential
        tmp= 0.5d0 *acc *qi*qj*diji *terfc
        if( j.le.natm ) then
          epi(i)= epi(i) +tmp
          epi(j)= epi(j) +tmp
          epotl = epotl +tmp +tmp
        else
          epi(i)= epi(i) +tmp
          epotl = epotl +tmp
        endif
!.....force
        texp = exp(-(dij/rhoij)**2)
        dedr= -acc *qi*qj*diji *(1d0*diji*terfc +2d0/rhoij /sqrt(pi) *texp)
        aa(1:3,i)= aa(1:3,i) -dxdi(1:3)*dedr
        aa(1:3,j)= aa(1:3,j) -dxdj(1:3)*dedr
!.....stress
        if( lstrs ) then
          do ixyz=1,3
            do jxyz=1,3
              strsl(jxyz,ixyz,i)= strsl(jxyz,ixyz,i) &
                   -0.5d0 *dedr*rij(ixyz)*(-dxdi(jxyz))
              strsl(jxyz,ixyz,j)= strsl(jxyz,ixyz,j) &
                   -0.5d0 *dedr*rij(ixyz)*(-dxdi(jxyz))
            enddo
          enddo
        endif
      enddo
    enddo

    if( lstrs ) then
!!$      call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
!!$           ,nn,mpi_md_world,strsl,9)
      strs(1:3,1:3,1:natm)= strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)
    endif

!-----gather epot
    call mpi_allreduce(epotl,epott,1,MPI_REAL8 &
         ,MPI_SUM,mpi_md_world,ierr)
    epot= epot +epott
!!$    write(6,'(a,es15.7)') ' screened Coulomb epott = ',epott

  end subroutine force_screened_Coulomb
!=======================================================================
  subroutine force_Ewald(namax,natm,tag,ra,nnmax,aa,strs &
       ,chg,chi,h,hi,tcom &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid,epi,epot,nismax,acon,lstrs,iprint,ifcoulomb &
       ,l1st,lcell_updated,lvc)
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
    integer,intent(in):: mpi_md_world,myid,ifcoulomb
    real(8),intent(in):: ra(3,namax),h(3,3),hi(3,3),rc &
         ,acon(nismax),tag(namax),sv(3,6),chi(namax)
    real(8),intent(inout):: chg(namax)
    real(8),intent(inout):: tcom
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: lstrs,l1st,lcell_updated,lvc

    integer:: i,j,ik,is,js,k1,k2,k3,ierr,jj,ixyz,jxyz
    real(8):: elrl,esrl,epotl,epott,qi,qj,tmp,ftmp &
         ,bdotr,terfc,diji,dij,ss2i,sgmsq2,rc2,q2tot,q2loc,bb2,sqpi &
         ,e0,q2,sgmi
    real(8),save:: eself,eselfl
    real(8),allocatable,save:: ri(:),bk(:),bk1(:),bk2(:),bk3(:)&
         ,bb(:),dxdi(:),dxdj(:),rij(:),xij(:),xj(:),xi(:)
    real(8),allocatable,save:: strsl(:,:,:)
    real(8),external:: sprod,absv

    if( l1st ) then
      if( .not.allocated(ri) ) then
        allocate(ri(3),bk(3),bk1(3),bk2(3),bk3(3),bb(3),dxdi(3) &
             ,dxdj(3),rij(3),xij(3),xj(3),xi(3))
      endif
      if( allocated(strsl) ) deallocate(strsl)
      allocate(strsl(3,3,namax))
    endif

    if( size(strsl).lt.3*3*namax ) then
      deallocate(strsl)
      allocate(strsl(3,3,namax))
    endif

    strsl(1:3,1:3,1:namax) = 0d0

    call Ewald_self(namax,natm,tag,chg,chi,epi,eselfl,iprint,lvc)

    call Ewald_short(namax,natm,tag,ra,nnmax,aa,strsl,chg,h,hi &
         ,lspr,epi,esrl,iprint,ifcoulomb,lstrs,rc)

    call Ewald_long(namax,natm,tag,ra,nnmax,aa,strsl,chg,h,hi,&
         lspr,epi,elrl,iprint,ifcoulomb,mpi_md_world,lstrs,lcell_updated)

    if( lstrs ) then
!!$      call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
!!$           ,nn,mpi_md_world,strsl,9)
      strs(1:3,1:3,1:natm)= strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)
    endif

    if( l1st .and. myid.eq.0 ) then
      print *,'Ewald energy term by term:'
      print '(a,f12.4)','   Self term         = ',eselfl
      print '(a,f12.4)','   Short-range term  = ',esrl
      print '(a,f12.4)','   Long-range term   = ',elrl
    endif

    epotl = esrl +elrl +eselfl
!.....Gather epot
    call mpi_allreduce(epotl,epott,1,mpi_real8 &
         ,mpi_sum,mpi_md_world,ierr)
    epot= epot +epott

  end subroutine force_Ewald
!=======================================================================
  subroutine force_Ewald_long(namax,natm,tag,ra,nnmax,aa,strs &
       ,chg,chi,h,hi,tcom &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid,epi,epot,nismax,acon,lstrs,iprint,ifcoulomb &
       ,l1st,lcell_updated,lvc)
!
!  Long-range only Coulomb potential with Gaussian distribution charges.
!  Currently this does not work in the parallel computation.
!  TODO: Parallelize long_Coulomb potential
!
    implicit none
    include "mpif.h"
    include "./params_unit.h"

    integer,intent(in):: namax,natm,nnmax,nismax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),lspr(0:nnmax,namax),nex(3)
    integer,intent(in):: mpi_md_world,myid,ifcoulomb
    real(8),intent(in):: ra(3,namax),h(3,3),hi(3,3),rc &
         ,acon(nismax),tag(namax),sv(3,6),chi(namax)
    real(8),intent(inout):: chg(namax)
    real(8),intent(inout):: tcom
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: lstrs,l1st,lcell_updated,lvc

    integer:: i,ierr,is
    real(8):: elrl,eselfl,epotl,epott,prefac,q2,sgmi,e0,tmp
    real(8),save:: eself
    real(8),allocatable,save:: strsl(:,:,:)

    if( l1st ) then
      if( allocated(strsl) ) deallocate(strsl)
      allocate(strsl(3,3,namax))
    endif

    if( size(strsl).lt.3*3*namax ) then
      deallocate(strsl)
      allocate(strsl(3,3,namax))
    endif

    strsl(1:3,1:3,1:namax) = 0d0

    call Ewald_self(namax,natm,tag,chg,chi,epi,eselfl,iprint,lvc)
    
    call Ewald_long(namax,natm,tag,ra,nnmax,aa,strsl,chg,h,hi,&
         lspr,epi,elrl,iprint,ifcoulomb,mpi_md_world,lstrs,lcell_updated)

    if( lstrs ) then
!!$      call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
!!$           ,nn,mpi_md_world,strsl,9)
      strs(1:3,1:3,1:natm)= strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)
    endif

    if( l1st .and. myid.eq.0 ) then
      print *,'Ewald energy term by term:'
      print '(a,f12.4)','   Self term         = ',eselfl
      print '(a,f12.4)','   Long-range term   = ',elrl
    endif

    epotl = elrl +eselfl
!.....Gather epot
    call mpi_allreduce(epotl,epott,1,mpi_real8 &
         ,mpi_sum,mpi_md_world,ierr)
    epot= epot +epott

    return
  end subroutine force_Ewald_long
!=======================================================================
  subroutine Ewald_short(namax,natm,tag,ra,nnmax,aa,strsl,chg,h,hi &
       ,lspr,epi,esrl,iprint,ifcoulomb,lstrs,rc)
    implicit none
    integer,intent(in):: namax,natm,nnmax,iprint,ifcoulomb, &
         lspr(0:nnmax,namax)
    real(8),intent(in)::tag(namax),ra(3,namax),chg(namax), &
         h(3,3),hi(3,3),rc
    logical,intent(in):: lstrs
    real(8),intent(inout):: aa(3,namax),strsl(3,3,namax), &
         epi(namax),esrl

    integer:: i,j,jj,is,js,ixyz,jxyz
    real(8):: rc2,sgmsq2,ss2i,sqpi,qi,qj,dij,diji,tmp,ftmp,terfc
    real(8):: xi(3),xj(3),xij(3),rij(3),dxdi(3),dxdj(3)

!.....Compute direct sum
    rc2 = rc*rc
    sgmsq2 = sqrt(2d0)*sgm_ew
    ss2i = 1d0 /sgmsq2
    sqpi = 1d0 /sqrt(pi)
    esrl = 0d0
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is= int(tag(i))
      qi = chg(i)
!!$      sgmi = sgm(is)
!!$      if( abs(qi).lt.qthd ) cycle
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
!!$        terfc = erfc(gmmij,dij)
!.....potential
        tmp = 0.5d0 *acc *qi*qj*diji *terfc
        if( j.le.natm ) then
          epi(i)= epi(i) +tmp
          epi(j)= epi(j) +tmp
          esrl = esrl +tmp +tmp
        else
          epi(i)= epi(i) +tmp
          esrl = esrl +tmp
        endif
!.....force
        ftmp = -acc *qj*qi*diji *( diji *terfc &
             +2d0 *sqpi *ss2i *exp(-(dij*ss2i)**2) )
        aa(1:3,i)= aa(1:3,i) -dxdi(1:3)*ftmp
        aa(1:3,j)= aa(1:3,j) -dxdj(1:3)*ftmp
!.....stress
        if( lstrs ) then
          do ixyz=1,3
            do jxyz=1,3
              strsl(jxyz,ixyz,i)= strsl(jxyz,ixyz,i) &
                   -0.5d0 *ftmp*rij(ixyz)*(-dxdi(jxyz))
              strsl(jxyz,ixyz,j)= strsl(jxyz,ixyz,j) &
                   -0.5d0 *ftmp*rij(ixyz)*(-dxdi(jxyz))
            enddo
          enddo
        endif
      enddo
    enddo
!!$    write(6,*) ' esrl = ',esrl

  end subroutine Ewald_short
!=======================================================================
  subroutine Ewald_long(namax,natm,tag,ra,nnmax,aa,strsl,chg,h,hi,&
       lspr,epi,elrl,iprint,ifcoulomb,mpi_md_world,lstrs,lcell_updated)
!
!  Long-range term of Ewald sum.
!  Not parallelized and not suitable for large system as it is too slow.
!
    implicit none
    include 'mpif.h'
    integer,intent(in):: namax,natm,nnmax,iprint,ifcoulomb, &
         mpi_md_world,lspr(0:nnmax,namax)
    real(8),intent(in):: tag(namax),ra(3,namax),chg(namax),&
         h(3,3),hi(3,3)
    logical,intent(in):: lstrs,lcell_updated
    real(8),intent(inout):: aa(3,namax),epi(namax),&
         elrl,strsl(3,3,namax)

    integer:: i,is,ik,k1,k2,k3,ierr,ixyz,jxyz
    real(8):: qi,xi(3),ri(3),bk1(3),bk2(3),bk3(3),bb(3),bdotr, &
         tmp,cs,sn,texp,bb2,bk
    real(8),external:: sprod,absv
    real(8):: emat(3,3)

!.....Compute reciprocal vectors
    if( lcell_updated ) call get_recip_vectors(h)
!.....Compute structure factor of the local processor
    call calc_qcos_qsin(namax,natm,tag,ra,chg,h,iprint,mpi_md_world)

!.....Compute long-range contribution to potential energy
    elrl = 0d0
!.....Long-range contribution to forces and stress
    emat(1:3,1) = (/ 1d0, 0d0, 0d0 /)
    emat(1:3,2) = (/ 0d0, 1d0, 0d0 /)
    emat(1:3,3) = (/ 0d0, 0d0, 1d0 /)
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is= int(tag(i))
      qi = chg(i)
!!$      if( abs(qi).lt.qthd ) cycle
      ri(1:3) = h(1:3,1)*xi(1) +h(1:3,2)*xi(2) +h(1:3,3)*xi(3)
      ik = 0
      do k1= -kmax1,kmax1
        bk1(1:3) = k1 *b1(1:3)
        do k2= -kmax2,kmax2
          bk2(1:3) = k2 *b2(1:3)
          do k3= -kmax3,kmax3
            if( .not. lkuse(k3,k2,k1) ) cycle
            ik= ik +1
            bk3(1:3) = k3 *b3(1:3)
            bb(1:3) = bk1(1:3) +bk2(1:3) +bk3(1:3)
            bdotr = sprod(3,bb,ri)
            cs = cos(bdotr)
            sn = sin(bdotr)
!.....Potential energy per atom
            tmp = 0.5d0 *acc /vol *qi *pflr(ik,is) &
                 *( cs*qcos(ik) +sn*qsin(ik) )
            epi(i) = epi(i) +tmp
            elrl = elrl +tmp
!.....Forces
            aa(1:3,i)= aa(1:3,i) -acc/vol *qi*bb(1:3) *pflr(ik,is) &
                 *( -sn*qcos(ik) +cs*qsin(ik) )
!.....Stress
            if( lstrs ) then
              bk = absv(3,bb)
!!$              tmp = 0.5d0 *acc/vol**2 &
!!$                   *pflr(ik,is) *(qcos(ik)**2 +qsin(ik)**2)
              do ixyz=1,3
                do jxyz=1,3
                  strsl(ixyz,jxyz,i) = strsl(ixyz,jxyz,i) +tmp &
                       *( bb(ixyz)*bb(jxyz)/bk*(bk *sgm(is)**2 +2d0/bk) &
                       -emat(ixyz,jxyz))
                enddo
              enddo

            endif
          enddo  ! k3
        enddo  ! k2
      enddo  ! k1
    enddo  ! ia

  end subroutine Ewald_long
!=======================================================================
  subroutine Ewald_self(namax,natm,tag,chg,chi,&
       epi,eselfl,iprint,lvc)
    implicit none
    integer,intent(in):: namax,natm,iprint
    logical,intent(in):: lvc
    real(8),intent(in):: tag(namax),chg(namax),chi(namax)
    real(8),intent(inout):: epi(namax),eselfl

    integer:: i,is
    real(8):: sgmi,qi,q2,e0,tmp

!.....Compute self term. If charge per atom is fixed, it is constant.
    if( lvc ) then  ! variable charge
      eselfl = 0d0
      do i=1,natm
        is = int(tag(i))
        sgmi = sgm(is)
        qi = chg(i)
        q2 = qi*qi
        e0 = vcg_e0(is)
        tmp = vcg_jii(is) -acc *sqrt(2d0/pi) /sgmi
!!$        tmp = vcg_jii(is)
        eselfl = eselfl +e0 +chi(i)*qi +0.5d0*tmp*q2
        epi(i) = epi(i) +e0 +chi(i)*qi +0.5d0*tmp*q2
!!$        write(6,'(a,2i4,3es12.4)') '  i,is,jii,chi*chg,tmp*q2/2=',&
!!$             i,is,vcg_jii(is), &
!!$             chi(i)*qi,0.5d0*tmp*q2
      enddo
    else  ! fixed charge
      eselfl = 0d0
      do i=1,natm
        is = int(tag(i))
        sgmi = sgm(is)
        q2 = chg(i)*chg(i)
        eselfl = eselfl -q2 /sgmi *acc/sqrt(2d0*pi)
        epi(i) = epi(i) -q2 /sgmi *acc/sqrt(2d0*pi)
      enddo
    endif

  end subroutine Ewald_self
!=======================================================================
  subroutine qforce_short(namax,natm,tag,ra,nnmax,chg,h &
       ,lspr,iprint,ifcoulomb,rc,fq,esr)
!
!  Compute q-force of short-range term of Ewald sum.
!
    implicit none
    integer,intent(in):: namax,natm,nnmax,iprint,ifcoulomb, &
         lspr(0:nnmax,namax)
    real(8),intent(in)::tag(namax),ra(3,namax),chg(namax), &
         h(3,3),rc
    real(8),intent(inout):: fq(namax),esr

    integer:: i,j,jj,is,js,ixyz,jxyz
    real(8):: rc2,sgmsq2,ss2i,sqpi,qi,qj,dij,diji,tmp,ftmp,terfc&
         ,sgmi,sgmj,gmmij
    real(8):: xi(3),xj(3),xij(3),rij(3),dxdi(3),dxdj(3)

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
        tmp = 0.5d0 *acc *diji *terfc
        if( j.le.natm ) then
          esr = esr +2d0*tmp*qi*qj
        else
          esr = esr +tmp*qi*qj
        endif
!.....Force on charge
        fq(i) = fq(i) -tmp*qj
        fq(j) = fq(j) -tmp*qi
      enddo
    enddo

    return
  end subroutine qforce_short
!=======================================================================
  subroutine qforce_long(namax,natm,tag,ra,chg,h, &
       tcom,mpi_md_world,myid,iprint,ifcoulomb,fq,elr)
!
!  Derivative of Ewald long-range term w.r.t. charges
!
    implicit none
    include 'mpif.h'
    integer,intent(in):: namax,natm,mpi_md_world,myid,iprint,ifcoulomb
    real(8),intent(in):: tag(namax),ra(3,namax),chg(namax),h(3,3)
    real(8),intent(inout):: tcom,fq(namax),elr

    integer:: i,ik,k1,k2,k3,is,ierr
    real(8):: prefac,bk1(3),bk2(3),bk3(3),bb(3),bb2,xi(3),ri(3), &
         sgmi,sgmi2,qi,bdotr,texp,cs,sn,elrl,eselfl,q2,epotl,tmp
    real(8),external:: sprod

!.....Compute reciprocal vectors
    call get_recip_vectors(h)
    prefac = 1d0 /(2d0*vol*eps0)
!!$    print *,'prefac=',prefac
!.....Compute structure factor
    call calc_qcos_qsin(namax,natm,tag,ra,chg,h,iprint,mpi_md_world)

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
          bb2 = sprod(3,bb,bb)
          do i=1,natm
            xi(1:3)= ra(1:3,i)
            is= int(tag(i))
            sgmi = sgm(is)
            sgmi2= sgmi*sgmi
            qi = chg(i)
            ri(1:3) = h(1:3,1)*xi(1) +h(1:3,2)*xi(2) +h(1:3,3)*xi(3)
            bdotr = sprod(3,bb,ri)
            texp = exp(-bb2*sgmi2/2)
            cs = cos(bdotr)
            sn = sin(bdotr)
!.....Potential energy per atom
            tmp = prefac/bb2*qi*texp &
                 *(cs*qcos(ik) +sn*qsin(ik))
            elr = elr +tmp
!.....Force on charge
            fq(i)= fq(i) -2d0 *prefac/bb2 *texp &
                 *(cs*qcos(ik) +sn*qsin(ik))
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
    enddo  ! ia

!!$    epotl = eselfl +elrl
!!$    call mpi_allreduce(epotl,epot,1,mpi_real8, &
!!$         mpi_sum,mpi_md_world,ierr)
    return
  end subroutine qforce_long
!=======================================================================
  subroutine qforce_self(namax,natm,tag,chg,chi,fq,eself)
!
!  Derivative of self term w.r.t. charges
!
    integer,intent(in):: namax,natm
    real(8),intent(in):: tag(namax),chi(namax),chg(namax)
    real(8),intent(inout):: fq(namax),eself

    integer:: i,is
    real(8):: qi,q2,sgmi,tmp

    eself = 0d0
    do i=1,natm
      is = int(tag(i))
      qi = chg(i)
      q2 = qi*qi
      sgmi = sgm(is)
      tmp = vcg_jii(is) -acc*sqrt(2d0/pi) /sgmi
!!$      tmp = vcg_jii(is)
      eself = eself +vcg_e0(is) +chi(i)*qi +0.5d0*tmp*q2
      fq(i) = fq(i) -(chi(i) +tmp*qi)
    enddo

  end subroutine qforce_self
!=======================================================================
  subroutine set_charge_BVS(natm,nb,tag,chg,myid,mpi_md_world,iprint)
!
! Reset actual charges of atoms using effective charge information
! in order to keep charge neutrality in the system.
! It is not necessary when performing screened Coulomb, but required
! for non-screened Coulomb such as Ewald sum method. So it may be
! natural to set neutral charges even when using screened Coulomb...
!
! This would be called only once at the beginning.
!
    include "mpif.h"
    integer,intent(in):: natm,nb,myid,mpi_md_world,iprint
    real(8),intent(in):: tag(natm+nb)
    real(8),intent(out):: chg(natm+nb)

    integer,allocatable:: nbvsl(:),nbvs(:)
    real(8),allocatable:: vc_bvs(:)
    integer:: i,is,ierr
    real(8):: sum_anion,sum_cation

!!$!.....BEGIN DEBUGGING
!!$    do i=1,natm+nb
!!$      is = int(tag(i))
!!$      if( is.eq.1 ) then
!!$        chg(i) = -vid_bvs(is)
!!$      else
!!$        chg(i) = vid_bvs(is)
!!$      endif
!!$      print *,' i,is,chg(i)=',i,is,chg(i)
!!$    enddo
!!$    return
!!$!.....END DEBUGGING

    allocate(nbvsl(msp),nbvs(msp),vc_bvs(msp))
    nbvsl(1:nsp) = 0
    nbvs(1:nsp) = 0
    do i=1,natm
      is = int(tag(i))
      nbvsl(is) = nbvsl(is) +1
    enddo
    call mpi_allreduce(nbvsl,nbvs,nsp,mpi_integer &
         ,mpi_sum,mpi_md_world,ierr)

    sum_anion = vid_bvs(1) *nbvs(1) /sqrt(dble(npq_bvs(1)))

    sum_cation = 0d0
    do is=2,nsp
      sum_cation = sum_cation +vid_bvs(is)*nbvs(is)/sqrt(dble(npq_bvs(is)))
    enddo

!.....Compute valence charges of species
    vc_bvs(1) = vid_bvs(1)/sqrt(dble(npq_bvs(1))) &
         *sqrt(sum_cation/sum_anion)
    do is=2,nsp
      vc_bvs(is)= vid_bvs(is)/sqrt(dble(npq_bvs(is))) &
           *sqrt(sum_anion/sum_cation)
    enddo

!!$!.....Overwrite charges for debugging...
!!$    vc_bvs(1:4) = (/ 1.27443, 0.78466, 1.10968, 3.20338/)

    if( myid.eq.0 .and. iprint.gt.0 ) then
      print *,''
      print *,'Charges fixed from ideal valences and composition:'
      do is=1,nsp
        write(6,'(a,i3,2f7.3)') '   isp, V_ideal, V_actual = '&
             ,is,vid_bvs(is),vc_bvs(is)
      enddo
    endif

    do i=1,natm+nb
      is= int(tag(i))
!.....Negative charge for anion and positive for cation
      if( is .eq. 1 ) then
        chg(i)= -vc_bvs(is)
      else
        chg(i)= vc_bvs(is)
      endif
!!$      if( i.le.natm ) then
!!$        write(6,'(a,2i5,f10.4)') ' i,is,chg=',i,is,chg(i)
!!$      endif
    enddo

    deallocate(nbvsl,nbvs,vc_bvs)

  end subroutine set_charge_BVS
!=======================================================================
  subroutine get_recip_vectors(h)
    implicit none
    real(8),intent(in):: h(3,3)

    real(8):: a1(3),a2(3),a3(3),a23(3),a12(3),a31(3),pi2
    real(8),external:: sprod

    a1(1:3) = h(1:3,1)
    a2(1:3) = h(1:3,2)
    a3(1:3) = h(1:3,3)
    call vprod(a2,a3,a23)
    call vprod(a3,a1,a31)
    call vprod(a1,a2,a12)
    vol = abs(sprod(3,a1,a23))
    pi2 = 2d0 *pi
    b1(1:3) = pi2 /vol *a23(1:3)
    b2(1:3) = pi2 /vol *a31(1:3)
    b3(1:3) = pi2 /vol *a12(1:3)
    return
  end subroutine get_recip_vectors
!=======================================================================
  subroutine calc_qcos_qsin(namax,natm,tag,ra,chg,h,iprint&
       ,mpi_md_world)
!
!  Compute qcos and qsin needed for calculation of Ewald long-range term.
!
    implicit none
    include 'mpif.h'
    integer,intent(in):: namax,natm,iprint,mpi_md_world
    real(8),intent(in):: tag(namax),ra(3,namax),chg(namax) &
         ,h(3,3)

    integer:: ik,k1,k2,k3,is,i,ierr
    real(8):: bk1(3),bk2(3),bk3(3),bb(3),xi(3),qi&
         ,ri(3),bdotr
    real(8),external:: sprod

!.....Compute structure factor of the local processor
    qcosl(1:nk) = 0d0
    qsinl(1:nk) = 0d0
    ik = 0
    do k1= -kmax1,kmax1
      bk1(1:3) = k1 *b1(1:3)
      do k2= -kmax2,kmax2
        bk2(1:3) = k2 *b2(1:3)
        do k3= -kmax3,kmax3
          if( .not. lkuse(k3,k2,k1) ) cycle
          ik= ik +1
          bk3(1:3) = k3 *b3(1:3)
          bb(1:3) = bk1(1:3) +bk2(1:3) +bk3(1:3)
          do i=1,natm
            xi(1:3)= ra(1:3,i)
            qi = chg(i)
            ri(1:3) = h(1:3,1)*xi(1) +h(1:3,2)*xi(2) +h(1:3,3)*xi(3)
            bdotr = sprod(3,bb,ri)
            qcosl(ik) = qcosl(ik) +qi*cos(bdotr)
            qsinl(ik) = qsinl(ik) +qi*sin(bdotr)
          enddo
        enddo
      enddo
    enddo  ! ia
!.....Allreduce qcos and qsin, which would be stupid and time consuming
    qcos(1:nk) = 0d0
    qsin(1:nk) = 0d0
    call mpi_allreduce(qcosl,qcos,nk,mpi_real8 &
         ,mpi_sum,mpi_md_world,ierr)
    call mpi_allreduce(qsinl,qsin,nk,mpi_real8 &
         ,mpi_sum,mpi_md_world,ierr)
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
    real(8),external:: absv

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
          bb2 = absv(3,bk)
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

    if( allocated(lkuse) ) deallocate(lkuse)
    allocate(lkuse(-kmax3:kmax3,-kmax2:kmax2,-kmax1:kmax1))

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
  subroutine cgopt_charge(natm,h,ra,tag,chg,chi)
!
!  Charge optimization for variable-charge potential
!  by means of conjugate gradient.
!  In this case, only Ewald_long part is used because the Gaussian
!  distribution is assumed for atomic charges.
!
    integer,intent(in):: natm
    real(8),intent(in):: h(3,3),ra(3,natm),tag(natm),chi(natm)
    real(8),intent(inout):: chg(natm)

    integer:: ia,ja,k1,k2,k3,is,js,ik,ierr
    real(8):: xi(3),xj(3),ri(3),rj(3),qi,qj,sgmi,sgmi2,sgmj,sgmj2, &
         bk1(3),bk2(3),bk3(3),bb(3),bdotri,bdotrj,bb2,djii,a0,sqpi,r2,&
         prefac
    real(8),allocatable,save:: amat(:,:),qvec(:),xvec(:)
    logical,save:: l1st = .true.
    real(8),external:: sprod

!.....TODO: not work in parallel
    if( l1st ) then
      if( allocated(amat) ) deallocate(amat,qvec,xvec)
      allocate(amat(natm+1,natm+1),qvec(natm+1),xvec(natm+1))
    endif

    amat(1:natm+1,1:natm+1) = 0d0
    qvec(natm+1) = 0d0

!.....Make A-matrix for linear equation, A*Q=X
!     Assuming all the atoms are in the root node (0),
!     and does not work in parallel mode.
!     TODO: And also the code is very primitive and could be faster.
    a0 = 1d0 /(2d0 *vol *eps0)
    do ia=1,natm
      xi(1:3) = ra(1:3,ia)
      ri(1:3) = h(1:3,1)*xi(1) +h(1:3,2)*xi(2) +h(1:3,3)*xi(3)
      is = int(tag(ia))
      sgmi = sgm(is)
      sgmi2= sgmi*sgmi
      do ja=ia,natm
        xj(1:3) = ra(1:3,ja)
        rj(1:3) = h(1:3,1)*xj(1) +h(1:3,2)*xj(2) +h(1:3,3)*xj(3)
        js = int(tag(ja))
        sgmj = sgm(js)
        sgmj2= sgmj*sgmj
        ik = 0
        do k1= -kmax1,kmax1
          bk1(1:3)= k1 *b1(1:3)
          do k2= -kmax2,kmax2
            bk2(1:3)= k2 *b2(1:3)
            do k3= -kmax3,kmax3
              if( .not. lkuse(k3,k2,k1) ) cycle
              ik = ik + 1
              bk3(1:3)= k3 *b3(1:3)
              bb(1:3) = bk1(1:3) +bk2(1:3) +bk3(1:3)
              bdotri = sprod(3,bb,ri)
              bdotrj = sprod(3,bb,rj)
              bb2 = sprod(3,bb,bb)

              amat(ia,ja) = amat(ia,ja) +2d0/bb2 &
                   *exp(-bb2/2*(sgmi2+sgmj2)) &
                   *( cos(bdotri)*cos(bdotrj) + sin(bdotri)*sin(bdotrj) )
            enddo  ! k3
          enddo  ! k2
        enddo  ! k1
        amat(ia,ja) = amat(ia,ja) *a0
      enddo  ! ja
    enddo  ! ia
    amat(1:natm,natm+1) = 1d0
!.....Self term
    prefac = 1d0/(4d0*pi*eps0)
    do ia=1,natm
      is = int(tag(ia))
      djii = vcg_jii(is)
      sgmi = sgm(is)
      amat(ia,ia) = amat(ia,ia) +djii !-2d0 *prefac /sqrt(pi) /sgmi
    enddo
!.....Symmetrize A-matrix
    do ia=1,natm
      do ja=ia+1,natm+1
        amat(ja,ia) = amat(ia,ja)
      enddo
    enddo

!.....Make Q vector
    do ia=1,natm
      qvec(ia) = chg(ia)
    enddo

!.....Make X vector
    do ia=1,natm
      is = int(tag(ia))
      xvec(ia) = -vcg_chi(is)
    enddo
    xvec(natm+1) = 0d0  ! charge neutrality, qtot = 0.0

!!$    print *,'amat:'
!!$    do ia=1,natm+1
!!$      write(6,'(10es10.2)') (amat(ia,ja),ja=1,natm+1)
!!$    enddo
!!$    print *,'qvec:'
!!$    write(6,'(a,10f8.4)') 'qvec=',(qvec(ja),ja=1,natm+1)
!!$    print *,'xvec:'
!!$    write(6,'(10es10.2)') (xvec(ja),ja=1,natm+1)

!.....Perform CG optimization to equilibrate the system and get charges
    call cg(natm+1,natm+2,1d-5,amat,xvec,qvec,ierr)
!!$    write(6,'(a,i6)') ' CG opt converged at ',ierr
!!$    print *,'qvec after CG opt:'
!!$    write(6,'(a,10f10.6)') 'q=',(qvec(ja),ja=1,natm)
!!$    print *,'A*Q:'
!!$    xvec = matmul(amat,qvec)
!!$    write(6,'(10es10.2)') (xvec(ja),ja=1,natm+1)

!.....Restore charge to chg
    do ia=1,natm
      chg(ia) = qvec(ia)
    enddo

    l1st = .false.
    return
  end subroutine cgopt_charge
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
end module Coulomb
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
