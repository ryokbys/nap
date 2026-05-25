module ShellModel
!-----------------------------------------------------------------------
!                     Last-modified: <2026-05-24 KOBAYASHI Ryo>
!-----------------------------------------------------------------------
!  Dick-Overhauser shell model: core-shell spring interaction.
!  Shell particles are treated as regular atoms with dedicated species
!  names (e.g. "O_s" for the oxygen shell).
!  This module computes only the spring potential between each
!  shell species and its paired core species:
!
!      V_spring = k2_s/2 * |r|^2 + k4_s/24 * |r|^4
!
!  Short-range (Buckingham) and long-range (Coulomb) interactions are
!  handled by their respective modules, using the shell species entries
!  in the normal specorder array.
!
!  Parameter file: in.params.ShellModel
!  Format (one entry per polarisable species pair):
!    # core_species  shell_species  k2_s[eV/Ang^2]  k4_s[eV/Ang^4]
!      O   O_s   74.92   0.0
!  (k4_s is optional; defaults to 0 if omitted)
!-----------------------------------------------------------------------
  use pmdmpi
  use mod_precision
  use pmdvars,only: nspmax,ntot,tag_itot,am
  use util,only: csp2isp
  use memory,only: accum_mem
  implicit none
  include "./const.h"
  save

  character(len=128):: paramsdir = '.'
  character(len=128),parameter:: cprmfname = 'in.params.ShellModel'
  integer,parameter:: ioprms = 20

  real(rp):: sm_k2s(nspmax)         ! 2nd-order spring constant [eV/Ang^2]
  real(rp):: sm_k4s(nspmax)         ! 4th-order spring constant [eV/Ang^4]
  integer :: sm_core_of(nspmax)     ! sm_core_of(ishell_sp) = icore_sp
  integer :: sm_shell_of(nspmax)    ! sm_shell_of(icore_sp) = ishell_sp
  logical :: is_shell_sp(nspmax)    ! .true. if species is a shell

  real(rp),parameter:: rc_coul_ex  = 1.5_rp
  real(rp),parameter:: rc_coul_ex2 = rc_coul_ex*rc_coul_ex

  logical:: lprmset_ShellModel = .false.

!.....Extended Lagrangian variables
  logical :: use_xl = .false.
  real(rp),allocatable,save:: xl_theta(:,:)   ! (3,namax) auxiliary shell positions (fractional)
  real(rp),allocatable,save:: xl_thdot(:,:)   ! (3,namax) auxiliary shell velocities
  real(rp),allocatable,save:: xl_thacc(:,:)   ! (3,namax) auxiliary shell accelerations
  real(rp):: xl_omega2 = 0.0_rp               ! = K/dt^2
!!$  real(rp),parameter:: xl_K = 0.05_rp         ! K = omega^2*dt^2; must satisfy K < 4*k2s/H_total
  real(rp),parameter:: xl_K = 2.0_rp         ! K = omega^2*dt^2; must satisfy K < 4*k2s/H_total
!  kappa = xl_step_factor/k2s; must satisfy xl_step_factor < k2s/H_total to avoid overshoot.
!  For strongly-ionic systems (BaTiO3) H_total >> k2s; use xl_step_factor << 1.
  real(rp),parameter:: xl_step_factor = 0.1_rp
  real(rp):: xl_kappa_sp(nspmax)              ! per-species gradient step: xl_step_factor/k2s_i

contains
!=======================================================================
  subroutine init_ShellModel()
    implicit none

    sm_k2s(:)      = 0.0_rp
    sm_k4s(:)      = 0.0_rp
    sm_core_of(:)  = 0
    sm_shell_of(:) = 0
    is_shell_sp(:) = .false.
    return
  end subroutine init_ShellModel
!=======================================================================
  subroutine force_ShellModel(namax,natm,tag_isp,ra,nnmax,aa,strs,h,hi &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid,epi,epot,nismax,lstrs,iprint,l1st)
!-----------------------------------------------------------------------
!  Compute core-shell spring forces and energies.
!  Uses the full-neighbour-list convention: force on centre atom only.
!  Each pair is processed twice (once with shell as centre, once with
!  core as centre), so the energy prefactor is 0.5 * V(d).
!
!  V(d) = k2s/2 * d^2 + k4s/24 * d^4
!  dV/dd = k2s * d + k4s/6 * d^3
!-----------------------------------------------------------------------
    implicit none
    include "./params_unit.h"
    integer,intent(in):: namax,natm,nnmax,nismax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),lspr(0:nnmax,namax),nex(3)
    integer,intent(in):: mpi_md_world,myid
    integer,intent(in):: tag_isp(namax)
    real(rp),intent(in):: ra(3,namax),h(3,3),hi(3,3),rc,sv(3,6)
    real(rp),intent(inout):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: l1st
    logical:: lstrs

    integer:: i,j,k,ierr,is,js,ixyz,jxyz,icore_sp,ishell_sp
    real(rp):: xi(3),xij(3),rij(3),dij,diji,dvdr,dij2,dij4 &
         ,dxdi(3),epotl,epott,k2s,k4s,tmp
    real(rp),allocatable,save:: strsl(:,:,:)

    if( l1st ) then
      if( allocated(strsl) ) then
        call accum_mem('force_ShellModel',-rp*size(strsl))
        deallocate(strsl)
      endif
      allocate(strsl(3,3,namax))
      call accum_mem('force_ShellModel',rp*size(strsl))
    endif

    if( size(strsl).lt.3*3*namax ) then
      call accum_mem('force_ShellModel',-rp*size(strsl))
      deallocate(strsl)
      allocate(strsl(3,3,namax))
      call accum_mem('force_ShellModel',rp*size(strsl))
    endif

    epotl = 0.0_rp
    strsl(1:3,1:3,1:namax) = 0.0_rp

!-----loop over resident atoms
    do i=1,natm
      xi(1:3) = ra(1:3,i)
      is = tag_isp(i)

!-----case 1: i is a shell atom -> find its core j (must be within rc_coul_ex)
      if( is_shell_sp(is) ) then
        icore_sp = sm_core_of(is)
        k2s = sm_k2s(is)
        k4s = sm_k4s(is)
        do k=1,lspr(0,i)
          j = lspr(k,i)
          if( j.eq.0 ) exit
          js = tag_isp(j)
          if( js.ne.icore_sp ) cycle
          xij(1:3) = ra(1:3,j) -xi(1:3)
          rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
          dij2 = rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
          if( dij2.gt.rc_coul_ex2 ) cycle  ! skip distant same-species atoms
          if( dij2.lt.1.0e-20_rp ) exit    ! coincident: spring force = 0
          dij4 = dij2*dij2
          dij  = sqrt(dij2)
          diji = 1.0_rp /dij
          dxdi(1:3) = -rij(1:3)*diji
          dvdr = k2s*dij + k4s*dij2*dij/6.0_rp
!---------force on shell i (toward core j)
          aa(1:3,i) = aa(1:3,i) -dxdi(1:3)*dvdr
!---------energy: 0.5*V(d) to avoid double-counting with case 2
          tmp = 0.5_rp *(0.5_rp*k2s*dij2 + k4s*dij4/24.0_rp)
          epi(i) = epi(i) +tmp
          epotl = epotl +tmp
!---------stress
          if( lstrs ) then
            do ixyz=1,3
              do jxyz=1,3
                strsl(jxyz,ixyz,i) = strsl(jxyz,ixyz,i) &
                     -0.5_rp*dvdr*rij(ixyz)*(-dxdi(jxyz))
              enddo
            enddo
          endif
          exit   ! each shell has exactly one core partner

!-----case 2: i is a core atom that has a shell -> find its shell j
        enddo
      else if( sm_shell_of(is).gt.0 ) then
        ishell_sp = sm_shell_of(is)
        k2s = sm_k2s(ishell_sp)
        k4s = sm_k4s(ishell_sp)
        do k=1,lspr(0,i)
          j = lspr(k,i)
          if( j.eq.0 ) exit
          js = tag_isp(j)
          if( js.ne.ishell_sp ) cycle
          xij(1:3) = ra(1:3,j) -xi(1:3)
          rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
          dij2 = rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
          if( dij2.gt.rc_coul_ex2 ) cycle  ! skip distant same-species atoms
          if( dij2.lt.1.0e-20_rp ) exit    ! coincident: spring force = 0
          dij4 = dij2*dij2
          dij  = sqrt(dij2)
          diji = 1.0_rp /dij
          dxdi(1:3) = -rij(1:3)*diji
          dvdr = k2s*dij + k4s*dij2*dij/6.0_rp
!---------force on core i (toward shell j)
          aa(1:3,i) = aa(1:3,i) -dxdi(1:3)*dvdr
!---------energy: 0.5*V(d) to avoid double-counting with case 1
          tmp = 0.5_rp *(0.5_rp*k2s*dij2 + k4s*dij4/24.0_rp)
          epi(i) = epi(i) +tmp
          epotl = epotl +tmp
!---------stress
          if( lstrs ) then
            do ixyz=1,3
              do jxyz=1,3
                strsl(jxyz,ixyz,i) = strsl(jxyz,ixyz,i) &
                     -0.5_rp*dvdr*rij(ixyz)*(-dxdi(jxyz))
              enddo
            enddo
          endif
          exit   ! each core has exactly one shell partner
        enddo
      endif
    enddo

    if( lstrs ) then
      strs(1:3,1:3,1:natm) = strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)
    endif

!-----gather epot
    epott = 0.0_rp
    call mpi_allreduce(epotl,epott,1,mpi_real_rp,mpi_sum,mpi_md_world,ierr)
    epot = epot +epott
    if( myid.eq.0 .and. iprint.ge.ipl_info ) &
         write(6,'(a,es15.7)') ' epot ShellModel = ',epott
    return
  end subroutine force_ShellModel
!=======================================================================
  subroutine xl_init(namax,natm,ra,va,tag_isp,dt)
!-----------------------------------------------------------------------
!  Initialise extended-Lagrangian auxiliary arrays for shell atoms.
!  Called once before the main VV loop when use_xl = .true.
!  Sets  θ = ra_shell,  θ̇ = 0,  θ̈ = 0.
!  Also computes ω² = 2/dt²  and  κ_i = 1/k2s_i per shell species.
!  Zeroes shell velocities va: in XL mode shells are propagated by the
!  auxiliary theta dynamics, not by Newton's law, so va_shell must be
!  zero to avoid a spurious constant contribution to ekin and temperature.
!
!  xl_theta/xl_thdot/xl_thacc are indexed by tag_itot(i) (global atom ID)
!  so they remain valid after any array rearrangement by bamove/bacopy.
!-----------------------------------------------------------------------
    implicit none
    integer,intent(in):: namax,natm,tag_isp(namax)
    real(rp),intent(in):: ra(3,namax),dt
    real(rp),intent(inout):: va(3,namax)

    integer:: i,is

    if( allocated(xl_theta) ) then
      call accum_mem('xl_init',-rp*size(xl_theta))
      deallocate(xl_theta,xl_thdot,xl_thacc)
    endif
    allocate(xl_theta(3,ntot),xl_thdot(3,ntot),xl_thacc(3,ntot))
    call accum_mem('xl_init',rp*size(xl_theta))

    xl_theta(:,:) = 0.0_rp
    xl_thdot(:,:) = 0.0_rp
    xl_thacc(:,:) = 0.0_rp
    do i=1,natm
      is = tag_isp(i)
      if( is.gt.0 .and. is_shell_sp(is) ) then
        xl_theta(1:3,tag_itot(i)) = ra(1:3,i)
!       zero shell velocity: position driven by theta dynamics, not Newton
        va(1:3,i) = 0.0_rp
      endif
    enddo

    xl_omega2 = xl_K / (dt*dt)

!   Per-species gradient step size: κ_i = xl_step_factor/k2s_i
!   xl_step_factor << 1 required for strongly-ionic systems where H_total >> k2s
    xl_kappa_sp(:) = 0.0_rp
    do is=1,nspmax
      if( is_shell_sp(is) .and. sm_k2s(is).gt.0.0_rp ) then
        xl_kappa_sp(is) = xl_step_factor / sm_k2s(is)
        write(6,'(a,i3,a,es12.4)') ' XL_INIT: shell species',is,' kappa =',xl_kappa_sp(is)
      endif
    enddo
    write(6,'(a,es12.4)') ' XL_INIT: xl_omega2 =',xl_omega2
    return
  end subroutine xl_init
!=======================================================================
  subroutine xl_predict(natm,tag_isp,ra,dt)
!-----------------------------------------------------------------------
!  Extended-Lagrangian predictor step for shell atoms.
!  Updates θ with velocity-Verlet half-step then full-step:
!    θ̇_half = θ̇ + 0.5*dt*θ̈
!    θ(t+dt) = θ(t) + dt*θ̇_half
!  Then sets ra_shell = θ(t+dt) for the upcoming force evaluation.
!-----------------------------------------------------------------------
    implicit none
    integer,intent(in):: natm,tag_isp(natm)
    real(rp),intent(inout):: ra(3,natm)
    real(rp),intent(in):: dt

    integer:: i,itot

    do i=1,natm
      if( .not.is_shell_sp(tag_isp(i)) ) cycle
      itot = tag_itot(i)
      xl_thdot(1:3,itot) = xl_thdot(1:3,itot) + 0.5_rp*dt*xl_thacc(1:3,itot)
      xl_theta(1:3,itot) = xl_theta(1:3,itot) + dt*xl_thdot(1:3,itot)
      ra(1:3,i) = xl_theta(1:3,itot)
    enddo
    return
  end subroutine xl_predict
!=======================================================================
  subroutine xl_gradient_step(natm,tag_isp,ra,aa,hi,h,dt,eaux)
!-----------------------------------------------------------------------
!  Extended-Lagrangian gradient step and corrector for shell atoms.
!  Uses the force aa (eV/Å, Cartesian) computed at (r_core, θ) to:
!    δra = hi*(κ * aa_shell)    [Cartesian force → fractional displacement]
!    ra_shell = θ + δra         [one Newton step toward equilibrium]
!    θ̈(t+dt) = ω² * δra
!    θ̇!  Also computes the XL auxiliary energy E_aux = E_coupling + E_kinetic:
!    E_coupling = 0.5 * κ * |F_cart|²                                [eV]
!    E_kinetic  = 0.5 / (ω² * κ) * |h*θ̇(t+dt)|²                       [eV]
!  which is derived from the fictitious mass mu = 1/(ω² * κ) [eV fs²/Å²].
!  The conserved quantity is H_XL = T_core + V(r_core,r_shell*) + E_aux.
!
!  Note: xl_thdot must already contain the half-step value from xl_predict.
!-----------------------------------------------------------------------
    implicit none
    include "./params_unit.h"
    integer,intent(in):: natm,tag_isp(natm)
    real(rp),intent(inout):: ra(3,natm)
    real(rp),intent(in):: aa(3,natm),hi(3,3),h(3,3),dt
    real(rp),intent(out):: eaux
 
    integer:: i,is,itot
    real(rp):: dfrac(3),kappa,vcart(3)
 
    eaux = 0.0_rp
    do i=1,natm
      is = tag_isp(i)
      if( .not.is_shell_sp(is) ) cycle
      itot = tag_itot(i)
      kappa = xl_kappa_sp(is)
!     δra in fractional coords = hi * (κ_i * F_cart) with species-dependent κ_i = 1/k2s_i
      dfrac(1) = hi(1,1)*aa(1,i) + hi(1,2)*aa(2,i) + hi(1,3)*aa(3,i)
      dfrac(2) = hi(2,1)*aa(1,i) + hi(2,2)*aa(2,i) + hi(2,3)*aa(3,i)
      dfrac(3) = hi(3,1)*aa(1,i) + hi(3,2)*aa(2,i) + hi(3,3)*aa(3,i)
      dfrac(1:3) = kappa * dfrac(1:3)
!     E_coupling = 0.5*κ*|F_cart|²
!     (Derived from Lagrangian with fictitious mass mu = 1/(omega^2 * kappa))
      eaux = eaux + 0.5_rp * kappa * (aa(1,i)**2 + aa(2,i)**2 + aa(3,i)**2)
!     update shell position (indexed by global atom ID, not array position)
      ra(1:3,i) = xl_theta(1:3,itot) + dfrac(1:3)
!     θ̈ corrector
      xl_thacc(1:3,itot) = xl_omega2 * dfrac(1:3)
!     complete θ̇ (xl_thdot currently holds half-step value from xl_predict)
      xl_thdot(1:3,itot) = xl_thdot(1:3,itot) + 0.5_rp*dt*xl_thacc(1:3,itot)
!     E_kinetic = 0.5*mu*|v_cart|² = 0.5 / (omega^2 * kappa) * |v_cart|²
      vcart(1) = h(1,1)*xl_thdot(1,itot)+h(1,2)*xl_thdot(2,itot)+h(1,3)*xl_thdot(3,itot)
      vcart(2) = h(2,1)*xl_thdot(1,itot)+h(2,2)*xl_thdot(2,itot)+h(2,3)*xl_thdot(3,itot)
      vcart(3) = h(3,1)*xl_thdot(1,itot)+h(3,2)*xl_thdot(2,itot)+h(3,3)*xl_thdot(3,itot)
      eaux = eaux + 0.5_rp / (xl_omega2 * kappa) &
           * (vcart(1)**2 + vcart(2)**2 + vcart(3)**2)
    enddo
    return
  end subroutine xl_gradient_step
!=======================================================================
  subroutine xl_sync_theta(natm,tag_isp,ra)
!-----------------------------------------------------------------------
!  After bamove() wraps shell atoms across periodic boundaries, xl_theta
!  must be updated to match the wrapped positions.  Without this, the
!  next xl_gradient_step would write  ra = xl_theta_OLD + kappa*F,
!  placing the shell on the WRONG side of the boundary, which makes GF2
!  use positions inconsistent with the pairlist (built after bamove).
!
!  We simply copy the current (bamove-corrected) fractional position
!  into xl_theta.  xl_thdot and xl_thacc are unchanged, because adding
!  a lattice vector to a fractional coordinate does not affect velocities
!  or accelerations.
!-----------------------------------------------------------------------
    implicit none
    integer,intent(in):: natm,tag_isp(natm)
    real(rp),intent(in):: ra(3,natm)
    integer:: i,itot
    do i=1,natm
      if( .not.is_shell_sp(tag_isp(i)) ) cycle
      itot = tag_itot(i)
      xl_theta(1:3,itot) = ra(1:3,i)
    enddo
    return
  end subroutine xl_sync_theta
!=======================================================================
  subroutine read_params_ShellModel(myid_md,mpi_md_world,iprint)
!
!  Read spring constants from in.params.ShellModel.
!  File format (one entry per polarisable species pair):
!-----------------------------------------------------------------------
!  #  ShellModel spring constants
!  #  core_species  shell_species  k2_s(eV/Ang^2)  k4_s(eV/Ang^4)
!     O   O_s   74.92   0.0
!     Zr  Zr_s  169.617  0.0
!  (k4_s column is optional; defaults to 0.0 if omitted)
!-----------------------------------------------------------------------
    implicit none
    integer,intent(in):: myid_md,mpi_md_world,iprint

    integer:: isp,jsp,ierr,ios
    real(rp):: k2s,k4s
    character(len=128):: cline,fname
    character(len=5):: cspi,cspj

    if( myid_md.eq.0 ) then
      fname = trim(paramsdir)//'/'//trim(cprmfname)
      open(ioprms,file=trim(fname),status='old')
      if( iprint.ge.ipl_basic ) write(6,'(/,a)') ' ShellModel parameters:'
      do while(.true.)
        read(ioprms,'(a)',end=10) cline
        if( cline(1:1).eq.'#' .or. cline(1:1).eq.'!' ) cycle
        if( len_trim(cline).eq.0 ) cycle
        k4s = 0.0_rp
        read(cline,*,iostat=ios) cspi,cspj,k2s,k4s
        if( ios.ne.0 ) read(cline,*) cspi,cspj,k2s
        isp = csp2isp(cspi)    ! core species index
        jsp = csp2isp(cspj)    ! shell species index
        if( isp.gt.0 .and. jsp.gt.0 ) then
          sm_k2s(jsp)     = k2s
          sm_k4s(jsp)     = k4s
          sm_core_of(jsp) = isp
          sm_shell_of(isp)= jsp
          is_shell_sp(jsp)= .true.
          if( iprint.ge.ipl_basic ) then
            write(6,'(a,2a6,2f12.4)') &
                 '   core,shell,k2_s,k4_s = ',trim(cspi),trim(cspj),k2s,k4s
          endif
        else
          if( iprint.ge.ipl_info ) then
            print *,' ShellModel parameter read but not used: ', &
                 trim(cspi),' ',trim(cspj)
          endif
        endif
      enddo
10    close(ioprms)
    endif

    call mpi_bcast(sm_k2s,nspmax,mpi_real_rp,0,mpi_md_world,ierr)
    call mpi_bcast(sm_k4s,nspmax,mpi_real_rp,0,mpi_md_world,ierr)
    call mpi_bcast(sm_core_of,nspmax,mpi_integer,0,mpi_md_world,ierr)
    call mpi_bcast(sm_shell_of,nspmax,mpi_integer,0,mpi_md_world,ierr)
    call mpi_bcast(is_shell_sp,nspmax,mpi_logical,0,mpi_md_world,ierr)

    lprmset_ShellModel = .true.
    return
  end subroutine read_params_ShellModel

end module ShellModel
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
