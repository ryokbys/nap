module ttm
!-----------------------------------------------------------------------
!                     Last-modified: <2021-01-15 14:57:29 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
!
! Module for two-temperature method (TTM).
!
! In the current implementation, it is assumed that the number of 
! parallel nodes are the common dividors of number of TTM meshes.
!
  implicit none
  save
  include 'mpif.h'
  include "./params_unit.h"
  
  character(len=128):: paramsdir = '.'
  character(len=128),parameter:: cfparams = 'in.ttm'
  character(len=128),parameter:: cTe_infile = 'in.Te'
  character(len=128),parameter:: cTe_outfile = 'out.Te'
  character(len=128),parameter:: cergio = 'out.eio_ttm'
  
  integer,parameter:: ioprms = 30
  integer,parameter:: ioTein = 31
  integer,parameter:: ioTeout = 32
  integer,parameter:: ioergio = 33

!!$  real(8),parameter:: pi = 3.14159265358979d0

  real(8):: t_ttm
!.....TTM mesh divisions
  integer:: nx,ny,nz,nxyz
!.....Mesh size in reduced unit [0:1)
  real(8):: dx,dy,dz,area,darea,dr2
!.....ODE solver; 1) Euler, 2) RK4 (4th order Runge-Kutta), 3) RK2
  character(len=20):: csolver = 'RK4'
!.....Time step in fs == dt in MD by default
  real(8):: dt_inner = -1.0
!.....Number of inner loop in TTM
  integer:: nstp_inner = 1
!.....Volume and area per mesh cell Ang^3 == dx*dy*dz
  real(8):: vcell
!.....Threshold kinetic energy in energy unit (or should be threshold velocity?)
  real(8):: ekth = 8.0d0
!.....Te of right edge, if negative, use free boundary
  real(8):: Te_right = -1d0
!.....Minimum Te, which limits the electron system goes too low energy
  real(8):: Te_min = -100d0
!.....Maximums in the system
  real(8):: Te_max = -1d0
  real(8):: alpha_max = -1d0
  real(8):: kappa_max, cete_min

!.....Ce dependence on Te: none(0), polynomial(1), tanh(2), linear(3)
  character(len=128):: Ce_Tdep = 'none'
  integer:: iCe_Tdep = 0
  real(8):: rho_e = 0.005
  real(8):: d_e = 20d0
  real(8):: c_0 = 1e-4
!.....Coefficients for polynomial Ce(Te)
  real(8):: a_0 = 0d0
  real(8):: a_1 = 0d0
  real(8):: a_2 = 0d0
  real(8):: a_3 = 0d0
  real(8):: a_4 = 0d0
  real(8):: A_exp = 0d0
!.....Coefficients for linear, (Ce(Te) = gmm_ce * Te ), in eV/(Ang^3*K^2)
  real(8):: gmm_ce = 6.648d-9
!.....Minimum Ce value for computational stability, in eV/K
  real(8):: ce_min = 1d-6
!.....Minimum atomic temperature for calculating kappa, in K
  real(8):: ta_min = 10.0

!.....Initial Te distribution: exp, homo, or read (from cTe_init)
  character(len=128):: cTe_init = 'homo'
  real(8):: Te_init = 300d0
!.....T-dependence of pulse shape: stepwise(1, default) or gaussian(2)
  character(len=128):: ctype_pulse = 'stepwise'
  integer:: itype_pulse = 1
!.....Function shape of thermal diffusivity, kappa: DCrho(1, default) or B2(2)
! B2: from Eq.(B2) in PRB68,064114 (2003)
  character(len=128):: ctype_kappa = 'DCrho'
  integer:: itype_kappa = 1
!.....Prefactor in case of kappa_type = B2, in eV/(fs*Ang*K)
  real(8):: kappa0 = 6.2422d-7

!.....Absobed laser fluence in eV/Ang^2 unit
  real(8):: fluence = 0d0
  real(8):: I_0
!.....Start time of laser injection
  real(8):: t0_laser = 0d0
!.....Pulse duration in fs
  real(8):: tau_pulse = 100d0
!.....Pulse sigma in case of Gaussian
  real(8):: sgm_pulse
!.....Surface skin length in Ang
  real(8):: lskin = 100d0
!.....Surface position (positive integer)
  integer:: lsurf = -1  ! left surface
  integer:: rsurf = -1  ! right surface
!.....Surface movement along x: no, plane (can move as a yz-plane surface)
  character(len=128):: surfmove = 'plane'
  integer:: nstp_surfmove = 10
!.....Threshold density: default = 0.025 Ang^{-3} (1/2 of ideal Si diamond density)
!     ex) 8 /5.427^3 = 0.05
  real(8):: dthresh = 0.025d0

!.....Temperature distribution
  real(8),allocatable:: te(:,:,:),tep(:,:,:),ta(:),tap(:),tex(:)
  integer,allocatable:: nac(:),nacp(:)
  real(8),allocatable:: eksum(:),ekpsum(:),vac(:,:),ekti(:)
!.....Atom to cell correspondance
  integer,allocatable:: a2c(:)

!.....Type of coupling constant: constant_gmmp(1) or constant_gp(2)
  character(len=128):: ctype_coupling = 'constant_gmmp'
  integer:: itype_coupling = 1
!.....e-ph coupling constant, g, in case of constant gp in eV/(fs*A^3*K)
!.....Parameter for Ni from Zhigilei et al., J.Phys.Chem. C 113 (2009)
  real(8):: e_ph_const = 2.247e-9
!.....Friction coefficient, gamma, in 1/fs (inverse of relaxation time)
  real(8):: gamma_p = 0.001
  real(8):: gamma_s = 0.1
  real(8),allocatable:: gmmp(:), gmms(:)
!.....Sigma of random force in Langevin thermostat
  real(8),allocatable:: sgm(:)
!.....gp, gs
  real(8),allocatable:: gp(:),gs(:)

  real(8),allocatable:: aai(:,:)

  real(8):: etot_e
!.....For energy difference
  real(8):: ein_e, ein_a, eout_e, eout_a
  real(8),allocatable:: dein(:),deout(:)

!.....DEBUGGING
!.....Cut interaction bewteen atom and electron systems
  logical:: lcut_interact = .false.

contains
!=======================================================================
  subroutine init_ttm(namax,natm,h,dtmd,lvardt,myid,mpi_world,iprint)
!
!  Read parameters for TTM from in.ttm and initialize
!
    integer,intent(in):: namax,natm,myid,mpi_world,iprint
    real(8),intent(in):: dtmd,h(3,3)
    logical,intent(in):: lvardt 

    integer:: ierr,ix,iy,iz,mem
    real(8):: t,t0,dtmax,tmp
    character(len=128):: c1st

    t_ttm = 0d0
    t0 = mpi_wtime()
!.....Read parameter file
    call read_ttm_params(myid,mpi_world,iprint)
    call sync_params(myid,mpi_world,iprint)

!.....Set some
    call set_inner_dt(dtmd)
    nxyz = nx*ny*nz
    dx = h(1,1)/nx
    dy = h(2,2)/ny
    dz = h(3,3)/nz
    vcell = dx*dy*dz
    area = h(2,2)*h(3,3)
    darea = dy*dz
!.....dr2 is used to compute time-step limit
    dr2 = 1d0/(1d0/dx**2 +1d0/dy**2 +1d0/dz**2)
!    print *,dr2

!.....Convert fluence to intensity in unit of energy not temperature, here.
!.....This is going to be converted to temperature by dividing by (rho*Ce)
    if( trim(ctype_pulse).eq.'stepwise' ) then
      itype_pulse = 1
      I_0 = fluence /lskin /tau_pulse *darea !*(2d0/(3d0*rho_e*vcell*fkb))
    else if( trim(ctype_pulse).eq.'gaussian' ) then
!!$      I_0 = fluence /lskin /tau_pulse *sqrt(4d0 *log(2d0) /pi) *darea &
!!$           *(2d0/(3d0*rho_e*vcell*fkb))
      itype_pulse = 2
      I_0 = fluence /lskin /tau_pulse *sqrt(4d0 *log(2d0) /pi) *darea
      sgm_pulse = tau_pulse /sqrt(8d0 *log(2d0))
    else
      if( myid.eq.0 ) then
        print *,'ERROR: pulse_type should be either stepwise or gaussian.'
      endif
      call mpi_finalize(ierr)
      stop
    endif

    if( trim(ctype_kappa).eq.'DCrho' ) then
      itype_kappa = 1
    else if( trim(ctype_kappa).eq.'B2' ) then
      itype_kappa = 2
    endif

!.....Conversion of Ce if needed
    if( trim(Ce_Tdep).eq.'none' ) then
      iCe_Tdep = 0
    else if( trim(Ce_Tdep).eq.'polynomial' ) then
      iCe_Tdep = 1
    else if( trim(Ce_Tdep).eq.'tanh' ) then
      iCe_Tdep = 2
    else if( trim(Ce_Tdep).eq.'linear' ) then
      iCe_Tdep = 3
      gmm_ce = gmm_ce !/rho_e
    endif

    if( myid.eq.0 .and. iprint.ne.0 ) then
      print *,''
      print *,'TTM parameters:'
      print '(a,es12.4,a,f8.4,a)','   Fluence = ',fluence,' eV/A^2, = ' &
           ,fluence*ev2j/(ang2m**2*10000),' J/cm^2'
      print '(a,f0.3,a)','   Pulse duration = ',tau_pulse,' fs'
      print '(a,es12.4,a)','   Intensity = ',I_0/darea,' eV/A^2/fs'
      print '(a,es12.4,a)','   Penetration depth = ',lskin,' A'
      print '(a,f0.1,a)','   Total incident energy = ',fluence*area,' eV'
      print '(a,f0.4,a)','   Electron density = ',rho_e,' e/A^3'
      print '(a,a)','   Diff. Eq. solver :  ',trim(csolver)
      print '(a,i5)','   inner_loop = ',nstp_inner
      if( .not. lvardt ) then
        print '(a,2es12.4)','   dtmd, dt = ',dtmd,dt_inner
      else
        print '(a,2es12.4,a)','   dtmd, dt = ',dtmd,dt_inner, &
             ' fs (but it is variable)'
      endif
      print '(a,3i5,i8)','   nx,ny,nz,nxyz = ',nx,ny,nz,nxyz
      print '(a,4es12.4)','   dx,dy,dz,vcell = ',dx,dy,dz,vcell
      print '(a,2es12.4)','   area,darea = ',area,darea
      print '(a,2i5)','   lsurf,rsurf = ',lsurf,rsurf
      if( trim(ctype_coupling).eq.'constant_gp' ) then
        print '(a,2es12.4)','   e_ph_const = ',e_ph_const
      else if( trim(ctype_coupling).eq.'constant_gmmp' ) then
        print '(a,2es12.4)','   gamma_p,gamma_s = ',gamma_p,gamma_s
      endif
      print '(a,a)','   Ce_Tdep = ',trim(Ce_Tdep)
      if( trim(Ce_Tdep).eq.'linear' ) then
        print '(a,2es12.4)','     gmm_ce,ce_min = ',gmm_ce, ce_min
      endif
      if( Te_min.gt.0d0 ) then
        print *,'    Note that the minimum Te is set:'
        print '(a,f0.1)','     Minimum Te = ',Te_min
      endif
      mem = 4 * 4*nxyz + 11 * 8*nxyz + 4 * 8*namax
      print '(a,f0.3,a)','   Memory for TTM = ',dble(mem)/1000/1000,' MByte'
    endif
    
!.....Error check
    if( nxyz.le.0 ) then
      if( myid.eq.0 ) then
        print *,'ERROR: wrong nxyz: ',nxyz
        print *,'  nx,ny,nz = ',nx,ny,nz
      endif
      goto 999
    else if( rsurf.lt.lsurf ) then
      if( myid.eq.0 ) then
        print *,'ERROR: rsurf < lsurf, which should not occur !'
        print *,'  lsurf,rsurf = ',lsurf,rsurf
      endif
      goto 999
    else if( rsurf.gt.nx ) then
      if( myid.eq.0 ) then
        print *,'ERROR: rsurf > nx, which should not occur !'
        print *,'  rsurf,nx = ',rsurf,nx
      endif
      goto 999
    endif

!!$!.....Check dt for TTM
!!$    t = 10d0
!!$    tmp = c_0 +(a_0 +a_1*t +a_2*t**2 +a_3*t**3 +a_4*t**4)&
!!$         *exp(-(A_exp*t)**2)
!!$    dtmax = 0.5d0/(D_e*rho_e*tmp) /(1d0/dx**2 +1d0/dy**2 +1d0/dz**2)
!!$    print *,'dt_ttm,dtmax = ',dt,dtmax
!!$    if( dt.gt.dtmax ) then
!!$      if( myid.eq.0 ) then
!!$        print *,'ERROR: dt > dtmax, you may have to change TTM mesh size.'
!!$        print *,'dt,dtmax = ',dt,dtmax
!!$      endif
!!$    endif
    
!.....Allocate initialize arrays
    allocate(nac(nxyz),nacp(nxyz),eksum(nxyz),ekpsum(nxyz), &
         sgm(nxyz),te(0:nx+1,0:ny+1,0:nz+1),tep(0:nx+1,0:ny+1,0:nz+1), &
         ta(nxyz),tap(nxyz),tex(nx),gp(nxyz),gs(nxyz),&
         gmmp(nxyz),gmms(nxyz),vac(3,nxyz))
    allocate(a2c(namax),aai(3,namax),ekti(namax))

    if( trim(ctype_coupling).eq.'constant_gmmp' ) then
      itype_coupling = 1
      gmmp(:) = gamma_p
      gmms(:) = gamma_s
    else if( trim(ctype_coupling).eq.'constant_gp' ) then
      itype_coupling = 2
      gp(:) = e_ph_const
      gs(:) = 0d0
    else 
      print *,'ERROR: no such coupling style: '//trim(ctype_coupling)
      goto 999
    endif

!.....Set initial Te distribution
    te(:,:,:) = 0d0
    if( trim(cTe_init).eq.'exp' ) then
      if( lsurf.le.0 ) then
        if( myid.eq.0 ) then
          print *,'ERROR: Initial Te distribution cannot be set with exp'&
               //' with lsurf <= 0.'
        endif
        goto 999
      endif
!!$      do ix=lsurf,rsurf
!!$        te(ix,1:ny,1:nz) = I_0 *exp(-dx*(ix-lsurf+1)/lskin)
!!$      enddo
    else if( cTe_init(1:4).eq.'homo' ) then  ! homogeneous Te
      if( myid.eq.0 ) then
        print '(a,f7.1,a)','   Initial Te =     ',Te_init,' K'
      endif
    else if( trim(cTe_init).eq.'read' ) then
      if( myid.eq.0 ) then
        print '(a)','   Initial Te is read from in.Te.'
        open(ioTein,file=trim(paramsdir)//'/'//trim(cTe_infile),status='old')
        do while(.true.)
          read(ioTein,*,end=10) c1st
          if( c1st(1:1).eq.'!' .or. c1st(1:1).eq.'#' ) cycle
          backspace(ioTein)
          read(ioTein,*) ix,iy,iz,t
          te(ix,iy,iz) = t
        enddo
10      close(ioTein)
      endif
!.....Broadcast Te distribution to all the nodes.
!.....There could be smarter way to reduce networking cost.
      call mpi_bcast(te,(nx+2)*(ny+2)*(nz+2),mpi_real8,0,mpi_world,ierr)
    endif

    t_ttm = t_ttm +mpi_wtime() -t0
    return

999 call mpi_finalize(ierr)
    stop
  end subroutine init_ttm
!=======================================================================
  subroutine set_inner_dt(dtmd)
    real(8),intent(in):: dtmd
    real(8):: tmp

!.....If alpha_max > 0, the upper limit of dt_inner can be determined.
    if( alpha_max.gt.0d0 ) then
      tmp = dr2 /(2d0*alpha_max)
!!$      if( dt_inner.gt.tmp ) dt_inner = tmp
      dt_inner = tmp
    endif

!.....If dt_inner is specified, nstp_inner and dt are determined from dt_inner and dtmd.
!.....Since the upper limit of dt is given by dx**2/(2*kappa),
!.....dt_inner should be determined regardless to dtmd and
!.....this would be appropriate especially in case of variable time-step.
    if( dt_inner .gt. 0d0 ) then
!!$      nstp_inner = max(int(dtmd /dt_inner),1)
      nstp_inner = int(dtmd /dt_inner) +1
    endif
    tmp = dtmd /nstp_inner
!.....Change dt_inner to match dt and use this dt_inner afterward
    dt_inner = tmp

!!$    print '(a,i4,5es12.4)','nstp_inner,dt_inner,alpha,Te,kappa,cete='&
!!$         ,nstp_inner,dt_inner,alpha_max,Te_max,kappa_max,cete_min
    
    return
  end subroutine set_inner_dt
!=======================================================================
  subroutine read_ttm_params(myid,mpi_world,iprint)
    use util, only: num_data
    integer,intent(in):: myid,mpi_world,iprint
!!$    integer,external:: num_data

    character(len=128):: c1st,fname
    
    if( myid.eq.0 ) then
      fname = trim(paramsdir)//'/'//trim(cfparams)
      open(ioprms,file=trim(fname),status='old')

      do while(.true.)
        read(ioprms,*,end=10) c1st
        if( num_data(c1st,' ').eq.0 ) cycle
        if( c1st(1:1).eq.'!' .or. c1st(1:1).eq.'#' ) cycle
        if( trim(c1st).eq.'mesh_size' ) then
          backspace(ioprms)
          read(ioprms,'(a)') c1st
          if( num_data(c1st,' ').ne.4 ) then
            print *,'Error: Wrong number of parameters for mesh_size !'
          endif
          backspace(ioprms)
          read(ioprms,*) c1st, nx, ny, nz
        else if( trim(c1st).eq.'inner_loop' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, nstp_inner
        else if( trim(c1st).eq.'inner_dt' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, dt_inner
        else if( trim(c1st).eq.'gamma_p' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, gamma_p
        else if( trim(c1st).eq.'gamma_s' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, gamma_s
        else if( trim(c1st).eq.'ekth' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, ekth
        else if( trim(c1st).eq.'rho_e' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, rho_e
        else if( trim(c1st).eq.'D_e' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, d_e
        else if( trim(c1st).eq.'kappa_type' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, ctype_kappa
        else if( trim(c1st).eq.'kappa_0' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, kappa0
        else if( trim(c1st).eq.'initial_Te' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, cTe_init
        else if( trim(c1st).eq.'Te_init' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, Te_init
        else if( trim(c1st).eq.'Te_right' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, Te_right
        else if( trim(c1st).eq.'t0_laser' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, t0_laser
        else if( trim(c1st).eq.'pulse_type' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, ctype_pulse
        else if( trim(c1st).eq.'pulse_duration' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, tau_pulse
        else if( trim(c1st).eq.'coupling_type' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, ctype_coupling
        else if( trim(c1st).eq.'coupling_constant' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, e_ph_const
        else if( trim(c1st).eq.'Ce_T-depend' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, ce_Tdep
        else if( trim(c1st).eq.'C_0' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, c_0
        else if( trim(c1st).eq.'Ce_poly_params' ) then
          backspace(ioprms)
          read(ioprms,'(a)') c1st
          if( num_data(c1st,' ').ne.7 ) then
            print *,'Error: Wrong number of parameters for Ce_poly_params !'
          endif
          backspace(ioprms)
          read(ioprms,*) c1st, a_0,a_1,a_2,a_3,a_4,A_exp
        else if( trim(c1st).eq.'Ce_linear_gamma' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, gmm_ce
        else if( trim(c1st).eq.'surface_move' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, surfmove
        else if( trim(c1st).eq.'threshold_density' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, dthresh
        else if( trim(c1st).eq.'laser_fluence' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, fluence
        else if( trim(c1st).eq.'Te_min' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, Te_min
        else if( trim(c1st).eq.'left_surface' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, lsurf
        else if( trim(c1st).eq.'right_surface' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, rsurf
        else if( trim(c1st).eq.'surface_skin_length' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, lskin
        else if( trim(c1st).eq.'cut_interaction' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, lcut_interact
        else if( trim(c1st).eq.'Ta_min' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, ta_min
        else if( trim(c1st).eq.'DE_solver' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, csolver
        else
          print *,'There is no TTM keyword: ',trim(c1st)
        endif
      enddo

10    close(ioprms)
    endif
  end subroutine read_ttm_params
!=======================================================================
  subroutine sync_params(myid,mpi_world,iprint)
    integer,intent(in):: myid,mpi_world,iprint
    integer:: ierr
    call mpi_bcast(nx,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(ny,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(nz,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(nstp_inner,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(dt_inner,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(gamma_p,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(gamma_s,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(ekth,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(rho_e,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(d_e,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(ctype_kappa,128,mpi_character,0,mpi_world,ierr)
    call mpi_bcast(kappa0,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(t0_laser,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(ctype_pulse,128,mpi_character,0,mpi_world,ierr)
    call mpi_bcast(tau_pulse,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(ctype_coupling,128,mpi_character,0,mpi_world,ierr)
    call mpi_bcast(e_ph_const,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(cTe_init,128,mpi_character,0,mpi_world,ierr)
    call mpi_bcast(Te_right,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(Te_init,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(ce_Tdep,128,mpi_character,0,mpi_world,ierr)
    call mpi_bcast(c_0,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(a_0,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(a_1,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(a_2,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(a_3,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(a_4,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(A_exp,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(gmm_ce,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(ce_min,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(ta_min,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(surfmove,1,mpi_logical,0,mpi_world,ierr)
    call mpi_bcast(fluence,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(Te_min,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(lsurf,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(rsurf,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(lskin,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(lcut_interact,1,mpi_logical,0,mpi_world,ierr)
    call mpi_bcast(csolver,20,mpi_character,0,mpi_world,ierr)
    return
  end subroutine sync_params
!=======================================================================
  subroutine assign_atom2cell(namax,natm,ra,sorg,boundary)
!
!  Assign atoms to TTM mesh cell
!
    integer,intent(in):: namax,natm
    real(8),intent(in):: ra(3,namax),sorg(3)
    character(len=3),intent(in):: boundary

    integer:: i,ix,iy,iz,ic,l
    real(8):: xi(3),udx,udy,udz,t0

    t0 = mpi_wtime()

    udx = 1d0/nx
    udy = 1d0/ny
    udz = 1d0/nz
    a2c(:) = 0
    do i=1,natm
      xi(1:3) = ra(1:3,i) +sorg(1:3)
      do l=1,3
        if( boundary(l:l).eq.'p' ) then
          xi(l) = mod(xi(l),1d0)
        endif
      enddo
      ix = int(xi(1)/udx) +1
      iy = int(xi(2)/udy) +1
      iz = int(xi(3)/udz) +1
      if( boundary(1:1).eq.'f' ) then
        ix = min(max(ix,1),nx)
      endif
      call ixyz2ic(ix,iy,iz,ic)
      a2c(i) = ic
    enddo

    t_ttm = t_ttm +mpi_wtime()-t0
    
  end subroutine assign_atom2cell
!=======================================================================
  subroutine calc_Ta(namax,natm,nspmax,h,tag,va,fmv,fekin &
       ,istp,myid,mpi_world,iprint)
!
!  Compute and set Ta and Tap array from atomic kinetic energies.
!
    use util,only: ifmvOf
    integer,intent(in):: namax,natm,nspmax,myid,mpi_world,istp,iprint
    real(8),intent(in):: tag(namax),fmv(3,0:9),va(3,namax),h(3,3) &
         ,fekin(nspmax)

    integer:: i,ic,ierr,is,ix,iy,iz,l,ifmv,idof
    real(8):: ek,t0,vat(3),vatr(3)
    integer,allocatable,save:: nacl(:),nacpl(:)
    real(8),allocatable,save:: eksuml(:),ekpsuml(:),vacl(:,:)
!!$    integer,external:: ifmvOf

    if( .not. allocated(nacl) ) then
      allocate(nacl(nxyz),nacpl(nxyz),eksuml(nxyz),ekpsuml(nxyz)&
           ,vacl(3,nxyz))
    endif

    t0 = mpi_wtime()

!.....First distinguish center of mass vectors of cells
    vacl(:,:) = 0d0
    nacl(:) = 0
    do i=1,natm
      ic = a2c(i)
      vacl(1:3,ic) = vacl(1:3,ic) +va(1:3,i)
      nacl(ic) = nacl(ic) + 1
    enddo
    vac(:,:) = 0d0
!.....NACL and NAC are temporal here, used just for normalization
    nac(:) = 0
    call mpi_reduce(vacl,vac,3*nxyz,mpi_real8,mpi_sum,0,mpi_world,ierr)
    call mpi_reduce(nacl,nac,nxyz,mpi_integer,mpi_sum,0,mpi_world,ierr)
    do ic=1,nxyz
      if( nac(ic).eq.0 ) cycle
      vac(1:3,ic) = vac(1:3,ic) /nac(ic)
    enddo
    call mpi_bcast(vac,3*nxyz,mpi_real8,0,mpi_world,ierr)

!.....Compute Ekin per atom using thermal part of velocities
    do i=1,natm
      ic = a2c(i)
      is = int(tag(i))
      vat(1:3) = va(1:3,i) -vac(1:3,ic)
      vatr(1:3) = h(1:3,1)*vat(1) +h(1:3,2)*vat(2) +h(1:3,3)*vat(3)
      ekti(i) = (vatr(1)**2 +vatr(2)**2 +vatr(3)**2) *fekin(is)
    enddo
    
    nacl(1:nxyz) = 0
    nacpl(1:nxyz) = 0
    eksuml(1:nxyz) = 0d0
    ekpsuml(1:nxyz) = 0d0
    do i=1,natm
      ic = a2c(i)
      ifmv = ifmvOf(tag(i))
      idof = 0
      do l=1,3
        idof = idof +nint(fmv(l,ifmv))
      enddo
      nacl(ic) = nacl(ic) + idof
      eksuml(ic) = eksuml(ic) +ekti(i)
      if( ek.gt.ekth ) then
        nacpl(ic) = nacpl(ic) +idof
        ekpsuml(ic) = ekpsuml(ic) +ekti(i)
      endif
    enddo
    nac(1:nxyz) = 0
    nacp(1:nxyz) = 0
    eksum(1:nxyz) = 0d0
    ekpsum(1:nxyz) = 0d0
    call mpi_reduce(nacl,nac,nxyz,mpi_integer,mpi_sum,0,mpi_world,ierr)
    call mpi_reduce(nacpl,nacp,nxyz,mpi_integer,mpi_sum,0,mpi_world,ierr)
    call mpi_reduce(eksuml,eksum,nxyz,mpi_real8,mpi_sum,0,mpi_world,ierr)
    call mpi_reduce(ekpsuml,ekpsum,nxyz,mpi_real8,mpi_sum,0,mpi_world,ierr)
!.....Compute Ta and Tap only at node-0
    if( myid.eq.0 ) then
      if( trim(ctype_coupling).eq.'constant_gmmp' ) then
        gp(:) = 0d0
        gs(:) = 0d0
        ta(:) = 0d0
        tap(:) = 0d0
        do ic=1,nxyz
          if( nac(ic).eq.0 ) cycle
!!$          call ic2ixyz(ic,ix,iy,iz)
!.....Degree of freedom per atom (3 in case of 3D) is included in nac
!.....CHECK: This factor 3 looks causing the difference of energy in/out between at/el systems.
          ta(ic) = eksum(ic) *2d0 /fkb /nac(ic)
          gp(ic) = nac(ic) *fkb *gmmp(ic) /vcell ! /3
          if( nacp(ic).eq.0 ) cycle
          tap(ic) = ekpsum(ic) *2d0 /fkb /nacp(ic)
          gs(ic) = nacp(ic) *fkb *gmms(ic) /vcell ! /3
        enddo
      else if( trim(ctype_coupling).eq.'constant_gp' ) then
!.....See Eq.(A5) in PRB 68 (2003) pp.064114
        gmmp(:) = 0d0
        gmms(:) = 0d0
        ta(:) = 0d0
        tap(:) = 0d0
        do ic=1,nxyz
          if( nac(ic).eq.0 ) cycle
          call ic2ixyz(ic,ix,iy,iz)
          ta(ic) = eksum(ic) *2d0 /fkb /nac(ic)
          if( tap(ic)*0d0 .ne. 0d0 ) then
            print *,'ERROR: tap==NaN !!!'
            print *,'   ic,ix,iy,iz=',ic,ix,iy,iz
            stop
          endif
!.....The definition from PRB 68 (2003) seems to have malfunctioning
!.....in case (Te-Ta) < 0.0, which causes sqrt(negative) for sigma calculation
!.....in Langevin_ttm,
!          gmmp(ic) = gp(ic)*vcell*(te(ix,iy,iz)-ta(ic))/2/eksum(ic)
!.....Here inverse of constant_gmmp will be used.
          gmmp(ic) = vcell*gp(ic) /fkb /nac(ic)
        enddo
      endif
    endif

    if( trim(surfmove).eq.'plane' .and. &
         (istp.eq.0 .or. mod(istp,nstp_surfmove).eq.0) ) then
      call update_surface_plane(myid,mpi_world,iprint)
    endif

    t_ttm = t_ttm +mpi_wtime() -t0
    return
  end subroutine calc_Ta
!=======================================================================
  subroutine model_dte(tnow,dtep,eitmp,eotmp,eptmp,iprint)
!
!  Create f(t,y) for ODE dy/dt = f(t,y), where y=Te(ix,iy,iz) here.
!  The ODE is now diffusion Eq. with two-temperature model.
!
    real(8),intent(in):: tnow
    integer,intent(in):: iprint
    real(8),intent(out):: dtep(0:nx+1,0:ny+1,0:nz+1),eitmp,eotmp,eptmp

    integer:: ic,ix,iy,iz
    real(8):: ce,dce,kappa,dkappa,pterm,sterm,dtemp,de,tmp&
         ,pulsefactor,xi
    
    dtep(:,:,:) = 0d0
    eitmp = 0d0
    eotmp = 0d0
    eptmp = 0d0
    call set_BC()
    do ic=1,nxyz
      call ic2ixyz(ic,ix,iy,iz)
      if( ix.lt.lsurf ) cycle
      if( ix.gt.rsurf ) cycle
      ce = cete(ix,iy,iz)
      dce = dcete(ix,iy,iz)
!!$      if( trim(ctype_kappa).eq.'DCrho' ) then
      if( itype_kappa.eq.1 ) then
        kappa = d_e *ce *rho_e
        dkappa = d_e *dce *rho_e
!!$      else if( trim(ctype_kappa).eq.'B2' ) then
      else if( itype_kappa.eq.2 ) then
        kappa = kappa0 *tep(ix,iy,iz) /max(ta(ic),ta_min)
        dkappa = kappa0 /max(ta(ic),ta_min)
      endif
      pterm = -gp(ic)*(tep(ix,iy,iz) -ta(ic))
      sterm = gs(ic)*tap(ic)
      if( lcut_interact ) then
        pterm = 0d0
        sterm = 0d0
      else
        eitmp = eitmp +(gp(ic)*ta(ic)+sterm) *vcell  !*dt
        eotmp = eotmp -gp(ic)*tep(ix,iy,iz) *vcell  !*dt
      endif
      dtemp = 1d0/((ce+tep(ix,iy,iz)*dce)*rho_e) &
           *( dkappa*dte2(ix,iy,iz) +kappa*d2te(ix,iy,iz) &
           +pterm +sterm )  ! *dt
      dtep(ix,iy,iz) = dtep(ix,iy,iz) +dtemp
    enddo  ! ic=1,nxyz
!!$    if( ctype_pulse(1:4).eq.'step' ) then
    if( itype_pulse.eq.1 ) then  ! step-wise pulse
      if( tnow.ge.t0_laser .and. &
           tnow.le.(t0_laser +tau_pulse) ) then
        do ic=1,nxyz
          call ic2ixyz(ic,ix,iy,iz)
          if( ix.lt.lsurf ) cycle
          if( ix.gt.rsurf ) cycle
          ce = cete(ix,iy,iz)
          dce = dcete(ix,iy,iz)
!.....To think the cell position is the center of the cell, add 0.5
          xi = (ix-lsurf+0.5d0)*dx
          tmp = 1d0 /((ce+tep(ix,iy,iz)*dce)*rho_e *vcell)
          de = I_0 *min(1d0,exp(-xi/lskin))*dx  !*dt
          dtep(ix,iy,iz) = dtep(ix,iy,iz) +de*tmp
          if( dtep(ix,iy,iz)*0d0 .ne. 0d0 ) then
            print *,'ERROR: tep==NaN !!!'
            print *,'  ic,ix,iy,iz=',ic,ix,iy,iz
            stop
          endif
          eptmp = eptmp +de
        enddo
      endif
!!$    else if( ctype_pulse(1:5).eq.'gauss' ) then
    else if( itype_pulse.eq.2 ) then  ! gaussian pulse
      if( tnow.ge.t0_laser .and. &
           tnow.lt.(t0_laser +tau_pulse*2) ) then
        pulsefactor = exp(-(tnow -(t0_laser+tau_pulse))**2 /(2d0*sgm_pulse**2))
        do ic=1,nxyz
          call ic2ixyz(ic,ix,iy,iz)
          if( ix.lt.lsurf ) cycle
          if( ix.gt.rsurf ) cycle
          ce = cete(ix,iy,iz)
          dce = dcete(ix,iy,iz)
          tmp = 1d0 /((ce+tep(ix,iy,iz)*dce)*rho_e *vcell)
!.....To think the cell position is the center of the cell, add 0.5
          xi = (ix-lsurf+0.5d0)*dx
          de = I_0 *min(1d0,exp(-xi/lskin))*dt_inner*dx *pulsefactor
          dtep(ix,iy,iz) = dtep(ix,iy,iz) +de*tmp
          eptmp = eptmp +de
        enddo
      endif
    endif

  end subroutine model_dte
!=======================================================================
  subroutine update_Te_old(tnow,dtmd,myid,mpi_world,iprint)
!
!  Update Te by solving the diffusion equation.
!  Currently use Euler method for time update.
!
    integer,intent(in):: myid,mpi_world,iprint
    real(8),intent(in):: tnow,dtmd

    integer:: ic,ix,iy,iz,ierr,istp
    real(8):: t0,ce,dce,xi,pterm,sterm,kappa,dkappa,pulsefactor&
         ,dtemp,tmp,de
    real(8):: de_surf, dte_surf
    real(8),save:: ein_pulse,dte_sum
    logical,save:: l1st = .true.

    t0 = mpi_wtime()

    if( l1st ) then
      ein_pulse = 0d0
      dte_sum = 0d0
!.....Get max Te
      if( myid.eq.0 ) then
        do ic=1,nxyz
          call ic2ixyz(ic,ix,iy,iz)
          if( te(ix,iy,iz).gt.Te_max ) Te_max = te(ix,iy,iz)
          if( itype_kappa.eq.1 ) then  ! DCrho
            kappa = d_e *cete(ix,iy,iz) *rho_e
          else if( itype_kappa.eq.2 ) then  ! B2
            kappa = kappa0 *te(ix,iy,iz)/max(ta(ic),ta_min)
          endif
          alpha_max = max(alpha_max,kappa/cete(ix,iy,iz)/rho_e)
        enddo
      endif
      l1st = .false.
    else
!.....Reset inner dt according to the updated alpha_max
      if( myid.eq.0 ) call set_inner_dt(dtmd)
    endif

    ein_e = 0d0
    eout_e = 0d0
    etot_e = 0d0
    if( myid.eq.0 ) then
      
10    continue
      do istp = 1,nstp_inner
        call set_BC()

        do ic=1,nxyz
          call ic2ixyz(ic,ix,iy,iz)
          tep(ix,iy,iz) = 0d0
          if( ix.lt.lsurf ) cycle
          if( ix.gt.rsurf ) cycle
          ce = cete(ix,iy,iz)
          dce = dcete(ix,iy,iz)
!!$          if( trim(ctype_kappa).eq.'DCrho' ) then
          if( itype_kappa.eq.1 ) then
            kappa = d_e *ce *rho_e
            dkappa = d_e *dce *rho_e
!!$          else if( trim(ctype_kappa).eq.'B2' ) then
          else if( itype_kappa.eq.2 ) then
            kappa = kappa0 *te(ix,iy,iz) /max(ta(ic),ta_min)
            dkappa = kappa0 /max(ta(ic),ta_min)
          endif
          pterm = -gp(ic)*(te(ix,iy,iz) -ta(ic))
          sterm = gs(ic)*tap(ic)
          if( lcut_interact ) then
            pterm = 0d0
            sterm = 0d0
          else
            ein_e = ein_e +(gp(ic)*ta(ic)+sterm)*vcell *dt_inner
            eout_e = eout_e -gp(ic)*te(ix,iy,iz)*vcell *dt_inner
          endif
          dtemp = dt_inner/((ce+te(ix,iy,iz)*dce)*rho_e) &
               *( dkappa*dte2(ix,iy,iz) +kappa*d2te(ix,iy,iz) &
               +pterm +sterm )
          tep(ix,iy,iz) = te(ix,iy,iz) +dtemp
          if( Te_min.gt.0d0 ) tep(ix,iy,iz) = max(tep(ix,iy,iz),Te_min)
          if( tep(ix,iy,iz).lt.0d0 ) then
            dt_inner = dt_inner/2
            call set_inner_dt(dtmd)
            print *,'ERROR: Since Te<0, change dt_inner and/or nstp_inner=',dt_inner,nstp_inner
!!$            goto 10  ! Go back and redo inner loop with the half dt_inner
            stop
          endif
        enddo  ! ic=1,nxyz
!!$        if( ctype_pulse(1:4).eq.'step' ) then
        if( itype_pulse.eq.1 ) then  ! step-wise pulse
          if( tnow.ge.t0_laser .and. &
               tnow.le.(t0_laser +tau_pulse) ) then
            do ic=1,nxyz
              call ic2ixyz(ic,ix,iy,iz)
              if( ix.lt.lsurf ) cycle
              if( ix.gt.rsurf ) cycle
              ce = cete(ix,iy,iz)
              dce = dcete(ix,iy,iz)
!.....To think the cell position is the center of the cell, add 0.5
              xi = (ix-lsurf+0.5d0)*dx
              tmp = 1d0 /((ce+te(ix,iy,iz)*dce)*rho_e *vcell)
              de = I_0 *min(1d0,exp(-xi/lskin))*dx *dt_inner
              tep(ix,iy,iz) = tep(ix,iy,iz) +de*tmp
              if( tep(ix,iy,iz)*0d0 .ne. 0d0 ) then
                print *,'ERROR: tep==NaN !!!'
                print *,'  ic,ix,iy,iz=',ic,ix,iy,iz
                stop
              endif
              ein_pulse = ein_pulse +de
            enddo
          endif
!!$        else if( ctype_pulse(1:5).eq.'gauss' ) then
        else if( itype_pulse.eq.2 ) then  ! gaussian shape pulse
          if( tnow.ge.t0_laser .and. &
               tnow.lt.(t0_laser +tau_pulse*2) ) then
            pulsefactor = exp(-(tnow -(t0_laser+tau_pulse))**2 /(2d0*sgm_pulse**2))
            do ic=1,nxyz
              call ic2ixyz(ic,ix,iy,iz)
              if( ix.lt.lsurf ) cycle
              if( ix.gt.rsurf ) cycle
              ce = cete(ix,iy,iz)
              dce = dcete(ix,iy,iz)
              tmp = 1d0 /((ce+te(ix,iy,iz)*dce)*rho_e *vcell)
!.....To think the cell position is the center of the cell, add 0.5
              xi = (ix-lsurf+0.5d0)*dx
              de = I_0 *min(1d0,exp(-xi/lskin))*dt_inner*dx *pulsefactor
              tep(ix,iy,iz) = tep(ix,iy,iz) +de*tmp
              ein_pulse = ein_pulse +de
            enddo
          endif
        endif
        te(:,:,:) = tep(:,:,:)
      enddo  ! istp=1,nstp_inner
!.....Update Te_max and alpha_max
      Te_max = -1d0
      alpha_max = -1d0
      cete_min = 1d30
      kappa_max = -1d0
      do ic=1,nxyz
        call ic2ixyz(ic,ix,iy,iz)
        ce = cete(ix,iy,iz)
        if( itype_kappa.eq.1 ) then  ! DCrho
          kappa = d_e *ce *rho_e
        else if( itype_kappa.eq.2 ) then  ! B2
          kappa = kappa0 *te(ix,iy,iz) /max(ta(ic),ta_min)
        endif
        kappa_max = max(kappa_max,kappa)
        cete_min = min(ce,cete_min)
        alpha_max = max(kappa/ce/rho_e,alpha_max)
        Te_max = max(Te_max,te(ix,iy,iz))
        etot_e = etot_e +te(ix,iy,iz)*fkb*rho_e*vcell*1.5d0
      enddo
!!$      print *,'cete_min,Te_max=',cete_min,Te_max
!.....Output
      if(iprint.gt.1) then
        if( ( itype_pulse.eq.2 .and. &  ! gaussian
             (tnow.ge.t0_laser .and. tnow.le.(t0_laser+tau_pulse*2)) ) &
             .or. ( itype_pulse.eq.1 .and. &  ! stepwise
             (tnow.ge.t0_laser .and. tnow.le.t0_laser+tau_pulse) ) ) then
          print '(a,2es13.4e3)',' tnow,ein_pulse=',tnow,ein_pulse
        endif
        do ic=1,nxyz
          call ic2ixyz(ic,ix,iy,iz)
          if( te(ix,iy,iz).lt.0d0 ) then
            print *,'ERROR: te(ix,iy,iz) < 0 !!'
            print *,'ic,ix,iy,iz,te=',ic,ix,iy,iz,te(ix,iy,iz)
            stop 1
          endif
        enddo
      endif
    endif
!.....Broadcast Te distribution to all the nodes.
!.....There could be smarter way to reduce networking cost.
    call mpi_bcast(te,(nx+2)*(ny+2)*(nz+2),mpi_real8,0,mpi_world,ierr)

    t_ttm = t_ttm +mpi_wtime()-t0

    return
  end subroutine update_Te_old
!=======================================================================
  subroutine update_Te(tnow,dtmd,myid,mpi_world,iprint)
!
!  Update Te by solving the diffusion equation.
!  Model calculation is separated as model_XXX routines,
!  and some ODE solvers are (to be) implemented.
!
    integer,intent(in):: myid,mpi_world,iprint
    real(8),intent(in):: tnow,dtmd

    integer:: ic,ix,iy,iz,ierr,istp
    real(8):: t0,ce,dce,xi,pterm,sterm,kappa,dkappa,pulsefactor&
         ,dtemp,tmp,de,eitmp,eotmp,eptmp
    real(8):: de_surf, dte_surf
    real(8),save:: ein_pulse,dte_sum
    logical,save:: l1st = .true.
    real(8),allocatable,save:: dtep(:,:,:)

    t0 = mpi_wtime()

    if( l1st ) then
      ein_pulse = 0d0
      dte_sum = 0d0
!.....Get max Te
      if( myid.eq.0 ) then
        do ic=1,nxyz
          call ic2ixyz(ic,ix,iy,iz)
          if( te(ix,iy,iz).gt.Te_max ) Te_max = te(ix,iy,iz)
          if( itype_kappa.eq.1 ) then  ! DCrho
            kappa = d_e *cete(ix,iy,iz) *rho_e
          else if( itype_kappa.eq.2 ) then  ! B2
            kappa = kappa0 *te(ix,iy,iz)/max(ta(ic),ta_min)
          endif
          alpha_max = max(alpha_max,kappa/cete(ix,iy,iz)/rho_e)
        enddo
        allocate(dtep(0:nx+1,0:ny+1,0:nz+1))
      endif
      l1st = .false.
    else
!.....Reset inner dt according to the updated alpha_max
      if( myid.eq.0 ) call set_inner_dt(dtmd)
    endif

    ein_e = 0d0
    eout_e = 0d0
    etot_e = 0d0
    if( myid.eq.0 ) then

      if( trim(csolver).eq.'Euler' ) then
        do istp = 1,nstp_inner
          tep(:,:,:) = te(:,:,:)
          call model_dte(tnow,dtep,eitmp,eotmp,eptmp,iprint)
          te(:,:,:) = te(:,:,:) +dtep(:,:,:)*dt_inner
          ein_e = ein_e +eitmp*dt_inner
          eout_e = eout_e +eotmp*dt_inner
          ein_pulse = ein_pulse +eptmp*dt_inner
        enddo  ! istp=1,nstp_inner
      else if( trim(csolver).eq.'RK4' ) then  ! 4th Runge-Kutta
        do istp=1,nstp_inner
          tep(:,:,:)= te(:,:,:)
!.....1st step
          call model_dte(tnow,dtep,eitmp,eotmp,eptmp,iprint)
          te(:,:,:) = te(:,:,:) +dtep(:,:,:)*dt_inner/6
          tep(:,:,:)= tep(:,:,:) +dtep(:,:,:)*dt_inner/2
          ein_e = ein_e +eitmp*dt_inner/6
          eout_e = eout_e +eotmp*dt_inner/6
          ein_pulse = ein_pulse +eptmp*dt_inner/6
!.....2nd step
          call model_dte(tnow,dtep,eitmp,eotmp,eptmp,iprint)
          te(:,:,:) = te(:,:,:) +dtep(:,:,:)*dt_inner/3
          tep(:,:,:)= tep(:,:,:) +dtep(:,:,:)*dt_inner/2
          ein_e = ein_e +eitmp*dt_inner/3
          eout_e = eout_e +eotmp*dt_inner/3
          ein_pulse = ein_pulse +eptmp*dt_inner/3
!.....3rd step
          call model_dte(tnow,dtep,eitmp,eotmp,eptmp,iprint)
          te(:,:,:) = te(:,:,:) +dtep(:,:,:)*dt_inner/3
          tep(:,:,:)= tep(:,:,:) +dtep(:,:,:)*dt_inner/2
          ein_e = ein_e +eitmp*dt_inner/3
          eout_e = eout_e +eotmp*dt_inner/3
          ein_pulse = ein_pulse +eptmp*dt_inner/3
!.....4th step
          call model_dte(tnow,dtep,eitmp,eotmp,eptmp,iprint)
          te(:,:,:) = te(:,:,:) +dtep(:,:,:)*dt_inner/6
          ein_e = ein_e +eitmp*dt_inner/6
          eout_e = eout_e +eotmp*dt_inner/6
          ein_pulse = ein_pulse +eptmp*dt_inner/6
        enddo
      else if( trim(csolver).eq.'RK2' ) then  ! 2th Runge-Kutta
        do istp=1,nstp_inner
          tep(:,:,:)= te(:,:,:)
!.....1st step
          call model_dte(tnow,dtep,eitmp,eotmp,eptmp,iprint)
          tep(:,:,:)= tep(:,:,:) +dtep(:,:,:)*dt_inner/2
!.....2nd step
          call model_dte(tnow,dtep,eitmp,eotmp,eptmp,iprint)
          te(:,:,:) = te(:,:,:) +dtep(:,:,:)*dt_inner
          ein_e = ein_e +eitmp*dt_inner
          eout_e = eout_e +eotmp*dt_inner
          ein_pulse = ein_pulse +eptmp*dt_inner
        enddo
      endif
      
!.....Update Te_max and alpha_max
      Te_max = -1d0
      alpha_max = -1d0
      cete_min = 1d30
      kappa_max = -1d0
      tep(:,:,:) = te(:,:,:)  ! need to use cete(...)
      do ic=1,nxyz
        call ic2ixyz(ic,ix,iy,iz)
        ce = cete(ix,iy,iz)
        if( itype_kappa.eq.1 ) then  ! DCrho
          kappa = d_e *ce *rho_e
        else if( itype_kappa.eq.2 ) then  ! B2
          kappa = kappa0 *te(ix,iy,iz) /max(ta(ic),ta_min)
        endif
        kappa_max = max(kappa_max,kappa)
        cete_min = min(ce,cete_min)
        alpha_max = max(kappa/ce/rho_e,alpha_max)
        Te_max = max(Te_max,te(ix,iy,iz))
        etot_e = etot_e +te(ix,iy,iz)*fkb*rho_e*vcell*1.5d0
      enddo

!.....Output
      if(iprint.gt.1) then
        if( ( itype_pulse.eq.2  &  ! gaussian
             .and. (tnow.ge.t0_laser .and. tnow.le.(t0_laser+tau_pulse*2)) ) &
             .or. ( itype_pulse.eq.1 & ! stepwise
             .and. (tnow.ge.t0_laser .and. tnow.le.t0_laser+tau_pulse) ) ) then
          print '(a,2es13.4e3)',' tnow,ein_pulse=',tnow,ein_pulse
        endif
        do ic=1,nxyz
          call ic2ixyz(ic,ix,iy,iz)
          if( te(ix,iy,iz).lt.0d0 ) then
            print *,'ERROR: te(ix,iy,iz) < 0 !!'
            print *,'ic,ix,iy,iz,te=',ic,ix,iy,iz,te(ix,iy,iz)
            stop 1
          endif
        enddo
      endif
    endif
!.....Broadcast Te distribution to all the nodes.
!.....There could be smarter way to reduce communication cost.
    call mpi_bcast(te,(nx+2)*(ny+2)*(nz+2),mpi_real8,0,mpi_world,ierr)

    t_ttm = t_ttm +mpi_wtime()-t0

    return
  end subroutine update_Te
!=======================================================================
  function cete(ix,iy,iz)
!
!  Ce(Te) at (ix,iy,iz)
!
    integer,intent(in):: ix,iy,iz

    real(8):: cete
    real(8):: t

    cete = 0d0
    if( iCe_Tdep.eq.0 ) then  ! none
      cete = c_0
    else if( iCe_Tdep.eq.1 ) then ! polynomial
      t = tep(ix,iy,iz)/1000
      cete = c_0 +(a_0 +a_1*t +a_2*t**2 +a_3*t**3 +a_4*t**4)&
           *exp(-(A_exp*t)**2) +ce_min
    else if( iCe_Tdep.eq.2 ) then ! tanh
      t = tep(ix,iy,iz)
      cete = 3d0 *tanh(2d-4 *t) +ce_min
    else if( iCe_Tdep.eq.3 ) then  ! linear
      cete = gmm_ce *tep(ix,iy,iz) +ce_min
    endif
    return
  end function cete
!=======================================================================
  function dcete(ix,iy,iz)
!
!  dCe(Te)/dTe at (ix,iy,iz)
!
    integer,intent(in):: ix,iy,iz

    real(8):: dcete
    real(8):: t,texp

    dcete = 0d0
    if( iCe_Tdep.eq.1 ) then  ! polynomial
      t = tep(ix,iy,iz)
      texp = exp(-(A_exp*t)**2)
      dcete = (a_1 +2d0*a_2*t +3d0*a_3*t**2 +4d0*a_4*t**3)*texp &
           -2d0*A_exp*t *(a_0 +a_1*t +a_2*t**2 +a_3*t**3 +a_4*t**4)*texp
    else if( iCe_Tdep.eq.2 ) then  ! tanh
      t = tep(ix,iy,iz)
      dcete = 3d0*2d-4 *(1d0 -tanh(2d-4 *t)**2)
    else if( iCe_Tdep.eq.3 ) then  ! linear
      dcete = gmm_ce
    endif
    return
    
  end function dcete
!=======================================================================
  function dte2(ix,iy,iz)
    integer,intent(in):: ix,iy,iz

    integer:: ixp,ixm,iyp,iym,izp,izm
    real(8):: dte2
    ixp = ix +1
    ixm = ix -1
    iyp = iy +1
    iym = iy -1
    izp = iz +1
    izm = iz -1

    dte2 = 0d0
    dte2 = (tep(ixp,iy,iz)**2 -2d0*tep(ixp,iy,iz)*tep(ixm,iy,iz) &
         +tep(ixm,iy,iz)**2 )/4/dx/dx &
         + (tep(ix,iyp,iz)**2 -2d0*tep(ix,iyp,iz)*tep(ix,iym,iz) &
         +tep(ix,iym,iz)**2 )/4/dy/dy &
         + (tep(ix,iy,izp)**2 -2d0*tep(ix,iy,izp)*tep(ix,iy,izm) &
         +tep(ix,iy,izm)**2 )/4/dz/dz
    return
  end function dte2
!=======================================================================
  function d2te(ix,iy,iz)
    integer,intent(in):: ix,iy,iz

    integer:: ixp,ixm,iyp,iym,izp,izm
    real(8):: d2te
    ixp = ix +1
    ixm = ix -1
    iyp = iy +1
    iym = iy -1
    izp = iz +1
    izm = iz -1

    d2te = (tep(ixp,iy,iz) -2d0*tep(ix,iy,iz) &
         +tep(ixm,iy,iz)) /dx/dx &
         + (tep(ix,iyp,iz) -2d0*tep(ix,iy,iz) &
         +tep(ix,iym,iz)) /dy/dy &
         + (tep(ix,iy,izp) -2d0*tep(ix,iy,iz) &
         +tep(ix,iy,izm)) /dz/dz
    return
  end function d2te
!=======================================================================
  subroutine langevin_ttm(namax,natm,va,aa,tag,am,h, &
       nspmax,fa2v,fekin,ediff,dtmd,myid,mpi_world,iprint)
!
!  Langevin thermostat for atomic system.
!
    use util,only: itotOf, ifmvOf
    include "params_unit.h"
    integer,intent(in):: namax,natm,nspmax,myid,mpi_world,iprint
    real(8),intent(in):: aa(3,namax),tag(namax),am(nspmax) &
         ,fa2v(nspmax),fekin(nspmax),dtmd,h(3,3)
    real(8),intent(inout):: va(3,namax),ediff(nspmax)

    integer:: ic,i,l,ifmv,ix,iy,iz,naccp,ierr,isp
    real(8):: hscl(3),sgmi,ami,ek,gmmi,vl(3),vi(3),aai(3),t0,vt(3)&
         ,aain(3),aaout(3),vin(3),vout(3),v0(3)
    real(8):: ediffl(nspmax),deinl(nspmax),deoutl(nspmax)
!!$    integer,external:: ifmvOf,itotOf
    real(8),external:: box_muller,sprod
    logical,save:: l1st = .true.

    if( l1st ) then
      allocate(dein(nspmax),deout(nspmax))
      l1st = .false.
    endif

    t0 = mpi_wtime()
    ein_a = 0d0
    eout_a = 0d0

    if( lcut_interact ) return

    do ic=1,nxyz
      call ic2ixyz(ic,ix,iy,iz)
      sgm(ic) = dsqrt(2d0*gmmp(ic)*te(ix,iy,iz)/dtmd *k2ue )
!!$      if( ic.lt.20 ) print *,'ic,sgm,Te,fkb*Te=',ic,sgm(ic)&
!!$           ,te(ix,iy,iz),fkb*te(ix,iy,iz)
    enddo

!.....Langevin thermostat with Mannella integrator
    ediffl(1:nspmax) = 0d0
    deinl(1:nspmax) = 0d0
    deoutl(1:nspmax) = 0d0
    hscl(1:3)= 0d0
    do l=1,3
      hscl(l)= dsqrt(h(1,l)**2 +h(2,l)**2 +h(3,l)**2)
    enddo
    do i=1,natm
      ifmv = ifmvOf(tag(i))
      isp = int(tag(i))
      if( ifmv.eq.0 ) then
        va(1:3,i)= 0d0
      else
        ic = a2c(i)
        call ic2ixyz(ic,ix,iy,iz)
        if( ix.lt.lsurf ) cycle
        vt(1:3) = va(1:3,i) -vac(1:3,ic)
        ami= am(isp)
        sgmi = sgm(ic) *dsqrt(ami)
        ek = ekti(i)
        gmmi = gmmp(ic)
        if( ek.gt.ekth ) gmmi = gmmp(ic) + gmms(ic)
        aai(1:3)= 0d0
        aain(1:3)= 0d0
        aaout(1:3)= 0d0
!.....SGMI should be [eV/Ang] whereas it is [ue/Ang]
!     and V0*GMMI*AMI is also [ue/Ang],
!     so need to multiply ue2ev
        do l=1,3
          aain(l) = sgmi*box_muller()/hscl(l) *ue2ev
          aaout(l) = -vt(l)*gmmi*ami *ue2ev
          aai(l) = aaout(l) +aain(l)
        enddo
!.....To compensate the factor 1/2 in fa2v, multiply 2 here.
        va(1:3,i)= va(1:3,i) +aai(1:3)*fa2v(isp)*dtmd *2d0
        if( va(1,i)*0d0.ne.0d0 .or. va(2,i)*0d0.ne.0d0 &
             .or. va(3,i)*0d0.ne.0d0 ) then
          if( myid.eq.0 ) then
            print *,'ERROR: va==NaN !!!'
            print *,'  ic,i,va(:)=',ic,i,va(1:3,i)
            print *,'  aain,aaout=',aain(1:3),aaout(1:3)
            print *,'  sgmi=',sgmi
            print *,'  gmmp(ic),te(ix,iy,iz)=',gmmp(ic),te(ix,iy,iz)
          endif
          stop
        endif
!.....accumulate energy difference
        vi(1:3)= h(1:3,1)*vt(1) &
             +h(1:3,2)*vt(2) &
             +h(1:3,3)*vt(3)
        vl(1:3)= h(1:3,1)*aai(1)*fa2v(isp)*dtmd *2d0 &
             +h(1:3,2)*aai(2)*fa2v(isp)*dtmd *2d0 &
             +h(1:3,3)*aai(3)*fa2v(isp)*dtmd *2d0
        vin(1:3)= h(1:3,1)*aain(1)*fa2v(isp)*dtmd *2d0 &
             +h(1:3,2)*aain(2)*fa2v(isp)*dtmd *2d0 &
             +h(1:3,3)*aain(3)*fa2v(isp)*dtmd *2d0
        vout(1:3)= h(1:3,1)*aaout(1)*fa2v(isp)*dtmd *2d0 &
             +h(1:3,2)*aaout(2)*fa2v(isp)*dtmd *2d0 &
             +h(1:3,3)*aaout(3)*fa2v(isp)*dtmd *2d0
        ediffl(isp)= ediffl(isp) +fekin(isp) &
             *(2d0*sprod(3,vi,vl)+sprod(3,vl,vl))
        deinl(isp)= deinl(isp) +fekin(isp) &
             *(2d0*sprod(3,vi,vin)+sprod(3,vin,vin))
        deoutl(isp)= deoutl(isp) +fekin(isp) &
             *(2d0*sprod(3,vi,vout)+sprod(3,vout,vout))
      endif
    enddo

    ediff(1:nspmax) = 0d0
    dein(1:nspmax) = 0d0
    deout(1:nspmax) = 0d0
    call mpi_reduce(ediffl,ediff,nspmax &
         ,mpi_real8,mpi_sum,0,mpi_world,ierr)
    call mpi_reduce(deinl,dein,nspmax &
         ,mpi_real8,mpi_sum,0,mpi_world,ierr)
    call mpi_reduce(deoutl,deout,nspmax &
         ,mpi_real8,mpi_sum,0,mpi_world,ierr)
    do isp=1,nspmax
      ein_a = ein_a +dein(isp)
      eout_a = eout_a +deout(isp)
    enddo

    t_ttm = t_ttm +mpi_wtime()-t0
    return
  end subroutine langevin_ttm
!=======================================================================
  subroutine ic2ixyz(ic,ix,iy,iz)
    integer,intent(in):: ic
    integer,intent(out):: ix,iy,iz

    iz = (ic-1) / (nx*ny) +1
    iy = mod(ic-1,nx*ny)/nx +1
    ix = mod(ic-1,nx) +1
    return
  end subroutine ic2ixyz
!=======================================================================
  subroutine ixyz2ic(ix,iy,iz,ic)
    integer,intent(in):: ix,iy,iz
    integer,intent(out):: ic

    ic = (iz-1)*nx*ny +(iy-1)*nx +ix
    return
  end subroutine ixyz2ic
!=======================================================================
  subroutine output_Te(istp,tnow,myid,iprint)
    integer,intent(in):: istp,myid,iprint
    real(8),intent(in):: tnow

    integer:: ix,iy,iz,ic,n
    real(8):: ave,eetot
    character(len=128):: cnum

    if( myid.eq.0 ) then
!!$!.....Average over y-z plane
!!$      tex(:) = 0d0
!!$      do ix=1,nx
!!$        do iy=1,ny
!!$          do iz=1,nz
!!$            tex(ix) = tex(ix) +te(ix,iy,iz)
!!$          enddo
!!$        enddo
!!$        tex(ix) = tex(ix) /(ny*nz)
!!$      enddo

!.....Output
      if( iprint.ne.0 ) then
        write(cnum,'(i0)') istp
        open(ioTeout,file=trim(cTe_outfile)//'_'//trim(cnum),status='replace')
        write(ioTeout,'(a)') '# ix,   iy,   iz,   te(ix,iy,iz),   ta(ic),'&
             //'  nac(ic)'
        do ix=1,nx
          do iy=1,ny
            do iz=1,nz
              call ixyz2ic(ix,iy,iz,ic)
              write(ioTeout,'(3i6,2es15.5,i6)') ix,iy,iz,te(ix,iy,iz) &
                   ,ta(ic),nac(ic)
            enddo
          enddo
        enddo
        close(ioTeout)
      endif

      if( iprint.gt.1 ) then
!.....Average Te
        eetot = 0d0
        ave = 0d0
        n = 0
        do ix=1,nx
          do iy=1,ny
            do iz=1,nz
              call ixyz2ic(ix,iy,iz,ic)
              if( nac(ic).eq.0 ) continue
              ave = ave +te(ix,iy,iz)
              eetot = eetot +3d0/2 *rho_e*vcell*te(ix,iy,iz)
              n = n + 1
!!$              write(ioTeout,'(3i6,2es15.5,i6)') ix,iy,iz,te(ix,iy,iz) &
!!$                   ,ta(ic),nac(ic)
            enddo
          enddo
        enddo
        ave = ave /n
        print *,'istp,tnow,te_ave,Ee_tot=',istp,tnow,ave,eetot*fkb
      endif
    endif
  end subroutine output_Te
!=======================================================================
  subroutine set_BC()
!
!  Set boundary condition
!
    
    if( lsurf.le.0 .and. rsurf.le.0 ) then
!.....Periodic for x,y,z
      tep(0,1:ny,1:nz) = tep(nx,1:ny,1:nz)
      tep(nx+1,1:ny,1:nz) = tep(1,1:ny,1:nz)
      tep(1:nx,0,1:nz) = tep(1:nx,ny,1:nz)
      tep(1:nx,ny+1,1:nz) = tep(1:nx,1,1:nz)
      tep(1:nx,1:ny,0) = tep(1:nx,1:ny,nz)
      tep(1:nx,1:ny,nz+1) = tep(1:nx,1:ny,1)
    else  ! Laser-ablation situation
!.....Free boundary for x
      tep(lsurf-1,1:ny,1:nz) = tep(lsurf,1:ny,1:nz)
      if( Te_right.lt.0d0 ) then
        tep(rsurf+1,1:ny,1:nz) = tep(rsurf,1:ny,1:nz)
      else
        tep(rsurf+1,1:ny,1:nz) = Te_right
      endif
!.....Periodic for y and z
      tep(lsurf:rsurf,0,1:nz)    = tep(lsurf:rsurf,ny,1:nz)
      tep(lsurf:rsurf,ny+1,1:nz) = tep(lsurf:rsurf,1,1:nz)
      tep(lsurf:rsurf,1:ny,0)    = tep(lsurf:rsurf,1:ny,nz)
      tep(lsurf:rsurf,1:ny,nz+1) = tep(lsurf:rsurf,1:ny,1)
    endif
  end subroutine set_BC
!=======================================================================
  subroutine update_surface_plane(myid,mpi_world,iprint)
!
!  Update left surface position according to number of atoms in the cells.
!  Criterion for lsurf:
!    - the average density over yz plane is less than a certain density (vacuum)
!    - right-most vacuum cell
!
    integer,intent(in):: myid,mpi_world,iprint

    integer:: ix,iy,iz,icl,ivac_right,ierr,lsurf_new,lsurf_true,imatt_right
    real(8):: tmp
    logical:: lupdate
    logical,allocatable,save:: lexists(:)
    real(8),allocatable,save:: densx(:)
    real(8),save:: volyz
    logical,save:: l1st = .true.

    if( l1st ) then
      allocate(lexists(nx),densx(nx))
      volyz = vcell*ny*nz
    endif

    if( myid.ne.0 ) goto 10

    lexists(:) = .false.
    densx(:) = 0d0
    do ix=nx,1,-1
      tmp = 0d0
      do iy=1,ny
        do iz=1,nz
          call ixyz2ic(ix,iy,iz,icl)
!.....Since nac is not a number of atoms, but DOF in a cell, so divide it by 3
          densx(ix) = densx(ix) +dble(nac(icl))/3
!.....Check true lsurf from given Te(:,:,:)
          tmp = tmp +te(ix,iy,iz)
        enddo
      enddo
      densx(ix) = densx(ix)/volyz
      if( densx(ix) .gt. dthresh ) lexists(ix) = .true.
!!$      print *,'ix,densx,lexists=',ix,densx(ix),lexists(ix)
      if( tmp.gt.0d0 ) lsurf_true = ix
    enddo

!.....Right-most layer of vacuum and right-most layer of material
    ivac_right = 0
    imatt_right = 0
    do ix=1,nx
      if( lexists(ix) ) imatt_right = ix
      if( .not. lexists(ix) ) ivac_right = ix
    enddo
!!$    print *,'ivac_right,imatt_right=',ivac_right,imatt_right

!.....Update lsurf
    lupdate = .false.
    lsurf_new = max(ivac_right+1, 1)  ! not to take 0
!!$    print *,'lsurf,new,true=',lsurf,lsurf_new,lsurf_true

!.....Update Te
    if( cTe_init(1:4).eq.'homo' .and. l1st ) then  ! the first call in case of homogeneous Te
      te(:,:,:) = 0d0
      do ix=lsurf_new,nx
        te(ix,:,:) = te_init
      enddo
    else if( lsurf_new .gt. lsurf_true ) then
      do ix=lsurf_true,lsurf_new-1
        te(ix,:,:) = 0d0
      enddo
    else if( lsurf_new .lt. lsurf_true ) then
!!$      print *,'lsurf_new.lt.lsurf'
      do ix=lsurf_new,lsurf_true-1
        te(ix,:,:) = te(lsurf_true,:,:)
      enddo
    endif
    lsurf = lsurf_new

    if( l1st ) then
      rsurf = imatt_right
!!$      do ix=1,nx
!!$        print *,'ix,te=',ix,te(ix,1,1)
!!$      enddo
    endif

10  continue
!.....Broadcast Te distribution to all the nodes.
!.....There could be smarter way to reduce networking cost.
    call mpi_bcast(te,(nx+2)*(ny+2)*(nz+2),mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(lsurf,1,mpi_integer,0,mpi_world,ierr)

    if( lsurf.ge.rsurf ) then
      if( myid.eq.0 ) print *,'Warning: lsurf.ge.rsurf, which should not happen!'
    endif
    if( myid.eq.0 .and. iprint.gt.1 ) then
      print '(a,2i4,5es12.4)', ' lsurf,ivac_right,densx= ' &
           ,lsurf,ivac_right,densx(max(1,lsurf-2):lsurf+2)
    endif
    l1st = .false.
    return
  end subroutine update_surface_plane
!=======================================================================
  subroutine output_energy_balance(istp,simtime,myid,iprint)
!
!  Output in/out of energy at electronic and atomic systems
!
    integer,intent(in):: istp, myid, iprint
    real(8),intent(in):: simtime
    logical:: lopen

    if( myid.eq.0 ) then
      inquire(unit=ioergio, opened=lopen)
      if( .not. lopen ) then
        open(ioergio,file=trim(cergio),status='replace')
        write(ioergio,'(a)') '# 1:istp, 2:time,        3:ein_el,' &
             //'       4:eout_el,        5:ein_at,      6:eout_at,' &
             //'       7:etot_el'
      endif
      write(ioergio,'(i8,6es15.7)') istp, simtime, ein_e, eout_e, ein_a &
           ,eout_a ,etot_e
      call flush(ioergio)
    endif
    
  end subroutine output_energy_balance
!=======================================================================
  subroutine remove_ablated_atoms(simtime,namax,natm,tag,ra,va,&
       chg,chi,h,sorg)
!
!  Remove atoms evapolated from left surface by laser ablation.
!  Store the following removed-atom data:
!    - ID, time, position (y,z), velocity (x,y,z)
!
    use pmdmpi
    use util,only: itotOf
    integer,intent(in):: namax
    integer,intent(inout):: natm
    real(8),intent(inout):: tag(namax),ra(3,namax),va(3,namax)&
         ,chg(namax),chi(namax)
    real(8),intent(in):: h(3,3),simtime,sorg(3)
    
    integer:: ix,iy,iz,inc,nabl,n,ia,ja,inc2,ierr,itot,iyz,n0
    real(8):: r(3),rt(3),v(3)
    logical:: lremove 
    integer,save:: ntmp = 1000
    integer,parameter:: nmpi = 4
    logical,save:: l1st = .true.
    integer,allocatable,save:: list(:)
    real(8),allocatable,save:: tagabl(:),rabl(:,:),vabl(:,:)
    integer:: istat(mpi_status_size),itag

!!$    integer,external:: itotOf

    if( l1st ) then
      if( myid_md.eq.0 ) then
        print '(a)','# ABLATED: ID, time, position (y,z), velocity (x,y,z)'
      endif
      allocate(tagabl(ntmp),rabl(3,ntmp),vabl(3,ntmp),list(ntmp))
      l1st = .false.
    endif

    call nid2xyz(myid_md,ix,iy,iz)

!.....Only left-most node is concerned
    if( ix.eq.0 ) then
      inc = 0
      do ia=1,natm
        if( ra(1,ia).lt.0d0 ) then
          inc = inc + 1
        endif
      enddo
      nabl = inc
      if( myid_md.eq.0 ) then
        do iyz = 1,ny*nz-1 ! node-ID of ix = 0
          itag = iyz
          call mpi_recv(n,1,mpi_integer,iyz,itag,mpi_md_world,istat,ierr)
          nabl = nabl + n
        enddo
      else  ! myid_md.ne.0
        itag = myid_md
        call mpi_send(inc,1,mpi_integer,0,itag,mpi_md_world,ierr)
      endif

      if( nabl.gt.ntmp ) then
        ntmp = nabl *1.5
        deallocate(tagabl,rabl,vabl,list)
        allocate(tagabl(ntmp),rabl(3,ntmp),vabl(3,ntmp),list(ntmp))
      endif

!.....Actually list up to-be-removed atoms
      inc = 0
      list(:) = 0
      do ia=1,natm
        if( ra(1,ia).lt.0d0 ) then
          inc = inc + 1
          list(inc) = ia
          tagabl(inc) = tag(ia)
          rt(1:3) = ra(1:3,ia) +sorg(1:3)
          r(1:3) = h(1:3,1)*rt(1) +h(1:3,2)*rt(2) +h(1:3,3)*rt(3)
          v(1:3) = h(1:3,1)*va(1,ia) +h(1:3,2)*va(2,ia) +h(1:3,3)*va(3,ia)
          rabl(1:3,inc) = r(1:3)
          vabl(1:3,inc) = v(1:3)
        endif
      enddo
!.....Gather to-be-removed atoms for writing out
      n0 = inc + 1
      if( myid_md.eq.0 ) then
        do iyz = 1,ny*nz-1
          itag = iyz*nmpi -nmpi
          call mpi_recv(n,1,mpi_integer,iyz,itag,mpi_md_world,istat,ierr)
          call mpi_recv(tagabl(n0),n,mpi_real8,iyz,itag+1,mpi_md_world,istat,ierr)
          call mpi_recv(rabl(1,n0),3*n,mpi_real8,iyz,itag+2,mpi_md_world,istat,ierr)
          call mpi_recv(vabl(1,n0),3*n,mpi_real8,iyz,itag+3,mpi_md_world,istat,ierr)
          n0 = n0 + n
        enddo
      else  ! myid_md.ne.0
        itag = myid_md*nmpi -nmpi
        call mpi_send(inc,1,mpi_integer,0,itag,mpi_md_world,ierr)
        call mpi_send(tagabl,inc,mpi_real8,0,itag+1,mpi_md_world,ierr)
        call mpi_send(rabl,3*inc,mpi_real8,0,itag+2,mpi_md_world,ierr)
        call mpi_send(vabl,3*inc,mpi_real8,0,itag+3,mpi_md_world,ierr)
      endif
!.....Write out only at node-0
      if( myid_md.eq.0 ) then
        do ia = 1,nabl
          itot = itotOf(tagabl(ia))
          print '(a,i10,6es13.5)','ABLATED: ',itot,simtime,rabl(2:3,ia),vabl(1:3,ia)
        enddo
      endif
!.....Remove atoms from the ra,va,tag,...
      inc2 = 0
      do ia=1,natm
        lremove = .false.
        do ja=1,inc
          if( ia.eq.list(ja) ) then
            lremove = .true.
!!$            print *,' remove itot,myid=',itotOf(tag(ia)),myid_md
            exit
          endif
        enddo
        if( .not. lremove ) then
          inc2 = inc2 + 1
          ra(1:3,inc2) = ra(1:3,ia)
          va(1:3,inc2) = va(1:3,ia)
          tag(inc2) = tag(ia)
          chg(inc2) = tag(ia)
          chi(inc2) = tag(ia)
        endif
      enddo
      natm = inc2
    endif
    
  end subroutine remove_ablated_atoms
!=======================================================================
  subroutine te2tei(namax,natm,tei)
!
!  Get electronic temperatures of atoms from electronic temperatures on sites.
!
    integer,intent(in):: namax,natm
    real(8),intent(out):: tei(namax)

    integer:: ic,ia,ix,iy,iz

    do ia=1,natm
      ic = a2c(ia)
      call ic2ixyz(ic,ix,iy,iz)
      if( ix.lt.lsurf ) cycle
      tei(ia) = te(ix,iy,iz)
    enddo
    return
  end subroutine te2tei
end module ttm
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
