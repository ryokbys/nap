module ttm
!-----------------------------------------------------------------------
!                     Last-modified: <2021-02-01 17:51:38 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
!
! Module for two(or three?)-temperature method (TTM).
!
! In the current implementation, it is assumed that the number of 
! parallel nodes are the common dividors of number of TTM meshes.
!
  implicit none
  save
  include 'mpif.h'
  include "./params_unit.h"
  
  character(len=128),parameter:: cfparams = 'in.ttm'
  character(len=128),parameter:: cin_ts3d = 'in.Ts3d'
  character(len=128),parameter:: cin_ts1d = 'in.Ts1d'
  character(len=128),parameter:: cTe_outfile = 'out.Te'
  character(len=128),parameter:: cergio = 'out.eio_ttm'
  character(len=128),parameter:: dname_ttm = 'outs_ttm/'
  character(len=128),parameter:: cout_ts3d = 'out.Ts3d'
  character(len=128),parameter:: cout_ts1d = 'out.Ts1d'
  
  integer,parameter:: ioprms = 30
  integer,parameter:: ioTein = 31
  integer,parameter:: ioTeout = 32
  integer,parameter:: ioergio = 33
  integer,parameter:: iots3d = 34
  integer,parameter:: iots1d = 35

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
!.....Coefficients for polynomial Ce(Te) [eV/(K.electron)]
!.....See Jay et al., IEEE Trans. Nucl. Sci. 64 (2017)
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

!.....1D continume TTM and boundary condition parameters
  integer:: nd1d = 1000  ! Num of nodes for 1D TTM
  real(8):: dx1d = 10d0  ! dx (Ang) of 1D TTM
  real(8):: cl1d = 1d-4 ! specific heat of lattice system [eV/K/atom]
  integer:: ibc1d,ibc3d
  real(8),allocatable:: te1d(:),tep1d(:),tl1d(:),tlp1d(:), &
       gp1d(:),gmmp1d(:)
  real(8):: rho_bulk
!.....Non-reflecting boundary condition
  real(8):: dnr = 10d0   ! NRBC region length [default: 10 Ang]
  real(8):: xrmd
  real(8):: ssound = 8433d-5  ! speed of sound [default: 8433d-5 Ang/fs] for Si

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
    real(8):: t,t0,t1,t2,dtmax,tmp
    character(len=128):: c1st

    t_ttm = 0d0
    t0 = mpi_wtime()
!.....Read parameter file
    call read_ttm_params(myid,mpi_world,iprint)
    call sync_params(myid,mpi_world,iprint)

!.....Make directory for storing TTM temperature files
    if( myid.eq.0 ) call system('mkdir -p '//trim(dname_ttm))

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
    if( trim(csolver).eq.'RK4' ) then
      dr2 = dr2*4
    else if( trim(csolver).eq.'RK2') then
      dr2 = dr2*2
    endif

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
        print *,'ERROR@init_ttm: pulse_type should be either stepwise or gaussian.'
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
      print '(a,2i5)','   lsurf = ',lsurf
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
        print *,'ERROR @init_ttm: wrong nxyz: ',nxyz
        print *,'  nx,ny,nz = ',nx,ny,nz
      endif
      goto 999
    else if( nd1d*dx1d.le.h(1,1) ) then
      if( myid.eq.0 ) then
        print *,'ERROR @init_ttm: nd1d*dx1d < h(1,1), which means 1D-TTM system is '&
             //'shorter than 3D-TTM, which should not happen.'
        print *,'  Increase num_node_ttm1d and/or dx_ttm1d.'
      endif
      goto 999
    else if( dx1d.lt.dx ) then
      if( myid.eq.0 ) then
        print *,'ERROR @init_ttm: dx1d < dx, which means 1D-TTM mesh size is '&
             //'smaller than that of 3D-TTM, which is not available.'
        print '(a,es12.4)','   Set dx_ttm1d greater than ',dx
      endif
      goto 999
    else if( cl1d.lt.0d0 ) then
      if( myid.eq.0 ) then
        print *,'ERROR @init_ttm: Cl_ttm1d is not set. Need to set a positive value.'
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
!.....1D TTM related arrays
    allocate(te1d(nd1d+1),tep1d(nd1d+1), &
         tl1d(nd1d+1),tlp1d(nd1d+1), &
         gp1d(nd1d),gmmp1d(nd1d))

    if( trim(ctype_coupling).eq.'constant_gmmp' ) then
      itype_coupling = 1
      gmmp(:) = gamma_p
      gmms(:) = gamma_s
      gmmp1d(:) = gamma_p
    else if( trim(ctype_coupling).eq.'constant_gp' ) then
      itype_coupling = 2
      gp(:) = e_ph_const
      gs(:) = 0d0
      gp1d(:) = e_ph_const
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
    else if( cTe_init(1:4).eq.'homo' ) then  ! homogeneous Te
      if( myid.eq.0 ) then
        print '(a,f7.1,a)','   Initial Te =     ',Te_init,' K'
      endif
      te(:,:,:) = Te_init
      te1d(:) = Te_init
      tl1d(:) = Te_init  ! Initial Temp of 1D-TTM lattice sys == electron sys
    else if( trim(cTe_init).eq.'read' ) then
      if( myid.eq.0 ) then
!  Initial temperatures are read from in.ts3d (Te in TTM3D region).
!  The format of in.ts3d should be like the following.
!-----------------------------------------------------------------------
!  #    ix,  iy,  iz,  temp
!        1    1    1    300.0
!        1    1    2    300.0
!  ...
!-----------------------------------------------------------------------
        print '(a)','   Initial Te for 3D TTM-MMD region is read from in.ts3d'
        open(iots3d,file=trim(cin_ts3d),status='old')
        do while(.true.)
          read(iots3d,*,end=10) c1st
          if( c1st(1:1).eq.'!' .or. c1st(1:1).eq.'#' ) cycle
          backspace(iots3d)
          read(iots3d,*) ix,iy,iz,t
          te(ix,iy,iz) = t
        enddo
10      close(iots3d)
!-----------------------------------------------------------------------
!  Initial temperatures in 1D TTM region are read from in.ts1d
!  The format of in.ts1d should be like the following.
!-----------------------------------------------------------------------
!  #  dx: ????
!  #    ix,  te,  tl
!        1   300.0  301.0
!        2   291.3  298.0
!  ...
!-----------------------------------------------------------------------
        print '(a)','   Initial Te and Tl in 1D TTM region are read from in.ts1d'
        open(iots1d,file=trim(cin_ts1d),status='old')
        do while(.true.)
          read(iots1d,*,end=20) c1st
          if( c1st(1:1).eq.'!' .or. c1st(1:1).eq.'#' ) cycle
          backspace(iots1d)
          read(iots1d,*) ix, t1, t2
          te1d(ix) = t1
          tl1d(ix) = t2
        enddo
20      close(iots1d)
      endif
!.....Broadcast Te distribution to all the nodes.
!.....There could be smarter way to reduce networking cost.
      call mpi_bcast(te,(nx+2)*(ny+2)*(nz+2),mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(te1d,nd1d+1,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(tl1d,nd1d+1,mpi_real8,0,mpi_world,ierr)
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
      dt_inner = dr2 /(2d0*alpha_max)
!!$      print *,'dr2,alpha_max,dt_inner,dtmd=',dr2,alpha_max,dt_inner,dtmd
    endif

!.....If dt_inner is specified, nstp_inner and dt are determined from dt_inner and dtmd.
!.....Since the upper limit of dt is given by dx**2/(2*kappa),
!.....dt_inner should be determined regardless to dtmd and
!.....this would be appropriate especially in case of variable time-step.
    if( dt_inner .gt. 0d0 ) then
!!$      nstp_inner = max(int(dtmd /dt_inner),1)
      nstp_inner = int(dtmd /dt_inner) +1
    endif
!.....Change dt_inner to match dt and use this dt_inner afterward
    dt_inner = dtmd /nstp_inner

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
      fname = trim(cfparams)
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
          print *,'WARNING: right_surface is deprecated. Not used in this version.'
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
        else if( trim(c1st).eq.'num_node_ttm1d' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, nd1d
        else if( trim(c1st).eq.'dx_ttm1d' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, dx1d
        else if( trim(c1st).eq.'Cl_ttm1d' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, cl1d
        else if( trim(c1st).eq.'NRBC_length' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, dnr
        else if( trim(c1st).eq.'speed_sound' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, ssound
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
    call mpi_bcast(rsurf,1,mpi_integer,0,mpi_world,ierr)  ! deprecated
    call mpi_bcast(lskin,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(lcut_interact,1,mpi_logical,0,mpi_world,ierr)
    call mpi_bcast(csolver,20,mpi_character,0,mpi_world,ierr)
    call mpi_bcast(nd1d,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(dx1d,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(dnr,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(ssound,1,mpi_real8,0,mpi_world,ierr)
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
        rho_bulk = 0d0
        do ic=1,nxyz
          if( nac(ic).eq.0 ) cycle
!!$          call ic2ixyz(ic,ix,iy,iz)
!.....Degree of freedom per atom (3 in case of 3D) is included in nac
!.....CHECK: This factor 3 looks causing the difference of energy in/out between at/el systems.
          ta(ic) = eksum(ic) *2d0 /fkb /nac(ic)
          gp(ic) = nac(ic) *fkb *gmmp(ic) /vcell ! /3
          rho_bulk = max(dble(nac(ic))/3/vcell,rho_bulk)
          if( nacp(ic).eq.0 ) cycle
          tap(ic) = ekpsum(ic) *2d0 /fkb /nacp(ic)
          gs(ic) = nacp(ic) *fkb *gmms(ic) /vcell ! /3
        enddo
        gp1d(:) = fkb *gmmp1d(:)*(rho_bulk*3)
      else if( trim(ctype_coupling).eq.'constant_gp' ) then
!.....See Eq.(A5) in PRB 68 (2003) pp.064114
        gmmp(:) = 0d0
        gmms(:) = 0d0
        ta(:) = 0d0
        tap(:) = 0d0
        rho_bulk = 0d0
        do ic=1,nxyz
          if( nac(ic).eq.0 ) cycle
          call ic2ixyz(ic,ix,iy,iz)
          ta(ic) = eksum(ic) *2d0 /fkb /nac(ic)
          rho_bulk = max(dble(nac(ic))/3/vcell,rho_bulk)
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
        gmmp1d(:) = gp1d(:)/fkb /(rho_bulk*3)
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
  subroutine update_ttm(tnow,dtmd,natm,ra,h,sorg,myid,mpi_world,iprint)
!
!  Wrapper routine for updating 3D-TTM and 1D-TTM systems.
!
    integer,intent(in):: myid,mpi_world,iprint,natm
    real(8),intent(in):: tnow,dtmd,ra(3,natm),h(3,3),sorg(3)

    call set_3d1d_bc_pos(natm,ra,h,sorg,myid,mpi_world,iprint)
    call couple_3d1d(myid,mpi_world,iprint)
    call update_2tm3d(tnow,dtmd,myid,mpi_world,iprint)
    call update_2tm1d(tnow,myid,mpi_world,iprint)
    
  end subroutine update_ttm
!=======================================================================
  subroutine update_2tm3d(tnow,dtmd,myid,mpi_world,iprint)
!
!  Update Te by solving the diffusion equation.
!  Model calculation is separated as model_XXX routines,
!  and some ODE solvers are (to be) implemented.
!
    integer,intent(in):: myid,mpi_world,iprint
    real(8),intent(in):: tnow,dtmd

    integer:: ic,ix,iy,iz,ierr,istp,ix0,ix1
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

      ix0 = 0
      ix1 = ibc3d

      if( trim(csolver).eq.'Euler' ) then
        do istp = 1,nstp_inner
          tep(ix0:ix1,:,:) = te(ix0:ix1,:,:)
          call model_2tm3d(tnow,dtep,eitmp,eotmp,eptmp,iprint)
          te(ix0:ix1,:,:) = te(ix0:ix1,:,:) +dtep(ix0:ix1,:,:)*dt_inner
          ein_e = ein_e +eitmp*dt_inner
          eout_e = eout_e +eotmp*dt_inner
          ein_pulse = ein_pulse +eptmp*dt_inner
        enddo  ! istp=1,nstp_inner
      else if( trim(csolver).eq.'RK4' ) then  ! 4th Runge-Kutta
        do istp=1,nstp_inner
          tep(:,:,:)= te(:,:,:)
!.....1st step
          call model_2tm3d(tnow,dtep,eitmp,eotmp,eptmp,iprint)
          te(ix0:ix1,:,:) = te(ix0:ix1,:,:) +dtep(ix0:ix1,:,:)*dt_inner/6
          tep(ix0:ix1,:,:)= tep(ix0:ix1,:,:) +dtep(ix0:ix1,:,:)*dt_inner/2
          ein_e = ein_e +eitmp*dt_inner/6
          eout_e = eout_e +eotmp*dt_inner/6
          ein_pulse = ein_pulse +eptmp*dt_inner/6
!.....2nd step
          call model_2tm3d(tnow,dtep,eitmp,eotmp,eptmp,iprint)
          te(ix0:ix1,:,:) = te(ix0:ix1,:,:) +dtep(ix0:ix1,:,:)*dt_inner/3
          tep(ix0:ix1,:,:)= tep(ix0:ix1,:,:) +dtep(ix0:ix1,:,:)*dt_inner/2
          ein_e = ein_e +eitmp*dt_inner/3
          eout_e = eout_e +eotmp*dt_inner/3
          ein_pulse = ein_pulse +eptmp*dt_inner/3
!.....3rd step
          call model_2tm3d(tnow,dtep,eitmp,eotmp,eptmp,iprint)
          te(ix0:ix1,:,:) = te(ix0:ix1,:,:) +dtep(ix0:ix1,:,:)*dt_inner/3
          tep(ix0:ix1,:,:)= tep(ix0:ix1,:,:) +dtep(ix0:ix1,:,:)*dt_inner/2
          ein_e = ein_e +eitmp*dt_inner/3
          eout_e = eout_e +eotmp*dt_inner/3
          ein_pulse = ein_pulse +eptmp*dt_inner/3
!.....4th step
          call model_2tm3d(tnow,dtep,eitmp,eotmp,eptmp,iprint)
          te(ix0:ix1,:,:) = te(ix0:ix1,:,:) +dtep(ix0:ix1,:,:)*dt_inner/6
          ein_e = ein_e +eitmp*dt_inner/6
          eout_e = eout_e +eotmp*dt_inner/6
          ein_pulse = ein_pulse +eptmp*dt_inner/6
        enddo
      else if( trim(csolver).eq.'RK2' ) then  ! 2th Runge-Kutta
        do istp=1,nstp_inner
          tep(ix0:ix1,:,:)= te(ix0:ix1,:,:)
!.....1st step
          call model_2tm3d(tnow,dtep,eitmp,eotmp,eptmp,iprint)
          tep(ix0:ix1,:,:)= tep(ix0:ix1,:,:) +dtep(ix0:ix1,:,:)*dt_inner/2
!.....2nd step
          call model_2tm3d(tnow,dtep,eitmp,eotmp,eptmp,iprint)
          te(ix0:ix1,:,:) = te(ix0:ix1,:,:) +dtep(ix0:ix1,:,:)*dt_inner
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
      tep(ix0:ix1,:,:) = te(ix0:ix1,:,:)  ! need to use cete(...)
      do ic=1,nxyz
        call ic2ixyz(ic,ix,iy,iz)
        if( ix.gt.ibc3d ) cycle
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
          if( ix.gt.ibc3d ) cycle
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
  end subroutine update_2tm3d
!=======================================================================
  subroutine model_2tm3d(tnow,dtep,eitmp,eotmp,eptmp,iprint)
!
!  Create f(t,y) for ODE dy/dt = f(t,y), where y=Te(ix,iy,iz) here.
!  The ODE is now diffusion Eq. with two-temperature model.
!
    real(8),intent(in):: tnow
    integer,intent(in):: iprint
    real(8),intent(out):: dtep(0:nx+1,0:ny+1,0:nz+1),eitmp,eotmp,eptmp

    integer:: ic,ix,iy,iz
    real(8):: ce,dce,kappa,dkappa,pterm,sterm,dtemp,de,tmp&
         ,pulsefactor,xi,denom
    
    dtep(:,:,:) = 0d0
    eitmp = 0d0
    eotmp = 0d0
    eptmp = 0d0
    call set_bc_2tm3d()
    do ic=1,nxyz
      call ic2ixyz(ic,ix,iy,iz)
      if( ix.lt.lsurf ) cycle
      if( ix.gt.ibc3d ) cycle
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
      denom = (ce+tep(ix,iy,iz)*dce)*rho_e
!!$      if( ix.eq.ibc3d ) then
!!$        print *,'3D: denom,gp,dT,pterm= ',denom,gp(ic),(tep(ix,iy,iz)-ta(ic)),pterm
!!$        print *,'3D: 1st,2nd,3rd,4th= ',dkappa*dte2(ix,iy,iz)/denom &
!!$             ,kappa*d2te(ix,iy,iz)/denom,pterm/denom,sterm/denom
!!$      endif
      dtemp = ( dkappa*dte2(ix,iy,iz) +kappa*d2te(ix,iy,iz) &
           +pterm +sterm ) /denom ! *dt
      dtep(ix,iy,iz) = dtep(ix,iy,iz) +dtemp
    enddo  ! ic=1,nxyz
!!$    if( ctype_pulse(1:4).eq.'step' ) then
    if( itype_pulse.eq.1 ) then  ! step-wise pulse
      if( tnow.ge.t0_laser .and. &
           tnow.le.(t0_laser +tau_pulse) ) then
        do ic=1,nxyz
          call ic2ixyz(ic,ix,iy,iz)
          if( ix.lt.lsurf ) cycle
          if( ix.gt.ibc3d ) cycle
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
          if( ix.gt.ibc3d ) cycle
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

  end subroutine model_2tm3d
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
    real(8):: t,texp,x

    dcete = 0d0
    if( iCe_Tdep.eq.1 ) then  ! polynomial
      x = tep(ix,iy,iz)/1000
      texp = exp(-(A_exp*x)**2)
      dcete = (a_1 +2d0*a_2*x +3d0*a_3*x**2 +4d0*a_4*x**3)*texp &
           -2d0*A_exp*x *(a_0 +a_1*x +a_2*x**2 +a_3*x**3 +a_4*x**4)*texp
      dcete = dcete/1000
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
  subroutine non_reflecting_bc(natm,tag,ra,va,h,sorg,dtmd &
       ,nspmax,am,fa2v,myid,mpi_world,iprint)
!
!  Langevin Non-reflecting_bc for atomic system.
!  See Shugaev, et al., PRB96 (2017) for details.
!
    integer,intent(in):: natm,myid,mpi_world,iprint,nspmax
    real(8),intent(in):: ra(3,natm),h(3,3),sorg(3),tag(natm) &
         ,dtmd,am(nspmax),fa2v(nspmax)
    real(8),intent(inout):: va(3,natm)

    integer:: i,nabcl,nabc,ierr,isp,ic
    real(8):: xdnr,area,xi,hxi,areatom,vx,ami,zimp,sgmi,axi
    real(8),external:: box_muller
    
    xdnr = dnr/h(1,1)
    area = h(2,2)*h(3,3)
!.....1st, count atoms to which LNRBC is applied
    nabcl = 0
    do i=1,natm
      xi = ra(1,i) +sorg(1)
      if( xi.lt.xrmd-xdnr ) cycle
      nabcl = nabcl +1
    enddo
    call mpi_allreduce(nabcl,nabc,1,mpi_integer,mpi_sum,mpi_world,ierr)
    areatom = area/nabc

    hxi = 1d0/h(1,1)
    do i=1,natm
      xi = ra(1,i) +sorg(1)
      if( xi.lt.xrmd-xdnr ) cycle
      vx = va(1,i)
      isp = int(tag(i))
      ami = am(isp)
      zimp = ami *rho_bulk *ssound
      ic = a2c(i)
      sgmi = dsqrt(2d0*zimp*areatom*ta(ic)/dtmd*k2ue)
      axi = (-vx*zimp*areatom +sgmi*box_muller()*hxi)*ue2ev
!.....To compensate the factor 1/2 in fa2v, multiply 2 here.
      va(1,i) = va(1,i) +axi*fa2v(isp)*dtmd *2d0
    enddo
    
  end subroutine non_reflecting_bc
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
  subroutine output_ttm(istp,tnow,myid,iprint)
!
!  Wrapper function for output 3D-TTM and 1D-TTM systems.
!
    integer,intent(in):: istp,myid,iprint
    real(8),intent(in):: tnow
    
    call output_ttm3d(istp,tnow,myid,iprint)
    call output_ttm1d(istp,tnow,myid,iprint)
  end subroutine output_ttm
!=======================================================================
  subroutine output_ttm3d(istp,tnow,myid,iprint)
!
!  Output Te data in 3D TTM-MD (ttm3d) region.
!
    integer,intent(in):: istp,myid,iprint
    real(8),intent(in):: tnow

    integer:: ix,iy,iz,ic,n
    real(8):: ave,eetot
    character(len=128):: cnum

    if( myid.eq.0 ) then

!.....Output
      if( iprint.ne.0 ) then
        write(cnum,'(i0)') istp
        open(iots3d,file=trim(dname_ttm)// &
             trim(cout_ts3d)//'_'//trim(cnum),status='replace')
        write(iots3d,'(a,2es15.7,i6)') '# tnow,dx,ibc3d: ',tnow,dx,ibc3d
        write(iots3d,'(a)') '# ix,   iy,   iz,   te(ix,iy,iz),   ta(ic),'&
             //'  nac(ic)'
        do ix=1,nx
          do iy=1,ny
            do iz=1,nz
              call ixyz2ic(ix,iy,iz,ic)
              write(iots3d,'(3i6,2es15.5,i6)') ix,iy,iz,te(ix,iy,iz) &
                   ,ta(ic),nac(ic)
            enddo
          enddo
        enddo
        close(iots3d)
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
  end subroutine output_ttm3d
!=======================================================================
  subroutine set_bc_2tm3d()
!
!  Set BC for 3D-TTM system
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
!!$      if( Te_right.lt.0d0 ) then
!!$        tep(rsurf+1,1:ny,1:nz) = tep(rsurf,1:ny,1:nz)
!!$      else
!!$        tep(rsurf+1,1:ny,1:nz) = Te_right
!!$      endif
!.....Periodic for y and z
      tep(lsurf:nx,0,1:nz)    = tep(lsurf:nx,ny,1:nz)
      tep(lsurf:nx,ny+1,1:nz) = tep(lsurf:nx,1,1:nz)
      tep(lsurf:nx,1:ny,0)    = tep(lsurf:nx,1:ny,nz)
      tep(lsurf:nx,1:ny,nz+1) = tep(lsurf:nx,1:ny,1)
    endif
  end subroutine set_bc_2tm3d
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

!!$    if( l1st ) then
!!$      rsurf = imatt_right
!!$    endif

10  continue
!.....Broadcast Te distribution to all the nodes.
!.....There could be smarter way to reduce networking cost.
    call mpi_bcast(te,(nx+2)*(ny+2)*(nz+2),mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(lsurf,1,mpi_integer,0,mpi_world,ierr)

    if( lsurf.ge.ibc3d ) then
      if( myid.eq.0 ) print *,'Warning: lsurf.ge.ibc3d, which should not happen!'
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
!=======================================================================
  subroutine update_2tm1d(tnow,myid,mpi_world,iprint)
!
!  Solve 1D continume TTM attached to right-most end of TTM-MD region.
!
    integer,intent(in):: myid,mpi_world,iprint
    real(8),intent(in):: tnow
    
    integer:: istp,ix0,ix1
    real(8),allocatable,save:: dtep1d(:),dtlp1d(:)
    logical,save:: l1st = .true.

    if( l1st ) then
      allocate(dtep1d(nd1d),dtlp1d(nd1d))
      l1st = .false.
    endif

!.....Only node-0 solve 1D-TTM 
    if( myid.eq.0 ) then
      ix0 = ibc1d+1
      ix1 = nd1d

      if( trim(csolver).eq.'Euler' ) then
        do istp=1,nstp_inner
          tep1d(:) = te1d(:)
          tlp1d(:) = tl1d(:)
          call model_2tm1d(tnow,dtep1d,dtlp1d,iprint)
          te1d(ix0:ix1) = te1d(ix0:ix1) +dtep1d(ix0:ix1)*dt_inner
          tl1d(ix0:ix1) = tl1d(ix0:ix1) +dtlp1d(ix0:ix1)*dt_inner
        enddo
      else if( trim(csolver).eq.'RK4' ) then
!!$        print *,'nstp_inner=',nstp_inner
        do istp=1,nstp_inner
          tep1d(:) = te1d(:)
          tlp1d(:) = tl1d(:)
!!$          print '(a,10es11.3)','0,te1d=',te1d(ix0-1:ix0+5)
!!$          print '(a,10es11.3)','0,tl1d=',tl1d(ix0-1:ix0+5)
!.....1st step
          call model_2tm1d(tnow,dtep1d,dtlp1d,iprint)
          te1d(ix0:ix1) = te1d(ix0:ix1) +dtep1d(ix0:ix1)*dt_inner/6
          tl1d(ix0:ix1) = tl1d(ix0:ix1) +dtlp1d(ix0:ix1)*dt_inner/6
          tep1d(ix0:ix1) = tep1d(ix0:ix1) +dtep1d(ix0:ix1)*dt_inner/2
          tlp1d(ix0:ix1) = tlp1d(ix0:ix1) +dtlp1d(ix0:ix1)*dt_inner/2
!!$          print '(a,10es11.3)','1,te1d=',te1d(ix0-1:ix0+5)
!!$          print '(a,10es11.3)','1,tl1d=',tl1d(ix0-1:ix0+5)
!.....2nd step
          call model_2tm1d(tnow,dtep1d,dtlp1d,iprint)
          te1d(ix0:ix1) = te1d(ix0:ix1) +dtep1d(ix0:ix1)*dt_inner/3
          tl1d(ix0:ix1) = tl1d(ix0:ix1) +dtlp1d(ix0:ix1)*dt_inner/3
          tep1d(ix0:ix1) = tep1d(ix0:ix1) +dtep1d(ix0:ix1)*dt_inner/2
          tlp1d(ix0:ix1) = tlp1d(ix0:ix1) +dtlp1d(ix0:ix1)*dt_inner/2
!!$          print '(a,10es11.3)','2,te1d=',te1d(ix0-1:ix0+5)
!!$          print '(a,10es11.3)','2,tl1d=',tl1d(ix0-1:ix0+5)
!.....3rd step
          call model_2tm1d(tnow,dtep1d,dtlp1d,iprint)
          te1d(ix0:ix1) = te1d(ix0:ix1) +dtep1d(ix0:ix1)*dt_inner/3
          tl1d(ix0:ix1) = tl1d(ix0:ix1) +dtlp1d(ix0:ix1)*dt_inner/3
          tep1d(ix0:ix1) = tep1d(ix0:ix1) +dtep1d(ix0:ix1)*dt_inner/2
          tlp1d(ix0:ix1) = tlp1d(ix0:ix1) +dtlp1d(ix0:ix1)*dt_inner/2
!!$          print '(a,10es11.3)','3,te1d=',te1d(ix0-1:ix0+5)
!!$          print '(a,10es11.3)','3,tl1d=',tl1d(ix0-1:ix0+5)
!.....4th step
          call model_2tm1d(tnow,dtep1d,dtlp1d,iprint)
          te1d(ix0:ix1) = te1d(ix0:ix1) +dtep1d(ix0:ix1)*dt_inner/6
          tl1d(ix0:ix1) = tl1d(ix0:ix1) +dtlp1d(ix0:ix1)*dt_inner/6
!!$          print '(a,10es11.3)','4,te1d=',te1d(ix0-1:ix0+5)
!!$          print '(a,10es11.3)','4,tl1d=',tl1d(ix0-1:ix0+5)
        enddo
      else if( trim(csolver).eq.'RK2' ) then
        do istp=1,nstp_inner
          tep1d(:) = te1d(:)
          tlp1d(:) = tl1d(:)
!.....1st step
          call model_2tm1d(tnow,dtep1d,dtlp1d,iprint)
          tep1d(ix0:ix1) = tep1d(ix0:ix1) +dtep1d(ix0:ix1)*dt_inner/2
          tlp1d(ix0:ix1) = tlp1d(ix0:ix1) +dtlp1d(ix0:ix1)*dt_inner/2
!.....2nd step
          call model_2tm1d(tnow,dtep1d,dtlp1d,iprint)
          te1d(ix0:ix1) = te1d(ix0:ix1) +dtep1d(ix0:ix1)*dt_inner
          tl1d(ix0:ix1) = tl1d(ix0:ix1) +dtlp1d(ix0:ix1)*dt_inner
        enddo
      endif
    endif

  end subroutine update_2tm1d
!=======================================================================
  subroutine model_2tm1d(tnow,dtep,dtlp,iprint)
!
!  Create difference of Te and Tl according to 1D TTM.
!  The TTM model of Ref.[1] is used for both Te and Tl.
!  
!   [1] Zhigilei, et al., J. Phys. Chem. C 113, 11892–11906 (2009)
!
    integer,intent(in):: iprint
    real(8),intent(in):: tnow
    real(8),intent(out):: dtep(nd1d),dtlp(nd1d)

    integer:: ix
    real(8):: ce,dce,kappa,dkappa,pterm,dtemp,tmp,xi,de,pulsefactor
    real(8):: denom
    real(8),parameter:: kappa_lat = 8.125d-7  ! kappa for Si lattice in eV/(fs.Ang.K)
    
    dtep(:) = 0d0
    dtlp(:) = 0d0

    do ix=ibc1d+1,nd1d
      ce = cete1d(ix)
      dce = dcete1d(ix)
      if( itype_kappa.eq.1 ) then  ! DCrho
        kappa = d_e *ce *rho_e
        dkappa = d_e *dce *rho_e
      else if( itype_kappa.eq.2 ) then  ! B2
        kappa = kappa0 *tep1d(ix) /max(tlp1d(ix),ta_min)
        dkappa = kappa0 /max(tlp1d(ix),ta_min)
      endif
      pterm = -gp1d(ix) *(tep1d(ix) -tlp1d(ix))
      denom = (ce +tep1d(ix)*dce) *rho_e
      dtemp = ( dkappa*dte21d(ix) +kappa*d2te1d(ix) +pterm ) /denom
      dtep(ix) = dtep(ix) +dtemp
!.....1/cl in lattice system may not be correct, since if cl=cl(Tl) and d(Tl)/dx!=0,
!.....as in the electronic system, the derivative d(cl)/dx!=0...
      dtlp(ix) = dtlp(ix) +(kappa*d2tl1d(ix) -pterm)/(cl1d*rho_bulk)
!!$      dtlp(ix) = dtlp(ix) -pterm/(cl1d*rho_bulk)
    enddo
!.....Laser pulse
    if( itype_pulse.eq.1 ) then  ! stepwise pulse
      if( tnow.ge.t0_laser .and. &
           tnow.le.(t0_laser +tau_pulse) ) then
        do ix=ibc1d+1,nd1d
          ce = cete1d(ix)
          dce = dcete1d(ix)
          tmp = 1d0 /((ce+tep1d(ix)*dce)*rho_e *vcell)
!.....Shift 0.5 since the mesh points are at the center of cells
          xi = (ix -1)*dx1d
          de = I_0 *min(1d0, exp(-xi/lskin))*dx
          dtep(ix) = dtep(ix) +de*tmp
        enddo
      endif
    else if( itype_pulse.eq.2 ) then  ! Gaussian pulse
      if( tnow.ge.t0_laser .and. &
           tnow.lt.(t0_laser +tau_pulse*2) ) then
        pulsefactor = exp(-(tnow -(t0_laser+tau_pulse))**2 /(2d0*sgm_pulse**2))
        do ix=ibc1d+1,nd1d
          ce = cete1d(ix)
          dce = dcete1d(ix)
          tmp = 1d0 /((ce+tep1d(ix)*dce)*rho_e *vcell)
!.....Shift 0.5 since the mesh points are at the center of cells
          xi = (ix -1)*dx1d
          de = I_0 *min(1d0, exp(-xi/lskin)) *dx *pulsefactor
          dtep(ix) = dtep(ix) +de*tmp
        enddo
      endif
    endif
  end subroutine model_2tm1d
!=======================================================================
  function cete1d(ix)
!
!  Ce(Te) at ix
!
    integer,intent(in):: ix
    real(8):: cete1d
    real(8):: t

    cete1d = 0d0
    if( iCe_Tdep.eq.0 ) then  ! none
      cete1d = c_0
    else if( iCe_Tdep.eq.1 ) then ! polynomial
      t = tep1d(ix)/1000
      cete1d = c_0 +(a_0 +a_1*t +a_2*t**2 +a_3*t**3 +a_4*t**4)&
           *exp(-(A_exp*t)**2) +ce_min
    else if( iCe_Tdep.eq.2 ) then ! tanh
      t = tep1d(ix)
      cete1d = 3d0 *tanh(2d-4 *t) +ce_min
    else if( iCe_Tdep.eq.3 ) then  ! linear
      cete1d = gmm_ce *tep1d(ix) +ce_min
    endif
    return
  end function cete1d
!=======================================================================
  function dcete1d(ix)
!
!  dCe(Te)/dTe at ix
!
    integer,intent(in):: ix

    real(8):: dcete1d
    real(8):: t,texp,x
    
    dcete1d = 0d0
    if( iCe_Tdep.eq.1 ) then  ! polynomial
      x = tep1d(ix)/1000
      texp = exp(-(A_exp*x)**2)
      dcete1d = (a_1 +2d0*a_2*x +3d0*a_3*x**2 +4d0*a_4*x**3)*texp &
           -2d0*A_exp*x *(a_0 +a_1*x +a_2*x**2 +a_3*x**3 +a_4*x**4)*texp
      dcete1d = dcete1d /1000
    else if( iCe_Tdep.eq.2 ) then  ! tanh
      t = tep1d(ix)
      dcete1d = 3d0*2d-4 *(1d0 -tanh(2d-4 *t)**2)
    else if( iCe_Tdep.eq.3 ) then  ! linear
      dcete1d = gmm_ce
    endif
    return
  end function dcete1d
!=======================================================================
  function dte21d(ix)
!
!  (dT/dx)^2 = [(T(ix+1)-T(ix-1))/2dx]^2
!            = [(T(ix+1)^2 -2*T(ix+1)*T(ix-1) -T(ix-1)^2)/(2dx)^2]
!
    integer,intent(in):: ix
    real(8):: dte21d
    real(8):: t,tp,tm

    t = tep1d(ix)
    tp= tep1d(ix+1)
    tm= tep1d(ix-1)
    dte21d = (tp**2 -2d0*tp*tm +tm**2)/4/dx1d**2
    return
  end function dte21d
!=======================================================================
  function d2te1d(ix)
!
!  (d^2/dx^2)Te
!
    integer,intent(in):: ix
    real(8):: d2te1d
    integer:: t,tp,tm

    t = tep1d(ix)
    tp= tep1d(ix+1)
    tm= tep1d(ix-1)
    d2te1d = (tp -2d0*t +tm)/dx1d**2
    return
  end function d2te1d
!=======================================================================
  function d2tl1d(ix)
!
!  (d^2/dx^2)Tl
!
    integer,intent(in):: ix
    real(8):: d2tl1d
    integer:: t,tp,tm

    t = tlp1d(ix)
    tp= tlp1d(ix+1)
    tm= tlp1d(ix-1)
    d2tl1d = (tp -2d0*t +tm)/dx1d**2
    return
  end function d2tl1d
!=======================================================================
  subroutine output_ttm1d(istp,tnow,myid,iprint)
!
!  Output temperatures (Te and Tl) in 1D TTM region.
!
    integer,intent(in):: istp,myid,iprint
    real(8),intent(in):: tnow

    integer:: ix
    character(len=128):: cnum
    
    if( myid.eq.0 ) then
      if( iprint.ne.0 ) then
        write(cnum,'(i0)') istp
        open(iots1d,file=trim(dname_ttm)// &
             trim(cout_ts1d)//'_'//trim(cnum),status='replace')
        write(iots1d,'(a,2es15.7,i6)') '#  tnow,dx,ibc1d: ',tnow,dx1d,ibc1d
        write(iots1d,'(a)') '#  ix,  te(ix),   tl(ix)'
        do ix=1,nd1d
          write(iots1d,'(2x,i6,2es15.5)') ix,te1d(ix),tl1d(ix)
        enddo
        close(iots1d)
      endif
    endif
    return
  end subroutine output_ttm1d
!=======================================================================
  subroutine couple_3d1d(myid,mpi_world,iprint)
!
!  Couple 3D-TTM-MD and 1D-TTM regions, which gives BCs to two systems..
!  See RK's note "1D-TTM for heat conducting medium" for details.
!
    integer,intent(in):: myid,mpi_world,iprint

    integer:: ix,iy,iz,mx,mxp,jx,ic,icp
    real(8):: tebc3d,tebcp3d,tabc3d,tabcp3d,x3d,x3dp,x1d,x1dp,tet3d,tat3d

!.....Get two ix3d points that sandwich the ibc1d point
    mxp = 0
    x1d = dx1d*(ibc1d-1)
    do ix= ibc3d-int(dx1d/dx)-1,ibc3d+1
      x3d = dx*(ix-1+0.5d0)
      if( x3d.gt.x1d ) then
        mxp = ix
        x3dp = x3d
        exit
      endif
    enddo
    mx = mxp -1
    x3d = dx*(mx-1+0.5d0)
!!$    print *,'mx,mxp,x3d,x3dp=',mx,mxp,x3d,x3dp

!.....Take averages of Te and Ta at mx and mxp over y and z
    tebc3d = 0d0
    tebcp3d = 0d0
    tabc3d = 0d0
    tabcp3d = 0d0
    do iy=1,ny
      do iz=1,nz
        tebc3d = tebc3d +te(mx,iy,iz)
        tebcp3d = tebcp3d +te(mxp,iy,iz)
        call ixyz2ic(mx,iy,iz,ic)
        tabc3d = tabc3d +ta(ic)
        call ixyz2ic(mxp,iy,iz,icp)
        tabcp3d = tabcp3d +ta(icp)
      enddo
    enddo
    tebc3d = tebc3d/(ny*nz)
    tebcp3d= tebcp3d/(ny*nz)
    tabc3d = tabc3d/(ny*nz)
    tabcp3d= tabcp3d/(ny*nz)
!!$    print '(a,2es12.4)','tebc3d,tabc3d=',tebc3d,tabc3d

!.....BC for 1D-TTM given from 3D-TTM system
!!$    te1d(ibc1d) = tebc3d
!!$    tl1d(ibc1d) = tabc3d
    te1d(ibc1d) = tebc3d +(tebcp3d-tebc3d)*(x1d-x3d)/(x3dp-x3d)
    tl1d(ibc1d) = tabc3d +(tabcp3d-tabc3d)*(x1d-x3d)/(x3dp-x3d)
    if( Te_right.lt.0d0 ) then
      te1d(nd1d+1) = te1d(nd1d)
      tl1d(nd1d+1) = tl1d(nd1d)
    else
      te1d(nd1d+1) = Te_right
      tl1d(nd1d+1) = Te_right
    endif
!!$    print *,'ibc1d,ibc3d=',ibc1d,ibc3d
!!$    print *,'te1d(ibc1d-2:ibc1d+2)=',te1d(ibc1d-2:ibc1d+2)
!!$    print *,'tl1d(ibc1d-2:ibc1d+2)=',tl1d(ibc1d-2:ibc1d+2)

!.....BC for 3D-TTM given from 1D-TTM system
    do ix=ibc3d+1,nx+1
      x3d = dx*(ix -1 +0.5d0)
      mx = 0
      do jx=ibc1d,nd1d
        x1d = dx1d*(jx-1)
        if( x1d.ge.x3d ) exit
        mx = jx
      enddo
      mxp = mx +1
      if( mx.eq.0 ) stop 'ERROR@couple_3d1d: could not find mx !!'
!.....Calc Te and Ta by linear interpolation from neighboring 1D-TTM node values
      x1d = dx1d*(mx-1)
      x1dp= dx1d*(mxp-1)
      tet3d = te1d(mx) +(te1d(mxp)-te1d(mx))*(x3d-x1d)/(x1dp-x1d)
      tat3d = tl1d(mx) +(tl1d(mxp)-tl1d(mx))*(x3d-x1d)/(x1dp-x1d)
!!$      print '(a,2i5,6es12.4)','ix,mx,te1d,,tet3d,tl1d,,tat3d= ',ix,mx &
!!$           ,te1d(mxp),te1d(mx),tet3d &
!!$           ,tat3d,tl1d(mxp),tl1d(mx)
      do iy=0,ny+1
        do iz=0,nz+1
          te(ix,iy,iz) = tet3d
          call ixyz2ic(ix,iy,iz,ic)
          ta(ic) = tat3d
        enddo
      enddo
    enddo
    
  end subroutine couple_3d1d
!=======================================================================
  subroutine set_3d1d_bc_pos(natm,ra,h,sorg,myid,mpi_world,iprint)
!
!  Determine boundary position between 3D-TTM and 1D-TTM systems.
!  See RK's note "1D-TTM for heat conducting medium" for details.
!
    integer,intent(in):: natm,myid,mpi_world,iprint
    real(8),intent(in):: ra(3,natm),h(3,3),sorg(3)

    integer:: i,ierr
    real(8):: xrl,xdnr,x1d,x3d

!.....Get right-most x-pos of atoms, xr
    xrl = 0d0
    do i=1,natm
      xrl = max(ra(1,i)+sorg(1),xrl)
    enddo
    call mpi_allreduce(xrl,xrmd,1,mpi_real8,mpi_max,mpi_world,ierr)

!.....Get x-index of boundary cell in 3D-TTM system (ibc3d) from (xr - dnr), 
!.....where dnr is the length of non-reflecting boundary region
    xdnr = dnr/h(1,1)
    ibc3d = int(nx*(xrmd-xdnr)) +1
    if( ibc3d.ge.nx ) ibc3d = nx-1  ! boundary position must be .lt. nx
!.....x3d as a right-edge position of boundary cell which is a bit different from the note
    x3d = dx*ibc3d
!!$    print *,'xdnr,ibc3d,x3d = ', xdnr, ibc3d, x3d
    
!.....Get x-index of 1D-TTM system (ibc1d) from x_{ibc3d}
    ibc1d = 0
    do i=1,nd1d
      x1d = dx1d*(i-1)
      if( x1d.gt.x3d ) exit
      ibc1d = i
    enddo
!!$    print *,'ibc1d,x1d = ', ibc1d, dx1d*(ibc1d-1)

  end subroutine set_3d1d_bc_pos
end module ttm
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
