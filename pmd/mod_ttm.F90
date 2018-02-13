module ttm
!-----------------------------------------------------------------------
!                     Last-modified: <2018-02-08 17:37:27 Ryo KOBAYASHI>
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
  integer,parameter:: ioprms = 30
  integer,parameter:: ioTein = 31
  integer,parameter:: ioTeout = 32

  real(8):: t_ttm
!.....TTM mesh divisions
  integer:: nx,ny,nz,nxyz
!.....Mesh size in reduced unit [0:1)
  real(8):: dx,dy,dz,area,darea
!.....Time step in fs == dt in MD by default
  real(8):: dt
!.....Number of inner loop in TTM
  integer:: nstp_inner = 1
!.....Volume per mesh cell Ang^3 == dx*dy*dz
  real(8):: vcell
!.....Threshold kinetic energy in energy unit (or should be threshold velocity?)
  real(8):: ekth = 8.0d0
!.....Te of right edge, if negative, use free boundary
  real(8):: Te_right = -1.0

!.....Ce dependence on Te: none, polynomial or tanh
  character(len=128):: Ce_Tdep = 'none'
  real(8):: rho_e = 0.005
  real(8):: d_e = 20000d0
  real(8):: c_0 = 1e-4
!.....Coefficients for polynomial Ce(Te)
  real(8):: a_0 = 0d0
  real(8):: a_1 = 0d0
  real(8):: a_2 = 0d0
  real(8):: a_3 = 0d0
  real(8):: a_4 = 0d0
  real(8):: A_exp = 0d0

!.....Pulse shape: exp or read (from cTe_init)
  character(len=128):: cTe_init = 'exp'
!.....Laser power at the surface top in eV/(fs*Ang^2) unit
  real(8):: I_0 = 5d+4
!.....Pulse duration in fs
  real(8):: tau_pulse = 100d0
!.....Surface skin length in Ang
  real(8):: lskin = 100d0
!.....Surface position (positive integer)
  integer:: lsurf = -1  ! left surface
  integer:: rsurf = -1  ! right surface
!.....Surface movement along x: no, plane (can move as a yz-plane surface)
  character(len=128):: surfmove = 'plane'
  integer:: nstp_surfmove = 10
!.....Minimum Te when surface moves
  real(8):: Te_min = 300d0

!.....Temperature distribution
  real(8),allocatable:: te(:,:,:),tep(:,:,:),ta(:),tap(:),tex(:)
  integer,allocatable:: nac(:),nacp(:)
  real(8),allocatable:: eksum(:),ekpsum(:)
!.....Atom to cell correspondance
  integer,allocatable:: a2c(:)

!.....Sigma of random force in Langevin thermostat
  real(8),allocatable:: sgm(:)
!.....Gammas
  real(8):: gmmp, gmms
!.....gp, gs
  real(8),allocatable:: gp(:),gs(:)

  real(8),allocatable:: aai(:,:)

contains
!=======================================================================
  subroutine init_ttm(namax,natm,h,dtmd,myid,mpi_world,iprint)
!
!  Read parameters for TTM from in.ttm and initialize
!
    integer,intent(in):: namax,natm,myid,mpi_world,iprint
    real(8),intent(in):: dtmd,h(3,3)

    integer:: ierr,ix,iy,iz,mem
    real(8):: t,t0,dtmax,tmp
    character(len=128):: c1st

    t_ttm = 0d0
    t0 = mpi_wtime()
!.....Read parameter file
    call read_ttm_params(myid,mpi_world,iprint)
    call sync_params(myid,mpi_world,iprint)

!.....Convert fluence to energy density
    I_0 = I_0 /lskin

!.....Set some
    nxyz = nx*ny*nz
    dt = dtmd /nstp_inner
    dx = h(1,1)/nx
    dy = h(2,2)/ny
    dz = h(3,3)/nz
    vcell = dx*dy*dz
    area = h(2,2)*h(3,3)
    darea = dy*dz

    if( myid.eq.0 .and. iprint.ne.0 ) then
      print *,''
      print *,'TTM parameters:'
      print '(a,i5)','   inner_loop = ',nstp_inner
      print '(a,3i5,i8)','   nx,ny,nz,nxyz = ',nx,ny,nz,nxyz
      print '(a,4es12.4)','   dx,dy,dz,vcell = ',dx,dy,dz,vcell
      print '(a,2es12.4)','   area,darea = ',area,darea
      print '(a,2i5,es12.4)','   lsurf,rsurf = ',lsurf,rsurf,lskin
      print '(a,2es12.4)','   gmmp,gmms = ',gmmp,gmms
      print '(a,a)','   Ce_Tdep = ',trim(Ce_Tdep)
      mem = 4 * 4*nxyz + 11 * 8*nxyz + 4 * 8*namax
      print '(a,f0.3,a)',' Memory for TTM = ',dble(mem)/1000/1000,' MByte'
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
         ta(nxyz),tap(nxyz),tex(nx),gp(nxyz),gs(nxyz))
    allocate(a2c(namax),aai(3,namax))

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
    else if( trim(cTe_init).eq.'read' ) then
      if( myid.eq.0 ) then
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
    endif

    t_ttm = t_ttm +mpi_wtime() -t0
    return

999 call mpi_finalize(ierr)
    stop
  end subroutine init_ttm
!=======================================================================
  subroutine read_ttm_params(myid,mpi_world,iprint)
    integer,intent(in):: myid,mpi_world,iprint
    integer,external:: num_data

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
        else if( trim(c1st).eq.'gamma_p' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, gmmp
        else if( trim(c1st).eq.'gamma_s' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, gmms
        else if( trim(c1st).eq.'ekth' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, ekth
        else if( trim(c1st).eq.'rho_e' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, rho_e
        else if( trim(c1st).eq.'D_e' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, d_e
        else if( trim(c1st).eq.'Te_init' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, cTe_init
        else if( trim(c1st).eq.'Te_right' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, Te_right
        else if( trim(c1st).eq.'pulse_duration' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, tau_pulse
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
        else if( trim(c1st).eq.'surface_move' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, surfmove
        else if( trim(c1st).eq.'laser_power' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, I_0
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
    call mpi_bcast(gmmp,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(gmms,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(ekth,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(rho_e,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(d_e,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(tau_pulse,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(cTe_init,128,mpi_character,0,mpi_world,ierr)
    call mpi_bcast(Te_right,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(ce_Tdep,128,mpi_character,0,mpi_world,ierr)
    call mpi_bcast(c_0,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(a_0,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(a_1,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(a_2,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(a_3,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(a_4,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(A_exp,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(surfmove,1,mpi_logical,0,mpi_world,ierr)
    call mpi_bcast(I_0,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(Te_min,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(lsurf,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(rsurf,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(lskin,1,mpi_real8,0,mpi_world,ierr)
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
  subroutine calc_Ta(namax,natm,eki,tag,fmv,istp,myid,mpi_world)
!
!  Compute and set Ta and Tap array from atomic kinetic energies.
!
    integer,intent(in):: namax,natm,myid,mpi_world,istp
    real(8),intent(in):: eki(3,3,namax),tag(namax),fmv(3,0:9)

    integer:: i,ic,ierr,ix,iy,iz,l,ifmv,idof
    real(8):: ek,t0
    integer,allocatable,save:: nacl(:),nacpl(:)
    real(8),allocatable,save:: eksuml(:),ekpsuml(:)
    integer,external:: ifmvOf

    if( .not. allocated(nacl) ) then
      allocate(nacl(nxyz),nacpl(nxyz),eksuml(nxyz),ekpsuml(nxyz))
    endif
    
    t0 = mpi_wtime()

    nacl(1:nxyz) = 0
    nacpl(1:nxyz) = 0
    eksuml(1:nxyz) = 0d0
    ekpsuml(1:nxyz) = 0d0
    do i=1,natm
      ek = eki(1,1,i) +eki(2,2,i) +eki(3,3,i)
      ic = a2c(i)
      ifmv = ifmvOf(tag(i))
      idof = 0
      do l=1,3
        idof = idof +nint(fmv(l,ifmv))
      enddo
      nacl(ic) = nacl(ic) + idof
      eksuml(ic) = eksuml(ic) +ek
!!$      if( ic.eq.4 ) print *,'i,ek,eksum=',i,ek,eksuml(ic)
      if( ek.gt.ekth ) then
        nacpl(ic) = nacpl(ic) +idof
        ekpsuml(ic) = ekpsuml(ic) +ek
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
      gp(:) = 0d0
      gs(:) = 0d0
      ta(:) = 0d0
      tap(:) = 0d0
      do ic=1,nxyz
        if( nac(ic).eq.0 ) cycle
        call ic2ixyz(ic,ix,iy,iz)
!!$        if( iy.eq.1 .and. iz.eq.1 ) then
!!$          print *,'ix,nac(ic)=',ix,nac(ic)
!!$        endif
!.....Degree of freedom per atom (3 in case of 3D) is included in nac
        ta(ic) = eksum(ic) *2d0 /fkb /nac(ic)
        gp(ic) = nac(ic) *fkb *gmmp /vcell
!!$        if( ic.eq.4 ) print *,'Ta: ic,eksum,nac,ta,gp=',ic,eksum(ic),nac(ic),ta(ic),gp(ic)
        if( nacp(ic).eq.0 ) cycle
        tap(ic) = ekpsum(ic) *2d0 /fkb /nacp(ic)
        gp(ic) = nacp(ic) *fkb *gmms /vcell
      enddo
!!$      do ix=1,nx
!!$        iy = 1
!!$        iz = 1
!!$        call ixyz2ic(ix,iy,iz,ic)
!!$        print *,'ix,iy,iz,ic,ta,eksum,eksuml=',ix,iy,iz,ic,ta(ic) &
!!$             ,eksum(ic),eksuml(ic)
!!$      enddo
    endif

    if( trim(surfmove).eq.'plane' .and. &
         mod(istp,nstp_surfmove).eq.0) then
      call update_surface_plane(myid,mpi_world)
    endif

    t_ttm = t_ttm +mpi_wtime() -t0
    return
  end subroutine calc_Ta
!=======================================================================
  subroutine update_Te(tnow,myid,mpi_world)
!
!  Update Te by solving the diffusion equation.
!
    integer,intent(in):: myid,mpi_world
    real(8),intent(in):: tnow

    integer:: ic,ix,iy,iz,ierr,istp
    real(8):: t0,ce,xi

    t0 = mpi_wtime()

    if( myid.eq.0 ) then
      do istp = 1,nstp_inner
        call set_BC()

        do ic=1,nxyz
          call ic2ixyz(ic,ix,iy,iz)
          tep(ix,iy,iz) = 0d0
          if( ix.lt.lsurf ) cycle
          if( ix.gt.rsurf ) cycle
          ce = cete(ix,iy,iz)
          tep(ix,iy,iz) = te(ix,iy,iz) &
               +dt/ce/rho_e &
               *(rho_e*d_e *( dcete(ix,iy,iz)*dte2(ix,iy,iz) &
               +ce*d2te(ix,iy,iz) ) &
               -gp(ic)*(te(ix,iy,iz) -ta(ic)) +gs(ic)*tap(ic) )

        enddo
        if( tnow.lt.tau_pulse ) then
          do ic=1,nxyz
            call ic2ixyz(ic,ix,iy,iz)
            if( ix.lt.lsurf ) cycle
            if( ix.gt.rsurf ) cycle
            ce = cete(ix,iy,iz)
            xi = (ix-lsurf)*dx
!.....TODO: check consistency of units
            tep(ix,iy,iz) = tep(ix,iy,iz) &
                 +I_0 *min(1d0,exp(-xi/lskin))/ce/rho_e*dt
!!$            if( ic.eq.126 ) then
!!$              print *,'ic,I_0,ce,other=',ic,I_0,ce &
!!$                   ,min(1d0,exp(-xi/lskin))/ce/rho_e*dt,1d0/ce,1d0/rho_e
!!$            endif
          enddo
        endif
        te(:,:,:) = tep(:,:,:)
      enddo
    endif
!.....Broadcast Te distribution to all the nodes.
!.....There could be smarter way to reduce networking cost.
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
    if( trim(Ce_Tdep).eq.'none' ) then
      cete = c_0
    else if( trim(Ce_Tdep).eq.'polynomial' ) then
      t = te(ix,iy,iz)/1000
      cete = c_0 +(a_0 +a_1*t +a_2*t**2 +a_3*t**3 +a_4*t**4)&
           *exp(-(A_exp*t)**2)
    else if( trim(Ce_Tdep).eq.'tanh' ) then
      cete = 3d0 *tanh(2d-4 *t)
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
    if( trim(Ce_Tdep).eq.'polynomial' ) then
      t = te(ix,iy,iz)
      texp = exp(-(A_exp*t)**2)
      dcete = (a_1 +2d0*a_2*t +3d0*a_3*t**2 +4d0*a_4*t**3)*texp &
           -2d0*A_exp*t *(a_0 +a_1*t +a_2*t**2 +a_3*t**3 +a_4*t**4)*texp
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
    dte2 = (te(ixp,iy,iz)**2 -2d0*te(ixp,iy,iz)*te(ixm,iy,iz) &
         +te(ixm,iy,iz))/4/dx/dx &
         + (te(ix,iyp,iz)**2 -2d0*te(ix,iyp,iz)*te(ix,iym,iz) &
         +te(ix,iym,iz))/4/dy/dy &
         + (te(ix,iy,izp)**2 -2d0*te(ix,iy,izp)*te(ix,iy,izm) &
         +te(ix,iy,izm))/4/dz/dz
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

    d2te = (te(ixp,iy,iz) -2d0*te(ix,iy,iz) &
         +te(ixm,iy,iz)) /dx/dx &
         + (te(ix,iyp,iz) -2d0*te(ix,iy,iz) &
         +te(ix,iym,iz)) /dy/dy &
         + (te(ix,iy,izp) -2d0*te(ix,iy,iz) &
         +te(ix,iy,izm)) /dz/dz
    return
  end function d2te
!=======================================================================
  subroutine langevin_ttm(namax,natm,va,aa,tag,eki,am,h, &
       nspmax,fa2v,fekin,ediff,dtmd)
!
!  Langevin thermostat for atomic system.
!
    integer,intent(in):: namax,natm,nspmax
    real(8),intent(in):: aa(3,namax),tag(namax),am(nspmax) &
         ,fa2v(nspmax),fekin(nspmax),dtmd,h(3,3),eki(3,3,namax)
    real(8),intent(inout):: va(3,namax),ediff(nspmax)

    integer:: ic,i,l,is,ifmv,ix,iy,iz,naccp
    real(8):: hscl(3),sgmi,ami,ek,gmmi,vl(3),vi(3),aai(3),t0
    integer,external:: ifmvOf
    real(8),external:: box_muller,sprod

    t0 = mpi_wtime()

    do ic=1,nxyz
      call ic2ixyz(ic,ix,iy,iz)
      sgm(ic) = dsqrt(2d0*gmmp*fkb*te(ix,iy,iz)/dtmd)
!!$      if( ic.eq.4 ) print *,'ic,sgm,gmmp,fkb*Te=',ic,sgm(ic)&
!!$           ,gmmp,fkb*te(ix,iy,iz)
    enddo
    
!.....Langevin thermostat with Mannella integrator
    hscl(1:3)= 0d0
    do l=1,3
      hscl(l)= dsqrt(h(1,l)**2 +h(2,l)**2 +h(3,l)**2)
    enddo
    do i=1,natm
      ifmv = ifmvOf(tag(i))
      is = int(tag(i))
      if( ifmv.eq.0 ) then
        va(1:3,i)= 0d0
      else
        va(1:3,i)=va(1:3,i) +aa(1:3,i)*fa2v(is)*dtmd
        ami= am(is)
        ic = a2c(i)
        sgmi = sgm(ic) *dsqrt(ami)
!!$        if( ic.eq.4 ) then
!!$          print *,'ia,ic,ami,sgmi=',i,ic,ami,sgmi
!!$        endif
        ek = eki(1,1,i) +eki(2,2,i) +eki(3,3,i)
        gmmi = gmmp
        if( ek.gt.ekth ) gmmi = gmmi + gmms
        aai(1:3)= 0d0
        do l=1,3
          aai(l)= -va(l,i)*gmmi*ami +sgmi*box_muller()/hscl(l)
        enddo
!.....To compensate the factor 1/2 in fa2v, multiply 2 here.
        va(1:3,i)= va(1:3,i) +aai(1:3)*fa2v(is)*dtmd *2d0
!.....accumulate energy difference
        vi(1:3)= h(1:3,1)*va(1,i) &
             +h(1:3,2)*va(2,i) &
             +h(1:3,3)*va(3,i)
        vl(1:3)= h(1:3,1)*aai(1) &
             +h(1:3,2)*aai(2) &
             +h(1:3,3)*aai(3)
        ediff(ifmv)= ediff(ifmv) +fekin(is) &
             *(2d0*sprod(3,vi,vl)+sprod(3,vl,vl))
      endif
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
  subroutine output_Te(istp,myid,iprint)
    integer,intent(in):: istp,myid,iprint

    integer:: ix,iy,iz,ic
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
    endif
  end subroutine output_Te
!=======================================================================
  subroutine set_BC()
!
!  Set boundary condition
!
    if( lsurf.le.0 .and. rsurf.le.0 ) then
!.....Periodic for x,y,z
      te(0,1:ny,1:nz) = te(nx,1:ny,1:nz)
      te(nx+1,1:ny,1:nz) = te(1,1:ny,1:nz)
      te(1:nx,0,1:nz) = te(1:nx,ny,1:nz)
      te(1:nx,ny+1,1:nz) = te(1:nx,1,1:nz)
      te(1:nx,1:ny,0) = te(1:nx,1:ny,nz)
      te(1:nx,1:ny,nz+1) = te(1:nx,1:ny,1)
    else  ! Laser-ablation situation
!.....Free boundary for x
      te(lsurf-1,1:ny,1:nz) = te(lsurf,1:ny,1:nz)
      if( Te_right.lt.0d0 ) then
        te(rsurf+1,1:ny,1:nz) = te(rsurf,1:ny,1:nz)
      else
        te(rsurf+1,1:ny,1:nz) = Te_right
      endif
!.....Periodic for y and z
      te(lsurf:rsurf,0,1:nz) = te(lsurf:rsurf,ny,1:nz)
      te(lsurf:rsurf,ny+1,1:nz) = te(lsurf:rsurf,1,1:nz)
      te(lsurf:rsurf,1:ny,0) = te(lsurf:rsurf,1:ny,nz)
      te(lsurf:rsurf,1:ny,nz+1) = te(lsurf:rsurf,1:ny,1)
    endif
  end subroutine set_BC
!=======================================================================
  subroutine update_surface_plane(myid,mpi_world)
!
!  Update surface position according to number of atoms in the cells.
!
    integer,intent(in):: myid,mpi_world
    integer:: ic,icl,icr,ix,iy,iz,ierr
    integer,allocatable,save:: nacl(:,:,:),nacr(:,:,:)
    logical:: lupdm,lupdp,rupdm,rupdp

    if( .not. allocated(nacl) ) then
      allocate(nacl(nx,ny,nz),nacr(nx,ny,nz))
    endif

    if( myid.eq.0 ) then

!.....Count number of atoms in the cells next to lsurf/rsurf
      nacl(:,:,:) = 0
      nacr(:,:,:) = 0
      do iz=1,nz
        do iy=1,ny
          if( lsurf-1.ge.1 ) then
            do ix=lsurf-1,lsurf
              call ixyz2ic(ix,iy,iz,icl)
              nacl(ix,iy,iz) = nac(icl)
            enddo
          endif
          if( rsurf+1.le.nx ) then
            do ix=rsurf,rsurf+1
              call ixyz2ic(ix,iy,iz,icr)
              nacr(ix,iy,iz) = nac(icr)
            enddo
          endif
        enddo
      enddo
!.....Check whether all the cells in y-z plane contains at least 2 atoms
      lupdm = .true.
      lupdp = .true.
      rupdm = .true.
      rupdp = .true.
      do iz=1,nz
        do iy=1,ny
          if( nacl(lsurf-1,iy,iz) .lt. 2 ) lupdm = .false.
          if( nacl(lsurf  ,iy,iz) .ge. 2 ) lupdp = .false.
          if( nacr(rsurf+1,iy,iz) .lt. 2 ) rupdp = .false.
          if( nacr(rsurf  ,iy,iz) .ge. 2 ) rupdm = .false.
        enddo
      enddo
!.....Update surface if needed
      if( lupdm ) then
        lsurf = lsurf - 1
        te(lsurf,:,:) = te(lsurf+1,:,:)
      else if( lupdp ) then
        lsurf = lsurf + 1
        te(lsurf-1,:,:) = 0d0
      endif
      if( rupdp ) then
        rsurf = rsurf + 1
        te(rsurf,:,:) = te(rsurf-1,:,:)
      else if( lupdm ) then
        rsurf = rsurf - 1
        te(rsurf+1,:,:) = 0d0
      endif
    endif
!.....Broadcast Te distribution to all the nodes.
!.....There could be smarter way to reduce networking cost.
    call mpi_bcast(te,(nx+2)*(ny+2)*(nz+2),mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(rsurf,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(lskin,1,mpi_real8,0,mpi_world,ierr)
  end subroutine update_surface_plane
end module ttm
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
