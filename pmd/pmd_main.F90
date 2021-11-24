program pmd
!-----------------------------------------------------------------------
!                     Last-modified: <2021-11-24 21:36:23 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
! Spatial decomposition parallel molecular dynamics program.
! Core part is separated to pmd_core.F.
!-----------------------------------------------------------------------
! INPUT FILES:
! ------------
!   in.pmd:     Main input file
!   pmdini:     Cell info and atom coordinations
!
! OUTPUT FILES:
! -------------
!   out.erg:    Total, kinetic, and potential energy
!   pmd_####:   Cell info and atom coordinations of a certain steps
!               in the MD run if required by "flag_out_pmd" in in.pmd.
!   out.stress: Stress component normal to z-upper surface of nanorod,
!               and z-strain of the nanorod.
!-----------------------------------------------------------------------
  use pmdvars
  use pmdio,only: read_pmdtot_bin, read_pmdtot_ascii, write_pmdtot_bin, &
       write_pmdtot_ascii,get_ntot_ascii,get_ntot_bin
  use force
  use Coulomb, only: cterms
  use util, only: itotOf, cell_info, iauxof, make_cdumpauxarr, spcs_info
  use time, only: time_stamp, accum_time, report_time
  use memory, only: accum_mem, report_mem
  use element
  use random,only: urnd, set_seed
  use clrchg,only: lclrchg,init_clrchg
  use localflux,only: lflux,init_lflux,final_lflux
  use pdens,only: lpdens,init_pdens,final_pdens
!$  use omp_lib
  implicit none
  include "mpif.h"
  include "./params_unit.h"
  include "./const.h"

#ifdef __DISL__
!.....Epot threshold for disl core extraction [Hartree]
  real(8),parameter:: epith = -0.1410d0
#endif

  real(8):: hunit,hmat(3,3,0:1)
  integer:: ntot0,ntot
  real(8),allocatable:: tagtot(:),rtot(:,:),vtot(:,:),atot(:,:)
  real(8),allocatable:: stot(:,:,:),epitot(:),ekitot(:,:,:)
  real(8),allocatable:: auxtot(:,:)

  integer:: i,j,k,l,m,n,ia,ib,is,ifmv,nave,nspl,i_conv,inc
  integer:: mpicolor,mpikey,ierr,jerr,itmp,nprocs,nnmax_est,mem
  real(8):: tmp,hscl(3),aai(3),ami,dt2,tave,vi(3),vl(3),rmin
  real(8):: epot,ekin,stnsr(3,3)
  real(8):: t0,t1
  character(len=3):: csp
  type(atom):: elem

!-----initialize the MPI environment
  call mpi_init(ierr)
!-----total number of MD-nodes
  call mpi_comm_size(MPI_COMM_WORLD,nprocs,ierr)
!-----my rank in MD-nodes
  call mpi_comm_rank(MPI_COMM_WORLD,myid_md,ierr)
  call mpi_comm_dup(MPI_COMM_WORLD,mpicomm,ierr)
  mpi_md_world = mpicomm
  t0 = mpi_wtime()

  call init_element()

!.....Set fmv as default value before reading 'in.pmd'
!!$  call set_fmv(fmv)

  if( myid_md.eq.0 ) then
    call write_headline()

!.....Read atom configuration file 1st
    if( trim(ciofmt).eq.'bin' .or. trim(ciofmt).eq.'binary' ) then
      write(6,*) 'Read pmdini in binary mode.'
      ntot0 = get_ntot_bin(20,trim(cpmdini))
      allocate(tagtot(ntot0),rtot(3,ntot0),vtot(3,ntot0),epitot(ntot0) &
           ,ekitot(3,3,ntot0),stot(3,3,ntot0),atot(3,ntot0))
      call read_pmdtot_bin(20,trim(cpmdini),ntot0,hunit,hmat,tagtot,rtot,vtot)
    else if( trim(ciofmt).eq.'ascii' ) then
      write(6,*) 'Read pmdini in ascii mode.'
      ntot0 = get_ntot_ascii(20,trim(cpmdini))
      allocate(tagtot(ntot0),rtot(3,ntot0),vtot(3,ntot0),epitot(ntot0) &
           ,ekitot(3,3,ntot0),stot(3,3,ntot0),atot(3,ntot0))
      call read_pmdtot_ascii(20,trim(cpmdini),ntot0,hunit,hmat,tagtot,rtot,vtot)
    else
      write(6,*) 'Error: io_format must be either ascii, ' &
           //'bin or binary.'
      stop
    endif
!.....Memory assessment
    mem = 8*ntot0*(1 +3 +3 +1 +3*3 +3*3 +3)
    call accum_mem('main',mem)

!.....Set mass of species if specorder is already set.
!.....This could be overwritten by mass entry in in.pmd
    if( has_specorder ) then
      do is=1,nspmax
        csp = specorder(is)
        if( trim(csp).ne.'x' ) then
          elem = get_element(trim(csp))
          am(is) = elem%mass
        endif
      enddo
    endif

    call cell_info(h)
    call spcs_info(ntot0,tagtot)
    write(6,*) ''
    write(6,'(a,i0)') ' Num of MPI processes = ',nprocs
!.....Read in.pmd after reading the atom configuration file.
    call read_inpmd(10,trim(cinpmd))
!$    if( nomp.gt.0 ) then
!$      call omp_set_num_threads(nomp)
!$    endif
!$omp parallel
!$omp single
!$    write(6,'(a,i0)') ' Num of OpenMP processes = ',omp_get_num_threads()
!$omp end single
!$omp end parallel

    call check_cmin()
    if( ifpmd.eq.2 ) then ! if dump output
      call make_cdumpauxarr()
    endif
    call write_initial_setting()
!        call write_inpmd(10,trim(cinpmd))
    if( num_forces.eq.0 ) stop ' ERROR: no force-field specified'

!.....Initialize random seeds in the function urnd
    if( rseed.ge.0d0 ) call set_seed(rseed+myid_md)
    tmp = urnd()
    
    if( trim(ctctl).eq.'ttm' ) then
      print *,''
      print *,'NOTICE: Since using the two-temperature model (TTM) MD:'
!.....Set x-boundary free if TTM
      if( boundary(1:1).ne.'f' ) then
        print *,'  - Free boundary condition is set' &
             //' for x direction.'
      endif
      boundary(1:1) = 'f'
!.....Damping off if TTM
      print *,'  - Set damping off.'
      ifdmp = 0
!.....Set nrmtrans to 0
      print *,'  - Set removal of translation off.'
      nrmtrans = 0
    endif

!.....Some FF requires other FFs
    if( use_force('BMH') ) then
      if( .not. use_force('Coulomb') ) then
        num_forces = num_forces +1
        force_list(num_forces) = 'Coulomb'
      endif
      if( .not. use_force('dipole') ) then
        num_forces = num_forces +1
        force_list(num_forces) = 'dipole'
      endif
    endif
    if( use_force('fpc') ) then
      if( .not. use_force('Coulomb') ) then
        num_forces = num_forces +1
        force_list(num_forces) = 'Coulomb'
      endif
!!$      if( .not. use_force('LJ_repul') ) then
!!$        num_forces = num_forces +1
!!$        force_list(num_forces) = 'LJ_repul'
!!$      endif
    endif

!.....Check whether localflux is used with color charge NEMD
    if( lflux .and. .not. lclrchg ) then
      print *,'ERROR: local flux must be used with color-charge NEMD !'
      stop
    endif
    if( lflux .or. lpdens ) then
!.....Set cutoff_buffer to zero, since it can affect local flux results
      rbuf = 0d0
      print *,''
      print *,'NOTICE: cutoff_buffer is reset to 0, since it can affect ' &
           //'the following:'
      if( lflux ) print *,'  - local flux'
      if( lpdens ) print *,'  - probability density'
    endif

!.....Correct nnmax if the given nnmax is too small compared to
!.....the estimated one
!.....Now assume the minimum interatomic distance is about 2.0A
    rmin = 2.0d0
    nnmax_est = (rc+rbuf)**3 /rmin**3 !*dsqrt(2d0)*pi/6
    if( nnmax.lt.nnmax_est ) then
      print *,'NNMAX is replaced since it is too small w.r.t. '// &
           'given cutoff radius.'
      print '(a,2(2x,i0))','    nnmax_orig, nnmax_new = ',nnmax,nnmax_est
      nnmax = nnmax_est
    endif
  endif  ! end of myid.eq.0

  call bcast_params()
  call mpi_bcast(hunit,1,mpi_real8,0,mpi_md_world,ierr)
  call mpi_bcast(hmat,9*2,mpi_real8,0,mpi_md_world,ierr)

!.....Before allocating auxiliary array, set naux (num of auxiliary data)
  call set_use_charge()
  call set_use_elec_temp()
  naux = 0
  if( luse_charge ) then
    naux = naux +1  ! chg
  endif
  if( luse_elec_temp .or. trim(ctctl).eq.'ttm' ) then
    naux = naux +1
  endif
  if( lclrchg ) then
    naux = naux +1
  endif
  allocate(cauxarr(naux))
  inc = 0
  if( luse_charge ) then
    inc = inc +1
    cauxarr(inc) = 'chg'
  endif
  if( luse_elec_temp .or. trim(ctctl).eq.'ttm' ) then
    inc = inc +1
    cauxarr(inc) = 'tei'
  endif
  if( lclrchg ) then
    inc = inc +1
    cauxarr(inc) = 'clr'
  endif

!.....Now allcoate the auxiliary array
  if( myid_md.eq.0 ) then
    if( .not. allocated(auxtot) ) then
      allocate(auxtot(naux,ntot0))
    else
      if( size(auxtot).ne.naux*ntot0 ) deallocate(auxtot)
      allocate(auxtot(naux,ntot0))
    endif
    auxtot(:,:) = 0d0
!.....Memory assessment
    mem = 8*naux*ntot0
    call accum_mem('main',mem)

!.....Determine nx,ny,nz using rc and hmat info
    if( .not. (nx.gt.0 .and. ny.gt.0 .and. nz.gt.0 ) ) then
      call determine_division(hmat,myid_md,nprocs,rc,nx,ny,nz,iprint)
    endif

!.....Make ntot and ?tot() not null in nodes myid_md != 0
  else
    ntot0 = 1
    allocate(tagtot(ntot0),rtot(3,ntot0),vtot(3,ntot0),epitot(ntot0) &
         ,ekitot(3,3,ntot0),stot(3,3,ntot0),atot(3,ntot0) )
    allocate(auxtot(naux,ntot0))
  endif
  ntot = ntot0

!.....Broadcast species data read from pmdini  
  call mpi_bcast(specorder,3*nspmax,mpi_character,0,mpicomm,ierr)
  call mpi_bcast(has_specorder,1,mpi_logical,0,mpicomm,ierr)

!.....Broadcast determined nx,ny,nz to reset MPI communicator if needed.
  call mpi_bcast(nx,1,mpi_integer,0,mpicomm,ierr)
  call mpi_bcast(ny,1,mpi_integer,0,mpicomm,ierr)
  call mpi_bcast(nz,1,mpi_integer,0,mpicomm,ierr)
  if( nx*ny*nz .gt. nprocs ) then
    if( myid_md.eq.0 ) then
      print '(a,2i5)',' ERROR: nxyz > nprocs, which should not happen!, nxyz,nprocs = '&
           ,nx*ny*nz,nprocs
    endif
  else if( nx*ny*nz .lt. nprocs ) then
    nodes_md = nx*ny*nz
    if( myid_md.eq.0 .and. iprint.ge.ipl_basic ) then
      print '(a)',' WARNING: Since nxyz < nprocs, use less processes than prepared.'
      print '(a,5(2x,i0))','          nx,ny,nz,nxyz,nprocs=',nx,ny,nz,nodes_md,nprocs
    endif
    mpikey = myid_md
    if( myid_md.lt.nodes_md ) then
      mpicolor = 0
    else
      mpicolor = 1
    endif
    call mpi_comm_split(mpicomm,mpicolor,mpikey,mpi_md_world,jerr)
    call mpi_comm_size(mpi_md_world,itmp,ierr)
    if( myid_md.lt.nodes_md .and. itmp.ne.nodes_md ) &
         stop 'ERROR: itmp.ne.nodes_md in mpi_comm_size'
    call mpi_comm_rank(mpi_md_world,itmp,ierr)
    if( myid_md.lt.nodes_md .and. itmp.ne.myid_md ) &
         stop 'ERROR: itmp.ne.myid in mpi_comm_rank'
!!$    print *,' myid,itmp,jerr,color,mpi_md_world,mpicomm=',myid_md,itmp &
!!$         ,jerr,mpicolor,mpi_md_world,mpicomm
  else
    nodes_md = nprocs
    mpi_md_world = mpicomm  ! it is already set
  endif

!.....If the node is not used, skip to the end and wait the other nodes
  if( myid_md .ge. nodes_md ) goto 1

!.....Check ensemble
  if( myid_md.eq.0 ) call check_ensemble()

!.....Set Coulomb flag here
  if( use_force('screened_Coulomb') ) then
    ifcoulomb = 1
    cterms = 'screened'
  else if( use_force('Ewald') ) then
    ifcoulomb = 2
    cterms = 'full'
  else if( use_force('Ewald_long') ) then
    ifcoulomb = 3
    cterms = 'long'
  endif

!.....Initial settting for color charge NEMD
  if( lclrchg ) call init_clrchg(specorder,ntot0,auxtot(iauxof('clr'),:),tagtot &
       ,myid_md,iprint)
!.....Init for local flux
  if( lflux ) call init_lflux(myid_md,nx,ny,nz,lclrchg &
       ,nstp,mpi_md_world,iprint)
  if( lpdens ) call init_pdens(myid_md,hmat,mpi_md_world,iprint)

!.....Add PKA velocity to some atom
  if( pka_energy .gt. 0d0 ) then
    call add_pka_velocity(ntot0,hmat,tagtot,rtot,vtot)
  endif

  call accum_time('overhead',mpi_wtime()-t0)
!.....Call pmd_core to perform MD; all the arguments are in pmdvars module
  call pmd_core(hunit,hmat,ntot0,ntot,tagtot,rtot,vtot,atot,stot, &
       ekitot,epitot,auxtot,epot,ekin,stnsr)

  if( myid_md.eq.0 ) then
    tmp = mpi_wtime()
    if( trim(ciofmt).eq.'bin' .or. trim(ciofmt).eq.'binary' ) then
      call write_pmdtot_bin(20,cpmdfin,ntot,hunit,hmat, &
             tagtot,rtot,vtot)
    elseif( trim(ciofmt).eq.'ascii' ) then
      call write_pmdtot_ascii(20,cpmdfin,ntot,hunit,hmat, &
             tagtot,rtot,vtot)
    endif
    call accum_time('write_xxx',mpi_wtime()-tmp)
  endif

  if( lflux ) call final_lflux(myid_md)
  if( lpdens ) call final_pdens(myid_md,mpi_md_world,hmat)

!.....write energy, forces and stresses only for fitpot
  if( myid_md.eq.0 ) then
    call write_force(21,'.pmd',hmat,epot,ntot,tagtot,atot,stnsr)
    call accum_time('total',mpi_wtime()-t0)
    call report_time(6,iprint)
    call report_mem(6,iprint)
    print *,''
    call time_stamp(' Job finished')
  endif

1 continue
  call mpi_barrier(mpicomm,ierr)
  call mpi_comm_free(mpi_md_world,ierr)
  deallocate(tagtot,rtot,vtot,epitot,ekitot,stot,atot)
  deallocate(auxtot)
  call mpi_finalize(ierr)

end program pmd
!=======================================================================
subroutine write_headline()
!
! Write out headline info for pmd users.
! Assuming that this is called only at 0-th node.
!
  use time, only: time_stamp
  use version
  
  write(6,*) ''
  write(6,'(a)') ' pmd --- Parallel Molecular Dynamics ---'
  write(6,*) ''
  call write_version()
  call write_authors()
  write(6,*) ''
  call time_stamp(' Job started')
  write(6,*) ''

end subroutine write_headline
!=======================================================================
subroutine set_fmv(fmv)
!
! Set default fmv values which might be override
!
  implicit none
  real(8),intent(out):: fmv(3,0:9)

!-----set fmv(1:3,ifmv) to be multiplied to the velocities
  fmv(1:3,0)= (/ 0d0, 0d0, 0d0 /) ! fix
  fmv(1:3,1)= (/ 1d0, 1d0, 1d0 /) ! free move
  fmv(1:3,2)= (/ 1d0, 1d0, 1d0 /) ! free move
  fmv(1:3,3)= (/ 1d0, 1d0, 1d0 /) ! free move
  fmv(1:3,4)= (/ 1d0, 1d0, 1d0 /) ! free move
  fmv(1:3,5)= (/ 1d0, 1d0, 1d0 /) ! free move
  fmv(1:3,6)= (/ 1d0, 1d0, 1d0 /) ! free move
  fmv(1:3,7)= (/ 1d0, 1d0, 1d0 /) ! free move
  fmv(1:3,8)= (/ 1d0, 1d0, 1d0 /) ! free move
  fmv(1:3,9)= (/ 1d0, 1d0, 1d0 /) ! free move

end subroutine set_fmv
!=======================================================================
subroutine write_initial_setting()
  use pmdvars
  use pmdmpi
  use force
  use clrchg,only: lclrchg, cspc_clrchg, clrfield, clr_init
  use localflux,only: lflux,nlx,nly,nlz,noutlflux
  use pdens,only: lpdens,cspc_pdens,npx,npy,npz,orig_pdens,hmat_pdens
  implicit none 
  integer:: i

  write(6,*) ''
  write(6,'(a)') '---------------------------------' &
       //'---------------------------------------'
  write(6,'(a)') '           Initial setting   '
  write(6,'(a)') '---------------------------------' &
       //'---------------------------------------'
  write(6,'(2x,a,5x,i0)')   'num_nodes_x',nx
  write(6,'(2x,a,5x,i0)')   'num_nodes_y',ny
  write(6,'(2x,a,5x,i0)')   'num_nodes_z',nz
  if( nomp.gt.0 ) write(6,'(2x,a,5x,i0)')   'num_omp_threads',nomp
  write(6,'(2x,a)') ''
  write(6,'(2x,a,5x,a)') 'io_format',ciofmt
  write(6,'(2x,a,5x,i0)') 'print_level',iprint
  write(6,'(2x,a)') ''
  write(6,'(2x,a,5x,f0.3)') 'time_interval',dt
  write(6,'(2x,a,5x,i0)')   'num_iteration',nstp
  write(6,'(2x,a,5x,i0)')   'num_out_energy',nerg
  write(6,'(2x,a)') ''
  write(6,'(2x,a,5x,i0)')   'flag_out_pmd',ifpmd
  write(6,'(2x,a,5x,i0)')   'num_out_pmd',npmd
  if( ifpmd.eq.2 ) then  ! if dump output
    write(6,'(2x,a,3x)',advance='no') 'dump_aux_order'
    do i=1,ndumpaux
      write(6,'(x,a)',advance='no') trim(cdumpauxarr(i))
    enddo
    print *,''
  endif
  write(6,'(2x,a)') ''
  write(6,'(2x,a,10(2x,a))') 'force_type   ', &
       (trim(force_list(i)),i=1,num_forces)
  write(6,'(2x,a,5x,f6.3)') 'cutoff_radius',rc
  write(6,'(2x,a,5x,f6.3)') 'cutoff_buffer',rbuf
  write(6,'(2x,a)') ''
  if( cmin.ne.'' ) then
    write(6,'(2x,a,5x,a)') 'minimization',trim(cmin)
  endif
  write(6,'(2x,a,5x,i0)') 'flag_damping',ifdmp
  write(6,'(2x,a,5x,f0.4)') 'damping_coeff',dmp
  write(6,'(2x,a,5x,f0.4)') 'converge_eps',eps_conv
  write(6,'(2x,a,5x,i0)')   'converge_num',n_conv
  write(6,'(2x,a,5x,i0)')   'min_iteration',minstp
  write(6,'(2x,a)') ''
!.....temperature control
  write(6,'(2x,a,5x,a)') 'temperature_control',trim(ctctl)
  write(6,'(2x,a,5x,f8.2)') 'initial_temperature',tinit
  if( trim(ctctl).eq.'Berendsen' .or. &
       trim(ctctl).eq.'Langevin' )  then
    if( tfin.ge.0d0 ) then
      write(6,'(2x,a,5x,f8.2)') 'final_temperature  ',tfin
    else
      do i=1,9
        write(6,'(2x,a,i3,f8.2)') 'temperature_target',i,ttgt(i)
      enddo
    endif
    write(6,'(2x,a,5x,f0.1)') 'temperature_relax_time',trlx
    if( ltdst ) then
!.....temperature distribution
      write(6,'(2x,a)') ''
      write(6,'(2x,a,5x,l2)') 'flag_temp_dist',ltdst
      write(6,'(2x,a,5x,i0)') 'num_temp_dist',ntdst
    endif
  endif
  write(6,'(2x,a,5x,i0)') 'remove_translation',nrmtrans
  write(6,'(2x,a)') ''
!.....pressure control
  write(6,'(2x,a,5x,a)') 'stress_control',trim(cpctl)
  if( trim(cpctl).eq.'Berendsen' .or. &
       trim(cpctl).eq.'vc-Berendsen' ) then
    write(6,'(2x,a)') 'stress_target'
    write(6,'(5x,3es11.3)') stgt(1,1:3)
    write(6,'(5x,3es11.3)') stgt(2,1:3)
    write(6,'(5x,3es11.3)') stgt(3,1:3)
    write(6,'(2x,a)') 'cell_fix'
    write(6,'(5x,3(2x,l))') lcellfix(1,1:3)
    write(6,'(5x,3(2x,l))') lcellfix(2,1:3)
    write(6,'(5x,3(2x,l))') lcellfix(3,1:3)
    write(6,'(2x,a)') ''
    
  else if( trim(cpctl).eq.'vv-Berendsen' ) then
    write(6,'(2x,a,5x,f0.3)') 'pressure_target',ptgt
    write(6,'(2x,a)') 'cell_fix'
    write(6,'(5x,3(2x,l))') lcellfix(1,1:3)
    write(6,'(5x,3(2x,l))') lcellfix(2,1:3)
    write(6,'(5x,3(2x,l))') lcellfix(3,1:3)
  endif
  write(6,*) ''
!.....strain control
  write(6,'(2x,a,5x,a)') 'zload_type',trim(czload_type)
  if( trim(czload_type).ne.'none' ) then
    write(6,'(2x,a,5x,f0.3)') 'zload_skin_width',zskin_width
    write(6,'(2x,a,5x,f0.3)') 'zload_shear_angle',zshear_angle
    write(6,'(2x,a,5x,f0.3)') 'final_strain',strfin
  endif
  write(6,'(2x,a)') ''
!.....velocity multiplying factor
  write(6,'(2x,a)') 'factor_direction'
  do i=0,9
    write(6,'(5x,i3,3f8.3)') i,fmv(1:3,i)
  enddo
  write(6,'(2x,a)') ''
!.....Mass
  write(6,'(2x,a)') 'mass'
  do i=1,nspmax
    if( trim(specorder(i)).ne.'x' ) write(6,'(5x,i3,a4,f10.4)') i,specorder(i),am(i)
  enddo
  write(6,'(2x,a)') ''
!.....Boundary condition
  write(6,'(2x,a,5x,a)') 'boundary',trim(boundary)
  write(6,'(2x,a)') ''
!.....Color charge NEMD
  if( lclrchg ) then
    write(6,'(2x,a,5x,l)') 'flag_clrchg',lclrchg
    write(6,'(2x,a,5x,a)') 'clr_init',trim(clr_init)
    write(6,'(2x,a,5x,a)') 'spcs_clrchg',trim(cspc_clrchg)
    write(6,'(2x,a,3(2x,f0.4))') 'clrfield',clrfield(1:3)
    write(6,'(2x,a)') ''
  endif
!.....Local flux
  if( lflux ) then
    write(6,'(2x,a,5x,l)') 'flag_lflux',lclrchg
    write(6,'(2x,a,2x,i0)') 'num_out_lflux',noutlflux
    write(6,'(2x,a,3(2x,i0))') 'ndiv_lflux',nlx,nly,nlz
    write(6,'(2x,a)') ''
  endif
!.....Probability density
  if( lpdens ) then
    write(6,'(2x,a,5x,l)') 'flag_pdens',lpdens
    write(6,'(2x,a,5x,a)') 'spcs_pdens',trim(cspc_pdens)
    write(6,'(2x,a,3x,3(2x,i0))') 'ndiv_pdens',npx,npy,npz
    write(6,'(2x,a,3x,3(2x,f8.2))') 'orig_pdens',orig_pdens(1:3)
    write(6,'(2x,a)') 'hmat_pdens'
    do i=1,3
      write(6,'(5x,3(2x,f8.2))') hmat_pdens(1:3,i)
    enddo
    write(6,'(2x,a)') ''
  endif
!.....Charge
!!$  write(6,'(2x,a)') 'charge'
!!$  do i=1,nspmax
!!$    if( trim(specorder(i)).ne.'x' ) write(6,'(5x,i3,f10.5)') i,schg(i)
!!$  enddo
!!$  write(6,'(2x,a,5x,a)') 'fix_charge',trim(chgfix)
!!$  write(6,'(2x,a)') ''

  if( lmetaD ) write(6,'(2x,a,5x,l2)') 'metadynamics',lmetaD
  if( lconst ) write(6,'(2x,a,5x,l2)') 'constraints',lconst

  write(6,'(a)') '---------------------------------' &
       //'---------------------------------------'

end subroutine write_initial_setting
!=======================================================================
subroutine write_inpmd(ionum,cfname)
  integer,intent(in):: ionum
  character(len=*),intent(in):: cfname

  character(len=128):: cline

  write(6,*) ''
  write(6,'(a)') '---------------------------------' &
       //'---------------------------------------'
  write(6,'(a)') '                           in.pmd'
  write(6,'(a)') '---------------------------------' &
       //'---------------------------------------'

  open(ionum,file=trim(cfname),status='old')
  do while(.true.)
    read(ionum,'(a)',end=10) cline
    write(6,'(a)') trim(cline)
  enddo
10 close(ionum)
  write(6,'(a)') '---------------------------------' &
       //'---------------------------------------'

end subroutine write_inpmd
!=======================================================================
subroutine bcast_params()
  use pmdvars
  use force
  use extforce,only: lextfrc,cspc_extfrc,extfrc
  use clrchg,only: lclrchg,cspc_clrchg,clr_init,clrfield
  use localflux,only: lflux,nlx,nly,nlz,noutlflux
  use pdens,only: lpdens,npx,npy,npz,cspc_pdens,orig_pdens,hmat_pdens
  implicit none
  include 'mpif.h'

  integer:: ierr

  if( myid_md.eq.0 ) write(6,'(/,a)') ' Broadcast data to be shared' &
       //' with all the nodes.'
!-----Broadcast input parameters to all nodes
  call mpi_bcast(namax,1,mpi_integer,0,mpicomm,ierr)
  call mpi_bcast(nbmax,1,mpi_integer,0,mpicomm,ierr)
  call mpi_bcast(nnmax,1,mpi_integer,0,mpicomm,ierr)
  call mpi_bcast(nstp,1,MPI_INTEGER,0,mpicomm,ierr)
  call mpi_bcast(minstp,1,MPI_INTEGER,0,mpicomm,ierr)
  call mpi_bcast(dt,1,MPI_REAL8,0,mpicomm,ierr)
  call mpi_bcast(vardt_len,1,MPI_REAL8,0,mpicomm,ierr)
  call mpi_bcast(rc,1,MPI_REAL8,0,mpicomm,ierr)
  call mpi_bcast(rc1nn,1,MPI_REAL8,0,mpicomm,ierr)
  call mpi_bcast(rbuf,1,MPI_REAL8,0,mpicomm,ierr)
  call mpi_bcast(cmin,20,mpi_character,0,mpicomm,ierr)
  call mpi_bcast(dmp,1,MPI_REAL8,0,mpicomm,ierr)
  call mpi_bcast(eps_conv,1,MPI_REAL8,0,mpicomm,ierr)
  call mpi_bcast(n_conv,1,mpi_integer,0,mpicomm,ierr)
  call mpi_bcast(tinit,1,mpi_real8,0,mpicomm,ierr)
  call mpi_bcast(tfin,1,mpi_real8,0,mpicomm,ierr)
  call mpi_bcast(ctctl,20,mpi_character,0,mpicomm,ierr)
  call mpi_bcast(ttgt,9,MPI_REAL8,0,mpicomm,ierr)
  call mpi_bcast(trlx,1,MPI_REAL8,0,mpicomm,ierr)
  call mpi_bcast(rseed,1,MPI_REAL8,0,mpicomm,ierr)
  call mpi_bcast(nerg,1,MPI_INTEGER,0,mpicomm,ierr)
  call mpi_bcast(ifpmd,1,MPI_INTEGER,0,mpicomm,ierr)
  call mpi_bcast(npmd,1,MPI_INTEGER,0,mpicomm,ierr)
  call mpi_bcast(ifsort,1,mpi_integer,0,mpicomm,ierr)
  call mpi_bcast(ifdmp,1,MPI_INTEGER,0,mpicomm,ierr)
  call mpi_bcast(iprint,1,mpi_integer,0,mpicomm,ierr)
  call mpi_bcast(fmv,30,mpi_real8,0,mpicomm,ierr)
  call mpi_bcast(shrst,1,mpi_real8,0,mpicomm,ierr)
  call mpi_bcast(cpctl,20,mpi_character,0,mpicomm,ierr)
  call mpi_bcast(ptgt,1,mpi_real8,0,mpicomm,ierr)
  call mpi_bcast(pini,1,mpi_real8,0,mpicomm,ierr)
  call mpi_bcast(pfin,1,mpi_real8,0,mpicomm,ierr)
  call mpi_bcast(srlx,1,mpi_real8,0,mpicomm,ierr)
  call mpi_bcast(stgt,9,mpi_real8,0,mpicomm,ierr)
  call mpi_bcast(lcellfix,9,mpi_logical,0,mpicomm,ierr)
  call mpi_bcast(czload_type,128,mpi_character,0,mpicomm,ierr)
  call mpi_bcast(zskin_width,1,mpi_real8,0,mpicomm,ierr)
  call mpi_bcast(zshear_angle,1,mpi_real8,0,mpicomm,ierr)
  call mpi_bcast(strfin,1,mpi_real8,0,mpicomm,ierr)
  call mpi_bcast(am,nspmax,mpi_real8,0,mpicomm,ierr)
  call mpi_bcast(specorder,nspmax*3,mpi_character,0,mpicomm,ierr)
  call mpi_bcast(ciofmt,6,mpi_character,0,mpicomm,ierr)
  call mpi_bcast(nrmtrans,6,mpi_integer,0,mpicomm,ierr)
  call mpi_bcast(lstrs0,1,mpi_logical,0,mpicomm,ierr)
  call mpi_bcast(boundary,3,mpi_character,0,mpicomm,ierr)
  call mpi_bcast(pka_energy,1,mpi_real8,0,mpicomm,ierr)
  call mpi_bcast(nomp,1,mpi_integer,0,mpicomm,ierr)

!.....Charge related
  call mpi_bcast(chgfix,20,mpi_character,0,mpicomm,ierr)
!.....Force-fields
  call mpi_bcast(cforce,20,mpi_character,0,mpicomm,ierr)
  call bcast_force(mpicomm)
!!$  call mpi_bcast(num_forces,1,mpi_integer,0,mpicomm,ierr)
!!$  if( num_forces.eq.0 ) then
!!$    if( myid_md.eq.0 ) write(6,'(a)') &
!!$         ' Error: no force-field specified'
!!$    call mpi_finalize(ierr)
!!$    stop
!!$  endif
!!$  if( myid_md.ne.0 ) then
!!$    allocate(force_list(num_forces))
!!$  endif
!!$  call mpi_bcast(force_list,128*num_forces,mpi_character &
!!$       ,0,mpicomm,ierr)
  call mpi_bcast(ifcoulomb,1,mpi_integer,0,mpicomm,ierr)
  call mpi_bcast(schg,nspmax,mpi_real8,0,mpicomm,ierr)
!.....NEMD
  call mpi_bcast(ltdst,1,mpi_logical,0,mpicomm,ierr)
  if( ltdst ) then
    call mpi_bcast(ntdst,1,mpi_integer,0,mpicomm,ierr)
  endif
!.....Color charge NEMD
  call mpi_bcast(lclrchg,1,mpi_logical,0,mpicomm,ierr)
  if( lclrchg ) then
    call mpi_bcast(cspc_clrchg,3,mpi_character,0,mpicomm,ierr)
    call mpi_bcast(clr_init,20,mpi_character,0,mpicomm,ierr)
    call mpi_bcast(clrfield,3,mpi_real8,0,mpicomm,ierr)
  endif
!.....Local flux
  call mpi_bcast(lflux,1,mpi_logical,0,mpicomm,ierr)
  if( lflux ) then
    call mpi_bcast(noutlflux,1,mpi_integer,0,mpicomm,ierr)
    call mpi_bcast(nlx,1,mpi_integer,0,mpicomm,ierr)
    call mpi_bcast(nly,1,mpi_integer,0,mpicomm,ierr)
    call mpi_bcast(nlz,1,mpi_integer,0,mpicomm,ierr)
  endif
!.....External force
  call mpi_bcast(lextfrc,1,mpi_logical,0,mpicomm,ierr)
  if( lextfrc ) then
    call mpi_bcast(cspc_extfrc,3,mpi_character,0,mpicomm,ierr)
    call mpi_bcast(extfrc,3,mpi_real8,0,mpicomm,ierr)
  endif
!.....Probability density
  call mpi_bcast(lpdens,1,mpi_logical,0,mpicomm,ierr)
  if( lpdens ) then
    call mpi_bcast(cspc_pdens,3,mpi_character,0,mpicomm,ierr)
    call mpi_bcast(npx,1,mpi_integer,0,mpicomm,ierr)
    call mpi_bcast(npy,1,mpi_integer,0,mpicomm,ierr)
    call mpi_bcast(npz,1,mpi_integer,0,mpicomm,ierr)
    call mpi_bcast(orig_pdens,3,mpi_real8,0,mpicomm,ierr)
    call mpi_bcast(hmat_pdens,3*3,mpi_real8,0,mpicomm,ierr)
  endif
!.....Metadynamics
  call mpi_bcast(lmetaD,1,mpi_logical,0,mpicomm,ierr)
!.....Constratins
  call mpi_bcast(lconst,1,mpi_logical,0,mpicomm,ierr)
!.....Reduced force
  call mpi_bcast(lrdcfrc,1,mpi_logical,0,mpicomm,ierr)
!.....Linked-cell reordering
  call mpi_bcast(lreorder,1,mpi_logical,0,mpicomm,ierr)
!.....Deformation
  call mpi_bcast(cdeform,128,mpi_character,0,mpicomm,ierr)
  call mpi_bcast(dhratio,9,mpi_real8,0,mpicomm,ierr)
!.....Structure analysis
  call mpi_bcast(cstruct,128,mpi_character,0,mpicomm,ierr)
  call mpi_bcast(istruct,1,mpi_integer,0,mpicomm,ierr)

end subroutine bcast_params
!=======================================================================
subroutine write_force(ionum,cpostfix,h,epot,ntot,tagtot,atot,stnsr)
  use util,only: itotOf
  implicit none
  include "./params_unit.h"
  integer,intent(in):: ionum,ntot
  character(len=*),intent(in):: cpostfix
  real(8),intent(in):: h(3,3),epot,stnsr(3,3)
  real(8),intent(inout):: tagtot(ntot),atot(3,ntot)

  integer:: i,n0,ixyz,ierr
  real(8):: at(3),ptmp(6)
  integer,parameter:: nmpi = 2
!!$  integer,external:: itotOf


!.....Write out forces
  open(ionum,file='frc'//trim(cpostfix),status='replace')
  write(ionum,'(i10)') ntot
  do i=1,ntot
    write(ionum,'(3es15.7,i8)') atot(1:3,i),itotOf(tagtot(i))
  enddo
  close(ionum)

!.....Write out energy
  open(ionum+1,file='erg'//trim(cpostfix),status='replace')
  write(ionum+1,'(es23.14e3)') epot
  close(ionum+1)

!.....Write stress tensor, negative as compressive, positive as tensile
  ptmp(1) = -stnsr(1,1)
  ptmp(2) = -stnsr(2,2)
  ptmp(3) = -stnsr(3,3)
  ptmp(4) = -stnsr(3,2)
  ptmp(5) = -stnsr(1,3)
  ptmp(6) = -stnsr(1,2)
  open(ionum+2,file='strs'//trim(cpostfix),status='replace')
  write(ionum+2,'(6f12.4)') ptmp(1:6)
  close(ionum+2)

end subroutine write_force
!=======================================================================
subroutine set_atomic_charges(ntot,chg,tag,nspmax,chgfix,schg,myid,iprint)
  implicit none 
  integer,intent(in):: ntot,nspmax,myid,iprint
  real(8),intent(in):: tag(ntot),schg(nspmax)
  character(len=*),intent(in):: chgfix
  real(8),intent(out):: chg(ntot)

  integer:: i,is

  if( trim(chgfix).eq.'input' ) then
    if( myid.eq.0 .and. iprint.ne.1 ) then
      print *,'Charges are set from input.'
    endif
    do i=1,ntot
      is = int(tag(i))
      chg(i)= schg(is)
    enddo
  endif

end subroutine set_atomic_charges
!=======================================================================
subroutine check_ensemble()
  use pmdvars
  implicit none

  logical:: l_temp
  character:: c_2*1,c_3*1,c_ensemble*3

  c_2 = 'V'
  c_3 = 'E'

  l_temp = .false.

  if(  trim(ctctl).eq.'Langevin' .or. &
       trim(ctctl).eq.'Berendsen' .or. &
       trim(ctctl).eq.'ttm' ) then
    l_temp = .true.
    c_3 = 'T'
  endif

  if(  trim(cpctl).eq.'Berendsen' .or. &
       trim(cpctl).eq.'vc-Berendsen' ) then
    c_2 = 'p'
    if( .not. l_temp ) c_3 = 'H'
  else if( trim(cpctl).eq.'vv-Berendsen' ) then
    c_2 = 'P'
    if( .not. l_temp ) c_3 = 'H'
  endif

  c_ensemble = 'N'//c_2//c_3
  print *,''
  print *,'Ensemble = ',trim(c_ensemble)

end subroutine check_ensemble
!=======================================================================
subroutine add_pka_velocity(ntot0,hmat,tagtot,rtot,vtot)
! 
! Add PKA velocity to some atom from a given PKA energy
! 
  use pmdvars
  use random,only: urnd
  implicit none
  include './params_unit.h'
  integer,intent(in):: ntot0
  real(8),intent(in):: tagtot(ntot0),rtot(3,ntot0),hmat(3,3,0:1)
  real(8),intent(inout):: vtot(3,ntot0)

  integer:: i,icntr,is
  real(8):: dmin,d,theta,phi,rx,ry,rz,vel,vx,vy,vz

  if( myid_md.eq.0 ) then
    if( iatom_pka.le.0 ) then
!.....Find a center atom ICNTR
      icntr = 0
      dmin = 1d+10
      do i=1,ntot0
        d = (rtot(1,i)-0.5d0)**2 &
             +(rtot(2,i)-0.5d0)**2 &
             +(rtot(3,i)-0.5d0)**2
        if( d .lt. dmin ) then
          icntr = i
          dmin = d
        endif
      enddo
    else
      if( iatom_pka.gt.ntot0 ) then
        print *,'ERROR: iatom_pka.gt.ntot0 !'
        stop
      endif
      icntr = iatom_pka
    endif

!.....Add PKA velocity at ICNTR atom
!.....Get random theta and phi
    theta = 90d0 *urnd() /180d0 *pi
    phi = 90d0 *urnd() /180d0 *pi
    rx = sin(theta)*cos(phi)
    ry = sin(theta)*sin(phi)
    rz = cos(theta)
    is = int(tagtot(icntr))
!.....[eV] to [Ang/fs]
    vel = sqrt(pka_energy*ev2j *2d0 /(am(is)*ump2kg)) *m2ang /s2fs
    vx = rx*vel
    vy = ry*vel
    vz = rz*vel
    print *,''
    print '(a)',' Primary knock-on atom: '
    print '(a,i0,2x,3es12.4)',      '   atom-id,rtot = ',icntr,rtot(1:3,icntr)
    print '(a,es12.4,a)','   energy  = ',pka_energy,' eV'
    print '(a,2f9.3)',   '   theta,phi = ',theta,phi
    print '(a,4es12.4)', '   vel,vx,vy,vz = ',vel,vx,vy,vz
    
!.....Assume an orthogonal simulation box
    vx = vx /hmat(1,1,0)
    vy = vy /hmat(2,2,0)
    vz = vz /hmat(3,3,0)
    print '(a,3es12.4)','   vtot before = ',vtot(1:3,icntr)
    vtot(1,icntr) = vtot(1,icntr) +vx
    vtot(2,icntr) = vtot(2,icntr) +vy
    vtot(3,icntr) = vtot(3,icntr) +vz
    print '(a,3es12.4)','   vtot after  = ',vtot(1:3,icntr)

  endif
end subroutine add_pka_velocity
!=======================================================================
subroutine determine_division(h,myid,nnode,rc,nx,ny,nz,iprint)
!
!  Determine parallel division along a1,a2,a3.
!
  use vector,only: dot
  implicit none
  include "./const.h"
  integer,intent(in):: nnode,myid,iprint
  real(8),intent(in):: h(3,3),rc
  integer,intent(inout):: nx,ny,nz

  integer:: imax,nd(3),ndnew(3),nnew,iorder(3),i,j,k,l,m,n,ncycle
  real(8):: al0(3),al(3)

!.....If serial run, NX,NY,NZ should all be 1.
  if( nnode.eq.1 ) then
    nd(1) = 1
    nd(2) = 1
    nd(3) = 1 
    goto 10
  endif

!.....Lengths of each axis
  al0(1) = sqrt(dot(h(1:3,1),h(1:3,1)))
  al0(2) = sqrt(dot(h(1:3,2),h(1:3,2)))
  al0(3) = sqrt(dot(h(1:3,3),h(1:3,3)))

  nd(1) = 1
  nd(2) = 1
  nd(3) = 1 
  al(1:3) = al0(1:3) /nd(1:3)
  if( al(1).lt.rc .and. al(2).lt.rc .and. al(3).lt.rc ) then
    goto 10
  endif

!.....Loop until total num of division exceeds num of MPI parallel nodes
  do
!.....Increase num of division of the axis of which divided length is the longest
    iorder(1:3) = (/1,2,3/)
    do i=2,1,-1
      do j=1,i
        l = iorder(j)
        m = iorder(j+1)
        if( al(l).lt.al(m) ) then
          iorder(j) = m
          iorder(j+1) = l
        endif
      enddo
    enddo

    ncycle = 0
    do n=1,3
      imax = iorder(n)
      ndnew(1:3) = nd(1:3)
      ndnew(imax) = ndnew(imax) +1

!.....New total num of divisions
      nnew = ndnew(1)*ndnew(2)*ndnew(3)
      al(1:3) = al0(1:3) /ndnew(1:3)

!.....Check whether the total num exceeds num of MPI parallel
      if( nnew.gt.nnode ) then
        ncycle = ncycle +1
        cycle
      else if( nnew.eq.nnode ) then
        nd(1:3) = ndnew(1:3)  ! use current ndnew(:)
        goto 10
!.....Check whether the minimum divided length is shorter than cutoff radius         
      else if( (ndnew(1).ne.1.and.al(1).lt.rc) .or. &
         (ndnew(2).ne.1.and.al(2).lt.rc) .or. &
         (ndnew(3).ne.1.and.al(3).lt.rc) ) then
        goto 10  ! use previous nd(:)
      else
        exit  ! exit do n=1,3
      endif
    enddo

    if( ncycle.ge.3 ) goto 10
    nd(1:3) = ndnew(1:3)
  enddo

10 continue
  nx = nd(1)
  ny = nd(2)
  nz = nd(3)

 if( iprint.ge.ipl_basic ) then
    print '(a,3(1x,i0))',' Number of spatial divisions ' &
        //'automatically set, NX,NY,NZ=',nx,ny,nz
  endif
  return

end subroutine determine_division
!=======================================================================
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd lib"
!     End:
