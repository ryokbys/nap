program pmd
!-----------------------------------------------------------------------
!                     Last-modified: <2019-07-24 11:14:02 Ryo KOBAYASHI>
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
  use pmdio
  use pmdmpi
  use version
  use force
  use Coulomb, only: cterms
  use util, only: time_stamp, itotOf
  use element
  implicit none
  include "mpif.h"
  include "./params_unit.h"

#ifdef __DISL__
!.....Epot threshold for disl core extraction [Hartree]
  real(8),parameter:: epith = -0.1410d0
#endif

  integer:: i,j,k,l,m,n,ia,ib,is,ifmv,nave,nspl,i_conv,nstp_done
  integer:: mpicolor,mpikey,ierr,jerr,itmp,nprocs,nnmax_est
  real(8):: tmp,hscl(3),aai(3),ami,dt2,tave,vi(3),vl(3),epot,ekin,rmin
  character(len=3):: csp
  type(atom):: elem
!!$  integer,external:: itotOf

!-----initialize the MPI environment
  call mpi_init(ierr)
!-----total number of MD-nodes
  call mpi_comm_size(MPI_COMM_WORLD,nprocs,ierr)
!-----my rank in MD-nodes
  call mpi_comm_rank(MPI_COMM_WORLD,myid_md,ierr)
  call mpi_comm_dup(MPI_COMM_WORLD,mpicomm,ierr)
  mpi_md_world = mpicomm

  call init_element()

!.....Set fmv as default value before reading 'in.pmd'
  call set_fmv(fmv)

  if( myid_md.eq.0 ) then
    write(6,'(a)') '=================================' &
         //'======================================='
    write(6,'(a)') ' PMD: A Parallel Molecular Dynamics program '
    write(6,*) ''
    call write_version()
    call write_authors()
    write(6,'(a)') '=================================' &
         //'======================================='
    write(6,*) ''
    call time_stamp(' Job started')
    write(6,*) ''

!.....Read atom configuration file 1st
    if( trim(ciofmt).eq.'bin' .or. trim(ciofmt).eq.'binary' ) then
      write(6,*) 'Read pmdini in binary mode.'
      call read_pmdtot_bin(20,trim(cpmdini))
    else if( trim(ciofmt).eq.'ascii' ) then
      write(6,*) 'Read pmdini in ascii mode.'
      call read_pmdtot_ascii(20,trim(cpmdini))
    else
      write(6,*) 'Error: io_format must be either ascii, ' &
           //'bin or binary.'
      stop
    endif

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
    call spcs_info()

    write(6,*) ''
    write(6,'(a,i0)') ' Number of processes in MPI = ',nprocs
!.....Read in.pmd after reading the atom configuration file.
    call read_input(10,trim(cinpmd))
    call check_cmin(cmin,ifdmp)
    if( ifpmd.eq.2 ) then ! if dump output
      call make_cdumpauxarr()
    endif
    call write_initial_setting()
!        call write_inpmd(10,trim(cinpmd))
    if( num_forces.eq.0 ) stop ' ERROR: no force-field specified'
    if( trim(ctctl).eq.'ttm' ) then
      print *,''
      print *,'Since the two-temperature model (TTM) MD...'
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
  endif

  call bcast_params()


  if( myid_md.eq.0 ) then
    allocate(chgtot(ntot0),chitot(ntot0))
    chitot(1:ntot0) = 0d0
    chgtot(1:ntot0) = 0d0
!!$    call set_atomic_charges(ntot0,chgtot,tagtot,nspmax &
!!$         ,chgfix,schg,myid_md,iprint)

!.....Determine nx,ny,nz using rc and hmat info
    if( .not. (nx.gt.0 .and. ny.gt.0 .and. nz.gt.0 ) ) then
      call determine_division(h,myid_md,nprocs,rc,nx,ny,nz,iprint)
    endif

!.....Make ntot and ?tot() not null in nodes myid_md != 0
  else
    ntot0 = 1
    allocate(tagtot(ntot0),rtot(3,ntot0),vtot(3,ntot0),epitot(ntot0) &
         ,ekitot(3,3,ntot0),stot(3,3,ntot0),atot(3,ntot0) &
         ,chgtot(ntot0),chitot(ntot0))
    chitot(1:ntot0) = 0d0
    chgtot(1:ntot0) = 0d0
  endif

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
    if( myid_md.eq.0 .and. iprint.gt.0 ) then
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

!.....Add PKA velocity to some atom
  if( pka_energy .gt. 0d0 ) then
    call add_pka_velocity(myid_md)
  endif

!.....call pmd_core to perfom MD
  call pmd_core(hunit,h,ntot0,tagtot,rtot,vtot,atot,stot &
       ,ekitot,epitot,chgtot,chitot,nstp,nerg,npmd &
       ,myid_md,mpi_md_world,nodes_md,nx,ny,nz,specorder &
       ,am,dt,vardt_len,ciofmt,ifpmd,rc,rbuf,rc1nn,ifdmp,dmp &
       ,minstp,tinit,tfin,ctctl,ttgt,trlx,ltdst,ntdst,nrmtrans,cpctl &
       ,stgt,ptgt,pini,pfin,srlx,stbeta,strfin,lstrs0,lcellfix,fmv &
       ,stnsr,epot,ekin,n_conv,ifcoulomb,czload_type,zskin_width &
       ,zshear_angle,eps_conv,ifsort,iprint,nstp_done,lvc,boundary &
       ,lmetaD,lconst,lrdcfrc,cstruct,istruct,cdeform,dhratio)

  if( myid_md.eq.0 ) then
    if( trim(ciofmt).eq.'bin' .or. trim(ciofmt).eq.'binary' ) then
      call write_pmdtot_bin(20,cpmdfin)
    elseif( trim(ciofmt).eq.'ascii' ) then
      call write_pmdtot_ascii(20,cpmdfin)
    endif
  endif

!.....write energy, forces and stresses only for fitpot
  if( myid_md.eq.0 ) then
    call write_force(21,'.pmd',h,epot,ntot,tagtot,atot,stnsr)
    print *,''
    call time_stamp(' Job finished')
  endif

1 continue
  call mpi_barrier(mpicomm,ierr)
  call mpi_comm_free(mpi_md_world,ierr)
  deallocate(tagtot,rtot,vtot,epitot,ekitot,stot,atot)
  call mpi_finalize(ierr)

end program pmd
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
  fmv(1:3,2)= (/ 1d0, 1d0, 0d0 /) ! xy-only
  fmv(1:3,3)= (/ 1d0, 1d0, 1d0 /) ! free move
  fmv(1:3,4)= (/ 1d0, 1d0, 1d0 /) ! free move
  fmv(1:3,5)= (/ 1d0, 1d0, 1d0 /) ! free move
  fmv(1:3,6)= (/ 1d0, 1d0, 1d0 /) ! free move
  fmv(1:3,7)= (/ 1d0, 1d0, 1d0 /) ! free move
  fmv(1:3,8)= (/ 1d0, 1d0, 1d0 /) ! free move
  fmv(1:3,9)= (/ 1d0, 1d0, 1d0 /) ! free move

end subroutine set_fmv
!=======================================================================
subroutine check_cmin(cmin,ifdmp)
  implicit none
  character,intent(in):: cmin 
  integer,intent(inout):: ifdmp

  if( cmin.ne.'' ) then
    if( cmin.eq.'none' ) then
      ifdmp = 0
    else if( cmin.eq.'damp' ) then
      ifdmp = 1
    else if( cmin.eq.'FIRE' ) then
      ifdmp = 2
    else
      write(6,'(a)') ' [Warning] There is no minimization' &
           //' method: '//cmin
      write(6,'(a)') '           So ifdmp is set 0.'
      ifdmp = 0
    endif
  endif

end subroutine check_cmin
!=======================================================================
subroutine write_initial_setting()
  use pmdio
  use pmdmpi
  use force
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
  if( lvc ) then
    write(6,'(2x,a,5x,l2)') 'charge_optimize',lvc
  endif
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
  write(6,'(2x,a,5x,f0.4)') 'initial_temperature',tinit
  write(6,'(2x,a,5x,a)') 'temperature_control',trim(ctctl)
  if( trim(ctctl).eq.'Berendsen' .or. &
       trim(ctctl).eq.'Langevin' )  then
    if( tfin.ge.0d0 ) then
      write(6,'(2x,a,5x,f0.4)') 'final_temperature',tfin
    else
      do i=1,9
        write(6,'(2x,a,i3,f8.1)') 'temperature_target',i,ttgt(i)
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
!.....Charge
!!$  write(6,'(2x,a)') 'charge'
!!$  do i=1,nspmax
!!$    if( trim(specorder(i)).ne.'x' ) write(6,'(5x,i3,f10.5)') i,schg(i)
!!$  enddo
!!$  write(6,'(2x,a,5x,a)') 'fix_charge',trim(chgfix)
!!$  write(6,'(2x,a,5x,l2)') 'variable_charge',lvc
  write(6,'(2x,a)') ''

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
  use pmdio
  use pmdmpi
  use force
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
!.....Charge related
  call mpi_bcast(lvc,1,mpi_logical,0,mpicomm,ierr)
  call mpi_bcast(chgfix,20,mpi_character,0,mpicomm,ierr)
!.....Force-fields
  call mpi_bcast(cforce,20,mpi_character,0,mpicomm,ierr)
  call bcast_force()
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
!.....Metadynamics
  call mpi_bcast(lmetaD,1,mpi_logical,0,mpicomm,ierr)
!.....Constratins
  call mpi_bcast(lconst,1,mpi_logical,0,mpicomm,ierr)
!.....Reduced force
  call mpi_bcast(lrdcfrc,1,mpi_logical,0,mpicomm,ierr)
!.....Deformation
  call mpi_bcast(cdeform,128,mpi_character,0,mpicomm,ierr)
  call mpi_bcast(dhratio,9,mpi_real8,0,mpicomm,ierr)
!.....Structure analysis
  call mpi_bcast(cstruct,128,mpi_character,0,mpicomm,ierr)
  call mpi_bcast(istruct,1,mpi_integer,0,mpicomm,ierr)

end subroutine bcast_params
!=======================================================================
subroutine write_force(ionum,cpostfix,h,epot,ntot &
     ,tagtot,atot,stnsr)
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
  ptmp(1) = stnsr(1,1)*up2gpa*(-1d0)
  ptmp(2) = stnsr(2,2)*up2gpa*(-1d0)
  ptmp(3) = stnsr(3,3)*up2gpa*(-1d0)
  ptmp(4) = stnsr(3,2)*up2gpa*(-1d0)
  ptmp(5) = stnsr(1,3)*up2gpa*(-1d0)
  ptmp(6) = stnsr(1,2)*up2gpa*(-1d0)
  open(ionum+2,file='strs'//trim(cpostfix),status='replace')
  write(ionum+2,'(6f12.4)') ptmp(1:6)
  close(ionum+2)

end subroutine write_force
!=======================================================================
subroutine set_atomic_charges(ntot,chg,tag,nspmax &
     ,chgfix,schg,myid,iprint)
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
  use pmdio
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
subroutine add_pka_velocity(myid_md)
! 
! Add PKA velocity to some atom from a given PKA energy
! 
  use pmdio
  implicit none
  include './params_unit.h'
  integer,intent(in):: myid_md

  integer:: i,icntr,is
  real(8):: dmin,d,theta,phi,rx,ry,rz,vel,vx,vy,vz
  interface
    function urnd(dseed0)
      real(8),intent(in),optional:: dseed0
      real(8):: urnd
    end function urnd
  end interface
!      real(8),external:: urnd

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
    vx = vx /h(1,1,0)
    vy = vy /h(2,2,0)
    vz = vz /h(3,3,0)
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

 if( iprint.gt.0 ) then
    print '(a,3(1x,i0))',' Number of spatial divisions ' &
        //'automatically set, NX,NY,NZ=',nx,ny,nz
  endif
  return

end subroutine determine_division
!=======================================================================
subroutine spcs_info()
  use pmdio

  integer:: nsps(nspmax),ia,is

  nsps(:) = 0
  do ia=1,ntot
    is = int(tagtot(ia))
    nsps(is) = nsps(is) +1
  enddo
  write(6,'(a)') ' Number of each species in the initial configuration'
  do is=1,nspmax
    if( trim(specorder(is)).eq.'x' ) cycle
    write(6,'(a,i0)') '   '//specorder(is)//':  ',nsps(is)
  enddo
  
end subroutine spcs_info
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
