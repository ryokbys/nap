module pmdvars
!-----------------------------------------------------------------------
!                    Last modified: <2024-07-09 15:42:03 KOBAYASHI Ryo>
!-----------------------------------------------------------------------
  implicit none
!=======================================================================
! Parameters or constants
!=======================================================================

!.....output #
  integer,parameter:: ioerg = 11
  integer,parameter:: iotemp= 12
  integer,parameter:: iostrs= 13
  integer,parameter:: iotdst= 14
  integer,parameter:: iozload= 15

  character(len=20),parameter:: cinpmd='in.pmd'

  character(len=20),parameter:: cpmdini = 'pmdini'
  character(len=20),parameter:: cpmdfin = 'pmdfin'

!.....Max num of category
!.....In pmd, any categorical value is limited from 0 to 9
  integer,parameter:: maxcat = 9
!.....Max num of species
  integer,parameter:: nspmax = 9
!.....Max num of temperatures
  integer,parameter:: maxntemps = 9

!-----------------------------------------------------------------------
! Global variables to be given by input file or wrapper function
!-----------------------------------------------------------------------

!!$!.....Data of total system
!!$  real(8):: hunit,h(3,3,0:1)
!!$  integer:: ntot0,ntot,naux
!!$  real(8),allocatable:: rtot(:,:),vtot(:,:),stot(:,:,:),epitot(:) &
!!$       ,ekitot(:,:,:),tagtot(:),atot(:,:)
!!$  real(8),allocatable:: auxtot(:,:)


!.....max. num. of atoms in a node
  integer:: namax = 200
!.....max. num. of boundary-particles
  integer:: nbmax = 100
!.....max. num. of neighbors
  integer:: nnmax = 50

!.....MPI variables 
  integer:: mpicomm
  integer:: myid_md,nodes_md,mpi_md_world
!.....NX,NY,NZ: if either one of these is negative, they are automatically
!     determined using rc and hmat information.
  integer:: nx = -1
  integer:: ny = -1
  integer:: nz = -1
  integer:: nxyz
!.....OpenMP
  integer:: nomp = -1

  integer:: nstp = 0
  integer:: minstp = 0
  integer:: nerg = 1000
  integer:: ifpmd= 2    ! 0:none, 1:pmd, 2:dump (default)
  integer:: npmd = 10
  logical:: lcomb_pos = .false.
  integer:: ifsort= 1
  real(8):: dt = 1d0
  real(8):: vardt_len = 0.1d0  ! Length criterion for variable time-step
  real(8):: rc = 5.0d0
  real(8):: rc1nn = 2.5d0
  real(8):: rbuf= 0d0
  integer:: ifdmp= 0 ! 0:none, 1:damped-MD, 2:FIRE
  character(len=20):: cmin= ''
  real(8):: dmp  = 0.9d0
  real(8):: eps_conv = 1d-4
  integer:: n_conv = 1
!.....Thermostat
  character(len=20):: ctctl='none'
  integer:: iftctl= 0
  real(8):: tinit= -1d0
  real(8):: tfin = -1d0
  real(8):: ttgt(maxntemps)
  data ttgt / 300d0, 300d0, 300d0, 300d0, 300d0, 300d0, &
       300d0, 300d0, 300d0 /
  real(8):: trlx = 100d0
  real(8):: tlimit = 1.0d+5
!.....Multiple temperatures
  logical:: lmultemps = .false.
  integer:: ntemps = 1
  real(8):: tfac(maxntemps),temps(maxntemps),eks(maxntemps)
!.....Random seed
!  If positve, use (RSEED+MYID) as the seed for each process
!  If negative, use the same random seeds for all the parallel process
  real(8):: rseed = 12345d0
!.....Remove translational motion:
!     N< 0: not to remove translation
!     N==0: remove translation only at the beginning
!     N> 1: remove translation at the begining and every N step.
  integer:: nrmtrans = 0
!.....temperature distribution on x-direction
  logical:: ltdst= .false.
  integer:: ntdst= 1
!.....shear stress
  real(8):: shrst = 0.0d0
!.....factors on each moving direction
  real(8):: fmv(3,0:9)
  data fmv &
       / 0d0, 0d0, 0d0, & ! 0
       1d0, 1d0, 1d0, & ! 1
       1d0, 1d0, 1d0, & ! 2
       1d0, 1d0, 1d0, & ! 3
       1d0, 1d0, 1d0, & ! 4
       1d0, 1d0, 1d0, & ! 5
       1d0, 1d0, 1d0, & ! 6
       1d0, 1d0, 1d0, & ! 7
       1d0, 1d0, 1d0, & ! 8
       1d0, 1d0, 1d0  & ! 9
       /
!.....Whether compute stress or not
  logical:: lstrs0 = .true.
  logical:: lstrs = .false.
!.....Barostat: (vv|vc)-Berendsen, (vv|vc)-Langevin
!.....  vv --- variable volume
!.....  vc --- variable cell
  character(len=20):: cpctl='none'
  real(8):: ptgt   = 0d0   ! target pressure [GPa]
  real(8):: pini   = 0d0
  real(8):: pfin   = 0d0
  real(8):: srlx   = 100d0  ! relaxation time [fs]
  real(8):: stbeta = 1d-2   ! 1/B where B is the bulk modulus [GPa]
  real(8):: strfin = 0.0d0
  real(8):: stgt(1:3,1:3)= 0d0  ! target stress tensor [GPa]
  logical:: lcellfix(1:3,1:3)= .false.
  logical:: lhydrostatic = .false.  ! enforce only hydrostatic pressure
!.....charge optimize or variable charge
  logical:: lvc = .false.
!.....Charge setting
  character(len=20):: chgfix='input'
!.....print level
!  0:quiet, 1:normal,
!  >10:fitpot data
  integer:: iprint= 1
  
  character(len=6):: ciofmt='ascii '
  character(len=20):: cforce='none'

!.....Auxiliary data order for dump output
  integer:: naux
  logical:: ldumpaux_changed = .false.
  character(len=256):: cdumpaux = 'ekin epot sxx syy szz syz sxz sxy'
  integer:: ndumpaux
  character(len=6),allocatable:: cdumpauxarr(:)
!.....Auxiliary data
  character(len=6),allocatable:: cauxarr(:)
!.....Indices of aux array
  integer:: iaux_chg,iaux_q,iaux_vq,iaux_tei,iaux_clr,iaux_edsp

!.....mass
  real(8):: am(1:nspmax)= 12.0d0
!.....charges
  real(8):: schg(1:nspmax)= 0d0
!.....species name
  character(len=3):: specorder(nspmax) = 'x'
  logical:: has_specorder = .false.
!.....forces in in.pmd [default: false]
  logical:: has_forces = .false.

!.....Boundary condition: p = periodic, f = free, w = wall
  character(len=3):: boundary = 'ppp'

!.....nnmax update ratio
  real(8):: ratio_nnmax_update = 1.1d0

!-----------------------------------------------------------------------
!  Global variables used in pmd
!-----------------------------------------------------------------------

  integer:: nstp_done
  integer:: nouterg,noutpmd,istp,iocntpmd,iocnterg
  integer:: natm,nb,nsp,nalmax,ntot
  integer:: maxnn = 0
  integer:: maxnb = 0
  real(8):: tcpu,tcpu0,tcpu1,tcpu2,tlspr
  real(8):: epot0,vmaxold,vmax,simtime
  real(8):: tgmm,cgmm,cmass
  real(8):: ediff(nspmax),ediff0(nspmax)
  integer:: ndof(nspmax)
  real(8):: ttgt_lang
!!$  integer,allocatable:: ndof(:)
  integer:: nxmlt
!.....Search time and expiration time
  real(8):: ts,te
  integer:: istpe
!.....simulation box
  real(8):: h(3,3,0:1)
!$acc declare create(h)
  real(8):: hi(3,3),vol,sgm(3,3),al(3),avol
  real(8):: ht(3,3,0:1),hti(3,3),dh
!.....positions, velocities, and accelerations
  real(8),allocatable:: ra(:,:),va(:,:),aa(:,:),ra0(:,:),strs(:,:,:),stt(:,:,:)
!$acc declare create(ra)
!.....real*8 identifier which includes species, index of FMV, total id
  real(8),allocatable:: tag(:)
  integer,allocatable:: lspr(:,:)
!$acc declare create(tag,lspr)
!.....potential and kinetic energy per atoms
  real(8),allocatable:: epi(:),eki(:,:,:),stp(:,:,:)
!$acc declare create(epi)
!.....mass, prefactors
  real(8),allocatable:: acon(:),fack(:)
!.....Factors for ekin
  real(8):: fekin(nspmax),fa2v(nspmax)
!.....atomic strain
!!$  real(8),allocatable:: stn(:,:,:)
!.....Auxiliary data
  real(8),allocatable:: aux(:,:)

  logical:: lcell_updated = .true.

!.....Reallocation
  logical:: lrealloc = .true.

!.....Shear stress
  real(8):: shrfx

!.....Barostat
  real(8):: phyd,ah(3,3),aht(3,3),g(3,3,0:1),gt(3,3,0:1),gi(3,3),gg(3,3)

!.....FIRE parameters
  integer:: nmin_fire = 5
  real(8):: finc_fire = 1.1
  real(8):: fdec_fire = 0.5
  real(8):: alp0_fire = 0.1
  real(8):: falp_fire = 0.99
  real(8):: dtmax_fire = 10.0
! factor to be multiplied to dt to get dtmax_fire
  real(8):: dtmfctr_fire = 10.0

!.....temperature distribution along x
  real(8),allocatable:: tdst(:)
  integer,allocatable:: nadst(:)

!.....space decomposition
  integer,allocatable:: lsb(:,:),lsex(:,:)
  integer:: nn(6),myparity(3),lsrc(6),nex(3),myx,myy,myz
  real(8):: sv(3,6),sorg(3),anxi,anyi,anzi

!.....Number of ifmv
  integer:: nfmv = 1
!.....Number of groups (1-4), but each group takes the value (0-9)
  integer,parameter:: max_group = 4
!  integer:: ngrp = 1

!.....Variable time-step
  logical:: lvardt = .false.

!.....Metadynamics
  logical:: lmetaD = .false. 
!.....Constraints
  logical:: lconst = .false. 
!.....Reduced force
  logical:: lrdcfrc = .false. 
!.....Sort arrays according to linked cell; BUGS with charges
  logical:: lsrt_arrs = .false. 

!.....zload type: zload or shear
  character(len=128):: czload_type= 'none'
!.....top and bottom skin width in which atoms are fixed and/or controlled
  real(8):: zskin_width = 5.0d0
!.....Shear angle from x in degree, shear direction is on xy-plane
  real(8):: zshear_angle = 0d0

!.....PKA for radiation damage
  integer:: iatom_pka = -1
  real(8):: pka_energy = -1.d0 ! in eV
  real(8):: pka_theta = 0.d0  ! in degree
  real(8):: pka_phi = 0.d0    ! in degree
  
!.....Structure analysis: CNA, a-CNA
  character(len=128):: cstruct = 'none'
  integer:: istruct = 1
  real(8):: rc_struct = 2.5d0

contains
!=======================================================================
  subroutine check_cmin()

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
  subroutine get_lspr(nao,nno,lspro)
!
!   Accessor method for lspr.
!
    integer,intent(out):: nao,nno
    integer,allocatable,intent(out):: lspro(:,:)

    nao = namax
    nno = nnmax
    if( allocated(lspro) ) then
      if( size(lspro).ne.(nnmax+1)*namax) then
        deallocate(lspro)
        allocate(lspro(0:nno,nao))
      endif
    else
      allocate(lspro(0:nno,nao))
    endif
    
  end subroutine get_lspr
!=======================================================================
end module pmdvars
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd lib"
!     End:
