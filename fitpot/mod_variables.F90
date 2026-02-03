module variables
  use pmdvars,only: nspmax
  use descriptor,only: desc
  implicit none
  save
  integer:: nsmpl
  integer:: niter= 1
  integer:: niter_eval= 1
  character(len=128):: cfmethod= 'BFGS'
  character(len=128):: cparfile= 'in.vars.fitpot'
  character(len=128):: crunmode= 'serial'
  character(len=128):: cevaltype= 'absolute' ! (absolute|relative)
!.....Sample file
  character(len=128):: csmplfile= ''  ! used since 250412
  character(len=128):: cdatasetdir= 'dataset'  ! obsolete since 250412
  character(len=128):: csmplistfile= ''
  character(len=128):: csmplftype= 'dir'  ! (dir|file)
  integer:: nprcs= 1
  real(8):: epse= 1d-4
  real(8):: epsf= 1d-4
  real(8):: xtol= 1d-4
  real(8):: gtol= 1d-5
  real(8):: ftol= 1d-6
!.....This specorder is the one for whole fitpot process not for each sample.
  character(len=3),dimension(nspmax):: specorder = &
       (/'x','x','x','x','x','x','x','x','x'/)
  integer:: nsp = 0  ! number of species used
  real(8):: eatom(nspmax)
  logical:: lematch= .false.   ! energy matching
  logical:: lfmatch= .false.  ! force matching
  logical:: lsmatch= .false.  ! stress matching
!.....Forces of species listed in it are to be neglected
  character(len=3):: cspcs_neglect(nspmax)
  integer:: nspcs_neglect = 0
  integer:: nsmpl_outfrc = 20000
  character(len=128):: cnormalize= 'std'
  logical:: lnormalize  = .false.  ! whether of not to normalize vars
  logical:: lnormalized = .false.  ! whether or not vars already normalized
  character(len=128):: cpot = ''
  character(len=128):: cpotlow = ''
!.....Cutoff for FF that is to be optimized [default: 5.0 Ang]
  real(8):: rcut   = 6.0d0
  real(8):: rc3    = 3.0d0
!.....Cutoff for other FFs that are subtracted [default: 5.0 Ang]
  real(8):: rc_other = 5.0
!!$  logical:: lgrad  = .true.
!!$  logical:: lgscale= .false.
  real(8):: gscl   = 1d0
  integer:: nfpsmpl= -10
  real(8):: seqcoef= 1d-2
  integer:: iprint = 1
!!$  character(len=128):: csgdupdate= 'adadelta'
!!$  integer:: nsgdbsize = 1   ! batch size per process for SGD
!!$  integer,allocatable:: ismplsgd(:)
!!$  real(8):: r0sgd = 1.0
!.....training or test
  real(8):: ratio_test= 0.1d0
  logical:: test_assigned = .false.  ! already assigned test set?
!.....initializing parameters: read or gauss/gaussian
  character(len=128):: cinitv= 'read'
  real(8):: vinitsgm = 1d0
  real(8):: vinitmu  = 0d0
  real(8):: vinitrs  = 12345d0
!.....Characteristic force scale used for limit or decay length
  real(8):: f_scale = 1d0
!.....Characteristic stress scale used for limit or decay length
  real(8):: s_scale = 1d0
!.....Loss function type: LS (least-square), Huber
  character(len=128):: ctype_loss = 'LS'
!.....Scaling type of force|stress in loss function:
!.....  one (default) | abs | rel | exp
  character(len=128):: cfrc_scale = 'one' !
  character(len=128):: cstrs_scale = 'one' !
!.....Gaussian density weight
  logical:: lgdw  = .false.  ! flag for GDW
  logical:: lgdwed = .false.  ! whether compuation of GDW is finished
  real(8):: gdsgm = 1.0d0  ! sigma in GDF, assuming that GDF is normalized

!.....Max value of species-ID
  integer:: maxisp
  real(8):: ebase(nspmax)
  real(8):: swgt2trn,swgt2tst
  logical:: interact(nspmax,nspmax)
  logical:: interact3(nspmax,nspmax,nspmax)

!.....Max num of atoms among reference samples
  integer:: maxna = 0
!.....Max num of atoms among nodes
  integer:: maxnin = 0
!.....Total num of atoms among reference samples
  integer:: natot = 0
!.....Max num of atoms whose forces are used for force-matching
  integer:: maxnf = 0

  integer:: nwgtindiv = 0
  character(len=128),allocatable:: cdirlist(:),cwgtindiv(:),csmplist(:)
  integer,allocatable:: nalist(:),iclist(:)
  real(8),allocatable:: wgtindiv(:)

!.....Relative weights bewteen energy, force, and stress
  real(8):: wgte = 0.5d0
  real(8):: wgtf = 0.3d0
  real(8):: wgts = 0.2d0

!.....sample error
  integer:: nserr = 0
  character(len=128),allocatable:: cserr(:)
  real(8),allocatable:: seerr(:),sferr(:),sserr(:)

!.....sample weights
  integer:: nswgt = 0
  character(len=128),allocatable:: cswgt(:)
  real(8),allocatable:: swgt0(:)

!.....Limit number of forces to be evaluated
  real(8):: rate_eval_frc = 1.0d0

  type mdsys
    character(len=128):: cdirname
    character(len=128):: csmplname
    integer:: natm,nfcal
    real(8):: h0,h(3,3,0:1),epot,eref,wgt,esub
    real(8):: eerr = 1.0d-3  ! in eV
    real(8):: ferr = 0.1d0   ! in eV/A
    real(8):: serr = 0.1d0   ! in GPa
    real(8),allocatable:: tag(:)
    real(8),allocatable:: ra(:,:),fa(:,:),fref(:,:)
    real(8):: strs(3,3),sref(3,3),ssub(3,3)
    integer,allocatable:: ifcal(:)
    real(8),allocatable:: fabs(:)
    real(8),allocatable:: va(:,:),strsi(:,:,:),eki(:,:,:),epi(:)&
         ,fsub(:,:),eatm(:)
    real(8),allocatable:: aux(:,:)
!.....Gradients w.r.t. parameters to avoid recalculation of gradw
!     in the case of linear regression (uf3 as well).
!     To save the memory usage, use single-precision variable
    real(4),allocatable:: gwe(:), gwf(:,:,:),gws(:,:)
!.....Check flag for {e,f,s}ref data
    logical:: leref_given = .false.
    logical:: lfref_given = .false.
    logical:: lsref_given = .false.
!.....This specorder is for this sample
    character(len=3),dimension(nspmax):: specorder  &
         = (/'x','x','x','x','x','x','x','x','x'/)
    integer:: ispmax       ! Max isp in the sample
    integer:: naps(nspmax)  ! num of atoms per species
    integer:: iclass       ! 1: training,  2: test
    logical:: charge_set = .false.
!.....Limit number of forces to be evaluated
    logical,allocatable:: lfrc_eval(:)
!.....Related to descriptors
    integer:: nsf,nal,nnl
    real(8),allocatable:: gsf(:,:),gsfo(:,:),dgsf(:,:,:,:)
    integer(2),allocatable:: igsf(:,:,:)
!.....Gaussian density functions (GDF) for atoms and weights using the GDF
    real(8),allocatable:: gdf(:), gdw(:)
!.....Specific to NN
    real(8),allocatable:: hl1(:,:)
!.....Design-matrix for force-matching
    real(8),allocatable:: dgsfa(:,:,:)
!.....flag for skipping gradw, will be set .true. when linear regression-type forces
    logical:: lgrad_done = .false.
  end type mdsys
  real(8):: erefmin
  real(8):: gsfmean,gsfvar
  real(8),allocatable:: gsfms(:),gsfvs(:),sgms(:),sgmis(:)
  real(8),allocatable:: gsfcorr(:,:),gsfss(:)
  real(8),parameter:: sgm_min = 1d-2
  real(8),parameter:: sq_min = 1d-2

  integer:: nsmpl_trn,nsmpl_tst
  type(mdsys),allocatable:: samples(:)
  integer:: nvars
  real(8),allocatable:: vars(:),vranges(:,:),gvar(:),dvar(:)
  integer:: ibest
  real(8):: fbest
  real(8),allocatable:: vbest(:)

  integer:: netrn,netst,nftrn,nftst,nstrn,nstst
  real(8):: evtrn,evtst,fvtrn,fvtst,svtrn,svtst
  real(8):: etrndnm,etstdnm,ftrndnm,ftstdnm,strndnm,ststdnm

  real(8):: time0,tcomm,tfunc,tgrad,twait
  real(8):: terg, tfrc, tstrs
  integer:: nfunc,ngrad
  integer:: iflag

  logical,allocatable:: lfevall(:,:),lfevalg(:,:)
  real(8),allocatable:: erefl(:),erefg(:),epotl(:),epotg(:),esubl(:),esubg(:)
  real(8),allocatable:: frefl(:,:,:),frefg(:,:,:),fal(:,:,:)&
       ,fag(:,:,:),fsubl(:,:,:),fsubg(:,:,:)
  real(8),allocatable:: srefl(:,:,:),srefg(:,:,:),strsl(:,:,:), &
       strsg(:,:,:),ssubl(:,:,:),ssubg(:,:,:)
  real(8),allocatable:: swgtl(:),swgtg(:)
  real(8),allocatable:: eerrl(:),eerrg(:)
  real(8),allocatable:: ferrl(:),ferrg(:)
  real(8),allocatable:: serrl(:),serrg(:)
  logical,allocatable:: lfcall(:,:),lfcalg(:,:)
  real(8),allocatable:: gdwl(:,:), gdwg(:,:)

!.....Force-fields which are subtracted from reference values
!.....and whose parameters to be fitted
  integer:: nsubff,nff
  character(len=20),allocatable:: csubffs(:),cffs(:)

!.....For the purpose of setting the reference energy from a structure
!.....within the whole sample structures.
  character(len=128):: crefstrct= ''
  integer:: myidrefsub = -1
  integer:: isidrefsub = -1
  real(8):: erefsub = 0d0

!.....Constants for NN2 potential
!!$  integer,parameter:: nn_nlmax = 2
  integer:: nn_nl = 0
  integer,allocatable:: nn_nhl(:)
!!$  integer:: nn_nhl(0:nn_nlmax+1)
  integer:: nn_sigtype = 2
  real(8):: nn_asig = 0.01d0

  integer:: mem = 0
  real(8):: dmem = 0d0

!.....For descriptors except Chebyshev
  type(desc),allocatable:: descs(:)
  integer:: nsp_desc,nsf_desc,nsf2_desc,nsf3_desc,nsff_desc
!.....List of isf's for each pair (ilsf2) and angle (ilsf3)
  integer,allocatable:: ilsf2(:,:,:),ilsf3(:,:,:,:)
  logical:: lcheby = .false.
  real(8),allocatable:: wgtsp_desc(:)
!.....Function types and num of constatns for types
  integer,parameter:: max_ncnst = 2
  integer:: ncnst_type(200)
  integer:: ncomb_type(200)
  real(8):: cnst(max_ncnst)

!.....ZBL potential-related parameters
  logical:: zbl_interact(nspmax,nspmax)
  real(8):: zbl_qnucl(nspmax)
  real(8):: zbl_rc
  real(8):: zbl_ri(1:nspmax) = -1d0
  real(8):: zbl_ro(1:nspmax) = -1d0

!.....Penalty is moved from minimize module
  character(len=128):: cpenalty = 'none'
  real(8):: penalty = 1d-15
!.....Repulsion correction for short distances in UF3 potential
  logical:: l_correct_short = .false.
  real(8):: short_radii(nspmax,nspmax)
  real(8):: pwgt2b = 0d0
  real(8):: pwgt2bd= 0d0
  real(8):: pwgt2bs= 0d0
  real(8):: pwgt3b = 0d0
  real(8):: pwgt3bd= 0d0

  real(8):: pwgt_repul = 1d-15  ! penalty for repulsion
  real(8):: valence_chgs(nspmax), core_chgs(nspmax)
  real(8):: repul_radii(nspmax,nspmax)
  integer:: n_repul_pnts = 0  ! # of points incl. 0 (orig.) and end
  real(8),allocatable:: drepul_tbl(:,:,:)  ! Repulsion gradient table

!------- MINIMIZATION PARAMETERS ---------------------------
!.....number of convergence criteria achieved
  integer:: numtol = 1

!.....penalty: lasso or ridge or smooth
  character(len=128):: cpena= 'none'
  character(len=128):: clinmin= 'backtrack'
  character(len=128):: cfsmode= 'grad'  ! [grad,grad0corr,df0corr]
  real(8):: pwgt = 1d-15

!.....SGD parameters
  integer:: nsgdbsize = 1
  integer:: nsgdbsnode = 1
  integer,allocatable:: ismask(:)
  character(len=128):: csgdupdate = 'normal'
  real(8):: sgd_rate_ini = 0.001d0
  real(8):: sgd_rate_fin = -0.001d0
  real(8):: sgd_eps = 1.0d-8
!.....Parameters for ADAM and AdaBound
  real(8):: adam_b1 = 0.9d0
  real(8):: adam_b2 = 0.999d0

!.....Group FS inner loop
  integer:: ninnergfs=100
  character(len=128):: cread_fsmask = ''
  character(len=128):: cfs_xrefresh = 'random' ! [zero, random, none]
  integer:: maxfsrefresh = 2

!.....Max iteration for line minimization
  integer:: niter_linmin   = 15
!.....Decreasing factor, should be < 1.0
  real(8):: fac_dec        = 0.2d0
!.....Increasing factor, should be > 1.0
  real(8):: fac_inc        = 5.0d0
!.....Armijo parameters
  real(8):: armijo_xi      = 1.0d-4
  real(8):: armijo_tau     = 0.5d0
  integer:: armijo_maxiter = 15

!.....CG
  integer:: icgbtype = 1 ! 1:FR, 2:PRP, 3:HS, 4:DY

!.....L-BFGS
  integer:: m_lbfgs   = 10

  real(8):: fupper_lim = 1d+5

!-----DEBUGGING------------------------------------------------------
  integer:: id_check_grad = 0

contains
  subroutine init_variables()

    interact(:,:) = .true.
    cspcs_neglect(:) = 'x'
    short_radii(:,:) = 0d0
    nsubff = 0

  end subroutine init_variables
!=======================================================================
  function csp_in_neglect(csp)
!
!  Return whether or not the given csp is in the neglect_list.
!
    character(len=*),intent(in):: csp

    logical:: csp_in_neglect
    integer:: i

    csp_in_neglect = .false.
    do i=1,nspmax
      if( trim(csp).eq.trim(cspcs_neglect(i)) ) then
        csp_in_neglect = .true.
        return
      endif
    enddo
    return
  end function csp_in_neglect
!=======================================================================
  function num_interact(ndim)
    integer,intent(in):: ndim

    integer:: num_interact
    integer:: i,j,k

    num_interact = 0
    if( ndim.eq.2 ) then  ! pairs
      do i=1,nspmax
        do j=i,nspmax
          if( .not. interact(i,j) ) cycle
          num_interact = num_interact +1
        enddo
      enddo
      return
    else if( ndim.eq.3 ) then  ! triplets
      do i=1,nspmax
        do j=1,nspmax
          do k=j,nspmax
            if( .not. interact3(i,j,k) ) cycle
            num_interact = num_interact +1
          enddo
        enddo
      enddo
      return
    else
      print *,' NDIM must be either 2 or 3.'
      stop 1
    endif
    return
  end function num_interact
!=======================================================================
  function csp2isp(csp) result(isp)
!
!  Convert csp to isp.
!  If not found, return -1.
!
    character(len=*),intent(in):: csp
    integer:: isp

    integer:: i

    isp = -1
    do i=1,nspmax
      if( trim(csp).eq.trim(specorder(i)) ) then
        isp = i
        return
      endif
    enddo
    return
  end function csp2isp
!=======================================================================
end module variables
