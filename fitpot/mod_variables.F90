module variables
  implicit none 
  save
  integer:: nsmpl
  integer:: niter= 1
  integer:: niter_eval= 1
  character(len=128):: cfmethod= 'BFGS'
  character(len=128):: cmaindir= 'data_set'
  character(len=128):: cparfile= 'in.vars.fitpot'
  character(len=128):: crunmode= 'serial'
  character(len=128):: cevaltype= 'absolute' ! (absolute|relative)
  character(len=128):: csmplist= ''
  integer:: nprcs= 1
  real(8):: epse= 1d-4
  real(8):: epsf= 1d-4
  real(8):: xtol= 1d-4
  real(8):: gtol= 1d-5
  real(8):: ftol= 1d-5
  integer,parameter:: maxnsp= 4
  real(8):: eatom(maxnsp)
  logical:: lematch= .true.   ! energy matching
  logical:: lfmatch= .false.  ! force matching
  logical:: lsmatch= .false.  ! stress matching
  integer:: nsmpl_outfrc = 20000
  character(len=128):: cnormalize= 'variance'
  logical:: lnormalized = .false.
  character(len=128):: cpot= 'NN'
  logical:: lgrad  = .true.
  logical:: lgscale= .false.
  real(8):: gscl   = 1d0
  integer:: nfpsmpl= -10
  real(8):: seqcoef= 1d-2
  integer:: iprint = 1
  character(len=128):: csgdupdate= 'adadelta'
  integer:: nsgdbsize = 1   ! batch size per process for SGD
  integer,allocatable:: ismplsgd(:)
  real(8):: r0sgd = 1.0
!.....training or test
  real(8):: ratio_test= 0.1d0
  logical:: test_assigned = .false.  ! already assigned test set?
!.....initializing parameters
  character(len=128):: cinitv= 'read'
  real(8):: vinitsgm = 1d0
  real(8):: vinitmu  = 0d0
  real(8):: vinitrs  = 12345d0
!.....Atomic forces over this value will not be used for fitting
  real(8):: force_limit = 100d0

!.....Max limit of number of species
  integer,parameter:: mspcs = 9
!.....Max value of species-ID
  integer:: maxisp
  real(8):: ebase(mspcs)
  real(8):: swgt2trn,swgt2tst

!.....max num of atoms among reference data
  integer:: maxna = 0

  integer:: nwgtindiv = 0
  character(len=128),allocatable:: cdirlist(:),cwgtindiv(:)
  integer,allocatable:: nalist(:),iclist(:)
  real(8),allocatable:: wgtindiv(:)

!.....sample error
  integer:: nserr = 0
  character(len=128),allocatable:: cserr(:)
  real(8),allocatable:: seerr(:),sferr(:),sserr(:)
  
!.....sample weights
  integer:: nswgt = 0
  character(len=128),allocatable:: cswgt(:)
  real(8),allocatable:: swerg0(:),swdenom(:)
  
  type mdsys
    character(len=128):: cdirname
    integer:: natm,nfcal
    real(8):: h0,h(3,3),epot,eref,wgt,esub
    real(8):: eerr = 1.0d-3  ! in eV
    real(8):: ferr = 0.1d0   ! in eV/A
    real(8):: serr = 0.1d0   ! in GPa
    real(8),allocatable:: tag(:)
    real(8),allocatable:: ra(:,:),fa(:,:),fref(:,:)
    real(8):: strs(3,3),sref(3,3),ssub(3,3)
    integer,allocatable:: ifcal(:)
    real(8),allocatable:: fabs(:)
    real(8),allocatable:: va(:,:),strsi(:,:,:),eki(:,:,:),epi(:)&
         ,chg(:),chi(:),fsub(:,:),eatm(:)
    real(8),allocatable:: gwe(:),gwf(:,:,:),gws(:,:)
    character(len=2),allocatable:: symbols(:)
    integer:: ispmax       ! Max isp in the sample
    integer:: naps(mspcs)  ! num of atoms per species
    integer:: iclass       ! 1: training,  2: test
    logical:: charge_set = .false.
!.....Related to descriptors
    integer:: nsf,nal,nnl
    real(8),allocatable:: gsf(:,:),gsfo(:,:) &
         ,dgsf(:,:,:,:),igsf(:,:,:)
!.....Specific to NN
    real(8),allocatable:: hl1(:,:)
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
  real(8):: rcut,rc3
!.....Cutoff for other FFs that are subtracted [default: 5.0 Ang]
  real(8):: rc_other = 5.0

  real(8):: time0,tcomm,tfunc,tgrad
  integer:: nfunc,ngrad
  integer:: iflag

  real(8),allocatable:: erefl(:),erefg(:),epotl(:),epotg(:),esubl(:),esubg(:)
  real(8),allocatable:: frefl(:,:,:),frefg(:,:,:),fal(:,:,:)&
       ,fag(:,:,:),fsubl(:,:,:),fsubg(:,:,:)
  real(8),allocatable:: srefl(:,:,:),srefg(:,:,:),strsl(:,:,:), &
       strsg(:,:,:),ssubl(:,:,:),ssubg(:,:,:)
  real(8),allocatable:: swgtl(:),swgtg(:)
  real(8),allocatable:: eerrl(:),eerrg(:)
  real(8),allocatable:: ferrl(:),ferrg(:)
  real(8),allocatable:: serrl(:),serrg(:)

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
  integer,parameter:: nn_nlmax = 2
  integer:: nn_nl = 0
  integer:: nn_nhl(0:nn_nlmax)

end module variables
