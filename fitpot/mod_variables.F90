module variables

  save
  integer:: nsmpl
  integer:: niter= 1
  integer:: niter_eval= 1
  character(len=128):: cfmethod= 'BFGS'
  character(len=128):: cmaindir= 'data_set'
  character(len=128):: cparfile= 'in.params.NN'
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
  logical:: lfmatch= .false.
  character(len=128):: cnormalize= 'variance'
  character(len=128):: cpot= 'NN'
  logical:: lgrad  = .true.
  logical:: lgscale= .false.
  real(8):: gscl   = 0.1d0
  real(8):: fred   = -1d0
  integer:: nfpsmpl= -10
  logical:: lswgt  = .false.
  real(8):: swerg = 1d0
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

!.....max num of species
  integer,parameter:: mspcs = 4
  real(8):: ebase(mspcs)
  real(8):: swgt2trn,swgt2tst

  integer:: nwgtindiv = 0
  character(len=128),allocatable:: cdirlist(:),cwgtindiv(:)
  integer,allocatable:: nalist(:),iclist(:)
  real(8),allocatable:: wgtindiv(:)

!.....sample error
  integer:: nserr = 0
  character(len=128),allocatable:: cserr(:)
  real(8),allocatable:: seerr(:),sferr(:)
  
  type mdsys
    character(len=128):: cdirname
    integer:: natm,nfcal
    real(8):: h0,h(3,3),epot,eref,wgt,esub
    real(8):: eerr = 1.0d-3  ! in eV
    real(8):: ferr = 0.1d0   ! in eV/A
    real(8),allocatable:: tag(:)
    real(8),allocatable:: ra(:,:),fa(:,:),fref(:,:)
    integer,allocatable:: ifcal(:)
    real(8),allocatable:: fabs(:)
    real(8),allocatable:: va(:,:),strsi(:,:,:),eki(:,:,:),epi(:)&
         ,chg(:),chi(:),fsub(:,:)
    integer:: naps(mspcs)  ! num of atoms per species
    integer:: iclass       ! 1: training,  2: test
  end type mdsys
  real(8):: erefmin

  integer:: nsmpl_trn,nsmpl_tst
  type(mdsys),allocatable:: samples(:)
  integer:: nvars
  real(8),allocatable:: vars(:),vranges(:,:),gvar(:),dvar(:)
  real(8):: rcut,rc3

  real(8):: time0,tcomm,tfunc,tgrad
  integer:: nfunc,ngrad
  integer:: iflag

  real(8),allocatable:: erefl(:),erefg(:),epotl(:),epotg(:)
  real(8),allocatable:: frefl(:,:,:),frefg(:,:,:),fal(:,:,:)&
       ,fag(:,:,:)
  real(8),allocatable:: swgtl(:),swgtg(:)
  real(8),allocatable:: eerrl(:),eerrg(:)
  real(8),allocatable:: ferrl(:),ferrg(:)

!.....Force-fields which are subtracted from reference values
  integer:: numff
  character(len=20),allocatable:: cffs(:)

end module variables
