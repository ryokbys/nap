module variables

  integer:: nsmpl
  integer:: niter= 1
  integer:: niter_eval= 1
  character(len=128):: cfmethod= 'BFGS'
  character(len=128):: cmaindir= 'learning_set'
  character(len=128):: cparfile= 'in.params.NN'
  character(len=128):: crunmode= 'serial'
  character(len=128):: cevaltype= 'absolute' ! (absolute|relative)
  integer:: nprcs= 1
  real(8):: epse= 1d-4
  real(8):: epsf= 1d-4
  real(8):: xtol= 1d-5
  real(8):: gtol= 1d-5
  real(8):: ftol= 1d-5
  integer,parameter:: maxnsp= 4
  real(8):: eatom(maxnsp)
  logical:: lfmatch= .true.
  logical:: lreg   = .false.
  character(len=128):: cpot= 'NN'
  logical:: lgrad  = .true.
  logical:: lgscale= .false.
  real(8):: gscl   = 0.1d0
  logical:: lfscale= .false.
  real(8):: fscl   = 1d0
  real(8):: fred   = -1d0
  integer:: nfpsmpl= 10
  logical:: lswgt  = .false.
  real(8):: swerg = 1d0
  real(8):: seqcoef= 1d-2
  integer:: iprint = 1
  real(8),save:: rseed  = 12345d0
!.....training or test
  real(8):: ratio_test= 0.1d0

!.....max num of species
  integer,parameter:: mspcs = 4
  real(8):: ebase(mspcs)
  real(8):: swgt2

  integer:: nwgtindiv = 0
  character(len=128),allocatable,save:: cdirlist(:),cwgtindiv(:)
  integer,allocatable,save:: nalist(:),iclist(:)
  real(8),allocatable,save:: wgtindiv(:)

!.....sample error
  integer:: nserr = 0
  character(len=128),allocatable,save:: cserr(:)
  real(8),allocatable,save:: seerr(:),sferr(:)
  
  type mdsys
    character(len=128):: cdirname
    integer:: natm
    real(8):: h0,h(3,3),epot,eref,wgt
    real(8):: eerr = 1.0d-3 ! in eV
    real(8):: ferr = 0.1d0  ! in eV/A
    real(8),allocatable:: tag(:)
    real(8),allocatable:: ra(:,:),fa(:,:),fref(:,:)
    integer,allocatable:: ifcal(:)
    real(8),allocatable:: fabs(:)
    integer:: naps(mspcs) ! num of atoms per species
    integer:: iclass
  end type mdsys
  real(8):: erefmin

  integer:: nsmpl_trn,nsmpl_tst
  type(mdsys),save,allocatable:: samples(:)
  integer,save:: nvars
  real(8),save,allocatable:: vars(:),vranges(:,:),gvar(:),dvar(:)
  real(8),save:: rcut,rc3

  real(8),save:: time0,tcomm,tfunc,tgrad
  integer,save:: nfunc,ngrad
  integer,save:: iflag

  real(8),allocatable,save:: erefl(:),erefg(:),epotl(:),epotg(:)
  real(8),allocatable,save:: frefl(:,:,:),frefg(:,:,:),fal(:,:,:)&
       ,fag(:,:,:)
  real(8),allocatable,save:: swgtl(:),swgtg(:)
  real(8),allocatable,save:: eerrl(:),eerrg(:)
  real(8),allocatable,save:: ferrl(:),ferrg(:)

end module variables
