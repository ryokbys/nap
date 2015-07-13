module variables

  integer:: nsmpl
  integer:: niter= 1
  integer:: niter_eval= 1
  character(len=128):: cfmethod= 'BFGS'
  character(len=128):: cmaindir= 'learning_set'
  character(len=128):: cparfile= 'in.params.NN'
  character(len=128):: crunmode= 'serial'
  integer:: nprcs= 1
  real(8):: eps= 1d-8
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
  logical:: lswgt  = .false.
  real(8):: swbeta = 1d0
  real(8):: seqcoef= 1d-2
  integer:: iprint = 1
!.....training or test
  real(8):: ratio_test= 0.1d0

  character(len=5),allocatable,save:: cdirlist(:)
  integer,allocatable,save:: nalist(:),iclist(:)
  
  type mdsys
    character(len=5):: cdirname
    integer:: natm
    real(8):: h0,h(3,3),epot,eref
    real(8),allocatable:: tag(:)
    real(8),allocatable:: ra(:,:),fa(:,:),fref(:,:)
    integer:: iclass
  end type mdsys
  real(8):: erefmin

  integer:: nsmpl_trn,nsmpl_tst
  type(mdsys),save,allocatable:: samples(:)
  integer,save:: nvars
  real(8),save,allocatable:: vars(:),vranges(:,:),gvar(:),dvar(:)
  real(8),save:: rcut

  real(8),save:: time0,tcomm,tfunc,tgrad
  integer,save:: nfunc,ngrad
  integer,save:: iflag

  real(8),allocatable,save:: erefl(:),erefg(:),epotl(:),epotg(:)
  real(8),allocatable,save:: frefl(:,:,:),frefg(:,:,:),fal(:,:,:)&
       ,fag(:,:,:)

end module variables
