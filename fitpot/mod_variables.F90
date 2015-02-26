module variables

  integer:: nsmpl
  integer:: nstp= 1
  character(len=128):: cfmethod= 'BFGS'
  character(len=128):: cmaindir= 'learning_set'
  character(len=128):: cparfile= 'in.params.NN'
  character(len=128):: crunmode= 'serial'
  integer:: nprcs= 1
  real(8):: eps= 1d-8
  real(8):: xtol= 1d-5
  integer,parameter:: maxnsp= 3
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
  logical:: lpena  = .false.
  real(8):: pwgt   = 1d0

  character(len=5),allocatable:: cdirlist(:)

  type mdsys
    character(len=5):: cdirname
    integer:: natm
    real(8):: h0,h(3,3),epot,eref
    real(8),allocatable:: tag(:)
    real(8),allocatable:: ra(:,:),fa(:,:),fref(:,:)
  end type mdsys

  type(mdsys),save,allocatable:: samples(:)
  integer:: nvars
  real(8),allocatable:: vars(:),vranges(:,:)
  real(8):: rcut

  real(4),save:: timef,timeg

end module variables
