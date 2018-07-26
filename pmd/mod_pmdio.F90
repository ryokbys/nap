module pmdio
!-----------------------------------------------------------------------
!                     Last modified: <2018-07-26 18:00:52 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
  implicit none
  save

  character(len=20),parameter:: cinpmd='in.pmd'

  integer:: ntot0,ntot
!.....data of total system
  real(8),allocatable:: rtot(:,:),vtot(:,:),stot(:,:,:),epitot(:) &
       ,ekitot(:,:,:),tagtot(:),atot(:,:),chgtot(:),chitot(:)
  real(8):: hunit,h(3,3,0:1)
  
  integer:: nstp = 0
  integer:: minstp = 0
  integer:: nerg = 1000
  integer:: ifpmd= 1
  integer:: npmd = 10
  integer:: ifsort= 1
  real(8):: dt = 1d0
  real(8):: vardt_len = 0.1d0  ! Length criterion for variable time-step
  real(8):: rc = 5.0d0
  real(8):: rbuf= 0d0
  integer:: ifdmp= 0 ! 0:none, 1:damped-MD, 2:FIRE
  character(len=20):: cmin= ''
  real(8):: dmp  = 0.9d0
  real(8):: eps_conv = 1d-4
  integer:: n_conv = 1
!.....temperature
  character(len=20):: ctctl='none'
  integer:: iftctl= 0
  real(8):: tinit= -1d0
  real(8):: tfin = -1d0
  real(8):: ttgt(9)
  data ttgt / 300d0, 300d0, 300d0, 300d0, 300d0, 300d0, &
       300d0, 300d0, 300d0 /
  real(8):: trlx = 100d0
!.....Remove translational motion every NRMTRANS steps,
!.....if it is less than 1, not to remove translational motion.
  integer:: nrmtrans = 1
!.....Coulomb system?
  integer:: ifcoulomb = 0
!.....temperature distribution on x-direction
  logical:: ltdst= .false.
  integer:: ntdst= 1
!.....shear stress
  real(8):: shrst = 0.0d0
!.....factors on each moving direction
  real(8):: fmv(3,0:9)
  data fmv &
       / 0d0, 0d0, 0d0, &
         1d0, 1d0, 1d0, &
         1d0, 1d0, 1d0, &
         1d0, 1d0, 1d0, &
         1d0, 1d0, 1d0, &
         1d0, 1d0, 1d0, &
         1d0, 1d0, 1d0, &
         1d0, 1d0, 1d0, &
         1d0, 1d0, 1d0, &
         1d0, 1d0, 1d0  &
         /
!.....whether compute stress or not
  logical:: lstrs0 = .true.
!.....barostat
  character(len=20):: cpctl='none'
  real(8):: ptgt   = 0d0
  real(8):: srlx   = 100d0
  real(8):: stbeta = 1d0
  real(8):: strfin = 0.0d0
  real(8):: stnsr(3,3)
  real(8):: stgt(1:3,1:3)= 0d0
  logical:: lcellfix(1:3,1:3)= .false.
!.....charge optimize or variable charge
  logical:: lvc = .false.
!.....Charge setting
  character(len=20):: chgfix='input'
!.....print level
!  0:quiet, 1:normal,
!  >10:fitpot data
  integer:: iprint= 1

  character(len=128):: cpmdini = 'pmdini'
  character(len=128):: cpmdfin = 'pmdfin'
  character(len=6):: ciofmt='ascii '
  character(len=20):: cforce='none'
!!$  integer:: numff = 0 ! number of force-fields
!!$  character(len=20),allocatable:: cffs(:)  ! force-fields
!.....max. num. of species
  integer,parameter:: nspmax= 9
!.....mass
  real(8):: am(1:nspmax)= 12.0d0
!.....charges
  real(8):: schg(1:nspmax)= 0d0
!.....zload type
  character(len=5):: czload_type= 'no'

!.....Boundary condition: p = periodic, f = free, w = wall
  character(len=3):: boundary = 'ppp'

!.....PKA for radiation damage
  real(8):: pka_energy = -1.d0 ! in eV
  real(8):: pka_theta = 0.d0  ! in degree
  real(8):: pka_phi = 0.d0    ! in degree
end module pmdio
