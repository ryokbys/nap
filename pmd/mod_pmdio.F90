module pmdio
!-----------------------------------------------------------------------
!                     Last modified: <2017-03-18 21:00:17 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
  implicit none
  save

  integer:: ntot
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
!.....charged system?
  integer:: ifchg = 0
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
  character(len=20):: cpctl='non'
  real(8):: ptgt   = 0d0
  real(8):: srlx   = 100d0
  real(8):: stbeta = 1d-1
  real(8):: strfin = 0.0d0
  real(8):: ptnsr(3,3)
  real(8):: stgt(1:3,1:3)= 0d0

!.....print level
!  0:quiet, 1:normal,
!  >10:fitpot data
  integer:: iprint= 1

  character(len=6):: ciofmt='ascii '
  character(len=20):: cforce='LJ_Ar'
!.....max. num. of species
  integer,parameter:: nismax= 9
!.....mass
  real(8):: am(1:nismax)= 12.0d0
!.....zload type
  character(len=5):: czload_type= 'no'
end module pmdio
