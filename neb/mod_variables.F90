  module variables
    implicit none
    save
!.....Constants
    integer,parameter:: nismax = 10
    real(8),parameter:: aut   = 2.41889d-17
!.....Energy
    real(8),parameter:: ehrt  = 27.2116d0
    real(8),parameter:: ehrt2j= 4.359744d-18
    real(8),parameter:: j2ehrt= 1d0 /ehrt2j
    real(8),parameter:: ev2hrt= 1d0 /ehrt
!.....Bohr radius
    real(8),parameter:: aume  = 9.1093897d-31
    real(8),parameter:: aump  = 1.67262171d-27
    real(8),parameter:: amu   = aump/aume
    real(8),parameter:: bohr  = 0.5291772d-10
    real(8),parameter:: aa2bohr= 1d-10 /bohr
!.....Boltzmann factor in atomic unit
    real(8),parameter:: fkb   = 1.3806503d-23 *aut**2 /aume /bohr**2

!.....Variables to be given via "in.neb"................................
    integer:: nslc
    integer:: nstp = 100
    integer:: iclmb = 0
    character(len=6):: ciofmt = 'ascii '
    character(len=20):: cminimize = 'velocity_damping    '
    character(len=6):: cfprg = 'pmd   '
    character(len=6):: cmethod = 'nudged'
    real(8):: dmp     = 0.8d0
    logical:: lcsdmp  = .true.
    real(8):: tinit   = -100d0
    logical:: ltctl   = .false.
    real(8):: treq    = 300d0
    real(8):: trlx    = 1d-13
    logical:: lcnvg   = .true.
    real(8):: feps    = 1d-3
    real(8):: deps    = 1d-3
    real(8):: scnst   = 1d-2
    real(8):: dt      = 1d-15
    real(8):: fmv(1:3,0:9)= &
         (/ 0.0d0, 0.0d0, 0.0d0, &
            1.0d0, 1.0d0, 1.0d0, &
            1.0d0, 1.0d0, 1.0d0, &
            1.0d0, 1.0d0, 1.0d0, &
            1.0d0, 1.0d0, 1.0d0, &
            1.0d0, 1.0d0, 1.0d0, &
            1.0d0, 1.0d0, 1.0d0, &
            1.0d0, 1.0d0, 1.0d0, &
            1.0d0, 1.0d0, 1.0d0, &
            1.0d0, 1.0d0, 1.0d0  &
            /)

!.....Variables used in each slice......................................
!  natm: num of atoms
!  nis: num of species
    integer:: natm,nis
    real(8):: h(3,3,0:1),hunit
    real(8):: am(1:nismax)= 1d0
    real(8),allocatable:: ra(:,:,:),va(:,:,:),fa(:,:,:),tag(:,:)

  end module variables
