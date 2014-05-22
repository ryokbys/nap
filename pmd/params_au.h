!
!  Note:
!    atomic unit of length = 0.5292Angstrom = 5.292e-11 meter
!    atomic unit of time = 2.418e-17 second
!    atomic unit of energy = 27.2114eV = 4.36e-18 Joule
!    atomic unit of mass = 9.109e-31 kg
!    atomic unit of force = 27.21eV/0.5292Angstrom
!
!-----time
      real(8),parameter:: aut   = 2.41889d-17
!-----energy
      real(8),parameter:: ehrt  = 27.2114d0
      real(8),parameter:: ehrt2j= 4.359744d-18
      real(8),parameter:: j2ehrt= 1d0 /ehrt2j
      real(8),parameter:: ev2hrt= 1d0 /ehrt
!-----Bohr radius
      real(8),parameter:: bohr  = 0.5291772d-10
      real(8),parameter:: aa2bohr= 1d-10 /bohr
      real(8),parameter:: bohr2aa= 1d0 /aa2bohr
!.....pressure
      real(8),parameter:: aup= ehrt2j/bohr**3
      real(8),parameter:: aup2gpa= 1d-9*aup
      real(8),parameter:: gpa2aup= 1d0/aup2gpa
!-----mass of electron
      real(8),parameter:: aume  = 9.1093897d-31
!-----mass of proton
      real(8),parameter:: aump  = 1.67262171d-27
!-----atomic mass unit
      real(8),parameter:: amu   = aump/aume
!-----Boltzmann factor in atomic unit
      real(8),parameter:: fkb   = 1.3806503d-23 *aut**2 /aume /bohr**2
