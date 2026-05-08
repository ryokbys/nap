!
!  Note:
!    atomic unit of length = 0.5292Angstrom = 5.292e-11 meter
!    atomic unit of time = 2.418e-17 second
!    atomic unit of energy = 27.2114eV = 4.36e-18 Joule
!    atomic unit of mass = 9.109e-31 kg
!    atomic unit of force = 27.21eV/0.5292Angstrom
!
!-----time
      real(rp),parameter:: aut   = real(2.41889d-17, rp)
!-----energy
      real(rp),parameter:: ehrt  = real(27.2114d0, rp)
      real(rp),parameter:: ehrt2j= real(4.359744d-18, rp)
      real(rp),parameter:: j2ehrt= real(1d0, rp) /ehrt2j
      real(rp),parameter:: ev2hrt= real(1d0, rp) /ehrt
!-----Bohr radius
      real(rp),parameter:: bohr  = real(0.5291772d-10, rp)
      real(rp),parameter:: aa2bohr= real(1d-10, rp) /bohr
      real(rp),parameter:: bohr2aa= real(1d0, rp) /aa2bohr
!.....pressure
      real(rp),parameter:: aup= ehrt2j/bohr**3
      real(rp),parameter:: aup2gpa= real(1d-9, rp)*aup
      real(rp),parameter:: gpa2aup= real(1d0, rp)/aup2gpa
!-----mass of electron
      real(rp),parameter:: aume  = real(9.1093897d-31, rp)
!-----mass of proton
      real(rp),parameter:: aump  = real(1.67262171d-27, rp)
!-----atomic mass unit
      real(rp),parameter:: amu   = aump/aume
!-----Boltzmann factor in atomic unit
      real(rp),parameter:: fkb   = real(1.3806503d-23, rp) *aut**2 /aume /bohr**2
