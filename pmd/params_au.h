!
!  Note:
!    atomic unit of length = 0.5292Angstrom = 5.292e-11 meter
!    atomic unit of time = 2.418e-17 second
!    atomic unit of energy = 27.2114eV = 4.36e-18 Joule
!    atomic unit of mass = 9.109e-31 kg
!    atomic unit of force = 27.21eV/0.5292Angstrom
!
!-----time
      real(rp),parameter:: aut   = 2.41889e-17_rp
!-----energy
      real(rp),parameter:: ehrt  = 27.2114_rp
      real(rp),parameter:: ehrt2j= 4.359744e-18_rp
      real(rp),parameter:: j2ehrt= 1.0_rp /ehrt2j
      real(rp),parameter:: ev2hrt= 1.0_rp /ehrt
!-----Bohr radius
      real(rp),parameter:: bohr  = 0.5291772e-10_rp
      real(rp),parameter:: aa2bohr= 1e-10_rp /bohr
      real(rp),parameter:: bohr2aa= 1.0_rp /aa2bohr
!.....pressure
      real(rp),parameter:: aup= ehrt2j/bohr**3
      real(rp),parameter:: aup2gpa= 1e-9_rp*aup
      real(rp),parameter:: gpa2aup= 1.0_rp/aup2gpa
!-----mass of electron
      real(rp),parameter:: aume  = 9.1093897e-31_rp
!-----mass of proton
      real(rp),parameter:: aump  = 1.67262171e-27_rp
!-----atomic mass unit
      real(rp),parameter:: amu   = aump/aume
!-----Boltzmann factor in atomic unit
      real(rp),parameter:: fkb   = 1.3806503e-23_rp *aut**2 /aume /bohr**2
