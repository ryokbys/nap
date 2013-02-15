c
c  Note:
c    atomic unit of length = 0.5292Angstrom = 5.292e-11 meter
c    atomic unit of time = 2.418e-17 second
c    atomic unit of energy = 27.2114eV = 4.36e-18 Joule
c    atomic unit of mass = 9.109e-31 kg
c    atomic unit of force = 27.21eV/0.5292Angstrom
c
c-----time
      real(8),parameter:: aut   = 2.41889d-17
c-----energy
      real(8),parameter:: ehrt  = 27.2114d0
      real(8),parameter:: ehrt2j= 4.359744d-18
      real(8),parameter:: j2ehrt= 1d0 /ehrt2j
      real(8),parameter:: ev2hrt= 1d0 /ehrt
c-----Bohr radius
      real(8),parameter:: bohr  = 0.5291772d-10
      real(8),parameter:: aa2bohr= 1d-10 /bohr
      real(8),parameter:: bohr2aa= 1d0 /aa2bohr
c.....pressure
      real(8),parameter:: aup= ehrt2j/bohr**3
      real(8),parameter:: aup2gpa= 1d-9*aup
      real(8),parameter:: gpa2aup= 1d0/aup2gpa
c-----mass of electron
      real(8),parameter:: aume  = 9.1093897d-31
c-----mass of proton
      real(8),parameter:: aump  = 1.67262171d-27
c-----atomic mass unit
      real(8),parameter:: amu   = aump/aume
c-----Boltzmann factor in atomic unit
      real(8),parameter:: fkb   = 1.3806503d-23 *aut**2 /aume /bohr**2
