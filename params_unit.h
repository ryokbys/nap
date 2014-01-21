!
!  Note:
!    unit of length = 1 Angstrom = 1e-10 meter
!    unit of time = 1e-15 second
!    unit of energy = 1 eV = 1.60217657e-19 Joule
!    unit of mass = unified atomic mass unit
!    unit of force = 1 eV/1Angstrom
!
!-----time
      real(8),parameter:: ut   = 1d-15
!-----energy
      real(8),parameter:: ev2j = 1.60217657d-19
      real(8),parameter:: j2ev = 1d0/ev2j
!-----Angstrom
      real(8),parameter:: ang  = 1.0d-10
!.....pressure
      real(8),parameter:: uprs   = ev2j/ang**3
      real(8),parameter:: up2gpa = 1d-9 *uprs
      real(8),parameter:: gpa2up = 1d0 /up2gpa
!-----mass of electron
      real(8),parameter:: ume  = 9.1093897d-31
!-----mass of 1/12 of Carbon
      real(8),parameter:: ump  = 1.660538921d-27
!-----unified atomic mass unit
      real(8),parameter:: umass= ump
!-----Boltzmann factor in eV/K
      real(8),parameter:: fkb  = 8.6173324d-5
