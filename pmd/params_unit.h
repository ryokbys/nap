!
!  Note:
!    unit of length = 1 Angstrom = 1e-10 meter
!    unit of time = 1e-15 second
!    unit of energy = 1 eV = 1.60217657e-19 Joule
!    unit of mass = unified atomic mass unit
!    unit of force = 1 eV/1Angstrom
!.....Pi
      real(8),parameter:: pi = 3.14159265358979d0
!-----time
      real(8),parameter:: ut   = 1d-15
      real(8),parameter:: fs2s = ut
      real(8),parameter:: s2fs = 1d0/ut
      real(8),parameter:: atu  = 2.418884326505d-17 ! atomic time unit
      real(8),parameter:: s2atu = 1d0/atu
      real(8),parameter:: fs2atu = ut/atu
      real(8),parameter:: atu2fs = 1d0/fs2atu
!-----energy
      real(8),parameter:: ev2j = 1.60217657d-19
      real(8),parameter:: j2ev = 1d0/ev2j
      real(8),parameter:: ht2ev = 27.2114d0
      real(8),parameter:: ev2ht = 1d0/ht2ev
!-----length
      real(8),parameter:: ang  = 1.0d-10
      real(8),parameter:: ang2m= ang
      real(8),parameter:: m2ang= 1d0/ang
      real(8),parameter:: bohr = 0.52917721092d-10
      real(8),parameter:: ang2bohr = ang/bohr
      real(8),parameter:: bohr2ang = 1d0/ang2bohr
!.....pressure
      real(8),parameter:: uprs   = ev2j/ang**3
      real(8),parameter:: up2gpa = 1d-9 *uprs
      real(8),parameter:: gpa2up = 1d0 /up2gpa
!-----mass of electron
      real(8),parameter:: ume  = 9.1093897d-31
!-----mass of 1/12 of Carbon
      real(8),parameter:: ump  = 1.660538921d-27
      real(8),parameter:: ume2ump = ume/ump
      real(8),parameter:: ump2ume = ump/ume
      real(8),parameter:: ump2kg = ump
      real(8),parameter:: kg2ump = 1d0/ump
!-----unified atomic mass unit
      real(8),parameter:: umass= ump
      real(8),parameter:: amu2kg= umass
      real(8),parameter:: kg2amu= 1d0/umass
!-----Boltzmann factor in eV/K
      real(8),parameter:: fkb  = 8.6173324d-5
      real(8),parameter:: k2ev = fkb
      real(8),parameter:: ev2k = 1d0 /k2ev
      real(8),parameter:: k2j  = k2ev *ev2j
      real(8),parameter:: j2k  = j2ev *ev2k
!-----special energy unit = [ump*Ang**2/fs**2]
      real(8),parameter:: j2ue  = kg2ump*m2ang**2/s2fs**2
      real(8),parameter:: ue2j  = 1d0/j2ue
      real(8),parameter:: ev2ue = ev2j *j2ue
      real(8),parameter:: ue2ev = 1d0/ev2ue
      real(8),parameter:: k2ue  = k2j*j2ue
      real(8),parameter:: ue2k  = ue2j*j2k
!.....Plank constant in eV*fs (6.62607e-34 /1602e-19 /1e-15)
      real(8),parameter:: plankh = 4.135667733d0
