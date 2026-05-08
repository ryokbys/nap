!
!  Note:
!    unit of length = 1 Angstrom = 1e-10 meter
!    unit of time = 1e-15 second
!    unit of energy = 1 eV = 1.60217657e-19 Joule
!    unit of mass = unified atomic mass unit
!    unit of force = 1 eV/1Angstrom
!.....Pi
      real(rp),parameter:: pi = real(3.14159265358979d0, rp)
!-----time
      real(rp),parameter:: ut   = real(1d-15, rp)
      real(rp),parameter:: fs2s = ut
      real(rp),parameter:: s2fs = real(1d0, rp)/ut
      real(rp),parameter:: atu  = real(2.418884326505d-17, rp) ! atomic time unit
      real(rp),parameter:: s2atu = real(1d0, rp)/atu
      real(rp),parameter:: fs2atu = ut/atu
      real(rp),parameter:: atu2fs = real(1d0, rp)/fs2atu
!-----energy
      real(rp),parameter:: ev2j = real(1.60217657d-19, rp)
      real(rp),parameter:: j2ev = real(1d0, rp)/ev2j
      real(rp),parameter:: ht2ev = real(27.2114d0, rp)
      real(rp),parameter:: ev2ht = real(1d0, rp)/ht2ev
!-----length
      real(rp),parameter:: ang  = real(1.0d-10, rp)
      real(rp),parameter:: ang2m= ang
      real(rp),parameter:: m2ang= real(1d0, rp)/ang
      real(rp),parameter:: bohr = real(0.52917721092d-10, rp)
      real(rp),parameter:: ang2bohr = ang/bohr
      real(rp),parameter:: bohr2ang = real(1d0, rp)/ang2bohr
!.....pressure
      real(rp),parameter:: uprs   = ev2j/ang**3
      real(rp),parameter:: up2gpa = real(1d-9, rp) *uprs
      real(rp),parameter:: gpa2up = real(1d0, rp) /up2gpa
!-----mass of electron
      real(rp),parameter:: ume  = real(9.1093897d-31, rp)
!-----mass of 1/12 of Carbon
      real(rp),parameter:: ump  = real(1.660538921d-27, rp)
      real(rp),parameter:: ume2ump = ume/ump
      real(rp),parameter:: ump2ume = ump/ume
      real(rp),parameter:: ump2kg = ump
      real(rp),parameter:: kg2ump = real(1d0, rp)/ump
!-----unified atomic mass unit
      real(rp),parameter:: umass= ump
      real(rp),parameter:: amu2kg= umass
      real(rp),parameter:: kg2amu= real(1d0, rp)/umass
!-----Boltzmann factor in eV/K
      real(rp),parameter:: fkb  = real(8.6173324d-5, rp)
      real(rp),parameter:: k2ev = fkb
      real(rp),parameter:: ev2k = real(1d0, rp) /k2ev
      real(rp),parameter:: k2j  = k2ev *ev2j
      real(rp),parameter:: j2k  = j2ev *ev2k
!-----special energy unit = [ump*Ang**2/fs**2]
      real(rp),parameter:: j2ue  = kg2ump*m2ang**2/s2fs**2
      real(rp),parameter:: ue2j  = real(1d0, rp)/j2ue
      real(rp),parameter:: ev2ue = ev2j *j2ue
      real(rp),parameter:: ue2ev = real(1d0, rp)/ev2ue
      real(rp),parameter:: k2ue  = k2j*j2ue
      real(rp),parameter:: ue2k  = ue2j*j2k
!.....Plank constant in eV*fs (6.62607e-34 /1602e-19 /1e-15)
      real(rp),parameter:: plankh = real(4.135667733d0, rp)
