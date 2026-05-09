!
!  Note:
!    unit of length = 1 Angstrom = 1e-10 meter
!    unit of time = 1e-15 second
!    unit of energy = 1 eV = 1.60217657e-19 Joule
!    unit of mass = unified atomic mass unit
!    unit of force = 1 eV/1Angstrom
!.....Pi
      real(rp),parameter:: pi = 3.14159265358979_rp
!-----time
      real(rp),parameter:: ut   = 1e-15_rp
      real(rp),parameter:: fs2s = ut
      real(rp),parameter:: s2fs = 1.0_rp/ut
      real(rp),parameter:: atu  = 2.418884326505e-17_rp ! atomic time unit
      real(rp),parameter:: s2atu = 1.0_rp/atu
      real(rp),parameter:: fs2atu = ut/atu
      real(rp),parameter:: atu2fs = 1.0_rp/fs2atu
!-----energy
      real(rp),parameter:: ev2j = 1.60217657e-19_rp
      real(rp),parameter:: j2ev = 1.0_rp/ev2j
      real(rp),parameter:: ht2ev = 27.2114_rp
      real(rp),parameter:: ev2ht = 1.0_rp/ht2ev
!-----length
      real(rp),parameter:: ang  = 1.0e-10_rp
      real(rp),parameter:: ang2m= ang
      real(rp),parameter:: m2ang= 1.0_rp/ang
      real(rp),parameter:: bohr = 0.52917721092e-10_rp
      real(rp),parameter:: ang2bohr = ang/bohr
      real(rp),parameter:: bohr2ang = 1.0_rp/ang2bohr
!.....pressure
      real(rp),parameter:: uprs   = ev2j/ang**3
      real(rp),parameter:: up2gpa = 1e-9_rp *uprs
      real(rp),parameter:: gpa2up = 1.0_rp /up2gpa
!-----mass of electron
      real(rp),parameter:: ume  = 9.1093897e-31_rp
!-----mass of 1/12 of Carbon
      real(rp),parameter:: ump  = 1.660538921e-27_rp
      real(rp),parameter:: ume2ump = ume/ump
      real(rp),parameter:: ump2ume = ump/ume
      real(rp),parameter:: ump2kg = ump
      real(rp),parameter:: kg2ump = 1.0_rp/ump
!-----unified atomic mass unit
      real(rp),parameter:: umass= ump
      real(rp),parameter:: amu2kg= umass
      real(rp),parameter:: kg2amu= 1.0_rp/umass
!-----Boltzmann factor in eV/K
      real(rp),parameter:: fkb  = 8.6173324e-5_rp
      real(rp),parameter:: k2ev = fkb
      real(rp),parameter:: ev2k = 1.0_rp /k2ev
      real(rp),parameter:: k2j  = k2ev *ev2j
      real(rp),parameter:: j2k  = j2ev *ev2k
!-----special energy unit = [ump*Ang**2/fs**2]
      real(rp),parameter:: j2ue  = kg2ump*(m2ang/s2fs)**2
      real(rp),parameter:: ue2j  = 1.0_rp/j2ue
      real(rp),parameter:: ev2ue = ev2j *j2ue
      real(rp),parameter:: ue2ev = 1.0_rp/ev2ue
      real(rp),parameter:: k2ue  = k2j*j2ue
      real(rp),parameter:: ue2k  = ue2j*j2k
!.....Plank constant in eV*fs (6.62607e-34 /1602e-19 /1e-15)
      real(rp),parameter:: plankh = 4.135667733_rp
