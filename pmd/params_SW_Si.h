!-----Si mass (to be multiplied by umass)
      real(rp),parameter:: am_si = 28.0855_rp
!-----SW unit energy in eV
      real(rp),parameter:: swe   = 2.1678_rp
!.....length scaling factor for matching this potential to VASP
      real(rp),parameter:: sfac  = 1.0062662_rp
!-----SW unit length in Ang
      real(rp),parameter:: swl   = 2.0951_rp*sfac
!-----si-si
      real(rp),parameter:: swa   = 7.049556277_rp
      real(rp),parameter:: swb   = 0.6022245584_rp
      real(rp),parameter:: swp   = 4._rp
      real(rp),parameter:: swq   = 0._rp
      real(rp),parameter:: swc   = 1._rp
      real(rp),parameter:: swrc  = 1.8_rp
!-----si-si-si
      real(rp),parameter:: sws   = 21._rp
      real(rp),parameter:: swt   = 1.2_rp
