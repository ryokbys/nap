!-----Si mass (to be multiplied by umass)
      real(rp),parameter:: am_si = real(28.0855d0, rp)
!-----SW unit energy in eV
      real(rp),parameter:: swe   = real(2.1678d0, rp)
!.....length scaling factor for matching this potential to VASP
      real(rp),parameter:: sfac  = real(1.0062662d0, rp)
!-----SW unit length in Ang
      real(rp),parameter:: swl   = real(2.0951d0, rp)*sfac
!-----si-si
      real(rp),parameter:: swa   = real(7.049556277d0, rp)
      real(rp),parameter:: swb   = real(0.6022245584d0, rp)
      real(rp),parameter:: swp   = real(4.d0, rp)
      real(rp),parameter:: swq   = real(0.d0, rp)
      real(rp),parameter:: swc   = real(1.d0, rp)
      real(rp),parameter:: swrc  = real(1.8d0, rp)
!-----si-si-si
      real(rp),parameter:: sws   = real(21.d0, rp)
      real(rp),parameter:: swt   = real(1.2d0, rp)
