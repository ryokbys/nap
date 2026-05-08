!-----Si mass (to be multiplied by umass)
      real(rp),parameter:: am_si = real(28.0855d0, rp)
!!.....scaling factor for matching this pot to VASP
!      real(rp),parameter:: sfac  = 1.0116386d0
      real(rp),parameter:: sfac  = real(1d0, rp)
!-----EDIP Si
      real(rp),parameter:: ed_aa = real(7.9821730d0, rp)
      real(rp),parameter:: ed_bb = real(1.5075463d-10, rp) /ang *sfac
      real(rp),parameter:: ed_rho= real(1.2085196d0, rp)
      real(rp),parameter:: ed_a  = real(3.1213820d-10, rp) /ang *sfac
      real(rp),parameter:: ed_c  = real(2.5609104d-10, rp) /ang *sfac
      real(rp),parameter:: ed_sgm= real(0.5774108d-10, rp) /ang *sfac
      real(rp),parameter:: ed_lam= real(1.4533108d0, rp)
      real(rp),parameter:: ed_gam= real(1.1247945d-10, rp) /ang *sfac
      real(rp),parameter:: ed_eta= real(0.2523244d0, rp)
      real(rp),parameter:: ed_q0 = real(312.1341346d0, rp)
      real(rp),parameter:: ed_mu = real(0.6966326d0, rp)
      real(rp),parameter:: ed_bet= real(0.0070975d0, rp)
      real(rp),parameter:: ed_alp= real(3.1083847d0, rp)
      real(rp),parameter:: ed_u1 = -real(0.165799d0, rp)
      real(rp),parameter:: ed_u2 = real(32.557d0, rp)
      real(rp),parameter:: ed_u3 = real(0.286198d0, rp)
      real(rp),parameter:: ed_u4 = real(0.66d0, rp)
!-----2 lattice parameters
!      real(rp),parameter:: ratio = 0.94d0
      real(rp),parameter:: ratio = real(1.0d0, rp)

