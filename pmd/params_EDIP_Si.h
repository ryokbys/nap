!-----Si mass (to be multiplied by umass)
      real(rp),parameter:: am_si = 28.0855_rp
!!.....scaling factor for matching this pot to VASP
!      real(rp),parameter:: sfac  = 1.0116386d0
      real(rp),parameter:: sfac  = 1.0_rp
!-----EDIP Si
      real(rp),parameter:: ed_aa = 7.9821730_rp
      real(rp),parameter:: ed_bb = 1.5075463e-10_rp /ang *sfac
      real(rp),parameter:: ed_rho= 1.2085196_rp
      real(rp),parameter:: ed_a  = 3.1213820e-10_rp /ang *sfac
      real(rp),parameter:: ed_c  = 2.5609104e-10_rp /ang *sfac
      real(rp),parameter:: ed_sgm= 0.5774108e-10_rp /ang *sfac
      real(rp),parameter:: ed_lam= 1.4533108_rp
      real(rp),parameter:: ed_gam= 1.1247945e-10_rp /ang *sfac
      real(rp),parameter:: ed_eta= 0.2523244_rp
      real(rp),parameter:: ed_q0 = 312.1341346_rp
      real(rp),parameter:: ed_mu = 0.6966326_rp
      real(rp),parameter:: ed_bet= 0.0070975_rp
      real(rp),parameter:: ed_alp= 3.1083847_rp
      real(rp),parameter:: ed_u1 = -0.165799_rp
      real(rp),parameter:: ed_u2 = 32.557_rp
      real(rp),parameter:: ed_u3 = 0.286198_rp
      real(rp),parameter:: ed_u4 = 0.66_rp
!-----2 lattice parameters
!      real(rp),parameter:: ratio = 0.94d0
      real(rp),parameter:: ratio = 1.0_rp

