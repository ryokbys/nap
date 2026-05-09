c.....Mass
      real(rp),parameter:: am_W   = 183.84_rp
      real(rp),parameter:: am_He  =   4.0026_rp

c.....energy scaling
      real(rp),parameter:: sfac  = 1.0_rp !/2

c.....length scaling for hybrid QMCL calculation
      real(rp),parameter:: slen  = 1.0_rp/3.204_rp *3.172_rp

c.....2-body
c.....W-W
      real(rp),parameter:: p_WW_c    = 3.25_rp *slen
      real(rp),parameter:: p_WW_c0   = 47.042512_rp /slen**2
      real(rp),parameter:: p_WW_c1   = -33.987390_rp /slen**3
      real(rp),parameter:: p_WW_c2   = 6.319348_rp /slen**4
      real(rp),parameter:: p_WW_B    = 90.300865_rp /slen**3
      real(rp),parameter:: p_WW_alpha= 0.997985_rp /slen
      real(rp),parameter:: p_WW_b0   = 2.951811_rp *slen
c.....W-He
      real(rp),parameter:: p_WHe_rc = 4.0_rp
      real(rp),parameter:: p_WHe_rs = 2.2_rp
      real(rp),parameter:: p_WHe_cs = p_WHe_rc -p_WHe_rs
      real(rp),parameter:: p_WHe_z1 = 1488.31_rp
      real(rp),parameter:: p_WHe_a1 = 4.129_rp
      real(rp),parameter:: p_WHe_b0 = -0.0033_rp
c.....He-He
      real(rp),parameter:: p_HeHe_rc = 2.4_rp
      real(rp),parameter:: p_HeHe_rs = 1.1_rp
      real(rp),parameter:: p_HeHe_cs = p_HeHe_rc -p_HeHe_rs
      real(rp),parameter:: p_HeHe_z1 = 87.8831_rp
      real(rp),parameter:: p_HeHe_a1 = 3.01517_rp
      real(rp),parameter:: p_HeHe_b0 = 0.0_rp

c.....Many-body
c.....W
      real(rp),parameter:: p_W_w    = 1.786367_rp /slen
      real(rp),parameter:: p_W_d    = 4.40_rp *slen
      real(rp),parameter:: p_W_beta = -0.23433_rp
c.....He
      real(rp),parameter:: p_He_A = -0.154896_rp
      real(rp),parameter:: p_He_w = 0.812281_rp
      real(rp),parameter:: p_He_d = p_W_d
      real(rp),parameter:: p_He_beta = p_W_beta

