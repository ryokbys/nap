c.....Mass
      real(rp),parameter:: am_W   = real(183.84d0, rp)
      real(rp),parameter:: am_He  =   real(4.0026d0, rp)

c.....energy scaling
      real(rp),parameter:: sfac  = real(1d0, rp) !/2

c.....length scaling for hybrid QMCL calculation
      real(rp),parameter:: slen  = real(1d0, rp)/real(3.204d0, rp) *real(3.172d0, rp)

c.....2-body
c.....W-W
      real(rp),parameter:: p_WW_c    = real(3.25d0, rp) *slen
      real(rp),parameter:: p_WW_c0   = real(47.042512d0, rp) /slen**2
      real(rp),parameter:: p_WW_c1   = -real(33.987390d0, rp) /slen**3
      real(rp),parameter:: p_WW_c2   = real(6.319348d0, rp) /slen**4
      real(rp),parameter:: p_WW_B    = real(90.300865d0, rp) /slen**3
      real(rp),parameter:: p_WW_alpha= real(0.997985d0, rp) /slen
      real(rp),parameter:: p_WW_b0   = real(2.951811d0, rp) *slen
c.....W-He
      real(rp),parameter:: p_WHe_rc = real(4.0d0, rp)
      real(rp),parameter:: p_WHe_rs = real(2.2d0, rp)
      real(rp),parameter:: p_WHe_cs = p_WHe_rc -p_WHe_rs
      real(rp),parameter:: p_WHe_z1 = real(1488.31d0, rp)
      real(rp),parameter:: p_WHe_a1 = real(4.129d0, rp)
      real(rp),parameter:: p_WHe_b0 = -real(0.0033d0, rp)
c.....He-He
      real(rp),parameter:: p_HeHe_rc = real(2.4d0, rp)
      real(rp),parameter:: p_HeHe_rs = real(1.1d0, rp)
      real(rp),parameter:: p_HeHe_cs = p_HeHe_rc -p_HeHe_rs
      real(rp),parameter:: p_HeHe_z1 = real(87.8831d0, rp)
      real(rp),parameter:: p_HeHe_a1 = real(3.01517d0, rp)
      real(rp),parameter:: p_HeHe_b0 = real(0.0d0, rp)

c.....Many-body
c.....W
      real(rp),parameter:: p_W_w    = real(1.786367d0, rp) /slen
      real(rp),parameter:: p_W_d    = real(4.40d0, rp) *slen
      real(rp),parameter:: p_W_beta = -real(0.23433d0, rp)
c.....He
      real(rp),parameter:: p_He_A = -real(0.154896d0, rp)
      real(rp),parameter:: p_He_w = real(0.812281d0, rp)
      real(rp),parameter:: p_He_d = p_W_d
      real(rp),parameter:: p_He_beta = p_W_beta

