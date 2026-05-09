!-----C
      real(rp),parameter:: am_c  = 12.0107
!-----interatomic length in graphene
      real(rp),parameter:: dcc   = 1.44e-10_rp /ang
!c-----interatomic length in diamond
!      real(rp),parameter:: dcc   = 1.605d-10 /ang
!-----Brenner Carbon (in eV,Ang)
      real(rp),parameter:: br_dd = 6.325_rp
      real(rp),parameter:: br_s  = 1.29_rp
      real(rp),parameter:: br_bet= 1.5e+10_rp *ang 
      real(rp),parameter:: br_re = 1.315e-10_rp /ang
      real(rp),parameter:: br_r1 = 1.7e-10_rp /ang
      real(rp),parameter:: br_r2 = 2.0e-10_rp /ang
      real(rp),parameter:: br_del= 0.80469_rp
      real(rp),parameter:: br_a0 = 0.011304_rp
      real(rp),parameter:: br_c0 = 19._rp
      real(rp),parameter:: br_d0 = 2.5_rp
!-----interlayer potential by A.J.Heim et.al
      real(rp),parameter:: evdw  = 5e-3_rp
      real(rp),parameter:: cvdw  = 200.0_rp
      real(rp),parameter:: p3vdw =-2290.707325617024_rp
      real(rp),parameter:: p4vdw = 3603.929410034915_rp
      real(rp),parameter:: p5vdw =-1513.930039501751_rp
      real(rp),parameter:: sgmvdw= 2.988_rp *2.0_rp**(1.0_rp/6.0_rp)
      real(rp),parameter:: r2vdw = 2.0_rp**(-1.0_rp/6.0_rp) *sgmvdw
      real(rp),parameter:: r1vdw = 0.683_rp *r2vdw
      real(rp),parameter:: r3vdw = (13.0_rp/7.0_rp)**(1.0_rp/6.0_rp) *sgmvdw
      real(rp),parameter:: r4vdw = 67.0_rp/48.0_rp *r3vdw
      real(rp),parameter:: c0vdw = 2.575275778983429_rp
      real(rp),parameter:: c1vdw =-4.316677142326428_rp
      real(rp),parameter:: c2vdw = 1.376573417835169_rp
      real(rp),parameter:: c3vdw =-0.12340088128894569_rp
