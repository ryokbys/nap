!-----C
      real(rp),parameter:: am_c  = 12.0107
!-----interatomic length in graphene
      real(rp),parameter:: dcc   = real(1.44d-10, rp) /ang
!c-----interatomic length in diamond
!      real(rp),parameter:: dcc   = 1.605d-10 /ang
!-----Brenner Carbon (in eV,Ang)
      real(rp),parameter:: br_dd = real(6.325d0, rp)
      real(rp),parameter:: br_s  = real(1.29d0, rp)
      real(rp),parameter:: br_bet= real(1.5d+10, rp) *ang 
      real(rp),parameter:: br_re = real(1.315d-10, rp) /ang
      real(rp),parameter:: br_r1 = real(1.7d-10, rp) /ang
      real(rp),parameter:: br_r2 = real(2.0d-10, rp) /ang
      real(rp),parameter:: br_del= real(0.80469d0, rp)
      real(rp),parameter:: br_a0 = real(0.011304d0, rp)
      real(rp),parameter:: br_c0 = real(19.d0, rp)
      real(rp),parameter:: br_d0 = real(2.5d0, rp)
!-----interlayer potential by A.J.Heim et.al
      real(rp),parameter:: evdw  = real(5d-3, rp)
      real(rp),parameter:: cvdw  = real(200d0, rp)
      real(rp),parameter:: p3vdw =-real(2290.707325617024d0, rp)
      real(rp),parameter:: p4vdw = real(3603.929410034915d0, rp)
      real(rp),parameter:: p5vdw =-real(1513.930039501751d0, rp)
      real(rp),parameter:: sgmvdw= real(2.988d0, rp) *real(2d0, rp)**(real(1d0, rp)/real(6d0, rp))
      real(rp),parameter:: r2vdw = real(2d0, rp)**(-real(1d0, rp)/real(6d0, rp)) *sgmvdw
      real(rp),parameter:: r1vdw = real(0.683d0, rp) *r2vdw
      real(rp),parameter:: r3vdw = (real(13d0, rp)/real(7d0, rp))**(real(1d0, rp)/real(6d0, rp)) *sgmvdw
      real(rp),parameter:: r4vdw = real(67d0, rp)/real(48d0, rp) *r3vdw
      real(rp),parameter:: c0vdw = real(2.575275778983429d0, rp)
      real(rp),parameter:: c1vdw =-real(4.316677142326428d0, rp)
      real(rp),parameter:: c2vdw = real(1.376573417835169d0, rp)
      real(rp),parameter:: c3vdw =-real(0.12340088128894569d0, rp)
