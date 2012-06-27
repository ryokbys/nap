c-----C
      real(8),parameter:: am_c  = 12.0107 *aump/aume
c-----interatomic length in graphene
      real(8),parameter:: dcc   = 1.44d-10 /bohr
cc-----interatomic length in diamond
c      real(8),parameter:: dcc   = 1.605d-10 /bohr
c-----Brenner Carbon (in eV,AA)
      real(8),parameter:: br_dd = 6.325d0 /ehrt
      real(8),parameter:: br_s  = 1.29d0
      real(8),parameter:: br_bet= 1.5d+10 *bohr 
      real(8),parameter:: br_re = 1.315d-10 /bohr
      real(8),parameter:: br_r1 = 1.7d-10 /bohr
      real(8),parameter:: br_r2 = 2.0d-10 /bohr
      real(8),parameter:: br_del= 0.80469d0
      real(8),parameter:: br_a0 = 0.011304d0
      real(8),parameter:: br_c0 = 19.d0
      real(8),parameter:: br_d0 = 2.5d0
c-----interlayer potential by A.J.Heim et.al
      real(8),parameter:: evdw  = 5d-3 /ehrt
      real(8),parameter:: cvdw  = 200d0
      real(8),parameter:: p3vdw =-2290.707325617024d0 *(bohr*1d10)**3
      real(8),parameter:: p4vdw = 3603.929410034915d0 *(bohr*1d10)**4
      real(8),parameter:: p5vdw =-1513.930039501751d0 *(bohr*1d10)**5
      real(8),parameter:: sgmvdw= 2.988d-10 *2d0**(1d0/6d0) /bohr
      real(8),parameter:: r2vdw = 2d0**(-1d0/6d0) *sgmvdw
      real(8),parameter:: r1vdw = 0.683d0 *r2vdw
      real(8),parameter:: r3vdw = (13d0/7d0)**(1d0/6d0) *sgmvdw
      real(8),parameter:: r4vdw = 67d0/48d0 *r3vdw
      real(8),parameter:: c0vdw = 2.575275778983429d0
      real(8),parameter:: c1vdw =-4.316677142326428d0 *(bohr*1d10)
      real(8),parameter:: c2vdw = 1.376573417835169d0 *(bohr*1d10)**2
      real(8),parameter:: c3vdw =-0.12340088128894569d0  *(bohr*1d10)**3
