c-----Al and N mass (to be multiplied by amu)
      real(8),parameter:: am_Al = 26.9815d0 *amu
      real(8),parameter:: am_N  = 14.0067d0 *amu
c-----RK potential parameters for AlN
      real(8),parameter:: rk_z(1:2)=
     &     (/ 1.d0,
     &       -1.d0 /)
c-----energy-scale base
      real(8),parameter:: rk_eps= 0.01d0
c-----length-scale base
      real(8),parameter:: rk_sgm= 4d0
c-----2-body
      real(8),parameter:: rk_l1 = 0.25d0 *rk_sgm
      real(8),parameter:: rk_l2 = 2d0 *rk_sgm
c-----3-body
      real(8),parameter:: rk_b  = 540d0 *rk_eps
      real(8),parameter:: rk_l3 = 0.5d0 *rk_sgm
      real(8),parameter:: rk_rc = 1d0 *rk_sgm
      real(8),parameter:: rk_c  = 30d0
c-----num of division for 2-body term table
      integer,parameter:: nd2b = 2048
