!-----Al,N mass (to be multiplied by umass)
      real(rp),parameter:: am_Al = 26.9815_rp *umass
      real(rp),parameter:: am_N  = 14.0067_rp *umass
!-----J*m/C**2 ---> eV*Ang/e**2
!      real(rp),parameter:: cfct= 8.9876d+9 *j2ev /ang
!     &     /(6.2414d+18)**2
!      real(rp),parameter:: cfct= 1d0
!-----Vashishta potential parameters for AlN
      real(rp),parameter:: v_z(1:2)= &
           (/ 1.0366_rp, &
             -1.0366_rp /)
      real(rp),parameter:: v_alp(1:2)= (/ 0.0_rp, &
                                         3.0_rp /)
      real(rp),parameter:: v_n(1:2,1:2)= &
           reshape([ 7.0_rp, 9.0_rp, 9.0_rp, 7.0_rp ],[2,2])
      real(rp),parameter:: v_h(1:2,1:2)= reshape([&
           8.7162e-17_rp *j2ev *1.0_rp**v_n(1,1), &
           2.6824e-18_rp *j2ev *1.0_rp**v_n(2,1), &
           2.6824e-18_rp *j2ev *1.0_rp**v_n(1,2), &
!     &        2.6824d-17 *j2ev *1d0**v_n(2,1),
!     &        2.6824d-17 *j2ev *1d0**v_n(1,2), &
           2.1790e-17_rp *j2ev *1.0_rp**v_n(2,2)],[2,2])
      real(rp),parameter:: v_w(1:2,1:2)= reshape([ &
           0.0_rp, &
           9.7901e-18_rp *j2ev *1.0_rp**6, &
           9.7901e-18_rp *j2ev *1.0_rp**6, &
           0.0_rp],[2,2])
      real(rp),parameter:: v_r1s(1:2,1:2)= reshape([&
           5.0_rp, &
           5.0_rp, &
           5.0_rp, &
           5.0_rp],[2,2])
      real(rp),parameter:: v_r4s(1:2,1:2)= reshape([ &
           3.75_rp, &
           3.75_rp, &
           3.75_rp, &
           3.75_rp],[2,2])
      real(rp),parameter:: v_b  = 40e-19_rp *j2ev
      real(rp),parameter:: v_c  = 30.0_rp
      real(rp),parameter:: v_xi = 1.0_rp
      real(rp),parameter:: v_r0 = 2.6_rp
      real(rp),parameter:: v_tht= 109.47122_rp /180.0_rp *3.14159265358979_rp
!-----num of division for 2-body term table
      integer,parameter:: nd2b = 2048
