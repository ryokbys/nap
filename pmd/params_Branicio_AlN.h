!-----Al,N mass (to be multiplied by umass)
      real(8),parameter:: am_Al = 26.9815d0 *umass
      real(8),parameter:: am_N  = 14.0067d0 *umass
!-----J*m/C**2 ---> eV*Ang/e**2
!      real(8),parameter:: cfct= 8.9876d+9 *j2ev /ang
!     &     /(6.2414d+18)**2
!      real(8),parameter:: cfct= 1d0
!-----Vashishta potential parameters for AlN
      real(8),parameter:: v_z(1:2)= &
           (/ 1.0366d0, &
             -1.0366d0 /)
      real(8),parameter:: v_alp(1:2)= (/ 0.0d0, &
                                         3.0d0 /)
      real(8),parameter:: v_n(1:2,1:2)= &
           reshape([ 7d0, 9d0, 9d0, 7d0 ],[2,2])
      real(8),parameter:: v_h(1:2,1:2)= reshape([&
           8.7162d-17 *j2ev *1d0**v_n(1,1), &
           2.6824d-18 *j2ev *1d0**v_n(2,1), &
           2.6824d-18 *j2ev *1d0**v_n(1,2), &
!     &        2.6824d-17 *j2ev *1d0**v_n(2,1),
!     &        2.6824d-17 *j2ev *1d0**v_n(1,2), &
           2.1790d-17 *j2ev *1d0**v_n(2,2)],[2,2])
      real(8),parameter:: v_w(1:2,1:2)= reshape([ &
           0d0, &
           9.7901d-18 *j2ev *1d0**6, &
           9.7901d-18 *j2ev *1d0**6, &
           0d0],[2,2])
      real(8),parameter:: v_r1s(1:2,1:2)= reshape([&
           5d0, &
           5d0, &
           5d0, &
           5d0],[2,2])
      real(8),parameter:: v_r4s(1:2,1:2)= reshape([ &
           3.75d0, &
           3.75d0, &
           3.75d0, &
           3.75d0],[2,2])
      real(8),parameter:: v_b  = 40d-19 *j2ev
      real(8),parameter:: v_c  = 30d0
      real(8),parameter:: v_xi = 1d0
      real(8),parameter:: v_r0 = 2.6d0
      real(8),parameter:: v_tht= 109.47122d0 /180d0 *3.14159265358979d0
!-----num of division for 2-body term table
      integer,parameter:: nd2b = 2048
