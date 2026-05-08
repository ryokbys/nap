!-----Al,N mass (to be multiplied by umass)
      real(rp),parameter:: am_Al = real(26.9815d0, rp) *umass
      real(rp),parameter:: am_N  = real(14.0067d0, rp) *umass
!-----J*m/C**2 ---> eV*Ang/e**2
!      real(rp),parameter:: cfct= 8.9876d+9 *j2ev /ang
!     &     /(6.2414d+18)**2
!      real(rp),parameter:: cfct= 1d0
!-----Vashishta potential parameters for AlN
      real(rp),parameter:: v_z(1:2)= &
           (/ real(1.0366d0, rp), &
             -real(1.0366d0, rp) /)
      real(rp),parameter:: v_alp(1:2)= (/ real(0.0d0, rp), &
                                         real(3.0d0, rp) /)
      real(rp),parameter:: v_n(1:2,1:2)= &
           reshape([ real(7d0, rp), real(9d0, rp), real(9d0, rp), real(7d0, rp) ],[2,2])
      real(rp),parameter:: v_h(1:2,1:2)= reshape([&
           real(8.7162d-17, rp) *j2ev *real(1d0, rp)**v_n(1,1), &
           real(2.6824d-18, rp) *j2ev *real(1d0, rp)**v_n(2,1), &
           real(2.6824d-18, rp) *j2ev *real(1d0, rp)**v_n(1,2), &
!     &        2.6824d-17 *j2ev *1d0**v_n(2,1),
!     &        2.6824d-17 *j2ev *1d0**v_n(1,2), &
           real(2.1790d-17, rp) *j2ev *real(1d0, rp)**v_n(2,2)],[2,2])
      real(rp),parameter:: v_w(1:2,1:2)= reshape([ &
           real(0d0, rp), &
           real(9.7901d-18, rp) *j2ev *real(1d0, rp)**6, &
           real(9.7901d-18, rp) *j2ev *real(1d0, rp)**6, &
           real(0d0, rp)],[2,2])
      real(rp),parameter:: v_r1s(1:2,1:2)= reshape([&
           real(5d0, rp), &
           real(5d0, rp), &
           real(5d0, rp), &
           real(5d0, rp)],[2,2])
      real(rp),parameter:: v_r4s(1:2,1:2)= reshape([ &
           real(3.75d0, rp), &
           real(3.75d0, rp), &
           real(3.75d0, rp), &
           real(3.75d0, rp)],[2,2])
      real(rp),parameter:: v_b  = real(40d-19, rp) *j2ev
      real(rp),parameter:: v_c  = real(30d0, rp)
      real(rp),parameter:: v_xi = real(1d0, rp)
      real(rp),parameter:: v_r0 = real(2.6d0, rp)
      real(rp),parameter:: v_tht= real(109.47122d0, rp) /real(180d0, rp) *real(3.14159265358979d0, rp)
!-----num of division for 2-body term table
      integer,parameter:: nd2b = 2048
