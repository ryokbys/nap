c
c     Vector analysis routines
c
      function sprod(a,b)
      implicit none
      real(8),intent(in):: a(3),b(3)
      real(8):: sprod
      
      sprod= a(1)*b(1) +a(2)*b(2) +a(3)*b(3)
      return
      end function sprod
c=======================================================================
      function absv(a)
      implicit none
      real(8),intent(in):: a(3)
      real(8):: absv

      absv= sqrt(a(1)*a(1) +a(2)*a(2) +a(3)*a(3))
      return
      end function absv
c=======================================================================
      subroutine vprod(a,b,ab)
c  vector product ab of vectors a and b, (ab = a x b)
      implicit none 
      real(8),intent(in):: a(3),b(3)
      real(8),intent(out):: ab(3)
      
      ab(1)= a(2)*b(3) -a(3)*b(2)
      ab(2)= a(3)*b(1) -a(1)*b(3)
      ab(3)= a(1)*b(2) -a(2)*b(1)
      return
      end subroutine vprod
