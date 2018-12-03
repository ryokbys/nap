module vector
  implicit none

contains
!=======================================================================
  function abc2cart(h,abc) result(cart)
    real(8),intent(in):: h(3,3),abc(3)
    real(8):: cart(3)

    cart(1:3) = h(1:3,1)*abc(1) &
         +h(1:3,2)*abc(2) &
         +h(1:3,3)*abc(3)
    return
  end function abc2cart
!=======================================================================
  function cart2abc(hi,cart) result(abc)
    real(8),intent(in):: hi(3,3),cart(3)
    real(8):: abc(3)

    abc(1:3) = hi(1:3,1)*cart(1) &
         +hi(1:3,2)*cart(2) &
         +hi(1:3,3)*cart(3)
    return
  end function cart2abc
!=======================================================================
  function dot(a,b)
    real(8):: a(:),b(:)
    real(8):: dot
    integer:: i
    dot = 0.0
    do i=1,ubound(a,1)
      dot = dot +a(i)*b(i)
    enddo
    return
  end function dot
!=======================================================================
  function cross(a,b) result(c)
    real(8):: a(3),b(3)
    real(8):: c(3)

    c(1) = a(2)*b(3) -b(2)*a(3)
    c(2) = a(3)*b(1) -b(3)*a(1)
    c(3) = a(1)*b(2) -b(1)*a(2)
    return
  end function cross
!=======================================================================
  function norm(r)
    real(8):: r(:)
    real(8):: norm
    integer:: i
    do i=1,ubound(r,1)
      norm = norm +r(i)*r(i)
    enddo
    return
  end function norm
end module vector
