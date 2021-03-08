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
  function matxvec3(mat,vec) result(vecout)
    real(8),intent(in):: mat(3,3),vec(3)
    real(8):: vecout(3)

    vecout(1:3) = mat(1:3,1)*vec(1) +mat(1:3,2)*vec(2) +mat(1:3,3)*vec(3)
    return
  end function matxvec3
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
  function norm(v)
    real(8):: v(:)
    real(8):: norm
    integer:: i

    norm = 0d0
    do i=1,ubound(v,1)
      norm = norm +v(i)*v(i)
    enddo
    norm = sqrt(norm)
    return
  end function norm
!=======================================================================
  function norm2(r)
    real(8):: r(:)
    real(8):: norm2
    integer:: i

    norm2 = 0d0
    do i=1,ubound(r,1)
      norm2 = norm2 +r(i)*r(i)
    enddo
    return
  end function norm2
!=======================================================================
  function matinv3(mat) result(mati)
    real(8),intent(in):: mat(3,3)
    real(8):: mati(3,3)

    integer:: i,j,im,ip,jm,jp
    real(8):: vol,sgm(3,3),matit(3,3)

!.....cofactor matrix, SGM
    do j=1,3
      jm=mod(j+1,3)+1
      jp=mod(j,  3)+1
      do i=1,3
        im=mod(i+1,3)+1
        ip=mod(i,  3)+1
        sgm(i,j)=mat(ip,jp)*mat(im,jm)-mat(im,jp)*mat(ip,jm)
      enddo
    enddo
!.....volume
    vol=mat(1,1)*sgm(1,1)+mat(2,1)*sgm(2,1)+mat(3,1)*sgm(3,1)
    do j=1,3
      do i=1,3
        matit(i,j)= sgm(i,j)/vol
      enddo
    enddo
!....transpose
    do j=1,3
      do i=1,3
        mati(i,j)= matit(j,i)
      enddo
    enddo
    return
  end function matinv3
end module vector
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd lib"
!     End:

