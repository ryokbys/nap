!
!  Vector analysis routines
!
function sprod(ndim,a,b)
  implicit none
  integer,intent(in):: ndim
  real(8),intent(in):: a(ndim),b(ndim)
  real(8):: sprod
  integer:: i

  sprod = 0d0
  do i=1,ndim
    sprod = sprod +a(i)*b(i)
  enddo
  return
end function sprod
!=======================================================================
function absv(ndim,a)
  implicit none
  integer,intent(in):: ndim
  real(8),intent(in):: a(ndim)
  real(8):: absv
  integer:: i

  absv = 0d0
  do i=1,ndim
    absv = absv +a(i)*a(i)
  enddo
  absv = sqrt(absv)

  return
end function absv
!=======================================================================
subroutine vprod(a,b,ab)
!
!  vector product ab of vectors a and b, (ab = a x b)
!
  implicit none 
  real(8),intent(in):: a(3),b(3)
  real(8),intent(out):: ab(3)

  ab(1)= a(2)*b(3) -a(3)*b(2)
  ab(2)= a(3)*b(1) -a(1)*b(3)
  ab(3)= a(1)*b(2) -a(2)*b(1)
  return
end subroutine vprod

