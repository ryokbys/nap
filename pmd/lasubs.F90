!=======================================================================
!     SUBROUTINE ludc_inv
!       Calculate INversion of the MATrix using LU decomposition.
!=======================================================================
SUBROUTINE ludc_inv( n,a,c )
  IMPLICIT NONE
!-----------------------------------------------------------------------
!     n:  dimension of the matrix A
!     a:  the matrix A
!-----------------------------------------------------------------------
!-----arguments
  integer,intent(in):: n
  real*8,intent(in):: a(n,n)
!-----return 
  real*8,intent(out):: c(n,n)
!-----locals
  integer i,j,k,pivot
  real*8 amax,tmp
  real*8,allocatable:: d(:,:)
!-----function
!      real*8 vbtime

!      write(6,'(A)') "---ludc_inv"
!      write(6,'(A,F12.3)') "TIME =",vbtime(4)

!     Initialize?
  allocate(d(n,n))
  DO i=1,n
    DO j=1,n
      d(j,i)=a(j,i)
    ENDDO
!     Initialize.
    DO j=1,n
      IF( j.EQ.i ) c(j,i)=1.D0
      IF( j.NE.i ) c(j,i)=0.D0
    ENDDO
  ENDDO

!     LU decompose of the matrix. O(N^3)
  DO i=1,n
!     Pivoting...
    amax=0
    pivot=i
    DO j=i,n
      IF( ABS(d(j,i)) .GT. amax ) THEN
        amax=ABS(d(j,i))
        pivot=j
      ENDIF
    ENDDO
    IF( pivot .NE. i ) THEN
      DO j=1,n
        tmp=d(i,j)
        d(i,j)=d(pivot,j)
        d(pivot,j)=tmp
        tmp=c(i,j)
        c(i,j)=c(pivot,j)
        c(pivot,j)=tmp
      ENDDO
    ENDIF
!     Gauss-Jordan elimination.
    DO j=i+1,n
      d(j,i) = d(j,i)/d(i,i)
      DO k=i+1,n
        d(j,k) = d(j,k) - d(i,k)*d(j,i)
      ENDDO
    ENDDO
  ENDDO
!      write(6,*) "d:"
!      do i=1,n
!        write(6,'(30D11.3)') (d(i,j),j=1,n)
!      enddo
!      write(6,'(A,F12.3)') "TIME =",vbtime(4)

!     Calculate inversion of the matrix. O(N^3)
  DO k=1,n
!     Using L matrix...
    DO i=1,n
      DO j=i+1,n
        c(j,k)=c(j,k)-c(i,k)*d(j,i)
      ENDDO
    ENDDO
    DO i=n,1,-1
      DO j=i+1,n
        c(i,k)=c(i,k)-d(i,j)*c(j,k)
      ENDDO
      c(i,k)=c(i,k)/d(i,i)
    ENDDO
  ENDDO
!      write(6,'(A,F12.3)') "TIME =",vbtime(4)
  deallocate(d)
  RETURN
END subroutine ludc_inv
!=======================================================================
subroutine chckinv(n,a,ai)
!-----------------------------------------------------------------------
!     subroutine chckinv
!-----------------------------------------------------------------------
  implicit none
  integer,intent(in):: n
  real*8,intent(in):: a(n,n),ai(n,n)
  integer i,j,k
  real*8 tmp
  real*8,allocatable:: e(:,:),dum(:,:)

  write(6,'(A)') "---chckinv"
  allocate(e(n,n),dum(n,n))
  e(1:n,1:n)=0.d0
  dum(1:n,1:n)=0.d0
  do i=1,n
    e(i,i)=1.d0
  enddo

  write(6,*) "a:"
  do j=1,min(n,6)
    write(6,'(6E11.3)') (a(i,j),i=1,min(n,6))
  enddo
  write(6,*) "ai:"
  do j=1,min(n,6)
    write(6,'(6E11.3)') (ai(i,j),i=1,min(n,6))
  enddo
  tmp=0.d0
  do i=1,n
    do j=1,n
      do k=1,n
        dum(j,i)=dum(j,i)+a(k,j)*ai(i,k)
      enddo
      if(abs(e(j,i)-dum(j,i)).gt.tmp) tmp=abs(e(j,i)-dum(j,i))
    enddo
  enddo
  write(6,*) "dum:"
  do j=1,min(n,6)
    write(6,'(6E11.3)') (dum(i,j),i=1,min(n,6))
  enddo
  write(6,'(A,E12.3)') " max diff =",tmp
  deallocate(e,dum)
end subroutine chckinv
!=======================================================================
subroutine choldc_inv(n,a,c)
!-----------------------------------------------------------------------
!     SUBROUTINE choldc_inv
!       Inversion routine using Cholesky decomposition with cholsl,
!       choldc, from Numerical Recipes in Fortran.
!-----------------------------------------------------------------------
  implicit real*8 (a-h,o-z)
  real*8,intent(in):: a(n,n)
  real*8,allocatable:: p(:),b(:),x(:),ar(:,:)
!      real*8 p(n),b(n),x(n),ar(n,n)
!     return
  real*8,intent(out) :: c(n,n)
!     functions
  real*8 vbtime

!      write(6,'(A)') "---choldc_inv"

  allocate(p(n),b(n),x(n),ar(n,n),stat=istat)
  c(1:n,1:n)=a(1:n,1:n)
  call choldc(c,n,p)

!      write(6,*) "c:"
!      do i=1,n
!        write(6,'(30D11.3)') (c(i,j),j=1,n)
!      enddo

  do i=1,n
    b(1:n)=0d0
    b(i)=1d0
    call cholsl(c,n,p,b,x)
    ar(1:n,i)=x(1:n)
  enddo

  c(1:n,1:n)=ar(1:n,1:n)
!      do i=1,n
!        do j=1,n
!          c(j,i)=ar(j,i)
!        enddo
!      enddo

  deallocate(p,b,x,ar)
  return
end subroutine choldc_inv
!=======================================================================
subroutine choldc(a,n,p)
!-----------------------------------------------------------------------
!     subroutine choldc: 
!       From Numerical Recipes in Fortran.
!       Lower triangle matrix is returned in a.
!-----------------------------------------------------------------------
  integer n
  real*8 a(n,n),p(n)
  integer i,j,k,l
  real*8 sum

  do i=1,n
    do j=i,n
      sum=a(i,j)
      do k=i-1,1,-1
!            if(i.eq.j) write(6,'(3I5,2E11.3)') i,j,k,sum,a(i,k)*a(j,k)
        sum=sum-a(i,k)*a(j,k)
      enddo
      if(i.eq.j)then
        if(sum.le.0d0) then
          write(6,'(A,I5,E11.3)') "Error (choldc): i,sum=",i,sum
          stop
        endif
        p(i)=sqrt(sum)
      else
        a(j,i)=sum/p(i)
      endif
!          if(i.eq.j) write(6,'(A,I5,E11.3)') "i,sum=",i,sum
    enddo
  enddo
  return
end subroutine choldc
!=======================================================================
subroutine cholsl(a,n,p,b,x)
  integer n
  real*8 a(n,n),b(n),p(n),x(n)
  integer i,k
  real*8 sum

  do i=1,n
    sum=b(i)
    do k=i-1,1,-1
      sum=sum-a(i,k)*x(k)
    enddo
    x(i)=sum/p(i)
  enddo
  do i=n,1,-1
    sum=x(i)
    do k=i+1,n
      sum=sum-a(k,i)*x(k)
    enddo
    x(i)=sum/p(i)
  enddo
  return
end subroutine cholsl
!=======================================================================
