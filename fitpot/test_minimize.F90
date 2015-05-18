  program test_linmin
    implicit none
    integer,parameter:: iprint = 2
    real(8):: x0(10),g(10),a,b,c,fa,fb,fc,f
    real(8),external:: f1,f2
    real(8),parameter:: xtol = 1d-4
    real(8),parameter:: gtol = 1d-5
    real(8),parameter:: ftol = 1d-6
    interface
      function rosenbrock(n,x)
        integer,intent(in):: n
        real(8),intent(in):: x(n)
        real(8):: rosenbrock
      end function rosenbrock
      function drosen(n,x)
        integer,intent(in):: n
        real(8),intent(in):: x(n)
        real(8):: drosen(n)
      end function drosen
    end interface

    x0(1)= -2.0d0
    x0(2)=  1.0d0
!!$    f= rosenbrock(2,x0)
!!$    g(1:2)= drosen(2,x0)
!!$    print *,' x0=',x0(1:2)
!!$    print *,' f =',f
!!$    print *,' g =',g(1:2)
    call steepest_descent(2,x0,f,xtol,gtol,ftol,2000 &
         ,iprint,rosenbrock,drosen)
!!$    call bfgs(2,x0,f,xtol,gtol,ftol,2000 &
!!$         ,iprint,rosenbrock,drosen)
!!$    call cg(2,x0,f,xtol,gtol,ftol,2000 &
!!$         ,iprint,rosenbrock,drosen)
!!$    call bfgs(2,x0,f,xtol,gtol,ftol,f3,df3)
!!$    call quad_interpolate(1,x0,g,f,a,f1)
!!$    print *,' x0,a,f=',x0(1:2),a,f
  end program test_linmin
!=======================================================================
  function f1(n,x)
    implicit none
    integer,intent(in):: n
    real(8),intent(in):: x(n)
    real(8):: f1

    f1= x(1)*x(1) +2d0*exp(-x(1))
    return
  end function f1
!=======================================================================
  function f2(n,x)
    implicit none
    integer,intent(in):: n
    real(8),intent(in):: x(n)
    real(8):: f2

    f2= -x(1)*cos(x(1))
    return
  end function f2
!=======================================================================
  function rosenbrock(n,x)
    implicit none
    integer,intent(in):: n
    real(8),intent(in):: x(n)
    real(8):: rosenbrock
    
    rosenbrock= 100d0*(x(2)-x(1)*x(1))**2 +(1d0-x(1))*(1d0-x(1))
    return
  end function rosenbrock
!=======================================================================
  function drosen(n,x)
    implicit none
    integer,intent(in):: n
    real(8),intent(in):: x(n)
    real(8):: drosen(n)

    drosen(1)= -400d0 *x(1) *(x(2)-x(1)*x(1)) -2d0*(1d0-x(1))
    drosen(2)= 200d0 *(x(2)-x(1)*x(1))
    return
  end function drosen
!=======================================================================
  function f3(n,x)
    implicit none
    integer,intent(in):: n
    real(8),intent(in):: x(n)
    real(8):: f3
    f3= (x(1)+2*x(2)+7d0)**2 +(2*x(1)+x(2)-5)**2
    return
  end function f3
!=======================================================================
  function df3(n,x)
    implicit none
    integer,intent(in):: n
    real(8),intent(in):: x(n)
    real(8):: df3(n)
    df3(1)= 2d0*(x(1)+2*x(2)+7d0) +4d0*(2*x(1)+x(2)-5)
    df3(2)= 4d0*(x(1)+2*x(2)+7d0) +2d0*(2*x(1)+x(2)-5)
    return
  end function df3
!=======================================================================

