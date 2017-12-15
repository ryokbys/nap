module minimize
  save
!.....number of convergence criteria achieved
  integer:: numtol = 1
  
!.....penalty: lasso or ridge
  character(len=128):: cpena= 'none'
  character(len=128):: clinmin= 'armijo'
  real(8):: pwgt = 1d-15
!.....group lasso
  integer:: ngl
  integer,allocatable:: iglid(:)
  real(8),allocatable:: glval(:)
!.....group fs and mask
  integer,allocatable:: mskgfs(:),msktmp(:)
  integer:: nitergfs=100

!.....Armijo parameters
  real(8):: armijo_xi      = 1.0d-4
  real(8):: armijo_tau     = 0.5d0
  integer:: armijo_maxiter = 15

!.....Simulated annealing parameters
  real(8):: sa_temp0 = 1d0
  real(8):: sa_tau   = 10d0
  real(8):: sa_xw0   = 1d-3
  real(8):: sa_fctr  = 0.5d0
  real(8),allocatable:: sa_xws(:)
!.....T control method
!.....  - linear
!.....  - exp
!.....  - best
!.....  - constant
  character(len=128):: sa_tctrl = 'best'
  real(8):: sa_div_best = 10d0

!.....Metadynamics
  real(8):: md_height = 1d0
  real(8):: md_sigma  = 1d0
!.....Max num of gaussian potentials
  integer:: md_ng = 1000
  real(8),allocatable:: md_gp(:,:)

!.....CG
  integer:: icgbtype = 1 ! 1:FR, 2:PRP, 3:HS, 4:DY

!.....L-BFGS
  integer:: mstore   = 10

!.....Genetic Algorithm variables.......................................
  real(8):: ga_rate_mutate = 0.1
  integer:: ga_nindivs = 10
  integer:: ga_nbits = 16
  integer:: ga_noffsp = 0
  real(8):: ga_temp = 1d0
  integer:: ga_ngenes
  character(len=128):: ga_fitness = 'inv'
  
  type gene  ! A Gene corresponds to a parameter/variable to be optimized
    integer(2),allocatable:: bits(:)
    real(8):: val
    real(8):: vmax,vmin
  end type gene

  type individual  ! An individual is a set of parameters
    integer:: iid  ! ID for this individual
    type(gene),allocatable:: genes(:)
    real(8):: fvalue  ! Loss function value of the individual
    real(8):: fitness ! Fittness value of the individual
    real(8),allocatable:: vel(:)
  end type individual

!.....Differential evolution variables..................................
  character(len=128):: de_fitness = 'inv'
  character(len=128):: de_algo = 'local_neigh'
  integer:: de_nindivs = 10
  real(8):: de_frac    = 1d0
  real(8):: de_lambda  = -1d0
  real(8):: de_cross_rate = 0.5d0
  real(8):: de_wmin    = 0.4d0
  real(8):: de_wmax    = 0.8d0

!.....Particle swarm optimization variables.............................
  integer:: pso_nindivs = 10
  real(8):: pso_w       = 0.99d0
  real(8):: pso_c1      = 2d0
  real(8):: pso_c2      = 2d0
  real(8):: pso_vinimax = 0.1d0

  real(8):: fupper_lim = 1d+5
  real(8),allocatable:: ranges(:,:)

contains
!=======================================================================
  subroutine set_ranges(ndim,xranges)
    implicit none
    integer,intent(in):: ndim
    real(8),intent(in):: xranges(2,ndim)

    if( .not. allocated(ranges) ) allocate(ranges(2,ndim))
    if( size(ranges).ne.2*ndim ) then
      deallocate(ranges)
      allocate(ranges(2,ndim))
    endif
    return
  end subroutine set_ranges
!=======================================================================
  subroutine wrap_ranges(ndim,x)
    implicit none
    integer,intent(in):: ndim
    real(8),intent(inout):: x(ndim)

    integer:: i

    if( .not.allocated(ranges) ) then
      print *,'Error: ranges is not allocated yet...'
      stop
    endif

    do i=1,ndim
      if( x(i).lt.ranges(1,ndim) ) then
        x(i) = ranges(1,ndim)
      else if( x(i).gt.ranges(2,ndim) ) then
        x(i) = ranges(2,ndim)
      endif
    enddo
    return
  end subroutine wrap_ranges
!=======================================================================
  subroutine steepest_descent(ndim,x,f,g,d,xranges,xtol,gtol,ftol,maxiter &
       ,iprint,iflag,myid,func,grad)
    implicit none
    integer,intent(in):: ndim,maxiter,iprint,myid
    integer,intent(inout):: iflag
    real(8),intent(in):: xtol,gtol,ftol,xranges(2,ndim)
    real(8),intent(inout):: f,x(ndim),g(ndim),d(ndim)
!!$    real(8):: func,grad
    interface
      subroutine func(n,x,ftrn,ftst)
        integer,intent(in):: n
        real(8),intent(in):: x(n)
        real(8),intent(out):: ftrn,ftst
      end subroutine func
      subroutine grad(n,x,gtrn)
        integer,intent(in):: n
        real(8),intent(in):: x(n)
        real(8),intent(out):: gtrn(n)
      end subroutine grad
    end interface

    integer:: iter,i,niter
    real(8):: alpha,fp,gnorm,ftst

    if( .not.allocated(ranges) ) call set_ranges(ndim,xranges)

    iter= 0
    call wrap_ranges(ndim,x)
    call func(ndim,x,f,ftst)
    if( trim(cpena).eq.'lasso' .or. trim(cpena).eq.'LASSO' ) then
      do i=1,ndim
        f= f +pwgt*abs(x(i))
      enddo
    else if( trim(cpena).eq.'ridge' ) then
      do i=1,ndim
        f= f +pwgt*x(i)*x(i)
      enddo
    endif
    call grad(ndim,x,g)
    gnorm= sqrt(sprod(ndim,g,g))
    d(1:ndim)= -g(1:ndim)
!!$    g(1:ndim)= -g(1:ndim)/gnorm
!!$    gnorm= gnorm/ndim
    if( myid.eq.0 ) then
      if( iprint.eq.1 ) then
        write(6,'(a,i8,10es15.7)') ' iter,f,gnorm=' &
             ,iter,f,gnorm
        flush(6)
      else if( iprint.ge.2 ) then
        write(6,'(a,i8,10es15.7)') ' iter,x,f,gnorm=' &
             ,iter,x(1:ndim),f,gnorm
        flush(6)
      endif
    endif

    do iter=1,maxiter
      fp= f
!.....line minimization
      if( trim(clinmin).eq.'quadratic' ) then
        call quad_interpolate(ndim,x,d,f,ftst,xtol,gtol,ftol,alpha &
             ,iprint,iflag,myid,func)
!.....if quad interpolation failed, perform golden section
        if( iflag/100.ne.0 ) then
          iflag= iflag -(iflag/100)*100
          if(myid.eq.0) then
            print *,'since quad_interpolate failed, call golden_section.'
          endif
          call golden_section(ndim,x,d,f,ftst,xtol,gtol,ftol,alpha &
               ,iprint,iflag,myid,func)
        endif
      else if ( trim(clinmin).eq.'golden') then
        call golden_section(ndim,x,d,f,ftst,xtol,gtol,ftol,alpha &
             ,iprint,iflag,myid,func)
      else ! armijo (default)
        alpha= 1d0
        call armijo_search(ndim,x,d,f,ftst,g,alpha,iprint &
             ,iflag,myid,func,niter)
      endif
      if( iflag/100.ne.0 ) return
      x(1:ndim)= x(1:ndim) +alpha*d(1:ndim)
      call wrap_ranges(ndim,x)
      call func(ndim,x,f,ftst)
      if( trim(cpena).eq.'lasso' .or. trim(cpena).eq.'LASSO' ) then
        do i=1,ndim
          f= f +pwgt*abs(x(i))
        enddo
      else if( trim(cpena).eq.'ridge' ) then
        do i=1,ndim
          f= f +pwgt*x(i)*x(i)
        enddo
      endif
      call grad(ndim,x,g)
      gnorm= sqrt(sprod(ndim,g,g))
      d(1:ndim)= -g(1:ndim)
!!$      g(1:ndim)= -g(1:ndim)/gnorm
!!$      gnorm= gnorm/ndim
      if( myid.eq.0 ) then
        if( iprint.eq.1 ) then
          write(6,'(a,i8,10es15.7)') ' iter,f,gnorm=' &
               ,iter,f,gnorm
        else if( iprint.ge.2 ) then
          write(6,'(a,i8,10es15.7)') ' iter,x,f,gnorm=' &
               ,iter,x(1:ndim),f,gnorm
        endif
      endif
!.....check convergence 
!!$      if( abs(alpha).lt.xtol ) then
!!$        if( myid.eq.0 ) then
!!$          print *,'>>> SD converged wrt xtol'
!!$          write(6,'(a,2es15.7)') '   alpha,xtol=',alpha,xtol
!!$        endif
!!$        iflag= iflag +1
!!$        return
      if( gnorm.lt.gtol ) then
        if( myid.eq.0 ) then
          print *,'>>> SD converged wrt gtol'
          write(6,'(a,2es15.7)') '   gnorm,gtol=',gnorm,gtol
        endif
        iflag= iflag +2
        return
!!$      else if( abs(f-fp)/abs(fp).lt.ftol ) then
!!$        if( myid.eq.0 ) then
!!$          print *,'>>> Sd converged wrt ftol'
!!$          write(6,'(a,2es15.7)') '   f-fp,ftol=',abs(f-fp)/abs(fp),ftol
!!$        endif
!!$        iflag= iflag +3
!!$        return
      endif
    enddo
    
    return
  end subroutine steepest_descent
!=======================================================================
  subroutine cg(ndim,x,f,g,u,xranges,xtol,gtol,ftol,maxiter,iprint,iflag,myid &
       ,func,grad,cfmethod,niter_eval,sub_eval)
!
!  Conjugate gradient minimization
!
    implicit none
    integer,intent(in):: ndim,maxiter,iprint,myid,niter_eval
    integer,intent(inout):: iflag
    real(8),intent(in):: xtol,gtol,ftol,xranges(2,ndim)
    real(8),intent(inout):: f,x(ndim),g(ndim),u(ndim)
    character(len=*),intent(in):: cfmethod
    interface
      subroutine func(n,x,ftrn,ftst)
        integer,intent(in):: n
        real(8),intent(in):: x(n)
        real(8),intent(out):: ftrn,ftst
      end subroutine func
      subroutine grad(n,x,gtrn)
        integer,intent(in):: n
        real(8),intent(in):: x(n)
        real(8),intent(out):: gtrn(n)
      end subroutine grad
      subroutine sub_eval(iter)
        integer,intent(in):: iter
      end subroutine sub_eval
    end interface
    real(8),parameter:: xtiny  = 1d-14
    real(8),parameter:: t_DL   = 0.2d0
    logical:: ltwice = .false.

    integer:: i,iter,nftol,nxtol,niter
    real(8):: alpha,fp,gnorm,gnormp,beta,pval,sgnorm,ftst,dxnorm,unorm
    real(8),save,allocatable:: gpena(:),gp(:),y(:),xp(:),s(:),dx(:),uu(:)

    if( myid.eq.0 ) then
      print *,'entering CG routine...'
    endif

    if( .not.allocated(gpena) ) allocate(gpena(ndim),gp(ndim)&
         ,y(ndim),xp(ndim),s(ndim),dx(ndim),uu(ndim))

    if( .not.allocated(ranges) ) call set_ranges(ndim,xranges)

    iter= 0
    nftol= 0
    nxtol= 0
    gpena(1:ndim)= 0d0
    call wrap_ranges(ndim,x)
    call func(ndim,x,f,ftst)
    call grad(ndim,x,g)
!.....penalty
    call penalty(cpena,pwgt,ndim,f,g,pval,gpena,x)
    g(1:ndim)= g(1:ndim) +gpena(1:ndim)
    gnorm= sprod(ndim,g,g)
    sgnorm= sqrt(gnorm)
!!$    g(1:ndim)= g(1:ndim)/sqrt(gnorm)
!!$    gnorm= gnorm/ndim
    if( myid.eq.0 ) then
      if( iprint.eq.1 ) then
        write(6,'(a,i8,10es15.7)') ' iter,f,gnorm=' &
             ,iter,f,sgnorm
      else if( iprint.ge.2 ) then
        write(6,'(a,i8,10es15.7)') ' iter,x,f,gnorm=' &
             ,iter,f,sgnorm,x(1:5)
      endif
    endif
    u(1:ndim)= -g(1:ndim)

    call sub_eval(0)
    do iter=1,maxiter
      fp= f
      xp(1:ndim)= x(1:ndim)
!.....normalize u-vector only for line search
      unorm = sqrt(sprod(ndim,u,u))
      uu(1:ndim) = u(1:ndim)/unorm
!.....evaluate statistics at every niter_eval
      if( mod(iter,niter_eval).eq.0 ) &
           call sub_eval(iter)
!.....line minimization
      if( trim(clinmin).eq.'quadratic' ) then
        call quad_interpolate(ndim,x,uu,f,ftst,xtol,gtol,ftol &
             ,alpha,iprint,iflag,myid,func)
!.....if quad interpolation failed, perform golden section
        if( iflag/100.ne.0 ) then
          iflag= iflag -(iflag/100)*100
          if(myid.eq.0) then
            print *,'since quad_interpolate failed, call golden_section.'
          endif
          call golden_section(ndim,x,uu,f,ftst,xtol,gtol,ftol &
               ,alpha,iprint,iflag,myid,func)
        endif
      else if ( trim(clinmin).eq.'golden') then
        call golden_section(ndim,x,uu,f,ftst,xtol,gtol,ftol,alpha &
             ,iprint,iflag,myid,func)
      else ! armijo (default)
        alpha= 1d0
        call armijo_search(ndim,x,uu,f,ftst,g,alpha,iprint &
             ,iflag,myid,func,niter)
      endif

      if( iflag/100.ne.0 ) then
        if( ltwice ) then
          if(myid.eq.0) then
            print *,'>>> Line search failed twice continuously.'
          endif
          return
        else
          ltwice= .true.
          if(myid.eq.0) then
            print *,'>>> gg initialized because alpha was not found.'
          endif
          u(1:ndim)= -g(1:ndim)
          f= fp
          iflag= iflag -100*(iflag/100)
          cycle
        endif
      else
        ltwice=.false.
      endif

      if( trim(cpena).eq.'lasso' .or. trim(cpena).eq.'LASSO' ) then
        call soft_threshold(ndim,x,uu,alpha)
        call wrap_ranges(ndim,x)
      else if( trim(cpena).eq.'ridge' ) then
        x(1:ndim)= x(1:ndim) +alpha*uu(1:ndim)
        call wrap_ranges(ndim,x)
        do i=1,ndim
          pval= pval +pwgt*x(i)*x(i)
          gpena(i)= 2d0*pwgt*x(i)
        enddo
      else
        x(1:ndim)= x(1:ndim) +alpha*uu(1:ndim)
        call wrap_ranges(ndim,x)
      endif

      dx(1:ndim)= x(1:ndim) -xp(1:ndim)
      gnormp= gnorm
      gp(1:ndim)= g(1:ndim)
      call grad(ndim,x,g)
      g(1:ndim)= g(1:ndim) +gpena(1:ndim)
      gnorm= sprod(ndim,g,g)
      sgnorm= sqrt(gnorm)
      dxnorm= sqrt(sprod(ndim,dx,dx))
      if( icgbtype.eq.2 ) then
!.....Polak-Ribiere-Polyak (PRP)
        y(1:ndim)= g(1:ndim) -gp(1:ndim)
        beta= sprod(ndim,g,y)/gnormp
      else if( icgbtype.eq.3 ) then
!.....Hestenes-Stiefel (HS)
        y(1:ndim)= g(1:ndim) -gp(1:ndim)
        beta= sprod(ndim,g,y)/sprod(ndim,u,y)
      else if( icgbtype.eq.4 ) then
!.....Dai-Yuan (DY)
        y(1:ndim)= g(1:ndim) -gp(1:ndim)
        beta= gnorm/sprod(ndim,u,y)
      else if( icgbtype.eq.5 ) then
!.....Dai-Liao (DL)
        s(1:ndim)= x(1:ndim) -xp(1:ndim)
        y(1:ndim)= g(1:ndim) -gp(1:ndim)
        beta= (sprod(ndim,g,y) -t_DL*sprod(ndim,g,s))&
             /sprod(ndim,u,y)
      else ! including icgbtype == 1
!.....Fletcher-Reeves (FR)
        beta= gnorm/gnormp
      endif
      u(1:ndim)= -g(1:ndim) +beta*u(1:ndim)
      if( myid.eq.0 ) then
        if( iprint.eq.1 ) then
          write(6,'(a,i8,10es15.7)') ' iter,f,gnorm,beta=' &
               ,iter,f,sgnorm,beta
        else if( iprint.ge.2 ) then
          write(6,'(a,i8,10es15.7)') ' iter,f,gnorm,beta,x=' &
               ,iter,f,sgnorm,beta,x(1:5)
        endif
      endif
!.....check convergence 
      if( dxnorm.lt.xtol ) then
        if( myid.eq.0 ) then
          print *,'>>> CG converged wrt xtol'
          write(6,'(a,2es15.7)') '   dxnorm,xtol=',dxnorm,xtol
        endif
        iflag= iflag +2
        return
      else if( abs(f-fp).lt.ftol ) then
        nftol= nftol +1
        if( nftol.gt.numtol ) then
          if( myid.eq.0 ) then
            print *,'>>> CG converged wrt ftol'
            write(6,'(a,2es15.7)') '   f-fp,ftol=' &
                 ,abs(f-fp),ftol
          endif
          iflag= iflag +3
          return
        endif
      endif
    enddo
    
    return
  end subroutine cg
!=======================================================================
  subroutine qn(ndim,x0,f,g,u,xranges,xtol,gtol,ftol,maxiter &
       ,iprint,iflag,myid,func,grad,cfmethod,niter_eval,sub_eval)
!
!  Broyden-Fletcher-Goldfarb-Shanno type of Quasi-Newton method.
!
    implicit none
    integer,intent(in):: ndim,iprint,myid,maxiter,niter_eval
    integer,intent(inout):: iflag
    real(8),intent(in):: xtol,gtol,ftol,xranges(2,ndim)
    real(8),intent(inout):: f,x0(ndim),g(ndim),u(ndim)
    character(len=*),intent(in):: cfmethod
    interface
      subroutine func(n,x,ftrn,ftst)
        integer,intent(in):: n
        real(8),intent(in):: x(n)
        real(8),intent(out):: ftrn,ftst
      end subroutine func
      subroutine grad(n,x,gtrn)
        integer,intent(in):: n
        real(8),intent(in):: x(n)
        real(8),intent(out):: gtrn(n)
      end subroutine grad
      subroutine sub_eval(iter)
        integer,intent(in):: iter
      end subroutine sub_eval
    end interface
    real(8),parameter:: xtiny  = 1d-14
    logical:: ltwice = .false.
    real(8),save,allocatable:: gg(:,:),x(:),s(:),y(:),gp(:) &
         ,ggy(:),gpena(:),dx(:)
    real(8):: tmp1,tmp2,b,sy,syi,fp,alpha,gnorm,ynorm,vnorm,pval &
         ,sgnx,absx,estmem,ftst,unorm,dxnorm
    integer:: i,j,iter,nftol,ngtol,nxtol,ig,mem,niter

    if( .not.allocated(gg) ) then
      if(myid.eq.0) then
        print *,''
        print *, '******************************* QN(BFGS) '&
             //'*******************************'
        estmem = (ndim*ndim +ndim*6)*8
        mem= estmem/1000/1000
        if( mem.eq.0 ) then
          mem= estmem/1000
          write(6,'(a,i6,a)') ' memory for BFGS = ' &
               ,int(estmem/1000),' kB'
        else
          write(6,'(a,i6,a)') ' memory for BFGS = ' &
               ,int(estmem/1000/1000),' MB'
        endif
      endif
      allocate(gg(ndim,ndim),x(ndim),dx(ndim) &
         ,s(ndim),y(ndim),gp(ndim),ggy(ndim),gpena(ndim))
    endif

    if( .not.allocated(ranges) ) call set_ranges(ndim,xranges)


!.....initialize alpha (line minimization factor)
    alpha = 1d0

    nftol= 0
    ngtol= 0
    nxtol= 0
!.....initial G = I
    gg(1:ndim,1:ndim)= 0d0
    do i=1,ndim
      gg(i,i)= 1d0
    enddo
    call wrap_ranges(ndim,x0)
    call func(ndim,x0,f,ftst)
    call grad(ndim,x0,g)
!.....penalty
    pval= 0d0
    gpena(1:ndim)= 0d0
    if( trim(cpena).eq.'lasso' ) then
      do i=1,ndim
        absx= abs(x0(i))
        pval= pval +pwgt*absx
        sgnx= sign(1d0,x0(i))
        if( absx.gt.xtiny ) gpena(i)= pwgt*sgnx
      enddo
    else if( trim(cpena).eq.'glasso' ) then
      glval(0:ngl)= 0d0
      do i=1,ndim
        ig= iglid(i)
        if(ig.gt.0) glval(ig)= glval(ig) +x0(i)*x0(i)
      enddo
      glval(0)= 1d0
      do ig=1,ngl
        glval(ig)= sqrt(glval(ig))
        pval= pval +pwgt*glval(ig)
      enddo
      do i=1,ndim
        ig= iglid(i)
        if( ig.eq.0 ) then ! i is not in a group
          absx= abs(x0(i))
          sgnx= sign(1d0,x0(i))
          if( absx.gt.xtiny ) gpena(i)= pwgt*sgnx
          pval= pval +pwgt*absx
        else if( ig.gt.0 ) then ! i is in a group
          if( glval(ig).gt.xtiny) gpena(i)= pwgt*x0(i)/glval(ig)
        endif
      enddo
    else if( trim(cpena).eq.'ridge' ) then
      do i=1,ndim
        pval= pval +pwgt*x0(i)*x0(i)
        gpena(i)= 2d0*pwgt*x0(i)
      enddo
    endif
!!$    if( myid.eq.0 ) then
!!$      do i=1,ndim
!!$        write(6,'(a,i6,2es15.7)') 'i,g,gpena=',i,g(i),gpena(i)
!!$      enddo
!!$    endif
    g(1:ndim)= g(1:ndim) +gpena(1:ndim)

    gnorm= sqrt(sprod(ndim,g,g))
    x(1:ndim)= x0(1:ndim)
    call wrap_ranges(ndim,x)
    vnorm= sqrt(sprod(ndim,x,x))
    dxnorm = 0d0

    iter= 0
    if( myid.eq.0 ) then
      if( iprint.eq.1 ) then
        if( trim(cpena).eq.'lasso' .or. trim(cpena).eq.'glasso' &
             .or. trim(cpena).eq.'ridge' ) then
          write(6,'(a,i8,7es13.5)') &
               ' iter,ftrn,ftst,p,vnorm,gnorm,dxnorm,f-fp=',iter,f,ftst &
               ,pval,vnorm,gnorm,dxnorm,f
        else
          write(6,'(a,i8,7es13.5)') &
               ' iter,ftrn,ftst,vnorm,gnorm,dxnorm,f-fp=' &
               ,iter,f,ftst,vnorm,gnorm,dxnorm,f
        endif
        call flush(6)
      else if( iprint.ge.2 ) then
        if( trim(cpena).eq.'lasso' .or. trim(cpena).eq.'glasso' &
             .or. trim(cpena).eq.'ridge' ) then
          write(6,'(a,i8,15es13.5)') &
               ' iter,ftrn,ftst,p,vnorm,gnorm,dxnorm,x(1:5)=' &
               ,iter,f,ftst,pval,vnorm,gnorm,dxnorm,x(1:5)
        else
          write(6,'(a,i8,15es13.5)') &
               ' iter,frn,ftst,vnorm,gnorm,dxnorm,x(1:5)=' &
               ,iter,f,ftst,vnorm,gnorm,dxnorm,x(1:5)
        endif
        call flush(6)
      endif
    endif

    call sub_eval(0)
    do iter=1,maxiter
      u(1:ndim)= 0d0
      do i=1,ndim
        u(1:ndim)= u(1:ndim) -gg(1:ndim,i)*g(i)
      enddo
      unorm = sqrt(sprod(ndim,u,u))
      u(1:ndim) = u(1:ndim) /unorm
!!$      unorm = sqrt(sprod(ndim,u,u))
!!$      print *,'qn: iter,gnorm,unorm = ',iter,gnorm,unorm
!!$      print *,' u =',u(1:10)
!!$      print *,' g =',g(1:10)
!!$      print *,' gg=',gg(1:5,1:5)
!.....store previous func and grad values
      fp= f
      gp(1:ndim)= g(1:ndim)
!.....line minimization
      if( trim(clinmin).eq.'quadratic' ) then
        call quad_interpolate(ndim,x,u,f,ftst,xtol,gtol,ftol,alpha &
             ,iprint,iflag,myid,func)
!.....if quad interpolation failed, perform golden section
        if( iflag/100.ne.0 ) then
          iflag= iflag -(iflag/100)*100
          if(myid.eq.0) then
            print *,'since quad_interpolate failed, call golden_section.'
          endif
          call golden_section(ndim,x,u,f,ftst,xtol,gtol,ftol,alpha &
               ,iprint,iflag,myid,func)
        endif
      else if ( trim(clinmin).eq.'golden') then
        call golden_section(ndim,x,u,f,ftst,xtol,gtol,ftol,alpha &
             ,iprint,iflag,myid,func)
      else ! armijo (default)
!.....To enhance the convergence in Armijo search,
!.....use the history of previous alpha by multiplying 2
!.....avoiding constant decrease, but alpha should not be greater than 1.
        alpha = min(max(alpha,xtol*2d0)*2d0, 1d0)
        call armijo_search(ndim,x,u,f,ftst,g,alpha,iprint &
             ,iflag,myid,func,niter)
      endif
!!$      if(myid.eq.0) print *,'armijo steps, alpha=',niter,alpha
      if( iflag/100.ne.0 ) then
        if( ltwice ) then
          x0(1:ndim)= x(1:ndim)
          if(myid.eq.0) then
            print *,'>>> line_search failed twice continuously.'
          endif
          return
        else
          ltwice= .true.
          if(myid.eq.0) then
            print *,'>>> gg initialized because alpha was not found.'
          endif
          alpha= 1d0  ! reset alpha to 1
          gg(1:ndim,1:ndim)= 0d0
          do i=1,ndim
            gg(i,i)= 1d0
          enddo
          f= fp
          iflag= iflag -100*(iflag/100)
          cycle
        endif
      else
        ltwice=.false.
      endif
!.....evaluate statistics at every niter_eval
      if( mod(iter,niter_eval).eq.0 ) &
           call sub_eval(iter)
      pval= 0d0
      gpena(1:ndim)= 0d0
      if( trim(cpena).eq.'lasso' ) then
        call soft_threshold(ndim,x,u,alpha)
        call wrap_ranges(ndim,x)
        do i=1,ndim
          absx= abs(x(i))
          pval= pval +pwgt*absx
          sgnx= sign(1d0,x(i))
          if( absx.gt.xtiny ) gpena(i)= pwgt*sgnx
        enddo
      else if( trim(cpena).eq.'glasso' ) then
        call soft_threshold(ndim,x,u,alpha)
        call wrap_ranges(ndim,x)
        glval(0:ngl)= 0d0
        do i=1,ndim
          ig= iglid(i)
          if( ig.gt.0 ) glval(ig)= glval(ig) +x(i)*x(i)
        enddo
        glval(0)= 1d0
        do ig=1,ngl
          glval(ig)= sqrt(glval(ig))
          pval= pval +pwgt*glval(ig)
        enddo
        do i=1,ndim
          ig= iglid(i)
          if( ig.eq.0 ) then ! i is not in a group
            absx= abs(x(i))
            sgnx= sign(1d0,x(i))
            if( absx.gt.xtiny ) gpena(i)= pwgt*sgnx
            pval= pval +pwgt*absx
          else if( ig.gt.0 ) then ! i is in a group
            if( glval(ig).gt.xtiny) gpena(i)= pwgt*x(i)/glval(ig)
          endif
        enddo
      else if( trim(cpena).eq.'ridge' ) then
        x(1:ndim)= x(1:ndim) +alpha*u(1:ndim)
        call wrap_ranges(ndim,x)
        do i=1,ndim
          pval= pval +pwgt*x(i)*x(i)
          gpena(i)= 2d0*pwgt*x(i)
        enddo
      else
        x(1:ndim)= x(1:ndim) +alpha*u(1:ndim)
        call wrap_ranges(ndim,x)
      endif
      dx(1:ndim)= x(1:ndim) -x0(1:ndim)
      x0(1:ndim)= x(1:ndim)
      call grad(ndim,x,g)
      g(1:ndim)= g(1:ndim) +gpena(1:ndim)
      gnorm= sqrt(sprod(ndim,g,g))
      vnorm= sqrt(sprod(ndim,x,x))
      dxnorm= sqrt(sprod(ndim,dx,dx))
!!$      g(1:ndim)= g(1:ndim)/sqrt(gnorm)
!!$      gnorm= gnorm/ndim
      if( myid.eq.0 ) then
        if( iprint.eq.1 ) then
          if( trim(cpena).eq.'lasso' .or. trim(cpena).eq.'glasso' &
               .or.trim(cpena).eq.'ridge' ) then
            write(6,'(a,i8,7es13.5)') &
                 ' iter,ftrn,ftst,p,vnorm,gnorm,dxnorm,f-fp=',&
                 iter,f,ftst &
                 ,pval,vnorm,gnorm,dxnorm,f-fp
          else
            write(6,'(a,i8,6es13.5)') &
                 ' iter,ftrn,ftst,vnorm,gnorm,dxnorm,f-fp=' &
                 ,iter,f,ftst &
                 ,vnorm,gnorm,dxnorm,f-fp
          endif
          call flush(6)
        else if( iprint.ge.2 ) then
          if( trim(cpena).eq.'lasso' .or. trim(cpena).eq.'glasso' &
               .or. trim(cpena).eq.'ridge' ) then
            write(6,'(a,i8,15es13.5)') &
                 ' iter,ftrn,ftst,p,vnorm,gnorm,dxnorm,f-fp,x(1:5)=' &
                 ,iter,f,ftst,pval,vnorm,gnorm,dxnorm,f-fp,x(1:5)
          else
            write(6,'(a,i8,15es13.5)') &
                 ' iter,ftrn,ftst,vnorm,gnorm,dxnorm,f-fp,x(1:5)=' &
                 ,iter,f,ftst,vnorm,gnorm,dxnorm,f-fp,x(1:5)
          endif
          call flush(6)
        endif
      endif
!.....check convergence 
      if( dxnorm.lt.xtol ) then
        nxtol = nxtol +1
        ngtol = 0
        nftol = 0
        if( nxtol.ge.numtol ) then
          if( myid.eq.0 ) then
            print *,'>>> QN converged because of xtol over ' &
                 ,numtol,' times.'
            write(6,'(a,2es13.5)') '   dxnorm,xtol=',dxnorm,xtol
          endif
          x0(1:ndim)= x(1:ndim)
          iflag= iflag +1
          return
        endif
      else if( gnorm.lt.gtol ) then
        ngtol = ngtol +1
        nxtol = 0
        nftol = 0
        if( ngtol.ge.numtol ) then
          if( myid.eq.0 ) then
            print *,'>>> QN converged because of gtol over ' &
                 ,numtol,' times.'
            write(6,'(a,2es13.5)') '   gnorm,gtol=',gnorm,gtol
          endif
          x0(1:ndim)= x(1:ndim)
          iflag= iflag +2
          return
        endif
      else if( abs(f-fp).lt.ftol) then
        nftol= nftol +1
        nxtol = 0
        ngtol = 0
        if( nftol.ge.numtol ) then
          if( myid.eq.0 ) then
            print *,'>>> QN converged because of ftol over ' &
                 ,numtol,' times.'
            write(6,'(a,2es13.5)') '   abs(f-fp),ftol=',abs(f-fp), ftol
          endif
          x0(1:ndim)= x(1:ndim)
          iflag= iflag +3
          return
        endif
      else
        nxtol = 0
        ngtol = 0
        nftol = 0
      endif
      
      s(1:ndim)= alpha *u(1:ndim)
      y(1:ndim)= g(1:ndim) -gp(1:ndim)
      ynorm= sprod(ndim,y,y)
      if( ynorm.lt.1d-14 .or. dxnorm.lt.xtol .or. gnorm.lt.gtol &
           .or. abs(f-fp).lt.ftol ) then
        if(myid.eq.0) then
          print *,'>>> gg initialized'
        endif
        gg(1:ndim,1:ndim)= 0d0
        do i=1,ndim
          gg(i,i)= 1d0
        enddo
        cycle
      endif

!.....update matrix gg
      sy= sprod(ndim,s,y)
      syi= 1d0/sy
      do i=1,ndim
        tmp1= 0d0
        tmp2= 0d0
        do j=1,ndim
          tmp1= tmp1 +gg(j,i)*y(j)
        enddo
        ggy(i)= tmp1 *syi
      enddo
      b= 1d0
      do i=1,ndim
        b=b +y(i)*ggy(i)
      enddo
      b= b*syi
!.....without temporary matrix aa
      do j=1,ndim
        do i=1,ndim
          gg(i,j)=gg(i,j) +s(j)*s(i)*b &
               -(s(i)*ggy(j) +ggy(i)*s(j))
        enddo
      enddo
    enddo

    
!!$    if( myid.eq.0 ) print *,'maxiter exceeded in qn'
!!$    iflag= iflag +10
    x0(1:ndim)= x(1:ndim)
    return
  end subroutine qn
!=======================================================================
  subroutine lbfgs(ndim,x0,f,g,u,xranges,xtol,gtol,ftol,maxiter &
       ,iprint,iflag,myid,func,grad,cfmethod,niter_eval,sub_eval)
!
!  Limited memory BFGS(Broyden-Fletcher-Goldfarb-Shanno).
!  See, https://en.wikipedia.org/wiki/Limited-memory_BFGS
!
    implicit none
    integer,intent(in):: ndim,iprint,myid,maxiter,niter_eval
    integer,intent(inout):: iflag
    real(8),intent(in):: xtol,gtol,ftol,xranges(2,ndim)
    real(8),intent(inout):: f,x0(ndim),g(ndim),u(ndim)
    character(len=*),intent(in):: cfmethod
    interface
      subroutine func(n,x,ftrn,ftst)
        integer,intent(in):: n
        real(8),intent(in):: x(n)
        real(8),intent(out):: ftrn,ftst
      end subroutine func
      subroutine grad(n,x,gtrn)
        integer,intent(in):: n
        real(8),intent(in):: x(n)
        real(8),intent(out):: gtrn(n)
      end subroutine grad
      subroutine sub_eval(iter)
        integer,intent(in):: iter
      end subroutine sub_eval
    end interface
    real(8),parameter:: xtiny  = 1d-14
    logical:: ltwice = .false.
!!$    real(8),external:: sprod
    real(8),save,allocatable:: x(:),s(:,:),y(:,:)&
         ,gp(:),xp(:),gpena(:),a(:),rho(:)
    real(8):: tmp1,tmp2,dsy,dsyi,fp,alpha,gnorm,ynorm,pval,sgnx,absx&
         ,beta,estmem,ftst
    integer:: i,j,k,l,m,n,iter,nftol,ig,mem,niter

    if( .not.allocated(x) ) then
      if(myid.eq.0) then
        print *, ''
        print *, '**** L-BFGS starts ****'
        print *,'history length in L-BFGS =',mstore
        estmem= (ndim*(mstore+1)*2 +ndim*4 +(mstore+1)*2)*8
        mem= estmem/1000/1000
        if( mem.eq.0 ) then
          mem= estmem/1000
          print *,'memory for L-BFGS =',mem,'kB'
        else
          print *,'memory for L-BFGS =',mem,'MB'
        endif
      endif
      allocate(x(ndim) &
         ,s(ndim,0:mstore),y(ndim,0:mstore) &
         ,gp(ndim),gpena(ndim),xp(ndim) &
         ,a(0:mstore),rho(0:mstore))
    endif

    if( .not.allocated(ranges) ) call set_ranges(ndim,xranges)

    s(1:ndim,0:mstore)= 0d0
    y(1:ndim,0:mstore)= 0d0
    rho(0:mstore)= 0d0
    a(0:mstore)= 0d0

    nftol= 0
    call wrap_ranges(ndim,x0)
!.....initial G = I
    call func(ndim,x0,f,ftst)
    call grad(ndim,x0,g)
!.....penalty
    call penalty(cpena,pwgt,ndim,f,g,pval,gpena,x0)
    g(1:ndim)= g(1:ndim) +gpena(1:ndim)

    gnorm= sqrt(sprod(ndim,g,g))
    x(1:ndim)= x0(1:ndim)

    iter= 0
    if( myid.eq.0 ) then
      if( iprint.eq.1 ) then
        if( trim(cpena).eq.'lasso' .or. trim(cpena).eq.'glasso' &
             .or. trim(cpena).eq.'ridge' ) then
          write(6,'(a,i8,3es13.5)') ' iter,f,p,gnorm=',iter,f &
               ,pval,gnorm
        else
          write(6,'(a,i8,2es13.5)') ' iter,f,gnorm=',iter,f,gnorm
        endif
        call flush(6)
      else if( iprint.ge.2 ) then
        if( trim(cpena).eq.'lasso' .or. trim(cpena).eq.'glasso' &
             .or. trim(cpena).eq.'ridge' ) then
          write(6,'(a,i8,12es13.5)') ' iter,f,p,gnorm,x(1:5)=' &
               ,iter,f,pval,gnorm,x(1:5)
        else
          write(6,'(a,i8,12es13.5)') ' iter,f,gnorm,x(1:5)=' &
               ,iter,f,gnorm,x(1:5)
        endif
        call flush(6)
      endif
    endif

    call sub_eval(0)
    u(1:ndim)= -g(1:ndim)
    do iter=1,maxiter
!.....store previous func and grad values
      fp= f
      gp(1:ndim)= g(1:ndim)
      xp(1:ndim)= x(1:ndim)
!.....line minimization
      if( trim(clinmin).eq.'quadratic' ) then
        call quad_interpolate(ndim,x,u,f,ftst,xtol,gtol,ftol,alpha &
             ,iprint,iflag,myid,func)
!.....if quad interpolation failed, perform golden section
        if( iflag/100.ne.0 ) then
          iflag= iflag -(iflag/100)*100
          if(myid.eq.0) then
            print *,'since quad_interpolate failed, call golden_section.'
          endif
          call golden_section(ndim,x,u,f,ftst,xtol,gtol,ftol,alpha &
               ,iprint,iflag,myid,func)
        endif
      else if ( trim(clinmin).eq.'golden') then
        call golden_section(ndim,x,u,f,ftst,xtol,gtol,ftol,alpha &
             ,iprint,iflag,myid,func)
      else ! armijo (default)
        alpha= 1d0
!!$        write(6,'(a,10es11.3)') 'u before armijo=',u(1:10)
        call armijo_search(ndim,x,u,f,ftst,g,alpha,iprint &
             ,iflag,myid,func,niter)
      endif

!!$      if(myid.eq.0) print *,'alpha=',alpha
      if( iflag/100.ne.0 ) then
        if( ltwice ) then
          x0(1:ndim)= x(1:ndim)
          if(myid.eq.0) then
            print *,'   line_search failed twice continuously...'
          endif
          return
        else
          ltwice= .true.
          if(myid.eq.0) then
            print *,'   gg initialized because alpha was not found.'
          endif
          u(1:ndim)= -g(1:ndim)
          f= fp
          iflag= iflag -100*(iflag/100)
          cycle
        endif
      else
        ltwice=.false.
      endif
!.....evaluate statistics at every niter_eval
      if( mod(iter,niter_eval).eq.0 ) &
           call sub_eval(iter)
      pval= 0d0
      gpena(1:ndim)= 0d0
      if( trim(cpena).eq.'lasso' ) then
        call soft_threshold(ndim,x,u,alpha)
        call wrap_ranges(ndim,x)
        do i=1,ndim
          absx= abs(x(i))
          pval= pval +pwgt*absx
          sgnx= sign(1d0,x(i))
          if( absx.gt.xtiny ) gpena(i)= pwgt*sgnx
        enddo
      else if( trim(cpena).eq.'glasso' ) then
        call soft_threshold(ndim,x,u,alpha)
        call wrap_ranges(ndim,x)
        glval(0:ngl)= 0d0
        do i=1,ndim
          ig= iglid(i)
          if( ig.gt.0 ) glval(ig)= glval(ig) +x(i)*x(i)
        enddo
        glval(0)= 1d0
        do ig=1,ngl
          glval(ig)= sqrt(glval(ig))
          pval= pval +pwgt*glval(ig)
        enddo
        do i=1,ndim
          ig= iglid(i)
          if( ig.eq.0 ) then ! i is not in a group
            absx= abs(x(i))
            sgnx= sign(1d0,x(i))
            if( absx.gt.xtiny ) gpena(i)= pwgt*sgnx
            pval= pval +pwgt*absx
          else if( ig.gt.0 ) then ! i is in a group
            if( glval(ig).gt.xtiny) gpena(i)= pwgt*x(i)/glval(ig)
          endif
        enddo
      else if( trim(cpena).eq.'ridge' ) then
        x(1:ndim)= x(1:ndim) +alpha*u(1:ndim)
        call wrap_ranges(ndim,x)
        do i=1,ndim
          pval= pval +pwgt*x(i)*x(i)
          gpena(i)= 2d0*pwgt*x(i)
        enddo
      else
        x(1:ndim)= x(1:ndim) +alpha*u(1:ndim)
        call wrap_ranges(ndim,x)
      endif
      x0(1:ndim)= x(1:ndim)
      call grad(ndim,x,g)
      g(1:ndim)= g(1:ndim) +gpena(1:ndim)
      gnorm= sqrt(sprod(ndim,g,g))
      if( myid.eq.0 ) then
        if( iprint.eq.1 ) then
          if( trim(cpena).eq.'lasso' .or. trim(cpena).eq.'glasso' &
               .or.trim(cpena).eq.'ridge' ) then
            write(6,'(a,i8,4es13.5)') ' iter,f,p,gnorm,f-fp=',iter,f &
                 ,pval,gnorm,f-fp
          else
            write(6,'(a,i8,3es13.5)') ' iter,f,gnorm,f-fp=',iter,f &
                 ,gnorm,f-fp
          endif
          call flush(6)
        else if( iprint.ge.2 ) then
          if( trim(cpena).eq.'lasso' .or. trim(cpena).eq.'glasso' &
               .or. trim(cpena).eq.'ridge' ) then
            write(6,'(a,i8,13es13.5)') ' iter,f,p,gnorm,f-fp,x(1:5)=' &
                 ,iter,f,pval,gnorm,f-fp,x(1:5)
          else
            write(6,'(a,i8,13es13.5)') ' iter,f,gnorm,f-fp,x(1:5)=' &
                 ,iter,f,gnorm,f-fp,x(1:5)
          endif
          call flush(6)
        endif
      endif
!.....check convergence 
      if( gnorm.lt.gtol ) then
        if( myid.eq.0 ) then
          print *,'>>> QN converged wrt gtol'
          write(6,'(a,2es13.5)') '   gnorm,gtol=',gnorm,gtol
        endif
        x0(1:ndim)= x(1:ndim)
        iflag= iflag +2
        return
      else if( abs(f-fp)/abs(fp).lt.ftol) then
        nftol= nftol +1
        if( nftol.gt.10 ) then
          if( myid.eq.0 ) then
            print *,'>>> QN may be converged because of ftol ' // &
                 'over 10 times.'
!!$            write(6,'(a,2es13.5)') '   f-fp/fp,ftol=' &
!!$                 ,abs(f-fp)/abs(fp),ftol
          endif
          x0(1:ndim)= x(1:ndim)
          iflag= iflag +3
          return
        endif
      endif
      
!.....limited-memory BFGS procedure
!.....compute rho,s,y
      s(1:ndim,0)= x(1:ndim) -xp(1:ndim)
      y(1:ndim,0)= g(1:ndim) -gp(1:ndim)
      dsy= sprod(ndim,s(1,0),y(1,0))
      rho(0)= 1d0/dsy
      ynorm= sprod(ndim,y(1,0),y(1,0))
!.....shift the history values
      do m=mstore-1,0,-1
        rho(m+1)= rho(m)
        s(1:ndim,m+1)= s(1:ndim,m)
        y(1:ndim,m+1)= y(1:ndim,m)
      enddo
      if( ynorm.lt.1d-14 ) then
        if(myid.eq.0) then
          print *,'>>> gg initialized because y*y < 1d-14'
        endif
        u(1:ndim)= -g(1:ndim)
        cycle
      endif

      if( iter.eq.1 ) then
        u(1:ndim)= -g(1:ndim)
        cycle
      endif
      u(1:ndim)= -g(1:ndim)
      do m=1,mstore
        a(m)= rho(m)*sprod(ndim,s(1,m),u)
        u(1:ndim)= u(1:ndim) -a(m)*y(1:ndim,m)
      enddo
      m=min(mstore,iter)
      dsy= sprod(ndim,y(1,m),s(1,m))
      ynorm= sprod(ndim,y(1,m),y(1,m))
      u(1:ndim)= dsy/ynorm *u(1:ndim)
      do m=mstore,1,-1
        beta= rho(m)*sprod(ndim,y(1,m),u)
        u(1:ndim)= u(1:ndim) +s(1:ndim,m)*(a(m)-beta)
      enddo

    enddo
    
!!$    if( myid.eq.0 ) print *,'maxiter exceeded in qn'
!!$    iflag= iflag +10
    x0(1:ndim)= x(1:ndim)
    return
  end subroutine lbfgs
!=======================================================================
  subroutine get_bracket(ndim,x0,d,a,b,c,fa,fb,fc,fta,ftb,ftc &
       ,iprint,iflag,myid,func)
    implicit none
    integer,intent(in):: ndim,iprint,myid
    integer,intent(inout):: iflag
    real(8),intent(in):: x0(ndim),d(ndim)
    real(8),intent(inout):: a,b,fa,fb,fta,ftb
    real(8),intent(out):: c,fc,ftc
    
    interface
      subroutine func(n,x,ftrn,ftst)
        integer,intent(in):: n
        real(8),intent(in):: x(n)
        real(8),intent(out):: ftrn,ftst
      end subroutine func
    end interface

    real(8),parameter:: RATIO = 1.61803398875d0
    real(8),parameter:: RATIOI= 1d0/RATIO
    real(8),parameter:: TINY= 1d-12
    real(8),parameter:: GLIMIT= 100d0
    real(8),parameter:: MAXITER= 50
    real(8):: dum,r,q,u,ulim,fu,ftst
    integer:: iter

    call func(ndim,x0+a*d,fa,fta)
    call func(ndim,x0+b*d,fb,ftb)
    iter= 0
10  continue
    iter= iter +1
    if( iter.gt.MAXITER ) then
      if( myid.eq.0 ) then
        print *,'[Error] iter.gt.MAXITER in get_bracket'
        print *,'  Search direction may not be a descent direction.'
      endif
      iflag= iflag +1000
      return
!!$      stop
    endif
    if( abs(b-a).lt.1d-12) then
      if( myid.eq.0 ) then
        print *,'[Error] a and b is too close in get_bracket'
        print *,'  Search direction may not be a descent direction.'
      endif
      iflag= iflag +2000
      return
!!$      stop
    endif
    if( fa.lt.fb ) then
      c= a +RATIOI*(b-a)
      call func(ndim,x0+c*d,fc,ftc)
      call exchange(c,b)
      call exchange(fc,fb)
      call exchange(ftc,ftb)
      if( iprint.eq.3 .and. myid.eq.0 ) then
        write(6,'(a,2(1x,3es12.4))') ' a,b,c,fa,fb,fc=',a,b,c,fa,fb,fc
      endif
      goto 10
    else
      c= a +RATIO*(b-a)
      call func(ndim,x0+c*d,fc,ftc)
      if( fb.gt.fc ) then
        b= a +RATIO*(c-a)
        call func(ndim,x0+b*d,fb,ftb)
        call exchange(b,c)
        call exchange(fb,fc)
        call exchange(ftb,ftc)
        if( iprint.eq.3 .and. myid.eq.0 ) then
          write(6,'(a,2(1x,3es12.4))') ' a,b,c,fa,fb,fc=',a,b,c,fa,fb,fc
        endif
        goto 10
      endif
!!$      write(6,'(a,2(1x,3es12.4))') ' a,b,c,fa,fb,fc=',a,b,c,fa,fb,fc
    endif

    return
  end subroutine get_bracket
!=======================================================================
  subroutine exchange(a,b)
    implicit none
    real(8),intent(inout):: a,b
    real(8):: tmp
    tmp= a
    a= b
    b= tmp
    return
  end subroutine exchange
!=======================================================================
  subroutine quad_interpolate(ndim,x0,g,f,ftst,xtol,gtol,ftol,a,iprint &
       ,iflag,myid,func)
    implicit none
    integer,intent(in):: ndim,iprint,myid
    integer,intent(inout):: iflag

    real(8),intent(in):: x0(ndim),xtol,gtol,ftol,g(ndim)
    real(8),intent(out):: f,a,ftst

    real(8),parameter:: STP0    = 1d-1
    real(8),parameter:: STPMAX  = 1d+1
    real(8),parameter:: TINY    = 1d-15
    integer,parameter:: MAXITER = 100

    interface
      subroutine func(n,x,ftrn,ftst)
        integer,intent(in):: n
        real(8),intent(in):: x(n)
        real(8),intent(out):: ftrn,ftst
      end subroutine func
    end interface

    integer:: iter,imin,imax,ix
    real(8):: r,q,fmin,fmax,dmin,dmax,d,xmin
    real(8),save,allocatable:: xi(:),fi(:),fti(:)

    if( .not. allocated(xi) ) allocate(xi(4),fi(4),fti(4))
    
    xi(1)= 0d0
    xi(2)= xi(1) +STP0
    call get_bracket(ndim,x0,g,xi(1),xi(2),xi(3),fi(1),fi(2),fi(3) &
         ,fti(1),fti(2),fti(3) &
         ,iprint,iflag,myid,func)
    if( iflag/1000.ne.0 ) return
    
    iter= 0
10  continue
    iter= iter +1
    if( iter.gt.MAXITER ) then
      if( myid.eq.0 ) then
        print *,' [Error] iter.gt.MAXITER in quad_interpolate !!!'
        print *,'   iter,MAXITER= ',iter,MAXITER
      endif
      iflag= iflag +100
      return
!!$      stop
    endif
    !.....step 3; compute turning point
    r= (xi(2)-xi(1))*(fi(2)-fi(3))
    q= (xi(2)-xi(3))*(fi(2)-fi(1))
    xi(4)= xi(2) -((xi(2)-xi(3))*q -(xi(2)-xi(1))*r) &
         /(2d0*sign(max(abs(q-r),TINY),q-r))
    call func(ndim,x0+xi(4)*g,fi(4),fti(4))
!!$    write(6,'(a,2(2x,4f11.2))') ' xi,fi=',xi(1:4),fi(1:4)

    !.....step4
    fmin= min(fi(1),fi(2),fi(3))
    fmax= max(fi(1),fi(2),fi(3))
    dmin= min(abs(xi(4)-xi(1)),abs(xi(4)-xi(2)),abs(xi(4)-xi(3)))
    if( fi(4).lt.fmin .and. dmin.gt.STPMAX ) then
!!$      print *,' 01'
      imax= 0
      dmax= 0d0
      do ix=1,3
        d= abs(xi(4)-xi(ix))
        if( dmax.lt.d ) then
          dmax= d
          imax= ix
        endif
      enddo
      !.....eliminate max function point and set point at STPMAX
      imax= 0
      fmax= -1d+30
      do ix=1,3
        if( fmax.lt.fi(ix) ) then
          fmax= fi(ix)
          imax= ix
        endif
      enddo
      do ix=imax+1,3
        xi(ix-1)= xi(ix)
        fi(ix-1)= fi(ix)
        fti(ix-1)= fti(ix)
      enddo
      if( fi(2).gt.fi(1) ) then
        xi(3)= xi(1) +STPMAX
      else
        xi(3)= xi(2) +STPMAX
      endif
      call func(ndim,x0+xi(3)*g,fi(3),fti(3))
      goto 10
    else if( fi(4).gt.fmax ) then ! fi(4) is maximum
!!$      print *,' 02'
      imin= 0
      dmin= 1d+30
      do ix=1,3
        d= abs(xi(4)-xi(ix))
        if( dmin.gt.d ) then
          dmin= d
          xmin= xi(ix)
          imin= ix
        endif
      enddo
      !.....eliminate nearest point to xi(4) and add a new point
      do ix=imin+1,3
        xi(ix-1)= xi(ix)
        fi(ix-1)= fi(ix)
        fti(ix-1)= fti(ix)
      enddo
      xi(3)= (xmin +xi(4))*0.5
      call func(ndim,x0+xi(3)*g,fi(3),fti(3))
      goto 10
    endif

    !.....step 5: check convergence
    if( dmin.lt.xtol ) then
      imin= 0
      fmin= 1d+30
      do ix=1,4
        if( fmin.gt.fi(ix) ) then
          imin= ix
          fmin= fi(ix)
        endif
      enddo
      f= fi(imin)
      a= xi(imin)
      return
    endif

!.....step 6: discard point of highest f value and replace it by xi(4)
!!$    print *,' 03'
    imax= 0
    fmax= -1d+30
    do ix=1,3
      if( fmax.lt.fi(ix) ) then
        imax= ix
        fmax= fi(ix)
      endif
    enddo
    xi(imax)= xi(4)
    fi(imax)= fi(4)
    fti(imax)= fti(4)
    goto 10

  end subroutine quad_interpolate
!=======================================================================
  subroutine golden_section(ndim,x0,g,f,ftst,xtol,gtol,ftol,alpha,iprint &
       ,iflag,myid,func)
    implicit none
    integer,intent(in):: ndim,iprint,myid
    integer,intent(inout):: iflag
    real(8),intent(in):: xtol,gtol,ftol,x0(ndim),g(ndim)
    real(8),intent(inout):: f,alpha,ftst
    interface
      subroutine func(n,x,ftrn,ftst)
        integer,intent(in):: n
        real(8),intent(in):: x(n)
        real(8),intent(out):: ftrn,ftst
      end subroutine func
    end interface

    real(8),parameter:: STP0 = 1d-1
    real(8),parameter:: GR   = 0.61803398875d0
    real(8),parameter:: GR2  = 1d0 -GR
    integer,parameter:: MAXITER= 100

    integer:: iter
    real(8):: a,b1,b2,c,fa,fb1,fb2,fc,xl
    real(8):: ftb1,ftb2,fta,ftc

    a= 0d0
    b1= STP0
    call get_bracket(ndim,x0,g,a,b1,c,fa,fb1,fc,fta,ftb1,ftc,&
         iprint,iflag,myid,func)
    if( iflag/1000.ne.0 ) return
    xl= (c-a)
    b1= a +GR2*xl
    b2= a +GR *xl
    call func(ndim,x0+b1*g,fb1,ftb1)
    call func(ndim,x0+b2*g,fb2,ftb2)

    iter= 0
10  continue
    iter= iter +1
    if( iter.gt.MAXITER ) then
      if( myid.eq.0 ) then
        print *,'[Error] iter.gt.MAXITER in golden_section.'
        print *,'  iter,MAXITER = ',iter,MAXITER
      endif
      iflag= iflag +100
      return
!!$      stop
    endif
!!$    write(6,'(a,2(2x,4es11.3))') ' a,b1,b2,c,fa,fb1,fb2,fc=' &
!!$         ,a,b1,b2,c,fa,fb1,fb2,fc
    if( fb1.gt.fb2 ) then
      a= b1
      fa= fb1
      fta= ftb1
      b1= b2
      fb1= fb2
      ftb1= ftb2
      xl= (c-a)
      b2= a +GR*xl
      call func(ndim,x0+b2*g,fb2,ftb2)
    else
      c= b2
      fc= fb2
      ftc= ftb2
      b2= b1
      fb2= fb1
      ftb2= ftb1
      xl= (c-a)
      b1= a +GR2*xl
      call func(ndim,x0+b1*g,fb1,ftb1)
    endif
!!$    print *,' xl,c,a,xtol=',xl,c,a,xtol
    if( xl.lt.xtol ) then
      if( fb1.gt.fb2 ) then
        alpha= b2
        f= fb2
        ftst= ftb2
      else
        alpha= b1
        f= fb1
        ftst= ftb1
      endif
      return
    endif
    goto 10

  end subroutine golden_section
!=======================================================================
  subroutine armijo_search(ndim,x0,d,f,ftst,g,alpha,iprint &
       ,iflag,myid,func,niter)
!  
!  1D search using Armijo rule.
!    
    implicit none
    integer,intent(in):: ndim,iprint,myid
    integer,intent(inout):: iflag,niter
    real(8),intent(in):: x0(ndim),g(ndim),d(ndim)
    real(8),intent(inout):: f,alpha,ftst
    interface
      subroutine func(n,x,ftrn,ftst)
        integer,intent(in):: n
        real(8),intent(in):: x(n)
        real(8),intent(out):: ftrn,ftst
      end subroutine func
    end interface

!!$  real(8),external:: sprod
    real(8),parameter:: xtiny  = 1d-14
    integer:: iter,i,ig
    real(8):: alphai,xigd,f0,fi,sgnx,pval,pval0,absx,fp,pvalp,alphap,ftsti
    real(8),allocatable,dimension(:):: x1(:),gpena(:)
    logical,save:: l1st = .true.

    if( l1st ) then
      if( myid.eq.0 .and. iprint.gt.0 ) then
        write(6,'(a)') ' Armijo rule parameters:'
        write(6,'(a,es12.4)') '   c       = ',armijo_xi
        write(6,'(a,f10.4)') '   tau     = ',armijo_tau
        write(6,'(a,i5)')   '   maxiter = ',armijo_maxiter
      endif
      l1st = .false.
    endif

    if( .not. allocated(x1)) allocate(x1(ndim),gpena(ndim))
    xigd= sprod(ndim,g,d)*armijo_xi
    if( xigd.gt.0d0 ) then
      iflag= iflag + 100
      if( myid.eq.0 .and. iprint.gt.0 ) print *,' WARNING: g*d > 0.0'
      return
    endif
    alphai= alpha
    pval0= 0d0
    gpena(1:ndim)= 0d0
    if( trim(cpena).eq.'lasso' ) then
      do i=1,ndim
        absx= abs(x0(i))
        sgnx= sign(1d0,x0(i))
        if( absx.gt.xtiny ) gpena(i)= pwgt*sgnx
        pval0= pval0 +pwgt*absx
      enddo
    else if( trim(cpena).eq.'glasso' ) then
      glval(0:ngl)= 0d0
      do i=1,ndim
        ig= iglid(i)
        if( ig.gt.0 ) glval(ig)= glval(ig) +x0(i)*x0(i)
      enddo
      glval(0)= 1d0
      do ig=1,ngl
        glval(ig)= sqrt(glval(ig))
        pval0= pval0 +pwgt*glval(ig)
      enddo
      do i=1,ndim
        ig= iglid(i)
        if( ig.eq.0 ) then ! i is not in a group
          absx= abs(x0(i))
          sgnx= sign(1d0,x0(i))
          if( absx.gt.xtiny ) gpena(i)= pwgt*sgnx
          pval0= pval0 +pwgt*absx
        else if( ig.gt.0 ) then ! i is in a group
          if( glval(ig).gt.xtiny) gpena(i)= pwgt*x0(i)/glval(ig)
        endif
      enddo
    else if( trim(cpena).eq.'ridge' ) then
      do i=1,ndim
        pval0= pval0 +pwgt*x0(i)*x0(i)
        gpena(i)= 2d0*pwgt*x0(i)
      enddo
    endif

    f0= f
    do iter=1,armijo_maxiter
      x1(1:ndim)= x0(1:ndim)
      if( trim(cpena).eq.'lasso' .or.trim(cpena).eq.'glasso') then
        call soft_threshold(ndim,x1,d,alphai)
        call wrap_ranges(ndim,x1)
      else
        x1(1:ndim)= x1(1:ndim) +alphai*d(1:ndim)
        call wrap_ranges(ndim,x1)
      endif
      call func(ndim,x1,fi,ftsti)
!!$      if( myid.eq.0 ) print *,'iter,alphai,fi=',iter,alphai,fi
      pval= 0d0
      if( trim(cpena).eq.'lasso' ) then
        do i=1,ndim
          pval= pval +pwgt*abs(x1(i))
        enddo
      else if( trim(cpena).eq.'glasso' ) then
        pval= 0d0
        glval(0:ngl)= 0d0
        do i=1,ndim
          ig= iglid(i)
          if( ig.gt.0 ) glval(ig)= glval(ig) +x1(i)*x1(i)
          if( ig.eq.0 ) pval= pval +pwgt*abs(x1(i))
        enddo
        do ig=1,ngl
          glval(ig)= sqrt(glval(ig))
          pval= pval +pwgt*glval(ig)
        enddo
      else if( trim(cpena).eq.'ridge' ) then
        do i=1,ndim
          pval= pval +pwgt*x1(i)*x1(i)
        enddo
      endif
      if( myid.eq.0 .and. iprint.ge.10 ) write(6,'(a,i5,3es15.7)') &
           ' armijo: iter,fi+pval-(f0+pval0),xigd*alphai,alphai=',&
           iter,fi+pval-(f0+pval0),xigd*alphai,alphai
      if( fi+pval-(f0+pval0).le.xigd*alphai ) then
        f= fi
        alpha= alphai
        ftst= ftsti
!!$        if(myid.eq.0) print *,'armijo finishes with iter=',iter
        niter = iter
        return
      endif
      fp= fi
      pvalp= pval
      alphap= alphai
      alphai= alphai*armijo_tau
    enddo

    if(myid.eq.0 .and. iprint.gt.0 ) &
         print *,'[Error] iter.gt.MAXITER in armijo_search.'
    iflag= iflag +100
    niter= iter
    if( myid.eq.0 .and. iprint.gt.0 ) then
      write(6,'(a,es13.5)') '  alphai   = ',alphai
      write(6,'(a,es13.5)') '  xigd    = ',xigd
      write(6,'(a,es13.5)') '  norm(g) = ',sqrt(sprod(ndim,g,g))
      if( trim(cpena).eq.'lasso' .or. trim(cpena).eq.'glasso' .or. &
           trim(cpena).eq.'ridge' ) then
        write(6,'(a,es13.5)') '  pval    = ',pval
      endif
    endif
    return

  end subroutine armijo_search
!=======================================================================
  function sprod(n,a,b)
    implicit none
    integer,intent(in):: n
    real(8),intent(in):: a(n),b(n)
    real(8):: sprod

    integer:: i
    sprod= 0d0
    do i=1,n
      sprod= sprod +a(i)*b(i)
    enddo
    return
  end function sprod
!=======================================================================
  subroutine soft_threshold(ndim,x,d,alpha)
!
!  Estimate next weight value using soft threshold
!
    implicit none
    integer,intent(in):: ndim
    real(8),intent(in):: d(ndim),alpha
    real(8),intent(inout):: x(ndim)

    integer:: i
    real(8):: xad,sgn,val,xt

    do i=1,ndim
      xad= x(i) +alpha*d(i)
      sgn=sign(1d0,xad) 
      val= max(abs(xad)-alpha*pwgt,0d0)
      x(i)= sgn*val
    enddo
    return
  end subroutine soft_threshold
!=======================================================================
  subroutine fs(ndim,x,f,g,d,xtol,gtol,ftol,maxiter &
       ,iprint,iflag,myid,func,grad)
!
!  Forward Stagewise (FS) regression
!
    implicit none
    integer,intent(in):: ndim,maxiter,iprint,myid
    integer,intent(inout):: iflag
    real(8),intent(in):: xtol,gtol,ftol
    real(8),intent(inout):: f,x(ndim),g(ndim),d(ndim)
!!$    real(8):: func,grad
    interface
      subroutine func(n,x,ftrn,ftst)
        integer,intent(in):: n
        real(8),intent(in):: x(n)
        real(8),intent(out):: ftrn,ftst
      end subroutine func
      subroutine grad(n,x,gtrn)
        integer,intent(in):: n
        real(8),intent(in):: x(n)
        real(8),intent(out):: gtrn(n)
      end subroutine grad
    end interface

    real(8),parameter:: eps = 1d0
    real(8),parameter:: xtiny= 1d-14
    integer:: iter,i,imax,ig
    real(8):: alpha,gnorm,gmax,absg,sgnx,xad,val,absx,pval,fp,ftst
    real(8),allocatable,dimension(:):: xt,gpena,grpg

    if( trim(cpena).ne.'lasso' .and. trim(cpena).ne.'glasso' ) then
      if(myid.eq.0) then
        print *,'>>> fs works only with lasso or glasso.'
      endif
      iflag= iflag +100
      return
    endif

    if( .not.allocated(xt) ) allocate(xt(ndim) &
         ,gpena(ndim),grpg(0:ngl))

    xt(1:ndim)= x(1:ndim)
!!$    xt(1:ndim)= 1d-6

    do iter=1,maxiter
      fp= f
!.....find maximum contribution in g
      call func(ndim,xt,f,ftst)
      call grad(ndim,xt,g)
      pval= 0d0
      gpena(1:ndim)= 0d0
      if( trim(cpena).eq.'lasso' ) then
        do i=1,ndim
          sgnx= sign(1d0,xt(i))
          absx= abs(xt(i))
          pval= pval +pwgt*absx
          if(absx.gt.xtiny) g(i)= g(i) +pwgt*sgnx
        enddo
      else if( trim(cpena).eq.'glasso' ) then
        glval(0:ngl)= 0d0
        do i=1,ndim
          ig= iglid(i)
          if( ig.gt.0 ) glval(ig)= glval(ig) +xt(i)*xt(i)
        enddo
        glval(0)= 1d0
        do ig=1,ngl
          glval(ig)= sqrt(glval(ig))
          pval= pval +pwgt*glval(ig)
        enddo
        do i=1,ndim
          ig= iglid(i)
          if( ig.eq.0 ) then ! i is not in a group
            absx= abs(xt(i))
            sgnx= sign(1d0,xt(i))
            if( absx.gt.xtiny ) gpena(i)= pwgt*sgnx
            pval= pval +pwgt*absx
          else if( ig.gt.0 ) then ! i is in a group
            if( glval(ig).gt.xtiny) gpena(i)= pwgt*xt(i)/glval(ig)
          endif
        enddo
      endif
      g(1:ndim)= g(1:ndim) +gpena(1:ndim)
      gnorm= sqrt(sprod(ndim,g,g))
      if( myid.eq.0 ) then
!!$        if( iprint.eq.1 .and. mod(iter,ndim).eq.1 ) then
        if( iprint.eq.1 ) then
          write(6,'(a,i8,4es13.5)') ' iter,f,p,gnorm,f-fp=',iter,f &
               ,pval,gnorm,f-fp
          call flush(6)
        else if( iprint.ge.2 ) then
          write(6,'(a,i8,12es13.5)') ' iter,f,p,gnorm,f-fp,x(1:5)=' &
               ,iter,f,pval,gnorm,f-fp,x(1:5)
          call flush(6)
        endif
      endif
      if( gnorm.lt.gtol ) then
        if( myid.eq.0 ) then
          print *,'>>> FS converged wrt gtol'
          write(6,'(a,2es13.5)') '   gnorm,gtol=',gnorm,gtol
        endif
        x(1:ndim)= xt(1:ndim)
        iflag= iflag +2
        return
      endif

      if( trim(cpena).eq.'lasso' ) then
!.....set 0 except the maximum contribution
        imax= 0
        gmax= 0d0
        do i=1,ndim
          absg= abs(g(i))
          if( gmax.lt.absg ) then
            gmax= absg
            imax= i
          endif
        enddo
        alpha= eps
        xad= xt(imax) -alpha*g(imax)
        sgnx= sign(1d0,xad)
        val= max(abs(xad)-alpha*pwgt,0d0)
        xt(imax)= sgnx*val
      else if( trim(cpena).eq.'glasso' ) then
        grpg(0:ngl)= 0d0
        do i=1,ndim
          ig= iglid(i)
          if( ig.gt.0 ) grpg(ig)= grpg(ig) +g(i)*g(i)
        enddo
        imax= 0
        gmax= 0d0
        do ig=1,ngl
          if( gmax.lt.grpg(ig) ) then
            gmax= grpg(ig)
            imax= ig
          endif
        enddo
!!$        print *,'myid,imax,gmax=',myid,imax,gmax
        alpha= eps
        do i=1,ndim
          ig= iglid(i)
          if( ig.eq.0 .or. ig.eq.imax ) then
            xad= xt(i) -alpha*g(i)
            sgnx= sign(1d0,xad)
            val= max(abs(xad)-alpha*pwgt,0d0)
            xt(i)= sgnx*val
          endif
        enddo
      endif
    enddo

    if( myid.eq.0 ) print *,'maxiter exceeded in fs'
    iflag= iflag +10
    x(1:ndim)= xt(1:ndim)
    return
    
  end subroutine fs
!=======================================================================
  subroutine gfs(ndim,x,f,g,d,xtol,gtol,ftol,maxiter &
       ,iprint,iflag,myid,func,grad,cfmethod,niter_eval &
       ,sub_eval,analyze)
!
!  Grouped Forward Stepwise (grouped FS) regression
!
    implicit none
    integer,intent(in):: ndim,maxiter,iprint,myid,niter_eval
    integer,intent(inout):: iflag
    real(8),intent(in):: xtol,gtol,ftol
    real(8),intent(inout):: f,x(ndim),g(ndim),d(ndim)
    character(len=*),intent(in):: cfmethod
!!$    real(8):: func,grad
    interface
      subroutine func(n,x,ftrn,ftst)
        integer,intent(in):: n
        real(8),intent(in):: x(n)
        real(8),intent(out):: ftrn,ftst
      end subroutine func
      subroutine grad(n,x,gtrn)
        integer,intent(in):: n
        real(8),intent(in):: x(n)
        real(8),intent(out):: gtrn(n)
      end subroutine grad
      subroutine sub_eval(iter)
        integer,intent(in):: iter
      end subroutine sub_eval
      subroutine analyze(num)
        integer,intent(in):: num
      end subroutine analyze
    end interface

    integer:: iter,i,imax,ig,itmp,j,igmm,itergfs,niter
    real(8):: alpha,gnorm,gmax,absg,sgnx,xad,val,absx,pval,fp,tmp,gmm,ftst
    real(8),allocatable,save:: xt(:),gmaxgl(:),u(:)
    real(8),save,allocatable:: gg(:,:),v(:),y(:),gp(:) &
         ,ggy(:),ygg(:),aa(:,:),cc(:,:)
    integer:: nmsks,imsk,nftol,nbases,nbasesp
    real(8):: ynorm,svy,svyi,tmp1,tmp2,b

    if( trim(cpena).eq.'lasso' .and. trim(cpena).eq.'glasso' ) then
      if(myid.eq.0) then
        print *,'>>> gfs does not work with lasso or glasso.'
        print *,'>>> so gfs neglects lasso and glasso.'
      endif
    endif

    if( .not.allocated(xt) ) allocate(xt(ndim),u(ndim))
    if( .not.allocated(gg) ) allocate(gg(ndim,ndim) &
         ,v(ndim),y(ndim),gp(ndim),ggy(ndim),ygg(ndim) &
         ,aa(ndim,ndim),cc(ndim,ndim))
    if( .not.allocated(mskgfs) ) then
      allocate(mskgfs(ngl),msktmp(ngl),gmaxgl(ngl))
      mskgfs(1:ngl)= 1
      msktmp(1:ngl)= mskgfs(1:ngl)
    endif
    
    xt(1:ndim)= x(1:ndim)
    do i=1,ndim
      ig= iglid(i)
      if( ig.gt.0 ) xt(i)= 1d-14
    enddo

    nmsks= 0
    do ig=1,ngl
      if( mskgfs(ig).eq.0 ) cycle
      nmsks= nmsks +1
    enddo
    nbasesp= ngl -nmsks

!.....do loop until the conversion criterion is achieved
    iter= 0
    do while(.true.)
!.....First, calc of gradient needs to be done with no masks
!     because it is used to find another new basis
      msktmp(1:ngl)= mskgfs(1:ngl)
      mskgfs(1:ngl)= 0
      call func(ndim,xt,f,ftst)
      call grad(ndim,xt,g)
      mskgfs(1:ngl)= msktmp(1:ngl)
      gnorm= sqrt(sprod(ndim,g,g))
      if( myid.eq.0 ) then
        if( iprint.eq.1 ) then
          write(6,'(a,i8,2es13.5)') ' itergfs,f,gnorm=',itergfs,f,gnorm
        else if( iprint.eq.2 ) then
          write(6,'(a,i8,12es13.5)') ' itergfs,f,gnorm,x(1:5)=' &
               ,itergfs,f,gnorm,xt(1:5)
        endif
        call flush(6)
      endif

      if( nmsks.eq.0 ) then
        if( myid.eq.0 ) then
          print *,'ngl=',ngl
          print *,'nmsks is already 0, and while loop should be finished.'
        endif
        exit
      endif
!.....Find bases with the largest gradient
      gmaxgl(1:ngl)= 0d0
      do i=1,ndim
        ig= iglid(i)
        if( ig.gt.0 ) then
          gmaxgl(ig)= gmaxgl(ig) &
               +g(i)*g(i)
        endif
      enddo
      gmm= 0d0
      igmm= 0
      do ig=1,ngl
!.....Do not take mskgfs==2 into account !
        if( mskgfs(ig).eq.1 .and. gmaxgl(ig).gt.gmm ) then
          gmm= gmaxgl(ig)
          igmm= ig
        endif
      enddo
      if( igmm.eq.0 ) then
        if(myid.eq.0) then
          print *,'igmm.eq.0 !!!'
          print *,'Nothing to do here, and going out from FS.'
        endif
        x(1:ndim)= xt(1:ndim)
        return
      endif
!.....remove mask of bases with large variations
      mskgfs(igmm)= 0
      if(myid.eq.0) print '(a,i5,es12.4,100i2)',' igmm,gmm,mskgfs= ' &
           ,igmm,gmm,mskgfs(1:min(ngl,100))
      nmsks= 0
      do ig=1,ngl
        if( mskgfs(ig).eq.0 ) cycle
        nmsks= nmsks +1
      enddo
      nbases= ngl -nmsks
      if( myid.eq.0 ) print '(a,4i8)','iter,ngl,nmsks,nbases=' &
           ,iter,ngl,nmsks,nbases
      if( nbases > nbasesp .and. nbases.gt.1 ) then
        call analyze(nbasesp)
        nbasesp= nbases
      endif

!.....preparation for BFGS
      gg(1:ndim,1:ndim)= 0d0
      do i=1,ndim
        gg(i,i)= 1d0
      enddo
!.....mask some g that have small contributions
      do i=1,ndim
        ig= iglid(i)
        if( ig.le.0 ) cycle
        if( mskgfs(ig).ne.0 ) g(i)= 0d0
      enddo
      call cap_grad(ndim,g)
      gnorm= sqrt(sprod(ndim,g,g))
      if( myid.eq.0 ) then
        if( iprint.eq.1 ) then
          write(6,'(a,i8,2es13.5)') ' iter,f,gnorm=',iter,f,gnorm
        else if( iprint.eq.2 ) then
          write(6,'(a,i8,12es13.5)') ' iter,f,gnorm,x(1:5)=' &
               ,iter,f,gnorm,x(1:5)
        endif
        call flush(6)
      endif
      nftol= 0
      iflag= 0
      if( mod(iter,niter_eval).eq.0 ) &
           call sub_eval(iter)
!.....BFGS loop begins
      do itergfs=1,nitergfs
        iter= iter +1
        if( iter.gt.maxiter ) exit
        u(1:ndim)= 0d0
        do i=1,ndim
          u(1:ndim)= u(1:ndim) -gg(1:ndim,i)*g(i)
        enddo
!.....mask some u that have small contributions
        do i=1,ndim
          ig= iglid(i)
          if( ig.le.0 ) cycle
          if( mskgfs(ig).ne.0 ) then
            g(i)= 0d0
            u(i)= 0d0
          endif
        enddo
        fp= f
        gp(1:ndim)= g(1:ndim)
        if( mod(iter,niter_eval).eq.0 ) &
             call sub_eval(iter)
!.....line minimization
        if( trim(clinmin).eq.'quadratic' ) then
          call quad_interpolate(ndim,xt,u,f,ftst,xtol,gtol,ftol,alpha &
               ,iprint,iflag,myid,func)
!.....if quad interpolation failed, perform golden section
          if( iflag/100.ne.0 ) then
            iflag= iflag -(iflag/100)*100
            if(myid.eq.0) then
              print *,'Since quad_interpolate failed, call golden_section.'
            endif
            call golden_section(ndim,xt,u,f,ftst,xtol,gtol,ftol,alpha &
                 ,iprint,iflag,myid,func)
          endif
        else if ( trim(clinmin).eq.'golden') then
          call golden_section(ndim,xt,u,f,ftst,xtol,gtol,ftol,alpha &
               ,iprint,iflag,myid,func)
        else ! armijo (default)
          alpha= 1d0
          call armijo_search(ndim,xt,u,f,ftst,g,alpha,iprint &
               ,iflag,myid,func,niter)
!.....if something wrong with armijo search, try opposite direction
          if( iflag/100.ne.0 ) then
!.....Reset iflag when we do opposite direction search.
!            iflag= iflag -(iflag/100)*100
            alpha= -1d0
            if(myid.eq.0) print *,'trying opposite direction...'
            call armijo_search(ndim,xt,u,f,ftst,g,alpha,iprint &
                 ,iflag,myid,func,niter)
          endif
        endif
!.....get out of bfgs loop
        if( iflag/100.ne.0 ) then
          if( itergfs.eq.1 ) then
!!$            if( myid.eq.0 ) then
!!$              print *,'Armijo failed at 1st step of FS, so going out of FS.'
!!$            endif
!!$            x(1:ndim)= xt(1:ndim)
!!$            return
            if( myid.eq.0 ) then
              print *,"Armijo failed at 1st step of FS,"&
                   //" but keep going forward..."
            endif
!!$!.....Set mask as 2, which means this basis will be not included
!!$!     and not taken into consideration anymore.
!!$            mskgfs(igmm)= 2
          endif
          exit
        endif
        xt(1:ndim)= xt(1:ndim) +alpha*u(1:ndim)
        call grad(ndim,xt,g)
        do i=1,ndim
          ig= iglid(i)
          if( ig.le.0 ) cycle
          if( mskgfs(ig).ne.0 ) then
            g(i)= 0d0
            u(i)= 0d0
          endif
        enddo
        call cap_grad(ndim,g)
        gnorm= sqrt(sprod(ndim,g,g))
        if( myid.eq.0 ) then
          if( iprint.eq.1 ) then
            write(6,'(a,i8,2es13.5)') ' itergfs,f,gnorm=',itergfs,f,gnorm
          else if( iprint.eq.2 ) then
            write(6,'(a,i8,12es13.5)') ' itergfs,f,gnorm,x(1:5)=' &
                 ,itergfs,f,gnorm,xt(1:5)
          endif
          call flush(6)
        endif

!.....check convergence
        if( gnorm.lt.gtol ) then
          if( myid.eq.0 ) then
            print *,'>>> QM in gFS converged wrt gtol'
            write(6,'(a,2es13.5)') '   gnorm,gtol=',gnorm,gtol
          endif
          x(1:ndim)= xt(1:ndim)
          iflag= iflag +2
          exit
        else if( abs(f-fp)/abs(fp).lt.ftol) then
          nftol= nftol +1
          if( nftol.gt.numtol ) then
            if( myid.eq.0 ) then
              print *,'>>> gFS may be converged because of ftol ' // &
                   'over 10 times.'
            endif
            x(1:ndim)= xt(1:ndim)
            iflag= iflag +3
            exit
          else
            if( myid.eq.0 ) then
              print *,'>>> gg initialized because |f-fp|/|fp|<ftol '
            endif
            gg(1:ndim,1:ndim)= 0d0
            do i=1,ndim
              gg(i,i)= 1d0
            enddo
            cycle
          endif
        endif

        v(1:ndim)= alpha *u(1:ndim)
        y(1:ndim)= g(1:ndim) -gp(1:ndim)
        ynorm= sprod(ndim,y,y)
        if( ynorm.lt.1d-14 ) then
          if(myid.eq.0) then
            print *,'>>> gg initialized because y*y < 1d-14'
            print *,'  ynorm=',ynorm
            print *,'  alpha=',alpha
          endif
          gg(1:ndim,1:ndim)= 0d0
          do i=1,ndim
            gg(i,i)= 1d0
          enddo
          cycle
        endif

!.....update G matrix, gg, according to BFGS
        svy= sprod(ndim,v,y)
        svyi= 1d0/svy
        do i=1,ndim
          tmp1= 0d0
          tmp2= 0d0
          do j=1,ndim
            aa(j,i)= v(j)*v(i) *svyi
            tmp1= tmp1 +gg(i,j)*y(j)
            tmp2= tmp2 +y(j)*gg(j,i)
          enddo
          ggy(i)= tmp1
          ygg(i)= tmp2
        enddo
        cc(1:ndim,1:ndim)= 0d0
        do j=1,ndim
          do i=1,ndim
            cc(i,j)=cc(i,j) +(v(i)*ygg(j) +ggy(i)*v(j)) *svyi
          enddo
        enddo
        b= 1d0
        do i=1,ndim
          b=b +y(i)*ggy(i) *svyi
        enddo
        aa(1:ndim,1:ndim)= aa(1:ndim,1:ndim) *b
        gg(1:ndim,1:ndim)=gg(1:ndim,1:ndim) +aa(1:ndim,1:ndim) &
             -cc(1:ndim,1:ndim)
      enddo
      if( iter.gt.maxiter ) exit
    enddo

999 if( myid.eq.0 ) print *,'maxiter exceeded in gfs'
    iflag= iflag +10
    x(1:ndim)= xt(1:ndim)
    return

  end subroutine gfs
!=======================================================================
  subroutine cap_grad(ndim,g)
!
!  Set the ceiling of gradient to avoid too large gnorm.
!
    implicit none
    integer,intent(in):: ndim
    real(8),intent(inout):: g(ndim)
    real(8),parameter:: gmax= 1.0d0
    real(8):: gnorm
    
    gnorm= sqrt(sprod(ndim,g,g))

    if( gnorm.gt.gmax ) then
      g(1:ndim)= g(1:ndim) *gmax /gnorm
    endif
    return
  end subroutine cap_grad
!=======================================================================
  subroutine sa(ndim,xbest,fbest,xranges,xtol,gtol,ftol,maxiter &
       ,iprint,iflag,myid,func,cfmethod,niter_eval,sub_eval)
!
! Simulated Annealing
!
    use random
    implicit none
    integer,intent(in):: ndim,iprint,myid,maxiter,niter_eval
    integer,intent(inout):: iflag
    real(8),intent(in):: xtol,gtol,ftol,xranges(2,ndim)
    real(8),intent(inout):: fbest,xbest(ndim)
    character(len=*),intent(in):: cfmethod
    interface
      subroutine func(n,x,ftrn,ftst)
        integer,intent(in):: n
        real(8),intent(in):: x(n)
        real(8),intent(out):: ftrn,ftst
      end subroutine func
      subroutine sub_eval(iter)
        integer,intent(in):: iter
      end subroutine sub_eval
    end interface

    integer:: iter,idim,nadpt,i,l
    real(8):: f,ft,temp,xw,dx,p,pt,ptrans,ftst,tau,xmin,xmax
    real(8),allocatable:: x(:),xt(:)
    logical,save:: l1st = .true.

    if( l1st ) then
      if(.not.allocated(sa_xws) ) allocate(sa_xws(ndim))
      do i=1,ndim
!!$        sa_xws(i) = max(xbest(i)*sa_xw0,1d-2)
        sa_xws(i) = 1d-2
      enddo
      tau = max(sa_tau,1d0)
      if(myid.eq.0 .and. iprint.ne.0) then
        write(6,*) ''
        write(6,'(a)') '------------------------------------------------------------------------'
        write(6,'(a)') '                    Simulated annealing'
        write(6,'(a)') '------------------------------------------------------------------------'
        write(6,'(a,a)') ' Temperature control method = ',trim(sa_tctrl)
        if( sa_tctrl(1:3).eq.'exp' ) then
          write(6,'(a,f10.4)') ' Initial temperature = ',sa_temp0
          write(6,'(a,f10.4)') ' Relaxation iteration (tau) = ',tau
        else if( sa_tctrl(1:3).eq.'bes' ) then
          write(6,'(a,f8.1)') ' Division of fbest = ',sa_div_best
        else
          write(6,'(a,f10.4)') ' Initial temperature = ',sa_temp0
        endif
        print *,''
      endif
      
    endif

    if( .not.allocated(x) ) allocate(x(ndim),xt(ndim))

!.....Initialize
    x(1:ndim)= xbest(1:ndim)
    do idim=1,ndim
      if( x(idim).lt.xranges(1,idim) ) x(idim) = xranges(1,idim)
      if( x(idim).gt.xranges(2,idim) ) x(idim) = xranges(2,idim)
    enddo
    call func(ndim,x,f,ftst)
!!$    p = exp(-f/temp)
    if( f*0d0 .ne. 0d0 ) f = 1d+10
    fbest= f
    if( sa_tctrl(1:3).eq.'bes' ) then
      temp = fbest/sa_div_best
    else
      temp= max(sa_temp0,1d-8)
    endif
    xw= sa_xw0
    nadpt= 0

    iter= 0
    ft = 0d0
    idim = 0
    if( myid.eq.0 .and. iprint.ne.0 ) then
      write(6,'(a,2i10,4es13.5,2f9.5)')&
           ' iter,idim,temp,f,ft,fbest,ptrans,radpt='&
           ,iter,idim,temp,f,ft,fbest,ptrans,0d0
    endif
    
    call sub_eval(0)
!.....Main loop of random displacements
    do iter=1,maxiter

!.....Choose a parameter to be displaced
      idim= urnd()*ndim +1

!.....Compute the displacement using a uniform random number
      xt(1:ndim)= x(1:ndim)
      dx= (urnd()-0.5d0)*sa_xws(idim)
      xt(idim)= xt(idim) +dx
      if( xt(idim).lt.xranges(1,idim) ) xt(idim) = xranges(1,idim)
      if( xt(idim).gt.xranges(2,idim) ) xt(idim) = xranges(2,idim)

!.....Compute function value
      call func(ndim,xt,ft,ftst)
!.....Detect NaN and skip this trial
      if( ft*0d0 .ne. 0d0 ) then
        if( myid.eq.0 .and. iprint.ne.0 ) then
          write(6,'(a,2i10,es12.4,3es13.5,2f9.5)')&
               ' [ft.eq.NaN] iter,idim,sa_xws(idim)=' &
               ,iter,idim,sa_xws(idim)
        endif
!.....Decrease the width of deviation
        sa_xws(idim) = sa_xws(idim) *sa_fctr
        goto 10
      endif
!!$      pt = exp(-ft/temp)

!.....Compute probability of taking the displacement
      ptrans= min(1d0,exp(-(ft-f)/temp))
!!$      ptrans = min(1d0,pt/p)

!.....Store the best one
      if( ft.lt.fbest ) then
        fbest= ft
        xbest(1:ndim)= xt(1:ndim)
        if( sa_tctrl(1:3).eq.'bes' ) temp = fbest/sa_div_best
      endif

      if( mod(iter,niter_eval).eq.0 ) then
        call sub_eval(iter)
      endif

      if( myid.eq.0 .and. iprint.ne.0 ) then
        write(6,'(a,2i10,4es13.5,2f9.5)')&
             ' iter,idim,temp,f,ft,fbest,ptrans,radpt='&
             ,iter,idim,temp,f,ft,fbest,ptrans,dble(nadpt)/iter
        flush(6)
      endif
      
!.....Update the parameter if needed
      if( urnd().lt.ptrans ) then
        x(idim)= xt(idim)
        f= ft
        nadpt= nadpt +1
!.....Increase the width of deviation
        sa_xws(idim) = sa_xws(idim) /sa_fctr
      else
!.....Decrease the width of deviation
        sa_xws(idim) = sa_xws(idim) *sa_fctr
      endif

10    continue
      if( sa_tctrl(1:3).eq.'lin' ) then
!.....Update temperature (linear)
        temp= dble(maxiter-iter)/maxiter *sa_temp0
      else if( sa_tctrl(1:3).eq.'exp' ) then
!.....Update temperature (exponential)
        temp= sa_temp0 *exp(-dble(iter)/tau)
      endif
    enddo

    if( myid.eq.0 ) then
      write(6,'(a,i10,a,i10)') ' Num of adoption in SA='&
           ,nadpt,'/',maxiter
!!$      print *,'sa_xws:'
!!$      do i=1,ndim
!!$        write(6,'(i6,es15.7)') i,sa_xws(i)
!!$      enddo
    endif

!.....Finally compute the function value of the best candidate
    call func(ndim,xbest,f,ftst)

    l1st = .false.
    
  end subroutine sa
!=======================================================================
  subroutine penalty(cpena,pwgt,ndim,f,g,fp,gp,x)
!
! Calculate penalty term and its derivative.
! lasso and ridge are available.
!
    implicit none
    character(len=*),intent(in):: cpena
    integer,intent(in):: ndim
    real(8),intent(in):: pwgt,f,g(ndim),x(ndim)
    real(8),intent(out):: fp,gp(ndim)

    integer:: i,ig
    real(8):: absx,sgnx
    real(8),parameter:: xtiny  = 1d-14

    fp= 0d0
    gp(1:ndim)= 0d0
    if( trim(cpena).eq.'lasso' ) then
      do i=1,ndim
        absx= abs(x(i))
        fp= fp +pwgt*absx
        sgnx= sign(1d0,x(i))
        if( absx.gt.xtiny ) gp(i)= pwgt*sgnx
      enddo
    else if( trim(cpena).eq.'glasso' ) then
      glval(0:ngl)= 0d0
      do i=1,ndim
        ig= iglid(i)
        if(ig.gt.0) glval(ig)= glval(ig) +x(i)*x(i)
      enddo
      glval(0)= 1d0
      do ig=1,ngl
        glval(ig)= sqrt(glval(ig))
        fp= fp +pwgt*glval(ig)
      enddo
      do i=1,ndim
        ig= iglid(i)
        if( ig.eq.0 ) then ! i is not in a group
          absx= abs(x(i))
          sgnx= sign(1d0,x(i))
          if( absx.gt.xtiny ) gp(i)= pwgt*sgnx
          fp= fp +pwgt*absx
        else if( ig.gt.0 ) then ! i is in a group
          if( glval(ig).gt.xtiny) gp(i)= pwgt*x(i)/glval(ig)
        endif
      enddo
    else if( trim(cpena).eq.'ridge' ) then
      do i=1,ndim
        fp= fp +pwgt*x(i)*x(i)
        gp(i)= 2d0*pwgt*x(i)
      enddo
    endif

    return
  end subroutine penalty
!=======================================================================
  subroutine random_search(ndim,xbest,fbest,xranges,xtol,gtol,ftol,maxiter &
       ,iprint,iflag,myid,func,cfmethod,niter_eval,sub_eval)
!
!  Pure random search with variable range
!
    use random
    implicit none
    integer,intent(in):: ndim,iprint,myid,maxiter,niter_eval
    integer,intent(inout):: iflag
    real(8),intent(in):: xtol,gtol,ftol,xranges(2,ndim)
    real(8),intent(inout):: fbest,xbest(ndim)
    character(len=*),intent(in):: cfmethod
    interface
      subroutine func(n,x,ftrn,ftst)
        integer,intent(in):: n
        real(8),intent(in):: x(n)
        real(8),intent(out):: ftrn,ftst
      end subroutine func
      subroutine sub_eval(iter)
        integer,intent(in):: iter
      end subroutine sub_eval
    end interface

    integer:: iter,idim,i
    real(8):: f,ftst,fmin,xmin,xmax,xi
    real(8),allocatable:: x(:)
    logical,save:: l1st = .true.

    if( l1st ) then
      if( .not. allocated(x) ) allocate(x(ndim))
      l1st = .false.
    endif

    x(1:ndim)= xbest(1:ndim)
    call func(ndim,x,f,ftst)
    if( fmin*0d0.ne.0d0 ) then  ! NaN
      fmin = 1.0d+10
    endif
    fmin = f
    if( myid.eq.0 .and. iprint.gt.0 ) &
         write(6,'(a,i8,2es15.7,100f7.3)') ' iter,f,fmin = ',0,f,fmin
    
    do iter=1,maxiter
      do idim=1,ndim
        xmin = xranges(1,idim)
        xmax = xranges(2,idim)
        xi = (xmax-xmin)*urnd() +xmin
        x(idim) = xi
      enddo

      call func(ndim,x,f,ftst)
      if( f*0d0.ne.0d0 ) then  !NaN
        f = 1.0d+10
      endif
      if( myid.eq.0 .and. iprint.gt.0 ) &
           write(6,'(a,i8,2es15.7,100f7.3)') ' iter,f,fmin = ',iter,f,fmin
      if( f.lt.fmin ) then
        if( myid.eq.0 .and. iprint.gt.0 ) then
          write(6,'(a,i8,2es15.7)') ' fmin is updated: iter,fmin,df= ' &
               ,iter,f,abs(f-fmin)
        endif
        fmin = f
        fbest = fmin
        xbest(1:ndim) = x(1:ndim)
        call sub_eval(iter)
      endif
    enddo

  end subroutine random_search
!=======================================================================
  subroutine metadynamics(ndim,xbest,fbest,xranges,xtol,gtol,ftol,maxiter &
       ,iprint,iflag,myid,func,cfmethod,niter_eval,sub_eval)
!
!  Metadynamics for minimum search
!  Use simulated annealing-like search for local minimum search
!  and once the local minimum is found, add a gaussian potential to
!  the minimum to make it possible to escape from the minimum.
!
    use random
    implicit none
    integer,intent(in):: ndim,iprint,myid,maxiter,niter_eval
    integer,intent(inout):: iflag
    real(8),intent(in):: xtol,gtol,ftol,xranges(2,ndim)
    real(8),intent(inout):: fbest,xbest(ndim)
    character(len=*),intent(in):: cfmethod
    interface
      subroutine func(n,x,ftrn,ftst)
        integer,intent(in):: n
        real(8),intent(in):: x(n)
        real(8),intent(out):: ftrn,ftst
      end subroutine func
      subroutine sub_eval(iter)
        integer,intent(in):: iter
      end subroutine sub_eval
    end interface

    integer:: iter,idim,nadpt,i,ng,ig,interval
    real(8):: f,ft,temp,xw,dx,p,pt,ptrans,ftst,tau,xmin,xmax,adx,&
         fg,fgt,tmp,gval
    real(8),allocatable:: x(:),xt(:),xxt(:)
    logical,save:: l1st = .true.

    if( l1st ) then
      interval = max(maxiter/md_ng,1)
      if(.not.allocated(md_gp) ) allocate(md_gp(ndim,md_ng))
      md_gp(:,:) = 0d0
      if(.not.allocated(sa_xws) ) allocate(sa_xws(ndim))
      do i=1,ndim
!!$        sa_xws(i) = max(xbest(i)*sa_xw0,1d-2)
        sa_xws(i) = 1d-2
      enddo
      tau = max(sa_tau,1d0)
      if(myid.eq.0 .and. iprint.ne.0) then
        print *,''
        print *,'******************** Metadynamics ********************'
        print *,''
        write(6,'(a,i5)') ' Number of gaussians to be added = ',md_ng
        write(6,'(a,f10.4)') ' Gaussian height         = ',md_height
        write(6,'(a,f10.4)') ' Gaussian width (sigma)  = ',md_sigma
      endif
      l1st = .false.
    endif

    if( .not.allocated(x) ) allocate(x(ndim),xt(ndim),xxt(ndim))

!.....Initialize
    x(1:ndim)= xbest(1:ndim)
    do idim=1,ndim
      if( x(idim).lt.xranges(1,idim) ) x(idim) = xranges(1,idim)
      if( x(idim).gt.xranges(2,idim) ) x(idim) = xranges(2,idim)
    enddo
    call func(ndim,x,f,ftst)
!!$    p = exp(-f/temp)
    if( f*0d0.ne.0d0 ) f = 1.0d+10
    fbest= f
    fg = f
!!$    temp= sa_temp0
    temp= fbest/sa_div_best
    xw= sa_xw0
    nadpt= 0
    ng = 0

    call sub_eval(0)
!.....Main loop of random displacements
    do iter=1,maxiter

!.....Choose a parameter to be displaced
      idim= urnd()*ndim +1

!.....Compute the displacement using a uniform random number
      xt(1:ndim)= x(1:ndim)
      dx= (urnd()-0.5d0)*sa_xws(idim)
      xt(idim)= xt(idim) +dx
      if( xt(idim).lt.xranges(1,idim) ) xt(idim) = xranges(1,idim)
      if( xt(idim).gt.xranges(2,idim) ) xt(idim) = xranges(2,idim)

!.....Compute function value
      call func(ndim,xt,ft,ftst)
!.....Detect NaN and skip this trial
      if( ft*0d0 .ne. 0d0 ) then
        if( myid.eq.0 .and. iprint.ne.0 ) then
          write(6,'(a,2i10,es12.4,3es13.5,2f9.5)')&
               ' [ft.eq.NaN] iter,idim,sa_xws(idim)=' &
               ,iter,idim,sa_xws(idim)
        endif
!.....Decrease the width of deviation
        sa_xws(idim) = sa_xws(idim) *sa_fctr
        goto 10
      endif
!.....Add gaussian potentials to ft
      fgt = ft
      do ig=1,ng
        xxt(1:ndim) = xt(1:ndim)-md_gp(1:ndim,ig)
        adx = dot_product(xxt(1:ndim),xxt(1:ndim))
        tmp = md_height*exp(-adx/2/md_sigma**2)
        fgt = fgt +tmp
      enddo
      fg = f
      do ig=1,ng
        xxt(1:ndim) = x(1:ndim)-md_gp(1:ndim,ig)
        adx = dot_product(xxt(1:ndim),xxt(1:ndim))
        tmp = md_height*exp(-adx/2/md_sigma**2)
        fg = fg +tmp
      enddo
!!$      pt = exp(-ft/temp)

!.....Store the best one
      if( ft.lt.fbest ) then
        fbest= ft
        xbest(1:ndim)= xt(1:ndim)
        temp = fbest/sa_div_best
      endif

!.....Compute probability of taking the displacement
      ptrans= min(1d0,exp(-(fgt-fg)/temp))
!!$      ptrans = min(1d0,pt/p)

      if( mod(iter,niter_eval).eq.0 ) then
        call sub_eval(iter)
      endif

      if( myid.eq.0 .and. iprint.ne.0 ) then
        write(6,'(a,2i10,6es13.5,2f9.5,i5,es13.5)')&
             ' iter,idim,temp,f,ft,fbest,fg,fgt,ptrans,radpt,ng='&
             ,iter,idim,temp,f,ft,fbest,fg,fgt,ptrans,dble(nadpt)/iter,ng
      endif
      
!.....Update the parameter if adopted
      if( urnd().lt.ptrans ) then
        x(idim)= xt(idim)
        f= ft
        fg= fgt
        nadpt= nadpt +1
!.....Increase the width of deviation
        sa_xws(idim) = sa_xws(idim) /sa_fctr
      else
!.....Decrease the width of deviation
        sa_xws(idim) = sa_xws(idim) *sa_fctr
      endif

!.....Put a gaussian potential to the current variable position
      if( mod(iter,interval).eq.0 ) then
        ng = ng + 1
        if( ng.gt.md_ng ) then
          if( myid.eq.0 .and. iprint.gt.0 ) then
            write(6,'(a,i8)') ' Number of gaussian exceeds md_ng = ',md_ng
          endif
          exit
        endif
        md_gp(1:ndim,ng) = x(1:ndim)
      endif

10    continue
!!$!.....Update temperature (linear)
!!$      temp= dble(maxiter-iter)/maxiter *sa_temp0
!!$!.....Update temperature (exponential)
!!$      temp= sa_temp0 *exp(-dble(iter)/tau)

    enddo

    if( myid.eq.0 ) then
      write(6,'(a,i10,a,i10)') ' Num of adoption in SA in Metadynamics='&
           ,nadpt,'/',maxiter
    endif

!.....Finally compute the function value of the best candidate
    call func(ndim,xbest,f,ftst)

    
  end subroutine metadynamics
!=======================================================================
  subroutine ga(ndim,xbest,fbest,xranges,xtol,gtol,ftol,maxiter &
       ,iprint,iflag,myid,func,cfmethod,niter_eval,sub_eval,sub_ergrel)
!
! Genetic algorithm (GA) which does not use gradient information.
! GA itself is a serial code, but the function evaluation can be parallel.
!
    use random
    implicit none
    integer,intent(in):: ndim,iprint,myid,maxiter,niter_eval
    integer,intent(inout):: iflag
    real(8),intent(in):: xtol,gtol,ftol,xranges(2,ndim)
    real(8),intent(inout):: fbest,xbest(ndim)
    character(len=*),intent(in):: cfmethod
    interface
      subroutine func(n,x,ftrn,ftst)
        integer,intent(in):: n
        real(8),intent(in):: x(n)
        real(8),intent(out):: ftrn,ftst
      end subroutine func
      subroutine sub_eval(iter)
        integer,intent(in):: iter
      end subroutine sub_eval
      subroutine sub_ergrel(cadd)
        character(len=*),intent(in):: cadd
      end subroutine sub_ergrel
    end interface

    integer:: i,j,k,l,m,n,iter,i1,i2
    real(8):: ftrn,ftst
    logical,save:: l1st = .true.
    integer:: iid,iidbest
    type(individual),allocatable:: indivs(:),offsprings(:)
    real(8),allocatable:: xtmp(:)
    character(len=128):: cadd

    integer,parameter:: io_indivs = 30
    character(len=128),parameter:: cf_indivs = 'out.ga.individuals'
    integer,parameter:: io_steps = 31
    character(len=128),parameter:: cf_steps = 'out.ga.generations'

    if( l1st ) then
!.....Initialize
      ga_ngenes = ndim
      if( ga_noffsp.le.0 ) then
        ga_noffsp = ga_nindivs
      endif
!.....Allocate necessary memory spaces
      allocate(indivs(ga_nindivs),offsprings(ga_noffsp))
      allocate(xtmp(ndim))
      do i=1,ga_noffsp
        allocate(offsprings(i)%genes(ga_ngenes))
        do j=1,ga_ngenes
          allocate(offsprings(i)%genes(j)%bits(ga_nbits))
        enddo
      enddo
      do i=1,ga_nindivs
        allocate(indivs(i)%genes(ga_ngenes))
        do j=1,ga_ngenes
          allocate(indivs(i)%genes(j)%bits(ga_nbits))
        enddo
      enddo
      if( myid.eq.0 .and. iprint.ne.0 ) then
        write(6,*) ''
        write(6,'(a)') '------------------------------------------------------------------------'
        write(6,'(a)') '                    Genetic Algorithm'
        write(6,'(a)') '------------------------------------------------------------------------'
        print '(a,i4)',' Number of individuals = ',ga_nindivs
        print '(a,i4)',' Number of genes       = ',ga_ngenes
        print '(a,i4)',' Number of bits        = ',ga_nbits
        print '(a,i4)',' Number of offsprings  = ',ga_noffsp
        print *,''
      endif
      l1st = .false.
    endif

    if( myid.eq.0 ) then
      open(io_indivs,file=cf_indivs,status='replace')
      write(io_indivs,'(a)') '# iid, fval, vars...'
      open(io_steps,file=cf_steps,status='replace')
      write(io_steps,'(a)') '# iter, iid, ftrn'
    endif
100 format(i6,es14.6,100f9.4)

    iter = 0
    
!.....Create population that includes some individuals
    iid = 0
    fbest = 1d+30
    do i=1,ga_nindivs
      do j=1,ga_ngenes
        indivs(i)%genes(j)%vmin = xranges(1,j)
        indivs(i)%genes(j)%vmax = xranges(2,j)
      enddo
      do j=1,ga_ngenes
        call var2gene(xbest(j),indivs(i)%genes(j))
      enddo
!.....Use a given X in case of i==1, otherwise mutate it with high rate.
      if( i.ne.1 ) then
        do j=1,ga_ngenes
          call mutate(indivs(i)%genes(j),0.25d0)
        enddo
      endif
      call indiv2vars(indivs(i),ndim,xtmp)
      call func(ndim,xtmp,ftrn,ftst)
      iid = iid + 1
!.....Detect NaN and replace it with 1d+10
      if( ftrn*0d0 .ne. 0d0 ) then
        if( myid.eq.0 .and. iprint.ne.0 ) then
          write(6,'(a,2i10)')&
               ' [ftrn.eq.NaN] iter,iid = ',iter,iid
        endif
        ftrn = fupper_lim
      else if( ftrn.gt.fupper_lim ) then
        ftrn = fupper_lim
      endif
      indivs(i)%fvalue = ftrn
      indivs(i)%fitness = 1d0/ftrn
      indivs(i)%iid = iid
      if( myid.eq.0 ) write(io_indivs,100) iid,ftrn,xtmp(1:min(ndim,100))
      if( ftrn.lt.fbest ) then
        fbest = ftrn
        iidbest = iid
        xbest(1:ndim) = xtmp(1:ndim)
      endif
      if( i.eq.1 ) call sub_eval(iter)
    enddo

    if( myid.eq.0 ) then
      write(6,'(a,i8,1x,100es12.4)') &
           " iter,fbest,fvals= ",&
           iter,fbest,(indivs(i)%fvalue,i=1,min(ndim,10))
      do i=1,ga_nindivs
        write(io_steps,'(2i8,es15.7)') iter, indivs(i)%iid, indivs(i)%fvalue
      enddo
    endif

!.....GA loop starts....................................................
    do iter=1,maxiter

!.....Give birth some offsprings by crossover
!!$      print *,'giving birth...'
      do i=1,ga_noffsp
        i1 = int(urnd()*ga_nindivs) +1
10      i2 = int(urnd()*ga_nindivs) +1
        if( i1.eq.i2 ) goto 10
        call crossover(indivs(i1),indivs(i2),offsprings(i))
!.....Mutation of new-born babies
        do j=1,ga_ngenes
          call mutate(offsprings(i)%genes(j),ga_rate_mutate)
        enddo
!.....Evaluate the value of each new-born babies
        call indiv2vars(offsprings(i),ndim,xtmp)
        call func(ndim,xtmp,ftrn,ftst)
        iid = iid + 1
!.....Detect NaN and replace it with 1d+10
        if( ftrn*0d0 .ne. 0d0 ) then
          if( myid.eq.0 .and. iprint.ne.0 ) then
            write(6,'(a,2i10)')&
                 ' [ftrn.eq.NaN] iter,iid = ',iter,iid
          endif
          ftrn = fupper_lim
        endif
        if( ftrn.gt.fupper_lim ) then
          ftrn = fupper_lim
        endif
        offsprings(i)%fvalue = ftrn
        offsprings(i)%fitness = 1d0/ftrn
        offsprings(i)%iid = iid
        if( myid.eq.0 ) write(io_indivs,100) iid,ftrn,xtmp(1:min(ndim,100))
        if( ftrn.lt.fbest ) then
          fbest = ftrn
          iidbest = iid
          xbest(1:ndim) = xtmp(1:ndim)
          if( iprint.ge.2 ) then
            write(cadd,'(i0)') iid
            call sub_ergrel(cadd)
          endif
          call sub_eval(iid)
        endif
      enddo

!.....Selection
!!$      print *,'selecting...'
      call roulette_selection(ga_nindivs,indivs,ga_noffsp,offsprings,fbest)

      if( myid.eq.0 ) then
        write(6,'(a,i8,1x,100es12.4)') &
           " iter,fbest,fvals= ",&
           iter,fbest,(indivs(i)%fvalue,i=1,min(ndim,10))
        do i=1,ga_nindivs
          write(io_steps,'(2i8,es15.7)') iter, indivs(i)%iid, indivs(i)%fvalue
        enddo
        flush(io_indivs)
        flush(io_steps)
      endif
    end do
!.....END of GA loop....................................................

!.....Output information of the best
    if( myid.eq.0 ) then
      write(6,*) ''
      write(6,'(a)') ' The best one in this GA simulation run:'
      write(6,'(a,i0,2x,f0.4)') '   ID, f-value: ',iidbest,fbest
!!$      write(6,'(a,100f7.3)')  '   Variables: ',xbest(1:ndim)
      write(6,*) ''
    endif

!.....Just for outputing out.erg.fin.1 of the best one     
    call func(ndim,xbest,ftrn,ftst)
!!$    if( myid.eq.0 ) then
!!$      write(6,*) 'best one re-caluclated = ', ftrn
!!$    endif

    close(io_indivs)
    close(io_steps)
    return
  end subroutine ga
!=======================================================================
  subroutine bin2dec(nbit,bin,dec)
    integer,intent(in):: nbit
    integer(2),intent(in):: bin(nbit)
    integer,intent(out):: dec
    integer:: i

    dec = 0
    do i=1,nbit
      dec = dec +bin(i)*2**(i-1)
    end do
    return
  end subroutine bin2dec
!=======================================================================
  subroutine dec2bin(nbit,dec,bin)
    integer,intent(in):: nbit
    integer,intent(in):: dec
    integer(2),intent(out):: bin(nbit)
    integer:: i,idec

    idec = dec
    bin(1:nbit) = 0
    do i=1,nbit
      bin(i) = mod(idec,2)
      idec = idec /2
    enddo
    return
  end subroutine dec2bin
!=======================================================================
  subroutine make_pairs(num,pairs)
!
!  Make random pairs from NUM elements.
!
    use random
    integer,intent(in):: num
    integer,intent(out):: pairs(2,num/2)

    integer:: i,j,k,l,m,n,ival,jval
    integer:: chosen(num),navail

    chosen(1:num) = 0
    
    navail = num
    do n=1,num/2
      i = int(urnd()*navail) +1
10    j = int(urnd()*navail) +1
      if( j.eq.i ) goto 10
      l=0
      ival = 0
      jval = 0
      do m=1,num
        if( chosen(m).eq.0 ) then
          l=l+1
          if( l.eq.i ) then
            ival = m
            chosen(m) = 1
          else if( l.eq.j ) then
            jval = m
            chosen(m) = 1
          endif
          if( ival.ne.0 .and. jval.ne.0 ) then
            exit
          endif
        endif
      enddo
      pairs(1,n) = ival
      pairs(2,n) = jval
      navail = 0
      do m=1,num
        if( chosen(m).eq.0 ) then
          navail = navail + 1
        endif
      enddo
    enddo
    return
  end subroutine make_pairs
!=======================================================================
  subroutine mutate(g,rate)
    use random
    type(gene),intent(inout):: g
    real(8),intent(in):: rate
    integer:: i

    do i=1,ga_nbits
      if( urnd().lt.rate ) then
        g%bits(i) = mod(g%bits(i)+1,2)
      endif
    enddo
    return
  end subroutine mutate
!=======================================================================
  subroutine gene2var(g,v)
    type(gene),intent(in):: g
    real(8),intent(out):: v
    integer:: i,dec

    call bin2dec(ga_nbits,g%bits,dec)
    v = g%vmin +dble(dec)*(g%vmax-g%vmin)/(2**ga_nbits-1)
    return
  end subroutine gene2var
!=======================================================================
  subroutine var2gene(v,g)
    real(8),intent(in):: v
    type(gene),intent(inout):: g
    integer:: dec
    
    dec = int((v -g%vmin)/(g%vmax-g%vmin)*(2**ga_nbits-1))
    call dec2bin(ga_nbits,dec,g%bits)
    return
  end subroutine var2gene
!=======================================================================
  subroutine indiv2vars(indiv,ndim,vars)
    type(individual),intent(in):: indiv
    integer,intent(in):: ndim
    real(8),intent(out):: vars(ndim)
    integer:: i
    
    do i=1,ndim
      call gene2var(indiv%genes(i),vars(i))
    enddo
    return
  end subroutine indiv2vars
!=======================================================================
  subroutine crossover(ind1,ind2,offspring)
!
!  Homogeneous crossover of two individuals to create an offspring
!  that has some similarities to the parents.
!
    use random
    type(individual),intent(in):: ind1,ind2
    type(individual),intent(inout):: offspring
    
    integer:: i
    type(gene):: g1,g2
    real(8):: v1,v2,r1
    
    do i=1,ga_ngenes
      g1 = ind1%genes(i)
      g2 = ind2%genes(i)
      v1 = ind1%fvalue
      v2 = ind2%fvalue
      r1 = log(v2+1d0)/(log(v1+1d0)+log(v2+1d0))
      do j=1,ga_nbits
        offspring%genes(i)%bits(j) = g1%bits(j)
        offspring%genes(i)%vmin = g1%vmin
        offspring%genes(i)%vmax = g1%vmax
        if( g1%bits(j).ne.g2%bits(j) .and. urnd() .gt. r1 ) then
          offspring%genes(i)%bits(j) = g2%bits(j)
        end if
      end do
!!$      print '(a,i4,2x,100i1)','i,bits=',i,&
!!$           (mod(offspring%genes(i)%bits(j),10),j=1,ga_nbits)
    end do
    return
  end subroutine crossover
!=======================================================================
  subroutine roulette_selection(nindivs,indivs,noffsp,offsprings,fbest)
!
!  Select individuals that are alive in the next generation according to
!  their evaulation values.
!  Selected ones are returned as an INDIVS array.
!  The best one is always selected at first.
!
    use random
    integer,intent(in):: nindivs,noffsp
    type(individual),intent(inout):: indivs(nindivs)
    type(individual),intent(in):: offsprings(noffsp)
    real(8),intent(in):: fbest

    integer:: i,j,k,l,m,n,ibest
    integer:: islct(nindivs)
    real(8):: fbestl

    integer,save:: nall
    real(8),save,allocatable:: probs(:)
    logical,save:: l1st = .true.
    type(individual),save,allocatable:: tmp_indivs(:)
    real(8),parameter:: pmax = 1d+10

    if( l1st ) then
      nall = nindivs + noffsp
      allocate(probs(nall),tmp_indivs(nindivs))
      l1st = .false.
    endif

!.....Compute all the probabilities using func values and temperature
    n = 0
    ibest = 0
    fbestl = 1d+30
    do i=1,nindivs
      n = n + 1
      if( indivs(i)%fvalue.lt.fbestl ) then
        fbestl = indivs(i)%fvalue
        ibest = n
      endif
      probs(n) = indivs(i)%fitness
!!$      if( trim(ga_fitness).eq.'exp' ) then
!!$        probs(n) = exp(-(indivs(i)%fvalue-fbest)/ga_temp)
!!$      else if( trim(ga_fitness).eq.'inv' ) then
!!$        probs(n) = 1d0/indivs(i)%fvalue
!!$      endif
!!$      print *,'n,fvalue,prob=',n,indivs(i)%fvalue,probs(n)
    enddo
    do i=1,noffsp
      n = n + 1
      if( offsprings(i)%fvalue.lt.fbestl ) then
        fbestl = offsprings(i)%fvalue
        ibest = n
      endif
      probs(n) = indivs(i)%fitness
!!$      if( trim(ga_fitness).eq.'exp' ) then
!!$        probs(n) = exp(-(offsprings(i)%fvalue-fbest)/ga_temp)
!!$      else if( trim(ga_fitness).eq.'inv' ) then
!!$        probs(n) = 1d0/offsprings(i)%fvalue
!!$      endif
!!$      print *,'n,fvalue,prob=',n,offsprings(i)%fvalue,probs(n)
    enddo

!.....Select individuals
    islct(1) = ibest
    probs(ibest) = 0d0
    do i=2,nindivs
      ptot = 0d0
      do j=1,nall
        ptot = ptot + probs(j)
      enddo
      prnd = urnd()*ptot
      ptot = 0d0
      do j=1,nall
        ptot = ptot +probs(j)
        if( prnd.lt.ptot ) then
          islct(i) = j
          probs(j) = 0d0
          exit
        endif
      enddo
    enddo

!!$    print *,'islct:'
!!$    do i=1,ga_nindivs
!!$      print *,'i,islct(i)=',i,islct(i)
!!$    enddo

!.....Replace indivs elements with selected ones
    do i=1,nindivs
      j = islct(i)
      if( j.le.nindivs ) then
        tmp_indivs(i) = indivs(j)
      else
        j = j - nindivs
        tmp_indivs(i) = offsprings(j)
      endif
    enddo
    do i=1,nindivs
      indivs(i) = tmp_indivs(i)
    enddo
    return
  end subroutine roulette_selection
!=======================================================================
  subroutine de(ndim,xbest,fbest,xranges,xtol,gtol,ftol,maxiter &
       ,iprint,iflag,myid,func,cfmethod,niter_eval,sub_eval,sub_ergrel)
!
! Differential evolution (DE) which does not use gradient information.
! DE itself is a serial code, but the function evaluation can be parallel.
!
    use random
    implicit none
    integer,intent(in):: ndim,iprint,myid,maxiter,niter_eval
    integer,intent(inout):: iflag
    real(8),intent(in):: xtol,gtol,ftol,xranges(2,ndim)
    real(8),intent(inout):: fbest,xbest(ndim)
    character(len=*),intent(in):: cfmethod
    interface
      subroutine func(n,x,ftrn,ftst)
        integer,intent(in):: n
        real(8),intent(in):: x(n)
        real(8),intent(out):: ftrn,ftst
      end subroutine func
      subroutine sub_eval(iter)
        integer,intent(in):: iter
      end subroutine sub_eval
      subroutine sub_ergrel(cadd)
        character(len=*),intent(in):: cadd
      end subroutine sub_ergrel
    end interface

    integer:: i,j,k,l,m,n,iter,i1,i2,ip,iq,ir,is
    real(8):: ftrn,ftst,fracl,fracg,lmdl,lmdg,w
    logical,save:: l1st = .true.
    integer:: iid,iidbest
    type(individual),allocatable:: indivs(:),offsprings(:)
    real(8),allocatable,dimension(:):: xtmp,xi,xp,xq,xr,xs,xbestl,xbestg&
         ,xl,xg,xd
    real(8),allocatable:: xpbest(:,:)
    character(len=128):: cadd

    integer,parameter:: io_indivs = 30
    character(len=128),parameter:: cf_indivs = 'out.de.individuals'
    integer,parameter:: io_steps = 31
    character(len=128),parameter:: cf_steps = 'out.de.generations'

    if( l1st ) then
!.....Allocate necessary memory spaces
      allocate(indivs(de_nindivs),offsprings(de_nindivs))
      allocate(xtmp(ndim),xi(ndim),xp(ndim),xq(ndim),xr(ndim),xs(ndim)&
           ,xbestl(ndim),xbestg(ndim),xl(ndim),xg(ndim),xd(ndim))
      allocate(xpbest(ndim,de_nindivs))
      do i=1,de_nindivs
        allocate(indivs(i)%genes(ndim),offsprings(i)%genes(ndim))
      enddo
!.....Initialize
      fracg = de_frac
      fracl = fracg
      if( de_lambda.le.0d0 ) de_lambda = de_frac
      lmdg = de_lambda
      lmdl = lmdg
      if( myid.eq.0 .and. iprint.ne.0 ) then
        write(6,*) ''
        write(6,'(a)') '------------------------------------------------------------------------'
        write(6,'(a)') '                          Differential Evolution'
        write(6,'(a)') '------------------------------------------------------------------------'
        print '(a,i4)',' Number of individuals = ',de_nindivs
        print '(a,f8.4)',' Fraction              =',de_frac
        print '(a,f8.4)',' Crossover rate        =',de_cross_rate
        print '(a,f8.4)',' wmin                  =',de_wmin
        print '(a,f8.4)',' wmax                  =',de_wmax
        print '(a,f8.4)',' fracg                 =',fracg
        print '(a,f8.4)',' fracl                 =',fracl
        print '(a,f8.4)',' lmdg                  =',lmdg
        print '(a,f8.4)',' lmdl                  =',lmdl

        print *,''
      endif
      l1st = .false.
    endif

    if( myid.eq.0 ) then
      open(io_indivs,file=cf_indivs,status='replace')
      write(io_indivs,'(a)') '# iid, fval, vars...'
      open(io_steps,file=cf_steps,status='replace')
      write(io_steps,'(a)') '# iter, iid, ftrn'
    endif
10  format(i6,es14.6,100f9.4)

    iter = 0

!.....Create population that includes some individuals
    iid = 0
    fbest = 1d+30
    do i=1,de_nindivs
      do j=1,ndim
        indivs(i)%genes(j)%vmin = xranges(1,j)
        indivs(i)%genes(j)%vmax = xranges(2,j)
        if( i.eq.1 ) then
          indivs(i)%genes(j)%val = xbest(j)
        else
          indivs(i)%genes(j)%val = xranges(1,j) + &
               (xranges(2,j)-xranges(1,j))*urnd()
        endif
      enddo
      do j=1,ndim
        xtmp(j) = indivs(i)%genes(j)%val
      enddo
      call func(ndim,xtmp,ftrn,ftst)
      iid = iid + 1
!.....Detect NaN and replace it with 1d+10
      if( ftrn*0d0 .ne. 0d0 ) then
        if( myid.eq.0 .and. iprint.ne.0 ) then
          write(6,'(a,2i10)')&
               ' [ftrn.eq.NaN] iter,iid = ',iter,iid
        endif
        ftrn = fupper_lim
      endif
      if( ftrn.gt.fupper_lim ) then
        ftrn = fupper_lim
      endif
      indivs(i)%fvalue = ftrn
      if( ftrn*0d0.ne.0d0 ) then
        indivs(i)%fitness = 0d0
      else
        indivs(i)%fitness = 1d0/ftrn
      endif
      indivs(i)%iid = iid
      if( myid.eq.0 ) write(io_indivs,10) iid,ftrn,xtmp(1:min(ndim,100))
      if( ftrn.lt.fbest ) then
        fbest = ftrn
        iidbest = iid
        xbest(1:ndim) = xtmp(1:ndim)
      endif

      if( iprint.ge.2 ) then
        write(cadd,'(i0)') iid
        call sub_ergrel(cadd)
      endif
      if( i.eq.1 ) call sub_eval(iter)
    enddo
    w = de_wmin + (de_wmax -de_wmin)*dble(iter)/maxiter
    if( maxiter.eq.0 ) w = de_wmin

    if( myid.eq.0 ) then
      write(6,'(a,i8,es12.4,f5.2,1x,100es12.4)') &
           " iter,fbest,w,fvals= ",&
           iter,fbest,w,(indivs(i)%fvalue,i=1,min(de_nindivs,10))
      do i=1,de_nindivs
        write(io_steps,'(2i8,es15.7)') iter, indivs(i)%iid, indivs(i)%fvalue
      enddo
    endif

!.....DE loop starts....................................................
    do iter=1,maxiter

      w = de_wmin + (de_wmax -de_wmin)*dble(iter)/maxiter
      call make_global_best(de_nindivs,indivs,ndim,xbestg)

!.....Loop for individuals
      do i=1,de_nindivs

        do j=1,ndim
          xi(j) = indivs(i)%genes(j)%val
        enddo

        if( trim(de_algo).eq.'local_neighbor' ) then

!.....Create a local vector
          ip = i+1
          if( ip.gt.de_nindivs ) ip = ip - de_nindivs
          iq = i-1
          if( iq.le.0 ) iq = iq + de_nindivs
          do j=1,ndim
            xp(j) = indivs(ip)%genes(j)%val
            xq(j) = indivs(iq)%genes(j)%val
          enddo
          xbestl(1:ndim) = 0d0
          call make_local_best(indivs(i),indivs(ip),indivs(iq),ndim,xbestl)
!!$        xl(1:ndim) = xi(1:ndim) +lmdl*(xbestl(1:ndim)-xi(1:ndim)) &
!!$             +fracl*(xp(1:ndim)-xq(1:ndim))
          xl(1:ndim) = xi(1:ndim) +lmdl*(xbestl(1:ndim)-xi(1:ndim)) &
               +urnd()*(xp(1:ndim)-xq(1:ndim))
!!$        print '(a,8es12.4)','xi,xp,xq,xbestl,fi,fp,fq,xl=' &
!!$             ,xi(1),xp(1),xq(1),xbestl(1)&
!!$             ,indivs(i)%fvalue,indivs(ip)%fvalue,indivs(iq)%fvalue,xl(1)
!.....Create a global vector
100       ir = int(urnd() *de_nindivs) +1
          if( ir.eq.i ) goto 100
110       is = int(urnd() *de_nindivs) +1
          if( is.eq.i .or. is.eq.ir ) goto 110
          do j=1,ndim
            xr(j) = indivs(ir)%genes(j)%val
            xs(j) = indivs(is)%genes(j)%val
          enddo
!!$        xg(1:ndim) = xi(1:ndim) +lmdg*(xbestg(1:ndim)-xi(1:ndim)) &
!!$             +fracg*(xr(1:ndim)-xs(1:ndim))
          xg(1:ndim) = xi(1:ndim) +lmdg*(xbestg(1:ndim)-xi(1:ndim)) &
               +urnd()*(xr(1:ndim)-xs(1:ndim))
!!$        print '(a,100f7.3)','xg(1:ndim)=',xg(1:ndim)
!.....Make the donor vector from the local and global vectors
          xd(1:ndim) = w*xg(1:ndim) +(1d0-w)*xl(1:ndim)

!.....Classical DE
        else
200       ip = int(urnd()*de_nindivs) +1
          if( ip.eq.i ) goto 200
210       ir = int(urnd()*de_nindivs) +1
          if( ir.eq.i .or. ir.eq.ip ) goto 210
220       is = int(urnd()*de_nindivs) +1
          if( is.eq.i .or. is.eq.ip .or. is.eq.ir )  goto 220
          do j=1,ndim
            xp(j) = indivs(ip)%genes(j)%val
            xr(j) = indivs(ir)%genes(j)%val
            xs(j) = indivs(is)%genes(j)%val
          enddo
!!$        xg(1:ndim) = xi(1:ndim) +lmdg*(xbestg(1:ndim)-xi(1:ndim)) &
!!$             +fracg*(xr(1:ndim)-xs(1:ndim))
!!$          xd(1:ndim) = xp(1:ndim) +de_frac *(xr(1:ndim)-xs(1:ndim))
          xd(1:ndim) = xp(1:ndim) +urnd()*de_frac *(xr(1:ndim)-xs(1:ndim))
        endif  ! de_algo

        iid = iid + 1
!.....Make a new candidate by the crossover of xd and xi
        do j=1,ndim
          if( urnd().lt.de_cross_rate ) then
            xtmp(j) = xd(j)
          else
            xtmp(j) = xi(j)
          endif
        enddo
        do j=1,ndim
          xtmp(j) = max(xtmp(j),indivs(i)%genes(j)%vmin)
          xtmp(j) = min(xtmp(j),indivs(i)%genes(j)%vmax)
!!$          offsprings(i)%genes(j)%val = xtmp(j)
        enddo
!!$        offsprings(i)%iid = iid
        call func(ndim,xtmp,ftrn,ftst)
!!$        offsprings(i)%fvalue = ftrn
!.....Detect NaN and replace it with fupper_lim
!!$        if( ftrn*0d0.ne.0d0 ) then
!!$          offsprings(i)%fitness = 0d0
!!$        else
!!$          offsprings(i)%fitness = 1d0/ftrn
!!$        endif
!!$        print *,'i,fvalue,ftrn=',i,indivs(i)%fvalue,ftrn
        if( ftrn.le.indivs(i)%fvalue .or. indivs(i)%fvalue*0d0.ne.0d0 ) then
          do j=1,ndim
            indivs(i)%genes(j)%val = xtmp(j)
          enddo
          if( ftrn.gt.fupper_lim ) ftrn = fupper_lim
          indivs(i)%fvalue = ftrn
          indivs(i)%fitness = 1d0/ftrn
        else
          cycle
        endif
        if( myid.eq.0 ) write(io_indivs,10) iid,ftrn,xtmp(1:min(ndim,100))
        if( ftrn.lt.fbest ) then
          fbest = ftrn
          iidbest = iid
          xbest(1:ndim) = xtmp(1:ndim)
          if( iprint.ge.2 ) then
            write(cadd,'(i0)') iid
            call sub_ergrel(cadd)
          endif
          call sub_eval(iid)
        endif
      enddo  ! loop over individuals

      if( myid.eq.0 ) then
        write(6,'(a,i8,es12.4,f5.2,1x,100es12.4)') &
             " iter,fbest,w,fvals= ",&
             iter,fbest,w,(indivs(i)%fvalue,i=1,min(de_nindivs,10))
        do i=1,de_nindivs
          write(io_steps,'(2i8,es15.7)') iter, indivs(i)%iid, indivs(i)%fvalue
        enddo
        flush(io_indivs)
        flush(io_steps)
      endif
    enddo
!.....DE loop ends......................................................

!.....Output information of the best
    if( myid.eq.0 ) then
      write(6,*) ''
      write(6,'(a)') ' The best one in this DE simulation run:'
      write(6,'(a,i8,1x,f0.4)') '   ID, f-value: ',iidbest,fbest
!!$      write(6,'(a,100f7.3)') '   Variables: ', xbest(1:min(ndim,100))
      write(6,*) ''
    endif

!.....Just for outputing out.erg.fin.1 of the best one     
    call func(ndim,xbest,ftrn,ftst)

    close(io_indivs)
    close(io_steps)
    return
  end subroutine de
!=======================================================================
  subroutine make_local_best(ind0,ind1,ind2,ndim,xbestl)
    type(individual),intent(in):: ind0,ind1,ind2
    integer,intent(in):: ndim
    real(8),intent(out):: xbestl(ndim)

    integer:: i
    real(8):: allf
    
    xbestl(1:ndim) = 0d0
    allf = ind0%fitness +ind1%fitness +ind2%fitness
    do i=1,ndim
      xbestl(i) = xbestl(i) + (&
           ind0%genes(i)%val *ind0%fitness  &
           +ind1%genes(i)%val *ind1%fitness  &
           +ind2%genes(i)%val *ind2%fitness ) /allf
    enddo
    return
  end subroutine make_local_best
!=======================================================================
  subroutine make_global_best(nindivs,indivs,ndim,xbestg)
    integer,intent(in):: nindivs,ndim
    type(individual),intent(in):: indivs(nindivs)
    real(8),intent(out):: xbestg(ndim)

    integer:: i,j
    real(8):: allf

    xbestg(1:ndim) = 0d0
    allf = 0d0
    do i=1,nindivs
      allf = allf +indivs(i)%fitness
    enddo
    do i=1,ndim
      do j=1,nindivs
        xbestg(i) = xbestg(i) +indivs(j)%fitness *indivs(j)%genes(i)%val
      enddo
      xbestg(i) = xbestg(i) /allf
    enddo
    return
  end subroutine make_global_best
!=======================================================================
  subroutine pso(ndim,xbest,fbest,xranges,xtol,gtol,ftol,maxiter &
       ,iprint,iflag,myid,func,cfmethod,niter_eval,sub_eval)
!
! Particle Swarm Optimization (PSO).
! DE itself is a serial code, but the function evaluation can be parallel.
!
    use random
    implicit none
    integer,intent(in):: ndim,iprint,myid,maxiter,niter_eval
    integer,intent(inout):: iflag
    real(8),intent(in):: xtol,gtol,ftol,xranges(2,ndim)
    real(8),intent(inout):: fbest,xbest(ndim)
    character(len=*),intent(in):: cfmethod
    interface
      subroutine func(n,x,ftrn,ftst)
        integer,intent(in):: n
        real(8),intent(in):: x(n)
        real(8),intent(out):: ftrn,ftst
      end subroutine func
      subroutine sub_eval(iter)
        integer,intent(in):: iter
      end subroutine sub_eval
    end interface

    integer:: i,j,k,l,m,n,iter
    real(8):: ftrn,ftst,tmp,xj,vj,r1,r2
    logical,save:: l1st = .true.
    integer:: iid,iidbest
    type(individual),allocatable:: indivs(:)
    real(8),allocatable:: xtmp(:),xpbest(:,:),fpbest(:)

    integer,parameter:: io_indivs = 30
    character(len=128),parameter:: cf_indivs = 'out.pso.individuals'
    integer,parameter:: io_fvalues = 31
    character(len=128),parameter:: cf_fvalues = 'out.pso.fvalues'
    integer,parameter:: io_steps = 32
    character(len=128),parameter:: cf_steps = 'out.pso.generations'

    if( l1st ) then
!.....Allocate necessary memory spaces
      allocate(indivs(pso_nindivs))
      allocate(xtmp(ndim),fpbest(pso_nindivs),xpbest(ndim,pso_nindivs))
      do i=1,pso_nindivs
        allocate(indivs(i)%genes(ndim),indivs(i)%vel(ndim))
      enddo
      if( myid.eq.0 .and. iprint.ne.0 ) then
        write(6,*) ''
        write(6,'(a)') '------------------------------------------------------------------------'
        write(6,'(a)') '                   Particle Swarm Optimization (PSO)'
        write(6,'(a)') '------------------------------------------------------------------------'
        print '(a,i4)',  ' Number of individuals = ',pso_nindivs
        print '(a,f8.4)',' w                     =',pso_w
        print '(a,f8.4)',' C1                    =',pso_c1
        print '(a,f8.4)',' C2                    =',pso_c2
        print *,''
      endif
      l1st = .false.
    endif

    if( myid.eq.0 ) then
      open(io_indivs,file=cf_indivs,status='replace')
      write(io_indivs,'(a)') '# iid, fval, vars...'
      open(io_fvalues,file=cf_fvalues,status='replace')
      write(io_fvalues,'(a)') '# iter, iid, ftrn'
      open(io_steps,file=cf_steps,status='replace')
      write(io_steps,'(a)') '# iter, iid, ftrn'
    endif
10  format(i6,es14.6,100f9.4)

    iter = 0

!.....Create population that includes some individuals
    iid = 0
    fbest = 1d+30
    do i=1,pso_nindivs
      do j=1,ndim
        indivs(i)%genes(j)%vmin = xranges(1,j)
        indivs(i)%genes(j)%vmax = xranges(2,j)
        if( i.eq.1 ) then
          indivs(i)%genes(j)%val = xbest(j)
        else
          indivs(i)%genes(j)%val = xranges(1,j) + &
               (xranges(2,j)-xranges(1,j))*urnd()
        endif
        indivs(i)%vel(j) = pso_vinimax*(urnd()-0.5d0)
      enddo
      do j=1,ndim
        xtmp(j) = indivs(i)%genes(j)%val
      enddo
      call func(ndim,xtmp,ftrn,ftst)
      iid = iid + 1
!.....Detect NaN and replace it with 1d+10
      if( ftrn*0d0 .ne. 0d0 ) then
        if( myid.eq.0 .and. iprint.ne.0 ) then
          write(6,'(a,2i10)')&
               ' [ftrn.eq.NaN] iter,iid = ',iter,iid
        endif
        ftrn = fupper_lim
      endif
      if( ftrn.gt.fupper_lim ) then
        ftrn = fupper_lim
      endif
      indivs(i)%fvalue = ftrn
      fpbest(i) = ftrn
      xpbest(1:ndim,i) = xtmp(1:ndim)
      if( ftrn*0d0.ne.0d0 ) then
        indivs(i)%fitness = 0d0
      else
        indivs(i)%fitness = 1d0/ftrn
      endif
      indivs(i)%iid = iid
      if( myid.eq.0 ) write(io_indivs,10) iid,ftrn,xtmp(1:min(ndim,100))
      if( ftrn.lt.fbest ) then
        fbest = ftrn
        iidbest = iid
        xbest(1:ndim) = xtmp(1:ndim)
      endif
      if( i.eq.1 ) call sub_eval(iter)
    enddo
    if( myid.eq.0 ) then
      write(6,'(a,i8,es12.4,1x,100es12.4)') &
           " iter,fbest,fvals= ",&
           iter,fbest,(indivs(i)%fvalue,i=1,min(ndim,10))
      write(io_fvalues,'(2i8, 100es12.4)') iter, indivs(i)%iid, &
           (indivs(i)%fvalue,i=1,pso_nindivs)
      do i=1,pso_nindivs
        write(io_steps,'(2i8,es15.7)') iter, indivs(i)%iid, indivs(i)%fvalue
      enddo
    endif

!.....PSO loop starts....................................................
    do iter=1,maxiter

!.....Loop for individuals
      do i=1,pso_nindivs

!.....Update velocity and position
        r1 = urnd()
        r2 = urnd()
        do j=1,ndim
          xj = indivs(i)%genes(j)%val
          vj = indivs(i)%vel(j)
          indivs(i)%vel(j) = vj &
               +pso_c1*r1*( xpbest(j,i) -xj ) &
               +pso_c2*r2*( xbest(j) -xj )
          indivs(i)%genes(j)%val = xj +indivs(i)%vel(j)
          xtmp(j) = indivs(i)%genes(j)%val
!.....Make sure the range of xtmp
          xtmp(j) = max(xtmp(j),indivs(i)%genes(j)%vmin)
          xtmp(j) = min(xtmp(j),indivs(i)%genes(j)%vmax)
        enddo
        call func(ndim,xtmp,ftrn,ftst)
        iid = iid + 1
!.....Detect NaN and replace it with 1d+10
        if( ftrn*0d0 .ne. 0d0 ) then
          ftrn = fupper_lim
        endif
        if( ftrn.gt.fupper_lim ) then
          ftrn = fupper_lim
        endif
        indivs(i)%fvalue = ftrn
        indivs(i)%iid = iid
        if( myid.eq.0 ) write(io_indivs,10) iid,ftrn,xtmp(1:min(ndim,100))
!.....Check the global best and update if needed
        if( ftrn.lt.fbest ) then
          fbest = ftrn
          iidbest = iid
          xbest(1:ndim) = xtmp(1:ndim)
          call sub_eval(iid)
        endif
!.....Check the particle best and update if needed
        if( ftrn.lt.fpbest(i) ) then
          fpbest(i)= ftrn
          xpbest(1:ndim,i) = xtmp(1:ndim)
        endif

      enddo  ! loop over individuals

      if( myid.eq.0 ) then
        write(6,'(a,i8,es12.4,1x,100es12.4)') &
             " iter,fbest,fvals= ",&
             iter,fbest,(indivs(i)%fvalue,i=1,min(ndim,10))
        write(io_fvalues,'(2i8, 100es12.4)') iter, indivs(i)%iid, &
             (indivs(i)%fvalue,i=1,pso_nindivs)
        do i=1,pso_nindivs
          write(io_steps,'(2i8,es15.7)') iter, indivs(i)%iid, indivs(i)%fvalue
        enddo
        flush(io_indivs)
        flush(io_steps)
        flush(io_fvalues)
      endif
    enddo
!.....DE loop ends......................................................

!.....Output information of the best
    if( myid.eq.0 ) then
      write(6,*) ''
      write(6,'(a)') ' The best one in this PSO simulation run:'
      write(6,'(a,i8,1x,f0.4)') '   ID, f-value: ',iidbest,fbest
!!$      write(6,'(a,100f7.3)')  '   Variables: ',xbest(1:ndim)
      write(6,*) ''

!!$      do i=1,pso_nindivs
!!$        print '(a,i4,f10.1,1x,100f9.4)', ' i,fpbest,xpbest =' ,i,fpbest(i),&
!!$             xpbest(1:ndim,i)
!!$      enddo
      print *,''
    endif

!.....Just for outputing out.erg.fin.1 of the best one     
    call func(ndim,xbest,ftrn,ftst)

    close(io_indivs)
    close(io_fvalues)
    close(io_steps)
    return
  end subroutine pso
!=======================================================================
  
end module
