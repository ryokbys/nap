module minimize
  use variables,only: dmem
  implicit none 
  save
!.....number of convergence criteria achieved
  integer:: numtol = 1

!.....penalty: lasso or ridge or smooth
  character(len=128):: cpena= 'none'
  character(len=128):: clinmin= 'backtrack'
  character(len=128):: cfsmode= 'grad'  ! [grad,grad0corr,df0corr]
  real(8):: pwgt = 1d-15

!.....SGD parameters
  integer:: nsgdbsize = 1
  integer:: nsgdbsnode = 1
  integer,allocatable:: ismask(:)
  character(len=128):: csgdupdate = 'normal'
  real(8):: sgd_rate_ini = 0.001d0
  real(8):: sgd_rate_fin = -0.001d0
  real(8):: sgd_eps = 1.0d-8
!.....Parameters for ADAM and AdaBound
  real(8):: adam_b1 = 0.9d0
  real(8):: adam_b2 = 0.999d0

!.....Group FS inner loop
  integer:: ninnergfs=100
  character(len=128):: cread_fsmask = ''
  character(len=128):: cfs_xrefresh = 'random' ! [zero, random, none]
  integer:: maxfsrefresh = 2

!.....Max iteration for line minimization
  integer:: niter_linmin   = 15
!.....Decreasing factor, should be < 1.0
  real(8):: fac_dec        = 0.2d0
!.....Increasing factor, should be > 1.0
  real(8):: fac_inc        = 5.0d0
!.....Armijo parameters
  real(8):: armijo_xi      = 1.0d-4
  real(8):: armijo_tau     = 0.5d0
  integer:: armijo_maxiter = 15

!.....CG
  integer:: icgbtype = 1 ! 1:FR, 2:PRP, 3:HS, 4:DY

!.....L-BFGS
  integer:: m_lbfgs   = 10

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
    
    ranges(1:2,1:ndim) = xranges(1:2,1:ndim)
    return
  end subroutine set_ranges
!=======================================================================
  subroutine wrap_ranges(ndim,x,xranges)
    implicit none
    integer,intent(in):: ndim
    real(8),intent(in):: xranges(2,ndim)
    real(8),intent(inout):: x(ndim)

    integer:: i

!!$    if( .not.allocated(ranges) ) then
!!$      print *,'ERROR: ranges is not allocated yet...'
!!$      stop
!!$    endif

    do i=1,ndim
      if( x(i).lt.xranges(1,i) ) then
        x(i) = xranges(1,i)
      else if( x(i).gt.xranges(2,i) ) then
        x(i) = xranges(2,i)
      endif
    enddo
    return
  end subroutine wrap_ranges
!=======================================================================
  subroutine write_status(ionum,myid,iprint,cpena,iter,ninner &
       ,ftrn,ftst,pval,xnorm,gnorm,dxnorm,fprev)
    integer,intent(in):: ionum,myid,iprint,iter,ninner
    character(len=128),intent(in):: cpena
    real(8),intent(in)::ftrn,ftst,pval,xnorm,gnorm,dxnorm,fprev
    
    if( myid.eq.0 ) then
      if( iprint.ge.1 ) then
        if( trim(cpena).eq.'ridge' ) then
          write(6,'(a,i5,i4,7es11.3)') &
               ' iter,ninner,ftrn,ftst,penalty,|x|,|g|,|dx|,|df|=' &
               ,iter,ninner,ftrn-pval,ftst &
               ,pval,xnorm,gnorm,dxnorm,abs(ftrn-fprev)
        else
          write(6,'(a,i5,i4,7es11.3)') &
               ' iter,ninner,ftrn,ftst,|x|,|g|,|dx|,|df|=' &
               ,iter,ninner,ftrn,ftst,xnorm,gnorm,dxnorm,abs(ftrn-fprev)
        endif
        call flush(6)
      endif
    endif
  end subroutine write_status
!=======================================================================
  subroutine check_converge(myid,iprint,cmin,xtol,gtol,ftol &
       ,dxnorm,gnorm,fdiff,nxtol,ngtol,nftol,iflag,lconverged)
    integer,intent(in):: myid,iprint
    real(8),intent(in):: xtol,gtol,ftol,dxnorm,gnorm,fdiff
    integer,intent(inout):: nxtol,ngtol,nftol,iflag
    logical,intent(out):: lconverged
    character(len=*):: cmin

    lconverged = .false.
    if( myid.eq.0 .and. iprint.gt.1 ) then
      print '(a,3es12.4)','  dxnorm,gnorm,fdiff=',dxnorm,gnorm,fdiff
    endif
    if( dxnorm.lt.xtol ) then
      nxtol = nxtol +1
      ngtol = 0
      nftol = 0
      if( nxtol.ge.numtol ) then
        if( myid.eq.0 .and. iprint.gt.0 ) then
          print '(a,i0,a)',' >>> '//trim(cmin)//' converged' &
               //' because xdiff < xtol over ',numtol,' times.'
          write(6,'(a,2es13.5)') '   dxnorm,xtol=',dxnorm,xtol
        endif
        iflag= iflag +1
        lconverged = .true.
!!$        x0(1:ndim)= x(1:ndim)
!!$        maxiter = iter
        return
      endif
    else if( gnorm.lt.gtol ) then
      ngtol = ngtol +1
      nxtol = 0
      nftol = 0
      if( ngtol.ge.numtol ) then
        if( myid.eq.0 ) then
          print '(a,i0,a)',' >>> '//trim(cmin)//' converged' &
               //' because gdiff < gtol over ',numtol,' times.'
          write(6,'(a,2es13.5)') '   gnorm,gtol=',gnorm,gtol
        endif
        iflag= iflag +2
        lconverged = .true.
!!$        x0(1:ndim)= x(1:ndim)
!!$        maxiter = iter
        return
      endif
    else if( fdiff.lt.ftol) then
      nftol= nftol +1
      nxtol = 0
      ngtol = 0
      if( nftol.ge.numtol ) then
        if( myid.eq.0 ) then
          print '(a,i0,a)',' >>> '//trim(cmin)//' converged' &
               //' because fdiff < ftol over ',numtol,' times.'
          write(6,'(a,2es13.5)') '   fdiff,ftol=',fdiff, ftol
        endif
        iflag= iflag +3
        lconverged = .true.
!!$        x0(1:ndim)= x(1:ndim)
!!$        maxiter = iter
        return
      endif
    else
      nxtol = 0
      ngtol = 0
      nftol = 0
    endif
    return
  end subroutine check_converge
!=======================================================================
  subroutine steepest_descent(ndim,x0,f,g,d,xranges,xtol,gtol,ftol,maxiter &
       ,iprint,iflag,myid,func,grad,cfmethod,niter_eval,sub_eval)
    implicit none
    integer,intent(in):: ndim,iprint,myid,niter_eval
    integer,intent(inout):: iflag,maxiter
    real(8),intent(in):: xtol,gtol,ftol,xranges(2,ndim)
    real(8),intent(inout):: f,x0(ndim),g(ndim),d(ndim)
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

    integer:: iter,niter,nxtol,ngtol,nftol
    real(8):: alpha,fp,gnorm,dxnorm,vnorm,ftst,pval
    real(8),save,allocatable:: s(:),y(:),x(:),gpena(:),gp(:),xp(:),dx(:)
    logical:: lconverged = .false. 

    if( myid.eq.0 ) then
      print *,''
      print *, '********************** Steepest Descent (SD) '&
           //  '**********************'
    endif

    if( .not.allocated(gpena) ) then
      allocate(gpena(ndim),dx(ndim),xp(ndim),x(ndim) &
           ,s(ndim),y(ndim),gp(ndim))
    endif

    iter= 0
    niter = 0
    x(:) = x0(:)
    call wrap_ranges(ndim,x,xranges)
    call func(ndim,x,f,ftst)
    call grad(ndim,x,g)
    gnorm= sqrt(sprod(ndim,g,g))
    vnorm= sqrt(sprod(ndim,x,x))
    d(1:ndim)= -g(1:ndim)
    call write_status(6,myid,iprint,cpena,iter,niter &
         ,f,ftst,pval,vnorm,gnorm,dxnorm,f)

    alpha = 1d0
    do iter=1,maxiter
      fp= f
      xp(:) = x(:)
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
      else if( trim(clinmin).eq.'two-point' ) then
!.....Like backtrack, but once alpha is defined, not to perform line minimization
!     no matter the func value gets larger.
!.....Currently this does NOT work well...
        if( iter.le.1 ) then  ! only 1st call, perform line minimization a bit
          call backtrack(ndim,x,xranges,d,f,ftst,alpha,iprint &
               ,iflag,myid,func,niter)
        else
          alpha = sprod(ndim,s,s)/sprod(ndim,s,y)
          alpha = max(min(alpha,1d0),xtol)
        endif
      else if( trim(clinmin).eq.'armijo' ) then
!!$        alpha = min(max(alpha,xtol*2d0)*2d0, 1d0)
!!$        alpha = max(alpha,xtol/gnorm)*2d0
        alpha = alpha *fac_inc
        call armijo_search(ndim,x,xranges,d,f,ftst,g,alpha,iprint &
             ,iflag,myid,func,niter)
      else ! Default = backtrack
!.....Increase alpha a bit every step,
!.....alpha is to be decreased in subroutine backtrack to decrease func value.
!!$        alpha = min(max(alpha,xtol*2d0)*2d0, 1d0)
!!$        alpha = max(alpha,xtol/gnorm)*2d0
        alpha = alpha *fac_inc
        call backtrack(ndim,x,xranges,d,f,ftst,alpha,iprint &
             ,iflag,myid,func,niter)
      endif
      if( iflag/100.ne.0 ) then
        if( myid.eq.0 ) then
          print *,'ERROR: iflag/100.ne.0 in SD !!!'
        endif
        return
      endif
!.....evaluate statistics at every niter_eval
      if( mod(iter,niter_eval).eq.0 ) &
           call sub_eval(iter)
!.....Update x
      x(1:ndim)= x(1:ndim) +alpha*d(1:ndim)
      s(1:ndim) = alpha*d(1:ndim)
      dx(:) = x(:) -xp(:)
      xp(:) = x(:)
      dxnorm = sqrt(sprod(ndim,dx,dx))
      call wrap_ranges(ndim,x,xranges)
      gp(:) = g(:)
      call grad(ndim,x,g)
      y(:) = g(:) -gp(:)
      gnorm= sqrt(sprod(ndim,g,g))
      d(1:ndim)= -g(1:ndim)
!!$      g(1:ndim)= -g(1:ndim)/gnorm
!!$      gnorm= gnorm/ndim
      vnorm= sqrt(sprod(ndim,x,x))
      call write_status(6,myid,iprint,cpena,iter,niter &
           ,f,ftst,pval,vnorm,gnorm,dxnorm,fp)
      call check_converge(myid,iprint,'SD',xtol,gtol,ftol &
           ,dxnorm,gnorm,abs(f-fp),nxtol,ngtol,nftol,iflag,lconverged)
      if( lconverged ) then
        x0(:) = x(:)
        maxiter = iter
        return
      endif
    enddo

    x0(:) = x(:)
    return
  end subroutine steepest_descent
!=======================================================================
  subroutine sgd(ndim,x0,f,g,u,xranges,xtol,gtol,ftol,maxiter,iprint &
       ,iflag,myid,mpi_world,mynsmpl,myntrn,isid0,isid1 &
       ,func,grad,cfmethod,niter_eval,sub_eval)
!
! Stochastic gradient decent (SGD)
!
    integer,intent(in):: ndim,iprint,niter_eval,myid,mpi_world, &
         mynsmpl,myntrn,isid0,isid1
    integer,intent(inout):: iflag,maxiter
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
    include 'mpif.h'
    real(8),parameter:: tiny  = 1d-14

    integer:: i,ismpl,iter,niter,nftol,ngtol,nxtol
    integer:: ninnerstp,innerstp,ierr,itmp
    real(8):: gnorm,xnorm,dxnorm,pval,sgd_rate,sgd_ratei,pvaltmp
    real(8):: fp,ftmp,ftst,ftsttmp
    real(8):: rate_upper, rate_lower
    real(8),allocatable:: x(:),dx(:),rm(:),rmh(:),gpena(:),gtmp(:) &
         ,gp(:),v(:),vh(:),xp(:),gpenatmp(:)
    integer,allocatable:: imaskarr(:)
    logical:: lconverged

    if( .not.allocated(x) ) then
      nsgdbsize = min(nsgdbsnode,myntrn)
      if(myid.eq.0) then
        print *,''
        print *, '************************ Stochastic gradient descent (SGD) '&
             //'************************'
        print *,'   Update method: ',trim(csgdupdate)
        if( trim(csgdupdate).ne.'normal' .and. trim(csgdupdate).ne.'adam' &
             .and. trim(csgdupdate).ne.'adabound' ) then
          if( myid.eq.0 ) then
            print *,'   WARNING: update method '//trim(csgdupdate)//' is not available.'
            print *,'            Use normal sgd instead...'
          endif
        endif
        print '(a,i0)','    SGD batch size per node = ',nsgdbsize
        if( trim(csgdupdate).eq.'adam' ) then
          print '(a,es11.3)','    epsilon value = ',sgd_eps
          print '(a,es11.3)','    beta1 = ',adam_b1
          print '(a,es11.3)','    beta2 = ',adam_b2
        else if( trim(csgdupdate).eq.'adabound' ) then
          print '(a,es11.3)','    epsilon value = ',sgd_eps
          print '(a,es11.3)','    beta1 = ',adam_b1
          print '(a,es11.3)','    beta2 = ',adam_b2
        endif
        if( sgd_rate_fin.le.0d0 ) then
          print '(a,es11.3)','    learning rate = ',sgd_rate_ini
        else
          print '(a,es11.3)','    initial learning rate = ',sgd_rate_ini
          print '(a,es11.3)','    final learning rate  = ',sgd_rate_fin
        endif
        print *,''
      endif
      allocate(x(ndim),dx(ndim),rm(ndim),rmh(ndim),gpena(ndim) &
           ,gp(ndim),gtmp(ndim),v(ndim),vh(ndim),xp(ndim),gpenatmp(ndim))
      allocate(ismask(isid0:isid1))
      sgd_rate = sgd_rate_ini
    endif

!.....Initialization
    nftol= 0
    ngtol= 0
    nxtol= 0
    rm(:) = 0d0
    v(:) = 0d0
    x(1:ndim)= x0(1:ndim)
    allocate(imaskarr(mynsmpl))
    ninnerstp = myntrn /nsgdbsize
    if( mod(myntrn,nsgdbsize) .ne. 0 ) ninnerstp = ninnerstp + 1
    itmp = ninnerstp
    call mpi_allreduce(itmp,ninnerstp,1,mpi_integer,mpi_max,mpi_world,ierr)

!.....Unset mask to compute all the samples at the first evaluation
    ismask(:) = 0
    call wrap_ranges(ndim,x0,xranges)
    call func(ndim,x0,f,ftst)
    call grad(ndim,x0,g)
    gnorm= sqrt(sprod(ndim,g,g))
    xnorm= sqrt(sprod(ndim,x,x))
    dxnorm = 0d0

    iter= 0
    niter = 0
    call write_status(6,myid,iprint,cpena,iter,ninnerstp &
         ,f,ftst,pval,xnorm,gnorm,dxnorm,f)

    call sub_eval(0)

!.....One iteration includes evaluation of all the training data.
    do iter=1,maxiter
      fp = f
      gp(:) = g(:)
      xp(:) = x(:)
      
!.....Gradual increasing/descreasing learning rate if sgd_rate_fin is given
      if( sgd_rate_fin .gt. 0d0 ) then
        sgd_rate = sgd_rate_ini +(sgd_rate_fin -sgd_rate_ini)/maxiter *(iter-1)
      endif

!.....Make imaskarr that contains the order of samples to be computed in each innerstp
      imaskarr(:) = 0
      call get_order_iarr(myntrn,nsgdbsize,imaskarr)
!!$      print *,'myid,iter,imaskarr(:)=',myid,iter,imaskarr(:)

!.....Inner loop for batch process
      do innerstp = 1,ninnerstp
!.....Unmask only the samples whose imaskarr(i)==innerstp
        ismask(:) = 1
        do i=1,myntrn  ! All the test samples remain masked
          ismpl = isid0 + i -1
          if( imaskarr(i).eq.innerstp ) then
            ismask(ismpl) = 0
!!$            print *,'myid,ismpl= ',myid,ismpl
          endif
        enddo

        call wrap_ranges(ndim,x,xranges)
        call func(ndim,x,f,ftst)
        call grad(ndim,x,g)
        gnorm= sqrt(sprod(ndim,g,g))
!!$        print *,'myid,innerstp,gnorm=',myid,innerstp,gnorm

!.....Compute step size of x
        if( trim(csgdupdate).eq.'adam' .or. trim(csgdupdate).eq.'Adam' ) then
          rm(:) = adam_b1*rm(:) +(1d0 -adam_b1)*g(:)
          v(:) = adam_b2*v(:) +(1d0 -adam_b2)*g(:)*g(:)
          rmh(:) = rm(:)/(1d0-adam_b1**iter)
          vh(:) = v(:)/(1d0-adam_b2**iter)
          dx(:) = -sgd_rate*rmh(:)/(sqrt(vh(:)) +sgd_eps)
        else if( trim(csgdupdate).eq.'adabound' ) then
          rm(:) = adam_b1*rm(:) +(1d0 -adam_b1/iter)*g(:)
          v(:) = adam_b2*v(:) +(1d0 -adam_b2)*g(:)*g(:)
          rmh(:) = rm(:)/(1d0-adam_b1**iter)
          vh(:) = v(:)/(1d0-adam_b2**iter)
          rate_lower = sgd_rate *(1d0 -1d0/((1d0 -adam_b2)*iter+1d0))
          rate_upper = sgd_rate *(1d0 +1d0/((1d0 -adam_b2)*iter))
          if( iprint.gt.1 ) print '(a,i6,2es12.4)','iter,lower,upper=',iter,rate_lower,rate_upper
          do i=1,ndim
            sgd_ratei = sgd_rate/(sqrt(vh(i)) +sgd_eps)
            sgd_ratei = min(max(rate_lower,sgd_ratei),rate_upper) /sqrt(dble(iter))
            dx(i) = -sgd_ratei*rmh(i)
          enddo
        else  ! normal SGD
          if( gnorm/xnorm .gt. 1d0 ) g(:) = g(:) /gnorm *xnorm
          dx(:) = -sgd_rate *g(:)
        endif

!.....Update x
!!$        if( trim(cpena).eq.'ridge' ) dx(:) = dx(:) -gpena(:)
!!$        if( trim(cpena).eq.'lasso' .or. trim(cpena).eq.'glasso' ) then
!!$          call soft_threshold(ndim,x,dx,1d0)
!!$        else
          x(1:ndim)= x(1:ndim) +dx(1:ndim)
!!$        endif
      enddo  ! innerstp

      call wrap_ranges(ndim,x,xranges)
      dx(:) = x(:) -xp(:)
      dxnorm = sqrt(sprod(ndim,dx,dx))
      xnorm= sqrt(sprod(ndim,x,x))
      gnorm= sqrt(sprod(ndim,g,g))
!!$      print *,'myid,iter,gnorm=',myid,iter,gnorm

!.....Evaluate statistics at every niter_eval.
      if( mod(iter,niter_eval).eq.0 ) then
!.....Before the output, compute all the remaining samples
        do ismpl=isid0,isid1
          ismask(ismpl) = mod(ismask(ismpl)+1,2)
        enddo
        call wrap_ranges(ndim,x,xranges)
        call func(ndim,x,ftmp,ftst)
        call grad(ndim,x,gtmp)
        gnorm= sqrt(sprod(ndim,gtmp,gtmp))
        call sub_eval(iter)
      endif
      call write_status(6,myid,iprint,cpena,iter,ninnerstp &
           ,f,ftst,pval,xnorm,gnorm,dxnorm,fp)
      call check_converge(myid,iprint,'SGD',xtol,gtol,ftol &
           ,dxnorm,gnorm,abs(f-fp),nxtol,ngtol,nftol,iflag,lconverged)
      if( lconverged ) then
        x0(:) = x(:)
        maxiter = iter
        deallocate(x,dx,rm,rmh,gpena)
        return
      endif
    enddo  ! iter

    x0(1:ndim)= x(1:ndim)
    deallocate(x,dx,rm,rmh,gpena,imaskarr)
    return
  end subroutine sgd
!=======================================================================
  subroutine get_order_iarr(ndim,nbatch,iarr)
!
!  Create an array of length ndim that has the order of samples to be computed.
!
    use random
    integer,intent(in):: ndim,nbatch
    integer,intent(out):: iarr(ndim)
    integer:: jarr(ndim)
    integer:: i,j,k,l,m,inc,idx

    do i=1,ndim
      jarr(i) = i
    enddo
    inc = 0
    l=ndim
    do while(.true.)
      inc = inc + 1
      if( l.le.nbatch ) then
        do k=1,l
          idx = jarr(k)
          iarr(idx) = inc
        enddo
        exit
      else  ! if not l.le.nbatch
        do j=1,nbatch
          k = l*urnd() +1
          idx = jarr(k)
          iarr(idx) = inc
          do m=k,l-1
            jarr(m) = jarr(m+1)
          enddo
          l= l-1
        enddo
      endif
    end do
    return
  end subroutine get_order_iarr
!=======================================================================
  subroutine get_uniq_iarr(n,m,iarr)
!
! Create an array with length m which includes a unique set of integers
! randomly chosen from from 1 to n.
!
    use random
    integer,intent(in):: n,m
    integer,intent(out):: iarr(m)
    integer:: jarr(n)
    integer:: i,j,k,l

    do i=1,n
      jarr(i) = i
    enddo
    l=n
    do i=1,m
      j = l*urnd()+1
      iarr(i)= jarr(j)
      do k=j,l-1
        jarr(k)= jarr(k+1)
      enddo
      l= l-1
    enddo
    return
  end subroutine get_uniq_iarr
!=======================================================================
  subroutine cg(ndim,x0,f,g,u,xranges,xtol,gtol,ftol,maxiter,iprint,iflag,myid &
       ,func,grad,cfmethod,niter_eval,sub_eval)
!
!  Conjugate gradient minimization
!
    integer,intent(in):: ndim,iprint,myid,niter_eval
    integer,intent(inout):: iflag,maxiter
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
    real(8),parameter:: t_DL   = 0.2d0
    logical:: ltwice = .false.

    integer:: iter,nftol,ngtol,nxtol,niter
    real(8):: alpha,fp,gnorm,gnormp,vnorm,beta,pval,sgnorm,ftst,dxnorm
    real(8),save,allocatable:: gpena(:),gp(:),y(:),xp(:),s(:),dx(:),uu(:),x(:)
    logical:: lconverged = .false.
    integer:: mem

    if( myid.eq.0 ) then
      print *,''
      print *, '********************** Conjugate Gradient (CG) '&
           //  '**********************'
      mem = ndim*8*8
      dmem = dmem +mem
      if( mem > 1000000000 ) then
        print '(a,f0.3,a)',' Memory for CG = ',float(mem)/1000000000,' GB'
      else
        print '(a,f0.3,a)',' Memory for CG = ',float(mem)/1000000,' MB'
      endif
    endif

    if( .not.allocated(gpena) ) allocate(gpena(ndim),gp(ndim)&
         ,y(ndim),xp(ndim),s(ndim),dx(ndim),uu(ndim),x(ndim))

!.....Initialize alpha (line minimization factor)
    alpha = 1d0
    
    iter= 0
    niter = 0
    nftol= 0
    nxtol= 0
    gpena(1:ndim)= 0d0
    x(:) = x0(:)
    call wrap_ranges(ndim,x,xranges)
    call func(ndim,x,f,ftst)
    call grad(ndim,x,g)
    gnorm= sprod(ndim,g,g)
    sgnorm= sqrt(gnorm)
    vnorm= sqrt(sprod(ndim,x,x))
    call write_status(6,myid,iprint,cpena,iter,niter &
         ,f,ftst,pval,vnorm,sgnorm,dxnorm,f)
    u(1:ndim)= -g(1:ndim)

    call sub_eval(0)
    do iter=1,maxiter
      fp= f
      xp(1:ndim)= x(1:ndim)
!!$!.....normalize u-vector only for line search
!!$      unorm = sqrt(sprod(ndim,u,u))
!!$      uu(1:ndim) = u(1:ndim)/unorm
!.....line minimization
      if( trim(clinmin).eq.'quadratic' ) then
        call quad_interpolate(ndim,x,u,f,ftst,xtol,gtol,ftol &
             ,alpha,iprint,iflag,myid,func)
!.....if quad interpolation failed, perform golden section
        if( iflag/100.ne.0 ) then
          iflag= iflag -(iflag/100)*100
          if(myid.eq.0) then
            print *,'since quad_interpolate failed, call golden_section.'
          endif
          call golden_section(ndim,x,u,f,ftst,xtol,gtol,ftol &
               ,alpha,iprint,iflag,myid,func)
        endif
      else if ( trim(clinmin).eq.'golden') then
        call golden_section(ndim,x,u,f,ftst,xtol,gtol,ftol,alpha &
             ,iprint,iflag,myid,func)
      else if( trim(clinmin).eq.'armijo' ) then
!.....To enhance the convergence in Armijo search,
!.....use the history of previous alpha by multiplying 2
!.....avoiding constant decrease, but alpha should not be greater than 1.
!!$        alpha = min(max(alpha,xtol*2d0)*2d0, 1d0)
!!$        alpha = max(alpha,xtol/gnorm)*2d0
        alpha = alpha *fac_inc
        call armijo_search(ndim,x,xranges,u,f,ftst,g,alpha,iprint &
             ,iflag,myid,func,niter)
      else ! backtrack (default)
!!$        alpha = min(max(alpha,xtol*2d0)*2d0, 1d0)
!!$        alpha = max(alpha,xtol/gnorm)*2d0
        alpha = alpha *fac_inc
        call backtrack(ndim,x,xranges,u,f,ftst,alpha,iprint &
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
            print *,'>>> Initialize gg because alpha was not found.'
          endif
          alpha = 1d0
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

!.....Update x
      x(1:ndim)= x(1:ndim) +alpha*u(1:ndim)
      call wrap_ranges(ndim,x,xranges)

      dx(1:ndim)= x(1:ndim) -xp(1:ndim)
      gnormp= gnorm
      gp(1:ndim)= g(1:ndim)
      call grad(ndim,x,g)
!!$      if( trim(cpena).eq.'ridge' ) g(1:ndim)= g(1:ndim) +gpena(1:ndim)
      gnorm= sprod(ndim,g,g)
      sgnorm= sqrt(gnorm)
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
      vnorm= sqrt(sprod(ndim,x,x))
      dxnorm= sqrt(sprod(ndim,dx,dx))
      call write_status(6,myid,iprint,cpena,iter,niter &
           ,f,ftst,pval,vnorm,sgnorm,dxnorm,fp)
      call check_converge(myid,iprint,'CG',xtol,gtol,ftol &
           ,dxnorm,sgnorm,abs(f-fp),nxtol,ngtol,nftol,iflag,lconverged)
      if( lconverged ) then
        x0(:) = x(:)
        maxiter = iter
        return
      endif
    enddo

    x0(:) = x(:)
    return
  end subroutine cg
!=======================================================================
  subroutine qn(ndim,x0,f,g,u,xranges,xtol,gtol,ftol,maxiter, &
       iprint,iflag,myid,func,grad,cfmethod,niter_eval,sub_eval)
!
!  Limited-memory Broyden-Fletcher-Goldfarb-Shanno type of Quasi-Newton method.
!  Since BFGS is a bit too heavy for high-dimension parameters (> 10,000),
!  adopt L-BFGS and avoid allocating ndim*ndim array.
!
    integer,intent(in):: ndim,iprint,myid,niter_eval
    integer,intent(inout):: iflag,maxiter
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
    real(8),save,allocatable:: x(:),gp(:),s(:),y(:),dx(:)
    real(8),save,allocatable:: gg(:,:),ggy(:)
    real(8),save,allocatable:: si(:,:),yi(:,:),q(:), &
         rho(:),ai(:),bi(:)
    real(8):: tmp1,tmp2,b,sy,syi,fp,alpha,gnorm,ynorm,vnorm,pval &
         ,estmem,ftst,dxnorm,yy
    integer:: i,j,iter,nftol,ngtol,nxtol,mem,niter
    logical:: lconverged = .false.
    logical:: limited_mem = .true.  ! L-BFGS

    if( .not.allocated(x) ) then
      allocate(x(ndim),dx(ndim),gp(ndim))
      estmem = 8*3*ndim
      if( trim(cfmethod).eq.'bfgs' .or. trim(cfmethod).eq.'BFGS' ) then ! BFGS
        limited_mem = .false.
        allocate(gg(ndim,ndim),s(ndim),y(ndim),ggy(ndim))
        estmem = estmem +8*(ndim**2 +3*ndim)
      else  ! L-BFGS
        limited_mem = .true.
        allocate(rho(m_lbfgs),si(ndim,m_lbfgs),yi(ndim,m_lbfgs), &
             q(ndim),ai(m_lbfgs),bi(m_lbfgs))
        estmem = estmem +8*(3*m_lbfgs +ndim +2*ndim*m_lbfgs)
        yi(:,:) = 0d0
        si(:,:) = 0d0
        rho(:) = 0d0
      endif
      if(myid.eq.0) then
        print *,''
        print *, '******************************* QN(BFGS) '&
             //'*******************************'
!!$        estmem = (ndim*ndim +ndim*6)*8
!!$        estmem = 8*(6*ndim +2*ndim*m_lbfgs +3*m_lbfgs)
        dmem = dmem +estmem
        mem= estmem/1000/1000
        if( mem.eq.0 ) then
          mem= estmem/1000
          write(6,'(a,i6,a)') ' Memory for BFGS = ' &
               ,int(estmem/1000),' kB'
        else
          write(6,'(a,i6,a)') ' Memory for BFGS = ' &
               ,int(estmem/1000/1000),' MB'
        endif
      endif
    endif

!.....initialize alpha (line minimization factor)
    alpha = 1d0

    nftol= 0
    ngtol= 0
    nxtol= 0
    if( .not. limited_mem ) then
!.....initial G = I
      gg(1:ndim,1:ndim)= 0d0
      do i=1,ndim
        gg(i,i)= 1d0
      enddo
    endif

    call wrap_ranges(ndim,x0,xranges)
    call func(ndim,x0,f,ftst)
    call grad(ndim,x0,g)
    gnorm= sqrt(sprod(ndim,g,g))
    x(1:ndim)= x0(1:ndim)
    vnorm= sqrt(sprod(ndim,x,x))
    dxnorm = 0d0

    iter= 0
    niter = 0
    call write_status(6,myid,iprint,cpena,iter,niter &
         ,f,ftst,pval,vnorm,gnorm,dxnorm,f)

    call sub_eval(0)
!.....Initial direction assuming H_0 == I in L-BFGS
    if( limited_mem ) u(:) = -g(:)
    do iter=1,maxiter
      if( .not.limited_mem ) then
        u(1:ndim)= 0d0
        do i=1,ndim
          u(1:ndim)= u(1:ndim) -gg(1:ndim,i)*g(i)
        enddo
      endif
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
      else if( trim(clinmin).eq.'golden') then
        call golden_section(ndim,x,u,f,ftst,xtol,gtol,ftol,alpha &
             ,iprint,iflag,myid,func)
      else if( trim(clinmin).eq.'armijo' ) then
!.....To enhance the convergence in Armijo search,
!.....use the history of previous alpha by multiplying 2
!.....avoiding constant decrease, but alpha should not be greater than 1.
!!$        alpha = min(max(alpha,xtol*2d0)*2d0, 1d0)
!!$        alpha = max(alpha,xtol)*2d0
        alpha = alpha *fac_inc
        call armijo_search(ndim,x,xranges,u,f,ftst,g,alpha,iprint &
             ,iflag,myid,func,niter)
      else ! backtrack (default)
!!$        alpha = min(max(alpha,xtol*2d0)*2d0, 1d0)
!        alpha = max(alpha,xtol)*2d0
        alpha = alpha *fac_inc
        call backtrack(ndim,x,xranges,u,f,ftst,alpha,iprint &
             ,iflag,myid,func,niter)
      endif
!!$      if(myid.eq.0) print *,'armijo steps, alpha=',niter,alpha
      if( iflag/100.ne.0 ) then
        if( ltwice ) then
          x0(1:ndim)= x(1:ndim)
          if(myid.eq.0) then
            print *,'>>> Line_minimization failed twice continuously.'
          endif
          return
        else
          ltwice= .true.
          if(myid.eq.0) then
            print *,'>>> Initialize gg because alpha was not found.'
          endif
          alpha= 1d0  ! reset alpha to 1
          if( .not.limited_mem ) then
            gg(1:ndim,1:ndim)= 0d0
            do i=1,ndim
              gg(i,i)= 1d0
            enddo
          endif
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
!.....Update x
      x(1:ndim)= x(1:ndim) +alpha*u(1:ndim)
      call wrap_ranges(ndim,x,xranges)

      dx(1:ndim)= x(1:ndim) -x0(1:ndim)
      x0(1:ndim)= x(1:ndim)
      call grad(ndim,x,g)
      gnorm= sqrt(sprod(ndim,g,g))
      vnorm= sqrt(sprod(ndim,x,x))
      dxnorm= sqrt(sprod(ndim,dx,dx))
      call write_status(6,myid,iprint,cpena,iter,niter &
           ,f,ftst,pval,vnorm,gnorm,dxnorm,fp)
      call check_converge(myid,iprint,'QN',xtol,gtol,ftol &
           ,dxnorm,gnorm,abs(f-fp),nxtol,ngtol,nftol,iflag,lconverged)
      if( lconverged ) then
        x0(:) = x(:)
        maxiter = iter
        return
      endif

      
      if( limited_mem ) then
!.....L-BFGS algorithm, see wikipedia page for detail
        
!.....Newest data is at s(:,m_lbfgs), oldest data is at s(:,1)
        do i=2,m_lbfgs
          si(:,i-1) = si(:,i)
          yi(:,i-1) = yi(:,i)
          rho(i-1) = rho(i)
        enddo
        si(:,m_lbfgs)= alpha *u(:)
        yi(:,m_lbfgs)= g(:) -gp(:)
!.....Compute current rho, s*y, y*y
        rho(m_lbfgs) = 0d0
        yy = 0d0
        sy = 0d0
        rho(m_lbfgs) = 1d0/sprod(ndim,yi(1,m_lbfgs),si(1,m_lbfgs))
        sy = sprod(ndim,si(1,m_lbfgs-1),yi(1,m_lbfgs-1))
        yy = sprod(ndim,yi(1,m_lbfgs-1),yi(1,m_lbfgs-1))
        if( iter.eq.1 ) then
          u(:) = -g(:)
          cycle
        endif
!.....L-BFGS update
        q(:) = g(:)
        do i= m_lbfgs-1,1,-1
          ai(i) = rho(i) *sprod(ndim,si(1,i),q)
          q(:) = q(:) -ai(i)*yi(:,i)
        enddo
!.....H_k^0
        q(:) = sy/yy *q(:)
        do i= 1,m_lbfgs-1
          bi(i) = rho(i) *sprod(ndim,yi(1,i),q)
          q(:) = q(:) +(ai(i)-bi(i))*si(:,i)
        enddo
        u(:) = -q(:)
        
      else  ! BFGS
        s(:)= alpha *u(:)
        y(:)= g(:) -gp(:)
        ynorm= sprod(ndim,y,y)
        if( ynorm.lt.1d-14 .or. dxnorm.lt.xtol .or. gnorm.lt.gtol &
             .or. abs(f-fp).lt.ftol ) then
          if(myid.eq.0) then
            print *,'>>> Initialize gg'
          endif
          gg(1:ndim,1:ndim)= 0d0
          do i=1,ndim
            gg(i,i)= 1d0
          enddo
          cycle
        endif

!.....update matrix gg in BFGS
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
      endif

    enddo  ! iter

    x0(1:ndim)= x(1:ndim)
    return
  end subroutine qn
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
    integer:: iter

    call func(ndim,x0+a*d,fa,fta)
    call func(ndim,x0+b*d,fb,ftb)
    iter= 0
10  continue
    iter= iter +1
    if( iter.gt.MAXITER ) then
      if( myid.eq.0 ) then
        print *,'WARNING: iter.gt.MAXITER in get_bracket'
        print *,'  Search direction may not be a descent direction.'
      endif
      iflag= iflag +1000
      return
!!$      stop
    endif
    if( abs(b-a).lt.1d-12) then
      if( myid.eq.0 ) then
        print *,'WARNING: a and b is too close in get_bracket'
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

    interface
      subroutine func(n,x,ftrn,ftst)
        integer,intent(in):: n
        real(8),intent(in):: x(n)
        real(8),intent(out):: ftrn,ftst
      end subroutine func
    end interface

    integer:: iter,imin,imax,ix
    real(8):: r,q,fmin,fmax,dmin,dmax,d,xmin
    real(8):: xi(4),fi(4),fti(4)
    
    xi(1)= 0d0
    xi(2)= xi(1) +STP0
    call get_bracket(ndim,x0,g,xi(1),xi(2),xi(3),fi(1),fi(2),fi(3) &
         ,fti(1),fti(2),fti(3) &
         ,iprint,iflag,myid,func)
    if( iflag/1000.ne.0 ) return
    
    iter= 0
10  continue
    iter= iter +1
    if( iter.gt.niter_linmin ) then
      if( myid.eq.0 ) then
        print *,'WARNING: iter.gt.NITER_LINMIN in quad_interpolate !!!'
        print *,'  iter,niter_linmin= ',iter,niter_linmin
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
    if( iter.gt.niter_linmin ) then
      if( myid.eq.0 ) then
        print *,'WARNING: iter.gt.NITER_LINMIN in golden_section.'
        print *,'  iter,niter_linmin = ',iter,niter_linmin
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
  subroutine armijo_search(ndim,x0,xranges,d,f,ftst,g,alpha,iprint &
       ,iflag,myid,func,niter)
!  
!  1D search using Armijo rule.
!
    implicit none
    integer,intent(in):: ndim,iprint,myid
    integer,intent(inout):: iflag,niter
    real(8),intent(in):: x0(ndim),g(ndim),d(ndim),xranges(2,ndim)
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
    integer:: iter
    real(8):: alphai,xigd,f0,fi,pval,fp,pvalp,alphap,ftsti
    real(8),allocatable,dimension(:):: x1(:),gpena(:)
    logical,save:: l1st = .true.

    if( l1st ) then
      if( myid.eq.0 .and. iprint.gt.0 ) then
        write(6,'(a)') ' Armijo rule parameters:'
        write(6,'(a,es12.4)') '   c       = ',armijo_xi
        write(6,'(a,f10.4)') '   tau     = ',armijo_tau
        write(6,'(a,i5)')   '   maxiter = ',niter_linmin
      endif
      l1st = .false.
    endif

    if( .not. allocated(x1)) allocate(x1(ndim),gpena(ndim))
    xigd= sprod(ndim,g,d)*armijo_xi
    if( xigd.gt.0d0 ) then
      iflag= iflag + 100
      if( myid.eq.0 .and. iprint.gt.0 ) print *,'WARNING: g*d > 0.0'
      return
    endif
    alphai= alpha

    f0= f
    fp= f0
    do iter=1,niter_linmin
      x1(1:ndim)= x0(1:ndim)
      x1(1:ndim)= x1(1:ndim) +alphai*d(1:ndim)
      call wrap_ranges(ndim,x1,xranges)
      call func(ndim,x1,fi,ftsti)
      if( myid.eq.0 .and. iprint.gt.2 ) write(6,'(a,i5,5es12.4)') &
           ' armijo: iter,fi,fi-f0,fi-fp,xigd*alphai,alphai=',&
           iter,fi,fi-fp,xigd*alphai,alphai
      if( fi-fp.le.xigd*alphai ) then
        f= fi
        alpha= alphai
        ftst= ftsti
        niter = iter
        return
      endif
      fp= fi
      pvalp= pval
      alphap= alphai
      alphai= alphai*armijo_tau
    enddo

    iflag= iflag +100
    niter= iter
    if( myid.eq.0 .and. iprint.gt.0 ) then
      print *,'WARNING: iter.gt.NITER_LINMIN in armijo_search.'
!!$      write(6,'(a,es13.5)') '   alphai   = ',alphai
!!$      write(6,'(a,es13.5)') '   xigd    = ',xigd
!!$      write(6,'(a,es13.5)') '   norm(g) = ',sqrt(sprod(ndim,g,g))
!!$      if( trim(cpena).eq.'lasso' .or. trim(cpena).eq.'glasso' .or. &
!!$           trim(cpena).eq.'ridge' ) then
!!$        write(6,'(a,es13.5)') '   pval    = ',pval
!!$      endif
    endif
    return

  end subroutine armijo_search
!=======================================================================
  subroutine backtrack(ndim,x0,xranges,d,f,ftst,alpha,iprint &
       ,iflag,myid,func,niter)
!
!  Simply move onestep towards current direction with max length
!  that can decreases function value.
!
    implicit none
    integer,intent(in):: ndim,iprint,myid
    integer,intent(inout):: iflag,niter
    real(8),intent(in):: x0(ndim),d(ndim),xranges(2,ndim)
    real(8),intent(inout):: f,alpha,ftst
    interface
      subroutine func(n,x,ftrn,ftst)
        integer,intent(in):: n
        real(8),intent(in):: x(n)
        real(8),intent(out):: ftrn,ftst
      end subroutine func
    end interface

!.....Precision
    real(8),parameter:: tiny = 1d-15
    integer:: iter,iterp
    real(8):: alphai,alphap,f0,fi,fp,ftsti,fpi,fti
    real(8),save,allocatable:: x1(:),gpena(:)
    logical,save:: l1st = .true.

    if( l1st ) then
      l1st = .false.
      if( myid.eq.0 .and. iprint.gt.1 ) &
           print *,'backtrack, alpha=',alpha
    endif
    if( .not.allocated(x1) ) allocate(x1(ndim),gpena(ndim))
    f0 = f
    fp = f0
    alphai = alpha
    do iter=1,niter_linmin
      x1(:) = x0(:)
      x1(1:ndim) = x1(1:ndim) +alphai*d(1:ndim)
      call wrap_ranges(ndim,x1,xranges)
      call func(ndim,x1,fi,ftsti)
      if( myid.eq.0 .and. iprint.gt.2 ) then
        print '(a,i8,5es12.4)','   iter,alphai,fi,fti,fi-f0,fi-fp = ' &
             ,iter,alphai,fi,fti,fi-f0,fi-fp
      endif
      if( fi.lt.f0 ) then
        f = fi
        alpha = alphai
        ftst = ftsti
        niter = iter
        return
      else  ! if fi > f0, decrease alpha
!!$        fp = min(fi,f0)
!!$        ftstp = ftsti
        fp = fi
        alphap = alphai
        alphai = alphai *fac_dec
        iterp = iter
        if( alphai.lt.tiny ) then
          if( myid.eq.0 .and. iprint.gt.0 ) then
            print *,'WARNING: alpha < tiny in backtrack,'
!!$            print *,'         The search direction would be wrong.'
!!$            print *,'   iter,alphai,fi=',iter,alphai,fi
          endif
          iflag = iflag + 100
          niter = iter
          return
        endif
      endif
    enddo

    iflag = iflag + 100
    niter = iter
    if( myid.eq.0 .and. iprint.gt.0 ) then
      print *, 'WARNING: iter exceeds NITER_LINMIN in backtrack.'
!!$      print *, '         The search direction would be wrong.'
!!$      write(6,'(a,es13.5)') '   alphai = ',alphai
    endif
    return
  end subroutine backtrack
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
    real(8):: xad,sgn,val

    do i=1,ndim
      xad= x(i) +alpha*d(i)
      sgn= sign(1d0,xad)   ! sign(a,b) = |a| * (sign of b)
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
    use descriptor,only: glval,ngl,iglid
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

!!$    if( trim(cpena).ne.'lasso' .and. trim(cpena).ne.'glasso' ) then
!!$      if(myid.eq.0) then
!!$        print *,'>>> fs works only with lasso or glasso.'
!!$      endif
!!$      iflag= iflag +100
!!$      return
!!$    endif

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
  subroutine gfs(ndim,x,f,ftst,g,d,xtol,gtol,ftol,xranges,maxiter &
       ,iprint,iflag,myid,func,grad,cfmethod,niter_eval &
       ,sub_eval)
!
!  Grouped Forward Stepwise (grouped FS) regression
!
    use descriptor,only: ngl,mskgfs,msktmp,iglid
    use variables,only: gsfcorr
    use random
    implicit none
    integer,intent(in):: ndim,iprint,myid,niter_eval
    integer,intent(inout):: iflag,maxiter
    real(8),intent(in):: xtol,gtol,ftol,xranges(2,ndim)
    real(8),intent(inout):: f,ftst,x(ndim),g(ndim),d(ndim)
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
    end interface

    integer:: iter,i,ig,jg,j,igmm,itergfs,niter,inc,jnc&
         ,nfailinmin
    real(8):: alpha,gnorm,pval,fp,f0,gmm &
         ,ftstp,fpgfs,frefb
    real(8),allocatable,save:: xt(:),xtb(:),gmaxgl(:),u(:),gmaxgl0(:)
    real(8),save,allocatable:: gg(:,:),y(:),gp(:),rg(:) &
         ,ggy(:),ygg(:),s(:),g0(:),gpena(:)  !,aa(:,:),cc(:,:),v(:)
    logical,allocatable,save:: lexclude(:)
    integer:: nmsks,nftol,nbases,nvar,nrefresh
    real(8):: ynorm,tmp1,b,sy,syi  !,svy,svyi
    character(len=128):: cnum
    real(8),parameter:: sgm = 1.0d-1

    if( trim(cpena).eq.'lasso' .and. trim(cpena).eq.'glasso' ) then
      if(myid.eq.0) then
        print *,'>>> gfs does not work with lasso or glasso.'
        print *,'>>> so gfs neglects lasso and glasso.'
      endif
      cpena = 'none'
    endif

    if( .not. allocated(gsfcorr) ) then
      if( myid.eq.0 ) then
        print *,'ERROR@gfs: gsfcorr not allocated!'
        print *,'  normalize_input should be norm/variance.'
      endif
      return
    endif

    if( .not.allocated(xt) ) allocate(xt(ndim),xtb(ndim),u(ndim),rg(ngl),g0(ndim))
    if( .not.allocated(gg) ) allocate(gg(ndim,ndim) &
         ,y(ndim),gp(ndim),ggy(ndim),ygg(ndim) &
         ,s(ndim),gpena(ndim))  !,v(ndim),aa(ndim,ndim),cc(ndim,ndim)
    if( .not. allocated(gmaxgl) ) then
      allocate(gmaxgl(ngl),gmaxgl0(ngl),lexclude(ngl))
      lexclude(:) = .false.
    endif
    if( .not.allocated(mskgfs) ) then
      allocate(mskgfs(ngl),msktmp(ngl))
    endif

    if( trim(cread_fsmask).eq.'' ) then
      mskgfs(1:ngl)= 1
    else
      call read_fsmask(cread_fsmask)
    endif
    msktmp(1:ngl)= mskgfs(1:ngl)

    xt(1:ndim)= x(1:ndim)
    do i=1,ndim
      ig= iglid(i)
      if( ig.lt.0 ) cycle
      if( mskgfs(ig).ne.0 ) xt(i)= 0d0
    enddo
    x(1:ndim)= xt(1:ndim)

    nmsks= 0
    do ig=1,ngl
      if( mskgfs(ig).eq.0 ) cycle
      nmsks= nmsks +1
    enddo
    nbases = ngl -nmsks

!.....cfsmode==df0corr, loss-func decrease of each descriptor
    if( index(cfsmode,'df0').ne.0 ) then
      call func(ndim,xt,f,ftst)
      call grad(ndim,xt,g)
      fp = f
      gp(:) = g(:)
!!$      print *,'initial f = ',fp
      do ig=1,ngl
        u(:) = -gp(:)
        g(:) = gp(:)
        do i=1,ndim
          if( ig.eq.i ) cycle
          u(i) = 0d0
          g(i) = 0d0
        enddo
        f = fp
!!$      alpha = 1d0
!!$      call backtrack(ndim,xt,xranges,u,f,ftst,alpha,iprint &
!!$           ,iflag,myid,func,niter)
        alpha = min(max(alpha,xtol*10d0)*2d0, 1d0)
        call armijo_search(ndim,xt,xranges,u,f,ftst,g,alpha,iprint &
             ,iflag,myid,func,niter)
!!$        print *,'ig,f,f-fp=',ig,f,f-fp
        gmaxgl0(ig) = -(f-fp)
      enddo
    else if( index(cfsmode,'grad0').ne.0 ) then
      call func(ndim,xt,f,ftst)
      call grad(ndim,xt,g)
      f0 = f
      g0(:) = g(:)
    endif

    if( myid.eq.0 .and. iprint.gt.0 ) then
      print '(a,i6,2es15.6e3,i6)',' iter,ftrn,ftst,nbases=',0,f,ftst,nbases
    endif
    call sub_eval(0)

!.....Start gFS loop here
    fpgfs = f
    ftstp = ftst
    do iter=1,maxiter
      if( index(cfsmode,'grad0').ne.0 ) then
        g(:) = g0(:)
      else if( index(cfsmode,'grad').ne.0 ) then
!.....First, calc of gradient needs to be done with no masks
!     because it is used to find another new basis
        msktmp(1:ngl)= mskgfs(1:ngl)
        mskgfs(1:ngl)= 0
        call func(ndim,xt,f,ftst)
        call grad(ndim,xt,g)
!.....Restore mask
        mskgfs(1:ngl)= msktmp(1:ngl)
      endif
      gnorm= sqrt(sprod(ndim,g,g))

      if( nmsks.eq.0 ) then
        if( myid.eq.0 ) then
          print *,'ngl=',ngl
          print *,'nmsks is already 0, and finsh the WHILE loop.'
        endif
        exit
      endif
!.....Find bases with the largest gradient
      gmaxgl(1:ngl)= 0d0
      do i=1,ndim
        ig= iglid(i)
        if( ig.le.0 ) cycle
        gmaxgl(ig)= gmaxgl(ig) +abs(g(i))
      enddo
      if( trim(cfsmode).eq.'grad0corr' ) then
        if( iter.eq.1 ) then
          gmaxgl0(:) = gmaxgl(:)
        else
          gmaxgl(:) = gmaxgl0(:)
        endif
      else if( trim(cfsmode).eq.'df0corr' ) then
        gmaxgl(:) = gmaxgl0(:)
      endif
      if( index(cfsmode,'corr').ne.0  .and. nbases.gt.0 ) then
        rg(:) = 0d0
        do ig=1,ngl
          if( mskgfs(ig).ne.0 ) then
            do jg=1,ngl
              if( jg.eq.ig ) cycle
              if( mskgfs(jg).ne.0 ) cycle
!.....Adopt maximum correlation as rg
              rg(ig) = max(rg(ig),abs(gsfcorr(ig,jg)))
            enddo
          endif
          if( myid.eq.0 .and. iprint.gt.1 ) then
            print '(a,i5,3es12.3e3)','  ig,rg,gmax,gmax*(1-rg)=',ig,rg(ig),gmaxgl(ig) &
                 ,gmaxgl(ig)*(1d0 -rg(ig))
          endif
!.....Scale gmaxgl w.r.t. rg
          gmaxgl(ig) = gmaxgl(ig) *(1d0- rg(ig))
        enddo
      else
        if( myid.eq.0 .and. iprint.gt.1 ) then
          do ig=1,ngl
            print '(a,i5,es12.4)','  ig,gmaxgl=',ig,gmaxgl(ig)
          enddo
        endif
      endif
      gmm= 0d0
      igmm= 0
      do ig=1,ngl
!.....Do not take mskgfs==2 into account !
!!$        print *,'myid,ig,mskgfs,gmaxgl=',myid,ig,mskgfs(ig),gmaxgl(ig)
        if( mskgfs(ig).eq.1 .and. gmaxgl(ig).gt.gmm .and. &
             .not. lexclude(ig) ) then
          gmm= gmaxgl(ig)
          igmm= ig
        endif
      enddo
      if( igmm.eq.0 ) then
        if( myid.eq.0 .and. iprint.gt.0 ) then
          print *,'igmm.eq.0 !!!'
          print *,'Nothing to do here, and going out from FS.'
        endif
        return
      endif
!.....remove mask of bases with large variations
      mskgfs(igmm)= 0
      if( myid.eq.0 ) then
        if( iprint.gt.0 ) then
          print '(a,2i5,es12.4)',' iter,igmm,gmm= ',iter,igmm,gmm
          call flush(6)
        endif
      endif
      nmsks= 0
      do ig=1,ngl
        if( mskgfs(ig).eq.0 ) cycle
        nmsks= nmsks +1
      enddo
      nbases= ngl -nmsks
      write(cnum,'(i0)') nbases
      call write_fsmask('out.fsmask.'//trim(cnum))

      nrefresh = 0
      frefb = 1d+30
      xtb(:) = xt(:)
10    continue
      nrefresh = nrefresh + 1
      x(1:ndim)= xt(1:ndim)
!.....Reset xt before going into the BFGS.
!     Because it can easily get stuck at the local minimum by starting BFGS
!     from the minimum of the previous BFGS even if another variable is added.
      if( trim(cfs_xrefresh).eq.'zero' ) then
!.....Reset all the xt connected to bases to 0.0 
        do i=1,ndim
          ig = iglid(i)
          if( ig.le.0 ) then
            xt(i) = sgm *polarbm()
          else
            xt(i)= 0d0
          endif
!!$          print *,'i,ig,xt=',i,ig,xt(i)
        enddo
      else if( trim(cfs_xrefresh).eq.'random' ) then
!.....Reset xt randomly
        do i=1,ndim
          xt(i) = sgm *polarbm()
          ig = iglid(i)
          if( ig.le.0 ) cycle
          if( mskgfs(ig).ne.0 ) xt(i) = 0d0
        enddo
      else
!.....Do nothing, use current variable values.
      endif
      call func(ndim,xt,f,ftst)
      call grad(ndim,xt,g)
!!$!.....Penalty
!!$      if( trim(cpena).eq.'ridge' ) then
!!$        call penalty(cpena,ndim,pval,gpena,xt)
!!$        f = f +pval
!!$        g(1:ndim) = g(1:ndim) +gpena(1:ndim)
!!$      endif
!.....preparation for BFGS
      gg(1:ndim,1:ndim)= 0d0
      do i=1,ndim
        gg(i,i)= 1d0
      enddo
!.....Mask some g's
      do i=1,ndim
        ig= iglid(i)
        if( ig.le.0 ) cycle
        if( mskgfs(ig).ne.0 ) g(i)= 0d0
      enddo

      gnorm= sqrt(sprod(ndim,g,g))
      nftol= 0
      iflag= 0
      nfailinmin = 0
!.....BFGS loop begins
      do itergfs=1,ninnergfs
        
        u(1:ndim)= 0d0
        inc = 0
        do i=1,ndim
          ig = iglid(i)
          if( ig.gt.0 .and. mskgfs(ig).ne.0 ) cycle
          inc = inc + 1
          jnc = 0
          do j=1,ndim
            jg = iglid(j)
            if( jg.gt.0 .and. mskgfs(jg).ne.0 ) cycle
            jnc = jnc + 1
            u(i)= u(i) -gg(inc,jnc)*g(j)
          enddo
        enddo
!.....mask u
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
!.....line minimization
        if( trim(clinmin).eq.'armijo' ) then
          alpha = min(max(alpha,xtol*2d0)*2d0, 1d0)
          call armijo_search(ndim,xt,xranges,u,f,ftst,g,alpha,iprint &
               ,iflag,myid,func,niter)
!.....if something wrong with armijo search, try opposite direction
          if( iflag/100.ne.0 ) then
            alpha= -1d0
            if(myid.eq.0) print *,'trying opposite direction...'
            call armijo_search(ndim,xt,xranges,u,f,ftst,g,alpha,iprint &
                 ,iflag,myid,func,niter)
          endif
        else   ! backtrack (default)
          alpha = min(max(alpha,xtol*2d0)*2d0, 1d0)
          call backtrack(ndim,xt,xranges,u,f,ftst,alpha,iprint &
               ,iflag,myid,func,niter)
        endif
!.....get out of bfgs loop
!!$        print '(a,3i4,2es12.4)','iter,itergfs,iflag,ftrn,fpgfs=' &
!!$             ,iter,itergfs,iflag,f,fpgfs
        if( itergfs.lt.ninnergfs .and. iflag/100.ne.0 ) then
          if( itergfs.eq.1 ) then
            if( myid.eq.0 ) then
              print *,"Line minimizaitonat failed 1st step of FS,"&
                   //" but keep going forward..."
            endif
!!$!.....Set mask as 2, which means this basis will be not included
!!$!     and not taken into consideration anymore.
!!$            mskgfs(igmm)= 2
          else ! itergsf.ne.1
            nfailinmin = nfailinmin +1
            if( nfailinmin.eq.1 ) then
              if( myid.eq.0 .and. iprint.gt.1 ) then
                print *,'>>> Initialize gg because alpha was not found.'
              endif
              gg(:,:) = 0d0
              do i=1,ndim
                gg(i,i) = 1d0
              enddo
              f = fp
              iflag = iflag -100*(iflag/100)
              alpha = 1d0
              cycle
            else ! nfailnmin > 1
              if( myid.eq.0 ) then
                print *,'>>> Line minimization failed twice consecutively.'
              endif
              if( f.gt.fpgfs .and. trim(cfs_xrefresh).eq.'random' &
                   .and. nrefresh.le.maxfsrefresh ) then
                if( myid.eq.0 .and. iprint.gt.0 ) then
                  print *,'>>> But refresh x and go back to BFGS, ' &
                       //'because current f > fpgfs.'
!!$                  print *,'f,fpgfs=',f,fpgfs
                endif
                goto 10
              else if( f.gt.fpgfs ) then
!.....If f is greater than the previous value,
!.....selected symmetry function is probably not a good one, so exclude it
!.....in the following gFS steps.
                lexclude(igmm) = .true.
                mskgfs(igmm) = 1
                if( myid.eq.0 .and. iprint.gt.0 ) &
                     print *,'>>> Exclude the basis since f > fpgfs, igmm=',igmm
              endif
              exit
            endif
          endif
        else ! iflag/100.eq.0
          nfailinmin = 0
        endif
        xt(1:ndim)= xt(1:ndim) +alpha*u(1:ndim)
        call wrap_ranges(ndim,xt,xranges)
        call grad(ndim,xt,g)
        do i=1,ndim
          ig= iglid(i)
          if( ig.le.0 ) cycle
          if( mskgfs(ig).ne.0 ) then
            g(i)= 0d0
          endif
        enddo

        gnorm= sqrt(sprod(ndim,g,g))
        if( myid.eq.0 ) then
          if( iprint.gt.1 ) then
            write(6,'(a,i8,2es13.5)') ' itergfs,f,gnorm=',itergfs,f,gnorm
            call flush(6)
          endif
        endif

!.....check convergence
        if( gnorm.lt.gtol ) then
          if( myid.eq.0 .and. iprint.gt.0 ) then
            print '(a,2es13.5)',' >>> QN in gFS converged; gnorm,gtol= ' &
                 ,gnorm,gtol
          endif
          if( f.gt.fpgfs .and. trim(cfs_xrefresh).eq.'random' &
               .and. nrefresh.le.maxfsrefresh ) then
            if( myid.eq.0 .and. iprint.gt.0 ) then
              print *,'>>> But refresh x and go back to BFGS, because current f > fpgfs.'
!!$              print *,'f,fpgfs=',f,fpgfs
            endif
            goto 10
          else if( f.gt.fpgfs ) then
!.....If f is greater than the previous value,
!.....selected symmetry function is probably not a good one, so exclude it
!.....in the following gFS steps.
            lexclude(igmm) = .true.
            mskgfs(igmm) = 1
            if( myid.eq.0 .and. iprint.gt.0 ) &
                 print *,'>>> Exclude the basis since f > fpgfs, igmm=',igmm
          endif
          iflag= iflag +2
          exit
        else if( abs(f-fp)/abs(fp).lt.ftol) then
          nftol= nftol +1
          if( nftol.ge.numtol ) then
            if( myid.eq.0 .and. iprint.gt.0 ) then
              print '(a,i2,a)',' >>> QN in gFS is probably converged because ' // &
                   'of ftol ',numtol,' times consecutively.'
            endif
            if( f.gt.fpgfs .and. trim(cfs_xrefresh).eq.'random' &
                 .and. nrefresh.le.maxfsrefresh ) then
              if( myid.eq.0 .and. iprint.gt.0 ) then
                print *,'>>> But refresh x and go back to BFGS, because current f > fpgfs.'
!!$                print *,'f,fpgfs=',f,fpgfs
              endif
              goto 10
            else if( f.gt.fpgfs ) then
!.....If f is greater than the previous value,
!.....selected symmetry function is probably not a good one, so exclude it
!.....in the following gFS steps.
              lexclude(igmm) = .true.
              mskgfs(igmm) = 1
              if( myid.eq.0 .and. iprint.gt.0 ) &
                   print *,'>>> Exclude the basis since f > fpgfs, igmm=',igmm
            endif
            iflag= iflag +3
            exit
          else if( itergfs.lt.ninnergfs ) then
            if( myid.eq.0 .and. iprint.gt.1 ) then
              print *,'>>> gg initialized because |f-fp|/|fp|<ftol '
            endif
            gg(1:ndim,1:ndim)= 0d0
            do i=1,ndim
              gg(i,i)= 1d0
            enddo
            cycle
          endif
        endif
        if( itergfs.ge.ninnergfs ) then
          if( myid.eq.0 .and. iprint.gt.0 ) then
            print *,'>>> itergfs exceeds fs_num_inner in QN in gFS.'
          endif
          if( f.gt.fpgfs .and. trim(cfs_xrefresh).eq.'random' &
               .and. nrefresh.le.maxfsrefresh ) then
            if( myid.eq.0 .and. iprint.gt.0 ) then
              print *,'>>> But refresh x and go back to BFGS, because current f > fpgfs.'
!!$              print *,'f,fpgfs=',f,fpgfs
            endif
            goto 10
          else if( f.gt.fpgfs ) then
!.....If f is greater than the previous value,
!.....selected symmetry function is probably not a good one, so exclude it
!.....in the following gFS steps.
            lexclude(igmm) = .true.
            mskgfs(igmm) = 1
            if( myid.eq.0 .and. iprint.gt.0 ) &
                 print *,'>>> Exclude the basis since f > fpgfs, igmm=',igmm
          endif
          exit
        endif
        nftol = 0

!.....BFGS treatment with reduced number of variables
        inc = 0
        s(:) = 0d0
        y(:) = 0d0
        do i=1,ndim
          ig = iglid(i)
          if( ig.gt.0 .and. mskgfs(ig).ne.0 ) cycle
          inc = inc + 1
          s(inc) = alpha *u(i)
          y(inc) = g(i) -gp(i)
        enddo
        nvar = inc
!!$        s(1:ndim)= alpha *u(1:ndim)
!!$        y(1:ndim)= g(1:ndim) -gp(1:ndim)
        ynorm= sprod(ndim,y,y)
        if( ynorm.lt.1d-14 ) then
          if( myid.eq.0 .and. iprint.gt.1 ) then
            print *,'>>> Initialize gg because y*y < 1d-14'
          endif
          gg(1:ndim,1:ndim)= 0d0
          do i=1,ndim
            gg(i,i)= 1d0
          enddo
          cycle
        endif

!.....update matrix gg
!!$        sy= sprod(ndim,s,y)
        sy= sprod(nvar,s,y)
        syi= 1d0/sy
!!$        do i=1,ndim
        do i=1,nvar
          tmp1= 0d0
!!$          tmp2= 0d0
!!$          do j=1,ndim
          do j=1,nvar
            tmp1= tmp1 +gg(j,i)*y(j)
          enddo
          ggy(i)= tmp1 *syi
        enddo
        b= 1d0
!!$        do i=1,ndim
        do i=1,nvar
          b=b +y(i)*ggy(i)
        enddo
        b= b*syi
!.....without temporary matrix aa
!!$        do j=1,ndim
!!$          do i=1,ndim
        do j=1,nvar
          do i=1,nvar
            gg(i,j)=gg(i,j) +s(j)*s(i)*b &
                 -(s(i)*ggy(j) +ggy(i)*s(j))
          enddo
        enddo
!!$        print *,'End of inner BFGS'
      enddo  ! End of inner BFGS

!!$      print '(a,i5,2es12.4)','iter,ftrn,fpgfs=',iter,f,fpgfs
      if( f.gt.fpgfs .and. .not.lexclude(igmm) ) then
        print *,'SOMETHING WRONG!!!'
        print *,'itergfs,ninnergfs=',itergfs,ninnergfs
      endif
      if( .not. lexclude(igmm) ) then
        fpgfs = f
        ftstp = ftst
        x(1:ndim)= xt(1:ndim)
        if( myid.eq.0 .and. iprint.gt.0 ) then
          print '(a,i6,2es15.6e3,i6)',' iter,ftrn,ftst,nbases=' &
               ,iter,f,ftst,nbases
        endif
        call sub_eval(iter)
      endif
    enddo

    if( myid.eq.0 .and. iprint.gt.0 ) &
         print *,'iter exceeds maxiter in gFS.'
    iflag= iflag +10
    x(1:ndim)= xt(1:ndim)
    return

  end subroutine gfs
!=======================================================================
  subroutine read_fsmask(cfname)
!
!  Read mask for gfs from file
!
    use descriptor,only: ngl,mskgfs
    use parallel
    implicit none
    character(len=*),intent(in):: cfname

    integer,parameter:: iomask = 40
    integer:: ndata,i,ierrcode,itmp,ig

    ierrcode = 0
    if( myid.eq.0 ) then
      print *,'Read mskgfs from '//trim(cfname)
      open(iomask,file=trim(cfname),status='old')
      read(iomask,*) ndata
      if( ndata.ne.ngl ) then
        ierrcode = 1
        close(iomask)
        goto 999
      endif
      do i=1,ndata
        read(iomask,*,err=99) ig, mskgfs(ig)
      enddo
      close(iomask)
    endif

999 continue
    itmp = 0
    call mpi_allreduce(ierrcode,itmp,1,mpi_integer,mpi_max &
         ,mpi_world,ierr)
    if( itmp.gt.0 ) then
      if( myid.eq.0 ) then
        if( itmp.eq.1 ) print *,'ERROR@read_fsmask: ndata.ne.ngl !'
        if( itmp.eq.2 ) print *,'ERROR@read_fsmask: something wrong with data !'
      endif
      call mpi_finalize(ierr)
    endif
    call mpi_bcast(mskgfs,ngl,mpi_integer,0,mpi_world,ierr)
    return

99  ierrcode = 2
    close(iomask)
    goto 999
  end subroutine read_fsmask
!=======================================================================
  subroutine write_fsmask(cfname)
    use parallel
    use descriptor,only: ngl,mskgfs
    implicit none 
    character(len=*),intent(in):: cfname

    integer,parameter:: iomask = 41
    integer:: i

    if( myid.eq.0 ) then
      open(iomask,file=trim(cfname),status='replace')
      write(iomask,'(i6)') ngl
      do i=1,ngl
        write(iomask,'(i7,i3)') i, mskgfs(i)
      enddo
      close(iomask)
    endif
    call mpi_barrier(mpi_world,ierr)
    return
  end subroutine write_fsmask
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
!!$!=======================================================================
!!$  subroutine penalty(cpena,ndim,fp,gp,x)
!!$!
!!$! Calculate penalty term and its derivative.
!!$! lasso and ridge are available.
!!$!
!!$    use descriptor,only: glval,ngl,iglid
!!$    implicit none
!!$    character(len=*),intent(in):: cpena
!!$    integer,intent(in):: ndim
!!$    real(8),intent(in):: x(ndim)
!!$    real(8),intent(out):: fp,gp(ndim)
!!$
!!$    integer:: i,ig
!!$    real(8):: absx,sgnx
!!$    real(8),parameter:: xtiny  = 1d-14
!!$
!!$    fp= 0d0
!!$    gp(1:ndim)= 0d0
!!$    if( trim(cpena).eq.'lasso' ) then
!!$      do i=1,ndim
!!$        absx= abs(x(i))
!!$        fp= fp +pwgt*absx
!!$        sgnx= sign(1d0,x(i))
!!$        if( absx.gt.xtiny ) gp(i)= pwgt*sgnx
!!$      enddo
!!$    else if( trim(cpena).eq.'glasso' ) then
!!$      glval(0:ngl)= 0d0
!!$      do i=1,ndim
!!$        ig= iglid(i)
!!$        if(ig.gt.0) glval(ig)= glval(ig) +x(i)*x(i)
!!$      enddo
!!$      glval(0)= 1d0
!!$      do ig=1,ngl
!!$        glval(ig)= sqrt(glval(ig))
!!$        fp= fp +pwgt*glval(ig)
!!$      enddo
!!$      do i=1,ndim
!!$        ig= iglid(i)
!!$        if( ig.eq.0 ) then ! i is not in a group
!!$          absx= abs(x(i))
!!$          sgnx= sign(1d0,x(i))
!!$          if( absx.gt.xtiny ) gp(i)= pwgt*sgnx
!!$          fp= fp +pwgt*absx
!!$        else if( ig.gt.0 ) then ! i is in a group
!!$          if( glval(ig).gt.xtiny) gp(i)= pwgt*x(i)/glval(ig)
!!$        endif
!!$      enddo
!!$    else if( trim(cpena).eq.'ridge' ) then
!!$      do i=1,ndim
!!$        fp= fp +pwgt*x(i)*x(i)
!!$        gp(i)= 2d0*pwgt*x(i)
!!$      enddo
!!$    endif
!!$
!!$    return
!!$  end subroutine penalty
!=======================================================================
  
end module
!-----------------------------------------------------------------------
! Local Variables:
! compile-command: "make fitpot"
! End:
