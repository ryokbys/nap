module linreg
!.....parameter file name
  character(128),parameter:: cpfname= 'in.params.linreg'
  character(128),parameter:: ccfname='in.const.linreg'
!.....parameters
  integer:: ncoeff
  real(8),allocatable:: coeff(:)
!.....constants
  integer:: nelem
  integer,allocatable:: itype(:)
  real(8),allocatable:: cnst(:,:)
!.....function types and num of constatns for types
  integer,parameter:: max_ncnst= 10
  integer,parameter:: ncnst_type(1:3)= &
       (/ 2, &  ! Gaussian
          1, &  ! cosine
          5 /)  ! polynomial
!.....max exponent of the basis function
  integer:: max_nexp
  
contains
  subroutine force_linreg(namax,natm,tag,ra,nnmax,aa,strs,h,hi,tcom &
       ,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_world,myid,epi,epot,nismax,acon,avol)
!-----------------------------------------------------------------------
!  Parallel implementation of linear regression potential for pmd
!    - 2014.06.11 by R.K.
!      1st implementation
!-----------------------------------------------------------------------
    implicit none
    include "mpif.h"
    include "./params_unit.h"
    integer,intent(in):: namax,natm,nnmax,nismax
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_world,myid,lspr(0:nnmax,namax)
    real(8),intent(in):: ra(3,namax),tag(namax),acon(nismax) &
         ,h(3,3),hi(3,3),sv(3,6)
    real(8),intent(inout):: tcom,avol,rc
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)

!.....local
    integer:: i,j,k,l,m,n,ixyz,jxyz,is,js,ks,ierr,nbl,ia,nexp,ielem
    real(8):: rcin,b_na,at(3),epotl
    real(8),allocatable:: fat(:,:)
!.....1st call
    logical,save:: l1st=.true.

    if( l1st ) then
!.....read in.params.linreg
      call read_params(myid,mpi_world,rcin)
!.....reset rc
      if( myid.eq.0 ) then
        write(6,'(a,f10.5,a,f10.5)') &
             ' Cutoff radius rc may have been changed from '&
             ,rc,' to ',rcin
      endif
      rc= rcin
      allocate(fat(3,namax))
      l1st= .false.
    endif

    epotl= 0d0
    epi(1:natm+nb)= 0d0
    strs(1:3,1:3,1:namax)= 0d0
    aa(1:3,1:namax)= 0d0

    do ia=1,natm
      do nexp=1,max_nexp
        do ielem=1,nelem
          b_na= 0d0
          fat(1:3,1:natm+nb)= 0d0
          call bfunc(ia,natm,namax,nnmax,ra,lspr,h,tag,fat,rc &
               ,ielem,nexp,b_na)
          epotl=epotl +b_na
          epi(ia)= epi(ia) +b_na
          aa(1:3,1:natm+nb)= aa(1:3,1:natm+nb) +fat(1:3,1:natm+nb)
        enddo
      enddo
    enddo

    call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,lsrc,myparity &
         ,nn,mpi_world,aa,3)
!-----reduced force
    do i=1,natm
      at(1:3)= aa(1:3,i)
      aa(1:3,i)= hi(1:3,1)*at(1) +hi(1:3,2)*at(2) +hi(1:3,3)*at(3)
    enddo
!-----multiply 0.5d0*dt**2/am(i)
    do i=1,natm
      is= int(tag(i))
      aa(1:3,i)= acon(is)*aa(1:3,i)
    enddo

!-----gather epot
    epot= 0d0
    call mpi_allreduce(epotl,epot,1,mpi_double_precision &
         ,mpi_sum,mpi_world,ierr)
    return
  end subroutine force_linreg
!=======================================================================
  function gauss(r,b,a)
    implicit none
    real(8),intent(in):: a,b,r
    real(8):: gauss

    gauss= r**b *exp(-a*r*r)
    return
  end function gauss
!=======================================================================
  function dgauss(r,b,a)
    implicit none
    real(8),intent(in):: a,b,r
    real(8):: dgauss

    dgauss= (b-2d0*a*r*r) *r**(b-1) *exp(-a*r*r)
    return
  end function dgauss
!=======================================================================
  function poly(r,a0,a1,a2,a3,a4)
    implicit none
    real(8),intent(in):: r,a0,a1,a2,a3,a4
    real(8):: poly,r2

    r2= r*r
    poly= a0 +a1*r +a2*r2 +a3*r*r2 +a4*r2*r2
    return
  end function poly
!=======================================================================
  function dpoly(r,a0,a1,a2,a3,a4)
    implicit none
    real(8),intent(in):: r,a0,a1,a2,a3,a4
    real(8):: dpoly,r2

    r2=r*r
    dpoly= a1 +2d0*a2*r +3d0*a3*r2 +4d0*a4*r*r2
    return
  end function dpoly
!=======================================================================
  function fc(r,rc)
    implicit none
    real(8),intent(in):: r,rc
    real(8):: fc
    real(8),parameter:: pi= 3.14159265358979d0

    fc= 0.5d0 *(cos(r/rc*pi)+1d0)
    return
  end function fc
!=======================================================================
  function dfc(r,rc)
    implicit none
    real(8),intent(in):: r,rc
    real(8):: dfc
    real(8),parameter:: pi= 3.14159265358979d0

    dfc= -pi/2/rc *sin(r/rc*pi)
    return
  end function dfc
!=======================================================================
  subroutine bfunc(ia,natm,namax,nnmax,ra,lspr,h,tag,fat,rc &
       ,ielem,nexp,b_na)
!
!  basis function in the linear regression potetnial
!
    implicit none
    integer,intent(in):: ia,natm,namax,nnmax,lspr(0:nnmax,natm) &
         ,ielem,nexp
    real(8),intent(in):: ra(3,namax),h(3,3),tag(namax),rc
    real(8),intent(out):: b_na,fat(3,namax)

    integer:: ja,n,nn
    real(8):: xi(3),xj(3),xij(3),rij(3),r,dirij(3),djrij(3),tmp

    xi(1:3)= ra(1:3,ia)
    nn= lspr(0,ia)
    do n=1,nn
      ja= lspr(n,ia)
      xj(1:3)= ra(1:3,ja)
      xij(1:3)= xj(1:3)-xi(1:3)
      rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) *h(1:3,3)*xij(3)
      r= sqrt(rij(1)**2 +rij(2)**2 +rij(3)**2)
      b_na= b_na +func(r,ielem)*fc(r,rc)
    enddo
    !.....b_na will be used in the following loop.
    !.....Therefore these two loops cannot be merged.
    do n=1,nn
      ja= lspr(n,ia)
      xj(1:3)= ra(1:3,ja)
      xij(1:3)= xj(1:3)-xi(1:3)
      rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) *h(1:3,3)*xij(3)
      r= sqrt(rij(1)**2 +rij(2)**2 +rij(3)**2)
      dirij(1:3)= -rij(1:3)/r
      djrij(1:3)= -dirij(1:3)
      tmp= dfunc(r,ielem)*fc(r,rc) +func(r,ielem)*dfc(r,rc)
      fat(1:3,ia)= fat(1:3,ia) -dirij(1:3)*tmp *nexp*b_na**(nexp-1)
      fat(1:3,ja)= fat(1:3,ja) -djrij(1:3)*tmp *nexp*b_na**(nexp-1)
    enddo
    b_na= b_na**nexp
    
  end subroutine bfunc
!=======================================================================
  function func(rij,ielem)
!
!  Calculate selected basis function specified by ielem.
!
!  ***If you add a function, you have to change ncnst_type parameter
!     in this module header.
!
    implicit none
    integer,intent(in):: ielem
    real(8),intent(in):: rij
    real(8):: func

    func= 0d0
    if( itype(ielem).eq.1 ) then ! Gaussian-type
      func= gauss(rij,cnst(ielem,1),cnst(ielem,2))

    elseif( itype(ielem).eq.2 ) then ! cosine-type
      func= cos(rij*cnst(ielem,1))

    elseif( itype(ielem).eq.3 ) then ! polynomial
      func= poly(rij,cnst(ielem,1),cnst(ielem,2) &
           ,cnst(ielem,3),cnst(ielem,4),cnst(ielem,5))

    endif
    
  end function func
!=======================================================================
  function dfunc(rij,ielem)
    implicit none
    integer,intent(in):: ielem
    real(8),intent(in):: rij
    real(8):: dfunc

    dfunc= 0d0
    if( itype(ielem).eq.1 ) then ! Gaussian-type
      dfunc= dgauss(rij,cnst(ielem,1),cnst(ielem,2))

    elseif( itype(ielem).eq.2 ) then ! cosine-type
      dfunc= -cnst(ielem,1)*sin(rij*cnst(ielem,1))

    elseif( itype(ielem).eq.3 ) then ! polynomial
      dfunc= dpoly(rij,cnst(ielem,1),cnst(ielem,2) &
           ,cnst(ielem,3),cnst(ielem,4),cnst(ielem,5))

    endif
    
  end function dfunc
!=======================================================================
  subroutine read_params(myid,mpi_world,rcin)
    implicit none
    include 'mpif.h'

    integer,intent(in):: myid,mpi_world
    real(8),intent(out):: rcin
    integer:: itmp,ierr,i,j
    logical:: lexist

!.....read parameters at the 1st call
    inquire(file=trim(cpfname),exist=lexist)
    if( .not. lexist ) then
      if( myid.eq.0 ) then
        write(6,'(a)') ' [Error] '//cpfname//' does not exist !!!.'
        write(6,'(a)') '   The linreg potential needs '//cpfname//'.'
      endif
      call mpi_finalize(ierr)
      stop
    endif
    open(50,file=trim(cpfname),status='old')
    read(50,*) ncoeff,rcin
    allocate(coeff(ncoeff))
    do i=1,ncoeff
      read(50,*) coeff(i)
    enddo
    close(50)
!.....read constants at the 1st call
    inquire(file=trim(ccfname),exist=lexist)
    if( .not. lexist ) then
      if( myid.eq.0 ) then
        write(6,'(a)') ' [Error] '//ccfname//' does not exist !!!.'
        write(6,'(a)') '   The linreg potential needs '//ccfname//'.'
      endif
      call mpi_finalize(ierr)
      stop
    endif
    open(51,file=trim(ccfname),status='old')
    read(51,*) nelem,max_nexp
    allocate(itype(nelem),cnst(nelem,max_ncnst))
    do i=1,nelem
      read(51,*) itype(i),(cnst(i,j),j=1,ncnst_type(itype(i)))
    enddo
    close(51)
!.....check whether the num of parameters is correct
    if( .not. ncoeff .eq. nelem*max_nexp) then
      write(6,'(a)') ' [Error] num of parameters is not correct !!!'
      write(6,'(a,i10)') ' ncoeff=',ncoeff
      write(6,'(a,2i10)') ' nelem,max_nexp=',nelem,max_nexp
      write(6,'(a,i10)') ' ncoeff should be ',nelem*max_nexp
      stop
    endif
    return
  end subroutine read_params

end module linreg
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
