module NN1
!.....parameter file name
  character(128),parameter:: cpfname= 'in.params.NN1'
  character(128),parameter:: ccfname='in.const.NN2'
!.....parameters
  integer:: nwgt1,nwgt2
  real(8),allocatable:: wgt1(:,:),wgt2(:)
!.....constants
  integer:: nsf,nhl1
  integer,allocatable:: itype(:)
  real(8),allocatable:: cnst(:,:)
!.....function types and num of constatns for types
  integer,parameter:: max_ncnst= 2
  integer,parameter:: ncnst_type(1:2)= &
       (/ 2, &  ! Gaussian
          1 /) ! angular
!.....max exponent of the basis function
  integer:: max_nexp
  
contains
  subroutine force_NN1(namax,natm,tag,ra,nnmax,aa,strs,h,hi,tcom &
       ,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_world,myid,epi,epot,nismax,acon,avol)
!-----------------------------------------------------------------------
!  Parallel implementation of neural-network potential of only one
!  hidden layer.
!    - 2014.07.01 by R.K.
!      start coding
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
    integer:: i,j,k,l,m,n,ixyz,jxyz,is,js,ks,ierr,nbl,ia,nexp,isf &
         ,icoeff,ihl1
    real(8):: rcin,b_na,at(3),epotl,wgt,hl1i,tmp2,tmp
    real(8),save,allocatable:: fat(:,:),gsf(:,:),dgsf(:,:,:)
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
      allocate(fat(3,namax),gsf(nsf,natm),dgsf(3,natm+nb,nsf))
      l1st= .false.
    endif

    epotl= 0d0
    epi(1:natm+nb)= 0d0
    strs(1:3,1:3,1:namax)= 0d0
    aa(1:3,1:namax)= 0d0
    fat(1:3,1:natm+nb)= 0d0

!.....first, calculate all the symmetry functions
    call eval_sf(nsf,namax,natm,nb,nnmax,h,tag,ra,lspr,gsf,dgsf,rc)

#ifdef __FITPOT__
    open(80,file='out.gsf-hl1',status='replace')
    write(80,'(2i10)') nsf,nhl1
    do ia=1,natm
      do isf=1,nsf
        write(80,'(2i8,es22.14)') ia,isf,gsf(isf,ia)
      enddo
    enddo
#endif

    do ia=1,natm
!.....second, sum up according to NN with one hidden layer
      epotl= epotl +wgt2(0)
      epi(ia)= epi(ia) + wgt2(0)
      tmp= wgt2(0)
      do ihl1=1,nhl1
        hl1i= wgt1(ihl1,0)
        do isf=1,nsf
          hl1i= hl1i +wgt1(ihl1,isf)*gsf(isf,ia)
        enddo
#ifdef __FITPOT__
        write(80,'(2i8,es22.14)') ia,ihl1,hl1i
#endif
        tmp2= tmp +wgt2(0)*sigmoid(hl1i)
      enddo
      epotl=epotl +tmp2
      epi(ia)= epi(ia) +tmp2
      aa(1:3,1:natm+nb)= aa(1:3,1:natm+nb) +fat(1:3,1:natm+nb)
    enddo

#ifdef __FITPOT__
    close(80)
#endif

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
  end subroutine force_NN1
!=======================================================================
  subroutine eval_sf(nsf,namax,natm,nb,nnmax,h,tag,ra,lspr,gsf,dgsf,rc)
!
!  Evaluate symmetry functions and derivatives of them.
!
    implicit none
    integer,intent(in):: nsf,namax,natm,nb,nnmax,lspr(0:nnmax,namax)
    real(8),intent(in):: h(3,3),tag(namax),ra(3,namax),rc
    real(8),intent(out):: gsf(nsf,natm),dgsf(3,natm+nb,nsf)

    integer:: isf,ia,jj,ja,kk,ka
    real(8):: xi(3),xj(3),xij(3),rij(3),dij,fcij,eta,rs,texp,driji(3), &
         dfcij,drijj(3),dgdr,xk(3),xik(3),rik(3),dik,fcik,dfcik, &
         driki(3),drikk(3),almbd,spijk,cs,t1,t2,dgdij,dgdik,dgcs, &
         dcsdj(3),dcsdk(3),dcsdi(3)

    real(8),external:: sprod

    do isf=1,nsf
      if( itype(isf).eq.1 ) then ! Gaussian (2-body)
        do ia=1,natm
          xi(1:3)= ra(1:3,ia)
          do jj=1,lspr(0,ia)
            ja= lspr(jj,ia)
            if( ja.eq.ia ) cycle
            xj(1:3)= ra(1:3,ja)
            xij(1:3)= xj(1:3)-xi(1:3)
            rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
            dij= sqrt(rij(1)**2 +rij(2)**2 +rij(3)**2)
            if( dij.ge.rc ) cycle
            fcij= fc(dij,rc)
            eta= cnst(1,isf)
            rs=  cnst(2,isf)
            !.....function value
            texp= exp(-eta*(dij-rs))
            gsf(isf,ia)= gsf(isf,ia) +texp*fcij
            !.....derivative
            driji(1:3)= -rij(1:3)/dij
            drijj(1:3)= -driji(1:3)
            dgdr= -2d0*eta*(dij-rs)*texp*fcij +texp*dfc(dij,rc)
            dgsf(1:3,ia,isf)= dgsf(1:3,ia,isf) +driji(1:3)*dgdr
            dgsf(1:3,ja,isf)= dgsf(1:3,ja,isf) +drijj(1:3)*dgdr
          enddo
        enddo

      else if( itype(isf).eq.2 ) then ! angular (3-body)
        do ia=1,natm
          xi(1:3)= ra(1:3,ia)
          do jj=1,lspr(0,ia)
            ja= lspr(jj,ia)
            if( ja.eq.ia ) cycle
            xj(1:3)= ra(1:3,ja)
            xij(1:3)= xj(1:3)-xi(1:3)
            rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
            dij= sqrt(rij(1)**2 +rij(2)**2 +rij(3)**2)
            if( dij.ge.rc ) cycle
            fcij= fc(dij,rc)
            dfcij= dfc(dij,rc)
            driji(1:3)= -rij(1:3)/dij
            drijj(1:3)= -driji(1:3)
            do kk=1,lspr(0,ia)
              ka= lspr(kk,ia)
              if( ka.eq.ia .or. ka.le.ja ) cycle
              xk(1:3)= ra(1:3,ka)
              xik(1:3)= xk(1:3)-xi(1:3)
              rik(1:3)= h(1:3,1)*xik(1) +h(1:3,2)*xik(2) +h(1:3,3)*xik(3)
              dik= sqrt(rik(1)**2 +rik(2)**2 +rik(3)**2)
              if( dik.ge.rc ) cycle
              fcik= fc(dik,rc)
              dfcik= dfc(dik,rc)
              driki(1:3)= -rik(1:3)/dik
              drikk(1:3)= -driki(1:3)
              almbd= cnst(1,isf)
              !.....function value
              spijk= sprod(rij,rik)
              cs= spijk/dij/dik
              t1= (almbd +cs)**2
              t2= (abs(almbd)+1d0)**2
              gsf(isf,ia)= gsf(isf,ia) +1d0/t2*t1*fcij*fcik
              !.....derivative
              dgdij= t1/t2 *dfcij*fcik
              dgdik= t1/t2 *fcij*dfcik
              dgcs= 2d0*(almbd+cs)/t2 *fcij*fcik
              dgsf(1:3,ia,isf)= dgsf(1:3,ia,isf) +dgdij*driji(1:3) &
                   +dgdik*driki(1:3)
              dgsf(1:3,ja,isf)= dgsf(1:3,ja,isf) +dgdij*drijj(1:3)
              dgsf(1:3,ka,isf)= dgsf(1:3,ka,isf) +dgdik*drikk(1:3)
              dcsdj(1:3)= rij(1:3)*(1d0-spijk/dij/dij)/dij/dik
              dcsdk(1:3)= rik(1:3)*(1d0-spijk/dik/dik)/dij/dik
              dcsdi(1:3)= -dcsdj(1:3) -dcsdk(1:3)
              dgsf(1:3,ia,isf)= dgsf(1:3,ia,isf) +dgcs*dcsdi(1:3)
              dgsf(1:3,ja,isf)= dgsf(1:3,ja,isf) +dgcs*dcsdj(1:3)
              dgsf(1:3,ka,isf)= dgsf(1:3,ka,isf) +dgcs*dcsdk(1:3)
            enddo
          enddo
        enddo
      endif
    enddo

  end subroutine eval_sf
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
  function sigmoid(x)
    implicit none
    real(8),intent(in):: x
    real(8):: sigmoid

    sigmoid= 1d0/(1d0 +exp(-x))
    return
  end function sigmoid
!=======================================================================
  function dsigmoid(x)
    implicit none
    real(8),intent(in):: x
    real(8):: dsigmoid,sx

    sx= sigmoid(x)
!    dsigmoid= -exp(-x)/(1d0+exp(-x))**2
    dsigmoid= sx*(1d0-sx)
    return
  end function dsigmoid
!=======================================================================
  subroutine read_params(myid,mpi_world,rcin)
    implicit none
    include 'mpif.h'

    integer,intent(in):: myid,mpi_world
    real(8),intent(out):: rcin
    integer:: itmp,ierr,i,j,nc,ncoeff,isf,ihl1
    logical:: lexist

!.....read constants at the 1st call
    inquire(file=trim(ccfname),exist=lexist)
    if( .not. lexist ) then
      if( myid.eq.0 ) then
        write(6,'(a)') ' [Error] '//ccfname//' does not exist !!!.'
        write(6,'(a)') '   The NN1 potential needs '//ccfname//'.'
      endif
      call mpi_finalize(ierr)
      stop
    endif
    open(51,file=trim(ccfname),status='old')
!.....num of symmetry functions, num of node in 1st hidden layer
    read(51,*) nsf,nhl1
    allocate(itype(nsf),cnst(max_ncnst,nsf))
    do i=1,nsf
      read(51,*) itype(i),(cnst(j,i),j=1,ncnst_type(itype(i)))
    enddo
    close(51)

!.....calc number of weights taking basis node into account
    nwgt1= (nsf+1)*nhl1
    nwgt2= (nhl1+1)

!.....read parameters at the 1st call
    inquire(file=trim(cpfname),exist=lexist)
    if( .not. lexist ) then
      if( myid.eq.0 ) then
        write(6,'(a)') ' [Error] '//cpfname//' does not exist !!!.'
        write(6,'(a)') '   The NN1 potential needs '//cpfname//'.'
      endif
      call mpi_finalize(ierr)
      stop
    endif
    open(50,file=trim(cpfname),status='old')
    read(50,*) ncoeff,rcin
!.....check whether the num of parameters is correct
    nc= nwgt1 +nwgt2
    if( ncoeff .ne. nc ) then
      write(6,'(a)') ' [Error] num of parameters is not correct !!!'
      write(6,'(a,i10)') ' ncoeff=',ncoeff
      write(6,'(a,2i10)') ' nsf,nhl1=',nsf,nhl1
      write(6,'(a,i10)') ' ncoeff should be ',nc
      stop
    endif
    allocate(wgt1(nhl1,0:nsf),wgt2(0:nhl1))
    do isf=0,nsf
      do ihl1=1,nhl1
        read(50,*) wgt1(ihl1,isf)
      enddo
    enddo
    do ihl1=0,nhl1
      read(50,*) wgt2(ihl1)
    enddo
    close(50)


    return
  end subroutine read_params

end module NN1
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
