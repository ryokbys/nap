module NN1
!.....parameter file name
  character(128),parameter:: cpfname= 'in.params.NN1'
  character(128),parameter:: ccfname='in.const.NN1'
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
!      started coding
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
    integer:: i,j,k,l,m,n,ixyz,jxyz,is,js,ks,ierr,nbl,ia,ja,nexp,isf &
         ,icoeff,ihl1,jj,jsf
    real(8):: rcin,b_na,at(3),epotl,wgt,hl1i,tmp2,tmp,tmp3(3)
    real(8),save,allocatable:: gsf(:,:),dgsf(:,:,:,:),hl1(:,:)
    real(8),allocatable:: aml(:,:,:,:),bml(:,:,:,:)
!.....1st call
    logical,save:: l1st=.true.

    if( l1st ) then
!.....read in.params.NN1
      call read_params(myid,mpi_world,rcin)
!.....reset rc
      if( myid.eq.0 ) then
        write(6,'(a,f10.5,a,f10.5)') &
             ' Cutoff radius rc may have been changed from '&
             ,rc,' to ',rcin
      endif
      rc= rcin
      allocate(gsf(0:nsf,namax),dgsf(3,nsf,0:nnmax,namax) &
           ,hl1(0:nhl1,natm))
      gsf(0:nsf,1:natm+nb)= 0d0
      l1st= .false.
    endif

!.....first, calculate all the symmetry functions
    call eval_sf(nsf,namax,natm,nb,nnmax,h,tag,ra,lspr,gsf,dgsf,rc)

#ifdef __FITPOT__
    open(80,file='out.NN1.gsf',status='replace')
    write(80,'(2i10)') nsf
    do ia=1,natm
      do isf=0,nsf
        write(80,'(2i8,es23.14e3)') ia,isf,gsf(isf,ia)
      enddo
    enddo
    close(80)
    open(81,file='out.NN1.hl1')
    write(81,'(2i10)') nhl1
#endif

!.....2nd, calculate the node values by summing contributions from
!.....  symmetry functions
    hl1(0:nhl1,1:natm)= 0d0
    do ia=1,natm
      hl1(0,ia)= 1d0
      do ihl1=1,nhl1
        tmp= 0d0
        do isf=0,nsf
          tmp= tmp +wgt1(ihl1,isf) *gsf(isf,ia)
        enddo
        hl1(ihl1,ia)= sigmoid(tmp)
      enddo
    enddo

#ifdef __FITPOT__
    do ia=1,natm
      do ihl1=0,nhl1
        write(81,'(2i8,es23.14e3)') ia,ihl1,hl1(ihl1,ia)
      enddo
    enddo
#endif

!.....then calculate the energy of atom by summing up the node values
    epotl= 0d0
    do ia=1,natm
      epi(ia)= 0d0
      do ihl1=0,nhl1
        hl1i= hl1(ihl1,ia)
        epi(ia)= epi(ia) +wgt2(ihl1)*hl1i
      enddo
      epotl=epotl +epi(ia)
    enddo

!.....sum up for forces
    strs(1:3,1:3,1:natm)= 0d0
    aa(1:3,1:natm+nb)= 0d0
    do ia=1,natm
      do ihl1=1,nhl1 ! there is no longer a bias node, 0
        hl1i= hl1(ihl1,ia)
        tmp= wgt2(ihl1)*hl1i*(1d0-hl1i)
        do jj=1,lspr(0,ia)
          ja= lspr(jj,ia)
          do isf=1,nsf ! there is no longer a bias node, 0
            aa(1:3,ja)=aa(1:3,ja) &
                 -tmp*wgt1(ihl1,isf)*dgsf(1:3,isf,jj,ia)
          enddo
        enddo
        !.....atom ia
        do isf= 1,nsf
          aa(1:3,ia)=aa(1:3,ia) &
               -tmp*wgt1(ihl1,isf)*dgsf(1:3,isf,0,ia)
        enddo
      enddo
    enddo

#ifdef __FITPOT__
    close(81)
    allocate( aml(3,nsf,nhl1,natm+nb),bml(3,nsf,nhl1,natm+nb) )
    call copy_dba_fwd(tcom,namax,natm,nbmax,nb,lsb,lsrc,myparity &
         ,nn,mpi_world,gsf,nsf)
!.....make lumps of terms
    aml(1:3,1:nsf,1:nhl1,1:natm+nb)= 0d0
    bml(1:3,1:nsf,1:nhl1,1:natm+nb)= 0d0
    do ia=1,natm
      do jj=1,lspr(0,ia)
        ja= lspr(jj,ia)
        do ihl1=1,nhl1 ! no need of a bias node, 0
          hl1i= hl1(ihl1,ia)
          tmp= hl1i*(1d0-hl1i)
          tmp2= hl1i*(1d0-hl1i)*(1d0-2d0*hl1i)
          do isf=1,nsf ! no need of a bias node, 0
            aml(1:3,isf,ihl1,ja)= aml(1:3,isf,ihl1,ja) &
                 +tmp*dgsf(1:3,isf,jj,ia)
            tmp3(1:3)= 0d0
            do jsf=1,nsf
              tmp3(1:3)= tmp3(1:3) +dgsf(1:3,jsf,jj,ia)*wgt1(ihl1,jsf)
            enddo
            bml(1:3,isf,ihl1,ja)= bml(1:3,isf,ihl1,ja) &
                 +tmp2*gsf(isf,ja)*tmp3(1:3)
            if( abs(bml(1,isf,ihl1,ja)).gt.1d+10 ) then
              write(6,'(a,5i6)') ' ia,jj,ja,ihl1,isf=',ia,jj,ja,ihl1,isf
              write(6,'(a,5es12.4)') ' tmp2,gsf,tmp3=',tmp2,gsf(isf,ja),tmp3(1:3)
            endif
          enddo
        enddo
      enddo
    enddo
!.....copy back Ms of buffer atoms
    call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,lsrc,myparity &
         ,nn,mpi_world,aml,3*nsf*nhl1)
    call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,lsrc,myparity &
         ,nn,mpi_world,bml,3*nsf*nhl1)
!.....write aml
    open(82,file='out.NN1.aml')
    write(82,'(3i10)') natm, nsf, nhl1
    do ia=1,natm
      do ihl1=1,nhl1
        do isf=1,nsf
          write(82,'(i6,2i4,3es23.14e3)') ia,ihl1,isf &
               ,aml(1:3,isf,ihl1,ia)
        enddo
      enddo
    enddo
    close(82)
!.....write bml
    open(83,file='out.NN1.bml')
    write(83,'(3i10)') natm, nsf, nhl1
    do ia=1,natm
      do ihl1=1,nhl1
        do isf=1,nsf
          write(83,'(i6,2i4,3es23.14e3)') ia,ihl1,isf &
               ,bml(1:3,isf,ihl1,ia)
        enddo
      enddo
    enddo
    close(83)
    deallocate(aml,bml)
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
    real(8),intent(out):: gsf(0:nsf,natm),dgsf(3,nsf,0:nnmax,namax)

    integer:: isf,ia,jj,ja,kk,ka
    real(8):: xi(3),xj(3),xij(3),rij(3),dij,fcij,eta,rs,texp,driji(3), &
         dfcij,drijj(3),dgdr,xk(3),xik(3),rik(3),dik,fcik,dfcik, &
         driki(3),drikk(3),almbd,spijk,cs,t1,t2,dgdij,dgdik,dgcs, &
         dcsdj(3),dcsdk(3),dcsdi(3)

    real(8),external:: sprod

    gsf(0:nsf,1:natm)= 0d0
    gsf(0,1:natm)= 1d0
    dgsf(1:3,1:nsf,0:nnmax,1:natm)= 0d0

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
            texp= exp(-eta*(dij-rs)**2)
            gsf(isf,ia)= gsf(isf,ia) +texp*fcij
            !.....derivative
            driji(1:3)= -rij(1:3)/dij
            drijj(1:3)= -driji(1:3)
            dgdr= -2d0*eta*(dij-rs)*texp*fcij +texp*dfc(dij,rc)
            dgsf(1:3,isf,0,ia)= dgsf(1:3,isf,0,ia) +driji(1:3)*dgdr
            dgsf(1:3,isf,jj,ia)= dgsf(1:3,isf,jj,ia) +drijj(1:3)*dgdr
          enddo
        enddo

      else if( itype(isf).eq.2 ) then ! angular (3-body)
        almbd= cnst(1,isf)
        t2= (abs(almbd)+1d0)**2
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
              !.....function value
              spijk= rij(1)*rik(1) +rij(2)*rik(2) +rij(3)*rik(3)
              cs= spijk/dij/dik
              t1= (almbd +cs)**2
              gsf(isf,ia)= gsf(isf,ia) +t1/t2 *fcij*fcik 
              !.....derivative
              dgdij= dfcij *fcik *t1/t2
              dgdik= fcij *dfcik *t1/t2
              dgsf(1:3,isf,0,ia)= dgsf(1:3,isf,0,ia) &
                   +dgdij*driji(1:3) +dgdik*driki(1:3)
              dgsf(1:3,isf,jj,ia)= dgsf(1:3,isf,jj,ia) +dgdij*drijj(1:3)
              dgsf(1:3,isf,kk,ia)= dgsf(1:3,isf,kk,ia) +dgdik*drikk(1:3)
              dgcs= 2d0*(almbd+cs)/t2 *fcij*fcik
              dcsdj(1:3)= rik(1:3)/dij/dik -rij(1:3)*spijk/dij**3/dik
              dcsdk(1:3)= rij(1:3)/dij/dik -rik(1:3)*spijk/dik**3/dij
              dcsdi(1:3)= -dcsdj(1:3) -dcsdk(1:3)
              dgsf(1:3,isf,0,ia)= dgsf(1:3,isf,0,ia) +dgcs*dcsdi(1:3)
              dgsf(1:3,isf,jj,ia)= dgsf(1:3,isf,jj,ia) +dgcs*dcsdj(1:3)
              dgsf(1:3,isf,kk,ia)= dgsf(1:3,isf,kk,ia) +dgcs*dcsdk(1:3)
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

!.....calc number of weights taking into account the bias nodes
    nwgt1= (nsf+1)*nhl1
    nwgt2= (nhl1+1)
    if( myid.eq.0 ) then
      write(6,'(a,4i6)') ' nfs, nhl1, nwgt1, nwgt2 =', &
           nsf,nhl1,nwgt1,nwgt2
    endif

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
!=======================================================================
  subroutine copy_dba_fwd(tcom,namax,natm,nbmax,nb,lsb,lsrc,myparity &
       ,nn,mpi_world,x,ndim)
!   Send forward data of resident atomd to buffer atoms
    implicit none 
    include 'mpif.h'
    integer,intent(in):: namax,natm,nbmax,nb,mpi_world,ndim
    integer,intent(in):: lsb(0:nbmax,6),lsrc(6),myparity(3),nn(6)
    real(8),intent(inout):: x(ndim,namax),tcom

    integer:: status(mpi_status_size)
    integer:: i,j,k,nbnew,kd,kdd,ku,inode,nsd,nrc,ierr
    real(8):: tcom1,tcom2
    real(8),allocatable,dimension(:,:):: dbuf,dbufr
    
    allocate(dbuf(ndim,nbmax),dbufr(ndim,nbmax))

    nbnew= 0
    do kd=1,3
      tcom1= mpi_wtime()

      do kdd=-1,0
        ku= 1*kd+kdd
        inode= nn(ku)
        nsd= lsb(0,ku)

        call mespasi(inode,myparity(kd),nsd,nrc,1,1,10 &
             ,mpi_world)
        
        do i=1,nsd
          j= lsb(i,ku)
          dbuf(1:ndim,i)= x(1:ndim,j)
        enddo

        call mespasd(inode,myparity(kd),dbuf,dbufr,nsd*ndim,nrc*ndim &
             ,21,mpi_world)

        do i=1,nrc
          x(1:ndim,natm+nbnew+i)= dbufr(1:ndim,i)
        enddo
        call mpi_barrier(mpi_world,ierr)
        nbnew= nbnew +nrc
      enddo

      tcom2= mpi_wtime()
      tcom= tcom +tcom2-tcom1
    enddo
    deallocate(dbuf,dbufr)
  end subroutine copy_dba_fwd
end module NN1
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
