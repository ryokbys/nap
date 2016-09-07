module NN
!-----------------------------------------------------------------------
!                        Time-stamp: <2016-09-06 16:34:13 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
!  Parallel implementation of neural-network potential with 1 hidden
!  layer. It is available for plural number of species.
!-----------------------------------------------------------------------
!.....parameter file name
  character(128),parameter:: cpfname= 'in.params.NN'
  character(128),parameter:: ccfname='in.const.NN'
!.....parameters
  integer:: nwgt1,nwgt2
  real(8),allocatable:: wgt11(:,:),wgt12(:)
  real(8),allocatable:: wgt21(:,:),wgt22(:,:),wgt23(:)
!.....constants
  integer,parameter:: nlmax= 2
  integer:: nsfc,nsfc1,nsfc2,nc1,nc2,nsp,nl,nhl(0:nlmax+1)
  integer,allocatable:: itype(:)
  real(8),allocatable:: cnst(:,:)
  real(8),allocatable:: hl1(:,:),hl2(:,:)
  real(8),allocatable:: gsf(:,:),dgsf(:,:,:,:)
  integer,allocatable:: icmb2(:,:),icmb3(:,:,:)
  integer,allocatable:: iaddr2(:,:,:),iaddr3(:,:,:,:)
!.....function types and num of constatns for types
  integer,parameter:: max_ncnst= 2
  integer:: ncnst_type(200)
  integer:: ncomb_type(200)

!.....max exponent of the basis function
  integer:: max_nexp

!.....cutoff region width ratio to rc
  real(8):: rcw = 0.9d0

!.....num of atoms and neighbors for dgsf array
  integer,save:: nal, nnl, nalmax,nnlmax,nnltmp
  logical:: lrealloc = .false.
  
contains
  subroutine force_NN(namax,natm,tag,ra,nnmax,aa,strs,h,hi,tcom &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_world,myid,epi,epot,nismax,acon,lstrs,iprint)
    implicit none
    include "mpif.h"
    include "./params_unit.h"
    integer,intent(in):: namax,natm,nnmax,nismax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_world,myid,lspr(0:nnmax,namax),nex(3)
    real(8),intent(in):: ra(3,namax),tag(namax),acon(nismax) &
         ,h(3,3),hi(3,3),sv(3,6)
    real(8),intent(inout):: tcom,rc
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical:: lstrs

!.....local
    integer:: i,j,k,l,m,n,ixyz,jxyz,is,js,ks,ierr,nbl,ia,ja,nexp,isf &
         ,icoeff,ihl0,ihl1,ihl2,jj,jsf
    real(8):: rcin,b_na,at(3),epotl,wgt,hl1i,hl2i,tmp2,tmp1,tmp,tmp3(3)
    real(8),save:: rc3
!    real(8),allocatable:: aml(:,:,:,:),bml(:,:,:,:)
!.....1st call
    logical,save:: l1st=.true.

    if( l1st ) then
!.....read in.params.NN
      call read_params(myid,mpi_world,rcin,rc3)
!.....reset rc
      if( myid.le.0 .and. rc .lt. rcin- 1d-8) then
        write(6,'(a,f10.5,a,f10.5)') &
             ' Error: Cutoff radius rc should be corrected from '&
             ,rc,' to ',rcin
        if( myid.ge.0 ) then
          call mpi_finalize(ierr)
          stop
        else
          stop
        endif
      endif
      rc= rcin

!  To reduce the memory usage, compute num of atoms and num of neighbors,
!  and add some margin for those numbers because they can change during
!  the simulation.
      nal = int(natm*1.1)
      nnltmp = 0
      do i=1,natm
        nnltmp = max(nnltmp,lspr(0,i))
      enddo
      nnl = int(nnltmp*1.1)
      if( nal .gt. namax ) then
        write(6,'(a)') ' [Error] nal .gt.namax'
        write(6,'(a,3i10)') '   myid,nal,namax = ',myid,nal,namax
        stop
      endif
      if( nnl.gt.nnmax ) then
        write(6,'(a)') ' [Error] nnl.gt.nnmax'
        write(6,'(a,3i10)') '   myid,nnl,nnmax = ',myid,nnl,nnmax
        stop
      endif
      if( myid.ge.0 ) then
        call mpi_reduce(nal,nalmax,1,mpi_integer,mpi_max,0,mpi_world,ierr)
        call mpi_reduce(nnl,nnlmax,1,mpi_integer,mpi_max,0,mpi_world,ierr)
      else
        nalmax = nal
        nnlmax = nnl
      endif
      if( myid.le.0 ) then
        write(6,'(a,2i10)') ' max num of (local atoms *1.1) = ',nalmax
        write(6,'(a,2i10)') ' max num of (neighbors *1.1)   = ',nnlmax
        write(6,'(a,i10,a)') ' gsf size  = ', &
             nhl(0)*nalmax*8/1000/1000,' MB'
        write(6,'(a,i10,a)') ' dgsf size = ', &
             int(3*nhl(0),8)*(nnlmax+1)*nalmax*8/1000/1000,' MB'
      endif
      allocate( gsf(nhl(0),nal),dgsf(3,nhl(0),0:nnl,nal) )

      lrealloc = .false.
      if( nl.eq.1 ) then
        allocate( hl1(nhl(1),nal) )
      else if( nl.eq.2 ) then
        allocate( hl1(nhl(1),nal), hl2(nhl(2),nal) )
      endif
      l1st= .false.
    endif

!  Since natm and nn can change every step of MD,
!  if natm/nnltmp becomes nal/nnl, they should be updated and
!  gsf/dgsf as well.
    if( natm.gt.nal ) then
      nal = int(natm*1.1)
      if( nal .gt. namax ) then
        write(6,'(a)') ' [Error] nal .gt.namax'
        write(6,'(a,3i10)') '   myid,nal,namax = ',myid,nal,namax
        stop
      endif
      lrealloc=.true.
    endif
    nnltmp = 0
    do i=1,natm
      nnltmp = max(nnltmp,lspr(0,i))
    enddo
    if( nnltmp.gt.nnl ) then
      nnl = int(nnltmp*1.1)
      if( nnlmax.gt.nnmax ) then
        write(6,'(a)') ' [Error] nnl.gt.nnmax'
        write(6,'(a,3i10)') '   myid,nnl,nnmax = ',myid,nnl,nnmax
        stop
      endif
      lrealloc=.true.
    endif
    if( allocated(dgsf).and.lrealloc ) then
      deallocate( gsf,dgsf )
      allocate( gsf(nhl(0),nal),dgsf(3,nhl(0),0:nnl,nal) )
      if( nl.eq.1 ) then
        deallocate( hl1 )
        allocate( hl1(nhl(1),nal) )
      else if( nl.eq.2 ) then
        deallocate( hl1,hl2 )
        allocate( hl1(nhl(1),nal), hl2(nhl(2),nal) )
      endif
      lrealloc=.false.
    endif

!.....first, calculate all the symmetry functions
    call eval_sf(nhl(0),namax,natm,nb,nnmax,h,tag,ra &
         ,lspr,rc,rc3)

    if( mod(iprint,100)/10.eq.1 .and. myid.le.0 ) then
      open(80,file='out.NN.gsf',status='replace',form='unformatted')
      write(80) nhl(0)
      do ia=1,natm
        write(80) (gsf(ihl0,ia),ihl0=1,nhl(0))
      enddo
      close(80)
      call write_dgsf(84,natm,namax,nnmax,lspr,tag,nhl(0))
    endif

!.....2nd, calculate the node values by summing contributions from
!.....  symmetry functions
    if( nl.eq.1 ) then
      hl1(1:nhl(1),1:natm)= 0d0
      do ia=1,natm
!.....debug
        if( iprint.ge.100 ) then
          print *,"ia=",ia
          do ihl1=1,nhl(1)
            tmp= 0d0
            do ihl0=1,nhl(0)
              tmp= tmp +wgt11(ihl0,ihl1) *gsf(ihl0,ia)
              write(6,'(a,i4,3es15.7)') "ihl0,wgt1*gsf=",ihl0 &
                   ,wgt11(ihl0,ihl1),gsf(ihl0,ia) &
                   ,wgt11(ihl0,ihl1) *gsf(ihl0,ia)
            enddo
            hl1(ihl1,ia)= sigmoid(tmp)
            print *,"ihl1,hl1=",ihl1,hl1(ihl1,ia)
          enddo
        endif !end of debug
        do ihl1=1,nhl(1)
          tmp= 0d0
          do ihl0=1,nhl(0)
            tmp= tmp +wgt11(ihl0,ihl1) *gsf(ihl0,ia)
          enddo
          hl1(ihl1,ia)= sigmoid(tmp)
        enddo
      enddo
    else if( nl.eq.2 ) then
      hl1(1:nhl(1),1:natm)= 0d0
      hl2(1:nhl(2),1:natm)= 0d0
      do ia=1,natm
        do ihl1=1,nhl(1)
          tmp= 0d0
          do ihl0=1,nhl(0)
            tmp= tmp +wgt21(ihl0,ihl1) *gsf(ihl0,ia)
          enddo
          hl1(ihl1,ia)= sigmoid(tmp)
        enddo
        do ihl2=1,nhl(2)
          tmp= 0d0
          do ihl1=1,nhl(1)
            tmp= tmp +wgt22(ihl1,ihl2) *(hl1(ihl1,ia)-0.5d0)
          enddo
          hl2(ihl2,ia)= sigmoid(tmp)
        enddo
      enddo

    endif

!.....then calculate the energy of atom by summing up the node values
    epotl= 0d0
    if( nl.eq.1 ) then
      do ia=1,natm
        epi(ia)= 0d0
        do ihl1=1,nhl(1)
          epi(ia)= epi(ia) +wgt12(ihl1) *(hl1(ihl1,ia)-0.5d0)
        enddo
        epotl=epotl +epi(ia)
#ifdef __3BODY__
        write(6,'(a,i8,es22.14)') ' 3-body term:',ia,epi(ia)
#endif
      enddo
    else if( nl.eq.2 ) then
      do ia=1,natm
        epi(ia)= 0d0
        do ihl2=1,nhl(2)
          epi(ia)= epi(ia) +wgt23(ihl2) *(hl2(ihl2,ia)-0.5d0)
        enddo
        epotl=epotl +epi(ia)
#ifdef __3BODY__
        write(6,'(a,i8,es22.14)') ' 3-body term:',ia,epi(ia)
#endif
      enddo
    endif

!.....sum up for forces
    aa(1:3,1:natm+nb)= 0d0
    if( nl.eq.1 ) then
      do ia=1,natm
        do ihl1=1,nhl(1)
          hl1i= hl1(ihl1,ia)
          tmp= wgt12(ihl1)*hl1i*(1d0-hl1i)
          do jj=1,lspr(0,ia)
            ja= lspr(jj,ia)
            do ihl0=1,nhl(0)
              aa(1:3,ja)=aa(1:3,ja) &
                   -tmp*wgt11(ihl0,ihl1)*dgsf(1:3,ihl0,jj,ia)
            enddo
          enddo
          !.....atom ia
          do ihl0= 1,nhl(0)
            aa(1:3,ia)=aa(1:3,ia) &
                 -tmp*wgt11(ihl0,ihl1)*dgsf(1:3,ihl0,0,ia)
          enddo
        enddo
      enddo
    else if( nl.eq.2 ) then
      do ia=1,natm
        do ihl2=1,nhl(2)
          hl2i= hl2(ihl2,ia)
          tmp2= wgt23(ihl2) *hl2i*(1d0-hl2i)
          do ihl1=1,nhl(1)
            hl1i= hl1(ihl1,ia)
            tmp1= wgt22(ihl1,ihl2) *hl1i*(1d0-hl1i)
            do jj=1,lspr(0,ia)
              ja= lspr(jj,ia)
              do ihl0=1,nhl(0)
                aa(1:3,ja)=aa(1:3,ja) &
                     -tmp2 *tmp1 &
                     *wgt21(ihl0,ihl1)*dgsf(1:3,ihl0,jj,ia)
              enddo
            enddo
!.....atom ia
            do ihl0= 1,nhl(0)
              aa(1:3,ia)=aa(1:3,ia) &
                   -tmp2*tmp1*wgt21(ihl0,ihl1)*dgsf(1:3,ihl0,0,ia)
            enddo
          enddo
        enddo
      enddo
    endif

    call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,aa,3)
!!$    if( myid.ge.0 ) then
!!$      call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,lsrc,myparity &
!!$           ,nn,mpi_world,aa,3)
!!$    else
!!$      call reduce_dba_bk(natm,namax,tag,aa,3)
!!$    endif

!-----reduced force
    do i=1,natm
      at(1:3)= aa(1:3,i)
      aa(1:3,i)= hi(1:3,1)*at(1) +hi(1:3,2)*at(2) +hi(1:3,3)*at(3)
!      aa(1:3,i)= hi(1,1:3)*at(1) +hi(2,1:3)*at(2) +hi(3,1:3)*at(3)
    enddo
!-----multiply 0.5d0*dt**2/am(i)
    do i=1,natm
      is= int(tag(i))
      aa(1:3,i)= acon(is)*aa(1:3,i)
    enddo


    if( lstrs ) then
      call compute_stress(namax,natm,tag,ra,nnmax,strs,h &
           ,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nn,rc,lspr &
           ,mpi_world,myid)
    endif

!-----gather epot
    epot= 0d0
    if( myid.ge.0 ) then
      call mpi_allreduce(epotl,epot,1,mpi_double_precision &
           ,mpi_sum,mpi_world,ierr)
    else
      epot= epotl
    endif
    return
  end subroutine force_NN
!=======================================================================
  subroutine eval_sf(nsf,namax,natm,nb,nnmax,h,tag,ra,lspr,rc,rc3)
!
!  Evaluate symmetry functions and derivatives for multi-species system.
!
    implicit none
    integer,intent(in):: nsf,namax,natm,nb,nnmax,lspr(0:nnmax,namax)
    real(8),intent(in):: h(3,3),tag(namax),ra(3,namax),rc,rc3
!    real(8),intent(out):: gsf(nsf,nal),dgsf(3,nsf,0:nnl,nal)

    integer:: isf,isfc,ia,jj,ja,kk,ka,is,js,ks,isfc1,isfc2
    real(8):: xi(3),xj(3),xij(3),rij(3),dij,fcij,eta,rs,texp,driji(3), &
         dfcij,drijj(3),dgdr,xk(3),xik(3),rik(3),dik,fcik,dfcik, &
         driki(3),drikk(3),almbd,spijk,cs,t1,t2,dgdij,dgdik,dgcs, &
         dcsdj(3),dcsdk(3),dcsdi(3),tcos,tpoly,a1,a2,tmorse

    real(8),external:: sprod

    gsf(1:nsf,1:nal)= 0d0
    dgsf(1:3,1:nsf,0:nnl,1:nal)= 0d0
    do ia=1,natm
      xi(1:3)= ra(1:3,ia)
      is= int(tag(ia))
      do jj=1,lspr(0,ia)
        ja= lspr(jj,ia)
        if( ja.eq.ia ) cycle
        xj(1:3)= ra(1:3,ja)
        xij(1:3)= xj(1:3)-xi(1:3)
        rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij= sqrt(rij(1)**2 +rij(2)**2 +rij(3)**2)
        if( dij.ge.rc ) cycle
        js= int(tag(ja))
        isfc=0
        driji(1:3)= -rij(1:3)/dij
        drijj(1:3)= -driji(1:3)
        fcij= fc(dij,rc)
        dfcij= dfc(dij,rc)
        do isf=iaddr2(1,is,js),iaddr2(2,is,js)
!!$          isfc= isfc+1
!!$          isf= (icmb2(is,js)-1)*nsfc1 +isfc1
          if( itype(isf).eq.1 ) then ! Gaussian
            eta= cnst(1,isf)
            rs=  cnst(2,isf)
            !.....function value
            texp= exp(-eta*(dij-rs)**2)
            gsf(isf,ia)= gsf(isf,ia) +texp*fcij
            !.....derivative
            ! dgsf(ixyz,isf,jj,ia): derivative of isf-th basis of atom-ia
            ! by ixyz coordinate of atom-jj.
            ! jj=0 means derivative by atom-ia.
            dgdr= -2d0*eta*(dij-rs)*texp*fcij +texp*dfcij
            dgsf(1:3,isf,0,ia)= dgsf(1:3,isf,0,ia) +driji(1:3)*dgdr
            dgsf(1:3,isf,jj,ia)= dgsf(1:3,isf,jj,ia) +drijj(1:3)*dgdr
          else if( itype(isf).eq.2 ) then ! cosine
            a1= cnst(1,isf)
            !.....func value
            tcos= (1d0+cos(dij*a1))
            gsf(isf,ia)= gsf(isf,ia) +tcos*fcij
            !.....derivative
            dgdr= -a1*sin(dij*a1)*fcij +tcos*dfc(dij,rc)
            dgsf(1:3,isf,0,ia)= dgsf(1:3,isf,0,ia) +driji(1:3)*dgdr
            dgsf(1:3,isf,jj,ia)= dgsf(1:3,isf,jj,ia) +drijj(1:3)*dgdr
          else if( itype(isf).eq.3 ) then ! polynomial
            a1= cnst(1,isf)
            !.....func value
            tpoly= 1d0*dij**(-a1)
            gsf(isf,ia)= gsf(isf,ia) +tpoly*fcij
            !.....derivative
            dgdr= -a1*dij**(-a1-1d0)*fcij +tpoly*dfc(dij,rc)
            dgsf(1:3,isf,0,ia)= dgsf(1:3,isf,0,ia) +driji(1:3)*dgdr
            dgsf(1:3,isf,jj,ia)= dgsf(1:3,isf,jj,ia) +drijj(1:3)*dgdr
          else if( itype(isf).eq.4 ) then ! Morse-type
            a1= cnst(1,isf)
            a2= cnst(2,isf)
            !.....func value
            texp= exp(-a1*(dij-a2))
            tmorse= ((1d0-texp)**2 -1d0)
            gsf(isf,ia)= gsf(isf,ia) +tmorse*fcij
            !.....derivative
            dgdr= 2d0*a1*(1d0-texp)*texp*fcij +tmorse*dfcij
            dgsf(1:3,isf,0,ia)= dgsf(1:3,isf,0,ia) +driji(1:3)*dgdr
            dgsf(1:3,isf,jj,ia)= dgsf(1:3,isf,jj,ia) +drijj(1:3)*dgdr
          endif
        enddo

!!$        fcij= fc(dij,rc)
!!$        dfcij= dfc(dij,rc)
!!$        driji(1:3)= -rij(1:3)/dij
!!$        drijj(1:3)= -driji(1:3)
        if( dij.gt.rc3 ) cycle
        do kk=1,lspr(0,ia)
          ka= lspr(kk,ia)
          ks= int(tag(ka))
          if( iaddr3(1,is,js,ks).lt.0 ) cycle
          if( ka.eq.ia .or. ka.le.ja ) cycle
          xk(1:3)= ra(1:3,ka)
          xik(1:3)= xk(1:3)-xi(1:3)
          rik(1:3)= h(1:3,1)*xik(1) +h(1:3,2)*xik(2) +h(1:3,3)*xik(3)
          dik= sqrt(rik(1)**2 +rik(2)**2 +rik(3)**2)
          if( dik.ge.rc3 ) cycle
          do isf=iaddr3(1,is,js,ks),iaddr3(2,is,js,ks)
!!$            isf= nsfc1*nc1 +(icmb3(is,js,ks)-1)*nsfc2 +isfc2
            almbd= cnst(1,isf)
            t2= (abs(almbd)+1d0)**2
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
    enddo

  end subroutine eval_sf
!=======================================================================
  function fc(r,rc)
    implicit none
    real(8),intent(in):: r,rc
    real(8):: fc,rs
    real(8),parameter:: pi= 3.14159265358979d0
    rs= rc*rcw
    if( r.le.rs ) then
      fc= 1d0
    else if( r.gt.rs .and. r.le.rc ) then
      fc= 0.5d0 *(cos((r-rs)/(rc-rs)*pi)+1d0)
    else
      fc= 0d0
    endif
    return
  end function fc
!=======================================================================
  function dfc(r,rc)
    implicit none
    real(8),intent(in):: r,rc
    real(8):: dfc,rs
    real(8),parameter:: pi= 3.14159265358979d0
    rs= rc*rcw
    if( r.le.rs ) then
      dfc= 0d0
    else if( r.gt.rs .and. r.le.rc ) then
      dfc= -pi/2/(rc-rs) *sin((r-rs)/(rc-rs)*pi)
    else
      dfc= 0d0
    endif
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
  subroutine read_params(myid,mpi_world,rcin,rc3)
    implicit none
    include 'mpif.h'

    integer,intent(in):: myid,mpi_world
    real(8),intent(out):: rcin,rc3
    integer:: itmp,ierr,i,j,k,nc,ncoeff,is,js,ks &
         ,n,ihl0,ihl1,ihl2,icmb(3),nsf,nsf1,nsf2,iap,jap,kap
    integer,allocatable:: nwgt(:)
    logical:: lexist

!.....initialize some
    ncnst_type(1)= 2   ! Gaussian
    ncnst_type(2)= 1   ! cosine
    ncnst_type(3)= 1   ! polynomial
    ncnst_type(4)= 2   ! Morse
    ncnst_type(101)= 1 ! angular

    ncomb_type(1:100)= 2    ! pair
    ncomb_type(101:200)= 3  ! triplet

!.....read constants at the 1st call
    inquire(file=trim(ccfname),exist=lexist)
    if( .not. lexist ) then
      if( myid.ge.0 ) then
        if( myid.eq.0 ) then
          write(6,'(a)') ' [Error] '//ccfname//' does not exist !!!.'
          write(6,'(a)') '   The NN potential needs '//ccfname//'.'
        endif
        call mpi_finalize(ierr)
        stop
      else
        write(6,'(a)') ' [Error] '//ccfname//' does not exist !!!.'
        write(6,'(a)') '   The NN potential needs '//ccfname//'.'
        stop
      endif
    endif
    open(51,file=trim(ccfname),status='old')
!.....num of symmetry functions, num of node in 1st hidden layer
    read(51,*) nl,nsp,(nhl(i),i=0,nl)
!    print *,' nl,nsp,(nhl(i),i=0,nl)=',nl,nsp,(nhl(i),i=0,nl)
    if( nl.gt.nlmax ) then
      if( myid.ge.0 ) then
        if( myid.eq.0 ) then
          print *, '[Error] nl.gt.nlmax '
          print *, '  nl,nlmax=',nl,nlmax
        endif
        call mpi_finalize(ierr)
        stop
      else
        print *, '[Error] nl.gt.nlmax '
        print *, '  nl,nlmax=',nl,nlmax
        stop
      endif
    endif
    nsf= nhl(0)
    nhl(nl+1)= 1
    allocate(itype(nsf),cnst(max_ncnst,nsf))
    allocate(iaddr2(2,nsp,nsp),iaddr3(2,nsp,nsp,nsp))
    iaddr2(1:2,1:nsp,1:nsp)= -1
    iaddr3(1:2,1:nsp,1:nsp,1:nsp)= -1
    nsf1= 0
    nsf2= 0
    iap= 0
    jap= 0
    kap= 0
    do i=1,nsf
      read(51,*) itype(i),(icmb(k),k=1,ncomb_type(itype(i))) &
           ,(cnst(j,i),j=1,ncnst_type(itype(i)))
      if( itype(i).le.100 ) then
        if( icmb(1).ne.iap .or. icmb(2).ne.jap ) then
          iaddr2(1,icmb(1),icmb(2))= i
          iaddr2(1,icmb(2),icmb(1))= i
        endif
        iaddr2(2,icmb(1),icmb(2))= i
        iaddr2(2,icmb(2),icmb(1))= i
        nsf1= nsf1 +1
        iap= icmb(1)
        jap= icmb(2)
      else if( itype(i).le.200 ) then
        if( icmb(1).ne.iap .or. icmb(2).ne.jap .or. &
             icmb(3).ne.kap ) then
          iaddr3(1,icmb(1),icmb(2),icmb(3))= i
          iaddr3(1,icmb(1),icmb(3),icmb(2))= i
        endif
        iaddr3(2,icmb(1),icmb(2),icmb(3))= i
        iaddr3(2,icmb(1),icmb(3),icmb(2))= i
        nsf2= nsf2 +1
        iap= icmb(1)
        jap= icmb(2)
        kap= icmb(3)
      endif
    enddo
    if( nsf.ne.nsf1+nsf2 ) then
      if(myid.eq.0) then
        print *,'[Error] nsf.ne.nsf1+nsf2 !!!'
      endif
      call mpi_finalize(ierr)
      stop
    endif
    close(51)

!.....calc number of weights
!!$    nc1= nsp +factorial(nsp,2)/2
!!$    nc2= nsp*nc1
!!$    nhl(0)= nsfc1*nc1 +nsfc2*nc2
!!$    nwgt1= nsf*nhl1
!!$    nwgt2= nhl1
    allocate(nwgt(nl+1))
    do i=1,nl+1
      nwgt(i)= nhl(i-1)*nhl(i)
!      print *,' i,nhl(i-1),nhl(i),nwgt(i)=',i,nhl(i-1),nhl(i),nwgt(i)
    enddo
    if( myid.le.0 ) then
      print *, 'num of basis funcs =',nhl(0)
      do i=1,nl
        print *, 'ihl, nhl(ihl)  =',i,nhl(i)
      enddo
      do i=1,nl+1
        print *, 'ihl, nwgt(ihl)  =',i,nwgt(i)
      enddo
    endif

!.....read parameters at the 1st call
    inquire(file=trim(cpfname),exist=lexist)
    if( .not. lexist ) then
      if( myid.ge.0 ) then
        if( myid.eq.0 ) then
          write(6,'(a)') ' [Error] '//cpfname//' does not exist !!!.'
          write(6,'(a)') '   The NN potential needs '//cpfname//'.'
        endif
        call mpi_finalize(ierr)
        stop
      else
        write(6,'(a)') ' [Error] '//cpfname//' does not exist !!!.'
        write(6,'(a)') '   The NN potential needs '//cpfname//'.'
        stop
      endif
    endif
    open(50,file=trim(cpfname),status='old')
    read(50,*) ncoeff,rcin,rc3
!.....check whether the num of parameters is correct
    if( rc3.gt.rcin ) then
      rc3= rcin
      if( myid.le.0 ) then
        write(6,*) ' rc3 was corrected to rcin = ',rcin
        write(6,*) ' because input rc3 > rc, which should not happen.'
      endif
    endif
    
    nc= 0
    do i=1,nl+1
      nc= nc +nwgt(i)
    enddo
    if( ncoeff .ne. nc ) then
      write(6,'(a)') ' [Error] num of parameters is not correct !!!'
      write(6,'(a,i10)')  '   ncoeff=',ncoeff
      write(6,'(a,i10)')  '   ncoeff should be ',nc
      stop
    endif
!.....different number of weights for different number of layers
    if( nl.eq.1 ) then
      allocate(wgt11(nhl(0),nhl(1)),wgt12(nhl(1)))
      wgt11(1:nhl(0),1:nhl(1)) = 0d0
      wgt12(1:nhl(1)) = 0d0
    else if( nl.eq.2 ) then
      allocate(wgt21(nhl(0),nhl(1)),wgt22(nhl(1),nhl(2)),wgt23(nhl(2)))
      wgt21(1:nhl(0),1:nhl(1)) = 0d0
      wgt22(1:nhl(1),1:nhl(2)) = 0d0
      wgt23(1:nhl(2)) = 0d0
    endif
    if( nl.eq.1 ) then
      do ihl0=1,nhl(0)
        do ihl1=1,nhl(1)
          read(50,*) wgt11(ihl0,ihl1)
        enddo
      enddo
      do ihl1=1,nhl(1)
        read(50,*) wgt12(ihl1)
      enddo
    else if( nl.eq.2 ) then
      do ihl0=1,nhl(0)
        do ihl1=1,nhl(1)
          read(50,*) wgt21(ihl0,ihl1)
        enddo
      enddo
      do ihl1=1,nhl(1)
        do ihl2=1,nhl(2)
          read(50,*) wgt22(ihl1,ihl2)
        enddo
      enddo
      do ihl2=1,nhl(2)
        read(50,*) wgt23(ihl2)
      enddo
    endif
    close(50)

#ifdef __DEBUG__
    if(myid.le.0) then
      write(6,'(a)') ' DEBUG: ihl0,ihl1,wgt11(ihl0,ihl1)'
      do ihl0=1,nhl(0)
        do ihl1=1,nhl(1)
          write(6,'(2i5,es15.7)') ,ihl0,ihl1,wgt11(ihl0,ihl1)
        enddo
      enddo
    endif
#endif

!!$!.....read in.comb.NN
!!$    allocate(icmb2(nsp,nsp),icmb3(nsp,nsp,nsp))
!!$    inquire(file=trim(cmbfname),exist=lexist)
!!$    if( .not.lexist ) then
!!$      if( myid.ge.0 ) then
!!$        if( myid.eq.0 ) then
!!$          write(6,'(a)') ' [Error] '//cmbfname//' does not exist !!!.'
!!$          write(6,'(a)') '   The NN potential needs '//cmbfname//'.'
!!$        endif
!!$        call mpi_finalize(ierr)
!!$        stop
!!$      else
!!$        write(6,'(a)') ' [Error] '//cmbfname//' does not exist !!!.'
!!$        write(6,'(a)') '   The NN potential needs '//cmbfname//'.'
!!$        stop
!!$      endif
!!$    else
!!$      open(52,file=trim(cmbfname),status='old')
!!$!.....read pairs
!!$      do n=1,nc1
!!$        read(52,*) i,j,icmb2(i,j)
!!$        icmb2(j,i)= icmb2(i,j)
!!$      enddo
!!$!.....read triplets
!!$      do n=1,nc2
!!$        read(52,*) i,j,k,icmb3(i,j,k)
!!$        icmb3(i,k,j)= icmb3(i,j,k)
!!$      enddo
!!$      close(52)
!!$    endif

    deallocate(nwgt)
    return
  end subroutine read_params
!=======================================================================
  function factorial(n,m)
!  compute factorial of n, m-times.
    implicit none
    integer,intent(in):: n,m
    real(8):: factorial

    integer:: i

    factorial= 1
    do i=0,m-1
      factorial= factorial*(n-i)
    enddo
    return
  end function factorial
!=======================================================================
  subroutine write_dgsf(ionum,natm,namax,nnmax,lspr,tag,nsf)
!   Write out dgsf data.
!   Buffer atom indices are replaced to resident atom ones.
    implicit none
    integer,intent(in):: ionum
    integer,intent(in):: natm,namax,nnmax,nsf,lspr(0:nnmax,namax)
!    real(8),intent(in):: dgsf(3,nsf,0:nnl,nal),tag(namax)
    real(8),intent(in):: tag(namax)
    integer:: ia,jj,ja,jra,isf
    real(8),allocatable:: dgsfo(:,:,:,:)
    integer,external:: itotOf

    allocate(dgsfo(3,natm,nsf,natm))
!.....reduce dgsf data of buffer atoms to those of resident atoms
    dgsfo(1:3,1:natm,1:nsf,1:natm)= 0d0
    do ia=1,natm
      do jj=0,lspr(0,ia)
        if( jj.eq.0 ) then
          ja= ia
        else
          ja= lspr(jj,ia)
        endif
        jra= itotOf(tag(ja))
        do isf=1,nsf
          dgsfo(1:3,jra,isf,ia)= dgsfo(1:3,jra,isf,ia) &
               +dgsf(1:3,isf,jj,ia)
        enddo
      enddo
    enddo
!.....write
    open(ionum,file='out.NN.dgsf',status='replace',form='unformatted')
    do ia=1,natm
      do isf=1,nsf
        write(ionum) (dgsfo(1:3,jra,isf,ia),jra=1,natm)
      enddo
    enddo
    close(ionum)

    deallocate(dgsfo)
  end subroutine write_dgsf
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
        ku= 2*kd+kdd
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
!=======================================================================
  subroutine compute_stress(namax,natm,tag,ra,nnmax,strs,h &
       ,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nn,rc,lspr &
       ,mpi_world,myid)
    implicit none
    integer,intent(in):: namax,natm,nnmax,nb,nbmax,lsb(0:nbmax,6)&
         ,lsrc(6),myparity(3),nn(6),mpi_world,myid,lspr(0:nnmax,namax)&
         ,nex(3)
    real(8),intent(in):: ra(3,namax),tag(namax),h(3,3),rc,tcom
    real(8),intent(out):: strs(3,3,namax)

    integer:: ia,ja,ixyz,jxyz,ihl0,ihl1,ihl2,jj
    real(8):: xi(3),xj(3),xji(3),rij(3),rji(3),dji,sji,sii&
         ,hl2i,hl2j,tmp2i,tmp2j,hl1i,hl1j,tmp1i,tmp1j

    strs(1:3,1:3,1:namax) = 0d0
    if( nl.eq.1 ) then
      do ia=1,natm
        xi(1:3)= ra(1:3,ia)
        do jj=1,lspr(0,ia)
          ja= lspr(jj,ia)
          xj(1:3)= ra(1:3,ja)
          xji(1:3)= xj(1:3)-xi(1:3)
          rji(1:3)= h(1:3,1)*xji(1) +h(1:3,2)*xji(2) +h(1:3,3)*xji(3)
          rij(1:3)= -rji(1:3)
          dji= sqrt(rji(1)**2 +rji(2)**2 +rji(3)**2)
          if( dji.ge.rc ) cycle
          do ihl1=1,nhl(1)
            hl1i= hl1(ihl1,ia)
            tmp1i= wgt12(ihl1)*hl1i*(1d0-hl1i)
            do ihl0=1,nhl(0)
              do ixyz=1,3
                do jxyz=1,3
! derivative of gsf of atom-i by atom-j
                  sji= -tmp1i*wgt11(ihl0,ihl1)*dgsf(jxyz,ihl0,jj,ia) &
                       *rji(ixyz)
! counter contribution
                  sii= tmp1i*wgt11(ihl0,ihl1)*dgsf(jxyz,ihl0,jj,ia) &
                       *rij(ixyz)
                  strs(ixyz,jxyz,ja) = strs(ixyz,jxyz,ja) +sji
                  strs(ixyz,jxyz,ia) = strs(ixyz,jxyz,ia) +sii
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    else if( nl.eq.2 ) then
      do ia=1,natm
        xi(1:3)= ra(1:3,ia)
        do jj=1,lspr(0,ia)
          ja= lspr(jj,ia)
          xj(1:3)= ra(1:3,ja)
          xji(1:3)= xj(1:3)-xi(1:3)
          rji(1:3)= h(1:3,1)*xji(1) +h(1:3,2)*xji(2) +h(1:3,3)*xji(3)
          rij(1:3)= -rji(1:3)
          dji= sqrt(rji(1)**2 +rji(2)**2 +rji(3)**2)
          if( dji.ge.rc ) cycle
          do ihl2=1,nhl(2)
            hl2i= hl2(ihl2,ia)
            tmp2i= wgt23(ihl2) *hl2i*(1d0-hl2i)
            do ihl1=1,nhl(1)
              hl1i= hl1(ihl1,ia)
              tmp1i= wgt22(ihl1,ihl2) *hl1i*(1d0-hl1i)
              do ihl0=1,nhl(0)
!......derivative of gsf of atom-j by atom-i
                sji= -tmp2i *tmp1i &
                     *wgt21(ihl0,ihl1) *dgsf(jxyz,ihl0,jj,ia) &
                     *rji(ixyz)
!.....derivative of gsf of atom-i by atom-i
                sii= tmp2i *tmp1i &
                     *wgt21(ihl0,ihl1) *dgsf(jxyz,ihl0,jj,ia) &
                     *rij(ixyz)
                strs(ixyz,jxyz,ja) = strs(ixyz,jxyz,ja) +sji
                strs(ixyz,jxyz,ia) = strs(ixyz,jxyz,ia) +sii
              enddo
            enddo
          enddo
        enddo
      enddo
    endif

!-----send back (3-body)forces, stresses, and potentials on immigrants
    call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,strs,9)
    strs(1:3,1:3,1:natm) = strs(1:3,1:3,1:natm)*0.5d0
!!$    if( myid.ge.0 ) then
!!$      call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,lsrc,myparity &
!!$           ,nn,mpi_world,strs,9)
!!$    else
!!$      call reduce_dba_bk(natm,namax,tag,strs,9)
!!$    endif

  end subroutine compute_stress
end module NN
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
