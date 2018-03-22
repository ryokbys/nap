module NN
!-----------------------------------------------------------------------
!                     Last modified: <2018-03-22 17:53:34 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
!  Parallel implementation of neural-network potential with 1 hidden
!  layer. It is available for plural number of species.
!-----------------------------------------------------------------------
  save
  character(len=128):: paramsdir = '.'

!.....parameter file name
  character(128),parameter:: cpfname = 'in.params.NN'
  character(128),parameter:: ccfname = 'in.const.NN'

!.....logical flag for bias
  logical:: lbias = .false.
!.....logical flag for charge
  logical:: lcharge = .false.
!.....logical flag for electron temperature
  logical:: letemp = .false.

!.....Max num of species
  integer,parameter:: msp = 9
  logical:: interact(msp,msp)
  
!.....parameters
  integer:: nwgt1,nwgt2
  real(8),allocatable:: wgt11(:,:),wgt12(:)
  real(8),allocatable:: wgt21(:,:),wgt22(:,:),wgt23(:)
!.....constants
  integer,parameter:: nlmax= 2
  integer:: nsfc,nsfc1,nsfc2,nc1,nc2,nsp,nl
!.....number of nodes in each layer
!.....  nhl includes bias nodes whereas mhl does not
  integer:: nhl(0:nlmax+1),mhl(0:nlmax+1)
  integer,allocatable:: itype(:)
  real(8),allocatable:: cnst(:,:)
  real(8),allocatable:: hl1(:,:),hl2(:,:)
!.....symmetry function values and their derivatives
  real(8),allocatable:: gsf(:,:),dgsf(:,:,:,:)
!.....start and end points of symmetry function IDs for each species pair
  integer,allocatable:: iaddr2(:,:,:),iaddr3(:,:,:,:)
!.....symmetry function IDs for each pair
  integer(2),allocatable:: igsf(:,:,:)
!.....function types and num of constatns for types
  integer,parameter:: max_ncnst= 2
  integer:: ncnst_type(200)
  integer:: ncomb_type(200)

!.....max exponent of the basis function
  integer:: max_nexp

!.....cutoff region width ratio to rc
  real(8):: rcw = 0.9d0
  real(8):: rcnn = 4.0d0
  real(8):: rc3nn = 3.0d0

!.....num of atoms and neighbors for dgsf array
  integer:: nal, nnl, nalmax,nnlmax,nnltmp
  logical:: lrealloc = .false.

!.....parameters given from outside (fitpot)
  integer:: nprms
  real(8),allocatable:: prms(:)
  logical:: lprmset_NN = .false.

  
contains
  subroutine force_NN(namax,natm,tag,ra,nnmax,aa,strs,h,hi,tcom &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_world,myid,epi,epot,nismax,acon,lstrs,iprint,l1st)
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
    logical,intent(in):: l1st
    logical:: lstrs

!.....local
    integer:: i,j,k,l,m,n,is,js,ierr,ia,ja &
         ,ihl0,ihl1,ihl2,jj
    real(8):: at(3),epotl,epott,hl1i,hl2i,tmp2,tmp1,tmp
    real(8),allocatable,save:: strsl(:,:,:),aal(:,:)

!    real(8),allocatable:: aml(:,:,:,:),bml(:,:,:,:)

    if( l1st ) then
!!$      call read_const_NN(myid,mpi_world,iprint)
!!$      call read_params_NN(myid,iprint)
!.....reset rc
      if( myid.le.0 .and. rc .lt. rcnn- 1d-8) then
        write(6,'(a,f10.5,a,f10.5)') &
             ' Error: Cutoff radius rc should be corrected from '&
             ,rc,' to greater than ',rcnn
        if( myid.ge.0 ) then
          call mpi_finalize(ierr)
          stop
        else
          stop
        endif
      endif

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
      if( myid.le.0 .and. iprint.ne.0 ) then
        write(6,'(a,2i10)') ' max num of (local atoms *1.1) = ',nalmax
        write(6,'(a,2i10)') ' max num of (neighbors *1.1)   = ',nnlmax
        write(6,'(a,f10.3,a)') ' gsf size  = ', &
             dble(nhl(0)*nalmax*8)/1000/1000,' MB'
        write(6,'(a,f10.3,a)') ' dgsf size = ', &
             dble(3*nhl(0)*(nnlmax+1)*nalmax*8)/1000/1000,' MB'
        write(6,'(a,f10.3,a)') ' igsf size = ', &
             dble(nhl(0)*(nnlmax+1)*nalmax*2)/1000/1000,' MB'
      endif
      if( allocated(gsf) ) deallocate(gsf,dgsf,igsf)
      allocate( gsf(nhl(0),nal),dgsf(3,nhl(0),0:nnl,nal) &
           ,igsf(nhl(0),0:nnl,nal) )

      lrealloc = .false.
      if( nl.eq.1 ) then
        if( allocated(hl1) ) deallocate(hl1)
        allocate( hl1(nhl(1),nal) )
      else if( nl.eq.2 ) then
        if( allocated(hl1) ) deallocate(hl1,hl2)
        allocate( hl1(nhl(1),nal), hl2(nhl(2),nal) )
      endif

      if( allocated(strsl) ) deallocate(strsl)
      allocate(strsl(3,3,namax))
      if( allocated(aal) ) deallocate(aal)
      allocate(aal(3,namax))
    endif

    if( size(strsl).lt.3*3*namax ) then
      deallocate(strsl)
      allocate(strsl(3,3,namax))
    endif
    if( size(aal).lt.3*namax ) then
      deallocate(aal)
      allocate(aal(3,namax))
    endif
    strsl(1:3,1:3,1:namax) = 0d0
    aal(1:3,1:namax) = 0d0

!  Since natm and nn can change every step of MD,
!  if natm/nnltmp becomes nal/nnl, they should be updated and
!  gsf/dgsf as well.
    if( natm.gt.nal ) then
      nal = int(natm*1.1)
      if( nal .gt. namax ) then
        write(6,'(a)') ' [Error] nal .gt.namax'
        write(6,'(a,3i10)') '   myid,nal,namax = ',myid,nal,namax
        stop
      else
        write(6,*) ' Since natm.gt.nal, nal was updated at myid =',myid
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
      else
        write(6,*) ' Since nnltmp.gt.nnl, nnl was updated at myid =',myid
      endif
      lrealloc=.true.
    endif
    if( allocated(dgsf).and.lrealloc ) then
      deallocate( gsf,dgsf,igsf )
      allocate( gsf(nhl(0),nal),dgsf(3,nhl(0),0:nnl,nal) &
           ,igsf(nhl(0),0:nnl,nal))
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
         ,lspr,rcnn,rc3nn)
    if( lbias ) then
!.....set bias node to 1
      gsf(nhl(0),1:natm) = 1d0
      dgsf(1:3,nhl(0),:,:) = 0d0
    endif

    if( iprint.gt.10 .and. mod(iprint,10).eq.1 .and. myid.le.0 ) then
      open(80,file='out.NN.gsf',status='replace',form='unformatted')
      write(80) nhl(0)
      do ia=1,natm
        write(80) (gsf(ihl0,ia),ihl0=1,nhl(0))
      enddo
      close(80)
      call write_dgsf(84,natm,namax,nnmax,lspr,tag,nhl(0))
    endif

!.....initialize hidden-layer node values
    if( nl.eq.1 ) then
      hl1(1:nhl(1),1:natm)= 0d0
      if( lbias ) then
        hl1(nhl(1),1:natm)= 1d0
      endif
    else if( nl.eq.2 ) then
      hl1(1:nhl(1),1:natm)= 0d0
      hl2(1:nhl(2),1:natm)= 0d0
      if( lbias ) then
        hl1(nhl(1),1:natm)= 1d0
        hl2(nhl(2),1:natm)= 1d0
      endif
    endif

!.....2nd, calculate the node values by summing contributions from
!.....  symmetry functions
    if( nl.eq.1 ) then
      do ia=1,natm
!.....debug
        if( iprint.ge.100 ) then
          print *,"ia=",ia
          do ihl1=1,mhl(1)
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
        endif  ! end of debug
        do ihl1=1,mhl(1)
          tmp= 0d0
          do ihl0=1,nhl(0)
            tmp= tmp +wgt11(ihl0,ihl1) *gsf(ihl0,ia)
          enddo
          hl1(ihl1,ia)= sigmoid(tmp)
        enddo
      enddo
    else if( nl.eq.2 ) then
      do ia=1,natm
        do ihl1=1,mhl(1)
          tmp= 0d0
          do ihl0=1,nhl(0)
            tmp= tmp +wgt21(ihl0,ihl1) *gsf(ihl0,ia)
          enddo
          hl1(ihl1,ia)= sigmoid(tmp)
        enddo
        do ihl2=1,mhl(2)
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
        tmp = 0d0
        do ihl1=1,nhl(1)
          tmp = tmp +wgt12(ihl1) *(hl1(ihl1,ia)-0.5d0)
        enddo
        epi(ia) = epi(ia) +tmp
        epotl=epotl +tmp
      enddo
    else if( nl.eq.2 ) then
      do ia=1,natm
        tmp = 0d0
        do ihl2=1,nhl(2)
          tmp = tmp +wgt23(ihl2) *(hl2(ihl2,ia)-0.5d0)
        enddo
        epi(ia)= epi(ia) +tmp
        epotl=epotl +tmp
      enddo
    endif

!.....sum up for forces
    if( nl.eq.1 ) then
!.....loop over every energy per atom-i
      do ia=1,natm
        do ihl1=1,mhl(1)
          hl1i= hl1(ihl1,ia)
          tmp= wgt12(ihl1)*hl1i*(1d0-hl1i)
          do jj=1,lspr(0,ia)
            ja= lspr(jj,ia)
            js = int(tag(ja))
            do ihl0=1,nhl(0)
              if( igsf(ihl0,jj,ia).eq.0 ) cycle
              aal(1:3,ja)=aal(1:3,ja) &
                   -tmp*wgt11(ihl0,ihl1)*dgsf(1:3,ihl0,jj,ia)
            enddo
          enddo
          !.....atom ia
          do ihl0= 1,nhl(0)
            aal(1:3,ia)=aal(1:3,ia) &
                 -tmp*wgt11(ihl0,ihl1)*dgsf(1:3,ihl0,0,ia)
          enddo
        enddo
      enddo
    else if( nl.eq.2 ) then
      do ia=1,natm
        do ihl2=1,mhl(2)
          hl2i= hl2(ihl2,ia)
          tmp2= wgt23(ihl2) *hl2i*(1d0-hl2i)
          do ihl1=1,mhl(1)
            hl1i= hl1(ihl1,ia)
            tmp1= wgt22(ihl1,ihl2) *hl1i*(1d0-hl1i)
            do jj=1,lspr(0,ia)
              ja= lspr(jj,ia)
              do ihl0=1,nhl(0)
                if( igsf(ihl0,jj,ia).eq.0 ) cycle
                aal(1:3,ja)=aal(1:3,ja) &
                     -tmp2 *tmp1 &
                     *wgt21(ihl0,ihl1)*dgsf(1:3,ihl0,jj,ia)
              enddo
            enddo
!.....atom ia
            do ihl0= 1,nhl(0)
              aal(1:3,ia)=aal(1:3,ia) &
                   -tmp2*tmp1*wgt21(ihl0,ihl1)*dgsf(1:3,ihl0,0,ia)
            enddo
          enddo
        enddo
      enddo
    endif

    call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,aal,3)
    aa(1:3,1:natm) = aa(1:3,1:natm) +aal(1:3,1:natm)

    if( lstrs ) then
      call compute_stress(namax,natm,tag,ra,nnmax,strsl,h &
           ,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nn,rcnn,lspr &
           ,mpi_world,myid)
      strs(1:3,1:3,1:natm) = strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)
    endif

!-----gather epot
    if( myid.ge.0 ) then
      call mpi_allreduce(epotl,epott,1,mpi_double_precision &
           ,mpi_sum,mpi_world,ierr)
      if( iprint.gt.2 ) print *,'NN epot = ',epott
      epot= epot +epott
    else
      epot= epot +epotl
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

    integer:: isf,ia,jj,ja,kk,ka,is,js,ks
    real(8):: xi(3),xj(3),xij(3),rij(3),dij,fcij,eta,rs,texp,driji(3), &
         dfcij,drijj(3),dgdr,xk(3),xik(3),rik(3),dik,fcik,dfcik, &
         driki(3),drikk(3),almbd,spijk,cs,t1,t2,dgdij,dgdik,dgcs, &
         dcsdj(3),dcsdk(3),dcsdi(3),tcos,tpoly,a1,a2,tmorse,rc2

    rc2 = rc*rc
    gsf(1:nsf,1:nal)= 0d0
    dgsf(1:3,1:nsf,0:nnl,1:nal)= 0d0
    igsf(1:nsf,0:nnl,1:nal) = 0
    do ia=1,natm
      xi(1:3)= ra(1:3,ia)
      is= int(tag(ia))
      do jj=1,lspr(0,ia)
        ja= lspr(jj,ia)
        if( ja.eq.ia ) cycle
        xj(1:3)= ra(1:3,ja)
        xij(1:3)= xj(1:3)-xi(1:3)
        rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij= rij(1)**2 +rij(2)**2 +rij(3)**2
        if( dij.ge.rc2 ) cycle
        dij = sqrt(dij)
        js= int(tag(ja))
        driji(1:3)= -rij(1:3)/dij
        drijj(1:3)= -driji(1:3)
        fcij= fc(dij,rc)
        dfcij= dfc(dij,rc)
        do isf=iaddr2(1,is,js),iaddr2(2,is,js)
!!$          print *,'ia,is,ja,js,isf,itype=',ia,is,ja,js,isf,itype(isf)
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
            igsf(isf,0,ia) = 1
            igsf(isf,jj,ia) = 1
          else if( itype(isf).eq.2 ) then ! cosine
            a1= cnst(1,isf)
            !.....func value
            tcos= (1d0+cos(dij*a1))
            gsf(isf,ia)= gsf(isf,ia) +tcos*fcij
            !.....derivative
            dgdr= -a1*sin(dij*a1)*fcij +tcos*dfc(dij,rc)
            dgsf(1:3,isf,0,ia)= dgsf(1:3,isf,0,ia) +driji(1:3)*dgdr
            dgsf(1:3,isf,jj,ia)= dgsf(1:3,isf,jj,ia) +drijj(1:3)*dgdr
            igsf(isf,0,ia) = 1
            igsf(isf,jj,ia) = 1
          else if( itype(isf).eq.3 ) then ! polynomial
            a1= cnst(1,isf)
            !.....func value
            tpoly= 1d0*dij**(-a1)
            gsf(isf,ia)= gsf(isf,ia) +tpoly*fcij
            !.....derivative
            dgdr= -a1*dij**(-a1-1d0)*fcij +tpoly*dfc(dij,rc)
            dgsf(1:3,isf,0,ia)= dgsf(1:3,isf,0,ia) +driji(1:3)*dgdr
            dgsf(1:3,isf,jj,ia)= dgsf(1:3,isf,jj,ia) +drijj(1:3)*dgdr
            igsf(isf,0,ia) = 1
            igsf(isf,jj,ia) = 1
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
            igsf(isf,0,ia) = 1
            igsf(isf,jj,ia) = 1
          endif
        enddo

        if( dij.gt.rc3 ) cycle
        fcij= fc(dij,rc3)
        dfcij= dfc(dij,rc3)
        do kk=1,lspr(0,ia)
          ka= lspr(kk,ia)
          ks= int(tag(ka))
!!$          if( .not.interact(is,ks) ) cycle
          if( iaddr3(1,is,js,ks).lt.0 ) cycle
          if( ka.eq.ia .or. ka.le.ja ) cycle
          xk(1:3)= ra(1:3,ka)
          xik(1:3)= xk(1:3)-xi(1:3)
          rik(1:3)= h(1:3,1)*xik(1) +h(1:3,2)*xik(2) +h(1:3,3)*xik(3)
          dik= sqrt(rik(1)**2 +rik(2)**2 +rik(3)**2)
          if( dik.ge.rc3 ) cycle
          do isf=iaddr3(1,is,js,ks),iaddr3(2,is,js,ks)
            almbd= cnst(1,isf)
            t2= (abs(almbd)+1d0)**2
            fcik= fc(dik,rc3)
            dfcik= dfc(dik,rc3)
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
            dcsdj(1:3)= rik(1:3)/dij/dik -rij(1:3)*cs/dij**2
            dcsdk(1:3)= rij(1:3)/dij/dik -rik(1:3)*cs/dik**2
            dcsdi(1:3)= -dcsdj(1:3) -dcsdk(1:3)
            dgsf(1:3,isf,0,ia)= dgsf(1:3,isf,0,ia) +dgcs*dcsdi(1:3)
            dgsf(1:3,isf,jj,ia)= dgsf(1:3,isf,jj,ia) +dgcs*dcsdj(1:3)
            dgsf(1:3,isf,kk,ia)= dgsf(1:3,isf,kk,ia) +dgcs*dcsdk(1:3)
            igsf(isf,0,ia) = 1
            igsf(isf,jj,ia) = 1
            igsf(isf,kk,ia) = 1
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
  subroutine read_params(myid,mpi_world,iprint)
    implicit none
    include 'mpif.h'

    integer,intent(in):: myid,mpi_world,iprint

    integer:: ierr,i,j,k,nc,ncoeff &
         ,ihl0,ihl1,ihl2,icmb(3),nsf,nsf1,nsf2,iap,jap,kap,ndat
    integer,allocatable:: nwgt(:)
    logical:: lexist
    character:: ctmp*128,fname*128
    integer,external:: num_data

!.....initialize some
    ncnst_type(1)= 2   ! Gaussian
    ncnst_type(2)= 1   ! cosine
    ncnst_type(3)= 1   ! polynomial
    ncnst_type(4)= 2   ! Morse
    ncnst_type(101)= 1 ! angular

    ncomb_type(1:100)= 2    ! pair
    ncomb_type(101:200)= 3  ! triplet

!.....read constants at the 1st call
    fname = trim(paramsdir)//'/'//trim(ccfname)
    inquire(file=trim(fname),exist=lexist)
    if( .not. lexist ) then
      if( myid.ge.0 ) then
        if( myid.eq.0 ) then
          write(6,'(a)') ' [Error] '//trim(fname)//' does not exist !!!.'
          write(6,'(a)') '   The NN potential needs '//trim(fname)//'.'
        endif
        call mpi_finalize(ierr)
        stop
      else
        write(6,'(a)') ' [Error] '//trim(fname)//' does not exist !!!.'
        write(6,'(a)') '   The NN potential needs '//trim(fname)//'.'
        stop
      endif
    endif
    open(51,file=trim(fname),status='old')
!.....num of symmetry functions, num of node in 1st hidden layer
10  read(51,'(a)') ctmp
    if( ctmp(1:1).eq.'!' .or. ctmp(1:1).eq.'#' ) then
      call parse_option(ctmp,iprint,ierr)
      goto 10
    else
      backspace(51)
    endif
    read(51,*) nl,nsp,(nhl(i),i=0,nl)
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

    call mpi_bcast(lbias,1,mpi_logical,0,mpi_world,ierr)
    call mpi_bcast(lcharge,1,mpi_logical,0,mpi_world,ierr)
    call mpi_bcast(letemp,1,mpi_logical,0,mpi_world,ierr)

!.....Determine num of symmetry functions and nodes
    nsf= nhl(0)
    nhl(nl+1)= 1  ! only one output node, an energy
    mhl(0:nl+1)= nhl(0:nl+1)
    if( lbias ) then  ! bias node
      nhl(0:nl) = nhl(0:nl) +1
    endif
    if( letemp ) nhl(0) = nhl(0) +1  ! T_e
    if( myid.eq.0 .and. iprint.ne.0 ) then
      print *,'lbias  = ',lbias
      print *,'lcharge= ',lcharge
      print *,'letemp = ',letemp
      print *,'nhl = ',nhl(0:nl+1)
      print *,'mhl = ',mhl(0:nl+1)
    endif

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
      if(myid.eq.0 ) then
        print *,'[Error] nsf.ne.nsf1+nsf2 !!!'
      endif
      call mpi_finalize(ierr)
      stop
    endif
    close(51)

!.....read parameters at the 1st call
    fname = trim(paramsdir)//'/'//trim(cpfname)
    inquire(file=trim(fname),exist=lexist)
    if( .not. lexist ) then
      if( myid.ge.0 ) then
        if( myid.eq.0 ) then
          write(6,'(a)') ' [Error] '//trim(fname)//' does not exist !!!.'
          write(6,'(a)') '   The NN potential needs '//trim(fname)//'.'
        endif
        call mpi_finalize(ierr)
        stop
      else
        write(6,'(a)') ' [Error] '//trim(fname)//' does not exist !!!.'
        write(6,'(a)') '   The NN potential needs '//trim(fname)//'.'
        stop
      endif
    endif
    open(50,file=trim(fname),status='old')
    read(50,*) ncoeff,rcnn,rc3nn
!.....check whether the num of parameters is correct
    if( rc3nn.gt.rcnn ) then
      rc3nn= rcnn
      if( myid.le.0 .and. iprint.ne.0 ) then
        write(6,*) ' rc3nn was corrected to rcnn = ',rcnn
        write(6,*) ' because input rc3nn > rcnn, which should not happen.'
      endif
    endif

!.....calc number of weights
    allocate(nwgt(nl+1))
    nwgt(1:nl+1) = 0
    do i=1,nl+1
      nwgt(i)= nhl(i-1)*mhl(i)
    enddo
    if( myid.le.0 .and. iprint.ne.0 ) then
      print *, 'num of basis funcs =',nhl(0)
      do i=1,nl
        print *, 'ihl, nhl(ihl)  =',i,nhl(i)
      enddo
      do i=1,nl+1
        print *, 'ihl, nwgt(ihl)  =',i,nwgt(i)
      enddo
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
!.....different number of weights and number of layers
    if( nl.eq.1 ) then
      allocate(wgt11(nhl(0),mhl(1)),wgt12(nhl(1)))
      wgt11(1:nhl(0),1:mhl(1)) = 0d0
      wgt12(1:nhl(1)) = 0d0
    else if( nl.eq.2 ) then
      allocate(wgt21(nhl(0),mhl(1)),wgt22(nhl(1),mhl(2)),wgt23(nhl(2)))
      wgt21(1:nhl(0),1:mhl(1)) = 0d0
      wgt22(1:nhl(1),1:mhl(2)) = 0d0
      wgt23(1:nhl(2)) = 0d0
    endif
    if( nl.eq.1 ) then
      do ihl0=1,nhl(0)
        do ihl1=1,mhl(1)
          read(50,*) wgt11(ihl0,ihl1)
        enddo
      enddo
      do ihl1=1,nhl(1)
        read(50,*) wgt12(ihl1)
      enddo
    else if( nl.eq.2 ) then
      do ihl0=1,nhl(0)
        do ihl1=1,mhl(1)
          read(50,*) wgt21(ihl0,ihl1)
        enddo
      enddo
      do ihl1=1,nhl(1)
        do ihl2=1,mhl(2)
          read(50,*) wgt22(ihl1,ihl2)
        enddo
      enddo
      do ihl2=1,nhl(2)
        read(50,*) wgt23(ihl2)
      enddo
    endif
    close(50)


    deallocate(nwgt)
    return
  end subroutine read_params
!=======================================================================
  subroutine read_const_NN(myid,mpi_world,iprint)
    implicit none
    include 'mpif.h'

    integer,intent(in):: myid,mpi_world,iprint

    integer:: ierr,i,j,k,nc,ncoeff &
         ,ihl0,ihl1,ihl2,icmb(3),nsf,nsf1,nsf2,iap,jap,kap,ndat
    logical:: lexist
    character:: ctmp*128,fname*128
    integer,external:: num_data

!.....initialize some
    ncnst_type(1)= 2   ! Gaussian
    ncnst_type(2)= 1   ! cosine
    ncnst_type(3)= 1   ! polynomial
    ncnst_type(4)= 2   ! Morse
    ncnst_type(101)= 1 ! angular

    ncomb_type(1:100)= 2    ! pair
    ncomb_type(101:200)= 3  ! triplet

!.....read constants at the 1st call
    fname = trim(paramsdir)//'/'//trim(ccfname)
    inquire(file=trim(fname),exist=lexist)
    if( .not. lexist ) then
      if( myid.ge.0 ) then
        if( myid.eq.0 ) then
          write(6,'(a)') ' [Error] '//trim(fname)//' does not exist !!!.'
          write(6,'(a)') '   The NN potential needs '//trim(fname)//'.'
        endif
        call mpi_finalize(ierr)
        stop
      else
        write(6,'(a)') ' [Error] '//trim(fname)//' does not exist !!!.'
        write(6,'(a)') '   The NN potential needs '//trim(fname)//'.'
        stop
      endif
    endif
    open(51,file=trim(fname),status='old')
!.....num of symmetry functions, num of node in 1st hidden layer
10  read(51,'(a)') ctmp
    if( ctmp(1:1).eq.'!' .or. ctmp(1:1).eq.'#' ) then
      call parse_option(ctmp,iprint,ierr)
      goto 10
    else
      backspace(51)
    endif
    read(51,*) nl,nsp,(nhl(i),i=0,nl)
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

    call mpi_bcast(lbias,1,mpi_logical,0,mpi_world,ierr)
    call mpi_bcast(lcharge,1,mpi_logical,0,mpi_world,ierr)
    call mpi_bcast(letemp,1,mpi_logical,0,mpi_world,ierr)

!.....Determine num of symmetry functions and nodes
    nsf= nhl(0)
    nhl(nl+1)= 1  ! only one output node, an energy
    mhl(0:nl+1)= nhl(0:nl+1)
    if( lbias ) then  ! bias node
      nhl(0:nl) = nhl(0:nl) +1
    endif
    if( letemp ) nhl(0) = nhl(0) +1  ! T_e
    if( myid.eq.0 .and. iprint.ne.0 ) then
      print *,'lbias  = ',lbias
      print *,'lcharge= ',lcharge
      print *,'letemp = ',letemp
      print *,'nhl = ',nhl(0:nl+1)
      print *,'mhl = ',mhl(0:nl+1)
    endif

    if( .not.allocated(itype) ) then
      allocate(itype(nsf),cnst(max_ncnst,nsf))
      allocate(iaddr2(2,nsp,nsp),iaddr3(2,nsp,nsp,nsp))
    endif
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
      if(myid.eq.0 ) then
        print *,'[Error] nsf.ne.nsf1+nsf2 !!!'
      endif
      call mpi_finalize(ierr)
      stop
    endif
    close(51)

    return
  end subroutine read_const_NN
!=======================================================================
  subroutine read_params_NN(myid,iprint)
    implicit none
    include 'mpif.h'

    integer,intent(in):: myid,iprint

    integer:: ierr,i,j,k,nc,ncoeff &
         ,ihl0,ihl1,ihl2,icmb(3),nsf,nsf1,nsf2,iap,jap,kap,ndat
    integer,allocatable:: nwgt(:)
    logical:: lexist
    character:: ctmp*128,fname*128
    integer,external:: num_data

!.....read parameters at the 1st call
    fname = trim(paramsdir)//'/'//trim(cpfname)
    inquire(file=trim(fname),exist=lexist)
    if( .not. lexist ) then
      if( myid.ge.0 ) then
        if( myid.eq.0 ) then
          write(6,'(a)') ' [Error] '//trim(fname)//' does not exist !!!.'
          write(6,'(a)') '   The NN potential needs '//trim(fname)//'.'
        endif
        call mpi_finalize(ierr)
        stop
      else
        write(6,'(a)') ' [Error] '//trim(fname)//' does not exist !!!.'
        write(6,'(a)') '   The NN potential needs '//trim(fname)//'.'
        stop
      endif
    endif
    open(50,file=trim(fname),status='old')
    read(50,*) ncoeff,rcnn,rc3nn
!.....check whether the num of parameters is correct
    if( rc3nn.gt.rcnn ) then
      rc3nn= rcnn
      if( myid.le.0 .and. iprint.ne.0 ) then
        write(6,*) ' rc3nn was corrected to rcnn = ',rcnn
        write(6,*) ' because input rc3nn > rcnn, which should not happen.'
      endif
    endif

!.....calc number of weights
    allocate(nwgt(nl+1))
    nwgt(1:nl+1) = 0
    do i=1,nl+1
      nwgt(i)= nhl(i-1)*mhl(i)
    enddo
    if( myid.le.0 .and. iprint.ne.0 ) then
      print *, 'num of basis funcs =',nhl(0)
      do i=1,nl
        print *, 'ihl, nhl(ihl)  =',i,nhl(i)
      enddo
      do i=1,nl+1
        print *, 'ihl, nwgt(ihl)  =',i,nwgt(i)
      enddo
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
!.....different number of weights and number of layers
    if( nl.eq.1 ) then
      allocate(wgt11(nhl(0),mhl(1)),wgt12(nhl(1)))
      wgt11(1:nhl(0),1:mhl(1)) = 0d0
      wgt12(1:nhl(1)) = 0d0
    else if( nl.eq.2 ) then
      allocate(wgt21(nhl(0),mhl(1)),wgt22(nhl(1),mhl(2)),wgt23(nhl(2)))
      wgt21(1:nhl(0),1:mhl(1)) = 0d0
      wgt22(1:nhl(1),1:mhl(2)) = 0d0
      wgt23(1:nhl(2)) = 0d0
    endif
    if( nl.eq.1 ) then
      do ihl0=1,nhl(0)
        do ihl1=1,mhl(1)
          read(50,*) wgt11(ihl0,ihl1)
        enddo
      enddo
      do ihl1=1,nhl(1)
        read(50,*) wgt12(ihl1)
      enddo
    else if( nl.eq.2 ) then
      do ihl0=1,nhl(0)
        do ihl1=1,mhl(1)
          read(50,*) wgt21(ihl0,ihl1)
        enddo
      enddo
      do ihl1=1,nhl(1)
        do ihl2=1,mhl(2)
          read(50,*) wgt22(ihl1,ihl2)
        enddo
      enddo
      do ihl2=1,nhl(2)
        read(50,*) wgt23(ihl2)
      enddo
    endif
    close(50)

    deallocate(nwgt)
    return
  end subroutine read_params_NN
!=======================================================================
  subroutine set_params_NN(nprms_in,prms_in,rc,rc3)
!
!  Accessor routine to set NN parameters from outside.
!  Curretnly this routine is supposed to be called only on serial run.
!
    implicit none 
    integer,intent(in):: nprms_in
    real(8),intent(in):: prms_in(nprms_in),rc,rc3

    nprms = nprms_in
    if( .not.allocated(prms) ) allocate(prms(nprms))
    prms(1:nprms) = prms_in(1:nprms_in)

    rcnn = rc
    rc3nn= rc3
    
    lprmset_NN = .true.
    return
  end subroutine set_params_NN
!=======================================================================
  subroutine check_consistency(myid,iprint)
    implicit none
    integer,intent(in):: myid,iprint

    integer:: i,nc
    integer,allocatable:: nwgt(:)

!.....calc number of weights
    allocate(nwgt(nl+1))
    nwgt(1:nl+1) = 0
    do i=1,nl+1
      nwgt(i)= nhl(i-1)*mhl(i)
    enddo
    
    nc= 0
    do i=1,nl+1
      nc= nc +nwgt(i)
    enddo
    if( nprms .ne. nc ) then
      write(6,'(a)') ' [Error] Number of parameters is not correct !!!'
      write(6,'(a,i0)')  '   nprms = ',nprms
      write(6,'(a,i0)')  '   where nprms should be ',nc
      stop
    endif

    deallocate(nwgt)
    return
  end subroutine check_consistency
!=======================================================================
  subroutine update_params_NN()
!
!  Update NN parameters by taking parameter values from params array.
!  This routine would be called only from externally from fitpot.
!
    integer:: i,inc,is,js,ihl0,ihl1,ihl2

    if( .not.lprmset_NN ) then
      print *,'ERROR: params have not been set yet.'
      stop
    endif
    
    if( nl.eq.1 ) then
      if( .not.allocated(wgt11) ) &
           allocate(wgt11(nhl(0),mhl(1)),wgt12(nhl(1)))
      wgt11(1:nhl(0),1:mhl(1)) = 0d0
      wgt12(1:nhl(1)) = 0d0
    else if( nl.eq.2 ) then
      if( .not.allocated(wgt21) ) &
           allocate(wgt21(nhl(0),mhl(1)), &
           wgt22(nhl(1),mhl(2)), &
           wgt23(nhl(2)))
      wgt21(1:nhl(0),1:mhl(1)) = 0d0
      wgt22(1:nhl(1),1:mhl(2)) = 0d0
      wgt23(1:nhl(2)) = 0d0
    endif

    inc = 0
    if( nl.eq.1 ) then
      do ihl0=1,nhl(0)
        do ihl1=1,mhl(1)
          inc = inc +1
          wgt11(ihl0,ihl1) = prms(inc)
        enddo
      enddo
      do ihl1=1,nhl(1)
        inc = inc +1
        wgt12(ihl1) = prms(inc)
      enddo
    else if( nl.eq.2 ) then
      do ihl0=1,nhl(0)
        do ihl1=1,mhl(1)
          inc = inc +1
          wgt21(ihl0,ihl1) = prms(inc)
        enddo
      enddo
      do ihl1=1,nhl(1)
        do ihl2=1,mhl(2)
          inc = inc +1
          wgt22(ihl1,ihl2) = prms(inc)
        enddo
      enddo
      do ihl2=1,nhl(2)
        inc = inc +1
        wgt23(ihl2) = prms(inc)
      enddo
    endif
    
    return
  end subroutine update_params_NN
!=======================================================================
  subroutine parse_option(cline,iprint,ierr)
!
!  Parse options from a comment line.
!  Lines starting from ! or # are treated as comment lines,
!  but options can be given at the comment lines.
!  The option words should be put after these comment characters with
!  one or more spaces between them for example,
!
!  bias:  .true.
!
!  Currently available options are:
!    - "bias:" with an argument .true. (T) or .false. (F)
!    - "charge:" with an argument .true. (T) or .false. (F)
!
    implicit none
    character(len=*),intent(in):: cline
    integer,intent(in):: iprint
    integer,intent(out):: ierr

    character(len=10):: c1,copt
    logical:: lopt
    integer,external:: num_data

    ierr = 0
    if( index(cline,'bias:').ne.0 ) then
      read(cline,*) c1,copt,lopt
      if( trim(copt).ne.'bias:' ) then
        print *, 'Error: copt is not "bias:" !!!'
        ierr = 1
      endif
      lbias = lopt
    else if( index(cline,'charge:').ne.0 ) then
      read(cline,*) c1,copt,lopt
      if( trim(copt).ne.'charge:' ) then
        print *, 'Error: copt is not "charge:" !!!'
        ierr = 2
      endif
      lcharge = lopt
    endif
    
  end subroutine parse_option
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
    integer(2),allocatable:: igsfo(:,:,:)
    integer,external:: itotOf

    allocate(dgsfo(3,natm,nsf,natm),igsfo(nsf,natm,natm))
!.....reduce d(i)gsf data of buffer atoms to those of resident atoms
    dgsfo(1:3,1:natm,1:nsf,1:natm)= 0d0
    igsfo(1:nsf,1:natm,1:natm)= 0
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
          igsfo(isf,jra,ia) = igsf(isf,jj,ia)
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
!.....write igsf to ionum+1
    open(ionum+1,file='out.NN.igsf',status='replace',form='unformatted')
    do ia=1,natm
      do ja=1,natm
        write(ionum+1) (igsfo(isf,ja,ia),isf=1,nsf)
      enddo
    enddo
    close(ionum+1)

    deallocate(dgsfo,igsfo)
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
    real(8),intent(in):: ra(3,namax),tag(namax),h(3,3),rc
    real(8),intent(inout):: tcom
    real(8),intent(out):: strs(3,3,namax)

    integer:: ia,ja,ixyz,jxyz,ihl0,ihl1,ihl2,jj,is,js
    real(8):: xi(3),xj(3),xji(3),rij(3),rji(3),dji,sji,sii&
         ,hl2i,tmp2i,hl1i,tmp1i,stmp(3,3)

    if( nl.eq.1 ) then
      do ia=1,natm
        xi(1:3)= ra(1:3,ia)
        is = int(tag(ia))
        stmp(1:3,1:3)= 0d0
        do jj=1,lspr(0,ia)
          ja= lspr(jj,ia)
          xj(1:3)= ra(1:3,ja)
          xji(1:3)= xj(1:3)-xi(1:3)
          rji(1:3)= h(1:3,1)*xji(1) +h(1:3,2)*xji(2) +h(1:3,3)*xji(3)
          rij(1:3)= -rji(1:3)
          dji= sqrt(rji(1)**2 +rji(2)**2 +rji(3)**2)
          if( dji.ge.rc ) cycle
          js = int(tag(ja))
          do ihl1=1,mhl(1)
            hl1i= hl1(ihl1,ia)
            tmp1i= wgt12(ihl1)*hl1i*(1d0-hl1i)
!!$            do ihl0=iaddr2(1,is,js),iaddr2(2,is,js)
            do ihl0=1,nhl(0)
              if( igsf(ihl0,jj,ia).eq.0 ) cycle
              do ixyz=1,3
                do jxyz=1,3
! derivative of gsf of atom-i by atom-j
                  sji= -tmp1i*wgt11(ihl0,ihl1)*dgsf(jxyz,ihl0,jj,ia) &
                       *rji(ixyz)
!  Since counter contribution seems the same as the above,
!  no need to perform redundant calculations.
!!$! counter contribution
!!$                  sii= tmp1i*wgt11(ihl0,ihl1)*dgsf(jxyz,ihl0,jj,ia) &
!!$                       *rij(ixyz)
                  strs(ixyz,jxyz,ja) = strs(ixyz,jxyz,ja) +sji
                  strs(ixyz,jxyz,ia) = strs(ixyz,jxyz,ia) +sji
                  stmp(ixyz,jxyz) = stmp(ixyz,jxyz) +sji
                enddo  ! jxyz
              enddo  ! ixyz
            enddo  ! ihl0
          enddo  ! ihl1
        enddo  ! ja
      enddo  ! ia
    else if( nl.eq.2 ) then
      do ia=1,natm
        xi(1:3)= ra(1:3,ia)
        is = int(tag(ia))
        do jj=1,lspr(0,ia)
          ja= lspr(jj,ia)
          xj(1:3)= ra(1:3,ja)
          xji(1:3)= xj(1:3)-xi(1:3)
          rji(1:3)= h(1:3,1)*xji(1) +h(1:3,2)*xji(2) +h(1:3,3)*xji(3)
          rij(1:3)= -rji(1:3)
          dji= sqrt(rji(1)**2 +rji(2)**2 +rji(3)**2)
          if( dji.ge.rc ) cycle
          js = int(tag(ja))
          do ihl2=1,mhl(2)
            hl2i= hl2(ihl2,ia)
            tmp2i= wgt23(ihl2) *hl2i*(1d0-hl2i)
            do ihl1=1,mhl(1)
              hl1i= hl1(ihl1,ia)
              tmp1i= wgt22(ihl1,ihl2) *hl1i*(1d0-hl1i)
!!$              do ihl0=iaddr2(1,is,js),iaddr2(2,is,js)
              do ihl0=1,nhl(0)
                if( igsf(ihl0,jj,ia).eq.0 ) cycle
                do ixyz=1,3
                  do jxyz=1,3
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
        enddo
      enddo
    endif

!-----send back (3-body)forces, stresses, and potentials on immigrants
    call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,strs,9)
    strs(1:3,1:3,1:natm) = strs(1:3,1:3,1:natm)*0.5d0
    return
  end subroutine compute_stress
!=======================================================================
  subroutine set_paramsdir_NN(dname)
!
!  Accessor routine for setting paramsdir
!
    implicit none
    character(len=*),intent(in):: dname

    paramsdir = trim(dname)
    return
  end subroutine set_paramsdir_NN
!=======================================================================
  subroutine pderiv_NN()
!
!  Derivative w.r.t. NN parameters, {w}
!
    implicit none
    
  end subroutine pderiv_NN
end module NN
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
