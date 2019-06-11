module NN2
!-----------------------------------------------------------------------
!                     Last modified: <2019-06-11 16:39:48 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
!  Parallel implementation of neural-network potential with upto 2
!  hidden layers. It is available for plural number of species.
!
!  To separate the symmetry function calculations in descriptor.F90,
!  the module was re-created as NN2.
!-----------------------------------------------------------------------
  implicit none 
  save
  character(len=128):: paramsdir = '.'

!.....parameter file name
  character(128),parameter:: cpfname = 'in.params.NN2'

  real(8),parameter:: pi= 3.14159265358979d0

  integer:: mem
  real(8):: time
  
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
  integer:: nwgt1,nwgt2,nwtot
  real(8),allocatable:: wgt11(:,:),wgt12(:)
  real(8),allocatable:: wgt21(:,:),wgt22(:,:),wgt23(:)
!.....constants
  integer,parameter:: nlmax= 2
  integer:: nsfc,nsfc1,nsfc2,nc1,nc2,nsp,nl
!.....number of nodes in each layer
!.....  nhl includes bias nodes whereas mhl does not
  integer:: nhl(0:nlmax+1),mhl(0:nlmax+1)
  real(8),allocatable:: hl1(:,:),hl2(:,:)
  real(8),allocatable:: zl1(:,:),zl2(:,:)

!.....parameters given from outside (fitpot)
  integer:: nprms
  real(8),allocatable:: prms(:)
  logical:: lprmset_NN2 = .false.

!.....Sigmoid function types:
!     1) 1/(1+exp(-x))
!     2) 1/(1+exp(-x)) +asig*x
  integer:: itypesig = 1
!.....Coefficient of additional term in sigmoid2
  real(8):: asig = 0.01d0
  
contains
  subroutine force_NN2(namax,natm,tag,ra,nnmax,aa,strs,h,hi,tcom &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rcin,lspr &
       ,mpi_world,myid,epi,epot,nismax,lstrs,iprint,l1st)
    use descriptor,only: gsf,dgsf,igsf,nsf,nal,calc_desc,make_gsf_arrays
    use util,only: itotOf
    implicit none
    include "mpif.h"
    include "./params_unit.h"
    integer,intent(in):: namax,natm,nnmax,nismax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_world,myid,lspr(0:nnmax,namax),nex(3)
    real(8),intent(in):: ra(3,namax),tag(namax) &
         ,h(3,3),hi(3,3),sv(3,6)
    real(8),intent(inout):: tcom,rcin
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: l1st
    logical:: lstrs

!.....local
    integer:: i,j,k,l,m,n,is,js,ierr,ia,ja &
         ,ihl0,ihl1,ihl2,jj
    real(8):: at(3),epotl,epott,hl1i,hl2i,tmp2,tmp1,tmp,zl1i,zl2i,time0
    real(8),allocatable,save:: strsl(:,:,:),aal(:,:)

    integer:: itot
!!$    integer,external:: itotOf
    character(len=8):: cnum

    call make_gsf_arrays(l1st,namax,natm &
         ,tag,nnmax,lspr,myid,mpi_world,iprint)

    if( l1st ) then

      if( nl.eq.1 ) then
        if( allocated(hl1) ) deallocate(hl1,zl1)
        allocate( hl1(nhl(1),nal),zl1(nhl(1),nal) )
      else if( nl.eq.2 ) then
        if( allocated(hl1) ) deallocate(hl1,hl2,zl1,zl2)
        allocate( hl1(nhl(1),nal), hl2(nhl(2),nal) &
             ,zl1(nhl(1),nal), zl2(nhl(2),nal) )
      endif

      if( allocated(strsl) ) deallocate(strsl,aal)
      allocate(strsl(3,3,namax),aal(3,namax))

    endif ! l1st

    if( allocated(hl1) .and. size(hl1).lt.nhl(1)*nal ) then
      deallocate(hl1,zl1)
      allocate(hl1(nhl(1),nal), zl1(nhl(1),nal))
    endif
    if( allocated(hl2) .and. size(hl2).lt.nhl(1)*nal ) then
      deallocate(hl1,hl2,zl1,zl2)
      allocate(hl1(nhl(1),nal), hl2(nhl(2),nal) &
           ,zl1(nhl(1),nal), zl2(nhl(2),nal))
    endif
    if( size(strsl).ne.3*3*namax ) then
      deallocate(strsl,aal)
      allocate(strsl(3,3,namax),aal(3,namax))
    endif

    strsl(1:3,1:3,1:namax) = 0d0
    aal(1:3,1:namax) = 0d0

    call calc_desc(namax,natm,nb,nnmax,h &
         ,tag,ra,lspr,rcin,myid,mpi_world,l1st,iprint)

    time0 = mpi_wtime()

    if( lbias ) then
!.....set bias node to 1
      gsf(nhl(0),1:natm) = 1d0
      dgsf(1:3,nhl(0),:,:) = 0d0
    endif


!.....initialize hidden-layer node values
    if( nl.eq.1 ) then
      hl1(1:nhl(1),1:natm)= 0d0
      zl1(:,:) = 0d0
      if( lbias ) then
        hl1(nhl(1),1:natm)= 1d0
        zl1(nhl(1),1:natm)= 1d0
      endif
    else if( nl.eq.2 ) then
      hl1(1:nhl(1),1:natm)= 0d0
      hl2(1:nhl(2),1:natm)= 0d0
      zl1(:,:) = 0d0
      zl2(:,:) = 0d0
      if( lbias ) then
        hl1(nhl(1),1:natm)= 1d0
        hl2(nhl(2),1:natm)= 1d0
      endif
    endif
!.....2nd, calculate the node values by summing contributions from
!.....  symmetry functions
    if( nl.eq.1 ) then
      do ia=1,natm
        do ihl1=1,mhl(1)
          tmp= 0d0
          do ihl0=1,nhl(0)
            tmp= tmp +wgt11(ihl0,ihl1) *gsf(ihl0,ia)
          enddo
          zl1(ihl1,ia)= tmp
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
          zl1(ihl1,ia)= tmp
          hl1(ihl1,ia)= sigmoid(tmp)
        enddo
        do ihl2=1,mhl(2)
          tmp= 0d0
          do ihl1=1,nhl(1)
            tmp= tmp +wgt22(ihl1,ihl2) *(hl1(ihl1,ia)-0.5d0)
          enddo
          zl2(ihl2,ia)= tmp
          hl2(ihl2,ia)= sigmoid(tmp)
        enddo
      enddo
    endif

!.....Calculate the energy of atom by summing up the node values
    epotl= 0d0
    if( nl.eq.1 ) then
      do ia=1,natm
        tmp = 0d0
        do ihl1=1,nhl(1)
!!$          tmp = tmp +wgt12(ihl1) *(hl1(ihl1,ia)-0.5d0)
          tmp = tmp +wgt12(ihl1) *(sigmoid(zl1(ihl1,ia))-0.5d0)
        enddo
        epi(ia) = epi(ia) +tmp
        epotl=epotl +tmp
      enddo
    else if( nl.eq.2 ) then
      do ia=1,natm
        tmp = 0d0
        do ihl2=1,nhl(2)
!!$          tmp = tmp +wgt23(ihl2) *(hl2(ihl2,ia)-0.5d0)
          tmp = tmp +wgt23(ihl2) *(sigmoid(zl2(ihl2,ia))-0.5d0)
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
          zl1i= zl1(ihl1,ia)
!!$          tmp= wgt12(ihl1)*hl1i*(1d0-hl1i)
          tmp = wgt12(ihl1) *dsigmoid(zl1i,hl1i)
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
          zl2i= zl2(ihl2,ia)
!!$          tmp2= wgt23(ihl2) *hl2i*(1d0-hl2i)
          tmp2= wgt23(ihl2) *dsigmoid(zl2i,hl2i)
          do ihl1=1,mhl(1)
            hl1i= hl1(ihl1,ia)
            zl1i= zl1(ihl1,ia)
!!$            tmp1= wgt22(ihl1,ihl2) *hl1i*(1d0-hl1i)
            tmp1= wgt22(ihl1,ihl2) *dsigmoid(zl1i,hl1i)
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
           ,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nn,rcin,lspr &
           ,mpi_world,myid)
      strs(1:3,1:3,1:natm) = strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)
    endif

!-----gather epot
    call mpi_allreduce(epotl,epott,1,mpi_double_precision &
         ,mpi_sum,mpi_world,ierr)
    if( iprint.gt.2 ) print *,'NN2 epot = ',epott
    epot= epot +epott

    time = time +(mpi_wtime() -time0)

    return
  end subroutine force_NN2
!=======================================================================
  subroutine force_NN2_overlay_pot(namax,natm,tag,ra,nnmax,aa,strs &
       ,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rcin,lspr &
       ,mpi_world,myid,epi,epot,nismax,lstrs,iprint,l1st)
!
!  NN2 potential overlaying with ZBL by potential
!
    use descriptor,only: gsf,dgsf,igsf,nsf,nal,calc_desc,make_gsf_arrays
    use util,only: itotOf
    use force,only: ol_ranges,ol_alphas,ol_dalphas,ol_pair,ol_type
    implicit none
    include "mpif.h"
    include "./params_unit.h"
    integer,intent(in):: namax,natm,nnmax,nismax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_world,myid,lspr(0:nnmax,namax),nex(3)
    real(8),intent(in):: ra(3,namax),tag(namax) &
         ,h(3,3),hi(3,3),sv(3,6)
    real(8),intent(inout):: tcom,rcin
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: l1st
    logical:: lstrs

!.....local
    integer:: i,j,k,l,m,n,is,js,ierr,ia,ja &
         ,ihl0,ihl1,ihl2,jj,ixyz,jxyz
    real(8):: at(3),epotl,epott,hl1i,hl2i,tmp2,tmp1,tmp,zl1i,zl2i,time0
    real(8),allocatable,save:: strsl(:,:,:),aal(:,:),epit(:)
    real(8):: xij(3),rij(3),dij2,dij,xi(3),alpi,drdxi(3),dtmp,ri,ro

    integer:: itot
!!$    integer,external:: itotOf
    character(len=8):: cnum

    call make_gsf_arrays(l1st,namax,natm &
         ,tag,nnmax,lspr,myid,mpi_world,iprint)

    if( l1st ) then

      if( nl.eq.1 ) then
        if( allocated(hl1) ) deallocate(hl1,zl1)
        allocate( hl1(nhl(1),nal),zl1(nhl(1),nal) )
      else if( nl.eq.2 ) then
        if( allocated(hl1) ) deallocate(hl1,hl2,zl1,zl2)
        allocate( hl1(nhl(1),nal), hl2(nhl(2),nal) &
             ,zl1(nhl(1),nal), zl2(nhl(2),nal) )
      endif

      if( allocated(strsl) ) deallocate(strsl,aal,epit)
      allocate(strsl(3,3,namax),aal(3,namax),epit(namax))

    endif ! l1st

    if( allocated(hl1) .and. size(hl1).lt.nhl(1)*nal ) then
      deallocate(hl1,zl1)
      allocate(hl1(nhl(1),nal), zl1(nhl(1),nal))
    endif
    if( allocated(hl2) .and. size(hl2).lt.nhl(1)*nal ) then
      deallocate(hl1,hl2,zl1,zl2)
      allocate(hl1(nhl(1),nal), hl2(nhl(2),nal) &
           ,zl1(nhl(1),nal), zl2(nhl(2),nal))
    endif
    if( size(strsl).ne.3*3*namax ) then
      deallocate(strsl,aal,epit)
      allocate(strsl(3,3,namax),aal(3,namax),epit(namax))
    endif

    strsl(1:3,1:3,1:natm+nb) = 0d0
    aal(1:3,1:natm+nb) = 0d0
    epit(1:natm+nb) = 0d0

!!$    return

    call calc_desc(namax,natm,nb,nnmax,h &
         ,tag,ra,lspr,rcin,myid,mpi_world,l1st,iprint)

    time0 = mpi_wtime()

    if( lbias ) then
!.....set bias node to 1
      gsf(nhl(0),1:natm) = 1d0
      dgsf(1:3,nhl(0),:,:) = 0d0
    endif


!.....initialize hidden-layer node values
    if( nl.eq.1 ) then
      hl1(1:nhl(1),1:natm)= 0d0
      zl1(:,:) = 0d0
      if( lbias ) then
        hl1(nhl(1),1:natm)= 1d0
        zl1(nhl(1),1:natm)= 1d0
      endif
    else if( nl.eq.2 ) then
      hl1(1:nhl(1),1:natm)= 0d0
      hl2(1:nhl(2),1:natm)= 0d0
      zl1(:,:) = 0d0
      zl2(:,:) = 0d0
      if( lbias ) then
        hl1(nhl(1),1:natm)= 1d0
        hl2(nhl(2),1:natm)= 1d0
      endif
    endif
!.....2nd, calculate the node values by summing contributions from
!.....  symmetry functions
    if( nl.eq.1 ) then
      do ia=1,natm
        do ihl1=1,mhl(1)
          tmp= 0d0
          do ihl0=1,nhl(0)
            tmp= tmp +wgt11(ihl0,ihl1) *gsf(ihl0,ia)
          enddo
          zl1(ihl1,ia)= tmp
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
          zl1(ihl1,ia)= tmp
          hl1(ihl1,ia)= sigmoid(tmp)
        enddo
        do ihl2=1,mhl(2)
          tmp= 0d0
          do ihl1=1,nhl(1)
            tmp= tmp +wgt22(ihl1,ihl2) *(hl1(ihl1,ia)-0.5d0)
          enddo
          zl2(ihl2,ia)= tmp
          hl2(ihl2,ia)= sigmoid(tmp)
        enddo
      enddo
    endif

!.....Calculate the energy of atom by summing up the node values
    epotl= 0d0
    if( nl.eq.1 ) then
      do ia=1,natm
        tmp = 0d0
        alpi = ol_alphas(0,ia)
        do ihl1=1,nhl(1)
!!$          tmp = tmp +wgt12(ihl1) *(hl1(ihl1,ia)-0.5d0)
          tmp = tmp +wgt12(ihl1) *(sigmoid(zl1(ihl1,ia))-0.5d0)
        enddo
        epit(ia) = tmp
        epi(ia) = epi(ia) +tmp *alpi
        epotl=epotl +tmp *alpi
      enddo
    else if( nl.eq.2 ) then
      do ia=1,natm
        tmp = 0d0
        alpi = ol_alphas(0,ia)
        do ihl2=1,nhl(2)
!!$          tmp = tmp +wgt23(ihl2) *(hl2(ihl2,ia)-0.5d0)
          tmp = tmp +wgt23(ihl2) *(sigmoid(zl2(ihl2,ia))-0.5d0)
        enddo
        epit(ia)= tmp
        epi(ia) = epi(ia) +tmp*alpi
        epotl=epotl +tmp*alpi
      enddo
    endif

!.....sum up for forces
    if( nl.eq.1 ) then
!.....loop over every energy per atom-i
      do ia=1,natm
        alpi = ol_alphas(0,ia)
        do ihl1=1,mhl(1)
          hl1i= hl1(ihl1,ia)
          zl1i= zl1(ihl1,ia)
!!$          tmp= wgt12(ihl1)*hl1i*(1d0-hl1i)
          tmp = wgt12(ihl1) *dsigmoid(zl1i,hl1i) *alpi
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
        alpi = ol_alphas(0,ia)
        do ihl2=1,mhl(2)
          hl2i= hl2(ihl2,ia)
          zl2i= zl2(ihl2,ia)
!!$          tmp2= wgt23(ihl2) *hl2i*(1d0-hl2i)
          tmp2= wgt23(ihl2) *dsigmoid(zl2i,hl2i) *alpi
          do ihl1=1,mhl(1)
            hl1i= hl1(ihl1,ia)
            zl1i= zl1(ihl1,ia)
!!$            tmp1= wgt22(ihl1,ihl2) *hl1i*(1d0-hl1i)
            tmp1= wgt22(ihl1,ihl2) *dsigmoid(zl1i,hl1i)
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

!.....Derivative wrt.alphas, which requires energy per atom, epit
    do ia=1,natm
      is = int(tag(ia))
      xi(1:3) = ra(1:3,ia)
      alpi = ol_alphas(0,ia)
      do jj=1,lspr(0,ia)
        ja = lspr(jj,ia)
        js = int(tag(ja))
        ri = ol_pair(1,is,js)
        ro = ol_pair(2,is,js)
        xij(1:3) = ra(1:3,ja) -xi(1:3)
        rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij2 = rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
        if( dij2.le.ri*ri .or. dij2.gt.ro*ro ) cycle
        dij = dsqrt(dij2)
        drdxi(1:3)= -rij(1:3)/dij
!.....Force
        dtmp = epit(ia) *alpi/ol_alphas(jj,ia) *ol_dalphas(jj,ia)
        aal(1:3,ia)= aal(1:3,ia) -dtmp*drdxi(1:3)
        aal(1:3,ja)= aal(1:3,ja) +dtmp*drdxi(1:3)
!.....Stress
        if( .not.lstrs ) cycle
        do ixyz=1,3
          do jxyz=1,3
            strsl(jxyz,ixyz,i)=strsl(jxyz,ixyz,i) &
                 -0.5d0*dtmp*rij(ixyz)*(-drdxi(jxyz))
            strsl(jxyz,ixyz,j)=strsl(jxyz,ixyz,j) &
                 -0.5d0*dtmp*rij(ixyz)*(-drdxi(jxyz))
          enddo
        enddo
      enddo
    enddo

    call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex &
         ,lsrc,myparity,nn,mpi_world,aal,3)
    aa(1:3,1:natm) = aa(1:3,1:natm) +aal(1:3,1:natm)

    if( lstrs ) then
      call compute_stress(namax,natm,tag,ra,nnmax,strsl,h &
           ,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nn,rcin,lspr &
           ,mpi_world,myid)
      strs(1:3,1:3,1:natm) = strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)
    endif

!!$    epi(1:natm) = epi(1:natm) +epit(1:natm)

!-----gather epot
    call mpi_allreduce(epotl,epott,1,mpi_double_precision &
         ,mpi_sum,mpi_world,ierr)
    if( iprint.gt.2 ) print *,'epot NN2_overlay = ',epott
    epot= epot +epott

    time = time +(mpi_wtime() -time0)

    return
  end subroutine force_NN2_overlay_pot
!=======================================================================
  subroutine force_NN2_overlay_frc(namax,natm,tag,ra,nnmax,aa,strs &
       ,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rcin,lspr &
       ,mpi_world,myid,epi,epot,nismax,lstrs,iprint,l1st)
!
!  NN2 potential overlaying with ZBL by force
!
    use descriptor,only: gsf,dgsf,igsf,nsf,nal,calc_desc,make_gsf_arrays &
         ,rcmax,rcmax2,descs
    use util,only: itotOf
    use force,only: ol_ranges,ol_alphas,ol_dalphas,ol_pair,ol_type
    use ZBL,only: dvnucl
    implicit none
    include "mpif.h"
    include "./params_unit.h"
    integer,intent(in):: namax,natm,nnmax,nismax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_world,myid,lspr(0:nnmax,namax),nex(3)
    real(8),intent(in):: ra(3,namax),tag(namax) &
         ,h(3,3),hi(3,3),sv(3,6)
    real(8),intent(inout):: tcom,rcin
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: l1st
    logical:: lstrs

!.....local
    integer:: i,j,k,l,m,n,is,js,ierr,ia,ja &
         ,ihl0,ihl1,ihl2,jj,ixyz,jxyz
    real(8):: at(3),epotl,epott,hl1i,hl2i,tmp2,tmp1,tmp,zl1i,zl2i
    real(8),allocatable,save:: strsl(:,:,:),aal(:,:)
    real(8):: xij(3),rij(3),dij2,dij,xi(3),xj(3),alpij,drdxi(3) &
         ,dtmp,ri,ro,tmpij(3),time0,sij

    integer:: itot
!!$    integer,external:: itotOf
    character(len=8):: cnum

    call make_gsf_arrays(l1st,namax,natm &
         ,tag,nnmax,lspr,myid,mpi_world,iprint)

    if( l1st ) then

      if( nl.eq.1 ) then
        if( allocated(hl1) ) deallocate(hl1,zl1)
        allocate( hl1(nhl(1),nal),zl1(nhl(1),nal) )
      else if( nl.eq.2 ) then
        if( allocated(hl1) ) deallocate(hl1,hl2,zl1,zl2)
        allocate( hl1(nhl(1),nal), hl2(nhl(2),nal) &
             ,zl1(nhl(1),nal), zl2(nhl(2),nal) )
      endif

      if( allocated(strsl) ) deallocate(strsl,aal)
      allocate(strsl(3,3,namax),aal(3,namax))

    endif ! l1st

    if( allocated(hl1) .and. size(hl1).lt.nhl(1)*nal ) then
      deallocate(hl1,zl1)
      allocate(hl1(nhl(1),nal), zl1(nhl(1),nal))
    endif
    if( allocated(hl2) .and. size(hl2).lt.nhl(1)*nal ) then
      deallocate(hl1,hl2,zl1,zl2)
      allocate(hl1(nhl(1),nal), hl2(nhl(2),nal) &
           ,zl1(nhl(1),nal), zl2(nhl(2),nal))
    endif
    if( size(strsl).ne.3*3*namax ) then
      deallocate(strsl,aal)
      allocate(strsl(3,3,namax),aal(3,namax))
    endif

    strsl(1:3,1:3,1:natm+nb) = 0d0
    aal(1:3,1:natm+nb) = 0d0

!!$    return

    call calc_desc(namax,natm,nb,nnmax,h &
         ,tag,ra,lspr,rcin,myid,mpi_world,l1st,iprint)

    time0 = mpi_wtime()

    if( lbias ) then
!.....set bias node to 1
      gsf(nhl(0),1:natm) = 1d0
      dgsf(1:3,nhl(0),:,:) = 0d0
    endif


!.....initialize hidden-layer node values
    if( nl.eq.1 ) then
      hl1(1:nhl(1),1:natm)= 0d0
      zl1(:,:) = 0d0
      if( lbias ) then
        hl1(nhl(1),1:natm)= 1d0
        zl1(nhl(1),1:natm)= 1d0
      endif
    else if( nl.eq.2 ) then
      hl1(1:nhl(1),1:natm)= 0d0
      hl2(1:nhl(2),1:natm)= 0d0
      zl1(:,:) = 0d0
      zl2(:,:) = 0d0
      if( lbias ) then
        hl1(nhl(1),1:natm)= 1d0
        hl2(nhl(2),1:natm)= 1d0
      endif
    endif
!.....2nd, calculate the node values by summing contributions from
!.....  symmetry functions
    if( nl.eq.1 ) then
      do ia=1,natm
        do ihl1=1,mhl(1)
          tmp= 0d0
          do ihl0=1,nhl(0)
            tmp= tmp +wgt11(ihl0,ihl1) *gsf(ihl0,ia)
          enddo
          zl1(ihl1,ia)= tmp
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
          zl1(ihl1,ia)= tmp
          hl1(ihl1,ia)= sigmoid(tmp)
        enddo
        do ihl2=1,mhl(2)
          tmp= 0d0
          do ihl1=1,nhl(1)
            tmp= tmp +wgt22(ihl1,ihl2) *(hl1(ihl1,ia)-0.5d0)
          enddo
          zl2(ihl2,ia)= tmp
          hl2(ihl2,ia)= sigmoid(tmp)
        enddo
      enddo
    endif

!.....Calculate the energy of atom by summing up the node values
    epotl= 0d0
    if( nl.eq.1 ) then
      do ia=1,natm
        tmp = 0d0
        do ihl1=1,nhl(1)
          tmp = tmp +wgt12(ihl1) *(sigmoid(zl1(ihl1,ia))-0.5d0)
        enddo
        epi(ia) = epi(ia) +tmp
        epotl=epotl +tmp
      enddo
    else if( nl.eq.2 ) then
      do ia=1,natm
        tmp = 0d0
        do ihl2=1,nhl(2)
          tmp = tmp +wgt23(ihl2) *(sigmoid(zl2(ihl2,ia))-0.5d0)
        enddo
        epi(ia) = epi(ia) +tmp
        epotl=epotl +tmp
      enddo
    endif

!.....sum up for forces
    if( nl.eq.1 ) then
!.....loop over every energy per atom-i
      do ia=1,natm
        is = int(tag(ia))
        xi(1:3) = ra(1:3,ia)
        do jj=1,lspr(0,ia)
          ja= lspr(jj,ia)
          js = int(tag(ja))
          xj(1:3) = ra(1:3,ja)
          xij(1:3)= xj(1:3) -xi(1:3)
          rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
          dij2 = rij(1)**2 +rij(2)**2 +rij(3)**2
          if( dij2.gt.rcmax2 ) cycle
          dij = dsqrt(dij2)
          ri = ol_pair(1,is,js)
          ro = ol_pair(2,is,js)
          alpij = ol_alphas(jj,ia)
          if( dij.ge.ro ) then
            tmpij(1:3) = 0d0
            do ihl1=1,mhl(1)
              hl1i= hl1(ihl1,ia)
              zl1i= zl1(ihl1,ia)
              tmp = wgt12(ihl1) *dsigmoid(zl1i,hl1i)
              do ihl0=1,nhl(0)
                if( dij.ge.descs(ihl0)%rcut ) cycle
                if( igsf(ihl0,jj,ia).eq.0 ) cycle
                tmpij(1:3)= tmpij(1:3) +tmp*wgt11(ihl0,ihl1)*dgsf(1:3,ihl0,jj,ia)
              enddo
            enddo
          else if( dij.ge.ri .and. dij.lt.ro ) then
            drdxi(1:3) = -rij(1:3)/dij
            tmpij(1:3)= -0.5d0*dvnucl(is,js,dij)*drdxi(1:3)*(1d0-alpij)
            do ihl1=1,mhl(1)
              hl1i= hl1(ihl1,ia)
              zl1i= zl1(ihl1,ia)
              tmp = wgt12(ihl1) *dsigmoid(zl1i,hl1i)
              do ihl0=1,nhl(0)
                if( dij.ge.descs(ihl0)%rcut ) cycle
                if( igsf(ihl0,jj,ia).eq.0 ) cycle
                tmpij(1:3)=tmpij(1:3) +alpij*tmp*wgt11(ihl0,ihl1)*dgsf(1:3,ihl0,jj,ia)
              enddo
            enddo
          else if( dij.lt.ri ) then
            drdxi(1:3) = -rij(1:3)/dij
            tmpij(1:3) = -0.5d0*dvnucl(is,js,dij)*drdxi(1:3)
          endif
          aal(1:3,ja)=aal(1:3,ja) -tmpij(1:3)
          aal(1:3,ia)=aal(1:3,ia) +tmpij(1:3)
          do ixyz=1,3
            do jxyz=1,3
              sij = -tmpij(jxyz)*rij(ixyz)
              strsl(ixyz,jxyz,ja) = strsl(ixyz,jxyz,ja) +sij
              strsl(ixyz,jxyz,ia) = strsl(ixyz,jxyz,ia) +sij
            enddo
          enddo
!!$!.....atom ia
!!$          do ihl0= 1,nhl(0)
!!$            aal(1:3,ia)=aal(1:3,ia) &
!!$                 -tmp*wgt11(ihl0,ihl1)*dgsf(1:3,ihl0,0,ia)
!!$          enddo
        enddo
      enddo
    else if( nl.eq.2 ) then
      do ia=1,natm
        is = int(tag(ia))
        xi(1:3) = ra(1:3,ia)
        do ihl2=1,mhl(2)
          hl2i= hl2(ihl2,ia)
          zl2i= zl2(ihl2,ia)
!!$          tmp2= wgt23(ihl2) *hl2i*(1d0-hl2i)
          tmp2= wgt23(ihl2) *dsigmoid(zl2i,hl2i)
          do ihl1=1,mhl(1)
            hl1i= hl1(ihl1,ia)
            zl1i= zl1(ihl1,ia)
!!$            tmp1= wgt22(ihl1,ihl2) *hl1i*(1d0-hl1i)
            tmp1= wgt22(ihl1,ihl2) *dsigmoid(zl1i,hl1i)
            do jj=1,lspr(0,ia)
              ja= lspr(jj,ia)
              js = int(tag(ja))
              xj(1:3) = ra(1:3,ja)
              xij(1:3)= xj(1:3) -xi(1:3)
              rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
              dij2 = rij(1)**2 +rij(2)**2 +rij(3)**2
              if( dij2.gt.rcmax2 ) cycle
              dij = dsqrt(dij2)
              ri = ol_pair(1,is,js)
              ro = ol_pair(2,is,js)
              alpij = ol_alphas(jj,ia)
              if( dij.ge.ro ) then
                do ihl0=1,nhl(0)
                  if( dij.ge.descs(ihl0)%rcut ) cycle
                  if( igsf(ihl0,jj,ia).eq.0 ) cycle
                  aal(1:3,ja)=aal(1:3,ja) -tmp2 *tmp1 *wgt21(ihl0,ihl1)*dgsf(1:3,ihl0,jj,ia)
                  aal(1:3,ia)=aal(1:3,ia) +tmp2 *tmp1 *wgt21(ihl0,ihl1)*dgsf(1:3,ihl0,jj,ia)
                enddo
              else if( dij.ge.ri .and. dij.lt.ro ) then
                drdxi(1:3) = -xij(1:3)/dij
                tmp1 = dvnucl(is,js,dij)
                aal(1:3,ia)=aal(1:3,ia) -tmp1*drdxi(1:3)*(1d0-alpij)
                aal(1:3,ja)=aal(1:3,ja) +tmp1*drdxi(1:3)*(1d0-alpij)
                do ihl0=1,nhl(0)
                  if( dij.ge.descs(ihl0)%rcut ) cycle
                  if( igsf(ihl0,jj,ia).eq.0 ) cycle
                  aal(1:3,ja)=aal(1:3,ja) -alpij *tmp2 *tmp1 *wgt21(ihl0,ihl1)*dgsf(1:3,ihl0,jj,ia)
                  aal(1:3,ia)=aal(1:3,ia) +alpij *tmp2 *tmp1 *wgt21(ihl0,ihl1)*dgsf(1:3,ihl0,jj,ia)
                enddo
              else if( dij.lt.ri ) then
                drdxi(1:3) = -xij(1:3)/dij
                tmp1 = dvnucl(is,js,dij)
                aal(1:3,ia)=aal(1:3,ia) -tmp1*drdxi(1:3)
                aal(1:3,ja)=aal(1:3,ja) +tmp1*drdxi(1:3)
              endif
            enddo
!!$!.....atom ia
!!$            do ihl0= 1,nhl(0)
!!$              aal(1:3,ia)=aal(1:3,ia) &
!!$                   -tmp2*tmp1*wgt21(ihl0,ihl1)*dgsf(1:3,ihl0,0,ia)
!!$            enddo
          enddo
        enddo
      enddo
    endif

    call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex &
         ,lsrc,myparity,nn,mpi_world,aal,3)
    aa(1:3,1:natm) = aa(1:3,1:natm) +aal(1:3,1:natm)

!-----send back stresses
    call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex &
         ,lsrc,myparity,nn,mpi_world,strsl,9)
    strs(:,:,1:natm) = strs(:,:,1:natm) +strsl(:,:,1:natm)*0.5d0

!!$    if( lstrs ) then
!!$      call compute_stress(namax,natm,tag,ra,nnmax,strsl,h &
!!$           ,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nn,rcin,lspr &
!!$           ,mpi_world,myid)
!!$      strs(1:3,1:3,1:natm) = strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)
!!$    endif

!-----gather epot
    call mpi_allreduce(epotl,epott,1,mpi_real8,mpi_sum,mpi_world,ierr)
    if( iprint.gt.2 ) print *,'epot NN2_overlay = ',epott
    epot= epot +epott

    time = time +(mpi_wtime() -time0)

    return
  end subroutine force_NN2_overlay_frc
!=======================================================================
  subroutine set_sigtype_NN2(itype)
!
!  Set sigmoid function type
!
    integer,intent(in):: itype

    itypesig = itype
    return
  end subroutine set_sigtype_NN2
!=======================================================================
  function sigmoid(x)
    real(8),intent(in):: x
    real(8):: sigmoid

    select case(itypesig)
    case(1)
      sigmoid = 1d0/(1d0 +exp(-x))
    case(2)
      sigmoid = 1d0/(1d0 +exp(-x)) +asig*x
    case default
      sigmoid = 0d0
    end select
    return
  end function sigmoid
!=======================================================================
  function dsigmoid(x,sx)
    real(8),intent(in):: x,sx
    real(8):: dsigmoid
    real(8):: sxt
    
    select case(itypesig)
    case(1)
      dsigmoid = sx*(1d0-sx)
    case(2)
      sxt = 1d0/(1d0 +exp(-x))
      dsigmoid = sxt *(1d0-sxt) +asig
    case default
      dsigmoid = 0d0
    end select
    return
  end function dsigmoid
!=======================================================================
  function ddsigmoid(x,sx)
    real(8),intent(in):: x,sx
    real(8):: ddsigmoid
    real(8):: sxt
    
    select case(itypesig)
    case(1)
      ddsigmoid = sx*(1d0-sx)*(1d0-2d0*sx)
    case(2)
      sxt = 1d0 /(1d0 +exp(-x))
      ddsigmoid = sxt *(1d0 -sxt) *(1d0 -2d0*sxt)
    case default
      ddsigmoid = 0d0
    end select
    return
  end function ddsigmoid
!=======================================================================
  subroutine read_params_NN2(myid,mpi_world,iprint)
!-----------------------------------------------------------------------
!  Assume that the descriptor information is already read.
!  Input file format is as follows:
!-----------------------------------------------------------------------
!  ! comments or options
!  NL, NHL(0:NL)
!  WGT(i),  LOWER(i),  UPPER(i)
!  ...
!-----------------------------------------------------------------------
!  - NL: number of hidden-layer (1 or 2)
!  - NHL: number of nodes in 0th, 1st, and 2nd layers
!  - WGT(i): weight of the i-th edge
!  - LOWER,UPPER: lower and upper limit of the weight
!-----------------------------------------------------------------------
    use descriptor,only: nsf,iglid
    use util, only: num_data
    implicit none
    include 'mpif.h'

    integer,intent(in):: myid,mpi_world,iprint

    integer:: ierr,i,j,k,nc &
         ,ihl0,ihl1,ihl2,icmb(3)
    integer,allocatable:: nwgt(:)
    logical:: lexist
    character:: ctmp*128,fname*128
!!$    integer,external:: num_data

!.....Check if the file exists
    fname = trim(paramsdir)//'/'//trim(cpfname)
    inquire(file=trim(fname),exist=lexist)
    if( .not. lexist ) then
      if( myid.eq.0 ) then
        write(6,'(a)') ' [Error] '//trim(fname)//' does not exist !!!.'
        write(6,'(a)') '   The NN potential needs '//trim(fname)//'.'
      endif
      call mpi_finalize(ierr)
      stop
    endif

!.....Read parameters at the 1st call
    if( myid.eq.0 ) then
      open(50,file=trim(fname),status='old')
!.....num of symmetry functions, num of node in 1st hidden layer
10    read(50,'(a)') ctmp
      if( ctmp(1:1).eq.'!' .or. ctmp(1:1).eq.'#' ) then
        call parse_option(ctmp,iprint,ierr)
        goto 10
      else
        backspace(50)
      endif
      read(50,*) nl, (nhl(i),i=0,nl)
    endif

    call mpi_bcast(nl,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(nhl,nlmax+2,mpi_integer,0,mpi_world,ierr)

    if( nhl(0).ne.nsf ) then
      if( myid.eq.0 ) then
        print *,'ERROR: nhl(0).ne.nsf !'
        print *,'  Check the consistency in in.params.desc and in.params.NN'
      endif
      call mpi_finalize(ierr)
      stop
    endif

    call mpi_bcast(lbias,1,mpi_logical,0,mpi_world,ierr)
    call mpi_bcast(lcharge,1,mpi_logical,0,mpi_world,ierr)
    call mpi_bcast(letemp,1,mpi_logical,0,mpi_world,ierr)
    call mpi_bcast(itypesig,1,mpi_integer,0,mpi_world,ierr)

!.....Determine num of weights
    nhl(nl+1) = 1
    mhl(0:nl+1) = nhl(0:nl+1)
    if( lbias ) then
      nhl(0:nl) = nhl(0:nl) + 1
    endif
    if( letemp ) nhl(0) = nhl(0) + 1
    if( myid.eq.0 .and. iprint.ne.0 ) then
      print *,''
      print *,'NN potential parameters:'
      print *,'  lbias  = ',lbias
      print *,'  lcharge= ',lcharge
      print *,'  letemp = ',letemp
      print *,'  nhl = ',nhl(0:nl+1)
      print *,'  mhl = ',mhl(0:nl+1)
      if( itypesig.eq.1 ) then
        print *,'  activation function: 1) 1/(1+exp(-x))'
      else if( itypesig.eq.2 ) then
        print *,'  activation function: 2) 1/(1+exp(-x)) +asig*x'
        print '(a,f7.4)','                          asig = ',asig
      else
        print *,'  activation function: unknown'
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
    nwtot = nc

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

    if( myid.eq.0 ) then
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
    endif
    if( nl.eq.1 ) then
      call mpi_bcast(wgt11,nhl(0)*mhl(1),mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(wgt12,nhl(1),mpi_real8,0,mpi_world,ierr)
    else if( nl.eq.2 ) then
      call mpi_bcast(wgt21,nhl(0)*mhl(1),mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(wgt22,nhl(1)*mhl(2),mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(wgt23,nhl(2),mpi_real8,0,mpi_world,ierr)
    endif

!.....Allocate Group-LASSO/FS related variable, which is not used in pmd
    if( .not. allocated(iglid) ) allocate(iglid(nwtot))

    deallocate(nwgt)
    return
  end subroutine read_params_NN2
!=======================================================================
  subroutine set_params_NN2(nprms_in,prms_in,nl_in,nhl_in)
!
!  Accessor routine to set NN parameters from outside.
!  Curretnly this routine is supposed to be called only on serial run.
!
    use descriptor,only: iglid
    implicit none 
    integer,intent(in):: nprms_in,nl_in,nhl_in(0:nl_in)
    real(8),intent(in):: prms_in(nprms_in)

    integer:: i

    nl = nl_in
    if( nl.eq.0 ) then
      print *,'ERROR: nl==0 which should not happen.'
      print *,'  Probably NN_num_layers and NN_num_nodes are not set in in.fitpot.'
      print *,'  Those should be consistent with the number of NN weigts to be optimized.'
      stop
    endif
    nhl(0:nl) = nhl_in(0:nl_in)
    nhl(nl+1) = 1
    mhl(0:nl+1) = nhl(0:nl+1)

    nwtot = 0
    do i=1,nl+1
      nwtot = nwtot +mhl(i)*nhl(i-1)
    enddo

    if( nwtot.ne.nprms_in ) then
      print *,'ERROR: nl_in,nhl_in,nprms_in not consistent !!'
      print *,'  Check in.vars.fitpot or in.params.NN2 and'&
           //' NN_num_nodes parameters in in.fitpot.'
      stop
    endif
    
    nprms = nprms_in
    if( .not.allocated(prms) ) allocate(prms(nprms))
    prms(1:nprms) = prms_in(1:nprms_in)

!.....Allocate Group-LASSO/FS related variable, which is not used in pmd
    if( .not. allocated(iglid) .or. size(iglid).ne.nwtot ) &
         allocate(iglid(nwtot))

    lprmset_NN2 = .true.
    return
  end subroutine set_params_NN2
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
  subroutine update_params_NN2()
!
!  Update NN parameters by taking parameter values from params array.
!  This routine would be called only from externally from fitpot.
!
    integer:: i,inc,is,js,ihl0,ihl1,ihl2

    if( .not.lprmset_NN2 ) then
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
  end subroutine update_params_NN2
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
!    - "sigtype:" sigmoid type: 1 or 2.
!
    use util, only: num_data
    implicit none
    character(len=*),intent(in):: cline
    integer,intent(in):: iprint
    integer,intent(out):: ierr

    real(8):: ropt
    character(len=10):: c1,copt
    logical:: lopt
!!$    integer,external:: num_data
    integer:: iopt

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
    else if( index(cline,'sigtype:').ne.0 ) then
      read(cline,*) c1,copt,iopt
      if( trim(copt).ne.'sigtype:' ) then
        print *, 'Error: copt is not "sigtype:" !!!'
        ierr = 2
      endif
      itypesig = iopt
    endif
    
  end subroutine parse_option
!=======================================================================
  subroutine compute_stress(namax,natm,tag,ra,nnmax,strs,h &
       ,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nn,rc,lspr &
       ,mpi_world,myid)
    use descriptor,only: igsf,dgsf
    implicit none
    integer,intent(in):: namax,natm,nnmax,nb,nbmax,lsb(0:nbmax,6)&
         ,lsrc(6),myparity(3),nn(6),mpi_world,myid,lspr(0:nnmax,namax)&
         ,nex(3)
    real(8),intent(in):: ra(3,namax),tag(namax),h(3,3),rc
    real(8),intent(inout):: tcom
    real(8),intent(out):: strs(3,3,namax)

    integer:: ia,ja,ixyz,jxyz,ihl0,ihl1,ihl2,jj,is,js
    real(8):: xi(3),xj(3),xji(3),rij(3),rji(3),dji,dji2,sji,sii&
         ,hl2i,tmp2i,hl1i,tmp1i,stmp(3,3),zl1i,zl2i

    real(8),save:: rcmax2
    logical,save:: l1st = .true.

    if( l1st ) then
      rcmax2 = rc*rc
      l1st = .false.
    endif

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
          dji2= rji(1)**2 +rji(2)**2 +rji(3)**2
          if( dji2.ge.rcmax2 ) cycle
          dji = dsqrt(dji2)
          js = int(tag(ja))
          do ihl1=1,mhl(1)
            hl1i= hl1(ihl1,ia)
            zl1i= zl1(ihl1,ia)
!!$            tmp1i= wgt12(ihl1)*hl1i*(1d0-hl1i)
            tmp1i= wgt12(ihl1)* dsigmoid(zl1i,hl1i)
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
          dji2= rji(1)**2 +rji(2)**2 +rji(3)**2
          if( dji2.ge.rcmax2 ) cycle
          dji= dsqrt(dji2)
          js = int(tag(ja))
          do ihl2=1,mhl(2)
            hl2i= hl2(ihl2,ia)
            zl2i= zl2(ihl2,ia)
!!$            tmp2i= wgt23(ihl2) *hl2i*(1d0-hl2i)
            tmp2i= wgt23(ihl2) *dsigmoid(zl2i,hl2i)
            do ihl1=1,mhl(1)
              hl1i= hl1(ihl1,ia)
              zl1i= zl1(ihl1,ia)
!!$              tmp1i= wgt22(ihl1,ihl2) *hl1i*(1d0-hl1i)
              tmp1i= wgt22(ihl1,ihl2) *dsigmoid(zl1i,hl1i)
!!$              do ihl0=iaddr2(1,is,js),iaddr2(2,is,js)
              do ihl0=1,nhl(0)
                if( igsf(ihl0,jj,ia).eq.0 ) cycle
                do ixyz=1,3
                  do jxyz=1,3
!......derivative of gsf of atom-j by atom-i
                    sji= -tmp2i *tmp1i &
                         *wgt21(ihl0,ihl1) *dgsf(jxyz,ihl0,jj,ia) &
                         *rji(ixyz)
!!$!.....derivative of gsf of atom-i by atom-i
!!$                    sii= tmp2i *tmp1i &
!!$                         *wgt21(ihl0,ihl1) *dgsf(jxyz,ihl0,jj,ia) &
!!$                         *rij(ixyz)
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
  subroutine set_paramsdir_NN2(dname)
!
!  Accessor routine for setting paramsdir
!
    implicit none
    character(len=*),intent(in):: dname

    paramsdir = trim(dname)
    return
  end subroutine set_paramsdir_NN2
!=======================================================================
  subroutine gradw_NN2(namax,natm,tag,ra,nnmax &
       ,h,rc,lspr,iprint,ndimp,gwe,gwf,gws &
       ,lematch,lfmatch,lsmatch,iprm0)
!=======================================================================
!  Gradient w.r.t. NN weights, {w}
!  Note: This routine is always called in single run,
!  thus no need of parallel implementation.
!  Currently only 1 hidden layer is implemented.
!=======================================================================
    use descriptor,only: gsf,dgsf,igsf,nsf,nnl,nal,mskgfs
    use util,only: itotOf
    implicit none
    integer,intent(in):: namax,natm,nnmax,iprint,iprm0
    integer,intent(in):: lspr(0:nnmax,namax)
    real(8),intent(in):: ra(3,namax),h(3,3),rc,tag(namax)
    integer,intent(in):: ndimp
    real(8),intent(inout):: gwe(ndimp),gwf(ndimp,3,natm),gws(ndimp,6)
    logical,intent(in):: lematch,lfmatch,lsmatch

    integer:: iv,ia,ihl0,ihl1,jj,ja,jra
    real(8):: g,h1,z1,tmp,w2,w1,ds,dds,ftmp(3),xi(3),xj(3),xij(3),rij(3)
!!$    integer,external:: itotOf
    real(8),allocatable:: dgsf2(:,:,:,:)

!.....TO CHECK: Need to make hl1 every time?
    if( size(hl1).ne.nhl(1)*nal ) then
      deallocate(hl1,zl1)
      allocate(hl1(nhl(1),nal), zl1(nhl(1),nal))
    endif
    do ia=1,natm
      do ihl1=1,mhl(1)
        tmp= 0d0
        do ihl0=1,nhl(0)
          tmp= tmp +wgt11(ihl0,ihl1) *gsf(ihl0,ia)
        enddo
        zl1(ihl1,ia)= tmp
        hl1(ihl1,ia)= sigmoid(tmp)
      enddo
    enddo

    if( lematch ) then
      do ia=1,natm
        iv = iprm0
        do ihl0=1,nhl(0)
          g = gsf(ihl0,ia)
          do ihl1=1,mhl(1)
            w2 = wgt12(ihl1)
            h1 = hl1(ihl1,ia)
            z1 = zl1(ihl1,ia)
            iv = iv + 1
            gwe(iv) = gwe(iv) +w2 *g *dsigmoid(z1,h1)
          enddo
        enddo
        do ihl1=1,nhl(1)
          h1 = hl1(ihl1,ia)
          iv = iv + 1
          gwe(iv) = gwe(iv) + (h1 -0.5d0)
        enddo
      enddo
    endif

    if( lfmatch .or. lsmatch ) then
!.....Make dgsf2 array
      if( .not.allocated(dgsf2) ) then
        allocate(dgsf2(3,0:nnl,mhl(1),nal))
      else if( size(dgsf2).ne.3*nnl*mhl(1)*nal ) then
        deallocate(dgsf2)
        allocate(dgsf2(3,0:nnl,mhl(1),nal))
      endif
      dgsf2(:,:,:,:) = 0d0
      do ia=1,natm
        do ihl0=1,nhl(0)
          do ihl1=1,mhl(1)
            w1 = wgt11(ihl0,ihl1)
            do jj=0,lspr(0,ia)  ! Notice: from 0 (ja==ia) not 1
              dgsf2(1:3,jj,ihl1,ia) = dgsf2(1:3,jj,ihl1,ia) &
                   +w1 *dgsf(1:3,ihl0,jj,ia)
            enddo
          enddo
        enddo
      enddo
    endif

    if( lfmatch .and. .not.lsmatch ) then
!.....Copy back derivatives of forces on atoms outside of the node
!.....Compute derivative of forces w.r.t. weights
      do ia=1,natm
        iv = iprm0
!.....Weights between layer 0 and 1
        do ihl0=1,nhl(0)
          g = gsf(ihl0,ia)
          do ihl1=1,mhl(1)
            w1 = wgt11(ihl0,ihl1)
            w2 = wgt12(ihl1)
            h1 = hl1(ihl1,ia)
            z1 = zl1(ihl1,ia)
            ds = dsigmoid(z1,h1)
            dds= ddsigmoid(z1,h1)
            iv = iv +1
            do jj=0,lspr(0,ia)  ! Notice: from 0 not 1
              if( jj.eq.0 ) then
                ja = ia
              else
                ja = lspr(jj,ia)
              endif
              jra = itotOf(tag(ja))
              ftmp(1:3) = -w2 *(dds *gsf(ihl0,ia) *dgsf2(1:3,jj,ihl1,ia) &
                   +ds*dgsf(1:3,ihl0,jj,ia) )
!.....Derivative of forces wrt weights
              gwf(iv,1:3,jra) = gwf(iv,1:3,jra) +ftmp(1:3)
            enddo  ! jj=
          enddo  ! ihl1=
!!$          endif
        enddo  ! ihl0=
!.....Weights between layer-1 and output
        do ihl1=1,mhl(1)
          w2 = wgt12(ihl1)
          h1 = hl1(ihl1,ia)
          z1 = zl1(ihl1,ia)
          tmp = dsigmoid(z1,h1)
          iv = iv +1
          do jj=0,lspr(0,ia)
            if( jj.eq.0 ) then
              ja = ia
            else
              ja = lspr(jj,ia)
            endif
            jra = itotOf(tag(ja))
            do ihl0=1,nhl(0)
              w1 = wgt11(ihl0,ihl1)
              ftmp(1:3) = -w1 *tmp *dgsf(1:3,ihl0,jj,ia)
!.....Derivative of force
              gwf(iv,1:3,jra) = gwf(iv,1:3,jra) +ftmp(1:3)
            enddo  ! ihl0=
          enddo  ! jj=
        enddo  ! ihl1=
      enddo  ! ia=
    endif

! Since stress-matching causes considerable increase of computational cost,
! only compute stress contribution when stress-matching is ON.
! But even if force-matching is OFF, force contribution is computed,
! because the addictional cost for force contribution is not very much.
    if( lsmatch ) then
!.....Compute derivative of forces w.r.t. weights
      do ia=1,natm
        iv = iprm0
        xi(1:3) = ra(1:3,ia)
!.....Weights between layer 0 and 1
        do ihl0=1,nhl(0)
          g = gsf(ihl0,ia)
          do ihl1=1,mhl(1)
            w1 = wgt11(ihl0,ihl1)
            w2 = wgt12(ihl1)
            h1 = hl1(ihl1,ia)
            z1 = zl1(ihl1,ia)
            ds = dsigmoid(z1,h1)
            dds= ddsigmoid(z1,h1)
            iv = iv +1
            do jj=0,lspr(0,ia)  ! Notice: from 0 not 1
              if( jj.eq.0 ) then
                ja = ia
              else
                ja = lspr(jj,ia)
                xj(1:3) = ra(1:3,ja)
                xij(1:3) = xj(1:3) -xi(1:3)
                rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
              endif
              jra = itotOf(tag(ja))
              ftmp(1:3) = -w2 *(dds *gsf(ihl0,ia) *dgsf2(1:3,jj,ihl1,ia) &
                   +ds*dgsf(1:3,ihl0,jj,ia) )
!.....Derivative of forces wrt weights
              gwf(iv,1:3,jra) = gwf(iv,1:3,jra) +ftmp(1:3)
!.....Derivative of stress wrt weights
              if( jj.eq.0 ) cycle
              gws(iv,1) = gws(iv,1) +rij(1)*ftmp(1)
              gws(iv,2) = gws(iv,2) +rij(2)*ftmp(2)
              gws(iv,3) = gws(iv,3) +rij(3)*ftmp(3)
              gws(iv,4) = gws(iv,4) +rij(2)*ftmp(3)
              gws(iv,5) = gws(iv,5) +rij(1)*ftmp(3)
              gws(iv,6) = gws(iv,6) +rij(1)*ftmp(2)
            enddo  ! jj=
          enddo  ! ihl1=
!!$          endif
        enddo  ! ihl0=
!.....Weights between layer-1 and output
        do ihl1=1,mhl(1)
          w2 = wgt12(ihl1)
          h1 = hl1(ihl1,ia)
          z1 = zl1(ihl1,ia)
          tmp = dsigmoid(z1,h1)
          iv = iv +1
          do jj=0,lspr(0,ia)
            if( jj.eq.0 ) then
              ja = ia
            else
              ja = lspr(jj,ia)
              xj(1:3) = ra(1:3,ja)
              xij(1:3) = xj(1:3) -xi(1:3)
              rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
            endif
            jra = itotOf(tag(ja))
            if( jj.eq.0 ) then
              do ihl0=1,nhl(0)
                w1 = wgt11(ihl0,ihl1)
                ftmp(1:3) = -w1 *tmp *dgsf(1:3,ihl0,jj,ia)
!.....Derivative of force
                gwf(iv,1:3,jra) = gwf(iv,1:3,jra) +ftmp(1:3)
              enddo  ! ihl0=
            else
              do ihl0=1,nhl(0)
                w1 = wgt11(ihl0,ihl1)
                ftmp(1:3) = -w1 *tmp *dgsf(1:3,ihl0,jj,ia)
!.....Derivative of force
                gwf(iv,1:3,jra) = gwf(iv,1:3,jra) +ftmp(1:3)
!.....Derivative of stress
                gws(iv,1) = gws(iv,1) +rij(1)*ftmp(1)
                gws(iv,2) = gws(iv,2) +rij(2)*ftmp(2)
                gws(iv,3) = gws(iv,3) +rij(3)*ftmp(3)
                gws(iv,4) = gws(iv,4) +rij(2)*ftmp(3)
                gws(iv,5) = gws(iv,5) +rij(1)*ftmp(3)
                gws(iv,6) = gws(iv,6) +rij(1)*ftmp(2)
              enddo  ! ihl0=
            endif
          enddo  ! jj=
        enddo  ! ihl1=
      enddo  ! ia=

    endif

    return
  end subroutine gradw_NN2
!=======================================================================
  subroutine get_NN2_hl1(hl1o)
!
!  Access to hl1 from outside
!
    use descriptor,only: nal
    real(8),intent(out):: hl1o(nhl(1),nal)

    hl1o(1:nhl(1),1:nal) = hl1(1:nhl(1),1:nal)
    return
  end subroutine get_NN2_hl1
!=======================================================================
  subroutine set_NN2_hl1(hl1o)
!
!  Set hl1 from outside
!
    use descriptor,only: nal
    real(8),intent(in):: hl1o(nhl(1),nal)

    integer:: ihl1

    if( size(hl1o).ne.size(hl1) ) then
      print *,'ERROR: size of hl1o and hl1 different !'
      print *,'  size(hl1o),size(hl1)=',size(hl1o),size(hl1)
      stop
    endif
    hl1(:,:) = hl1o(:,:)
!!$    print *,'hl1 after set, nhl(1),mhl(1)=',nhl(1),mhl(1)
!!$    do ihl1=1,nhl(1)
!!$      print *,'ihl1,hl1(ihl1,1)=',ihl1,hl1(ihl1,1)
!!$    enddo
    return
  end subroutine set_NN2_hl1
!=======================================================================
  subroutine set_iglid_NN2(cpena,cfmethod)
!
!  Initialize some only required for fitpot.
!
    use descriptor,only: ngl,glval,iglid
    character(len=*),intent(in):: cpena,cfmethod

    integer:: i,ihl0,ihl1
    
!.....Make groups for group LASSO/FS
    if( trim(cpena).eq.'glasso' &
         .or. trim(cfmethod).eq.'gfs') then
      if( .not.allocated(iglid) ) allocate(iglid(nwtot))
      iglid(1:nwtot)= 0
      i= 0
      do ihl0=1,nhl(0)
        do ihl1=1,nhl(1)
          i= i +1
          iglid(i)= ihl0
        enddo
      enddo
!.....weights not connected to symmetry functions are not penalized in g-lasso
      do i=nhl(0)*nhl(1)+1,nwtot
        iglid(i)= -1
      enddo
    endif
    
  end subroutine set_iglid_NN2
!=======================================================================
  function mem_NN2()
    integer:: mem_NN2

    mem_NN2 = mem
    return
  end function mem_NN2
!=======================================================================
  function time_NN2()
    real(8):: time_NN2

    time_NN2 = time
    return
  end function time_NN2

end module NN2
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
