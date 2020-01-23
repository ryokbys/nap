module DNN
!-----------------------------------------------------------------------
!                     Last modified: <2020-01-22 18:11:09 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
!  Parallel implementation of deep neural-network potential.
!  See RK's memo 2020-01-21 for formulation details.
!  To separate the symmetry function calculations in descriptor.F90.
!-----------------------------------------------------------------------
#ifdef _DEBUG
#define DEBUG_PRINT(a) write(*,*) (a)
#else
#define DEBUG_PRINT(a)
#endif
  implicit none 
  save
  character(len=128):: paramsdir = '.'

!.....parameter file name
  character(128),parameter:: cpfname = 'in.params.DNN'

  real(8),parameter:: pi= 3.14159265358979d0

  integer:: mem
  real(8):: time
  
!.....logical flag for bias
  logical:: lbias = .true.

!.....Max num of species
  integer,parameter:: msp = 9
  logical:: interact(msp,msp)

!.....Parameters
  integer:: nlayer, maxnnode, nwtot
  integer,allocatable:: nhl(:), mhl(:),iactf(:)
  real(8),allocatable:: hls(:,:,:),gls(:,:,:),wgts(:,:,:)
  real(8),allocatable:: zls(:,:,:),sgm1(:,:,:),sgm2(:,:,:)

!.....parameters given from outside (fitpot)
  integer:: nprms
  real(8),allocatable:: prms(:)
  logical:: lprmset_DNN = .false.
  
!.....Sigmoid function types:
!     1) 1/(1+exp(-x))
!     2) 1/(1+exp(-x)) +asig*x
  integer:: itypesig = 2
!.....Coefficient of additional term in sigmoid2
  real(8):: asig = 0.01d0

contains
!=======================================================================
  subroutine force_DNN(namax,natm,tag,ra,nnmax,aa,strs,h,hi,tcom &
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
    integer:: i,j,k,l,m,n,is,js,ierr,ia,ja,ml0,ml1,ml2,jj
    real(8):: at(3),epotl,epott,hl1i,hl2i,tmp2,tmp1,tmp,zl1i,zl2i,time0
    real(8),allocatable,save:: strsl(:,:,:),aal(:,:),gw(:)

    integer:: itot
    character(len=8):: cnum

    call make_gsf_arrays(l1st,namax,natm &
         ,tag,nnmax,lspr,myid,mpi_world,iprint)

    if( l1st ) then

      if( allocated(hls) ) deallocate(hls,gls,zls,sgm1,sgm2,gw)
      allocate( hls(0:maxnnode,0:nlayer,nal), gls(maxnnode,nlayer,nal), &
           zls(maxnnode,nlayer,nal), sgm1(0:maxnnode,nlayer,nal), &
           sgm2(0:maxnnode,nlayer,nal), gw(0:nhl(0)))

      if( allocated(strsl) ) deallocate(strsl,aal)
      allocate(strsl(3,3,namax),aal(3,namax))

    endif ! l1st

    if( allocated(hls) .and. ubound(hls,3).lt.nal ) then
      deallocate(hls,gls,zls,sgm1,sgm2)
      allocate( hls(0:maxnnode,0:nlayer,nal), gls(maxnnode,nlayer,nal), &
           zls(maxnnode,nlayer,nal), sgm1(0:maxnnode,nlayer,nal), &
           sgm2(0:maxnnode,nlayer,nal))
    endif

    if( size(strsl).ne.3*3*namax ) then
      deallocate(strsl,aal)
      allocate(strsl(3,3,namax),aal(3,namax))
    endif

    call calc_desc(namax,natm,nb,nnmax,h &
         ,tag,ra,lspr,rcin,myid,mpi_world,l1st,iprint)

    time0 = mpi_wtime()

    call comp_nodes(natm)

!.....Energies of atoms by summing up (N-1)-th node values
    epotl= 0d0
    do ia=1,natm
      epi(ia) = hls(1,nlayer,ia)
      epotl = epotl +epi(ia)
    enddo

!.....Forces of atoms
    strsl(1:3,1:3,1:namax) = 0d0
    aal(1:3,1:namax) = 0d0
    do ia=1,natm
      is = int(tag(ia))
      do ml0=0,nhl(0)
        gw(ml0) = 0d0
        do ml1=1,nhl(1)
          gw(ml0)= gw(ml0) +gls(ml1,1,ia)*wgts(ml0,ml1,1)
        enddo
      enddo
!.....Derivative of SF of atom-i w.r.t. atom-j 
      do jj=1,lspr(0,ia)
        ja = lspr(jj,ia)
        js = int(tag(ja))
        do ml0=1,nhl(0)  ! no bias contribution to force
!!$          if( igsf(ml0,jj,ia).eq.0 ) cycle  ! either way, dgsf=0d0
          aal(1:3,ja) = aal(1:3,ja) -gw(ml0)*dgsf(1:3,ml0,jj,ia)
        enddo
      enddo
!.....Derivative of SF of atom-i w.r.t. atom-i
      do ml0=1,nhl(0)
!!$        if( igsf(ml0,0,ia).eq.0 ) cycle  ! either way, dgsf=0d0
        aal(1:3,ia) = aal(1:3,ia) -gw(ml0)*dgsf(1:3,ml0,0,ia)
      enddo
    enddo

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
    if( iprint.gt.2 ) print *,'DNN epot = ',epott
    epot= epot +epott

    time = time +(mpi_wtime() -time0)

    return
  end subroutine force_DNN
!=======================================================================
  subroutine comp_nodes(natm)
!
!  Compute node values hls and zls by the forward propagation,
!  and gls by the back propgation.
!
    use descriptor,only: gsf,dgsf
    integer,intent(in):: natm

    integer:: ia,il,ml0,ml1,ml2,itype
    real(8):: z,sgmz,y

!.....Initialize
    hls(:,:,:) = 0d0
    zls(:,:,:) = 0d0
    hls(1:nhl(0),0,1:natm) = gsf(1:nhl(0),1:natm)
!!$    gsf(0,1:natm) = 1d0
!!$    dgsf(1:3,0,:,1:natm) = 0d0
!.....NOTE: 0-th node in layer is the bias.
    hls(0,0:nlayer,1:natm) = 1d0
    sgm1(0,0:nlayer,1:natm) = 0d0
    sgm2(0,0:nlayer,1:natm) = 0d0

!.....Compute the node values by forward propagation.
    do ia=1,natm
      do il=1,nlayer
        itype = iactf(il)
        do ml1=1,nhl(il)
          z = 0d0
          do ml0=0,nhl(il-1)
            z = z +wgts(ml0,ml1,il)*hls(ml0,il-1,ia)
          enddo
          zls(ml1,il,ia) = z
          sgmz = actf(itype,z)
          hls(ml1,il,ia) = sgmz
          sgm1(ml1,il,ia) = dactf(itype,z,sgmz)
          sgm2(ml1,il,ia) = ddactf(itype,z,sgmz)
        enddo ! ml1=...
      enddo ! il=...
!.....Compute gls by back propagation.
!.....NOTE: gls = prod(W*G) in RK's memo.
      gls(:,:,ia) = 0d0
      gls(1,nlayer,ia) = 1d0
      do il=nlayer-1,1,-1
        do ml1=1,nhl(il) ! Since sgm1(0,:,:)=0d0, no need to take ml1=0
          z = 0d0
          do ml2=1,nhl(il+1)
            z = z + gls(ml2,il+1,ia)*wgts(ml1,ml2,il+1)
          enddo
          gls(ml1,il,ia) = z*sgm1(ml1,il,ia)
        enddo ! ml0=...
      enddo ! il=...
      
    enddo
    
    return
  end subroutine comp_nodes
!=======================================================================
  subroutine set_sigtype_DNN(itype)
!
!  Set sigmoid function type
!
    integer,intent(in):: itype

    itypesig = itype
    return
  end subroutine set_sigtype_DNN
!=======================================================================
  function actf(itype,x)
!
!  Activation function, which could be DIFFERENT from sigmoid.
!
    integer,intent(in):: itype
    real(8),intent(in):: x
    real(8):: actf

    select case(itype)
    case(0)
      actf = x
    case(1)
      actf = 1d0/(1d0 +exp(-x))
    case(2)
      actf = 1d0/(1d0 +exp(-x)) +asig*x
    case default
      actf = 0d0
    end select
    return
  end function actf
!=======================================================================
  function dactf(itype,x,sx)
    integer,intent(in):: itype
    real(8),intent(in):: x,sx
    real(8):: dactf
    real(8):: sxt
    
    select case(itype)
    case(0)
      dactf = 1d0
    case(1)
      dactf = sx*(1d0-sx)
    case(2)
      sxt = 1d0/(1d0 +exp(-x))
      dactf = sxt *(1d0-sxt) +asig
    case default
      dactf = 0d0
    end select
    return
  end function dactf
!=======================================================================
  function ddactf(itype,x,sx)
    integer,intent(in):: itype
    real(8),intent(in):: x,sx
    real(8):: ddactf
    real(8):: sxt
    
    select case(itype)
    case(0)
      ddactf = 0d0
    case(1)
      ddactf = sx*(1d0-sx)*(1d0-2d0*sx)
    case(2)
      sxt = 1d0 /(1d0 +exp(-x))
      ddactf = sxt *(1d0 -sxt) *(1d0 -2d0*sxt)
    case default
      ddactf = 0d0
    end select
    return
  end function ddactf
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

    integer:: ia,ja,ixyz,jxyz,jj,is,js,ml0,ml1
    real(8):: xi(3),xj(3),xji(3),rij(3),rji(3),dji,dji2,sji
    real(8),allocatable,save:: gw(:)

    real(8),save:: rcmax2
    logical,save:: l1st = .true.

    if( l1st ) then
      rcmax2 = rc*rc
      l1st = .false.
      if( .not. allocated(gw) ) allocate(gw(0:nhl(0)))
    endif

    do ia=1,natm
      is = int(tag(ia))
      xi(1:3) = ra(1:3,ia)
      do ml0=0,nhl(0)
        gw(ml0) = 0d0
        do ml1=1,nhl(1)
          gw(ml0)= gw(ml0) +gls(ml1,1,ia)*wgts(ml0,ml1,1)
        enddo
      enddo
      do jj=1,lspr(0,ia)
        ja = lspr(jj,ia)
        xj(1:3)= ra(1:3,ja)
        xji(1:3)= xj(1:3)-xi(1:3)
        rji(1:3)= h(1:3,1)*xji(1) +h(1:3,2)*xji(2) +h(1:3,3)*xji(3)
        rij(1:3)= -rji(1:3)
        dji2= rji(1)**2 +rji(2)**2 +rji(3)**2
        if( dji2.ge.rcmax2 ) cycle
        dji = dsqrt(dji2)
        js = int(tag(ja))
        do ml0=1,nhl(0) ! no bias contribution to force/stress
          do ixyz=1,3
            do jxyz=1,3
              sji = -gw(ml0)*dgsf(jxyz,ml0,jj,ia)*rji(ixyz)
              strs(ixyz,jxyz,ja) = strs(ixyz,jxyz,ja) +sji
              strs(ixyz,jxyz,ia) = strs(ixyz,jxyz,ia) +sji
            enddo
          enddo
        enddo
      enddo
    enddo

!-----send back (3-body)forces, stresses, and potentials on immigrants
    call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,strs,9)
    strs(1:3,1:3,1:natm) = strs(1:3,1:3,1:natm)*0.5d0
    return
  end subroutine compute_stress
!=======================================================================
  subroutine gradw_DNN(namax,natm,tag,ra,nnmax &
       ,h,rc,lspr,iprint,ndimp,gwe,gwf,gws &
       ,lematch,lfmatch,lsmatch,iprm0)
!=======================================================================
!  Gradient w.r.t. NN weights, {w}
!  Note: This routine is always called in single run,
!  thus no need of parallel implementation.
!  See RK's memo 2020-01-21.
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
    integer:: ihl2
    real(8):: h2,z2,w3
!!$    integer,external:: itotOf
    real(8),allocatable:: dgsf2(:,:,:,:)

!.....TO CHECK: Need to make hl1 every time?
    if( size(hl1).ne.nhl(1)*nal ) then
      deallocate(hl1,zl1)
      allocate(hl1(nhl(1),nal), zl1(nhl(1),nal))
    endif
    if( nl.eq.2 .and. size(hl2).ne.nhl(2)*nal ) then
      deallocate(hl2,zl2)
      allocate(hl2(nhl(2),nal),zl2(nhl(2),nal))
    endif
    
    call comp_nodes(natm)

    if( lematch ) then
      if( nl.eq.1 ) then
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
      else if( nl.eq.2 ) then
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
            do ihl2=1,mhl(2)
              w3 = wgt23(ihl2)
              h2 = hl2(ihl2,ia)
              z2 = zl2(ihl2,ia)
              iv = iv + 1
              gwe(iv) = gwe(iv) +w3 *h1 *dsigmoid(z2,h2)
            enddo
          enddo
          do ihl2=1,nhl(2)
            h2 = hl2(ihl2,ia)
            iv = iv + 1
            gwe(iv) = gwe(iv) + (h2 -0.5d0)
          enddo
        enddo
      endif
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
      if( nl.eq.1 ) then
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
!!$      else if( nl.eq.2 ) then
        
      endif  ! nl=
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
  end subroutine gradw_DNN
!=======================================================================
  subroutine set_paramsdir_DNN(dname)
!
!  Accessor routine for setting paramsdir
!
    implicit none
    character(len=*),intent(in):: dname

    paramsdir = trim(dname)
    return
  end subroutine set_paramsdir_DNN
!=======================================================================
  subroutine read_params_DNN(myid,mpi_world,iprint)
!-----------------------------------------------------------------------
!  Assume that the descriptor information is already read.
!  Output of the NN is an energy of an atom-i, thus the number of nodes
!  in the output layer is 1.
!  Definition of "num of hidden layers", NL, is such that:
!    NL==1 if there are 1 input layer, 1 hidden layer and 1 output layer,
!    which means that NL does not count the input and outputlayer.
!  On the other hand, NLAYER = NL + 1, which counts an output layer.
!  Input file format is as follows:
!-----------------------------------------------------------------------
!  ! comments or options
!  NL, NHL(0:NL)
!  WGT(i),  LOWER(i),  UPPER(i)
!  ...
!-----------------------------------------------------------------------
!  - NL: number of layers
!  - NHL: number of nodes (0th: input, NL-th: output)
!  - WGT(i): weight of the i-th edge
!  - LOWER,UPPER: lower and upper limit of the weight
!-----------------------------------------------------------------------
    use descriptor,only: nsf
    use util, only: num_data
    implicit none
    include 'mpif.h'

    integer,intent(in):: myid,mpi_world,iprint

    integer:: ierr,i,j,k,nc,ndat,icmb(3),itmp,nl,il &
         ,ml,ml0,ml1,istart
    integer,allocatable:: nwgt(:)
    logical:: lexist
    character:: ctmp*128,fname*128

!.....Check if the file exists
    fname = trim(paramsdir)//'/'//trim(cpfname)
    inquire(file=trim(fname),exist=lexist)
    if( .not. lexist ) then
      if( myid.eq.0 ) then
        write(6,'(a)') ' [Error] '//trim(fname)//' does not exist !!!.'
        write(6,'(a)') '   The DNN potential needs '//trim(fname)//'.'
      endif
      call mpi_finalize(ierr)
      stop
    endif

!.....Read parameters at the 1st call
    if( myid.eq.0 ) then
      open(50,file=trim(fname),status='old')
!.....num of symmetry functions, num of nodes in hidden layers
10    read(50,'(a)') ctmp
      if( ctmp(1:1).eq.'!' .or. ctmp(1:1).eq.'#' ) then
        call parse_option(ctmp,iprint,ierr)
        goto 10
      else
        read(ctmp,*) nl
        backspace(50)
        nlayer = nl + 1
      endif
    endif

    call mpi_bcast(lbias,1,mpi_logical,0,mpi_world,ierr)
    call mpi_bcast(itypesig,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(nlayer,1,mpi_integer,0,mpi_world,ierr)
    allocate(nhl(0:nlayer),iactf(nlayer))
    if( myid.eq.0 ) read(50,*) itmp, (nhl(i),i=0,nlayer-1)
    nhl(nlayer) = 1
    iactf(1:nlayer-1) = itypesig
    iactf(nlayer) = 0
    call mpi_bcast(nhl,nlayer+1,mpi_integer,0,mpi_world,ierr)

    if( nhl(0).ne.nsf ) then
      print '(a,3(1x,i0))',' ERROR: nhl(0).ne.nsf!, myid,nhl(0),nsf=',myid,nhl(0),nsf
      if( myid.eq.0 ) then
        print *,'  Check the consistency in in.params.desc and in.params.DNN'
      endif
      call mpi_finalize(ierr)
      stop
    endif

!.....calc number of weights
    allocate(nwgt(nlayer))
    nwgt(1:nlayer) = 0
    maxnnode = nhl(0)
    do i=1,nlayer
!.....NOTE: including bias
      nwgt(i)= (nhl(i-1)+1)*(nhl(i))
      maxnnode = max(maxnnode,nhl(i))
    enddo
    nc= 0
    do i=1,nlayer
      nc= nc +nwgt(i)
    enddo
    nwtot = nc

    if( myid.le.0 .and. iprint.ne.0 ) then
      print *,''
      print *,'DNN potential parameters:'
      print '(a,i0,a)','   Num of layers = ',nlayer, '   (input & hidden, not incl. output)'
      print '(a,100(1x,i0))','   nhl(:) = ',nhl(0:nlayer)
      if( itypesig.eq.1 ) then
        print *,'  Activation function: 1) 1/(1+exp(-x))'
      else if( itypesig.eq.2 ) then
        print *,'  Activation function: 2) 1/(1+exp(-x)) +asig*x'
        print '(a,f7.4)','                          asig = ',asig
      else
        print *,'  Activation function: unknown'
      endif
      print '(a,i3,i5)',      '   ml, nhl(ml)           = ',0,nhl(0)
      do ml=1,nlayer
        print '(a,i3,i5,i6)', '   ml, nhl(ml), nwgt(ml) = ',ml,nhl(ml),nwgt(ml)
      enddo
      print '(a,i0)','   Max num of nodes in a layer = ',maxnnode
      print '(a,i0)','   Total num of weights = ',nwtot
      print *,''
    endif

!.....NOTE: W_{i,j,il} == wgts(j,i,il), row and column positions are swapped
!.....      and only the column(j) has bias(0-th) components.
    allocate(wgts(0:maxnnode,maxnnode,nlayer))
    if( myid.eq.0 ) then
      wgts(:,:,:) = 0d0
      istart = 0
      if( .not. lbias ) istart = 1
      do il=1,nlayer
        do ml0=istart,nhl(il-1)
          do ml1=1,nhl(il)
            read(50,*) wgts(ml0,ml1,il)
          enddo
        enddo
      enddo
    endif
    nc = (maxnnode+1)*maxnnode*nlayer
    call mpi_bcast(wgts,nc,mpi_real8,0,mpi_world,ierr)

    deallocate(nwgt)
    return
  end subroutine read_params_DNN
!=======================================================================
  subroutine set_params_DNN(nprms_in,prms_in,nl_in,nhl_in)
!
!  Accessor routine to set DNN parameters from outside.
!  Curretnly this routine is supposed to be called only on serial run.
!  NOTE that the definition of NL_IN and NLAYER could be different such as:
!    NL_IN does not count output layer, whereas NLAYER does.
!
    implicit none 
    integer,intent(in):: nprms_in,nl_in,nhl_in(0:nl_in)
    real(8),intent(in):: prms_in(nprms_in)

    integer:: i

    nlayer = nl_in +1
    if( nl_in.eq.0 ) then
      print *,'ERROR: NL_IN==0 which should not happen.'
      print *,'  Probably NN_num_layers and NN_num_nodes are not set in in.fitpot.'
      print *,'  Those should be consistent with the number of NN weigts to be optimized.'
      stop
    endif
    nhl(0:nlayer-1) = nhl_in(0:nl_in)
    nhl(nlayer) = 1

    nwtot = 0
    do i=1,nlayer
      nwtot = nwtot +(nhl(i-1)+1)*nhl(i)
    enddo

    if( nwtot.ne.nprms_in ) then
      print *,'ERROR: nl_in,nhl_in,nprms_in not consistent !!'
      print *,'  Check in.vars.fitpot or in.params.DNN and'&
           //' NN_num_nodes parameters in in.fitpot.'
      stop
    endif
    
    nprms = nprms_in
    if( .not.allocated(prms) ) allocate(prms(nprms))
    prms(1:nprms) = prms_in(1:nprms_in)

    lprmset_DNN = .true.
    return
  end subroutine set_params_DNN
!=======================================================================
  subroutine check_consistency(myid,iprint)
    implicit none
    integer,intent(in):: myid,iprint

    integer:: i,nc
    integer,allocatable:: nwgt(:)

!.....calc number of weights
    allocate(nwgt(nlayer))
    nwgt(1:nlayer) = 0
    do i=1,nlayer
      nwgt(i)= (nhl(i-1)+1)*nhl(i)
    enddo
    
    nc= 0
    do i=1,nlayer
      nc= nc +nwgt(i)
    enddo
    if( nprms .ne. nc ) then
      write(6,'(a)') ' ERROR: Number of parameters is not correct !!!'
      write(6,'(a,i0)')  '   nprms = ',nprms
      write(6,'(a,i0)')  '   where nprms should be ',nc
      stop
    endif

    deallocate(nwgt)
    return
  end subroutine check_consistency
!=======================================================================
  subroutine update_params_DNN()
!
!  Update DNN parameters by taking parameter values from params array.
!  This routine would be called only from externally from fitpot.
!
    integer:: i,inc,il,ml0,ml1

    if( .not.lprmset_DNN ) then
      print *,'ERROR: params have not been set yet.'
      stop
    endif

    maxnnode = nhl(0)
    do i=1,nlayer
      maxnnode = max(maxnnode,nhl(i))
    enddo
    if( .not.allocated(wgts) ) allocate(wgts(0:maxnnode,maxnnode,nlayer))

    inc = 0
    do il=1,nlayer
      do ml1=1,nhl(il)
        do ml0=0,nhl(il-1)
          inc = inc + 1
          wgts(ml0,ml1,il) = prms(inc)
        enddo
      enddo
    enddo
    
    return
  end subroutine update_params_DNN
!=======================================================================
  subroutine parse_option(cline,iprint,ierr)
!
!  Parse options from a comment line.
!  Lines starting from ! or # are treated as comment lines,
!  but options can be given at the comment lines.
!  The option words should be put after these comment characters with
!  one or more spaces between them for example,
!
!  sigtype:  .true.
!
!  Currently available options are:
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
    if( index(cline,'sigtype:').ne.0 ) then
      read(cline,*) c1,copt,iopt
      if( trim(copt).ne.'sigtype:' ) then
        print *, 'Error: copt is not "sigtype:" !!!'
        ierr = 2
      endif
      itypesig = iopt
    else if( index(cline,'bias:').ne.0 ) then
      read(cline,*) c1,copt,lopt
      if( trim(copt).ne.'bias:' ) then
        print *, 'Error: copt is not "bias:" !!!'
        ierr = 1
      endif
      lbias = lopt
    endif
    
  end subroutine parse_option
!=======================================================================
  function mem_DNN()
    integer:: mem_DNN

    mem_DNN = mem
    return
  end function mem_DNN
!=======================================================================
  function time_DNN()
    real(8):: time_DNN

    time_DNN = time
    return
  end function time_DNN
end module DNN
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
