module DNN
!-----------------------------------------------------------------------
!                     Last modified: <2020-01-26 12:37:03 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
!  Parallel implementation of deep neural-network potential.
!  See RK's memo 2020-01-21 for formulation details.
!  To separate the symmetry function calculations in descriptor.F90.
!-----------------------------------------------------------------------
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
  integer,allocatable:: nhl(:), mhl(:),iactf(:),nwgt(:)
  real(8),allocatable:: hls(:,:),gls(:,:),wgts(:,:,:)
  real(8),allocatable:: zls(:,:),sgm1(:,:),sgm2(:,:)

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

!.....Whether or not called from fitpot [default: .false.]
  logical:: lfitpot = .false.

contains
!=======================================================================
  subroutine force_DNN(namax,natm,tag,ra,nnmax,aa,strs,h,hi,tcom &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rcin,lspr &
       ,mpi_world,myid,epi,epot,nismax,lstrs,iprint,l1st)
    use descriptor,only: gsf,dgsf,igsf,nsf,nal,calc_desc,make_gsf_arrays, &
         lupdate_gsf,gsfi,dgsfi,igsfi,calc_desci,prepare_desci,nnl
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
    integer:: i,j,k,l,m,n,is,js,ierr,ia,ja,ml,ml0,ml1,ml2,jj,ixyz,jxyz
    real(8):: at(3),epotl,epott,hl1i,hl2i,tmp2,tmp1,tmp,zl1i,zl2i,time0
    real(8):: xi(3),xj(3),xji(3),rji(3),rij(3),dji,dji2,sji
    real(8),allocatable,save:: strsl(:,:,:),aal(:,:),gw(:)
    real(8),save:: rcmax2
    integer:: itot
    character(len=8):: cnum

    call make_gsf_arrays(l1st,namax,natm &
         ,tag,nnmax,lspr,myid,mpi_world,iprint)

    if( l1st ) then

      if( allocated(hls) ) deallocate(hls,gls,zls,sgm1,sgm2,gw)
      allocate( hls(0:maxnnode,0:nlayer), gls(maxnnode,nlayer+1), &
           zls(maxnnode,nlayer), sgm1(0:maxnnode,nlayer), &
           sgm2(0:maxnnode,nlayer), gw(0:maxnnode))

      if( allocated(strsl) ) deallocate(strsl,aal)
      allocate(strsl(3,3,namax),aal(3,namax))

!.....Set activation function type here
      iactf(1:nlayer-1) = itypesig
      iactf(nlayer) = 0
      
!.....Compute memory used
      mem = mem +8*size(hls) +8*size(gls) +8*size(zls) +8*size(sgm1) &
           +8*size(sgm2) +8*size(gw) +8*size(strsl) +8*size(aal)

      rcmax2 = rcin*rcin

      if( myid.le.0 .and. iprint.ne.0 ) then
        print *,''
        print *,'DNN potential parameters:'
        print '(a,i0,a)','   Num of layers = ',nlayer, '   (input & hidden, not incl. output)'
        print '(a,100(1x,i0))','   Num of nodes in each layer = ',nhl(0:nlayer)
        if( itypesig.eq.1 ) then
          print *,'   Activation function: 1) 1/(1+exp(-x))'
        else if( itypesig.eq.2 ) then
          print '(a,f7.4)','   Activation function: 2) 1/(1+exp(-x))+asig*x' &
               //', w/ asig=',asig
        else
          print *,'   Activation function: unknown'
        endif
        print '(a,i3,i5)',      '   ml, nhl(ml)           = ',0,nhl(0)
        do ml=1,nlayer
          print '(a,i3,i5,i6)', '   ml, nhl(ml), nwgt(ml) = ',ml,nhl(ml),nwgt(ml)
        enddo
        print '(a,i0)','   Max num of nodes in a layer = ',maxnnode
        print '(a,i0)','   Total num of weights = ',nwtot
        print '(a,f10.3,a)','   Memory in force_DNN = ', &
             dble(mem)/1000/1000, ' MB'
      endif

      call prepare_desci(myid,iprint,rcin)

    endif ! l1st

    if( allocated(hls) .and. size(hls).eq.(maxnnode+1)*(nlayer+1) ) then
      deallocate(hls,gls,zls,sgm1,sgm2,gw)
      allocate( hls(0:maxnnode,0:nlayer), gls(maxnnode,nlayer+1), &
           zls(maxnnode,nlayer), sgm1(0:maxnnode,nlayer), &
           sgm2(0:maxnnode,nlayer), gw(0:maxnnode))
    endif

    if( size(strsl).ne.3*3*namax ) then
      deallocate(strsl,aal)
      allocate(strsl(3,3,namax),aal(3,namax))
    endif

    if( iprint.gt.1 ) then
      print *,'nal, nnl    =',nal,nnl
      print *,'shape(gsf)  =',shape(gsf)
      print *,'shape(dgsf) =',shape(dgsf)
      print *,'shape(gsfi) =',shape(gsfi)
      print *,'shape(dgsfi)=',shape(dgsfi)
    end if

!!$    call calc_desc(namax,natm,nb,nnmax,h &
!!$         ,tag,ra,lspr,rcin,myid,mpi_world,l1st,iprint)

    time0 = mpi_wtime()

!.....Eenergies, forces and stresses of atoms
    epotl= 0d0
    strsl(1:3,1:3,1:namax) = 0d0
    aal(1:3,1:namax) = 0d0
    do ia=1,natm
      call calc_desci(ia,namax,natm,nnmax,h &
           ,tag,ra,lspr,rcin,iprint)
      is = int(tag(ia))
      xi(1:3) = ra(1:3,ia)
      call comp_nodes_of(ia)
      epi(ia) = hls(1,nlayer)
      epotl = epotl +epi(ia)
!.....Pre-compute series of gl*Wl
      do ml0=0,nhl(0)
        gw(ml0) = 0d0
        do ml1=1,nhl(1)
          gw(ml0)= gw(ml0) +gls(ml1,1)*wgts(ml0,ml1,1)
        enddo
      enddo
!.....Derivative of SF of atom-i w.r.t. atom-j 
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
        do ml0=1,nhl(0)  ! no bias contribution to force
!!$          aal(1:3,ja) = aal(1:3,ja) -gw(ml0)*dgsf(1:3,ml0,jj,ia)
          aal(1:3,ja) = aal(1:3,ja) -gw(ml0)*dgsfi(1:3,ml0,jj)
!.....Stress
          do ixyz=1,3
            do jxyz=1,3
!!$              sji = -gw(ml0)*dgsf(jxyz,ml0,jj,ia)*rji(ixyz)
              sji = -gw(ml0)*dgsfi(jxyz,ml0,jj)*rji(ixyz)
              strs(ixyz,jxyz,ja) = strs(ixyz,jxyz,ja) +sji
              strs(ixyz,jxyz,ia) = strs(ixyz,jxyz,ia) +sji
            enddo
          enddo
        enddo ! ml0=
      enddo
!.....Derivative of SF of atom-i w.r.t. atom-i
      do ml0=1,nhl(0)
!!$        aal(1:3,ia) = aal(1:3,ia) -gw(ml0)*dgsf(1:3,ml0,0,ia)
        aal(1:3,ia) = aal(1:3,ia) -gw(ml0)*dgsfi(1:3,ml0,0)
      enddo
    enddo

!.....Send back forces on immigrant atoms
    call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,aal,3)
    aa(1:3,1:natm) = aa(1:3,1:natm) +aal(1:3,1:natm)

!.....Send back stresses on immigrant atoms
    call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,strsl,9)
    strs(1:3,1:3,1:natm) = strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)*0.5d0

!.....Gather epot
    call mpi_allreduce(epotl,epott,1,mpi_double_precision &
         ,mpi_sum,mpi_world,ierr)
    if( iprint.gt.2 ) print *,'DNN epot = ',epott
    epot= epot +epott

    time = time +(mpi_wtime() -time0)

    return
  end subroutine force_DNN
!=======================================================================
  subroutine gradw_DNN(namax,natm,tag,ra,nnmax &
       ,h,rc,lspr,iprint,ndimp,gwe,gwf,gws &
       ,lematch,lfmatch,lsmatch,iprm0)
!=======================================================================
!  Gradient w.r.t. NN weights, {w}
!  Note: This routine is always called in single run,
!  thus no need of parallel implementation.
!=======================================================================
    use descriptor,only: gsf,dgsf,igsf,nsf,nnl,nal,mskgfs, &
         gsfi,dgsfi,igsfi,calc_desci
    use util,only: itotOf
    implicit none
    integer,intent(in):: namax,natm,nnmax,iprint,iprm0
    integer,intent(in):: lspr(0:nnmax,namax)
    real(8),intent(in):: ra(3,namax),h(3,3),rc,tag(namax)
    integer,intent(in):: ndimp
    real(8),intent(inout):: gwe(ndimp),gwf(ndimp,3,natm),gws(ndimp,6)
    logical,intent(in):: lematch,lfmatch,lsmatch

    integer:: iv,ia,jj,il,ml1,ml0,ml2,ml,nni,n,mn0,mn1,l,memg
    real(8):: tmp
    real(8),allocatable,save:: fls(:,:,:,:),wfgw(:,:,:,:),wsgm1(:,:),gw(:)

    if( .not. allocated(fls) ) then
      allocate(fls(3,0:nnmax,0:maxnnode,0:nlayer), gw(0:maxnnode), &
           wfgw(3,0:nnmax,maxnnode,nlayer), wsgm1(maxnnode,maxnnode))

      if( iprint.ne.0 ) then
        memg = 8*size(fls) +8*size(gw) +8*size(wfgw) +8*size(wsgm1)
        print '(a,f10.3,a)',' Memnory in gradw_DNN = ',dble(memg)/1000/1000,' MB'
      endif
    endif

    do ia=1,natm
      call calc_desci(ia,namax,natm,nnmax,h,tag,ra,lspr,rc,iprint)
      call comp_nodes_of(ia)
      if( lematch ) then
        iv = 0
        do il=1,nlayer
          do ml1=1,nhl(il)
            do ml0=0,nhl(il-1)
              iv = iv + 1
              gwe(iv) = gwe(iv) +gls(ml1,il)*hls(ml0,il-1)
            enddo
          enddo
        enddo
      endif  ! lematch

!!$      if( lfmatch .or. lsmatch ) then
!!$!.....1st, create fls(:,:,:,:) by forward propagation that are
!!$!.....required in force and stress matching
!!$        nni = lspr(0,ia)
!!$!.....NOTE: fls of bias is 0; fls(:,:,0,:) = 0d0
!!$        fls(:,:,:,:) = 0d0
!!$        do jj=0,nni
!!$          do ml0=1,nhl(0)
!!$            fls(1:3,jj,ml0,0) = dgsfi(1:3,ml0,jj)
!!$          enddo
!!$        enddo
!!$        do il=1,nlayer
!!$          do ml1=1,nhl(il)
!!$            do ml0=1,nhl(il-1)
!!$              fls(1:3,0:nni,ml1,il) = fls(1:3,0:nni,ml1,il) &
!!$                   + wgts(ml0,ml1,il)*fls(1:3,0:nni,ml0,il-1)
!!$            enddo
!!$            fls(1:3,0:nni,ml1,il) = fls(1:3,0:nni,ml1,il)*sgm1(ml1,il)
!!$          enddo
!!$        enddo ! il=...
!!$        wfgw(:,:,:,:) = 0d0
!!$        do il=1,nlayer
!!$          do ml1=0,nhl(il)
!!$            gw(ml1) = 0d0
!!$            do ml2=1,nhl(il+1)
!!$              gw(ml1)= gw(ml1) +gls(ml2,il+1)*wgts(ml1,ml2,il+1)
!!$            enddo
!!$          enddo
!!$          do ml1=1,nhl(il)
!!$            wfgw(1:3,0:nni,ml1,il) = 0d0
!!$            do ml0=0,nhl(il-1)
!!$              wfgw(1:3,0:nni,ml1,il) = wfgw(1:3,0:nni,ml1,il) &
!!$                   +wgts(ml0,ml1,il)*fls(1:3,0:nni,ml0,il-1)
!!$            enddo
!!$            wfgw(1:3,0:nni,ml1,il) = wfgw(1:3,0:nni,ml1,il)*gw(ml1)
!!$          enddo
!!$        enddo
!!$      endif ! lfmatch .or. lsmatch
!!$      
!!$      if( lfmatch ) then
!!$!.....Direct derivative of force term w.r.t. W_l
!!$        iv = 0
!!$        do il=1,nlayer
!!$          do ml1=1,nhl(il)
!!$            do ml0=0,nhl(il-1)
!!$              iv = iv + 1
!!$              do jj=0,nni
!!$                gwf(iv,1:3,ia) = gwf(iv,1:3,ia) &
!!$                     +gls(ml1,il)*fls(1:3,jj,ml0,il-1)
!!$              enddo
!!$            enddo
!!$          enddo
!!$        enddo
!!$!.....Indirect derivative of force term w.r.t. W_l
!!$        iv = 0
!!$        do n=1,nlayer
!!$          do mn1=1,nhl(n)
!!$            do mn0=0,nhl(n-1)
!!$              iv = iv +1
!!$              do l=n,nlayer
!!$                wsgm1 = wxs(l,n)
!!$                tmp = 0d0
!!$                do ml=1,nhl(l)
!!$                  tmp = wsgm1(mn1,ml)*sgm2(ml,l)*hls(mn0,n)
!!$                  do jj=0,nni
!!$                    gwf(iv,1:3,ia) = gwf(iv,1:3,ia) &
!!$                         +tmp*wfgw(1:3,jj,ml,l)
!!$                  enddo
!!$                enddo
!!$!                if( ia.eq.1 .and. n.eq.nlayer .and. mn0.eq.0 ) then
!!$!                  print '(a,6i5)','ia,n,mn1,mn0,iv,l=',ia,n,mn1,mn0,iv,l
!!$!                  print '(a,i5,6es11.3)','nhl(l),tmp,wsgm1(1,1),sgm2(1,2)=' &
!!$!                       ,nhl(l),tmp,wsgm1(1,1),sgm2(1,2),gwf(iv,1:3,ia)
!!$!                endif
!!$              enddo ! l=...
!!$            enddo ! ml0=...
!!$          enddo ! ml1=...
!!$        enddo ! n=...
!!$      endif ! lfmatch
    enddo  ! ia=...

    return
  end subroutine gradw_DNN
!=======================================================================
  subroutine comp_nodes_of(ia)
!
!  Compute node values hls and zls of atom-IA by the forward propagation,
!  and gls by the back propgation.
!
    use descriptor,only: gsf,dgsf,gsfi,dgsfi
    integer,intent(in):: ia

    integer:: il,ml0,ml1,ml2,itype
    real(8):: z,sgmz,y

!.....Initialize
    hls(:,:) = 0d0
    zls(:,:) = 0d0
!!$    hls(1:nhl(0),0) = gsf(1:nhl(0),ia)
    hls(1:nhl(0),0) = gsfi(1:nhl(0))
!.....NOTE: 0-th node in layer is the bias.
    hls(0,0:nlayer) = 1d0
    sgm1(0,nlayer) = 0d0
    sgm2(0,nlayer) = 0d0

!.....Compute the node values by forward propagation.
    do il=1,nlayer
      itype = iactf(il)
      do ml1=1,nhl(il)
        z = 0d0
        do ml0=0,nhl(il-1)
          z = z +wgts(ml0,ml1,il)*hls(ml0,il-1)
        enddo
        zls(ml1,il) = z
        sgmz = actf(itype,z)
        hls(ml1,il) = sgmz
        sgm1(ml1,il) = dactf(itype,z,sgmz)
        sgm2(ml1,il) = ddactf(itype,z,sgmz)
      enddo ! ml1=...
    enddo ! il=...
!.....Compute gls by back propagation.
!.....NOTE: gls = prod(W*G) in RK's memo.
    gls(:,:) = 0d0
    gls(:,nlayer:nlayer+1) = 1d0
    do il=nlayer-1,1,-1
      do ml1=1,nhl(il) ! Since sgm1(0,:,:)=0d0, no need to take ml1=0
        z = 0d0
        do ml2=1,nhl(il+1)
          z = z + gls(ml2,il+1)*wgts(ml1,ml2,il+1)
        enddo
        gls(ml1,il) = z*sgm1(ml1,il)
      enddo ! ml0=...
    enddo ! il=...


    return
  end subroutine comp_nodes_of
!=======================================================================
  function wxs(l,n)
!
!  Compute Prod_{k=1}^{l-n} W_{l-k+1}G_{l-k} ==> wxs
!
    integer,intent(in):: l,n
    real(8):: wxs(maxnnode,maxnnode)

    integer:: i,k,ml,ml0,ml1
    real(8):: ymat(0:maxnnode,maxnnode),tmp(maxnnode,maxnnode)

    wxs(:,:) = 0d0
    do i=1,maxnnode
      wxs(i,i) = 1d0
    enddo
    do k=1,l-n
      do ml1=1,nhl(l-k+1)
        do ml0=0,nhl(l-k)
          ymat(ml0,ml1) = wgts(ml0,ml1,l-k+1)*sgm1(ml0,l-k)
        end do
      enddo
      do ml1=1,nhl(l)
        do ml0=1,nhl(l-k)
          tmp(ml0,ml1) = 0d0
          do ml=1,nhl(l-k+1)
            tmp(ml0,ml1) = tmp(ml0,ml1) +wxs(ml,ml1)*ymat(ml0,ml)
          end do
        enddo
      enddo
      wxs(1:nhl(l-k),1:nhl(l)) = tmp(1:nhl(l-k),1:nhl(l))
    enddo
    return
  end function wxs
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
    if( .not. allocated(nhl)) allocate(nhl(0:nlayer+1),iactf(nlayer)&
         ,nwgt(nlayer))
    if( myid.eq.0 ) read(50,*) itmp, (nhl(i),i=0,nlayer-1)
    nhl(nlayer) = 1
    nhl(nlayer+1) = 1
    call mpi_bcast(nhl,nlayer+2,mpi_integer,0,mpi_world,ierr)

    if( nhl(0).ne.nsf ) then
      print '(a,3(1x,i0))',' ERROR: nhl(0).ne.nsf!, myid,nhl(0),nsf=',myid,nhl(0),nsf
      if( myid.eq.0 ) then
        print *,'  Check the consistency in in.params.desc and in.params.DNN'
      endif
      call mpi_finalize(ierr)
      stop
    endif

!.....calc number of weights
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

!.....NOTE: W_{i,j,il} == wgts(j,i,il), row and column positions are swapped
!.....      and only the column(j) has bias(0-th) components.
    allocate(wgts(0:maxnnode,maxnnode,nlayer+1))
    if( myid.eq.0 ) then
      wgts(:,:,:) = 0d0
      wgts(:,:,nlayer+1) = 1d0
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
    nc = (maxnnode+1)*maxnnode*(nlayer+1)
    call mpi_bcast(wgts,nc,mpi_real8,0,mpi_world,ierr)

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
    if( .not. allocated(nhl) ) allocate(nhl(0:nlayer+1),iactf(nlayer)&
         ,nwgt(nlayer))
    nhl(0:nlayer-1) = nhl_in(0:nl_in)
    nhl(nlayer) = 1
    nhl(nlayer+1) = 1

    nwtot = 0
    maxnnode = nhl(0)
    nwgt(1:nlayer) = 0
    do i=1,nlayer
      nwgt(i)= (nhl(i-1)+1)*(nhl(i))
      nwtot = nwtot +(nhl(i-1)+1)*nhl(i)
      maxnnode = max(maxnnode,nhl(i))
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

!.....calc number of weights
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

    if( .not.allocated(wgts) ) allocate(wgts(0:maxnnode,maxnnode,nlayer+1))
    wgts(:,:,:) = 0d0
    wgts(:,:,nlayer+1) = 1d0

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