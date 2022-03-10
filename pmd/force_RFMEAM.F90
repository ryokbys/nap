module RFMEAM
!-----------------------------------------------------------------------
!                     Last modified: <2022-03-10 17:37:40 KOBAYASHI Ryo>
!-----------------------------------------------------------------------
!  Parallel implementation of the RF-MEAM pontential.
!  Ref:
!    [1] Lazic and Thijsse, Comp. Mat. Sci. 53, (2012) 483-492
!    [2] Timonova and Thijsse, MSMSE 19 (2011) 015003
!-----------------------------------------------------------------------
  use pmdvars, only: nspmax
  use memory,only: accum_mem
  use util,only: csp2isp
  implicit none
  include "./const.h"
  private
  save

  public:: read_params_RFMEAM, force_RFMEAM, lprmset_RFMEAM
  
  character(len=128):: paramsdir = '.'
  character(len=128),parameter:: paramsfname = 'in.params.RFMEAM'

  integer,parameter:: ioprms = 20
  integer,parameter:: lmax = 3
  real(8),parameter:: tiny = 1d-14
  real(8),parameter:: sgm = 1d-5
  integer,parameter:: allspcs = 1000
  
  logical:: lprmset_RFMEAM = .false.

  
  real(8):: rcmax,rc2max
  real(8):: rcij(nspmax,nspmax),rsij(nspmax,nspmax), &
       rcij2(nspmax,nspmax),trcij2(nspmax,nspmax),trcij(nspmax,nspmax), &
       cmax(nspmax,nspmax,nspmax), cmin(nspmax,nspmax,nspmax), &
       epij(nspmax,nspmax), alpij(nspmax,nspmax), c2ij(nspmax,nspmax), &
       c3ij(nspmax,nspmax), rpij(nspmax,nspmax), pj(0:lmax,nspmax), &
       qj(0:lmax,nspmax), e0(nspmax), e1(nspmax), e2(nspmax), &
       ni(nspmax), ti(lmax,nspmax), rcfac2(nspmax,nspmax)

  logical:: interact(nspmax,nspmax), &
       interact3(nspmax,nspmax,nspmax)


contains
!=======================================================================
  subroutine read_params_RFMEAM(myid_md,mpi_md_world,iprint,specorder)
!
!  Read parameters from file
!--parameter file format------------------------------------------------
!   pairtype   Timonova
!   atomic   cspi  ni  e0  e1  e2  pj(0:lmax)  qj(0:lmax)  ti(lmax)
!   pair     cspi cspj rcij rsij epij alpij c2ij c3ij rpij
!   triple   cspk cspi cspj  cmin cmax
!--or-------------------------------------------------------------------
!   pairtype   Srinivasan
!   atomic   cspi  ni  e0  e1  e2  pj(0:lmax)  qj(0:lmax)  ti(lmax)
!   pair     cspi cspj rcij rsij epij alpij bij(1:mmax) sij(1:mmax)
!   triple   cspk cspi cspj cmin cmax
!-----------------------------------------------------------------------
    include 'mpif.h'
    include './const.h'
    integer,intent(in):: myid_md,mpi_md_world,iprint
    character(len=3),intent(in):: specorder(nspmax)

    integer:: is,js,ks,ierr,l,i,j,k
    integer:: imask(nspmax),jmask(nspmax),kmask(nspmax)
    real(8):: rc,rs,cmaxi,cmini,ep,alp,c2,c3,rp,e0i,e1i,e2i,n, &
         pjl(0:lmax),qjl(0:lmax),til(lmax)
    character(len=128):: cline,fname,c1
    character(len=3):: cspi,cspj,cspk
    logical:: latomic(nspmax)

    if( myid_md.eq.0 ) then
      rcij(:,:) = 0d0
      rsij(:,:) = 0d0
      cmax(:,:,:) = 0d0
      cmin(:,:,:) = 0d0
      epij(:,:) = 0d0
      alpij(:,:) = 0d0
      c2ij(:,:) = 0d0
      c3ij(:,:) = 0d0
      rpij(:,:) = 0d0
      pj(:,:) = 0d0
      qj(:,:) = 0d0
      e0(:) = 0d0
      e1(:) = 0d0
      e2(:) = 0d0
      ni(:) = 0d0
      ti(:,:) = 0d0
      interact(:,:) = .false.
      interact3(:,:,:) = .false.
      latomic(:) = .false.

      fname = trim(paramsdir)//'/'//trim(paramsfname)
      open(ioprms,file=trim(fname),status='old')
      if( iprint.ge.ipl_basic ) write(6,'(/a)') ' RFMEAM parameters:'
      do while(.true.)
        read(ioprms,'(a)',end=10) cline
!.....Check comment lines
        if( cline(1:1).eq.'#' .or. cline(1:1).eq.'!' ) cycle
        if( len_trim(cline).eq.0 ) cycle
        backspace(ioprms)
!.....The 1st entry specifies ATOMIC or PAIR parameters
        read(ioprms,*) c1
        if( trim(c1).eq.'atomic' .or. trim(c1).eq.'element' ) then
          backspace(ioprms)
          read(ioprms,*) c1,cspi, n,e0i,e1i,e2i,pjl(0:lmax),qjl(0:lmax),til(1:lmax)
          is = csp2isp(trim(cspi))
          if( is.lt.1 .and. trim(cspi).eq.'*' ) then  ! '*' means all.
            is = allspcs
          else if( is.gt.nspmax .or. is.lt.1 ) then
            print *,'Warning @read_params_RFMEAM: IS not in range [1,nspmax]: ',&
                 trim(cspi)//'=',is
            cycle
          endif
          if( is.eq.allspcs ) then
            latomic(:) = .true.
            ni(:) = n
            e0(:) = e0i
            e1(:) = e1i
            e2(:) = e2i
            do i=1,nspmax
              pj(:,i) = pjl(:)
              qj(:,i) = qjl(:)
              ti(:,i) = til(:)
            enddo
          else
            latomic(is) = .true.
            ni(is) = n
            e0(is) = e0i
            e1(is) = e1i
            e2(is) = e2i
            pj(:,is) = pjl(:)
            qj(:,is) = qjl(:)
            ti(:,is) = til(:)
          endif

!.....Otherwise the entry is for pair parameter
        else if( trim(c1).eq.'pair' ) then
          backspace(ioprms)
          read(ioprms,*) c1,cspi,cspj, rc,rs,ep,alp,c2,c3,rp
          is = csp2isp(trim(cspi))
          js = csp2isp(trim(cspj))
          if( is.gt.nspmax .or. (is.lt.1.and.trim(cspi).ne.'*') &
               .or. js.gt.nspmax .or. (js.lt.1.and.trim(cspj).ne.'*') ) then
            print *,'Warning @read_params_RFMEAM: IS/JS is not in range [1,nspmax]: ',&
                 trim(cspi)//','//trim(cspj)//'=',is,js
            cycle
          endif
!.....Set imask,jmask
          imask(:) = 1
          jmask(:) = 1
          do i=1,nspmax
            if( i.eq.is ) imask(i) = 0  ! unmask is
            if( i.eq.js ) jmask(i) = 0  ! unmask js
          enddo
          if( trim(cspi).eq.'*' ) imask(:) = 0  ! unmask all
          if( trim(cspj).eq.'*' ) jmask(:) = 0  ! unmask all
!.....Set params using the masks
          do is=1,nspmax
            if( imask(is).eq.1 ) cycle
            do js=1,nspmax
              if( jmask(js).eq.1 ) cycle
              interact(is,js) = .true.
              rcij(is,js) = rc
              rsij(is,js) = rs
              epij(is,js) = ep
              alpij(is,js) = alp
              c2ij(is,js) = c2
              c3ij(is,js) = c3
              rpij(is,js) = rp
!.....Symmetrize parameters
              interact(js,is) = .true.
              rcij(js,is) = rc
              rsij(js,is) = rs
              epij(js,is) = ep
              alpij(js,is) = alp
              c2ij(js,is) = c2
              c3ij(js,is) = c3
              rpij(js,is) = rp
            enddo
          enddo
        else if( c1(1:6).eq.'triple' ) then
          backspace(ioprms)
          read(ioprms,*) c1,cspk,cspi,cspj, cmini,cmaxi
          is = csp2isp(trim(cspi))
          js = csp2isp(trim(cspj))
          ks = csp2isp(trim(cspk))
          if( is.gt.nspmax .or. (is.lt.1.and.trim(cspi).ne.'*') .or. &
               js.gt.nspmax .or. (js.lt.1.and.trim(cspj).ne.'*') .or. &
               ks.gt.nspmax .or. (ks.lt.1.and.trim(cspk).ne.'*') ) then
            print *,'Warning @read_params_RFMEAM: KS/IS/JS is not in range [1,nspmax]: ',&
                 trim(cspk)//' '//trim(cspi)//' '//trim(cspj)
            cycle
          endif
!.....Set masks
          imask(:) = 1
          jmask(:) = 1
          kmask(:) = 1
          do i=1,nspmax
            if( i.eq.is ) imask(i) = 0  ! unmask is
            if( i.eq.js ) jmask(i) = 0  ! unmask js
            if( i.eq.ks ) kmask(i) = 0  ! unmask ks
          enddo
          if( trim(cspi).eq.'*' ) imask(:) = 0  ! unmask all
          if( trim(cspj).eq.'*' ) jmask(:) = 0  ! unmask all
          if( trim(cspk).eq.'*' ) kmask(:) = 0  ! unmask all
!.....Set params using the masks
          do is=1,nspmax
            if( imask(is).eq.1 ) cycle
            do js=1,nspmax
              if( jmask(js).eq.1 ) cycle
              do ks=1,nspmax
                if( kmask(ks).eq.1 ) cycle
                interact3(ks,is,js) = .true.
                cmax(ks,is,js) = cmaxi
                cmin(ks,is,js) = cmini
!.....Symmetrize
                interact3(ks,js,is) = .true.
                cmax(ks,js,is) = cmaxi
                cmin(ks,js,is) = cmini
              enddo
            enddo
          enddo
        endif  ! c1 
      enddo
10    continue
      close(ioprms)

!.....TODO: Check consistency

!.....If is-js parameters are not specified even though is-is and js-js
!     parameters are given, is-js parameters are computed as averages.
      do is=1,nspmax
        if( .not. interact(is,is) ) cycle
        do js=is+1,nspmax
          if( .not. interact(js,js) .or. interact(is,js) ) cycle
          rcij(is,js) = (rcij(is,is)+rcij(js,js))/2
          rsij(is,js) = (rsij(is,is)+rsij(js,js))/2
          epij(is,js) = (epij(is,is)+epij(js,js))/2
          alpij(is,js)= (alpij(is,is)+alpij(js,js))/2
          c2ij(is,js) = (c2ij(is,is)+c2ij(js,js))/2
          c3ij(is,js) = (c3ij(is,is)+c3ij(js,js))/2
!.....Symmetrize
          rcij(js,is) = rcij(is,js) 
          rsij(js,is) = rsij(is,js) 
          epij(js,is) = epij(is,js) 
          alpij(js,is)= alpij(is,js)
          c2ij(js,is) = c2ij(is,js) 
          c3ij(js,is) = c3ij(is,js)
        enddo
      enddo

      if( iprint.ge.ipl_basic ) then
        do is=1,nspmax
          if( latomic(is) .and. trim(specorder(is)).ne.'x' ) then
            cspi = trim(specorder(is))
            write(6,'(a,1x,a3,4(2x,4f8.3))') '   atomic  ',trim(cspi), &
                 ni(is),e0(is),e1(is),e2(is), &
                 pj(0:lmax,is),qj(0:lmax,is),ti(1:lmax,is)
          endif
        enddo
        do is=1,nspmax
          cspi = trim(specorder(is))
          if( cspi.eq.'x' ) cycle
          do js=is,nspmax
            cspj = trim(specorder(js))
            if( cspj.eq.'x' ) cycle
            if( interact(is,js) ) then
              write(6,'(a,2(1x,a3),7f8.3)') '   pair  ',trim(cspi),trim(cspj), &
                   rcij(is,js),rsij(is,js),epij(is,js),alpij(is,js), &
                   c2ij(is,js),c3ij(is,js),rpij(is,js)
            endif
          enddo
        enddo
        do is=1,nspmax
          cspi = trim(specorder(is))
          if( cspi.eq.'x' ) cycle
          do js=is,nspmax
            cspj = trim(specorder(js))
            if( cspj.eq.'x' ) cycle
            do ks=1,nspmax
              cspk = trim(specorder(ks))
              if( cspk.eq.'x' ) cycle
              if( .not. interact3(ks,is,js) ) cycle
              write(6,'(a,a2,"-",a2,"-",a2,2f7.3)') '   triple  ', &
                   trim(cspk),trim(cspi),trim(cspj), &
                   cmin(ks,is,js),cmax(ks,is,js)
            enddo
          enddo
        enddo
      endif
    endif  ! myid_md == 0

    call mpi_bcast(rcij,nspmax**2,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(rsij,nspmax**2,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(epij,nspmax**2,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(alpij,nspmax**2,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(c2ij,nspmax**2,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(c3ij,nspmax**2,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(rpij,nspmax**2,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(pj,4*nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(qj,4*nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(ti,3*nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(ni,nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(e0,nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(e1,nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(e2,nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(cmax,nspmax**3,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(cmin,nspmax**3,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(interact,nspmax**2,mpi_logical,0,mpi_md_world,ierr)
    call mpi_bcast(interact3,nspmax**3,mpi_logical,0,mpi_md_world,ierr)
    call mpi_bcast(lmax,1,mpi_integer,0,mpi_md_world,ierr)

    if( myid_md.eq.0 .and. iprint.gt.ipl_basic ) then
      write(6,'(a)') ' Finished reading '//trim(fname)
    endif

    return
  end subroutine read_params_RFMEAM
!=======================================================================
  subroutine force_RFMEAM(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rcg,lspr,d2lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,lstrs,iprint,l1st)
    implicit none
    include "mpif.h"
    include "./params_unit.h"
    integer,intent(in):: namax,natm,nnmax,nismax,lspr(0:nnmax,namax)&
         ,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_md_world,myid_md,nex(3)
    real(8),intent(in):: ra(3,namax),h(3,3),hi(3,3),sv(3,6) &
         ,rcg,tag(namax),d2lspr(nnmax,namax)
    real(8),intent(inout):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: l1st
    logical:: lstrs

    integer:: i,j,k,l,m,n,ierr,is,js,ks,ixyz,jxyz,jj,kk,nni
    real(8):: xi(3),xij(3),rij(3),xj(3),dij2,dij,ep,alpha,c2,c3, &
         dp,eta,expeta,phi,rs,rc,zij,z2,fcij,dfcij,tmp,epotl, &
         driji(3),drijj(3),dphi,dtmp,phifc,eqjr,dfcdr, &
         rik(3),dik2,dik,fcik,cs,cs2,plcs(0:lmax),dplcs(0:lmax),sfcjk, &
         rhoi2(0:lmax),gam,egam,ggam,dgdgam,rhoi,yi,gyi,fyi,frhoi,rhoi0, &
         pcs,pcsi,dcsdij(3),dcsdik(3),strho2,dgdy,dfdy,atmp(3), &
         cmaxkij,cmaxmax,truerc,epott,epot2l,epotml,epot2,epotm
    real(8),allocatable,save:: aal(:,:),strsl(:,:,:),epil(:)
    real(8),allocatable,save:: sij(:),dsij(:,:),sfc(:),fl(:,:), &
         dfl(:,:,:),dsfc(:,:,:),drhoi2(:,:,:),drhoi0(:,:),dstrho2(:,:), &
         dgam(:,:),drho(:,:),rijs(:,:)

    if( l1st ) then
      if( allocated(aal) ) then
        deallocate(aal,strsl,epil,sij,dsij,sfc,fl,dfl,dsfc, &
             drhoi2,drhoi0,dstrho2,dgam,drho,rijs)
      endif
      allocate(aal(3,namax),strsl(3,3,namax),epil(namax),sij(nnmax), &
           dsij(3,nnmax),sfc(nnmax),fl(0:lmax,nnmax), &
           dfl(3,0:lmax,nnmax),dsfc(3,nnmax,nnmax),drhoi2(3,nnmax,lmax), &
           drhoi0(3,nnmax),dstrho2(3,nnmax),dgam(3,nnmax),drho(3,nnmax), &
           rijs(5,nnmax))
      call accum_mem('force_RFMEAM',8*(size(aal) +size(strsl) +size(epil) &
           +size(sij) +size(dsij) +size(sfc) +size(dsfc) &
           +size(fl) +size(dfl) +size(drhoi2) +size(drhoi0) +size(dstrho2) &
           +size(dgam) +size(drho) +size(rijs)))
!.....True rcut = max(rc, rc*Cmax/(2*sqrt(Cmax -1)))
      cmaxmax = 0d0
      rcmax = 0d0
      rcij2(:,:) = 0d0
      trcij2(:,:) = 0d0
      rcfac2(:,:) = 1d0
      do is=1,nspmax
        do js=1,nspmax
          if( .not.interact(is,js) ) cycle
          rc = rcij(is,js)
          do ks=1,nspmax
            cmaxkij = cmax(ks,is,js)
            if( cmaxkij.le.1d0 ) cycle
            truerc = max(rc, rc*cmaxkij/(2d0*sqrt(cmaxkij -1d0)))
            rcfac2(is,js) = max(rcfac2(is,js), &
                 cmaxkij**2/(4d0*(cmaxkij -1d0)))
          enddo
          trcij(is,js) = truerc
          rcij2(is,js) = rc*rc
          trcij2(is,js) = truerc*truerc
          rcmax = max(rcmax,truerc)
        enddo
      enddo
      if( rcg.lt.rcmax ) then
        if( myid_md.eq.0 ) then
          write(6,*) " Error: Global rc is smaller than one of RFMEAM rc's."
        endif
        call mpi_finalize(ierr)
        stop
      endif
      rc2max = rcmax*rcmax
    endif

    if( size(aal).lt.namax ) then
      deallocate(aal,strsl,epil,sij,dsij,sfc,fl,dfl,dsfc, &
           drhoi2,drhoi0,dstrho2,dgam,drho,rijs)
      allocate(aal(3,namax),strsl(3,3,namax),epil(namax),sij(nnmax), &
           dsij(3,nnmax),sfc(nnmax),fl(0:lmax,nnmax), &
           dfl(3,0:lmax,nnmax),dsfc(3,nnmax,nnmax),drhoi2(3,nnmax,lmax), &
           drhoi0(3,nnmax),dstrho2(3,nnmax),dgam(3,nnmax),drho(3,nnmax), &
           rijs(5,nnmax))
    endif

    epotl= 0d0
    epot2l= 0d0
    epotml= 0d0
    epil(:) = 0d0
    aal(:,:) = 0d0
    strsl(:,:,:) = 0d0

    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is = int(tag(i))
      sij(:) = 0d0
      dsij(:,:) = 0d0
      dsfc(:,:,:) = 0d0
      sfc(:) = 0d0
      fl(:,:) = 0d0
      dfl(:,:,:) = 0d0
!.....Create ij vector and distances and store them for after heavy use
      rijs(:,:) = 0d0
      nni = lspr(0,i)
      do jj=1,nni
!!$        if( d2lspr(jj,i).ge.rc2max ) cycle
        j = lspr(jj,i)
        xj(1:3) = ra(1:3,j)
        xij(1:3) = xj(1:3) -xi(1:3)
        rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij2 = rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
        dij = sqrt(dij2)
!.....Store
        rijs(1:3,jj) = rij(1:3)
        rijs(4,jj) = dij2
        rijs(5,jj) = dij
      enddo
!.....Sij and pair potential
      do jj=1,nni
        j = lspr(jj,i)
        js = int(tag(j))
        if( .not. interact(is,js) ) cycle
        dij2 = rijs(4,jj)
        if( dij2.gt.rcij2(is,js) ) cycle
        rij(1:3) = rijs(1:3,jj)
        dij = rijs(5,jj)
        call compute_sij(i,j,jj,is,js,rijs, &
             namax,natm,nnmax,tag,lspr,sij(jj),dsij(:,:))
!.....Skip computing i-j interaction if it is completely screened (Sij==0).
        if( sij(jj).lt.tiny ) cycle
!-----Pair potential part ------------------------------------
        ep = epij(is,js)
        alpha = alpij(is,js)
        c2 = c2ij(is,js)
        c3 = c3ij(is,js)
        dp = rpij(is,js)
!.....Cutoff function, fc
        fcij = fcut(is,js,dij)
        dfcij = dfcut(is,js,dij)
        if( fcij.lt.tiny ) cycle
        eta = alpha *(dij/dp -1d0)
        expeta = exp(-eta)
        phi = -ep*(1d0 +eta +c2*eta**2 +c3*eta**3) *expeta
        tmp = 0.5d0 *phi*fcij*sij(jj)
        epil(i) = epil(i) +0.5d0 *tmp
        epil(j) = epil(j) +0.5d0 *tmp
        epot2l = epot2l +tmp
!.....Forces related to pair potential, which also depends on k != i,j
        driji(1:3) = -rij(1:3)/dij
        drijj(1:3) = -driji(1:3)
        dphi = alpha/dp *ep*expeta *(eta*(1d0 -2d0*c2) &
             +eta**2 *(c2 -3d0*c3) +c3*eta**3)
        dtmp = 0.5d0 *(dphi*sij(jj)*fcij +phi*sij(jj)*dfcij)
        aal(1:3,j) = aal(1:3,j) -dtmp*drijj(1:3)
        aal(1:3,i) = aal(1:3,i) -dtmp*driji(1:3)
        phifc = 0.5d0 *phi *fcij
        do kk=1,nni
          k = lspr(kk,i)
          aal(1:3,k) = aal(1:3,k) -phifc*dsij(1:3,kk)
          aal(1:3,i) = aal(1:3,i) +phifc*dsij(1:3,kk)
        enddo
!.....Atomic stress for pair part
        if( lstrs ) then
          do ixyz=1,3
            do jxyz=1,3
              strsl(jxyz,ixyz,i)=strsl(jxyz,ixyz,i) &
                   -0.5d0*dtmp*rij(ixyz)*(-driji(jxyz))
              strsl(jxyz,ixyz,j)=strsl(jxyz,ixyz,j) &
                   -0.5d0*dtmp*rij(ixyz)*(-driji(jxyz))
            enddo
          enddo
          do kk=1,nni
            k = lspr(kk,i)
            do ixyz=1,3
              do jxyz=1,3
                strsl(jxyz,ixyz,i)=strsl(jxyz,ixyz,i) &
                     -0.5d0*phifc*dsij(jxyz,kk)
                strsl(jxyz,ixyz,k)=strsl(jxyz,ixyz,k) &
                     -0.5d0*phifc*dsij(jxyz,kk)
              enddo
            enddo
          enddo  ! kk-loop
        endif  ! lstrs
!.....Compute Sij*fcij,flij and their derivatives for embedded part
        sfc(jj) = sij(jj)*fcij
        do l=0,lmax
          eqjr = exp(-qj(l,js)*dij)
          fl(l,jj) = pj(l,js) *eqjr
          dfl(1:3,l,jj)= pj(l,js)*(-qj(l,js))*eqjr *rij(1:3)/dij
        enddo
!.....dsfc wrt non-j
        do kk=1,nni
          if( kk.eq.jj) cycle
          k = lspr(kk,i)
          ks = int(tag(k))
          if( .not.interact3(ks,is,js) ) cycle
          dik2 = rijs(4,kk)
          if( dik2.gt.trcij2(is,js) ) cycle
          dsfc(1:3,kk,jj) = dsfc(1:3,kk,jj) +dsij(1:3,kk)*fcij
        enddo
!.....dsfc wrt j
        dsfc(1:3,jj,jj) = dsfc(1:3,jj,jj) +dsij(1:3,jj)*fcij &
             +sij(jj)*dfcij*rij(1:3)/dij
      enddo ! jj-loop

      rhoi2(:) = 0d0
      do jj=1,nni
        if( sfc(jj).lt.tiny ) cycle
        j = lspr(jj,i)
        js = int(tag(j))
        rij(1:3) = rijs(1:3,jj)
        dij2 = rijs(4,jj)
        dij = rijs(5,jj)
!.....Compute rhoi2(l)
        do kk=1,nni
          if( sfc(kk).lt.tiny ) cycle
          k = lspr(kk,i)
          ks = int(tag(k))
          if( .not.interact3(ks,is,js) ) cycle
          dik2 = rijs(4,kk)
          if( dik2.gt.rcij2(is,ks) ) cycle
          rik(1:3) = rijs(1:3,kk)
          dik = rijs(5,kk)
          cs = (rij(1)*rik(1) +rij(2)*rik(2) +rij(3)*rik(3))/dij/dik
          cs2 = cs*cs
          plcs(0) = 1d0
          plcs(1) = cs
          plcs(2) = (3d0*cs2 -1d0)/2
          plcs(3) = (5d0*cs*cs2 -3d0*cs)/2
          sfcjk = sfc(jj)*sfc(kk)
          do l=0,lmax
            rhoi2(l) = rhoi2(l) +sfcjk*fl(l,jj)*fl(l,kk)*plcs(l)
          enddo
        enddo  ! kk-loop
      enddo  ! jj-loop

!.....F[rhoi]
      rhoi0 = sqrt(rhoi2(0))
      gam = 0d0
      do l=1,lmax
        gam = gam +ti(l,is)*rhoi2(l)
      enddo
      gam = gam /rhoi2(0)
      egam = exp(-gam)
      ggam = 2d0 /(1d0 +egam)
      rhoi = rhoi0 *ggam
      yi = rhoi/ni(is)
      gyi = 1d0 -exp(-yi*yi /2d0/sgm**2)
      fyi = e0(is)*yi*log(yi) +e1(is)*yi +e2(is)*yi*yi
      frhoi = fyi*gyi
      epil(i)= epil(i) +frhoi
      epotml = epotml +frhoi

!.....d(rhoi(l)^2)/dr_{il} needed for dF[rhoi]/dr_{il}
      drhoi2(:,:,:) = 0d0
      drhoi0(:,:) = 0d0
      do jj=1,nni
        if( sij(jj).lt.tiny ) cycle
        j = lspr(jj,i)
        js = int(tag(j))
        if( .not.interact(is,js) ) cycle
        dij2 = rijs(4,jj)
        if( dij2.gt.rcij2(is,js) ) cycle
        rij(1:3) = rijs(1:3,jj)
        dij = rijs(5,jj)
        fcij = fcut(is,js,dij)
        do kk=1,nni
          k = lspr(kk,i)
          ks = int(tag(k))
          if( .not.interact3(ks,is,js) ) cycle
          dik2 = rijs(4,kk)
          if( dik2.gt.rcij2(is,ks) ) cycle  ! NOTE: cutoff for triplet is betw is-ks
          rik(1:3) = rijs(1:3,kk)
          dik = rijs(5,kk)
          fcik = fcut(is,ks,dik)
          pcs = dij*dik
          pcsi = 1d0/pcs
          cs = (rij(1)*rik(1) +rij(2)*rik(2) +rij(3)*rik(3)) *pcsi
          dcsdij(1:3) = (rik(1:3) -cs*dik/dij*rij(1:3)) *pcsi
          dcsdik(1:3) = (rij(1:3) -cs*dij/dik*rik(1:3)) *pcsi
          sfcjk = sfc(jj)*sfc(kk)
          cs2 = cs*cs
          plcs(1) = cs
          plcs(2) = (3d0*cs2 -1d0)/2
          plcs(3) = (5d0*cs*cs2 -3d0*cs)/2
          dplcs(1)= 1d0
          dplcs(2)= 3d0*cs
          dplcs(3)= (15d0*cs2 -3d0)/2
          do l=1,lmax
            drhoi2(1:3,1:nni,l) = drhoi2(1:3,1:nni,l) +fl(l,jj)*fl(l,kk)*plcs(l) &
                 *(sfc(jj)*dsfc(1:3,1:nni,kk) +sfc(kk)*dsfc(1:3,1:nni,jj))
            drhoi2(1:3,jj,l) = drhoi2(1:3,jj,l) +sfcjk*fl(l,kk) &
                 *(plcs(l)*dfl(1:3,l,jj) +fl(l,jj)*dcsdij(1:3)*dplcs(l))
            drhoi2(1:3,kk,l) = drhoi2(1:3,kk,l) +sfcjk*fl(l,jj) &
                 *(plcs(l)*dfl(1:3,l,kk) +fl(l,kk)*dcsdik(1:3)*dplcs(l))
          enddo
          drhoi0(1:3,1:nni) = drhoi0(1:3,1:nni) +0.5d0/rhoi0 &
               *fl(0,jj)*fl(0,kk) &
               *(dsfc(1:3,1:nni,jj)*sfc(kk) +sfc(jj)*dsfc(1:3,1:nni,kk))
          drhoi0(1:3,jj) = drhoi0(1:3,jj) +0.5d0/rhoi0*dfl(1:3,0,jj)*fl(0,kk)*sfcjk
          drhoi0(1:3,kk) = drhoi0(1:3,kk) +0.5d0/rhoi0*dfl(1:3,0,kk)*fl(0,jj)*sfcjk
        enddo  ! kk-loop
      enddo  ! jj-loop

!.....dF[rhoi]
      strho2 = 0d0
      dstrho2(:,:) = 0d0
      do l=1,lmax
        strho2 = strho2 +ti(l,is)*rhoi2(l)
        dstrho2(1:3,1:nni)= dstrho2(1:3,1:nni) +ti(l,is)*drhoi2(1:3,1:nni,l)
      enddo
      dgam(1:3,1:nni) = -2d0/rhoi0**3 *strho2 *drhoi0(1:3,1:nni) &
           +1d0/rhoi2(0) *dstrho2(1:3,1:nni)
      dgdgam = 2d0*egam/(1d0+egam)**2
      drho(1:3,1:nni) = drhoi0(1:3,1:nni)*ggam &
           +rhoi0 *dgdgam *dgam(1:3,1:nni)
      dgdy = yi /sgm**2 *exp(-yi*yi/2/sgm**2)
      dfdy = (e0(is)*log(yi) +e0(is) +e1(is) +2d0*e2(is)*yi)*gyi +dgdy*fyi
      do jj=1,nni
!!$        if( sij(jj).lt.tiny ) cycle
        j = lspr(jj,i)
        js = int(tag(j))
        if( .not. interact(is,js) ) cycle
        dij2 = rijs(4,jj)
        if( dij2.gt.trcij2(is,js) ) cycle
        rij(1:3) = rijs(1:3,jj)
        atmp(1:3) = dfdy /ni(is) *drho(1:3,jj)
        aal(1:3,j) = aal(1:3,j) -atmp(1:3)
        aal(1:3,i) = aal(1:3,i) +atmp(1:3)
        if( lstrs ) then
          do ixyz=1,3
            do jxyz=1,3
              strsl(jxyz,ixyz,i)=strsl(jxyz,ixyz,i) &
                   -0.5d0*rij(ixyz)*atmp(jxyz)
              strsl(jxyz,ixyz,j)=strsl(jxyz,ixyz,j) &
                   -0.5d0*rij(ixyz)*atmp(jxyz)
            enddo
          enddo
        endif
      enddo ! jj-loop
!!$      stop
    enddo  ! i-loop

!.....Send back epi on immigrants
    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_md_world,epil,1)
    epi(1:natm) = epi(1:natm) +epil(1:natm)
    
!.....Send back forces on immigrants
    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_md_world,aal,3)
    aa(1:3,1:natm)= aa(1:3,1:natm) +aal(1:3,1:natm)
    
    if( lstrs ) then
!.....Send back stresses on immigrants
      call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
           ,nn,mpi_md_world,strsl,9)
      strs(1:3,1:3,1:natm)= strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)
    endif

!.....Gather epot
    epotl = epot2l +epotml
    call mpi_allreduce(epotl,epott,1,mpi_real8,mpi_sum,mpi_md_world,ierr)
    epot= epot +epott
    if( iprint.ge.ipl_info ) then
      call mpi_allreduce(epot2l,epot2,1,mpi_real8,mpi_sum,mpi_md_world,ierr)
      call mpi_allreduce(epotml,epotm,1,mpi_real8,mpi_sum,mpi_md_world,ierr)
      if( myid_md.eq.0 ) &
           write(6,'(a,3es15.7)') ' epot RFMEAM (2-body,embed,total = ', &
           epot2,epotm,epott
    endif
    return
  end subroutine force_RFMEAM
!=======================================================================
  subroutine compute_sij(i,j,jj,is,js,rijs, &
       namax,natm,nnmax,tag,lspr,sij,dsij)
!
!  Compute Sij and its derivatives wrt r_(i,x).
!  
!  Sij = prod_k B(y_{jik}), and y_{jik} = (C_jik -C_min)/(C_max -C_min)
!  Ckij = (1-cos**2)/((rij/rik -cos)*cos)
!  Use above Ckij definition (using cos) instead of the def with rij,rik,rjk,
!  since that causes j-k action-reaction that makes it complicated when
!  implementing stress calculation.
!
    implicit none
    integer,intent(in):: i,j,jj,is,js,namax,natm,nnmax
    integer,intent(in):: lspr(0:nnmax,namax)
    real(8),intent(in):: rijs(5,nnmax)
    real(8),intent(in):: tag(namax)
    real(8),intent(out):: sij,dsij(3,nnmax)

    logical,save:: l1st = .true.
    real(8),allocatable,save:: skij(:)

    integer:: kk,k,ks
    real(8):: dij4,cmaxkij,cminkij,dij,dij2,rij(3)
    real(8):: xk(3),xik(3),rik(3),dik2,dik
    real(8):: denom,numer,ckij,y,denom2
    real(8):: dc,driki(3),drikk(3),pcs,pcsi,qcs,cs,sn, &
         ddij,ddik,dcsdij(3),dcsdik(3), &
         dcdij(3),dcdik(3),dbdy,tmp,dnumdcs,ddendcs,dcdcs, &
         dydc,sijperkij

    if( l1st ) then
      allocate(skij(nnmax))
      l1st = .false.
    endif

!!$    print *,' compute_...'
    rij(1:3) = rijs(1:3,jj)
    dij2 = rijs(4,jj)
    dij = rijs(5,jj)
    dij4 = dij2*dij2
    sij = 1d0
    skij(:) = 0d0
    dsij(:,:) = 0d0
!.....Firstly compute all the Skij
    do kk=1,lspr(0,i)
      dik2 = rijs(4,kk)
      if( dik2.gt.dij2*rcfac2(is,js) ) cycle  ! NOTE: cutoff wrt dij
      k= lspr(kk,i)
      ks = int(tag(k))
      if( k.eq.j ) cycle
      if( .not. interact3(ks,is,js) ) cycle
      cmaxkij = cmax(ks,is,js)
      cminkij = cmin(ks,is,js)
      rik(1:3) = rijs(1:3,kk)
      dik = rijs(5,kk)
      cs = (rij(1)*rik(1) +rij(2)*rik(2) +rij(3)*rik(3))/dij/dik
      numer = 1d0 -cs*cs
      denom = (dij/dik -cs)*cs
      ckij = numer/denom
      y = (ckij -cminkij)/(cmaxkij -cminkij)
      if( denom.le.tiny .or. y.ge.1d0 ) then
        skij(kk) = 1d0
      else if( y.le.tiny ) then
        skij(kk) = 0d0
        sij = 0d0
        exit
      else  ! case: 0 < y < 1
        skij(kk) = (1d0 -(1d0-y)**4)**2
      endif
      sij = sij *skij(kk)
    enddo  ! kk
!.....If sij==0d0 (skij==0d0 exists), skip derivative calculation
    if( sij.lt.tiny ) return
!.....Secondly, compute dSkij/dr_(i,x), dSkij/dr_(j,y), ...
    do kk=1,lspr(0,i)
      if( skij(kk).gt.1d0-tiny .or. skij(kk).lt.tiny ) cycle
      k= lspr(kk,i)
      if( k.eq.j ) cycle
      ks = int(tag(k))
      if( .not. interact3(ks,is,js) ) cycle
      dik2 = rijs(4,kk)
      if( dik2.gt.dij2*rcfac2(is,js) ) cycle  ! NOTE: cutoff wrt dij
! Only the cases ( 0 <= y <= 1 ) survive hereafter
      cmaxkij = cmax(ks,is,js)
      cminkij = cmin(ks,is,js)
      rik(1:3) = rijs(1:3,kk)
      dik = rijs(5,kk)
      pcs = dij*dik
      pcsi= 1d0/pcs
      qcs = rij(1)*rik(1) +rij(2)*rik(2) +rij(3)*rik(3)
      cs = qcs *pcsi
      numer = 1d0 -cs*cs
      denom = (dij/dik -cs)*cs
      ckij = numer/denom
      denom2 = denom*denom
      dcsdij(1:3) = (rik(1:3) -cs*dik/dij*rij(1:3)) *pcsi
      dcsdik(1:3) = (rij(1:3) -cs*dij/dik*rik(1:3)) *pcsi
      dcdij(1:3) = -(1d0-cs*cs)/denom2 *cs /dik *rij(1:3)/dij
      dcdik(1:3) =  (1d0-cs*cs)/denom2 *cs /dik/dik *dij *rik(1:3)/dik
      dnumdcs = -2d0 *cs
      ddendcs = dij/dik -2d0*cs
      dcdcs = (dnumdcs*denom -numer*ddendcs)/denom2
      dcdij(1:3) = dcdij(1:3) +dcdcs*dcsdij(1:3)
      dcdik(1:3) = dcdik(1:3) +dcdcs*dcsdik(1:3)
      y = (ckij -cminkij)/(cmaxkij -cminkij)
      dydc = 1d0 /(cmaxkij -cminkij)
      dbdy = 8d0*(1d0 -(1d0-y)**4) *(1d0-y)**3
      sijperkij = sij/skij(kk)
      tmp = dbdy *dydc *sijperkij
      dsij(1:3,jj) = dsij(1:3,jj) +tmp*dcdij(1:3)
      dsij(1:3,kk) = dsij(1:3,kk) +tmp*dcdik(1:3)
    enddo  ! kk
    return
  end subroutine compute_sij
!=======================================================================
  function fcut(is,js,rij)
    integer,intent(in):: is,js
    real(8),intent(in):: rij
    real(8):: fcut
    real(8):: zij,rs,rc,z2
    
    rs = rsij(is,js)
    rc = rcij(is,js)
    zij = (rij-rs)/(rc-rs)
    if(zij.lt.0d0) then
      fcut = 1d0
    else if(zij.ge.1d0 ) then
      fcut = 0d0
    else
      z2 = zij*zij
      fcut = 1d0 -z2*zij*(6d0*z2 -15d0*zij +10d0)
    endif
    return
  end function fcut
!=======================================================================
  function dfcut(is,js,rij)
    integer,intent(in):: is,js
    real(8),intent(in):: rij
    real(8):: dfcut
    real(8):: zij,rs,rc,z2
    
    rs = rsij(is,js)
    rc = rcij(is,js)
    zij = (rij-rs)/(rc-rs)
    if(zij.lt.0d0) then
      dfcut = 0d0
    else if(zij.ge.1d0 ) then
      dfcut = 0d0
    else
      z2 = zij*zij
      dfcut = -30d0*(z2*z2 -2d0*z2*zij +z2)/(rc -rs)
    endif
    return
  end function dfcut
!=======================================================================
end module RFMEAM
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd lib"
!     End:
