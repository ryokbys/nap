module tersoff
!-----------------------------------------------------------------------
!                     Last modified: <2021-12-23 14:05:18 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
! Ref:
!   [1] Tersoff, Physical Review B, 38(14), 9902–9905 (1988).
!   [2] Kumagai, et al., Comput. Mater. Sci. 39, 457 (2007).
!   [3] Murty, et al., Physical Review B, 51(8), 4889–4893 (1995).
!   [4] https://lammps.sandia.gov/doc/pair_tersoff_mod.html#tersoff-12
!   [5] Shokeen, & Schelling, (1995).IEEE Transactions on Microwave Theory and Techniques, 43(8), 1826–1833
!-----------------------------------------------------------------------
  use pmdvars,only: nspmax, nsp
  use util,only: csp2isp
  implicit none
  include "./const.h"
  save
  integer,parameter:: ioprms = 50
  character(len=128):: paramsdir = '.'
  character(len=128),parameter:: paramsfname = 'in.params.Tersoff'

!.....global variables
  integer:: mem = 0
  real(8):: time = 0d0
  real(8):: ts_epot

!.....Tersoff type: (default) Kumagai's modified Tersoff, (Te-dependent) Te-dependent[5]
  character(len=12):: ts_type = 'default'

!.....Cutoff type: (1) Tersoff,  (2) Kumagai, (3) Murty (default)
  integer:: ts_fc_type = 3

!.....Original parameters from [2]
  real(8):: ts_a(nspmax,nspmax), ts_b(nspmax,nspmax), ts_lmbd1(nspmax,nspmax), &
       ts_lmbd2(nspmax,nspmax), ts_eta(nspmax), ts_delta(nspmax), &
       ts_alpha(nspmax,nspmax), ts_beta(nspmax), ts_c1(nspmax), &
       ts_c2(nspmax), ts_c3(nspmax), ts_c4(nspmax), ts_c5(nspmax), ts_h(nspmax), &
       ts_r1(nspmax,nspmax), ts_r2(nspmax,nspmax), ts_f0(nspmax), &
       ts_r(nspmax,nspmax), ts_d(nspmax,nspmax), &
       ts_rc(nspmax,nspmax), ts_rc2(nspmax,nspmax)

!.....Te-dependent parameters
  integer:: ntemp
  real(8),allocatable:: ted_te(:),ted_a(:),ted_b(:),ted_lmbd1(:),ted_lmbd2(:) &
       ,ted_eta(:),ted_delta(:),ted_alpha(:),ted_beta(:) &
       ,ted_h(:),ted_r1(:),ted_r2(:),ted_f0(:) &
       ,ted_c1(:),ted_c2(:),ted_c3(:),ted_c4(:),ted_c5(:)

  
contains
  subroutine init_tersoff(myid,mpi_world,iprint,specorder)
    integer,intent(in):: myid,mpi_world,iprint
    character(len=3),intent(in):: specorder(nspmax)

    integer:: isp,jsp

    call read_params_tersoff(myid,mpi_world,iprint,specorder)

!.....Lorentz-Berthelot rule for inter-species parameters
    do isp=1,nspmax
      do jsp=1,nspmax
        ts_a(isp,jsp) = sqrt(ts_a(isp,isp)*ts_a(jsp,jsp))
        ts_b(isp,jsp) = sqrt(ts_b(isp,isp)*ts_b(jsp,jsp))
        ts_lmbd1(isp,jsp)= (ts_lmbd1(isp,isp)+ts_lmbd1(jsp,jsp))/2
        ts_lmbd2(isp,jsp)= (ts_lmbd2(isp,isp)+ts_lmbd2(jsp,jsp))/2
        ts_r1(isp,jsp)= (ts_r1(isp,isp)+ts_r1(jsp,jsp))/2
        ts_r2(isp,jsp)= (ts_r2(isp,isp)+ts_r2(jsp,jsp))/2
!.....R, D, Rc
        ts_r(isp,jsp)= (ts_r1(isp,jsp)+ts_r2(isp,jsp))/2
        ts_d(isp,jsp)= (ts_r2(isp,jsp)-ts_r1(isp,jsp))/2
        ts_rc(isp,jsp)= ts_r(isp,jsp) +ts_d(isp,jsp)
        ts_rc2(isp,jsp)= ts_rc(isp,jsp)**2
      enddo
    enddo
!!$    ts_r = (ts_r1 +ts_r2)/2
!!$    ts_d = (ts_r2 -ts_r1)/2
!!$    ts_rc = ts_r +ts_d
!!$    ts_rc2= ts_rc**2
  end subroutine init_tersoff
!=======================================================================
  subroutine force_tersoff(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr,d2lspr &
       ,mpi_world,myid,epi,epot,nismax,specorder,lstrs,iprint &
       ,tei)
    use vector,only: dot
    include 'mpif.h'
    include './params_unit.h'
!.....arguments
    integer,intent(in):: namax,natm,nnmax,nismax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_world,myid,lspr(0:nnmax,namax),nex(3)
    real(8),intent(in):: ra(3,namax),h(3,3),hi(3,3),sv(3,6) &
         ,tag(namax),rc,d2lspr(nnmax,namax)
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: lstrs
    character(len=3),intent(in):: specorder(nspmax)
    real(8),intent(in),optional:: tei(namax)

!.....local variables
    integer:: ia,ja,ka,jj,kk,ierr,ixyz,is,js,ks
    real(8):: xi(3),xij(3),rij(3),xik(3),rik(3),dij2,dij,diji &
         ,dik2,dik,diki,fc,dfc,fcij,dfcij,fcik,dfcik,tmp,dvdr &
         ,texp2ij,faij,zeta,texp3,zijk,gzijk,bij,dfaij,dvdij &
         ,dbijpref,sumj(3),dirij(3),djrij(3),dirik(3),dkrik(3) &
         ,dexp3ij,dexp3ik,cs,dgzijk,dics(3),djcs(3),dkcs(3) &
         ,diz(3),djz(3),dkz(3),dzdcs,texp,tmpk(3),pref
    real(8):: epotl,epotl1,epotl2,epott,t0
    real(8):: rcmax
    
    logical,save:: l1st = .true.
    real(8),allocatable,save:: aa1(:,:),aa2(:,:),strsl(:,:,:)

    t0 = mpi_wtime()

!.....only at 1st call
    if ( l1st ) then
      rcmax = 0d0
      do is=1,nsp
        rcmax = max(rcmax,ts_rc(is,is))
      enddo
      if( myid.eq.0 .and. iprint.ge.ipl_basic ) then
        print *,'Tersoff potential:'
        write(6,'(a,f7.3,a)') '   max cutoff radius = ',rcmax,' Ang.'
      endif
      if( rc.lt.rcmax ) then
        if( myid.eq.0 ) then
          print *,'ERROR: Cutoff radius too short !!!'
          write(6,'(a,2es12.4)') '  rc, rcmax = ',rc,rcmax
        endif
        call mpi_finalize(ierr)
        stop
      endif
!!$      if( ts_fc_type.eq.2 ) then
!!$        ts_rc = rc
!!$        ts_rc2= rc*rc
!!$        if( myid.eq.0 ) then
!!$          print *,'Use the given cutoff radius for this potential, ' &
!!$               //'since the cutoff function type is tanh.'
!!$        endif
!!$      endif
      if( ts_type(1:2).eq.'Te' .and. .not. present(tei) ) then
        if( myid.eq.0 ) print *,'ERROR: tei is not present.'
        stop 1
      endif
      
      allocate(aa1(3,namax),aa2(3,namax),strsl(3,3,namax))
      time = 0d0
      l1st = .false.
    endif

    if( size(aa1).ne.3*namax ) then
      deallocate(aa1,aa2,strsl)
      allocate(aa1(3,namax),aa2(3,namax),strsl(3,3,namax))
    endif

!.....initialize
    epotl = 0d0
    epi(1:natm+nb) = 0d0
    strsl(1:3,1:3,1:natm+nb) = 0d0

!.....Repulsive term which contains only 2-body
    epotl1 = 0d0
    aa1(:,1:natm+nb) = 0d0
    do ia=1,natm
      xi(1:3) = ra(1:3,ia)
      is = int(tag(ia))
      if( ts_type(1:2).eq.'Te' ) call set_ted_params(tei(ia),is)
      do jj=1,lspr(0,ia)
        ja = lspr(jj,ia)
        js = int(tag(ja))
        if( d2lspr(jj,ia).ge.ts_rc2(is,js) ) cycle
        if( ja.le.ia ) cycle
        xij(1:3) = ra(1:3,ja) -xi(1:3)
        rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij2 = rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
!!$        if( dij2.gt.ts_rc2 ) cycle
        dij = dsqrt(dij2)
!.....Potential
        fc = f_c(dij,is,js)
        dfc = df_c(dij,is,js)
        texp = exp(-ts_lmbd1(is,js)*dij)
        tmp = 0.5d0 *fc *ts_a(is,js) *texp
        epi(ia)= epi(ia) +tmp
        epotl1 = epotl1 +tmp
        if( ja.le.natm ) then
          epi(ja)= epi(ja) +tmp
          epotl1 = epotl1 +tmp
        endif
!.....Force
        diji = 1d0/dij
        dirij(1:3)= -rij(1:3)*diji
        djrij(1:3)= -dirij(1:3)
        dvdr= ts_a(is,js)*texp *(dfc -ts_lmbd1(is,js)*fc)
        aa1(1:3,ia)= aa1(1:3,ia) -dvdr*dirij(1:3)
        aa1(1:3,ja)= aa1(1:3,ja) -dvdr*djrij(1:3)
!.....Stress
        if( .not.lstrs ) cycle
        if( ja.le.natm) then
          do ixyz=1,3
            strsl(1:3,ixyz,ia)= strsl(1:3,ixyz,ia) &
                 -0.5d0 *rij(ixyz) *dvdr *djrij(1:3)
            strsl(1:3,ixyz,ja)= strsl(1:3,ixyz,ja) &
                 -0.5d0 *rij(ixyz) *dvdr *djrij(1:3)
          enddo
        else
          do ixyz=1,3
            strsl(1:3,ixyz,ia)= strsl(1:3,ixyz,ia) &
                 -0.5d0 *rij(ixyz) *dvdr *djrij(1:3)
          enddo
        endif
      enddo
      epi(ia) = epi(ia) +ts_f0(is)
      epotl1 = epotl1 +ts_f0(is)
    enddo
!!$    print *,'epotl1 =',epotl1

!.....Attractive term which contains 3-body as well
    epotl2 = 0d0
    aa2(:,1:natm+nb) = 0d0
    do ia=1,natm
      is = int(tag(ia))
      if( ts_type(1:2).eq.'Te' ) call set_ted_params(tei(ia),is)
      xi(1:3) = ra(1:3,ia)
      do jj=1,lspr(0,ia)
        ja = lspr(jj,ia)
        js = int(tag(ja))
        if( d2lspr(jj,ia).ge.ts_rc2(is,js) ) cycle
        if( ja.eq.ia ) cycle
        xij(1:3) = ra(1:3,ja) -xi(1:3)
        rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij2 = rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
!!$        if( dij2.gt.ts_rc2 ) cycle
        dij = dsqrt(dij2)
        diji = 1d0/dij
        dirij(1:3)= -rij(1:3)*diji
        djrij(1:3)= -dirij(1:3)
        fcij = f_c(dij,is,js)
        texp2ij = exp(-ts_lmbd2(is,js)*dij)
        faij = -ts_b(is,js) *texp2ij
!.....Precompute zeta which is necessary for derivative calculation
        zeta = 0d0
        do kk=1,lspr(0,ia)
          ka = lspr(kk,ia)
          ks = int(tag(ka))
          if( d2lspr(kk,ia).ge.ts_rc2(is,ks) ) cycle
          if( ka.eq.ja ) cycle
          xik(1:3) = ra(1:3,ka) -xi(1:3)
          rik(1:3) = h(1:3,1)*xik(1) +h(1:3,2)*xik(2) +h(1:3,3)*xik(3)
          dik2 = rik(1)*rik(1) +rik(2)*rik(2) +rik(3)*rik(3)
!!$          if( dik2.gt.ts_rc2 ) cycle
          dik = dsqrt(dik2)
          fcik = f_c(dik,is,ks)
          texp3 = exp(ts_alpha(is,js)*(dij-dik)**ts_beta(is))
          cs = dot(rij,rik) /dik *diji
          zijk = (ts_h(is) -cs)*(ts_h(is) -cs)
          gzijk = gz(zijk,is)
          zeta = zeta +fcik*gzijk*texp3
        enddo
!.....Potential
        bij = (1d0 +zeta**ts_eta(is))**(-ts_delta(is))
        tmp = 0.25d0 *fcij*bij*faij
        epi(ia)= epi(ia) +tmp
        epi(ja)= epi(ja) +tmp
        epotl2 = epotl2 +tmp +tmp
!.....Force except derivative of bij
        dfcij = df_c(dij,is,js)
        dfaij = ts_b(is,js) *ts_lmbd2(is,js) *texp2ij
        dvdij = 0.5d0 *bij*(dfcij*faij +dfaij*fcij)
!.....Since bij contains the angle around atom-i that is not symmetric to ia and ja,
!     the contribution to ja should be considered as well, even if the summations
!     of ia and ja are taken fully.
        aa2(1:3,ia)= aa2(1:3,ia) -dvdij*dirij(1:3)
        aa2(1:3,ja)= aa2(1:3,ja) -dvdij*djrij(1:3)
!.....Stress
        if( lstrs ) then
          do ixyz=1,3
            strsl(1:3,ixyz,ja)= strsl(1:3,ixyz,ja) &
                 -0.5d0 *rij(ixyz) *dvdij *djrij(1:3)
            strsl(1:3,ixyz,ia)= strsl(1:3,ixyz,ia) &
                 -0.5d0 *rij(ixyz) *dvdij *djrij(1:3)
          enddo
        endif
!.....Force w.r.t. derivative of bij
        dbijpref = (-ts_delta(is))*(1d0 +zeta**ts_eta(is))**(-(ts_delta(is)+1d0))  &
             *ts_eta(is) *zeta**(ts_eta(is)-1d0)
        pref = 0.5d0 *dbijpref *fcij*faij
        sumj(1:3) = 0d0
        do kk=1,lspr(0,ia)
          ka = lspr(kk,ia)
          ks = int(tag(ka))
          if( d2lspr(kk,ia).ge.ts_rc2(is,ks) ) cycle
          if( ka.eq.ja ) cycle
          xik(1:3) = ra(1:3,ka) -xi(1:3)
          rik(1:3) = h(1:3,1)*xik(1) +h(1:3,2)*xik(2) +h(1:3,3)*xik(3)
          dik2 = rik(1)*rik(1) +rik(2)*rik(2) +rik(3)*rik(3)
!!$          if( dik2.gt.ts_rc2 ) cycle
          dik = dsqrt(dik2)
          diki= 1d0/dik
          dirik(1:3)= -rik(1:3)*diki
          dkrik(1:3)= -dirik(1:3)
          fcik = f_c(dik,is,ks)
          dfcik = df_c(dik,is,ks)
          texp3 = exp(ts_alpha(is,js)*(dij-dik)**ts_beta(is))
          dexp3ij = ts_alpha(is,js) *(dij-dik)**(ts_beta(is) -1d0) *texp3
          dexp3ik = -dexp3ij
          cs = dot(rij,rik) *diji *diki
          zijk = (ts_h(is) -cs)*(ts_h(is) -cs)
          dzdcs = -2d0*(ts_h(is) -cs)
          gzijk = gz(zijk,is)
          dgzijk = dgz(zijk,is)
          djcs(1:3)= diji*diki*rik(1:3) -cs*djrij(1:3)*diji
          dkcs(1:3)= diji*diki*rij(1:3) -cs*dkrik(1:3)*diki
          dics(1:3)= -djcs(1:3) -dkcs(1:3)
          diz(1:3) = dzdcs*dics(1:3)
          djz(1:3) = dzdcs*djcs(1:3)
          dkz(1:3) = dzdcs*dkcs(1:3)
          sumj(1:3)= sumj(1:3) +fcik*(dgzijk*djz(1:3)*texp3 &
               +gzijk*dexp3ij*djrij(1:3))
          tmpk(1:3)= dfcik*gzijk*texp3*dkrik(1:3) &
               +dgzijk*fcik*texp3*dkz(1:3) +fcik*gzijk*dexp3ik*dkrik(1:3)
          aa2(1:3,ka)= aa2(1:3,ka) -tmpk(1:3) *pref
          aa2(1:3,ia)= aa2(1:3,ia) +tmpk(1:3) *pref
!.....Stress
          if( .not.lstrs ) cycle
          do ixyz=1,3
            strsl(1:3,ixyz,ka)= strsl(1:3,ixyz,ka) &
                 -0.5d0 *rik(ixyz) *pref *tmpk(1:3)
            strsl(1:3,ixyz,ia)= strsl(1:3,ixyz,ia) &
                 -0.5d0 *rik(ixyz) *pref *tmpk(1:3)
          enddo
        enddo
        aa2(1:3,ja) = aa2(1:3,ja) -pref*sumj(1:3)
        aa2(1:3,ia) = aa2(1:3,ia) +pref*sumj(1:3)
!.....Stress
        if( .not.lstrs ) cycle
        do ixyz=1,3
          strsl(1:3,ixyz,ja)= strsl(1:3,ixyz,ja) &
               -0.5d0 *rij(ixyz) *pref *sumj(1:3)
          strsl(1:3,ixyz,ia)= strsl(1:3,ixyz,ia) &
               -0.5d0 *rij(ixyz) *pref *sumj(1:3)
        enddo
      enddo  ! jj=...
    enddo  ! ia=...

!.....send back forces and potentials on immigrants
    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,aa2,3)
    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,epi,1)
    aa(1:3,1:natm)= aa(1:3,1:natm) +aa1(1:3,1:natm) +aa2(1:3,1:natm)

    if( lstrs ) then
      call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
           ,nn,mpi_world,strsl,9)
      strs(1:3,1:3,1:natm) = strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)
    endif

!.....gather epot
    epotl = epotl1 +epotl2
    call mpi_allreduce(epotl,epott,1,mpi_real8,mpi_sum,mpi_world,ierr)
    ts_epot = epott
    epot= epot +epott
    time = time + (mpi_wtime() -t0)

    if( iprint.ge.ipl_info .and. myid.eq.0 ) print *,'Tersoff epot = ',ts_epot
    return
  end subroutine force_tersoff
!=======================================================================
  function f_c(r,is,js)
    integer,intent(in):: is,js
    real(8),intent(in):: r
    real(8):: f_c

    if( ts_fc_type.eq.1 ) then
      f_c = f_c1(r,is,js)
      continue
    else if( ts_fc_type.eq.2 ) then
      f_c = f_c2(r,is,js)
      continue
    else ! default
      f_c = f_c3(r,is,js)
    endif
    return
  end function f_c
!=======================================================================
  function df_c(r,is,js)
    integer,intent(in):: is,js
    real(8),intent(in):: r
    real(8):: df_c

    df_c = 0d0
    if( ts_fc_type.eq.1 ) then
      df_c = df_c1(r,is,js)
      continue
    else if( ts_fc_type.eq.2 ) then
      df_c = df_c2(r,is,js)
      continue
    else ! default
      df_c = df_c3(r,is,js)
    endif
    return
  end function df_c
!=======================================================================
  function f_c1(r,is,js)
    integer,intent(in):: is,js
    real(8),intent(in):: r
    real(8):: f_c1
    include './params_unit.h'

    
    if( r.lt.ts_r(is,js) -ts_d(is,js) ) then
      f_c1 = 1d0
    else if( r.ge.ts_r(is,js) -ts_d(is,js) .and. &
         r.lt.ts_r(is,js) +ts_d(is,js) ) then
!!$      f_c1 =  0.5d0 -9d0/16d0 *sin(0.5d0*pi*(r-ts_r)/ts_d) &
!!$           -1d0/16d0*sin(1.5d0*pi*(r-ts_r)/ts_d)
      f_c1 =  0.5d0 *(1d0 -sin(0.5d0*pi*(r -ts_r(is,js))/ts_d(is,js)))
    else
      f_c1 = 0d0
    endif
    return
  end function f_c1
!=======================================================================
  function df_c1(r,is,js)
    integer,intent(in):: is,js
    real(8),intent(in):: r
    real(8):: df_c1
    include './params_unit.h'
    real(8):: p

    if( r.lt.ts_r(is,js) -ts_d(is,js) ) then
      df_c1 = 0d0
    else if( r.ge.ts_r(is,js) -ts_d(is,js) .and. &
         r.lt.ts_r(is,js) +ts_d(is,js) ) then
      p = 0.5d0 *pi /ts_d(is,js)
      df_c1 = -0.5d0 *p *cos(p*(r -ts_r(is,js)))
    else
      df_c1 = 0d0
    endif
    return
  end function df_c1
!=======================================================================
  function f_c2(r,is,js)
    integer,intent(in):: is,js
    real(8),intent(in):: r
    real(8):: f_c2
    include './params_unit.h'
    
    f_c2 = 0.5d0 *(1d0 -tanh(0.5d0*pi*(r -ts_r(is,js))/ts_d(is,js)))
    return
  end function f_c2
!=======================================================================
  function df_c2(r,is,js)
    integer,intent(in):: is,js
    real(8),intent(in):: r
    real(8):: df_c2
    include './params_unit.h'
    real(8):: p

    p = 0.5d0 *pi /ts_d(is,js)
    df_c2 = -0.5d0 *p /cosh(p*(r -ts_r(is,js)))**2
    return
  end function df_c2
!=======================================================================
  function f_c3(r,is,js)
    integer,intent(in):: is,js
    real(8),intent(in):: r
    real(8):: f_c3
    include './params_unit.h'
    
    if( r.lt.ts_r(is,js) -ts_d(is,js) ) then
      f_c3 = 1d0
    else if( r.ge.ts_r(is,js) -ts_d(is,js) .and. &
         r.lt.ts_r(is,js) +ts_d(is,js) ) then
      f_c3 =  0.5d0 -9d0/16d0 *sin(0.5d0*pi*(r-ts_r(is,js))/ts_d(is,js)) &
           -1d0/16d0*sin(1.5d0*pi*(r-ts_r(is,js))/ts_d(is,js))
    else
      f_c3 = 0d0
    endif
    return
  end function f_c3
!=======================================================================
  function df_c3(r,is,js)
    integer,intent(in):: is,js
    real(8),intent(in):: r
    real(8):: df_c3
    include './params_unit.h'
    real(8):: p

    if( r.lt.ts_r(is,js) -ts_d(is,js) ) then
      df_c3 = 0d0
    else if( r.ge.ts_r(is,js) -ts_d(is,js) .and. &
         r.lt.ts_r(is,js) +ts_d(is,js) ) then
      p = pi/2d0/ts_d(is,js)
      df_c3 = -9d0/16d0*p *cos(p*(r-ts_r(is,js))) &
           -1d0/16d0*3d0*p*cos(3d0*p*(r-ts_r(is,js)))
    else
      df_c3 = 0d0
    endif
    return
  end function df_c3
!=======================================================================
  function gz(z,is)
    integer,intent(in):: is
    real(8),intent(in):: z
    real(8):: gz
    real(8):: go,ga
    
    go = ts_c2(is)*z /(ts_c3(is) +z)
    ga = 1d0 +ts_c4(is)*exp(-ts_c5(is)*z)
    gz = ts_c1(is) +go*ga
    return
  end function gz
!=======================================================================
  function dgz(z,is)
    integer,intent(in):: is
    real(8),intent(in):: z
    real(8):: dgz
    real(8):: dgo,dga,go,ga,texp
    
    go = ts_c2(is)*z /(ts_c3(is) +z)
    dgo = ts_c2(is)*ts_c3(is)/(ts_c3(is) +z)**2
    texp = exp(-ts_c5(is)*z)
    ga = 1d0 +ts_c4(is)*texp
    dga = -ts_c4(is)*ts_c5(is)*texp

    dgz = dgo*ga +go*dga
    return
  end function dgz
!=======================================================================
  subroutine g_dg(z,g,dg,is)
    integer,intent(in):: is
    real(8),intent(in):: z
    real(8),intent(out):: g,dg

    real(8):: texp,go,dgo,ga,dga

!.....go term
    go = ts_c2(is)*z /(ts_c3(is) +z)
    dgo = ts_c2(is)*ts_c3(is)/(ts_c3(is) +z)**2
!.....ga term
    texp = exp(-ts_c5(is)*z)
    ga = 1d0 +ts_c4(is)*texp
    dga = -ts_c4(is)*ts_c5(is)*texp

    g = ts_c1(is) +go*ga
    dg = dgo*ga +go*dga
    return
  end subroutine g_dg
!=======================================================================
  subroutine read_params_tersoff(myid,mpi_world,iprint,specorder)
    use util, only: num_data
    implicit none
    include 'mpif.h'
    integer,intent(in):: myid,mpi_world,iprint
    character(len=3),intent(in):: specorder(nspmax)

    integer:: ite,ierr,isp,jsp,ksp,nd
    real(8):: te,a,b,lmbd1,lmbd2,eta,delta,alpha,beta,h,r1,r2 &
         ,c1,c2,c3,c4,c5,f0
    logical:: lexist
    character(len=128):: cfname,ctmp,cline
    character(len=3):: cspi,cspj,cspk

    if( myid.eq.0 ) then
!.....Check whether the file exists      
      cfname = trim(paramsdir)//'/'//trim(paramsfname)
      inquire(file=cfname,exist=lexist)
      if( .not. lexist ) then
        if( iprint.ge.ipl_basic ) then
          write(6,'(a)') ' WARNING: in.params.Tersoff does not exist !!!.'
          write(6,'(a)') '          Use default parameters.'
        endif
        call set_orig_params()
        goto 20
      endif

!.....Detect ts_type first to know the format of input file
      open(ioprms,file=cfname,status='old')
      do while(.true.)
        read(ioprms,'(a)',end=1) cline
        nd = num_data(cline,' ')
        if( nd.eq.0 ) cycle
        if( cline(1:1).eq.'!' .or. cline(1:1).eq.'#' ) cycle
        if( index(cline,'Tersoff_type').eq.0 ) cycle
        if( nd.lt.3 ) then
          print *,'ERROR: entry for Tersoff_type should have at least two arguments.'
          print *,'     e.g.)  Tersoff_type  Te-dependent  5'
          print *,'  The second argument is the number of electronic temperature points,' &
               //' which includes zero temperature.'
          stop 1
        endif
        backspace(ioprms)
        read(ioprms,*) ctmp, ts_type, ntemp
        allocate(ted_te(ntemp),ted_a(ntemp),ted_b(ntemp),ted_lmbd1(ntemp),ted_lmbd2(ntemp) &
             ,ted_eta(ntemp),ted_delta(ntemp),ted_alpha(ntemp),ted_beta(ntemp) &
             ,ted_h(ntemp),ted_r1(ntemp),ted_r2(ntemp),ted_c1(ntemp),ted_c2(ntemp) &
             ,ted_c3(ntemp),ted_c4(ntemp),ted_c5(ntemp),ted_f0(ntemp))
        exit
      enddo
1     close(ioprms)
      if( iprint.ge.ipl_basic ) print *,'ts_type = ',trim(ts_type)

!.....Read input file according to ts_type
      if( iprint.ge.ipl_basic ) write(6,'(/,a)') ' Tersoff parameters from file:'
      open(ioprms,file=cfname,status='old')

!.....Te-dependent) Te-dependent
!     Format of a single entry (one or more lines)
!       element, Te, A, B, lambda1, lambda2, eta, delta,
!       alpha, beta, h, R1, R2,
!       c1, c2, c3, c4, c5
      if( ts_type(1:2).eq.'Te' ) then
        ite = 0
        do while(.true.)
          read(ioprms,'(a)',end=10) cline
          nd = num_data(cline,' ')
          if( nd.eq.0 ) cycle
          if( cline(1:1).eq.'!' .or. cline(1:1).eq.'#' ) cycle
          if( index(cline,'Tersoff_type').ne.0 ) cycle

          ite = ite +1
          backspace(ioprms)
          read(ioprms,*) cspi,te,a,b,lmbd1,lmbd2,eta,delta,alpha,beta,h,r1,r2, &
               c1,c2,c3,c4,c5,f0
          if( iprint.ge.ipl_basic ) then
            print '(a5, f6.2, 2f11.5, 2f11.7)',trim(cspi),te,a,b,lmbd1,lmbd2
            print '(10x,f4.1,2f10.6,f4.1,f8.3,3f4.1)',eta,delta,alpha,beta,h,r1,r2
            print '(10x,f10.6,f11.2,f10.1,f6.3,f5.1,f6.3)',c1,c2,c3,c4,c5,f0
          endif

          if( iprint.ge.ipl_basic .and. trim(cspi).ne.'Si' ) then
            print *,'WARNING: Te-dependent is now available for only Si.'
          endif
          ted_te(ite) = te
          if( ite.eq.1 ) ted_te(ite) = 0d0  ! 1st temperature should be 0
          ted_a(ite) = a
          ted_b(ite) = b
          ted_lmbd1(ite) = lmbd1
          ted_lmbd2(ite) = lmbd2
          ted_eta(ite) = eta
          ted_delta(ite) = delta
          ted_alpha(ite) = alpha
          ted_beta(ite) = beta
          ted_h(ite) = h
          ted_r1(ite) = r1
          ted_r2(ite) = r2
          ted_c1(ite) = c1
          ted_c2(ite) = c2
          ted_c3(ite) = c3
          ted_c4(ite) = c4
          ted_c5(ite) = c5
          ted_f0(ite) = f0
          if( ite.eq.ntemp ) exit
        enddo
        if( ite.ne.ntemp ) then
          print *,'ERROR: ite.ne.ntemp !'
          stop 1
        endif

!.....default) Kumagai's modified Tersoff potential
      else
        do while(.true.)
          read(ioprms,'(a)',end=10) cline
          nd = num_data(cline,' ')
          if( nd.eq.0 ) cycle
          if( cline(1:1).eq.'!' .or. cline(1:1).eq.'#' ) cycle
          if( index(cline,'Tersoff_type').ne.0 ) cycle

          backspace(ioprms)
          read(ioprms,*) cspi, a,b,lmbd1,lmbd2,eta,delta &
               ,alpha,beta,h,r1,r2,c1,c2,c3,c4,c5
          isp = csp2isp(cspi)
          ts_a(isp,isp) = a
          ts_b(isp,isp) = b
          ts_lmbd1(isp,isp) = lmbd1
          ts_lmbd2(isp,isp) = lmbd2
          ts_eta(isp) = eta
          ts_delta(isp) = delta
          ts_alpha(isp,isp) = alpha
          ts_beta(isp) = beta
          ts_h(isp) = h
          ts_r1(isp,isp) = r1
          ts_r2(isp,isp) = r2
          ts_c1(isp) = c1
          ts_c2(isp) = c2
          ts_c3(isp) = c3
          ts_c4(isp) = c4
          ts_c5(isp) = c5

          if( iprint.ge.ipl_basic ) then
            write(6,'(a5,4es12.4)') trim(cspi), a, b, lmbd1,lmbd2
            write(6,'(15x,7es12.4)') eta,delta,alpha,beta,h,r1,r2
            write(6,'(15x,5es12.4)') c1,c2,c3,c4,c5
          endif
          goto 10
        enddo
      endif
10    close(ioprms)

    endif  !  myid.eq.0

20  continue

    call mpi_bcast(ts_type,12,mpi_character,0,mpi_world,ierr)

    if( ts_type(1:2).eq.'Te' ) then
      call mpi_bcast(ntemp,1,mpi_integer,0,mpi_world,ierr)
      if( .not. allocated(ted_te) ) allocate(ted_te(ntemp),ted_a(ntemp) &
           ,ted_b(ntemp),ted_lmbd1(ntemp),ted_lmbd2(ntemp) &
           ,ted_eta(ntemp),ted_delta(ntemp),ted_alpha(ntemp),ted_beta(ntemp) &
           ,ted_h(ntemp),ted_r1(ntemp),ted_r2(ntemp),ted_c1(ntemp),ted_c2(ntemp) &
           ,ted_c3(ntemp),ted_c4(ntemp),ted_c5(ntemp),ted_f0(ntemp))
      call mpi_bcast(ted_te,ntemp,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ted_a,ntemp,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ted_b,ntemp,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ted_lmbd1,ntemp,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ted_lmbd2,ntemp,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ted_eta,ntemp,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ted_delta,ntemp,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ted_alpha,ntemp,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ted_beta,ntemp,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ted_h,ntemp,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ted_r1,ntemp,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ted_r2,ntemp,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ted_c1,ntemp,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ted_c2,ntemp,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ted_c3,ntemp,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ted_c4,ntemp,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ted_c5,ntemp,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ted_f0,ntemp,mpi_real8,0,mpi_world,ierr)
      ts_a(:,:) = ted_a(1)
      ts_b(:,:) = ted_b(1)
      ts_lmbd1(:,:) = ted_lmbd1(1)
      ts_lmbd2(:,:) = ted_lmbd2(1)
      ts_eta(:) = ted_eta(1)
      ts_delta(:) = ted_delta(1)
      ts_alpha(:,:) = ted_alpha(1)
      ts_beta(:) = ted_beta(1)
      ts_h(:) = ted_h(1)
      ts_r1(:,:) = ted_r1(1)
      ts_r2(:,:) = ted_r2(1)
      ts_c1(:) = ted_c1(1)
      ts_c2(:) = ted_c2(1)
      ts_c3(:) = ted_c3(1)
      ts_c4(:) = ted_c4(1)
      ts_c5(:) = ted_c5(1)
      ts_f0(:) = ted_f0(1)
    else
      call mpi_bcast(ts_a,nspmax*nspmax,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ts_b,nspmax*nspmax,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ts_lmbd1,nspmax*nspmax,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ts_lmbd2,nspmax*nspmax,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ts_eta,nspmax,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ts_delta,nspmax,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ts_alpha,nspmax*nspmax,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ts_beta,nspmax,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ts_h,nspmax,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ts_r1,nspmax*nspmax,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ts_r2,nspmax*nspmax,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ts_c1,nspmax,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ts_c2,nspmax,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ts_c3,nspmax,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ts_c4,nspmax,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ts_c5,nspmax,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ts_f0,nspmax,mpi_real8,0,mpi_world,ierr)
    endif
    return
  end subroutine read_params_tersoff
!=======================================================================
  subroutine set_ted_params(te,is)
!
!  Set Te dependent parameters by linear interpolation using given Te
!
    include './params_unit.h'
    integer,intent(in):: is
    real(8),intent(in):: te

    integer:: ite,ite0,ite1
    real(8):: x,tev

!.....Since the ted_te is in eV, unit conversion is required
    tev = te*k2ev
    do ite=ntemp,1,-1
      if( tev.ge.ted_te(ite) ) then
        if( ite.eq.ntemp ) then
          ite0 = ite-1
          ite1 = ite
          x = 1d0
        else
          ite0 = ite
          ite1 = ite +1
          x = (tev -ted_te(ite))/(ted_te(ite1) -ted_te(ite))
        endif
        exit
      endif
    enddo
    
!.....Interpolation
    ts_a(is,is) = ted_a(ite0) +(ted_a(ite1) -ted_a(ite0)) *x
    ts_b(is,is) = ted_b(ite0) +(ted_b(ite1) -ted_b(ite0)) *x
    ts_lmbd1(is,is) = ted_lmbd1(ite0) +(ted_lmbd1(ite1) -ted_lmbd1(ite0)) *x
    ts_lmbd2(is,is) = ted_lmbd2(ite0) +(ted_lmbd2(ite1) -ted_lmbd2(ite0)) *x
    ts_eta(is) = ted_eta(ite0) +(ted_eta(ite1) -ted_eta(ite0)) *x
    ts_delta(is) = ted_delta(ite0) +(ted_delta(ite1) -ted_delta(ite0)) *x
    ts_alpha(is,is) = ted_alpha(ite0) +(ted_alpha(ite1) -ted_alpha(ite0)) *x
    ts_beta(is) = ted_beta(ite0) +(ted_beta(ite1) -ted_beta(ite0)) *x
    ts_h(is) = ted_h(ite0) +(ted_h(ite1) -ted_h(ite0)) *x
    ts_r1(is,is) = ted_r1(ite0) +(ted_r1(ite1) -ted_r1(ite0)) *x
    ts_r2(is,is) = ted_r2(ite0) +(ted_r2(ite1) -ted_r2(ite0)) *x
    ts_c1(is) = ted_c1(ite0) +(ted_c1(ite1) -ted_c1(ite0)) *x
    ts_c2(is) = ted_c2(ite0) +(ted_c2(ite1) -ted_c2(ite0)) *x
    ts_c3(is) = ted_c3(ite0) +(ted_c3(ite1) -ted_c3(ite0)) *x
    ts_c4(is) = ted_c4(ite0) +(ted_c4(ite1) -ted_c4(ite0)) *x
    ts_c5(is) = ted_c5(ite0) +(ted_c5(ite1) -ted_c5(ite0)) *x
    ts_f0(is) = ted_f0(ite0) +(ted_f0(ite1) -ted_f0(ite0)) *x

    ts_r(is,is) = (ts_r1(is,is) +ts_r2(is,is))/2
    ts_d(is,is) = (ts_r2(is,is) -ts_r1(is,is))/2
    ts_rc(is,is) = ts_r(is,is) +ts_d(is,is)
    ts_rc2(is,is)= ts_rc(is,is)**2
    return
  end subroutine set_ted_params
!=======================================================================
  subroutine set_orig_params()
!
!  Set original parameter set taken from Kumagai et al.[2], if not given.
!

    ts_type = 'default'
!.....Cutoff type: (1) Tersoff,  (2) Kumagai, (3) Murty (default)
    ts_fc_type = 3
    
!.....Original parameters from [2]
    ts_a(:,:) = 3281.5905d0
    ts_b(:,:) = 121.00047d0
    ts_lmbd1(:,:) = 3.2300135d0
    ts_lmbd2(:,:) = 1.345797d0
    ts_eta(:) = 1.0d0
    ts_delta(:) = 0.53298909d0
    ts_alpha(:,:) = 2.3890327d0
    ts_beta(:) = 1d0
!!$    ts_beta_ters = 1d0  ! beta of original Tersoff [1], must be 1.0 for [2]
    ts_c1(:) = 0.20173476d0
    ts_c2(:) = 730418.72d0
    ts_c3(:) = 1000000d0
    ts_c4(:) = 1.0d0
    ts_c5(:) = 26.0d0
    ts_h(:) = -0.365
    ts_r1(:,:) = 2.7d0
    ts_r2(:,:) = 3.3d0
    ts_f0(:) = 0d0
    return
  end subroutine set_orig_params
end module tersoff
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd lib"
!     End:
