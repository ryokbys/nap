module tersoff
!-----------------------------------------------------------------------
!                     Last modified: <2023-03-30 14:31:10 KOBAYASHI Ryo>
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

!.....Cutoff type: (1) Tersoff, (3) Murty (default)
  integer:: ts_fc_type = 3

!.....Original parameters from [2]
  real(8):: ts_a(nspmax,nspmax), ts_b(nspmax,nspmax), ts_lmbd1(nspmax,nspmax), &
       ts_lmbd2(nspmax,nspmax), ts_eta(nspmax,nspmax), ts_delta(nspmax,nspmax), &
       ts_alpha(nspmax,nspmax), ts_beta(nspmax,nspmax), ts_c1(nspmax,nspmax), &
       ts_c2(nspmax,nspmax), ts_c3(nspmax,nspmax), ts_c4(nspmax,nspmax),&
       ts_c5(nspmax,nspmax), ts_h(nspmax,nspmax), &
       ts_rc2in(nspmax,nspmax), ts_rc2out(nspmax,nspmax), ts_f0(nspmax,nspmax), &
       ts_rc3in(nspmax,nspmax), ts_rc3out(nspmax,nspmax), &
       ts_r(nspmax,nspmax), ts_d(nspmax,nspmax), &
       ts_rc(nspmax,nspmax), ts_rc2(nspmax,nspmax)

!$omp threadprivate(ts_a,ts_b,ts_lmbd1,ts_lmbd2,ts_eta,ts_delta,ts_alpha,ts_beta,ts_c1, &
!$omp               ts_c2,ts_c3,ts_c4,ts_c5,ts_h,ts_rc2in,ts_rc2out,ts_f0,&
!$omp               ts_rc3in,ts_rc3out,ts_r,ts_d,ts_rc,ts_rc2 )

  logical:: interact(nspmax,nspmax)

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

    ts_rc(:,:) = ts_rc2out(:,:)
    ts_rc2(:,:)= ts_rc(:,:)**2

!!$!.....Lorentz-Berthelot rule for inter-species parameters
!!$    do isp=1,nspmax
!!$      do jsp=1,nspmax
!!$        ts_a(isp,jsp) = sqrt(ts_a(isp,isp)*ts_a(jsp,jsp))
!!$        ts_b(isp,jsp) = sqrt(ts_b(isp,isp)*ts_b(jsp,jsp))
!!$        ts_lmbd1(isp,jsp)= (ts_lmbd1(isp,isp)+ts_lmbd1(jsp,jsp))/2
!!$        ts_lmbd2(isp,jsp)= (ts_lmbd2(isp,isp)+ts_lmbd2(jsp,jsp))/2
!!$        ts_rc2in(isp,jsp)= (ts_rc2in(isp,isp)+ts_rc2in(jsp,jsp))/2
!!$        ts_rc2out(isp,jsp)= (ts_rc2out(isp,isp)+ts_rc2out(jsp,jsp))/2
!!$        ts_rc3in(isp,jsp)= (ts_rc3in(isp,isp)+ts_rc3in(jsp,jsp))/2
!!$        ts_rc3out(isp,jsp)= (ts_rc3out(isp,isp)+ts_rc3out(jsp,jsp))/2
!!$        ts_rc(isp,jsp)= ts_rc2out(isp,jsp)
!!$        ts_rc2(isp,jsp)= ts_rc(isp,jsp)**2
!!$      enddo
!!$    enddo
!!$    ts_r = (ts_r1 +ts_r2)/2
!!$    ts_d = (ts_r2 -ts_r1)/2
!!$    ts_rc = ts_r +ts_d
!!$    ts_rc2= ts_rc**2
  end subroutine init_tersoff
!=======================================================================
  subroutine force_tersoff(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
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
         ,tag(namax),rc
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: lstrs
    character(len=3),intent(in):: specorder(nspmax)
    real(8),intent(in),optional:: tei(namax)

!.....local variables
    integer:: ia,ja,ka,jj,kk,ierr,ixyz,jxyz,is,js,ks
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
        do js=is,nsp
          rcmax = max(rcmax,ts_rc(is,js))
        enddo
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
      if( ts_type(1:2).eq.'Te' .and. .not. present(tei) ) then
        if( myid.eq.0 ) print *,'ERROR: tei is not present.'
        stop 1
      endif
      
      allocate(aa1(3,namax),aa2(3,namax),strsl(3,3,namax))
      time = 0d0
      l1st = .false.
    endif

    if( size(aa1).lt.3*namax ) then
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
!$omp parallel copyin(ts_a,ts_b,ts_lmbd1,ts_lmbd2,ts_eta,ts_delta,ts_alpha,ts_beta,ts_c1, &
!$omp                 ts_c2,ts_c3,ts_c4,ts_c5,ts_h,ts_rc2in,ts_rc2out,ts_f0,&
!$omp                 ts_rc3in,ts_rc3out,ts_r,ts_d,ts_rc,ts_rc2)

!$omp do private(ia,xi,is,jj,ja,js,xij,rij,dij2,dij,fc,dfc,texp,tmp, &
!$omp            diji,dirij,djrij,dvdr,ixyz,jxyz ) &
!$omp    reduction(+:epotl1)
    do ia=1,natm
      xi(1:3) = ra(1:3,ia)
      is = int(tag(ia))
      if( ts_type(1:2).eq.'Te' ) call set_ted_params(tei(ia),is)
      do jj=1,lspr(0,ia)
        ja = lspr(jj,ia)
        js = int(tag(ja))
        if( .not.interact(is,js) ) cycle
        if( ja.le.ia ) cycle
        xij(1:3) = ra(1:3,ja) -xi(1:3)
        rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij2 = rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
        if( dij2.ge.ts_rc2(is,js) ) cycle
        dij = dsqrt(dij2)
!.....Potential
        fc = f_c(dij,ts_rc2in(is,js),ts_rc2out(is,js))
        dfc = df_c(dij,ts_rc2in(is,js),ts_rc2out(is,js))
        texp = exp(-ts_lmbd1(is,js)*dij)
        tmp = 0.5d0 *fc *ts_a(is,js) *texp
!omp atomic
        epi(ia)= epi(ia) +tmp
        epotl1 = epotl1 +tmp
        if( ja.le.natm ) then
!omp atomic
          epi(ja)= epi(ja) +tmp
          epotl1 = epotl1 +tmp
        endif
!.....Force
        diji = 1d0/dij
        dirij(1:3)= -rij(1:3)*diji
        djrij(1:3)= -dirij(1:3)
        dvdr= ts_a(is,js)*texp *(dfc -ts_lmbd1(is,js)*fc)
        do ixyz=1,3
!$omp atomic
          aa1(ixyz,ia)= aa1(ixyz,ia) -dvdr*dirij(ixyz)
!$omp atomic
          aa1(ixyz,ja)= aa1(ixyz,ja) -dvdr*djrij(ixyz)
        enddo
!.....Stress
        if( ja.le.natm) then
          do ixyz=1,3
            do jxyz=1,3
!$omp atomic
              strsl(jxyz,ixyz,ia)= strsl(jxyz,ixyz,ia) &
                   -0.5d0 *rij(ixyz) *dvdr *djrij(jxyz)
!$omp atomic
              strsl(jxyz,ixyz,ja)= strsl(jxyz,ixyz,ja) &
                   -0.5d0 *rij(ixyz) *dvdr *djrij(jxyz)
            enddo
          enddo
        else
          do ixyz=1,3
            do jxyz=1,3
!$omp atomic
              strsl(jxyz,ixyz,ia)= strsl(jxyz,ixyz,ia) &
                   -0.5d0 *rij(ixyz) *dvdr *djrij(jxyz)
            enddo
          enddo
        endif
      enddo
!$omp atomic
      epi(ia) = epi(ia) +ts_f0(is,is)
      epotl1 = epotl1 +ts_f0(is,is)
    enddo
!$omp end do

!.....Attractive term which contains 3-body as well
!$omp single
    epotl2 = 0d0
    aa2(:,1:natm+nb) = 0d0
!$omp end single
    
!$omp do private(ia,is,xi,jj,ja,js,xij,rij,dij2,dij,diji,dirij,djrij,fcij, &
!$omp            texp2ij,faij,zeta,kk,ka,ks,xik,rik,dik2,dik,fcik,texp3,ixyz,jxyz, &
!$omp            cs,zijk,gzijk,bij,tmp,dfcij,dfaij,dvdij,dbijpref,pref,sumj, &
!$omp            diki,dirik,dkrik,dfcik, &
!$omp            dexp3ij,dexp3ik,dzdcs,dgzijk,djcs,dkcs,dics,diz,djz,dkz,tmpk ) &
!$omp    reduction(+:epotl2)
    do ia=1,natm
      is = int(tag(ia))
      if( ts_type(1:2).eq.'Te' ) call set_ted_params(tei(ia),is)
      xi(1:3) = ra(1:3,ia)
      do jj=1,lspr(0,ia)
        ja = lspr(jj,ia)
        js = int(tag(ja))
        if( .not.interact(is,js) ) cycle
        if( ja.eq.ia ) cycle
        xij(1:3) = ra(1:3,ja) -xi(1:3)
        rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij2 = rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
        if( dij2.gt.ts_rc2(is,js) ) cycle
        dij = dsqrt(dij2)
        diji = 1d0/dij
        dirij(1:3)= -rij(1:3)*diji
        djrij(1:3)= -dirij(1:3)
        fcij = f_c(dij,ts_rc2in(is,js),ts_rc2out(is,js))
        texp2ij = exp(-ts_lmbd2(is,js)*dij)
        faij = -ts_b(is,js) *texp2ij
!.....Precompute zeta which is necessary for derivative calculation
        zeta = 0d0
        do kk=1,lspr(0,ia)
          ka = lspr(kk,ia)
          ks = int(tag(ka))
!.....Currently Tersoff potential is only available fro js==ks
          if( ks.ne.js ) cycle
          if( ka.eq.ja ) cycle
          xik(1:3) = ra(1:3,ka) -xi(1:3)
          rik(1:3) = h(1:3,1)*xik(1) +h(1:3,2)*xik(2) +h(1:3,3)*xik(3)
          dik2 = rik(1)*rik(1) +rik(2)*rik(2) +rik(3)*rik(3)
          if( dik2.gt.ts_rc2(is,ks) ) cycle
          dik = dsqrt(dik2)
          fcik = f_c(dik,ts_rc3in(is,ks),ts_rc3out(is,ks))
          texp3 = exp(ts_alpha(is,ks)*(dij-dik)**ts_beta(is,ks))
          cs = dot(rij,rik) /dik *diji
          zijk = (ts_h(is,ks) -cs)*(ts_h(is,ks) -cs)
          gzijk = gz(zijk,is,js,ks)
          zeta = zeta +fcik*gzijk*texp3
        enddo
!.....Potential
        bij = (1d0 +zeta**ts_eta(is,js))**(-ts_delta(is,js))
        tmp = 0.25d0 *fcij*bij*faij
!$omp atomic
        epi(ia)= epi(ia) +tmp
!$omp atomic
        epi(ja)= epi(ja) +tmp
        epotl2 = epotl2 +tmp +tmp
!.....Force except derivative of bij
        dfcij = df_c(dij,ts_rc2in(is,js),ts_rc3in(is,js))
        dfaij = ts_b(is,js) *ts_lmbd2(is,js) *texp2ij
        dvdij = 0.5d0 *bij*(dfcij*faij +dfaij*fcij)
!.....Since bij contains the angle around atom-i that is not symmetric to ia and ja,
!     the contribution to ja should be considered as well, even if the summations
!     of ia and ja are taken fully.
        do ixyz=1,3
!$omp atomic
          aa2(ixyz,ia)= aa2(ixyz,ia) -dvdij*dirij(ixyz)
!$omp atomic
          aa2(ixyz,ja)= aa2(ixyz,ja) -dvdij*djrij(ixyz)
        enddo
!.....Stress
        do ixyz=1,3
          do jxyz=1,3
!$omp atomic
            strsl(jxyz,ixyz,ja)= strsl(jxyz,ixyz,ja) &
                 -0.5d0 *rij(ixyz) *dvdij *djrij(jxyz)
!$omp atomic
            strsl(jxyz,ixyz,ia)= strsl(jxyz,ixyz,ia) &
                 -0.5d0 *rij(ixyz) *dvdij *djrij(jxyz)
          enddo
        enddo
!.....Force w.r.t. derivative of bij
        dbijpref = (-ts_delta(is,js))*(1d0 +zeta**ts_eta(is,js))**(-(ts_delta(is,js)+1d0))  &
             *ts_eta(is,js) *zeta**(ts_eta(is,js)-1d0)
        pref = 0.5d0 *dbijpref *fcij*faij
        sumj(1:3) = 0d0
        do kk=1,lspr(0,ia)
          ka = lspr(kk,ia)
          if( ka.eq.ja ) cycle
          ks = int(tag(ka))
!.....Currently Tersoff potential is only available fro js==ks
          if( ks.ne.js ) cycle
          xik(1:3) = ra(1:3,ka) -xi(1:3)
          rik(1:3) = h(1:3,1)*xik(1) +h(1:3,2)*xik(2) +h(1:3,3)*xik(3)
          dik2 = rik(1)*rik(1) +rik(2)*rik(2) +rik(3)*rik(3)
          if( dik2.ge.ts_rc2(is,ks) ) cycle
          dik = dsqrt(dik2)
          diki= 1d0/dik
          dirik(1:3)= -rik(1:3)*diki
          dkrik(1:3)= -dirik(1:3)
          fcik = f_c(dik,ts_rc3in(is,ks),ts_rc3out(is,ks))
          dfcik = df_c(dik,ts_rc3in(is,ks),ts_rc3out(is,ks))
          texp3 = exp(ts_alpha(is,js)*(dij-dik)**ts_beta(is,js))
          dexp3ij = ts_alpha(is,js) *(dij-dik)**(ts_beta(is,js) -1d0) *texp3
          dexp3ik = -dexp3ij
          cs = dot(rij,rik) *diji *diki
          zijk = (ts_h(is,ks) -cs)*(ts_h(is,ks) -cs)
          dzdcs = -2d0*(ts_h(is,ks) -cs)
          gzijk = gz(zijk,is,js,ks)
          dgzijk = dgz(zijk,is,js,ks)
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
          do ixyz=1,3
!$omp atomic
            aa2(ixyz,ka)= aa2(ixyz,ka) -tmpk(ixyz) *pref
!$omp atomic
            aa2(ixyz,ia)= aa2(ixyz,ia) +tmpk(ixyz) *pref
          enddo
!.....Stress
          do ixyz=1,3
            do jxyz=1,3
!$omp atomic
              strsl(jxyz,ixyz,ka)= strsl(jxyz,ixyz,ka) &
                   -0.5d0 *rik(ixyz) *pref *tmpk(jxyz)
!$omp atomic
              strsl(jxyz,ixyz,ia)= strsl(jxyz,ixyz,ia) &
                   -0.5d0 *rik(ixyz) *pref *tmpk(jxyz)
            enddo
          enddo
        enddo
        do ixyz=1,3
!$omp atomic
          aa2(ixyz,ja) = aa2(ixyz,ja) -pref*sumj(ixyz)
!$omp atomic
          aa2(ixyz,ia) = aa2(ixyz,ia) +pref*sumj(ixyz)
        enddo
!.....Stress
        do ixyz=1,3
          do jxyz=1,3
!$omp atomic
            strsl(jxyz,ixyz,ja)= strsl(jxyz,ixyz,ja) &
                 -0.5d0 *rij(ixyz) *pref *sumj(jxyz)
!$omp atomic
            strsl(jxyz,ixyz,ia)= strsl(jxyz,ixyz,ia) &
                 -0.5d0 *rij(ixyz) *pref *sumj(jxyz)
          enddo
        enddo
      enddo  ! jj=...
    enddo  ! ia=...
!$omp end do    
!$omp end parallel

!.....send back forces and potentials on immigrants
    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,aa2,3)
    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,epi,1)
    aa(1:3,1:natm)= aa(1:3,1:natm) +aa1(1:3,1:natm) +aa2(1:3,1:natm)
    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,strsl,9)
    strs(1:3,1:3,1:natm) = strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)

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
  function f_c(r,rcin,rcout)
    real(8),intent(in):: r,rcin,rcout
    real(8):: f_c

    if( ts_fc_type.eq.1 ) then
      f_c = f_c1(r,rcin,rcout)
      continue
!!$    else if( ts_fc_type.eq.2 ) then
!!$      f_c = f_c2(r,rcin,rcout)
!!$      continue
    else ! default
      f_c = f_c3(r,rcin,rcout)
    endif
    return
  end function f_c
!=======================================================================
  function df_c(r,rcin,rcout)
    real(8),intent(in):: r,rcin,rcout
    real(8):: df_c

    df_c = 0d0
    if( ts_fc_type.eq.1 ) then
      df_c = df_c1(r,rcin,rcout)
!!$      continue
!!$    else if( ts_fc_type.eq.2 ) then
!!$      df_c = df_c2(r,rcin,rcout)
!!$      continue
    else ! default
      df_c = df_c3(r,rcin,rcout)
    endif
    return
  end function df_c
!=======================================================================
  function f_c1(r,rcin,rcout)
    real(8),intent(in):: r,rcin,rcout
    real(8):: f_c1
    include './params_unit.h'

    
    if( r.lt.rcin ) then
      f_c1 = 1d0
    else if( r.ge.rcin .and. r.lt.rcout ) then
      f_c1 =  0.5d0 *(1d0 +cos(pi*(r -rcin)/(rcout-rcin)))
    else
      f_c1 = 0d0
    endif
    return
  end function f_c1
!=======================================================================
  function df_c1(r,rcin,rcout)
    real(8),intent(in):: r,rcin,rcout
    real(8):: df_c1
    include './params_unit.h'
    real(8):: p

    if( r.lt.rcin ) then
      df_c1 = 0d0
    else if( r.ge.rcin .and. r.lt.rcout ) then
      p = pi /(rcout-rcin)
      df_c1 = -0.5d0 *p *sin(p*(r -rcin))
    else
      df_c1 = 0d0
    endif
    return
  end function df_c1
!=======================================================================
  function f_c3(r,rcin,rcout)
    real(8),intent(in):: r,rcin,rcout
    real(8):: f_c3
    include './params_unit.h'
    real(8):: p
    
    if( r.lt.rcin ) then
      f_c3 = 1d0
    else if( r.ge.rcin .and. r.lt.rcout ) then
      p = pi *(r-rcin)/(rcout-rcin)
      f_c3 =  0.5d0 +9d0/16d0 *cos(p) -1d0/16d0*cos(3d0*p)
    else
      f_c3 = 0d0
    endif
    return
  end function f_c3
!=======================================================================
  function df_c3(r,rcin,rcout)
    real(8),intent(in):: r,rcin,rcout
    real(8):: df_c3
    include './params_unit.h'
    real(8):: p

    if( r.lt.rcin ) then
      df_c3 = 0d0
    else if( r.ge.rcin .and. r.lt.rcout ) then
      p = pi/(rcout-rcin)
      df_c3 = -9d0/16d0*p *sin(p*(r-rcin)) +1d0/16d0*3d0*p*sin(3d0*p*(r-rcin))
    else
      df_c3 = 0d0
    endif
    return
  end function df_c3
!=======================================================================
  function gz(z,is,js,ks)
    integer,intent(in):: is,js,ks
    real(8),intent(in):: z
    real(8):: gz
    real(8):: go,ga
    
    go = ts_c2(is,js)*z /(ts_c3(is,js) +z)
    ga = 1d0 +ts_c4(is,js)*exp(-ts_c5(is,js)*z)
    gz = ts_c1(is,js) +go*ga
    return
  end function gz
!=======================================================================
  function dgz(z,is,js,ks)
    integer,intent(in):: is,js,ks
    real(8),intent(in):: z
    real(8):: dgz
    real(8):: dgo,dga,go,ga,texp
    
    go = ts_c2(is,js)*z /(ts_c3(is,js) +z)
    dgo = ts_c2(is,js)*ts_c3(is,js)/(ts_c3(is,js) +z)**2
    texp = exp(-ts_c5(is,js)*z)
    ga = 1d0 +ts_c4(is,js)*texp
    dga = -ts_c4(is,js)*ts_c5(is,js)*texp

    dgz = dgo*ga +go*dga
    return
  end function dgz
!=======================================================================
  subroutine g_dg(z,g,dg,is,js,ks)
    integer,intent(in):: is,js,ks
    real(8),intent(in):: z
    real(8),intent(out):: g,dg

    real(8):: texp,go,dgo,ga,dga

!.....go term
    go = ts_c2(is,js)*z /(ts_c3(is,js) +z)
    dgo = ts_c2(is,js)*ts_c3(is,js)/(ts_c3(is,js) +z)**2
!.....ga term
    texp = exp(-ts_c5(is,js)*z)
    ga = 1d0 +ts_c4(is,js)*texp
    dga = -ts_c4(is,js)*ts_c5(is,js)*texp

    g = ts_c1(is,js) +go*ga
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
    real(8):: te,a,b,lmbd1,lmbd2,eta,delta,alpha,beta,h, &
         c1,c2,c3,c4,c5,f0, rc2in,rc2out,rc3in,rc3out
    logical:: lexist
    character(len=128):: cfname,ctmp,cline
    character(len=3):: cspi,cspj,cspk

    if( myid.eq.0 ) then
      interact(:,:) = .false.
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

!.....Read input file according to ts_type
      if( iprint.ge.ipl_basic ) then
        write(6,'(/,a)') ' Tersoff parameters from file:'
        print *,'  ts_type = ',trim(ts_type)
      endif
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
          read(ioprms,*) cspi,te,a,b,lmbd1,lmbd2,eta,delta,alpha,beta,h,rc2in,rc2out, &
               c1,c2,c3,c4,c5,f0

          isp = csp2isp(cspi)
          if( isp.le.0 ) then
            print *,' Tersoff parameter read but not used for ',trim(cspi)
            cycle
          endif
          if( iprint.ge.ipl_basic ) then
            print '(a5, f6.2, 2f11.5, 2f11.7)',trim(cspi),te,a,b,lmbd1,lmbd2
            print '(10x,f4.1,2f10.6,f4.1,f8.3,3f4.1)',eta,delta,alpha,beta,h,rc2in,rc2out
            print '(10x,f10.6,f11.2,f10.1,f6.3,f5.1,f6.3)',c1,c2,c3,c4,c5,f0
          endif

          if( iprint.ge.ipl_basic .and. trim(cspi).ne.'Si' ) then
            print *,'WARNING: Te-dependent is now available for only Si.'
            cycle
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
          ted_r1(ite) = rc2in
          ted_r2(ite) = rc2out
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
!     Format of a single entry (one or more lines)
!       element1, element2, A, B, lambda1, lambda2, eta, delta,
!       alpha, beta, h, R1, R2,
!       c1, c2, c3, c4, c5
      else
        do while(.true.)
          read(ioprms,'(a)',end=10) cline
          nd = num_data(cline,' ')
          if( nd.eq.0 ) cycle
          if( cline(1:1).eq.'!' .or. cline(1:1).eq.'#' ) cycle
          if( index(cline,'Tersoff_type').ne.0 ) cycle

          backspace(ioprms)
!.....TODO: make it applicable to multi-species...
          read(ioprms,*) cspi, cspj, a,b,lmbd1,lmbd2, rc2in, rc2out, &
               eta,delta,alpha,beta,h,c1,c2,c3,c4,c5,rc3in,rc3out
          isp = csp2isp(cspi)
          jsp = csp2isp(cspj)
          if( isp.le.0 .or. jsp.le.0 ) then
            print '(a,2a4)',' Tersoff parameter read but no such species in the system: ', &
                 trim(cspi),trim(cspj)
            cycle
          endif

          ts_a(isp,jsp) = a
          ts_b(isp,jsp) = b
          ts_lmbd1(isp,jsp) = lmbd1
          ts_lmbd2(isp,jsp) = lmbd2
          ts_rc2in(isp,jsp) = rc2in
          ts_rc2out(isp,jsp) = rc2out
          ts_eta(isp,jsp) = eta
          ts_delta(isp,jsp) = delta
          ts_alpha(isp,jsp) = alpha
          ts_beta(isp,jsp) = beta
          ts_h(isp,jsp) = h
          ts_c1(isp,jsp) = c1
          ts_c2(isp,jsp) = c2
          ts_c3(isp,jsp) = c3
          ts_c4(isp,jsp) = c4
          ts_c5(isp,jsp) = c5
          ts_rc3in(isp,jsp) = rc3in
          ts_rc3out(isp,jsp) = rc3out
          interact(isp,jsp) = .true.

!.....Symmetrize
          if( isp.ne.jsp ) then
            ts_a(jsp,isp) = a
            ts_b(jsp,isp) = b
            ts_lmbd1(jsp,isp) = lmbd1
            ts_lmbd2(jsp,isp) = lmbd2
            ts_rc2in(jsp,isp) = rc2in
            ts_rc2out(jsp,isp) = rc2out
            ts_eta(jsp,isp) = eta
            ts_delta(jsp,isp) = delta
            ts_alpha(jsp,isp) = alpha
            ts_beta(jsp,isp) = beta
            ts_h(jsp,isp) = h
            ts_c1(jsp,isp) = c1
            ts_c2(jsp,isp) = c2
            ts_c3(jsp,isp) = c3
            ts_c4(jsp,isp) = c4
            ts_c5(jsp,isp) = c5
            ts_rc3in(jsp,isp) = rc3in
            ts_rc3out(jsp,isp) = rc3out
            interact(jsp,isp) = .true.
          endif

          if( iprint.ge.ipl_basic ) then
            write(6,'(2a5,4es12.4,2f6.3)') trim(cspi),trim(cspj), a, b, lmbd1,lmbd2,rc2in,rc2out
            write(6,'(15x,7es12.4)') eta,delta,alpha,beta,h
            write(6,'(15x,5es12.4,2f6.3)') c1,c2,c3,c4,c5,rc3in,rc3out
          endif
!!$          goto 10
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
      ts_eta(:,:) = ted_eta(1)
      ts_delta(:,:) = ted_delta(1)
      ts_alpha(:,:) = ted_alpha(1)
      ts_beta(:,:) = ted_beta(1)
      ts_h(:,:) = ted_h(1)
      ts_rc2in(:,:) = ted_r1(1)
      ts_rc2out(:,:) = ted_r2(1)
      ts_rc3in(:,:) = ts_rc2in(:,:)
      ts_rc3out(:,:) = ts_rc2out(:,:)
      ts_c1(:,:) = ted_c1(1)
      ts_c2(:,:) = ted_c2(1)
      ts_c3(:,:) = ted_c3(1)
      ts_c4(:,:) = ted_c4(1)
      ts_c5(:,:) = ted_c5(1)
      ts_f0(:,:) = ted_f0(1)
      interact(:,:) = .true.

    else
      call mpi_bcast(ts_a,nspmax*nspmax,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ts_b,nspmax*nspmax,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ts_lmbd1,nspmax*nspmax,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ts_lmbd2,nspmax*nspmax,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ts_eta,nspmax*nspmax,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ts_delta,nspmax*nspmax,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ts_alpha,nspmax*nspmax,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ts_beta,nspmax*nspmax,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ts_h,nspmax*nspmax,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ts_rc2in,nspmax*nspmax,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ts_rc2out,nspmax*nspmax,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ts_rc3in,nspmax*nspmax,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ts_rc3out,nspmax*nspmax,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ts_c1,nspmax*nspmax,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ts_c2,nspmax*nspmax,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ts_c3,nspmax*nspmax,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ts_c4,nspmax*nspmax,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ts_c5,nspmax*nspmax,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(ts_f0,nspmax*nspmax,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(interact,nspmax*nspmax,mpi_logical,0,mpi_world,ierr)
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
    ts_eta(is,is) = ted_eta(ite0) +(ted_eta(ite1) -ted_eta(ite0)) *x
    ts_delta(is,is) = ted_delta(ite0) +(ted_delta(ite1) -ted_delta(ite0)) *x
    ts_alpha(is,is) = ted_alpha(ite0) +(ted_alpha(ite1) -ted_alpha(ite0)) *x
    ts_beta(is,is) = ted_beta(ite0) +(ted_beta(ite1) -ted_beta(ite0)) *x
    ts_h(is,is) = ted_h(ite0) +(ted_h(ite1) -ted_h(ite0)) *x
    ts_c1(is,is) = ted_c1(ite0) +(ted_c1(ite1) -ted_c1(ite0)) *x
    ts_c2(is,is) = ted_c2(ite0) +(ted_c2(ite1) -ted_c2(ite0)) *x
    ts_c3(is,is) = ted_c3(ite0) +(ted_c3(ite1) -ted_c3(ite0)) *x
    ts_c4(is,is) = ted_c4(ite0) +(ted_c4(ite1) -ted_c4(ite0)) *x
    ts_c5(is,is) = ted_c5(ite0) +(ted_c5(ite1) -ted_c5(ite0)) *x
    ts_f0(is,is) = ted_f0(ite0) +(ted_f0(ite1) -ted_f0(ite0)) *x
    ts_rc2in(is,is) = ted_r1(ite0) +(ted_r1(ite1) -ted_r1(ite0)) *x
    ts_rc2out(is,is) = ted_r2(ite0) +(ted_r2(ite1) -ted_r2(ite0)) *x

!!$    ts_r(is,is) = (ts_r1(is,is) +ts_r2(is,is))/2
!!$    ts_d(is,is) = (ts_r2(is,is) -ts_r1(is,is))/2
    ts_rc(is,is) = ts_rc2out(is,is)
    ts_rc2(is,is)= ts_rc(is,is)**2
    return
  end subroutine set_ted_params
!=======================================================================
  subroutine set_orig_params()
!
!  Set original parameter set taken from Kumagai et al.[2], if not given.
!

    ts_type = 'default'
!.....Cutoff type: (1) Tersoff,  (3) Murty (default)
    ts_fc_type = 3
    
!.....Original parameters from [2]
    ts_a(:,:) = 3281.5905d0
    ts_b(:,:) = 121.00047d0
    ts_lmbd1(:,:) = 3.2300135d0
    ts_lmbd2(:,:) = 1.345797d0
    ts_eta(:,:) = 1.0d0
    ts_delta(:,:) = 0.53298909d0
    ts_alpha(:,:) = 2.3890327d0
    ts_beta(:,:) = 1d0
!!$    ts_beta_ters = 1d0  ! beta of original Tersoff [1], must be 1.0 for [2]
    ts_c1(:,:) = 0.20173476d0
    ts_c2(:,:) = 730418.72d0
    ts_c3(:,:) = 1000000d0
    ts_c4(:,:) = 1.0d0
    ts_c5(:,:) = 26.0d0
    ts_h(:,:) = -0.365
    ts_rc2in(:,:) = 2.7d0
    ts_rc2out(:,:) = 3.3d0
    ts_rc3in(:,:) = ts_rc2in(:,:)
    ts_rc3out(:,:) = ts_rc2out(:,:)
    ts_f0(:,:) = 0d0
    interact(:,:) = .true.
    return
  end subroutine set_orig_params
end module tersoff
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd lib"
!     End:
