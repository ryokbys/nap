module UF3
!-----------------------------------------------------------------------
!                     Last modified: <2024-10-23 14:07:17 KOBAYASHI Ryo>
!-----------------------------------------------------------------------
!  Parallel implementation of Ultra-Fast Force-Field (UF3) for pmd
!    - 2024.09.02 by R.K., start to implement
!-----------------------------------------------------------------------
  use pmdvars,only: nspmax
  use util,only: csp2isp, num_data
  use memory,only: accum_mem
  use vector,only: dot
  implicit none
  include "mpif.h"
  include './const.h'
  include "./params_unit.h"
  save

  character(len=128):: paramsdir = '.'
!.....parameter file name
  character(128),parameter:: cpfname= 'in.params.uf3'
  character(128),parameter:: ccfname='in.const.uf3'
  integer,parameter:: ioprms = 50

  logical:: lprmset_uf3 = .false.

!.....uf2 parameters
  type prm2
    character(2):: cb, csi, csj, cknot
!  cb: NA, 2B or 3B
!  csi,csj,csk: species name
!  cknot: nk (non-uniform knot spacing) or uk (uniform knot spacing)
    integer:: nklead, nktrail
    integer:: nknot, ncoef
    real(8):: rc,rc2
    real(8),allocatable:: knots(:), coefs(:)
  end type prm2

!.....uf3 parameters
  type prm3
    character(2):: cb, csi, csj, csk, cknot
!  cb: NA, 2B or 3B
!  csi,csj,csk: species name
!  cknot: nk (non-uniform knot spacing) or uk (uniform knot spacing)
    integer:: nklead, nktrail
    integer:: nknij, nknik, nknjk, ncfij, ncfik, ncfjk
    real(8):: rcij, rcik, rcjk, rcij2, rcik2, rcjk2
    real(8),allocatable:: knij(:), knik(:), knjk(:), coefs(:,:,:)
  end type prm3

  integer:: n2b, n3b
  type(prm2),allocatable:: prm2s(:)
  type(prm3),allocatable:: prm3s(:)
  logical:: has_trios = .false.
  real(8):: rcmax = 0.0d0
  
!.....Map of pairs (trios) to parameter set id
  integer:: interact2(nspmax,nspmax), interact3(nspmax,nspmax,nspmax)
!.....Cutoffs
  real(8):: rc2_3b(nspmax,nspmax)
  
!.....constants
  integer:: nelem,nexp,nsp

contains
  subroutine read_params_uf3(myid,mpi_world,iprint)
!
!  Read parameters of uf3 potential from in.params.uf3 file that is given
!  by uf3/lammps_plugin/scripts/generate_uf3_lammps_pots.py, 
!-----------------------------------------------------------------------
!  #UF3 POT UNITS: metal DATE: 2024-09-12 17:12:36 AUTHOR: RK CITATION:
!  2B W W 0 3 nk
!  5.5  22
!  0.001 0.001 0.001 0.001 0.3676 0.7342 ... 5.1334 5.5 5.5 5.5 5.5
!  18
!  36.82 36.82 26.69 ... -0.032 0 0 0
!  #
!  #UF3 POT UNITS: metal DATE: 2024-09-12 17:12:36 AUTHOR: RK CITATION:
!  3B W W W 0 3 nk
!  7.0 3.5 3.5 19 13 13
!  1.5 1.5 1.5 1.96 2.42 ... 6.54 7 7 7
!  1.5 1.5 1.5 1.83 2.17 ... 3.17 3.5 3.5 3.5
!  1.5 1.5 1.5 1.83 2.17 ... 3.17 3.5 3.5 3.5
!  9 9 15
!  0.0 0.0 0.0 0.0 ... 0.0 0.0 0.0
!  ...
!  #
!-----------------------------------------------------------------------
!  Each pair or trio is sandwitched with lines beginning with "#".
!  For trios, in the case of knot information the order is jk,ik,ij,
!  but for coefficient information the order is ij,ik,jk.
!  And there is a limitation in the current lammps implementation version
!  such that r_max_jk = 2*r_max_ij = 2_r_max_ik.
!
    implicit none
    integer,intent(in):: myid,mpi_world,iprint
    
    integer:: itmp,ierr,i,j,i2b,i3b
    integer:: nklead, nktrail
!  nklead, nktrail: num of leading or trailing knots
    logical:: lexist
    character:: fname*128, cmode*4, cb*2, csi*2, csj*2, csk*2, &
         cknot*2, ctmp*128, cline*128
!  cmode: none or read
    
    if( myid == 0 ) then
      if( iprint >= ipl_basic ) print '(/,a)',' Read UF3 parameters...'
      fname = trim(paramsdir)//'/'//trim(cpfname)
!.....read parameters at the 1st call
      inquire(file=trim(fname),exist=lexist)
      if( .not. lexist ) then
        write(6,'(a)') '   [Error] '//trim(fname)//' does not exist !!!.'
!!$        call mpi_finalize(ierr)
        stop
      endif
      cmode = 'none'
      n2b = 0
      n3b = 0
      interact2(:,:) = -1
      interact3(:,:,:) = -1
      open(ioprms,file=trim(fname),status='old')
!.....1st, count number of 2B and 3B entries to allocate type objects
      do while(.true.)
        read(ioprms,'(a)',end=10) cline
        if( cline(1:1) == '#' ) then
          if( cline(1:3) == '#UF' ) then
            cmode = 'read'
            cb = 'NA'
          else
            cmode = 'none'
          endif
          cycle
        endif
        if( trim(cmode) == 'none' ) cycle
        if( cb == 'NA' ) then  ! if cb is not read yet
          read(cline,*,iostat=ierr) cb
          if( cb == '2B' ) then
            n2b = n2b + 1
          else if( cb == '3B' ) then
            n3b = n3b + 1
            has_trios = .true.
          endif
! do nothing if cb is already set...
        endif
!        if( ierr /= 0 ) cycle
      end do  ! finished counting 2B & 3B entries
10    continue
      if( iprint >= ipl_debug ) then
        print *,'  n2b, n3b = ',n2b,n3b
        print *,'  has_trios= ',has_trios
      endif
!.....allocate lists of uf2prms and uf3prms
      if( allocated(prm2s) ) deallocate(prm2s)
      if( allocated(prm3s) ) deallocate(prm3s)
      if( n2b /= 0 ) allocate(prm2s(n2b))
      if( n3b /= 0 ) allocate(prm3s(n3b))
      rewind(ioprms)
      i2b = 0
      i3b = 0
      do while(.true.)
        read(ioprms,'(a)',end=20) cline
        if( cline(1:1) == '#') then
          if( cline(1:3) == '#UF') then
            cmode = 'read'
            cb = 'NA'
          else
            cmode = 'none'
          endif
          cycle
        endif
        if( trim(cmode) == 'none' ) cycle
        if( cb == 'NA' ) then  ! if cb is not read yet
          read(cline,*,iostat=ierr) cb
          if( cb == '2B' ) then
            backspace(ioprms)
            i2b = i2b +1
            call read_2b(prm2s(i2b),i2b)
            if( iprint >= ipl_basic ) call print_2b(prm2s(i2b))
          else if( cb == '3B' ) then
            backspace(ioprms)
            i3b = i3b +1
            call read_3b(prm3s(i3b),i3b)
            if( iprint >= ipl_basic ) call print_3b(prm3s(i3b))
          endif
        endif
!        if( ierr /= 0 ) cycle

      enddo  ! while(.true.)
20    continue  ! when the file reached the end
    endif

    call bcast_uf3_params(mpi_world,myid)

    return
  end subroutine read_params_uf3
!=======================================================================
  subroutine read_2b(ps,i2b)
!
!  Read 2B part of file via ioprms.
!  Assuming that the starting line of IOPRMS is the 1st line after
!  the line starting with '#UF3'.
!
    type(prm2),intent(out):: ps
    integer,intent(in):: i2b
    integer:: i, isp, jsp
    
    read(ioprms,*) ps%cb, ps%csi, ps%csj, ps%nklead, ps%nktrail, ps%cknot
    if( ps%cb /= '2B' ) stop 'ERROR: CB should be 2B.'
    read(ioprms,*) ps%rc, ps%nknot
    ps%rc2 = ps%rc**2
    if( allocated(ps%knots) ) deallocate(ps%knots)
    allocate(ps%knots(ps%nknot))
    read(ioprms,*) (ps%knots(i), i=1,ps%nknot)
    read(ioprms,*) ps%ncoef
    if( allocated(ps%coefs) ) deallocate(ps%coefs)
    allocate(ps%coefs(ps%ncoef))
    read(ioprms,*) (ps%coefs(i), i=1,ps%ncoef)

    isp = csp2isp(ps%csi)
    jsp = csp2isp(ps%csj)
    interact2(isp,jsp) = i2b
    interact2(jsp,isp) = i2b
    rcmax = max(ps%rc,rcmax)
  end subroutine read_2b
!=======================================================================
  subroutine read_3b(ps,i3b)
!
!  Read 3B part of file via ioprms
!  Assuming that the starting line of IOPRMS is the 1st line after
!  the line starting with '#UF3'.
!
    type(prm3),intent(out):: ps
    integer,intent(in):: i3b
    integer:: i,j,k,isp,jsp,ksp
    
    read(ioprms,*) ps%cb, ps%csi, ps%csj, ps%csk, ps%nklead, ps%nktrail, ps%cknot
    if( ps%cb /= '3B' ) stop 'ERROR: CB should be 3B.'
    read(ioprms,*) ps%rcjk, ps%rcij, ps%rcik, ps%nknjk, ps%nknij, ps%nknik
    ps%rcjk2 = ps%rcjk**2
    ps%rcik2 = ps%rcik**2
    ps%rcij2 = ps%rcij**2
    if( allocated(ps%knij) ) deallocate(ps%knij, ps%knik, ps%knjk)
    allocate(ps%knij(ps%nknij), ps%knik(ps%nknik), ps%knjk(ps%nknjk))
    read(ioprms,*) (ps%knjk(i), i=1,ps%nknjk)
    read(ioprms,*) (ps%knij(i), i=1,ps%nknij)
    read(ioprms,*) (ps%knik(i), i=1,ps%nknik)
    read(ioprms,*) ps%ncfij, ps%ncfik, ps%ncfjk
    if( allocated(ps%coefs) ) deallocate(ps%coefs)
    allocate(ps%coefs(ps%ncfij, ps%ncfik, ps%ncfjk))
    do i=1, ps%ncfij
      do j=1, ps%ncfik
        read(ioprms,*) (ps%coefs(i,j,k), k=1,ps%ncfjk)
      enddo
    enddo

    isp = csp2isp(ps%csi)
    jsp = csp2isp(ps%csj)
    ksp = csp2isp(ps%csk)
    interact3(isp,jsp,ksp) = i3b
    interact3(isp,ksp,jsp) = i3b
    rcmax = max(ps%rcij,rcmax)
    rcmax = max(ps%rcik,rcmax)
    rc2_3b(isp,jsp) = max(ps%rcij2, rc2_3b(isp,jsp))
    rc2_3b(isp,ksp) = max(ps%rcik2, rc2_3b(isp,ksp))
    rc2_3b(jsp,isp) = rc2_3b(isp,jsp)
    rc2_3b(ksp,isp) = rc2_3b(isp,ksp)
  end subroutine read_3b
!=======================================================================
  subroutine bcast_uf3_params(mpi_world,myid)
!
!  Broadcast 2B & 3B parameters.
!
    integer,intent(in):: mpi_world, myid
    integer:: i, i2b, i3b, ierr
    type(prm2):: p2
    type(prm3):: p3

    call mpi_bcast(n2b, 1, mpi_integer, 0,mpi_world,ierr)
    call mpi_bcast(n3b, 1, mpi_integer, 0,mpi_world,ierr)
    if( .not. allocated(prm2s) ) allocate(prm2s(n2b))
    if( .not. allocated(prm3s) ) allocate(prm3s(n3b))
    
    do i2b=1,n2b
      p2 = prm2s(i2b)
      call mpi_bcast(p2%cb,2,mpi_character,0,mpi_world,ierr)
      call mpi_bcast(p2%csi,2,mpi_character,0,mpi_world,ierr)
      call mpi_bcast(p2%csj,2,mpi_character,0,mpi_world,ierr)
      call mpi_bcast(p2%cknot,2,mpi_character,0,mpi_world,ierr)
      call mpi_bcast(p2%nklead,1,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(p2%nktrail,1,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(p2%nknot,1,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(p2%ncoef,1,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(p2%rc,1,mpi_real8,0,mpi_world,ierr)
      if( .not.allocated(p2%knots) ) allocate(p2%knots(p2%nknot))
      if( .not.allocated(p2%coefs) ) allocate(p2%coefs(p2%ncoef))
      call mpi_bcast(p2%knots,p2%nknot,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(p2%coefs,p2%ncoef,mpi_real8,0,mpi_world,ierr)
    enddo

    do i3b=1,n3b
      p3 = prm3s(i3b)
      call mpi_bcast(p3%cb,2,mpi_character,0,mpi_world,ierr)
      call mpi_bcast(p3%csi,2,mpi_character,0,mpi_world,ierr)
      call mpi_bcast(p3%csj,2,mpi_character,0,mpi_world,ierr)
      call mpi_bcast(p3%csk,2,mpi_character,0,mpi_world,ierr)
      call mpi_bcast(p3%cknot,2,mpi_character,0,mpi_world,ierr)
      call mpi_bcast(p3%nklead,1,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(p3%nktrail,1,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(p3%nknij,1,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(p3%nknik,1,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(p3%nknjk,1,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(p3%ncfij,1,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(p3%ncfik,1,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(p3%ncfjk,1,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(p3%rcij,1,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(p3%rcik,1,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(p3%rcjk,1,mpi_real8,0,mpi_world,ierr)
      if( .not.allocated(p3%knij) ) allocate(p3%knij(p3%nknij))
      if( .not.allocated(p3%knik) ) allocate(p3%knik(p3%nknik))
      if( .not.allocated(p3%knjk) ) allocate(p3%knjk(p3%nknjk))
      if( .not.allocated(p3%coefs)) allocate(p3%coefs(p3%ncfij, p3%ncfik, p3%ncfjk))
      call mpi_bcast(p3%knij,p3%nknij,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(p3%knik,p3%nknik,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(p3%knjk,p3%nknjk,mpi_real8,0,mpi_world,ierr)
      call mpi_bcast(p3%coefs,p3%ncfjk*p3%ncfik*p3%ncfij,mpi_real8,0,mpi_world,ierr)
    enddo

    call mpi_bcast(has_trios, 1, mpi_logical, 0,mpi_world,ierr)
    call mpi_bcast(interact2, nspmax*nspmax, mpi_integer, 0,mpi_world,ierr)
    call mpi_bcast(interact3, nspmax**3, mpi_integer, 0,mpi_world,ierr)
    
  end subroutine bcast_uf3_params
!=======================================================================
  subroutine force_uf3(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rcin,lspr &
       ,mpi_world,myid,epi,epot,lstrs,iprint,l1st)
!
!  UF3 implementation without using recursive function of b-spline.
!  
!  TODO: More efficient B-spline implementation is available (see, lammps src/ML-UF3/pair_uf3.cpp).
!        But that is a bit complicated, use simple b_spl() routine for now.
!
    use util, only: itotOf
    implicit none
    integer,intent(in):: namax,natm,nnmax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_world,myid,lspr(0:nnmax,namax),nex(3)
    real(8),intent(in):: ra(3,namax),tag(namax) &
         ,h(3,3),hi(3,3),sv(3,6)
    real(8),intent(inout):: rcin
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: l1st 
    logical:: lstrs

!.....local
    integer:: ia,ja,ka,jj,kk,l,is,nr2,n,nij3,inc,nik3,njk3,&
         nik,njk,nij,itot,jtot,ktot,i2b,i3b,js,ks,jsp,ksp,ierr, &
         ixyz,jxyz,lij,lik,ljk
    real(8):: epotl2,epotl3,epot2,epot3,tmp,tmp2,bij(-3:0),dbij(-3:0), &
         bij3(-3:0),dbij3(-3:0),bik3(-3:0),dbik3(-3:0),bjk3(-3:0),dbjk3(-3:0), &
         c2t,c3t
    real(8):: xi(3),xj(3),xk(3),xij(3),xik(3),xjk(3),rij(3),rik(3),&
         rjk(3),dij2,dij,dik2,dik,djk2,djk,drijj(3),drikk(3),&
         drjkk(3),tmpij(3),tmpik(3),tmpjk(3)
    real(8),save,allocatable:: aal2(:,:),aal3(:,:),strsl(:,:,:)
    real(8),save:: rcin2
    integer,save,allocatable:: ls3b(:)

    type(prm2):: p2
    type(prm3):: p3

    if( l1st ) then
      if( allocated(aal2) ) deallocate(aal2,aal3,strsl,ls3b)
      allocate(aal2(3,namax),aal3(3,namax),strsl(3,3,namax),ls3b(0:nnmax))
      if( rcin < rcmax ) then
        if( myid == 0 ) then
          write(6,'(1x,a)') "ERROR: Cutoff radius is not appropriate !!!"
          write(6,'(1x,a,f0.3)') "  rc should be longer than ", rcmax
        endif
        call mpi_finalize(ierr)
        stop
      endif
      rcin2 = rcin*rcin
    endif

    if( size(aal2) < 3*namax ) then
      deallocate(aal2,aal3,strsl)
      allocate(aal2(3,namax),aal3(3,namax),strsl(3,3,namax))
    endif

    if( size(ls3b) < nnmax+1 ) then
      deallocate(ls3b)
      allocate(ls3b(0:nnmax))
    endif

    aal2(:,:) = 0d0
    aal3(:,:) = 0d0
    strsl(:,:,:) = 0d0
    epotl2 = 0d0
    epotl3 = 0d0
!$omp parallel
!$omp do private(ia,is,xi,jj,ja,js,i2b,p2,xj,xij,rij,dij2,dij,drijj,nr2, &
!$omp      n,bij,c2t,tmp,dbij,tmp2,ixyz,jxyz, &
!$omp      jsp,ksp,i3b,p3,nij3,inc,bij3,dbij3, &
!$omp      kk,ka,ks,xk,xik,xjk,rik,dik2,dik,rjk,djk2,djk, &
!$omp      drikk,drjkk,nik3,njk3,nik,bik3,dbik3,njk,bjk3,dbjk3, &
!$omp      l,nij,c3t,tmpij,tmpik,tmpjk) &
!$omp      reduction(+:epotl2,epotl3)
    do ia=1,natm
      is = int(tag(ia))
      xi(1:3) = ra(1:3,ia)
      ls3b(0) = 0
      do jj=1,lspr(0,ia)
        ja = lspr(jj,ia)
!!$        if( ja <= ia ) cycle
        js = int(tag(ja))
!.....Pair terms
        i2b = interact2(is,js)
        if( i2b <= 0 ) cycle
        p2 = prm2s(i2b)
        xj(1:3) = ra(1:3,ja)
        xij(1:3) = xj(1:3) -xi(1:3)
        rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij2 = rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
        if( dij2 > p2%rc2 ) cycle
        dij = sqrt(dij2)
!.....Make short-distance pair-list for 3-body term
        if( has_trios ) then
          if( dij2 < rc2_3b(is,js) ) then
            ls3b(0) = ls3b(0) +1
            ls3b(ls3b(0)) = ja
          endif
        endif
        drijj(1:3) = rij(1:3)/dij
        call b_spl(dij,p2%knots,p2%nknot,nr2,bij,dbij)
        do lij = -3,0
          n = nr2 +lij
          if( n < 1 .or. n > p2%nknot-4 ) cycle
          c2t = p2%coefs(n)
          tmp = c2t *bij(lij)
!.....Energy
          epi(ia) = epi(ia) +tmp
          epotl2 = epotl2 +tmp
!.....Forces
          tmp2 = c2t *dbij(lij)
          do ixyz=1,3
!$omp atomic
            aal2(ixyz,ia) = aal2(ixyz,ia) +drijj(ixyz)*tmp2
!$omp atomic
            aal2(ixyz,ja) = aal2(ixyz,ja) -drijj(ixyz)*tmp2
          enddo
!.....Stresses
          do ixyz=1,3
            do jxyz=1,3
!$omp atomic
              strsl(jxyz,ixyz,ia)= strsl(jxyz,ixyz,ia) &
                   -0.5d0 *tmp2*rij(ixyz)*drijj(jxyz)
!$omp atomic
              strsl(jxyz,ixyz,ja)= strsl(jxyz,ixyz,ja) &
                   -0.5d0 *tmp2*rij(ixyz)*drijj(jxyz)
            enddo
          enddo
        enddo  ! lij
      enddo

!.....Trio part is separated from pair part,
!.....which may be slower because of double computation of dij,
!.....but this code is a bit simpler.
      if( .not.has_trios ) cycle
      do jsp=1,nspmax
        do ksp=jsp,nspmax
          i3b = interact3(is,jsp,ksp)
          if( i3b <= 0 ) cycle
          p3 = prm3s(i3b)
!!$          do jj=1,lspr(0,ia)
!!$            ja = lspr(jj,ia)
          do jj=1,ls3b(0)
            ja = ls3b(jj)
            js = int(tag(ja))
            if( js /= jsp ) cycle
            xj(1:3) = ra(1:3,ja)
            xij(1:3) = xj(1:3) -xi(1:3)
            rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
            dij2 = rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
            if( dij2 > p3%rcij2 ) cycle
            dij = sqrt(dij2)
            drijj(1:3) = rij(1:3)/dij
            call b_spl(dij, p3%knij, p3%nknij, nij3, bij3, dbij3)
            
!!$            do kk=1,lspr(0,ia)
!!$              ka = lspr(kk,ia)
            do kk=1,ls3b(0)
              ka = ls3b(kk)
              if( jsp == ksp .and. ka <= ja ) cycle
              ks = int(tag(ka))
              if( ks /= ksp ) cycle
              xk(1:3) = ra(1:3,ka)
              xik(1:3) = xk(1:3) -xi(1:3)
              xjk(1:3) = xk(1:3) -xj(1:3)
              rik(1:3) = h(1:3,1)*xik(1) +h(1:3,2)*xik(2) +h(1:3,3)*xik(3)
              dik2 = rik(1)*rik(1) +rik(2)*rik(2) +rik(3)*rik(3)
              if( dik2 > p3%rcik2 ) cycle
              dik = sqrt(dik2)
              rjk(1:3) = h(1:3,1)*xjk(1) +h(1:3,2)*xjk(2) +h(1:3,3)*xjk(3)
              djk2 = rjk(1)*rjk(1) +rjk(2)*rjk(2) +rjk(3)*rjk(3)
              if( djk2 > p3%rcjk2 ) cycle
              djk = sqrt(djk2)
              drikk(1:3) = rik(1:3)/dik
              drjkk(1:3) = rjk(1:3)/djk
!.....B-spline part
              call b_spl(dik, p3%knik, p3%nknik, nik3, bik3, dbik3)
              call b_spl(djk, p3%knjk, p3%nknjk, njk3, bjk3, dbjk3)
              do lik = -3,0
                nik = nik3 +lik
                if( nik < 1 .or. nik > p3%nknik-4 ) cycle
                do ljk = -3,0
                  njk = njk3 +ljk
                  if( njk < 1 .or. njk > p3%nknjk-4 ) cycle
                  do lij = -3,0
                    nij = nij3 +lij
                    if( nij < 1 .or. nij > p3%nknij-4 ) cycle
!.....Energy
                    c3t = p3%coefs(nij,nik,njk)
                    tmp = c3t*bij3(lij)*bik3(lik)*bjk3(ljk)
                    epi(ia) = epi(ia) +tmp
                    epotl3 = epotl3 +tmp
!.....Force
                    tmpij(1:3) = dbij3(lij)*bik3(lik)*bjk3(ljk)*drijj(1:3)
                    tmpik(1:3) = bij3(lij)*dbik3(lik)*bjk3(ljk)*drikk(1:3)
                    tmpjk(1:3) = bij3(lij)*bik3(lik)*dbjk3(ljk)*drjkk(1:3)
                    do ixyz=1,3
!$omp atomic
                      aal3(ixyz,ia)= aal3(ixyz,ia) +c3t*(tmpij(ixyz) +tmpik(ixyz))
!$omp atomic
                      aal3(ixyz,ja)= aal3(ixyz,ja) +c3t*(-tmpij(ixyz) +tmpjk(ixyz))
!$omp atomic
                      aal3(ixyz,ka)= aal3(ixyz,ka) +c3t*(-tmpik(ixyz) -tmpjk(ixyz))
                    enddo
!.....Stresses
                    do jxyz=1,3
                      do ixyz=1,3
!$omp atomic
                        strsl(ixyz,jxyz,ia)= strsl(ixyz,jxyz,ia) &
                             -0.5d0 *c3t *(xij(jxyz)*tmpij(ixyz) +xik(jxyz)*tmpik(ixyz))
!$omp atomic
                        strsl(ixyz,jxyz,ja)= strsl(ixyz,jxyz,ja) &
                             -0.5d0 *c3t *(xij(jxyz)*tmpij(ixyz) +xjk(jxyz)*tmpjk(ixyz))
!$omp atomic
                        strsl(ixyz,jxyz,ka)= strsl(ixyz,jxyz,ka) &
                             -0.5d0 *c3t *(xik(jxyz)*tmpik(ixyz) +xjk(jxyz)*tmpjk(ixyz))
                      enddo
                    enddo
                    
                  enddo  ! lij
                enddo  ! ljk
              enddo  ! lik

            enddo  ! kk
          enddo  ! jj
        enddo  ! ksp
      enddo  ! jsp

    enddo ! ia
!$omp end do
!$omp end parallel

    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,aal2,3)
    if( has_trios ) call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,aal3,3)
    aa(1:3,1:natm) = aa(1:3,1:natm) +aal2(1:3,1:natm) +aal3(1:3,1:natm)

    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,strsl,9)
    strs(1:3,1:3,1:natm) = strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)

!-----gather epot
    epot2 = 0d0
    epot3 = 0d0
    call mpi_allreduce(epotl2,epot2,1,mpi_real8,mpi_sum,mpi_world,ierr)
    if( has_trios ) call mpi_allreduce(epotl3,epot3,1,mpi_real8,mpi_sum,mpi_world,ierr)
    epot= epot +epot2 +epot3
    if( myid == 0 .and. iprint > 2 ) &
         print '(a,2es12.4)',' force_uf3 epot2,epot3 = ',epot2,epot3

    return
  end subroutine force_uf3
!=======================================================================
  subroutine force_uf3_rec(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rcin,lspr &
       ,mpi_world,myid,epi,epot,lstrs,iprint,l1st)
!
!  Recursive implementation of force_uf3.
!
    use util, only: itotOf
    implicit none
    integer,intent(in):: namax,natm,nnmax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_world,myid,lspr(0:nnmax,namax),nex(3)
    real(8),intent(in):: ra(3,namax),tag(namax) &
         ,h(3,3),hi(3,3),sv(3,6)
    real(8),intent(inout):: rcin
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: l1st 
    logical:: lstrs

!.....local
    integer:: ia,ja,ka,jj,kk,l,is,nr2,n,nij3,inc,nik3,njk3,&
         nik,njk,nij,itot,jtot,ktot,i2b,i3b,js,ks,jsp,ksp,ierr, &
         ixyz,jxyz,lij,lik,ljk
    real(8):: epotl2,epotl3,epot2,epot3,tmp,tmp2,bij,dbij, &
         bij3(4),dbij3(4),bik3,dbik3,bjk3,dbjk3, &
         c2t,c3t
    real(8):: xi(3),xj(3),xk(3),xij(3),xik(3),xjk(3),rij(3),rik(3),&
         rjk(3),dij2,dij,dik2,dik,djk2,djk,drijj(3),drikk(3),&
         drjkk(3),tmpij(3),tmpik(3),tmpjk(3)
    real(8),save,allocatable:: aal2(:,:),aal3(:,:),strsl(:,:,:)
    real(8),save:: rcin2

    type(prm2):: p2
    type(prm3):: p3

    if( l1st ) then
      if( allocated(aal2) ) deallocate(aal2,strsl)
      allocate(aal2(3,namax),aal3(3,namax),strsl(3,3,namax))
      if( rcin < rcmax ) then
        if( myid == 0 ) then
          write(6,'(1x,a)') "ERROR: Cutoff radius is not appropriate !!!"
          write(6,'(1x,a,f0.3)') "  rc should be longer than ", rcmax
        endif
        call mpi_finalize(ierr)
        stop
      endif
      rcin2 = rcin*rcin
    endif

    if( size(aal2) < 3*namax ) then
      deallocate(aal2,aal3,strsl)
      allocate(aal2(3,namax),aal3(3,namax),strsl(3,3,namax))
    endif

    aal2(:,:) = 0d0
    aal3(:,:) = 0d0
    strsl(1:3,1:3,1:namax) = 0d0
    epotl2 = 0d0
    epotl3 = 0d0
!$omp parallel
!$omp do private(ia,is,xi,jj,ja,js,i2b,p2,xj,xij,rij,dij2,dij,drijj,nr2, &
!$omp      n,bij,c2t,tmp,dbij,tmp2,ixyz,jxyz, &
!$omp      jsp,ksp,i3b,p3,nij3,inc,bij3,dbij3, &
!$omp      kk,ka,ks,xk,xik,xjk,rik,dik2,dik,rjk,djk2,djk, &
!$omp      drikk,drjkk,nik3,njk3,nik,bik3,dbik3,njk,bjk3,dbjk3, &
!$omp      l,nij,c3t,tmpij,tmpik,tmpjk) &
!$omp      reduction(+:epotl2,epotl3)
    do ia=1,natm
      is = int(tag(ia))
      xi(1:3) = ra(1:3,ia)
      do jj=1,lspr(0,ia)
        ja = lspr(jj,ia)
!!$        if( ja <= ia ) cycle
        js = int(tag(ja))
!.....Pair terms
        i2b = interact2(is,js)
        if( i2b <= 0 ) cycle
        p2 = prm2s(i2b)
        xj(1:3) = ra(1:3,ja)
        xij(1:3) = xj(1:3) -xi(1:3)
        rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij2 = rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
        if( dij2 > p2%rc2 ) cycle
        dij = sqrt(dij2)
        drijj(1:3) = rij(1:3)/dij
!.....Look for nr2 from knots data
        nr2 = knot_index(dij, p2%nknot, p2%knots)
        do n=max(1,nr2-3),min(nr2,p2%nknot-4)
          bij = b_spl_rec(n,3,dij,p2%knots,p2%nknot)
          c2t = p2%coefs(n)
          tmp = c2t *bij
          epi(ia) = epi(ia) +tmp
          epotl2 = epotl2 + tmp
!.....Forces
          dbij = db_spl(n,3,dij,p2%knots,p2%nknot)
          tmp2 = dbij*c2t
          do ixyz=1,3
!$omp atomic
            aal2(ixyz,ia) = aal2(ixyz,ia) +drijj(ixyz)*tmp2
!$omp atomic
            aal2(ixyz,ja) = aal2(ixyz,ja) -drijj(ixyz)*tmp2
          enddo
!.....Stresses
          do ixyz=1,3
            do jxyz=1,3
!$omp atomic
              strsl(jxyz,ixyz,ia)= strsl(jxyz,ixyz,ia) &
                   -0.5d0 *tmp2*rij(ixyz)*drijj(jxyz)
!$omp atomic
              strsl(jxyz,ixyz,ja)= strsl(jxyz,ixyz,ja) &
                   -0.5d0 *tmp2*rij(ixyz)*drijj(jxyz)
            enddo
          enddo
        enddo
      enddo

!.....Trio part is separated from pair part,
!.....which may be slower because of double computation of dij,
!.....but this code is a bit simpler.
      if( .not.has_trios ) cycle
      do jsp=1,nspmax
        do ksp=jsp,nspmax
          i3b = interact3(is,jsp,ksp)
          if( i3b <= 0 ) cycle
          p3 = prm3s(i3b)
          do jj=1,lspr(0,ia)
            ja = lspr(jj,ia)
            js = int(tag(ja))
            if( js /= jsp ) cycle
            xj(1:3) = ra(1:3,ja)
            xij(1:3) = xj(1:3) -xi(1:3)
            rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
            dij2 = rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
            if( dij2 > p3%rcij2 ) cycle
            dij = sqrt(dij2)
            drijj(1:3) = rij(1:3)/dij
!!$            call b_spl(dij, p3%knij, p3%nknij, nij3, bij3, dbij3)
            nij3 = knot_index(dij, p3%nknij, p3%knij)
            inc = 0
            bij3(:) = 0d0
            dbij3(:) = 0d0
            do n=nij3-3,nij3
              inc = inc +1
              if( n < 1 .or. n > p3%nknij-4 ) cycle
              bij3(inc) = b_spl_rec(n,3,dij, p3%knij, p3%nknij)
              dbij3(inc) = db_spl(n,3,dij, p3%knij, p3%nknij)
            enddo
              
            do kk=1,lspr(0,ia)
              ka = lspr(kk,ia)
              if( jsp == ksp .and. ka <= ja ) cycle
              ks = int(tag(ka))
              if( ks /= ksp ) cycle
              xk(1:3) = ra(1:3,ka)
              xik(1:3) = xk(1:3) -xi(1:3)
              xjk(1:3) = xk(1:3) -xj(1:3)
              rik(1:3) = h(1:3,1)*xik(1) +h(1:3,2)*xik(2) +h(1:3,3)*xik(3)
              dik2 = rik(1)*rik(1) +rik(2)*rik(2) +rik(3)*rik(3)
              if( dik2 > p3%rcik2 ) cycle
              dik = sqrt(dik2)
              rjk(1:3) = h(1:3,1)*xjk(1) +h(1:3,2)*xjk(2) +h(1:3,3)*xjk(3)
              djk2 = rjk(1)*rjk(1) +rjk(2)*rjk(2) +rjk(3)*rjk(3)
              if( djk2 > p3%rcjk2 ) cycle
              djk = sqrt(djk2)
              nik3 = knot_index(dik, p3%nknik, p3%knik)
              njk3 = knot_index(djk, p3%nknjk, p3%knjk)
!.....B-spline part
              do nik=nik3-3,nik3
                if( nik < 1 .or. nik > p3%nknik-4 ) cycle
                bik3 = b_spl_rec(nik,3,dik, p3%knik, p3%nknik)
                dbik3 = db_spl(nik,3,dik, p3%knik, p3%nknik)
                do njk=njk3-3,njk3
                  if( njk < 1 .or. njk > p3%nknjk-4 ) cycle
                  bjk3 = b_spl_rec(njk,3,djk, p3%knjk, p3%nknjk)
                  dbjk3 = db_spl(njk,3,djk, p3%knjk, p3%nknjk)
                  l = 0
                  do nij=nij3-3,nij3
                    l = l +1
                    if( nij < 1 .or. nij > p3%nknij-4 ) cycle
!.....Energy
                    c3t = p3%coefs(nij,nik,njk)
                    tmp = c3t*bij3(l)*bik3*bjk3
                    epi(ia) = epi(ia) +tmp
                    epotl3 = epotl3 +tmp
!.....Force
                    tmpij(1:3) = dbij3(l)*bik3*bjk3*drijj(1:3)
                    tmpik(1:3) = bij3(l)*dbik3*bjk3*drikk(1:3)
                    tmpjk(1:3) = bij3(l)*bik3*dbjk3*drjkk(1:3)
                    do ixyz=1,3
!$omp atomic
                      aal3(ixyz,ia)= aal3(ixyz,ia) +c3t*(tmpij(ixyz) +tmpik(ixyz))
!$omp atomic
                      aal3(ixyz,ja)= aal3(ixyz,ja) +c3t*(-tmpij(ixyz) +tmpjk(ixyz))
!$omp atomic
                      aal3(ixyz,ka)= aal3(ixyz,ka) +c3t*(-tmpik(ixyz) -tmpjk(ixyz))
                    enddo
!.....Stresses
                    do jxyz=1,3
                      do ixyz=1,3
!$omp atomic
                        strsl(ixyz,jxyz,ia)= strsl(ixyz,jxyz,ia) &
                             -0.5d0 *c3t *(xij(jxyz)*tmpij(ixyz) +xik(jxyz)*tmpik(ixyz))
!$omp atomic
                        strsl(ixyz,jxyz,ja)= strsl(ixyz,jxyz,ja) &
                             -0.5d0 *c3t *(xij(jxyz)*tmpij(ixyz) +xjk(jxyz)*tmpjk(ixyz))
!$omp atomic
                        strsl(ixyz,jxyz,ka)= strsl(ixyz,jxyz,ka) &
                             -0.5d0 *c3t *(xik(jxyz)*tmpik(ixyz) +xjk(jxyz)*tmpjk(ixyz))
                      enddo
                    enddo
                  enddo ! nij
                enddo ! njk
              enddo ! nik

            enddo  ! kk
          enddo  ! jj
        enddo  ! ksp
      enddo  ! jsp

    enddo ! ia
!$omp end do
!$omp end parallel

    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,aal2,3)
    if( has_trios ) call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,aal3,3)
    aa(1:3,1:natm) = aa(1:3,1:natm) +aal2(1:3,1:natm) +aal3(1:3,1:natm)

    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,strsl,9)
    strs(1:3,1:3,1:natm) = strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)

!-----gather epot
    epot2 = 0d0
    epot3 = 0d0
    call mpi_allreduce(epotl2,epot2,1,mpi_real8,mpi_sum,mpi_world,ierr)
    if( has_trios ) call mpi_allreduce(epotl3,epot3,1,mpi_real8,mpi_sum,mpi_world,ierr)
    epot= epot +epot2 +epot3
    if( myid == 0 .and. iprint > 2 ) &
         print '(a,2es12.4)',' force_uf3 epot2,epot3 = ',epot2,epot3

    return
  end subroutine force_uf3_rec
!=======================================================================
  subroutine b_spl(r,ts,nmax,nr,b,db)
!
!  Non-recursive implementation of B-spline function at R.
!  Assuming maximum order (d) = 3.
!
!  Args:
!    r --- position to be evaluated
!    ts --- a list of knots {t_n}
!    nmax --- length of ts
!
!  Return:
!    nr: index n in ts of position r
!    b: B array
!    db: dB array (derivative of B)
!
    integer,intent(in):: nmax
    real(8),intent(in):: r,ts(nmax)
    integer,intent(out):: nr
    real(8),intent(out):: b(-3:0), db(-3:0)
!.....local variables
    real(8):: btmp(-3:+1,0:3)  ! Temporal B(n,d) array with n in (-3,+1)
    real(8):: dbtmp(-3:0)
    integer:: id, n, l
    real(8):: tn0,tn1,tn2,tn3,tn4,dt1,dt2,tmp1,tmp2
    real(8),parameter:: teps = 1d-8  ! epsilon for neighboring knot distance

!...index in the knot where ts(nr) <= r < ts(nr+1)
    nr = knot_index(r,nmax,ts)

!.....Compute B(n,d)
    btmp(:,:) = 0d0
    btmp(0,0) = 1d0
    do id = 1,3
      do l = -id,0
        n = nr +l
        tn0 = ts(n)
        tn1 = ts(n+id)
        dt1 = tn1 -tn0
        tmp1 = 0d0
        if( abs(dt1) > teps ) tmp1 = (r-tn0)/dt1 *btmp(l,id-1)
        tn2 = ts(n+1)
        tn3 = ts(n+id+1)
        dt2 = tn3 -tn2
        tmp2 = 0d0
        if( abs(dt2) > teps ) tmp2 = (tn3-r)/dt2 *btmp(l+1,id-1)
        btmp(l,id) = tmp1 + tmp2
      enddo
    enddo

!.....Compute dB(n) where d=3 is fixed
    dbtmp(:) = 0d0
    do l = -3,0
      n = nr +l
      tn0 = ts(n)
      tn1 = ts(n+1)
      tn3 = ts(n+3)
      tn4 = ts(n+1+3)
      tmp1 = 0d0
      if( abs(tn3-tn0) > teps ) tmp1 = btmp(l,2)/(tn3 -tn0)
      tmp2 = 0d0
      if( abs(tn4-tn1) > teps ) tmp2 = btmp(l+1,2)/(tn4 -tn1)
      dbtmp(l) = 3d0 *(tmp1 -tmp2)
    enddo

    b(-3:0) = btmp(-3:0,3)
    db(-3:0) = dbtmp(-3:0)

    return
  end subroutine b_spl
!=======================================================================
  recursive function b_spl_rec(n,d,r,ts,nmax) result(val)
!
!  TODO: check the efficiency of recursive func
!  
!  Recursive implementation of B-spline function with N and D as indices
!  and R as an argument.
!  TS --- a list of {t_n}
!
    integer,intent(in):: n,d,nmax
    real(8),intent(in):: r,ts(nmax)
    real(8):: val
    real(8):: denom1,denom2
    
    if( d == 0 ) then
      if( r >= ts(n) .and. r < ts(n+1) ) then
        val = 1d0
        return
      else
        val = 0d0
        return
      endif
    else
      val = 0d0
      denom1 = ts(n+d)-ts(n)
      if( abs(denom1) > 1d-8 ) then
        val = val +(r-ts(n))/denom1 *b_spl_rec(n,d-1,r,ts,nmax)
      endif
      denom2 = ts(n+d+1)-ts(n+1)
      if( abs(denom2) > 1d-8 ) then
        val = val +(ts(n+d+1)-r)/denom2 *b_spl_rec(n+1,d-1,r,ts,nmax)
      endif
    endif
    return
  end function b_spl_rec
!=======================================================================
  function db_spl(n,d,r,ts,nmax)
    integer,intent(in):: n,d,nmax
    real(8),intent(in):: r,ts(nmax)
    real(8):: db_spl
    real(8):: denom1, denom2

    db_spl = 0d0
    denom1 = ts(n+d)-ts(n)
    if( denom1 > 1d-8 ) then
      db_spl = db_spl +d*b_spl_rec(n,d-1,r,ts,nmax) /denom1
    endif
    denom2 = ts(n+d+1)-ts(n+1)
    if( denom2 > 1d-8 ) then
      db_spl = db_spl -d*b_spl_rec(n+1,d-1,r,ts,nmax)/denom2
    endif
    return
  end function db_spl
!=======================================================================
  function knot_index(r, nknot, knots) result(n)
    real(8),intent(in):: r
    integer,intent(in):: nknot
    real(8),intent(in):: knots(nknot)
    integer:: n, i

!.....TODO: use faster algorithm
    n = 0
    do i=1,nknot
      if( r > knots(i) ) then
        n = i
      else
        return
      endif
    enddo
    return
  end function knot_index
!=======================================================================
  subroutine set_paramsdir_uf3(dname)
!
!  Accessor routine to set paramsdir.
!
    implicit none
    character(len=*),intent(in):: dname

    paramsdir = trim(dname)
    return
  end subroutine set_paramsdir_uf3
!=======================================================================
  subroutine print_2b(ps)
    type(prm2),intent(in):: ps
    integer:: i
    character:: c*2
    print '(/,a)','   UF3 parameters of 2B for '&
         //trim(ps%csi)//'-'//trim(ps%csj)
    print '(5(a,1x))','     cb,csi,csj,cknot = ',ps%cb,ps%csi,ps%csj,ps%cknot
    print '(a,2i3)','     nklead,nktrail = ',ps%nklead,ps%nktrail
    print '(a,2i3)','     nknot,ncoef = ',ps%nknot, ps%ncoef
    print '(a,f6.3)','     rc = ',ps%rc
    write(c,'(i2)') size(ps%knots)
    print '(a,'//c//'(1x,f7.3))','     knots =',(ps%knots(i),i=1,size(ps%knots))
    write(c,'(i2)') size(ps%coefs)
    print '(a,'//trim(c)//'(1x,f7.3))','     coefs =',(ps%coefs(i),i=1,size(ps%coefs))
  end subroutine print_2b
!=======================================================================
  subroutine print_3b(ps)
    type(prm3),intent(in):: ps
    integer:: i
    character:: c*2
    print '(/,a)','   UF3 parameters of 3B for '&
         //trim(ps%csi)//'-'//trim(ps%csj)//'-'//trim(ps%csk)
    print '(6(a,1x))','     cb,csi,csj,csk,cknot = ',ps%cb,ps%csi,ps%csj,ps%csk,ps%cknot
    print '(a,2i3)','     nklead,nktrail = ',ps%nklead,ps%nktrail
    print '(a,6i3)','     nknij,nknik,nknjk,ncfij,ncfik,ncfjk = ', &
         ps%nknij, ps%nknik, ps%nknjk, ps%ncfij, ps%ncfik, ps%ncfjk
    print '(a,3f6.3)','     rcij,rcik,rcjk = ',ps%rcij,ps%rcik,ps%rcjk
    write(c,'(i2)') ps%nknij
    print '(a,'//c//'(1x,f7.3))','     knij =',(ps%knij(i),i=1,ps%nknij)
    write(c,'(i2)') ps%nknik
    print '(a,'//c//'(1x,f7.3))','     knik =',(ps%knik(i),i=1,ps%nknik)
    write(c,'(i2)') ps%nknjk
    print '(a,'//c//'(1x,f7.3))','     knjk =',(ps%knjk(i),i=1,ps%nknjk)
!!$    write(c,'(i2)') size(ps%coefs)
!!$    print '(a,'//trim(c)//'es11.2)','     coefs =',(ps%coefs(i),i=1,size(ps%coefs))
  end subroutine print_3b
end module UF3
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
