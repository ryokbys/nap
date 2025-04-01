module angular
!-----------------------------------------------------------------------
!                     Last modified: <2025-03-29 00:03:31 KOBAYASHI Ryo>
!-----------------------------------------------------------------------
  use pmdvars,only: nspmax,nsp
  use util,only: csp2isp
  use memory,only: accum_mem
  include "./const.h"
  save
  
  integer,parameter:: ioprms = 50
  character(len=128):: paramsdir = '.'
  character(len=128),parameter:: paramsfname = 'in.params.angular'

!.....Small enough value for some criterion
  real(8),parameter:: eps = 1d-10
  integer,parameter:: msp = nspmax

  real(8):: rc3s(nspmax,nspmax,nspmax)
  real(8):: alps(nspmax,nspmax,nspmax)
  real(8):: bets(nspmax,nspmax,nspmax)
  real(8):: gmms(nspmax,nspmax,nspmax)
  real(8):: shfts(nspmax,nspmax,nspmax)
  logical:: interact3(nspmax,nspmax,nspmax)

!.....For fitpot
  logical:: lprmset_angular = .false.
  integer:: nprms
  real(8),allocatable:: params(:)
  real(8),allocatable:: ge_alp(:,:,:),ge_bet(:,:,:),ge_gmm(:,:,:)
  real(8),allocatable:: gf_alp(:,:,:,:,:),gf_bet(:,:,:,:,:),gf_gmm(:,:,:,:,:)
  real(8),allocatable:: gs_alp(:,:,:,:),gs_bet(:,:,:,:),gs_gmm(:,:,:,:)
  
contains
  subroutine force_angular(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_world,myid,epi,epot,nismax,specorder,lstrs,iprint,l1st)
!-----------------------------------------------------------------------
!  Parallel implementation of SW-like angular force-field.
!  Currently, only type-1 and -2 are available, and the difference
!  between them is only the shift of energy and thus no need to separate
!  the code for them.
!-----------------------------------------------------------------------
    implicit none
    include "mpif.h"
    include "./params_unit.h"
    integer,intent(in):: namax,natm,nnmax,nismax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_world,myid,lspr(0:nnmax,namax),nex(3)
    real(8),intent(in):: ra(3,namax),tag(namax) &
         ,h(3,3),hi(3,3),sv(3,6),rc
    real(8),intent(inout):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    character(len=3),intent(in):: specorder(msp)
    logical,intent(in):: lstrs, l1st

!-----local
    integer:: i,j,k,l,m,n,ixyz,jxyz,is,js,ks,ierr
    real(8):: rij,rik,riji,riki,rij2,rik2,rc2,rc3,alp,bet,gmm,shft
    real(8):: tmp,tmpj(3),tmpk(3),vexp,df2,csn,tcsn,tcsn2,dhrij,dhrik &
         ,dhcsn,vol,voli,volj,volk,drijj(3),rcmax
    real(8):: drikk(3),dcsni(3),dcsnj(3),dcsnk(3),drijc,drikc,x,y,z,bl &
         ,xi(3),xj(3),xk(3),xij(3),xik(3),at(3)
    real(8):: epotl,epotl2,epotl3,epott
    real(8),save:: rcmax2
    real(8),allocatable,save:: aa3(:,:)
    real(8),allocatable,save:: strsl(:,:,:)

!-----only at 1st call
    if( l1st ) then
!!$      call read_params_angular(myid,mpi_world,iprint,specorder)
      if( allocated(aa3) ) then
        call accum_mem('force_angular',-8*(size(aa3)+size(strsl)))
        deallocate(aa3,strsl)
      endif
      allocate(aa3(3,namax),strsl(3,3,namax))
      call accum_mem('force_angular',8*(size(aa3)+size(strsl)))
!-------check rc
      rcmax = 0d0
      do is=1,nspmax
        do js=1,nspmax
          do ks=1,nspmax
            rcmax = max(rcmax,rc3s(is,js,ks))
          enddo
        enddo
      enddo
      rcmax2= rcmax*rcmax
      if( myid.eq.0 .and. iprint.ge.ipl_basic ) then
        print *,''
        print *,'Angular potential:'
        write(6,'(a,es12.4)') '   rc from in.pmd           =',rc
        write(6,'(a,es12.4)') '   rcmax for this potential =',rcmax
      endif
      if( rc .lt. rcmax ) then
!!$      if( int(rc*100d0) &
!!$           .ne.int(swrc*swl*100d0) ) then
        if( myid.eq.0 ) then
          write(6,'(1x,a)') "ERROR: Cutoff radius is not appropriate !!!"
          write(6,'(1x,a,es12.4)') "  rc should be longer than ", rcmax
        endif
        call mpi_finalize(ierr)
        stop
      endif
    endif

    if( size(aa3).lt.3*namax ) then
      call accum_mem('force_angular',-8*(size(aa3)+size(strsl)))
      deallocate(aa3,strsl)
      allocate(aa3(3,namax),strsl(3,3,namax))
      call accum_mem('force_angular',8*(size(aa3)+size(strsl)))
    endif

    strsl(1:3,1:3,1:namax)= 0d0
    epotl3= 0d0
    aa3(1:3,1:namax)= 0d0
!.....Loop over i
!$omp parallel
!$omp do reduction(+:epotl3) &
!$omp    private(i,xi,is,n,j,js,xj,x,y,z,xij,rij2,rij,riji,drijj,m,k, &
!$omp            ks,rc3,xk,xik,rik2,rik,riki,drijc,drikc,alp,bet,gmm, &
!$omp            shft,csn,tcsn,tcsn2,vexp,tmp,dhrij,dhrik,dhcsn,drikk, &
!$omp            dcsnj,dcsnk,dcsni,tmpj,tmpk,ixyz)
    do i=1,natm
      xi(1:3)=ra(1:3,i)
      is= int(tag(i))
!.....Loop over j
      do n=1,lspr(0,i)-1
        j=lspr(n,i)
        js= int(tag(j))
        xj(1:3)= ra(1:3,j)
        x = xj(1) -xi(1)
        y = xj(2) -xi(2)
        z = xj(3) -xi(3)
        xij(1:3)= (h(1:3,1)*x +h(1:3,2)*y +h(1:3,3)*z)
        rij2 = xij(1)*xij(1) +xij(2)*xij(2) +xij(3)*xij(3)
        if( rij2.gt.rcmax2 ) cycle
        rij = dsqrt(rij2)
        riji= 1d0/rij
        drijj(1:3)= xij(1:3)*riji
!.....Loop over k
        do m=n+1,lspr(0,i)
          k=lspr(m,i)
          ks= int(tag(k))
          if( .not. interact3(is,js,ks) ) cycle
          rc3 = rc3s(is,js,ks)
          if( rij.ge.rc3 ) cycle
          xk(1:3)= ra(1:3,k)
          x = xk(1) -xi(1)
          y = xk(2) -xi(2)
          z = xk(3) -xi(3)
          xik(1:3)= (h(1:3,1)*x +h(1:3,2)*y +h(1:3,3)*z)
          rik2 = xik(1)*xik(1) +xik(2)*xik(2) +xik(3)*xik(3)
          if( rik2.ge.rc3*rc3 ) cycle
          rik = dsqrt(rik2)
          riki= 1d0/rik
          drijc= 1d0/(rij-rc3)
          drikc= 1d0/(rik-rc3)
!.....Parameters
          alp = alps(is,js,ks)
          bet = bets(is,js,ks)
          gmm = gmms(is,js,ks)
          shft= shfts(is,js,ks)
!.....Common terms
          csn=(xij(1)*xik(1) +xij(2)*xik(2) +xij(3)*xik(3)) *(riji*riki)
          tcsn = csn -gmm
          tcsn2= (tcsn*tcsn +shft)
          vexp= dexp(bet*drijc +bet*drikc)
!.....Potential
          tmp= alp *vexp *tcsn2
          epi(i)= epi(i) +tmp
          epotl3= epotl3 +tmp
!.....Force
          dhrij= -drijc*drijc *bet *tmp
          dhrik= -drikc*drikc *bet *tmp
          dhcsn= 2d0 *alp *vexp *tcsn 
          drikk(1:3)= xik(1:3)*riki
          dcsnj(1:3)= (-xij(1:3)*csn*(riji*riji) +xik(1:3)*(riji*riki))
          dcsnk(1:3)= (-xik(1:3)*csn*(riki*riki) +xij(1:3)*(riji*riki))
          dcsni(1:3)= -dcsnj(1:3) -dcsnk(1:3)
          tmpj(1:3)= dhcsn*dcsnj(1:3) +dhrij*drijj(1:3)
          tmpk(1:3)= dhcsn*dcsnk(1:3) +dhrik*drikk(1:3)
!...Use omp atomic instead of reduction(+:aa3) for better parallel performace.
          do ixyz=1,3
!$omp atomic
            aa3(ixyz,i)= aa3(ixyz,i) +(tmpj(ixyz)+tmpk(ixyz))
!$omp atomic
            aa3(ixyz,j)= aa3(ixyz,j) -tmpj(ixyz)
!$omp atomic
            aa3(ixyz,k)= aa3(ixyz,k) -tmpk(ixyz)
          enddo
!.....Stress
          do jxyz=1,3
            do ixyz=1,3
!$omp atomic
              strsl(ixyz,jxyz,i)= strsl(ixyz,jxyz,i) + &
                   (-0.5d0*xij(jxyz)*tmpj(ixyz) &
                   -0.5d0*xik(jxyz)*tmpk(ixyz))
!$omp atomic
              strsl(ixyz,jxyz,j)=strsl(ixyz,jxyz,j) &
                   -0.5d0*xij(jxyz)*tmpj(ixyz) !*volj
!$omp atomic
              strsl(ixyz,jxyz,k)=strsl(ixyz,jxyz,k) &
                   -0.5d0*xik(jxyz)*tmpk(ixyz) !*volk
            enddo
          enddo ! jxyz

        enddo
      enddo
    enddo
!$omp end do
!$omp end parallel

    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,aa3,3)
    aa(1:3,1:natm)= aa(1:3,1:natm) +aa3(1:3,1:natm)

    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,strsl,9)
    strs(1:3,1:3,1:natm) = strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)

!.....Gather epot
    epott= 0d0
    epotl= epotl3
    call mpi_allreduce(epotl,epott,1,mpi_real8,mpi_sum,mpi_world,ierr)
    epot= epot +epott
    if( iprint.ge.ipl_info ) print *,'epot angular = ',epott

    return
  end subroutine force_angular
!=======================================================================
  subroutine read_params_angular(myid,mpi_world,iprint,specorder)
    use util, only: num_data, is_numeric
    implicit none
    include 'mpif.h'
    integer,intent(in):: myid,mpi_world,iprint
    character(len=3),intent(in):: specorder(msp)

    integer:: itmp,ierr,isp,jsp,ksp,nd,itype
    real(8):: alp,bet,gmm,shft,rc3
    logical:: lexist
    character(len=128):: cfname,ctmp,cline,ctype
    character(len=3):: cspi,cspj,cspk

!.....read parameters at the 1st call
    if( myid.eq.0 ) then
      interact3(:,:,:) = .false.
!.....Initialize parameters
      rc3s(:,:,:) = 0d0
      gmms(:,:,:) = -1d0/3
      alps(:,:,:) = 1d0
      bets(:,:,:) = 1d0
      shfts(:,:,:) = 0d0
!.....Check whether the file exists      
      cfname = trim(paramsdir)//'/'//trim(paramsfname)
      inquire(file=cfname,exist=lexist)
      if( .not. lexist ) then
        if( iprint.ge.ipl_warn ) then
          write(6,'(a)') ' WARNING: in.params.angular does not exist !!!.'
          write(6,'(a)') '          Default parameters will be used.'
        endif
        goto 20
      endif
      
!.....Read file if exists
      if( iprint.ge.ipl_basic ) write(6,'(/,a)') ' Angular parameters read from file:'
      open(ioprms,file=cfname,status='old')
      do while(.true.)
        read(ioprms,'(a)',end=10) cline
        nd = num_data(cline,' ')
        if( nd.eq.0 ) cycle
        if( cline(1:1).eq.'!' .or. cline(1:1).eq.'#' ) cycle
        read(cline,*) ctype  ! read angular-type
        if( is_numeric(ctype) ) then
          read(ctype,*) itype
        else
          itype = itype_from(ctype)
        endif
        if( itype.eq.1 ) then
          if( nd.ne.8 .and. iprint.ge.ipl_basic ) then
            print *,'WARNING@read_params_angular: # of entry wrong for type-1,' &
                 //' so skip the line.'
            print *,'   '//trim(cline)
          endif
          read(cline,*) ctmp, cspi,cspj,cspk,rc3,alp,bet,gmm
          isp = csp2isp(cspi)
          jsp = csp2isp(cspj)
          ksp = csp2isp(cspk)
          if( iprint.ge.ipl_basic ) print '(a,3(a3,1x),4f7.3)', &
               '   cspi,cspj,cspk,rc3,alp,bet,gmm=', &
               trim(cspi),trim(cspj),trim(cspk),rc3,alp,bet,gmm
          if( isp.gt.0 .and. jsp.gt.0 .and. ksp.gt.0 ) then
            interact3(isp,jsp,ksp) = .true.
            rc3s(isp,jsp,ksp) = rc3
            alps(isp,jsp,ksp) = alp
            bets(isp,jsp,ksp) = bet
            gmms(isp,jsp,ksp) = gmm
!.....Symmetrize
            interact3(isp,ksp,jsp) = .true.
            rc3s(isp,ksp,jsp) = rc3
            alps(isp,ksp,jsp) = alp
            bets(isp,ksp,jsp) = bet
            gmms(isp,ksp,jsp) = gmm
          else
            if( iprint.ge.ipl_info ) then
              print *,' Angular parameter read but not used: cspi,cspj,cspk='&
                   ,cspi,cspj,cspk
            endif
          endif

!.....Read type-2 parameters
        else if( itype.eq.2 ) then
          if( nd.ne.9 .and. iprint.ge.ipl_basic ) then
            print *,'WARNING@read_params_angular: # of entry wrong for type-1,' &
                 //' so skip the line.'
            print *,'   '//trim(cline)
          endif
          read(cline,*) ctmp, cspi,cspj,cspk,rc3,alp,bet,gmm,shft
          isp = csp2isp(cspi)
          jsp = csp2isp(cspj)
          ksp = csp2isp(cspk)
          if( iprint.ge.ipl_basic ) print '(a,3(a3,1x),5f7.3)', &
               '   cspi,cspj,cspk,rc3,alp,bet,gmm,shft=', &
               trim(cspi),trim(cspj),trim(cspk),rc3,alp,bet,gmm,shft
          if( isp.gt.0 .and. jsp.gt.0 .and. ksp.gt.0 ) then
            interact3(isp,jsp,ksp) = .true.
            rc3s(isp,jsp,ksp) = rc3
            alps(isp,jsp,ksp) = alp
            bets(isp,jsp,ksp) = bet
            gmms(isp,jsp,ksp) = gmm
            shfts(isp,jsp,ksp) = shft
!.....Symmetrize
            interact3(isp,ksp,jsp) = .true.
            rc3s(isp,ksp,jsp) = rc3
            alps(isp,ksp,jsp) = alp
            bets(isp,ksp,jsp) = bet
            gmms(isp,ksp,jsp) = gmm
            shfts(isp,ksp,jsp) = shft
          else
            if( iprint.ge.ipl_info ) then
              print *,' Angular parameter read but not used: cspi,cspj,cspk='&
                   ,cspi,cspj,cspk
            endif
          endif
          
        endif  ! itype
      enddo
10    close(ioprms)
      
    endif

20  continue

    call mpi_bcast(interact3,nspmax*nspmax*nspmax,mpi_logical,0,mpi_world,ierr)
    call mpi_bcast(rc3s,nspmax*nspmax*nspmax,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(alps,nspmax*nspmax*nspmax,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(bets,nspmax*nspmax*nspmax,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(gmms,nspmax*nspmax*nspmax,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(shfts,nspmax*nspmax*nspmax,mpi_real8,0,mpi_world,ierr)

    return
  end subroutine read_params_angular
!=======================================================================
  subroutine set_paramsdir_angular(dname)
!
!  Accessor routine to set paramsdir.
!
    implicit none
    character(len=*),intent(in):: dname

    paramsdir = trim(dname)
    return
  end subroutine set_paramsdir_angular
!=======================================================================
  subroutine set_params_angular(ndimp,params_in,ctype,rc3,interact3_in)
!
!  Accessor routine to set angular parameters from outside.
!  Curretnly this routine is supposed to be called only on serial run.
!
    integer,intent(in):: ndimp
    real(8),intent(in):: params_in(ndimp),rc3
    character(len=*),intent(in):: ctype
    logical,intent(in):: interact3_in(nspmax,nspmax,nspmax)

    integer:: i,j,k,inc,itmp,nspt,nint
    real(8):: alp,bet,gmm

    nprms = ndimp
    if( .not.allocated(params) ) allocate(params(nprms))
    params(1:nprms) = params_in(1:ndimp)
    lprmset_angular = .true.

    alps(:,:,:) = 0d0
    bets(:,:,:) = 0d0
    gmms(:,:,:) = 0d0
    rc3s(:,:,:) = 0d0
    interact3(:,:,:) = interact3_in(:,:,:)
    nint = 0
    do i=1,nspmax
      do j=1,nspmax
        do k=j,nspmax
          if( .not.interact3(i,j,k) ) cycle
          rc3s(i,j,k) = rc3
          nint = nint +1
        enddo
      enddo
    enddo

    if( index(ctype,'angular1').ne.0 ) then
      if( nint*3.ne.nprms ) then
        print *,'ERROR@set_params_angular: nint*3.ne.nprms !!!'
        print *,'  nint, nprms = ', nint, nprms
        print *,'Probably you need to set interactions3 correctly in in.fitpot.'
        stop 1
      endif

      inc = 0
      do i=1,nspmax
        do j=1,nspmax
          do k=j,nspmax
            if( .not.interact3(i,j,k) ) cycle
            inc = inc +1
            alps(i,j,k) = params(inc)
            alps(i,k,j) = alps(i,j,k)
            inc = inc +1
            bets(i,j,k) = params(inc)
            bets(i,k,j) = bets(i,j,k)
            inc = inc +1
            gmms(i,j,k) = params(inc)
            gmms(i,k,j) = gmms(i,j,k)
!!$            print *,'i,j,k,alp,bet,gmm=',i,j,k,alps(i,j,k),bets(i,j,k),gmms(i,j,k)
          enddo
        enddo
      enddo

    endif
    

    return
  end subroutine set_params_angular
!=======================================================================
  function itype_from(str)
!
!  Convert string type such as angular1 to integer type.
!
    character(len=*),intent(in):: str
    integer:: itype_from

    itype_from = 0
    if( trim(str).eq.'angular1' ) then
      itype_from = 1
    else if( trim(str).eq.'angular2' ) then
      itype_from = 2
    endif
  end function itype_from
!=======================================================================
  subroutine gradw_angular(namax,natm,tag,ra,nnmax,h,rcin,lspr, &
       iprint,ndimp,gwe,gwf,gws,lematch,lfmatch,lsmatch,iprm0, &
       nfcal,lfrc_eval)
!
!  Gradient of angular wrt parameters.
!  Note: This routine is always called in single run,
!  thus no need of parallel implementation.
!
    use util, only: itotOf
    implicit none
    integer,intent(in):: namax,natm,nnmax,iprint,iprm0
    integer,intent(in):: lspr(0:nnmax,namax)
    real(8),intent(in):: ra(3,namax),tag(namax),h(3,3)
    real(8),intent(inout):: rcin
    integer,intent(in):: ndimp
    integer,intent(in):: nfcal
    logical,intent(in):: lfrc_eval(natm)
    real(8),intent(inout):: gwe(ndimp),gwf(3,ndimp,nfcal),gws(6,ndimp)
    logical,intent(in):: lematch,lfmatch,lsmatch

!.....local
    integer:: i,j,k,l,m,n,ixyz,jxyz,is,js,ks,ierr,ia,ka,ja,jra,kra
    integer:: ifcal,jfcal,kfcal,ip
    real(8):: rij,rik,riji,riki,rij2,rik2,rc2,rc3,alp,bet,gmm,shft
    real(8):: tmp,tmpj(3),tmpk(3),vexp,df2,csn,tcsn,tcsn2,dhrij,dhrik &
         ,dhcsn,vol,voli,volj,volk,drijj(3),rcmax,rcmax2
    real(8):: drikk(3),dcsni(3),dcsnj(3),dcsnk(3),drijc,drikc,x,y,z,bl &
         ,xi(3),xj(3),xk(3),xij(3),xik(3),at(3)
    real(8):: tmpja(3),tmpjb(3),tmpjg(3),tmpka(3),tmpkb(3),tmpkg(3)
    real(8),save:: rcin2 = -1d0
    integer,save,allocatable:: ia2ifcal(:)

    if( .not.allocated(ia2ifcal) ) allocate(ia2ifcal(namax))
    if( rcin2 < 0d0 ) then
      if( rcin < rcmax ) then
        call mpi_finalize(ierr)
        stop
      endif
      rcin2 = rcin*rcin
    endif
    rcmax2 = rcmax*rcmax

    if( lematch .and. .not.allocated(ge_alp) ) then
      allocate(ge_alp(nspmax,nspmax,nspmax),ge_bet(nspmax,nspmax,nspmax),&
           ge_gmm(nspmax,nspmax,nspmax))
      call accum_mem('force_angular', &
           8*(size(ge_alp)+size(ge_bet)+size(ge_gmm)))
    endif

    if( lfmatch .and. (.not.allocated(gf_alp) &
         .or. size(gf_alp).ne.3*natm *nspmax**3) ) then
      if( allocated(gf_alp) ) then
        call accum_mem('force_angular', &
             -8*(size(gf_alp)+size(gf_bet)+size(gf_gmm)))
        deallocate(gf_alp, gf_bet, gf_gmm)
      endif
      allocate(gf_alp(3,nspmax,nspmax,nspmax,natm), &
           gf_bet(3,nspmax,nspmax,nspmax,natm), &
           gf_gmm(3,nspmax,nspmax,nspmax,natm) )
      call accum_mem('force_angular', &
           8*(size(gf_alp)+size(gf_bet)+size(gf_gmm)))
    endif

    if( lsmatch .and. .not.allocated(gs_alp) ) then
      allocate(gs_alp(6,nspmax,nspmax,nspmax), &
           gs_bet(6,nspmax,nspmax,nspmax), &
           gs_gmm(6,nspmax,nspmax,nspmax))
      call accum_mem('force_angular', &
           8*(size(gs_alp)+size(gs_bet)+size(gs_gmm)))
    endif

!.....Set nsp by max isp of atoms in the system
    nsp = 0
    do i=1,natm
      nsp = max(nsp,int(tag(i)))
    enddo

    if( lematch ) then
      ge_alp(:,:,:) = 0d0
      ge_bet(:,:,:) = 0d0
      ge_gmm(:,:,:) = 0d0
    endif
    if( lfmatch ) then
      gf_alp(:,:,:,:,:) = 0d0
      gf_bet(:,:,:,:,:) = 0d0
      gf_gmm(:,:,:,:,:) = 0d0
    endif
    if( lsmatch ) then
      gs_alp(:,:,:,:) = 0d0
      gs_bet(:,:,:,:) = 0d0
      gs_gmm(:,:,:,:) = 0d0
    endif

!.....gradw_uf3 hereafter.....

    if( size(ia2ifcal) < namax ) then
      deallocate(ia2ifcal)
      allocate(ia2ifcal(namax))
    endif
!.....Construct ia2ifcal
    ifcal = 0
    ia2ifcal(:) = 0
    do ia=1,natm
      if( lfrc_eval(ia) ) then
        ifcal = ifcal +1
        ia2ifcal(ia) = ifcal
      endif
    enddo

    do i=1,natm
      xi(1:3)=ra(1:3,i)
      is= int(tag(i))
!.....Loop over j
      do n=1,lspr(0,i)-1
        j=lspr(n,i)
        js= int(tag(j))
        xj(1:3)= ra(1:3,j)
        x = xj(1) -xi(1)
        y = xj(2) -xi(2)
        z = xj(3) -xi(3)
        xij(1:3)= (h(1:3,1)*x +h(1:3,2)*y +h(1:3,3)*z)
        rij2 = xij(1)*xij(1) +xij(2)*xij(2) +xij(3)*xij(3)
        if( rij2.gt.rcmax2 ) cycle
        jra = itotOf(tag(j))
        rij = dsqrt(rij2)
        riji= 1d0/rij
        drijj(1:3)= xij(1:3)*riji
!.....Loop over k
        do m=n+1,lspr(0,i)
          k=lspr(m,i)
          ks= int(tag(k))
          if( .not. interact3(is,js,ks) ) cycle
          rc3 = rc3s(is,js,ks)
          if( rij.ge.rc3 ) cycle
          xk(1:3)= ra(1:3,k)
          x = xk(1) -xi(1)
          y = xk(2) -xi(2)
          z = xk(3) -xi(3)
          xik(1:3)= (h(1:3,1)*x +h(1:3,2)*y +h(1:3,3)*z)
          rik2 = xik(1)*xik(1) +xik(2)*xik(2) +xik(3)*xik(3)
          if( rik2.ge.rc3*rc3 ) cycle
          kra = itotOf(tag(k))
          rik = dsqrt(rik2)
          riki= 1d0/rik
          drijc= 1d0/(rij-rc3)
          drikc= 1d0/(rik-rc3)
!.....Parameters
          alp = alps(is,js,ks)
          bet = bets(is,js,ks)
          gmm = gmms(is,js,ks)
          shft= shfts(is,js,ks)
!.....Common terms
          csn=(xij(1)*xik(1) +xij(2)*xik(2) +xij(3)*xik(3)) *(riji*riki)
          tcsn = csn -gmm
          tcsn2= (tcsn*tcsn +shft)
          vexp= dexp(bet*drijc +bet*drikc)
!.....Potential
          tmp= alp *vexp *tcsn2
!!$          epi(i)= epi(i) +tmp
!!$          epotl3= epotl3 +tmp
          if( lematch ) then
            ge_alp(is,js,ks)= ge_alp(is,js,ks) +vexp*tcsn2
            ge_bet(is,js,ks)= ge_bet(is,js,ks) +(drijc+drikc)*tmp
            ge_gmm(is,js,ks)= ge_gmm(is,js,ks) -2d0*alp*vexp*tcsn
          endif
!.....Force
          if( .not.lfmatch .and. .not.lsmatch ) cycle
          dhrij= -drijc*drijc *bet *tmp
          dhrik= -drikc*drikc *bet *tmp
          dhcsn= 2d0 *alp *vexp *tcsn 
          drikk(1:3)= xik(1:3)*riki
          dcsnj(1:3)= (-xij(1:3)*csn*(riji*riji) +xik(1:3)*(riji*riki))
          dcsnk(1:3)= (-xik(1:3)*csn*(riki*riki) +xij(1:3)*(riji*riki))
          dcsni(1:3)= -dcsnj(1:3) -dcsnk(1:3)
!!$          tmpj(1:3)= dhcsn*dcsnj(1:3) +dhrij*drijj(1:3)
!!$          tmpk(1:3)= dhcsn*dcsnk(1:3) +dhrik*drikk(1:3)
!.....w.r.t. alp
          tmpja(:)= drijj(:)*(-bet)*drijc**2*vexp*tcsn2 +dcsnj(:)*2d0*vexp*tcsn
          tmpka(:)= drikk(:)*(-bet)*drikc**2*vexp*tcsn2 +dcsnk(:)*2d0*vexp*tcsn
!.....w.r.t. bet
          tmpjb(:)= drijj(:)*(-alp)*drijc**2*vexp*tcsn2*(1d0+bet*drijc+bet*drikc) &
               +dcsnj(:)*2d0*alp*(drijc+drikc)*vexp*tcsn
          tmpkb(:)= drikk(:)*(-alp)*drikc**2*vexp*tcsn2*(1d0+bet*drijc+bet*drikc) &
               +dcsnk(:)*2d0*alp*(drijc+drikc)*vexp*tcsn
!.....w.r.t. gmm
          tmpjg(:)= drijj(:)*2d0*alp*bet*drijc**2*vexp*tcsn -dcsnj(:)*2d0*alp*vexp
          tmpkg(:)= drikk(:)*2d0*alp*bet*drikc**2*vexp*tcsn -dcsnk(:)*2d0*alp*vexp
          if( lfmatch ) then
            ifcal = ia2ifcal(ia)
            jfcal = ia2ifcal(jra)
            kfcal = ia2ifcal(kra)
            if( ifcal.ne.0 ) then
              gf_alp(:,is,js,ks,ifcal)= gf_alp(:,is,js,ks,ifcal) +tmpja(:) +tmpka(:)
              gf_bet(:,is,js,ks,ifcal)= gf_bet(:,is,js,ks,ifcal) +tmpjb(:) +tmpkb(:)
              gf_gmm(:,is,js,ks,ifcal)= gf_gmm(:,is,js,ks,ifcal) +tmpjg(:) +tmpkg(:)
            endif
            if( jfcal.ne.0 ) then
              gf_alp(:,is,js,ks,jfcal)= gf_alp(:,is,js,ks,jfcal) -tmpja(:)
              gf_bet(:,is,js,ks,jfcal)= gf_bet(:,is,js,ks,jfcal) -tmpjb(:)
              gf_gmm(:,is,js,ks,jfcal)= gf_gmm(:,is,js,ks,jfcal) -tmpjg(:)
            endif
            if( kfcal.ne.0 ) then
              gf_alp(:,is,js,ks,kfcal)= gf_alp(:,is,js,ks,kfcal) -tmpka(:)
              gf_bet(:,is,js,ks,kfcal)= gf_bet(:,is,js,ks,kfcal) -tmpkb(:)
              gf_gmm(:,is,js,ks,kfcal)= gf_gmm(:,is,js,ks,kfcal) -tmpkg(:)
            endif
          endif
!!$!...Use omp atomic instead of reduction(+:aa3) for better parallel performace.
!!$          do ixyz=1,3
!!$            aa3(ixyz,i)= aa3(ixyz,i) +(tmpj(ixyz)+tmpk(ixyz))
!!$            aa3(ixyz,j)= aa3(ixyz,j) -tmpj(ixyz)
!!$            aa3(ixyz,k)= aa3(ixyz,k) -tmpk(ixyz)
!!$          enddo
!.....Stress
          if( lsmatch ) then
            gs_alp(1,is,js,ks)= gs_alp(1,is,js,ks) +xij(1)*tmpja(1) +xik(1)*tmpka(1)
            gs_alp(2,is,js,ks)= gs_alp(2,is,js,ks) +xij(2)*tmpja(2) +xik(2)*tmpka(2)
            gs_alp(3,is,js,ks)= gs_alp(3,is,js,ks) +xij(3)*tmpja(3) +xik(3)*tmpka(3)
            gs_alp(4,is,js,ks)= gs_alp(4,is,js,ks) +xij(2)*tmpja(3) +xik(2)*tmpka(3)
            gs_alp(5,is,js,ks)= gs_alp(5,is,js,ks) +xij(1)*tmpja(3) +xik(1)*tmpka(3)
            gs_alp(6,is,js,ks)= gs_alp(6,is,js,ks) +xij(1)*tmpja(2) +xik(1)*tmpka(2)

            gs_bet(1,is,js,ks)= gs_bet(1,is,js,ks) +xij(1)*tmpjb(1) +xik(1)*tmpkb(1)
            gs_bet(2,is,js,ks)= gs_bet(2,is,js,ks) +xij(2)*tmpjb(2) +xik(2)*tmpkb(2)
            gs_bet(3,is,js,ks)= gs_bet(3,is,js,ks) +xij(3)*tmpjb(3) +xik(3)*tmpkb(3)
            gs_bet(4,is,js,ks)= gs_bet(4,is,js,ks) +xij(2)*tmpjb(3) +xik(2)*tmpkb(3)
            gs_bet(5,is,js,ks)= gs_bet(5,is,js,ks) +xij(1)*tmpjb(3) +xik(1)*tmpkb(3)
            gs_bet(6,is,js,ks)= gs_bet(6,is,js,ks) +xij(1)*tmpjb(2) +xik(1)*tmpkb(2)

            gs_gmm(1,is,js,ks)= gs_gmm(1,is,js,ks) +xij(1)*tmpjg(1) +xik(1)*tmpkg(1)
            gs_gmm(2,is,js,ks)= gs_gmm(2,is,js,ks) +xij(2)*tmpjg(2) +xik(2)*tmpkg(2)
            gs_gmm(3,is,js,ks)= gs_gmm(3,is,js,ks) +xij(3)*tmpjg(3) +xik(3)*tmpkg(3)
            gs_gmm(4,is,js,ks)= gs_gmm(4,is,js,ks) +xij(2)*tmpjg(3) +xik(2)*tmpkg(3)
            gs_gmm(5,is,js,ks)= gs_gmm(5,is,js,ks) +xij(1)*tmpjg(3) +xik(1)*tmpkg(3)
            gs_gmm(6,is,js,ks)= gs_gmm(6,is,js,ks) +xij(1)*tmpjg(2) +xik(1)*tmpkg(2)
          endif
!!$          do jxyz=1,3
!!$            do ixyz=1,3
!!$              strsl(ixyz,jxyz,i)= strsl(ixyz,jxyz,i) + &
!!$                   (-0.5d0*xij(jxyz)*tmpj(ixyz) &
!!$                   -0.5d0*xik(jxyz)*tmpk(ixyz))
!!$              strsl(ixyz,jxyz,j)=strsl(ixyz,jxyz,j) &
!!$                   -0.5d0*xij(jxyz)*tmpj(ixyz) !*volj
!!$              strsl(ixyz,jxyz,k)=strsl(ixyz,jxyz,k) &
!!$                   -0.5d0*xik(jxyz)*tmpk(ixyz) !*volk
!!$            enddo
!!$          enddo ! jxyz

        enddo
      enddo
    enddo


!.....Tidy up gradient arrays
    if( lematch ) then  ! energy matching
      ip = iprm0
      do is=1,nspmax
        do js=1,nspmax
          do ks=js,nspmax
            if( .not.interact3(is,js,ks) ) cycle
            ip = ip +1
            gwe(ip)= gwe(ip) +ge_alp(is,js,ks)
            ip = ip +1
            gwe(ip)= gwe(ip) +ge_bet(is,js,ks)
            ip = ip +1
            gwe(ip)= gwe(ip) +ge_gmm(is,js,ks)
          enddo
        enddo
      enddo
    endif

    if( lfmatch ) then  ! force matching
      do ia=1,natm
        if( .not. lfrc_eval(ia) ) cycle
        ifcal = ia2ifcal(ia)
        is = int(tag(ia))
        ip = iprm0
        do is=1,nspmax
          do js=1,nspmax
            do ks=js,nspmax
              if( .not.interact3(is,js,ks) ) cycle
              ip = ip +1
              gwf(1:3,ip,ifcal)= gwf(1:3,ip,ifcal) +gf_alp(1:3,is,js,ks,ifcal)
              ip = ip +1
              gwf(1:3,ip,ifcal)= gwf(1:3,ip,ifcal) +gf_bet(1:3,is,js,ks,ifcal)
              ip = ip +1
              gwf(1:3,ip,ifcal)= gwf(1:3,ip,ifcal) +gf_gmm(1:3,is,js,ks,ifcal)
            enddo
          enddo
        enddo
      enddo ! ia
    endif

    if( lsmatch ) then  ! stress matching
      ip = iprm0
      do is=1,nspmax
        do js=1,nspmax
          do ks=js,nspmax
            if( .not.interact3(is,js,ks) ) cycle
            ip = ip +1
            gws(1:6,ip)= gws(1:6,ip) +gs_alp(1:6,is,js,ks)
            ip = ip +1
            gws(1:6,ip)= gws(1:6,ip) +gs_bet(1:6,is,js,ks)
            ip = ip +1
            gws(1:6,ip)= gws(1:6,ip) +gs_gmm(1:6,is,js,ks)
          enddo
        enddo
      enddo
    endif

    return
  end subroutine gradw_angular
!=======================================================================
  
end module angular
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd lib"
!     End:
