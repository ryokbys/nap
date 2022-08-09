module angular
!-----------------------------------------------------------------------
!                     Last modified: <2022-08-05 09:21:06 KOBAYASHI Ryo>
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

    if( size(aa3).ne.3*namax ) then
      call accum_mem('force_angular',-8*(size(aa3)+size(strsl)))
      deallocate(aa3,strsl)
      allocate(aa3(3,namax),strsl(3,3,namax))
      call accum_mem('force_angular',8*(size(aa3)+size(strsl)))
    endif

    strsl(1:3,1:3,1:natm+nb)= 0d0
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
          aa3(1:3,i)= aa3(1:3,i) +(tmpj(1:3)+tmpk(1:3))
!...Use omp atomic instead of reduction(+:aa3) for better parallel performace.
          do ixyz=1,3
!$omp atomic
            aa3(ixyz,j)= aa3(ixyz,j) -tmpj(ixyz)
          enddo
          do ixyz=1,3
!$omp atomic
            aa3(ixyz,k)= aa3(ixyz,k) -tmpk(ixyz)
          enddo
!.....Stress
          do jxyz=1,3
            strsl(1:3,jxyz,i)=strsl(1:3,jxyz,i) &
                 -0.5d0*xij(jxyz)*tmpj(1:3) & !*volj &
                 -0.5d0*xik(jxyz)*tmpk(1:3) !*volk
            do ixyz=1,3
!$omp atomic
              strsl(ixyz,jxyz,j)=strsl(ixyz,jxyz,j) &
                   -0.5d0*xij(jxyz)*tmpj(ixyz) !*volj
            enddo
            do ixyz=1,3
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
        read(cline,*) ctype
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
end module angular
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd lib"
!     End:
