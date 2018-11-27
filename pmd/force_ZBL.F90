module ZBL
!-----------------------------------------------------------------------
!                     Last modified: <2018-11-27 16:52:23 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
!  Parallel implementation of ZBL repulsive potential with switching
!  function zeta(x).
!  See:
!   [1] Ziegler, Biersack, Littmark (ZBL) The Stopping and Range of Ions
!       in Matter (SRIM), 1985.
!   [2] M.Z. Hossain, J.B. Freund, and H.T. Johnson, Nuclear Inst.
!       and Methods in Physics Research, B 267, 1061 (2009).
!   [3] G. Bonny et al., JAP 121, 165107 (2017) for switching function.
!-----------------------------------------------------------------------
  implicit none
  save

  character(len=128):: paramsdir = '.'
  character(len=128),parameter:: paramsfname = 'in.params.ZBL'
  integer,parameter:: ioprms = 20
  
!.....Max num of species
  integer,parameter:: msp = 9

!.....Coulomb's constant, acc = 1.0/(4*pi*epsilon0)
  real(8),parameter:: acc  = 14.3998554737d0
!.....permittivity of vacuum
  real(8),parameter:: eps0 = 0.00552634939836d0  ! e^2 /Ang /eV

  logical:: interact(1:msp,1:msp) = .false.
  real(8):: qnucl(msp)
  real(8):: zbl_rc
  real(8):: r_inner(1:msp) = -1d0
  real(8):: r_outer(1:msp) = -1d0

!.....ZBL parameters
  real(8):: zbl_aa = 0.4683766d0
  real(8):: zbl_gamma = 0.23d0
  real(8):: zbl_alpha(1:4) = (/ &
       0.18180d0, &
       0.50990d0, &
       0.28020d0, &
       0.02817d0 /)
  real(8):: zbl_beta(1:4) = (/ &
       3.20d0, &
       0.94230d0, &
       0.40290d0, &
       0.20160d0 /)

contains
!=======================================================================
  subroutine read_params_ZBL(myid,mpi_world,iprint)
    use force, only: loverlay
    implicit none
    include "mpif.h"
    integer,intent(in):: myid,mpi_world,iprint

    integer:: isp,jsp,ierr
    real(8):: qnucli,ri,ro
    character(len=128):: cline,fname,cmode,ctmp
    real(8),parameter:: qtiny = 1d-10
    integer,external:: num_data

    if( myid.eq.0 ) then
      cmode = ''
      fname = trim(paramsdir)//'/'//trim(paramsfname)
      open(ioprms,file=trim(fname),status='old')
      interact(1:msp,1:msp) = .true.
      qnucl(1:msp) = 0d0
      zbl_rc = 0d0
      
      if( iprint.ne.0 ) write(6,'(/,a)') ' ZBL parameters:'
      do while(.true.)
        read(ioprms,*,end=10) cline
        if( num_data(cline,' ').eq.0 ) cycle
        if( cline(1:1).eq.'#' .or. cline(1:1).eq.'!' ) cycle
!.....Mode detection
        if( trim(cline).eq.'parameters' ) then
          cmode = trim(cline)
          cycle
        else if( trim(cline).eq.'interactions' ) then
          cmode = trim(cline)
          interact(1:msp,1:msp) = .false.
          cycle
        else if( trim(cline).eq.'overlay' ) then
          cmode = trim(cline)
          backspace(ioprms)
          read(ioprms,*) ctmp, loverlay
        endif
!.....Read parameters depending on the mode
        if( trim(cmode).eq.'parameters' ) then
          backspace(ioprms)
          read(ioprms,*) isp, qnucli, ri, ro
          if( isp.gt.msp ) then
            write(6,*) ' Warning @read_params: since isp is greater than msp,'&
                 //' skip reading the line.'
          endif
          qnucl(isp) = qnucli
          r_inner(isp) = ri
          r_outer(isp) = ro
          zbl_rc = max(zbl_rc,ro)
          if( iprint.ne.0 ) then
            write(6,'(a,i4,3(2x,f7.3))') &
                 '   isp,qnucl,ri,ro = ',isp,qnucli,ri,ro
          endif
        else if( trim(cmode).eq.'interactions' ) then
          backspace(ioprms)
          read(ioprms,*) isp, jsp
          interact(isp,jsp) = .true.
          interact(jsp,isp) = .true.
        endif
      enddo
10    close(ioprms)
!!$      do isp=1,msp
!!$        if( abs(qnucl(isp)).lt.qtiny ) cycle
!!$        do jsp=isp,msp
!!$          if( abs(qnucl(jsp)).lt.qtiny ) cycle
!!$          interact(isp,jsp) = .true.
!!$          interact(jsp,isp) = .true.
!!$        enddo
!!$      enddo
      if( iprint.gt.1 ) then
        do isp=1,msp
          do jsp=isp,msp
            if( interact(isp,jsp) ) then
              write(6,'(a,2i3,l2)') '   isp,jsp,interact = ',isp,jsp,interact(isp,jsp)
            endif
          enddo
        enddo
      endif
    endif

    call mpi_bcast(zbl_rc,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(qnucl,msp,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(r_inner,msp,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(r_outer,msp,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(zbl_rc,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(interact,msp*msp,mpi_logical,0,mpi_world,ierr)
    call mpi_bcast(loverlay,1,mpi_logical,0,mpi_world,ierr)
    return
  end subroutine read_params_ZBL
!=======================================================================
  subroutine force_ZBL(namax,natm,tag,ra,nnmax,aa,strs,h,hi,tcom &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,lstrs,iprint,l1st)
    implicit none
    include "mpif.h"
    include "./params_unit.h"
    integer,intent(in):: namax,natm,nnmax,nismax,lspr(0:nnmax,namax)&
         ,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_md_world,myid_md,nex(3)
    real(8),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),sv(3,6) &
         ,rc,tag(namax)
    real(8),intent(inout):: tcom
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: l1st
    logical:: lstrs

    integer:: i,j,k,l,m,n,ierr,is,js,ixyz,jxyz
    real(8):: xij(3),rij,rcij,dfi,dfj,drdxi(3),drdxj(3),r,at(3)
    real(8):: x,y,z,xi(3),epotl,epott,tmp,tmp2,dtmp
    real(8),allocatable,save:: strsl(:,:,:)
    real(8),external:: fcut1,dfcut1

    real(8),save:: zbl_rc2

    if( l1st ) then
      if( rc.lt.zbl_rc ) then
        if( myid_md.eq.0 .and. iprint.gt.0 ) then
          print '(/,a)',' Input cutoff radius is smaller than rc of ZBL potential.'
          print '(a,f0.3)', '   Input rc     = ',rc
          print '(a,f0.3)', '   Potential rc = ',zbl_rc
        endif
        call mpi_finalize(ierr)
        stop
      endif
      if( allocated(strsl) ) deallocate(strsl)
      allocate(strsl(3,3,namax))
      zbl_rc2 = zbl_rc*zbl_rc
    endif

    if( size(strsl).lt.3*3*namax ) then
      deallocate(strsl)
      allocate(strsl(3,3,namax))
    endif
    
    epotl= 0d0
    strsl(1:3,1:3,1:namax) = 0d0

!-----dE/dr_i
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is = int(tag(i))
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
        if(j.le.i) cycle
        js = int(tag(j))
        if( .not. interact(is,js) ) cycle
        x= ra(1,j) -xi(1)
        y= ra(2,j) -xi(2)
        z= ra(3,j) -xi(3)
        xij(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
        rij= xij(1)**2+ xij(2)**2 +xij(3)**2
        if( rij.gt.zbl_rc2 ) cycle
        rij = sqrt(rij)
        drdxi(1:3)= -xij(1:3)/rij
!.....2-body term
        tmp = vij(is,js,rij)
        tmp2 = 0.5d0 *tmp
        if(j.le.natm) then
          epi(i)= epi(i) +tmp2
          epi(j)= epi(j) +tmp2
          epotl=epotl +tmp2 +tmp2
        else
          epi(i)= epi(i) +tmp2
          epotl=epotl +tmp2
        endif
!!$        dtmp = dvij(is,js,rij)*fcut1(rij,zbl_rc) &
!!$             +tmp *dfcut1(rij,zbl_rc)
        dtmp = dvij(is,js,rij)
        aa(1:3,i)=aa(1:3,i) -dtmp*drdxi(1:3)
        aa(1:3,j)=aa(1:3,j) +dtmp*drdxi(1:3)
!.....Atomic stress for 2-body terms
        if( lstrs ) then
          do ixyz=1,3
            do jxyz=1,3
              strsl(jxyz,ixyz,i)=strsl(jxyz,ixyz,i) &
                   -0.5d0*dtmp*xij(ixyz)*(-drdxi(jxyz))
              strsl(jxyz,ixyz,j)=strsl(jxyz,ixyz,j) &
                   -0.5d0*dtmp*xij(ixyz)*(-drdxi(jxyz))
            enddo
          enddo
        endif
      enddo
    enddo

    if( lstrs ) then
      strs(1:3,1:3,1:natm)= strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)
    endif

!-----gather epot
    call mpi_allreduce(epotl,epott,1,mpi_real8 &
         ,mpi_sum,mpi_md_world,ierr)
    if( iprint.gt.2 ) print *,'ZBL epot = ',epott
    epot= epot +epott

  end subroutine force_ZBL
!=======================================================================
  function vij(is,js,rij)
!
! Main two-body function.
!
    implicit none
    integer,intent(in):: is,js
    real(8),intent(in):: rij
    real(8):: vij
    real(8):: ri,ro,veqt

    ri = (r_inner(is)+r_inner(js))/2
    ro = (r_outer(is)+r_outer(js))/2
    vij = 0d0
    if( rij.lt.ri ) then
      vij = vnucl(is,js,rij)
    else if( rij.ge.ri .and. rij.lt.ro ) then
      vij = zeta((ro+ri-2d0*rij)/(ro-ri)) &
           *vnucl(is,js,rij)
    endif
    return
  end function vij
!=======================================================================
  function dvij(is,js,rij)
!
!  Derivative of the main two-body function.
!
    implicit none
    integer,intent(in):: is,js
    real(8),intent(in):: rij
    real(8):: dvij
    real(8):: ri,ro,x

    ri = (r_inner(is)+r_inner(js))/2
    ro = (r_outer(is)+r_outer(js))/2
    dvij = 0d0
    if( rij.lt.ri ) then
      dvij = dvnucl(is,js,rij)
    else if( rij.ge.ri .and. rij.lt.ro ) then
      x = (ro+ri-2d0*rij)/(ro-ri)
      dvij = dzeta(x)*(-2d0/(ro-ri)) &
           *vnucl(is,js,rij) &
           +zeta(x)*dvnucl(is,js,rij)
    endif
    return
  end function dvij
!=======================================================================
  function vnucl(is,js,rij)
!
!  Repulsive potential between nuclei
!
    implicit none
    integer,intent(in):: is,js
    real(8),intent(in):: rij
    real(8):: vnucl
    real(8):: rs,qi,qj

    qi = qnucl(is)
    qj = qnucl(js)
    rs = zbl_aa  /(qi**zbl_gamma +qj**zbl_gamma)
    vnucl = qi*qj/rij *xi(rij/rs)
    return
  end function vnucl
!=======================================================================
  function dvnucl(is,js,rij)
!
!  Derivative of the repulsive potential between nuclei
!
    implicit none
    integer,intent(in):: is,js
    real(8),intent(in):: rij
    real(8):: dvnucl
    real(8):: rs,qi,qj

    qi = qnucl(is)
    qj = qnucl(js)
    rs = zbl_aa  /(qi**zbl_gamma +qj**zbl_gamma)
    dvnucl = -qi*qj/rij* ( 1d0/rij*xi(rij/rs) -dxi(rij/rs)/rs )
    return
  end function dvnucl
!=======================================================================
  function xi(x)
    implicit none
    real(8),intent(in):: x
    real(8):: xi
    integer:: i

    xi = 0d0
    do i=1,4
      xi = xi +zbl_alpha(i)*exp(-zbl_beta(i)*x)
    enddo
!!$    xi= 0.1818d0*exp(-3.2d0*x) &
!!$         +0.5099d0*exp(-0.9423d0*x) &
!!$         +0.2802d0*exp(-0.4029d0*x) &
!!$         +0.02817d0*exp(-0.2016d0*x)
    return
  end function xi
!=======================================================================
  function dxi(x)
    implicit none
    real(8),intent(in):: x
    real(8):: dxi
    integer:: i

    dxi = 0d0
    do i=1,4
      dxi = dxi -zbl_beta(i)*zbl_alpha(i)*exp(-zbl_beta(i)*x)
    enddo
!!$    dxi= -0.58176d0*exp(-3.2d0*x) &
!!$         -0.48047877d0*exp(-0.9423d0*x) &
!!$         -0.11289258d0*exp(-0.4029d0*x) &
!!$         -0.005679072d0*exp(-0.2016d0*x)
    return
  end function dxi
!=======================================================================
  function zeta(x)
    implicit none
    real(8),intent(in):: x
    real(8):: zeta

    zeta = (3d0*x**5 -10d0*x**3 +15d0*x +8d0)/16d0
    return
  end function zeta
!=======================================================================
  function dzeta(x)
    implicit none
    real(8),intent(in):: x
    real(8):: dzeta

    dzeta = (15d0*x**4 -30d0*x**2 +15d0)/16d0
    return
  end function dzeta
!=======================================================================
  subroutine set_paramsdir_ZBL(dname)
!
!  Accessor routine to set paramsdir.
!
    implicit none
    character(len=*),intent(in):: dname

    paramsdir = trim(dname)
    return
  end subroutine set_paramsdir_ZBL
!=======================================================================
  
end module ZBL
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
