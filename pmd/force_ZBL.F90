module ZBL
!-----------------------------------------------------------------------
!                     Last modified: <2024-11-23 22:06:16 KOBAYASHI Ryo>
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
  use pmdvars,only: nspmax
  use util,only: csp2isp
  use force,only: loverlay
  implicit none
  include "./const.h"
  save

  character(len=128):: paramsdir = '.'
  character(len=128),parameter:: paramsfname = 'in.params.ZBL'
  integer,parameter:: ioprms = 20
  logical:: lprmset_zbl = .false.
  
!.....Max num of species
  integer,parameter:: msp = nspmax

!.....Coulomb's constant, acc = 1.0/(4*pi*epsilon0)
  real(8),parameter:: acc  = 14.3998554737d0
!.....permittivity of vacuum
  real(8),parameter:: eps0 = 0.00552634939836d0  ! e^2 /Ang /eV

  logical:: interact(1:msp,1:msp) = .true.
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
  subroutine init_ZBL(iprint)
    use pmdvars,only: specorder
    use element,only: nelem,elmts
    integer,intent(in):: iprint

    integer:: isp,iz
    character(len=3):: csp

!.....Set qnucl from specorder
    do isp=1,nspmax
      csp = trim(specorder(isp))
      if( trim(csp).eq.'x' ) cycle
      do iz=1,nelem
        if( trim(csp).eq.trim(elmts(iz)%symbol) ) then
          qnucl(isp) = dble(iz)
          exit
        endif
      enddo
    enddo
    return
  end subroutine init_ZBL
!=======================================================================
  subroutine read_params_ZBL(myid,mpi_world,iprint,specorder)
    use force, only: loverlay
    use util, only: num_data
    implicit none
    include "mpif.h"
    integer,intent(in):: myid,mpi_world,iprint
    character(len=3),intent(in):: specorder(nspmax)

    integer:: isp,jsp,ierr
    real(8):: qnucli,ri,ro
    character(len=128):: cline,fname,cmode,ctmp
    character(len=5):: cspi,cspj
    real(8),parameter:: qtiny = 1d-10
!!$    integer,external:: num_data

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
        endif
!.....Read parameters depending on the mode
        if( trim(cmode).eq.'parameters' ) then
          backspace(ioprms)
          read(ioprms,*) cspi, qnucli, ri, ro
          isp = csp2isp(trim(cspi))
          if( isp.le.0 ) cycle
          qnucl(isp) = qnucli
          r_inner(isp) = ri
          r_outer(isp) = ro
          zbl_rc = max(zbl_rc,ro)
          if( iprint.ne.0 ) then
            write(6,'(a,a3,3(2x,f7.3))') &
                 '   csp,qnucl,ri,ro = ',trim(cspi),qnucli,ri,ro
          endif
        else if( trim(cmode).eq.'interactions' ) then
          backspace(ioprms)
          read(ioprms,*) cspi, cspj
          isp = csp2isp(trim(cspi))
          jsp = csp2isp(trim(cspj))
          if( isp.gt.0 .and. jsp.gt.0 ) then
            interact(isp,jsp) = .true.
            interact(jsp,isp) = .true.
          else
            print *,'  interacion read but not used: ',isp,jsp
          endif
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
      if( iprint.ge.ipl_info ) then
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
    call mpi_bcast(interact,msp*msp,mpi_logical,0,mpi_world,ierr)
    lprmset_zbl = .true.
    return
  end subroutine read_params_ZBL
!=======================================================================
  subroutine force_ZBL(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
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
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: l1st
    logical:: lstrs

    integer:: i,j,k,l,m,n,ierr,is,js,ixyz,jxyz
    real(8):: xij(3),rij,rij2,rcij,dfi,dfj,drdxi(3),drdxj(3),r,at(3)
    real(8):: x,y,z,xi(3),epotl,epott,tmp,tmp2,dtmp
    real(8),allocatable,save:: strsl(:,:,:)
    real(8),external:: fcut1,dfcut1

    real(8),save:: zbl_rc2

    if( l1st ) then
      if( rc.lt.zbl_rc ) then
        if( myid_md.eq.0 .and. iprint.ge.ipl_basic ) then
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
        rij2= xij(1)**2+ xij(2)**2 +xij(3)**2
        if( rij2.gt.zbl_rc2 ) cycle
        rij = sqrt(rij2)
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
        dtmp = dvij(is,js,rij)
        aa(1:3,i)=aa(1:3,i) -dtmp*drdxi(1:3)
        aa(1:3,j)=aa(1:3,j) +dtmp*drdxi(1:3)
!.....Atomic stress
        if( .not.lstrs ) cycle
        do ixyz=1,3
          do jxyz=1,3
            strsl(jxyz,ixyz,i)=strsl(jxyz,ixyz,i) &
                 -0.5d0*dtmp*xij(ixyz)*(-drdxi(jxyz))
            strsl(jxyz,ixyz,j)=strsl(jxyz,ixyz,j) &
                 -0.5d0*dtmp*xij(ixyz)*(-drdxi(jxyz))
          enddo
        enddo
      enddo
    enddo

    if( lstrs ) then
      strs(1:3,1:3,1:natm)= strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)
    endif

!-----gather epot
    call mpi_allreduce(epotl,epott,1,mpi_real8 &
         ,mpi_sum,mpi_md_world,ierr)
    if( iprint.ge.ipl_info ) print *,'epot ZBL = ',epott
    epot= epot +epott

  end subroutine force_ZBL
!=======================================================================
  subroutine force_ZBL_overlay(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,lstrs,iprint,l1st)
!
!  ZBL potential used as overlaying potential for close ion-ion distance.
!  
!  Since the energy of atom-i is defined using alpha_i, and the force-ji
!  and force-ij is not symmetry (I think), the loop over i and j has to be
!  fully taken.
!
    use force,only: ol_ranges,ol_alphas,ol_dalphas,ol_pair,ol_type
    implicit none
    include "mpif.h"
    include "./params_unit.h"
    integer,intent(in):: namax,natm,nnmax,nismax,lspr(0:nnmax,namax)&
         ,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_md_world,myid_md,nex(3)
    real(8),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),sv(3,6) &
         ,rc,tag(namax)
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: l1st
    logical:: lstrs

    integer:: i,j,k,l,m,n,ierr,is,js,ixyz,jxyz
    real(8):: xij(3),rij,rij2,rcij,dfi,dfj,drdxi(3),drdxj(3),r,at(3)
    real(8):: x,y,z,xi(3),epotl,epott,tmp,tmp2,dtmp,tmp2i,tmp2j &
         ,dtmpi,dtmpj,alpi,ri,ro,dij,dij2
    real(8),allocatable,save:: strsl(:,:,:),epit(:),aal(:,:)
    real(8),external:: fcut1,dfcut1
    real(8),save:: zbl_rc2

    if( l1st ) then
      if( allocated(strsl) ) deallocate(strsl,epit,aal)
      allocate(strsl(3,3,namax),epit(namax),aal(3,namax))
      zbl_rc2 = rc*rc
    endif

    if( size(strsl).lt.3*3*namax ) then
      deallocate(strsl,epit,aal)
      allocate(strsl(3,3,namax),epit(namax),aal(3,namax))
    endif

    epotl= 0d0
    epit(1:natm+nb) = 0d0
    strsl(1:3,1:3,1:natm+nb) = 0d0
    aal(:,1:natm+nb)= 0d0

    do i=1,natm
      xi(1:3) = ra(1:3,i)
      is = int(tag(i))
      alpi = ol_alphas(0,i)
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
        js = int(tag(j))
!!$        if( .not. interact(is,js) ) cycle
        x= ra(1,j) -xi(1)
        y= ra(2,j) -xi(2)
        z= ra(3,j) -xi(3)
        xij(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
        dij2= xij(1)*xij(1)+ xij(2)*xij(2) +xij(3)*xij(3)
        if( dij2.gt.zbl_rc2 ) cycle
        dij = sqrt(dij2)
        drdxi(1:3)= -xij(1:3)/dij
!.....2-body term, taking overlay into account
        tmp = 0.5d0 *vnucl(is,js,dij) 
        epit(i)= epit(i) +tmp
        epotl=epotl +tmp *(1d0-alpi)
!.....Force
        dtmp = 0.5d0 *dvnucl(is,js,dij) *(1d0 -alpi)
        aal(1:3,i)=aal(1:3,i) -dtmp*drdxi(1:3)
        aal(1:3,j)=aal(1:3,j) +dtmp*drdxi(1:3)
!.....Stress
        if( .not.lstrs ) cycle
        do ixyz=1,3
          do jxyz=1,3
            strsl(jxyz,ixyz,i)=strsl(jxyz,ixyz,i) &
                 -0.5d0*dtmp*xij(ixyz)*(-drdxi(jxyz))
            strsl(jxyz,ixyz,j)=strsl(jxyz,ixyz,j) &
                 -0.5d0*dtmp*xij(ixyz)*(-drdxi(jxyz))
          enddo
        enddo
      enddo  ! k=
      epi(i) = epi(i) +epit(i) *(1d0-alpi)
    enddo

!.....Derivatives wrt alphas, which requires energy per atom, epit
    do i=1,natm
      xi(1:3) = ra(1:3,i)
      is = int(tag(i))
      alpi = ol_alphas(0,i)
      do k=1,lspr(0,i)
        j= lspr(k,i)
        js = int(tag(j))
!!$        if( .not. interact(is,js) ) cycle
        ri = ol_pair(1,is,js)
        ro = ol_pair(2,is,js)
        x= ra(1,j) -xi(1)
        y= ra(2,j) -xi(2)
        z= ra(3,j) -xi(3)
        xij(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
        dij2= xij(1)*xij(1)+ xij(2)*xij(2) +xij(3)*xij(3)
        if( dij2.gt.ro*ro .or. dij2.le.ri*ri ) cycle
        dij = dsqrt(dij2)
        drdxi(1:3)= -xij(1:3)/dij
!.....Force
        dtmp = -epit(i) *alpi/ol_alphas(k,i) *ol_dalphas(k,i)
        aal(1:3,i)=aal(1:3,i) -dtmp*drdxi(1:3)
        aal(1:3,j)=aal(1:3,j) +dtmp*drdxi(1:3)
!.....Stress
        if( .not.lstrs ) cycle
        do ixyz=1,3
          do jxyz=1,3
            strsl(jxyz,ixyz,i)=strsl(jxyz,ixyz,i) &
                 -0.5d0*dtmp*xij(ixyz)*(-drdxi(jxyz))
            strsl(jxyz,ixyz,j)=strsl(jxyz,ixyz,j) &
                 -0.5d0*dtmp*xij(ixyz)*(-drdxi(jxyz))
          enddo
        enddo
      enddo
    enddo  ! i=

    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex &
         ,lsrc,myparity,nn,mpi_md_world,aal,3)
    aa(1:3,1:natm) = aa(1:3,1:natm) +aal(1:3,1:natm)

    if( lstrs ) then
      strs(1:3,1:3,1:natm)= strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)
    endif

!-----gather epot
    call mpi_allreduce(epotl,epott,1,mpi_real8,mpi_sum,mpi_md_world,ierr)
    if( iprint.ge.ipl_info ) print *,'epot ZBL_overlay = ',epott
    epot= epot +epott

  end subroutine force_ZBL_overlay
!=======================================================================
  function vij(is,js,rij)
!
!   Main two-body function.
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
    vnucl = acc *qi*qj/rij *xi(rij/rs)
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
    dvnucl = -acc *qi*qj/rij* ( 1d0/rij*xi(rij/rs) -dxi(rij/rs)/rs )
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
  subroutine set_params_ZBL(rcin,qnuclin,riin,roin,interactin)
!
!  Accessor routine to set ZBL parameters from outside (fitpot).
!  Curretnly this routine is supposed to be called only on serial run.
!
    real(8),intent(in):: rcin,qnuclin(nspmax),riin(nspmax),roin(nspmax)
    logical,intent(in):: interactin(nspmax,nspmax)

    if( lprmset_zbl ) return

    zbl_rc = rcin
    qnucl(:) = qnuclin(:)
    r_inner(:) = riin(:)
    r_outer(:) = roin(:)
    interact(:,:) = interactin(:,:)

    lprmset_zbl = .true.
    return
  end subroutine set_params_ZBL
!=======================================================================
  subroutine set_qnucl_ZBL(qnuclin,interactin)
!
!  Accessor routine to set ZBL parameters from outside (fitpot).
!  Curretnly this routine is supposed to be called only on serial run.
!
    real(8),intent(in):: qnuclin(nspmax)
    logical,intent(in):: interactin(nspmax,nspmax)

    qnucl(:) = qnuclin(:)
    interact(:,:) = interactin(:,:)
  end subroutine set_qnucl_ZBL
!=======================================================================
  subroutine get_dvnucl(is,js,npnts,rijs,dvnucls)
!
!  Accessor routine to set ZBL parameters from outside (fitpot).
!  Curretnly this routine is supposed to be called only on serial run.
!
    integer,intent(in):: is,js,npnts
    real(8),intent(in):: rijs(npnts)
    real(8),intent(out):: dvnucls(npnts)
    integer:: i
    real(8):: qi,qj,rs,rij

    qi = qnucl(is)
    qj = qnucl(js)

    do i=1,npnts
      rs = zbl_aa  /(qi**zbl_gamma +qj**zbl_gamma)
      rij = rijs(i)
      dvnucls(i) = -acc *qi*qj/rij &
           *( 1d0/rij*xi(rij/rs) -dxi(rij/rs)/rs )
    enddo
    return
  end subroutine get_dvnucl
end module ZBL
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
