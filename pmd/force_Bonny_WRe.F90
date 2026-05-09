module Bonny_WRe
!-----------------------------------------------------------------------
!                     Last modified: <2023-01-23 17:09:20 KOBAYASHI Ryo>
!-----------------------------------------------------------------------
!  Parallel implementation of EAM poetntial of Bonney et al.
!  See G. Bonny et al., J. Appl. Phys. 121, 165107 (2017).
!-----------------------------------------------------------------------
!  Species 1 should be W and 2 should be Re.
!-----------------------------------------------------------------------
  use pmdmpi
  use mod_precision
  implicit none
  include "./const.h"
  character(len=128):: paramsdir = ''

  real(rp),external:: hvsd

!.....Max num of species
  integer,parameter:: msp = 9

  logical:: interact(msp,msp)

!.....Coulomb's constant, acc = 1.0/(4*pi*epsilon0)
  real(rp),parameter:: acc  = 14.3998554737_rp
!.....permittivity of vacuum
  real(rp),parameter:: eps0 = 0.00552634939836_rp  ! e^2 /Ang /eV

  real(rp):: qnucl(1:2) = (/ &
       74.0_rp, & ! 1, W
       75.0_rp  & ! 2, Re
       /)

  real(rp):: r_inner(1:2,1:2) = reshape( (/ 1.0_rp, 1.0_rp, 1.0_rp, 1.0_rp /), &
       shape(r_inner) )
  real(rp):: r_outer(1:2,1:2) = reshape( (/ 2.0_rp, 2.0_rp, 2.0_rp, 2.0_rp /), &
       shape(r_outer) )

  real(rp),parameter:: bonny_rc(1:2,1:2) = reshape( &
       (/ 5.460_rp, 3.825_rp , &  ! 1-1 W-W, 1-2 W-Re
          3.825_rp, 5.460_rp /) &  ! 2-1 Re-W, 2-2 Re-Re
          , shape(bonny_rc) )

!.....Pure W parameters
  real(rp),parameter:: gauge_C = 1.848055990_rp
  real(rp),parameter:: gauge_S = 0.2232322602_rp

  real(rp),parameter:: A0_W  = -5.524855802_rp
  real(rp),parameter:: A1_W  =  2.317313103e-1_rp
  real(rp),parameter:: A2_W  = -3.665345949e-2_rp
  real(rp),parameter:: A3_W  =  8.989367404e-3_rp
  real(rp),parameter:: rhoi_W = 1.359141225_rp
  
  real(rp),parameter:: rhospln_rc = 2.002970124727_rp
  integer,parameter:: n_rhospln = 4
  real(rp),parameter:: rhospln_a(1:4) = (/&
       -0.420429107805055e+1_rp, &
        0.518217702261442_rp,  &
        0.562720834534370e-1_rp, &
        0.344164178842340e-1_rp &
        /)
  real(rp),parameter:: rhospln_r(1:4) = (/&
       2.5_rp, &
       3.1_rp, &
       3.5_rp, &
       4.9_rp &
       /)
  real(rp),parameter:: fspln_a1 = -5.946454472402710_rp
  real(rp),parameter:: fspln_a2 = -0.049477376935239_rp
  integer,parameter:: n_vspln = 15
  real(rp),parameter:: vspln_a(1:15) = (/ &
        0.960851701343041e+2_rp, &
       -0.184410923895214e+3_rp, &
        0.935784079613550e+2_rp, &
       -0.798358265041677e+1_rp, &
        0.747034092936229e+1_rp, &
       -0.152756043708453e+1_rp, &
        0.125205932634393e+1_rp, &
        0.163082162159425e+1_rp, &
       -0.141854775352260e+1_rp, &
       -0.819936046256149_rp,  &
        0.198013514305908e1_rp,  &
       -0.696430179520267_rp,  &
        0.304546909722160e-1_rp, &
       -0.163131143161660e+1_rp, &
        0.138409896486177e+1_rp  &
        /)
  real(rp),parameter:: vspln_r(1:15) = (/ &
       2.5648975_rp, &
       2.6297950_rp, &
       2.6946925_rp, &
       2.8663175_rp, &
       2.9730450_rp, &
       3.0797725_rp, &
       3.5164725_rp, &
       3.8464450_rp, &
       4.1764175_rp, &
       4.7008450_rp, &
       4.8953000_rp, &
       5.0897550_rp, &
       5.3429525_rp, &
       5.4016950_rp, &
       5.4600000_rp  &  ! shortened from original value 5.4604375
       /)
  
!.....Pure Re parameters
  real(rp),parameter:: A_Re = -7.046791948_rp
  real(rp),parameter:: B_Re =  1.236584720_rp
  real(rp),parameter:: C_Re =  1.143405627_rp
  real(rp),parameter:: C0_Re=  3.704045964e-3_rp

  integer,parameter:: n_veq_ReRe = 6
  real(rp),parameter:: veq_ReRe_a(1:6) = (/ &
        6.726805309_rp, &
        3.217593889_rp, &
       -6.545857587e-1_rp, &
        1.453229484e-1_rp, &
       -2.063629464e-1_rp, &
        6.114909116e-2_rp  &
        /)
  real(rp),parameter:: veq_ReRe_r(1:6) = (/ &
        2.700_rp, &
        3.252_rp, &
        3.804_rp, &
        4.356_rp, &
        4.908_rp, &
        5.460_rp  &
        /)

!.....W-Re parameters
  integer,parameter:: n_veq_WRe = 6
  real(rp),parameter:: veq_WRe_a(1:6) = (/ &
       -2.335000000e+1_rp,  &
        2.456959229e+1_rp,  &
       -2.585878138_rp,   &
        3.457586051_rp,   &
       -7.013105493e-1_rp,  &
       -0.25133241_rp     &
       /)
  real(rp),parameter:: veq_WRe_r(1:6) = (/ &
       2.650_rp, &
       2.700_rp, &
       3.075_rp, &
       3.450_rp, &
       3.825_rp, &
       4.200_rp  &
       /)
  
contains
  subroutine force_Bonny_WRe(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,lstrs,iprint,l1st)
    implicit none
    include "./params_unit.h"
    integer,intent(in):: namax,natm,nnmax,nismax,lspr(0:nnmax,namax)&
         ,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_md_world,myid_md,nex(3)
    real(rp),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),sv(3,6) &
         ,rc,tag(namax)
    real(rp),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: l1st
    logical:: lstrs

    integer:: i,j,k,l,m,n,ierr,is,js,ixyz,jxyz
    real(rp):: xij(3),rij(3),dij,dij2,rcij,dfi,dfj,drdxi(3),drdxj(3),r,at(3) &
         ,xi(3),xj(3),epotl,epott,tmp,dtmp,drhoi,drhoj,d
    real(rp),allocatable,save:: rho(:)
    real(rp),allocatable,save:: strsl(:,:,:)

    real(rp),save:: rcmax,rcmax2

    if( l1st ) then
      if( rc.lt.bonny_rc(1,1) ) then
        if( myid_md.eq.0 .and. iprint.ge.ipl_basic ) then
          print '(/,a)',' Input cutoff radius is smaller than rc of Bonny potential.'
          print '(a,f0.3)', '   Input rc     = ',rc
          print '(a,f0.3)', '   Potential rc = ',bonny_rc(1,1)
        endif
        call mpi_finalize(ierr)
        stop
      endif
      interact(1:msp,1:msp) = .false.
      interact(1:2,1:2) = .true.
      rcmax = 0.0_rp
      do is=1,2
        do js=1,2
          rcmax = max(rcmax, bonny_rc(is,js))
        enddo
      enddo
      rcmax2 = rcmax*rcmax
      if( myid_md.eq.0 .and. iprint.ge.ipl_basic )  then
        print '(/,a,f0.3)', ' Max cutoff in Bonny potential = ',rcmax
      endif

      if( allocated(rho) ) deallocate(rho)
      allocate(rho(namax))
      if( allocated(strsl) ) deallocate(strsl)
      allocate(strsl(3,3,namax))
    endif

    if( size(strsl).lt.3*3*namax ) then
      deallocate(strsl)
      allocate(strsl(3,3,namax))
    endif
    if( size(rho).lt.namax ) then
      deallocate(rho)
      allocate(rho(namax))
    endif
    
    epotl= 0.0_rp
    rho(1:namax)= 0.0_rp
    strsl(1:3,1:3,1:namax) = 0.0_rp

!-----rho(i)
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is = int(tag(i))
      if( is.gt.2 ) cycle
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
        js = int(tag(j))
        if( .not. interact(is,js) ) cycle
        xj(1:3) = ra(1:3,j)
        xij(1:3)= xj(1:3) -xi(1:3)
        rij(1:3)= h(1:3,1,0)*xij(1) +h(1:3,2,0)*xij(2) +h(1:3,3,0)*xij(3)
        dij2=rij(1)*rij(1)+ rij(2)*rij(2) +rij(3)*rij(3)
!!$        rij = dlspr(0,k,i)
        if( dij2.gt.rcmax2 ) cycle
        dij= sqrt(dij2)
        rho(i) = rho(i) +rhoij(js,dij)
      enddo
    enddo

!-----copy rho of boundary atoms
    call copy_dba_fwd(namax,natm,nb,nbmax,lsb,nex,&
         lsrc,myparity,nn,sv,mpi_md_world,rho,1)

!-----dE/dr_i
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is = int(tag(i))
      if( is.gt.2 ) cycle
      dfi = dfrho(is,rho(i))
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
        if(j.le.i) cycle
        js = int(tag(j))
        if( .not. interact(is,js) ) cycle
        xj(1:3)= ra(1:3,j)
        xij(1:3)= xj(1:3) -xi(1:3)
        rij(1:3)= h(1:3,1,0)*xij(1) +h(1:3,2,0)*xij(2) +h(1:3,3,0)*xij(3)
        dij2=rij(1)*rij(1)+ rij(2)*rij(2) +rij(3)*rij(3)
        if( dij2.gt.rcmax2 ) cycle
        dij= sqrt(dij2)
!!$        rij=dsqrt(xij(1)**2+ xij(2)**2 +xij(3)**2)
!!$        rij = dlspr(0,k,i)
!!$        if( rij.gt.rcmax ) cycle
!!$        xij(1:3) = dlspr(1:3,k,i)
        drdxi(1:3)= -rij(1:3)/dij
!.....2-body term
        tmp = 0.5_rp *vij(is,js,dij)
        epi(i)= epi(i) +tmp
        epi(j)= epi(j) +tmp
        if(j.le.natm) then
          epotl=epotl +tmp +tmp
        else
          epotl=epotl +tmp
        endif
        dtmp = dvij(is,js,dij)
        aa(1:3,i)=aa(1:3,i) -dtmp*drdxi(1:3)
        aa(1:3,j)=aa(1:3,j) +dtmp*drdxi(1:3)
!.....Atomic stress for 2-body terms
        if( lstrs ) then
          do ixyz=1,3
            do jxyz=1,3
              strsl(jxyz,ixyz,i)=strsl(jxyz,ixyz,i) &
                   -0.5_rp*dtmp*rij(ixyz)*(-drdxi(jxyz))
              strsl(jxyz,ixyz,j)=strsl(jxyz,ixyz,j) &
                   -0.5_rp*dtmp*rij(ixyz)*(-drdxi(jxyz))
            enddo
          enddo
        endif
!.....Embedded term
        if( rho(i).lt.1e-10_rp .or. rho(j).lt.1e-10_rp ) cycle
        drhoi = drhoij(is,dij)
        drhoj = drhoij(js,dij)
        dfj = dfrho(js,rho(j))
        tmp = dfi*drhoj + dfj*drhoi
        aa(1:3,i)=aa(1:3,i) -tmp*drdxi(1:3)
        aa(1:3,j)=aa(1:3,j) +tmp*drdxi(1:3)
!.....Atomic stress of many-body contributions
        if( lstrs ) then
          do ixyz=1,3
            do jxyz=1,3
              strsl(jxyz,ixyz,i)=strsl(jxyz,ixyz,i) &
                   -0.5_rp*tmp*rij(ixyz)*(-drdxi(jxyz))
              strsl(jxyz,ixyz,j)=strsl(jxyz,ixyz,j) &
                   -0.5_rp*tmp*rij(ixyz)*(-drdxi(jxyz))
            enddo
          enddo
        endif
      enddo
      tmp = frho(is,rho(i))
      epi(i)=epi(i) +tmp
      epotl=epotl +tmp
    enddo

    if( lstrs ) then
!!$      call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
!!$           ,nn,mpi_md_world,strsl,9)
      strs(1:3,1:3,1:natm)= strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)
    endif

!-----gather epot
    call mpi_allreduce(epotl,epott,1,mpi_real_rp &
         ,mpi_sum,mpi_md_world,ierr)
    epot= epot +epott

  end subroutine force_Bonny_WRe
!=======================================================================
  function rhoij(js,rij)
!
! Calculate rho of atom j at distance rij.
!
    real(rp),intent(in):: rij
    integer,intent(in):: js
    real(rp):: rhoij

    rhoij = 0.0_rp
    if( js.eq.1 ) then  ! Only in case of W, scaling with S
      rhoij = gauge_S *rhospln(rij)
    else if( js.eq.2 ) then
      rhoij = C0_Re *(bonny_rc(2,2) -rij)**3 *hvsd(bonny_rc(2,2) -rij)
    endif
    return
  end function rhoij
!=======================================================================
  function drhoij(js,rij)
    implicit none
    real(rp),intent(in):: rij
    integer,intent(in):: js
    real(rp):: drhoij

    drhoij = 0.0_rp
    if( js.eq.1 ) then
      drhoij = gauge_S *drhospln(rij)
    else if( js.eq.2 ) then
      drhoij = -3.0_rp *C0_Re *(bonny_rc(2,2) -rij)**2 *hvsd(bonny_rc(2,2) -rij)
    endif
    return
  end function drhoij
!=======================================================================
  function rhospln(rij)
!
!  rho_j(rij) from Marinica et al., JAP 121 (2017)
!
    implicit none
    real(rp),intent(in):: rij
    real(rp):: rhospln,ri
    integer:: i

    rhospln = 0.0_rp
    if( rij.le.rhospln_rc ) then
      do i=1,n_rhospln
        ri = rhospln_r(i)
        rhospln = rhospln +rhospln_a(i)*(ri -rhospln_rc)**3 &
             *hvsd(ri -rhospln_rc)
      enddo
    else
      do i=1,n_rhospln
        ri = rhospln_r(i)
        rhospln = rhospln +rhospln_a(i)*(ri -rij)**3 &
             *hvsd(ri -rij)
      enddo
    endif
    return
  end function rhospln
!=======================================================================
  function drhospln(rij)
!
!  rho_j(rij) from Marinica et al., JAP 121 (2017)
!
    implicit none
    real(rp),intent(in):: rij
    real(rp):: drhospln,ri
    integer:: i

    drhospln = 0.0_rp
    if( rij.le.rhospln_rc ) then
      drhospln = 0.0_rp
    else
      do i=1,n_rhospln
        ri = rhospln_r(i)
        drhospln = drhospln -rhospln_a(i)*(ri -rij)**2 &
             *hvsd(ri -rij)
      enddo
      drhospln = drhospln*3.0_rp
    endif
    return
  end function drhospln
!=======================================================================
  function frho(is,rho)
    implicit none
    integer,intent(in):: is
    real(rp),intent(in):: rho
    real(rp):: frho

    if( is.eq.1 ) then  ! W
      if( rho.le.rhoi_W ) then
        frho = feff(rho)
      else
        frho = A0_W +A1_W*rho &
             +A2_W*rho*rho +A3_W*rho*rho*rho
      endif
    else if( is.eq.2 ) then  ! Re
      frho = A_Re*sqrt(rho) +B_Re*rho +C_Re*rho*rho
    endif
    
  end function frho
!=======================================================================
  function dfrho(is,rho)
    implicit none
    integer,intent(in):: is
    real(rp),intent(in):: rho
    real(rp):: dfrho

    if( is.eq.1 ) then  ! W
      if( rho.le.rhoi_W ) then
        dfrho = dfeff(rho)
      else
        dfrho = A1_W +2.0_rp*A2_W*rho +3.0_rp*A3_W*rho*rho
      endif
    else if( is.eq.2 ) then  ! Re
      dfrho = 0.5_rp*A_Re/sqrt(rho) +B_Re +2.0_rp*C_Re*rho
    endif
    
  end function dfrho
!=======================================================================
  function feff(rho)
!
!  F^{eff} for pure W
!
    implicit none
    real(rp),intent(in):: rho
    real(rp):: feff

    feff = fspln(rho/gauge_S) +gauge_C/gauge_S*rho
    return
  end function feff
!=======================================================================
  function dfeff(rho)
!
!  Derivative of F^{eff} for pure W
!
    implicit none
    real(rp),intent(in):: rho
    real(rp):: dfeff

    dfeff = dfspln(rho/gauge_S)/gauge_S +gauge_C/gauge_S
    return
  end function dfeff
!=======================================================================
  function fspln(rho)
!
!  F[rho] function from Marinica, JAP 121, 165107 (2017)
!
    implicit none
    real(rp),intent(in):: rho
    real(rp):: fspln

    fspln = fspln_a1*sqrt(rho) +fspln_a2*rho*rho
    return
  end function fspln
!=======================================================================
  function dfspln(rho)
!
!  Derivative of F[rho] function from Marinica et al.
!
    implicit none
    real(rp),intent(in):: rho
    real(rp):: dfspln

    dfspln = 0.5_rp*fspln_a1/sqrt(rho) +2.0_rp*fspln_a2*rho
    return
  end function dfspln
!=======================================================================
  function vij(is,js,rij)
!
! Main two-body function.
!
    implicit none
    integer,intent(in):: is,js
    real(rp),intent(in):: rij
    real(rp):: vij
    real(rp):: ri,ro,veqt

    ri = r_inner(is,js)
    ro = r_outer(is,js)
    if( rij.lt.ri ) then
      vij = vnucl(is,js,rij)
    else if( rij.ge.ri .and. rij.lt.ro ) then
      veqt = veq(is,js,rij)
      vij = veqt +zeta((ro+ri-2.0_rp*rij)/(ro-ri)) &
           *(vnucl(is,js,rij) -veqt)
    else if( rij.ge.ro ) then
      vij = veq(is,js,rij)
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
    real(rp),intent(in):: rij
    real(rp):: dvij
    real(rp):: ri,ro

    ri = r_inner(is,js)
    ro = r_outer(is,js)
    if( rij.lt.ri ) then
      dvij = dvnucl(is,js,rij)
    else if( rij.ge.ri .and. rij.lt.ro ) then
      dvij = dveq(is,js,rij) +dzeta((ro+ri-2.0_rp*rij)/(ro-ri))*(-2.0_rp/(ro-ri)) &
           *(vnucl(is,js,rij) -veq(is,js,rij)) &
           +zeta((ro+ri-2.0_rp*rij)/(ro-ri))*(dvnucl(is,js,rij) -dveq(is,js,rij))
    else if( rij.ge.ro ) then
      dvij = dveq(is,js,rij)
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
    real(rp),intent(in):: rij
    real(rp):: vnucl
    real(rp):: rs,qi,qj

    qi = qnucl(is)
    qj = qnucl(js)
!!$    rs = 0.4683766d0  /sqrt(qi**(2d0/3) +qj**(2d0/3))
!!$    rs = 0.4683766d0  /(qi**(2d0/3) +qj**(2d0/3))
    rs = 0.4683766_rp  /(qi**(0.23_rp) +qj**(0.23_rp))
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
    real(rp),intent(in):: rij
    real(rp):: dvnucl
    real(rp):: rs,qi,qj

    qi = qnucl(is)
    qj = qnucl(js)
!!$    rs = 0.4683766d0  /sqrt(qi**(2d0/3) +qj**(2d0/3))
!!$    rs = 0.4683766d0  /(qi**(2d0/3) +qj**(2d0/3))
    rs = 0.4683766_rp  /(qi**(0.23_rp) +qj**(0.23_rp))
    dvnucl = acc* qi*qj/rij* ( -1.0_rp/rij*xi(rij/rs) &
         +dxi(rij/rs)/rs )
    return
  end function dvnucl
!=======================================================================
  function veq(is,js,rij)
    implicit none
    integer,intent(in):: is,js
    real(rp),intent(in):: rij
    real(rp):: veq
    real(rp):: rk
    integer:: i

    veq = 0.0_rp
    if( rij.gt.bonny_rc(is,js) ) return
    if( is.eq.1 .and. js.eq.1 ) then  ! W-W
      veq = vspln(rij) -2.0_rp*gauge_C*rhospln(rij)
    else if( (is.eq.1 .and. js.eq.2) .or.&
         (is.eq.2 .and. js.eq.1) ) then  ! W-Re, Re-W
      do i=1,n_veq_WRe
        rk = veq_WRe_r(i)
        veq = veq +veq_WRe_a(i)*(rk -rij)**3 &
             *hvsd(rk -rij)
      enddo
    else if( is.eq.2 .and. js.eq.2 ) then  ! Re-Re
      do i=1,n_veq_ReRe
        rk = veq_ReRe_r(i)
        veq = veq +veq_ReRe_a(i)*(rk -rij)**3 &
             *hvsd(rk -rij)
      enddo
    endif
    return
  end function veq
!=======================================================================
  function dveq(is,js,rij)
    implicit none
    integer,intent(in):: is,js
    real(rp),intent(in):: rij
    real(rp):: dveq
    real(rp):: rk
    integer:: i

    dveq = 0.0_rp
    if( rij.gt.bonny_rc(is,js) ) return
    if( is.eq.1 .and. js.eq.1 ) then  ! W-W
      dveq = dvspln(rij) -2.0_rp*gauge_C*drhospln(rij)
    else if( (is.eq.1 .and. js.eq.2) .or. &
         (is.eq.2 .and. js.eq.1) ) then  ! W-Re, Re-W
      do i=1,n_veq_WRe
        rk = veq_WRe_r(i)
        dveq = dveq -veq_WRe_a(i)*(rk -rij)**2 &
             *hvsd(rk -rij)
      enddo
      dveq = dveq *3.0_rp
    else if (is.eq.2 .and. js.eq.2 ) then ! Re-Re
      do i=1,n_veq_ReRe
        rk = veq_ReRe_r(i)
        dveq = dveq -veq_ReRe_a(i)*(rk -rij)**2 &
             *hvsd(rk -rij)
      enddo
      dveq = dveq *3.0_rp
    endif
    return
  end function dveq
!=======================================================================
  function vspln(rij)
    implicit none
    real(rp),intent(in):: rij
    real(rp):: vspln,ri
    integer:: i

    vspln = 0.0_rp
    do i=1,n_vspln
      ri = vspln_r(i)
      vspln = vspln +vspln_a(i)*(ri -rij)**3 &
           *hvsd(ri -rij)
    enddo
    return
  end function vspln
!=======================================================================
  function dvspln(rij)
    implicit none
    real(rp),intent(in):: rij
    real(rp):: dvspln,ri
    integer:: i

    dvspln = 0.0_rp
    do i=1,n_vspln
      ri = vspln_r(i)
      dvspln = dvspln -vspln_a(i)*(ri -rij)**2 &
           *hvsd(ri -rij)
    enddo
    dvspln = dvspln*3.0_rp
    return
  end function dvspln
!=======================================================================
  function xi(x)
    implicit none
    real(rp),intent(in):: x
    real(rp):: xi

    xi= 0.1818_rp*exp(-3.2_rp*x) &
         +0.5099_rp*exp(-0.9423_rp*x) &
         +0.2802_rp*exp(-0.4029_rp*x) &
         +0.02817_rp*exp(-0.2016_rp*x)
    return
  end function xi
!=======================================================================
  function dxi(x)
    implicit none
    real(rp),intent(in):: x
    real(rp):: dxi

    dxi= -0.58176_rp*exp(-3.2_rp*x) &
         -0.48047877_rp*exp(-0.9423_rp*x) &
         -0.11289258_rp*exp(-0.4029_rp*x) &
         -0.005679072_rp*exp(-0.2016_rp*x)
    return
  end function dxi
!=======================================================================
  function zeta(x)
    implicit none
    real(rp),intent(in):: x
    real(rp):: zeta

    zeta = (3.0_rp*x**5 -10.0_rp*x**3 +15.0_rp*x +8.0_rp)/16.0_rp
    return
  end function zeta
!=======================================================================
  function dzeta(x)
    implicit none
    real(rp),intent(in):: x
    real(rp):: dzeta

    dzeta = (15.0_rp*x**4 -30.0_rp*x**2 +15.0_rp)/16.0_rp
    return
  end function dzeta
!=======================================================================
  subroutine set_paramsdir_Bonny(dname)
!
!  Accessor routine to set paramsdir.
!
    implicit none
    character(len=*),intent(in):: dname

    paramsdir = trim(dname)
    return
  end subroutine set_paramsdir_Bonny
!=======================================================================
  
end module Bonny_WRe
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
