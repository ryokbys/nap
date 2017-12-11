module Bonny_WRe
!-----------------------------------------------------------------------
!                     Last modified: <2017-12-11 17:03:38 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
!  Parallel implementation of EAM poetntial of Bonney et al.
!  See G. Bonny et al., J. Appl. Phys. 121, 165107 (2017).
!-----------------------------------------------------------------------
!  Species 1 should be W and 2 should be Re.
!-----------------------------------------------------------------------
  implicit none

!.....Max num of species
  integer,parameter:: msp = 9

  logical:: interact(msp,msp)

!.....Coulomb's constant, acc = 1.0/(4*pi*epsilon0)
  real(8),parameter:: acc  = 14.3998554737d0
!.....permittivity of vacuum
  real(8),parameter:: eps0 = 0.00552634939836d0  ! e^2 /Ang /eV

  real(8):: qnucl(1:2) = (/ &
       74.0d0, & ! 1, W
       75.0d0  & ! 2, Re
       /)

  real(8):: r_inner(1:2,1:2) = reshape( (/ 1d0, 1d0, 1d0, 1d0 /), &
       shape(r_inner) )
  real(8):: r_outer(1:2,1:2) = reshape( (/ 2d0, 2d0, 2d0, 2d0 /), &
       shape(r_outer) )

  real(8),parameter:: bonny_rc(1:2,1:2) = reshape( &
       (/ 5.460d0, 3.825d0 , &  ! 1-1 W-W, 1-2 W-Re
          3.825d0, 5.460d0 /) &  ! 2-1 Re-W, 2-2 Re-Re
       , shape(bonny_rc) )

!.....Pure W parameters
  real(8),parameter:: gauge_C = 1.848055990d0
  real(8),parameter:: gauge_S = 0.2232322602d0

  real(8),parameter:: A0_W  = -5.524855802d0
  real(8),parameter:: A1_W  =  2.317313103d-1
  real(8),parameter:: A2_W  = -3.665345949d-2
  real(8),parameter:: A3_W  =  8.989367404d-3
  real(8),parameter:: rhoi_W = 1.359141225d0
  
  real(8),parameter:: rhospln_rc = 2.002970124727d0
  integer,parameter:: n_rhospln = 4
  real(8),parameter:: rhospln_a(1:4) = (/&
       -0.420429107805055d+1, &
        0.518217702261442d0,  &
        0.562720834534370d-1, &
        0.344164178842340d-1 &
        /)
  real(8),parameter:: rhospln_r(1:4) = (/&
       2.5d0, &
       3.1d0, &
       3.5d0, &
       4.9d0 &
       /)
  real(8),parameter:: fspln_a1 = -5.946454472402710d0
  real(8),parameter:: fspln_a2 = -0.049477376935239d0
  integer,parameter:: n_vspln = 15
  real(8),parameter:: vspln_a(1:15) = (/ &
        0.960851701343041d+2, &
       -0.184410923895214d+3, &
        0.935784079613550d+2, &
       -0.798358265041677d+1, &
        0.747034092936229d+1, &
       -0.152756043708453d+1, &
        0.125205932634393d+1, &
        0.163082162159425d+1, &
       -0.141854775352260d+1, &
       -0.819936046256149d0,  &
        0.198013514305908d1,  &
       -0.696430179520267d0,  &
        0.304546909722160d-1, &
       -0.163131143161660d+1, &
        0.138409896486177d+1  &
        /)
  real(8),parameter:: vspln_r(1:15) = (/ &
       2.5648975d0, &
       2.6297950d0, &
       2.6946925d0, &
       2.8663175d0, &
       2.9730450d0, &
       3.0797725d0, &
       3.5164725d0, &
       3.8464450d0, &
       4.1764175d0, &
       4.7008450d0, &
       4.8953000d0, &
       5.0897550d0, &
       5.3429525d0, &
       5.4016950d0, &
       5.4600000d0  &  ! shortened from original value 5.4604375
       /)
  
!.....Pure Re parameters
  real(8),parameter:: A_Re = -7.046791948d0
  real(8),parameter:: B_Re =  1.236584720d0
  real(8),parameter:: C_Re =  1.143405627d0
  real(8),parameter:: C0_Re=  3.704045964d-3

  integer,parameter:: n_veq_ReRe = 6
  real(8),parameter:: veq_ReRe_a(1:6) = (/ &
        6.726805309d0, &
        3.217593889d0, &
       -6.545857587d-1, &
        1.453229484d-1, &
       -2.063629464d-1, &
        6.114909116d-2  &
        /)
  real(8),parameter:: veq_ReRe_r(1:6) = (/ &
        2.700d0, &
        3.252d0, &
        3.804d0, &
        4.356d0, &
        4.908d0, &
        5.460d0  &
        /)

!.....W-Re parameters
  integer,parameter:: n_veq_WRe = 5
  real(8),parameter:: veq_WRe_a(1:5) = (/ &
       -2.335000000d+1,  &
        2.456959229d+1,  &
       -2.585878138d0,   &
        3.457586051d0,   &
       -7.013105493d-1   &
       /)
  real(8),parameter:: veq_WRe_r(1:5) = (/ &
       2.650d0, &
       2.700d0, &
       3.075d0, &
       3.450d0, &
       3.825d0  &
       /)
  
contains
  subroutine force_Bonny_WRe(namax,natm,tag,ra,nnmax,aa,strs,h,hi,tcom &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint,l1st)
    use force, only: copy_dba_fwd
    implicit none
    include "mpif.h"
    include "./params_unit.h"
    integer,intent(in):: namax,natm,nnmax,nismax,lspr(0:nnmax,namax)&
         ,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_md_world,myid_md,nex(3)
    real(8),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),sv(3,6) &
         ,acon(nismax),rc,tag(namax)
    real(8),intent(inout):: tcom
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: l1st
    logical:: lstrs

    integer:: i,j,k,l,m,n,ierr,is,js,ixyz,jxyz
    real(8):: xij(3),rij,rcij,dfi,dfj,drdxi(3),drdxj(3),r,at(3)
    real(8):: x,y,z,xi(3),epotl,epott,tmp,dtmp,drhoi,drhoj
    real(8),allocatable,save:: rho(:)
    real(8),allocatable,save:: strsl(:,:,:)

    real(8),save:: rcmax

    if( l1st ) then
      if( rc.lt.bonny_rc(1,1) ) then
        if( myid_md.eq.0 .and. iprint.gt.0 ) then
          print '(/,a)',' Input cutoff radius is smaller than rc of Bonny potential.'
          print '(a,f0.3)', '   Input rc     = ',rc
          print '(a,f0.3)', '   Potential rc = ',bonny_rc(1,1)
        endif
        call mpi_finalize(ierr)
        stop
      endif
      interact(1:msp,1:msp) = .false.
      interact(1:2,1:2) = .true.
      rcmax = 0d0
      do is=1,2
        do js=1,2
          rcmax = max(rcmax, bonny_rc(is,js))
        enddo
      enddo
      if( myid_md.eq.0 .and. iprint.gt.0 )  then
        print '(/,a,f0.3)', ' Max cutoff in Bonny potential = ',rcmax
      endif

      if( allocated(rho) ) deallocate(rho)
      allocate(rho(namax+nbmax))
      if( allocated(strsl) ) deallocate(strsl)
      allocate(strsl(3,3,namax))
    endif

    if( size(strsl).lt.3*3*namax ) then
      deallocate(strsl)
      allocate(strsl(3,3,namax))
    endif
    if( size(rho).lt.namax+nbmax ) then
      deallocate(rho)
      allocate(rho(namax+nbmax))
    endif
    
    epotl= 0d0
    rho(1:namax)= 0d0
    strsl(1:3,1:3,1:namax) = 0d0

!-----rho(i)
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is = int(tag(i))
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
        js = int(tag(j))
        if( .not. interact(is,js) ) cycle
        x= ra(1,j) -xi(1)
        y= ra(2,j) -xi(2)
        z= ra(3,j) -xi(3)
        xij(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
        rij=sqrt(xij(1)*xij(1)+ xij(2)*xij(2) +xij(3)*xij(3))
        if( rij.gt.rcmax ) cycle
        rho(i) = rho(i) +rhoij(js,rij)
      enddo
    enddo

!-----copy rho of boundary atoms
    call copy_dba_fwd(tcom,namax,natm,nb,nbmax,lsb,nex,&
         lsrc,myparity,nn,sv,mpi_md_world,rho,1)

!-----dE/dr_i
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is = int(tag(i))
      dfi = dfrho(is,rho(i))
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
        rij=sqrt(xij(1)**2+ xij(2)**2 +xij(3)**2)
        if( rij.gt.rcmax ) cycle
        drdxi(1:3)= -xij(1:3)/rij
!.....2-body term
        tmp = 0.5d0 *vij(is,js,rij)
        epi(i)= epi(i) +tmp
        epi(j)= epi(j) +tmp
        if(j.le.natm) then
          epotl=epotl +tmp +tmp
        else
          epotl=epotl +tmp
        endif
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
!.....Embedded term
        drhoi = drhoij(is,rij)
        drhoj = drhoij(js,rij)
        dfj = dfrho(js,rho(j))
        tmp = dfi*drhoj + dfj*drhoi
        aa(1:3,i)=aa(1:3,i) -tmp*drdxi(1:3)
        aa(1:3,j)=aa(1:3,j) +tmp*drdxi(1:3)
!.....Atomic stress of many-body contributions
        if( lstrs ) then
          do ixyz=1,3
            do jxyz=1,3
              strsl(jxyz,ixyz,i)=strsl(jxyz,ixyz,i) &
                   -0.5d0*tmp*xij(ixyz)*(-drdxi(jxyz))
              strsl(jxyz,ixyz,j)=strsl(jxyz,ixyz,j) &
                   -0.5d0*tmp*xij(ixyz)*(-drdxi(jxyz))
            enddo
          enddo
        endif
      enddo
      tmp = frho(is,rho(i))
      epi(i)=epi(i) +tmp
      epotl=epotl +tmp
    enddo

    if( lstrs ) then
!!$      call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
!!$           ,nn,mpi_md_world,strsl,9)
      strs(1:3,1:3,1:natm)= strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)
    endif

!-----gather epot
    call mpi_allreduce(epotl,epott,1,mpi_real8 &
         ,mpi_sum,mpi_md_world,ierr)
    epot= epot +epott

  end subroutine force_Bonny_WRe
!=======================================================================
  function rhoij(js,rij)
!
! Calculate rho of atom j at distance rij.
!
    use force, only: hvsd
    real(8),intent(in):: rij
    integer,intent(in):: js
    real(8):: rhoij

    rhoij = 0d0
    if( js.eq.1 ) then  ! Only in case of W, scaling with S
      rhoij = gauge_S *rhospln(rij)
    else if( js.eq.2 ) then
      rhoij = C0_Re *(bonny_rc(2,2) -rij)**3 *hvsd(bonny_rc(2,2) -rij)
    endif
    return
  end function rhoij
!=======================================================================
  function drhoij(js,rij)
    use force, only: hvsd
    implicit none
    real(8),intent(in):: rij
    integer,intent(in):: js
    real(8):: drhoij

    drhoij = 0d0
    if( js.eq.1 ) then
      drhoij = gauge_S *drhospln(rij)
    else if( js.eq.2 ) then
      drhoij = -3d0 *C0_Re *(bonny_rc(2,2) -rij)**2 *hvsd(bonny_rc(2,2) -rij)
    endif
    return
  end function drhoij
!=======================================================================
  function rhospln(rij)
!
!  rho_j(rij) from Marinica et al., JAP 121 (2017)
!
    use force, only: hvsd
    implicit none
    real(8),intent(in):: rij
    real(8):: rhospln,ri
    integer:: i

    rhospln = 0d0
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
    use force, only: hvsd
    implicit none
    real(8),intent(in):: rij
    real(8):: drhospln,ri
    integer:: i

    drhospln = 0d0
    if( rij.le.rhospln_rc ) then
      drhospln = 0d0
    else
      do i=1,n_rhospln
        ri = rhospln_r(i)
        drhospln = drhospln -rhospln_a(i)*(ri -rij)**2 &
             *hvsd(ri -rij)
      enddo
      drhospln = drhospln*3d0
    endif
    return
  end function drhospln
!=======================================================================
  function frho(is,rho)
    implicit none
    integer,intent(in):: is
    real(8),intent(in):: rho
    real(8):: frho

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
    real(8),intent(in):: rho
    real(8):: dfrho

    if( is.eq.1 ) then  ! W
      if( rho.le.rhoi_W ) then
        dfrho = dfeff(rho)
      else
        dfrho = A1_W +2d0*A2_W*rho +3d0*A3_W*rho*rho
      endif
    else if( is.eq.2 ) then  ! Re
      dfrho = 0.5d0*A_Re/sqrt(rho) +B_Re +2d0*C_Re*rho
    endif
    
  end function dfrho
!=======================================================================
  function feff(rho)
!
!  F^{eff} for pure W
!
    implicit none
    real(8),intent(in):: rho
    real(8):: feff

    feff = fspln(rho/gauge_S) +gauge_C/gauge_S*rho
    return
  end function feff
!=======================================================================
  function dfeff(rho)
!
!  Derivative of F^{eff} for pure W
!
    implicit none
    real(8),intent(in):: rho
    real(8):: dfeff

    dfeff = dfspln(rho/gauge_S)/gauge_S +gauge_C/gauge_S
    return
  end function dfeff
!=======================================================================
  function fspln(rho)
!
!  F[rho] function from Marinica, JAP 121, 165107 (2017)
!
    implicit none
    real(8),intent(in):: rho
    real(8):: fspln

    fspln = fspln_a1*sqrt(rho) +fspln_a2*rho*rho
    return
  end function fspln
!=======================================================================
  function dfspln(rho)
!
!  Derivative of F[rho] function from Marinica et al.
!
    implicit none
    real(8),intent(in):: rho
    real(8):: dfspln

    dfspln = 0.5d0*fspln_a1/sqrt(rho) +2d0*fspln_a2*rho
    return
  end function dfspln
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

    ri = r_inner(is,js)
    ro = r_outer(is,js)
    if( rij.lt.ri ) then
      vij = vnucl(is,js,rij)
    else if( rij.ge.ri .and. rij.lt.ro ) then
      veqt = veq(is,js,rij)
      vij = veqt +zeta((ro+ri-2d0*rij)/(ro-ri)) &
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
    real(8),intent(in):: rij
    real(8):: dvij
    real(8):: ri,ro

    ri = r_inner(is,js)
    ro = r_outer(is,js)
    if( rij.lt.ri ) then
      dvij = dvnucl(is,js,rij)
    else if( rij.ge.ri .and. rij.lt.ro ) then
      dvij = dveq(is,js,rij) +dzeta((ro+ri-2d0*rij)/(ro-ri))*(-2d0/(ro-ri)) &
           *(vnucl(is,js,rij) -veq(is,js,rij)) &
           +zeta((ro+ri-2d0*rij)/(ro-ri))*(dvnucl(is,js,rij) -dveq(is,js,rij))
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
    real(8),intent(in):: rij
    real(8):: vnucl
    real(8):: rs,qi,qj

    qi = qnucl(is)
    qj = qnucl(js)
!!$    rs = 0.4683766d0  /sqrt(qi**(2d0/3) +qj**(2d0/3))
    rs = 0.4683766d0  /(qi**(0.23d0) +qj**(0.23d0))
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
!!$    rs = 0.4683766d0  /sqrt(qi**(2d0/3) +qj**(2d0/3))
    rs = 0.4683766d0  /(qi**(0.23d0) +qj**(0.23d0))
    dvnucl = acc* qi*qj/rij* ( -1d0/rij*xi(rij/rs) &
         +dxi(rij/rs)/rs )
    return
  end function dvnucl
!=======================================================================
  function veq(is,js,rij)
    use force, only: hvsd
    implicit none
    integer,intent(in):: is,js
    real(8),intent(in):: rij
    real(8):: veq
    real(8):: rk
    integer:: i

    veq = 0d0
    if( rij.gt.bonny_rc(is,js) ) return
    if( is.eq.1 .and. js.eq.1 ) then  ! W-W
      veq = vspln(rij) -2d0*gauge_C*rhospln(rij)
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
    use force, only: hvsd
    implicit none
    integer,intent(in):: is,js
    real(8),intent(in):: rij
    real(8):: dveq
    real(8):: rk
    integer:: i

    dveq = 0d0
    if( rij.gt.bonny_rc(is,js) ) return
    if( is.eq.1 .and. js.eq.1 ) then  ! W-W
      dveq = dvspln(rij) -2d0*gauge_C*drhospln(rij)
    else if( (is.eq.1 .and. js.eq.2) .or. &
         (is.eq.2 .and. js.eq.1) ) then  ! W-Re, Re-W
      do i=1,n_veq_WRe
        rk = veq_WRe_r(i)
        dveq = dveq -veq_WRe_a(i)*(rk -rij)**2 &
             *hvsd(rk -rij)
      enddo
      dveq = dveq *3d0
    else if (is.eq.2 .and. js.eq.2 ) then ! Re-Re
      do i=1,n_veq_ReRe
        rk = veq_ReRe_r(i)
        dveq = dveq -veq_ReRe_a(i)*(rk -rij)**2 &
             *hvsd(rk -rij)
      enddo
      dveq = dveq *3d0
    endif
    return
  end function dveq
!=======================================================================
  function vspln(rij)
    use force, only: hvsd
    implicit none
    real(8),intent(in):: rij
    real(8):: vspln,ri
    integer:: i

    vspln = 0d0
    do i=1,n_vspln
      ri = vspln_r(i)
      vspln = vspln +vspln_a(i)*(ri -rij)**3 &
           *hvsd(ri -rij)
    enddo
    return
  end function vspln
!=======================================================================
  function dvspln(rij)
    use force, only: hvsd
    implicit none
    real(8),intent(in):: rij
    real(8):: dvspln,ri
    integer:: i

    dvspln = 0d0
    do i=1,n_vspln
      ri = vspln_r(i)
      dvspln = dvspln -vspln_a(i)*(ri -rij)**2 &
           *hvsd(ri -rij)
    enddo
    dvspln = dvspln*3d0
    return
  end function dvspln
!=======================================================================
  function xi(x)
    implicit none
    real(8),intent(in):: x
    real(8):: xi

    xi= 0.1818d0*exp(-3.2d0*x) &
         +0.5099d0*exp(-0.9423d0*x) &
         +0.2802d0*exp(-0.4029d0*x) &
         +0.02817d0*exp(-0.2016d0*x)
    return
  end function xi
!=======================================================================
  function dxi(x)
    implicit none
    real(8),intent(in):: x
    real(8):: dxi

    dxi= -0.58176d0*exp(-3.2d0*x) &
         -0.48047877d0*exp(-0.9423d0*x) &
         -0.11289258d0*exp(-0.4029d0*x) &
         -0.005679072d0*exp(-0.2016d0*x)
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
  
end module Bonny_WRe
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
