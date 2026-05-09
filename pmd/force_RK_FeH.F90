module RK_FeH

  use pmdmpi
use mod_precision

contains
  subroutine force_RK_FeH(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,lstrs,iprint)
!-----------------------------------------------------------------------
!  Parallel implementation of EAM Ackland model for Fe (iron) and H.
!    - See Philos. Mag. 83(35) (2003) 3977--3994
!    - See also PRB 79, 174101 (2009) for Fe-H (Potential A)
!    - rho of boundary atoms are sent to the neighbor nodes
!    - only force on i is calculated, not necessary to send-back
!-----------------------------------------------------------------------
    implicit none
    include "./params_unit.h"
    include "params_RK_FeH.h"
    integer,intent(in):: namax,natm,nnmax,nismax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_md_world,myid_md,nex(3)
    integer,intent(in):: lspr(0:nnmax,namax)
    real(rp),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),sv(3,6) &
         ,rc,tag(namax)
    real(rp),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical:: lstrs

    integer:: i,j,k,l,m,n,ierr,is,js,ixyz,jxyz
    real(rp):: xij(3),rij,dfi,dfj,drhoij,drdxi(3),drdxj(3),at(3)
    real(rp):: x,y,z,xi(3),epotl,epott,tmp,t1,t2,t3,vemb

    logical,save:: l1st=.true.
    real(rp),allocatable,save:: rho(:)
    real(rp),save:: rs,rs_feh

    if( l1st ) then
      allocate(rho(namax))
      rs= a_rs /sqrt(2.0_rp)/z_fe**(1.0_rp/3)
      rs_feh= a_rs /sqrt(z_fe**(2.0_rp/3)+z_h**(2.0_rp/3))
!!$!.....assuming fixed (constant) atomic volume (BCC)
!!$      avol= alcfe**3 /2
!!$      if(myid_md.eq.0) write(6,'(a,es12.4)') ' avol =',avol
      l1st=.false.
!.....check cutoff radius
      if( myid_md.eq.0 ) then
        write(6,'(a,es22.14)') ' rc of input    =',rc
        write(6,'(a,es22.14)') ' rc of this pot =',rc_vphi
      endif
      if( rc.lt.rc_vphi ) then
        if( myid_md.eq.0 ) write(6,'(a)') &
             ' [get_force] rc.lt.rc_vphi !!!'
        call mpi_finalize(ierr)
        stop
      endif
    endif

    if( size(rho).lt.namax ) then
      deallocate(rho)
      allocate(rho(namax))
    endif

    epotl= 0.0_rp
    rho(:)= 0.0_rp

!.....rho(i)
    do i=1,natm
      is= int(tag(i))
      xi(1:3)= ra(1:3,i)
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
        js= int(tag(j))
        x= ra(1,j) -xi(1)
        y= ra(2,j) -xi(2)
        z= ra(3,j) -xi(3)
        xij(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
        rij=sqrt(xij(1)*xij(1)+ xij(2)*xij(2) +xij(3)*xij(3))
        if( rij.gt.rc_rho ) cycle ! Fortunately, rc_rho is common
        if( is.eq.1 .and. js.eq.1 ) then ! Fe-Fe
          rho(i)=rho(i) +fpsi(rij)
        else if( is.eq.1 .and. js.eq.2 ) then ! Fe-H
          rho(i)=rho(i) +rho_hfe(rij)
        else if( is.eq.2 .and. js.eq.1 ) then ! H-Fe
          rho(i)=rho(i) +rho_feh(rij)
        else if( is.eq.2 .and. js.eq.2 ) then ! H-H
          rho(i)=rho(i) +rho_hh(rij)
        endif
      enddo
    enddo

    call copy_dba_fwd(namax,natm,nb,nbmax,lsb,nex,&
         lsrc,myparity,nn,sv,mpi_md_world,rho,1)
!!$    if( myid_md.ge.0 ) then
!!$!.....copy rho of boundary atoms
!!$      call copy_rho_ba(namax,natm,nb,nbmax,lsb &
!!$           ,lsrc,myparity,nn,sv,mpi_md_world,rho)
!!$    else
!!$      call distribute_dba(natm,namax,tag,rho,1)
!!$    endif

!.....dE/dr_i
    do i=1,natm
      is= int(tag(i))
      xi(1:3)= ra(1:3,i)
      vemb=0.0_rp
      dfi= 0.0_rp
      if( is.eq.1 ) then ! F_{Fe}
        vemb= femb(rho(i))
        dfi= dfemb(rho(i))
      else if( is.eq.2 ) then ! Fe_{H}
        vemb= fh(rho(i))
        dfi= dfh(rho(i))
      endif
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
        if(j.le.i) cycle
        js= int(tag(j))
        x= ra(1,j) -xi(1)
        y= ra(2,j) -xi(2)
        z= ra(3,j) -xi(3)
        xij(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
        rij=sqrt(xij(1)**2+ xij(2)**2 +xij(3)**2)
        if( rij.gt.rc_vphi ) cycle
        drdxi(1:3)= -xij(1:3)/rij
!.....2-body term
        if( is.eq.1 .and. js.eq.1 ) then ! Fe-Fe
          t1= 0.5_rp *fvphi(rij,rs)
          t2= dfvphi(rij,rs)
        else if( is.ne.js ) then ! Fe-H
          if( rij.ge.r_phi_feh(7) ) cycle
          t1= 0.5_rp *fvphi_feh(rij,rs_feh)
          t2= dfvphi_feh(rij,rs_feh)
        else if( is.eq.2 .and. js.eq.2 ) then ! H-H
          if( rij.ge.rc_phi_hh ) cycle
          t1= 0.5_rp *fvphi_hh(rij)
          t2= dfvphi_hh(rij)
        endif
        epi(i)= epi(i) +t1
        epi(j)= epi(j) +t1
        if( j.le.natm ) then
          epotl=epotl +t1 +t1
        else
          epotl=epotl +t1
        endif
        aa(1:3,i)=aa(1:3,i) -drdxi(1:3)*t2
        aa(1:3,j)=aa(1:3,j) +drdxi(1:3)*t2
!.....atomic stress for 2-body terms
        do ixyz=1,3
          do jxyz=1,3
            strs(jxyz,ixyz,i)=strs(jxyz,ixyz,i) &
                 -0.5_rp*t2*xij(ixyz)*(-drdxi(jxyz))
            strs(jxyz,ixyz,j)=strs(jxyz,ixyz,j) &
                 -0.5_rp*t2*xij(ixyz)*(-drdxi(jxyz))
          enddo
        enddo
!.....Embedded term
        if( rij.gt.rc_rho ) cycle
        dfj= 0.0_rp
        tmp= 0.0_rp
        if( is.eq.1 .and. js.eq.1 ) then ! F_{Fe}, Fe-Fe
          dfj= dfemb(rho(j))
          tmp= (dfi+dfj) *dfpsi(rij)
        else if( is.eq.1 .and. js.eq.2 ) then ! F_{Fe}, Fe-H
          dfj= dfh(rho(j))
          tmp= dfi*drho_hfe(rij) +dfj*drho_feh(rij)
        else if( is.eq.2 .and. js.eq.1 ) then ! F_{H}, H-Fe
          dfj= dfemb(rho(j))
          tmp= dfi*drho_feh(rij) +dfj*drho_hfe(rij)
        else if( is.eq.2 .and. js.eq.2 ) then ! F_{H}, H-H
          dfj= dfh(rho(j))
          tmp= (dfi+dfj) *drho_hh(rij)
        endif
        aa(1:3,i)=aa(1:3,i) -drdxi(1:3)*tmp
        aa(1:3,j)=aa(1:3,j) +drdxi(1:3)*tmp
!.....atomic stress of many-body contributions
        do ixyz=1,3
          do jxyz=1,3
            strs(jxyz,ixyz,i)=strs(jxyz,ixyz,i) &
                 -0.5_rp*tmp*xij(ixyz)*(-drdxi(jxyz))
            strs(jxyz,ixyz,j)=strs(jxyz,ixyz,j) &
                 -0.5_rp*tmp*xij(ixyz)*(-drdxi(jxyz))
          enddo
        enddo
      enddo
      epi(i)=epi(i) +vemb
      epotl=epotl +vemb
    enddo

    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_md_world,strs,9)
!!$    if( myid_md.ge.0 ) then
!!$!.....copy strs of boundary atoms
!!$      call copy_dba_bk(namax,natm,nbmax,nb,lsb,lsrc,myparity &
!!$           ,nn,mpi_world,strs,9)
!!$    else
!!$      call reduce_dba_bk(natm,namax,tag,strs,9)
!!$    endif

!!$!.....atomic level stress in [eV/Ang^3] assuming 1 Ang thick
!!$    do i=1,natm
!!$      strs(1:3,1:3,i)= strs(1:3,1:3,i) /avol
!!$    enddo

!.....gather epot
    if( myid_md.ge.0 ) then
      call mpi_allreduce(epotl,epott,1,mpi_real_rp &
           ,MPI_SUM,mpi_md_world,ierr)
      epot= epot +epott
    else
      epot= epot +epotl
    endif

  end subroutine force_RK_FeH
!=======================================================================
  function fphi(x)
!
!  Screening function for r < r1
!
    implicit none
    real(rp),intent(in):: x
    real(rp):: fphi

    fphi= 0.1818_rp*exp(-3.2_rp*x) &
         +0.5099_rp*exp(-0.9423_rp*x) &
         +0.2802_rp*exp(-0.4029_rp*x) &
         +0.02817_rp*exp(-0.2016_rp*x)
    return
  end function fphi
!=======================================================================
  function dfphi(x)
!
!  1st derivative of the screening function for r < r1
!
    implicit none 
    real(rp),intent(in):: x
    real(rp):: dfphi

    dfphi= -0.58176_rp*exp(-3.2_rp*x) &
         -0.48047877_rp*exp(-0.9423_rp*x) &
         -0.11289258_rp*exp(-0.4029_rp*x) &
         -0.005679072_rp*exp(-0.2016_rp*x)
    return
  end function dfphi
!=======================================================================
  function fvphi(r,rs)
    implicit none 
    include './params_unit.h'
    include './params_RK_FeH.h'
    real(rp),intent(in):: r,rs
    real(rp):: fvphi
    integer:: i
    real(rp),external:: hvsd

    fvphi=0.0_rp
    if( r.le.r1 ) then
      fvphi=fvphi +z2q2/r*fphi(r/rs)
    else if( r.le.r2 ) then
      fvphi=fvphi +exp(b0 +b1*r +b2*r**2 +b3*r**3)
    else
      do i=2,14
        fvphi=fvphi +a_vphi(i)*(r_vphi(i)-r)**3 &
             *hvsd(r_vphi(i)-r)
      enddo
    endif
    return
  end function fvphi
!=======================================================================
  function dfvphi(r,rs)
    implicit none 
    include './params_unit.h'
    include './params_RK_FeH.h'
    real(rp),intent(in):: r,rs
    real(rp):: dfvphi
    integer:: i
    real(rp),external:: hvsd

    dfvphi= 0.0_rp
    if( r.le.r1 ) then
      dfvphi=dfvphi -z2q2/r**2*fphi(r/rs) &
           +z2q2/r/rs*dfphi(r/rs)
    else if( r.le.r2 ) then
      dfvphi=dfvphi +(b1 +2.0_rp*b2*r +3.0_rp*b3*r**2) &
           *exp(b0 +b1*r +b2*r**2 +b3*r**3)
    else
      do i=2,14
        dfvphi=dfvphi -a_vphi(i)*(r_vphi(i)-r)**2 &
             *hvsd(r_vphi(i)-r)
      enddo
      dfvphi=dfvphi*3.0_rp
    endif
    return
  end function dfvphi
!=======================================================================
  function fpsi(r)
!
!  Cubic spline function for calculating rho
!
    implicit none 
    include './params_unit.h'
    include './params_RK_FeH.h'
    real(rp),intent(in):: r
    real(rp):: fpsi
    integer:: i
    real(rp),external:: hvsd

    fpsi=0.0_rp
    do i=1,3
      fpsi=fpsi +a_psi(i)*(r_psi(i)-r)**3 &
           *hvsd(r_psi(i)-r)
    enddo
    return
  end function fpsi
!=======================================================================
  function dfpsi(r)
!
!  1st derivative of the cubic spline func of calculation of rho
!
    implicit none 
    include './params_unit.h'
    include './params_RK_FeH.h'
    real(rp),intent(in):: r
    real(rp):: dfpsi
    integer:: i
    real(rp),external:: hvsd

    dfpsi=0.0_rp
    do i=1,3
      dfpsi=dfpsi -a_psi(i)*(r_psi(i)-r)**2 &
           *hvsd(r_psi(i)-r)
    enddo
    dfpsi=dfpsi*3.0_rp
    return
  end function dfpsi
!=======================================================================
  function femb(rho)
!
!  Embedding function
!
    implicit none
    include './params_unit.h'
    include './params_RK_FeH.h'
    real(rp),intent(in):: rho
    real(rp):: femb

    femb= (-sqrt(rho) +a_emb*rho*rho)
    return
  end function femb
!=======================================================================
  function dfemb(rho)
!
!  1st derivative of the embedding function
!
    implicit none 
    include './params_unit.h'
    include './params_RK_FeH.h'
    real(rp),intent(in):: rho
    real(rp):: dfemb

    dfemb= (-0.5_rp/sqrt(rho)+2.0_rp*a_emb*rho)
    return
  end function dfemb
!=======================================================================
  function rho_feh(r)
    implicit none
    include './params_unit.h'
    include './params_RK_FeH.h'
    real(rp),intent(in):: r
    real(rp):: rho_feh
    real(rp),external:: hvsd

    integer:: i

    rho_feh= 0.0_rp
    do i=1,6
      rho_feh= rho_feh +a_rho_feh(i)*(r_rho_feh(i)-r)**3 &
           *hvsd(r_rho_feh(i)-r)
    enddo
    return
  end function rho_feh
!=======================================================================
  function drho_feh(r)
    implicit none
    include './params_unit.h'
    include './params_RK_FeH.h'
    real(rp),intent(in):: r
    real(rp):: drho_feh
    integer:: i
    real(rp),external:: hvsd

    drho_feh=0.0_rp
    do i=1,6
      drho_feh=drho_feh -a_rho_feh(i)*(r_rho_feh(i)-r)**2 &
           *hvsd(r_rho_feh(i)-r)
    enddo
    drho_feh=drho_feh*3.0_rp
    return
  end function drho_feh
!=======================================================================
  function rho_hfe(r)
    implicit none
    include './params_unit.h'
    include './params_RK_FeH.h'
    real(rp),intent(in):: r
    real(rp):: rho_hfe
    real(rp),external:: hvsd

    integer:: i

    rho_hfe= 0.0_rp
    do i=1,5
      rho_hfe= rho_hfe +a_rho_hfe(i)*(r_rho_hfe(i)-r)**3 &
           *hvsd(r_rho_hfe(i)-r)
    enddo
    return
  end function rho_hfe
!=======================================================================
  function drho_hfe(r)
    implicit none
    include './params_unit.h'
    include './params_RK_FeH.h'
    real(rp),intent(in):: r
    real(rp):: drho_hfe
    integer:: i
    real(rp),external:: hvsd

    drho_hfe=0.0_rp
    do i=1,5
      drho_hfe=drho_hfe -a_rho_hfe(i)*(r_rho_hfe(i)-r)**2 &
           *hvsd(r_rho_hfe(i)-r)
    enddo
    drho_hfe=drho_hfe*3.0_rp
    return
  end function drho_hfe
!=======================================================================
  function rho_hh(r)
    implicit none
    include './params_unit.h'
    include './params_RK_FeH.h'
    real(rp),intent(in):: r
    real(rp):: rho_hh
    
    rho_hh=0.0_rp
    if( r.ge.rc_phi_hh ) return
    rho_hh= c_rho_hh *r**2 *exp(-2.0_rp*r) *fcut(r)
    return
  end function rho_hh
!=======================================================================
  function drho_hh(r)
    implicit none
    include './params_unit.h'
    include './params_RK_FeH.h'
    real(rp),intent(in):: r
    real(rp):: drho_hh

    real(rp):: e1,fc

    drho_hh=0.0_rp
    if( r.ge.rc_phi_hh ) return
    e1= exp(-2.0_rp*r)
    fc= fcut(r)
    drho_hh= c_rho_hh *(2.0_rp*r*e1*fc -2.0_rp*r**2*e1*fc +r**2*e1*dfcut(r))
    return
  end function drho_hh
!=======================================================================
  function fcut(r)
    implicit none
    include './params_unit.h'
    include './params_RK_FeH.h'
    real(rp),intent(in):: r
    real(rp):: fcut

    fcut= exp(1.0_rp/(r-rc_phi_hh))
    return
  end function fcut
!=======================================================================
  function dfcut(r)
    implicit none
    include './params_unit.h'
    include './params_RK_FeH.h'
    real(rp),intent(in):: r
    real(rp):: dfcut

    dfcut= -1.0_rp/(r-rc_phi_hh)**2 *fcut(r)
    return
  end function dfcut
!=======================================================================
  function fh(rho)
    implicit none
    include './params_unit.h'
    include './params_RK_FeH.h'
    real(rp),intent(in):: rho
    real(rp):: fh
    integer:: i

    fh= 0.0_rp
    do i=1,6
      fh= fh +a_f(i) *rho**i
    enddo
    return
  end function fh
!=======================================================================
  function dfh(rho)
    implicit none
    include './params_unit.h'
    include './params_RK_FeH.h'
    real(rp),intent(in):: rho
    real(rp):: dfh
    integer:: i

    dfh= 0.0_rp
    do i=1,6
      dfh=dfh +a_f(i) *rho**(i-1) *i
    enddo
    return
  end function dfh
!=======================================================================
  function fvphi_hh(r)
    implicit none
    include './params_unit.h'
    include './params_RK_FeH.h'
    real(rp),intent(in):: r
    real(rp):: fvphi_hh

    real(rp):: rho,emol,adag,sr,x,ex,exm

    fvphi_hh= 0.0_rp
!.....correction term to avoid H-H clustering
    if( r.ge.r0_hh_corr .and. r.lt.r1_hh_corr ) then
      x=(r-r0_hh_corr)/lmbd_hh_corr
      ex= exp( -x**k_hh_corr )
      fvphi_hh=fvphi_hh +c0_hh_corr*x**(k_hh_corr-1.0_rp)*ex
    else if( r.ge.r1_hh_corr ) then
      x=(r-r0_hh_corr)/lmbd_hh_corr
      ex= exp( -x**k_hh_corr -b0_hh_corr*(r-r0_hh_corr)**2 )
      fvphi_hh=fvphi_hh +c0_hh_corr*x**(k_hh_corr-1.0_rp)*ex
    endif

    if( r.gt.rc_phi_hh ) return
    x=a_tanh_hh*(r-r_tanh_hh)
    ex= exp(x)
    exm=1.0_rp/ex
    sr= 0.5_rp *(1.0_rp -(ex-exm)/(ex+exm))
!      sr= 0.5d0 *(1d0 -tanh(a_tanh_hh*(r-r_tanh_hh)))
    adag= (r-r0_hh)/(r0_hh*almbd_hh)
    emol= -2.0_rp *eb_hh *(1.0_rp+adag) *exp(-adag)
    rho= rho_hh(r)

    fvphi_hh= sr *(emol -2.0_rp*fh(rho))
!      fvphi_hh= sr *(-2d0*fh(rho))

    return
  end function fvphi_hh
!=======================================================================
  function dfvphi_hh(r)
    implicit none
    include './params_unit.h'
    include './params_RK_FeH.h'
    real(rp),intent(in):: r
    real(rp):: dfvphi_hh

    real(rp):: emol,rho,adag,sr,x,ex,exm,xk1,xk2
    real(rp):: dsr,demolr,dfhr,drhor

    dfvphi_hh= 0.0_rp

!.....correction term to avoid H-H clustering
    if( r.ge.r0_hh_corr .and. r.lt.r1_hh_corr ) then
      x=(r-r0_hh_corr)/lmbd_hh_corr
      xk1= x**(k_hh_corr-1.0_rp)
      xk2= x**(k_hh_corr-2.0_rp)
      ex= exp( -x**k_hh_corr )
      dfvphi_hh=dfvphi_hh &
           +c0_hh_corr/lmbd_hh_corr*(k_hh_corr-1.0_rp) &
            *xk2 *ex &
           +c0_hh_corr*xk1*ex &
            *(-k_hh_corr/lmbd_hh_corr *xk1)
    else if( r.ge.r1_hh_corr ) then
      x=(r-r0_hh_corr)/lmbd_hh_corr
      xk1= x**(k_hh_corr-1.0_rp)
      xk2= x**(k_hh_corr-2.0_rp)
      ex= exp( -x**k_hh_corr -b0_hh_corr*(r-r0_hh_corr)**2 )
      dfvphi_hh=dfvphi_hh &
           +c0_hh_corr/lmbd_hh_corr*(k_hh_corr-1.0_rp) &
            *xk2 *ex &
           +c0_hh_corr*xk1*ex &
            *(-k_hh_corr/lmbd_hh_corr *xk1 &
            -2.0_rp*b0_hh_corr*(r-r1_hh_corr))
    endif

    if( r.gt.rc_phi_hh ) return
!      sr= 0.5d0 *(1d0 -tanh(a_tanh_hh*(r-r_tanh_hh)))
    x=a_tanh_hh*(r-r_tanh_hh)
    ex=exp(x)
    exm=1.0_rp/ex
    sr= 0.5_rp *(1.0_rp -(ex-exm)/(ex+exm))
    adag= (r-r0_hh)/(r0_hh*almbd_hh)
    emol= -2.0_rp *eb_hh *(1.0_rp+adag) *exp(-adag)
    rho= rho_hh(r)
    drhor= drho_hh(r)

!      dsr= -0.5d0 *a_tanh_hh *2d0 /cosh(a_tanh_hh*(r-r_tanh_hh))**2
!      dsr= -a_tanh_hh /cosh(a_tanh_hh*(r-r_tanh_hh))**2
    dsr= -a_tanh_hh *2.0_rp/(ex+exm)**2
    demolr= 2.0_rp *eb_hh *adag *exp(-adag) /(r0_hh*almbd_hh)
    dfhr= drhor *dfh(rho)

    dfvphi_hh= dsr*(emol-2.0_rp*fh(rho)) +sr*(demolr-2.0_rp*dfhr)
!      dfvphi_hh= dsr*(-2d0*fh(rho)) +sr*(-2d0*dfhr)
    return
  end function dfvphi_hh
!=======================================================================
  function s_hh(r)
    implicit none
    include './params_unit.h'
    include './params_RK_FeH.h'
    real(rp),intent(in):: r
    real(rp):: s_hh

    s_hh= 0.5_rp *(1.0_rp -tanh(a_tanh_hh*(r-r_tanh_hh)))
    return
  end function s_hh
!=======================================================================
  function fvphi_feh(r,rs)
    implicit none 
    include './params_unit.h'
    include './params_RK_FeH.h'
    real(rp),intent(in):: r,rs
    real(rp):: fvphi_feh
    integer:: i
    real(rp):: rr,r4
    real(rp),external:: hvsd

    fvphi_feh=0.0_rp
    if( r.le.r1_feh ) then
      fvphi_feh= z2q2_feh/r*fphi(r/rs)
    else if( r.le.r2_feh ) then
      rr=r*r
      r4=rr*rr
      fvphi_feh= (b0_feh +b1_feh*r +b2_feh*rr &
           +b3_feh*rr*r +b4_feh*r4 +b5_feh*r4*r)
!!$      fvphi_feh= exp(b0_feh +b1_feh*r +b2_feh*rr &
!!$           +b3_feh*rr*r +b4_feh*r4 +b5_feh*r4*r)
    else
      do i=1,7
        fvphi_feh=fvphi_feh +a_phi_feh(i)*(r_phi_feh(i)-r)**3 &
             *hvsd(r_phi_feh(i)-r)
      enddo
    endif
    return

  end function fvphi_feh
!=======================================================================
  function dfvphi_feh(r,rs)
    implicit none 
    include './params_unit.h'
    include './params_RK_FeH.h'
    real(rp),intent(in):: r,rs
    real(rp):: dfvphi_feh
    integer:: i
    real(rp):: r4,r3,rr
    real(rp),external:: hvsd

    dfvphi_feh= 0.0_rp
    if( r.le.r1_feh ) then
      dfvphi_feh= -z2q2_feh/r**2*fphi(r/rs) &
           +z2q2_feh/r/rs*dfphi(r/rs)
    else if( r.le.r2_feh ) then
      rr=r**2
      r3=rr*r
      r4=r3*r
      dfvphi_feh= (b1_feh +2.0_rp*b2_feh*r +3.0_rp*b3_feh*rr &
           +4.0_rp*b4_feh*r3 +5.0_rp*b5_feh*r4)
!!$      dfvphi_feh= (b1_feh +2d0*b2_feh*r +3d0*b3_feh*rr &
!!$           +4d0*b4_feh*r3 +5d0*b5_feh*r4) &
!!$           *exp(b0_feh +b1_feh*r +b2_feh*rr +b3_feh*r3 &
!!$           +b4_feh*r4 +b5_feh*r4*r)
    else
      do i=1,7
        dfvphi_feh=dfvphi_feh -a_phi_feh(i)*(r_phi_feh(i)-r)**2 &
             *hvsd(r_phi_feh(i)-r)
      enddo
      dfvphi_feh=dfvphi_feh*3.0_rp
    endif
    return
  end function dfvphi_feh
end module RK_FeH
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
