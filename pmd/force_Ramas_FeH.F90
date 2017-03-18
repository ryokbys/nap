module Ramas_FeH

contains
  subroutine force_Ramas_FeH(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
       ,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint) 
!-----------------------------------------------------------------------
!  Parallel implementation of EAM Ackland model for Fe (iron) and H.
!    - See Philos. Mag. 83(35) (2003) 3977--3994
!    - See also PRB 79, 174101 (2009) for Fe-H (Potential A)
!    - rho of boundary atoms are sent to the neighbor nodes
!    - only force on i is calculated, not necessary to send-back
!-----------------------------------------------------------------------
    implicit none
    include "mpif.h"
    include "./params_unit.h"
    include "params_Ramas_FeH.h"
    integer,intent(in):: namax,natm,nnmax,nismax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_md_world,myid_md,nex(3)
    integer,intent(in):: lspr(0:nnmax,namax)
    real(8),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),sv(3,6) &
         ,acon(nismax),rc,tag(namax)
    real(8),intent(inout):: tcom
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical:: lstrs

    integer:: i,j,k,l,m,n,ierr,is,js,ixyz,jxyz
    real(8):: xij(3),rij,dfi,dfj,drhoij,drdxi(3),drdxj(3),at(3)
    real(8):: x,y,z,xi(3),epotl,tmp,t1,t2,t3,vemb

    logical,save:: l1st=.true.
    real(8),allocatable,save:: rho(:)
    real(8),save:: rs,rs_feh

    if( l1st ) then
      allocate(rho(namax+nbmax))
      rs    = a_rs*a0 /(z_fe**(2d0/3) +z_fe**(2d0/3))
      rs_feh= a_rs*a0 /(z_fe**(2d0/3) +z_h**(2d0/3))
!.....assuming fixed (constant) atomic volume (BCC)
!!$      avol= alcfe**3 /2
!!$      if(myid_md.eq.0) write(6,'(a,es12.4)') ' avol =',avol
      l1st=.false.
!.....check cutoff radius
      if( myid_md.eq.0 ) then
        write(6,'(a,es22.14)') ' rc of input    =',rc
        write(6,'(a,es22.14)') ' rc of this pot =',rc_vphi
      endif
      if( rc.lt.rc_vphi ) then
        if( myid_md.eq.0 ) write(6,'(a)') ' [get_force] rc.lt.rc_vphi !!!'
        call mpi_finalize(ierr)
        stop
      endif
    endif

    aa(1:3,1:natm)=0d0
    epi(1:natm+nb)= 0d0
    epotl= 0d0
    rho(1:natm+nb)= 0d0
    strs(1:3,1:3,1:natm+nb)= 0d0

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
        rij=dsqrt(xij(1)*xij(1)+ xij(2)*xij(2) +xij(3)*xij(3))
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

    call copy_dba_fwd(tcom,namax,natm,nb,nbmax,lsb,nex,&
         lsrc,myparity,nn,sv,mpi_md_world,rho,1)
!!$    if( myid_md.ge.0 ) then
!!$!.....copy rho of boundary atoms
!!$      call copy_rho_ba(tcom,namax,natm,nb,nbmax,lsb &
!!$           ,lsrc,myparity,nn,sv,mpi_md_world,rho)
!!$    else
!!$      call distribute_dba(natm,namax,tag,rho,1)
!!$    endif

!.....dE/dr_i
    do i=1,natm
      is= int(tag(i))
      xi(1:3)= ra(1:3,i)
      vemb= 0d0
      dfi= 0d0
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
          t1= 0.5d0 *fvphi(rij,rs)
          t2= dfvphi(rij,rs)
        else if( is.ne.js ) then ! Fe-H
          if( rij.ge.r_phi_feh(7) ) cycle
          t1= 0.5d0 *fvphi_feh(rij,rs_feh)
          t2= dfvphi_feh(rij,rs_feh)
        else if( is.eq.2 .and. js.eq.2 ) then ! H-H
          if( rij.ge.rc_phi_hh ) cycle
          t1= 0.5d0 *fvphi_hh(rij)
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
                 +0.5d0*t2*xij(ixyz)*(-drdxi(jxyz))
            strs(jxyz,ixyz,j)=strs(jxyz,ixyz,j) &
                 +0.5d0*t2*xij(ixyz)*(-drdxi(jxyz))
          enddo
        enddo
!.....Embedded term
        if( rij.gt.rc_rho ) cycle
        dfj= 0d0
        tmp= 0d0
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
                 +0.5d0*tmp*xij(ixyz)*(-drdxi(jxyz))
            strs(jxyz,ixyz,j)=strs(jxyz,ixyz,j) &
                 +0.5d0*tmp*xij(ixyz)*(-drdxi(jxyz))
          enddo
        enddo
      enddo
      epi(i)=epi(i) +vemb
      epotl=epotl +vemb
    enddo

    call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_md_world,strs,9)
!!$    if( myid_md.ge.0 ) then
!!$!.....copy strs of boundary atoms
!!$      call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,lsrc,myparity &
!!$           ,nn,mpi_md_world,strs,9)
!!$    else
!!$      call reduce_dba_bk(natm,namax,tag,strs,9)
!!$    endif

!!$!.....atomic level stress in [eV/Ang^3] assuming 1 Ang thick
!!$    do i=1,natm
!!$      strs(1:3,1:3,i)= strs(1:3,1:3,i) /avol
!!$    enddo

!.....gather epot
    epot= 0d0
    if( myid_md.ge.0 ) then
      call mpi_allreduce(epotl,epot,1,MPI_DOUBLE_PRECISION &
           ,MPI_SUM,mpi_md_world,ierr)
    else
      epot= epotl
    endif

  end subroutine force_Ramas_FeH
!=======================================================================
  subroutine force_Ackland_Fe(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
       ,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr,mpi_md_world &
       ,myid_md,epi,epot,nismax,acon,lstrs,iprint) 
!-----------------------------------------------------------------------
!  Parallel implementation of EAM Ackland model for Fe (iron).
!    - See Philos. Mag. 83(35) (2003) 3977--3994
!    - rho of boundary atoms are sent to the neighbor nodes
!    - only force on i is calculated, not necessary to send-back
!-----------------------------------------------------------------------
    implicit none
    include "mpif.h"
    include "./params_unit.h"
    include "params_Ramas_FeH.h"
    integer,intent(in):: namax,natm,nnmax,nismax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_md_world,myid_md,nex(3)
    integer,intent(in):: lspr(0:nnmax,namax)
    real(8),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),sv(3,6) &
         ,acon(nismax),rc,tag(namax)
    real(8),intent(inout):: tcom
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical:: lstrs

    integer:: i,j,k,l,m,n,ierr,is,ixyz,jxyz
    real(8):: xij(3),rij,dfi,dfj,drhoij,drdxi(3),drdxj(3),at(3)
    real(8):: x,y,z,xi(3),epotl,tmp

    logical,save:: l1st=.true.
    real(8),allocatable,save:: rho(:)
    real(8),save:: rs

    if( l1st ) then
      allocate(rho(namax+nbmax))
      rs=  a_rs /dsqrt(2d0)/z_fe**(1d0/3)
!.....assuming fixed (constant) atomic volume (BCC)
!!$      avol= alcfe**3 /2
!!$      if(myid_md.eq.0) write(6,'(a,es12.4)') ' avol =',avol
      l1st=.false.
!.....check cutoff radius
      if( myid_md.eq.0 ) then
        write(6,'(a,es22.14)') ' rc of input    =',rc
        write(6,'(a,es22.14)') ' rc of this pot =',rc_vphi
      endif
      if( rc.lt.rc_vphi ) then
        if( myid_md.eq.0 ) write(6,'(a)') ' [get_force] rc.lt.rc_vphi !!!'
        call mpi_finalize(ierr)
        stop
      endif
    endif

    aa(1:3,1:natm)=0d0
    epi(1:natm)= 0d0
    epotl= 0d0
    rho(1:natm)= 0d0
    strs(1:3,1:3,1:natm+nb)= 0d0

!.....rho(i)
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
        x= ra(1,j) -xi(1)
        y= ra(2,j) -xi(2)
        z= ra(3,j) -xi(3)
        xij(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
        rij=dsqrt(xij(1)*xij(1)+ xij(2)*xij(2) +xij(3)*xij(3))
        if( rij.gt.rc_rho ) cycle
        rho(i)=rho(i) +fpsi(rij)
      enddo
    enddo

!.....copy rho of boundary atoms
    call copy_dba_fwd(tcom,namax,natm,nb,nbmax,lsb,nex,&
         lsrc,myparity,nn,sv,mpi_md_world,rho,1)
!!$    call copy_rho_ba(tcom,namax,natm,nb,nbmax,lsb &
!!$         ,lsrc,myparity,nn,sv,mpi_md_world,rho)

!.....dE/dr_i
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      dfi= dfemb(rho(i))
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
        if(j.le.i) cycle
        x= ra(1,j) -xi(1)
        y= ra(2,j) -xi(2)
        z= ra(3,j) -xi(3)
        xij(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
        rij=sqrt(xij(1)**2+ xij(2)**2 +xij(3)**2)
        if( rij.gt.rc_vphi ) cycle
        drdxi(1:3)= -xij(1:3)/rij
!.....2-body term
        tmp= 0.5d0*fvphi(rij,rs)
        epi(i)= epi(i) +tmp
        epi(j)= epi(j) +tmp
        if( j.le.natm ) then
          epotl=epotl +tmp +tmp
        else
          epotl=epotl +tmp
        endif
        tmp= dfvphi(rij,rs)
        aa(1:3,i)=aa(1:3,i) -drdxi(1:3)*tmp
        aa(1:3,j)=aa(1:3,j) +drdxi(1:3)*tmp
!.....atomic stress for 2-body terms
        do ixyz=1,3
          do jxyz=1,3
            strs(jxyz,ixyz,i)=strs(jxyz,ixyz,i) &
                 -0.5d0*tmp*xij(ixyz)*(-drdxi(jxyz))
            strs(jxyz,ixyz,j)=strs(jxyz,ixyz,j) &
                 -0.5d0*tmp*xij(ixyz)*(-drdxi(jxyz))
          enddo
        enddo
!.....Embedded term
        if( rij.gt.rc_rho ) cycle
        dfj= dfemb(rho(j))
        tmp= (dfi+dfj)*dfpsi(rij)
        aa(1:3,i)=aa(1:3,i) -drdxi(1:3)*tmp
        aa(1:3,j)=aa(1:3,j) +drdxi(1:3)*tmp
!.....atomic stress of many-body contributions
        do ixyz=1,3
          do jxyz=1,3
            strs(jxyz,ixyz,i)=strs(jxyz,ixyz,i) &
                 -0.5d0*tmp*xij(ixyz)*(-drdxi(jxyz))
            strs(jxyz,ixyz,j)=strs(jxyz,ixyz,j) &
                 -0.5d0*tmp*xij(ixyz)*(-drdxi(jxyz))
          enddo
        enddo
      enddo
      tmp= femb(rho(i))
      epi(i)=epi(i) +tmp
      epotl=epotl +tmp
    enddo

    call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_md_world,strs,9)
!!$    if( myid_md.ge.0 ) then
!!$!-----copy strs of boundary atoms
!!$      call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,lsrc,myparity &
!!$           ,nn,mpi_world,strs,9)
!!$    else
!!$      call reduce_dba_bk(natm,namax,tag,strs,9)
!!$    endif
!!$    call copy_strs_ba(tcom,namax,natm,nb,nbmax,lsb &
!!$         ,lsrc,myparity,nn,sv,mpi_md_world,strs)
!!$!.....atomic level stress in [eV/Ang^3] assuming 1 Ang thick
!!$    do i=1,natm
!!$      strs(1:3,1:3,i)= strs(1:3,1:3,i) /avol
!!$    enddo

!.....reduced force
    do i=1,natm
      at(1:3)= aa(1:3,i)
      aa(1:3,i)= hi(1:3,1)*at(1) +hi(1:3,2)*at(2) +hi(1:3,3)*at(3)
    enddo
!.....multiply 0.5d0*dt**2/am(i)
    do i=1,natm
      is= int(tag(i))
      aa(1:3,i)= acon(is)*aa(1:3,i)
    enddo

!.....gather epot
    epot= 0d0
    call mpi_allreduce(epotl,epot,1,MPI_DOUBLE_PRECISION &
         ,MPI_SUM,mpi_md_world,ierr)

  end subroutine force_Ackland_Fe
!=======================================================================
  function fphi(x)
!
!  Screening function for r < r1
!
    implicit none
    real(8),intent(in):: x
    real(8):: fphi

    fphi= 0.1818d0*exp(-3.2d0*x) &
         +0.5099d0*exp(-0.9423d0*x) &
         +0.2802d0*exp(-0.4029d0*x) &
         +0.02817d0*exp(-0.2016d0*x)
    return
  end function fphi
!=======================================================================
  function dfphi(x)
!
!  1st derivative of the screening function for r < r1
!
    implicit none 
    real(8),intent(in):: x
    real(8):: dfphi

    dfphi= -0.58176d0*exp(-3.2d0*x) &
         -0.48047877d0*exp(-0.9423d0*x) &
         -0.11289258d0*exp(-0.4029d0*x) &
         -0.005679072d0*exp(-0.2016d0*x)
    return
  end function dfphi
!=======================================================================
  function fvphi(r,rs)
    implicit none 
    include './params_unit.h'
    include './params_Ramas_FeH.h'
    real(8),intent(in):: r,rs
    real(8):: fvphi
    integer:: i
    real(8),external:: hvsd

    fvphi=0d0
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
    include './params_Ramas_FeH.h'
    real(8),intent(in):: r,rs
    real(8):: dfvphi
    integer:: i
    real(8),external:: hvsd

    dfvphi= 0d0
    if( r.le.r1 ) then
      dfvphi=dfvphi -z2q2/r**2*fphi(r/rs) +z2q2/r/rs*dfphi(r/rs)
    else if( r.le.r2 ) then
      dfvphi=dfvphi +(b1 +2d0*b2*r +3d0*b3*r**2) &
           *exp(b0 +b1*r +b2*r**2 +b3*r**3)
    else
      do i=2,14
        dfvphi=dfvphi -a_vphi(i)*(r_vphi(i)-r)**2 &
             *hvsd(r_vphi(i)-r)
      enddo
      dfvphi=dfvphi*3d0
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
    include './params_Ramas_FeH.h'
    real(8),intent(in):: r
    real(8):: fpsi
    integer:: i
    real(8),external:: hvsd

    fpsi=0d0
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
    include './params_Ramas_FeH.h'
    real(8),intent(in):: r
    real(8):: dfpsi
    integer:: i
    real(8),external:: hvsd

    dfpsi=0d0
    do i=1,3
      dfpsi=dfpsi -a_psi(i)*(r_psi(i)-r)**2 &
           *hvsd(r_psi(i)-r) 
    enddo
    dfpsi=dfpsi*3d0
    return
  end function dfpsi
!=======================================================================
  function femb(rho)
!
!  Embedding function
!
    implicit none
    include './params_unit.h'
    include './params_Ramas_FeH.h'
    real(8),intent(in):: rho
    real(8):: femb

    femb= (-dsqrt(rho) +a_emb*rho*rho)
    return
  end function femb
!=======================================================================
  function dfemb(rho)
!
!  1st derivative of the embedding function
!
    implicit none 
    include './params_unit.h'
    include './params_Ramas_FeH.h'
    real(8),intent(in):: rho
    real(8):: dfemb

    dfemb= (-0.5d0/dsqrt(rho)+2d0*a_emb*rho)
    return
  end function dfemb
!=======================================================================
  function rho_feh(r)
    implicit none
    include './params_unit.h'
    include './params_Ramas_FeH.h'
    real(8),intent(in):: r
    real(8):: rho_feh

    integer:: i
    real(8),external:: hvsd

    rho_feh= 0d0
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
    include './params_Ramas_FeH.h'
    real(8),intent(in):: r
    real(8):: drho_feh
    integer:: i
    real(8),external:: hvsd

    drho_feh=0d0
    do i=1,6
      drho_feh=drho_feh -a_rho_feh(i)*(r_rho_feh(i)-r)**2 &
           *hvsd(r_rho_feh(i)-r)
    enddo
    drho_feh=drho_feh*3d0
    return
  end function drho_feh
!=======================================================================
  function rho_hfe(r)
    implicit none
    include './params_unit.h'
    include './params_Ramas_FeH.h'
    real(8),intent(in):: r
    real(8):: rho_hfe

    integer:: i
    real(8),external:: hvsd

    rho_hfe= 0d0
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
    include './params_Ramas_FeH.h'
    real(8),intent(in):: r
    real(8):: drho_hfe
    integer:: i
    real(8),external:: hvsd

    drho_hfe=0d0
    do i=1,5
      drho_hfe=drho_hfe -a_rho_hfe(i)*(r_rho_hfe(i)-r)**2 &
           *hvsd(r_rho_hfe(i)-r)
    enddo
    drho_hfe=drho_hfe*3d0
    return
  end function drho_hfe
!=======================================================================
  function rho_hh(r)
    implicit none
    include './params_unit.h'
    include './params_Ramas_FeH.h'
    real(8),intent(in):: r
    real(8):: rho_hh

    rho_hh= 0d0
    if( r.ge.rc_phi_hh ) return
    rho_hh= c_rho_hh *r**2 *exp(-2d0*r/a0) *fcut(r)
    return
  end function rho_hh
!=======================================================================
  function drho_hh(r)
    implicit none
    include './params_unit.h'
    include './params_Ramas_FeH.h'
    real(8),intent(in):: r
    real(8):: drho_hh

    real(8):: e1,fc

    drho_hh= 0d0
    if( r.ge.rc_phi_hh ) return
    e1= exp(-2d0*r/a0)
    fc= fcut(r)
    drho_hh= c_rho_hh *(2d0*r*e1*fc &
         -2d0/a0*r**2*e1*fc &
         +r**2*e1*dfcut(r))
    return
  end function drho_hh
!=======================================================================
  function fcut(r)
    implicit none
    include './params_unit.h'
    include './params_Ramas_FeH.h'
    real(8),intent(in):: r
    real(8):: fcut

    fcut= exp(1d0/(r-rc_phi_hh))
    return
  end function fcut
!=======================================================================
  function dfcut(r)
    implicit none
    include './params_unit.h'
    include './params_Ramas_FeH.h'
    real(8),intent(in):: r
    real(8):: dfcut

    dfcut= -1d0/(r-rc_phi_hh)**2 *fcut(r)
    return
  end function dfcut
!=======================================================================
  function fh(rho)
    implicit none
    include './params_unit.h'
    include './params_Ramas_FeH.h'
    real(8),intent(in):: rho
    real(8):: fh
    integer:: i

    fh= 0d0
    do i=1,6
      fh= fh +a_f(i) *rho**i
    enddo
    return
  end function fh
!=======================================================================
  function dfh(rho)
    implicit none
    include './params_unit.h'
    include './params_Ramas_FeH.h'
    real(8),intent(in):: rho
    real(8):: dfh
    integer:: i

    dfh= 0d0
    do i=1,6
      dfh=dfh +a_f(i) *rho**(i-1) *i
    enddo
    return
  end function dfh
!=======================================================================
  function fvphi_hh(r)
    implicit none
    include './params_unit.h'
    include './params_Ramas_FeH.h'
    real(8),intent(in):: r
    real(8):: fvphi_hh
    real(8):: rho,emol,adag,sr,x,ex,exm

    fvphi_hh= 0d0

!!$!.....correction term to avoid H-H clustering
!!$    if( r.ge.r0_hh_corr .and. r.lt.r1_hh_corr ) then
!!$      x=(r-r0_hh_corr)/lmbd_hh_corr
!!$      ex= exp( -x**k_hh_corr )
!!$      fvphi_hh=fvphi_hh +c0_hh_corr*x**(k_hh_corr-1d0)*ex
!!$    else if( r.ge.r1_hh_corr ) then
!!$      x=(r-r0_hh_corr)/lmbd_hh_corr
!!$      ex= exp( -x**k_hh_corr -b0_hh_corr*(r-r0_hh_corr)**2 )
!!$      fvphi_hh=fvphi_hh +c0_hh_corr*x**(k_hh_corr-1d0)*ex
!!$    endif

    if( r.gt.rc_phi_hh ) return
    x=a_tanh_hh*(r-r_tanh_hh)
    ex= exp(x)
    exm=1d0/ex
    sr= 0.5d0 *(1d0 -(ex-exm)/(ex+exm))
!      sr= 0.5d0 *(1d0 -tanh(a_tanh_hh*(r-r_tanh_hh)))
    adag= (r-r0_hh)/(r0_hh*almbd_hh)
    emol= -2d0 *eb_hh *(1d0+adag) *exp(-adag)
    rho= rho_hh(r)

!.....because C1,C2 are 0, not calculate their term
    fvphi_hh= fvphi_hh +sr *(emol -2d0*fh(rho))

    return
  end function fvphi_hh
!=======================================================================
  function dfvphi_hh(r)
    implicit none
    include './params_unit.h'
    include './params_Ramas_FeH.h'
    real(8),intent(in):: r
    real(8):: dfvphi_hh

    real(8):: emol,rho,adag,sr,x,ex,exm,xk1,xk2
    real(8):: dsr,demolr,dfhr,drhor

    dfvphi_hh= 0d0

!!$!.....correction term to avoid H-H clustering
!!$    if( r.ge.r0_hh_corr .and. r.lt.r1_hh_corr ) then
!!$      x=(r-r0_hh_corr)/lmbd_hh_corr
!!$      xk1= x**(k_hh_corr-1d0)
!!$      xk2= x**(k_hh_corr-2d0)
!!$      ex= exp( -x**k_hh_corr )
!!$      dfvphi_hh=dfvphi_hh &
!!$           +c0_hh_corr/lmbd_hh_corr*(k_hh_corr-1d0) &
!!$            *xk2 *ex &
!!$           +c0_hh_corr*xk1*ex &
!!$            *(-k_hh_corr/lmbd_hh_corr *xk1)
!!$    else if( r.ge.r1_hh_corr ) then
!!$      x=(r-r0_hh_corr)/lmbd_hh_corr
!!$      xk1= x**(k_hh_corr-1d0)
!!$      xk2= x**(k_hh_corr-2d0)
!!$      ex= exp( -x**k_hh_corr -b0_hh_corr*(r-r0_hh_corr)**2 )
!!$      dfvphi_hh=dfvphi_hh &
!!$           +c0_hh_corr/lmbd_hh_corr*(k_hh_corr-1d0) &
!!$            *xk2 *ex &
!!$           +c0_hh_corr*xk1*ex &
!!$            *(-k_hh_corr/lmbd_hh_corr *xk1 &
!!$            -2d0*b0_hh_corr*(r-r1_hh_corr))
!!$    endif


    if( r.gt.rc_phi_hh ) return
!      sr= 0.5d0 *(1d0 -tanh(a_tanh_hh*(r-r_tanh_hh)))
    x=a_tanh_hh*(r-r_tanh_hh)
    ex=exp(x)
    exm=1d0/ex
    sr= 0.5d0 *(1d0 -(ex-exm)/(ex+exm))
    adag= (r-r0_hh)/(r0_hh*almbd_hh)
    emol= -2d0 *eb_hh *(1d0+adag) *exp(-adag)
    rho= rho_hh(r)
    drhor= drho_hh(r)

!      dsr= -0.5d0 *a_tanh_hh *2d0 /cosh(a_tanh_hh*(r-r_tanh_hh))**2
!      dsr= -a_tanh_hh /cosh(a_tanh_hh*(r-r_tanh_hh))**2
    dsr= -a_tanh_hh *2d0/(ex+exm)**2
    demolr= 2d0 *eb_hh *adag *exp(-adag) /(r0_hh*almbd_hh)
    dfhr= drhor *dfh(rho)

    dfvphi_hh=dfvphi_hh +dsr*(emol-2d0*fh(rho)) +sr*(demolr-2d0*dfhr)
!      dfvphi_hh= dsr*(-2d0*fh(rho)) +sr*(-2d0*dfhr)
    return
  end function dfvphi_hh
!=======================================================================
  function s_hh(r)
    implicit none
    include './params_unit.h'
    include './params_Ramas_FeH.h'
    real(8),intent(in):: r
    real(8):: s_hh

    s_hh= 0.5d0 *(1d0 -tanh(a_tanh_hh*(r-r_tanh_hh)))
    return
  end function s_hh
!=======================================================================
  function fvphi_feh(r,rs)
    implicit none 
    include './params_unit.h'
    include './params_Ramas_FeH.h'
    real(8),intent(in):: r,rs
    real(8):: fvphi_feh
    integer:: i
    real(8),external:: hvsd
    real(8):: rr,r4

    fvphi_feh=0d0
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
    include './params_Ramas_FeH.h'
    real(8),intent(in):: r,rs
    real(8):: dfvphi_feh
    integer:: i
    real(8),external:: hvsd
    real(8):: r4,r3,rr

    dfvphi_feh= 0d0
    if( r.le.r1_feh ) then
      dfvphi_feh= -z2q2_feh/r**2*fphi(r/rs) &
           +z2q2_feh/r/rs*dfphi(r/rs)
    else if( r.le.r2_feh ) then
      rr=r**2
      r3=rr*r
      r4=r3*r
      dfvphi_feh= (b1_feh +2d0*b2_feh*r +3d0*b3_feh*rr &
           +4d0*b4_feh*r3 +5d0*b5_feh*r4)
!!$      dfvphi_feh= (b1_feh +2d0*b2_feh*r +3d0*b3_feh*rr &
!!$           +4d0*b4_feh*r3 +5d0*b5_feh*r4) &
!!$           *exp(b0_feh +b1_feh*r +b2_feh*rr +b3_feh*r3 &
!!$           +b4_feh*r4 +b5_feh*r4*r)
    else
      do i=1,7
        dfvphi_feh=dfvphi_feh -a_phi_feh(i)*(r_phi_feh(i)-r)**2 &
             *hvsd(r_phi_feh(i)-r)
      enddo
      dfvphi_feh=dfvphi_feh*3d0
    endif
    return
  end function dfvphi_feh
end module Ramas_FeH
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
