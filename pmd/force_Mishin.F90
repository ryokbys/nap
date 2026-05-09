module Mishin
!-----------------------------------------------------------------------
!  Parallel implementation of force and potential energy calculation
!  of EAM for Al and Ni by Mishin et al.
!    - rho of boundary atoms are sent to the neighbor nodes
!    - only force on i is calculated, not necessary to send-back
!-----------------------------------------------------------------------
!  See Mishin et al. PRB 59(5) (1999) 3393--3407.
!-----------------------------------------------------------------------
  use pmdmpi
  use mod_precision
  implicit none
  save

  real(rp):: rc_eam,rcmax2
  integer,parameter:: neamd = 25

  real(rp):: rtbl(neamd)
  real(rp):: rhotbl(neamd)
  real(rp):: rhoprm(4,neamd-1)
  real(rp):: vprm(4,neamd-1)
  real(rp):: fprm(4,neamd-1)
  
contains
  subroutine force_Mishin_Al(namax,natm,tag_isp,ra,nnmax,aa,strs,h,hi &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,lstrs,iprint,l1st)
    implicit none
    include "./params_unit.h"
    include "params_Mishin_Al.h"
    integer,intent(in):: namax,natm,nnmax,nismax,lspr(0:nnmax,namax)&
         ,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_md_world,myid_md,nex(3)
    integer,intent(in):: tag_isp(namax)
    real(rp),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),sv(3,6) &
         ,rc
    real(rp),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: lstrs,l1st

    if( l1st ) then
      rtbl(1:neamd) = rtbl_Al(1:neamd)
      rhotbl(1:neamd) = rhotbl_Al(1:neamd)
      rhoprm(1:4,1:neamd-1) = rhoprm_Al(1:4,1:neamd-1)
      vprm(1:4,1:neamd-1) = vprm_Al(1:4,1:neamd-1)
      fprm(1:4,1:neamd-1) = fprm_Al(1:4,1:neamd-1)
      rc_eam = rc_eam_Al
      rcmax2 = rc_eam*rc_eam
      if( myid_md.eq.0 ) then
        print *,''
        write(6,'(a)') ' force_Mishin_Al:'
        write(6,'(a,es12.4)') '   rc of input    =',rc
        write(6,'(a,es12.4)') '   rc of this pot =',rc_eam
      endif
      if( rc.lt.rc_eam ) then
        if( myid_md.eq.0 ) print *,'ERROR: rc < rc_eam !!!'
        stop 1
      endif
    endif

    call force_Mishin(namax,natm,tag_isp,ra,nnmax,aa,strs,h,hi &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,lstrs,iprint,l1st)
    return
  end subroutine force_Mishin_Al
!=======================================================================
  subroutine force_Mishin_Ni(namax,natm,tag_isp,ra,nnmax,aa,strs,h,hi &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,lstrs,iprint,l1st)
    implicit none
    include "./params_unit.h"
    include "params_Mishin_Ni.h"
    integer,intent(in):: namax,natm,nnmax,nismax,lspr(0:nnmax,namax)&
         ,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_md_world,myid_md,nex(3)
    integer,intent(in):: tag_isp(namax)
    real(rp),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),sv(3,6) &
         ,rc
    real(rp),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: lstrs,l1st

    if( l1st ) then
      rtbl(1:neamd) = rtbl_Ni(1:neamd)
      rhotbl(1:neamd) = rhotbl_Ni(1:neamd)
      rhoprm(1:4,1:neamd-1) = rhoprm_Ni(1:4,1:neamd-1)
      vprm(1:4,1:neamd-1) = vprm_Ni(1:4,1:neamd-1)
      fprm(1:4,1:neamd-1) = fprm_Ni(1:4,1:neamd-1)
      rc_eam = rc_eam_Ni
      rcmax2 = rc_eam*rc_eam
      if( myid_md.eq.0 ) then
        print *,''
        write(6,'(a)') ' force_Mishin_Ni:'
        write(6,'(a,es12.4)') '   rc of input    =',rc
        write(6,'(a,es12.4)') '   rc of this pot =',rc_eam
      endif
      if( rc.lt.rc_eam ) then
        if( myid_md.eq.0 ) print *,'ERROR: rc < rc_eam !!!'
        stop 1
      endif
    endif

    call force_Mishin(namax,natm,tag_isp,ra,nnmax,aa,strs,h,hi &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,lstrs,iprint,l1st)
    return
  end subroutine force_Mishin_Ni
!=======================================================================
  subroutine force_Mishin(namax,natm,tag_isp,ra,nnmax,aa,strs,h,hi &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,lstrs,iprint,l1st)
    implicit none
    include "./params_unit.h"
    include "params_Mishin_Al.h"
    integer,intent(in):: namax,natm,nnmax,nismax,lspr(0:nnmax,namax)&
         ,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_md_world,myid_md,nex(3)
    integer,intent(in):: tag_isp(namax)
    real(rp),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),sv(3,6) &
         ,rc
    real(rp),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: lstrs,l1st

    integer:: i,j,k,l,m,n,ierr,is,ixyz,jxyz
    real(rp):: xij(3),rij(3),dij,dij2,dfi,dfj,drhoij,drdxi(3),drdxj(3) &
         ,r,dphi,at(3),xi(3),xj(3),epotl,epott,tmp,sji
!.....Saved variables
    real(rp),allocatable,save:: rho(:)
    real(rp),allocatable,save:: strsl(:,:,:)

    if( l1st ) then
      if( allocated(rho) ) deallocate(rho)
      allocate(rho(namax))
      if( allocated(strsl) ) deallocate(strsl)
      allocate(strsl(3,3,namax))
    endif

    if( size(rho).lt.namax ) then
      deallocate(rho,strsl)
      allocate(rho(namax),strsl(3,3,namax))
    endif

    epotl= 0.0_rp
    rho(1:natm)= 0.0_rp
    strsl(1:3,1:3,1:natm) = 0.0_rp

!-----rho(i)
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
        xj(1:3)= ra(1:3,j)
        xij(1:3)= xj(1:3) -xi(1:3)
        rij(1:3)= h(1:3,1,0)*xij(1) +h(1:3,2,0)*xij(2) +h(1:3,3,0)*xij(3)
        dij2= rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
        if( dij2.gt.rcmax2 ) cycle
        dij= sqrt(dij2)
        rho(i)= rho(i) +calc_rho(dij)
      enddo
    enddo

    call copy_dba_fwd(namax,natm,nb,nbmax,lsb,nex,&
         lsrc,myparity,nn,sv,mpi_md_world,rho,1)

!-----dE/dr_i
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      dfi= calc_df(rho(i))
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
        if(j.le.i) cycle
        xj(1:3)= ra(1:3,j)
        xij(1:3)= xj(1:3) -xi(1:3)
        rij(1:3)= h(1:3,1,0)*xij(1) +h(1:3,2,0)*xij(2) +h(1:3,3,0)*xij(3)
        dij2= rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
        if( dij2.gt.rcmax2 ) cycle
        dij= sqrt(dij2)
        drdxi(1:3)= -rij(1:3)/dij
!.....2-body term
        tmp= 0.5_rp *calc_v(dij)
        epi(i)= epi(i) +tmp
        epi(j)= epi(j) +tmp
        if(j.le.natm) then
          epotl=epotl +tmp +tmp
        else
          epotl=epotl +tmp
        endif
        dphi= calc_dv(dij)
        aa(1:3,i)=aa(1:3,i) -dphi*drdxi(1:3)
        aa(1:3,j)=aa(1:3,j) +dphi*drdxi(1:3)
!.....Atomic stress of 2-body terms
        do ixyz=1,3
          do jxyz=1,3
            sji = -dphi*drdxi(ixyz)*(-rij(jxyz)) *0.5_rp
            strsl(jxyz,ixyz,i)=strsl(jxyz,ixyz,i) +sji
            strsl(jxyz,ixyz,j)=strsl(jxyz,ixyz,j) +sji
          enddo
        enddo
!.....Embeded term
        drhoij= calc_drho(dij)
        dfj= calc_df(rho(j))
        tmp = (dfi+dfj)*drhoij
        aa(1:3,i)=aa(1:3,i) -tmp*drdxi(1:3)
        aa(1:3,j)=aa(1:3,j) +tmp*drdxi(1:3)
!.....Atomic stress of many-body contributions
        do ixyz=1,3
          do jxyz=1,3
            sji = -tmp*drdxi(ixyz)*(-rij(jxyz)) *0.5_rp
            strsl(jxyz,ixyz,i)=strsl(jxyz,ixyz,i) +sji
            strsl(jxyz,ixyz,j)=strsl(jxyz,ixyz,j) +sji
          enddo
        enddo
      enddo
      tmp= calc_f(rho(i))
      epi(i)=epi(i) +tmp
      epotl=epotl +tmp
    enddo

    strs(1:3,1:3,1:natm)= strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)

!-----gather epot
    if( myid_md.ge.0 ) then
      call mpi_allreduce(epotl,epott,1,mpi_real_rp &
           ,MPI_SUM,mpi_md_world,ierr)
      epot= epot +epott
    else
      epot= epot +epotl
    endif

!      deallocate(sqrho)
  end subroutine force_Mishin
!=======================================================================
  function calc_rho(r)
    implicit none
    include './params_unit.h'
    include 'params_Mishin_Al.h'
    real(rp),intent(in):: r
    real(rp):: calc_rho

    integer:: i
    real(rp):: a(4),r0,rho0,drho0

!!$    call check_range(r,neamd,rtbl,'calc_rho')
    if( r.ge.rtbl(neamd) ) then
      calc_rho= 0.0_rp
      return
    else if( r.le.rtbl(1) ) then
      a(1:4)= rhoprm(1:4,1)
      r0 = rtbl(1)
      rho0 = a(1) +a(2)*r0 +a(3)*r0*r0 +a(4)*r0*r0*r0
      drho0 = a(2) +2.0_rp*a(3)*r0 +3.0_rp*a(4)*r0*r0
      calc_rho = rho0 +(r-r0)*drho0
      return
    endif
    do i=1,neamd-1
      if( r.lt.rtbl(i+1) ) then
        a(1:4)= rhoprm(1:4,i)
        exit
      endif
    enddo

    calc_rho= a(1) +a(2)*r +a(3)*r*r +a(4)*r*r*r
    return
  end function calc_rho
!=======================================================================
  function calc_drho(r)
    implicit none
    include './params_unit.h'
    include 'params_Mishin_Al.h'
    real(rp),intent(in):: r
    real(rp):: calc_drho

    integer:: i
    real(rp):: a(4),r0

    if( r.ge.rtbl(neamd) ) then
      calc_drho= 0.0_rp
      return
    else if( r.le.rtbl(1) ) then
      a(1:4)= rhoprm(1:4,1)
      r0 = rtbl(1)
      calc_drho = a(2) +2.0_rp*a(3)*r0 +3.0_rp*a(4)*r0*r0
      return
    endif
    do i=1,neamd-1
      if( r.lt.rtbl(i+1) ) then
        a(1:4)= rhoprm(1:4,i)
        exit
      endif
    enddo

    calc_drho= a(2) +2.0_rp*a(3)*r +3.0_rp*a(4)*r*r
    return
  end function calc_drho
!=======================================================================
  function calc_v(r)
    implicit none
    include './params_unit.h'
    include 'params_Mishin_Al.h'
    real(rp),intent(in):: r
    real(rp):: calc_v

    integer:: i
    real(rp):: a(4),r0,v0,dv0

    if( r.ge.rtbl(neamd) ) then
      calc_v= 0.0_rp
      return
    else if( r.le.rtbl(1) ) then
      a(1:4)= vprm(1:4,1)
      r0 = rtbl(1)
      v0 = a(1) +a(2)*r0 +a(3)*r0*r0 +a(4)*r0*r0*r0
      dv0 = a(2) +2.0_rp*a(3)*r0 +3.0_rp*a(4)*r0*r0
      calc_v= v0 +(r -r0)*dv0
      return
    endif
    do i=1,neamd-1
      if( r.lt.rtbl(i+1) ) then
        a(1:4)= vprm(1:4,i)
        exit
      endif
    enddo

    calc_v= a(1) +a(2)*r +a(3)*r*r +a(4)*r*r*r
    return
  end function calc_v
!=======================================================================
  function calc_dv(r)
    implicit none
    include './params_unit.h'
    include 'params_Mishin_Al.h'
    real(rp),intent(in):: r
    real(rp):: calc_dv

    integer:: i
    real(rp):: a(4),r0

    if( r.ge.rtbl(neamd) ) then
      calc_dv= 0.0_rp
      return
    else if( r.le.rtbl(1) ) then
      a(1:4)= vprm(1:4,1)
      r0 = rtbl(1)
      calc_dv = a(2) +2.0_rp*a(3)*r0 +3.0_rp*a(4)*r0*r0
      return
    endif
    do i=1,neamd-1
      if( r.lt.rtbl(i+1) ) then
        a(1:4)= vprm(1:4,i)
        exit
      endif
    enddo

    calc_dv= a(2) +2.0_rp*a(3)*r +3.0_rp*a(4)*r*r
    return
  end function calc_dv
!=======================================================================
  function calc_f(rho)
    implicit none
    include './params_unit.h'
    include 'params_Mishin_Al.h'
    real(rp),intent(in):: rho
    real(rp):: calc_f

    integer:: i
    real(rp):: a(4),rho0,f0,df0

    if( rho.ge.rhotbl(neamd) ) then
      a(1:4) = fprm(1:4,neamd-1)
    else if( rho.le.rhotbl(1) ) then
      a(1:4)= fprm(1:4,1)
      rho0 = rhotbl(1)
      f0 = a(1) +a(2)*rho0 +a(3)*rho0*rho0 +a(4)*rho0*rho0*rho0
      df0 = a(2) +2.0_rp*a(3)*rho0 +3.0_rp*a(4)*rho0*rho0
      calc_f = f0 +(rho -rho0)*df0
    else
      do i=1,neamd-1
        if( rho.lt.rhotbl(i+1) ) then
          a(1:4)= fprm(1:4,i)
          exit
        endif
      enddo
    endif

    calc_f= a(1) +a(2)*rho +a(3)*rho*rho +a(4)*rho*rho*rho
    return
  end function calc_f
!=======================================================================
  function calc_df(rho)
    implicit none
    include './params_unit.h'
    include 'params_Mishin_Al.h'
    real(rp),intent(in):: rho
    real(rp):: calc_df

    integer:: i
    real(rp):: a(4),rho0

    if( rho.ge.rhotbl(neamd) ) then
      a(1:4)= fprm(1:4,neamd-1)
    else if( rho.le.rhotbl(1) ) then
      a(1:4)= fprm(1:4,1)
      rho0 = rhotbl(1)
      calc_df = a(2) +2.0_rp*a(3)*rho0 +3.0_rp*a(4)*rho0*rho0
    else
      do i=1,neamd-1
        if( rho.lt.rhotbl(i+1) ) then
          a(1:4)= fprm(1:4,i)
          exit
        endif
      enddo
    endif

    calc_df= a(2) +2.0_rp*a(3)*rho +3.0_rp*a(4)*rho*rho
    return
  end function calc_df
!=======================================================================
  subroutine check_range(x,n,tbl,c)
    integer,intent(in):: n
    real(rp),intent(in):: x,tbl(n)
    character(len=*),intent(in):: c

    if( x.lt.tbl(1) ) then
      write(6,'(2es12.4)') x,tbl(1)
      write(6,'(a)') ' ['//trim(c)//'] x.lt.tbl(1)'
      stop
    endif

  end subroutine check_range
!=======================================================================
end module Mishin
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
