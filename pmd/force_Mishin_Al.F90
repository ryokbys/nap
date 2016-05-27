module Mishin_Al
contains
  subroutine force_Mishin_Al(namax,natm,tag,ra,nnmax,aa,strs,h,hi,tcom &
       ,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs)
!-----------------------------------------------------------------------
!  Parallel implementation of force and potential energy calculation
!  of EAM for Al by Mishin et al.
!    - rho of boundary atoms are sent to the neighbor nodes
!    - only force on i is calculated, not necessary to send-back
!-----------------------------------------------------------------------
!  See Mishin et al. PRB 59(5) (1999) 3393--3407.
!-----------------------------------------------------------------------
    implicit none
    include "mpif.h"
    include "./params_unit.h"
    include "params_Mishin_Al.h"
    integer,intent(in):: namax,natm,nnmax,nismax,lspr(0:nnmax,namax)
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_md_world,myid_md
    real(8),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),sv(3,6) &
         ,acon(nismax),rc,tag(namax)
    real(8),intent(inout):: tcom
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical:: lstrs

    integer:: i,j,k,l,m,n,ierr,is
    real(8):: xij(3),rij,dfi,dfj,drhoij,drdxi(3),drdxj(3),r,dphi,at(3)
    real(8):: x,y,z,xi(3),epotl,tmp
!.....Saved variables
    logical,save:: l1st=.true.
    real(8),allocatable,save:: rho(:)

    if( l1st ) then
      allocate(rho(namax+nbmax))
      l1st=.false.
    endif

    aa(1:3,1:natm)=0d0
    epi(1:natm)= 0d0
    epotl= 0d0
    rho(1:natm)= 0d0
    strs(1:3,1:3,natm+nb)= 0d0

!-----rho(i)
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
        x= ra(1,j) -xi(1)
        y= ra(2,j) -xi(2)
        z= ra(3,j) -xi(3)
        xij(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
        rij=sqrt(xij(1)*xij(1)+ xij(2)*xij(2) +xij(3)*xij(3))
        if( rij.gt.rc_eam ) cycle
        rho(i)= rho(i) +calc_rho(rij)
      enddo
    enddo

    if( myid_md.ge.0 ) then
!-----copy rho of boundary atoms
      call copy_rho_ba(tcom,namax,natm,nb,nbmax,lsb,lsrc,myparity,nn,sv &
           ,mpi_md_world,rho)
    else
      call distribute_dba(natm,namax,tag,rho,1)
    endif

!-----dE/dr_i
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      dfi= calc_df(rho(i))
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
        if(j.le.i) cycle
        x= ra(1,j) -xi(1)
        y= ra(2,j) -xi(2)
        z= ra(3,j) -xi(3)
        xij(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
        rij=sqrt(xij(1)**2+ xij(2)**2 +xij(3)**2)
        if( rij.gt.rc_eam ) cycle
        drdxi(1:3)= -xij(1:3)/rij
!.....2-body term
        tmp= 0.5d0 *calc_v(rij)
        epi(i)= epi(i) +tmp
        epi(j)= epi(j) +tmp
        if(j.le.natm) then
          epotl=epotl +tmp +tmp
        else
          epotl=epotl +tmp
        endif
        dphi= calc_dv(rij)
        aa(1:3,i)=aa(1:3,i) -dphi*drdxi(1:3)
        aa(1:3,j)=aa(1:3,j) +dphi*drdxi(1:3)
!.....Embeded term
        drhoij= calc_drho(rij)
        dfj= calc_df(rho(j))
        aa(1:3,i)=aa(1:3,i) -(dfi+dfj)*drhoij*drdxi(1:3)
        aa(1:3,j)=aa(1:3,j) +(dfi+dfj)*drhoij*drdxi(1:3)
      enddo
      tmp= calc_f(rho(i))
      epi(i)=epi(i) +tmp
      epotl=epotl +tmp
    enddo

!-----reduced force
    do i=1,natm
      at(1:3)= aa(1:3,i)
      aa(1:3,i)= hi(1:3,1)*at(1) +hi(1:3,2)*at(2) +hi(1:3,3)*at(3)
    enddo
!-----multiply 0.5d0*dt**2/am(i)
    do i=1,natm
      is= int(tag(i))
      aa(1:3,i)= acon(is)*aa(1:3,i)
    enddo

!-----gather epot
    epot= 0d0
    if( myid_md.ge.0 ) then
      call mpi_allreduce(epotl,epot,1,MPI_DOUBLE_PRECISION &
           ,MPI_SUM,mpi_md_world,ierr)
    else
      epot= epotl
    endif

!      deallocate(sqrho)
  end subroutine force_Mishin_Al
!=======================================================================
  function calc_rho(r)
    implicit none
    include './params_unit.h'
    include 'params_Mishin_Al.h'
    real(8),intent(in):: r
    real(8):: calc_rho

    integer:: i
    real(8):: a(4)

    call check_range(r,neamd,rtbl,'calc_rho')
    if( r.ge.rtbl(neamd) ) then
      calc_rho= 0d0
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
    real(8),intent(in):: r
    real(8):: calc_drho

    integer:: i
    real(8):: a(4)

    if( r.ge.rtbl(neamd) ) then
      calc_drho= 0d0
      return
    endif
    do i=1,neamd-1
      if( r.lt.rtbl(i+1) ) then
        a(1:4)= rhoprm(1:4,i)
        exit
      endif
    enddo

    calc_drho= a(2) +2d0*a(3)*r +3d0*a(4)*r*r
    return
  end function calc_drho
!=======================================================================
  function calc_v(r)
    implicit none
    include './params_unit.h'
    include 'params_Mishin_Al.h'
    real(8),intent(in):: r
    real(8):: calc_v

    integer:: i
    real(8):: a(4)

    if( r.ge.rtbl(neamd) ) then
      calc_v= 0d0
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
    real(8),intent(in):: r
    real(8):: calc_dv

    integer:: i
    real(8):: a(4)

    if( r.ge.rtbl(neamd) ) then
      calc_dv= 0d0
      return
    endif
    do i=1,neamd-1
      if( r.lt.rtbl(i+1) ) then
        a(1:4)= vprm(1:4,i)
        exit
      endif
    enddo

    calc_dv= a(2) +2d0*a(3)*r +3d0*a(4)*r*r
    return
  end function calc_dv
!=======================================================================
  function calc_f(rho)
    implicit none
    include './params_unit.h'
    include 'params_Mishin_Al.h'
    real(8),intent(in):: rho
    real(8):: calc_f

    integer:: i
    real(8):: a(4)

    if( rho.ge.rhotbl(neamd) ) then
      calc_f= 0d0
      return
    endif
    do i=1,neamd-1
      if( rho.lt.rhotbl(i+1) ) then
        a(1:4)= fprm(1:4,i)
        exit
      endif
    enddo

    calc_f= a(1) +a(2)*rho +a(3)*rho*rho +a(4)*rho*rho*rho
    return
  end function calc_f
!=======================================================================
  function calc_df(rho)
    implicit none
    include './params_unit.h'
    include 'params_Mishin_Al.h'
    real(8),intent(in):: rho
    real(8):: calc_df

    integer:: i
    real(8):: a(4)

    if( rho.ge.rhotbl(neamd) ) then
      calc_df= 0d0
      return
    endif
    do i=1,neamd-1
      if( rho.lt.rhotbl(i+1) ) then
        a(1:4)= fprm(1:4,i)
        exit
      endif
    enddo

    calc_df= a(2) +2d0*a(3)*rho +3d0*a(4)*rho*rho
    return
  end function calc_df
!=======================================================================
  subroutine check_range(x,n,tbl,c)
    integer,intent(in):: n
    real(8),intent(in):: x,tbl(n)
    character(len=*),intent(in):: c

    if( x.lt.tbl(1) ) then
      write(6,'(2es12.4)') x,tbl(1)
      write(6,'(a)') ' ['//trim(c)//'] x.lt.tbl(1)'
      stop
    endif

  end subroutine check_range
!=======================================================================
end module Mishin_Al
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
