module Branicio_AlN
contains
  subroutine force_Branicio_AlN(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_world,myid,epi,epot,nismax,lstrs,iprint)
!-----------------------------------------------------------------------
!  Parallel implementation of Branicio potential for AlN
!    - Branicio potential (Vashishta group)
!      Ref: Branicio et al., J. Mech. Phys. Solids 56 (2008) pp.1955.
!    - 2010.04.02 by R.K.
!-----------------------------------------------------------------------
    implicit none
    include "mpif.h"
    include "./params_unit.h"
    include "params_Branicio_AlN.h"
    integer,intent(in):: namax,natm,nnmax,nismax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_world,myid,lspr(0:nnmax,namax),nex(3)
    real(8),intent(in):: ra(3,namax),tag(namax) &
         ,h(3,3),hi(3,3),sv(3,6),rc
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical:: lstrs

    integer:: i,j,k,l,m,n,ierr,is,js,ks,ir
    real(8):: f2rc,df2rc,ri,cst,cs
    real(8):: epotl,epotl2,epotl3,t1,t2,tmp1,tmp2,vol,epott
    real(8):: rij,riji,rik,riki,d,v2,dv2,expij,expik,v3,dt1j,dt1k,dt2
    real(8):: xi(3),xj(3),xij(3),xik(3),at(3),xx(3),drij(3),drik(3) &
         ,dcsj(3),dcsk(3),dcsi(3)

!-----saved values
    real(8),save:: rmin,rmax,dr
!-----saved allocatable arrays
    real(8),allocatable,save:: tblf2(:,:,:),tbldf2(:,:,:)
    real(8),allocatable,save:: aa2(:,:),aa3(:,:)
!-----1st call
    logical,save:: l1st=.true.

!-----only at 1st call
    if( l1st ) then
!-------allocate 2-body force table at the 1st call (do not deallocate!)
      allocate(tblf2(nd2b,2,2),tbldf2(nd2b,2,2))
!-------make 2-body (smoothed) force table
      rmin= 0.5d0
      rmax= rc
      dr= (rmax-rmin)/(nd2b-1)
      do is=1,2
        do js=1,2
          f2rc = f2_r(rc,is,js)
          df2rc=df2_r(rc,is,js)
          do i=1,nd2b
            ri= rmin +dr*(i-1)
            tblf2(i,is,js) = f2_r(ri,is,js) -f2rc -(ri-rc)*df2rc
            tbldf2(i,is,js)=df2_r(ri,is,js) -df2rc
          enddo
        enddo
      enddo
      do i=1,nd2b
        ri= rmin +dr*(i-1)
        write(90,'(10es12.4)') ri,tblf2(i,1:2,1:2) &
             ,tbldf2(i,1:2,1:2)
      enddo
!        write(6,'(a,2es12.4)') " cfct=",cfct
      write(6,'(a,2es12.4)') " alp =",v_alp(1:2)
      write(6,'(a,4es12.4)') " n   =",v_n(1:2,1:2)
      write(6,'(a,4es12.4)') " h   =",v_h(1:2,1:2)
      write(6,'(a,4es12.4)') " w   =",v_w(1:2,1:2)
      write(6,'(a,4es12.4)') " r1s =",v_r1s(1:2,1:2)
      write(6,'(a,4es12.4)') " r4s =",v_r4s(1:2,1:2)

      allocate(aa2(3,namax),aa3(3,namax))

!-------finally set l1st
      l1st=.false.
    endif

    epotl= 0d0
    aa2(1:3,1:namax)= 0d0
    aa3(1:3,1:namax)= 0d0
    epotl2= 0d0
    epotl3= 0d0

!-----2-body term
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is= int(tag(i))
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if( j.le.i ) cycle
        js= int(tag(j))
        xx(1:3)= ra(1:3,j) -xi(1:3)
        xij(1:3)= h(1:3,1)*xx(1) +h(1:3,2)*xx(2) +h(1:3,3)*xx(3)
        rij= dsqrt(xij(1)*xij(1) +xij(2)*xij(2) +xij(3)*xij(3))
        ir= int( (rij-rmin)/dr +1 )
        d = (rij-rmin)/dr -(ir-1)
!---------potential
        v2= tblf2(ir,is,js) +(tblf2(ir+1,is,js)-tblf2(ir,is,js))*d
        v2= v2 /2
        epi(i)= epi(i) +v2
        epotl2= epotl2 +v2
        if( j.le.natm ) then
          epi(j)= epi(j) +v2
          epotl2= epotl2 +v2
        endif
!---------force
        drij(1:3)= -xij(1:3)/rij
        dv2= tbldf2(ir,is,js) &
             +(tbldf2(ir+1,is,js)-tbldf2(ir,is,js))*d
        aa2(1:3,i)= aa2(1:3,i) +drij(1:3)*dv2
        aa2(1:3,j)= aa2(1:3,j) -drij(1:3)*dv2
      enddo
    enddo

!-----3-body term
    cst= cos(v_tht)
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is= int(tag(i))
      do m=1,lspr(0,i)
        j=lspr(m,i)
        js= int(tag(j))
        xx(1:3)= ra(1:3,j) -xi(1:3)
        xij(1:3)= h(1:3,1)*xx(1) +h(1:3,2)*xx(2) +h(1:3,3)*xx(3)
        rij= dsqrt(xij(1)*xij(1) +xij(2)*xij(2) +xij(3)*xij(3))
        riji= 1d0/rij
        drij(1:3)= -xij(1:3)*riji
        expij= exp(v_xi/(rij-v_r0))
        do n=1,lspr(0,i)
          k=lspr(n,i)
          if( k.le.j ) cycle
          xx(1:3)= ra(1:3,k) -xi(1:3)
          xik(1:3)= h(1:3,1)*xx(1) +h(1:3,2)*xx(2) +h(1:3,3)*xx(3)
          rik=dsqrt(xik(1)*xik(1)+xik(2)*xik(2)+xik(3)*xik(3))
          riki= 1d0/rik
          drik(1:3)= -xik(1:3)*riki
          expik= exp(v_xi/(rik-v_r0))
!-----------potential
          t1= v_b*expij*expik
          cs= (xij(1)*xik(1)+xij(2)*xik(2)+xij(3)*xik(3))*riji*riki
          t2= (cs-cst)**2 /(1d0+v_c*(cs-cst)**2)
          v3= t1*t2
          epi(i)= epi(i) +v3
          epotl3= epotl3 +v3
!-----------force
          dt1j= v_b*expij*expik*(-v_xi/(rij-v_r0)**2)
          dt1k= v_b*expij*expik*(-v_xi/(rik-v_r0)**2)
          dt2= 2d0*(cs-cst) /( 1d0 +v_c*(cs-cst)**2 )**2
          dcsj(1:3)= xik(1:3)*riji*riki -xij(1:3)*cs*riji**2
          dcsk(1:3)= xij(1:3)*riji*riki -xik(1:3)*cs*riki**2
          dcsi(1:3)= -dcsj(1:3) -dcsk(1:3)
          aa3(1:3,i)=aa3(1:3,i) +t2*( dt1j*drij(1:3)+dt1k*drik(1:3) ) &
               +t1*dt2*dcsi(1:3)
          aa3(1:3,j)=aa3(1:3,j) +t2*dt1j*(-drij(1:3)) &
               +t1*dt2*dcsj(1:3)
          aa3(1:3,k)=aa3(1:3,k) +t2*dt1k*(-drik(1:3)) &
               +t1*dt2*dcsk(1:3)
        enddo
      enddo
    enddo

!-----send back (3-body) forces and potentials on immigrants
    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,aa3,3)
    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_world,epi,1)
!!$    if( myid.ge.0 ) then
!!$!-----send back (3-body) forces and potentials on immigrants
!!$      call copy_dba_bk(namax,natm,nbmax,nb,lsb,lsrc,myparity &
!!$           ,nn,mpi_world,aa3,3)
!!$      call copy_dba_bk(namax,natm,nbmax,nb,lsb,lsrc,myparity &
!!$           ,nn,mpi_world,epi,1)
!!$    else
!!$      call reduce_dba_bk(natm,namax,tag,aa3,3)
!!$      call reduce_dba_bk(natm,namax,tag,epi,1)
!!$    endif

!-----sum
    aa(1:3,1:natm)= aa(1:3,1:natm) -aa2(1:3,1:natm) -aa3(1:3,1:natm)

!-----gather epot
    epotl= epotl2 +epotl3
    if( myid.ge.0 ) then
      epott = 0d0
      call mpi_allreduce(epotl,epott,1,MPI_DOUBLE_PRECISION &
           ,MPI_SUM,mpi_world,ierr)
      epot= epot +epott
    else
      epot= epot +epotl
    endif

!!$!-----get min bond length
!!$    tmp1= 1d10
!!$    do k=1,lspr(0,1)
!!$      j=lspr(k,1)
!!$      xi(1:3)= ra(1:3,j) -ra(1:3,1)
!!$      xij(1:3)= h(1:3,1)*xi(1) +h(1:3,2)*xi(2) +h(1:3,3)*xi(3)
!!$      rij= sqrt(xij(1)**2 +xij(2)**2 +xij(3)**2)
!!$      tmp1= min(tmp1,rij)
!!$    enddo
!!$!-----output lattice constant
!!$    vol= h(1,1)*h(2,2)*h(3,3) *0.529177d0**3
!!$    write(92,'(10es12.4)') tmp1*0.529177d0,vol/natm,epot/natm

  end subroutine force_Branicio_AlN
!=======================================================================
  function f2_r(r,is,js)
!-----2-body force
    implicit none
    include "./params_unit.h"
    include "params_Branicio_AlN.h"
    real(8),intent(in):: r
    integer,intent(in):: is,js

    real(8):: dij
!-----value
    real(8):: f2_r

    dij= 0.5d0 *( v_alp(is)*v_z(is)**2 +v_alp(js)*v_z(js)**2 )
    f2_r= &
         v_h(is,js)/r**v_n(is,js) &
         +v_z(is)*v_z(js)/r *exp(-r/v_r1s(is,js)) &
         -dij/r**4 *exp(-r/v_r4s(is,js)) &
         -v_w(is,js)/r**6

    if(is.eq.1 .and. js.eq.2) then
      write(91,'(2i6,5es12.4)') is,js,r &
           ,v_h(is,js)/r**v_n(is,js) &
           ,v_z(is)*v_z(js)/r *exp(-r/v_r1s(is,js)) &
           ,-dij/r**4 *exp(-r/v_r4s(is,js)) &
           ,-v_w(is,js)/r**6
    endif

    return
  end function f2_r
!=======================================================================
  function df2_r(r,is,js)
!-----Derivative of 2-body term
    implicit none 
    include "./params_unit.h"
    include "params_Branicio_AlN.h"
    real(8),intent(in):: r
    integer,intent(in):: is,js

    real(8):: df2_r,dij

    dij= 0.5d0 *( v_alp(is)*v_z(is)**2 +v_alp(js)*v_z(js)**2 )
    df2_r= &
         -v_n(is,js)*v_h(is,js)/r**(v_n(is,js)+1d0) &
         -v_z(is)*v_z(js)/r *exp(-r/v_r1s(is,js)) &
         *(1d0/r +1d0/v_r1s(is,js)) &
         +dij/r**4 *exp(-r/v_r4s(is,js)) &
         *(4d0/r +1d0/v_r4s(is,js)) &
         +6d0*v_w(is,js)/r**7d0

    return
  end function df2_r
end module Branicio_AlN
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
