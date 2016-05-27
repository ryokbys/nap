module EDIP_Si
contains
  subroutine force_EDIP_Si(namax,natm,tag,ra,nnmax,aa,strs,h,hi,tcom &
       ,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_world,myid,epi,epot,nismax,acon,lstrs)
!-----------------------------------------------------------------------
!  Parallel implementation of EDIP(Si) force calculation for pmd
!    - Environment Dependent Interatomic Potential (EDIP) for Si
!      Ref: PRB 58, 2539 (1998), J.F.Just et al.
!    - 2010.03.31 by R.K.
!-----------------------------------------------------------------------
    implicit none
    include "mpif.h"
    include "./params_unit.h"
    include "params_EDIP_Si.h"
    integer,intent(in):: namax,natm,nnmax,nismax
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_world,myid,lspr(0:nnmax,namax)
    real(8),intent(in):: ra(3,namax),tag(namax),acon(nismax) &
         ,h(3,3),hi(3,3),sv(3,6),rc
    real(8),intent(inout):: tcom
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical:: lstrs

!-----local
    integer:: i,j,k,l,m,n,ixyz,is,js,ks,ierr,jj,kk,ll
    real(8):: epotl,epotl2,epotl3
    real(8):: rij,rij2,riji,rim,rimi,v2,t1,t2,t3,dft,aexp,gg &
         ,eda2,eda,edb,edc,edg,eds
    real(8):: v3,qz,dqz,tz,dtz,grij,dgrij,grik,dgrik,rik,rik2,riki,h1 &
         ,h2,cs,hijk
    real(8),save,allocatable,dimension(:):: z,pz,xi,xj,xij,dxi,dxj,at &
         ,xx,dixij,djxij,dixik,dkxik,xk &
         ,xik,dcsi,dcsj,dcsk,dit1,dlt1
    real(8),save,allocatable:: dz(:,:,:)
    real(8),save,allocatable:: aa2(:,:),aa3(:,:)
    real(8),save,allocatable:: teda(:,:),tedb(:,:),tedc(:,:) &
         ,tedg(:,:),teds(:,:)
!-----1st call
    logical,save:: l1st=.true.

!-----only at 1st call
    if( l1st ) then
      allocate(aa2(3,namax),aa3(3,namax))
      allocate(z(namax),pz(namax),xi(3),xj(3),xij(3),dxi(3),dxj(3) &
           ,at(3),xx(3))
      allocate(dz(3,0:nnmax,namax))
      allocate(dixij(3),djxij(3),dixik(3),dkxik(3),xk(3),xik(3) &
           ,dcsi(3),dcsj(3),dcsk(3),dit1(3),dlt1(3))
      allocate(teda(2,2),tedb(2,2),tedc(2,2),teds(2,2),tedg(2,2))
      teda(1,1)= ed_a
      teda(2,2)= ed_a*ratio
      teda(1,2)= (teda(1,1)+teda(2,2))*0.5d0
      teda(2,1)= (teda(1,1)+teda(2,2))*0.5d0
      tedb(1,1)= ed_bb
      tedb(2,2)= ed_bb*ratio
      tedb(1,2)= (tedb(1,1)+tedb(2,2))*0.5d0
      tedb(2,1)= (tedb(1,1)+tedb(2,2))*0.5d0
      tedc(1,1)= ed_c
      tedc(2,2)= ed_c*ratio
      tedc(1,2)= (tedc(1,1)+tedc(2,2))*0.5d0
      tedc(2,1)= (tedc(1,1)+tedc(2,2))*0.5d0
      tedg(1,1)= ed_gam
      tedg(2,2)= ed_gam*ratio
      tedg(1,2)= (tedg(1,1)+tedg(2,2))*0.5d0
      tedg(2,1)= (tedg(1,1)+tedg(2,2))*0.5d0
      teds(1,1)= ed_sgm
      teds(2,2)= ed_sgm*ratio
      teds(1,2)= (teds(1,1)+teds(2,2))*0.5d0
      teds(2,1)= (teds(1,1)+teds(2,2))*0.5d0
!-------check rc
      if( int(rc*100d0).ne.int(max(teda(1,1),teda(2,2))*100d0) ) then
        if(myid.eq.0) then
          write(6,'(1x,a)') "!!! Cutoff radius is not appropriate !!!"
          write(6,'(1x,a,es12.4)') "rc should be" &
               ,max(teda(1,1),teda(2,2))
        endif
        call mpi_finalize(ierr)
        stop
      endif
!-------finally set l1st
      l1st=.false.
    endif

    epotl= 0d0
    epi(1:namax)= 0d0
    aa2(1:3,1:namax)= 0d0
    aa3(1:3,1:namax)= 0d0
    epotl2= 0d0
    epotl3= 0d0

!-----set Z
    z(1:natm)= 0d0
    dz(1:3,0:nnmax,natm)= 0d0
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is= int(tag(i))
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
        js= int(tag(j))
        xx(1:3)= ra(1:3,j)-xi(1:3)
        xij(1:3)= h(1:3,1)*xx(1) +h(1:3,2)*xx(2) +h(1:3,3)*xx(3)
        rij2=xij(1)*xij(1) +xij(2)*xij(2) +xij(3)*xij(3)
        eda= teda(is,js)
        eda2= eda*eda
        if( rij2.ge.eda2 ) cycle
        rij= dsqrt(rij2)
        edc= tedc(is,js)
        z(i)= z(i) +f_r(eda,edc,rij,ed_alp)
        dz(1:3,k,i)= df_r(eda,edc,rij,ed_alp)*xij(1:3)/rij
        dz(1:3,0,i)= dz(1:3,0,i) -dz(1:3,k,i)
      enddo
      pz(i)= dexp(-ed_bet*z(i)*z(i))
    enddo

!-----2-body term
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is= int(tag(i))
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
        js= int(tag(j))
        xx(1:3)= ra(1:3,j)-xi(1:3)
        xij(1:3)= h(1:3,1)*xx(1) +h(1:3,2)*xx(2) +h(1:3,3)*xx(3)
        rij2=xij(1)*xij(1) +xij(2)*xij(2) +xij(3)*xij(3)
        eda= teda(is,js)
        eda2= eda*eda
        if( rij2.ge.eda2 ) cycle
        rij= dsqrt(rij2)
        riji= 1d0/rij
        eds= teds(is,js)
        aexp= exp(eds/(rij-eda))
!---------potential
        edb= tedb(is,js)
        v2= ed_aa*((edb*riji)**ed_rho-pz(i))*aexp
        epi(i)= epi(i) +v2
        epotl2= epotl2 +v2
!---------force
        t1= -ed_rho*(edb*riji)**ed_rho*riji *ed_aa*aexp
        t2= -2d0*ed_bet*z(i)*pz(i) *ed_aa*aexp
        t3= -eds/(rij-eda)**2 *v2
        dxi(1:3)= -xij(1:3)*riji
        dxj(1:3)=  xij(1:3)*riji
        aa2(1:3,i)= aa2(1:3,i) +dxi(1:3)*(t1+t3) -t2*dz(1:3,0,i)
        aa2(1:3,j)= aa2(1:3,j) +dxj(1:3)*(t1+t3)
!---------This code works because dz(1:3,l,i)=0d0
!--------- if l-th neighbor is outside the cutoff
        do l=1,lspr(0,i)
          m=lspr(l,i)
          aa2(1:3,m)= aa2(1:3,m) -t2*dz(1:3,l,i)
        enddo
      enddo
    enddo

!-----3-body term
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is= int(tag(i))
      qz=ed_q0*dexp(-ed_mu*z(i))
      dqz= -ed_mu*qz
      tz= ed_u1 +ed_u2*(ed_u3*dexp(-ed_u4*z(i))-dexp(-2d0*ed_u4*z(i)))
      dtz= -ed_u2*ed_u4*(ed_u3*dexp(-ed_u4*z(i)) &
           -2d0*dexp(-2d0*ed_u4*z(i)))
      do jj=1,lspr(0,i)
        j=lspr(jj,i)
        if(j.eq.0) exit
        js= int(tag(j))
        xx(1:3)= ra(1:3,j)-xi(1:3)
        xij(1:3)= ( h(1:3,1)*xx(1) +h(1:3,2)*xx(2) +h(1:3,3)*xx(3) )
        rij2=xij(1)*xij(1) +xij(2)*xij(2) +xij(3)*xij(3)
        eda= teda(is,js)
        eda2= eda*eda
        if( rij2.ge.eda2 ) cycle
        rij= dsqrt(rij2)
        edg= tedg(is,js)
        grij= dexp(edg/(rij-eda))
        riji= 1d0/rij
        dgrij= -edg/(rij-eda)**2*grij
        dixij(1:3)= -xij(1:3)*riji
        djxij(1:3)=  xij(1:3)*riji
        do kk=1,lspr(0,i)
          k=lspr(kk,i)
          if(k.eq.0) exit
          if( k.le.j ) cycle
          ks= int(tag(k))
          xx(1:3)= ra(1:3,k)-xi(1:3)
          xik(1:3)= ( h(1:3,1)*xx(1) +h(1:3,2)*xx(2) +h(1:3,3)*xx(3) )
          rik2= xik(1)*xik(1)+xik(2)*xik(2)+xik(3)*xik(3)
          eda= teda(is,ks)
          eda2= eda*eda
          if( rik2.ge.eda2 ) cycle
          rik=dsqrt(rik2)
          edg= tedg(is,ks)
          grik= dexp(edg/(rik-eda))
          dgrik= -edg/(rik-eda)**2 *grik
          riki= 1d0/rik
          dixik(1:3)= -xik(1:3)*riki
          dkxik(1:3)=  xik(1:3)*riki
          cs=(xij(1)*xik(1)+xij(2)*xik(2)+xij(3)*xik(3))*riji*riki
          t1= qz*(cs+tz)**2
          aexp= dexp(-t1)
          h1= ed_lam*(1d0-aexp)
          h2= ed_lam*ed_eta*t1
          hijk= h1 +h2
!-----------potential
          gg= grij*grik
          v3= gg*hijk
          epotl3=epotl3 +v3
!-----------force calc.
          t2=(xij(1)*xik(1)+xij(2)*xik(2)+xij(3)*xik(3)) &
               *riji*riji *riki*riki
          dcsi(1:3)= -(xij(1:3)+xik(1:3))*riji*riki &
               -t2*(rik*dixij(1:3)+rij*dixik(1:3))
          dcsj(1:3)= xik(1:3)*riji*riki -t2*rik*djxij(1:3)
          dcsk(1:3)= xij(1:3)*riji*riki -t2*rij*dkxik(1:3)
!-----------deriv. of grij, grik
          aa3(1:3,i)= aa3(1:3,i) +dixij(1:3)*dgrij*grik*hijk &
               +dixik(1:3)*dgrik*grij*hijk
          aa3(1:3,j)= aa3(1:3,j) +djxij(1:3)*dgrij*grik*hijk
          aa3(1:3,k)= aa3(1:3,k) +dkxik(1:3)*dgrik*grij*hijk
!-----------deriv. of t1 about i
          dit1(1:3)=dz(1:3,0,i)*dqz*(cs+tz)**2 &
               +2d0*qz*(cs+tz)*(dcsi(1:3)+dz(1:3,0,i)*dtz)
!-----------deriv. of h2 about i
          aa3(1:3,i)=aa3(1:3,i) +dit1(1:3)*ed_lam*ed_eta*gg
!-----------deriv. of h1 about i
          aa3(1:3,i)=aa3(1:3,i) +dit1(1:3)*ed_lam*aexp*gg
!-----------deriv. of h2 about j,k only cs part
          aa3(1:3,j)=aa3(1:3,j) &
               +2d0*qz*(cs+tz)*dcsj(1:3)*ed_lam*ed_eta*gg
          aa3(1:3,k)=aa3(1:3,k) &
               +2d0*qz*(cs+tz)*dcsk(1:3)*ed_lam*ed_eta*gg
!-----------deriv. of h1 about j,k only cs part
          aa3(1:3,j)=aa3(1:3,j) &
               +2d0*qz*(cs+tz)*dcsj(1:3)*ed_lam*aexp*gg
          aa3(1:3,k)=aa3(1:3,k) &
               +2d0*qz*(cs+tz)*dcsk(1:3)*ed_lam*aexp*gg
          do ll=1,lspr(0,i)
            l=lspr(ll,i)
!-------------deriv. of t1 about l except cs part
!              dlt1(1:3)= dz(1:3,ll,i)*dqz*(cs+tz)**2
!     &             +2d0*qz*(cs+tz)*dz(1:3,ll,i)*dtz
            dlt1(1:3)= dz(1:3,ll,i)*(cs+tz) &
                 *( dqz*(cs+tz) +2d0*qz*dtz )
!-------------deriv. of h2 about l except cs part
            aa3(1:3,l)=aa3(1:3,l)+dlt1(1:3)*ed_lam*ed_eta*gg
!-------------deriv. of h1 about l except cs part
            aa3(1:3,l)=aa3(1:3,l)+dlt1(1:3)*ed_lam*aexp*gg
          enddo
        enddo
      enddo
    enddo

    if( myid.ge.0 ) then
!-----send back (3-body) forces and potentials on immigrants
      call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,lsrc,myparity &
           ,nn,mpi_world,aa2,3)
      call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,lsrc,myparity &
           ,nn,mpi_world,aa3,3)
      call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,lsrc,myparity &
           ,nn,mpi_world,epi,1)
    else
      call reduce_dba_bk(natm,namax,tag,aa3,3)
      call reduce_dba_bk(natm,namax,tag,aa2,3)
      call reduce_dba_bk(natm,namax,tag,epi,1)
    endif

!-----sum
    aa(1:3,1:natm)= -aa2(1:3,1:natm) -aa3(1:3,1:natm)

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
    epotl= epotl2 +epotl3
    if( myid.ge.0 ) then
      call mpi_allreduce(epotl,epot,1,MPI_DOUBLE_PRECISION &
           ,MPI_SUM,mpi_world,ierr)
    else
      epot= epotl
    endif

  end subroutine force_EDIP_Si
!=======================================================================
  function f_r(a,c,r,alpha)
!-----Weighting function of EDIP
    implicit none
    real(8),intent(in):: a,c,r,alpha
    real(8):: f_r,x

    f_r= 0d0
    if(r.lt.c) then
      f_r= 1d0
    elseif(c.le.r .and. r.lt.a) then
      x= (r-c)/(a-c)
      f_r= exp(alpha/(1d0-1d0/x**3))
    endif
    return
  end function f_r
!=======================================================================
  function df_r(a,c,r,alpha)
!-----Derivative: df(r)/dr
    implicit none 
    real(8),intent(in):: a,c,r,alpha
    real(8):: df_r,x,t1,t2,t3

    df_r= 0d0
    if(c.le.r .and. r.lt.a) then
      x= (r-c)/(a-c)
      t1= 1d0/(a-c)
      t2= exp(alpha/(1d0-1d0/x**3))
      t3= -3d0/x**4*alpha/(1d0 -1d0/x**3)**2
      df_r= t1 *t2 *t3
    endif
    return
  end function df_r
end module EDIP_Si
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
