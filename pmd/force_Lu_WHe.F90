module Lu_WHe
  use pmdmpi
  use mod_precision
  use vector,only: dot

contains
  subroutine force_Lu_WHe(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,lstrs,iprint)
!-----------------------------------------------------------------------
! Parallel implementation of the G.-H. Lu potential and forces
!                                                     2013.07.08 by R.K.
!   - See, X.-C. Li et al., J. Nuclear Mater. 426 (2012) 31--37
!-----------------------------------------------------------------------
    implicit none
    include "./params_unit.h"
    include "./params_Lu_WHe.h"
    integer,intent(in):: namax,natm,nnmax,nismax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),mpi_md_world,myid_md,lspr(0:nnmax,namax),nex(3)
    real(rp),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),sv(3,6) &
         ,tag(namax),rc
    real(rp),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical:: lstrs

!-----locals
    integer:: i,j,jj,k,kk,ierr,is,js,ks
    real(rp):: rij,riji,rik,riki,rjk,rjki,rexp,aexp,dvrdr,cs,gc,tk,t1, &
         va,dgc,frik,dfrik,frjk,dfrjk,bij,bji,b,dvadr,frij,dfrij, &
         gfi,gfj,tmp,x,y,z,epotl,epott,vx
    real(rp):: R1,D1,R1ij,D1ij,R1ik,D1ik,R1jk,D1jk,beta,s,D0,r0 &
         ,gij,cij,dij,hij,x2i,x6i
    real(rp):: xi(3),xj(3),xk(3),xij(3),xji(3),xik(3),xjk(3),dixij(3) &
         ,djxij(3),dixji(3),djxji(3),fi(3),fj(3),dixik(3),djxjk(3) &
         ,dkxik(3),dkxjk(3),dics(3),djcs(3),dkcs(3),at(3)
!-----1st call
    logical,save:: l1st=.true.
    real(rp),save:: vc_HeHe,dvc_HeHe,xc
    real(rp),allocatable,dimension(:,:),save:: aa1,aa2

!-----only at 1st call
    if( l1st ) then
      allocate(aa1(3,namax),aa2(3,namax))
!.....cut-off for He-He term
      xc= p_HeHe_rc /p_HeHe_rm
      x2i= 1.0_rp/xc/xc
      x6i= x2i*x2i*x2i
      rexp= exp(-p_HeHe_alpha*xc +p_HeHe_beta*xc*xc)
      tmp=  p_HeHe_c6*x6i &
           +p_HeHe_c8*x6i*x2i &
           +p_HeHe_c10*x6i*x6i/x2i
      vc_HeHe= p_HeHe_A*rexp -f_x(xc,p_HeHe_D)*tmp
      dvc_HeHe= p_HeHe_A*rexp*(-p_HeHe_alpha +2.0_rp*p_HeHe_beta*xc) &
           -df_x(xc,p_HeHe_D)*tmp &
           -f_x(xc,p_HeHe_D)*( -6.0_rp*p_HeHe_c6*x6i &
           -8.0_rp*p_HeHe_c8*x6i*x2i &
           -10.0_rp*p_HeHe_c10*x6i*x6i/x2i )/xc
!.....If you do not want to apply cutoff, set vc_HeHe and dvc_HeHe as 0
!        vc_HeHe=0d0
!        dvc_HeHe= 0d0
!-------finally set l1st
      l1st=.false.
    endif

    if( size(aa1).lt.3*namax ) then
      deallocate(aa1,aa2)
      allocate(aa1(3,namax),aa2(3,namax))
    endif

!-----initialize
    aa1(1:3,1:namax)= 0.0_rp
    aa2(1:3,1:namax)= 0.0_rp
    epotl= 0.0_rp

!-----Repulsive term: V_R
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is= int(tag(i))
      do jj=1,lspr(0,i)
        j= lspr(jj,i)
        if(j.le.i) cycle
        js= int(tag(j))
        x= ra(1,j) -xi(1)
        y= ra(2,j) -xi(2)
        z= ra(3,j) -xi(3)
        xij(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
        rij= sqrt(dot(xij,xij))
!.....He-He
        if( is.eq.2 .and. js.eq.2 ) then
          if( rij.gt.p_HeHe_rc ) cycle
          x= rij/p_HeHe_rm
          x2i= 1.0_rp/x/x
          x6i= x2i*x2i*x2i
!.....potential
          rexp= exp(-p_HeHe_alpha*x +p_HeHe_beta*x*x)
          tmp=  p_HeHe_c6*x6i &
               +p_HeHe_c8*x6i*x2i &
               +p_HeHe_c10*x6i*x6i/x2i
          vx= p_HeHe_A *rexp -f_x(x,p_HeHe_D)*tmp &
               -vc_HeHe &
               -(x-xc)*dvc_HeHe
          vx= vx *p_HeHe_eps *0.5_rp
          epi(i)= epi(i) +vx
          if(j.le.natm) then
            epi(j)=epi(j)+vx
            epotl=epotl +vx +vx
          else
            epotl=epotl +vx
          endif
!.....force
          riji= 1.0_rp/rij
          dixij(1:3)= -xij(1:3)*riji
          djxij(1:3)=  xij(1:3)*riji
          dvrdr= p_HeHe_A*rexp*(-p_HeHe_alpha +2.0_rp*p_HeHe_beta*x) &
               -df_x(x,p_HeHe_D)*tmp &
               -f_x(x,p_HeHe_D)*( -6.0_rp*p_HeHe_c6*x6i &
               -8.0_rp*p_HeHe_c8*x6i*x2i &
               -10.0_rp*p_HeHe_c10*x6i*x6i/x2i )/x &
               -dvc_HeHe
          dvrdr= dvrdr *p_HeHe_eps /p_HeHe_rm
          aa1(1:3,i)= aa1(1:3,i) +dvrdr*dixij(1:3)
          aa1(1:3,j)= aa1(1:3,j) +dvrdr*djxij(1:3)

!.....other pairs
        else
          R1= p_R1(is,js)
          D1= p_D1(is,js)
!---------cutoff judgement
          if(rij.gt.R1+D1) cycle
          beta= p_beta(is,js)
          s= p_s(is,js)
          D0= p_D0(is,js)
          r0= p_r0(is,js)
!---------potential
          rexp= exp(-beta*sqrt(2.0_rp*s)*(rij-r0))
          tmp=0.5_rp*fc_r(rij,R1,D1)*D0/(s-1.0_rp)*rexp
          epi(i)=epi(i)+tmp
          if(j.le.natm) then
            epi(j)=epi(j)+tmp
            epotl=epotl +tmp +tmp
          else
            epotl=epotl +tmp
          endif
!---------force
          riji= 1.0_rp/rij
          dixij(1:3)= -xij(1:3)*riji
          djxij(1:3)=  xij(1:3)*riji
          dvrdr= D0/(s-1.0_rp)*rexp*( dfc_r(rij,R1,D1) &
               -beta*sqrt(2.0_rp*s)*fc_r(rij,R1,D1) )
          aa1(1:3,i)= aa1(1:3,i) +dvrdr*dixij(1:3)
          aa1(1:3,j)= aa1(1:3,j) +dvrdr*djxij(1:3)
        endif
      enddo
    enddo

!-----Attractive term: -Bij*V_A
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is=int(tag(i))
      do jj=1,lspr(0,i)
        j= lspr(jj,i)
        if(j.le.i) cycle
        js= int(tag(j))
!.....Skip He-He pair
        if( is.eq.2 .and. js.eq.2 ) cycle
        xj(1:3)= ra(1:3,j)
        x= xj(1) -xi(1)
        y= xj(2) -xi(2)
        z= xj(3) -xi(3)
        xij(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
        xji(1:3)= -xij(1:3)
        rij= sqrt(dot(xij,xij))
        R1ij= p_R1(is,js)
        D1ij= p_D1(is,js)
!---------cutoff judgement
        if(rij.gt.R1ij+D1ij) cycle
        beta= p_beta(is,js)
        s= p_s(is,js)
        D0= p_D0(is,js)
        r0= p_r0(is,js)
        riji= 1.0_rp/rij
        dixij(1:3)= -xij(1:3)*riji
        djxij(1:3)=  xij(1:3)*riji
        djxji(1:3)= -xji(1:3)*riji
        dixji(1:3)=  xji(1:3)*riji
        aexp= exp(-beta*sqrt(2.0_rp/s)*(rij-r0))
        va= fc_r(rij,R1ij,D1ij)*D0*s/(s-1.0_rp)*aexp
        dvadr= D0*s/(s-1.0_rp)*aexp*( dfc_r(rij,R1ij,D1ij) &
             -beta*sqrt(2.0_rp/s)*fc_r(rij,R1ij,D1ij) ) 
!.....parameters for g(theta)
        gij= p_gamma(is,js)
        cij= p_c(is,js)
        dij= p_d(is,js)
        hij= p_h(is,js)
!---------make gfi
        tk= 0.0_rp
        do kk=1,lspr(0,i)
          k= lspr(kk,i)
          if(k.eq.j) cycle
          ks= int(tag(k))
          x= ra(1,k) -xi(1)
          y= ra(2,k) -xi(2)
          z= ra(3,k) -xi(3)
          xik(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
          rik=sqrt(dot(xik,xik))
          R1ik= p_R1(is,ks)
          D1ik= p_D1(is,ks)
!-----------cutoff judgement
          if(rik.gt.R1ik+D1ik) cycle
          riki= 1.0_rp/rik
          cs= dot(xij,xik)*riji*riki
          gc= gij*(1.0_rp +cij**2/dij**2 &
               -(cij**2/(dij**2 +(hij+cs)**2)))
!.....alpha is always 0d0 in cases of W-W and W-He
          aexp= 1.0_rp
          tk=tk +gc*fc_r(rik,R1ik,D1ik) *aexp
        enddo
        gfi= 1.0_rp +tk
!---------scan k around i to get bij
        tk= 0.0_rp
        fi(1:3)= 0.0_rp
        fj(1:3)= 0.0_rp
        t1= -0.5_rp*gfi**(-1.5_rp)
        do kk=1,lspr(0,i)
          k= lspr(kk,i)
          if(k.eq.j) cycle
          ks=int(tag(k))
          x= ra(1,k) -xi(1)
          y= ra(2,k) -xi(2)
          z= ra(3,k) -xi(3)
          xik(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
          rik=sqrt(dot(xik,xik))
          R1ik= p_R1(is,ks)
          D1ik= p_D1(is,ks)
!-----------cutoff judgement
          if(rik.gt.R1ik+D1ik) cycle
          riki= 1.0_rp/rik
          dixik(1:3)= -xik(1:3)*riki
          dkxik(1:3)=  xik(1:3)*riki
          cs= dot(xij,xik)*riji*riki
          dics(1:3)= -riji*riki*(xij(1:3)+xik(1:3)) &
               -cs*(dixij(1:3)*riji+dixik(1:3)*riki)
          djcs(1:3)= riji*riki*xik(1:3) -cs*djxij(1:3)*riji
          dkcs(1:3)= riji*riki*xij(1:3) -cs*dkxik(1:3)*riki
          gc= gij*(1.0_rp +cij**2/dij**2 -(cij**2/(dij**2 &
               +(hij+cs)**2)))
          dgc= 2.0_rp*(hij+cs)*gij*cij**2/(dij**2+(hij+cs)**2)**2
          frik = fc_r(rik,R1ik,D1ik)
          dfrik= dfc_r(rik,R1ik,D1ik)
          tk=tk +gc*frik
          fi(1:3)=fi(1:3) +dics(1:3)*dgc*frik &
               +gc*dfrik*dixik(1:3)
          fj(1:3)=fj(1:3) +djcs(1:3)*dgc*frik
          aa2(1:3,k)=aa2(1:3,k) -0.5_rp*va*t1*( dkcs(1:3)*dgc*frik &
               +gc*dfrik*dkxik(1:3) )
        enddo
!---------derivative of bij part
        aa2(1:3,i)=aa2(1:3,i) -0.5_rp*va*t1*fi(1:3)
        aa2(1:3,j)=aa2(1:3,j) -0.5_rp*va*t1*fj(1:3)
!---------potential
        bij= gfi**(-0.5_rp)
        epi(i)=epi(i) -0.5_rp*0.5_rp*bij*va
        epi(j)=epi(j) -0.5_rp*0.5_rp*bij*va
        epotl=epotl -0.5_rp*bij*va
!---------derivative of va part
        aa2(1:3,i)=aa2(1:3,i) -0.5_rp*bij*dvadr*dixij(1:3)
        aa2(1:3,j)=aa2(1:3,j) -0.5_rp*bij*dvadr*djxij(1:3)
!---------if j.gt.natm, no need to calc the term around j
        if(j.gt.natm) cycle
!---------make gfj
        tk= 0.0_rp
        do kk=1,lspr(0,j)
          k= lspr(kk,j)
          if(k.eq.i) cycle
          ks=int(tag(k))
          x= ra(1,k) -xj(1)
          y= ra(2,k) -xj(2)
          z= ra(3,k) -xj(3)
          xjk(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
          rjk=sqrt(dot(xjk,xjk))
          R1jk= p_R1(js,ks)
          D1jk= p_D1(js,ks)
!-----------cutoff judgement
          if(rjk.gt.R1jk+D1jk) cycle
          rjki= 1.0_rp/rjk
          cs= dot(xji,xjk)*riji*rjki
          gc= gij*(1.0_rp +cij**2/dij**2 &
               -(cij**2/(dij**2 +(hij+cs)**2)))
!.....alpha is always 0d0 in cases of W-W and W-He
          aexp= 1.0_rp
          tk=tk +gc*fc_r(rjk,R1jk,D1jk) *aexp
        enddo
        gfj= 1.0_rp +tk
!---------scan k around j to get bji
        tk= 0.0_rp
        fi(1:3)= 0.0_rp
        fj(1:3)= 0.0_rp
        t1= -0.5_rp*gfj**(-1.5_rp)
        do kk=1,lspr(0,j)
          k= lspr(kk,j)
          if(k.eq.i) cycle
          ks=int(tag(k))
          x= ra(1,k) -xj(1)
          y= ra(2,k) -xj(2)
          z= ra(3,k) -xj(3)
          xjk(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
          rjk=sqrt(dot(xjk,xjk))
          R1jk= p_R1(js,ks)
          D1jk= p_D1(js,ks)
!-----------cutoff judgement
          if(rjk.gt.R1jk+D1jk) cycle
          rjki= 1.0_rp/rjk
          djxjk(1:3)= -xjk(1:3)*rjki
          dkxjk(1:3)=  xjk(1:3)*rjki
          cs= dot(xji,xjk)*riji*rjki
          djcs(1:3)= -riji*rjki*(xji(1:3)+xjk(1:3)) &
               -cs*(djxji(1:3)*riji+djxjk(1:3)*rjki)
          dics(1:3)= riji*rjki*xjk(1:3) -cs*dixji(1:3)*riji
          dkcs(1:3)= riji*rjki*xji(1:3) -cs*dkxjk(1:3)*rjki
          gc= gij*(1.0_rp +cij**2/dij**2 -(cij**2/(dij**2 &
               +(hij+cs)**2)))
          dgc= 2.0_rp*(hij+cs)*gij*cij**2/(dij**2+(hij+cs)**2)**2
          frjk = fc_r(rjk,R1jk,D1jk)
          dfrjk= dfc_r(rjk,R1jk,D1jk)
          tk=tk +gc*frjk
          fj(1:3)=fj(1:3) +djcs(1:3)*dgc*frjk &
               +gc*dfrjk*djxjk(1:3)
          fi(1:3)=fi(1:3) +dics(1:3)*dgc*frjk
          aa2(1:3,k)=aa2(1:3,k) -0.5_rp*va*t1*( dkcs(1:3)*dgc*frjk &
               +gc*dfrjk*dkxjk(1:3))
        enddo
!---------derivative of bji part
        aa2(1:3,i)=aa2(1:3,i) -0.5_rp*va*t1*fi(1:3)
        aa2(1:3,j)=aa2(1:3,j) -0.5_rp*va*t1*fj(1:3)
!---------potential
        bji= gfj**(-0.5_rp)
        epi(i)=epi(i) -0.5_rp*0.5_rp*bji*va
        epi(j)=epi(j) -0.5_rp*0.5_rp*bji*va
        epotl=epotl -0.5_rp*bji*va
!---------derivative of va part
        aa2(1:3,i)=aa2(1:3,i) -0.5_rp*bji*dvadr*dixij(1:3)
        aa2(1:3,j)=aa2(1:3,j) -0.5_rp*bji*dvadr*djxij(1:3)
      enddo
    enddo

    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_md_world,aa2,3)
    call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
         ,nn,mpi_md_world,epi,1)
!!$    if( myid_md.ge.0 ) then
!!$!-----send back forces and potentials on immigrants
!!$      call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
!!$           ,nn,mpi_md_world,aa2,3)
!!$      call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
!!$           ,nn,mpi_md_world,epi,1)
!!$    else
!!$      call reduce_dba_bk(natm,namax,tag,aa2,3)
!!$      call reduce_dba_bk(natm,namax,tag,epi,1)
!!$    endif

!-----sum
    aa(1:3,1:natm)= aa(1:3,1:natm) -aa1(1:3,1:natm) -aa2(1:3,1:natm)

!-----gather epot
    call mpi_allreduce(epotl,epott,1,mpi_real_rp &
         ,MPI_SUM,mpi_md_world,ierr)
    epot= epot +epott
    
  end subroutine force_Lu_WHe
!=======================================================================
  function fc_r(r,r1,d1)
    implicit none
    real(rp),intent(in):: r,r1,d1
    real(rp):: fc_r
    real(rp),parameter:: pi = 3.14159265358979_rp

    fc_r= 0.0_rp
    if(r.lt.r1-d1) then
      fc_r= 1.0_rp
    elseif(r1-d1.le.r .and. r.lt.r1+d1) then
      fc_r= 0.5_rp*( 1.0_rp -sin(0.5_rp*pi*(r-r1)/d1) )
    endif
    return
  end function fc_r
!=======================================================================
  function dfc_r(r,r1,d1)
!-----Derivative: df(r)/dr
    implicit none 
    real(rp),intent(in):: r,r1,d1
    real(rp):: dfc_r
    real(rp),parameter:: pi = 3.14159265358979_rp

    dfc_r= 0.0_rp
    if(r1-d1.le.r .and. r.lt.r1+d1) then
      dfc_r= -0.5_rp*cos(0.5_rp*pi*(r-r1)/d1)*0.5_rp*pi/d1
    endif
    return
  end function dfc_r
!=======================================================================
  function f_x(x,d)
    implicit none
    real(rp),intent(in):: x,d
    real(rp):: f_x,tmp

    if( x.lt.d ) then
      tmp= d/x -1.0_rp
      f_x= exp(-tmp*tmp)
    else
      f_x= 1.0_rp
    endif

    return
  end function f_x
!=======================================================================
  function df_x(x,d)
    implicit none
    real(rp),intent(in):: x,d
    real(rp):: df_x,tmp

    if( x.lt.d ) then
      tmp= d/x -1.0_rp
      df_x= 2.0_rp*d/x/x *tmp *exp(-tmp*tmp)
    else
      df_x= 0.0_rp
    endif
    return
  end function df_x
end module Lu_WHe
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
