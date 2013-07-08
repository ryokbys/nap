      subroutine get_force(namax,natm,tag,ra,nnmax,aa,strs,h,hi,tcom
     &     ,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
     &     ,mpi_md_world,myid_md,epi,epot,nismax,acon)
c-----------------------------------------------------------------------
c Parallel implementation of the G.-H. Lu's potential and forces
c                                                     2013.07.08 by R.K.
c   - See, X.-C. Li et al., J. Nuclear Mater. 426 (2012) 31--37
c-----------------------------------------------------------------------
      implicit none
      include "mpif.h"
      include "./params_au.h"
      include "./params_Lu_W-He.h"
      integer,intent(in):: namax,natm,nnmax,nismax
      integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3)
     &     ,nn(6),mpi_md_world,myid_md,lspr(0:nnmax,namax)
      real(8),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),sv(3,6)
     &     ,acon(nismax),tag(namax),rc
      real(8),intent(inout):: tcom
      real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)

c-----locals
      integer:: i,j,jj,k,kk,ierr,is,js,ks
      real(8):: rij,riji,rik,riki,rjk,rjki,rexp,aexp,dvrdr,cs,gc,tk,t1,
     &     va,dgc,frik,dfrik,frjk,dfrjk,bij,bji,b,dvadr,frij,dfrij,
     &     gfi,gfj,tmp,x,y,z,epotl,vx
      real(8):: R1,D1,R1ij,D1ij,R1ik,D1ik,R1jk,D1jk,beta,s,D0,r0
     &     ,gij,cij,dij,hij,x2i,x6i
      real(8):: xi(3),xj(3),xk(3),xij(3),xji(3),xik(3),xjk(3),dixij(3)
     &     ,djxij(3),dixji(3),djxji(3),fi(3),fj(3),dixik(3),djxjk(3)
     &     ,dkxik(3),dkxjk(3),dics(3),djcs(3),dkcs(3),at(3)
c-----functions
      real(8),external:: fc_r,dfc_r,sprod,f_x,df_x
c-----1st call
      logical,save:: l1st=.true.
      real(8),save:: vc_HeHe,dvc_HeHe,xc
      real(8),allocatable,dimension(:,:),save:: aa1,aa2

c-----only at 1st call
      if( l1st ) then
        allocate(aa1(3,namax),aa2(3,namax))
c.....cut-off for He-He term
        xc= p_HeHe_rc /p_HeHe_rm
        x2i= 1d0/xc/xc
        x6i= x2i*x2i*x2i
        rexp= exp(-p_HeHe_alpha*xc +p_HeHe_beta*xc*xc)
        tmp=  p_HeHe_c6*x6i
     &       +p_HeHe_c8*x6i*x2i
     &       +p_HeHe_c10*x6i*x6i/x2i
        vc_HeHe= p_HeHe_A*rexp -f_x(xc,p_HeHe_D)*tmp
        dvc_HeHe= p_HeHe_A*rexp*(-p_HeHe_alpha +2d0*p_HeHe_beta*xc)
     &       -df_x(xc,p_HeHe_D)*tmp
     &       -f_x(xc,p_HeHe_D)*( -6d0*p_HeHe_c6*x6i
     &       -8d0*p_HeHe_c8*x6i*x2i
     &       -10d0*p_HeHe_c10*x6i*x6i/x2i )/xc
c.....If you do not want to apply cutoff, set vc_HeHe and dvc_HeHe as 0
c        vc_HeHe=0d0
c        dvc_HeHe= 0d0
c-------finally set l1st
        l1st=.false.
      endif

c-----initialize
      aa1(1:3,1:namax)= 0d0
      aa2(1:3,1:namax)= 0d0
      epotl= 0d0
      epi(1:namax)= 0d0
      strs(1:3,1:3,1:namax)= 0d0

c-----Repulsive term: V_R
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
          rij= dsqrt(sprod(xij,xij))
c.....He-He
          if( is.eq.2 .and. js.eq.2 ) then
            if( rij.gt.p_HeHe_rc ) cycle
            x= rij/p_HeHe_rm
            x2i= 1d0/x/x
            x6i= x2i*x2i*x2i
c.....potential
            rexp= exp(-p_HeHe_alpha*x +p_HeHe_beta*x*x)
            tmp=  p_HeHe_c6*x6i
     &           +p_HeHe_c8*x6i*x2i
     &           +p_HeHe_c10*x6i*x6i/x2i
            vx= p_HeHe_A *rexp -f_x(x,p_HeHe_D)*tmp
     &           -vc_HeHe
     &           -(x-xc)*dvc_HeHe
            vx= vx *p_HeHe_eps *0.5d0
            epi(i)= epi(i) +vx
            if(j.le.natm) then
              epi(j)=epi(j)+vx
              epotl=epotl +vx +vx
            else
              epotl=epotl +vx
            endif
c.....force
            riji= 1d0/rij
            dixij(1:3)= -xij(1:3)*riji
            djxij(1:3)=  xij(1:3)*riji
            dvrdr= p_HeHe_A*rexp*(-p_HeHe_alpha +2d0*p_HeHe_beta*x)
     &           -df_x(x,p_HeHe_D)*tmp
     &           -f_x(x,p_HeHe_D)*( -6d0*p_HeHe_c6*x6i
     &           -8d0*p_HeHe_c8*x6i*x2i
     &           -10d0*p_HeHe_c10*x6i*x6i/x2i )/x
     &           -dvc_HeHe
            dvrdr= dvrdr *p_HeHe_eps /p_HeHe_rm
            aa1(1:3,i)= aa1(1:3,i) +dvrdr*dixij(1:3)
            aa1(1:3,j)= aa1(1:3,j) +dvrdr*djxij(1:3)

c.....other pairs
          else
            R1= p_R1(is,js)
            D1= p_D1(is,js)
c---------cutoff judgement
            if(rij.gt.R1+D1) cycle
            beta= p_beta(is,js)
            s= p_s(is,js)
            D0= p_D0(is,js)
            r0= p_r0(is,js)
c---------potential
            rexp= dexp(-beta*dsqrt(2d0*s)*(rij-r0))
            tmp=0.5d0*fc_r(rij,R1,D1)*D0/(s-1d0)*rexp
            epi(i)=epi(i)+tmp
            if(j.le.natm) then
              epi(j)=epi(j)+tmp
              epotl=epotl +tmp +tmp
            else
              epotl=epotl +tmp
            endif
c---------force
            riji= 1d0/rij
            dixij(1:3)= -xij(1:3)*riji
            djxij(1:3)=  xij(1:3)*riji
            dvrdr= D0/(s-1d0)*rexp*( dfc_r(rij,R1,D1)
     &           -beta*dsqrt(2d0*s)*fc_r(rij,R1,D1) )
            aa1(1:3,i)= aa1(1:3,i) +dvrdr*dixij(1:3)
            aa1(1:3,j)= aa1(1:3,j) +dvrdr*djxij(1:3)
          endif
        enddo
      enddo

c-----Attractive term: -Bij*V_A
      do i=1,natm
        xi(1:3)= ra(1:3,i)
        is=int(tag(i))
        do jj=1,lspr(0,i)
          j= lspr(jj,i)
          if(j.le.i) cycle
          js= int(tag(j))
c.....Skip He-He pair
          if( is.eq.2 .and. js.eq.2 ) cycle
          xj(1:3)= ra(1:3,j)
          x= xj(1) -xi(1)
          y= xj(2) -xi(2)
          z= xj(3) -xi(3)
          xij(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
          xji(1:3)= -xij(1:3)
          rij= dsqrt(sprod(xij,xij))
          R1ij= p_R1(is,js)
          D1ij= p_D1(is,js)
c---------cutoff judgement
          if(rij.gt.R1ij+D1ij) cycle
          beta= p_beta(is,js)
          s= p_s(is,js)
          D0= p_D0(is,js)
          r0= p_r0(is,js)
          riji= 1d0/rij
          dixij(1:3)= -xij(1:3)*riji
          djxij(1:3)=  xij(1:3)*riji
          djxji(1:3)= -xji(1:3)*riji
          dixji(1:3)=  xji(1:3)*riji
          aexp= dexp(-beta*dsqrt(2d0/s)*(rij-r0))
          va= fc_r(rij,R1ij,D1ij)*D0*s/(s-1d0)*aexp
          dvadr= D0*s/(s-1d0)*aexp*( dfc_r(rij,R1ij,D1ij)
     &         -beta*dsqrt(2d0/s)*fc_r(rij,R1ij,D1ij) ) 
c.....parameters for g(theta)
          gij= p_gamma(is,js)
          cij= p_c(is,js)
          dij= p_d(is,js)
          hij= p_h(is,js)
c---------make gfi
          tk= 0d0
          do kk=1,lspr(0,i)
            k= lspr(kk,i)
            if(k.eq.j) cycle
            ks= int(tag(k))
            x= ra(1,k) -xi(1)
            y= ra(2,k) -xi(2)
            z= ra(3,k) -xi(3)
            xik(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
            rik=dsqrt(sprod(xik,xik))
            R1ik= p_R1(is,ks)
            D1ik= p_D1(is,ks)
c-----------cutoff judgement
            if(rik.gt.R1ik+D1ik) cycle
            riki= 1d0/rik
            cs= sprod(xij,xik)*riji*riki
            gc= gij*(1d0 +cij**2/dij**2
     &           -(cij**2/(dij**2 +(hij+cs)**2)))
c.....alpha is always 0d0 in cases of W-W and W-He
            aexp= 1d0
            tk=tk +gc*fc_r(rik,R1ik,D1ik) *aexp
          enddo
          gfi= 1d0 +tk
c---------scan k around i to get bij
          tk= 0d0
          fi(1:3)= 0d0
          fj(1:3)= 0d0
          t1= -0.5d0*gfi**(-1.5d0)
          do kk=1,lspr(0,i)
            k= lspr(kk,i)
            if(k.eq.j) cycle
            ks=int(tag(k))
            x= ra(1,k) -xi(1)
            y= ra(2,k) -xi(2)
            z= ra(3,k) -xi(3)
            xik(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
            rik=dsqrt(sprod(xik,xik))
            R1ik= p_R1(is,ks)
            D1ik= p_D1(is,ks)
c-----------cutoff judgement
            if(rik.gt.R1ik+D1ik) cycle
            riki= 1d0/rik
            dixik(1:3)= -xik(1:3)*riki
            dkxik(1:3)=  xik(1:3)*riki
            cs= sprod(xij,xik)*riji*riki
            dics(1:3)= -riji*riki*(xij(1:3)+xik(1:3))
     &           -cs*(dixij(1:3)*riji+dixik(1:3)*riki)
            djcs(1:3)= riji*riki*xik(1:3) -cs*djxij(1:3)*riji
            dkcs(1:3)= riji*riki*xij(1:3) -cs*dkxik(1:3)*riki
            gc= gij*(1d0 +cij**2/dij**2 -(cij**2/(dij**2
     &           +(hij+cs)**2)))
            dgc= 2d0*(hij+cs)*gij*cij**2/(dij**2+(hij+cs)**2)**2
            frik = fc_r(rik,R1ik,D1ik)
            dfrik= dfc_r(rik,R1ik,D1ik)
            tk=tk +gc*frik
            fi(1:3)=fi(1:3) +dics(1:3)*dgc*frik
     &           +gc*dfrik*dixik(1:3)
            fj(1:3)=fj(1:3) +djcs(1:3)*dgc*frik
            aa2(1:3,k)=aa2(1:3,k) -0.5d0*va*t1*( dkcs(1:3)*dgc*frik
     &           +gc*dfrik*dkxik(1:3) )
          enddo
c---------derivative of bij part
          aa2(1:3,i)=aa2(1:3,i) -0.5d0*va*t1*fi(1:3)
          aa2(1:3,j)=aa2(1:3,j) -0.5d0*va*t1*fj(1:3)
c---------potential
          bij= gfi**(-0.5d0)
          epi(i)=epi(i) -0.5d0*0.5d0*bij*va
          epi(j)=epi(j) -0.5d0*0.5d0*bij*va
          epotl=epotl -0.5d0*bij*va
c---------derivative of va part
          aa2(1:3,i)=aa2(1:3,i) -0.5d0*bij*dvadr*dixij(1:3)
          aa2(1:3,j)=aa2(1:3,j) -0.5d0*bij*dvadr*djxij(1:3)
c---------if j.gt.natm, no need to calc the term around j
          if(j.gt.natm) cycle
c---------make gfj
          tk= 0d0
          do kk=1,lspr(0,j)
            k= lspr(kk,j)
            if(k.eq.i) cycle
            ks=int(tag(k))
            x= ra(1,k) -xj(1)
            y= ra(2,k) -xj(2)
            z= ra(3,k) -xj(3)
            xjk(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
            rjk=dsqrt(sprod(xjk,xjk))
            R1jk= p_R1(js,ks)
            D1jk= p_D1(js,ks)
c-----------cutoff judgement
            if(rjk.gt.R1jk+D1jk) cycle
            rjki= 1d0/rjk
            cs= sprod(xji,xjk)*riji*rjki
            gc= gij*(1d0 +cij**2/dij**2
     &           -(cij**2/(dij**2 +(hij+cs)**2)))
c.....alpha is always 0d0 in cases of W-W and W-He
            aexp= 1d0
            tk=tk +gc*fc_r(rjk,R1jk,D1jk) *aexp
          enddo
          gfj= 1d0 +tk
c---------scan k around j to get bji
          tk= 0d0
          fi(1:3)= 0d0
          fj(1:3)= 0d0
          t1= -0.5d0*gfj**(-1.5d0)
          do kk=1,lspr(0,j)
            k= lspr(kk,j)
            if(k.eq.i) cycle
            ks=int(tag(k))
            x= ra(1,k) -xj(1)
            y= ra(2,k) -xj(2)
            z= ra(3,k) -xj(3)
            xjk(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
            rjk=dsqrt(sprod(xjk,xjk))
            R1jk= p_R1(js,ks)
            D1jk= p_D1(js,ks)
c-----------cutoff judgement
            if(rjk.gt.R1jk+D1jk) cycle
            rjki= 1d0/rjk
            djxjk(1:3)= -xjk(1:3)*rjki
            dkxjk(1:3)=  xjk(1:3)*rjki
            cs= sprod(xji,xjk)*riji*rjki
            djcs(1:3)= -riji*rjki*(xji(1:3)+xjk(1:3))
     &           -cs*(djxji(1:3)*riji+djxjk(1:3)*rjki)
            dics(1:3)= riji*rjki*xjk(1:3) -cs*dixji(1:3)*riji
            dkcs(1:3)= riji*rjki*xji(1:3) -cs*dkxjk(1:3)*rjki
            gc= gij*(1d0 +cij**2/dij**2 -(cij**2/(dij**2
     &           +(hij+cs)**2)))
            dgc= 2d0*(hij+cs)*gij*cij**2/(dij**2+(hij+cs)**2)**2
            frjk = fc_r(rjk,R1jk,D1jk)
            dfrjk= dfc_r(rjk,R1jk,D1jk)
            tk=tk +gc*frjk
            fj(1:3)=fj(1:3) +djcs(1:3)*dgc*frjk
     &           +gc*dfrjk*djxjk(1:3)
            fi(1:3)=fi(1:3) +dics(1:3)*dgc*frjk
            aa2(1:3,k)=aa2(1:3,k) -0.5d0*va*t1*( dkcs(1:3)*dgc*frjk
     &           +gc*dfrjk*dkxjk(1:3))
          enddo
c---------derivative of bji part
          aa2(1:3,i)=aa2(1:3,i) -0.5d0*va*t1*fi(1:3)
          aa2(1:3,j)=aa2(1:3,j) -0.5d0*va*t1*fj(1:3)
c---------potential
          bji= gfj**(-0.5d0)
          epi(i)=epi(i) -0.5d0*0.5d0*bji*va
          epi(j)=epi(j) -0.5d0*0.5d0*bji*va
          epotl=epotl -0.5d0*bji*va
c---------derivative of va part
          aa2(1:3,i)=aa2(1:3,i) -0.5d0*bji*dvadr*dixij(1:3)
          aa2(1:3,j)=aa2(1:3,j) -0.5d0*bji*dvadr*djxij(1:3)
        enddo
      enddo

c-----send back forces and potentials on immigrants
      call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,lsrc,myparity
     &     ,nn,mpi_md_world,aa2,3)
      call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,lsrc,myparity
     &     ,nn,mpi_md_world,epi,1)

c-----sum
      aa(1:3,1:natm)= -aa1(1:3,1:natm) -aa2(1:3,1:natm)

c-----reduced force
      do i=1,natm
        at(1:3)= aa(1:3,i)
        aa(1:3,i)= hi(1:3,1)*at(1) +hi(1:3,2)*at(2) +hi(1:3,3)*at(3)
      enddo
c-----multiply 0.5d0*dt**2/am(i)
      do i=1,natm
        is= int(tag(i))
        aa(1:3,i)= acon(is)*aa(1:3,i)
      enddo

c-----gather epot
      epot= 0d0
      call mpi_allreduce(epotl,epot,1,MPI_DOUBLE_PRECISION
     &     ,MPI_SUM,mpi_md_world,ierr)

      end subroutine get_force
c=======================================================================
      function fc_r(r,r1,d1)
      implicit none
      real(8),intent(in):: r,r1,d1
      real(8):: fc_r
      real(8),parameter:: pi = 3.14159265358979d0
      
      fc_r= 0d0
      if(r.lt.r1-d1) then
        fc_r= 1d0
      elseif(r1-d1.le.r .and. r.lt.r1+d1) then
        fc_r= 0.5d0*( 1d0 -sin(0.5d0*pi*(r-r1)/d1) )
      endif
      return
      end function fc_r
c=======================================================================
      function dfc_r(r,r1,d1)
c-----Derivative: df(r)/dr
      implicit none 
      real(8),intent(in):: r,r1,d1
      real(8):: dfc_r
      real(8),parameter:: pi = 3.14159265358979d0
      
      dfc_r= 0d0
      if(r1-d1.le.r .and. r.lt.r1+d1) then
        dfc_r= -0.5d0*cos(0.5d0*pi*(r-r1)/d1)*0.5d0*pi/d1
      endif
      return
      end function dfc_r
c=======================================================================
      function f_x(x,d)
      implicit none
      real(8),intent(in):: x,d
      real(8):: f_x,tmp
      
      if( x.lt.d ) then
        tmp= d/x -1d0
        f_x= dexp(-tmp*tmp)
      else
        f_x= 1d0
      endif

      return
      end function f_x
c=======================================================================
      function df_x(x,d)
      implicit none
      real(8),intent(in):: x,d
      real(8):: df_x,tmp

      if( x.lt.d ) then
        tmp= d/x -1d0
        df_x= 2d0*d/x/x *tmp *dexp(-tmp*tmp)
      else
        df_x= 0d0
      endif
      return
      end function df_x
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make pmd"
c     End:
