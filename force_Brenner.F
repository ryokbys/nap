      subroutine get_force(namax,natm,tag,ra,nnmax,aa,h,hi,tcom
     &     ,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,anxi,anyi,anzi
     &     ,mpi_md_world,myid_md,epi,epot,nismax,acon)
c-----------------------------------------------------------------------
c Parallel implementation of the Brenner potential and forces
c                                               since 2009.04.27 by R.K.
c   - Use 1st generation Brenner potential,
c       see, D.W.Brenner, PRB 42 (1990) pp.9458-9471
c-----------------------------------------------------------------------
      implicit none
      include "mpif.h"
      include "./params_au.h"
      include "./params_Brenner.h"
      integer,intent(in):: namax,natm,nnmax,nismax
      integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3)
     &     ,nn(6),mpi_md_world,myid_md
      real(8),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),sv(3,6)
     &     ,acon(nismax),tag(namax),rc,anxi,anyi,anzi
      real(8),intent(inout):: tcom
      real(8),intent(out):: aa(3,namax),epi(namax),epot

c-----locals
      integer:: i,j,jj,k,kk,ierr,is
      real(8):: rij,riji,rik,riki,rjk,rjki,rexp,aexp,dvrdr,cs,gc,tk,t1,
     &     va,dgc,frik,dfrik,frjk,dfrjk,bij,bji,b,dvadr,frij,dfrij,
     &     gfi,gfj,tmp,x,y,z,epotl
cc-----allocatables
c      real(8),allocatable,dimension(:):: xi,xj,xk,xij,xji,xik,xjk,dixij
c     &     ,djxij,dixji,djxji,fi,fj,dixik,dkxik,djxjk,dkxjk,dics,djcs
c     &     ,dkcs
      real(8):: xi(3),xj(3),xk(3),xij(3),xji(3),xik(3),xjk(3),dixij(3)
     &     ,djxij(3),dixji(3),djxji(3),fi(3),fj(3),dixik(3),djxjk(3)
     &     ,dkxik(3),dkxjk(3),dics(3),djcs(3),dkcs(3),at(3)
      real(8),allocatable,dimension(:,:):: aa1,aa2
c-----functions
      real(8),external:: f_r,df_r,sprod
c-----1st call
      logical,save:: l1st=.true.

c-----only at 1st call
      if( l1st ) then
c-------check rc
        if( int(rc).ne.int(br_r2) ) then
          write(6,'(1x,a)') "!!! Cutoff radius is not appropriate !!!"
          write(6,'(1x,a,es12.4)') "rc should be",br_r2*bohr
          call mpi_finalize(ierr)
          stop
        endif
c-------finally set l1st
        l1st=.false.
      endif

c      allocate(xi(3),xj(3),xk(3),xij(3),xji(3),xik(3),xjk(3),dixij(3)
c     &     ,djxij(3),dixji(3),djxji(3),fi(3),fj(3),dixik(3),djxjk(3)
c     &     ,dkxik(3),dkxjk(3),dics(3),djcs(3),dkcs(3))
      allocate(aa1(3,namax),aa2(3,namax),lspr(0:nnmax,natm))

c-----initialize
      aa1(1:3,1:namax)= 0d0
      aa2(1:3,1:namax)= 0d0
      epotl= 0d0
      epi(1:namax)= 0d0

c-----make pair list for 2-body term
      call mk_lspr(namax,natm,nb,nnmax,tag,ra,rc,h,hi
     &     ,anxi,anyi,anzi,lspr)

c-----Repulsive term: V_R
      do i=1,natm
        xi(1:3)= ra(1:3,i)
        do jj=1,lspr(0,i)
          j= lspr(jj,i)
          if(j.le.i) cycle
          x= ra(1,j) -xi(1)
          y= ra(2,j) -xi(2)
          z= ra(3,j) -xi(3)
          xij(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
          rij= dsqrt(sprod(xij,xij))
cc---------cutoff judgement
c          if(rij.gt.br_r2) cycle
c---------potential
          rexp= dexp(-br_bet*dsqrt(2d0*br_s)*(rij-br_re))
          tmp=0.5d0*f_r(rij,br_r1,br_r2)*br_dd/(br_s-1d0)*rexp
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
          dvrdr= br_dd/(br_s-1d0)*rexp*( df_r(rij,br_r1,br_r2)
     &         -br_bet*dsqrt(2d0*br_s)*f_r(rij,br_r1,br_r2) )
          aa1(1:3,i)= aa1(1:3,i) +dvrdr*dixij(1:3)
          aa1(1:3,j)= aa1(1:3,j) +dvrdr*djxij(1:3)
        enddo
      enddo

c-----Attractive term: -Bij*V_A
      do i=1,natm
        xi(1:3)= ra(1:3,i)
        do jj=1,lspr(0,i)
          j= lspr(jj,i)
          if(j.le.i) cycle
          xj(1:3)= ra(1:3,j)
          x= xj(1) -xi(1)
          y= xj(2) -xi(2)
          z= xj(3) -xi(3)
          xij(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
          xji(1:3)= -xij(1:3)
          rij= dsqrt(sprod(xij,xij))
cc---------cutoff judgement
c          if(rij.gt.br_r2) cycle
          riji= 1d0/rij
          dixij(1:3)= -xij(1:3)*riji
          djxij(1:3)=  xij(1:3)*riji
          djxji(1:3)= -xji(1:3)*riji
          dixji(1:3)=  xji(1:3)*riji
          aexp= dexp(-br_bet*dsqrt(2d0/br_s)*(rij-br_re))
          va= f_r(rij,br_r1,br_r2)*br_dd*br_s/(br_s-1d0)*aexp
          dvadr= br_dd*br_s/(br_s-1d0)*aexp*( df_r(rij,br_r1,br_r2)
     &         -br_bet*dsqrt(2d0/br_s)*f_r(rij,br_r1,br_r2) ) 
c---------make gfi
          tk= 0d0
          do kk=1,lspr(0,i)
            k= lspr(kk,i)
            if(k.eq.j) cycle
            x= ra(1,k) -xi(1)
            y= ra(2,k) -xi(2)
            z= ra(3,k) -xi(3)
            xik(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
            rik=dsqrt(sprod(xik,xik))
cc-----------cutoff judgement
c            if(rik.gt.br_r2) cycle
            riki= 1d0/rik
            cs= sprod(xij,xik)*riji*riki
            gc= br_a0*(1d0 +br_c0**2/br_d0**2
     &           -(br_c0**2/(br_d0**2 +(1d0+cs)**2)))
            tk=tk +gc*f_r(rik,br_r1,br_r2)
          enddo
          gfi= 1d0 +tk
c---------scan k around i to get bij
          tk= 0d0
          fi(1:3)= 0d0
          fj(1:3)= 0d0
          t1= -br_del*gfi**(-br_del-1d0)
          do kk=1,lspr(0,i)
            k= lspr(kk,i)
            if(k.eq.j) cycle
            x= ra(1,k) -xi(1)
            y= ra(2,k) -xi(2)
            z= ra(3,k) -xi(3)
            xik(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
            rik=dsqrt(sprod(xik,xik))
cc-----------cutoff judgement
c            if(rik.gt.br_r2) cycle
            riki= 1d0/rik
            dixik(1:3)= -xik(1:3)*riki
            dkxik(1:3)=  xik(1:3)*riki
            cs= sprod(xij,xik)*riji*riki
            dics(1:3)= -riji*riki*(xij(1:3)+xik(1:3))
     &           -cs*(dixij(1:3)*riji+dixik(1:3)*riki)
            djcs(1:3)= riji*riki*xik(1:3) -cs*djxij(1:3)*riji
            dkcs(1:3)= riji*riki*xij(1:3) -cs*dkxik(1:3)*riki
            gc= br_a0*(1d0 +br_c0**2/br_d0**2 -(br_c0**2/(br_d0**2
     &           +(1d0+cs)**2)))
            dgc= 2d0*(1d0+cs)*br_a0*br_c0**2/(br_d0**2+(1d0+cs)**2)**2
            frik = f_r(rik,br_r1,br_r2)
            dfrik= df_r(rik,br_r1,br_r2)
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
          bij= gfi**(-br_del)
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
            x= ra(1,k) -xj(1)
            y= ra(2,k) -xj(2)
            z= ra(3,k) -xj(3)
            xjk(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
            rjk=dsqrt(sprod(xjk,xjk))
cc-----------cutoff judgement
c            if(rjk.gt.br_r2) cycle
            rjki= 1d0/rjk
            cs= sprod(xji,xjk)*riji*rjki
            gc= br_a0*(1d0 +br_c0**2/br_d0**2
     &           -(br_c0**2/(br_d0**2 +(1d0+cs)**2)))
            tk=tk +gc*f_r(rjk,br_r1,br_r2)
          enddo
          gfj= 1d0 +tk
c---------scan k around j to get bji
          tk= 0d0
          fi(1:3)= 0d0
          fj(1:3)= 0d0
          t1= -br_del*gfj**(-br_del-1d0)
          do kk=1,lspr(0,j)
            k= lspr(kk,j)
            if(k.eq.i) cycle
            x= ra(1,k) -xj(1)
            y= ra(2,k) -xj(2)
            z= ra(3,k) -xj(3)
            xjk(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
            rjk=dsqrt(sprod(xjk,xjk))
cc-----------cutoff judgement
c            if(rjk.gt.br_r2) cycle
            rjki= 1d0/rjk
            djxjk(1:3)= -xjk(1:3)*rjki
            dkxjk(1:3)=  xjk(1:3)*rjki
            cs= sprod(xji,xjk)*riji*rjki
            djcs(1:3)= -riji*rjki*(xji(1:3)+xjk(1:3))
     &           -cs*(djxji(1:3)*riji+djxjk(1:3)*rjki)
            dics(1:3)= riji*rjki*xjk(1:3) -cs*dixji(1:3)*riji
            dkcs(1:3)= riji*rjki*xji(1:3) -cs*dkxjk(1:3)*rjki
            gc= br_a0*(1d0 +br_c0**2/br_d0**2 -(br_c0**2/(br_d0**2
     &           +(1d0+cs)**2)))
            dgc= 2d0*(1d0+cs)*br_a0*br_c0**2/(br_d0**2+(1d0+cs)**2)**2
            frjk = f_r(rjk,br_r1,br_r2)
            dfrjk= df_r(rjk,br_r1,br_r2)
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
          bji= gfj**(-br_del)
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

c      deallocate(xi,xj,xk,xij,xji,xik,xjk,dixij,djxij,dixji,djxji,fi,fj
c     &     ,dixik,djxjk,dkxik,dkxjk,dics,djcs,dkcs)
      deallocate(aa1,aa2)
      end subroutine get_force
c=======================================================================
      subroutine force_brenner_vdW(namax,natm,tag,ra,nnmax,lspr,aa
     &     ,h,hi,tcom,nb,nbmax,lsb,lsrc,myparity,nn,sv,mpi_md_world
     &     ,myid_md,epi,epot,nismax,acon)
c-----------------------------------------------------------------------
c Parallel implementation of the Brenner potential and forces
c                                               since 2009.04.27 by R.K.
c   - Use 1st generation Brenner potential,
c       see, D.W.Brenner, PRB 42 (1990) pp.9458-9471
c   - Include van der Waals potential for interlayer interaction,
c       see, A.J.Heim et.al, cond-mat (25 Jun 2008)
c   - Note: this mixed potential returns wrong C-C bond length
c-----------------------------------------------------------------------
      implicit none
      include "mpif.h"
      include "./params_au.h"
      include "./params_Brenner.h"
      integer,intent(in):: namax,natm,nnmax,lspr(0:nnmax,namax),nismax
      integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3)
     &     ,nn(6),mpi_md_world,myid_md
      real(8),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),sv(3,6)
     &     ,acon(nismax),tag(namax)
      real(8),intent(inout):: tcom
      real(8),intent(out):: aa(3,namax),epi(namax),epot

c-----max num of neighbors for Brenner potential
      integer,parameter:: nnmaxb = 20
c-----locals
      integer:: i,j,jj,k,kk,ierr,is
      real(8):: rij,riji,rik,riki,rjk,rjki,rexp,aexp,dvrdr,cs,gc,tk,t1,
     &     va,dgc,frik,dfrik,frjk,dfrjk,bij,bji,b,dvadr,frij,dfrij,
     &     gfi,gfj,tmp,x,y,z,epotl
cc-----allocatables
c      real(8),allocatable,dimension(:):: xi,xj,xk,xij,xji,xik,xjk,dixij
c     &     ,djxij,dixji,djxji,fi,fj,dixik,dkxik,djxjk,dkxjk,dics,djcs
c     &     ,dkcs
      real(8):: xi(3),xj(3),xk(3),xij(3),xji(3),xik(3),xjk(3),dixij(3)
     &     ,djxij(3),dixji(3),djxji(3),fi(3),fj(3),dixik(3),djxjk(3)
     &     ,dkxik(3),dkxjk(3),dics(3),djcs(3),dkcs(3),at(3)
      real(8),allocatable,dimension(:,:):: aa1,aa2
      integer,allocatable:: lsprb(:,:)
c-----functions
      real(8),external:: f_r,df_r,sprod
      real(8),external:: vvdw,dvvdw

c      allocate(xi(3),xj(3),xk(3),xij(3),xji(3),xik(3),xjk(3),dixij(3)
c     &     ,djxij(3),dixji(3),djxji(3),fi(3),fj(3),dixik(3),djxjk(3)
c     &     ,dkxik(3),dkxjk(3),dics(3),djcs(3),dkcs(3))
      allocate(aa1(3,namax),aa2(3,namax),lsprb(0:nnmaxb,namax))

c-----initialize
      aa1(1:3,1:namax)= 0d0
      aa2(1:3,1:namax)= 0d0
      epotl= 0d0
      epi(1:namax)= 0d0
      lsprb(0:nnmaxb,namax)= 0

c-----van der Waals potential with making LSPRB
      do i=1,natm
        xi(1:3)= ra(1:3,i)
        do jj=1,lspr(0,i)
          j= lspr(jj,i)
          if(j.le.i) cycle
          x= ra(1,j) -xi(1)
          y= ra(2,j) -xi(2)
          z= ra(3,j) -xi(3)
          xij(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
          rij= dsqrt(sprod(xij,xij))
          riji= 1d0/rij
          dixij(1:3)= -xij(1:3)*riji
          djxij(1:3)=  xij(1:3)*riji
c---------potential
          tmp= 0.5d0*vvdw(rij,riji)
          epi(i)=epi(i) +tmp
          if(j.le.natm) then
            epi(j)=epi(j) +tmp
            epotl=epotl +tmp +tmp
          else
            epotl=epotl +tmp
          endif
c---------force
          tmp= dvvdw(rij,riji)
          aa1(1:3,i)=aa1(1:3,i) +dixij(1:3)*tmp
          aa1(1:3,j)=aa1(1:3,j) +djxij(1:3)*tmp
c---------make LSPRB
          if(rij.le.br_r2) then
            lsprb(0,i)=lsprb(0,i)+1
            lsprb(0,j)=lsprb(0,j)+1
c-----------error trap
            if(lsprb(0,i).gt.nnmaxb
     &           .or.lsprb(0,j).gt.nnmaxb) then
              write(6,'(a)') "!!!lsprb(0,i).gt.nnmaxb!!!"
              call mpi_finalize(ierr)
              stop
            endif
            lsprb(lsprb(0,i),i)= j
            lsprb(lsprb(0,j),j)= i
          endif
        enddo
      enddo

c-----Repulsive term: V_R
      do i=1,natm
        xi(1:3)= ra(1:3,i)
        do jj=1,lsprb(0,i)
          j= lsprb(jj,i)
          if(j.le.i) cycle
          x= ra(1,j) -xi(1)
          y= ra(2,j) -xi(2)
          z= ra(3,j) -xi(3)
          xij(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
          rij= dsqrt(sprod(xij,xij))
cc---------cutoff judgement
c          if(rij.gt.br_r2) cycle
c---------potential
          rexp= dexp(-br_bet*dsqrt(2d0*br_s)*(rij-br_re))
          tmp=0.5d0*f_r(rij,br_r1,br_r2)*br_dd/(br_s-1d0)*rexp
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
          dvrdr= br_dd/(br_s-1d0)*rexp*( df_r(rij,br_r1,br_r2)
     &         -br_bet*dsqrt(2d0*br_s)*f_r(rij,br_r1,br_r2) )
          aa1(1:3,i)= aa1(1:3,i) +dvrdr*dixij(1:3)
          aa1(1:3,j)= aa1(1:3,j) +dvrdr*djxij(1:3)
        enddo
      enddo

c-----Attractive term: -Bij*V_A
      do i=1,natm
        xi(1:3)= ra(1:3,i)
        do jj=1,lsprb(0,i)
          j= lsprb(jj,i)
          if(j.le.i) cycle
          xj(1:3)= ra(1:3,j)
          x= xj(1) -xi(1)
          y= xj(2) -xi(2)
          z= xj(3) -xi(3)
          xij(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
          xji(1:3)= -xij(1:3)
          rij= dsqrt(sprod(xij,xij))
cc---------cutoff judgement
c          if(rij.gt.br_r2) cycle
          riji= 1d0/rij
          dixij(1:3)= -xij(1:3)*riji
          djxij(1:3)=  xij(1:3)*riji
          djxji(1:3)= -xji(1:3)*riji
          dixji(1:3)=  xji(1:3)*riji
          aexp= dexp(-br_bet*dsqrt(2d0/br_s)*(rij-br_re))
          va= f_r(rij,br_r1,br_r2)*br_dd*br_s/(br_s-1d0)*aexp
          dvadr= br_dd*br_s/(br_s-1d0)*aexp*( df_r(rij,br_r1,br_r2)
     &         -br_bet*dsqrt(2d0/br_s)*f_r(rij,br_r1,br_r2) ) 
c---------make gfi
          tk= 0d0
          do kk=1,lsprb(0,i)
            k= lsprb(kk,i)
            if(k.eq.j) cycle
            x= ra(1,k) -xi(1)
            y= ra(2,k) -xi(2)
            z= ra(3,k) -xi(3)
            xik(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
            rik=dsqrt(sprod(xik,xik))
cc-----------cutoff judgement
c            if(rik.gt.br_r2) cycle
            riki= 1d0/rik
            cs= sprod(xij,xik)*riji*riki
            gc= br_a0*(1d0 +br_c0**2/br_d0**2
     &           -(br_c0**2/(br_d0**2 +(1d0+cs)**2)))
            tk=tk +gc*f_r(rik,br_r1,br_r2)
          enddo
          gfi= 1d0 +tk
c---------scan k around i to get bij
          tk= 0d0
          fi(1:3)= 0d0
          fj(1:3)= 0d0
          t1= -br_del*gfi**(-br_del-1d0)
          do kk=1,lsprb(0,i)
            k= lsprb(kk,i)
            if(k.eq.j) cycle
            x= ra(1,k) -xi(1)
            y= ra(2,k) -xi(2)
            z= ra(3,k) -xi(3)
            xik(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
            rik=dsqrt(sprod(xik,xik))
cc-----------cutoff judgement
c            if(rik.gt.br_r2) cycle
            riki= 1d0/rik
            dixik(1:3)= -xik(1:3)*riki
            dkxik(1:3)=  xik(1:3)*riki
            cs= sprod(xij,xik)*riji*riki
            dics(1:3)= -riji*riki*(xij(1:3)+xik(1:3))
     &           -cs*(dixij(1:3)*riji+dixik(1:3)*riki)
            djcs(1:3)= riji*riki*xik(1:3) -cs*djxij(1:3)*riji
            dkcs(1:3)= riji*riki*xij(1:3) -cs*dkxik(1:3)*riki
            gc= br_a0*(1d0 +br_c0**2/br_d0**2 -(br_c0**2/(br_d0**2
     &           +(1d0+cs)**2)))
            dgc= 2d0*(1d0+cs)*br_a0*br_c0**2/(br_d0**2+(1d0+cs)**2)**2
            frik = f_r(rik,br_r1,br_r2)
            dfrik= df_r(rik,br_r1,br_r2)
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
          bij= gfi**(-br_del)
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
          do kk=1,lsprb(0,j)
            k= lsprb(kk,j)
            if(k.eq.i) cycle
            x= ra(1,k) -xj(1)
            y= ra(2,k) -xj(2)
            z= ra(3,k) -xj(3)
            xjk(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
            rjk=dsqrt(sprod(xjk,xjk))
cc-----------cutoff judgement
c            if(rjk.gt.br_r2) cycle
            rjki= 1d0/rjk
            cs= sprod(xji,xjk)*riji*rjki
            gc= br_a0*(1d0 +br_c0**2/br_d0**2
     &           -(br_c0**2/(br_d0**2 +(1d0+cs)**2)))
            tk=tk +gc*f_r(rjk,br_r1,br_r2)
          enddo
          gfj= 1d0 +tk
c---------scan k around j to get bji
          tk= 0d0
          fi(1:3)= 0d0
          fj(1:3)= 0d0
          t1= -br_del*gfj**(-br_del-1d0)
          do kk=1,lsprb(0,j)
            k= lsprb(kk,j)
            if(k.eq.i) cycle
            x= ra(1,k) -xj(1)
            y= ra(2,k) -xj(2)
            z= ra(3,k) -xj(3)
            xjk(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
            rjk=dsqrt(sprod(xjk,xjk))
cc-----------cutoff judgement
c            if(rjk.gt.br_r2) cycle
            rjki= 1d0/rjk
            djxjk(1:3)= -xjk(1:3)*rjki
            dkxjk(1:3)=  xjk(1:3)*rjki
            cs= sprod(xji,xjk)*riji*rjki
            djcs(1:3)= -riji*rjki*(xji(1:3)+xjk(1:3))
     &           -cs*(djxji(1:3)*riji+djxjk(1:3)*rjki)
            dics(1:3)= riji*rjki*xjk(1:3) -cs*dixji(1:3)*riji
            dkcs(1:3)= riji*rjki*xji(1:3) -cs*dkxjk(1:3)*rjki
            gc= br_a0*(1d0 +br_c0**2/br_d0**2 -(br_c0**2/(br_d0**2
     &           +(1d0+cs)**2)))
            dgc= 2d0*(1d0+cs)*br_a0*br_c0**2/(br_d0**2+(1d0+cs)**2)**2
            frjk = f_r(rjk,br_r1,br_r2)
            dfrjk= df_r(rjk,br_r1,br_r2)
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
          bji= gfj**(-br_del)
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

c      deallocate(xi,xj,xk,xij,xji,xik,xjk,dixij,djxij,dixji,djxji,fi,fj
c     &     ,dixik,djxjk,dkxik,dkxjk,dics,djcs,dkcs)
      deallocate(aa1,aa2,lsprb)
      end subroutine force_brenner_vdW
c=======================================================================
      function f_r(r,r1,r2)
      implicit none
      real(8),intent(in):: r,r1,r2
      real(8):: f_r
      real(8),parameter:: pi = 3.14159265358979d0
      
      f_r= 0d0
      if(r.lt.r1) then
        f_r= 1d0
      elseif(r1.le.r .and. r.lt.r2) then
        f_r= 0.5d0*(1d0 +cos(pi*(r-r1)/(r2-r1)))
      endif
      return
      end function f_r
c=======================================================================
      function df_r(r,r1,r2)
c-----Derivative: df(r)/dr
      implicit none 
      real(8),intent(in):: r,r1,r2
      real(8):: df_r
      real(8),parameter:: pi = 3.14159265358979d0
      
      df_r= 0d0
      if(r1.le.r .and. r.lt.r2) then
        df_r= -0.5d0*sin(pi*(r-r1)/(r2-r1))*pi/(r2-r1)
      endif
      return
      end function df_r
c=======================================================================
      function vvdw(r,ri)
c-----------------------------------------------------------------------
c Graphene interlayer potential
c   See, A.J.Heim et.al, cond-mat, (25 Jun 2008)
c-----------------------------------------------------------------------
      implicit none
      include "./params_au.h"
      include "./params_Brenner.h"
      real(8),intent(in):: r,ri
      real(8):: vvdw

      vvdw= 0d0
      if(r.lt.r1vdw) then
        vvdw= evdw*cvdw
      elseif(r.ge.r1vdw .and. r.lt.r2vdw) then
        vvdw= evdw*(cvdw +p3vdw*(r-r1vdw)**3
     &       +p4vdw*(r-r1vdw)**4 +p5vdw*(r-r1vdw)**5)
      elseif(r.ge.r2vdw .and. r.lt.r3vdw) then
        vvdw= evdw*( (sgmvdw*ri)**12 -2d0*(sgmvdw*ri)**6 )
      elseif(r.ge.r3vdw .and. r.lt.r4vdw) then
        vvdw= evdw*( c0vdw +c1vdw*r +c2vdw*r**2 +c3vdw*r**3 )
      endif
      return
      end function vvdw
c=======================================================================
      function dvvdw(r,ri)
c-----------------------------------------------------------------------
c Derivative of the graphene interlayer potential
c-----------------------------------------------------------------------
      implicit none
      include "./params_au.h"
      include "./params_Brenner.h"
      real(8),intent(in):: r,ri
      real(8):: dvvdw

      dvvdw= 0d0
      if(r.ge.r1vdw .and. r.lt.r2vdw) then
        dvvdw= evdw*(3d0*p3vdw*(r-r1vdw)**2
     &       +4d0*p4vdw*(r-r1vdw)**3 +5d0*p5vdw*(r-r1vdw)**4)
      elseif(r.ge.r2vdw .and. r.lt.r3vdw) then
        dvvdw= -12d0*evdw*( (sgmvdw*ri)**12 -(sgmvdw*ri)**6 ) *ri
      elseif(r.ge.r3vdw .and. r.lt.r4vdw) then
        dvvdw= evdw*( c1vdw +2d0*c2vdw*r +3d0*c3vdw*r**2 )
      endif
      return
      end function dvvdw
