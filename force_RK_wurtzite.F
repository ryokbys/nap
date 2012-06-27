      subroutine setup(nismax,am,acon,fack,dt)
      implicit none
      include "params_au.h"
      include "params_RK_wurtzite.h"
      integer,intent(in):: nismax
      real(8),intent(in):: dt
      real(8),intent(out):: am(nismax),acon(nismax),fack(nismax)

      integer:: i

c-----atomic masses
      am(1)= am_Al
      am(2)= am_N
      am(3:nismax)= am(1)

c-----prefactors for normalized accelerration
      do i=1,nismax
        acon(i)= 0.5d0*dt**2/am(i)
      enddo

c-----prefactors for kinetic energy, FACK
      do i=1,nismax
        fack(i)= 0.5d0*am(i)/dt**2
      enddo

cc-----check
c      write(6,'(a,10es12.4)') " am  =",am(1:nismax)
c      write(6,'(a,10es12.4)') " acon=",acon(1:nismax)
c      write(6,'(a,10es12.4)') " fack=",fack(1:nismax)

      end subroutine setup
c=======================================================================
      subroutine get_force(namax,natm,tag,ra,nnmax,aa,h,hi,tcom
     &     ,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,anxi,anyi,anzi
     &     ,mpi_world,myid,epi,epot,nismax,acon)
c-----------------------------------------------------------------------
c  Parallel implementation of RK potential for wurtzite structure
c    - 2010.04.07 by R.K.
c-----------------------------------------------------------------------
      implicit none
      include "mpif.h"
      include "params_au.h"
      include "params_RK_wurtzite.h"
      integer,intent(in):: namax,natm,nnmax,nismax
      integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3)
     &     ,nn(6),mpi_world,myid
      real(8),intent(in):: ra(3,namax),tag(namax),acon(nismax)
     &     ,h(3,3),hi(3,3),sv(3,6),rc,anxi,anyi,anzi
      real(8),intent(inout):: tcom
      real(8),intent(out):: aa(3,namax),epi(namax),epot

      integer:: i,j,k,l,m,n,ierr,is,js,ks,ir
      real(8):: f2rc,df2rc,ri,cst,cs
      real(8):: epotl,epotl2,epotl3,t1,t2,tmp1,tmp2,vol
      real(8):: rij,riji,rik,riki,d,v2,dv2,expij,expik,v3,dt1j,dt1k,dt2
      integer,allocatable:: lspr(:,:)
      real(8),allocatable,dimension(:):: xi,xj,xij,xik,at,xx,drij,drik
     &     ,dcsj,dcsk,dcsi
      real(8),allocatable:: aa2(:,:),aa3(:,:)

c-----functions
      real(8),external:: f2_r,df2_r
c-----saved values
      real(8),save:: rmin,rmax,dr
c-----saved allocatable arrays
      real(8),allocatable,save:: tblf2(:,:,:),tbldf2(:,:,:)
c-----1st call
      logical,save:: l1st=.true.

c-----only at 1st call
      if( l1st ) then
c-------allocate 2-body force table at the 1st call (do not deallocate!)
        allocate(tblf2(nd2b,2,2),tbldf2(nd2b,2,2))
c-------make 2-body (smoothed) force table
        rmin= 0.5d-10 /bohr
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
            write(91,'(a)') ""
          enddo
        enddo
        do i=1,nd2b
          ri= rmin +dr*(i-1)
          write(90,'(10es12.4)') ri*bohr,tblf2(i,1:2,1:2)
     &         ,tbldf2(i,1:2,1:2)
        enddo
        write(6,'(a,2es12.4)') " eps =",rk_eps
        write(6,'(a,4es12.4)') " b   =",rk_b
        write(6,'(a,4es12.4)') " sgm =",rk_sgm
        write(6,'(a,4es12.4)') " l1  =",rk_l1
        write(6,'(a,4es12.4)') " l2  =",rk_l2
        write(6,'(a,4es12.4)') " l3  =",rk_l3
        write(6,'(a,4es12.4)') " rc  =",rk_rc
        write(6,'(a,i8)') " nb =",nb

c-------finally set l1st
        l1st=.false.
      endif
      
      allocate(lspr(0:nnmax,natm))
      allocate(aa2(3,namax),aa3(3,namax))
      allocate(xi(3),xij(3),at(3),xx(3),xik(3),drij(3),drik(3)
     &     ,dcsj(3),dcsk(3),dcsi(3))

      epotl= 0d0
      epi(1:namax)= 0d0
      aa2(1:3,1:namax)= 0d0
      aa3(1:3,1:namax)= 0d0
      epotl2= 0d0
      epotl3= 0d0

c-----make pair list for 2-body term
      call mk_lspr(namax,natm,nb,nnmax,tag,ra,rc,h,hi
     &     ,anxi,anyi,anzi,lspr)

c-----2-body term
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
c---------potential
          v2= tblf2(ir,is,js) +(tblf2(ir+1,is,js)-tblf2(ir,is,js))*d
          v2= v2 /2
          epi(i)= epi(i) +v2
          epotl2= epotl2 +v2
          if( j.le.natm ) then
            epi(j)= epi(j) +v2
            epotl2= epotl2 +v2
          endif
c---------force
          drij(1:3)= -xij(1:3)/rij
          dv2= tbldf2(ir,is,js)
     &         +(tbldf2(ir+1,is,js)-tbldf2(ir,is,js))*d
          aa2(1:3,i)= aa2(1:3,i) +drij(1:3)*dv2
          aa2(1:3,j)= aa2(1:3,j) -drij(1:3)*dv2
        enddo
      enddo

c-----make pair list for 3-body term
      call mk_lspr(namax,natm,nb,nnmax,tag,ra,rk_rc,h,hi
     &     ,anxi,anyi,anzi,lspr)

c-----3-body term
      cst= -1d0/3d0
      do i=1,natm
        xi(1:3)= ra(1:3,i)
        is= int(tag(i))
        tmp1= 0d0
        tmp2= 0d0
        do m=1,lspr(0,i)
          j=lspr(m,i)
          js= int(tag(j))
          xx(1:3)= ra(1:3,j) -xi(1:3)
          xij(1:3)= h(1:3,1)*xx(1) +h(1:3,2)*xx(2) +h(1:3,3)*xx(3)
          rij= dsqrt(xij(1)*xij(1) +xij(2)*xij(2) +xij(3)*xij(3))
          riji= 1d0/rij
          drij(1:3)= -xij(1:3)*riji
          expij= exp(rk_l3/(rij-rk_rc))
          do n=1,lspr(0,i)
            k=lspr(n,i)
            if( k.le.j ) cycle
            xx(1:3)= ra(1:3,k) -xi(1:3)
            xik(1:3)= h(1:3,1)*xx(1) +h(1:3,2)*xx(2) +h(1:3,3)*xx(3)
            rik=dsqrt(xik(1)*xik(1)+xik(2)*xik(2)+xik(3)*xik(3))
            riki= 1d0/rik
            drik(1:3)= -xik(1:3)*riki
            expik= exp(rk_l3/(rik-rk_rc))
c-----------potential
            t1= rk_b*expij*expik
            cs= (xij(1)*xik(1)+xij(2)*xik(2)+xij(3)*xik(3))*riji*riki
            t2= (cs-cst)**2 /(1d0+rk_c*(cs-cst)**2)
            tmp2= tmp2 +t2
            tmp1= tmp1 +t1
            v3= t1*t2
            epi(i)= epi(i) +v3
            epotl3= epotl3 +v3
c            write(6,'(a,3i6,4es13.4e3)') "i,j,k,t1,t2,v3=",i,j,k,t1,t2
c     &           ,v3,epotl3
c-----------force
            dt1j= rk_b*expij*expik*(-rk_l3/(rij-rk_rc)**2)
            dt1k= rk_b*expij*expik*(-rk_l3/(rik-rk_rc)**2)
            dt2= 2d0*(cs-cst) /( 1d0 +rk_c*(cs-cst)**2 )**2
            dcsj(1:3)= xik(1:3)*riji*riki -xij(1:3)*cs*riji**2
            dcsk(1:3)= xij(1:3)*riji*riki -xik(1:3)*cs*riki**2
            dcsi(1:3)= -dcsj(1:3) -dcsk(1:3)
            aa3(1:3,i)=aa3(1:3,i) +t2*( dt1j*drij(1:3)+dt1k*drik(1:3) )
     &           +t1*dt2*dcsi(1:3)
            aa3(1:3,j)=aa3(1:3,j) +t2*dt1j*(-drij(1:3))
     &           +t1*dt2*dcsj(1:3)
            aa3(1:3,k)=aa3(1:3,k) +t2*dt1k*(-drik(1:3))
     &           +t1*dt2*dcsk(1:3)
          enddo
        enddo
c        write(6,'(a,i8,3es12.4)') " i,tmp1,tmp2=",i,tmp1,tmp2,epotl3
      enddo
c      stop

c-----send back (3-body) forces and potentials on immigrants
      call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,lsrc,myparity
     &     ,nn,mpi_world,aa3,3)
      call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,lsrc,myparity
     &     ,nn,mpi_world,epi,1)

c-----sum
      aa(1:3,1:natm)= -aa2(1:3,1:natm) -aa3(1:3,1:natm)

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
      epotl= epotl2 +epotl3
      call mpi_allreduce(epotl,epot,1,MPI_DOUBLE_PRECISION
     &     ,MPI_SUM,mpi_world,ierr)
c      write(6,'(a,2es13.4e3)') " epot,epotl3=",epot,epotl3

c-----get min bond length
      tmp1= 1d10
      do k=1,lspr(0,1)
        j=lspr(k,1)
        xi(1:3)= ra(1:3,j) -ra(1:3,1)
        xij(1:3)= h(1:3,1)*xi(1) +h(1:3,2)*xi(2) +h(1:3,3)*xi(3)
        rij= sqrt(xij(1)**2 +xij(2)**2 +xij(3)**2)
        tmp1= min(tmp1,rij)
      enddo
c-----output lattice constant
      vol= h(1,1)*h(2,2)*h(3,3) *0.529177d0**3
      write(92,'(10es12.4)') tmp1*0.529177d0,vol/natm,epot/natm

      deallocate(aa2,aa3,lspr,xi,xij,xik,at,xx,drij,drik,dcsj,dcsk,dcsi)
      end subroutine get_force
c=======================================================================
      function f2_r(r,is,js)
c-----2-body force
      implicit none
      include "params_au.h"
      include "params_RK_wurtzite.h"
      real(8),intent(in):: r
      integer,intent(in):: is,js
      
c-----value
      real(8):: f2_r

      f2_r=
     &     rk_eps*( exp(-2d0*(r-rk_sgm)/rk_l1)
     &       -2d0*exp(-(r-rk_sgm)/rk_l1) )
     &     +rk_z(is)*rk_z(js)/r *exp(-r/rk_l2)

      write(91,'(2i6,5es12.4)') is,js,r*bohr
     &     ,rk_eps*(exp(-2d0*(r-rk_sgm)/rk_l1)
     &       -2d0*exp(-(r-rk_sgm)/rk_l1))
     &     ,rk_z(is)*rk_z(js)/r *exp(-r/rk_l2)

      return
      end function f2_r
c=======================================================================
      function df2_r(r,is,js)
c-----Derivative of 2-body term
      implicit none 
      include "params_au.h"
      include "params_RK_wurtzite.h"
      real(8),intent(in):: r
      integer,intent(in):: is,js

      real(8):: df2_r,dij
      
      df2_r=
     &     -2d0*rk_eps/rk_l1
     &      *( exp(-2d0*(r-rk_sgm)/rk_l1)
     &       -exp(-(r-rk_sgm)/rk_l1) )
     &     -rk_z(is)*rk_z(js)/r *exp(-r/rk_l2)
     &       *(1d0/r +1d0/rk_l2)

      return
      end function df2_r
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make pmd"
c     End:
