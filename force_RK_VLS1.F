      subroutine get_force(namax,natm,tag,ra,nnmax,aa,h,hi,tcom
     &     ,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,anxi,anyi,anzi
     &     ,mpi_md_world,myid_md,epi,epot,nismax,acon)
c-----------------------------------------------------------------------
c  Parallel implementation of 2 species LJ force calculation
c    - only force on i is considered, no need to send back
c-----------------------------------------------------------------------
      implicit none
      include "mpif.h"
      include "./params_unit.h"
      include "params_RK_VLS1.h"
      integer,intent(in):: namax,natm,nnmax,nismax
      integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3)
     &     ,nn(6),mpi_md_world,myid_md
      real(8),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),rc
     &     ,acon(nismax),tag(namax),sv(3,6),anxi,anyi,anzi
      real(8),intent(inout):: tcom
      real(8),intent(out):: aa(3,namax),epi(namax),epot

      integer:: i,j,k,l,m,n,ierr,is,js,ir
      real(8):: vrc,dvdrc,xi(3),xx(3),xij(3),drij(3),rij,ri
     &     ,epotl,at(3),tmp,d,v2,dv2,f2rc,df2rc
      integer,allocatable:: lspr(:,:)
c-----2-body force table (allocated at 1st call and reuse them)
      real(8),allocatable,save:: tf(:,:,:),tdf(:,:,:)
      real(8),save:: rmin,rmax,dr
c-----1st call
      logical,save:: l1st=.true.
c-----function
      real(8),external:: f2_rk,df2_rk

c-----only at 1st call
      if( l1st ) then
c-------allocate 2-body force table
        allocate(tf(nd_rk,2,2),tdf(nd_rk,2,2))
c-------make 2-body (smoothed) force table
        rmin= 0.5d0
        rmax= rc
        dr= (rmax-rmin)/(nd_rk-1)
        do is=1,2
          do js=1,2
            f2rc = f2_rk(rc,is,js)
            df2rc=df2_rk(rc,is,js)
            do i=1,nd_rk
              ri= rmin +dr*(i-1)
              tf(i,is,js) = f2_rk(ri,is,js) -f2rc -(ri-rc)*df2rc
              tdf(i,is,js)=df2_rk(ri,is,js) -df2rc
            enddo
          enddo
        enddo
c-------finally set l1st as false
        l1st= .false.
      endif

      aa(1:3,1:namax)=0d0
      epi(1:namax)= 0d0
      epotl= 0d0

      allocate(lspr(0:nnmax,natm))

c-----make pair list for 2-body term
      call mk_lspr(namax,natm,nb,nnmax,tag,ra,rc,h,hi
     &     ,anxi,anyi,anzi,lspr)

c-----loop over resident atoms
      do i=1,natm
        xi(1:3)= ra(1:3,i)
        is= int(tag(i))
        do k=1,lspr(0,i)
          j=lspr(k,i)
          if(j.eq.0) exit
          if(j.le.i) cycle
          js= int(tag(j))
          xx(1:3)= ra(1:3,j) -xi(1:3)
          xij(1:3)= h(1:3,1,0)*xx(1) +h(1:3,2,0)*xx(2) +h(1:3,3,0)*xx(3)
          rij= sqrt(xij(1)**2+ xij(2)**2 +xij(3)**2)
          ir= int( (rij-rmin)/dr +1 )
          d = (rij-rmin)/dr -(ir-1)
c---------potential
          v2= tf(ir,is,js) +(tf(ir+1,is,js)-tf(ir,is,js))*d
          v2= v2 /2
          epi(i)= epi(i) +v2
          epotl= epotl +v2
          if( j.le.natm ) then
            epi(j)= epi(j) +v2
            epotl = epotl +v2
          endif
c---------force
          drij(1:3)= -xij(1:3)/rij
          dv2= tdf(ir,is,js) +(tdf(ir+1,is,js)-tdf(ir,is,js))*d
          aa(1:3,i)= aa(1:3,i) -drij(1:3)*dv2
          aa(1:3,j)= aa(1:3,j) +drij(1:3)*dv2
        enddo
      enddo

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

      deallocate(lspr)
      end subroutine get_force
c=======================================================================
      function f2_rk(r,is,js)
c-----2-body force
      implicit none
      include "./params_unit.h"
      include "params_RK_VLS1.h"
      real(8),intent(in):: r
      integer,intent(in):: is,js
      
c-----value
      real(8):: f2_rk

      f2_rk= 4d0 *eps_rk(is,js) *( (sgm_rk(is,js)/r)**12
     &     -(sgm_rk(is,js)/r)**6 )

      return
      end function f2_rk
c=======================================================================
      function df2_rk(r,is,js)
c-----Derivative of 2-body term
      implicit none 
      include "./params_unit.h"
      include "params_RK_VLS1.h"
      real(8),intent(in):: r
      integer,intent(in):: is,js

      real(8):: df2_rk,dij
      
      df2_rk= -24d0 *eps_rk(is,js) *( 2d0*sgm_rk(is,js)**12/r**13
     &     -sgm_rk(is,js)**6/r**7 )

      return
      end function df2_rk
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make pmd"
c     End:
