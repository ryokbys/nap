      subroutine get_force(namax,natm,tag,ra,nnmax,aa,strs,h,hi,tcom
     &     ,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
     &     ,mpi_md_world,myid_md,epi,epot,nismax,acon)
c-----------------------------------------------------------------------
c  Parallel implementation of EAM force calculation for Al
c    - smoothing is applied to both 2- and many-body terms
c    - rho of boundary atoms are sent to the neighbor nodes
c    - only force on i is calculated, not necessary to send-back
c-----------------------------------------------------------------------
c  See Streitz and Mintmire, PRB 50(16), 11996 (1994).
c  This is implementation of only Al-Al potential in the literature,
c  and non-variable charge implementation.
c-----------------------------------------------------------------------
      implicit none
      include "mpif.h"
      include "./params_au.h"
      include "params_EAM_Al.h"
      integer,intent(in):: namax,natm,nnmax,nismax
      integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3)
     &     ,nn(6),mpi_md_world,myid_md
      integer,intent(inout):: lspr(0:nnmax,namax)
      real(8),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),sv(3,6)
     &     ,acon(nismax),rc,tag(namax)
      real(8),intent(inout):: tcom
      real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)

      integer:: i,j,k,l,m,n,ierr,is
      real(8):: xij(3),rij,dfi,dfj,drhoij,drdxi(3),drdxj(3),r,dphi,at(3)
      real(8):: x,y,z,xi(3),epotl,tmp
      logical,save:: l1st=.true.
      real(8),allocatable,save:: sqrho(:)
      real(8),save:: exrc,phic,dphic

      if( l1st ) then
        allocate(sqrho(namax+nbmax))
c-----smoothing embeded term
        exrc= exp(-ea_bt*(rc-ea_re))
c-----smoothing 2-body term
        r= rc -ea_re
        phic= 2d0*ea_b*exp(-0.5d0*ea_bt*r)
     &       -ea_c*(1d0+ea_al*r)*exp(-ea_al*r)
        dphic= -ea_bt*ea_b*exp(-0.5d0*ea_bt*r) 
     &       +ea_c*ea_al*ea_al*r*exp(-ea_al*r)
        l1st=.false.
      endif

      aa(1:3,1:natm)=0d0
      epi(1:natm)= 0d0
      epotl= 0d0
      sqrho(1:natm)= 0d0
      strs(1:3,1:3,natm+nb)= 0d0

c-----rho(i)
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
          if( rij.gt.rc ) cycle
          sqrho(i)=sqrho(i) +exp(-ea_bt*(rij-ea_re))
     &         -exrc -(rij-rc)*(-ea_bt)*exrc
        enddo
        sqrho(i)= dsqrt(sqrho(i))
      enddo

c-----copy rho of boundary atoms
      call copy_rho_ba(tcom,namax,natm,nb,nbmax,lsb,lsrc,myparity,nn,sv
     &     ,mpi_md_world,sqrho)

c-----dE/dr_i
      do i=1,natm
        xi(1:3)= ra(1:3,i)
        dfi= -0.5d0*ea_a/sqrho(i)
        do k=1,lspr(0,i)
          j=lspr(k,i)
          if(j.eq.0) exit
          if(j.le.i) cycle
          x= ra(1,j) -xi(1)
          y= ra(2,j) -xi(2)
          z= ra(3,j) -xi(3)
          xij(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
          rij=sqrt(xij(1)**2+ xij(2)**2 +xij(3)**2)
          if( rij.gt.rc ) cycle
          drdxi(1:3)= -xij(1:3)/rij
          r= rij -ea_re
c---------2-body term
          tmp= 0.5d0 *( 2d0*ea_b*exp(-0.5d0*ea_bt*r)
     &         -ea_c*(1d0+ea_al*r)*exp(-ea_al*r)
     &         -phic -(rij-rc)*dphic )
          epi(i)= epi(i) +tmp
          epi(j)= epi(j) +tmp
          if(j.le.natm) then
            epotl=epotl +tmp +tmp
          else
            epotl=epotl +tmp
          endif
          dphi= -ea_bt*ea_b*exp(-0.5d0*ea_bt*r) 
     &         +ea_c*ea_al*ea_al*r*exp(-ea_al*r)
     &         -dphic
          aa(1:3,i)=aa(1:3,i) -dphi*drdxi(1:3)
          aa(1:3,j)=aa(1:3,j) +dphi*drdxi(1:3)
c---------embedded term
          drhoij= -ea_bt*exp(-ea_bt*r) +ea_bt*exrc
          dfj= -0.5d0*ea_a/sqrho(j)
          aa(1:3,i)=aa(1:3,i) -(dfi+dfj)*drhoij*drdxi(1:3)
          aa(1:3,j)=aa(1:3,j) +(dfi+dfj)*drhoij*drdxi(1:3)
        enddo
        epi(i)=epi(i) -ea_a*sqrho(i)
        epotl=epotl -ea_a*sqrho(i)
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

c      deallocate(sqrho)
      end subroutine get_force
c=======================================================================
      subroutine copy_rho_ba(tcom,namax,natm,nb,nbmax,lsb
     &     ,lsrc,myparity,nn,sv,mpi_md_world,rho)
c-----------------------------------------------------------------------
c     Exchanges boundary-atom data among neighbor nodes
c-----------------------------------------------------------------------
      implicit none
      include "mpif.h"
      integer:: status(MPI_STATUS_SIZE)
c-----in
      integer,intent(in):: namax,natm,nb,nbmax,mpi_md_world
      integer,intent(in):: lsb(0:nbmax,6),lsrc(6),myparity(3),nn(6)
      real(8),intent(in):: sv(3,6)
c-----out
      real(8),intent(inout):: rho(natm+nb),tcom

c-----locals
      integer:: i,j,k,l,m,n,kd,kdd,ku,inode,nsd,nsd3,nrc,nrc3,nbnew,ierr
      real(8):: tcom1,tcom2
      logical,save:: l1st=.true.
      real(8),allocatable,save:: dbuf(:),dbufr(:)

      if( l1st ) then
        allocate(dbuf(nbmax),dbufr(nbmax))
        l1st=.false.
      endif

      nbnew= 0

c-----loop over z, y, & x directions
      do kd=1,3
        tcom1= mpi_wtime()
        do kdd=-1,0
          ku= 2*kd +kdd
          inode= nn(ku)
c---------num. of to-be-sent particles
          nsd= lsb(0,ku)
c---------num. of to-be-recieved particles
          nrc= lsrc(ku)

c---------exchange x
          do i=1,nsd
            j=lsb(i,ku)
            dbuf(i)= rho(j)
          enddo
          call mespasd(inode,myparity(kd),dbuf,dbufr,nsd,nrc,21
     &         ,mpi_md_world)
          do i=1,nrc
            rho(natm+nbnew+i)= dbufr(i)
          enddo

c---------mpi barrier
          call mpi_barrier(mpi_md_world,ierr)
c---------accumulate num. of boundary particles
c          write(6,'(a,2i8)') "nbnew,nrc=",nbnew,nrc
          nbnew=nbnew +nrc
        enddo
        tcom2= mpi_wtime()
        tcom= tcom +tcom2-tcom1
      enddo

      if(nbnew.ne.nb) then
        write(6,'(a,2i8)') "nbnew,(natm+nb)=",nbnew,natm+nb
        stop "error: nbnew.ne.(natm+nb)!!"
      endif
      
      end subroutine copy_rho_ba
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make pmd"
c     End:
