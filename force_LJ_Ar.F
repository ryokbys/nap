      subroutine setup(nismax,am,acon,fack,dt)
      implicit none
      include "params_au.h"
      include "params_LJ_Ar.h"
      integer,intent(in):: nismax
      real(8),intent(in):: dt
      real(8),intent(out):: am(nismax),acon(nismax),fack(nismax)

      integer:: i

c-----atomic masses
      am(1:nismax)= am_ar

c-----prefactors for normalized accelerration
      do i=1,nismax
        acon(i)= 0.5d0*dt**2/am(i)
      enddo

c-----prefactors for kinetic energy, FACK
      do i=1,nismax
        fack(i)= 0.5d0*am(i)/dt**2
      enddo

      end subroutine setup
c=======================================================================
      subroutine get_force(namax,natm,tag,ra,nnmax,aa,strs,h,hi,tcom
     &     ,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
     &     ,mpi_md_world,myid_md,epi,epot,nismax,acon)
c-----------------------------------------------------------------------
c  Parallel implementation of LJ force calculation
c    - only force on i is considered, no need to send back
c-----------------------------------------------------------------------
      implicit none
      include "mpif.h"
      include "params_au.h"
      include "params_LJ_Ar.h"
      integer,intent(in):: namax,natm,nnmax,nismax
      integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3)
     &     ,nn(6),mpi_md_world,myid_md,lspr(0:nnmax,namax)
      real(8),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),rc
     &     ,acon(nismax),tag(namax),sv(3,6)
      real(8),intent(inout):: tcom
      real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)

      integer:: i,j,k,l,m,n,ierr,is,ixyz,jxyz
      real(8):: xi(3),xij(3),rij,riji,dvdr
     &     ,dxdi(3),dxdj(3),x,y,z,epotl,at(3),tmp

      logical,save:: l1st=.true.
      real(8),save:: vrc,dvdrc,avol

      if( l1st ) then
c.....prefactors
        vrc= 4d0 *epslj *((sgmlj/rc)**12 -(sgmlj/rc)**6)
        dvdrc=-24.d0 *epslj *( 2.d0*sgmlj**12/(rc**13)
     &       -sgmlj**6/(rc**7) )
c.....assuming fixed (constant) atomic volume
        avol= (2d0**(1d0/6) *sgmlj)**2 *sqrt(3d0) /2
        l1st=.false.
      endif

      aa(1:3,1:namax)=0d0
      epi(1:namax)= 0d0
      epotl= 0d0
      strs(1:3,1:3,1:natm+nb)= 0d0

c-----loop over resident atoms
      do i=1,natm
        xi(1:3)= ra(1:3,i)
        do k=1,lspr(0,i)
          j=lspr(k,i)
          if(j.eq.0) exit
          if(j.le.i) cycle
          x= ra(1,j) -xi(1)
          y= ra(2,j) -xi(2)
          z= ra(3,j) -xi(3)
          xij(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
          rij= sqrt(xij(1)**2+ xij(2)**2 +xij(3)**2)
          riji= 1d0/rij
          dxdi(1:3)= -xij(1:3)*riji
          dxdj(1:3)=  xij(1:3)*riji
          dvdr=(-24.d0*epslj)*(2.d0*(sgmlj*riji)**12*riji
     &         -(sgmlj*riji)**6*riji)
     &         -dvdrc
c---------force
          aa(1:3,i)=aa(1:3,i) -dxdi(1:3)*dvdr
          aa(1:3,j)=aa(1:3,j) -dxdj(1:3)*dvdr
c---------potential
          tmp= 0.5d0*( 4d0*epslj*((sgmlj*riji)**12
     &         -(sgmlj*riji)**6) -vrc -dvdrc*(rij-rc) )
          if(j.le.natm) then
            epi(i)=epi(i) +tmp
            epi(j)=epi(j) +tmp
            epotl= epotl +tmp +tmp
          else
            epi(i)=epi(i) +tmp
            epotl= epotl +tmp
          endif
c---------stress
          do ixyz=1,3
            do jxyz=1,3
              strs(jxyz,ixyz,i)=strs(jxyz,ixyz,i)
     &             -0.5d0*dvdr*xij(ixyz)*(-dxdi(jxyz))
              strs(jxyz,ixyz,j)=strs(jxyz,ixyz,j)
     &             -0.5d0*dvdr*xij(ixyz)*(-dxdi(jxyz))
            enddo
          enddo
        enddo
      enddo

c-----copy strs of boundary atoms
      call copy_strs_ba(tcom,namax,natm,nb,nbmax,lsb
     &     ,lsrc,myparity,nn,sv,mpi_md_world,strs)
c-----atomic level stress in [Hartree/Bohr^3] assuming 1 Bohr thick
      do i=1,natm
        strs(1:3,1:3,i)= strs(1:3,1:3,i) /avol
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

      end subroutine get_force
c=======================================================================
      subroutine copy_strs_ba(tcom,namax,natm,nb,nbmax,lsb
     &     ,lsrc,myparity,nn,sv,mpi_md_world,strs)
c-----------------------------------------------------------------------
c  Exchanges boundary-atom data among neighbor nodes
c-----------------------------------------------------------------------
      implicit none
      include "mpif.h"
      integer:: status(MPI_STATUS_SIZE)
c-----in
      integer,intent(in):: namax,natm,nb,nbmax,mpi_md_world
      integer,intent(in):: lsb(0:nbmax,6),lsrc(6),myparity(3),nn(6)
      real(8),intent(in):: sv(3,6)
c-----out
      real(8),intent(inout):: strs(9,natm+nb),tcom

c-----locals
      integer:: i,j,k,l,m,n,kd,kdd,ku,inode,nsd,nrc,nbnew,ierr
      real(8):: tcom1,tcom2
      real(8),allocatable:: dbuf(:,:),dbufr(:,:)

      allocate(dbuf(9,nbmax),dbufr(9,nbmax))

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

c---------exchange strs
          do i=1,nsd
            j=lsb(i,ku)
            dbuf(1:9,i)= strs(1:9,j)
          enddo
          call mespasd(inode,myparity(kd),dbuf,dbufr,9*nsd,9*nrc,21
     &         ,mpi_md_world)
          do i=1,nrc
            strs(1:9,natm+nbnew+i)= dbufr(1:9,i)
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
      
      deallocate(dbuf,dbufr)
      end subroutine copy_strs_ba
c=======================================================================
      subroutine set_avol()
      use variables
      implicit none
      include './params_au.h'
      include './params_LJ_Ar.h'

      avol= alcar**3 /4
      if(myid_md.eq.0) write(6,'(a,es12.4)') ' avol =',avol
      return
      end subroutine set_avol
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make pmd"
c     End:
