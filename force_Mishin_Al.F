      subroutine get_force(namax,natm,tag,ra,nnmax,aa,strs,h,hi,tcom
     &     ,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
     &     ,mpi_md_world,myid_md,epi,epot,nismax,acon)
c-----------------------------------------------------------------------
c  Parallel implementation of force and potential energy calculation
c  of EAM for Al by Mishin et al.
c    - rho of boundary atoms are sent to the neighbor nodes
c    - only force on i is calculated, not necessary to send-back
c-----------------------------------------------------------------------
c  See Mishin et al. PRB 59(5) (1999) 3393--3407.
c-----------------------------------------------------------------------
      implicit none
      include "mpif.h"
      include "./params_au.h"
      include "params_Mishin_Al.h"
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
c.....Saved variables
      real(8),save:: exrc
      logical,save:: l1st=.true.
      real(8),allocatable,save:: rho(:)
c.....Functions
      real(8),external:: calc_rho,calc_v,calc_f
     &     ,calc_drho,calc_dv,calc_df

      if( l1st ) then
        allocate(rho(namax+nbmax))
        l1st=.false.
      endif

      aa(1:3,1:natm)=0d0
      epi(1:natm)= 0d0
      epotl= 0d0
      rho(1:natm)= 0d0
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
          if( rij.gt.rc_eam ) cycle
          rho(i)= rho(i) +calc_rho(rij*bohr2aa)
        enddo
      enddo

c-----copy rho of boundary atoms
      call copy_rho_ba(tcom,namax,natm,nb,nbmax,lsb,lsrc,myparity,nn,sv
     &     ,mpi_md_world,rho)

c-----dE/dr_i
      do i=1,natm
        xi(1:3)= ra(1:3,i)
        dfi= calc_df(rho(i)) *ev2hrt/aa2bohr
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
c.....2-body term
          tmp= 0.5d0 *calc_v(rij*bohr2aa) *ev2hrt
          epi(i)= epi(i) +tmp
          epi(j)= epi(j) +tmp
          if(j.le.natm) then
            epotl=epotl +tmp +tmp
          else
            epotl=epotl +tmp
          endif
          dphi= calc_dv(rij*bohr2aa) *ev2hrt/aa2bohr
          aa(1:3,i)=aa(1:3,i) -dphi*drdxi(1:3)
          aa(1:3,j)=aa(1:3,j) +dphi*drdxi(1:3)
c.....Embeded term
          drhoij= calc_drho(rij*bohr2aa)
          dfj= calc_df(rho(j)) *ev2hrt/aa2bohr
          aa(1:3,i)=aa(1:3,i) -(dfi+dfj)*drhoij*drdxi(1:3)
          aa(1:3,j)=aa(1:3,j) +(dfi+dfj)*drhoij*drdxi(1:3)
        enddo
        tmp= calc_f(rho(i)) *ev2hrt
        epi(i)=epi(i) +tmp
        epotl=epotl +tmp
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
      function calc_rho(r)
      implicit none
      include './params_au.h'
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
c=======================================================================
      function calc_drho(r)
      implicit none
      include './params_au.h'
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
c=======================================================================
      function calc_v(r)
      implicit none
      include './params_au.h'
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
c=======================================================================
      function calc_dv(r)
      implicit none
      include './params_au.h'
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
c=======================================================================
      function calc_f(rho)
      implicit none
      include './params_au.h'
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
c=======================================================================
      function calc_df(rho)
      implicit none
      include './params_au.h'
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
c=======================================================================
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
