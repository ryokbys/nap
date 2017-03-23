      subroutine get_force(namax,natm,tag,ra,nnmax,aa,strs,chg,chi
     &     ,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr
     &     ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs
     &     ,cforce,ifchg,iprint)
c-----------------------------------------------------------------------
c  Wrapper routine for force calculations.
c  Each force calculation routine is called from this subroutine.
c  Any new force routine should be implemented in this subroutine.
c-----------------------------------------------------------------------
      use RK_FeH,only:force_RK_FeH
      use Ramas_FeH,only:force_Ramas_FeH,force_Ackland_Fe
      use RK_WHe,only:force_RK_WHe
      use Ito3_WHe,only:force_Ito3_WHe
      use LJ_Ar,only:force_LJ_Ar
      use SW_Si,only:force_SW_Si
      use EDIP_Si,only:force_EDIP_Si
      use Brenner,only:force_brenner,force_brenner_vdW
      use Lu_WHe,only:force_Lu_WHe
      use Branicio_AlN,only:force_Branicio_AlN
      use Mishin_Al,only:force_Mishin_Al
      use AFS_W,only:force_AFS_W
      use SC_Fe,only:force_SC_Fe
      use SM_Al,only:force_SM_Al
      use linreg,only:force_linreg
      use NN,only:force_NN
      implicit none
      integer,intent(in):: namax,natm,nnmax,nismax,iprint,ifchg
      integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3)
     &     ,nnn(6),mpi_md_world,myid_md,nex(3)
      integer,intent(in):: lspr(0:nnmax,namax)
      real(8),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),sv(3,6)
     &     ,acon(nismax),tag(namax)
      real(8),intent(inout):: tcom,rc
      real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
     &     ,chg(namax),chi(namax)
      character(len=*),intent(in):: cforce
      logical:: lstrs
      
      integer:: ierr,is,i
      real(8):: at(3)

      epot = 0d0
      aa(1:3,1:natm+nb)=0d0
      epi(1:natm+nb)= 0d0
      strs(1:3,1:3,1:natm+nb)= 0d0
      
      select case (trim(cforce))
      case ('LJ_Ar')
        call force_LJ_Ar(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &       ,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr
     &       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
      case ('Ito3_WHe')
        call force_Ito3_WHe(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &       ,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr
     &       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
      case ('RK_WHe')
        call force_RK_WHe(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &       ,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr
     &       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
      case ('RK_FeH')
        call force_RK_FeH(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &       ,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr
     &       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
      case ('Ramas_FeH')
        call force_Ramas_FeH(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &       ,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr
     &       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
      case ('Ackland_Fe')
        call force_Ackland_Fe(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &       ,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr
     &       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
      case ('SW_Si')
        call force_SW_Si(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &       ,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr
     &       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
      case ('EDIP_Si')
        call force_EDIP_Si(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &       ,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr
     &       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
      case ('Brenner')
        call force_Brenner(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &       ,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr
     &       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
      case ('Brenner_vdW')
        call force_Brenner_vdW(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &       ,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr
     &       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
      case ('Lu_WHe')
        call force_Lu_Whe(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &       ,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr
     &       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
      case ('Branicio_AlN')
        call force_Branicio_AlN(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &       ,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr
     &       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
      case ('Mishin_Al')
        call force_Mishin_Al(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &       ,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr
     &       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
      case ('AFS_W')
        call force_AFS_W(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &       ,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr
     &       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
      case ('SC_Fe')
        call force_SC_Fe(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &       ,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr
     &       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
      case ('SM_Al')
        call force_SM_Al(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &       ,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr
     &       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
      case ('linreg')
        call force_linreg(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &       ,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr
     &       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
      case ('NN')
        call force_NN(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &       ,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr
     &       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
      case default
        write(6,'(a)') ' [get_force] unknown force: '
     &       //trim(cforce)//' !!'
        call mpi_finalize(ierr)
        stop
      end select

      if( ifchg .eq. 1 ) then   ! screened Coulomb

      endif

c.....convert forces from hmat-coordinates to Cartesian coordinates
      do i=1,natm
        at(1:3)= aa(1:3,i)
        aa(1:3,i)= hi(1:3,1)*at(1) +hi(1:3,2)*at(2) +hi(1:3,3)*at(3)
      enddo
c-----multiply 0.5d0*dt**2/am(i)
      do i=1,natm
        is= int(tag(i))
        aa(1:3,i)= acon(is)*aa(1:3,i)
      enddo

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
      
      logical,save:: l1st=.true.
      real(8),save,allocatable:: dbuf(:,:),dbufr(:,:)

      if( l1st ) then
        allocate(dbuf(9,nbmax),dbufr(9,nbmax))
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
      
      end subroutine copy_strs_ba
c=======================================================================
      subroutine copy_dba_fwd(tcom,namax,natm,nb,nbmax,lsb,nex
     &     ,lsrc,myparity,nn,sv,mpi_md_world,x,ndim)
c-----------------------------------------------------------------------
c     Exchanges boundary-atom data among neighbor nodes
c-----------------------------------------------------------------------
      implicit none
      include "mpif.h"
      integer:: status(MPI_STATUS_SIZE)
c-----in
      integer,intent(in):: namax,natm,nb,nbmax,mpi_md_world,ndim
      integer,intent(in):: lsb(0:nbmax,6),lsrc(6),myparity(3),nn(6)
     &     ,nex(3)
      real(8),intent(in):: sv(3,6)
c-----out
      real(8),intent(inout):: x(ndim,namax),tcom

c-----locals
      integer:: i,j,k,l,m,n,kd,kdd,ku,inode,nsd,nsd3,nrc,nrc3,nbnew,ierr
      real(8):: tcom1,tcom2
      logical,save:: l1st=.true.
      real(8),allocatable,save:: dbuf(:,:),dbufr(:,:)
      integer,save:: mdim

      if( l1st ) then
        mdim= ndim
        allocate(dbuf(mdim,nbmax),dbufr(mdim,nbmax))
        l1st=.false.
      endif

      if( ndim.gt.mdim ) then
        deallocate(dbuf,dbufr)
        mdim= ndim
        allocate(dbuf(mdim,nbmax),dbufr(mdim,nbmax))
      endif

      nbnew= 0

c-----loop over z, y, & x directions
      do kd=1,3

        if( nex(kd).gt.1 ) then

          do kdd=-1,0
            ku= 2*kd +kdd
            nrc= lsb(0,ku)
            do i=1,nrc
              j= lsb(i,ku)
              x(1:ndim,natm+nbnew+i)= x(1:ndim,j)
            enddo
            nbnew= nbnew +nrc
          enddo
        else

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
              dbuf(1:ndim,i)= x(1:ndim,j)
            enddo
            call mespasd(inode,myparity(kd),dbuf,dbufr,ndim*nsd
     &           ,ndim*nrc,21,mpi_md_world)
            do i=1,nrc
              x(1:ndim,natm+nbnew+i)= dbufr(1:ndim,i)
            enddo
c---------mpi barrier
            call mpi_barrier(mpi_md_world,ierr)
c---------accumulate num. of boundary particles
c          write(6,'(a,2i8)') "nbnew,nrc=",nbnew,nrc
            nbnew=nbnew +nrc
          enddo
          tcom2= mpi_wtime()
          tcom= tcom +tcom2-tcom1
        endif
      enddo

      if(nbnew.ne.nb) then
        write(6,'(a,2i8)') "nbnew,(natm+nb)=",nbnew,natm+nb
        stop "error: nbnew.ne.(natm+nb)!!"
      endif
      
      end subroutine copy_dba_fwd
c=======================================================================
      subroutine copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex
     &     ,lsrc,myparity,nn,mpi_md_world,x,ndim)
c-----------------------------------------------------------------------
c     Send-back & receive reaction on cached-atoms
c-----------------------------------------------------------------------
      implicit none
      include "mpif.h"
      integer,intent(in):: namax,natm,nbmax,nb,mpi_md_world,ndim
      integer,intent(in):: lsb(0:nbmax,6),lsrc(6),myparity(3),nn(6)
     &     ,nex(3)
      real(8),intent(inout):: x(ndim,namax),tcom

      integer:: status(MPI_STATUS_SIZE)
      integer:: i,j,k,l,m,n,kd,kdd,ku,kuc,ibkwd,nsd,nsd3,nrc,nrc3,nsdbk
     &     ,ierr,natmx
      real(8):: tcom1,tcom2
      real(8),save,allocatable:: dbuf(:,:),dbufr(:,:)
      logical,save:: l1st=.true.
      integer,save:: mdim

      if( l1st ) then
        mdim= ndim
        allocate(dbuf(mdim,nbmax),dbufr(mdim,nbmax))
        l1st=.false.
      endif

      if( ndim.gt.mdim ) then
        deallocate(dbuf,dbufr)
        mdim= ndim
        allocate(dbuf(mdim,nbmax),dbufr(mdim,nbmax))
      endif

c-----natmx
      natmx= natm +nb

c-----num. of sent-back reactions
      nsdbk= 0

c-----send-back reactions in the reverse orer, z, y, & x
      do kd=3,1,-1

        if( nex(kd).gt.1 ) then
          do kdd=0,-1,-1
            ku= 2*kd +kdd
            nsd= lsb(0,ku)
            nrc= nsd
            do i=1,nrc
              j= lsb(i,ku)
              x(1:ndim,j)= x(1:ndim,j) +x(1:ndim,natmx-nsdbk-nsd+i)
            enddo
c---------accumulate num. of already sent-back-particles
            nsdbk=nsdbk +nsd
          enddo
        else
c-------To calculate the communication time
          tcom1=MPI_WTIME()

c-------higher & lower directions
          do kdd=0,-1,-1
            ku= 2*kd +kdd
            if(mod(ku,2).eq.0) then
              kuc= ku-1
            else
              kuc= ku+1
            endif
            ibkwd= nn(kuc)
c---------num. of to-be-sent particles
            nsd= lsrc(ku)
c          nsd3= ndim*nsd
            nsd3= mdim*nsd
c---------num. of to-be-recieved particles
            nrc= lsb(0,ku)
c          nrc3= ndim*nrc
            nrc3= mdim*nrc

c---------to-be-sent-back particles
            do i=1,nsd
              dbuf(1:ndim,i)= x(1:ndim,natmx-nsdbk-nsd+i)
            enddo
            call mespasd(ibkwd,myparity(kd),dbuf,dbufr,nsd3,nrc3,500
     &           ,mpi_md_world)
            do k=1,nrc
              i=lsb(k,ku)
              x(1:ndim,i)= x(1:ndim,i) +dbufr(1:ndim,k)
            enddo

c---------mpi barrier
            call mpi_barrier(mpi_md_world,ierr)
c---------accumulate num. of already sent-back-particles
            nsdbk=nsdbk +nsd
          enddo

c-------Add the communication time to COMT
          tcom2=MPI_WTIME()
          tcom=tcom+tcom2-tcom1
        endif

      enddo

c-----check
      if(nsdbk.ne.nb) then
        write(6,'(a,2i8)') "nsdbk,nb=",nsdbk,nb
        stop "error: nsdbk.ne.nb!!"
      endif

c      deallocate(dbuf,dbufr)
      end subroutine copy_dba_bk
c=======================================================================
      subroutine reduce_dba_bk(natm,namax,tag,x,ndim)
c-----------------------------------------------------------------------
c  Send-back or reduce reaction on cached-atoms.
c  This routine works only on small MD not on parallel version.
c-----------------------------------------------------------------------
      implicit none
      integer,intent(in):: namax,natm,ndim
      real(8),intent(in):: tag(namax)
      real(8),intent(inout):: x(ndim,namax)
      integer,external:: itotOf
      integer:: ia,ja

      do ia=natm+1,namax
        if( tag(ia).lt.1d0 ) cycle
        ja= itotOf(tag(ia))
        x(1:ndim,ja)= x(1:ndim,ja) +x(1:ndim,ia)
      enddo

      end subroutine reduce_dba_bk
c=======================================================================
      subroutine distribute_dba(natm,namax,tag,x,ndim)
c-----------------------------------------------------------------------
c  Distribute some values to the cached (boundary) atoms.
c  This routine works only on small MD not on parallel version.
c-----------------------------------------------------------------------
      implicit none
      integer,intent(in):: namax,natm,ndim
      real(8),intent(in):: tag(namax)
      real(8),intent(inout):: x(ndim,namax)
      integer,external:: itotOf
      integer:: ia,ja

      do ia=natm+1,namax
        if( tag(ia).lt.1d0 ) cycle
        ja= itotOf(tag(ia))
        x(1:ndim,ia)= x(1:ndim,ja)
      enddo

      end subroutine distribute_dba
c=======================================================================
      function hvsd(x)
c
c  Heaviside's stepwise function
c
      implicit none
      real(8),intent(in):: x
      real(8):: hvsd

      hvsd= 0d0
      if( x.ge.0 ) then
        hvsd= 1d0
        return
      endif
      return 

      end function hvsd
c=======================================================================
