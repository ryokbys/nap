module Morse
!-----------------------------------------------------------------------
!                     Last modified: <2017-03-30 17:51:07 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
!  Parallel implementation of Morse pontential.
!    - For BVS, see Adams & Rao, Phys. Status Solidi A 208, No.8 (2011)
!    - Currently no cutoff tail treatment is done. (170310)
!-----------------------------------------------------------------------
  implicit none
  character(len=128),parameter:: paramsfname = 'in.params.Morse'
  character(len=128),parameter:: configfname = 'in.config.Morse'

  integer,parameter:: ioprms = 20
  integer,parameter:: iocnfg = 21

!.....Number of species
  integer:: nsp
!.....Morse parameters
  real(8),allocatable:: alp(:,:),d0(:,:),rmin(:,:)
  logical,allocatable:: interact(:,:)
  
contains
  subroutine force_Morse(namax,natm,tag,ra,nnmax,aa,strs,h,hi,tcom &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid,epi,epot,nismax,acon,lstrs,iprint)
    implicit none
    include "mpif.h"
    include "./params_unit.h"
!!$    include "params_BVS_Morse.h"
    integer,intent(in):: namax,natm,nnmax,nismax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),lspr(0:nnmax,namax),nex(3)
    integer,intent(in):: mpi_md_world,myid
    real(8),intent(in):: ra(3,namax),h(3,3),hi(3,3),rc &
         ,acon(nismax),tag(namax),sv(3,6)
    real(8),intent(inout):: tcom
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical:: lstrs

    integer:: i,j,k,l,m,n,ierr,is,js,ixyz,jxyz
    real(8):: xi(3),xj(3),xij(3),rij(3),dij,diji(3),dedr,epott &
         ,dxdi(3),dxdj(3),x,y,z,epotl,at(3),tmp,texp,d0ij,alpij,rminij
    real(8),allocatable,save:: strsl(:,:,:)

    logical,save:: l1st=.true.

    if( l1st ) then
      call initialize(natm,tag,mpi_md_world)
      call read_params(myid,mpi_md_world,iprint)
      if( .not. allocated(strsl) ) then
        allocate(strsl(3,3,namax))
      endif
      l1st=.false.
    endif

    epotl= 0d0
    strsl(1:3,1:3,1:namax) = 0d0

!.....Loop over resident atoms
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is=int(tag(i))
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
        if(j.le.i) cycle
        js= int(tag(j))
!.....Check if these two species interact
        if( .not. interact(is,js) ) cycle
        xj(1:3)= ra(1:3,j)
        xij(1:3)= xj(1:3)-xi(1:3)
        rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij= sqrt(rij(1)**2 +rij(2)**2 +rij(3)**2)
        if( dij.gt.rc ) cycle
        diji= 1d0/dij
        dxdi(1:3)= -rij(1:3)*diji
        dxdj(1:3)=  rij(1:3)*diji
        d0ij = d0(is,js)
        alpij= alp(is,js)
        rminij=rmin(is,js)
        texp = exp(alpij*(rminij-dij))
!.....potential
        tmp= 0.5d0 * d0ij*((texp-1d0)**2 -1d0)
        if( j.le.natm ) then
          epi(i)= epi(i) +tmp
          epi(j)= epi(j) +tmp
          epotl = epotl +tmp +tmp
        else
          epi(i)= epi(i) +tmp
          epotl = epotl +tmp
        endif
!.....force
        dedr= 2d0 *alpij *d0ij *texp *(1d0 -texp)
        aa(1:3,i)= aa(1:3,i) -dxdi(1:3)*dedr
        aa(1:3,j)= aa(1:3,j) -dxdj(1:3)*dedr
!.....stress
        if( lstrs ) then
          do ixyz=1,3
            do jxyz=1,3
              strsl(jxyz,ixyz,i)= strsl(jxyz,ixyz,i) &
                   -0.5d0 *dedr*rij(ixyz)*(-dxdi(jxyz))
              strsl(jxyz,ixyz,j)= strsl(jxyz,ixyz,j) &
                   -0.5d0 *dedr*rij(ixyz)*(-dxdi(jxyz))
            enddo
          enddo
        endif
      enddo
    enddo

    if( lstrs ) then
      call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
           ,nn,mpi_md_world,strsl,9)
      strs(1:3,1:3,1:natm)= strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)
    endif

    
!-----gather epot
    epott= 0d0
    call mpi_allreduce(epotl,epott,1,MPI_DOUBLE_PRECISION &
         ,MPI_SUM,mpi_md_world,ierr)
    epot= epot +epott
 
  end subroutine force_Morse
!=======================================================================
  subroutine initialize(natm,tag,mpi_md_world)
!
!  Allocate and initialize parameters to be used.
!
    include 'mpif.h'
    integer,intent(in):: natm,mpi_md_world
    real(8),intent(in):: tag(natm)
    integer:: i,nspl,ierr

!.....Get umber of species
    nspl = 0
    do i=1,natm
      nspl = max(int(tag(i)),nspl)
    enddo
    call mpi_allreduce(nspl,nsp,1,mpi_integer,mpi_max &
         ,mpi_md_world,ierr)

!.....Allocate parameter arrays
    allocate(alp(nsp,nsp),d0(nsp,nsp),rmin(nsp,nsp),interact(nsp,nsp))
    
  end subroutine initialize
!=======================================================================
  subroutine read_params(myid_md,mpi_md_world,iprint)
!
!  Read pair parameters from file
!
    include 'mpif.h'
    integer,intent(in):: myid_md,mpi_md_world,iprint
    integer:: i,j,isp,jsp,id,ierr
    character(len=128):: cline
    real(8):: d,r,a

    if( myid_md.eq.0 ) then
      open(ioprms,file=trim(paramsfname),status='old')
      interact(1:nsp,1:nsp) = .false.
      d0(1:nsp,1:nsp)= 0d0
      rmin(1:nsp,1:nsp)= 0d0
      alp(1:nsp,1:nsp)= 0d0
      do while(.true.)
        read(ioprms,*,end=10) cline
        if( cline(1:1).eq.'#' .or. cline(1:1).eq.'!' ) cycle
        backspace(ioprms)
        read(ioprms,*) isp,jsp,d,r,a
        if( isp.gt.nsp .or. jsp.gt.nsp ) then
          write(6,*) ' Warning @read_params: since isp/jsp is greater than nsp,'&
               //' skip reading the line.'
          cycle
        endif
        d0(isp,jsp) = d
        rmin(isp,jsp) = r
        alp(isp,jsp) = a
        interact(isp,jsp) = .true.
!.....Symmetrize parameters
        d0(jsp,isp) = d0(isp,jsp)
        rmin(jsp,isp)= rmin(isp,jsp)
        alp(jsp,isp)= alp(isp,jsp)
        interact(jsp,isp)= interact(isp,jsp)
      enddo

10    close(ioprms)
      if( iprint.ne.0 ) then
        write(6,'(a)') ' finish reading '//trim(paramsfname)
        write(6,*) ''
      endif
    endif

    call mpi_bcast(d0,nsp*nsp,mpi_double_precision,0,mpi_md_world,ierr)
    call mpi_bcast(rmin,nsp*nsp,mpi_double_precision,0,mpi_md_world,ierr)
    call mpi_bcast(alp,nsp*nsp,mpi_double_precision,0,mpi_md_world,ierr)
    call mpi_bcast(interact,nsp*nsp,mpi_logical,0,mpi_md_world,ierr)
    
  end subroutine read_params
end module Morse
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
