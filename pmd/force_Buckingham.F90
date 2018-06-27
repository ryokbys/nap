module Buckingham
!-----------------------------------------------------------------------
!                     Last modified: <2018-06-27 17:42:31 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
!  Parallel implementation of Buckingham calculation
!    - only force on i is considered, no need to send back
!-----------------------------------------------------------------------
  implicit none
  character(len=128):: paramsdir = '.'
  character(len=128),parameter:: cprmfname = 'in.params.Buckingham'

  integer,parameter:: ioprms = 20

!.....Max number of species available in the potential
  integer,parameter:: msp = 9
  integer:: nspcs
  real(8),allocatable:: buck_a(:,:), buck_rho(:,:), buck_c(:,:)
  logical,allocatable:: interact(:,:)

  logical:: lprmset_Buckingham
  
contains
  subroutine force_Buckingham(namax,natm,tag,ra,nnmax,aa,strs,h,hi,tcom &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid,epi,epot,nismax,lstrs,iprint,l1st)
    implicit none
    include "mpif.h"
    include "./params_unit.h"
    integer,intent(in):: namax,natm,nnmax,nismax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),lspr(0:nnmax,namax),nex(3)
    integer,intent(in):: mpi_md_world,myid
    real(8),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),rc &
         ,tag(namax),sv(3,6)
    real(8),intent(inout):: tcom
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: l1st 
    logical:: lstrs

    integer:: i,j,k,l,m,n,ierr,is,js,ixyz,jxyz
    real(8):: xi(3),xij(3),rij,riji,dvdr,rij2,rij6 &
         ,dxdi(3),dxdj(3),x,y,z,epotl,epott,at(3),tmp &
         ,a,rho,c
    real(8),allocatable,save:: strsl(:,:,:)

    integer,external:: itotOf

    if( l1st ) then
      if( allocated(strsl) ) deallocate(strsl)
      allocate(strsl(3,3,namax))
    endif

    if( size(strsl).lt.3*3*namax ) then
      deallocate(strsl)
      allocate(strsl(3,3,namax))
    endif

    epotl= 0d0
    strsl(1:3,1:3,1:namax) = 0d0

!-----loop over resident atoms
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is = int(tag(i))
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
        if(j.le.i) cycle
        js = int(tag(j))
        if( .not.interact(is,js) ) cycle
        x= ra(1,j) -xi(1)
        y= ra(2,j) -xi(2)
        z= ra(3,j) -xi(3)
        xij(1:3)= h(1:3,1,0)*x +h(1:3,2,0)*y +h(1:3,3,0)*z
        rij= sqrt(xij(1)**2+ xij(2)**2 +xij(3)**2)
        riji= 1d0/rij
        dxdi(1:3)= -xij(1:3)*riji
        dxdj(1:3)=  xij(1:3)*riji
        a = buck_a(is,js)
        rho = buck_rho(is,js)
        c = buck_c(is,js)
        rij2 = rij*rij
        rij6 = rij2*rij2*rij2
        dvdr= -a/rho *exp(-rij/rho) +6d0 *c/rij6 /rij
!---------force
        aa(1:3,i)=aa(1:3,i) -dxdi(1:3)*dvdr
        aa(1:3,j)=aa(1:3,j) -dxdj(1:3)*dvdr
!---------potential
        tmp= 0.5d0 *(a*exp(-rij/rho) -c/rij6)
        if(j.le.natm) then
          epi(i)=epi(i) +tmp
          epi(j)=epi(j) +tmp
          epotl= epotl +tmp +tmp
        else
          epi(i)=epi(i) +tmp
          epotl= epotl +tmp
        endif
!.....Stress
        if( lstrs ) then
          do ixyz=1,3
            do jxyz=1,3
              strsl(jxyz,ixyz,i)=strsl(jxyz,ixyz,i) &
                   -0.5d0*dvdr*xij(ixyz)*(-dxdi(jxyz))
              strsl(jxyz,ixyz,j)=strsl(jxyz,ixyz,j) &
                   -0.5d0*dvdr*xij(ixyz)*(-dxdi(jxyz))
            enddo
          enddo
        endif
      enddo
    enddo

    if( lstrs ) then
!!$      call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
!!$           ,nn,mpi_md_world,strsl,9)
      strs(1:3,1:3,1:natm)= strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)
    endif

!-----gather epot
    epott = 0d0
    call mpi_allreduce(epotl,epott,1,mpi_real8 &
         ,mpi_sum,mpi_md_world,ierr)
    epot= epot +epott
    
  end subroutine force_Buckingham
!=======================================================================
  subroutine init_Buckingham()
    implicit none

    if( .not.allocated(buck_a) ) then
      allocate(buck_a(msp,msp),buck_rho(msp,msp),buck_c(msp,msp)&
           ,interact(msp,msp))
    endif
    
  end subroutine init_Buckingham
!=======================================================================
  subroutine read_params_Buckingham(myid_md,mpi_md_world,iprint)
    implicit none
    include 'mpif.h'
    integer,intent(in):: myid_md,mpi_md_world,iprint

    integer:: isp,jsp,ierr
    real(8):: a,rho,c
    character(len=128):: cline,fname

    if( myid_md.eq.0 ) then
      fname = trim(paramsdir)//'/'//trim(cprmfname)
      open(ioprms,file=trim(fname),status='old')
      interact(1:msp,1:msp) = .false.
      buck_a(1:msp,1:msp) = 0d0
      buck_rho(1:msp,1:msp) = 0d0
      buck_c(1:msp,1:msp) = 0d0
      write(6,'(/,a)') ' Buckingham parameters:'
      do while(.true.)
        read(ioprms,*,end=10) cline
        if( cline(1:1).eq.'#' .or. cline(1:1).eq.'!' ) cycle
        backspace(ioprms)
        read(ioprms,*) isp,jsp,a,rho,c
!!$        print *,'isp,jsp,a,rho,c=',isp,jsp,a,rho,c
        if( isp.gt.msp .or. jsp.gt.msp ) then
          write(6,*) ' Warning @read_params: since isp/jsp is greater than msp,'&
               //' skip reading the line.'
          cycle
        endif
        buck_a(isp,jsp) = a
        buck_rho(isp,jsp) = rho
        buck_c(isp,jsp) = c
        interact(isp,jsp) = .true.
        write(6,'(a,2i3,3f10.3)') '   is,js,A,rho,C = ',isp,jsp,a,rho,c
!.....Symmetrize parameters
        buck_a(jsp,isp) = buck_a(isp,jsp)
        buck_rho(jsp,isp)= buck_rho(isp,jsp)
        buck_c(jsp,isp)= buck_c(isp,jsp)
        interact(jsp,isp)= interact(isp,jsp)
      enddo
10    close(ioprms)
    endif

    call mpi_bcast(buck_a,msp*msp,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(buck_rho,msp*msp,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(buck_c,msp*msp,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(interact,msp*msp,mpi_logical,0,mpi_md_world,ierr)
    
  end subroutine read_params_Buckingham
end module Buckingham
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
