module LJ
  implicit none
  include "mpif.h"
  include "./params_unit.h"
  integer,parameter:: ioprms = 50
  character(len=128):: paramsdir = '.'
  character(len=128),parameter:: paramsfname = 'in.params.LJ_repul'
  character(len=128),parameter:: paramsfname_repul = 'in.params.LJ_repul'

!-----LJ parameters for Argon
  real(8),parameter:: epslj = 120d0 *fkb
  real(8),parameter:: sgmlj = 3.41d-10 /ang
  real(8),parameter:: am_ar = 39.948d0
  real(8),parameter:: alcar = 2d0**(1d0/6)*sgmlj *1.41421356d0*0.996d0

!.....Max num of species
  integer,parameter:: msp = 9
  integer:: nsp

  real(8),allocatable:: strsl(:,:,:)
  real(8),allocatable:: eps(:,:), sgm(:,:)
  logical:: interact(msp,msp)
  real(8):: repul(msp,msp),  rclj(msp,msp)
  
contains
  subroutine force_LJ(namax,natm,tag,ra,nnmax,aa,strs,h,hi,tcom &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid,epi,epot,nismax,lstrs,iprint)
!-----------------------------------------------------------------------
!  Parallel implementation of LJ force calculation
!    - only force on i is considered, no need to send back
!-----------------------------------------------------------------------
    implicit none
!!$    include "./params_unit.h"
!!$    include "params_LJ.h"
    integer,intent(in):: namax,natm,nnmax,nismax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),lspr(0:nnmax,namax),nex(3)
    integer,intent(in):: mpi_md_world,myid
    real(8),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),rc &
         ,tag(namax),sv(3,6)
    real(8),intent(inout):: tcom
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical:: lstrs

    integer:: i,j,k,l,m,n,ierr,is,ixyz,jxyz
    real(8):: xi(3),xij(3),rij,riji,dvdr &
         ,dxdi(3),dxdj(3),x,y,z,epotl,epott,at(3),tmp

    logical,save:: l1st=.true.
    real(8),save:: vrc,dvdrc

    if( l1st ) then
!!$!.....assuming fixed atomic volume
!!$      avol= alcar**3 /4
!!$      if(myid.eq.0) write(6,'(a,es12.4)') ' avol =',avol
!.....prefactors
      vrc= 4d0 *epslj *((sgmlj/rc)**12 -(sgmlj/rc)**6)
      dvdrc=-24.d0 *epslj *( 2.d0*sgmlj**12/(rc**13) &
           -sgmlj**6/(rc**7) )
!!$!.....assuming fixed (constant) atomic volume
!!$      avol= (2d0**(1d0/6) *sgmlj)**2 *sqrt(3d0) /2
      l1st=.false.
    endif

    epotl= 0d0
!    print *, ' force_LJ_Ar 1'

!-----loop over resident atoms
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
        dvdr=(-24.d0*epslj)*(2.d0*(sgmlj*riji)**12*riji &
             -(sgmlj*riji)**6*riji) &
             -dvdrc
!---------force
        aa(1:3,i)=aa(1:3,i) -dxdi(1:3)*dvdr
        aa(1:3,j)=aa(1:3,j) -dxdj(1:3)*dvdr
!---------potential
        tmp= 0.5d0*( 4d0*epslj*((sgmlj*riji)**12 &
             -(sgmlj*riji)**6) -vrc -dvdrc*(rij-rc) )
        if(j.le.natm) then
          epi(i)=epi(i) +tmp
          epi(j)=epi(j) +tmp
          epotl= epotl +tmp +tmp
        else
          epi(i)=epi(i) +tmp
          epotl= epotl +tmp
        endif
!---------stress
        if( lstrs ) then
          do ixyz=1,3
            do jxyz=1,3
              strs(jxyz,ixyz,i)=strs(jxyz,ixyz,i) &
                   -0.5d0*dvdr*xij(ixyz)*(-dxdi(jxyz))
              strs(jxyz,ixyz,j)=strs(jxyz,ixyz,j) &
                   -0.5d0*dvdr*xij(ixyz)*(-dxdi(jxyz))
            enddo
          enddo
        endif
      enddo
    enddo

    if( lstrs ) then
      call copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
           ,nn,mpi_md_world,strs,9)
    endif

!-----gather epot
    if( myid.ge.0 ) then
      call mpi_allreduce(epotl,epott,1,MPI_DOUBLE_PRECISION &
           ,MPI_SUM,mpi_md_world,ierr)
      epot= epot +epott
    else
      epot= epot +epotl
    endif
  end subroutine force_LJ
!=======================================================================
  subroutine force_LJ_repul(namax,natm,tag,ra,nnmax,aa,strs,h,hi,tcom &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid,epi,epot,nismax,lstrs,iprint,l1st)
!
!  LJ potential of only repulsive term
!
    implicit none
!!$    include "./params_unit.h"
!!$    include "params_LJ.h"
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

    integer:: i,j,k,l,m,n,ierr,is,js,ixyz,jxyz,jj
    real(8):: xi(3),xij(3),rij(3),dij,dij2,diji,diji6,diji12,dvdr &
         ,dxdi(3),dxdj(3),x,y,z,epotl,epott,at(3),tmp,rcij,repij
    real(8),external:: fcut1,dfcut1

    if( l1st ) then
!!$      call init_Morse(natm,tag,mpi_md_world)
!!$      call read_params_Morse(myid,mpi_md_world,iprint)
      if( allocated(strsl) ) deallocate(strsl)
      allocate(strsl(3,3,namax))
    endif

    if( size(strsl).lt.3*3*namax ) then
      deallocate(strsl)
      allocate(strsl(3,3,namax))
    endif

    epotl = 0d0
    strsl(1:3,1:3,1:namax) = 0d0

    do i=1,natm
      xi(1:3) = ra(1:3,i)
      is = int(tag(i))
      do jj=1,lspr(0,i)
        j = lspr(jj,i)
        if( j.eq.0 ) exit
        if( j.le.i ) cycle
        js = int(tag(j))
        rcij = rclj(is,js)
        rij(1:3) = ra(1:3,j) -xi(1:3)
        xij(1:3) = h(1:3,1,0)*rij(1) +h(1:3,2,0)*rij(2) +h(1:3,3,0)*rij(3)
        dij2 = xij(1)*xij(1) +xij(2)*xij(2) +xij(3)*xij(3)
        if( dij2.gt.rcij*rcij ) cycle
        repij = repul(is,js)
        dij = sqrt(dij2)
        diji = 1d0/dij
        dxdi(1:3) = -xij(1:3)*diji
        dxdj(1:3) =  xij(1:3)*diji
!.....Potential
        diji6 = 1d0/(dij2*dij2*dij2)
        diji12 = diji6*diji6
        tmp = 0.5d0 *repij *diji12 *fcut1(dij,0d0,rcij)
        if( j.le.natm ) then
          epi(i) = epi(i) +tmp
          epi(j) = epi(j) +tmp
          epotl = epotl +tmp +tmp
        else
          epi(i) = epi(i)
          epotl = epotl +tmp
        endif
!.....Force
        dvdr = -12d0*repij*diji12*diji *fcut1(dij,0d0,rcij) &
             -repij*diji12 *dfcut1(dij,0d0,rcij)
        aa(1:3,i) = aa(1:3,i) -dxdi(1:3)*dvdr
        aa(1:3,j) = aa(1:3,j) -dxdj(1:3)*dvdr
!.....Stress
        if( lstrs ) then
          do ixyz=1,3
            do jxyz=1,3
              strsl(jxyz,ixyz,i) = strsl(jxyz,ixyz,i) &
                   -0.5d0*dvdr*xij(ixyz)*(-dxdi(jxyz))
              strsl(jxyz,ixyz,j) = strsl(jxyz,ixyz,j) &
                   -0.5d0*dvdr*xij(ixyz)*(-dxdi(jxyz))
            enddo
          enddo
        endif
      enddo
    enddo

    if( lstrs ) then
      strs(1:3,1:3,1:natm)= strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)
    endif
    
!-----gather epot
    epott= 0d0
    call mpi_allreduce(epotl,epott,1,mpi_real8,mpi_sum,mpi_md_world,ierr)
    epot= epot +epott
    if( iprint.gt.2 ) print *,'LJ_repul epot = ',epott

  end subroutine force_LJ_repul
!=======================================================================
  subroutine read_params_LJ_repul(myid,mpi_world,iprint)
!
!  Read input for LJ_repul force
!
    integer,intent(in):: myid,mpi_world,iprint
    character(len=128):: cline,cfname
    integer:: isp,jsp,nd,ierr
    logical:: lexist
    real(8):: repij,rcij
    integer,external:: num_data

    if( myid.eq.0 ) then
!.....Initialization
      interact(1:msp,1:msp) = .false.
!.....Check whether the file exists      
      cfname = trim(paramsdir)//'/'//trim(paramsfname_repul)
      inquire(file=cfname,exist=lexist)
      if( .not. lexist ) then
        if( iprint.gt.0 ) then
          write(6,'(a)') ' WARNING: in.params.LJ_repul does not exist !!!.'
          write(6,'(a)') '           Default parameters will be used.'
        endif
        interact(1,1) = .true.
        repul(1,1) = 4d0*epslj*sgmlj**12
        return
      endif
!.....Read file if exits
      if( iprint.ne.0 ) write(6,'(/,a)') ' LJ_repul parameters read from file:'
      open(ioprms,file=cfname,status='old')
      do while(.true.)
        read(ioprms,'(a)',end=10) cline
        nd = num_data(cline,' ')
        if( nd.eq.0 ) cycle
        if( cline(1:1).eq.'!' .or. cline(1:1).eq.'#' ) cycle
        if( nd.eq. 4 ) then
          backspace(ioprms)
          read(ioprms,*) isp,jsp, repij, rcij
          if( iprint.ne.0 ) write(6,'(a,2i4,f8.4,f6.2)') &
               '   isp,jsp,repij,rcij = ',isp,jsp,repij,rcij
          interact(isp,jsp) = .true.
          repul(isp,jsp) = repij
          rclj(isp,jsp) = rcij
!.....Symmetrize
          interact(jsp,isp) = interact(isp,jsp)
          repul(jsp,isp) = repul(isp,jsp)
          rclj(jsp,isp) = rclj(isp,jsp)
        else
          if( iprint.ne.0 ) then
            write(6,*) 'WARNING@read_params_LJ_repul: number of entry wrong, ' &
                 //'so skip the line.'
          endif
          cycle
        endif
      enddo
10    close(ioprms)
      
    endif

    call mpi_bcast(interact,msp*msp,mpi_logical,0,mpi_world,ierr)
    call mpi_bcast(repul,msp*msp,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(rclj,msp*msp,mpi_real8,0,mpi_world,ierr)
    return
  end subroutine read_params_LJ_repul
!=======================================================================
  subroutine set_paramsdir_LJ(dname)
!
!  Accessor routine to set paramsdir.
!
    implicit none
    character(len=*),intent(in):: dname

    paramsdir = trim(dname)
    return
  end subroutine set_paramsdir_LJ
end module LJ
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
