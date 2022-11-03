module Buckingham
!-----------------------------------------------------------------------
!                     Last modified: <2022-11-03 14:21:14 KOBAYASHI Ryo>
!-----------------------------------------------------------------------
!  Parallel implementation of Buckingham calculation
!    - only force on i is considered, no need to send back
!-----------------------------------------------------------------------
  use pmdvars,only: nspmax
  use util,only: csp2isp
  use memory,only: accum_mem
  implicit none
  include "./const.h"
  save
  character(len=128):: paramsdir = '.'
  character(len=128),parameter:: cprmfname = 'in.params.Buckingham'

  integer,parameter:: ioprms = 20

!.....Max number of species available in the potential
  integer:: nspcs
  real(8):: buck_a(nspmax,nspmax), buck_rho(nspmax,nspmax), buck_c(nspmax,nspmax)
  logical:: interact(nspmax,nspmax)
!.....Smooth cutoff
  real(8):: vrcs(nspmax,nspmax), dvdrcs(nspmax,nspmax)

  logical:: lprmset_Buckingham
  
contains
!=======================================================================
  subroutine init_Buckingham()
    implicit none

    return
  end subroutine init_Buckingham
!=======================================================================
  subroutine force_Buckingham(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid,epi,epot,nismax,lstrs,iprint,l1st)
    use util,only: itotOf
    implicit none
    include "mpif.h"
    include "./params_unit.h"
    integer,intent(in):: namax,natm,nnmax,nismax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),lspr(0:nnmax,namax),nex(3)
    integer,intent(in):: mpi_md_world,myid
    real(8),intent(in):: ra(3,namax),h(3,3),hi(3,3),rc &
         ,tag(namax),sv(3,6)
    real(8),intent(inout):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: l1st 
    logical:: lstrs

    integer:: i,j,k,l,m,n,ierr,is,js,ixyz,jxyz
    real(8):: xi(3),xij(3),rij(3),dij,diji,dvdr,dij2,dij6 &
         ,dxdi(3),dxdj(3),epotl,epott,at(3),tmp &
         ,aij,rhoij,cij,texp,vrc,dvdrc
    real(8),allocatable,save:: strsl(:,:,:)

    real(8),save:: rcmax2

!!$    integer,external:: itotOf

    if( l1st ) then
      if( allocated(strsl) ) then
        call accum_mem('force_Buckingham',-8*size(strsl))
        deallocate(strsl)
      endif
      allocate(strsl(3,3,namax))
      call accum_mem('force_Buckingham',8*size(strsl))
      rcmax2 = rc*rc
!.....Initialize smooth cutoff
      vrcs(:,:) = 0d0
      dvdrcs(:,:) = 0d0
      do is=1,nspmax
        do js=1,nspmax
          aij = buck_a(is,js)
          rhoij = buck_rho(is,js)
          cij = buck_c(is,js)
          vrc = aij*exp(-rc/rhoij) -cij/(rc**6)
          vrcs(is,js) = vrc
          vrcs(js,is) = vrc
          dvdrc = -aij/rhoij*exp(-rc/rhoij) +6d0*cij/(rc**7)
          dvdrcs(is,js) = dvdrc
          dvdrcs(js,is) = dvdrc
        enddo
      enddo
    endif

    if( size(strsl).lt.3*3*namax ) then
      call accum_mem('force_Buckingham',-8*size(strsl))
      deallocate(strsl)
      allocate(strsl(3,3,namax))
      call accum_mem('force_Buckingham',8*size(strsl))
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
        js = int(tag(j))
        if( .not.interact(is,js) ) cycle
        xij(1:3) = ra(1:3,j) -xi(1:3)
        rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij2= rij(1)*rij(1)+ rij(2)*rij(2) +rij(3)*rij(3)
        if( dij2.ge.rcmax2 ) cycle
        dij = dsqrt(dij2)
        diji= 1d0/dij
        dxdi(1:3)= -rij(1:3)*diji
        aij = buck_a(is,js)
        rhoij = buck_rho(is,js)
        cij = buck_c(is,js)
        dij2 = dij*dij
        dij6 = dij2*dij2*dij2
        vrc = vrcs(is,js)
        dvdrc = dvdrcs(is,js)
        texp = exp(-dij/rhoij)
        dvdr= -aij/rhoij *texp +6d0 *cij/dij6 *diji -dvdrc
!---------force
        aa(1:3,i)=aa(1:3,i) -dxdi(1:3)*dvdr
!---------potential
        tmp= 0.5d0 *(aij*texp -cij/dij6 -vrc -dvdrc*(dij-rc))
        epi(i)= epi(i) +tmp
        epotl = epotl +tmp
!.....Stress
        if( lstrs ) then
          do ixyz=1,3
            do jxyz=1,3
              strsl(jxyz,ixyz,i)=strsl(jxyz,ixyz,i) &
                   -0.5d0*dvdr*rij(ixyz)*(-dxdi(jxyz))
            enddo
          enddo
        endif
      enddo
    enddo

    if( lstrs ) then
      strs(1:3,1:3,1:natm)= strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)
    endif

!-----gather epot
    epott = 0d0
    call mpi_allreduce(epotl,epott,1,mpi_real8,mpi_sum,mpi_md_world,ierr)
    epot= epot +epott
    if( myid.eq.0 .and. iprint.ge.ipl_info ) &
         write(6,'(a,es15.7)') ' epot Buckingham = ',epott
    return
  end subroutine force_Buckingham
!=======================================================================
  subroutine read_params_Buckingham(myid_md,mpi_md_world,iprint)
!
!  Read potential parameters from in.params.Buckingham.
!  Potential function shape is written as,
!  
!      phi(r) = A*exp( -r/rho) -C/r^6,    for r < rc
!  
!  The file format is like below:
!-----------------------------------------------------------------------
!  #  Buckingham potential parameters for garnet LLZ
!  #  cspi, cspj, Aij(eV), rhoij(Ang), cij(eV/Ang^6)
!     Li    O      876.86      0.2433      0.0
!     La    O    14509.63      0.2438     30.83
!     Zr    O     1366.09      0.3181      0.0
!     O     O     4869.99      0.2402     27.22
!-----------------------------------------------------------------------
    implicit none
    include 'mpif.h'
    integer,intent(in):: myid_md,mpi_md_world,iprint

    integer:: isp,jsp,ierr
    real(8):: a,rho,c
    character(len=128):: cline,fname
    character(len=3):: cspi,cspj

    if( myid_md.eq.0 ) then
      fname = trim(paramsdir)//'/'//trim(cprmfname)
      open(ioprms,file=trim(fname),status='old')
      interact(1:nspmax,1:nspmax) = .false.
      buck_a(1:nspmax,1:nspmax) = 0d0
      buck_rho(1:nspmax,1:nspmax) = 0d0
      buck_c(1:nspmax,1:nspmax) = 0d0
      if( iprint.ge.ipl_basic ) write(6,'(/,a)') ' Buckingham parameters:'
      do while(.true.)
        read(ioprms,*,end=10) cline
        if( cline(1:1).eq.'#' .or. cline(1:1).eq.'!' ) cycle
        backspace(ioprms)
        read(ioprms,*) cspi,cspj,a,rho,c
        isp = csp2isp(cspi)
        jsp = csp2isp(cspj)
!!$        print *,'isp,jsp,a,rho,c=',isp,jsp,a,rho,c
        if( isp.gt.0 .and. jsp.gt.0 ) then
          buck_a(isp,jsp) = a
          buck_rho(isp,jsp) = rho
          buck_c(isp,jsp) = c
          interact(isp,jsp) = .true.
          if( iprint.ge.ipl_basic ) then
            write(6,'(a,2a4,3f10.3)') '   cspi,cspj,A,rho,C = ',trim(cspi),trim(cspj),a,rho,c
          endif
!.....Symmetrize parameters
          buck_a(jsp,isp) = buck_a(isp,jsp)
          buck_rho(jsp,isp)= buck_rho(isp,jsp)
          buck_c(jsp,isp)= buck_c(isp,jsp)
          interact(jsp,isp)= interact(isp,jsp)
        else
          if( iprint.ge.ipl_info ) then
            print *,' Buckingham parameter read but not used: cspi,cspj=',cspi,cspj
          endif
        endif
      enddo
10    close(ioprms)
      if( iprint.ge.ipl_debug ) then
        write(6,'(a)') ' Finished reading '//trim(fname)
        write(6,*) ''
      endif
    endif

    call mpi_bcast(buck_a,nspmax*nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(buck_rho,nspmax*nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(buck_c,nspmax*nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(interact,nspmax*nspmax,mpi_logical,0,mpi_md_world,ierr)
    
  end subroutine read_params_Buckingham
end module Buckingham
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
