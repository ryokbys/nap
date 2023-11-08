module Morse
!-----------------------------------------------------------------------
!                     Last modified: <2023-11-01 12:29:04 KOBAYASHI Ryo>
!-----------------------------------------------------------------------
!  Parallel implementation of Morse pontential.
!    - For BVS, see Adams & Rao, Phys. Status Solidi A 208, No.8 (2011)
!    - Currently no cutoff tail treatment is done. (170310)
!-----------------------------------------------------------------------
  use pmdvars,only: nspmax,rc
  use util,only: csp2isp, num_data
  use memory,only: accum_mem
  use vector,only: dot
  implicit none
  include "./const.h"
  save
  character(len=128):: paramsdir = '.'
  character(len=128),parameter:: paramsfname = 'in.params.Morse'
  character(len=128),parameter:: configfname = 'in.config.Morse'
  character(len=128),parameter:: descfname   = 'in.desc.vcMorse'
  character(len=128),parameter:: paramsvcfname  = 'in.params.vcMorse'

  integer,parameter:: ioprms = 20
  integer,parameter:: iocnfg = 21
  integer,parameter:: iodesc = 22
  integer,parameter:: ioprmsvc = 23

!.....Max number of species available in this potential
!!$  integer,parameter:: nspmax = 9
  integer:: nsp
!.....Morse parameters
  real(8):: alp(nspmax,nspmax),d0(nspmax,nspmax),rmin(nspmax,nspmax), &
       rcs(nspmax,nspmax),rc2s(nspmax,nspmax)
  logical:: interact(nspmax,nspmax)
!.....fixed_bond Morse parameters
  real(8):: rc_fb(nspmax,nspmax)
  logical:: lfixbond(nspmax,nspmax)

!.....Smooth cutoff
  real(8):: vrcs(nspmax,nspmax), dvdrcs(nspmax,nspmax)

  integer,parameter:: ivoigt(3,3)= &
       reshape((/ 1, 6, 5, 6, 2, 4, 5, 4, 3 /),shape(ivoigt))

  real(8),allocatable:: strsl(:,:,:)
  real(8),allocatable:: ge_alp(:,:),ge_d0(:,:),ge_rmin(:,:)
  real(8),allocatable:: gf_alp(:,:,:,:),gf_d0(:,:,:,:),gf_rmin(:,:,:,:)
  real(8),allocatable:: gs_alp(:,:,:),gs_d0(:,:,:),gs_rmin(:,:,:)
  
!.....Atomic descriptor
  type atdesc
    integer:: na            ! atomic number
    character(len=3):: csym ! Symbol
    real(8):: eion1, eion2  ! 1st and 2nd ionization energies (eV)
    real(8):: eaff          ! electron affinity (eV)
    real(8):: atrad         ! atomic radius (Ang)
    real(8):: enpaul        ! Pauling's electronegativity (Pauling unit)
  end type atdesc

  type(atdesc),allocatable:: atdescs(:)

!.....Varible-charge Morse parameters
!     Currently, descriptors are:
!       - charge
!       - IE1: 1st ionization energy
!       - IE2: 2nd
!       - EA: electron affinity
!       - raidus: atomic radius
!       - EN_Pauling: electronegativity of Pauling's scale
!       - na: atomic number
  integer,parameter:: ndesc = 7
  integer,parameter:: nprm  = ndesc*2  ! 2 for each descriptor
  real(8):: walp(0:nprm),wd(0:nprm),wrmin(0:nprm),&
       pdij(0:nprm)
!.....Derivatives w.r.t. weights
  real(8):: gwalp(0:nprm),gwd(0:nprm),gwrmin(0:nprm)

!.....beta: ratio to D_ij such that the absolute value of Eij
!     should be smaller than beta*Dij
  real(8),parameter:: beta = 0.1d0
  real(8):: prefbeta

  logical:: lprmset_Morse = .false.

!.....params
  integer:: nprms
  real(8),allocatable:: params(:)

!.....Limit of exponential term to avoid Inf and NaN...
!.....The term larger than this is replaced by square of r
  real(8),parameter:: ecore = 2.d0
  real(8),parameter:: ln_ecore = log(ecore)

contains
  subroutine force_Morse(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rct,lspr &
       ,mpi_md_world,myid,epi,epot,nismax,lstrs,iprint,l1st)
    use util,only: itotOf
    implicit none
    include "mpif.h"
    include "./params_unit.h"
!!$    include "params_BVS_Morse.h"
    integer,intent(in):: namax,natm,nnmax,nismax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),lspr(0:nnmax,namax),nex(3)
    integer,intent(in):: mpi_md_world,myid
    real(8),intent(in):: ra(3,namax),h(3,3),hi(3,3),rct &
         ,tag(namax),sv(3,6)
    real(8),intent(inout):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: l1st
    logical:: lstrs

    integer:: i,j,k,l,m,n,ierr,is,js,ixyz,jxyz
    real(8):: xi(3),xj(3),xij(3),rij(3),dij,diji,dedr,epott &
         ,dxdi(3),dxdj(3),x,y,z,epotl,at(3),tmp,tmp2,texp &
         ,d0ij,alpij,rminij,dij2,vrc,dvdrc,rc,rc2
!!$    real(8),save:: rc2
    real(8),external:: fcut1,dfcut1

    if( l1st ) then
!!$      call init_Morse(natm,tag,mpi_md_world)
!!$      call read_params_Morse(myid,mpi_md_world,iprint)
      if( allocated(strsl) ) then
        call accum_mem('force_Morse',-8*size(strsl))
        deallocate(strsl)
      endif
      allocate(strsl(3,3,namax))
      call accum_mem('force_Morse',8*size(strsl))
!!$      rc2 = rc*rc
!.....Initialize smooth cutoff
      vrcs(:,:) = 0d0
      dvdrcs(:,:) = 0d0
      do is=1,nspmax
        do js=is,nspmax
          rc = rcs(is,js)
          rc2s(is,js) = rc**2
          rc2s(js,is) = rc**2
          rminij = rmin(is,js)
          if( rmin(is,js).lt.0d0 ) cycle
          alpij = alp(is,js)
          d0ij = d0(is,js)
          texp = exp( alpij*(rminij -rc))
          vrc = d0ij*( (texp-1d0)**2 -1d0 )
          vrcs(is,js) = vrc
          vrcs(js,is) = vrc
          dvdrc = 2d0 *alpij *d0ij *texp *(1d0-texp)
          dvdrcs(is,js) = dvdrc
          dvdrcs(js,is) = dvdrc
        enddo
      enddo
    endif

    if( .not.allocated(strsl) ) then
      allocate(strsl(3,3,namax))
      call accum_mem('force_Morse',8*size(strsl))
    else if( size(strsl).lt.3*3*namax ) then
      call accum_mem('force_Morse',-8*size(strsl))
      deallocate(strsl)
      allocate(strsl(3,3,namax))
      call accum_mem('force_Morse',8*size(strsl))
    endif

    epotl= 0d0
    strsl(1:3,1:3,1:namax) = 0d0

!.....Loop over resident atoms
!$omp parallel
!$omp do private(i,xi,is,k,j,js,xj,xij,rij,dij2,dij,diji,dxdi, &
!$omp     d0ij,alpij,rminij,rc2,rc,vrc,dvdrc,texp,tmp,tmp2,dedr,ixyz,jxyz) &
!$omp     reduction(+:epotl)
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is=int(tag(i))
      do k=1,lspr(0,i)
        j=lspr(k,i)
!!$        if( j.le.i ) cycle
        js= int(tag(j))
!.....Check if these two species interact
        if( .not. interact(is,js) ) cycle
        xj(1:3)= ra(1:3,j)
        xij(1:3)= xj(1:3)-xi(1:3)
        rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij2 = rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
        rc2 = rc2s(is,js)
        if( dij2.gt.rc2 ) cycle
        dij= sqrt(dij2)
        diji= 1d0/dij
        dxdi(1:3)= -rij(1:3)*diji
        d0ij = d0(is,js)
        alpij= alp(is,js)
        rminij=rmin(is,js)
        vrc = vrcs(is,js)
        rc = rcs(is,js)
        dvdrc = dvdrcs(is,js)
        texp = exp(alpij*(rminij-dij))
!.....potential
        tmp= d0ij*((texp-1d0)**2 -1d0)
        tmp2 = 0.5d0 *(tmp -vrc -dvdrc*(dij-rc))
        epi(i)= epi(i) +tmp2
        epotl= epotl +tmp2
!!$        if( j.le.natm ) then
!!$          epi(i)= epi(i) +tmp2
!!$          epotl= epotl +tmp2 +tmp2
!!$!$omp atomic
!!$          epi(j)= epi(j) +tmp2
!!$        else
!!$          epi(i)= epi(i) +tmp2
!!$          epotl= epotl +tmp2
!!$        endif
!.....force
        dedr= 2d0 *alpij *d0ij *texp *(1d0 -texp) -dvdrc
        aa(1:3,i)= aa(1:3,i) -dxdi(1:3)*dedr
!!$        do ixyz=1,3
!!$!$omp atomic
!!$          aa(ixyz,i)= aa(ixyz,i) +dxdi(ixyz)*dedr
!!$        enddo
!.....stress
        do ixyz=1,3
          do jxyz=1,3
            strsl(jxyz,ixyz,i)= strsl(jxyz,ixyz,i) &
                 -0.5d0 *dedr*rij(ixyz)*(-dxdi(jxyz))
!!$!$omp atomic
!!$            strsl(jxyz,ixyz,j)= strsl(jxyz,ixyz,j) &
!!$                 -0.5d0 *dedr*rij(ixyz)*(-dxdi(jxyz))
          enddo
        enddo
      enddo
    enddo
!$omp end do
!$omp end parallel

    strs(1:3,1:3,1:natm)= strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)

!-----gather epot
    epott= 0d0
    call mpi_allreduce(epotl,epott,1,MPI_REAL8 &
         ,MPI_SUM,mpi_md_world,ierr)
    epot= epot +epott
    if( iprint.ge.ipl_info ) print *,'epot Morse = ',epott
    return
  end subroutine force_Morse
!=======================================================================
  subroutine force_Morse_repul(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid,epi,epot,nismax,lstrs,iprint,l1st)
!
!  Repulsive-only Morse potential
!
    implicit none
    include "mpif.h"
    include "./params_unit.h"
    integer,intent(in):: namax,natm,nnmax,nismax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),lspr(0:nnmax,namax),nex(3)
    integer,intent(in):: mpi_md_world,myid
    real(8),intent(in):: ra(3,namax),h(3,3),hi(3,3),rc &
         ,tag(namax),sv(3,6)
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: l1st
    logical:: lstrs

    integer:: i,j,k,l,m,n,ierr,is,js,ixyz,jxyz
    real(8):: xi(3),xj(3),xij(3),rij(3),dij,diji,dedr,epott &
         ,dxdi(3),dxdj(3),x,y,z,epotl,at(3),tmp,tmp2,texp &
         ,d0ij,alpij,rminij
    real(8),external:: fcut1,dfcut1
!!$    real(8),allocatable,save:: strsl(:,:,:)

    if( l1st ) then
!!$      call init_Morse(natm,tag,mpi_md_world)
!!$      call read_params_Morse(myid,mpi_md_world,iprint)
      if( allocated(strsl) ) then
        call accum_mem('force_Morse',-8*size(strsl))
        deallocate(strsl)
      endif
      allocate(strsl(3,3,namax))
      call accum_mem('force_Morse',8*size(strsl))
    endif

    if( .not.allocated(strsl) ) then
      allocate(strsl(3,3,namax))
      call accum_mem('force_Morse',8*size(strsl))
    else if( size(strsl).lt.3*3*namax ) then
      call accum_mem('force_Morse',-8*size(strsl))
      deallocate(strsl)
      allocate(strsl(3,3,namax))
      call accum_mem('force_Morse',8*size(strsl))
    endif

    epotl= 0d0
    strsl(1:3,1:3,1:namax) = 0d0

!!$    do i=1,natm+nb
!!$      write(6,'(a,i5,3f10.5)') 'i,ra(1:3,i)=',i,ra(1:3,i)
!!$    enddo

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
        texp = exp(2d0*alpij*(rminij-dij))
!.....potential
        tmp= d0ij*texp
        tmp2 = 0.5d0 *tmp *fcut1(dij,0d0,rc)
        if( j.le.natm ) then
          epi(i)= epi(i) +tmp2
          epi(j)= epi(j) +tmp2
          epotl = epotl +tmp2 +tmp2
        else
          epi(i)= epi(i) +tmp2
          epotl = epotl +tmp2
        endif
!!$        if( i.eq.1 .and. iprint.ne.0 ) then
!!$          write(6,'(a,4i6,10es12.4)') 'i,j,is,js,dij,d0ij,alpij,rminij,texp,tmp2='&
!!$               ,i,j,is,js,dij,d0ij,alpij,rminij,texp,tmp2
!!$        endif
!.....force
        dedr= -2d0 *alpij *d0ij *texp *fcut1(dij,0d0,rc) &
             + tmp*dfcut1(dij,0d0,rc)
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
!!$      call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
!!$           ,nn,mpi_md_world,strsl,9)
      strs(1:3,1:3,1:natm)= strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)
    endif
    
!-----gather epot
    epott= 0d0
    call mpi_allreduce(epotl,epott,1,MPI_REAL8 &
         ,MPI_SUM,mpi_md_world,ierr)
    epot= epot +epott
!!$    write(6,'(a,es15.7)') ' Morse repul epott = ',epott
 
  end subroutine force_Morse_repul
!=======================================================================
  subroutine force_fbMorse(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rct,lspr &
       ,mpi_md_world,myid,epi,epot,nismax,lstrs,iprint,l1st)
!
!  Fixed-bond (fb) Morse potential
!  
!  To use fbMorse, specify "fbMorse" for potential in in.pmd file.
!  Atraction forces on specified bonds are replaced with spring force
!  whose spring-constant is given by 2.0 *alpha**2 *D. (see the note on 2023-10-13).
!  The bond pairs are specified in "in.params.Morse" file.
!
    use util,only: itotOf
    implicit none
    include "mpif.h"
    include "./params_unit.h"
!!$    include "params_BVS_Morse.h"
    integer,intent(in):: namax,natm,nnmax,nismax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),lspr(0:nnmax,namax),nex(3)
    integer,intent(in):: mpi_md_world,myid
    real(8),intent(in):: ra(3,namax),h(3,3),hi(3,3),rct &
         ,tag(namax),sv(3,6)
    real(8),intent(inout):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: l1st
    logical:: lstrs

    integer:: i,j,k,l,m,n,ierr,is,js,ixyz,jxyz
    real(8):: xi(3),xj(3),xij(3),rij(3),dij,diji,dedr,epott &
         ,dxdi(3),dxdj(3),x,y,z,epotl,at(3),tmp,tmp2,texp &
         ,d0ij,alpij,rminij,dij2,vrc,dvdrc,rc,rc2
    real(8),external:: fcut1,dfcut1

    if( l1st ) then
!!$      call init_Morse(natm,tag,mpi_md_world)
!!$      call read_params_Morse(myid,mpi_md_world,iprint)
      if( allocated(strsl) ) then
        call accum_mem('force_fbMorse',-8*size(strsl))
        deallocate(strsl)
      endif
      allocate(strsl(3,3,namax))
      call accum_mem('force_fbMorse',8*size(strsl))
!!$      rc2 = rc*rc
!.....Initialize smooth cutoff
      vrcs(:,:) = 0d0
      dvdrcs(:,:) = 0d0
      do is=1,nspmax
        do js=is,nspmax
          rc = rcs(is,js)
          rc2s(is,js) = rc**2
          rc2s(js,is) = rc**2
          rminij = rmin(is,js)
          if( rmin(is,js).lt.0d0 ) cycle
          alpij = alp(is,js)
          d0ij = d0(is,js)
          texp = exp( alpij*(rminij -rc))
          vrc = d0ij*( (texp-1d0)**2 -1d0 )
          vrcs(is,js) = vrc
          vrcs(js,is) = vrc
          dvdrc = 2d0 *alpij *d0ij *texp *(1d0-texp)
          dvdrcs(is,js) = dvdrc
          dvdrcs(js,is) = dvdrc
        enddo
      enddo
    endif

    if( .not.allocated(strsl) ) then
      allocate(strsl(3,3,namax))
      call accum_mem('force_fbMorse',8*size(strsl))
    else if( size(strsl).lt.3*3*namax ) then
      call accum_mem('force_fbMorse',-8*size(strsl))
      deallocate(strsl)
      allocate(strsl(3,3,namax))
      call accum_mem('force_fbMorse',8*size(strsl))
    endif

    epotl= 0d0
    strsl(1:3,1:3,1:namax) = 0d0

!.....Loop over resident atoms
!$omp parallel
!$omp do private(i,xi,is,k,j,js,xj,xij,rij,dij2,dij,diji,dxdi, &
!$omp     d0ij,alpij,rminij,vrc,dvdrc,texp,tmp,tmp2,dedr,ixyz,jxyz) &
!$omp     reduction(+:epotl)
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is=int(tag(i))
      do k=1,lspr(0,i)
        j=lspr(k,i)
!!$        if( j.le.i ) cycle
        js= int(tag(j))
!.....Check if these two species interact
        if( .not. interact(is,js) ) cycle
        xj(1:3)= ra(1:3,j)
        xij(1:3)= xj(1:3)-xi(1:3)
        rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij2 = rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
        rc2 = rc2s(is,js)
        if( dij2.gt.rc2 ) cycle
        dij= sqrt(dij2)
        diji= 1d0/dij
        dxdi(1:3)= -rij(1:3)*diji
        d0ij = d0(is,js)
        alpij= alp(is,js)
        rminij=rmin(is,js)
        rc = rcs(is,js)
        if( lfixbond(is,js) .and. dij.lt.rc_fb(is,js) .and. dij.gt.rminij ) then
          tmp2 = 0.25d0 *2d0 *alpij**2 *d0ij *(dij-rminij)**2  &
               -0.5d0*(d0ij -vrc -dvdrc*(dij-rc))
          dedr = 2d0 *alpij**2 *d0ij *(dij-rminij)
        else
          vrc = vrcs(is,js)
          dvdrc = dvdrcs(is,js)
          texp = exp(alpij*(rminij-dij))
!.....potential
          tmp= d0ij*((texp-1d0)**2 -1d0)
          tmp2 = 0.5d0 *(tmp -vrc -dvdrc*(dij-rc))
!.....force
          dedr= 2d0 *alpij *d0ij *texp *(1d0 -texp) -dvdrc
        endif
        epi(i)= epi(i) +tmp2
        epotl= epotl +tmp2
        aa(1:3,i)= aa(1:3,i) -dxdi(1:3)*dedr
!.....stress
        do ixyz=1,3
          do jxyz=1,3
            strsl(jxyz,ixyz,i)= strsl(jxyz,ixyz,i) &
                 -0.5d0 *dedr*rij(ixyz)*(-dxdi(jxyz))
          enddo
        enddo
      enddo
    enddo
!$omp end do
!$omp end parallel

    strs(1:3,1:3,1:natm)= strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)

!-----gather epot
    epott= 0d0
    call mpi_allreduce(epotl,epott,1,MPI_REAL8 &
         ,MPI_SUM,mpi_md_world,ierr)
    epot= epot +epott
    if( iprint.ge.ipl_info ) print *,'epot Morse = ',epott
    return
  end subroutine force_fbMorse
!=======================================================================
  subroutine force_vcMorse(namax,natm,tag,ra,nnmax,aa,strs,chg &
       ,h,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid,epi,epot,nismax,lstrs,iprint,l1st)
!
!  Variable-charge Morse potential.
!  Morse parameters depend on atomic charges which vary time to time.
!
    implicit none
    include "mpif.h"
    include "./params_unit.h"
    integer,intent(in):: namax,natm,nnmax,nismax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),lspr(0:nnmax,namax),nex(3)
    integer,intent(in):: mpi_md_world,myid
    real(8),intent(in):: ra(3,namax),h(3,3),hi(3,3),rc &
         ,tag(namax),sv(3,6),chg(namax)
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: lstrs,l1st

    real(8),external:: fcut1,dfcut1

    integer:: i,j,k,l,m,n,ierr,is,js,ixyz,jxyz
    real(8):: dij,dedr,epott,x,y,z,epotl,tmp,texp,d0ij,alpij,rminij &
         ,chgi,chgj,tmp2,diji,rcore,dij2
!!$    real(8),allocatable,save:: strsl(:,:,:)
    type(atdesc):: atdi,atdj
    real(8),save:: rc2
    real(8),save,allocatable:: xi(:),xj(:),xij(:),rij(:)&
         ,dxdi(:),dxdj(:),at(:)
    if( .not.allocated(xi) ) allocate(xi(3),xj(3),xij(3),rij(3),&
         dxdi(3),dxdj(3),at(3) )

    if( l1st ) then
!!$      call init_vcMorse(natm,tag,mpi_md_world)
!!$      call read_params_vcMorse(myid,mpi_md_world,iprint)
      if( allocated(strsl) ) then
        call accum_mem('force_Morse',-8*size(strsl))
        deallocate(strsl)
      endif
      allocate(strsl(3,3,namax))
      call accum_mem('force_Morse',8*size(strsl))
      prefbeta = log(1d0 -sqrt(1d0 -beta))
      rc2 = rc*rc
    endif

    if( .not.allocated(strsl) ) then
      allocate(strsl(3,3,namax))
      call accum_mem('force_Morse',8*size(strsl))
    else if( size(strsl).lt.3*3*namax ) then
      call accum_mem('force_Morse',-8*size(strsl))
      deallocate(strsl)
      allocate(strsl(3,3,namax))
      call accum_mem('force_Morse',8*size(strsl))
    endif

    epotl= 0d0
    strsl(1:3,1:3,1:namax) = 0d0

!!$    do is = 1,nsp
!!$      atdi = atdescs(is)
!!$      write(6,'(a,2i5,a,5es11.3)') 'is,na,csym,eion1,eion2,eaff,atrad,enpaul =',&
!!$           is,atdi%na,atdi%csym,atdi%eion1,atdi%eion2,atdi%eaff,&
!!$           atdi%atrad,atdi%enpaul
!!$    enddo

!!$    write(6,'(a,30es10.2)') ' chg @force = ',chg(1:natm)

!.....Loop over resident atoms
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is=int(tag(i))
      chgi = chg(i)
      atdi = atdescs(is)
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
        if(j.le.i) cycle
        js= int(tag(j))
        xj(1:3)= ra(1:3,j)
        xij(1:3)= xj(1:3)-xi(1:3)
        rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij2 = rij(1)**2 +rij(2)**2 +rij(3)**2
        if( dij2.gt.rc2 ) cycle
        dij= sqrt(dij2)
!!$        if( dij.gt.rc ) cycle
        diji= 1d0/dij
        dxdi(1:3)= -rij(1:3)*diji
        dxdj(1:3)=  rij(1:3)*diji
        chgj = chg(j)
        atdj = atdescs(js)
        call make_pair_desc(chgi,chgj,atdi,atdj,pdij)
!!$        write(6,'(a,20es11.3)') 'pdij = ',pdij(0:nprm)
!.....Create Morse parameters that depend on current atom charges
        d0ij = dot(wd,pdij)
        alpij= dot(walp,pdij)
        rminij= dot(wrmin,pdij)
        d0ij = max(d0ij, 0d0)
        rminij= max(rminij, (atdi%atrad +atdj%atrad)/2)
        alpij= max(alpij, prefbeta/(rminij -rc))
!!$        if( j.le.natm ) then
!!$          write(6,*) ' i,is,j,js,d0ij,alpij,rminij= ', &
!!$               i,is,j,js,d0ij,alpij,rminij
!!$        endif
        texp = exp(alpij*(rminij-dij))
        if( texp.gt.ecore ) then
          rcore = rminij -ln_ecore/alpij
          texp = -dij*dij +rcore*rcore +ecore
          tmp= d0ij*((texp-1d0)**2 -1d0)
          dedr= -4d0 *dij *d0ij *(texp-1d0) *fcut1(dij,0d0,rc) &
               + tmp *dfcut1(dij,0d0,rc)
        else
          tmp= d0ij*((texp-1d0)**2 -1d0)
          dedr= -2d0 *alpij *d0ij *texp *(texp-1d0) *fcut1(dij,0d0,rc) &
               + tmp *dfcut1(dij,0d0,rc)
        endif
!.....potential
        tmp2 = 0.5d0 *tmp *fcut1(dij,0d0,rc)
        if( j.le.natm ) then
          epi(i)= epi(i) +tmp2
          epi(j)= epi(j) +tmp2
          epotl = epotl +tmp2 +tmp2
        else
          epi(i)= epi(i) +tmp2
          epotl = epotl +tmp2
        endif
!.....force
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
!!$      call copy_dba_bk(namax,natm,nbmax,nb,lsb,nex,lsrc,myparity &
!!$           ,nn,mpi_md_world,strsl,9)
      strs(1:3,1:3,1:natm)= strs(1:3,1:3,1:natm) +strsl(1:3,1:3,1:natm)
    endif

    
!-----gather epot
    epott= 0d0
    call mpi_allreduce(epotl,epott,1,MPI_REAL8 &
         ,MPI_SUM,mpi_md_world,ierr)
    epot= epot +epott
 
  end subroutine force_vcMorse
!=======================================================================
  subroutine qforce_vcMorse(namax,natm,tag,ra,fq,nnmax,chg &
       ,h,rc,lspr,mpi_md_world,myid,epot,iprint,l1st)
!
!  Variable-charge Morse potential.
!  Morse parameters depend on atomic charges which vary time to time.
!
    implicit none
    include "mpif.h"
    include "./params_unit.h"
    integer,intent(in):: namax,natm,nnmax,iprint
    integer,intent(in):: lspr(0:nnmax,namax)
    integer,intent(in):: mpi_md_world,myid
    real(8),intent(in):: ra(3,namax),h(3,3),rc &
         ,tag(namax),chg(namax)
    real(8),intent(out):: epot,fq(namax)
    logical,intent(in):: l1st

    integer:: i,j,k,l,m,n,ierr,is,js,ixyz,jxyz
    real(8):: dij,dedr,epott,x,y,z,epotl,tmp,texp,d0ij,alpij,rminij,rcore &
         ,chgi,chgj,dd0dq,dalpdq,drmindq,dedd0,dedalp,dedrmin,tmp2,diji
    type(atdesc):: atdi,atdj
    real(8),external:: fcut1,dfcut1
    real(8),save,allocatable:: xi(:),xj(:),xij(:),rij(:) &
         ,dxdi(:),dxdj(:),at(:)

    if( .not.allocated(xi) ) then
      allocate(xi(3),xj(3),xij(3),rij(3), &
           dxdi(3),dxdj(3),at(3) )
      prefbeta = log(1d0 -sqrt(1d0 -beta))
    endif

    epotl= 0d0

!!$    write(6,'(a,9f10.6)') 'h(1:3,1:3) in qforce_vcMorse=',h(1:3,1:3)
!!$    write(6,'(a,30es10.2)') 'chg in qforce_Morse =',chg(1:natm)

!.....Loop over resident atoms
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is=int(tag(i))
      chgi = chg(i)
      atdi = atdescs(is)
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
        if(j.le.i) cycle
        js= int(tag(j))
        xj(1:3)= ra(1:3,j)
        xij(1:3)= xj(1:3)-xi(1:3)
        rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij= sqrt(rij(1)**2 +rij(2)**2 +rij(3)**2)
        if( dij.gt.rc ) cycle
        diji= 1d0/dij
        dxdi(1:3)= -rij(1:3)*diji
        dxdj(1:3)=  rij(1:3)*diji
        chgj = chg(j)
        atdj = atdescs(js)
        call make_pair_desc(chgi,chgj,atdi,atdj,pdij)
!.....Create Morse parameters that depend on current atom charges
        d0ij = dot(wd,pdij)
        rminij= dot(wrmin,pdij)
        alpij= dot(walp,pdij)
        d0ij = max(d0ij, 0d0)
        rminij= max(rminij, (atdi%atrad +atdj%atrad)/2)
        rminij= min(rminij,rc*0.9d0)
        alpij= max(alpij, prefbeta/(rminij -rc))
        texp = exp(alpij*(rminij-dij))
!.....force on charge
        dd0dq = wd(1) +wd(2)*chgj
        dalpdq = walp(1) +walp(2)*chgj
        drmindq = wrmin(1) +wrmin(2)*chgj
        if( texp.gt.ecore ) then
          rcore = rminij -ln_ecore/alpij
          texp = -dij*dij +rcore*rcore +ecore
          dedalp = 2d0*d0ij*(texp-1d0)*ecore*(rminij-rcore)
          dedrmin = 2d0*d0ij*(texp-1d0)*ecore*alpij
!!$          print '(a,2i5,10es11.3)', 'i,j,alpij,d0ij,rminij,rcore,texp,dedalp,dedrmin=' &
!!$               ,i,j,alpij,d0ij,rminij,rcore,texp,dedalp,dedrmin
!!$          print '(a,3es12.4)', 'alpij,d0ij,rminij = ',alpij,d0ij,rminij
!!$          print '(a,3es12.4)', 'rcore,texp = ',rcore,texp
!!$          print '(a,3es12.4)', 'dedalp,dedrmin = ',dedalp,dedrmin
        else
          dedalp = 2d0*d0ij*(texp-1d0)*texp*(rminij-dij)
          dedrmin = 2d0*d0ij*(texp-1d0)*texp*alpij
        endif
        dedd0 = (texp -1d0)**2 -1d0
        fq(i) = fq(i) +(dd0dq*dedd0 +dalpdq*dedalp +drmindq*dedrmin) &
             *fcut1(dij,0d0,rc)
!!$        print '(a,4i5,7es11.3)','i,is,j,js,dd0dq,dedd0,dalpdq,dedalp,drmindq,dedrmin,fqi='&
!!$             ,i,is,j,js,dd0dq,dedd0,dalpdq,dedalp,drmindq,dedrmin,fq(i)
!.....potential
        tmp= 0.5d0 * d0ij*((texp-1d0)**2 -1d0)
        tmp2 = tmp *fcut1(dij,0d0,rc)
        if( j.le.natm ) then
          epotl = epotl +tmp2 +tmp2
        else
          epotl = epotl +tmp2
        endif
      enddo
!!$      if( l1st ) write(6,'(a,i5,2es10.2)') 'i,chgi,fq(i)=',i,chgi,fq(i)
    enddo
    
!-----gather epot
    epott= 0d0
    call mpi_allreduce(epotl,epott,1,mpi_real8 &
         ,mpi_sum,mpi_md_world,ierr)
    epot= epot +epott

  end subroutine qforce_vcMorse
!=======================================================================
  subroutine read_params_Morse(myid_md,mpi_md_world,iprint,specorder)
!
!  Read pair parameters for Morse potential from file
!
    include 'mpif.h'
    integer,intent(in):: myid_md,mpi_md_world,iprint
    character(len=3),intent(in):: specorder(nspmax)
    integer:: i,j,isp,jsp,id,ierr,jerr,ndat
    character(len=128):: cline,fname,ctmp,cerr
    character(len=3):: cspi,cspj
    real(8):: d,r,a,rct

    jerr = 0
    if( myid_md.eq.0 ) then
      fname = trim(paramsdir)//'/'//trim(paramsfname)
      open(ioprms,file=trim(fname),status='old')
      interact(1:nspmax,1:nspmax) = .false.
      d0(1:nspmax,1:nspmax)= 0d0
      rmin(1:nspmax,1:nspmax)= 0d0
      alp(1:nspmax,1:nspmax)= 0d0
      rcs(1:nspmax,1:nspmax) = rc
      rc_fb(1:nspmax,1:nspmax) = -1.0
      lfixbond(1:nspmax,1:nspmax) = .false.
      if( iprint.ne.0 ) write(6,'(/,a)') ' Morse parameters:'
      do while(.true.)
        read(ioprms,'(a)',end=10) cline
        if( cline(1:1).eq.'#' .or. cline(1:1).eq.'!' ) cycle
!!$        backspace(ioprms)
!!$        read(ioprms,*) isp,jsp,d,a,r
        read(cline,*,iostat=ierr) ctmp
        if( ierr .ne. 0 ) cycle
        if( ctmp.eq.'fix_bond' ) then
          read(cline,*) ctmp, cspi,cspj,rct
          isp = csp2isp(cspi)
          jsp = csp2isp(cspj)
          if( isp.gt.0 .and. jsp.gt.0 ) then
            rc_fb(isp,jsp) = rct
            rc_fb(jsp,isp) = rct
            lfixbond(isp,jsp) = .true.
            lfixbond(jsp,isp) = .true.
            if( iprint.ge.ipl_basic ) then
              write(6,'(a,2a4,f7.3)') '   fix_bond: cspi,cspj,rc_fb = ', &
                   trim(cspi),trim(cspj),rct
            endif
          endif
        else if( ctmp.eq.'' ) then
          cycle
        else ! not fix_bond, normal entry of Morse potential parameters
          ndat = num_data(cline, ' ')
          if( ndat.eq.5 ) then
            read(cline,*) cspi,cspj,d,a,r
          else if( ndat.eq.6 ) then
            read(cline,*) cspi,cspj,d,a,r,rct
          else
            jerr = 80
            cerr = 'Number of values in a line must be 5 or 6 in in.params.Morse.'
            exit
          endif
          isp = csp2isp(cspi)
          jsp = csp2isp(cspj)
          if( isp.gt.0 .and. jsp.gt.0 ) then
            d0(isp,jsp) = d
            rmin(isp,jsp) = r
            alp(isp,jsp) = a
            interact(isp,jsp) = .true.
!.....Symmetrize parameters
            d0(jsp,isp) = d0(isp,jsp)
            rmin(jsp,isp)= rmin(isp,jsp)
            alp(jsp,isp)= alp(isp,jsp)
            interact(jsp,isp)= interact(isp,jsp)
            if( ndat.eq.6 ) then
              rcs(isp,jsp) = min(rct,rc)
              rcs(jsp,isp) = min(rct,rc)
            endif
            if( iprint.ge.ipl_basic ) then
              write(6,'(a,2a4,4f7.3)') '   cspi,cspj,D,alpha,rmin,rc = ',trim(cspi),trim(cspj),d,a,r,rcs(isp,jsp)
            endif
          else
            if( iprint.ge.ipl_info ) then
              print *,' Morse parameter read but not used: cspi,cspj=',cspi,cspj
            endif
          endif
        endif
      enddo  ! do while(.true.)

10    close(ioprms)
      if( iprint.ge.ipl_debug .and. jerr.eq.0 ) then
        write(6,'(a)') ' Finished reading '//trim(fname)
        write(6,*) ''
      endif
    endif

    call mpi_bcast(jerr,1,mpi_integer,0,mpi_md_world,ierr)
    if( jerr.gt.0 ) then
      if( myid_md.eq.0 ) print *,'Error: ',trim(cerr)
      stop
    endif

    call mpi_bcast(d0,nspmax*nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(rmin,nspmax*nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(alp,nspmax*nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(interact,nspmax*nspmax,mpi_logical,0,mpi_md_world,ierr)
    call mpi_bcast(rcs,nspmax*nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(rc_fb,nspmax*nspmax,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(lfixbond,nspmax*nspmax,mpi_logical,0,mpi_md_world,ierr)

  end subroutine read_params_Morse
!=======================================================================
  subroutine read_params_vcMorse(myid,mpi_md_world,iprint)
!
!  Read parameters for VC-Morse potential from file, in.params.vcMorse
!
    include 'mpif.h'
    integer,intent(in):: myid,mpi_md_world,iprint

    integer:: i,ierr
    character(len=128) :: cline,fname

    if( myid .eq.0 ) then
      fname = trim(paramsdir)//'/'//trim(paramsvcfname)
      open(ioprmsvc,file=trim(fname),status='old')
      read(ioprmsvc,'(a)') cline
      if( cline(1:1).eq.'!' .or. cline(1:1).eq.'#' ) then
!.....Read comment line, if needed, it can contain options
        
      else
        backspace(ioprmsvc)
      endif
!.....Read parameters hereafter
      do i=0,nprm
        read(ioprmsvc,*) walp(i),wd(i),wrmin(i)
      enddo
      close(ioprmsvc)

!!$      if( iprint.ne.0 ) then
!!$        print *, 'Finished reading '//trim(fname)
!!$        print *, ''
!!$      endif
    endif

    call mpi_bcast(walp,nprm+1,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(wd,nprm+1,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(wrmin,nprm+1,mpi_real8,0,mpi_md_world,ierr)

    lprmset_Morse = .true.
    
  end subroutine read_params_vcMorse
!=======================================================================
  subroutine set_paramsdir_Morse(dname)
!
!  Accessor routine to set paramsdir.
!
    implicit none
    character(len=*),intent(in):: dname

    paramsdir = trim(dname)
    return
  end subroutine set_paramsdir_Morse
!=======================================================================
  subroutine set_params_Morse(ndimp,params_in,ctype,interact_in)
!
!  Accessor routine to set Morse parameters from outside.
!  Curretnly this routine is supposed to be called only on serial run.
!
    use Coulomb,only: vid_bvs, npq_bvs, acc
    integer,intent(in):: ndimp
    real(8),intent(in):: params_in(ndimp)
    character(len=*),intent(in):: ctype
    logical,intent(in):: interact_in(nspmax,nspmax)

    integer:: i,j,inc,itmp,nspt,nint
    real(8):: c

    nprms = ndimp
    if( .not.allocated(params) ) allocate(params(nprms))
    params(1:nprms) = params_in(1:ndimp)
    lprmset_Morse = .true.

    interact(:,:) = interact_in(:,:)
    nint = 0
    do i=1,nspmax
      do j=i,nspmax
        if( .not.interact(i,j) ) cycle
        nint = nint +1
      enddo
    enddo

    if( index(ctype,'BVS').ne.0 ) then

      d0(:,:)= 0d0
      alp(:,:)= 0d0
      rmin(:,:)= -1d0
    
!.....In case of BVS3, only alpha and Rmin are to be optimized and
!     and D0 is computed using alpha, Rmin, vid, npq.
      if( ctype(4:4).eq.'3' ) then
        if( nint*2.ne.nprms ) then
          write(6,*) 'ERROR @set_params_Morse: nint*2.ne.nprms !!!'
          write(6,*) '  nint,nprms= ',nint,nprms
          write(6,*) 'Probably you need to set interactions correctly...'
          stop
        endif

        inc = 0
        do i=1,nspmax
          do j=i,nspmax
            if( .not.interact(i,j) ) cycle
            inc= inc +1
            alp(i,j) = params(inc)
            alp(j,i) = alp(i,j)
            inc= inc +1
            rmin(i,j) = params(inc)
            rmin(j,i) = rmin(i,j)
          enddo
        enddo

!.....Determine D0 using alpha, Rmin,vid,npq according to
!     Adams & Rao, Phys. Status Solidi A 208, No.8 (2011)
!.....c value should be determined according to anion species,
!     but in most cases c=1 (anion is s or p element)
        c = 1d0
        do i=1,nspmax
          do j=1,nspmax
            if( .not.interact(i,j) ) cycle
            d0(i,j) = c*acc*(abs(vid_bvs(i)*vid_bvs(j)))**(1d0/c)/rmin(i,j) &
                 /sqrt(dble(npq_bvs(i)*npq_bvs(j))) /(2d0*alp(i,j)**2)
            d0(j,i) = d0(i,j)
          enddo
        enddo

!.....Except BVS3, 3 parameters per pair
      else
        if( nint*3.ne.nprms ) then
          write(6,*) 'ERROR @set_params_Morse: nint*3.ne.nprms !!!'
          write(6,*) '  nint,nprms= ',nint,nprms
          write(6,*) 'Probably you need to set interactions correctly...'
          write(6,*) 'interact:'
          do i=1,nspmax
            do j=i,nspmax
              write(6,'(a,2i5,l4)') '  i,j,interact(i,j)=',i,j,interact(i,j)
            enddo
          enddo
          stop
        endif

        inc = 0
        do i=1,nspmax
          do j=i,nspmax
            if( .not.interact(i,j) ) cycle
            inc= inc +1
            d0(i,j) = params(inc)
            d0(j,i) = d0(i,j)
            inc= inc +1
            alp(i,j) = params(inc)
            alp(j,i) = alp(i,j)
            inc= inc +1
            rmin(i,j) = params(inc)
            rmin(j,i) = rmin(i,j)
          enddo
        enddo
      endif  ! ctype(4)
    endif  ! index(ctype,'BVS')

! !.....Different operations for different potential type
! !.....for example, only O-X interactions in BVS potential,
! !.....whereas all the pair interactions for normal Morse potential
!     if( trim(ctype).eq.'bvs' .or. trim(ctype).eq.'BVS' ) then
!       nspt = nprms/3 +1
! 
!       d0(1:nspt,1:nspt)= 0d0
!       alp(1:nspt,1:nspt)= 0d0
!       rmin(1:nspt,1:nspt)= 0d0
!       interact(1:nspt,1:nspt) = .false.
! 
!       inc = 0
!       do i=2,nspt
!         inc= inc +1
!         d0(1,i) = params(inc)
!         d0(i,1) = d0(1,i)
!         inc= inc +1
!         alp(1,i) = params(inc)
!         alp(i,1) = alp(1,i)
!         inc= inc +1
!         rmin(1,i) = params(inc)
!         rmin(i,1) = rmin(1,i)
!         interact(1,i) = .true.
!         interact(i,1) = .true.
!       enddo
! 
!     else  ! All the pair interactions for normal Morse potential
! !.....Number of pairs should be, 1 or 3, 6, 10, 15, 21, 28, 36, 45
!       itmp = nprms/3
!       if( itmp.eq.1 ) then
!         nspt = 1
!       else if( itmp.eq.3 ) then
!         nspt = 2
!       else if( itmp.eq.6 ) then
!         nspt = 3
!       else if( itmp.eq.10 ) then
!         nspt = 4
!       else if( itmp.eq.15 ) then
!         nspt = 5
!       else if( itmp.eq.21 ) then
!         nspt = 6
!       else if( itmp.eq.28 ) then
!         nspt = 7
!       else if( itmp.eq.36 ) then
!         nspt = 8
!       else if( itmp.eq.45 ) then
!         nspt = 9
!       else
!         print *,'ERROR: number of pairs wrong.'
!         print *,'  number of pairs extracted from nprms = ',itmp
!         stop
!       endif
! 
!       d0(1:nspt,1:nspt)= 0d0
!       alp(1:nspt,1:nspt)= 0d0
!       rmin(1:nspt,1:nspt)= 0d0
!       interact(1:nspt,1:nspt) = .false.
! 
!       inc = 0
!       do i=1,nspt
!         do j=i,nspt
!           inc= inc +1
!           d0(i,j) = params(inc)
!           d0(j,i) = d0(i,j)
!           inc= inc +1
!           alp(i,j) = params(inc)
!           alp(j,i) = alp(i,j)
!           inc= inc +1
!           rmin(i,j) = params(inc)
!           rmin(j,i) = rmin(i,j)
! !!$          write(6,'(a,2i5,3f8.4)') ' isp,jsp,d0,alp,rmin=', &
! !!$               i,j,d0(i,j),alp(i,j),rmin(i,j)
!           interact(i,j) = .true.
!           interact(j,i) = .true.
!         enddo
!       enddo
!     endif

    return
  end subroutine set_params_Morse
!=======================================================================
  subroutine set_params_vcMorse(ndimp,params)
!
! Accessor routine to set vcMorse parameters from outside.
! Curretnly this routine is supposed to be called only on serial run.
! So no need of treating this as parallel code.
!
    integer,intent(in):: ndimp
    real(8),intent(in):: params(ndimp)

    integer:: i,inc

    if( ndimp.ne.3*(nprm+1) ) then
      print *,'Error: ndimp.ne.3*(nprm+1) !!!'
      stop
    endif

    inc = 0
    do i=0,ndesc*2
      inc = inc + 1
      walp(i) = params(inc)
    enddo
    do i=0,ndesc*2
      inc = inc + 1
      wd(i) = params(inc)
    enddo
    do i=0,ndesc*2
      inc = inc + 1
      wrmin(i) = params(inc)
    enddo
    
    lprmset_Morse = .true.
    return
  end subroutine set_params_vcMorse
!=======================================================================
  subroutine update_params_Morse(ctype)
!
!  Update Morse parameters by taking parameter values from params array.
!  This routine would be called only from fitpot externally.
!
    character(len=*),intent(in):: ctype
    integer:: i,j,inc, nspt, itmp

    if( .not.lprmset_Morse ) then
      print *,'ERROR: params have not been set yet.'
      stop
    endif

!.....Different operations for different potential type
!.....for example, only O-X interactions in BVS potential,
!.....whereas all the pair interactions for normal Morse potential
    if( trim(ctype).eq.'bvs' .or. trim(ctype).eq.'BVS' ) then
!!$    if( nprms.ne.3*(nsp-1) ) then
!!$      print *,'ERROR: nprms.ne.3*(nsp-1), nprms,nsp=',nprms,nsp
!!$      stop
!!$    endif
      nspt = nprms/3 +1

      d0(1:nspt,1:nspt)= 0d0
      alp(1:nspt,1:nspt)= 0d0
      rmin(1:nspt,1:nspt)= 0d0
      interact(1:nspt,1:nspt) = .false.

      inc = 0
      do i=2,nspt
        inc= inc +1
        d0(1,i) = params(inc)
        d0(i,1) = d0(1,i)
        inc= inc +1
        alp(1,i) = params(inc)
        alp(i,1) = alp(1,i)
        inc= inc +1
        rmin(1,i) = params(inc)
        rmin(i,1) = rmin(1,i)
!!$        write(6,'(a,2i5,3f8.4)') ' isp,jsp,d0,alp,rmin=', &
!!$             1,i,d0(1,i),alp(1,i),rmin(1,i)
        interact(1,i) = .true.
        interact(i,1) = .true.
      enddo
    else  ! All the pair interactions for normal Morse potential
!.....Number of pairs should be, 1 or 3, 6, 10, 15, 21, 28, 36, 45
      itmp = nprms/3
      if( itmp.eq.1 ) then
        nspt = 1
      else if( itmp.eq.3 ) then
        nspt = 2
      else if( itmp.eq.6 ) then
        nspt = 3
      else if( itmp.eq.10 ) then
        nspt = 4
      else if( itmp.eq.15 ) then
        nspt = 5
      else if( itmp.eq.21 ) then
        nspt = 6
      else if( itmp.eq.28 ) then
        nspt = 7
      else if( itmp.eq.36 ) then
        nspt = 8
      else if( itmp.eq.45 ) then
        nspt = 9
      else
        print *,'ERROR: number of pairs wrong.'
        print *,'  number of pairs extracted from nprms = ',itmp
        stop
      endif

      d0(1:nspt,1:nspt)= 0d0
      alp(1:nspt,1:nspt)= 0d0
      rmin(1:nspt,1:nspt)= 0d0
      interact(1:nspt,1:nspt) = .false.

      inc = 0
      do i=1,nspt
        do j=i,nspt
          inc= inc +1
          d0(i,j) = params(inc)
          d0(j,i) = d0(i,j)
          inc= inc +1
          alp(i,j) = params(inc)
          alp(j,i) = alp(i,j)
          inc= inc +1
          rmin(i,j) = params(inc)
          rmin(j,i) = rmin(i,j)
!!$          write(6,'(a,2i5,3f8.4)') ' isp,jsp,d0,alp,rmin=', &
!!$               i,j,d0(i,j),alp(i,j),rmin(i,j)
          interact(i,j) = .true.
          interact(j,i) = .true.
        enddo
      enddo
    endif


    return
  end subroutine update_params_Morse
!=======================================================================
  subroutine read_element_descriptors(myid_md,mpi_md_world,iprint)
!
!  Read descriptors for elements from file.
!
    integer,intent(in):: myid_md,mpi_md_world,iprint

    integer:: itmp,isp,nspt
    character(len=5):: ctmp
    character(len=128):: cline,fname
    type(atdesc):: atd
    real(8):: eion1,eion2,eaff,atrad,enpaul

    if( allocated(atdescs) ) deallocate(atdescs)
    allocate(atdescs(nspmax))
    
    if( myid_md.eq.0 ) then
      fname = trim(paramsdir)//'/'//trim(descfname)
      open(iodesc,file=trim(fname),status='old')
      isp = 0
      do while(.true.)
        read(iodesc,*,end=10) cline
        if( cline(1:1).eq.'#' .or. cline(1:1).eq.'!' ) then
!.....Comment or command line
        else
          backspace(iodesc)
          read(iodesc,*) itmp, ctmp, eion1, eion2, eaff, atrad, enpaul
          isp = isp + 1
          if( isp.gt.nspmax ) then
            if( iprint.ne.0 ) then
              write(6,*) trim(fname)//' has more entries than NSPMAX.' &
                   //' So skip reading it.'
              write(6,*) '  NSPMAX = ',nspmax
            endif
            goto 10
          endif
          atdescs(isp)%na = itmp
          atdescs(isp)%csym = trim(ctmp)
          atdescs(isp)%eion1 = eion1
          atdescs(isp)%eion2 = eion2
          atdescs(isp)%eaff = eaff
          atdescs(isp)%atrad = atrad
          atdescs(isp)%enpaul= enpaul
        endif
      enddo
10    close(iodesc)
      nspt = isp
      if( iprint.ne.0 ) then
        write(6,*) 'Atomic descriptors were loaded from '//trim(fname)
        write(6,*) 'atomic number, symbol,    IE1,    IE2,    EA,'//&
             '  radius,  EN_Pauling'
        do isp=1,nspt
          atd = atdescs(isp)
          write(6,'(i14,5x,a3,5f8.3)') atd%na, &
               trim(atd%csym), &
               atd%eion1, atd%eion2, atd%eaff, atd%atrad, atd%enpaul
        enddo
        print *,''
      endif
    endif

!.....Share the descriptor information with all nodes
    call bcast_atdescs(nspt,myid_md,mpi_md_world)

  end subroutine read_element_descriptors
!=======================================================================
  subroutine bcast_atdescs(nspt,myid_md,mpi_md_world)
!
!  Broadcast type atdesc
!
    implicit none 
    include 'mpif.h'
    integer,intent(in):: nspt,myid_md,mpi_md_world

    integer,allocatable:: nas(:)
    character(len=3),allocatable:: csyms(:)
    real(8),allocatable:: eion1s(:),eion2s(:),eaffs(:),enpauls(:),atrads(:)
    integer:: isp,ierr

    allocate(eion1s(nspt),eion2s(nspt),eaffs(nspt),enpauls(nspt),&
         atrads(nspt),nas(nspt),csyms(nspt))

    if( myid_md.eq.0 ) then
      do isp=1,nspt
        nas(isp) = atdescs(isp)%na
        csyms(isp) = atdescs(isp)%csym
        eion1s(isp) = atdescs(isp)%eion1
        eion2s(isp) = atdescs(isp)%eion2
        eaffs(isp) = atdescs(isp)%eaff
        atrads(isp) = atdescs(isp)%atrad
        enpauls(isp) = atdescs(isp)%enpaul
      enddo
    endif

    call mpi_bcast(nas,nspt,mpi_integer,0,mpi_md_world,ierr)
    call mpi_bcast(csyms,nspt,mpi_character,0,mpi_md_world,ierr)
    call mpi_bcast(eion1s,nspt,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(eion2s,nspt,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(eaffs,nspt,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(atrads,nspt,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(enpauls,nspt,mpi_real8,0,mpi_md_world,ierr)
    
    deallocate(eion1s,eion2s,eaffs,enpauls,atrads,nas,csyms)
  end subroutine bcast_atdescs
!=======================================================================
  subroutine make_pair_desc(chgi,chgj,atdi,atdj,pdij)
!
!  Make a pair descriptor vector from atomic descriptros
!
    real(8),intent(in):: chgi,chgj
    type(atdesc),intent(in):: atdi,atdj
    real(8),intent(inout):: pdij(0:nprm)

    integer:: i
    real(8):: di(ndesc),dj(ndesc)

    di(1) = chgi
    di(2) = atdi%eion1
    di(3) = atdi%eion2
    di(4) = atdi%eaff
    di(5) = atdi%atrad
    di(6) = atdi%enpaul
    di(7) = atdi%na

    dj(1) = chgj
    dj(2) = atdj%eion1
    dj(3) = atdj%eion2
    dj(4) = atdj%eaff
    dj(5) = atdj%atrad
    dj(6) = atdj%enpaul
    dj(7) = atdj%na

    pdij(0) = 1d0  ! bias
    do i=1,ndesc
      pdij(2*i-1) = di(i)+dj(i)
      pdij(2*i)   = di(i)*dj(i)
    enddo
    return
    
  end subroutine make_pair_desc
!=======================================================================
  subroutine gradw_Morse(namax,natm,tag,ra,nnmax &
       ,h,rc,lspr,epot,iprint,ndimp,gwe,gwf,gws &
       ,lematch,lfmatch,lsmatch,iprm0)
!
!  Gradient w.r.t. weights.
!  Note: This routine is always called in single run,
!  thus no need of parallel implementation.
!  - iprm0: The starting point -1 in parameter array for this FF.
!
    implicit none
    include "./params_unit.h"
    integer,intent(in):: namax,natm,nnmax,iprint,iprm0
    integer,intent(in):: lspr(0:nnmax,namax)
    real(8),intent(in):: ra(3,namax),h(3,3),rc,tag(namax)
    real(8),intent(inout):: epot
    integer,intent(in):: ndimp
    real(8),intent(inout):: gwe(ndimp),gwf(3,ndimp,natm),gws(6,ndimp)
    logical,intent(in):: lematch,lfmatch,lsmatch

    integer:: i,j,k,l,m,n,jj,ierr,is,js,ixyz,jxyz,inc,nspt &
         ,ne,nf,ns
    real(8):: dij,dedr,rc2,dij2 &
         ,x,y,z,epotl,tmp,texp,d0ij,alpij,rminij &
         ,dd0dq,dalpdq,drmindq,dedd0,dedalp,dedrmin,tmp2 &
         ,diji,factor,fc1,dfc1,dr,dedrd0,dedralp,dedrrmin
    real(8),allocatable,save:: xi(:),xj(:),xij(:),rij(:),dxdi(:),dxdj(:)
    real(8),external:: fcut1,dfcut1

    if( .not. allocated(ge_alp) ) then
      allocate(ge_alp(nspmax,nspmax),ge_d0(nspmax,nspmax),ge_rmin(nspmax,nspmax))
      allocate(gs_alp(nspmax,nspmax,6),gs_d0(nspmax,nspmax,6), &
           gs_rmin(nspmax,nspmax,6))
    endif
    if( .not.allocated(xi) ) then
      allocate(xi(3),xj(3),xij(3),rij(3),dxdi(3),dxdj(3))
    endif

    if( .not.allocated(gf_alp) &
         .or. size(gf_alp).ne.nspmax*nspmax*3*natm ) then
      if( allocated(gf_alp) ) then
        call accum_mem('force_Morse',-8*(size(gf_alp)+size(gf_d0)+size(gf_rmin)))
        deallocate(gf_alp,gf_d0,gf_rmin)
      endif
      allocate(gf_alp(nspmax,nspmax,3,natm),gf_d0(nspmax,nspmax,3,natm), &
           gf_rmin(nspmax,nspmax,3,natm))
      call accum_mem('force_Morse',8*(size(gf_alp)+size(gf_d0)+size(gf_rmin)))
    endif

!.....Set nsp by max isp of atoms in the system
    nsp = 0
    do i=1,natm
      nsp = max(nsp,int(tag(i)))
    enddo

    rc2 = rc*rc

    ge_alp(1:nspmax,1:nspmax) = 0d0
    ge_d0(1:nspmax,1:nspmax) = 0d0
    ge_rmin(1:nspmax,1:nspmax) = 0d0
    gf_alp(1:nspmax,1:nspmax,1:3,1:natm) = 0d0
    gf_d0(1:nspmax,1:nspmax,1:3,1:natm) = 0d0
    gf_rmin(1:nspmax,1:nspmax,1:3,1:natm) = 0d0
    gs_alp(1:nspmax,1:nspmax,1:6) = 0d0
    gs_d0(1:nspmax,1:nspmax,1:6) = 0d0
    gs_rmin(1:nspmax,1:nspmax,1:6) = 0d0
    
!.....Loop over resident atoms
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is=int(tag(i))
      do jj=1,lspr(0,i)
        j=lspr(jj,i)
        if(j.eq.0) exit
        if(j.le.i) cycle
        js= int(tag(j))
!.....Check if these two species interact
        if( .not. interact(is,js) ) cycle
        xj(1:3)= ra(1:3,j)
        xij(1:3)= xj(1:3)-xi(1:3)
        rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij2= rij(1)**2 +rij(2)**2 +rij(3)**2
        if( dij2.gt.rc2 ) cycle
        dij= sqrt(dij2)
        diji = 1d0/dij
        dxdi(1:3)= -rij(1:3)*diji
        dxdj(1:3)=  rij(1:3)*diji
        d0ij = d0(is,js)
        alpij= alp(is,js)
        rminij=rmin(is,js)
        dr = rminij-dij
        texp = exp(alpij*dr)
!.....Potential
        fc1 = fcut1(dij,0d0,rc)
        dfc1= dfcut1(dij,0d0,rc)
!!$        tmp= 0.5d0 * d0ij*((texp-1d0)**2 -1d0)
        if( j.le.natm ) then
          factor = 1d0
        else
          factor = 0.5d0
        endif
        if( lematch ) then
!!$        tmp2 = tmp *fc1
!!$          if( j.le.natm ) then
!!$            epotl = epotl +tmp2 +tmp2
!!$          else
!!$            epotl = epotl +tmp2
!!$          endif
!.....Derivative of potential energy w.r.t. {w}
          dedd0 = ((texp -1d0)**2 -1d0) *fc1 *factor 
          dedalp = 2d0*d0ij*(texp-1d0)*texp*dr *fc1*factor
          dedrmin = 2d0*d0ij*(texp-1d0)*texp*alpij *fc1*factor
          ge_d0(is,js) = ge_d0(is,js) +dedd0
          ge_alp(is,js) = ge_alp(is,js) +dedalp
          ge_rmin(is,js) = ge_rmin(is,js) +dedrmin
        endif
!.....Pre-compute some factors required in force and stress derivatives
        if( lfmatch .or. lsmatch ) then
          dedr= 2d0 *alpij *d0ij *texp *(1d0 -texp) *fc1 &
               + dfc1 *d0ij*((texp-1d0)**2-1d0)
          dedrd0  = dedr/d0ij
          dedralp = 2d0*d0ij*texp*((1d0-texp) +alpij*dr*(1d0-2d0*texp))*fc1 &
               +dfc1*(2d0*d0ij*(texp-1d0)*texp*dr)
          dedrrmin= 2d0*d0ij*texp*alpij*alpij*( 1d0 -2d0*texp)*fc1 &
               +dfc1*(2d0*d0ij*(texp-1d0)*texp*alpij)
        endif
        if( lfmatch ) then
!.....Force
!!$          aa(1:3,i)= aa(1:3,i) -dxdi(1:3)*dedr
!!$          aa(1:3,j)= aa(1:3,j) -dxdj(1:3)*dedr
!.....Derivative of forces
          gf_d0(is,js,1:3,i)= gf_d0(is,js,1:3,i) -dxdi(1:3) *dedrd0
          gf_alp(is,js,1:3,i)= gf_alp(is,js,1:3,i) -dxdi(1:3) *dedralp
          gf_rmin(is,js,1:3,i)= gf_rmin(is,js,1:3,i) -dxdi(1:3) *dedrrmin
          if( j.le.natm ) then
            gf_d0(is,js,1:3,j)= gf_d0(is,js,1:3,j) -dxdj(1:3)*dedrd0
            gf_alp(is,js,1:3,j)= gf_alp(is,js,1:3,j) -dxdj(1:3) *dedralp
            gf_rmin(is,js,1:3,j)= gf_rmin(is,js,1:3,j) -dxdj(1:3) *dedrrmin
          endif
        endif
!.....Stress
        if( lsmatch ) then
!.....Derivative of stress
          do ixyz=1,3
            do jxyz=1,3
              k = ivoigt(ixyz,jxyz)
!!$              strsl(jxyz,ixyz,i)= strsl(jxyz,ixyz,i) &
!!$                   -0.5d0 *dedr*rij(ixyz)*(-dxdi(jxyz))
!!$              strsl(jxyz,ixyz,j)= strsl(jxyz,ixyz,j) &
!!$                   -0.5d0 *dedr*rij(ixyz)*(-dxdi(jxyz))
              gs_d0(is,js,k)= gs_d0(is,js,k) &
                   -factor *rij(ixyz)*(-dxdi(jxyz)) *dedrd0
              gs_alp(is,js,k)= gs_alp(is,js,k) &
                   -factor *rij(ixyz) *(-dxdi(jxyz)) *dedralp
              gs_rmin(is,js,k)= gs_rmin(is,js,k) &
                   -factor *rij(ixyz) *(-dxdi(jxyz)) *dedrrmin
!!$              if( j.le.natm ) then
!!$                gs_d0(is,js,k)= gs_d0(is,js,k) &
!!$                     -0.5d0 *rij(ixyz)*(-dxdi(jxyz)) *dedrd0
!!$                gs_alp(is,js,k)= gs_alp(is,js,k) &
!!$                     -0.5d0 *rij(ixyz) *(-dxdi(jxyz)) *dedralp
!!$                gs_rmin(is,js,k)= gs_rmin(is,js,k) &
!!$                     -0.5d0 *rij(ixyz) *(-dxdi(jxyz)) *dedrrmin
!!$              endif
            enddo
          enddo
        endif
      enddo
    enddo

!!$    do is=1,nsp
!!$      do js=is,nsp
!!$        if( .not.interact(is,js) ) cycle
!!$        write(6,'(a,2i4,3es12.4)') ' is,js,ge_alp,d0,rmin=', &
!!$             is,js,ge_alp(is,js),ge_d0(is,js),ge_rmin(is,js)
!!$      enddo
!!$    enddo

!.....Tidy up gradient arrays
    if( lematch ) then
      ne = iprm0
      do is=1,nsp
        do js=is,nsp
          if( .not. interact(is,js) ) cycle
          ne = ne + 1
          gwe(ne) = gwe(ne) +ge_d0(is,js)
          if( is.ne.js ) gwe(ne) = gwe(ne) +ge_d0(js,is)
!!$          print '(a,2i4,3es11.3)','is,js,ge_d0(is,js),(js,is),gwe='&
!!$               ,is,js,ge_d0(is,js),ge_d0(js,is),gwe(ne)
          ne = ne + 1
          gwe(ne) = gwe(ne) +ge_alp(is,js)
          if( is.ne.js ) gwe(ne) = gwe(ne) +ge_alp(js,is)
          ne = ne + 1
          gwe(ne) = gwe(ne) +ge_rmin(is,js)
          if( is.ne.js ) gwe(ne) = gwe(ne) +ge_rmin(js,is)
        enddo
      enddo
    endif
    if( lfmatch ) then
      do i=1,natm
        nf = iprm0
        do is=1,nsp
          do js=is,nsp
            if( .not. interact(is,js) ) cycle
            nf = nf + 1
            gwf(1:3,nf,i)=gwf(1:3,nf,i) +gf_d0(is,js,1:3,i)
            if( is.ne.js ) gwf(1:3,nf,i)=gwf(1:3,nf,i) +gf_d0(js,is,1:3,i)
            nf = nf + 1
            gwf(1:3,nf,i)=gwf(1:3,nf,i) +gf_alp(is,js,1:3,i)
            if( is.ne.js ) gwf(1:3,nf,i)=gwf(1:3,nf,i) +gf_alp(js,is,1:3,i)
            nf = nf + 1
            gwf(1:3,nf,i)=gwf(1:3,nf,i) +gf_rmin(is,js,1:3,i)
            if( is.ne.js ) gwf(1:3,nf,i)=gwf(1:3,nf,i) +gf_rmin(js,is,1:3,i)
          enddo
        enddo
      enddo
    endif
    if( lsmatch ) then
      do i=1,natm
        ns = iprm0
        do is=1,nsp
          do js=is,nsp
            if( .not. interact(is,js) ) cycle
            ns = ns + 1
            do k=1,6
              gws(k,ns)=gws(k,ns) +gs_d0(is,js,k)
              if( is.ne.js ) gws(k,ns)=gws(k,ns) +gs_d0(js,is,k)
            enddo
            ns = ns + 1
            do k=1,6
              gws(k,ns)=gws(k,ns) +gs_alp(is,js,k)
              if( is.ne.js ) gws(k,ns)=gws(k,ns) +gs_alp(js,is,k)
            enddo
            ns = ns + 1
            do k=1,6
              gws(k,ns)=gws(k,ns) +gs_rmin(is,js,k)
              if( is.ne.js ) gws(k,ns)=gws(k,ns) +gs_rmin(js,is,k)
            enddo
          enddo
        enddo
      enddo
    endif

    return
  end subroutine gradw_Morse
!=======================================================================
  subroutine gradw_vcMorse(namax,natm,tag,ra,nnmax,chg &
       ,h,rc,lspr,epot,iprint,ndimp,gwe,gwf,gws)
!
!  Derivative w.r.t.parameters {w}.
!  Note: This routine is always called in single run,
!  thus no need of parallel implementation.
!
    implicit none
    include "./params_unit.h"
    integer,intent(in):: namax,natm,nnmax,iprint
    integer,intent(in):: lspr(0:nnmax,namax)
    real(8),intent(in):: ra(3,namax),h(3,3),rc &
         ,tag(namax),chg(namax)
    real(8),intent(inout):: epot
    integer,intent(in):: ndimp
    real(8),intent(inout):: gwe(ndimp),gwf(3,ndimp,natm),gws(6,ndimp)

    integer:: i,j,k,l,m,n,ierr,is,js,ixyz,jxyz,inc
    real(8):: xi(3),xj(3),xij(3),rij(3),dij,dedr &
         ,x,y,z,epotl,tmp,texp,d0ij,alpij,rminij &
         ,chgi,chgj,dd0dq,dalpdq,drmindq,dedd0,dedalp,dedrmin,tmp2 &
         ,dr
    type(atdesc):: atdi,atdj
    real(8),external:: fcut1,dfcut1

    epotl= 0d0

    gwalp(0:nprm) = 0d0
    gwd(0:nprm) = 0d0
    gwrmin(0:nprm) = 0d0

!!$    write(6,'(a,30es16.8)') ' chg @pderiv = ',chg(1:natm)
    
!.....Loop over resident atoms
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is=int(tag(i))
      chgi = chg(i)
      atdi = atdescs(is)
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
!!$        if(j.le.i) cycle
        js= int(tag(j))
        xj(1:3)= ra(1:3,j)
        xij(1:3)= xj(1:3)-xi(1:3)
        rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij= sqrt(rij(1)**2 +rij(2)**2 +rij(3)**2)
        if( dij.gt.rc ) cycle
        chgj = chg(j)
        atdj = atdescs(js)
        call make_pair_desc(chgi,chgj,atdi,atdj,pdij)
!.....Create Morse parameters that depend on current atom charges
        d0ij = dot(wd,pdij)
        alpij= dot(walp,pdij)
        rminij= dot(wrmin,pdij)
        d0ij = max(d0ij, 0d0)
        rminij= max(rminij, (atdi%atrad +atdj%atrad)/2)
        alpij= max(alpij, prefbeta/(rminij -rc))
        dr = rminij-dij
        texp = exp(alpij*dr)
!!$        write(6,'(a,4i5,4es15.7)') 'i,is,j,js,d0ij,alpij,rminij,texp=',&
!!$             i,is,j,js,d0ij,alpij,rminij,texp
!!$        write(6,'(a,13f8.3)') 'pdij=',pdij(0:nprm)
!.....Potential
        tmp= 0.5d0 * d0ij*((texp-1d0)**2 -1d0)
        tmp2 = tmp *fcut1(dij,0d0,rc)
        epotl = epotl +tmp2
!.....Derivative of potential energy w.r.t. {w}
        dedd0 = ((texp -1d0)**2 -1d0)
        dedalp = 2d0*d0ij*(texp-1d0)*texp*dr
        dedrmin = 2d0*d0ij*(texp-1d0)*texp*alpij
        gwd(0:nprm) = gwd(0:nprm) +0.5d0 *dedd0 *pdij(0:nprm)
        gwalp(0:nprm) = gwalp(0:nprm) +0.5d0 *dedalp *pdij(0:nprm)
        gwrmin(0:nprm) = gwrmin(0:nprm) +0.5d0 *dedrmin *pdij(0:nprm)
!!$        write(6,'(a,2i5,3es15.7)') 'i,j,gwalp(6),dedalp,pdij(6)=', &
!!$             i,j,gwalp(6),dedalp,pdij(6)
      enddo
    enddo

!!$!.....Not parallel
!!$    epot= epotl

!!$!.....gwalp,gwd,gwrmin to pderiv
!!$    inc = 0
!!$    do i=0,ndesc*2
!!$      inc = inc + 1
!!$      pderiv(inc) = gwalp(i)
!!$    enddo
!!$    do i=0,ndesc*2
!!$      inc = inc + 1
!!$      pderiv(inc) = gwd(i)
!!$    enddo
!!$    do i=0,ndesc*2
!!$      inc = inc + 1
!!$      pderiv(inc) = gwrmin(i)
!!$    enddo
    return
  end subroutine gradw_vcMorse
!=======================================================================

end module Morse
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd lib"
!     End:
