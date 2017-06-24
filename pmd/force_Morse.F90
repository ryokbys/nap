module Morse
!-----------------------------------------------------------------------
!                     Last modified: <2017-06-24 21:49:19 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
!  Parallel implementation of Morse pontential.
!    - For BVS, see Adams & Rao, Phys. Status Solidi A 208, No.8 (2011)
!    - Currently no cutoff tail treatment is done. (170310)
!-----------------------------------------------------------------------
  implicit none
  character(len=128):: paramsdir = '.'
  character(len=128),parameter:: paramsfname = 'in.params.Morse'
  character(len=128),parameter:: configfname = 'in.config.Morse'
  character(len=128),parameter:: descfname   = 'in.desc.vcMorse'
  character(len=128),parameter:: paramsvcfname  = 'in.params.vcMorse'

  integer,parameter:: ioprms = 20
  integer,parameter:: iocnfg = 21
  integer,parameter:: iodesc = 22
  integer,parameter:: ioprmsvc = 23

!.....Number of species
  integer:: nsp
!.....Morse parameters
  real(8),allocatable:: alp(:,:),d0(:,:),rmin(:,:)
  logical,allocatable:: interact(:,:)

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
  integer,parameter:: ndesc = 6
  integer,parameter:: nprm  = ndesc*2  ! 2 for each descriptor
  real(8):: walp(0:nprm),wd(0:nprm),wrmin(0:nprm),&
       pdij(0:nprm)
!.....Derivatives w.r.t. weights
  real(8):: gwalp(0:nprm),gwd(0:nprm),gwrmin(0:nprm)

  logical:: lprmset = .false.
  
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
      call init_Morse(natm,tag,mpi_md_world)
      call read_params_Morse(myid,mpi_md_world,iprint)
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
  subroutine force_vcMorse(namax,natm,tag,ra,nnmax,aa,strs,chg &
       ,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid,epi,epot,nismax,acon,lstrs,iprint,l1st)
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
         ,acon(nismax),tag(namax),sv(3,6),chg(namax)
    real(8),intent(inout):: tcom
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical,intent(in):: lstrs,l1st

    integer:: i,j,k,l,m,n,ierr,is,js,ixyz,jxyz
    real(8):: dij,dedr,epott,x,y,z,epotl,tmp,texp,d0ij,alpij,rminij &
         ,chgi,chgj,tmp2
    real(8),allocatable,save:: strsl(:,:,:)
    type(atdesc):: atdi,atdj
    real(8),external:: fcut1,dfcut1,sprod
    real(8),save,allocatable:: xi(:),xj(:),xij(:),rij(:),diji(:)&
         ,dxdi(:),dxdj(:),at(:)
    if( .not.allocated(xi) ) allocate(xi(3),xj(3),xij(3),rij(3),diji(3),&
         dxdi(3),dxdj(3),at(3) )

    if( l1st ) then
!!$      call init_vcMorse(natm,tag,mpi_md_world)
!!$      call read_params_vcMorse(myid,mpi_md_world,iprint)
      if( allocated(strsl) ) deallocate(strsl)
      allocate(strsl(3,3,namax))
    endif

    epotl= 0d0
    strsl(1:3,1:3,1:namax) = 0d0

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
        d0ij = sprod(wd,pdij)
        alpij= sprod(walp,pdij)
        rminij= sprod(wrmin,pdij)
        texp = exp(alpij*(rminij-dij))
!.....potential
        tmp= d0ij*((texp-1d0)**2 -1d0)
        tmp2 = 0.5d0 *tmp *fcut1(dij,rc)
        if( j.le.natm ) then
          epi(i)= epi(i) +tmp2
          epi(j)= epi(j) +tmp2
          epotl = epotl +tmp2 +tmp2
        else
          epi(i)= epi(i) +tmp2
          epotl = epotl +tmp2
        endif
!.....force
        dedr= -2d0 *alpij *d0ij *texp *(texp-1d0) *fcut1(dij,rc) &
             + tmp *dfcut1(dij,rc)
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
 
  end subroutine force_vcMorse
!=======================================================================
  subroutine qforce_vcMorse(namax,natm,tag,ra,fq,nnmax,chg &
       ,h,tcom,rc,lspr,mpi_md_world,myid,epot,iprint,l1st)
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
    real(8),intent(inout):: tcom
    real(8),intent(out):: epot,fq(namax)
    logical,intent(in):: l1st 

    integer:: i,j,k,l,m,n,ierr,is,js,ixyz,jxyz
    real(8):: dij,dedr,epott,x,y,z,epotl,tmp,texp,d0ij,alpij,rminij &
         ,chgi,chgj,dd0dq,dalpdq,drmindq,dedd0,dedalp,dedrmin,tmp2
    type(atdesc):: atdi,atdj
    real(8),external:: fcut1,sprod
    real(8),save,allocatable:: xi(:),xj(:),xij(:),rij(:),diji(:) &
         ,dxdi(:),dxdj(:),at(:)

    if( .not.allocated(xi) ) allocate(xi(3),xj(3),xij(3),rij(3),diji(3), &
         dxdi(3),dxdj(3),at(3) )

    epotl= 0d0

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
        d0ij = sprod(wd,pdij)
        alpij= sprod(walp,pdij)
        rminij= sprod(wrmin,pdij)
        texp = exp(alpij*(rminij-dij))
!.....potential
        tmp= 0.5d0 * d0ij*((texp-1d0)**2 -1d0)
        tmp2 = tmp *fcut1(dij,rc)
        if( j.le.natm ) then
          epotl = epotl +tmp2 +tmp2
        else
          epotl = epotl +tmp2
        endif
!.....force on charge
        dd0dq = wd(1) +wd(2)*chgj
        dalpdq = walp(1) +walp(2)*chgj
        drmindq = wrmin(1) +wrmin(2)*chgj
        dedd0 = (texp -1d0)**2 -1d0
        dedalp = 2d0*d0ij*(texp-1d0)*texp*(rminij-dij)
        dedrmin = -2d0*d0ij*(texp-1d0)*texp*alpij
        fq(i) = fq(i) +(dd0dq*dedd0 +dalpdq*dedalp +drmindq*dedrmin) &
             *fcut1(dij,rc)
!!$        if( l1st ) then
!!$          write(6,'(a,2i5,30es10.2)') 'i,j,rij,chgi,chgj,d0ij,alpij,rminij'&
!!$               //',epotl,fq(i)=',i,j,rij(1:3),chgi,chgj,d0ij,alpij,rminij &
!!$               ,epotl,fq(i)
!!$        endif
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
  subroutine init_Morse(natm,tag,mpi_md_world)
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
    
  end subroutine init_Morse
!=======================================================================
  subroutine init_vcMorse(natm,tag,mpi_md_world)
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

!!$!.....Allocate parameter arrays
!!$    if( .not.allocated(walp) ) then
!!$      allocate(walp(0:nprm),wd(0:nprm),wrmin(0:nprm),pdij(0:nprm))
!!$    endif

  end subroutine init_vcMorse
!=======================================================================
  subroutine read_params_Morse(myid_md,mpi_md_world,iprint)
!
!  Read pair parameters for Morse potential from file
!
    include 'mpif.h'
    integer,intent(in):: myid_md,mpi_md_world,iprint
    integer:: i,j,isp,jsp,id,ierr
    character(len=128):: cline,fname
    real(8):: d,r,a

    if( myid_md.eq.0 ) then
      fname = trim(paramsdir)//'/'//trim(paramsfname)
      open(ioprms,file=trim(fname),status='old')
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
        write(6,'(a)') ' finish reading '//trim(fname)
        write(6,*) ''
      endif
    endif

    call mpi_bcast(d0,nsp*nsp,mpi_double_precision,0,mpi_md_world,ierr)
    call mpi_bcast(rmin,nsp*nsp,mpi_double_precision,0,mpi_md_world,ierr)
    call mpi_bcast(alp,nsp*nsp,mpi_double_precision,0,mpi_md_world,ierr)
    call mpi_bcast(interact,nsp*nsp,mpi_logical,0,mpi_md_world,ierr)
    
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

      if( iprint.ne.0 ) then
        print *, 'Finished reading '//trim(fname)
        print *, ''
      endif
    endif

    call mpi_bcast(walp,nprm+1,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(wd,nprm+1,mpi_real8,0,mpi_md_world,ierr)
    call mpi_bcast(wrmin,nprm+1,mpi_real8,0,mpi_md_world,ierr)

    lprmset = .true.
    
  end subroutine read_params_vcMorse
!=======================================================================
  subroutine set_params_vcMorse(ndimp,params,ierr)
!
! Accessor routine to set vcMorse parameters from outside.
! Curretnly this routine is supposed to be called only on serial run.
! So no need of treating this as parallel code.
!
    integer,intent(in):: ndimp
    real(8),intent(in):: params(ndimp)
    integer,intent(inout):: ierr

    integer:: i,inc

    ierr = 0
    if( ndimp.ne.3*(nprm+1) ) then
      ierr = 10
      print *,'Error: ndimp.ne.3*(nprm+1) !!!'
      stop
      return
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
    
    lprmset = .true.
    return
  end subroutine set_params_vcMorse
!=======================================================================
  subroutine read_element_descriptors(myid_md,mpi_md_world,iprint)
!
!  Read descriptors for elements from file.
!
    integer,intent(in):: myid_md,mpi_md_world,iprint

    integer:: itmp,isp
    character(len=5):: ctmp
    character(len=128):: cline,fname
    type(atdesc):: atd
    real(8):: eion1,eion2,eaff,atrad,enpaul

    if( allocated(atdescs) ) deallocate(atdescs)
    allocate(atdescs(nsp))
    
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
          if( isp.gt.nsp ) then
            if( iprint.ne.0 ) then
              write(6,*) trim(fname)//' has more entries than nsp.' &
                   //' So skip reading it.'
              write(6,*) '  NSP = ',nsp
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
      if( iprint.ne.0 ) then
        write(6,*) 'Atomic descriptors were loaded from '//trim(fname)
        write(6,*) 'atomic number, symbol,    IE1,    IE2,    EA,'//&
             '  radius,  EN_Pauling'
        do isp=1,nsp
          atd = atdescs(isp)
          write(6,'(i14,5x,a3,5f8.3)') atd%na, &
               trim(atd%csym), &
               atd%eion1, atd%eion2, atd%eaff, atd%atrad, atd%enpaul
        enddo
        print *,''
      endif
    endif

!.....Share the descriptor information with all nodes
    call bcast_atdescs(nsp,myid_md,mpi_md_world)

  end subroutine read_element_descriptors
!=======================================================================
  subroutine bcast_atdescs(nsp,myid_md,mpi_md_world)
!
!  Broadcast type atdesc
!
    implicit none 
    include 'mpif.h'
    integer,intent(in):: nsp,myid_md,mpi_md_world

    integer,allocatable:: nas(:)
    character(len=3),allocatable:: csyms(:)
    real(8),allocatable:: eion1s(:),eion2s(:),eaffs(:),enpauls(:),atrads(:)
    integer:: isp,ierr

    allocate(eion1s(nsp),eion2s(nsp),eaffs(nsp),enpauls(nsp),&
         atrads(nsp),nas(nsp),csyms(nsp))

    if( myid_md.eq.0 ) then
      do isp=1,nsp
        nas(isp) = atdescs(isp)%na
        csyms(isp) = atdescs(isp)%csym
        eion1s(isp) = atdescs(isp)%eion1
        eion2s(isp) = atdescs(isp)%eion2
        eaffs(isp) = atdescs(isp)%eaff
        atrads(isp) = atdescs(isp)%atrad
        enpauls(isp) = atdescs(isp)%enpaul
      enddo
    endif

    call mpi_bcast(nas,nsp,mpi_integer,0,mpi_md_world,ierr)
    call mpi_bcast(csyms,nsp,mpi_character,0,mpi_md_world,ierr)
    call mpi_bcast(eion1s,nsp,mpi_double_precision,0,mpi_md_world,ierr)
    call mpi_bcast(eion2s,nsp,mpi_double_precision,0,mpi_md_world,ierr)
    call mpi_bcast(eaffs,nsp,mpi_double_precision,0,mpi_md_world,ierr)
    call mpi_bcast(atrads,nsp,mpi_double_precision,0,mpi_md_world,ierr)
    call mpi_bcast(enpauls,nsp,mpi_double_precision,0,mpi_md_world,ierr)
    
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

    dj(1) = chgj
    dj(2) = atdj%eion1
    dj(3) = atdj%eion2
    dj(4) = atdj%eaff
    dj(5) = atdj%atrad
    dj(6) = atdj%enpaul

    pdij(0) = 1d0  ! bias
    do i=1,ndesc
      pdij(2*i-1) = di(i)+dj(i)
      pdij(2*i)   = di(i)*dj(i)
    enddo
    return
    
  end subroutine make_pair_desc
!=======================================================================
  function get_coeff(prms,chgi,chgj,atdi,atdj)
!
!  Determine D which depends on charges of atom-i and atom-j.
!
    real(8),intent(in):: chgi,chgj,prms(0:nprm)
    type(atdesc),intent(in):: atdi,atdj
    real(8):: get_coeff

    integer:: i
    real(8):: di(ndesc),dj(ndesc)

    di(1) = chgi
    di(2) = atdi%eion1
    di(3) = atdi%eion2
    di(4) = atdi%eaff
    di(5) = atdi%atrad
    di(6) = atdi%enpaul

    dj(1) = chgj
    dj(2) = atdj%eion1
    dj(3) = atdj%eion2
    dj(4) = atdj%eaff
    dj(5) = atdj%atrad
    dj(6) = atdj%enpaul

    get_coeff = prms(0)
    do i=1,ndesc
      get_coeff = get_coeff +prms(2*i-1)*(di(i)+dj(i)) &
           +prms(2*i)*( di(i)*dj(i) )
    enddo
    return
    
  end function get_coeff
!=======================================================================
  subroutine pderiv_vcMorse(namax,natm,tag,ra,nnmax,chg &
       ,h,rc,lspr,epot,iprint,ndimp,pderiv)
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
    real(8),intent(out):: epot
    integer,intent(in):: ndimp
    real(8),intent(inout):: pderiv(ndimp)

    integer:: i,j,k,l,m,n,ierr,is,js,ixyz,jxyz,inc
    real(8):: xi(3),xj(3),xij(3),rij(3),dij,diji(3),dedr,epott &
         ,dxdi(3),dxdj(3),x,y,z,epotl,at(3),tmp,texp,d0ij,alpij,rminij &
         ,chgi,chgj,dd0dq,dalpdq,drmindq,dedd0,dedalp,dedrmin,tmp2
    type(atdesc):: atdi,atdj
    real(8),external:: fcut1,sprod

!!$    if( .not.allocated(gwalp) ) allocate(gwalp(0:nprm), &
!!$         gwd(0:nprm),gwrmin(0:nprm))

    epotl= 0d0

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
        diji= 1d0/dij
        dxdi(1:3)= -rij(1:3)*diji
        dxdj(1:3)=  rij(1:3)*diji
        chgj = chg(j)
        atdj = atdescs(js)
        call make_pair_desc(chgi,chgj,atdi,atdj,pdij)
!.....Create Morse parameters that depend on current atom charges
        d0ij = sprod(wd,pdij)
        alpij= sprod(walp,pdij)
        rminij= sprod(wrmin,pdij)
        texp = exp(alpij*(rminij-dij))
!.....potential
        tmp= 0.5d0 * d0ij*((texp-1d0)**2 -1d0)
        tmp2 = tmp *fcut1(dij,rc)
        epotl = epotl +tmp2
!.....Derivative of potential energy w.r.t. {w}
        dedd0 = (texp -1d0)**2 -1d0
        dedalp = 2d0*d0ij*(texp-1d0)*texp*(rminij-dij)
        dedrmin = -2d0*d0ij*(texp-1d0)*texp*alpij
        gwd(0:nprm) = gwd(0:nprm) +0.5d0 *dedd0 *pdij(0:nprm)
        gwalp(0:nprm) = gwalp(0:nprm) +0.5d0 *dedalp *pdij(0:nprm)
        gwrmin(0:nprm) = gwrmin(0:nprm) +0.5d0 *dedrmin *pdij(0:nprm)
      enddo
    enddo
    
    epot= epot +epotl

!.....gwalp,gwd,gwrmin to pderiv
    inc = 0
    do i=0,ndesc*2
      inc = inc + 1
      pderiv(inc) = gwalp(i)
    enddo
    do i=0,ndesc*2
      inc = inc + 1
      pderiv(inc) = gwd(i)
    enddo
    do i=0,ndesc*2
      inc = inc + 1
      pderiv(inc) = gwrmin(i)
    enddo
    return
  end subroutine pderiv_vcMorse
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

end module Morse
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
