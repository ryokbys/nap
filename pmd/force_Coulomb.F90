module Coulomb
!-----------------------------------------------------------------------
!                     Last modified: <2017-04-07 14:01:00 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
!  Parallel implementation of Coulomb potential
!  For screened Coulomb potential:  
!    - Ref. Adams & Rao, Phys. Status Solidi A 208, No.8 (2011)
!    - No cutoff treatment
!  For Ewald Coulomb potential:
!    - ...
!-----------------------------------------------------------------------
  implicit none
  save

  character(len=128),parameter:: paramsfname = 'in.params.Coulomb'
  real(8),parameter:: pi = 3.14159265398979d0

  integer,parameter:: ioprms = 20
!.....Coulomb's constant, acc = 1.0/(4*pi*epsilon0)
  real(8),parameter:: acc  = 14.3998554737d0

  logical,allocatable:: interact(:,:)

  integer:: nsp
!.....ideal valence charges of species
  real(8),allocatable:: vid_bvs(:)
!.....principal quantum numbers of species
  integer,allocatable:: npq_bvs(:)
!.....covalent radius
  real(8),allocatable:: rad_bvs(:)
!.....screening length
  real(8),allocatable:: rho_bvs(:,:)
  real(8),parameter:: fbvs = 0.74d0

!.....charge threshold for Coulomb interaction [default: 0.01]
  real(8),parameter:: qthd = 1d-8

!.....Gaussian width of Ewald sum [default: 1.0 (Angstrom)]
  real(8):: sgm_ew = 1d0
!.....Accuracy controlling parameter for Ewald sum
!.....See, http://www.jncasr.ac.in/ccms/sbs2007/lecturenotes/5day10nov/SBS_Ewald.pdf
!.....Exp(-pacc) = 1e-7 when pacc= 18.0
  real(8),parameter:: pacc   = 18d0
!.....real-space cell volume
  real(8):: avol
!.....k-space variables
  real(8):: b1(3),b2(3),b3(3)
!.....Initial kmax = 20 is hard coded, which has no meaning.
  integer,parameter:: kmaxini = 20
  real(8):: bkmax
  integer:: kmax1,kmax2,kmax3,nk
  logical,allocatable:: lkuse(:,:,:)
!.....kmax threshold
  real(8),parameter:: threshold_kmax = 1d-4

contains
  subroutine force_screened_Coulomb(namax,natm,tag,ra,nnmax,aa,strs &
       ,chg,h,hi,tcom &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid,epi,epot,nismax,acon,lstrs,iprint,ifcoulomb)
    implicit none
    include "mpif.h"
    include "./params_unit.h"
!!$    include "params_BVS_Morse.h"
    integer,intent(in):: namax,natm,nnmax,nismax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),lspr(0:nnmax,namax),nex(3)
    integer,intent(in):: mpi_md_world,myid,ifcoulomb
    real(8),intent(in):: ra(3,namax),h(3,3),hi(3,3),rc &
         ,acon(nismax),tag(namax),sv(3,6)
    real(8),intent(inout):: chg(namax)
    real(8),intent(inout):: tcom
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical:: lstrs

    
    integer:: i,j,k,l,m,n,ierr,is,js,ixyz,jxyz,nconnect(4)
    real(8):: xi(3),xj(3),xij(3),rij(3),dij,diji,dedr &
         ,dxdi(3),dxdj(3),x,y,z,epotl,epott,at(3),tmp &
         ,qi,qj,radi,radj,rhoij,terfc,texp
    real(8),allocatable,save:: strsl(:,:,:)
    logical,save:: l1st=.true.

    if( l1st ) then
      if( myid.eq.0 ) write(6,'(a)') ' use force_screened_Coulomb'
      call initialize(natm,tag,myid,mpi_md_world,ifcoulomb,rc)
      call read_params(myid,mpi_md_world,ifcoulomb,iprint)
      call set_charge_BVS(natm,nb,tag,chg,myid,mpi_md_world)
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
      is= int(tag(i))
      nconnect(i) = 0
      qi= chg(i)
      if( abs(qi).lt.qthd ) cycle
      do k=1,lspr(0,i)
        j=lspr(k,i)
        if(j.eq.0) exit
        if(j.le.i) cycle
        js= int(tag(j))
        if( .not.interact(is,js) ) cycle
        qj= chg(j)
        if( abs(qj).lt.qthd ) cycle
        xj(1:3)= ra(1:3,j)
        xij(1:3)= xj(1:3)-xi(1:3)
        rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij= sqrt(rij(1)**2 +rij(2)**2 +rij(3)**2)
        if( dij.gt.rc ) cycle
        nconnect(i)= nconnect(i) +1
        diji= 1d0/dij
        dxdi(1:3)= -rij(1:3)*diji
        dxdj(1:3)=  rij(1:3)*diji
        rhoij = rho_bvs(is,js)
        terfc = erfc(dij/rhoij)
!.....potential
        tmp= 0.5d0 *acc *qi*qj*diji *terfc
        if( j.le.natm ) then
          epi(i)= epi(i) +tmp
          epi(j)= epi(j) +tmp
          epotl = epotl +tmp +tmp
        else
          epi(i)= epi(i) +tmp
          epotl = epotl +tmp
        endif
!.....force
        texp = exp(-(dij/rhoij)**2)
        dedr= -acc *qi*qj*diji *(1d0*diji*terfc +2d0/rhoij /sqrt(pi) *texp)
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
    call mpi_allreduce(epotl,epott,1,MPI_DOUBLE_PRECISION &
         ,MPI_SUM,mpi_md_world,ierr)
    epot= epot +epott
  end subroutine force_screened_Coulomb
!=======================================================================
  subroutine force_Ewald_Coulomb(namax,natm,tag,ra,nnmax,aa,strs &
       ,chg,h,hi,tcom &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid,epi,epot,nismax,acon,lstrs,iprint,ifcoulomb)
!
!  Coulomb potential and force computation using Ewald sum method.
!  Currently, the efficiency when parallel computation is not taken into account.
!  So it is not appropriate to use this routine for large system.
!
    implicit none
    include "mpif.h"
    include "./params_unit.h"

    integer,intent(in):: namax,natm,nnmax,nismax,iprint
    integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
         ,nn(6),lspr(0:nnmax,namax),nex(3)
    integer,intent(in):: mpi_md_world,myid,ifcoulomb
    real(8),intent(in):: ra(3,namax),h(3,3),hi(3,3),rc &
         ,acon(nismax),tag(namax),sv(3,6)
    real(8),intent(inout):: chg(namax)
    real(8),intent(inout):: tcom
    real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
    logical:: lstrs

    integer:: i,j,ik,is,js,k1,k2,k3,ierr,jj,ixyz,jxyz
    real(8):: elrl,esrl,epot_self,epotl,epott,qi,qj,tmp,ftmp &
         ,bdotr,terfc,diji,dij,ss2i,sgmsq2,rc2,q2tot,q2loc,bb2,sqpi
    real(8),allocatable,save:: qcosl(:),qcos(:),qsinl(:),qsin(:),pflr(:)&
         ,ri(:),bk(:),bk1(:),bk2(:),bk3(:),bb(:),dxdi(:),dxdj(:),rij(:)&
         ,xij(:),xj(:),xi(:)
    real(8),allocatable,save:: strsl(:,:,:)
    logical,save:: l1st= .true.
    real(8),external:: sprod,absv

    if( l1st ) then
      if( myid.eq.0 ) write(6,'(a)') ' Use force_Ewald_Coulomb'
      allocate(ri(3),bk(3),bk1(3),bk2(3),bk3(3),bb(3),dxdi(3),dxdj(3),rij(3)&
         ,xij(3),xj(3),xi(3))
      call initialize(natm,tag,myid,mpi_md_world,ifcoulomb,rc)
      call read_params(myid,mpi_md_world,ifcoulomb,iprint)
      call get_recip_vectors(h)
      if( myid.eq.0 ) then
        write(6,'(a)') ' Reciprocal vectors:'
        write(6,'(a,3es12.3)') '   b1 = ',b1(1:3)
        write(6,'(a,3es12.3)') '   b2 = ',b2(1:3)
        write(6,'(a,3es12.3)') '   b3 = ',b3(1:3)
      endif
!.....kmax# is constant during MD run even if h-matrix can change...
      call setup_kspace()
      if( myid.eq.0 ) then
        write(6,'(a)') ' Number of k-points for Ewald sum:'
        write(6,'(a,i8)') '   kmax1 = ',kmax1
        write(6,'(a,i8)') '   kmax2 = ',kmax2
        write(6,'(a,i8)') '   kmax3 = ',kmax3
        write(6,'(a,i8)') '   total = ',nk
        write(6,*) ''
      endif
      if( .not.allocated(qcos) ) then
        allocate(qcosl(nk),qcos(nk),qsinl(nk),qsin(nk),pflr(nk))
      endif
!.....prefactor for long-range term
      ik = 0
      do k1= -kmax1,kmax1
        bk1(1:3) = k1*b1(1:3)
        do k2= -kmax2,kmax2
          bk2(1:3) = k2*b2(1:3)
          do k3= -kmax3,kmax3
            if( .not. lkuse(k3,k2,k1) ) cycle
            ik= ik+1
            bk3(1:3) = k3*b3(1:3)
            bk(1:3) = bk1(1:3) +bk2(1:3) +bk3(1:3)
            bb2 = absv(bk)
            bb2 = bb2*bb2
            pflr(ik)= 4d0 *pi /bb2 *exp(-0.5d0 *sgm_ew**2 *bb2)
          enddo
        enddo
      enddo
      if( .not. allocated(strsl) ) then
        allocate(strsl(3,3,namax))
      endif
      l1st=.false.
    endif

!.....Compute self term, if fixed charge per atom, it is constant
    q2loc = 0d0
    do i=1,natm
      q2loc = q2loc +chg(i)*chg(i)
    enddo
    call mpi_allreduce(q2loc,q2tot,1,mpi_double_precision &
         ,mpi_sum,mpi_md_world,ierr)
    epot_self = q2tot *acc /sqrt(2d0*pi) /sgm_ew
!!$    write(6,*) ' epot_self = ',epot_self

!.....Compute direct sum
    rc2 = rc*rc
    sgmsq2 = sqrt(2d0)*sgm_ew
    ss2i = 1d0 /sgmsq2
    sqpi = 1d0 /sqrt(pi)
    esrl = 0d0
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is= int(tag(i))
      qi = chg(i)
      if( abs(qi).lt.qthd ) cycle
      do jj=1,lspr(0,i)
        j = lspr(jj,i)
        if( j.eq.0 ) exit
        if( j.le.i ) cycle
        js = int(tag(j))
        qj = chg(j)
        if( abs(qj).lt.qthd ) cycle
        xj(1:3) = ra(1:3,j)
        xij(1:3)= xj(1:3)-xi(1:3)
        rij(1:3)= h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij = rij(1)**2 +rij(2)**2 +rij(3)**2
        if( dij.gt.rc2 ) cycle
        dij = sqrt(dij)
        diji = 1d0/dij
        dxdi(1:3)= -rij(1:3)*diji
        dxdj(1:3)=  rij(1:3)*diji
        terfc = erfc(dij*ss2i)
!.....potential
        tmp = 0.5d0 *acc *qi*qj*diji *terfc
        if( j.le.natm ) then
          epi(i)= epi(i) +tmp
          epi(j)= epi(j) +tmp
          esrl = esrl +tmp +tmp
        else
          epi(i)= epi(i) +tmp
          esrl = esrl +tmp
        endif
!.....force
        ftmp = -acc *qj*qi*diji *( diji *terfc &
             +2d0 *sqpi *ss2i *exp(-(dij*ss2i)**2) )
        aa(1:3,i)= aa(1:3,i) -dxdi(1:3)*ftmp
        aa(1:3,j)= aa(1:3,j) -dxdj(1:3)*ftmp
!.....stress
        if( lstrs ) then
          do ixyz=1,3
            do jxyz=1,3
              strsl(jxyz,ixyz,i)= strsl(jxyz,ixyz,i) &
                   -0.5d0 *ftmp*rij(ixyz)*(-dxdi(jxyz))
              strsl(jxyz,ixyz,j)= strsl(jxyz,ixyz,j) &
                   -0.5d0 *ftmp*rij(ixyz)*(-dxdi(jxyz))
            enddo
          enddo
        endif
      enddo
    enddo
!!$    write(6,*) ' esrl = ',esrl

!.....Compute reciprocal vectors
    call get_recip_vectors(h)
!.....Compute structure factor of the local processor
    qcosl(1:nk)= 0d0
    qsinl(1:nk)= 0d0
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is= int(tag(i))
      qi = chg(i)
      if( abs(qi).lt.qthd ) cycle
      ri(1:3) = h(1:3,1)*xi(1) +h(1:3,2)*xi(2) +h(1:3,3)*xi(3)
      ik = 0
      do k1= -kmax1,kmax1
        bk1(1:3) = k1 *b1(1:3)
        do k2= -kmax2,kmax2
          bk2(1:3) = k2 *b2(1:3)
          do k3= -kmax3,kmax3
            if( .not. lkuse(k3,k2,k1) ) cycle
            ik= ik +1
            bk3(1:3) = k3 *b3(1:3)
            bb(1:3) = bk1(1:3) +bk2(1:3) +bk3(1:3)
            bdotr = sprod(bb,ri)
            qcosl(ik) = qcosl(ik) +qi*cos(bdotr)
            qsinl(ik) = qsinl(ik) +qi*sin(bdotr)
          enddo
        enddo
      enddo
    enddo
!.....allreduce qcos and qsin, which would be stupid and time consuming
    call mpi_allreduce(qcosl,qcos,nk,mpi_double_precision &
         ,mpi_sum,mpi_md_world,ierr)
    call mpi_allreduce(qsinl,qsin,nk,mpi_double_precision &
         ,mpi_sum,mpi_md_world,ierr)
!.....Compute long-range contribution to potential energy
    elrl = 0d0
    ik = 0
    do k1= -kmax1,kmax1
      do k2= -kmax2,kmax2
        do k3= -kmax3,kmax3
          if( .not. lkuse(k3,k2,k1) ) cycle
          ik= ik+1
          elrl= elrl +pflr(ik)*(qcos(ik)*qcos(ik) +qsin(ik)*qsin(ik))
        enddo
      enddo
    enddo
    elrl = elrl /avol *acc /2
!.....Long-range contribution to forces
    do i=1,natm
      xi(1:3)= ra(1:3,i)
      is= int(tag(i))
      qi = chg(i)
      if( abs(qi).lt.qthd ) cycle
      ri(1:3) = h(1:3,1)*xi(1) +h(1:3,2)*xi(2) +h(1:3,3)*xi(3)
      ik = 0
      do k1= -kmax1,kmax1
        bk1(1:3) = k1 *b1(1:3)
        do k2= -kmax2,kmax2
          bk2(1:3) = k2 *b2(1:3)
          do k3= -kmax3,kmax3
            if( .not. lkuse(k3,k2,k1) ) cycle
            ik= ik +1
            bk3(1:3) = k3 *b3(1:3)
            bb(1:3) = bk1(1:3) +bk2(1:3) +bk3(1:3)
            bdotr = sprod(bb,ri)
            aa(1:3,i)= aa(1:3,i) -acc/avol *qi*bb(1:3) *pflr(ik) &
                 *(-sin(bdotr)*qcos(ik) +cos(bdotr)*qsin(ik))
          enddo
        enddo
      enddo
    enddo
    
    
    epotl = elrl +esrl
!.....Gather epot
    call mpi_allreduce(epotl,epott,1,MPI_DOUBLE_PRECISION &
         ,MPI_SUM,mpi_md_world,ierr)
    epot= epot +epott -epot_self
    

  end subroutine force_Ewald_Coulomb
!=======================================================================
  subroutine initialize(natm,tag,myid,mpi_md_world,ifcoulomb,rc)
!
!  Allocate and initialize parameters to be used.
!
    include "mpif.h"
    integer,intent(in):: myid,mpi_md_world,natm,ifcoulomb
    real(8),intent(in):: tag(natm),rc
    integer:: i,nspl,ierr

!.....Get umber of species
    nspl = 0
    do i=1,natm
      nspl = max(int(tag(i)),nspl)
    enddo
    call mpi_allreduce(nspl,nsp,1,mpi_integer,mpi_max &
         ,mpi_md_world,ierr)

    if( ifcoulomb.eq.1 ) then  ! screened Coulomb
      if( .not.allocated(interact) ) allocate(interact(nsp,nsp))
    else if( ifcoulomb.eq.2 ) then  ! Ewald Coulomb
      sgm_ew = rc/sqrt(2d0*pacc)
      bkmax  = 2d0*pacc /rc
      if( myid.eq.0 ) then
        write(6,'(a,es12.3)') '  accuracy parameter p = ', pacc
        write(6,'(a,es12.3)') '  Gaussian width sgm   = ', sgm_ew
        write(6,'(a,es12.3)') '  K-space cutoff       = ', bkmax
      endif
    endif

  end subroutine initialize
!=======================================================================
  subroutine read_params(myid_md,mpi_md_world,ifcoulomb,iprint)
!
!  Read parameters
!
    include "mpif.h"
    integer,intent(in):: myid_md,mpi_md_world,ifcoulomb,iprint
    character(len=128):: cline,c1st
    character(len=5):: cname
    real(8):: vid,rad
    integer:: npq,isp,jsp,ierr,mode
    
    if( myid_md.eq.0 .and. ifcoulomb.eq.1 ) then
      open(ioprms,file=trim(paramsfname),status='old')
      mode = 1
      interact(1:nsp,1:nsp) = .false.
!.....1st line for check the Coulomb computation
      read(ioprms,*) c1st
      if( trim(c1st).ne.'screened_bvs' ) then
        write(6,*) ' Error@read_params: ifcoulomb does not match '&
             //'with Coulomb type: '//trim(c1st)
        stop
      endif
      allocate(rad_bvs(nsp),npq_bvs(nsp),vid_bvs(nsp) &
           ,rho_bvs(nsp,nsp))
      vid_bvs(1:nsp)= 0d0
      do while(.true.)
        read(ioprms,*,end=10) cline
        if( cline(1:1).eq.'#' .or. cline(1:1).eq.'!' ) cycle
        if( trim(cline).eq.'interactions' ) then
          mode = 2
          cycle
        endif
        if( mode.eq.1 ) then
          backspace(ioprms)
          read(ioprms,*) isp, cname, vid, rad, npq
          if( isp.gt.nsp ) then
            write(6,*) ' Warning @read_params: since isp is greater than nsp,'&
                 //' skip reading the line.'
            cycle
          endif
          if( iprint.ne.0 ) then
            write(6,'(a,i3,a5,2f7.3,i4)') ' isp,cname,vid,rad,npq =' &
                 ,isp,trim(cname),vid,rad,npq
          endif
          vid_bvs(isp) = vid
          rad_bvs(isp) = rad
          npq_bvs(isp) = npq
        else if( mode.eq.2 ) then
          backspace(ioprms)
          read(ioprms,*) isp, jsp
          interact(isp,jsp) = .true.
          interact(jsp,isp) = interact(isp,jsp)
        endif
      enddo
10    close(ioprms)
      if( iprint.ne.0 ) then
        write(6,'(a)') ' Finish reading '//trim(paramsfname)
        write(6,*) ''
      endif

      if( c1st.eq.'screened_bvs' ) then
!.....Set screening length
        do isp=1,nsp
          do jsp=1,nsp
            rho_bvs(isp,jsp) = fbvs*(rad_bvs(isp)+rad_bvs(jsp))
          enddo
        enddo
      endif
    endif

    if( ifcoulomb.eq.1 ) then
      call mpi_bcast(vid_bvs,nsp,mpi_double_precision,0,mpi_md_world,ierr)
      call mpi_bcast(rad_bvs,nsp,mpi_double_precision,0,mpi_md_world,ierr)
      call mpi_bcast(npq_bvs,nsp,mpi_integer,0,mpi_md_world,ierr)
      call mpi_bcast(rho_bvs,nsp*nsp,mpi_double_precision &
           ,0,mpi_md_world,ierr)
      call mpi_bcast(interact,nsp*nsp,mpi_logical,0,mpi_md_world,ierr)
    endif
    
  end subroutine read_params
!=======================================================================
  subroutine set_charge_BVS(natm,nb,tag,chg,myid,mpi_md_world)
!
! Reset actual charges of atoms using effective charge information
! in order to keep charge neutrality in the system.
! It is not necessary when performing screened Coulomb, but required
! for non-screened Coulomb such as Ewald sum method. So it may be
! natural to set neutral charges even when using screened Coulomb...
!
! This would be called only once at the beginning.
!
    include "mpif.h"
    integer,intent(in):: natm,nb,myid,mpi_md_world
    real(8),intent(in):: tag(natm+nb)
    real(8),intent(out):: chg(natm+nb)
    
    integer,allocatable:: nbvsl(:),nbvs(:)
    real(8),allocatable:: vc_bvs(:)
    integer:: i,is,ierr
    real(8):: sum_anion,sum_cation

!!$!.....BEGIN DEBUGGING
!!$    do i=1,natm+nb
!!$      is = int(tag(i))
!!$      if( is.eq.1 ) then
!!$        chg(i) = -vid_bvs(is)
!!$      else
!!$        chg(i) = vid_bvs(is)
!!$      endif
!!$      print *,' i,is,chg(i)=',i,is,chg(i)
!!$    enddo
!!$    return
!!$!.....END DEBUGGING

    allocate(nbvsl(nsp),nbvs(nsp),vc_bvs(nsp))
    nbvsl(1:nsp) = 0
    nbvs(1:nsp) = 0
    do i=1,natm
      is = int(tag(i))
      nbvsl(is) = nbvsl(is) +1
    enddo

    call mpi_allreduce(nbvsl,nbvs,nsp,mpi_integer &
         ,mpi_sum,mpi_md_world,ierr)

    sum_anion = vid_bvs(1) *nbvs(1) /sqrt(dble(npq_bvs(1)))

    sum_cation = 0d0
    do is=2,nsp
      sum_cation = sum_cation +vid_bvs(is)*nbvs(is)/sqrt(dble(npq_bvs(is)))
    enddo

!.....Compute valence charges of species
    vc_bvs(1) = vid_bvs(1)/sqrt(dble(npq_bvs(1))) &
         *sqrt(sum_cation/sum_anion)
    do is=2,nsp
      vc_bvs(is)= vid_bvs(is)/sqrt(dble(npq_bvs(is))) &
           *sqrt(sum_anion/sum_cation)
    enddo

    if( myid.eq.0 ) then
      do is=1,nsp
        write(6,'(a,i3,2f7.3)') ' isp, V_ideal, V_actual = '&
             ,is,vid_bvs(is),vc_bvs(is)
      enddo
    endif

    do i=1,natm+nb
      is= int(tag(i))
!.....Negative charge for anion and positive for cation
      if( is .eq. 1 ) then
        chg(i)= -vc_bvs(is)
      else
        chg(i)= vc_bvs(is)
      endif
    enddo
    
    deallocate(nbvsl,nbvs,vc_bvs)
    
  end subroutine set_charge_BVS
!=======================================================================
  subroutine equilibrate_charges(natm,tag,chg,chi,myid,mpi_world)
!
!  Equilbrate atomic charges using atomic electronegativities.
!
    integer,intent(in):: natm,myid,mpi_world
    real(8),intent(in):: tag(natm),chi(natm)
    real(8),intent(out):: chg(natm)

    
    
  end subroutine equilibrate_charges
!=======================================================================
  subroutine get_recip_vectors(h)
    implicit none
    real(8),intent(in):: h(3,3)

    real(8):: a1(3),a2(3),a3(3),a23(3),a12(3),a31(3),pi2
    real(8),external:: sprod,absv

    a1(1:3) = h(1:3,1)
    a2(1:3) = h(1:3,2)
    a3(1:3) = h(1:3,3)
    call vprod(a2,a3,a23)
    call vprod(a3,a1,a31)
    call vprod(a1,a2,a12)
    avol = sprod(a1,absv(a23))
    pi2 = 2d0 *pi
    b1(1:3) = pi2 /avol *a23(1:3)
    b2(1:3) = pi2 /avol *a31(1:3)
    b3(1:3) = pi2 /avol *a12(1:3)
    return
  end subroutine get_recip_vectors
!=======================================================================
  subroutine setup_kspace()
!
!  Estimate kmax# so that exp(-sgm^2*k^2/2)/k^2 outside kmax# is small enough.
!  And set flags to k-points to be used.
!
    real(8),allocatable:: bbs(:,:,:)
    integer:: k1,k2,k3
    real(8):: bk1(3),bk2(3),bk3(3),bk(3),bb,bb2,bbmax
    real(8),external:: absv

    allocate(bbs(-kmaxini:kmaxini,-kmaxini:kmaxini,-kmaxini:kmaxini))

    bbmax = 0d0
    do k1= -kmaxini,kmaxini
      bk1(1:3) = k1*b1(1:3)
      do k2= -kmaxini,kmaxini
        bk2(1:3) = k2*b2(1:3)
        do k3= -kmaxini,kmaxini
          if( k1.eq.0 .and. k2.eq.0 .and. k3.eq.0 ) cycle
          bk3(1:3) = k3*b3(1:3)
          bk(1:3) = bk1(1:3) +bk2(1:3) +bk3(1:3)
          bb2 = absv(bk)
          bbs(k3,k2,k1) = bb2
!!$          bb2= bb2*bb2
!!$          bbs(k3,k2,k1) = exp(-0.5d0 *sgm_ew**2 *bb2)/bb2
!!$          bbmax = max(bbmax,bbs(k3,k2,k1))
        enddo
      enddo
    enddo

    kmax1= 1
    kmax2= 1
    kmax3= 1
    do k1= 1,kmaxini
      do k2= 1,kmaxini
        do k3= 1,kmaxini
!.....threshold factor = 1d-4 is hard coded, which has no meaning
!!$          if( bbs(k3,k2,k1).gt.bbmax*threshold_kmax ) then
          if( bbs(k3,k2,k1).le.bkmax ) then
            kmax1= max(kmax1,k1)
            kmax2= max(kmax2,k2)
            kmax3= max(kmax3,k3)
          endif
        enddo
      enddo
    enddo

    if( .not. allocated(lkuse) ) then
      allocate(lkuse(-kmax3:kmax3,-kmax2:kmax2,-kmax1:kmax1))
    endif

    lkuse(:,:,:) = .false.
    nk = 0
    do k1= -kmax1,kmax1
      do k2= -kmax2,kmax2
        do k3= -kmax3,kmax3
          if( k1.eq.0 .and. k2.eq.0 .and. k3.eq.0 ) cycle
!!$          print *,'k1,k2,k3,bbs(k3,k2,k1),bkmax =',k1,k2,k3,bbs(k3,k2,k1),bkmax
!!$          if( bbs(k3,k2,k1).gt.bbmax*threshold_kmax ) then
          if( bbs(k3,k2,k1).le.bkmax ) then
            lkuse(k3,k2,k1) = .true.
            nk = nk +1
          endif
        enddo
      enddo
    enddo
    
    deallocate(bbs)
  end subroutine setup_kspace
end module Coulomb
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
