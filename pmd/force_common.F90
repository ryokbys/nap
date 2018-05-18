subroutine get_force(namax,natm,tag,ra,nnmax,aa,strs,chg,chi,stnsr &
     ,h,hi,tcom,nb,nbmax,lsb,lsex,nex,lsrc,myparity,nnn,sv,rc,lspr &
     ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs &
     ,ifcoulomb,iprint,l1st &
     ,lvc,lcell_updated,boundary)
!-----------------------------------------------------------------------
!  Wrapper routine for force calculations.
!  Every force calculation routine is called from this subroutine and
!  new force routine should also be implemented in this subroutine.
!-----------------------------------------------------------------------
  use force
  use RK_FeH,only:force_RK_FeH
  use Ramas_FeH,only:force_Ramas_FeH,force_Ackland_Fe
  use RK_WHe,only:force_RK_WHe
  use Ito3_WHe,only:force_Ito3_WHe
  use LJ,only:force_LJ,force_LJ_repul
  use SW,only:force_SW
  use EDIP_Si,only:force_EDIP_Si
  use Brenner,only:force_brenner,force_brenner_vdW
  use Lu_WHe,only:force_Lu_WHe
  use Branicio_AlN,only:force_Branicio_AlN
  use Mishin_Al,only:force_Mishin_Al
  use AFS_W,only:force_AFS_W
  use SC_Fe,only:force_SC_Fe
  use SM_Al,only:force_SM_Al
  use EAM,only:force_EAM
  use linreg,only:force_linreg
  use NN,only:force_NN
  use NN2,only: force_NN2
  use Coulomb, only: force_screened_Coulomb, force_Ewald &
       ,initialize_coulomb, force_Ewald_long, force_Coulomb
  use Morse, only: force_Morse, force_Morse_repul, force_vcMorse
  use Buckingham,only:force_Buckingham
  use Bonny_WRe,only: force_Bonny_WRe
  use ZBL,only: force_ZBL
  implicit none
  integer,intent(in):: namax,natm,nnmax,nismax,iprint
  integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsex(nbmax,6),lsrc(6) &
       ,myparity(3),nnn(6),mpi_md_world,myid_md,nex(3)
  integer,intent(in):: lspr(0:nnmax,namax) !,numff
  integer,intent(inout):: ifcoulomb
  real(8),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),sv(3,6) &
       ,acon(nismax),tag(namax)
  real(8),intent(inout):: tcom,rc
  real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax) &
       ,chg(namax),chi(namax),stnsr(3,3)
!!$    character(len=20),intent(in):: cffs(numff)
  logical,intent(in):: l1st,lstrs,lcell_updated
  logical,intent(inout):: lvc
  character(len=3),intent(in):: boundary

  integer:: ierr,is,i
  real(8):: at(3),tmp

  epot = 0d0
  aa(1:3,1:namax)=0d0
  epi(1:namax)= 0d0
  strs(1:3,1:3,1:namax)= 0d0
  stnsr(1:3,1:3) = 0d0

!.....If varaible charge, optimize charges before any force calc
  if( lvc ) then
    if( l1st .and. myid_md.eq.0 .and. iprint.ne.0 ) then
      write(6,'(/a)') ' Charges are to be equilibrated.'
    endif
    call dampopt_charge(namax,natm,tag,h,ra,chg,chi,nnmax,lspr,rc, &
         lsb,lsex,nbmax,nb,nnn,myparity,lsrc,nex,&
         tcom,myid_md,mpi_md_world,iprint,l1st,boundary)
    if( l1st .and. myid_md.eq.0 .and. iprint.ge.20 ) then
      write(6,'(/a)') ' Charges:'
      tmp = 0d0
      do i=1,natm
        tmp = tmp +chg(i)
        if( i.gt.100 ) cycle
        write(6,'(a,i5,i3,f8.3)') '   i,is,chg(i) = ',i,int(tag(i)),chg(i)
      enddo
      write(6,'(a,f0.3)') ' Total charge = ',tmp
    endif
  endif

!.....Non-exclusive (additive) choice of force-fields
  if( use_force('LJ') ) call force_LJ(namax,natm,tag,ra,nnmax,aa,strs,h &
       ,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
  if( use_force('LJ_repul') ) call force_LJ_repul(namax,natm,tag,ra,nnmax &
       ,aa,strs,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint,l1st)
  if( use_force('Ito3_WHe') ) call force_Ito3_WHe(namax,natm,tag,ra,nnmax,aa &
       ,strs,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
  if( use_force('RK_WHe') ) call force_RK_WHe(namax,natm,tag,ra,nnmax,aa &
       ,strs,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
  if( use_force('RK_FeH') ) call force_RK_FeH(namax,natm,tag,ra,nnmax,aa &
       ,strs,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
  if( use_force('Ramas_FeH') ) call force_Ramas_FeH(namax,natm,tag,ra,nnmax &
       ,aa,strs,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc &
       ,lspr,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
  if( use_force('Ackland_Fe') ) call force_Ackland_Fe(namax,natm,tag,ra &
       ,nnmax,aa,strs,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn &
       ,sv,rc,lspr,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs &
       ,iprint)
  if( use_force('SW') ) call force_SW(namax,natm,tag,ra,nnmax,aa,strs &
       ,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
  if( use_force('EDIP_Si') ) call force_EDIP_Si(namax,natm,tag,ra,nnmax,aa &
       ,strs,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
  if( use_force('Brenner') ) call force_Brenner(namax,natm,tag,ra,nnmax,aa &
       ,strs,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
  if( use_force('Brenner_vdW') ) call force_Brenner_vdW(namax,natm,tag,ra &
       ,nnmax,aa,strs,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn &
       ,sv,rc,lspr,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs &
       ,iprint)
  if( use_force('Lu_WHe') ) call force_Lu_Whe(namax,natm,tag,ra,nnmax,aa &
       ,strs,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
  if( use_force('Branicio_AlN') ) call force_Branicio_AlN(namax,natm,tag,ra &
       ,nnmax,aa,strs,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn &
       ,sv,rc,lspr,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs &
       ,iprint)
  if( use_force('Mishin_Al') ) call force_Mishin_Al(namax,natm,tag,ra,nnmax &
       ,aa,strs,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc &
       ,lspr,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
  if( use_force('AFS_W') ) call force_AFS_W(namax,natm,tag,ra,nnmax,aa,strs &
       ,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
  if( use_force('SC_Fe') ) call force_SC_Fe(namax,natm,tag,ra,nnmax,aa,strs &
       ,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
  if( use_force('SM_Al') ) call force_SM_Al(namax,natm,tag,ra,nnmax,aa,strs,h &
       ,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
  if( use_force('EAM') ) call force_EAM(namax,natm,tag,ra,nnmax,aa,strs,h &
       ,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint,l1st)
  if( use_force('linreg') ) call force_linreg(namax,natm,tag,ra,nnmax,aa &
       ,strs,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint,l1st)
  if( use_force('NN') ) call force_NN(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
       ,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint,l1st)
  if( use_force('NN2') ) call force_NN2(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
       ,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint,l1st)
  if( use_force('Morse') ) call force_Morse(namax,natm,tag,ra,nnmax,aa,strs &
       ,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint,l1st)
  if( use_force('Morse_repul') ) call force_Morse_repul(namax,natm,tag,ra,nnmax &
       ,aa,strs,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint,l1st)
  if( use_force('vcMorse') ) call force_vcMorse(namax,natm,tag,ra,nnmax,aa,strs &
       ,chg,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint,l1st)
  if( use_force('Buckingham') ) call force_Buckingham(namax,natm,tag,ra,nnmax,aa,strs &
       ,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint,l1st)
  if( use_force('Bonny_WRe') ) call force_Bonny_WRe(namax,natm,tag,ra,nnmax,aa,strs &
       ,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint,l1st)
  if( use_force('ZBL') ) call force_ZBL(namax,natm,tag,ra,nnmax,aa,strs &
       ,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint,l1st)

!.....Exclusive choice of different Coulomb force-fields
  if( use_force('screened_Coulomb') ) then ! screened Coulomb
    call force_screened_Coulomb(namax,natm,tag,ra,nnmax,aa,strs &
         ,chg,h,hi,tcom &
         ,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
         ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint &
         ,l1st)
  else if( use_force('Ewald') ) then  ! Ewald Coulomb
    call force_Ewald(namax,natm,tag,ra,nnmax,aa,strs &
         ,chg,chi,h,hi,tcom &
         ,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
         ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint &
         ,l1st,lcell_updated,lvc)
  else if( use_force('Ewald_long') ) then ! long-range Coulomb
    call force_Ewald_long(namax,natm,tag,ra,nnmax,aa,strs &
         ,chg,chi,h,hi,tcom &
         ,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
         ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint &
         ,l1st,lcell_updated,lvc)
  else if( use_force('Coulomb') ) then  ! Coulomb
    call force_Coulomb(namax,natm,tag,ra,nnmax,aa,strs &
         ,chg,chi,h,hi,tcom &
         ,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
         ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint &
         ,l1st,lcell_updated,lvc)
  endif

!.....convert forces from hmat-coordinates to Cartesian coordinates
  do i=1,natm
    at(1:3)= aa(1:3,i)
    aa(1:3,i)= hi(1:3,1)*at(1) +hi(1:3,2)*at(2) +hi(1:3,3)*at(3)
  enddo

end subroutine get_force
!=======================================================================
subroutine init_force(namax,natm,nsp,tag,chg,chi,myid_md,mpi_md_world, &
     iprint,h,rc,lvc,ifcoulomb)
!
!  Initialization routine is separated from main get_force routine.
!
  use force
  use Coulomb, only: initialize_coulomb, initialize_coulombx
  use Morse, only: read_params_vcMorse, lprmset_Morse, &
       read_element_descriptors,read_params_Morse,&
       update_params_Morse
  use EAM, only: init_EAM, read_params_EAM, update_params_EAM, lprmset_EAM
  use NN, only: read_const_NN, read_params_NN, update_params_NN, lprmset_NN
  use Buckingham, only: init_Buckingham, read_params_Buckingham, lprmset_Buckingham
  use ZBL, only: read_params_ZBL
  use LJ, only: read_params_LJ_repul
  use linreg, only: read_params_linreg,lprmset_linreg
  use descriptor, only: read_params_desc,init_desc
  use NN2, only: read_params_NN2,lprmset_NN2
  implicit none
  integer,intent(in):: namax,natm,nsp,myid_md,mpi_md_world,iprint !,numff
  real(8),intent(in):: tag(namax),h(3,3),rc
!!$    character(len=20),intent(in):: cffs(numff)
  integer,intent(inout):: ifcoulomb
  real(8),intent(inout):: chg(namax),chi(namax)
  logical,intent(inout):: lvc

  integer:: i

  if( iprint.ne.0 ) call write_forces(myid_md)

!.....vcMorse requires charge optimization, 
!.....everywhen atomic positions or potential parameters change
  if( use_force('vcMorse') ) then
    if( use_force('Ewald_long') ) then
      if( myid_md.eq.0 .and. iprint.ne.0 ) print *,'Use Ewald_long,' &
           //' because vcMorse is chosen.'
      ifcoulomb = 3
    endif
    lvc = .true.
  endif

!.....Coulomb interaction
  if( use_force('screened_Coulomb') .or. &
       use_force('Ewald') .or. &
       use_force('Ewald_long') ) then
    call initialize_coulomb(natm,nsp,tag,chg,chi,myid_md &
         ,mpi_md_world,ifcoulomb,iprint,h,rc,lvc)
  else if( use_force('Coulomb') ) then
    call initialize_coulombx(natm,nsp,tag,chg,chi,myid_md &
         ,mpi_md_world,ifcoulomb,iprint,h,rc,lvc)
  endif

!.....vcMorse
  if( use_force('vcMorse') ) then
    if( .not. lprmset_Morse ) then
      call read_params_vcMorse(myid_md,mpi_md_world,iprint)
    endif
    call read_element_descriptors(myid_md,mpi_md_world,iprint)
  endif
!.....Morse
  if( use_force('Morse') .or. use_force('Morse_repul') ) then
    if( .not.lprmset_Morse ) then
      call read_params_Morse(myid_md,mpi_md_world,iprint)
!!$    else
!!$!.....This code is not parallelized, and only for fitpot
!!$      if( ifcoulomb.eq.1 ) then
!!$        call update_params_Morse('BVS')
!!$      else
!!$        call update_params_Morse('full_Morse')
!!$      endif
    endif
  endif
!.....EAM
  if( use_force('EAM') ) then
    call init_EAM()
    if( .not.lprmset_EAM ) then
      call read_params_EAM(myid_md,mpi_md_world,iprint)
    else
!.....This code is not parallelized, and only for fitpot
      call update_params_EAM()
    endif
  endif
!.....NN
  if( use_force('NN') ) then
    call read_const_NN(myid_md,mpi_md_world,iprint)
    if( .not.lprmset_NN ) then
      call read_params_NN(myid_md,iprint)
    else
!.....This code is not parallelized, and only for fitpot
      call update_params_NN()
    endif
  else if( use_force('NN2') ) then
    call init_desc()
    if( .not.lprmset_NN2 ) then
!.....Read both in.params.desc and in.params.linreg
      call read_params_desc(myid_md,mpi_md_world,iprint)
      call read_params_NN2(myid_md,mpi_md_world,iprint)
    else
!.....Read only in.params.desc
      call read_params_desc(myid_md,mpi_md_world,iprint)
    endif
  endif
!.....Buckingham
  if( use_force('Buckingham') ) then
    call init_Buckingham()
    if( .not.lprmset_Buckingham ) then
      call read_params_Buckingham(myid_md,mpi_md_world,iprint)
    endif
  endif
!.....ZBL
  if( use_force('ZBL') ) then
    call read_params_ZBL(myid_md,mpi_md_world,iprint)
  endif
!.....LJ_repul
  if( use_force('LJ_repul') ) then
    call read_params_LJ_repul(myid_md,mpi_md_world,iprint)
  endif
!.....Linear regression
  if( use_force('linreg') ) then
    call init_desc()
    if( .not.lprmset_linreg ) then
!.....Read both in.params.desc and in.params.linreg
      call read_params_desc(myid_md,mpi_md_world,iprint)
      call read_params_linreg(myid_md,mpi_md_world,iprint)
    else
!.....Read only in.params.desc
      call read_params_desc(myid_md,mpi_md_world,iprint)
    endif
  endif

end subroutine init_force
!=======================================================================
subroutine copy_rho_ba(tcom,namax,natm,nb,nbmax,lsb &
     ,lsrc,myparity,nn,sv,mpi_md_world,rho)
!-----------------------------------------------------------------------
!     Exchanges boundary-atom data among neighbor nodes
!-----------------------------------------------------------------------
  implicit none
  include "mpif.h"
  integer:: status(MPI_STATUS_SIZE)
!-----in
  integer,intent(in):: namax,natm,nb,nbmax,mpi_md_world
  integer,intent(in):: lsb(0:nbmax,6),lsrc(6),myparity(3),nn(6)
  real(8),intent(in):: sv(3,6)
!-----out
  real(8),intent(inout):: rho(natm+nb),tcom

!-----locals
  integer:: i,j,k,l,m,n,kd,kdd,ku,inode,nsd,nsd3,nrc,nrc3,nbnew,ierr
  real(8):: tcom1,tcom2
  logical,save:: l1st=.true.
  real(8),allocatable,save:: dbuf(:),dbufr(:)

  if( l1st ) then
    allocate(dbuf(nbmax),dbufr(nbmax))
    l1st=.false.
  endif

  nbnew= 0

!-----loop over z, y, & x directions
  do kd=1,3
    tcom1= mpi_wtime()
    do kdd=-1,0
      ku= 2*kd +kdd
      inode= nn(ku)
!---------num. of to-be-sent particles
      nsd= lsb(0,ku)
!---------num. of to-be-recieved particles
      nrc= lsrc(ku)

!---------exchange x
      do i=1,nsd
        j=lsb(i,ku)
        dbuf(i)= rho(j)
      enddo
      call mespasd(inode,myparity(kd),dbuf,dbufr,nsd,nrc,21 &
           ,mpi_md_world)
      do i=1,nrc
        rho(natm+nbnew+i)= dbufr(i)
      enddo

!---------mpi barrier
      call mpi_barrier(mpi_md_world,ierr)
!---------accumulate num. of boundary particles
!          write(6,'(a,2i8)') "nbnew,nrc=",nbnew,nrc
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
!=======================================================================
subroutine copy_strs_ba(tcom,namax,natm,nb,nbmax,lsb &
     ,lsrc,myparity,nn,sv,mpi_md_world,strs)
!-----------------------------------------------------------------------
!  Exchanges boundary-atom data among neighbor nodes
!-----------------------------------------------------------------------
  implicit none
  include "mpif.h"
  integer:: status(MPI_STATUS_SIZE)
!-----in
  integer,intent(in):: namax,natm,nb,nbmax,mpi_md_world
  integer,intent(in):: lsb(0:nbmax,6),lsrc(6),myparity(3),nn(6)
  real(8),intent(in):: sv(3,6)
!-----out
  real(8),intent(inout):: strs(9,natm+nb),tcom

!-----locals
  integer:: i,j,k,l,m,n,kd,kdd,ku,inode,nsd,nrc,nbnew,ierr
  real(8):: tcom1,tcom2

  logical,save:: l1st=.true.
  real(8),save,allocatable:: dbuf(:,:),dbufr(:,:)

  if( l1st ) then
    allocate(dbuf(9,nbmax),dbufr(9,nbmax))
    l1st=.false.
  endif

  nbnew= 0

!-----loop over z, y, & x directions
  do kd=1,3
    tcom1= mpi_wtime()
    do kdd=-1,0
      ku= 2*kd +kdd
      inode= nn(ku)
!---------num. of to-be-sent particles
      nsd= lsb(0,ku)
!---------num. of to-be-recieved particles
      nrc= lsrc(ku)

!---------exchange strs
      do i=1,nsd
        j=lsb(i,ku)
        dbuf(1:9,i)= strs(1:9,j)
      enddo
      call mespasd(inode,myparity(kd),dbuf,dbufr,9*nsd,9*nrc,21&
           ,mpi_md_world)
      do i=1,nrc
        strs(1:9,natm+nbnew+i)= dbufr(1:9,i)
      enddo

!---------mpi barrier
      call mpi_barrier(mpi_md_world,ierr)
!---------accumulate num. of boundary particles
!          write(6,'(a,2i8)') "nbnew,nrc=",nbnew,nrc
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
!=======================================================================
subroutine copy_dba_fwd(tcom,namax,natm,nb,nbmax,lsb,nex &
     ,lsrc,myparity,nn,sv,mpi_md_world,x,ndim)
!-----------------------------------------------------------------------
!     Exchanges boundary-atom data among neighbor nodes
!-----------------------------------------------------------------------
  implicit none
  include "mpif.h"
  integer:: status(MPI_STATUS_SIZE)
!-----in
  integer,intent(in):: namax,natm,nb,nbmax,mpi_md_world,ndim
  integer,intent(in):: lsb(0:nbmax,6),lsrc(6),myparity(3),nn(6) &
       ,nex(3)
  real(8),intent(in):: sv(3,6)
!-----out
  real(8),intent(inout):: x(ndim,namax),tcom

!-----locals
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

!-----loop over z, y, & x directions
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
!---------num. of to-be-sent particles
        nsd= lsb(0,ku)
!---------num. of to-be-recieved particles
        nrc= lsrc(ku)

!---------exchange x
        do i=1,nsd
          j=lsb(i,ku)
          dbuf(1:ndim,i)= x(1:ndim,j)
        enddo
        call mespasd(inode,myparity(kd),dbuf,dbufr,ndim*nsd &
             ,ndim*nrc,21,mpi_md_world)
        do i=1,nrc
          x(1:ndim,natm+nbnew+i)= dbufr(1:ndim,i)
        enddo
!---------mpi barrier
        call mpi_barrier(mpi_md_world,ierr)
!---------accumulate num. of boundary particles
!          write(6,'(a,2i8)') "nbnew,nrc=",nbnew,nrc
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
!=======================================================================
subroutine copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb,nex &
     ,lsrc,myparity,nn,mpi_md_world,x,ndim)
!-----------------------------------------------------------------------
!     Send-back & receive reaction on cached-atoms
!-----------------------------------------------------------------------
  implicit none
  include "mpif.h"
  integer,intent(in):: namax,natm,nbmax,nb,mpi_md_world,ndim
  integer,intent(in):: lsb(0:nbmax,6),lsrc(6),myparity(3),nn(6) &
       ,nex(3)
  real(8),intent(inout):: x(ndim,namax),tcom

  integer:: status(MPI_STATUS_SIZE)
  integer:: i,j,k,l,m,n,kd,kdd,ku,kuc,ibkwd,nsd,nsd3,nrc,nrc3,nsdbk &
       ,ierr,natmx
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

!-----natmx
  natmx= natm +nb

!-----num. of sent-back reactions
  nsdbk= 0

!-----send-back reactions in the reverse orer, z, y, & x
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
!---------accumulate num. of already sent-back-particles
        nsdbk=nsdbk +nsd
      enddo
    else
!-------To calculate the communication time
      tcom1=MPI_WTIME()

!-------higher & lower directions
      do kdd=0,-1,-1
        ku= 2*kd +kdd
        if(mod(ku,2).eq.0) then
          kuc= ku-1
        else
          kuc= ku+1
        endif
        ibkwd= nn(kuc)
!---------num. of to-be-sent particles
        nsd= lsrc(ku)
!          nsd3= ndim*nsd
        nsd3= mdim*nsd
!---------num. of to-be-recieved particles
        nrc= lsb(0,ku)
!          nrc3= ndim*nrc
        nrc3= mdim*nrc

!---------to-be-sent-back particles
        do i=1,nsd
          dbuf(1:ndim,i)= x(1:ndim,natmx-nsdbk-nsd+i)
        enddo
        call mespasd(ibkwd,myparity(kd),dbuf,dbufr,nsd3,nrc3,500 &
             ,mpi_md_world)
        do k=1,nrc
          i=lsb(k,ku)
          x(1:ndim,i)= x(1:ndim,i) +dbufr(1:ndim,k)
        enddo

!---------mpi barrier
        call mpi_barrier(mpi_md_world,ierr)
!---------accumulate num. of already sent-back-particles
        nsdbk=nsdbk +nsd
      enddo

!-------Add the communication time to COMT
      tcom2=MPI_WTIME()
      tcom=tcom+tcom2-tcom1
    endif

  enddo

!-----check
  if(nsdbk.ne.nb) then
    write(6,'(a,2i8)') "nsdbk,nb=",nsdbk,nb
    stop "error: nsdbk.ne.nb!!"
  endif

!      deallocate(dbuf,dbufr)
end subroutine copy_dba_bk
!=======================================================================
subroutine reduce_dba_bk(natm,namax,tag,x,ndim)
!-----------------------------------------------------------------------
!  Send-back or reduce reaction on cached-atoms.
!  This routine works only on small MD not on parallel version.
!-----------------------------------------------------------------------
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
!=======================================================================
subroutine distribute_dba(natm,namax,tag,x,ndim)
!-----------------------------------------------------------------------
!  Distribute some values to the cached (boundary) atoms.
!  This routine works only on small MD not on parallel version.
!-----------------------------------------------------------------------
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
!=======================================================================
function hvsd(x)
!
!  Heaviside's stepwise function
!
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
!=======================================================================
function fcut1(r,rc)
!
!     Cutof function type-1
!     f(r) = 1                                  for r <= rs
!          = 1/2 *[1+cos(pi*(r-rs)/(rc-rs))]    for rs < r <= rc
!          = 0                                  for rc < r
!
  real(8),intent(in):: r,rc
  real(8),parameter:: pi = 3.14159265358979d0
  real(8),parameter:: rsr= 0.9d0  ! rs ratio to rc
  real(8):: fcut1
  real(8):: rs

  rs = rc*rsr
  if( r.le.rs ) then
    fcut1 = 1d0
  else if( rs.lt.r .and. r.le.rc ) then
    fcut1 = 0.5d0 *(1d0 +cos(pi*(r-rs)/(rc-rs)))
  else
    fcut1 = 0d0
  endif
  return
end function fcut1
!=======================================================================
function dfcut1(r,rc)
!
!     Derivative of the cutoff function type-1
!
  real(8),intent(in):: r,rc
  real(8),parameter:: pi = 3.14159265358979d0
  real(8),parameter:: rsr= 0.9d0  ! rs ratio to rc
  real(8):: dfcut1,rs

  rs = rc *rsr
  if( r.le.rs ) then
    dfcut1 = 0d0
  else if( rs.lt.r .and. r.le.rc ) then
    dfcut1 = -0.5d0 *pi/(rc-rs) *sin(pi*(r-rs)/(rc-rs))
  else
    dfcut1 = 0d0
  endif
  return
end function dfcut1
!=======================================================================
function fcut2(r,rc)
!
!     Cutoff function type-2
!     f(r) = (1-(r/rc)^2)^3      for r <= rc
!          = 0d0                 for rc < r
!
  real(8),intent(in):: r,rc
  real(8):: fcut2

  if( r.le.rc ) then
    fcut2 = (1d0 -(r/rc)**2)**3
  else
    fcut2 = 0d0
  endif
end function fcut2
!=======================================================================
function dfcut2(r,rc)
!
!     Derivative of the cutoff function type-2
!     f(r) = (1-(r/rc)^2)^3      for r <= rc
!          = 0d0                 for rc < r
!
  real(8),intent(in):: r,rc
  real(8):: dfcut2
  real(8):: rrc

  if( r.le.rc ) then
    rrc = r/rc
    dfcut2 = 3d0 *(1d0-rrc**2)**2  *(-2d0*rrc) /rc
  else
    dfcut2 = 0d0
  endif
end function dfcut2
!=======================================================================
function force_on(force_name,numff,cffs)
  implicit none
  character(len=*),intent(in):: force_name
  integer,intent(in):: numff
  character(len=20),intent(in):: cffs(numff)
  logical:: force_on
  integer:: i

  force_on = .false.
  do i=1,numff
    if( trim(force_name).eq.trim(cffs(i)) ) then
      force_on = .true.
      exit
    endif
  enddo
  return

end function force_on
!=======================================================================
subroutine dampopt_charge(namax,natm,tag,h,ra,chg,chi,nnmax,lspr,rc, &
     lsb,lsex,nbmax,nb,nnn,myparity,lsrc,nex,&
     tcom,myid,mpi_md_world,iprint,l1st,boundary)
!
!  Charge optimization/equilibration by damped dynamics.
!  Since not only Coulomb interaction but also other force-fields can
!  depend on atomic charges, this routine is separated out from the
!  force_Coulomb.F90.
!
  use force
  use Coulomb, only: qforce_long,qforce_short,qforce_self,qlower,qupper&
       ,cterms,avmu,conv_eps
  use Morse, only: qforce_vcMorse
  implicit none
  include "mpif.h"
  integer,intent(in):: namax,natm,myid,mpi_md_world,iprint &
       ,nnmax,lspr(0:nnmax,namax)
  real(8),intent(in):: chi(namax),h(3,3),tag(namax),ra(3,namax),rc
  real(8),intent(inout):: chg(namax),tcom
  logical,intent(in):: l1st
  integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
       ,nnn(6),nex(3),lsex(nbmax,6)
  character(len=3),intent(in):: boundary

  integer:: istp,i,istp_pos,istp_conv,is,ntot,ierr
  real(8):: eclong,epot,epotp,eMorse,fqnorm,vqnorm,dt,alpha,p &
       ,ecoul,ecshort,eself
  real(8),save,allocatable:: vq(:),fq(:)
  logical,save,allocatable:: lqfix(:) 

!!$  character(len=128),parameter:: algorithm = 'damping'
  character(len=128),parameter:: algorithm = 'FIRE'

  integer,parameter:: nstp_dampopt = 1000
  real(8),parameter:: dt_dampopt = 0.001  ! 0.005 fs
!!$  real(8),parameter:: ecrit_dampopt = 1.0d-4
  real(8),parameter:: amassq = 0.002
  real(8),parameter:: qtot = 0d0
!.....Simple velocity damping
  real(8),parameter:: fdamp = 0.95
!.....FIRE parameters 
  integer,parameter:: nstp_min = 3
  integer,parameter:: nstp_conv = 3
  real(8),parameter:: finc = 1.1d0
  real(8),parameter:: fdec = 0.5d0
  real(8),parameter:: alpha0 = 0.1d0
  real(8),parameter:: falpha = 0.99d0
  real(8),parameter:: dtmax  = dt_dampopt*10

  integer,external:: count_nonfixed

  if( l1st ) then
    if( allocated(vq) ) deallocate(vq,fq,lqfix)
    allocate(vq(namax),fq(namax),lqfix(namax))
  endif

  call mpi_allreduce(natm,ntot,1,mpi_integer, &
       mpi_sum,mpi_md_world,ierr)
  
!.....Gather forces on charges
  fq(1:namax) = 0d0
  ecshort = 0d0
  eclong = 0d0
  eself = 0d0
  istp = 0

  lqfix(1:namax) = .false.
  call impose_qtot(namax,natm,chg,qtot,lqfix,myid,mpi_md_world)
  do i=1,natm
    is = int(tag(i))
    if( chg(i).gt.qupper(is) ) then
      chg(i) = qupper(is)
      lqfix(i) = .true.
    else if( chg(i).lt.qlower(is) ) then
      chg(i) = qlower(is)
      lqfix(i) = .true.
    endif
  enddo
  if( count_nonfixed(namax,natm,lqfix,myid,mpi_md_world).eq.0 ) then
    if( myid.eq.0 .and. iprint.gt.1 ) &
         print *,'WARNING: exited from dampopt_charge because all Qs are '&
         //'fixed at upper/lower bounds.'
    return
  endif
!!$  if( use_force('Ewald_long') ) then
  if( trim(cterms).eq.'long' ) then
    call qforce_self(namax,natm,tag,chg,chi,fq,eself)
    call qforce_long(namax,natm,tag,ra,chg,h,tcom,mpi_md_world, &
         myid,iprint,fq,eclong)
!!$  else if( use_force('Ewald') ) then
  else if( trim(cterms).eq.'full' ) then
    call qforce_self(namax,natm,tag,chg,chi,fq,eself)
    call qforce_short(namax,natm,tag,ra,nnmax,chg,h,lspr,iprint &
         ,rc,fq,ecshort)
    call qforce_long(namax,natm,tag,ra,chg,h,tcom,mpi_md_world, &
         myid,iprint,fq,eclong)
  else if( trim(cterms).eq.'short' .or. trim(cterms).eq.'screened' ) then
    call qforce_self(namax,natm,tag,chg,chi,fq,eself)
    call qforce_short(namax,natm,tag,ra,nnmax,chg,h,lspr,iprint &
         ,rc,fq,ecshort)
  endif
  epot = eself +ecshort + eclong

  if( use_force('vcMorse') ) then
    eMorse = 0d0
    call qforce_vcMorse(namax,natm,tag,ra,fq,nnmax,chg &
         ,h,tcom,rc,lspr,mpi_md_world,myid,eMorse,iprint,.true.)
    epot = epot + eMorse
  endif
!.....Set fq(i)=0, if lqfix(i)
  do i=1,natm
    if( lqfix(i) ) then
      fq(i) = 0d0
      vq(i) = 0d0
    endif
  enddo
  call suppress_fq(namax,natm,fq,myid,mpi_md_world)
  call get_average_fq(namax,natm,fq,avmu,lqfix,myid,mpi_md_world)
  if( myid.eq.0 .and. iprint.ge.20 ) then
    write(6,'(a,i6,4es12.4)') ' istp,eself,eshort,elong,epot= ',0,eself,ecshort,eclong,epot
    write(6,'(a,i5,20es11.3)') ' istp,avmu,fqs = ',istp,avmu,fq(1:min(natm,10))
  endif

  if( avmu*0d0.ne.0d0 ) then
    if( myid.eq.0 ) print *,'WARNING: exited from dampopt_charge because avmu == NaN'
    return   ! NaN
  endif
  do i=1,natm
    if( lqfix(i) ) cycle
    fq(i) = fq(i) -avmu
  enddo
  fq(1:natm) = fq(1:natm) -avmu

  vq(1:natm) = 0d0
  istp_pos = 0
  istp_conv = 0
  alpha = alpha0
  dt = dt_dampopt
  do istp=1,nstp_dampopt
    epotp = epot
!.....first update of velocity
    vq(1:natm) = vq(1:natm) +0.5d0*dt/amassq*fq(1:natm)

    if( trim(algorithm).eq.'FIRE' ) then
      p = dot_product(fq(1:natm),vq(1:natm))
      vqnorm = sqrt(dot_product(vq(1:natm),vq(1:natm)))
      fqnorm = sqrt(dot_product(fq(1:natm),fq(1:natm)))
      vq(1:natm) = (1d0-alpha)*vq(1:natm) -alpha*vqnorm*fq(1:natm)/fqnorm
      if( p.gt.0d0 ) then
        istp_pos = istp_pos +1
        if( istp_pos.gt.nstp_min ) then  ! start FIRE after nstp_min
          dt = min(dt*finc,dtmax)
          alpha = alpha *falpha
        endif
      else if( p.le.0d0 ) then
        istp_pos = 0
        dt = dt*fdec
        alpha = alpha0
        vq(1:natm) = 0d0
      endif
    else  ! simple velocity damping
      vq(1:natm) = vq(1:natm)*fdamp
    endif
!.....update charges
    chg(1:natm)= chg(1:natm) +vq(1:natm)*dt
    call bacopy_chg_fixed(tcom,lsb,lsex,nbmax,namax &
         ,natm,nb,nnn,myid,myparity,lsrc &
         ,nex,mpi_md_world,chg,boundary)

!.....get new qforces
    eclong = 0d0
    ecshort = 0d0
    eself = 0d0
    fq(1:namax) = 0d0
    call impose_qtot(namax,natm,chg,qtot,lqfix,myid,mpi_md_world)
    do i=1,natm
      if( lqfix(i) ) cycle
      is = int(tag(i))
      if( chg(i).gt.qupper(is) ) then
        chg(i) = qupper(is)
        lqfix(i) = .true.
      else if( chg(i).lt.qlower(is) ) then
        chg(i) = qlower(is)
        lqfix(i) = .true.
      endif
    enddo
    if( count_nonfixed(namax,natm,lqfix,myid,mpi_md_world).eq.0 ) then
      if( myid.eq.0 .and. iprint.gt.1 ) &
           print *,'WARNING: exited from dampopt_charge because all Qs are '&
           //'fixed at upper/lower bounds.'
      return
    endif
!!$    if( use_force('Ewald_long') ) then
    if( trim(cterms).eq.'long' ) then
      call qforce_self(namax,natm,tag,chg,chi,fq,eself)
      call qforce_long(namax,natm,tag,ra,chg,h,tcom,mpi_md_world, &
           myid,iprint,fq,eclong)
!!$    else if( use_force('Ewald') ) then
    else if( trim(cterms).eq.'full' ) then
      call qforce_self(namax,natm,tag,chg,chi,fq,eself)
      call qforce_short(namax,natm,tag,ra,nnmax,chg,h,lspr,iprint &
           ,rc,fq,ecshort)
      call qforce_long(namax,natm,tag,ra,chg,h,tcom,mpi_md_world, &
           myid,iprint,fq,eclong)
    else if( trim(cterms).eq.'short' .or. trim(cterms).eq.'screened' ) then
      call qforce_self(namax,natm,tag,chg,chi,fq,eself)
      call qforce_short(namax,natm,tag,ra,nnmax,chg,h,lspr,iprint &
           ,rc,fq,ecshort)
    endif
    epot = eself +ecshort + eclong
    if( use_force('vcMorse') ) then
      eMorse = 0d0
      call qforce_vcMorse(namax,natm,tag,ra,fq,nnmax,chg &
           ,h,tcom,rc,lspr,mpi_md_world,myid,eMorse,iprint,.false.)
      epot = epot + eMorse
    endif

!.....Set fq(i)=0, if lqfix(i)
    do i=1,natm
      if( lqfix(i) ) then
        fq(i) = 0d0
        vq(i) = 0d0
      endif
    enddo
    call suppress_fq(namax,natm,fq,myid,mpi_md_world)
    call get_average_fq(namax,natm,fq,avmu,lqfix,myid,mpi_md_world)
    if( avmu*0d0.ne.0d0 ) then
      if( myid.eq.0 ) print *,'WARNING: exited from dampopt_charge because avmu == NaN'
      return   ! NaN
    endif
    do i=1,natm
      if( lqfix(i) ) cycle
      fq(i) = fq(i) -avmu
    enddo
    if( myid.eq.0 .and. iprint.ge.20 ) then
      write(6,'(a,i6,4es12.4)') ' istp,eself,eshort,elong,epot= ',istp,eself,ecshort,eclong,epot
      write(6,'(a,i5,2es12.4,20es11.3)') ' istp,epot,de,avmu,fqs= ',istp &
           ,epot,abs(epot-epotp),avmu,fq(1:min(natm,10))
    endif

!.....second update of velocity
    vq(1:natm) = vq(1:natm) +0.5d0*dt/amassq*fq(1:natm)

!.....check convergence
    if( istp.gt.nstp_min .and. &
         abs(epot-epotp)/ntot.lt.conv_eps ) then
      istp_conv = istp_conv +1
      if( istp_conv.gt.nstp_conv ) then
        if( myid.eq.0 .and. iprint.ge.1 ) then
          write(6,'(a,i0,a)') ' Dampopt_charge converged at ', &
               istp,' steps.'
        endif
        exit
      endif
    else
      istp_conv = 0
    endif
  enddo

  return
end subroutine dampopt_charge
!=======================================================================
subroutine get_average_fq(namax,natm,fq,avmu,lqfix,myid,mpi_md_world)
!
!  Get an average FQ that is an average chemical potential.
!
  implicit none
  include "mpif.h"
  integer,intent(in):: namax,natm,myid,mpi_md_world
  real(8),intent(in):: fq(namax)
  logical,intent(in):: lqfix(namax) 
  real(8),intent(out):: avmu

  integer:: i,ierr
  integer:: nonfix,nonfixl
  real(8):: avmul
  logical,save:: l1st = .true.

  nonfixl = 0
  avmul = 0d0
  do i=1,natm
    if( lqfix(i) ) cycle
    nonfixl = nonfixl +1
    avmul = avmul +fq(i)
  enddo
  nonfix = 0
  avmu = 0d0
  call mpi_allreduce(nonfixl,nonfix,1,mpi_integer, &
       mpi_sum,mpi_md_world,ierr)
  call mpi_allreduce(avmul,avmu,1,mpi_real8,mpi_sum,mpi_md_world,ierr)
  avmu = avmu/nonfix

  return
end subroutine get_average_fq
!=======================================================================
subroutine impose_qtot(namax,natm,chg,qtot,lqfix,myid,mpi_md_world)
!
!  Impose total charge to qtot.
!
  implicit none
  include 'mpif.h'
  integer,intent(in):: namax,natm,myid,mpi_md_world
  real(8),intent(in):: qtot
  real(8),intent(inout):: chg(namax)
  logical,intent(in):: lqfix(namax)

  integer:: i,nonfixl,nonfix,ierr
  real(8):: ql,qg,dq,qdist

  ql = 0d0
  nonfixl = 0
  do i=1,natm
    ql = ql +chg(i)
    if( lqfix(i) ) cycle
    nonfixl = nonfixl + 1
  enddo
  nonfix = 0
  call mpi_allreduce(nonfixl,nonfix,1,mpi_integer, &
       mpi_sum,mpi_md_world,ierr)
  if( nonfix.eq.0 ) return  ! no atom for charge to be distributed
  qg = 0d0
  call mpi_allreduce(ql,qg,1,mpi_real8, &
       mpi_sum,mpi_md_world,ierr)
  dq = qg -qtot
  qdist = dq /nonfix
!!$  call mpi_bcast(qdist,1,mpi_real8,0,mpi_md_world,ierr)

  do i=1,natm
    if( lqfix(i) ) cycle
    chg(i) = chg(i) -qdist
  enddo
  
end subroutine impose_qtot
!=======================================================================
function count_nonfixed(namax,natm,lqfix,myid,mpi_md_world)
  implicit none
  include 'mpif.h'
  integer,intent(in):: namax,natm,myid,mpi_md_world
  logical,intent(in):: lqfix(namax)
  integer:: count_nonfixed
  integer:: i,nonfixl,nonfix,ierr

  nonfixl = 0
  do i=1,natm
    if( lqfix(i) ) cycle
    nonfixl = nonfixl + 1
  enddo
  nonfix = 0
  call mpi_allreduce(nonfixl,nonfix,1,mpi_integer, &
       mpi_sum,mpi_md_world,ierr)
  count_nonfixed = nonfix
  return
end function count_nonfixed
!=======================================================================
subroutine suppress_fq(namax,natm,fq,myid,mpi_md_world)
!
! Scale fqs in order to suppress too large fq values
!
! Method == scale:
!   To keep ratio between fqs, all the fqs are multiplied by the same factor
!
  implicit none
  include "mpif.h"
  integer,intent(in):: namax,natm,myid,mpi_md_world
  real(8),intent(inout):: fq(namax)

  integer:: i,ierr
  real(8):: fqmaxl,fqmax,fqfactor
  real(8),parameter:: fqlimit = 10d0
  character(len=128),parameter:: method = 'scale'
!!$  character(len=128),parameter:: method = 'ceiling'

  if( trim(method).eq.'scale' ) then
    fqmaxl = 0d0
    do i=1,natm
      fqmaxl = max(fqmaxl,abs(fq(i)))
    enddo
    fqmax = 0d0
    call mpi_allreduce(fqmaxl,fqmax,1,mpi_real8,mpi_max,mpi_md_world,ierr)

    fqfactor = fqlimit /fqmax
    do i=1,natm
      fq(i) = fq(i)*fqfactor
    enddo
  else if( trim(method).eq.'ceiling' ) then
    do i=1,natm
      fq(i) = sign(1d0,fq(i)) *min(abs(fq(i)),fqlimit)
    enddo
  endif

  return
end subroutine suppress_fq
!=======================================================================
subroutine get_gradw(namax,natm,tag,ra,nnmax,aa,strs,chg,chi &
     ,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
     ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs &
     ,ifcoulomb,iprint,l1st,lvc &
     ,ndimp,gwe,gwf,gws)
!
!  Compute derivative of potential energy (and forces) 
!  w.r.t. potential parameters.
!
  use force
  use Morse,only: gradw_vcMorse, gradw_Morse
  implicit none
  integer,intent(in):: namax,natm,nnmax,nismax,iprint
  integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
       ,nnn(6),mpi_md_world,myid_md,nex(3)
  integer,intent(in):: lspr(0:nnmax,namax) !,numff
  integer,intent(inout):: ifcoulomb
  real(8),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),sv(3,6) &
       ,acon(nismax),tag(namax)
  real(8),intent(inout):: tcom,rc
  real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax) &
       ,chg(namax),chi(namax)
!!$    character(len=20),intent(in):: cffs(numff)
  logical,intent(in):: l1st
  logical,intent(inout):: lvc
  logical,intent(in):: lstrs
  integer,intent(in):: ndimp
  real(8),intent(out):: gwe(ndimp),gwf(ndimp,3,natm),gws(ndimp,6)

  if( use_force('vcMorse') ) call gradw_vcMorse(namax,natm,tag,ra,nnmax,chg &
       ,h,rc,lspr,epot,iprint,ndimp,gwe,gwf,gws)

end subroutine get_gradw
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
