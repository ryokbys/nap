subroutine get_force(namax,natm,tag,ra,nnmax,aa,strs,chg,chi,stnsr &
     ,h,hi,tcom,nb,nbmax,lsb,lsex,nex,lsrc,myparity,nnn,sv,rc,lspr &
     ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs &
     ,numff,cffs,ifcoulomb,iprint,l1st, &
     luse_LJ,luse_Ito3_WHe,luse_RK_WHe,luse_RK_FeH,  &
     luse_Ramas_FeH,luse_Ackland_Fe,luse_SW_Si,luse_EDIP_Si, &
     luse_Brenner,luse_Brenner_vdW,luse_Lu_WHe,luse_Branicio_AlN, &
     luse_Mishin_Al,luse_AFS_W,luse_SC_Fe,luse_SM_Al,luse_EAM, &
     luse_linreg,luse_NN,luse_Morse,luse_Morse_repul,luse_vcMorse,lvc,&
     luse_Buckingham,luse_Bonny_WRe,lcell_updated)
!-----------------------------------------------------------------------
!  Wrapper routine for force calculations.
!  Every force calculation routine is called from this subroutine and
!  new force routine should also be implemented in this subroutine.
!-----------------------------------------------------------------------
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
  use EAM,only:force_EAM
  use linreg,only:force_linreg
  use NN,only:force_NN
  use Coulomb, only: force_screened_Coulomb, force_Ewald_Coulomb &
       ,initialize_coulomb, force_long_Coulomb
  use Morse, only: force_Morse, force_Morse_repul, force_vcMorse
  use Buckingham,only:force_Buckingham
  use Bonny_WRe,only: force_Bonny_WRe
  implicit none
  integer,intent(in):: namax,natm,nnmax,nismax,iprint
  integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsex(nbmax,6),lsrc(6) &
       ,myparity(3),nnn(6),mpi_md_world,myid_md,nex(3)
  integer,intent(in):: lspr(0:nnmax,namax),numff
  integer,intent(inout):: ifcoulomb
  real(8),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),sv(3,6) &
       ,acon(nismax),tag(namax)
  real(8),intent(inout):: tcom,rc
  real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax) &
       ,chg(namax),chi(namax),stnsr(3,3)
  character(len=20),intent(in):: cffs(numff)
  logical,intent(in):: l1st,lstrs,lcell_updated
  logical,intent(inout):: luse_LJ,luse_Ito3_WHe,luse_RK_WHe,luse_RK_FeH, &
       luse_Ramas_FeH,luse_Ackland_Fe,luse_SW_Si,luse_EDIP_Si, &
       luse_Brenner,luse_Brenner_vdW,luse_Lu_WHe,luse_Branicio_AlN, &
       luse_Mishin_Al,luse_AFS_W,luse_SC_Fe,luse_SM_Al,luse_EAM, &
       luse_linreg,luse_NN,luse_Morse,luse_Morse_repul,luse_vcMorse, &
       luse_Buckingham,luse_Bonny_WRe
  logical,intent(inout):: lvc

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
         tcom,myid_md,mpi_md_world,iprint,ifcoulomb, &
         luse_vcMorse,l1st)
  endif

!.....Exclusive choice of different Coulomb force-fields
  if( ifcoulomb.eq.1 ) then ! screened Coulomb
    call force_screened_Coulomb(namax,natm,tag,ra,nnmax,aa,strs &
         ,chg,h,hi,tcom &
         ,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
         ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint &
         ,ifcoulomb,l1st)
  else if( ifcoulomb.eq.2 ) then  ! Ewald Coulomb
    call force_Ewald_Coulomb(namax,natm,tag,ra,nnmax,aa,strs &
         ,chg,h,hi,tcom &
         ,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
         ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint &
         ,ifcoulomb,l1st,lcell_updated)
  else if( ifcoulomb.eq.3 ) then ! long-range Coulomb
    call force_long_Coulomb(namax,natm,tag,ra,nnmax,aa,strs &
         ,chg,chi,h,hi,tcom &
         ,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
         ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint &
         ,ifcoulomb,l1st,lcell_updated)
  endif

  if( luse_LJ ) call force_LJ_Ar(namax,natm,tag,ra,nnmax,aa,strs,h &
       ,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
  if( luse_Ito3_WHe ) call force_Ito3_WHe(namax,natm,tag,ra,nnmax,aa &
       ,strs,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
  if( luse_RK_WHe ) call force_RK_WHe(namax,natm,tag,ra,nnmax,aa &
       ,strs,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
  if( luse_RK_FeH ) call force_RK_FeH(namax,natm,tag,ra,nnmax,aa &
       ,strs,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
  if( luse_Ramas_FeH ) call force_Ramas_FeH(namax,natm,tag,ra,nnmax &
       ,aa,strs,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc &
       ,lspr,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
  if( luse_Ackland_Fe ) call force_Ackland_Fe(namax,natm,tag,ra &
       ,nnmax,aa,strs,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn &
       ,sv,rc,lspr,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs &
       ,iprint)
  if( luse_SW_Si ) call force_SW_Si(namax,natm,tag,ra,nnmax,aa,strs &
       ,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
  if( luse_EDIP_Si ) call force_EDIP_Si(namax,natm,tag,ra,nnmax,aa &
       ,strs,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
  if( luse_Brenner ) call force_Brenner(namax,natm,tag,ra,nnmax,aa &
       ,strs,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
  if( luse_Brenner_vdW ) call force_Brenner_vdW(namax,natm,tag,ra &
       ,nnmax,aa,strs,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn &
       ,sv,rc,lspr,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs &
       ,iprint)
  if( luse_Lu_WHe ) call force_Lu_Whe(namax,natm,tag,ra,nnmax,aa &
       ,strs,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
  if( luse_Branicio_AlN ) call force_Branicio_AlN(namax,natm,tag,ra &
       ,nnmax,aa,strs,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn &
       ,sv,rc,lspr,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs &
       ,iprint)
  if( luse_Mishin_Al ) call force_Mishin_Al(namax,natm,tag,ra,nnmax &
       ,aa,strs,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc &
       ,lspr,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
  if( luse_AFS_W ) call force_AFS_W(namax,natm,tag,ra,nnmax,aa,strs &
       ,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
  if( luse_SC_Fe ) call force_SC_Fe(namax,natm,tag,ra,nnmax,aa,strs &
       ,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
  if( luse_SM_Al) call force_SM_Al(namax,natm,tag,ra,nnmax,aa,strs,h &
       ,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
  if( luse_EAM) call force_EAM(namax,natm,tag,ra,nnmax,aa,strs,h &
       ,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint,l1st)
  if( luse_linreg ) call force_linreg(namax,natm,tag,ra,nnmax,aa &
       ,strs,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint)
  if( luse_NN ) call force_NN(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
       ,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint,l1st)
  if( luse_Morse ) call force_Morse(namax,natm,tag,ra,nnmax,aa,strs &
       ,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint,l1st)
  if( luse_Morse_repul ) call force_Morse_repul(namax,natm,tag,ra,nnmax &
       ,aa,strs,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint,l1st)
  if( luse_vcMorse ) call force_vcMorse(namax,natm,tag,ra,nnmax,aa,strs &
       ,chg,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint,l1st)
  if( luse_Buckingham ) call force_Buckingham(namax,natm,tag,ra,nnmax,aa,strs &
       ,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint,l1st)
  if( luse_Bonny_WRe ) call force_Bonny_WRe(namax,natm,tag,ra,nnmax,aa,strs &
       ,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs,iprint,l1st)

!.....convert forces from hmat-coordinates to Cartesian coordinates
  do i=1,natm
    at(1:3)= aa(1:3,i)
    aa(1:3,i)= hi(1:3,1)*at(1) +hi(1:3,2)*at(2) +hi(1:3,3)*at(3)
  enddo

end subroutine get_force
!=======================================================================
subroutine init_force(namax,natm,tag,chg,chi,myid_md,mpi_md_world, &
     iprint,h,rc,numff,cffs, &
     luse_LJ,luse_Ito3_WHe,luse_RK_WHe,luse_RK_FeH,  &
     luse_Ramas_FeH,luse_Ackland_Fe,luse_SW_Si,luse_EDIP_Si, &
     luse_Brenner,luse_Brenner_vdW,luse_Lu_WHe,luse_Branicio_AlN, &
     luse_Mishin_Al,luse_AFS_W,luse_SC_Fe,luse_SM_Al,luse_EAM, &
     luse_linreg,luse_NN,luse_Morse,luse_Morse_repul, &
     luse_vcMorse,lvc,ifcoulomb,luse_Buckingham,luse_Bonny_WRe)
!
!  Initialization routine is separated from main get_force routine,
!  mainly because the initialization is not necessary in case of fitpot.
!
  use Coulomb, only: initialize_coulomb
  use Morse, only: init_vcMorse, read_params_vcMorse, lprmset, &
       read_element_descriptors,init_Morse,read_params_Morse,&
       update_params_Morse
  use EAM, only: init_EAM, read_params_EAM, update_params_EAM, lprmset_EAM
  use NN, only: read_const_NN, read_params_NN, update_params_NN, lprmset_NN
  use Buckingham, only: init_Buckingham, read_params_Buckingham, lprmset_Buckingham
  implicit none
  integer,intent(in):: namax,natm,myid_md,mpi_md_world,iprint,numff
  real(8),intent(in):: tag(namax),h(3,3),rc
  character(len=20),intent(in):: cffs(numff)
  integer,intent(inout):: ifcoulomb
  real(8),intent(inout):: chg(namax),chi(namax)
  logical,intent(inout):: luse_LJ,luse_Ito3_WHe,luse_RK_WHe,luse_RK_FeH, &
       luse_Ramas_FeH,luse_Ackland_Fe,luse_SW_Si,luse_EDIP_Si, &
       luse_Brenner,luse_Brenner_vdW,luse_Lu_WHe,luse_Branicio_AlN, &
       luse_Mishin_Al,luse_AFS_W,luse_SC_Fe,luse_SM_Al,luse_EAM, &
       luse_linreg,luse_NN,luse_Morse,luse_Morse_repul,luse_vcMorse, &
       luse_Buckingham,luse_Bonny_WRe
  logical,intent(inout):: lvc

  call set_force_flags(luse_LJ,luse_Ito3_WHe,luse_RK_WHe, &
       luse_RK_FeH,luse_Ramas_FeH,luse_Ackland_Fe,luse_SW_Si, &
       luse_EDIP_Si,luse_Brenner,luse_Brenner_vdW,luse_Lu_WHe, &
       luse_Branicio_AlN,luse_Mishin_Al,luse_AFS_W,luse_SC_Fe, &
       luse_SM_Al,luse_EAM,luse_linreg,luse_NN,luse_Morse,luse_Morse_repul, &
       luse_vcMorse,luse_Buckingham,luse_Bonny_WRe, &
       ifcoulomb,numff,cffs,myid_md,iprint)

!.....vcMorse requires charge optimization, 
!.....everywhen atomic positions or potential parameters change
  if( luse_vcMorse ) then
    if( ifcoulomb.ne.3 ) then
      if( myid_md.eq.0 .and. iprint.ne.0 ) print *,'ifcoulomb is set to 3,' &
           //' because vcMorse is chosen.'
      ifcoulomb = 3
    endif
    lvc = .true.
  endif

!.....Coulomb interaction
  if( ifcoulomb.eq.1 .or. ifcoulomb.eq.2 .or. ifcoulomb.eq.3 ) then
    call initialize_coulomb(natm,tag,chg,chi,myid_md &
         ,mpi_md_world,ifcoulomb,iprint,h,rc,lvc)
  endif

!.....vcMorse
  if( luse_vcMorse ) then
    call init_vcMorse(natm,tag,mpi_md_world)
    if( .not. lprmset ) then
      call read_params_vcMorse(myid_md,mpi_md_world,iprint)
    endif
    call read_element_descriptors(myid_md,mpi_md_world,iprint)
  endif
!.....Morse
  if( luse_Morse .or. luse_Morse_repul ) then
    call init_Morse(natm,tag,mpi_md_world)
    if( .not.lprmset ) then
      call read_params_Morse(myid_md,mpi_md_world,iprint)
    else
!.....This code is not parallelized, and only for fitpot
      if( ifcoulomb.eq.1 ) then
        call update_params_Morse('BVS')
      else
        call update_params_Morse('full_Morse')
      endif
    endif
  endif
!.....EAM
  if( luse_EAM ) then
    call init_EAM()
    if( .not.lprmset_EAM ) then
      call read_params_EAM(myid_md,mpi_md_world,iprint)
    else
!.....This code is not parallelized, and only for fitpot
      call update_params_EAM()
    endif
  endif
!.....NN
  if( luse_NN ) then
    call read_const_NN(myid_md,mpi_md_world,iprint)
    if( .not.lprmset_NN ) then
      call read_params_NN(myid_md,iprint)
    else
!.....This code is not parallelized, and only for fitpot
      call update_params_NN()
    endif
  endif
!.....Buckingham
  if( luse_Buckingham ) then
    call init_Buckingham()
    if( .not.lprmset_Buckingham ) then
      call read_params_Buckingham(myid_md,mpi_md_world,iprint)
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
subroutine set_force_flags(luse_LJ,luse_Ito3_WHe,luse_RK_WHe, &
     luse_RK_FeH,luse_Ramas_FeH,luse_Ackland_Fe,luse_SW_Si, &
     luse_EDIP_Si,luse_Brenner,luse_Brenner_vdW,luse_Lu_WHe, &
     luse_Branicio_AlN,luse_Mishin_Al,luse_AFS_W,luse_SC_Fe, &
     luse_SM_Al,luse_EAM,luse_linreg,luse_NN,luse_Morse,luse_Morse_repul, &
     luse_vcMorse,luse_Buckingham,luse_Bonny_WRe, &
     ifcoulomb,numff,cffs,myid,iprint)
!     
!     Set flags for forces whether or not using them.
!
  implicit none
  integer,intent(in):: myid,numff,iprint
  integer,intent(inout):: ifcoulomb
  character(len=20),intent(in):: cffs(numff)
  logical,intent(inout):: luse_LJ,luse_Ito3_WHe,luse_RK_WHe &
       ,luse_RK_FeH,luse_Ramas_FeH,luse_Ackland_Fe,luse_SW_Si &
       ,luse_EDIP_Si,luse_Brenner,luse_Brenner_vdW,luse_Lu_WHe &
       ,luse_Branicio_AlN,luse_Mishin_Al,luse_AFS_W,luse_SC_Fe &
       ,luse_SM_Al,luse_EAM,luse_linreg,luse_NN,luse_Morse &
       ,luse_Morse_repul,luse_vcMorse,luse_Buckingham &
       ,luse_Bonny_WRe

  logical,external:: force_on

  luse_LJ = .false.
  luse_Ito3_WHe = .false.
  luse_RK_WHe = .false.
  luse_RK_FeH = .false.
  luse_Ramas_FeH = .false.
  luse_Ackland_Fe = .false.
  luse_SW_Si = .false.
  luse_EDIP_Si = .false.
  luse_Brenner = .false.
  luse_Brenner_vdW = .false.
  luse_Lu_WHe = .false.
  luse_Branicio_AlN = .false.
  luse_Mishin_Al = .false.
  luse_AFS_W = .false.
  luse_SC_Fe = .false.
  luse_SM_Al = .false.
  luse_EAM = .false.
  luse_linreg = .false.
  luse_NN = .false.
  luse_Morse = .false.
  luse_Morse_repul = .false.
  luse_vcMorse = .false.
  luse_Buckingham = .false.
  luse_Bonny_WRe = .false.
  if( force_on('LJ_Ar',numff,cffs) ) luse_LJ = .true.
  if( force_on('Ito3_WHe',numff,cffs) ) luse_Ito3_WHe = .true.
  if( force_on('RK_WHe',numff,cffs) ) luse_RK_WHe = .true.
  if( force_on('RK_FeH',numff,cffs) ) luse_RK_FeH = .true.
  if( force_on('Ramas_FeH',numff,cffs) ) luse_Ramas_FeH = .true.
  if( force_on('Ackland_Fe',numff,cffs) ) luse_Ackland_Fe = .true.
  if( force_on('SW_Si',numff,cffs) ) luse_SW_Si = .true.
  if( force_on('EDIP_Si',numff,cffs) ) luse_EDIP_Si = .true.
  if( force_on('Brenner',numff,cffs) ) luse_Brenner = .true.
  if( force_on('Brenner_vdW',numff,cffs) ) luse_Brenner_vdW = .true.
  if( force_on('Lu_WHe',numff,cffs) ) luse_Lu_WHe = .true.
  if( force_on('Branicio_AlN',numff,cffs) ) luse_Branicio_AlN=.true.
  if( force_on('Mishin_Al',numff,cffs) ) luse_Mishin_Al = .true.
  if( force_on('AFS_W',numff,cffs) ) luse_AFS_W = .true.
  if( force_on('SC_Fe',numff,cffs) ) luse_SC_Fe = .true.
  if( force_on('SM_Al',numff,cffs) ) luse_SM_Al = .true.
  if( force_on('EAM',numff,cffs) ) luse_EAM = .true.
  if( force_on('linreg',numff,cffs) ) luse_linreg = .true.
  if( force_on('NN',numff,cffs) ) luse_NN = .true.
  if( force_on('Morse',numff,cffs) ) luse_Morse = .true.
  if( force_on('Morse_repul',numff,cffs) ) luse_Morse_repul = .true.
  if( force_on('vcMorse',numff,cffs) ) luse_vcMorse = .true.
  if( force_on('Buckingham',numff,cffs) ) luse_Buckingham = .true.
  if( force_on('Bonny_WRe',numff,cffs) ) luse_Bonny_WRe = .true.
!.....Coulomb forces should be exclusive each other
  if( force_on('screened_Coulomb',numff,cffs) ) then
    ifcoulomb = 1
  else if( force_on('Ewald_Coulomb',numff,cffs) ) then
    ifcoulomb = 2
  else if( force_on('long_Coulomb',numff,cffs) ) then
    ifcoulomb = 3
  endif

  if( myid.eq.0 .and. iprint.ne.0 ) then
    write(6,'(/,a)') ' Use the following force-fields:'
    if( luse_LJ ) print *,'  LJ_Ar'
    if( luse_Ito3_WHe ) print *,'  Ito3_WHe'
    if( luse_RK_WHe ) print *,'  RK_WHe'
    if( luse_RK_FeH ) print *,'  RK_FeH'
    if( luse_Ramas_FeH ) print *,'  Ramas_FeH'
    if( luse_Ackland_Fe ) print *,'  Ackland_Fe'
    if( luse_SW_Si ) print *,'  SW_Si'
    if( luse_EDIP_Si ) print *,'  EDIP_Si'
    if( luse_Brenner ) print *,'  Brenner'
    if( luse_Brenner_vdW ) print *,'  Brenner_vdW'
    if( luse_Lu_WHe ) print *,'  Lu_WHe'
    if( luse_Branicio_AlN ) print *,'  Branicio_AlN'
    if( luse_Mishin_Al ) print *,'  Mishin_Al'
    if( luse_AFS_W ) print *,'  AFS_W'
    if( luse_SC_Fe ) print *,'  SC_Fe'
    if( luse_SM_Al ) print *,'  SM_Al'
    if( luse_EAM ) print *,'  EAM'
    if( luse_linreg ) print *,'  linreg'
    if( luse_NN ) print *,'  NN'
    if( luse_Morse ) print *,'  Morse'
    if( luse_Morse_repul ) print *,'  Morse_repul'
    if( luse_vcMorse ) print *,'  vcMorse'
    if( luse_Buckingham ) print *,'  Buckingham'
    if( luse_Bonny_WRe ) print *,'  Bonny_WRe'
!.....Coulomb forces should be exclusive each other
    if( ifcoulomb.eq.1 ) then
      print *,'  screened_Coulomb'
    else if( ifcoulomb.eq.2 ) then
      print *,'  Ewald_Coulomb'
    else if( ifcoulomb.eq.3 ) then
      print *,'  long_Coulomb'
    endif
  endif

end subroutine set_force_flags
!=======================================================================
subroutine dampopt_charge(namax,natm,tag,h,ra,chg,chi,nnmax,lspr,rc, &
     lsb,lsex,nbmax,nb,nnn,myparity,lsrc,nex,&
     tcom,myid,mpi_md_world,iprint,ifcoulomb,luse_vcMorse,l1st)
!
!  Charge optimization/equilibration by damped dynamics.
!  Since not only Coulomb interaction but also other force-fields can
!  depend on atomic charges, this routine is separated out from the
!  force_Coulomb.F90.
!
  use Coulomb, only: qforce_long
  use Morse, only: qforce_vcMorse
  implicit none
  integer,intent(in):: namax,natm,myid,mpi_md_world,iprint,ifcoulomb &
       ,nnmax,lspr(0:nnmax,namax)
  real(8),intent(in):: chi(namax),h(3,3),tag(namax),ra(3,namax),rc
  real(8),intent(inout):: chg(namax),tcom
  logical,intent(in):: luse_vcMorse,l1st
  integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
       ,nnn(6),nex(3),lsex(nbmax,6)

  integer:: istp,i,istp_pos,istp_conv
  real(8):: eclong,epot,epotp,afq,eMorse,fqnorm,vqnorm,dt,alpha,p
  real(8),save,allocatable:: vq(:),fq(:)

  character(len=128),parameter:: algorithm = 'damping'
!!$  character(len=128),parameter:: algorithm = 'FIRE'
  
  integer,parameter:: nstp_dampopt = 1000
  real(8),parameter:: dt_dampopt = 0.001  ! 0.005 fs
  real(8),parameter:: ecrit_dampopt = 1.0d-4
  real(8),parameter:: amassq = 0.002
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

  if( l1st ) then
    if( allocated(vq) ) deallocate(vq,fq)
    allocate(vq(namax),fq(namax))
  endif

!.....Gather forces on charges
  fq(1:natm) = 0d0
  eclong = 0d0
  istp = 0

  call qforce_long(namax,natm,tag,ra,chg,chi,h,tcom,mpi_md_world, &
       myid,iprint,ifcoulomb,fq,eclong)
  epot = eclong

  if( myid.eq.0 .and. iprint.ge.20 ) then
      write(6,'(a)') ' After qforce_long:'
    write(6,'(a,i5,20es11.3)') ' istp,fqs = ',istp,fq(1:min(natm,20))
  endif
  
  if( luse_vcMorse ) then
    eMorse = 0d0
    call qforce_vcMorse(namax,natm,tag,ra,fq,nnmax,chg &
         ,h,tcom,rc,lspr,mpi_md_world,myid,eMorse,iprint,.true.)
    epot = epot + eMorse
!!$    if( myid.eq.0 .and. iprint.ge.20 ) then
!!$      write(6,'(a)') ' After qforce_vcMorse:'
!!$      write(6,'(a,i5,20es11.3)') ' istp,fqs = ',istp,fq(1:min(natm,20))
!!$    endif
  endif

  call suppress_fq(namax,natm,fq,myid,mpi_md_world)
!!$  if( myid.eq.0 .and. iprint.ge.20 ) then
!!$    write(6,'(a)') ' After log(fq):'
!!$    write(6,'(a,i5,20es11.3)') ' istp,fqs = ',istp,fq(1:min(natm,20))
!!$  endif
  call get_average_fq(namax,natm,fq,afq,myid,mpi_md_world)
  if( afq*0d0.ne.0d0 ) return   ! NaN
  fq(1:natm) = fq(1:natm) -afq

  vq(1:natm) = 0d0
  istp_pos = 0
  istp_conv = 0
  alpha = alpha0
  dt = dt_dampopt
  do istp=1,nstp_dampopt
    epotp = epot
    if( myid.eq.0 .and. iprint.ge.20 ) then
      write(6,'(a,i5,20es11.3)') ' istp,chgs = ',istp,chg(1:min(natm,5))&
           ,fq(1:min(natm,5)),vq(1:min(natm,5))
    endif
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
            ,nex,mpi_md_world,chg)

!.....get new qforces
    fq(1:natm) = 0d0
    call qforce_long(namax,natm,tag,ra,chg,chi,h,tcom,mpi_md_world, &
         myid,iprint,ifcoulomb,fq,eclong)
    epot = eclong
    if( myid.eq.0 .and. iprint.ge.20 ) then
      write(6,'(a)') ' After qforce_long:'
      write(6,'(a,i5,20es11.3)') ' istp,fqs = ',istp,fq(1:min(natm,20))
    endif
    if( luse_vcMorse ) then
      eMorse = 0d0
      call qforce_vcMorse(namax,natm,tag,ra,fq,nnmax,chg &
           ,h,tcom,rc,lspr,mpi_md_world,myid,eMorse,iprint,.false.)
      epot = epot + eMorse
      if( myid.eq.0 .and. iprint.ge.20 ) then
        write(6,'(a)') ' After qforce_vcMorse:'
        write(6,'(a,i5,20es11.3)') ' istp,fqs = ',istp,fq(1:min(natm,20))
      endif
    endif

    call suppress_fq(namax,natm,fq,myid,mpi_md_world)
    call get_average_fq(namax,natm,fq,afq,myid,mpi_md_world)
    if( afq*0d0.ne.0d0 ) return   ! NaN
    fq(1:natm) = fq(1:natm) -afq

!.....second update of velocity
    vq(1:natm) = vq(1:natm) +0.5d0*dt/amassq*fq(1:natm)

!.....check convergence
    if( istp.gt.nstp_min .and. &
         abs(epot-epotp).lt.ecrit_dampopt ) then
      istp_conv = istp_conv +1
      if( istp_conv.gt.nstp_conv ) then
        if( myid.eq.0 .and. iprint.ge.10 ) then
          write(6,'(a,i4,a)') ' Dampopt_charge converged with ', &
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
subroutine get_average_fq(namax,natm,fq,afq,myid,mpi_md_world)
!
!  Get an average FQ that is an average chemical potential.
!
  implicit none
  include "mpif.h"
  integer,intent(in):: namax,natm,myid,mpi_md_world
  real(8),intent(in):: fq(namax)
  real(8),intent(out):: afq

  integer:: i,ierr
  integer,save:: ntot
  real(8):: afql
  logical,save:: l1st = .true.

  if( l1st ) then
    ntot = 0
    call mpi_allreduce(natm,ntot,1,mpi_integer, &
         mpi_sum,mpi_md_world,ierr)
  endif

  afql = 0d0
  do i=1,natm
    afql = afql +fq(i)
  enddo
  afq = 0d0
  call mpi_allreduce(afql,afq,1,mpi_real8,mpi_sum,mpi_md_world,ierr)
  afq = afq/ntot
  
  return
end subroutine get_average_fq
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
subroutine get_pderiv(namax,natm,tag,ra,nnmax,aa,strs,chg,chi &
     ,h,hi,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
     ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs &
     ,numff,cffs,ifcoulomb,iprint,l1st,luse_vcMorse,lvc &
     ,ndimp,pderiv)
!
!  Compute derivative of potential energy (and forces) 
!  w.r.t. potential parameters.
!
  use Morse,only: pderiv_vcMorse
  implicit none
  integer,intent(in):: namax,natm,nnmax,nismax,iprint
  integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
       ,nnn(6),mpi_md_world,myid_md,nex(3)
  integer,intent(in):: lspr(0:nnmax,namax),numff
  integer,intent(inout):: ifcoulomb
  real(8),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),sv(3,6) &
       ,acon(nismax),tag(namax)
  real(8),intent(inout):: tcom,rc
  real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax) &
       ,chg(namax),chi(namax)
  character(len=20),intent(in):: cffs(numff)
  logical,intent(in):: l1st
  logical,intent(inout):: luse_vcMorse
  logical,intent(inout):: lvc
  logical,intent(in):: lstrs
  integer,intent(in):: ndimp
  real(8),intent(out):: pderiv(ndimp)

  if( luse_vcMorse ) call pderiv_vcMorse(namax,natm,tag,ra,nnmax,chg &
       ,h,rc,lspr,epot,iprint,ndimp,pderiv)
  
  
end subroutine get_pderiv
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
