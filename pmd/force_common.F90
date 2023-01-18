subroutine get_force(l1st,epot,stnsr)
!-----------------------------------------------------------------------
!  Wrapper routine for force calculations.
!  Every force calculation routine is called from this subroutine and
!  new force routine should also be implemented in this subroutine.
!-----------------------------------------------------------------------
  use force
  use pmdvars,only: namax,natm,nb,tag,ra,nnmax,aa,strs,aux,naux,nspmax, &
       h,hi,nb,nbmax,lsb,lsex,nex,lsrc,myparity,nn,sv,rc,lspr, &
       sorg,mpi_md_world,myid_md,epi,specorder,lstrs, &
       iprint,lvc,lcell_updated,boundary, &
       iaux_chg, iaux_tei, iaux_q, iaux_vq
  use util,only: iauxof
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
  use Mishin,only:force_Mishin_Al,force_Mishin_Ni
  use AFS_W,only:force_AFS_W
  use SC_Fe,only:force_SC_Fe
  use SM_Al,only:force_SM_Al
  use EAM,only:force_EAM
  use linreg,only:force_linreg
!!$  use NN,only:force_NN
!!$  use NN2,only: force_NN2,force_NN2_overlay_pot, force_NN2_overlay_frc
  use DNN,only: force_DNN
  use Coulomb, only: force_Coulomb, init_for_Ewald, &
       chgopt_method, bacopy_auxq_fixed, qtop, qbot
  use Morse, only: force_Morse, force_Morse_repul, force_vcMorse
  use Buckingham,only:force_Buckingham
  use Bonny_WRe,only: force_Bonny_WRe
  use ZBL,only: force_ZBL,force_ZBL_overlay,r_inner
  use cspline,only: force_cspline
  use tersoff,only: force_tersoff, ts_type
  use Abell,only: force_Abell
  use BMH,only: force_BMH
  use dipole,only: force_dipole
  use fpc,only: force_fpc
  use angular,only: force_angular
  use RFMEAM,only: force_RFMEAM
  use Pellenq,only: force_Pellenq
  use repel,only: force_repel
  use time, only: accum_time
  implicit none
  include "mpif.h"
  include "./const.h"
!!$  integer,intent(in):: namax,natm,nnmax,iprint
!!$  integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsex(nbmax,6),lsrc(6) &
!!$       ,myparity(3),nnn(6),mpi_md_world,myid_md,nex(3)
!!$  integer,intent(in):: lspr(0:nnmax,namax) !,numff
!!$  real(8),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),sv(3,6) &
!!$       ,tag(namax),sorg(3)
!!$  real(8),intent(inout):: rc
!!$  real(8),intent(out):: aa(3,namax),epi(namax),strs(3,3,namax)
!!$  integer,intent(in):: naux
!!$  real(8),intent(inout):: aux(naux,namax)
!!$  logical,intent(in):: lstrs,lcell_updated
!!$  logical,intent(inout):: lvc
!!$  character(len=3),intent(in):: boundary, specorder(nspmax)
  logical,intent(in):: l1st
  real(8),intent(out):: epot,stnsr(3,3)

  integer:: ierr,is,i,ichg,ia
  real(8):: at(3),tmp

  epot = 0d0
  aa(1:3,1:namax)=0d0
  epi(1:natm)= 0d0
  strs(1:3,1:3,1:namax)= 0d0
  stnsr(1:3,1:3) = 0d0

!.....Compute overlay coefficient first
  if( loverlay ) then
    call calc_overlay(namax,natm,nb,nnmax,h,tag,ra,lspr,l1st,iprint)
  endif

!.....init_for_Ewald must be called before chgopt_damping
  if( use_force('Coulomb') ) then
    call init_for_Ewald(h,rc,myid_md,mpi_md_world,iprint,l1st)
  endif

!.....If varaible charge, optimize charges before any charge-dependent potential
  if( lvc ) then
    tmp = mpi_wtime()
    if( l1st .and. myid_md.eq.0 .and. iprint.ge.ipl_basic ) then
      write(6,'(/a)') ' Charges are to be equilibrated.'
    endif
!.....Set initial atomic charges
    if( l1st ) then
      do ia=1,natm+nb
        is = int(tag(ia))
        if( abs(qtop(is)).gt.abs(qbot(is)) ) then
          aux(iaux_chg,ia) = qtop(is)
        else
          aux(iaux_chg,ia) = qbot(is)
        endif
      enddo
    endif
    if( chgopt_method(1:4).eq.'damp' .or. chgopt_method(1:4).eq.'FIRE' ) then
      call chgopt_damping(aux(iaux_chg,:),l1st)
      call accum_time('chgopt_damping',mpi_wtime() -tmp)
    else if( chgopt_method(1:4).eq.'xlag' ) then
      if( l1st ) then
        call chgopt_damping(aux(iaux_chg,:),l1st)
        aux(iaux_q,:) = aux(iaux_chg,:)
      else
        call bacopy_auxq_fixed(aux(iaux_q,:),aux(iaux_vq,:))
        call chgopt_xlag(aux(iaux_chg,:),aux(iaux_q,:),ierr)
        if( ierr.ne.0 ) then
          call chgopt_damping(aux(iaux_chg,:),.false.)
!!$          aux(iaux_q,:) = aux(iaux_chg,:)
        endif
      endif
      call accum_time('chgopt_xlag',mpi_wtime() -tmp)
    else
      if( myid_md.eq.0 ) print *,'Non-available chgopt_method: '//trim(chgopt_method)
      stop 1
    endif
  endif

!.....Non-exclusive (additive) choice of force-fields
  if( use_force('LJ') ) then
    tmp = mpi_wtime()
    call force_LJ(namax,natm,tag,ra,nnmax,aa,strs,h &
       ,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nspmax,lstrs,iprint)
    call accum_time('force_LJ',mpi_wtime() -tmp)
  endif
  if( use_force('LJ_repul') ) then
    tmp = mpi_wtime()
    call force_LJ_repul(namax,natm,tag,ra,nnmax &
       ,aa,strs,h,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nspmax,lstrs,iprint,l1st)
    call accum_time('force_LJ',mpi_wtime() -tmp)
  endif
  if( use_force('Ito3_WHe') ) then
    tmp = mpi_wtime()
    call force_Ito3_WHe(namax,natm,tag,ra,nnmax,aa &
       ,strs,h,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nspmax,lstrs,iprint)
    call accum_time('force_Ito3_WHe',mpi_wtime() -tmp)
  endif
  if( use_force('RK_WHe') ) then
    tmp = mpi_wtime()
    call force_RK_WHe(namax,natm,tag,ra,nnmax,aa &
       ,strs,h,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nspmax,lstrs,iprint)
    call accum_time('force_RK_WHe',mpi_wtime() -tmp)
  endif
  if( use_force('RK_FeH') ) then
    tmp = mpi_wtime()
    call force_RK_FeH(namax,natm,tag,ra,nnmax,aa &
       ,strs,h,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nspmax,lstrs,iprint)
    call accum_time('force_RK_FeH',mpi_wtime() -tmp)
  endif
  if( use_force('Ramas_FeH') ) then
    tmp = mpi_wtime()
    call force_Ramas_FeH(namax,natm,tag,ra,nnmax &
       ,aa,strs,h,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc &
       ,lspr,mpi_md_world,myid_md,epi,epot,nspmax,lstrs,iprint)
    call accum_time('force_Ramas_FeH',mpi_wtime() -tmp)
  endif
  if( use_force('Ackland_Fe') ) then
    tmp = mpi_wtime()
    call force_Ackland_Fe(namax,natm,tag,ra &
       ,nnmax,aa,strs,h,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn &
       ,sv,rc,lspr,mpi_md_world,myid_md,epi,epot,nspmax,lstrs &
       ,iprint)
    call accum_time('force_Ackland_Fe',mpi_wtime() -tmp)
  endif
  if( use_force('SW') ) then
    tmp = mpi_wtime()
    call force_SW(namax,natm,tag,ra,nnmax,aa,strs &
       ,h,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nspmax,specorder,lstrs,iprint)
    call accum_time('force_SW',mpi_wtime() -tmp)
  endif
  if( use_force('EDIP_Si') ) then
    tmp = mpi_wtime()
    call force_EDIP_Si(namax,natm,tag,ra,nnmax,aa &
       ,strs,h,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nspmax,lstrs,iprint)
    call accum_time('force_EDIP',mpi_wtime() -tmp)
  endif
  if( use_force('Tersoff') ) then
    tmp = mpi_wtime()
    call force_tersoff(namax,natm,tag,ra,nnmax,aa &
       ,strs,h,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nspmax,specorder,lstrs,iprint &
       ,aux(iaux_tei,:))
    call accum_time('force_Tersoff',mpi_wtime() -tmp)
  endif
  if( use_force('Brenner') ) then
    tmp = mpi_wtime()
    call force_Brenner(namax,natm,tag,ra,nnmax,aa &
       ,strs,h,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nspmax,lstrs,iprint)
    call accum_time('force_Brenner',mpi_wtime() -tmp)
  endif
  if( use_force('Brenner_vdW') ) then
    tmp = mpi_wtime()
    call force_Brenner_vdW(namax,natm,tag,ra &
       ,nnmax,aa,strs,h,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn &
       ,sv,rc,lspr,mpi_md_world,myid_md,epi,epot,nspmax,lstrs &
       ,iprint)
    call accum_time('force_Brenner_vdW',mpi_wtime() -tmp)
  endif
  if( use_force('Lu_WHe') ) then
    tmp = mpi_wtime()
    call force_Lu_Whe(namax,natm,tag,ra,nnmax,aa &
       ,strs,h,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nspmax,lstrs,iprint)
    call accum_time('force_Lu_WHe',mpi_wtime() -tmp)
  endif
  if( use_force('Branicio_AlN') ) then
    tmp = mpi_wtime()
    call force_Branicio_AlN(namax,natm,tag,ra &
       ,nnmax,aa,strs,h,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn &
       ,sv,rc,lspr,mpi_md_world,myid_md,epi,epot,nspmax,lstrs &
       ,iprint)
    call accum_time('force_Branicio_AlN',mpi_wtime() -tmp)
  endif
  if( use_force('Mishin_Al') ) then
    tmp = mpi_wtime()
    call force_Mishin_Al(namax,natm,tag,ra,nnmax &
       ,aa,strs,h,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc &
       ,lspr,mpi_md_world,myid_md,epi,epot,nspmax,lstrs,iprint,l1st)
    call accum_time('force_Mishin_Al',mpi_wtime() -tmp)
  endif
  if( use_force('Mishin_Ni') ) then
    tmp = mpi_wtime()
    call force_Mishin_Ni(namax,natm,tag,ra,nnmax &
       ,aa,strs,h,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc &
       ,lspr,mpi_md_world,myid_md,epi,epot,nspmax,lstrs,iprint,l1st)
    call accum_time('force_Mishin_Ni',mpi_wtime() -tmp)
  endif
  if( use_force('AFS_W') ) then
    tmp = mpi_wtime()
    call force_AFS_W(namax,natm,tag,ra,nnmax,aa,strs &
       ,h,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nspmax,lstrs,iprint)
    call accum_time('force_AFS_W',mpi_wtime() -tmp)
  endif
  if( use_force('SC_Fe') ) then
    tmp = mpi_wtime()
    call force_SC_Fe(namax,natm,tag,ra,nnmax,aa,strs &
       ,h,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nspmax,lstrs,iprint)
    call accum_time('force_SC_Fe',mpi_wtime() -tmp)
  endif
  if( use_force('SM_Al') ) then
    tmp = mpi_wtime()
    call force_SM_Al(namax,natm,tag,ra,nnmax,aa,strs,h &
       ,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nspmax,lstrs,iprint)
    call accum_time('force_SM_Al',mpi_wtime() -tmp)
  endif
  if( use_force('EAM') ) then
    tmp = mpi_wtime()
    call force_EAM(namax,natm,tag,ra,nnmax,aa,strs,h &
       ,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nspmax,lstrs,iprint,l1st)
    call accum_time('force_EAM',mpi_wtime() -tmp)
  endif
  if( use_force('RFMEAM') ) then
    tmp = mpi_wtime()
    call force_RFMEAM(namax,natm,tag,ra,nnmax,aa,strs,h &
       ,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nspmax,lstrs,iprint,l1st)
    call accum_time('force_RFMEAM',mpi_wtime() -tmp)
  endif
  if( use_force('Pellenq') ) then
    tmp = mpi_wtime()
    call force_Pellenq(namax,natm,tag,ra,nnmax,aa,strs,h &
       ,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nspmax,lstrs,iprint,l1st)
    call accum_time('force_Pellenq',mpi_wtime() -tmp)
  endif
  if( use_force('repel') ) then
    tmp = mpi_wtime()
    call force_repel(namax,natm,tag,ra,nnmax,aa,strs,h &
       ,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nspmax,lstrs,iprint,l1st)
    call accum_time('force_repel',mpi_wtime() -tmp)
  endif
  if( use_force('linreg') ) then
    tmp = mpi_wtime()
    call force_linreg(namax,natm,tag,ra,nnmax,aa &
       ,strs,h,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nspmax,lstrs,iprint,l1st)
    call accum_time('force_linreg',mpi_wtime() -tmp)
  endif
!!$  if( use_force('NN') ) call force_NN(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
!!$       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
!!$       ,mpi_md_world,myid_md,epi,epot,nspmax,lstrs,iprint,l1st)
!!$  if( use_force('NN2') ) then
!!$    if( loverlay ) then
!!$      if( ol_type(1:3).eq.'pot' ) then
!!$        call force_NN2_overlay_pot(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
!!$             ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
!!$             ,mpi_md_world,myid_md,epi,epot,nspmax,lstrs,iprint,l1st)
!!$      else if( ol_type(1:5).eq.'force' ) then
!!$        call force_NN2_overlay_frc(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
!!$             ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
!!$             ,mpi_md_world,myid_md,epi,epot,nspmax,lstrs,iprint,l1st)
!!$      endif
!!$    else
!!$      call force_NN2(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
!!$           ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
!!$           ,mpi_md_world,myid_md,epi,epot,nspmax,lstrs,iprint,l1st)
!!$    endif
!!$  endif
  if( use_force('DNN') ) then
    tmp = mpi_wtime()
    call force_DNN(namax,natm,tag,ra,nnmax,aa,strs,h,hi &
       ,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nspmax,lstrs,iprint,l1st)
    call accum_time('force_DNN',mpi_wtime() -tmp)
  endif
  if( use_force('Morse') ) then
    tmp = mpi_wtime()
    call force_Morse(namax,natm,tag,ra,nnmax,aa,strs &
       ,h,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nspmax,lstrs,iprint,l1st)
    call accum_time('force_Morse',mpi_wtime() -tmp)
  endif
  if( use_force('Morse_repul') ) then
    tmp = mpi_wtime()
    call force_Morse_repul(namax,natm,tag,ra,nnmax &
       ,aa,strs,h,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nspmax,lstrs,iprint,l1st)
    call accum_time('force_Morse_repul',mpi_wtime() -tmp)
  endif
!!$  if( use_force('vcMorse') ) then
!!$    tmp = mpi_wtime()
!!$    call force_vcMorse(namax,natm,tag,ra,nnmax,aa,strs &
!!$       ,aux(iaux_chg,:),h,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
!!$       ,mpi_md_world,myid_md,epi,epot,nspmax,lstrs,iprint,l1st)
!!$    call accum_time('force_vcMorse',mpi_wtime() -tmp)
!!$  endif
  if( use_force('Buckingham') ) then
    tmp = mpi_wtime()
    call force_Buckingham(namax,natm,tag,ra,nnmax,aa,strs &
       ,h,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nspmax,lstrs,iprint,l1st)
    call accum_time('force_Buckingham',mpi_wtime() -tmp)
  endif
  if( use_force('Bonny_WRe') ) then
    tmp = mpi_wtime()
    call force_Bonny_WRe(namax,natm,tag,ra,nnmax,aa,strs &
       ,h,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nspmax,lstrs,iprint,l1st)
    call accum_time('force_Bonny_WRe',mpi_wtime() -tmp)
  endif
  if( use_force('ZBL') ) then
    tmp = mpi_wtime()
    call force_ZBL(namax,natm,tag,ra,nnmax,aa,strs &
         ,h,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
         ,mpi_md_world,myid_md,epi,epot,nspmax,lstrs,iprint,l1st)
    call accum_time('force_ZBL',mpi_wtime() -tmp)
  endif
  if( loverlay .and. trim(ol_force).eq.'ZBL' .and. ol_type(1:3).eq.'pot' ) then
    tmp = mpi_wtime()
    call force_ZBL_overlay(namax,natm,tag,ra,nnmax,aa,strs &
         ,h,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
         ,mpi_md_world,myid_md,epi,epot,nspmax,lstrs,iprint,l1st)
    call accum_time('force_ZBL',mpi_wtime() -tmp)
  endif
  if( use_force('dipole') ) then
    call force_dipole(namax,natm,tag,ra,nnmax,aa,strs &
       ,h,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nspmax,specorder,lstrs,iprint,l1st)
    call accum_time('force_dipole',mpi_wtime() -tmp)
  endif
  if( use_force('cspline') ) then
    tmp = mpi_wtime()
    call force_cspline(namax,natm,tag,ra,nnmax,aa,strs &
       ,h,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nspmax,specorder,lstrs,iprint,l1st)
    call accum_time('force_cspline',mpi_wtime() -tmp)
  endif
  if( use_force('BMH') ) then
    tmp = mpi_wtime()
    call force_BMH(namax,natm,tag,ra,nnmax,aa,strs &
       ,h,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nspmax,specorder,lstrs,iprint,l1st)
    call accum_time('force_BMH',mpi_wtime() -tmp)
  endif
  if( use_force('Abell') ) then
    tmp = mpi_wtime()
    call force_Abell(namax,natm,tag,ra,nnmax,aa,strs &
       ,h,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nspmax,specorder,lstrs,iprint,l1st)
    call accum_time('force_Abell',mpi_wtime() -tmp)
  endif
  if( use_force('fpc') ) then
    tmp = mpi_wtime()
    call force_fpc(namax,natm,tag,ra,nnmax,aa,strs &
       ,h,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nspmax,specorder,lstrs,iprint,l1st)
    call accum_time('force_fpc',mpi_wtime() -tmp)
  endif
  if( use_force('angular') ) then
    tmp = mpi_wtime()
    call force_angular(namax,natm,tag,ra,nnmax,aa,strs &
       ,h,hi,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr &
       ,mpi_md_world,myid_md,epi,epot,nspmax,specorder,lstrs,iprint,l1st)
    call accum_time('force_angular',mpi_wtime() -tmp)
  endif
  
  if( use_force('Coulomb') ) then
    tmp = mpi_wtime()
    if( lvc .and. chgopt_method(1:4).eq.'xlag' .and. .not.l1st ) then
      call force_Coulomb(namax,natm,tag,ra,nnmax,aa,strs &
           ,aux(iaux_q,:),h,hi,nb,nbmax &
           ,lsb,nex,lsrc,myparity,nn,sv,rc,lspr,sorg &
           ,mpi_md_world,myid_md,epi,epot,nspmax,lstrs,iprint &
           ,l1st,lcell_updated,lvc,specorder)
    else
      call force_Coulomb(namax,natm,tag,ra,nnmax,aa,strs &
           ,aux(iaux_chg,:),h,hi,nb,nbmax &
           ,lsb,nex,lsrc,myparity,nn,sv,rc,lspr,sorg &
           ,mpi_md_world,myid_md,epi,epot,nspmax,lstrs,iprint &
           ,l1st,lcell_updated,lvc,specorder)
    endif
    call accum_time('force_Coulomb',mpi_wtime() -tmp)
  endif

!!$!.....convert forces from hmat-coordinates to Cartesian coordinates
!!$  do i=1,natm
!!$    at(1:3)= aa(1:3,i)
!!$    aa(1:3,i)= hi(1:3,1)*at(1) +hi(1:3,2)*at(2) +hi(1:3,3)*at(3)
!!$  enddo
  return
end subroutine get_force
!=======================================================================
subroutine init_force(linit)
!
!  Initialization routine is separated from main get_force routine.
!  And this routine is called from pmd_main not from pmd_core, so h and ra
!  are not yet determined but hmat and rtot are.
!
  use pmdvars,only: namax,nspmax,nsp,myid_md,mpi_md_world,iprint, &
       specorder,rc,lvc,am
  use force
  use Coulomb, only: init_coulomb, lprmset_Coulomb
  use Morse, only: read_params_vcMorse, lprmset_Morse, &
       read_element_descriptors,read_params_Morse,&
       update_params_Morse
  use EAM, only: init_EAM, read_params_EAM, lprmset_EAM
!!$  use NN, only: read_const_NN, read_params_NN, update_params_NN, lprmset_NN
  use Buckingham, only: init_Buckingham, read_params_Buckingham, lprmset_Buckingham
  use ZBL, only: read_params_ZBL, init_ZBL
  use LJ, only: read_params_LJ_repul
  use linreg, only: read_params_linreg,lprmset_linreg
  use descriptor, only: read_params_desc,init_desc,lprmset_desc
  use dspring, only: ldspring
!!$  use NN2, only: read_params_NN2,lprmset_NN2,update_params_NN2
  use DNN, only: read_params_DNN,lprmset_DNN,update_params_DNN
  use tersoff,only: init_tersoff
  use dipole,only: read_params_dipole
  use Abell,only: read_params_Abell, lprmset_Abell
  use BMH,only: read_params_BMH, lprmset_BMH
  use fpc,only: read_params_fpc, lprmset_fpc
  use angular,only: read_params_angular, lprmset_angular
  use RFMEAM, only: read_params_RFMEAM, lprmset_RFMEAM
  use Pellenq,only: read_params_Pellenq, lprmset_Pellenq
  use repel,only: read_params_repel, lprmset_repel
  implicit none
  include "./const.h"
  
  logical,intent(in):: linit

  integer:: i,j
  real(8):: ri,ro
  character(len=3):: cspi,cspj

  if( .not. linit ) return

  if( iprint.ne.0 ) call write_forces(myid_md)

  if( loverlay ) then
    if( myid_md.eq.0 .and. iprint.ge.ipl_basic ) then
      print *,'Overlay parameters of pairs:'
      do i=1,nspmax
        cspi = specorder(i)
        if( trim(cspi).eq.'x' ) cycle
        do j=i,nspmax
          cspj = specorder(j)
          if( trim(cspj).eq.'x' ) cycle
          ri = ol_pair(1,i,j)
          ro = ol_pair(2,i,j)
          print '(3x,a,2f7.3)',trim(cspi)//'-'//trim(cspj),ri,ro
        enddo
      enddo
    endif
  endif

!!$!.....vcMorse requires charge optimization, 
!!$!.....everywhen atomic positions or potential parameters change
!!$  if( use_force('vcMorse') ) then
!!$    if( use_force('Ewald_long') ) then
!!$      if( myid_md.eq.0 .and. iprint.ne.0 ) print *,'Use Ewald_long,' &
!!$           //' because vcMorse is chosen.'
!!$    endif
!!$    lvc = .true.
!!$  endif

!.....Coulomb interaction
  if( use_force('Coulomb') ) then
    if( .not. lprmset_Coulomb ) then
      call init_coulomb(myid_md,mpi_md_world,iprint,lvc,specorder)
    endif
  endif

!!$!.....vcMorse
!!$  if( use_force('vcMorse') ) then
!!$    if( .not. lprmset_Morse ) then
!!$      call read_params_vcMorse(myid_md,mpi_md_world,iprint)
!!$    endif
!!$    call read_element_descriptors(myid_md,mpi_md_world,iprint)
!!$  endif
!.....Morse
  if( use_force('Morse') .or. use_force('Morse_repul') ) then
    if( .not.lprmset_Morse ) then
      if( myid_md.eq.0 .and. iprint.ge.ipl_debug ) print*,'read_params_Morse...'
      call read_params_Morse(myid_md,mpi_md_world,iprint,specorder)
    endif
  endif
!.....dipole
  if( use_force('dipole') ) then
    call read_params_dipole(myid_md,mpi_md_world,iprint &
         ,specorder,am)
  endif
!.....Abell
  if( use_force('Abell') ) then
    if( .not.lprmset_Abell ) then
      call read_params_Abell(myid_md,mpi_md_world,iprint,specorder)
    endif
  endif
!.....fpc
  if( use_force('fpc') ) then
    if( .not.lprmset_fpc ) then
      call read_params_fpc(myid_md,mpi_md_world,iprint,specorder)
    endif
  endif
!.....BMH
  if( use_force('BMH') ) then
    if( .not.lprmset_BMH ) then
      call read_params_BMH(myid_md,mpi_md_world,iprint,specorder)
    endif
  endif
!.....EAM
  if( use_force('EAM') ) then
    call init_EAM()
    if( .not.lprmset_EAM ) then
      call read_params_EAM(myid_md,mpi_md_world,iprint,specorder)
    endif
  endif
!.....RF-MEAM
  if( use_force('RFMEAM') ) then
    if( .not.lprmset_RFMEAM ) then
      call read_params_RFMEAM(myid_md,mpi_md_world,iprint,specorder)
    endif
  endif
!.....Pellenq
  if( use_force('Pellenq') ) then
    if( .not.lprmset_Pellenq ) then
      call read_params_Pellenq(myid_md,mpi_md_world,iprint,specorder)
    endif
  endif
!.....Pellenq
  if( use_force('repel') ) then
    if( .not.lprmset_repel ) then
      call read_params_repel(myid_md,mpi_md_world,iprint,specorder)
    endif
  endif
  
!.....Need to set descriptors before NN or linreg
  if( use_force('DNN') .or. use_force('linreg') .or. ldspring ) then
!.....If descs are already set, no need to read descs from file.
!.....This happens when descs are set from fitpot and re-used for all the samples.
    if( .not.lprmset_desc ) then
      call init_desc()
      call read_params_desc(myid_md,mpi_md_world, &
           iprint,specorder)
    endif
  endif
!!$  if( use_force('NN2') ) then
!!$    if( .not.lprmset_NN2 ) then
!!$!.....Read both in.params.desc and in.params.NN2
!!$      call read_params_NN2(myid_md,mpi_md_world,iprint)
!!$    else
!!$!.....Read only in.params.desc
!!$      call update_params_NN2()
!!$    endif
!!$  endif
!.....Deep NN
  if( use_force('DNN') ) then
    if( .not.lprmset_DNN ) then
!.....Read both in.params.desc and in.params.DNN
      call read_params_DNN(myid_md,mpi_md_world,iprint)
    else
!.....Read only in.params.desc
      call update_params_DNN()
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
    call read_params_ZBL(myid_md,mpi_md_world,iprint,specorder)
  else if( loverlay .and. trim(ol_force).eq.'ZBL' ) then
    call init_ZBL(iprint)
  endif
!.....LJ_repul
  if( use_force('LJ_repul') ) then
    call read_params_LJ_repul(myid_md,mpi_md_world,iprint,specorder)
  endif
!.....Linear regression
  if( use_force('linreg') ) then
    if( .not.lprmset_linreg ) then
!.....Read both in.params.desc and in.params.linreg
      call read_params_linreg(myid_md,mpi_md_world,iprint)
    endif
  endif

!.....Tersoff
  if( use_force('Tersoff') ) then
    call init_tersoff(myid_md,mpi_md_world,iprint,specorder)
  endif
  
!.....angular
  if( use_force('angular') ) then
    if( .not.lprmset_angular ) then
      if( myid_md.eq.0 .and. iprint.ge.ipl_debug ) print*,'read_params_angular...'
      call read_params_angular(myid_md,mpi_md_world,iprint,specorder)
    endif
  endif


end subroutine init_force
!=======================================================================
subroutine copy_rho_ba(namax,natm,nb,nbmax,lsb &
     ,lsrc,myparity,nn,sv,mpi_md_world,rho)
!-----------------------------------------------------------------------
!     Exchanges boundary-atom data among neighbor nodes
!-----------------------------------------------------------------------
  use memory, only: accum_mem
  implicit none
  include "mpif.h"
  integer:: status(MPI_STATUS_SIZE)
!-----in
  integer,intent(in):: namax,natm,nb,nbmax,mpi_md_world
  integer,intent(in):: lsb(0:nbmax,6),lsrc(6),myparity(3),nn(6)
  real(8),intent(in):: sv(3,6)
!-----out
  real(8),intent(inout):: rho(natm+nb)

!-----locals
  integer:: i,j,k,l,m,n,kd,kdd,ku,inode,nsd,nsd3,nrc,nrc3,nbnew,ierr,mem
  logical,save:: l1st=.true.
  real(8),allocatable,save:: dbuf(:),dbufr(:)

  if( l1st ) then
    call accum_mem('force_common',8*nbmax*2)
    allocate(dbuf(nbmax),dbufr(nbmax))
    l1st=.false.
  endif

  if( .not.allocated(dbuf) ) then
    call accum_mem('force_common',8*nbmax*2)
    allocate(dbuf(nbmax),dbufr(nbmax))
  else if( size(dbuf).ne.nbmax ) then
    call accum_mem('force_common',-8*size(dbuf)*2 +8*nbmax*2)
    deallocate(dbuf,dbufr)
    allocate(dbuf(nbmax),dbufr(nbmax))
  endif

  nbnew= 0

!-----loop over z, y, & x directions
  do kd=1,3
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
  enddo

  if(nbnew.ne.nb) then
    write(6,'(a,2i8)') "nbnew,(natm+nb)=",nbnew,natm+nb
    stop "error: nbnew.ne.(natm+nb)!!"
  endif

end subroutine copy_rho_ba
!=======================================================================
subroutine copy_strs_ba(namax,natm,nb,nbmax,lsb &
     ,lsrc,myparity,nn,sv,mpi_md_world,strs)
!-----------------------------------------------------------------------
!  Exchanges boundary-atom data among neighbor nodes
!-----------------------------------------------------------------------
  use memory, only: accum_mem
  implicit none
  include "mpif.h"
  integer:: status(MPI_STATUS_SIZE)
!-----in
  integer,intent(in):: namax,natm,nb,nbmax,mpi_md_world
  integer,intent(in):: lsb(0:nbmax,6),lsrc(6),myparity(3),nn(6)
  real(8),intent(in):: sv(3,6)
!-----out
  real(8),intent(inout):: strs(9,natm+nb)

!-----locals
  integer:: i,j,k,l,m,n,kd,kdd,ku,inode,nsd,nrc,nbnew,ierr,mem

  logical,save:: l1st=.true.
  real(8),save,allocatable:: dbuf(:,:),dbufr(:,:)
  integer,save:: narrsize

  if( l1st ) then
    narrsize = 9*nbmax
    call accum_mem('force_common',8*9*nbmax*2)
    allocate(dbuf(9,nbmax),dbufr(9,nbmax))
    l1st=.false.
  endif

  if( .not. allocated(dbuf) ) then
    call accum_mem('force_common',8*9*nbmax*2)
    allocate(dbuf(9,nbmax),dbufr(9,nbmax))
  else if( size(dbuf).ne.narrsize ) then
    deallocate(dbuf,dbufr)
    call accum_mem('force_common',8*9*nbmax*2 -8*size(dbuf)*2)
    narrsize = 9*nbmax
    allocate(dbuf(9,nbmax),dbufr(9,nbmax))
  endif

  nbnew= 0

!-----loop over z, y, & x directions
  do kd=1,3
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
  enddo

  if(nbnew.ne.nb) then
    write(6,'(a,2i8)') "nbnew,(natm+nb)=",nbnew,natm+nb
    stop "error: nbnew.ne.(natm+nb)!!"
  endif

end subroutine copy_strs_ba
!=======================================================================
subroutine copy_dba_fwd(namax,natm,nb,nbmax,lsb,nex &
     ,lsrc,myparity,nn,sv,mpi_md_world,x,ndim)
!-----------------------------------------------------------------------
!     Exchanges boundary-atom data among neighbor nodes
!-----------------------------------------------------------------------
  use memory, only: accum_mem
  implicit none
  include "mpif.h"
  integer:: status(MPI_STATUS_SIZE)
!-----in
  integer,intent(in):: namax,natm,nb,nbmax,mpi_md_world,ndim
  integer,intent(in):: lsb(0:nbmax,6),lsrc(6),myparity(3),nn(6) &
       ,nex(3)
  real(8),intent(in):: sv(3,6)
!-----out
  real(8),intent(inout):: x(ndim,namax)

!-----locals
  integer:: i,j,k,l,m,n,kd,kdd,ku,inode,nsd,nsd3,nrc,nrc3,nbnew,ierr
  logical,save:: l1st=.true.
  real(8),allocatable,save:: dbuf(:,:),dbufr(:,:)
  integer,save:: maxdim
  integer,save:: maxbmax

  if( l1st ) then
    maxdim = ndim
    maxbmax = nbmax
    call accum_mem('force_common',8*maxbmax*maxdim*2)
    allocate(dbuf(maxdim,maxbmax),dbufr(ndim,maxbmax))
    l1st=.false.
  endif

  if( ndim.gt.maxdim .or. nbmax.gt.maxbmax ) then
    call accum_mem('force_common', &
         -8*size(dbuf)-8*size(dbufr) +8*maxbmax*maxdim)
    maxdim = max(ndim,maxdim)
    maxbmax = max(nbmax,maxbmax)
    if( allocated(dbuf) ) then
      call accum_mem('force_common', -8*size(dbuf)-8*size(dbufr) )
      deallocate(dbuf,dbufr)
    endif
    allocate(dbuf(maxdim,maxbmax),dbufr(maxdim,maxbmax))
    call accum_mem('force_common', 8*size(dbuf) +8*size(dbufr) )
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
        call mespasd(inode,myparity(kd),dbuf,dbufr,maxdim*nsd &
             ,maxdim*nrc,21,mpi_md_world)
        do i=1,nrc
          x(1:ndim,natm+nbnew+i)= dbufr(1:ndim,i)
        enddo
!---------mpi barrier
        call mpi_barrier(mpi_md_world,ierr)
!---------accumulate num. of boundary particles
!          write(6,'(a,2i8)') "nbnew,nrc=",nbnew,nrc
        nbnew=nbnew +nrc
      enddo
    endif
  enddo

  if(nbnew.ne.nb) then
    write(6,'(a,2i8)') "nbnew,(natm+nb)=",nbnew,natm+nb
    stop "error: nbnew.ne.(natm+nb)!!"
  endif

end subroutine copy_dba_fwd
!=======================================================================
subroutine copy_dba_bk(namax,natm,nbmax,nb,lsb,nex &
     ,lsrc,myparity,nn,mpi_md_world,x,ndim)
!-----------------------------------------------------------------------
!     Send-back & receive reaction on cached-atoms
!-----------------------------------------------------------------------
  use memory, only: accum_mem
  implicit none
  include "mpif.h"
  integer,intent(in):: namax,natm,nbmax,nb,mpi_md_world,ndim
  integer,intent(in):: lsb(0:nbmax,6),lsrc(6),myparity(3),nn(6) &
       ,nex(3)
  real(8),intent(inout):: x(ndim,namax)

  integer:: status(MPI_STATUS_SIZE)
  integer:: i,j,k,l,m,n,kd,kdd,ku,kuc,ibkwd,nsd,nsd3,nrc,nrc3,nsdbk &
       ,ierr,natmx
  real(8),save,allocatable:: dbuf(:,:),dbufr(:,:)
  logical,save:: l1st=.true.
  integer,save:: maxdim
  integer,save:: maxbmax

  if( l1st ) then
    maxdim = ndim
    maxbmax = nbmax
    call accum_mem('force_common',8*maxdim*maxbmax*2)
    allocate(dbuf(maxdim,maxbmax),dbufr(ndim,maxbmax))
    l1st=.false.
  endif

  if( ndim.gt.maxdim .or. nbmax.gt.maxbmax ) then
    maxdim = max(ndim,maxdim)
    maxbmax = max(nbmax,maxbmax)
    if( allocated(dbuf) ) then
      call accum_mem('force_common', -8*size(dbuf) -8*size(dbufr))
      deallocate(dbuf,dbufr)
    endif
    call accum_mem('force_common', 8*maxdim*maxbmax*2)
    allocate(dbuf(maxdim,maxbmax),dbufr(maxdim,maxbmax))
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
        nsd3= maxdim*nsd
!---------num. of to-be-recieved particles
        nrc= lsb(0,ku)
!          nrc3= ndim*nrc
        nrc3= maxdim*nrc

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

    endif

  enddo

!-----check
  if(nsdbk.ne.nb) then
    write(6,'(a,2i8)') "nsdbk,nb=",nsdbk,nb
    stop "ERROR: nsdbk.ne.nb!!"
  endif

!      deallocate(dbuf,dbufr)
end subroutine copy_dba_bk
!=======================================================================
subroutine reduce_dba_bk(natm,namax,tag,x,ndim)
!-----------------------------------------------------------------------
!  Send-back or reduce reaction on cached-atoms.
!  This routine works only on small MD not on parallel version.
!-----------------------------------------------------------------------
  use util,only: itotOf
  implicit none
  integer,intent(in):: namax,natm,ndim
  real(8),intent(in):: tag(namax)
  real(8),intent(inout):: x(ndim,namax)
!!$  integer,external:: itotOf
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
  use util,only: itotOf
  implicit none
  integer,intent(in):: namax,natm,ndim
  real(8),intent(in):: tag(namax)
  real(8),intent(inout):: x(ndim,namax)
!!$  integer,external:: itotOf
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
function fcut1(r,rin,rout)
!
!     Cutof function type-1
!     f(r) = 1                                  for r <= rin
!          = 1/2 *[1+cos(pi*(r-rs)/(rc-rs))]    for rin < r <= rout
!          = 0                                  for rout < r
!
  real(8),intent(in):: r,rin,rout
  real(8),parameter:: pi = 3.14159265358979d0
!!$  real(8),parameter:: rsr= 0.9d0  ! rs ratio to rc
  real(8):: fcut1

!  rs = rc*rsr
  if( r.le.rin ) then
    fcut1 = 1d0
  else if( rin.lt.r .and. r.le.rout ) then
    fcut1 = 0.5d0 *(1d0 +cos(pi*(r-rin)/(rout-rin)))
  else
    fcut1 = 0d0
  endif
  return
end function fcut1
!=======================================================================
function dfcut1(r,rin,rout)
!
!     Derivative of the cutoff function type-1
!
  real(8),intent(in):: r,rin,rout
  real(8),parameter:: pi = 3.14159265358979d0
!!$  real(8),parameter:: rsr= 0.9d0  ! rs ratio to rc
  real(8):: dfcut1

!!$  rs = rc *rsr
  if( r.le.rin ) then
    dfcut1 = 0d0
  else if( rin.lt.r .and. r.le.rout ) then
    dfcut1 = -0.5d0 *pi/(rout-rin) *sin(pi*(r-rin)/(rout-rin))
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
subroutine chgopt_damping(chg,l1st)
!
!  Charge optimization/equilibration by damped dynamics.
!  Since not only Coulomb interaction but also other force-fields can
!  depend on atomic charges, this routine is separated out from the
!  force_Coulomb.F90.
!  Parameters related to chgopt_damping have "codmp" pre/postfix in their variable names.
!  Choices of damping method, given by codmp_method, are:
!    - damping --- simple, stable, but not very efficient.
!    - FIRE --- it is said that robust and fast, but sometimes unstable, do not know why.
!
  use pmdvars,only: namax,natm,tag,h,ra,nnmax,lspr,rc, &
       lsb,lsex,nbmax,nb,nn,myparity,lsrc,nex,sorg,myid_md, &
       mpi_md_world,iprint,boundary,ntot
  use force
  use Coulomb, only: qtop,qbot, get_qforce, &
       cterms,qforce_screened_cut, conv_eps_qeq, chgopt_method, &
       nstp_qeq, dt_codmp, qmass, minstp_qeq, &
       minstp_conv_qeq, finc_codmp, fdec_codmp, alpha0_codmp, falpha_codmp, &
       fdamp_codmp, dfdamp_codmp, impose_qtot
!!$  use Morse, only: qforce_vcMorse
  use memory, only: accum_mem
  implicit none
  include "mpif.h"
  include "./const.h"
  logical,intent(in):: l1st
!.....Output
  real(8),intent(inout):: chg(namax)

  integer:: istp,i,istp_pos,istp_conv,is,ierr,iflag
  real(8):: epot,epotp,epotpp,eMorse,fqnorm,vqnorm,dt,alpha,p &
       ,ecoul,dtmax,fdamp, dqmax
  real(8),save,allocatable:: vq(:),fq(:),dq(:)
!.....small value for checking range
  real(8),parameter:: qeps   = 1d-10

  integer,external:: count_nonfixed

  if( l1st ) then
    if( allocated(vq) ) then
      call accum_mem('force_common',-8*(size(vq)+size(fq)+size(dq)))
      deallocate(vq,fq,dq)
    endif
    allocate(vq(namax),fq(namax),dq(namax))
    call accum_mem('force_common',8*(size(vq) +size(fq)+size(dq)))
  endif

!!$  call mpi_allreduce(natm,ntot,1,mpi_integer,mpi_sum,mpi_md_world,ierr)

!.....Gather forces on charges
  istp = 0
  call impose_qtot(chg)
  call get_fq_wrapper(chg,fq,epot,l1st)

  vq(:) = 0d0
  fdamp = fdamp_codmp
  istp_pos = 0
  istp_conv = 0
  alpha = alpha0_codmp
  dt = dt_codmp
  dtmax = dt_codmp *10
  epotp = 0d0
  do istp=1,nstp_qeq
    epotpp = epotp
    epotp = epot
!.....velocity update
    vq(1:natm) = vq(1:natm) +dt/qmass*fq(1:natm)

    if( trim(chgopt_method).eq.'FIRE' ) then
      p = dot_product(fq(1:natm),vq(1:natm))
      vqnorm = sqrt(dot_product(vq(1:natm),vq(1:natm)))
      fqnorm = sqrt(dot_product(fq(1:natm),fq(1:natm)))
      vq(1:natm) = (1d0-alpha)*vq(1:natm) -alpha*vqnorm*fq(1:natm)/fqnorm
      if( p.gt.0d0 ) then
        istp_pos = istp_pos +1
        if( istp_pos.gt.minstp_qeq ) then
          dt = min(dt*finc_codmp,dtmax)
          alpha = alpha *falpha_codmp
          istp_pos = 0
        endif
      else if( p.le.0d0 ) then
        istp_pos = 0
        dt = dt*fdec_codmp
        alpha = alpha0_codmp
        vq(1:natm) = 0d0
      endif
    else  ! simple velocity damping
      do i=1,natm
        if( vq(i)*fq(i).lt.0d0 ) then
          vq(i) = dt/qmass *fq(i)
        else
          vq(i) = vq(i)*fdamp
        endif
      enddo
    endif

!.....update charges
    dq(1:natm) = vq(1:natm)*dt
    call suppress_dq(dq,0.2d0,dqmax)
    chg(1:natm)= chg(1:natm) +dq(1:natm)
!.....Check charges
    if( iprint.ge.ipl_debug ) then
      do i=1,natm
        if( chg(i)*0d0 .ne. 0d0 ) then
          print *,'ERROR: chg is NaN !?'
          print '(a,i8,es12.3)',' i,chg(i)= ',i,chg(i)
          stop
        endif
      enddo
    endif

    call impose_qtot(chg)
!!$    print '(a,200es11.3)',' chg=',chg(1:natm)
    call get_fq_wrapper(chg,fq,epot,.false.)
    if( myid_md.eq.0 .and. iprint.ge.ipl_debug ) then
      write(6,'(a,i5,2es16.8,4es11.3)') ' istp,epot/ntot,de/ntot,dqmax= ',istp &
           ,epot/ntot,abs(epot-epotp)/ntot,dqmax,chg(1),chg(5)
    endif

!.....check convergence
    if( istp.gt.minstp_qeq ) then
      if( abs(epot-epotp)/ntot.lt.conv_eps_qeq ) then
        istp_conv = istp_conv +1
        if( istp_conv.ge.minstp_conv_qeq ) then
          if( myid_md.eq.0 .and. iprint.ge.ipl_info ) then
            write(6,'(a,i0,a)') ' chgopt_damping converged at ', &
                 istp,' steps.'
          endif
          exit
        endif
      else
        istp_conv = 0
      endif
    endif
  enddo

  if( istp.ge.nstp_qeq .and. myid_md.eq.0 .and. iprint.ge.ipl_info ) then
    write(6,'(a,i0,a)') ' chgopt_damping not converged within ', istp,' steps.'
  endif

  return
end subroutine chgopt_damping
!=======================================================================
subroutine chgopt_xlag(chg,auxq,iflag)
!
!  Charge optimization via extended Lagrangian proposed by Nomura[1].
!  Charges are optimized only once (line minimization) along the forces
!  on charges. But the initial guess of charges are given by auxiliary
!  variables corresponding to charges.
!
!  [1]: Nomura, et al., Comput. Phys. Commun. 192, 91–96 (2015)
!
  use pmdvars,only: namax,natm,iprint,myid_md,mpi_md_world,ntot
  use Coulomb,only: qmass,conv_eps_qeq
  use memory, only: accum_mem
  implicit none
  include "mpif.h"
  include "./const.h"
  real(8),intent(in):: auxq(namax)
  real(8),intent(inout):: chg(namax)
  integer,intent(out):: iflag

  logical,save:: l1st = .true.
  real(8),save,allocatable:: fq(:)
  integer:: i,imax,ierr
  real(8):: alpha,epot,dqmax,dq

  if( l1st ) then
    if( allocated(fq) ) then
      call accum_mem('force_common', -8*size(fq))
      deallocate(fq)
    endif
    allocate(fq(namax))
    call accum_mem('force_common',+8*size(fq))
    l1st = .false.
  endif

  dqmax = -1d0
  imax = -1
  do i=1,natm
    dq = abs(chg(i) -auxq(i))
!!$    print '(a,i5,3f12.2)','i,chg,auxq,dq=',i,chg(i),auxq(i),dq
    if( dqmax.lt.dq ) then
      dqmax = dq
      imax = i
    endif
  enddo

!.....Set initial guess of charges
  chg(:) = auxq(:)
!.....Get forces on charges
  call get_fq_wrapper(chg,fq,epot,.false.)
!!$  print *,'epot/ntot = ',epot/ntot

  iflag = 0
  call linmin_chg(chg,fq,conv_eps_qeq,alpha,epot,iflag)
  if( iflag.ne.0 ) then
    if( myid_md.eq.0 ) then
      print *,'something wrong with linmin_chg? iflag = ',iflag
    endif
  else
    chg(1:natm) = chg(1:natm) +fq(1:natm)*alpha
  endif
  return
end subroutine chgopt_xlag
!=======================================================================
subroutine get_fq_wrapper(chg,fq,epot,l1st)
!
!  Wrapper routine for getting forces on charges.
!
  use pmdvars,only: namax,myid_md
  use Coulomb,only: get_qforce,bacopy_chg_fixed
  real(8),intent(inout):: chg(namax)
  real(8),intent(out):: fq(namax),epot
  logical,intent(in):: l1st

  real(8):: epotb

  call bacopy_chg_fixed(chg)
  call get_qforce(chg,fq,epot,l1st)
  call get_bounding_fq(chg,fq,epotb)
  epot = epot +epotb
  call set_fqtot_zero(fq)
  return
end subroutine get_fq_wrapper
!=======================================================================
subroutine set_fqtot_zero(fq)
!
!  Set total force on q's as zero, which means the conservation of total charge.
!
  use pmdvars,only: namax,natm,myid_md,mpi_md_world,ntot
  implicit none
  include "mpif.h"
  real(8),intent(inout):: fq(namax)

  integer:: i,ierr
  real(8):: avmul,avmu

  avmul = 0d0
  do i=1,natm
    avmul = avmul +fq(i)
  enddo
  avmu = 0d0
  call mpi_allreduce(avmul,avmu,1,mpi_real8,mpi_sum,mpi_md_world,ierr)
  avmu = avmu/ntot
  do i=1,natm
    fq(i) = fq(i) -avmu
  enddo
  
  return
end subroutine set_fqtot_zero
!=======================================================================
subroutine get_bounding_fq(chg,fq,epot)
!
!  Compute forces on qs that bound qs inside [qbot,qtop].
!
  use pmdvars,only: namax,natm,iprint,tag,mpi_md_world
  use Coulomb,only: bound_k2,bound_k4,qtop,qbot
  implicit none
  include 'mpif.h'
  include './const.h'
  real(8),intent(in):: chg(namax)
  real(8),intent(inout):: fq(namax)
  real(8),intent(out):: epot

  integer:: i,is,ierr
  real(8):: epotl

  epotl = 0d0
  
  do i=1,natm
    is = int(tag(i))
    if( chg(i).ge.qtop(is) ) then
      epotl = epotl +bound_k2*(chg(i)-qtop(is))**2 &
           +bound_k4*(chg(i)-qtop(is))**4
      fq(i) = fq(i) -2d0*bound_k2*(chg(i)-qtop(is)) &
           -4d0*bound_k4*(chg(i)-qtop(is))**3
    else if( chg(i).le.qbot(is) ) then
      epotl = epotl +bound_k2*(chg(i)-qbot(is))**2 &
           +bound_k4*(chg(i)-qbot(is))**4
      fq(i) = fq(i) -2d0*bound_k2*(chg(i)-qbot(is)) &
           -4d0*bound_k4*(chg(i)-qbot(is))**3
    endif
  enddo
  epot = 0d0
  call mpi_allreduce(epotl,epot,1,mpi_real8,mpi_sum,mpi_md_world,ierr)
  return
end subroutine get_bounding_fq
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
subroutine suppress_dq(dq,dqlim,dqmax)
!
!  Suppress too large change of charge within dq < dqlim.
!
  use pmdvars,only: namax,natm,tag,iprint,mpi_md_world,myid_md
  implicit none
  include "mpif.h"
  real(8),intent(in):: dqlim
  real(8),intent(inout):: dq(namax)
  real(8),intent(out):: dqmax

  integer:: i,ierr
  real(8):: dqmaxl,fac
  
  dqmaxl = 0d0
  do i=1,natm
    dqmaxl = max(dqmaxl,abs(dq(i)))
  enddo
  dqmax = 0d0
  call mpi_allreduce(dqmaxl,dqmax,1,mpi_real8,mpi_max,mpi_md_world,ierr)
!!$  if( myid_md.eq.0 ) print *,'dqmax,dqlim=',dqmax,dqlim
  if( dqmax.gt.dqlim ) then
    fac = dqlim/dqmax
    do i=1,natm
      dq(i) = dq(i)*fac
    end do
    dqmax = dqlim
  end if
  return
end subroutine suppress_dq
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

    fqfactor = 1d0
    if( fqmax.gt.1d-8 ) fqfactor = fqlimit /fqmax
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
subroutine get_gradw(namax,natm,tag,ra,nnmax,aa,strs,chg &
     ,h,hi,nb,nbmax,lsb,nex,lsrc,myparity,nnn,sv,rc,lspr &
     ,mpi_md_world,myid_md,epi,epot,nismax,lstrs &
     ,iprint,l1st,lvc &
     ,ndimp,gwe,gwf,gws)
!
!  Compute derivative of potential energy (and forces) 
!  w.r.t. potential parameters.
!
  use force
!!$  use Morse,only: gradw_vcMorse, gradw_Morse
  implicit none
  integer,intent(in):: namax,natm,nnmax,nismax,iprint
  integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3) &
       ,nnn(6),mpi_md_world,myid_md,nex(3)
  integer,intent(in):: lspr(0:nnmax,namax) !,numff
  real(8),intent(inout):: ra(3,namax),h(3,3,0:1),hi(3,3),sv(3,6) &
       ,tag(namax)
  real(8),intent(inout):: rc
  real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax) &
       ,chg(namax)
!!$    character(len=20),intent(in):: cffs(numff)
  logical,intent(in):: l1st
  logical,intent(inout):: lvc
  logical,intent(in):: lstrs
  integer,intent(in):: ndimp
  real(8),intent(out):: gwe(ndimp),gwf(ndimp,3,natm),gws(ndimp,6)

!!$  if( use_force('vcMorse') ) call gradw_vcMorse(namax,natm,tag,ra,nnmax,chg &
!!$       ,h,rc,lspr,epot,iprint,ndimp,gwe,gwf,gws)

end subroutine get_gradw
!=======================================================================
subroutine linmin_chg(chg0,dchg,ftol,alpha,falpha,iflag)
!
!  Search for coeff, alpha, that are multiplied to fq to get the charge
!  distribution of minimum energy.
!
  use pmdvars, only: namax,natm,nnmax,tag,ra,h,lspr,rc, &
       lsb,lsex,nbmax,nb,nn,myparity,lsrc,nex,boundary,sorg, &
       myid_md,mpi_md_world,iprint,ntot
  use Coulomb, only: get_qforce,qbot,qtop
  implicit none
  include 'mpif.h'
  include './const.h'
!.....input
  real(8),intent(in):: chg0(namax),dchg(namax)
  real(8),intent(in):: ftol
!.....output
  real(8),intent(out):: alpha,falpha
  integer,intent(out):: iflag

  real(8),parameter:: stp0 = 0.1d0
  real(8),parameter:: gr = 0.61803398875d0
  real(8),parameter:: gr2 = 1d0 -gr
  integer,parameter:: nstp = 100
  real(8),save,allocatable:: chgt(:),fqt(:)

  integer:: istp
  real(8):: a,b1,b2,xl,c,fa,fb1,fb2,fc

  if( .not.allocated(chgt) ) then
    allocate( chgt(namax),fqt(namax) )
  else if( size(chgt).ne.namax ) then
    deallocate(chgt,fqt)
    allocate( chgt(namax),fqt(namax) )
  endif

  a = 0d0
  b1 = stp0
  call get_range(chg0,dchg,a,b1,c,fa,fb1,fc,iflag)
  if( iflag/1000.ne.0 ) return
  xl = (c-a)
  b1 = a +gr2*xl
  b2 = a +gr *xl

  chgt(:) = chg0(:)+b1*dchg(:)
  call get_fq_wrapper(chgt,fqt,fb1,.false.)
  chgt(:) = chg0(:)+b2*dchg(:)
  call get_fq_wrapper(chgt,fqt,fb2,.false.)

  istp = 0
  if( iprint.ge.ipl_debug .and. myid_md.eq.0 ) then
    print '(a,i5,4f12.8,4es16.8)','   linmin_chg: istp,a,b1,b2,c,fa,fb1,fb2,fc=', &
         istp,a,b1,b2,c,fa,fb1,fb2,fc
  endif

10 continue
  istp = istp +1
  if( istp.gt.nstp ) then
    if( myid_md.eq.0 .and. iprint.ge.ipl_info ) then
      print *,'WARNING @linmin: istp.gt.nstp !!'
    endif
    iflag = iflag +100
    return
  endif

  if( fb1.gt.fb2 ) then
    a = b1
    fa = fb1
    b1 = b2
    fb1 = fb2
    xl = (c-a)
    b2 = a +gr*xl
    chgt(:) = chg0(:)+b2*dchg(:)
    call get_fq_wrapper(chgt,fqt,fb2,.false.)
    if( iprint.ge.ipl_debug .and. myid_md.eq.0 ) then
      print '(a,i5,4f12.8,4es16.8)','   linmin_chg: istp,a,b1,b2,c,fa,fb1,fb2,fc=', &
           istp,a,b1,b2,c,fa,fb1,fb2,fc
    endif
  else
    c = b2
    fc = fb2
    b2 = b1
    fb2 = fb1
    xl = (c-a)
    b1 = a +gr2*xl
    chgt(:) = chg0(:)+b1*dchg(:)
    call get_fq_wrapper(chgt,fqt,fb1,.false.)
    if( iprint.ge.ipl_debug .and. myid_md.eq.0 ) then
      print '(a,i5,4f12.8,4es16.8)','   linmin_chg: istp,a,b1,b2,c,fa,fb1,fb2,fc=', &
           istp,a,b1,b2,c,fa,fb1,fb2,fc
    endif
  endif

  if( abs(fb1-fb2).lt.ftol ) then
    iflag = 0
    if( fb1.gt.fb2 ) then
      alpha = b2
      falpha = fb2
    else
      alpha = b1
      falpha = fb1
    endif
    return
  endif
  goto 10

end subroutine linmin_chg
!=======================================================================
subroutine get_range(chg0,dchg,a,b,c,fa,fb,fc,iflag)
!
!  Get a range of factor of line minimization in QEq.
!
  use pmdvars,only: namax,natm,nnmax,tag,ra,h,lspr,rc,sorg, &
       myid_md,mpi_md_world,iprint,ntot
  use Coulomb,only: get_qforce
  include './const.h'
  real(8),intent(in):: chg0(namax),dchg(namax)
  real(8),intent(inout):: a,b,c,fa,fb,fc
  integer,intent(inout):: iflag
  
  real(8),parameter:: gr = 1.61803398875d0  ! golden ratio
  real(8),parameter:: gri = 1d0/gr
  real(8),parameter:: eps = 1d-12
  real(8),parameter:: nstp = 100
  integer:: istp
  real(8),save,allocatable:: fqt(:),chgt(:)

  if( .not.allocated(fqt) ) then
    allocate(fqt(namax),chgt(namax))
  else if( size(fqt).ne.namax ) then
    deallocate(fqt,chgt)
    allocate(fqt(namax),chgt(namax))
  endif

  chgt(:) = chg0(:)+a*dchg(:)
  call get_fq_wrapper(chgt,fqt,fa,.false.)
  chgt(:) = chg0(:)+b*dchg(:)
  call get_fq_wrapper(chgt,fqt,fb,.false.)

  istp = 0

  if( iprint.ge.ipl_debug .and. myid_md.eq.0 ) then
    print '(a,i5,2f12.8,2es16.8)','   get_range: istp,a,b,fa,fb=', &
         istp,a,b,fa,fb
  endif
  
10 continue
  istp = istp +1
  if( istp.gt.nstp ) then
    if( myid.eq.0 .and. iprint.ge.ipl_info ) then
      print *,'WARNING @get_bracket: istp exceeds nstp.'
    endif
    iflag = iflag +1000
    return
  endif
  if( abs(b-a).lt.eps ) then
    if( myid.eq.0 .and. iprint.ge.ipl_info ) then
      print *,'WARNING @get_bracket: a and b are too close.'
      print *,'  Search direction may not be a descent direction.'
    endif
    iflag = iflag +2000
    return
  endif
  if( fa.lt.fb ) then
    c = a +gri*(b-a)
    chgt(:) = chg0(:)+c*dchg(:)
    call get_fq_wrapper(chgt,fqt,fc,.false.)
    call swap_vals(c,b)
    call swap_vals(fc,fb)
    if( iprint.ge.ipl_debug .and. myid_md.eq.0 ) then
      print '(a,i5,2f12.8,2es16.8)','   get_range: istp,a,b,fa,fb=', &
           istp,a,b,fa,fb
    endif
    goto 10
  else
    c = a +gr*(b-a)
    chgt(:) = chg0(:)+c*dchg(:)
    call get_fq_wrapper(chgt,fqt,fc,.false.)
    if( fb.gt.fc ) then
      b = a +gr*(c-a)
      chgt(:) = chg0(:)+b*dchg(:)
      call get_fq_wrapper(chgt,fqt,fb,.false.)
      call swap_vals(b,c)
      call swap_vals(fb,fc)
      if( iprint.ge.ipl_debug .and. myid_md.eq.0 ) then
        print '(a,i5,2f12.8,2es16.8)','   get_range: istp,a,b,fa,fb=', &
             istp,a,b,fa,fb
      endif
      goto 10
    endif
  endif
  return
end subroutine get_range
!=======================================================================
subroutine swap_vals(a,b)
  real(8),intent(inout):: a,b
  real(8):: tmp
  tmp = b
  b = a
  a = tmp
  return
end subroutine swap_vals
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd lib"
!     End:
