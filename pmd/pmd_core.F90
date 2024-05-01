!-----------------------------------------------------------------------
!                     Last-modified: <2024-04-10 16:48:46 KOBAYASHI Ryo>
!-----------------------------------------------------------------------
! Core subroutines/functions needed for pmd.
!-----------------------------------------------------------------------
subroutine pmd_core(hunit,hmat,ntot0,tagtot,rtot,vtot,atot,stot &
     ,ekitot,epitot,auxtot,epot,ekin,stnsr)
!.....All the arguments are in pmdvars module
  use pmdio,only: write_pmdtot_ascii, write_pmdtot_bin, write_dump
  use util,only: iauxof, calc_nfmv, cell_info
  use pmdvars
  use zload
  use force
  use vector,only: dot
  use random,only: box_muller
  use ttm,only: init_ttm,langevin_ttm,output_ttm, &
       calc_Ta,update_ttm,assign_atom2cell,output_energy_balance, &
       remove_ablated_atoms,set_inner_dt, te2tei, non_reflecting_bc, &
       set_3d1d_bc_pos, lcouple_3d1d, couple_3d1d, compute_nac
  use pmdmpi,only: nid2xyz,xyz2nid
  use metadynamics,only: init_metaD,update_metaD,force_metaD &
       ,write_metaD_potential
  use constraints,only: init_const, update_const, update_const_vel &
       ,update_const_pos
  use rdcfrc,only: init_rdcfrc, reduce_forces, finalize_rdcfrc
  use structure,only: cna,acna
  use deform,only: init_deform, apply_deform, cdeform
  use util,only: itotOf, ifmvOf, ithOf
  use extforce,only: lextfrc,rm_trans_extfrc,add_extfrc
  use clrchg,only: lclrchg,clrchg_force,rm_trans_clrchg
  use localflux,only: lflux,accum_lflux
  use pdens,only: lpdens,accum_pdens
  use time, only: sec2hms, accum_time
  use pairlist, only: mk_lspr_para
  use Coulomb,only: chgopt_method, update_auxq, update_vauxq, get_aauxq
  use isostat,only: setup_langevin, vel_update_langevin, vel_update_berendsen, &
       cell_update_berendsen,cell_force_berendsen, setup_cell_langevin, &
       cell_update_langevin, cvel_update_langevin, setup_cell_berendsen, &
       setup_cell_langevin
  use descriptor,only: write_desc,lout_desc
  use dspring, only: ldspring, init_dspring, final_dspring, force_dspring

  implicit none
  include "mpif.h"
  include "./params_unit.h"
  include "./const.h"
  integer,intent(inout):: ntot0
  real(8),intent(in):: hunit
  real(8),intent(inout):: tagtot(ntot0),rtot(3,ntot0),vtot(3,ntot0), &
       hmat(3,3,0:1)
  real(8),intent(out):: atot(3,ntot0),stot(3,3,ntot0), &
       ekitot(3,3,ntot0),epitot(ntot0),auxtot(naux,ntot0)
  real(8),intent(out):: epot,ekin,stnsr(3,3)

  integer:: i,j,k,l,m,n,ia,ib,is,ifmv,itemp,nave,nspl,i_conv,ierr
  integer:: ihour,imin,isec
  real(8):: tmp,hscl(3),aai(3),ami,tave,vi(3),vl(3),epotp, &
       htmp(3,3),prss,dtmax,vmaxt,rbufres,tnow,sth(3,3)
  logical:: l1st
  logical:: lconverged = .false.
!.....FIRE variables
  real(8):: alp_fire,fnorm,vnorm,fdotv
  integer:: num_fire
!.....For tensile test of Al nanorod
  real(8):: strnow,ftop,fbot
!-----output file names
  character:: cnum*128, ctmp*128
  logical:: ltot_updated = .true.
!.....Formats for output
  character(len=20):: cfistp  = 'i10' !or larger
  character(len=20):: cfstime = 'es18.10'
  character(len=20):: cfetime = 'f12.2' ! for elapsed time
  character(len=20):: cftave  = 'f12.2' ! or 'es12.4' for high-T

  tcpu0= mpi_wtime()
  h(:,:,:) = hmat(:,:,:)  ! use pmdvars variable h instead of hmat
  ntot = ntot0
  call mpi_bcast(ntot,1,mpi_integer,0,mpi_md_world,ierr)

  if( nstp.le.0 ) then
    cfistp = 'i2'
  else
    write(cfistp(2:3),'(i2.2)') int(log10(dble(nstp)))+2
  endif

!.....Variable time-step
  if( dt.lt.0d0 ) then
    dtmax = abs(dt)
    lvardt = .true.
    if( myid_md.eq.0 .and. iprint.ge.ipl_basic ) then
      print *,''
      print '(a,f8.3,a)',' Use variable time-step: dtmax = ' &
           ,dtmax,' fs'
    endif
    dt = dtmax
  endif

!.....Set 0d0 to z of fmv(9) in case of z-loading
  if( trim(czload_type).eq.'atoms' .or. &
       trim(czload_type).eq.'box' ) then
    fmv(3,9) = 0d0
    if( myid_md.eq.0 .and. iprint.ge.ipl_basic ) then
      print *,''
      write(6,'(a)') 'Set fmv(3,9)=0d0, because czload_type=' &
           //trim(czload_type)
    endif
  else if( trim(czload_type).eq.'shear' ) then
    fmv(1,9) = 0d0
    fmv(3,9) = 0d0
    if( myid_md.eq.0 .and. iprint.ge.ipl_basic ) then
      print *,''
      write(6,'(a)') 'Set fmv(1,9) and fmv(3,9) = 0d0,' &
           //' because czload_type='//trim(czload_type)
    endif
  endif

!-----parallel configuration
  nxyz= nx*ny*nz
  anxi= 1d0/nx
  anyi= 1d0/ny
  anzi= 1d0/nz
!-----error trap
  if(nodes_md.ne.nxyz) then
    write(6,'(a)') " error: nodes_md .ne. nxyz!!"
    write(6,'(a,3i5)') ' myid_md,nodes_md,nxyz=' &
         ,myid_md,nodes_md,nxyz
    print *,'nx,ny,nz =',nx,ny,nz
    call mpi_finalize(ierr)
    stop
  endif
!-----vector node indices: range [0:nx-1]
  call nid2xyz(myid_md,myx,myy,myz)
!-----reduced node origin
  sorg(1)= anxi*myx
  sorg(2)= anyi*myy
  sorg(3)= anzi*myz

!     check number of parallel node and vc_gaussian
  if( nxyz.gt.1 .and. use_force('long_Coulomb') ) then
    write(6,'(a)') ' [Error] vc_gaussian potential is ' &
         //'not available in parallel.'
    call mpi_finalize(ierr)
    stop
  endif

!-----setup
  call setup(nspmax,am,fekin,fa2v)
!-----set HI and SGM
  call boxmat(h,hi,ht,g,gi,gt,vol,sgm)
!-----ntset
  call ntset(myx,myy,myz,nx,ny,nz,nn,sv,myparity,anxi,anyi,anzi)

!!$  print *,'Time at 0 = ',mpi_wtime() -tcpu0
!-----output every these steps, NOUTERG, NOUTPMD
  if( nerg.gt.0 ) then
    nouterg = max(nstp/nerg,1)
  else
    nouterg = nstp +1
  endif
  if( npmd.ne.0 ) then
    noutpmd = max(nstp/npmd,1)
  else
    noutpmd = nstp +1
  endif
!.....perform space decomposition after reading atomic configuration
  tmp = mpi_wtime()
  call space_decomp(ntot0,tagtot,rtot,vtot,auxtot)
  call accum_time('space_decomp',mpi_wtime()-tmp)
!.....Some conversions
  do i=1,natm
    ra(1:3,i)= ra(1:3,i) -sorg(1:3)
  enddo

#ifdef __FITPOT__
!.....check whether order of atoms and total-id of atoms match
  do ia=1,natm
    if( itotOf(tag(ia)).ne.ia ) then
      print *, '[Error] itotOf(tag(ia)).ne.ia !!!'
      print *, '  In case of FITPOT mode, the order of atom'// &
           ' must be the same of itot of the atom-tag.'
      stop
    endif
  enddo
#endif

!!$  print *,'Time at 1 = ',mpi_wtime() -tcpu0

!.....NEMD setting
  if( ltdst ) then
    if( mod(ntdst,nx).ne.0 ) &
         call error_mpi_stop('mod(ntdst,nx).ne.0')
    allocate(tdst(ntdst),nadst(ntdst))
  endif

!.....Deformation setup
  if( trim(cdeform).ne.'none' ) then
    call init_deform(h,myid_md,iprint)
  endif

!.....get_num_dof is called once in a MD run
!!$  call get_num_dof(natm,tag,fmv,ndof,myid_md,mpi_md_world,iprint)
  call get_num_dof()

!.....Call init_metaD after setting namax and nsp
  if( lmetaD ) then
    call init_metaD(nstp,nsp,namax,nodes_md &
         ,myid_md,mpi_md_world,iprint)
  endif
!.....Initialize constraints if needed
  if( lconst ) call init_const(myid_md,mpi_md_world,nodes_md &
       ,iprint,h)
!.....Initialize rdcfrc if needed
  if( lrdcfrc ) call init_rdcfrc(myid_md,mpi_md_world,nodes_md &
       ,iprint)

!.....Setup for FIRE
  if( ifdmp.eq.2 ) then
    alp_fire = alp0_fire
!        am(1:nspmax) = 20.d0
    num_fire = 0
    dtmax_fire = dtmfctr_fire * dt
    if( myid_md.eq.0 .and. iprint.ne.0 ) then
      write(6,'(/,a)') ' FIRE relaxation parameters:'
      write(6,'(a,f0.2,a)') '   Max dt = ',dtmax_fire,' fs'
      write(6,'(a,f0.2)') '   Factor of dt increment = ' &
           ,finc_fire
      write(6,'(a,f0.2)') '   Default alpha = ',alp0_fire
    endif
  endif

!.....Set ifmv of top and bottom atoms for z-loading if needed
  if( trim(czload_type) .eq. 'atoms' ) then
    call set_zload_atoms(natm,ra,tag,h,fmv,sorg,strfin,nstp &
         ,zskin_width,myid_md,mpi_md_world,iprint)
  else if( trim(czload_type).eq.'shear' ) then
    call set_shear(natm,ra,tag,h,fmv,sorg,strfin,nstp &
         ,zskin_width,zshear_angle,myid_md,mpi_md_world,iprint)
  endif

!.....Set initial temperature if needed
  if( tinit.gt.1d-5 ) then
    call setv(h,hi,natm,tag,va,nspmax,am,tinit,dt)
  elseif( abs(tinit).le.1d-5 ) then
    va(1:3,1:natm)= 0d0
  endif
  if( lclrchg ) then  ! special treatment for translational momentum
    call rm_trans_clrchg(natm,tag,va,am,mpi_md_world,myid_md,iprint)
  elseif( nrmtrans.ge.0 ) then
    call rm_trans_motion(natm,tag,va,nspmax,am,mpi_md_world,myid_md,iprint)
  endif

  if( ifdmp.eq.2 ) then
    va(1:3,1:natm) = 0d0
    h(1:3,1:3,1) = 0d0
  endif

!-----calc kinetic energy
  call get_ekin(namax,natm,va,tag,h,nspmax,fekin,ekin,eki,eks &
       ,vmax,mpi_md_world)
  vmaxold=vmax

  if( index(ctctl,'langevin').ne.0 ) then
    call setup_langevin(myid_md,iprint)
  else if( index(ctctl,'ttm').ne.0 ) then
    tmp = mpi_wtime()
    call init_ttm(namax,natm,ra,h,sorg,dt,lvardt, &
         boundary,myid_md,mpi_md_world,iprint)
!!$    call assign_atom2cell(namax,natm,ra,sorg,boundary)
    call calc_Ta(namax,natm,nspmax,h,tag,va,fmv,fekin &
         ,0,myid_md,mpi_md_world,iprint)
    call te2tei(namax,natm,aux(iaux_tei,:))
    if( lcouple_3d1d ) then
      call couple_3d1d(myid_md,mpi_md_world,iprint)
    endif
    call output_ttm(0,simtime,myid_md,iprint)
    call accum_time('ttm',mpi_wtime()-tmp)
  endif

  tcpu1= mpi_wtime()

!!$  call init_force(.true.)
!-----copy RA of boundary atoms
  call check_size_and_parallel(sgm,vol,rc,anxi,anyi,anzi &
       ,nx,ny,nz,myid_md)
  l1st = .true.
  tmp = mpi_wtime()
  call bacopy(.true.)
  call accum_time('ba_xxx',mpi_wtime()-tmp)
!-----Make pair list
  tmp = mpi_wtime()
!$acc update device(ra,h)
!!$  call mk_lspr_para(namax,natm,nbmax,nb,nnmax,tag,ra,va,rc+rbuf &
!!$       ,h,hi,anxi,anyi,anzi,lspr,iprint,l1st)
  call mk_lspr_para(l1st)
  call accum_time('lspr',mpi_wtime()-tmp)

  if( iprint.gt.0 .and. myid_md.eq.0 ) then
    print '(/a,i5)', ' Max num of neighbors = ',maxnn
  endif

  if( chgopt_method(1:4).eq.'xlag' ) then
    aux(iaux_vq,:) = 0d0
  endif

!.....Output descriptor and stop
  if( lout_desc ) then
    if( nodes_md.ne.1 ) stop 'ERROR: write_desc cannot be done in parallel.'
    call write_desc(namax,natm,nnmax,lspr,h,tag,ra,rc, &
         myid_md,mpi_md_world,iprint)
  endif

  if( ldspring ) call init_dspring(myid_md,mpi_md_world,iprint)

!!$  print *,'Time at 2 = ',mpi_wtime() -tcpu0
!.....Calc forces
  lstrs = lstrs0 .or. (index(cpctl,'beren').ne.0)
!.....Cell is new at the first call of get_force
  lcell_updated = .true.
  tmp = mpi_wtime()
  call get_force(.true.,epot0,stnsr)
  if( chgopt_method(1:4).eq.'xlag' ) then
    call get_aauxq(aux(iaux_chg,:),aux(iaux_q,:))
  endif
  call accum_time('get_force',mpi_wtime()-tmp)
  lcell_updated = .false.
  lstrs = .false.
  epot= epot0
  epotp = 0d0

!.....Descriptor spring
  if( ldspring ) call force_dspring(namax,natm,nnmax,lspr,rc,h,hi,tag,ra, &
       aa,aux(iaux_edsp,:), &
       nb,nbmax,lsb,nex,lsrc,myparity,nn,myid_md,mpi_md_world,iprint,.true.)
!.....Structure analysis
  if( trim(cstruct).eq.'CNA' ) then
    call cna(namax,natm,h,ra,nb,nnmax,lspr,rc_struct)
  else if( trim(cstruct).eq.'a-CNA' ) then
    call acna(namax,natm,h,ra,nb,nnmax,lspr,rc_struct)
  endif
!.....Color charge NEMD
  if( lclrchg ) call clrchg_force(namax,natm,tag,aa,aux,hi,specorder &
       ,myid_md,iprint)
!.....External force
  if( lextfrc ) call add_extfrc(natm,tag,aa,hi,specorder,myid_md,iprint)
!.....Constraints
  if( lconst ) then
    call update_const(namax,natm,tag,ra,h,0,nstp)
  endif
!.....Force_modify
  if( lrdcfrc ) then
    call reduce_forces(namax,natm,aa,tag,ra,h,nnmax,lspr)
  endif

#ifdef __DISL__
  call perf_disl_pos_by_pot(epith,natm,ra,h,epi,sorg &
       ,nodes_md,myid_md,mpi_md_world,0,21)
#endif

  call sa2stnsr(natm,strs,eki,stnsr,vol,mpi_md_world)
  if( index(cpctl,'beren').ne.0 ) then
    call setup_cell_berendsen(myid_md,iprint)
    call cell_force_berendsen(stnsr,ah,mpi_md_world)
  else if( index(cpctl,'lange').ne.0 ) then
    call setup_cell_langevin(myid_md,iprint)
    call cvel_update_langevin(stnsr,h,mpi_md_world,2)
  endif

  istp= 0
  simtime = 0d0
  rbufres = rbuf

  if(myid_md.eq.0 .and. iprint.ne.0 ) then
    write(6,*) ''
    write(6,'(1x,a)') "Initial values:"
    write(6,'(1x,a,f16.5,a,f10.3,a)') "  Kinetic energy  = ",ekin &
         ,' eV = ',ekin/ntot0,' eV/atom'
    write(6,'(1x,a,f16.5,a,f10.3,a)') "  Potential energy= ",epot0 &
         ,' eV = ',epot0/ntot0,' eV/atom'
    nave= 0
    tave= 0d0
    do itemp=1,ntemps
      if( ndof(itemp).eq.0 ) cycle
      temps(itemp)= eks(itemp) /max(ndof(itemp)-3,3) /fkb *2d0
      nave= nave +ndof(itemp)
      if( lmultemps ) then
        write(6,'(1x,a,i1,a,f16.5,a)') "  Temperature ",itemp &
             ,"   = ",temps(itemp),' K'
      endif
      tave= tave +temps(itemp)*max(ndof(itemp)-3,3)
    enddo
    tave= tave/(nave-3)
    write(6,'(1x,a,f16.5,a)') "  Temperature     = ",tave,' K'

!.....Human-friendly stress unit GPa
    sth(:,:) = stnsr(:,:) *up2gpa
    prss = (sth(1,1)+sth(2,2)+sth(3,3))/3
    if( prss.lt.0d0 ) then
      ctmp = '(tensile)'
    else
      ctmp = '(compressive)'
    endif
    write(6,'(1x,a,f16.5,a)') "  Pressure        = ", &
         prss,' GPa '//trim(ctmp)
    write(6,'(1x,a,6f10.3)') '  Stress tensor   =', &
         sth(1,1),sth(2,2),sth(3,3), &
         sth(2,3),sth(3,1),sth(1,2)
    write(6,*) ''

    if( tave.gt.100000d0 ) cftave = 'es12.4'
    tcpu = mpi_wtime() -tcpu0
    write(6,'(a,'//cfistp//','//cfetime//','//cftave &
         //',es13.4,2es11.3)') &
         " istp,etime,temp,epot,vol,prss=" &
         ,istp,tcpu,tave,epot,vol,prss
  endif

  call sanity_check(ekin,epot,stnsr,tave,myid_md,mpi_md_world)

!!$  print *,'Time at 3 = ',mpi_wtime() -tcpu0
!.....output initial configuration including epi, eki, and strs
!      write(cnum(1:4),'(i4.4)') 0
  write(cnum,'(i0)') 0
  tmp = mpi_wtime()
  call space_comp(ntot0,tagtot,rtot,vtot,atot,stot, &
       ekitot,epitot,auxtot)
  call accum_time('space_comp',mpi_wtime()-tmp)
  if( ifpmd.gt.0 .and. myid_md.eq.0 ) then
    if( ifsort.gt.0 ) then
      tmp = mpi_wtime()
      call sort_by_tag(ntot,tagtot,rtot,vtot &
           ,atot,ekitot,epitot,stot,auxtot,naux,ifsort)
      call accum_time('sort_by_tag',mpi_wtime()-tmp)
    endif
    tmp = mpi_wtime()
    if( ifpmd.eq.1 ) then  ! pmd format
      if( trim(ciofmt).eq.'bin' .or. trim(ciofmt).eq.'binary' ) &
           then
        call write_pmdtot_bin(20,"pmd_"//trim(cnum),ntot,hunit,h, &
             tagtot,rtot,vtot)
      elseif( trim(ciofmt).eq.'ascii' ) then
        if( lcomb_pos ) then
          call write_pmdtot_ascii(20,"pmdsnap",ntot,hunit,h, &
               tagtot,rtot,vtot,atot,epot,ekin,sth,.false.,0)
        else
          call write_pmdtot_ascii(20,"pmd_"//trim(cnum),ntot,hunit,h, &
               tagtot,rtot,vtot,atot,epot,ekin,sth,.false.,0)
        endif
      endif
    else if( ifpmd.eq.2 ) then ! LAMMPS-dump format
      if( lcomb_pos ) then  ! if combined, filename is dump.
        call write_dump(20,'dump',ntot,hunit,h,tagtot, &
             rtot,vtot,atot,stot,ekitot,epitot,naux,auxtot,0)
      else  ! otherwise, filename contains timestep
        call write_dump(20,'dump_'//trim(cnum),ntot,hunit,h,tagtot, &
             rtot,vtot,atot,stot,ekitot,epitot,naux,auxtot,0)
      endif
    endif
    call accum_time('write_xxx',mpi_wtime() -tmp)
  endif

!-----initialize the counter for output
  iocntpmd=0
  iocnterg=0

  if( myid_md.eq.0 ) then
    if( nerg.gt.0 .and. iprint.ge.ipl_basic ) then
!.....write out energies
      open(ioerg,file="out.erg",status='replace')
      write(ioerg,'(a)') '# 1:istp, 2:simtime[fs],' &
           //'   3:etot[eV],  4:ekin,' &
           //'  5:epot,  6:temp[K],  7:vol[Ang^3],  8:pressure[GPa]'
      write(ioerg,'(a,es16.7e3,a)') '#  Epot0 =',epot0,' [eV]'
      if( tave.gt.10000d0) cftave = 'es12.4'
      write(ioerg,'('//cfistp//','//cfstime//',3es16.7e3' &
           //','//cftave//',2es16.7e3)') istp &
           ,simtime,ekin+epot0,ekin,epot0,tave,vol,prss
      call flush(ioerg)
!.....Write stress components
      open(iostrs,file="out.strs",status='replace')
      write(iostrs,'(a)') '# 1:istp, 2:simtime[fs],' &
           //'  3:sxx[GPa],  4:syy,  5:szz,  6:syz,  7:sxz,  8:sxy'
      call flush(iostrs)
    endif
!c.....write out temperatures
!        open(iotemp,file='out.temperature',status='replace')
!        write(iotemp,'(a)') '# istp, temperature[0-9]'
!        ediff(1:9)= 0d0
!        write(iotemp,'(i10,18es16.7e3)') istp,temps(1:9),ediff(1:9)
!        call flush(iotemp)

!.....write out tensile forces
    if(trim(czload_type).eq.'atoms' .or. trim(czload_type).eq.'box' &
         .or. trim(czload_type).eq.'shear' ) then
      open(iozload,file='out.zload',status='replace')
      write(iozload,'(a)') '#  istp,   strnow,      ftop,     fbot'
    endif
  endif

  if( trim(czload_type).eq.'atoms' .or. &
       trim(czload_type).eq.'shear' ) then
    call get_forces_on_base(natm,ra,aa,tag,h,ftop &
         ,fbot,sorg,myid_md,mpi_md_world,iprint,czload_type)
    if( myid_md.eq.0 ) then
      write(iozload,'(i8,3es15.7)') 0,0.0,ftop,fbot
      call flush(iozload)
    endif
  endif

!.....temperature distribution along x
  if( ltdst ) then
    if( myid_md.eq.0 ) then
      open(iotdst,file='out.temp-dist')
      write(iotdst,'(a)') '#   x-pos,   temp   '
    endif
    call calc_temp_dist(iotdst,ntdst,tdst,nadst,natm,ra,eki &
         ,istp,nouterg,myid_md,mpi_md_world,sorg)
  endif

!.....Set al
  al(1)= h(1,1,0)
  al(2)= h(2,2,0)
  al(3)= h(3,3,0)

  if( lflux ) call accum_lflux(namax,natm,h,ra,va,aux(iaux_clr,:),istp,dt &
       ,myid_md,mpi_md_world,nxyz)
  if( lpdens ) call accum_pdens(namax,natm,tag,ra,sorg)

!!$  print *,'Time at 4 = ',mpi_wtime() -tcpu0
  i_conv = 0
  lconverged = .false.
!-----velocity-Verlet loop starts---------------------------------------
  do istp=1,nstp

!.....Metadynamics
    if( lmetaD ) then
      call update_metaD(istp,namax,natm,nsp,tag,ra,h,nnmax,lspr &
           ,myid_md,mpi_md_world,iprint)
    endif

!.....In case of isobaric MD, lstrs has to be always TRUE.
    if( index(cpctl,'beren').ne.0 .or. index(cpctl,'lange').ne.0 ) then
      lstrs = .true.
    else
      lstrs = lstrs0
    endif

!.....Variable time-step
    if( lvardt ) then
      call get_vmax(namax,natm,va,h,vmax,mpi_md_world)
      dt = min( dtmax, vardt_len/vmax )
      if( myid_md.eq.0 ) then
        if( iprint.ge.ipl_info ) then
          print '(a,f8.5,a)',' update dt to ',dt,' fs'
        endif
      endif
!.....Update dt-related values
      if( index(ctctl,'lange').ne.0 ) then
        tgmm = 1d0/trlx
        do itemp=1,ntemps
          if( ttgt(itemp).lt.0d0 ) then
            tfac(itemp)= -1d0
          else
            tfac(itemp)= dsqrt(2d0*tgmm*ttgt(itemp)/dt *k2ue)
          endif
        enddo
      else if( index(ctctl,'ttm').ne.0 ) then
        tmp = mpi_wtime()
        call set_inner_dt(dt,myid_md,iprint)
        call accum_time('ttm',mpi_wtime()-tmp)
      endif
    endif

!-------first kick of velocities
    do i=1,natm
      is = int(tag(i))
      va(1:3,i)=va(1:3,i) +aa(1:3,i)*fa2v(is)*dt
    enddo
    if( chgopt_method(1:4).eq.'xlag' ) call update_vauxq(aux(iaux_vq,:))
    if( index(cpctl,'lange').ne.0 ) call cvel_update_langevin(stnsr,h,mpi_md_world,1)

    if( ifdmp.eq.2 ) then
      call vfire(num_fire,alp0_fire,alp_fire,falp_fire,dtmax_fire &
           ,finc_fire,fdec_fire,nmin_fire &
           ,natm,va,aa,myid_md,mpi_md_world,dt,iprint)
!.....Restrict dt (use vardt_len here)
      call get_vmax(namax,natm,va,h,vmax,mpi_md_world)
      dt = min( dt, vardt_len/vmax )
    else if( ifdmp.eq.1 ) then
      do i=1,natm
        va(1:3,i) = va(1:3,i) *dmp
      enddo
      h(:,:,1) = h(:,:,1) *dmp
    endif

!-------multiply fmv
    do i=1,natm
      l= int(mod(tag(i)*10,10d0))
      va(1:3,i)=va(1:3,i) *fmv(1:3,l)
    enddo

    if( lclrchg ) then  ! special treatment for translational momentum
      call rm_trans_clrchg(natm,tag,va,am,mpi_md_world,myid_md,iprint)
    elseif( nrmtrans.gt.0 .and. mod(istp,nrmtrans).eq.0 ) then
      call rm_trans_motion(natm,tag,va,nspmax,am &
           ,mpi_md_world,myid_md,iprint)
    endif

!!$    print *,'Time at 5 = ',mpi_wtime() -tcpu0
!.....Update positions
    if( lconst ) then
      call update_const_pos(namax,natm,h,hi,tag,ra,va,dt,nspmax,am)
    else
      do i=1,natm
        ra(1:3,i)=ra(1:3,i) +(hi(1:3,1)*va(1,i) &
             +hi(1:3,2)*va(2,i) +hi(1:3,3)*va(3,i) )*dt
      enddo
    endif
    ltot_updated = .false.
    if( chgopt_method(1:4).eq.'xlag' ) &
         call update_auxq(aux(iaux_q,:),aux(iaux_vq,:))

    if( trim(czload_type).eq.'atoms' ) then
      call zload_atoms(natm,ra,tag,nstp,strfin,strnow &
           ,sorg,myid_md,mpi_md_world)
    else if( trim(czload_type).eq.'box' ) then
      call zload_box(natm,nstp,istp,dt,strfin,strnow,h,myid_md)
    else if( trim(czload_type).eq.'shear' ) then
      call shear_atoms(natm,ra,tag,nstp,strfin,strnow &
           ,sorg,myid_md,mpi_md_world)
    endif

    if( index(cpctl,'beren').ne.0 ) then
      call cell_update_berendsen(ah,h,lcellfix,lcell_updated)
    else if( index(cpctl,'lange').ne.0 ) then
      call cell_update_langevin(h,lcellfix,lcell_updated)
    endif

    if( trim(cdeform).ne.'none' ) then
      call apply_deform(h,dt,simtime)
      lcell_updated = .true.
    endif

    simtime = simtime +dt

!.....Reset matrices
    if( lcell_updated ) call boxmat(h,hi,ht,g,gi,gt,vol,sgm)

!.....Update rbufres considering change of vmax
    rbufres = rbufres -vmax*dt*2d0

!.....Update pair list and boundary atoms
!.....if making new pair list is needed.
    if( rbufres.le.0d0 .or. &
         (ifpmd.gt.0.and.mod(istp,noutpmd).eq.0) ) then
      if( iprint.ge.ipl_info .and. myid_md.eq.0 ) then
        print *,'Update boundary atoms and thus pair-list, too.'
        if( rbufres.le.0d0 ) then
          print *,'  because of rbufres < 0.'
        else
          print *,'  because of writing out positions.'
        endif
      endif
!.....Move atoms that cross the boundary
      tmp = mpi_wtime()
      call bamove()
      l1st = .false.
!.....Copy RA of boundary atoms
      call bacopy(.false.)
      call accum_time('ba_xxx',mpi_wtime()-tmp)
!.....Make pair list
      tmp = mpi_wtime()
!$acc update device(ra,h)
      call mk_lspr_para(l1st)
      call accum_time('lspr',mpi_wtime()-tmp)
      rbufres = rbuf
!!$!$acc update device(lspr)
    else
!.....Copy RA of boundary atoms determined by 'bacopy'
      tmp = mpi_wtime()
      call bacopy_fixed()
      call accum_time('ba_xxx',mpi_wtime()-tmp)
!$acc update device(ra,h)
    endif
!!$    print *,'Time at 6 = ',mpi_wtime() -tcpu0

    if(ifpmd.gt.0.and. mod(istp,noutpmd).eq.0 )then
      lstrs = lstrs0
    endif
!-------Calc forces
    tmp = mpi_wtime()
    call get_force(.false.,epot,stnsr)
    if( chgopt_method(1:4).eq.'xlag' ) then
      call get_aauxq(aux(iaux_chg,:),aux(iaux_q,:))
    endif
    call accum_time('get_force',mpi_wtime()-tmp)
    lcell_updated = .false.
    lstrs = lstrs0
!.....Descriptor spring
    if( ldspring ) call force_dspring(namax,natm,nnmax,lspr,rc,h,hi,tag,ra, &
         aa,aux(iaux_edsp,:), &
         nb,nbmax,lsb,nex,lsrc,myparity,nn,myid_md,mpi_md_world,iprint,.false.)
!.....Structure analysis
    if( trim(cstruct).eq.'CNA' &
         .and. mod(istp,istruct).eq.0 ) then
      call cna(namax,natm,h,ra,nb,nnmax,lspr,rc_struct)
    else if( trim(cstruct).eq.'a-CNA' &
         .and. mod(istp,istruct).eq.0 ) then
      call acna(namax,natm,h,ra,nb,nnmax,lspr,rc_struct)
    endif
!.....Color charge NEMD
    if( lclrchg ) call clrchg_force(namax,natm,tag,aa,aux,hi,specorder &
         ,myid_md,iprint)
!.....External force
    if( lextfrc ) call add_extfrc(natm,tag,aa,hi,specorder,myid_md,iprint)
!.....Force from metadynamics
    if( lmetaD ) then
      call force_metaD(istp,namax,natm,tag,ra,aa,h,hi,epot &
           ,nnmax,lspr,myid_md,mpi_md_world,iprint)
    endif
!.....Update some variables in constraints if needed
    if( lconst ) then
      call update_const(namax,natm,tag,ra,h,istp,nstp)
    endif
!.....Force_modify
    if( lrdcfrc ) call reduce_forces(namax,natm,aa,tag,ra &
         ,h,nnmax,lspr)
!!$    print *,'Time at 7 = ',mpi_wtime() -tcpu0

!.....Second kick of velocities
    if( index(ctctl,'lange').ne.0 ) then
      call vel_update_langevin(natm,tag,va,aa)
    else
      do i=1,natm
        is = int(tag(i))
        va(1:3,i)=va(1:3,i) +aa(1:3,i)*fa2v(is)*dt
      enddo
    endif
    if( chgopt_method(1:4).eq.'xlag' ) call update_vauxq(aux(iaux_vq,:))

!.....Constraints for velocities
    if( lconst ) then
      call update_const_vel(namax,natm,h,hi,tag,va,dt,nspmax,am)
    endif

!.....Calc kinetic energy
    vmaxold= vmax
    call get_ekin(namax,natm,va,tag,h,nspmax,fekin,ekin,eki,eks &
         ,vmax,mpi_md_world)
!!$    print *,'Time at 8 = ',mpi_wtime() -tcpu0

!.....Some thermostats come after get_ekin, since they require eks values
    if( index(ctctl,'beren').ne.0 ) then
      call vel_update_berendsen(natm,tag,va)
    else if( index(ctctl,'ttm').ne.0 ) then
      tmp = mpi_wtime()
      call assign_atom2cell(namax,natm,ra,sorg,boundary)
      call compute_nac(natm,myid_md,mpi_md_world,iprint)
      call calc_Ta(namax,natm,nspmax,h,tag,va,fmv,fekin &
           ,istp,myid_md,mpi_md_world,iprint)
      call langevin_ttm(namax,natm,va,aa,tag,am,h &
           ,nspmax,fa2v,fekin,ediff,dt,myid_md,mpi_md_world,iprint)
      call update_ttm(simtime,dt,natm,ra,h,sorg,myid_md,mpi_md_world,iprint)
      call non_reflecting_bc(natm,tag,ra,va,h,sorg,dt,nspmax,am,fa2v &
           ,myid_md,mpi_md_world,iprint)
      call te2tei(namax,natm,aux(iaux_tei,:))
      if( lcouple_3d1d ) then
        call set_3d1d_bc_pos(natm,ra,h,sorg,myid_md,mpi_md_world,iprint)
        call couple_3d1d(myid_md,mpi_md_world,iprint)
      endif
      call accum_time('ttm',mpi_wtime()-tmp)
    endif

    if( abs(pini-pfin).gt.0.1d0 ) then
      ptgt = ( pini +(pfin-pini)*istp/nstp ) *gpa2up
    endif
    call sa2stnsr(natm,strs,eki,stnsr,vol,mpi_md_world)
    if( index(cpctl,'beren').ne.0 ) then
      call cell_force_berendsen(stnsr,ah,mpi_md_world)
    else if( index(cpctl,'lange').ne.0 ) then
      call cvel_update_langevin(stnsr,h,mpi_md_world,2)
    endif
    sth(:,:) = stnsr(:,:) *up2gpa
    prss = (sth(1,1)+sth(2,2)+sth(3,3))/3

!.....temperature distribution along x
    if( ltdst ) call calc_temp_dist(iotdst,ntdst,tdst,nadst,natm,ra &
         ,eki,istp,nouterg,myid_md,mpi_md_world,sorg)

!c.....Eliminate translational motion
!        if( mod(istp,100).eq.0 ) then
!          call rm_trans_motion()
!        endif

!-------output energies
    if( mod(istp,nouterg).eq.0 .or. &
         (ifpmd.gt.0.and. mod(istp,noutpmd).eq.0) ) then
      iocnterg=iocnterg+1
      nave= 0
      tave= 0d0
      do itemp=1,ntemps
        if( ndof(itemp).le.0 ) cycle
        temps(itemp)= eks(itemp) *2d0 /fkb /max(ndof(itemp)-3,3)
        nave= nave +ndof(itemp)
        tave= tave +temps(itemp)*max(ndof(itemp)-3,3)
      enddo
      tave= tave/(nave-3)

      if( trim(czload_type).eq.'atoms' &
           .or. trim(czload_type).eq.'shear' ) then
        call get_forces_on_base(natm,ra,aa,tag,h,ftop &
             ,fbot,sorg,myid_md,mpi_md_world,iprint,czload_type)
      endif

      if( iprint.ne.0 ) then
        ediff(1:9)= ediff(1:9) /nouterg
        ediff0(1:9)= 0d0
        call mpi_reduce(ediff,ediff0,9 &
             ,mpi_double_precision,mpi_sum,0,mpi_md_world,ierr)
      endif

      if( myid_md.eq.0 .and. iprint.ge.ipl_basic ) then
!.....write energies
        if( tave.gt.10000d0) cftave = 'es12.4'
        write(ioerg,'('//cfistp//','//cfstime//',3es16.7e3' &
             //','//cftave//',2es16.7e3)') istp &
             ,simtime,ekin+epot,ekin,epot,tave,vol,prss
        call flush(ioerg)
!.....write stresses
        write(iostrs,'('//cfistp//','//cfstime//',6es12.3e3)') istp &
             ,simtime, sth(1,1), sth(2,2), sth(3,3), sth(2,3), sth(1,3), sth(1,2)
        call flush(iostrs)
!.....write temperature
        ediff(1:9)= 0d0

        if( trim(czload_type).eq.'atoms' .or. &
             trim(czload_type).eq.'shear' ) then
          write(iozload,'(i8,3es15.7)') istp,strnow,ftop,fbot
          call flush(iozload)
        endif
      endif   ! myid_md.eq.0

!---------output step, time, and temperature
      tcpu= mpi_wtime() -tcpu1
      if( myid_md.eq.0 .and. iprint.ne.0 ) then
        if( tave.gt.10000d0 ) cftave = 'es12.4'
        tcpu = mpi_wtime() -tcpu0
        write(6,'(a,'//cfistp//','//cfetime//','//cftave &
             //',es13.4,2es11.3)') &
             " istp,etime,temp,epot,vol,prss=" &
             ,istp,tcpu,tave,epot,vol,prss
        if( (index(cpctl,'beren').ne.0 .or. index(cpctl,'lange').ne.0 ) &
             .and. iprint.ne.0 ) then
          write(6,'(a)') ' Lattice vectors:' !,h(1:3,1:3,0)
          write(6,'(a,"[ ",3f12.3," ]")') '   a = ',h(1:3,1,0)
          write(6,'(a,"[ ",3f12.3," ]")') '   b = ',h(1:3,2,0)
          write(6,'(a,"[ ",3f12.3," ]")') '   c = ',h(1:3,3,0)

          write(6,'(a,6f10.4)') ' Stress (GPa):', &
               sth(1,1),sth(2,2),sth(3,3), &
               sth(2,3),sth(1,3),sth(1,2)
        endif
        call flush(6)
      endif

      if( index(ctctl,'ttm').ne.0 ) then
        call output_ttm(istp,simtime,myid_md,iprint)
        call output_energy_balance(istp,simtime,myid_md,iprint)
      endif

#ifdef __DISL__
!.....Output disl core pos
      call perf_disl_pos_by_pot(epith,natm,ra,h,epi,sorg &
           ,nodes_md,myid_md,mpi_md_world,iocnterg,21)
#endif
    endif ! energy output

    call sanity_check(ekin,epot,stnsr,tave,myid_md,mpi_md_world)

!.....check convergence criteria if it is dumped MD
    if( ifdmp.gt.0 .and. epot-epotp.le.0d0 .and. n_conv.gt.0 .and. &
         abs(epot-epotp).lt.eps_conv .and. istp.gt.minstp ) then
      i_conv = i_conv + 1
      if( i_conv.ge.n_conv ) then
        if( myid_md.eq.0 .and. iprint.ne.0 ) then
          write(6,'(/,a,i6,a,es10.3,a,i3,a)') &
               ' Damped MD converged with ',istp, &
               ' steps, since dE < ', &
               eps_conv,', ',n_conv,' times'
          write(6,'(a,3es20.10)') '  epot,epotp,epot-epotp = ' &
               ,epot,epotp,epot-epotp
        endif
        lconverged = .true.
      endif
    else
      epotp = epot
      i_conv = 0
    endif

!!$    if( lclrchg ) call clrchg_force(namax,natm,tag,aa,aux,hi,specorder &
!!$         ,myid_md,iprint)
    if( lflux ) call accum_lflux(namax,natm,h,ra,va,aux(iaux_clr,:),istp,dt &
         ,myid_md,mpi_md_world,nxyz)
    if( lpdens ) call accum_pdens(namax,natm,tag,ra,sorg)

!-------write the particle positions
    if(ifpmd.gt.0.and. &
         (mod(istp,noutpmd).eq.0.or.lconverged) )then
      tmp = mpi_wtime()
!---------decide pmd-file name
      iocntpmd=iocntpmd+1
      write(cnum,'(i0)') istp
      call space_comp(ntot0,tagtot,rtot,vtot,atot,stot, &
           ekitot,epitot,auxtot)
      call accum_time('space_comp',mpi_wtime()-tmp)
      ltot_updated = .true.
      if( ifpmd.gt.0 .and. myid_md.eq.0 ) then
        if( ifsort.gt.0 ) then
          tmp = mpi_wtime()
          call sort_by_tag(ntot0,tagtot,rtot,vtot &
               ,atot,ekitot,epitot,stot,auxtot,naux,ifsort)
          call accum_time('sort_by_tag',mpi_wtime()-tmp)
        endif
        tmp = mpi_wtime()
        if( ifpmd.eq.1 ) then  ! pmd format
          if( trim(ciofmt).eq.'bin' .or. trim(ciofmt).eq.'binary' ) &
               then
            call write_pmdtot_bin(20,"pmd_"//trim(cnum),ntot,hunit,h, &
                 tagtot,rtot,vtot)
          elseif( trim(ciofmt).eq.'ascii' ) then
            if( lcomb_pos ) then
              call write_pmdtot_ascii(20,"pmdsnap",ntot,hunit,h, &
                   tagtot,rtot,vtot,atot,epot,ekin,sth,.false.,istp)
            else
              call write_pmdtot_ascii(20,"pmd_"//trim(cnum),ntot,hunit,h, &
                   tagtot,rtot,vtot,atot,epot,ekin,sth,.false.,istp)
            endif
          endif
        else if( ifpmd.eq.2 ) then  ! LAMMPS-dump format
          if( lcomb_pos ) then
            call write_dump(20,'dump',ntot,hunit,h,tagtot, &
                 rtot,vtot,atot,stot,ekitot,epitot,naux,auxtot,istp)
          else
            call write_dump(20,'dump_'//trim(cnum),ntot,hunit,h,tagtot, &
                 rtot,vtot,atot,stot,ekitot,epitot,naux,auxtot,istp)
          endif
        endif
        call accum_time('write_xxx',mpi_wtime() -tmp)
      endif
    endif

    if( lconverged ) exit
!.....End of velocity-verlet loop    
  enddo

!.....Close out_pos file in the end, if combined out_pos.
  if( lcomb_pos ) close(20)

  if( .not. ltot_updated ) then
    tmp = mpi_wtime()
    call space_comp(ntot0,tagtot,rtot,vtot,atot,stot,ekitot,epitot, &
         auxtot)
    call accum_time('space_comp',mpi_wtime()-tmp)
    if( myid_md.eq.0 ) then
      tmp = mpi_wtime()
      call sort_by_tag(ntot0,tagtot,rtot,vtot &
           ,atot,ekitot,epitot,stot,auxtot,naux,ifsort)
      call accum_time('sort_by_tag',mpi_wtime()-tmp)
    endif
  endif

  nstp_done = istp
  tcpu2= mpi_wtime()

  if(myid_md.eq.0) then
    if( nerg.gt.0 .and. iprint.ge.ipl_basic ) then
      close(ioerg)
      close(iostrs)
    endif
!        close(iotemp)
    if(czload_type.eq.'atoms' .or. czload_type.eq.'box') then
      close(iozload)
    endif
    if( ltdst ) close(iotdst)
  endif

  tcpu= tcpu2 -tcpu1
  if( myid_md.eq.0 .and. iprint.ne. 0 ) then
    write(6,*) ''
    write(6,'(1x,a)') "Final values:"
    write(6,'(1x,a,f16.5,a,f10.3,a)') "  Kinetic energy  = ",ekin &
         ,' eV = ',ekin/ntot0,' eV/atom'
    write(6,'(1x,a,f16.5,a,f10.3,a)') "  Potential energy= ",epot &
         ,' eV = ',epot/ntot0,' eV/atom'
    nave= 0
    tave= 0d0
    do itemp=1,ntemps
      if( ndof(itemp).le.0 ) cycle
      temps(itemp)= eks(itemp) *2d0 /fkb /max(ndof(itemp)-3,3)
      nave= nave +ndof(itemp)
      if( lmultemps ) then
        write(6,'(1x,a,i1,a,f16.5,a)') "  Temperature ",itemp &
             ,"   = ",temps(itemp),' K'
      endif
      tave= tave +temps(itemp)*max(ndof(itemp)-3,3)
    enddo
    tave= tave/(nave-3)
    write(6,'(1x,a,f16.5,a)') "  Temperature     = ",tave,' K'
    if( prss.lt.0d0 ) then
      ctmp = '(tensile)'
    else
      ctmp = '(compressive)'
    endif
    write(6,'(1x,a,f16.5,a)') "  Pressure        = ", &
         prss,' GPa '//trim(ctmp)
    if( index(cpctl,'beren').ne.0 ) then
!!$      if( trim(cpctl).eq.'berendsen' .or. &
!!$         trim(cpctl).eq.'vc-berendsen' .or. &
!!$         trim(cpctl).eq.'vv-berendsen' ) then
      call cell_info(h)
    endif
    write(6,*) ''
    write(6,'(1x,a,2(2x,i0))') "Max num of neighbors during MD and nnmax = ",maxnn,nnmax
    write(6,'(1x,a,i0)') "Max num of boundary atoms during MD = ",maxnb
    write(6,*) ''
  endif

  hmat(:,:,:) = h(:,:,:)  ! Return h-matrix as hmat

!.....Return stnsr in human-friendly unit, GPa
  stnsr(:,:) = stnsr(:,:) *up2gpa

!.....Output metadynamics potential
  if( lmetaD ) then
    call write_metaD_potential(nstp,nsp,myid_md,iprint)
  endif

  if( lrdcfrc ) then
    call finalize_rdcfrc()
  endif

  if( ldspring ) call final_dspring(myid_md)

!.....deallocate all the arrays allocated in pmd_core
  if( ltdst ) then
    deallocate(tdst,nadst)
  endif
  deallocate(ra,va,aa,ra0,strs,tag,lspr &
       ,epi,eki,lsb,lsex)
  deallocate(aux)
end subroutine pmd_core
!=======================================================================
subroutine oneshot(hunit,hmat,ntot0,tagtot,rtot,vtot,atot,stot, &
     ekitot,epitot,auxtot,ekin,epot,stnsr,linit)
!
!  In case that only one shot force calculation is required,
!  especially called from fitpot.
!
  use pmdvars
  use force
  use pairlist,only: mk_lspr_para
  implicit none
  include "mpif.h"
  include "./params_unit.h"
  include "./const.h"
  integer,intent(in):: ntot0
  real(8),intent(in):: hunit,hmat(3,3,0:1)
  real(8),intent(in):: tagtot(ntot0),rtot(3,ntot0),vtot(3,ntot0)
  real(8),intent(out):: atot(3,ntot0),stot(3,3,ntot0),auxtot(naux,ntot0)
  real(8),intent(out):: ekitot(3,3,ntot0),epitot(ntot0),ekin,epot,stnsr(3,3)
  logical,intent(in):: linit

  integer:: i,ierr,is,nspl,iprm0
  real(8):: aai(3),epott
  logical:: l1st
  character(len=3):: csp

  ntot = ntot0
  h(:,:,:) = hmat(:,:,:)

  call set_domain_vars()
!.....perform space decomposition after reading atomic configuration
  if( iprint.ge.ipl_basic ) then
    write(6,'(a)',advance='no') ' Species order: '
    do is=1,nspmax
      csp = specorder(is)
      if( trim(csp).ne.'x' ) write(6,'(1x,3a)',advance='no') csp
    enddo
    print *,''
    write(6,'(a,i8)') ' Number of total atoms = ',ntot0
    write(6,'(a)') " Lattice vectors:"
    write(6,'(a,"[ ",3f12.3," ]")') '   a = ',h(1:3,1,0)
    write(6,'(a,"[ ",3f12.3," ]")') '   b = ',h(1:3,2,0)
    write(6,'(a,"[ ",3f12.3," ]")') '   c = ',h(1:3,3,0)
  endif
  call space_decomp(ntot0,tagtot,rtot,vtot,auxtot)

!.....Some conversions
  nsp= 0
  do i=1,natm
!-------species of atom-i
    nsp= max(int(tag(i)),nsp)
  enddo
!-----get total number of species
  nspl = nsp
  call mpi_allreduce(nspl,nsp,1,mpi_integer,mpi_max &
       ,mpi_md_world,ierr)
!-----setup
  call setup(nspmax,am,fekin,fa2v)
!-----set HI and SGM
  call boxmat(h,hi,ht,g,gi,gt,vol,sgm)
!-----ntset
  call ntset(myx,myy,myz,nx,ny,nz,nn,sv,myparity,anxi,anyi,anzi)

  tcpu1= mpi_wtime()

!!$  call init_force(linit)

!-----copy RA of boundary atoms
  call check_size_and_parallel(sgm,vol,rc,anxi,anyi,anzi &
       ,nx,ny,nz,myid_md)
  call bacopy(.true.)
!-----Make pair list
  l1st = .true.
!!$  call mk_lspr_para(namax,natm,nbmax,nb,nnmax,tag,ra,va,rc+rbuf &
!!$       ,h,hi,anxi,anyi,anzi,lspr,iprint,l1st)
  call mk_lspr_para(l1st)
  lstrs = .true.

  if( iprint.ge.ipl_info ) print *,'get_force...'
  call get_force(.true.,epot,stnsr)

!      print *,'one_shot: 07'
  if( iprint.ge.ipl_info ) print *,'sa2stnsr...'
  call sa2stnsr(natm,strs,eki,stnsr,vol,mpi_md_world)
  stnsr(:,:) = stnsr(:,:) *up2gpa

  if( iprint.ge.ipl_info ) print *,'space_comp...'
  call space_comp(ntot0,tagtot,rtot,vtot,atot,stot,ekitot,epitot, &
       auxtot)

!.....revert forces to the unit eV/A before going out 
  if( myid_md.eq.0 ) then
    do i=1,ntot0
      is= int(tagtot(i))
      aai(1:3)= h(1:3,1,0)*atot(1,i) &
           +h(1:3,2,0)*atot(2,i) &
           +h(1:3,3,0)*atot(3,i)
      atot(1:3,i) = aai(1:3)
    enddo
  endif

  return
end subroutine oneshot
!=======================================================================
subroutine oneshot4fitpot(hunit,hmat,ntot0,tagtot,rtot,vtot,atot,stot, &
     ekitot,epitot,auxtot,ekin,epot,stnsr,lcalcgrad,ndimp,maxisp, &
     gwe,gwf,gws,lematch,lfmatch,lsmatch)
!
!  In case that only one shot force calculation is required,
!  especially called from fitpot.
!
  use util,only: iauxof
  use pmdvars
  use force
  use Morse,only: gradw_Morse,gradw_vcMorse
  use Coulomb,only: gradw_Coulomb
  use linreg,only: gradw_linreg
  use DNN,only: gradw_DNN
  use pairlist,only: mk_lspr_para
  implicit none
  include "mpif.h"
  include "./params_unit.h"
  include "./const.h"
  integer,intent(in):: ntot0
  real(8),intent(in):: hunit,hmat(3,3,0:1)
  real(8),intent(in):: tagtot(ntot0),rtot(3,ntot0),vtot(3,ntot0)
  real(8),intent(out):: atot(3,ntot0),stot(3,3,ntot0),auxtot(naux,ntot0)
  real(8),intent(out):: ekitot(3,3,ntot0),epitot(ntot0),ekin,epot,stnsr(3,3)
  logical,intent(in):: lcalcgrad
  integer,intent(in):: ndimp,maxisp
  real(8),intent(inout):: gwe(ndimp),gwf(3,ndimp,ntot0),gws(6,ndimp)
  logical,intent(in):: lematch,lfmatch,lsmatch

  integer:: i,ierr,is,nspl,iprm0
  real(8):: aai(3),epott
  logical:: l1st
  character(len=3):: csp

  ntot = ntot0
  h(:,:,:) = hmat(:,:,:)

  call set_domain_vars()
!.....perform space decomposition after reading atomic configuration
  if( iprint.ge.ipl_basic ) then
    write(6,'(a)',advance='no') ' Species order: '
    do is=1,nspmax
      csp = specorder(is)
      if( trim(csp).ne.'x' ) write(6,'(1x,3a)',advance='no') csp
    enddo
    print *,''
    write(6,'(a,i8)') ' Number of total atoms = ',ntot0
    write(6,'(a)') " Lattice vectors:"
    write(6,'(a,"[ ",3f12.3," ]")') '   a = ',h(1:3,1,0)
    write(6,'(a,"[ ",3f12.3," ]")') '   b = ',h(1:3,2,0)
    write(6,'(a,"[ ",3f12.3," ]")') '   c = ',h(1:3,3,0)
  endif
  call boxmat(h,hi,ht,g,gi,gt,vol,sgm)
  call space_decomp(ntot0,tagtot,rtot,vtot,auxtot)

!.....Some conversions
  nsp= 0
  do i=1,natm
!-------species of atom-i
    nsp= max(int(tag(i)),nsp)
  enddo
!-----get total number of species
  nspl = nsp
  call mpi_allreduce(nspl,nsp,1,mpi_integer,mpi_max &
       ,mpi_md_world,ierr)
!-----setup
  call setup(nspmax,am,fekin,fa2v)
!-----set HI and SGM
  call boxmat(h,hi,ht,g,gi,gt,vol,sgm)
!-----ntset
  call ntset(myx,myy,myz,nx,ny,nz,nn,sv,myparity,anxi,anyi,anzi)

  tcpu1= mpi_wtime()

!!$  call init_force(.true.)

!-----copy RA of boundary atoms
  call check_size_and_parallel(sgm,vol,rc,anxi,anyi,anzi &
       ,nx,ny,nz,myid_md)
  call bacopy(.true.)
!-----Make pair list
  l1st = .true.
!!$  call mk_lspr_para(namax,natm,nbmax,nb,nnmax,tag,ra,va,rc+rbuf &
!!$       ,h,hi,anxi,anyi,anzi,lspr,iprint,l1st)
  call mk_lspr_para(l1st)
  lstrs = .true.

  if( .not.lcalcgrad ) then
    if( iprint.ge.ipl_basic ) print *,'get_force...'
    call get_force(.true.,epot,stnsr)
    if( iprint.ge.ipl_basic ) print '(a,es15.7)',' Potential energy = ',epot
  else  ! lcalcgrad = .true.
    if( iprint.ge.ipl_basic ) print *,'gradw_xxxx...'
    epot = 0d0
    gwe(1:ndimp) = 0d0
    gwf(1:3,1:ndimp,1:natm) = 0d0
    gws(1:6,1:ndimp) = 0d0
    if( use_force('Morse') &
         .and. use_force('screened_Coulomb') ) then
      iprm0 = 0
      call gradw_Coulomb(namax,natm,nb,tag,ra,aux(iauxof('chg'),:), &
           nnmax,h,rc,lspr,epott,iprint,ndimp,gwe,gwf,gws, &
           lematch,lfmatch,lsmatch,iprm0,myid_md,mpi_md_world,specorder)
      iprm0 = maxisp
      call gradw_Morse(namax,natm,tag, &
           ra,nnmax,h,rc,lspr,epott,iprint,ndimp,gwe,gwf,gws, &
           lematch,lfmatch,lsmatch,iprm0)
    else if( use_force('Morse') ) then
      iprm0 = 0
      call gradw_Morse(namax,natm,tag, &
           ra,nnmax,h,rc,lspr,epott,iprint,ndimp,gwe,gwf,gws, &
           lematch,lfmatch,lsmatch,iprm0)
    else if( use_force('linreg') ) then
      iprm0 = 0
      call gradw_linreg(namax,natm,tag,ra,nnmax,h,rc,lspr, &
           iprint,ndimp,gwe,gwf,gws,lematch,lfmatch,lsmatch,iprm0)
    else if( use_force('DNN') ) then
      iprm0 = 0
      call gradw_DNN(namax,natm,tag,ra,nnmax,h,rc,lspr, &
           iprint,ndimp,gwe,gwf,gws,lematch,lfmatch,lsmatch,iprm0)
    endif
!.....Derivative of stress should be divided by the cell volume
    gws(:,:) = gws(:,:) /vol
  endif

!      print *,'one_shot: 07'
  if( iprint.ge.ipl_basic ) print *,'sa2stnsr...'
  call sa2stnsr(natm,strs,eki,stnsr,vol,mpi_md_world)
  stnsr(:,:) = stnsr(:,:) *up2gpa

  if( iprint.ge.ipl_basic ) print *,'space_comp...'
  call space_comp(ntot0,tagtot,rtot,vtot,atot,stot,ekitot,epitot, &
       auxtot)

!.....revert forces to the unit eV/A before going out 
  if( myid_md.eq.0 ) then
    do i=1,ntot0
      is= int(tagtot(i))
      aai(1:3)= h(1:3,1,0)*atot(1,i) &
           +h(1:3,2,0)*atot(2,i) &
           +h(1:3,3,0)*atot(3,i)
      atot(1:3,i) = aai(1:3)
    enddo
  endif

  return
end subroutine oneshot4fitpot
!=======================================================================
subroutine min_core(hunit,hmat,ntot0,tagtot,rtot,vtot,atot,stot &
     ,ekitot,epitot,auxtot,epot,stnsr)
!
!  Minimization/relaxation routine.
!  Since MD and minimization could be very different, use different
!  subroutine for minimization.
!
!.....All the arguments are in pmdvars module
  use pmdio,only: write_pmdtot_ascii, write_pmdtot_bin, write_dump
  use util,only: iauxof, cell_info
  use pmdvars
  use force
  use vector,only: dot
  use random,only: box_muller
  use pmdmpi,only: nid2xyz,xyz2nid
  use util,only: itotOf, ifmvOf
  use time, only: sec2hms, accum_time
  use pairlist, only: mk_lspr_para
  use Coulomb,only: chgopt_method, update_auxq, update_vauxq, get_aauxq
  use isostat,only: setup_cell_min

  implicit none
  include "mpif.h"
  include "./params_unit.h"
  include "./const.h"
  integer,intent(inout):: ntot0
  real(8),intent(in):: hunit
  real(8),intent(inout):: tagtot(ntot0),rtot(3,ntot0),vtot(3,ntot0), &
       hmat(3,3,0:1)
  real(8),intent(out):: atot(3,ntot0),stot(3,3,ntot0), &
       ekitot(3,3,ntot0),epitot(ntot0),auxtot(naux,ntot0)
  real(8),intent(out):: epot,stnsr(3,3)

  integer:: ifmv,ierr,i_conv,i
  real(8):: tmp,tave,prss,epotp,sth(3,3),ekin
  logical:: l1st
  logical:: lconverged = .false.
  character:: cnum*128, ctmp*128
!.....Formats for output
  character(len=20):: cfistp  = 'i10' !or larger
  character(len=20):: cfstime = 'es18.10'
  character(len=20):: cfetime = 'f12.2' ! for elapsed time
  character(len=20):: cftave  = 'f12.2' ! or 'es12.4' for high-T

  if( nodes_md.ne.1 ) then
    write(6,'(a)') 'Error: CG minimization is not available in parallel.'
    write(6,'(a)') '  Use damping MD instead for large systems.'
    stop
  endif

  tcpu0= mpi_wtime()
  h(:,:,:) = hmat(:,:,:)  ! use pmdvars variable h instead of hmat
  ntot = ntot0

  call mpi_bcast(ntot,1,mpi_integer,0,mpi_md_world,ierr)

!-----parallel configuration
  nxyz= nx*ny*nz
  anxi= 1d0/nx
  anyi= 1d0/ny
  anzi= 1d0/nz
!-----error trap
  if(nodes_md.ne.nxyz) then
    write(6,'(a)') " error: nodes_md .ne. nxyz!!"
    write(6,'(a,3i5)') ' myid_md,nodes_md,nxyz=' &
         ,myid_md,nodes_md,nxyz
    print *,'nx,ny,nz =',nx,ny,nz
    call mpi_finalize(ierr)
    stop
  endif
!-----vector node indices: range [0:nx-1]
  call nid2xyz(myid_md,myx,myy,myz)
!-----reduced node origin
  sorg(1)= anxi*myx
  sorg(2)= anyi*myy
  sorg(3)= anzi*myz

!.....perform space decomposition after reading atomic configuration
  tmp = mpi_wtime()
  call space_decomp(ntot0,tagtot,rtot,vtot,auxtot)
  call accum_time('space_decomp',mpi_wtime()-tmp)
!.....Some conversions
  do i=1,natm
    ra(1:3,i)= ra(1:3,i) -sorg(1:3)
  enddo

!-----setup
  call setup(nspmax,am,fekin,fa2v)
!-----set HI and SGM
  call boxmat(h,hi,ht,g,gi,gt,vol,sgm)
!-----ntset
  call ntset(myx,myy,myz,nx,ny,nz,nn,sv,myparity,anxi,anyi,anzi)

!.....get_num_dof is called once in a MD run
!!$  call get_num_dof(natm,tag,fmv,ndof,myid_md,mpi_md_world,iprint)
  call get_num_dof()

!.....Set initial velocity as zero
  va(1:3,1:natm) = 0d0

  tcpu1= mpi_wtime()

!-----copy RA of boundary atoms
  call check_size_and_parallel(sgm,vol,rc,anxi,anyi,anzi &
       ,nx,ny,nz,myid_md)
  l1st = .true.

  if( myid_md.eq.0 ) then
    if( nerg.gt.0 .and. iprint.ge.ipl_basic ) then
!.....write out energies
      open(ioerg,file="out.erg",status='replace')
      write(ioerg,'(a)') '# 1:istp, 2:simtime[fs],' &
           //'   3:etot[eV],  4:ekin,' &
           //'  5:epot,  6:temp[K],  7:vol[Ang^3],  8:pressure[GPa]'
      write(ioerg,'(a,es16.7e3,a)') '#  Epot0 =',epot0,' [eV]'
      call flush(ioerg)
!.....Write stress components
      open(iostrs,file="out.strs",status='replace')
      write(ioerg,'(a)') '# 1:istp, 2:simtime[fs],' &
           //'  3:sxx[GPa],  4:syy,  5:szz,  6:syz,  7:sxz,  8:sxy'
      call flush(iostrs)
    endif
  endif

  i_conv = 0
  lconverged = .false.

!.....Minimization loop starts
  do istp=1,nstp

!.....In case of variable-cell/volume, lstrs has to be always TRUE.
    if( index(cpctl,'vv').ne.0 .or. index(cpctl,'vc').ne.0 ) then
      lstrs = .true.
    else
      lstrs = lstrs0
    endif

    tmp = mpi_wtime()
    call bacopy(.true.)
    call accum_time('ba_xxx',mpi_wtime()-tmp)
!-----Make pair list
    tmp = mpi_wtime()
!!$    call mk_lspr_para(namax,natm,nbmax,nb,nnmax,tag,ra,va,rc+rbuf &
!!$         ,h,hi,anxi,anyi,anzi,lspr,iprint,l1st)
    call mk_lspr_para(l1st)
    call accum_time('lspr',mpi_wtime()-tmp)

    if( iprint.gt.0 .and. myid_md.eq.0 ) then
      print '(/a,i5)', ' Max num of neighbors = ',maxnn
    endif

!.....Calc forces
    lstrs = lstrs0 .or. (index(cpctl,'beren').ne.0)
!.....Cell is new at the first call of get_force
    lcell_updated = .true.
    tmp = mpi_wtime()
    call get_force(.true.,epot0,stnsr)
    call accum_time('get_force',mpi_wtime()-tmp)
    lcell_updated = .false.
    lstrs = .false.
    epot= epot0
    epotp = 0d0

    call sa2stnsr(natm,strs,eki,stnsr,vol,mpi_md_world)
    call setup_cell_min(myid_md,iprint)

    simtime = 0d0

    if(myid_md.eq.0 .and. iprint.ne.0 ) then

!.....Human-friendly stress unit GPa
      sth(:,:) = stnsr(:,:) *up2gpa
      prss = (sth(1,1)+sth(2,2)+sth(3,3))/3
      if( prss.lt.0d0 ) then
        ctmp = '(tensile)'
      else
        ctmp = '(compressive)'
      endif
      write(6,'(1x,a,f16.5,a)') "  Pressure        = ", &
           prss,' GPa '//trim(ctmp)
      write(6,'(1x,a,6f10.3)') '  Stress tensor   =', &
           sth(1,1),sth(2,2),sth(3,3), &
           sth(2,3),sth(3,1),sth(1,2)
      write(6,*) ''

      tcpu = mpi_wtime() -tcpu0
      tave = 0d0
      write(6,'(a,'//cfistp//','//cfetime//','//cftave &
           //',es13.4,2es11.3)') &
           " istp,etime,temp,epot,vol,prss=" &
           ,istp,tcpu,tave,epot,vol,prss
    endif

!.....output initial configuration including epi, eki, and strs
    write(cnum,'(i0)') istp
    tmp = mpi_wtime()
    call space_comp(ntot0,tagtot,rtot,vtot,atot,stot, &
         ekitot,epitot,auxtot)
    call accum_time('space_comp',mpi_wtime()-tmp)
    if( ifpmd.gt.0 .and. myid_md.eq.0 ) then
      if( ifsort.gt.0 ) then
        tmp = mpi_wtime()
        call sort_by_tag(ntot,tagtot,rtot,vtot &
             ,atot,ekitot,epitot,stot,auxtot,naux,ifsort)
        call accum_time('sort_by_tag',mpi_wtime()-tmp)
      endif
      tmp = mpi_wtime()
      if( ifpmd.eq.1 ) then  ! pmd format
        if( trim(ciofmt).eq.'bin' .or. trim(ciofmt).eq.'binary' ) then
          call write_pmdtot_bin(20,"pmd_"//trim(cnum),ntot,hunit,h, &
               tagtot,rtot,vtot)
        elseif( trim(ciofmt).eq.'ascii' ) then
          sth(:,:) = stnsr(:,:)*up2gpa
          if( lcomb_pos ) then
            call write_pmdtot_ascii(20,"pmdsnap",ntot,hunit,h, &
                 tagtot,rtot,vtot,atot,epot,ekin,sth,.false.,istp)
          else            
            call write_pmdtot_ascii(20,"pmd_"//trim(cnum),ntot,hunit,h, &
                 tagtot,rtot,vtot,atot,epot,ekin,sth,.false.,istp)
          endif
        endif
      else if( ifpmd.eq.2 ) then ! LAMMPS-dump format
        if( lcomb_pos ) then
          call write_dump(20,'dump',ntot,hunit,h,tagtot, &
               rtot,vtot,atot,stot,ekitot,epitot,naux,auxtot,istp)
        else
          call write_dump(20,'dump_'//trim(cnum),ntot,hunit,h,tagtot, &
               rtot,vtot,atot,stot,ekitot,epitot,naux,auxtot,istp)
        endif
      endif
      call accum_time('write_xxx',mpi_wtime() -tmp)
    endif

    if( myid_md.eq.0 ) then
      if( nerg.gt.0 .and. iprint.ge.ipl_basic ) then
!.....write out energies
        tave = 0d0
        write(ioerg,'('//cfistp//','//cfstime//',3es16.7e3' &
             //','//cftave//',2es16.7e3)') istp &
             ,simtime,ekin+epot0,ekin,epot0,tave,vol,prss
        call flush(ioerg)
!.....write stresses
        write(iostrs,'('//cfistp//','//cfstime//',6es11.3e3)') istp &
             ,simtime, sth(1,1), sth(2,2), sth(3,3), sth(2,3), sth(1,3), sth(1,2)
        call flush(iostrs)
      endif
    endif

!!! TODO: complete CG minimization code here...    

  end do  ! istp

  nstp_done = istp
  tcpu2= mpi_wtime()

  if(myid_md.eq.0) then
    if( nerg.gt.0 .and. iprint.ge.ipl_basic ) then
      close(ioerg)
      close(iostrs)
    endif
  endif

  tcpu= tcpu2 -tcpu1
  if( myid_md.eq.0 .and. iprint.ne. 0 ) then
    write(6,*) ''
    write(6,'(1x,a)') "Final values:"
    write(6,'(1x,a,f16.5,a,f10.3,a)') "  Potential energy= ",epot &
         ,' eV = ',epot/ntot0,' eV/atom'
    if( prss.lt.0d0 ) then
      ctmp = '(tensile)'
    else
      ctmp = '(compressive)'
    endif
    write(6,'(1x,a,f16.5,a)') "  Pressure        = ", &
         prss,' GPa '//trim(ctmp)
    if( index(cpctl,'vv').ne.0 .or. &
         index(cpctl,'vc').ne.0  ) then
      call cell_info(h)
    endif
    write(6,*) ''
  endif
  
  hmat(:,:,:) = h(:,:,:)  ! Return h-matrix as hmat

!.....Return stnsr in human-friendly unit, GPa
  stnsr(:,:) = stnsr(:,:) *up2gpa
  
  return
end subroutine min_core
!=======================================================================
subroutine set_fmv(fmv)
!
! Set default fmv values which might be override
!
  implicit none
  real(8),intent(out):: fmv(3,0:9)

!-----set fmv(1:3,ifmv) to be multiplied to the velocities
  fmv(1:3,0)= (/ 0d0, 0d0, 0d0 /) ! fix
  fmv(1:3,1)= (/ 1d0, 1d0, 1d0 /) ! free move
  fmv(1:3,2)= (/ 1d0, 1d0, 1d0 /) ! free move
  fmv(1:3,3)= (/ 1d0, 1d0, 1d0 /) ! free move
  fmv(1:3,4)= (/ 1d0, 1d0, 1d0 /) ! free move
  fmv(1:3,5)= (/ 1d0, 1d0, 1d0 /) ! free move
  fmv(1:3,6)= (/ 1d0, 1d0, 1d0 /) ! free move
  fmv(1:3,7)= (/ 1d0, 1d0, 1d0 /) ! free move
  fmv(1:3,8)= (/ 1d0, 1d0, 1d0 /) ! free move
  fmv(1:3,9)= (/ 1d0, 1d0, 1d0 /) ! free move

end subroutine set_fmv
!=======================================================================
subroutine set_nsp(ntot,tagtot)
!
!  Set nsp from specorder.
!
  use pmdvars,only: nsp,nspmax,specorder,myid_md,mpi_md_world,iprint
  implicit none 
  integer,intent(in):: ntot
  real(8),intent(in):: tagtot(ntot)
  include './const.h'
  include 'mpif.h'

  integer:: i,ierr

  if( myid_md.eq.0 ) then
    nsp = 0
    do i=1,ntot
      nsp = max(int(tagtot(i)),nsp)
    enddo
  
    if( nsp.gt.nspmax ) then
      print *,'ERROR: nsp.gt.nspmax !!!'
      stop 1
    endif

    if( iprint.ge.ipl_basic ) then
      write(6,*) ''
      write(6,'(a,i0)') " Number of total atoms = ",ntot
      write(6,'(a,i0)') " Number of species     = ",nsp
    endif
  endif
  call mpi_bcast(nsp,1,mpi_integer,0,mpi_md_world,ierr)
  return
end subroutine set_nsp
!=======================================================================
subroutine set_domain_vars()
  use pmdvars
  implicit none
  
  nxyz = nx*ny*nz
  anxi= 1d0/nx
  anyi= 1d0/ny
  anzi= 1d0/nz
!-----vector node indices: range [0:nx-1]
  myx=myid_md/(ny*nz)
  myy=mod(myid_md/nz,ny)
  myz=mod(myid_md,nz)
!-----reduced node origin
  sorg(1)= anxi*myx
  sorg(2)= anyi*myy
  sorg(3)= anzi*myz

  return
end subroutine set_domain_vars
!=======================================================================
subroutine setup(nspmax,am,fekin,fa2v)
  implicit none
  include "params_unit.h"
  integer,intent(in):: nspmax
  real(8),intent(in):: am(nspmax)
  real(8),intent(out):: fekin(nspmax),fa2v(nspmax)

  integer:: is

!.....Factors for the following quantities
!.....Force: eV/Ang
!.....Acceleration: Ang/fs/fs
!.....Velocity: Ang/fs
!.....Position: Ang
!.....Energy: eV = 1.0/1.602e-19 [J = kg*m**2/s**2]
  do is=1,nspmax
    fa2v(is) = 0.5d0 *ev2j *(kg2ump *m2ang**2/s2fs**2) /am(is)
    fekin(is)= 0.5d0 *(am(is)*ump2kg) *ang2m**2 /fs2s**2 *j2ev
  enddo

end subroutine setup
!=======================================================================
subroutine boxmat(h,hi,ht,g,gi,gt,vol,sgm)
!-----------------------------------------------------------------------
!  setup matrices of MD-box
!    H:   MD-box matrix
!    HI:  inverse MD-box matrix
!    SGM: cofactor matrix
!-----------------------------------------------------------------------
  implicit none
  real(8),intent(in):: h(3,3,0:1)
  real(8),intent(out):: vol,sgm(3,3),hi(3,3),ht(3,3,0:1) &
       ,g(3,3,0:1),gi(3,3),gt(3,3,0:1)

  real(8):: hit(3,3)
  integer:: i,j,k,im,ip,jm,jp

!-----cofactor matrix, SGM
  do j=1,3
    jm=mod(j+1,3)+1
    jp=mod(j,  3)+1
    do i=1,3
      im=mod(i+1,3)+1
      ip=mod(i,  3)+1
      sgm(i,j)=h(ip,jp,0)*h(im,jm,0)-h(im,jp,0)*h(ip,jm,0)
    enddo
  enddo
!-----MD-box volume
  vol=h(1,1,0)*sgm(1,1)+h(2,1,0)*sgm(2,1)+h(3,1,0)*sgm(3,1)
  do j=1,3
    do i=1,3
      hit(i,j)= sgm(i,j)/vol
    enddo
  enddo
!-----transpose
  do j=1,3
    do i=1,3
      hi(i,j)= hit(j,i)
    enddo
  enddo

!.....Set transpose
  do j=1,3
    do i=1,3
      ht(i,j,0:1)= h(j,i,0:1)
    enddo
  enddo

!.....Set G-matrix
  g(1:3,1:3,0:1)= 0d0
  do j=1,3
    do i=1,3
      do k=1,3
        g(i,j,0)=g(i,j,0) +ht(i,k,0)*h(k,j,0)
        g(i,j,1)=g(i,j,1) +ht(i,k,1)*h(k,j,0) &
             +ht(i,k,0)*h(k,j,1)
      enddo
    enddo
  enddo
!.....Transpose of G
  do j=1,3
    do i=1,3
      gt(i,j,0:1)= g(j,i,0:1)
    enddo
  enddo
!.....Inverse of G
  call ludc_inv(3,g(1,1,0),gi)

  return
end subroutine boxmat
!=======================================================================
subroutine ntset(myx,myy,myz,nx,ny,nz,nn,sv,myparity,anxi,anyi,anzi)
!-----------------------------------------------------------------------
!  Preparation for network related properties
!-----------------------------------------------------------------------
  implicit none
  integer,intent(in):: myx,myy,myz,nx,ny,nz
  real(8),intent(in):: anxi,anyi,anzi
  integer,intent(out):: nn(6),myparity(3)
  real(8),intent(out):: sv(3,6)
  integer:: iv(3,6),ku,k1x,k1y,k1z

  iv(1:3,1)= (/ -1, 0, 0 /)
  iv(1:3,2)= (/  1, 0, 0 /)
  iv(1:3,3)= (/  0,-1, 0 /)
  iv(1:3,4)= (/  0, 1, 0 /)
  iv(1:3,5)= (/  0, 0,-1 /)
  iv(1:3,6)= (/  0, 0, 1 /)

  do ku=1,6
    k1x=mod(myx+iv(1,ku)+nx,nx)
    k1y=mod(myy+iv(2,ku)+ny,ny)
    k1z=mod(myz+iv(3,ku)+nz,nz)
!-------scalar neighbor ID, nn
    nn(ku)=k1x*(ny*nz)+k1y*nz+k1z
!-------shift vector for exchnaging ra
    sv(1,ku)= anxi*iv(1,ku)
    sv(2,ku)= anyi*iv(2,ku)
    sv(3,ku)= anzi*iv(3,ku)
  enddo

!-----Set up the node parity table
  if (nx.eq.1) then
    myparity(1)=2
  else
    myparity(1)=mod(myx,2)
  endif

  if (ny.eq.1) then
    myparity(2)=2
  else
    myparity(2)=mod(myy,2)
  endif

  if (nz.eq.1) then
    myparity(3)=2
  else
    myparity(3)=mod(myz,2)
  endif

  return
end subroutine ntset
!=======================================================================
subroutine get_ekin(namax,natm,va,tag,h,nspmax,fekin,ekin,eki,eks &
     ,vmax,mpi_md_world)
  use util,only: ithOf
  use time,only: accum_time
  implicit none 
  include "mpif.h"
  integer,intent(in):: namax,natm,mpi_md_world,nspmax
  real(8),intent(in):: va(3,namax),h(3,3),fekin(nspmax) &
       ,tag(namax)
  real(8),intent(out):: ekin,eki(3,3,namax),vmax,eks(nspmax)
!-----locals
  integer:: i,ierr,is,ixyz,jxyz,imax,igrp,itemp
  real(8):: ekinl,x,y,z,v(3),v2,vmaxl,eksl(nspmax),tmp

  ekinl=0d0
  eki(1:3,1:3,1:natm)= 0d0
  eksl(:)= 0d0
  vmaxl= 0d0

  igrp = 1  ! temperature category is Group-#1 (same as fmv)
  do i=1,natm
    is= int(tag(i))
    itemp = ithOf(tag(i),igrp)
    if( itemp.eq.0 ) cycle
!.....Tensor form eki
    do jxyz=1,3
      do ixyz=1,3
        eki(ixyz,jxyz,i)= va(ixyz,i)*va(jxyz,i)
      enddo
    enddo
    v2= eki(1,1,i) +eki(2,2,i) +eki(3,3,i)
    eki(1:3,1:3,i)=eki(1:3,1:3,i)*fekin(is)
    ekinl=ekinl +eki(1,1,i) +eki(2,2,i) +eki(3,3,i)
!.....ekin of each itemp
    eksl(itemp)= eksl(itemp) +eki(1,1,i) +eki(2,2,i) +eki(3,3,i)
!.....Find max speed
!!$    if( v2.gt.vmaxl ) imax=i
    vmaxl=max(vmaxl,v2)
  enddo
!      print *,'imax,vmax,ekin=',imax,vmaxl,ekinl
!      print *,'va(1:3,imax)=',va(1:3,imax)

  tmp = mpi_wtime()  
  ekin= 0d0
  call mpi_allreduce(ekinl,ekin,1,mpi_real8 &
       ,mpi_sum,mpi_md_world,ierr)
  eks(:)= 0d0
  call mpi_allreduce(eksl,eks,nspmax,mpi_real8 &
       ,mpi_sum,mpi_md_world,ierr)
  vmax= 0d0
  call mpi_allreduce(vmaxl,vmax,1,mpi_real8 &
       ,mpi_max,mpi_md_world,ierr)
  vmax=dsqrt(vmax)
  call accum_time('mpi_allreduce',mpi_wtime()-tmp)

end subroutine get_ekin
!=======================================================================
subroutine get_vmax(namax,natm,va,h,vmax,mpi_md_world)
  use time,only: accum_time
  implicit none 
  include "mpif.h"
  integer,intent(in):: namax,natm,mpi_md_world
  real(8),intent(in):: va(3,namax),h(3,3)
  real(8),intent(out):: vmax

  integer:: i,ierr
  real(8):: vx,vy,vz,v(3),v2,vmaxl,tmp

  vmaxl = 0d0
  do i=1,natm
!!$    vx = va(1,i)
!!$    vy = va(2,i)
!!$    vz = va(3,i)
!!$    v(1:3) = h(1:3,1)*vx +h(1:3,2)*vy +h(1:3,3)*vz
!!$    v2 = v(1)*v(1) +v(2)*v(2) +v(3)*v(3)
    v2 = va(1,i)**2 +va(2,i)**2 +va(3,i)**2
    vmaxl = max(vmaxl,v2)
  enddo
  tmp = mpi_wtime()
  vmax= 0d0
  call mpi_allreduce(vmaxl,vmax,1,mpi_real8 &
       ,mpi_max,mpi_md_world,ierr)
  vmax=dsqrt(vmax)
  call accum_time('mpi_allreduce',mpi_wtime()-tmp)
  return
end subroutine get_vmax
!=======================================================================
subroutine calc_temp_dist(iotdst,ntdst,tdst,nadst,natm,ra,eki &
     ,istp,nouterg,myid_md,mpi_md_world,sorg)
  implicit none
  include 'mpif.h'
  include './params_unit.h'
  integer,intent(in):: iotdst,ntdst,natm,istp,nouterg,myid_md &
       ,mpi_md_world
  real(8),intent(in):: ra(3,natm),eki(3,3,natm),sorg(3)
  integer,intent(out):: nadst(ntdst)
  real(8),intent(out):: tdst(ntdst)

  integer:: i,ix,ierr
  real(8):: dx,xi
  real(8),allocatable,save:: tdl(:)
  integer,allocatable,save:: ndl(:)
  logical,save:: l1st=.true.

  if( l1st ) then
    allocate(tdl(ntdst),ndl(ntdst))
    l1st= .false.
    tdst(1:ntdst)= 0d0
    nadst(1:ntdst)= 0
  endif

  dx= 1d0/ntdst
  do i=1,natm
    xi= ra(1,i)+sorg(1)
    ix= int(xi/dx)+1
    if( ix.le.0 ) ix= 1
    if( ix.gt.ntdst ) ix= ntdst
    tdst(ix)= tdst(ix) +(eki(1,1,i)+eki(2,2,i)+eki(3,3,i))
    nadst(ix)= nadst(ix) +1
  enddo

  if( mod(istp,nouterg).eq.0 ) then
!.....before writing to file, accumurate data of other nodes
    tdl(1:ntdst)= 0d0
    ndl(1:ntdst)= 0
    call mpi_reduce(nadst,ndl,ntdst &
         ,mpi_integer,mpi_sum,0,mpi_md_world,ierr)
    call mpi_reduce(tdst,tdl,ntdst &
         ,mpi_double_precision,mpi_sum,0,mpi_md_world,ierr)
    do i=1,ntdst
      if( ndl(i).ne.0 ) tdl(i)= tdl(i)*2d0/3/ndl(i)/fkb
      if( myid_md.eq.0 ) then
        write(iotdst,'(f10.5,f15.3)') ((i-1)+0.5)/ntdst, tdl(i)
      endif
    enddo
    if( myid_md.eq.0 ) write(iotdst,*) ''
!.....initialize tdst
    tdst(1:ntdst)= 0d0
    nadst(1:ntdst)= 0
  endif

end subroutine calc_temp_dist
!=======================================================================
subroutine get_num_dof()
  use pmdvars,only: natm,tag,fmv,nfmv,ndof, &
       myid_md,mpi_md_world,iprint
  use time,only: accum_time
  use util,only: ithOf
  implicit none
  include 'mpif.h'
  include "./const.h"
!!$  integer,intent(in):: natm,myid_md,mpi_md_world,iprint
!!$!.....TODO: hard coding number 9 should be replaced...
!!$  real(8),intent(in):: tag(natm),fmv(3,0:9)
!!$  integer,intent(out):: ndof(9)

  integer:: i,l,k,ndofl(9),ierr,igrp,ifmv
  real(8):: tmp
  real(8),parameter:: deps= 1d-12

  igrp = 1  ! Group-#1 for multiple-temperature (same as fmv)
  ndofl(1:9)= 0
  do i=1,natm
!!$    l= int(mod(tag(i)*10,10d0))
    ifmv= ithOf(tag(i),igrp)
    do k=1,3
      if( abs(fmv(k,ifmv)).lt.0.5d0 ) cycle
      ndofl(ifmv)= ndofl(ifmv) +1
    enddo
  enddo

  tmp = mpi_wtime()
  ndof(1:9)= 0
  call mpi_allreduce(ndofl,ndof,9,mpi_integer,mpi_sum &
       ,mpi_md_world,ierr)
  call accum_time('mpi_allreduce',mpi_wtime()-tmp)
  if(myid_md.eq.0 .and. iprint.ge.ipl_basic ) &
       write(6,'(/a,9(2x,i0,") ",i0))') &
       ' Degrees of freedom for each ifmv =',(i,ndof(i),i=1,nfmv)
  return
end subroutine get_num_dof
!=======================================================================
subroutine check_size_and_parallel(sgm,vol,rc,anxi,anyi,anzi &
     ,nx,ny,nz,myid_md)
  implicit none
  include 'mpif.h'
  integer,intent(in):: nx,ny,nz,myid_md
  real(8),intent(in):: sgm(3,3),vol,rc,anxi,anyi,anzi

  integer:: ierr
  real(8):: vala,valb,valc,rca,rcb,rcc

!-----calculate the cut-off lengths
  vala=dsqrt(sgm(1,1)**2+sgm(2,1)**2+sgm(3,1)**2)/vol
  valb=dsqrt(sgm(1,2)**2+sgm(2,2)**2+sgm(3,2)**2)/vol
  valc=dsqrt(sgm(1,3)**2+sgm(2,3)**2+sgm(3,3)**2)/vol
  rca=rc*vala
  rcb=rc*valb
  rcc=rc*valc

!.....check whether the size of the system is large enough
  if(  (rca.gt.anxi .and. nx.gt.1) .or. &
       (rcb.gt.anyi .and. ny.gt.1) .or. &
       (rcc.gt.anzi .and. nz.gt.1) ) then
    if( myid_md.le.0 ) then
      write(6,'(a)') " Error: about size and parallelization"
      write(6,'(a,2f15.3,i5)') "rca,anxi,nx = ",rca,anxi,nx
      write(6,'(a,2f15.3,i5)') "rcb,anyi,ny = ",rcb,anyi,ny
      write(6,'(a,2f15.3,i5)') "rcc,anzi,nz = ",rcc,anzi,nz
      write(6,'(a)') " * All the size (xyz) per node " &
           //"with more than 1 node should be larger " &
           //"than the cutoff radius."
      write(6,*) ""
    endif
    call mpi_finalize(ierr)
    stop
  endif

end subroutine check_size_and_parallel
!=======================================================================
subroutine bacopy(l1st)
!-----------------------------------------------------------------------
!  Exchanges boundary-atom data among neighbor nodes: tag and ra
!  Normal data passing and repeating the cell are mixed, and
!  automatically selected according to the size of the system.
!  
!  Note: parallelized to smaller than rcut should not happen.
!-----------------------------------------------------------------------
  use pmdvars
  use force
  use pmdmpi,only: nid2xyz
  use clrchg,only: lclrchg
  use time,only: accum_time
  implicit none
  include 'mpif.h'
  include './const.h'
  logical,intent(in):: l1st

!.....integer:: status(MPI_STATUS_SIZE)
  integer:: i,j,m,kd,kdd,kul,kuh,ku,ierr,iex,ix,iy,iz,itmp,istatus,iaux
  integer:: nav,maxna,maxb,inode,nsd,nrc,nbnew
  real(8):: xi(3),rcv(3),asgm,tmp
  logical,external:: bbd
  real(8),save,allocatable:: dbuf(:,:),dbufr(:,:)
  logical:: lshort(3)

  integer,save:: ndimbuf = 4

  if( l1st ) then
    if( allocated(dbuf) ) deallocate(dbuf,dbufr)
    ndimbuf = 4 +naux
    allocate(dbuf(ndimbuf,nbmax),dbufr(ndimbuf,nbmax))
  endif

  if( .not.allocated(dbuf) ) then
    allocate(dbuf(ndimbuf,nbmax),dbufr(ndimbuf,nbmax))
  else if( size(dbuf).lt.ndimbuf*nbmax ) then
    deallocate(dbuf,dbufr)
    allocate(dbuf(ndimbuf,nbmax),dbufr(ndimbuf,nbmax))
  endif

  call nid2xyz(myid_md,ix,iy,iz)

!-----reset the num of "received" boundary atoms
  nbnew=0

!-----calculate the cut-off lengths
  do kd=1,3
!.....ASGM is like a inverse of the lengh of cell vector
    asgm= dsqrt(sgm(1,kd)**2 +sgm(2,kd)**2 +sgm(3,kd)**2)
!.....RCV is like rc/(length of cell vector)
    rcv(kd)= rc*asgm/vol
!.....NEX > 1 if the cell vector is shorter than the rc,
!     meaning that the cell is so small that we need to consider multiple copies
!     along a certain direction.
    nex(kd)= int(rcv(kd)) +1
  enddo
  if( l1st .and. myid_md.eq.0 .and. iprint.ge.ipl_info ) then
    print '(a)',' bacopy info:'
    write(6,'(a,3f10.3)') '   rcv = ',rcv(1:3)
    write(6,'(a,3i10)')   '   nex = ',nex(1:3)
  endif

!.....Update nbmax if nex(i)>1,
  if( nex(1).gt.1 .or. nex(2).gt.1 .or.nex(3).gt.1 ) then
    maxb = ((2*nex(1)+1)*(2*nex(2)+1)*(2*nex(3)+1)-1)*natm
    if( maxb.gt.nbmax ) then
      if (myid_md.eq.0 .and. iprint.ne.0 ) then
        print *,'Updated namax and array since nbmax changed' &
             //' from ',nbmax,' to ',int(maxb*1.2), &
             ' because of nex(i)>1.'
      endif
      call realloc_namax_related(namax-nbmax,int(maxb*1.2))
    endif
  endif

!-----loop over x, y, & z directions
  do kd=1,3

!-------No. of to-be-copied atoms, LSB(0,)
    do kdd= -1,0
      lsb(0,2*kd+kdd)= 0
    enddo

!.....Dependent on the unit vector size, different treatment is required
    if( nex(kd).gt.1 ) then
      kul=2*kd-1
      kuh=2*kd
      do iex=-nex(kd),nex(kd)
        if( iex.eq.0 ) cycle
        if( iex.lt.0 ) then
          ku= 2*kd-1
        else
          ku= 2*kd
        endif
        do i=1,natm+nbnew
          lsb(0,ku)= lsb(0,ku) +1
          lsb(lsb(0,ku),ku)= i
          lsex(lsb(0,ku),ku)= iex
        enddo
      enddo
    else   ! long enough for normal boundary-atom copy
!-------Scan all the residents & copies
      do i=1,natm+nbnew
        xi(1:3)= ra(1:3,i)
!---------For low & high directions
        kul=2*kd-1
        kuh=2*kd
!---------Lower neighbor
        if (bbd(xi(1),xi(2),xi(3),rcv(1),rcv(2),rcv(3) &
             ,kul,anxi,anyi,anzi)) then
          lsb(0,kul)=lsb(0,kul)+1
          lsb(lsb(0,kul),kul)=i
        endif
!---------Higher neighbor
        if(bbd(xi(1),xi(2),xi(3),rcv(1),rcv(2),rcv(3) &
             ,kuh,anxi,anyi,anzi)) then
          lsb(0,kuh)=lsb(0,kuh)+1
          lsb(lsb(0,kuh),kuh)=i
        endif
      enddo
    endif    ! if(nex.gt.1)

!.....If BC is not periodic (p),
!.....number of to-be-sent atoms is set 0.
    if( boundary(kd:kd).ne.'p' )  then
      kul = 2*kd -1
      kuh = 2*kd
      if( kd.eq.1 ) then
        if( ix.eq.0 ) then
          lsb(0,kul) = 0
        endif
        if (ix.eq.nx-1 ) then
          lsb(0,kuh) = 0
        endif
      else if( kd.eq.2 ) then
        if( iy.eq.0 ) then
          lsb(0,kul) = 0
        endif
        if (iy.eq.ny-1 ) then
          lsb(0,kuh) = 0
        endif
      else if( kd.eq.3 ) then
        if( iz.eq.0 ) then
          lsb(0,kul) = 0
        endif
        if ( iz.eq.nz-1 ) then
          lsb(0,kuh) = 0
        endif
      endif
    endif

!-------Error trap
    itmp = nbnew
    do kdd= -1,0
      ku=2*kd+kdd
      itmp = itmp +lsb(0,ku)
    enddo
    tmp = mpi_wtime()
    call mpi_allreduce(itmp,maxb,1,mpi_integer,mpi_max, &
         mpi_md_world,ierr)
    call accum_time('mpi_allreduce',mpi_wtime()-tmp)
    if(maxb.gt.nbmax) then
      if( lrealloc ) then
        if (myid_md.eq.0 .and. iprint.ne.0 ) then
          print *,'Updated namax and array since nbmax changed' &
               //' from ',nbmax,' to ',int(maxb*1.2)
        endif
        call realloc_namax_related(namax-nbmax,int(maxb*1.2))
      else
        if( myid_md.eq.0 ) then
          print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          print *,'Exit pmd because maxb > nbmax.'
          print *,'If you do not want to stop, set allow_reallocation T in in.pmd.'
          print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        endif
        call mpi_finalize(ierr)
        stop
      endif
    endif

    if( nex(kd).gt.1 ) then
      do kdd= -1,0
        ku= 2*kd+kdd
        do i=1,lsb(0,ku)
          j= lsb(i,ku)
          iex= lsex(i,ku)
          ra(1:3,natm+nbnew+i)= ra(1:3,j)
          ra(kd,natm+nbnew+i)= ra(kd,natm+nbnew+i) +iex
          tag(natm+nbnew+i)= tag(j)
          aux(1:naux,natm+nbnew+i) = aux(1:naux,j)
        enddo
        nbnew= nbnew +lsb(0,ku)
      enddo
    else   ! long enough compared to rc for normal boundary-atom copy
      do kdd= -1,0
        ku=2*kd+kdd
        inode=nn(ku)
        nsd=lsb(0,ku)
        call mespasi(inode,myparity(kd),nsd,nrc,1,1,10 &
             ,mpi_md_world)
!---------Store the # of received boundary atoms in LSRC
        lsrc(ku)=nrc

!---------Exchange ra and tag
        do i=1,nsd
          j= lsb(i,ku)
          dbuf(1:3,i)= ra(1:3,j) -sv(1:3,ku)
          dbuf(4,i)  = tag(j)
          do iaux=1,naux
            dbuf(4+iaux,i) = aux(iaux,j)
          enddo
        enddo
        call mespasd(inode,myparity(kd),dbuf,dbufr,nsd*ndimbuf &
             ,nrc*ndimbuf,21,mpi_md_world)
        do i=1,nrc
          ra(1:3,natm+nbnew+i)= dbufr(1:3,i)
          tag(natm+nbnew+i)   = dbufr(4,i)
          do iaux=1,naux
            aux(iaux,natm+nbnew+i) = dbufr(4+iaux,i)
          enddo
        enddo

        call mpi_barrier(mpi_md_world,ierr)
!---------increase the # of received boundary atoms
        nbnew=nbnew+nrc

200     continue
      enddo ! kdd= -1,0
    endif

!-------Error trap
    nav=natm+nbnew
    tmp = mpi_wtime()
    call mpi_allreduce(nav,maxna,1,mpi_integer,mpi_max &
         ,mpi_md_world,ierr)
    call accum_time('mpi_allreduce',mpi_wtime()-tmp)
    if (maxna.gt.namax) then
      if (myid_md.eq.0) then
        write(*,*)'NAMAX overflowed at bacopy'
        write(*,*)'N+NB NAMAX = ',maxna,namax
      endif
      call mpi_finalize(ierr)
      stop
    endif

100 continue
  enddo ! kd=1,3

!-----num. of received boundary atoms
  nb=nbnew
  maxnb = max(maxnb,nb)

end subroutine bacopy
!=======================================================================
subroutine bacopy_fixed()
!-----------------------------------------------------------------------
!  Exchanges boundary-atom data among neighbor nodes: tag and ra
!  This does not search using position, just send & recv data of atoms
!    which were listed by 'bacopy'.
!  Different number of data are copied depending on whether 
!    using atomic charges or not.
!-----------------------------------------------------------------------
  use pmdvars
  use force
  use pmdmpi,only: nid2xyz
  use clrchg,only: lclrchg
  implicit none
  include 'mpif.h'

  integer:: i,j,m,kd,kdd,ku,ierr,iex,ix,iy,iz,iaux
  integer:: inode,nsd,nrc,nbnew
  real(8),save,allocatable:: dbuf(:,:),dbufr(:,:)
  logical,save:: l1st=.true.
  integer,save:: ndimbuf = 4

  if( l1st ) then
    if( allocated(dbuf) ) deallocate(dbuf,dbufr)
    ndimbuf = 4 +naux
    allocate(dbuf(ndimbuf,nbmax),dbufr(ndimbuf,nbmax))    
    l1st=.false.
  endif

  if( .not. allocated(dbuf) ) then
    allocate(dbuf(ndimbuf,nbmax),dbufr(ndimbuf,nbmax))
  else if( size(dbuf).lt.ndimbuf*nbmax ) then
    deallocate(dbuf,dbufr)
    allocate(dbuf(ndimbuf,nbmax),dbufr(ndimbuf,nbmax))
  endif

  call nid2xyz(myid_md,ix,iy,iz)

  nbnew= 0

!-----loop over x, y, & z directions
  do kd=1,3

    if( nex(kd).gt.1 ) then
      do kdd= -1,0
        ku= 2*kd+kdd
        do i=1,lsb(0,ku)
          j= lsb(i,ku)
          iex= lsex(i,ku)
          ra(1:3,natm+nbnew+i)= ra(1:3,j)
          ra(kd,natm+nbnew+i)= ra(kd,natm+nbnew+i) +iex
          tag(natm+nbnew+i)= tag(j)
          aux(1:naux,natm+nbnew+i) = aux(1:naux,j)
        enddo
        nbnew= nbnew +lsb(0,ku)
      enddo
    else  ! long enough compared to rc for normal boundary-atom copy
      do kdd= -1,0
        ku=2*kd+kdd
        inode=nn(ku)
        nsd=lsb(0,ku)
        nrc=lsrc(ku)

!---------Exchange ra and tag
        do i=1,nsd
          j= lsb(i,ku)
          dbuf(1:3,i)= ra(1:3,j) -sv(1:3,ku)
          dbuf(4,i)  = tag(j)
          do iaux=1,naux
            dbuf(4+iaux,i) = aux(iaux,j)
          enddo
        enddo
        call mespasd(inode,myparity(kd),dbuf,dbufr,nsd*ndimbuf &
             ,nrc*ndimbuf,21,mpi_md_world)
        do i=1,nrc
          ra(1:3,natm+nbnew+i)= dbufr(1:3,i)
          tag(natm+nbnew+i)   = dbufr(4,i)
          do iaux=1,naux
            aux(iaux,natm+nbnew+i) = dbufr(4+iaux,i)
          enddo
        enddo

        call mpi_barrier(mpi_md_world,ierr)
        nbnew=nbnew +nrc
200     continue
      enddo
    endif

100 continue
  enddo

end subroutine bacopy_fixed
!=======================================================================
subroutine bamove()
!-----------------------------------------------------------------------
!  Exchange atoms between neighbor nodes and atomic data as well.
!
!  MVQUE(0:NBMAX,6):
!    MVQUE(0,ku) is the # of to-be-moved atoms to neighbor ku;
!    MVQUE(i,ku) is the adress, in IS, of atom i
!-----------------------------------------------------------------------
  use pmdvars
  use force
  use pmdmpi,only: nid2xyz
  use clrchg,only: lclrchg
  use time,only: accum_time
  implicit none
  include 'mpif.h'

  integer:: i,j,m,ku,kd,kdd,kul,kuh,inode,nsd,nrc,ipt,ierr,is,ix,iy,iz,iaux
  integer:: mvque(0:nbmax,6),newim,maxa,itmp,mvql,mvqg
  real(8):: xi(3),tmp
  logical,external:: bmv
  real(8),save,allocatable:: dbuf(:,:),dbufr(:,:)
  logical,save:: l1st=.true.
  integer,save:: ndimbuf = 7

  if( l1st ) then
    ndimbuf = 7 +naux
    if( allocated(dbuf) ) deallocate(dbuf,dbufr)
    allocate(dbuf(ndimbuf,nbmax),dbufr(ndimbuf,nbmax))
    l1st=.false.
  endif

  if( .not.allocated(dbuf) ) then
    allocate(dbuf(ndimbuf,nbmax),dbufr(ndimbuf,nbmax))
  else if( size(dbuf).ne.ndimbuf*nbmax ) then
    deallocate(dbuf,dbufr)
    allocate(dbuf(ndimbuf,nbmax),dbufr(ndimbuf,nbmax))
  endif

  call nid2xyz(myid_md,ix,iy,iz)

!-----newim: num. of new immigrants
  newim= 0
  mvque(0,1:6)= 0

  do kd=1,3
    kul=2*kd-1
    kuh=2*kd

!-------num of to-be-moved atoms
    do i=1,natm+newim
      xi(1:3)= ra(1:3,i)
      is= int(tag(i))
      if (is.gt.0) then
        if (bmv(xi(1),xi(2),xi(3),kul,anxi,anyi,anzi)) then
          mvque(0,kul)=mvque(0,kul)+1
          mvque(mvque(0,kul),kul)=i
        else if (bmv(xi(1),xi(2),xi(3),kuh,anxi,anyi,anzi)) then
          mvque(0,kuh)=mvque(0,kuh)+1
          mvque(mvque(0,kuh),kuh)=i
        endif
      endif
    enddo

!.....Depending on boundary condition and cell position,
!.....number of to-be-moved atoms is set 0.
    if( boundary(kd:kd).ne.'p' )  then
      if( kd.eq.1 ) then
        if( ix.eq.0 ) then
          mvque(0,kul) = 0
        endif
        if (ix.eq.nx-1 ) then
          mvque(0,kuh) = 0
        endif
      else if( kd.eq.2 ) then
        if( iy.eq.0 ) then
          mvque(0,kul) = 0
        endif
        if (iy.eq.ny-1 ) then
          mvque(0,kuh) = 0
        endif
      else if( kd.eq.3 ) then
        if( iz.eq.0 ) then
          mvque(0,kul) = 0
        endif
        if ( iz.eq.nz-1 ) then
          mvque(0,kuh) = 0
        endif
      endif
    endif

!-------Error trap
    mvql = max(mvque(0,kul),mvque(0,kuh))
!!$    if (mvque(0,kul).gt.nbmax) then
!!$      jerr = mvque(0,kul)
!!$      print *,'Buffer overflowed at bamove node',myid_md
!!$      print *,'# in MVQUE=',mvque(0,kul)
!!$      stop
!!$    else if (mvque(0,kuh).gt.nbmax) then
!!$      jerr = mvque(0,kuh)
!!$      print *,'Buffer overflowed at bamove node',myid_md
!!$      print *,'# in MVQUE=',mvque(0,kuh)
!!$      stop
!!$    endif
    tmp = mpi_wtime()
    call mpi_allreduce(mvql,mvqg,1,mpi_integer,mpi_max, &
         mpi_md_world,ierr)
    call accum_time('mpi_allreduce',mpi_wtime()-tmp)
    if( mvqg.gt.nbmax ) then
      if( lrealloc ) then
        if( myid_md.eq.0 .and. iprint.ne.0 ) then
          print *,'Update nbmax from ',nbmax,' to ',nbmax*2
        endif
        call realloc_namax_related(namax-nbmax,nbmax*2)
      else
        if( myid_md.eq.0 ) then
          print *,'Exit pmd because mvqg > nbmax.'
          print *,'If you do not want to stop, set allow_reallocation T in in.pmd.'
        endif
        call mpi_finalize(ierr)
        stop
      endif
    endif
    
!.....Error trap for reallocation like bacopy
    itmp = natm + newim
    do kdd= -1,0
      ku = 2*kd +kdd
      itmp = itmp +mvque(0,ku)
    enddo
    tmp = mpi_wtime()
    call mpi_allreduce(itmp,maxa,1,mpi_integer,mpi_max, &
         mpi_md_world,ierr)
    call accum_time('mpi_allreduce',mpi_wtime()-tmp)
    if( maxa.gt.namax-nbmax ) then
      if( lrealloc ) then
        if( myid_md.eq.0 .and. iprint.ne.0 ) then
          print *,'Update namax from ',namax,' to ',maxa*2+nbmax
        endif
        call realloc_namax_related(maxa*2,nbmax)
      else
        if( myid_md.eq.0 ) then
          print *,'Exit pmd because maxa > nalmax.'
          print *,'If you do not want to stop, set allow_reallocation T in in.pmd.'
        endif
        call mpi_finalize(ierr)
        stop
      endif
    endif

    do kdd= -1,0

      ku=2*kd+kdd
      inode=nn(ku)
      nsd=mvque(0,ku)

      call mespasi(inode,myparity(kd),nsd,nrc,1,1,60, &
           mpi_md_world)

!---------move atom positions, RA
      do i=1,nsd
        j= mvque(i,ku)
        dbuf(1:3,i)= ra(1:3,j) -sv(1:3,ku)
        dbuf(4:6,i)= va(1:3,j)
        dbuf(7,i)  = tag(j)
!-----------Eliminate the record of a moved-out atom
        tag(j)= 0d0
        m = 7
        do iaux=1,naux
          dbuf(7+iaux,i) = aux(iaux,j)
        enddo
      enddo
      call mespasd(inode,myparity(kd),dbuf,dbufr,ndimbuf*nsd, &
           ndimbuf*nrc,71,mpi_md_world)
      do i=1,nrc
        ra(1:3,natm+newim+i)= dbufr(1:3,i)
        va(1:3,natm+newim+i)= dbufr(4:6,i)
        tag(natm+newim+i)   = dbufr(7,i)
        m = 7
        do iaux=1,naux
          aux(iaux,natm+newim+i) = dbufr(7+iaux,i)
        enddo
      enddo

      newim=newim+nrc
      call mpi_barrier(mpi_md_world,ierr)

    enddo

  enddo

!-----Compression
  ipt=0
  do i=1,natm+newim
    is= int(tag(i))
    if(is.ne.0) then
      ipt=ipt+1
      ra(1:3,ipt)= ra(1:3,i)
      va(1:3,ipt)= va(1:3,i)
      tag(ipt)   = tag(i)
      aux(1:naux,ipt) = aux(1:naux,i)
    endif
  enddo
!-----Update # of resident atoms
  natm=ipt

  return
end subroutine bamove
!=======================================================================
function bbd(xv,yv,zv,rcav,rcbv,rccv,ku,anxi,anyi,anzi)
!-----------------------------------------------------------------------
!  BBD = .true. if the coordinates are in the boundary to neighbor ku
!-----------------------------------------------------------------------
  implicit none
  real(8),intent(in):: xv,yv,zv,rcav,rcbv,rccv,anxi,anyi,anzi
  integer,intent(in):: ku
  logical:: bbd
  
  bbd = .false.
  if (ku.eq.1) then
    bbd = xv.lt.rcav
  else if (ku.eq.2) then
    bbd = anxi-rcav.lt.xv
  else if (ku.eq.3) then
    bbd = yv.lt.rcbv
  else if (ku.eq.4) then
    bbd = anyi-rcbv.lt.yv
  else if (ku.eq.5) then
    bbd = zv.lt.rccv
  else if (ku.eq.6) then
    bbd = anzi-rccv.lt.zv
  else
    write(6,'(a)') ' ERROR: BBD call is out of range'
    stop
  endif
  return
end function bbd
!=======================================================================
function bmv(xv,yv,zv,ku,anxi,anyi,anzi)
!-----------------------------------------------------------------------
!  BMV = .true. if the coordinates should belong to neighbor ku
!------------------------------------------------------------------------
  implicit none
  real(8),intent(in):: xv,yv,zv,anxi,anyi,anzi
  integer,intent(in):: ku
  
  logical:: bmv
  bmv = .false.
  if (ku.eq.1) then
    bmv = xv.lt.0d0
  else if (ku.eq.2) then
    bmv = anxi.lt.xv
  else if (ku.eq.3) then
    bmv = yv.lt.0d0
  else if (ku.eq.4) then
    bmv = anyi.lt.yv
  else if (ku.eq.5) then
    bmv = zv.lt.0d0
  else if (ku.eq.6) then
    bmv = anzi.lt.zv
  else
    write(6,'(a)')' ERROR: BMV call is out of range'
    stop
  endif
  return
end function bmv
!=======================================================================
subroutine sa2stnsr(natm,strs,eki,stnsr,vol,mpi_md_world)
!      
!  Take sum of atomic stresses to compute cell stress
!
  use time,only: accum_time
  implicit none
  include "mpif.h"
  integer,intent(in):: natm,mpi_md_world
  real(8),intent(in):: eki(3,3,natm),strs(3,3,natm),vol
  real(8),intent(out):: stnsr(3,3)

  integer:: i,ixyz,jxyz,ierr
  real(8):: stp(3,3),stk(3,3),stl(3,3),stg(3,3),tmp

  stp(1:3,1:3)= 0d0
  stk(1:3,1:3)= 0d0
  do i=1,natm
    do jxyz=1,3
      do ixyz=1,3
        stk(ixyz,jxyz)=stk(ixyz,jxyz) +2d0*eki(ixyz,jxyz,i)  ! eki in eV
        stp(ixyz,jxyz)=stp(ixyz,jxyz) +strs(ixyz,jxyz,i)  ! strs as rij*fij in eV
      enddo
    enddo
  enddo

  stl(1:3,1:3) = stk(1:3,1:3) +stp(1:3,1:3)
!!$  stl(1:3,1:3) = stp(1:3,1:3)
  stg(1:3,1:3)= 0d0
  tmp = mpi_wtime()
  call mpi_allreduce(stl,stg,9,mpi_real8,mpi_sum &
       ,mpi_md_world,ierr)
  call accum_time('mpi_allreduce',mpi_wtime()-tmp)
  stnsr(1:3,1:3) = stg(1:3,1:3)/vol
  return
end subroutine sa2stnsr
!=======================================================================
subroutine setv(h,hi,natm,tag,va,nspmax,am,tinit,dt)
  use random,only: box_muller
  implicit none
  include 'mpif.h'
  include 'params_unit.h'
  integer,intent(in):: natm,nspmax
  real(8),intent(in):: tag(natm),am(nspmax),tinit,dt &
       ,h(3,3,0:1),hi(3,3)
  real(8),intent(out):: va(3,natm)
  integer:: i,l,is
  real(8):: dseed,sumvx,sumvy,sumvz,tmp,facv(nspmax),vt(3)

!      facv(1:nspmax)=dsqrt(tinit*fkb*ev2j /(am(1:nspmax)*amu2kg))
!     &     *m2ang /s2fs
  facv(1:nspmax) = dsqrt(tinit *k2j &
       *(kg2ump* m2ang**2/s2fs**2) /am(1:nspmax))

!-----velocities in Maxwell-Boltzmann distribution
  do i=1,natm
    is= int(tag(i))
    do l=1,3
      va(l,i)=facv(is) *box_muller()
    enddo
!        write(6,'(i3,i5,i3,es24.14,3es15.7)') myid_md,i,is,tag(i)
!     &       ,va(1:3,i)
  enddo

!!$  do i=1,natm
!!$    vt(1:3) = va(1:3,i)
!!$    va(1:3,i) = hi(1:3,1)*vt(1) +hi(1:3,2)*vt(2) +hi(1:3,3)*vt(3)
!!$  enddo

end subroutine setv
!=======================================================================
subroutine rm_trans_motion(natm,tag,va,nspmax,am &
     ,mpi_md_world,myid_md,iprint)
  use time,only: accum_time
  implicit none
  include 'mpif.h'
  include "./const.h"
  integer,intent(in):: natm,nspmax,mpi_md_world,myid_md,iprint
  real(8),intent(in):: tag(natm),am(nspmax)
  real(8),intent(out):: va(3,natm)

  integer:: i,is,ierr
  real(8):: sumpx,sumpy,sumpz,amss,amtot,tmp,ttmp

!-----set center of mass motion to zero
  sumpx=0d0
  sumpy=0d0
  sumpz=0d0
  amtot=0d0
  do i=1,natm
    is= int(tag(i))
    amss= am(is)
    sumpx=sumpx+amss*va(1,i)
    sumpy=sumpy+amss*va(2,i)
    sumpz=sumpz+amss*va(3,i)
    amtot= amtot +amss
  enddo
  ttmp = mpi_wtime()
  tmp= sumpx
  call mpi_allreduce(tmp,sumpx,1,mpi_double_precision,mpi_sum &
       ,mpi_md_world,ierr)
  tmp= sumpy
  call mpi_allreduce(tmp,sumpy,1,mpi_double_precision,mpi_sum &
       ,mpi_md_world,ierr)
  tmp= sumpz
  call mpi_allreduce(tmp,sumpz,1,mpi_double_precision,mpi_sum &
       ,mpi_md_world,ierr)
  tmp= amtot
  call mpi_allreduce(tmp,amtot,1,mpi_double_precision,mpi_sum &
       ,mpi_md_world,ierr)
  call accum_time('mpi_allreduce',mpi_wtime()-ttmp)
  do i=1,natm
    va(1,i)=va(1,i)-sumpx/amtot
    va(2,i)=va(2,i)-sumpy/amtot
    va(3,i)=va(3,i)-sumpz/amtot
  enddo

  if( myid_md.eq.0 .and. iprint.ge.ipl_info ) then
    write(6,'(a,3es12.4)') ' sumpx,y,z/amtot=' &
         ,sumpx/amtot,sumpy/amtot,sumpz/amtot
  endif

end subroutine rm_trans_motion
!=======================================================================
subroutine vfire(num_fire,alp0_fire,alp_fire,falp_fire,dtmax_fire &
     ,finc_fire,fdec_fire,nmin_fire &
     ,natm,va,aa,myid_md,mpi_md_world,dt,iprint)
  use time,only: accum_time
  implicit none
  include 'mpif.h'
  include './const.h'
  integer,intent(in):: natm,myid_md,mpi_md_world,nmin_fire &
       ,iprint
  real(8),intent(in):: aa(3,natm),falp_fire,alp0_fire,dtmax_fire &
       ,finc_fire,fdec_fire
  integer,intent(inout):: num_fire
  real(8),intent(inout):: alp_fire,va(3,natm),dt

  integer:: i,ixyz,ierr
  real(8):: fdotv,vnorm,fnorm,tmp
  real(8):: fdotvl,vnorml,fnorml

  fdotvl = 0d0
  vnorml = 0d0
  fnorml = 0d0
  do i=1,natm
    do ixyz=1,3
      fdotvl = fdotvl +aa(ixyz,i)*va(ixyz,i)
      vnorml = vnorml +va(ixyz,i)*va(ixyz,i)
      fnorml = fnorml +aa(ixyz,i)*aa(ixyz,i)
    enddo
  enddo
  vnorm = 0d0
  fnorm = 0d0
  fdotv = 0d0
  tmp = mpi_wtime()
  call mpi_allreduce(vnorml,vnorm,1,mpi_real8 &
       ,mpi_sum,mpi_md_world,ierr)
  call mpi_allreduce(fnorml,fnorm,1,mpi_real8 &
       ,mpi_sum,mpi_md_world,ierr)
  call mpi_allreduce(fdotvl,fdotv,1,mpi_real8 &
       ,mpi_sum,mpi_md_world,ierr)
  call accum_time('mpi_allreduce',mpi_wtime()-tmp)
  vnorm = dsqrt(vnorm)
  fnorm = dsqrt(fnorm)
  do i=1,natm
    do ixyz=1,3
      va(ixyz,i) = (1d0 -alp_fire)*va(ixyz,i) &
           +alp_fire*aa(ixyz,i)/fnorm *vnorm
    enddo
  enddo
!      print '(a,2e14.5)','vnorm,fnorm =',vnorm,fnorm

  if( fdotv.gt.0d0  ) then
    num_fire = num_fire + 1
    if( num_fire.gt.nmin_fire ) then
      dt = min(dtmax_fire,dt*finc_fire)
      alp_fire = alp_fire *falp_fire
      if( iprint.ge.ipl_info .and. myid_md.eq.0 ) then
        write(6,'(a,f10.3,f10.5,es12.4)') ' dt,alp_fire,fdotv = ' &
             ,dt,alp_fire,fdotv
      endif
    endif
  else
    dt = min(dtmax_fire,dt*fdec_fire)
    va(1:3,1:natm) = 0d0
    alp_fire = alp0_fire
    num_fire = 0
    if( iprint.ge.ipl_info .and. myid_md.eq.0 ) then
      write(6,'(a,f10.3,f10.5,es12.4)') ' dt,alp_fire,fdotv = ' &
           ,dt,alp_fire,fdotv
    endif
  endif

!      call mpi_bcast(dt,1,mpi_real8,0,mpi_md_world,ierr)
!      call mpi_bcast(alp_fire,1,mpi_real8,0,mpi_md_world
!     &     ,ierr)
end subroutine vfire
!=======================================================================
subroutine space_decomp(ntot0,tagtot,rtot,vtot,auxtot)
!
!  Decompose the system and scatter atoms to every process.
!
  use pmdvars
  use clrchg,only: lclrchg
  implicit none
  include 'mpif.h'
  integer,intent(in):: ntot0
  real(8),intent(inout):: rtot(3,ntot0),tagtot(ntot0)
  real(8),intent(in):: vtot(3,ntot0)
  real(8),intent(in):: auxtot(naux,ntot0)
!!$  real(8),intent(in):: hunit,h(3,3,0:1)
!!$  real(8),intent(inout):: rtot(3,ntot0),tagtot(ntot0)
!!$  integer,intent(in):: ntot0,naux
!!$  real(8),intent(in):: vtot(3,ntot0),rcut,rbuf
!!$  real(8),intent(in):: auxtot(naux,ntot0)

  integer:: istat(mpi_status_size)
  integer:: i,j,ixyz,n,ierr,nacc,ir
  integer:: myxt,myyt,myzt,nmin
  real(8):: sxogt,syogt,szogt
  real(8):: t0
  logical,save:: l1st = .true.

  t0 = mpi_wtime()

!.....at 1st call
  if( .not. allocated(ra) ) then
    if( myid_md.eq.0 ) then
!.....wrap atoms into [0:1)
      do i =1,ntot0
        do ixyz=1,3
          if( rtot(ixyz,i).lt.0d0 ) then
            rtot(ixyz,i) = rtot(ixyz,i) +1d0
          else if( rtot(ixyz,i).ge.1d0 ) then
            rtot(ixyz,i) = rtot(ixyz,i) -1d0
          endif
        enddo
      enddo

!.....count max number of atoms in a node
      nalmax = 0
      nmin = 1000000000
      nacc = 0
      do ixyz=0,nxyz-1
        myxt = ixyz/(ny*nz)
        myyt = mod(ixyz/nz,ny)
        myzt = mod(ixyz,nz)
        sxogt = dble(myxt)/nx
        syogt = dble(myyt)/ny
        szogt = dble(myzt)/nz
        n = 0
        do i=1,ntot0
          if(  rtot(1,i).ge.sxogt .and. &
               rtot(1,i).le.sxogt+1d0/nx .and. &
               rtot(2,i).ge.syogt .and. &
               rtot(2,i).le.syogt+1d0/ny .and. &
               rtot(3,i).ge.szogt .and. &
               rtot(3,i).le.szogt+1d0/nz .and. &
               tagtot(i).gt.0d0) then
            n = n+1
!.....Set the tag negative if the atom is assigned to a certain cell
            tagtot(i) = -tagtot(i)
          endif
        enddo
        nalmax = max(nalmax,n)
        nmin = min(nmin,n)
        nacc = nacc + n
      enddo
      namax = max(int(nalmax*1.2),200)
!          nbmax = max(namax*27,nbmax)
      call estimate_nbmax(nalmax,h,nx,ny,nz,vol,rc,rbuf,nbmax,boundary)
      namax = namax +nbmax
      if( iprint.ne.0 .and. l1st ) then
        l1st = .false.
        print *,''
        print '(a)', ' space_decomp:'
        print '(a,2f6.3)','   rcut, rbuf = ',rc,rbuf
        write(6,'(a,i10)') '   Min number of local atoms = ',nmin
        write(6,'(a,i10)') '   Max number of local atoms = ',nalmax
        write(6,'(a,i10)')   '     nbmax = ',nbmax
        write(6,'(a,i10)')   '     namax = nalmax*1.2 + nbmax  = ' &
             ,namax
      endif
!.....Reset the tags positive
      do i=1,ntot0
        tagtot(i) = abs(tagtot(i))
      enddo
    endif ! myid.eq.0
    call mpi_bcast(namax,1,mpi_integer,0,mpi_md_world,ierr)
    call mpi_bcast(nbmax,1,mpi_integer,0,mpi_md_world,ierr)
    call alloc_namax_related()
    eki(1:3,1:3,1:namax) = 0d0
  endif  ! .not. allocated(ra)

!!$  call mpi_bcast(hunit,1,mpi_real8,0,mpi_md_world,ierr)
!!$  call mpi_bcast(h,9*2,mpi_real8,0,mpi_md_world,ierr)

  ixyz= 0
  if( myid_md.eq.0 ) then
    do ixyz=nxyz-1,0,-1
      myxt = ixyz/(ny*nz)
      myyt = mod(ixyz/nz,ny)
      myzt = mod(ixyz,nz)
      sxogt = dble(myxt)/nx
      syogt = dble(myyt)/ny
      szogt = dble(myzt)/nz
      natm = 0
      do i=1,ntot0
        if(  rtot(1,i).ge.sxogt .and. &
             rtot(1,i).le.sxogt+1d0/nx .and. &
             rtot(2,i).ge.syogt .and. &
             rtot(2,i).le.syogt+1d0/ny .and. &
             rtot(3,i).ge.szogt .and. &
             rtot(3,i).le.szogt+1d0/nz .and. &
             tagtot(i).gt.0d0 ) then
          natm = natm +1
          tag(natm)= tagtot(i)
          ra(1:3,natm)= rtot(1:3,i)
!!$          va(1:3,natm)= vtot(1:3,i)
!.....va is in real unit, whereas vtot is in scaled unit
          va(1:3,natm)= h(1:3,1,0)*vtot(1,i) +h(1:3,2,0)*vtot(2,i) &
               +h(1:3,3,0)*vtot(3,i)
          tagtot(i) = -tagtot(i)
          if( naux.gt.0 ) aux(1:naux,natm) = auxtot(1:naux,i)
        endif
      enddo
      if( ixyz.ne.0 ) then
        call mpi_send(natm,1,mpi_integer,ixyz,ixyz+1 &
             ,mpi_md_world,ierr)
        call mpi_send(tag,natm,mpi_real8,ixyz,ixyz+2 &
             ,mpi_md_world,ierr)
        call mpi_send(ra,3*natm,mpi_real8,ixyz,ixyz+3 &
             ,mpi_md_world,ierr)
        call mpi_send(va,3*natm,mpi_real8,ixyz,ixyz+4 &
             ,mpi_md_world,ierr)
        if( naux.gt.0 ) then
          call mpi_send(aux,natm*naux,mpi_real8,ixyz,ixyz+8 &
               ,mpi_md_world,ierr)
        endif
      endif
    enddo
!        write(6,'(a,f10.3)') ' time space_decomp = ',mpi_wtime() -t0
!.....Reset the tags positive
    do i=1,ntot0
      tagtot(i) = abs(tagtot(i))
    enddo
  else ! myid_md.ne.0
    call mpi_recv(natm,1,mpi_integer,0,myid_md+1 &
         ,mpi_md_world,istat,ierr)
    call mpi_recv(tag,natm,mpi_real8,0,myid_md+2 &
         ,mpi_md_world,istat,ierr)
    call mpi_recv(ra,3*natm,mpi_real8,0,myid_md+3 &
         ,mpi_md_world,istat,ierr)
    call mpi_recv(va,3*natm,mpi_real8,0,myid_md+4 &
         ,mpi_md_world,istat,ierr)
    if( naux.gt.0 ) then
      call mpi_recv(aux,natm*naux,mpi_real8,0,myid_md+8 &
           ,mpi_md_world,istat,ierr)
    endif
  endif
  call mpi_barrier(mpi_md_world,ierr)

end subroutine space_decomp
!=======================================================================
subroutine space_comp(ntot0,tagtot,rtot,vtot,atot,stot, &
     ekitot,epitot,auxtot)
!
!  Opposite to space_decomp, gather atoms from every process
!  to create the total system for output.
!
  use pmdvars
  use util,only: itotOf
  implicit none
  include 'mpif.h'
  integer,intent(in):: ntot0
  real(8),intent(out):: tagtot(ntot0),rtot(3,ntot0),vtot(3,ntot0) &
       ,atot(3,ntot0),epitot(ntot0),ekitot(3,3,ntot0),stot(3,3,ntot0)
  real(8),intent(out):: auxtot(naux,ntot0)
  integer,parameter:: nmpi = 10
  integer:: n0,ixyz,natmt,i,ierr,ntott
  integer:: istat(mpi_status_size),itag
  real(8):: t0
  real(8),allocatable,save:: ratmp(:,:),vatmp(:,:),aatmp(:,:)

  if( .not. allocated(ratmp) ) then
    allocate(ratmp(3,natm),vatmp(3,natm),aatmp(3,natm))
  else if( size(ratmp).lt.3*natm ) then
    deallocate(ratmp,vatmp,aatmp)
    allocate(ratmp(3,natm),vatmp(3,natm),aatmp(3,natm))
  endif

  if( myid_md.eq.0 ) then
    n0 = natm
    tagtot(1:natm) = tag(1:natm)
    rtot(1:3,1:natm) = ra(1:3,1:natm)
    do i=1,natm
      vtot(1:3,i) = hi(1:3,1)*va(1,i) +hi(1:3,2)*va(2,i) +hi(1:3,3)*va(3,i)
      atot(1:3,i) = hi(1:3,1)*aa(1,i) +hi(1:3,2)*aa(2,i) +hi(1:3,3)*aa(3,i)
    enddo
!!$    vtot(1:3,1:natm) = va(1:3,1:natm)
!!$    atot(1:3,1:natm) = aa(1:3,1:natm)
    epitot(1:natm) = epi(1:natm)
    ekitot(1:3,1:3,1:natm) = eki(1:3,1:3,1:natm)
    stot(1:3,1:3,1:natm) = strs(1:3,1:3,1:natm)
    if( naux.gt.0 ) then
      auxtot(:,1:natm) = aux(:,1:natm)
    endif
    ntott = natm
    n0 = n0 +1
    do ixyz=1,nxyz-1
      itag = ixyz*nmpi -nmpi
      call mpi_recv(natmt,1,mpi_integer,ixyz,itag &
           ,mpi_md_world,istat,ierr)
      if( natmt.eq.0 ) cycle
      ntott = ntott + natmt
      call mpi_recv(tagtot(n0),natmt,mpi_real8 &
           ,ixyz,itag+1,mpi_md_world,istat,ierr)
      call mpi_recv(rtot(1,n0),3*natmt,mpi_real8 &
           ,ixyz,itag+2,mpi_md_world,istat,ierr)
      call mpi_recv(vtot(1,n0),3*natmt,mpi_real8 &
           ,ixyz,itag+3,mpi_md_world,istat,ierr)
      call mpi_recv(epitot(n0),natmt,mpi_real8 &
           ,ixyz,itag+4,mpi_md_world,istat,ierr)
      call mpi_recv(ekitot(1,1,n0),3*3*natmt,mpi_real8 &
           ,ixyz,itag+5,mpi_md_world,istat,ierr)
      call mpi_recv(stot(1,1,n0),3*3*natmt,mpi_real8 &
           ,ixyz,itag+6,mpi_md_world,istat,ierr)
      call mpi_recv(atot(1,n0),3*natmt,mpi_real8 &
           ,ixyz,itag+7,mpi_md_world,istat,ierr)
      if( naux.gt.0 ) then
        call mpi_recv(auxtot(1,n0),natmt*naux,mpi_real8 &
             ,ixyz,itag+11,mpi_md_world,istat,ierr)
      endif
      n0 = n0 + natmt
    enddo
!.....Update ntot
    ntot = ntott
    if( ntot.gt.ntot0 ) then
      print *,'ERROR: ntot is greater than ntot0,' &
           //' which should not happen !'
      stop
    endif
  else ! myid_md.ne.0
    itag = myid_md*nmpi -nmpi
    call mpi_send(natm,1,mpi_integer,0,itag &
         ,mpi_md_world,ierr)
    if( natm.gt.0 ) then
      call mpi_send(tag,natm,mpi_real8,0,itag+1 &
           ,mpi_md_world,ierr)
!.....Positions should be shifted, velocities and forces should be converted to scaled ones
      do i=1,natm
        ratmp(1:3,i) = ra(1:3,i) + sorg(1:3)
        vatmp(1:3,i) = hi(1:3,1)*va(1,i) +hi(1:3,2)*va(2,i) +hi(1:3,3)*va(3,i)
        aatmp(1:3,i) = hi(1:3,1)*aa(1,i) +hi(1:3,2)*aa(2,i) +hi(1:3,3)*aa(3,i)
      enddo
      call mpi_send(ratmp,3*natm,mpi_real8,0,itag+2 &
           ,mpi_md_world,ierr)
      call mpi_send(vatmp,3*natm,mpi_real8,0,itag+3 &
           ,mpi_md_world,ierr)
      call mpi_send(epi,natm,mpi_real8,0,itag+4 &
           ,mpi_md_world,ierr)
      call mpi_send(eki,3*3*natm,mpi_real8,0,itag+5 &
           ,mpi_md_world,ierr)
      call mpi_send(strs,3*3*natm,mpi_real8,0,itag+6 &
           ,mpi_md_world,ierr)
      call mpi_send(aatmp,3*natm,mpi_real8,0,itag+7 &
           ,mpi_md_world,ierr)
      if( naux.gt.0 ) then
        call mpi_send(aux,natm*naux,mpi_real8,0,itag+11 &
             ,mpi_md_world,ierr)
      endif
    endif
  endif
  call mpi_barrier(mpi_md_world,ierr)
  call mpi_bcast(ntot,1,mpi_integer,0,mpi_md_world,ierr)

end subroutine space_comp
!=======================================================================
subroutine sort_by_tag(natm,tag,ra,va,aa,eki,epi,strs,aux,naux,ifsort)
!
!  Sort by tag for output.
!  - ifsort
!      1: quick sort
!      2: heap sort
!  
  use util,only: itotOf
  use time,only: accum_time
!!$  use force,only: luse_charge, luse_elec_temp
  implicit none
  include "mpif.h"
  integer,intent(in):: natm,ifsort
  real(8),intent(inout):: ra(3,natm),va(3,natm),aa(3,natm) &
       ,eki(3,3,natm),epi(natm),strs(3,3,natm),tag(natm)
  integer,intent(in):: naux
  real(8),intent(inout):: aux(naux,natm)

  integer,allocatable,save:: itag(:),idxarr(:)
  real(8),allocatable,save:: buf(:,:)
  integer:: i,j,k
  integer,save:: nsave = 0
  integer,save:: ndata
  real(8):: tmp
!!$  integer,external:: itotOf

!!$  ndata = 33
  ndata = 9

  if( .not. allocated(itag) .or. natm.gt.nsave ) then
    if( allocated(itag) ) deallocate(itag,idxarr,buf)
    nsave = natm
    allocate(itag(natm),idxarr(natm),buf(ndata,natm))
  endif

  do i=1,natm
    itag(i)= itotOf(tag(i))
  enddo

!!$  tmp = mpi_wtime()
  if( ifsort.eq.1 ) then
    call arg_heapsort_iarr(natm,natm,itag,idxarr)
  else  ! default 2
    call arg_qsort_iarr(natm,1,natm,itag,idxarr)
  endif
!!$  call accum_time('sorting',mpi_wtime()-tmp)

  buf(1:3,1:natm) = ra(1:3,1:natm)
  buf(4:6,1:natm) = va(1:3,1:natm)
  buf(7:9,1:natm) = aa(1:3,1:natm)
  do i=1,natm
    j = idxarr(i)
    ra(1:3,i) = buf(1:3,j)
    va(1:3,i) = buf(4:6,j)
    aa(1:3,i) = buf(7:9,j)
  enddo
  
  buf(1:3,1:natm) = strs(1:3,1,1:natm)
  buf(4:6,1:natm) = strs(1:3,2,1:natm)
  buf(7:9,1:natm) = strs(1:3,3,1:natm)
  do i=1,natm
    j = idxarr(i)
    strs(1:3,1,i) = buf(1:3,j)
    strs(1:3,2,i) = buf(4:6,j)
    strs(1:3,3,i) = buf(7:9,j)
  enddo
  
  buf(1:3,1:natm) = eki(1:3,1,1:natm)
  buf(4:6,1:natm) = eki(1:3,2,1:natm)
  buf(7:9,1:natm) = eki(1:3,3,1:natm)
  do i=1,natm
    j = idxarr(i)
    eki(1:3,1,i) = buf(1:3,j)
    eki(1:3,2,i) = buf(4:6,j)
    eki(1:3,3,i) = buf(7:9,j)
  enddo
  
  buf(1,1:natm) = tag(1:natm)
  buf(2,1:natm) = epi(1:natm)
!.....This part may cause cache mishits and harm the efficiency
  buf(2+1:2+naux,1:natm) = aux(1:naux,1:natm)
!!$  do k=1,naux
!!$    buf(2+k,1:natm) = aux(k,1:natm)
!!$  enddo
  do i=1,natm
    j = idxarr(i)
    tag(i) = buf(1,j)
    epi(i) = buf(2,j)
!!$    do k=1,naux
!!$      aux(k,i) = buf(2+k,j)
!!$    enddo
    aux(1:naux,i) = buf(2+1:2+naux,j)
  enddo

end subroutine sort_by_tag
!=======================================================================
subroutine error_mpi_stop(cerrmsg)
  use pmdvars,only: myid_md
  implicit none
  character(len=*):: cerrmsg

  integer:: ierr

  if(myid_md.eq.0) then
    write(6,'(a)') &
         '=================================================='
    write(6,'(a)') ' !!! ERROR !!!'
    write(6,'(a)') '   ' &
         //trim(cerrmsg)
    write(6,'(a)') &
         '=================================================='
  endif
  call mpi_finalize(ierr)
  stop
end subroutine error_mpi_stop
!=======================================================================
subroutine estimate_nbmax(nalmax,h,nx,ny,nz,vol,rcut,rbuf,nbmax,boundary)
  implicit none
  integer,intent(in):: nalmax,nx,ny,nz
  real(8),intent(in):: vol,rcut,rbuf,h(3,3)
  integer,intent(inout):: nbmax
  character(len=3),intent(in):: boundary

  integer:: nest
  real(8):: alx,aly,alz,area,dens_2A3,densl,dens,rc,volex,dvol

  rc = rcut +rbuf
  alx = dsqrt(h(1,1)**2 +h(2,1)**2 +h(3,1)**2)/nx +2*rc
  aly = dsqrt(h(1,2)**2 +h(2,2)**2 +h(3,2)**2)/ny +2*rc
  alz = dsqrt(h(1,3)**2 +h(2,3)**2 +h(3,3)**2)/nz +2*rc
!.....Estimated volume of extended system including immigrants region
  volex = alx*aly*alz
  dvol = volex -vol/(nx*ny*nz)
!.....Density cannot be obtained from nalmax because sometimes the system contains vacuum.
!.....Density = 0.1 means 1 atom in 10 A^3 which is dense enough for most cases.
  dens_2A3 = 0.1d0
!.....But for very dense cases, used the higher density.
  densl = dble(nalmax)/(vol/(nx*ny*nz))
  dens = max(dens_2A3,densl)
!!$  area = 0.0d0
!!$  if( boundary(1:1).eq.'p' ) area = area + aly*alz*2
!!$  if( boundary(2:2).eq.'p' ) area = area + alx*alz*2
!!$  if( boundary(3:3).eq.'p' ) area = area + alx*aly*2
!.....Estimated number of atoms in the margin region.
!!$  nest = density *area*(rcut+rbuf)
  nest = dens *dvol
!!$  print *, volex,dvol,dens_2A3,densl,dens,nest
  nbmax = max(nbmax,nest)

end subroutine estimate_nbmax
!=======================================================================
subroutine alloc_namax_related()
!     
!     Allocated arrays related to NAMAX.
!
  use pmdvars
  use memory, only: accum_mem
  implicit none

  integer:: mem

  if( allocated(ra) ) deallocate(ra)
  if( allocated(va) ) deallocate(va)
  if( allocated(aa) ) deallocate(aa)
  if( allocated(ra0) ) deallocate(ra0)
  if( allocated(strs) ) deallocate(strs)
  if( allocated(tag) ) deallocate(tag)
  if( allocated(epi) ) deallocate(epi)
  if( allocated(eki) ) deallocate(eki)
  if( allocated(aux) ) deallocate(aux)
  if( allocated(lsb) ) deallocate(lsb)
  if( allocated(lsex) ) deallocate(lsex)
  allocate(ra(3,namax),va(3,namax),aa(3,namax),ra0(3,namax) &
       ,strs(3,3,namax),tag(namax) &
       ,epi(namax),eki(3,3,namax) &
       ,lsb(0:nbmax,6),lsex(nbmax,6))
  allocate(aux(naux,namax))
  mem = 8*namax*(3 +3 +3 +3 +9 +1 +1 +9 +naux)
  mem = 4*6*(nbmax+1) +4*6*nbmax
  call accum_mem('pmd',mem)
  
  return
end subroutine alloc_namax_related
!=======================================================================
subroutine realloc_namax_related(newnalmax,newnbmax)
!     
!     Reallocated namax-related arrays everywhen namax needed to be
!     updated
!
  use pmdvars
  use memory, only: accum_mem
  implicit none
  integer,intent(in):: newnalmax,newnbmax

  integer:: ierr,newnamax,ndim,l,m,inc,mem
  real(8),allocatable:: arr(:)
  integer,allocatable:: iarr(:)

  if( .not. allocated(ra) ) then
    print *,'Error: Arrays are not allocated yet!'
    call mpi_finalize(ierr)
    stop
  endif

  newnamax = newnalmax + newnbmax
!.....No need of updating arrays
  if( newnamax .lt. namax ) return

!.....Since the arrays already have data,
!.....those data must be restored after reallocation.
!.....The hard coding here is too messy but I do not know how to avoid..

  mem = 0
!.....ra
  ndim = size(ra)
  allocate(arr(ndim))
  call copy_arr(ndim,ra,arr)
  deallocate(ra)
  allocate(ra(3,newnamax))
  call copy_arr(ndim,arr,ra)
  deallocate(arr)
  mem = mem -8*ndim +3*8*newnamax

!.....va
  ndim = size(va)
  allocate(arr(ndim))
  call copy_arr(ndim,va,arr)
  deallocate(va)
  allocate(va(3,newnamax))
  call copy_arr(ndim,arr,va)
  deallocate(arr)
  mem = mem -8*ndim +3*8*newnamax

!.....aa
  ndim = size(aa)
  allocate(arr(ndim))
  call copy_arr(ndim,aa,arr)
  deallocate(aa)
  allocate(aa(3,newnamax))
  call copy_arr(ndim,arr,aa)
  deallocate(arr)
  mem = mem -8*ndim +3*8*newnamax

!.....ra0
  ndim = size(ra0)
  allocate(arr(ndim))
  call copy_arr(ndim,ra0,arr)
  deallocate(ra0)
  allocate(ra0(3,newnamax))
  call copy_arr(ndim,arr,ra0)
  deallocate(arr)
  mem = mem -8*ndim +3*8*newnamax

!.....strs
  ndim = size(strs)
  allocate(arr(ndim))
  call copy_arr(ndim,strs,arr)
  deallocate(strs)
  allocate(strs(3,3,newnamax))
  call copy_arr(ndim,arr,strs)
  deallocate(arr)
  mem = mem -8*ndim +8*9*newnamax

!!$!.....stt
!!$  ndim = size(stt)
!!$  allocate(arr(ndim))
!!$  call copy_arr(ndim,stt,arr)
!!$  deallocate(stt)
!!$  allocate(stt(3,3,newnamax))
!!$  call copy_arr(ndim,arr,stt)
!!$  deallocate(arr)

!.....tag
  ndim = size(tag)
  allocate(arr(ndim))
  call copy_arr(ndim,tag,arr)
  deallocate(tag)
  allocate(tag(newnamax))
  call copy_arr(ndim,arr,tag)
  deallocate(arr)
  mem = mem -8*ndim +8*newnamax

!!$!.....lspr
!!$  ndim = size(lspr)
!!$  allocate(iarr(ndim))
!!$  call copy_iarr(ndim,lspr,iarr)
!!$  deallocate(lspr)
!!$  allocate(lspr(0:nnmax,newnamax))
!!$  call copy_iarr(ndim,iarr,lspr)
!!$  deallocate(iarr)
!!$  mem = mem -4*ndim +4*(nnmax+1)*newnamax

!.....epi
  ndim = size(epi)
  allocate(arr(ndim))
  call copy_arr(ndim,epi,arr)
  deallocate(epi)
  allocate(epi(newnamax))
  call copy_arr(ndim,arr,epi)
  deallocate(arr)
  mem = mem -8*ndim +8*newnamax

!.....eki
  ndim = size(eki)
  allocate(arr(ndim))
  call copy_arr(ndim,eki,arr)
  deallocate(eki)
  allocate(eki(3,3,newnamax))
  call copy_arr(ndim,arr,eki)
  deallocate(arr)
  mem = mem -8*ndim +8*9*newnamax

!!$!.....stp
!!$  ndim = size(stp)
!!$  allocate(arr(ndim))
!!$  call copy_arr(ndim,stp,arr)
!!$  deallocate(stp)
!!$  allocate(stp(3,3,newnamax))
!!$  call copy_arr(ndim,arr,stp)
!!$  deallocate(arr)

!!$!.....stn
!!$  ndim = size(stn)
!!$  allocate(arr(ndim))
!!$  call copy_arr(ndim,stn,arr)
!!$  deallocate(stn)
!!$  allocate(stn(3,3,newnamax))
!!$  call copy_arr(ndim,arr,stn)
!!$  deallocate(arr)

!.....aux
  ndim = size(aux)
  allocate(arr(ndim))
  call copy_arr(ndim,aux,arr)
  deallocate(aux)
  allocate(aux(naux,newnamax))
  call copy_arr(ndim,arr,aux)
  deallocate(arr)
  mem = mem -8*ndim +8*naux*newnamax

!.....lsb
  ndim = size(lsb)
  allocate(iarr(ndim))
  call copy_iarr(ndim,lsb,iarr)
  deallocate(lsb)
  allocate(lsb(0:newnbmax,6))
  mem = mem -4*ndim +4*6*(newnbmax+1)
  inc = 0
  do l=1,6
    do m=0,nbmax
      inc=inc+1
      lsb(m,l) = iarr(inc)
    enddo
  enddo
!      call copy_iarr(ndim,iarr,lsb)
  deallocate(iarr)

!.....lsex
  ndim = size(lsex)
  allocate(iarr(ndim))
  call copy_iarr(ndim,lsex,iarr)
  deallocate(lsex)
  allocate(lsex(newnbmax,6))
  mem = mem -4*ndim +4*6*newnbmax
  inc = 0
  do l=1,6
    do m=1,nbmax
      inc=inc+1
      lsex(m,l) = iarr(inc)
    enddo
  enddo
!      call copy_iarr(ndim,iarr,lsex)
  deallocate(iarr)

!.....Update NAMAX and NBMAX
  namax = newnamax
  nalmax= newnalmax
  nbmax = newnbmax

  call accum_mem('pmd',mem)

  return
end subroutine realloc_namax_related
!=======================================================================
subroutine copy_arr(ndim,srcarr,destarr)
  implicit none 
  integer,intent(in):: ndim
  real(8),intent(in):: srcarr(ndim)
  real(8),intent(out):: destarr(ndim)

  destarr(:) = srcarr(:)
  return
end subroutine copy_arr
!=======================================================================
subroutine copy_iarr(ndim,srcarr,destarr)
  implicit none 
  integer,intent(in):: ndim
  integer,intent(in):: srcarr(ndim)
  integer,intent(out):: destarr(ndim)

  destarr(:) = srcarr(:)
  return
end subroutine copy_iarr
!=======================================================================
subroutine sanity_check(ekin,epot,stnsr,tave,myid,mpi_world)
  use pmdvars,only: tlimit
  implicit none 
  real(8),intent(in):: ekin,epot,stnsr(3,3),tave
  integer,intent(in):: myid,mpi_world
  include "mpif.h"

  integer:: i,j,ierr
  character(len=128):: msg
  
  msg = ''

  if( epot*0d0 .ne. 0d0 ) then
    msg = 'ERROR: epot == NaN !'
    goto 10
  endif
  if( ekin*0d0 .ne. 0d0 ) then
    msg = 'ERROR: ekin == NaN !'
    goto 10
  endif
  do j=1,3
    do i=1,3
      if( stnsr(i,j)*0d0 .ne. 0d0 ) then
        msg = 'ERROR: stnsr(i,j) == NaN !'
        goto 10
      endif
    enddo
  enddo
  if( tave .gt. tlimit ) then
    msg = 'ERROR: too high temperature (tave > temperature_limit) !'
    goto 10
  endif

10 continue
  call mpi_bcast(msg,128,mpi_character,0,mpi_world,ierr)
  if( trim(msg).ne.'' ) then
    if( myid.eq.0 ) then
      print *,'Exit pmd, because of '//trim(msg)
    endif
    call mpi_finalize(ierr)
    stop
  endif
  return

end subroutine sanity_check
!=======================================================================
subroutine set_cauxarr()
  use pmdvars,only: cauxarr,naux, iaux_chg, iaux_q, iaux_vq, iaux_tei,&
       iaux_clr, ctctl, iaux_edsp
  use force,only: set_use_charge, set_use_elec_temp, &
       luse_charge, luse_elec_temp
  use Coulomb,only: chgopt_method
  use clrchg,only: lclrchg
  use dspring,only: ldspring

  integer:: inc

  call set_use_charge()
  call set_use_elec_temp()
  naux = 0
  if( luse_charge ) then
    naux = naux +1  ! chg
  endif
  if( chgopt_method(1:4).eq.'xlag' ) then
    naux = naux +2  ! auxq, vauxq
  endif
  if( luse_elec_temp .or. trim(ctctl).eq.'ttm' ) then
    naux = naux +1
  endif
  if( lclrchg ) then
    naux = naux +1
  endif
  if( ldspring ) then
    naux = naux +1
  endif
  if( allocated(cauxarr) ) then
    if( size(cauxarr).ne.naux ) deallocate(cauxarr)
  endif
  if( .not.allocated(cauxarr) ) allocate(cauxarr(naux))
  inc = 0
  if( luse_charge ) then
    inc = inc +1
    cauxarr(inc) = 'chg'
    iaux_chg = inc
  endif
  if( chgopt_method(1:4).eq.'xlag' ) then
    inc = inc +1
    cauxarr(inc) = 'auxq'
    iaux_q = inc
    inc = inc +1
    cauxarr(inc) = 'vauxq'
    iaux_vq = inc
  endif
  if( luse_elec_temp .or. trim(ctctl).eq.'ttm' ) then
    inc = inc +1
    cauxarr(inc) = 'tei'
    iaux_tei = inc
  endif
  if( lclrchg ) then
    inc = inc +1
    cauxarr(inc) = 'clr'
    iaux_clr = inc
  endif
  if( ldspring ) then
    inc = inc +1
    cauxarr(inc) = 'edsp'
    iaux_edsp = inc
  endif
  
end subroutine set_cauxarr
!=======================================================================
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd lib"
!     End:
