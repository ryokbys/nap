!-----------------------------------------------------------------------
!                     Last-modified: <2019-06-07 13:05:07 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
! Core subroutines/functions needed for pmd.
!-----------------------------------------------------------------------
subroutine pmd_core(hunit,h,ntot0,tagtot,rtot,vtot,atot,stot &
     ,ekitot,epitot,chgtot,chitot,maxstp,nerg,npmd &
     ,myid_md,mpi_md_world,nodes_md,nx,ny,nz,specorder &
     ,am,dt,vardt_len,ciofmt,ifpmd,rc,rbuf,rc1nn,ifdmp,dmp &
     ,minstp,tinit,tfin,ctctl,ttgt,trlx,ltdst,ntdst,nrmtrans,cpctl &
     ,stgt,ptgt,pini,pfin,srlx,stbeta,strfin,lstrs0,lcellfix,fmv &
     ,stnsr,epot,ekin,n_conv,ifcoulomb,czload_type,zskin_width &
     ,zshear_angle,eps_conv,ifsort,iprint,nstp_done,lvc,boundary &
     ,lmetaD,lconst,lrdcfrc,cstruct,istruct,cdeform,dhratio)
  use pmdio,only: write_pmdtot_ascii, write_pmdtot_bin, write_dump &
       ,namax,nbmax,nnmax,nspmax,alptot
  use pmdvars
  use zload
  use force
  use ttm,only: init_ttm,langevin_ttm,output_Te,t_ttm, &
       calc_Ta,update_Te,assign_atom2cell,output_energy_balance, &
       remove_ablated_atoms,set_inner_dt
  use pmdmpi,only: nid2xyz,xyz2nid
  use metadynamics,only: init_metaD,update_metaD,force_metaD &
       ,write_metaD_potential
  use constraints,only: init_const, update_const, update_const_vel &
       ,update_const_pos
  use rdcfrc,only: init_rdcfrc, reduce_forces, finalize_rdcfrc
  use structure,only: cna,acna
  use deform,only: init_deform, apply_deform
  use util,only: itotOf, ifmvOf
  implicit none
  include "mpif.h"
  include "./params_unit.h"
  integer,intent(in):: ntot0,maxstp,nerg,npmd,myid_md,mpi_md_world &
       ,ifpmd,ifdmp,minstp,ntdst,ifsort & !,numff &
       ,iprint,nodes_md,nx,ny,nz,n_conv,nrmtrans,istruct
  integer,intent(inout):: ifcoulomb
  integer,intent(out):: nstp_done
  real(8),intent(in):: hunit,tinit,tfin &
       ,trlx,srlx,stbeta,strfin,dmp,eps_conv,vardt_len &
       ,zskin_width,zshear_angle,dhratio(3,3)
  real(8),intent(inout):: tagtot(ntot0),rtot(3,ntot0),vtot(3,ntot0) &
       ,atot(3,ntot0),stot(3,3,ntot0),ekitot(3,3,ntot0) &
       ,epitot(ntot0),dt,rc,rbuf,rc1nn,h(3,3,0:1),stnsr(3,3) &
       ,fmv(3,0:9),epot,ekin,am(nspmax),stgt(3,3),ptgt,pini &
       ,pfin,ttgt(9),chgtot(ntot0),chitot(ntot0)
  character,intent(in):: ciofmt*6, cpctl*20, ctctl*20 &
       ,boundary*3
  character(len=3),intent(in):: specorder(nspmax) 
  character(len=*),intent(in):: czload_type,cstruct,cdeform
!      character(len=20),intent(in):: cffs(numff)
  logical,intent(in):: ltdst,lstrs0,lcellfix(3,3),lmetaD,lconst &
       ,lrdcfrc
  logical,intent(inout):: lvc

  integer:: i,j,k,l,m,n,ia,ib,is,ifmv,nave,nspl,i_conv,ierr,nxyz
  integer:: ihour,imin,isec
  real(8):: tmp,hscl(3),aai(3),ami,tave,vi(3),vl(3),epotp, &
       htmp(3,3),prss,dtmax,vmaxt,rbufres,tnow
  real(8),external:: box_muller,sprod
  logical:: lconverged = .false.
  logical:: lstrs = .false.
  logical:: lcell_updated = .true.
!.....FIRE variables
  real(8):: alp_fire,fnorm,vnorm,fdotv
  integer:: num_fire
!.....For tensile test of Al nanorod
  real(8):: strnow,ftop,fbot
!-----output file names
  character:: cnum*128, ctmp*128
!!$  integer,external:: itotOf,ifmvOf
  logical:: ltot_updated = .true.
!.....Formats for output
  character(len=20):: cfistp  = 'i10' !or larger
  character(len=20):: cfstime = 'es18.10'
  character(len=20):: cfetime = 'f12.2' ! for elapsed time
  character(len=20):: cftave  = 'f12.2' ! or 'es12.4' for high-T

  tcpu0= mpi_wtime()
  call initialize_pmdvars(nspmax)
  call calc_nfmv(ntot0,tagtot,myid_md,mpi_md_world)
  call set_use_charge()

  if( maxstp.le.0 ) then
    cfistp = 'i2'
  else
    write(cfistp(2:3),'(i2.2)') int(log10(dble(maxstp)))+2
  endif

!.....Variable time-step
  if( dt.lt.0d0 ) then
    dtmax = abs(dt)
    lvardt = .true.
    if( myid_md.eq.0 .and. iprint.gt.0 ) then
      print *,''
      print '(a,f8.3,a)',' Use variable time-step: dtmax = ' &
           ,dtmax,' fs'
    endif
    dt = dtmax
  endif

!.....multiply damping factor 1st to avoid doing every step
  if( ifdmp.eq.1 ) then
    fmv(1:3,0:9)= fmv(1:3,0:9) *dmp
  endif

!.....Set 0d0 to z of fmv(9) in case of z-loading
  if( trim(czload_type).eq.'atoms' .or. &
       trim(czload_type).eq.'box' ) then
    fmv(3,9) = 0d0
    if( myid_md.eq.0 .and. iprint.gt.0 ) then
      print *,''
      write(6,'(a)') 'Set fmv(3,9)=0d0, because czload_type=' &
           //trim(czload_type)
    endif
  else if( trim(czload_type).eq.'shear' ) then
    fmv(1,9) = 0d0
    fmv(3,9) = 0d0
    if( myid_md.eq.0 .and. iprint.gt.0 ) then
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

!-----output every these steps, NOUTERG, NOUTPMD
  if( nerg.ne.0 ) then
    nouterg = max(maxstp/nerg,1)
  else
    nouterg = maxstp +1
  endif
  if( npmd.ne.0 ) then
    noutpmd = max(maxstp/npmd,1)
  else
    noutpmd = maxstp +1
  endif
!.....perform space decomposition after reading atomic configuration
  call space_decomp(hunit,h,ntot0,tagtot,rtot,vtot,chgtot,chitot &
       ,myid_md,mpi_md_world,nx,ny,nz,nxyz,rc,rbuf,iprint)
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

!.....NEMD setting
  if( ltdst ) then
    if( mod(ntdst,nx).ne.0 ) &
         call error_mpi_stop('mod(ntdst,nx).ne.0')
    allocate(tdst(ntdst),nadst(ntdst))
  endif


!-----setup
  allocate(fekin(nspmax),fa2v(nspmax))
  call setup(nspmax,am,fekin,fa2v)
!-----set HI and SGM
  call boxmat(h,hi,ht,g,gi,gt,vol,sgm)
  if( myid_md.eq.0 ) then
    write(6,'(a,f0.2,a)') ' Cell volume = ',vol,' Ang^3'
  endif
!-----ntset
  call ntset(myx,myy,myz,nx,ny,nz,nn,sv,myparity,anxi,anyi,anzi)

!.....Deformation setup
  if( trim(cdeform).ne.'none' ) then
    call init_deform(maxstp,h,dhratio,myid_md,iprint)
  endif

!-----get total number of species
  nspl = 0
  do i=1,natm
    nspl = max(int(tag(i)),nspl)
  enddo
  call mpi_allreduce(nspl,nsp,1,mpi_integer,mpi_max &
       ,mpi_md_world,ierr)
  if( nsp.gt.nspmax ) then
    if( myid_md.eq.0 ) then
      print *,'ERROR: nsp.gt.nspmax !!!'
    endif
    call mpi_finalize(ierr)
    stop
  endif
  if(myid_md.eq.0 .and. iprint.ne.0 ) then
    write(6,*) ''
    write(6,'(a,i0)') " Number of total atoms = ",ntot0
    write(6,'(a,i0)') " Number of species     = ",nsp
  endif
!.....get_num_dof is called once in a MD run
  call get_num_dof(natm,tag,fmv,ndof,myid_md,mpi_md_world,iprint)

!.....Call init_metaD after setting namax and nsp
  if( lmetaD ) then
    call init_metaD(maxstp,nsp,namax,nodes_md &
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
    call set_zload_atoms(natm,ra,tag,h,fmv,sorg,strfin,maxstp &
         ,zskin_width,myid_md,mpi_md_world,iprint)
  else if( trim(czload_type).eq.'shear' ) then
    call set_shear(natm,ra,tag,h,fmv,sorg,strfin,maxstp &
         ,zskin_width,zshear_angle,myid_md,mpi_md_world,iprint)
  endif

!.....Set initial temperature if needed
  if( tinit.gt.1d-5 ) then
    call setv(h,natm,tag,va,nspmax,am,tinit,dt)
  elseif( abs(tinit).le.1d-5 ) then
    va(1:3,1:natm)= 0d0
  endif
  if( nrmtrans.ge.0 ) call rm_trans_motion(natm,tag,va,nspmax,am &
       ,mpi_md_world,myid_md,iprint)

  if( ifdmp.eq.2 ) then
    va(1:3,1:natm) = 0d0
  endif

  if( trim(ctctl).eq.'Langevin' ) then
    tgmm= 1d0/trlx
    do ifmv=1,9
      if( ttgt(ifmv).lt.0d0 ) then
        tfac(ifmv)= -1d0
      else
!.....TFAC should have [sqrt(ue/fs**2)] = [ue/Ang/sqrt(ump)] unit
!            tfac(ifmv)= sqrt(2d0*tgmm*(fkb*ev2j)*ttgt(ifmv)/dt)
!     &           *m2ang/s2fs
!            tfac(ifmv)= dsqrt(2d0*tgmm*fkb*ttgt(ifmv)/dt)
        tfac(ifmv)= dsqrt(2d0*tgmm*ttgt(ifmv)/dt *k2ue)
      endif
    enddo
    if( myid_md.eq.0 .and. iprint.gt.0 ) then
      print *,''
      print *,'Langevin thermostat parameters:'
      print '(a,f10.2)','   Relaxation time = ',trlx
      do ifmv=1,nfmv
        print '(a,i3,f10.2,es15.4)','   ifmv,ttgt,tfac = ' &
             ,ifmv,ttgt(ifmv),tfac(ifmv)
      enddo
    endif
  else if( trim(ctctl).eq.'ttm' ) then
    call init_ttm(namax,natm,h,dt,lvardt,myid_md,mpi_md_world &
         ,iprint)
  endif

  tcpu1= mpi_wtime()
  tcom = 0d0
  tlspr = 0d0
  tdump = 0d0

  call init_force(namax,natm,nsp,tag,chg,chi,myid_md,mpi_md_world, &
       iprint,h,rc,lvc,ifcoulomb,specorder)
!-----copy RA of boundary atoms
  call check_size_and_parallel(sgm,vol,rc,anxi,anyi,anzi &
       ,nx,ny,nz,myid_md)
  call bacopy(rc,myid_md,mpi_md_world,iprint,ifcoulomb &
       ,.true.,boundary)
!-----Make pair list
  tmp = mpi_wtime()
  call mk_lspr_para(namax,natm,nbmax,nb,nnmax,tag,ra,rc+rbuf,rc1nn &
       ,h,hi,anxi,anyi,anzi,lspr,ls1nn,iprint,.true.)
  tlspr = tlspr +(mpi_wtime() -tmp)

!.....Calc forces
  lstrs = lstrs0 .or. &
       ( trim(cpctl).eq.'Berendsen' .or. &
       trim(cpctl).eq.'vc-Berendsen' .or. &
       trim(cpctl).eq.'vv-Berendsen')
!.....Cell is new at the first call of get_force
  lcell_updated = .true.
  call get_force(namax,natm,tag,ra,nnmax,aa,strs,chg,chi,stnsr &
       ,h,hi,tcom,nb,nbmax,lsb,lsex,nex,lsrc,myparity,nn,sv,rc &
       ,lspr,sorg,mpi_md_world,myid_md,epi,epot0,nspmax,specorder,lstrs &
       ,ifcoulomb,iprint,.true.,lvc,lcell_updated,boundary)
  lcell_updated = .false.
  lstrs = .false.
  epot= epot0
  epotp = 0d0
!      write(6,'(a,f12.4,200f10.4)') 'erg,chg after get_force = '
!     &     ,epot,chg(1:natm)

!.....Structure analysis
  if( trim(cstruct).eq.'CNA' ) then
    if( rc1nn.gt.(rc+rbuf)/2 .and. myid_md.eq.0 ) then
      print *,' WARNING: rc1nn > (rc+rbuf)/2, ' &
           //'which may cause problem in CNA.'
      print *,'   rc,rbuf,rc1nn =',rc,rbuf,rc1nn
    endif
    call cna(namax,natm,nb,nnmax,ls1nn,tag)
  else if( trim(cstruct).eq.'a-CNA' ) then
    call acna(namax,natm,nb,nnmax,lspr,h,ra,tag)
  endif
!.....Constraints
  if( lconst ) then
    call update_const(namax,natm,tag,ra,h,0,maxstp)
  endif
!.....Force_modify
  if( lrdcfrc ) then
    call reduce_forces(namax,natm,aa,tag,ra,h,nnmax,lspr)
  endif

#ifdef __DISL__
  call perf_disl_pos_by_pot(epith,natm,ra,h,epi,sorg &
       ,nodes_md,myid_md,mpi_md_world,0,21)
#endif

!-----calc kinetic energy
  call get_ekin(namax,natm,va,tag,h,nspmax,fekin,ekin,eki,ekl &
       ,vmax,mpi_md_world)
  vmaxold=vmax

  if( trim(cpctl).eq.'Berendsen' .or. &
       trim(cpctl).eq.'vc-Berendsen' ) then
    if(myid_md.eq.0 .and. iprint.ne.0 ) then
      write(6,*) ''
      write(6,'(a)') ' Barostat: variable-cell Berendsen'
      write(6,'(a,6f10.3)') '   Target stress [GPa]: ' &
           ,stgt(1,1),stgt(2,2),stgt(3,3) &
           ,stgt(2,3),stgt(3,1),stgt(1,2)
    endif
    stgt(1:3,1:3)= stgt(1:3,1:3) *gpa2up
  else if( trim(cpctl).eq.'vv-Berendsen' ) then
    if( abs(pini-pfin).gt. 0.1d0 ) then
      if(myid_md.eq.0 .and. iprint.ne.0 ) then
        write(6,*) ''
        write(6,'(a)') ' Barostat: variable-volume Berendsen'
        write(6,'(a,f0.3,a,f0.3,a)') &
             '   Target pressure from = ',pini,' GPa to ' &
             ,pfin,' GPa'
      endif
      ptgt = pini *gpa2up
    else
      if(myid_md.eq.0 .and. iprint.ne.0 ) then
        write(6,*) ''
        write(6,'(a)') ' Barostat: variable-volume Berendsen'
        write(6,'(a,f0.3,a)') &
             '   Target pressure = ',ptgt,' GPa'
      endif
      ptgt = ptgt *gpa2up
    endif
  endif

  call force_isobaric(stgt,ptgt,ah,natm,eki,strs,sgm &
       ,dt,srlx,stbeta,vol,stnsr,mpi_md_world,cpctl)
  prss = (stnsr(1,1)+stnsr(2,2)+stnsr(3,3))/3*up2gpa

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
    do ifmv=1,9
      if( ndof(ifmv).eq.0 ) cycle
      temp(ifmv)= ekl(ifmv) /ndof(ifmv) /fkb *2d0
!          print *,' ifmv,ekl,temp = ',ifmv,ekl(ifmv),temp(ifmv)
      nave= nave +ndof(ifmv)
      write(6,'(1x,a,i1,a,f16.5,a)') "  Temperature ",ifmv &
           ,"   = ",temp(ifmv),' K'
      tave= tave +temp(ifmv)*ndof(ifmv)
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
    write(6,'(1x,a,6(1x,f0.3))') '  Stress tensor   =' &
         ,stnsr(1,1)*up2gpa,stnsr(2,2)*up2gpa,stnsr(3,3)*up2gpa &
         ,stnsr(2,3)*up2gpa,stnsr(3,1)*up2gpa,stnsr(1,2)*up2gpa
    write(6,*) ''

!!$    print '(a,20f8.5)',' alphas=',ol_alphas(0,1:natm)
    if( tave.gt.10000d0 ) cftave = 'es12.4'
    tcpu = mpi_wtime() -tcpu0
    write(6,'(a,'//cfistp//','//cfetime//','//cftave &
         //',es13.4,2es11.3)') &
         " istp,etime,temp,epot,vol,prss=" &
         ,istp,tcpu,tave,epot,vol,prss
  endif

!.....output initial configuration including epi, eki, and strs
!      write(cnum(1:4),'(i4.4)') 0
  write(cnum,'(i0)') 0
  tmp = mpi_wtime()
  call space_comp(ntot0,tagtot,rtot,vtot,atot,epitot,ekitot &
       ,stot,chgtot,chitot,natm,tag,ra,va,aa,epi,eki,strs &
       ,chg,chi,sorg,nxyz,myid_md,mpi_md_world,tspdcmp)
  if( myid_md.eq.0 ) then
    if( ifsort.eq.1 ) call sort_by_tag(ntot0,tagtot,rtot,vtot &
         ,atot,ekitot,epitot,stot,chgtot,chitot)
    if( ifpmd.eq.1 ) then  ! pmd format
      if( trim(ciofmt).eq.'bin' .or. trim(ciofmt).eq.'binary' ) &
           then
        call write_pmdtot_bin(20,"pmd_"//trim(cnum))
      elseif( trim(ciofmt).eq.'ascii' ) then
        call write_pmdtot_ascii(20,"pmd_"//trim(cnum))
      endif
    else if( ifpmd.eq.2 ) then ! LAMMPS-dump format
      call write_dump(20,'dump_'//trim(cnum))
    endif
  endif
  tdump = tdump +(mpi_wtime() -tmp)

  if( trim(ctctl).eq.'ttm' ) then
    call assign_atom2cell(namax,natm,ra,sorg,boundary)
    call calc_Ta(namax,natm,nspmax,h,tag,va,fmv,fekin &
         ,istp,myid_md,mpi_md_world,iprint)
    call output_Te(0,simtime,myid_md,iprint)
  endif

!-----initialize the counter for output
  iocntpmd=0
  iocnterg=0

  if(myid_md.eq.0) then
!.....write out energies
    open(ioerg,file="out.erg",status='replace')
    write(ioerg,'(a)') '# 1:istp, 2:simtime[fs],' &
         //'   3:etot[eV],  4:ekin,' &
         //'  5:epot,  6:temp[K],  7:vol[Ang^3],  8:pressure[GPa]'
    write(ioerg,'(a,es16.7e3,a)') '#  Epot0 =',epot0,' [eV]'
    if( tave.gt.10000d0) cftave = 'es12.4'
    write(ioerg,'('//cfistp//','//cfstime//',3es16.7e3' &
         //','//cftave//',2es16.7e3)') istp &
         ,simtime,ekin+epot0,ekin,epot0,tave &
         ,vol &
         ,(stnsr(1,1)+stnsr(2,2)+stnsr(3,3))/3*up2gpa
    call flush(ioerg)
!c.....write out temperatures
!        open(iotemp,file='out.temperature',status='replace')
!        write(iotemp,'(a)') '# istp, temperature[0-9]'
!        ediff(1:9)= 0d0
!        write(iotemp,'(i10,18es16.7e3)') istp,temp(1:9),ediff(1:9)
!        call flush(iotemp)

!.....write out tensile forces
    if(trim(czload_type).eq.'atoms' .or. trim(czload_type).eq.'box' &
         .or. trim(czload_type).eq.'shear' ) then
      open(iostrs,file='out.zload',status='replace')
      write(iostrs,'(a)') '#  istp,   strnow,      ftop,     fbot'
    endif
  endif

  if( trim(czload_type).eq.'atoms' .or. &
       trim(czload_type).eq.'shear' ) then
    call get_forces_on_base(natm,ra,aa,tag,h,ftop &
         ,fbot,sorg,myid_md,mpi_md_world,iprint,czload_type)
    if( myid_md.eq.0 ) then
      write(iostrs,'(i8,3es15.7)') 0,0.0,ftop,fbot
      call flush(iostrs)
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


  i_conv = 0
  lconverged = .false.
!-----velocity-Verlet loop starts---------------------------------------
  do istp=1,maxstp

!.....Metadynamics
    if( lmetaD ) then
      call update_metaD(istp,namax,natm,nsp,tag,ra,h,nnmax,lspr &
           ,myid_md,mpi_md_world,iprint)
    endif

!.....In case of isobaric MD, lstrs has to be always TRUE.
    if( trim(cpctl).eq.'Berendsen' .or. &
         trim(cpctl).eq.'vc-Berendsen' .or. &
         trim(cpctl).eq.'vv-Berendsen' .or. lstrs0 ) then
      lstrs = .true.
    else
      lstrs = .false.
    endif

!.....Variable time-step
    if( lvardt ) then
      call get_vmax(namax,natm,va,h,vmax,mpi_md_world)
      dt = min( dtmax, vardt_len/vmax )
      if( myid_md.eq.0 ) then
        if( iprint.gt.10 ) then
          print *,'dt,vmax,vardt_len/vmax = ',dt,vmax,vardt_len/vmax
        else if( iprint.gt.1 ) then
          print '(a,f8.5,a)',' update dt to ',dt,' fs'
        endif
      endif
!.....Update dt-related values
      if( trim(ctctl).eq.'Langevin' ) then
        tgmm = 1d0/trlx
        do ifmv=1,9
          if( ttgt(ifmv).lt.0d0 ) then
            tfac(ifmv)= -1d0
          else
            tfac(ifmv)= dsqrt(2d0*tgmm*ttgt(ifmv)/dt *k2ue)
          endif
        enddo
      else if( trim(ctctl).eq.'ttm' ) then
        call set_inner_dt(dt)
      endif
    endif

!-------first kick of velocities
    do i=1,natm
      is = int(tag(i))
      va(1:3,i)=va(1:3,i) +aa(1:3,i)*fa2v(is)*dt
    enddo

    if( ifdmp.eq.2 ) then
      call vfire(num_fire,alp0_fire,alp_fire,falp_fire,dtmax_fire &
           ,finc_fire,fdec_fire,nmin_fire &
           ,natm,va,aa,myid_md,mpi_md_world,dt,iprint)
    endif

!-------multiply fmv or damping
    do i=1,natm
      l= int(mod(tag(i)*10,10d0))
      va(1:3,i)=va(1:3,i) *fmv(1:3,l)
    enddo
!        call get_vmax(namax,natm,va,h,vmaxt,mpi_md_world)
!        print '(a,i5,2es12.4)','myid,vmax,vmaxt='
!     &       ,myid_md,vmax,vmaxt

!.....Update positions
    if( lconst ) then
      call update_const_pos(namax,natm,h,hi,tag,ra,va,dt,nspmax,am)
    else
      do i=1,natm
        ra(1:3,i)=ra(1:3,i) +va(1:3,i)*dt
      enddo
    endif
    ltot_updated = .false.

    if( trim(czload_type).eq.'atoms' ) then
      call zload_atoms(natm,ra,tag,maxstp,strfin,strnow &
           ,sorg,myid_md,mpi_md_world)
    else if( trim(czload_type).eq.'box' ) then
      call zload_box(natm,maxstp,istp,dt,strfin,strnow,h,myid_md)
    else if( trim(czload_type).eq.'shear' ) then
      call shear_atoms(natm,ra,tag,maxstp,strfin,strnow &
           ,sorg,myid_md,mpi_md_world)
    endif

    if( trim(cpctl).eq.'Berendsen' .or. &
         trim(cpctl).eq.'vc-Berendsen' .or. &
         trim(cpctl).eq.'vv-Berendsen' ) then
!     h(1:3,1:3,0)= matmul(ah,h(1:3,1:3,0))
      htmp(1:3,1:3) = matmul(ah,h(1:3,1:3,0))
      do i=1,3
        do j=1,3
!.....NOTE: The definition of lcellfix from input and
!.....that of h-matrix actually used are in relation of transpose.
          if( .not. lcellfix(j,i) ) then
            h(i,j,0) = htmp(i,j)
          else
            h(i,j,0) = h(i,j,0)
          endif
        enddo
      enddo
      lcell_updated = .true.
    endif

    if( trim(cdeform).ne.'none' ) then
      call apply_deform(h)
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
      if( trim(ctctl).eq.'ttm' ) then
        call remove_ablated_atoms(simtime,namax,natm &
             ,tag,ra,va,chg,chi,h,sorg)
      endif
!.....Move atoms that cross the boundary
      call bamove(tcom,namax,nbmax,natm,ra,va,tag,chg,chi &
           ,anxi,anyi,anzi,myid_md,nn,sv,myparity,mpi_md_world &
           ,boundary)
!.....Copy RA of boundary atoms
      call bacopy(rc,myid_md,mpi_md_world,iprint,ifcoulomb &
           ,.false.,boundary)
!.....Make pair list
      tmp = mpi_wtime()
      call mk_lspr_para(namax,natm,nbmax,nb,nnmax,tag,ra,rc+rbuf &
           ,rc1nn,h,hi,anxi,anyi,anzi,lspr,ls1nn,iprint,.false.)
      tlspr = tlspr +(mpi_wtime() -tmp)
      rbufres = rbuf
    else
!.....Copy RA of boundary atoms determined by 'bacopy'
      call bacopy_fixed(tcom,sgm,vol,lsb,lsex,nbmax,ra,namax &
           ,natm,nb,anxi,anyi,anzi,nn,tag,rc,myid_md,myparity,lsrc &
           ,sv,nex,mpi_md_world,ifcoulomb,chg,chi,boundary)
    endif

    if(ifpmd.gt.0.and. mod(istp,noutpmd).eq.0 )then
      lstrs = lstrs0
    endif
!-------Calc forces
    call get_force(namax,natm,tag,ra,nnmax,aa,strs,chg,chi,stnsr &
         ,h,hi,tcom,nb,nbmax,lsb,lsex,nex,lsrc,myparity,nn,sv,rc &
         ,lspr,sorg,mpi_md_world,myid_md,epi,epot,nspmax,specorder,lstrs &
         ,ifcoulomb,iprint,.false.,lvc,lcell_updated,boundary)
    lcell_updated = .false.
    lstrs = .false.
!.....Structure analysis
    if( trim(cstruct).eq.'CNA' &
         .and. mod(istp,istruct).eq.0 ) then
      call cna(namax,natm,nb,nnmax,ls1nn,tag)
    else if( trim(cstruct).eq.'a-CNA' &
         .and. mod(istp,istruct).eq.0 ) then
      call acna(namax,natm,nb,nnmax,lspr,h,ra,tag)
    endif
!.....Force from metadynamics
    if( lmetaD ) then
      call force_metaD(istp,namax,natm,tag,ra,aa,h,hi,epot &
           ,nnmax,lspr,myid_md,mpi_md_world,iprint)
    endif
!.....Update some variables in constraints if needed
    if( lconst ) then
      call update_const(namax,natm,tag,ra,h,istp,maxstp)
    endif
!.....Force_modify
    if( lrdcfrc ) call reduce_forces(namax,natm,aa,tag,ra &
         ,h,nnmax,lspr)

!.....Second kick of velocities
    do i=1,natm
      is = int(tag(i))
      va(1:3,i)=va(1:3,i) +aa(1:3,i)*fa2v(is)*dt
    enddo

!.....Constraints for velocities
    if( lconst ) then
      call update_const_vel(namax,natm,h,hi,tag,va,dt,nspmax,am)
    endif

!.....Calc kinetic energy
    vmaxold= vmax
    call get_ekin(namax,natm,va,tag,h,nspmax,fekin,ekin,eki,ekl &
         ,vmax,mpi_md_world)

    if( trim(ctctl).eq.'Langevin' ) then
!.....Langevin thermostat with Mannella integrator
      hscl(1:3)= 0d0
      do l=1,3
        hscl(l)= dsqrt(h(1,l,0)**2 +h(2,l,0)**2 +h(3,l,0)**2)
      enddo
!.....If final temperature is assigned,
!.....target temperatures are forced to be set intermediate temperatures
      if( tfin.gt.0d0 ) then
        ttgt(1:9) = tinit +(tfin-tinit)*istp/maxstp
        tfac(1:9) = dsqrt(2d0*tgmm*ttgt(1:9)/dt *k2ue)
      endif
      do i=1,natm
!            ifmv= int(mod(tag(i)*10,10d0))
        ifmv = ifmvOf(tag(i))
        is = int(tag(i))
        if( .not. (ifmv.eq.0 .or. tfac(ifmv).lt.0d0) ) then
          ami= am(is)
          aai(1:3)= 0d0
!.....Here unit of TMP should be [eV/Ang],
!     whereas TFAC is [eu/Ang/sqrt(ump)], so need to multiply ue2ev
          tmp= tfac(ifmv)*dsqrt(ami)*ue2ev
!.....Here the unit of va*tgmm*ami is [ump*(Ang)/fs**2 = ue/(Ang)],
!.....where (Ang) is actually scaled to unitless,
!.....but this should be [eV/Ang], so need to multiply ue2ev.
          do l=1,3
            aai(l)= -va(l,i)*tgmm*ami*ue2ev &
                 +tmp*box_muller()/hscl(l)
          enddo
          if( iprint.gt.1 ) then
!.....accumulate energy difference
            vi(1:3)= h(1:3,1,0)*va(1,i) &
                 +h(1:3,2,0)*va(2,i) &
                 +h(1:3,3,0)*va(3,i)
            vl(1:3)= h(1:3,1,0)*aai(1)*fa2v(is)*dt *2d0 &
                 +h(1:3,2,0)*aai(2)*fa2v(is)*dt *2d0 &
                 +h(1:3,3,0)*aai(3)*fa2v(is)*dt *2d0
            ediff(ifmv)= ediff(ifmv) +fekin(is) &
                 *(2d0*sprod(3,vi,vl)+sprod(3,vl,vl))
          endif
!.....To compensate the factor 1/2 in fa2v, multiply 2 here.
          va(1:3,i)= va(1:3,i) +aai(1:3)*fa2v(is)*dt *2d0
        endif
      enddo
!-------temperature control by Berendsen thermostat
    else if( trim(ctctl).eq.'Berendsen' ) then
      tfac(1:9)= 0d0
!.....if final temperature is assigned,
!.....target temperatures are forced to be set intermediate temperatures
      if( tfin.gt.0d0 ) then
        ttgt(1:9) = tinit +(tfin-tinit)*istp/maxstp
      endif
      do ifmv=1,9
        if(ndof(ifmv).le.0 .or. ttgt(ifmv).lt.0d0 ) cycle
        temp(ifmv)= ekl(ifmv) *2d0 /fkb /ndof(ifmv)
        if( abs(ttgt(ifmv)-temp(ifmv))/temp(ifmv).gt.100d0 ) then
          tfac(ifmv)= dsqrt(1d0 +dt/trlx*100d0 )
        else
          tfac(ifmv)= dsqrt(1d0 +dt/trlx*(ttgt(ifmv)-temp(ifmv)) &
               /temp(ifmv))
        endif
      enddo
      do i=1,natm
!            ifmv= int(mod(tag(i)*10,10d0))
        ifmv = ifmvOf(tag(i))
        if( ifmv.le.0 .or. ttgt(ifmv).lt.0d0 ) cycle
        va(1:3,i)= va(1:3,i) *tfac(ifmv)
!.....accumulate energy difference
        is= int(tag(i))
        ediff(ifmv)= ediff(ifmv) +(tfac(ifmv)**2-1d0)*fekin(is)
      enddo
    else if( trim(ctctl).eq.'ttm' ) then
      call assign_atom2cell(namax,natm,ra,sorg,boundary)
      call calc_Ta(namax,natm,nspmax,h,tag,va,fmv,fekin &
           ,istp,myid_md,mpi_md_world,iprint)
      call langevin_ttm(namax,natm,va,aa,tag,am,h &
           ,nspmax,fa2v,fekin,ediff,dt,myid_md,mpi_md_world,iprint)
      call update_Te(simtime,dt,myid_md,mpi_md_world,iprint)
    endif  ! end of thermostat


    if( abs(pini-pfin).gt.0.1d0 ) then
      ptgt = ( pini +(pfin-pini)*istp/maxstp ) *gpa2up
    endif
    call force_isobaric(stgt,ptgt,ah,natm,eki,strs,sgm &
         ,dt,srlx,stbeta,vol,stnsr,mpi_md_world,cpctl)
    prss = (stnsr(1,1)+stnsr(2,2)+stnsr(3,3))/3*up2gpa

    if( nrmtrans.gt.0 .and. mod(istp,nrmtrans).eq.0 ) then
      call rm_trans_motion(natm,tag,va,nspmax,am &
           ,mpi_md_world,myid_md,iprint)
    endif


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
      do ifmv=1,9
        if( ndof(ifmv).le.0 ) cycle
        temp(ifmv)= ekl(ifmv) *2d0 /fkb /ndof(ifmv)
        nave= nave +ndof(ifmv)
        tave= tave +temp(ifmv)*ndof(ifmv)
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

      if( myid_md.eq.0 .and. iprint.gt.0 ) then
!.....write energies
        if( tave.gt.10000d0) cftave = 'es12.4'
        write(ioerg,'('//cfistp//','//cfstime//',3es16.7e3' &
             //','//cftave//',2es16.7e3)') istp &
             ,simtime,ekin+epot,ekin,epot,tave &
             ,vol &
             ,(stnsr(1,1)+stnsr(2,2)+stnsr(3,3))/3*up2gpa
        call flush(ioerg)
!.....write temperature
!            write(iotemp,'(i10,18es16.7e3)') istp,temp(1:9),ediff0(1:9)
!            call flush(iotemp)
        ediff(1:9)= 0d0

        if( trim(czload_type).eq.'atoms' .or. &
             trim(czload_type).eq.'shear' ) then
          write(iostrs,'(i8,3es15.7)') istp,strnow,ftop,fbot
          call flush(iostrs)
        endif
      endif   ! myid_md.eq.0

!---------output step, time, and temperature
      tcpu= mpi_wtime() -tcpu1
      if( myid_md.eq.0 .and. iprint.gt.0 ) then
!!$        print '(a,20f8.5)',' alphas=',ol_alphas(0,1:natm)
        if( tave.gt.10000d0 ) cftave = 'es12.4'
        tcpu = mpi_wtime() -tcpu0
        write(6,'(a,'//cfistp//','//cfetime//','//cftave &
             //',es13.4,2es11.3)') &
             " istp,etime,temp,epot,vol,prss=" &
             ,istp,tcpu,tave,epot,vol,prss
        if( trim(cpctl).eq.'Berendsen' .or. &
             trim(cpctl).eq.'vc-Berendsen' .or. &
             trim(cpctl).eq.'vv-Berendsen' .and. &
             iprint.ne.0 ) then
!!$          write(6,'(a)') ' Cell-matrix:' !,h(1:3,1:3,0)
!!$          write(6,'(3f15.7)') h(1,1:3,0)
!!$          write(6,'(3f15.7)') h(2,1:3,0)
!!$          write(6,'(3f15.7)') h(3,1:3,0)
          write(6,'(a)') ' Lattice vectors:' !,h(1:3,1:3,0)
          write(6,'(a,"[ ",3f12.3," ]")') '   a = ',h(1:3,1,0)
          write(6,'(a,"[ ",3f12.3," ]")') '   b = ',h(1:3,2,0)
          write(6,'(a,"[ ",3f12.3," ]")') '   c = ',h(1:3,3,0)

          write(6,'(a,6f10.3)') ' Stress (GPa):' &
               ,stnsr(1,1)*up2gpa ,stnsr(2,2)*up2gpa &
               ,stnsr(3,3)*up2gpa ,stnsr(2,3)*up2gpa &
               ,stnsr(1,3)*up2gpa ,stnsr(1,2)*up2gpa
        endif
        call flush(6)
      endif

      if( trim(ctctl).eq.'ttm' ) then
        call output_Te(istp,simtime,myid_md,iprint)
        call output_energy_balance(istp,simtime,myid_md,iprint)
      endif

#ifdef __DISL__
!.....Output disl core pos
      call perf_disl_pos_by_pot(epith,natm,ra,h,epi,sorg &
           ,nodes_md,myid_md,mpi_md_world,iocnterg,21)
#endif
    endif ! energy output

!.....check convergence criteria if it is dumped MD
    if( ifdmp.gt.0 .and. epot-epotp.le.0d0 .and. n_conv.gt.0 .and. &
         abs(epot-epotp).lt.eps_conv .and. istp.gt.minstp ) then
!          print *,'ifdmp.gt.0 = ',ifdmp.gt.0
!          print *,'epot-epotp.le.0d0 = ',epot-epotp.le.0d0
!          print *,'abs(epot-epotp).lt.eps_conv = ',abs(epot-epotp)
!     &         .lt.eps_conv
!          print *,'istp.gt.minstp = ',istp.gt.minstp
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

!-------write the particle positions
    if(ifpmd.gt.0.and. &
         (mod(istp,noutpmd).eq.0.or.lconverged) )then
      tmp = mpi_wtime()
!---------decide pmd-file name
      iocntpmd=iocntpmd+1
!          write(cnum(1:4),'(i4.4)') iocntpmd
      write(cnum,'(i0)') istp
!          call system("mkdir -p "//cnum)
      call space_comp(ntot0,tagtot,rtot,vtot,atot,epitot,ekitot &
           ,stot,chgtot,chitot,natm,tag,ra,va,aa,epi,eki,strs &
           ,chg,chi,sorg,nxyz,myid_md,mpi_md_world,tspdcmp)
      ltot_updated = .true.
      if( myid_md.eq.0 ) then
        if( ifsort.eq.1 ) call sort_by_tag(ntot0,tagtot,rtot,vtot &
             ,atot,ekitot,epitot,stot,chgtot,chitot)
        if( ifpmd.eq.1 ) then  ! pmd format
          if( trim(ciofmt).eq.'bin' .or. trim(ciofmt).eq.'binary' ) &
               then
            call write_pmdtot_bin(20,"pmd_"//trim(cnum))
          elseif( trim(ciofmt).eq.'ascii' ) then
            call write_pmdtot_ascii(20,"pmd_"//trim(cnum))
          endif
        else if( ifpmd.eq.2 ) then  ! LAMMPS-dump format
          call write_dump(20,'dump_'//trim(cnum))
        endif
      endif
      tdump = tdump +(mpi_wtime() -tmp)
    endif

    if( lconverged ) exit
  enddo ! end of istp

  if( .not. ltot_updated ) then
    call space_comp(ntot0,tagtot,rtot,vtot,atot,epitot,ekitot &
         ,stot,chgtot,chitot,natm,tag,ra,va,aa,epi,eki,strs &
         ,chg,chi,sorg,nxyz,myid_md,mpi_md_world,tspdcmp)
  endif

  nstp_done = istp
  tcpu2= mpi_wtime()

  if(myid_md.eq.0) then
    close(ioerg)
!        close(iotemp)
    if(czload_type.eq.'atoms' .or. czload_type.eq.'box') then
      close(iostrs)
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
    do ifmv=1,9
      if( ndof(ifmv).le.0 ) cycle
      temp(ifmv)= ekl(ifmv) *2d0 /fkb /ndof(ifmv)
      nave= nave +ndof(ifmv)
      write(6,'(1x,a,i1,a,f16.5,a)') "  Temperature ",ifmv &
           ,"   = ",temp(ifmv),' K'
      tave= tave +temp(ifmv)*ndof(ifmv)
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
    if( trim(cpctl).eq.'Berendsen' .or. &
         trim(cpctl).eq.'vc-Berendsen' .or. &
         trim(cpctl).eq.'vv-Berendsen' ) then
      call cell_info(h)
      write(6,'(a,f16.5,a)')  '  Cell volume     = ',vol,' Ang^3'
    endif
    write(6,*) ''
    if( iprint.gt.1 ) then
      write(6,'(1x,a,f10.2)') "Time for space decomp = ",tspdcmp
      write(6,'(1x,a,f10.2)') "Time for comm         = ",tcom
      write(6,'(1x,a,f10.2)') "Time for neighbor     = ",tlspr
      write(6,'(1x,a,f10.2)') "Time for dump         = ",tdump
    endif
    if( trim(ctctl).eq.'ttm' ) then
      write(6,'(1x,a,f10.2)') "Time for TTM          = ",t_ttm
    endif

    ihour = int(tcpu/3600)
    imin  = int((tcpu-ihour*3600)/60)
    isec  = int(tcpu -ihour*3600 -imin*60)
    write(6,'(1x,a,f10.2,a,i3,"h",i2.2,"m",i2.2,"s")') &
         "Time                  = ",tcpu, &
         " sec  = ",ihour,imin,isec
  endif

!.....revert forces to the unit eV/A before going out pmd_core
  if( myid_md.eq.0 ) then
    do i=1,ntot0
      is= int(tagtot(i))
      aai(1:3)= h(1:3,1,0)*atot(1,i) &
           +h(1:3,2,0)*atot(2,i) &
           +h(1:3,3,0)*atot(3,i)
      atot(1:3,i) = aai(1:3)
    enddo
  endif

!.....Output metadynamics potential
  if( lmetaD ) then
    call write_metaD_potential(maxstp,nsp,myid_md,iprint)
  endif

  if( lrdcfrc ) then
    call finalize_rdcfrc()
  endif

!.....deallocate all the arrays allocated in pmd_core
  if( ltdst ) then
    deallocate(tdst,nadst)
  endif
  deallocate(fekin,fa2v)
  deallocate(ra,va,aa,ra0,strs,stt,tag,lspr,ls1nn &
       ,epi,eki,stp,stn,lsb,lsex,chg,chi)

end subroutine pmd_core
!=======================================================================
subroutine one_shot(hunit,h,ntot0,tagtot,rtot,vtot,atot,stot &
     ,ekitot,epitot,chgtot,chitot &
     ,myid_md,mpi_md_world,nodes_md,nx,ny,nz &
     ,specorder,am,dt,rc,rbuf,rc1nn,stnsr,epot &
     ,ekin,ifcoulomb,lvc,iprint,lcalcgrad,ndimp,maxisp &
     ,gwe,gwf,gws,lematch,lfmatch,lsmatch,boundary)
!
!  In case that only one shot force calculation is required,
!  especially called from fitpot.
!
  use pmdio,only: namax,nbmax,nnmax,nspmax
  use pmdvars
  use force
  use Morse,only: gradw_Morse,gradw_vcMorse
  use Coulomb,only: gradw_Coulomb
  use linreg,only: gradw_linreg
  use NN2,only: gradw_NN2
  implicit none
  include "mpif.h"
  include "./params_unit.h"
  integer,intent(in):: ntot0,myid_md,mpi_md_world &
       ,iprint,nodes_md,nx,ny,nz
  integer,intent(inout):: ifcoulomb
  real(8),intent(in):: hunit
  real(8),intent(inout):: tagtot(ntot0),rtot(3,ntot0),vtot(3,ntot0) &
       ,atot(3,ntot0),stot(3,3,ntot0),ekitot(3,3,ntot0) &
       ,epitot(ntot0),rc,rbuf,rc1nn,h(3,3,0:1),stnsr(3,3) &
       ,epot,ekin,am(nspmax),dt &
       ,chgtot(ntot0),chitot(ntot0)
!      character(len=20),intent(in):: cffs(numff)
  logical,intent(in):: lcalcgrad
  integer,intent(in):: ndimp,maxisp
  real(8),intent(inout):: gwe(ndimp),gwf(ndimp,3,ntot0),gws(ndimp,6)
  logical,intent(inout):: lvc
  logical,intent(in):: lematch,lfmatch,lsmatch
  character(len=3),intent(in):: boundary,specorder(nspmax)

  integer:: i,ierr,is,nspl,nxyz,iprm0
  real(8):: aai(3),epott
  logical:: lstrs = .false.
  logical:: lcell_updated = .false.
  character(len=3):: csp

!      print *,'one_shot: 01'
  call initialize_pmdvars(nspmax)
  call set_use_charge()

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
!.....perform space decomposition after reading atomic configuration
!      write(6,'(a,30es10.2)') 'chgtot before space_decomp = ',
!     &     chgtot(1:ntot0)
  if( iprint.gt.0 ) then
    print *,''
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
!      print *,'one_shot: 02'
  call space_decomp(hunit,h,ntot0,tagtot,rtot,vtot,chgtot,chitot &
       ,myid_md,mpi_md_world,nx,ny,nz,nxyz,rc,rbuf,iprint)
!      write(6,'(a,200es10.2)') 'chg after space_decomp = ',chg(1:natm)

!      print *,'one_shot: 03'
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
  allocate(fekin(nspmax),fa2v(nspmax))
  call setup(nspmax,am,fekin,fa2v)
!-----set HI and SGM
  call boxmat(h,hi,ht,g,gi,gt,vol,sgm)
!-----ntset
  call ntset(myx,myy,myz,nx,ny,nz,nn,sv,myparity,anxi,anyi,anzi)

  tcpu1= mpi_wtime()
  tcom = 0d0

!      print *,'one_shot: 04'
!      if( iprint.ge.10 ) print *,'init_force,myid_md=',myid_md
  call init_force(namax,natm,nsp,tag,chg,chi,myid_md,mpi_md_world, &
       iprint,h,rc,lvc,ifcoulomb,specorder)

!      print *,'one_shot: 05'
!-----copy RA of boundary atoms
  call check_size_and_parallel(sgm,vol,rc,anxi,anyi,anzi &
       ,nx,ny,nz,myid_md)
  call bacopy(rc,myid_md,mpi_md_world,iprint,ifcoulomb &
       ,.true.,boundary)
!-----Make pair list
  call mk_lspr_para(namax,natm,nbmax,nb,nnmax,tag,ra,rc+rbuf &
       ,rc1nn,h,hi,anxi,anyi,anzi,lspr,ls1nn,iprint,.true.)
  lstrs = .true.

!      print *,'one_shot: 06'
!      if( iprint.ge.10 ) print *,'get_force,myid_md,lcalcgrad='
!     &     ,myid_md,lcalcgrad
  if( .not.lcalcgrad ) then
    call get_force(namax,natm,tag,ra,nnmax,aa,strs,chg,chi,stnsr &
         ,h,hi,tcom,nb,nbmax,lsb,lsex,nex,lsrc,myparity,nn,sv,rc &
         ,lspr,sorg,mpi_md_world,myid_md,epi,epot,nspmax,specorder,lstrs &
         ,ifcoulomb,iprint,.true.,lvc,lcell_updated,boundary)
    if( iprint.gt.0 ) then
      write(6,'(a,es15.7)') ' Potential energy = ',epot
    endif
  else  ! lcalcgrad = .true.
    if( iprint.gt.0 ) print *,'Computing gradient...'
    epot = 0d0
    gwe(1:ndimp) = 0d0
    gwf(1:ndimp,1:3,1:natm) = 0d0
    gws(1:ndimp,1:6) = 0d0
    if( use_force('Morse') &
         .and. use_force('screened_Coulomb') ) then
      iprm0 = 0
      call gradw_Coulomb(namax,natm,nb,tag,ra,chg, &
           nnmax,h,rc,lspr,epott,iprint,ndimp,gwe,gwf,gws, &
           lematch,lfmatch,lsmatch,iprm0,myid_md,mpi_md_world)
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
    else if( use_force('NN2') ) then
      iprm0 = 0
      call gradw_NN2(namax,natm,tag,ra,nnmax,h,rc,lspr, &
           iprint,ndimp,gwe,gwf,gws,lematch,lfmatch,lsmatch,iprm0)
    endif
!        if( use_force('vcMorse') ) call gradw_vcMorse(namax,natm,tag,ra
!     &       ,nnmax,chg,h,rc,lspr,epott,iprint,ndimp,gwe,gwf,gws)
!.....Derivative of stress should be divided by the cell volume
    gws(1:ndimp,1:6) = gws(1:ndimp,1:6) /vol
  endif

!      print *,'one_shot: 07'
  if( iprint.gt.0 ) print *,'Compute stresses...'
  call sa2stnsr(natm,strs,eki,stnsr,vol,mpi_md_world)

  call space_comp(ntot0,tagtot,rtot,vtot,atot,epitot,ekitot &
       ,stot,chgtot,chitot,natm,tag,ra,va,aa,epi,eki,strs &
       ,chg,chi,sorg,nxyz,myid_md,mpi_md_world,tspdcmp)
  if( iprint.gt.0 ) print *,'Compute stresses done'

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

  if( iprint.ge.20 .and. myid_md.eq.0 .and. lvc ) then
    print *,'Atomic charges:'
    do i=1,ntot0
      print '(a,i5,es11.3,4f7.3)','   ia,epi,chg,frc=',i &
           ,epitot(i),chgtot(i),atot(1:3,i)
    enddo
  endif

!.....deallocate all the arrays allocated 
  deallocate(fekin,fa2v)
!      deallocate(ra,va,aa,ra0,strs,stt,tag,lspr
!     &     ,epi,eki,stp,stn,lsb,lsex,chg,chi)

  if( iprint.ge.10 ) then
    write(6,'(a)') ' one_shot done.'
    write(6,*) ''
  endif
  return
end subroutine one_shot
!=======================================================================
subroutine setup(nspmax,am,fekin,fa2v)
  implicit none
  include "params_unit.h"
  integer,intent(in):: nspmax
  real(8),intent(in):: am(nspmax)
  real(8),intent(out):: fekin(nspmax),fa2v(nspmax)

  integer:: i

!.....Factors for the following quantities
!.....Force: eV/Ang
!.....Acceleration: Ang/fs/fs
!.....Velocity: Ang/fs
!.....Position: Ang
!.....Energy: eV = 1.0/1.602e-19 [J = kg*m**2/s**2]
  do i=1,nspmax
    fa2v(i) = 0.5d0 *ev2j *(kg2ump *m2ang**2/s2fs**2) /am(i)
    fekin(i)= 0.5d0 *(am(i)*ump2kg) *ang2m**2 /fs2s**2 *j2ev
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
!      write(6,'(a)') ' g:'
!      write(6,'(3es12.4)') g(1,1:3,0)
!      write(6,'(3es12.4)') g(2,1:3,0)
!      write(6,'(3es12.4)') g(3,1:3,0)
!      write(6,'(a)') ' gi:'
!      write(6,'(3es12.4)') gi(1,1:3)
!      write(6,'(3es12.4)') gi(2,1:3)
!      write(6,'(3es12.4)') gi(3,1:3)

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
subroutine get_ekin(namax,natm,va,tag,h,nspmax,fekin,ekin,eki,ekl &
     ,vmax,mpi_md_world)
  use util,only: ifmvOf
  implicit none 
  include "mpif.h"
  integer,intent(in):: namax,natm,mpi_md_world,nspmax
  real(8),intent(in):: va(3,namax),h(3,3),fekin(nspmax) &
       ,tag(namax)
  real(8),intent(out):: ekin,eki(3,3,namax),vmax,ekl(9)
!-----locals
  integer:: i,ierr,is,ixyz,jxyz,imax,ifmv
  real(8):: ekinl,x,y,z,v(3),v2,vmaxl,ekll(9)
!!$  integer,external:: ifmvOf

  ekinl=0d0
  eki(1:3,1:3,1:natm)= 0d0
  ekll(1:9)= 0d0
  vmaxl= 0d0

  do i=1,natm
    is= int(tag(i))
!.....ifmv= int(mod(tag(i)*10,10d0))
    ifmv = ifmvOf(tag(i))
    if( ifmv.eq.0 ) cycle
    x= va(1,i)
    y= va(2,i)
    z= va(3,i)
    v(1:3)= h(1:3,1)*x +h(1:3,2)*y +h(1:3,3)*z
!        print *,'i,is,va,v=',i,is,x,y,z,v(1:3)
!.....Tensor form eki
    do jxyz=1,3
      do ixyz=1,3
        eki(ixyz,jxyz,i)= v(ixyz)*v(jxyz)
      enddo
    enddo
    v2= eki(1,1,i) +eki(2,2,i) +eki(3,3,i)
    eki(1:3,1:3,i)=eki(1:3,1:3,i)*fekin(is)
    ekinl=ekinl +eki(1,1,i) +eki(2,2,i) +eki(3,3,i)
!.....ekin of each ifmv
    ekll(ifmv)= ekll(ifmv) +eki(1,1,i) +eki(2,2,i) +eki(3,3,i)
!.....Find max speed
    if( v2.gt.vmaxl ) imax=i
    vmaxl=max(vmaxl,v2)
  enddo
!      print *,'imax,vmax,ekin=',imax,vmaxl,ekinl
!      print *,'va(1:3,imax)=',va(1:3,imax)

  ekin= 0d0
  call mpi_allreduce(ekinl,ekin,1,mpi_real8 &
       ,MPI_SUM,mpi_md_world,ierr)
  ekl(1:9)= 0d0
  call mpi_allreduce(ekll,ekl,9,mpi_real8 &
       ,mpi_sum,mpi_md_world,ierr)
  vmax= 0d0
  call mpi_allreduce(vmaxl,vmax,1,mpi_real8 &
       ,mpi_max,mpi_md_world,ierr)
  vmax=dsqrt(vmax)

end subroutine get_ekin
!=======================================================================
subroutine get_vmax(namax,natm,va,h,vmax,mpi_md_world)
  implicit none 
  include "mpif.h"
  integer,intent(in):: namax,natm,mpi_md_world
  real(8),intent(in):: va(3,namax),h(3,3)
  real(8),intent(out):: vmax

  integer:: i,ierr
  real(8):: vx,vy,vz,v(3),v2,vmaxl

  vmaxl = 0d0
  do i=1,natm
    vx = va(1,i)
    vy = va(2,i)
    vz = va(3,i)
    v(1:3) = h(1:3,1)*vx +h(1:3,2)*vy +h(1:3,3)*vz
    v2 = v(1)*v(1) +v(2)*v(2) +v(3)*v(3)
    vmaxl = max(vmaxl,v2)
  enddo
  vmax= 0d0
  call mpi_allreduce(vmaxl,vmax,1,mpi_real8 &
       ,mpi_max,mpi_md_world,ierr)
  vmax=dsqrt(vmax)
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
subroutine get_num_dof(natm,tag,fmv,ndof,myid_md,mpi_md_world &
     ,iprint)
  implicit none
  include 'mpif.h'
  integer,intent(in):: natm,myid_md,mpi_md_world,iprint
!.....TODO: hard coding number 9 should be replaced...
  real(8),intent(in):: tag(natm),fmv(3,0:9)
  integer,intent(out):: ndof(9)

  integer:: i,l,k,ndofl(9),ierr
  real(8),parameter:: deps= 1d-12

  ndofl(1:9)= 0
  do i=1,natm
    l= int(mod(tag(i)*10,10d0))
    do k=1,3
      if( abs(fmv(k,l)).lt.0.5d0 ) cycle
      ndofl(l)= ndofl(l) +1
    enddo
  enddo

  ndof(1:9)= 0
  call mpi_allreduce(ndofl,ndof,9,mpi_integer,mpi_sum &
       ,mpi_md_world,ierr)
  if(myid_md.eq.0 .and. iprint.ne.0 ) &
       write(6,'(/,a,9(2x,i0))') &
       ' Degrees of freedom for each ifmv =',ndof(1:9)
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
subroutine bacopy(rc,myid,mpi_md_world,iprint,ifcoulomb,l1st &
     ,boundary)
!-----------------------------------------------------------------------
!  Exchanges boundary-atom data among neighbor nodes: tag and ra
!  Normal data passing and repeating the cell are mixed, and
!  automatically selected according to the size of the system.
!  
!  Note: parallelized to smaller than rcut should not happen.
!-----------------------------------------------------------------------
  use pmdio,only: namax,nbmax
  use pmdvars
  use force
  use pmdmpi,only: nx,ny,nz,nid2xyz
  implicit none
  include 'mpif.h'
  integer,intent(in):: myid,mpi_md_world,iprint
  real(8),intent(in):: rc
  integer,intent(in):: ifcoulomb
  logical,intent(in):: l1st
  character(len=3):: boundary

!.....integer:: status(MPI_STATUS_SIZE)
  integer:: i,j,kd,kdd,kul,kuh,ku,ierr,iex,ix,iy,iz,itmp,istatus
  integer:: nav,maxna,maxb,inode,nsd,nrc,nbnew
  real(8):: tcom1,tcom2,xi(3),rcv(3),asgm
  logical,external:: bbd
  real(8),save,allocatable:: dbuf(:,:),dbufr(:,:)
  logical:: lshort(3)

  if( l1st ) then
    if( allocated(dbuf) ) deallocate(dbuf,dbufr)
    if( luse_charge ) then
      allocate(dbuf(6,nbmax),dbufr(6,nbmax))
    else
      allocate(dbuf(4,nbmax),dbufr(4,nbmax))
    endif
  endif

  if( luse_charge ) then
    if( size(dbuf).lt.6*nbmax ) then
      deallocate(dbuf,dbufr)
      allocate(dbuf(6,nbmax),dbufr(6,nbmax))
    endif
  else
    if( size(dbuf).lt.4*nbmax ) then
      deallocate(dbuf,dbufr)
      allocate(dbuf(4,nbmax),dbufr(4,nbmax))
    endif
  endif

  call nid2xyz(myid,ix,iy,iz)

!-----reset the num of "received" boundary atoms
  nbnew=0

!-----calculate the cut-off lengths
  do kd=1,3
    asgm= dsqrt(sgm(1,kd)**2 +sgm(2,kd)**2 +sgm(3,kd)**2)
    rcv(kd)= rc*asgm/vol
    nex(kd)= int(rcv(kd)) +1
  enddo
  if( l1st .and. myid.eq.0 .and. iprint.ge.2 ) then
    write(6,'(a,3f10.3)') ' rcv = ',rcv(1:3)
    write(6,'(a,3i10)')   ' nex = ',nex(1:3)
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
    else                    ! long enough for normal boundary-atom copy
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
    endif                   ! if(nex.gt.1)

!.....Depending on boundary condition and cell position,
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
    call mpi_allreduce(itmp,maxb,1,mpi_integer,mpi_max, &
         mpi_md_world,ierr)
    if(maxb.gt.nbmax) then
      if (myid.eq.0 .and. iprint.ne.0 ) then
        print *,'Updated namax and array since nbmax changed' &
             //' from ',nbmax,' to ',int(maxb*1.5)
!            write(*,*)'Buffer or list overflowed at bacopy'
!            write(*,*)'LSB(0) NBMAX =',maxb,nbmax
      endif
!          call mpi_finalize(ierr)
!          stop
      call realloc_namax_related(nalmax,int(maxb*1.5),iprint)
    endif

    if( nex(kd).gt.1 ) then
      do kdd= -1,0
        ku= 2*kd+kdd
        if( .not. luse_charge ) then
          do i=1,lsb(0,ku)
            j= lsb(i,ku)
            iex= lsex(i,ku)
            ra(1:3,natm+nbnew+i)= ra(1:3,j)
            ra(kd,natm+nbnew+i)= ra(kd,natm+nbnew+i) +iex
            tag(natm+nbnew+i)= tag(j)
          enddo
        else
          do i=1,lsb(0,ku)
            j= lsb(i,ku)
            iex= lsex(i,ku)
            ra(1:3,natm+nbnew+i)= ra(1:3,j)
            ra(kd,natm+nbnew+i)= ra(kd,natm+nbnew+i) +iex
            tag(natm+nbnew+i)= tag(j)
            chg(natm+nbnew+i)= chg(j)
            chi(natm+nbnew+i)= chi(j)
          enddo
        endif
        nbnew= nbnew +lsb(0,ku)
      enddo
    else
      tcom1=mpi_wtime()
      do kdd= -1,0
        ku=2*kd+kdd
        inode=nn(ku)
        nsd=lsb(0,ku)
        call mespasi(inode,myparity(kd),nsd,nrc,1,1,10 &
             ,mpi_md_world)
!---------Store the # of received boundary atoms in LSRC
        lsrc(ku)=nrc

!---------Exchange ra and tag
        if( .not. luse_charge ) then
          do i=1,nsd
            j= lsb(i,ku)
            dbuf(1:3,i)= ra(1:3,j) -sv(1:3,ku)
            dbuf(4,i)  = tag(j)
          enddo
          call mespasd(inode,myparity(kd),dbuf,dbufr,nsd*4,nrc*4,21 &
               ,mpi_md_world)
          do i=1,nrc
            ra(1:3,natm+nbnew+i)= dbufr(1:3,i)
            tag(natm+nbnew+i)   = dbufr(4,i)
          enddo
        else
!.....Exchange ra, tag, chg and chi in case of using charge.
          do i=1,nsd
            j= lsb(i,ku)
            dbuf(1:3,i)= ra(1:3,j) -sv(1:3,ku)
            dbuf(4,i)  = tag(j)
            dbuf(5,i)  = chg(j)
            dbuf(6,i)  = chi(j)
          enddo
          call mespasd(inode,myparity(kd),dbuf,dbufr,nsd*6,nrc*6,21 &
               ,mpi_md_world)
          do i=1,nrc
            ra(1:3,natm+nbnew+i)= dbufr(1:3,i)
            tag(natm+nbnew+i)   = dbufr(4,i)
            chg(natm+nbnew+i)   = dbufr(5,i)
            chi(natm+nbnew+i)   = dbufr(6,i)
          enddo
        endif

        call mpi_barrier(mpi_md_world,ierr)
!---------increase the # of received boundary atoms
        nbnew=nbnew+nrc

200     continue
      enddo
!-------Add the communication time to COMT
      tcom2=MPI_WTIME()
      tcom=tcom+tcom2-tcom1
    endif

!-------Error trap
    nav=natm+nbnew
    call mpi_allreduce(nav,maxna,1,mpi_integer,mpi_max &
         ,mpi_md_world,ierr)
    if (maxna.gt.namax) then
      if (myid.eq.0) then
        write(*,*)'NAMAX overflowed at bacopy'
        write(*,*)'N+NB NAMAX = ',maxna,namax
      endif
      call mpi_finalize(ierr)
      stop
    endif

100 continue
  enddo

!-----num. of received boundary atoms
  nb=nbnew

end subroutine bacopy
!=======================================================================
subroutine bacopy_fixed(tcom,sgm,vol,lsb,lsex,nbmax,ra,namax &
     ,natm,nb,anxi,anyi,anzi,nn,tag,rc,myid_md,myparity,lsrc,sv &
     ,nex &
     ,mpi_md_world,ifcoulomb,chg,chi,boundary)
!-----------------------------------------------------------------------
!  Exchanges boundary-atom data among neighbor nodes: tag and ra
!  This doesnt search using position, just send & recv data of atoms
!    which were listed by 'bacopy'.
!  Different number of data are copied depending on whether 
!    using atomic charges or not.
!-----------------------------------------------------------------------
  use force
  use pmdmpi,only: nx,ny,nz,nid2xyz
  implicit none
  include 'mpif.h'
  integer,intent(in):: namax,nbmax,nn(6),natm,nb &
       ,myid_md,myparity(3),mpi_md_world,nex(3)
  integer,intent(inout):: lsb(0:nbmax,6),lsex(nbmax,6),lsrc(6)
  real(8),intent(in):: sv(3,6),sgm(3,3),vol,anxi,anyi,anzi,rc
  real(8),intent(inout):: tag(namax),ra(3,namax),tcom
  integer,intent(in):: ifcoulomb
  real(8),intent(inout):: chg(namax),chi(namax)
  character(len=3),intent(in):: boundary

!      integer:: status(MPI_STATUS_SIZE)
  integer:: i,j,kd,kdd,ku,ierr,iex,ix,iy,iz
  integer:: inode,nsd,nrc,nbnew
  real(8):: tcom1,tcom2 !,rcv(3),asgm
  real(8),save,allocatable:: dbuf(:,:),dbufr(:,:)
  logical,save:: l1st=.true.

  if( l1st ) then
    if( allocated(dbuf) ) deallocate(dbuf,dbufr)
    if( .not. luse_charge ) then
      allocate(dbuf(4,nbmax),dbufr(4,nbmax))
    else
      allocate(dbuf(6,nbmax),dbufr(6,nbmax))
    endif
    l1st=.false.
  endif

  if( .not. luse_charge ) then
    if( size(dbuf).lt.4*nbmax ) then
      deallocate(dbuf,dbufr)
      allocate(dbuf(4,nbmax),dbufr(4,nbmax))
    endif
  else
    if( size(dbuf).lt.6*nbmax ) then
      deallocate(dbuf,dbufr)
      allocate(dbuf(6,nbmax),dbufr(6,nbmax))
    endif
  endif

  call nid2xyz(myid_md,ix,iy,iz)

  nbnew= 0

!c-----calculate the cut-off lengths
!      do kd=1,3
!        asgm= dsqrt(sgm(1,i)**2 +sgm(2,i)**2 +sgm(3,i)**2)
!        rcv(kd)= rc*asgm/vol
!        nex(kd)= int(rcv(kd)) +1
!      enddo

!-----loop over x, y, & z directions
  do kd=1,3

!-------To calculate the communication time
    tcom1=MPI_WTIME()

    if( nex(kd).gt.1 ) then
      do kdd= -1,0
        ku= 2*kd+kdd
        if( .not. luse_charge ) then
          do i=1,lsb(0,ku)
            j= lsb(i,ku)
            iex= lsex(i,ku)
            ra(1:3,natm+nbnew+i)= ra(1:3,j)
            ra(kd,natm+nbnew+i)= ra(kd,natm+nbnew+i) +iex
            tag(natm+nbnew+i)= tag(j)
          enddo
        else
          do i=1,lsb(0,ku)
            j= lsb(i,ku)
            iex= lsex(i,ku)
            ra(1:3,natm+nbnew+i)= ra(1:3,j)
            ra(kd,natm+nbnew+i)= ra(kd,natm+nbnew+i) +iex
            tag(natm+nbnew+i)= tag(j)
            chg(natm+nbnew+i)= chg(j)
            chi(natm+nbnew+i)= chi(j)
          enddo
        endif
        nbnew= nbnew +lsb(0,ku)
      enddo
    else
      do kdd= -1,0
        ku=2*kd+kdd
        inode=nn(ku)
        nsd=lsb(0,ku)
        nrc=lsrc(ku)

!---------Exchange ra and tag
        if( .not. luse_charge ) then  ! in case of no charge
          do i=1,nsd
            j= lsb(i,ku)
            dbuf(1:3,i)= ra(1:3,j) -sv(1:3,ku)
            dbuf(4,i)  = tag(j)
          enddo
          call mespasd(inode,myparity(kd),dbuf,dbufr,nsd*4,nrc*4,21 &
               ,mpi_md_world)
          do i=1,nrc
            ra(1:3,natm+nbnew+i)= dbufr(1:3,i)
            tag(natm+nbnew+i)   = dbufr(4,i)
          enddo
        else  ! in case of using atomic charge
          do i=1,nsd
            j= lsb(i,ku)
            dbuf(1:3,i)= ra(1:3,j) -sv(1:3,ku)
            dbuf(4,i)  = tag(j)
            dbuf(5,i)  = chg(j)
            dbuf(6,i)  = chi(j)
          enddo
          call mespasd(inode,myparity(kd),dbuf,dbufr,nsd*6,nrc*6,21 &
               ,mpi_md_world)
          do i=1,nrc
            ra(1:3,natm+nbnew+i)= dbufr(1:3,i)
            tag(natm+nbnew+i)   = dbufr(4,i)
            chg(natm+nbnew+i)   = dbufr(5,i)
            chi(natm+nbnew+i)   = dbufr(6,i)
          enddo
        endif

        call MPI_BARRIER(mpi_md_world,ierr)
        nbnew=nbnew +nrc
200     continue
      enddo
    endif

!-------Add the communication time to COMT
    tcom2=MPI_WTIME()
    tcom=tcom+tcom2-tcom1

100 continue
  enddo

end subroutine bacopy_fixed
!=======================================================================
subroutine bacopy_chg_fixed(tcom,lsb,lsex,nbmax,namax &
     ,natm,nb,nn,myid_md,myparity,lsrc &
     ,nex,mpi_md_world,chg,boundary)
!-----------------------------------------------------------------------
!  Exchange only chg of boundary-atom
!  This doesnt search using position, just send & recv data of atoms
!    which were listed by 'bacopy'.
!-----------------------------------------------------------------------
  use pmdmpi,only: nx,ny,nz,nid2xyz
  implicit none
  include 'mpif.h'
  integer,intent(in):: namax,nbmax,nn(6),natm,nb &
       ,myid_md,myparity(3),mpi_md_world,nex(3)
  integer,intent(inout):: lsb(0:nbmax,6),lsex(nbmax,6),lsrc(6)
  real(8),intent(inout):: tcom
  real(8),intent(inout):: chg(namax)
  character(len=3),intent(in):: boundary

  integer:: i,j,kd,kdd,ku,ierr,iex,ix,iy,iz
  integer:: inode,nsd,nrc,nbnew
  real(8):: tcom1,tcom2 !,rcv(3),asgm
  real(8),save,allocatable:: dbuf(:),dbufr(:)
  logical,save:: l1st=.true.

  if( l1st ) then
    allocate(dbuf(nbmax),dbufr(nbmax))
    l1st=.false.
  endif

  call nid2xyz(myid_md,ix,iy,iz)

  nbnew= 0

!c-----calculate the cut-off lengths
!      do kd=1,3
!        asgm= dsqrt(sgm(1,i)**2 +sgm(2,i)**2 +sgm(3,i)**2)
!        rcv(kd)= rc*asgm/vol
!        nex(kd)= int(rcv(kd)) +1
!      enddo

!-----loop over x, y, & z directions
  do kd=1,3

!-------To calculate the communication time
    tcom1=MPI_WTIME()

    if( nex(kd).gt.1 ) then
      do kdd= -1,0
        ku= 2*kd+kdd
        do i=1,lsb(0,ku)
          j= lsb(i,ku)
          iex= lsex(i,ku)
          chg(natm+nbnew+i)= chg(j)
        enddo
        nbnew= nbnew +lsb(0,ku)
      enddo
    else
      do kdd= -1,0
        ku=2*kd+kdd
        inode=nn(ku)
        nsd=lsb(0,ku)
        nrc=lsrc(ku)

!---------Exchange ra and tag
        do i=1,nsd
          j= lsb(i,ku)
          dbuf(i)  = chg(j)
        enddo
        call mespasd(inode,myparity(kd),dbuf,dbufr,nsd,nrc,23 &
             ,mpi_md_world)
        do i=1,nrc
          chg(natm+nbnew+i)   = dbufr(i)
        enddo

        call mpi_barrier(mpi_md_world,ierr)
        nbnew=nbnew +nrc
200     continue
      enddo
    endif

!-------Add the communication time to COMT
    tcom2=MPI_WTIME()
    tcom=tcom+tcom2-tcom1

100 continue
  enddo

end subroutine bacopy_chg_fixed
!=======================================================================
subroutine bamove(tcom,namax,nbmax,natm,ra,va,tag,chg,chi &
     ,anxi,anyi,anzi,myid_md,nn,sv,myparity,mpi_md_world &
     ,boundary)
!-----------------------------------------------------------------------
!  Exchange atoms between neighbor nodes and atomic data as well.
!
!  MVQUE(0:NBMAX,6):
!    MVQUE(0,ku) is the # of to-be-moved atoms to neighbor ku;
!    MVQUE(i,ku) is the adress, in IS, of atom i
!-----------------------------------------------------------------------
  use pmdmpi,only: nid2xyz,nx,ny,nz
  implicit none
  include 'mpif.h'
  integer,intent(in):: namax,nbmax
  integer,intent(in):: myid_md,mpi_md_world,myparity(3),nn(6)
  real(8),intent(in):: anxi,anyi,anzi,sv(3,6)
  integer,intent(inout):: natm
  real(8),intent(inout):: tcom,ra(3,namax),va(3,namax),tag(namax) &
       ,chg(namax),chi(namax)
  character(len=3),intent(in):: boundary

!      integer:: status(MPI_STATUS_SIZE)
  integer:: i,j,ku,kd,kdd,kul,kuh,inode,nsd,nrc,ipt,ierr,is,ix,iy,iz
  integer:: mvque(0:nbmax,6),newim
  real(8):: tcom1,tcom2,xi(3)
  logical,external:: bmv
  real(8),save,allocatable:: dbuf(:,:),dbufr(:,:)
  logical,save:: l1st=.true.

  if( l1st ) then
    allocate(dbuf(9,nbmax),dbufr(9,nbmax))
    l1st=.false.
  endif

  if( size(dbuf).ne.9*nbmax ) then
    deallocate(dbuf,dbufr)
    allocate(dbuf(9,nbmax),dbufr(9,nbmax))
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
!
!        if( kd.eq.1 ) then
!          if( ix.eq.0 .and. mvque(0,kul).gt.0 ) then
!            print *,' myid,kd,ix,mvque(0,kul)=',myid_md
!     &           ,kd,ix,mvque(0,kul)
!          else if( ix.eq.nx-1 .and. mvque(0,kuh).gt.0 ) then
!            print *,' myid,kd,ix,mvque(0,kuh)=',myid_md
!     &           ,kd,ix,mvque(0,kuh)
!          endif
!        endif

!-------Error trap
    if (mvque(0,kul).gt.nbmax) then
      print *,'Buffer overflowed at bamove node',myid_md
      print *,'# in MVQUE=',mvque(0,kul)
      stop
    else if (mvque(0,kuh).gt.nbmax) then
      print *,'Buffer overflowed at bamove node',myid_md
      print *,'# in MVQUE=',mvque(0,kuh)
      stop
    endif

!-------To calculate the communacation time
    tcom1=MPI_WTIME()

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
        dbuf(8,i)  = chg(j)
        dbuf(9,i)  = chi(j)
!-----------Eliminate the record of a moved-out atom
        tag(j)= 0d0
        chg(j)= 0d0
        chi(j)= 0d0
      enddo
      call mespasd(inode,myparity(kd),dbuf,dbufr,9*nsd,9*nrc,71 &
           ,mpi_md_world)
      do i=1,nrc
        ra(1:3,natm+newim+i)= dbufr(1:3,i)
        va(1:3,natm+newim+i)= dbufr(4:6,i)
        tag(natm+newim+i)   = dbufr(7,i)
        chg(natm+newim+i)   = dbufr(8,i)
        chi(natm+newim+i)   = dbufr(9,i)
      enddo

      newim=newim+nrc
      call MPI_BARRIER(mpi_md_world,ierr)

    enddo

    tcom2=MPI_WTIME()
    tcom=tcom+tcom2-tcom1 

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
      chg(ipt)   = chg(i)
      chi(ipt)   = chi(i)
    endif
  enddo
!-----Update # of resident atoms
  natm=ipt

!.....Update max natm globally
!      nalmax = max(natm,nalmax)

  return
end subroutine bamove
!=======================================================================
function bbd(xv,yv,zv,rcav,rcbv,rccv,ku,anxi,anyi,anzi)
!-----------------------------------------------------------------------
!  BBD = .true. if the coordinates are in the boundary to neighbor ku
!-----------------------------------------------------------------------
  implicit real*8(a-h,o-z)
  logical:: bbd

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
    write(*,*)'BBD call is out of range'
  endif
  return
end function bbd
!=======================================================================
function bmv(xv,yv,zv,ku,anxi,anyi,anzi)
!-----------------------------------------------------------------------
!  BMV = .true. if the coordinates should belong to neighbor ku
!------------------------------------------------------------------------
  implicit real*8(a-h,o-z)
  logical bmv

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
    write(*,*)'BMV call is out of range'
  endif
  return
end function bmv
!=======================================================================
subroutine force_isobaric(stgt,ptgt,ah,natm,eki,strs,sgm &
     ,dt,srlx,stbeta,vol,stnsr,mpi_md_world,cpctl)
!
!  Calc. acceralation of h-matrix used in isobaric MD
!  using given stress tensor stgt.
!
  implicit none
  include "mpif.h"
  integer,intent(in):: natm,mpi_md_world
  real(8),intent(in):: stgt(3,3),eki(3,3,natm),strs(3,3,natm),ptgt &
       ,sgm(3,3),dt,vol,srlx,stbeta
  real(8),intent(out):: ah(3,3),stnsr(3,3)
  character(len=*):: cpctl

!.....Max change rate (2%)
  real(8),parameter:: RMAX = 0.02d0

  integer:: i,ixyz,jxyz,ierr
  real(8):: prss,fac,tmp

  call sa2stnsr(natm,strs,eki,stnsr,vol,mpi_md_world)
!!$  print '(a,6f8.3)','stnsr= ',stnsr(1,1),stnsr(2,2),stnsr(3,3) &
!!$       ,stnsr(3,2),stnsr(1,3),stnsr(1,2)

!.....Berendsen for variable-cell
  if( trim(cpctl).eq.'Berendsen' .or. &
       trim(cpctl).eq.'vc-Berendsen' ) then
!.....now ah is scaling factor for h-mat
    ah(1:3,1:3)= 0d0
    ah(1,1)= 1d0
    ah(2,2)= 1d0
    ah(3,3)= 1d0
!.....Limit change rate of h (ah) to RMAX
    do ixyz=1,3
      do jxyz=1,3
        tmp = stbeta*dt/3/srlx*( stgt(ixyz,jxyz)-stnsr(ixyz,jxyz) )
        tmp = min(max(tmp,-RMAX),RMAX)
        ah(ixyz,jxyz) = ah(ixyz,jxyz) -tmp
      enddo
    enddo
!        ah(1:3,1:3)= ah(1:3,1:3)
!     &       -stbeta*dt/3/srlx*( stgt(1:3,1:3)-stnsr(1:3,1:3) )
!.....Berendsen for variable-volume not variable-cell
  else if( trim(cpctl).eq.'vv-Berendsen' ) then
    ah(1:3,1:3)= 0d0
    ah(1,1)= 1d0
    ah(2,2)= 1d0
    ah(3,3)= 1d0
    prss = (stnsr(1,1)+stnsr(2,2)+stnsr(3,3))/3
    fac = 1.0 -stbeta/dt/3/srlx *(ptgt-prss)
!.....Limit change rate of h (ah) to RMAX
    fac = min(max(fac,1d0-RMAX),1d0+RMAX)
    ah(1:3,1:3) = ah(1:3,1:3)*fac
  endif

end subroutine force_isobaric
!=======================================================================
subroutine sa2stnsr(natm,strs,eki,stnsr,vol,mpi_md_world)
!      
!  Take sum of atomic stresses to compute cell stress
!
  implicit none
  include "mpif.h"
  integer,intent(in):: natm,mpi_md_world
  real(8),intent(in):: eki(3,3,natm),strs(3,3,natm),vol
  real(8),intent(out):: stnsr(3,3)

  integer:: i,ixyz,jxyz,ierr
  real(8):: stp(3,3),stk(3,3),stl(3,3),stg(3,3)

  stp(1:3,1:3)= 0d0
  stk(1:3,1:3)= 0d0
  do i=1,natm
    do jxyz=1,3
      do ixyz=1,3
        stk(ixyz,jxyz)=stk(ixyz,jxyz) +2d0*eki(ixyz,jxyz,i)
        stp(ixyz,jxyz)=stp(ixyz,jxyz) +strs(ixyz,jxyz,i)
      enddo
    enddo
  enddo

  stl(1:3,1:3) = stk(1:3,1:3) +stp(1:3,1:3)
!!$  stl(1:3,1:3) = stp(1:3,1:3)
  stg(1:3,1:3)= 0d0
  call mpi_allreduce(stl,stg,9,mpi_real8,mpi_sum &
       ,mpi_md_world,ierr)
  stnsr(1:3,1:3) = stg(1:3,1:3)/vol
  return
end subroutine sa2stnsr
!=======================================================================
subroutine setv(h,natm,tag,va,nspmax,am,tinit,dt)
  implicit none
  include 'mpif.h'
  include 'params_unit.h'
  integer,intent(in):: natm,nspmax
  real(8),intent(in):: tag(natm),am(nspmax),tinit,dt,h(3,3,0:1)
  real(8),intent(out):: va(3,natm)
  integer:: i,l,is
  real(8):: dseed,sumvx,sumvy,sumvz,tmp,facv(nspmax)
  real(8),external:: box_muller

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

  do i=1,natm
    va(1,i)= va(1,i) /h(1,1,0) !*dt
    va(2,i)= va(2,i) /h(2,2,0) !*dt
    va(3,i)= va(3,i) /h(3,3,0) !*dt
  enddo

end subroutine setv
!=======================================================================
subroutine rm_trans_motion(natm,tag,va,nspmax,am &
     ,mpi_md_world,myid_md,iprint)
  implicit none
  include 'mpif.h'
  integer,intent(in):: natm,nspmax,mpi_md_world,myid_md,iprint
  real(8),intent(in):: tag(natm),am(nspmax)
  real(8),intent(out):: va(3,natm)

  integer:: i,is,ierr
  real(8):: sumpx,sumpy,sumpz,amss,amtot,tmp

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
  do i=1,natm
    va(1,i)=va(1,i)-sumpx/amtot
    va(2,i)=va(2,i)-sumpy/amtot
    va(3,i)=va(3,i)-sumpz/amtot
  enddo

  if( myid_md.eq.0 .and. iprint.gt.2 ) then
    write(6,'(a,3es12.4)') ' sumpx,y,z/amtot=' &
         ,sumpx/amtot,sumpy/amtot,sumpz/amtot
  endif

end subroutine rm_trans_motion
!=======================================================================
subroutine vfire(num_fire,alp0_fire,alp_fire,falp_fire,dtmax_fire &
     ,finc_fire,fdec_fire,nmin_fire &
     ,natm,va,aa,myid_md,mpi_md_world,dt,iprint)
  implicit none
  include 'mpif.h'
  integer,intent(in):: natm,myid_md,mpi_md_world,nmin_fire &
       ,iprint
  real(8),intent(in):: aa(3,natm),falp_fire,alp0_fire,dtmax_fire &
       ,finc_fire,fdec_fire
  integer,intent(inout):: num_fire
  real(8),intent(inout):: alp_fire,va(3,natm),dt

  integer:: i,ixyz,ierr
  real(8):: fdotv,vnorm,fnorm
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
  call mpi_allreduce(vnorml,vnorm,1,mpi_real8 &
       ,mpi_sum,mpi_md_world,ierr)
  call mpi_allreduce(fnorml,fnorm,1,mpi_real8 &
       ,mpi_sum,mpi_md_world,ierr)
  call mpi_allreduce(fdotvl,fdotv,1,mpi_real8 &
       ,mpi_sum,mpi_md_world,ierr)
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
      if( iprint.ne.0 .and. myid_md.eq.0 ) then
        write(6,'(a,f10.3,f10.5,es12.4)') ' dt,alp_fire,fdotv = ' &
             ,dt,alp_fire,fdotv
      endif
    endif
  else
    dt = min(dtmax_fire,dt*fdec_fire)
    va(1:3,1:natm) = 0d0
    alp_fire = alp0_fire
    num_fire = 0
    if( iprint.ne.0 .and. myid_md.eq.0 ) then
      write(6,'(a,f10.3,f10.5,es12.4)') ' dt,alp_fire,fdotv = ' &
           ,dt,alp_fire,fdotv
    endif
  endif

!      call mpi_bcast(dt,1,mpi_real8,0,mpi_md_world,ierr)
!      call mpi_bcast(alp_fire,1,mpi_real8,0,mpi_md_world
!     &     ,ierr)
end subroutine vfire
!=======================================================================
subroutine space_decomp(hunit,h,ntot0,tagtot,rtot,vtot &
     ,chgtot,chitot,myid_md,mpi_md_world,nx,ny,nz,nxyz &
     ,rcut,rbuf,iprint)
!
!  Decompose the system and scatter atoms to every process.
!
  use pmdio,only: namax,nbmax,nnmax
  use pmdvars
  implicit none
  include 'mpif.h'
  integer,intent(in):: ntot0,myid_md,mpi_md_world,nx,ny,nz,nxyz &
       ,iprint
  real(8),intent(in):: vtot(3,ntot0),chgtot(ntot0) &
       ,chitot(ntot0),rcut,rbuf
  real(8),intent(in):: hunit,h(3,3,0:1)
  real(8),intent(inout):: rtot(3,ntot0),tagtot(ntot0)

  integer:: istat(mpi_status_size)
  integer:: i,j,ixyz,n,ierr
  integer:: myxt,myyt,myzt,nmin
  real(8):: sxogt,syogt,szogt
  real(8):: t0

  t0 = mpi_wtime()

!.....at 1st call
  if( .not. allocated(ra) ) then
    if( myid_md.eq.0 ) then
!.....wrap atom position inside [0:1)
      do i =1,ntot0
        if( rtot(1,i).lt.0d0 ) then
          rtot(1,i) = rtot(1,i) +1d0
        else if( rtot(1,i).ge.1d0 ) then
          rtot(1,i) = rtot(1,i) -1d0
        endif
        if( rtot(2,i).lt.0d0 ) then
          rtot(2,i) = rtot(2,i) +1d0
        else if( rtot(2,i).ge.1d0 ) then
          rtot(2,i) = rtot(2,i) -1d0
        endif
        if( rtot(3,i).lt.0d0 ) then
          rtot(3,i) = rtot(3,i) +1d0
        else if( rtot(3,i).ge.1d0 ) then
          rtot(3,i) = rtot(3,i) -1d0
        endif
      enddo

!.....count max number of atoms in a node
      nalmax = 0
      nmin = 1000000000
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
               rtot(1,i).lt.sxogt+1d0/nx .and. &
               rtot(2,i).ge.syogt .and. &
               rtot(2,i).lt.syogt+1d0/ny .and. &
               rtot(3,i).ge.szogt .and. &
               rtot(3,i).lt.szogt+1d0/nz .and. &
               tagtot(i).gt.0d0) then
            n = n+1
!.....Set the tag negative if the atom is assigned to a certain cell
            tagtot(i) = -tagtot(i)
          endif
        enddo
        nalmax = max(nalmax,n)
        nmin = min(nmin,n)
      enddo
      namax = max(int(nalmax*1.2),200)
!          nbmax = max(namax*27,nbmax)
      call estimate_nbmax(nalmax,h,nx,ny,nz,rcut,rbuf,nbmax)
      namax = namax +nbmax
      if( iprint.ne.0 ) then
        write(6,'(a,i10)') ' Min number of local atoms = ',nmin
        write(6,'(a,i10)') ' Max number of local atoms = ',nalmax
        write(6,'(a,i10)')   '   nbmax = ',nbmax
        write(6,'(a,i10)')   '   namax = nalmax*1.2 + nbmax  = ' &
             ,namax
        write(6,'(a,i10)')   '   nnmax = ',nnmax
        write(6,'(a,i10,a)') ' Memory for main routine   = ' &
             ,int(dble(namax)*(3*4 +3*3*5 +1 +nnmax/2) &
             *8/1000/1000),' MByte'
      endif
!.....Reset the tags positive
      do i=1,ntot0
        tagtot(i) = abs(tagtot(i))
      enddo
    endif ! myid.eq.0
    call mpi_bcast(namax,1,mpi_integer,0,mpi_md_world,ierr)
    call mpi_bcast(nbmax,1,mpi_integer,0,mpi_md_world,ierr)
!c.....allocate memory for local arrays of length namax
!        allocate(ra(3,namax),va(3,namax),aa(3,namax),ra0(3,namax)
!     &       ,strs(3,3,namax),stt(3,3,namax),tag(namax)
!     &       ,lspr(0:nnmax,namax),epi(namax),eki(3,3,namax)
!     &       ,stp(3,3,namax),stn(3,3,namax)
!     &       ,chg(namax),chi(namax)
!     &       ,lsb(0:nbmax,6),lsex(nbmax,6))
    call alloc_namax_related()
    eki(1:3,1:3,1:namax) = 0d0
  endif

  call mpi_bcast(hunit,1,mpi_real8,0,mpi_md_world,ierr)
  call mpi_bcast(h,9*2,mpi_real8,0,mpi_md_world,ierr)

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
             rtot(1,i).lt.sxogt+1d0/nx .and. &
             rtot(2,i).ge.syogt .and. &
             rtot(2,i).lt.syogt+1d0/ny .and. &
             rtot(3,i).ge.szogt .and. &
             rtot(3,i).lt.szogt+1d0/nz .and. &
             tagtot(i).gt.0d0 ) then
          natm = natm +1
          tag(natm)= tagtot(i)
          ra(1:3,natm)= rtot(1:3,i)
          va(1:3,natm)= vtot(1:3,i)
          chg(natm)= chgtot(i)
          chi(natm)= chitot(i)
          tagtot(i) = -tagtot(i)
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
        call mpi_send(chg,natm,mpi_real8,ixyz,ixyz+5 &
             ,mpi_md_world,ierr)
        call mpi_send(chi,natm,mpi_real8,ixyz,ixyz+6 &
             ,mpi_md_world,ierr)
      endif
    enddo
!        write(6,'(a,f10.3)') ' time space_decomp = ',mpi_wtime() -t0
    tspdcmp = mpi_wtime() -t0
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
    call mpi_recv(chg,natm,mpi_real8,0,myid_md+5 &
         ,mpi_md_world,istat,ierr)
    call mpi_recv(chi,natm,mpi_real8,0,myid_md+6 &
         ,mpi_md_world,istat,ierr)
  endif
  call mpi_barrier(mpi_md_world,ierr)

end subroutine space_decomp
!=======================================================================
subroutine space_comp(ntot0,tagtot,rtot,vtot,atot,epitot &
     ,ekitot,stot,chgtot,chitot,natm,tag,ra,va,aa,epi,eki,strs,chg &
     ,chi,sorg,nxyz,myid_md,mpi_md_world,tspdcmp)
!
!  Opposite to space_decomp, gather atoms from every process
!  to create the total system for output.
!
  use pmdio,only: ntot,alptot
  use util,only: itotOf
  use force,only: loverlay,ol_alphas
  implicit none
  include 'mpif.h'
  integer,intent(in):: ntot0,natm,nxyz,myid_md,mpi_md_world
  real(8),intent(in):: va(3,natm),aa(3,natm),epi(natm),eki(3,3,natm) &
       ,strs(3,3,natm),tag(natm),sorg(3),chg(natm),chi(natm)
  real(8),intent(inout):: ra(3,natm),tspdcmp
  real(8),intent(out):: tagtot(ntot0),rtot(3,ntot0),vtot(3,ntot0) &
       ,atot(3,ntot0),epitot(ntot0),ekitot(3,3,ntot0) &
       ,stot(3,3,ntot0),chgtot(ntot0),chitot(ntot0)
  integer,parameter:: nmpi = 10
  integer:: n0,ixyz,natmt,i,ierr,ntott
  integer:: istat(mpi_status_size),itag
  real(8):: t0
  real(8),allocatable,save:: ratmp(:,:),alptmp(:)
!!$  integer,external:: itotOf

  t0 = mpi_wtime()

  if( .not. allocated(ratmp) ) then
    allocate(ratmp(3,natm),alptmp(natm))
  else if( size(ratmp).lt.3*natm ) then
    deallocate(ratmp,alptmp)
    allocate(ratmp(3,natm),alptmp(natm))
  endif

  if( loverlay ) then
    alptmp(1:natm) = ol_alphas(0,1:natm)
  endif

  if( myid_md.eq.0 ) then
    n0 = natm
    tagtot(1:natm) = tag(1:natm)
    rtot(1:3,1:natm) = ra(1:3,1:natm)
    vtot(1:3,1:natm) = va(1:3,1:natm)
    atot(1:3,1:natm) = aa(1:3,1:natm)
    epitot(1:natm) = epi(1:natm)
    ekitot(1:3,1:3,1:natm) = eki(1:3,1:3,1:natm)
    stot(1:3,1:3,1:natm) = strs(1:3,1:3,1:natm)
    chgtot(1:natm) = chg(1:natm)
    chitot(1:natm) = chi(1:natm)
    if( loverlay ) alptot(1:natm) = alptmp(1:natm)
    ntott = natm
    n0 = n0 +1
    do ixyz=1,nxyz-1
      itag = ixyz*nmpi -nmpi
      call mpi_recv(natmt,1,mpi_integer,ixyz,itag &
           ,mpi_md_world,istat,ierr)
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
      call mpi_recv(chgtot(n0),natmt,mpi_real8 &
           ,ixyz,itag+8,mpi_md_world,istat,ierr)
      call mpi_recv(chitot(n0),natmt,mpi_real8 &
           ,ixyz,itag+9,mpi_md_world,istat,ierr)
      if( loverlay ) then
        call mpi_recv(alptot(n0),natmt,mpi_real8 &
             ,ixyz,itag+10,mpi_md_world,istat,ierr)
      endif
      n0 = n0 + natmt
    enddo
    tspdcmp = tspdcmp +(mpi_wtime()-t0)
!.....Update ntot
    ntot = ntott
    if( ntot.gt.ntot0 ) then
      print *,'ERROR: ntot is greater than ntot0,' &
           //' which should not happen !'
      stop
    endif
!        write(6,'(a,f10.3)') ' time space_comp = ',mpi_wtime() -t0
  else ! myid_md.ne.0
    itag = myid_md*nmpi -nmpi
    call mpi_send(natm,1,mpi_integer,0,itag &
         ,mpi_md_world,ierr)
    call mpi_send(tag,natm,mpi_real8,0,itag+1 &
         ,mpi_md_world,ierr)
!.....positions should be shifted before passing data
    do i=1,natm
      ratmp(1:3,i) = ra(1:3,i) + sorg(1:3)
    enddo
    call mpi_send(ratmp,3*natm,mpi_real8,0,itag+2 &
         ,mpi_md_world,ierr)
    call mpi_send(va,3*natm,mpi_real8,0,itag+3 &
         ,mpi_md_world,ierr)
    call mpi_send(epi,natm,mpi_real8,0,itag+4 &
         ,mpi_md_world,ierr)
    call mpi_send(eki,3*3*natm,mpi_real8,0,itag+5 &
         ,mpi_md_world,ierr)
    call mpi_send(strs,3*3*natm,mpi_real8,0,itag+6 &
         ,mpi_md_world,ierr)
    call mpi_send(aa,3*natm,mpi_real8,0,itag+7 &
         ,mpi_md_world,ierr)
    call mpi_send(chg,natm,mpi_real8,0,itag+8 &
         ,mpi_md_world,ierr)
    call mpi_send(chi,natm,mpi_real8,0,itag+9 &
         ,mpi_md_world,ierr)
    if( loverlay ) then
      call mpi_send(alptmp,natm,mpi_real8,0,itag+10 &
           ,mpi_md_world,ierr)
    endif
  endif
  call mpi_barrier(mpi_md_world,ierr)

end subroutine space_comp
!=======================================================================
subroutine sort_by_tag(natm,tag,ra,va,aa,eki,epi,strs,chg,chi)
!
!  Sort by tag for output.
!  
  use util,only: itotOf
  use pmdio,only: alptot
  use force,only: loverlay
  implicit none
  integer,intent(in):: natm
  real(8),intent(inout):: ra(3,natm),va(3,natm),aa(3,natm) &
       ,eki(3,3,natm),epi(natm),strs(3,3,natm),tag(natm) &
       ,chg(natm),chi(natm)

  integer,allocatable:: itag(:)
  real(8),allocatable:: buf(:,:)
  integer:: i
  integer,save:: nsave = 0
  integer,save:: ndata
!!$  integer,external:: itotOf

  ndata = 31
  if( loverlay ) ndata = ndata +1

  if( .not. allocated(itag) .or. natm.gt.nsave ) then
    if( allocated(itag) ) deallocate(itag,buf)
    nsave = natm
    allocate(itag(natm),buf(ndata,natm))
  endif

  do i=1,natm
    buf(1:3,i)= ra(1:3,i)
    buf(4:6,i)= va(1:3,i)
    buf(7,i)= epi(i)
    buf(8:10,i)= eki(1:3,1,i)
    buf(11:13,i)= eki(1:3,2,i)
    buf(14:16,i)= eki(1:3,3,i)
    buf(17:19,i) = strs(1:3,1,i)
    buf(20:22,i)= strs(1:3,2,i)
    buf(23:25,i)= strs(1:3,3,i)
    buf(26,i)= tag(i)
    buf(27,i)= chg(i)
    buf(28,i)= chi(i)
    buf(29:31,i)= aa(1:3,i)
    if( loverlay ) buf(ndata,i)= alptot(i)
    itag(i)= itotOf(tag(i))
  enddo

!!$  call heapsort_itag(natm,natm,itag,31,buf)
  call qsort_itag(natm,1,natm,itag,ndata,buf)

  do i=1,natm
    ra(1:3,i)= buf(1:3,i)
    va(1:3,i)= buf(4:6,i)
    epi(i)= buf(7,i)
    eki(1:3,1,i)= buf(8:10,i)
    eki(1:3,2,i)= buf(11:13,i)
    eki(1:3,3,i)= buf(14:16,i)
    strs(1:3,1,i)= buf(17:19,i)
    strs(1:3,2,i)= buf(20:22,i)
    strs(1:3,3,i)= buf(23:25,i)
    tag(i)= buf(26,i)
    chg(i)= buf(27,i)
    chi(i)= buf(28,i)
    aa(1:3,i)= buf(29:31,i)
    if( loverlay ) alptot(i)= buf(ndata,i)
  enddo

end subroutine sort_by_tag
!=======================================================================
subroutine error_mpi_stop(cerrmsg)
  use pmdmpi
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
subroutine estimate_nbmax(nalmax,h,nx,ny,nz,rcut,rbuf,nbmax)
  implicit none
  integer,intent(in):: nalmax,nx,ny,nz
  real(8),intent(in):: rcut,rbuf,h(3,3)
  integer,intent(inout):: nbmax

  integer:: nest
  real(8):: alx,aly,alz,area,vol,density

  alx = dsqrt(h(1,1)**2 +h(2,1)**2 +h(3,1)**2)/nx
  aly = dsqrt(h(1,2)**2 +h(2,2)**2 +h(3,2)**2)/ny
  alz = dsqrt(h(1,3)**2 +h(2,3)**2 +h(3,3)**2)/nz
!.....Estimated volume and area, which is not correct in case of
!.....non-orthogonal cell
  vol = alx*aly*alz
!.....Density cannot be obtained from nalmax
!.....because sometimes the system contains vacuum.
!.....Density = 0.1 means 1 atom in 10 Ang^3 which is dense enough.
  density = 0.1d0
  area = alx*aly*2 +alx*alz*2 +aly*alz*2
!.....Estimated number of atoms in the margin region.
!.....Factor 2 is just for a surplus
  nest = density *area*(rcut+rbuf) *2
  nbmax = max(nbmax,nest)

end subroutine estimate_nbmax
!=======================================================================
subroutine alloc_namax_related()
!     
!     Allocated arrays related to NAMAX.
!
  use pmdio,only: namax,nbmax,nnmax
  use pmdvars
  implicit none

  if( allocated(ra) ) deallocate(ra)
  if( allocated(va) ) deallocate(va)
  if( allocated(aa) ) deallocate(aa)
  if( allocated(ra0) ) deallocate(ra0)
  if( allocated(strs) ) deallocate(strs)
  if( allocated(stt) ) deallocate(stt)
  if( allocated(tag) ) deallocate(tag)
  if( allocated(lspr) ) deallocate(lspr)
  if( allocated(ls1nn) ) deallocate(ls1nn)
  if( allocated(epi) ) deallocate(epi)
  if( allocated(eki) ) deallocate(eki)
  if( allocated(stp) ) deallocate(stp)
  if( allocated(stn) ) deallocate(stn)
  if( allocated(chg) ) deallocate(chg)
  if( allocated(chi) ) deallocate(chi)
  if( allocated(lsb) ) deallocate(lsb)
  if( allocated(lsex) ) deallocate(lsex)
  allocate(ra(3,namax),va(3,namax),aa(3,namax),ra0(3,namax) &
       ,strs(3,3,namax),stt(3,3,namax),tag(namax) &
       ,lspr(0:nnmax,namax),ls1nn(0:nnmax,namax) &
       ,epi(namax),eki(3,3,namax) &
       ,stp(3,3,namax),stn(3,3,namax) &
       ,chg(namax),chi(namax) &
       ,lsb(0:nbmax,6),lsex(nbmax,6))
  return
end subroutine alloc_namax_related
!=======================================================================
subroutine realloc_namax_related(newnalmax,newnbmax,iprint)
!     
!     Reallocated namax-related arrays everywhen namax needed to be
!     updated
!
  use pmdio,only: namax,nbmax,nnmax
  use pmdvars
  implicit none
  integer,intent(in):: iprint,newnalmax,newnbmax

  integer:: ierr,newnamax,ndim,l,m,inc
  real(8),allocatable:: arr(:)
  integer,allocatable:: iarr(:)

  if( .not. allocated(ra) ) then
    print *,'Error: Arrays are not allocated yet!'
    call mpi_finalize(ierr)
    stop
  endif

  newnamax = int(newnalmax*1.2) + newnbmax
!.....No need of updating arrays
  if( newnamax .lt. namax ) return

!.....Since the arrays already have data,
!.....those data must be restored after reallocation.
!.....The hard coding here is too messy but I do not know how to avoid..

!.....ra
  ndim = size(ra)
  allocate(arr(ndim))
  call copy_arr(ndim,ra,arr)
  deallocate(ra)
  allocate(ra(3,newnamax))
  call copy_arr(ndim,arr,ra)
  deallocate(arr)

!.....va
  ndim = size(va)
  allocate(arr(ndim))
  call copy_arr(ndim,va,arr)
  deallocate(va)
  allocate(va(3,newnamax))
  call copy_arr(ndim,arr,va)
  deallocate(arr)

!.....aa
  ndim = size(aa)
  allocate(arr(ndim))
  call copy_arr(ndim,aa,arr)
  deallocate(aa)
  allocate(aa(3,newnamax))
  call copy_arr(ndim,arr,aa)
  deallocate(arr)

!.....ra0
  ndim = size(ra0)
  allocate(arr(ndim))
  call copy_arr(ndim,ra0,arr)
  deallocate(ra0)
  allocate(ra0(3,newnamax))
  call copy_arr(ndim,arr,ra0)
  deallocate(arr)

!.....strs
  ndim = size(strs)
  allocate(arr(ndim))
  call copy_arr(ndim,strs,arr)
  deallocate(strs)
  allocate(strs(3,3,newnamax))
  call copy_arr(ndim,arr,strs)
  deallocate(arr)

!.....stt
  ndim = size(stt)
  allocate(arr(ndim))
  call copy_arr(ndim,stt,arr)
  deallocate(stt)
  allocate(stt(3,3,newnamax))
  call copy_arr(ndim,arr,stt)
  deallocate(arr)

!.....tag
  ndim = size(tag)
  allocate(arr(ndim))
  call copy_arr(ndim,tag,arr)
  deallocate(tag)
  allocate(tag(newnamax))
  call copy_arr(ndim,arr,tag)
  deallocate(arr)

!.....lspr
  ndim = size(lspr)
  allocate(iarr(ndim))
  call copy_iarr(ndim,lspr,iarr)
  deallocate(lspr)
  allocate(lspr(0:nnmax,newnamax))
  call copy_iarr(ndim,iarr,lspr)
  deallocate(iarr)


!.....ls1nn
  ndim = size(ls1nn)
  allocate(iarr(ndim))
  call copy_iarr(ndim,ls1nn,iarr)
  deallocate(ls1nn)
  allocate(ls1nn(0:nnmax,newnamax))
  call copy_iarr(ndim,iarr,ls1nn)
  deallocate(iarr)

!.....epi
  ndim = size(epi)
  allocate(arr(ndim))
  call copy_arr(ndim,epi,arr)
  deallocate(epi)
  allocate(epi(newnamax))
  call copy_arr(ndim,arr,epi)
  deallocate(arr)

!.....eki
  ndim = size(eki)
  allocate(arr(ndim))
  call copy_arr(ndim,eki,arr)
  deallocate(eki)
  allocate(eki(3,3,newnamax))
  call copy_arr(ndim,arr,eki)
  deallocate(arr)

!.....stp
  ndim = size(stp)
  allocate(arr(ndim))
  call copy_arr(ndim,stp,arr)
  deallocate(stp)
  allocate(stp(3,3,newnamax))
  call copy_arr(ndim,arr,stp)
  deallocate(arr)

!.....stn
  ndim = size(stn)
  allocate(arr(ndim))
  call copy_arr(ndim,stn,arr)
  deallocate(stn)
  allocate(stn(3,3,newnamax))
  call copy_arr(ndim,arr,stn)
  deallocate(arr)

!.....chg
  ndim = size(chg)
  allocate(arr(ndim))
  call copy_arr(ndim,chg,arr)
  deallocate(chg)
  allocate(chg(newnamax))
  call copy_arr(ndim,arr,chg)
  deallocate(arr)

!.....chi
  ndim = size(chi)
  allocate(arr(ndim))
  call copy_arr(ndim,chi,arr)
  deallocate(chi)
  allocate(chi(newnamax))
  call copy_arr(ndim,arr,chi)
  deallocate(arr)

!.....lsb
  ndim = size(lsb)
  allocate(iarr(ndim))
  call copy_iarr(ndim,lsb,iarr)
  deallocate(lsb)
  allocate(lsb(0:newnbmax,6))
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

  return
end subroutine realloc_namax_related
!=======================================================================
subroutine copy_arr(ndim,srcarr,destarr)
  integer,intent(in):: ndim
  real(8),intent(in):: srcarr(ndim)
  real(8),intent(out):: destarr(ndim)

  destarr(:) = srcarr(:)
  return
end subroutine copy_arr
!=======================================================================
subroutine copy_iarr(ndim,srcarr,destarr)
  integer,intent(in):: ndim
  integer,intent(in):: srcarr(ndim)
  integer,intent(out):: destarr(ndim)

  destarr(:) = srcarr(:)
  return
end subroutine copy_iarr
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:
