c-----------------------------------------------------------------------
c                     Last-modified: <2018-02-05 21:12:50 Ryo KOBAYASHI>
c-----------------------------------------------------------------------
c This file includes core subroutines/functions needed for pmd.
c-----------------------------------------------------------------------
      subroutine pmd_core(hunit,h,ntot,tagtot,rtot,vtot,atot,stot
     &     ,ekitot,epitot,chgtot,chitot,maxstp,nerg,npmd
     &     ,myid_md,mpi_md_world,nodes_md,nx,ny,nz
     &     ,nspmax,am,dt,ciofmt,ifpmd,rc,rbuf,ifdmp,dmp
     &     ,minstp,tinit,tfin,ctctl,ttgt,trlx,ltdst,ntdst,nrmtrans,cpctl
     &     ,stgt,ptgt,srlx,stbeta,strfin,lstrs0,lcellfix,fmv,stnsr,epot
     &     ,ekin,n_conv,ifcoulomb,czload_type,eps_conv,ifsort,iprint
     &     ,nstp_done,lvc,boundary)
      use pmdvars
      use zload
      use force
      use ttm,only: init_ttm,langevin_ttm,output_Te,t_ttm,
     &     calc_Ta,update_Te,assign_atom2cell
      use pmdmpi,only: nid2xyz,xyz2nid
      implicit none
      include "mpif.h"
      include "./params_unit.h"
      integer,intent(in):: ntot,maxstp,nerg,npmd,myid_md,mpi_md_world
     &     ,nspmax,ifpmd,ifdmp,minstp,ntdst,ifsort !,numff
     &     ,iprint,nodes_md,nx,ny,nz,n_conv,nrmtrans
      integer,intent(inout):: ifcoulomb
      integer,intent(out):: nstp_done
      real(8),intent(in):: hunit,tinit,tfin
     &     ,trlx,srlx,stbeta,strfin,dmp,eps_conv
      real(8),intent(inout):: tagtot(ntot),rtot(3,ntot),vtot(3,ntot)
     &     ,atot(3,ntot),stot(3,3,ntot),ekitot(3,3,ntot)
     &     ,epitot(ntot),dt,rc,rbuf,h(3,3,0:1),stnsr(3,3),fmv(3,0:9)
     &     ,epot,ekin,am(nspmax),stgt(3,3),ptgt,ttgt(9)
     &     ,chgtot(ntot),chitot(ntot)
      character,intent(in):: ciofmt*6, cpctl*20, ctctl*20
     &     ,czload_type*5,boundary*3
!      character(len=20),intent(in):: cffs(numff)
      logical,intent(in):: ltdst,lstrs0,lcellfix(3,3)
      logical,intent(inout):: lvc

      integer:: i,j,k,l,m,n,ia,ib,is,ifmv,nave,nspl,i_conv,ierr,nxyz
      integer:: ihour,imin,isec
      real(8):: tmp,hscl(3),aai(3),ami,tave,vi(3),vl(3),epotp,
     &     htmp(3,3),prss
      real(8),external:: box_muller,sprod
      logical:: lconverged = .false.
      logical:: lstrs = .false.
      logical:: lcell_updated = .true.
c.....FIRE variables
      real(8):: alp_fire,fnorm,vnorm,fdotv
      integer:: num_fire
c.....For tensile test of Al nanorod
      real(8):: strnow,ftop,fbot
c-----output file names
      character:: cnum*128, ctmp*128
      integer,external:: itotOf,ifmvOf
      logical:: ltot_updated = .true.

c      logical:: luse_LJ,luse_Ito3_WHe,luse_RK_WHe,luse_RK_FeH, 
c     &     luse_Ramas_FeH,luse_Ackland_Fe,luse_SW_Si,luse_EDIP_Si,
c     &     luse_Brenner,luse_Brenner_vdW,luse_Lu_WHe,luse_Branicio_AlN,
c     &     luse_Mishin_Al,luse_AFS_W,luse_SC_Fe,luse_SM_Al,luse_EAM,
c     &     luse_linreg,luse_NN,luse_Morse,luse_Morse_repul,luse_vcMorse,
c     &     luse_Buckingham,luse_Bonny_WRe,luse_SRIM

      call initialize_pmdvars(nspmax)
      call calc_nfmv(ntot,tagtot,myid_md,mpi_md_world)

c.....multiply damping factor 1st to avoid doing every step
      if( ifdmp.eq.1 ) then
        fmv(1:3,0:9)= fmv(1:3,0:9) *dmp
      endif

c.....Set 0d0 to z of fmv(9) in case of z-loading
      if( czload_type.eq.'atoms' .or. czload_type.eq.'box' ) then
        fmv(3,9) = 0d0
        if( myid_md.eq.0 ) then
          write(6,'(a)') 'Set fmv(3,9)=0d0, because czload_type='
     &         //czload_type
        endif
      endif
      
c-----parallel configuration
      nxyz= nx*ny*nz
      anxi= 1d0/nx
      anyi= 1d0/ny
      anzi= 1d0/nz
c-----error trap
      if(nodes_md.ne.nxyz) then
        write(6,'(a)') " error: nodes_md .ne. nxyz!!"
        write(6,'(a,3i5)') ' myid_md,nodes_md,nxyz='
     &       ,myid_md,nodes_md,nxyz
        print *,'nx,ny,nz =',nx,ny,nz
        call mpi_finalize(ierr)
        stop
      endif
c-----vector node indices: range [0:nx-1]
c      myx=myid_md/(ny*nz)
c      myy=mod(myid_md/nz,ny)
c      myz=mod(myid_md,nz)
      call nid2xyz(myid_md,myx,myy,myz)
c-----reduced node origin
      sorg(1)= anxi*myx
      sorg(2)= anyi*myy
      sorg(3)= anzi*myz

c     check number of parallel node and vc_gaussian
      if( nxyz.gt.1 .and. use_force('long_Coulomb') ) then
        write(6,'(a)') ' [Error] vc_gaussian potential is '
     &       //'not available in parallel.'
        call mpi_finalize(ierr)
        stop
      endif

c-----output every these steps, NOUTERG, NOUTPMD
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
c      do i=1,ntot
c        rtot(1:3,i) = rtot(1:3,i) +0.1d0
c      enddo
c.....perform space decomposition after reading atomic configuration
      call space_decomp(hunit,h,ntot,tagtot,rtot,vtot,chgtot,chitot
     &     ,myid_md,mpi_md_world,nx,ny,nz,nxyz,iprint)
      if( myid_md.eq.0 .and. iprint.ne.0 ) then
        call cell_info(h)
      endif
      print *,'myid,natm=',myid_md,natm
c.....Some conversions
      nsp= 0
      do i=1,natm
        ra(1:3,i)= ra(1:3,i) -sorg(1:3)
c-------scaling
        va(1:3,i)= va(1:3,i) !*dt
c-------species of atom-i
        nsp= max(int(tag(i)),nsp)
        if(nsp.gt.nspmax) then
          write(6,'(a)') " Error: nsp.gt.nspmax!!"
          call mpi_finalize(ierr)
          stop
        endif
      enddo

#ifdef __FITPOT__
c.....check whether order of atoms and total-id of atoms match
      do ia=1,natm
        if( itotOf(tag(ia)).ne.ia ) then
          print *, '[Error] itotOf(tag(ia)).ne.ia !!!'
          print *, '  In case of FITPOT mode, the order of atom'//
     &         ' must be the same of itot of the atom-tag.'
          stop
        endif
      enddo
#endif

c.....NEMD setting
      if( ltdst ) then
        if( mod(ntdst,nx).ne.0 )
     &       call error_mpi_stop('mod(ntdst,nx).ne.0')
        allocate(tdst(ntdst),nadst(ntdst))
      endif


c-----setup
      allocate(acon(nspmax),fack(nspmax))
c      call setup(nspmax,am,acon,fack,dt)
      allocate(fekin(nspmax),fa2v(nspmax))
      call setup2(nspmax,am,fekin,fa2v,dt)
c-----set HI and SGM
      call boxmat(h,hi,ht,g,gi,gt,vol,sgm)
      if( myid_md.eq.0 ) then
        write(6,'(a,f0.2,a)') ' Cell volume = ',vol,' Ang^3'
      endif
c-----ntset
      call ntset(myx,myy,myz,nx,ny,nz,nn,sv,myparity,anxi,anyi,anzi)

c-----get total number of particles
      nspl = nsp
      call mpi_allreduce(nspl,nsp,1,MPI_INTEGER,mpi_max
     &     ,mpi_md_world,ierr)
      if(myid_md.eq.0 .and. iprint.ne.0 ) then
        write(6,*) ''
        write(6,'(a,i0)') " Number of total atoms = ",ntot
        write(6,'(a,i0)') " Number of species     = ",nsp
      endif
c.....get_num_dof is called once in a MD run
      call get_num_dof(natm,tag,fmv,ndof,myid_md,mpi_md_world,iprint)

c.....Setup for FIRE
      if( ifdmp.eq.2 ) then
        alp_fire = alp0_fire
c        am(1:nspmax) = 20.d0
        num_fire = 0
        dtmax_fire = dtmfctr_fire * dt
        if( myid_md.eq.0 .and. iprint.ne.0 ) then
          write(6,'(/,a)') ' FIRE relaxation parameters:'
          write(6,'(a,f0.2,a)') '   Max dt = ',dtmax_fire,' fs'
          write(6,'(a,f0.2)') '   Factor of dt increment = '
     &         ,finc_fire
          write(6,'(a,f0.2)') '   Default alpha = ',alp0_fire
        endif
      endif

c.....Set ifmv of top and bottom atoms for z-loading if needed
      if( czload_type .eq. 'atoms' ) then
        call set_zload_atoms(natm,ra,tag,h,fmv,sorg,strfin,maxstp
     &       ,myid_md,mpi_md_world)
      endif

c.....Set initial temperature if needed
      if( tinit.gt.1d-5 ) then
        call setv(h,natm,tag,va,nspmax,am,tinit,dt)
        if( nrmtrans.gt.0 ) call rm_trans_motion(natm,tag,va,nspmax,am
     &       ,mpi_md_world,myid_md,iprint)
      elseif( abs(tinit).le.1d-5 ) then
        va(1:3,1:natm)= 0d0
      endif

      if( ifdmp.eq.2 ) then
        va(1:3,1:natm) = 0d0
      endif

      if( trim(ctctl).eq.'Langevin' ) then
        tgmm= 1d0/trlx
        do ifmv=1,9
          if( ttgt(ifmv).lt.0d0 ) then
            tfac(ifmv)= -1d0
          else
c            tfac(ifmv)= sqrt(2d0*tgmm*(fkb*ev2j)*ttgt(ifmv)/dt)
c     &           *m2ang/s2fs
            tfac(ifmv)= dsqrt(2d0*tgmm*fkb*ttgt(ifmv)/dt)
c            print *,'ifmv,ttgt,fkb,dt,tfac='
c     &           ,ifmv,ttgt(ifmv),fkb,dt,tfac(ifmv)
          endif
        enddo
        if( myid_md.eq.0 ) then
          print *,''
          print *,'Langevin thermostat parameters:'
          print '(a,f10.2)','   Relaxation time = ',trlx
          do ifmv=1,nfmv
            print '(a,i3,f10.2,es15.4)','   ifmv,ttgt,tfac = '
     &           ,ifmv,ttgt(ifmv),tfac(ifmv)
          enddo
        endif
      else if( trim(ctctl).eq.'ttm' ) then
        call init_ttm(namax,natm,h,dt,myid_md,mpi_md_world,iprint)
      endif

      tcpu1= mpi_wtime()
      tcom = 0d0

      call init_force(namax,natm,nsp,tag,chg,chi,myid_md,mpi_md_world,
     &     iprint,h,rc,lvc,ifcoulomb)
c-----copy RA of boundary atoms
      call check_size_and_parallel(sgm,vol,rc,anxi,anyi,anzi
     &     ,nx,ny,nz,myid_md)
      call bacopy(tcom,sgm,vol,lsb,lsex,nbmax,ra,namax
     &     ,natm,nb,anxi,anyi,anzi,nn,tag,rc,myid_md,myparity
     &     ,lsrc,sv,nex,mpi_md_world,iprint,ifcoulomb,chg,chi,.true.
     &     ,boundary)
c-----Make pair list
      call mk_lspr_para(namax,natm,nbmax,nb,nnmax,tag,ra,rc+rbuf,h,hi
     &     ,anxi,anyi,anzi,lspr,iprint,.true.)
      
c.....Calc forces
      lstrs = lstrs0 .or.
     &     ( trim(cpctl).eq.'Berendsen' .or.
     &       trim(cpctl).eq.'vc-Berendsen' .or.
     &     trim(cpctl).eq.'vv-Berendsen')
      lcell_updated = ( trim(cpctl).eq.'Berendsen' .or.
     &       trim(cpctl).eq.'vc-Berendsen' .or.
     &     trim(cpctl).eq.'vv-Berendsen')
      call get_force(namax,natm,tag,ra,nnmax,aa,strs,chg,chi,stnsr
     &     ,h,hi,tcom,nb,nbmax,lsb,lsex,nex,lsrc,myparity,nn,sv,rc
     &     ,lspr,mpi_md_world,myid_md,epi,epot0,nspmax,acon,lstrs
     &     ,ifcoulomb,iprint,.true.,lvc,lcell_updated,boundary)
      lstrs = .false.
      epot= epot0
      epotp = 0d0
c      write(6,'(a,f12.4,200f10.4)') 'erg,chg after get_force = '
c     &     ,epot,chg(1:natm)

#ifdef __DISL__
      call perf_disl_pos_by_pot(epith,natm,ra,h,epi,sorg
     &     ,nodes_md,myid_md,mpi_md_world,0,21)
#endif

c-----calc kinetic energy
      call get_ekin(namax,natm,va,tag,h,nspmax,fekin,ekin,eki,ekl
     &     ,vmax,mpi_md_world)
      vmaxold=vmax

      if( trim(cpctl).eq.'Berendsen' .or.
     &     trim(cpctl).eq.'vc-Berendsen' ) then
        if(myid_md.eq.0 .and. iprint.ne.0 ) then
          write(6,*) ''
          write(6,'(a)') ' Variable-cell Berendsen is chosen'
          write(6,'(a,6f10.3)') '   Target stress [GPa]: '
     &         ,stgt(1,1),stgt(2,2),stgt(3,3)
     &         ,stgt(2,3),stgt(3,1),stgt(1,2)
        endif
        stgt(1:3,1:3)= stgt(1:3,1:3) *gpa2up
      else if( trim(cpctl).eq.'vv-Berendsen' ) then
        if(myid_md.eq.0 .and. iprint.ne.0 ) then
          write(6,*) ''
          write(6,'(a)') ' Variable-volume Berendsen is chosen'
          write(6,'(a,f0.3,a)')
     &       '   Target pressure = ',ptgt,' GPa'
        endif
        ptgt = ptgt *gpa2up
      endif

      call force_isobaric(stgt,ptgt,ah,natm,eki,strs,sgm
     &     ,dt,srlx,stbeta,vol,stnsr,mpi_md_world,cpctl)
      prss = (stnsr(1,1)+stnsr(2,2)+stnsr(3,3))/3*up2gpa

      istp= 0
      istps= istp
      ts= istps*dt
      te= rbuf/(2d0*vmax/dt) +istp*dt
      istpe= int(te/dt)
c.....Considering vmax==0d0, make istpe small only for the first step
      if( istpe.gt.100 ) istpe=100
c      if(myid_md.eq.0) write(6,'(a,i10)') " istpe=",istpe

      if(myid_md.eq.0 .and. iprint.ne.0 ) then
        write(6,*) ''
        write(6,'(1x,a)') "Initial values:"
        write(6,'(1x,a,f16.5,a)') "  Kinetic energy  = ",ekin,' eV'
        write(6,'(1x,a,f16.5,a)') "  Potential energy= ",epot0,' eV'
        nave= 0
        tave= 0d0
        do ifmv=1,9
          if( ndof(ifmv).eq.0 ) cycle
          temp(ifmv)= ekl(ifmv) /ndof(ifmv) /fkb *2d0
c          print *,' ifmv,ekl,temp = ',ifmv,ekl(ifmv),temp(ifmv)
          nave= nave +ndof(ifmv)
          write(6,'(1x,a,i1,a,f16.5,a)') "  Temperature ",ifmv
     &         ,"   = ",temp(ifmv),' K'
          tave= tave +temp(ifmv)*ndof(ifmv)
        enddo
        tave= tave/(nave-3)
        write(6,'(1x,a,f16.5,a)') "  Temperature     = ",tave,' K'
        if( prss.ge.0d0 ) then
          ctmp = '(tensile)'
        else
          ctmp = '(compressive)'
        endif
        write(6,'(1x,a,f16.5,a)') "  Pressure        = ",
     &       prss,' GPa '//trim(ctmp)
        write(6,'(1x,a,6(1x,f0.3))') '  Stress tensor   ='
     &       ,stnsr(1,1)*up2gpa,stnsr(2,2)*up2gpa,stnsr(3,3)*up2gpa
     &       ,stnsr(2,3)*up2gpa,stnsr(3,1)*up2gpa,stnsr(1,2)*up2gpa
        write(6,*) ''
        if( trim(cpctl).eq.'Berendsen' .or.
     &       trim(cpctl).eq.'vc-Berendsen' .or.
     &       trim(cpctl).eq.'vv-Berendsen') then
          write(6,'(a,i8,2f10.2,es15.6,2es11.3)')
     &         " istp,time,temp,vol,prss="
     &         ,istp,tcpu,tave,epot,vol
     &         ,prss
        else
          if( tave.ge.10000.0 ) then
            write(6,'(a,i8,f11.2,es11.3,es15.6,es11.3)')
     &           " istp,time,temp,epot,prss="
     &           ,istp,tcpu,tave,epot,prss
          else
            write(6,'(a,i8,2f11.2,es15.6,es11.3)')
     &           " istp,time,temp,epot,prss="
     &           ,istp,tcpu,tave,epot,prss
          endif
        endif
      endif

c.....output initial configuration including epi, eki, and strs
c      write(cnum(1:4),'(i4.4)') 0
      write(cnum,'(i0)') 0
      call space_comp(ntot,tagtot,rtot,vtot,atot,epitot,ekitot
     &     ,stot,chgtot,chitot,natm,tag,ra,va,aa,epi,eki,strs
     &     ,chg,chi,sorg,nxyz,myid_md,mpi_md_world,tspdcmp)
      if( myid_md.eq.0 ) then
        if( ifsort.eq.1 ) call sort_by_tag(ntot,tagtot,rtot,vtot
     &       ,ekitot,epitot,stot,chgtot,chitot)
        if( ifpmd.eq.1 ) then  ! pmd format
          if( trim(ciofmt).eq.'bin' .or. trim(ciofmt).eq.'binary' )
     &         then
            call write_pmdtot_bin(20,"pmd_"//trim(cnum))
          elseif( trim(ciofmt).eq.'ascii' ) then
            call write_pmdtot_ascii(20,"pmd_"//trim(cnum))
          endif
        else if( ifpmd.eq.2 ) then ! LAMMPS-dump format
          call write_dump(20,'dump_'//trim(cnum))
        endif
      endif
      
      if( trim(ctctl).eq.'ttm' ) then
        call assign_atom2cell(namax,natm,ra,sorg)
        call calc_Ta(namax,natm,eki,istp,myid_md,mpi_md_world)
        call output_Te(0,myid_md,iprint)
      endif

c-----initialize the counter for output
      iocntpmd=0
      iocnterg=0

      if(myid_md.eq.0) then
c.....write out energies
        open(ioerg,file="out.erg",status='replace')
        write(ioerg,'(a)') '# istp,  etot-epot0[eV],  ekin,'
     &       //'  epot-epot0,  temp[K],  vol[Ang^3],  pressure[GPa]'
        write(ioerg,'(a,es16.7e3,a)') '#  Epot0 =',epot0,' [eV]'
        if( tave.gt.10000d0) then
          write(ioerg,'(i10,3es16.7e3,es11.3,2es16.7e3)') istp
     &         ,ekin+epot0-epot0,ekin,epot0-epot0,tave
     &         ,vol
     &         ,(stnsr(1,1)+stnsr(2,2)+stnsr(3,3))/3*up2gpa
        else
          write(ioerg,'(i10,3es16.7e3,f11.2,2es16.7e3)') istp
     &         ,ekin+epot0-epot0,ekin,epot0-epot0,tave
     &         ,vol
     &         ,(stnsr(1,1)+stnsr(2,2)+stnsr(3,3))/3*up2gpa
        endif
        call flush(ioerg)
c.....write out temperatures
        open(iotemp,file='out.temperature',status='replace')
        write(iotemp,'(a)') '# istp, temperature[0-9]'
        ediff(1:9)= 0d0
        write(iotemp,'(i10,18es16.7e3)') istp,temp(1:9),ediff(1:9)
        call flush(iotemp)
        
c.....write out tensile forces
        if(czload_type.eq.'atoms' .or. czload_type.eq.'box') then
          open(iostrs,file='out.zload',status='replace')
          write(iostrs,'(a)') '#  istp,   strnow,      ftop,     fbot'
        endif
      endif

c.....temperature distribution along x
      if( ltdst ) then
        if( myid_md.eq.0 ) then
          open(iotdst,file='out.temp-dist')
          write(iotdst,'(a)') '#   x-pos,   temp   '
        endif
        call calc_temp_dist(iotdst,ntdst,tdst,nadst,natm,ra,eki
     &     ,istp,nouterg,myid_md,mpi_md_world,sorg)
      endif

c.....Set al
      al(1)= h(1,1,0)
      al(2)= h(2,2,0)
      al(3)= h(3,3,0)

      
      simtime = 0d0
      i_conv = 0
      lconverged = .false.
c-----velocity-Verlet loop starts---------------------------------------
      do istp=1,maxstp

c.....In case of isobaric MD, lstrs has to be always TRUE.
        if( trim(cpctl).eq.'Berendsen' .or.
     &       trim(cpctl).eq.'vc-Berendsen' .or.
     &       trim(cpctl).eq.'vv-Berendsen' .or. lstrs0 ) then
          lstrs = .true.
        else
          lstrs = .false.
        endif

c-------first kick of velocities
        do i=1,natm
          is = int(tag(i))
          va(1:3,i)=va(1:3,i) +aa(1:3,i)*fa2v(is)*dt
        enddo

        if( ifdmp.eq.2 ) then
          call vfire(num_fire,alp0_fire,alp_fire,falp_fire,dtmax_fire
     &         ,finc_fire,fdec_fire,nmin_fire
     &         ,natm,va,aa,myid_md,mpi_md_world,dt,iprint)
c.....Update dt-related variables
c          call setup(nspmax,am,acon,fack,dt)
          call setup2(nspmax,am,fekin,fa2v,dt)
        endif

c-------multiply fmv or damping
        do i=1,natm
          l= int(mod(tag(i)*10,10d0))
          va(1:3,i)=va(1:3,i) *fmv(1:3,l)
        enddo

c-------update positions
        do i=1,natm
          is = int(tag(i))
          ra(1:3,i)=ra(1:3,i) +va(1:3,i)*dt
        enddo
        ltot_updated = .false.

        if( czload_type.eq.'atoms' ) then
          call zload_atoms(natm,ra,tag,maxstp,strfin,strnow
     &         ,sorg,myid_md,mpi_md_world)
        else if( czload_type.eq.'box' ) then
          call zload_box(natm,maxstp,istp,dt,strfin,strnow,h,myid_md)
        endif

        if( trim(cpctl).eq.'Berendsen' .or.
     &       trim(cpctl).eq.'vc-Berendsen' .or.
     &       trim(cpctl).eq.'vv-Berendsen' ) then
c     h(1:3,1:3,0)= matmul(ah,h(1:3,1:3,0))
          htmp(1:3,1:3) = matmul(ah,h(1:3,1:3,0))
          do i=1,3
            do j=1,3
              if( .not. lcellfix(i,j) ) then
                h(i,j,0) = htmp(i,j)
              else
                h(i,j,0) = h(i,j,0)
              endif
            enddo
          enddo
c.....Reset matrices
          call boxmat(h,hi,ht,g,gi,gt,vol,sgm)
        endif
        simtime = simtime +dt

c.....Update te considering change of vmax
        te=(te-ts)*vmaxold/vmax +ts
        istpe= int(te/dt)

c.....Update pair list and boundary atoms
c.....if making new pair list is needed.
        if( istp.ge.istpe .or.
     &       (ifpmd.gt.0.and.mod(istp,noutpmd).eq.0) ) then
c.....Move atoms that cross the boundary
          call bamove(tcom,namax,nbmax,natm,ra,va,tag,chg,chi
     &         ,anxi,anyi,anzi,myid_md,nn,sv,myparity,mpi_md_world
     &         ,boundary)
c.....Copy RA of boundary atoms
          call bacopy(tcom,sgm,vol,lsb,lsex,nbmax,ra,namax,natm,nb
     &         ,anxi,anyi,anzi,nn,tag,rc,myid_md,myparity,lsrc,sv,nex
     &         ,mpi_md_world,iprint,ifcoulomb,chg,chi,.false.,boundary)
c.....Make pair list
          call mk_lspr_para(namax,natm,nbmax,nb,nnmax,tag,ra,rc+rbuf,h
     &         ,hi,anxi,anyi,anzi,lspr,iprint,.false.)
          ts= istp*dt
          te=rbuf/(2d0*vmax/dt) +ts
          istpe= int(te/dt)
c          if(myid_md.eq.0) write(6,'(a,i10)')
c     &         " Next update step will be ",istpe
        else
c.....Copy RA of boundary atoms determined by 'bacopy'
          call bacopy_fixed(tcom,sgm,vol,lsb,lsex,nbmax,ra,namax
     &         ,natm,nb,anxi,anyi,anzi,nn,tag,rc,myid_md,myparity,lsrc
     &         ,sv,nex,mpi_md_world,ifcoulomb,chg,chi,boundary)
        endif
        
        if(ifpmd.gt.0.and. mod(istp,noutpmd).eq.0 )then
          lstrs = lstrs0
        endif
c-------Calc forces
        call get_force(namax,natm,tag,ra,nnmax,aa,strs,chg,chi,stnsr
     &       ,h,hi,tcom,nb,nbmax,lsb,lsex,nex,lsrc,myparity,nn,sv,rc
     &       ,lspr,mpi_md_world,myid_md,epi,epot,nspmax,acon,lstrs
     &       ,ifcoulomb,iprint,.false.,lvc,lcell_updated,boundary)
        lstrs = .false.

c.....Second kick of velocities
        if( trim(ctctl).eq.'Langevin' ) then
c.....Langevin thermostat with Mannella integrator
          hscl(1:3)= 0d0
          do l=1,3
            hscl(l)= dsqrt(h(1,l,0)**2 +h(2,l,0)**2 +h(3,l,0)**2)
          enddo
          do i=1,natm
c            ifmv= int(mod(tag(i)*10,10d0))
            ifmv = ifmvOf(tag(i))
            is = int(tag(i))
            if( ifmv.eq.0 .or. tfac(ifmv).lt.0d0 ) then
              va(1:3,i)=va(1:3,i) +aa(1:3,i)*fa2v(is)*dt
            else
              va(1:3,i)=va(1:3,i) +aa(1:3,i)*fa2v(is)*dt
              ami= am(is)
              aai(1:3)= 0d0
              tmp= tfac(ifmv)*dsqrt(ami) !/sqrt(ami*amu2kg) !*dt
c              if( i.eq.1 ) print *,'ifmv,tfac,ami,tgmm,tmp='
c     &             ,ifmv,tfac(ifmv),ami,tgmm,tmp
              do l=1,3
                aai(l)= -va(l,i)*tgmm*ami +tmp*box_muller()/hscl(l)
              enddo
c.....To compensate the factor 1/2 in fa2v, multiply 2 here.
              va(1:3,i)= va(1:3,i) +aai(1:3)*fa2v(is)*dt *2d0
c.....accumulate energy difference
              vi(1:3)= h(1:3,1,0)*va(1,i)
     &             +h(1:3,2,0)*va(2,i)
     &             +h(1:3,3,0)*va(3,i)
              vl(1:3)= h(1:3,1,0)*aai(1)
     &             +h(1:3,2,0)*aai(2)
     &             +h(1:3,3,0)*aai(3)
              ediff(ifmv)= ediff(ifmv) +fekin(is)
     &             *(2d0*sprod(3,vi,vl)+sprod(3,vl,vl))
            endif
          enddo
          if( mod(istp,nrmtrans).eq.0 ) then
            call rm_trans_motion(natm,tag,va,nspmax,am
     &           ,mpi_md_world,myid_md,0)
          endif
        else if( trim(ctctl).eq.'ttm' ) then
          call assign_atom2cell(namax,natm,ra,sorg)
          call langevin_ttm(namax,natm,va,aa,tag,eki,am,h
     &         ,nspmax,fa2v,fekin,ediff,dt)
        else
          do i=1,natm
            is = int(tag(i))
            va(1:3,i)=va(1:3,i) +aa(1:3,i)*fa2v(is)*dt
c            if( i.eq.52 ) then
c              print *,'i,is,aa(1:3,i),va(1:3,i)=',
c     &             i,is,aa(1:3,i),va(1:3,i)
c            endif
          enddo
        endif

        

c.....Calc kinetic energy
        vmaxold= vmax
        call get_ekin(namax,natm,va,tag,h,nspmax,fekin,ekin,eki,ekl
     &       ,vmax,mpi_md_world)

        if( trim(ctctl).eq.'ttm' ) then
          call calc_Ta(namax,natm,eki,istp,myid_md,mpi_md_world)
          call update_Te(istp*dt,myid_md,mpi_md_world)
        endif

        call force_isobaric(stgt,ptgt,ah,natm,eki,strs,sgm
     &       ,dt,srlx,stbeta,vol,stnsr,mpi_md_world,cpctl)
        prss = (stnsr(1,1)+stnsr(2,2)+stnsr(3,3))/3*up2gpa
c-------temperature control by Berendsen thermostat
        if( trim(ctctl).eq.'Berendsen' ) then
          tfac(1:9)= 0d0
c.....if final temperature is assigned,
c.....target temperatures are forced to be set intermediate temperatures
          if( tfin.gt.0d0 ) then
            ttgt(1:9) = tinit +(tfin-tinit)*istp/maxstp
          endif
          do ifmv=1,9
            if(ndof(ifmv).le.0 .or. ttgt(ifmv).lt.0d0 ) cycle
            temp(ifmv)= ekl(ifmv) *2d0 /fkb /ndof(ifmv)
            if( abs(ttgt(ifmv)-temp(ifmv))/temp(ifmv).gt.100d0 ) then
              tfac(ifmv)= dsqrt(1d0 +dt/trlx*100d0 )
            else
              tfac(ifmv)= dsqrt(1d0 +dt/trlx*(ttgt(ifmv)-temp(ifmv))
     &             /temp(ifmv))
            endif
          enddo
          do i=1,natm
c            ifmv= int(mod(tag(i)*10,10d0))
            ifmv = ifmvOf(tag(i))
            if( ifmv.le.0 .or. ttgt(ifmv).lt.0d0 ) cycle
            va(1:3,i)= va(1:3,i) *tfac(ifmv)
c.....accumulate energy difference
            is= int(tag(i))
            ediff(ifmv)= ediff(ifmv) +(tfac(ifmv)**2-1d0)*fekin(is)
          enddo
          if( mod(istp,nrmtrans).eq.0 ) then
            call rm_trans_motion(natm,tag,va,nspmax,am
     &           ,mpi_md_world,myid_md,0)
          endif
        endif

c.....temperature distribution along x
        if( ltdst ) call calc_temp_dist(iotdst,ntdst,tdst,nadst,natm,ra
     &       ,eki,istp,nouterg,myid_md,mpi_md_world,sorg)

cc.....Eliminate translational motion
c        if( mod(istp,100).eq.0 ) then
c          call rm_trans_motion()
c        endif

c-------output energies
        if(mod(istp,nouterg).eq.0) then
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

          if(czload_type.eq.'atoms') then
            call get_forces_on_base(natm,ra,aa,tag,h,acon(1),ftop,fbot
     &           ,sorg,myid_md,mpi_md_world)
          endif

          ediff(1:9)= ediff(1:9) /nouterg
          ediff0(1:9)= 0d0
          call mpi_reduce(ediff,ediff0,9
     &         ,mpi_double_precision,mpi_sum,0,mpi_md_world,ierr)
          
          if(myid_md.eq.0) then
c.....write energies
            if( tave.gt.10000d0 ) then
              write(ioerg,'(i10,3es16.7e3,es11.3,2es16.7e3)') istp
     &             ,ekin+epot-epot0,ekin,epot-epot0,tave
     &             ,vol,prss
            else
              write(ioerg,'(i10,3es16.7e3,f11.2,2es16.7e3)') istp
     &             ,ekin+epot-epot0,ekin,epot-epot0,tave
     &             ,vol,prss
            endif
            call flush(ioerg)
c.....write temperature
            write(iotemp,'(i10,18es16.7e3)') istp,temp(1:9),ediff0(1:9)
            call flush(iotemp)
            ediff(1:9)= 0d0

            if(czload_type.eq.'atoms') then
              write(iostrs,'(i8,3es15.7)') istp,strnow,ftop,fbot
              call flush(iostrs)
            endif
          endif   ! myid_md.eq.0
          
c---------output step, time, and temperature
          tcpu= mpi_wtime() -tcpu1
          if(myid_md.eq.0) then
            if( trim(cpctl).eq.'Berendsen' .or.
     &           trim(cpctl).eq.'vc-Berendsen' .or.
     &           trim(cpctl).eq.'vv-Berendsen' .and.
     &           iprint.ne.0 ) then
              write(6,'(a,i8,2f10.2,es15.6,2es11.3)')
     &             " istp,time,temp,epot,vol,prss="
     &             ,istp,tcpu,tave,epot,vol,prss
              write(6,'(a)') ' H-matrix:' !,h(1:3,1:3,0)
              write(6,'(3f15.7)') h(1:3,1:3,0)
              
              write(6,'(a,6f10.3)') ' Stress (GPa):'
     &             ,stnsr(1,1)*up2gpa ,stnsr(2,2)*up2gpa
     &             ,stnsr(3,3)*up2gpa ,stnsr(2,3)*up2gpa
     &             ,stnsr(1,3)*up2gpa ,stnsr(1,2)*up2gpa
            elseif( iprint.ne.0 ) then
              if( tave.ge.10000.0 ) then
                write(6,'(a,i8,f11.2,es11.3,es15.6,es11.3)')
     &               " istp,time,temp,epot,prss="
     &               ,istp,tcpu,tave,epot,prss
              else
                write(6,'(a,i8,2f11.2,es15.6,es11.3)')
     &               " istp,time,temp,epot,prss="
     &               ,istp,tcpu,tave,epot,prss
              endif
            endif
          endif

          if( trim(ctctl).eq.'ttm' ) then
            call output_Te(istp,myid_md,iprint)
          endif
            
#ifdef __DISL__
c.....Output disl core pos
          call perf_disl_pos_by_pot(epith,natm,ra,h,epi,sorg
     &         ,nodes_md,myid_md,mpi_md_world,iocnterg,21)
#endif
        endif ! energy output

c.....check convergence criteria if it is dumped MD
        if( ifdmp.gt.0 .and. epot-epotp.le.0d0 .and.
     &       abs(epot-epotp).lt.eps_conv .and. istp.gt.minstp ) then
c          print *,'ifdmp.gt.0 = ',ifdmp.gt.0
c          print *,'epot-epotp.le.0d0 = ',epot-epotp.le.0d0
c          print *,'abs(epot-epotp).lt.eps_conv = ',abs(epot-epotp)
c     &         .lt.eps_conv
c          print *,'istp.gt.minstp = ',istp.gt.minstp
          i_conv = i_conv + 1
          if( i_conv.ge.n_conv ) then
            if( myid_md.eq.0 .and. iprint.ne.0 ) then
              write(6,'(/,a,i6,a,es10.3,a,i3,a)')
     &             ' Damped MD converged with ',istp,
     &             ' steps, since ediff < ',
     &             eps_conv,', ',n_conv,' times'
              write(6,'(a,3es20.10)') '  epot,epotp,epot-epotp = '
     &             ,epot,epotp,epot-epotp
            endif
            lconverged = .true.
          endif
        else
          epotp = epot
          i_conv = 0
        endif

c-------write the particle positions
        if(ifpmd.gt.0.and.
     &       (mod(istp,noutpmd).eq.0.or.lconverged) )then
c---------decide pmd-file name
          iocntpmd=iocntpmd+1
c          write(cnum(1:4),'(i4.4)') iocntpmd
          write(cnum,'(i0)') istp
c          call system("mkdir -p "//cnum)
          call space_comp(ntot,tagtot,rtot,vtot,atot,epitot,ekitot
     &         ,stot,chgtot,chitot,natm,tag,ra,va,aa,epi,eki,strs
     &         ,chg,chi,sorg,nxyz,myid_md,mpi_md_world,tspdcmp)
          ltot_updated = .true.
          if( myid_md.eq.0 ) then
            if( ifsort.eq.1 ) call sort_by_tag(ntot,tagtot,rtot,vtot
     &           ,ekitot,epitot,stot,chgtot,chitot)
            if( ifpmd.eq.1 ) then  ! pmd format
              if( trim(ciofmt).eq.'bin' .or. trim(ciofmt).eq.'binary' )
     &             then
                call write_pmdtot_bin(20,"pmd_"//trim(cnum))
              elseif( trim(ciofmt).eq.'ascii' ) then
                call write_pmdtot_ascii(20,"pmd_"//trim(cnum))
              endif
            else if( ifpmd.eq.2 ) then  ! LAMMPS-dump format
              call write_dump(20,'dump_'//trim(cnum))
            endif
          endif
        endif

        if( lconverged ) exit
      enddo ! end of istp

      if( .not. ltot_updated ) then
        call space_comp(ntot,tagtot,rtot,vtot,atot,epitot,ekitot
     &       ,stot,chgtot,chitot,natm,tag,ra,va,aa,epi,eki,strs
     &       ,chg,chi,sorg,nxyz,myid_md,mpi_md_world,tspdcmp)
      endif
      
      nstp_done = istp
      tcpu2= mpi_wtime()

      if(myid_md.eq.0) then
        close(ioerg)
        close(iotemp)
        if(czload_type.eq.'atoms' .or. czload_type.eq.'box') then
          close(iostrs)
        endif
        if( ltdst ) close(iotdst)
      endif
      
      tcpu= tcpu2 -tcpu1
      if( myid_md.eq.0 .and. iprint.ne. 0 ) then
        write(6,*) ''
        write(6,'(1x,a)') "Final values:"
        write(6,'(1x,a,f16.5,a)') "  Kinetic energy  = ",ekin,' eV'
        write(6,'(1x,a,f16.5,a)') "  Potential energy= ",epot,' eV'
        nave= 0
        tave= 0d0
        do ifmv=1,9
          if( ndof(ifmv).le.0 ) cycle
          temp(ifmv)= ekl(ifmv) *2d0 /fkb /ndof(ifmv)
          nave= nave +ndof(ifmv)
          write(6,'(1x,a,i1,a,f16.5,a)') "  Temperature ",ifmv
     &         ,"   = ",temp(ifmv),' K'
          tave= tave +temp(ifmv)*ndof(ifmv)
        enddo
        tave= tave/(nave-3)
        write(6,'(1x,a,f16.5,a)') "  Temperature     = ",tave,' K'
        if( prss.ge.0d0 ) then
          ctmp = '(tensile)'
        else
          ctmp = '(compressive)'
        endif
        write(6,'(1x,a,f16.5,a)') "  Pressure        = ",
     &       prss,' GPa '//trim(ctmp)
        if( trim(cpctl).eq.'Berendsen' .or.
     &       trim(cpctl).eq.'vc-Berendsen' .or.
     &       trim(cpctl).eq.'vv-Berendsen' ) then
          call cell_info(h)
          write(6,'(a,f16.5,a)')  '  Cell volume     = ',vol,' Ang^3'
        endif
        write(6,*) ''
        write(6,'(1x,a,f10.2)') "Time for space decomp = ",tspdcmp
        write(6,'(1x,a,f10.2)') "Time for comm         = ",tcom
        if( trim(ctctl).eq.'ttm' ) then
          write(6,'(1x,a,f10.2)') "Time for ttm          = ",t_ttm
        endif

        ihour = int(tcpu/3600)
        imin  = int((tcpu-ihour*3600)/60)
        isec  = int(tcpu -ihour*3600 -imin*60)
        write(6,'(1x,a,f10.2,a,i3,"h",i2.2,"m",i2.2,"s")')
     &       "Time                  = ",tcpu,
     &       " sec  = ",ihour,imin,isec
      endif

c.....revert forces to the unit eV/A before going out pmd_core
      if( myid_md.eq.0 ) then
        do i=1,ntot
          is= int(tagtot(i))
          aai(1:3)= h(1:3,1,0)*atot(1,i)
     &         +h(1:3,2,0)*atot(2,i)
     &         +h(1:3,3,0)*atot(3,i)
c          atot(1:3,i) = aai(1:3)/acon(is)
        enddo
      endif

c.....deallocate all the arrays allocated in pmd_core
      if( ltdst ) then
        deallocate(tdst,nadst)
      endif
      deallocate(acon,fack)
      deallocate(fekin,fa2v)
      deallocate(ra,va,aa,ra0,strs,stt,tag,lspr
     &     ,epi,eki,stp,stn,lsb,lsex,chg,chi)

      end subroutine pmd_core
c=======================================================================
      subroutine one_shot(hunit,h,ntot,tagtot,rtot,vtot,atot,stot
     &     ,ekitot,epitot,chgtot,chitot
     &     ,myid_md,mpi_md_world,nodes_md,nx,ny,nz
     &     ,nspmax,am,dt,rc,rbuf,stnsr,epot
     &     ,ekin,ifcoulomb,lvc,iprint,lcalcgrad,ndimp,gwe
     &     ,gwf,gws,lematch,lfmatch,lsmatch,boundary)
c
c  In case that only one shot force calculation is required,
c  especially called from fitpot.
c
      use pmdvars
      use force
      use Morse,only: gradw_Morse,gradw_vcMorse
      implicit none
      include "mpif.h"
      include "./params_unit.h"
      integer,intent(in):: ntot,myid_md,mpi_md_world
     &     ,nspmax  !,numff
     &     ,iprint,nodes_md,nx,ny,nz
      integer,intent(inout):: ifcoulomb
      real(8),intent(in):: hunit
      real(8),intent(inout):: tagtot(ntot),rtot(3,ntot),vtot(3,ntot)
     &     ,atot(3,ntot),stot(3,3,ntot),ekitot(3,3,ntot)
     &     ,epitot(ntot),rc,rbuf,h(3,3,0:1),stnsr(3,3)
     &     ,epot,ekin,am(nspmax),dt
     &     ,chgtot(ntot),chitot(ntot)
c      character(len=20),intent(in):: cffs(numff)
      logical,intent(in):: lcalcgrad
      integer,intent(in):: ndimp
      real(8),intent(inout):: gwe(ndimp),gwf(ndimp,3,ntot),gws(ndimp,6)
      logical,intent(inout):: lvc
      logical,intent(in):: lematch,lfmatch,lsmatch
      character(len=3),intent(in):: boundary

      integer:: i,ierr,is,nspl,nxyz
      real(8):: aai(3),epott
      logical:: lstrs = .false.
      logical:: lcell_updated = .false.

c      print *,'one_shot: 01'
      call initialize_pmdvars(nspmax)

      nxyz = nx*ny*nz
      anxi= 1d0/nx
      anyi= 1d0/ny
      anzi= 1d0/nz
c-----vector node indices: range [0:nx-1]
      myx=myid_md/(ny*nz)
      myy=mod(myid_md/nz,ny)
      myz=mod(myid_md,nz)
c-----reduced node origin
      sorg(1)= anxi*myx
      sorg(2)= anyi*myy
      sorg(3)= anzi*myz
c.....perform space decomposition after reading atomic configuration
c      write(6,'(a,30es10.2)') 'chgtot before space_decomp = ',
c     &     chgtot(1:ntot)
      if( iprint.ne.0 ) then
        write(6,'(a,i8)') ' number of total atoms = ',ntot
        write(6,'(a)') " h-matrix:"
        write(6,'("   | ",3f12.3," |")') h(1,1:3,0)
        write(6,'("   | ",3f12.3," |")') h(2,1:3,0)
        write(6,'("   | ",3f12.3," |")') h(3,1:3,0)
        write(6,'(a)') " which means the following lattice vectors:"
        write(6,'(a,"[ ",3f12.3," ]")') '   a = ',h(1:3,1,0)
        write(6,'(a,"[ ",3f12.3," ]")') '   b = ',h(1:3,2,0)
        write(6,'(a,"[ ",3f12.3," ]")') '   c = ',h(1:3,3,0)
        write(6,*) ''
      endif
c      print *,'one_shot: 02'
      call space_decomp(hunit,h,ntot,tagtot,rtot,vtot,chgtot,chitot
     &     ,myid_md,mpi_md_world,nx,ny,nz,nxyz,iprint)
c      write(6,'(a,200es10.2)') 'chg after space_decomp = ',chg(1:natm)

c      print *,'one_shot: 03'
c.....Some conversions
      nsp= 0
      do i=1,natm
c-------species of atom-i
        nsp= max(int(tag(i)),nsp)
      enddo
c-----setup
      allocate(acon(nspmax),fack(nspmax))
c      call setup(nspmax,am,acon,fack,dt)
      allocate(fekin(nspmax),fa2v(nspmax))
      call setup2(nspmax,am,fekin,fa2v,dt)
c-----set HI and SGM
      call boxmat(h,hi,ht,g,gi,gt,vol,sgm)
c-----ntset
      call ntset(myx,myy,myz,nx,ny,nz,nn,sv,myparity,anxi,anyi,anzi)
c-----get total number of species
      nspl = nsp
      call mpi_allreduce(nspl,nsp,1,mpi_integer,mpi_max
     &     ,mpi_md_world,ierr)

      tcpu1= mpi_wtime()
      tcom = 0d0

c      print *,'one_shot: 04'
      if( iprint.ge.10 ) print *,'init_force,myid_md=',myid_md
      call init_force(namax,natm,nsp,tag,chg,chi,myid_md,mpi_md_world,
     &     iprint,h,rc,lvc,ifcoulomb)

c      print *,'one_shot: 05'
c-----copy RA of boundary atoms
      call check_size_and_parallel(sgm,vol,rc,anxi,anyi,anzi
     &     ,nx,ny,nz,myid_md)
      call bacopy(tcom,sgm,vol,lsb,lsex,nbmax,ra,namax
     &     ,natm,nb,anxi,anyi,anzi,nn,tag,rc,myid_md,myparity
     &     ,lsrc,sv,nex,mpi_md_world,iprint,ifcoulomb,chg,chi,.true.
     &     ,boundary)
c-----Make pair list
      call mk_lspr_para(namax,natm,nbmax,nb,nnmax,tag,ra,rc+rbuf,h,hi
     &     ,anxi,anyi,anzi,lspr,iprint,.true.)
      lstrs = .true.

c      print *,'one_shot: 06'
      if( iprint.ge.10 ) print *,'get_force,myid_md,lcalcgrad='
     &     ,myid_md,lcalcgrad
      if( .not.lcalcgrad ) then
        call get_force(namax,natm,tag,ra,nnmax,aa,strs,chg,chi,stnsr
     &       ,h,hi,tcom,nb,nbmax,lsb,lsex,nex,lsrc,myparity,nn,sv,rc
     &       ,lspr,mpi_md_world,myid_md,epi,epot,nspmax,acon,lstrs
     &       ,ifcoulomb,iprint,.true.,lvc,lcell_updated,boundary)
      else  ! lcalcgrad = .true.
        epot = 0d0
        if( use_force('Morse') ) call gradw_Morse(namax,natm,tag,
     &       ra,nnmax,h,rc,lspr,epott,iprint,ndimp,gwe,gwf,gws,
     &       lematch,lfmatch,lsmatch)
c        if( use_force('vcMorse') ) call gradw_vcMorse(namax,natm,tag,ra
c     &       ,nnmax,chg,h,rc,lspr,epott,iprint,ndimp,gwe,gwf,gws)
c.....Derivative of stress should be divided by the cell volume
        gws(1:ndimp,1:6) = gws(1:ndimp,1:6) /vol
      endif

      if( iprint.ne.0 ) then
        write(6,'(a,es15.7)') ' potential energy = ',epot
      endif

c      print *,'one_shot: 07'
      call sa2stnsr(natm,strs,eki,stnsr,vol,mpi_md_world)

      call space_comp(ntot,tagtot,rtot,vtot,atot,epitot,ekitot
     &     ,stot,chgtot,chitot,natm,tag,ra,va,aa,epi,eki,strs
     &     ,chg,chi,sorg,nxyz,myid_md,mpi_md_world,tspdcmp)

c.....revert forces to the unit eV/A before going out 
      if( myid_md.eq.0 ) then
        do i=1,ntot
          is= int(tagtot(i))
          aai(1:3)= h(1:3,1,0)*atot(1,i)
     &         +h(1:3,2,0)*atot(2,i)
     &         +h(1:3,3,0)*atot(3,i)
          atot(1:3,i) = aai(1:3)
c          atot(1:3,i) = aai(1:3)/acon(is)
        enddo
      endif

      if( iprint.ne.20 .and. myid_md.eq.0 .and. lvc ) then
        print *,'Atomic charges:'
        do i=1,ntot
          print '(a,i5,es11.3,4f7.3)','   ia,epi,chg,frc=',i
     &         ,epitot(i),chgtot(i),atot(1:3,i)
        enddo
      endif
      
c.....deallocate all the arrays allocated 
      deallocate(acon,fack)
      deallocate(fekin,fa2v)
      deallocate(ra,va,aa,ra0,strs,stt,tag,lspr
     &     ,epi,eki,stp,stn,lsb,lsex,chg,chi)

      if( iprint.ge.10 ) then
        write(6,'(a)') ' one_shot done.'
        write(6,*) ''
      endif
      return
      end subroutine one_shot
c=======================================================================
      subroutine setup(nspmax,am,acon,fack,dt)
      implicit none
      include "params_unit.h"
      integer,intent(in):: nspmax
      real(8),intent(in):: dt
      real(8),intent(in):: am(nspmax)
      real(8),intent(out):: acon(nspmax),fack(nspmax)

      integer:: i
      real(8):: tmp

c      am(1:nspmax)= am(1:nspmax)

c.....umass*Ang^2/fs^2 ---> kg*m^2/sec^2 (=Joule)
      tmp= amu2kg *ang2m**2 /fs2s**2
c      write(6,'(a,es12.4,a,es12.4,a)') ' umass*Ang^2/fs^2 ='
c     &     ,tmp,' Joule =',tmp*j2ev,' eV'

c-----prefactors for normalized accelerration
      do i=1,nspmax
        acon(i)=  ev2j *0.5d0 *(dt*fs2s)**2 /(am(i)*amu2kg) *m2ang**2
c        acon(i)=  ev2j *0.5d0 *(dt*fs2s)**2 /(am(i)*amu2kg)
c        print *,'i,acon(i)=',i,acon(i)
      enddo

c-----prefactors for kinetic energy in eV
      do i=1,nspmax
        fack(i)= 0.5d0*(am(i)*amu2kg)/(dt*fs2s)**2 *ang2m**2 *j2ev
      enddo

      end subroutine setup
c=======================================================================
      subroutine setup2(nspmax,am,fekin,fa2v,dt)
      implicit none
      include "params_unit.h"
      integer,intent(in):: nspmax
      real(8),intent(in):: dt
      real(8),intent(in):: am(nspmax)
      real(8),intent(out):: fekin(nspmax),fa2v(nspmax)

      integer:: i
      real(8):: tmp

c      am(1:nspmax)= am(1:nspmax)

c.....umass*Ang^2/fs^2 ---> kg*m^2/sec^2 (=Joule)
      tmp= amu2kg *ang2m**2 /fs2s**2
c      write(6,'(a,es12.4,a,es12.4,a)') ' umass*Ang^2/fs^2 ='
c     &     ,tmp,' Joule =',tmp*j2ev,' eV'

c.....Factors for the following quantities
c.....Force: eV/Ang = (me*bohr**2/atu**2)/Ang
c.....Acceleration: Ang/fs/fs
c.....Velocity: Ang/fs
c.....Position: Ang
c.....Energy: eV = (me*bohr*bohr/atu/atu)
      do i=1,nspmax
        fa2v(i) = 0.5d0 *ume2ump *bohr2ang**2 /atu2fs**2 /am(i)
        fekin(i)=  0.5d0 *(am(i)*ump2ume) *ang2bohr**2 /fs2atu**2
c        acon(i)=  ev2j *0.5d0 *(dt*fs2s)**2 /(am(i)*amu2kg)
c        print *,'i,fekin(i),fa2v(i)=',i,fekin(i),fa2v(i)
      enddo

      end subroutine setup2
c=======================================================================
      subroutine boxmat(h,hi,ht,g,gi,gt,vol,sgm)
c-----------------------------------------------------------------------
c  setup matrices of MD-box
c    H:   MD-box matrix
c    HI:  inverse MD-box matrix
c    SGM: cofactor matrix
c-----------------------------------------------------------------------
      implicit none
      real(8),intent(in):: h(3,3,0:1)
      real(8),intent(out):: vol,sgm(3,3),hi(3,3),ht(3,3,0:1)
     &     ,g(3,3,0:1),gi(3,3),gt(3,3,0:1)

      real(8):: hit(3,3)
      integer:: i,j,k,im,ip,jm,jp

c-----cofactor matrix, SGM
      do j=1,3
        jm=mod(j+1,3)+1
        jp=mod(j,  3)+1
        do i=1,3
          im=mod(i+1,3)+1
          ip=mod(i,  3)+1
          sgm(i,j)=h(ip,jp,0)*h(im,jm,0)-h(im,jp,0)*h(ip,jm,0)
        enddo
      enddo
c-----MD-box volume
      vol=h(1,1,0)*sgm(1,1)+h(2,1,0)*sgm(2,1)+h(3,1,0)*sgm(3,1)
      do j=1,3
        do i=1,3
          hit(i,j)= sgm(i,j)/vol
        enddo
      enddo
c-----transpose
      do j=1,3
        do i=1,3
          hi(i,j)= hit(j,i)
        enddo
      enddo

c.....Set transpose
      do j=1,3
        do i=1,3
          ht(i,j,0:1)= h(j,i,0:1)
        enddo
      enddo

c.....Set G-matrix
      g(1:3,1:3,0:1)= 0d0
      do j=1,3
        do i=1,3
          do k=1,3
            g(i,j,0)=g(i,j,0) +ht(i,k,0)*h(k,j,0)
            g(i,j,1)=g(i,j,1) +ht(i,k,1)*h(k,j,0)
     &           +ht(i,k,0)*h(k,j,1)
          enddo
        enddo
      enddo
c.....Transpose of G
      do j=1,3
        do i=1,3
          gt(i,j,0:1)= g(j,i,0:1)
        enddo
      enddo
c.....Inverse of G
      call ludc_inv(3,g(1,1,0),gi)
c      write(6,'(a)') ' g:'
c      write(6,'(3es12.4)') g(1,1:3,0)
c      write(6,'(3es12.4)') g(2,1:3,0)
c      write(6,'(3es12.4)') g(3,1:3,0)
c      write(6,'(a)') ' gi:'
c      write(6,'(3es12.4)') gi(1,1:3)
c      write(6,'(3es12.4)') gi(2,1:3)
c      write(6,'(3es12.4)') gi(3,1:3)

      return
      end subroutine boxmat
c=======================================================================
      subroutine ntset(myx,myy,myz,nx,ny,nz,nn,sv,myparity,anxi,anyi
     &     ,anzi)
c-----------------------------------------------------------------------
c  Preparation for network related properties
c-----------------------------------------------------------------------
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
c-------scalar neighbor ID, nn
        nn(ku)=k1x*(ny*nz)+k1y*nz+k1z
c-------shift vector for exchnaging ra
        sv(1,ku)= anxi*iv(1,ku)
        sv(2,ku)= anyi*iv(2,ku)
        sv(3,ku)= anzi*iv(3,ku)
      enddo
      
c-----Set up the node parity table
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
c=======================================================================
      subroutine get_ekin(namax,natm,va,tag,h,nspmax,fekin,ekin,eki,ekl
     &     ,vmax,mpi_md_world)
      implicit none 
      include "mpif.h"
      integer,intent(in):: namax,natm,mpi_md_world,nspmax
      real(8),intent(in):: va(3,namax),h(3,3),fekin(nspmax)
     &     ,tag(namax)
      real(8),intent(out):: ekin,eki(3,3,namax),vmax,ekl(9)
c-----locals
      integer:: i,ierr,is,ixyz,jxyz,imax,ifmv
      real(8):: ekinl,x,y,z,v(3),v2,vmaxl,ekll(9)
      integer,external:: ifmvOf

      ekinl=0d0
      eki(1:3,1:3,1:natm)= 0d0
      ekll(1:9)= 0d0
      vmaxl= 0d0

      do i=1,natm
        is= int(tag(i))
c.....ifmv= int(mod(tag(i)*10,10d0))
        ifmv = ifmvOf(tag(i))
        if( ifmv.eq.0 ) cycle
        x= va(1,i)
        y= va(2,i)
        z= va(3,i)
        v(1:3)= h(1:3,1)*x +h(1:3,2)*y +h(1:3,3)*z
c        print *,'i,is,va,v=',i,is,x,y,z,v(1:3)
c.....Tensor form eki
        do jxyz=1,3
          do ixyz=1,3
            eki(ixyz,jxyz,i)= v(ixyz)*v(jxyz)
          enddo
        enddo
        v2= eki(1,1,i) +eki(2,2,i) +eki(3,3,i)
c        eki(1:3,1:3,i)=eki(1:3,1:3,i)*fack(is)
        eki(1:3,1:3,i)=eki(1:3,1:3,i)*fekin(is)
        ekinl=ekinl +eki(1,1,i) +eki(2,2,i) +eki(3,3,i)
c.....ekin of each ifmv
        ekll(ifmv)= ekll(ifmv) +eki(1,1,i) +eki(2,2,i) +eki(3,3,i)
c.....Find max speed
        if( v2.gt.vmaxl ) imax=i
        vmaxl=max(vmaxl,v2)
      enddo
c      print *,'imax,vmax,ekin=',imax,vmaxl,ekinl
c      print *,'va(1:3,imax)=',va(1:3,imax)

      ekin= 0d0
      call mpi_allreduce(ekinl,ekin,1,mpi_real8
     &     ,MPI_SUM,mpi_md_world,ierr)
      ekl(1:9)= 0d0
      call mpi_allreduce(ekll,ekl,9,mpi_real8
     &     ,mpi_sum,mpi_md_world,ierr)
      vmax= 0d0
      call mpi_allreduce(vmaxl,vmax,1,mpi_real8
     &     ,mpi_max,mpi_md_world,ierr)
      vmax=dsqrt(vmax)

      end subroutine get_ekin
c=======================================================================
      subroutine calc_temp_dist(iotdst,ntdst,tdst,nadst,natm,ra,eki
     &     ,istp,nouterg,myid_md,mpi_md_world,sorg)
      implicit none
      include 'mpif.h'
      include './params_unit.h'
      integer,intent(in):: iotdst,ntdst,natm,istp,nouterg,myid_md
     &     ,mpi_md_world
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
c.....before writing to file, accumurate data of other nodes
        tdl(1:ntdst)= 0d0
        ndl(1:ntdst)= 0
        call mpi_reduce(nadst,ndl,ntdst
     &       ,mpi_integer,mpi_sum,0,mpi_md_world,ierr)
        call mpi_reduce(tdst,tdl,ntdst
     &       ,mpi_double_precision,mpi_sum,0,mpi_md_world,ierr)
        do i=1,ntdst
          if( ndl(i).ne.0 ) tdl(i)= tdl(i)*2d0/3/ndl(i)/fkb
          if( myid_md.eq.0 ) then
            write(iotdst,'(f10.5,f15.3)') ((i-1)+0.5)/ntdst, tdl(i)
          endif
        enddo
        if( myid_md.eq.0 ) write(iotdst,*) ''
c.....initialize tdst
        tdst(1:ntdst)= 0d0
        nadst(1:ntdst)= 0
      endif

      end subroutine calc_temp_dist
c=======================================================================
      subroutine get_num_dof(natm,tag,fmv,ndof,myid_md,mpi_md_world
     &     ,iprint)
      implicit none
      include 'mpif.h'
      integer,intent(in):: natm,myid_md,mpi_md_world,iprint
c.....TODO: hard coding number 9 should be replaced...
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
      call mpi_allreduce(ndofl,ndof,9,mpi_integer,mpi_sum
     &     ,mpi_md_world,ierr)
      if(myid_md.eq.0 .and. iprint.ne.0 )
     &     write(6,'(/,a,9(2x,i0))')
     &     ' Degrees of freedom for each ifmv =',ndof(1:9)
      return
      end subroutine get_num_dof
c=======================================================================
      subroutine check_size_and_parallel(sgm,vol,rc,anxi,anyi,anzi
     &     ,nx,ny,nz,myid_md)
      implicit none
      include 'mpif.h'
      integer,intent(in):: nx,ny,nz,myid_md
      real(8),intent(in):: sgm(3,3),vol,rc,anxi,anyi,anzi

      integer:: ierr
      real(8):: vala,valb,valc,rca,rcb,rcc
      
c-----calculate the cut-off lengths
      vala=dsqrt(sgm(1,1)**2+sgm(2,1)**2+sgm(3,1)**2)/vol
      valb=dsqrt(sgm(1,2)**2+sgm(2,2)**2+sgm(3,2)**2)/vol
      valc=dsqrt(sgm(1,3)**2+sgm(2,3)**2+sgm(3,3)**2)/vol
      rca=rc*vala
      rcb=rc*valb
      rcc=rc*valc

c.....check whether the size of the system is large enough
      if(  (rca.gt.anxi .and. nx.gt.1) .or.
     &     (rcb.gt.anyi .and. ny.gt.1) .or.
     &     (rcc.gt.anzi .and. nz.gt.1) ) then
        if( myid_md.le.0 ) then
          write(6,'(a)') " Error: about size and parallelization"
          write(6,'(a,2f15.3,i5)') "rca,anxi,nx = ",rca,anxi,nx
          write(6,'(a,2f15.3,i5)') "rcb,anyi,ny = ",rcb,anyi,ny
          write(6,'(a,2f15.3,i5)') "rcc,anzi,nz = ",rcc,anzi,nz
          write(6,'(a)') " * All the size (xyz) per node "
     &         //"with more than 1 node should be larger "
     &         //"than the cutoff radius."
          write(6,*) ""
        endif
        call mpi_finalize(ierr)
        stop
      endif

      end subroutine check_size_and_parallel
c=======================================================================
      subroutine bacopy(tcom,sgm,vol,lsb,lsex,nbmax,ra,namax
     &     ,natm,nb,anxi,anyi,anzi,nn,tag,rc,myid,myparity
     &     ,lsrc,sv,nex,mpi_md_world,iprint,ifcoulomb,chg,chi,l1st
     &     ,boundary)
c-----------------------------------------------------------------------
c  Exchanges boundary-atom data among neighbor nodes: tag and ra
c  Normal data passing and repeating the cell are mixed, and
c  automatically selected according to the size of the system.
c  
c  Note: parallelized to smaller than rcut should not happen.
c-----------------------------------------------------------------------
      use force
      use pmdmpi,only: nx,ny,nz,nid2xyz
      implicit none
      include 'mpif.h'
      integer,intent(in):: namax,nbmax,nn(6)
     &     ,myid,myparity(3),mpi_md_world,iprint
      integer,intent(inout):: natm,nb,lsb(0:nbmax,6),lsex(nbmax,6)
     &     ,lsrc(6),nex(3)
      real(8),intent(in):: sv(3,6),sgm(3,3),vol,anxi,anyi,anzi,rc
      real(8),intent(inout):: tag(namax),ra(3,namax),tcom
      integer,intent(in):: ifcoulomb
      real(8),intent(inout):: chg(namax),chi(namax)
      logical,intent(in):: l1st
      character(len=3):: boundary

c.....integer:: status(MPI_STATUS_SIZE)
      integer:: i,j,kd,kdd,kul,kuh,ku,ierr,iex,ix,iy,iz
      integer:: nav,maxna,maxb,inode,nsd,nrc,nbnew
      real(8):: tcom1,tcom2,xi(3),rcv(3),asgm
      logical,external:: bbd
      real(8),save,allocatable:: dbuf(:,:),dbufr(:,:)
      logical:: lshort(3)

      if( l1st ) then
        if( allocated(dbuf) ) deallocate(dbuf,dbufr)
        if( use_charge() ) then
          allocate(dbuf(6,nbmax),dbufr(6,nbmax))
        else
          allocate(dbuf(4,nbmax),dbufr(4,nbmax))
        endif
      endif

      if( use_charge() ) then
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
      
c-----reset the num of "received" boundary atoms
      nbnew=0

c-----calculate the cut-off lengths
      do kd=1,3
        asgm= dsqrt(sgm(1,kd)**2 +sgm(2,kd)**2 +sgm(3,kd)**2)
        rcv(kd)= rc*asgm/vol
        nex(kd)= int(rcv(kd)) +1
      enddo
      if( l1st .and. myid.eq.0 .and. iprint.ge.2 ) then
        write(6,'(a,3f10.3)') ' rcv = ',rcv(1:3)
        write(6,'(a,3i10)')   ' nex = ',nex(1:3)
      endif

c-----loop over x, y, & z directions
      do 100 kd=1,3

c-------No. of to-be-copied atoms, LSB(0,)
        do kdd= -1,0
          lsb(0,2*kd+kdd)= 0
        enddo

c.....Dependent on the unit vector size, different treatment is required
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
c-------Scan all the residents & copies
          do i=1,natm+nbnew
            xi(1:3)= ra(1:3,i)
c---------For low & high directions
            kul=2*kd-1
            kuh=2*kd
c---------Lower neighbor
            if (bbd(xi(1),xi(2),xi(3),rcv(1),rcv(2),rcv(3)
     &           ,kul,anxi,anyi,anzi)) then
              lsb(0,kul)=lsb(0,kul)+1
              lsb(lsb(0,kul),kul)=i
            endif
c---------Higher neighbor
            if(bbd(xi(1),xi(2),xi(3),rcv(1),rcv(2),rcv(3)
     &           ,kuh,anxi,anyi,anzi)) then
              lsb(0,kuh)=lsb(0,kuh)+1
              lsb(lsb(0,kuh),kuh)=i
            endif
          enddo
        endif                   ! if(nex.gt.1)

c.....Depending on boundary condition and cell position,
c.....number of to-be-sent atoms is set 0.
        if( boundary(kd:kd).ne.'p' )  then
          kul = 2*kd -1
          kuh = 2*kd
          if( kd.eq.1 ) then
            if( ix.eq.0 ) then
              lsb(0,kul) = 0
            else if (ix.eq.nx-1 ) then
              lsb(0,kuh) = 0
            endif
          else if( kd.eq.2 ) then
            if( iy.eq.0 ) then
              lsb(0,kul) = 0
            else if (iy.eq.ny-1 ) then
              lsb(0,kuh) = 0
            endif
          else if( kd.eq.3 ) then
            if( iz.eq.0 ) then
              lsb(0,kul) = 0
            else if ( iz.eq.nz-1 ) then
              lsb(0,kuh) = 0
            endif
          endif
        endif

c-------Error trap
        do kdd= -1,0
          ku=2*kd+kdd
          call MPI_ALLREDUCE(lsb(0,ku),maxb,1,MPI_INTEGER,MPI_MAX,
     &         mpi_md_world,ierr)
          if (maxb.gt.nbmax) then
            if (myid.eq.0) then
              write(*,*)'Buffer or list overflowed at bacopy'
              write(*,*)'LSB(0) NBMAX =',maxb,nbmax
            endif
            call MPI_FINALIZE(ierr)
            stop
          endif
        enddo

        if( nex(kd).gt.1 ) then
          do kdd= -1,0
            ku= 2*kd+kdd
            if( .not. use_charge() ) then
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
c-------To calculate the communication time
          tcom1=MPI_WTIME()
          do 200 kdd= -1,0
            ku=2*kd+kdd
            inode=nn(ku)
            nsd=lsb(0,ku)
            call mespasi(inode,myparity(kd),nsd,nrc,1,1,10
     &           ,mpi_md_world)
c---------Store the # of received boundary atoms in LSRC
            lsrc(ku)=nrc

c---------Exchange ra and tag
            if( .not. use_charge() ) then
              do i=1,nsd
                j= lsb(i,ku)
                dbuf(1:3,i)= ra(1:3,j) -sv(1:3,ku)
                dbuf(4,i)  = tag(j)
              enddo
              call mespasd(inode,myparity(kd),dbuf,dbufr,nsd*4,nrc*4,21
     &             ,mpi_md_world)
              do i=1,nrc
                ra(1:3,natm+nbnew+i)= dbufr(1:3,i)
                tag(natm+nbnew+i)   = dbufr(4,i)
              enddo
            else
              do i=1,nsd
                j= lsb(i,ku)
                dbuf(1:3,i)= ra(1:3,j) -sv(1:3,ku)
                dbuf(4,i)  = tag(j)
                dbuf(5,i)  = chg(j)
                dbuf(6,i)  = chi(j)
              enddo
              call mespasd(inode,myparity(kd),dbuf,dbufr,nsd*6,nrc*6,21
     &             ,mpi_md_world)
              do i=1,nrc
                ra(1:3,natm+nbnew+i)= dbufr(1:3,i)
                tag(natm+nbnew+i)   = dbufr(4,i)
                chg(natm+nbnew+i)   = dbufr(5,i)
                chi(natm+nbnew+i)   = dbufr(6,i)
              enddo
            endif

            call MPI_BARRIER(mpi_md_world,ierr)
c---------increase the # of received boundary atoms
            nbnew=nbnew+nrc

 200      continue
c-------Add the communication time to COMT
          tcom2=MPI_WTIME()
          tcom=tcom+tcom2-tcom1
        endif

c-------Error trap
        nav=natm+nbnew
        call MPI_ALLREDUCE(nav,maxna,1,MPI_INTEGER,MPI_MAX
     &       ,mpi_md_world,ierr)
        if (maxna.gt.namax) then
          if (myid.eq.0) then
            write(*,*)'NAMAX overflowed at bacopy'
            write(*,*)'N+NB NAMAX = ',maxna,namax
          endif
          call MPI_FINALIZE(ierr)
          stop
        endif

 100  continue

c-----num. of received boundary atoms
      nb=nbnew

      end subroutine bacopy
c=======================================================================
      subroutine bacopy_fixed(tcom,sgm,vol,lsb,lsex,nbmax,ra,namax
     &     ,natm,nb,anxi,anyi,anzi,nn,tag,rc,myid_md,myparity,lsrc,sv
     &     ,nex
     &     ,mpi_md_world,ifcoulomb,chg,chi,boundary)
c-----------------------------------------------------------------------
c  Exchanges boundary-atom data among neighbor nodes: tag and ra
c  This doesnt search using position, just send & recv data of atoms
c    which were listed by 'bacopy'.
c  Different number of data are copied depending on whether 
c    using atomic charges or not.
c-----------------------------------------------------------------------
      use force
      use pmdmpi,only: nx,ny,nz,nid2xyz
      implicit none
      include 'mpif.h'
      integer,intent(in):: namax,nbmax,nn(6),natm,nb
     &     ,myid_md,myparity(3),mpi_md_world,nex(3)
      integer,intent(inout):: lsb(0:nbmax,6),lsex(nbmax,6),lsrc(6)
      real(8),intent(in):: sv(3,6),sgm(3,3),vol,anxi,anyi,anzi,rc
      real(8),intent(inout):: tag(namax),ra(3,namax),tcom
      integer,intent(in):: ifcoulomb
      real(8),intent(inout):: chg(namax),chi(namax)
      character(len=3),intent(in):: boundary

c      integer:: status(MPI_STATUS_SIZE)
      integer:: i,j,kd,kdd,ku,ierr,iex,ix,iy,iz
      integer:: inode,nsd,nrc,nbnew
      real(8):: tcom1,tcom2 !,rcv(3),asgm
      real(8),save,allocatable:: dbuf(:,:),dbufr(:,:)
      logical,save:: l1st=.true.

      if( l1st ) then
        if( allocated(dbuf) ) deallocate(dbuf,dbufr)
        if( .not. use_charge() ) then
          allocate(dbuf(4,nbmax),dbufr(4,nbmax))
        else
          allocate(dbuf(6,nbmax),dbufr(6,nbmax))
        endif
        l1st=.false.
      endif

      if( .not. use_charge() ) then
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

cc-----calculate the cut-off lengths
c      do kd=1,3
c        asgm= dsqrt(sgm(1,i)**2 +sgm(2,i)**2 +sgm(3,i)**2)
c        rcv(kd)= rc*asgm/vol
c        nex(kd)= int(rcv(kd)) +1
c      enddo

c-----loop over x, y, & z directions
      do 100 kd=1,3

c-------To calculate the communication time
        tcom1=MPI_WTIME()

        if( nex(kd).gt.1 ) then
          do kdd= -1,0
            ku= 2*kd+kdd
            if( .not. use_charge() ) then
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
          do 200 kdd= -1,0
            ku=2*kd+kdd
            inode=nn(ku)
            nsd=lsb(0,ku)
            nrc=lsrc(ku)

c---------Exchange ra and tag
            if( .not. use_charge() ) then  ! in case of no charge
              do i=1,nsd
                j= lsb(i,ku)
                dbuf(1:3,i)= ra(1:3,j) -sv(1:3,ku)
                dbuf(4,i)  = tag(j)
              enddo
              call mespasd(inode,myparity(kd),dbuf,dbufr,nsd*4,nrc*4,21
     &             ,mpi_md_world)
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
              call mespasd(inode,myparity(kd),dbuf,dbufr,nsd*6,nrc*6,21
     &             ,mpi_md_world)
              do i=1,nrc
                ra(1:3,natm+nbnew+i)= dbufr(1:3,i)
                tag(natm+nbnew+i)   = dbufr(4,i)
                chg(natm+nbnew+i)   = dbufr(5,i)
                chi(natm+nbnew+i)   = dbufr(6,i)
              enddo
            endif

            call MPI_BARRIER(mpi_md_world,ierr)
            nbnew=nbnew +nrc
 200      continue
        endif

c-------Add the communication time to COMT
        tcom2=MPI_WTIME()
        tcom=tcom+tcom2-tcom1

100   continue

      end subroutine bacopy_fixed
c=======================================================================
      subroutine bacopy_chg_fixed(tcom,lsb,lsex,nbmax,namax
     &     ,natm,nb,nn,myid_md,myparity,lsrc
     &     ,nex,mpi_md_world,chg,boundary)
c-----------------------------------------------------------------------
c  Exchange only chg of boundary-atom
c  This doesnt search using position, just send & recv data of atoms
c    which were listed by 'bacopy'.
c-----------------------------------------------------------------------
      use pmdmpi,only: nx,ny,nz,nid2xyz
      implicit none
      include 'mpif.h'
      integer,intent(in):: namax,nbmax,nn(6),natm,nb
     &     ,myid_md,myparity(3),mpi_md_world,nex(3)
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

cc-----calculate the cut-off lengths
c      do kd=1,3
c        asgm= dsqrt(sgm(1,i)**2 +sgm(2,i)**2 +sgm(3,i)**2)
c        rcv(kd)= rc*asgm/vol
c        nex(kd)= int(rcv(kd)) +1
c      enddo

c-----loop over x, y, & z directions
      do 100 kd=1,3

c-------To calculate the communication time
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
          do 200 kdd= -1,0
            ku=2*kd+kdd
            inode=nn(ku)
            nsd=lsb(0,ku)
            nrc=lsrc(ku)

c---------Exchange ra and tag
            do i=1,nsd
              j= lsb(i,ku)
              dbuf(i)  = chg(j)
            enddo
            call mespasd(inode,myparity(kd),dbuf,dbufr,nsd,nrc,23
     &           ,mpi_md_world)
            do i=1,nrc
              chg(natm+nbnew+i)   = dbufr(i)
            enddo

            call mpi_barrier(mpi_md_world,ierr)
            nbnew=nbnew +nrc
 200      continue
        endif

c-------Add the communication time to COMT
        tcom2=MPI_WTIME()
        tcom=tcom+tcom2-tcom1

100   continue

      end subroutine bacopy_chg_fixed
c=======================================================================
      subroutine bamove(tcom,namax,nbmax,natm,ra,va,tag,chg,chi
     &     ,anxi,anyi,anzi,myid_md,nn,sv,myparity,mpi_md_world
     &     ,boundary)
c-----------------------------------------------------------------------
c  Exchange atoms between neighbor nodes and atomic data as well.
c
c  MVQUE(0:NBMAX,6):
c    MVQUE(0,ku) is the # of to-be-moved atoms to neighbor ku;
c    MVQUE(i,ku) is the adress, in IS, of atom i
c----------------------------------------------------------------------
      use pmdmpi,only: nid2xyz,nx,ny,nz
      implicit none
      include 'mpif.h'
      integer,intent(in):: namax,nbmax
      integer,intent(in):: myid_md,mpi_md_world,myparity(3),nn(6)
      real(8),intent(in):: anxi,anyi,anzi,sv(3,6)
      integer,intent(inout):: natm
      real(8),intent(inout):: tcom,ra(3,namax),va(3,namax),tag(namax)
     &     ,chg(namax),chi(namax)
      character(len=3),intent(in):: boundary

c      integer:: status(MPI_STATUS_SIZE)
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

      call nid2xyz(myid_md,ix,iy,iz)

c-----newim: num. of new immigrants
      newim= 0
      mvque(0,1:6)= 0

      do 100 kd=1,3
        kul=2*kd-1
        kuh=2*kd

c-------num of to-be-moved atoms
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

c.....Depending on boundary condition and cell position,
c.....number of to-be-moved atoms is set 0.
        if( boundary(kd:kd).ne.'p' )  then
          if( kd.eq.1 ) then
            if( ix.eq.0 ) then
              mvque(0,kul) = 0
            else if (ix.eq.nx-1 ) then
              mvque(0,kuh) = 0
            endif
          else if( kd.eq.2 ) then
            if( iy.eq.0 ) then
              mvque(0,kul) = 0
            else if (iy.eq.ny-1 ) then
              mvque(0,kuh) = 0
            endif
          else if( kd.eq.3 ) then
            if( iz.eq.0 ) then
              mvque(0,kul) = 0
            else if ( iz.eq.nz-1 ) then
              mvque(0,kuh) = 0
            endif
          endif
        endif
c
c        if( kd.eq.1 ) then
c          if( ix.eq.0 .and. mvque(0,kul).gt.0 ) then
c            print *,' myid,kd,ix,mvque(0,kul)=',myid_md
c     &           ,kd,ix,mvque(0,kul)
c          else if( ix.eq.nx-1 .and. mvque(0,kuh).gt.0 ) then
c            print *,' myid,kd,ix,mvque(0,kuh)=',myid_md
c     &           ,kd,ix,mvque(0,kuh)
c          endif
c        endif

c-------Error trap
        if (mvque(0,kul).gt.nbmax) then
          print *,'Buffer overflowed at bamove node',myid_md
          print *,'# in MVQUE=',mvque(0,kul)
          stop
        else if (mvque(0,kuh).gt.nbmax) then
          print *,'Buffer overflowed at bamove node',myid_md
          print *,'# in MVQUE=',mvque(0,kuh)
          stop
        endif

c-------To calculate the communacation time
        tcom1=MPI_WTIME()

        do 200 kdd= -1,0

          ku=2*kd+kdd
          inode=nn(ku)
          nsd=mvque(0,ku)
          
          call mespasi(inode,myparity(kd),nsd,nrc,1,1,60,
     &         mpi_md_world)

c---------move atom positions, RA
          do i=1,nsd
            j= mvque(i,ku)
            dbuf(1:3,i)= ra(1:3,j) -sv(1:3,ku)
            dbuf(4:6,i)= va(1:3,j)
            dbuf(7,i)  = tag(j)
            dbuf(8,i)  = chg(j)
            dbuf(9,i)  = chi(j)
c-----------Eliminate the record of a moved-out atom
            tag(j)= 0d0
            chg(j)= 0d0
            chi(j)= 0d0
          enddo
          call mespasd(inode,myparity(kd),dbuf,dbufr,9*nsd,9*nrc,71
     &         ,mpi_md_world)
          do i=1,nrc
            ra(1:3,natm+newim+i)= dbufr(1:3,i)
            va(1:3,natm+newim+i)= dbufr(4:6,i)
            tag(natm+newim+i)   = dbufr(7,i)
            chg(natm+newim+i)   = dbufr(8,i)
            chi(natm+newim+i)   = dbufr(9,i)
          enddo

          newim=newim+nrc
          call MPI_BARRIER(mpi_md_world,ierr)

200     continue

        tcom2=MPI_WTIME()
        tcom=tcom+tcom2-tcom1 

100   continue

c-----Compression
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
c-----Update # of resident atoms
      natm=ipt

      return
      end subroutine bamove
c=======================================================================
      function bbd(xv,yv,zv,rcav,rcbv,rccv,ku,anxi,anyi,anzi)
c-----------------------------------------------------------------------
c  BBD = .true. if the coordinates are in the boundary to neighbor ku
c-----------------------------------------------------------------------
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
c=======================================================================
      function bmv(xv,yv,zv,ku,anxi,anyi,anzi)
c-----------------------------------------------------------------------
c  BMV = .true. if the coordinates should belong to neighbor ku
c------------------------------------------------------------------------
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
c=======================================================================
      subroutine force_isobaric(stgt,ptgt,ah,natm,eki,strs,sgm
     &     ,dt,srlx,stbeta,vol,stnsr,mpi_md_world,cpctl)
c
c  Calc. acceralation of h-matrix used in isobaric MD
c  using given stress tensor stgt.
c
      implicit none
      include "mpif.h"
      integer,intent(in):: natm,mpi_md_world
      real(8),intent(in):: stgt(3,3),eki(3,3,natm),strs(3,3,natm),ptgt
     &     ,sgm(3,3),dt,vol,srlx,stbeta
      real(8),intent(out):: ah(3,3),stnsr(3,3)
      character(len=*):: cpctl
      
      integer:: i,ixyz,jxyz,ierr
      real(8):: prss,fac

      call sa2stnsr(natm,strs,eki,stnsr,vol,mpi_md_world)

c      print '(a,9f8.3)','stnsr = ',stnsr(1:3,1:3)

c.....Berendsen for variable-cell
      if( trim(cpctl).eq.'Berendsen' .or.
     &     trim(cpctl).eq.'vc-Berendsen' ) then
c.....now ah is scaling factor for h-mat
        ah(1:3,1:3)= 0d0
        ah(1,1)= 1d0
        ah(2,2)= 1d0
        ah(3,3)= 1d0
        ah(1:3,1:3)= ah(1:3,1:3)
     &       -stbeta*dt/3/srlx*( stgt(1:3,1:3)-stnsr(1:3,1:3) )
c.....Berendsen for variable-volume not variable-cell
      else if( trim(cpctl).eq.'vv-Berendsen' ) then
        ah(1:3,1:3)= 0d0
        ah(1,1)= 1d0
        ah(2,2)= 1d0
        ah(3,3)= 1d0
        prss = (stnsr(1,1)+stnsr(2,2)+stnsr(3,3))/3
        fac = 1.0 -stbeta/dt/3/srlx *(ptgt-prss)
        ah(1:3,1:3) = ah(1:3,1:3)*fac
      endif

      end subroutine force_isobaric
c=======================================================================
      subroutine sa2stnsr(natm,strs,eki,stnsr,vol,mpi_md_world)
c      
c  Take sum of atomic stresses to compute cell stress
c
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
      stg(1:3,1:3)= 0d0
      call mpi_allreduce(stl,stg,9,mpi_real8,mpi_sum
     &     ,mpi_md_world,ierr)
      stnsr(1:3,1:3) = stnsr(1:3,1:3) +stg(1:3,1:3)/vol
      return
      end subroutine sa2stnsr
c=======================================================================
      subroutine setv(h,natm,tag,va,nspmax,am,tinit,dt)
      implicit none
      include 'mpif.h'
      include 'params_unit.h'
      integer,intent(in):: natm,nspmax
      real(8),intent(in):: tag(natm),am(nspmax),tinit,dt,h(3,3,0:1)
      real(8),intent(out):: va(3,natm)
      integer:: i,l,is
      real(8):: dseed,sumvx,sumvy,sumvz,tmp,facv(nspmax)
      real(8),parameter:: pi = 3.14159265358979d0
      real(8),external:: box_muller

c      facv(1:nspmax)=dsqrt(tinit*fkb*ev2j /(am(1:nspmax)*amu2kg))
c     &     *m2ang /s2fs
      facv(1:nspmax) = dsqrt(tinit*fkb
     &     *(ume2ump* bohr2ang**2/atu2fs**2) /am(1:nspmax))
      
c-----velocities in Maxwell-Boltzmann distribution
      do i=1,natm
        is= int(tag(i))
        do l=1,3
          va(l,i)=facv(is) *box_muller()
        enddo
c        write(6,'(i3,i5,i3,es24.14,3es15.7)') myid_md,i,is,tag(i)
c     &       ,va(1:3,i)
      enddo

      do i=1,natm
        va(1,i)= va(1,i) /h(1,1,0) !*dt
        va(2,i)= va(2,i) /h(2,2,0) !*dt
        va(3,i)= va(3,i) /h(3,3,0) !*dt
      enddo

      end subroutine setv
c=======================================================================
      subroutine rm_trans_motion(natm,tag,va,nspmax,am
     &     ,mpi_md_world,myid_md,iprint)
      implicit none
      include 'mpif.h'
      integer,intent(in):: natm,nspmax,mpi_md_world,myid_md,iprint
      real(8),intent(in):: tag(natm),am(nspmax)
      real(8),intent(out):: va(3,natm)
      
      integer:: i,is,ierr
      real(8):: sumpx,sumpy,sumpz,amss,amtot,tmp

c-----set center of mass motion to zero
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
      call mpi_allreduce(tmp,sumpx,1,mpi_double_precision,mpi_sum
     &     ,mpi_md_world,ierr)
      tmp= sumpy
      call mpi_allreduce(tmp,sumpy,1,mpi_double_precision,mpi_sum
     &     ,mpi_md_world,ierr)
      tmp= sumpz
      call mpi_allreduce(tmp,sumpz,1,mpi_double_precision,mpi_sum
     &     ,mpi_md_world,ierr)
      tmp= amtot
      call mpi_allreduce(tmp,amtot,1,mpi_double_precision,mpi_sum
     &     ,mpi_md_world,ierr)
      do i=1,natm
        va(1,i)=va(1,i)-sumpx/amtot
        va(2,i)=va(2,i)-sumpy/amtot
        va(3,i)=va(3,i)-sumpz/amtot
      enddo

c      if( myid_md.eq.0 .and. iprint.ne.0 ) then
c        write(6,'(a,3es12.4)') ' sumpx,y,z/amtot='
c     &       ,sumpx/amtot,sumpy/amtot,sumpz/amtot
c      endif

      end subroutine rm_trans_motion
c=======================================================================
      subroutine vfire(num_fire,alp0_fire,alp_fire,falp_fire,dtmax_fire
     &     ,finc_fire,fdec_fire,nmin_fire
     &     ,natm,va,aa,myid_md,mpi_md_world,dt,iprint)
      implicit none
      include 'mpif.h'
      integer,intent(in):: natm,myid_md,mpi_md_world,nmin_fire
     &     ,iprint
      real(8),intent(in):: aa(3,natm),falp_fire,alp0_fire,dtmax_fire
     &     ,finc_fire,fdec_fire
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
      call mpi_allreduce(vnorml,vnorm,1,mpi_real8
     &     ,mpi_sum,mpi_md_world,ierr)
      call mpi_allreduce(fnorml,fnorm,1,mpi_real8
     &     ,mpi_sum,mpi_md_world,ierr)
      call mpi_allreduce(fdotvl,fdotv,1,mpi_real8
     &     ,mpi_sum,mpi_md_world,ierr)
      vnorm = dsqrt(vnorm)
      fnorm = dsqrt(fnorm)
      do i=1,natm
        do ixyz=1,3
          va(ixyz,i) = (1d0 -alp_fire)*va(ixyz,i)
     &         +alp_fire*aa(ixyz,i)/fnorm *vnorm
        enddo
      enddo
c      print '(a,2e14.5)','vnorm,fnorm =',vnorm,fnorm
      
      if( fdotv.gt.0d0  ) then
        num_fire = num_fire + 1
        if( num_fire.gt.nmin_fire ) then
          dt = min(dtmax_fire,dt*finc_fire)
          alp_fire = alp_fire *falp_fire
          if( iprint.ne.0 .and. myid_md.eq.0 ) then
            write(6,'(a,f10.3,f10.5,es12.4)') ' dt,alp_fire,fdotv = '
     &           ,dt,alp_fire,fdotv
          endif
        endif
      else
        dt = min(dtmax_fire,dt*fdec_fire)
        va(1:3,1:natm) = 0d0
        alp_fire = alp0_fire
        num_fire = 0
        if( iprint.ne.0 .and. myid_md.eq.0 ) then
          write(6,'(a,f10.3,f10.5,es12.4)') ' dt,alp_fire,fdotv = '
     &         ,dt,alp_fire,fdotv
        endif
      endif

c      call mpi_bcast(dt,1,mpi_real8,0,mpi_md_world,ierr)
c      call mpi_bcast(alp_fire,1,mpi_real8,0,mpi_md_world
c     &     ,ierr)
      end subroutine vfire
c=======================================================================
      subroutine space_decomp(hunit,h,ntot,tagtot,rtot,vtot
     &     ,chgtot,chitot,myid_md,mpi_md_world,nx,ny,nz,nxyz,iprint)
c
c  Decompose the system and scatter atoms to every process.
c
      use pmdvars
      implicit none
      include 'mpif.h'
      integer,intent(in):: ntot,myid_md,mpi_md_world,nx,ny,nz,nxyz
     &     ,iprint
      real(8),intent(in):: vtot(3,ntot),chgtot(ntot)
     &     ,chitot(ntot)
      real(8),intent(in):: hunit,h(3,3,0:1)
      real(8),intent(inout):: rtot(3,ntot),tagtot(ntot)

      integer:: istat(mpi_status_size)
      integer:: i,j,ixyz,n,ierr
      integer:: myxt,myyt,myzt,nmax,nmin
      real(8):: sxogt,syogt,szogt
      real(8):: t0

      t0 = mpi_wtime()

c.....at 1st call
      if( .not. allocated(ra) ) then
        if( myid_md.eq.0 ) then
c.....wrap atom position inside [0:1)
          do i =1,ntot
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

c.....count max number of atoms in a node
          nmax = 0
          nmin = 1000000000
          do ixyz=0,nxyz-1
            myxt = ixyz/(ny*nz)
            myyt = mod(ixyz/nz,ny)
            myzt = mod(ixyz,nz)
            sxogt = dble(myxt)/nx
            syogt = dble(myyt)/ny
            szogt = dble(myzt)/nz
            n = 0
            do i=1,ntot
              if(  rtot(1,i).ge.sxogt .and.
     &             rtot(1,i).lt.sxogt+1d0/nx .and.
     &             rtot(2,i).ge.syogt .and.
     &             rtot(2,i).lt.syogt+1d0/ny .and.
     &             rtot(3,i).ge.szogt .and.
     &             rtot(3,i).lt.szogt+1d0/nz .and.
     &             tagtot(i).gt.0d0) then
                n = n+1
c.....Set the tag negative if the atom is assigned to a certain cell
                tagtot(i) = -tagtot(i)
              endif
            enddo
            nmax = max(nmax,n)
            nmin = min(nmin,n)
            print '(a,4i6,3es12.4,3i6)'
     &           ,'ixyz,myxt,myyt,myzt,sxogt,syogt,szogt,n,nmax,nmin='
     &           ,ixyz,myxt,myyt,myzt,sxogt,syogt,szogt,n,nmax,nmin
          enddo
          namax = max(int(nmax*1.2),200)
          nbmax = max(namax*27,nbmax)
          namax = namax +nbmax
          if( iprint.ne.0 ) then
            write(6,'(a,i10)') ' Min number of local atoms = ',nmin
            write(6,'(a,i10)') ' Max number of local atoms = ',nmax
            write(6,'(a,i10)')   '   nbmax = max(nmax*1.2*27,nbmax) = '
     &           ,nbmax
            write(6,'(a,i10)')   '   namax = nmax*1.2 + nbmax  = ',namax
            write(6,'(a,i10,a)') ' Memory for main routine   = '
     &           ,namax*(3*4 +3*3*5 +1 +nnmax/2)*8/1000/1000,' MByte'
          endif
c.....Reset the tags positive
          do i=1,ntot
            tagtot(i) = abs(tagtot(i))
          enddo
        endif ! myid.eq.0
        call mpi_bcast(namax,1,mpi_integer,0,mpi_md_world,ierr)
        call mpi_bcast(nbmax,1,mpi_integer,0,mpi_md_world,ierr)
c.....allocate memory for local arrays of length namax
        allocate(ra(3,namax),va(3,namax),aa(3,namax),ra0(3,namax)
     &       ,strs(3,3,namax),stt(3,3,namax),tag(namax)
     &       ,lspr(0:nnmax,namax),epi(namax),eki(3,3,namax)
     &       ,stp(3,3,namax),stn(3,3,namax)
     &       ,chg(namax),chi(namax)
     &       ,lsb(0:nbmax,6),lsex(nbmax,6))
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
          do i=1,ntot
            if(  rtot(1,i).ge.sxogt .and.
     &           rtot(1,i).lt.sxogt+1d0/nx .and.
     &           rtot(2,i).ge.syogt .and.
     &           rtot(2,i).lt.syogt+1d0/ny .and.
     &           rtot(3,i).ge.szogt .and.
     &           rtot(3,i).lt.szogt+1d0/nz .and.
     &           tagtot(i).gt.0d0 ) then
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
            call mpi_send(natm,1,mpi_integer,ixyz,ixyz+1
     &           ,mpi_md_world,ierr)
            call mpi_send(tag,natm,mpi_real8,ixyz,ixyz+2
     &           ,mpi_md_world,ierr)
            call mpi_send(ra,3*natm,mpi_real8,ixyz,ixyz+3
     &           ,mpi_md_world,ierr)
            call mpi_send(va,3*natm,mpi_real8,ixyz,ixyz+4
     &           ,mpi_md_world,ierr)
            call mpi_send(chg,natm,mpi_real8,ixyz,ixyz+5
     &           ,mpi_md_world,ierr)
            call mpi_send(chi,natm,mpi_real8,ixyz,ixyz+6
     &           ,mpi_md_world,ierr)
          endif
        enddo
c        write(6,'(a,f10.3)') ' time space_decomp = ',mpi_wtime() -t0
        tspdcmp = mpi_wtime() -t0
c.....Reset the tags positive
        do i=1,ntot
          tagtot(i) = abs(tagtot(i))
        enddo
      else ! myid_md.ne.0
        call mpi_recv(natm,1,mpi_integer,0,myid_md+1
     &       ,mpi_md_world,istat,ierr)
        call mpi_recv(tag,natm,mpi_real8,0,myid_md+2
     &       ,mpi_md_world,istat,ierr)
        call mpi_recv(ra,3*natm,mpi_real8,0,myid_md+3
     &       ,mpi_md_world,istat,ierr)
        call mpi_recv(va,3*natm,mpi_real8,0,myid_md+4
     &       ,mpi_md_world,istat,ierr)
        call mpi_recv(chg,natm,mpi_real8,0,myid_md+5
     &       ,mpi_md_world,istat,ierr)
        call mpi_recv(chi,natm,mpi_real8,0,myid_md+6
     &       ,mpi_md_world,istat,ierr)
      endif
      call mpi_barrier(mpi_md_world,ierr)

      end subroutine space_decomp
c=======================================================================
      subroutine space_comp(ntot,tagtot,rtot,vtot,atot,epitot,ekitot
     &     ,stot,chgtot,chitot,natm,tag,ra,va,aa,epi,eki,strs,chg,chi
     &     ,sorg,nxyz,myid_md,mpi_md_world,tspdcmp)
c
c  Opposite to space_decomp, gather atoms from every process
c  to create the total system for output.
c
      implicit none
      include 'mpif.h'
      integer,intent(in):: ntot,natm,nxyz,myid_md,mpi_md_world
      real(8),intent(in):: va(3,natm),aa(3,natm),epi(natm),eki(3,3,natm)
     &     ,strs(3,3,natm),tag(natm),sorg(3),chg(natm),chi(natm)
      real(8),intent(inout):: ra(3,natm),tspdcmp
      real(8),intent(out):: tagtot(ntot),rtot(3,ntot),vtot(3,ntot)
     &     ,atot(3,ntot),epitot(ntot),ekitot(3,3,ntot),stot(3,3,ntot)
     &     ,chgtot(ntot),chitot(ntot)
      integer,parameter:: nmpi = 8
      integer:: n0,ixyz,natmt,i,ierr
      integer:: istat(mpi_status_size),itag
      real(8):: t0

      t0 = mpi_wtime()

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
        n0 = n0 +1
        do ixyz=1,nxyz-1
          itag = ixyz*nmpi -nmpi
          call mpi_recv(natmt,1,mpi_integer,ixyz,itag
     &         ,mpi_md_world,istat,ierr)
          call mpi_recv(tagtot(n0),natmt,mpi_real8
     &         ,ixyz,itag+1,mpi_md_world,istat,ierr)
          call mpi_recv(rtot(1,n0),3*natmt,mpi_real8
     &         ,ixyz,itag+2,mpi_md_world,istat,ierr)
          call mpi_recv(vtot(1,n0),3*natmt,mpi_real8
     &         ,ixyz,itag+3,mpi_md_world,istat,ierr)
          call mpi_recv(epitot(n0),natmt,mpi_real8
     &         ,ixyz,itag+4,mpi_md_world,istat,ierr)
          call mpi_recv(ekitot(1,1,n0),3*3*natmt,mpi_real8
     &         ,ixyz,itag+5,mpi_md_world,istat,ierr)
          call mpi_recv(stot(1,1,n0),3*3*natmt,mpi_real8
     &         ,ixyz,itag+6,mpi_md_world,istat,ierr)
          call mpi_recv(atot(1,n0),3*natmt,mpi_real8
     &         ,ixyz,itag+7,mpi_md_world,istat,ierr)
          call mpi_recv(chgtot(n0),natmt,mpi_real8
     &         ,ixyz,itag+8,mpi_md_world,istat,ierr)
          call mpi_recv(chitot(n0),natmt,mpi_real8
     &         ,ixyz,itag+9,mpi_md_world,istat,ierr)
          n0 = n0 + natmt
        enddo
        tspdcmp = tspdcmp +(mpi_wtime()-t0)
c        write(6,'(a,f10.3)') ' time space_comp = ',mpi_wtime() -t0
      else ! myid_md.ne.0
        itag = myid_md*nmpi -nmpi
        call mpi_send(natm,1,mpi_integer,0,itag
     &       ,mpi_md_world,ierr)
        call mpi_send(tag,natm,mpi_real8,0,itag+1
     &       ,mpi_md_world,ierr)
c.....positions should be shifted before passing data
        do i=1,natm
          ra(1:3,i) = ra(1:3,i) + sorg(1:3)
        enddo
        call mpi_send(ra,3*natm,mpi_real8,0,itag+2
     &       ,mpi_md_world,ierr)
c.....positions should be reverted after passing data
        do i=1,natm
          ra(1:3,i) = ra(1:3,i) - sorg(1:3)
        enddo
        call mpi_send(va,3*natm,mpi_real8,0,itag+3
     &       ,mpi_md_world,ierr)
        call mpi_send(epi,natm,mpi_real8,0,itag+4
     &       ,mpi_md_world,ierr)
        call mpi_send(eki,3*3*natm,mpi_real8,0,itag+5
     &       ,mpi_md_world,ierr)
        call mpi_send(strs,3*3*natm,mpi_real8,0,itag+6
     &       ,mpi_md_world,ierr)
        call mpi_send(aa,3*natm,mpi_real8,0,itag+7
     &       ,mpi_md_world,ierr)
        call mpi_send(chg,natm,mpi_real8,0,itag+8
     &       ,mpi_md_world,ierr)
        call mpi_send(chi,natm,mpi_real8,0,itag+9
     &       ,mpi_md_world,ierr)
      endif

      end subroutine space_comp
c=======================================================================
      subroutine sort_by_tag(natm,tag,ra,va,eki,epi,strs,chg,chi)
c
c  Sort by tag for output.
c  
      implicit none
      integer,intent(in):: natm
      real(8),intent(inout):: ra(3,natm),va(3,natm),eki(3,3,natm)
     &     ,epi(natm),strs(3,3,natm),tag(natm),chg(natm),chi(natm)

      integer,allocatable:: itag(:)
      real(8),allocatable:: buf(:,:)
      integer:: i
      integer,save:: nsave = 0
      integer,external:: itotOf

      if( .not. allocated(itag) .or. natm.gt.nsave ) then
        if( allocated(itag) ) deallocate(itag,buf)
        nsave = natm
        allocate(itag(natm),buf(20,natm))
      endif

      do i=1,natm
        buf(1:3,i)= ra(1:3,i)
        buf(4:6,i)= va(1:3,i)
        buf(7,i)= eki(1,1,i)
        buf(8,i)= epi(i)
        buf(9:11,i) = strs(1:3,1,i)
        buf(12:14,i)= strs(1:3,2,i)
        buf(15:17,i)= strs(1:3,3,i)
        buf(18,i)= tag(i)
        buf(19,i)= chg(i)
        buf(20,i)= chi(i)
        itag(i)= itotOf(tag(i))
      enddo

      call heapsort_itag(natm,natm,itag,20,buf)

      do i=1,natm
        ra(1:3,i)= buf(1:3,i)
        va(1:3,i)= buf(4:6,i)
        eki(1,1,i)= buf(7,i)
        epi(i)= buf(8,i)
        strs(1:3,1,i)= buf(9:11,i)
        strs(1:3,2,i)= buf(12:14,i)
        strs(1:3,3,i)= buf(15:17,i)
        tag(i)= buf(18,i)
        chg(i)= buf(19,i)
        chi(i)= buf(20,i)
      enddo

      end subroutine sort_by_tag
c=======================================================================
      subroutine error_mpi_stop(cerrmsg)
      use pmdmpi
      implicit none
      character(len=*):: cerrmsg

      integer:: ierr
      
      if(myid_md.eq.0) then
        write(6,'(a)')
     &       '=================================================='
        write(6,'(a)') ' !!! ERROR !!!'
        write(6,'(a)') '   '
     &       //trim(cerrmsg)
        write(6,'(a)')
     &       '=================================================='
      endif
      call mpi_finalize(ierr)
      stop
      end subroutine error_mpi_stop
c=======================================================================
      subroutine cell_info(h)
      implicit none
      real(8),intent(in):: h(3,3)

      real(8):: a,b,c,alpha,beta,gamma
      real(8),parameter:: pi = 3.14159265358979d0
      real(8),external:: sprod

      write(6,*) ''
      write(6,'(a)') " H-matrix:"
      write(6,'("   | ",3f12.3," |")') h(1,1:3)
      write(6,'("   | ",3f12.3," |")') h(2,1:3)
      write(6,'("   | ",3f12.3," |")') h(3,1:3)
      write(6,'(a)') " which means the following lattice vectors:"
      write(6,'(a,"[ ",3f12.3," ]")') '   a = ',h(1:3,1)
      write(6,'(a,"[ ",3f12.3," ]")') '   b = ',h(1:3,2)
      write(6,'(a,"[ ",3f12.3," ]")') '   c = ',h(1:3,3)

      a = dsqrt(sprod(3,h(1:3,1),h(1:3,1)))
      b = dsqrt(sprod(3,h(1:3,2),h(1:3,2)))
      c = dsqrt(sprod(3,h(1:3,3),h(1:3,3)))
      alpha = acos(sprod(3,h(1:3,2),h(1:3,3))/b/c) /pi *180d0
      beta  = acos(sprod(3,h(1:3,1),h(1:3,3))/a/c) /pi *180d0
      gamma = acos(sprod(3,h(1:3,1),h(1:3,2))/a/b) /pi *180d0

      write(6,'(a)') ' Lattice parameters:'
      write(6,'(a,f10.3,a,f7.2,a)') '   a = ',a,' Ang.,  alpha = '
     &     ,alpha,' deg.'
      write(6,'(a,f10.3,a,f7.2,a)') '   b = ',b,' Ang.,  beta  = '
     &     ,beta,' deg.'
      write(6,'(a,f10.3,a,f7.2,a)') '   c = ',c,' Ang.,  gamma = '
     &     ,gamma,' deg.'

      end subroutine cell_info
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make pmd"
c     End:
