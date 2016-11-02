      program parallel_md
c-----------------------------------------------------------------------
c                     Last-modified: <2016-10-30 07:27:59 Ryo KOBAYASHI>
c-----------------------------------------------------------------------
c Spatial decomposition parallel molecular dynamics program.
c
c-----------------------------------------------------------------------
c INPUT FILES:
c ------------
c   in.pmd:     Main input file
c   pmd0000:    Cell info and atom coordinations
c
c OUTPUT FILES:
c -------------
c   out.erg:    Total, kinetic, and potential energy
c   pmd####:    Cell info and atom coordinations of a certain steps
c               in the MD run if required by "flag_out_pmd" in in.pmd.
c   out.stress: Stress component normal to z-upper surface of nanorod,
c               and z-strain of the nanorod.
c-----------------------------------------------------------------------
      use variables
      implicit none
      include "mpif.h"
      include "./params_unit.h"

#ifdef __DISL__
c.....Epot threshold for disl core extraction [Hartree]
      real(8),parameter:: epith = -0.1410d0
#endif


      integer:: i,j,k,l,m,n,ia,ib,is,ifmv,nave,nisl,i_conv
      integer:: ihour,imin,isec
      real(8):: tmp,hscl(3),aai(3),ami,dt2,tave,vi(3),vl(3),epotp
      integer,external:: ispOf,ifmvOf,itotOf
      real(8),external:: box_muller,sprod
      logical:: lconverged = .false.
      logical:: lstrs = .false.
c-----output file names
      character:: cnode*5,cnum*4
c.....FIRE variables
      real(8):: alp_fire,fnorm,vnorm,fdotv
      integer:: num_fire

c.....For tensile test of Al nanorod
      real(8):: strnow,ftop,fbot

c-----initialize the MPI environment
      call mpi_init(ierr)
c-----total number of MD-nodes
      call mpi_comm_size(MPI_COMM_WORLD, nodes_md, ierr)
c-----my rank in MD-nodes
      call mpi_comm_rank(MPI_COMM_WORLD, myid_md, ierr)
      mpi_md_world= MPI_COMM_WORLD


c.....Set fmv as default value before reading 'in.pmd'
      call set_fmv(fmv)

      if( myid_md.eq.0 ) then
        call read_input(10,'in.pmd')
        call check_cmin()
        call write_initial_setting()
c.....Multiply damping factor to avoid multiplication of every step
        if( ifdmp.eq.1 ) then
          fmv(1:3,0:9)= fmv(1:3,0:9) *dmp
        endif
      endif

      call bcast_params()


c-----output every these steps, NOUTERG, NOUTPMD
      if( nerg.ne.0 ) then
        nouterg = max(nstp/nerg,1)
      else
        nouterg = nstp +1
      endif
      if( npmd.ne.0 ) then
        noutpmd = max(nstp/npmd,1)
      else
        noutpmd = nstp +1
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
      myx=myid_md/(ny*nz)
      myy=mod(myid_md/nz,ny)
      myz=mod(myid_md,nz)
c-----reduced node origin
      sorg(1)= anxi*myx
      sorg(2)= anyi*myy
      sorg(3)= anzi*myz

cc.....Each node opens one MD-configuration file
c      write(cnode(1:5),'(i5.5)') myid_md
c      if( trim(ciofmt).eq.'bin' .or. trim(ciofmt).eq.'binary' ) then
c        call read_pmd_bin(20,"0000/pmd"//cnode
c     &       ,namax,natm,h,hunit,tag,ra,va,eki,epi,strs)
c      else if( trim(ciofmt).eq.'ascii' ) then
c        call read_pmd_ascii(20,"0000/pmd"//cnode
c     &     ,namax,natm,h,hunit,tag,ra,va,eki,epi,strs)
c      endif
c.....only 0th-node reads pmd0000 file
      if( myid_md.eq.0 ) then
        if( trim(ciofmt).eq.'bin' .or. trim(ciofmt).eq.'binary' ) then
          call read_pmdtot_bin(20,"pmd0000")
        else if( trim(ciofmt).eq.'ascii' ) then
          call read_pmdtot_ascii(20,"pmd0000")
        endif
      endif
c.....perform space decomposition after reading atomic configuration
      call space_decomp()
      if( myid_md.eq.0 ) then
        write(6,'(a)') " H-matrix:"
        write(6,'(3es15.7)') h(1:3,1,0)
        write(6,'(3es15.7)') h(1:3,2,0)
        write(6,'(3es15.7)') h(1:3,3,0)
      endif
c.....Some conversions
      nis= 0
      do i=1,natm
        ra(1:3,i)= ra(1:3,i) -sorg(1:3)
c-------scaling
        va(1:3,i)= va(1:3,i) *dt
c-------species of atom-i
        nis= max(int(tag(i)),nis)
        if(nis.gt.nismax) then
          write(6,'(a)') " nis.gt.nismax!!"
          stop
        endif
      enddo
c      do i=1,natm
c        write(6,'(i3,i5,3es12.3)') myid_md,i,ra(1:3,i)
c      enddo
c      write(6,'(a,2i8)') " myid,natm=",myid_md,natm

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

c.....Setup for FIRE
      if( ifdmp.eq.2 ) then
        alp_fire = alp0_fire
        am(1:nismax) = 10.d0
        num_fire = 0
        dtmax_fire = dtmfctr_fire * dt
      endif

c-----setup
      call setup(nismax,am,acon,fack,dt)
c-----set HI and SGM
      call boxmat(h,hi,ht,g,gi,gt,vol,sgm)
c-----ntset
      call ntset(myx,myy,myz,nx,ny,nz,nn,sv,myparity,anxi,anyi,anzi)

c-----get total number of particles
c      call mpi_allreduce(natm,ntot,1,MPI_INTEGER,MPI_SUM
c     &     ,mpi_md_world,ierr)
      nisl = nis
      call mpi_allreduce(nisl,nis,1,MPI_INTEGER,mpi_max
     &     ,mpi_md_world,ierr)
      if(myid_md.eq.0) then
        write(6,'(a,i10)') " ntot=",ntot
        write(6,'(a,i10)') " nis =",nis
      endif

c.....Set initial temperature if needed
      if( tinit.gt.1d-5 ) then
        call setv()
        call rm_trans_motion()
      elseif( abs(tinit).le.1d-5 ) then
        va(1:3,1:natm)= 0d0
      endif

      if( ifdmp.eq.2 ) then
        va(1:3,1:natm) = 0d0
      endif

      if( trim(ctctl).eq.'Langevin' ) then
        dt2= dt/2
        tgmm= 1d0/trlx
        do ifmv=1,9
          if( ttgt(ifmv).lt.0d0 ) then
            tfac(ifmv)= -1d0
          else
            tfac(ifmv)= sqrt(2d0*tgmm*dt*(fkb*ev2j)*ttgt(ifmv))
     &           *m2ang/s2fs
          endif
        enddo
      endif

#ifdef __SHEAR__
      call init_shear_stress(shrst,shrfx,natm,tag,ra,h,sorg,acon(1)
     &     ,myid_md,mpi_md_world)
#endif

      tcpu1= mpi_wtime()
      tcom = 0d0

c-----copy RA of boundary atoms
      call check_size_and_parallel()
      call bacopy_new(tcom,sgm,vol,lsb,lsex,nbmax,ra,namax
     &     ,natm,nb,anxi,anyi,anzi,nn,tag,rc,myid_md,myparity
     &     ,lsrc,sv,nex,mpi_md_world)
c-----Make pair list
      call mk_lspr_para(namax,natm,nbmax,nb,nnmax,tag,ra,rc+rbuf,h,hi
     &     ,anxi,anyi,anzi,lspr)
c.....Calc forces
      lstrs = .true.
      call get_force(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &     ,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr
     &     ,mpi_md_world,myid_md,epi,epot0,nismax,acon,lstrs
     &     ,cforce,iprint)
      lstrs = .false.
      epot= epot0
      epotp = 0d0

#ifdef __SHEAR__
      call add_shear_stress(natm,ra,aa,tag,sorg,shrfx)
#endif

#ifdef __DISL__
      call perf_disl_pos_by_pot(epith,natm,ra,h,epi,sorg
     &     ,nodes_md,myid_md,mpi_md_world,0,21)
#endif

c-----calc kinetic energy
      call get_ekin(namax,natm,va,tag,h,nismax,fack,ekin,eki,ekl
     &     ,vmax,mpi_md_world)
      vmaxold=vmax
c.....get_num_dof is called once in a MD run
      call get_num_dof()

      if( trim(cpctl).eq.'Berendsen' .or.
     &     trim(cpctl).eq.'vc-Berendsen' ) then
        if(myid_md.eq.0) then
          write(6,*) ''
          write(6,'(a)') ' variable-cell Berendsen is chosen'
          write(6,'(a)') ' target stress (GPa):'
          write(6,'(3f8.3)') stgt(1,1:3)
          write(6,'(3f8.3)') stgt(2,1:3)
          write(6,'(3f8.3)') stgt(3,1:3)
        endif
        stgt(1:3,1:3)= stgt(1:3,1:3) *gpa2up
      else if( trim(cpctl).eq.'vv-Berendsen' ) then
        if(myid_md.eq.0) then
          write(6,*) ''
          write(6,'(a)') ' variable-volume Berendsen is chosen'
          write(6,'(a,f10.2,a,es12.4,a)')
     &       ' target pressure = ',ptgt,' GPa'
        endif
        ptgt = ptgt *gpa2up
      endif
      call force_isobaric(stgt,ptgt,ah,natm,eki,strs,sgm
     &     ,dt,srlx,stbeta,vol,avol,ptnsr,mpi_md_world,cpctl)
c      ptnsr(1:3,1:3)= ptnsr(1:3,1:3)/vol
      if( myid_md.eq.0 ) then
        write(6,'(a,6f10.3)')' stress tensor [GPa]: '
     &       ,ptnsr(1,1)*up2gpa,ptnsr(2,2)*up2gpa,ptnsr(3,3)*up2gpa
     &       ,ptnsr(2,3)*up2gpa,ptnsr(3,1)*up2gpa,ptnsr(1,2)*up2gpa
      endif

      istp= 0
      istps= istp
      ts= istps*dt
      te= rbuf/(2d0*vmax/dt) +istp*dt
      istpe= te/dt
c.....Considering vmax==0d0, make istpe small only for the first step
      if( istpe.gt.100 ) istpe=100
c      if(myid_md.eq.0) write(6,'(a,i10)') " istpe=",istpe

      if(myid_md.eq.0) then
        write(6,*) ''
        write(6,'(1x,a)') "initial values:"
        write(6,'(1x,a,es22.14,a)') " kinetic energy  =",ekin,' eV'
        write(6,'(1x,a,es22.14,a)') " potential energy=",epot0,' eV'
        nave= 0
        tave= 0d0
        do ifmv=1,9
          if( ndof(ifmv).eq.0 ) cycle
          temp(ifmv)= ekl(ifmv) /ndof(ifmv) /fkb *2d0
          nave= nave +ndof(ifmv)
          write(6,'(1x,a,i1,a,es15.7,a)') " temperature ",ifmv
     &         ,"   =",temp(ifmv),' K'
          tave= tave +temp(ifmv)*ndof(ifmv)
        enddo
        tave= tave/nave
        write(6,'(1x,a,es15.7,a)') " temperature     =",tave,' K'
        write(6,'(1x,a,f8.3,a)') " stress          =",
     &       (ptnsr(1,1)+ptnsr(2,2)+ptnsr(3,3))/3*up2gpa,' GPa'
        

        if( trim(cpctl).eq.'Berendsen' .or.
     &       trim(cpctl).eq.'vc-Berendsen' .or.
     &       trim(cpctl).eq.'vv-Berendsen') then
          write(6,'(a,i8,2f10.2,es15.7,2es11.3)')
     &         " istp,time,temp,vol,prss="
     &         ,istp,tcpu,tave,epot,vol
     &         ,(ptnsr(1,1)+ptnsr(2,2)+ptnsr(3,3))*up2gpa/3
        else
          write(6,'(a,i8,2f10.2,es15.7)')
     &         " istp,time,epot,temp="
     &         ,istp,tcpu,tave,epot
        endif
        write(6,'(a,6f10.3)')' stress tensor [GPa]: '
     &       ,ptnsr(1,1)*up2gpa,ptnsr(2,2)*up2gpa,ptnsr(3,3)*up2gpa
     &       ,ptnsr(2,3)*up2gpa,ptnsr(3,1)*up2gpa,ptnsr(1,2)*up2gpa
        write(6,*) ''
      endif

c-----initialize the counter for output
      iocntpmd=0
      iocnterg=0
c.....I do not know why I need to write initial configuration,
c.....since initial configuration is already given as pmd0000.
cc-----output initial configuration
c      if( ifpmd.eq.1 ) then
c        call system("mkdir -p "//cnum)
c        if( trim(ciofmt).eq.'bin' .or. trim(ciofmt).eq.'binary' ) then
c          call write_pmd_bin(20,cnum//"/pmd"//cnode
c     &         ,natm,h,hunit,tag,ra,va,eki,epi,strs,sorg,dt)
c        elseif( trim(ciofmt).eq.'ascii' ) then
c          call write_pmd_ascii(20,cnum//"/pmd"//cnode
c     &         ,natm,h,hunit,tag,ra,va,eki,epi,strs,sorg,dt)
c        endif
c      endif

      if(myid_md.eq.0) then
c.....write out energies
        open(ioerg,file="out.erg",status='replace')
        write(ioerg,'(a)') '# istp,etot-epot0[eV],ekin,epot-epot0'
     &       //',temp[K],vol[Ang^3],pressure[GPa]'
        write(ioerg,'(a,es16.7e3,a)') '# epot0 =',epot0,' [eV]'
        write(ioerg,'(i10,3es16.7e3,f10.2,2es16.7e3)') istp
     &       ,ekin+epot0-epot0,ekin,epot0-epot0,tave
     &       ,vol
     &       ,(ptnsr(1,1)+ptnsr(2,2)+ptnsr(3,3))/3*up2gpa
        call flush(ioerg)
c.....write out temperatures
        open(iotemp,file='out.temperature',status='replace')
        write(iotemp,'(a)') '# istp, temperature[0-9]'
        ediff(1:9)= 0d0
        write(iotemp,'(i10,18es16.7e3)') istp,temp(1:9),ediff(1:9)
        call flush(iotemp)
        
c.....write out stresses
        if(czload_type.eq.'atoms' .or. czload_type.eq.'box') then
          open(iostrs,file='out.stress',status='replace')
          write(iostrs,'(a)') '#  istp,   strnow,      ftop,     fbot'
        endif
      endif

c.....temperature distribution along x
      if( ltdst ) then
        if( myid_md.eq.0 ) then
          open(iotdst,file='out.temp-dist')
          write(iotdst,'(a)') '#   x-pos,   temp   '
        endif
        call calc_temp_dist()
      endif

c.....Set al
      al(1)= h(1,1,0)
      al(2)= h(2,2,0)
      al(3)= h(3,3,0)

      i_conv = 0

c-----MD loop with velocity-Verlet--------------------------------------
      do istp=1,nstp

c.....In case of isobaric MD, lstrs has to be always TRUE.
        if( trim(cpctl).eq.'Berendsen' .or.
     &       trim(cpctl).eq.'vc-Berendsen' .or.
     &       trim(cpctl).eq.'vv-Berendsen') then
          lstrs = .true.
        else
          lstrs = .false.
        endif

c-------first kick of velocities
        va(1:3,1:natm)=va(1:3,1:natm) +aa(1:3,1:natm)

        if( ifdmp.eq.2 ) then
          call vfire(num_fire,alp_fire)
        endif

c-------multiply fmv or damping
        do i=1,natm
          l= int(mod(tag(i)*10,10d0))
          va(1:3,i)=va(1:3,i) *fmv(1:3,l)
        enddo


c-------update positions
        ra(1:3,1:natm)=ra(1:3,1:natm) +va(1:3,1:natm)

        if( czload_type.eq.'atoms' ) then
          call zload_atoms(natm,ra,tag,nstp,strfin,strnow
     &         ,sorg,myid_md,mpi_md_world)
        else if( czload_type.eq.'box' ) then
          call zload_box(natm,nstp,istp,dt,strfin,strnow,h,myid_md)
        endif

        if( trim(cpctl).eq.'Berendsen' .or.
     &       trim(cpctl).eq.'vc-Berendsen' .or.
     &       trim(cpctl).eq.'vv-Berendsen' ) then
          h(1:3,1:3,0)= matmul(ah,h(1:3,1:3,0))
c.....Reset matrices
          call boxmat(h,hi,ht,g,gi,gt,vol,sgm)
        endif

c.....Update te considering change of vmax
        te=(te-ts)*vmaxold/vmax +ts
        istpe= te/dt

c.....Update pair list and boundary atoms
c.....if making new pair list is needed.
        if( istp.ge.istpe .or.
     &       (ifpmd.eq.1.and.mod(istp,noutpmd).eq.0) ) then
c.....Move atoms that cross the boundary
          call bamove_stress(tcom,namax,nbmax,natm,ra,va,tag,strs
     &         ,anxi,anyi,anzi,myid_md,nn,sv,myparity,mpi_md_world)
c.....Copy RA of boundary atoms
          call bacopy_new(tcom,sgm,vol,lsb,lsex,nbmax,ra,namax,natm,nb
     &         ,anxi,anyi,anzi,nn,tag,rc,myid_md,myparity,lsrc,sv,nex
     &         ,mpi_md_world)
c.....Make pair list
          call mk_lspr_para(namax,natm,nbmax,nb,nnmax,tag,ra,rc+rbuf,h
     &         ,hi,anxi,anyi,anzi,lspr)
          ts= istp*dt
          te=rbuf/(2d0*vmax/dt) +ts
          istpe= te/dt
c          if(myid_md.eq.0) write(6,'(a,i10)')
c     &         " Next update step will be ",istpe
        else
c.....Copy RA of boundary atoms determined by 'bacopy'
          call bacopy_fixed_new(tcom,sgm,vol,lsb,lsex,nbmax,ra,namax
     &         ,natm,nb,anxi,anyi,anzi,nn,tag,rc,myid_md,myparity,lsrc
     &         ,sv,nex,mpi_md_world)
        endif
        
        if(ifpmd.eq.1.and. mod(istp,noutpmd).eq.0 )then
          lstrs = .true.
        endif
c-------Calc forces
        call get_force(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &       ,tcom,nb,nbmax,lsb,nex,lsrc,myparity,nn,sv,rc,lspr
     &       ,mpi_md_world,myid_md,epi,epot,nismax,acon,lstrs
     &       ,cforce,iprint)
        lstrs = .false.

#ifdef __SHEAR__
        call add_shear_stress(natm,ra,aa,tag,sorg,shrfx)
#endif

c.....Second kick of velocities
        if( ctctl.eq.'Langevin' ) then
c.....Langevin thermostat with Mannella integrator
          hscl(1:3)= 0d0
          do l=1,3
            hscl(l)= sqrt(h(1,l,0)**2 +h(2,l,0)**2 +h(3,l,0)**2)
          enddo
          do i=1,natm
            ifmv= int(mod(tag(i)*10,10d0))
            if( ifmv.eq.0 .or. tfac(ifmv).lt.0d0 ) then
              va(1:3,i)=va(1:3,i) +aa(1:3,i)
              cycle
            else
              is= int(tag(i))
              ami= am(is)
              aai(1:3)= 0d0
              tmp= tfac(ifmv) /sqrt(2d0*ami*amu2kg)*dt
              do l=1,3
                aai(l)= -va(l,i)*tgmm*dt2 +tmp*box_muller()/hscl(l)
              enddo
              va(1:3,i)=va(1:3,i) +aa(1:3,i)
c              aa(1:3,i)= aa(1:3,i) +aai(1:3)
c.....accumulate energy difference
              vi(1:3)= h(1:3,1,0)*va(1,i)
     &             +h(1:3,2,0)*va(2,i)
     &             +h(1:3,3,0)*va(3,i)
              vl(1:3)= h(1:3,1,0)*aai(1)
     &             +h(1:3,2,0)*aai(2)
     &             +h(1:3,3,0)*aai(3)
              ediff(ifmv)= ediff(ifmv) +fack(is)
     &             *(2d0*sprod(vi,vl)+sprod(vl,vl))
              va(1:3,i)= va(1:3,i) +aai(1:3)
            endif
          enddo
        else
          va(1:3,1:natm)=va(1:3,1:natm) +aa(1:3,1:natm)
        endif

        

c.....Calc kinetic energy
        vmaxold= vmax
        call get_ekin(namax,natm,va,tag,h,nismax,fack,ekin,eki,ekl
     &       ,vmax,mpi_md_world)

        call force_isobaric(stgt,ptgt,ah,natm,eki,strs,sgm
     &       ,dt,srlx,stbeta,vol,avol,ptnsr,mpi_md_world,cpctl)

c-------temperature control by Berendsen thermostat
        if( trim(ctctl).eq.'Berendsen' ) then
          tfac(1:9)= 0d0
c.....if final temperature is assigned,
c.....target temperatures are forced to be set intermediate temperatures
          if( tfin.gt.0d0 ) then
            ttgt(1:9) = tinit +(tfin-tinit)*istp/nstp
          endif
          do ifmv=1,9
            if(ndof(ifmv).le.0 .or. ttgt(ifmv).lt.0d0 ) cycle
            temp(ifmv)= ekl(ifmv) *2d0 /fkb /ndof(ifmv)
            if( abs(ttgt(ifmv)-temp(ifmv))/temp(ifmv).gt.100d0 ) then
              tfac(ifmv)= sqrt(1d0 +dt/trlx*100d0 )
            else
              tfac(ifmv)= sqrt(1d0 +dt/trlx*(ttgt(ifmv)-temp(ifmv))
     &             /temp(ifmv))
            endif
          enddo
          do i=1,natm
            ifmv= int(mod(tag(i)*10,10d0))
            if( ifmv.le.0 .or. ttgt(ifmv).lt.0d0 ) cycle
            va(1:3,i)= va(1:3,i) *tfac(ifmv)
c.....accumulate energy difference
            is= int(tag(i))
            ediff(ifmv)= ediff(ifmv) +(tfac(ifmv)**2-1d0)*fack(is)
          enddo
        endif

c.....temperature distribution along x
        if( ltdst ) call calc_temp_dist()

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
          tave= tave/nave

          if(czload_type.eq.'atoms') then
            call get_stress_on_base(natm,ra,aa,tag,h,acon(1),ftop,fbot
     &           ,sorg,myid_md,mpi_md_world)
          endif

          ediff(1:9)= ediff(1:9) /nouterg
          ediff0(1:9)= 0d0
          call mpi_reduce(ediff,ediff0,9
     &         ,mpi_double_precision,mpi_sum,0,mpi_md_world,ierr)
          
          if(myid_md.eq.0) then
c.....write energies
            write(ioerg,'(i10,3es16.7e3,f10.2,2es16.7e3)') istp
     &           ,ekin+epot-epot0,ekin,epot-epot0,tave
     &           ,vol
     &           ,(ptnsr(1,1)+ptnsr(2,2)+ptnsr(3,3))/3*up2gpa
            call flush(ioerg)
c.....write temperature
            write(iotemp,'(i10,18es16.7e3)') istp,temp(1:9),ediff0(1:9)
            call flush(iotemp)
            ediff(1:9)= 0d0

            if(czload_type.eq.'atoms') then
              write(iostrs,'(i8,3es15.7)') istp,strnow,ftop,fbot
              call flush(iostrs)
            endif

          endif
c---------output step, time, and temperature
          tcpu= mpi_wtime() -tcpu1
          if(myid_md.eq.0) then
            if( trim(cpctl).eq.'Berendsen' .or.
     &           trim(cpctl).eq.'vc-Berendsen' .or.
     &           trim(cpctl).eq.'vv-Berendsen' ) then
              write(6,'(a,i8,2f10.2,es15.7,2es11.3)')
     &             " istp,time,temp,epot,vol,prss="
     &             ,istp,tcpu,tave,epot,vol
     &             ,(ptnsr(1,1)+ptnsr(2,2)+ptnsr(3,3))/3*up2gpa
              write(6,'(a)') ' H-matrix:' !,h(1:3,1:3,0)
              write(6,'(3f15.7)') h(1:3,1:3,0)
              
              write(6,'(a,6f8.3)') ' stress (GPa):'
     &             ,ptnsr(1,1)*up2gpa ,ptnsr(2,2)*up2gpa
     &             ,ptnsr(3,3)*up2gpa ,ptnsr(2,3)*up2gpa
     &             ,ptnsr(1,3)*up2gpa ,ptnsr(1,2)*up2gpa
            else
              write(6,'(a,i8,2f10.2,es15.7)')
     &             " istp,time,temp,epot="
     &             ,istp,tcpu,tave,epot
            endif
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
          i_conv = i_conv + 1
          if( i_conv.ge.n_conv ) then
            if( myid_md.eq.0 ) then
              write(6,'(a,i6,a,es10.3,a,i3,a)')
     &             ' damped MD converged with ',istp,
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
        if(ifpmd.eq.1.and.
     &       (mod(istp,noutpmd).eq.0.or.lconverged) )then
c---------decide pmd-file name
          iocntpmd=iocntpmd+1
          write(cnum(1:4),'(i4.4)') iocntpmd
c          call system("mkdir -p "//cnum)
          call space_comp()
          if( myid_md.eq.0 ) then
            if( ifsort.eq.1 ) call sort_by_tag(ntot,tagtot,rtot,vtot
     &           ,ekitot,epitot,stot)
            if( trim(ciofmt).eq.'bin' .or. trim(ciofmt).eq.'binary' )
     &           then
              call write_pmdtot_bin(20,"pmd"//cnum)
            elseif( trim(ciofmt).eq.'ascii' ) then
              call write_pmdtot_ascii(20,"pmd"//cnum)
            endif
          endif
c          if( trim(ciofmt).eq.'bin' .or. trim(ciofmt).eq.'binary' ) then
c            call write_pmd_bin(20,cnum//"/pmd"//cnode
c     &           ,natm,h,hunit,tag,ra,va,eki,epi,strs,sorg,dt)
c          elseif( trim(ciofmt).eq.'ascii' ) then
c            call write_pmd_ascii(20,cnum//"/pmd"//cnode
c     &           ,natm,h,hunit,tag,ra,va,eki,epi,strs,sorg,dt)
c          endif
        endif

        if( lconverged ) exit
      enddo

      tcpu2= mpi_wtime()

      if(myid_md.eq.0) then
        close(ioerg)
        close(iotemp)
        if(czload_type.eq.'atoms' .or. czload_type.eq.'box') then
          close(iostrs)
        endif
        if( ltdst ) close(iotdst)
      endif

c.....write energy and forces only for fitpot (2nd digit of iprint==1)
      call write_force(21,'.pmd',natm,h,tag,aa,acon
     &     ,nismax,epot,myid_md,mpi_md_world,nxyz,ntot,tagtot,atot
     &     ,ptnsr)
      
      tcpu= tcpu2 -tcpu1
      if(myid_md.eq.0) then
        write(6,*) ''
        write(6,'(1x,a)') "final values:"
        write(6,'(1x,a,es22.14,a)') " kinetic energy  =",ekin,' eV'
        write(6,'(1x,a,es22.14,a)') " potential energy=",epot,' eV'
        nave= 0
        tave= 0d0
        do ifmv=1,9
          if( ndof(ifmv).le.0 ) cycle
          temp(ifmv)= ekl(ifmv) *2d0 /fkb /ndof(ifmv)
          nave= nave +ndof(ifmv)
          write(6,'(1x,a,i1,a,es15.7,a)') " temperature ",ifmv
     &         ,"   =",temp(ifmv),' K'
          tave= tave +temp(ifmv)*ndof(ifmv)
        enddo
        tave= tave/nave
        write(6,'(1x,a,es15.7,a)') " temperature     =",tave,' K'
        write(6,'(1x,a,f8.3,a)') " stress          =",
     &       (ptnsr(1,1)+ptnsr(2,2)+ptnsr(3,3))/3*up2gpa,' GPa'
        write(6,*) ''
        write(6,'(1x,a,f10.2)') "time space decomp =",tspdcmp
        write(6,'(1x,a,f10.2)') "time comm         =",tcom

        ihour = int(tcpu/3600)
        imin  = int((tcpu-ihour*3600)/60)
        isec  = int(tcpu -ihour*3600 -imin*60)
        write(6,'(1x,a,f10.2,a,i3,"h",i2.2,"m",i2.2,"s")')
     &       "time              =",tcpu,
     &       " sec  = ",ihour,imin,isec
      endif

      call mpi_finalize(ierr)

      end program parallel_md
c=======================================================================
      subroutine write_initial_setting()
      use variables
      implicit none 
      integer:: i

      write(6,'(a)') '========== INITIAL SETTING ==========='
      write(6,'(2x,a25,2x,a)') 'io_format',ciofmt
      write(6,'(2x,a)') ''
      write(6,'(2x,a25,i8)')   'num_nodes_x',nx
      write(6,'(2x,a25,i8)')   'num_nodes_y',ny
      write(6,'(2x,a25,i8)')   'num_nodes_z',nz
      write(6,'(2x,a)') ''
      write(6,'(2x,a25,es11.3)') 'time_interval',dt
      write(6,'(2x,a25,i8)')   'num_iteration',nstp
      write(6,'(2x,a25,i8)')   'num_out_energy',nerg
      write(6,'(2x,a)') ''
      write(6,'(2x,a25,i8)')   'flag_out_pmd',ifpmd
      write(6,'(2x,a25,i8)')   'num_out_pmd',npmd
      write(6,'(2x,a)') ''
      write(6,'(2x,a25,2x,a)') 'force_type',cforce
      write(6,'(2x,a25,es15.7)') 'cutoff_radius',rc
      write(6,'(2x,a25,es15.7)') 'cutoff_buffer',rbuf
      write(6,'(2x,a)') ''
      if( cmin.ne.'' ) then
        write(6,'(2x,a25,2x,a)') 'minimization',cmin
      endif
      write(6,'(2x,a25,i8)') 'flag_damping',ifdmp
      write(6,'(2x,a25,es11.3)') 'damping_coeff',dmp
      write(6,'(2x,a25,i8)')   'min_iteration',minstp
      write(6,'(2x,a)') ''
c.....temperature control
      write(6,'(2x,a25,f12.4)') 'initial_temperature',tinit
      if( tfin.ge.0d0 ) then
        write(6,'(2x,a25,f12.4)') 'final_temperature',tfin
      endif
      write(6,'(2x,a25,2x,a)') 'temperature_control',ctctl
      do i=1,9
        write(6,'(2x,a25,i3,f12.4)') 'temperature_target',i,ttgt(i)
      enddo
      write(6,'(2x,a25,f12.4)') 'temperature_relax_time',trlx
      write(6,'(2x,a)') ''
c.....temperature distribution
      write(6,'(2x,a25,l8)') 'flag_temp_dist',ltdst
      write(6,'(2x,a25,i8)') 'num_temp_dist',ntdst
      write(6,'(2x,a)') ''
c.....pressure control
      write(6,'(2x,a25,2x,a)') 'stress_control',cpctl
      if( trim(cpctl).eq.'Berendsen' .or.
     &     trim(cpctl).eq.'vc-Berendsen' ) then
        write(6,'(2x,a)') 'stress_target'
        write(6,'(5x,3es11.3)') stgt(1,1:3)
        write(6,'(5x,3es11.3)') stgt(2,1:3)
        write(6,'(5x,3es11.3)') stgt(3,1:3)
        write(6,'(2x,a)') ''
      else if( trim(cpctl).eq.'vv-Berendsen' ) then
        write(6,'(2x,a25,f12.3)') 'pressure_target',ptgt
      endif
c.....strain control
      write(6,'(2x,a25,f12.4)') 'final_strain',strfin
      write(6,'(2x,a)') ''
c.....velocity multiplying factor
      write(6,'(2x,a)') 'factor_direction'
      do i=0,9
        write(6,'(4x,i2,3es11.3)') i,fmv(1:3,i)
      enddo
      write(6,'(2x,a)') ''
      write(6,'(2x,a)') 'mass'
      do i=1,nismax
        write(6,'(4x,i2,f10.3)') i,am(i)
      enddo
      write(6,'(a)') '======================================'

      end subroutine write_initial_setting
c=======================================================================
      subroutine bcast_params()
      use variables
      implicit none
      include 'mpif.h'

c-----Broadcast input parameters to all nodes
      call mpi_bcast(nx,1,MPI_INTEGER,0,mpi_md_world,ierr)
      call mpi_bcast(ny,1,MPI_INTEGER,0,mpi_md_world,ierr)
      call mpi_bcast(nz,1,MPI_INTEGER,0,mpi_md_world,ierr)
      call mpi_bcast(nstp,1,MPI_INTEGER,0,mpi_md_world,ierr)
      call mpi_bcast(minstp,1,MPI_INTEGER,0,mpi_md_world,ierr)
      call mpi_bcast(dt,1,MPI_DOUBLE_PRECISION,0,mpi_md_world,ierr)
      call mpi_bcast(rc,1,MPI_DOUBLE_PRECISION,0,mpi_md_world,ierr)
      call mpi_bcast(rbuf,1,MPI_DOUBLE_PRECISION,0,mpi_md_world,ierr)
      call mpi_bcast(dmp,1,MPI_DOUBLE_PRECISION,0,mpi_md_world,ierr)
      call mpi_bcast(eps_conv,1,MPI_DOUBLE_PRECISION,0,
     &     mpi_md_world,ierr)
      call mpi_bcast(n_conv,1,mpi_integer,0,mpi_md_world,ierr)
      call mpi_bcast(tinit,1,mpi_double_precision,0,mpi_md_world,ierr)
      call mpi_bcast(tfin,1,mpi_double_precision,0,mpi_md_world,ierr)
      call mpi_bcast(ctctl,20,mpi_character,0,mpi_md_world,ierr)
      call mpi_bcast(ttgt,9,MPI_DOUBLE_PRECISION,0,mpi_md_world,ierr)
      call mpi_bcast(trlx,1,MPI_DOUBLE_PRECISION,0,mpi_md_world,ierr)
      call mpi_bcast(nerg,1,MPI_INTEGER,0,mpi_md_world,ierr)
      call mpi_bcast(ifpmd,1,MPI_INTEGER,0,mpi_md_world,ierr)
      call mpi_bcast(npmd,1,MPI_INTEGER,0,mpi_md_world,ierr)
      call mpi_bcast(ifsort,1,mpi_integer,0,mpi_md_world,ierr)
      call mpi_bcast(ifdmp,1,MPI_INTEGER,0,mpi_md_world,ierr)
      call mpi_bcast(iprint,1,mpi_integer,0,mpi_md_world,ierr)
      call mpi_bcast(fmv,30,mpi_double_precision,0,mpi_md_world,ierr)
      call mpi_bcast(shrst,1,mpi_double_precision,0,mpi_md_world,ierr)
      call mpi_bcast(cpctl,20,mpi_character,0,mpi_md_world,ierr)
      call mpi_bcast(ptgt,1,mpi_double_precision,0,mpi_md_world,ierr)
      call mpi_bcast(srlx,1,mpi_double_precision,0,mpi_md_world,ierr)
      call mpi_bcast(stgt,9,mpi_double_precision,0,mpi_md_world,ierr)
      call mpi_bcast(czload_type,5,mpi_character,0,mpi_md_world,ierr)
      call mpi_bcast(strfin,1,mpi_double_precision,0,mpi_md_world,ierr)
      call mpi_bcast(am,nismax,mpi_double_precision,0,mpi_md_world,ierr)
      call mpi_bcast(ciofmt,6,mpi_character,0,mpi_md_world,ierr)
      call mpi_bcast(cforce,20,mpi_character,0,mpi_md_world,ierr)
c.....NEMD
      call mpi_bcast(ltdst,1,mpi_logical,0,mpi_md_world,ierr)
      if( ltdst ) then
        call mpi_bcast(ntdst,1,mpi_integer,0,mpi_md_world,ierr)
        if( mod(ntdst,nx).ne.0 )
     &       call error_mpi_stop('mod(ntdst,nx).ne.0')
        allocate(tdst(ntdst),nadst(ntdst))
      endif

      end subroutine bcast_params
c=======================================================================
      subroutine error_mpi_stop(cerrmsg)
      use variables
      implicit none
      character(len=*):: cerrmsg
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
      subroutine setup(nismax,am,acon,fack,dt)
      implicit none
      include "params_unit.h"
      integer,intent(in):: nismax
      real(8),intent(in):: dt
      real(8),intent(out):: acon(nismax),fack(nismax),am(nismax)

      integer:: i
      real(8):: tmp

      am(1:nismax)= am(1:nismax)

c.....umass*Ang^2/fs^2 ---> kg*m^2/sec^2 (=Joule)
      tmp= amu2kg *ang2m**2 /fs2s**2
c      write(6,'(a,es12.4,a,es12.4,a)') ' umass*Ang^2/fs^2 ='
c     &     ,tmp,' Joule =',tmp*j2ev,' eV'

c-----prefactors for normalized accelerration
      do i=1,nismax
        acon(i)=  ev2j/ang2m *0.5d0 *(dt*fs2s)**2 /(am(i)*amu2kg) *m2ang
      enddo

c-----prefactors for kinetic energy in eV
      do i=1,nismax
        fack(i)= 0.5d0*(am(i)*amu2kg)/(dt*fs2s)**2 *ang2m**2 *j2ev
      enddo

      end subroutine setup
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
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine set_fmv(fmv)
c
c Set default fmv values which might be override
c
      implicit none
      real(8),intent(out):: fmv(3,0:9)

c-----set fmv(1:3,ifmv) to be multiplied to the velocities
      fmv(1:3,0)= (/ 0d0, 0d0, 0d0 /) ! fix
      fmv(1:3,1)= (/ 1d0, 1d0, 1d0 /) ! free move
      fmv(1:3,2)= (/ 1d0, 1d0, 0d0 /) ! xy-only
      fmv(1:3,3)= (/ 1d0, 1d0, 1d0 /) ! free move
      fmv(1:3,4)= (/ 1d0, 1d0, 1d0 /) ! free move
      fmv(1:3,5)= (/ 1d0, 1d0, 1d0 /) ! free move
      fmv(1:3,6)= (/ 1d0, 1d0, 1d0 /) ! free move
      fmv(1:3,7)= (/ 1d0, 1d0, 1d0 /) ! free move
      fmv(1:3,8)= (/ 1d0, 1d0, 1d0 /) ! free move
      fmv(1:3,9)= (/ 1d0, 1d0, 1d0 /) ! free move

      end subroutine set_fmv
c=======================================================================
      subroutine get_ekin(namax,natm,va,tag,h,nismax,fack,ekin,eki,ekl
     &     ,vmax,mpi_md_world)
      implicit none 
      include "mpif.h"
      integer,intent(in):: namax,natm,mpi_md_world,nismax
      real(8),intent(in):: va(3,namax),h(3,3),fack(nismax)
     &     ,tag(namax)
      real(8),intent(out):: ekin,eki(3,3,namax),vmax,ekl(9)
c-----locals
      integer:: i,ierr,is,ixyz,jxyz,imax,ifmv
      real(8):: ekinl,x,y,z,v(3),v2,vmaxl,ekll(9)

      ekinl=0d0
      eki(1:3,1:3,1:natm)= 0d0
      ekll(1:9)= 0d0
      vmaxl= 0d0

      do i=1,natm
        is= int(tag(i))
        ifmv= int(mod(tag(i)*10,10d0))
        if( ifmv.eq.0 ) cycle
        x= va(1,i)
        y= va(2,i)
        z= va(3,i)
        v(1:3)= h(1:3,1)*x +h(1:3,2)*y +h(1:3,3)*z
c.....Tensor form eki
        do jxyz=1,3
          do ixyz=1,3
            eki(ixyz,jxyz,i)= v(ixyz)*v(jxyz)
          enddo
        enddo
        v2= eki(1,1,i) +eki(2,2,i) +eki(3,3,i)
        eki(1:3,1:3,i)=eki(1:3,1:3,i)*fack(is)
        ekinl=ekinl +eki(1,1,i) +eki(2,2,i) +eki(3,3,i)
c.....ekin of each ifmv
        ekll(ifmv)= ekll(ifmv) +eki(1,1,i) +eki(2,2,i) +eki(3,3,i)
c.....Find max speed
        if( v2.gt.vmaxl ) imax=i
        vmaxl=max(vmaxl,v2)
      enddo

      ekin= 0d0
      call mpi_allreduce(ekinl,ekin,1,MPI_DOUBLE_PRECISION
     &     ,MPI_SUM,mpi_md_world,ierr)
      ekl(1:9)= 0d0
      call mpi_allreduce(ekll,ekl,9,mpi_double_precision
     &     ,mpi_sum,mpi_md_world,ierr)
      vmax= 0d0
      call mpi_allreduce(vmaxl,vmax,1,mpi_double_precision
     &     ,mpi_max,mpi_md_world,ierr)
      vmax=sqrt(vmax)

      end subroutine get_ekin
c=======================================================================
      subroutine calc_temp_dist()
      use variables
      implicit none
      include 'mpif.h'
      include './params_unit.h'
      integer:: i,ix,itmp
      real(8):: dx,xi,tmp
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
      subroutine get_num_dof()
      use variables
      implicit none
      include 'mpif.h'
      integer:: i,l,k,ndofl(9)
      real(8),parameter:: deps= 1d-12

      ndofl(1:9)= 0d0
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
      if(myid_md.eq.0) write(6,'(a,9i9)') ' ndof=',ndof(1:9)
c      write(6,'(a,i3,a,9i6)') ' myid=',myid_md,': ndof=',ndof(1:9)
      return
      end subroutine get_num_dof
c=======================================================================
      subroutine check_size_and_parallel()
      use variables
      implicit none
      include 'mpif.h'

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
      subroutine bacopy(tcom,sgm,vol,lsb,nbmax,ra,namax,natm,nb
     &     ,anxi,anyi,anzi,nn,tag,rc,myid,myparity,lsrc,sv
     &     ,mpi_md_world)
c-----------------------------------------------------------------------
c     Exchanges boundary-atom data among neighbor nodes: tag and ra
c-----------------------------------------------------------------------
      implicit none
      include 'mpif.h'
      integer,intent(in):: namax,nbmax,nn(6)
     &     ,myid,myparity(3),mpi_md_world
      integer,intent(inout):: natm,nb,lsb(0:nbmax,6),lsrc(6)
      real(8),intent(in):: sv(3,6),sgm(3,3),vol,anxi,anyi,anzi,rc
      real(8),intent(inout):: tag(namax),ra(3,namax),tcom

c      integer:: status(MPI_STATUS_SIZE)
      integer:: i,j,kd,kdd,kul,kuh,ku,ierr
      integer:: nav,maxna,maxb,inode,nsd,nrc,nbnew
      real(8):: tcom1,tcom2,vala,valb,valc,rca,rcb,rcc,xi(3)
      logical,external:: bbd
      real(8),save,allocatable:: dbuf(:,:),dbufr(:,:)
      logical,save:: l1st=.true.

      if( l1st ) then
        allocate(dbuf(4,nbmax),dbufr(4,nbmax))
        l1st=.false.
      endif

c-----reset the num of "received" boundary atoms
      nbnew=0

c-----calculate the cut-off lengths
      vala=dsqrt(sgm(1,1)**2+sgm(2,1)**2+sgm(3,1)**2)/vol
      valb=dsqrt(sgm(1,2)**2+sgm(2,2)**2+sgm(3,2)**2)/vol
      valc=dsqrt(sgm(1,3)**2+sgm(2,3)**2+sgm(3,3)**2)/vol
      rca=rc*vala
      rcb=rc*valb
      rcc=rc*valc

c-----loop over x, y, & z directions
      do 100 kd=1,3

c-------No. of to-be-copied atoms, LSB(0,)
        do kdd= -1,0
          lsb(0,2*kd+kdd)=0
        enddo
c-------Scan all the residents & copies
        do i=1,natm+nbnew
          xi(1:3)= ra(1:3,i)
c---------For low & high directions
          kul=2*kd-1
          kuh=2*kd
c---------Lower neighbor
          if (bbd(xi(1),xi(2),xi(3),rca,rcb,rcc
     &         ,kul,anxi,anyi,anzi)) then
            lsb(0,kul)=lsb(0,kul)+1
            lsb(lsb(0,kul),kul)=i
          endif
c---------Higher neighbor
          if(bbd(xi(1),xi(2),xi(3),rca,rcb,rcc
     &           ,kuh,anxi,anyi,anzi)) then
            lsb(0,kuh)=lsb(0,kuh)+1
            lsb(lsb(0,kuh),kuh)=i
          endif
        enddo
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

c-------To calculate the communication time
        tcom1=MPI_WTIME()

        do 200 kdd= -1,0
          ku=2*kd+kdd
          inode=nn(ku)
          nsd=lsb(0,ku)

          call mespasi(inode,myparity(kd),nsd,nrc,1,1,10
     &         ,mpi_md_world)
c---------Store the # of received boundary atoms in LSRC
          lsrc(ku)=nrc

c---------Exchange ra and tag
          do i=1,nsd
            j= lsb(i,ku)
            dbuf(1:3,i)= ra(1:3,j) -sv(1:3,ku)
            dbuf(4,i)  = tag(j)
          enddo
          call mespasd(inode,myparity(kd),dbuf,dbufr,nsd*4,nrc*4,21
     &         ,mpi_md_world)
          do i=1,nrc
            ra(1:3,natm+nbnew+i)= dbufr(1:3,i)
            tag(natm+nbnew+i)   = dbufr(4,i)
          enddo

          call MPI_BARRIER(mpi_md_world,ierr)
c---------increase the # of received boundary atoms
          nbnew=nbnew+nrc
c          write(6,'(a,2i8)') " nbnew,nrc=",nbnew,nrc
200     continue

c-------Add the communication time to COMT
        tcom2=MPI_WTIME()
        tcom=tcom+tcom2-tcom1

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

100   continue

c-----num. of received boundary atoms
      nb=nbnew

c      print *, 'bacopy: nb=',nb

      end subroutine bacopy
c=======================================================================
      subroutine bacopy_new(tcom,sgm,vol,lsb,lsex,nbmax,ra,namax
     &     ,natm,nb,anxi,anyi,anzi,nn,tag,rc,myid,myparity
     &     ,lsrc,sv,nex,mpi_md_world)
c-----------------------------------------------------------------------
c Exchanges boundary-atom data among neighbor nodes: tag and ra
c Normal data passing and repeating the cell are mixed, and
c automatically selected according to the size of the system.
c
c Note: parallelized to smaller than rcut should not happen.
c-----------------------------------------------------------------------
      implicit none
      include 'mpif.h'
      integer,intent(in):: namax,nbmax,nn(6)
     &     ,myid,myparity(3),mpi_md_world
      integer,intent(inout):: natm,nb,lsb(0:nbmax,6),lsex(nbmax,6)
     &     ,lsrc(6),nex(3)
      real(8),intent(in):: sv(3,6),sgm(3,3),vol,anxi,anyi,anzi,rc
      real(8),intent(inout):: tag(namax),ra(3,namax),tcom

c      integer:: status(MPI_STATUS_SIZE)
      integer:: i,j,kd,kdd,kul,kuh,ku,ierr,iex
      integer:: nav,maxna,maxb,inode,nsd,nrc,nbnew
      real(8):: tcom1,tcom2,xi(3),rcv(3),asgm
      logical,external:: bbd
      real(8),save,allocatable:: dbuf(:,:),dbufr(:,:)
      logical,save:: l1st=.true.
      logical:: lshort(3)

      if( l1st ) then
        allocate(dbuf(4,nbmax),dbufr(4,nbmax))
      endif

c-----reset the num of "received" boundary atoms
      nbnew=0

c-----calculate the cut-off lengths
      do kd=1,3
        asgm= dsqrt(sgm(1,kd)**2 +sgm(2,kd)**2 +sgm(3,kd)**2)
        rcv(kd)= rc*asgm/vol
        nex(kd)= int(rcv(kd)) +1
      enddo
      if( l1st .and. myid.eq.0 ) then
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
c          write(6,'(a,i5,i6)') ' ku,lsb(0,ku) = ',kul,lsb(0,kul)
c          write(6,'(a,i5,i6)') ' ku,lsb(0,ku) = ',kuh,lsb(0,kuh)
        else ! long enough for normal boundary-atom copy
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
        endif ! if(nex.gt.1)

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
            do i=1,lsb(0,ku)
              j= lsb(i,ku)
              iex= lsex(i,ku)
              ra(1:3,natm+nbnew+i)= ra(1:3,j)
              ra(kd,natm+nbnew+i)= ra(kd,natm+nbnew+i) +iex
c              write(6,*) 'kd,iex,ku,i,ra=',kd,iex,ku
c     &             ,i,ra(1:3,natm+nbnew+i)
              tag(natm+nbnew+i)= tag(j)
            enddo
            nbnew= nbnew +lsb(0,ku)
          enddo
        else
c-------To calculate the communication time
          tcom1=MPI_WTIME()
          do 200 kdd= -1,0
            ku=2*kd+kdd
            inode=nn(ku)
            nsd=lsb(0,ku)
c            write(6,*) ' kd,kdd,ku,nsd =',kd,kdd,ku,nsd
            call mespasi(inode,myparity(kd),nsd,nrc,1,1,10
     &           ,mpi_md_world)
c---------Store the # of received boundary atoms in LSRC
            lsrc(ku)=nrc

c---------Exchange ra and tag
            do i=1,nsd
              j= lsb(i,ku)
              dbuf(1:3,i)= ra(1:3,j) -sv(1:3,ku)
              dbuf(4,i)  = tag(j)
            enddo
            call mespasd(inode,myparity(kd),dbuf,dbufr,nsd*4,nrc*4,21
     &           ,mpi_md_world)
            do i=1,nrc
              ra(1:3,natm+nbnew+i)= dbufr(1:3,i)
              tag(natm+nbnew+i)   = dbufr(4,i)
            enddo

            call MPI_BARRIER(mpi_md_world,ierr)
c---------increase the # of received boundary atoms
            nbnew=nbnew+nrc
c            write(6,'(a,2i8)') " nbnew,nrc=",nbnew,nrc
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

100   continue

c-----num. of received boundary atoms
      nb=nbnew

c      print *, 'bacopy: nb=',nb
c      do i=1,natm+nb
c        write(6,'(i5,3f15.7)') i,ra(1:3,i)
c      enddo

      l1st=.false.
      end subroutine bacopy_new
c=======================================================================
      subroutine bacopy_fixed(tcom,sgm,vol,lsb,nbmax,ra,namax,natm,nb
     &     ,anxi,anyi,anzi,nn,tag,rc,myid_md,myparity,lsrc,sv
     &     ,mpi_md_world)
c-----------------------------------------------------------------------
c  Exchanges boundary-atom data among neighbor nodes: tag and ra
c  This doesnt search using position, just send & recv data of atoms
c    which were listed by 'bacopy'.
c-----------------------------------------------------------------------
      implicit none
      include 'mpif.h'
      integer,intent(in):: namax,nbmax,nn(6),natm,nb
     &     ,myid_md,myparity(3),mpi_md_world
      integer,intent(inout):: lsb(0:nbmax,6),lsrc(6)
      real(8),intent(in):: sv(3,6),sgm(3,3),vol,anxi,anyi,anzi,rc
      real(8),intent(inout):: tag(namax),ra(3,namax),tcom

c      integer:: status(MPI_STATUS_SIZE)
      integer:: i,j,kd,kdd,ku,ierr
      integer:: inode,nsd,nrc,nbnew
      real(8):: tcom1,tcom2
      real(8),save,allocatable:: dbuf(:,:),dbufr(:,:)
      logical,save:: l1st=.true.

      if( l1st ) then
        allocate(dbuf(4,nbmax),dbufr(4,nbmax))
        l1st=.false.
      endif

      nbnew= 0

c-----loop over x, y, & z directions
      do 100 kd=1,3

c-------To calculate the communication time
        tcom1=MPI_WTIME()

        do 200 kdd= -1,0
          ku=2*kd+kdd
          inode=nn(ku)
          nsd=lsb(0,ku)
          nrc=lsrc(ku)

c---------Exchange ra and tag
          do i=1,nsd
            j= lsb(i,ku)
            dbuf(1:3,i)= ra(1:3,j) -sv(1:3,ku)
            dbuf(4,i)  = tag(j)
          enddo
          call mespasd(inode,myparity(kd),dbuf,dbufr,nsd*4,nrc*4,21
     &         ,mpi_md_world)
          do i=1,nrc
            ra(1:3,natm+nbnew+i)= dbufr(1:3,i)
            tag(natm+nbnew+i)   = dbufr(4,i)
          enddo

          call MPI_BARRIER(mpi_md_world,ierr)
          nbnew=nbnew +nrc
200     continue

c-------Add the communication time to COMT
        tcom2=MPI_WTIME()
        tcom=tcom+tcom2-tcom1

100   continue

      end subroutine bacopy_fixed
c=======================================================================
      subroutine bacopy_fixed_new(tcom,sgm,vol,lsb,lsex,nbmax,ra,namax
     &     ,natm,nb,anxi,anyi,anzi,nn,tag,rc,myid_md,myparity,lsrc,sv
     &     ,nex
     &     ,mpi_md_world)
c-----------------------------------------------------------------------
c  Exchanges boundary-atom data among neighbor nodes: tag and ra
c  This doesnt search using position, just send & recv data of atoms
c    which were listed by 'bacopy'.
c-----------------------------------------------------------------------
      implicit none
      include 'mpif.h'
      integer,intent(in):: namax,nbmax,nn(6),natm,nb
     &     ,myid_md,myparity(3),mpi_md_world,nex(3)
      integer,intent(inout):: lsb(0:nbmax,6),lsex(nbmax,6),lsrc(6)
      real(8),intent(in):: sv(3,6),sgm(3,3),vol,anxi,anyi,anzi,rc
      real(8),intent(inout):: tag(namax),ra(3,namax),tcom

c      integer:: status(MPI_STATUS_SIZE)
      integer:: i,j,kd,kdd,ku,ierr,iex
      integer:: inode,nsd,nrc,nbnew
      real(8):: tcom1,tcom2 !,rcv(3),asgm
      real(8),save,allocatable:: dbuf(:,:),dbufr(:,:)
      logical,save:: l1st=.true.

      if( l1st ) then
        allocate(dbuf(4,nbmax),dbufr(4,nbmax))
        l1st=.false.
      endif

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
              ra(1:3,natm+nbnew+i)= ra(1:3,j)
              ra(kd,natm+nbnew+i)= ra(kd,natm+nbnew+i) +iex
              tag(natm+nbnew+i)= tag(j)
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
              dbuf(1:3,i)= ra(1:3,j) -sv(1:3,ku)
              dbuf(4,i)  = tag(j)
            enddo
            call mespasd(inode,myparity(kd),dbuf,dbufr,nsd*4,nrc*4,21
     &           ,mpi_md_world)
            do i=1,nrc
              ra(1:3,natm+nbnew+i)= dbufr(1:3,i)
              tag(natm+nbnew+i)   = dbufr(4,i)
            enddo

            call MPI_BARRIER(mpi_md_world,ierr)
            nbnew=nbnew +nrc
 200      continue
        endif

c-------Add the communication time to COMT
        tcom2=MPI_WTIME()
        tcom=tcom+tcom2-tcom1

100   continue

      end subroutine bacopy_fixed_new
c=======================================================================
      subroutine bamove(tcom,namax,nbmax,natm,ra,va,tag
     &     ,anxi,anyi,anzi,myid_md,nn,sv,myparity,mpi_md_world)
c-----------------------------------------------------------------------
c  Exchange atoms between neighbor nodes.
c
c  MVQUE(0:NBMAX,6):
c    MVQUE(0,ku) is the # of to-be-moved atoms to neighbor ku;
c    MVQUE(i,ku) is the adress, in IS, of atom i
c----------------------------------------------------------------------
      implicit none
      include 'mpif.h'
      integer,intent(in):: namax,nbmax
      integer,intent(in):: myid_md,mpi_md_world,myparity(3),nn(6)
      real(8),intent(in):: anxi,anyi,anzi,sv(3,6)
      integer,intent(inout):: natm
      real(8),intent(inout):: tcom,ra(3,namax),va(3,namax),tag(namax)

c      integer:: status(MPI_STATUS_SIZE)
      integer:: i,j,ku,kd,kdd,kul,kuh,inode,nsd,nrc,ipt,ierr,is
      integer:: mvque(0:nbmax,6),newim
      real(8):: tcom1,tcom2,xi(3)
      logical,external:: bmv
      real(8),save,allocatable:: dbuf(:,:),dbufr(:,:)
      logical,save:: l1st=.true.

      if( l1st ) then
        allocate(dbuf(7,nbmax),dbufr(7,nbmax))
        l1st=.false.
      endif

c-----newim: num. of new immigrants
      newim= 0
      mvque(0,1:6)= 0

      do 100 kd=1,3

c-------num of to-be-moved atoms
        do i=1,natm+newim
          xi(1:3)= ra(1:3,i)
          kul=2*kd-1
          kuh=2*kd
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
c-----------Eliminate the record of a moved-out atom
            tag(j)= 0d0
          enddo
          call mespasd(inode,myparity(kd),dbuf,dbufr,7*nsd,7*nrc,71
     &         ,mpi_md_world)
          do i=1,nrc
            ra(1:3,natm+newim+i)= dbufr(1:3,i)
            va(1:3,natm+newim+i)= dbufr(4:6,i)
            tag(natm+newim+i)   = dbufr(7,i)
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
        endif
      enddo
c-----Update # of resident atoms
      natm=ipt

      return
      end subroutine bamove
c=======================================================================
      subroutine bamove_stress(tcom,namax,nbmax,natm,ra,va,tag,strs
     &     ,anxi,anyi,anzi,myid_md,nn,sv,myparity,mpi_md_world)
c-----------------------------------------------------------------------
c  Exchange atoms between neighbor nodes.
c
c  MVQUE(0:NBMAX,6):
c    MVQUE(0,ku) is the # of to-be-moved atoms to neighbor ku;
c    MVQUE(i,ku) is the adress, in IS, of atom i
c----------------------------------------------------------------------
      implicit none
      include 'mpif.h'
      integer,intent(in):: namax,nbmax
      integer,intent(in):: myid_md,mpi_md_world,myparity(3),nn(6)
      real(8),intent(in):: anxi,anyi,anzi,sv(3,6)
      integer,intent(inout):: natm
      real(8),intent(inout):: tcom,ra(3,namax),va(3,namax),tag(namax)
     &     ,strs(9,namax)

c      integer:: status(MPI_STATUS_SIZE)
      integer:: i,j,ku,kd,kdd,kul,kuh,inode,nsd,nrc,ipt,ierr,is
      integer:: mvque(0:nbmax,6),newim
      real(8):: tcom1,tcom2,xi(3)
      logical,external:: bmv
      real(8),save,allocatable:: dbuf(:,:),dbufr(:,:)
      logical,save:: l1st=.true.

      if( l1st ) then
        allocate(dbuf(16,nbmax),dbufr(16,nbmax))
        l1st=.false.
      endif

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
            dbuf(8:16,i)= strs(1:9,j)
c-----------Eliminate the record of a moved-out atom
            tag(j)= 0d0
            strs(1:9,j)= 0d0
          enddo
          call mespasd(inode,myparity(kd),dbuf,dbufr,16*nsd,16*nrc,71
     &         ,mpi_md_world)
          do i=1,nrc
            ra(1:3,natm+newim+i)= dbufr(1:3,i)
            va(1:3,natm+newim+i)= dbufr(4:6,i)
            tag(natm+newim+i)   = dbufr(7,i)
            strs(1:9,natm+newim+i)= dbufr(8:16,i)
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
          strs(1:9,ipt)= strs(1:9,i)
        endif
      enddo
c-----Update # of resident atoms
      natm=ipt

      return
      end subroutine bamove_stress
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
     &     ,dt,srlx,stbeta,vol,avol,ptnsr,mpi_md_world,cpctl)
c
c  Calc. acceralation of h-matrix used in isobaric MD
c  using given stress tensor stgt.
c
      implicit none
      include "mpif.h"
      integer,intent(in):: natm,mpi_md_world
      real(8),intent(in):: stgt(3,3),eki(3,3,natm),strs(3,3,natm),ptgt
     &     ,sgm(3,3),dt,vol,avol,srlx,stbeta
      real(8),intent(out):: ah(3,3),ptnsr(3,3)
      character(len=*):: cpctl
      
      integer:: i,ixyz,jxyz,kxyz,ierr
      real(8):: pt(3,3),prss,fac

      pt(1:3,1:3)= 0d0
      do i=1,natm
        do jxyz=1,3
          do ixyz=1,3
c            pt(ixyz,jxyz)=pt(ixyz,jxyz)
c     &           +2d0*eki(ixyz,jxyz,i) +strs(ixyz,jxyz,i)*avol
            pt(ixyz,jxyz)=pt(ixyz,jxyz)
     &           +2d0*eki(ixyz,jxyz,i) +strs(ixyz,jxyz,i) !*avol
          enddo
        enddo
      enddo

      ptnsr(1:3,1:3)= 0d0
      call mpi_allreduce(pt,ptnsr,9,mpi_double_precision,mpi_sum
     &     ,mpi_md_world,ierr)
      ptnsr(1:3,1:3)=ptnsr(1:3,1:3) /vol

c.....Berendsen for variable-cell
      if( trim(cpctl).eq.'Berendsen' .or.
     &     trim(cpctl).eq.'vc-Berendsen') then
c.....now ah is scaling factor for h-mat
        ah(1:3,1:3)= 0d0
        ah(1,1)= 1d0
        ah(2,2)= 1d0
        ah(3,3)= 1d0
        ah(1:3,1:3)= ah(1:3,1:3)
     &       -stbeta*dt/3/srlx*( stgt(1:3,1:3)-ptnsr(1:3,1:3) )
c.....Berendsen for variable-volume not variable-cell
      else if( trim(cpctl).eq.'vv-Berendsen' ) then
        ah(1:3,1:3)= 0d0
        ah(1,1)= 1d0
        ah(2,2)= 1d0
        ah(3,3)= 1d0
        prss = (ptnsr(1,1)+ptnsr(2,2)+ptnsr(3,3))/3
        fac = 1.0 -stbeta/dt/3/srlx *(ptgt-prss)
        ah(1:3,1:3) = ah(1:3,1:3)*fac
      endif

      end subroutine force_isobaric
c=======================================================================
      subroutine write_force(ionum,cpostfix,natm,h,tag,aa,acon,nismax
     &     ,epot,myid_md,mpi_md_world,nxyz,ntot,tagtot,atot,ptnsr)
      implicit none
      include 'mpif.h'
      integer,intent(in):: natm,nismax,ionum,nxyz
     &     ,ntot,myid_md,mpi_md_world
      character(len=*),intent(in):: cpostfix
      real(8),intent(in):: h(3,3),acon(nismax),aa(3,natm),tag(natm)
     &     ,epot,ptnsr(3,3)
      real(8),intent(inout):: tagtot(ntot),atot(3,ntot)

      integer:: i,is,n0,ixyz,itag,natmt,ierr
      integer:: istat(mpi_status_size)
      real(8):: at(3)
      integer,parameter:: nmpi = 2

c.....Gather acceralations
c  Since this routine is mainly for fitpot which requires only
c  very small number of atoms, usually one does not need to think about
c  parallelization.
      if( myid_md.eq.0 ) then
        n0 = natm
        atot(1:3,1:natm)= aa(1:3,1:natm)
        n0 = n0 +1
        do ixyz=1,nxyz-1
          itag = ixyz*nmpi -nmpi
          call mpi_recv(natmt,1,mpi_integer,ixyz,itag
     &         ,mpi_md_world,istat,ierr)
          call mpi_recv(tagtot(n0),natmt,mpi_double_precision
     &         ,ixyz,itag+1,mpi_md_world,istat,ierr)
          call mpi_recv(atot(1,n0),3*natmt,mpi_double_precision
     &         ,ixyz,itag+2,mpi_md_world,istat,ierr)
          n0 = n0 + natmt
        enddo
      else
        itag = myid_md*nmpi -nmpi
        call mpi_send(natm,1,mpi_integer,0,itag
     &       ,mpi_md_world,ierr)
        call mpi_send(tag,natm,mpi_double_precision,0,itag+1
     &       ,mpi_md_world,ierr)
        call mpi_send(aa,3*natm,mpi_double_precision,0,itag+2
     &       ,mpi_md_world,ierr)
      endif

      if( myid_md.eq.0 ) then
c.....Write out forces
        open(ionum,file='frc'//trim(cpostfix),status='replace')
        write(ionum,'(i10)') ntot
        do i=1,ntot
          is= int(tagtot(i))
          at(1:3)= h(1:3,1)*atot(1,i)
     &         +h(1:3,2)*atot(2,i)
     &         +h(1:3,3)*atot(3,i)
          at(1:3)= at(1:3)/acon(is)
          write(ionum,'(3es23.14e3)') at(1:3)
        enddo
        close(ionum)

c.....Write out energy
        open(ionum+1,file='erg'//trim(cpostfix),status='replace')
        write(ionum+1,'(es23.14e3)') epot
        close(ionum+1)

c.....Write stress tensor
        open(ionum+2,file='str'//trim(cpostfix),status='replace')
        write(ionum+2,'(6es16.7e3)') ptnsr(1,1),ptnsr(2,2),ptnsr(3,3),
     &       ptnsr(2,3),ptnsr(1,3),ptnsr(1,2)
        close(ionum+2)
      endif

      end subroutine write_force
c=======================================================================
      subroutine setv()
      use variables
      implicit none
      include 'mpif.h'
      include 'params_unit.h'
      integer:: i,l,is
      real(8):: dseed,sumvx,sumvy,sumvz,tmp,facv(nismax)
      real(8),parameter:: pi = 3.14159265358979d0
      real(8),external:: box_muller

      facv(1:nismax)=dsqrt(tinit*fkb*ev2j /(am(1:nismax)*amu2kg))
     &     *m2ang /s2fs

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
        va(1,i)= va(1,i) /h(1,1,0) *dt
        va(2,i)= va(2,i) /h(2,2,0) *dt
        va(3,i)= va(3,i) /h(3,3,0) *dt
      enddo

      end subroutine setv
c=======================================================================
      subroutine rm_trans_motion()
      use variables
      implicit none
      include 'mpif.h'
      integer:: i,is
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

      if( myid_md.eq.0 ) then
        write(6,'(a,3es12.4)') ' sumpx,y,z/amtot='
     &       ,sumpx/amtot,sumpy/amtot,sumpz/amtot
      endif

      end subroutine rm_trans_motion
c=======================================================================
      subroutine check_cmin()
      use variables
      implicit none 

      if( cmin.ne.'' ) then
        if( cmin.eq.'none' ) then
          ifdmp = 0
        else if( cmin.eq.'damp' ) then
          ifdmp = 1
        else if( cmin.eq.'FIRE' ) then
          ifdmp = 2
        else
          write(6,'(a)') ' [Warning] There is no minimization'
     &         //' method: '//cmin
          write(6,'(a)') '           So ifdmp is set 0.'
          ifdmp = 0
        endif
      endif
      
      end subroutine check_cmin
c=======================================================================
      subroutine vfire(num_fire,alp_fire)
      use variables
      implicit none
      include 'mpif.h'
      integer,intent(inout):: num_fire
      real(8),intent(inout):: alp_fire

      integer:: i,ixyz
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
      call mpi_allreduce(vnorml,vnorm,1,mpi_double_precision
     &     ,mpi_sum,mpi_md_world,ierr)
      call mpi_allreduce(fnorml,fnorm,1,mpi_double_precision
     &     ,mpi_sum,mpi_md_world,ierr)
      call mpi_allreduce(fdotvl,fdotv,1,mpi_double_precision
     &     ,mpi_sum,mpi_md_world,ierr)
      vnorm = sqrt(vnorm)
      fnorm = sqrt(fnorm)
      do i=1,natm
        do ixyz=1,3
          va(ixyz,i) = (1d0 -alp_fire)*va(ixyz,i)
     &         +alp_fire*aa(ixyz,i)/fnorm *vnorm
        enddo
      enddo
      
      if( myid_md.eq.0 ) then
        if( fdotv.gt.0d0  ) then
          num_fire = num_fire + 1
          if( num_fire.gt.nmin_fire ) then
            dt = min(dtmax_fire,dt*finc_fire)
            alp_fire = alp_fire*falp_fire
            write(6,'(a,2f10.3,es12.4)') ' dt,alp_fire,fdotv = '
     &           ,dt,alp_fire,fdotv
          endif
        else
          dt = min(dtmax_fire,dt*fdec_fire)
          va(1:3,1:natm) = 0d0
          alp_fire = alp0_fire
          num_fire = 0
          write(6,'(a,2f10.3,es12.4)') ' dt,alp_fire,fdotv = '
     &         ,dt,alp_fire,fdotv
        endif
      endif
      call mpi_bcast(dt,1,mpi_double_precision,0,mpi_md_world,ierr)
      call mpi_bcast(alp_fire,1,mpi_double_precision,0,mpi_md_world
     &     ,ierr)
      end subroutine vfire
c=======================================================================
      subroutine space_decomp()
c
c  Decompose the system and scatter atoms to every process.
c  This would be called only the first time, but the routine is
c  written as if it is called several times.
c
      use variables
      implicit none
      include 'mpif.h'
      integer:: istat(mpi_status_size)

      integer:: i,j,ixyz,n
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
     &             rtot(3,i).lt.szogt+1d0/nz ) then
                n = n+1
              endif
            enddo
            nmax = max(nmax,n)
            nmin = min(nmin,n)
          enddo
          write(6,'(a,i10)') ' min number of local atoms = ',nmin
          write(6,'(a,i10)') ' max number of local atoms = ',nmax
          namax = max(int(nmax*1.2),200)
          nbmax = namax*27
          namax = namax +nbmax
          write(6,'(a,i10)')   ' nbmax = nmax*1.2*27       = ',nbmax
          write(6,'(a,i10)')   ' namax = nmax*1.2 + nbmax  = ',namax
          write(6,'(a,i10,a)') ' memory for main routine   = '
     &         ,namax*(3*4 +3*3*5 +1 +nnmax/2)*8/1000/1000,' MB'
        endif
        call mpi_bcast(namax,1,mpi_integer,0,mpi_md_world,ierr)
        call mpi_bcast(nbmax,1,mpi_integer,0,mpi_md_world,ierr)
c.....allocate memory for local arrays of length namax
        allocate(ra(3,namax),va(3,namax),aa(3,namax),ra0(3,namax)
     &       ,strs(3,3,namax),stt(3,3,namax),tag(namax)
     &       ,lspr(0:nnmax,namax),epi(namax),eki(3,3,namax)
     &       ,stp(3,3,namax),stn(3,3,namax),lsb(0:nbmax,6)
     &       ,lsex(nbmax,6))
      endif

      call mpi_bcast(hunit,1,mpi_double_precision,0,mpi_md_world,ierr)
      call mpi_bcast(h,9,mpi_double_precision,0,mpi_md_world,ierr)

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
     &           rtot(3,i).lt.szogt+1d0/nz ) then
              natm = natm +1
              tag(natm)= tagtot(i)
              ra(1:3,natm)= rtot(1:3,i)
              va(1:3,natm)= vtot(1:3,i)
            endif
          enddo
          if( ixyz.ne.0 ) then
            call mpi_send(natm,1,mpi_integer,ixyz,ixyz
     &           ,mpi_md_world,ierr)
            call mpi_send(tag,natm,mpi_double_precision,ixyz
     &           ,3*(ixyz+1)-3,mpi_md_world,ierr)
            call mpi_send(ra,3*natm,mpi_double_precision,ixyz
     &           ,3*(ixyz+1)-2,mpi_md_world,ierr)
            call mpi_send(va,3*natm,mpi_double_precision,ixyz
     &           ,3*(ixyz+1)-1,mpi_md_world,ierr)
          endif
        enddo
c        write(6,'(a,f10.3)') ' time space_decomp = ',mpi_wtime() -t0
        tspdcmp = mpi_wtime() -t0
      else ! myid_md.ne.0
        call mpi_recv(natm,1,mpi_integer,0,myid_md
     &       ,mpi_md_world,istat,ierr)
        call mpi_recv(tag,natm,mpi_double_precision,0
     &       ,3*(myid_md+1)-3,mpi_md_world,istat,ierr)
        call mpi_recv(ra,3*natm,mpi_double_precision,0
     &       ,3*(myid_md+1)-2,mpi_md_world,istat,ierr)
        call mpi_recv(va,3*natm,mpi_double_precision,0
     &       ,3*(myid_md+1)-1,mpi_md_world,istat,ierr)
      endif
      call mpi_barrier(mpi_md_world,ierr)

      end subroutine space_decomp
c=======================================================================
      subroutine space_comp()
c
c  Opposite to space_decomp, gather atoms from every process
c  to create the total system for output.
c
      use variables
      implicit none
      include 'mpif.h'
      integer,parameter:: nmpi = 7
      integer:: n0,ixyz,natmt,i
      integer:: istat(mpi_status_size),itag
      real(8):: t0

      t0 = mpi_wtime()

      if( myid_md.eq.0 ) then
        n0 = natm
        tagtot(1:natm) = tag(1:natm)
        rtot(1:3,1:natm) = ra(1:3,1:natm)
        vtot(1:3,1:natm) = va(1:3,1:natm)
        epitot(1:natm) = epi(1:natm)
        ekitot(1:3,1:3,1:natm) = eki(1:3,1:3,1:natm)
        stot(1:3,1:3,1:natm) = strs(1:3,1:3,1:natm)
        n0 = n0 +1
        do ixyz=1,nxyz-1
          itag = ixyz*nmpi -nmpi
          call mpi_recv(natmt,1,mpi_integer,ixyz,itag
     &         ,mpi_md_world,istat,ierr)
          call mpi_recv(tagtot(n0),natmt,mpi_double_precision
     &         ,ixyz,itag+1,mpi_md_world,istat,ierr)
          call mpi_recv(rtot(1,n0),3*natmt,mpi_double_precision
     &         ,ixyz,itag+2,mpi_md_world,istat,ierr)
          call mpi_recv(vtot(1,n0),3*natmt,mpi_double_precision
     &         ,ixyz,itag+3,mpi_md_world,istat,ierr)
          call mpi_recv(epitot(n0),natmt,mpi_double_precision
     &         ,ixyz,itag+4,mpi_md_world,istat,ierr)
          call mpi_recv(ekitot(1,1,n0),3*3*natmt,mpi_double_precision
     &         ,ixyz,itag+5,mpi_md_world,istat,ierr)
          call mpi_recv(stot(1,1,n0),3*3*natmt,mpi_double_precision
     &         ,ixyz,itag+6,mpi_md_world,istat,ierr)
          n0 = n0 + natmt
        enddo
        tspdcmp = tspdcmp +(mpi_wtime()-t0)
c        write(6,'(a,f10.3)') ' time space_comp = ',mpi_wtime() -t0
      else ! myid_md.ne.0
        itag = myid_md*nmpi -nmpi
        call mpi_send(natm,1,mpi_integer,0,itag
     &       ,mpi_md_world,ierr)
        call mpi_send(tag,natm,mpi_double_precision,0,itag+1
     &       ,mpi_md_world,ierr)
        do i=1,natm
          ra(1:3,i) = ra(1:3,i) + sorg(1:3)
        enddo
        call mpi_send(ra,3*natm,mpi_double_precision,0,itag+2
     &       ,mpi_md_world,ierr)
        do i=1,natm
          ra(1:3,i) = ra(1:3,i) - sorg(1:3)
        enddo
        call mpi_send(va,3*natm,mpi_double_precision,0,itag+3
     &       ,mpi_md_world,ierr)
        call mpi_send(epi,natm,mpi_double_precision,0,itag+4
     &       ,mpi_md_world,ierr)
        call mpi_send(eki,3*3*natm,mpi_double_precision,0,itag+5
     &       ,mpi_md_world,ierr)
        call mpi_send(strs,3*3*natm,mpi_double_precision,0,itag+6
     &       ,mpi_md_world,ierr)
      endif

      end subroutine space_comp
c=======================================================================
      subroutine sort_by_tag(natm,tag,ra,va,eki,epi,strs)
c
c  Sort by tag for output.
c  
      implicit none
      integer,intent(in):: natm
      real(8),intent(inout):: ra(3,natm),va(3,natm),eki(3,3,natm)
     &     ,epi(natm),strs(3,3,natm),tag(natm)

      integer,allocatable:: itag(:)
      real(8),allocatable:: buf(:,:)
      integer:: i
      integer,save:: nsave = 0
      integer,external:: itotOf

      if( .not. allocated(itag) .or. natm.gt.nsave ) then
        if( allocated(itag) ) deallocate(itag,buf)
        nsave = natm
        allocate(itag(natm),buf(18,natm))
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
        itag(i)= itotOf(tag(i))
      enddo

      call heapsort_itag(natm,natm,itag,18,buf)

      do i=1,natm
        ra(1:3,i)= buf(1:3,i)
        va(1:3,i)= buf(4:6,i)
        eki(1,1,i)= buf(7,i)
        epi(i)= buf(8,i)
        strs(1:3,1,i)= buf(9:11,i)
        strs(1:3,2,i)= buf(12:14,i)
        strs(1:3,3,i)= buf(15:17,i)
        tag(i)= buf(18,i)
      enddo

      end subroutine sort_by_tag
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make pmd"
c     End:
