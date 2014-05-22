      program parallel_md
c-----------------------------------------------------------------------
c Spatial decomposition parallel molecular dynamics program.
c
c-----------------------------------------------------------------------
c INPUT FILES:
c ------------
c   in.pmd:     Main input file
c   pmd###-000: MD-coordination files, ### is the node number
c
c OUTPUT FILES:
c -------------
c   out.erg:    Total, kinetic, and potential energy
c   out.stress: Stress component normal to z-upper surface of nanorod,
c               and z-strain of the nanorod.
c   pmd###-???: MD-coordination file of each MD-step,
c               ### is the node number and ??? is the frame number
c-----------------------------------------------------------------------
      use variables
#ifdef __WALL__
      use wall
#endif
      implicit none
      include "mpif.h"
      include "./params_unit.h"

#ifdef __DISL__
c.....Epot threshold for disl core extraction [Hartree]
      real(8),parameter:: epith = -0.1410d0
#endif

      integer:: i,j,k,l,m,n,ia,ib,is
      real(8):: tmp
      integer,external:: num_dof
c-----output file names
      character:: cnode*5,cnum*4

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
        call write_initial_setting()
c.....Multiply damping factor to avoid multiplication of every step
        if( ifdmp.eq.1 ) then
          fmv(1:3,0:9)= fmv(1:3,0:9) *dmp
        endif
      endif

c-----Broadcast input parameters to all nodes
      call mpi_bcast(nx,1,MPI_INTEGER,0,mpi_md_world,ierr)
      call mpi_bcast(ny,1,MPI_INTEGER,0,mpi_md_world,ierr)
      call mpi_bcast(nz,1,MPI_INTEGER,0,mpi_md_world,ierr)
      call mpi_bcast(nstp,1,MPI_INTEGER,0,mpi_md_world,ierr)
      call mpi_bcast(dt,1,MPI_DOUBLE_PRECISION,0,mpi_md_world,ierr)
      call mpi_bcast(rc,1,MPI_DOUBLE_PRECISION,0,mpi_md_world,ierr)
      call mpi_bcast(rbuf,1,MPI_DOUBLE_PRECISION,0,mpi_md_world,ierr)
      call mpi_bcast(dmp,1,MPI_DOUBLE_PRECISION,0,mpi_md_world,ierr)
      call mpi_bcast(tinit,1,mpi_double_precision,0,mpi_md_world,ierr)
      call mpi_bcast(iftctl,1,MPI_INTEGER,0,mpi_md_world,ierr)
      call mpi_bcast(treq,1,MPI_DOUBLE_PRECISION,0,mpi_md_world,ierr)
      call mpi_bcast(trlx,1,MPI_DOUBLE_PRECISION,0,mpi_md_world,ierr)
      call mpi_bcast(nerg,1,MPI_INTEGER,0,mpi_md_world,ierr)
      call mpi_bcast(ifpmd,1,MPI_INTEGER,0,mpi_md_world,ierr)
      call mpi_bcast(npmd,1,MPI_INTEGER,0,mpi_md_world,ierr)
      call mpi_bcast(ifdmp,1,MPI_INTEGER,0,mpi_md_world,ierr)
      call mpi_bcast(fmv,30,mpi_double_precision,0,mpi_md_world,ierr)
      call mpi_bcast(shrst,1,mpi_double_precision,0,mpi_md_world,ierr)
      call mpi_bcast(ifpctl,1,MPI_INTEGER,0,mpi_md_world,ierr)
      call mpi_bcast(ptgt,1,mpi_double_precision,0,mpi_md_world,ierr)
      call mpi_bcast(vmcoeff,1,mpi_double_precision,0,mpi_md_world,ierr)
      call mpi_bcast(voldmp,1,mpi_double_precision,0,mpi_md_world,ierr)
      call mpi_bcast(strfin,1,mpi_double_precision,0,mpi_md_world,ierr)
      call mpi_bcast(am,nismax,mpi_double_precision,0,mpi_md_world,ierr)
      call mpi_bcast(ciofmt,6,mpi_character,0,mpi_md_world,ierr)
      call mpi_bcast(cforce,20,mpi_character,0,mpi_md_world,ierr)
#ifdef __WALL__
      call mpi_bcast(wtop,1,mpi_double_precision,0,mpi_md_world,ierr)
      call mpi_bcast(wbot,1,mpi_double_precision,0,mpi_md_world,ierr)
      call mpi_bcast(ptgt_wall,1,mpi_double_precision,0,mpi_md_world
     &     ,ierr)
      call mpi_bcast(trlx_wall,1,mpi_double_precision,0,mpi_md_world
     &     ,ierr)
      call mpi_bcast(nout_wall,1,mpi_integer,0,mpi_md_world,ierr)
#endif

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
        call mpi_finalize(ierr)
        stop
      endif
#ifdef __WALL__
      trlx_wall= trlx_wall
      dt_wall= dt*nout_wall
c.....GPa --> eV/AA^3
      ptgt_wall= ptgt_wall *gpa2up
#endif
      
c-----vector node indices: range [0:nx-1]
      myx=myid_md/(ny*nz)
      myy=mod(myid_md/nz,ny)
      myz=mod(myid_md,nz)
c-----reduced node origin
      sorg(1)= anxi*myx
      sorg(2)= anyi*myy
      sorg(3)= anzi*myz

c.....Each node opens one MD-configuration file
      write(cnode(1:5),'(i5.5)') myid_md
      if( trim(ciofmt).eq.'bin' .or. trim(ciofmt).eq.'binary' ) then
        call read_pmd_bin(20,"0000/pmd"//cnode
     &       ,namax,natm,h,tag,ra,va,eki,epi,strs)
      else if( trim(ciofmt).eq.'ascii' ) then
        call read_pmd_ascii(20,"0000/pmd"//cnode
     &     ,namax,natm,h,tag,ra,va,eki,epi,strs)
      endif
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
      write(6,'(a,2i8)') " myid,natm=",myid_md,natm


c-----setup
      call setup(nismax,am,acon,fack,dt)
c-----set HI and SGM
      call boxmat(h,hi,ht,g,gi,gt,vol,sgm)
c-----ntset
      call ntset(myx,myy,myz,nx,ny,nz,nn,sv,myparity,anxi,anyi,anzi)

c-----get total number of particles
      call mpi_allreduce(natm,ntot,1,MPI_INTEGER,MPI_SUM
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

#ifdef __SHEAR__
      call init_shear_stress(shrst,shrfx,natm,tag,ra,h,sorg,acon(1)
     &     ,myid_md,mpi_md_world)
#endif

      tcpu1= mpi_wtime()
      tcom = 0d0

c-----copy RA of boundary atoms
      call bacopy(tcom,sgm,vol,lsb,nbmax,ra,namax,natm,nb
     &     ,anxi,anyi,anzi,nn,tag,rc,myid_md,myparity,lsrc,sv
     &     ,mpi_md_world)
c-----Make pair list
      call mk_lspr(namax,natm,nbmax,nb,nnmax,tag,ra,rc+rbuf,h,hi
     &     ,anxi,anyi,anzi,lspr)
c.....Calc forces
      call get_force(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &     ,tcom,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
     &     ,mpi_md_world,myid_md,epi,epot0,nismax,acon,avol
     &     ,cforce)
      epot= epot0

#ifdef __SHEAR__
      call add_shear_stress(natm,ra,aa,tag,sorg,shrfx)
#endif

#ifdef __DISL__
      call perf_disl_pos_by_pot(epith,natm,ra,h,epi,sorg
     &     ,nodes_md,myid_md,mpi_md_world,0,21)
#endif

c-----calc kinetic energy
      call get_ekin(namax,natm,va,tag,h,nismax,fack,ekin,eki
     &     ,vmax,mpi_md_world)
      vmaxold=vmax
c-----in case of 3D, <ekin> = 3/2 *NkT
      temp= ekin /ntot/fkb *2d0/3d0
cc-----in case of 2D, <ekin> = 2/2 *NkT
c      temp= ekin /ntot/fkb *2d0/2d0

      if( ifpctl.eq.1 ) then
        if(myid_md.eq.0) write(6,'(a,f10.2,a,es12.4,a)')
     &       ' Target Pressure = ',ptgt,' GPa, ',ptgt*gpa2up,' a.u.'
        ptgt= ptgt *gpa2up
c.....Parrinello-Rahman Isobaric
        vm= am(1)*vmcoeff
        call force_isobaric(ptgt,ah,natm,eki,strs,sgm,vm,vol,avol
     &       ,ptnsr,mpi_md_world)
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
        write(6,'(1x,a)') "initial values:"
        write(6,'(1x,a,es22.14,a)') " kinetic energy  =",ekin,' eV'
        write(6,'(1x,a,es22.14,a)') " potential energy=",epot0,' eV'
        write(6,'(1x,a,es15.7,a)') " temperature     =",temp,' K'

        if( ifpctl .eq. 1 ) then
          write(6,'(a,i8,2f10.2,2es11.3)')
     &         " istp,time,temp,vol,prss="
     &         ,istp,tcpu,temp,vol
     &         ,(ptnsr(1,1)+ptnsr(2,2)+ptnsr(3,3))/3*up2gpa
        else
          write(6,'(a,i8,2f10.2)')
     &         " istp,time,temp="
     &         ,istp,tcpu,temp
        endif
      endif

c-----initialize the counter for output
      iocntpmd=0
      iocnterg=0
      write(cnum(1:4),'(i4.4)') iocntpmd
      if(myid_md.eq.0) call write_force(21,cnum,natm,h,tag,aa,acon
     &     ,nismax,epot0)
c-----output initial configuration
      if( ifpmd.eq.1 ) then
        call system("mkdir -p "//cnum)
        if( trim(ciofmt).eq.'bin' .or. trim(ciofmt).eq.'binary' ) then
          call write_pmd_bin(20,cnum//"/pmd"//cnode
     &         ,natm,h,tag,ra,va,eki,epi,strs,sorg,dt)
        elseif( trim(ciofmt).eq.'ascii' ) then
          call write_pmd_ascii(20,cnum//"/pmd"//cnode
     &         ,natm,h,tag,ra,va,eki,epi,strs,sorg,dt)
        endif
      endif

#ifdef __WALL__
      call update_wall()
#endif

c-----Energy output
      if(myid_md.eq.0) then
        open(11,file="out.erg",status='replace')
        write(11,'(a)') '# istp,etot-epot0[eV],ekin,epot-epot0,temp[K]'
     &       //'vol[Ang^3],pressure[GPa]'
        write(11,'(a,es15.7,a)') '# epot0 =',epot0,' [eV]'
        write(11,'(i8,3es15.7,f10.2,2es15.7)') istp
     &       ,ekin+epot0-epot0,ekin,epot0-epot0,temp
     &       ,vol
     &       ,(ptnsr(1,1)+ptnsr(2,2)+ptnsr(3,3))/3*up2gpa
        call flush(11)
#ifdef __KOUKYU__
        open(80,file='out.stress',status='replace')
        write(80,'(a)') '#  istp,   strnow,      ftop,     fbot'
#endif
#ifdef __WALL__
        open(81,file='out.wall',status='replace')
        write(81,'(a)') '#   istp,nitop,nibot,ptop(GPa),pbot(GPa)'
     &       //',wtop,wbot'
        write(81,'(3i8,4es15.7)')
     &       istp,nitop,nibot,ptop*aup2gpa,pbot*aup2gpa,wtop,wbot
#endif
      endif

c.....Set al
      al(1)= h(1,1,0)
      al(2)= h(2,2,0)
      al(3)= h(3,3,0)

c-----MD loop with velocity-Verlet--------------------------------------
      do istp=1,nstp

c-------first kick of velocities
        va(1:3,1:natm)=va(1:3,1:natm) +aa(1:3,1:natm)

c-------multiply fmv or damping
        do i=1,natm
          l= int(mod(tag(i)*10,10d0))
          va(1:3,i)=va(1:3,i) *fmv(1:3,l)
        enddo

c.....Parrinello-Rahman isobaric
        if( ifpctl.eq.1 ) then
          h(1:3,1:3,1)=h(1:3,1:3,1) +ah(1:3,1:3)*dt/2
          h(1:3,1:3,1)=h(1:3,1:3,1)*voldmp
c.....Only orthogonal change
          do i=1,3
            do j=1,3
              if( j.eq.i ) cycle
              h(j,i,1)= 0d0
            enddo
          enddo
        endif

c-------update positions
        ra(1:3,1:natm)=ra(1:3,1:natm) +va(1:3,1:natm)

#ifdef __KOUKYU__
        call tensile_loading(natm,ra,tag,nstp,strfin,strnow
     &       ,sorg,mpi_md_world)
#endif

        if( ifpctl.eq.1 ) then
          h(1:3,1:3,0)=h(1:3,1:3,0) +h(1:3,1:3,1)*dt
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
          call bacopy(tcom,sgm,vol,lsb,nbmax,ra,namax,natm,nb
     &         ,anxi,anyi,anzi,nn,tag,rc,myid_md,myparity,lsrc,sv
     &         ,mpi_md_world)
c.....Make pair list
          call mk_lspr(namax,natm,nbmax,nb,nnmax,tag,ra,rc+rbuf,h,hi
     &         ,anxi,anyi,anzi,lspr)
          ts= istp*dt
          te=rbuf/(2d0*vmax/dt) +ts
          istpe= te/dt
c          if(myid_md.eq.0) write(6,'(a,i10)')
c     &         " Next update step will be ",istpe
        else
c.....Copy RA of boundary atoms determined by 'bacopy'
          call bacopy_fixed(tcom,sgm,vol,lsb,nbmax,ra,namax,natm,nb
     &         ,anxi,anyi,anzi,nn,tag,rc,myid_md,myparity,lsrc,sv
     &         ,mpi_md_world)
        endif

c-------Calc forces
        call get_force(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &       ,tcom,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
     &       ,mpi_md_world,myid_md,epi,epot,nismax,acon,avol
     &       ,cforce)

#ifdef __SHEAR__
        call add_shear_stress(natm,ra,aa,tag,sorg,shrfx)
#endif

c.....Second kick of velocities
        va(1:3,1:natm)=va(1:3,1:natm) +aa(1:3,1:natm)

c.....Calc kinetic energy
        vmaxold= vmax
        call get_ekin(namax,natm,va,tag,h,nismax,fack,ekin,eki
     &       ,vmax,mpi_md_world)

c        write(6,'(a,10es12.4)') '251,va(1:3),eki(251),epi(251)='
c     &       ,va(1:3,251),eki(1,1,251)+eki(2,2,251)+eki(3,3,251)
c     &       ,epi(251)

        if( ifpctl.eq.1 ) then
          call force_isobaric(ptgt,ah,natm,eki,strs,sgm
     &         ,vm,vol,avol,ptnsr,mpi_md_world)
          h(1:3,1:3,1)=h(1:3,1:3,1) +ah(1:3,1:3)*dt/2
          gg(1:3,1:3)= 0d0
          do j=1,3
            do i=1,3
              do k=1,3
                gg(i,j)=gg(i,j) +gi(i,k)*g(k,j,1)
              enddo
            enddo
          enddo
          do i=1,natm
            do j=1,3
              do k=1,3
                aa(j,i)=aa(j,i) -gg(j,k)*va(k,i)
              enddo
            enddo
          enddo
        endif

c-------temperature control by velocity scaling
        if(iftctl.eq.1) then
          ndof= num_dof()
          if( ndof.ne.0 ) then
!            write(6,'(a,2i6)') ' myid,ndof=',myid_md,ndof
c.....for the case that there are fixed atoms..., <ekin>= 1/2 *ndof*kT
            temp= ekin *2d0 /fkb /ndof
c---------in case of 3D, <ekin> = 3/2 *NkT
c          temp= ekin /ntot/fkb *2d0/3d0
cc---------in case of 2D, <ekin> = 2/2 *NkT
c          temp= ekin /ntot/fkb *2d0/2d0
c---------Berendsen thermostat
            if( (treq-temp)/temp.gt.100d0 ) then
              tmp= sqrt(1d0 +dt/trlx*100d0 )
            else
              tmp= sqrt(1d0 +dt/trlx*(treq-temp)/temp)
            endif
            va(1:3,1:natm)= va(1:3,1:natm) *tmp
          endif
        endif

cc.....Eliminate translational motion
c        if( mod(istp,100).eq.0 ) then
c          call rm_trans_motion()
c        endif

c-------output energies
        if(mod(istp,nouterg).eq.0) then
          iocnterg=iocnterg+1
cc---------calc kinetic energy
c          if(iftctl.ne.1) then
c            vmaxold= vmax
c            call get_ekin(namax,natm,va,tag,h,nismax,fack,ekin,eki
c     &           ,vmax,mpi_md_world)
c          endif
          if( iftctl.ne.1 ) then
c---------in case of 3D, <ekin> = 3/2 *NkT
            temp= ekin /ntot/fkb *2d0/3d0
cc---------in case of 2D, <ekin> = 2/2 *NkT
c          temp= ekin /ntot/fkb *2d0/2d0
          endif
#ifdef __KOUKYU__
          call get_stress_on_base(natm,ra,aa,tag,h,acon(1),ftop,fbot
     &         ,sorg,mpi_md_world)
#endif
          if(myid_md.eq.0) then
            write(11,'(i8,3es15.7,f10.2,2es15.7)') istp
     &           ,ekin+epot-epot0,ekin,epot-epot0,temp
     &           ,vol
     &           ,(ptnsr(1,1)+ptnsr(2,2)+ptnsr(3,3))/3*up2gpa
            call flush(11)
#ifdef __KOUKYU__
            write(80,'(i8,3es15.7)') istp,strnow,ftop,fbot
            call flush(80)
#endif
          endif
c---------output step, time, and temperature
          tcpu= mpi_wtime() -tcpu1
          if(myid_md.eq.0) then
            if( ifpctl .eq. 1 ) then
              write(6,'(a,i8,2f10.2,2es11.3)')
     &             " istp,time,temp,vol,prss="
     &             ,istp,tcpu,temp,vol
     &             ,(ptnsr(1,1)+ptnsr(2,2)+ptnsr(3,3))/3*up2gpa
            else
              write(6,'(a,i8,2f10.2)')
     &             " istp,time,temp="
     &             ,istp,tcpu,temp
            endif
          endif
            
#ifdef __DISL__
c.....Output disl core pos
          call perf_disl_pos_by_pot(epith,natm,ra,h,epi,sorg
     &         ,nodes_md,myid_md,mpi_md_world,iocnterg,21)
#endif
        endif

c-------write the particle positions
        if(ifpmd.eq.1.and.mod(istp,noutpmd).eq.0)then
c---------decide pmd-file name
          iocntpmd=iocntpmd+1
          write(cnum(1:4),'(i4.4)') iocntpmd
          call system("mkdir -p "//cnum)
c---------write pmd file for Akira
          if( trim(ciofmt).eq.'bin' .or. trim(ciofmt).eq.'binary' ) then
            call write_pmd_bin(20,cnum//"/pmd"//cnode
     &           ,natm,h,tag,ra,va,eki,epi,strs,sorg,dt)
          elseif( trim(ciofmt).eq.'ascii' ) then
            call write_pmd_ascii(20,cnum//"/pmd"//cnode
     &           ,natm,h,tag,ra,va,eki,epi,strs,sorg,dt)
          endif
        endif

#ifdef __WALL__
        call wall_reflect()
        if( mod(istp,nout_wall).eq.0 ) then
          call update_wall()
          if( myid_md.eq.0 ) then
            write(81,'(3i8,4es15.7)')
     &           istp,nitop,nibot,ptop*aup2gpa,pbot*aup2gpa,wtop,wbot
          endif
c.....reset pressure
          ptop= 0d0
          pbot= 0d0
          nitop= 0
          nibot= 0
        endif
        wtop=wtop +dwtop
        wbot=wbot +dwbot
        wtop= min(wtop,0.99d0)
        wtop= max(wtop,btop)
        wbot= max(wbot,0.01d0)
        wbot= min(wbot,bbot)
#endif

      enddo

      tcpu2= mpi_wtime()

      if(myid_md.eq.0) then
        close(11)
#ifdef __KOUKYU__
        close(80)
#endif
#ifdef __WALL__
        close(81)
#endif
      endif

      tcpu= tcpu2 -tcpu1
      if(myid_md.eq.0) then
        write(6,'(1x,a)') "final values:"
        write(6,'(1x,a,es22.14,a)') " kinetic energy  =",ekin,' eV'
        write(6,'(1x,a,es22.14,a)') " potential energy=",epot,' eV'
        write(6,'(1x,a,es15.7,a)') " temperature     =",temp,' K'
        write(6,'(1x,a,2f10.2)') "time cpu & com=",tcpu,tcom
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
      write(6,'(2x,a25,i8)') 'flag_damping',ifdmp
      write(6,'(2x,a25,es11.3)') 'damping_coeff',dmp
      write(6,'(2x,a)') ''
      write(6,'(2x,a25,f12.4)') 'initial_temperature',tinit
      write(6,'(2x,a25,i8)') 'flag_temperature',iftctl
      write(6,'(2x,a25,f12.4)') 'temperature_target',treq
      write(6,'(2x,a25,f12.4)') 'temperature_relax_time',trlx
      write(6,'(2x,a)') ''
      write(6,'(2x,a25,i8)') 'flag_isobaric',ifpctl
      write(6,'(2x,a25,es11.3)') 'pressure_target',ptgt
      write(6,'(2x,a25,es11.3)') 'vol_mass_coeff',vmcoeff
      write(6,'(2x,a25,f12.4)') 'vol_change_damping',voldmp
      write(6,'(2x,a)') ''
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
      write(6,'(a,es12.4,a,es12.4,a)') ' umass*Ang^2/fs^2 ='
     &     ,tmp,' Joule =',tmp*j2ev,' eV'

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
      subroutine get_ekin(namax,natm,va,tag,h,nismax,fack,ekin,eki
     &     ,vmax,mpi_md_world)
      implicit none 
      include "mpif.h"
      integer,intent(in):: namax,natm,mpi_md_world,nismax
      real(8),intent(in):: va(3,namax),h(3,3),fack(nismax)
     &     ,tag(namax)
      real(8),intent(out):: ekin,eki(3,3,namax),vmax
c-----locals
      integer:: i,ierr,is,ixyz,jxyz,imax
      real(8):: ekinl,x,y,z,v(3),v2,vmaxl

      ekinl=0d0
      eki(1:3,1:3,1:natm)= 0d0
      vmaxl= 0d0

      do i=1,natm
        is= int(tag(i))
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
c.....Find max speed
        if( v2.gt.vmaxl ) imax=i
        vmaxl=max(vmaxl,v2)
      enddo

      ekin= 0d0
      call mpi_allreduce(ekinl,ekin,1,MPI_DOUBLE_PRECISION
     &     ,MPI_SUM,mpi_md_world,ierr)
      vmax= 0d0
      call mpi_allreduce(vmaxl,vmax,1,mpi_double_precision
     &     ,MPI_MAX,mpi_md_world,ierr)
      vmax=sqrt(vmax)

      end subroutine get_ekin
c=======================================================================
      function num_dof()
      use variables
      implicit none
      include 'mpif.h'
      integer:: num_dof
      integer:: i,l,k,ndofl
      real(8),parameter:: deps= 1d-12
      
      ndofl= 0
      do i=1,natm
        l= int(mod(tag(i)*10,10d0))
        do k=1,3
          if( abs(fmv(k,l)).lt.1d-12 ) cycle
          ndofl= ndofl +1
        enddo
      enddo
      
      num_dof= 0
      call mpi_allreduce(ndofl,num_dof,1,mpi_integer,mpi_sum
     &     ,mpi_md_world,ierr)
      return
      end function num_dof
c=======================================================================
      subroutine bacopy(tcom,sgm,vol,lsb,nbmax,ra,namax,natm,nb
     &     ,anxi,anyi,anzi,nn,tag,rc,myid_md,myparity,lsrc,sv
     &     ,mpi_md_world)
c-----------------------------------------------------------------------
c     Exchanges boundary-atom data among neighbor nodes: tag and ra
c-----------------------------------------------------------------------
      implicit none
      include 'mpif.h'
      integer,intent(in):: namax,nbmax,nn(6)
     &     ,myid_md,myparity(3),mpi_md_world
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
c      write(6,'(a,3es12.4)') "rca,rcb,rcc=",rca,rcb,rcc

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
            if (myid_md.eq.0) then
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
          if (myid_md.eq.0) then
            write(*,*)'NAMAX overflowed at bacopy'
            write(*,*)'N+NB NAMAX = ',maxna,namax
          endif
          call MPI_FINALIZE(ierr)
          stop
        endif

100   continue

c-----num. of received boundary atoms
      nb=nbnew

      end subroutine bacopy
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
      subroutine copy_dba_bk(tcom,namax,natm,nbmax,nb,lsb
     &     ,lsrc,myparity,nn,mpi_md_world,x,ndim)
c-----------------------------------------------------------------------
c     Send-back & receive reaction on cached-atoms
c-----------------------------------------------------------------------
      implicit none
      include "mpif.h"
      integer,intent(in):: namax,natm,nbmax,nb,mpi_md_world,ndim
      integer,intent(in):: lsb(0:nbmax,6),lsrc(6),myparity(3),nn(6)
      real(8),intent(inout):: x(ndim,namax),tcom

      integer:: status(MPI_STATUS_SIZE)
      integer:: i,j,k,l,m,n,kd,kdd,ku,kuc,ibkwd,nsd,nsd3,nrc,nrc3,nsdbk
     &     ,ierr,natmx
      real(8):: tcom1,tcom2
      real(8),save,allocatable:: dbuf(:,:),dbufr(:,:)
      logical,save:: l1st=.true.

c      if( l1st ) then
        allocate(dbuf(ndim,nbmax),dbufr(ndim,nbmax))
c        l1st=.false.
c      endif

c-----natmx
      natmx= natm +nb

c-----num. of sent-back reactions
      nsdbk= 0

c-----send-back reactions in the reverse orer, z, y, & x
      do kd=3,1,-1

c-------To calculate the communication time
        tcom1=MPI_WTIME()

c-------higher & lower directions
        do kdd=0,-1,-1
          ku= 2*kd +kdd
          if(mod(ku,2).eq.0) then
            kuc= ku-1
          else
            kuc= ku+1
          endif
          ibkwd= nn(kuc)
c---------num. of to-be-sent particles
          nsd= lsrc(ku)
          nsd3= ndim*nsd
c---------num. of to-be-recieved particles
          nrc= lsb(0,ku)
          nrc3= ndim*nrc

c---------to-be-sent-back particles
          do i=1,nsd
            dbuf(1:ndim,i)= x(1:ndim,natmx-nsdbk-nsd+i)
          enddo
          call mespasd(ibkwd,myparity(kd),dbuf,dbufr,nsd3,nrc3,500
     &         ,mpi_md_world)
          do k=1,nrc
            i=lsb(k,ku)
            x(1:ndim,i)= x(1:ndim,i) +dbufr(1:ndim,k)
          enddo

c---------mpi barrier
          call mpi_barrier(mpi_md_world,ierr)
c---------accumulate num. of already sent-back-particles
          nsdbk=nsdbk +nsd
        enddo

c-------Add the communication time to COMT
        tcom2=MPI_WTIME()
        tcom=tcom+tcom2-tcom1

      enddo

c-----check
      if(nsdbk.ne.nb) then
        write(6,'(a,2i8)') "nsdbk,nb=",nsdbk,nb
        stop "error: nsdbk.ne.nb!!"
      endif

      deallocate(dbuf,dbufr)
      end subroutine copy_dba_bk
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
      subroutine mk_lspr(namax,natm,nbmax,nb,nnmax,tag,ra,rc
     &     ,h,hi,anxi,anyi,anzi,lspr)
      implicit none
      integer,intent(in):: namax,natm,nbmax,nb,nnmax
      integer,intent(out):: lspr(0:nnmax,natm)
      real(8),intent(in):: ra(3,namax),rc,anxi,anyi,anzi
     &     ,hi(3,3),h(3,3),tag(namax)

      integer:: i,j,k,l,m,n
      integer:: mx,my,mz,kux,kuy,kuz,m1x,m1y,m1z,m1,ic,jc,ierr
      real(8):: xi(3),xij(3),rij(3),rij2
      
      logical,save:: l1st=.true.
      integer,allocatable,save:: lscl(:),lshd(:)
      real(8),save:: rc2,rcx,rcy,rcz,rcxi,rcyi,rczi
      integer,save:: lcx,lcy,lcz,lcx2,lcy2,lcz2,lcyz2,lcxyz2

      if( l1st ) then
        rc2= rc**2
c-----make a linked cell list, LSCL
        lcx=anxi/dsqrt(hi(1,1)**2+hi(1,2)**2+hi(1,3)**2)/rc
        lcy=anyi/dsqrt(hi(2,1)**2+hi(2,2)**2+hi(2,3)**2)/rc
        lcz=anzi/dsqrt(hi(3,1)**2+hi(3,2)**2+hi(3,3)**2)/rc
c.....In case that system is thinner than rc, modify lc?.
c.....but notice this modification does not correct results.
        if( lcx.eq.0 ) lcx=1
        if( lcy.eq.0 ) lcy=1
        if( lcz.eq.0 ) lcz=1
        lcx2= lcx +2
        lcy2= lcy +2
        lcz2= lcz +2
        lcyz2=lcy2*lcz2
        lcxyz2=lcx2*lcyz2
        rcx= anxi/lcx
        rcy= anyi/lcy
        rcz= anzi/lcz
        rcxi=1d0/rcx
        rcyi=1d0/rcy
        rczi=1d0/rcz
c        write(6,'(a,3i8)') ' lcx,lcy,lcz=',lcx,lcy,lcz
c        write(6,'(a,3i8)') ' lcx2,lcy2,lcz2=',lcx2,lcy2,lcz2
c        write(6,'(a,i8)') ' lcxyz2=',lcxyz2
c        write(6,'(a,3es12.4)') ' rcx,rcy,rcz=',rcx,rcy,rcz
c-----allocate LSCL & LSHD after obtaining lcxyz2
        allocate(lscl(namax+nbmax),lshd(lcxyz2))
        l1st=.false.
      endif

c-----reset pair list, LSPR
      lspr(0,1:natm)= 0

c-----reset headers
      lshd(1:lcxyz2)= 0


c-----construct a linked-cell list, LSCL, & a header list, LSHD
      do i=1,natm+nb
c-------assign a vector cell index
        mx=(ra(1,i)+rcx)*rcxi
        my=(ra(2,i)+rcy)*rcyi
        mz=(ra(3,i)+rcz)*rczi
c-------classify residents in inner cells even if they are not
        if(i.le.natm) then
          mx= min(max(mx,1),lcx)
          my= min(max(my,1),lcy)
          mz= min(max(mz,1),lcz)
c-------copied atoms are either in inner or surface cells
        else
          mx= min(max(mx,0),lcx+1)
          my= min(max(my,0),lcy+1)
          mz= min(max(mz,0),lcz+1)
        endif
        m= mx*lcyz2 +my*lcz2 +mz +1
        lscl(i)= lshd(m)
c-------the last one goes to the header
        lshd(m)= i
      enddo
c      write(6,'(a)') ' lscl,lshd done'

c-----make a pair list, LSPR
c-----Scan resident cells
      do 100 mz=1,lcz
      do 100 my=1,lcy
      do 100 mx=1,lcx
        m= mx*lcyz2 +my*lcz2 +mz +1
        if (lshd(m).eq.0) goto 100
        do 102 kuz= -1,1
        do 102 kuy= -1,1
        do 102 kux= -1,1
          m1x= mx +kux
          m1y= my +kuy
          m1z= mz +kuz
          m1=m1x*lcyz2 +m1y*lcz2 +m1z +1
          if (lshd(m1).eq.0) goto 102

          i=lshd(m)
1         continue
          if (natm.lt.i) goto 4

c          ic=is(i)
          ic= int(tag(i))
          xi(1:3)= ra(1:3,i)

          j=lshd(m1)

2         continue
          if (j.eq.i) goto 3
c          jc=is(j)
          jc= int(tag(i))
          xij(1:3)= ra(1:3,j) -xi(1:3)
          rij(1)= h(1,1)*xij(1) +h(1,2)*xij(2) +h(1,3)*xij(3)
          rij(2)= h(2,1)*xij(1) +h(2,2)*xij(2) +h(2,3)*xij(3)
          rij(3)= h(3,1)*xij(1) +h(3,2)*xij(2) +h(3,3)*xij(3)
          rij2= rij(1)**2 +rij(2)**2 +rij(3)**2

c          if (rij2.lt.rcij2(ic,jc)) then
          if(rij2.lt.rc2) then
            lspr(0,i)=lspr(0,i)+1
            if(lspr(0,i).gt.nnmax) then
              write(6,'(a)') " !!!lspr(0,i).gt.nnmax"
              call mpi_finalize(ierr)
              stop
            endif
            lspr(lspr(0,i),i)=j
          endif

c---------Continue until j= 0
3         j=lscl(j)
          if (j.gt.0) goto 2

c---------Continue until i= 0
4         i=lscl(i)
          if (i.gt.0) goto 1

102     continue
100   continue

      end subroutine mk_lspr
c=======================================================================
      subroutine mespasi(inode,parity,ibufs,ibufr,nsd,nrc,tag,
     &     mpi_md_world)
c-----------------------------------------------------------------------
c     Integer message passing.  
c-----------------------------------------------------------------------
      include 'mpif.h'
      integer,intent(in):: inode,parity,nsd,nrc,tag
      integer,intent(in):: ibufs(nsd)
      integer,intent(out):: ibufr(nrc)
c-----locals
      integer status(MPI_STATUS_SIZE),ierr

c-----Even: send & recv
      if (parity.eq.0) then
        call MPI_SEND(ibufs,nsd,MPI_INTEGER,inode,tag,
     &       mpi_md_world,ierr)
        call MPI_RECV(ibufr,nrc,MPI_INTEGER,MPI_ANY_SOURCE,tag,
     &       mpi_md_world,status,ierr)
c-----Odd: recv & send
      else if (parity.eq.1) then
        call MPI_RECV(ibufr,nrc,MPI_INTEGER,MPI_ANY_SOURCE,tag,
     &       mpi_md_world,status,ierr)
        call MPI_SEND(ibufs,nsd,MPI_INTEGER,inode,tag,
     &       mpi_md_world,ierr)
c-----Exchange information with myself
      else
        do i=1,nrc
          ibufr(i)=ibufs(i)
        enddo
      endif
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mespasd(inode,parity,bufs,bufr,nsd,nrc,tag,
     & mpi_md_world)
c-----------------------------------------------------------------------
c     Real*8 message passing.
c-----------------------------------------------------------------------
      include 'mpif.h'
      integer,intent(in):: inode,parity,nsd,nrc,tag
      real(8),intent(in):: bufs(nsd)
      real(8),intent(out):: bufr(nrc)
      integer:: status(MPI_STATUS_SIZE),ierr

c-----Even: send & recv
      if (parity.eq.0) then
        call MPI_SEND(bufs,nsd,MPI_DOUBLE_PRECISION,inode,tag,
     &       mpi_md_world,ierr) 
        call MPI_RECV(bufr,nrc,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,tag,
     &       mpi_md_world,status,ierr) 
c-----Odd: recv & send
      else if (parity.eq.1) then
        call MPI_RECV(bufr,nrc,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,tag,
     &       mpi_md_world,status,ierr) 
        call MPI_SEND(bufs,nsd,MPI_DOUBLE_PRECISION,inode,tag,
     &       mpi_md_world,ierr)
c-----Exchange information with myself
      else
        do i=1,nrc
          bufr(i)=bufs(i)
        enddo
      endif
      return
      end
c=======================================================================
      subroutine force_isobaric(ptgt,ah,natm,eki,strs,sgm
     &     ,vm,vol,avol,ptnsr,mpi_md_world)
c
c  Calc. acceralation of h-matrix used in isobaric MD
c
      implicit none
      include "mpif.h"
      integer,intent(in):: natm,mpi_md_world
      real(8),intent(in):: ptgt,eki(3,3,natm),strs(3,3,natm)
     &     ,sgm(3,3),vm,vol,avol
      real(8),intent(out):: ah(3,3),ptnsr(3,3)
      
      integer:: i,ixyz,jxyz,kxyz,ierr
      real(8):: pt(3,3)

      pt(1:3,1:3)= 0d0
      do i=1,natm
        do jxyz=1,3
          do ixyz=1,3
            pt(ixyz,jxyz)=pt(ixyz,jxyz)
     &           +2d0*eki(ixyz,jxyz,i) +strs(ixyz,jxyz,i)*avol
          enddo
        enddo
      enddo

      ptnsr(1:3,1:3)= 0d0
      call mpi_allreduce(pt,ptnsr,9,mpi_double_precision,mpi_sum
     &     ,mpi_md_world,ierr)
      ptnsr(1:3,1:3)=ptnsr(1:3,1:3)/vol

c.....Acceralation of h-matrix
      ah(1:3,1:3)= 0d0
      ptnsr(1,1)= ptnsr(1,1)-ptgt
      ptnsr(2,2)= ptnsr(2,2)-ptgt
      ptnsr(3,3)= ptnsr(3,3)-ptgt
      do jxyz=1,3
        do ixyz=1,3
          do kxyz=1,3
            ah(ixyz,jxyz)=ah(ixyz,jxyz)
     &           +ptnsr(ixyz,kxyz)*sgm(kxyz,jxyz)
          enddo
        enddo
      enddo
      ah(1:3,1:3)=ah(1:3,1:3)/vm
c.....Restore ptnsr
      ptnsr(1,1)= ptnsr(1,1)+ptgt
      ptnsr(2,2)= ptnsr(2,2)+ptgt
      ptnsr(3,3)= ptnsr(3,3)+ptgt

      end subroutine force_isobaric
c=======================================================================
      subroutine write_force(ionum,cnum,natm,h,tag,aa,acon,nismax
     &     ,epot)
      implicit none
      integer,intent(in):: natm,nismax,ionum
      character(len=*),intent(in):: cnum 
      real(8),intent(in):: h(3,3),acon(nismax),aa(3,natm),tag(natm)
     &     ,epot

      integer:: i,is
      real(8):: at(3)

c.....Write out forces
      open(ionum,file='frc'//trim(cnum),status='replace')
      write(ionum,'(i10)') natm
      do i=1,natm
        is= int(tag(i))
        at(1:3)= h(1:3,1)*aa(1,i) +h(1:3,2)*aa(2,i) +h(1:3,3)*aa(3,i)
        at(1:3)= at(1:3)/acon(is)
        write(ionum,'(3f20.7)') at(1:3)
      enddo
      close(ionum)

c.....Write out energy
      open(ionum+1,file='erg'//trim(cnum),status='replace')
      write(ionum+1,'(f20.7)') epot
      close(ionum+1)

      end subroutine write_force
c=======================================================================
      subroutine setv()
      use variables
      implicit none
      include 'mpif.h'
      include 'params_unit.h'
      integer:: i,l,is
      real(8):: dseed,sumvx,sumvy,sumvz,rnd1,rnd2,tmp,facv(nismax)
      real(8),parameter:: pi = 3.14159265358979d0

      facv(1:nismax)=dsqrt(2d0*tinit*fkb*ev2j /(am(1:nismax)*amu2kg))
     &     *m2ang /s2fs

c-----velocities in Maxwell-Boltzmann distribution
      dseed=12345
      do i=1,natm
        is= int(tag(i))
        do l=1,3
          call myrnd(rnd1,dseed)
          call myrnd(rnd2,dseed)
          va(l,i)=facv(is)*dsqrt(-dlog(rnd1))*dcos(2d0*pi*rnd2)
        enddo
      enddo

      do i=1,natm
        va(1,i)= va(1,i) /h(1,1,0) *dt
        va(2,i)= va(2,i) /h(2,2,0) *dt
        va(3,i)= va(3,i) /h(3,3,0) *dt
      enddo

      end subroutine setv
c=======================================================================
      subroutine myrnd(rnd,dseed)
      real*8 rnd,dseed
      real*8 d2p31m,d2p31
      save d2p31m,d2p31
      data d2p31m/2147483647d0/
      data d2p31 /2147483648d0/
      
      dseed=dmod(16807d0*dseed,d2p31m)
      rnd=dseed/d2p31
      return
      end subroutine myrnd
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
      subroutine get_force(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &     ,tcom,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
     &     ,mpi_md_world,myid_md,epi,epot,nismax,acon,avol
     &     ,cforce)
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
      implicit none
      integer,intent(in):: namax,natm,nnmax,nismax
      integer,intent(in):: nb,nbmax,lsb(0:nbmax,6),lsrc(6),myparity(3)
     &   ,nn(6),mpi_md_world,myid_md
      integer,intent(in):: lspr(0:nnmax,namax)
      real(8),intent(in):: ra(3,namax),h(3,3,0:1),hi(3,3),sv(3,6)
     &     ,acon(nismax),rc,tag(namax)
      real(8),intent(inout):: tcom,avol
      real(8),intent(out):: aa(3,namax),epi(namax),epot,strs(3,3,namax)
      character(len=*),intent(in):: cforce
      
      integer:: ierr
      
      select case (trim(cforce))
        case ('LJ_Ar')
          call force_LJ_Ar(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &         ,tcom,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
     &         ,mpi_md_world,myid_md,epi,epot,nismax,acon,avol)
        case ('Ito3_WHe')
          call force_Ito3_WHe(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &         ,tcom,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
     &         ,mpi_md_world,myid_md,epi,epot,nismax,acon,avol)
        case ('RK_WHe')
          call force_RK_WHe(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &         ,tcom,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
     &         ,mpi_md_world,myid_md,epi,epot,nismax,acon,avol)
        case ('RK_FeH')
          call force_RK_FeH(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &         ,tcom,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
     &         ,mpi_md_world,myid_md,epi,epot,nismax,acon,avol)
        case ('Ramas_FeH')
          call force_Ramas_FeH(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &         ,tcom,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
     &         ,mpi_md_world,myid_md,epi,epot,nismax,acon,avol)
        case ('Ackland_Fe')
          call force_Ackland_Fe(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &         ,tcom,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
     &         ,mpi_md_world,myid_md,epi,epot,nismax,acon,avol)
        case ('SW_Si')
          call force_SW_Si(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &         ,tcom,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
     &         ,mpi_md_world,myid_md,epi,epot,nismax,acon,avol)
        case ('EDIP_Si')
          call force_EDIP_Si(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &         ,tcom,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
     &         ,mpi_md_world,myid_md,epi,epot,nismax,acon,avol)
        case ('Brenner')
          call force_Brenner(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &         ,tcom,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
     &         ,mpi_md_world,myid_md,epi,epot,nismax,acon,avol)
        case ('Brenner_vdW')
          call force_Brenner_vdW(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &         ,tcom,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
     &         ,mpi_md_world,myid_md,epi,epot,nismax,acon,avol)
        case ('Lu_WHe')
          call force_Lu_Whe(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &         ,tcom,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
     &         ,mpi_md_world,myid_md,epi,epot,nismax,acon,avol)
        case ('Branicio_AlN')
          call force_Branicio_AlN(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &         ,tcom,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
     &         ,mpi_md_world,myid_md,epi,epot,nismax,acon,avol)
        case ('Mishin_Al')
          call force_Mishin_Al(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &         ,tcom,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
     &         ,mpi_md_world,myid_md,epi,epot,nismax,acon,avol)
        case ('AFS_W')
          call force_AFS_W(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &         ,tcom,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
     &         ,mpi_md_world,myid_md,epi,epot,nismax,acon,avol)
        case ('SC_Fe')
          call force_SC_Fe(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &         ,tcom,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
     &         ,mpi_md_world,myid_md,epi,epot,nismax,acon,avol)
        case ('SM_Al')
          call force_SM_Al(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &         ,tcom,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
     &         ,mpi_md_world,myid_md,epi,epot,nismax,acon,avol)
        case default
          write(6,'(a)') ' [get_force] unknown force: '
     &         //trim(cforce)//' !!'
          call mpi_finalize(ierr)
          stop
      end select

      end subroutine get_force
c=======================================================================
#ifdef __WALL__
      subroutine wall_reflect()
      use variables
      use wall
      implicit none
      include "mpif.h"
      integer:: i,is
      real(8):: zi
      
      do i=1,natm
        is= int(tag(i))
        if( is.ne.2 ) cycle
        zi= (ra(3,i)+sorg(3))
c.....if atom over wtop (top wall), go downward
        if( zi.ge.wtop ) then
          va(3,i)=-abs(va(3,i))
          va(3,i)= min(va(3,i),dwtop)
          ptop=ptop +am(is)*2d0*abs(va(3,i))*h(3,3,0)/dt
          nitop=nitop +1
        endif
c.....if atom under wbot (bottom wall), go upward
        if( zi.le.wbot ) then
          va(3,i)= abs(va(3,i))
          va(3,i)= max(va(3,i),dwbot)
          pbot=pbot +am(is)*2d0*abs(va(3,i))*h(3,3,0)/dt
          nibot=nibot +1
        endif
      enddo
      end subroutine wall_reflect
c=======================================================================
      subroutine update_wall()
c
c  Change of wall with relaxation time.
c  See RK's note on 2013.02.13
c
      use variables
      use wall
      implicit none
      include "mpif.h"
      integer:: i,itmp,is
      real(8):: zi,tmp,vtop,vbot,dvtop,dvbot
      logical,save:: l1st=.true.

      if(l1st) then
        if( wtop.lt.wbot ) stop ' [Error] wtop.lt.wbot !!!'
c.....btop,bbot: top/bottom of bulk material are fixed
        btop=0d0
        bbot=1d0
        do i=1,natm
          is= int(tag(i))
          if( is.ne.1 ) cycle
          zi= (ra(3,i)+sorg(3))
          btop= max(btop,zi)
          bbot= min(bbot,zi)
        enddo
        tmp= btop
        call mpi_allreduce(tmp,btop,1,mpi_double_precision,mpi_max
     &       ,mpi_md_world,ierr)
        tmp= bbot
        call mpi_allreduce(tmp,bbot,1,mpi_double_precision,mpi_min
     &       ,mpi_md_world,ierr)
c.....area of wall
        area_wall= h(2,2,0)*h(1,1,0)
        if( myid_md.eq.0 ) then
          write(6,'(a,es12.4)') ' btop =',btop
          write(6,'(a,es12.4)') ' bbot =',bbot
          write(6,'(a,es12.4)') ' area =',area_wall
        endif
c.....reset pressure
        ptop= 0d0
        pbot= 0d0
        nitop= 0
        nibot= 0
        l1st=.false.
c.....return at 1st call
        return
      endif

      itmp= nitop
      call mpi_reduce(itmp,nitop,1,mpi_integer,mpi_sum,0,
     &     mpi_md_world,ierr)
      itmp= nibot
      call mpi_reduce(itmp,nibot,1,mpi_integer,mpi_sum,0,
     &     mpi_md_world,ierr)
      tmp= ptop
      call mpi_reduce(tmp,ptop,1,mpi_double_precision,mpi_sum,0,
     &     mpi_md_world,ierr)
      ptop= ptop /nodes_md
      tmp= pbot
      call mpi_reduce(tmp,pbot,1,mpi_double_precision,mpi_sum,0,
     &     mpi_md_world,ierr)
      pbot= pbot /nodes_md

      if( myid_md.eq.0 ) then
        if( nitop.ne.0 ) then
          ptop= ptop /area_wall /dt_wall
c.....current volume
          vtop= abs(wtop-btop)*h(3,3,0) *area_wall
c.....change of volume
          dvtop= -vtop*(ptgt_wall -ptop)/ptop *dt_wall/trlx_wall
c.....speed of wall
          dwtop= abs(wtop-btop) *dvtop/vtop /nout_wall
cc.....change of wall position
c          wtop= btop +abs(wtop-btop)*(1d0+dvtop/vtop)
c          wtop= min(wtop,0.99d0)
        endif
        if( nibot.ne.0 ) then
          pbot= pbot /area_wall /dt_wall
          vbot= abs(wbot-bbot)*h(3,3,0) *area_wall
          dvbot= -vbot*(ptgt_wall -pbot)/pbot *dt_wall/trlx_wall
          dwbot= -abs(wbot-bbot) *dvbot/vbot /nout_wall
c          wbot= bbot -abs(wbot-bbot)*(1d0+dvbot/vbot)
c          wbot= max(wbot,0.01d0)
        endif
      endif

c      call mpi_bcast(wtop,1,mpi_double_precision,0,mpi_md_world,ierr)
c      call mpi_bcast(wbot,1,mpi_double_precision,0,mpi_md_world,ierr)
      call mpi_bcast(dwtop,1,mpi_double_precision,0,mpi_md_world,ierr)
      call mpi_bcast(dwbot,1,mpi_double_precision,0,mpi_md_world,ierr)
      
      end subroutine update_wall
#endif
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make pmd"
c     End:
