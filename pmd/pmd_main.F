      program pmd
c-----------------------------------------------------------------------
c                     Last-modified: <2016-10-30 08:26:50 Ryo KOBAYASHI>
c-----------------------------------------------------------------------
c Spatial decomposition parallel molecular dynamics program.
c Core part is separated to pmd_core.F.
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

c.....only 0th-node reads pmd0000 file
      if( myid_md.eq.0 ) then
        if( trim(ciofmt).eq.'bin' .or. trim(ciofmt).eq.'binary' ) then
          call read_pmdtot_bin(20,"pmd0000")
        else if( trim(ciofmt).eq.'ascii' ) then
          call read_pmdtot_ascii(20,"pmd0000")
        endif
      endif

c.....call pmd_core to perfom MD
      call pmd_core()

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

      end program pmd
c=======================================================================
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
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make pmd2"
c     End:
