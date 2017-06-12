      program pmd
c-----------------------------------------------------------------------
c                     Last-modified: <2017-06-13 07:32:10 Ryo KOBAYASHI>
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
      use pmdio
      use pmdmpi
      implicit none
      include "mpif.h"
      include "./params_unit.h"

#ifdef __DISL__
c.....Epot threshold for disl core extraction [Hartree]
      real(8),parameter:: epith = -0.1410d0
#endif

      integer:: i,j,k,l,m,n,ia,ib,is,ifmv,nave,nspl,i_conv,ierr
     &     ,nstp_done
      real(8):: tmp,hscl(3),aai(3),ami,dt2,tave,vi(3),vl(3),epot,ekin
      logical,external:: force_on

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
        call check_cmin(cmin,ifdmp)
        call write_initial_setting()
      endif

      print *,' numff = ',numff
      print *,' cffs  = ',(cffs(i),i=1,numff)

      call bcast_params()

c.....Set Coulomb flag here
      if( force_on('screened_Coulomb') ) then
        ifcoulomb = 1
      else if( force_on('Ewald_Coulomb') ) then
        ifcoulomb = 2
      else if( force_on('vc_gaussian') ) then
        ifcoulomb = 3
      endif

c.....only 0th-node reads pmdini file
      if( myid_md.eq.0 ) then
c        call system('cp pmd0000 pmd0000.orig')
        if( trim(ciofmt).eq.'bin' .or. trim(ciofmt).eq.'binary' ) then
          call read_pmdtot_bin(20,trim(cpmdini))
        else if( trim(ciofmt).eq.'ascii' ) then
          call read_pmdtot_ascii(20,trim(cpmdini))
        else
          write(6,*) 'Error: io_format must be either ascii, '
     &         //'bin or binary.'
          stop
        endif

c.....if needed, read atomic charges from chg0000 file
c.....but even if not needed, allocate ntot long array
c.....for the simplicty of code
        allocate(chgtot(ntot),chitot(ntot))
        chitot(1:ntot) = 0d0
        call set_atomic_charges(ntot,chgtot,tagtot,nspmax,schg)

c.....make ntot and ?tot() not null in nodes myid_md != 0
      else
        ntot = 1
        allocate(tagtot(ntot),rtot(3,ntot),vtot(3,ntot),epitot(ntot)
     &       ,ekitot(3,3,ntot),stot(3,3,ntot),atot(3,ntot)
     &       ,chgtot(ntot),chitot(ntot))
        chitot(1:ntot) = 0d0
      endif

c.....call pmd_core to perfom MD
      call pmd_core(hunit,h,ntot,tagtot,rtot,vtot,atot,stot
     &     ,ekitot,epitot,chgtot,chitot,nstp,nerg,npmd
     &     ,myid_md,mpi_md_world,nodes_md,nx,ny,nz
     &     ,nspmax,am,dt,ciofmt,ifpmd,numff,cffs,rc,rbuf,ifdmp,dmp
     &     ,minstp,tinit,tfin,ctctl,ttgt,trlx,ltdst,ntdst,cpctl,stgt
     &     ,ptgt,srlx,stbeta,strfin,lstrs0,lcellfix,fmv,ptnsr,epot,ekin
     &     ,n_conv,ifcoulomb,czload_type,eps_conv,ifsort,iprint
     &     ,nstp_done)

      if( myid_md.eq.0 ) then
        if( trim(ciofmt).eq.'bin' .or. trim(ciofmt).eq.'binary' ) then
          call write_pmdtot_bin(20,cpmdfin)
        elseif( trim(ciofmt).eq.'ascii' ) then
          call write_pmdtot_ascii(20,cpmdfin)
        endif
      endif

c.....write energy, forces and stresses only for fitpot
      if( myid_md.eq.0 ) then
        call write_force(21,'.pmd',h,epot,ntot,tagtot,atot,ptnsr)
      endif

      deallocate(tagtot,rtot,vtot,epitot,ekitot,stot,atot,cffs)
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
      subroutine check_cmin(cmin,ifdmp)
      implicit none
      character,intent(in):: cmin 
      integer,intent(inout):: ifdmp

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
      use pmdio
      use pmdmpi
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
      write(6,'(2x,a25,10(2x,a))') 'force_type',(cffs(i),i=1,numff)
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
      do i=1,nspmax
        write(6,'(4x,i2,f10.3)') i,am(i)
      enddo
      write(6,'(a)') '======================================'

      end subroutine write_initial_setting
c=======================================================================
      subroutine bcast_params()
      use pmdio
      use pmdmpi
      implicit none
      include 'mpif.h'

      integer:: ierr

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
      call mpi_bcast(lcellfix,9,mpi_logical,0,mpi_md_world,ierr)
      call mpi_bcast(czload_type,5,mpi_character,0,mpi_md_world,ierr)
      call mpi_bcast(strfin,1,mpi_double_precision,0,mpi_md_world,ierr)
      call mpi_bcast(am,nspmax,mpi_double_precision,0,mpi_md_world,ierr)
      call mpi_bcast(ciofmt,6,mpi_character,0,mpi_md_world,ierr)
c.....Force-fields
      call mpi_bcast(cforce,20,mpi_character,0,mpi_md_world,ierr)
      call mpi_bcast(numff,1,mpi_integer,0,mpi_md_world,ierr)
      if( myid_md.ne.0 ) then
        allocate(cffs(numff))
      endif
      call mpi_bcast(cffs,20*numff,mpi_character,0,mpi_md_world,ierr)
      call mpi_bcast(ifcoulomb,1,mpi_integer,0,mpi_md_world,ierr)
      call mpi_bcast(schg,nspmax,mpi_double_precision,0,mpi_md_world
     &     ,ierr)
c.....NEMD
      call mpi_bcast(ltdst,1,mpi_logical,0,mpi_md_world,ierr)
      if( ltdst ) then
        call mpi_bcast(ntdst,1,mpi_integer,0,mpi_md_world,ierr)
      endif

      end subroutine bcast_params
c=======================================================================
      subroutine write_force(ionum,cpostfix,h,epot,ntot
     &     ,tagtot,atot,ptnsr)
      implicit none
      include "./params_unit.h"
      integer,intent(in):: ionum,ntot
      character(len=*),intent(in):: cpostfix
      real(8),intent(in):: h(3,3),epot,ptnsr(3,3)
      real(8),intent(inout):: tagtot(ntot),atot(3,ntot)

      integer:: i,n0,ixyz,ierr
      real(8):: at(3),ptmp(6)
      integer,parameter:: nmpi = 2


c.....Write out forces
      open(ionum,file='frc'//trim(cpostfix),status='replace')
      write(ionum,'(i10)') ntot
      do i=1,ntot
        write(ionum,'(3f12.6)') atot(1:3,i)
      enddo
      close(ionum)

c.....Write out energy
      open(ionum+1,file='erg'//trim(cpostfix),status='replace')
      write(ionum+1,'(es23.14e3)') epot
      close(ionum+1)

c.....Write stress tensor, negative as compressive, positive as tensile
      ptmp(1) = ptnsr(1,1)*up2gpa*(-1d0)
      ptmp(2) = ptnsr(2,2)*up2gpa*(-1d0)
      ptmp(3) = ptnsr(3,3)*up2gpa*(-1d0)
      ptmp(4) = ptnsr(3,2)*up2gpa*(-1d0)
      ptmp(5) = ptnsr(1,3)*up2gpa*(-1d0)
      ptmp(6) = ptnsr(1,2)*up2gpa*(-1d0)
      open(ionum+2,file='strs'//trim(cpostfix),status='replace')
      write(ionum+2,'(6f12.4)') ptmp(1:6)
      close(ionum+2)

      end subroutine write_force
c=======================================================================
      subroutine set_atomic_charges(ntot,chg,tag,nspmax,schg)
      implicit none 
      integer,intent(in):: ntot,nspmax
      real(8),intent(in):: tag(ntot),schg(nspmax)
      real(8),intent(out):: chg(ntot)

      integer:: i,is

      do i=1,ntot
        is = int(tag(i))
        chg(i)= schg(is)
      enddo
      
      end subroutine set_atomic_charges
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make pmd"
c     End:
