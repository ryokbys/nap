      program pmd
c-----------------------------------------------------------------------
c                     Last-modified: <2019-05-08 15:08:11 Ryo KOBAYASHI>
c-----------------------------------------------------------------------
c Spatial decomposition parallel molecular dynamics program.
c Core part is separated to pmd_core.F.
c-----------------------------------------------------------------------
c INPUT FILES:
c ------------
c   in.pmd:     Main input file
c   pmdini:     Cell info and atom coordinations
c
c OUTPUT FILES:
c -------------
c   out.erg:    Total, kinetic, and potential energy
c   pmd_####:   Cell info and atom coordinations of a certain steps
c               in the MD run if required by "flag_out_pmd" in in.pmd.
c   out.stress: Stress component normal to z-upper surface of nanorod,
c               and z-strain of the nanorod.
c-----------------------------------------------------------------------
      use pmdio
      use pmdmpi
      use version
      use force
      use Coulomb, only: cterms
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
      integer,external:: itotOf

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
        write(6,'(a)') '================================='
     &       //'======================================='
        write(6,'(a)') ' PMD: A Parallel Molecular Dynamics program '
        write(6,*) ''
        call write_version()
        call write_authors()
        write(6,'(a)') '================================='
     &       //'======================================='
        write(6,*) ''
        call time_stamp(' Job started')
        write(6,*) ''
        write(6,'(a,i0)') ' Number of processes in MPI = ',nodes_md
        call read_input(10,trim(cinpmd))
        call check_cmin(cmin,ifdmp)
        call write_initial_setting()
c        call write_inpmd(10,trim(cinpmd))
        if( trim(ctctl).eq.'ttm' ) then
c.....Set x-boundary free if TTM
          if( boundary(1:1).ne.'f' ) then
            print *,''
            print *,'Since the two-temperature model (TTM) MD...'
            print *,'  - Free boundary condition is set'
     &           //' for x direction.'
          endif
          boundary(1:1) = 'f'
c.....Damping off if TTM
          print *,'  - Set damping off.'
          ifdmp = 0
c.....Set nrmtrans to 0
          print *,'  - Set removal of translation off.'
          nrmtrans = 0
        endif
      endif

      call bcast_params()

c.....Check ensemble
      if( myid_md.eq.0 ) call check_ensemble()

c.....Set Coulomb flag here
      if( use_force('screened_Coulomb') ) then
        ifcoulomb = 1
        cterms = 'screened'
      else if( use_force('Ewald') ) then
        ifcoulomb = 2
        cterms = 'full'
      else if( use_force('Ewald_long') ) then
        ifcoulomb = 3
        cterms = 'long'
      endif

c.....only 0th-node reads pmdini file
      if( myid_md.eq.0 ) then
c        call system('cp pmd0000 pmd0000.orig')
        if( trim(ciofmt).eq.'bin' .or. trim(ciofmt).eq.'binary' ) then
          write(6,*) 'Read pmdini in binary mode.'
          call read_pmdtot_bin(20,trim(cpmdini))
        else if( trim(ciofmt).eq.'ascii' ) then
          write(6,*) 'Read pmdini in ascii mode.'
          call read_pmdtot_ascii(20,trim(cpmdini))
        else
          write(6,*) 'Error: io_format must be either ascii, '
     &         //'bin or binary.'
          stop
        endif

        allocate(chgtot(ntot0),chitot(ntot0))
        chitot(1:ntot0) = 0d0
        call set_atomic_charges(ntot0,chgtot,tagtot,nspmax
     &       ,chgfix,schg,myid_md,iprint)

c.....Make ntot and ?tot() not null in nodes myid_md != 0
      else
        ntot0 = 1
        allocate(tagtot(ntot0),rtot(3,ntot0),vtot(3,ntot0),epitot(ntot0)
     &       ,ekitot(3,3,ntot0),stot(3,3,ntot0),atot(3,ntot0)
     &       ,chgtot(ntot0),chitot(ntot0))
        chitot(1:ntot0) = 0d0
      endif

c.....Add PKA velocity to some atom
      if( pka_energy .gt. 0d0 ) then
        call add_pka_velocity(myid_md)
      endif

c.....call pmd_core to perfom MD
      call pmd_core(hunit,h,ntot0,tagtot,rtot,vtot,atot,stot
     &     ,ekitot,epitot,chgtot,chitot,nstp,nerg,npmd
     &     ,myid_md,mpi_md_world,nodes_md,nx,ny,nz
     &     ,nspmax,am,dt,vardt_len,ciofmt,ifpmd,rc,rbuf,rc1nn,ifdmp,dmp
     &     ,minstp,tinit,tfin,ctctl,ttgt,trlx,ltdst,ntdst,nrmtrans,cpctl
     &     ,stgt,ptgt,pini,pfin,srlx,stbeta,strfin,lstrs0,lcellfix,fmv
     &     ,stnsr,epot,ekin,n_conv,ifcoulomb,czload_type,zskin_width
     &     ,zshear_angle,eps_conv,ifsort,iprint,nstp_done,lvc,boundary
     &     ,lmetaD,lconst,lrdcfrc,lcna,iter_cna)

      if( myid_md.eq.0 ) then
        if( trim(ciofmt).eq.'bin' .or. trim(ciofmt).eq.'binary' ) then
          call write_pmdtot_bin(20,cpmdfin)
        elseif( trim(ciofmt).eq.'ascii' ) then
          call write_pmdtot_ascii(20,cpmdfin)
        endif
      endif

c.....write energy, forces and stresses only for fitpot
      if( myid_md.eq.0 ) then
        call write_force(21,'.pmd',h,epot,ntot,tagtot,atot,stnsr)
        print *,''
        call time_stamp(' Job finished')
      endif

      deallocate(tagtot,rtot,vtot,epitot,ekitot,stot,atot)
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
      use force
      implicit none 
      integer:: i

      write(6,*) ''
      write(6,'(a)') '---------------------------------'
     &       //'---------------------------------------'
      write(6,'(a)') '           Initial setting   '
      write(6,'(a)') '---------------------------------'
     &       //'---------------------------------------'
      write(6,'(2x,a,5x,i0)')   'num_nodes_x',nx
      write(6,'(2x,a,5x,i0)')   'num_nodes_y',ny
      write(6,'(2x,a,5x,i0)')   'num_nodes_z',nz
      write(6,'(2x,a)') ''
      write(6,'(2x,a,5x,a)') 'io_format',ciofmt
      write(6,'(2x,a,5x,i0)') 'print_level',iprint
      write(6,'(2x,a)') ''
      write(6,'(2x,a,5x,f0.3)') 'time_interval',dt
      write(6,'(2x,a,5x,i0)')   'num_iteration',nstp
      write(6,'(2x,a,5x,i0)')   'num_out_energy',nerg
      write(6,'(2x,a)') ''
      write(6,'(2x,a,5x,i0)')   'flag_out_pmd',ifpmd
      write(6,'(2x,a,5x,i0)')   'num_out_pmd',npmd
      write(6,'(2x,a)') ''
      write(6,'(2x,a,10(2x,a))') 'force_type   ',
     &     (trim(force_list(i)),i=1,num_forces)
      if( lvc ) then
        write(6,'(2x,a,5x,l2)') 'charge_optimize',lvc
      endif
      write(6,'(2x,a,5x,f6.3)') 'cutoff_radius',rc
      write(6,'(2x,a,5x,f6.3)') 'cutoff_buffer',rbuf
      write(6,'(2x,a)') ''
      if( cmin.ne.'' ) then
        write(6,'(2x,a,5x,a)') 'minimization',trim(cmin)
      endif
      write(6,'(2x,a,5x,i0)') 'flag_damping',ifdmp
      write(6,'(2x,a,5x,f0.4)') 'damping_coeff',dmp
      write(6,'(2x,a,5x,f0.4)') 'converge_eps',eps_conv
      write(6,'(2x,a,5x,i0)')   'converge_num',n_conv
      write(6,'(2x,a,5x,i0)')   'min_iteration',minstp
      write(6,'(2x,a)') ''
c.....temperature control
      write(6,'(2x,a,5x,f0.4)') 'initial_temperature',tinit
      write(6,'(2x,a,5x,a)') 'temperature_control',trim(ctctl)
      if( trim(ctctl).eq.'Berendsen' .or.
     &     trim(ctctl).eq.'Langevin' )  then
        if( tfin.ge.0d0 ) then
          write(6,'(2x,a,5x,f0.4)') 'final_temperature',tfin
        else
          do i=1,9
            write(6,'(2x,a,i3,f8.1)') 'temperature_target',i,ttgt(i)
          enddo
        endif
        write(6,'(2x,a,5x,f0.1)') 'temperature_relax_time',trlx
        if( ltdst ) then
c.....temperature distribution
          write(6,'(2x,a)') ''
          write(6,'(2x,a,5x,l2)') 'flag_temp_dist',ltdst
          write(6,'(2x,a,5x,i0)') 'num_temp_dist',ntdst
        endif
      endif
      write(6,'(2x,a,5x,i0)') 'remove_translation',nrmtrans
      write(6,'(2x,a)') ''
c.....pressure control
      write(6,'(2x,a,5x,a)') 'stress_control',trim(cpctl)
      if( trim(cpctl).eq.'Berendsen' .or.
     &     trim(cpctl).eq.'vc-Berendsen' ) then
        write(6,'(2x,a)') 'stress_target'
        write(6,'(5x,3es11.3)') stgt(1,1:3)
        write(6,'(5x,3es11.3)') stgt(2,1:3)
        write(6,'(5x,3es11.3)') stgt(3,1:3)
        write(6,'(2x,a)') ''
      else if( trim(cpctl).eq.'vv-Berendsen' ) then
        write(6,'(2x,a,5x,f0.3)') 'pressure_target',ptgt
      endif
      write(6,*) ''
c.....strain control
      write(6,'(2x,a,5x,a)') 'zload_type',trim(czload_type)
      if( trim(czload_type).ne.'none' ) then
        write(6,'(2x,a,5x,f0.3)') 'zload_skin_width',zskin_width
        write(6,'(2x,a,5x,f0.3)') 'zload_shear_angle',zshear_angle
        write(6,'(2x,a,5x,f0.3)') 'final_strain',strfin
      endif
      write(6,'(2x,a)') ''
c.....velocity multiplying factor
      write(6,'(2x,a)') 'factor_direction'
      do i=0,9
        write(6,'(5x,i3,3f8.3)') i,fmv(1:3,i)
      enddo
      write(6,'(2x,a)') ''
c.....Mass
      write(6,'(2x,a)') 'mass'
      do i=1,nspmax
        write(6,'(5x,i3,f10.4)') i,am(i)
      enddo
      write(6,'(2x,a)') ''
c.....Boundary condition
      write(6,'(2x,a,5x,a)') 'boundary',trim(boundary)
      write(6,'(2x,a)') ''
c.....Charge
      write(6,'(2x,a)') 'charge'
      do i=1,nspmax
        write(6,'(5x,i3,f10.5)') i,schg(i)
      enddo
      write(6,'(2x,a,5x,a)') 'fix_charge',trim(chgfix)
      write(6,'(2x,a,5x,l2)') 'variable_charge',lvc
      write(6,'(2x,a)') ''

      if( lmetaD ) write(6,'(2x,a,5x,l2)') 'metadynamics',lmetaD
      if( lconst ) write(6,'(2x,a,5x,l2)') 'constraints',lconst

      write(6,'(a)') '---------------------------------'
     &       //'---------------------------------------'

      end subroutine write_initial_setting
c=======================================================================
      subroutine write_inpmd(ionum,cfname)
      integer,intent(in):: ionum
      character(len=*),intent(in):: cfname

      character(len=128):: cline
      
      write(6,*) ''
      write(6,'(a)') '---------------------------------'
     &       //'---------------------------------------'
      write(6,'(a)') '                           in.pmd'
      write(6,'(a)') '---------------------------------'
     &       //'---------------------------------------'

      open(ionum,file=trim(cfname),status='old')
      do while(.true.)
        read(ionum,'(a)',end=10) cline
        write(6,'(a)') trim(cline)
      enddo
 10   close(ionum)
      write(6,'(a)') '---------------------------------'
     &       //'---------------------------------------'
      
      end subroutine write_inpmd
c=======================================================================
      subroutine bcast_params()
      use pmdio
      use pmdmpi
      use force
      implicit none
      include 'mpif.h'

      integer:: ierr

      if( myid_md.eq.0 ) write(6,'(/,a)') ' Broadcast data to be shared'
     &     //' with all the nodes.'
c-----Broadcast input parameters to all nodes
      call mpi_bcast(nx,1,MPI_INTEGER,0,mpi_md_world,ierr)
      call mpi_bcast(ny,1,MPI_INTEGER,0,mpi_md_world,ierr)
      call mpi_bcast(nz,1,MPI_INTEGER,0,mpi_md_world,ierr)
      call mpi_bcast(nstp,1,MPI_INTEGER,0,mpi_md_world,ierr)
      call mpi_bcast(minstp,1,MPI_INTEGER,0,mpi_md_world,ierr)
      call mpi_bcast(dt,1,MPI_REAL8,0,mpi_md_world,ierr)
      call mpi_bcast(vardt_len,1,MPI_REAL8,0,mpi_md_world,ierr)
      call mpi_bcast(rc,1,MPI_REAL8,0,mpi_md_world,ierr)
      call mpi_bcast(rc1nn,1,MPI_REAL8,0,mpi_md_world,ierr)
      call mpi_bcast(rbuf,1,MPI_REAL8,0,mpi_md_world,ierr)
      call mpi_bcast(cmin,20,mpi_character,0,mpi_md_world,ierr)
      call mpi_bcast(dmp,1,MPI_REAL8,0,mpi_md_world,ierr)
      call mpi_bcast(eps_conv,1,MPI_REAL8,0,
     &     mpi_md_world,ierr)
      call mpi_bcast(n_conv,1,mpi_integer,0,mpi_md_world,ierr)
      call mpi_bcast(tinit,1,mpi_real8,0,mpi_md_world,ierr)
      call mpi_bcast(tfin,1,mpi_real8,0,mpi_md_world,ierr)
      call mpi_bcast(ctctl,20,mpi_character,0,mpi_md_world,ierr)
      call mpi_bcast(ttgt,9,MPI_REAL8,0,mpi_md_world,ierr)
      call mpi_bcast(trlx,1,MPI_REAL8,0,mpi_md_world,ierr)
      call mpi_bcast(nerg,1,MPI_INTEGER,0,mpi_md_world,ierr)
      call mpi_bcast(ifpmd,1,MPI_INTEGER,0,mpi_md_world,ierr)
      call mpi_bcast(npmd,1,MPI_INTEGER,0,mpi_md_world,ierr)
      call mpi_bcast(ifsort,1,mpi_integer,0,mpi_md_world,ierr)
      call mpi_bcast(ifdmp,1,MPI_INTEGER,0,mpi_md_world,ierr)
      call mpi_bcast(iprint,1,mpi_integer,0,mpi_md_world,ierr)
      call mpi_bcast(fmv,30,mpi_real8,0,mpi_md_world,ierr)
      call mpi_bcast(shrst,1,mpi_real8,0,mpi_md_world,ierr)
      call mpi_bcast(cpctl,20,mpi_character,0,mpi_md_world,ierr)
      call mpi_bcast(ptgt,1,mpi_real8,0,mpi_md_world,ierr)
      call mpi_bcast(pini,1,mpi_real8,0,mpi_md_world,ierr)
      call mpi_bcast(pfin,1,mpi_real8,0,mpi_md_world,ierr)
      call mpi_bcast(srlx,1,mpi_real8,0,mpi_md_world,ierr)
      call mpi_bcast(stgt,9,mpi_real8,0,mpi_md_world,ierr)
      call mpi_bcast(lcellfix,9,mpi_logical,0,mpi_md_world,ierr)
      call mpi_bcast(czload_type,128,mpi_character,0,mpi_md_world,ierr)
      call mpi_bcast(zskin_width,1,mpi_real8,0,mpi_md_world,ierr)
      call mpi_bcast(zshear_angle,1,mpi_real8,0,mpi_md_world,ierr)
      call mpi_bcast(strfin,1,mpi_real8,0,mpi_md_world,ierr)
      call mpi_bcast(am,nspmax,mpi_real8,0,mpi_md_world,ierr)
      call mpi_bcast(ciofmt,6,mpi_character,0,mpi_md_world,ierr)
      call mpi_bcast(nrmtrans,6,mpi_integer,0,mpi_md_world,ierr)
      call mpi_bcast(lstrs0,1,mpi_logical,0,mpi_md_world,ierr)
      call mpi_bcast(boundary,3,mpi_character,0,mpi_md_world,ierr)
      call mpi_bcast(pka_energy,1,mpi_real8,0,mpi_md_world,ierr)
c.....Charge related
      call mpi_bcast(lvc,1,mpi_logical,0,mpi_md_world,ierr)
      call mpi_bcast(chgfix,20,mpi_character,0,mpi_md_world,ierr)
c.....Force-fields
      call mpi_bcast(cforce,20,mpi_character,0,mpi_md_world,ierr)
      call mpi_bcast(num_forces,1,mpi_integer,0,mpi_md_world,ierr)
      if( num_forces.eq.0 ) then
        if( myid_md.eq.0 ) write(6,'(a)')
     &       ' Error: no force-field specified'
        call mpi_finalize(ierr)
        stop
      endif
      if( myid_md.ne.0 ) then
        allocate(force_list(num_forces))
      endif
      call mpi_bcast(force_list,128*num_forces,mpi_character
     &     ,0,mpi_md_world,ierr)
      call mpi_bcast(ifcoulomb,1,mpi_integer,0,mpi_md_world,ierr)
      call mpi_bcast(schg,nspmax,mpi_real8,0,mpi_md_world
     &     ,ierr)
c.....NEMD
      call mpi_bcast(ltdst,1,mpi_logical,0,mpi_md_world,ierr)
      if( ltdst ) then
        call mpi_bcast(ntdst,1,mpi_integer,0,mpi_md_world,ierr)
      endif
c.....Metadynamics
      call mpi_bcast(lmetaD,1,mpi_logical,0,mpi_md_world,ierr)
c.....Constratins
      call mpi_bcast(lconst,1,mpi_logical,0,mpi_md_world,ierr)
c.....Reduced force
      call mpi_bcast(lrdcfrc,1,mpi_logical,0,mpi_md_world,ierr)

c.....Structure analysis
      call mpi_bcast(lcna,1,mpi_logical,0,mpi_md_world,ierr)
      call mpi_bcast(iter_cna,1,mpi_integer,0,mpi_md_world,ierr)

      end subroutine bcast_params
c=======================================================================
      subroutine write_force(ionum,cpostfix,h,epot,ntot
     &     ,tagtot,atot,stnsr)
      implicit none
      include "./params_unit.h"
      integer,intent(in):: ionum,ntot
      character(len=*),intent(in):: cpostfix
      real(8),intent(in):: h(3,3),epot,stnsr(3,3)
      real(8),intent(inout):: tagtot(ntot),atot(3,ntot)

      integer:: i,n0,ixyz,ierr
      real(8):: at(3),ptmp(6)
      integer,parameter:: nmpi = 2
      integer,external:: itotOf


c.....Write out forces
      open(ionum,file='frc'//trim(cpostfix),status='replace')
      write(ionum,'(i10)') ntot
      do i=1,ntot
        write(ionum,'(3f12.6,i8)') atot(1:3,i),itotOf(tagtot(i))
      enddo
      close(ionum)

c.....Write out energy
      open(ionum+1,file='erg'//trim(cpostfix),status='replace')
      write(ionum+1,'(es23.14e3)') epot
      close(ionum+1)

c.....Write stress tensor, negative as compressive, positive as tensile
      ptmp(1) = stnsr(1,1)*up2gpa*(-1d0)
      ptmp(2) = stnsr(2,2)*up2gpa*(-1d0)
      ptmp(3) = stnsr(3,3)*up2gpa*(-1d0)
      ptmp(4) = stnsr(3,2)*up2gpa*(-1d0)
      ptmp(5) = stnsr(1,3)*up2gpa*(-1d0)
      ptmp(6) = stnsr(1,2)*up2gpa*(-1d0)
      open(ionum+2,file='strs'//trim(cpostfix),status='replace')
      write(ionum+2,'(6f12.4)') ptmp(1:6)
      close(ionum+2)

      end subroutine write_force
c=======================================================================
      subroutine set_atomic_charges(ntot,chg,tag,nspmax
     &     ,chgfix,schg,myid,iprint)
      implicit none 
      integer,intent(in):: ntot,nspmax,myid,iprint
      real(8),intent(in):: tag(ntot),schg(nspmax)
      character(len=*),intent(in):: chgfix
      real(8),intent(out):: chg(ntot)

      integer:: i,is

      if( trim(chgfix).eq.'input' ) then
        if( myid.eq.0 .and. iprint.ne.0 ) then
          print *,'Charges are set from input.'
        endif
        do i=1,ntot
          is = int(tag(i))
          chg(i)= schg(is)
        enddo
      endif
      
      end subroutine set_atomic_charges
c=======================================================================
      subroutine check_ensemble()
      use pmdio
      implicit none

      logical:: l_temp
      character:: c_2*1,c_3*1,c_ensemble*3

      c_2 = 'V'
      c_3 = 'E'

      l_temp = .false.
      
      if(  trim(ctctl).eq.'Langeven' .or.
     &     trim(ctctl).eq.'Berendsen' .or.
     &     trim(ctctl).eq.'ttm' ) then
        l_temp = .true.
        c_3 = 'T'
      endif

      if(  trim(cpctl).eq.'Berendsen' .or.
     &     trim(cpctl).eq.'vc-Berendsen' ) then
        c_2 = 'p'
        if( .not. l_temp ) c_3 = 'H'
      else if( trim(cpctl).eq.'vv-Berendsen' ) then
        c_2 = 'P'
        if( .not. l_temp ) c_3 = 'H'
      endif

      c_ensemble = 'N'//c_2//c_3
      print *,''
      print *,'Ensemble = ',trim(c_ensemble)
      
      end subroutine check_ensemble
c=======================================================================
      subroutine add_pka_velocity(myid_md)
c.....
c.....Add PKA velocity to some atom from a given PKA energy
c.....
      use pmdio
      implicit none
      include './params_unit.h'
      integer,intent(in):: myid_md
      real(8),parameter:: pi= 3.14159265358979d0

      integer:: i,icntr,is
      real(8):: dmin,d,theta,phi,rx,ry,rz,vel,vx,vy,vz
      interface
        function urnd(dseed0)
        real(8),intent(in),optional:: dseed0
        real(8):: urnd
        end function urnd
      end interface
c      real(8),external:: urnd

      if( myid_md.eq.0 ) then
        if( iatom_pka.le.0 ) then
c.....Find a center atom ICNTR
          icntr = 0
          dmin = 1d+10
          do i=1,ntot0
            d = (rtot(1,i)-0.5d0)**2
     &           +(rtot(2,i)-0.5d0)**2
     &           +(rtot(3,i)-0.5d0)**2
            if( d .lt. dmin ) then
              icntr = i
              dmin = d
            endif
          enddo
        else
          if( iatom_pka.gt.ntot0 ) then
            print *,'ERROR: iatom_pka.gt.ntot0 !'
            stop
          endif
          icntr = iatom_pka
        endif

c.....Add PKA velocity at ICNTR atom
c.....Get random theta and phi
        theta = 90d0 *urnd() /180d0 *pi
        phi = 90d0 *urnd() /180d0 *pi
        rx = sin(theta)*cos(phi)
        ry = sin(theta)*sin(phi)
        rz = cos(theta)
        is = int(tagtot(icntr))
c.....[eV] to [Ang/fs]
        vel = sqrt(pka_energy*ev2j *2d0 /(am(is)*ump2kg)) *m2ang /s2fs
        vx = rx*vel
        vy = ry*vel
        vz = rz*vel
c.....Assume an orthogonal simulation box
        vx = vx /h(1,1,0)
        vy = vy /h(2,2,0)
        vz = vz /h(3,3,0)
        vtot(1,icntr) = vtot(1,icntr) +vx
        vtot(2,icntr) = vtot(2,icntr) +vy
        vtot(3,icntr) = vtot(3,icntr) +vz

      endif
      end subroutine add_pka_velocity
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make pmd"
c     End:
