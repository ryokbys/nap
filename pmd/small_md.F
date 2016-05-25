      program small_md_main
c-----------------------------------------------------------------------
c                        Time-stamp: <2016-05-25 14:50:58 Ryo KOBAYASHI>
c-----------------------------------------------------------------------
c  Molecular dynamics program for systems small enough not to use 
c  cell-list method and parallelization.
c-----------------------------------------------------------------------
      use variables
      implicit none
      include "./params_unit.h"

      integer:: i,ia
      real(8):: tmp
      character:: cnum*4 
      integer,external:: itotOf

      call read_input(10,'in.smd')
      call check_cmin()
      call write_initial_setting()
c.....Multiply damping factor to avoid multiplication of every step
      if( ifdmp.eq.1 ) then
        fmv(1:3,0:9)= fmv(1:3,0:9) *dmp
      endif

      if( trim(ciofmt).ne.'ascii' ) then
        write(6,*) ' [Error] io_format must be ascii !!!'
        stop
      endif
      namax = 26000
      nbmax = namax
      allocate(tag(namax),ra(3,namax),va(3,namax),epi(namax)
     &     ,eki(3,3,namax),strs(3,3,namax),aa(3,namax))
      call read_pmd_ascii(20,"smd0000",namax,natm,h,hunit,tag,ra,va
     &     ,eki,epi,strs)
      ntot= natm
      print *, ' ntot =',ntot
      write(6,'(a)') " H-matrix:"
      write(6,'(3es15.7)') h(1:3,1,0)
      write(6,'(3es15.7)') h(1:3,2,0)
      write(6,'(3es15.7)') h(1:3,3,0)
      nis= 0
      do i=1,natm
c        ra(1:3,i)= ra(1:3,i) -sorg(1:3)
c-------scaling
        va(1:3,i)= va(1:3,i) *dt
c-------species of atom-i
        nis= max(int(tag(i)),nis)
        if(nis.gt.nismax) then
          write(6,'(a)') " nis.gt.nismax!!"
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

      call smd_core(namax,natm,h,hunit,tag,ra,va,eki,epi
     &     ,nis,nismax,am,dt,vol,tinit,tfin,ttgt,tfac,trlx
     &     ,ctctl,ndof,nbmax,nb,nnmax,aa,strs,nstp,minstp,nerg,ifpmd
     &     ,npmd,cforce,rc,rbuf,ifdmp,dmp,eps_conv,n_conv
     &     ,fmv,ekin,epot,temp,1
     &     ,nmin_fire,finc_fire,fdec_fire
     &     ,alp0_fire,falp_fire,dtmax_fire)

      call write_force(20,natm,aa)
      call write_energy(21,epot)

      write(6,'(a)') " program smd finished correctly."

      end program small_md_main
c=======================================================================
      subroutine write_initial_setting()
      use variables
      implicit none 
      integer:: i

      write(6,'(a)') '========== INITIAL SETTING ==========='
      write(6,'(2x,a25,2x,a)') 'io_format',ciofmt
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
      write(6,'(2x,a)') ''
      write(6,'(2x,a25,f12.4)') 'initial_temperature',tinit
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
      subroutine write_force(ionum,natm,aa)
      implicit none
      include "./params_unit.h"
      integer,intent(in):: ionum,natm
      real(8),intent(in):: aa(3,natm)
      
      integer:: i

c.....Write out forces
      open(ionum,file='frc.smd',status='replace')
      write(ionum,'(i10)') natm
      do i=1,natm
        write(ionum,'(3es23.14e3)') aa(1:3,i)
      enddo
      close(ionum)

      end subroutine write_force
c=======================================================================
      subroutine write_energy(ionum,erg)
      implicit none
      include "./params_unit.h"
      integer,intent(in):: ionum
      real(8),intent(in):: erg
      
c.....Write out energy
      open(ionum,file='erg.smd',status='replace')
      write(ionum,'(es23.14e3)') erg
      close(ionum)

      end subroutine write_energy
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
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make smd"
c     End:
