      program small_md
c-----------------------------------------------------------------------
c  Molecular dynamics program for systems small enough not to use 
c  cell-list method and parallelization.
c-----------------------------------------------------------------------
      use variables
      implicit none
      include "./params_unit.h"

      real(8):: tmp
      character:: cnum*4 

      call read_input(10,'in.smd')
      call write_initial_setting()
      call convert_units()

      if( trim(ciofmt).eq.'ascii' ) then
        call read_pmd_ascii(20,"smd000",namax,natm,h,tag,ra,va
     &       ,eki,epi,strs)
      else
        call read_POSCAR(20,"POSCAR",namax,natm,h,tag,ra,va)
      endif

c-----setup
      call setup(nismax,am,acon,fack,dt)
c-----set HI and SGM
      call boxmat(h,hi,ht,g,gi,gt,vol,sgm)

c.....Set initial temperature if needed
      if( tinit.gt.1d-5 ) then
        call setv()
        call rm_trans_motion()
      elseif( abs(tinit).le.1d-5 ) then
        va(1:3,1:natm)= 0d0
      endif

      call mk_lspr_brute(namax,natm,nbmax,nb,nnmax,tag,ra,rc+rbuf,h,hi
     &     ,lspr)
      call get_force(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &     ,tcom,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
     &     ,mpi_md_world,myid_md,epi,epot,nismax,acon,avol
     &     ,cforce)

      call write_force(12)
      call write_energy(13,epot0)

c-----calc kinetic energy
      call get_ekin(namax,natm,va,tag,h,nismax,fack,ekin,eki)
c-----in case of 3D, <ekin> = 3/2 *NkT
      temp= ekin /ntot/fkb *2d0/3d0
cc-----in case of 2D, <ekin> = 2/2 *NkT
c      temp= ekin /ntot/fkb *2d0/2d0

      write(6,'(1x,a)') "initial values:"
      write(6,'(1x,a,es22.14,a)') " kinetic energy  ="
     &     ,ekin*ehrt," eV"
      write(6,'(1x,a,es22.14,a)') " potential energy="
     &     ,epot0*ehrt," eV"
      write(6,'(1x,a,es15.7,a)') " temperature     =",temp," K"

      open(11,file="out.erg",status='replace')
      write(11,'(a)') '# istp,etot-epot0,ekin,epot-epot0,temp'
     &     //'vol(A^3),pressure(GPa)'
      write(11,'(a,es15.7)') '# epot0 =',epot0
      write(11,'(i8,3es15.7,f10.2,2es15.7)') istp
     &     ,(ekin+epot0-epot0)*ehrt,ekin*ehrt,(epot0-epot0)*ehrt
     &     ,temp
     &     ,vol
     &     ,(ptnsr(1,1)+ptnsr(2,2)+ptnsr(3,3))/3*aup2gpa
      call flush(11)

c.....make MD loop starts
      do istp=1,nstp
        va(1:3,1:natm)= va(1:3,1:natm) +aa(1:3,1:natm)
        ra(1:3,1:natm)= ra(1:3,1:natm) +va(1:3,1:natm)

        call mk_lspr_brute(namax,natm,nbmax,nb,nnmax,tag,ra,rc+rbuf,h,hi
     &     ,lspr)
        call get_force(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &       ,tcom,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
     &       ,mpi_md_world,myid_md,epi,epot,nismax,acon,avol
     &       ,cforce)
        va(1:3,1:natm)= va(1:3,1:natm) +aa(1:3,1:natm)

        call get_ekin(namax,natm,va,tag,h,nismax,fack,ekin,eki)
c.....temperature control by velocity scaling
        if(iftctl.eq.1) then
c.....in case of 3D, <ekin> = 3/2 *NkT
          temp= ekin /ntot/fkb *2d0/3d0
c---------Berendsen thermostat
          if( (treq-temp)/temp.gt.100d0 ) then
            tmp= sqrt(1d0 +dt/trlx*100d0 )
          else
            tmp= sqrt(1d0 +dt/trlx*(treq-temp)/temp)
          endif
          va(1:3,1:natm)= va(1:3,1:natm) *tmp
        endif

        write(11,'(i8,3es15.7,f10.2,2es15.7)') istp
     &       ,(ekin+epot-epot0)*ehrt,ekin*ehrt,(epot-epot0)*ehrt
     &       ,temp
     &       ,vol
     &       ,(ptnsr(1,1)+ptnsr(2,2)+ptnsr(3,3))/3*aup2gpa
        call flush(11)

c-------write the particle positions
        if(ifpmd.eq.1.and.mod(istp,noutpmd).eq.0)then
c---------decide pmd-file name
          iocntpmd=iocntpmd+1
          write(cnum(1:4),'(i4.4)') iocntpmd
c---------write pmd file for Akira
          if( trim(ciofmt).eq.'ascii' ) then
            call write_pmd_ascii(20,"smd"//cnum,namax,natm,h,tag,ra,va
     &           ,eki,epi,strs)
          else
            call read_POSCAR(20,"POSCAR"//cnum,namax,natm,h,tag,ra,va)
          endif
        endif

      enddo ! end of istp=1,nstp

      write(6,'(1x,a)') "final values:"
      write(6,'(1x,a,es22.14)') " kinetic energy  ="
     &     ,ekin*ehrt," eV"
      write(6,'(1x,a,es22.14)') " potential energy="
     &     ,epot*ehrt," eV"
      write(6,'(1x,a,es15.7)') " temperature     ="
     &     ,temp," K"
      write(6,'(a)') " program finished correctly."

      end program small_md
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
      subroutine convert_units()
      use variables
      implicit none
      include "./params_unit.h"

c-----sec. --> a.u.
      dt= dt /ut
      trlx= trlx /ut
c-----meter --> bohr
      rc= rc /bohr
      rbuf= rbuf /bohr
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

      end subroutine convert_units
c=======================================================================
      subroutine setup(nismax,am,acon,fack,dt)
      implicit none
      include "params_unit.h"
      integer,intent(in):: nismax
      real(8),intent(in):: dt
      real(8),intent(out):: acon(nismax),fack(nismax),am(nismax)

      integer:: i

      am(1:nismax)= am(1:nismax) *ump/ume

c-----prefactors for normalized accelerration
      do i=1,nismax
        acon(i)= 0.5d0*dt**2/am(i)
      enddo

c-----prefactors for kinetic energy, FACK
      do i=1,nismax
        fack(i)= 0.5d0*am(i)/dt**2
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

      return
      end subroutine boxmat
c=======================================================================
      subroutine setv()
      use variables
      implicit none
      include 'params_unit.h'
      integer:: i,l,is
      real(8):: dseed,sumvx,sumvy,sumvz,rnd1,rnd2,tmp,facv(nismax)
      real(8),parameter:: pi = 3.14159265358979d0

      facv(1:nismax)=dsqrt(2d0*tinit*fkb/am(1:nismax))

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
      do i=1,natm
        va(1,i)=va(1,i)-sumpx/amtot
        va(2,i)=va(2,i)-sumpy/amtot
        va(3,i)=va(3,i)-sumpz/amtot
      enddo
      write(6,'(a,3es12.4)') ' sumpx,y,z/amtot='
     &     ,sumpx/amtot,sumpy/amtot,sumpz/amtot

      end subroutine rm_trans_motion
c=======================================================================
      subroutine get_ekin(namax,natm,va,tag,h,nismax,fack,ekin,eki)
      implicit none 
      integer,intent(in):: namax,natm,nismax
      real(8),intent(in):: va(3,namax),h(3,3),fack(nismax)
     &     ,tag(namax)
      real(8),intent(out):: ekin,eki(3,3,namax)
c.....locals
      integer:: i,ierr,is,ixyz,jxyz,imax
      real(8):: x,y,z,v(3),v2

      ekin=0d0
      eki(1:3,1:3,1:natm)= 0d0

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
        ekin=ekin +eki(1,1,i) +eki(2,2,i) +eki(3,3,i)
      enddo

      end subroutine get_ekin
c=======================================================================
      subroutine write_force(ionum)
      use variables
      implicit none
      include "./params_unit.h"
      integer,intent(in):: ionum
      
      integer:: i,is
      real(8):: at(3)

c.....Write out forces
      open(ionum,file='frc.smd',status='replace')
      write(ionum,'(i10)') natm
      do i=1,natm
        is= int(tag(i))
        at(1:3)= h(1:3,1,0)*aa(1,i)
     &       +h(1:3,2,0)*aa(2,i)
     &       +h(1:3,3,0)*aa(3,i)
        at(1:3)= at(1:3)/acon(is)
        write(ionum,'(3es22.14)') at(1:3)*ehrt/bohr2aa
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
      write(ionum,'(es22.14)') erg*ehrt
      close(ionum)

      end subroutine write_energy
c=======================================================================
      subroutine get_force(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &     ,tcom,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
     &     ,mpi_md_world,myid_md,epi,epot,nismax,acon,avol
     &     ,cforce)
      use LJ_Ar,only:force_LJ_Ar
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
        case default
          write(6,'(a)') ' [get_force] unknown force: '
     &         //trim(cforce)//' !!'
          stop
      end select

      end subroutine get_force
c=======================================================================
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make smd"
c     End:
