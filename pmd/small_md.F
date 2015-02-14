      program small_md
c-----------------------------------------------------------------------
c                        Time-stamp: <2015-02-12 10:44:04 Ryo KOBAYASHI>
c-----------------------------------------------------------------------
c  Molecular dynamics program for systems small enough not to use 
c  cell-list method and parallelization.
c  Brute-force searching is used to make a neighbor list with copying
c  the original cell to make image cells. Thus it requires 26 times
c  large number of atoms than the original cell.
c-----------------------------------------------------------------------
      use variables
      implicit none
      include "./params_unit.h"

      integer:: i,ia
      real(8):: tmp
      character:: cnum*4 
      integer,external:: itotOf

c.....set myid to indicate that it is not a parallel computation
      myid_md= -1

      call read_input(10,'in.smd')
      call write_initial_setting()
      call convert_units()

      if( trim(ciofmt).ne.'ascii' ) then
        write(6,*) ' [Error] io_format must be ascii !!!'
        stop
      endif
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

      call mk_lspr_brute(namax,natm,nbmax,nb,nnmax,tag,ra,rc+rbuf
     &     ,h,hi,sgm,lspr)
      call get_force(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &     ,tcom,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
     &     ,mpi_md_world,myid_md,epi,epot0,nismax,acon,avol
     &     ,cforce)
      epot= epot0

      call write_force(12)
      call write_energy(13,epot0)

c-----calc kinetic energy
      call get_ekin(namax,natm,va,tag,h,nismax,fack,ekin,eki)
c-----in case of 3D, <ekin> = 3/2 *NkT
      temp(1)= ekin /ntot/fkb *2d0/3d0
cc-----in case of 2D, <ekin> = 2/2 *NkT
c      temp= ekin /ntot/fkb *2d0/2d0

      write(6,'(1x,a)') "initial values:"
      write(6,'(1x,a,es22.14,a)') " kinetic energy  ="
     &     ,ekin," eV"
      write(6,'(1x,a,es22.14,a)') " potential energy="
     &     ,epot0," eV"
      write(6,'(1x,a,es15.7,a)') " temperature     =",temp(1)," K"

      open(ioerg,file="out.erg",status='replace')
      write(ioerg,'(a)') '# istp,etot-epot0,ekin,epot-epot0,temp'
     &     //'vol(A^3),pressure(GPa)'
      write(ioerg,'(a,es15.7)') '# epot0 =',epot0
      write(ioerg,'(i8,3es15.7,f10.2,2es15.7)') istp
     &     ,(ekin+epot0-epot0),ekin,(epot0-epot0)
     &     ,temp(1)
     &     ,vol
     &     ,(ptnsr(1,1)+ptnsr(2,2)+ptnsr(3,3))/3*up2gpa
      call flush(ioerg)


c.....make MD loop starts
      do istp=1,nstp
        va(1:3,1:natm)= va(1:3,1:natm) +aa(1:3,1:natm)
        ra(1:3,1:natm)= ra(1:3,1:natm) +va(1:3,1:natm)

        call mk_lspr_brute(namax,natm,nbmax,nb,nnmax,tag,ra,rc+rbuf
     &       ,h,hi,sgm,lspr)
        call get_force(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &       ,tcom,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
     &       ,mpi_md_world,myid_md,epi,epot,nismax,acon,avol
     &       ,cforce)
        va(1:3,1:natm)= va(1:3,1:natm) +aa(1:3,1:natm)

        call get_ekin(namax,natm,va,tag,h,nismax,fack,ekin,eki)
cc.....temperature control by velocity scaling
c        if(iftctl.eq.1) then
cc.....in case of 3D, <ekin> = 3/2 *NkT
c          temp= ekin /ntot/fkb *2d0/3d0
cc---------Berendsen thermostat
c          if( (ttgt()-temp)/temp.gt.100d0 ) then
c            tmp= sqrt(1d0 +dt/trlx*100d0 )
c          else
c            tmp= sqrt(1d0 +dt/trlx*(treq-temp)/temp)
c          endif
c          va(1:3,1:natm)= va(1:3,1:natm) *tmp
c        endif
        temp(1)= ekin /ntot/fkb *2d0/3d0
        write(6,'(a,i8,f10.2)') ' istp,temp = ',istp,temp(1)
        write(ioerg,'(i8,3es15.7,f10.2,2es15.7)') istp
     &       ,(ekin+epot-epot0),ekin,(epot-epot0)
     &       ,temp(1)
     &       ,vol
     &       ,(ptnsr(1,1)+ptnsr(2,2)+ptnsr(3,3))/3*up2gpa
        call flush(ioerg)

c-------write the particle positions
        if(ifpmd.eq.1.and.mod(istp,noutpmd).eq.0)then
c---------decide pmd-file name
          iocntpmd=iocntpmd+1
          write(cnum(1:4),'(i4.4)') iocntpmd
c---------write pmd file for Akira
          call write_pmd_ascii(20,"smd"//cnum,natm,h,hunit
     &         ,tag,ra,va,eki,epi,strs,sorg,dt)
        endif

      enddo ! end of istp=1,nstp

      write(6,'(1x,a)') "final values:"
      write(6,'(1x,a,es22.14,a)') " kinetic energy  ="
     &     ,ekin," eV"
      write(6,'(1x,a,es22.14,a)') " potential energy="
     &     ,epot," eV"
      write(6,'(1x,a,es15.7,a)') " temperature     ="
     &     ,temp(1)," K"
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
      write(6,'(2x,a25,2x,a)') 'temperature_control',ctctl
      do i=1,9
        write(6,'(2x,a25,i3,f12.4)') 'temperature_target',i,ttgt(i)
      enddo
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

cc-----sec. --> a.u.
c      dt= dt /ut
c      trlx= trlx /ut
cc-----meter --> bohr
c      rc= rc /bohr
c      rbuf= rbuf /bohr
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
      real(8):: tmp

      am(1:nismax)= am(1:nismax) !*ump/ume

c.....umass*Ang^2/fs^2 ---> kg*m^2/sec^2 (=Joule)
      tmp= amu2kg *ang2m**2 /fs2s**2

c-----prefactors for normalized accelerration
      do i=1,nismax
c        acon(i)= 0.5d0*dt**2/am(i)
        acon(i)= ev2j/ang2m *0.5d0 *(dt*fs2s)**2 /(am(i)*amu2kg) *m2ang
      enddo

c-----prefactors for kinetic energy, FACK
      do i=1,nismax
c        fack(i)= 0.5d0*am(i)/dt**2
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

      return
      end subroutine boxmat
c=======================================================================
      subroutine setv()
      use variables
      implicit none
      include 'params_unit.h'
      integer:: i,l,is
      real(8):: sumvx,sumvy,sumvz,rnd1,rnd2,tmp,facv(nismax)
      real(8),parameter:: pi = 3.14159265358979d0
      real(8),external:: box_muller

      facv(1:nismax)=dsqrt(tinit*fkb*ev2j/(am(1:nismax)*amu2kg))
     &     *m2ang /s2fs

c-----velocities in Maxwell-Boltzmann distribution
      do i=1,natm
        is= int(tag(i))
        do l=1,3
          va(l,i)=facv(is) *box_muller()
        enddo
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
      real(8):: x,y,z,v(3)

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
        write(ionum,'(3es23.14e3)') at(1:3)
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
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make smd"
c     End:
