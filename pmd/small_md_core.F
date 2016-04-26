      subroutine smd_core(namax,natm,h,hunit,tag,ra,va,eki,epi
     &     ,nis,nismax,am,dt,vol,tinit,tfin,ttgt,tfac,trlx
     &     ,ctctl,ndof,nbmax,nb,nnmax,aa,strs,nstp,nerg,ifpmd,npmd
     &     ,cforce,rc,rbuf,ifdmp,dmp,eps_conv,n_conv
     &     ,fmv,ekin,epot,temp,iolevel)
c-----------------------------------------------------------------------
c                        Time-stamp: <2016-04-25 14:33:07 Ryo KOBAYASHI>
c-----------------------------------------------------------------------
c  Molecular dynamics core routine for systems small enough not to use 
c  cell-list method and parallelization.
c  Brute-force searching is used to make a neighbor list with copying
c  the original cell to make image cells. Thus it requires 26 times
c  large number of atoms than the original cell.
c-----------------------------------------------------------------------
c  iolevel: [0] for no print-out to stdout
c-----------------------------------------------------------------------
      implicit none
      include "./params_unit.h"

c.....arguments
      integer,intent(in):: namax,natm,nis,nismax,nstp,nerg,ifdmp
     &     ,npmd,nbmax,nnmax,iolevel,n_conv
      integer,intent(inout):: nb,ifpmd,ndof(9)
      real(8),intent(in):: rbuf,dmp,h(3,3,0:1),hunit,dt,fmv(3,0:9)
     &     ,am(nismax),tinit,tfin,eps_conv,trlx
      real(8),intent(inout):: tag(namax),ra(3,namax)
     &     ,va(3,namax),eki(3,3,namax),epi(namax),aa(3,namax)
     &     ,strs(3,3,namax),rc,tfac(9),ttgt(9)
      real(8),intent(out):: ekin,epot,temp(9),vol
      character(len=*),intent(in):: cforce,ctctl

c.....local variables
      integer:: i,j,k,l,m,n
      integer:: ia,is,ioerg,nouterg,noutpmd,iocntpmd,istp,ifmv
      real(8):: tmp,at(3)
      character:: cnum*4
      integer,external:: itotOf
      real(8):: acon(nismax),fack(nismax),tcom,avol,epot0,epotp
     &     ,g(3,3,0:1),gi(3,3),gt(3,3,0:1),hi(3,3),ht(3,3,0:1),sgm(3,3)
     &     ,sv(3,6),sorg(3)
      integer:: nn(6),myparity(3),lsrc(6)
     &     ,myid_md,nodes_md,mpi_md_world,myx,myy,myz,ierr
      integer:: i_conv
      integer,allocatable:: lspr(:,:),lsb(:,:)
      logical:: lconverged = .false.

      allocate(lspr(0:nnmax,natm),lsb(0:nbmax,6))

      ioerg= 10

c.....set myid to indicate it is not a parallel computation.
      myid_md= -1

c.....If iolevel is 0, force ifpmd to be 0 as well.
      if( iolevel.eq.0 .and. ifpmd.ne.0 ) then
        ifpmd = 0
      endif

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
      iocntpmd= 0
      
      call setup(nismax,am,acon,fack,dt)
      call boxmat(h,hi,ht,g,gi,gt,vol,sgm)

c.....Set initial temperature if needed
      if( tinit.gt.1d-5 ) then
        call setv(natm,tag,va,h,dt,tinit,nismax,am)
        call rm_trans_motion(natm,tag,va,nismax,am)
      elseif( abs(tinit).le.1d-5 ) then
        va(1:3,1:natm)= 0d0
      endif

c      print *,' hi:'
c      print *, hi(1:3,1)
c      print *, hi(1:3,2)
c      print *, hi(1:3,3)
c      print *,' sgm:'
c      print *, sgm(1:3,1)
c      print *, sgm(1:3,2)
c      print *, sgm(1:3,3)

      call mk_lspr_brute(namax,natm,nbmax,nb,nnmax,tag,ra,rc+rbuf
     &     ,h,hi,sgm,lspr)
c      print *,'lspr:'
c      do ia=1,natm
c        write(6,'(i3,"(",i2,")",100i4)') ia,lspr(0,ia)
c     &       ,(lspr(i,ia),i=1,lspr(0,ia))
c      enddo
      call get_force(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &     ,tcom,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
     &     ,mpi_md_world,myid_md,epi,epot0,nismax,acon,avol
     &     ,cforce)
      epot= epot0
      epotp= epot

c      print *,'aa:'
c      do ia=1,natm
c        write(6,'(i5,3es15.7)') ia, aa(1:3,ia)
c      enddo

      call get_ekin(namax,natm,va,tag,h,nismax,fack,ekin,eki)
      call get_num_dof(natm,tag,fmv,ndof)
c.....only use ifmv(1) for temperature
      temp(1)= ekin /fkb *2d0 /ndof(1)

      if( iolevel.gt.0 ) then
        write(6,'(1x,a)') "initial values:"
        write(6,'(1x,a,es22.14,a)') " kinetic energy  ="
     &       ,ekin," eV"
        write(6,'(1x,a,es22.14,a)') " potential energy="
     &       ,epot0," eV"
        write(6,'(1x,a,es15.7,a)') " temperature     =",temp(1)," K"

        open(ioerg,file="out.erg",status='replace')
        write(ioerg,'(a)') '# istp,etot-epot0,ekin,epot-epot0,temp(K)'
        write(ioerg,'(a,es15.7)') '# epot0 =',epot0
        istp= 0
        write(ioerg,'(i8,3es15.7,f10.2,2es15.7)') istp
     &       ,(ekin+epot0-epot0),ekin,(epot0-epot0)
     &       ,temp(1)
        call flush(ioerg)
      endif

      i_conv = 0

c.....make MD loop starts
      do istp=1,nstp
        va(1:3,1:natm)= va(1:3,1:natm) +aa(1:3,1:natm)
        ra(1:3,1:natm)= ra(1:3,1:natm) +va(1:3,1:natm)
c-------multiply fmv or damping
        do i=1,natm
          l= int(mod(tag(i)*10,10d0))
          va(1:3,i)=va(1:3,i) *fmv(1:3,l)
        enddo

        call mk_lspr_brute(namax,natm,nbmax,nb,nnmax,tag,ra,rc+rbuf
     &       ,h,hi,sgm,lspr)
        call get_force(namax,natm,tag,ra,nnmax,aa,strs,h,hi
     &       ,tcom,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
     &       ,mpi_md_world,myid_md,epi,epot,nismax,acon,avol
     &       ,cforce)
        va(1:3,1:natm)= va(1:3,1:natm) +aa(1:3,1:natm)

        call get_ekin(namax,natm,va,tag,h,nismax,fack,ekin,eki)
        temp(1)= ekin /fkb *2d0 /ndof(1)
        if( iolevel.gt.0 ) then
          write(6,'(a,i8,f10.2)') ' istp,temp = ',istp,temp(1)
          write(ioerg,'(i8,3es15.7,f10.2,2es15.7)') istp
     &         ,(ekin+epot-epot0),ekin,(epot-epot0)
     &         ,temp(1)
          call flush(ioerg)
        endif

c.....temperature control by velocity scaling
c.....note: only 1st element of array is used
        if( trim(ctctl).eq.'Berendsen' ) then
          ifmv = 1
          tfac(ifmv)= 0d0
          if( tfin.gt.0d0 ) then
            ttgt(ifmv) = tinit +(tfin-tinit)*istp/nstp
          endif
          if(ndof(ifmv).le.0 .or. ttgt(ifmv).lt.0d0 ) cycle
          temp(ifmv)= ekin *2d0 /fkb /ndof(ifmv)
          if( abs(ttgt(ifmv)-temp(ifmv))/temp(ifmv).gt.100d0 ) then
            tfac(ifmv)= sqrt(1d0 +dt/trlx*100d0 )
          else
            tfac(ifmv)= sqrt(1d0 +dt/trlx*(ttgt(ifmv)-temp(ifmv))
     &           /temp(ifmv))
          endif
c          write(6,*) ttgt(ifmv),ndof(ifmv),temp(ifmv),tfac(ifmv)
          do i=1,natm
            ifmv= int(mod(tag(i)*10,10d0))
            if( ifmv.ne.1 .or. ttgt(ifmv).lt.0d0 ) cycle
            va(1:3,i)= va(1:3,i) *tfac(ifmv)
          enddo
          
        endif

c.....check convergence criteria if it is dumping simulation
        if( ifdmp.eq.1 .and. epot.lt.epotp .and.
     &       abs((epot-epotp)/epotp).lt.eps_conv ) then
          i_conv = i_conv + 1
          if( i_conv.eq.n_conv ) then
            write(6,'(a,es10.3,a,i3,a)')
     &           'Dumping MD converged since ediff < ',
     &           eps_conv,', ',n_conv,' times'
            write(6,'(a,2es20.10)') '  epot,epotp = ',epot,epotp
            lconverged = .true.
          endif
        else
          epotp = epot
          i_conv = 0
        endif

c-------write the particle positions
        if(ifpmd.eq.1.and.
     &       (mod(istp,noutpmd).eq.0 .or. lconverged))then
c---------decide pmd-file name
          iocntpmd=iocntpmd+1
          write(cnum(1:4),'(i4.4)') iocntpmd
c---------write pmd file for Akira
          call write_pmd_ascii(20,"smd"//cnum,natm,h,hunit
     &         ,tag,ra,va,eki,epi,strs,sorg,dt)
        endif

        if( lconverged ) exit

      enddo ! end of istp=1,nstp

c.....convert accelerations to the real unit
      do ia=1,natm
        is= int(tag(ia))
        at(1:3)= h(1:3,1,0)*aa(1,ia)
     &       +h(1:3,2,0)*aa(2,ia)
     &       +h(1:3,3,0)*aa(3,ia)
        aa(1:3,ia)= at(1:3)/acon(is)
      enddo

c      print *,'aa:'
c      do ia=1,natm
c        write(6,'(i5,3es15.7)') ia, aa(1:3,ia)
c      enddo

      if( iolevel.gt.0 ) then
        write(6,'(1x,a)') "final values:"
        write(6,'(1x,a,es22.14,a)') " kinetic energy  ="
     &       ,ekin," eV"
        write(6,'(1x,a,es22.14,a)') " potential energy="
     &       ,epot," eV"
        write(6,'(1x,a,es15.7,a)') " temperature     ="
     &       ,temp(1)," K"
      endif

      end subroutine smd_core
c=======================================================================
      subroutine setup(nismax,am,acon,fack,dt)
      implicit none
      include "params_unit.h"
      integer,intent(in):: nismax
      real(8),intent(in):: dt,am(nismax)
      real(8),intent(out):: acon(nismax),fack(nismax)

      integer:: i
      real(8):: tmp

c      am(1:nismax)= am(1:nismax) !*ump/ume

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
      subroutine setv(natm,tag,va,h,dt,tinit,nismax,am)
      implicit none
      include 'params_unit.h'
c.....arguments
      integer,intent(in):: natm,nismax
      real(8),intent(in):: tag(natm),h(3,3,0:1),dt,tinit,am(nismax)
      real(8),intent(out):: va(3,natm)
c.....local variables
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
      subroutine rm_trans_motion(natm,tag,va,nismax,am)
      implicit none
c.....arguments
      integer,intent(in):: natm,nismax
      real(8),intent(in):: tag(natm),am(nismax)
      real(8),intent(inout):: va(3,natm)
c.....local variables
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
      subroutine get_num_dof(natm,tag,fmv,ndof)
      implicit none
      integer,intent(in):: natm
      real(8),intent(in):: fmv(3,0:9),tag(natm)
      integer,intent(inout):: ndof(9)
      integer:: i,l,k
      real(8),parameter:: deps= 1d-12

      ndof(1:9)= 0d0
      do i=1,natm
        l= int(mod(tag(i)*10,10d0))
        do k=1,3
          if( abs(fmv(k,l)).lt.0.5d0 ) cycle
          ndof(l)= ndof(l) +1
        enddo
      enddo
      
      write(6,'(a,9i9)') ' ndof=',ndof(1:9)
      return
      end subroutine get_num_dof
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make smd"
c     End:
