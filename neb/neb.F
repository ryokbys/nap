      program neb
c-----------------------------------------------------------------------
c  USAGE:
c    $ ./neb > out.neb
c-----------------------------------------------------------------------
c  You need following files in the working directory:
c    - in.neb
c    - neb###
c    - run_force.sh
c       + in.pmd
c       + 30pmdrun.king
c-----------------------------------------------------------------------
      use variables
      implicit none

      integer:: i,ierr,ionum,istp,is
      integer:: islc,islcp,islcm
      character(len=3):: cslc
      real(8):: te
      real(8),external:: get_etime

      te= get_etime() 
      
      call read_input(10,'in.neb')
      allocate(ra(3,natm,nslc),va(3,natm,nslc),fa(3,natm,nslc)
     &     ,tag(natm,nslc))
cc.....eV/A^2 ---> Hartree/Bohr^2
c      scnst=scnst *ev2hrt /aa2bohr**2
cc.....eV/A ---> Hartree/Bohr
c      feps= feps *ev2hrt /aa2bohr
cc.....A ---> Bohr
c      deps= deps *aa2bohr

c.....Read atom configuration of each slice
      do islc=1,nslc
        ionum= 10+islc
        write(cslc(1:3),'(i3.3)') islc
        call read_atom_conf(ionum,'neb'//cslc,natm
     &       ,h,hunit,ra(1,1,islc),va(1,1,islc),tag(1,islc))
      enddo

      if( trim(cminimize).eq.'velocity_damping' ) then
        call damping_md()
      else if( trim(cminimize).eq.'conjugate_gradient' ) then
        call conjugate_gradient()
      else
        write(6,'(a)') ' [Error] There is no minimization method: '
     &       //cminimize
        stop
      endif

c.....Write final atom configuration
      do islc=1,nslc
        ionum= 10+islc
        write(cslc(1:3),'(i3.3)') islc
        call write_atom_conf(ionum,'neb'//cslc,natm
     &       ,h,hunit,ra(1,1,islc),va(1,1,islc),tag(1,islc))
      enddo
      
      write(6,'(a,f10.1)') ' Total elapsed time=',get_etime()

      end program neb
c=======================================================================
      subroutine read_atom_conf(ionum,cfname,natm,h,hunit,ra,va,tag)
      implicit none
      integer,intent(in):: ionum,natm
      character(len=*),intent(in):: cfname
      real(8),intent(out):: ra(3,natm),tag(natm),h(3,3,0:1),va(3,natm)
     &     ,hunit

      integer:: ia,ib,l,i,itmp

      open(ionum,file=trim(cfname),status='old')
      read(ionum,*) hunit
      read(ionum,*) (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
      h(1:3,1:3,0:1)= h(1:3,1:3,0:1) *hunit
      read(ionum,*) itmp
      if( itmp.ne.natm) then
        write(6,'(a)') ' [Error] itmp.ne.natm !!!'
        write(6,'(a,2i6)') '  itmp,natm=',itmp,natm
        stop
      endif
      do i=1,natm
        read(ionum,*) tag(i),ra(1:3,i),va(1:3,i)
      enddo
      close(ionum)
      
      end subroutine read_atom_conf
c=======================================================================
      subroutine write_atom_conf(ionum,cfname,natm,h,hunit,ra,va,tag)
      implicit none
      integer,intent(in):: natm,ionum
      character(len=*),intent(in):: cfname
      real(8),intent(in):: h(3,3,0:1),ra(3,natm),tag(natm),va(3,natm)
     &     ,hunit

      integer:: ia,ib,l,i

      open(ionum,file=trim(cfname),status='replace')
      write(ionum,'(es23.14e3)') hunit
      write(ionum,'(3es23.14e3)') (((h(ia,ib,l)/hunit,ia=1,3)
     &     ,ib=1,3),l=0,1)
      write(ionum,'(i10)') natm
      do i=1,natm
        write(ionum,'(7es23.14e3,8es12.4)') tag(i),ra(1:3,i)
     &       ,va(1:3,i)
     &       ,0d0,0d0           ! eki(i),epi(i)
     &       ,0d0,0d0,0d0 ,0d0,0d0,0d0 ! strs1,2,3,4,5,6
      enddo
      close(ionum)
      
      end subroutine write_atom_conf
c=======================================================================
      subroutine damping_md()
c
c  Perform MD by means of velocity-Verlet method with velocity damping
c  depending on the cosine of velocity and force.
c
      use variables
      implicit none

      integer:: i,ierr,ionum,istp,is,ifmv
      integer:: islc,islcp,islcm,imin
      integer,external:: ifmv_tag
      logical,external:: lconverged
      real(8),external:: get_etime
      real(8):: epot0,epot,dr(3),di,dmax,vf,vv,ff,cs,dt2,ekin,temp,tmp
     &     ,epotmin,estr0,estr,fi(3),vi(3)
      real(8),allocatable:: rat(:,:,:)
      character:: cslc*3,cstp*4

      allocate(rat(3,natm,nslc))
      rat(1:3,1:natm,1:nslc)= 0d0

      dt= dt /aut
      am(1:nismax)= am(1:nismax) *amu
      dt2= dt/2
      trlx= trlx /aut
      write(6,'(a,es12.4)') '   dt   =',dt
      write(6,'(a,es12.4)') '   dt2  =',dt2
      write(6,'(a,es12.4)') '   trlx =',trlx
      write(6,'(a,es12.4)') '   am(1)=',am(1)
      write(6,'(a,es12.4)') '   am(2)=',am(2)

c.....Set initial temperature
      if( tinit.gt.0d0 ) then
        write(6,'(a,f10.3)') ' Setting initial temperature to ',tinit
        do islc=2,nslc-1
          write(6,'(a,i5)') '  image #=',islc
          call setv(natm,va(1,1,islc),tag(1,islc),tinit)
        enddo
      elseif( tinit.lt.0d0) then
        ! do nothing
        write(6,'(a,f10.3)') ' Not to change initial temperature'
      else
        ! set temperature to zero
        write(6,'(a,f10.3)') ' Setting initial temperature to ',tinit
        va(1:3,1:natm,1:nslc)= 0d0
      endif

c.....First calc of forces
      call get_force(nslc,natm,h,ra,va,tag,fa,epot0,estr0
     &     ,scnst,deps,iclmb)
      write(6,'(a,2es22.14)') ' Initial epot0,estr0= ',epot0,estr0
      
      istp= 0
      call output_energies(istp,nslc)

      do islc=1,nslc
        ionum=10+islc
        write(cslc,'(i3.3)') islc
        write(cstp,'(i4.4)') istp
        call write_atom_conf(ionum,'neb'//cslc//'-'//cstp,natm,h
     &       ,hunit,ra(1,1,islc),va(1,1,islc),tag(1,islc))
      enddo

c.....Check force convergence
      if( lcnvg ) then
        if( lconverged(nslc,natm,h,tag,ra,rat,va,fa,deps,feps,0
     &       ,scnst,iclmb,epot0-epot0) ) return
      else
          call get_ekin(nslc,natm,tag,va,nismax,am,ekin)
          temp= ekin /(natm*(nslc-2))/fkb *2d0/3d0
          write(6,'(a,i6,f10.1,4es22.14,f10.3)')
     &         ' istp,etime,ekin,epot,estr,etot,temp=',
     &         istp,get_etime(),ekin,epot0-epot0,estr0-estr0
     &         ,ekin+epot0-epot0,temp
          imin= 0
          epotmin=epot0+estr0
          rat(1:3,1:natm,1:nslc)= ra(1:3,1:natm,1:nslc)
      endif

c.....Start minimization loop
      do istp=1,nstp

        if( lcsdmp ) then
c.....Get cosine of velocity and force
          vf= 0d0
          vv= 0d0
          ff= 0d0
          do islc=2,nslc-1
            do i=1,natm
              ifmv= ifmv_tag(tag(i,islc))
              fi(1:3)=fa(1:3,i,islc)*fmv(1:3,ifmv)
              vi(1:3)=va(1:3,i,islc)*fmv(1:3,ifmv)
              vf= vf +vi(1)*fi(1) +vi(2)*fi(2) +vi(3)*fi(3)
              vv= vv +vi(1)*vi(1) +vi(2)*vi(2) +vi(3)*vi(3)
              ff= ff +fi(1)*fi(1) +fi(2)*fi(2) +fi(3)*fi(3)
            enddo
          enddo
        else
          vv= 0d0
        endif
        if( vv.lt.1d-15 ) then
          cs= 1d0
        else
          cs= vf/sqrt(vv)/sqrt(ff)
          write(6,'(a,es15.7)') ' cos(theta)=',cs
        endif

c.....First update of velocities
        do islc=2,nslc-1
          do i=1,natm
            is= int(tag(i,islc))
            ifmv= ifmv_tag(tag(i,islc))
            va(1:3,i,islc)=va(1:3,i,islc) *cs
            va(1:3,i,islc)=va(1:3,i,islc) +fa(1:3,i,islc)/am(is) *dt2
            va(1:3,i,islc)=va(1:3,i,islc) *fmv(1:3,ifmv)*dmp
          enddo
        enddo

c.....Update of positions
        if( lcnvg ) rat(1:3,1:natm,1:nslc)= ra(1:3,1:natm,1:nslc)
        do islc=2,nslc-1
          do i=1,natm
            dr(1:3)= va(1:3,i,islc)*dt
            ra(1,i,islc)= ra(1,i,islc) +dr(1)/h(1,1,0)
            ra(2,i,islc)= ra(2,i,islc) +dr(2)/h(2,2,0)
            ra(3,i,islc)= ra(3,i,islc) +dr(3)/h(3,3,0)
            if( ra(1,i,islc).lt.0d0 ) ra(1,i,islc)=ra(1,i,islc) +1d0
            if( ra(2,i,islc).lt.0d0 ) ra(2,i,islc)=ra(2,i,islc) +1d0
            if( ra(3,i,islc).lt.0d0 ) ra(3,i,islc)=ra(3,i,islc) +1d0
            if( ra(1,i,islc).ge.1d0 ) ra(1,i,islc)=ra(1,i,islc) -1d0
            if( ra(2,i,islc).ge.1d0 ) ra(2,i,islc)=ra(2,i,islc) -1d0
            if( ra(3,i,islc).ge.1d0 ) ra(3,i,islc)=ra(3,i,islc) -1d0
          enddo
        enddo

        call get_force(nslc,natm,h,ra,va,tag,fa,epot,estr
     &       ,scnst,deps,iclmb)

c.....Second update of velocities
        do islc=2,nslc-1
          do i=1,natm
            is= int(tag(i,islc))
            va(1:3,i,islc)=va(1:3,i,islc) +fa(1:3,i,islc)/am(is) *dt2
          enddo
        enddo

c-------temperature control by velocity scaling
        if( ltctl ) then
          call get_ekin(nslc,natm,tag,va,nismax,am,ekin)
c---------in case of 3D, <ekin> = 3/2 *NkT
          temp= ekin /(natm*(nslc-2))/fkb *2d0/3d0
c---------Berendsen thermostat
          if( (treq-temp)/temp.gt.100d0 ) then
            tmp= sqrt(1d0 +dt/trlx*100d0 )
          else
            tmp= sqrt(1d0 +dt/trlx*(treq-temp)/temp)
          endif
          va(1:3,1:natm,2:nslc-1)= va(1:3,1:natm,2:nslc-1) *tmp
          write(6,'(a,es12.4)') ' temperature scaling factor=',tmp
        endif

        call output_energies(istp,nslc)

        do islc=2,nslc-1
          ionum=10+islc
          write(cslc,'(i3.3)') islc
          write(cstp,'(i4.4)') istp
          call write_atom_conf(ionum,'neb'//cslc//'-'//cstp,natm,h
     &         ,hunit,ra(1,1,islc),va(1,1,islc),tag(1,islc))
        enddo

c.....Check convergence
        if( lcnvg ) then
          if( lconverged(nslc,natm,h,tag,ra,rat,va,fa,deps,feps,istp
     &         ,scnst,iclmb,epot+estr-epot0-estr0) ) return
c.....Not to check convergence === (NVE or NVT) MD
        else
          if(.not.ltctl) then
            call get_ekin(nslc,natm,tag,va,nismax,am,ekin)
            temp= ekin /(natm*(nslc-2))/fkb *2d0/3d0
          endif
          write(6,'(a,i6,f10.1,4es22.14,f10.3)')
     &         ' istp,etime,ekin,epot,estr,etot,temp=',
     &         istp,get_etime(),ekin,epot-epot0,estr-estr0
     &         ,ekin+epot+estr-epot0-estr0,temp
          if( epot+estr.lt.epotmin ) then
            imin= istp
            epotmin=epot+estr
            rat(1:3,1:natm,1:nslc)= ra(1:3,1:natm,1:nslc)
            write(6,'(a,i6,es22.14)')
     &           ' Writing minimum config.: imin,epotmin=',imin,epotmin
c.....Write minimum atom configuration
            do islc=1,nslc
              ionum= 30+islc
              write(cslc(1:3),'(i3.3)') islc
              call write_atom_conf(ionum,'min'//cslc,natm
     &             ,h,hunit,ra(1,1,islc),va(1,1,islc),tag(1,islc))
            enddo
          endif
        endif

      enddo ! end of istp loop

      end subroutine damping_md
c=======================================================================
      subroutine rel2abs(vec,h)
c
c  Convert vector in relative representation in h-matrix 
c    to absolute value.
c
      implicit none
      real(8),intent(in):: h(3,3)
      real(8),intent(inout):: vec(3)
      real(8):: vt(3)
      
      vt(1:3)= vec(1:3)
      vec(1:3)= h(1:3,1)*vt(1) +h(1:3,2)*vt(2) +h(1:3,3)*vt(3)
      
      end subroutine rel2abs
c=======================================================================
      subroutine get_ekin(nslc,natm,tag,va,nismax,am,ekin)
      implicit none
      integer,intent(in):: nslc,natm,nismax
      real(8),intent(in):: tag(natm,nslc),va(3,natm,nslc),am(nismax)
      real(8),intent(out):: ekin
      
      integer:: i,islc,is

      ekin= 0d0
      do islc=2,nslc-1
        do i=1,natm
          is= int(tag(i,islc))
          ekin= ekin +( va(1,i,islc)*va(1,i,islc)
     &         +va(2,i,islc)*va(2,i,islc)
     &         +va(3,i,islc)*va(3,i,islc) )*am(is)
        enddo
      enddo
      ekin=ekin/2

      end subroutine get_ekin
c=======================================================================
      subroutine get_force(nslc,natm,h,ra,va,tag,fa
     &     ,ergp,ergs,scnst,deps,iclmb)
c  
c  Calculate forces of the current atomic configuration
c  using an external program.
c  
      use variables, only : cmethod, ev2hrt, aa2bohr, hunit
      implicit none
      integer,intent(in):: nslc,natm,iclmb
      real(8),intent(in):: h(3,3),ra(3,natm,nslc),tag(natm,nslc)
     &     ,va(3,natm,nslc),scnst,deps
      real(8),intent(out):: fa(3,natm,nslc),ergp,ergs

      real(8),parameter:: pi = 3.14159265358979d0

      integer:: i,j,itmp,islc,islcp,islcm,ionum,islcini,islcfin
     &     ,nfix
      real(8):: xpm(3),rpm,ti(3),tti(3,3),ttim(3,3),xp(3),xm(3)
     &     ,rxp,rxm,cosphi,fphi,fs(3),ft(3),ergt,dr(3),tmp
      real(8):: ftau,fstau
      character:: cslc*3, carg*128
      real(8),external:: anint_rk

      logical,save:: l1st=.true.
      real(8),allocatable,save:: vxpm(:,:),vtau(:,:),vfs(:,:)

      if( l1st ) then
        islcini= 1
        islcfin= nslc
        write(6,'(a,i6)') 'natm=',natm
        allocate(vxpm(3,natm),vtau(3,natm),vfs(3,natm))
      else
        islcini= 2
        islcfin= nslc-1
      endif
      
c.....Write atom config file
      do islc=islcini,islcfin
        ionum=10+islc
        write(cslc,'(i3.3)') islc
        call write_atom_conf(ionum,'neb'//cslc,natm,h
     &       ,hunit,ra(1,1,islc),va(1,1,islc),tag(1,islc))
      enddo
c.....write(6,'(a)') ' wrote files...'

c.....Calc forces using an external program
      write(carg,'(1x,2i3)') islcini,islcfin
c.....write(6,'(a)') ' ./run_force.sh '//trim(carg)
      call system('./run_force.sh '//trim(carg))

c.....Read forces given by the external program
      ergp= 0d0
      fa(1:3,1:natm,1:nslc)= 0d0
      do islc=1,nslc
        ionum=20+islc
        write(cslc,'(i3.3)') islc
        open(ionum,file='frc'//cslc,status='old')
        read(ionum,*) itmp
        if( itmp.ne.natm ) stop ' [Error] itmp.ne.natm !!!'
        do i=1,natm
          read(ionum,*) fa(1:3,i,islc)
c.....In case of vasp, we need to convert from eV/A to Ht/Bohr
c          fa(1:3,i,islc)=fa(1:3,i,islc) *ev2hrt /aa2bohr
          fa(1:3,i,islc)=fa(1:3,i,islc)
        enddo
        close(ionum)
c.....Energies, too
        open(ionum+10,file='erg'//cslc,status='old')
        read(ionum+10,*) ergt
c.....In case of vasp, we need to convert from eV to Ht
        ergt= ergt *ev2hrt
        close(ionum+10)
        ergp=ergp +ergt
      enddo

c.....Another definition of the energy
c      ergt= 0d0
c      ergp = 0d0
c      do islc=1,nslc-1
c        do i=1,natm
c          dr(1:3)= ra(1:3,i,islc+1) -ra(1:3,i,islc)
c          dr(1)= dr(1) -anint_rk(dr(1))
c          dr(2)= dr(2) -anint_rk(dr(2))
c          dr(3)= dr(3) -anint_rk(dr(3))
c          dr(1:3)= h(1:3,1)*dr(1) +h(1:3,2)*dr(2) +h(1:3,3)*dr(3)
c          ft(1:3)= (fa(1:3,i,islc)+fa(1:3,i,islc+1))/2
c          ergt= ergt -ft(1)*dr(1)
c     &         -ft(2)*dr(2)
c     &         -ft(3)*dr(3)
c        enddo
c        ergp= ergp +ergt
c      enddo

      ergs= 0d0
      nfix= 0
c.....Nudged elastic band
      if( cmethod(1:1).eq.'n' .or. cmethod(1:1).eq.'N' ) then
        do islc=1,nslc
          islcp=islc+1
          islcm=islc-1
          if( islcp.gt.nslc ) islcp=nslc
          if( islcm.lt.1    ) islcm=1
c.....NEB forces using tangential spring forces
          rpm= 0d0
          do i=1,natm
            vxpm(1:3,i)= ra(1:3,i,islcp)-ra(1:3,i,islcm)
            vxpm(1,i)= vxpm(1,i) -anint_rk(vxpm(1,i))
            vxpm(2,i)= vxpm(2,i) -anint_rk(vxpm(2,i))
            vxpm(3,i)= vxpm(3,i) -anint_rk(vxpm(3,i))
            call rel2abs(vxpm(1,i),h)
            rpm=rpm +sqrt(vxpm(1,i)**2 +vxpm(2,i)**2 +vxpm(3,i)**2)
          enddo
c            write(6,'(a,2i8,es12.4)') ' islc,i,rpm/deps='
c     &           ,islc,i,rpm/deps
c.....stop if there is no change between images
          if( rpm.lt.deps ) then
            write(6,'(a)') ' [Error] something wrong with images!'
            write(6,'(a)') '    rpm.lt.eps'
            write(6,'(a)') ' which means there is no change between'//
     &           ' neighboring images...'
            stop
          endif
c.....tangential vector (tau) and scalar product (f*tau)
          ftau= 0d0
          do i=1,natm
            vtau(1:3,i)= vxpm(1:3,i)/rpm
            ftau= ftau +fa(1,i,islc)*vtau(1,i) +fa(2,i,islc)*vtau(2,i)
     &           +fa(3,i,islc)*vtau(3,i)
          enddo

c.....switching function
          cosphi= 0d0
          rxp=0d0
          rxm=0d0
          do i=1,natm
            xp(1:3)= ra(1:3,i,islcp) -ra(1:3,i,islc)
            xp(1)=xp(1) -anint_rk(xp(1))
            xp(2)=xp(2) -anint_rk(xp(2))
            xp(3)=xp(3) -anint_rk(xp(3))
            xm(1:3)= ra(1:3,i,islcm) -ra(1:3,i,islc)
            xm(1)=xm(1) -anint_rk(xm(1))
            xm(2)=xm(2) -anint_rk(xm(2))
            xm(3)=xm(3) -anint_rk(xm(3))
            call rel2abs(xp,h)
            call rel2abs(xm,h)
            rxp=rxp +sqrt(xp(1)**2 +xp(2)**2 +xp(3)**2)
            rxm=rxm +sqrt(xm(1)**2 +xm(2)**2 +xm(3)**2)
            cosphi=cosphi +xp(1)*xm(1) +xp(2)*xm(2) +xp(3)*xm(3)
            vfs(1:3,i)= scnst*xp(1:3) +scnst*xm(1:3)
            do j=1,3
              ergs= ergs +0.25d0*scnst*(xp(j)**2 +xm(j)**2)
            enddo
          enddo
          if( islc.ne.1 .and. islc.ne.nslc ) then
            cosphi=cosphi/rxp/rxm
            fphi= 0.5d0*(1d0 +cos(pi*cosphi))
          else
            fphi= 1d0
          endif
c.....Spring force
          fstau= 0d0
          do i=1,natm
            fstau=fstau +vfs(1,i)*vtau(1,i) +vfs(2,i)*vtau(2,i)
     &           +vfs(3,i)*vtau(3,i)
          enddo
c.....Forces on atoms
          if( islc.eq.iclmb ) then
            do i=1,natm
              fa(1:3,i,islc)= fa(1:3,i,islc)
     &             -2d0*ftau*vtau(1:3,i)
            enddo
          else
c.....Forces normal to tangential vector
            do i=1,natm
              fa(1:3,i,islc)= fa(1:3,i,islc)
     &             -ftau*vtau(1:3,i) +fstau*vtau(1:3,i)
c.....Correction using switching function
              fa(1:3,i,islc)= fa(1:3,i,islc)
     &             +fphi*(vfs(1:3,i)-fstau*vtau(1:3,i))
            enddo
          endif
        enddo ! end of islc loop

      else if( cmethod(1:1).eq.'a' .or. cmethod(1:1).eq.'A' ) then
c.....atom individual version of NEB
        do islc=1,nslc
          islcp=islc+1
          islcm=islc-1
          if( islcp.gt.nslc ) islcp=nslc
          if( islcm.lt.1    ) islcm=1
c.....NEB forces using tangential spring forces
          do i=1,natm
            xpm(1:3)= ra(1:3,i,islcp)-ra(1:3,i,islcm)
            xpm(1)= xpm(1) -anint_rk(xpm(1))
            xpm(2)= xpm(2) -anint_rk(xpm(2))
            xpm(3)= xpm(3) -anint_rk(xpm(3))
            call rel2abs(xpm,h)
            rpm= sqrt(xpm(1)**2 +xpm(2)**2 +xpm(3)**2)
c            write(6,'(a,2i8,es12.4)') ' islc,i,rpm/deps='
c     &           ,islc,i,rpm/deps
c.....Skip if atom-i has almost no displacement
            if( rpm.lt.deps ) then
              fa(1:3,i,islc)= 0d0
              nfix= nfix +1
              cycle
            endif
c.....Tangential vector, tensor
            ti(1:3)= xpm(1:3)/rpm
            tti(1,1:3)= (/ ti(1)*ti(1), ti(1)*ti(2), ti(1)*ti(3) /)
            tti(2,1:3)= (/ ti(2)*ti(1), ti(2)*ti(2), ti(2)*ti(3) /)
            tti(3,1:3)= (/ ti(3)*ti(1), ti(3)*ti(2), ti(3)*ti(3) /)
            ttim(1:3,1:3)= -tti(1:3,1:3)
            ttim(1,1)= ttim(1,1) +1d0
            ttim(2,2)= ttim(2,2) +1d0
            ttim(3,3)= ttim(3,3) +1d0
            xp(1:3)= ra(1:3,i,islcp) -ra(1:3,i,islc)
            xp(1)=xp(1) -anint_rk(xp(1))
            xp(2)=xp(2) -anint_rk(xp(2))
            xp(3)=xp(3) -anint_rk(xp(3))
            xm(1:3)= ra(1:3,i,islcm) -ra(1:3,i,islc)
            xm(1)=xm(1) -anint_rk(xm(1))
            xm(2)=xm(2) -anint_rk(xm(2))
            xm(3)=xm(3) -anint_rk(xm(3))
            call rel2abs(xp,h)
            call rel2abs(xm,h)
c.....Switching function
            rxp= sqrt(xp(1)**2 +xp(2)**2 +xp(3)**2)
            rxm= sqrt(xm(1)**2 +xm(2)**2 +xm(3)**2)
            if( islc.ne.1 .and. islc.ne.nslc ) then
              cosphi= (xp(1)*xm(1) +xp(2)*xm(2) +xp(3)*xm(3)) /rxp /rxm
              fphi= 0.5d0 *(1d0 +cos(pi*cosphi))
            else
              fphi= 1d0
            endif
c.....Spring force
            fs(1:3)= scnst*xp(1:3) +scnst*xm(1:3)
            do j=1,3
              ergs= ergs +0.25d0*scnst*(xp(j)**2 +xm(j)**2)
            enddo
c.....Forces on atoms
            ft(1:3)= fa(1:3,i,islc)
            if( islc.eq.iclmb ) then
              fa(1:3,i,islc)= ft(1:3)
     &             -2d0*tti(1:3,1)*ft(1)
     &             -2d0*tti(1:3,2)*ft(2)
     &             -2d0*tti(1:3,3)*ft(3)
            else
c.....Forces normal to tangential vector
              fa(1:3,i,islc)= ft(1)*ttim(1:3,1)
     &             +ft(2)*ttim(1:3,2)
     &             +ft(3)*ttim(1:3,3)
     &             +fs(1)*tti(1:3,1)
     &             +fs(2)*tti(1:3,2)
     &             +fs(3)*tti(1:3,3)
c.....Correction using switching function
              fa(1:3,i,islc)= fa(1:3,i,islc) +fphi*( fs(1:3)
     &             -fs(1)*tti(1:3,1)
     &             -fs(2)*tti(1:3,2)
     &             -fs(3)*tti(1:3,3) )
            endif
          enddo
        enddo ! end of islc loop

c.....Plain elastic method
      else if( cmethod(1:1).eq.'p' .or. cmethod(1:1).eq.'P' ) then
        do islc=1,nslc
          islcp=islc+1
          islcm=islc-1
          if( islcp.gt.nslc ) islcp=nslc
          if( islcm.lt.1    ) islcm=1
c.....NEB forces using tangential spring forces
          do i=1,natm
            xpm(1:3)= ra(1:3,i,islcp)-ra(1:3,i,islcm)
            xpm(1)= xpm(1) -anint(xpm(1))
            xpm(2)= xpm(2) -anint(xpm(2))
            xpm(3)= xpm(3) -anint(xpm(3))
            call rel2abs(xpm,h)
            rpm= sqrt(xpm(1)**2 +xpm(2)**2 +xpm(3)**2)
c.....Skip if atom-i has almost no displacement
            if( rpm.lt.deps ) then
              fa(1:3,i,islc)= 0d0
              cycle
            endif
c.....Tangential vector, tensor
            ti(1:3)= xpm(1:3)/rpm
            tti(1,1:3)= (/ ti(1)*ti(1), ti(1)*ti(2), ti(1)*ti(3) /)
            tti(2,1:3)= (/ ti(2)*ti(1), ti(2)*ti(2), ti(2)*ti(3) /)
            tti(3,1:3)= (/ ti(3)*ti(1), ti(3)*ti(2), ti(3)*ti(3) /)
            xp(1:3)= ra(1:3,i,islcp) -ra(1:3,i,islc)
            xp(1)=xp(1) -anint_rk(xp(1))
            xp(2)=xp(2) -anint_rk(xp(2))
            xp(3)=xp(3) -anint_rk(xp(3))
            xm(1:3)= ra(1:3,i,islcm) -ra(1:3,i,islc)
            xm(1)=xm(1) -anint_rk(xm(1))
            xm(2)=xm(2) -anint_rk(xm(2))
            xm(3)=xm(3) -anint_rk(xm(3))
            call rel2abs(xp,h)
            call rel2abs(xm,h)
c.....Spring force
            fs(1:3)= scnst*xp(1:3) +scnst*xm(1:3)
            do j=1,3
              ergs= ergs +0.25d0*scnst*(xp(j)**2 +xm(j)**2)
            enddo
c.....Forces on atoms
            ft(1:3)= fa(1:3,i,islc)
            if( islc.eq.iclmb ) then
              fa(1:3,i,islc)= ft(1:3)
     &             -2d0*tti(1:3,1)*ft(1)
     &             -2d0*tti(1:3,2)*ft(2)
     &             -2d0*tti(1:3,3)*ft(3)
            else
              fa(1:3,i,islc)=ft(1:3) +fs(1:3)
            endif
          enddo
        enddo ! end of islc loop
      endif

c      write(6,'(a,2i8)') ' natm,nfix=',natm,nfix

 9999 l1st=.false.
      end subroutine get_force
c=======================================================================
      function lconverged(nslc,natm,h,tag,ra,rat,va,fa,deps,feps
     &     ,istp,scnst,iclmb,erg)
      use variables, only: cminimize,fmv
      implicit none
      integer,intent(in):: nslc,natm,istp,iclmb
      real(8),intent(in):: tag(natm,nslc),ra(3,natm,nslc),erg
     &     ,rat(3,natm,nslc),feps,deps,h(3,3),scnst
     &     ,va(3,natm,nslc)
      real(8),intent(inout)::fa(3,natm,nslc) 

      logical:: lconverged

      integer:: i,ierr,islc,islcmax,imax,ifmv
      real(8):: fmax,fi,ft(3),dmax,di,dr(3),epot,estr,dimax
      real(8),external:: get_etime,anint_rk
      integer,external:: ifmv_tag

      lconverged= .false.

      dmax= 0d0
      fmax= 0d0
c      write(6,'(a,2i6)') ' nslc,natm=',nslc,natm
      do islc=2,nslc-1
        do i=1,natm
          dr(1:3)= ra(1:3,i,islc)-rat(1:3,i,islc)
          dr(1)= dr(1) -anint_rk(dr(1))
          dr(2)= dr(2) -anint_rk(dr(2))
          dr(3)= dr(3) -anint_rk(dr(3))
          dr(1:3)= h(1:3,1)*dr(1)
     &         +h(1:3,2)*dr(2)
     &         +h(1:3,3)*dr(3)
c          dr(1:3)= h(1:3,1)*(ra(1,i,islc)-rat(1,i,islc))
c     &         +h(1:3,2)*(ra(2,i,islc)-rat(2,i,islc))
c     &         +h(1:3,3)*(ra(3,i,islc)-rat(3,i,islc))
          di= sqrt(dr(1)**2 +dr(2)**2 +dr(3)**2)
          dmax= max(di,dmax)
          ifmv= ifmv_tag(tag(i,islc))
          ft(1:3)= fa(1:3,i,islc)*fmv(1:3,ifmv)
          fi= sqrt(ft(1)**2 +ft(2)**2 +ft(3)**2)
c          write(6,'(a,2i4,4es12.4)') 'islc,i,fi,ft(1:3)='
c     &         ,islc,i,fi,ft(1:3)
          if( fi .gt. fmax ) then
            imax= i
            islcmax= islc
            dimax= di
          endif
          fmax= max(fi,fmax)
        enddo
      enddo
      write(6,'(a,i6,f10.1,5es22.14)')
     &     ' istp,etime,dE,fmax,feps,dmax,deps='
     &     ,istp,get_etime(),erg,fmax,feps,dmax,deps
c      write(6,'(a,2i6,es12.4)') ' imax,islcmax,dimax='
c     &     ,imax,islcmax,dimax
c      if( dmax.lt.deps .or. fmax.lt.feps ) then
      if( fmax.lt.feps ) then
c        if( fmax.lt.feps ) then
          write(6,'(a)') ' Convergence achieved for force.'
          write(6,'(a,2es15.7)') '  fmax, feps=',fmax,feps
c        else
c          write(6,'(a)') ' Convergence achieved for displacement.'
c          write(6,'(a,2es15.7)') '  dmax, deps=',dmax,deps
c        endif

c.....Get force of the final configuration in case of CG
        if( cminimize(1:1).eq.'c' ) then
          write(6,'(a)') ' Computing forces of final configuration...'
          call get_force(nslc,natm,h,ra,va,tag,fa,epot,estr
     &         ,scnst,deps,iclmb)
        endif
        lconverged= .true.
      endif

      return
      end function lconverged
c=======================================================================
      subroutine conjugate_gradient()
      use variables
      implicit none
      integer:: istp,islc,i,ifmv
      integer,external:: ifmv_tag
      logical,external:: lconverged
      real(8):: dmax,rt(3),ft(3),di,epot,epot0,estr,estr0,fnrm,fnrm0
     &     ,erg,erg0
      real(8),allocatable:: rat(:,:,:),fa0(:,:,:)

      allocate(rat(3,natm,nslc),fa0(3,natm,nslc))

      fnrm= 1d0
      fa(1:3,1:natm,1:nslc)= 0d0
      va(1:3,1:natm,1:nslc)= 0d0
      fa0(1:3,1:natm,1:nslc)= fa(1:3,1:natm,1:nslc)

c.....First calc of forces
      call get_force(nslc,natm,h,ra,va,tag,fa,epot0,estr0
     &     ,scnst,deps,iclmb)
      erg0= epot0+estr0
      erg= erg0

      call output_energies(0,nslc)

      do istp=1,nstp
        fnrm0= fnrm
        rat(1:3,1:natm,1:nslc)= ra(1:3,1:natm,1:nslc)

c.....Force norm
        fnrm= 0d0
        do islc=2,nslc-1
          do i=1,natm
            ft(1:3)= fa(1:3,i,islc)
            fnrm= fnrm +ft(1)**2 +ft(2)**2 +ft(3)**2
          enddo
        enddo

c.....Fletcher-Reeves CG
        fa(1:3,1:natm,2:nslc-1)= fa(1:3,1:natm,2:nslc-1)
     &       +fnrm/fnrm0 *fa0(1:3,1:natm,2:nslc-1)

c.....Multiplying fmv
        do islc=2,nslc-1
          do i=1,natm
            ifmv= ifmv_tag(tag(i,islc))
            fa(1:3,i,islc)=fa(1:3,i,islc)*fmv(1:3,ifmv)
          enddo
        enddo

c.....Line minimization in fa direction to get updated ra
        call line_minimize(nslc,natm,h,tag,ra,va,fa,erg)
        
        call output_energies(istp,nslc)

c.....Check force convergence
        if( lconverged(nslc,natm,h,tag,ra,rat,va,fa,deps,feps,istp
     &       ,scnst,iclmb,erg-erg0) ) return

        fa0(1:3,1:natm,1:nslc)= fa(1:3,1:natm,1:nslc)
        call get_force(nslc,natm,h,ra,va,tag,fa,epot,estr
     &       ,scnst,deps,iclmb)
        erg= epot +estr
      enddo

      write(6,'(a,i6)') ' Convergence not achieved within ',nstp
      return
      end subroutine conjugate_gradient
c=======================================================================
      subroutine line_minimize(nslc,natm,h,tag,ra,va,fa,erg)
      use variables, only : deps, aa2bohr, scnst, iclmb
      implicit none
      integer,intent(in):: nslc,natm
      real(8),intent(in):: h(3,3),fa(3,natm,nslc),tag(natm,nslc)
     &     ,va(3,natm,nslc)
      real(8),intent(inout):: ra(3,natm,nslc),erg
      
      real(8),parameter:: xfac= 2d0
      integer,parameter:: itmax= 100

      integer:: islc,i,it
      real(8):: fi,fmax,ft(3),ds,s0,s1,s2,s3,s4,erg0,erg1,erg2,erg3,erg4
     &     ,tau,d,hmax,epot,estr

      logical,save:: l1st=.true.
      real(8),allocatable,save:: ra1(:,:,:),ra2(:,:,:)
     &     ,ra3(:,:,:),ra4(:,:,:),fat(:,:,:)

      if( l1st ) then
        allocate(ra1(3,natm,nslc),ra2(3,natm,nslc)
     &       ,ra3(3,natm,nslc),ra4(3,natm,nslc)
     &       ,fat(3,natm,nslc))
        ra1(1:3,1:natm,1)= ra(1:3,1:natm,1)
        ra2(1:3,1:natm,1)= ra(1:3,1:natm,1)
        ra3(1:3,1:natm,1)= ra(1:3,1:natm,1)
        ra4(1:3,1:natm,1)= ra(1:3,1:natm,1)
        ra1(1:3,1:natm,nslc)= ra(1:3,1:natm,nslc)
        ra2(1:3,1:natm,nslc)= ra(1:3,1:natm,nslc)
        ra3(1:3,1:natm,nslc)= ra(1:3,1:natm,nslc)
        ra4(1:3,1:natm,nslc)= ra(1:3,1:natm,nslc)
        fat(1:3,1:natm,1)= fa(1:3,1:natm,1)
        fat(1:3,1:natm,nslc)= fa(1:3,1:natm,nslc)
        l1st= .false.
      endif

c.....Get max force
      fmax=0d0
      do islc=2,nslc-1
        do i=1,natm
          ft(1:3)=fa(1:3,i,islc)
          fmax= max(fmax,abs(ft(1)))
          fmax= max(fmax,abs(ft(2)))
          fmax= max(fmax,abs(ft(3)))
c          fi= sqrt(ft(1)**2 +ft(2)**2 +ft(3)**2)
c          fmax= max(fi,fmax)
        enddo
      enddo
c.....First step distance: 0.1 AA
c      ds= 0.01d0 *aa2bohr  /fmax
      ds= 0.01d0 /fmax
      hmax= max(h(1,1),h(2,2))
      hmax= max(hmax,h(3,3))
      ds= ds/hmax

c      write(6,'(a,3es15.7)') ' fmax,ds,fmax*ds=',fmax,ds,fmax*ds
c      write(6,'(a,es15.7)') ' fmax*ds*h(1,1)=',fmax*ds*h(1,1)
c      write(6,'(a,es15.7)') ' fmax*ds*h(2,2)=',fmax*ds*h(2,2)
c      write(6,'(a,es15.7)') ' fmax*ds*h(3,3)=',fmax*ds*h(3,3)

c.....Determine range; within [x1,x2], there must be minimum.
      s0= 0d0
      s1= 0d0
      s2= ds
      erg1= erg
      do it=1,itmax
        ra2(1:3,1:natm,2:nslc-1)= ra(1:3,1:natm,2:nslc-1)
     &       +s2*fa(1:3,1:natm,2:nslc-1)
        call get_force(nslc,natm,h,ra2,va,tag,fat,epot,estr
     &       ,scnst,deps,iclmb)
        erg2=epot+estr
        write(6,'(a,i6,2es22.14)') '   it,erg1,erg2=',it,erg1,erg2
        if( erg2.gt.erg1 ) then
          if( it.gt.1 ) then
            erg1=erg0
            s1= s0
          endif
          exit
        else
          s0= s1
          s1= s2
          s2= s2*xfac
          erg0= erg1
          erg1= erg2
        endif
      enddo
c      stop
      if( it.gt.itmax ) then 
        write(6,'(a)') ' [Error] Could not find miniimum within a range'
        write(6,'(a,es12.4)') '  ds= ',ds
        stop
      endif

c.....Golden section search
      ra1(1:3,1:natm,2:nslc-1)= ra(1:3,1:natm,2:nslc-1)
     &     +s1*fa(1:3,1:natm,2:nslc-1)
      tau= (1d0 +sqrt(5d0))/2
      do it=1,itmax
        write(6,'(a,i6,2es22.14)') '   LINMIN steps=',it,erg1,erg2
        d= s2-s1
c.....Convergence checking in bohr unit
        if( d*fmax.lt.deps ) exit
        s3= (tau-1d0)/tau *d
        s4= 1d0/tau *d
        ra3(1:3,1:natm,2:nslc-1)= ra1(1:3,1:natm,2:nslc-1)
     &       +s3*fa(1:3,1:natm,2:nslc-1)
        ra4(1:3,1:natm,2:nslc-1)= ra1(1:3,1:natm,2:nslc-1)
     &       +s4*fa(1:3,1:natm,2:nslc-1)
        call get_force(nslc,natm,h,ra3,va,tag,fat,epot,estr
     &       ,scnst,deps,iclmb)
        erg3=epot+estr
        call get_force(nslc,natm,h,ra4,va,tag,fat,epot,estr
     &       ,scnst,deps,iclmb)
        erg4=epot+estr
c.....Update x1 or x4 dependent on relation between erg3 and erg4
        if( erg3.ge.erg4 ) then
          s1= s3
          ra1(1:3,1:natm,2:nslc-1)= ra3(1:3,1:natm,2:nslc-1)
          erg1= erg3
        else
          s2= s4
          ra2(1:3,1:natm,2:nslc-1)= ra4(1:3,1:natm,2:nslc-1)
          erg2= erg4
        endif
      enddo

      if( erg1.lt.erg2 ) then
        ra(1:3,1:natm,2:nslc-1)= ra1(1:3,1:natm,2:nslc-1)
        erg= erg1
      else
        ra(1:3,1:natm,2:nslc-1)= ra2(1:3,1:natm,2:nslc-1)
        erg= erg2
      endif

      return
      end subroutine line_minimize
c=======================================================================
      subroutine output_energies(istp,nslc)
      integer,intent(in):: istp,nslc

      integer:: islc
      character:: cfouterg*11,cslc*3

c.....Output energies of every slices every steps
      cfouterg='out.erg0000'
      write(cfouterg(8:11),'(i4.4)') istp
      call system('if [ -f '//cfouterg//' ]; then rm '
     &     //cfouterg//'; fi')
      do islc=1,nslc
        cslc='000'
        write(cslc,'(i3.3)') islc
        call system('cat erg'//cslc//' >> '//cfouterg)
      enddo

      end subroutine output_energies
c=======================================================================
      function ifmv_tag(tag)
      implicit none
      real(8),intent(in):: tag
      integer:: ifmv_tag
      
      ifmv_tag= int((tag-int(tag))*10)
      return
      end function ifmv_tag
c=======================================================================
      function anint_rk(x)
      real(8),intent(in):: x
      real(8):: anint_rk
      anint_rk= 0d0
      if( x.gt.0.5d0 ) then
        anint_rk= 1.d0
      else if( x.lt.-0.5d0 ) then
        anint_rk= -1.d0
      endif
      return
      end function anint_rk
c=======================================================================
      function get_etime()
c-----------------------------------------------------------------------
c  Time measurement routine which provides elapsed time
c  from the 1st call.
c-----------------------------------------------------------------------
      implicit none
      real(8):: get_etime
      integer:: i1,irate,imax
      logical,save:: l1st=.true.
      integer,save:: i0

      if( l1st ) then
        call system_clock(i0)
        get_etime= 0d0
        l1st=.false.
      else
        call system_clock(i1,irate,imax)
        if( i1 < i0 ) then
          get_etime= i1+imax -i0
        else
          get_etime= i1 -i0
        endif
        get_etime= get_etime/irate
      endif

      return
      end function get_etime
c=======================================================================
      subroutine setv(natm,va,tag,tempini)
      use variables, only: fkb,amu
      implicit none
      integer,intent(in):: natm
      real(8),intent(in):: tempini,tag(natm)
      real(8),intent(out):: va(3,natm)

      integer:: i,l,is
      real(8):: dseed,sumvx,sumvy,sumvz,rnd1,rnd2,tmp,facv,am
      real(8),parameter:: pi = 3.14159265358979d0
      real(8),parameter:: am_fe = 55.847d0 *amu

      facv=dsqrt(2d0*tempini*fkb/am_fe)
      am= am_fe

c-----velocities in Maxwell-Boltzmann distribution
      dseed=12345
      do i=1,natm
        is= int(tag(i))
        do l=1,3
          call myrnd(rnd1,dseed)
          call myrnd(rnd2,dseed)
          va(l,i)=facv*dsqrt(-dlog(rnd1))*dcos(2d0*pi*rnd2)
        enddo
      enddo
c-----set center of mass motion to zero
      sumvx=0d0
      sumvy=0d0
      sumvz=0d0
      do i=1,natm
        sumvx=sumvx+va(1,i)
        sumvy=sumvy+va(2,i)
        sumvz=sumvz+va(3,i)
      enddo
      do i=1,natm
        va(1,i)=va(1,i)-sumvx/dble(natm)
        va(2,i)=va(2,i)-sumvy/dble(natm)
        va(3,i)=va(3,i)-sumvz/dble(natm)
      enddo

      tmp=0d0
      do i=1,natm
        is= int(tag(i))
        tmp= tmp +0.5d0*am*(va(1,i)**2 +va(2,i)**2 +va(3,i)**2)
      enddo
      write(6,'(a,es12.4)') " ekin=",tmp
      write(6,'(a,es12.4)') " temp.=",tmp*2d0/3d0/fkb/natm

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
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make neb"
c     End:
