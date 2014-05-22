      program qmcl
c-----------------------------------------------------------------------
c  Hybrid QM-CL program which uses external MD and DFT codes
c    to calculate classical and quantum forces, respectively,
c    via a shell script.
c-----------------------------------------------------------------------
      use variables
      implicit none
      include 'params_unit.h'
      integer:: i
      real(8):: t0
      real(8),external:: get_etime

      t0= get_etime()

      call initialize()

      call read_input(10,'in.qmcl')
c.....Multiply damping factor to avoid multiplication of every step
      if( ldmp ) then
        fmv(1:3,0:9)= fmv(1:3,0:9) *dmp
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

      call setup()

c.....Check input values
      if( rclst.lt.rcqm ) then
        write(6,'(a)') ' [qmcl] rclst.lt.rcqm !!!'
        write(6,'(a,2es12.4)') '  rclst,rcqm=',rclst,rcqm
        stop
      endif

c.....check force_type
      if( trim(cfrct).ne.'qm' .and. trim(cfrct).ne.'QM' .and.
     &     trim(cfrct).ne.'cl' .and. trim(cfrct).ne.'CL' .and.
     &     trim(cfrct).ne.'qmcl' .and. trim(cfrct).ne.'QMCL' ) then
        write(6,'(a)') ' [Error] Force type should be qm/cl/qmcl !!!'
        write(6,'(a)') '   Currentforce type = '//trim(cfrct)
        stop
      endif
        
c.....check run_type
      if( crtype(1:1).ne.'t' .and. crtype(1:1).ne.'T' .and.
     &     crtype(1:1).ne.'m' .and. crtype(1:1).ne.'M' .and.
     &     crtype(1:1).ne.'c' .and. crtype(1:1).ne.'C' ) then
        write(6,'(a)') ' [Error] Force type should be test/md/cg !!!'
        write(6,'(a)') '   Currentforce type = '//trim(crtype)
        stop
      endif
        
c      call read_pmd_ascii(20,"qmclref"
c     &     ,namax,natm,h0,tag,ra0,va,eki,epi,strs)
      call read_pmd_ascii(20,"qmcl000"
     &     ,namax,natm,h,hunit,tag,ra,va,eki,epi,strs)
c.....Some conversions
      do i=1,natm
c.....Species of atom-i
        nis= max(nis,int(tag(i)))
        if( nis.gt.nismax ) then
          write(6,'(a)') ' [qmcl] nis.gt.nismax !'
          stop
        endif
      enddo
      write(6,'(a,i8)') ' Number of atoms   =',natm
      write(6,'(a,i8)') ' Number of species =',nis

      if( namax.lt.natm+nbmax ) then
        stop ' [qmcl] namax.lt.natm+nbmax !!!'
      endif

c.....set initial temperature if needed
      if( tinit.gt.1d-1 ) then
        call setv()
        call rm_trans_motion()
      elseif( abs(tinit).le.1d-1 ) then
        va(1:3,1:natm)= 0d0
      endif

      call boxmat(h,hi,vol,sgm)
c      write(6,'(a)') ' h0-mat:'
c      write(6,'(3es12.4)') h0(1,1:3,0)
c      write(6,'(3es12.4)') h0(2,1:3,0)
c      write(6,'(3es12.4)') h0(3,1:3,0)

      if( lconst ) then
        call read_frcav(20,'frcav000')
        call read_vnorm(21,cfnorm)
        call set_vnorm(3*natm,vn)
      endif

      if( crtype(1:1).eq.'t' .or. crtype(1:1).eq.'T' ) then
        call velocity_verlet()
      else if( crtype(1:1).eq.'m' .or. crtype(1:1).eq.'M' ) then
        call velocity_verlet()
      else if( crtype(1:1).eq.'c' .or. crtype(1:1).eq.'C' ) then
        call conjugate_gradient()
      endif

      end program qmcl
c=======================================================================
      subroutine boxmat(h,hi,vol,sgm)
c-----------------------------------------------------------------------
c  setup matrices of MD-box
c    H:   MD-box matrix
c    HI:  inverse MD-box matrix
c    SGM: cofactor matrix
c-----------------------------------------------------------------------
      implicit none
      real(8),intent(in):: h(3,3,0:1)
      real(8),intent(out):: vol,sgm(3,3),hi(3,3)

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

      write(6,'(a)') ' h-mat:'
      write(6,'(3es12.4)') h(1,1:3,0)
      write(6,'(3es12.4)') h(2,1:3,0)
      write(6,'(3es12.4)') h(3,1:3,0)
      write(6,'(a)') ' hi-mat:'
      write(6,'(3es12.4)') hi(1,1:3)
      write(6,'(3es12.4)') hi(2,1:3)
      write(6,'(3es12.4)') hi(3,1:3)

      return
      end subroutine boxmat
c=======================================================================
      subroutine initialize()
      use variables
      implicit none
      include 'params_unit.h'
      integer:: i
      real(8),parameter:: am_al = 26.9815d0

c.....Set all mass to Al
      am(1:nismax)= am_al

      end subroutine initialize
c=======================================================================
      subroutine setup()
      use variables
      implicit none
      include 'params_unit.h'
      integer:: i

c.....prefactors for normalized accelerration
      write(6,'(a)') ' acon:'
      do i=1,nismax
        acon(i)= ev2j/ang2m *0.5d0 *(dt*fs2s)**2 /(am(i)*amu2kg)*m2ang
        write(6,'(i5,es12.4)') i,acon(i)
      enddo

c.....prefactors for kinetic energy in eV
      write(6,'(a)') ' fack:'
      do i=1,nismax
        fack(i)= 0.5d0*(am(i)*amu2kg)/(dt*fs2s)**2 *ang2m**2 *j2ev
        write(6,'(i5,es12.4)') i,fack(i)
      enddo

      end subroutine setup
c=======================================================================
      subroutine velocity_verlet()
      use variables
      implicit none
      include 'params_unit.h'
      integer:: i,is,l
      real(8):: dr(3),tmp,cs,ff,vv,vf,fi(3),vi(3)
      character:: cnum*3
      real(8),external:: get_etime,sprod
      logical,external:: lconverged

      call get_force(natm,h,ra,tag,va,fa,epi,eki,epot0,0)

      call get_ekin(natm,tag,va,eki,ekin,nismax,fack)
c-----in case of 3D, <ekin> = 3/2 *NkT
      temp= ekin /natm/fkb *2d0/3d0

      write(6,'(a)') " Initial values: [Hartree], [K]"
      write(6,'(a,es15.7)') " Boltzmann factor=",fkb
      write(6,'(a,es15.7)') " Kinetic energy  =",ekin
      write(6,'(a,es15.7)') " Potential energy=",epot0
      write(6,'(a,es15.7)') " Temperature     =",temp

      write(cnum(1:3),'(i3.3)') 0
      call write_force(21,cnum,natm,h,tag,fa,epot0)

c.....initialize the counter for output
      iocntpmd=0
      iocnterg=0
c.....output initial configuration
      if( lpmd ) then
        write(cnum(1:3),'(i3.3)') iocntpmd
        call write_pmd_ascii(30,"qmcl"//cnum
     &       ,natm,h,hunit,tag,ra,va,eki,epi,stn)
      endif

c.....Energy output
      istp= 0
      open(ioerg,file=cferg,status='replace')
      call output_energy(ioerg,0,ekin,epot,epot0,temp)

c.....Check convergence
      if( lconv ) then
        if( lconverged(natm,h,tag,fa,feps,istp,fmv,epot0-epot0) ) return
      endif

c.....MD loop with velocity-Verlet......................................
      do istp=1,nstp

        if( lcsdmp ) then
c.....Get cosine of velocity and force
          vf= 0d0
          vv= 0d0
          ff= 0d0
          do i=1,natm
            l= ifmv_tag(tag(i))
            fi(1:3)=fa(1:3,i)*fmv(1:3,l)
            vi(1:3)=va(1:3,i)*fmv(1:3,l)
            vf= vf +vi(1)*fi(1) +vi(2)*fi(2) +vi(3)*fi(3)
            vv= vv +vi(1)*vi(1) +vi(2)*vi(2) +vi(3)*vi(3)
            ff= ff +fi(1)*fi(1) +fi(2)*fi(2) +fi(3)*fi(3)
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


c.....1st kick of velocities
        do i=1,natm
          is= int(tag(i))
          l= int(mod(tag(i)*10,10d0))
          va(1:3,i)= va(1:3,i) *cs
c          va(1:3,i)= va(1:3,i) +0.5d0*dt*fa(1:3,i)/am(is)
          va(1:3,i)= va(1:3,i) +fa(1:3,i)*acon(is)
          va(1:3,i)= va(1:3,i) *fmv(1:3,l)
        enddo

c.....apply constraint to velocities if needed
c.....note: this is not applicable to non-othogonal systems
        if( lconst ) then
c          tmp= sprod(3*natm,va,vn)
          do i=1,natm
            tmp= sprod(3,va(1,i),vn(1,i))
            va(1:3,i)= va(1:3,i) -tmp*vn(1:3,i)
          enddo
        endif

c.....Update positions
        do i=1,natm
          dr(1:3)= va(1:3,i) !*dt
          ra(1:3,i)= ra(1:3,i) +hi(1:3,1)*dr(1)
     &         +hi(1:3,2)*dr(2) +hi(1:3,3)*dr(3)
        enddo

c.....PBC
        do i=1,natm
          if( ra(1,i).lt.0d0 ) ra(1,i)=ra(1,i) +1d0
          if( ra(2,i).lt.0d0 ) ra(2,i)=ra(2,i) +1d0
          if( ra(3,i).lt.0d0 ) ra(3,i)=ra(3,i) +1d0
          if( ra(1,i).ge.1d0 ) ra(1,i)=ra(1,i) -1d0
          if( ra(2,i).ge.1d0 ) ra(2,i)=ra(2,i) -1d0
          if( ra(3,i).ge.1d0 ) ra(3,i)=ra(3,i) -1d0
        enddo

        call get_force(natm,h,ra,tag,va,fa,epi,eki,epot,istp)

c.....accumulate force for averaging
        if( lconst ) then
c          tmp= sprod(3*natm,fa,vn)
          do i=1,natm
            facc(i)=facc(i) +sprod(3,fa(1,i),vn(1,i))
          enddo
        endif

c.....2nd kick of velocities
        do i=1,natm
          is= int(tag(i))
c          va(1:3,i)= va(1:3,i) +0.5d0*dt*fa(1:3,i)/am(is)
          va(1:3,i)= va(1:3,i) +fa(1:3,i)*acon(is)
        enddo

c.....Calc kinetic energy
        call get_ekin(natm,tag,va,eki,ekin,nismax,fack)

c.....temperature control by velocity scaling
        if( ltctl ) then
c.....in case of 3D, <ekin> = 3/2 *NkT
          temp= ekin /natm/fkb *2d0/3d0
cc.....in case of 2D, <ekin> = 2/2 *NkT
c          temp= ekin /natm/fkb *2d0/2d0
c.....Berendsen thermostat
          if( (treq-temp)/temp.gt.100d0 ) then
            tmp= sqrt(1d0 +dt/trlx*100d0 )
          else
            tmp= sqrt(1d0 +dt/trlx*(treq-temp)/temp)
          endif
          va(1:3,1:natm)= va(1:3,1:natm) *tmp
        endif

c.....Output energies
        if(mod(istp,nouterg).eq.0) then
          iocnterg=iocnterg+1
c.....In case of 3D, <ekin> = 3/2 *NkT
          temp= ekin /natm/fkb *2d0/3d0
cc---------in case of 2D, <ekin> = 2/2 *NkT
c          temp= ekin /natm/fkb *2d0/2d0
          call output_energy(ioerg,istp,ekin,epot,epot0,temp)
        endif

c.....write the particle positions
        if( lpmd .and. mod(istp,noutpmd).eq.0 )then
c.....decide qmcl-file name
          iocntpmd=iocntpmd+1
          write(cnum(1:3),'(i3.3)') iocntpmd
cc.....Calc atomic strain
c          call get_atomic_strain(namax,natm,h,ra,stn)
c.....write qmcl file for Akira
          call write_pmd_ascii(30,"qmcl"//cnum
     &         ,natm,h,hunit,tag,ra,va,eki,epi,stn)
          if( lconst ) then
            call write_frcav(21,'frcav'//cnum)
          endif

        endif

c.....Check convergence
        if( lconv ) then
          if( lconverged(natm,h,tag,fa,feps,istp,fmv,epot-epot0) )
     &         return
        endif

      enddo ! istp
      
      end subroutine velocity_verlet
c=======================================================================
      subroutine get_ekin(natm,tag,va,eki,ekin,nismax,fack)
      use variables,only: am
      implicit none
      integer,intent(in):: natm,nismax
      real(8),intent(in):: tag(natm),va(3,natm),fack(nismax)
      real(8),intent(out):: eki(natm),ekin

      integer:: i,is
      
      eki(1:natm)= 0d0
      ekin= 0d0

      do i=1,natm
        is= int(tag(i))
c        eki(i)= 0.5d0*am(is) *(va(1,i)**2 +va(2,i)**2 +va(3,i)**2)
        eki(i)= (va(1,i)**2 +va(2,i)**2 +va(3,i)**2) *fack(is)
        ekin= ekin +eki(i)
      enddo
      
      return
      end subroutine get_ekin
c=======================================================================
      subroutine get_force(natm,h,ra,tag,va,fa,epi,eki,erg,istp)
      use variables, only:nacmax,nac,hc,rac,tagc,xc,rclst,rcqm,wvac
     &     ,idclst,fac,naqm,idqm,strs,nis,is_tag,rcut,cfrct,lterm,nhc
     &     ,nhcmax,rhc,crtype,hunit
      implicit none
      integer,intent(in):: natm,istp
      real(8),intent(in):: h(3,3,0:1),ra(3,natm)
     &     ,epi(natm),eki(natm)
      real(8),intent(inout):: va(3,natm)
      real(8),intent(out):: fa(3,natm),erg,tag(natm)

      integer:: ic,i,iqm,is,ncp,ix,iy,iz,m
      real(8):: ht(3,3,0:1),dcpx,dcpy,dcpz,ergcl,ergqm
      character:: cstp*5
      logical,save:: l1st= .true.
      integer,save:: ncpx,ncpy,ncpz,natmt
      real(8),allocatable,save:: rat(:,:),vat(:,:),fat(:,:),tagt(:)
     &     ,ekit(:),epit(:),strst(:,:,:)
      logical:: l_exists_vasp_dir,l_exists_pmd_dir

      if( l1st ) then
cc.....Create a cluster only at the first call
c        call get_qm_cluster_id(xc,rclst**2,rcqm**2,h,natm,ra,tag,nacmax
c     &       ,nac,idclst,naqm,idqm)
c.....read cluster data from files out.idclst and out.idqm
        call read_qm_cluster(nacmax,nac,naqm,idclst,idqm)
c.....Compute the num of copies along the direction shorter than rcut
        ncpx=1
        ncpy=1
        ncpz=1
        if( h(1,1,0).lt.2*rcut ) then
          ncpx= 1+int(2*rcut/h(1,1,0))
        elseif( h(2,2,0).lt.2*rcut ) then
          ncpy= 1+int(2*rcut/h(2,2,0))
        elseif( h(3,3,0).lt.2*rcut ) then
          ncpz= 1+int(2*rcut/h(3,3,0))
        endif
        write(6,'(a,3i5)') ' ncpx,ncpy,ncpz=',ncpx,ncpy,ncpz
        natmt= natm*ncpx*ncpy*ncpz
        write(6,'(a,2i10)') ' natm,natmt    =',natm,natmt
        ht(1,1,0)= h(1,1,0)*ncpx
        ht(2,2,0)= h(2,2,0)*ncpy
        ht(3,3,0)= h(3,3,0)*ncpz
        allocate(rat(3,natmt),vat(3,natmt),fat(3,natmt),tagt(natmt)
     &       ,ekit(natmt),epit(natmt),strst(3,3,natmt))
        m=0
        do iz=1,ncpz
          dcpz= dble(iz-1)/ncpz
          do iy=1,ncpy
            dcpy= dble(iy-1)/ncpy
            do ix=1,ncpx
              dcpx= dble(ix-1)/ncpx
              do i=1,natm
                m=m+1
                tagt(m)= tag(i)
              enddo
            enddo
          enddo
        enddo
cc.....Set tags of QM atoms
c        do iqm=1,naqm
c          ic= idqm(iqm)
c          i= idclst(ic)
c          is= is_tag(tag(i))
c          tag(i)= tag(i) -is
c          if( is.gt.nis ) is=is-nis
c          tag(i)= tag(i) +(nis+is)
c        enddo
        l1st= .false.
      endif

c.....Create qm configuration
      call create_qm_conf(h,natm,ra,tag,nac,nacmax,idclst
     &     ,hc,rac,tagc,wvac,xc)

c.....If the system is thinner than 2*rcut, enlarge it for pmd
      m=0
      do iz=1,ncpz
        dcpz= dble(iz-1)/ncpz
        do iy=1,ncpy
          dcpy= dble(iy-1)/ncpy
          do ix=1,ncpx
            dcpx= dble(ix-1)/ncpx
            do i=1,natm
              m=m+1
              rat(1,m)= ra(1,i)/ncpx +dcpx
              rat(2,m)= ra(2,i)/ncpy +dcpy
              rat(3,m)= ra(3,i)/ncpz +dcpz
              vat(1:3,m)= va(1:3,i)
            enddo
          enddo
        enddo
      enddo

c.....Write atom config file for pmd
      call write_pmd_ascii(20,'CL_input',natmt,ht,hunit,tagt,rat,vat
     &     ,ekit,epit,strst)
c.....Write atom config file for vasp
      call write_vasp(21,'QM_input',nac+nhc,hc,rac,tagc)

      if( crtype(1:1).eq.'t' .or. crtype(1:1).eq.'T' ) then
        write(6,'(a)') ' CL_input and QM_input have been written.'
        write(6,'(a)') ' Test run finished.'
        stop
      endif

c.....Run shellscript which executes external program to get forces
      write(cstp,'(i5)') istp
      ergcl= 0d0
      ergqm= 0d0
      if( trim(cfrct).eq.'qmcl' .or. trim(cfrct).eq.'QMCL' ) then
        call system('./run_CL+QM.sh '//cstp)
        call read_force_CL(22,'frc.CL',natmt,fat)
        call read_force_QM(23,'frc.QM',nac+nhc,fac)
        call read_erg(24,'erg.CL',ergcl)
        call read_erg(25,'erg.QM',ergqm)
        fa(1:3,1:natm)= fat(1:3,1:natm)
c.....Replace forces of 'QM'-atoms
        do iqm=1,naqm
          ic= idqm(iqm)
          i= idclst(ic)
          fa(1:3,i)= fac(1:3,ic)
        enddo

      else if( trim(cfrct).eq.'qm' .or. trim(cfrct).eq.'QM' ) then
        call system('./run_QM.sh '//cstp)
        call read_force_QM(23,'frc.QM',nac+nhc,fac)
        call read_erg(25,'erg.QM',ergqm)
        fa(1:3,1:natm)= fat(1:3,1:natm)
c.....Replace forces of 'QM'-atoms
        do iqm=1,naqm
          ic= idqm(iqm)
          i= idclst(ic)
          fa(1:3,i)= fac(1:3,ic)
        enddo

      else if( trim(cfrct).eq.'cl' .or. trim(cfrct).eq.'CL' ) then
        call system('./run_CL.sh '//cstp)
        call read_force_CL(22,'frc.CL',natmt,fat)
        call read_erg(24,'erg.CL',ergcl)
        fa(1:3,1:natm)= fat(1:3,1:natm)
      endif


c.....Calc potential energy
      erg= ergcl +ergqm

      return
      end subroutine get_force
c=======================================================================
      subroutine get_qm_cluster_id(xc,rc2,rcqm2,h,natm,ra,tag
     &     ,nacmax,nac,idclst,naqm,idqm)
      use variables, only: ra0
      implicit none
      integer,intent(in):: natm,nacmax
      integer,intent(out):: nac,naqm,idclst(nacmax),idqm(nacmax)
      real(8),intent(in):: h(3,3),ra(3,natm),tag(natm),rc2,xc(3),rcqm2

      integer:: i
      real(8):: xt(3),xr(3),r2

c.....From the center position and radius, create the cluster
      nac= 0
      naqm= 0
      do i=1,natm
        xt(1:3)= ra(1:3,i)-xc(1:3) -anint(ra(1:3,i)-xc(1:3))
c        xt(1:3)= ra0(1:3,i)-xc(1:3) -anint(ra0(1:3,i)-xc(1:3))
        xr(1:3)= h(1:3,1)*xt(1) +h(1:3,2)*xt(2) +h(1:3,3)*xt(3)
c.....If xc(?) is minus, no length needed in ?-direction
        if( xc(1).lt.0d0 ) xr(1)= 0d0
        if( xc(2).lt.0d0 ) xr(2)= 0d0
        if( xc(3).lt.0d0 ) xr(3)= 0d0
        r2= xr(1)**2 +xr(2)**2 +xr(3)**2
        if( r2.ge.rc2 ) cycle
        nac=nac +1
        idclst(nac)= i
c        write(6,'(a,2i6)') ' nac,idclst=',nac,idclst(nac)
        if( r2.ge.rcqm2 ) cycle
        naqm=naqm +1
        idqm(naqm)= nac
c        write(6,'(a,2i6)') ' naqm,idqm(naqm)=',naqm,idqm(naqm)
      enddo
      write(6,'(a,i6)') ' num of cluster atoms=',nac
      write(6,'(a,i6)') ' num of QM atoms     =',naqm

      open(90,file='out.idclst')
      write(90,'(i6)') nac
      do i=1,nac
        write(90,'(i8)') idclst(i)
      enddo
      close(90)

      open(91,file='out.idqm')
      write(91,'(i6)') naqm
      do i=1,naqm
        write(91,'(i8)') idqm(i)
      enddo
      close(91)
      
      return
      end subroutine get_qm_cluster_id
c=======================================================================
      subroutine read_qm_cluster(nacmax,nac,naqm,idclst,idqm)
      implicit none
      integer,intent(in):: nacmax
      integer,intent(out):: nac,naqm,idclst(nacmax),idqm(nacmax)

      integer:: i
      
      open(90,file='out.idclst',status='old')
      read(90,*) nac
      if( nac.gt.nacmax ) stop ' [Error] nac.gt.nacmax !!!'
      do i=1,nac
        read(90,*) idclst(i)
      enddo
      close(90)

      open(91,file='out.idqm',status='old')
      read(91,*) naqm
      do i=1,naqm
        read(91,*) idqm(i)
      enddo
      close(91)

      end subroutine read_qm_cluster
c=======================================================================
      subroutine create_qm_conf(h,natm,ra,tag,nac,nacmax,idclst
     &     ,hc,rac,tagc,wvac,xc)
      use variables,only: lterm,rhc,nhc,nhcmax
      implicit none
      integer,intent(in):: natm,nac,nacmax,idclst(nacmax)
      real(8),intent(in):: ra(3,natm),tag(natm),h(3,3),wvac(3),xc(3)
      real(8),intent(out):: rac(3,nacmax),tagc(nacmax),hc(3,3)
      
      integer:: i,ic
      real(8):: xmin(3),xmax(3),xt(3)
      logical,save:: l1st=.true.

      if( l1st ) then
        if(  (xc(1).lt.0d0 .and. nint(wvac(1)*1d+2).ne.0) .or.
     &       (xc(2).lt.0d0 .and. nint(wvac(2)*1d+2).ne.0) .or.
     &       (xc(3).lt.0d0 .and. nint(wvac(3)*1d+2).ne.0)  ) then
          write(6,'(a)') ' [create_qm_cell] xc(?).lt.0d0 '
     &         //'.and. wvac(?).ne.0d0 !!!'
          stop
        endif
      endif
      
c.....Copy atom positions of 'QM'-atoms
      do ic=1,nac
        i= idclst(ic)
        rac(1:3,ic)= ra(1:3,i)
        tagc(ic)= tag(i)
      enddo

c.....put termination H atoms if needed
      nhc=0
      if( lterm ) then
        call add_termination_H()
      endif

c.....Determine cell size: hc
      xmin(1:3)= 1d10
      xmax(1:3)= -1d0
      do ic=1,nac+nhc
        xt(1:3)= rac(1:3,ic)
        if( xc(1).ge.0d0 ) xt(1)= xt(1) -anint(xt(1)-xc(1))
        if( xc(2).ge.0d0 ) xt(2)= xt(2) -anint(xt(2)-xc(2))
        if( xc(3).ge.0d0 ) xt(3)= xt(3) -anint(xt(3)-xc(3))
        rac(1:3,ic)= h(1:3,1)*xt(1) +h(1:3,2)*xt(2) +h(1:3,3)*xt(3)
        xmin(1:3)= min(xmin(1:3),rac(1:3,ic))
        xmax(1:3)= max(xmax(1:3),rac(1:3,ic))
      enddo
c      write(6,'(a,3es12.4)') ' xmin(1:3)=',xmin(1:3)
c      write(6,'(a,3es12.4)') ' xmax(1:3)=',xmax(1:3)
c.....Assume orthorhombic cell, and add 2*(vacuum width)
      if( xc(1).lt.0d0 ) then
        xmin(1)= 0d0
        xmax(1)= h(1,1)
      endif
      if( xc(2).lt.0d0 ) then
        xmin(2)= 0d0
        xmax(2)= h(2,2)
      endif
      if( xc(3).lt.0d0 ) then
        xmin(3)= 0d0
        xmax(3)= h(3,3)
      endif
c      write(6,'(a,3es12.4)') ' xmin(1:3)=',xmin(1:3)
c      write(6,'(a,3es12.4)') ' xmax(1:3)=',xmax(1:3)
      if( l1st ) then
        hc(1:3,1:3)= 0d0
        hc(1,1)= xmax(1)-xmin(1) +2d0*wvac(1)
        hc(2,2)= xmax(2)-xmin(2) +2d0*wvac(2)
        hc(3,3)= xmax(3)-xmin(3) +2d0*wvac(3)
      endif
c      write(6,'(a,3es12.4)') ' hc(1,1:3)=',hc(1,1:3)
c      write(6,'(a,3es12.4)') ' hc(2,1:3)=',hc(2,1:3)
c      write(6,'(a,3es12.4)') ' hc(3,1:3)=',hc(3,1:3)
c      write(6,'(a)') 'rac:'
c      do ic=1,nac
c        write(6,'(i6,3es12.4)') ic,rac(1:3,ic)
c      enddo

c.....Convert h-rep. to hc-rep.
      do ic=1,nac+nhc
c.....Shift
        rac(1:3,ic)= rac(1:3,ic) -xmin(1:3) +wvac(1:3)
c.....Normalize
        rac(1,ic)= rac(1,ic)/hc(1,1)
        rac(2,ic)= rac(2,ic)/hc(2,2)
        rac(3,ic)= rac(3,ic)/hc(3,3)
      enddo
c      write(6,'(a)') 'rac:'
c      do ic=1,nac
c        write(6,'(i6,3es12.4)') ic,rac(1:3,ic)
c      enddo

      l1st=.false.
      return
      end subroutine create_qm_conf
c=======================================================================
      subroutine add_termination_H()
      use variables
      implicit none
      include './params_unit.h'
c.....cutoff radius for Si
      real(8),parameter:: rcsi  = 3.00d0
c.....the ratio used to shorten Si-Si bond to get Si-H bond
      real(8),parameter:: bratio= 0.7d0
c.....Ideal Si-H bond length
      real(8),parameter:: sih_ratio  = 1.5d0/2.45d0

      integer:: ic,i,j,k,l,m,n,maxn,maxi
      real(8):: xi(3),xj(3),xij(3),xnj(3),rij,yij(3)
      logical,external:: exists_in_cluster
      real(8),external:: pbc
      logical,save:: l1st= .true.
      
      if( l1st ) then
        write(6,'(a,f10.4)') ' Si-H/Si-Si ratio=',sih_ratio

c.....Make lspr for looking at neighbors using rcut from user input.
c.....lspr should include only 1st N.N.
        ra0(1:3,1:natm)= ra(1:3,1:natm)
        call mklspr(rcsi)

c.....Check num of 1st N.N. atoms
        maxn= 0
        maxi= 0
        do ic=1,nac
          i= idclst(ic)
          if( lspr(0,i).gt.maxn ) then
            maxn= lspr(0,i)
            maxi= i
          endif
        enddo
c      if( maxn.gt.4 ) then ! depends on the system 
        if( maxn.gt.8 ) then
          write(6,'(a,i5)') ' [Error] Num of neighbors is too large:'
     &         ,maxn
          write(6,'(a)') '   You should decrease cutoff_radius.'
          write(6,'(a,es12.4)') '   Currnet cutoff_radius = ',rcsi
          write(6,'(a,i6)') '   maxi=',maxi
          write(6,'(a,20i6)') '   lspr=',lspr(1:lspr(0,maxi),maxi)
          write(6,'(3x,i6,a,3es12.4)') maxi,': ',ra0(1:3,maxi)
          do m=1,lspr(0,maxi)
            j= lspr(m,maxi)
            write(6,'(3x,i6,a,3es12.4)') j,': ',ra0(1:3,j)
          enddo
          stop
        endif
      endif
      
c.....For all cluster atoms, search neighbors that are not included
c.....in the cluster. And add H atom on the line connecting these
c.....atoms included and not included.
      nhc= 0
      do ic=1,nac
        i= idclst(ic)
        xi(1:3)= ra0(1:3,i)
        do m=1,lspr(0,i)
          j=lspr(m,i)
          xj(1:3)= ra0(1:3,j)
          if( exists_in_cluster(j,nac,idclst,xj,ra0,h,namax) ) cycle
          xij(1:3)= xj(1:3)-xi(1:3)
          yij(1:3)= h(1:3,1,0)*xij(1)
     &         +h(1:3,2,0)*xij(2)
     &         +h(1:3,3,0)*xij(3)
          rij= sqrt(yij(1)**2 +yij(2)**2 +yij(3)**2)
          xnj(1:3)= ra0(1:3,i) +xij(1:3)*sih_ratio
          xnj(1)= pbc(xnj(1))
          xnj(2)= pbc(xnj(2))
          xnj(3)= pbc(xnj(3))
          nhc= nhc +1
          if( nac+nhc.gt.nacmax ) stop ' [Error] nac+nhc.gt.nacmax !!!'
          rac(1:3,nac+nhc)= xnj(1:3)
          tagc(nac+nhc)= 1d0*(nis+1) +1d-1*1
        enddo
      enddo
      if( l1st ) write(6,'(a,i6)') ' num of termination H atoms=',nhc

c.....Check if these atoms apart enough from each other.
      do i=1,nhc-1
        xi(1:3)= rac(1:3,nac+i)
        do j=i+1,nhc
          xj(1:3)= rac(1:3,nac+j)
          xij(1:3)= xj(1:3)-xi(1:3) -anint(xj(1:3)-xi(1:3))
          xnj(1:3)= h(1:3,1,0)*xij(1)
     &         +h(1:3,2,0)*xij(2)
     &         +h(1:3,3,0)*xij(3)
          rij= sqrt(xnj(1)**2 +xnj(2)**2 +xnj(3)**2)
          if( rij.lt.0.5d0 ) then
            write(6,'(a)') ' [Error] H-H is too short ???'
            write(6,'(a,i3,a,i3,a,es12.4,a)') '  H-H:',i,'-',j
     &           ,'=',rij,' Ang'
            write(6,'(a,i6,3es12.4)') ' i,rac(1:3,nac+i)='
     &           ,i,rac(1:3,nac+i)
            write(6,'(a,i6,3es12.4)') ' j,rac(1:3,nac+j)='
     &           ,j,rac(1:3,nac+j)
            stop
          endif
        enddo
      enddo

      l1st= .false.
      end subroutine add_termination_H
c=======================================================================
      function exists_in_cluster(ia,nac,idclst,x,ra,h,namax)
      implicit none
      integer,intent(in):: ia,nac,idclst(nac),namax
      real(8),intent(in):: x(3),ra(3,namax),h(3,3,0:1)
      logical:: exists_in_cluster

      integer:: ic,i
      real(8):: xi(3),xj(3),xt(3)
      real(8),external:: pbc

      exists_in_cluster= .false.

      do ic=1,nac
        i= idclst(ic)
        xi(1:3)= h(1:3,1,0)*ra(1,i)
     &       +h(1:3,2,0)*ra(2,i)
     &       +h(1:3,3,0)*ra(3,i)
        xt(1)= pbc(x(1))
        xt(2)= pbc(x(2))
        xt(3)= pbc(x(3))
        xj(1:3)= h(1:3,1,0)*xt(1)
     &       +h(1:3,2,0)*xt(2)
     &       +h(1:3,3,0)*xt(3)
        if(  int((xj(1)-xi(1))*100).eq.0 .and.
     &       int((xj(2)-xi(2))*100).eq.0 .and.
     &       int((xj(3)-xi(3))*100).eq.0 ) then
          exists_in_cluster= .true.
          return
        endif
      enddo
      return
      end function exists_in_cluster
c=======================================================================
      subroutine write_vasp(ionum,cfname,nac,hc,hunit,rac,tagc)
c
c  Write POSCAR file for VASP run
c
      implicit none
      include 'params_unit.h'
      integer,intent(in):: ionum,nac
      real(8),intent(in):: hunit,hc(3,3),rac(3,nac),tagc(nac)
      character(len=*),intent(in):: cfname

      integer:: i,j,is
      logical:: lexist
      integer,parameter:: nsmax= 10
      logical,save:: l1st=.true.
      integer,save:: ns,lsis(nsmax),nis(nsmax)

      if( l1st ) then
c.....check ns
        lsis(1:nsmax)= 0
        nis(1:nsmax)= 0
        ns=0
        do i=1,nac
          is= int(tagc(i))
          lexist=.false.
          do j=1,ns
            if( is.eq.lsis(j) ) then
              lexist=.true.
              nis(j)=nis(j)+1
              exit
            endif
          enddo
          if( .not.lexist ) then
            ns=ns+1
            if( ns.gt.nsmax ) stop ' Error: ns.gt.nsmax !!!'
            lsis(ns)= is
            nis(ns)= 1
          endif
        enddo

        write(6,'(a,i5)') ' num of species in QM cluster=',ns
        do i=1,ns
          write(6,'("is=",i3,":",i5)') lsis(i),nis(i)
        enddo
        l1st=.false.
      endif
      
      open(ionum,file=trim(cfname),status='replace')
      write(ionum,'(a)') ' QM atoms created by QMCL' ! comment line
      write(ionum,'(f10.2)') hunit
      write(ionum,'(3f12.5)') hc(1,1:3)/hunit
      write(ionum,'(3f12.5)') hc(2,1:3)/hunit
      write(ionum,'(3f12.5)') hc(3,1:3)/hunit
      write(ionum,'(10i4)') (nis(i),i=1,ns)
c      write(ionum,'(2i6)') nac-1,1
      write(ionum,'(a)') 'direct'
      do j=1,ns
        is= lsis(j)
        do i=1,nac
          if( int(tagc(i)).ne.is ) cycle
          write(ionum,'(3f12.7)') rac(1:3,i)
        enddo
      enddo
      close(ionum)

      return
      end subroutine write_vasp
c=======================================================================
      subroutine read_force_CL(ionum,cfname,natm,fa)
      implicit none
      integer,intent(in):: ionum,natm
      character(len=*),intent(in):: cfname
      real(8),intent(out):: fa(3,natm)

      integer:: i,itmp

      open(ionum,file=trim(cfname),status='old')
      read(ionum,*) itmp
      if( itmp.ne.natm ) stop ' [read_force_CL] itmp.ne.natm !!!'
      do i=1,natm
        read(ionum,*) fa(1:3,i)
      enddo
      close(ionum)
      
      end subroutine read_force_CL
c=======================================================================
      subroutine read_erg(ionum,cfname,erg)
      implicit none
      integer,intent(in):: ionum
      character(len=*),intent(in):: cfname
      real(8),intent(out):: erg

      open(ionum,file=trim(cfname),status='old')
      read(ionum,*) erg
      close(ionum)
      
      end subroutine read_erg
c=======================================================================
      subroutine read_force_QM(ionum,cfname,natm,fa)
      implicit none
      include 'params_unit.h'
      integer,intent(in):: ionum,natm
      character(len=*),intent(in):: cfname
      real(8),intent(out):: fa(3,natm)

      integer:: i,itmp

      open(ionum,file=trim(cfname),status='old')
      read(ionum,*) itmp
      if( itmp.ne.natm ) stop ' [read_force_QM] itmp.ne.natm !!!'
      do i=1,natm
        read(ionum,*) fa(1:3,i)
      enddo
      close(ionum)

      end subroutine read_force_QM
c=======================================================================
      subroutine get_atomic_strain(namax,na,h,ra,stn)
c-----------------------------------------------------------------------
c  Calculate strain tensor according to
c    F.Shimizu, et al., Mater. Trans. 48 (2007) pp.2923-2927.
c-----------------------------------------------------------------------
      use variables,only: ra0,h0,nnmax,lspr,mklspr
      implicit none
      include 'params_unit.h'
      integer,intent(in):: namax,na
      real(8),intent(in):: ra(3,namax),h(3,3)
      real(8),intent(out):: stn(3,3,namax)

      integer:: i,j,l,ixyz,jxyz
      real(8):: x,y,z,rij2
      logical,save:: l1st=.true.
      real(8),save:: rc2
      real(8),allocatable,save:: xi(:),xi0(:),xij(:),xij0(:),ddt(:,:)
     &     ,dd0i(:,:)

      if( l1st ) then
        allocate(xi(3),xi0(3),xij(3),xij0(3),ddt(3,3),dd0i(3,3))
c.....rc2: cutoff for 1st nearest neighbors in FCC (Al)
        rc2= 4.0425d0 *(1d0 +1d0/sqrt(2d0))/2
        rc2= rc2*rc2
cc.....rc2: cutoff for 1st nearest neighbors in BCC-Fe
c        rc2= 2.8553d0 *(1d0 +sqrt(3d0)/2) /2
c        rc2= rc2*rc2
        l1st=.false.
      endif

      call mklspr(sqrt(rc2))

c-----initialize strain tensor
      stn(1:3,1:3,1:na)= 0d0

      do i=1,na
        xi(1:3)=ra(1:3,i)
        xi0(1:3)=ra0(1:3,i)
        ddt(1:3,1:3)= 0d0
        dd0i(1:3,1:3)= 0d0
        do l=1,lspr(0,i)
          j=lspr(l,i)
          xij(1)= (ra(1,j)-xi(1)-anint(ra(1,j)-xi(1)))*h(1,1)
          xij(2)= (ra(2,j)-xi(2)-anint(ra(2,j)-xi(2)))*h(2,2)
          xij(3)= (ra(3,j)-xi(3)-anint(ra(3,j)-xi(3)))*h(3,3)
          rij2= xij(1)**2 +xij(2)**2 +xij(3)**2
c---------only 1st nearest neighbors
          if(rij2.gt.rc2) cycle
          xij0(1)= (ra0(1,j)-xi0(1))*h0(1,1,0)
          xij0(2)= (ra0(2,j)-xi0(2))*h0(2,2,0)
          xij0(3)= (ra0(3,j)-xi0(3))*h0(3,3,0)
          do jxyz=1,3
            do ixyz=1,3
              stn(ixyz,jxyz,i)= stn(ixyz,jxyz,i) +xij0(ixyz)*xij(jxyz)
              ddt(ixyz,jxyz)= ddt(ixyz,jxyz) +xij0(ixyz)*xij0(jxyz)
            enddo
          enddo
        enddo
c        if( i.eq.1 ) then
c          write(6,'(a,9es11.3)') ' stn=',stn(1:3,1:3,i)
c          write(6,'(a,9es11.3)') ' ddt=',ddt(1:3,1:3)
c        endif
c-------inverse of ddt ---> dd0i
c        call ludc_inv(2,ddt(1:2,1:2),dd0i(1:2,1:2))
        call ludc_inv(3,ddt,dd0i)
c        if( i.eq.1 ) then
c          write(6,'(a,9es11.3)') ' dd0i=',dd0i(1:3,1:3)
c        endif
c-------J =[dd0i]*dd =[dd0i]*stn
        ddt(1:3,1:3)= stn(1:3,1:3,i)
        do jxyz=1,3
          do ixyz=1,3
            stn(ixyz,jxyz,i)= dd0i(1,ixyz)*ddt(jxyz,1)
     &           +dd0i(2,ixyz)*ddt(jxyz,2) +dd0i(3,ixyz)*ddt(jxyz,3)
          enddo
        enddo
c        if( i.eq.1 ) then
c          write(6,'(a,9es11.3)') ' stn=',stn(1:3,1:3,i)
c        endif
c-------J^{T}J
        ddt(1:3,1:3)= stn(1:3,1:3,i)
        do jxyz=1,3
          do ixyz=1,3
            stn(ixyz,jxyz,i)= ddt(ixyz,1)*ddt(jxyz,1)
     &           +ddt(ixyz,2)*ddt(jxyz,2) +ddt(ixyz,3)*ddt(jxyz,3)
          enddo
        enddo
c-------Lagrange strain tensor
        stn(1,1,i)= stn(1,1,i) -1d0
        stn(2,2,i)= stn(2,2,i) -1d0
        stn(3,3,i)= stn(3,3,i) -1d0
        stn(1:3,1:3,i)= 0.5d0 *stn(1:3,1:3,i)
c        if( i.eq.1 ) then
c          write(6,'(a,9es11.3)') ' stn=',stn(1:3,1:3,i)
c        endif
      enddo

      end subroutine get_atomic_strain
c=======================================================================
      subroutine write_force(ionum,cnum,natm,h,tag,fa,epot)
      implicit none
      integer,intent(in):: ionum,natm
      character(len=*),intent(in):: cnum
      real(8),intent(in):: h(3,3),tag(natm),fa(3,natm),epot

      integer:: i

c.....Write forces
      open(ionum,file='frc'//trim(cnum),status='replace')
      write(ionum,'(i10)') natm
      do i=1,natm
        write(ionum,'(3es22.14)') fa(1:3,i)
      enddo
      close(ionum)

c.....Write energy
      open(ionum+1,file='erg'//trim(cnum),status='replace')
      write(ionum+1,'(es22.14)') epot
      close(ionum+1)

      end subroutine write_force
c=======================================================================
      subroutine write_frcav(ionum,cfname)
c
c write 'frcav???' file which includes accumulated and averaged forces
c
      use variables
      implicit none
      integer,intent(in):: ionum
      character(len=*),intent(in):: cfname

      integer:: i

      open(ionum,file=trim(cfname),status='replace')
      write(ionum,'(2i10)') natm,istp+istp0
      do i=1,natm
        write(ionum,'(2es22.14)') facc(i)
      enddo
      close(ionum)

      end subroutine write_frcav
c=======================================================================
      subroutine read_frcav(ionum,cfname)
      use variables
      implicit none
      integer,intent(in):: ionum
      character(len=*),intent(in):: cfname

      integer:: itmp,i
      real(8):: tmp

      facc(1:natm)= 0d0
      open(ionum,file=trim(cfname),status='old',err=100)
      read(ionum,*) itmp,istp0
      if( itmp.ne.natm ) stop ' [read_frcav] itmp.ne.natm !!!'
      do i=1,natm
        read(ionum,*) facc(i)
      enddo
      close(ionum)
      return

 100  write(6,'(a)') ' [read_frcav] '//trim(cfname)//' is not read !'
      return
      end subroutine read_frcav
c=======================================================================
      subroutine read_vnorm(ionum,cfname)
      use variables
      implicit none
      integer,intent(in):: ionum
      character(len=*),intent(in):: cfname

      integer:: itmp,i

      open(ionum,file=trim(cfname),status='old')
      read(ionum,*) itmp
      if( itmp.ne.natm ) stop ' [read_vnorm] itmp.ne.natm !!!'
      do i=1,natm
        read(ionum,*) vn(1:3,i)
      enddo
      close(ionum)
      write(6,'(a)') ' read '//trim(cfname)

      end subroutine read_vnorm
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
      function lconverged(natm,h,tag,fa,feps,istp,fmv,erg)
      use variables, only: ifmv_tag
      implicit none
      integer,intent(in):: natm,istp
      real(8),intent(in):: tag(natm),erg
     &     ,fa(3,natm),feps,h(3,3),fmv(3,0:9)

      logical:: lconverged

      integer:: i,ierr,imax,ifmv
      real(8):: fmax,fi,ft(3),epot,estr
      real(8),external:: get_etime

      lconverged= .false.

      fmax= 0d0
      do i=1,natm
        ifmv= ifmv_tag(tag(i))
        ft(1:3)= fa(1:3,i)*fmv(1:3,ifmv)
        fi= sqrt(ft(1)**2 +ft(2)**2 +ft(3)**2)
        fmax= max(fi,fmax)
      enddo
      write(6,'(a,i6,f10.1,5es22.14)')
     &     ' istp,etime,dE,fmax,feps='
     &     ,istp,get_etime(),erg,fmax,feps
      if( fmax.lt.feps ) then
        write(6,'(a)') ' Convergence achieved for force.'
        write(6,'(a,2es15.7)') '  fmax, feps=',fmax,feps
        lconverged= .true.
      endif

      return
      end function lconverged
c=======================================================================
      function sprod(ndim,v1,v2)
      implicit none
      integer,intent(in):: ndim
      real(8),intent(in):: v1(ndim),v2(ndim)
      real(8):: sprod

      integer:: i

      sprod= 0d0
      do i=1,ndim
        sprod= sprod +v1(i)*v2(i)
      enddo
      return
      end function sprod
c=======================================================================
      subroutine set_vnorm(ndim,vn)
      implicit none
      integer,intent(in):: ndim
      real(8),intent(inout):: vn(ndim)

      real(8):: tmp
      real(8),external:: sprod

      tmp= sprod(ndim,vn,vn)
      vn(1:ndim)= vn(1:ndim)/tmp
      
      end subroutine set_vnorm
c=======================================================================
      subroutine setv()
      use variables
      implicit none
      include 'params_unit.h'
      integer:: i,l,is
      real(8):: dseed,sumvx,sumvy,sumvz,rnd1,rnd2,tmp,facv(nismax)
      real(8),parameter:: pi = 3.14159265358979d0

      facv(1:nismax)=dsqrt(2d0*tinit*fkb/am(1:nismax))

c.....velocities in Maxwell-Boltzmann distribution
      dseed=12345d0
      do i=1,natm
        is= int(tag(i))
        do l=1,3
          call myrnd(rnd1,dseed)
          call myrnd(rnd2,dseed)
          va(l,i)=facv(is)*dsqrt(-dlog(rnd1))*dcos(2d0*pi*rnd2)
        enddo
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

      end subroutine rm_trans_motion
c=======================================================================
      function pbc(x)
      implicit none
      real(8),intent(in):: x
      real(8):: pbc
      
      pbc= x
      if( x.lt.0d0 ) pbc= x +1d0
      if( x.ge.1d0 ) pbc= x -1d0
      return
      end function pbc
c=======================================================================
      subroutine conjugate_gradient()
      use variables
      implicit none
      integer:: i,ifmv
      logical,external:: lconverged
      real(8):: dmax,rt(3),ft(3),di,fnrm,fnrm0,erg0
      real(8),allocatable:: rat(:,:),fa0(:,:)

      open(ioerg,file=cferg,status='replace')

      allocate(rat(3,natm),fa0(3,natm))

      fnrm= 1d0
      fa(1:3,1:natm)= 0d0
      va(1:3,1:natm)= 0d0
      fa0(1:3,1:natm)= fa(1:3,1:natm)
      ekin= 0d0
      temp= 0d0

c.....First calc of forces
      call get_force(natm,h,ra,tag,va,fa,epi,eki,epot0,0)
      epot= epot0

      call output_energy(ioerg,0,ekin,epot,epot0,temp)

c.....Check force convergence
      if( lconverged(natm,h,tag,fa,feps,0,fmv,epot-epot0) )
     &     return

      do istp=1,nstp
        fnrm0= fnrm
        rat(1:3,1:natm)= ra(1:3,1:natm)

c.....Force norm
        fnrm= 0d0
        do i=1,natm
          ft(1:3)= fa(1:3,i)
          fnrm= fnrm +(ft(1)**2 +ft(2)**2 +ft(3)**2)
        enddo

c.....Fletcher-Reeves CG
        fa(1:3,1:natm)= fa(1:3,1:natm)
     &       +fnrm/fnrm0 *fa0(1:3,1:natm)

c.....Multiplying fmv
        do i=1,natm
          ifmv= ifmv_tag(tag(i))
          fa(1:3,i)=fa(1:3,i)*fmv(1:3,ifmv)
        enddo

c.....Line minimization in fa direction to get updated ra
        call line_minimize(natm,h,tag,ra,va,fa,epot)
         
        call output_energy(ioerg,istp,ekin,epot,epot0,temp)

c.....Check force convergence
        if( lconverged(natm,h,tag,fa,feps,istp,fmv,epot-epot0) )
     &         return

        fa0(1:3,1:natm)= fa(1:3,1:natm)
        call get_force(natm,h,ra,tag,va,fa,epi,eki,epot,istp)
      enddo

      write(6,'(a,i8,a)') ' !!! Convergence not achieved within ',nstp
     &     ,' !!!'
      return
      end subroutine conjugate_gradient
c=======================================================================
      subroutine line_minimize(natm,h,tag,ra,va,fa,erg)
      use variables, only : deps,eki,epi
      implicit none
      integer,intent(in):: natm
      real(8),intent(in):: h(3,3),fa(3,natm),tag(natm)
     &     ,va(3,natm)
      real(8),intent(inout):: ra(3,natm),erg
      
      real(8),parameter:: xfac= 2d0
      integer,parameter:: itmax= 100

      integer:: islc,i,it
      real(8):: fi,fmax,ft(3),ds,s0,s1,s2,s3,s4,erg0,erg1,erg2,erg3,erg4
     &     ,tau,d,hmax,epot,estr

      logical,save:: l1st=.true.
      real(8),allocatable,save:: ra1(:,:),ra2(:,:)
     &     ,ra3(:,:),ra4(:,:),fat(:,:)

      if( l1st ) then
        allocate(ra1(3,natm),ra2(3,natm)
     &       ,ra3(3,natm),ra4(3,natm)
     &       ,fat(3,natm))
        ra1(1:3,1:natm)= ra(1:3,1:natm)
        ra2(1:3,1:natm)= ra(1:3,1:natm)
        ra3(1:3,1:natm)= ra(1:3,1:natm)
        ra4(1:3,1:natm)= ra(1:3,1:natm)
        ra1(1:3,1:natm)= ra(1:3,1:natm)
        ra2(1:3,1:natm)= ra(1:3,1:natm)
        ra3(1:3,1:natm)= ra(1:3,1:natm)
        ra4(1:3,1:natm)= ra(1:3,1:natm)
        fat(1:3,1:natm)= fa(1:3,1:natm)
        fat(1:3,1:natm)= fa(1:3,1:natm)
        l1st= .false.
      endif

c.....Get max force
      fmax=0d0
      do i=1,natm
        ft(1:3)=fa(1:3,i)
        fmax= max(fmax,abs(ft(1)))
        fmax= max(fmax,abs(ft(2)))
        fmax= max(fmax,abs(ft(3)))
c          fi= sqrt(ft(1)**2 +ft(2)**2 +ft(3)**2)
c          fmax= max(fi,fmax)
      enddo
c.....First step distance: 0.1 AA
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
        ra2(1:3,1:natm)= ra(1:3,1:natm)
     &       +s2*fa(1:3,1:natm)
        call get_force(natm,h,ra2,tag,va,fat,epi,eki,erg2,it)
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
      ra1(1:3,1:natm)= ra(1:3,1:natm) +s1*fa(1:3,1:natm)
      tau= (1d0 +sqrt(5d0))/2
      do it=1,itmax
        write(6,'(a,i6,2es22.14)') '   LINMIN steps=',it,erg1,erg2
        d= s2-s1
c.....Convergence checking
        if( d*fmax.lt.deps ) exit
        s3= (tau-1d0)/tau *d
        s4= 1d0/tau *d
        ra3(1:3,1:natm)= ra1(1:3,1:natm) +s3*fa(1:3,1:natm)
        ra4(1:3,1:natm)= ra1(1:3,1:natm) +s4*fa(1:3,1:natm)
        call get_force(natm,h,ra3,tag,va,fat,epi,eki,erg3,it)
        call get_force(natm,h,ra4,tag,va,fat,epi,eki,erg4,it)
c.....Update x1 or x4 dependent on relation between erg3 and erg4
        if( erg3.ge.erg4 ) then
          s1= s3
          ra1(1:3,1:natm)= ra3(1:3,1:natm)
          erg1= erg3
        else
          s2= s4
          ra2(1:3,1:natm)= ra4(1:3,1:natm)
          erg2= erg4
        endif
      enddo

      if( erg1.lt.erg2 ) then
        ra(1:3,1:natm)= ra1(1:3,1:natm)
        erg= erg1
      else
        ra(1:3,1:natm)= ra2(1:3,1:natm)
        erg= erg2
      endif

      return
      end subroutine line_minimize
c=======================================================================
      subroutine output_energy(ionum,istp,ekin,epot,epot0,temp)
      implicit none
      integer,intent(in):: istp,ionum
      real(8),intent(in):: ekin,epot,epot0,temp
      real(8),external:: get_etime

      if( istp.eq.0 ) then
        write(ionum,'(a,es15.7)') '# epot0 =',epot0
      endif

      write(ionum,'(i8,3es22.14,f12.3)') istp
     &     ,ekin+epot-epot0,ekin,epot-epot0,temp
      call flush(ionum)

      write(6,'(a,i8,f10.1,2es15.7,f10.3)')
     &     ' istp,etime,ekin,epot,temp='
     &     ,istp,get_etime(),ekin,epot,temp

      end subroutine output_energy
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make qmcl"
c     End:
