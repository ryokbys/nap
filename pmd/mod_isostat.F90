module isostat
!-----------------------------------------------------------------------
!                     Last modified: <2024-07-25 22:30:41 KOBAYASHI Ryo>
!-----------------------------------------------------------------------
! Isothermal and/or isobaric ensemble.
! Note that some variables used in this module are defined in pmdvars not here.
!
! Isothermal algorithms:
!   - Berendsen
!   - Langevin
!     - Mannella algorithm
!     - GJF (Gronbech-Jensen & Farago)
! Isobaric algorithms:
!   - Berendsen
!   - Langevin
!     - GJF
!-----------------------------------------------------------------------
! References:
!   1. Gronbech-Jensen and Farago, J Chem Phys 141, 194108 (2014)
!   2. PhD thesis by David Quigley (2005)
!-----------------------------------------------------------------------
  use pmdvars,only: tinit,tfin,trlx,nstp,istp,tgmm,ttgt,dt,am,fa2v,tfac, &
       ndof,eks,temps,nfmv,cmass,cgmm,nspmax,srlx,ttgt_lang,stgt,ptgt, &
       pini,pfin,cpctl,ifdmp,strs,eki,sgm,dt,stbeta,vol,natm,lhydrostatic, &
       ntemps, lmultemps
  use util,only: ithOf
  use random,only: box_muller
  implicit none
  include "./params_unit.h"
  include "./const.h"
  private
  save

!.....Max strain rate (1%)
  real(8):: sratemax = 0.01d0
  
  public:: setup_langevin, setup_cell_langevin, setup_cell_berendsen, &
       setup_cell_min, &
       vel_update_langevin, vel_update_berendsen, &
       cell_update_berendsen, cell_force_berendsen, &
       cell_update_langevin, cvel_update_langevin, sratemax
       
contains
!=======================================================================
  subroutine setup_langevin(myid,iprint)
    implicit none 
    integer,intent(in):: myid,iprint

    integer:: itemp
    
    tgmm= 1d0/trlx
    do itemp=1,ntemps
      if( ttgt(itemp).lt.0d0 ) then
        tfac(itemp)= -1d0
      else
!.....TFAC should have [sqrt(ue/fs**2)] = [ue/Ang/sqrt(ump)] unit
!            tfac(itemp)= sqrt(2d0*tgmm*(fkb*ev2j)*ttgt(itemp)/dt)
!     &           *m2ang/s2fs
!            tfac(itemp)= dsqrt(2d0*tgmm*fkb*ttgt(itemp)/dt)
        tfac(itemp)= dsqrt(2d0*tgmm*ttgt(itemp)/dt *k2ue)
      endif
    enddo
    if( myid.eq.0 .and. iprint.ge.ipl_basic ) then
      print *,''
      print *,'Langevin thermostat parameters:'
      print '(a,f10.2)','   Relaxation time = ',trlx
      do itemp=1,nfmv
        print '(a,i3,f10.2,es15.4)','   itemp,ttgt,tfac = ' &
             ,itemp,ttgt(itemp),tfac(itemp)
      enddo
    endif

  end subroutine setup_langevin
!=======================================================================
  subroutine setup_cell_langevin(myid,iprint)
!
!  Set the parameters for Langevin barostat.
!
!!$    use pmdvars,only: cmass,cgmm,nspmax,srlx,ttgt_lang,stgt,ptgt, &
!!$         pini,pfin,cpctl,ifdmp
    implicit none
    integer,intent(in):: myid,iprint
    real(8),parameter:: temp_min = 10d0
    integer:: is,ndoftot
    real(8):: ttgtmax,temp,j2
    
    ndoftot = 0
    ttgtmax = 0d0
    do is=1,nspmax
      ndoftot = ndoftot +ndof(is)
      ttgtmax = max(ttgt(is),ttgtmax)
    enddo
    ttgt_lang = max(temp_min, ttgtmax)
    cgmm = 1d0/srlx
!.....T [K] * k2ev ==> [eV], so [cmass] == [eV*fs**2]
    cmass = ndoftot *ttgt_lang /3 /cgmm**2 *k2ev
!.....In case of damping, set T as 0, but mass as finite value.
    if( ifdmp.gt.0 ) then
      cmass = ndoftot *temp_min /3 /cgmm**2 *k2ev
      ttgt_lang = 0d0
    endif


    if( index(cpctl,'vv-').ne.0 .and. abs(pini-pfin).gt.0.1d0 ) ptgt = pini
    if( myid.eq.0 .and. iprint.ge.ipl_basic ) then
      print *,''
      print *,'Langevin barostat parameters:'
      print '(a,f11.2,a)', '   Relaxation time = ',srlx,' fs'
      print '(a,es11.3,a)','   Cell mass       = ',cmass,' ump*ang**2'
      print '(a,f11.2,a)', '   Temperature     = ',ttgt_lang,' K'
      if( index(cpctl,'vv-').ne.0 ) then
        write(6,'(a)') '   Variable-volume to'
        if( abs(pini-pfin).gt.0.1d0 ) then
          write(6,'(a,f0.3,a,f0.3,a)') &
               '   target pressure from = ',pini,' GPa to ',pfin,' GPa'
        else
          write(6,'(a,f11.2,a)') &
               '   target pressure = ',ptgt,' GPa'          
        endif
      else
        write(6,'(a)') '   Variable-cell to'
        write(6,'(a,6f10.2)') '   target stress [GPa]: ' &
             ,stgt(1,1),stgt(2,2),stgt(3,3) &
             ,stgt(2,3),stgt(3,1),stgt(1,2)
      endif
      print '(a,f11.2,a)', '   Max strain rate = ',sratemax
    endif
    ptgt = ptgt *gpa2up
    stgt(:,:) = stgt(:,:)*gpa2up
    return
  end subroutine setup_cell_langevin
!=======================================================================
  subroutine setup_cell_berendsen(myid,iprint)
!!$    use pmdvars,only: cpctl,ptgt,stgt,pini,pfin
    implicit none 
    integer,intent(in):: myid,iprint

    if( index(cpctl,'vv-').ne.0 ) then
      if( abs(pini-pfin).gt. 0.1d0 ) then
        if(myid.eq.0 .and. iprint.ne.0 ) then
          write(6,*) ''
          write(6,'(a)') ' Barostat: variable-volume Berendsen'
          write(6,'(a,f0.3,a,f0.3,a)') &
               '   Target pressure from = ',pini,' GPa to ' &
               ,pfin,' GPa'
        endif
        ptgt = pini *gpa2up
      else
        if(myid.eq.0 .and. iprint.ne.0 ) then
          write(6,*) ''
          write(6,'(a)') ' Barostat: variable-volume Berendsen'
          write(6,'(a,f0.3,a)') &
               '   Target pressure = ',ptgt,' GPa'
        endif
        ptgt = ptgt *gpa2up
      endif
    else
      if(myid.eq.0 .and. iprint.ne.0 ) then
        write(6,*) ''
        write(6,'(a)') ' Barostat: variable-cell Berendsen'
        if( lhydrostatic ) then
          write(6,'(a,l)') '   Hydrostatic pressure: TRUE'
          write(6,'(a,3f10.3)') '   Target stress [GPa]: ' &
               ,stgt(1,1),stgt(2,2),stgt(3,3)
        else
          write(6,'(a,6f10.3)') '   Target stress [GPa]: ' &
               ,stgt(1,1),stgt(2,2),stgt(3,3) &
               ,stgt(2,3),stgt(3,1),stgt(1,2)
        endif
      endif
      stgt(1:3,1:3)= stgt(1:3,1:3) *gpa2up
    endif
    return
  end subroutine setup_cell_berendsen
!=======================================================================
  subroutine setup_cell_min(myid,iprint)
    implicit none
    integer,intent(in):: myid,iprint

    if(myid.eq.0 .and. iprint.ne.0 ) then
      write(6,*) ''
      if( index(cpctl,'vc').ne.0 ) then
        write(6,'(a)') ' Variable-cell'
      else if( index(cpctl,'vv').ne.0 ) then
        write(6,'(a)') ' Variable-volume'
      endif
      write(6,'(a,6f10.3)') '   Target stress [GPa]: ' &
           ,stgt(1,1),stgt(2,2),stgt(3,3) &
           ,stgt(2,3),stgt(3,1),stgt(1,2)
    endif
    stgt(1:3,1:3)= stgt(1:3,1:3) *gpa2up

  end subroutine setup_cell_min
!=======================================================================
  subroutine vel_update_langevin(natm,tag,va,aa)
!
!  Update of velocities in Langevin thermostat using G-JF algorithm.
!
    implicit none
    integer,intent(in):: natm
    real(8),intent(in):: tag(natm),aa(3,natm)
    real(8),intent(inout):: va(3,natm)

    integer:: i,is,itemp,l
    real(8):: ami,tmp,eta,bfac,afac,aai(3)

!.....For G-JF algorithm
    eta = tgmm*dt/2
    bfac = 1d0 /(1d0 +eta)
    afac = bfac *(1d0 -eta)

!.....If final temperature is assigned,
!.....target temperatures are forced to be set intermediate temperatures
    if( tfin.gt.0d0 ) then
      ttgt(1:9) = tinit +(tfin-tinit)*istp/nstp
      tfac(1:9) = dsqrt(2d0*tgmm*ttgt(1:9)/dt *k2ue)
    endif

    do i=1,natm
      itemp = ithOf(tag(i),1)  ! Group-ID for itemp == 1
      is = int(tag(i))
      if( itemp.eq.0 ) cycle
      if( tfac(itemp).lt.0d0 ) cycle
      ami= am(is)
!.....Here unit of TMP should be [eV/Ang],
!     whereas TFAC is [eu/Ang/sqrt(ump)], so need to multiply ue2ev
      tmp= tfac(itemp)*dsqrt(ami)*ue2ev
!!$!.....Here the unit of va*tgmm*ami is [ump*(Ang)/fs**2 = ue/(Ang)],
!!$!.....where (Ang) is actually scaled to unitless,
!!$!.....but this should be [eV/Ang], so need to multiply ue2ev.
!!$      do l=1,3
!!$        aai(l)= -va(l,i)*tgmm*ami*ue2ev &
!!$             +tmp*box_muller() !/hscl(l)
!!$      enddo
!!$!.....To compensate the factor 1/2 in fa2v, multiply 2 here.
!!$      va(1:3,i)= va(1:3,i) +aai(1:3)*fa2v(is)*dt *2d0
!.....G-JF algorithm
      do l=1,3
        aai(l) = tmp*box_muller()
      enddo
      va(1:3,i) = afac*va(1:3,i) +fa2v(is)*dt*(aa(1:3,i) +bfac*aai(1:3)*2d0 )
    enddo
  end subroutine vel_update_langevin
!=======================================================================
  subroutine vel_update_berendsen(natm,tag,va)
!
!  2nd update of velocities in Berendsen thermostat.
!
    implicit none 
    integer,intent(in):: natm
    real(8),intent(in):: tag(natm)
    real(8),intent(inout):: va(3,natm)

    integer:: i,is,itemp
    
    tfac(1:9)= 0d0
!.....if final temperature is assigned,
!.....target temperatures are forced to be set intermediate temperatures
    if( tfin.gt.0d0 ) then
      ttgt(1:9) = tinit +(tfin-tinit)*istp/nstp
    endif
    do itemp=1,ntemps
      if(ndof(itemp).le.0 .or. ttgt(itemp).lt.0d0 ) cycle
      temps(itemp)= eks(itemp) *2d0 /fkb /max(ndof(itemp)-3,3)
      if( abs(ttgt(itemp)-temps(itemp))/temps(itemp).gt.100d0 ) then
        tfac(itemp)= dsqrt(1d0 +dt/trlx*100d0 )
      else
        tfac(itemp)= dsqrt(1d0 +dt/trlx*(ttgt(itemp)-temps(itemp)) &
             /temps(itemp))
      endif
    enddo
    do i=1,natm
      itemp = ithOf(tag(i),1)  ! Group-ID for itemp == 1
      if( itemp.le.0 .or. ttgt(itemp).lt.0d0 ) cycle
      va(1:3,i)= va(1:3,i) *tfac(itemp)
    enddo
  end subroutine vel_update_berendsen
!=======================================================================
  subroutine cell_force_berendsen(stnsr,ah,mpi_md_world)
!
!  Calc. deviation of h-matrix used in isobaric MD
!  using given stress tensor stgt.
!
    implicit none
    include "mpif.h"
    integer,intent(in):: mpi_md_world
    real(8),intent(in):: stnsr(3,3)
    real(8),intent(out):: ah(3,3)


    integer:: i,ixyz,jxyz,ierr,l,jm,jp,im,ip
    real(8):: prss,fac,tmp,bxc(3),cxa(3),axb(3),sgmnrm
    real(8):: detah,ahcof(3,3)

!.....now ah is scaling factor for h-mat
      ah(1:3,1:3)= 0d0
      ah(1,1)= 1d0
      ah(2,2)= 1d0
      ah(3,3)= 1d0
!.....Berendsen for variable-cell
    if( trim(cpctl).eq.'berendsen' .or. &
         trim(cpctl).eq.'vc-berendsen' .or. &
         trim(cpctl).eq.'cv-berendsen' ) then
!.....Limit change rate of h (ah) to sratemax
      do jxyz=1,3
        sgmnrm = sqrt(sgm(1,jxyz)**2 +sgm(2,jxyz)**2 +sgm(3,jxyz)**2)
        do ixyz=1,3
          tmp = 0d0
          if( lhydrostatic ) then
            tmp = tmp + ( stgt(ixyz,ixyz)-stnsr(ixyz,ixyz) ) *sgm(ixyz,jxyz)
          else
            do l=1,3
              tmp = tmp + ( stgt(ixyz,l)-stnsr(ixyz,l) ) *sgm(l,jxyz)
            enddo
          endif
          tmp = tmp *(stbeta/gpa2up)*dt/3/srlx /sgmnrm
          tmp = min(max(tmp,-sratemax),sratemax)
          ah(ixyz,jxyz) = ah(ixyz,jxyz) -tmp
        enddo
      enddo
!.....Conserving (constant)-volume
      if( cpctl(1:3).eq.'cv-' ) then
        do jxyz=1,3
          jm=mod(jxyz+1,3)+1
          jp=mod(jxyz,  3)+1
          do ixyz=1,3
            im=mod(ixyz+1,3)+1
            ip=mod(ixyz,  3)+1
            ahcof(ixyz,jxyz)= ah(ip,jp)*ah(im,jm) -ah(im,jp)*ah(ip,jm)
          enddo
        enddo
        detah = ah(1,1)*ahcof(1,1) +ah(2,1)*ahcof(2,1) +ah(3,1)*ahcof(3,1)
        ah(:,:) = ah(:,:) /detah**(1d0/3)
      endif
!.....Berendsen for variable-volume not variable-cell
    else if( trim(cpctl).eq.'vv-berendsen' ) then
      prss = (stnsr(1,1)+stnsr(2,2)+stnsr(3,3))/3
      fac = 1.0 -(stbeta/gpa2up)*dt/3/srlx *(ptgt-prss)
!.....Limit change rate of h (ah) to RMAX
      fac = min(max(fac,1d0-sratemax),1d0+sratemax)
      ah(1:3,1:3) = ah(1:3,1:3)*fac
    endif
    return
  end subroutine cell_force_berendsen
!=======================================================================
  subroutine cell_update_berendsen(ah,h,lcellfix,lcell_updated)
!
!  Update cell vector, h, by Berendsen barostat
!
    implicit none 
    real(8),intent(in):: ah(3,3)
    logical,intent(in):: lcellfix(3,3)
    real(8),intent(inout):: h(3,3,0:1)
    logical,intent(inout):: lcell_updated

    integer:: i,j
    real(8):: htmp(3,3)
    
    htmp(1:3,1:3) = matmul(ah,h(1:3,1:3,0))
    do i=1,3
      do j=1,3
!.....NOTE: The definition of lcellfix from input and
!.....that of h-matrix actually used are in relation of transpose.
        if( .not. lcellfix(j,i) ) then
          h(i,j,0) = htmp(i,j)
        endif
      enddo
    enddo
    lcell_updated = .true.

  end subroutine cell_update_berendsen
!=======================================================================
  subroutine cvel_update_langevin(stnsr,h,mpi_md_world,ikick)
!
!  Update cell velocity in Langevin scheme
!  Parameters:
!    - ikick: 1) 1st-kick, 2) 2nd-kick
!
!!$    use pmdvars,only: cpctl,sgm,dt,srlx,vol,stgt,ptgt,cgmm,cmass,&
!!$         ttgt_lang
    implicit none
    include "mpif.h"
    integer,intent(in):: mpi_md_world,ikick
    real(8),intent(in):: stnsr(3,3)
    real(8),intent(inout):: h(3,3,0:1)

    integer:: ixyz,jxyz,l
    real(8):: tfac_lang,tmp,tmp2,eta,bfac,afac,prss, &
         sdiff(3,3),sgmnrm,dh

!.....For G-JF algorithm
    eta = cgmm *dt/2
    bfac = 1d0 /(1d0 +eta)
    afac = bfac *(1d0 -eta)
    if( ikick.eq.1 ) afac = 1d0
!.....[K] *k2ue ==> special energy unit (ue)
!.....cmass*cgmm/dt = [ump*ang**2 /fs**2] = [ue] *ue2ev = [eV]
    tfac_lang = dsqrt(2d0*cgmm*cmass*(ttgt_lang*k2ue)/dt) *ue2ev
    
    if( index(cpctl,'vv-').ne.0 ) then  ! variable-volume
      prss = (stnsr(1,1) +stnsr(2,2) +stnsr(3,3))/3
      dh = (vol*(prss -ptgt) +tfac_lang*box_muller())*dt/2/cmass
      do ixyz=1,3
        h(ixyz,ixyz,1) = afac*h(ixyz,ixyz,1) +dh
      enddo
    else  ! variable-cell
      sdiff(:,:) = vol*(stnsr(:,:) -stgt(:,:))
      do jxyz=1,3
        sgmnrm = sqrt(sgm(1,jxyz)**2 +sgm(2,jxyz)**2 +sgm(3,jxyz)**2)
        do ixyz=1,3
          tmp2 = 0d0
          do l=1,3
            tmp2 = tmp2 + sdiff(ixyz,l)*sgm(l,jxyz)
          enddo
          tmp2 = tmp2/sgmnrm
          h(ixyz,jxyz,1) = afac*h(ixyz,jxyz,1) &
               +(tmp2+tfac_lang*box_muller())*dt/2/cmass
        enddo
      enddo
    endif
    return
  end subroutine cvel_update_langevin
!=======================================================================
  subroutine cell_update_langevin(h,lcellfix,lcell_updated)
!
!  Update cell vectors, h, by Langevin barostat
!
!!$    use pmdvars,only: dt,cgmm,cmass
    implicit none
    logical,intent(in):: lcellfix(3,3)
    real(8),intent(inout):: h(3,3,0:1)
    logical,intent(inout):: lcell_updated

    integer:: ixyz,jxyz
    real(8):: eta,bfac,ah(3,3),htmp(3,3),elemax,fac
    
!.....For G-JF algorithm
    eta = cgmm *dt/2
    bfac = 1d0 /(1d0 +eta)
    ah(1:3,1:3) = 0d0
    elemax = 0d0
    do jxyz=1,3
      do ixyz=1,3
        ah(ixyz,jxyz) = bfac*dt*h(ixyz,jxyz,1) !/cmass
        elemax = max( elemax, abs(ah(ixyz,jxyz)) )
      enddo
    enddo
    if( elemax.gt.sratemax ) then
      fac = sratemax/elemax
      ah(:,:) = ah(:,:) *fac
      h(:,:,1) = h(:,:,1) *fac
    endif
    ah(1,1) = ah(1,1) +1d0
    ah(2,2) = ah(2,2) +1d0
    ah(3,3) = ah(3,3) +1d0
    
    htmp(1:3,1:3) = matmul(ah,h(1:3,1:3,0))
    do ixyz=1,3
      do jxyz=1,3
!.....NOTE: The definition of lcellfix from input and
!.....that of h-matrix actually used are in relation of transpose.
        if( .not. lcellfix(jxyz,ixyz) ) then
          h(ixyz,jxyz,0) = htmp(ixyz,jxyz)
        endif
      enddo
    enddo
    lcell_updated = .true.
    return
  end subroutine cell_update_langevin
end module isostat
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd lib"
!     End:

