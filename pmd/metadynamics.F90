module metadynamics
!-----------------------------------------------------------------------
!                     Last modified: <2019-05-17 13:26:58 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
!  Metadynamics module mainly for creating samples for fitpot.
!-----------------------------------------------------------------------
  implicit none
  save
  include 'mpif.h'

  character(len=128),parameter:: cfparams = 'in.metaD'
  character(len=128),parameter:: cfoutput = 'out.metaD.potential'
  integer,parameter:: ioprms = 40
  integer,parameter:: iopot  = 41

  real(8),parameter:: pi = 3.14159265358979d0

!.....CV type:
!       - species_pair: mainly for the purpose of sample collection for fitpot
!       - bonds: specific bond lengths
!       - bonds_from_atoms: all the bonds connected to given atoms
  character(len=128):: cvtype = 'species_pair'

  real(8):: resrcs
  
  integer:: nhistory
  integer:: nskip
  real(8):: gwidth = 0.2d0
  real(8):: gheight = 0.001d0  ! in eV per Gaussian
  real(8):: rc = 5.0
  real(8):: rc2, dr
  integer:: ndiv = 1024
  real(8),allocatable:: fpair(:,:,:)
!.....Specific to bonds
  integer:: nbond = 0
  real(8),allocatable:: fbonds1(:),fbonds2(:,:)
  integer,allocatable:: ibonds(:,:)
!.....Specific to bonds_from_atoms
  integer:: natm4bnd
  integer,allocatable:: iatm4bnd(:)
  real(8),allocatable:: dbonds(:,:,:)
  
contains
!=======================================================================
  subroutine init_metaD(nstp,nsp,namax,nodes,myid,mpi_world,iprint)
    integer,intent(in):: nstp,nsp,namax,myid,mpi_world,iprint,nodes

    integer:: ia,ib
    character:: cnum*4

    if( nodes.gt.1 ) then
      print *,'ERROR: metadynamics is not available in parallel mode.'
      stop
    endif

    call read_metaD(myid,mpi_world,iprint)
    call sync_params(myid,mpi_world,iprint)

    nskip = nstp/nhistory +1

    resrcs = 0d0

    if( trim(cvtype).eq.'species_pair' ) then
      allocate(fpair(ndiv,nsp,nsp))
      resrcs = resrcs +ndiv*nsp*nsp*8
      fpair(:,:,:) = 0d0
    else if( trim(cvtype).eq.'bonds' ) then
      if( nbond.gt.2 ) then
        print *,'NBOND is limited up to 2, but NBOND=',nbond
        stop
      else if( dble(ndiv)**nbond *8.gt.1d+10 ) then
        print *,'Dimension would be too large, > ' &
             ,dble(ndiv)**nbond*8 /1d+9,' GB'
        stop
      endif
      if( nbond.eq.1 ) then
        resrcs = resrcs +ndiv*8
        allocate(fbonds1(0:ndiv))
        fbonds1(:) = 0d0
      else if( nbond.eq.2 ) then
        resrcs = resrcs +ndiv*ndiv*8
        allocate(fbonds2(0:ndiv,0:ndiv))
        fbonds2(:,:) = 0d0
      endif
    else if( trim(cvtype).eq.'bonds_from_atoms' ) then
      resrcs = namax*nhistory*natm4bnd*8d0
      if( resrcs.gt.1d+10 ) then
        print *,'Dimension would be too large, > ' &
             ,resrcs /1d+9,' GB'
        stop
      endif
      allocate(dbonds(namax,natm4bnd,nhistory))
      dbonds(:,:,:) = 0d0
    endif

    rc2 = rc**2
    dr = rc/ndiv

    if( iprint.gt.0 ) then
      print *,''
      print *,'Metadynamics parameters:'
      print '(a,a)',  '   CV_type = ',trim(cvtype)
      print '(a,i0)', '   num_history = ',nhistory
      print '(a,f8.2)', '   gaussian_width  = ',gwidth
      print '(a,f8.4)', '   gaussian_height = ',gheight
      print '(a,i0)', '   num_division    = ',ndiv
      print '(a,f8.4,f8.5)', '   cutoff_radius,dr= ',rc,dr
      print '(a,f8.4,a)', '   Memory = ',resrcs/1d6,' MB'
      if( trim(cvtype).eq.'bonds' ) then
        print '(a,i0)', '   Num of bonds = ',nbond
        do ib=1,nbond
          print '(a,2(1x,i0))', '     ia,ja = ',ibonds(1,ib),ibonds(2,ib)
        enddo
      else if( trim(cvtype).eq.'bonds_from_atoms' ) then
        print '(a,i0)', '   Num of given atoms = ',natm4bnd
        write(cnum,'(i0)') natm4bnd 
        print '(a,'//trim(cnum)//'(1x,i0))', '    ia = ' &
             ,(iatm4bnd(ia),ia=1,natm4bnd)
      endif
    endif

    return
  end subroutine init_metaD
!=======================================================================
  subroutine read_metaD(myid,mpi_world,iprint)
    use util, only: num_data
    integer,intent(in):: myid,mpi_world,iprint
!!$    integer,external:: num_data

    character(len=128):: c1st, cmode, cline
    integer:: ia,ja,ib,nd
    
    if( myid.eq.0 ) then
      open(ioprms,file=trim(cfparams),status='old')

      do while(.true.)
        read(ioprms,*,end=10) c1st
        if( num_data(c1st,' ').eq.0 ) cycle
        if( c1st(1:1).eq.'!' .or. c1st(1:1).eq.'#' ) cycle
        if( trim(c1st).eq.'CV_type' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, cvtype
          if( trim(cvtype).eq.'bonds' ) then
            cmode = 'bonds'
            ib = 0
          else if( trim(cvtype).eq.'bonds_from_atoms' ) then
            cmode = 'bonds_from_atoms'
          endif
        else if( trim(c1st).eq.'gaussian_width' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, gwidth
        else if( trim(c1st).eq.'gaussian_height' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, gheight
        else if( trim(c1st).eq.'cutoff_radius' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, rc
        else if( trim(c1st).eq.'num_history' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, nhistory
        else if( trim(c1st).eq.'num_division' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, ndiv
        else
          if( trim(cmode).eq.'bonds' ) then
            backspace(ioprms)
            read(ioprms,'(a)') cline
            nd = num_data(cline,' ')
            if( nd.eq.1 ) then
              read(cline,*) nbond
              allocate(ibonds(2,nbond))
            else if( nd.eq.2 ) then
              read(cline,*) ia, ja
              ib = ib + 1
              if( ib.gt.nbond ) then
                if( iprint.gt.0 ) then
                  print *,'WARNING: number of given bonds exceeds NBOND !!'
                endif
              else
                ibonds(1,ib) = ia
                ibonds(2,ib) = ja
              endif
            else
              if( iprint.gt.0 ) then
                print *,'WARNING: wrong format for bonds in in.metaD !!'
              endif
            endif
          else if( trim(cmode).eq.'bonds_from_atoms' ) then
            backspace(ioprms)
            read(ioprms,'(a)') cline
            natm4bnd = num_data(cline,' ')
            allocate(iatm4bnd(natm4bnd))
            read(cline,*) (iatm4bnd(ia),ia=1,natm4bnd)
          endif
        endif
      enddo
10    continue
    endif
    return
  end subroutine read_metaD
!=======================================================================
  subroutine sync_params(myid,mpi_world,iprint)
    integer,intent(in):: myid,mpi_world,iprint
    integer:: ierr

    call mpi_bcast(cvtype,128,mpi_character,0,mpi_world,ierr)
    call mpi_bcast(gwidth,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(gheight,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(rc,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(nhistory,1,mpi_integer,0,mpi_world,ierr)
    if( trim(cvtype).eq.'bonds' ) then
      call mpi_bcast(nbond,1,mpi_integer,0,mpi_world,ierr)
      call mpi_bcast(ibonds,2*nbond,mpi_integer,0,mpi_world,ierr)
    endif
    
    return
  end subroutine sync_params
!=======================================================================
  subroutine update_metaD(istp,namax,natm,nsp,tag,ra,h,nnmax,lspr &
       ,myid,mpi_world,iprint)
!
!  Update the additive potential in metaD.
!
    integer,intent(in):: istp,namax,natm,nnmax,lspr(0:nnmax,namax),nsp
    integer,intent(in):: myid,mpi_world,iprint
    real(8),intent(in):: tag(natm),ra(3,natm),h(3,3)

    integer:: ihist

    if( mod(istp,nskip).ne.0 ) return
    ihist = istp/nskip
    if( trim(cvtype).eq.'species_pair' ) then
      call update_species_pair(namax,natm,nsp,tag,ra,h,nnmax,lspr &
           ,myid,mpi_world,iprint)
    else if( trim(cvtype).eq.'bonds' ) then
      call update_bonds(namax,natm,nsp,tag,ra,h &
           ,myid,mpi_world,iprint,ihist)
    else if( trim(cvtype).eq.'bonds_from_atoms' ) then
      call update_bonds_from_atoms(namax,natm,nsp,tag,ra,h &
           ,nnmax,lspr,myid,mpi_world,iprint,ihist)
    endif

  end subroutine update_metaD
!=======================================================================
  subroutine update_species_pair(namax,natm,nsp,tag,ra,h,nnmax,lspr &
       ,myid,mpi_world,iprint)
    integer,intent(in):: namax,natm,nnmax,lspr(0:nnmax,namax),nsp
    integer,intent(in):: myid,mpi_world,iprint
    real(8),intent(in):: tag(natm),ra(3,natm),h(3,3)

    integer:: ia,ja,is,js,isp,jsp,jj,idiv
    real(8):: xi(3),xj(3),xij(3),rij(3),dij,pref,r,fval,fc

    do ia=1,natm
      is = int(tag(ia))
      xi(1:3)= ra(1:3,ia)
      do jj=1,lspr(0,ia)
        ja = lspr(jj,ia)
        js = int(tag(ja))
        if( is.le.js ) then
          isp = is
          jsp = js
        else
          isp = js
          jsp = is
        endif
        xj(1:3) = ra(1:3,ja)
        xij(1:3) = xj(1:3) -xi(1:3)
        rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij = rij(1)**2 +rij(2)**2 +rij(3)**2
        if( dij.gt.rc2 ) cycle
        dij = sqrt(dij)
        pref = 1d0/(4d0*pi*dij*dij*dr)
        fc = fcut(dij,0d0,rc)
        do idiv=1,ndiv
          r = idiv*dr
          fval = gheight*exp(-(r-dij)**2/2/gwidth**2)*pref
          fpair(idiv,isp,jsp) = fpair(idiv,isp,jsp) +fval*fc
        enddo
      enddo
    enddo
!.....Symmetrize
    do isp=1,nsp-1
      do jsp=isp+1,nsp
        fpair(:,jsp,isp) = fpair(:,isp,jsp)
      enddo
    enddo
    return
  end subroutine update_species_pair
!=======================================================================
  function fcut(r,rin,rout)
!
!  Cutoff function working between rin and rout.
!
    real(8),intent(in):: r,rin,rout
    real(8):: fcut

    if( r.lt.rin ) then
      fcut = 1d0
    else if( rin.le.r .and. r.lt.rout ) then
      fcut = 0.5d0 *(1d0 +cos((r-rin)/(rout-rin)*pi))
    else
      fcut = 0d0
    endif
    return
  end function fcut
!=======================================================================
  subroutine update_bonds(namax,natm,nsp,tag,ra,h &
           ,myid,mpi_world,iprint,ihist)
    use util,only: itotOf
    integer,intent(in):: namax,natm,nsp
    integer,intent(in):: myid,mpi_world,iprint,ihist
    real(8),intent(in):: tag(natm),ra(3,natm),h(3,3)

    integer:: ia,ja,iat,jat,i,itot,ib,idiv1,idiv2
    real(8):: xi(3),xj(3),xij(3),rij(3),dij,rbonds(2),d1,d2,tmp
!!$    integer,external:: itotOf

    do ib=1,nbond
      ia = ibonds(1,ib)
      ja = ibonds(2,ib)
      iat = 0
      jat = 0
      do i=1,natm
        itot = itotOf(tag(i))
        if( itot.eq.ia ) iat = i
        if( itot.eq.ja ) jat = i
        if( iat.gt.0 .and. jat.gt.0 ) exit
      enddo
      xi(1:3) = ra(1:3,iat)
      xj(1:3) = ra(1:3,jat)
      xij(1:3) = xj(1:3)-xi(1:3) -anint(xj(1:3)-xi(1:3))
      rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
      dij = sqrt(rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3))
      rbonds(ib) = dij
!!$      print *,'ib,ia,iat,ja,jat,dij=',ib,ia,iat,ja,jat,dij
!!$      print *,'ia,iat,xi=',ia,iat,xi(1:3)
!!$      print *,'ja,jat,xi=',ja,jat,xj(1:3)
    enddo

    if( nbond.eq.1 ) then
      do idiv1=0,ndiv
        d1 = dr*idiv1
        fbonds1(idiv1) = fbonds1(idiv1) +gheight*exp(-(d1 -rbonds(1))**2 /2/gwidth**2)
      enddo
    else if( nbond.eq.2 ) then
      do idiv2=0,ndiv
        d2 = dr*idiv2
        tmp = -(d2 -rbonds(2))**2 /2 /gwidth
        do idiv1=0,ndiv
          d1 = dr*idiv1
          fbonds2(idiv1,idiv2) = fbonds2(idiv1,idiv2) &
               + gheight *exp( -(d1 -rbonds(1))**2 /2 /gwidth**2&
               +tmp )
        enddo
      enddo
    endif

    return
  end subroutine update_bonds
!=======================================================================
  subroutine update_bonds_from_atoms(namax,natm,nsp,tag,ra,h &
       ,nnmax,lspr,myid,mpi_world,iprint,ihist)
    use util,only: itotOf
    integer,intent(in):: namax,natm,nsp,nnmax,lspr(0:nnmax,namax)
    integer,intent(in):: myid,mpi_world,iprint,ihist
    real(8),intent(in):: tag(natm),ra(3,natm),h(3,3)

    integer:: i,j,ia,ja,jat,iat
    real(8):: xi(3),xj(3),xij(3),rij(3),dij
!!$    integer,external:: itotOf

    do i=1,natm4bnd
      iat = iatm4bnd(i)
      ia = ia_from_itot(iat,natm,tag)
      xi(1:3) = ra(1:3,ia)
      do j=1,lspr(0,ia)
        ja = lspr(j,ia)
        jat = itotOf(tag(ja))
        xj(1:3) = ra(1:3,jat)
        xij(1:3) = xj(1:3)-xi(1:3) -anint(xj(1:3) -xi(1:3))
        rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij = sqrt(rij(1)**2 +rij(2)**2 +rij(3)**2)
        dbonds(jat,ia,ihist) = dij
      enddo
    enddo
  end subroutine update_bonds_from_atoms
!=======================================================================
  function ia_from_itot(itot,natm,tag)
    use util,only: itotOf
    integer,intent(in):: itot,natm
    real(8),intent(in):: tag(natm)
    integer:: ia_from_itot
!!$    integer,external:: itotOf
    integer:: ia

    ia_from_itot = 0
    do ia=1,natm
      if( itot.eq.itotOf(tag(ia)) ) then
        ia_from_itot = ia
        return
      endif
    enddo
    return
  end function ia_from_itot
!=======================================================================
  subroutine force_metaD(istp,namax,natm,tag,ra,aa,h,hi,epot &
       ,nnmax,lspr,myid,mpi_world,iprint)
!
!  Wrapper for calling different force routine for given CV.
!
    integer,intent(in):: namax,natm,istp,myid,mpi_world,iprint,nnmax&
         ,lspr(0:nnmax,namax)
    real(8),intent(in):: tag(namax),ra(3,namax),h(3,3),hi(3,3)
    real(8),intent(inout):: aa(3,namax),epot

    integer:: ihist

    ihist = istp/nskip
    if( trim(cvtype).eq.'species_pair' ) then
      call force_species_pair(namax,natm,tag,ra,aa,h,hi,epot &
           ,nnmax,lspr,myid,mpi_world,iprint)
    else if( trim(cvtype).eq.'bonds' ) then
      call force_bonds(namax,natm,tag,ra,aa,h,hi,epot &
           ,nnmax,lspr,myid,mpi_world,iprint,ihist)
    else if( trim(cvtype).eq.'bonds_from_atoms' ) then
      call force_bonds_from_atoms(namax,natm,tag,ra,aa,h,hi,epot &
           ,nnmax,lspr,myid,mpi_world,iprint,ihist)
    endif
  end subroutine force_metaD
!=======================================================================
  subroutine force_species_pair(namax,natm,tag,ra,aa,h,hi,epot&
       ,nnmax,lspr,myid,mpi_world,iprint)
!
!  CV is pair_species.
!
!  See, J.E. Herr, K. Yao, R. McIntyre, D.W. Toth, and J. Parkhill,
!  J. Chem. Phys. 148, 241710 (2018).
!
    integer,intent(in):: namax,natm,myid,mpi_world,iprint,nnmax &
         ,lspr(0:nnmax,namax)
    real(8),intent(in):: tag(namax),ra(3,namax),h(3,3),hi(3,3)
    real(8),intent(inout):: aa(3,namax),epot

    integer:: ia,ja,is,js,isp,jsp,ierr,jj
    real(8):: xi(3),xj(3),xij(3),rij(3),dij,dxdi(3),dxdj(3),fval,dfval &
         ,epotl,at(3)
    real(8),save,allocatable:: aal(:,:)

    if( .not. allocated(aal) ) then
      allocate(aal(3,namax))
    endif
    if( size(aal).ne.3*namax ) then
      deallocate(aal)
      allocate(aal(3,namax))
    endif
    
    epotl = 0d0
    aal(:,:) = 0d0
    
    do ia=1,natm
      is = int(tag(ia))
      xi(1:3)= ra(1:3,ia)
      do jj=1,lspr(0,ia)
        ja = lspr(jj,ia)
        js = int(tag(ja))
        isp = is
        jsp = js
        if( isp.gt.jsp ) then
          isp = js
          jsp = is
        endif
        xj(1:3) = ra(1:3,ja)
        xij(1:3) = xj(1:3) -xi(1:3)
        rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij = rij(1)**2 +rij(2)**2 +rij(3)**2
        if( dij.gt.rc2 ) cycle
        dij = sqrt(dij)
        fval = fpair_at(dij,isp,jsp)
        if( ja.le.natm ) then
          epotl = epotl + fval
        else
          epotl = epotl +0.5d0
        endif
        dfval = dfpair_at(dij,isp,jsp)
        dxdi(1:3) = -rij(1:3)/dij
        dxdj(1:3) =  rij(1:3)/dij
        aal(1:3,ia) = aal(1:3,ia) -dxdi(1:3)*dfval
        aal(1:3,ja) = aal(1:3,ja) -dxdj(1:3)*dfval
      enddo
    enddo

    do ia=1,natm
      at(1:3)= aal(1:3,ia)
      aal(1:3,ia)= hi(1:3,1)*at(1) +hi(1:3,2)*at(2) +hi(1:3,3)*at(3)
      aa(1:3,ia) = aa(1:3,ia) + aal(1:3,ia)
    enddo
    
!!$!.....Gather epot
!!$    call mpi_allreduce(epotl,epot,1,mpi_real8 &
!!$         ,mpi_sum,mpi_world,ierr)
    
    return
  end subroutine force_species_pair
!=======================================================================
  function fpair_at(r,isp,jsp)
!
!  F value at r
!
    real(8),intent(in):: r
    integer,intent(in):: isp,jsp
    real(8):: fpair_at

    integer:: ir

    ir = int(r/dr)+1
    fpair_at = fpair(ir,isp,jsp)
    return
  end function fpair_at
!=======================================================================
  function dfpair_at(r,isp,jsp)
!
!  dF value at r
!
    real(8),intent(in):: r
    integer,intent(in):: isp,jsp
    real(8):: dfpair_at

    integer:: ir

    ir = int(r/dr)+1
    if( ir.eq.1 ) then
      dfpair_at = (fpair(ir+1,isp,jsp)-fpair(ir,isp,jsp))/dr
    else if( ir.eq.ndiv ) then
      dfpair_at = (fpair(ir,isp,jsp)-fpair(ir-1,isp,jsp))/dr
    else
      dfpair_at = (fpair(ir+1,isp,jsp)-fpair(ir-1,isp,jsp))/(2*dr)
    endif
    return
  end function dfpair_at
!=======================================================================
  subroutine force_bonds(namax,natm,tag,ra,aa,h,hi,epot&
       ,nnmax,lspr,myid,mpi_world,iprint,ihist)
!
!  CV is bonds.
!
    use util,only: itotOf
    integer,intent(in):: namax,natm,myid,mpi_world,iprint,nnmax &
         ,lspr(0:nnmax,namax),ihist
    real(8),intent(in):: tag(namax),ra(3,namax),h(3,3),hi(3,3)
    real(8),intent(inout):: aa(3,namax),epot

    integer:: ia,ja,is,js,isp,jsp,ierr,jj,i,iat,jat,itot,ib
    real(8):: xi(3),xj(3),xij(3),rij(3),dij,dxdi(3),dxdj(3),fval,dfval &
         ,epotl,at(3),rbonds(2)
    real(8),save,allocatable:: aal(:,:),dbdr(:)
!!$    integer,external:: itotOf

    if( .not.allocated(dbdr) ) allocate(dbdr(nbond))

    if( .not. allocated(aal) ) then
      allocate(aal(3,namax))
    endif
    if( size(aal).ne.3*namax ) then
      deallocate(aal)
      allocate(aal(3,namax))
    endif

    epotl = 0d0
    aal(:,:) = 0d0

    if( nbond.eq.1 ) then
      ib = 1
      ia = ibonds(1,ib)
      ja = ibonds(2,ib)
      iat = 0
      jat = 0
      do i=1,natm
        itot = itotOf(tag(i))
        if( itot.eq.ia ) iat = i
        if( itot.eq.ja ) jat = i
        if( iat.gt.0 .and. jat.gt.0 ) exit
      enddo
      xi(1:3) = ra(1:3,iat)
      xj(1:3) = ra(1:3,jat)
      xij(1:3) = xj(1:3)-xi(1:3) -anint(xj(1:3)-xi(1:3))
      rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
      dij = sqrt(rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3))
      dxdi(1:3) = -rij(1:3)/dij
      dxdj(1:3) =  rij(1:3)/dij
      epotl = epotl +fbonds1_at(dij)
      dbdr(1) = dfbonds1_at(dij)
      aal(1:3,iat) = -dbdr(1)*dxdi(1:3)
      aal(1:3,jat) = -dbdr(1)*dxdj(1:3)
    else
      do ib=1,nbond
        ia = ibonds(1,ib)
        ja = ibonds(2,ib)
        iat = 0
        jat = 0
        do i=1,natm
          itot = itotOf(tag(i))
          if( itot.eq.ia ) iat = i
          if( itot.eq.ja ) jat = i
          if( iat.gt.0 .and. jat.gt.0 ) exit
        enddo
        xi(1:3) = ra(1:3,iat)
        xj(1:3) = ra(1:3,jat)
        xij(1:3) = xj(1:3)-xi(1:3) -anint(xj(1:3)-xi(1:3))
        rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij = sqrt(rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3))
        rbonds(ib) = dij
      enddo
      epotl = epotl +fbonds2_at(rbonds(1),rbonds(2))
      call dfbonds2_at(rbonds(1),rbonds(2),dbdr(1),dbdr(2))
      do ib=1,nbond
        ia = ibonds(1,ib)
        ja = ibonds(2,ib)
        iat = 0
        jat = 0
        do i=1,natm
          itot = itotOf(tag(i))
          if( itot.eq.ia ) iat = i
          if( itot.eq.ja ) jat = i
          if( iat.gt.0 .and. jat.gt.0 ) exit
        enddo
        xi(1:3) = ra(1:3,iat)
        xj(1:3) = ra(1:3,jat)
        xij(1:3) = xj(1:3)-xi(1:3) -anint(xj(1:3)-xi(1:3))
        rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij = sqrt(rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3))
        dxdi(1:3) = -rij(1:3)/dij
        dxdj(1:3) =  rij(1:3)/dij
        aal(1:3,iat) = aal(1:3,iat) -dbdr(ib)*dxdi(1:3)
        aal(1:3,jat) = aal(1:3,jat) -dbdr(ib)*dxdj(1:3)
      enddo
    endif

    do ia=1,natm
      at(1:3)= aal(1:3,ia)
      aal(1:3,ia)= hi(1:3,1)*at(1) +hi(1:3,2)*at(2) +hi(1:3,3)*at(3)
      aa(1:3,ia) = aa(1:3,ia) + aal(1:3,ia)
    enddo

    return
  end subroutine force_bonds
!=======================================================================
  function fbonds1_at(r1)
    real(8),intent(in):: r1
    real(8):: fbonds1_at

    integer:: ir1

    ir1 = int(r1/dr) +1
    fbonds1_at = fbonds1(ir1)
    return
  end function fbonds1_at
!=======================================================================
  function fbonds2_at(r1,r2)
    real(8),intent(in):: r1,r2
    real(8):: fbonds2_at

    integer:: ir1,ir2

    ir1 = int(r1/dr) +1
    ir2 = int(r2/dr) +1
    fbonds2_at = fbonds2(ir1,ir2)
    return
  end function fbonds2_at
!=======================================================================
  function dfbonds1_at(r1)
    real(8),intent(in):: r1
    real(8):: dfbonds1_at

    integer:: ir1

    ir1 = int(r1/dr) +1

    if( ir1.eq.1 ) then
      dfbonds1_at = (fbonds1(ir1+1) -fbonds1(ir1))/dr
    else if( ir1.eq.ndiv ) then
      dfbonds1_at = (fbonds1(ir1) -fbonds1(ir1-1))/dr
    else
      dfbonds1_at = (fbonds1(ir1+1) -fbonds1(ir1-1))/(2*dr)
    endif
    return
  end function dfbonds1_at
!=======================================================================
  subroutine dfbonds2_at(r1,r2,df1,df2)
    real(8),intent(in):: r1,r2
    real(8),intent(out):: df1,df2

    integer:: ir1,ir2,ir1p,ir1m,ir2p,ir2m
    real(8):: dr1,dr2

    ir1 = int(r1/dr) +1
    ir2 = int(r2/dr) +1

    ir1p = ir1+1
    ir1m = ir1-1
    dr1 = dr *2
    if( ir1.eq.1 ) then
      ir1m = ir1
      dr1 = dr
    else if( ir1.eq.ndiv ) then
      ir1p = ir1
      dr1 = dr
    endif
    ir2p = ir2+1
    ir2m = ir2-1
    dr2 = dr *2
    if( ir2.eq.1 ) then
      ir2m = ir2
      dr2 = dr
    else if( ir2.eq.ndiv ) then
      ir2p = ir2
      dr2 = dr
    endif
    df1 = (fbonds2(ir1p,ir2) -fbonds2(ir1m,ir2))/dr1
    df2 = (fbonds2(ir1,ir2p) -fbonds2(ir1,ir2m))/dr2
    return
  end subroutine dfbonds2_at
!=======================================================================
  subroutine force_bonds_from_atoms(namax,natm,tag,ra,aa,h,hi,epot&
       ,nnmax,lspr,myid,mpi_world,iprint,ihist)
!
!  CV is bonds_from_atoms
!
    use util,only: itotOf
    integer,intent(in):: namax,natm,myid,mpi_world,iprint,nnmax &
         ,lspr(0:nnmax,namax),ihist
    real(8),intent(in):: tag(namax),ra(3,namax),h(3,3),hi(3,3)
    real(8),intent(inout):: aa(3,namax),epot

    integer:: ia,ja,is,js,isp,jsp,ierr,jj,i,iat,jat,itot,ib,j,jt,ih
    real(8):: xi(3),xj(3),xij(3),rij(3),dij,dxdi(3),dxdj(3) &
         ,epotl,at(3),texp,tmp
    logical:: l_jat_in_lspr 
    real(8),save,allocatable:: aal(:,:)
!!$    integer,external:: itotOf

    if( .not. allocated(aal) ) then
      allocate(aal(3,namax))
    endif
    if( size(aal).ne.3*namax ) then
      deallocate(aal)
      allocate(aal(3,namax))
    endif

    epotl = 0d0
    aal(:,:) = 0d0

    do ih=1,ihist
      do i=1,natm4bnd
        iat = iatm4bnd(i)
        ia = ia_from_itot(iat,natm,tag)
        xi(1:3) = ra(1:3,ia)
        tmp = 0d0
        do ja=1,natm
          jat = itotOf(tag(ja))
          l_jat_in_lspr = .false.
          do j=1,lspr(0,ia)
            jj = lspr(j,ia)
            jt = itotOf(tag(jj))
            if( jat.eq.jt ) then
              l_jat_in_lspr = .true.
              exit
            endif
          enddo
          if( .not. l_jat_in_lspr ) cycle
          xj(1:3) = ra(1:3,jat)
          xij(1:3) = xj(1:3)-xi(1:3) -anint(xj(1:3)-xi(1:3))
          rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
          dij = sqrt(rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3))
          tmp = tmp -(dij -dbonds(jat,ia,ih))**2 /2/ gwidth**2
        enddo
        texp = exp(tmp)
        do ja=1,natm
          jat = itotOf(tag(ja))
          l_jat_in_lspr = .false.
          do j=1,lspr(0,ia)
            jj = lspr(j,ia)
            jt = itotOf(tag(jj))
            if( jat.eq.jt ) then
              l_jat_in_lspr = .true.
              exit
            endif
          enddo
          if( .not. l_jat_in_lspr ) cycle
          xj(1:3) = ra(1:3,jat)
          xij(1:3) = xj(1:3)-xi(1:3) -anint(xj(1:3)-xi(1:3))
          rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
          dij = sqrt(rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3))
          dxdi(1:3) = -rij(1:3)/dij
          dxdj(1:3) =  rij(1:3)/dij
          tmp = -(dij -dbonds(jat,ia,ih))/gwidth**2 *texp
          aal(1:3,ia) = aal(1:3,ia) -dxdi(1:3)*tmp
          aal(1:3,ja) = aal(1:3,ja) -dxdj(1:3)*tmp
        enddo
      enddo
    enddo

    do ia=1,natm
      at(1:3)= aal(1:3,ia)
      aal(1:3,ia)= hi(1:3,1)*at(1) +hi(1:3,2)*at(2) +hi(1:3,3)*at(3)
      aa(1:3,ia) = aa(1:3,ia) + aal(1:3,ia)
    enddo

    return
  end subroutine force_bonds_from_atoms
!=======================================================================
  subroutine write_metaD_potential(istp,nsp,myid,iprint)
!
!  Write metaD potential to the file
!
    integer,intent(in):: nsp,myid,iprint,istp

    integer:: idiv,isp,jsp,idiv1,idiv2
    real(8):: r1,r2

    if( myid.eq.0 ) then
      open(iopot,file=trim(cfoutput),status='replace')
      if( trim(cvtype).eq.'species_pair' ) then
        write(iopot,'(a)') '# rad, fpair(:,:,:)'
        do idiv=1,ndiv
          write(iopot,'(f8.4,100es13.4e3)') idiv*dr,((fpair(idiv,isp,jsp) &
               ,jsp=isp,nsp),isp=1,nsp)
        enddo
      else if( trim(cvtype).eq.'bonds' ) then
        write(iopot,'(a)') '# r1, r2, fbonds(:,:)'
        if( nbond.eq.1 )  then
          do idiv1=1,ndiv
            r1 = idiv1*dr
            write(iopot,'(f8.4, es13.4e3)') r1,fbonds1(idiv1)
          enddo
        else if( nbond.eq.2 ) then
          do idiv1=1,ndiv
            r1 = idiv1*dr
            do idiv2=1,ndiv
              r2 = idiv2*dr
              write(iopot,'(f8.4, f8.4, es13.4e3)') r1,r2,fbonds2(idiv1,idiv2)
            enddo
            write(iopot,*) ''
          enddo
        endif
      endif
      close(iopot)
    endif
  end subroutine write_metaD_potential
end module metadynamics
