module constraints
!
!  Currently only one bond distant constraint is available.
!
  use vector,only: abc2cart, cart2abc, dot, norm2
  use util, only: num_data, itotOf, is_numeric
  implicit none
  save
  include 'mpif.h'

  character(len=128):: paramsdir = '.'
  character(len=128),parameter:: cfparams = 'in.constraints'
  integer,parameter:: ioprms = 66

!.....Constraint type: [bonds]
  character(len=128):: ctype_const = 'bonds'
  integer:: nconst = 0
  integer:: maxiter = 100
  real(8),allocatable:: dij(:),dij2(:),dini(:),dfin(:),dtol2(:)
  integer,allocatable:: idcs0(:,:),idcs(:,:),idxyz(:)
  real(8),allocatable:: rio(:,:),rjo(:,:)
  real(8),allocatable:: ris(:,:),rjs(:,:),vis(:,:),vjs(:,:)
  real(8),allocatable:: rijo(:,:)
  real(8),allocatable:: dxyzi(:), dxyzf(:),dxyz(:)
  character(len=20),allocatable:: constype(:)
  real(8):: tol = 1d-2
  real(8):: vtol = 1d-8
  
  
contains
!=======================================================================
  subroutine init_const(myid,mpi_world,nnode,iprint,h)
    integer,intent(in):: myid,mpi_world,nnode,iprint
    real(8),intent(in):: h(3,3)

    integer:: ierr,ic
    real(8):: alen(3),dlim

    if( nnode.gt.1 ) then
      print *,'ERROR: contraints does not work in parallel mode.'
      call mpi_finalize(ierr)
      stop
    endif

    if( myid.eq.0 .and. iprint.ne.0 ) then
      print *,''
      print *,'CONSTRAINTS parameters:'
    endif
    
    call read_const_params(myid,mpi_world,iprint)

    do ic=1,nconst
      if( constype(ic).ne.'bond' .and. constype(ic).ne.'dxyz' ) then
        print *,'ERROR: No such contraint available: '//trim(constype(ic))
        stop
      endif
    enddo

    alen(1) = norm2(h(:,1))
    alen(2) = norm2(h(:,2))
    alen(3) = norm2(h(:,3))
    dlim = min(alen(1),alen(2))
    dlim = min(dlim,alen(3))
    do ic=1,nconst
      if( trim(constype(ic)).eq.'bond' ) then
        if( dini(ic).gt.dlim/2 .or. dfin(ic).gt.dlim/2 ) then
          print *,'ERROR: A bond constraint too long w.r.t. '&
               //'simulation cell.'
          print *,'       Bond distance should be shorter than ' &
               //'a half of the cell size,'
          print '(a,3f7.3)','        where dini,dfin,dlim/2 = ', &
               dini(ic),dfin(ic),dlim/2
          stop 2
        endif
      else if( constype(ic) == 'dxyz' ) then
        if( abs(dxyzi(ic)) > alen(idxyz(ic))/2 &
             .or. abs(dxyzf(ic)) > alen(idxyz(ic))/2 ) then
          print *,'ERROR: A dxyz constraint too long w.r.t. cell vector,'
          print '(a,3f7.3)','        where dxyzi,dxyzf,alen/2 = ',  &
               dxyzi(ic),dxyzf(ic),alen(idxyz(ic))/2
          stop 2
        endif
      endif
    enddo

    if( myid.eq.0 .and. iprint.ne.0 ) then
      do ic=1,nconst
        if( constype(ic) == 'bond' ) then
          print '(3x,a,2(2x,i0),2(2x,f0.3))', 'bond: i, j, dini, dfin = ', &
               idcs0(1,ic),idcs0(2,ic) ,dini(ic),dfin(ic)
        else if( constype(ic) == 'dxyz' ) then
          print '(3x,a,3(2x,i0),2(2x,f0.3))', 'dxyz: i, j, ixyz, dini, dfin = ', &
               idcs0(1,ic),idcs0(2,ic), idxyz(ic), dxyzi(ic),dxyzf(ic)
        endif
      enddo
      print '(a,i0)','   maxiter = ',maxiter
    endif
    return
  end subroutine init_const
!=======================================================================
  subroutine read_const_params(myid,mpi_world,iprint)
    integer,intent(in):: myid,mpi_world,iprint
!!$    integer,external:: num_data

    integer:: i,j,ic,ixyz
    real(8):: d,di,df,xyzi,xyzf
    character(len=128):: c1st,cline,fname,ctmp
    character(len=20):: cxi,cyi,czi,cxf,cyf,czf
    
    if( myid.eq.0 ) then
      fname = trim(paramsdir)//'/'//trim(cfparams)
      open(ioprms,file=trim(fname),status='old')

!.....1st, count number of constraints
      nconst = 0
      do while(.true.)
        read(ioprms,*,end=10) c1st
        if( c1st(1:1).eq.'!' .or. c1st(1:1).eq.'#' ) cycle
        if( trim(c1st).eq.'const' ) nconst= nconst +1
      enddo
10    continue
      rewind(ioprms)
      print '(a,2x,i0)','   Number of constraints = ',nconst
      allocate(dij(nconst),dij2(nconst),dini(nconst),dfin(nconst) &
           ,idcs0(2,nconst),idcs(2,nconst),dtol2(nconst) &
           ,rio(3,nconst),rjo(3,nconst) &
           ,ris(3,nconst),rjs(3,nconst),vis(3,nconst),vjs(3,nconst) &
           ,rijo(3,nconst), dxyzi(nconst), dxyzf(nconst), &
           dxyz(nconst), idxyz(nconst), constype(nconst))

      ic = 0  ! initialize index of constraint
      do while(.true.)
        read(ioprms,'(a)',end=99) cline
        if( num_data(cline,' ').eq.0 ) cycle
        if( cline(1:1).eq.'!' .or. cline(1:1).eq.'#' ) cycle
        read(cline,*) c1st, ctype_const
        if( trim(c1st).eq.'const' ) then  ! entry of a constraint
          ic = ic + 1
          if( trim(ctype_const).eq.'bond' ) then
            constype(ic) = 'bond'
            read(cline,*) ctmp, ctmp, i, j, di, df
            idcs0(1,ic) = i
            idcs0(2,ic) = j
            dij(ic) = di
            dini(ic) = di
            dfin(ic) = df
          else if( trim(ctype_const).eq.'dxyz' ) then
            constype(ic) = 'dxyz'
            read(cline,*) ctmp, ctmp, i, j, ixyz, xyzi, xyzf
            idcs0(1,ic) = i
            idcs0(2,ic) = j
            idxyz(ic) = ixyz
            dxyzi(ic) = xyzi
            dxyzf(ic) = xyzf
          endif
        else if( trim(c1st).eq.'max_iteration' ) then
          read(cline,*) ctmp, maxiter
        else if( trim(c1st).eq.'tolerance' ) then
          read(cline,*) ctmp, tol
        else if( trim(c1st).eq.'tolerance_velocity' ) then
          read(cline,*) ctmp, vtol
        else
          print *,'There is no CONSTRAINTS keyword: ',trim(c1st)
        endif
      enddo
          
!!$!.....Num of bonds and bond pairs are read in the following lines
!!$          if( trim(ctype_const).eq.'bonds' ) then
!!$            read(ioprms,*,end=99) nconst
!!$            do ic=1,nconst
!!$              read(ioprms,'(a)') cline
!!$              if( num_data(cline,' ').eq.3 ) then  ! constant distance
!!$                backspace(ioprms)
!!$                read(ioprms,*)  i, j, d
!!$                idcs0(1,ic) = i
!!$                idcs0(2,ic) = j
!!$                dij(ic) = d
!!$                dini(ic) = d
!!$                dfin(ic) = d
!!$                dij2(ic) = d**2
!!$              else if( num_data(cline,' ').eq.4 ) then  ! linearly varying distance
!!$                backspace(ioprms)
!!$                read(ioprms,*) i, j, di, df
!!$                idcs0(1,ic) = i
!!$                idcs0(2,ic) = j
!!$                dij(ic) = di
!!$                dini(ic) = di
!!$                dfin(ic) = df
!!$              else
!!$                print *,'Error: Wrong number of parameters for bond !'
!!$              endif
!!$            enddo
!!$          endif
!!$        else if( trim(c1st).eq.'max_iteration' ) then
!!$          backspace(ioprms)
!!$          read(ioprms,*) c1st, maxiter
!!$        else if( trim(c1st).eq.'tolerance' ) then
!!$          backspace(ioprms)
!!$          read(ioprms,*) c1st, tol
!!$        else if( trim(c1st).eq.'tolerance_velocity' ) then
!!$          backspace(ioprms)
!!$          read(ioprms,*) c1st, vtol
!!$        else
!!$          print *,'There is no CONSTRAINTS keyword: ',trim(c1st)
!!$        endif
!!$      enddo

99    close(ioprms)
    endif

    return
  end subroutine read_const_params
!=======================================================================
  subroutine get_indice(namax,natm,tag)
!
!  Convert itot of constrained atoms (k1,k2) to current indice (i1,i2).
!  NOTE: This is not applicable to parallel MD...
!
    integer,intent(in):: namax,natm
    real(8),intent(in):: tag(namax)
!!$    integer,external:: itotOf
    integer:: ia,ic
    
!.....At first, find out two atoms that are constrained
    idcs(1:2,1:nconst) = -1
    do ia=1,natm
      do ic=1,nconst
        if( itotOf(tag(ia)).eq.idcs0(1,ic) ) idcs(1,ic) = ia
        if( itotOf(tag(ia)).eq.idcs0(2,ic) ) idcs(2,ic) = ia
      enddo
!.....If you want to skip searching once all the indices are set,
!     write a code for that hereafter...
    enddo
!!$    do ic=1,nconst
!!$      print *,'ic,idcs0,idcs=',ic,idcs0(1:2,ic),idcs(1:2,ic)
!!$    enddo
    return
  end subroutine get_indice
!=======================================================================
  subroutine update_const(namax,natm,tag,ra,h,istp,maxstp)
    integer,intent(in):: namax,natm,istp,maxstp
    real(8),intent(in):: tag(namax),ra(1:3,namax),h(3,3)
!!$    integer,external:: itotOf

    integer:: ic,i,j
    real(8):: ri(3),rj(3)

    call get_indice(namax,natm,tag)
    do ic=1,nconst
      if( constype(ic) == 'bond' ) then  ! bond
        if( maxstp.lt.1 ) then
          dij(ic) = dini(ic)
        else
          dij(ic) = dini(ic) +(dfin(ic)-dini(ic))*istp/(maxstp-1)
        endif
        dij2(ic)= dij(ic)**2
        dtol2(ic) = dij(ic)*tol
        dtol2(ic) = dtol2(ic)**2
        i = idcs(1,ic)
        j = idcs(2,ic)
        rio(:,ic) = ra(:,i)
        rjo(:,ic) = ra(:,j)
      else if( constype(ic) == 'dxyz' ) then  ! dxyz
        if( maxstp.lt.1 ) then
          dxyz(ic) = dxyzi(ic)
        else
          dxyz(ic) = dxyzi(ic) +(dxyzf(ic)-dxyzi(ic))*istp/(maxstp-1)
        endif
        i = idcs(1,ic)
        j = idcs(2,ic)
        rio(:,ic) = ra(:,i)
        rjo(:,ic) = ra(:,j)
      endif
    enddo

    return
  end subroutine update_const
!=======================================================================
  subroutine update_const_pos(namax,natm,h,hi,tag,ra,va,dt,nspmax,am)
!
!  Constraints on positions
!
    integer,intent(in):: namax,natm,nspmax
    real(8),intent(in):: h(3,3),hi(3,3),tag(namax),dt,am(nspmax)
    real(8),intent(inout):: ra(3,namax),va(3,namax)

    integer:: i,j,ia,ic,is,js,iter,ixyz,jxyz
    real(8):: rij(3),rijos(3),vij(3),ami,amj,amij,dd,gmk,dijtmp(3)
    logical:: not_conv

!.....Normal velocity Verlet update of positions
    do ia=1,natm
      ra(1:3,ia) = ra(1:3,ia) +(hi(1:3,1)*va(1,ia) &
           +hi(1:3,2)*va(2,ia) +hi(1:3,3)*va(3,ia))*dt
    enddo

!.....Constraints hereafter
    do ic=1,nconst
      i = idcs(1,ic)
      j = idcs(2,ic)
      rio(1:3,ic) = ra(1:3,i)
      rjo(1:3,ic) = ra(1:3,j)
      rij(1:3) = rjo(1:3,ic) -rio(1:3,ic)
      rij(1:3) = rij(1:3) -anint(rij(1:3))
      rijo(:,ic) = abc2cart(h,rij)
    enddo
    not_conv = .false. 
    do ic=1,nconst
      i = idcs(1,ic)
      j = idcs(2,ic)
      ris(1:3,ic) = ra(1:3,i)
      rjs(1:3,ic) = ra(1:3,j)
      vis(1:3,ic) = va(1:3,i)
      vjs(1:3,ic) = va(1:3,j)
      rij(1:3) = rjs(1:3,ic) -ris(1:3,ic)
      rij(1:3) = rij(1:3) -anint(rij(1:3))
      rij = abc2cart(h,rij)
      if( constype(ic) == 'bond' ) then
        dd = norm2(rij)
        if( abs(dd-dij2(ic)).gt.dtol2(ic) ) not_conv = .true. 
      else if( constype(ic) == 'dxyz' ) then
        if( abs(rij(idxyz(ic))-dxyz(ic)) > tol ) not_conv = .true.
      endif
    enddo
    if( .not. not_conv ) goto 10

    do iter=1,maxiter
!!$      if( .not. not_conv ) exit
      not_conv = .false.
      do ic=1,nconst
        i = idcs(1,ic)
        j = idcs(2,ic)
        is = int(tag(i))
        js = int(tag(j))
        ami = am(is)
        amj = am(js)
        amij = ami*amj/(ami+amj)
        rij(1:3) = rjs(1:3,ic) -ris(1:3,ic)
        rij(1:3) = rij(1:3) -anint(rij(1:3))
        rij = abc2cart(h,rij)
        dd = norm2(rij)
        rijos = cart2abc(hi,rijo)
        if( constype(ic) == 'bond' ) then
          gmk = amij *(dd-dij2(ic)) /dot(rij,rijo(:,ic))
          ris(1:3,ic) = ris(1:3,ic) +gmk/(2d0*ami)*rijos(1:3)
          rjs(1:3,ic) = rjs(1:3,ic) -gmk/(2d0*amj)*rijos(1:3)
          vis(1:3,ic) = vis(1:3,ic) +gmk/(2d0*ami)*rijo(1:3,ic) /dt
          vjs(1:3,ic) = vjs(1:3,ic) -gmk/(2d0*amj)*rijo(1:3,ic) /dt
!.....Check convergence
          rij(1:3) = rjs(1:3,ic) -ris(1:3,ic)
          rij(1:3) = rij(1:3) -anint(rij(1:3))
          rij = abc2cart(h,rij)
          dd = norm2(rij)
!!$        print '(a,i6,3es12.4)','iter,gmk,dd-dij2,dtol2=',iter,gmk,abs(dd-dij2(ic)),dtol2(ic)
          if( abs(dd-dij2(ic)).gt.dtol2(ic) ) not_conv = .true. 
        else if( constype(ic) == 'dxyz' ) then
          ixyz = idxyz(ic)
          gmk = 2d0*amij *(rij(ixyz)-dxyz(ic))
          dijtmp(:) = 0d0
          do jxyz=1,3
            dijtmp(jxyz) = dijtmp(jxyz) +hi(jxyz,ixyz)*gmk
          enddo
          ris(:,ic) = ris(:,ic) +dijtmp(:)/(2d0*ami)
          rjs(:,ic) = rjs(:,ic) -dijtmp(:)/(2d0*amj)
          vis(ixyz,ic) = vis(ixyz,ic) +gmk/(2d0*ami) /dt
          vjs(ixyz,ic) = vjs(ixyz,ic) -gmk/(2d0*amj) /dt
!.....Check convergence
          rij(1:3) = rjs(1:3,ic) -ris(1:3,ic)
          rij(1:3) = rij(1:3) -anint(rij(1:3))
          rij = abc2cart(h,rij)
          if( abs(rij(ixyz)-dxyz(ic)).gt.tol ) not_conv = .true. 
        endif
      enddo
      if( .not. not_conv ) then
!!$        print *,' converged at ',iter
        goto 10
      endif
    enddo
!!$    print *,' not converged'
10  continue
    do ic=1,nconst
      i = idcs(1,ic)
      j = idcs(2,ic)
      ra(1:3,i) = ris(1:3,ic)
      ra(1:3,j) = rjs(1:3,ic)
      va(1:3,i) = vis(1:3,ic)
      va(1:3,j) = vjs(1:3,ic)
    enddo

    return
  end subroutine update_const_pos
!=======================================================================
  subroutine update_const_vel(namax,natm,h,hi,tag,va,dt,nspmax,am)
!
!  Update velocities of atoms involved by the contraints
!
    integer,intent(in):: namax,natm,nspmax
    real(8),intent(in):: h(3,3),hi(3,3),tag(namax),dt,am(nspmax)
    real(8),intent(inout):: va(3,namax)

    integer:: ic,i,j,iter,is,js,ixyz
    real(8):: vij(3),rij(3),sgm,ami,amj,amij,gmk
    logical:: not_conv

    do ic=1,nconst
      i = idcs(1,ic)
      j = idcs(2,ic)
      vis(1:3,ic) = va(1:3,i)
      vjs(1:3,ic) = va(1:3,j)
    enddo

    do iter=1,maxiter
      not_conv = .false.
      do ic=1,nconst
        vij(1:3) = vjs(1:3,ic) -vis(1:3,ic)
        rij(1:3) = rjs(1:3,ic) -ris(1:3,ic)
!!$        rij(1:3) = rjo(1:3,ic) -rio(1:3,ic)
        rij(1:3) = rij(1:3) -anint(rij(1:3))
        rij = abc2cart(h,rij)
        if( constype(ic) == 'bond' ) then
          sgm = dot(vij,rij)
          if( abs(sgm).gt.vtol ) then
            not_conv = .true.
          else
            cycle
          endif
          i = idcs(1,ic)
          j = idcs(2,ic)
          is = int(tag(i))
          js = int(tag(j))
          ami = am(is)
          amj = am(js)
          amij = ami*amj /(ami+amj)
          gmk = amij /dij2(ic)*sgm
          vis(1:3,ic) = vis(1:3,ic) +gmk/ami*rij(1:3)
          vjs(1:3,ic) = vjs(1:3,ic) -gmk/amj*rij(1:3)
        else if( constype(ic) == 'dxyz' ) then
          ixyz = idxyz(ic)
          sgm = vij(ixyz)
          if( abs(sgm) > vtol ) then
            not_conv = .true.
          else
            cycle
          endif
          i = idcs(1,ic)
          j = idcs(2,ic)
          is = int(tag(i))
          js = int(tag(j))
          ami = am(is)
          amj = am(js)
          amij = ami*amj /(ami+amj)
          gmk = amij *vij(ixyz)
          vis(ixyz,ic) = vis(ixyz,ic) +gmk/ami
          vjs(ixyz,ic) = vjs(ixyz,ic) -gmk/amj
        endif
      enddo
      if( .not. not_conv ) exit
    enddo

    do ic=1,nconst
      i = idcs(1,ic)
      j = idcs(2,ic)
      va(1:3,i) = vis(1:3,ic)
      va(1:3,j) = vjs(1:3,ic)
    enddo

    return
  end subroutine update_const_vel
end module constraints
!-----------------------------------------------------------------------
!  Local Variables:
!  compile-command: "make pmd"
!  End:
