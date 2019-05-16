module constraints
!
!  Currently only one bond distant constraint is available.
!
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
  integer,allocatable:: idcs0(:,:),idcs(:,:)
  real(8),allocatable:: rio(:,:),rjo(:,:)
  real(8),allocatable:: ris(:,:),rjs(:,:),vis(:,:),vjs(:,:)
  real(8),allocatable:: rijo(:,:)
  real(8):: tol = 1d-2
  real(8):: vtol = 1d-8
  
  
contains
!=======================================================================
  subroutine init_const(myid,mpi_world,nnode,iprint,h)
    use vector
    integer,intent(in):: myid,mpi_world,nnode,iprint
    real(8),intent(in):: h(3,3)

    integer:: ierr,ic
    real(8):: alen(3),dlim

    if( nnode.gt.1 ) then
      print *,'ERROR: contraints does not work in parallel mode.'
      call mpi_finalize(ierr)
      stop
    endif

    call read_const_params(myid,mpi_world,iprint)

    if( trim(ctype_const).ne.'bonds' ) then
      print *,'ERROR: No such contraint available: '//trim(ctype_const)
      stop
    endif

    if( trim(ctype_const).eq.'bonds' ) then
      alen(1) = norm2(h(:,1))
      alen(2) = norm2(h(:,2))
      alen(3) = norm2(h(:,3))
      dlim = min(alen(1),alen(2))
      dlim = min(dlim,alen(3))
      do ic=1,nconst
        if( dfin(ic).gt.dlim/2 ) then
          print *,'ERROR: a bond constraint too long w.r.t. '&
               //'simulation cell.'
          print *,'       bond distance should be shorter than ' &
               //'a half of the cell size'
          stop 2
        endif
      enddo
    endif

    if( myid.eq.0 .and. iprint.ne.0 ) then
      print *,''
      print *,'CONSTRAINTS parameters:'
      print '(a)', '  Constraint type: '//trim(ctype_const)
      if( trim(ctype_const).eq.'bonds' ) then
        print '(a,i0)', '   Number of bonds = ',nconst
        print '(a)','    #i,  #j,  dini,  dfin'
        do ic=1,nconst
          print '(3x,2(2x,i0),2(2x,f0.3))', idcs0(1,ic),idcs0(2,ic) &
               ,dini(ic),dfin(ic)
        enddo
      endif
      print '(a,i0)','   maxiter = ',maxiter
    endif
    return
  end subroutine init_const
!=======================================================================
  subroutine read_const_params(myid,mpi_world,iprint)
    use util, only: num_data
    integer,intent(in):: myid,mpi_world,iprint
!!$    integer,external:: num_data

    integer:: i,j,ic
    real(8):: d,di,df
    character(len=128):: c1st,cline,fname
    
    if( myid.eq.0 ) then
      fname = trim(paramsdir)//'/'//trim(cfparams)
      open(ioprms,file=trim(fname),status='old')

      do while(.true.)
        read(ioprms,*,end=10) c1st
        if( num_data(c1st,' ').eq.0 ) cycle
        if( c1st(1:1).eq.'!' .or. c1st(1:1).eq.'#' ) cycle
        if( trim(c1st).eq.'const_type' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, ctype_const
!.....Num of bonds and bond pairs are read in the following lines
          if( trim(ctype_const).eq.'bonds' ) then
            read(ioprms,*,end=10) nconst
            allocate(dij(nconst),dij2(nconst),dini(nconst),dfin(nconst) &
                 ,idcs0(2,nconst),idcs(2,nconst),dtol2(nconst) &
                 ,rio(3,nconst),rjo(3,nconst) &
                 ,ris(3,nconst),rjs(3,nconst),vis(3,nconst),vjs(3,nconst) &
                 ,rijo(3,nconst))
            do ic=1,nconst
              read(ioprms,'(a)') cline
              if( num_data(cline,' ').eq.3 ) then  ! constant distance
                backspace(ioprms)
                read(ioprms,*)  i, j, d
                idcs0(1,ic) = i
                idcs0(2,ic) = j
                dij(ic) = d
                dini(ic) = d
                dfin(ic) = d
                dij2(ic) = d**2
              else if( num_data(cline,' ').eq.4 ) then  ! linearly varying distance
                backspace(ioprms)
                read(ioprms,*) i, j, di, df
                idcs0(1,ic) = i
                idcs0(2,ic) = j
                dij(ic) = di
                dini(ic) = di
                dfin(ic) = df
              else
                print *,'Error: Wrong number of parameters for bond !'
              endif
            enddo
          endif
        else if( trim(c1st).eq.'max_iteration' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, maxiter
        else if( trim(c1st).eq.'tolerance' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, tol
        else if( trim(c1st).eq.'tolerance_velocity' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, vtol
        else
          print *,'There is no CONSTRAINTS keyword: ',trim(c1st)
        endif
      enddo

10    close(ioprms)
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
    integer,external:: itotOf
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
    use vector
    integer,intent(in):: namax,natm,istp,maxstp
    real(8),intent(in):: tag(namax),ra(1:3,namax),h(3,3)
    integer,external:: itotOf

    integer:: ic,i,j
    real(8):: ri(3),rj(3)

    call get_indice(namax,natm,tag)
    do ic=1,nconst
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
      rio(:,ic) = ra(1:3,i)
      rjo(:,ic) = ra(1:3,j)
!!$      print '(a,5i4,6f7.3)','ic,i0,j0,i,j,dij,dij2=',ic,idcs0(1:2,ic),i,j &
!!$           ,dij(ic),dij2(ic)
    enddo

    return
  end subroutine update_const
!=======================================================================
  subroutine update_const_pos(namax,natm,h,hi,tag,ra,va,dt,nspmax,am)
!
!  Constraints on positions
!
    use vector
    integer,intent(in):: namax,natm,nspmax
    real(8),intent(in):: h(3,3),hi(3,3),tag(namax),dt,am(nspmax)
    real(8),intent(inout):: ra(3,namax),va(3,namax)

    integer:: i,j,ia,ic,is,js,iter
    real(8):: rij(3),rijos(3),vij(3),ami,amj,amij,dd,gmk
    logical:: not_conv
    integer,external:: itotOf

!.....Normal velocity Verlet update of positions
    do ia=1,natm
      ra(1:3,ia) = ra(1:3,ia) +va(1:3,ia) *dt
    enddo

!.....Constraints hereafter
    do ic=1,nconst
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
      vis(1:3,ic) = abc2cart(h,va(:,i))
      vjs(1:3,ic) = abc2cart(h,va(:,j))
      rij(1:3) = rjs(1:3,ic) -ris(1:3,ic)
      rij(1:3) = rij(1:3) -anint(rij(1:3))
      rij = abc2cart(h,rij)
      dd = norm2(rij)
      if( abs(dd-dij2(ic)).gt.dtol2(ic) ) not_conv = .true. 
!!$      print '(a,5es11.3,2x,3l)','rij,dd,dij2,not_conv = ',rij(1:3),sqrt(dd),sqrt(dij2(ic)),not_conv
    enddo
    if( .not. not_conv ) goto 10

    do iter=1,maxiter
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
        gmk = amij *(dd-dij2(ic)) /dot(rij,rijo(:,ic))
        rijos = cart2abc(hi,rijo)
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
      enddo
      if( .not. not_conv ) goto 10
    enddo

10  continue
    do ic=1,nconst
      i = idcs(1,ic)
      j = idcs(2,ic)
      ra(1:3,i) = ris(:,ic)
      ra(1:3,j) = rjs(:,ic)
      va(1:3,i) = cart2abc(hi,vis(:,ic))
      va(1:3,j) = cart2abc(hi,vjs(:,ic))
    enddo

    return
  end subroutine update_const_pos
!=======================================================================
  subroutine update_const_vel(namax,natm,h,hi,tag,va,dt,nspmax,am)
!
!  Update velocities of atoms involved by the contraints
!
    use vector
    integer,intent(in):: namax,natm,nspmax
    real(8),intent(in):: h(3,3),hi(3,3),tag(namax),dt,am(nspmax)
    real(8),intent(inout):: va(3,namax)

    integer:: ic,i,j,iter,is,js
    real(8):: vij(3),rij(3),sgm,ami,amj,amij,gmk
    logical:: not_conv

    do ic=1,nconst
      i = idcs(1,ic)
      j = idcs(2,ic)
      vis(1:3,ic) = abc2cart(h,va(:,i))
      vjs(1:3,ic) = abc2cart(h,va(:,j))
    enddo

    do iter=1,maxiter
      not_conv = .false.
      do ic=1,nconst
        vij(1:3) = vjs(1:3,ic) -vis(1:3,ic)
        rij(1:3) = rjo(1:3,ic) -rio(1:3,ic)
        rij(1:3) = rij(1:3) -anint(rij(1:3))
        rij = abc2cart(h,rij)
        sgm = dot(vij,rij)
!!$        print *,'iter,vij,rij,sgm,vtol=',iter,norm2(vij),norm2(rij) &
!!$             ,abs(sgm),vtol
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
      enddo
      if( .not. not_conv ) exit
    enddo

    do ic=1,nconst
      i = idcs(1,ic)
      j = idcs(2,ic)
      va(1:3,i) = cart2abc(hi,vis(:,ic))
      va(1:3,j) = cart2abc(hi,vjs(:,ic))
    enddo

    return
  end subroutine update_const_vel
!=======================================================================
  subroutine store_ra_va(namax,natm,h,ra,va)
!
!  Store ra and va before updating with forces.
!
    use vector
    integer,intent(in):: namax,natm
    real(8),intent(in):: ra(1:3,namax),va(1:3,namax),h(3,3)
    integer:: ic,i,j
    real(8):: r(3)

    do ic=1,nconst
      i = idcs(1,ic)
      r(1:3) = ra(1:3,i)
      rio(:,ic) = abc2cart(h,r)
      j = idcs(2,ic)
      r(1:3) = ra(1:3,j)
      rjo(:,ic) = abc2cart(h,r)
    enddo
    return
  end subroutine store_ra_va
end module constraints
!-----------------------------------------------------------------------
!  Local Variables:
!  compile-command: "make pmd"
!  End:
