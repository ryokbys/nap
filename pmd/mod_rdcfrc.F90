module rdcfrc
  implicit none
  save
  include 'mpif.h'

  character(len=128):: paramsdir = '.'
  character(len=128),parameter:: cfparams = 'in.rdcfrc'
  integer,parameter:: ioprms = 67

  real(8),parameter:: pi = 3.14159265358979d0

!.....Type of selection of atoms whose forces are to be reduced.
!.....Available types: species, ID, neigh_XXXXX
  character(len=128):: ctype_fmod = 'species' ! default
  
!.....alpha: factor to be multiplied to forces
  real(8):: alpha = 0.5d0

!.....Species whose forces are to be reduced.
  integer:: is_rdcfrc = -1 ! default
  integer,allocatable:: ids_fmod(:)
  integer:: nids = 0

!.....Number of neighbors per atom
  real(8),allocatable:: ann(:)
!.....Distance where neighbor counting starts to decrease
  real(8):: rin = -1d0
!.....Cutoff distance of neighbor counting
  real(8):: rout = -1d0
  real(8):: rout2
!.....Forces of atoms having this num of neighbors are not modified
  integer:: nn_intact = 8
!.....Threshold of difference in num of neighbors
  real(8):: dth_nn = 0.5d0
!.....Boundary in num of neighbors
  real(8):: nn_boundary = -1d0

!.....For the calculation of histogram of num of neighbors
  logical:: lhist = .false.
  integer:: nbins = 200
  real(8),allocatable:: hist_nn(:)
  real(8):: dhist
  
contains
!=======================================================================
  subroutine init_rdcfrc(myid,mpi_world,nnode,iprint)
    integer,intent(in):: myid,mpi_world,nnode,iprint

    integer:: ierr,i
    character(len=128):: cnum
    
    call read_rdcfrc_params(myid,mpi_world,iprint)

    call mpi_bcast(ctype_fmod,128,mpi_character,0,mpi_world,ierr)
    call mpi_bcast(is_rdcfrc,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(alpha,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(rin,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(rout,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(dth_nn,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(nn_boundary,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(nn_intact,1,mpi_integer,mpi_world,ierr)
    call mpi_bcast(nids,1,mpi_integer,0,mpi_world,ierr)
    if( nids.gt.0 ) then
      if( myid.ne.0 ) allocate(ids_fmod(nids))
      call mpi_bcast(ids_fmod,nids,mpi_integer,0,mpi_world,ierr)
    endif
    call mpi_bcast(lhist,1,mpi_logical,0,mpi_world,ierr)
    call mpi_bcast(nbins,1,mpi_integer,0,mpi_world,ierr)
    

!.....Error handling
    if( trim(ctype_fmod).eq.'species' .and. is_rdcfrc.lt.1 ) then
      print *,'ERROR: species is not specified in rdcfrc, is_rdcfrc=',is_rdcfrc
      stop
    else if( trim(ctype_fmod).eq.'id' .and. nids.lt.1 ) then
      print *,'ERROR: IDs should be specified.'
      stop
    else if( ctype_fmod(1:5).eq.'neigh' ) then
      if( rin.lt.0d0 .or. rout.lt.0d0 .or. rout.lt.rin ) then
        print *,'ERROR: rin and/or rout are wrong, rin,rout=',rin,rout
        stop
      endif
      if( ctype_fmod(7:11).eq.'under' .and. nn_boundary.lt.0d0 ) then
        print *,'ERROR: nn_boundary should be specified, nn_boundary = ',nn_boundary
        stop
      endif
    endif

    if( ctype_fmod(1:5).eq.'neigh' ) then
      rout2 = rout*rout
    endif

    if( lhist ) then
      if( nbins.lt.1 ) then
        print *,'ERROR: nbins < 1, nbins=',nbins
        stop
      endif
      allocate(hist_nn(nbins))
      hist_nn(:) = 0d0
      dhist = 13d0 /nbins
    endif

    if( myid.eq.0 ) then
      print *,'RDCFRC parameters:'
      print '(a,a)','  type:    ',trim(ctype_fmod)
      print '(a,f6.3)','  alpha:   ',alpha
      if( ctype_fmod(1:5).eq.'neigh' ) then
        print '(a,2f7.3)','  rin,rout:   ',rin,rout
        print '(a,f6.2)','  dth_nn:   ',dth_nn
        if( ctype_fmod(7:11).eq.'under' ) then
          print '(a,f7.3)','  nn_boundary:   ',nn_boundary
        endif
      else if( trim(ctype_fmod).eq.'ID' ) then
        print '(a)','  IDs of atoms whose forces are to be modified: '
        write(cnum,'(i0)') nids
        print '(a,'//trim(cnum)//'(x,i0))', '  ',(ids_fmod(i),i=1,nids)
      else  ! default: species
        print '(a,i3)','  species: ',is_rdcfrc
      endif
    endif
    return
  end subroutine init_rdcfrc
!=======================================================================
  subroutine read_rdcfrc_params(myid,mpi_world,iprint)
    integer,intent(in):: myid,mpi_world,iprint
    
    integer,external:: num_data
    character(len=128):: cfname,c1st
    integer:: i
    
    if( myid.eq.0 ) then
      cfname = trim(paramsdir)//'/'//trim(cfparams)
      open(ioprms,file=trim(cfname),status='old')
      do while(.true.)
        read(ioprms,*,end=10) c1st
        if( num_data(c1st,' ').eq.0 ) cycle
        if( c1st(1:1).eq.'!' .or. c1st(1:1).eq.'#' ) cycle
        if( trim(c1st).eq.'type' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, ctype_fmod
        else if( trim(c1st).eq.'species' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, is_rdcfrc
        else if( trim(c1st).eq.'IDs' ) then
          backspace(ioprms)
          read(ioprms,'(a)') c1st
          nids = num_data(c1st,' ') -1
          allocate(ids_fmod(nids))
          backspace(ioprms)
          read(ioprms,*) c1st, (ids_fmod(i),i=1,nids)
        else if( trim(c1st).eq.'alpha' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, alpha
        else if( trim(c1st).eq.'rin' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, rin
        else if( trim(c1st).eq.'rout' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, rout
        else if( trim(c1st).eq.'dnn_threshold' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, dth_nn
        else if( trim(c1st).eq.'nn_boundary' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, nn_boundary
        else if( trim(c1st).eq.'nn_intact' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, nn_intact
        else if( trim(c1st).eq.'histogram' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, lhist
        else if( trim(c1st).eq.'nbins' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, nbins
        else
          print *,'No such RDCFRC keyword: ',trim(c1st)
        endif
      enddo
      
10    close(ioprms)
    endif ! myid.eq.0
    return
  end subroutine read_rdcfrc_params
!=======================================================================
  subroutine reduce_forces(namax,natm,aa,tag,ra,h,nnmax,lspr)
!
!  Reduce forces on atoms that are selected in some way...
!
    integer,intent(in):: namax,natm,nnmax,lspr(0:nnmax,namax)
    real(8),intent(in):: ra(3,namax),tag(namax),h(3,3)
    real(8),intent(inout):: aa(3,namax)

    integer:: ia,is,i,itot
    real(8):: beta
    integer,external:: itotOf

    if( ctype_fmod(1:5).eq.'neigh' ) then
      if( .not. allocated(ann) ) allocate(ann(namax))
      call count_nn(namax,natm,ra,h,nnmax,lspr)
      if( lhist ) then
        do ia=1,natm
          i = int(ann(ia)/dhist)
          if( i.le.nbins ) then
            hist_nn(i) = hist_nn(i) + 1.0
          endif
        enddo
      endif
      if( ctype_fmod(7:11).eq.'under' ) then
        do ia=1,natm
          beta = 1d0
          if( ann(ia).lt. nn_boundary-dth_nn ) then
            beta = alpha
          else if( ann(ia).lt. nn_boundary+dth_nn ) then
            beta = alpha +(1d0-alpha)*0.5d0*(1d0 &
                 -cos(pi *(ann(ia)-(nn_boundary-dth_nn)) /(2*dth_nn)))
          endif
!!$        print *,'ia,ann(ia),alpha,beta=',ia,ann(ia),alpha,beta
!!$        beta = 1d0 -(1d0-alpha)*exp(-(ann(ia)-7.0)**2/dth_nn**2/2)
!!$        print '(a,i6,3f8.3)','ia,ann,beta,exp=',ia,ann(ia),beta&
!!$             ,exp(-(ann(ia)-7.0)**2/dth_nn**2/2)
          aa(1:3,ia) = aa(1:3,ia) *beta
!!$        aa(1:3,ia) = aa(1:3,ia)*(alpha  &
!!$             +(1d0-alpha)*exp(-(ann(ia)-nn_intact)**2/dth_nn**2/2))
!!$        if( abs(ann(ia) -nn_intact).gt.dth_nn ) then
!!$          aa(1:3,ia) = alpha *aa(1:3,ia)
!!$        endif
        enddo
      endif
    else if( trim(ctype_fmod).eq.'ID' ) then
      do ia=1,natm
        itot = itotOf(tag(ia))
        do is=1,nids
          if( itot.eq.ids_fmod(is) ) then
            aa(1:3,ia) = alpha *aa(1:3,ia)
            exit
          endif
        enddo
      enddo
    else ! default: species
      do ia=1,natm
        is = int(tag(ia))
        if( is.eq.is_rdcfrc ) aa(1:3,ia) = alpha *aa(1:3,ia)
      enddo
    endif
    return
  end subroutine reduce_forces
!=======================================================================
  subroutine count_nn(namax,natm,ra,h,nnmax,lspr)
    integer,intent(in):: namax,natm,nnmax,lspr(0:nnmax,namax)
    real(8),intent(in):: ra(3,namax),h(3,3)

    integer:: ia,jj,ja
    real(8):: xi(3),xj(3),xij(3),rij(3),dij2,dij

    ann(1:natm) = 0d0
    do ia=1,natm
      xi(1:3) = ra(1:3,ia)
      do jj=1,lspr(0,ia)
        ja = lspr(jj,ia)
        xj(1:3) = ra(1:3,ja)
        xij(1:3) = xj(1:3) -xi(1:3)
        rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij2 = rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
        if( dij2.gt.rout2 ) cycle
        dij = dsqrt(dij2)
        ann(ia) = ann(ia) +fnn(dij)
      enddo
    enddo
    
  end subroutine count_nn
!=======================================================================
  function fnn(rij)
    real(8),intent(in):: rij
    real(8):: fnn

    fnn = 0d0
    if( rij.le.rin ) then
      fnn = 1d0
    else if( rin.lt.rij .and. rij.le.rout ) then
      fnn = 0.5d0 *(1d0 +cos(pi*(rij-rin)/(rout-rin)))
    endif
    return
  end function fnn
!=======================================================================
  subroutine finalize_rdcfrc()
    integer:: i

    if( lhist ) then
!.....Output histogram
      open(90,file='out.hist_rdcfrc',status='replace')
      do i=1,nbins
        write(90,'(i5,2es15.7)') i,dhist*(dble(i)-0.5),hist_nn(i)
      enddo
      close(90)
    endif
  end subroutine finalize_rdcfrc
end module rdcfrc
!-----------------------------------------------------------------------
!  Local Variables:
!  compile-command: "make pmd"
!  End:
