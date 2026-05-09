module rdcfrc
  use pmdmpi
  use mod_precision
  implicit none
  save

  character(len=128):: paramsdir = '.'
  character(len=128),parameter:: cfparams = 'in.rdcfrc'
  integer,parameter:: ioprms = 67

  real(rp),parameter:: pi = 3.14159265358979_rp

!.....Type of selection of atoms whose forces are to be reduced.
!.....Available types: species, ID, neigh_XXXXX, cna
  character(len=128):: ctype_fmod = 'species' ! default
  
!.....alpha: factor to be multiplied to forces
  real(rp):: alpha = 0.5_rp

!.....Species whose forces are to be reduced.
  integer:: is_rdcfrc = -1 ! default
  integer,allocatable:: ids_fmod(:)
  integer:: nids = 0

!.....Number of neighbors per atom
  real(rp),allocatable:: ann(:)
!.....Distance where neighbor counting starts to decrease
  real(rp):: rin = -1.0_rp
!.....Cutoff distance of neighbor counting
  real(rp):: rout = -1.0_rp
  real(rp):: rout2
!.....Forces of atoms having this num of neighbors are not modified
  integer:: nn_intact = 8
!.....Threshold of difference in num of neighbors
  real(rp):: dth_nn = 0.5_rp
!.....Boundary in num of neighbors
  real(rp):: nn_boundary = -1.0_rp

!.....For the calculation of histogram of num of neighbors
  logical:: lhist = .false.
  integer:: nbins = 200
  real(rp),allocatable:: hist_nn(:)
  real(rp):: dhist
  
contains
!=======================================================================
  subroutine init_rdcfrc(myid,mpi_world,nnode,iprint)
    integer,intent(in):: myid,mpi_world,nnode,iprint

    integer:: ierr,i
    character(len=128):: cnum
    
    call read_rdcfrc_params(myid,mpi_world,iprint)

    call mpi_bcast(ctype_fmod,128,mpi_character,0,mpi_world,ierr)
    call mpi_bcast(is_rdcfrc,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(alpha,1,mpi_real_rp,0,mpi_world,ierr)
    call mpi_bcast(rin,1,mpi_real_rp,0,mpi_world,ierr)
    call mpi_bcast(rout,1,mpi_real_rp,0,mpi_world,ierr)
    call mpi_bcast(dth_nn,1,mpi_real_rp,0,mpi_world,ierr)
    call mpi_bcast(nn_boundary,1,mpi_real_rp,0,mpi_world,ierr)
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
      if( rin.lt.0.0_rp .or. rout.lt.0.0_rp .or. rout.lt.rin ) then
        print *,'ERROR: rin and/or rout are wrong, rin,rout=',rin,rout
        stop
      endif
      if( ctype_fmod(7:11).eq.'under' .and. nn_boundary.lt.0.0_rp ) then
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
      hist_nn(:) = 0.0_rp
      dhist = 13.0_rp /nbins
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
      else if( trim(ctype_fmod).eq.'CNA' ) then
        if( is_rdcfrc.ge.0 ) then
          print '(a,i0)','  CNA-ID whose forces are to be reduced: ', is_rdcfrc
        else
          print '(a,i0)','  Forces of all atoms are to be reduced, ' &
               //'except those of CNA-ID = ',abs(is_rdcfrc)
        endif
      else  ! default: species
        print '(a,i3)','  species: ',is_rdcfrc
      endif
    endif
    return
  end subroutine init_rdcfrc
!=======================================================================
  subroutine read_rdcfrc_params(myid,mpi_world,iprint)
    use util, only: num_data
    integer,intent(in):: myid,mpi_world,iprint
    
!!$    integer,external:: num_data
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
        else if( trim(c1st).eq.'cna_id' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, is_rdcfrc
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
  subroutine reduce_forces(namax,natm,aa,tag_isp,tag_itot,ra,h,nnmax,lspr)
!
!  Reduce forces on atoms that are selected in some way...
!
    use structure,only: idcna
    integer,intent(in):: namax,natm,nnmax,lspr(0:nnmax,namax)
    integer,intent(in):: tag_isp(namax),tag_itot(namax)
    real(rp),intent(in):: ra(3,namax),h(3,3)
    real(rp),intent(inout):: aa(3,namax)

    integer:: ia,is,i,itot,icna
    real(rp):: beta


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
          beta = 1.0_rp
          if( ann(ia).lt. nn_boundary-dth_nn ) then
            beta = alpha
          else if( ann(ia).lt. nn_boundary+dth_nn ) then
            beta = alpha +(1.0_rp-alpha)*0.5_rp*(1.0_rp &
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
        itot = tag_itot(ia)
        do is=1,nids
          if( itot.eq.ids_fmod(is) ) then
            aa(1:3,ia) = alpha *aa(1:3,ia)
            exit
          endif
        enddo
      enddo
    else if( trim(ctype_fmod).eq.'CNA' ) then
      if( .not. allocated(idcna) ) then
        print *,'ERROR @reduce_forces: idcna is not allocated !!!'
        stop 1
      endif
      if( is_rdcfrc.ge.0 ) then
        do ia=1,natm
          icna = idcna(ia)
          if( icna.eq.is_rdcfrc ) aa(1:3,ia) = alpha *aa(1:3,ia)
        enddo
      else
        do ia=1,natm
          icna = idcna(ia)
          if( icna.ne.abs(is_rdcfrc) ) aa(1:3,ia) = alpha *aa(1:3,ia)
        enddo
      endif
    else ! default: species
      do ia=1,natm
        is = tag_isp(ia)
        if( is.eq.is_rdcfrc ) aa(1:3,ia) = alpha *aa(1:3,ia)
      enddo
    endif
    return
  end subroutine reduce_forces
!=======================================================================
  subroutine count_nn(namax,natm,ra,h,nnmax,lspr)
    integer,intent(in):: namax,natm,nnmax,lspr(0:nnmax,namax)
    real(rp),intent(in):: ra(3,namax),h(3,3)

    integer:: ia,jj,ja
    real(rp):: xi(3),xj(3),xij(3),rij(3),dij2,dij

    ann(1:natm) = 0.0_rp
    do ia=1,natm
      xi(1:3) = ra(1:3,ia)
      do jj=1,lspr(0,ia)
        ja = lspr(jj,ia)
        xj(1:3) = ra(1:3,ja)
        xij(1:3) = xj(1:3) -xi(1:3)
        rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        dij2 = rij(1)*rij(1) +rij(2)*rij(2) +rij(3)*rij(3)
        if( dij2.gt.rout2 ) cycle
        dij = sqrt(dij2)
        ann(ia) = ann(ia) +fnn(dij)
      enddo
    enddo
    
  end subroutine count_nn
!=======================================================================
  function fnn(rij)
    real(rp),intent(in):: rij
    real(rp):: fnn

    fnn = 0.0_rp
    if( rij.le.rin ) then
      fnn = 1.0_rp
    else if( rin.lt.rij .and. rij.le.rout ) then
      fnn = 0.5_rp *(1.0_rp +cos(pi*(rij-rin)/(rout-rin)))
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
        write(90,'(i5,2es15.7)') i,dhist*(real(i, rp)-0.5),hist_nn(i)
      enddo
      close(90)
    endif
  end subroutine finalize_rdcfrc
end module rdcfrc
!-----------------------------------------------------------------------
!  Local Variables:
!  compile-command: "make pmd"
!  End:
