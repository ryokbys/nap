module frcmod
  implicit none
  save
  include 'mpif.h'

  character(len=128):: paramsdir = '.'
  character(len=128),parameter:: cfparams = 'in.frcmod'
  integer,parameter:: ioprms = 67

  real(8),parameter:: pi = 3.14159265358979d0

!.....Type of selection of atoms whose forces are to be reduced.
!.....Available types: species, neighbor
  character(len=128):: ctype_fmod = 'species' ! default
  
!.....alpha: factor to be multiplied to forces
  real(8):: alpha = 0.5d0

!.....Species whose forces are to be reduced.
  integer:: is_fmod = -1 ! default

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
  
contains
!=======================================================================
  subroutine init_frcmod(myid,mpi_world,nnode,iprint)
    integer,intent(in):: myid,mpi_world,nnode,iprint

    integer:: ierr
    
    call read_frcmod_params(myid,mpi_world,iprint)

    call mpi_bcast(ctype_fmod,128,mpi_character,0,mpi_world,ierr)
    call mpi_bcast(is_fmod,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(alpha,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(rin,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(rout,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(dth_nn,1,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(nn_intact,1,mpi_integer,mpi_world,ierr)

!.....Error handling
    if( trim(ctype_fmod).eq.'species' .and. is_fmod.lt.1 ) then
      print *,'ERROR: species is not specified in frcmod, is_fmod=',is_fmod
      stop
    else if( trim(ctype_fmod).eq.'num_neighbors' .and. &
         (rin.lt.0d0 .or. rout.lt.0d0 .or. rout.lt.rin) ) then
      print *,'ERROR: rin and/or rout are wrong, rin,rout=',rin,rout
      stop
    endif

    if( trim(ctype_fmod).eq.'num_neighbors' ) then
      rout2 = rout*rout
    endif

    if( myid.eq.0 ) then
      print *,'FRCMOD parameters:'
      print '(a,a)','  type:    ',trim(ctype_fmod)
      print '(a,f6.3)','  alpha:   ',alpha
      if( trim(ctype_fmod).eq.'num_neighbors' ) then
        print '(a,2f7.3)','  rin,rout:   ',rin,rout
        print '(a,i3)','  nn_intact:   ',nn_intact
        print '(a,f6.2)','  dth_nn:   ',dth_nn
      else  ! default: species
        print '(a,i3)','  species: ',is_fmod
      endif
    endif
    return
  end subroutine init_frcmod
!=======================================================================
  subroutine read_frcmod_params(myid,mpi_world,iprint)
    integer,intent(in):: myid,mpi_world,iprint
    
    integer,external:: num_data
    character(len=128):: cfname,c1st
    
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
          read(ioprms,*) c1st, is_fmod
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
        else if( trim(c1st).eq.'nn_intact' ) then
          backspace(ioprms)
          read(ioprms,*) c1st, nn_intact
        else
          print *,'No such FRCMOD keyword: ',trim(c1st)
        endif
      enddo
      
10    close(ioprms)
    endif ! myid.eq.0
    return
  end subroutine read_frcmod_params
!=======================================================================
  subroutine modify_forces(namax,natm,aa,tag,ra,h,nnmax,lspr)
!
!  Modify forces on atoms that are selected in some way...
!
    integer,intent(in):: namax,natm,nnmax,lspr(0:nnmax,namax)
    real(8),intent(in):: ra(3,namax),tag(namax),h(3,3)
    real(8),intent(inout):: aa(3,namax)

    integer:: ia,is
    real(8):: beta

    if( trim(ctype_fmod).eq.'num_neighbors' ) then
      if( .not. allocated(ann) ) allocate(ann(namax))
      call count_nn(namax,natm,ra,h,nnmax,lspr)
      do ia=1,natm
!!$        beta = 1d0
!!$        if( ann(ia).lt.nn_intact -(dth_nn+0.1d0) ) then
!!$          beta = alpha
!!$        else if( ann(ia).lt.nn_intact -(dth_nn-0.1d0) ) then
!!$          beta = alpha +(1d0-alpha)*0.5d0*(1d0 &
!!$               -cos(pi *(ann(ia)-(nn_intact-(dth_nn-0.1d0))) /(2*0.1d0)))
!!$        endif
        beta = 1d0 -(1d0-alpha)*exp(-(ann(ia)-7.0)**2/dth_nn**2/2)
!!$        print '(a,i6,3f8.3)','ia,ann,beta,exp=',ia,ann(ia),beta&
!!$             ,exp(-(ann(ia)-7.0)**2/dth_nn**2/2)
        aa(1:3,ia) = aa(1:3,ia) *beta
!!$        aa(1:3,ia) = aa(1:3,ia)*(alpha  &
!!$             +(1d0-alpha)*exp(-(ann(ia)-nn_intact)**2/dth_nn**2/2))
!!$        if( abs(ann(ia) -nn_intact).gt.dth_nn ) then
!!$          aa(1:3,ia) = alpha *aa(1:3,ia)
!!$        endif
      enddo
    else ! default: species
      do ia=1,natm
        is = int(tag(ia))
        if( is.eq.is_fmod ) aa(1:3,ia) = alpha *aa(1:3,ia)
      enddo
    endif
    return
  end subroutine modify_forces
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
end module frcmod
