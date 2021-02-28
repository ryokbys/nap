module composition
  use pmdvars,only: nspmax
  implicit none
  save

!.....Composition information should be ready if composition is used.
!  File format is like the following:
!  
!  $ cat in.data.compos
!  #  ID, c(1), c(2), ..., c(nsp), Emin per atom
!     30      <=== num of compositions
!     1    0.1   0.2  ...   0.0    -5.423
!     2    0.3   0.2  ...   0.0    -4.239
!     ...
!     -3.453  <=== at the end, Emin per atom of the other compositions
!
  character(len=128):: composfname = 'in.data.compos'

  integer:: io_compos = 80

  logical:: lwgt_compos = .false. 

  integer,parameter:: max_compos = 128
  integer:: num_compos = 0
  real(8):: data_compos(nspmax,0:max_compos)
  real(8):: emin_compos(0:max_compos)

  real(8),parameter:: eps = 1.0d-3

!.....Scale of weight in eV, (wgt = exp(-dE/escl_compos))
  real(8):: escl_compos = 1.0d0

contains
!=======================================================================
  subroutine read_compos_list()
!
!  Read the list of compositions of all the samples from the file.
!  The composition data should be ready beforehand.
!
    use variables
    use parallel
    
    integer:: ismpl,natms(nspmax),ic,i,itmp
    real(8):: dcompos(nspmax)
    character:: ctmp*128,cnum*3
    

!.....Read data only at the root node
    if( myid.eq.0 ) then
      open(io_compos,file=trim(composfname),status='old')
      read(io_compos,*) ctmp
      read(io_compos,*) num_compos
      do ic=1,num_compos
        read(io_compos,*) itmp,data_compos(1:nsp,ic), emin_compos(ic)
      enddo
      read(io_compos,*) emin_compos(0)
      close(io_compos)
    endif

!.....Broadcast
    call mpi_bcast(num_compos,1,mpi_integer,0,mpi_world,ierr)
    call mpi_bcast(data_compos,(max_compos+1)*nspmax,mpi_real8,0,mpi_world,ierr)
    call mpi_bcast(emin_compos,(max_compos+1),mpi_real8,0,mpi_world,ierr)

    if( myid.eq.0 .and. iprint.gt.1 ) then
      print *,''
      print *,'Composition list:'
      write(cnum,'(i0)') nsp
      do ic=1,num_compos
        print '(i5,'//trim(cnum)//'(2x,f5.3),es11.3)', &
             ic,(data_compos(i,ic),i=1,nsp),emin_compos(ic)
      enddo
      print '(a,es11.3)', '   others     ',emin_compos(0)
    endif

    return
  end subroutine read_compos_list
!=======================================================================
  subroutine assign_compos_weight()
!
!  Assign composition ID to all the samples.
!
    use variables
    use parallel

    integer:: ismpl,natm,natms(nspmax),ia,is,ic
    real(8):: dcompos(nspmax),wgt,de,eref
    logical:: lsame
    type(mdsys):: smpl

    do ismpl=isid0,isid1
      smpl = samples(ismpl)
      natm = smpl%natm
      eref = smpl%eref
      natms(:) = 0
      do ia=1,natm
        is = int(smpl%tag(ia))
        natms(is) = natms(is) +1
      enddo
      dcompos(:) = dble(natms(:)) /natm
      do ic=1,num_compos
        lsame = .true.
        do is=1,nsp
          if( abs(data_compos(is,ic)-dcompos(is)).gt.eps ) then
            lsame = .false.
            exit
          endif
        enddo
        if( lsame ) then
          de = max(eref/natm -emin_compos(ic), 0d0)
          wgt = exp(-de /escl_compos)
!!$          print '(a,2i6,2es12.3)','  ismpl,ic,de,wgt=',ismpl,ic,de,wgt
          samples(ismpl)%wgt = wgt
          exit
        endif
      enddo
    enddo
    
  end subroutine assign_compos_weight
end module composition
!-----------------------------------------------------------------------
! Local Variables:
! compile-command: "make fitpot"
! End:
