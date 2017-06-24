module fp_common
!
! Module that contains common functions/subroutines for fitpot.
!
  implicit none

  integer:: maxna = 0
  real(8),allocatable:: fdiff(:,:)


contains
!=======================================================================
  subroutine func_w_pmd(ndim,x,ftrn,ftst)
!
!  Evaluate loss function value using pmd (actually one_shot routine.)
!
    use variables,only:nsmpl,nsmpl_trn,samples,nprcs,tfunc &
         ,lfmatch,nfunc,tcomm,mdsys,erefmin &
         ,cmaindir,cevaltype,swgt2trn,swgt2tst,cpot
    use parallel
    use minimize
    use Coulomb,only: set_paramsdir_Coulomb
    use Morse,only: set_paramsdir_Morse
    implicit none
    
    integer,intent(in):: ndim
    real(8),intent(in):: x(ndim)
    real(8),intent(out):: ftrn,ftst

    integer:: ismpl,natm,ia,ixyz,idim
    real(8):: dn3i,ediff,fscale,eref,swgt,wgtidv
    real(8):: eerr,ferr,ferri
    real(8):: ftrnl,ftstl,ftmp
    real(8):: edenom,fdenom
    real(8):: tfl,tcl,tfg,tcg,tf0,tc0
    type(mdsys):: smpl
    logical:: l1st = .true.
    logical:: lcalcgrad = .false.
    real(8),save,allocatable:: gdummy(:)
    character(len=128):: cdirname

    nfunc= nfunc +1

    tc0= mpi_wtime()
    call mpi_bcast(x,ndim,mpi_real8,0,mpi_world,ierr)
    tcl= mpi_wtime() -tc0
    tf0= mpi_wtime()

    if( l1st ) then
      if( .not.allocated(gdummy) ) allocate(gdummy(ndim))
    endif
    
    ftrnl = 0d0
    ftstl = 0d0
    do ismpl=isid0,isid1
      smpl = samples(ismpl)
      cdirname = smpl%cdirname
      if( trim(cpot).eq.'vcMorse' ) then
        call set_paramsdir_Morse(cdirname)
        call set_paramsdir_Coulomb(cdirname)
      endif
      call run_pmd(smpl,lcalcgrad,ndim,gdummy)

      ftmp = 0d0
      natm= smpl%natm
      if( natm.gt.maxna ) then
        maxna = natm
        if( allocated(fdiff) ) deallocate(fdiff)
        allocate(fdiff(3,maxna))
      endif
      eref= smpl%eref
      eerr = smpl%eerr
      swgt = smpl%wgt
!.....Energy matching
      ediff= (smpl%epot -eref)/natm /eerr
      ediff= ediff*ediff
      ftmp= ftmp +ediff *swgt
!.....Force matching
      if( lfmatch .and. smpl%nfcal.ne.0 ) then
        ferr = smpl%ferr
        ferri = 1d0/ferr
        dn3i = 1d0/3/smpl%nfcal
        do ia=1,natm
          if( smpl%ifcal(ia).eq.0 ) cycle
          do ixyz=1,3
            fdiff(ixyz,ia)= (smpl%fa(ixyz,ia) &
                 -smpl%fref(ixyz,ia)) *ferri
            fdiff(ixyz,ia)= fdiff(ixyz,ia)*fdiff(ixyz,ia)
            ftmp= ftmp +fdiff(ixyz,ia) *dn3i *swgt
          enddo
        enddo
      endif
      if( smpl%iclass.eq.1 ) then
        ftrnl = ftrnl +ftmp
      else if( smpl%iclass.eq.2 ) then
        ftstl = ftstl +ftmp
      endif
    enddo  ! ismpl

!    tfunc= tfunc +mpi_wtime() -tf0
    tfl = mpi_wtime() -tf0

    tc0= mpi_wtime()
    ftrn= 0d0
    ftst = 0d0
    call mpi_allreduce(ftrnl,ftrn,1,mpi_real8,mpi_sum,mpi_world,ierr)
    call mpi_allreduce(ftstl,ftst,1,mpi_real8,mpi_sum,mpi_world,ierr)
    ftrn = ftrn /swgt2trn
    ftst = ftst /swgt2tst
    tcl = tcl + (mpi_wtime() -tc0)

!.....only the bottle-neck times are taken into account
    call mpi_reduce(tcl,tcg,1,mpi_real8,mpi_max,0,mpi_world,ierr)
    call mpi_reduce(tfl,tfg,1,mpi_real8,mpi_max,0,mpi_world,ierr)
    tcomm= tcomm +tcg
    tfunc= tfunc +tfg
    l1st = .false.

  end subroutine func_w_pmd
!=======================================================================
  subroutine grad_w_pmd(ndim,x,gtrn)
!
!  Evaluate the gradient of loss function value
!  using pmd (actually one_shot routine.)
!
    use variables,only: nsmpl,nsmpl_trn,tgrad,ngrad,tcomm &
         ,samples,mdsys,swgt2trn,swgt2tst,cpot
    use parallel
    use minimize
    use Coulomb,only: set_paramsdir_Coulomb
    use Morse,only: set_paramsdir_Morse
    implicit none
    integer,intent(in):: ndim
    real(8),intent(in):: x(ndim)
    real(8),intent(out):: gtrn(ndim)
    
    integer:: ismpl,i,idim,natm
    real(8),save,allocatable:: gs(:),gtrnl(:)
    real(8):: tcl,tgl,tcg,tgg,tc0,tg0
    real(8):: ediff,eerr,eref,swgt
    type(mdsys):: smpl
    logical:: lcalcgrad = .true. 
    character(len=128):: cdirname

    if( .not.allocated(gs) ) allocate(gs(ndim),gtrnl(ndim))

    ngrad= ngrad +1
    tg0= mpi_wtime()

    gtrnl(1:ndim) = 0d0

    do ismpl=isid0,isid1
      smpl = samples(ismpl)
      cdirname = smpl%cdirname
!.....Since g calc is time consuming,
!.....not calculate g for test set.
      if( smpl%iclass.ne.1 ) cycle
      if( trim(cpot).eq.'vcMorse' ) then
        call set_paramsdir_Morse(cdirname)
        call set_paramsdir_Coulomb(cdirname)
      endif
      call run_pmd(smpl,lcalcgrad,ndim,gs)

      natm= smpl%natm
      eref= smpl%eref
      eerr= smpl%eerr
      swgt= smpl%wgt
!.....Energy matching
      ediff= (smpl%epot -eref) /natm /eerr
      ediff= 2d0 *ediff /natm /eerr *swgt
      gtrnl(1:ndim)= gtrnl(1:ndim) +gs(1:ndim)*ediff
!.....TODO: force matching
      
    enddo  ! ismpl

    tgl= mpi_wtime() -tg0

    tc0= mpi_wtime()
    gtrn(1:ndim) = 0d0
    call mpi_allreduce(gtrnl,gtrn,ndim,mpi_real8,mpi_sum,mpi_world,ierr)
    tcl= mpi_wtime() -tc0

    gtrn(1:ndim)= gtrn(1:ndim) /swgt2trn

!.....only the bottle-neck times are taken into account
    call mpi_reduce(tcl,tcg,1,mpi_real8,mpi_max,0,mpi_world,ierr)
    call mpi_reduce(tgl,tgg,1,mpi_real8,mpi_max,0,mpi_world,ierr)
    tcomm= tcomm +tcg
    tgrad= tgrad +tgg
    
    return
    
  end subroutine grad_w_pmd
!=======================================================================
  function ndat_in_line(ionum,delim) result(ndat)
    implicit none
    integer,intent(in):: ionum
    character(len=1),intent(in):: delim
    integer:: ndat
    integer,external:: num_data
    character:: ctmp*128

    read(ionum,'(a)') ctmp
    ndat = num_data(trim(ctmp),delim)
    backspace(ionum)
    return

  end function ndat_in_line
  
end module fp_common
