subroutine read_vars()
  use variables
  use parallel
  use random
  implicit none
  integer:: i
  real(8):: rs0
  character(len=128):: fname

  if( index(cpot,'NN').ne.0 ) then
    call read_vars_NN()
  else
    call read_vars_fitpot()
  endif

end subroutine read_vars
!=======================================================================
subroutine write_vars(cadd)
  use variables
  use parallel
!!$  use NNd, only: NN_standardize, NN_restore_standard
  use fp_common,only: normalize, restore_normalize
  implicit none
  character(len=*),intent(in):: cadd
  integer:: i
  character(len=128):: cfname

  if( cnormalize(1:4).ne.'none' ) then
!!$    if( trim(cpot).eq.'NN' .and. .not. &
!!$         (trim(cfmethod).eq.'sa' .or. trim(cfmethod).eq.'SA' .or. &
!!$         trim(cfmethod).eq.'ga' .or. trim(cfmethod).eq.'GA' .or. &
!!$         trim(cfmethod).eq.'de' .or. trim(cfmethod).eq.'DE' .or. &
!!$         trim(cfmethod).eq.'pso' .or. trim(cfmethod).eq.'PSO') ) then
!!$      call NN_restore_standard()
!!$    else if( lnormalize ) then
    if( lnormalize ) then
      call restore_normalize()
     endif
  endif

!!$  cfname= trim(cmaindir)//'/'//trim(cparfile)//'.'//trim(cadd)
  cfname= trim(cparfile)//'.'//trim(cadd)

  if( index(cpot,'NN').ne.0 ) then
    call write_vars_NN(cfname)
  else
    call write_vars_fitpot(cfname)
  endif

  if( cnormalize(1:4).ne.'none' ) then
!!$    if( trim(cpot).eq.'NN' .and. .not. &
!!$         (trim(cfmethod).eq.'sa' .or. trim(cfmethod).eq.'SA' .or. &
!!$         trim(cfmethod).eq.'ga' .or. trim(cfmethod).eq.'GA' .or. &
!!$         trim(cfmethod).eq.'de' .or. trim(cfmethod).eq.'DE' .or. &
!!$         trim(cfmethod).eq.'pso' .or. trim(cfmethod).eq.'PSO') ) then
!!$      call NN_standardize()
!!$    else if( lnormalize ) then
    if( lnormalize ) then
      call normalize()
    endif
  endif

end subroutine write_vars
!=======================================================================
subroutine read_vars_fitpot()
!
!  Read fitpot original param-format file.
!  Linreg potenitla param file is also the same style.
!
  use variables
  use parallel
  use random
  implicit none
  integer,parameter:: ionum = 15
  integer:: i
  real(8):: rs0
  character(len=128):: fname

  if( myid.eq.0 ) then
    print *,'Read parameters to be optimized from '//trim(cparfile)
    open(ionum,file=trim(cparfile),status='old')
    read(ionum,*) nvars, rcut, rc3
    rcut = max(rcut,rc3)
  endif
  call mpi_bcast(nvars,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(rcut,1,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(rc3,1,mpi_real8,0,mpi_world,ierr)
  allocate(vars(nvars),vranges(2,nvars))
  if( myid.eq.0 ) then
    print '(a,i0)',' Number of variables to be optimized = ',nvars
    do i=1,nvars
      read(ionum,*) vars(i),vranges(1:2,i)
    enddo
    if( index(cinitv,'gauss').ne.0 ) then
      rs0 = get_seed()
      call set_seed(vinitrs)
      do i=1,nvars
        vars(i) = vinitsgm*(polarbm()-vinitmu)
      enddo
      call set_seed(rs0)
      write(6,'(a)') ' Potential parameters are shuffled'&
           //' to give normal distribution'
      write(6,'(a,2es10.2)') '   with mu and sgm =',vinitmu,vinitsgm
    else
      write(6,'(a)') ' Potential parameters are read from file: '//trim(cparfile)
    endif
    close(ionum)
  endif
  call mpi_bcast(vars,nvars,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(vranges,2*nvars,mpi_real8,0,mpi_world,ierr)

end subroutine read_vars_fitpot
!=======================================================================
subroutine write_vars_fitpot(cfname)
!
!  Write params in original fitpot format.
!
  use variables
  use parallel
  implicit none
  character(len=128),intent(in):: cfname
  integer,parameter:: ionum = 16
  integer:: i

  if( myid.eq.0 ) then
    open(ionum,file=trim(cfname),status='replace')
    write(ionum,'(i10,2es15.4)') nvars,rcut,rc3
    do i=1,nvars
      write(ionum,'(es23.14e3,2es12.4)') vars(i),vranges(1:2,i)
    enddo
    close(ionum)
!    print *, 'wrote '//trim(cfname)
  endif

end subroutine write_vars_fitpot
!=======================================================================
subroutine read_vars_NN()
!
!  Read in.params.NN2/DNN and convert data to vars array.
!
  use variables
  use parallel
  use random
  implicit none

  integer,parameter:: ionum = 15
  integer:: i
  real(8):: rs0
  character(len=128):: fname,ctmp

  if( myid.eq.0 ) then
    print *,'Read parameters to be optimized from '//trim(cparfile)
    open(ionum,file=trim(cparfile),status='old')
!.....num of symmetry functions, num of node in 1st hidden layer
10  read(ionum,'(a)') ctmp
    if( ctmp(1:1).eq.'!' .or. ctmp(1:1).eq.'#' ) then
      call parse_option_NN(ctmp)
      goto 10
    else
      read(ctmp,*) nn_nl
      backspace(ionum)
    endif
    allocate(nn_nhl(0:nn_nl+1))
    read(ionum,*) nn_nl, (nn_nhl(i),i=0,nn_nl)
    nn_nhl(nn_nl+1)= 1

!.....Bias treatment depends on NN version
    nvars= 0
    if( trim(cpot).eq.'NN2' ) then
      do i=1,nn_nl+1
        nvars= nvars +nn_nhl(i-1)*nn_nhl(i)
      enddo
    else if( trim(cpot).eq.'DNN' ) then
      do i=1,nn_nl+1
        nvars = nvars +(nn_nhl(i-1)+1)*nn_nhl(i)
      enddo
    endif
!!$!.....
!!$    read(ionum,*) nvars, rcut, rc3
  endif
  call mpi_bcast(nn_nl,1,mpi_integer,0,mpi_world,ierr)
  if( .not. allocated(nn_nhl) ) allocate(nn_nhl(0:nn_nl+1))
  call mpi_bcast(nn_nhl,nn_nl+2,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(nn_sigtype,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(nvars,1,mpi_integer,0,mpi_world,ierr)
  allocate(vars(nvars),vranges(2,nvars))
  if( myid.eq.0 ) then
    print '(a,i0)',' Number of variables to be optimized = ',nvars
    do i=1,nvars
      read(ionum,*) vars(i),vranges(1:2,i)
    enddo
    if( index(cinitv,'gauss').ne.0 ) then
      rs0 = get_seed()
      call set_seed(vinitrs)
      do i=1,nvars
        vars(i) = vinitsgm*(polarbm()-vinitmu)
      enddo
      call set_seed(rs0)
      write(6,'(a)') ' Potential parameters are shuffled'&
           //' to give normal distribution'
      write(6,'(a,2es10.2)') '   with mu and sgm =',vinitmu,vinitsgm
    else
      write(6,'(a)') ' Potential parameters are read from file: '//trim(cparfile)
    endif
    close(ionum)
  endif
  call mpi_bcast(vars,nvars,mpi_real8,0,mpi_world,ierr)
  call mpi_bcast(vranges,2*nvars,mpi_real8,0,mpi_world,ierr)

end subroutine read_vars_NN
!=======================================================================
subroutine write_vars_NN(cfname)
!
!  Write parameters in NN2/DNN format
!
  use variables
  use parallel
  implicit none
  character(len=*),intent(in):: cfname
  integer:: i
  integer,parameter:: ionum = 16
  
  if( myid.eq.0 ) then
    open(ionum,file=trim(cfname),status='replace')
    write(ionum,'(a,i2)') '!  sigtype: ',nn_sigtype
    write(ionum,'(a)') '!  '
    write(ionum,'(100i4)') nn_nl, (nn_nhl(i),i=0,nn_nl)
!!$    write(ionum,'(i10,2es15.4)') nvars,rcut,rc3
    do i=1,nvars
      write(ionum,'(es23.14e3,2es12.4)') vars(i),vranges(1:2,i)
    enddo
    close(ionum)
!    print *, 'wrote '//trim(cfname)
  endif

end subroutine write_vars_NN
!=======================================================================
subroutine parse_option_NN(cline)
!
!  Parse options from a comment line.
!  Lines starting from ! or # are treated as comment lines,
!  but options can be given at the comment lines.
!  The option words should be put after these comment characters with
!  one or more spaces between them for example,
!
!  bias:  .true.
!
!  Currently available options are:
!    - "sigtype:" sigmoid type: 1 or 2.
!
  use variables
  use parallel
  implicit none
  character(len=*),intent(in):: cline

  character(len=10):: c1,copt
  logical:: lopt
  integer:: iopt

  ierr = 0
  if( index(cline,'sigtype:').ne.0 ) then
    read(cline,*) c1,copt,iopt
    nn_sigtype = iopt
  endif

end subroutine parse_option_NN
!-----------------------------------------------------------------------
! Local Variables:
! compile-command: "make fitpot"
! End:
