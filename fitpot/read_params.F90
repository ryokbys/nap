subroutine read_params_fitpot()
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
    if( trim(cinitv).eq.'gaussian' .or. trim(cinitv).eq.'gauss' ) then
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

end subroutine read_params_fitpot
!=======================================================================
subroutine write_params_fitpot(cfname)
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

end subroutine write_params_fitpot
!=======================================================================
subroutine read_params_NN2()
!
!  Read in.params.NN2 and convert data to vars array.
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
      call parse_option_NN2(ctmp)
      goto 10
    else
      backspace(ionum)
    endif
    read(ionum,*) nn_nl, (nn_nhl(i),i=0,nn_nl)
    nn_nhl(nn_nl+1)= 1
    nvars= 0
    do i=1,nn_nl+1
      nvars= nvars +nn_nhl(i-1)*nn_nhl(i)
    enddo
!!$!.....
!!$    read(ionum,*) nvars, rcut, rc3
  endif
  call mpi_bcast(nn_nl,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(nn_nhl,nn_nlmax+1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(nn_sigtype,1,mpi_integer,0,mpi_world,ierr)
  call mpi_bcast(nvars,1,mpi_integer,0,mpi_world,ierr)
!!$  call mpi_bcast(rcut,1,mpi_real8,0,mpi_world,ierr)
!!$  call mpi_bcast(rc3,1,mpi_real8,0,mpi_world,ierr)
  allocate(vars(nvars),vranges(2,nvars))
  if( myid.eq.0 ) then
    print '(a,i0)',' Number of variables to be optimized = ',nvars
    do i=1,nvars
      read(ionum,*) vars(i),vranges(1:2,i)
    enddo
    if( trim(cinitv).eq.'gaussian' .or. trim(cinitv).eq.'gauss' ) then
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

end subroutine read_params_NN2
!=======================================================================
subroutine write_params_NN2(cfname)
!
!  Write parameters in NN2 format
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
    write(ionum,'(i4i5)') nn_nl, (nn_nhl(i),i=0,nn_nl)
!!$    write(ionum,'(i10,2es15.4)') nvars,rcut,rc3
    do i=1,nvars
      write(ionum,'(es23.14e3,2es12.4)') vars(i),vranges(1:2,i)
    enddo
    close(ionum)
!    print *, 'wrote '//trim(cfname)
  endif

end subroutine write_params_NN2
!=======================================================================
subroutine parse_option_NN2(cline)
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

end subroutine parse_option_NN2
!=======================================================================
!-----------------------------------------------------------------------
! Local Variables:
! compile-command: "make fitpot"
! End:
