module time
!-----------------------------------------------------------------------
!                     Last modified: <2021-02-06 17:45:15 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
! Module for time measurement.
!-----------------------------------------------------------------------
  implicit none
  save

!.....Maximum number of time measurements
  integer,parameter:: mtimes = 128
  integer:: ntimes = 0
  integer:: ncalls(mtimes)
  real(8):: etimes(mtimes)
  character(len=32):: ctnames(mtimes)

contains
!=======================================================================
  subroutine time_stamp(prefix)
    implicit none
    character(len=*),intent(in):: prefix
    character(len=10):: c1,c2,c3
    integer:: date(8)
    character(len=128):: cdate,ctime

    call date_and_time(c1,c2,c3,date)

    write(ctime,'(i0.2,a,i0.2,a,i0.2)') date(5),':',date(6) &
         ,':',date(7)
    write(cdate,'(i4,a,i0.2,a,i0.2)') date(1),'-',date(2),'-',date(3)
    write(6,'(a,1x,a,1x,a)') prefix, 'at '//trim(ctime), &
         'on '//trim(cdate)
    return
  end subroutine time_stamp
!=======================================================================
  function get_time_id(cname)
    character(len=*),intent(in):: cname
    integer:: i
    integer:: get_time_id

    get_time_id = 0

!.....Search for the name in the list
    do i=1,ntimes
      if( trim(ctnames(i)).eq.trim(cname) ) then
        get_time_id = i
        return
      endif
    enddo

!.....If not match, register the name to the list, and increase ntimes
    if( get_time_id.eq.0 ) then
      ntimes = ntimes +1
      ctnames(ntimes) = trim(cname)
      etimes(ntimes) = 0d0
      ncalls(ntimes) = 0
      get_time_id = ntimes
    endif
    return
    
  end function get_time_id
!=======================================================================
  subroutine accum_time(cname,etime)
    character(len=*),intent(in):: cname
    real(8),intent(in):: etime
    integer:: idtime

    idtime = get_time_id(cname)
    etimes(idtime) = etimes(idtime) +etime
    ncalls(idtime) = ncalls(idtime) +1
    return
  end subroutine accum_time
!=======================================================================
  subroutine report_time(ionum)
!
! Write a report about measured times.
!
    integer,intent(in):: ionum
    integer:: i,maxlen
    character(len=3):: clen
    real(8):: time_tot

!.....Max length of name for time
    maxlen = 0
    do i=1,ntimes
      maxlen = max(maxlen,len_trim(ctnames(i)))
    enddo
    write(clen,'(i0)') maxlen

!.....Total time
    time_tot = -1d0
    do i=1,ntimes
      if( trim(ctnames(i)).eq.'total' ) then
        time_tot = etimes(i)
        exit
      endif
    enddo

    write(ionum,*) ''
    write(ionum,'(a)') ' Report from time module:'
    do i=1,ntimes
      if( trim(ctnames(i)).ne.'' ) then
        if( time_tot.ge.0d0 ) then
          write(ionum,'(a,a'//trim(clen)//',a,f10.2,a,f7.1,a,i0)') &
               '   Time ',trim(ctnames(i)),' = ', etimes(i), &
               ' sec, ',etimes(i)/time_tot*100, &
               ' %, # of calls = ',ncalls(i)
        else
          write(ionum,'(a,a'//trim(clen)//',a,f10.2,a,i0)') &
               '   Time ',trim(ctnames(i)),' = ', etimes(i), &
               ' sec, # of calls = ',ncalls(i)
        endif
      endif
    enddo
  end subroutine report_time
!=======================================================================
  subroutine sec2hms(sec,h,m,s)
!
!  Convert second to hours,minutes,seconds
!
    real(8),intent(in):: sec
    integer,intent(out):: h,m,s
    h = int(sec/3600)
    m = int((sec-h*3600)/60)
    s = int(sec -h*3600 -m*60)
    return
  end subroutine sec2hms
end module time
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:

