module time
!-----------------------------------------------------------------------
!                     Last modified: <2021-07-01 14:55:28 Ryo KOBAYASHI>
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
  subroutine report_time(ionum,iprint)
!
! Write a report about measured times.
!
    include './const.h'
    integer,intent(in):: ionum,iprint
    integer:: i,maxlen,itot,nspace,ihour,imin,isec
    character(len=3):: clen,cspace
    character(len=128):: ctn
    real(8):: time_tot

!.....Total time
    time_tot = -1d0
    do i=1,ntimes
      if( trim(ctnames(i)).eq.'total' ) then
        time_tot = etimes(i)
        itot = i
        exit
      endif
    enddo

    if( iprint.ge.ipl_time ) then
!.....Max length of name for time
      maxlen = 0
      do i=1,ntimes
        maxlen = max(maxlen,len_trim(ctnames(i)))
      enddo
      write(clen,'(i0)') maxlen

      write(ionum,*) ''
      write(ionum,'(a)') ' Report from time module:'
      nspace = 3
      write(cspace,'(i0)') nspace
      do i=1,ntimes
        ctn = trim(ctnames(i))
        if( trim(ctn).ne.'' .and. trim(ctn).ne.'total' ) then
          if( time_tot.ge.0d0 ) then
            write(ionum,'('//trim(cspace)//'x,a,a'//trim(clen)//',a,f10.3,a,f7.1,a,i0)') &
                 'Time ',trim(ctnames(i)),' = ', etimes(i), &
                 ' sec, ',etimes(i)/time_tot*100, &
                 ' %, # of calls = ',ncalls(i)
          else
            write(ionum,'('//trim(cspace)//'x,a,a'//trim(clen)//',a,f10.3,a,i0)') &
                 'Time ',trim(ctnames(i)),' = ', etimes(i), &
                 ' sec, # of calls = ',ncalls(i)
          endif
        endif
      enddo
      
    else  ! iprint.lt.ipl_time
      nspace = 1
      write(cspace,'(i0)') nspace
      write(clen,'(i0)') 7
    endif
    call sec2hms(time_tot,ihour,imin,isec)
    write(ionum,'('//trim(cspace)//'x,a,a'//trim(clen)//',a,f10.3,a,i3,"h",i2.2,"m",i2.2,"s")') &
         'Time ','total',' = ', time_tot, ' sec = ',ihour,imin,isec
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

