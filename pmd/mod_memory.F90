module memory
!-----------------------------------------------------------------------
!                     Last modified: <2021-11-04 11:41:47 Ryo KOBAYASHI>
!-----------------------------------------------------------------------
! Module for memory measurement.
!-----------------------------------------------------------------------
  implicit none
  save

!.....Maximum number of memory measurements
  integer,parameter:: mmems = 128
  integer:: nmems = 0
  integer:: mems(mmems)   ! byte
  character(len=32):: cmnames(mmems)

contains
!=======================================================================
  function get_mem_id(cname)
    character(len=*),intent(in):: cname
    integer:: i
    integer:: get_mem_id

    get_mem_id = 0

!.....Search for the name in the list
    do i=1,nmems
      if( trim(cmnames(i)).eq.trim(cname) ) then
        get_mem_id = i
        return
      endif
    enddo

!.....If not match, register the name to the list, and increase nmems
    if( get_mem_id.eq.0 ) then
      nmems = nmems +1
      cmnames(nmems) = trim(cname)
      mems(nmems) = 0
      get_mem_id = nmems
    endif
    return
    
  end function get_mem_id
!=======================================================================
  subroutine accum_mem(cname,memin)
    character(len=*),intent(in):: cname
    integer,intent(in):: memin
    integer:: idmem

    logical,save:: l1st = .true.

    if( l1st ) then
      mems(:) = 0
      l1st = .false.
    endif

    idmem = get_mem_id(cname)
    mems(idmem) = mems(idmem) +memin
    return
  end subroutine accum_mem
!=======================================================================
  function get_mem_total()
    integer:: get_mem_total
    integer:: i
    
    get_mem_total = 0
    do i=1,nmems
      get_mem_total = get_mem_total +mems(i)
    enddo
    return
  end function get_mem_total
!=======================================================================
  subroutine report_mem(ionum,iprint)
!
! Write a report about measured memory.
!
    include './const.h'
    integer,intent(in):: ionum,iprint
    integer:: i,maxlen
    character(len=3):: clen
    integer:: mem_tot
    character(len=4):: cunit
    real(8):: rmem

    mem_tot = get_mem_total()

    if( iprint.ge.ipl_time ) then
!.....Max length of name for mem
      maxlen = 0
      do i=1,nmems
        maxlen = max(maxlen,len_trim(cmnames(i)))
      enddo
      write(clen,'(i0)') maxlen

      write(ionum,*) ''
      write(ionum,'(a)') ' Report from memory module:'
      do i=1,nmems
        if( trim(cmnames(i)).ne.'' ) then
          call mem_in_unit(mems(i),rmem,cunit)
          write(ionum,'(a,a'//trim(clen)//',a,f10.3,a)') &
               '   Memory ',trim(cmnames(i)),' = ', rmem,' '//cunit
        endif
      enddo
      call mem_in_unit(mem_tot,rmem,cunit)
      write(ionum,'(a,a'//trim(clen)//',a,f10.3,a)') &
           '   Memory ','total',' = ', rmem,' '//cunit
    else
      call mem_in_unit(mem_tot,rmem,cunit)
      write(ionum,'(a,f10.3,a)') &
           ' Memory per MPI-proc = ', rmem,' '//cunit
    endif
  end subroutine report_mem
!=======================================================================
  subroutine mem_in_unit(inmem,outmem,cout)
!
!  Convert the unit of memory.
!  Original unit is Byte, and returns the memory in Byte, kB, MB, or GB.
!
    integer,intent(in):: inmem
    real(8),intent(out):: outmem
    character(len=4):: cout

    if( inmem.lt.1000 ) then
      outmem = dble(inmem)
      cout = 'Byte'
    else if( inmem.lt.1000000 ) then
      outmem = dble(inmem)/1000
      cout = 'kB'
    else if( inmem.lt.1000000000 ) then
      outmem = dble(inmem)/1000000
      cout = 'MB'
    else
      outmem = dble(inmem)/1000000000
      cout = 'GB'
    endif
    
    return
  end subroutine mem_in_unit
end module memory
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd"
!     End:

