      subroutine write_akr(ionum,cfname
     &     ,natm,h,hunit,tag,ra,va)
      implicit none
      include '../pmd/params_unit.h'
      integer,intent(in):: ionum
      character(len=*),intent(in) :: cfname
      integer,intent(in):: natm
      real(8),intent(in):: h(3,3,0:1),tag(natm),ra(3,natm),va(3,natm)
     &     ,hunit

      integer:: ia,ib,l,i

      open(ionum,file=cfname,status='replace')
      write(ionum,'(es15.7)') hunit
      write(ionum,'(3es15.7)') ((h(ia,ib,0)/hunit,ia=1,3),ib=1,3)
      write(ionum,'(i10,3i5)') natm,3,0,0
      do i=1,natm
        write(ionum,'(i4,3f11.7,11f7.3)') int(tag(i)),ra(1:3,i)
     &       ,va(1:3,i)
      enddo
      close(ionum)

      end subroutine write_akr
