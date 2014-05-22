      subroutine write_akr(ionum,cfname
     &     ,natm,h,hunit,tag,ra,va)
      implicit none
      include './params_unit.h'
      integer,intent(in):: ionum
      character(len=*),intent(in) :: cfname
      integer,intent(in):: natm
      real(8),intent(in):: h(3,3,0:1),tag(natm),ra(3,natm),va(3,natm)
     &     ,eki(3,3,natm),epi(natm),hunit

      integer:: ia,ib,l,i

      open(ionum,file=cfname,status='replace')
      write(ionum,'(es15.7)') hunit
      write(ionum,'(3es15.7)') ((h(ia,ib,l)/hunit,ia=1,3),ib=1,3)
      write(ionum,'(i10,3i5)') natm,3,0,0
      write(ionum,'(7es23.14e3,11es22.14)') (int(tag(i)),ra(1:3,i)
     &     ,va(1:3,i),i=1,natm)
      close(ionum)

      end subroutine write_akr
