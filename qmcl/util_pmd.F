      subroutine write_pmd_ascii(ionum,cfname
     &     ,natm,h,hunit,tag,ra,va,eki,epi,strs)
      implicit none
      integer,intent(in):: ionum
      character(len=*),intent(in) :: cfname
      integer,intent(in):: natm
      real(8),intent(in):: hunit,h(3,3,0:1),tag(natm),ra(3,natm)
     &     ,va(3,natm),eki(natm),epi(natm),strs(3,3,natm)

      integer:: ia,ib,l,i

      open(ionum,file=trim(cfname),status='replace')
      write(ionum,'(es23.14e3)') hunit
      write(ionum,'(3es23.14e3)') (((h(ia,ib,l)/hunit,ia=1,3)
     &     ,ib=1,3),l=0,1)
      write(ionum,'(i10)') natm
      do i=1,natm
        write(ionum,'(7es23.14e3,11es12.4)') tag(i),ra(1:3,i)
     &       ,va(1:3,i),eki(i)
     &       ,epi(i),strs(1:3,1:3,i)!*2.94d+4
      enddo
      close(ionum)

      end subroutine write_pmd_ascii
c=======================================================================
      subroutine read_pmd_ascii(ionum,cfname
     &     ,namax,natm,h,hunit,tag,ra,va,eki,epi,strs)
      implicit none
      integer,intent(in):: ionum,namax
      character(len=*),intent(in):: cfname
      integer,intent(out):: natm
      real(8),intent(out):: hunit,h(3,3,0:1),tag(namax),ra(3,namax)
     &     ,va(3,namax),eki(namax),epi(namax),strs(3,3,namax)
      
      integer:: ia,ib,l,i

      open(ionum,file=trim(cfname),status='old')
c-----natm: num. of particles in this node
      read(ionum,*) hunit
      read(ionum,*) (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
      h(1:3,1:3,0:1)= h(1:3,1:3,0:1)*hunit
      read(ionum,*) natm
      read(ionum,*) (tag(i),ra(1:3,i),va(1:3,i)
     &     ,eki(i),epi(i),strs(1:3,1:3,i),i=1,natm)
      close(ionum)

      end subroutine read_pmd_ascii
c=======================================================================
