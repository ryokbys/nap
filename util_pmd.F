      subroutine write_pmd(ionum,cfname,cstat,natm,tag,ra,va,h
     &     ,eki,epi,strs)
      integer,intent(in):: ionum,natm
      real(8),intent(in):: tag(natm),ra(3,natm),va(3,natm)
     &     ,eki(natm),epi(natm),strs(3,3,natm),h(3,3,0:1)
      character(len=*),intent(in):: cfname,cstat

      integer:: ia,ib,l,i

      open(ionum,file=trim(cfname),form='unformatted'
     &     ,status=trim(cstat))
      write(ionum) natm
      write(ionum) (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
      write(ionum) (tag(i),ra(1:3,i),va(1:3,i)
     &     ,eki(i),epi(i),strs(1:3,1:3,i),i=1,natm)
      close(ionum)

      end subroutine write_pmd
