      subroutine write_pmd0_bin(ionum,cfname,cstat,natm,tag,ra,va,h
     &     ,hunit,eki,epi,strs)
      integer,intent(in):: ionum,natm
      real(8),intent(in):: tag(natm),ra(3,natm),va(3,natm)
     &     ,eki(natm),epi(natm),strs(3,3,natm),h(3,3,0:1),hunit
      character(len=*),intent(in):: cfname,cstat

      integer:: ia,ib,l,i

      open(ionum,file=trim(cfname),form='unformatted'
     &     ,status=trim(cstat))
      write(ionum) hunit
      write(ionum) (((h(ia,ib,l)/hunit,ia=1,3),ib=1,3),l=0,1)
      write(ionum) natm
      write(ionum) (tag(i),ra(1:3,i),va(1:3,i)
     &     ,eki(i),epi(i),strs(1:3,1:3,i),i=1,natm)
      close(ionum)

      end subroutine write_pmd0_bin
c=======================================================================
      subroutine write_pmd0_ascii(ionum,cfname,cstat,natm,tag,ra,va,h
     &     ,hunit,eki,epi,strs)
      integer,intent(in):: ionum,natm
      real(8),intent(in):: tag(natm),ra(3,natm),va(3,natm)
     &     ,eki(natm),epi(natm),strs(3,3,natm),h(3,3,0:1),hunit
      character(len=*),intent(in):: cfname,cstat

      integer:: ia,ib,l,i

      open(ionum,file=trim(cfname),status=trim(cstat))
      write(ionum,'(es23.14e3)') hunit
      write(ionum,'(3es23.14e3)') (((h(ia,ib,l)/hunit,ia=1,3)
     &     ,ib=1,3),l=0,1)
      write(ionum,'(i10)') natm
      write(ionum,'(7es23.14e3,11es12.4)') (tag(i),ra(1:3,i),va(1:3,i)
     &     ,eki(i),epi(i),strs(1:3,1:3,i),i=1,natm)
      close(ionum)

      end subroutine write_pmd0_ascii
c=======================================================================
      subroutine write_pmd_bin(ionum,cfname
     &     ,natm,h,hunit,tag,ra,va,eki,epi,strs,sorg,dt)
      implicit none
      include './params_unit.h'
      integer,intent(in):: ionum
      character(len=*),intent(in) :: cfname
      integer,intent(in):: natm
      real(8),intent(in):: h(3,3,0:1),tag(natm),ra(3,natm),va(3,natm)
     &     ,eki(3,3,natm),epi(natm),strs(3,3,natm),sorg(3),dt,hunit

      integer:: ia,ib,l,i

      open(ionum,file=cfname,form='unformatted'
     &     ,status='replace')
      write(ionum) hunit
      write(ionum) (((h(ia,ib,l)/hunit,ia=1,3),ib=1,3),l=0,1)
      write(ionum) natm
      write(ionum) (tag(i),ra(1:3,i)+sorg(1:3),va(1:3,i)/dt
     &     ,eki(1,1,i)+eki(2,2,i)+eki(3,3,i)
     &     ,epi(i),strs(1:3,1:3,i),i=1,natm)
      close(ionum)

      end subroutine write_pmd_bin
c=======================================================================
      subroutine read_pmd_bin(ionum,cfname
     &     ,namax,natm,h,hunit,tag,ra,va,eki,epi,strs)
      implicit none
      integer,intent(in):: ionum,namax
      character(len=*),intent(in):: cfname
      integer,intent(out):: natm
      real(8),intent(out):: hunit,h(3,3,0:1),tag(namax),ra(3,namax)
     &     ,va(3,namax),eki(3,3,namax),epi(namax),strs(3,3,namax)

      integer:: ia,ib,l,i
      
      open(ionum,file=trim(cfname),form='unformatted',status='old')
c-----natm: num. of particles in this node
      read(ionum) hunit
      read(ionum) (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
      h(1:3,1:3,0:1)= h(1:3,1:3,0:1)*hunit
      read(ionum) natm
      read(ionum) (tag(i),ra(1:3,i),va(1:3,i)
     &     ,eki(1,1,i),epi(i),strs(1:3,1:3,i),i=1,natm)
      close(ionum)

      end subroutine read_pmd_bin
c=======================================================================
      subroutine write_pmd_ascii(ionum,cfname
     &     ,natm,h,hunit,tag,ra,va,eki,epi,strs,sorg,dt)
      implicit none
      include './params_unit.h'
      integer,intent(in):: ionum
      character(len=*),intent(in) :: cfname
      integer,intent(in):: natm
      real(8),intent(in):: h(3,3,0:1),tag(natm),ra(3,natm),va(3,natm)
     &     ,eki(3,3,natm),epi(natm),strs(3,3,natm),sorg(3),dt,hunit

      integer:: ia,ib,l,i

      open(ionum,file=cfname,status='replace')
      write(ionum,'(es23.14e3)') hunit
      write(ionum,'(3es23.14e3)') (((h(ia,ib,l)/hunit,ia=1,3)
     &     ,ib=1,3),l=0,1)
      write(ionum,'(i10)') natm
      write(ionum,'(7es23.14e3,11es22.14)') (tag(i),ra(1:3,i)+sorg(1:3)
     &     ,va(1:3,i)/dt
     &     ,eki(1,1,i)+eki(2,2,i)+eki(3,3,i)
     &     ,epi(i),strs(1:3,1:3,i),i=1,natm)
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
     &     ,va(3,namax),eki(3,3,namax),epi(namax),strs(3,3,namax)
      
      integer:: ia,ib,l,i

      open(ionum,file=trim(cfname),status='old')
c-----natm: num. of particles in this node
      read(ionum,*) hunit
      read(ionum,*) (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
      h(1:3,1:3,0:1)= h(1:3,1:3,0:1)*hunit
      read(ionum,*) natm
      read(ionum,*) (tag(i),ra(1:3,i),va(1:3,i)
     &     ,eki(1,1,i),epi(i),strs(1:3,1:3,i),i=1,natm)
      close(ionum)

      end subroutine read_pmd_ascii
c=======================================================================
      subroutine read_pmdtot_ascii(ionum,cfname)
      use variables
      implicit none
      integer,intent(in):: ionum
      character(len=*),intent(in):: cfname

      integer:: ia,ib,l,i

      open(ionum,file=trim(cfname),status='old')
      read(ionum,*) hunit
      read(ionum,*) (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
      h(1:3,1:3,0:1)= h(1:3,1:3,0:1)*hunit
      read(ionum,*) ntot
      allocate(tagtot(ntot),rtot(3,ntot),vtot(3,ntot),epitot(ntot)
     &     ,ekitot(3,3,ntot),stot(3,3,ntot),atot(3,ntot))
      read(ionum,*) (tagtot(i),rtot(1:3,i),vtot(1:3,i)
     &     ,ekitot(1,1,i),epitot(i),stot(1:3,1:3,i),i=1,ntot)
      close(ionum)
      
      end subroutine read_pmdtot_ascii
c=======================================================================
      subroutine write_pmdtot_ascii(ionum,cfname)
      use variables
      implicit none
      include './params_unit.h'
      integer,intent(in):: ionum
      character(len=*),intent(in) :: cfname

      integer:: ia,ib,l,i

      open(ionum,file=cfname,status='replace')
      write(ionum,'(es23.14e3)') hunit
      write(ionum,'(3es23.14e3)') (((h(ia,ib,l)/hunit,ia=1,3)
     &     ,ib=1,3),l=0,1)
      write(ionum,'(i10)') ntot
      write(ionum,'(7es23.14e3,11es22.14)') (tagtot(i)
     &     ,rtot(1:3,i)
     &     ,vtot(1:3,i)/dt
     &     ,ekitot(1,1,i)+ekitot(2,2,i)+ekitot(3,3,i)
     &     ,epitot(i),stot(1:3,1:3,i),i=1,ntot)
      close(ionum)

      end subroutine write_pmdtot_ascii
c=======================================================================
      subroutine read_pmdtot_bin(ionum,cfname)
      use variables
      implicit none
      integer,intent(in):: ionum
      character(len=*),intent(in):: cfname

      integer:: ia,ib,l,i
      
      open(ionum,file=trim(cfname),form='unformatted',status='old')
c-----natm: num. of particles in this node
      read(ionum) hunit
      read(ionum) (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
      h(1:3,1:3,0:1)= h(1:3,1:3,0:1)*hunit
      read(ionum) ntot
      allocate(tagtot(ntot),rtot(3,ntot),atot(3,ntot)
     &     ,vtot(3,ntot),epitot(ntot)
     &     ,ekitot(3,3,ntot),stot(3,3,ntot))
      read(ionum) (tagtot(i),rtot(1:3,i),vtot(1:3,i)
     &     ,ekitot(1,1,i),epitot(i),stot(1:3,1:3,i),i=1,ntot)
      close(ionum)

      end subroutine read_pmdtot_bin
c=======================================================================
      subroutine write_pmdtot_bin(ionum,cfname)
      use variables
      implicit none
      include './params_unit.h'
      integer,intent(in):: ionum
      character(len=*),intent(in) :: cfname

      integer:: ia,ib,l,i

      open(ionum,file=cfname,form='unformatted'
     &     ,status='replace')
      write(ionum) hunit
      write(ionum) (((h(ia,ib,l)/hunit,ia=1,3),ib=1,3),l=0,1)
      write(ionum) ntot
      write(ionum) (tagtot(i),rtot(1:3,i),vtot(1:3,i)/dt
     &     ,ekitot(1,1,i)+ekitot(2,2,i)+ekitot(3,3,i)
     &     ,epitot(i),stot(1:3,1:3,i),i=1,ntot)
      close(ionum)

      end subroutine write_pmdtot_bin
c=======================================================================
