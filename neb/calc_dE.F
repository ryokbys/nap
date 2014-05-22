      program calc_dE
c-----------------------------------------------------------------------
c  Calculate the energy difference from the initial state 
c  by means of force integration.
c-----------------------------------------------------------------------
c  Usage:
c    $ ./calc_dE
c-----------------------------------------------------------------------
      use variables
      implicit none

      integer:: i,ierr,ionum,istp,is,itmp,l,ia,ib
      integer:: islc,islcp,islcm
      character(len=3):: cslc
      real(8):: dr(3),erg,ft(3)
      real(8),allocatable,dimension(:):: ergi
      real(8),external:: anint_rk

      call read_input(10,'in.neb')
      allocate(ra(3,natm,nslc),va(3,natm,nslc),fa(3,natm,nslc)
     &     ,tag(natm,nslc),ergi(natm))

c.....Read atom configuration and forces of each slice
      do islc=1,nslc
        ionum= 10+islc
        write(cslc(1:3),'(i3.3)') islc
c.....Atom config
        open(ionum,file='neb'//cslc,status='old')
        read(ionum,*) itmp
        if( itmp.ne.natm ) stop ' [Error] itmp.ne.natm !!!'
        read(ionum,*) (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
        do i=1,natm
          read(ionum,*) tag(i,islc),ra(1:3,i,islc)
        enddo
        close(ionum)
c.....Forces
        open(ionum,file='frc'//cslc,status='old')
        read(ionum,*) itmp
        if( itmp.ne.natm ) stop ' [Error] itmp.ne.natm !!!'
        do i=1,natm
          read(ionum,*) fa(1:3,i,islc)
        enddo
        close(ionum)
      enddo

      
c.....Integrate force*dr
      open(90,file='out.dE',status='replace')
      open(91,file='out.dEi',status='replace')
      erg= 0d0
      ergi(1:natm)= 0d0
      write(90,'(2es22.14)') erg,erg
      write(91,'(10000es12.4)')  (ergi(i),i=1,natm)
      do islc=1,nslc-1
        do i=1,natm
          dr(1:3)= ra(1:3,i,islc+1) -ra(1:3,i,islc)
          dr(1)= dr(1) -anint_rk(dr(1))
          dr(2)= dr(2) -anint_rk(dr(2))
          dr(3)= dr(3) -anint_rk(dr(3))
          dr(1:3)= h(1:3,1,0)*dr(1) +h(1:3,2,0)*dr(2) +h(1:3,3,0)*dr(3)
          ft(1:3)= (fa(1:3,i,islc)+fa(1:3,i,islc+1))/2
          erg= erg -ft(1)*dr(1)
     &         -ft(2)*dr(2)
     &         -ft(3)*dr(3)
          ergi(i)= ergi(i) -ft(1)*dr(1)
     &         -ft(2)*dr(2)
     &         -ft(3)*dr(3)
        enddo
        write(90,'(2es22.14)')  erg,ergi(natm)
        write(91,'(10000es12.4)')  (ergi(i),i=1,natm)
      enddo
      close(90)
      close(91)

      end program calc_dE
c=======================================================================
      function anint_rk(x)
      real(8),intent(in):: x
      real(8):: anint_rk
      anint_rk= 0d0
      if( x.gt.0.5d0 ) then
        anint_rk= 1.d0
      else if( x.lt.-0.5d0 ) then
        anint_rk= -1.d0
      endif
      return
      end function anint_rk
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make calc_dE"
c     End:
