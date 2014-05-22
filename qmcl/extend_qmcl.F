      program extend_pmd
      implicit none
      integer,parameter:: namax = 1000000
      integer:: natm,inc,ix,iy,iz,nex,ney,nez,m,is,ifmv
      real(8):: hunit,h(3,3,0:1),tag(namax),ra(3,namax),va(3,namax)
     &     ,eki(namax),epi(namax),strs(3,3,namax),rat(3,namax)
      character:: cifname*128,cofname*128

c.....ask num of extension
      write(6,'(a)') 'type input-file-name to be expanded:'
      read(5,*) cifname
      write(6,'(a)') 'number of extension for each axis, x,y,z:'
      read(5,*) nex,ney,nez
      write(6,'(a)') 'type output-file-name'
      read(5,*) cofname

      call read_pmd_ascii(10,trim(cifname),namax,natm,h,hunit,tag,ra,va
     &     ,eki,epi,strs)

      if( natm*nex*ney*nez .gt. namax ) then
        write(6,'(a)') '[error] natm*nex*ney*nez .gt. namax !!!'
        write(6,'(a,i10)') ' namax should be greater than '
     &       ,natm*nex*ney*nez
        stop
      endif

      h(1,1,0)= h(1,1,0)*nex
      h(1,2,0)= h(1,2,0)*nex
      h(1,3,0)= h(1,3,0)*nex
      h(2,1,0)= h(2,1,0)*ney
      h(2,2,0)= h(2,2,0)*ney
      h(2,3,0)= h(2,3,0)*ney
      h(3,1,0)= h(3,1,0)*nez
      h(3,2,0)= h(3,2,0)*nez
      h(3,3,0)= h(3,3,0)*nez
      
      inc= 0
      do ix=1,nex
        do iy=1,ney
          do iz=1,nez
            do m=1,natm
              inc=inc +1
              rat(1,inc)= (ra(1,m) +dble(ix-1))/nex
              rat(2,inc)= (ra(2,m) +dble(iy-1))/ney
              rat(3,inc)= (ra(3,m) +dble(iz-1))/nez
              va(1:3,inc)= va(1:3,m)
              eki(inc)= epi(m)
              epi(inc)= epi(m)
              strs(1:3,1:3,inc)= strs(1:3,1:3,m)
              is= int(tag(m))
              ifmv= int((tag(m)-is)*10)
              tag(inc)= 1d0*is +0.1d0*ifmv +1d-14*inc
            enddo
          enddo
        enddo
      enddo
      
      call write_pmd_ascii(11,trim(cofname),inc,h,hunit,tag,ra,va
     &     ,eki,epi,strs)

      write(6,'(a)') 'program done.'
      write(6,'(a)') 'check output file: '//trim(cofname)
      
      end program extend_pmd
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "make extpmd"
c     End:
