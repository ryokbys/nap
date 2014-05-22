      subroutine read_input(ionum,cfname)
c
c  Read frexible input format
c
      implicit none
      integer,intent(in):: ionum
      character(len=*),intent(in):: cfname
      character(len=128):: c1st

      open(ionum,file=trim(cfname))
      do
c.....Read 1st word in each line
        read(ionum,*,end=10) c1st
c.....Skip comment line
        if( c1st(1:1).eq.'!' .or.
     &       c1st(1:1).eq.'#' .or.
c.....Skip lines starting from digits or sign
     &       c1st(1:1).eq.'0' .or.
     &       c1st(1:1).eq.'1' .or.
     &       c1st(1:1).eq.'2' .or.
     &       c1st(1:1).eq.'3' .or.
     &       c1st(1:1).eq.'4' .or.
     &       c1st(1:1).eq.'5' .or.
     &       c1st(1:1).eq.'6' .or.
     &       c1st(1:1).eq.'7' .or.
     &       c1st(1:1).eq.'8' .or.
     &       c1st(1:1).eq.'9' .or.
     &       c1st(1:1).eq.'+' .or.
     &       c1st(1:1).eq.'-' ) cycle
c        write(6,'(a)') c1st
        call set_variable(ionum,c1st)
      enddo
      close(ionum)
 10   write(6,'(a)') " Finished reading "//trim(cfname)
      end subroutine read_input
c=======================================================================
      subroutine set_variable(ionum,cname)
      use variables
      implicit none
      integer,intent(in):: ionum
      character(len=*),intent(in):: cname

      character(len=128):: ctmp
      integer:: ndata,nrow,is
      
      if( trim(cname).eq.'time_interval' ) then
        call read_r1(ionum,dt)
        return
      elseif( trim(cname).eq.'num_iteration' ) then
        call read_i1(ionum,nstp)
        return
      elseif( trim(cname).eq.'num_out_energy' ) then
        call read_i1(ionum,nerg)
        return
      elseif( trim(cname).eq.'flag_out_qmcl' ) then
        call read_l1(ionum,lpmd)
        return
      elseif( trim(cname).eq.'num_out_qmcl' ) then
        call read_i1(ionum,npmd)
        return
      elseif( trim(cname).eq.'flag_cosine_damping' ) then
        call read_l1(ionum,lcsdmp)
        return
      elseif( trim(cname).eq.'flag_damping' ) then
        call read_l1(ionum,ldmp)
        return
      elseif( trim(cname).eq.'damping_coeff' ) then
        call read_r1(ionum,dmp)
        return
      elseif( trim(cname).eq.'flag_temperature' ) then
        call read_l1(ionum,ltctl)
        return
      elseif( trim(cname).eq.'initial_temperature' ) then
        call read_r1(ionum,tinit)
        return
      elseif( trim(cname).eq.'temperature_target' ) then
        call read_r1(ionum,treq)
        return
      elseif( trim(cname).eq.'temperature_relax_time' ) then
        call read_r1(ionum,trlx)
        return
      elseif( trim(cname).eq.'cutoff_radius' ) then
        call read_r1(ionum,rcut)
        return
      elseif( trim(cname).eq.'flag_convergence' ) then
        call read_l1(ionum,lconv)
        return
      elseif( trim(cname).eq.'convergence_force' ) then
        call read_r1(ionum,feps)
        return
      elseif( trim(cname).eq.'factor_direction' ) then
        backspace(ionum)
        read(ionum,*) ctmp,ndata,nrow
        if( ndata.ne.3 ) write(6,'(a)') ' [Error] ndata.ne.3 !!!'
        call read_rs(ionum,ctmp,3,nrow,fmv(1:3,1:nrow))
        return
      elseif( trim(cname).eq.'io_format' ) then
        call read_c1(ionum,ciofmt)
        return
      elseif( trim(cname).eq.'force_type' ) then
        call read_c1(ionum,cfrct)
        return
      elseif( trim(cname).eq.'cluster_shape' ) then
        call read_c1(ionum,cqmshape)
        return
      elseif( trim(cname).eq.'cluster_center_x' ) then
        call read_r1(ionum,xc(1))
        return
      elseif( trim(cname).eq.'cluster_center_y' ) then
        call read_r1(ionum,xc(2))
        return
      elseif( trim(cname).eq.'cluster_center_z' ) then
        call read_r1(ionum,xc(3))
        return
      elseif( trim(cname).eq.'cluster_radius' ) then
        call read_r1(ionum,rclst)
        return
      elseif( trim(cname).eq.'qm_region_radius' ) then
        call read_r1(ionum,rcqm)
        return
      elseif( trim(cname).eq.'cluster_vacuum_x' ) then
        call read_r1(ionum,wvac(1))
        return
      elseif( trim(cname).eq.'cluster_lx' ) then
        call read_r1(ionum,rclx)
        return
      elseif( trim(cname).eq.'cluster_ly' ) then
        call read_r1(ionum,rcly)
        return
      elseif( trim(cname).eq.'cluster_lz' ) then
        call read_r1(ionum,rclz)
        return
      elseif( trim(cname).eq.'qm_region_lx' ) then
        call read_r1(ionum,rcqmx)
        return
      elseif( trim(cname).eq.'qm_region_ly' ) then
        call read_r1(ionum,rcqmy)
        return
      elseif( trim(cname).eq.'qm_region_lz' ) then
        call read_r1(ionum,rcqmz)
        return
      elseif( trim(cname).eq.'cluster_vacuum_y' ) then
        call read_r1(ionum,wvac(2))
        return
      elseif( trim(cname).eq.'cluster_vacuum_z' ) then
        call read_r1(ionum,wvac(3))
        return
      elseif( trim(cname).eq.'cluster_termination' ) then
        call read_l1(ionum,lterm)
        return
      elseif( trim(cname).eq.'run_type' ) then
        call read_c1(ionum,crtype)
        return
      elseif( trim(cname).eq.'mass' ) then
        backspace(ionum)
        read(ionum,*) ctmp,is,am(is)
        return
c      elseif( trim(cname).eq.'' ) then
c        call read_i1(ionum,nz)
c        return
      endif
      
c      write(6,'(a)') " [Error] No match: "//trim(cname)//" !!!"
c      stop
      write(6,'(a)') ' [Warning] No match: '//trim(cname)//' !!!'
      return

      end subroutine set_variable
c=======================================================================
      subroutine read_r1(ionum,rval)
c
c  Read one read*8 parameter from the line
c
      integer,intent(in):: ionum
      real(8),intent(out):: rval
      character(len=128):: ctmp

      backspace(ionum)
      read(ionum,*) ctmp,rval
      write(6,'(1x,a,es15.3)') trim(ctmp),rval

      end subroutine read_r1
c=======================================================================
      subroutine read_rs(ionum,ctmp,ndata,nrow,rval)
c
c  Read several read*8 parameters
c
      integer,intent(in):: ionum,ndata,nrow
      real(8),intent(out):: rval(ndata,nrow)
      character(len=*),intent(in):: ctmp

      write(6,'(1x,a,2i8)') trim(ctmp), ndata, nrow
      do n=1,nrow
        read(ionum,*) (rval(i,n),i=1,ndata)
        write(6,'(1x,100es15.3)') (rval(i,n),i=1,ndata)
      enddo

      end subroutine read_rs
c=======================================================================
      subroutine read_i1(ionum,ival)
c
c  Read one integer parameter from the line
c
      integer,intent(in):: ionum
      integer,intent(out):: ival
      character(len=128):: ctmp

      backspace(ionum)
      read(ionum,*) ctmp,ival
      write(6,'(1x,a,i10)') trim(ctmp),ival

      end subroutine read_i1
c=======================================================================
      subroutine read_c1(ionum,cval)
c
c  Read one word from the line
c
      integer,intent(in):: ionum
      character(len=*),intent(out):: cval
      character(len=128):: ctmp

      backspace(ionum)
      read(ionum,*) ctmp,cval
      write(6,'(1x,a,5x,a)') trim(ctmp),trim(cval)

      end subroutine read_c1
c=======================================================================
      subroutine read_l1(ionum,lval)
c
c  Read logical variable
c
      integer,intent(in):: ionum
      logical,intent(out):: lval
      character(len=128):: ctmp

      backspace(ionum)
      read(ionum,*) ctmp,lval
      write(6,'(1x,a,5x,l1)') trim(ctmp),lval

      end subroutine read_l1
c=======================================================================
