      subroutine read_input(ionum,cfname)
c
c  Read frexible input format
c
      implicit none
      integer,intent(in):: ionum
      character(len=*),intent(in):: cfname
      character(len=128):: c1st

      write(6,'(a)') ' Start reading '//trim(cfname)

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
#ifdef __WALL__
      use wall
#endif
      implicit none
      integer,intent(in):: ionum
      character(len=*),intent(in):: cname

      character(len=128):: ctmp
      integer:: ndata,nrow,is,itmp
      
      if( trim(cname).eq.'time_interval' ) then
        call read_r1(ionum,dt)
        return
      elseif( trim(cname).eq.'num_iteration' .or.
     &       trim(cname).eq.'num_steps' ) then
        call read_i1(ionum,nstp)
        return
      elseif( trim(cname).eq.'num_out_energy' ) then
        call read_i1(ionum,nerg)
        return
      elseif( trim(cname).eq.'flag_out_pmd' ) then
        call read_i1(ionum,ifpmd)
        return
      elseif( trim(cname).eq.'num_out_pmd' ) then
        call read_i1(ionum,npmd)
        return
      elseif( trim(cname).eq.'cutoff_radius' ) then
        call read_r1(ionum,rc)
        return
      elseif( trim(cname).eq.'cutoff_buffer' ) then
        call read_r1(ionum,rbuf)
        return
      elseif( trim(cname).eq.'flag_damping' ) then
        call read_i1(ionum,ifdmp)
        return
      elseif( trim(cname).eq.'damping_coeff' ) then
        call read_r1(ionum,dmp)
        return
      elseif( trim(cname).eq.'initial_temperature' ) then
        call read_r1(ionum,tinit)
        return
      elseif( trim(cname).eq.'temperature_control' ) then
        call read_c1(ionum,ctctl)
        return
      elseif( trim(cname).eq.'temperature_target' ) then
        backspace(ionum)
        read(ionum,*) ctmp,itmp,ttgt(itmp)
c        call read_r1(ionum,ttgt)
        return
      elseif( trim(cname).eq.'temperature_relax_time' ) then
        call read_r1(ionum,trlx)
        return
c.....temperature distribution along x
      elseif( trim(cname).eq.'flag_temp_dist' ) then
        call read_l1(ionum,ltdst)
        return
      elseif( trim(cname).eq.'num_temp_dist' ) then
        call read_i1(ionum,ntdst)
        return
      elseif( trim(cname).eq.'num_nodes_x' ) then
        call read_i1(ionum,nx)
        return
      elseif( trim(cname).eq.'num_nodes_y' ) then
        call read_i1(ionum,ny)
        return
      elseif( trim(cname).eq.'num_nodes_z' ) then
        call read_i1(ionum,nz)
        return
      elseif( trim(cname).eq.'shear_stress' ) then
        call read_r1(ionum,shrst)
        return
      elseif( trim(cname).eq.'factor_direction' ) then
        backspace(ionum)
        read(ionum,*) ctmp,ndata,nrow
c        if( ndata.ne.3 ) write(6,'(a)') ' [Error] ndata.ne.3 !!!'
        call read_rs(ionum,ctmp,3,nrow,fmv(1:3,1:nrow))
        return
      elseif( trim(cname).eq.'pressure_target' ) then
        call read_r1(ionum,ptgt)
        return
      elseif( trim(cname).eq.'stress_target' ) then
        call read_rs(ionum,ctmp,3,3,stgt(1:3,1:3))
        return
      elseif( trim(cname).eq.'vol_mass_coeff' ) then
        call read_r1(ionum,vmcoeff)
        return
      elseif( trim(cname).eq.'vol_change_damping' ) then
        call read_r1(ionum,voldmp)
        return
      elseif( trim(cname).eq.'flag_isobaric' ) then
        call read_i1(ionum,ifpctl)
        return
      elseif( trim(cname).eq.'zload_type' ) then
        call read_c1(ionum,czload_type)
        return
      elseif( trim(cname).eq.'final_strain' ) then
        call read_r1(ionum,strfin)
        return
      elseif( trim(cname).eq.'mass' ) then
        backspace(ionum)
        read(ionum,*) ctmp,is,am(is)
        return
      elseif( trim(cname).eq.'io_format' ) then
        call read_c1(ionum,ciofmt)
        return
      elseif( trim(cname).eq.'force_type' ) then
        call read_c1(ionum,cforce)
        return
#ifdef __WALL__
      elseif( trim(cname).eq.'wall_pos_top' ) then
        call read_r1(ionum,wtop)
        return
      elseif( trim(cname).eq.'wall_pos_bottom' ) then
        call read_r1(ionum,wbot)
        return
      elseif( trim(cname).eq.'wall_target_pressure' ) then
        call read_r1(ionum,ptgt_wall)
        return
      elseif( trim(cname).eq.'wall_relax_time' ) then
        call read_r1(ionum,trlx_wall)
        return
      elseif( trim(cname).eq.'wall_nout' ) then
        call read_i1(ionum,nout_wall)
        return
#endif
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
c      write(6,'(1x,a,es15.3)') trim(ctmp),rval

      end subroutine read_r1
c=======================================================================
      subroutine read_rs(ionum,ctmp,ndata,nrow,rval)
c
c  Read several read*8 parameters
c
      integer,intent(in):: ionum,ndata,nrow
      real(8),intent(out):: rval(ndata,nrow)
      character(len=*),intent(in):: ctmp

c      write(6,'(1x,a,2i8)') trim(ctmp), ndata, nrow
      do n=1,nrow
        read(ionum,*) (rval(i,n),i=1,ndata)
c        write(6,'(1x,100es15.3)') (rval(i,n),i=1,ndata)
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
c      write(6,'(1x,a,i10)') trim(ctmp),ival

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
c      write(6,'(1x,2a)') trim(ctmp),trim(cval)

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

