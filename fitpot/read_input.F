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
        call set_variable(ionum,c1st)
      enddo
      close(ionum)
 10   write(6,'(a)') " Finished reading "//trim(cfname)
      end subroutine read_input
c=======================================================================
      subroutine set_variable(ionum,cname)
      use variables
      use minimize
      use random
      implicit none
      integer,intent(in):: ionum
      character(len=*),intent(in):: cname

      character(len=128):: ctmp
      integer:: ndata,nrow,is,itmp
      
      if( trim(cname).eq.'num_samples' ) then
        call read_i1(ionum,nsmpl)
        return
      elseif( trim(cname).eq.'num_iteration' ) then
        call read_i1(ionum,niter)
        return
      elseif( trim(cname).eq.'num_iter_eval' ) then
        call read_i1(ionum,niter_eval)
        return
      elseif( trim(cname).eq.'fitting_method' ) then
        call read_c1(ionum,cfmethod)
        return
      elseif( trim(cname).eq.'main_directory' ) then
        call read_c1(ionum,cmaindir)
        return
      elseif( trim(cname).eq.'param_file' ) then
        call read_c1(ionum,cparfile)
        return
      elseif( trim(cname).eq.'run_mode' ) then
        call read_c1(ionum,crunmode)
        return
      elseif( trim(cname).eq.'evaluation_type' ) then
        call read_c1(ionum,cevaltype)
        return
      elseif( trim(cname).eq.'xtol' ) then
        call read_r1(ionum,xtol)
        return
      elseif( trim(cname).eq.'ftol' ) then
        call read_r1(ionum,ftol)
        return
      elseif( trim(cname).eq.'gtol' ) then
        call read_r1(ionum,gtol)
        return
      elseif( trim(cname).eq.'atom_energy' ) then
        backspace(ionum)
        read(ionum,*) ctmp,itmp,eatom(itmp)
        return
      elseif( trim(cname).eq.'force_match' ) then
        call read_l1(ionum,lfmatch)
        return
      elseif( trim(cname).eq.'penalty' ) then
        call read_c1(ionum,cpena)
        return
      elseif( trim(cname).eq.'penalty_weight' ) then
        call read_r1(ionum,pwgt)
        return
      elseif( trim(cname).eq.'potential' ) then
        call read_c1(ionum,cpot)
        return
      elseif( trim(cname).eq.'gradient' ) then
        call read_l1(ionum,lgrad)
        return
      elseif( trim(cname).eq.'grad_scale' ) then
        call read_l1(ionum,lgscale)
        return
      elseif( trim(cname).eq.'gscale_factor' ) then
        call read_r1(ionum,gscl)
        return
      elseif( trim(cname).eq.'normalize_input' ) then
        call read_c1(ionum,cnormalize)
        return
      elseif( trim(cname).eq.'force_scale' ) then
        call read_l1(ionum,lfscale)
        return
      elseif( trim(cname).eq.'fscale_factor' ) then
        call read_r1(ionum,fscl)
        return
      elseif( trim(cname).eq.'freduce_threshold' ) then
        call read_r1(ionum,fred)
        return
      elseif( trim(cname).eq.'num_forces' ) then
        call read_i1(ionum,nfpsmpl)
        return
      elseif( trim(cname).eq.'sample_weight' ) then
        call read_l1(ionum,lswgt)
        return
      elseif( trim(cname).eq.'sample_weight_erg' ) then
        call read_r1(ionum,swerg)
        return
      elseif( trim(cname).eq.'num_multiprocess' ) then
        call read_i1(ionum,nprcs)
        return
      elseif( trim(cname).eq.'coeff_sequential' ) then
        call read_r1(ionum,seqcoef)
        return
      elseif( trim(cname).eq.'print_level' ) then
        call read_i1(ionum,iprint)
        return
      elseif( trim(cname).eq.'line_minimization' ) then
        call read_c1(ionum,clinmin)
        return
      elseif( trim(cname).eq.'test_ratio' ) then
        call read_r1(ionum,ratio_test)
        return
      elseif( trim(cname).eq.'num_iter_gfs' ) then
        call read_i1(ionum,nitergfs)
        return
      elseif( trim(cname).eq.'eps_energy' ) then
        call read_r1(ionum,epse)
        return
      elseif( trim(cname).eq.'eps_force' ) then
        call read_r1(ionum,epsf)
        return
      elseif( trim(cname).eq.'sa_temperature' ) then
        call read_r1(ionum,sa_temp0)
        return
      elseif( trim(cname).eq.'sa_dxwidth' ) then
        call read_r1(ionum,sa_xw0)
        return
      elseif( trim(cname).eq.'random_seed' ) then
        call read_r1(ionum,rseed)
        return
      elseif( trim(cname).eq.'sgd_update' ) then
        call read_c1(ionum,csgdupdate)
        return
      elseif( trim(cname).eq.'sgd_batch_size' ) then
        call read_i1(ionum,nsgdbsize)
        return
      elseif( trim(cname).eq.'sgd_rate0' ) then
        call read_r1(ionum,r0sgd)
        return
      elseif( trim(cname).eq.'init_params' ) then
        call read_c1(ionum,cinitv)
        return
      elseif( trim(cname).eq.'init_params_sgm' ) then
        call read_r1(ionum,vinitsgm)
        return
      elseif( trim(cname).eq.'init_params_mu' ) then
        call read_r1(ionum,vinitmu)
        return
      elseif( trim(cname).eq.'init_params_rs' ) then
        call read_r1(ionum,vinitrs)
        return
      elseif( trim(cname).eq.'individual_weight' ) then
        backspace(ionum)
        read(ionum,*) ctmp,nwgtindiv
        allocate(cwgtindiv(nwgtindiv),wgtindiv(nwgtindiv))
        call read_cr(ionum,nwgtindiv,cwgtindiv,wgtindiv)
        return
      elseif( trim(cname).eq.'sample_error' ) then
        backspace(ionum)
        read(ionum,*) ctmp,nserr
        allocate(cserr(nserr),seerr(nserr),sferr(nserr))
        call read_crr(ionum,nserr,cserr,seerr,sferr)
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
c      write(6,'(1x,a,es15.3)') trim(ctmp),rval

      end subroutine read_r1
c=======================================================================
      subroutine read_cr(ionum,nrow,cval,rval)
c
c  Read sets of (character, real*8)
c
      integer,intent(in):: ionum,nrow
      real(8),intent(out):: rval(nrow)
      character(len=*),intent(out):: cval(nrow)

      do irow=1,nrow
        read(ionum,*) cval(irow),rval(irow)
      enddo

      end subroutine read_cr
c=======================================================================
      subroutine read_crr(ionum,nrow,cval,rval1,rval2)
c
c  Read sets of (character, real*8, real*8)
c
      integer,intent(in):: ionum,nrow
      real(8),intent(out):: rval1(nrow),rval2(nrow)
      character(len=*),intent(out):: cval(nrow)

      do irow=1,nrow
        read(ionum,*) cval(irow),rval1(irow),rval2(irow)
      enddo

      end subroutine read_crr
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
c      write(6,'(1x,a,5x,l1)') trim(ctmp),lval

      end subroutine read_l1

