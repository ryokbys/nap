      subroutine read_input(ionum,cfname)
c
c  Read frexible input format
c
      implicit none
      integer,intent(in):: ionum
      character(len=*),intent(in):: cfname
      character(len=128):: c1st

c      write(6,'(a)') ' Start reading '//trim(cfname)
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
 10   close(ionum)
c 10   write(6,'(a)') " Finished reading "//trim(cfname)
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
      elseif( trim(cname).eq.'sample_list' ) then
        call read_c1(ionum,csmplist)
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
      elseif( trim(cname).eq.'main_directory' .or.
     &       trim(cname).eq.'sample_directory' ) then
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
      elseif( trim(cname).eq.'numtol' ) then
        call read_i1(ionum,numtol)
        return
      elseif( trim(cname).eq.'gtol' ) then
        call read_r1(ionum,gtol)
        return
      elseif( trim(cname).eq.'atom_energy' ) then
        backspace(ionum)
        read(ionum,*) ctmp,itmp,eatom(itmp)
        return
      elseif( trim(cname).eq.'reference_structure' ) then
        call read_c1(ionum,crefstrct)
        return
      elseif( trim(cname).eq.'energy_match' ) then
        call read_l1(ionum,lematch)
        return
      elseif( trim(cname).eq.'force_match' ) then
        call read_l1(ionum,lfmatch)
        return
      elseif( trim(cname).eq.'force_match_species' ) then
        call read_force_match_species(ionum)
        return
      elseif( trim(cname).eq.'stress_match' ) then
        call read_l1(ionum,lsmatch)
        return
      elseif( trim(cname).eq.'loss_function' ) then
        call read_c1(ionum,ctype_loss)
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
      elseif( trim(cname).eq.'subtract_force_field' .or.
     &       trim(cname).eq.'additive_force_field' ) then
        call read_force_field(ionum)
        return
      elseif( trim(cname).eq.'rcut_other_FF' ) then
        call read_r1(ionum,rc_other)
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
      elseif( trim(cname).eq.'num_forces' ) then
        call read_i1(ionum,nfpsmpl)
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
      elseif( trim(cname).eq.'fs_num_inner' ) then
        call read_i1(ionum,ninnergfs)
        return
      elseif( trim(cname).eq.'fs_mode' ) then
        call read_c1(ionum,cfsmode)
        return
      elseif( trim(cname).eq.'fs_read_mask' ) then
        call read_c1(ionum,cread_fsmask)
        return
      elseif( trim(cname).eq.'fs_xrefresh' ) then
        call read_c1(ionum,cfs_xrefresh)
        return
      elseif( trim(cname).eq.'fs_max_num_refresh' ) then
        call read_i1(ionum,maxfsrefresh)
        return
      elseif( trim(cname).eq.'eps_energy' ) then
        call read_r1(ionum,epse)
        return
      elseif( trim(cname).eq.'eps_force' ) then
        call read_r1(ionum,epsf)
        return
      elseif( trim(cname).eq.'armijo_xi' ) then
        call read_r1(ionum,armijo_xi)
        return
      elseif( trim(cname).eq.'armijo_tau' ) then
        call read_r1(ionum,armijo_tau)
        return
      elseif( trim(cname).eq.'armijo_maxiter' ) then
        call read_i1(ionum,armijo_maxiter)
        return
      elseif( trim(cname).eq.'niter_linmin' ) then
        call read_i1(ionum,niter_linmin)
        return
      elseif( trim(cname).eq.'sa_temp_control' ) then
        call read_c1(ionum,sa_tctrl)
        return
      elseif( trim(cname).eq.'sa_temperature' ) then
        call read_r1(ionum,sa_temp0)
        return
      elseif( trim(cname).eq.'sa_tau' ) then
        call read_r1(ionum,sa_tau)
        return
      elseif( trim(cname).eq.'sa_dxwidth' ) then
        call read_r1(ionum,sa_xw0)
        return
      elseif( trim(cname).eq.'sa_div_best' ) then
        call read_r1(ionum,sa_div_best)
        return
      elseif( trim(cname).eq.'md_height' ) then
        call read_r1(ionum,md_height)
        return
      elseif( trim(cname).eq.'md_sigma' ) then
        call read_r1(ionum,md_sigma)
        return
      elseif( trim(cname).eq.'md_num_gaussian' ) then
        call read_i1(ionum,md_ng)
        return
      elseif( trim(cname).eq.'ga_temperature' ) then
        call read_r1(ionum,ga_temp)
        return
      elseif( trim(cname).eq.'ga_num_bits' ) then
        call read_i1(ionum,ga_nbits)
        return
      elseif( trim(cname).eq.'ga_num_individuals' ) then
        call read_i1(ionum,ga_nindivs)
        return
      elseif( trim(cname).eq.'ga_num_offsprings' ) then
        call read_i1(ionum,ga_noffsp)
        return
      elseif( trim(cname).eq.'ga_mutation_rate' ) then
        call read_r1(ionum,ga_rate_mutate)
        return
      elseif( trim(cname).eq.'ga_fitness' ) then
        call read_c1(ionum,ga_fitness)
        return
      elseif( trim(cname).eq.'de_fitness' ) then
        call read_c1(ionum,de_fitness)
        return
      elseif( trim(cname).eq.'de_algorithm' ) then
        call read_c1(ionum,de_algo)
        return
      elseif( trim(cname).eq.'de_num_individuals' ) then
        call read_i1(ionum,de_nindivs)
        return
      elseif( trim(cname).eq.'de_fraction' ) then
        call read_r1(ionum,de_frac)
        return
      elseif( trim(cname).eq.'de_lambda' ) then
        call read_r1(ionum,de_lambda)
        return
      elseif( trim(cname).eq.'de_crossover_rate' ) then
        call read_r1(ionum,de_cross_rate)
        return
      elseif( trim(cname).eq.'de_temperature' ) then
        call read_r1(ionum,de_temp)
        return
      elseif( trim(cname).eq.'de_wmin' ) then
        call read_r1(ionum,de_wmin)
        return
      elseif( trim(cname).eq.'de_wmax' ) then
        call read_r1(ionum,de_wmax)
        return
      elseif( trim(cname).eq.'pso_num_individuals' ) then
        call read_i1(ionum,pso_nindivs)
        return
      elseif( trim(cname).eq.'pso_w' ) then
        call read_r1(ionum,pso_w)
        return
      elseif( trim(cname).eq.'pso_c1' ) then
        call read_r1(ionum,pso_c1)
        return
      elseif( trim(cname).eq.'pso_c2' ) then
        call read_r1(ionum,pso_c2)
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
      elseif( trim(cname).eq.'cg_beta_type' ) then
        call read_i1(ionum,icgbtype)
        return
      elseif( trim(cname).eq.'lbfgs_history' ) then
        call read_i1(ionum,mstore)
        return
      elseif( trim(cname).eq.'nsmpl_outfrc' ) then
        call read_i1(ionum,nsmpl_outfrc)
        return
      elseif( trim(cname).eq.'sample_error' ) then
        backspace(ionum)
        read(ionum,*) ctmp,nserr
        allocate(cserr(nserr),seerr(nserr),sferr(nserr),sserr(nserr))
        call read_smpl_err(ionum,nserr,cserr,seerr,sferr,sserr)
        return
      elseif( trim(cname).eq.'sample_weight' ) then
        backspace(ionum)
        read(ionum,*) ctmp,nswgt
        allocate(cswgt(nswgt),swerg0(nswgt),swdenom(nswgt))
        call read_smpl_wgt(ionum,nswgt,cswgt,swerg0,swdenom)
        return
      elseif( trim(cname).eq.'fval_upper_limit' ) then
        call read_r1(ionum,fupper_lim)
        return
      elseif( trim(cname).eq.'force_limit' ) then
        call read_r1(ionum,force_limit)
        return
      elseif( trim(cname).eq.'NN_num_layers' ) then
        call read_i1(ionum,nn_nl)
        return
      elseif( trim(cname).eq.'NN_num_nodes' ) then
        call read_nn_nhl(ionum)
        return
      elseif( trim(cname).eq.'NN_sigtype' ) then
        call read_i1(ionum,nn_sigtype)
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
      subroutine read_smpl_err(ionum,nrow,cval,eerr,ferr,serr)
c
c  Read sample errors
c
      implicit none
      integer,intent(in):: ionum,nrow
      real(8),intent(out):: eerr(nrow),ferr(nrow),serr(nrow)
      character(len=*),intent(out):: cval(nrow)

      integer,external:: num_data
      integer:: irow,ndat
      character(len=1024):: ctmp 
      
      read(ionum,'(a)') ctmp
      ndat = num_data(trim(ctmp),' ')

      backspace(ionum)
      if( ndat.eq.3 ) then
        do irow=1,nrow
          read(ionum,*) cval(irow),eerr(irow),ferr(irow)
          serr(irow) = 1d0
        enddo
      else if( ndat.eq.4 ) then
        do irow=1,nrow
          read(ionum,*) cval(irow),eerr(irow),ferr(irow),serr(irow)
        enddo
      endif

      end subroutine read_smpl_err
c=======================================================================
      subroutine read_smpl_wgt(ionum,nrow,cval,erg0,temp)
c
c  Read sample weights
c
      implicit none
      integer,intent(in):: ionum,nrow
      real(8),intent(out):: erg0(nrow),temp(nrow)
      character(len=*),intent(out):: cval(nrow)

      integer,external:: num_data
      integer:: irow,ndat
      character(len=1024):: ctmp

      do irow=1,nrow
        read(ionum,*) cval(irow),erg0(irow),temp(irow)
      enddo

      end subroutine read_smpl_wgt
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
c=======================================================================
      subroutine read_force_field(ionum)
c
c     Read forces
c     There is no limit of number of force-fields to be specified.
c
      use variables,only: nsubff,csubffs
      implicit none
      integer,intent(in):: ionum

      integer:: i,ndat
      character(len=1024):: ctmp
      integer,external:: num_data

      backspace(ionum)
      read(ionum,'(a)') ctmp
      ndat = num_data(trim(ctmp),' ')
      nsubff = ndat -1
      allocate(csubffs(nsubff))
      backspace(ionum)
      read(ionum,*) ctmp, (csubffs(i),i=1,nsubff)
      end subroutine read_force_field
c=======================================================================
      subroutine read_force_match_species(ionum)
      use variables,only: lsps_frc
      implicit none
      integer,intent(in):: ionum

      integer:: i,ndat,nsp_frc
      character(len=1024):: ctmp
      integer,external:: num_data
      integer,allocatable:: isps(:)

      backspace(ionum)
      read(ionum,'(a)') ctmp
      ndat = num_data(trim(ctmp),' ')
      nsp_frc = ndat -1
      allocate(isps(nsp_frc))
      backspace(ionum)
      read(ionum,*) ctmp, (isps(i),i=1,nsp_frc)

      lsps_frc(:) = .false.
      do i=1,nsp_frc
        lsps_frc(isps(i)) = .true.
      enddo
      deallocate(isps)
      end subroutine read_force_match_species
c=======================================================================
      subroutine read_nn_nhl(ionum)
      use variables,only: nn_nhl
      implicit none
      integer,intent(in):: ionum

      integer:: nl,ndat,i
      character(len=1024):: ctmp
      integer,external:: num_data

      backspace(ionum)
      read(ionum,'(a)') ctmp
      ndat = num_data(trim(ctmp),' ')
      nl = ndat -2
      backspace(ionum)
      read(ionum,*) ctmp, (nn_nhl(i),i=0,nl)
      
      return
      end subroutine read_nn_nhl

