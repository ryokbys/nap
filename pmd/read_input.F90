subroutine read_inpmd(ionum,cfname)
!
!  Read frexible input format
!
  implicit none
  integer,intent(in):: ionum
  character(len=*),intent(in):: cfname
  character(len=128):: c1st

  open(ionum,file=trim(cfname))
  do
!.....Read 1st word in each line
    read(ionum,*,end=10) c1st
!.....Skip comment line
    if( c1st(1:1).eq.'!' .or. &
         c1st(1:1).eq.'#' .or. &
!.....Skip lines starting from digits or sign 
         c1st(1:1).eq.'0' .or. &
         c1st(1:1).eq.'1' .or. &
         c1st(1:1).eq.'2' .or. &
         c1st(1:1).eq.'3' .or. &
         c1st(1:1).eq.'4' .or. &
         c1st(1:1).eq.'5' .or. &
         c1st(1:1).eq.'6' .or. &
         c1st(1:1).eq.'7' .or. &
         c1st(1:1).eq.'8' .or. &
         c1st(1:1).eq.'9' .or. &
         c1st(1:1).eq.'+' .or. &
         c1st(1:1).eq.'-' ) cycle
!        write(6,'(a)') c1st
    call set_variable(ionum,c1st)
  enddo
10 close(ionum)

end subroutine read_inpmd
!=======================================================================
subroutine set_variable(ionum,cname)
  use pmdvars
  use pmdmpi
  use util,only: csp2isp
  use force,only: ol_type, ol_force
  use extforce,only: lextfrc,cspc_extfrc,extfrc
  use clrchg,only: lclrchg,cspc_clrchg,clrfield,clr_init
  use localflux,only: lflux,nlx,nly,nlz,noutlflux
  use pdens,only: lpdens,cspc_pdens,npx,npy,npz,orig_pdens,hmat_pdens
#ifdef __WALL__
  use wall
#endif
  implicit none
  integer,intent(in):: ionum
  character(len=*),intent(in):: cname

  character(len=128):: ctmp,cval
  character(len=3):: csp
  integer:: ndata,nrow,is,itmp
  real(8):: tmp

  if( trim(cname).eq.'time_interval' ) then
    call read_r1(ionum,dt)
    return
  elseif( trim(cname).eq.'vardt_length_scale' .or. &
       trim(cname).eq.'fire_max_length' ) then
    call read_r1(ionum,vardt_len)
    return
  elseif( trim(cname).eq.'max_num_atoms' ) then
    call read_i1(ionum,namax)
    return
  elseif( trim(cname).eq.'max_num_boundary_atoms' ) then
    call read_i1(ionum,nbmax)
    return
  elseif( trim(cname).eq.'max_num_neighbors' ) then
    call read_i1(ionum,nnmax)
    return
  elseif( trim(cname).eq.'num_iteration' .or. &
       trim(cname).eq.'num_steps' ) then
    call read_i1(ionum,nstp)
    return
  elseif( trim(cname).eq.'min_iteration' .or. &
       trim(cname).eq.'min_steps' ) then
    call read_i1(ionum,minstp)
    return
  elseif( trim(cname).eq.'num_out_energy' ) then
    call read_i1(ionum,nerg)
    return
  elseif( trim(cname).eq.'flag_out_pmd' .or. &
       trim(cname).eq.'flag_out_pos' ) then
    call read_i1(ionum,ifpmd)
    return
  elseif( trim(cname).eq.'num_out_pmd' .or. &
       trim(cname).eq.'num_out_pos' ) then
    call read_i1(ionum,npmd)
    return
  elseif( trim(cname).eq.'dump_aux_order' ) then
    call read_dumpaux(ionum)
    return
  elseif( trim(cname).eq.'boundary' ) then
    backspace(ionum)
    read(ionum,*) ctmp, cval
    if( len(trim(cval)).ne.3 ) then
      print *,'WARNING: input format for boundary is wrong...'
      return
    endif
    boundary = trim(cval)
    return
  elseif( trim(cname).eq.'flag_sort' ) then
    call read_i1(ionum,ifsort)
    return
  elseif( trim(cname).eq.'cutoff_radius' ) then
    call read_r1(ionum,rc)
    return
  elseif( trim(cname).eq.'cutoff_radius_1nn' ) then
    call read_r1(ionum,rc1nn)
    return
  elseif( trim(cname).eq.'cutoff_buffer' ) then
    call read_r1(ionum,rbuf)
    return
  elseif( trim(cname).eq.'reorder_arrays' ) then
    call read_l1(ionum,lreorder)
    return
  elseif( trim(cname).eq.'flag_damping' ) then
    call read_i1(ionum,ifdmp)
    return
  elseif( trim(cname).eq.'minimization' ) then
    call read_c1(ionum,cmin)
    return
  elseif( trim(cname).eq.'damping_coeff' ) then
    call read_r1(ionum,dmp)
    return
  elseif( trim(cname).eq.'converge_eps' ) then
    call read_r1(ionum,eps_conv)
    return
  elseif( trim(cname).eq.'converge_num' ) then
    call read_i1(ionum,n_conv)
    return
  elseif( trim(cname).eq.'initial_temperature' ) then
    call read_r1(ionum,tinit)
    return
  elseif( trim(cname).eq.'final_temperature' ) then
    call read_r1(ionum,tfin)
    return
  elseif( trim(cname).eq.'temperature_control' ) then
    call read_c1(ionum,ctctl)
    return
  elseif( trim(cname).eq.'temperature_target' ) then
    backspace(ionum)
    read(ionum,*) ctmp,itmp,ttgt(itmp)
!        call read_r1(ionum,ttgt)
    return
  elseif( trim(cname).eq.'temperature_relax_time' ) then
    call read_r1(ionum,trlx)
    return
  elseif( trim(cname).eq.'remove_translation' ) then
    call read_i1(ionum,nrmtrans)
    return
  elseif( trim(cname).eq.'random_seed' ) then
    call read_r1(ionum,rseed)
    return
!.....temperature distribution along x
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
  elseif( trim(cname).eq.'num_omp_threads' ) then
    call read_i1(ionum,nomp)
    return
  elseif( trim(cname).eq.'shear_stress' ) then
    call read_r1(ionum,shrst)
    return
  elseif( trim(cname).eq.'factor_direction' ) then
    backspace(ionum)
    read(ionum,*) ctmp,ndata,nrow
!        if( ndata.ne.3 ) write(6,'(a)') ' [Error] ndata.ne.3 !!!'
    call read_rs(ionum,3,nrow,fmv(1:3,1:nrow))
    return
  elseif( trim(cname).eq.'pressure_target' ) then
    call read_r1(ionum,ptgt)
    return
  elseif( trim(cname).eq.'initial_pressure_target' ) then
    call read_r1(ionum,pini)
    return
  elseif( trim(cname).eq.'final_pressure_target' ) then
    call read_r1(ionum,pfin)
    return
  elseif( trim(cname).eq.'stress_target' ) then
    call read_rs(ionum,3,3,stgt(1:3,1:3))
    return
  elseif( trim(cname).eq.'stress_relax_time' .or. &
       trim(cname).eq.'pressure_relax_time' ) then
    call read_r1(ionum,srlx)
    return
  elseif( trim(cname).eq.'stress_control' ) then
    call read_c1(ionum,cpctl)
    return
  elseif( trim(cname).eq.'flag_compute_stress' ) then
    call read_l1(ionum,lstrs0)
    return
  elseif( trim(cname).eq.'cell_fix' ) then
    call read_ls(ionum,3,3,lcellfix)
    return
  elseif( trim(cname).eq.'deformation' ) then
    call read_c1(ionum,cdeform)
    return
  elseif( trim(cname).eq.'deform_ratio' ) then
    call read_rs(ionum,3,3,dhratio(1:3,1:3))
    return
  elseif( trim(cname).eq.'zload_type' ) then
    call read_c1(ionum,czload_type)
    return
  elseif( trim(cname).eq.'zload_skin_width' ) then
    call read_r1(ionum,zskin_width)
    return
  elseif( trim(cname).eq.'zload_shear_angle' ) then
    call read_r1(ionum,zshear_angle)
    return
  elseif( trim(cname).eq.'final_strain' ) then
    call read_r1(ionum,strfin)
    return
  elseif( trim(cname).eq.'mass' ) then
    backspace(ionum)
    read(ionum,*) ctmp,csp,tmp
    is = csp2isp(csp)
    if( is.gt.0 ) am(is) = tmp
    return
  elseif( trim(cname).eq.'charge' ) then
    backspace(ionum)
    read(ionum,*) ctmp,is,schg(is)
    return
  elseif( trim(cname).eq.'io_format' ) then
    call read_c1(ionum,ciofmt)
    return
  elseif( trim(cname).eq.'force_type' .or. &
       trim(cname).eq.'force_field' ) then
!        call read_c1(ionum,cforce)
    call read_force_field(ionum)
    return
  elseif( trim(cname).eq.'fix_charge' ) then
    call read_c1(ionum,chgfix)
    return
  elseif( trim(cname).eq.'charge_optimize' .or. &
       trim(cname).eq.'variable_charge' ) then
    call read_l1(ionum,lvc)
    return
  elseif( trim(cname).eq.'flag_coulomb' ) then
    call read_i1(ionum,ifcoulomb)
    return
  elseif( trim(cname).eq.'print_level' ) then
    call read_i1(ionum,iprint)
    return
  elseif( trim(cname).eq.'pka_atom') then
    call read_i1(ionum,iatom_pka)
    return
  elseif( trim(cname).eq.'pka_energy') then
    call read_r1(ionum,pka_energy)
    return
  elseif( trim(cname).eq.'pka_theta') then
    call read_r1(ionum,pka_theta)
    return
  elseif( trim(cname).eq.'pka_phi') then
    call read_r1(ionum,pka_phi)
    return
  elseif( trim(cname).eq.'metadynamics') then
    call read_l1(ionum,lmetaD)
    return
  elseif( trim(cname).eq.'constraints') then
    call read_l1(ionum,lconst)
    return
  elseif( trim(cname).eq.'reduced_force' .or. &
       trim(cname).eq.'reduce_force' ) then
    call read_l1(ionum,lrdcfrc)
    return
  elseif( trim(cname).eq.'structure_analysis') then
    backspace(ionum)
    read(ionum,*) ctmp, cstruct, istruct
    return
  elseif( trim(cname).eq.'structure_rcut') then
    call read_r1(ionum,rc_struct)
    return
  elseif( trim(cname).eq.'overlay') then
    call read_overlay(ionum)
    return
  elseif( trim(cname).eq.'overlay_type') then
    call read_c1(ionum,ol_type)
    return
  elseif( trim(cname).eq.'overlay_force') then
    call read_c1(ionum,ol_force)
    return
  elseif( trim(cname).eq.'flag_extfrc') then
    call read_l1(ionum,lextfrc)
    return
  elseif( trim(cname).eq.'spcs_extfrc') then
    call read_c1(ionum,cspc_extfrc)
    return
  elseif( trim(cname).eq.'extfrc') then
    backspace(ionum)
    read(ionum,*) ctmp,extfrc(1:3)
    return
!.....Color charge NEMD
  elseif( trim(cname).eq.'flag_clrchg') then
    call read_l1(ionum,lclrchg)
    return
  elseif( trim(cname).eq.'clr_init') then
    call read_c1(ionum,clr_init)
    return
  elseif( trim(cname).eq.'spcs_clrchg') then
    call read_c1(ionum,cspc_clrchg)
    return
  elseif( trim(cname).eq.'clrfield') then
    backspace(ionum)
    read(ionum,*) ctmp,clrfield(1:3)
    return
!.....Local flux
  elseif( trim(cname).eq.'flag_lflux') then
    call read_l1(ionum,lflux)
    return
  elseif( trim(cname).eq.'num_out_lflux') then
    call read_i1(ionum,noutlflux)
    return
  elseif( trim(cname).eq.'ndiv_lflux') then
    backspace(ionum)
    read(ionum,*) ctmp,nlx,nly,nlz
    return
!.....Probability density
  elseif( trim(cname).eq.'flag_pdens') then
    call read_l1(ionum,lpdens)
    return
  elseif( trim(cname).eq.'spcs_pdens') then
    call read_c1(ionum,cspc_pdens)
    return
  elseif( trim(cname).eq.'ndiv_pdens') then
    backspace(ionum)
    read(ionum,*) ctmp,npx,npy,npz
    return
  elseif( trim(cname).eq.'orig_pdens') then
    backspace(ionum)
    read(ionum,*) ctmp,orig_pdens(1:3)
    return
  elseif( trim(cname).eq.'hmat_pdens' ) then
    read(ionum,*) hmat_pdens(1:3,1)
    read(ionum,*) hmat_pdens(1:3,2)
    read(ionum,*) hmat_pdens(1:3,3)
    return
!.....Reallocation
  elseif( trim(cname).eq.'allow_reallocation') then
    call read_l1(ionum,lrealloc)
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
!      elseif( trim(cname).eq.'' ) then
!        call read_i1(ionum,nz)
!        return
  endif

!      write(6,'(a)') " [Error] No match: "//trim(cname)//" !!!"
!      stop
  write(6,'(a)') ' Warning: No such in.pmd entry, '//trim(cname)//' !!!'
  return

end subroutine set_variable
!=======================================================================
subroutine read_r1(ionum,rval)
!
!  Read one read*8 parameter from the line
!
  integer,intent(in):: ionum
  real(8),intent(out):: rval
  character(len=128):: ctmp

  backspace(ionum)
  read(ionum,*) ctmp,rval
!      write(6,'(1x,a,es15.3)') trim(ctmp),rval

end subroutine read_r1
!=======================================================================
subroutine read_rs(ionum,ndata,nrow,rval)
!
!  Read several read*8 parameters
!
  integer,intent(in):: ionum,ndata,nrow
  real(8),intent(out):: rval(ndata,nrow)

  do n=1,nrow
    read(ionum,*) (rval(i,n),i=1,ndata)
!        write(6,'(1x,100es15.3)') (rval(i,n),i=1,ndata)
  enddo

end subroutine read_rs
!=======================================================================
subroutine read_i1(ionum,ival)
!
!  Read one integer parameter from the line
!
  integer,intent(in):: ionum
  integer,intent(out):: ival
  character(len=128):: ctmp

  backspace(ionum)
  read(ionum,*) ctmp,ival
!      write(6,'(1x,a,i10)') trim(ctmp),ival

end subroutine read_i1
!=======================================================================
subroutine read_c1(ionum,cval)
!
!  Read one word from the line
!
  integer,intent(in):: ionum
  character(len=*),intent(out):: cval
  character(len=128):: ctmp

  backspace(ionum)
  read(ionum,*) ctmp,cval
!      write(6,'(1x,2a)') trim(ctmp),trim(cval)

end subroutine read_c1
!=======================================================================
subroutine read_l1(ionum,lval)
!
!  Read logical variable
!
  integer,intent(in):: ionum
  logical,intent(out):: lval
  character(len=128):: ctmp

  backspace(ionum)
  read(ionum,*) ctmp,lval
!      write(6,'(1x,a,5x,l1)') trim(ctmp),lval

end subroutine read_l1
!=======================================================================
subroutine read_ls(ionum,ndata,nrow,lval)
!
!  Read several logical parameters
!
  integer,intent(in):: ionum,ndata,nrow
  logical,intent(out):: lval(ndata,nrow)

  do n=1,nrow
    read(ionum,*) (lval(n,i),i=1,ndata)
!        write(6,'(1x,100es15.3)') (rval(i,n),i=1,ndata)
  enddo

end subroutine read_ls
!=======================================================================
subroutine read_force_field(ionum)
!
!     Read forces
!     There is no limit of number of force-fields to be specified.
!
!.....use pmdio
  use force, only: num_forces, force_list
  use util, only: num_data
  implicit none
  integer,intent(in):: ionum

  integer:: i,ndat
  character(len=1024):: ctmp
!!$  integer,external:: num_data

  backspace(ionum)
  read(ionum,'(a)') ctmp
  ndat = num_data(trim(ctmp),' ')
  if( ndat.lt.2 ) write(6,'(a)') 'There is no force-field' &
       //' specified.'
!      numff = ndat -1
  num_forces = ndat -1
!      allocate(cffs(numff))
!!$  allocate(force_list(num_forces))
  backspace(ionum)
!      read(ionum,*) ctmp, (cffs(i),i=1,numff)
  read(ionum,*) ctmp, (force_list(i),i=1,num_forces)
end subroutine read_force_field
!=======================================================================
subroutine read_overlay(ionum)
!
!  Read overlay of a given pair
!
  use pmdvars, only: specorder, nspmax
  use force, only: ol_ranges, loverlay
  use util, only: num_data, csp2isp
  implicit none 
  integer,intent(in):: ionum
  
  character(len=1024):: ctmp
  character(len=128):: ctmp1
  integer:: isp,ndat
  character(len=3):: cspi
  real(8):: rin, rout

  backspace(ionum)
  read(ionum,'(a)') ctmp
  ndat = num_data(trim(ctmp),' ')
  if( ndat.lt.4 ) stop 'ERROR: wrong format for overlay entry.'
  loverlay = .true.
  read(ctmp,*) ctmp1, cspi, rin, rout
  isp = csp2isp(cspi)
  if( isp.gt.0 ) then
    ol_ranges(1,isp) = rin
    ol_ranges(2,isp) = rout
  else
    print *,'Overlay for '//trim(cspi)//' is not set, '//&
         'because the specified species-pair is not found in the system.'
  endif
  
end subroutine read_overlay
!=======================================================================
subroutine read_dumpaux(ionum)
!
!  Read dump_aux_order entry
!
  use pmdvars, only: cdumpaux,ldumpaux_changed
  use util, only: num_data
  implicit none 
  integer,intent(in):: ionum
  
  character(len=1024):: ctmp
  character(len=20),allocatable:: ctmp1(:)
  integer:: ndat,i

  backspace(ionum)
  read(ionum,'(a)') ctmp
  ndat = num_data(trim(ctmp),' ')
  allocate(ctmp1(ndat))
  read(ctmp,*) (ctmp1(i),i=1,ndat)
!.....Create cdumpaux string without entry keyword, dump_aux_order
  write(cdumpaux,*) (trim(ctmp1(i))//' ',i=2,ndat)
!!$  ldumpaux_changed = .true.
  deallocate(ctmp1)
end subroutine read_dumpaux
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "make pmd lib"
!     End:
