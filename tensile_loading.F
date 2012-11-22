      subroutine tensile_loading(natm,ra,tag,nstp,strfin,strnow
     &     ,sorg,mpi_md_world)
c
c  To apply tensile loading to Al nanorod, move atoms of top and
c  bottom layers in opposite direction.
c  This routine is applicable only to non-parallel calculation.
c
      implicit none
      include 'mpif.h'
      integer,intent(in):: natm,nstp,mpi_md_world
      real(8),intent(in):: tag(natm),strfin,sorg(3)
      real(8),intent(inout):: ra(3,natm),strnow

      integer:: i,l,ierr
      real(8):: zmin,zmax,zmaxl,zminl
      logical,save:: l1st=.true.
      real(8),save:: zl0,dzlfin,dzl,zl

c.....At first call, determine velocity of top and bottom layers
      if( l1st ) then
c.....Get z-length
        zminl=1d0
        zmaxl=0d0
        do i=1,natm
          zminl= min(zminl,ra(3,i)+sorg(3))
          zmaxl= max(zmaxl,ra(3,i)+sorg(3))
        enddo
        zmax= 0d0
        zmin= 0d0
        call mpi_allreduce(zmaxl,zmax,1,mpi_double_precision,mpi_max
     &       ,mpi_md_world,ierr)
        call mpi_allreduce(zminl,zmin,1,mpi_double_precision,mpi_min
     &       ,mpi_md_world,ierr)
        zl0= zmax-zmin
        dzlfin= zl0 *strfin/100
        dzl= dzlfin/nstp /2
        zl= zl0
c        write(6,'(a,2es12.4)') ' zl0,dzl =',zl0,dzl
        l1st=.false.
      endif

      do i=1,natm
        l= int(mod(tag(i)*10,10d0))
        if( l.eq.2 .and. ra(3,i)+sorg(3).gt.0.5d0 ) then ! top
          ra(3,i)= ra(3,i) +dzl
        else if( l.eq.2 .and. ra(3,i)+sorg(3).le.0.5d0 ) then ! bottom
          ra(3,i)= ra(3,i) -dzl
        endif
      enddo
      zl= zl +2d0*dzl
      strnow= (zl-zl0)/zl0 *100d0
c      write(6,'(a,2es12.3e3)') ' zl,strnow=',zl,strnow
      
      end subroutine tensile_loading
c=======================================================================
      subroutine get_stress_on_base(natm,ra,aa,tag,h,acon,ftop,fbot
     &     ,sorg,mpi_md_world)
c
c  Compute forces on atoms at top and bottom bases,
c  and get stress on bases dividing by area.
c
      implicit none
      include 'mpif.h'
      integer,intent(in):: natm,mpi_md_world
      real(8),intent(in):: ra(3,natm),aa(3,natm),tag(natm),h(3,3,0:1)
     &     ,acon,sorg(3)
      real(8),intent(out):: ftop,fbot
      integer:: i,l,ierr
      real(8):: ftopl,fbotl

      real(8),parameter:: pi  = 3.14159265358979d0
      real(8),parameter:: rad = 1.0d-9 ! m
      real(8),parameter:: area= pi*rad*rad ! m^2

      ftopl= 0d0
      fbotl= 0d0
      do i=1,natm
        l= int(mod(tag(i)*10,10d0))
        if( l.eq.2 .and. ra(3,i)+sorg(3).gt.0.5d0 ) then ! top layer
          ftopl=ftopl +aa(3,i)
        else if( l.eq.2 .and. ra(3,i)+sorg(3).le.0.5d0 ) then ! bottom
          fbotl=fbotl +aa(3,i)
        endif
      enddo

      ftop= 0d0
      fbot= 0d0
      call mpi_allreduce(ftopl,ftop,1,mpi_double_precision
     &     ,mpi_sum,mpi_md_world,ierr)
      call mpi_allreduce(fbotl,fbot,1,mpi_double_precision
     &     ,mpi_sum,mpi_md_world,ierr)

c.....Convert unit to [N]
      ftop= ftop /acon *h(3,3,0) *82.387d-9
      fbot= fbot /acon *h(3,3,0) *82.387d-9
c.....Divide by area and get stress in [GPa]=[1d+9 Pa]=[1d+9 N/m^2]
      ftop= ftop /area *1d-9
      fbot= fbot /area *1d-9
      
      
      end subroutine get_stress_on_base
