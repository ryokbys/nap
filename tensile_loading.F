      subroutine tensile_loading(natm,ra,tag,nstp,strfin,strnow)
c
c  To apply tensile loading to Al nanorod, move atoms of top and
c  bottom layers in opposite direction.
c  This routine is applicable only to non-parallel calculation.
c
      implicit none
      integer,intent(in):: natm,nstp
      real(8),intent(in):: tag(natm),strfin
      real(8),intent(inout):: ra(3,natm),strnow

      integer:: i,is
      real(8):: zmin,zmax
      logical,save:: l1st=.true.
      real(8),save:: zl0,dzlfin,dzl,zl

c.....At first call, determine velocity of top and bottom layers
      if( l1st ) then
c.....Get z-length
        zmin=1d0
        zmax=0d0
        do i=1,natm
          zmin= min(zmin,ra(3,i))
          zmax= max(zmax,ra(3,i))
        enddo
        zl0= zmax-zmin
        dzlfin= zl0 *strfin/100
        dzl= dzlfin/nstp /2
        zl= zl0
        l1st=.false.
      endif

      do i=1,natm
        is= int(tag(i))
        if( is.eq.2 ) then       ! top layer
          ra(3,i)= ra(3,i) +dzl
        else if( is.eq.3 ) then  ! bottom layer
          ra(3,i)= ra(3,i) -dzl
        endif
      enddo
      zl= zl +2d0*dzl
      strnow= (zl-zl0)/zl0 *100d0
      
      end subroutine tensile_loading
c=======================================================================
      subroutine get_stress_on_base(natm,aa,tag,h,acon,ftop,fbot)
c
c  Compute forces on atoms at top and bottom bases,
c  and get stress on bases dividing by area.
c
      implicit none
      integer,intent(in):: natm
      real(8),intent(in):: aa(3,natm),tag(natm),h(3,3,0:1),acon
      real(8),intent(out):: ftop,fbot
      integer:: i,is

      real(8),parameter:: pi  = 3.14159265358979d0
      real(8),parameter:: rad = 2.0d-9 ! m
      real(8),parameter:: area= pi*rad*rad ! m^2

      ftop= 0d0
      fbot= 0d0
      do i=1,natm
        is= int(tag(i))
        if( is.eq.2 ) then      ! top layer
          ftop=ftop +aa(3,i)
        else if( is.eq.3 ) then ! bottom layer
          fbot=fbot +aa(3,i)
        endif
      enddo

c.....Convert unit to [N]
      ftop= ftop /acon *h(3,3,0) *82.387d-9
      fbot= fbot /acon *h(3,3,0) *82.387d-9
c.....Divide by area and get stress in [GPa]=[1d+9 Pa]=[1d+9 N/m^2]
      ftop= ftop /area *1d-9
      fbot= fbot /area *1d-9
      
      end subroutine get_stress_on_base
