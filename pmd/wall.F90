!=======================================================================
!                     Last modified: <2016-10-30 07:26:30 Ryo KOBAYASHI>
!=======================================================================
! Since these wall functions are obsolete and I do not understand
! what it does anymore. It was separated from the main routine.
!=======================================================================
subroutine wall_reflect()
  use variables
  use wall
  implicit none
  include "mpif.h"
  integer:: i,is
  real(8):: zi

  do i=1,natm
    is= int(tag(i))
    if( is.ne.2 ) cycle
    zi= (ra(3,i)+sorg(3))
!.....if atom over wtop (top wall), go downward
    if( zi.ge.wtop ) then
      va(3,i)=-abs(va(3,i))
      va(3,i)= min(va(3,i),dwtop)
      ptop=ptop +am(is)*2d0*abs(va(3,i))*h(3,3,0)/dt
      nitop=nitop +1
    endif
!.....if atom under wbot (bottom wall), go upward
    if( zi.le.wbot ) then
      va(3,i)= abs(va(3,i))
      va(3,i)= max(va(3,i),dwbot)
      pbot=pbot +am(is)*2d0*abs(va(3,i))*h(3,3,0)/dt
      nibot=nibot +1
    endif
  enddo
end subroutine wall_reflect
!=======================================================================
subroutine update_wall()
!
!  Change of wall with relaxation time.
!  See RK's note on 2013.02.13
!
  use variables
  use wall
  implicit none
  include "mpif.h"
  integer:: i,itmp,is
  real(8):: zi,tmp,vtop,vbot,dvtop,dvbot
  logical,save:: l1st=.true.

  if(l1st) then
    if( wtop.lt.wbot ) stop ' [Error] wtop.lt.wbot !!!'
!.....btop,bbot: top/bottom of bulk material are fixed
    btop=0d0
    bbot=1d0
    do i=1,natm
      is= int(tag(i))
      if( is.ne.1 ) cycle
      zi= (ra(3,i)+sorg(3))
      btop= max(btop,zi)
      bbot= min(bbot,zi)
    enddo
    tmp= btop
    call mpi_allreduce(tmp,btop,1,mpi_double_precision,mpi_max &
         ,mpi_md_world,ierr)
    tmp= bbot
    call mpi_allreduce(tmp,bbot,1,mpi_double_precision,mpi_min &
         ,mpi_md_world,ierr)
!.....area of wall
    area_wall= h(2,2,0)*h(1,1,0)
    if( myid_md.eq.0 ) then
      write(6,'(a,es12.4)') ' btop =',btop
      write(6,'(a,es12.4)') ' bbot =',bbot
      write(6,'(a,es12.4)') ' area =',area_wall
    endif
!.....reset pressure
    ptop= 0d0
    pbot= 0d0
    nitop= 0
    nibot= 0
    l1st=.false.
!.....return at 1st call
    return
  endif

  itmp= nitop
  call mpi_reduce(itmp,nitop,1,mpi_integer,mpi_sum,0, &
       mpi_md_world,ierr)
  itmp= nibot
  call mpi_reduce(itmp,nibot,1,mpi_integer,mpi_sum,0, &
       mpi_md_world,ierr)
  tmp= ptop
  call mpi_reduce(tmp,ptop,1,mpi_double_precision,mpi_sum,0, &
       mpi_md_world,ierr)
  ptop= ptop /nodes_md
  tmp= pbot
  call mpi_reduce(tmp,pbot,1,mpi_double_precision,mpi_sum,0, &
       mpi_md_world,ierr)
  pbot= pbot /nodes_md

  if( myid_md.eq.0 ) then
    if( nitop.ne.0 ) then
      ptop= ptop /area_wall /dt_wall
!.....current volume
      vtop= abs(wtop-btop)*h(3,3,0) *area_wall
!.....change of volume
      dvtop= -vtop*(ptgt_wall -ptop)/ptop *dt_wall/trlx_wall
!.....speed of wall
      dwtop= abs(wtop-btop) *dvtop/vtop /nout_wall
!.....change of wall position
!          wtop= btop +abs(wtop-btop)*(1d0+dvtop/vtop)
!          wtop= min(wtop,0.99d0)
    endif
    if( nibot.ne.0 ) then
      pbot= pbot /area_wall /dt_wall
      vbot= abs(wbot-bbot)*h(3,3,0) *area_wall
      dvbot= -vbot*(ptgt_wall -pbot)/pbot *dt_wall/trlx_wall
      dwbot= -abs(wbot-bbot) *dvbot/vbot /nout_wall
!          wbot= bbot -abs(wbot-bbot)*(1d0+dvbot/vbot)
!          wbot= max(wbot,0.01d0)
    endif
  endif

!      call mpi_bcast(wtop,1,mpi_double_precision,0,mpi_md_world,ierr)
!      call mpi_bcast(wbot,1,mpi_double_precision,0,mpi_md_world,ierr)
  call mpi_bcast(dwtop,1,mpi_double_precision,0,mpi_md_world,ierr)
  call mpi_bcast(dwbot,1,mpi_double_precision,0,mpi_md_world,ierr)

end subroutine update_wall

