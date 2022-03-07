module distfunc
!-----------------------------------------------------------------------
! Module for distribution functions such as RDF and ADF.
!-----------------------------------------------------------------------
  implicit none
  include './const.h'
  private
  save

  public:: calc_rdf, calc_adf

  real(8),parameter:: pi = 3.14159265358979d0
  
contains
!=======================================================================
  subroutine calc_rdf(natm,nnmax,tag,h,rmax,rmin,lspr,d2lspr, &
       iprint,l1st,lpairwise,msp,nbins,dists,rdfs)
!
!  Calculate RDF.
!
    implicit none
    integer,intent(in):: natm,nnmax,iprint,msp,nbins
    real(8),intent(in):: rmax,rmin,tag(natm),h(3,3)
    integer,intent(in):: lspr(0:nnmax,natm)
    real(8),intent(in):: d2lspr(nnmax,natm)
    logical,intent(in):: l1st,lpairwise
    real(8),intent(out):: dists(nbins)
    real(8),intent(out):: rdfs(nbins,0:msp,0:msp)

    integer:: ia,is,ja,js,jj,ib,ni,nj
    integer:: natms(msp)
    real(8):: dr,rc2,dij2,dij,rrdr,vol,tmp,r

    dr = (rmax-rmin)/nbins
    do ib=1,nbins
      dists(ib) = (rmin -dr/2) +dr*ib
    enddo
    
    rc2 = rmax*rmax
    rdfs(:,:,:) = 0d0
    do ia=1,natm
      is = int(tag(ia))
      do jj=1,lspr(0,ia)
        ja = lspr(jj,ia)
        js = int(tag(ja))
        dij2 = d2lspr(jj,ia)
        if( dij2.gt.rc2 ) cycle
        dij = dsqrt(dij2)
        rrdr = (dij-rmin)/dr
        if( rrdr.lt.0 ) cycle
        ib = min(int(rrdr)+1,nbins)
        rdfs(ib,js,is) = rdfs(ib,js,is) + 1d0
!!$      if( js.ne.is ) rdfs(js,is,ib) = rdfs(js,is,ib) + 1d0
      enddo
    enddo
    do is=1,msp
      do js=1,msp
        rdfs(:,0,0) = rdfs(:,0,0) +rdfs(:,js,is)
      enddo
    enddo

    vol = h(1,1)*(h(2,2)*h(3,3)-h(3,2)*h(2,3)) &
         +h(2,1)*(h(3,2)*h(1,3)-h(1,2)*h(3,3)) &
         +h(3,1)*(h(1,2)*h(2,3)-h(2,2)*h(1,3))

    if( lpairwise ) then
      natms(:) = 0
      do ia=1,natm
        is = int(tag(ia))
        natms(is) = natms(is) +1
      enddo
      tmp = 4d0 *pi *natm *(natm -1)/vol *dr
      do ib=1,nbins
        r = dists(ib)
        rdfs(ib,0,0) = rdfs(ib,0,0)/ (tmp*r*r)
      enddo
      do is=1,msp
        ni = natms(is)
        if( ni.eq.0 ) cycle
        do js=is,msp
          nj = natms(js)
          if( nj.eq.0 ) cycle
          tmp = 4d0*pi*dr /vol
          if( is.eq.js ) then
            if( ni.eq.1 ) cycle
            tmp = tmp *ni*(ni-1)
          else
            tmp = tmp *ni*nj
          endif
          do ib=1,nbins
            r = dists(ib)
            rdfs(ib,js,is) = rdfs(ib,js,is) /(tmp*r*r)
          enddo
        enddo
      enddo
    else
      tmp = 4d0 *pi *natm**2 /vol *dr
      do ib=1,nbins
        r = dists(ib)
        rdfs(ib,0,0) = rdfs(ib,0,0)/ (tmp*r*r)
      enddo
    endif
    
    return

  end subroutine calc_rdf
!=======================================================================
  subroutine calc_adf(natm,nnmax,tag,h,ra,rc,lspr,d2lspr,&
       ntrpl,itriples,dang,nang,angs,adfs)
!
!  Calculate the angular distribution function (ADF).
!
    implicit none
    integer,intent(in):: natm,nnmax,ntrpl,nang
    real(8),intent(in):: rc,tag(natm),h(3,3),ra(3,natm),dang
    integer,intent(in):: lspr(0:nnmax,natm),itriples(3,ntrpl)
    real(8),intent(in):: d2lspr(nnmax,natm)
    real(8),intent(out):: angs(nang)
    real(8),intent(out):: adfs(nang,ntrpl)

    integer:: ia,is,ja,js,jj,kk,ka,ks,itrpl,ijktrpl,iang
    real(8):: rc2,dij2,dik2,xi(3),xij(3),rij(3),xik(3),rik(3),dot,cs,rad,deg
    logical:: iexist,ijkexist

    do ia=1,nang
      angs(ia) = dang*(ia-1) +dang/2
    enddo
    
    rc2 = rc*rc
    adfs(:,:) = 0d0
    do ia=1,natm
      is = int(tag(ia))
      iexist = .false.
      do itrpl=1,ntrpl
        if( is.eq.itriples(1,itrpl) ) then
          iexist = .true.
          exit
        endif
      enddo
      if( .not. iexist ) cycle
      xi(1:3) = ra(1:3,ia)
      do jj=1,lspr(0,ia)
        ja = lspr(jj,ia)
        js = int(tag(ja))
        dij2 = d2lspr(jj,ia)
        if( dij2.gt.rc2 ) cycle
        xij(1:3) = ra(1:3,ja) -xi(1:3) -anint(ra(1:3,ja) -xi(1:3))
        rij(1:3) = h(1:3,1)*xij(1) +h(1:3,2)*xij(2) +h(1:3,3)*xij(3)
        do kk=jj+1,lspr(0,ia)
          ka = lspr(kk,ia)
          ks = int(tag(ka))
          dik2 = d2lspr(kk,ia)
          if( dik2.gt.rc2 ) cycle
          xik(1:3) = ra(1:3,ka) -xi(1:3) -anint(ra(1:3,ka) -xi(1:3))
          rik(1:3) = h(1:3,1)*xik(1) +h(1:3,2)*xik(2) +h(1:3,3)*xik(3)
          ijkexist = .false.
          do itrpl=1,ntrpl
            if( is.eq.itriples(1,itrpl) .and. &
                 ( (js.eq.itriples(2,itrpl).and.ks.eq.itriples(3,itrpl)) .or. &
                 (js.eq.itriples(3,itrpl).and.ks.eq.itriples(2,itrpl)) ) ) then
              ijkexist = .true.
              ijktrpl = itrpl
              exit
            endif
          enddo
          if( .not. ijkexist ) cycle
          dot = rij(1)*rik(1) +rij(2)*rik(2) +rij(3)*rik(3)
          cs = dot/sqrt(dij2)/sqrt(dik2)
          rad = acos(cs)
          deg = rad/pi *180d0
          iang = min(int(deg/dang)+1,nang)
          adfs(iang,ijktrpl) = adfs(iang,ijktrpl) +1d0
        enddo
      enddo
    enddo
    
    return
  end subroutine calc_adf
!=======================================================================
end module distfunc
!-----------------------------------------------------------------------
!     Local Variables:
!     compile-command: "gfortran -c mod_distfunc.F90"
!     End:
