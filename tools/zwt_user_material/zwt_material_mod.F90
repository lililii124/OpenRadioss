! SPDX-License-Identifier: AGPL-3.0-or-later
! Classical scalar ZWT with the explicit small-stretch tensor extension in README.md.
module zwt_material_mod
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  integer, parameter, public :: wp=real64, nparam=11, nstate=18
  public :: zwt_update, zwt_validate, zwt_defaults, tensor_to_voigt, voigt_to_tensor
contains
  pure function voigt_to_tensor(v) result(a)
    real(wp), intent(in) :: v(6)
    real(wp) :: a(3,3)
    a(1,:)=[v(1),v(4),v(6)]
    a(2,:)=[v(4),v(2),v(5)]
    a(3,:)=[v(6),v(5),v(3)]
  end function

  pure function tensor_to_voigt(a) result(v)
    real(wp), intent(in) :: a(3,3)
    real(wp) :: v(6)
    v=[a(1,1),a(2,2),a(3,3),a(1,2),a(2,3),a(1,3)]
  end function

  pure function det3(a) result(d)
    real(wp), intent(in) :: a(3,3)
    real(wp) :: d
    d=a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) &
     -a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)) &
     +a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
  end function

  pure function inverse3(a) result(b)
    real(wp), intent(in) :: a(3,3)
    real(wp) :: b(3,3)
    b(1,:)=[a(2,2)*a(3,3)-a(2,3)*a(3,2), &
            a(1,3)*a(3,2)-a(1,2)*a(3,3),a(1,2)*a(2,3)-a(1,3)*a(2,2)]
    b(2,:)=[a(2,3)*a(3,1)-a(2,1)*a(3,3), &
            a(1,1)*a(3,3)-a(1,3)*a(3,1),a(1,3)*a(2,1)-a(1,1)*a(2,3)]
    b(3,:)=[a(2,1)*a(3,2)-a(2,2)*a(3,1), &
            a(1,2)*a(3,1)-a(1,1)*a(3,2),a(1,1)*a(2,2)-a(1,2)*a(2,1)]
    b=b/det3(a)
  end function

  pure function cunit(e,nu) result(s)
    real(wp), intent(in) :: e(3,3),nu
    real(wp) :: s(3,3),p
    integer :: j
    s=e/(1.0_wp+nu)
    p=nu*(e(1,1)+e(2,2)+e(3,3))/((1.0_wp+nu)*(1.0_wp-2.0_wp*nu))
    do j=1,3
      s(j,j)=s(j,j)+p
    end do
  end function

  subroutine zwt_defaults(cm)
    real(wp), intent(out) :: cm(nparam)
    ! Synthetic SI verification parameters; not an experimental calibration.
    cm=[2.2e9_wp/3,1.1e9_wp,1.0e9_wp,-5.0e9_wp,2.0e10_wp, &
        4.0e8_wp,1.0e-2_wp,8.0e8_wp,1.0e-4_wp,0.0_wp,0.05_wp]
  end subroutine

  subroutine zwt_validate(cm,status)
    real(wp), intent(in) :: cm(nparam)
    integer, intent(out) :: status
    real(wp) :: nu,scale,young,cmin,cmax,a,b,limit,lower,vertex,margin
    status=1
    if (.not.all(ieee_is_finite(cm))) return
    status=2
    nu=cm(10)
    limit=cm(11)
    if (minval(cm(1:3))<=0.0_wp.or.cm(6)<0.0_wp.or.cm(8)<0.0_wp) return
    if (cm(7)<=0.0_wp.or.cm(9)<=0.0_wp.or.nu<=-1.0_wp.or.nu>=0.5_wp) return
    if (limit<=0.0_wp.or.limit>0.05_wp) return
    scale=maxval(abs(cm([1,2,3,4,5,6,8])))
    cmin=min(1.0_wp/(1.0_wp+nu),1.0_wp/(1.0_wp-2.0_wp*nu))
    cmax=max(1.0_wp/(1.0_wp+nu),1.0_wp/(1.0_wp-2.0_wp*nu))
    ! Reserve arithmetic headroom in both the stress and acoustic calculations.
    if (scale>huge(scale)/(1024.0_wp*max(1.0_wp,cmax))) return
    young=cm(3)/scale+cm(6)/scale+cm(8)/scale
    status=3
    if (abs(cm(1)/scale-young/(3*(1-2*nu)))>1.0e-8_wp*cm(1)/scale) return
    if (abs(cm(2)/scale-young/(2*(1+nu)))>1.0e-8_wp*cm(2)/scale) return
    ! A sufficient positive-definiteness certificate throughout ||e||F<=limit.
    ! Spectral divided differences of a*e^2+b*e^3 are bounded by min g'(x).
    a=cm(4)/scale
    b=cm(5)/scale
    lower=min(-2*a*limit+3*b*limit**2,2*a*limit+3*b*limit**2)
    if (b>0.0_wp) then
      if (abs(a)<3*b*limit) then
        vertex=-a/(3*b)
        lower=min(lower,2*a*vertex+3*b*vertex**2)
      end if
    end if
    margin=128*epsilon(scale)*max(cm(3)/scale*cmin,abs(lower))
    status=9
    if (cm(3)/scale*cmin+lower<=margin) return
    status=0
  end subroutine

  subroutine relaxation(dt,theta,decay,weight)
    real(wp), intent(in) :: dt,theta
    real(wp), intent(out) :: decay,weight
    real(wp) :: x
    ! Avoid overflow in dt/theta, including otherwise valid extreme time units.
    if (theta<dt/50.0_wp) then
      decay=0.0_wp
      weight=theta/dt
    else
      x=dt/theta
      decay=exp(-x)
      if (x<1.0e-4_wp) then
        weight=1+x*(-0.5_wp+x*(1.0_wp/6+x*(-1.0_wp/24+x/120)))
      else
        weight=(1-decay)/x
      end if
    end if
  end subroutine

  subroutine zwt_update(cm,dt,f,state,stress,acoustic,status)
    real(wp), intent(in) :: cm(nparam),dt,f(3,3)
    real(wp), intent(inout) :: state(nstate),stress(6),acoustic
    integer, intent(out) :: status
    real(wp) :: r(3,3),rn(3,3),u(3,3),e(3,3),de(3,3),eold(3,3)
    real(wp) :: total(3,3),e2(3,3),q(3,3),next(nstate),sig(6),modulus
    real(wp) :: limit,nu,err,jac,decay,weight,young,theta,cmax,bound,scale
    integer :: k,it,offset
    call zwt_validate(cm,status)
    if (status/=0) return
    status=1
    if (.not.all(ieee_is_finite(f)).or..not.all(ieee_is_finite(state))) return
    if (.not.ieee_is_finite(dt)) return
    status=2
    if (dt<0.0_wp) return
    nu=cm(10)
    limit=cm(11)
    cmax=max(1.0_wp/(1+nu),1.0_wp/(1-2*nu))
    status=10
    ! Reject invalid restart states before squaring or multiplying their values.
    if (maxval(abs(state(1:6)))>limit*(1+1.0e-8_wp)) return
    eold=voigt_to_tensor(state(1:6))
    if (sqrt(sum(eold*eold))>limit*(1+1.0e-8_wp)) return
    do k=1,2
      offset=6*k
      young=cm(2*k+4)
      if (young==0.0_wp) then
        if (any(state(offset+1:offset+6)/=0.0_wp)) return
      else
        ! q=Ei*C:(e-exponential_average(e)), for histories inside the strain ball.
        bound=2*young*cmax*limit
        if (maxval(abs(state(offset+1:offset+6)))>bound*(1+1.0e-8_wp)) return
      end if
    end do
    status=4
    ! Any admissible F=R*(I+e) has |Fij|<=1+limit; early bound prevents overflow.
    if (maxval(abs(f))>2.0_wp) return
    jac=det3(f)
    if (jac<(1-limit)**3*(1-1.0e-8_wp).or.jac>(1+limit)**3*(1+1.0e-8_wp)) return
    r=f
    do it=1,30
      rn=0.5_wp*(r+transpose(inverse3(r)))
      err=maxval(abs(rn-r))
      r=rn
      if (err<32*epsilon(err)) exit
    end do
    status=5
    if (it>30.or..not.all(ieee_is_finite(r))) return
    u=matmul(transpose(r),f)
    e=0.5_wp*(u+transpose(u))
    do k=1,3
      e(k,k)=e(k,k)-1
    end do
    status=6
    if (sqrt(sum(e*e))>limit*(1+1.0e-8_wp)) return
    de=e-eold
    status=7
    if (dt==0.0_wp) then
      if (maxval(abs(de))>256*epsilon(dt)) return
      ! Initialization/restart probes must never drift stored strain or memory.
      e=eold
      de=0.0_wp
    end if
    e2=matmul(e,e)
    total=cm(3)*cunit(e,nu)+cm(4)*e2+cm(5)*matmul(e2,e)
    next=state
    next(1:6)=tensor_to_voigt(e)
    do k=1,2
      offset=6*k
      young=cm(2*k+4)
      theta=cm(2*k+5)
      call relaxation(dt,theta,decay,weight)
      q=decay*voigt_to_tensor(state(offset+1:offset+6))+young*weight*cunit(de,nu)
      next(offset+1:offset+6)=tensor_to_voigt(q)
      total=total+q
    end do
    sig=tensor_to_voigt(matmul(matmul(r,total),transpose(r)))
    ! Upper bound for the small-strain tangent over the entire admitted strain ball.
    scale=maxval(abs(cm([1,2,3,4,5,6,8])))
    modulus=scale*((cm(3)/scale+cm(6)/scale+cm(8)/scale)*cmax &
                  +2*abs(cm(4)/scale)*limit+3*abs(cm(5)/scale)*limit**2)
    status=8
    if (.not.all(ieee_is_finite(sig)).or..not.all(ieee_is_finite(next))) return
    if (.not.ieee_is_finite(modulus).or.modulus<=0.0_wp) return
    if (dt>0.0_wp) state=next
    stress=sig
    acoustic=modulus
    status=0
  end subroutine
end module
