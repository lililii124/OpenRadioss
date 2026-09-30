! SPDX-License-Identifier: AGPL-3.0-or-later
subroutine lecmuser01(iin,iout,uparam,maxu,nu,ns,ifunc,maxf,nf,parmat,userbuf)
  use law_user, only: ulawbuf
  use zwt_material_mod, only: wp,nparam,nstate,zwt_validate
  implicit none
  integer, intent(in) :: iin,iout,maxu,maxf
  integer, intent(out) :: nu,ns,nf,ifunc(maxf)
  real(wp), intent(out) :: uparam(maxu),parmat(100)
  type(ulawbuf), intent(inout) :: userbuf
  integer :: ios,status
  nu=0
  ns=0
  nf=0
  ifunc=0
  parmat=0.0_wp
  if (maxu<nparam) then
    write(iout,*) 'ZWT USER01: insufficient parameter storage'
    call arret(2)
    return
  end if
  read(iin,*,iostat=ios) uparam(1:nparam)
  status=1
  if (ios==0) call zwt_validate(uparam(1:nparam),status)
  if (ios/=0.or.status/=0) then
    write(iout,*) 'ZWT USER01: invalid or unstable parameter card, material ',userbuf%id,' status ',status
    call arret(2)
    return
  end if
  nu=nparam
  ns=nstate
  parmat(1)=uparam(1)
  parmat(2)=uparam(3)+uparam(6)+uparam(8)
  parmat(3)=uparam(10)
  parmat(16)=2.0_wp
  parmat(17)=2*uparam(2)/(uparam(1)+4*uparam(2)/3)
  write(iout,*) 'ZWT USER01: small-stretch ZWT; solids require Ismstr=10, Iframe=1'
end subroutine
