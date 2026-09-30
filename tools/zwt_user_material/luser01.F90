! SPDX-License-Identifier: AGPL-3.0-or-later
subroutine luser01(nel,nu,ns,nf,ifunc,npf,tf,time,dt,uparam,rho,volume,eint, &
                   ngl,sound,visc,uvar,off,sigy,pla,userbuf)
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use law_userso, only: ulawintbuf
  use zwt_material_mod, only: wp,nparam,nstate,zwt_update
  implicit none
  integer, intent(in) :: nel,nu,ns,nf,ifunc(nf),npf(*),ngl(nel)
  real(wp), intent(in) :: tf(*),time,dt,uparam(nu),rho(nel),volume(nel),eint(nel)
  real(wp), intent(out) :: sound(nel),visc(nel),sigy(nel),pla(nel)
  real(wp), intent(inout) :: uvar(nel,ns),off(nel)
  type(ulawintbuf), intent(inout) :: userbuf
  real(wp) :: f(3,3),sig(6),state(nstate),acoustic
  integer :: i,status,nchar
  character(len=160) :: message
  if (nu/=nparam.or.ns/=nstate) then
    message='ZWT USER01: incompatible material or restart state size'
    nchar=len_trim(message)
    call write_iout(message,nchar)
    call arret(2)
    return
  end if
  do i=1,nel
    sound(i)=0.0_wp
    visc(i)=0.0_wp
    sigy(i)=0.0_wp
    pla(i)=0.0_wp
    userbuf%dpla(i)=0.0_wp
    sig=0.0_wp
    if (off(i)>0.0_wp) then
      status=11
      if (ieee_is_finite(rho(i)).and.rho(i)>0.0_wp) then
        f(1,:)=[userbuf%fpsxx(i),userbuf%fpsxy(i),userbuf%fpsxz(i)]
        f(2,:)=[userbuf%fpsyx(i),userbuf%fpsyy(i),userbuf%fpsyz(i)]
        f(3,:)=[userbuf%fpszx(i),userbuf%fpszy(i),userbuf%fpszz(i)]
        state=uvar(i,1:nstate)
        acoustic=0.0_wp
        call zwt_update(uparam(1:nparam),dt,f,state,sig,acoustic,status)
        if (status==0) then
          if (rho(i)<1.0_wp) then
            if (acoustic>(0.5_wp*huge(acoustic))*rho(i)) status=11
          end if
        end if
        if (status==0) then
          ! Separate square roots avoid overflow of acoustic/rho.
          sound(i)=sqrt(acoustic)/sqrt(rho(i))
          if (.not.ieee_is_finite(sound(i))) status=11
        end if
      end if
      if (status/=0) then
        write(message,'(A,I12,A,I4)') 'ZWT USER01: rejected state, element ',ngl(i),' status ',status
        nchar=len_trim(message)
        call write_iout(message,nchar)
        call arret(2)
        return
      end if
      uvar(i,1:nstate)=state
    end if
    userbuf%signxx(i)=sig(1)
    userbuf%signyy(i)=sig(2)
    userbuf%signzz(i)=sig(3)
    userbuf%signxy(i)=sig(4)
    userbuf%signyz(i)=sig(5)
    userbuf%signzx(i)=sig(6)
    userbuf%sigvxx(i)=0.0_wp
    userbuf%sigvyy(i)=0.0_wp
    userbuf%sigvzz(i)=0.0_wp
    userbuf%sigvxy(i)=0.0_wp
    userbuf%sigvyz(i)=0.0_wp
    userbuf%sigvzx(i)=0.0_wp
  end do
end subroutine
