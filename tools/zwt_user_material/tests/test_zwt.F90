! SPDX-License-Identifier: AGPL-3.0-or-later
program test_zwt
  use zwt_material_mod
  use law_userso, only: ulawintbuf
  use law_user, only: ulawbuf
  implicit none
  integer, parameter :: nel=4
  type(ulawintbuf) :: buf
  type(ulawbuf) :: stbuf
  real(wp) :: cm(nparam),readback(nparam),parmat(100),rho(nel),vol(nel),energy(nel)
  real(wp) :: sound(nel),visc(nel),history(nel,nstate),off(nel),sigy(nel),pla(nel)
  real(wp) :: expected(6,nel),hs(nstate,nel),acoustic(nel),f(3,3),r(3,3),u(3,3),v(6),tf(1)
  integer :: i,j,status,ngl(nel),ifunc(1),npf(1),nu,ns,nf,iu
  character(len=32) :: mode
  call check_core
  call get_command_argument(1,mode)
  call zwt_defaults(cm)
  if (trim(mode)=='invalid_card') then
    cm(4)=-2.0e10_wp
    cm(5)=0.0_wp
  end if
  stbuf%id=123
  open(newunit=iu,status='scratch',form='formatted')
  write(iu,*) cm
  rewind(iu)
  call lecmuser01(iu,6,readback,nparam,nu,ns,ifunc,1,nf,parmat,stbuf)
  close(iu)
  if (trim(mode)=='invalid_card') stop 90
  if (nu/=nparam.or.ns/=nstate.or.nf/=0) stop 1
  if (maxval(abs(readback-cm))>1.0e-10_wp) stop 2
  if (abs(parmat(1)/cm(1)-1)>1.0e-12_wp) stop 3
  if (abs(parmat(2)/(2*(1+parmat(3)))/cm(2)-1)>1.0e-12_wp) stop 4
  rho=1200.0_wp
  vol=1.0_wp
  energy=0.0_wp
  off=[1.0_wp,0.0_wp,1.0_wp,1.0_wp]
  history=0.0_wp
  history(2,:)=99.0_wp
  expected=0.0_wp
  hs=transpose(history)
  acoustic=0.0_wp
  tf=0.0_wp
  npf=0
  ifunc=0
  do i=1,nel
    ngl(i)=10*i
    u=voigt_to_tensor([0.001_wp*i,-0.0002_wp*i,0.0003_wp*i,0.0004_wp*i,-0.0005_wp*i,0.0006_wp*i])
    do j=1,3
      u(j,j)=u(j,j)+1
    end do
    r(1,:)=[cos(0.2_wp*i),-sin(0.2_wp*i),0.0_wp]
    r(2,:)=[sin(0.2_wp*i),cos(0.2_wp*i),0.0_wp]
    r(3,:)=[0.0_wp,0.0_wp,1.0_wp]
    f=matmul(r,u)
    if (off(i)>0) then
      call zwt_update(cm,0.001_wp,f,hs(:,i),expected(:,i),acoustic(i),status)
      if (status/=0) stop 5
    end if
    buf%fpsxx(i)=f(1,1)
    buf%fpsxy(i)=f(1,2)
    buf%fpsxz(i)=f(1,3)
    buf%fpsyx(i)=f(2,1)
    buf%fpsyy(i)=f(2,2)
    buf%fpsyz(i)=f(2,3)
    buf%fpszx(i)=f(3,1)
    buf%fpszy(i)=f(3,2)
    buf%fpszz(i)=f(3,3)
  end do
  if (trim(mode)=='missing_f') then
    buf%fpsxx(1)=0
    buf%fpsyy(1)=0
    buf%fpszz(1)=0
  end if
  if (trim(mode)=='bad_density') rho(1)=-1
  if (trim(mode)=='bad_state') history(1,1)=1
  if (trim(mode)=='wrong_size') ns=ns-1
  call luser01(nel,nparam,ns,0,ifunc,npf,tf,0.001_wp,0.001_wp,cm,rho,vol,energy, &
               ngl,sound,visc,history,off,sigy,pla,buf)
  if (len_trim(mode)>0) stop 91
  do i=1,nel
    v=[buf%signxx(i),buf%signyy(i),buf%signzz(i),buf%signxy(i),buf%signyz(i),buf%signzx(i)]
    if (maxval(abs(v-expected(:,i)))>1.0e-7_wp) stop 6
    if (maxval(abs(history(i,:)-hs(:,i)))>1.0e-7_wp) stop 7
    if (abs(sound(i)-sqrt(acoustic(i)/rho(i)))>1.0e-10_wp) stop 8
    if (pla(i)/=0.or.buf%dpla(i)/=0.or.sigy(i)/=0.or.visc(i)/=0) stop 9
    if (buf%sigvxx(i)/=0.or.buf%sigvxy(i)/=0) stop 10
  end do
  if (off(2)/=0.or.any(off([1,3,4])/=1)) stop 11
  print *, 'PASS: Starter, PARMAT, batch/inactive points, F ordering, rotation, sound, zero plastic outputs'
contains
  subroutine check_core
    real(wp) :: c(nparam),f0(3,3),h(nstate),saved(nstate),s(6),ac,saved_s(6),eq,q1,q2,dt0,e
    integer :: code,j
    call zwt_defaults(c)
    f0=0.0_wp
    do j=1,3
      f0(j,j)=1.0_wp
    end do
    e=0.01_wp
    f0(1,1)=1+e
    dt0=0.002_wp
    h=0.0_wp
    s=0.0_wp
    ac=0.0_wp
    eq=c(3)*e+c(4)*e**2+c(5)*e**3
    q1=c(6)*c(7)*(e/dt0)*(1-exp(-dt0/c(7)))
    q2=c(8)*c(9)*(e/dt0)*(1-exp(-dt0/c(9)))
    call zwt_update(c,dt0,f0,h,s,ac,code)
    if (code/=0.or.abs(s(1)/(eq+q1+q2)-1)>1.0e-12_wp) stop 20
    if (maxval(abs(s(2:6)))>1.0e-6_wp) stop 21
    call zwt_update(c,0.003_wp,f0,h,s,ac,code)
    if (code/=0.or.abs(s(1)/(eq+q1*exp(-.003_wp/c(7))+q2*exp(-.003_wp/c(9)))-1)>1.0e-12_wp) stop 22
    saved=h
    saved_s=s
    call zwt_update(c,0.0_wp,f0,h,s,ac,code)
    if (code/=0.or.any(h/=saved).or.any(s/=saved_s)) stop 23
    c(7)=1.0e-300_wp
    c(9)=1.0e-300_wp
    h=0.0_wp
    call zwt_update(c,1.0e300_wp,f0,h,s,ac,code)
    if (code/=0.or.abs(s(1)/eq-1)>1.0e-12_wp) stop 24
    c(7)=1.0e300_wp
    c(9)=1.0e300_wp
    h=0.0_wp
    call zwt_update(c,1.0e-300_wp,f0,h,s,ac,code)
    if (code/=0.or.abs(s(1)/(eq+(c(6)+c(8))*e)-1)>1.0e-12_wp) stop 25
    saved=h
    saved_s=s
    c(4)=-2.0e10_wp
    c(5)=0.0_wp
    call zwt_update(c,0.001_wp,f0,h,s,ac,code)
    if (code/=9.or.any(h/=saved).or.any(s/=saved_s)) stop 26
    call zwt_defaults(c)
    f0=1.0e300_wp
    call zwt_update(c,0.001_wp,f0,h,s,ac,code)
    if (code/=4.or.any(h/=saved).or.any(s/=saved_s)) stop 27
    c(5)=1.0e308_wp
    call zwt_validate(c,code)
    if (code==0) stop 28
    print *, 'PASS: ramp, relaxation, zero time, extreme ratios, stiffness certificate, transactional rejection'
  end subroutine
end program

subroutine arret(code)
  implicit none
  integer, intent(in) :: code
  print *, 'ARRET ',code
  stop 99
end subroutine

subroutine write_iout(message,nchar)
  implicit none
  integer, intent(in) :: nchar
  character(len=*), intent(in) :: message
  print *, message(1:nchar)
end subroutine
