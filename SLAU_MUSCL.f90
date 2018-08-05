module Mesh_Cell_PhysParam
  implicit none
  !!grid,cell parameter!!
  integer :: nxi !!cell number, not a grid number
  double precision, allocatable :: xgrid(:),xcell(:)

  !!physical parameter used in the governing equations
  double precision, allocatable :: rho(:),ux(:),pst(:),e_int(:)
  double precision, allocatable :: gamm(:),asou(:),tst(:)

  !!evaluation terms
  double precision, allocatable :: convx(:,:),dU(:,:),Q(:,:)

  !!time step parameter
  double precision :: dt
  integer :: nstep


  !!keizoku keisan option
  integer :: keizoku,timing,fileNo

contains

subroutine AllAllocate(nxi,xgrid,xcell,rho,ux,pst,e_int,gamm,asou,tst)
    implicit none
    integer :: nxi
    double precision, allocatable :: xgrid(:),xcell(:),rho(:),ux(:),pst(:)
    double precision, allocatable :: e_int(:),gamm(:),asou(:),tst(:)

    allocate(xgrid(nxi+1))
    allocate(xcell(nxi),rho(nxi),ux(nxi),pst(nxi))
    allocate(e_int(nxi),gamm(nxi),asou(nxi),tst(nxi))
end subroutine AllAllocate

subroutine TermsAllAllocate(nxi,convx,dU)
    implicit none
    integer :: nxi
    double precision, allocatable :: convx(:,:),dU(:,:)

    allocate(convx(2:nxi,3),dU(2:nxi-1,3),Q(2:nxi-1,3))
    
end subroutine TermsAllAllocate

subroutine ReadCond(nstep,timing,fileNo,dt,nxi,keizoku)
  implicit none
  integer, intent(out) :: nstep,timing,fileNo,nxi,keizoku
  double precision, intent(out) :: dt

  open(601,file='./cond.dat',status='old')
  read(601,*)nstep
  read(601,*)timing
  read(601,*)fileNo
  read(601,*)dt
  read(601,*)nxi
  read(601,*)keizoku
  close(601)

  print*, "nstep:",nstep,"timing:",timing
     
end subroutine ReadCond

subroutine ReadGrid(nxi,xgrid,xcell)
  implicit none
  integer, intent(in) :: nxi
  double precision, intent(inout) :: xgrid(nxi+1),xcell(nxi)
  integer i

  open(602,file='./grid.dat',status='old')

  do i=1,nxi+1
    read(602,*)xgrid(i)
  end do

  do i=1,nxi
    xcell(i)=0.5d0*(xgrid(i)+xgrid(i+1))
  end do

  close(602)
    
end subroutine ReadGrid

subroutine Initialize(nxi,keizoku,rho,ux,pst,e_int,gamm,asou,tst)
    implicit none
    integer, intent(in) :: nxi,keizoku
    double precision, intent(out) :: rho(nxi),ux(nxi),pst(nxi),e_int(nxi)
    double precision, intent(out) :: gamm(nxi),asou(nxi),tst(nxi)
    double precision :: R
    integer :: i

    R=208.2d0

    if (keizoku.eq.0) then
      do i=1,nxi/2
        gamm(i)=1.4d0
        rho(i)=1.0d0!1.0d0
        pst(i)=1.5d0!1.0d-1
        ux(i)=0.0d0
        ! tst(i)=pst(i)/(rho(i)*R)
        tst(i)=0.0d0
        asou(i)=sqrt(gamm(i)*pst(i)/rho(i))
        e_int(i)=pst(i)/(rho(i)*(gamm(i)-1.0d0))
      end do
      ! print *,i-1,asou(i-1),e_int(i-1)

      do i=nxi/2+1,nxi
        gamm(i)=1.4d0
        rho(i)=1.25d-1!1.25d-1
        pst(i)=1.0d-1!1.0d-1
        ux(i)=0.0d0
        ! tst(i)=pst(i)/(rho(i)*R)
        tst(i)=0.0d0
        asou(i)=sqrt(gamm(i)*pst(i)/rho(i))
        e_int(i)=pst(i)/(rho(i)*(gamm(i)-1.0d0))
      end do

      ! do i=1,nxi/3
      !   gamm(i)=1.4d0
      !   rho(i)=1.25d-1
      !   pst(i)=1.0d-1
      !   ux(i)=0.0d0
      !   ! tst(i)=pst(i)/(rho(i)*R)
      !   tst(i)=0.0d0
      !   asou(i)=sqrt(gamm(i)*pst(i)/rho(i))
      !   e_int(i)=pst(i)/(rho(i)*(gamm(i)-1.0d0))
      ! end do

      ! do i=nxi/3+1,nxi*2/3
      !   gamm(i)=1.4d0
      !   rho(i)=1.0d0
      !   pst(i)=1.0d0
      !   ux(i)=0.0d0
      !   ! tst(i)=pst(i)/(rho(i)*R)
      !   tst(i)=0.0d0
      !   asou(i)=sqrt(gamm(i)*pst(i)/rho(i))
      !   e_int(i)=pst(i)/(rho(i)*(gamm(i)-1.0d0))
      ! end do

      ! do i=nxi*2/3+1,nxi
      !   gamm(i)=1.4d0
      !   rho(i)=1.25d-1
      !   pst(i)=1.0d-1
      !   ux(i)=0.0d0
      !   ! tst(i)=pst(i)/(rho(i)*R)
      !   tst(i)=0.0d0
      !   asou(i)=sqrt(gamm(i)*pst(i)/rho(i))
      !   e_int(i)=pst(i)/(rho(i)*(gamm(i)-1.0d0))
      ! end do

      ! do i=1,nxi/3
      !   gamm(i)=1.4d0
      !   rho(i)=1.0d0
      !   pst(i)=1.0d0
      !   ux(i)=0.0d0
      !   ! tst(i)=pst(i)/(rho(i)*R)
      !   tst(i)=0.0d0
      !   asou(i)=sqrt(gamm(i)*pst(i)/rho(i))
      !   e_int(i)=pst(i)/(rho(i)*(gamm(i)-1.0d0))
      ! end do

      ! do i=nxi/3+1,nxi*2/3
      !   gamm(i)=1.4d0
      !   rho(i)=0.90d0
      !   pst(i)=0.90d0
      !   ux(i)=5.0d-1
      !   ! tst(i)=pst(i)/(rho(i)*R)
      !   tst(i)=0.0d0
      !   asou(i)=sqrt(gamm(i)*pst(i)/rho(i))
      !   e_int(i)=pst(i)/(rho(i)*(gamm(i)-1.0d0))
      ! end do

      ! do i=nxi*2/3+1,nxi
      !   gamm(i)=1.4d0
      !   rho(i)=1.0d0
      !   pst(i)=1.0d0
      !   ux(i)=0.0d0
      !   ! tst(i)=pst(i)/(rho(i)*R)
      !   tst(i)=0.0d0
      !   asou(i)=sqrt(gamm(i)*pst(i)/rho(i))
      !   e_int(i)=pst(i)/(rho(i)*(gamm(i)-1.0d0))
      ! end do
      ! print *,asou(nxi),e_int(nxi)
      
    endif
    
end subroutine Initialize

subroutine Ausm(nxi,rho,ux,pst,e_int,gamm,asou,convx)
  implicit none
  integer :: nxi
  double precision :: rho(nxi),ux(nxi),pst(nxi),e_int(nxi)
  double precision :: gamm(nxi),asou(nxi),convx(nxi,3)
  integer :: i
  double precision :: rhoL,rhoR,aL,aR,uL,uR,enthL,enthR,MachL,MachR
  double precision :: MLP,MRM,Mhalf,pL,pR,pLP,pRM

  do i=2,nxi
    rhoL=rho(i-1)
    rhoR=rho(i)
    aL=asou(i-1)
    aR=asou(i)
    uL=ux(i-1)
    uR=ux(i)
    enthL=e_int(i-1)+0.5d0*uL*uL+pst(i-1)/rhoL
    enthR=e_int(i)+0.5d0*uR*uR+pst(i)/rhoR
    MachL=uL/aL
    MachR=uR/aR
    pL=pst(i-1)
    pR=pst(i)
    
    MLP=0.5d0*(MachL+abs(MachL))
    pLP=pL*MLP/MachL
    if(abs(MachL).le.1.0d0)then
      MLP=0.25d0*(MachL+1.0d0)*(MachL+1.0d0)
      pLP=pL*MLP*(2.0d0-MachL)
    endif
    
    MRM=0.5d0*(MachR-abs(MachR))
    pRM=pR*MRM/MachR
    if(abs(MachR).le.1.0d0)then
      MRM=-0.25d0*(MachR-1.0d0)*(MachR-1.0d0)
      pRM=-pR*MRM*(2.0d0+MachR)
    endif

    Mhalf=MLP+MRM

    convx(i,1)=Mhalf*0.5d0*(rhoL*aL+rhoR*aR)&
             &-0.5d0*abs(Mhalf)*(rhoR*aR-rhoL*aL)
    convx(i,2)=Mhalf*0.5d0*(rhoL*aL*uL+rhoR*aR*uR)&
             &-0.5d0*abs(Mhalf)*(rhoR*aR*uR-rhoL*aL*uL)&
             &+pLP+pRM
    convx(i,3)=Mhalf*0.5d0*(rhoL*aL*enthL+rhoR*aR*enthR)&
             &-0.5d0*abs(Mhalf)*(rhoR*aR*enthR-rhoL*aL*enthL)
  
  end do
  ! print *, convx(499,2),convx(500,2),convx(501,2)

end subroutine Ausm

subroutine Ausm_w_MUSCL(nxi,rho,ux,pst,e_int,gamm,asou,convx)
  implicit none
  integer :: nxi
  double precision :: rho(nxi),ux(nxi),pst(nxi),e_int(nxi)
  double precision :: gamm(nxi),asou(nxi),convx(nxi,3)
  integer :: i
  double precision :: rhoL,rhoR,aL,aR,uL,uR,enth_m2,enth_m1,enth_0,enth_1,MachL,MachR
  double precision :: MLP,MRM,Mhalf,pL,pR,pLP,pRM,enthL,enthR

  do i=3,nxi-1
    rhoL=rho(i-1)+0.5d0*minmod((rho(i-1)-rho(i-2))/(rho(i)-rho(i-1)+1.0d-15))*(rho(i-1)-rho(i-2))
    rhoR=rho(i)-0.5d0*minmod((rho(i)-rho(i-1))/(rho(i+1)-rho(i)+1.0d-15))*(rho(i+1)-rho(i))
    aL=asou(i-1)+0.5d0*minmod((asou(i-1)-asou(i-2))/(asou(i)-asou(i-1)+1.0d-15))*(asou(i-1)-asou(i-2))
    aR=asou(i)-0.5d0*minmod((asou(i)-asou(i-1))/(asou(i+1)-asou(i)+1.0d-15))*(asou(i+1)-asou(i))
    uL=ux(i-1)+0.5d0*minmod((ux(i-1)-ux(i-2))/(ux(i)-ux(i-1)+1.0d-15))*(ux(i-1)-ux(i-2))
    uR=ux(i)-0.5d0*minmod((ux(i)-ux(i-1))/(ux(i+1)-ux(i)+1.0d-15))*(ux(i+1)-ux(i))
    enth_m2=e_int(i-2)+0.5d0*ux(i-2)*ux(i-2)+pst(i-2)/rho(i-2)
    enth_m1=e_int(i-1)+0.5d0*ux(i-1)*ux(i-1)+pst(i-1)/rho(i-1)
    enth_0=e_int(i)+0.5d0*ux(i)*ux(i)+pst(i)/rho(i)
    enth_1=e_int(i+1)+0.5d0*ux(i+1)*ux(i+1)+pst(i+1)/rho(i+1)
    enthL=enth_m1+0.5d0*minmod((enth_m1-enth_m2)/(enth_0-enth_m1+1.0d-15))*(enth_m1-enth_m2)
    enthR=enth_0-0.5d0*minmod((enth_0-enth_m1)/(enth_1-enth_0+1.0d-15))*(enth_1-enth_0)
    MachL=uL/aL
    MachR=uR/aR
    pL=pst(i-1)+0.5d0*minmod((pst(i-1)-pst(i-2))/(pst(i)-pst(i-1)+1.0d-15))*(pst(i-1)-pst(i-2))
    pR=pst(i)-0.5d0*minmod((pst(i)-pst(i-1))/(pst(i+1)-pst(i)+1.0d-15))*(pst(i+1)-pst(i))
    
    MLP=0.5d0*(MachL+abs(MachL))
    pLP=pL*MLP/MachL
    if(abs(MachL).le.1.0d0)then
      MLP=0.25d0*(MachL+1.0d0)*(MachL+1.0d0)
      pLP=pL*MLP*(2.0d0-MachL)
    endif
    
    MRM=0.5d0*(MachR-abs(MachR))
    pRM=pR*MRM/MachR
    if(abs(MachR).le.1.0d0)then
      MRM=-0.25d0*(MachR-1.0d0)*(MachR-1.0d0)
      pRM=-pR*MRM*(2.0d0+MachR)
    endif

    Mhalf=MLP+MRM

    convx(i,1)=Mhalf*0.5d0*(rhoL*aL+rhoR*aR)&
             &-0.5d0*abs(Mhalf)*(rhoR*aR-rhoL*aL)
    convx(i,2)=Mhalf*0.5d0*(rhoL*aL*uL+rhoR*aR*uR)&
             &-0.5d0*abs(Mhalf)*(rhoR*aR*uR-rhoL*aL*uL)&
             &+pLP+pRM
    convx(i,3)=Mhalf*0.5d0*(rhoL*aL*enthL+rhoR*aR*enthR)&
             &-0.5d0*abs(Mhalf)*(rhoR*aR*enthR-rhoL*aL*enthL)
  
  end do
  ! print *, convx(499,2),convx(500,2),convx(501,2)

end subroutine Ausm_w_MUSCL

subroutine SLAU(nxi,rho,ux,pst,e_int,gamm,asou,convx)
  implicit none
  integer :: nxi
  double precision :: rho(nxi),ux(nxi),pst(nxi),e_int(nxi)
  double precision :: gamm(nxi),asou(nxi),convx(nxi,3)
  integer :: i
  double precision :: rhoM,rhoP,uxM,uxP,rhouxM,rhouxP,hM,hP,pM,pP,betaM,betaP
  double precision :: enth_m2,enth_m1,enth_0,enth_1,enthP,enthM,aP,aM,aave,machP,machM
  double precision :: kai,machHat,funcfp,deltarho,deltap,funcg,Vnbar,mdot,ptilde
  double precision :: coeffP,coeffM

  ! E. Shima, K. Kitamura,
  ! Parameter-Free Simple Low-Dissipation AUSM-Family Scheme for All Speeds
  ! AIAA Journal, Vol. 49, No. 8, pp.1693--1709, 2011.


  ! M:minus, P:plus

  do i=3,nxi-1
    rhoP=rho(i-1)+0.5d0*minmod((rho(i-1)-rho(i-2))/(rho(i)-rho(i-1)+1.0d-15))*(rho(i-1)-rho(i-2))
    rhoM=rho(i)-0.5d0*minmod((rho(i)-rho(i-1))/(rho(i+1)-rho(i)+1.0d-15))*(rho(i+1)-rho(i))
    uxP=ux(i-1)+0.5d0*minmod((ux(i-1)-ux(i-2))/(ux(i)-ux(i-1)+1.0d-15))*(ux(i-1)-ux(i-2))
    uxM=ux(i)-0.5d0*minmod((ux(i)-ux(i-1))/(ux(i+1)-ux(i)+1.0d-15))*(ux(i+1)-ux(i))
    rhouxP=rho(i-1)*ux(i-1)+0.5d0*minmod((rho(i-1)*ux(i-1)-rho(i-2)*ux(i-2))&
         &/(rho(i)*ux(i)-rho(i-1)*ux(i-1)+1.0d-15))*(rho(i-1)*ux(i-1)-rho(i-2)*ux(i-2))
    rhouxM=rho(i)*ux(i)-0.5d0*minmod((rho(i)*ux(i)-rho(i-1)*ux(i-1))&
         &/(rho(i+1)*ux(i+1)-rho(i)*ux(i)+1.0d-15))*(rho(i+1)*ux(i+1)-rho(i)*ux(i))

    enth_m2=e_int(i-2)+0.5d0*ux(i-2)*ux(i-2)+pst(i-2)/rho(i-2)
    enth_m1=e_int(i-1)+0.5d0*ux(i-1)*ux(i-1)+pst(i-1)/rho(i-1)
    enth_0=e_int(i)+0.5d0*ux(i)*ux(i)+pst(i)/rho(i)
    enth_1=e_int(i+1)+0.5d0*ux(i+1)*ux(i+1)+pst(i+1)/rho(i+1)
    enthP=enth_m1+0.5d0*minmod((enth_m1-enth_m2)/(enth_0-enth_m1+1.0d-15))*(enth_m1-enth_m2)
    enthM=enth_0-0.5d0*minmod((enth_0-enth_m1)/(enth_1-enth_0+1.0d-15))*(enth_0-enth_1)
    
    aP=asou(i-1)+0.5d0*minmod((asou(i-1)-asou(i-2))/(asou(i)-asou(i-1)+1.0d-15))*(asou(i-1)-asou(i-2))
    aM=asou(i)-0.5d0*minmod((asou(i)-asou(i-1))/(asou(i+1)-asou(i)+1.0d-15))*(asou(i+1)-asou(i))
    aave=0.5d0*(aP+aM)
    machP=uxP/aave
    machM=uxM/aave

    pP=pst(i-1)+0.5d0*minmod((pst(i-1)-pst(i-2))/(pst(i)-pst(i-1)+1.0d-15))*(pst(i-1)-pst(i-2))
    pM=pst(i)-0.5d0*minmod((pst(i)-pst(i-1))/(pst(i+1)-pst(i)+1.0d-15))*(pst(i+1)-pst(i))

    betaP=0.25d0*(2.0d0-machP)*(machP+1.0d0)*(machP+1.0d0)
    betaM=0.25d0*(2.0d0+machM)*(machM-1.0d0)*(machM-1.0d0)
    if(abs(machP).ge.1.0d0)then
      betaP=0.5d0*(1.0d0+machP/abs(machP))
    end if

    if(abs(machM).ge.1.0d0)then
      betaM=0.5d0*(1.0d0-machM/abs(machM))
    end if

    machHat=min(1.0d0,sqrt(0.5d0*(uxP*uxP+uxM*uxM))/aave)
    kai=(1.0d0-machHat)*(1.0d0-machHat)
    funcfp=1.0d0-kai
    deltarho=rhoM-rhoP
    deltap=pM-pP

    funcg=-max(min(machP,0.0d0),-1.0d0)*min(max(machM,0.0d0),1.0d0)

    Vnbar=(rhoP*abs(uxP)+rhoM*abs(uxM))/(rhoP+rhoM)

    mdot=0.5d0*(rhouxP+rhouxM-Vnbar*deltarho)*(1.0d0-funcg)-kai*deltap/(2.0d0*aave)

    coeffP=0.5d0*(mdot+abs(mdot))
    coeffM=0.5d0*(mdot-abs(mdot))

    ptilde=0.5d0*(pP+pM)+0.5d0*(betaP-betaM)*(pP-pM)+funcfp*0.5d0*(betaP+betaM-1.0d0)*(pP+pM)


    convx(i,1)=coeffP+coeffM
    convx(i,2)=coeffP*uxP+coeffM*uxM+ptilde
    convx(i,3)=coeffP*enthP+coeffM*enthM

  end do
end subroutine SLAU

subroutine AusmDV(nxi,rho,ux,pst,e_int,gamm,asou,convx)
  implicit none
  integer :: nxi
  double precision :: rho(nxi),ux(nxi),pst(nxi),e_int(nxi)
  double precision :: gamm(nxi),asou(nxi),convx(nxi,3)
  integer :: i
  double precision :: rhoL,rhoR,am,uL,uR,pL,pR,enthL,enthR
  double precision :: MachL,MachR
  double precision :: p_rh_L,p_rh_R,alphaL,alphaR
  double precision :: uLP,uRM,pLP,pRM
  double precision :: rhoUHalf,sparam,Kparam

  Kparam=10.0d0

  do i=2,nxi
    rhoL=rho(i-1)
    rhoR=rho(i)
    am=max(asou(i-1),asou(i))
    uL=ux(i-1)
    uR=ux(i)
    pL=pst(i-1)
    pR=pst(i)
    enthL=e_int(i-1)+0.5d0*uL*uL+pL/rhoL
    enthR=e_int(i)+0.5d0*uR*uR+pR/rhoR
    MachL=uL/am
    MachR=uR/am

    p_rh_L=pL/rhoL
    p_rh_R=pR/rhoR

    alphaL=2.0d0*p_rh_L/(p_rh_L+p_rh_R)
    alphaR=2.0d0-alphaL
    ! alphaR=2.0d0*p_rh_R/(p_rh_L+p_rh_R)

    uLP=alphaL*((uL+am)*(uL+am)*0.25d0/am-0.5d0*(uL+abs(uL)))+0.5d0*(uL+abs(uL))
    pLP=0.25d0*pL*(MachL+1.0d0)*(MachL+1.0d0)*(2.0d0-MachL)
    if (abs(MachL).gt.1.0d0) then
      uLP=0.5d0*(uL+abs(uL))
      pLP=pL*0.5d0*(uL+abs(uL))/uL
    endif

    uRM=alphaR*(-(uR-am)*(uR-am)*0.25d0/am-0.5d0*(uR-abs(uR)))+0.5d0*(uR-abs(uR))
    pRM=0.25d0*pR*(MachR-1.0d0)*(MachR-1.0d0)*(2.0d0+MachR)
    if (abs(MachR).gt.1.0d0) then
      uRM=0.5d0*(uR-abs(uR))
      pRM=pR*0.5d0*(uR-abs(uR))/uR
    endif

    rhoUHalf=uLP*rhoL+uRM*rhoR

    sparam=0.5d0*min(1.0d0,Kparam*abs(pR-pL)/min(pL,pR))
    ! if(i.eq.500)print *, sparam
    convx(i,1)=rhoUHalf
    convx(i,2)=(0.5d0-sparam)*0.5d0*(rhoUHalf*(uL+uR)-abs(rhoUHalf)*(uR-uL))&
             &+(0.5d0+sparam)*(uLP*rhoL*uL+uRM*rhoR*uR)&
             &+pLP+pRM
    convx(i,3)=0.5d0*(rhoUHalf*(enthL+enthR)-abs(rhoUHalf)*(enthR-enthL))

  end do

end subroutine AusmDV

subroutine AusmDV_w_MUSCL(nxi,rho,ux,pst,e_int,gamm,asou,convx)
  implicit none
  integer :: nxi
  double precision :: rho(nxi),ux(nxi),pst(nxi),e_int(nxi)
  double precision :: gamm(nxi),asou(nxi),convx(nxi,3)
  integer :: i
  double precision :: rhoL,rhoR,am,uL,uR,pL,pR,enthL,enthR,enth_m2,enth_m1,enth_0,enth_1
  double precision :: MachL,MachR
  double precision :: p_rh_L,p_rh_R,alphaL,alphaR
  double precision :: uLP,uRM,pLP,pRM
  double precision :: rhoUHalf,sparam,Kparam

  Kparam=10.0d0

  do i=3,nxi-1
    rhoL=rho(i-1)+0.5d0*minmod((rho(i-1)-rho(i-2))/(rho(i)-rho(i-1)+1.0d-15))*(rho(i-1)-rho(i-2))
    rhoR=rho(i)-0.5d0*minmod((rho(i)-rho(i-1))/(rho(i+1)-rho(i)+1.0d-15))*(rho(i+1)-rho(i))
    am=max(asou(i-1)+0.5d0*minmod((asou(i-1)-asou(i-2))/(asou(i)-asou(i-1)+1.0d-15))*(asou(i-1)-asou(i-2)),&
      &    asou(i)-0.5d0*minmod((asou(i)-asou(i-1))/(asou(i+1)-asou(i)+1.0d-15))*(asou(i+1)-asou(i)))
    uL=ux(i-1)+0.5d0*minmod((ux(i-1)-ux(i-2))/(ux(i)-ux(i-1)+1.0d-15))*(ux(i-1)-ux(i-2))
    uR=ux(i)-0.5d0*minmod((ux(i)-ux(i-1))/(ux(i+1)-ux(i)+1.0d-15))*(ux(i+1)-ux(i))
    pL=pst(i-1)+0.5d0*minmod((pst(i-1)-pst(i-2))/(pst(i)-pst(i-1)+1.0d-15))*(pst(i-1)-pst(i-2))
    pR=pst(i)-0.5d0*minmod((pst(i)-pst(i-1))/(pst(i+1)-pst(i)+1.0d-15))*(pst(i+1)-pst(i))
    enth_m2=e_int(i-2)+0.5d0*ux(i-2)*ux(i-2)+pst(i-2)/rho(i-2)
    enth_m1=e_int(i-1)+0.5d0*ux(i-1)*ux(i-1)+pst(i-1)/rho(i-1)
    enth_0=e_int(i)+0.5d0*ux(i)*ux(i)+pst(i)/rho(i)
    enth_1=e_int(i+1)+0.5d0*ux(i+1)*ux(i+1)+pst(i+1)/rho(i+1)
    enthL=enth_m1+0.5d0*minmod((enth_m1-enth_m2)/(enth_0-enth_m1+1.0d-15))*(enth_m1-enth_m2)
    enthR=enth_0-0.5d0*minmod((enth_0-enth_m1)/(enth_1-enth_0+1.0d-15))*(enth_1-enth_0)
    MachL=uL/am
    MachR=uR/am

    p_rh_L=pL/rhoL
    p_rh_R=pR/rhoR

    alphaL=2.0d0*p_rh_L/(p_rh_L+p_rh_R)
    alphaR=2.0d0-alphaL
    ! alphaR=2.0d0*p_rh_R/(p_rh_L+p_rh_R)

    uLP=alphaL*((uL+am)*(uL+am)*0.25d0/am-0.5d0*(uL+abs(uL)))+0.5d0*(uL+abs(uL))
    pLP=0.25d0*pL*(MachL+1.0d0)*(MachL+1.0d0)*(2.0d0-MachL)
    if (abs(MachL).gt.1.0d0) then
      uLP=0.5d0*(uL+abs(uL))
      pLP=pL*0.5d0*(uL+abs(uL))/uL
    endif

    uRM=alphaR*(-(uR-am)*(uR-am)*0.25d0/am-0.5d0*(uR-abs(uR)))+0.5d0*(uR-abs(uR))
    pRM=0.25d0*pR*(MachR-1.0d0)*(MachR-1.0d0)*(2.0d0+MachR)
    if (abs(MachR).gt.1.0d0) then
      uRM=0.5d0*(uR-abs(uR))
      pRM=pR*0.5d0*(uR-abs(uR))/uR
    endif

    rhoUHalf=uLP*rhoL+uRM*rhoR

    sparam=0.5d0*min(1.0d0,Kparam*abs(pR-pL)/min(pL,pR))
    ! if(i.eq.500)print *, sparam
    convx(i,1)=rhoUHalf
    convx(i,2)=(0.5d0-sparam)*0.5d0*(rhoUHalf*(uL+uR)-abs(rhoUHalf)*(uR-uL))&
             &+(0.5d0+sparam)*(uLP*rhoL*uL+uRM*rhoR*uR)&
             &+pLP+pRM
    convx(i,3)=0.5d0*(rhoUHalf*(enthL+enthR)-abs(rhoUHalf)*(enthR-enthL))

  end do

end subroutine AusmDV_w_MUSCL

double precision function minmod(r)
    implicit none
    double precision :: r
    minmod=max(0.0d0,min(1.0d0,r))

end function minmod

double precision function superbee(r)
    implicit none
    double precision :: r
    superbee=max(0.0d0,min(2.0d0*r,1.0d0),min(r,2.0d0))

end function superbee

subroutine TimeIntEulerExplicit(nxi,xgrid,convx,dU,Q,dt)
    implicit none
    integer :: nxi
    double precision, intent(in) :: xgrid(nxi+1),convx(nxi,3)
    double precision, intent(out) ::dU(2:nxi-1,3),Q(2:nxi-1,3)
    double precision, intent(in) :: dt
    integer i,j
    double precision :: dx

    do i=2,nxi-1
      dx=xgrid(i+1)-xgrid(i)
      dU(i,1)=dt/dx*(convx(i,1)-convx(i+1,1))
      dU(i,2)=dt/dx*(convx(i,2)-convx(i+1,2))
      dU(i,3)=dt/dx*(convx(i,3)-convx(i+1,3))
    end do

    do i=2,nxi-1
      Q(i,1)=rho(i)+dU(i,1)
      Q(i,2)=rho(i)*ux(i)+dU(i,2)
      Q(i,3)=rho(i)*(e_int(i)+(0.5d0*ux(i)*ux(i)))+dU(i,3)
    end do
    ! print*, Q(500,2),dU(500,2),xcell(500)

    
end subroutine TimeIntEulerExplicit

subroutine TimeIntEulerExplicit_2nd(nxi,xgrid,convx,dU,Q,dt)
    implicit none
    integer :: nxi
    double precision, intent(in) :: xgrid(nxi+1),convx(nxi,3)
    double precision, intent(out) ::dU(2:nxi-1,3),Q(2:nxi-1,3)
    double precision, intent(in) :: dt
    integer i,j
    double precision :: dx

    do i=3,nxi-2
      dx=xgrid(i+1)-xgrid(i)
      dU(i,1)=dt/dx*(convx(i,1)-convx(i+1,1))
      dU(i,2)=dt/dx*(convx(i,2)-convx(i+1,2))
      dU(i,3)=dt/dx*(convx(i,3)-convx(i+1,3))
    end do

    do i=3,nxi-2
      Q(i,1)=rho(i)+dU(i,1)
      Q(i,2)=rho(i)*ux(i)+dU(i,2)
      Q(i,3)=rho(i)*(e_int(i)+(0.5d0*ux(i)*ux(i)))+dU(i,3)
    end do
    ! print*, Q(500,2),dU(500,2),xcell(500)

    
end subroutine TimeIntEulerExplicit_2nd

subroutine Result(nxi,rho,ux,pst,e_int,gamm,asou,tst,dU,Q)
  implicit none
  integer :: nxi
  double precision :: rho(nxi),ux(nxi),pst(nxi),e_int(nxi)
  double precision :: gamm(nxi),asou(nxi),tst(nxi)
  double precision :: dU(2:nxi-1,3),Q(2:nxi-1,3)
  integer :: i
  double precision :: R

  ! R=208.2d0

  do i=2,nxi-1
    rho(i)=Q(i,1)
    ux(i)=Q(i,2)/Q(i,1)
    e_int(i)=Q(i,3)/Q(i,1)-0.5d0*(ux(i)*ux(i))
    !!gamm doesn't change
    pst(i)=(gamm(i)-1.0d0)*rho(i)*e_int(i)
    ! tst(i)=pst(i)/(rho(i)*R)
    asou(i)=sqrt(gamm(i)*pst(i)/rho(i))

  end do

  !boundary
  ux(1)=0.0d0
  ux(nxi)=0.0d0
  rho(1)=rho(2)
  rho(nxi)=rho(nxi-1)
  pst(1)=pst(2)
  pst(nxi)=pst(nxi-1)
  e_int(1)=pst(1)/(rho(1)*(gamm(1)-1.0d0))
  e_int(nxi)=pst(nxi)/(rho(nxi)*(gamm(nxi)-1.0d0))
  ! tst(1)=pst(1)/(rho(1)*R)
  ! tst(nxi)=pst(nxi)/(rho(nxi)*R)
  asou(1)=sqrt(gamm(1)*pst(1)/rho(1))
  asou(nxi)=sqrt(gamm(nxi)*pst(nxi)/rho(nxi))


    
end subroutine Result

subroutine Result_2nd(nxi,rho,ux,pst,e_int,gamm,asou,tst,dU,Q)
  implicit none
  integer :: nxi
  double precision :: rho(nxi),ux(nxi),pst(nxi),e_int(nxi)
  double precision :: gamm(nxi),asou(nxi),tst(nxi)
  double precision :: dU(2:nxi-1,3),Q(2:nxi-1,3)
  integer :: i
  double precision :: R

  ! R=208.2d0

  do i=3,nxi-2
    rho(i)=Q(i,1)
    ux(i)=Q(i,2)/Q(i,1)
    e_int(i)=Q(i,3)/Q(i,1)-0.5d0*(ux(i)*ux(i))
    !!gamm doesn't change
    pst(i)=(gamm(i)-1.0d0)*rho(i)*e_int(i)
    ! tst(i)=pst(i)/(rho(i)*R)
    asou(i)=sqrt(gamm(i)*pst(i)/rho(i))

  end do

  !boundary
  ux(1)=0.0d0
  ux(2)=0.0d0
  ux(nxi-1)=0.0d0
  ux(nxi)=0.0d0
  rho(2)=rho(3)
  rho(1)=rho(3)
  rho(nxi-1)=rho(nxi-2)
  rho(nxi)=rho(nxi-2)
  pst(2)=pst(3)
  pst(1)=pst(3)
  pst(nxi-1)=pst(nxi-2)
  pst(nxi)=pst(nxi-2)
  e_int(1)=pst(1)/(rho(1)*(gamm(1)-1.0d0))
  e_int(2)=pst(2)/(rho(2)*(gamm(2)-1.0d0))
  e_int(nxi-1)=pst(nxi-1)/(rho(nxi-1)*(gamm(nxi-1)-1.0d0))
  e_int(nxi)=pst(nxi)/(rho(nxi)*(gamm(nxi)-1.0d0))
  ! tst(1)=pst(1)/(rho(1)*R)
  ! tst(nxi)=pst(nxi)/(rho(nxi)*R)
  asou(1)=sqrt(gamm(1)*pst(1)/rho(1))
  asou(2)=sqrt(gamm(2)*pst(2)/rho(2))
  asou(nxi)=sqrt(gamm(nxi)*pst(nxi)/rho(nxi))
  asou(nxi-1)=sqrt(gamm(nxi-1)*pst(nxi-1)/rho(nxi-1))

    
end subroutine Result_2nd

subroutine Output(nxi,xcell,rho,ux,pst,e_int,tst,fileNo)
  implicit none
  integer :: nxi
  integer :: fileNo
  double precision :: xcell(nxi),rho(nxi),ux(nxi)
  double precision :: pst(nxi),e_int(nxi),tst(nxi)
  integer i
  character(30) :: fileName

  write(fileName,*)fileNo

  open(301,file='./'//trim(adjustl(fileName))//'_2nd.dat',status='unknown')
  write(301,*)"1:x,          2:rho,       3:ux,        4:pst,       5:e_int,     6:Mach"

  do i=2,nxi-1
    write(301,'(1p6d13.5)')xcell(i),rho(i),ux(i),pst(i),e_int(i),abs(ux(i)/asou(i))
  end do

  fileNo=fileNo+1
    
end subroutine Output

end module Mesh_Cell_PhysParam

program CFD
  use Mesh_Cell_PhysParam
  implicit none
  integer, i

  !!!setup!!!
  !!!reading numerical conditions!!
  call ReadCond(nstep,timing,fileNo,dt,nxi,keizoku)
  !!!allocating mesh,cell,PhysParam
  call TermsAllAllocate(nxi,convx,dU)
  call AllAllocate(nxi,xgrid,xcell,rho,ux,pst,e_int,gamm,asou,tst)
  !!!reading mesh!!!
  call ReadGrid(nxi,xgrid,xcell)
  ! print *, "ReadGrid has finished."

  !!!initialize!!!
  call Initialize(nxi,keizoku,rho,ux,pst,e_int,gamm,asou,tst)
  
  do i=1,nstep

    ! call Ausm(nxi,rho,ux,pst,e_int,gamm,asou,convx)
    ! call Ausm_w_MUSCL(nxi,rho,ux,pst,e_int,gamm,asou,convx)
    call AusmDV_w_MUSCL(nxi,rho,ux,pst,e_int,gamm,asou,convx)
    ! call SLAU(nxi,rho,ux,pst,e_int,gamm,asou,convx)
    ! call AusmDV(nxi,rho,ux,pst,e_int,gamm,asou,convx)

    call TimeIntEulerExplicit_2nd(nxi,xgrid,convx,dU,Q,dt)
    ! call TimeIntEulerExplicit(nxi,xgrid,convx,dU,Q,dt)

    ! call Result(nxi,rho,ux,pst,e_int,gamm,asou,tst,dU,Q)
    call Result_2nd(nxi,rho,ux,pst,e_int,gamm,asou,tst,dU,Q)

    if(mod(i,timing).eq.0)then
      call Output(nxi,xcell,rho,ux,pst,e_int,tst,fileNo)
    endif

  end do

  ! call Output(nxi,xcell,rho,ux,pst,e_int,tst,fileNo)


  ! open(100,file='./grid.dat',status='unknown')
  ! do i=1,nxi+1
  !   write(100,"(1pd15.7)") dble(i-1)/dble(nxi)-0.5d0
  ! end do
  ! close(100)




end program CFD



