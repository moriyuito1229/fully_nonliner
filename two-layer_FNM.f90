
!  2D 2-layer flow past an obstacle
!  fully nonlinear models
!----------------------------------------------------------------------
      implicit none
!----- parameter ------------------------------------------------------
integer, parameter :: &
       imax=6000*4,  &
       itewri=5000*10,&
       itemn =5000*600
       ! itemn =20
double precision,parameter ::  &
       PI=3.14159,            &
       ! den1=0.1d0,           &
       d2=0.2d0,d1=1.0d0-d2  , &
       ! d2=0.3d0,d1=0.1d0  , &
       den1=0.1d0,           &
       ! h0=0.00
       ! h0=10d0**(-2d0),    &
       h0=0.005d0,    &
       ccc=0.3d0,    &
       Frc=1.0d0,             &          ! Frc=U/c
       dt=2.0d0*10d0**(-4d0),  &
       xmax=300.0d0,          &
       ! xmin=-100.0d0,         &
       xmin=-300.0d0,         &
       dx=(xmax-xmin)/dble(imax), &
       nva=1.0d0                   ! coefficient of numerical viscous term
!----------------------------------------------------------------------
double precision,dimension (-2:imax+1) :: &
       x,        &
       zeta,     &
       zeta0,zeta00,zeta000,&
       zeta0t,&
       f,    &  ! f(n+1)
       f0,f00,f000 ,  &  ! f(n)
       h,&
       hx,&
       hxx,&
       hxxx,&
       w10,w10x,&
       u1,&
       u10,u100,&
       u2,&
       u20,u200,&
       eta1,eta2,&
       eta10,eta20,&
       alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,&
       beta1,beta2,beta3,beta4,beta5,beta6,&
       r ,px,p2,G1,G2,F1,F2,K1,K2,&
       check1,check2
integer ::  &
   i,j,k,iten,itime,iet,ist,irt,&
   iet2,ist2,irt2,itmax,ith,itm,its
double precision:: fmax,Fd,c1,&
          zeta0x,zeta0xx,zeta0xxx,&
          zetax,&
          fx,fxx,fxxx,&
          f0x,f0xx,f0xxx,&
          ft,ft0,ft00,&
          ftt0,&
          f0xt,f0xxt,f0tt,f0xtt,&
          u10x,u10xx,u10xxx,u10t,u10xt,u10xxt,u20x,u20xx,u20xxx,u20t,u20xt,&
          u200x,&
          u10u10x,u20u20x,w10w10x,&
          eta10x,eta20x,&
          eta10xu10,eta20xu20,&
          r1,r2

character(7) citen

!----- matrix ---------------------------------------------------------
integer,parameter :: kl=3,ku=2,ldab=2*kl+ku+1
double precision, dimension(1:imax-1) ::   b1,d,du,dl, &
                                           AA,BB,CC,DD,EE,FF,&
                                           HH,II,JJ,KK,LL,MM,G,O
double precision, dimension(1:2*(imax-1)) :: b2
double precision, dimension(1:ldab,1:2*(imax-1)) ::  ab

integer , dimension(1:2*(imax-1)) :: ipiv
integer :: info

!----- Froude number --------------------------------------------------
      ! c1=sqrt((1.0d0-den1)/(den1/d1+1.0d0/d2))         ! c=c/sqrt(gD)
      c1=dsqrt(0.5d0*(1d0-dsqrt(1-4d0*d1*d2*(1-den1))))   !c-
      ! c1=dsqrt(0.5d0*(1d0+dsqrt(1d0-4d0*d1*d2*(1d0-den1))))   !c+

      Fd=c1*Frc                                        ! Fd=U/sqrt(gD)

!----- grid generation ------------------------------------------------
      do i=0,imax
        x(i)=xmin+dx*dble(i)
      end do

      x(-1)=2.0d0*x(0)-x(1)
      x(imax+1)=2.0d0*x(imax)-x(imax-1)

!----- printing of the parameter --------------------------------------
      write(6,2001) Frc,dt,dx,itemn,h0,ccc,d2,den1
 2001 format('Frc=',E10.2,/,'dt=',E10.2,/,'dx=',E10.2,/,&
            'itemn=',I10,/,'h0=',E10.2,/,'c=',E10.2,/, &
            'd2=',E10.2,/,'den1=',E10.2)

!----- initial condition ----------------------------------------------
      do i=-2,imax+1
        zeta(i)=0.0d0
        zeta0(i)=0.0d0
        zeta00(i)=0.0d0
        f(i)=0.0d0
        f0(i)=0.0d0
        f00(i)=0.0d0
        f000(i)=0.0d0
        ! ft(i)=0.0d0
        ! ft0(i)=0d0
        ! ft00(i)=0d0
        ! ftt0(i)=0d0
        ! h(i)=h0*(cosh(ccc*x(i))**(-2.0d0))
        h(i)=h0*(dcosh(ccc*x(i)))**(-2d0)
        ! u1(i)=1.0d0/(1.0d0-h(i))
        ! u2(i)=1.0d0/(1.0d0-h(i))
        u1(i)=1d0
        u2(i)=1d0
        u10(i)=u1(i)
        u20(i)=u2(i)
        u100(i)=u10(i)
        u200(i)=u20(i)
        eta10(i)=d1
        eta20(i)=d2-h(i)
        eta1(i)=d1
        eta2(i)=d2-h(i)
        w10(i)=0d0
        w10x(i)=0d0
      end do




      do i=1,imax
        hx(i)=(h(i+1)-h(i-1))/dx*0.5d0
        hxx(i)=(h(i+1)-2.0d0*h(i)+h(i-1))/(dx**2d0)
      end do

        hx(imax+1)=2d0*hx(imax)-hx(imax-1)
        hxx(imax+1)=2d0*hxx(imax)-hxx(imax-1)

      do i=1,imax-1
        hxxx(i)=(h(i+2)-2.0d0*h(i+1)         &
                +2.0d0*h(i-1)-h(i-2))/(2d0*dx**3d0)
      end do

      hxxx(imax)=2d0*hxxx(imax-1)-hxxx(imax-2)
      hxxx(imax+1)=2d0*hxxx(imax)-hxxx(imax-1)

      open(4,file='obstacle.dat')
      do i=1,imax-1
        write(4,*) x(i),h(i),hx(i)
      end do
      close(4)


!----- time development -----------------------------------------------
      call system_clock(ist)      !   monitor calculation time
      call system_clock(ist2)

      iten=0
      do while(iten<=itemn)

!----- calculation of eta1 -----------------------------------------------

!----- eta1 on the upstream boundary -------------------------------------
     eta1(0)=d1
     eta1(-1)=d1
     ! eta1(0)=2.0d0*eta1(1)-eta1(2)
     ! eta1(-1)=2.0d0*eta1(0)-eta1(1)
!----- eta1 on the downstream boundary -----------------------------------
     ! eta1(imax)=d1
     ! eta1(imax+1)=d1

     eta1(imax)=2.0d0*eta1(imax-1)-eta1(imax-2)
     eta1(imax+1)=2.0d0*eta1(imax)-eta1(imax-1)

!----- Explicit method --------------------------------------------

      do i=1,imax-1
        u10x=(u10(i+1)-u10(i-1))/dx*0.5d0
        if ( u10(i)>=0 ) then
          eta10xu10=u10(i)*(2d0*eta10(i-2)-10d0*eta10(i-1)+9d0*eta10(i)-2d0*eta10(i+1)+eta10(i+2))/(6d0*dx)
        else
          eta10xu10=u10(i)*(-eta10(i-2)+2d0*eta10(i-1)-9d0*eta10(i)+10d0*eta10(i+1)-2d0*eta10(i+2))/(6d0*dx)
        end if
        eta1(i)=eta10(i) &
               -dt*(eta10(i)*u10x+eta10xu10)

      end do

      !----- calculation of eta2-----------------------------------------------

      !----- eta1 on the upstream boundary -------------------------------------
           eta2(0)=d2-h(0)
           eta2(-1)=d2-h(-1)
           ! eta2(0)=2.0d0*eta2(1)-eta2(2)
           ! eta2(-1)=2.0d0*eta2(0)-eta2(1)
      !----- eta1 on the downstream boundary -----------------------------------
           ! eta2(imax)=d2-h(imax)
           ! eta2(imax+1)=d2-h(imax+1)

           eta2(imax)=2.0d0*eta2(imax-1)-eta2(imax-2)
           eta2(imax+1)=2.0d0*eta2(imax)-eta2(imax-1)
      !----- Explicit method --------------------------------------------

      do i=1,imax-1
        u20x=(u20(i+1)-u20(i-1))/dx*0.5d0
        if ( u20(i)>=0 ) then
          eta20xu20=u20(i)*(2d0*eta20(i-2)-10d0*eta20(i-1)+9d0*eta20(i)-2d0*eta20(i+1)+eta20(i+2))/(6d0*dx)
        else
          eta20xu20=u20(i)*(-eta20(i-2)+2d0*eta20(i-1)-9d0*eta20(i)+10d0*eta20(i+1)-2d0*eta20(i+2))/(6d0*dx)
        end if
        eta2(i)=eta20(i) &
               -dt*(eta20(i)*u20x+eta20xu20)

      end do
  !----- calculation of f,zeta-----------------------------------------------

            do i=1,imax-1
              f(i)=eta2(i)-d2+h(i)
              ! zeta(i)=eta1(i)-d1+f(i)
               zeta(i)=eta1(i)-d1+eta2(i)-d2+h(i)
            end do

            f(0)=0.0d0
            f(-1)=0.0d0
            ! f(imax)=0d0
            ! f(imax+1)=0d0

            ! f(0)=2.0d0*f(1)-f(2)
            ! f(-1)=2.0d0*f(0)-f(1)
            f(imax)=2.0d0*f(imax-1)-f(imax-2)
            f(imax+1)=2.0d0*f(imax)-f(imax-1)


            zeta(0)=0.0d0
            zeta(-1)=0.0d0
            ! zeta(0)=2.0d0*zeta(1)-zeta(2)
            ! zeta(-1)=2.0d0*zeta(0)-zeta(1)
            ! zeta(imax)=0d0
            ! zeta(imax+1)=0d0

            zeta(imax)=2.0d0*zeta(imax-1)-zeta(imax-2)
            zeta(imax+1)=2.0d0*zeta(imax)-zeta(imax-1)


            ! ! !----- calculation of ft -----------------------------------------------
            ! if(iten==0) then
            !
            !   !----- ft0 on the upstream boundary -------------------------------------
            !               ft0(0)=0.0d0
            !               ft0(-1)=0.0d0
            !               ! ft0(0)=2.0d0*ft0(1)-ft0(2)
            !               ! ft0(-1)=2.0d0*ft0(0)-ft0(1)
            !   !----- ft0 on the downstream boundary -----------------------------------
            !               ft0(imax)=2.0d0*ft0(imax-1)-ft0(imax-2)
            !               ft0(imax+1)=2.0d0*ft0(imax)-ft0(imax-1)
            !               ! ft0(imax)=0d0
            !               ! ft0(imax+1)=0d0
            !
            !   do i=1,imax-1
            !   !
            !   u20x=(u20(i+1)-u20(i-1))/dx*0.5d0
            !   ! if ( u20(i)>=0 ) then
            !   !   eta20xu20=u20(i)*(2d0*eta20(i-2)-10d0*eta20(i-1)+9d0*eta20(i)-2d0*eta20(i+1)+eta20(i+2))/(6d0*dx)
            !   ! else
            !   !   eta20xu20=u20(i)*(-eta20(i-2)+2d0*eta20(i-1)-9d0*eta20(i)+10d0*eta20(i+1)-2d0*eta20(i+2))/(6d0*dx)
            !   ! end if
            !
            !   eta20xu20 = (u20(i)*(-eta20(i+2)+8.0d0*eta20(i+1)   &
            !               -8.0d0*eta20(i-1)+eta20(i-2))                 &
            !               +nva*dabs(u20(i))*(eta20(i+2)-4.0d0*eta20(i+1)+6.0d0*eta20(i)  &
            !                            -4.0d0*eta20(i-1)+eta20(i-2)))/(12.0d0*dx)
            !
            !
            !   ft0(i)= -(eta20(i)*u20x+eta20xu20)
            !
            !   end do
            !   !  ft0(imax+1)=2.0d0*ft0(imax-1)-ft0(imax-2)
            ! end if
            ! !----- ft0 on the downstream boundary -----------------------------------
            !             ft0(imax)=2.0d0*ft0(imax-1)-ft0(imax-2)
            !             ft0(imax+1)=2.0d0*ft0(imax)-ft0(imax-1)
            !             ! ft0(imax)=0d0
            !             ! ft0(imax+1)=0d0
            !
            ! !----- ft on the upperstream boundary -----------------------------------
            ! !
            ! ft(0)=0.0d0
            ! ft(-1)=0.0d0
            ! ! ft(0)=2.0d0*ft(1)-ft(2)
            ! ! ft(-1)=2.0d0*ft(0)-ft(1)
            ! !----- ft on the downstream boundary -----------------------------------
            ! ft(imax)=2.0d0*ft(imax-1)-ft(imax-2)
            ! ft(imax+1)=2.0d0*ft(imax)-ft(imax-1)
            ! ! ft(imax)=0d0
            ! ! ft(imax+1)=0d0
            !
            !
            !   do i=1,imax-1
            !   !
            !   u20x=(u20(i+1)-u20(i-1))/dx*0.5d0
            !   ! if ( u20(i)>=0 ) then
            !   !   eta20xu20=u20(i)*(2d0*eta2(i-2)-10d0*eta2(i-1)+9d0*eta2(i)-2d0*eta2(i+1)+eta2(i+2))/(6d0*dx)
            !   ! else
            !   !   eta20xu20=u20(i)*(-eta2(i-2)+2d0*eta2(i-1)-9d0*eta2(i)+10d0*eta2(i+1)-2d0*eta2(i+2))/(6d0*dx)
            !   ! end if
            !
            !   eta20xu20 = (u20(i)*(-eta2(i+2)+8.0d0*eta2(i+1)   &
            !               -8.0d0*eta2(i-1)+eta2(i-2))                 &
            !               +nva*dabs(u20(i))*(eta2(i+2)-4.0d0*eta2(i+1)+6.0d0*eta2(i)  &
            !                            -4.0d0*eta2(i-1)+eta2(i-2)))/(12.0d0*dx)
            !
            !
            !   ft(i)= -(eta2(i)*u20x+eta20xu20)
            !
            !   end do
            !
            !   ft(imax)=2.0d0*ft(imax-1)-ft(imax-2)
            !   ft(imax+1)=2.0d0*ft(imax)-ft(imax-1)
            !   ! ft(imax)=0d0
            !   ! ft(imax+1)=0d0
            !


  ! b1(1)=b1(1)+u20(1)/dx*0.25d0*f(0)
!   b1(imax-1)=b1(imax-1)-u20(imax-1)/dx*0.25d0*ft(imax)
!
! !----- solve the equation (tridiagonal matrix) ------------------------
!   call dgtsv(imax-1,1,dl,d,du,b1,imax-1,info)
!
!   do i=1,imax-1
!     ft(i)=b1(i)
!   end do





  !----- calculation of ftt -----------------------------------------------
    ! do i=1,imax
    !
    !   f0x=(f0(i+1)-f0(i-1))/dx*0.5d0
    !   f0xt=(ft0(i+1)-ft0(i-1))/(2d0*dx)
    !   u20x=(u20(i+1)-u20(i-1))/dx*0.5d0
    !   u200x=(u200(i+1)-u200(i-1))/(2d0*dx)
    !
    !
    ! ftt0(i)= -f0xt*u20(i)-(f0x-hx(i))*(u20(i)-u200(i))/dt    &
    !          -ft0(i)*u20x-(d2+f0(i)-h(i))*(u20x-u200x)/dt
    !
    !
    ! end do
    !
    !
    !
    !
!---------------calculation of px---------------------------------------------
! do i=0,imax
!   u10x=(u10(i+1)-u10(i-1))/dx*0.5d0
!
!   f0x = (f0(i+1)-f0(i-1))/dx*0.5d0
!
!   w10(i)=-u10x*(d1+zeta0(i)-f0(i))+ft0(i)+u10(i)*f0x
!
! end do
!
! do i =1 ,imax-1
!   zeta0x = (zeta0(i+1)-zeta0(i-1))/dx*0.5d0
!   u10t=(u10(i)-u100(i))/(dt)
!
! px(i)=(-((w10(i+1)**2d0-w10(i-1)**2d0)+(u10(i+1)**2d0-u10(i-1)**2d0))/(4d0*dx) &
!     -1d0/(Fd*Fd)*zeta0x-u10t)*den1
! end do

!
! px(0)=0d0
! px(1)=0d0
! px(imax)=0d0
! px(imax+1)=0d0
!
!
! do i =0 ,imax
!     u10x=(u10(i+1)-u10(i-1))/dx*0.5d0
!     u10xx=   (u10(i+1)-2.0d0*u10(i)+u10(i-1))/dx**2d0
!     f0x = (f0(i+1)-f0(i-1))/dx*0.5d0
!     f0xx = (f0(i+1)-2.0d0*f0(i)+f0(i-1))/dx**2d0
!     f0xt=(ft0(i+1)-ft0(i-1))/(2d0*dx)
!     zeta0x = (zeta0(i+1)-zeta0(i-1))/dx*0.5d0
!
! w10(i)=-u10x*(d1+zeta0(i)-f0(i))+ft0(i)+u10(i)*f0x
! w10x(i)=-u10xx*(d1+zeta0(i)-f0(i))-u10x*(zeta0x-f0x)+f0xt+u10x*f0x+u10(i)*f0xx
! end do
!
! do i= 1,imax-2
!    u10t=(u10(i)-u100(i))/(dt)
!    zeta0x =  (zeta0(i+1)-zeta0(i-1))/dx*0.5d0
!    u10u10x = (u10(i)*(-u10(i+2)+8.0d0*u10(i+1)   &
!                -8.0d0*u10(i-1)+u10(i-2))                 &
!                +nva*dabs(u10(i))*(u10(i+2)-4.0d0*u10(i+1)+6.0d0*u10(i)  &
!                             -4.0d0*u10(i-1)+u10(i-2)))/(12.0d0*dx)
!
! px(i+1)=px(i)&
!         +den1*(u10t+u10u10x+w10(i)*w10x(i)+zeta0x/(Fd**2d0))
!
!         u10t=(u10(i+1)-u100(i+1))/(dt)
!         zeta0x  = (zeta0(i+2)-zeta0(i))/dx*0.5d0
!         u10u10x = (u10(i+1)*(-u10(i+3)+8.0d0*u10(i+2)                           &
!                     -8.0d0*u10(i)+u10(i-1))                                     &
!                     +nva*dabs(u10(i+1))*(u10(i+3)-4.0d0*u10(i+2)+6.0d0*u10(i+1)  &
!                                  -4.0d0*u10(i)+u10(i-1)))/(12.0d0*dx)
!
! px(i+1)=px(i+1)-den1*(u10t+u10u10x+w10(i+1)*w10x(i+1)+zeta0x/(Fd**2d0))
!
!
!
! end do


! px(-1)=0d0
! px(0)=0d0
! px(1)=0d0
! ! px(0)=2.0d0*px(1)-px(2)
! ! px(-1)=2.0d0*px(0)-px(1)
! ! px(imax)=0d0
! ! px(imax+1)=0d0
! px(imax)=2.0d0*px(imax-1)-px(imax-2)
! px(imax+1)=2.0d0*px(imax)-px(imax-1)
!
! do i =0 ,imax
!     u10x=(u10(i+1)-u10(i-1))/dx*0.5d0
!     u10xx=   (u10(i+1)-2.0d0*u10(i)+u10(i-1))/dx**2d0
!     f0x = (f0(i+1)-f0(i-1))/dx*0.5d0
!     f0xx = (f0(i+1)-2.0d0*f0(i)+f0(i-1))/dx**2d0
!     f0xt=(ft0(i+1)-ft0(i-1))/(2d0*dx)
!     zeta0x = (zeta0(i+1)-zeta0(i-1))/dx*0.5d0
!
! w10(i)=-u10x*(d1+zeta0(i)-f0(i))+ft0(i)+u10(i)*f0x
! ! w10x(i)=-u10xx*(d1+zeta0(i)-f0(i))-u10x*(zeta0x-f0x)+f0xt+u10x*f0x+u10(i)*f0xx
! end do
!
! ! do i= 1,imax-1
! ! w10x(i)=(w10(i+1)-w10(i-1))/dx*0.5d0
! ! end do
!
! do i= 1,imax-1
!    u10t=(u10(0)-u100(0))/(dt)
!    zeta0x =  (zeta0(1)-zeta0(-1))/dx*0.5d0
!    u10u10x = (u10(0)*(-u10(2)+8.0d0*u10(1)   &
!                  -8.0d0*u10(-1)+u10(-2))                 &
!                  +nva*dabs(u10(0))*(u10(2)-4.0d0*u10(1)+6.0d0*u10(0)  &
!                               -4.0d0*u10(-1)+u10(-2)))/(12.0d0*dx)
!
!    w10w10x = (w10(0)*(-w10(2)+8.0d0*w10(1)   &
!                     -8.0d0*w10(-1)+w10(-2))                 &
!                     +nva*dabs(w10(0))*(w10(2)-4.0d0*w10(1)+6.0d0*w10(0)  &
!                     -4.0d0*w10(-1)+w10(-2)))/(12.0d0*dx)
!
!
! px(i)=px(0)&
!          +den1*(u10t+u10u10x+w10w10x+zeta0x/(Fd**2d0))
!           ! +den1*(u10u10x+w10w10x+zeta0x/(Fd**2d0))
!
!         u10t=(u10(i)-u100(i))/(dt)
!         zeta0x  = (zeta0(i+1)-zeta0(i-1))/dx*0.5d0
!         u10u10x = (u10(i)*(-u10(i+2)+8.0d0*u10(i+1)   &
!                     -8.0d0*u10(i-1)+u10(i-2))                 &
!                     +nva*dabs(u10(i))*(u10(i+2)-4.0d0*u10(i+1)+6.0d0*u10(i)  &
!                                  -4.0d0*u10(i-1)+u10(i-2)))/(12.0d0*dx)
!
!         w10w10x = (w10(i)*(-w10(i+2)+8.0d0*w10(i+1)   &
!                    -8.0d0*w10(i-1)+w10(i-2))                 &
!                      +nva*dabs(w10(i))*(w10(i+2)-4.0d0*w10(i+1)+6.0d0*w10(i)  &
!                                   -4.0d0*w10(i-1)+w10(i-2)))/(12.0d0*dx)
!
! px(i)=px(i)  &
!         -den1*(u10t+u10u10x+w10w10x+zeta0x/(Fd**2d0))
!         ! -den1*(u10u10x+w10w10x+zeta0x/(Fd**2d0))
!
!
!
! end do

!


!----- calculation of u1,u2 -------------------------------------------

!----- u1,u2 on the upstream boundary ---------------------------------
      u1(0)=1.0d0
      u1(-1)=1.0d0
      u2(0)=1.0d0
      u2(-1)=1.0d0
      ! u1(0)=2.0d0*u1(1)-u1(2)
      ! u1(-1)=2.0d0*u1(0)-u1(1)
      ! u2(0)=2.0d0*u2(1)-u2(2)
      ! u2(-1)=2.0d0*u2(0)-u2(1)

!----- u1,u2 on the downstream boundary -------------------------------
      u1(imax)=u10(imax)-dt*u10(imax)*(u10(imax)-u10(imax-1))/dx         ! Sommerfeld radiation condition?
      u1(imax+1)=u10(imax+1)-dt*u10(imax+1)*(u10(imax+1)-u10(imax))/dx
      u2(imax)=u20(imax)-dt*u20(imax)*(u20(imax)-u20(imax-1))/dx
      u2(imax+1)=u20(imax+1)-dt*u20(imax+1)*(u20(imax+1)-u20(imax))/dx

      ! u1(imax)=1d0     ! Sommerfeld radiation condition?
      ! u1(imax+1)=1d0
      ! u2(imax)=1d0
      ! u2(imax+1)=1d0

      ! u1(imax)=0d0     ! Sommerfeld radiation condition?
      ! u1(imax+1)=0d0
      ! u2(imax)=0d0
      ! u2(imax+1)=0d0


!----- the right side of the equation ---------------------------------

!------------------------solve of u1-------------------------------------------
do i=1,imax-1


          u10x=    (u10(i+1)-u10(i-1))/dx*0.5d0
          u10xx=   (u10(i+1)-2.0d0*u10(i)+u10(i-1))/dx**2d0
          u10xxx=  (u10(i+2)-2.0d0*u10(i+1)  &
                  +2.0d0*u10(i-1)-u10(i-2))/(dx**3d0)*0.5d0

         u10t=  (u10(i)-u100(i))/dt
         u10xt=  ((u10(i+1)-u10(i-1))-(u100(i+1)-u100(i-1)))/(dx*dt)*0.5d0


          u20x =   (u20(i+1)-u20(i-1))/dx*0.5d0
          u20xx =  (u20(i+1)-2.0d0*u20(i)+u20(i-1))/dx**2d0
          u20xxx = (u20(i+2)-2.0d0*u20(i+1)   &
                  +2.0d0*u20(i-1)-u20(i-2))/(dx**3d0)*0.5d0

          u10u10x = (u10(i)*(-u10(i+2)+8.0d0*u10(i+1)   &
                      -8.0d0*u10(i-1)+u10(i-2))                 &
                      +nva*dabs(u10(i))*(u10(i+2)-4.0d0*u10(i+1)+6.0d0*u10(i)  &
                                   -4.0d0*u10(i-1)+u10(i-2)))/(12.0d0*dx)

          f0x = (f0(i+1)-f0(i-1))/dx*0.5d0
          f0xx = (f0(i+1)-2.0d0*f0(i)+f0(i-1))/dx**2d0
          f0xxx = (f0(i+2)-2.0d0*f0(i+1)   &
                 +2.0d0*f0(i-1)-f0(i-2))/(dx**3d0)*0.5d0


        ft0 =(f0(i)-f00(i))/dt
        f0xt=((f0(i+1)-f0(i-1))-(f00(i+1)-f00(i-1)))/(dt*dx)*0.5d0
        f0tt=(f0(i)-2d0*f00(i)+f000(i))/(dt*dt)
        f0xxt=((f0(i+1)-f00(i+1))-2d0*(f0(i)-f00(i))+(f0(i-1)-f00(i-1)))/(dx*dx*dt)
        f0xtt=((f0(i+1)-f0(i-1))-2d0*(f00(i+1)-f00(i-1))+(f000(i+1)-f000(i-1)))/(2d0*dt*dt*dx)

  zeta0x = (zeta0(i+1)-zeta0(i-1))/dx*0.5d0
  zeta0xx = (zeta0(i+1)-2.0d0*zeta0(i)+zeta0(i-1))/dx**2d0
  zeta0xxx = (zeta0(i+2)-2.0d0*zeta0(i+1)   &
         +2.0d0*zeta0(i-1)-zeta0(i-2))/(2d0*dx*dx*dx)

  alpha1(i)=(1.0d0+0.5d0*eta10(i)*f0xx+zeta0x*f0x)
  alpha2(i)=eta10(i)*(f0x-zeta0x)
  alpha3(i)=-(1d0/3d0)*eta10(i)**2d0

  AA(i)=-0.5d0*alpha2(i)/(dt*dx)+alpha3(i)/(dt*dx*dx)
  BB(i)=alpha1(i)/dt-2.0d0*alpha3(i)/(dt*dx*dx)
  CC(i)=0.5d0*alpha2(i)/(dt*dx)+alpha3(i)/(dt*dx*dx)




  b1(i)= (1d0/3d0)*(-u10x*u10xx+u10(i)*u10xxx)*eta10(i)**2d0                &
      + eta10(i)*(zeta0x-0.5d0*f0x)*(u10(i)*u10xx-u10x**2d0)         &
      -0.5d0*eta10(i)*((u10(i)*u10xx+u10x**2d0)*f0x+3d0*u10u10x*f0xx &
                        +f0xxx*u10(i)**2d0 +2d0*u10x*f0xt +2d0*u10(i)*f0xxt +f0xtt  )  &
     -zeta0x*(u10u10x*f0x+f0xx*u10(i)**2d0+2d0*u10(i)*f0xt+f0tt)&
     -u10u10x-(1d0/Fd**2d0)*zeta0x                         &
    +alpha1(i)*u10(i)/dt                                   &
    +alpha2(i)*(u10(i+1)-u10(i-1))/(dt*dx)*0.5d0           &
    +alpha3(i)*(u10(i+1)-2.0d0*u10(i)+u10(i-1))/(dt*dx*dx)
    ! -px(i)/den1

    G1(i)=u10xt+u10(i)*u10xx-u10x**2d0
    F1(i)=-((u10t+u10u10x)*f0x+f0xx*u10(i)**2d0+2d0*u10(i)*f0xt+f0tt)
    K1(i)=(1d0/3d0)*G1(i)*eta10(i)**3d0+0.5d0*eta10(i)**2d0*F1(i)


end do

b1(1)=b1(1)-AA(1)*u1(0)
b1(imax-1)=b1(imax-1)-CC(imax-1)*u1(imax)

do i=1,imax-1
  d(i)=BB(i)
end do

do i=1,imax-2
  du(i)=CC(i)
end do

do i=2,imax-1
  dl(i-1)=AA(i)
end do




!----- solve the u1 equation (tridiagonal matrix) ------------------------
call dgtsv(imax-1,1,dl,d,du,b1,imax-1,info)

do i=1,imax-1
  u1(i)=b1(i)
end do

!-----------------------------solve of u2---------------------------------------

do i=1,imax-1

  u10x=    (u10(i+1)-u10(i-1))/dx*0.5d0
  u10xx=   (u10(i+1)-2.0d0*u10(i)+u10(i-1))/dx**2d0
  u10xxx=  (u10(i+2)-2.0d0*u10(i+1)  &
          +2.0d0*u10(i-1)-u10(i-2))/(dx**3d0)*0.5d0
  u10t=  (u1(i)-u10(i))/dt
  u10xt=  ((u1(i+1)-u1(i-1))-(u10(i+1)-u10(i-1)))/(dx*dt)*0.5d0
  u10xxt=   ((u1(i+1)-2.0d0*u1(i)+u1(i-1))-(u10(i+1)-2.0d0*u10(i)+u10(i-1)))/(dx*dx*dt)

  u20x =   (u20(i+1)-u20(i-1))/dx*0.5d0
  u20xx =  (u20(i+1)-2.0d0*u20(i)+u20(i-1))/dx**2d0
  u20xxx = (u20(i+2)-2.0d0*u20(i+1)   &
          +2.0d0*u20(i-1)-u20(i-2))/(dx**3d0)*0.5d0

  u20t=  (u20(i)-u200(i))/dt
  u20xt=  ((u20(i+1)-u20(i-1))-(u200(i+1)-u200(i-1)))/(dx*dt)*0.5d0


  f0x = (f0(i+1)-f0(i-1))/dx*0.5d0
  f0xx = (f0(i+1)-2.0d0*f0(i)+f0(i-1))/dx**2d0
  f0xxx = (f0(i+2)-2.0d0*f0(i+1)   &
         +2.0d0*f0(i-1)-f0(i-2))/(dx**3d0)*0.5d0

 ft0 =(f0(i)-f00(i))/dt
 f0xt=((f0(i+1)-f0(i-1))-(f00(i+1)-f00(i-1)))/(dt*dx)*0.5d0
 f0tt=(f0(i)-2d0*f00(i)+f000(i))/(dt*dt)
 f0xxt=((f0(i+1)-f00(i+1))-2d0*(f0(i)-f00(i))+(f0(i-1)-f00(i-1)))/(dx*dx*dt)
 f0xtt=((f0(i+1)-f0(i-1))-2d0*(f00(i+1)-f00(i-1))+(f000(i+1)-f000(i-1)))/(2d0*dt*dt*dx)


  zeta0x = (zeta0(i+1)-zeta0(i-1))/dx*0.5d0
  zeta0xx = (zeta0(i+1)-2.0d0*zeta0(i)+zeta0(i-1))/dx**2d0
  zeta0xxx = (zeta0(i+2)-2.0d0*zeta0(i+1)   &
         +2.0d0*zeta0(i-1)-zeta0(i-2))/(2d0*dx*dx*dx)


  u10u10x = (u10(i)*(-u10(i+2)+8.0d0*u10(i+1)   &
      -8.0d0*u10(i-1)+u10(i-2))                 &
      +nva*dabs(u10(i))*(u10(i+2)-4.0d0*u10(i+1)+6.0d0*u10(i)  &
                   -4.0d0*u10(i-1)+u10(i-2)))/(12.0d0*dx)

  u20u20x = (u20(i)*(-u20(i+2)+8.0d0*u20(i+1)   &
      -8.0d0*u20(i-1)+u20(i-2))                 &
      +nva*dabs(u20(i))*(u20(i+2)-4.0d0*u20(i+1)+6.0d0*u20(i)  &
                    -4.0d0*u20(i-1)+u20(i-2)))/(12.0d0*dx)


  beta4(i)=(1d0+0.5d0*eta20(i)*hxx(i)+f0x*hx(i))
  beta5(i)=-(f0x-hx(i))*eta20(i)
  beta6(i)=-(1d0/3d0)*eta20(i)**2d0


  II(i)=-0.5d0*beta5(i)/(dt*dx)+beta6(i)/(dt*dx*dx)
  KK(i)=beta4(i)/dt-2.0d0*beta6(i)/(dt*dx*dx)
  MM(i)=0.5d0*beta5(i)/(dt*dx)+beta6(i)/(dt*dx*dx)


    b1(i)=(1d0/3d0)*(-u20x*u20xx+u20(i)*u20xxx)*eta20(i)**2d0           &
        +eta20(i)*(f0x-0.5d0*hx(i))*(u20(i)*u20xx-u20x**2d0)           &
        -0.5d0*eta20(i)*((u20(i)*u20xx+u20x**2d0)*hx(i)+3d0*u20u20x*hxx(i)+u20(i)**2d0*hxxx(i))  &
        -f0x*(u20u20x*hx(i)+hxx(i)*u20(i)**2d0)                     &
        +den1*(1d0/Fd**2d0)*(f0x-zeta0x)                              &
        +den1*(eta10(i)*(zeta0x-f0x)*(u10xt+u10(i)*u10xx-u10x**2d0)          &
                +0.5d0*(u10xxt-u10x*u10xx+u10(i)*u10xxx)*eta10(i)**2d0          &
                -eta10(i)*((u10xt+u10(i)*u10xx+u10x**2d0)*f0x+(u10t+3d0*u10u10x)*f0xx &
                          +f0xxx*u10(i)**2d0 +2d0*u10x*f0xt+2d0*u10(i)*f0xxt+f0xtt )&
              -(zeta0x-f0x)*((u10t+u10u10x)*f0x+f0xx*u10(i)**2d0+2d0*u10(i)*f0xt+f0tt)   )&
          -u20u20x-(1d0/Fd**2d0)*f0x                         &
          +beta4(i)*u20(i)/dt                                   &
          +beta5(i)*(u20(i+1)-u20(i-1))/(dt*dx)*0.5d0           &
          +beta6(i)*(u20(i+1)-2.0d0*u20(i)+u20(i-1))/(dt*dx*dx)
          ! -px(i)

          G2(i)=u20xt+u20(i)*u20xx-u20x**2d0
          F2(i)=-((u20t+u20u20x)*hx(i)+hxx(i)*u20(i)**2d0)
          p2(i)=-den1*((f0(i)-zeta0(i)-d1)/(Fd**2d0)+(0.5d0*G1(i)*eta10(i)**2d0+F1(i)*eta10(i)))
          K2(i)=(1d0/3d0)*G2(i)*eta20(i)**3d0+0.5d0*eta20(i)**2d0*F2(i)


end do

b1(1)=b1(1)-II(1)*u2(0)
b1(imax-1)=b1(imax-1)-MM(imax-1)*u2(imax)

do i=1,imax-1
  d(i)=KK(i)
end do

do i=1,imax-2
  du(i)=MM(i)
end do

do i=2,imax-1
  dl(i-1)=II(i)
end do



!----- solve the equation (tridiagonal matrix) ------------------------
call dgtsv(imax-1,1,dl,d,du,b1,imax-1,info)

do i=1,imax-1
  u2(i)=b1(i)
end do

!----- NaN check ------------------------------------------------------

  2002 format(' ','NaN output at iten=',I9)
  do i=-1,imax+1
  if(isnan(f(i)) .eqv. .true.) then
    write(6,*) 'f is NaN !'
    write(6,2002) iten

    stop   !do while exit

  end if
 end do
 do i=-1,imax+1

  if(isnan(zeta(i)) .eqv. .true.) then
    write(6,*) 'zeta is NaN !'
    write(6,2002) iten

    stop   !do while exit

  end if
end do

do j=-1,imax+1
  if(isnan(u1(j)) .eqv. .true.) then
    write(6,*) 'u1 is NaN !'

    write(6,2002) iten
    itime = int(iten*dt)
    write(citen,'(I7.7)') itime

     open(4,file=citen//'_error.dat')
       write(4,*) "#x f zeta u1 u2"
       fmax=0.0
       do i=-1,imax+1
           write(4,4000) x(i),f0(i),zeta0(i),u10(i),u20(i)
           if(fmax<f(i)) fmax=f(i)
       end do
  4000 format(F10.3,1X,E16.8e3,1x,E16.8e3,1x,E16.8e3,1x,E16.8e3)
       close(4)
  stop   !do while exit

  end if
end do
!
! itime = int(iten)
! write(citen,'(I7.7)') itime
!
!
! !----- writing of f,zeta,u1,u2 ---------------------------------------------
!
!  open(4,file=citen//'_FNM.dat')
!    write(4,*) "#x f zeta u1 u2 "
!    fmax=0.0
!    do i=-1,imax+1
!        write(4,4500) x(i),f(i),zeta(i),u1(i),u2(i)
!        if(fmax<f(i)) fmax=f(i)
!    end do
!
!  4100 format(F10.3,4e27.16E3)
!
!   open(4,file=citen//'_FNM_inexplcit.dat')
!     write(4,*) "#x f zeta u1 u2"
!     fmax=0.0
!     do i=-1,imax+1
!         write(4,4100) x(i),f0(i),zeta0(i),u10(i),u20(i)
!         if(fmax<f(i)) fmax=f(i)
!     end do
!
!     close(4)
!
!     open(10,file=citen//'_GF_inexplcit.dat')
!       write(10,*) "#x G1 G2 F1 F2 "
!       do i=-1,imax+1
!           write(10,4100) x(i),G1(i),G2(i),F1(i),F2(i)
!
!       end do
!
!       open(15,file=citen//'_Kp_inexplcit.dat')
!         write(15,*) "#x K1 K2 p2 "
!         do i=-1,imax+1
!             write(15,4100) x(i),K1(i),K2(i),p2(i)
!         end do
!
!       close(15)

!----- conservation for next time -------------------------------------
       do i=-1,imax+1
          ! zeta0(i)=zeta(i)
          ! zeta00(i)=zeta0(i)
          ! zeta000(i)=zeta00(i)
          ! f0(i)=f(i)
          ! f00(i)=f0(i)
          ! f000(i)=f00(i)
          ! ft0(i)=ft(i)
          ! ft00(i)=ft0(i)
          ! u10(i)=u1(i)
          ! u20(i)=u2(i)
          ! u100(i)=u10(i)
          ! u200(i)=u20(i)
          ! eta10(i)=eta1(i)
          ! eta20(i)=eta2(i)


          zeta000(i)=zeta00(i)
          zeta00(i)=zeta0(i)
          zeta0(i)=zeta(i)

          f000(i)=f00(i)
          f00(i)=f0(i)
          f0(i)=f(i)

          ! ft00(i)=ft0(i)
          ! ft0(i)=ft(i)

          u100(i)=u10(i)
          u200(i)=u20(i)
          u10(i)=u1(i)
          u20(i)=u2(i)

          eta10(i)=eta1(i)
          eta20(i)=eta2(i)
       end do

!-----------------------------------------------------------------------
! if(mod(iten,1000)==0) then
!   write(citen,'(I7.7)') iten
!
!     open(20,file=citen//'_xfzeta.dat')
!       write(20,*) "#x f zeta "
!       fmax=0.0
!       do i=-1,imax+1
!           write(20,4100) x(i),f(i),zeta(i),u1(i),u2(i)
!           if(fmax<f(i)) fmax=f(i)
!       end do
!  4100 format(F10.3,1X,E16.8e3,1x,E16.8e3,1x,E16.8e3,1x,E16.8e3)
!       close(20)
!
!
!   end if



!----- G2 on the upstream boundary ---------------------------------
G2(0)=0.0d0
G2(-1)=0.0d0
! G2(0)=2.0d0*G2(1)-G2(2)
! G2(-1)=2.0d0*G2(0)-G2(1)
!----- G2 on the downstream boundary -----------------------------------
G2(imax)=2.0d0*G2(imax-1)-G2(imax-2)
G2(imax+1)=2.0d0*G2(imax)-G2(imax-1)
!-----------------------------------------------------------------------

!----- F2 on the upstream boundary ---------------------------------
F2(0)=0.0d0
F2(-1)=0.0d0
! F2(0)=2.0d0*F2(1)-F2(2)
! F2(-1)=2.0d0*F2(0)-F2(1)
!----- F2 on the downstream boundary -----------------------------------
F2(imax)=2.0d0*F2(imax-1)-F2(imax-2)
F2(imax+1)=2.0d0*F2(imax)-F2(imax-1)
!-----------------------------------------------------------------------

!----- K2 on the upstream boundary ---------------------------------
K2(0)=0.0d0
K2(-1)=0.0d0
! K2(0)=2.0d0*K2(1)-K2(2)
! K2(-1)=2.0d0*K2(0)-K2(1)
!----- K2 on the downstream boundary -----------------------------------
K2(imax)=2.0d0*K2(imax-1)-K2(imax-2)
K2(imax+1)=2.0d0*K2(imax)-K2(imax-1)
!-----------------------------------------------------------------------

!----- p2 on the upstream boundary ---------------------------------
! p2(0)=0.0d0
! p2(-1)=0.0d0
p2(0)=2.0d0*p2(1)-p2(2)
p2(-1)=2.0d0*p2(0)-p2(1)
!----- p2 on the downstream boundary -----------------------------------
p2(imax)=2.0d0*p2(imax-1)-p2(imax-2)
p2(imax+1)=2.0d0*p2(imax)-p2(imax-1)
!-----------------------------------------------------------------------

! !----- judgement of iten ----------------------------------------------
if(mod(iten,itewri)==0) then

   itime = int(iten*dt)
   write(citen,'(I7.7)') itime


!----- writing of f,zeta,u1,u2 ---------------------------------------------

    open(4,file=citen//'_FNM.dat')
      write(4,*) "#x f zeta u1 u2 "
      fmax=0.0
      do i=-1,imax+1
          write(4,4500) x(i),f(i),zeta(i),u1(i),u2(i)
          if(fmax<f(i)) fmax=f(i)
      end do
 4500 format(F10.3,4e27.16E3)
      close(4)

      call SYSTEM_CLOCK(iet,irt)
      write(6,5000) iten,itime,(iet-ist)/irt
      ist=iet
 5000 format(' ','iten=',i8,1X,'itime=',i6,1X,'time=',i8,'s')
      write(6,*) 'fmax=', fmax

    end if


!   write(citen,'(I7.7)') iten
!
!     open(20,file=citen//'_xfzeta.dat')
!       write(20,*) "#x f zeta "
!       fmax=0.0
!       do i=-1,imax+1
!           write(20,4100) x(i),f(i),zeta(i),u1(i),u2(i)
!           if(fmax<f(i)) fmax=f(i)
!       end do
!  4100 format(F10.3,1X,E16.8e3,1x,E16.8e3,1x,E16.8e3,1x,E16.8e3)
!       close(20)
!
!



    iten=iten+1

  end do !end do while iten
!----- finish judgement -----------------------------------------------

!----- writing the end message ----------------------------------------

      write(6,2000) iten
 2000 format(' ','normal end at iten=',I9)




!----- finish output --------------------------------------------------


      open(70,file='finish_FNM.dat')
      write(70,*) "#x f0 zeta0 u10 u20"
      do i=-1,imax+1
          write(70,4000) x(i),f0(i),zeta0(i),u10(i),u20(i)
      end do
      close(70)

!----- finish time ----------------------------------------------------
      call system_clock(iet2,irt2,itmax)
      if(iet2<ist2) then
        its=(itmax-ist2+iet2)/irt2
      else
        its=(iet2-ist2)/irt2
      end if
      ith=its/3600
      itm=its/60-ith*60
      its=mod(its,60)
      open(60,file='finish.dat')
      write(6,5001) ith,itm,its
      write(60,5001) ith,itm,its
 5001 format(' ','all_time=',i2,'h',1X,i2,'m',1X,i2,'s')
      write(60,*) 'iten =', iten
      close(60)

!----------------------------------------------------------------------
      stop
      end
