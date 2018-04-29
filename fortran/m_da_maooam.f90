module m_da_maooam
  use m_mt, only : eigen, inv
  use mod_optimization, only : run_lbfgs, destroy_lbfgs
  implicit none

  private 
  public :: etkf
  public :: etkf2
  public :: tdvar
  public :: inctdvar

  !integer,parameter :: nx = 20
  !integer,parameter :: nyo = 20
  !integer,parameter :: nx = 36
  !integer,parameter :: nyo = 36
  !integer,parameter :: kitermax = 500
  real(8),parameter :: epsG   = 1.d-2

contains

subroutine inctdvar( nx, nyo, xb, B, lyo, yo, erro, xa, &
                     kitermax_out, kitermax_in )
  implicit none
! passed args
  integer,   intent(in   ) :: nx
  integer,   intent(in   ) :: nyo
  real(8),intent(in   ) :: xb(nx)
  real(8),intent(in   ) :: B(nx,nx)
  logical,   intent(in   ) :: lyo(nyo)
  real(8),intent(in   ) :: yo(nyo)
  real(8),intent(in   ) :: erro(nyo)
  real(8),intent(  out) :: xa(nx)
  integer,intent(in) :: kitermax_out
  integer,intent(in) :: kitermax_in
! local vars
  real(8) :: invB(nx,nx), invR(nyo)
  real(8) :: x(nx)
  real(8) :: Jc, Jo, Jb, Jc_old
  real(8) :: dJc(nx), dJo(nx), dJb(nx)

  real(8) :: yg(nyo)  ! H(xg)
  real(8) :: dyg(nyo) ! yobs - H(xg)
  real(8) :: xg(nx)   ! reference state
  real(8) :: dxg(nx)  ! xb -xg
  real(8) :: dx(nx)   
  real(8) :: Hdx(nx)

  integer :: kiter_out, kiter_in
  integer :: m, n

  integer :: ierr
  integer :: i, j
  logical :: lfound

  kiter_out = kitermax_out
  kiter_in  = kitermax_in

! get B^-1, and R^-1
  Call INV( nx, B, invB )

  Write(6,*) "B ="
  Do i = 1, nx
     Write(6,*)   ( B(i,j), j = 1, nx )
  Enddo

  Write(6,*) "inv(B) ="
  Do i = 1, nx
     Write(6,*)   ( invB(i,j), j = 1, nx )
  Enddo

  !If ( ierr /= 0 ) Stop "[lorenz63_3dvar] error: fail to get inverse of B"
  Do i = 1, nyo
     If ( lyo(i) ) Then
        invR(i) = 1.d0/( erro(i)*erro(i) )
     Endif
  Enddo


  xg(:) = xb(:)
  outer_loop: do m = 1, kiter_out

     ! 1. recalculate nonlinear h(x)


     ! 1.1 dxg = xb - xg
     do i = 1, nx
        dxg(i) = xb(i) - xg(i)
     enddo

     print*, "-----------------------------------------------"
     print*, "outerloop", m
     print*, "dxg=", dxg, xg
     !print*, "------------------"

     ! 1.2 dyg = yobs - h(xg)
     do i = 1, nyo
        !---------------------
        ! call hx(i,xg,yg(i))
        yg(i) = xg(i)
        !---------------------
        dyg(i) = yo(i) - yg(i) 
     enddo

     ! 2. Minimization loop
     !dx = dxg
     dx = 0.d0
     lfound = .false.

     inner_loop: do n = 1, kiter_in
        ! Calcualte Jb = 0.5 * (dx -dxg)* B^-1 * (dx-dxg), where dx = x-xg
        Jb = 0.0d0
        do i = 1, nx
           do j = 1, nx
              Jb = Jb + 0.5d0*( dx(i)-dxg(i) )*( dx(j)-dxg(j) )*invB(i,j) 
           enddo
        enddo
    
        ! Calculate Jo = 0.5 * ( dyg - Hdx ) * R^-1 * ( dyg - Hdx )
        Jo = 0.0d0
        do i = 1, nyo
           if ( lyo(i) ) then
              !-----------------------
              ! call hx_d(i,xb,xb_d, yb_d)
              !Hdx(i) = dx(i)
              !-----------------------
              Jo = Jo + 0.5d0*( dyg(i) - dx(i) )*( dyg(i) - dx(i) )*invR(i)  ! replace dx(i) with Hdx(i)
           endif
        enddo

        ! Calculate total cost function Jc = Jb + Jo
        Jc = Jb + Jo

        ! Calculate Gradient dJb = B^-1 * ( dx -dxg )
        dJb = 0.0d0
        do i = 1, nx
           do j = 1, nx
              dJb(i) = dJb(i) + invB(i,j) * ( dx(j)-dxg(j) ) 
           enddo
        enddo

        ! Calculate Gradient dJo = -H^T * R^-1 * ( dyg- Hdx )
        dJo = 0.0d0
        do i = 1, nx
           if ( lyo(i) ) then
              !-----------------------
              ! call hx_d(i,xb,xb_d, yb_d)
              !-----------------------
              dJo(i) = - invR(i) * ( dyg(i)-dx(i) ) ! replace dx(i) with Hdx(i)
           endif
        enddo

        ! Calculate Total gradient dJc = dJo + dJb
        dJc = dJo + dJb

        !print*, "out, in, J=", m, n, Jc, Jo, Jb

        !if (n>1 .and. abs(Jc_old-Jc)<1.0d-4*abs(Jc)) then
        !   print*, "no change of J: Jc_old, Jc=", Jc_old, Jc
        !   lfound = .true.
        !   call destroy_lbfgs()
        !   exit inner_loop
        !else
        !   Jc_old = Jc
        !endif

        call run_lbfgs( nx, dx, Jc, dJc, epsG, ierr )
        if ( ierr < 0 ) Stop "[lorenz63_inc3dvar] error:fail to run lbfgs."
        if ( ierr == 0 ) then
           lfound = .true.
           exit inner_loop
        endif

    enddo inner_loop

    if (.not.lfound) then
       write(*,*) "reach the max inner loop:", kitermax_in
       call destroy_lbfgs()
    endif

    print*, "before: xg=", xg
    xg(:) = xg(:) + dx(:)
    print*, "after: xg=", xg
    print*, "dx=", dx

  enddo outer_loop

  xa(:) = xg(:)
  print*, "finished"

endsubroutine



subroutine tdvar( nx, nyo, xb, B, lyo, yo, erro, xa, kitermax )
  implicit none
! passed args
  integer, intent(in) :: nx
  integer, intent(in) :: nyo
  real(8),intent(in   ) :: xb(nx)
  real(8),intent(in   ) :: B(nx,nx)
  logical,   intent(in   ) :: lyo(nyo)
  real(8),intent(in   ) :: yo(nyo)
  real(8),intent(in   ) :: erro(nyo)
  real(8),intent(  out) :: xa(nx)
  integer,intent(in) :: kitermax 
! local vars
  real(8) :: invB(nx,nx), invR(nyo)
  real(8) :: x(nx)
  real(8) :: Jc, Jo, Jb
  real(8) :: dJc(nx), dJo(nx), dJb(nx)

  integer :: kiter

  integer :: ierr
  integer :: i, j

! get B^-1, and R^-1
  Call INV( nx, B, invB )

  !Write(6,*) "B ="
  !Do i = 1, nx
  !   Write(6,*)   ( B(i,j), j = 1, nx )
  !Enddo

  !Write(6,*) "inv(B) ="
  !Do i = 1, nx
  !   Write(6,*)   ( invB(i,j), j = 1, nx )
  !Enddo

  !If ( ierr /= 0 ) Stop "[lorenz63_3dvar] error: fail to get inverse of B"
  Do i = 1, nyo
     If ( lyo(i) ) Then
        invR(i) = 1.d0/( erro(i)*erro(i) )
     Endif
  Enddo

! Set xa=xb at the inital step
  x = xb

  kiter = 0
! Minimization loop
  lp_mini: Do 

    ! Calcualte Jb = 0.5 * ( x - xb) * B^-1 * ( x - xb )
    Jb = 0.0d0
    Do i = 1, nx
       Do j = 1, nx
          Jb = Jb + 0.5d0*( x(i)-xb(i) )*( x(j)-xb(j) )*invB(i,j) 
       Enddo
    Enddo
    ! Calculate Jo = 0.5 * ( h[x] - yo ) * R^-1 * ( h[x] - yo )
    Jo = 0.0d0
    Do i = 1, nyo
       If ( lyo(i) ) Then
          Jo = Jo + 0.5d0*( x(i)-yo(i) )*( x(i)-yo(i) )*invR(i)
       Endif
    Enddo
    ! Calculate total cost function Jc = Jb + Jo
    Jc = Jb + Jo

    ! Calculate Gradient dJb = B^-1 * ( x - xb )
    dJb = 0.0d0
    Do i = 1, nx
       Do j = 1, nx
          dJb(i) = dJb(i) + invB(i,j) * ( x(j)-xb(j) ) 
       Enddo
    Enddo
    ! Calculate Gradient dJo = H^T * R^-1 * ( h[x] - yo )
    dJo = 0.0d0
    Do i = 1, nx
       If ( lyo(i) ) Then
          dJo(i) = invR(i) * ( x(i) - yo(i) ) 
       Endif
    Enddo
    ! Calculate Total gradient dJc = dJo + dJb
    dJc = dJo + dJb

    ! 
    Call run_lbfgs( nx, x, Jc, dJc, epsG, ierr )
    If ( ierr < 0 ) Stop "[lorenz63_3dvar] error:fail to run lbfgs."
    kiter = kiter + 1
    If ( kiter > kitermax ) Then
       Write(6,*) "[lorenz63_3dvar] warning: maximum iter reached. "
       call destroy_lbfgs()
       exit lp_mini
    Endif
    If ( ierr == 0 ) exit lp_mini

  Enddo lp_mini

  xa = x

endsubroutine


subroutine etkf( nx, nyo, nn, xb, lyo, yo, erro, infl, xa, xam, KH )
  implicit none
! passed args
  integer,    intent(in   ) :: nx
  integer,    intent(in   ) :: nyo
  integer,    intent(in   ) :: nn
  real(8), intent(in   ) :: xb(nx,nn)
  logical,    intent(in   ) :: lyo(nyo)
  real(8), intent(in   ) :: yo(nyo)
  real(8), intent(in   ) :: erro(nyo)
  real(8), intent(in   ) :: infl
  real(8), intent(  out) :: xa(nx,nn)
  real(8), intent(  out) :: xam(nx)
  real(8), intent(  out) :: KH(nx, nx)
! local vars
  integer :: nyo_use
  real(8),allocatable :: dyb(:,:)
  real(8),allocatable :: ybm(:)
  real(8),allocatable :: erro_use(:)
  real(8),allocatable :: yo_ybm(:)
  real(8) :: xbm(nx)
  real(8) :: deltaxb(nx,nn)
  real(8) :: invSPa(nn,nn), SPa(nn,nn) ! Pa in ensemble space
  real(8),allocatable :: C(:,:)
  real(8) :: eigval(nn), eigvect(nn,nn)
  integer :: np
  real(8) :: wam(nn), Wa(nn,nn)
  real(8) :: dxb(nn)
  integer :: ierr
  integer :: i, j, k

  nyo_use = 0
! calculate mean yb and perturbation matrix Yb
  do i = 1, nyo
     if ( lyo(i) ) then
        nyo_use = nyo_use + 1
     else
        write(6,*) "skip obs: iob =", i, ", yo(iob)=", yo(i)
     endif
  enddo
  allocate( dyb(nyo_use,nn) )
  allocate( ybm(nyo_use), yo_ybm(nyo_use), erro_use(nyo_use) )

  k = 0
  do i = 1, nyo
     if ( lyo(i) ) then
        k           = k + 1
        dyb(k,:)    = xb(i,:)
        ybm(k)      = SUM(dyb(k,:))/nn
        dyb(k,:)    = dyb(k,:) - ybm(k)
        yo_ybm(k)   = yo(i) - ybm(k)
        erro_use(k) = erro(i)
        !print*, "k, erro, yo, ybm=", k, erro_use(k), yo(i), ybm(k)
        !pause
     endif
  enddo

! calculate mean xb
  do i = 1, nx
     xbm(i) = SUM(xb(i,:))/nn
  enddo
 
  allocate( C(nn,nyo_use) )

! C=(Yb)^T * R^-1
  do j = 1, nyo_use
     C(:,j) = dyb(j,:)/(erro_use(j)**2)
  enddo
! invSPa = [(k-1)I+C*Yb]
  invSPa = MATMUL( C, dyb )
  do k = 1, nn
     invSPa(k,k) = invSPa(k,k)+(nn-1)/infl
  enddo
! 
! so A^-1=Q*V^-1*Q^T
  call eigen( nn, invSPa, eigvect, eigval)
  !if ( ierr/=0 .or. np/=nn ) then
  !   write(6,*) "[warning] lorenz63_letkf: fail to find nn (+) eigenvetor"
  !   xam = xbm
  !   xa = xb
  !   return
  !endif
  do k = 1, nn
     SPa(:,k) = eigvect(:,k)/eigval(k)
  enddo
  SPa = MATMUL( SPa, TRANSPOSE(eigvect) )
  ! Wa = [(k-1)SPa]^1/2
  !    = [(k-1)A^-1]^1/2
  !    = [ Q*(k-1)V^-1*Q^T ]^1/2
  !    = [ Q*sqrt((k-1)/V)*Q^T ]
  do k = 1, nn
     Wa(:,k) = eigvect(:,k)*SQRT( (nn-1)/eigval(k) )
  enddo
  Wa = MATMUL( Wa, TRANSPOSE(eigvect) )
  ! wam = SPa * C *( yo -ybm )
  wam = MATMUL( C, yo_ybm )
  wam = MATMUL( SPa, wam )  !!!
  ! Wa = Wa + wam
  do k = 1, nn
     Wa(:,k) = Wa(:,k) + wam(:)
  enddo

  ! Calculate gain

  do i = 1, nx
     deltaxb(i, :) = xb(i, :) - xbm(i)
  !   I_m(i, i) = 1.0
  end do

  !YRY = MATMUL(dyb, TRANSPOSE(C))

  KH = MATMUL(MATMUL(MATMUL(deltaxb, SPa), C), dyb)

  do i = 1, nx
     ! xam = xbm + Xb*wam
     dxb = xb(i,:)-xbm(i)
     xam(i) = xbm(i) + SUM(dxb(:)*wam(:))

     do k = 1, nn
        xa(i,k) = xbm(i) + SUM(dxb(:)*Wa(:,k))
     enddo
  enddo

endsubroutine


subroutine etkf2( nx, nyo, nn, xb, lyo, yo, erro, infl, xa, xam )
  implicit none
! passed args
  integer, intent(in) :: nx
  integer, intent(in) :: nyo
  integer,    intent(in   ) :: nn
  real(8), intent(in   ) :: xb(nx,nn)
  logical,    intent(in   ) :: lyo(nyo)
  real(8), intent(in   ) :: yo(nyo)
  real(8), intent(in   ) :: erro(nyo)
  real(8), intent(in   ) :: infl
  real(8), intent(  out) :: xa(nx,nn)
  real(8), intent(  out) :: xam(nx)
! local vars
  integer :: nyo_use
  real(8),allocatable :: dyb(:,:)
  real(8),allocatable :: ybm(:)
  real(8),allocatable :: erro_use(:)
  real(8),allocatable :: yo_ybm(:)
  real(8) :: xbm(nx)
  real(8) :: invSPa(nn,nn), SPa(nn,nn) ! Pa in ensemble space
  real(8),allocatable :: C(:,:)
  real(8) :: eigval(nn), eigvect(nn,nn)
  integer :: np
  real(8) :: wam(nn), Wa(nn,nn)
  real(8) :: dxb(nn)
  integer :: ierr
  integer :: i, j, k

  nyo_use = 0
! calculate mean yb and perturbation matrix Yb
  do i = 1, nyo
     if ( lyo(i) ) then
        nyo_use = nyo_use + 1
     else
        write(6,*) "skip obs: iob =", i, ", yo(iob)=", yo(i)
     endif
  enddo
  allocate( dyb(nyo_use,nn) )
  allocate( ybm(nyo_use), yo_ybm(nyo_use), erro_use(nyo_use) )

  k = 0
  do i = 1, nyo
     if ( lyo(i) ) then
        k           = k + 1
        dyb(k,:)    = xb(i,:)
        ybm(k)      = SUM(dyb(k,:))/nn
        dyb(k,:)    = dyb(k,:) - ybm(k)
        yo_ybm(k)   = yo(i) - ybm(k)
        erro_use(k) = erro(i)
        !print*, "k, erro, yo, ybm=", k, erro_use(k), yo(i), ybm(k)
        !pause
     endif
  enddo

! calculate mean xb
  do i = 1, nx
     xbm(i) = SUM(xb(i,:))/nn
  enddo
 
  allocate( C(nn,nyo_use) )

! C=(Yb)^T * R^-1
  do j = 1, nyo_use
     C(:,j) = dyb(j,:)/(erro_use(j)**2)
  enddo
! invSPa = [(k-1)I+C*Yb]
  invSPa = MATMUL( C, dyb )
  do k = 1, nn
     invSPa(k,k) = invSPa(k,k)+(nn-1)/infl
  enddo
! let A=invSPa, where A=Q*V*Q^T
! so A^-1=Q*V^-1*Q^T
  call eigen( nn, invSPa, eigvect, eigval)
  !if ( ierr/=0 .or. np/=nn ) then
  !   write(6,*) "[warning] lorenz63_letkf: fail to find nn (+) eigenvetor"
  !   xam = xbm
  !   xa = xb
  !   return
  !endif
  do k = 1, nn
     SPa(:,k) = eigvect(:,k)/eigval(k)
  enddo
  SPa = MATMUL( SPa, TRANSPOSE(eigvect) )
  ! Wa = [(k-1)SPa]^1/2
  !    = [(k-1)A^-1]^1/2
  !    = [ Q*(k-1)V^-1*Q^T ]^1/2
  !    = [ Q*sqrt((k-1)/V)*Q^T ]
  do k = 1, nn
     Wa(:,k) = eigvect(:,k)*SQRT( (nn-1)/eigval(k) )
  enddo
  Wa = MATMUL( Wa, TRANSPOSE(eigvect) )
  ! wam = SPa * C *( yo -ybm )
  wam = MATMUL( C, yo_ybm )
  wam = MATMUL( SPa, wam )  !!!
  ! Wa = Wa + wam
  do k = 1, nn
     Wa(:,k) = Wa(:,k) + wam(:)
  enddo


  do i = 1, nx
     ! xam = xbm + Xb*wam
     dxb = xb(i,:)-xbm(i)
     xam(i) = xbm(i) + SUM(dxb(:)*wam(:))

     do k = 1, nn
        xa(i,k) = xbm(i) + SUM(dxb(:)*Wa(:,k))
     enddo
  enddo

endsubroutine




endmodule
