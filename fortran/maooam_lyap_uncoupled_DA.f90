
!  maooam_lyap.f90
!
!> Fortran 90 implementation of the modular arbitrary-order ocean-atmosphere
!> model MAOOAM computing the Lyapunov spectrum.
!
!> @copyright
!> 2016 Lesley De Cruz, Sebastian Schubert & Jonathan Demaeyer.
!> See LICENSE.txt for license information.
!
!---------------------------------------------------------------------------!

PROGRAM maooam_lyap
  USE params, only: ndim, natm, noc, dt, tw, t_trans, t_run, writeout, rescaling_time
  USE aotensor_def, only: init_aotensor
  USE tl_ad_tensor, only: init_tltensor
  USE IC_def, only: load_IC, IC
  USE integrator, only: init_integrator,step
  USE tl_ad_integrator, only: init_tl_ad_integrator,prop_step
  USE lyap_vectors, only: lyapunov,lyapunov_atm,lyapunov_ocn,loclyap,&
      loclyap_atm,loclyap_ocn,init_lyap,init_lyap_atm,init_lyap_ocn,&
      multiply_prop,multiply_prop_atm,multiply_prop_ocn,benettin_step,&
      benettin_step_atm,benettin_step_ocn
  USE stat
  USE lyap_stat
  use m_maooam, only: get_res_maooam, init_maooam
  use m_io,        only : daio
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X, X_atm, X_ocn          !< State variable in the model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Xnew, Xnew_atm, Xnew_ocn       !< Updated state variable
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: prop_buf, prop_buf_atm, prop_buf_ocn !< Buffer for the propagator
  real(8),allocatable :: Pb(:,:), Pbatm(:,:), Pbocn(:,:)
  real(8),allocatable :: Rinv(:, :), R(:), eye(:, :), KH_atm(:, :), KH_ocn(:, :)
  REAL(KIND=8) :: t=0.D0                                !< Time variable
  REAL(KIND=8) :: t_up, gasdev
  INTEGER :: IndexBen,WRSTAT
  CHARACTER(LEN=19) :: FMTX
!  CHARACTER(LEN=3) :: sim_id
  INTEGER :: idum1
  INTEGER :: j, k, n

  integer :: nxa, nya, nxo, nyo

!  call get_res_maooam(nxa,nya,nxo,nyo)
!  call init_maooam(nxa,nya,nxo,nyo,ndim)

  idum1=-1254

  PRINT*, 'Model MAOOAM v1.0'
  PRINT*, '      - with computation of the Lyapunov spectrum'
  PRINT*, 'Loading information...'

!  CALL get_command_argument(1, sim_id)

  CALL init_aotensor!(sim_id)    ! Compute the tensors
  CALL init_tltensor
  CALL load_IC          ! Load the initial condition

  CALL init_integrator        ! Initialize the integrator
  CALL init_tl_ad_integrator  ! Initialize tangent linear integrator
  CALL init_lyap              ! Initialize Lyapunov computation
  CALL init_lyap_atm              ! Initialize Lyapunov computation
  CALL init_lyap_ocn              ! Initialize Lyapunov computation
  write(FMTX,'(A10,i3,A6)') '(F10.2,4x,',ndim,'E15.5)'
  t_up=dt/t_trans*100.D0

  IF (writeout) THEN
     OPEN(10,file='evol_field.dat')
     OPEN(11,file='lyapunov_exponents_both_3D.dat',status='replace',form='UNFORMATTED',access='DIRECT',recl=8*ndim)
     OPEN(14,file='lyapunov_exponents_atm_3D.dat',status='replace',form='UNFORMATTED',access='DIRECT',recl=8*ndim)
     OPEN(15,file='lyapunov_exponents_ocn_3D.dat',status='replace',form='UNFORMATTED',access='DIRECT',recl=8*ndim)
  END IF

  ALLOCATE(X(0:ndim),Xnew(0:ndim),X_atm(0:ndim),Xnew_atm(0:ndim),&
      X_ocn(0:ndim),Xnew_ocn(0:ndim),prop_buf(ndim,ndim),&
      prop_buf_atm(ndim,ndim),prop_buf_ocn(ndim,ndim))
  allocate(Pb(ndim,ndim))
  allocate(Pbatm(2*natm,2*natm))
  allocate(Pbocn(2*noc,2*noc))
  allocate(R(ndim))
  allocate(Rinv(ndim, ndim))
  allocate(eye(ndim, ndim))
  allocate(KH_atm(2*natm, 2*natm))
  allocate(KH_ocn(2*noc, 2*noc))


  ! load Pb
  do k = 1, ndim
     read(daio%B,"(10000(D24.17,1x))") (Pb(k,n),n=1,ndim)
  enddo
  Pbatm(1:2*natm,1:2*natm) = Pb(1:2*natm,1:2*natm)
  Pbocn(1:2*noc,1:2*noc) = Pb(1+2*natm:ndim,1+2*natm:ndim)

  ! load R
  Rinv = 0.D0
  eye = 0.D0
  do n = 1, ndim
     read(daio%R,"(10000(D24.17,1x))") R(n)
     Rinv(n, n) = 1/R(n)
     eye(n, n) = 1.D0
     print*, "n, R(n)=", n, Rinv(n, n)
  enddo

  KH_atm = Rinv(1:2*natm,1:2*natm) + inv(Pbatm)
  KH_atm = inv(KH_atm)
  KH_atm = MATMUL(KH_atm, Rinv(1:2*natm,1:2*natm))

  KH_ocn = Rinv(1+2*natm:ndim,1+2*natm:ndim) + inv(Pbocn)
  KH_ocn = inv(KH_ocn)
  KH_ocn = MATMUL(KH_ocn, Rinv(1+2*natm:ndim,1+2*natm:ndim))

  X = IC

  X_atm = IC
  X_ocn = IC

  PRINT*, 'Starting the transient time evolution... t_trans = ',t_trans

  DO WHILE (t<t_trans)
     CALL step(X,t,dt,Xnew)

     t = t - dt
     X_atm(2*natm+1:ndim) = X(2*natm+1:ndim)
     CALL step(X_atm, t, dt, Xnew_atm)
     X_atm = Xnew_atm

     t = t - dt
     X_ocn(1:2*natm) = X(1:2*natm)
     CALL step(X_ocn, t, dt, Xnew_ocn)
     X_ocn = Xnew_ocn

     X=Xnew

     IF (mod(t/t_trans*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_trans*100.D0,char(13)
  END DO

  PRINT*, 'Starting the time evolution... t_run = ',t_run

  CALL init_stat
  CALL lyap_init_stat
  t=0.D0
  IndexBen=0
  t_up=dt/t_run*100.D0

  DO WHILE (t<t_run)
     CALL prop_step(X,prop_buf,t,dt,Xnew,.false.) ! Obtains propagator prop_buf at X
     CALL multiply_prop(prop_buf) ! Multiplies prop_buf with prop

     X_atm(2*natm+1:ndim) = X(2*natm+1:ndim)
     t = t - dt
     CALL prop_step(X_atm,prop_buf_atm,t,dt,Xnew_atm,.false.) ! Obtains propagator prop_buf at X
     prop_buf_atm(2*natm+1:ndim, 2*natm+1:ndim) = 0.D0
     prop_buf_atm(2*natm+1:ndim, 1:2*natm) = 0.D0
     prop_buf_atm(1:2*natm, 2*natm+1:ndim) = 0.D0
     CALL multiply_prop_atm(prop_buf_atm) ! Multiplies prop_buf with prop
     X_atm = Xnew_atm

     X_ocn(1:2*natm) = X(1:2*natm)
     t = t - dt
     CALL prop_step(X_ocn,prop_buf_ocn,t,dt,Xnew_ocn,.false.) ! Obtains propagator prop_buf at X
     prop_buf_ocn(1:2*natm, 1:2*natm) = 0.D0
     prop_buf_ocn(2*natm+1:ndim, 1:2*natm) = 0.D0
     prop_buf_ocn(1:2*natm, 2*natm+1:ndim) = 0.D0
     CALL multiply_prop_ocn(prop_buf_ocn) ! Multiplies prop_buf with prop
     X_ocn = Xnew_ocn

     X=Xnew

     IF (mod(t,rescaling_time)<dt) THEN
        ! Multiply Kalman gain
        CALL multiply_prop_atm(eye(1:2*natm,1:2*natm) - KH_atm)
        CALL multiply_prop_ocn(eye(1+2*natm:ndim,1+2*natm:ndim) - KH_ocn)

        CALL  benettin_step ! Performs QR step with prop
        CALL lyap_acc(loclyap)

        CALL  benettin_step_atm ! Performs QR step with prop
        CALL lyap_acc_atm(loclyap_atm)

        CALL  benettin_step_ocn ! Performs QR step with prop
        CALL lyap_acc_ocn(loclyap_ocn)
     END IF
     IF (mod(t,tw)<dt) THEN
        CALL acc(X)
        IF (writeout) WRITE(10,FMTX) t,X(1:ndim)
        IndexBen=IndexBen+1
        IF (writeout) WRITE(11,rec=IndexBen,iostat=WRSTAT) loclyap
        IF (writeout) WRITE(14,rec=IndexBen,iostat=WRSTAT) loclyap_atm
        IF (writeout) WRITE(15,rec=IndexBen,iostat=WRSTAT) loclyap_ocn
     END IF
     IF (mod(t/t_run*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_run*100.D0,char(13)
  END DO
  PRINT*, 'Evolution finished.'

  IF (writeout) CLOSE(10)
  IF (writeout) CLOSE(11)
  IF (writeout) CLOSE(14)
  IF (writeout) CLOSE(15)

  IF (writeout) THEN
     OPEN(10,file='mean_lyapunov_both_3D.dat')
     lyapunov=lyap_mean()
     WRITE(10,*) 'mean',lyapunov(1:ndim)
     lyapunov=lyap_var()
     WRITE(10,*) 'var',lyapunov(1:ndim)
     CLOSE(10)
     OPEN(12,file='mean_lyapunov_atm_3D.dat')
     lyapunov_atm=lyap_mean_atm()
     WRITE(12,*) 'mean',lyapunov_atm(1:ndim)
     lyapunov_atm=lyap_var_atm()
     WRITE(12,*) 'var',lyapunov_atm(1:ndim)
     CLOSE(12)
     OPEN(13,file='mean_lyapunov_ocn_3D.dat')
     lyapunov_ocn=lyap_mean_ocn()
     WRITE(13,*) 'mean',lyapunov_ocn(1:ndim)
     lyapunov_ocn=lyap_var_ocn()
     WRITE(13,*) 'var',lyapunov_ocn(1:ndim)
  END IF

  IF (writeout) THEN
     OPEN(10,file='mean_field.dat')
     X=mean()
     WRITE(10,*) 'mean',X(1:ndim)
     X=var()
     WRITE(10,*) 'var',X(1:ndim)
  END IF

  contains

function inv(A) result(Ainv)
  real(8), dimension(:,:), intent(in) :: A
  real(8), dimension(size(A,1),size(A,2)) :: Ainv

  real(8), dimension(size(A,1)) :: work  ! work array for LAPACK
  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
  integer :: n, info

  ! External procedures defined in LAPACK
  external DGETRF
  external DGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A
  n = size(A,1)

  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call DGETRF(n, n, Ainv, n, ipiv, info)

  if (info /= 0) then
     stop 'Matrix is numerically singular!'
  end if

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call DGETRI(n, Ainv, n, ipiv, work, n, info)

  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if
end function inv

END PROGRAM maooam_lyap
