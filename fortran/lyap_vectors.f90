
! lyap_vectors.f90
!
!> Module for computation of Lyapunov exponents and vectors
!
!> @copyright                                                               
!> 2016 Sebastian Schubert.
!> See LICENSE.txt for license information.
!
!---------------------------------------------------------------------------!
!                                                                           
!>  @remark                                                                 
!>  This module contains the necessary tools to perform the Benettin
!>  steps to compute the lyapunov exponents. (Ginelli for CLV will be added later)
!>
!>  References :
!>  Benettin, G., Galgani, L., Giorgilli, A., & Strelcyn, J. M. (1980). Lyapunov
!>  characteristic exponents for smooth dynamical systems; a method for computing
!>  all of them. Part 2: Numerical application. \a Meccanica \a, 15, 21-30.
!                                                                           
!---------------------------------------------------------------------------


MODULE lyap_vectors

  !-----------------------------------------------------!
  !                                                     !
  ! Preamble and variables declaration                  !
  !                                                     !
  !-----------------------------------------------------!

  USE params, only: ndim,natm,noc,dt,rescaling_time
  USE util, only: init_one
  IMPLICIT NONE

  PRIVATE
  
  PUBLIC :: benettin_step,benettin_step_atm,benettin_step_ocn,loclyap,&
      loclyap_atm,loclyap_ocn,lyapunov,lyapunov_atm,lyapunov_ocn,ensemble,&
      ensemble_atm,ensemble_ocn,init_lyap,init_lyap_atm,init_lyap_ocn,&
      multiply_prop,multiply_prop_atm,multiply_prop_ocn,get_lyap_state,&
      get_lyap_state_atm,get_lyap_state_ocn
 
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: loclyap, loclyap_atm, loclyap_ocn    !< Buffer containing the local Lyapunov exponent
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: lyapunov, lyapunov_atm, lyapunov_ocn   !< Buffer containing the averaged Lyapunov exponent
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: ensemble, ensemble_atm, ensemble_ocn !< Buffer containing the QR decomposition of the ensemble
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: prop, prop_atm, prop_ocn     !< Buffer holding the propagator matrix
  
  INTEGER :: lwork
  REAL(kind=8), DIMENSION(:), ALLOCATABLE :: work       !< Temporary buffer for QR decomposition
  REAL(kind=8), DIMENSION(:), ALLOCATABLE :: work2      !< Temporary buffer for QR decomposition
  REAL(kind=8), DIMENSION(:), ALLOCATABLE :: tau, tau_atm, tau_ocn        !< Temporary buffer for QR decomposition
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: prop_buf, prop_buf_atm, prop_buf_ocn !< Buffer holding the local propagator matrix

  !-----------------------------------------------------!
  !                                                     !
  ! End of preamble                                     !
  !                                                     !
  !-----------------------------------------------------!


CONTAINS
 !-----------------------------------------------------!
 !                                                     !
 ! Function declarations                               !
 !                                                     !
 !-----------------------------------------------------!


  !> Initialize Lyapunov computation (possibly also vectors in later version)
  !> and initializes also a random orthogonal matrix for the matrix ensemble. 
  SUBROUTINE init_lyap
    INTEGER :: AllocStat,ilaenv,info
    lwork=ilaenv(1,"dgeqrf"," ",ndim,ndim,ndim,-1)
    lwork=ndim*lwork
    ALLOCATE(prop_buf(ndim,ndim),lyapunov(ndim),loclyap(ndim),ensemble(ndim,ndim),tau(ndim),prop(ndim,ndim), &
    & work2(ndim),work(lwork),STAT=AllocStat) 
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    lyapunov=0.0d0
    loclyap=0.0d0
    CALL init_one(prop)
    CALL random_number(ensemble)
    ensemble=2*(ensemble-0.5)
    CALL DGEQRF(ndim,ndim,ensemble,ndim,tau,work,lwork, info) ! qr decomposition
  END SUBROUTINE init_lyap

  SUBROUTINE init_lyap_atm
    INTEGER :: AllocStat,ilaenv,info
    lwork=ilaenv(1,"dgeqrf"," ",ndim,ndim,ndim,-1)
    lwork=ndim*lwork
    ALLOCATE(lyapunov_atm(ndim),loclyap_atm(ndim),ensemble_atm(ndim,ndim),&
        prop_atm(ndim,ndim),prop_buf_atm(ndim,ndim),tau_atm(ndim),STAT=AllocStat) 
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    lyapunov_atm=0.0d0
    loclyap_atm=0.0d0
    CALL init_one(prop_atm)
    CALL random_number(ensemble_atm)
    ensemble_atm=2*(ensemble_atm-0.5)
    CALL DGEQRF(ndim,ndim,ensemble_atm,ndim,tau_atm,work,lwork, info) ! qr decomposition
  END SUBROUTINE init_lyap_atm

  SUBROUTINE init_lyap_ocn
    INTEGER :: AllocStat,ilaenv,info
    lwork=ilaenv(1,"dgeqrf"," ",ndim,ndim,ndim,-1)
    lwork=ndim*lwork
    ALLOCATE(lyapunov_ocn(ndim),loclyap_ocn(ndim),ensemble_ocn(ndim,ndim),&
    prop_ocn(ndim,ndim),prop_buf_ocn(ndim,ndim),tau_ocn(ndim),STAT=AllocStat) 
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    lyapunov_ocn=0.0d0
    loclyap_ocn=0.0d0
    CALL init_one(prop_ocn)
    CALL random_number(ensemble_ocn)
    ensemble_ocn=2*(ensemble_ocn-0.5)
    CALL DGEQRF(ndim,ndim,ensemble_ocn,ndim,tau_ocn,work,lwork, info) ! qr decomposition
  END SUBROUTINE init_lyap_ocn

  !> Multiplies prop_mul from the left with the prop matrix defined in this
  !> module and saves the result to prop_mul
  !> @param prop_mul local propagator to multiply with the global one
  SUBROUTINE multiply_prop(prop_mul)
    REAL(KIND=8), DIMENSION(ndim,ndim),INTENT(IN) :: prop_mul
    prop_buf=prop    
    CALL DGEMM ('n', 'n', ndim, ndim, ndim, 1.0d0, prop_mul, ndim,prop_buf, ndim,0.0d0, prop, ndim)
  END SUBROUTINE multiply_prop

  SUBROUTINE multiply_prop_atm(prop_mul)
    REAL(KIND=8), DIMENSION(ndim,ndim),INTENT(IN) :: prop_mul
    prop_buf_atm=prop_atm    
    CALL DGEMM ('n', 'n', ndim, ndim, ndim, 1.0d0, prop_mul, ndim,prop_buf_atm, ndim,0.0d0, prop_atm, ndim)
  END SUBROUTINE multiply_prop_atm

  SUBROUTINE multiply_prop_ocn(prop_mul)
    REAL(KIND=8), DIMENSION(ndim,ndim),INTENT(IN) :: prop_mul
    prop_buf_ocn=prop_ocn    
    CALL DGEMM ('n', 'n', ndim, ndim, ndim, 1.0d0, prop_mul, ndim,prop_buf_ocn, ndim,0.0d0, prop_ocn, ndim)
  END SUBROUTINE multiply_prop_ocn

  !> Performs the benettin step in integration. Multiplies the aggregated
  !> propagators in prop with ensemble and performs QR decomposition (Gram-Schmidt
  !> orthogonalization gives Q and upper triangular matrix R). Computes also the
  !> Lyapunov exponents via the diagonal of R. WATCH OUT: prop is changed during
  !> the subroutine and restored to a unit matrix
  SUBROUTINE benettin_step
    INTEGER :: info,k

    ! Multiply the Propagator prop from the right side with the non transposed q matrix
    ! from the qr decomposition which is stored in ensemble.
    CALL DORM2R("r","n",ndim,ndim,ndim,ensemble,ndim,tau,prop,ndim,work2,info)
    ! prop contains prop*ensemble but QR decomposed(tau is needed for that as
    ! well !) => copy to ensemble 
    ensemble=prop

    ! From here on ensemble contains the new information prop*ensemble
    CALL DGEQRF(ndim,ndim,ensemble,ndim,tau,work,lwork, info) ! qr decomposition
    
    DO k=1,ndim
      loclyap(k)=log(abs(ensemble(k,k)))/rescaling_time
    END DO

    !
    ! Add here save for 
    !

    ! Initialise prop again with unit matrix
    CALL init_one(prop) 
    
   END SUBROUTINE benettin_step

  SUBROUTINE benettin_step_atm
    INTEGER :: info,k

    ! Multiply the Propagator prop from the right side with the non transposed q matrix
    ! from the qr decomposition which is stored in ensemble.
    CALL DORM2R("r","n",ndim,ndim,ndim,ensemble_atm,ndim,tau_atm,prop_atm,ndim,work2,info)
    ! prop contains prop*ensemble but QR decomposed(tau is needed for that as
    ! well !) => copy to ensemble 
     prop_atm(2*natm:ndim, 2*natm:ndim) = 0.D0
     prop_atm(2*natm:ndim, 1:2*natm) = 0.D0
     prop_atm(1:2*natm, 2*natm:ndim) = 0.D0
    ensemble_atm=prop_atm

    ! From here on ensemble contains the new information prop*ensemble
    CALL DGEQRF(ndim,ndim,ensemble_atm,ndim,tau_atm,work,lwork, info) ! qr decomposition
    
     ensemble_atm(2*natm:ndim, 2*natm:ndim) = 0.D0
     ensemble_atm(2*natm:ndim, 1:2*natm) = 0.D0
     ensemble_atm(1:2*natm, 2*natm:ndim) = 0.D0

    DO k=1,ndim
      loclyap_atm(k)=log(abs(ensemble_atm(k,k)))/rescaling_time
    END DO

    !
    ! Add here save for 
    !

    ! Initialise prop again with unit matrix
    CALL init_one(prop_atm) 
    
   END SUBROUTINE benettin_step_atm

  SUBROUTINE benettin_step_ocn
    INTEGER :: info,k

    ! Multiply the Propagator prop from the right side with the non transposed q matrix
    ! from the qr decomposition which is stored in ensemble.
    CALL DORM2R("r","n",ndim,ndim,ndim,ensemble_ocn,ndim,tau_ocn,prop_ocn,ndim,work2,info)
    ! prop contains prop*ensemble but QR decomposed(tau is needed for that as
    ! well !) => copy to ensemble 
     prop_ocn(1:2*natm, 1:2*natm) = 0.D0
     prop_ocn(2*natm:ndim, 1:2*natm) = 0.D0
     prop_ocn(1:2*natm, 2*natm:ndim) = 0.D0
    ensemble_ocn=prop_ocn

    ! From here on ensemble contains the new information prop*ensemble
    CALL DGEQRF(ndim,ndim,ensemble_ocn,ndim,tau_ocn,work,lwork, info) ! qr decomposition
     ensemble_ocn(1:2*natm, 1:2*natm) = 0.D0
     ensemble_ocn(2*natm:ndim, 1:2*natm) = 0.D0
     ensemble_ocn(1:2*natm, 2*natm:ndim) = 0.D0
   
    DO k=1,ndim
      loclyap_ocn(k)=log(abs(ensemble_ocn(k,k)))/rescaling_time
    END DO

    !
    ! Add here save for 
    !

    ! Initialise prop again with unit matrix
    CALL init_one(prop_ocn) 
    
   END SUBROUTINE benettin_step_ocn


   !> Routine that returns the current global propagator and ensemble of
   !> lyapunov vectors
   SUBROUTINE get_lyap_state(prop_ret,ensemble_ret)
     REAL(KIND=8), DIMENSION(ndim,ndim),INTENT(OUT) :: prop_ret,ensemble_ret
     prop_ret=prop
     ensemble_ret=ensemble
   END SUBROUTINE get_lyap_state

   SUBROUTINE get_lyap_state_atm(prop_ret,ensemble_ret)
     REAL(KIND=8), DIMENSION(ndim,ndim),INTENT(OUT) :: prop_ret,ensemble_ret
     prop_ret=prop_atm
     ensemble_ret=ensemble_atm
   END SUBROUTINE get_lyap_state_atm

   SUBROUTINE get_lyap_state_ocn(prop_ret,ensemble_ret)
     REAL(KIND=8), DIMENSION(ndim,ndim),INTENT(OUT) :: prop_ret,ensemble_ret
     prop_ret=prop_ocn
     ensemble_ret=ensemble_ocn
   END SUBROUTINE get_lyap_state_ocn

END MODULE lyap_vectors
     
