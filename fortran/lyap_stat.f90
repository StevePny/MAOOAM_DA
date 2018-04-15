
! lyap_stat.f90
!
!>  Statistics accumulators for the Lyapunov exponents
!
!> @copyright                                                               
!> 2015 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!



MODULE lyap_stat
  USE params, only: ndim
  IMPLICIT NONE

  PRIVATE
  
  INTEGER :: i=0 !< Number of stats accumulated
  
  ! Vectors holding the stats
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m, m_atm, m_ocn       !< Vector storing the inline mean
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: mprev, mprev_atm, mprev_ocn   !< Previous mean vector
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: v, v_atm, v_ocn       !< Vector storing the inline variance
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: mtmp, mtmp_atm, mtmp_ocn 


  PUBLIC :: lyap_acc,lyap_acc_atm,lyap_acc_ocn,lyap_init_stat,lyap_mean,&
      lyap_mean_atm,lyap_mean_ocn,lyap_var,lyap_var_atm,lyap_var_ocn,&
      lyap_iter,lyap_reset

  CONTAINS

    !> Initialise the accumulators
    SUBROUTINE lyap_init_stat
      INTEGER :: AllocStat
      
      ALLOCATE(m(0:ndim),m_atm(0:ndim),m_ocn(0:ndim),mprev(0:ndim),&
          mprev_atm(0:ndim),mprev_ocn(0:ndim),v(0:ndim),v_atm(0:ndim),&
          v_ocn(0:ndim),mtmp(0:ndim),mtmp_atm(0:ndim),mtmp_ocn(0:ndim),&
          STAT=AllocStat)
      IF (AllocStat /= 0) STOP '*** Not enough memory ***'
      m=0.D0
      m_atm=0.D0
      m_ocn=0.D0
      mprev=0.D0
      mprev_atm=0.D0
      mprev_ocn=0.D0
      v=0.D0
      v_atm=0.D0
      v_ocn=0.D0
      mtmp=0.D0
      mtmp_atm=0.D0
      mtmp_ocn=0.D0
      
    END SUBROUTINE lyap_init_stat

    !> Accumulate one state
    SUBROUTINE lyap_acc(x)
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: x
      i=i+1
      mprev=m+(x-m)/i
      mtmp=mprev
      mprev=m
      m=mtmp
      v=v+(x-mprev)*(x-m)
    END SUBROUTINE lyap_acc

    SUBROUTINE lyap_acc_atm(x)
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: x
      mprev_atm=m_atm+(x-m_atm)/i
      mtmp_atm=mprev_atm
      mprev_atm=m_atm
      m_atm=mtmp_atm
      v_atm=v_atm+(x-mprev_atm)*(x-m_atm)
    END SUBROUTINE lyap_acc_atm

    SUBROUTINE lyap_acc_ocn(x)
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: x
      mprev_ocn=m_ocn+(x-m_ocn)/i
      mtmp_ocn=mprev_ocn
      mprev_ocn=m_ocn
      m_ocn=mtmp_ocn
      v_ocn=v_ocn+(x-mprev_ocn)*(x-m_ocn)
    END SUBROUTINE lyap_acc_ocn

    !> Function returning the mean
    FUNCTION lyap_mean()
      REAL(KIND=8), DIMENSION(0:ndim) :: lyap_mean
      lyap_mean=m
    END FUNCTION lyap_mean

    FUNCTION lyap_mean_atm()
      REAL(KIND=8), DIMENSION(0:ndim) :: lyap_mean_atm
      lyap_mean_atm=m_atm
    END FUNCTION lyap_mean_atm

    FUNCTION lyap_mean_ocn()
      REAL(KIND=8), DIMENSION(0:ndim) :: lyap_mean_ocn
      lyap_mean_ocn=m_ocn
    END FUNCTION lyap_mean_ocn

    !> Function returning the variance
    FUNCTION lyap_var()
      REAL(KIND=8), DIMENSION(0:ndim) :: lyap_var
      lyap_var=v/(i-1)
    END FUNCTION lyap_var

    FUNCTION lyap_var_atm()
      REAL(KIND=8), DIMENSION(0:ndim) :: lyap_var_atm
      lyap_var_atm=v_atm/(i-1)
    END FUNCTION lyap_var_atm

    FUNCTION lyap_var_ocn()
      REAL(KIND=8), DIMENSION(0:ndim) :: lyap_var_ocn
      lyap_var_ocn=v_ocn/(i-1)
    END FUNCTION lyap_var_ocn

    !> Function returning the number of data accumulated
    FUNCTION lyap_iter()
      INTEGER :: lyap_iter
      lyap_iter=i
    END FUNCTION lyap_iter

    !> Routine resetting the accumulators
    SUBROUTINE lyap_reset
      m=0.D0
      m_atm = 0.D0
      m_ocn = 0.D0
      mprev=0.D0
      mprev_atm = 0.D0
      mprev_ocn = 0.D0
      v=0.D0
      v_atm = 0.D0
      v_ocn = 0.D0
      i=0
    END SUBROUTINE lyap_reset
      

  END MODULE lyap_stat
