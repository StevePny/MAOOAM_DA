
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
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X, X_atm, X_ocn          !< State variable in the model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: err_atm
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: err_ocn
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Xnew, Xnew_atm, Xnew_ocn       !< Updated state variable
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: prop_buf, prop_buf_atm, prop_buf_ocn !< Buffer for the propagator
  REAL(KIND=8) :: t=0.D0                                !< Time variable
  REAL(KIND=8) :: t_up, gasdev
  INTEGER :: IndexBen,WRSTAT
  CHARACTER(LEN=19) :: FMTX
  CHARACTER(LEN=3) :: sim_id
  INTEGER :: idum1
  INTEGER :: j

  idum1=-1254

  PRINT*, 'Model MAOOAM v1.0'
  PRINT*, '      - with computation of the Lyapunov spectrum'
  PRINT*, 'Loading information...'

  CALL get_command_argument(1, sim_id)

  CALL init_aotensor(sim_id)    ! Compute the tensors
  CALL init_tltensor
  CALL load_IC          ! Load the initial condition

  CALL init_integrator        ! Initialize the integrator
  CALL init_tl_ad_integrator  ! Initialize tangent linear integrator

  t_trans = 1.D3
  write(FMTX,'(A10,i3,A6)') '(F10.2,4x,',ndim,'E15.5)'
  t_up=dt/t_trans*100.D0

  IF (writeout) THEN
     OPEN(10,file='evol_field_' // sim_id // '.dat')
     OPEN(14,file='evol_field_atm_' // sim_id // '.dat')
     OPEN(15,file='evol_field_ocn_' // sim_id // '.dat')
  END IF

  ALLOCATE(X(0:ndim),Xnew(0:ndim),X_atm(0:ndim),Xnew_atm(0:ndim),&
      X_ocn(0:ndim),Xnew_ocn(0:ndim),prop_buf(ndim,ndim),&
      prop_buf_atm(ndim,ndim),prop_buf_ocn(ndim,ndim),err_atm(1:2*natm),&
      err_ocn(1:2*noc))
  X = IC


  DO j=1,2*natm
    err_atm(j) = 1D-1*gasdev(idum1)
  END DO
  DO j=1,2*noc
    err_ocn(j) = 1D-1*gasdev(idum1)
  END DO

  X_atm = IC
  X_ocn = IC

  X_ocn(2*natm+1:ndim) = X_ocn(2*natm+1:ndim) + err_ocn
  X_atm(1:2*natm) = X_atm(1:2*natm) + err_atm

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

    IF (writeout) WRITE(10,FMTX) t,X(1:ndim)
    IF (writeout) WRITE(14,FMTX) t,X_atm(1:ndim)
    IF (writeout) WRITE(15,FMTX) t,X_ocn(1:ndim)
  END DO
  PRINT*, 'Evolution finished.'

  IF (writeout) CLOSE(10)
  IF (writeout) CLOSE(11)
  IF (writeout) CLOSE(14)
  IF (writeout) CLOSE(15)

END PROGRAM maooam_lyap

FUNCTION gasdev(idum)
  INTEGER :: idum
  REAL(KIND=8) ::  gasdev,ran2
  !    USES ran2
  INTEGER :: iset
  REAL(KIND=8) :: fac,gset,rsq,v1,v2
  SAVE iset,gset
  DATA iset/0/
  if (idum.lt.0) iset=0
  if (iset.eq.0) then
1    v1=2.D0*ran2(idum)-1.
     v2=2.D0*ran2(idum)-1.
     rsq=v1**2+v2**2
     if (rsq.ge.1.D0.or.rsq.eq.0.D0) goto 1
     fac=sqrt(-2.*log(rsq)/rsq)
     gset=v1*fac
     gasdev=v2*fac
     iset=1
  else
     gasdev=gset
     iset=0
  endif
  return
END FUNCTION gasdev

FUNCTION ran2(idum)
  INTEGER :: idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
  REAL(KIND=8) :: ran2,AM,EPS,RNMX
  PARAMETER (IM1=2147483563,IM2=2147483399,AM=1.D0/IM1,IMM1=IM1-1&
       &,IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2&
       &=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2D-7,RNMX=1.D0-EPS)
  INTEGER :: idum2,j,k,iv(NTAB),iy
  SAVE iv,iy,idum2
  DATA idum2/123456789/, iv/NTAB*0/, iy/0/
  if (idum.le.0) then
     idum=max(-idum,1)
     idum2=idum
     do j=NTAB+8,1,-1
        k=idum/IQ1
        idum=IA1*(idum-k*IQ1)-k*IR1
        if (idum.lt.0) idum=idum+IM1
        if (j.le.NTAB) iv(j)=idum
     enddo
     iy=iv(1)
  endif
  k=idum/IQ1
  idum=IA1*(idum-k*IQ1)-k*IR1
  if (idum.lt.0) idum=idum+IM1
  k=idum2/IQ2
  idum2=IA2*(idum2-k*IQ2)-k*IR2
  if (idum2.lt.0) idum2=idum2+IM2
  j=1+iy/NDIV
  iy=iv(j)-idum2
  iv(j)=idum
  if (iy.lt.1) iy=iy+IMM1
  ran2=min(AM*iy,RNMX)
  return
END FUNCTION ran2
