
! params.f90                                                                
!
!>  The model parameters module. 
!
!> @copyright                                                               
!> 2015 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------
!                                                                           
!>  @remark                                                                 
!>  Once the init_params() subroutine is called, the parameters are loaded           
!>  globally in the main program and its subroutines and function           
!                                                                           
!---------------------------------------------------------------------------

MODULE params

  IMPLICIT NONE

  PUBLIC

  REAL(KIND=8) :: n         !< \f$n = 2 L_y / L_x\f$ - Aspect ratio
  REAL(KIND=8) :: phi0      !< Latitude in radian
  REAL(KIND=8) :: rra       !< Earth radius
  REAL(KIND=8) :: sig0      !< \f$\sigma_0\f$ - Non-dimensional static stability of the atmosphere.
  REAL(KIND=8) :: k         !< Bottom atmospheric friction coefficient.
  REAL(KIND=8) :: kp        !< \f$k'\f$ - Internal atmospheric friction coefficient.
  REAL(KIND=8) :: r         !< Frictional coefficient at the bottom of the ocean.
  REAL(KIND=8) :: d         !< Merchanical coupling parameter between the ocean and the atmosphere.
  REAL(KIND=8) :: f0        !< \f$f_0\f$ - Coriolis parameter
  REAL(KIND=8) :: gp        !< \f$g'\f$Reduced gravity
  REAL(KIND=8) :: H         !< Depth of the active water layer of the ocean.
  REAL(KIND=8) :: phi0_npi  !< Latitude exprimed in fraction of pi.

  REAL(KIND=8) :: lambda    !< \f$\lambda\f$ - Sensible + turbulent heat exchange between the ocean and the atmosphere.
  REAL(KIND=8) :: Co        !< \f$C_a\f$ - Constant short-wave radiation of the ocean.
  REAL(KIND=8) :: Go        !< \f$\gamma_o\f$ - Specific heat capacity of the ocean.
  REAL(KIND=8) :: Ca        !< \f$C_a\f$ - Constant short-wave radiation of the atmosphere.
  REAL(KIND=8) :: To0       !< \f$T_o^0\f$ -  Stationary solution for the 0-th order ocean temperature.
  REAL(KIND=8) :: Ta0       !< \f$T_a^0\f$ -  Stationary solution for the 0-th order atmospheric temperature.
  REAL(KIND=8) :: epsa      !< \f$\epsilon_a\f$ - Emissivity coefficient for the grey-body atmosphere.
  REAL(KIND=8) :: Ga        !< \f$\gamma_a\f$ - Specific heat capacity of the atmosphere.
  REAL(KIND=8) :: RR        !< \f$R\f$ - Gas constant of dry air

  REAL(KIND=8) :: scale     !< \f$L_y = L \, \pi\f$ - The characteristic space scale.
  REAL(KIND=8) :: pi        !< \f$\pi\f$
  REAL(KIND=8) :: LR        !< \f$L_R\f$ - Rossby deformation radius
  REAL(KIND=8) :: G         !< \f$\gamma\f$
  REAL(KIND=8) :: rp        !< \f$r'\f$ - Frictional coefficient at the bottom of the ocean.
  REAL(KIND=8) :: dp        !< \f$d'\f$ - Non-dimensional mechanical coupling parameter between the ocean and the atmosphere.
  REAL(KIND=8) :: kd        !< \f$k_d\f$ - Non-dimensional bottom atmospheric friction coefficient.
  REAL(KIND=8) :: kdp       !< \f$k'_d\f$ - Non-dimensional internal atmospheric friction coefficient.

  REAL(KIND=8) :: Cpo       !< \f$C'_a\f$ - Non-dimensional constant short-wave radiation of the ocean.
  REAL(KIND=8) :: Lpo       !< \f$\lambda'_o\f$ - Non-dimensional sensible + turbulent heat exchange from ocean to atmosphere.
  REAL(KIND=8) :: Cpa       !< \f$C'_a\f$ - Non-dimensional constant short-wave radiation of the atmosphere. @remark Cpa acts on psi1-psi3, not on theta.
  REAL(KIND=8) :: Lpa       !< \f$\lambda'_a\f$ - Non-dimensional sensible + turbulent heat exchange from atmosphere to ocean.
  REAL(KIND=8) :: sBpo      !< \f$\sigma'_{B,o}\f$ - Long wave radiation lost by ocean to atmosphere & space.
  REAL(KIND=8) :: sBpa      !< \f$\sigma'_{B,a}\f$ - Long wave radiation from atmosphere absorbed by ocean.
  REAL(KIND=8) :: LSBpo     !< \f$S'_{B,o}\f$ - Long wave radiation from ocean absorbed by atmosphere.
  REAL(KIND=8) :: LSBpa     !< \f$S'_{B,a}\f$ - Long wave radiation lost by atmosphere to space & ocean.
  REAL(KIND=8) :: L         !< \f$L\f$ - Domain length scale
  REAL(KIND=8) :: sc        !< Ratio of surface to atmosphere temperature.
  REAL(KIND=8) :: sB        !< Stefan–Boltzmann constant
  REAL(KIND=8) :: betp      !< \f$\beta'\f$ - Non-dimensional beta parameter

  REAL(KIND=8) :: nua=0.D0  !< Dissipation in the atmosphere
  REAL(KIND=8) :: nuo=0.D0  !< Dissipation in the ocean

  REAL(KIND=8) :: nuap      !< Non-dimensional dissipation in the atmosphere
  REAL(KIND=8) :: nuop      !< Non-dimensional dissipation in the ocean

  REAL(KIND=8) :: t_trans   !< Transient time period
  REAL(KIND=8) :: t_run     !< Effective intergration time (length of the generated trajectory)
  REAL(KIND=8) :: dt        !< Integration time step
  REAL(KIND=8) :: tw        !< Write all variables every tw time units
  LOGICAL :: writeout       !< Write to file boolean

  REAL(KIND=8) :: rescaling_time !< Rescaling time for the Lyapunov computation
  
  REAL(KIND=8) :: coupling_thermo
  REAL(KIND=8) :: coupling_motion
  INTEGER :: uncoupled

  INTEGER :: nboc   !< Number of atmospheric blocks
  INTEGER :: nbatm  !< Number of oceanic blocks
  INTEGER :: natm=0 !< Number of atmospheric basis functions
  INTEGER :: noc=0  !< Number of oceanic basis functions
  INTEGER :: ndim   !< Number of variables (dimension of the model)
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: oms   !< Ocean mode selection array
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: ams   !< Atmospheric mode selection array

  PRIVATE :: init_nml


CONTAINS

  !> copied from python version.
  !> generate AMS or OMS
  SUBROUTINE get_modes(nxmax,nymax,res)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: nxmax, nymax
    INTEGER,ALLOCATABLE :: res(:,:)

    INTEGER :: i, ix, iy

    ALLOCATE(res(nxmax*nymax,2))
    i = 0
    do ix = 1, nxmax
       do iy = 1, nymax
          i = i + 1
          res(i,1) = ix
          res(i,2) = iy
       enddo
    enddo

  ENDSUBROUTINE


  SUBROUTINE fstrout(nxa,nya,nxo,nyo,cout)  
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: nxa, nya, nxo, nyo
    CHARACTER(*),INTENT(INOUT) :: cout
    ! .ATMS002x002_OCN002x004
    write(cout,"(A,I3.3,A,I3.3,A,I3.3,A,I3.3)") "ATMS",nxa,"x",nya,"_OCN",nxo,"x",nyo

  ENDSUBROUTINE


  !< explicitly calculate resoluation-related vars
  !< CDA
  SUBROUTINE init_res(nxa,nya,nxo,nyo,lout)
    IMPLICIT NONE
    
    INTEGER,INTENT(IN) :: nxa, nya, nxo, nyo
    LOGICAL,OPTIONAL,INTENT(IN) :: lout

    CHARACTER(80) :: cout
    INTEGER :: i

    NAMELIST /modeselection/ oms,ams
    NAMELIST /numblocs/ nboc,nbatm
    NAMELIST /res/ nxa, nya, nxo, nyo

    nboc  = nxo * nyo
    nbatm = nxa * nya
    !ALLOCATE(oms(nboc,2),ams(nbatm,2), STAT=AllocStat)
    !IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
    call get_modes(nxa,nya,ams)
    call get_modes(nxo,nyo,oms)
    IF (ANY(shape(ams)/=(/nbatm,2/)) ) STOP "*** Error of AMS shape"
    IF (ANY(shape(oms)/=(/nboc,2/)) ) STOP "*** Error of OMS shape"
    
    !write(*,*) "nboc=", nboc
    !write(*,*) "nbatm=", nbatm
    !do i = 1, nbatm
    !   write(*,*) "AMS(",i, ",1:2)=",ams(i,1:2)
    !enddo
    !do i = 1, nboc
    !   write(*,*) "OMS(",i, ",1:2)=",oms(i,1:2)
    !enddo

    IF (PRESENT(lout).and.lout) THEN
       CALL fstrout(nxa,nya,nxo,nyo,cout)
       cout="mode."//TRIM(cout)//".nml"
       WRITE(*,*) "/modeselection/ & /numblocs/ written into file:", TRIM(cout)
       OPEN(18,FILE=TRIM(cout),ACTION="write")
       WRITE(18,NML=res)
       WRITE(18,NML=numblocs)
       WRITE(18,NML=modeselection)
       CLOSE(18)
    ENDIF

  ENDSUBROUTINE init_res


  !> Read the basic parameters and mode selection from the namelist.
  SUBROUTINE init_nml(sim_id)
    CHARACTER(len=3), OPTIONAL :: sim_id
    INTEGER :: AllocStat

    NAMELIST /aoscale/  scale,f0,n,rra,phi0_npi
    NAMELIST /oparams/  gp,r,H,d,nuo
    NAMELIST /aparams/  k,kp,sig0,nua
    NAMELIST /toparams/ Go,Co,To0
    NAMELIST /taparams/ Ga,Ca,epsa,Ta0
    NAMELIST /otparams/ sc,lambda,RR,sB
    NAMELIST /cparams/  coupling_thermo, coupling_motion, uncoupled

    NAMELIST /modeselection/ oms,ams
    NAMELIST /numblocs/ nboc,nbatm

    NAMELIST /int_params/ t_trans,t_run,dt,tw,writeout
    NAMELIST /lyap_params/ rescaling_time

    IF (PRESENT(sim_id)) THEN
        OPEN(8, file="params/params_" // sim_id // ".nml", status='OLD', recl=80, delim='APOSTROPHE')
    ELSE
        OPEN(8, file="params.nml", status='OLD', recl=80, delim='APOSTROPHE')
    END IF

    READ(8,nml=aoscale)
    READ(8,nml=oparams)
    READ(8,nml=aparams)
    READ(8,nml=toparams)
    READ(8,nml=taparams)
    READ(8,nml=otparams)
    READ(8,nml=cparams)

    CLOSE(8)

    OPEN(8, file="modeselection.nml", status='OLD', recl=80, delim='APOSTROPHE')
    READ(8,nml=numblocs)

    ALLOCATE(oms(nboc,2),ams(nbatm,2), STAT=AllocStat)
    !    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    READ(8,nml=modeselection)
    CLOSE(8)

    OPEN(8, file="int_params.nml", status='OLD', recl=80, delim='APOSTROPHE')
    READ(8,nml=int_params)
    READ(8,nml=lyap_params)

  END SUBROUTINE init_nml

  !> Parameters initialisation routine
  SUBROUTINE init_params(sim_id)
    CHARACTER(len=3), OPTIONAL :: sim_id
    INTEGER, DIMENSION(2) :: s
    INTEGER :: i

    IF (PRESENT(sim_id)) THEN
        CALL init_nml(sim_id)
    ELSE
        CALL init_nml
    END IF

    !---------------------------------------------------------!
    !                                                         !
    ! Computation of the dimension of the atmospheric         !
    ! and oceanic components                                  !
    !                                                         !
    !---------------------------------------------------------!

    natm=0
    DO i=1,nbatm
       IF (ams(i,1)==1) THEN
          natm=natm+3
       ELSE
          natm=natm+2
       ENDIF
    ENDDO
    s=shape(oms)
    noc=s(1)

    IF (uncoupled == 0) THEN
       ndim=2*natm+2*noc
    ELSE IF (uncoupled == 1) THEN
       ndim=4*natm+2*noc
    ELSE IF (uncoupled == 2) THEN
       ndim=2*natm+4*noc
    END IF

    !---------------------------------------------------------!
    !                                                         !
    ! Some general parameters (Domain, beta, gamma, coupling) !
    !                                                         !
    !---------------------------------------------------------!

    d=coupling_motion*d
    pi=dacos(-1.D0)
    L=scale/pi
    phi0=phi0_npi*pi
    LR=sqrt(gp*H)/f0
    G=-L**2/LR**2
    betp=L/rra*cos(phi0)/sin(phi0)
    rp=r/f0
    dp=d/f0
    kd=coupling_motion*k*2
    kdp=kp
    lambda=coupling_thermo*lambda

    !-----------------------------------------------------!
    !                                                     !
    ! DERIVED QUANTITIES                                  !
    !                                                     !
    !-----------------------------------------------------!

    Cpo=Co/(Go*f0) * RR/(f0**2*L**2)
    Lpo=lambda/(Go*f0)
    Cpa=Ca/(Ga*f0) * RR/(f0**2*L**2)/2 ! Cpa acts on psi1-psi3, not on theta
    Lpa=lambda/(Ga*f0)
    sBpo=coupling_thermo*4*sB*To0**3/(Go*f0) ! long wave radiation lost by ocean to atmosphere space
    sBpa=coupling_thermo*8*epsa*sB*Ta0**3/(Go*f0) ! long wave radiation from atmosphere absorbed by ocean
    LSBpo=coupling_thermo*2*epsa*sB*To0**3/(Ga*f0) ! long wave radiation from ocean absorbed by atmosphere
    LSBpa=coupling_thermo*8*epsa*sB*Ta0**3/(Ga*f0) ! long wave radiation lost by atmosphere to space & ocea
    nuap=nua/(f0*L**2)
    nuop=nuo/(f0*L**2)

  END SUBROUTINE init_params
END MODULE params
