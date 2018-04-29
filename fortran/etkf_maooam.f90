program main
  use params, only: natm, noc
  use m_maooam, only: dt, dt_day, t_trans, t_run, &
                      get_res_maooam, init_maooam, read_x0_maooam, write_le_maooam, &
                      step
  use m_da_maooam, only : etkf
  use m_io,        only : daio
  implicit none

  real(8),allocatable :: x(:,:)   ! trajectory (0:ndim,nt)
  real(8),allocatable :: x_comp(:,:)   ! trajectory (0:ndim,nt)
  integer :: nxa, nya, nxo, nyo, ndim
  namelist /maooam/ t_trans, t_run, dt, dt_day, nt, &
                    ndim, nxa, nya, nxo, nyo

  real(8),allocatable :: Xens(:,:) ! (0:ndim,nens)
  real(8),allocatable :: Xcomp_ens(:,:) ! (0:ndim,nens)
  real(8),allocatable :: Yens(:,:) ! (ndim,nens)
  integer :: nens
  real(8) :: infl
  real(8),allocatable :: KH(:, :)
  integer :: lout_etkf
  namelist /da_etkf/ nens, infl, lout_etkf

  real(8),allocatable :: yobs(:,:) ! (ndim,nt)
  real(8),allocatable :: R(:)
  integer :: nyobs
  logical,allocatable :: luse(:)


  real(8) :: t, cpu_t0, cpu_t1
  integer :: nt
  integer :: n, k, l


  CHARACTER(80) :: cin, cout

  CHARACTER(3) :: comp

  comp = 'atm'

  call get_res_maooam(nxa,nya,nxo,nyo)
  call init_maooam(nxa,nya,nxo,nyo,ndim)
  nt     = NINT(t_run/dt)
  write(*,nml=maooam)
  write(9,nml=maooam)

  allocate(x(0:ndim,nt))
  allocate(x_comp(0:ndim,nt))
  call read_x0_maooam(ndim,x(:,1))
  call read_x0_maooam(ndim,x_comp(:,1))

!----------------------------------------------------------------


  nens = 10
  infl = 1.00d0
  lout_etkf = 12
  read(daio%config,nml=da_etkf) 
  write(*,nml=da_etkf)
  pause "continue"

  CALL CPU_TIME(cpu_t0)
  nyobs = ndim
  allocate(Xens(0:ndim,nens))
  allocate(Xcomp_ens(0:ndim,nens))
  allocate(R(nyobs))
  allocate(luse(nyobs)); luse = .true.
  allocate(yobs(ndim,nt))
  allocate(Yens(nyobs,nens))
  allocate(KH(ndim,ndim))

! load obs
  do n = 1, nt
     read(daio%yobs,"(10000(D24.17,1x))") (yobs(k,n),k=1,nyobs)
  enddo
  print*, "y(:,end)=", yobs(:,nt)
  do n = 1, nyobs
     read(daio%R,"(10000(D24.17,1x))") R(n)
     print*, "n, R(n)=", n, R(n)
  enddo

  pause "set Ens"

! set initial ensembles
  read(daio%freerun,"(10000(D24.17,1x))") (x(k,1),k=1,ndim)  ! freerun
  x(0,1) = 1.d0
  do n = 1, nt-1
     call step(x(:,n),t,dt,x(:,n+1))
     if (comp == 'atm') then
         t = t - dt
         x_comp(2*natm+1:ndim, n) = x(2*natm+1:ndim, n)
         call step(x_comp(:,n), t, dt, x_comp(:,n+1))
     else if (comp == 'ocn') then
         t = t - dt
         x_comp(1:2*natm, n) = x(1:2*natm, n)
         call step(x_comp(:,n), t, dt, x_comp(:,n+1))
     end if
  enddo
  do n = 1, nens
     k = (nt/nens)*n
     Xens(1:ndim,n) = x(1:ndim,k)
     Xens(0,n)      = 1.d0
     print*, "ensemble: mem, step=", n, k, nt
  enddo

  call init_stat
  call lyap_init_stat

! run the DA cycle
  x(0,2:nt) = 0.0d0
  do n = 1, nt-1
     do k = 1, nens
        ! ensemble forecast
        call step(Xens(:,k), t, dt, Xens(:,k))

        if (comp == 'atm') then
            t = t - dt
            Xcomp_ens(2*natm+1:ndim, k) = Xens(2*natm+1:ndim, k)
            call step(Xcomp_ens(:,k), t, dt, Xcomp_ens(:,k))
        else if (comp == 'ocn') then
            t = t - dt
            Xcomp_ens(2*natm+1:ndim, k) = Xens(2*natm+1:ndim, k)
            call step(Xcomp_ens(:,k), t, dt, Xcomp_ens(:,k))
        end if

        ! H(x)
        Yens(:,k) = Xens(1:ndim,k)
     enddo

     ! luse
     ! (1:natm_maooam): atm streamfunction
     ! (natm_maooam+1:2*natm_maooam): atm temp
     ! (2*natm_maooam:2*natm_maooam+nocn_maooam): ocn streamfunction
     ! (2*natm_maooam+1:end): ocn temp
     !
     call etkf( ndim, &
                ndim, &
                nens, &
                Xens(1:ndim,:), &
                luse, &
                yobs(:,n+1), &
                R, &
                infl, &
                Xens(1:ndim,:), &
                x(1:ndim,n+1), &
                KH )

  enddo

  do n = 1, nt
     write(lout_etkf,"(10000(D24.17,1x))") (x(k,n),k=1,ndim)
  enddo




  PRINT*, 'Evolution finished.'
  CALL CPU_TIME(cpu_t1)
  PRINT*, "total run time equals", cpu_t1-cpu_t0, "sec."


END PROGRAM
