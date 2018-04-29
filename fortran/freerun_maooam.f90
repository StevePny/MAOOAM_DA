program main
  use m_maooam, only: dt, dt_day, t_trans, t_run, &
                      get_res_maooam, init_maooam, read_x0_maooam, write_le_maooam, &
                      step
  use m_mt,     only: randn
  use m_io,     only: daio
  implicit none

  real(8),allocatable :: x(:,:)   ! trajectory (0:ndim,0:nt)
  real(8),allocatable :: stdx(:) ! ndim
  integer :: nxa, nya, nxo, nyo, ndim
  namelist /maooam/ t_trans, t_run, dt, dt_day, nt, &
                    ndim, nxa, nya, nxo, nyo

  real(8) :: frac
  namelist /makeobs/ frac

  real(8) :: t, cpu_t0, cpu_t1
  integer :: nt
  integer :: n, k, l

  CHARACTER(80) :: cin, cout


  call get_res_maooam(nxa,nya,nxo,nyo)
  call init_maooam(nxa,nya,nxo,nyo,ndim)
  nt     = NINT(t_run/dt)
  write(*,nml=maooam)
  write(9,nml=maooam)

  allocate(x(0:ndim,nt))
  allocate(stdx(ndim))

!----------------------------------------------------------------
  CALL CPU_TIME(cpu_t0)

! read 1st step of nature run trajectory
  read(daio%truth,"(10000(D24.17,1x))") (x(k,1),k=1,ndim)
  x(0,1) = 1.d0
  
  frac = 0.1d0
  read(daio%config,nml=makeobs) 
  write(*,nml=makeobs)
  do n = 1, ndim
     read(daio%clim,"(10000(D24.17,1x))") stdx(n)
     x(n,1) = x(n,1) + stdx(n)*frac
  enddo


! fwd
  t = 0.d0
  do n = 1, nt-1
     call step(x(:,n),t,dt,x(:,n+1))
     if (mod(n,nt/10)==0) print*, 100.d0*n/nt,"%"
  enddo
  write(*,*) "write out trajectory of free-run to fort.", daio%freerun
  do n = 1, nt
     write(daio%freerun,"(10000(D24.17,1x))") (x(k,n),k=1,ndim)
  enddo

  PRINT*, 'Evolution finished.'
  CALL CPU_TIME(cpu_t1)
  PRINT*, "total run time equals", cpu_t1-cpu_t0, "sec."


END PROGRAM
