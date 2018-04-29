program main
  use m_maooam, only: dt, dt_day, t_trans, t_run, &
                      get_res_maooam, init_maooam, read_x0_maooam, write_le_maooam, &
                      step
  use m_mt,     only: randn
  use m_io,     only: daio
  implicit none

  real(8),allocatable :: x(:,:)   ! trajectory (0:ndim,0:nt)
  real(8),allocatable :: xm(:)   ! ndim
  real(8),allocatable :: stdx(:) ! ndim
  integer :: nxa, nya, nxo, nyo, ndim
  namelist /maooam/ t_trans, t_run, dt, dt_day, nt, &
                    ndim, nxa, nya, nxo, nyo

  real(8) :: frac
  namelist /makeobs/ frac
  real(8),allocatable :: R(:)
  real(8),allocatable :: err(:,:)

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
  allocate(xm(ndim))
  allocate(stdx(ndim))
  call read_x0_maooam(ndim,x(:,1))

!----------------------------------------------------------------
  CALL CPU_TIME(cpu_t0)

! nature run trajectory
  t = 0.d0
  do n = 1, nt-1
     call step(x(:,n),t,dt,x(:,n+1))
     if (mod(n,nt/10)==0) print*, 100.d0*n/nt,"%"
  enddo
  write(*,*) "write out trajectory of nature run to fort.", daio%truth
  do n = 1, nt
     write(daio%truth,"(10000(D24.17,1x))") (x(k,n),k=1,ndim)
  enddo

! climatology
  write(*,*) "write out climatology to fort.", daio%clim
  do n = 1, ndim
     xm(n) = sum( x(n,:) )/nt
     stdx(n) = sqrt(max(sum((x(n,:)-xm(n))**2)/(nt-1),0.d0))
  enddo
  do n = 1, ndim
     write(daio%clim,"(10000(D24.17,1x))") stdx(n), xm(n)
     write(*,*) "n, stdx(n), xm(n)=",n,stdx(n), xm(n)
  enddo

! obs
  frac = 0.01d0 ! error = 10% * clim_std
  print*,"daio%config=", daio%config
  read(daio%config,nml=makeobs)
  print*, "frac=", frac
  allocate(R(ndim))
  allocate(err(ndim,nt))
  R(:) = stdx(:) * frac
  write(*,*) "write out obs err (R) to fort.", daio%R
  do n = 1, ndim
     write(daio%R,"(10000(D24.17,1x))") R(n)
     write(*,*) "n, R(n), stdx(n)=",n,R(n), stdx(n)
  enddo
  call randn(ndim,nt,err)
  do n = 1, ndim
     err(n,1:nt) = err(n,1:nt)*R(n)
  enddo
  x(1:ndim,:) = x(1:ndim,:)+err(:,:)
  write(*,*) "write out obs to fort.", daio%yobs
  do n = 1, nt
     write(daio%yobs,"(10000(D24.17,1x))") (x(k,n),k=1,ndim)
  enddo


  PRINT*, 'Evolution finished.'
  CALL CPU_TIME(cpu_t1)
  PRINT*, "total run time equals", cpu_t1-cpu_t0, "sec."


END PROGRAM
