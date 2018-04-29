module m_io
  implicit none

  private
  public :: daio

  type t_daio
    integer :: truth = 200  ! Natural run (Index number) 
    integer :: clim  = 201  ! Natural variability
    integer :: R     = 202  ! R matrix
    integer :: freerun = 203  
    integer :: yobs   = 204  ! Observation
    integer :: B      = 205  ! B for total
    integer :: Batm   = 206  ! B for ATM only
    integer :: Bocn   = 207  ! B for Ocean only

    integer :: config = 300  ! All namelist of DA configuration
  endtype
  type(t_daio) :: daio


endmodule
