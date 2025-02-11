module machine

  use iso_fortran_env, only: int32, real32  ! Use single precision from ISO_FORTRAN_ENV

  implicit none

  integer, parameter :: kind_io4 = int32, kind_io8 = int32, kind_ior = int32
  integer, parameter :: kind_evod = real32, kind_dbl_prec = real32
  integer, parameter :: kind_qdt_prec = real32
  integer, parameter :: kind_rad = real32, kind_phys = real32, kind_taum = real32
  integer, parameter :: kind_grid = real32
  integer, parameter :: kind_REAL = real32  ! Used in cmp_comm
  integer, parameter :: kind_LOGICAL = int32
  integer, parameter :: kind_INTEGER = int32

  integer, parameter :: kind_dyn = real32

  real(kind=kind_evod), parameter :: mprec = 1.0e-6_real32  ! Machine precision for single precision
  real(kind=kind_evod), parameter :: grib_undef = 9.99e20_real32  ! GRIB undefined value

end module machine
