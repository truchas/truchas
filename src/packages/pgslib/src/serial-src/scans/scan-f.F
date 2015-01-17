!!
!! scan-f.F
!!
!! Serial-only dummy routines that correspond to the MPI parallel
!! C functions from par-src/scans.
!!
!! This is a modern rewrite of the original code by Robert Ferrell
!! using the C interoperability features of Fortran 2003.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! November 2014
!!

subroutine off_node_sum_prefix_int_c (Dest_Data, Dest_Seg, Src_Data, Src_Seg) bind(c)
  ! Since this is serial emulation, there is no contribution from other PEs
  use,intrinsic :: iso_c_binding, only: c_int
  integer(c_int) :: Dest_Data, Src_Data
  integer(c_int) :: Dest_Seg, Src_Seg
  Dest_Data = 0
  Dest_Seg  = 0
end subroutine

subroutine off_node_sum_prefix_float_c (Dest_Data, Dest_Seg, Src_Data, Src_Seg) bind(c)
  ! Since this is serial emulation, there is no contribution from other PEs
  use,intrinsic :: iso_c_binding, only: c_int, c_float
  real(c_float) :: Dest_Data, Src_Data
  integer(c_int) :: Dest_Seg, Src_Seg
  Dest_Data = 0
  Dest_Seg  = 0
end subroutine

subroutine off_node_sum_prefix_double_c (Dest_Data, Dest_Seg, Src_Data, Src_Seg) bind(c)
  ! Since this is serial emulation, there is no contribution from other PEs
  use,intrinsic :: iso_c_binding, only: c_int, c_double
  real(c_double) :: Dest_Data, Src_Data
  integer(c_int) :: Dest_Seg, Src_Seg
  Dest_Data = 0
  Dest_Seg  = 0
end subroutine

subroutine off_node_sum_suffix_int_c (Dest_Data, Dest_Seg, Src_Data, Src_Seg) bind(c)
  ! Since this is serial emulation, there is no contribution from other PEs
  use,intrinsic :: iso_c_binding, only: c_int
  integer(c_int) :: Dest_Data, Src_Data
  integer(c_int) :: Dest_Seg, Src_Seg
  Dest_Data = 0
  Dest_Seg  = 0
end subroutine

subroutine off_node_sum_suffix_float_c (Dest_Data, Dest_Seg, Src_Data, Src_Seg) bind(c)
  ! Since this is serial emulation, there is no contribution from other PEs
  use,intrinsic :: iso_c_binding, only: c_int, c_float
  real(c_float) :: Dest_Data, Src_Data
  integer(c_int) :: Dest_Seg, Src_Seg
  Dest_Data = 0
  Dest_Seg  = 0
end subroutine

subroutine off_node_sum_suffix_double_c (Dest_Data, Dest_Seg, Src_Data, Src_Seg) bind(c)
  ! Since this is serial emulation, there is no contribution from other PEs
  use,intrinsic :: iso_c_binding, only: c_int, c_double
  real(c_double) :: Dest_Data, Src_Data
  integer(c_int) :: Dest_Seg, Src_Seg
  Dest_Data = 0
  Dest_Seg  = 0
end subroutine
