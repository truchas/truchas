! -*- Mode: f90 -*-

#include "matvec_proto.F"
#include "precond_proto.F"

    optional :: ApplyPrecond
    real(real_type), intent(in), dimension(:) :: b
    type(control_type), intent(inout) :: control
    real(real_type), intent(inout), dimension(:), target :: x
    ! <arguments>

    ! Local:
    character(len=32) :: routine
    character(len=132) :: out_string
    logical(logical_type) :: outmode_iterates_, outmode_full_
    logical(logical_type) :: precond, stopping_uses_x, update_iterate
    integer(int_type) :: i, j, k, m, n, nunk, luerr_, luout_, norm_type_, status
    real(real_type), allocatable, dimension(:) :: c, s, r, rt, w
    real(real_type), allocatable, dimension(:,:) :: h
    real(real_type), allocatable, dimension(:), target :: xt
    type(vector_type) :: x_vec, xt_vec
