!!
!! FLOW_2D_OPERATORS_TYPE
!!
!! This module defines FLOW_2D_OPERATORS, geometry-dependent finite-volume
!! operators for two-dimensional flow on an unstructured mesh. The type has
!! no global state. Its first-order cell-to-face operators use a two-point
!! stencil in the face-normal direction. DIVERGENCE accumulates the signed
!! face fluxes in each on-process cell without volume scaling.
!!
!! NB: zjibben is working to implement the operators described in Ferzinger &
!! Peric's 2020 book, sections 9.7 and 9.8, including second-order
!! corrections for non-orthogonal meshes.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module flow_2d_operators_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_2d_mesh_type
  use bndry_func1_class
  use bndry_vfunc_class
  implicit none
  private

  type, public :: flow_2d_operators
    private
    type(unstr_2d_mesh), pointer :: mesh => null()  ! unowned reference
    real(r8), allocatable :: dx(:), dr(:,:,:), interpolation_factor(:)
  contains
    procedure :: init
    procedure :: derivative_cf_1r, derivative_cf_2r
    procedure :: interpolate_cf_1r, interpolate_cf_2r
    procedure :: gradient_cc
    procedure :: divergence
    procedure :: normal_distance
    generic :: derivative_cf => derivative_cf_1r, derivative_cf_2r
    generic :: interpolate_cf => interpolate_cf_1r, interpolate_cf_2r
  end type

contains

  subroutine init(this, mesh)
    class(flow_2d_operators), intent(out) :: this
    type(unstr_2d_mesh), target, intent(inout) :: mesh

    integer :: f, c1, c2
    real(r8) :: dcc(2), r(2,2)

    call mesh%init_cell_centroid()
    call mesh%init_face_centroid()
    this%mesh => mesh
    allocate(this%dx(mesh%nface), this%dr(2,2,mesh%nface_onP), &
        this%interpolation_factor(mesh%nface))
    this%dx = 0.0_r8
    this%interpolation_factor = 0.0_r8
    this%dr = 0.0_r8

    do f = 1, mesh%nface_onP
      c1 = mesh%fcell(1,f)
      c2 = mesh%fcell(2,f)
      if (c2 == 0) then
        this%dx(f) = dot_product(mesh%face_centroid(:,f) - mesh%cell_centroid(:,c1), &
            mesh%unit_normal(:,f))
        this%interpolation_factor(f) = 1.0_r8
        r(:,1) = mesh%face_centroid(:,f) - this%dx(f)*mesh%unit_normal(:,f)
        this%dr(:,1,f) = r(:,1) - mesh%cell_centroid(:,c1)
      else
        dcc = mesh%cell_centroid(:,c2) - mesh%cell_centroid(:,c1)
        this%dx(f) = 2.0_r8*min( &
            dot_product(mesh%face_centroid(:,f) - mesh%cell_centroid(:,c1), mesh%unit_normal(:,f)), &
            -dot_product(mesh%face_centroid(:,f) - mesh%cell_centroid(:,c2), mesh%unit_normal(:,f)))
        this%interpolation_factor(f) = 1.0_r8 - &
            dot_product(mesh%face_centroid(:,f) - mesh%cell_centroid(:,c1), dcc)/dot_product(dcc, dcc)
        r(:,1) = mesh%face_centroid(:,f) - 0.5_r8*this%dx(f)*mesh%unit_normal(:,f)
        r(:,2) = mesh%face_centroid(:,f) + 0.5_r8*this%dx(f)*mesh%unit_normal(:,f)
        this%dr(:,1,f) = r(:,1) - mesh%cell_centroid(:,c1)
        this%dr(:,2,f) = r(:,2) - mesh%cell_centroid(:,c2)
      end if
      ASSERT(this%dx(f) > 0.0_r8)
      ASSERT(this%interpolation_factor(f) >= 0.0_r8)
      ASSERT(this%interpolation_factor(f) <= 1.0_r8)
    end do
    call mesh%face_imap%gather_offp(this%dx)
    call mesh%face_imap%gather_offp(this%interpolation_factor)
  end subroutine


  !! Compute a cell-centered scalar gradient by a least-squares fit to the
  !! neighboring cell values. Pressure Dirichlet values contribute a fit row
  !! from the cell center to the face center; pressure Neumann values contribute
  !! a normal row of length NORMAL_DISTANCE. Boundary faces not listed by a
  !! condition are omitted from the fit.
  subroutine gradient_cc(this, field_cc, gradient_c, normal_flux_bc, dirichlet_bc, gravity_head)
    class(flow_2d_operators), intent(in) :: this
    real(r8), intent(in) :: field_cc(:)
    real(r8), intent(out) :: gradient_c(:,:)
    class(bndry_func1), optional, intent(in) :: normal_flux_bc, dirichlet_bc
    real(r8), optional, intent(in) :: gravity_head(:,:)

    integer :: c, i, f, neighbor, n
    real(r8) :: r(2), difference, matrix(2,2), rhs(2), determinant
    logical :: found

    ASSERT(size(field_cc) >= this%mesh%ncell)
    ASSERT(size(gradient_c,1) == 2)
    ASSERT(size(gradient_c,2) == this%mesh%ncell)
    gradient_c = 0.0_r8
    do c = 1, this%mesh%ncell_onP
      matrix = 0.0_r8
      rhs = 0.0_r8
      n = 0
      do i = this%mesh%cstart(c), this%mesh%cstart(c+1)-1
        f = this%mesh%cface(i)
        neighbor = this%mesh%cnhbr(i)
        if (neighbor > 0) then
          r = this%mesh%cell_centroid(:,neighbor) - this%mesh%cell_centroid(:,c)
          difference = field_cc(neighbor) - field_cc(c)
          if (present(gravity_head)) then
            if (c == this%mesh%fcell(1,f)) then
              difference = difference + gravity_head(2,f) - gravity_head(1,f)
            else
              difference = difference + gravity_head(1,f) - gravity_head(2,f)
            end if
          end if
        else
          found = .false.
          if (present(dirichlet_bc)) then
            call bndry_value(dirichlet_bc, f, difference, found)
            if (found) then
              r = this%mesh%face_centroid(:,f) - this%mesh%cell_centroid(:,c)
              if (present(gravity_head)) difference = difference - gravity_head(1,f)
              difference = difference - field_cc(c)
            end if
          end if
          if (.not.found .and. present(normal_flux_bc)) then
            call bndry_value(normal_flux_bc, f, difference, found)
            if (found) then
              r = this%normal_distance(f)*this%mesh%unit_normal(:,f)
              difference = difference*this%normal_distance(f)
            end if
          end if
          if (.not.found) cycle
        end if
        matrix = matrix + reshape([r(1)**2, r(1)*r(2), r(1)*r(2), r(2)**2], [2,2])
        rhs = rhs + r*difference
        n = n + 1
      end do
      ASSERT(n >= 2)
      determinant = matrix(1,1)*matrix(2,2) - matrix(1,2)*matrix(2,1)
      ASSERT(abs(determinant) > tiny(determinant))
      gradient_c(1,c) = (matrix(2,2)*rhs(1) - matrix(1,2)*rhs(2))/determinant
      gradient_c(2,c) = (matrix(1,1)*rhs(2) - matrix(2,1)*rhs(1))/determinant
    end do
    call this%mesh%cell_imap%gather_offp(gradient_c)
  end subroutine


  subroutine bndry_value(bndry, face, value, found)
    class(bndry_func1), intent(in) :: bndry
    integer, intent(in) :: face
    real(r8), intent(out) :: value
    logical, intent(out) :: found

    integer :: i

    i = findloc(bndry%index, face, dim=1)
    found = i > 0
    if (found) value = bndry%value(i)
  end subroutine


  function normal_distance(this, face) result(dx)
    class(flow_2d_operators), intent(in) :: this
    integer, intent(in) :: face
    real(r8) :: dx

    ASSERT(face >= 1 .and. face <= this%mesh%nface)
    dx = this%dx(face)
  end function


  !! Compute a first-order face-normal derivative from cell-centered scalar
  !! data. On a boundary face, NORMAL_FLUX_BC supplies the derivative directly
  !! and DIRICHLET_BC supplies the scalar value at the face center.
  subroutine derivative_cf_1r(this, field_cc, derivative_fn, normal_flux_bc, dirichlet_bc, &
      dirichlet_value, gravity_head)
    class(flow_2d_operators), intent(in) :: this
    real(r8), intent(in) :: field_cc(:)
    real(r8), intent(out) :: derivative_fn(:)
    class(bndry_func1), optional, intent(in) :: normal_flux_bc, dirichlet_bc
    real(r8), optional, intent(in) :: dirichlet_value(:)
    real(r8), optional, intent(in) :: gravity_head(:,:)

    integer :: f, c1, c2, i

    ASSERT(size(field_cc) >= this%mesh%ncell)
    ASSERT(size(derivative_fn) == this%mesh%nface)
    if (present(gravity_head)) then
      ASSERT(size(gravity_head,1) == 2)
      ASSERT(size(gravity_head,2) == this%mesh%nface)
    end if
    derivative_fn = 0.0_r8
    do f = 1, this%mesh%nface_onP
      c1 = this%mesh%fcell(1,f)
      c2 = this%mesh%fcell(2,f)
      if (c2 > 0) then
        derivative_fn(f) = field_cc(c2) - field_cc(c1)
        if (present(gravity_head)) derivative_fn(f) = derivative_fn(f) + &
            gravity_head(2,f) - gravity_head(1,f)
        derivative_fn(f) = derivative_fn(f)/this%dx(f)
      end if
    end do
    if (present(normal_flux_bc)) then
      do i = 1, size(normal_flux_bc%index)
        f = normal_flux_bc%index(i)
        if (f > this%mesh%nface_onP) cycle
        ASSERT(f >= 1)
        ASSERT(this%mesh%fcell(2,f) == 0)
        derivative_fn(f) = normal_flux_bc%value(i)
      end do
    end if
    if (present(dirichlet_bc)) then
      if (present(dirichlet_value)) then
        ASSERT(size(dirichlet_value) == size(dirichlet_bc%value))
      end if
      do i = 1, size(dirichlet_bc%index)
        f = dirichlet_bc%index(i)
        if (f > this%mesh%nface_onP) cycle
        ASSERT(f >= 1)
        c1 = this%mesh%fcell(1,f)
        ASSERT(this%mesh%fcell(2,f) == 0)
        if (present(dirichlet_value)) then
          derivative_fn(f) = dirichlet_value(i) - field_cc(c1)
        else
          derivative_fn(f) = dirichlet_bc%value(i) - field_cc(c1)
        end if
        if (present(gravity_head)) derivative_fn(f) = derivative_fn(f) - gravity_head(1,f)
        derivative_fn(f) = derivative_fn(f)/this%dx(f)
      end do
    end if
    call this%mesh%face_imap%gather_offp(derivative_fn)
  end subroutine


  !! Vector version of DERIVATIVE_CF_1R. The result is the face-normal
  !! derivative of each vector component. ZERO_NORMAL_BC imposes a zero
  !! normal velocity while retaining the adjacent cell's tangential velocity.
  subroutine derivative_cf_2r(this, field_cc, derivative_fn, zero_normal_bc, dirichlet_bc)
    class(flow_2d_operators), intent(in) :: this
    real(r8), intent(in) :: field_cc(:,:)
    real(r8), intent(out) :: derivative_fn(:,:)
    class(bndry_func1), optional, intent(in) :: zero_normal_bc
    class(bndry_vfunc), optional, intent(in) :: dirichlet_bc

    integer :: f, c1, c2, i
    real(r8) :: normal_velocity

    ASSERT(size(field_cc,1) == 2)
    ASSERT(size(field_cc,2) >= this%mesh%ncell)
    ASSERT(size(derivative_fn,1) == 2)
    ASSERT(size(derivative_fn,2) == this%mesh%nface)
    derivative_fn = 0.0_r8
    do f = 1, this%mesh%nface_onP
      c1 = this%mesh%fcell(1,f)
      c2 = this%mesh%fcell(2,f)
      if (c2 > 0) derivative_fn(:,f) = (field_cc(:,c2) - field_cc(:,c1))/this%dx(f)
    end do
    if (present(zero_normal_bc)) then
      do i = 1, size(zero_normal_bc%index)
        f = zero_normal_bc%index(i)
        if (f > this%mesh%nface_onP) cycle
        ASSERT(f >= 1)
        c1 = this%mesh%fcell(1,f)
        ASSERT(this%mesh%fcell(2,f) == 0)
        normal_velocity = dot_product(this%mesh%unit_normal(:,f), field_cc(:,c1))
        derivative_fn(:,f) = -normal_velocity*this%mesh%unit_normal(:,f)/this%dx(f)
      end do
    end if
    if (present(dirichlet_bc)) then
      do i = 1, size(dirichlet_bc%index)
        f = dirichlet_bc%index(i)
        if (f > this%mesh%nface_onP) cycle
        ASSERT(f >= 1)
        c1 = this%mesh%fcell(1,f)
        ASSERT(this%mesh%fcell(2,f) == 0)
        derivative_fn(:,f) = (dirichlet_bc%value(:,i) - field_cc(:,c1))/this%dx(f)
      end do
    end if
    call this%mesh%face_imap%gather_offp(derivative_fn)
  end subroutine


  !! Interpolate cell-centered scalar data to faces using the first-order
  !! method of Ferzinger & Peric (2020, Eq. 9.36). Boundary-face values are
  !! the adjacent cell values.
  subroutine interpolate_cf_1r(this, field_cc, field_f)
    class(flow_2d_operators), intent(in) :: this
    real(r8), intent(in) :: field_cc(:)
    real(r8), intent(out) :: field_f(:)

    integer :: f, c1, c2

    ASSERT(size(field_cc) >= this%mesh%ncell)
    ASSERT(size(field_f) == this%mesh%nface)
    do f = 1, this%mesh%nface_onP
      c1 = this%mesh%fcell(1,f)
      c2 = this%mesh%fcell(2,f)
      field_f(f) = this%interpolation_factor(f)*field_cc(c1)
      if (c2 > 0) field_f(f) = field_f(f) + (1.0_r8-this%interpolation_factor(f))*field_cc(c2)
    end do
    call this%mesh%face_imap%gather_offp(field_f)
  end subroutine


  !! Interpolate cell-centered velocity to its face-normal component using
  !! the same first-order interpolation as INTERPOLATE_CF_1R. DIRICHLET_BC
  !! and ZERO_NORMAL_BC override the values at their boundary faces.
  subroutine interpolate_cf_2r(this, field_cc, field_f, zero_normal_bc, dirichlet_bc)
    class(flow_2d_operators), intent(in) :: this
    real(r8), intent(in) :: field_cc(:,:)
    real(r8), intent(out) :: field_f(:)
    class(bndry_func1), optional, intent(in) :: zero_normal_bc
    class(bndry_vfunc), optional, intent(in) :: dirichlet_bc

    integer :: f, c1, c2, i

    ASSERT(size(field_cc,1) == 2)
    ASSERT(size(field_cc,2) >= this%mesh%ncell)
    ASSERT(size(field_f) == this%mesh%nface)
    do f = 1, this%mesh%nface_onP
      c1 = this%mesh%fcell(1,f)
      c2 = this%mesh%fcell(2,f)
      field_f(f) = this%interpolation_factor(f)*dot_product(this%mesh%unit_normal(:,f), field_cc(:,c1))
      if (c2 > 0) field_f(f) = field_f(f) + (1.0_r8-this%interpolation_factor(f))* &
          dot_product(this%mesh%unit_normal(:,f), field_cc(:,c2))
    end do
    if (present(zero_normal_bc)) then
      do i = 1, size(zero_normal_bc%index)
        f = zero_normal_bc%index(i)
        if (f > this%mesh%nface_onP) cycle
        ASSERT(f >= 1)
        ASSERT(this%mesh%fcell(2,f) == 0)
        field_f(f) = 0.0_r8
      end do
    end if
    if (present(dirichlet_bc)) then
      do i = 1, size(dirichlet_bc%index)
        f = dirichlet_bc%index(i)
        if (f > this%mesh%nface_onP) cycle
        ASSERT(f >= 1)
        ASSERT(this%mesh%fcell(2,f) == 0)
        field_f(f) = dot_product(this%mesh%unit_normal(:,f), dirichlet_bc%value(:,i))
      end do
    end if
    call this%mesh%face_imap%gather_offp(field_f)
  end subroutine


  !! Accumulate the signed face-normal fluxes in each on-process cell. The
  !! returned values are net fluxes, not volume-scaled divergence values.
  subroutine divergence(this, field_f, flux_c)
    class(flow_2d_operators), intent(in) :: this
    real(r8), intent(in) :: field_f(:)
    real(r8), intent(out) :: flux_c(:)

    integer :: c, i, f

    ASSERT(size(field_f) >= this%mesh%nface)
    ASSERT(size(flux_c) == this%mesh%ncell_onP)
    do c = 1, this%mesh%ncell_onP
      flux_c(c) = 0.0_r8
      do i = this%mesh%cstart(c), this%mesh%cstart(c+1)-1
        f = this%mesh%cface(i)
        if (btest(this%mesh%cfpar(c), i-this%mesh%cstart(c)+1)) then
          flux_c(c) = flux_c(c) - field_f(f)*this%mesh%area(f)
        else
          flux_c(c) = flux_c(c) + field_f(f)*this%mesh%area(f)
        end if
      end do
    end do
  end subroutine

end module flow_2d_operators_type
