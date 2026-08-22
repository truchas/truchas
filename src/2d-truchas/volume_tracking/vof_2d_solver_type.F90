!!
!! VOF_2D_SOLVER_TYPE
!!
!! This module defines VOF_2D_SOLVER, which advances a set of cell-centered
!! volume fractions with a prescribed two-dimensional velocity field.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module vof_2d_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_2d_mesh_type
  use vector_func_class
  use volume_tracker_2d_class
  implicit none
  private

  type, public :: vof_2d_solver
    private
    type(unstr_2d_mesh), pointer :: mesh => null()
    class(volume_tracker_2d), allocatable :: tracker
    real(r8), allocatable :: vfrac_in(:,:), vfrac_out(:,:), flux_volume(:,:)
    real(r8), allocatable :: flux_velocity(:), face_velocity(:), interface_normal(:,:,:)
  contains
    final :: delete
    procedure :: init
    procedure :: step
  end type

contains

  subroutine delete(this)
    type(vof_2d_solver), intent(inout) :: this
    nullify(this%mesh)
  end subroutine


  subroutine init(this, mesh, nmat, algorithm, axisymmetric, stat, errmsg)

    use simple_volume_tracker_type
    use geometric_volume_tracker_type

    class(vof_2d_solver), intent(out) :: this
    type(unstr_2d_mesh), target, intent(in) :: mesh
    integer, intent(in) :: nmat
    character(*), intent(in) :: algorithm
    logical, intent(in) :: axisymmetric
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    stat = 0
    if (nmat < 1) then
      stat = 1
      errmsg = 'require at least one material'
      return
    end if
    select case (algorithm)
    case ('simple')
      allocate(simple_volume_tracker :: this%tracker)
    case ('geometric')
      allocate(geometric_volume_tracker :: this%tracker)
    case default
      stat = 1
      errmsg = 'unknown volume-tracking algorithm: ' // algorithm
      return
    end select

    this%mesh => mesh
    allocate(this%vfrac_in(nmat,mesh%ncell), this%vfrac_out(nmat,mesh%ncell))
    allocate(this%flux_volume(nmat,size(mesh%cface)), this%flux_velocity(size(mesh%cface)))
    allocate(this%face_velocity(mesh%nface), this%interface_normal(2,nmat,mesh%ncell))
    call this%tracker%init(mesh, nmat, nmat, nmat, axisymmetric)
  end subroutine


  subroutine step(this, time, dt, velocity, vfrac)

    class(vof_2d_solver), intent(inout) :: this
    real(r8), intent(in) :: time, dt
    class(vector_func), intent(in) :: velocity
    real(r8), intent(inout) :: vfrac(:,:)

    integer :: c, f, j, f0, f1
    real(r8) :: v(2)

    ASSERT(associated(this%mesh))
    ASSERT(velocity%dim == 2)
    ASSERT(all(shape(vfrac) == [size(this%vfrac_in,1), this%mesh%ncell]))

    do f = 1, this%mesh%nface
      v = velocity%eval([time, this%mesh%face_centroid(:,f)])
      this%face_velocity(f) = dot_product(v, this%mesh%unit_normal(:,f))
    end do

    this%vfrac_in = vfrac
    this%flux_volume = 0.0_r8
    do c = 1, this%mesh%ncell
      f0 = this%mesh%cstart(c)
      f1 = this%mesh%cstart(c+1)-1
      do j = f0, f1
        f = this%mesh%cface(j)
        if (btest(this%mesh%cfpar(c), 1+j-f0)) then
          this%flux_velocity(j) = -this%face_velocity(f)
        else
          this%flux_velocity(j) = this%face_velocity(f)
        end if
      end do
    end do
    call this%tracker%flux_volumes(this%flux_velocity, this%vfrac_in, this%vfrac_out, &
        this%flux_volume, this%interface_normal, size(vfrac,1), 0, dt)
    vfrac = this%vfrac_out
  end subroutine

end module vof_2d_solver_type
