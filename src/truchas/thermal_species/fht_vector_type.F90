!!
!! FHT_VECTOR_TYPE
!!
!! This module defines the concrete vector type used by the flow-coupled heat
!! transport solver. It stores the cell and face temperature unknowns using wide
!! mesh-sized arrays whose trailing off-process entries provide persistent
!! workspace for model assembly routines. When view-factor radiation is active,
!! it also stores the owned radiative flux unknowns.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!
!! Notes
!!
!! The abstract vector operations act only on the on-process portions of the
!! cell and face arrays. Off-process entries are scratch values and must be
!! gathered explicitly by model code before they are used. The radiative flux
!! array has no on/off-process split; its entries are uniquely owned.
!!
!! Vector arguments use intent(inout) even when the operation is conceptually
!! data-only output. The vector object carries structural metadata, such as
!! mesh and component layout, that remains input to the operation while the
!! numeric component data may be overwritten.
!!

module fht_vector_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use vector_class
  use unstr_mesh_type
  use parallel_communication, only: global_sum, global_maxval
  implicit none
  private

  type, extends(vector), public :: fht_vector
    type(unstr_mesh), pointer :: mesh => null() ! unowned reference
    real(r8), allocatable :: tc(:), tf(:)
    real(r8), allocatable :: qrad(:)
  contains
    !! Deferred base class procedures
    procedure :: clone1
    procedure :: clone2
    procedure :: copy_
    procedure :: setval
    procedure :: scale
    procedure :: update1_
    procedure :: update2_
    procedure :: update3_
    procedure :: update4_
    procedure :: dot_
    procedure :: norm1 => norm1_
    procedure :: norm2 => norm2_
    procedure :: norm_max => norm_max_
    procedure :: checksum
    !! Additional procedures specific to this type
    generic :: init => init_mesh, init_mold
    procedure, private :: init_mesh, init_mold
    procedure :: gather_offp
    procedure :: set_cell_temp
    procedure :: set_face_temp
  end type

contains

  !! Specific subroutine for the generic INIT. Initialize a FHT_VECTOR object
  !! for the given unstructured MESH. The numeric elements are not initialized.

  subroutine init_mesh(this, mesh, qrad_size)
    class(fht_vector), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    integer, intent(in), optional :: qrad_size
    this%mesh => mesh
    allocate(this%tc(mesh%ncell), this%tf(mesh%nface))
    if (present(qrad_size)) allocate(this%qrad(qrad_size))
  end subroutine

  !! Specific subroutine for the generic INIT. Initialize a FHT_VECTOR object
  !! to be a clone of MOLD. The numeric elements are not initialized.

  subroutine init_mold(this, mold)
    class(fht_vector), intent(out) :: this
    class(fht_vector), intent(in) :: mold
    if (allocated(mold%qrad)) then
      call this%init(mold%mesh, size(mold%qrad))
    else
      call this%init(mold%mesh)
    end if
  end subroutine

  subroutine gather_offp(this)
    class(fht_vector), intent(inout) :: this
    call this%mesh%cell_imap%gather_offp(this%tc)
    call this%mesh%face_imap%gather_offp(this%tf)
  end subroutine

  subroutine set_cell_temp(this, cell, value)
    class(fht_vector), intent(inout) :: this
    integer, intent(in) :: cell
    real(r8), intent(in) :: value
    this%tc(cell) = value
  end subroutine

  subroutine set_face_temp(this, face, value)
    class(fht_vector), intent(inout) :: this
    integer, intent(in) :: face
    real(r8), intent(in) :: value
    this%tf(face) = value
  end subroutine

  subroutine clone1(this, clone)
    class(fht_vector), intent(in) :: this
    class(vector), allocatable, intent(out) :: clone
    allocate(fht_vector :: clone)
    select type (clone)
    class is (fht_vector)
      call clone%init(this)
    end select
  end subroutine

  subroutine clone2(this, clone, n)
    class(fht_vector), intent(in) :: this
    class(vector), allocatable, intent(out) :: clone(:)
    integer, intent(in) :: n
    integer :: j
    allocate(fht_vector :: clone(n))
    select type (clone)
    class is (fht_vector)
      do j = 1, n
        call clone(j)%init(this)
      end do
    end select
  end subroutine

  subroutine copy_(dest, src)
    class(fht_vector), intent(inout) :: dest
    class(vector), intent(in) :: src
    select type (src)
    class is (fht_vector)
      associate (ncell_onP => dest%mesh%ncell_onP, nface_onP => dest%mesh%nface_onP)
        dest%tc(:ncell_onP) = src%tc(:ncell_onP)
        dest%tf(:nface_onP) = src%tf(:nface_onP)
        if (allocated(dest%qrad)) dest%qrad(:) = src%qrad
      end associate
    end select
  end subroutine

  subroutine setval(this, val)
    class(fht_vector), intent(inout) :: this
    real(r8), intent(in) :: val
    associate (ncell_onP => this%mesh%ncell_onP, nface_onP => this%mesh%nface_onP)
      this%tc(:ncell_onP) = val
      this%tf(:nface_onP) = val
      if (allocated(this%qrad)) this%qrad(:) = val
    end associate
  end subroutine

  subroutine scale(this, a)
    class(fht_vector), intent(inout) :: this
    real(r8), intent(in) :: a
    associate (ncell_onP => this%mesh%ncell_onP, nface_onP => this%mesh%nface_onP)
      this%tc(:ncell_onP) = a * this%tc(:ncell_onP)
      this%tf(:nface_onP) = a * this%tf(:nface_onP)
      if (allocated(this%qrad)) this%qrad(:) = a * this%qrad
    end associate
  end subroutine

  subroutine update1_(this, a, x)
    class(fht_vector), intent(inout) :: this
    class(vector), intent(in) :: x
    real(r8), intent(in) :: a
    select type (x)
    class is (fht_vector)
      associate (ncell_onP => this%mesh%ncell_onP, nface_onP => this%mesh%nface_onP)
        this%tc(:ncell_onP) = a * x%tc(:ncell_onP) + this%tc(:ncell_onP)
        this%tf(:nface_onP) = a * x%tf(:nface_onP) + this%tf(:nface_onP)
        if (allocated(this%qrad)) this%qrad(:) = a * x%qrad + this%qrad
      end associate
    end select
  end subroutine

  subroutine update2_(this, a, x, b)
    class(fht_vector), intent(inout) :: this
    class(vector), intent(in) :: x
    real(r8), intent(in) :: a, b
    select type (x)
    class is (fht_vector)
      associate (ncell_onP => this%mesh%ncell_onP, nface_onP => this%mesh%nface_onP)
        this%tc(:ncell_onP) = a * x%tc(:ncell_onP) + b * this%tc(:ncell_onP)
        this%tf(:nface_onP) = a * x%tf(:nface_onP) + b * this%tf(:nface_onP)
        if (allocated(this%qrad)) this%qrad(:) = a * x%qrad + b * this%qrad
      end associate
    end select
  end subroutine

  subroutine update3_(this, a, x, b, y)
    class(fht_vector), intent(inout) :: this
    class(vector), intent(in) :: x, y
    real(r8), intent(in) :: a, b
    select type (x)
    class is (fht_vector)
      select type (y)
      class is (fht_vector)
        associate (ncell_onP => this%mesh%ncell_onP, nface_onP => this%mesh%nface_onP)
          this%tc(:ncell_onP) = a * x%tc(:ncell_onP) + b * y%tc(:ncell_onP) &
                              + this%tc(:ncell_onP)
          this%tf(:nface_onP) = a * x%tf(:nface_onP) + b * y%tf(:nface_onP) &
                              + this%tf(:nface_onP)
          if (allocated(this%qrad)) this%qrad(:) = a * x%qrad + b * y%qrad + this%qrad
        end associate
      end select
    end select
  end subroutine

  subroutine update4_(this, a, x, b, y, c)
    class(fht_vector), intent(inout) :: this
    class(vector), intent(in) :: x, y
    real(r8), intent(in) :: a, b, c
    select type (x)
    class is (fht_vector)
      select type (y)
      class is (fht_vector)
        associate (ncell_onP => this%mesh%ncell_onP, nface_onP => this%mesh%nface_onP)
          this%tc(:ncell_onP) = a * x%tc(:ncell_onP) + b * y%tc(:ncell_onP) &
                              + c * this%tc(:ncell_onP)
          this%tf(:nface_onP) = a * x%tf(:nface_onP) + b * y%tf(:nface_onP) &
                              + c * this%tf(:nface_onP)
          if (allocated(this%qrad)) this%qrad(:) = a * x%qrad + b * y%qrad + c * this%qrad
        end associate
      end select
    end select
  end subroutine

  function dot_(x, y) result(dp)
    class(fht_vector), intent(in) :: x
    class(vector), intent(in) :: y
    real(r8) :: dp
    select type (y)
    class is (fht_vector)
      associate (ncell_onP => x%mesh%ncell_onP, nface_onP => x%mesh%nface_onP)
        dp = dot_product(x%tc(:ncell_onP), y%tc(:ncell_onP)) &
           + dot_product(x%tf(:nface_onP), y%tf(:nface_onP))
      end associate
      if (allocated(x%qrad)) dp = dp + dot_product(x%qrad, y%qrad)
      dp = global_sum(dp)
    end select
  end function

  real(r8) function norm1_(this) result(norm)
    class(fht_vector), intent(in) :: this
    associate (ncell_onP => this%mesh%ncell_onP, nface_onP => this%mesh%nface_onP)
      norm = sum(abs(this%tc(:ncell_onP))) &
           + sum(abs(this%tf(:nface_onP)))
    end associate
    if (allocated(this%qrad)) norm = norm + sum(abs(this%qrad))
    norm = global_sum(norm)
  end function

  real(r8) function norm2_(this) result(norm)
    class(fht_vector), intent(in) :: this
    associate (ncell_onP => this%mesh%ncell_onP, nface_onP => this%mesh%nface_onP)
      norm = norm2(this%tc(:ncell_onP))**2 &
           + norm2(this%tf(:nface_onP))**2
    end associate
    if (allocated(this%qrad)) norm = norm + norm2(this%qrad)**2
    norm = sqrt(global_sum(norm))
  end function

  real(r8) function norm_max_(this) result(norm)
    class(fht_vector), intent(in) :: this
    associate (ncell_onP => this%mesh%ncell_onP, nface_onP => this%mesh%nface_onP)
      norm = max(0.0_r8, &
                 maxval(abs(this%tc(:ncell_onP))), &
                 maxval(abs(this%tf(:nface_onP))))
    end associate
    if (allocated(this%qrad)) norm = max(norm, maxval(abs(this%qrad)))
    norm = global_maxval(norm)
  end function

  function checksum(this, full) result(string)
    use md5_hash_type
    class(fht_vector), intent(in) :: this
    logical, intent(in), optional :: full
    character(:), allocatable :: string

    type(md5_hash) :: hash
    logical :: strict
    strict = .true.
    if (present(full)) strict = .not.full
    if (strict) then
      associate (ncell_onP => this%mesh%ncell_onP, nface_onP => this%mesh%nface_onP)
        call hash%update(this%tc(:ncell_onP))
        call hash%update(this%tf(:nface_onP))
      end associate
    else
      call hash%update(this%tc(:this%mesh%ncell))
      call hash%update(this%tf(:this%mesh%nface))
    end if
    if (allocated(this%qrad)) call hash%update(this%qrad)
    string = hash%hexdigest()
  end function

end module fht_vector_type
