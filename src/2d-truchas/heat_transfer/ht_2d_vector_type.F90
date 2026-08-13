!!
!! HT_2D_VECTOR_TYPE
!!
!! This module defines the concrete vector type used by the 2D heat-transfer
!! solver. It stores the cell enthalpy, cell temperature, and face temperature
!! unknowns using wide mesh-sized arrays whose trailing off-process entries
!! provide persistent workspace for model assembly routines.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!
!! Notes
!!
!! The abstract vector operations act only on the on-process portions of the
!! cell and face arrays. Off-process entries are scratch values and must be
!! gathered explicitly by model code before they are used.
!!
!! Vector arguments use intent(inout) even when the operation is conceptually
!! data-only output. The vector object carries structural metadata, such as
!! mesh and component layout, that remains input to the operation while the
!! numeric component data may be overwritten.
!!

module ht_2d_vector_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use vector_class
  use unstr_2d_mesh_type
  use parallel_communication, only: global_sum, global_maxval
  implicit none
  private

  type, extends(vector), public :: ht_2d_vector
    type(unstr_2d_mesh), pointer :: mesh => null() ! unowned reference
    real(r8), allocatable :: hc(:) ! cell enthalpy density
    real(r8), allocatable :: tc(:) ! cell temperature
    real(r8), allocatable :: tf(:) ! face temperature
  contains
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
    generic :: init => init_mesh, init_mold
    procedure, private :: init_mesh, init_mold
    procedure :: gather_offp
  end type

contains

  subroutine init_mesh(this, mesh)
    class(ht_2d_vector), intent(out) :: this
    type(unstr_2d_mesh), intent(in), target :: mesh
    this%mesh => mesh
    allocate(this%hc(mesh%ncell), this%tc(mesh%ncell), this%tf(mesh%nface))
  end subroutine

  subroutine init_mold(this, mold)
    class(ht_2d_vector), intent(out) :: this
    class(ht_2d_vector), intent(in) :: mold
    call this%init(mold%mesh)
  end subroutine

  subroutine gather_offp(this)
    class(ht_2d_vector), intent(inout) :: this
    call this%mesh%cell_imap%gather_offp(this%hc)
    call this%mesh%cell_imap%gather_offp(this%tc)
    call this%mesh%face_imap%gather_offp(this%tf)
  end subroutine

  subroutine clone1(this, clone)
    class(ht_2d_vector), intent(in) :: this
    class(vector), allocatable, intent(out) :: clone
    allocate(ht_2d_vector :: clone)
    select type (clone)
    class is (ht_2d_vector)
      call clone%init(this)
    end select
  end subroutine

  subroutine clone2(this, clone, n)
    class(ht_2d_vector), intent(in) :: this
    class(vector), allocatable, intent(out) :: clone(:)
    integer, intent(in) :: n
    integer :: j
    allocate(ht_2d_vector :: clone(n))
    select type (clone)
    class is (ht_2d_vector)
      do j = 1, n
        call clone(j)%init(this)
      end do
    end select
  end subroutine

  subroutine copy_(dest, src)
    class(ht_2d_vector), intent(inout) :: dest
    class(vector), intent(in) :: src
    select type (src)
    class is (ht_2d_vector)
      associate (ncell_onP => dest%mesh%ncell_onP, nface_onP => dest%mesh%nface_onP)
        dest%hc(:ncell_onP) = src%hc(:ncell_onP)
        dest%tc(:ncell_onP) = src%tc(:ncell_onP)
        dest%tf(:nface_onP) = src%tf(:nface_onP)
      end associate
    end select
  end subroutine

  subroutine setval(this, val)
    class(ht_2d_vector), intent(inout) :: this
    real(r8), intent(in) :: val
    associate (ncell_onP => this%mesh%ncell_onP, nface_onP => this%mesh%nface_onP)
      this%hc(:ncell_onP) = val
      this%tc(:ncell_onP) = val
      this%tf(:nface_onP) = val
    end associate
  end subroutine

  subroutine scale(this, a)
    class(ht_2d_vector), intent(inout) :: this
    real(r8), intent(in) :: a
    associate (ncell_onP => this%mesh%ncell_onP, nface_onP => this%mesh%nface_onP)
      this%hc(:ncell_onP) = a*this%hc(:ncell_onP)
      this%tc(:ncell_onP) = a*this%tc(:ncell_onP)
      this%tf(:nface_onP) = a*this%tf(:nface_onP)
    end associate
  end subroutine

  subroutine update1_(this, a, x)
    class(ht_2d_vector), intent(inout) :: this
    class(vector), intent(in) :: x
    real(r8), intent(in) :: a
    select type (x)
    class is (ht_2d_vector)
      associate (ncell_onP => this%mesh%ncell_onP, nface_onP => this%mesh%nface_onP)
        this%hc(:ncell_onP) = a*x%hc(:ncell_onP) + this%hc(:ncell_onP)
        this%tc(:ncell_onP) = a*x%tc(:ncell_onP) + this%tc(:ncell_onP)
        this%tf(:nface_onP) = a*x%tf(:nface_onP) + this%tf(:nface_onP)
      end associate
    end select
  end subroutine

  subroutine update2_(this, a, x, b)
    class(ht_2d_vector), intent(inout) :: this
    class(vector), intent(in) :: x
    real(r8), intent(in) :: a, b
    select type (x)
    class is (ht_2d_vector)
      associate (ncell_onP => this%mesh%ncell_onP, nface_onP => this%mesh%nface_onP)
        this%hc(:ncell_onP) = a*x%hc(:ncell_onP) + b*this%hc(:ncell_onP)
        this%tc(:ncell_onP) = a*x%tc(:ncell_onP) + b*this%tc(:ncell_onP)
        this%tf(:nface_onP) = a*x%tf(:nface_onP) + b*this%tf(:nface_onP)
      end associate
    end select
  end subroutine

  subroutine update3_(this, a, x, b, y)
    class(ht_2d_vector), intent(inout) :: this
    class(vector), intent(in) :: x, y
    real(r8), intent(in) :: a, b
    select type (x)
    class is (ht_2d_vector)
      select type (y)
      class is (ht_2d_vector)
        associate (ncell_onP => this%mesh%ncell_onP, nface_onP => this%mesh%nface_onP)
          this%hc(:ncell_onP) = a*x%hc(:ncell_onP) + b*y%hc(:ncell_onP) + this%hc(:ncell_onP)
          this%tc(:ncell_onP) = a*x%tc(:ncell_onP) + b*y%tc(:ncell_onP) + this%tc(:ncell_onP)
          this%tf(:nface_onP) = a*x%tf(:nface_onP) + b*y%tf(:nface_onP) + this%tf(:nface_onP)
        end associate
      end select
    end select
  end subroutine

  subroutine update4_(this, a, x, b, y, c)
    class(ht_2d_vector), intent(inout) :: this
    class(vector), intent(in) :: x, y
    real(r8), intent(in) :: a, b, c
    select type (x)
    class is (ht_2d_vector)
      select type (y)
      class is (ht_2d_vector)
        associate (ncell_onP => this%mesh%ncell_onP, nface_onP => this%mesh%nface_onP)
          this%hc(:ncell_onP) = a*x%hc(:ncell_onP) + b*y%hc(:ncell_onP) + c*this%hc(:ncell_onP)
          this%tc(:ncell_onP) = a*x%tc(:ncell_onP) + b*y%tc(:ncell_onP) + c*this%tc(:ncell_onP)
          this%tf(:nface_onP) = a*x%tf(:nface_onP) + b*y%tf(:nface_onP) + c*this%tf(:nface_onP)
        end associate
      end select
    end select
  end subroutine

  function dot_(x, y) result(dp)
    class(ht_2d_vector), intent(in) :: x
    class(vector), intent(in) :: y
    real(r8) :: dp
    select type (y)
    class is (ht_2d_vector)
      dp = dot_product(x%hc(:x%mesh%ncell_onP), y%hc(:x%mesh%ncell_onP)) &
         + dot_product(x%tc(:x%mesh%ncell_onP), y%tc(:x%mesh%ncell_onP)) &
         + dot_product(x%tf(:x%mesh%nface_onP), y%tf(:x%mesh%nface_onP))
      dp = global_sum(dp)
    end select
  end function

  real(r8) function norm1_(this) result(norm)
    class(ht_2d_vector), intent(in) :: this
    norm = sum(abs(this%hc(:this%mesh%ncell_onP))) + sum(abs(this%tc(:this%mesh%ncell_onP))) &
         + sum(abs(this%tf(:this%mesh%nface_onP)))
    norm = global_sum(norm)
  end function

  real(r8) function norm2_(this) result(norm)
    class(ht_2d_vector), intent(in) :: this
    norm = norm2(this%hc(:this%mesh%ncell_onP))**2 + norm2(this%tc(:this%mesh%ncell_onP))**2 &
         + norm2(this%tf(:this%mesh%nface_onP))**2
    norm = sqrt(global_sum(norm))
  end function

  real(r8) function norm_max_(this) result(norm)
    class(ht_2d_vector), intent(in) :: this
    norm = max(0.0_r8, maxval(abs(this%hc(:this%mesh%ncell_onP))), &
               maxval(abs(this%tc(:this%mesh%ncell_onP))), maxval(abs(this%tf(:this%mesh%nface_onP))))
    norm = global_maxval(norm)
  end function

  function checksum(this, full) result(string)
    use md5_hash_type
    class(ht_2d_vector), intent(in) :: this
    logical, intent(in), optional :: full
    character(:), allocatable :: string
    type(md5_hash) :: hash
    logical :: strict
    strict = .true.
    if (present(full)) strict = .not.full
    if (strict) then
      call hash%update(this%hc(:this%mesh%ncell_onP))
      call hash%update(this%tc(:this%mesh%ncell_onP))
      call hash%update(this%tf(:this%mesh%nface_onP))
    else
      call hash%update(this%hc)
      call hash%update(this%tc)
      call hash%update(this%tf)
    end if
    string = hash%hexdigest()
  end function

end module ht_2d_vector_type
