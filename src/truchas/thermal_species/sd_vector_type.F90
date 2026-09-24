!!
!! SD_VECTOR_TYPE
!!
!! This module defines the concrete vector type used by the species transport
!! solver. It stores the cell and face concentration unknowns for all species
!! components, using wide mesh-sized arrays whose trailing off-process entries
!! provide persistent workspace for model assembly routines.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!
!! Notes
!!
!! The abstract vector operations act only on the on-process portions of
!! those arrays. Off-process entries are scratch values and must be gathered
!! explicitly by model code before they are used.
!!
!! Vector arguments use intent(inout) even when the operation is conceptually
!! data-only output. The vector object carries structural metadata, such as
!! mesh and component layout, that remains input to the operation while the
!! numeric component data may be overwritten.
!!

module sd_vector_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use vector_class
  use unstr_mesh_type
  use parallel_communication, only: global_sum, global_maxval
  implicit none
  private

  type, extends(vector), public :: sd_vector
    type(unstr_mesh), pointer :: mesh => null() ! unowned reference
    real(r8), allocatable :: cc(:,:) ! cell concentration, including off-process cells
    real(r8), allocatable :: cf(:,:) ! face concentration, including off-process faces
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
  end type

contains

  !! Specific subroutine for the generic INIT. Initialize a SD_VECTOR object
  !! for the given unstructured MESH. The numeric elements are not initialized.

  subroutine init_mesh(this, mesh, num_comp)
    class(sd_vector), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    integer, intent(in) :: num_comp
    this%mesh => mesh
    allocate(this%cc(mesh%ncell,num_comp), this%cf(mesh%nface,num_comp))
  end subroutine

  !! Specific subroutine for the generic INIT. Initialize a SD_VECTOR object
  !! to be a clone of MOLD. The numeric elements are not initialized.

  subroutine init_mold(this, mold)
    class(sd_vector), intent(out) :: this
    class(sd_vector), intent(in) :: mold
    call this%init(mold%mesh, size(mold%cc,2))
  end subroutine

  subroutine gather_offp(this)
    class(sd_vector), intent(inout) :: this
    integer :: n
    do n = 1, size(this%cc,2)
      call this%mesh%cell_imap%gather_offp(this%cc(:,n))
      call this%mesh%face_imap%gather_offp(this%cf(:,n))
    end do
  end subroutine

  subroutine clone1(this, clone)
    class(sd_vector), intent(in) :: this
    class(vector), allocatable, intent(out) :: clone
    allocate(sd_vector :: clone)
    select type (clone)
    class is (sd_vector)
      call clone%init(this)
    end select
  end subroutine

  subroutine clone2(this, clone, n)
    class(sd_vector), intent(in) :: this
    class(vector), allocatable, intent(out) :: clone(:)
    integer, intent(in) :: n
    integer :: j
    allocate(sd_vector :: clone(n))
    select type (clone)
    class is (sd_vector)
      do j = 1, n
        call clone(j)%init(this)
      end do
    end select
  end subroutine

  subroutine copy_(dest, src)
    class(sd_vector), intent(inout) :: dest
    class(vector), intent(in) :: src
    select type (src)
    class is (sd_vector)
      associate (ncell_onP => dest%mesh%ncell_onP, nface_onP => dest%mesh%nface_onP)
        dest%cc(:ncell_onP,:) = src%cc(:ncell_onP,:)
        dest%cf(:nface_onP,:) = src%cf(:nface_onP,:)
      end associate
    end select
  end subroutine

  subroutine setval(this, val)
    class(sd_vector), intent(inout) :: this
    real(r8), intent(in) :: val
    associate (ncell_onP => this%mesh%ncell_onP, nface_onP => this%mesh%nface_onP)
      this%cc(:ncell_onP,:) = val
      this%cf(:nface_onP,:) = val
    end associate
  end subroutine

  subroutine scale(this, a)
    class(sd_vector), intent(inout) :: this
    real(r8), intent(in) :: a
    associate (ncell_onP => this%mesh%ncell_onP, nface_onP => this%mesh%nface_onP)
      this%cc(:ncell_onP,:) = a * this%cc(:ncell_onP,:)
      this%cf(:nface_onP,:) = a * this%cf(:nface_onP,:)
    end associate
  end subroutine

  subroutine update1_(this, a, x)
    class(sd_vector), intent(inout) :: this
    class(vector), intent(in) :: x
    real(r8), intent(in) :: a
    select type (x)
    class is (sd_vector)
      associate (ncell_onP => this%mesh%ncell_onP, nface_onP => this%mesh%nface_onP)
        this%cc(:ncell_onP,:) = a * x%cc(:ncell_onP,:) + this%cc(:ncell_onP,:)
        this%cf(:nface_onP,:) = a * x%cf(:nface_onP,:) + this%cf(:nface_onP,:)
      end associate
    end select
  end subroutine

  subroutine update2_(this, a, x, b)
    class(sd_vector), intent(inout) :: this
    class(vector), intent(in) :: x
    real(r8), intent(in) :: a, b
    select type (x)
    class is (sd_vector)
      associate (ncell_onP => this%mesh%ncell_onP, nface_onP => this%mesh%nface_onP)
        this%cc(:ncell_onP,:) = a * x%cc(:ncell_onP,:) + b * this%cc(:ncell_onP,:)
        this%cf(:nface_onP,:) = a * x%cf(:nface_onP,:) + b * this%cf(:nface_onP,:)
      end associate
    end select
  end subroutine

  subroutine update3_(this, a, x, b, y)
    class(sd_vector), intent(inout) :: this
    class(vector), intent(in) :: x, y
    real(r8), intent(in) :: a, b
    select type (x)
    class is (sd_vector)
      select type (y)
      class is (sd_vector)
        associate (ncell_onP => this%mesh%ncell_onP, nface_onP => this%mesh%nface_onP)
          this%cc(:ncell_onP,:) = a * x%cc(:ncell_onP,:) + b * y%cc(:ncell_onP,:) &
                                + this%cc(:ncell_onP,:)
          this%cf(:nface_onP,:) = a * x%cf(:nface_onP,:) + b * y%cf(:nface_onP,:) &
                                + this%cf(:nface_onP,:)
        end associate
      end select
    end select
  end subroutine

  subroutine update4_(this, a, x, b, y, c)
    class(sd_vector), intent(inout) :: this
    class(vector), intent(in) :: x, y
    real(r8), intent(in) :: a, b, c
    select type (x)
    class is (sd_vector)
      select type (y)
      class is (sd_vector)
        associate (ncell_onP => this%mesh%ncell_onP, nface_onP => this%mesh%nface_onP)
          this%cc(:ncell_onP,:) = a * x%cc(:ncell_onP,:) + b * y%cc(:ncell_onP,:) &
                                + c * this%cc(:ncell_onP,:)
          this%cf(:nface_onP,:) = a * x%cf(:nface_onP,:) + b * y%cf(:nface_onP,:) &
                                + c * this%cf(:nface_onP,:)
        end associate
      end select
    end select
  end subroutine

  function dot_(x, y) result(dp)
    class(sd_vector), intent(in) :: x
    class(vector), intent(in) :: y
    real(r8) :: dp
    integer :: n
    select type (y)
    class is (sd_vector)
      associate (ncell_onP => x%mesh%ncell_onP, nface_onP => x%mesh%nface_onP)
        dp = 0.0_r8
        do n = 1, size(x%cc,2)
          dp = dp + dot_product(x%cc(:ncell_onP,n), y%cc(:ncell_onP,n)) &
                  + dot_product(x%cf(:nface_onP,n), y%cf(:nface_onP,n))
        end do
      end associate
      dp = global_sum(dp)
    end select
  end function

  real(r8) function norm1_(this) result(norm)
    class(sd_vector), intent(in) :: this
    integer :: n
    associate (ncell_onP => this%mesh%ncell_onP, nface_onP => this%mesh%nface_onP)
      norm = 0.0_r8
      do n = 1, size(this%cc,2)
        norm = norm + sum(abs(this%cc(:ncell_onP,n))) &
                    + sum(abs(this%cf(:nface_onP,n)))
      end do
    end associate
    norm = global_sum(norm)
  end function

  real(r8) function norm2_(this) result(norm)
    class(sd_vector), intent(in) :: this
    integer :: n
    associate (ncell_onP => this%mesh%ncell_onP, nface_onP => this%mesh%nface_onP)
      norm = 0.0_r8
      do n = 1, size(this%cc,2)
        norm = norm + norm2(this%cc(:ncell_onP,n))**2 &
                    + norm2(this%cf(:nface_onP,n))**2
      end do
    end associate
    norm = sqrt(global_sum(norm))
  end function

  real(r8) function norm_max_(this) result(norm)
    class(sd_vector), intent(in) :: this
    integer :: n
    associate (ncell_onP => this%mesh%ncell_onP, nface_onP => this%mesh%nface_onP)
      norm = 0.0_r8
      do n = 1, size(this%cc,2)
        norm = max(norm, &
                   maxval(abs(this%cc(:ncell_onP,n))), &
                   maxval(abs(this%cf(:nface_onP,n))))
      end do
    end associate
    norm = global_maxval(norm)
  end function

  function checksum(this, full) result(string)
    use md5_hash_type
    class(sd_vector), intent(in) :: this
    logical, intent(in), optional :: full
    character(:), allocatable :: string
    integer :: n
    type(md5_hash) :: hash
    logical :: strict
    strict = .true.
    if (present(full)) strict = .not.full
    do n = 1, size(this%cc,2)
      if (strict) then
        call hash%update(this%cc(:this%mesh%ncell_onP,n))
        call hash%update(this%cf(:this%mesh%nface_onP,n))
      else
        call hash%update(this%cc(:,n))
        call hash%update(this%cf(:,n))
      end if
    end do
    string = hash%hexdigest()
  end function

end module sd_vector_type
