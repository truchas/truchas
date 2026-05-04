!!
!! VIZ_FIELD_REGISTRY_TYPE
!!
!! Registry of visualization fields offered by active providers.
!!
!! Providers register flat, user-visible field names here after they are known
!! to be active. The registry enforces a single global namespace, maps each
!! field name to the provider that owns it, resolves a stream's requested
!! names, and groups resolved names by provider so the manager can create
!! provider states with provider-local field selections.
!!
!! The provider id stored here is an opaque active-provider-array index supplied
!! by viz_manager_type during setup. It is only used to route selected fields
!! back to the provider that registered them.
!!
!! Missing requested fields are handled outside this type by higher-level
!! stream configuration code; duplicate registered names are fatal because
!! the current user input has no namespace qualification.
!!

#include "f90_assert.fpp"

module viz_field_registry_type

  implicit none
  private

  type :: viz_field_info
    character(:), allocatable :: name
    ! Active provider array index; not a stable provider identity.
    integer :: provider_id = 0
  end type

  type, public :: viz_provider_field_selection
    ! Setup-time grouped selection routed by active provider array index.
    integer :: provider_id = 0
    character(:), allocatable :: field_names(:)
  end type

  type, public :: viz_field_registry
    private
    type(viz_field_info), allocatable :: fields(:)
  contains
    procedure :: register_field
    procedure :: resolve_selection
    procedure :: group_fields
    procedure, private :: resolve
  end type

contains

  subroutine register_field(this, name, provider_id)

    class(viz_field_registry), intent(inout) :: this
    character(*), intent(in) :: name
    integer, intent(in) :: provider_id

    type(viz_field_info), allocatable :: tmp(:)
    integer :: n

    call this%resolve(name, n)
    INSIST(n == 0)

    if (allocated(this%fields)) then
      n = size(this%fields) + 1
      allocate(tmp(n))
      tmp(:n-1) = this%fields
      call move_alloc(tmp, this%fields)
    else
      n = 1
      allocate(this%fields(n))
    end if

    this%fields(n)%name = trim(name)
    this%fields(n)%provider_id = provider_id

  end subroutine register_field

  !! Return the field id for the given field name, or 0 if not found.
  subroutine resolve(this, name, field_id)

    class(viz_field_registry), intent(in) :: this
    character(*), intent(in) :: name
    integer, intent(out) :: field_id

    integer :: i

    field_id = 0
    if (allocated(this%fields)) then
      do i = 1, size(this%fields)
        if (this%fields(i)%name == trim(name)) then
          field_id = i
          return
        end if
      end do
    end if

  end subroutine resolve

  !! Resolve requested field names into registry ids and collect names that
  !! are unavailable in the current active-provider set.

  subroutine resolve_selection(this, names, field_ids, missing_names)

    class(viz_field_registry), intent(in) :: this
    character(*), intent(in) :: names(:)
    integer, allocatable, intent(out) :: field_ids(:)
    character(:), allocatable, intent(out) :: missing_names(:)

    integer, allocatable :: ids(:)
    integer :: i, nfield, nmissing, field_id
    character(max_field_name_len(names)), allocatable :: missing(:)

    allocate(ids(size(names)))
    allocate(missing(size(names)))
    nfield = 0
    nmissing = 0

    do i = 1, size(names)
      call this%resolve(names(i), field_id)
      if (field_id /= 0) then
        if (any(ids(:nfield) == field_id)) cycle
        nfield = nfield + 1
        ids(nfield) = field_id
      else
        nmissing = nmissing + 1
        missing(nmissing) = trim(names(i))
      end if
    end do

    allocate(field_ids(nfield))
    field_ids = ids(:nfield)
    allocate(character(len=max_field_name_len(names)) :: missing_names(nmissing))
    missing_names = missing(:nmissing)

  end subroutine resolve_selection

  !! Group resolved field ids by owning provider so stream setup can create
  !! one provider state per selected provider.

  subroutine group_fields(this, field_ids, groups)

    class(viz_field_registry), intent(in) :: this
    integer, intent(in) :: field_ids(:)
    type(viz_provider_field_selection), allocatable, intent(out) :: groups(:)

    integer, allocatable :: provider_ids(:), counts(:)
    integer :: i, n, provider_id, group_id

    allocate(provider_ids(size(field_ids)))
    allocate(counts(size(field_ids)))
    provider_ids = 0
    counts = 0
    n = 0

    do i = 1, size(field_ids)
      provider_id = this%fields(field_ids(i))%provider_id
      group_id = findloc(provider_ids(:n), provider_id, dim=1)
      if (group_id == 0) then
        n = n + 1
        provider_ids(n) = provider_id
        group_id = n
      end if
      counts(group_id) = counts(group_id) + 1
    end do

    allocate(groups(n))
    do i = 1, n
      groups(i)%provider_id = provider_ids(i)
      allocate(character(max_registered_field_name_len(this)) :: &
          groups(i)%field_names(counts(i)))
      counts(i) = 0
    end do

    do i = 1, size(field_ids)
      provider_id = this%fields(field_ids(i))%provider_id
      group_id = findloc(provider_ids(:n), provider_id, dim=1)
      counts(group_id) = counts(group_id) + 1
      groups(group_id)%field_names(counts(group_id)) = this%fields(field_ids(i))%name
    end do

  end subroutine group_fields

  pure integer function max_field_name_len(names) result(len)
    character(*), intent(in) :: names(:)
    integer :: i
    len = 0
    do i = 1, size(names)
      len = max(len, len_trim(names(i)))
    end do
  end function

  pure integer function max_registered_field_name_len(this) result(len)
    class(viz_field_registry), intent(in) :: this
    integer :: i
    len = 0
    if (.not.allocated(this%fields)) return
    do i = 1, size(this%fields)
      len = max(len, len_trim(this%fields(i)%name))
    end do
  end function

end module viz_field_registry_type
