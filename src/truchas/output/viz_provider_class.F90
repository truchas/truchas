!!
!! VIZ_PROVIDER_CLASS
!!
!! Abstract setup-time descriptor/factory for package participation in VTKHDF
!! output. A provider reports whether its package is active for the current
!! simulation, registers the selectable field names owned by that package,
!! and creates a stream-local viz_provider_state for either all fields or a
!! selected subset.
!!
!! The provider_id passed to register_fields() is the provider's current index
!! in the active provider array. Providers must pass it unchanged to
!! viz_field_registry%register_field(); it is not a persistent package id or
!! user-facing namespace.
!!
!! Providers are intentionally lightweight and should not hold open-file or
!! per-stream execution state; that belongs to the viz_provider_state objects
!! they create.
!!

module viz_provider_class

  use viz_field_registry_type
  use viz_provider_state_class
  implicit none
  private

  type, abstract, public :: viz_provider
  contains
    procedure(is_active_ifc), deferred :: is_active
    procedure(register_fields_ifc), deferred :: register_fields
    procedure(create_state_ifc), deferred :: create_state
  end type

  type, public :: viz_provider_box
    class(viz_provider), allocatable :: p
  end type

  abstract interface
    logical function is_active_ifc(this)
      import :: viz_provider
      class(viz_provider), intent(in) :: this
    end function

    subroutine register_fields_ifc(this, registry, provider_id)
      import :: viz_field_registry, viz_provider
      class(viz_provider), intent(in) :: this
      type(viz_field_registry), intent(inout) :: registry
      integer, intent(in) :: provider_id
    end subroutine

    subroutine create_state_ifc(this, nblock, state, field_names)
      import :: viz_provider, viz_provider_state
      class(viz_provider), intent(in) :: this
      integer, intent(in) :: nblock
      class(viz_provider_state), allocatable, intent(out) :: state
      character(*), intent(in), optional :: field_names(:)
    end subroutine
  end interface

end module viz_provider_class
