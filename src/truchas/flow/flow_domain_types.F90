!! In general, the flow mesh involves cells where the velocities and
!! pressure will be updated (termed `regular` cells) and cells where
!! these quantities will not be updated (either `void` or `solid` cells).
!!
!! This module defines three integer parameters which uniquely identify all
!! cells on the flow mesh.  The criteria for assigning a particualar type to a
!! cell requires:
!!   1 - A volume fraction describing the amount of fluid
!!   2 - A volume fraction describing the amount of void
!!   3 - A prescibed cutoff value (defaults to 0.01 in flow_props_type.F90)
!!
!! A cell is `regular_t` if #1 >= cutoff
!! A cell is `void_t` if #1 < cutoff and #2 >= cutoff
!! A cell is `solid_t` if #1 < cutoff and #2 < cutoff
!!
!! The type given to a face is the maxium of the neighboring cell types
!! regular_t/void_t -> void_t
!! regular_t/solid_t -> solid_t
!! void_t/solid_t -> solid_t
!!
!! These criteria are enforced in the flow_props_type
module flow_domain_types
  ! cell/face types for flow algorithm
  integer, parameter :: regular_t = 0
  integer, parameter :: void_t = 1
  integer, parameter :: solid_t = 2
end module flow_domain_types
