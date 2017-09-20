
#include "f90_assert.fpp"

module flow_model

  use kinds, only: r8
  use unstr_mesh_type
  use parameter_list_type
  use truchas_logging_services
  implicit none
  private

  type, public :: flow_model_t
    private
    type(unstr_mesh), pointer :: mesh => null()  ! reference only -- do not own
 end type flow_model_t

contains

end module flow_model
