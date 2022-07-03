#include "f90_assert.fpp"

module material_model_driver

  use material_model_type
  use material_database_type
  use truchas_logging_services
  implicit none
  private

  public :: init_material_model

  type(material_database), target :: matl_db
  type(material_model), public :: matl_model

  !! These are initialized through the PHYSICS namelist. They are for
  !! internal consumption and should not be used externally.
  integer, public :: nmat
  character(80), public :: materials(100)

contains

  subroutine init_material_model

    use material_factory,  only: load_material_database
    use material_namelist, only: params
    use legacy_matl_api, only: nmat_old => nmat

    integer :: stat
    character(:), allocatable :: errmsg

    call load_material_database(matl_db, params, stat, errmsg)
    if (stat /= 0) call TLS_fatal('error initializing material database: ' // errmsg)

    call matl_model%init(materials(1:nmat), matl_db, stat, errmsg)
    if (stat /= 0) call TLS_fatal('error initializing material model: ' // errmsg)

    nmat_old = matl_model%nphase !TODO: eliminate the need for this

  end subroutine init_material_model

end module material_model_driver
