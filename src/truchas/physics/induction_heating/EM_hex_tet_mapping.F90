!!
!! EM_HEX_TET_MAPPING
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! 22 Jun 2005
!!
!! This module provides procedures for performing the mapping of cell-based
!! fields between the EM tet mesh and the main Truchas hex mesh.  This is a
!! thin layer over GRID_MAPPING_MODULE where all the real work is performed.
!!

#include "f90_assert.fpp"

module EM_hex_tet_mapping

  use kinds, only : rk => r8
  use grid_mapping_module

  implicit none
  private

  public :: get_grid_mapping_data, map_hex_field, map_tet_field

  !! Objects from GRID_MAPPING_MODULE required by a user of EM_HEX_TET_MAPPING.
  public :: gm_mesh, destroy_gm_mesh, grid_int_vols, destroy_grid_int_vols

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! GET_GRID_MAPPING_DATA
 !!
 !! The subroutine returns the grid mapping data GMD required to map data
 !! between the main Truchas mesh HEXMESH and the EM mesh TETMESH.  It attempts
 !! to read this data from the file specified in the ALTMESH namelist.  If that
 !! isn't successful it has the mapping data computed and written to the file
 !! for use in future runs.  This is a serial subroutine and should only be
 !! called from the IO processor.
 !!

  subroutine get_grid_mapping_data (hexmesh, tetmesh, gmd)

    use altmesh_input, only: grid_transfer_file
    use parameter_module, only: string_len
    use truchas_logging_services
#ifdef SUPPORTS_NEWUNIT
    use truchas_env, only: output_dir
#else
    use truchas_env, only: output_dir, new_unit
#endif

    type(gm_mesh), intent(in) :: hexmesh, tetmesh
    type(grid_int_vols), intent(out) :: gmd

    integer :: lun, status
    logical :: file_exists, compute
    character(len=string_len) :: mapfile

    call TLS_info ('  Initializing the hex-tet grid mapping data ...')

    compute = .true.

    inquire(file=grid_transfer_file,exist=file_exists)
    if (file_exists) then
      call TLS_info ('   Reading the hex-tet grid mapping data from ' // trim(grid_transfer_file))
#ifdef SUPPORTS_NEWUNIT
      open(newunit=lun,file=grid_transfer_file,action='read',position='rewind',form='unformatted',iostat=status)
#else
      call new_unit (lun)
      open(lun,file=grid_transfer_file,action='read',position='rewind',form='unformatted',iostat=status)
#endif
      if (status == 0) then
        call read_int_volumes (gmd, lun, status)
        close(lun)
        if (status == 0) then
          compute = .not. right_int_volumes(hexmesh, tetmesh, gmd)
          if (compute) then
            call destroy_grid_int_vols (gmd)
            call TLS_info ('    Grid mapping data is not usable.')
          else
            call TLS_info ('    Using the grid mapping data.')
          endif
        else
          call destroy_grid_int_vols (gmd)
          call TLS_info ('    Error reading the file!')
        end if
      else
        call TLS_info ('    Unable to open the file for unformatted reading!')
      end if
    end if

    if (compute) then
      mapfile = 'altmesh_mapping_data.bin'
      call TLS_info ('   Computing the hex-tet grid mapping data.')
      call compute_int_volumes (hexmesh, tetmesh, gmd)
      call TLS_info ('   Writing the hex-tet grid mapping data to ' // trim(mapfile))
      mapfile = trim(output_dir) // trim(mapfile)
#ifdef SUPPORTS_NEWUNIT
      open(newunit=lun,file=trim(mapfile),action='write',position='rewind',form='unformatted',iostat=status)
#else
      call new_unit (lun)
      open(lun,file=trim(mapfile),action='write',position='rewind',form='unformatted',iostat=status)
#endif
      if (status == 0 ) then
        call write_int_volumes (gmd, lun)
        close(lun)
      else
        call TLS_info ('    Unable to open the file for unformatted writing!  Continuing anyway.')
      end if
    end if

    call TLS_info ('   Hex-tet grid mapping data initialized.')

  end subroutine get_grid_mapping_data

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! MAP_HEX_FIELD
 !!
 !! Maps the distributed, cell-based array HEX_FIELD on the main hex mesh to
 !! the distributed, cell-based array TET_FIELD on the EM tet mesh.  The
 !! mapping is constant-preserving and strict, and DEFVAL specifies the value
 !! to use for tets that do not intersect the hex mesh (see the documentation
 !! for GRID_MAPPING_MODULE).  It is assumed that the grid mapping data GMD was
 !! calculated with the hex mesh as MESH_A and the tet mesh as MESH_B (on the
 !! IO processor only.)  This is a parallel subroutine.
 !!

  subroutine map_hex_field (gmd, hex_field, tet_field, defval)
    type(grid_int_vols), intent(in) :: gmd
    real(rk), intent(in) :: defval
    real(rk), intent(in) :: hex_field(:)
    real(rk), intent(out) :: tet_field(:)
    call map_cell_field_aux (gmd, 1, defval, hex_field, tet_field)
  end subroutine map_hex_field

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! MAP_TET_FIELD
 !!
 !! Maps the distributed, cell-based array TET_FIELD on the EM tet mesh to
 !! the distributed, cell-based array HEX_FIELD on the main hex mesh.  The
 !! mapping is conservative and strict, and DEFVAL specifies the value to use
 !! for hexes that do not intersect the tet mesh (see the documentation for
 !! GRID_MAPPING_MODULE).  It is assumed that the grid mapping data GMD was
 !! calculated with the hex mesh as MESH_A and the tet mesh as MESH_B (on the
 !! IO processor only.)  This is a parallel subroutine.
 !!

  subroutine map_tet_field (gmd, tet_field, hex_field, defval)
    type(grid_int_vols), intent(in) :: gmd
    real(rk), intent(in) :: defval
    real(rk), intent(in) :: tet_field(:)
    real(rk), intent(out) :: hex_field(:)
    call map_cell_field_aux (gmd, 2, defval, tet_field, hex_field)
  end subroutine map_tet_field

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! MAP_CELL_FIELD_AUX
 !!
 !! This auxillary subroutine is a parallel wrapper around the serial
 !! MAP_CELL_FIELD procedure from the GRID_MAPPING_MODULE module.  The
 !! distributed source field is collated, mapped in serial, and the collated
 !! result distributed to the destination field.  Two forms of the mapping
 !! are possible: FORM=1  is the constant-preserving hex-to-tet mapping,
 !! and FORM=2 is the conservative tet-to-hex mapping.  It is assumed that
 !! the grid mapping data GMD was generated with the hex mesh as MESH_A and
 !! the tet mesh as MESH_B
 !!

  subroutine map_cell_field_aux (gmd, form, defval, src, dest)

    use parallel_communication

    type(grid_int_vols), intent(in) :: gmd
    integer,  intent(in) :: form
    real(rk), intent(in) :: defval
    real(rk), intent(in) :: src(:)
    real(rk), intent(out) :: dest(:)

    real(rk), pointer :: col_src(:), col_dest(:)

    call allocate_collated_array (col_src,  global_sum(size(src)))
    call allocate_collated_array (col_dest, global_sum(size(dest)))

    call collate (col_src, src)

    if (is_IOP) then
      select case (form)
      case (1)  ! hex-to-tet map; we want this to preserve constants.
        call map_cell_field (col_src, col_dest, gmd, preserve_constants=.true., strict=.true., defval=defval)
      case (2)  ! tet-to-hex map; we want this to conserve.
        call map_cell_field (col_src, col_dest, gmd, reverse_order=.true.,exactly_conservative=.true.,strict=.true.,defval=defval)
      case default
        ASSERT( .false. )
      end select
    end if

    call distribute (dest, col_dest)

    deallocate(col_src, col_dest)

  end subroutine map_cell_field_aux

end module EM_hex_tet_mapping
