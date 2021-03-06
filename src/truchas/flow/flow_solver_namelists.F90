!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module flow_solver_namelists

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parameter_list_type
  implicit none
  private

  public :: read_flow_pressure_solver_namelist, read_flow_viscous_solver_namelist

  real(r8) :: rel_tol, abs_tol, conv_rate_tol, amg_strong_threshold
  integer :: max_ds_iter, max_amg_iter, gmres_krylov_dim, cg_use_two_norm
  integer :: print_level, amg_max_levels, amg_coarsen_type, amg_coarsen_sweeps
  integer :: amg_smoothing_method, amg_smoothing_sweeps, amg_interp_method
  character(128) :: krylov_method, amg_coarsen_method

contains

  subroutine read_flow_pressure_solver_namelist(lun, plist)

    integer, intent(in) :: lun
    type(parameter_list), intent(inout) :: plist

    namelist /flow_pressure_solver/ rel_tol, &
        abs_tol, conv_rate_tol, max_ds_iter, max_amg_iter, &
        krylov_method, gmres_krylov_dim, cg_use_two_norm, &
        print_level, amg_strong_threshold, &
        amg_max_levels, amg_coarsen_method, amg_coarsen_type, &
        amg_smoothing_sweeps, amg_smoothing_method, amg_interp_method

    call read_hypre_namelist(lun, 'FLOW_PRESSURE_SOLVER', read_nml, plist)

  contains

    subroutine read_nml(lun, ios, iom)
      integer, intent(in) :: lun
      integer, intent(out) :: ios
      character(*), intent(inout) :: iom
      read(lun,nml=flow_pressure_solver,iostat=ios,iomsg=iom)
    end subroutine

  end subroutine read_flow_pressure_solver_namelist


  subroutine read_flow_viscous_solver_namelist(lun, plist)

    integer, intent(in) :: lun
    type(parameter_list), intent(inout) :: plist

    namelist /flow_viscous_solver/ rel_tol, &
        abs_tol, conv_rate_tol, max_ds_iter, max_amg_iter, &
        krylov_method, gmres_krylov_dim, cg_use_two_norm, &
        print_level, amg_strong_threshold, &
        amg_max_levels, amg_coarsen_method, amg_coarsen_type, &
        amg_smoothing_sweeps, amg_smoothing_method, amg_interp_method

    call read_hypre_namelist(lun, 'FLOW_VISCOUS_SOLVER', read_nml, plist)

  contains

    subroutine read_nml(lun, ios, iom)
      integer, intent(in) :: lun
      integer, intent(out) :: ios
      character(*), intent(inout) :: iom
      read(lun,nml=flow_viscous_solver,iostat=ios,iomsg=iom)
    end subroutine

  end subroutine read_flow_viscous_solver_namelist

  subroutine read_hypre_namelist(lun, name, read_nml, p)

    use string_utilities, only: i_to_c
    use parallel_communication, only: is_IOP, broadcast
    use flow_input_utils
    use truchas_logging_services

    integer, intent(in) :: lun
    character(*), intent(in) :: name
    type(parameter_list), intent(inout) :: p

    type(parameter_list), pointer :: pp
    integer :: ios
    logical :: found
    character(128) :: iom

    interface
      subroutine read_nml(lun, ios, iom)
        integer, intent(in) :: lun
        integer, intent(out) :: ios
        character(*), intent(inout) :: iom
      end subroutine
    end interface

    max_ds_iter = NULL_I
    max_amg_iter = NULL_I
    gmres_krylov_dim = NULL_I
    cg_use_two_norm = NULL_I
    print_level = NULL_I
    amg_max_levels = NULL_I
    amg_coarsen_type = NULL_I
    amg_coarsen_sweeps = NULL_I
    amg_smoothing_method = NULL_I
    amg_smoothing_sweeps = NULL_I
    amg_interp_method = NULL_I
    rel_tol = NULL_R
    abs_tol = NULL_R
    conv_rate_tol = NULL_R
    amg_strong_threshold = NULL_R
    krylov_method = NULL_C
    amg_coarsen_method = NULL_C

    if (is_IOP) then
      rewind lun
      call seek_to_namelist(lun, name, found, iostat=ios)
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))

    call broadcast(found)
    if (found) then
      call TLS_info('Reading ' // name // ' namelist ...')
      !! Read the namelist.
      if (is_IOP) then
        call read_nml(lun, ios, iom)
      end if
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading ' // name // ' namelist: ' // trim(iom))

      call broadcast(max_ds_iter)
      call broadcast(max_amg_iter)
      call broadcast(gmres_krylov_dim)
      call broadcast(cg_use_two_norm)
      call broadcast(print_level)
      call broadcast(amg_max_levels)
      call broadcast(amg_coarsen_type)
      call broadcast(amg_coarsen_sweeps)
      call broadcast(amg_smoothing_method)
      call broadcast(amg_interp_method)
      call broadcast(rel_tol)
      call broadcast(abs_tol)
      call broadcast(conv_rate_tol)
      call broadcast(amg_strong_threshold)
      call broadcast(krylov_method)
      call broadcast(amg_coarsen_method)

      pp => p%sublist("solver")
      call plist_set_if(pp, "rel-tol", rel_tol)
      call plist_set_if(pp, 'abs-tol', abs_tol)
      call plist_set_if(pp, 'conv-rate-tol', conv_rate_tol)
      call plist_set_if(pp, 'max-ds-iter', max_ds_iter)
      call plist_set_if(pp, 'max-amg-iter', max_amg_iter)
      call plist_set_if(pp, 'krylov-method', krylov_method)
      call plist_set_if(pp, 'gmres-krylov-dim', gmres_krylov_dim)
      call plist_set_if(pp, 'cg-use-two-norm', cg_use_two_norm /= 0)
      call plist_set_if(pp, 'print-level', print_level)
      call plist_set_if(pp, 'amg-strong-threshold', amg_strong_threshold)
      call plist_set_if(pp, 'amg-max-levels', amg_max_levels)
      call plist_set_if(pp, 'amg-coarsen-method', amg_coarsen_method)
      call plist_set_if(pp, 'amg-coarsen-type', amg_coarsen_type)
      call plist_set_if(pp, 'amg-smoothing-sweeps', amg_smoothing_sweeps)
      call plist_set_if(pp, 'amg-smoothing-method', amg_smoothing_method)
      call plist_set_if(pp, 'amg-interp-method', amg_interp_method)
    end if

  end subroutine read_hypre_namelist

end module flow_solver_namelists
