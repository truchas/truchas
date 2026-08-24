!!
!! HT_2D_MODEL_TYPE
!!
!! This module defines the 2D thermal-transport model. It owns the mimetic
!! finite difference discretization, material-dependent thermal properties,
!! boundary conditions, and heat sources, and supplies the cell-and-face
!! residual and  vector initialization required by the implicit thermal
!! integrator.
!!
!! David Neill-Asanza <davidhneill@gmail.com>, July 2020
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module ht_2d_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_2d_mesh_type
  use mfd_2d_disc_type
  use bndry_func1_class
  use bndry_func2_class
  use scalar_mesh_multifunc_type
  use new_mesh_func_class
  use cell_matl_prop_func_type
  use ht_2d_tofh_type
  use ht_2d_vector_type
  use material_composition_type
  use parallel_communication
  use parameter_list_type
  use simulation_environment_type
  implicit none
  private

  type, public :: ht_2d_model
    type(unstr_2d_mesh), pointer :: mesh => null() ! unowned reference
    type(mfd_2d_disc) :: disc
    !! Equation parameters
    class(new_mesh_func), allocatable :: conductivity
    class(new_mesh_func), allocatable :: H_of_T
    type(ht_2d_tofh) :: T_of_H            ! inverse enthalpy-temperature relation
    type(scalar_mesh_multifunc), allocatable :: source  ! external heat source
    real(r8), allocatable :: ext_enthalpy_rate(:)       ! cell-integrated external rate
    !! Boundary condition data
    class(bndry_func1), allocatable :: bc_dir   ! Dirichlet
    class(bndry_func1), allocatable :: bc_flux  ! Simple flux
    class(bndry_func2), allocatable :: bc_htc   ! External heat transfer coefficient
    class(bndry_func2), allocatable :: bc_rad   ! Simple radiation
  contains
    procedure :: init
    procedure :: init_vector
    procedure :: set_ext_enthalpy_rate
    procedure :: residual
  end type

contains

  subroutine init(this, env, mesh, matl_model, composition, params, stat, errmsg)

    use material_model_type
    use material_utilities

    class(ht_2d_model), intent(out), target :: this
    type(simulation_environment), intent(in) :: env
    type(unstr_2d_mesh), intent(in), target :: mesh
    type(material_model), intent(in) :: matl_model
    type(material_composition), pointer, intent(in) :: composition
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    character(:), allocatable :: context
    type(parameter_list), pointer :: plist
    real(r8) :: sigma, abszero

    context = 'processing ' // params%path() // ': '

    call this%disc%init(mesh)
    this%mesh => mesh
    allocate(this%ext_enthalpy_rate(mesh%ncell_onP), source=0.0_r8)

    call params%get('stefan-boltzmann', sigma, stat, errmsg, default=5.67e-8_r8)
    if (stat /= 0) return
    call params%get('absolute-zero', abszero, stat, errmsg, default=0.0_r8)
    if (stat /= 0) return

    !! Enthalpy density.
    call required_property_check(matl_model, 'enthalpy', stat, errmsg)
    if (stat /= 0) return
    allocate(cell_matl_prop_func :: this%H_of_T)
    select type (H_of_T => this%H_of_T)
    type is (cell_matl_prop_func)
      call H_of_T%init(matl_model, composition, 'enthalpy', stat, errmsg)
    end select
    if (global_any(stat /= 0)) then
      stat = -1
      errmsg = context // 'unexpected error defining H_of_T: ' // errmsg
      return
    end if

    !! Inverse of enthalpy-temperature relation.
    block
      real(r8) :: tol, step
      integer :: ntry
      plist => params%sublist('temperature-from-enthalpy')
      !! Absolute temperature root-finding tolerance, >= 0.
      call plist%get('temperature-tolerance', tol, stat, errmsg, default=0.0d0)
      if (stat /= 0) return
      !! Initial temperature increment for expanding a failed bracket, > 0.
      call plist%get('bracket-initial-step', step, stat, errmsg, default=1.0d-3)
      if (stat /= 0) return
      !! Maximum number of attempts to expand a failed bracket, >= 0.
      call plist%get('bracket-max-tries', ntry, stat, errmsg, default=50)
      if (stat /= 0) return
      call this%T_of_H%init(this%H_of_T, eps=tol, max_try=ntry, delta=step)
    end block

    !! Thermal conductivity.
    call required_property_check(matl_model, 'conductivity', stat, errmsg)
    if (stat /= 0) return
    allocate(cell_matl_prop_func :: this%conductivity)
    select type (conductivity => this%conductivity)
    type is (cell_matl_prop_func)
      call conductivity%init(matl_model, composition, 'conductivity', stat, errmsg)
    end select
    if (global_any(stat /= 0)) then
      stat = -1
      errmsg = context // 'unexpected error defining conductivity: ' // errmsg
      return
    end if

    !! Defines the boundary condition components
    if (params%is_sublist('bc')) then
      plist => params%sublist('bc')
      call init_bc(this, env, plist, sigma, abszero, stat, errmsg)
      if (stat /= 0) return
    else
      stat = 1
      errmsg = context // 'missing "bc" sublist parameter'
      return
    end if

    !! Defines the heat source
    plist => params%sublist('source')
    call init_source(this, env, plist, stat, errmsg)
    if (stat /= 0) return

  end subroutine init


  subroutine init_vector(this, vec)
    class(ht_2d_model), intent(in) :: this
    type(ht_2d_vector), intent(out) :: vec
    call vec%init(this%mesh)
  end subroutine


  !! Set a cell-integrated external enthalpy rate. This is normally supplied
  !! by a coupled physics model, such as explicit enthalpy advection.
  subroutine set_ext_enthalpy_rate(this, enthalpy_rate)
    class(ht_2d_model), intent(inout) :: this
    real(r8), intent(in) :: enthalpy_rate(:)

    ASSERT(size(enthalpy_rate) == this%mesh%ncell_onP)
    this%ext_enthalpy_rate = enthalpy_rate
  end subroutine

  subroutine init_bc(model, env, params, sigma, abszero, stat, errmsg)

    use bitfield_type
    use thermal_bc_factory_type
    use string_utilities, only: i_to_c

    class(ht_2d_model), intent(inout), target :: model
    type(simulation_environment), intent(in) :: env
    type(parameter_list), intent(inout), target :: params
    real(r8), intent(in) :: sigma, abszero
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(thermal_bc_factory) :: bc_fac
    type(bitfield) :: bitmask
    character(160) :: string
    logical, allocatable :: mask(:)
    integer, allocatable :: setids(:)
    integer :: j, n

    allocate(mask(model%mesh%nface), source=.false.)

    call bc_fac%init(model%mesh, sigma, abszero, params)

    !! Define the simple flux boundary conditions.
    call bc_fac%alloc_flux_bc(model%bc_flux, env, stat, errmsg)
    if (stat /= 0) return
    if (allocated(model%bc_flux)) then
      mask(model%bc_flux%index) = .true. ! mark the simple flux faces
    end if

    !! Define external heat-transfer-coefficient boundary conditions. HTC
    !! conditions are flux conditions and may be superimposed with simple
    !! flux conditions, but not with Dirichlet conditions.
    call bc_fac%alloc_htc_bc(model%bc_htc, env, stat, errmsg)
    if (stat /= 0) return
    if (allocated(model%bc_htc)) then
      mask(model%bc_htc%index) = .true. ! mark the HTC faces
    end if

    !! Define simple radiation boundary conditions. Radiation is a flux
    !! condition and may be superimposed with other flux conditions, but not
    !! with Dirichlet conditions.
    call bc_fac%alloc_rad_bc(model%bc_rad, env, stat, errmsg)
    if (stat /= 0) return
    if (allocated(model%bc_rad)) then
      mask(model%bc_rad%index) = .true. ! mark the radiation faces
    end if

    !! Define the Dirichlet boundary conditions.
    call bc_fac%alloc_dir_bc(model%bc_dir, env, stat, errmsg)
    if (stat /= 0) return
    if (allocated(model%bc_dir)) then
      if (global_any(mask(model%bc_dir%index))) then
        stat = -1
        errmsg = 'temperature dirichlet boundary condition overlaps with other conditions'
        return
      end if
      mask(model%bc_dir%index) = .true. ! mark the dirichlet faces
    end if

    !! Finally verify that a condition has been applied to every boundary face.
    mask = mask .neqv. btest(model%mesh%face_set_mask,0)
    if (global_any(mask)) then
      call model%mesh%get_face_set_ids(pack([(j,j=1,model%mesh%nface)], mask), setids)
      if (size(setids) == 0) then
        string = '(none)'
      else
        write(string,'(i0,*(:,", ",i0))') setids
      end if
      errmsg = 'incomplete temperature boundary specification;' // &
          ' remaining boundary faces belong to face sets ' // trim(string)
      bitmask = ibset(ZERO_BITFIELD, 0)
      mask = mask .and. (model%mesh%face_set_mask == bitmask)
      n = global_count(mask(:model%mesh%nface_onP))
      if (n > 0) errmsg = errmsg // '; ' // i_to_c(n) // &
          ' boundary faces belong to no face set'
      stat = -1
      return
    end if

  end subroutine init_bc

  subroutine init_source(model, env, params, stat, errmsg)

    use ht_2d_source_factory_type

    class(ht_2d_model), intent(inout), target :: model
    type(simulation_environment), intent(in) :: env
    type(parameter_list), intent(inout), target :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(ht_2d_source_factory) :: src_fac

    call src_fac%init(model%mesh, params)

    !! Allocated function-based source
    call src_fac%alloc_source_funcs(env, model%source, stat, errmsg)
    if (stat /= 0) return

  end subroutine init_source


  subroutine residual(this, t, u, udot, r)

    class(ht_2d_model), intent(inout) :: this
    real(r8), intent(in) :: t
    type(ht_2d_vector), intent(inout) :: u, udot
    type(ht_2d_vector), intent(inout) :: r

    real(r8), allocatable :: Tdir(:)
    real(r8) :: cval(this%mesh%ncell)

    !! The vector integrator guarantees only on-process entries at callback
    !! boundaries.  Gather the state needed by the MFD residual explicitly.
    call this%mesh%cell_imap%gather_offp(u%tc)
    call this%mesh%face_imap%gather_offp(u%tf)

    call impose_dirichlet
    call compute_thermal_residual
    call restore_dirichlet

  contains

    subroutine impose_dirichlet
      if (allocated(this%bc_dir)) then
        call this%bc_dir%compute(t)
        allocate(Tdir(size(this%bc_dir%index)))
        associate (index => this%bc_dir%index, value => this%bc_dir%value)
          Tdir = u%tf(index)
          u%tf(index) = value
        end associate
      end if
    end subroutine

    subroutine restore_dirichlet
      if (allocated(this%bc_dir)) then
        u%tf(this%bc_dir%index) = Tdir
        deallocate(Tdir)
      end if
    end subroutine

    subroutine compute_thermal_residual

      integer :: ncell_onP

      ncell_onP = this%mesh%ncell_onP

      !! Residual of the algebraic enthalpy-temperature relation.
      call this%H_of_T%compute_value(u%tc, cval(:ncell_onP))
      r%hc(:ncell_onP) = u%hc(:ncell_onP) - cval(:ncell_onP)

      call this%conductivity%compute_value(u%tc, cval(:ncell_onP))
      call this%mesh%cell_imap%gather_offp(cval)
      call this%disc%apply_diff(cval, u%tc, u%tf, r%tc, r%tf)
      r%tc(:ncell_onP) = r%tc(:ncell_onP) + this%mesh%volume(:ncell_onP)*udot%hc(:ncell_onP)

      if (allocated(this%ext_enthalpy_rate)) then
        r%tc(:ncell_onP) = r%tc(:ncell_onP) - this%ext_enthalpy_rate
      end if

      if (allocated(this%source)) then
        call this%source%compute(t, u%tc)
        r%tc(:ncell_onP) = r%tc(:ncell_onP) &
            - this%mesh%volume(:ncell_onP)*this%source%value(:ncell_onP)
      end if

      if (allocated(this%bc_dir)) r%tf(this%bc_dir%index) = Tdir - this%bc_dir%value

      if (allocated(this%bc_flux)) then
        call this%bc_flux%compute(t)
        associate (index => this%bc_flux%index, value => this%bc_flux%value)
          r%tf(index) = r%tf(index) + this%mesh%area(index)*value
        end associate
      end if

      if (allocated(this%bc_htc)) then
        call this%bc_htc%compute(t, u%tf)
        r%tf(this%bc_htc%index) = r%tf(this%bc_htc%index) + &
            this%bc_htc%value
      end if

      if (allocated(this%bc_rad)) then
        call this%bc_rad%compute(t, u%tf)
        r%tf(this%bc_rad%index) = r%tf(this%bc_rad%index) + &
            this%bc_rad%value
      end if

    end subroutine compute_thermal_residual

  end subroutine residual

end module ht_2d_model_type
