module flow_surface_tension_bc_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_mesh_type
  use bndry_vfunc_class
  private

  type, extends(bndry_vfunc), public :: surface_tension_bc
    private
    real(r8) :: dsig_dT

    type(unstr_mesh), pointer :: mesh => null()
    real(r8), pointer :: vof(:) => null(), temperature_fc(:) => null()
  contains
    procedure :: init
    procedure :: compute
  end type surface_tension_bc

contains

  subroutine init(this, params, mesh, vof, temperature_fc)

    use parameter_list_type
    use truchas_logging_services
    use bndry_face_group_builder_type
    use string_utilities, only: raise_case

    class(surface_tension_bc), intent(out) :: this
    type(parameter_list), intent(inout) :: params
    type(unstr_mesh), target, intent(in) :: mesh
    real(r8), target, intent(in) :: vof(:), temperature_fc(:)

    type(bndry_face_group_builder) :: builder
    type(parameter_list_iterator) :: piter
    type(parameter_list), pointer :: bc_params
    integer :: stat, ngroup
    integer, allocatable :: xgroup(:), setids(:)
    character(:), allocatable :: bc_type, errmsg
    logical :: found

    this%mesh => mesh
    this%vof => vof
    this%temperature_fc => temperature_fc

    ! initialize index list and get dsig_dT
    ! finds all face indices associated with face sets given
    ! by boundary conditions with "condition = 'surface tension'"
    ! NOTE: currently only one surface tension BC namelist is supported
    call builder%init(this%mesh)

    found = .false.
    piter = parameter_list_iterator(params, sublists_only=.true.)
    do while (.not.piter%at_end())
      bc_params => piter%sublist()
      call bc_params%get('condition', bc_type)
      if (raise_case(bc_type) == raise_case('surface tension')) then
        if (found) &
            call TLS_fatal('error: Found more than one surface tension BC namelist. ' // &
            'Currently only one supported.')
        call bc_params%get('data', this%dsig_dT)
        call bc_params%get('face sets', setids)
        call builder%add_face_group(setids, stat, errmsg)
        if (stat /= 0) call TLS_fatal('error generating boundary condition: ' // errmsg)
        found = .true.
      end if
      call piter%next()
    end do

    call builder%get_face_groups(ngroup, xgroup, this%index)

    allocate(this%value(3,size(this%index)))

  end subroutine init

  subroutine compute(this, t)

    class(surface_tension_bc), intent(inout) :: this
    real(r8), intent(in) :: t

    integer :: i, j, f, fi, fn
    real(r8) :: gradT(3), maxdist, dist(3)

    associate (faces => this%index, value => this%value)
      do i = 1, size(faces)
        fi = faces(i)
        j = this%mesh%fcell(1,fi)
        if (j > this%mesh%ncell_onP) cycle

        ! calculate the temperature gradient with green-gauss
        associate (facen => this%mesh%cface(this%mesh%xcface(j):this%mesh%xcface(j+1)-1))
          gradT = 0
          do f = 1,size(facen)
            fn = facen(f)

            if (btest(this%mesh%cfpar(j),pos=f)) then ! true if normal points inward
              gradT = gradT - this%temperature_fc(fn) * this%mesh%normal(:,fn)
            else
              gradT = gradT + this%temperature_fc(fn) * this%mesh%normal(:,fn)
            end if
          end do
          gradT = gradT / this%mesh%volume(j)

          ! get the component of the temperature gradient tangent to the boundary
          gradT = gradT &
              - dot_product(gradT, this%mesh%normal(:,fi)) * this%mesh%normal(:,fi) &
              / this%mesh%area(fi)**2

          ! the surface tension term here is a volume-averaged integral over the
          ! computational cell, but the term is only active on the boundary face,
          ! giving a face_area / volume component.
          value(:,i) = this%vof(j) * this%dsig_dT * gradT * this%mesh%area(fi)
        end associate
      end do
    end associate

  end subroutine compute

end module flow_surface_tension_bc_type
