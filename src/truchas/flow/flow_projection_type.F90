!!
!! FLOW_PROJECTION_TYPE
!!
!! This module provides a type that encapsulates the projection of the face
!! velocities onto a solenoidal field
!!
!! Peter Brady <ptb@lanl.gov>
!! Feb 2018
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  The module defines the derived type FLOW_PROJECTION that encapsulates the data
!!  and procedures for setting up and solving the pressure poisson system and
!!  projecting the face velocities onto a solenoidal field.
!!
!!  Interaction with this object occurs through its methods, generally via the flow_driver
!!
!!  READ_PARAMS(PARAMETER_LIST P) - Initializes the parameter list for this object.
!!    This must be call first for the object to behave sensibly.  P is passed on to
!!    a HYPRE_HYBRID instance.
!!
!!  INIT(FLOW_MESH M) - Allocates the data members for this object and initializes the
!!    internal solver
!!
!!  SETUP(FLOW_PROPS FP, REAL(R8) FACE_VELOCITY(:), REAL(R8) DT)
!!    prepares the lhs and rhs of the poisson system. FACE_VELOCITY is an array of outward
!!    face normal velocities in a "ragged" format indexed directly with xcface.
!!    All face data associated with on process cells must be valid
!!
!!  SOLVE(REAL(R8) SOLUTION(:))
!!    solves the system and returns the solution.  Only the on_process portion of
!!    SOLUTION is valid.
!!
!!  CORRECT_FACE_VELOCITY(REAL(R8) DP(:), REAL(R8) old_vel(:), REAL(R8) DT, REAL(R8) new_vel(:))
!!    Uses DP to project the predicted face velocity (old_vel) onto a divergence free
!!    field (new_vel)
!!
!!


module flow_projection_type

  use kinds, only: r8
  use truchas_logging_services
  use truchas_timers
  use flow_mesh_type
  use unstr_mesh_type
  use flow_props_type
  use fischer_guess_type
  use index_partitioning
  use hypre_hybrid_type
  use pcsr_matrix_type
  use parameter_list_type
  implicit none
  private

  public :: flow_projection

  type :: flow_projection
    private
    type(flow_mesh), pointer :: mesh ! unowned reference
    type(fischer_guess) :: fg
    type(hypre_hybrid) :: solver
    type(parameter_list) :: p
    real(r8), allocatable :: rhs(:)
    real(r8), allocatable :: dx_scaled(:,:) ! dxyz/||dxyz||^2
  contains
    procedure :: read_params
    procedure :: init
    procedure :: setup
    procedure :: solve
    procedure :: correct_face_velocity
    procedure :: divergence
  end type flow_projection

contains

  subroutine read_params(this, p)
    class(flow_projection), intent(inout) :: this
    type(parameter_list), intent(inout) :: p

    this%p = p
    call this%fg%read_params(p)
  end subroutine read_params

  subroutine init(this, mesh)
    class(flow_projection), intent(inout) :: this
    type(flow_mesh), pointer, intent(in) :: mesh
    !-
    integer :: j, i, fi, ni, ri, f0
    type(pcsr_graph), pointer :: g
    type(pcsr_matrix), pointer :: A
    type(ip_desc), pointer :: row_ip
    type(unstr_mesh), pointer :: m
    type(parameter_list), pointer :: p

    this%mesh => mesh
    m => mesh%mesh
    p = this%p

    allocate(this%rhs(m%ncell))
    allocate(this%dx_scaled(3, size(m%cface)))

    !! Create a CSR matrix graph for the pressure poisson system.
    allocate(g)
    row_ip => m%cell_ip
    call g%init (row_ip)
    do j = 1, m%ncell_onP
      call g%add_edge(j,j)
      associate (cn => m%cnhbr(m%xcnhbr(j):m%xcnhbr(j+1)-1))
        do i = 1, size(cn)
          if (cn(i) > 0) call g%add_edge(j, i)
        end do
      end associate
    end do
    call g%add_complete

    allocate(A)
    call A%init(g, take_graph=.true.)
    call this%solver%init(A, p)

    !! Create dx_scaled to compute face centered gradients of pressure
    this%dx_scaled = 0.0_r8

    do j = 1, m%ncell
      associate (cn => m%cnhbr(m%xcnhbr(j):m%xcnhbr(j+1)-1))
        f0 = m%xcface(j)-1

        do i = 1, size(cn)
          ri = f0+i ! ''ragged'' index
          fi = m%cface(ri) ! face index
          ni = cn(i) ! neighbor index

          if (ni > 0) then
            this%dx_scaled(:,ri) = this%mesh%cell_centroid(:,ni) - this%mesh%cell_centroid(:,j)
          else
            this%dx_scaled(:,ri) = this%mesh%face_centroid(:,fi) - this%mesh%cell_centroid(:,j)
          end if
          this%dx_scaled(:,ri) = this%dx_scaled(:,ri)/sum(this%dx_scaled(:,ri)**2)
        end do
      end associate
    end do

    call this%fg%init(mesh)

  end subroutine init


  subroutine setup(this, props, vel, dt)
    class(flow_projection), intent(inout) :: this
    type(flow_props), intent(in) :: props
    real(r8), intent(in) :: vel(:), dt
    !-
    type(pcsr_matrix), pointer :: A
    type(unstr_mesh), pointer :: m
    integer :: i, j, k, fi, ni, f0, f1, ri
    real(r8) :: coeff

    A = this%solver%matrix()
    call A%set_all(0.0_r8)

    m => this%mesh%mesh

    do j = 1, m%ncell_onP
      associate (cn => m%cnhbr(m%xcnhbr(j):m%xcnhbr(j+1)-1))
        f0 = m%xcface(j)-1

        do i = 1, size(cn)
          ri = f0+i ! ! ragged index
          fi = m%cface(ri) ! face index
          ni = cn(i) ! neighbor index

          !! FIXME: COME BACK AND HANDLE SOLID FACES, AND SOLID/VOID CELLS
          ! note that normal is already weighted by face area
          if (btest(m%cfpar(j), pos=i)) then
            coeff = dot_product(this%dx_scaled(:,ri), m%normal(:,fi))/props%rho_fc(fi)
          else
            coeff = -dot_product(this%dx_scaled(:,ri), m%normal(:,fi))/props%rho_fc(fi)
          end if

          if (ni > 0) then
            ! DOUBLE CHECK THE SIGN OF THESE COEFFICIENTS
            call A%add_to(j, j, coeff)
            call A%add_to(j, ni, -coeff)
          else
            call A%add_to(j, j, coeff)
          end if
        end do
      end associate
    end do

    call this%solver%setup()

    this%rhs = 0.0_r8
    !! FIXME: need correction for void
    do j = 1, m%ncell_onP
      f0 = m%xcface(j)
      f1 = m%xcface(j+1)-1
      do i = f0, f1
        k = m%cface(i)
        this%rhs(j) = this%rhs(j) + vel(i)*m%area(k)
      end do
      this%rhs(j) = this%rhs(j) / dt
    end do


  end subroutine setup


  subroutine solve(this, solution)
    class(flow_projection), intent(inout) :: this
    real(r8), intent(inout) :: solution(:)
    !-
    integer :: ierr

    call this%fg%guess(this%rhs, solution)
    call this%solver%solve(this%rhs, solution, ierr)
    if (ierr /= 0) call tls_error("projection solve unsuccessful")
    call this%fg%update(this%rhs, solution, this%solver%matrix())

  end subroutine solve


  subroutine correct_face_velocity(this, props, dp, vel_old, dt, vel)
    class(flow_projection), intent(inout) :: this
    type(flow_props) :: props
    real(r8), intent(in) :: dp(:), vel_old(:), dt
    real(r8), intent(out) :: vel(:)
    !-
    type(unstr_mesh), pointer :: m
    integer :: i, j, fi, ni, f0, f1, ri

    m => this%mesh%mesh

    do j = 1, m%ncell_onP
      associate (cn => m%cnhbr(m%xcnhbr(j):m%xcnhbr(j+1)-1))
        f0 = m%xcface(j)-1

        do i = 1, size(cn)
          ri = f0+i ! ragged index
          fi = m%cface(ri) ! face index
          ni = cn(i) ! neighbor index

          !! FIXME: COME BACK AND HANDLE SOLID FACES, AND SOLID/VOID CELLS
          if (btest(m%cfpar(j), pos=i)) then
            vel(ri) = vel_old(ri) + dt*dot_product(this%dx_scaled(:,ri),m%normal(:,fi))/&
                (props%rho_fc(fi)*m%area(fi))
          else
            vel(ri) = vel_old(ri) - dt*dot_product(this%dx_scaled(:,ri),m%normal(:,fi))/&
                (props%rho_fc(fi)*m%area(fi))
          end if
        end do
      end associate
    end do

  end subroutine correct_face_velocity


  subroutine divergence(this, vel, div)
    class(flow_projection), intent(inout) :: this
    real(r8), intent(in) :: vel(:)
    real(r8), intent(out) :: div(:)
    !-
    type(unstr_mesh), pointer :: m
    integer :: j, f0, f1, i, k, out

    div = 0.0_r8
    m => this%mesh%mesh

    !! FOR DEBUGGING PURPOSES
    open(file='divergence', newunit=out)
    do j = 1, m%ncell_onP
      f0 = m%xcface(j)
      f1 = m%xcface(j+1)-1
      do i = f0, f1
        k = m%cface(i)
        div(j) = div(j) + vel(i)*m%area(k)
      end do
      write(out, '(4es13.5)') this%mesh%cell_centroid(:,j), div(j)
    end do
    close(out)

    ! FIXEME: Truchas zeros out divergence in solid/void cells.  Seems like a hack
  end subroutine divergence

end module flow_projection_type
