!!
!! This module provides a type to calculate the thermomechanical
!! response of all relevant solid materials present in a given problem.
!!
!! References:
!!
!! - Truchas Physics & Algorithms handbook.
!!
!! - C. Bailey and M. Cross. A finite volume procedure to solve elastic solid
!! mechanics problems in three dimensions on an unstructured mesh. International
!! Journal for Numerical Methods in Engineering, 38:1757–1776, 1995.
!!
!! - Y.D. Fryer, C. Bailey, M. Cross, and C.H. Lai. A control volume procedure
!! for solving the elastic stress-strain equations on an unstructured mesh.
!! Applied Mathematical Modelling, 15:639–645, 1991.
!!
!! - P.S. Follansbee and U.F. Kocks. A constitutive description of the
!! deformation of copper based on the use of the mechanical threshold stress as
!! an internal state variable. Acta Metallurgica, 36:81–93, 1988.
!!
!! - O. C. Zienkiewicz. The Finite Element Method. McGraw-Hill, New York, NY,
!! 1977.
!!
!! Zach Jibben <zjibben@lanl.gov>
!! August 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module solid_mechanics_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_mesh_type
  use truchas_logging_services
  use truchas_timers
  use integration_geometry_type
  implicit none
  private

  type, public :: solid_mechanics
    private
    type(unstr_mesh), pointer :: mesh => null() ! unowned reference
    type(integration_geometry) :: ig

    !! TODO: These will probably vary across the mesh
    real(r8) :: poissons_ratio, youngs_modulus, thermal_expansion_coeff
    real(r8), allocatable :: density(:), delta_temperature(:)

    !! input parameters
    real(r8), allocatable :: body_force_density(:)
    real(r8) :: contact_distance, contact_normal_traction, contact_penalty, strain_limit

    integer, public :: thermoelastic_niter = 0 ! linear iteration count
    integer, public :: viscoplastic_niter = 0  ! nonlinear iteration count
  contains
    procedure :: init
    procedure :: step
    procedure, private :: compute_residual
    procedure, private :: compute_total_strain
  end type solid_mechanics

contains

  subroutine init(this, mesh, params)

    use parameter_list_type

    class(solid_mechanics), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(parameter_list), intent(inout) :: params

    type(parameter_list), pointer :: plist => null()

    call start_timer("solid mechanics")

    this%mesh => mesh
    call this%ig%init(mesh)

    call params%get('body-force-density', this%body_force_density, default=[0d0, 0d0, 0d0])
    call params%get('contact-distance', this%contact_distance, default=1e-7_r8)
    call params%get('contact-normal-traction', this%contact_normal_traction, default=1e4_r8)
    call params%get('contact-penalty', this%contact_penalty, default=1e3_r8)
    call params%get('strain-limit', this%strain_limit, default=1e-10_r8)

    ASSERT(size(this%body_force_density) == 3)

    plist => params%sublist("nonlinear-solver")
    ! call this%solver%init(plist)

    call stop_timer("solid mechanics")

  end subroutine init


  subroutine step(this)

    class(solid_mechanics), intent(inout) :: this

    call start_timer("solid mechanics")

    

    call stop_timer("solid mechanics")

  end subroutine step


  ! This evaluates the equations on page 1765 of Bailey & Cross 1995.
  subroutine compute_residual(this, displ, r)

    class(solid_mechanics), intent(in) :: this
    real(r8), intent(in) :: displ(:,:)
    real(r8), intent(out) :: r(:)

    integer :: n, xp, p
    real(r8) :: lhs(3*this%mesh%nnode_onP), rhs(3*this%mesh%nnode_onP), total_strain(6,this%ig%npt)
    real(r8) :: a, b, tmp, tmp2(3), Dmatr(3,6)

    ASSERT(size(displ,dim=1) == 3 .and. size(displ,dim=2) == this%mesh%nnode_onP)

    call this%compute_total_strain(displ, total_strain)

    do n = 1, this%mesh%nnode_onP
      associate (rhsn => rhs(3*(n-1)+1:3*(n-1)+3), &
          lhsn => lhs(3*(n-1)+1:3*(n-1)+3), &
          rn => r(3*(n-1)+1:3*(n-1)+3), &
          np => this%ig%npoint(this%ig%xnpoint(n):this%ig%xnpoint(n+1)-1))

        rhsn = this%body_force_density * this%density(n) * this%ig%volume(n)
        lhsn = 0
        do xp = 1, size(np)
          p = np(xp)

          ! right hand side
          tmp = this%youngs_modulus / (1-2*this%poissons_ratio) * this%thermal_expansion_coeff
          if (btest(this%ig%nppar(n),xp)) tmp = -tmp
          rhsn = rhsn + tmp * this%ig%n(:,p) * this%delta_temperature(p)

          ! left hand side
          tmp = this%youngs_modulus / (2*(1-this%poissons_ratio))
          if (btest(this%ig%nppar(n),xp)) tmp = -tmp
          
          ! Dmatr(1,1) = 2/(2-2*this%poissons_ratio) * ((1-this%poissons_ratio)*grad_displ(1,1,p) + this%poissons_ratio*(grad_displ(2,2,p)+grad_displ(3,3,p)))
          ! Dmatr(2,1) = grad_displ(1,2) + grad_displ(2,1) ! total_strain(4)
          ! Dmatr(3,1) = grad_displ(1,3) + grad_displ(3,1) ! total_strain(5)
          ! Dmatr(1,2) = grad_displ(2,1) + grad_displ(1,2) ! total_strain(4)
          ! Dmatr(2,2) = 2/(2-2*this%poissons_ratio) * ((1-this%poissons_ratio)*grad_displ(2,2,p) + this%poissons_ratio*(grad_displ(1,1,p)+grad_displ(3,3,p)))
          ! Dmatr(3,2) = grad_displ(2,3) + grad_displ(3,2) ! total_strain(6)
          ! Dmatr(1,3) = grad_displ(3,1) + grad_displ(1,3) ! total_strain(5)
          ! Dmatr(2,3) = grad_displ(3,2) + grad_displ(2,3) ! total_strain(6)
          ! Dmatr(3,3) = 2/(2-2*this%poissons_ratio) * ((1-this%poissons_ratio)*grad_displ(3,3,p) + this%poissons_ratio*(grad_displ(1,1,p)+grad_displ(2,2,p)))

          ! Dmatr(1,1) = 2/(2-2*this%poissons_ratio) * ((1-this%poissons_ratio)*total_strain(1,p) + this%poissons_ratio*(total_strain(2,p)+total_strain(3,p)))
          ! Dmatr(2,2) = 2/(2-2*this%poissons_ratio) * ((1-this%poissons_ratio)*total_strain(2,p) + this%poissons_ratio*(total_strain(1,p)+total_strain(3,p)))
          ! Dmatr(3,3) = 2/(2-2*this%poissons_ratio) * ((1-this%poissons_ratio)*total_strain(3,p) + this%poissons_ratio*(total_strain(1,p)+total_strain(2,p)))
          ! Dmatr(2,1) = total_strain(4)
          ! Dmatr(1,2) = total_strain(4)
          ! Dmatr(3,1) = total_strain(5)
          ! Dmatr(1,3) = total_strain(5)
          ! Dmatr(3,2) = total_strain(6)
          ! Dmatr(2,3) = total_strain(6)
          ! lhsn = lhsn + tmp * matmul(Dmatr, this%ig%n(:,p))

          Dmatr = 0
          a = 2/(1-2*this%poissons_ratio) * (1-this%poissons_ratio)
          b = 2/(1-2*this%poissons_ratio) * this%poissons_ratio
          Dmatr(:,1) = this%ig%n(:,p) * [a, b, b]
          Dmatr(:,2) = this%ig%n(:,p) * [b, a, b]
          Dmatr(:,3) = this%ig%n(:,p) * [b, b, a]
          Dmatr(1,4) = this%ig%n(2,p)
          Dmatr(1,5) = this%ig%n(3,p)
          Dmatr(2,4) = this%ig%n(1,p)
          Dmatr(2,6) = this%ig%n(3,p)
          Dmatr(3,5) = this%ig%n(1,p)
          Dmatr(3,6) = this%ig%n(2,p)
          lhsn = lhsn + tmp * matmul(Dmatr, total_strain(:,p))
        end do

        rn = lhsn - rhsn
      end associate
    end do
    
  end subroutine compute_residual


  !! From the node-centered displacement, compute the integration-point-centered total strain
  subroutine compute_total_strain(this, displ, total_strain)

    class(solid_mechanics), intent(in) :: this
    real(r8), intent(in) :: displ(:,:)
    real(r8), intent(out) :: total_strain(:,:)

    integer :: p, j, d
    real(r8) :: grad_displ(3,3)

    do j = 1, this%mesh%ncell
      associate (cn => this%mesh%cnode(this%mesh%xcnode(j):this%mesh%xcnode(j+1)-1))
        do p = this%ig%xcpoint(j), this%ig%xcpoint(j+1)-1
          !! Compute the derivative of a scalar phi in global coordinates using the
          !! formula on page 1765 of Bailey & Cross 1995--linear interpolation.
          !! The Jacobian inverse multiplication is already performed and stored
          !! in the grad_shape matrix.
          do d = 1, 3
            grad_displ(:,d) = matmul(this%ig%grad_shape(p)%p, displ(d,cn))
          end do

          total_strain(1,p) = grad_displ(1,1)
          total_strain(2,p) = grad_displ(2,2)
          total_strain(3,p) = grad_displ(3,3)
          total_strain(4,p) = grad_displ(1,2) + grad_displ(2,1)
          total_strain(5,p) = grad_displ(1,3) + grad_displ(3,1)
          total_strain(6,p) = grad_displ(2,3) + grad_displ(3,2)
        end do
      end associate
    end do

  end subroutine compute_total_strain


  ! !! Defines the temperature-dependent density function for each immobile
  ! !! material phase using the specified reference density and temperature
  ! !! values and the linear CTE function.
  ! subroutine define_tm_density_property(stat, errmsg)

  !   use scalar_func_class
  !   use tm_density, only: alloc_tm_density_func

  !   integer, intent(out) :: stat
  !   character(*), intent(out) :: errmsg
  !   character(:), allocatable :: errm

  !   integer :: m, p
  !   real(r8) :: dens0, temp0, state(0)
  !   class(scalar_func), allocatable :: cte_fun, rho_fun

  !   do p = 1, matl_model%nphase_real
  !     if (matl_model%is_fluid(p)) cycle
  !     dens0 = matl_model%const_phase_prop(p, 'tm-ref-density')
  !     temp0 = matl_model%const_phase_prop(p, 'tm-ref-temp')
  !     call matl_model%alloc_phase_prop(p, 'tm-linear-cte', cte_fun)
  !     ASSERT(allocated(cte_fun))

  !     !! Create the temperature-dependent density function.
  !     call alloc_tm_density_func(rho_fun, dens0, temp0, cte_fun, stat, errm)
  !     if (stat /= 0) then
  !       errmsg = 'problem with the "tm-linear-cte" property: ' // trim(errm)
  !       return
  !     end if

  !     !! Assign the density function as the TM density property.
  !     call matl_model%add_phase_prop(p, 'TM density', rho_fun)
  !   end do

  !   stat = 0
  !   errmsg = ''

  ! end subroutine define_tm_density_property


  ! !! This computes the specified property for the cells on the original Truchas
  ! !! mesh. The property is one given in a PHASE namelist PROPERTY_NAME variable,
  ! !! or one created internally by Truchas, and the associated property value is
  ! !! either constant or a function of temperature only.  The routine essentially
  ! !! handles the material averaging of the property over a cell using the volume
  ! !! fraction data from MATL.  The void material (one with a MATERIAL namelist
  ! !! density of zero) contributes nothing to the property value; e.g., zero is
  ! !! returned for an entirely void cell.
  ! subroutine compute_cell_property(prop, temp, mesh, value)

  !   use matl_module, only: gather_vof
  !   use scalar_func_class
  !   use scalar_func_tools, only: is_const

  !   character(*), intent(in) :: prop
  !   real(r8), intent(in) :: temp(:)
  !   type(unstr_mesh), intent(in) :: mesh
  !   real(r8), intent(out) :: value(:)

  !   integer :: m, j
  !   real(r8) ::vofm(mesh%ncell), state(1), mval
  !   class(scalar_func), allocatable :: prop_fun

  !   ASSERT(size(temp) == mesh%ncell)
  !   ASSERT(size(value) == mesh%ncell)

  !   value = 0.0_r8
  !   do m = 1, matl_model%nphase_real
  !     if (matl_model%is_fluid(m)) cycle
  !     call matl_model%alloc_phase_prop(m, prop, prop_fun)
  !     ASSERT(allocated(prop_fun))
  !     call gather_vof(m, vofm)
  !     if (is_const(prop_fun)) then
  !       mval = prop_fun%eval(state)  ! state is ignored, but required
  !       value = value + mval*vofm
  !     else
  !       do j = 1, mesh%ncell
  !         if (vofm(j) > 0.0_r8) then
  !           state(1) = temp(j)
  !           mval = prop_fun%eval(state)
  !           value(j) = value(j) + mval*vofm(j)
  !         end if
  !       end do
  !     end if
  !   end do

  ! end subroutine compute_cell_property

end module solid_mechanics_type
