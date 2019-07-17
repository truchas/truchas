!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module test_ht_2d_common

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use truchas_logging_services
  use unstr_2d_mesh_type
  use matl_mesh_func_type
  use material_database_type
  use material_model_type
  use scalar_func_class
  use mfd_2d_disc_type
  use HT_2d_model_type
  implicit none

contains

  !! Initializes material database and related objects needed by the HT types
  subroutine init_materials(mesh, matl_model, mmf)

    use parameter_list_type
    use parameter_list_json
    use material_factory, only: load_material_database
    use material_utilities, only: add_enthalpy_prop

    type(unstr_2d_mesh), target, intent(in) :: mesh
    type(material_model), intent(out) :: matl_model
    type(matl_mesh_func), intent(out) :: mmf

    type(parameter_list), pointer :: plist
    integer, allocatable :: matids(:)
    integer :: stat
    character(:), allocatable :: errmsg, string
    type(material_database), save :: matl_db

    string = '{"unobtanium": &
                {"properties":{"conductivity":1.0, &
                               "density":1.0,&
                               "specific-heat":1.0}}}'

    !! Initialize material database/model with single material
    call parameter_list_from_json_string(string, plist, errmsg)
    call load_material_database(matl_db, plist, stat, errmsg)
    if (stat /= 0) call TLS_FATAL(errmsg)
    call matl_model%init(['unobtanium'], matl_db, stat, errmsg)
    if (stat /= 0) call TLS_FATAL(errmsg)

    !! Initialize enthalpy
    call add_enthalpy_prop(matl_model, stat, errmsg)
    if (stat /= 0) call TLS_FATAL(errmsg)

    !! Layout material across the mesh
    matids = matl_model%matl_index(['unobtanium'])
    call mmf%init(mesh)
    call mmf%define_region(mesh%cell_set_id, matids, stat, errmsg)
    if (stat /= 0) call TLS_FATAL(errmsg)
    call mmf%define_complete(stat, errmsg)
    if (stat /= 0) call TLS_FATAL(errmsg)

  end subroutine init_materials


  !! Computes the average integral of a function over on-process faces and cells
  subroutine average_integral(disc, f, ucell, uface)

    type(mfd_2d_disc), intent(in) :: disc
    class(scalar_func), intent(in) :: f
    real(r8), intent(out) :: ucell(:), uface(:)

    real(r8) :: x0(2) = [sqrt(3.0_r8)/3.0_r8, -sqrt(3.0_r8)/3.0_r8]  ! Quadrature points
    real(r8) :: w(2) = [1.0_r8, 1.0_r8]   ! Quadrature weights
    ! real(r8) :: x0(3) = [sqrt(3.0_r8/5.0_r8), 0.0_r8, -sqrt(3.0_r8/5.0_r8)]  ! Quadrature points
    ! real(r8) :: w(3) = [5.0_r8/9.0_r8, 8.0_r8/9.0_r8, 5.0_r8/9.0_r8]   ! Quadrature weights
    real(r8) :: phi(4,size(x0),size(x0))  ! values of DG shape functions
    real(r8) :: N(4)  ! corner normals
    real(r8) :: val, jac, args(3)
    integer :: j, u, v

    ASSERT(size(ucell) == disc%mesh%ncell_onP)
    ASSERT(size(uface) == disc%mesh%nface_onP)

    args(1) = 0.0_r8  ! f is constant in time

    !! Pre-compute shape functions
    do u = 1, size(x0)
      do v = 1, size(x0)
        phi(1,u,v) = 0.25_r8*(1.0_r8-x0(u))*(1.0_r8-x0(v))
        phi(2,u,v) = 0.25_r8*(1.0_r8+x0(u))*(1.0_r8-x0(v))
        phi(3,u,v) = 0.25_r8*(1.0_r8+x0(u))*(1.0_r8+x0(v))
        phi(4,u,v) = 0.25_r8*(1.0_r8-x0(u))*(1.0_r8+x0(v))
      end do
    end do

    !! Compute cell average
    do j = 1, disc%mesh%ncell_onP
      associate (cface => disc%mesh%cface(disc%mesh%cstart(j):disc%mesh%cstart(j+1)-1), &
                 nodes => disc%mesh%x(:,disc%mesh%cnode(disc%mesh%cstart(j):disc%mesh%cstart(j+1)-1)))
        !! Corner normal lengths
        N(1) = cross_product_2D(nodes(:,2)-nodes(:,1), nodes(:,4)-nodes(:,1))
        N(2) = cross_product_2D(nodes(:,3)-nodes(:,2), nodes(:,1)-nodes(:,2))
        N(3) = cross_product_2D(nodes(:,4)-nodes(:,3), nodes(:,2)-nodes(:,3))
        N(4) = cross_product_2D(nodes(:,1)-nodes(:,4), nodes(:,3)-nodes(:,4))

        !! Gaussian quadrature
        val = 0.0_r8
        do u = 1, size(x0)
          do v = 1, size(x0)
            !! Jacobian
            jac = 0.25_r8*dot_product(N, phi(:,u,v))

            !! Quadrature point
            args(2) = dot_product(nodes(1,:), phi(:,u,v))
            args(3) = dot_product(nodes(2,:), phi(:,u,v))

            val = val + (f%eval(args) * w(u) * w(v) * jac)
          end do
        end do
        ucell(j) = val / disc%mesh%volume(j)
      end associate
    end do

    !! Compute face average
    do j = 1, disc%mesh%nface_onP
      associate (nodes => disc%mesh%x(:,disc%mesh%fnode(:,j)))
        !! Jacobian
        jac = 0.5_r8*disc%mesh%area(j)

        !! Gaussian quadrature
        val = 0.0_r8
        do u = 1, size(x0)
          !! Quadrature point
          args(2:) = 0.5_r8*(sum(nodes, dim=2) + (nodes(:,2)-nodes(:,1))*x0(u))
          val = val + (f%eval(args) * w(u) * jac)
        end do
        uface(j) = val / disc%mesh%area(j)
      end associate
    end do

  end subroutine average_integral


  !TODO: add to cell_geometry?
  !! Computes A x B
  pure function cross_product_2D (a, b) result (axb)
    real(r8), intent(in) :: a(2), b(2)
    real(r8) :: axb
    axb = a(1)*b(2) - a(2)*b(1)
  end function cross_product_2D

end module test_ht_2d_common
