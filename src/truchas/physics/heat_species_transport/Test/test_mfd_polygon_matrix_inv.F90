!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_mfd_polygon_matrix_inv
#ifdef NAGFOR
  use,intrinsic :: f90_unix, only: exit
#endif
  use kinds
  use mfd_disc_type
  use cell_geometry

  implicit none

  integer :: status = 0

  type(mfd_cell) :: hex_cell, pyr_cell, prism_cell

  call init_hex(hex_cell)
  call test_mass_matrix_inv(hex_cell)
!!$
!!$  call init_prism(prism_cell)
!!$  call test_mass_matrix_inv(prism_cell)

contains


  subroutine test_mass_matrix_inv(polygon_cell)

    type(mfd_cell) :: polygon_cell
    real(r8), allocatable :: matrix(:)
    real(r8) :: coef, diff, grad_norm, sur_int, sur_int_num
    real(r8) :: grad(3), u_center 
    real(r8), allocatable :: du(:), flux(:), flux_num(:)
    integer :: i, j, loc
    integer,parameter :: seed = 86456


    allocate(matrix(polygon_cell%nfaces*(1 + polygon_cell%nfaces)/2))
    coef = 1.

    call polygon_cell%compute_flux_matrix_inv(coef, matrix)

    print *, matrix(:)

    allocate(  du(polygon_cell%nfaces))
    allocate(flux(polygon_cell%nfaces))
    allocate(flux_num(polygon_cell%nfaces))

    !call srand(seed)
    call random_number(grad)
    
!!$!    grad(1) = 1.
!!$!    grad(2) = 6.
!!$!    grad(3) = -3.

    print *, "gradient ", grad
    grad_norm = 0
    do i=1,3
       grad_norm = grad_norm + grad(i)*grad(i)
    end do
    grad_norm = sqrt(grad_norm)

    u_center = 0
    do i=1, 3
       u_center = u_center + grad(i)*polygon_cell%cell_center(i)
    end do

    do i=1, polygon_cell%nfaces
       du(i) = -u_center
       do j=1, 3
          du(i) = du(i) + grad(j)*polygon_cell%face_centers(j, i)
       end do
       du(i) = du(i)*polygon_cell%face_area(i)

       flux(i) = 0.
       do j=1, 3
          flux(i) = flux(i) + coef*grad(j)*polygon_cell%face_normals(j,i)
       end do
    end do
   
    diff = 0.
    do i=1, polygon_cell%nfaces
       flux_num(i) = 0.
       do j=1, polygon_cell%nfaces
          if (i<j) then
             loc = i + j*(j-1)/2
          else
             loc = j + i*(i-1)/2
          end if
          flux_num(i) = flux_num(i) + matrix(loc)*du(j) 
       end do
       print *, "fl", flux(i)*polygon_cell%face_area(i), "fl_n", flux_num(i)
       diff = max( diff, abs(flux_num(i) - flux(i)*polygon_cell%face_area(i) ) )
    end do

    print *, "Error: ", diff

    sur_int = 0.;
    sur_int_num = 0.;
    do i=1, polygon_cell%nfaces
       sur_int = sur_int + flux(i)*polygon_cell%face_area(i)
       sur_int_num = sur_int_num + flux_num(i)*polygon_cell%face_area(i)
    end do

    print *, "Volume integral:"
    print *,  sur_int, sur_int_num

    
    if (diff < 1e-6*grad_norm) print *, "TEST: exit successfully ************"

  end subroutine  test_mass_matrix_inv



  subroutine init_hex(hex_cell)

    type(mfd_cell), intent(inout) :: hex_cell
    real(r8) :: vertices(3,8)
    integer :: cnodes(8), i

!!$    cnodes = [ (i, i = 1, 8) ]
!!$
!!$    vertices(1, 1) = 0
!!$    vertices(2, 1) = 0
!!$    vertices(3, 1) = 0
!!$
!!$    vertices(1, 2) = 2
!!$    vertices(2, 2) = 0
!!$    vertices(3, 2) = 0
!!$
!!$    vertices(1, 3) = 2
!!$    vertices(2, 3) = 2
!!$    vertices(3, 3) = 0
!!$
!!$    vertices(1, 4) = 0
!!$    vertices(2, 4) = 2
!!$    vertices(3, 4) = 0
!!$
!!$    vertices(1, 5) = 0
!!$    vertices(2, 5) = 0
!!$    vertices(3, 5) = 3
!!$
!!$    vertices(1, 6) = 2
!!$    vertices(2, 6) = 0
!!$    vertices(3, 6) = 3
!!$
!!$    vertices(1, 7) = 2
!!$    vertices(2, 7) = 2
!!$    vertices(3, 7) = 3
!!$
!!$    vertices(1, 8) = 0
!!$    vertices(2, 8) = 2
!!$    vertices(3, 8) = 3
!!$
!!$
!!$    print *, "Create: HEX cell"

!    call hex_cell%init(cnodes, vertices)

!    print *, cell_center(vertices)
!!$    print *, hex_cell%nfaces
!!$    print *, hex_cell%volume
!!$    print *, hex_cell%cell_center
!!$
!!$    return

    hex_cell%nfaces = 6

    allocate(hex_cell%face_area(hex_cell%nfaces))
    allocate(hex_cell%face_centers(3, hex_cell%nfaces))
    allocate(hex_cell%face_normals(3, hex_cell%nfaces))
    allocate(hex_cell%cell_center(3))


    hex_cell%volume = 8.00E-03
    hex_cell%cell_center(1) = 3.00E-01
    hex_cell%cell_center(2) = 0.
    hex_cell%cell_center(3) = 0.



    hex_cell%face_area = 4.00E-02
    hex_cell%face_centers(1,1) = 3.00E-01
    hex_cell%face_centers(2,1) = -1.00E-01
    hex_cell%face_centers(3,1) = 0

    hex_cell%face_centers(1,2) = 3.00E-01
    hex_cell%face_centers(2,2) = 0
    hex_cell%face_centers(3,2) = -1.00E-01 

    hex_cell%face_centers(1,3) = 3.00E-01
    hex_cell%face_centers(2,3) = 1.00E-01
    hex_cell%face_centers(3,3) = 0

    hex_cell%face_centers(1,4) = 3.00E-01
    hex_cell%face_centers(2,4) = 0
    hex_cell%face_centers(3,4) = 1.00E-01 

    hex_cell%face_centers(1,5) = 2.00E-01
    hex_cell%face_centers(2,5) = 0
    hex_cell%face_centers(3,5) = 0

    hex_cell%face_centers(1,6) = 4.00E-01
    hex_cell%face_centers(2,6) = 0
    hex_cell%face_centers(3,6) = 0


    hex_cell%face_normals(1,1) = 0
    hex_cell%face_normals(2,1) = -4.00E-02 
    hex_cell%face_normals(3,1) = 0

    hex_cell%face_normals(1,2) = 0
    hex_cell%face_normals(2,2) = 0
    hex_cell%face_normals(3,2) = -4.00E-02 

    hex_cell%face_normals(1,3) = 0
    hex_cell%face_normals(2,3) = 4.00E-02
    hex_cell%face_normals(3,3) = 0

    hex_cell%face_normals(1,4) = 0
    hex_cell%face_normals(2,4) = 0
    hex_cell%face_normals(3,4) = 4.00E-02 

    hex_cell%face_normals(1,5) = -4.00E-02
    hex_cell%face_normals(2,5) = 0
    hex_cell%face_normals(3,5) = 0

    hex_cell%face_normals(1,6) = 4.00E-02
    hex_cell%face_normals(2,6) = 0
    hex_cell%face_normals(3,6) = 0

!!$    do i = 1, 6
!!$       hex_cell%face_normals(:,i) = hex_cell%face_normals(:,i) /  hex_cell%face_area(i)
!!$    end do

  end subroutine init_hex


!!$  subroutine init_pyramid(pyr_cell)
!!$
!!$    type(mfd_polygon), intent(inout) :: pyr_cell
!!$
!!$    print *, "Create: PYRAMID cell"
!!$    
!!$    pyr_cell%nfaces = 5
!!$    pyr_cell%volume = 8.
!!$    pyr_cell%cell_center(1) = 1.
!!$    pyr_cell%cell_center(2) = 1.
!!$    pyr_cell%cell_center(3) = 1.
!!$
!!$    allocate(pyr_cell%face_area(pyr_cell%nfaces))
!!$    allocate(pyr_cell%face_centers(3, pyr_cell%nfaces))
!!$    allocate(pyr_cell%face_normals(3, pyr_cell%nfaces))
!!$
!!$    pyr_cell%face_area = 4.
!!$    pyr_cell%face_centers(1,1) = 1
!!$    pyr_cell%face_centers(2,1) = 0
!!$    pyr_cell%face_centers(3,1) = 1
!!$
!!$    pyr_cell%face_centers(1,2) = 2
!!$    pyr_cell%face_centers(2,2) = 1
!!$    pyr_cell%face_centers(3,2) = 1
!!$
!!$    pyr_cell%face_centers(1,3) = 1
!!$    pyr_cell%face_centers(2,3) = 2
!!$    pyr_cell%face_centers(3,3) = 1
!!$
!!$    pyr_cell%face_centers(1,4) = 0
!!$    pyr_cell%face_centers(2,4) = 1
!!$    pyr_cell%face_centers(3,4) = 1
!!$
!!$    pyr_cell%face_centers(1,5) = 1
!!$    pyr_cell%face_centers(2,5) = 1
!!$    pyr_cell%face_centers(3,5) = 0
!!$
!!$    pyr_cell%face_centers(1,6) = 1
!!$    pyr_cell%face_centers(2,6) = 1
!!$    pyr_cell%face_centers(3,6) = 2
!!$
!!$
!!$    pyr_cell%face_normals(1,1) = 0
!!$    pyr_cell%face_normals(2,1) = -1.
!!$    pyr_cell%face_normals(3,1) = 0
!!$
!!$    pyr_cell%face_normals(1,2) = 1
!!$    pyr_cell%face_normals(2,2) = 0
!!$    pyr_cell%face_normals(3,2) = 0
!!$
!!$    pyr_cell%face_normals(1,3) = 0
!!$    pyr_cell%face_normals(2,3) = 1
!!$    pyr_cell%face_normals(3,3) = 0
!!$
!!$    pyr_cell%face_normals(1,4) = -1
!!$    pyr_cell%face_normals(2,4) = 0
!!$    pyr_cell%face_normals(3,4) = 0
!!$
!!$    pyr_cell%face_normals(1,5) = 0
!!$    pyr_cell%face_normals(2,5) = 0
!!$    pyr_cell%face_normals(3,5) = -1
!!$
!!$    pyr_cell%face_normals(1,6) = 0
!!$    pyr_cell%face_normals(2,6) = 0
!!$    pyr_cell%face_normals(3,6) = 1
!!$
!!$
!!$  end subroutine init_pyramid
!!$
!!$  subroutine init_prism(prism_cell)
!!$
!!$    type(mfd_polygon), intent(inout) :: prism_cell
!!$
!!$    print *, "Create: PRISM cell"
!!$    
!!$    prism_cell%nfaces = 5
!!$    prism_cell%volume = 4.
!!$    prism_cell%cell_center(1) = 1
!!$    prism_cell%cell_center(2) = 2./3
!!$    prism_cell%cell_center(3) = 1.
!!$
!!$    allocate(prism_cell%face_area(prism_cell%nfaces))
!!$    allocate(prism_cell%face_centers(3, prism_cell%nfaces))
!!$    allocate(prism_cell%face_normals(3, prism_cell%nfaces))
!!$
!!$    prism_cell%face_area(1) = 2*sqrt(5.)
!!$    prism_cell%face_area(2) = 2*sqrt(5.)
!!$    prism_cell%face_area(3) = 4
!!$    prism_cell%face_area(4) = 2
!!$    prism_cell%face_area(5) = 2
!!$
!!$    prism_cell%face_centers(1,1) = 1.5
!!$    prism_cell%face_centers(2,1) = 1
!!$    prism_cell%face_centers(3,1) = 1
!!$
!!$    prism_cell%face_centers(1,2) = 0.5
!!$    prism_cell%face_centers(2,2) = 1
!!$    prism_cell%face_centers(3,2) = 1
!!$
!!$    prism_cell%face_centers(1,3) = 1
!!$    prism_cell%face_centers(2,3) = 0
!!$    prism_cell%face_centers(3,3) = 1
!!$
!!$    prism_cell%face_centers(1,4) = 1
!!$    prism_cell%face_centers(2,4) = 2./3.
!!$    prism_cell%face_centers(3,4) = 0
!!$
!!$    prism_cell%face_centers(1,5) = 1.
!!$    prism_cell%face_centers(2,5) = 2./3.
!!$    prism_cell%face_centers(3,5) = 2
!!$
!!$
!!$    prism_cell%face_normals(1,1) = 2./sqrt(5.)
!!$    prism_cell%face_normals(2,1) = 1./sqrt(5.)
!!$    prism_cell%face_normals(3,1) = 0
!!$
!!$    prism_cell%face_normals(1,2) = -2./sqrt(5.)
!!$    prism_cell%face_normals(2,2) =  1./sqrt(5.)
!!$    prism_cell%face_normals(3,2) = 0
!!$
!!$    prism_cell%face_normals(1,3) = 0
!!$    prism_cell%face_normals(2,3) = -1
!!$    prism_cell%face_normals(3,3) = 0
!!$
!!$    prism_cell%face_normals(1,4) = 0
!!$    prism_cell%face_normals(2,4) = 0
!!$    prism_cell%face_normals(3,4) = -1
!!$
!!$    prism_cell%face_normals(1,5) = 0
!!$    prism_cell%face_normals(2,5) = 0
!!$    prism_cell%face_normals(3,5) = 1
!!$
!!$
!!$  end subroutine init_prism


end program test_mfd_polygon_matrix_inv
