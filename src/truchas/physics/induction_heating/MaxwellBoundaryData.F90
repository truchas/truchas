!!
!! The MaxwellBoundaryData Module
!!
!! Neil N Carlson <nnc@newmexico.com> 10 Jul 2003
!! Last revised 19 Jul 2003
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module MaxwellBoundaryData

  use kinds, only: rk => r8
  use simpl_mesh_type
  implicit none
  private
  
  public :: create_boundary_data, destroy_boundary_data
  public :: set_Eb_function, set_Hb_function, set_Eb_values, get_Hb_source
  public :: zero_field, xhat_field, yhat_field, zhat_field
  
  public :: generate_bface, bit_mask
  
  type, public :: BoundaryData
    type(simpl_mesh), pointer :: mesh => null()
    !! Boundary conditions on nxE (essential)
    integer :: nebgroup = 0
    integer, pointer :: ebedge(:) => null()
    real(kind=rk), pointer :: eb(:) => null()
    !! Boundary conditions on nxH (natural)
    integer :: nhbgroup = 0
    integer, pointer :: hbface(:) => null()
    real(kind=rk), pointer :: hb(:,:) => null()
  end type BoundaryData

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! CREATE_BOUNDARY_DATA
 !!

  subroutine create_boundary_data (self, mesh, bface, ebgroup, hbgroup)
  
    use index_partitioning, only: gather_boundary
  
    type(BoundaryData), intent(out) :: self
    type(simpl_mesh), intent(in), target :: mesh
    integer, intent(in) :: bface(:)
    integer, intent(in) :: ebgroup(:)
    integer, intent(in), optional :: hbgroup(:)
    
    integer :: i, j, k
    
    ASSERT( size(bface) == mesh%nface )
    ASSERT( all(ebgroup > 0) )
    ! Fails in the case of a 0-sized bface
    !ASSERT( all(ebgroup <= maxval(abs(bface))) )
    
    self%mesh => mesh ! we'll need this for subsequent initialization steps
    
    allocate(self%ebedge(mesh%nedge), self%eb(mesh%nedge))
    self%ebedge = 0
    self%eb = 0.0_rk
    self%nebgroup = size(ebgroup)
    do i = 1, self%nebgroup
      do j = 1, mesh%nface
        if (abs(bface(j)) == ebgroup(i)) then
          do k = 1, 3
            if (self%ebedge(abs(mesh%fedge(k,j))) == 0) self%ebedge(abs(mesh%fedge(k,j))) = i
          end do
        end if
      end do
    end do
    
    !! Ensure consistency of values for edges on the partition boundary.
    call gather_boundary (mesh%edge_ip, self%ebedge)
    
    if (present(hbgroup)) then
      ASSERT( all(hbgroup > 0) )
      ! Assertion fails for 0-sized bface
      !ASSERT( all(hbgroup <= maxval(abs(bface))) )
      allocate(self%hbface(mesh%nface), self%hb(mesh%nedge,size(hbgroup)))
      self%hbface = 0
      self%nhbgroup = size(hbgroup)
      do i = 1, self%nhbgroup
        do j = 1, mesh%nface
          if (abs(bface(j)) == hbgroup(i)) then
            if (self%hbface(j) == 0) self%hbface(j) = sign(i, bface(j))
          end if
        end do
      end do
    end if
    
  end subroutine create_boundary_data
  
  subroutine destroy_boundary_data (self)
    type(BoundaryData), intent(inout) :: self
    if (associated(self%ebedge)) deallocate(self%ebedge)
    if (associated(self%eb)) deallocate(self%eb)
    if (associated(self%hbface)) deallocate(self%hbface)
    if (associated(self%hb)) deallocate(self%hb)
    self%mesh => null()
    self%nebgroup = 0
    self%nhbgroup = 0
  end subroutine destroy_boundary_data
  
  subroutine set_Eb_function (self, n, f, order)
  
    type(BoundaryData), intent(inout) :: self
    integer, intent(in) :: n
    integer, intent(in), optional :: order
    
    interface
      function f (x) result (fx)
        use kinds, only: rk => r8
        real(kind=rk), intent(in) :: x(:)
        real(kind=rk) :: fx(3)
      end function f
    end interface
    
    integer :: j, ord
    real(kind=rk), pointer :: x1(:), x2(:)
    
    ASSERT( associated(self%ebedge) )
    ! Assertion fails for 0-sized ebedge array
    !ASSERT( n > 0 .and. n <= maxval(self%ebedge) )
    
    ord = 1 ! use trapezoid rule quadrature by default
    if (present(order)) ord = order
    
    do j = 1, self%mesh%nedge
      if (self%ebedge(j) /= n) cycle
      x1 => self%mesh%x(:,self%mesh%enode(1,j))
      x2 => self%mesh%x(:,self%mesh%enode(2,j))
      self%eb(j) = dot_product(x2 - x1, vector_average(x1, x2, f, ord))
    end do
    
  end subroutine set_Eb_function
  
  subroutine set_Eb_values (self, coef, e)
  
    type(BoundaryData), intent(in) :: self
    real(kind=rk), intent(in) :: coef(:)
    real(kind=rk), intent(inout) :: e(:)
    
    integer :: j
    
    ASSERT( size(coef) == self%nebgroup)
    ASSERT( size(e) == size(self%eb) )
    
    do j = 1, size(e)
      if (self%ebedge(j) == 0) cycle
      e(j) = coef(self%ebedge(j)) * self%eb(j)
    end do
    
  end subroutine set_Eb_values
  
  subroutine set_Hb_function (self, n, f)
  
    type(BoundaryData), intent(inout) :: self
    integer, intent(in) :: n
    
    integer :: j, edge(3)
    real(kind=rk) :: h(3), c
    
    interface
      function f (x) result (fx)
        use kinds, only: rk => r8
        real(kind=rk), intent(in) :: x(:)
        real(kind=rk) :: fx(3)
      end function f
    end interface
    
    ASSERT( n > 0 .and. n <= self%nhbgroup )
    
    self%hb(:,n) = 0.0_rk
    do j = 1, self%mesh%nface
      if (abs(self%hbface(j)) /= n) cycle
      edge = self%mesh%fedge(:,j)
      c = sign(1, self%hbface(j)) / 6.0_rk
      call project_on_face_edges (self%mesh%x(:,self%mesh%fnode(:,j)), f, h)
      self%hb(edge(1),n) = self%hb(edge(1),n) + c * (h(3) - h(2))
      self%hb(edge(2),n) = self%hb(edge(2),n) - c * (h(1) - h(3))
      self%hb(edge(3),n) = self%hb(edge(3),n) + c * (h(2) - h(1))
    end do
    
  end subroutine set_Hb_function
  
  subroutine project_on_face_edges (x, f, p)
  
    real(kind=rk), intent(in) :: x(:,:)
    real(kind=rk), intent(out):: p(:)
    
    interface
      function f (x) result (fx)
        use kinds, only: rk => r8
        real(kind=rk), intent(in) :: x(:)
        real(kind=rk) :: fx(3)
      end function f
    end interface
    
    integer :: k
    real(kind=rk) :: fv(3,3)
    
    ASSERT( size(x,dim=1) == 3 )
    ASSERT( size(x,dim=2) == 3 )
    ASSERT( size(p) == 3 )
    
    !! Length-weighted edge projections using trapezoid rule
    do k = 1, 3
      fv(:,k) = f(x(:,k))
    end do
    p(1) = 0.5_rk * dot_product(x(:,3)-x(:,2), fv(:,3)+fv(:,2))
    p(2) = 0.5_rk * dot_product(x(:,1)-x(:,3), fv(:,1)+fv(:,3))
    p(3) = 0.5_rk * dot_product(x(:,2)-x(:,1), fv(:,2)+fv(:,1))
    
  end subroutine project_on_face_edges

  function vector_average (x1, x2, f, order) result (avg)
    real(kind=rk), intent(in) :: x1(:), x2(:)
    integer, intent(in) :: order
    real(kind=rk) :: avg(size(x1))
    interface
      function f (x) result (fx)
        use kinds, only: rk => r8
        real(kind=rk), intent(in) :: x(:)
        real(kind=rk) :: fx(3)
      end function f
    end interface
    select case (order)
    case (:0) ! midpoint rule
      avg = f(0.5_rk*(x1+x2))
    case (1)  ! trapezoid rule
      avg = 0.5_rk*(f(x1)+f(x2))
    case (2:) ! Simpson's rule
      avg = (f(x1)+4.0_rk*f(0.5_rk*(x1+x2))+f(x2))/6.0_rk
    end select
  end function vector_average
  
  
  subroutine get_Hb_source (self, coef, bsrc)
  
    type(BoundaryData), intent(in) :: self
    real(kind=rk), intent(in) :: coef(:)
    real(kind=rk), intent(out) :: bsrc(:)
    
    integer :: j, k
    
    ASSERT( size(coef) == size(self%hb,dim=2) )
    ASSERT( size(bsrc) == size(self%hb,dim=1) )
    
    bsrc = 0.0_rk
    do k = 1, size(coef)
      do j = 1, size(bsrc)
        bsrc(j) = bsrc(j) + coef(k) * self%hb(j,k)
      end do
    end do
    
    !! We need to suppress values on any edges where nxE is specified.
    where (self%ebedge /= 0) bsrc = 0.0_rk
    
  end subroutine get_Hb_source
  
  !!
  !! THESE TRIVIAL FIELDS ARE PROVIDED FOR THE USER'S CONVENIENCE
  !!
  
  function zero_field (x) result (v)
    real(kind=rk), intent(in) :: x(:)
    real(kind=rk) :: v(3)
    v = 0.0_rk
  end function zero_field
  
  function xhat_field (x) result (v)
    real(kind=rk), intent(in) :: x(:)
    real(kind=rk) :: v(3)
    v = (/ 1.0_rk, 0.0_rk, 0.0_rk /)
  end function xhat_field
  
  function yhat_field (x) result (v)
    real(kind=rk), intent(in) :: x(:)
    real(kind=rk) :: v(3)
    v = (/ 0.0_rk, 1.0_rk, 0.0_rk /)
  end function yhat_field
  
  function zhat_field (x) result (v)
    real(kind=rk), intent(in) :: x(:)
    real(kind=rk) :: v(3)
    v = (/ 0.0_rk, 0.0_rk, 1.0_rk /)
  end function zhat_field
  
      
  subroutine generate_bface (gm, mesh, filter, group, bface, status)
  
    use GeometricModeler
    use truchas_env, only: output_dir
    use parallel_communication, only: global_any
    use truchas_logging_services
    
    type(GeometricModel), intent(in) :: gm
    type(simpl_mesh), intent(in) :: mesh
    logical, intent(in) :: filter(:)
    integer, intent(in) :: group(:)
    integer, intent(out) :: bface(:)
    integer, intent(out), optional :: status
    
    integer :: i, j
    integer, pointer :: surf(:)
    logical :: mask(mesh%nface)
    character(:), allocatable :: vtk_file
    
    ASSERT( size(filter) == mesh%nface )
    ASSERT( size(bface) == mesh%nface )
    
    bface = 0
    do j = 1, mesh%nface
      if (.not.filter(j)) cycle ! ignore this face
      surf => OnSurface(gm, mesh%x(:,mesh%fnode(:,j)))
      select case (size(surf))
      case (1)
        if (surf(1) >= bit_size(group)) then
          print *, 'GENERATE_BFACE: PANIC! mask bit width overflow'
          stop
        else if (surf(1) < 1) then
          print *, 'GENERATE_BFACE: PANIC! mask bit width underflow'
          stop
        else
          do i = 1, size(group)
            if (btest(group(i), surf(1))) then
              if (same_normal(gm, surf(1), mesh%x(:,mesh%fnode(:,j)))) then
                bface(j) = i
              else
                bface(j) = -i
              end if
              exit
            end if
          end do
        end if
      case (2:)
        print *, 'GENERATE_BFACE: PANIC! face associated with more than one surface'
        print *, '  face:', mesh%x(:,mesh%fnode(:,j))
        print *, '  surfaces:', surf
        stop
      end select
      deallocate(surf)
    end do
    
    mask = ((bface == 0) .and. filter)
    if (global_any(mask)) then
      if (present(status)) then
        status = -1
        return
      else
        vtk_file = trim(output_dir) // 'em_no_bc.vtk'
        call mesh%write_faces_vtk(mask, vtk_file, 'Boundary faces without BC')
        call TLS_fatal('GENERATE_BFACE: no BC for some EM mesh boundary faces;&
            & problem faces written to ' // vtk_file)
      end if
    end if
    
    if (present(status)) status = 0
    
  end subroutine generate_bface
  
  function bit_mask (list) result (n)
  integer, intent(in) :: list(:)
  integer :: n, j
  n = 0
  do j = 1, size(list)
    n = ibset(n, list(j))
  end do
  end function bit_mask
  
end module MaxwellBoundaryData
