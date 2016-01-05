!!
!! CELL_GRAD_TYPE
!!
!! An object that computes a cell-based gradient of a cell-based scalar field
!! using a mimetic finite difference technique.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! May 2014
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module cell_grad_type

  use kinds, only: r8
  use unstr_mesh_type
  use pcsr_matrix_type
  use mfd_disc_type
  use hypre_hybrid_type
  use parameter_list_type
  use index_partitioning
  use truchas_logging_services
  implicit none
  private
  
  type, public :: cell_grad
    type(unstr_mesh), pointer :: mesh => null()  ! reference only -- do not own
    type(mfd_disc),  pointer :: disc => null()  ! reference only -- do not own
    type(pcsr_matrix) :: matrix
    type(hypre_hybrid) :: solver
    type(parameter_list), pointer :: params => null()    ! own
    logical, allocatable :: cell_mask(:)
    logical, allocatable :: active
  contains
    procedure :: init
    procedure :: compute
    final :: delete_cell_grad
  end type cell_grad

contains

  !! Final subroutine for CELL_GRAD objects.
  subroutine delete_cell_grad (this)
    type(cell_grad), intent(inout) :: this
    if (associated(this%params)) deallocate(this%params)
  end subroutine delete_cell_grad

  subroutine init (this, disc, mask, setids, stat, errmsg)

    class(cell_grad), intent(out) :: this
    type(mfd_disc), intent(in), target :: disc
    logical, intent(in) :: mask(:)
    integer, intent(in) :: setids(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    
    integer :: j, l, ic, ir
    logical, allocatable :: face_mask(:)
    type(pcsr_graph), pointer :: g
    type(ip_desc), pointer :: row_ip
    real(r8) :: c
    
    ASSERT(size(mask) == disc%mesh%ncell)
    
    this%disc => disc
    this%mesh => disc%mesh
    this%cell_mask = mask
    
    !! Define the face flux matrix for the subdomain defined by the active
    !! cell mask array.  Decoupled dummy equations are defined for inactive
    !! faces not associated with the subdomain, to yield a system over all
    !! mesh faces.  Only on-process face equations are significant.

    !! Create the CSR matrix graph for the matrix.
    allocate(g, face_mask(this%mesh%nface))
    face_mask = .true.  ! will tag inactive faces
    row_ip => this%mesh%face_ip
    call g%init (row_ip)
    do j = 1, this%mesh%ncell
      if (this%cell_mask(j)) then
        associate (cface => this%mesh%cface(this%mesh%xcface(j):this%mesh%xcface(j+1)-1))
          call g%add_clique (cface)
          face_mask(cface) = .false.  ! tag faces as active
        end associate
      end if
    end do
    do j = 1, this%mesh%nface  ! dummy equations for inactive faces
      if (face_mask(j)) call g%add_edge (j, j)
    end do
    call g%add_complete
    
    !! Create the face flux matrix.
    call this%matrix%init (g, take_graph=.true.)
    call this%matrix%set_all (0.0_r8)
    
    !! Matrix values for inactive face dummy equations.
    do j = 1, this%mesh%nface
      if (face_mask(j)) call this%matrix%set (j, j, 1.0_r8)
    end do
    
    !! Get the face mask identifying faces where passive BC will be applied.
    call get_bc_mask (this%mesh, this%cell_mask, setids, face_mask, stat, errmsg)
    if (stat /= 0) return
    
    !! Assemble the matrix elements corresponding to the active cells.
    do j = 1, this%mesh%ncell
      if (this%cell_mask(j)) then
        associate (index => this%mesh%cface(this%mesh%xcface(j):this%mesh%xcface(j+1)-1), &
                   minv => this%disc%minv(this%disc%xminv(j):this%disc%xminv(j+1)-1))
          l = 1
          do ic = 1, size(index)
            do ir = 1, ic-1
              call this%matrix%add_to (index(ir), index(ic), minv(l))
              call this%matrix%add_to (index(ic), index(ir), minv(l))
              l = l + 1
            end do
            call this%matrix%add_to (index(ic), index(ic), minv(l))
            l = l + 1
            !! Matrix modification for passive BC on face INDEX(IC).  Only
            !! affects that row of the matrix. NB: using IR to iterate over
            !! the elements in the row.
            if (face_mask(index(ic))) then
              !face_normals = hex_face_normals(this%mesh%x(:,this%mesh%cnode(:,j)))
              !do ir = 1, size(index)
              !  c = -dot_product(face_normals(:,ic),face_normals(:,ir))/this%mesh%volume(j)
              !  call this%matrix%add_to (index(ic), index(ir), c)
              !  if (dot_product(face_normals(:,ir),this%mesh%normal(:,index(ir))) > 0.0_r8 &
              !      .eqv. btest(this%mesh%cfpar(j),pos=ir)) then
              !    write(string,'(i0,3(1x,i0),1x,l1)') j, this%mesh%xcell(j), ir, index(ir), index(ir) > this%mesh%nface_onP!, this%mesh%face_set_mask(index(ir))
              !    call TLS_debug (string)
              !  end if
              !end do
              do ir = 1, size(index)
                c = -dot_product(this%mesh%normal(:,index(ic)),this%mesh%normal(:,index(ir)))/this%mesh%volume(j)
                if (btest(this%mesh%cfpar(j),pos=ir)) c = -c
                call this%matrix%add_to (index(ic), index(ir), c)
              end do
            end if
          end do
        end associate
      end if
    end do
    
    !! Setup the solver
    
    !! Hardwire the solver parameters for now; some need to be exposed.
    allocate(this%params)
    call this%params%set ('krylov-method', 'gmres')
    call this%params%set ('rel-tol', 1.0e-10_r8)
    call this%params%set ('gmres-krylov-dim', 5)
    call this%params%set ('max-ds-iter', 50)
    call this%params%set ('max-amg-iter', 20)
    if (TLS_VERBOSITY >= TLS_VERB_NOISY) then
      call this%params%set ('print-level', 1)
      call this%params%set ('logging-level', 1)
    end if
    
    call this%solver%init (this%matrix, this%params)
    call this%solver%setup
  
  end subroutine init
  
  subroutine compute (this, ucell, gradu, stat, errmsg)
  
    use upper_packed_matrix

    class(cell_grad), intent(inout) :: this
    real(r8), intent(in)  :: ucell(:)     ! u on cells (all)
    real(r8), intent(out) :: gradu(:,:)   ! gradient on cells (all or on-process)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    
    integer  :: j, num_itr, num_dscg_itr, num_pcg_itr
    !logical  :: bface(this%mesh%nface)
    real(r8) :: norm
    real(r8) :: uface(this%mesh%nface), rface(this%mesh%nface)
    real(r8), allocatable :: w(:)
    character(80) :: string
    
    ASSERT(size(ucell) == this%mesh%ncell)
    ASSERT(size(gradu,dim=1) == size(this%mesh%x,dim=1))
    ASSERT(size(gradu,dim=2) <= this%mesh%ncell)
    
    !! Compute the RHS for the system of face unknowns and an initial guess for
    !! their solution.  The initial guess is simply the average of neighboring
    !! cell unknowns.  The dummy equations for inactive faces have solution 0.
    !! NB: the initial guess for boundary faces is off by a factor of two.  This
    !! is easily fixed by accumulating a scale factor within the same loop, but
    !! testing didn't show any significant benefit to justify the expense.
    
    rface = 0.0_r8
    uface = 0.0_r8
    !bface = .false.
    do j = 1, this%mesh%ncell
      if (this%cell_mask(j)) then
        associate (cface => this%mesh%cface(this%mesh%xcface(j):this%mesh%xcface(j+1)-1), &
                   minv => this%disc%minv(this%disc%xminv(j):this%disc%xminv(j+1)-1))
          allocate(w(size(cface)))
          call upm_col_sum (minv, w)
          rface(cface) = rface(cface) + ucell(j)*w
          uface(cface) = uface(cface) + ucell(j)*0.5_r8
          !bface(cface) = .not.bface(cface)
          deallocate(w)
        end associate
      end if
    end do
    !where (bface) uface = 2*uface
    
    call this%solver%solve (rface, uface, stat)

    if (TLS_VERBOSITY >= TLS_VERB_NOISY) then
      call this%solver%get_metrics (num_itr, num_dscg_itr, num_pcg_itr, norm)
      write(string,'(3(a,i0),a,es9.2)') 'cell_grad%compute: num_itr = ', num_itr, &
          ' (', num_dscg_itr, ', ', num_pcg_itr, '), ||r||/||b|| = ', norm
      call TLS_info (string)
    end if
    
    if (stat /= 0) then
      call this%solver%get_metrics (num_itr, num_dscg_itr, num_pcg_itr, norm)
      write(string,'(3(a,i0),a,es9.2)') 'failed to converge: num_itr = ', num_itr, &
          ' (', num_dscg_itr, ', ', num_pcg_itr, '), ||r||/||b|| = ', norm
      errmsg = trim(string)
    endif
    
    call gather_boundary (this%mesh%face_ip, uface)
    call this%disc%compute_cell_grad (uface, this%cell_mask(:size(gradu,2)), gradu)

  end subroutine compute
    
  !! Identify boundary faces where BC need to be applied.  A passive BC will
  !! be applied at all active boundary faces NOT belonging to one of the
  !! specified boundary face sets, as well as interior faces separating
  !! active and inactive cells.  At remaining faces on the boundary of the
  !! active domain (i.e., those in the specified face sets) the equations
  !! will not be modified; this corresponds to a natural zero normal gradient
  !! condition.
  !!
  !! NB: This does the necessary communication to ensure a valid result for
  !! all faces.  Is this really necessary (we ignore off-process faces)?
  
  subroutine get_bc_mask (mesh, mask, setids, bc_mask, stat, errmsg)
  
    use bitfield_type
    use string_utilities, only: i_to_c
  
    type(unstr_mesh), intent(in) :: mesh
    logical, intent(in) :: mask(:)
    integer, intent(in) :: setids(:)
    logical, intent(out) :: bc_mask(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    
    integer :: i, j
    type(bitfield) :: bitmask

    !! Identify faces on the boundary of the active domain.  These are
    !! precisely those faces adjacent to exactly one active cell.
    bc_mask = .false.
    do j = 1, mesh%ncell
      if (mask(j)) then
        associate (cface => mesh%cface(mesh%xcface(j):mesh%xcface(j+1)-1))
          bc_mask(cface) = .not. bc_mask(cface)
        end associate
      end if
    end do
    call gather_boundary (mesh%face_ip, bc_mask)
    
    if (size(setids) > 0) then
      !! Create the bitmask corresponding to SETIDS.
      bitmask = ZERO_BITFIELD
      do i = 1, size(setids)
        do j = size(mesh%face_set_ID), 1, -1
          if (setids(i) == mesh%face_set_ID(j)) exit
        end do
        if (j == 0) then
          stat = -1
          errmsg = 'unknown face set ID: ' // i_to_c(setids(i))
          return
        end if
        bitmask = ibset(bitmask, j)
      end do
      !! Eliminate active boundary faces specified by SETIDS.
      where (btest(mesh%face_set_mask,pos=0))
        bc_mask = bc_mask .and. (iand(bitmask, mesh%face_set_mask) == ZERO_BITFIELD)
      end where
    end if
    
    stat = 0
    
  end subroutine get_bc_mask

end module cell_grad_type
