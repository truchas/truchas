!!
!! VOLUME_TRACKING_DRIVER
!!
!! This code is meant to be deleted and subsumed by high-level cycle driver.
!! At the moment, this driver serves as a mediator between the volume tracking
!! code with the desired api and the rest of truchas
!!
!! Peter Brady <ptb@lanl.gov>
!! 2017
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  This module provides procedures for driving the microstructure modeling
!!  kernel.  It serves as an adapter between top-level Truchas drivers (input,
!!  initialization, cycle, and output) and the microstructure modeling object.
!!  The existing top-level drivers manage no data; they merely provide for
!!  orchestrating the sequence of procedural steps.  Thus a major role of
!!  these subroutines is to assemble the data from various Truchas components
!!  that must be passed the the microstructure modeling object.
!!
!!  CALL READ_VOLUMETRACKING_NAMELIST (LUN) reads the first instance of the
!!    VOLUMETRACKING namelist from the file opened on logical unit LUN.  This
!!    is a collective procedure; the file is read on the I/O process and the
!!    data replicated to all other processes.  If an error occurs (I/O error
!!    or invalid data) a message is written and execution is halted.  The
!!    presence of this namelist serves to enable the volumetracker; it
!!    is permissible for there to be no instance of this namelist.
!!    NB: this should be called ONLY when flow is active.
!!
!!  CALL VTRACK_DRIVER_INIT (T) initializes the driver.  T is the initial time.
!!    This should be called only if heat transfer physics is enabled, and after
!!    its initialization.  If microstructure analysis has not been enabled this
!!    subroutine does nothing, and so may always be called.
!!
!!  CALL VTRACK_UPDATE (T) updates or advances the microstructure model to the
!!    next time level.  The subroutine handles collecting all the necessary
!!    mesh-based state arrays required by the model; only the current time T
!!    needs to be passed.
!!
!! NOTES
!!
!! The allocatation status of the private module variable THIS serves to
!! indicate whether the microstructure modeling kernel is enabled.  It is
!! allocated by READ_MICROSTRUCTURE_NAMELIST if the microstructure namelist
!! is present in the input file.
!!
#include "f90_assert.fpp"

module vtrack_driver
  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_mesh_type
  use volume_tracker_class
  use unsplit_volume_tracker_class  
  use parameter_list_type
  use truchas_logging_services
  use truchas_timers
  use index_partitioning
  use cell_tagged_mm_volumes_type, only : cell_tagged_mm_volumes
  implicit none
  private

  public :: vtrack_driver_init, vtrack_update, vtrack_driver_final
  public :: vtrack_enabled
  public :: vtrack_vof_view, vtrack_flux_vol_view, vtrack_unsplit_flux_vol_view
  public :: vtrack_liq_matid_view, vtrack_mat_band_view, vtrack_interface_band_view
  public :: get_vof_from_matl
  public :: vtrack_set_inflow_bc, vtrack_set_inflow_material
  public :: vtrack_velocity_overwrite
  private :: vtrack_update_mat_band
  public :: vtrack_open_interface_file, vtrack_close_interface_file, vtrack_write_interface

  integer, parameter, private :: band_map_width = 7 ! Size of band_map to create in +/- direction.

  !! Bundle up all the driver state data as a singleton THIS of private
  !! derived type.  All procedures use/modify this object.
  type :: vtrack_driver_data
    type(unstr_mesh), pointer :: mesh => null()  ! reference only -- do not own
    integer, allocatable :: liq_matid(:), sol_matid(:)
    class(volume_tracker), allocatable :: vt
    class(unsplit_volume_tracker), allocatable :: unsplit_vt
    real(r8), allocatable  :: vof(:,:) ! volume fractions for all materials
    real(r8), allocatable :: fvof_i(:,:) ! fluid/void volume fractions at start of update
    real(r8), allocatable :: fvof_o(:,:) ! fluid/void volume fractions at end of update
    real(r8), allocatable :: fvol_init(:) ! fluid/void volume fractions at start of simulation
    real(r8), allocatable :: flux_vol(:,:) ! flux volumes
    type(cell_tagged_mm_volumes), allocatable :: unsplit_flux_vol(:) ! Unsplit flux volumes
    real(r8), allocatable :: flux_vel(:) ! fluxing velocity - ragged form
    integer, allocatable :: mat_band(:,:) ! Signed material band for band_number = (material_id, cell_id)
    integer, allocatable :: mat_band_map(:,:) ! Unstructured mapping of band to material cell_id = (index, material)
    integer, allocatable :: xmat_band_map(:,:) ! Number of cells in mat_band_map for each band level.
    integer, allocatable :: interface_band(:) ! Signed band indicating if there is interface nearby
    integer, allocatable :: interface_band_map(:) ! Unstructued mapping of interface band to cell_id = (index)
    integer, allocatable :: xinterface_band_map(:) ! Number of cells in interface_band_map for each band level.    
    integer :: fluids ! number of fluids (not including void) to advance
    integer :: void ! 1 if simulation includes void else 0
    integer :: solid ! 1 if simulation includes solid else 0
    logical :: unsplit_advection
  end type vtrack_driver_data
  type(vtrack_driver_data), allocatable, target, save :: this
  
contains

  subroutine vtrack_driver_final
    if (allocated(this)) deallocate(this)
  end subroutine vtrack_driver_final

  logical function vtrack_enabled()
    vtrack_enabled = allocated(this)
  end function vtrack_enabled

  function vtrack_vof_view() result(p)
    real(r8), pointer :: p(:,:)
    ASSERT(vtrack_enabled())
    p => this%fvof_o
  end function vtrack_vof_view

  function vtrack_flux_vol_view() result(p)
    real(r8), pointer :: p(:,:)
    ASSERT(vtrack_enabled())
    p => this%flux_vol
  end function vtrack_flux_vol_view

  function vtrack_unsplit_flux_vol_view() result(p)
    type(cell_tagged_mm_volumes), pointer :: p(:)
    ASSERT(vtrack_enabled())
    p => this%unsplit_flux_vol
  end function vtrack_unsplit_flux_vol_view
  
  ! material ids corresponding to flux_vol columns
  function vtrack_liq_matid_view() result(p)
    integer, pointer :: p(:)
    p => this%liq_matid(:this%fluids)
  end function vtrack_liq_matid_view

  function vtrack_mat_band_view() result(p)
    integer, pointer :: p(:,:)
    ASSERT(vtrack_enabled())
    p => this%mat_band
  end function vtrack_mat_band_view

  function vtrack_interface_band_view() result(p)
    integer, pointer :: p(:)
    ASSERT(vtrack_enabled())
    p => this%interface_band
  end function vtrack_interface_band_view

  subroutine vtrack_open_interface_file()
    use truchas_phase_interface_output, only : TPIO_open
    if(vtrack_enabled()) then
      if(this%unsplit_advection) then
        call TPIO_open()
      end if
    end if
  end subroutine vtrack_open_interface_file
  
  subroutine vtrack_close_interface_file()
    use truchas_phase_interface_output, only : TPIO_close    
    if(vtrack_enabled()) then
      if(this%unsplit_advection) then
        call TPIO_close()
      end if
    end if
  end subroutine vtrack_close_interface_file

  subroutine vtrack_write_interface(t, dt, cycle_number)
    real(r8), intent(in) :: t
    real(r8), intent(in) :: dt
    integer, intent(in) :: cycle_number   

    if(vtrack_enabled()) then
      if(this%unsplit_advection) then
        call this%unsplit_vt%write_interface(t,dt,cycle_number)
      end if
    end if

  end subroutine vtrack_write_interface
  
  subroutine vtrack_driver_init(params)

    use mesh_manager, only: unstr_mesh_ptr
    use material_interop, only: void_material_index
    use property_data_module, only: isImmobile
    use parameter_module, only: nmat
    use geometric_volume_tracker_type
    use unsplit_geometric_volume_tracker_type
    use simple_volume_tracker_type
    use parallel_communication, only : global_sum

    type(parameter_list), intent(inout) :: params

    integer :: i,j,k
    logical :: track_interfaces

    allocate(this)

    call TLS_info('')
    call TLS_info('Configuring volume tracking ...')

    this%mesh => unstr_mesh_ptr('MAIN')
    INSIST(associated(this%mesh))

    this%solid = merge(1, 0, any(isImmobile(:nmat)))
    this%void = merge(1, 0, void_material_index > 0)
    this%fluids = count(.not.isImmobile(:nmat)) - this%void

    if (this%fluids == 0) then
      print *, 'no non-void fluids'
      return ! hmmm
    end if

    allocate(this%liq_matid(this%fluids+this%void))
    allocate(this%sol_matid(nmat-(this%fluids+this%void)))
    allocate(this%fvof_i(this%fluids+this%void+this%solid, this%mesh%ncell))
    allocate(this%fvof_o(this%fluids+this%void+this%solid, this%mesh%ncell))
    allocate(this%fvol_init(this%fluids+this%void+this%solid))    
    allocate(this%flux_vol(this%fluids+this%void,size(this%mesh%cface)))
    allocate(this%unsplit_flux_vol(size(this%mesh%cface)))    
    allocate(this%flux_vel(size(this%mesh%cface)))
    allocate(this%vof(nmat, this%mesh%ncell_onP))

    j = 1
    k = 1
    do i = 1, nmat
      if (i == void_material_index) then
        this%liq_matid(this%fluids+1) = i
      else if (.not.isImmobile(i)) then
        this%liq_matid(j) = i
        j = j + 1
      else
        this%sol_matid(k) = i
        k = k + 1
      end if
    end do

    ! Allocations for the band_map based on how near phases are
    allocate(this%mat_band(nmat, this%mesh%ncell))
    allocate(this%mat_band_map(this%mesh%ncell_onP, nmat))
    allocate(this%xmat_band_map(0:band_map_width+1, nmat))
    allocate(this%interface_band(this%mesh%ncell))
    allocate(this%interface_band_map(this%mesh%ncell_onP))
    allocate(this%xinterface_band_map(0:band_map_width+1))    

    call params%get('track_interfaces', track_interfaces)
    call params%get('unsplit_interface_advection', this%unsplit_advection, .true.)
    if (track_interfaces) then
      if(this%unsplit_advection) then
        allocate(unsplit_geometric_volume_tracker :: this%unsplit_vt)
      else
        allocate(geometric_volume_tracker :: this%vt)
      end if
    else
      allocate(simple_volume_tracker :: this%vt)
    end if

    call start_timer('Vof Initialization')

    call get_vof_from_matl(this%fvof_i)
    if(this%unsplit_advection) then
      call this%unsplit_vt%init(this%mesh, this%fluids, this%fluids+this%void, &
           this%fluids+this%void+this%solid, this%liq_matid, params)      
      ! Enforce VOF sums to 1
      do j = 1, this%mesh%ncell
        do k = 1, nmat
          if(this%fvof_i(k,j) < this%unsplit_vt%cutoff) then
            this%fvof_i(k,j) = 0.0_r8
          else if(this%fvof_i(k,j) > 1.0_r8 - this%unsplit_vt%cutoff) then
            this%fvof_i(k,j) = 1.0_r8
          end if
        end do
        this%fvof_i(:,j) = this%fvof_i(:,j) / sum(this%fvof_i(:,j))
      end do
      this%fvof_o = this%fvof_i
      call vtrack_update_mat_band()      
    else
      this%fvof_o = this%fvof_i      
      call this%vt%init(this%mesh, this%fluids, this%fluids+this%void, &
           this%fluids+this%void+this%solid, this%liq_matid, params)
    end if
    call stop_timer('Vof Initialization')

    ! Compute sum of volume of each phase in domain
    this%fvol_init = 0.0_r8
    do j = 1, this%mesh%ncell_onP
      this%fvol_init = this%fvol_init + this%fvof_o(:,j) * this%mesh%volume(j)
    end do
    do j = 1, nmat
      this%fvol_init(j) = global_sum(this%fvol_init(j))
    end do

    ! Open file for interface
    call vtrack_open_interface_file
    
    
  end subroutine vtrack_driver_init

  subroutine get_vof_from_matl(vof)

    use matl_utilities, only: matl_get_vof

    real(r8), intent(out) :: vof(:,:)

    integer :: i, n

    vof = 0.0_r8
    n = this%mesh%ncell_onP
    call matl_get_vof(this%vof)
    do i = 1, size(this%liq_matid)
      vof(i,:n) = this%vof(this%liq_matid(i),:)
    end do

    if (this%solid > 0) then
      do i = 1, n
        vof(this%fluids+this%void+this%solid,i) = sum(this%vof(this%sol_matid,i))
      end do
    end if

    call gather_boundary(this%mesh%cell_ip, vof)

  end subroutine get_vof_from_matl

  subroutine put_vof_into_matl()
    use matl_utilities, only: matl_set_vof
    integer :: i
    do i = 1, size(this%liq_matid)
      this%vof(this%liq_matid(i),:) = this%fvof_o(i,:this%mesh%ncell_onP)
    end do
    call matl_set_vof(this%vof)
  end subroutine put_vof_into_matl

  ! vel_fn is the outward oriented face-normal velocity
  ! vel_cc is the cell centered velocity
  subroutine vtrack_update(t, dt, vel_fn, vel_cc, initial)

    use parallel_communication, only : global_sum
    use parameter_module, only : string_len
    
    use constants_module
    real(r8), intent(in) :: t, dt
    real(r8), intent(in) :: vel_fn(:)
    real(r8), intent(in) :: vel_cc(:,:)
    logical, intent(in), optional :: initial

    integer :: i, j, k, f0, f1
    real(r8) :: vol_sum(size(this%fvof_i,1))

    real(r8) :: tmp
    character(string_len) :: message, myformat
    

    if (.not.allocated(this)) return
    if (this%fluids == 0) return

    call start_timer('Volumetracking')

    call get_vof_from_matl(this%fvof_i)

    if(this%unsplit_advection) then
      ! QUESTION : Where is correct place to put this? Do we ever stop changing VOF?
      ! Enforce VOF sums to 1
      do j = 1, this%mesh%ncell
        do k = 1, size(this%fvof_i,1)
          if(this%fvof_i(k,j) < this%unsplit_vt%cutoff) then
            this%fvof_i(k,j) = 0.0_r8
          else if(this%fvof_i(k,j) > 1.0_r8 - this%unsplit_vt%cutoff) then
            this%fvof_i(k,j) = 1.0_r8
          end if
        end do
        this%fvof_i(:,j) = this%fvof_i(:,j) / sum(this%fvof_i(:,j))
      end do 
      call vtrack_update_mat_band()      
      ! Unsplit advection written to use face velocities, operates on faces.
      ! Cell centered velocity also needed to form node velocities for projection.
      call this%unsplit_vt%flux_volumes(vel_fn, vel_cc, this%fvof_i, this%fvof_o, this%unsplit_flux_vol, &
           this%fluids, this%void, dt, this%interface_band)
    else
      ! Split advection works on a per-cell level.
      ! copy face velocities into cell-oriented array
      do i = 1,this%mesh%ncell
        f0 = this%mesh%xcface(i)
        f1 = this%mesh%xcface(i+1)-1
        do j = f0, f1
          k = this%mesh%cface(j)
          if (btest(this%mesh%cfpar(i),pos=1+j-f0)) then ! normal points inward
            this%flux_vel(j) = -vel_fn(k)
          else
            this%flux_vel(j) = vel_fn(k)
          end if          
        end do
      end do
    
      call this%vt%flux_volumes(this%flux_vel, vel_cc, this%fvof_i, this%fvof_o, this%flux_vol, &
          this%fluids, this%void, dt, this%interface_band)
    end if
   
    vol_sum = 0.0_r8
    do j = 1, this%mesh%ncell_onP
       vol_sum = vol_sum + this%fvof_o(:,j)*this%mesh%volume(j)
     end do
    do j = 1, size(vol_sum)
      vol_sum(j) = global_sum(vol_sum(j))
    end do
    write(myformat, '(a,i1,a)') '(a,',size(vol_sum),'es12.4)'
    write(message, trim(myformat)) 'Absolute volume change: ', vol_sum - this%fvol_init
    call TLS_info(message)
    write(message, trim(myformat)) 'Domain Normalized volume change: ',  (vol_sum - this%fvol_init)/sum(this%fvol_init)
    call TLS_info(message)
    write(message, trim(myformat)) 'Phase Normalized volume change: ',  (vol_sum - this%fvol_init)/this%fvol_init
    call TLS_info(message)    

    ! update the matl structure if this isn't the initial pass
    if (present(initial)) then
      if (.not.initial) &
          call put_vof_into_matl()
    else
      call put_vof_into_matl()
    end if

    call stop_timer('Volumetracking')

  end subroutine vtrack_update

  !! Sets the inflow material for the given list of boundary faces. MATID is
  !! the Truchas material number which must be a fluid. MATID can also be the
  !! special value 0 which indicates that materials should be fluxed into the
  !! adjacent cell in proportion to the material fractions currently present
  !! in the cell. The initial default for all boundary faces is the latter.

  subroutine vtrack_set_inflow_material(matid, faces)

    integer, intent(in) :: matid, faces(:)
    integer :: n

    if (matid == 0) then
      n = 0
    else
      do n = 1, size(this%liq_matid)
        if (this%liq_matid(n) == matid) exit
      end do
      ASSERT(n <= size(this%liq_matid))
    end if

    call this%vt%set_inflow_material(n, faces)

  end subroutine vtrack_set_inflow_material

  !! Sets the inflow material boundary conditions specified by the given
  !! parameter list. The structure of the list is a list of sublists. Only
  !! those that define an "inflow-material" parameter are significant. Its
  !! value and the corresponding value of the required "face-set-ids"
  !! parameter are used to set the inflow material on a portion of the
  !! boundary. Other sublists are ignored.
  subroutine vtrack_set_inflow_bc(params, stat, errmsg)

    use bndry_face_group_builder_type

    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list_iterator) :: piter
    type(parameter_list), pointer :: plist
    type(bndry_face_group_builder) :: builder
    integer, allocatable :: setids(:), mlist(:), xgroup(:), index(:)
    integer :: j, n, matid, ngroup

    call builder%init(this%mesh)

    piter = parameter_list_iterator(params, sublists_only=.true.)
    n = piter%count()
    allocate(mlist(n))
    n = 0
    do while (.not.piter%at_end())
      plist => piter%sublist()
      if (plist%is_parameter('inflow-material')) then
        call plist%get('inflow-material', matid, stat=stat, errmsg=errmsg)
        if (stat /= 0) return
        call plist%get('face-set-ids', setids, stat=stat, errmsg=errmsg)
        if (stat /= 0) return
        call builder%add_face_group(setids, stat, errmsg)
        if (stat /= 0) return
        n = n + 1
        mlist(n) = matid
      end if
      call piter%next
    end do

    call builder%get_face_groups(n, xgroup, index, omit_offp=.true.)

    do j = 1, n
      call vtrack_set_inflow_material(mlist(j), index(xgroup(j):xgroup(j+1)-1))
    end do

  end subroutine vtrack_set_inflow_bc

  subroutine vtrack_velocity_overwrite(a_time, a_face_velocity, a_cell_velocity)

    use vof_velocity_overwrite, only : velocity_overwrite_requested, &
                                       velocity_overwrite_case
    use physics_module, only : prescribed_flow
    use constants_module, only : pi

    real(r8), intent(in) :: a_time
    real(r8), intent(inout) :: a_face_velocity(:)
    real(r8), intent(inout) :: a_cell_velocity(:,:)

    integer :: j, f, n, e
    real(r8) :: edge(3,2)
    real(r8) :: node_location(3), full_face_velocity(3), vel(3)
    real(r8) :: projected_point_x, projected_point_y, rotation
    real(r8) :: time_multiplier

    if (.not.allocated(this)) return
    if(.not.velocity_overwrite_requested) return
    ASSERT(prescribed_flow)

    select case(trim(velocity_overwrite_case))
    case('Rotation2D')
       do f = 1,this%mesh%nface
         node_location = this%mesh%face_centroid(:,f)
         full_face_velocity(1) = -2.0_r8*pi*node_location(2)
         
         full_face_velocity(2) = 2.0_r8*pi*node_location(1)
         
         full_face_velocity(3) = 0.0_r8
         
         a_face_velocity(f) = dot_product(full_face_velocity, this%mesh%normal(:,f)/this%mesh%area(f))
         
       end do
       
       do j = 1,this%mesh%ncell
         node_location = this%mesh%cell_centroid(:,j)
         a_cell_velocity(1,j) = -2.0_r8*pi*node_location(2)
         
         a_cell_velocity(2,j) = 2.0_r8*pi*node_location(1)
         
         a_cell_velocity(3,j) = 0.0_r8
         
       end do
              
    case('Rotation3D')
       rotation = -0.25_r8*pi
       do f = 1,this%mesh%nface
         node_location = this%mesh%face_centroid(:,f)
         projected_point_x = dot_product(node_location, 1.0_r8/sqrt(2.0_r8)*[1.0_r8, 0.0_r8, 1.0_r8])
         projected_point_y = dot_product(node_location, [0.0_r8, 1.0_r8, 0.0_r8])         
         vel(1) = -2.0_r8*pi*projected_point_y         
         vel(2) = 2.0_r8*pi*projected_point_x
         vel(3) = 0.0_r8

         ! Rotate about y-axis
         full_face_velocity(1) = dot_product([cos(rotation), 0.0_r8, sin(rotation)], vel)
         full_face_velocity(2) = dot_product([0.0_r8, 1.0_r8, 0.0_r8], vel)
         full_face_velocity(3) = dot_product([-sin(rotation), 0.0_r8, cos(rotation)], vel)         
         a_face_velocity(f) = dot_product(full_face_velocity, this%mesh%normal(:,f)/this%mesh%area(f))
       end do

       do j = 1,this%mesh%ncell
         node_location = this%mesh%cell_centroid(:,j)
         projected_point_x = dot_product(node_location, 1.0_r8/sqrt(2.0_r8)*[1.0_r8, 0.0_r8, 1.0_r8])
         projected_point_y = dot_product(node_location, [0.0_r8, 1.0_r8, 0.0_r8])

         vel(1) = -2.0_r8*pi*projected_point_y         
         vel(2) = 2.0_r8*pi*projected_point_x
         vel(3) = 0.0_r8

         ! Rotate about y-axis
         a_cell_velocity(1,j) = dot_product([cos(rotation), 0.0_r8, sin(rotation)], vel)
         a_cell_velocity(2,j) = dot_product([0.0_r8, 1.0_r8, 0.0_r8], vel)
         a_cell_velocity(3,j) = dot_product([-sin(rotation), 0.0_r8, cos(rotation)], vel)         
         
       end do
       
    case('Deformation2D')
       ! Generally on -0.5,0.5 mesh, so add 0.5 to node position
       do f = 1,this%mesh%nface
          node_location = this%mesh%face_centroid(:,f)+[0.5_r8,0.5_r8,0.0_r8]
          full_face_velocity(1) = -2.0_r8*sin(pi*node_location(1))**2 &
                                 * sin(pi*node_location(2)) &
                                 * cos(pi*node_location(2)) &
                                 * cos(pi*a_time / 8.0_r8)
          
          full_face_velocity(2) = 2.0_r8*sin(pi*node_location(2))**2 &
                                 * sin(pi*node_location(1)) &
                                 * cos(pi*node_location(1)) &
                                 * cos(pi*a_time / 8.0_r8)          
          
          full_face_velocity(3) = 0.0_r8

          a_face_velocity(f) = dot_product(full_face_velocity, this%mesh%normal(:,f)/this%mesh%area(f))

       end do
       
       do j = 1,this%mesh%ncell
          node_location = this%mesh%cell_centroid(:,j)+[0.5_r8,0.5_r8,0.0_r8]
          a_cell_velocity(1,j) = -2.0_r8*sin(pi*node_location(1))**2 &
                                 * sin(pi*node_location(2)) &
                                 * cos(pi*node_location(2)) &
                                 * cos(pi*a_time / 8.0_r8)
          
          a_cell_velocity(2,j) = 2.0_r8*sin(pi*node_location(2))**2 &
                                 * sin(pi*node_location(1)) &
                                 * cos(pi*node_location(1)) &
                                 * cos(pi*a_time / 8.0_r8)
          
          a_cell_velocity(3,j) = 0.0_r8
          
       end do

    case('Deformation3D')
       ! Generally on -0.5,0.5 mesh, but should be on [0.1], so add 0.5 to node position
       time_multiplier = cos(pi*a_time/3.0_r8)
       do f = 1,this%mesh%nface

         a_face_velocity(f) = 0.0_r8
         associate(nodes => this%mesh%fnode(this%mesh%xfnode(f):this%mesh%xfnode(f+1)-1))           
           do n = 1, size(nodes) 
             edge(:,1) = this%mesh%x(:,nodes(n))+[0.5_r8, 0.5_r8, 0.5_r8]
             edge(:,2) = this%mesh%x(:,nodes(mod(n, size(nodes))+1))+[0.5_r8, 0.5_r8, 0.5_r8]
             vel = 0.0_r8
             do e = 1,2
               vel = vel + [0.25_r8/pi*(2.0_r8*sin(2.0_r8*pi*edge(1,e))*sin(pi*edge(2,e))**2*cos(2.0_r8*pi*edge(3,e)) &
                    -sin(2.0_r8*pi*edge(1,e))*cos(2.0_r8*pi*edge(2,e))) ,&
                    sin(pi*edge(1,e))**2*sin(2.0_r8*pi*edge(2,e))*cos(2.0_r8*pi*edge(3,e))/pi ,&
                    0.0_r8]               
             end do
             vel = 0.5_r8*vel
             a_face_velocity(f) = a_face_velocity(f) + dot_product(vel, edge(:,2)-edge(:,1))
           end do
         end associate
         a_face_velocity(f) = a_face_velocity(f) / this%mesh%area(f) * time_multiplier

       end do
       
       do j = 1,this%mesh%ncell
          node_location = this%mesh%cell_centroid(:,j)+[0.5_r8,0.5_r8,0.5_r8]
          a_cell_velocity(1,j) = 2.0_r8*sin(pi*node_location(1))**2 &
                                 * sin(2.0_r8*pi*node_location(2)) &
                                 * sin(2.0_r8*pi*node_location(3)) &
                                 * cos(pi*a_time / 3.0_r8)
          
          a_cell_velocity(2,j) =  -sin(2.0_r8*pi*node_location(1)) &
                                 * sin(pi*node_location(2))**2 &
                                 * sin(2.0_r8*pi*node_location(3)) &
                                 * cos(pi*a_time / 3.0_r8)
          
          a_cell_velocity(3,j) =  -sin(2.0_r8*pi*node_location(1)) &
                                 * sin(2.0_r8*pi*node_location(2)) &
                                 * sin(pi*node_location(3))**2 &
                                 * cos(pi*a_time / 3.0_r8)        
          
       end do


    case default
       call TLS_fatal('Unknown case provided for overwriting velocity to prescribed field')
    end select   

  end subroutine vtrack_velocity_overwrite

  ! NOTE: There is a lot of parallel communication in this routine.
  ! It could probably be rewritten to use less. Also, if a
  ! band is not required on a per-material basis, could save
  ! a factor of nmat communications by just constructing a
  ! "is there interface of any kind or not" band.
  subroutine vtrack_update_mat_band()    

    integer :: b, j, n, k, c
    integer :: total_phases
    logical :: has_interface
    integer :: smallest_band

    total_phases = this%fluids + this%void + this%solid

    ! QUESTION : How to incorporate knowledge of BC into band?
    ! If I have a full gas domain, and am injecting liquid, need to
    ! recognize that and make the BC internal to domain band(0)

    ! Seed the initial 0-band that has interface
    this%xmat_band_map(0,:) = 1
    this%xmat_band_map(1,:) = this%xmat_band_map(0,:)
    do j = 1,this%mesh%ncell_onP
      do k = 1,total_phases
        has_interface = (this%fvof_i(k,j) >  this%unsplit_vt%cutoff) .and. &
             (this%fvof_i(k,j) < 1.0_r8 - this%unsplit_vt%cutoff)
        if(.not. has_interface) then
           ! Need to also check condition where cell is full/empty and neighbor empty/full
           if(this%fvof_i(k,j) < this%unsplit_vt%cutoff) then ! Cell is empty
              ! Check if any neighbors are full
              associate(cneigh => this%mesh%cnhbr(this%mesh%xcnhbr(j):this%mesh%xcnhbr(j+1)-1))
                do n = 1, size(cneigh)
                   if(cneigh(n) /=0 ) then
                      if(this%fvof_i(k,cneigh(n)) > 1.0_r8 -  this%unsplit_vt%cutoff) then
                         has_interface = .true.
                         exit
                      end if                      
                   end if
                end do

              end associate              
           else
              ! Check if any neighbors are empty
              associate(cneigh => this%mesh%cnhbr(this%mesh%xcnhbr(j):this%mesh%xcnhbr(j+1)-1))
                do n = 1, size(cneigh)
                   if(cneigh(n) /=0 ) then
                      if(this%fvof_i(k,cneigh(n)) <  this%unsplit_vt%cutoff) then
                         has_interface = .true.
                         exit
                      end if                      
                   end if
                end do

              end associate              
           end if           
        end if        
        if(has_interface) then
           this%mat_band(k,j) = 0
           this%mat_band_map(this%xmat_band_map(1,k),k) = j            
           this%xmat_band_map(1,k) = this%xmat_band_map(1,k) + 1
        else
           ! Initialize band positive inside phase, negative outside.
           this%mat_band(k,j) = int(sign(real(band_map_width+1,r8), this%fvof_i(k,j)-0.5_r8))
        end if
      end do
    end do

    call gather_boundary(this%mesh%cell_ip, this%mat_band)

! print*,'ADVECTING EVERYWHERE RIGHT NOW - BAND TURNED OFF'
!this%mat_band(:,:) = 0
    ! Now loop through and fill out band to +/- band_map_width
    ! To fill in band b, loop though cells in b-1
    do b = 1, band_map_width
      do k = 1, total_phases
        this%xmat_band_map(b+1,k) = this%xmat_band_map(b,k)
        associate( cell_ind => this%mat_band_map(this%xmat_band_map(b-1,k):this%xmat_band_map(b,k)-1,k))
          do c = 1,size(cell_ind)
            j = cell_ind(c)   

             associate(cneigh => this%mesh%cnhbr(this%mesh%xcnhbr(j):this%mesh%xcnhbr(j+1)-1))
               do n = 1, size(cneigh)
                 if(cneigh(n) /= 0 .and. cneigh(n) <= this%mesh%ncell_onP) then
                   if(abs(this%mat_band(k,cneigh(n))) .eq. band_map_width+1) then
                     ! Not yet set
                     this%mat_band(k,cneigh(n)) = sign(b, this%mat_band(k,cneigh(n)))                     
                     this%mat_band_map(this%xmat_band_map(b+1,k),k) = cneigh(n)
                     this%xmat_band_map(b+1,k) =  this%xmat_band_map(b+1,k) + 1                     
                   end if
                 end if
               end do
           end associate          

          end do
        end associate

      end do

      call gather_boundary(this%mesh%cell_ip, this%mat_band)
      ! Also loop over halo to catch band updated on boundary
      do k = 1, total_phases
        do j = this%mesh%ncell_onP+1, this%mesh%ncell
          if(abs(this%mat_band(k,j)) == b-1) then
             associate(cneigh => this%mesh%cnhbr(this%mesh%xcnhbr(j):this%mesh%xcnhbr(j+1)-1))
               do n = 1, size(cneigh)
                 if(cneigh(n) /= 0 .and. cneigh(n) <= this%mesh%ncell_onP) then
                   if(abs(this%mat_band(k,cneigh(n))) == band_map_width+1) then
                     ! Not yet set
                     this%mat_band(k,cneigh(n)) = sign(b, this%mat_band(k,cneigh(n)))
                     this%mat_band_map(this%xmat_band_map(b+1,k),k) = cneigh(n)
                     this%xmat_band_map(b+1,k) =  this%xmat_band_map(b+1,k) + 1
                   end if
                 end if
               end do
             end associate
           end if
         end do
       end do
      
     end do

     ! Create interface band map that indicates presence of any interface
     ! Only has positive (or 0) band entries.
     this%interface_band = band_map_width+1
     this%xinterface_band_map(0) = 1
     do b = 0, band_map_width
         this%xinterface_band_map(b+1) = this%xinterface_band_map(b)
       do j = 1, this%mesh%ncell
         smallest_band = minval(abs(this%mat_band(:,j)))
         if(smallest_band == b) then
           ! Add to band
           this%interface_band(j) = smallest_band
           if(j <= this%mesh%ncell_onP) then
             this%interface_band_map(this%xinterface_band_map(b+1)) = j
             this%xinterface_band_map(b+1) =  this%xinterface_band_map(b+1) + 1
           end if
         end if
       end do
     end do
     

 end subroutine vtrack_update_mat_band

end module vtrack_driver
