#include "f90_assert.fpp"

module volume_tracker_type

  use kinds, only: r8
  use truchas_logging_services
  use truchas_timers
  use unstr_mesh_type
  use index_partitioning
  implicit none
  private

  public :: volume_tracker

  type :: volume_tracker
    type(unstr_mesh), pointer :: mesh ! unowned reference
    integer :: location_iter_max ! maximum number of iterations to use in fitting interface
    integer :: subcycles
    logical :: use_brents_method, nested_dissection
    real(r8) :: location_tol ! tolerance of plic fit
    real(r8) :: cutoff ! allow volume fraction {0,(cutoff,1]}
    real(r8), allocatable :: flux_vol_sub(:,:), normal(:,:,:)
    ! node/face/cell workspace
    real(r8), allocatable :: w_node(:,:), w_face(:,:)
    integer, allocatable :: priority(:)
    integer :: nrealfluid, nfluid, nmat ! # of non-void fluids, # of fluids incl. void, # of materials
  contains
    procedure :: init
    procedure :: flux_volumes
    procedure :: normals
    procedure :: donor_fluxes_nd
    procedure :: donor_fluxes_os
    procedure :: flux_renorm
    procedure :: flux_acceptor
    !procedure :: flux_bc
    procedure :: accumulate_volume
    procedure :: enforce_bounded_vof
  end type volume_tracker

contains

  subroutine init(this, mesh, nrealfluid, nfluid, nmat, liq_matid, params)

    use parameter_list_type
    use property_module, only: get_truchas_material_id
    use f08_intrinsics, only: findloc

    class(volume_tracker), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    integer, intent(in) :: nrealfluid, nfluid, nmat, liq_matid(:)
    type(parameter_list), intent(inout) :: params
    integer :: i

    this%mesh => mesh
    this%nrealfluid = nrealfluid
    this%nfluid = nfluid
    this%nmat = nmat

    call params%get('location_iter_max', this%location_iter_max, default=20)
    call params%get('location_tol', this%location_tol, default=1.0e-8_r8)
    call params%get('cutoff', this%cutoff, default=1.0e-8_r8)
    call params%get('subcycles', this%subcycles, default=2)
    call params%get('use_brents_method', this%use_brents_method, default=.true.)
    call params%get('nested_dissection', this%nested_dissection, default=.true.)

    ! convert user material ids to array index
    if (params%is_vector('material_priority')) then
      call params%get('material_priority', this%priority)
      do i = 1,size(this%priority)
        if (this%priority(i) < 1) cycle ! solid (-1) is handled later
        this%priority(i) = findloc(liq_matid, get_truchas_material_id(this%priority(i)))

        ! make sure we found a liquid material
        ! TODO: need better error handling here
        INSIST(this%priority(i) > 0)
      end do

      ! The current expectation is that a user will
      ! use a material number of -1 to indicate solid.
      ! Internally, if solid is present it is the last material
      if (this%nmat > this%nfluid) then
        where (this%priority == -1) this%priority = this%nmat
      end if
    else
      this%priority = [(i, i=1,this%nmat)]
    end if
    ASSERT(size(this%priority) == this%nmat)
    ASSERT(all(this%priority > 0) .and. all(this%priority <= this%nmat))
    ! TODO: assert that each material appears exactly once

    allocate(this%flux_vol_sub(this%nfluid,size(mesh%cface)))
    allocate(this%normal(3,this%nmat,mesh%ncell))
    allocate(this%w_node(2,mesh%nnode))
    allocate(this%w_face(this%nfluid,mesh%nface))

  end subroutine init


  ! flux volumes routine assuming vel/flux_vol is a cface-like array
  subroutine flux_volumes(this, vel, vof_n, vof, flux_vol, fluids, void, dt)
    class(volume_tracker), intent(inout) :: this
    real(r8), intent(in) :: vel(:), vof_n(:,:), dt
    real(r8), intent(out) :: flux_vol(:,:), vof(:,:)
    integer, intent(in) :: fluids, void

    integer :: i,ii,j,k,f0,f1,kk
    real(r8) :: sub_dt, disp(3), cent(3)

    flux_vol = 0.0_r8
    sub_dt = dt/real(this%subcycles, r8)
    vof = vof_n

    do i = 1, this%subcycles
      call this%normals(vof)

      if (this%nested_dissection) then
        call this%donor_fluxes_nd(vel, vof, sub_dt)
      else
        call this%donor_fluxes_os(vel, vof, sub_dt)
      end if


      call this%flux_renorm(vel, vof_n, flux_vol, sub_dt)

      call this%flux_acceptor()

      !call this%flux_bc(fluxing_velocity, vof_n, volume_flux_sub)

      call this%accumulate_volume(vof, flux_vol)

      call this%enforce_bounded_vof(vof, flux_vol, fluids, void)
    end do

  end subroutine flux_volumes


  subroutine normals(this, vof)

    use flow_operators, only: gradient_cc
    use f08_intrinsics, only: findloc
    intrinsic :: norm2

    class(volume_tracker), intent(inout) :: this
    real(r8), intent(in)  :: vof(:,:)

    real(r8) :: mag
    integer :: i,j,k,c
    logical :: hasvof(size(vof,dim=1))

    call start_timer('normals')
    do i = 1, this%nmat
      call gradient_cc(this%normal(1,i,:), this%normal(2,i,:), this%normal(3,i,:), &
          vof(i,:), this%w_node(1,:), this%w_node(2,:))
    end do

    do i = 1, this%mesh%ncell_onP
      hasvof = vof(:,i) > 0.0_r8
      c = count(hasvof)
!!$      if (c < 2) then
!!$        this%normal(:,:,i) = 1.0_r8
!!$        cycle
!!$      end if

      this%normal(:,:,i) = -this%normal(:,:,i)
      ! enforce consistency for two materials
      if (c == 2) then
        j = findloc(hasvof,.true.)
        k = findloc(hasvof,.true.,back=.true.)
        this%normal(:,k,i) = -this%normal(:,j,i)
      endif

      ! perform onion skin if requested
      ! normal vectors (gradients) include previous materials in the priority ordering
      if (.not.this%nested_dissection .and. c > 2) then
        do j = 2, this%nmat
          k = this%priority(j)
          this%normal(:,k,i) = this%normal(:,k,i) + this%normal(:,this%priority(j-1),i)
        end do
      end if

      ! normalize and remove smallish components due to robustness issues in nested disection
      do j = 1 , this%nmat
        !if (vof(j,i) <= this%cutoff) cycle ! should this be 0 or cutoff?
        ! remove small values
        do k = 1, 3
          if (abs(this%normal(k,j,i)) < epsilon(1.0_r8)) this%normal(k,j,i) = 0.0_r8
        end do
        ! normalize if possible
        mag = norm2(this%normal(:,j,i))
        if (mag > epsilon(1.0_r8)) this%normal(:,j,i) = this%normal(:,j,i)/mag
        ! remove slightly larger values
        do k = 1, 3
          if (abs(this%normal(k,j,i)) < 1.0e-6_r8) this%normal(k,j,i) = 0.0_r8
        end do
        ! normalize if possible
        mag = norm2(this%normal(:,j,i))
        if (mag > epsilon(1.0_r8)) then
          this%normal(:,j,i) = this%normal(:,j,i)/mag
        else
          this%normal(:,j,i) = 1.0_r8
        end if
      end do
    end do
    ! will need normals for vof reconstruction in ghost cells
    call gather_boundary(this%mesh%cell_ip, this%normal)

    call stop_timer('normals')

  end subroutine normals


  subroutine donor_fluxes_os(this, vel, vof, dt)

    use cell_topology
    use hex_types, only: cell_data

    class(volume_tracker), intent(inout) :: this
    real(r8), intent(in)  :: dt, vof(:,:), vel(:)

    real(r8) :: face_normal(3,6), flux_vols(size(vof,dim=1),6)
    integer :: i,j,k,ierr, face_vid(4,6)
    type(cell_data) :: cell

        ! calculate the flux volumes for each face
    call start_timer('reconstruct/advect')

    do i = 1, this%mesh%ncell
      associate (cn => this%mesh%cnode(this%mesh%xcnode(i):this%mesh%xcnode(i+1)-1), &
          fi => this%mesh%cface(this%mesh%xcface(i):this%mesh%xcface(i+1)-1))

        do j = 1, size(fi)
          k = fi(j)
          if (btest(this%mesh%cfpar(i),pos=j)) then
            face_normal(:,j) = -this%mesh%normal(:,k)/this%mesh%area(k)
          else
            face_normal(:,j) = this%mesh%normal(:,k)/this%mesh%area(k)
          end if
        end do

        select case (size(cn))
!!$        case (4)
!!$          call cell%init(this%mesh%x(:,cn), reshape(source=TET4_FACES,shape=[3,4]), &
!!$              TET4_EDGES, this%mesh%volume(i), face_normal(:,1:size(fi)))
!!$
!!$        case (5)
!!$          ! zero treated as sentinel value in multimat_cell procedures
!!$          face_vid = 0
!!$          do j = 1, size(fi)
!!$            face_vid(1:PYR5_FSIZE(j),j) = PYR5_FACES(PYR5_XFACE(j):PYR5_XFACE(j+1)-1)
!!$          end do
!!$          call cell%init(this%mesh%x(:,cn), face_vid, &
!!$              PYR5_EDGES, this%mesh%volume(i), face_normal(:,1:size(fi)))
!!$
!!$        case (6)
!!$          ! zero treated as sentinel value in multimat_cell procedures
!!$          face_vid = 0
!!$          do j = 1, size(fi)
!!$            face_vid(1:WED6_FSIZE(j),j) = WED6_FACES(WED6_XFACE(j):WED6_XFACE(j+1)-1)
!!$          end do
!!$          call cell%init(this%mesh%x(:,cn), face_vid, &
!!$              WED6_EDGES, this%mesh%volume(i), face_normal(:,1:size(fi)))

        case (8)
          call cell%init(this%mesh%x(:,cn), this%mesh%volume(i), this%mesh%area(fi), &
              face_normal(:,1:size(fi)))

        case default
          call TLS_fatal('unaccounted topology in donor_fluxes_nd')
        end select

!!$        if (ierr /= 0) call TLS_fatal('cell_outward_volflux failed: could not initialize cell')

!!$        call cell%partition(vof(:,i), this%normal(:,:,i), this%cutoff, this%location_tol, &
!!$            this%location_iter_max)

        this%flux_vol_sub(:,this%mesh%xcface(i):this%mesh%xcface(i+1)-1) = &
            cell_volume_flux(dt, cell, vof(:,i), this%normal(:,:,i), &
            vel(this%mesh%xcface(i):this%mesh%xcface(i+1)-1), this%cutoff, &
            this%priority, this%nfluid, this%nmat)
!!$        if (ierr /= 0) call TLS_fatal('cell_outward_volflux failed')
      end associate
    end do

    call stop_timer('reconstruct/advect')

  end subroutine donor_fluxes_os

    ! get the volume flux for every material in the given cell
  function cell_volume_flux (dt, cell, vof, int_norm, vel, cutoff, priority, nfluid, nmat)
    use hex_types,           only: cell_data
    use locate_plane_modulez, only: locate_plane_hex
    use flux_volume_modulez,  only: flux_vol_quantity
    use array_utils,         only: last_true_loc

    real(r8), intent(in) :: dt, int_norm(:,:), vof(:), vel(:), cutoff
    integer, intent(in) :: priority(:), nfluid, nmat
    type(cell_data), intent(in) :: cell
    real(r8) :: cell_volume_flux(nfluid,6)

    real(r8) :: Vofint, vp, dvol
    real(r8) :: flux_vol_sum(6)
    integer :: ni,f,locate_plane_niters,nlast, nint, ierr,nmat_in_cell
    logical :: is_mixed_donor_cell
    type(locate_plane_hex) :: plane_cell
    type(flux_vol_quantity) :: Flux_Vol

    !call start_timer ("ra_cell")

    cell_volume_flux = 0.0_r8
    flux_vol_sum = 0.0_r8
    nmat_in_cell = count(vof > 0.0_r8)
    !nint = count(vof > 0.0_r8)
    ! Here, I am not certain the conversion from pri_ptr to direct material indices worked properly.
    ! This will be clear when trying 3 or more materials. -zjibben

    ! Loop over the interfaces in priority order
    do ni = 1,nmat-1
      ! check if this is a mixed material cell
      ! First accumulate the volume fraction of this material and materials with lower priorities.
      ! Force 0.0 <= Vofint <= 1.0
      Vofint = min(max(sum(vof(priority(:ni))), 0.0_r8), 1.0_r8)
      is_mixed_donor_cell = cutoff < Vofint .and. Vofint < (1.0_r8 - cutoff)

      ! locate each interface plane by computing the plane constant
      if (is_mixed_donor_cell) then
        call plane_cell%init (int_norm(:,priority(ni)), vofint, cell%volume, cell%node)
        call plane_cell%locate_plane (locate_plane_niters)
      end if

      ! calculate delta advection volumes for this material at each donor face and accumulate the sum
      cell_volume_flux(priority(ni),:) =  material_volume_flux (flux_vol_sum, plane_cell, cell, &
           is_mixed_donor_cell, vel, dt, vof(priority(ni)), cutoff)
    end do ! interface loop

    ! get the id of the last material
    nlast = priority(last_true_loc(vof(priority) >= cutoff))

    !call start_timer ("last_loop")
    ! Compute the advection volume for the last material.
    do f = 1,6
      ! Recalculate the total flux volume for this face.
      Flux_Vol%Vol = dt*vel(f)*cell%face_area(f)
      if (abs(flux_vol%vol) > 0.5_r8 * cell%volume) then
        write(*,*) dt,flux_vol%vol,cell%volume,flux_vol%vol/cell%volume
        call TLS_fatal('advection timestep too large')
      end if
      if (Flux_Vol%Vol <= cutoff*cell%volume) cycle

      ! For donor cells containing only one material, assign the total flux.
      if (nmat_in_cell==1) then
        cell_volume_flux(nlast,f) = Flux_Vol%Vol
      else
        ! The volume flux of the last material shouldn't be less than
        ! zero nor greater than the volume of this material in the donor cell.
        dvol = min(max(abs(Flux_Vol%Vol - Flux_Vol_Sum(f)), 0.0_r8), Vof(nlast)*cell%volume)

        ! Store the last material's volume flux.
        if (dvol > cutoff*cell%volume) cell_volume_flux(nlast,f) = dvol
      end if
    end do ! face loop
    !call stop_timer ("last_loop")

    !call stop_timer ("ra_cell")

  end function cell_volume_flux

  ! calculate the flux of one material in a cell
  function material_volume_flux (flux_vol_sum, plane_cell, cell, is_mixed_donor_cell, vel, dt, vof, cutoff)
    use hex_types,              only: cell_data
    use truncate_volume_modulez, only: truncate_volume, face_param, truncvol_data
    use flux_volume_modulez,     only: flux_vol_quantity, flux_vol_vertices
    use locate_plane_modulez, only: locate_plane_hex

    real(r8),               intent(inout) :: flux_vol_sum(:)
    real(r8),               intent(in)    :: vel(:), dt, vof, cutoff
    type(locate_plane_hex), intent(in)    :: plane_cell
    type(cell_data),        intent(in)    :: cell
    logical,                intent(in)    :: is_mixed_donor_cell
    real(r8)                              :: material_volume_flux(6)

    integer                 :: f,ff, idbg
    real(r8)                :: vp, dvol
    type(truncvol_data)     :: trunc_vol(6)
    type(flux_vol_quantity) :: Flux_Vol

    !call start_timer ("material_flux_vol")

    material_volume_flux = 0.0_r8

    do f = 1,6
      ! Flux volumes
      Flux_Vol%Fc  = f
      Flux_Vol%Vol = dt*vel(f)*cell%face_area(f)
      if (Flux_Vol%Vol <= cutoff*cell%volume) cycle

      ! calculate the vertices describing the volume being truncated through the face
      call flux_vol_vertices (f, cell, is_mixed_donor_cell, vel(f)*dt, Flux_Vol, cutoff)

      if (is_mixed_donor_cell) then
        ! Now compute the volume truncated by interface planes in each flux volumes.
        do ff = 1,6
          trunc_vol(ff) = face_param (plane_cell, 'flux_cell', ff, flux_vol%xv)
        end do

        ! For mixed donor cells, the face flux is in Int_Flux%Advection_Volume.
        Vp = truncate_volume(plane_cell, trunc_vol)
      else
        ! For clean donor cells, the entire Flux volume goes to the single donor material.
        if (vof >= (1.0_r8-cutoff)) then
          Vp = abs(flux_vol%vol)
        else
          Vp = 0.0_r8
        end if
      end if

      ! If Vp is close to 0 set it to 0.  If it is close
      ! to 1 set it to 1. This will avoid numerical round-off.
      if (Vp > (1.0_r8-cutoff)*abs(Flux_Vol%Vol)) Vp = abs(Flux_Vol%Vol)

      ! Make sure that the current material-integrated advection
      ! volume hasn't decreased from its previous value.  (This
      ! can happen if the interface significantly changed its
      ! orientation and now crosses previous interfaces.)  Also
      ! limit Volume_Flux_Sub to take no more than the material
      ! occupied volume in the donor cell.
      dvol = min(max(Vp - Flux_Vol_Sum(f), 0.0_r8), Vof*cell%volume)

      ! Now gather advected volume information into Volume_Flux_Sub
      material_volume_flux(f) = dvol
      Flux_Vol_Sum(f) = Flux_Vol_Sum(f) + dvol
    end do ! face loop

    !call stop_timer ("material_flux_vol")

  end function material_volume_flux


  subroutine donor_fluxes_nd(this, vel, vof, dt)

    use cell_topology
    use multimat_cell_type

    class(volume_tracker), intent(inout) :: this
    real(r8), intent(in)  :: dt, vof(:,:), vel(:)

    real(r8) :: face_normal(3,6), flux_vols(size(vof,dim=1),6)
    integer :: i,j,k,ierr, face_vid(4,6)
    type(multimat_cell) :: cell

    ! calculate the flux volumes for each face
    call start_timer('reconstruct/advect')

    do i = 1, this%mesh%ncell
      associate (cn => this%mesh%cnode(this%mesh%xcnode(i):this%mesh%xcnode(i+1)-1), &
          fi => this%mesh%cface(this%mesh%xcface(i):this%mesh%xcface(i+1)-1))

        do j = 1, size(fi)
          k = fi(j)
          if (btest(this%mesh%cfpar(i),pos=j)) then
            face_normal(:,j) = -this%mesh%normal(:,k)/this%mesh%area(k)
          else
            face_normal(:,j) = this%mesh%normal(:,k)/this%mesh%area(k)
          end if
        end do

        select case (size(cn))
        case (4)
          call cell%init(ierr, this%mesh%x(:,cn), reshape(source=TET4_FACES,shape=[3,4]), &
              TET4_EDGES, this%mesh%volume(i), face_normal(:,1:size(fi)))

        case (5)
          ! zero treated as sentinel value in multimat_cell procedures
          face_vid = 0
          do j = 1, size(fi)
            face_vid(1:PYR5_FSIZE(j),j) = PYR5_FACES(PYR5_XFACE(j):PYR5_XFACE(j+1)-1)
          end do
          call cell%init(ierr, this%mesh%x(:,cn), face_vid, &
              PYR5_EDGES, this%mesh%volume(i), face_normal(:,1:size(fi)))

        case (6)
          ! zero treated as sentinel value in multimat_cell procedures
          face_vid = 0
          do j = 1, size(fi)
            face_vid(1:WED6_FSIZE(j),j) = WED6_FACES(WED6_XFACE(j):WED6_XFACE(j+1)-1)
          end do
          call cell%init(ierr, this%mesh%x(:,cn), face_vid, &
              WED6_EDGES, this%mesh%volume(i), face_normal(:,1:size(fi)))

        case (8)
          call cell%init(ierr, this%mesh%x(:,cn), reshape(source=HEX8_FACES,shape=[4,6]), &
              HEX8_EDGES, this%mesh%volume(i), face_normal(:,1:size(fi)))

        case default
          call TLS_fatal('unaccounted topology in donor_fluxes_nd')
        end select

        if (ierr /= 0) call TLS_fatal('cell_outward_volflux failed: could not initialize cell')

        call cell%partition(vof(:,i), this%normal(:,:,i), this%cutoff, this%priority, &
            this%location_tol, this%location_iter_max)

        this%flux_vol_sub(:,this%mesh%xcface(i):this%mesh%xcface(i+1)-1) = &
            cell%outward_volflux(dt, vel(this%mesh%xcface(i):this%mesh%xcface(i+1)-1),&
                                 this%mesh%area(fi), this%cutoff, this%nfluid, ierr)
        if (ierr /= 0) call TLS_fatal('cell_outward_volflux failed')
      end associate
    end do

    call stop_timer('reconstruct/advect')

  end subroutine donor_fluxes_nd


  ! Scan all faces with an outward flux and determine if any material is over-exhausted from this
  ! cell.  If so lower the fluxes until the material is just exhausted. Then loop over the faces
  ! and balance the individual material fluxes with the total face flux. The sum of the material
  ! volume fluxes (Volume_Flux_Sub) for each face should sum to the total volume flux for that
  ! face.  This balancing is an iterative adjust procedure such that by the end, two criteria are
  ! satisfied: 1) the sum of material fluxes from each face of the cell equals the total face
  ! volume flux and 2) The cumulative sum of individual material fluxes (from this and previous
  ! volume_track_subcycles) does not exceed the volume of a particular material originally within a
  ! cell.  To make this happen, we decrease some fluxes to equal the volume of material still
  ! available to be fluxed, and increase other fluxes appropriately.
  !
  ! This routine assumes that the total fluxed volume is dt*v*face_area.  This is obviously not
  ! true for partially solid/void cells.  In these cases, the fluxes will be balanced upwards
  ! towards a ficticiously large flux volume.  This routine also permits the creation of volume
  ! fluxes for materials which may not be present in the cell.
  subroutine flux_renorm(this, vel, vof_n, flux_vol, dt)

    use near_zero_function

    class(volume_tracker), intent(inout) :: this
    real(r8), intent(in) :: vel(:), vof_n(:,:), flux_vol(:,:), dt

    integer :: i,j,o,m,f0,f1,navail,nmat
    logical  :: adjust_fluxes, avail(size(vof_n,dim=1))
    real(r8) :: mat_flux_cur, mat_flux_acc, mat_avail
    real(r8) :: face_flux_fixed, face_flux_adjustable, face_flux

    nmat = this%nrealfluid

    do i = 1, this%mesh%ncell
      avail = .true.
      f0 = this%mesh%xcface(i)
      f1 = this%mesh%xcface(i+1)-1

      do o = 1, nmat+1
        adjust_fluxes = .false.

        do m = 1, nmat
          mat_flux_cur = sum(this%flux_vol_sub(m,f0:f1))
          if (near_zero(mat_flux_cur)) cycle

          mat_flux_acc = 0.0_r8
          do j = f0, f1
            if (flux_vol(m,j) > 0.0_r8) mat_flux_acc = mat_flux_acc + flux_vol(m,j)
          end do
          mat_avail = vof_n(m,i)*this%mesh%volume(i)
          ! check/correct for overflow
          if ((mat_flux_acc + mat_flux_cur) > mat_avail) then
            adjust_fluxes = .true.
            avail(m) = .false.
            ! rescale so that mat_flux_acc + mat_flux_cur == mat_avail
            this%flux_vol_sub(m,f0:f1) = this%flux_vol_sub(m,f0:f1) * &
                (mat_avail-mat_flux_acc)/mat_flux_cur
          end if
        end do

        if (.not.adjust_fluxes) exit

        do j = f0, f1
          face_flux = dt*vel(j)*this%mesh%area(this%mesh%cface(j))
          if (face_flux > this%cutoff*this%mesh%volume(i)) then
            face_flux_fixed = 0.0_r8
            face_flux_adjustable = 0.0_r8
            do m = 1, nmat
              if (avail(m)) then
                face_flux_adjustable = face_flux_adjustable + this%flux_vol_sub(m,j)
              else
                face_flux_fixed = face_flux_fixed + this%flux_vol_sub(m,j)
              end if
            end do

            if (face_flux_adjustable > 0.0_r8) then
              ! rescale so face_flux_adj+face_flux_fixed == face_flux
              do m = 1, nmat
                if (avail(m)) this%flux_vol_sub(m,j) = this%flux_vol_sub(m,j) * &
                    (face_flux - face_flux_fixed)/face_flux_adjustable
              end do
            else
              navail = count(avail)
              if (navail == 0) then
                call TLS_fatal ('FLUX_RENORM: cannot reassign face flux to any other material')
              end if

              ! arbitrarily add volume flux to potentially non-existent material in cell to balance
              ! equations.  This seems really, really, really, really, really, really bad.
              do m = 1, nmat
                if (avail(m)) this%flux_vol_sub(m,j) = (face_flux-face_flux_fixed)/real(navail,r8)
              end do
            end if
          end if
        end do
      end do
    end do

  end subroutine flux_renorm


  ! On entrance, this%flux_vol_sub only contains outward (i.e. positive) flux volumes.  On exit,
  ! this%flux_vol_sub is made consistent with appropriate negative entries
  subroutine flux_acceptor(this)

    class(volume_tracker), intent(inout) :: this

    integer :: m,j,i,f0,f1,nmat

    nmat = this%nfluid
    this%w_face(1:nmat,:) = 0.0_r8

    do i = 1, this%mesh%ncell
      f0 = this%mesh%xcface(i)
      f1 = this%mesh%xcface(i+1)-1
      do j = f0, f1
        do m = 1, nmat
          if (this%flux_vol_sub(m,j) > 0.0_r8) &
              this%w_face(m,this%mesh%cface(j)) = this%flux_vol_sub(m,j)
        end do
      end do
    end do

    do i = 1, this%mesh%ncell_onP
      f0 = this%mesh%xcface(i)
      f1 = this%mesh%xcface(i+1)-1
      do j = f0, f1
        do m = 1, nmat
          if (this%flux_vol_sub(m,j) /= this%w_face(m,this%mesh%cface(j))) &
              this%flux_vol_sub(m,j) = -this%w_face(m,this%mesh%cface(j))
        end do
      end do
    end do

  end subroutine flux_acceptor


  subroutine accumulate_volume(this, vof, flux_vol)

    class(volume_tracker), intent(inout) :: this
    real(r8), intent(inout) :: flux_vol(:,:), vof(:,:)

    integer :: i,j,f0,f1,m

    do i = 1, this%mesh%ncell_onP
      f0 = this%mesh%xcface(i)
      f1 = this%mesh%xcface(i+1)-1

      do j = f0, f1
        do m = 1, this%nrealfluid
          flux_vol(m,j) = flux_vol(m,j) + this%flux_vol_sub(m,j)
          vof(m,i) = vof(m,i) - this%flux_vol_sub(m,j)/this%mesh%volume(i)
        end do
      end do
    end do

  end subroutine accumulate_volume


  ! Enforce boundedness by allowing _inconsistent_ material flux volumes at faces
  subroutine enforce_bounded_vof(this, vof, flux_vol, fluids, void)

    class(volume_tracker), intent(inout) :: this
    real(r8), intent(inout) :: vof(:,:), flux_vol(:,:)
    integer, intent(in) :: fluids, void

    integer :: i,m,f0,f1
    real(r8) :: a1, v, q, excess

    a1 = 1.0_r8 - this%cutoff

    do i = 1, this%mesh%ncell_onP
      f0 = this%mesh%xcface(i)
      f1 = this%mesh%xcface(i+1)-1

      do m = 1, fluids
        if (vof(m,i) > a1 .and. vof(m,i) /= 1.0_r8) then
          call adjust_flux_matl(flux_vol(m,f0:f1), this%mesh%volume(i)*(1.0_r8-vof(m,i)))
          vof(m,i) = 1.0_r8
        else if (vof(m,i) < this%cutoff .and. vof(m,i) /= 0) then
          call adjust_flux_matl(flux_vol(m,f0:f1), -this%mesh%volume(i)*vof(m,i))
          vof(m,i) = 0.0_r8
        end if
      end do

      do m = fluids+1, fluids+void
        if (vof(m,i) > a1 .and. vof(m,i) /= 1.0_r8) then
          vof(m,i) = 1.0_r8
        else if (vof(m,i) < this%cutoff .and. vof(m,i) /= 0) then
          vof(m,i) = 0.0_r8
        end if
      end do

      v = sum(vof(:,i))
      excess = v - 1.0_r8
      if (excess == 0.0_r8) cycle

      q = sum(vof(fluids+1:fluids+void,i))

      if (q > 0.0_r8) then
         ! we can add or remove enough void from cell
        if (excess < 0.0_r8 .or. (excess > 0.0_r8 .and. q >= excess)) then
          do m = fluids+1, fluids+void
            vof(m,i) = vof(m,i) * (1.0_r8 - excess/q)
          end do
        else ! we cannot remove enough void from cell
          do m = fluids+1, fluids+void
            vof(m,i) = 0.0_r8
          end do
          call adjust_flux_all(flux_vol(:,f0:f1), vof(:,i), q-excess, this%mesh%volume(i), fluids)
        end if
        cycle
      end if

      ! There is no void to adjust in this cell.
      call adjust_flux_all(flux_vol(:,f0:f1), vof(:,i), -excess, this%mesh%volume(i), fluids)
    end do

    call gather_boundary(this%mesh%cell_ip, vof)

  end subroutine enforce_bounded_vof


  ! Adjust the material volume fluxes on the faces of a single cell to match the evaluated material
  ! volume to a target value. vof_delta is desired_vof-actual_vof
  ! IGNORES BC's FOR NOW
  ! The result of this process is an _inconsitent_ view of material flux volumes
  subroutine adjust_flux_matl(flux_vol, vol_delta)

    real(r8), intent(inout) :: flux_vol(:)
    real(r8), intent(in) :: vol_delta

    integer        :: i
    real(r8)       :: flux_vol_adj, r

    ! both inflow and outflow will be rescaled so take abs
    flux_vol_adj = sum(abs(flux_vol))

    if (flux_vol_adj > 0.0_r8) then
      r = vol_delta/flux_vol_adj

      do i = 1, size(flux_vol)
        if (flux_vol(i) > 0.0_r8) then
          flux_vol(i) = (1.0_r8-r)*flux_vol(i)
        else
          flux_vol(i) = (1.0_r8+r)*flux_vol(i)
        end if
      end do

    else if (vol_delta < 0.0_r8) then
      ! jms Note:   If the material is to be removed from the cell
      ! look for a face that doesn't have incoming material, and
      ! flux it out through that face
      do i = 1, size(flux_vol)
        if (flux_vol(i) >= 0.0_r8) then
          flux_vol(i) = flux_vol(i) - vol_delta
          exit
        end if
      end do
      !
      ! WHY IS THERE NO CHECK THAT THE FLUXES HAVE ACTUALLY BEEN ADJUSTED?
      !
    else
      !
      ! WHY NOT ARBITRARILY FLUX VOLUME IN AS WELL IF NEED BE?
      !
    end if

  end subroutine adjust_flux_matl

  ! Adjust all the material volume fluxes on the faces of a single cell to match the evaluated
  ! material volume to a target value. vof_delta is desired_vof-actual_vof
  ! IGNORES BC's FOR NOW The
  ! result of this process is an _inconsitent_ view of material flux volumes
  subroutine adjust_flux_all(flux_vol, vof, vof_delta, vol, fluids)

    real(r8), intent(inout) :: flux_vol(:,:), vof(:)
    real(r8), intent(in) :: vof_delta, vol
    integer, intent(in) :: fluids

    integer :: i, j
    real(r8) :: flux_vol_adj, excess, r

    ! both inflow and outflow will be rescaled so take abs. Ignore BCs for now
    flux_vol_adj = sum(abs(flux_vol(1:fluids,:)))

    if (flux_vol_adj > 0.0_r8) then
      r = vof_delta/flux_vol_adj

      do j = 1, size(flux_vol,dim=2)
        do i = 1, fluids
          if (flux_vol(i,j) > 0.0_r8) then
            vof(i) = vof(i) + r*flux_vol(i,j)
            flux_vol(i,j) = (1.0_r8-r*vol)*flux_vol(i,j)
          else
            vof(i) = vof(i) - r*flux_vol(i,j)
            flux_vol(i,j) = (1.0_r8+r*vol)*flux_vol(i,j)
          end if
        end do
      end do
    else
      ! truchas prints a warning... why is the matl routine different?
      return
    end if

  end subroutine adjust_flux_all

end module volume_tracker_type
