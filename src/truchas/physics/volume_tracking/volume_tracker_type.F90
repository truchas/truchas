module volume_tracker_type

  use truchas_logging_services
  use timing_tree
  use unstr_mesh_type
  implicit none
  private

  public :: volume_tracker

  type :: volume_tracker
    integer :: location_iter_max ! maximum number of iterations to use in fitting interface
    integer :: subcyles
    logical :: use_brents_method
    real(r8) :: location_tol ! tolerance of plic fit
    real(r8) :: cutoff ! allow volume fraction {0,(cutoff,1]}
    real(r8), allocatable :: flux_vol_sub(:,:), normal(:,:,:)
    type(unstr_mesh), pointer :: mesh ! unowned reference
    ! node/face/cell workspace
    real(r8), allocatable :: w_node(:,:), w_face(:,:), w_cell(:,:)
  contains
    procedure :: read_params
    procedure :: init
    procedure :: flux_volumes
    procedure :: advect
  end type volume_tracker


contains

  subroutine read_params(this, p)
    use parameter_list_type
    type(volume_tracker), intent(inout) :: this
    type(parameter_list), intent(in) :: p
    call p%get('location_iter_max', this%location_iter_max)
    call p%get('location_tol', this%location_tol)
    call p%get('cutoff', this%cutoff)
    call p%get('subcycles', this%subcyles)
    call p%get('use_brents_method', this%use_brents_method)
  end subroutine read_params

  subroutine init(this, m, fluids_and_void)
    type(volume_tracker), intent(inout) :: this
    type(unstr_mesh), intent(in) :: m
    integer, intent(in) :: fluids_and_void

    this%mesh = m
    allocate(this%flux_vol_sub(fluid_and_void,size(m%cface)))
    allocate(this%normal(fluid_and_void,m%ncell))
    ! workspace
    allocate(w_node(?,m%nnode))
    allocate(w_face(?,m%nface))
    allocate(w_cell(?,m%ncell))
  end subroutine init


  ! flux volumes routine assuming vel/flux_vol is a cface-like array
  subroutine flux_volumes(this, vel, vof_n, vof, flux_vol, fluids, void, dt)
    type(volume_tracker), intent(in) :: this
    real(r8), intent(in) :: vel(:), vof_n(:,:), dt
    real(r8), intent(out) :: flux_vol(:,:), vof(:,:)
    integer, intent(in) :: fluids, void

    integer :: i
    real(r8) :: sub_dt

    flux_vol = 0.0_r8
    sub_dt = dt/real(this%subcycle, r8)

    do i = 1, this%subcycles
      call this%normals(vof, fluids, void)

      call this%donor_fluxes_nd(vel, vof, sub_dt)

      call this%flux_renorm(vel, vof_n, flux_vol, fluids, void, sub_dt)

      call this%flux_acceptor(volume_flux_sub)

      call this%flux_bc(fluxing_velocity, vof_n, volume_flux_sub)

      ! add the volume fluxes from this subcycle (volume_flux_sub) to the
      ! total flux array (volume_flux_tot), and update the volume fraction
      ! array (vof).
      call this%volume_advance(volume_flux_sub, volume_flux_tot, vof)
    end do

  end subroutine flux_volumes


  subroutine normals(this, vof, fluids, void)

    use flow_operators, only: gradient_cc
    intrinsic :: norm2, findloc

    class(volume_tracker), intent(inout) :: this
    real(r8), intent(in)  :: vof(:,:)
    integer, intent(in) :: fluids, void

    integer :: i,j,k,c
    logical :: hasvof(size(vof,dim=1))

    call start_timer('normals')
    do i = 1 , fluids+void
      call gradient_cc(this%mesh, this%normal(1,i,:), this%normal(2,i,:), this%normal(3,i,:), &
          vof(i,:), this%w_node(1,:), this%w_node(2,:), this%w_face(1,:))
    end do

    do i = 1, this%mesh%ncell_onP
      hasvof = vof(:,i) > 0.0_r8
      c = count(hasvof)
      if (c < 2) then
        this%normal(:,:,i) = 1.0_r8
        cycle
      end if

      this%normal(:,:,i) = -this%normal(:,:,i)
      ! enforce consistency for two materials
      if (c == 2) then
        j = findloc(hasvof,.true.)
        k = findloc(hasvof(j+1:),.true.)
        this%normal(:,k,i) = -this%normal(:,j,i)
      endif

      ! normalize and remove smallish components due to robustness issues in nested disection
      do j = 1 , fluids+void
        if (vof(j,i) <= 0.0_r8) cycle ! should this be cutoff?

        this%normal(:,j,i) = this%normal(:,j,i)/norm2(this%normal(:,j,i))
        do k = 1, 3
          if (abs(this%normal(k,j,i)) < 1.0e-6_r8) this%normal(k,j,i) = 0.0_r8
        end do
        this%normal(:,j,i) = this%normal(:,j,i)/norm2(this%normal(:,j,i))
      end do
    end do
    call end_timer('normals')

  end subroutine normals


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

    do i = 1, this%mesh%ncell_onP
      associate (cn => this%mesh%cnode(this%mesh%xcnode(i):this%mesh%xcnode(i+1)-1), &
          fi => this%mesh%cface(this%mesh%xcface(i):this%mesh%xcface(i+1)-1))

        do j = 1, size(fi)
          k = fi(j)
          if (btest(this%mesh%cfpar(i),pos=j)) then
            face_normal(:,j) = -this%mesh%normal(:,k)
          else
            face_normal(:,j) = this%mesh%normal(:,k)
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

        call cell%partition(vof(:,i), this%normal(:,:,i))

        this%flux_vol_sub(:,this%mesh%xcface(i):this%mesh%xcface(i+1)-1) = &
            cell%outward_volflux(dt, vel(this%mesh%xcface(i):this%mesh%xcface(i+1)-1),&
                                 this%mesh%area(fi), ierr)
        if (ierr /= 0) call TLS_fatal('cell_outward_volflux failed')
      end associate
    end do

    call stop_timer ('reconstruct/advect')

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
  subroutine flux_renorm (this, vel, vof_n, flux_vol, dt)

    use near_zero_function

    class(volume_tracker), intent(inout) :: this
    real(r8), intent(in) :: vel(:), vof_n(:,:), flux_vol(:,:), dt

    integer :: i,o,m,f0,f1,navail,nmat
    logical  :: adjust_fluxes, avail(size(vof_n,dim=1))
    real(r8) :: mat_flux_cur, mat_flux_acc, mat_avail
    real(r8) :: face_flux_fixed, face_flux_adjustable, face_flux

    nmat = size(vof_n,dim=1)

    do i = 1, this%mesh%ncell
      avail = .true.
      f0 = this%mesh%xcface(i)
      f1 = this%mesh%xcface(i+1)-1

      do o = 1, nmat+1
        adjust_fluxes = .false.

        do m = 1, nmat
          mat_flux_cur = sum(this%flux_vol_sub(m,f0:f1))
          if (near_zero(mat_flux_cur)) cycle

          mat_flux_acc = sum(flux_vol(m,f0:f1))
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
            do m = 1, fluids+void
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


  ! renorm_cell: ensure no materials are over-exhausted in a given cell
  !
  subroutine renorm_cell (volume_flux_sub, volume_flux_tot, vof_n, &
       total_face_flux, cell_volume, matl_is_fluid, ierr)



    real(r8), intent(inout) :: volume_flux_sub(:,:)
    real(r8), intent(in)    :: volume_flux_tot(:,:), vof_n(:), &
         total_face_flux(:), cell_volume
    logical, intent(in) :: matl_is_fluid(:)
    integer,  intent(out)   :: ierr

    real(r8) :: Ratio, total_material_flux, cumulative_outward_material_flux, &
        total_flux_through_face, total_flux_through_face_not_maxed
    integer  :: norm_iter, f, m, number_not_maxed
    logical  :: flux_reduced, maxed(size(vof_n))

    maxed = .false.
    ierr = 0

    ! Loop over the renorm_loop a maximum of nmat - 1 times, to resolve all instances
    ! of fluxing more material than was in the cell at the beginning of the timestep
    do norm_iter = 1, size(vof_n)+1
      flux_reduced = .false.

      ! see note 1

      ! The first step is to determine if any material is being over-exhausted from
      ! a cell.  If so mark it as MAXED and lower the Volume_Flux_Sub's so that the
      ! material is just exhausted.
      do m = 1,size(vof_n)
        ! volume of material m attempting to leave the cell in this volume_track_subcycle
        total_material_flux = sum(Volume_Flux_Sub(m,:))
        if (near_zero(total_material_flux)) cycle

        ! If the cumulative_outward_material_flux exceeds the amount of material
        ! originally in the cell (from vof_n), calculate the 'Ratio' of fluid material volume
        ! still allowed to be fluxed to the flux volume, and note that we're not 'Done'


        ! cumulative volume of material m that left the cell in previous subcycles
        ! why is the mask checking for volume_flux_tot > 0?  Why does that matter?
        cumulative_outward_material_flux = sum(volume_flux_tot(m,:), mask=volume_flux_tot(m,:) > 0)

        ! ratio between the volume of original (from beginning of flow cycle)
        ! material remaining in the cell (material that has entered is disregarded)
        ! and the volume of material m attempting to leave the cell in this volume track cycle
        if (matl_is_fluid(m)) then
          ratio = (vof_n(m)*cell_volume - cumulative_outward_material_flux) / total_material_flux
        else
          ratio = 0.0_r8
        end if

        if (Ratio < 1) then
          ! lower the fluxes to match the material volume within the cell,
          ! and flag the material number as maxed
          flux_reduced = .true.
          maxed(m) = .true.
          volume_flux_sub(m,:) = ratio * volume_flux_sub(m,:)
        end if
      end do ! material loop

      if (.not.flux_reduced) exit

      ! This cell had one/more fluxes reduced.  For each of the faces, if the sum
      ! of material fluxes is less than Total_Face_Flux, multiply all non-maxed
      ! fluxes by another 'Ratio' (this time > 1) that restores the flux balance.
      ! This may in turn over-exhaust one or more of these materials, and so from
      ! the bottom of this loop, we head back to the top.
      do f = 1,NFC
        ! Calculate the total flux volume through the cell face (is this already
        ! available elsewhere?), and if the flux volume is greater than zero,
        ! then concern ourselves with adjusting individual material fluxes.
        if (Total_Face_Flux(f) > cutvof*cell_volume) then
          ! calculate the sum of material fluxes at a face, and the sum of un-maxed material fluxes.
          total_flux_through_face           = sum(Volume_Flux_Sub(:,f))
          total_flux_through_face_not_maxed = sum(volume_flux_sub(:,f), mask=.not.maxed)

          ! Ratio as defined below, when used to multiply the non-maxed fluxes at
          ! a face, will restore the flux balance.
          if (total_flux_through_face_not_maxed > 0) then
            Ratio = (Total_Face_Flux(f) + total_flux_through_face_not_maxed &
                - total_flux_through_face) &
                / total_flux_through_face_not_maxed
            where (.not.maxed) Volume_Flux_Sub(:,f) = Ratio * Volume_Flux_Sub(:,f)
          else
            number_not_maxed = count(.not.Maxed.and.matl_is_fluid)
            if (number_not_maxed == 0) then
              ierr = 1
              return
            end if

            where (.not.Maxed.and.matl_is_fluid) &
                Volume_Flux_Sub(:,f) = (Total_Face_Flux(f) - total_flux_through_face) &
                / number_not_maxed
          end if

        end if
      end do ! face loop

    end do ! renorm loop

  end subroutine renorm_cell

end module volume_tracker_type
