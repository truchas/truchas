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
  subroutine flux_volumes(this, vel, vof_n, vof, flux_vol, fluids, void, dt, svof)
    type(volume_tracker), intent(in) :: this
    real(r8), intent(in) :: vel(:), vof_n(:,:), dt, svof(:)
    real(r8), intent(out) :: flux_vol(:,:), vof(:,:)
    integer, intent(in) :: fluids, void

    integer :: i
    real(r8) :: sub_dt

    flux_vol = 0.0_r8
    sub_dt = dt/real(this%subcycle, r8)

    do i = 1, this%subcycles
      call this%normals(vof, fluids, void)

      call this%donor_fluxes_nd(vel, vof, sub_dt)

      call this%flux_renorm(vel, vof_n, flux_vol, sub_dt)

      call this%flux_acceptor()

      call this%flux_bc(fluxing_velocity, vof_n, volume_flux_sub)

      call this%accumulate_volume(vof, flux_vol)

      call this%enforce_bounded_vof(vof, flux_vol, fluids, void, svof)
    end do

  end subroutine flux_volumes


  subroutine normals(this, vof)

    use flow_operators, only: gradient_cc
    intrinsic :: norm2, findloc

    class(volume_tracker), intent(inout) :: this
    real(r8), intent(in)  :: vof(:,:)

    integer :: i,j,k,c,nmat
    logical :: hasvof(size(vof,dim=1))

    nmat = size(vof,dim=1)

    call start_timer('normals')
    do i = 1 , nmat
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
      do j = 1 , nmat
        if (vof(j,i) <= 0.0_r8) cycle ! should this be cutoff?

        this%normal(:,j,i) = this%normal(:,j,i)/norm2(this%normal(:,j,i))
        do k = 1, 3
          if (abs(this%normal(k,j,i)) < 1.0e-6_r8) this%normal(k,j,i) = 0.0_r8
        end do
        this%normal(:,j,i) = this%normal(:,j,i)/norm2(this%normal(:,j,i))
      end do
    end do
    ! will need normals for vof reconstruction in ghost cells
    call gather_boundary(this%mesh%cell_ip, this%normal)
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

    do i = 1, this%mesh%ncell
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


  ! On entrance, this%flux_vol_sub only contains outward (i.e. positive) flux volumes.  On exit,
  ! this%flux_vol_sub is made consistent with appropriate negative entries
  subroutine flux_acceptor(this)

    class(volume_tracker), intent(inout) :: this

    integer :: m,j,i,f0,f1,nmat

    nmat = size(this%flux_vol_sub,dim=1)
    this%w_face(1:nmat,:) = 0.0_r8

    do i = 1, this%mesh%ncell
      f0 = this%mesh%xcface(i)
      f1 = this%mesh%xcface(i+1)-1
      do j = f0, f1
        do m = 1, nmat
          if (this%flux_vol_sub(m,j,i) > 0.0_r8) &
              this%w_face(m,this%mesh%cface(j)) = this%flux_vol_sub(m,j,i)
        end do
      end do
    end do

    do i = 1, this%mesh%ncell_onP
      f0 = this%mesh%xcface(i)
      f1 = this%mesh%xcface(i+1)-1
      do j = f0, f1
        do m = 1, nmat
          if (this%flux_vol_sub(m,j,i) /= this%w_face(m,this%mesh%cface(j))) &
              this%flux_vol_sub(m,j,i) = -this%w_face(m,this%mesh%cface(j))
        end do
      end do
    end do

  end subroutine flux_acceptor


  subroutine accumulate_volume(this, vof, flux_vol)

    class(volume_tracker), intent(inout) :: this
    real(r8), intent(inout) :: flux_vol(:,:), vof(:,:)

    integer :: i,nmat,j,f0,f1

    nmat = size(vof,dim=1)

    do i = 1, this%mesh%ncell_onP
      f0 = this%mesh%xcface(i)
      f1 = this%mesh%xcface(i+1)-1

      do j = f0, f1
        do m = 1, nmat
          flux_vol(m,j) = flux_vol(m,j) + this%flux_vol_sub(m,j)
          vof(m,j) = vof(m,j) - this%flux_vol_sub(m,j)/this%mesh%volume(i)
      end do
    end do

  end subroutine accumulate_volume


  ! Enforce boundedness by allowing _inconsistent_ material flux volumes at faces
  subroutine enforce_bounded_vof(this, vof, flux_vol, fluids, void, svof)

    class(volume_tracker), intent(inout) :: this
    real(r8), intent(inout) :: vof(:,:), flux_vol(:,:)
    real(r8), intent(in) :: svof(:)
    integer, intent(in) :: fluids, void

    integer :: i,m,f0,f1
    real(r8) :: a1, v, q, excess

    a1 = 1.0_r8 - this%cutoff

    do i = 1, this%mesh%ncell_onP
      f0 = this%mesh%xcface(i)
      f1 = this%mesh%xcface(i+1)-1

      do m = 1, fluids
        if (vof(m,i) > a1 .and. vof(m,i) /= 1.0_r8) then
          call adjust_matl_flux(flux_vol(m,f0:f1), this%mesh%volume(i)*(1.0_r8-vof(m,i)))
          vof(m,i) = 1.0_r8
        else if (vof(m,i) < this%cutoff .and. vof(m,i) /= 0) then
          call adjust_matl_flux(flux_vol(m,f0:f1), -this%mesh%volume(i)*vof(m,i))
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

      v = sum(vof(:,i)) + svof(i)
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
          call adjust_flux_total(flux_vol(:,f0:f1), vof(:,i), q-excess, this%mesh%volume(i), fluids)
        end if
        cycle
      end if

      ! There is no void to adjust in this cell.
      call adjust_flux_total(flux_vol(:,f0:f1), vof(:,i), -excess, this%mesh%volume(i), fluids)
    end do

  end subroutine enforce_bounded_vof


  ! Adjust the material volume fluxes on the faces of a single cell to match the evaluated material
  ! volume to a target value. vof_delta is desired_vof-actual_vof
  ! IGNORES BC's FOR NOW
  ! The result of this process is an _inconsitent_ view of material flux volumes
  subroutine adjust_flux_matl(flux_vol, vol_delta)

    real(r8), intent(inout) :: flux_vol(:)
    real(r8), intent(in) :: vol_delta

    integer        :: i
    real(r8)       :: flux_vol_adj

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
    real(r8) :: flux_vol_adj, excess

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
            vof(i) = vof(i) - r*flux_vof(i,j)
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
