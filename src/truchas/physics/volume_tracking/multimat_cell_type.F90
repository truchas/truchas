!!
!! MULTIMAT_CELL_TYPE
!!
!! This module provides a cell type that
!! describes internal material geometries
!! as child arbitrary polyhedra
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! March 2017
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module multimat_cell_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use polygon_type
  use polyhedron_type
  use truchas_logging_services
  implicit none
  private

  ! a multimat_cell is a polyhedron itself, describing
  ! the cell geometry, and also contains an array
  ! of polyhedra each describing the geometry of a
  ! particular material
  type, public :: multimat_cell
    integer                       :: nmat,m ! number of materials actually present in cell
    type(polyhedron) :: geom
    type(polyhedron), allocatable :: mat_poly(:)
  contains
    procedure, private :: init_from_polyhedron
    procedure, private :: init_tet
    generic :: init => init_from_polyhedron, init_tet
    procedure :: partition
    procedure, private :: volumes_behind_plane
    procedure :: outward_volflux
    procedure :: interface_polygon
  end type multimat_cell

contains

  subroutine init_from_polyhedron (this, ierr, x, face_v, edge_v, face_normal, vol, tesselate)

    class(multimat_cell), intent(out) :: this
    integer, intent(out) :: ierr
    real(r8), intent(in)  :: x(:,:)
    integer, intent(in)  :: face_v(:,:), edge_v(:,:)
    real(r8), optional, intent(in)  :: face_normal(:,:), vol
    logical, optional, intent(in) :: tesselate

    call this%geom%init(ierr, x, face_v, edge_v, face_normal, vol, tesselate)

  end subroutine init_from_polyhedron

  subroutine init_tet (this, ierr, x, face_normal, vol, set_face_normals)
    class(multimat_cell), intent(out) :: this
    integer, intent(out) :: ierr
    real(r8), intent(in) :: x(:,:)
    real(r8), optional, intent(in) :: face_normal(:,:), vol
    logical, optional, intent(in) :: set_face_normals
    call this%geom%init(ierr, x, face_normal, vol, set_face_normals)
  end subroutine init_tet

  ! given a set of VoFs, normals, and an order,
  ! create child polyhedra for each material

  subroutine partition (this, vof, norm, cutvof, priority, max_reconstruction_iterations)

    use near_zero_function
    use plane_type
    use locate_plane_nd_module

    class(multimat_cell), intent(inout) :: this
    real(r8),             intent(in)    :: vof(:), norm(:,:), cutvof
    integer, intent(in) :: priority(:), max_reconstruction_iterations

    type(plane)      :: interface_plane
    type(polyhedron) :: tmp(2),remainder
    integer          :: m,p,nm,ierr

    ierr = 0
    if (allocated(this%mat_poly)) deallocate(this%mat_poly)
    allocate(this%mat_poly(size(vof)), stat=ierr)
    this%mat_poly(:)%parent%nVerts = 0
    this%m = 0

    remainder = this%geom

    this%nmat = count(vof > cutvof)
    nm = 0

    do p = 1,size(vof)
      m = priority(p)
      if (vof(m) < cutvof) cycle
      nm = nm+1 ! update the counter of how many materials we've seen thus far

      ! reconstruct the plane from the remaining free space
      ! use the plane to generate the polyhedron for this material,
      ! and update the free-space polyhedron
      if (nm==this%nmat .or. (1-cutvof)*remainder%volume() < vof(m)*this%geom%volume() .or. &
          near_zero(remainder%volume() / this%geom%volume())) then
        ! if this is the final material in the cell,
        ! or its volume is within a cutvof of the remaining volume,
        ! it gets the entire remainder of the polyhedron
        this%mat_poly(m) = remainder
        if (this%nmat==1) this%m = m ! if this is the only material, store its ID
        exit
      else
        ! if this is not the final material in the cell, split the cell
        interface_plane = locate_plane_nd(remainder, norm(:,m), vof(m)*this%geom%volume(), &
            this%geom%volume(), cutvof, max_reconstruction_iterations)
        call remainder%split (interface_plane,tmp,ierr)

        ! this check ensures the partitions give their vofs within the requested cutvof
        ! it will fail if the Brent's iterations did not converge within the maximum
        ! number of allowed iterations
        ! if (.not.near_zero(tmp(1)%volume()/this%volume() - (1.0_r8 - sum(vof(:m))), cutvof)) ierr=1

        if (ierr/=0) call partitionError ()

        remainder = tmp(1)
        this%mat_poly(m) = tmp(2)
      end if
    end do

  contains

    subroutine partitionError ()

      integer :: i
      real(r8) :: rem

      write(*,*)
      write(*,*) 'partition error!'

      write(*,*) 'cell:'
      call this%geom%print_data ()

      write(*,*) 'vof: ',vof

      write(*,*) 'norms:'
      do i = 1,size(norm, dim=2)
        write(*,'(3es20.10)') norm(:,i)
      end do

      write(*,*) 'previous vofs:'
      rem = 1.0_r8
      do i = 1,m-1
        write(*,'(i3,a,2es20.10)') i,': ',this%mat_poly(i)%volume() / this%geom%volume()
        rem = rem - this%mat_poly(i)%volume() / this%geom%volume()
      end do
      write(*,*) 'remaining vof: ',rem
      rem = rem - tmp(2)%volume() / this%geom%volume()
      write(*,*) 'remaining after last cutout vof: ',rem

      write(*,*)
      write(*,*) 'previous remainder polyhedron data: '
      call remainder%print_data ()

      write(*,*)
      write(*,*) 'material ',m,' cutout: '
      call tmp(2)%print_data ()

      write(*,*)
      write(*,*) 'remainder polyhedron data: '
      call tmp(1)%print_data ()

      write(*,*) 'interface plane:'
      call interface_plane%print_data ()

      write(*,*)
      ! write(*,*) 'remainder polyhedron volume does not match remaining volume from vof'
      ! write(*,'(a,es20.10)') 'remainder polyhedron vof: ',tmp(1)%volume() / this%geom%volume()
      ! write(*,'(a,es20.10)') 'exact remaining vof:      ',1.0_r8 - sum(vof(:m))
      ! write(*,'(a,es20.10)') 'error:                    ',&
      !     abs(tmp(1)%volume() / this%geom%volume() - (1.0_r8 - sum(vof(:m))))

      call TLS_fatal ("partition error")
    end subroutine partitionError

  end subroutine partition

  ! given a plane, find the volumes of each
  ! material behind that plane (flux volumes)
  !
  ! note: this will loop through the number of elements
  !       in the input 'vol' array, so that some last-priority
  !       materials may be skipped (such as solid).

  subroutine volumes_behind_plane (this, P, vol, ierr)

    use plane_type

    class(multimat_cell), intent(inout) :: this
    class(plane),         intent(in) :: P
    real(r8), intent(out) :: vol(:)
    integer,              intent(out) :: ierr

    integer                          :: m

    ASSERT(size(vol) <= size(this%mat_poly))

    ierr = 0
    do m = 1,size(vol)
      vol(m) = this%mat_poly(m)%volume_behind_plane (P,ierr)
      if (ierr /= 0) then
        write(*,*) 'volumes_behind_plane failed on material ',m
        return
      end if
    end do

  end subroutine volumes_behind_plane

  function outward_volflux (this, adv_dt, fluxing_velocity, face_area, cutvof, nfluid, ierr)

    use plane_type

    class(multimat_cell), intent(inout) :: this !inout because of call to volume
    real(r8),             intent(in)    :: adv_dt, fluxing_velocity(:), face_area(:), cutvof
    integer, intent(in) :: nfluid
    integer, intent(out) :: ierr
    real(r8)                            :: outward_volflux(nfluid,size(face_area))

    type(plane) :: flux_plane
    real(r8)    :: xf(3)
    integer     :: f,nV,m

    ierr = 0
    do f = 1,size(face_area)
      if (fluxing_velocity(f)*adv_dt*face_area(f) < cutvof*this%geom%volume()) then
        ! if we would be fluxing very very little, don't flux anything
        outward_volflux(:,f) = 0.0_r8
      else
        ! find the plane equation for the back end of the flux volume
        ! WARNING: in general, this could be non-planar, just like cell faces
        flux_plane%normal = -this%geom%parent%face_normal(:,f)

        nV = count(this%geom%parent%face_vid(:,f) /= 0) ! number of vertices on this face
        xf = sum(this%geom%parent%x(:,this%geom%parent%face_vid(1:nV,f)),dim=2) / nV ! face center

        flux_plane%rho  = sum(xf*flux_plane%normal) + adv_dt * fluxing_velocity(f)

        ! find the volume of the materials behind the flux plane
        if (this%nmat==1) then
          ! pure hex cells are easy, don't cut up polyhedrons to calculate the flux volume
          outward_volflux(:,f) = 0.0_r8
          outward_volflux(this%m,f) = adv_dt * fluxing_velocity(f) * face_area(f)
        else
          ! calculate the volumes of materials behind the flux plane,
          ! considering only the first nfluid materials (indicated by
          ! the size of the first dim of outward_volflux). Any material
          ! after index nfluid is solid, by convention.
          call this%volumes_behind_plane (flux_plane, outward_volflux(:,f), ierr)
          if (ierr /= 0) then
            write(*,*) 'outward_volflux failed on face ',f
            return
          end if
        end if

        ! make sure we have a valid outward volume flux
        if (any(outward_volflux(:,f) < 0.0_r8)) then
          write(*,*)
          write(*,'(a,i6,4es14.4)') 'f,volflux: ',f,outward_volflux(:,f)
          write(*,'(a,es14.4)') 'correct tot volflux: ', adv_dt * fluxing_velocity(f) * face_area(f)

          write(*,'(a,4es20.10)') 'flux plane n,p: ',flux_plane%normal, flux_plane%rho
          write(*,*) 'nmat ',this%nmat, this%m

          call this%geom%print_data ()

          call TLS_fatal ("in nested dissection outward_volflux: negative fluxes calculated!")
        end if
      end if
    end do

  end function outward_volflux

  ! returns a polygon for a desired interface
  type(polygon) function interface_polygon (this,m)

    use polygon_type

    class(multimat_cell), intent(in) :: this
    integer,              intent(in) :: m

    integer :: interface_face_id,nVerts

    interface_polygon%nVerts = 0
    ! if this polyhedron doesn't exist, or the polyhedron describes a pure cell,
    ! there is no interface to find
    if (m > size(this%mat_poly) .or. this%nmat<2 .or. this%mat_poly(m)%parent%nVerts<4) return

    ! by the convention set in polyhedron_type%polyhedron_on_side_of_plane,
    ! the face corresponding to the phase interface is the last face in the polyhedron
    interface_face_id = this%mat_poly(m)%parent%nFaces

    ! count how many real vertices are listed for this face (0s represent non-existent vertices)
    nVerts = count(this%mat_poly(m)%parent%face_vid(:,interface_face_id)/=0)

    ! initialize the polyhedron with the vertices used by the interface face
    ! and the corresponding normal vector
    call interface_polygon%init(&
        this%mat_poly(m)%parent%x(:,this%mat_poly(m)%parent%face_vid(1:nVerts,interface_face_id)))
    !this%mat_poly(m)%face_normal(:,interface_face_id))

  end function interface_polygon

end module multimat_cell_type
