#include "f90_assert.fpp"

module vof_2d_test_driver

#ifdef NAGFOR
  use,intrinsic :: f90_unix, only: exit
#endif

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use truchas_logging_services
  use unstr_2d_mesh_type
  use volume_tracker_2d_class
  use xdmf_file_type
  implicit none

  !! Bundle up all the driver state data as a singleton THIS of private
  !! derived type.  All procedures use/modify this object.
  type :: driver_data
    type(unstr_2d_mesh), pointer :: mesh => null()  ! reference only -- do not own
    class(volume_tracker_2d), allocatable :: vt
    real(r8), allocatable :: fvof_i(:,:) ! fluid/void volume fractions at start of update
    real(r8), allocatable :: fvof_o(:,:) ! fluid/void volume fractions at end of update
    real(r8), allocatable :: flux_vol(:,:) ! flux volumes
    real(r8), allocatable :: flux_vel(:) ! fluxing velocity - ragged form
  end type driver_data
  type(driver_data), allocatable, target :: this

contains

  subroutine timestep_driver(tsmax, dt, mesh, vel_fn, nmat, nvtrack, problem_vel, vof, &
    outfile, int_normal, axisym, myproc)

    use simple_volume_tracker_type
    use geometric_volume_tracker_type
    integer, intent(in) :: tsmax, nmat, nvtrack
    type(unstr_2d_mesh), intent(in), target :: mesh
    logical, intent(in) :: axisym
    procedure(), pointer :: problem_vel
    real(r8), intent(inout) :: vof(:,:)
    real(r8), intent(inout) :: dt, vel_fn(mesh%nface)
    real(r8), intent(out) :: int_normal(:,:,:)
    real(r8), intent(in) :: myproc(:)
    type(xdmf_file), intent(inout) :: outfile

    integer :: i, j, k, f0, f1, its, void
    real(r8) :: tot_t

    allocate(this)

    this%mesh => mesh
    void = 0
    allocate(this%fvof_i(nmat, mesh%ncell))
    allocate(this%fvof_o(nmat, mesh%ncell))
    allocate(this%flux_vol(nmat, size(mesh%cface)))
    allocate(this%flux_vel(size(mesh%cface)))

    !! allocate required volume tracker derived class
    if (nvtrack == 1) then
      call TLS_info('  Simple volume tracker:', TLS_VERB_NORMAL)
      allocate(simple_volume_tracker :: this%vt)
    else if (nvtrack == 2) then
      call TLS_info('  Geometric volume tracker:', TLS_VERB_NORMAL)
      allocate(geometric_volume_tracker :: this%vt)
    else
      write(*,*) "Incorrect vtrack option: ", nvtrack
      call exit
    end if

    call this%vt%init(this%mesh, nmat, nmat, nmat, axisym)

    !! Time-stepping loop
    tot_t = 0.0_r8
    do its = 1, tsmax
      this%flux_vol = 0.0_r8

      !! get_vof
      this%fvof_i = vof

      !! Initialize the face-based normal-velocity field
      call problem_vel(mesh, tot_t, vel_fn)

      !! copy face velocities into cell-oriented array for VOF
      do i = 1, mesh%ncell
        f0 = mesh%cstart(i)
        f1 = mesh%cstart(i+1)-1
        do j = f0, f1
          k = mesh%cface(j)
          if (btest(mesh%cfpar(i), 1+j-f0)) then ! normal points inward
            this%flux_vel(j) = -vel_fn(k)
          else
            this%flux_vel(j) = vel_fn(k)
          end if
        end do
      end do

      call this%vt%flux_volumes(this%flux_vel, this%fvof_i, this%fvof_o, this%flux_vol, &
                                int_normal, nmat, void, dt)

      !! put_vof
      vof = this%fvof_o

      !! flow time
      tot_t = tot_t+dt

      !if (mod(its,100)==0 .or. its==1 .or. its==tsmax) write(*,*) its, tot_t

      if (mod(its,100)==0 .or. its==1 .or. its==tsmax) then
        !! Write P and V fields at a time snapshot
        call outfile%begin_variables(tot_t)
        call outfile%write_cell_var(vof(1,:), 'VOF1')
        call outfile%write_cell_var(vof(2,:), 'VOF2')
        call outfile%write_cell_var(myproc(:), 'mype')
        call outfile%write_cell_var(int_normal(1,1,:), 'x-normal')
        call outfile%write_cell_var(int_normal(2,1,:), 'y-normal')
        call outfile%end_variables
      end if
    end do !its

  end subroutine timestep_driver

  !! Face-based normal-velocity field with constant uniform velocity
  subroutine constant_vel(mesh, t, vel_fn)
    type(unstr_2d_mesh), intent(in), target :: mesh
    real(r8), intent(in) :: t
    real(r8), intent(inout) :: vel_fn(mesh%nface)
    integer :: j
    do j = 1, mesh%nface
      vel_fn(j) = dot_product([2.0_r8, 1.0_r8], mesh%normal(:,j)/mesh%area(j))
    end do
  end subroutine constant_vel

  !! Face-based normal-velocity field with vortical velocity
  subroutine vortex_vel(mesh, t, vel_fn)

    type(unstr_2d_mesh), intent(in), target :: mesh
    real(r8), intent(in) :: t
    real(r8), intent(inout) :: vel_fn(mesh%nface)

    integer :: j
    real(r8) :: pi, T_period, xf, yf, ux, uy

    pi = 4.0_r8 * atan(1.0_r8)

    T_period = 8.0_r8

    do j = 1, mesh%nface
      xf = 0.5_r8 * (mesh%x(1,mesh%fnode(1,j)) + mesh%x(1,mesh%fnode(2,j)))
      yf = 0.5_r8 * (mesh%x(2,mesh%fnode(1,j)) + mesh%x(2,mesh%fnode(2,j)))
      ux = -2.0_r8 * (sin(pi*xf))**2 * sin(pi*yf) * cos(pi*yf) *cos(pi*t/T_period)
      uy =  2.0_r8 * (sin(pi*yf))**2 * sin(pi*xf) * cos(pi*xf) *cos(pi*t/T_period)
      vel_fn(j) = dot_product([ux, uy], mesh%normal(:,j)/mesh%area(j))
    end do

  end subroutine vortex_vel

  !! Face-based normal-velocity field with axisymmetric velocity
  subroutine axisymmetric_vel(mesh, t, vel_fn)

    type(unstr_2d_mesh), intent(in), target :: mesh
    real(r8), intent(in) :: t
    real(r8), intent(inout) :: vel_fn(mesh%nface)

    integer :: j
    real(r8) :: pi, T_period, rf, zf, ur, uz

    pi = 4.0_r8 * atan(1.0_r8)

    T_period = 1.0_r8

    do j = 1, mesh%nface
      rf = 0.5_r8 * (mesh%x(1,mesh%fnode(1,j)) + mesh%x(1,mesh%fnode(2,j)))
      ! this mesh velocity field can be used only for radii that are non-zero; since at a
      ! radius of zero, it becomes infinite.
      INSIST(abs(rf) > 1e-12_r8)
      zf = 0.5_r8 * (mesh%x(2,mesh%fnode(1,j)) + mesh%x(2,mesh%fnode(2,j)))
      ur = rf/2.0_r8 * cos(pi*t/T_period)
      uz = -zf * cos(pi*t/T_period)
      vel_fn(j) = dot_product([ur, uz], mesh%normal(:,j)/mesh%area(j))
    end do

  end subroutine axisymmetric_vel

end module vof_2d_test_driver
