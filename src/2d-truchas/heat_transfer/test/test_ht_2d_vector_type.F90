program test_ht_2d_vector_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mpi_f08, only: MPI_COMM_WORLD, MPI_Comm_rank, MPI_Comm_size
  use parallel_communication
  use simulation_environment_type
  use unstr_2d_mesh_type
  use unstr_2d_mesh_factory
  use vector_class
  use ht_2d_vector_type
  implicit none

  type(unstr_2d_mesh), pointer :: mesh
  real(r8) :: xmin(2), xmax(2)
  integer :: nx(2)
  integer :: status, stat
  character(:), allocatable :: errmsg
  type(simulation_environment) :: env

  call init_parallel_communication
  env%comm = MPI_COMM_WORLD
  call MPI_Comm_rank(env%comm, env%rank)
  call MPI_Comm_size(env%comm, env%nproc)
  call env%simlog%init(env%comm, 'test_ht_2d_vector_type.log', stat, errmsg, terminal_output=.false.)
  if (stat /= 0) then
    if (is_IOP) print '(2a)', 'FAIL: ', errmsg
    call halt_parallel_communication
    stop 1
  end if
  xmin = [0.0_r8, 0.0_r8]
  xmax = [1.0_r8, 1.0_r8]
  nx = [8, 8]
  mesh => new_unstr_2d_mesh(env, xmin, xmax, nx, 0.0_r8)

  status = 0
  call test_storage_and_gather(mesh)
  call test_vector_operations(mesh)

  call halt_parallel_communication
  stop status

contains

  subroutine test_storage_and_gather(mesh)
    type(unstr_2d_mesh), target, intent(in) :: mesh

    type(ht_2d_vector) :: u
    real(r8), allocatable :: hc_offp(:), tc_offp(:), tf_offp(:)
    integer :: j

    call u%init(mesh)
    call require(size(u%hc) == mesh%ncell, 'incorrect cell enthalpy storage size')
    call require(size(u%tc) == mesh%ncell, 'incorrect cell temperature storage size')
    call require(size(u%tf) == mesh%nface, 'incorrect face temperature storage size')

    do j = 1, mesh%ncell_onP
      u%hc(j) = real(mesh%cell_imap%global_index(j), r8)
      u%tc(j) = -real(mesh%cell_imap%global_index(j), r8)
    end do
    do j = 1, mesh%nface_onP
      u%tf(j) = 2.0_r8*real(mesh%face_imap%global_index(j), r8)
    end do
    call u%gather_offp

    call require(all(u%hc == real(mesh%cell_imap%global_index([(j, j=1,mesh%ncell)]), r8)), &
                 'cell enthalpy gather failed')
    call require(all(u%tc == -real(mesh%cell_imap%global_index([(j, j=1,mesh%ncell)]), r8)), &
                 'cell temperature gather failed')
    call require(all(u%tf == 2.0_r8*real(mesh%face_imap%global_index([(j, j=1,mesh%nface)]), r8)), &
                 'face temperature gather failed')

    if (mesh%ncell_onP < mesh%ncell) then
      hc_offp = u%hc(mesh%ncell_onP+1:)
      tc_offp = u%tc(mesh%ncell_onP+1:)
    end if
    if (mesh%nface_onP < mesh%nface) tf_offp = u%tf(mesh%nface_onP+1:)
    call u%setval(7.0_r8)

    call require(all(u%hc(:mesh%ncell_onP) == 7.0_r8), 'setval failed for cell enthalpy')
    call require(all(u%tc(:mesh%ncell_onP) == 7.0_r8), 'setval failed for cell temperature')
    call require(all(u%tf(:mesh%nface_onP) == 7.0_r8), 'setval failed for face temperature')
    if (mesh%ncell_onP < mesh%ncell) then
      call require(all(u%hc(mesh%ncell_onP+1:) == hc_offp), 'setval modified off-process enthalpy')
      call require(all(u%tc(mesh%ncell_onP+1:) == tc_offp), 'setval modified off-process temperature')
    end if
    if (mesh%nface_onP < mesh%nface) &
        call require(all(u%tf(mesh%nface_onP+1:) == tf_offp), 'setval modified off-process face temperature')
  end subroutine test_storage_and_gather


  subroutine test_vector_operations(mesh)
    type(unstr_2d_mesh), target, intent(in) :: mesh

    type(ht_2d_vector) :: x, y, z
    class(vector), allocatable :: clone
    integer :: n

    call x%init(mesh)
    call y%init(mesh)
    call z%init(mesh)
    call x%setval(2.0_r8)
    call y%setval(3.0_r8)
    call y%update(4.0_r8, x)
    call require_owned_value(y, 11.0_r8, 'single-vector update failed')

    call y%update(2.0_r8, x, 3.0_r8)
    call require_owned_value(y, 37.0_r8, 'two-term update failed')

    call z%setval(5.0_r8)
    call z%update(2.0_r8, x, 3.0_r8, y)
    call require_owned_value(z, 120.0_r8, 'three-term update failed')

    call z%update(2.0_r8, x, 3.0_r8, y, 4.0_r8)
    call require_owned_value(z, 595.0_r8, 'four-term update failed')

    n = global_sum(2*mesh%ncell_onP + mesh%nface_onP)
    call require(abs(x%dot(x) - 4.0_r8*real(n, r8)) == 0.0_r8, 'incorrect dot product')
    call require(abs(x%norm1() - 2.0_r8*real(n, r8)) == 0.0_r8, 'incorrect one norm')
    call require(abs(x%norm2() - 2.0_r8*sqrt(real(n, r8))) < 1.0e-12_r8, 'incorrect two norm')
    call require(x%norm_max() == 2.0_r8, 'incorrect maximum norm')

    call x%clone(clone)
    call clone%copy(x)
    select type (clone)
    class is (ht_2d_vector)
      call require_owned_value(clone, 2.0_r8, 'clone or copy failed')
    class default
      call require(.false., 'incorrect clone dynamic type')
    end select
  end subroutine test_vector_operations


  subroutine require_owned_value(u, value, message)
    type(ht_2d_vector), intent(in) :: u
    real(r8), intent(in) :: value
    character(*), intent(in) :: message

    call require(all(u%hc(:u%mesh%ncell_onP) == value) .and. &
                 all(u%tc(:u%mesh%ncell_onP) == value) .and. &
                 all(u%tf(:u%mesh%nface_onP) == value), message)
  end subroutine require_owned_value


  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message

    if (global_any(.not.condition)) then
      if (is_IOP) print '("ERROR: ",a)', message
      status = 1
    end if
  end subroutine require

end program test_ht_2d_vector_type
