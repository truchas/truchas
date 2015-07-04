#include "f90_assert.fpp"

module test_cell_grad_type_tools

  use kinds, only: r8
  implicit none
  private

  type, abstract, public :: test_func
  contains
    procedure(test_func_evaluate), deferred :: evaluate
  end type

  abstract interface
    subroutine test_func_evaluate (this, x, u, gradu)
      import test_func, r8
      class(test_func), intent(in) :: this
      real(r8), intent(in)  :: x(:)
      real(r8), intent(out) :: u, gradu(:)
    end subroutine
  end interface

  type, extends(test_func), public :: linear_func
    real(r8) :: c(3)
  contains
    procedure :: evaluate => linear_func_evaluate
  end type

  type, extends(test_func), public :: quad_func
    real(r8) :: c(9)
  contains
    procedure :: evaluate => quad_func_evaluate
  end type

  type, extends(test_func), public :: spher_func
    real(r8) :: r, w
  contains
    procedure :: evaluate => spher_func_evaluate
  end type

contains

  subroutine linear_func_evaluate (this, x, u, gradu)
    class(linear_func), intent(in) :: this
    real(r8), intent(in)  :: x(:)
    real(r8), intent(out) :: u, gradu(:)
    u = dot_product(this%c, x)
    gradu = this%c
  end subroutine

  subroutine quad_func_evaluate (this, x, u, gradu)
    class(quad_func), intent(in) :: this
    real(r8), intent(in)  :: x(:)
    real(r8), intent(out) :: u, gradu(:)
    u = x(1)*this%c(1) + x(2)*this%c(2) + x(3)*this%c(3) &
      + x(1)*(this%c(4)*x(1) + this%c(5)*x(2) + this%c(6)*x(3)) &
      + x(2)*(this%c(5)*x(1) + this%c(7)*x(2) + this%c(8)*x(3)) &
      + x(3)*(this%c(6)*x(1) + this%c(8)*x(2) + this%c(9)*x(3))
    gradu(1) = this%c(1) + 2.0_r8*(this%c(4)*x(1) + this%c(5)*x(2) + this%c(6)*x(3))
    gradu(2) = this%c(2) + 2.0_r8*(this%c(5)*x(1) + this%c(7)*x(2) + this%c(8)*x(3))
    gradu(3) = this%c(3) + 2.0_r8*(this%c(6)*x(1) + this%c(8)*x(2) + this%c(9)*x(3))
  end subroutine

  subroutine spher_func_evaluate (this, x, u, gradu)
    class(spher_func), intent(in) :: this
    real(r8), intent(in)  :: x(:)
    real(r8), intent(out) :: u, gradu(:)
    real(r8) :: z
    z = (sum(x**2)-this%r**2)/this%w
    u = tanh(z)
    gradu = (2.0_r8 / (this%w * cosh(z)**2)) * x
  end subroutine

end module test_cell_grad_type_tools


program test_linear_consistency

#ifdef NAGFOR
  use,intrinsic :: f90_unix, only: exit
#endif
  use kinds, only: r8
  use pgslib_module, only: PGSLib_CL_MAX_TOKEN_LENGTH, pgslib_finalize
  use parallel_util_module, only: parallel_init
  use parallel_communication
  use truchas_env, only: prefix
  use mesh_broker
  use truchas_logging_services
  use parameter_list_type
  use mfd_disc_type
  use cell_grad_type
  use unstr_mesh_type
  use test_cell_grad_type_tools
  implicit none

  character(PGSLib_CL_MAX_TOKEN_LENGTH), pointer :: argv(:) => null()
  integer :: status = 0

  character(:), allocatable :: indir

  call parallel_init (argv)
  call init_parallel_communication
  prefix='run'  ! TLS will write to 'run.log'
  call TLS_initialize
  call TLS_set_verbosity (TLS_VERB_NOISY)

  call process_command_line (indir)

  call make_meshes (indir)

  call test1_unif
  call test1_rand
  call test1_pave

  call test2_unif
  !call test2_pave

  !call test3_unif

  call pgslib_finalize
  call exit (status)

contains

  !! Regular mesh on [-1,1]^2. 1 cell thick in z with symmetry conditions on
  !! the z=const surfaces.  We expect to recover the exact cell gradient
  !! for linear functions (in x,y) with the imposed passive BC.  For higher
  !! order polynomials the passive BC is inexact, leading to larger errors
  !! on boundary cells.

  subroutine test1_unif

    logical, allocatable :: mask(:)
    integer, allocatable :: setids(:)
    character(:), allocatable :: errmsg
    class(test_func), allocatable :: func
    type(unstr_mesh), pointer :: mesh
    type(mfd_disc),  pointer :: disc
    type(cell_grad) :: grad
    integer :: stat

    call TLS_info ('')
    call TLS_info ('TEST1_UNIF: 2D square/uniform mesh')

    mesh => unstr_mesh_ptr('mesh1-unif')
    INSIST(associated(mesh))

    allocate(disc)
    call disc%init (mesh, use_new_mfd=.false.)

    allocate(mask(mesh%ncell), setids(2))
    mask = .true.
    setids = [5, 6]

    call grad%init (disc, mask, setids, stat, errmsg)
    if (stat /= 0) then
      call TLS_error ('from CELL_GRAD%INIT: ' // errmsg)
      status = 1
      return
    end if

    call TLS_info ('Recovering linear function gradient...')
    !F2008: func = linear_func([real(r8)::1.1,0.9,0])
    allocate(func, source=linear_func([real(r8)::1.1,0.9,0]))
    call run_test (mesh, grad, func, 1.0e-13_r8, 'test1-unif-lin.gmv')

    !! Here normal second derivative is 0 and passive BC is exact.
    call TLS_info ('Recovering quadratic function gradient...')
    !F2008: func = quad_func([real(r8)::1,2,0,0,1,0,0,0,0])
    deallocate(func)
    allocate(func, source=quad_func([real(r8)::1,2,0,0,1,0,0,0,0]))
    call run_test (mesh, grad, func, 1.0e-13_r8, 'test1-unif-quad1.gmv')

    !! Here normal second derivative is non-zero and passive BC is inexact.
    !! The error is exactly 0.05 on boundary cells and 0 on others.
    call TLS_info ('Recovering another quadratic function gradient...')
    !F2008: func = quad_func([real(r8)::1,2,0,1,0,0,1,0,0])
    deallocate(func)
    allocate(func, source=quad_func([real(r8)::1,2,0,1,0,0,1,0,0]))
    call run_test (mesh, grad, func, 5.1e-2_r8, 'test1-unif-quad2.gmv')

  end subroutine test1_unif

  !! Same problem as test1_unif but with the mesh node positions perturbed,
  !! yielding a non-orthogonal mesh.  Would still hope to be exact for
  !! linear functions if we were interpolating at cell centroids (center of
  !! mass), but we are using simple cell centers (average of cell vertices)
  !! which is not quite the same.  NB: cells do not have planar faces because
  !! the top surface grid is perturbed differently from the bottom surface.
  !! This probably wouldn't give exactness; the mesh could be fixed though.

  subroutine test1_rand

    logical, allocatable :: mask(:)
    integer, allocatable :: setids(:)
    character(:), allocatable :: errmsg
    class(test_func), allocatable :: func
    type(unstr_mesh), pointer :: mesh
    type(mfd_disc),  pointer :: disc
    type(cell_grad) :: grad
    integer :: stat

    call TLS_info ('')
    call TLS_info ('TEST1_RAND: 2D square/randomized mesh')

    mesh => unstr_mesh_ptr('mesh1-rand2')
    INSIST(associated(mesh))

    allocate(disc)
    call disc%init (mesh, use_new_mfd=.false.)

    allocate(mask(mesh%ncell), setids(2))
    mask = .true.
    setids = [5, 6]

    call grad%init (disc, mask, setids, stat, errmsg)
    if (stat /= 0) then
      call TLS_error ('from CELL_GRAD%INIT: ' // errmsg)
      status = 1
      return
    end if

    call TLS_info ('Recovering linear function gradient...')
    !F2008: func = linear_func([real(r8)::1.1,0.9,0])
    allocate(func, source=linear_func([real(r8)::1.1,0.9,0]))
    call run_test (mesh, grad, func, 1.0e-2_r8, 'test1-rand-lin.gmv')

  end subroutine test1_rand

  !! Same problem as test1_unif but using an unstructured, non-orthogonal
  !! paved mesh.  Comments for test1_rand apply here too.  Note that the
  !! cells are better shaped here (so that the center is closer to the
  !! centroid) and this seems to be reflected in a much more "accurate"
  !! solution for linear functions.

  subroutine test1_pave

    logical, allocatable :: mask(:)
    integer, allocatable :: setids(:)
    character(:), allocatable :: errmsg
    class(test_func), allocatable :: func
    type(unstr_mesh), pointer :: mesh
    type(mfd_disc),  pointer :: disc
    type(cell_grad) :: grad
    integer :: stat

    call TLS_info ('')
    call TLS_info ('TEST1_PAVE: 2D square/paved mesh')

    mesh => unstr_mesh_ptr('mesh1-pave')
    INSIST(associated(mesh))

    allocate(disc)
    call disc%init (mesh, use_new_mfd=.false.)

    allocate(mask(mesh%ncell), setids(2))
    mask = .true.
    setids = [5, 6]

    call grad%init (disc, mask, setids, stat, errmsg)
    if (stat /= 0) then
      call TLS_error ('from CELL_GRAD%INIT: ' // errmsg)
      status = 1
      return
    end if

    call TLS_info ('Recovering linear function gradient...')
    !F2008: func = linear_func([real(r8)::1.1,0.9,0])
    allocate(func, source=linear_func([real(r8)::1.1,0.9,0]))
    call run_test (mesh, grad, func, 7.0e-5_r8, 'test1-pave-lin.gmv')

  end subroutine test1_pave

  !! Regular mesh on [-1.2,1.2]^2. 1 cell thick in z with symmetry conditions
  !! on the z=const surfaces.  The outer 0.2-wide border is void and does not
  !! participate in the calculation. We expect to recover the exact cell gradient
  !! for linear functions (in x,y) with the imposed passive BC.  For higher
  !! order polynomials the passive BC is inexact, leading to larger errors
  !! on boundary cells.

  subroutine test2_unif

    logical, allocatable :: mask(:)
    integer, allocatable :: setids(:)
    character(:), allocatable :: errmsg
    class(test_func), allocatable :: func
    type(unstr_mesh), pointer :: mesh
    type(mfd_disc),  pointer :: disc
    type(cell_grad) :: grad
    integer :: stat

    call TLS_info ('')
    call TLS_info ('TEST2_UNIF: 2D square/uniform mesh with void border')

    mesh => unstr_mesh_ptr('mesh2-unif')
    INSIST(associated(mesh))

    allocate(disc)
    call disc%init (mesh, use_new_mfd=.false.)

    allocate(mask(mesh%ncell), setids(2))
    mask = (mesh%cblock == 1)
    setids = [5, 6]

    call grad%init (disc, mask, setids, stat, errmsg)
    if (stat /= 0) then
      call TLS_error ('from CELL_GRAD%INIT: ' // errmsg)
      status = 1
      return
    end if

    call TLS_info ('Recovering linear function gradient...')
    !F2008: func = linear_func([real(r8)::1.1,0.9,0])
    allocate(func, source=linear_func([real(r8)::1.1,0.9,0]))
    call run_test (mesh, grad, func, 1.0e-13_r8, 'test2-unif-lin.gmv')

    !! Here normal second derivative is 0 and passive BC is exact.
    call TLS_info ('Recovering quadratic function gradient...')
    !F2008: func = quad_func([real(r8)::1,2,0,0,1,0,0,0,0])
    deallocate(func)
    allocate(func, source=quad_func([real(r8)::1,2,0,0,1,0,0,0,0]))
    call run_test (mesh, grad, func, 1.0e-13_r8, 'test2-unif-quad1.gmv')

    !! Here normal second derivative is non-zero and passive BC is inexact.
    !! The error is exactly 0.05 on boundary cells and 0 on others.
    call TLS_info ('Recovering another quadratic function gradient...')
    !F2008: func = quad_func([real(r8)::1,2,0,1,0,0,1,0,0])
    deallocate(func)
    allocate(func, source=quad_func([real(r8)::1,2,0,1,0,0,1,0,0]))
    call run_test (mesh, grad, func, 5.1e-2_r8, 'test2-unif-quad2.gmv')

  end subroutine test2_unif


  subroutine run_test (mesh, grad, func, tol, filename)

    type(unstr_mesh) :: mesh
    type(cell_grad) :: grad
    class(test_func) :: func
    real(r8) :: tol
    character(*) :: filename  ! for gmv output

    integer :: j, stat
    real(r8) :: xc(3), maxerr
    real(r8) :: ucell(mesh%ncell), gradu(3,mesh%ncell), error(3,mesh%ncell)
    character(:), allocatable :: errmsg
    character(80) :: string
    
    !call centroid (reshape([real(r8)::0,0,1,0,2,2,0,1],shape=[2,4]),xc(1:2))
    !print *,xc
    !stop

    xc(3) = 0.0_r8
    do j = 1, mesh%ncell
      if (grad%cell_mask(j)) then
        !xc = sum(mesh%x(:,mesh%cnode(:,j)),dim=2) / size(mesh%cnode,dim=1)
        !print *, mesh%cnode(:,j)
        !print '(8f6.2)', mesh%x(1,mesh%cnode(:,j))
        !print '(8f6.2)', mesh%x(2,mesh%cnode(:,j))
        !print '(8f6.2)', mesh%x(3,mesh%cnode(:,j))
        !stop
        !print *, mesh%cnode([1,4,3,2],j)
        !print '(8f6.2)', mesh%x(1,mesh%cnode([1,4,3,2],j))
        !print '(8f6.2)', mesh%x(2,mesh%cnode([1,4,3,2],j))
        !print '(8f6.2)', mesh%x(3,mesh%cnode([1,4,3,2],j))
        call centroid (mesh%x(1:2,mesh%cnode([1,4,3,2],j)), xc(1:2))
        call func%evaluate (xc, ucell(j), gradu=error(:,j))
      else
        ucell(j) = 0.0_r8
        error(:,j) = 0.0_r8
      end if
    end do

    call grad%compute (ucell, gradu, stat, errmsg)
    if (stat /= 0) then
      status = 1  ! signal overall unit test failure
      call TLS_error ('from CELL_GRAD%COMPUTE: ' // errmsg)
    end if

    error = gradu - error
    call write_gmv_output (mesh, ucell, gradu, error, filename)

    maxerr = global_maxval(abs(error))
    write(string,'(a,es9.2)') 'max error = ', maxerr
    call TLS_info (string)

    if (maxerr > tol) then
      status = 1  ! signal overall unit test failure
      write(string,'(a,es10.2)') 'error exceeds tolerance:', tol
      call TLS_error (string)
    end if

  end subroutine run_test
  
  subroutine centroid (x, xc)
    real(r8), intent(in)  :: x(:,:)
    real(r8), intent(out) :: xc(:)
    real(r8) :: a1, a2, a3, a4, a
    !print *, x(1,:)
    !print *, x(2,:)
    a1 = 0.5_r8 * ((x(1,2)-x(1,1))*(x(2,4)-x(2,1)) - (x(2,2)-x(2,1))*(x(1,4)-x(1,1)))
    a2 = 0.5_r8 * ((x(1,3)-x(1,2))*(x(2,1)-x(2,2)) - (x(2,3)-x(2,2))*(x(1,1)-x(1,2)))
    a3 = 0.5_r8 * ((x(1,4)-x(1,3))*(x(2,2)-x(2,3)) - (x(2,4)-x(2,3))*(x(1,2)-x(1,3)))
    a4 = 0.5_r8 * ((x(1,1)-x(1,4))*(x(2,3)-x(2,4)) - (x(2,1)-x(2,4))*(x(1,3)-x(1,4)))
    a = 0.5_r8 * (a1 + a2 + a3 + a4)
    print *, abs((a1+a3)-(a2+a4))/a
    xc = 0.25_r8 * sum(x,dim=2) !+ ((a3-a1)*(x(:,3)-x(:,1)) + (a4-a2)*(x(:,4)-x(:,2)))/(12.0_r8*a)
    !print *, a1, a2, a3, a4
    !stop
  end subroutine centroid

  subroutine write_gmv_output (mesh, ucell, gradu, error, filename)
    use unstr_mesh_gmv
    type(unstr_mesh), intent(in) :: mesh
    real(r8), intent(in) :: ucell(:), gradu(:,:), error(:,:)
    character(*), intent(in) :: filename
    call gmv_open (filename)
    call gmv_write_unstr_mesh (mesh)
    call gmv_begin_variables
    call gmv_write_dist_cell_var (mesh, ucell(:mesh%ncell_onP), 'u')
    call gmv_write_dist_cell_var (mesh, gradu(1,:mesh%ncell_onP), 'dudx')
    call gmv_write_dist_cell_var (mesh, gradu(2,:mesh%ncell_onP), 'dudy')
    call gmv_write_dist_cell_var (mesh, gradu(3,:mesh%ncell_onP), 'dudz')
    call gmv_write_dist_cell_var (mesh, error(1,:mesh%ncell_onP), 'error-x')
    call gmv_write_dist_cell_var (mesh, error(2,:mesh%ncell_onP), 'error-y')
    call gmv_write_dist_cell_var (mesh, error(3,:mesh%ncell_onP), 'error-z')
    call gmv_end_variables
    call gmv_close
  end subroutine write_gmv_output

  !! Handle a single optional argument which is the directory where the mesh
  !! files are located.  This is needed for cmake/ctest which puts and runs
  !! the executable from a different directory that this source code and mesh
  !! files.  Also handle '-h' or '--help' as an option to write out the usage.

  subroutine process_command_line (indir)

    use,intrinsic :: iso_fortran_env, only: output_unit

    character(:), allocatable, intent(out) :: indir

    character(:), allocatable :: prog
    character(256) :: arg
    integer :: n, stat

    call get_command_argument (0, arg)
    n = scan(arg, '/', back=.true.) ! remove the leading path component, if any
    prog = trim(arg(n+1:))

    stat = 0
    select case (command_argument_count())
    case (0)
      indir = '.'
    case (1)
      call get_command_argument (1, arg)
      select case (arg)
      case ('-h','--help')
        stat = 1
      case default
        indir = trim(arg)
      end select
    case default
      stat = 1
    end select

    if (stat /= 0) then
      if (is_IOP) then
        write(output_unit,'(a)') 'Usage: ' // prog // ' [INDIR]'
        write(output_unit,'(a)') 'INDIR is the directory containing the input files (default ".")'
      end if
      call exit (stat) ! do not want this to count as a successful test
    end if

  end subroutine process_command_line

  !! Instantiate the meshes that the tests will use

  subroutine make_meshes (indir)

    character(*), intent(in) :: indir

    type(parameter_list) :: plist
    type(parameter_list), pointer :: params

    params => plist%sublist('mesh1-unif')
    call params%set ('mesh-file', indir // '/mesh1-unif.gen')

    params => plist%sublist('mesh1-rand')
    call params%set ('mesh-file', indir // '/mesh1-rand.gen')

    params => plist%sublist('mesh1-rand2')
    call params%set ('mesh-file', indir // '/mesh1-rand2.gen')

    params => plist%sublist('mesh1-pave')
    call params%set ('mesh-file', indir // '/mesh1-pave.gen')

    params => plist%sublist('mesh2-unif')
    call params%set ('mesh-file', indir // '/mesh2-unif.gen')
    call params%set ('interface-side-sets', [2])

!    params => plist%sublist('mesh2-pave')
!    call params%set ('mesh-file', indir // '/mesh2-pave.gen')
!    call params%set ('interface-side-sets', [2])
!
!    params => plist%sublist('mesh3-unif')
!    call params%set ('mesh-file', indir // '/mesh3-unif.gen')

    call init_mesh_broker (plist)

  end subroutine make_meshes

end program test_linear_consistency
