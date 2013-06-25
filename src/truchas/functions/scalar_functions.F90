!!
!! SCALAR_FUNCTIONS
!!
!! This module provides a derived type for describing a general scalar-valued
!! function of a vector argument, and an associated evaluation method.
!!
!! PROGRAMMING INTERFACE
!!
!! The derived type SCAFUN serves as a container for one of serveral different
!! types of scalar function objects.  It should be thought of as an abstract
!! type or class, and instances as polymorphic variables with a specific
!! dynamic type, which are created through a call to one of the following
!! subroutines.
!!
!!  CALL CREATE_SCAFUN_CONST (THIS, C) creates a constant-valued SCAFUN
!!    object THIS whose evaluation will return the value C.
!!
!!  CALL CREATE_SCAFUN_USER (THIS, INDEX[, P]) creates a SCAFUN object THIS
!!    whose evaluation invokes a user-customized function named USER_SCAFUN
!!    having the following interface:
!!
!!      function user_scafun (index, x, p) result (f)
!!        use kinds, only: r8
!!        integer,  intent(in) :: index
!!        real(r8), intent(in) :: x(:), p(:)
!!        real(r8) :: f
!!      end function
!!
!!    The specified INDEX value is passed as the first function argument, and
!!    the values of the rank-1 real array P are passed as the third function
!!    argument it it is specified, otherwise a 0-sized array is passed.
!!    Within the function, the INDEX value is intended to be used as the
!!    selector in a select-case construct to choose between the code for
!!    multiple functions of variables whose values are passed in the X array
!!    and which may depend on parameter values passed in the P array.
!!    Note that the SCAFUN object THIS holds a copy of the values in the
!!    specified P array; subsequently changing the values in the P array has
!!    no effect on THIS, so that the function parameters are essentially fixed.
!!
!!  CALL CREATE_SCAFUN_DLL (THIS, LIB, SYM, P) creates a SCAFUN object THIS
!!    whose evaluation invokes a function from a dynamically linked library.
!!    The character strings LIB and SYM specify the path to the shared library
!!    and the symbol name within the library, respectively.  The function
!!    interface must be interoperable with the Fortran 77 style interface
!!
!!      function f(x, p)
!!        double precision, intent(in) :: x(*), p(*)
!!        double precision :: f
!!      end function
!!
!!    assuming that Fortran real kind R8 is equivalent to double precision.
!!    The values of the rank-1 real array P are passed as the second function
!!    argument if P is specified, otherwise a 0-sized array is passed.  These
!!    are fixed parameter values that the function may depend on.  The function
!!    argument values are passed in the array X.
!!
!!  CALL CREATE_SCAFUN_POLY (THIS, C, E[, X0]) creates a SCAFUN object THIS that
!!    describes the polynomial function $\sum_{j=1}^n c_j (x-x_0)^{e_j}$. The
!!    coefficients are passed in the real array C and the corresponding exponents
!!    in the integer array E.  The optional real argument X0 specifies the
!!    reference point, which defaults to 0.  Note that SCAFUN objects describe
!!    functions of a vector argument, and that this is a polynomial of a single
!!    variable whose value is taken to be the first component of the vector
!!    argument.
!!
!!  CALL CREATE_SCAFUN_MPOLY (THIS, C, E[, X0]) creates a SCAFUN object THIS that
!!    describes the multivariate polynomial function
!!    $\sum_{j=1}^n c_j (x - x_0)^{e_j}$, where $x=(x_1,\dots,x_d)$ and
!!    $e_j=(e_{j,1},\dots,e_{j,d})$ and $x^e=x_1^{e_1}\dots x_d^{e_d}$.
!!    The coefficients are passed in the real array C.  The exponents are passed
!!    is the rank-2 integer array E, where E(:,j) are the variable exponents
!!    corressponding to C(j).  The reference point is passed in the array X0,
!!    but defaults to 0 (vector) if it is not specified.
!!
!!  CALL CREATE_SCAFUN_TM_DENSITY (THIS, REF_DENS, REF_TEMP, CTE, STAT, ERRMSG)
!!    creates a SCAFUN object THIS that describes a temperature-dependent
!!    density function used by thermo-mechanics.  This function has the form
!!      \[ \rho(T) = \rho_0 \exp(-3 \int_{T_0}^{T} \alpha(T) dT) \]
!!    where $\rho_0$ and $T_0$ are given by the real arguments REF_DENS and
!!    REF_TEMP, respectively, and the linear coefficient of thermal expansion
!!    $\alpha(T)$ is given by the SCAFUN object CTE.  CTE must be of a form
!!    that CREATE_SCAFUN_ANTIDERIV can be used with, namely a constant function
!!    or a single-variable polynomial.  STAT returns a nonzero value if the
!!    the antiderivative of CTE cannot be created, and an expanatory error
!!    message is returned in ERRMSG.
!!
!! This module provides a defined assignment for SCAFUN variables that does a
!! deep copy so that the lhs of the assignment is an identical but independent
!! copy of the rhs;  subsequently modifying or destroying the rhs has no effect
!! on the value of the lhs.
!!
!! The primary method provided by the module is the evaluation of the function
!! described by a SCAFUN object:
!!
!!  EVAL(THIS, X) evaluates the function described by the SCAFUN object THIS
!!    at the vector argument values X.
!!
!! Other methods are:
!!
!!  DESTROY(THIS) deallocates all storage associated with the SCAFUN object
!!    THIS and returns it to its default initialization state.
!!
!!  DEFINED(THIS) returns true if the SCAFUN object is well-defined.  Mostly
!!    intended for assertion checks and debugging.
!!
!!  IS_CONST_SCAFUN(THIS) returns true if the SCAFUN function THIS is of
!!    constant type; that is, it was defined by the CREATE_SCAFUN_CONST method.
!!
!!  CALL CREATE_SCAFUN_ANTIDERIV (F, X0, G0, G, STAT, ERRMSG) creates the
!!    SCAFUN object G that is the antiderivative of the SCAFUN object F and
!!    that satisfies EVAL(G, X0) == G0.  The antiderivative is with respect to
!!    the first variable.  Only certain SCAFUN functions can be integrated in
!!    this way, namely constant functions and single variable polynomial
!!    functions (without a 1/x term).  STAT returns a nonzero value if the
!!    subroutine is unable to create the antiderivative function, and an
!!    explanatory error message is returned in ERRMSG.
!!
!!  CALL CREATE_SCAFUN_PRODUCT (F, G, FG, STAT, ERRMSG) creates the SCAFUN
!!    object FG that is the product of the SCAFUN objects F and G.  The
!!    product of only certain types of SCAFUN functions can be formed.  STAT
!!    returns a nonzero value if the subroutine is unable to create the
!!    product function, and an explanatory error message is assigned to ERRMSG.
!!    Currently F must be a constant function and G a constant, polynomial, or
!!    multivariate polynomial function.  The returned function is the same
!!    type as G.
!!

#include "f90_assert.fpp"

module scalar_functions

  use kinds
  use dynamic_linking_loader
  implicit none
  private

  public :: scafun, assignment(=)  ! polymorphic type
  public :: eval, destroy, defined
  public :: create_scafun_const
  public :: create_scafun_user
  public :: create_scafun_dll
  public :: create_scafun_poly
  public :: create_scafun_mpoly
  public :: create_scafun_sstep
  public :: create_scafun_tm_density
  public :: create_scafun_antideriv, is_const_scafun, create_scafun_product
  public :: scafun_list, scafun_vector_list, append_to_list, convert_list_to_array

  public :: user_scafun, dll_scafun ! external interface definitions

  !! Abstract base type (poor man's polymorphism).
  !! Deferred type-bound procedures: eval, destroy, defined
  type :: scafun
    private
    integer :: dyn_type = 0
    type(scafun_const), pointer :: const => null()
    type(scafun_user),  pointer :: user  => null()
    type(scafun_dll),   pointer :: dll   => null()
    type(scafun_poly),  pointer :: poly  => null()
    type(scafun_mpoly), pointer :: mpoly => null()
    type(scafun_sstep), pointer :: sstep => null()
    type(scafun_tm_density), pointer :: tm_density => null()
  end type scafun

  !! Extended type: constant function
  type :: scafun_const
    private
    real(r8) :: c = 0.0_r8
  end type
  integer, parameter :: TYPE_IS_CONST = 1

  !! Extended type: statically-linked indexed function (user customized).
  type :: scafun_user
    private
    integer :: index = 0
    real(r8), pointer :: p(:) => null()
  end type scafun_user
  integer, parameter :: TYPE_IS_USER = 2

  !! Extended type: dynamically-linked user-provided function.
  type :: scafun_dll
    private
    integer(C_PTR_KIND) :: so_handle = 0
    integer(C_PTR_KIND) :: f_addr = 0
    real(r8), pointer :: p(:) => null()
  end type scafun_dll
  integer, parameter :: TYPE_IS_DLL  = 3

  !! Extended type: single-variable polynomial function.
  type :: scafun_poly
    private
    integer  :: emin = 0    ! minimum exponent
    integer  :: emax = 0    ! maximum exponent
    real(r8) :: x0 = 0.0_r8 ! reference point
    real(r8), pointer :: c(:) => null() ! array of coefficients
  end type scafun_poly
  integer, parameter :: TYPE_IS_POLY = 4

  !! Extended type: multivariate polynomial function.
  type :: scafun_mpoly
    private
    real(r8), pointer :: x0(:) => null()      ! reference point
    integer,  pointer :: expon(:,:) => null() ! array of exponents
    real(r8), pointer :: coef(:) => null()    ! array of coefficients
  end type scafun_mpoly
  integer, parameter :: TYPE_IS_MPOLY = 5

  !! Extended type: smooth (single-variable) step function
  type :: scafun_sstep
    private
    real(r8) :: x0 = 0.0_r8, y0 = 0.0_r8, x1 = 1.0_r8, y1 = 0.0_r8
  end type scafun_sstep
  integer, parameter :: TYPE_IS_SSTEP = 6

  !! Extended type: phase density function for thermo-mechanics
  type :: scafun_tm_density
    private
    real(r8) :: d0 = 0.0_r8 ! reference density
    type(scafun) :: int_cte ! indefinite integral of the linear CTE
  end type scafun_tm_density
  integer, parameter :: TYPE_IS_TM_DENSITY = 7

  interface eval
    module procedure eval_scafun
  end interface

  interface destroy
    module procedure destroy_scafun
  end interface

  interface defined
    module procedure defined_scafun
  end interface

  interface assignment(=)
    module procedure scafun_copy
  end interface

  interface ! for the customizable statically-linked indexed function.
    function user_scafun (index, x, p) result (f)
      use kinds, only: r8
      integer,  intent(in) :: index
      real(r8), intent(in) :: x(:), p(:)
      real(r8) :: f
    end function
  end interface

  interface ! for the dynamically-linked function.
    function dll_scafun (x, p) result (f)
      use kinds, only: r8
      real(r8), intent(in) :: x(*), p(*)
      real(r8) :: f
    end function
  end interface

  interface ! for the procedure to invoke the dynamically-linked function.
    function call_dll_scafun (addr, x, p) result(f)
      use kinds, only: r8
      use dynamic_linking_loader, only: c_ptr_kind
#if (defined(COMPAQ_COMPILER) || defined(PATHSCALE_COMPILER))
      integer(c_ptr_kind) :: addr
#else
      integer(c_ptr_kind), value :: addr
#endif
      real(r8), intent(in) :: x(*), p(*)
      real(r8) :: f
    end function
  end interface

  !!
  !! Types and generic procedures dealing with lists of scafun-type variables.
  !!

  type :: scafun_list
    private
    type(scafun) :: f
    type(scafun_list), pointer :: next => null()
  end type scafun_list

  type :: scafun_vector_list
    private
    type(scafun), pointer :: f(:) => null()
    type(scafun_vector_list), pointer :: next => null()
  end type scafun_vector_list

  interface append_to_list
    module procedure append_to_list0, append_to_list1
  end interface

  interface convert_list_to_array
    module procedure convert_list_to_array0, convert_list_to_array1
  end interface

  interface destroy
    module procedure destroy_scafun_list, destroy_scafun_vector_list
  end interface

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! PROCEDURES TO CREATE SPECIFIC DYNAMIC TYPES OF A SCAFUN POLYMORPHIC OBJECT
 !!

  subroutine create_scafun_const (this, c)
    type(scafun), intent(out) :: this
    real(r8), intent(in) :: c
    this%dyn_type = TYPE_IS_CONST
    allocate(this%const)
    call define_scafun_const (this%const, c)
  end subroutine create_scafun_const

  subroutine create_scafun_user (this, index, p)
    type(scafun), intent(out) :: this
    integer, intent(in) :: index
    real(r8), intent(in), optional :: p(:)
    this%dyn_type = TYPE_IS_USER
    allocate(this%user)
    call define_scafun_user (this%user, index, p)
  end subroutine create_scafun_user

  subroutine create_scafun_dll (this, lib, sym, p)
    type(scafun), intent(out) :: this
    character(*), intent(in) :: lib, sym
    real(r8), intent(in), optional :: p(:)
    this%dyn_type = TYPE_IS_DLL
    allocate(this%dll)
    call define_scafun_dll (this%dll, lib, sym, p)
  end subroutine create_scafun_dll

  subroutine create_scafun_poly (this, c, e, x0)
    type(scafun), intent(out) :: this
    real(r8), intent(in) :: c(:)
    integer, intent(in) :: e(:)
    real(r8), intent(in), optional :: x0
    this%dyn_type = TYPE_IS_POLY
    allocate(this%poly)
    call define_scafun_poly (this%poly, c, e, x0)
  end subroutine create_scafun_poly

  subroutine create_scafun_mpoly (this, c, e, x0)
    type(scafun), intent(out) :: this
    real(r8), intent(in) :: c(:)
    integer, intent(in) :: e(:,:)
    real(r8), intent(in), optional :: x0(:)
    this%dyn_type = TYPE_IS_MPOLY
    allocate(this%mpoly)
    call define_scafun_mpoly (this%mpoly, c, e, x0)
  end subroutine create_scafun_mpoly

  subroutine create_scafun_sstep (this, x0, y0, x1, y1)
    type(scafun), intent(out) :: this
    real(r8), intent(in) :: x0, y0, x1, y1
    this%dyn_type = TYPE_IS_SSTEP
    allocate(this%sstep)
    call define_scafun_sstep (this%sstep, x0, y0, x1, y1)
  end subroutine create_scafun_sstep

  subroutine create_scafun_tm_density (this, ref_dens, ref_temp, cte, stat, errmsg)
    type(scafun), intent(out) :: this
    real(r8), intent(in) :: ref_dens, ref_temp  ! reference density and temperature
    type(scafun), intent(in) :: cte ! linear coefficient of thermal expansion
    integer, intent(out) :: stat
    character(*), intent(out) :: errmsg
    this%dyn_type = TYPE_IS_TM_DENSITY
    allocate(this%tm_density)
    call define_scafun_tm_density (this%tm_density, ref_dens, ref_temp, cte, stat, errmsg)
  end subroutine create_scafun_tm_density

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! DYNAMIC DISPATCH PROCEDURES FOR POLYMORPHIC SCAFUN OBJECTS
 !!

  elemental logical function defined_scafun (this)
    type(scafun), intent(in) :: this
    defined_scafun = defined_scafun_aux(this)
  end function defined_scafun

  pure recursive logical function defined_scafun_aux (this)
    type(scafun), intent(in) :: this
    defined_scafun_aux = .false.
    select case (this%dyn_type)
    case (TYPE_IS_CONST)
      if (associated(this%const)) defined_scafun_aux = defined_scafun_const(this%const)
    case (TYPE_IS_USER)
      if (associated(this%user))  defined_scafun_aux = defined_scafun_user(this%user)
    case (TYPE_IS_DLL)
      if (associated(this%dll))   defined_scafun_aux = defined_scafun_dll(this%dll)
    case (TYPE_IS_POLY)
      if (associated(this%poly))  defined_scafun_aux = defined_scafun_poly(this%poly)
    case (TYPE_IS_MPOLY)
      if (associated(this%mpoly)) defined_scafun_aux = defined_scafun_mpoly(this%mpoly)
    case (TYPE_IS_SSTEP)
      if (associated(this%sstep)) defined_scafun_aux = defined_scafun_sstep(this%sstep)
    case (TYPE_IS_TM_DENSITY)
      if (associated(this%tm_density)) defined_scafun_aux = defined_scafun_tm_density(this%tm_density)
    end select
  end function defined_scafun_aux

  elemental subroutine destroy_scafun (this)
    type(scafun), intent(inout) :: this
    call destroy_scafun_aux (this)
  end subroutine destroy_scafun
  
  pure recursive subroutine destroy_scafun_aux (this)
    type(scafun), intent(inout) :: this
    type(scafun) :: default
    if (associated(this%const)) then
      call destroy_scafun_const (this%const)
      deallocate(this%const)
    end if
    if (associated(this%user)) then
      call destroy_scafun_user (this%user)
      deallocate(this%user)
    end if
    if (associated(this%dll)) then
      call destroy_scafun_dll (this%dll)
      deallocate(this%dll)
    end if
    if (associated(this%poly)) then
      call destroy_scafun_poly (this%poly)
      deallocate(this%poly)
    end if
    if (associated(this%mpoly)) then
      call destroy_scafun_mpoly (this%mpoly)
      deallocate(this%mpoly)
    end if
    if (associated(this%sstep)) then
      call destroy_scafun_sstep (this%sstep)
      deallocate(this%sstep)
    end if
    if (associated(this%tm_density)) then
      call destroy_scafun_tm_density (this%tm_density)
      deallocate(this%tm_density)
    end if
    this = default  ! assign default initialization values
  end subroutine destroy_scafun_aux
  
  recursive function eval_scafun (this, x) result (f)
    type(scafun), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8) :: f
    ASSERT( defined(this) )
    select case (this%dyn_type)
    case (TYPE_IS_CONST)
      f = eval_scafun_const(this%const, x)
    case (TYPE_IS_USER)
      f = eval_scafun_user(this%user, x)
    case (TYPE_IS_DLL)
      f = eval_scafun_dll(this%dll, x)
    case (TYPE_IS_POLY)
      f = eval_scafun_poly(this%poly, x)
    case (TYPE_IS_MPOLY)
      f = eval_scafun_mpoly(this%mpoly, x)
    case (TYPE_IS_SSTEP)
      f = eval_scafun_sstep(this%sstep, x)
    case (TYPE_IS_TM_DENSITY)
      f = eval_scafun_tm_density(this%tm_density, x)
    end select
  end function eval_scafun
  
  elemental subroutine scafun_copy (dest, src)
    type(scafun), intent(out) :: dest
    type(scafun), intent(in)  :: src
    call scafun_copy_aux (dest, src)
  end subroutine scafun_copy

  pure recursive subroutine scafun_copy_aux (dest, src)
    type(scafun), intent(out) :: dest
    type(scafun), intent(in)  :: src
    dest%dyn_type = src%dyn_type
    select case (src%dyn_type)
    case (TYPE_IS_CONST)
      if (associated(src%const)) then
        allocate(dest%const)
        dest%const = src%const
      end if
    case (TYPE_IS_USER)
      if (associated(src%user)) then
        allocate(dest%user)
        call scafun_user_copy (dest%user, src%user)
      end if
    case (TYPE_IS_DLL)
      if (associated(src%dll)) then
        allocate(dest%dll)
        call scafun_dll_copy (dest%dll, src%dll)
      end if
    case (TYPE_IS_POLY)
      if (associated(src%poly)) then
        allocate(dest%poly)
        call scafun_poly_copy (dest%poly, src%poly)
      end if
    case (TYPE_IS_MPOLY)
      if (associated(src%mpoly)) then
        allocate(dest%mpoly)
        call scafun_mpoly_copy (dest%mpoly, src%mpoly)
      end if
    case (TYPE_IS_SSTEP)
      if (associated(src%sstep)) then
        allocate(dest%sstep)
        dest%sstep = src%sstep
      end if
    case (TYPE_IS_TM_DENSITY)
      if (associated(src%tm_density)) then
        allocate(dest%tm_density)
        call scafun_tm_density_copy (dest%tm_density, src%tm_density)
      end if
    end select
  end subroutine scafun_copy_aux

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SPECIFIC PROCEDURES FOR CONSTANT FUNCTIONS
 !!
 !! No defined assignment is needed; all components have intrinsic type.
 !!

  subroutine define_scafun_const (this, c)
    type(scafun_const), intent(out) :: this
    real(r8), intent(in) :: c
    this%c = c
  end subroutine define_scafun_const

  pure logical function defined_scafun_const (this)
    type(scafun_const), intent(in) :: this
    defined_scafun_const = .true.
  end function defined_scafun_const

  pure subroutine destroy_scafun_const (this)
    type(scafun_const), intent(inout) :: this
    type(scafun_const) :: default
    this = default  ! assign default initialization values
  end subroutine destroy_scafun_const

  function eval_scafun_const (this, x) result(f)
    type(scafun_const), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8) :: f
    f =  this%c
  end function eval_scafun_const

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SPECIFIC PROCEDURES FOR STATICALLY-LINKED, INDEXED FUNCTIONS
 !!
 !! No defined assignment is needed; all components have intrinsic type.
 !!

  subroutine define_scafun_user (this, index, p)
    type(scafun_user), intent(out) :: this
    integer, intent(in) :: index
    real(r8), intent(in), optional :: p(:)
    INSIST( index > 0 )
    this%index = index
    if (present(p)) then
      allocate(this%p(size(p)))
      this%p = p
    else
      allocate(this%p(0))
    end if
  end subroutine define_scafun_user

  pure logical function defined_scafun_user (this)
    type(scafun_user), intent(in) :: this
    defined_scafun_user = (this%index > 0) .and. associated(this%p)
  end function defined_scafun_user

  pure subroutine destroy_scafun_user (this)
    type(scafun_user), intent(inout) :: this
    type(scafun_user) :: default
    if (associated(this%p)) deallocate(this%p)
    this = default  ! assign default initialization values
  end subroutine destroy_scafun_user

  function eval_scafun_user(this, x) result(f)
    type(scafun_user), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8) :: f
    f =  user_scafun(this%index, x, this%p)
  end function eval_scafun_user

  pure subroutine scafun_user_copy (dest, src)
    type(scafun_user), intent(out) :: dest
    type(scafun_user), intent(in)  :: src
    dest%index = src%index
    if (associated(src%p)) then
      allocate(dest%p(size(src%p)))
      dest%p = src%p
    end if
  end subroutine scafun_user_copy

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SPECIFIC PROCEDURES FOR USER-PROVIDED, DYNAMICALLY-LOADED FUNCTIONS
 !!
 !! Ideally the shared library should be unloaded when the object is destroyed
 !! but this spawns several problems.  We need to be able to make copies of an
 !! instance of the object.  Because the type only has components of instrinsic
 !! type, we could use intrinsic assignment of the type to effect copies except
 !! that we need would need to increment the reference count in the shared
 !! object handle so that the library would not be prematurely unloaded (due to
 !! a destroy of the first instance, for example).  The only way to increase
 !! this reference count is to call dll_open again as part of a defined
 !! assignment, which would require the path of the shared library to be stored
 !! as a component of the structure.  To do this we have to impose a maximum
 !! string length (fixed in F2003), which is something that is otherwise not
 !! needed and rather distasteful.  It may be possible to maintain an
 !! independent reference count (for objects created by assignment), but for
 !! now we choose not to unload the shared library.
 !!

  subroutine define_scafun_dll (this, lib, sym, p)
    type(scafun_dll), intent(out) :: this
    character(*), intent(in) :: lib, sym
    real(r8), intent(in), optional :: p(:)
    if (scan(lib, '/') == 0) then
      call dll_open ('./'//lib, RTLD_NOW, this%so_handle)
    else
      call dll_open (lib, RTLD_NOW, this%so_handle)
    end if
    call dll_symbol (this%so_handle, sym, this%f_addr)
    if (present(p)) then
      allocate(this%p(size(p)))
      this%p = p
    else
      allocate(this%p(0))
    end if
  end subroutine define_scafun_dll

  pure logical function defined_scafun_dll (this)
    type(scafun_dll), intent(in) :: this
    defined_scafun_dll = (this%so_handle /= 0) .and. (this%f_addr /= 0) .and. associated(this%p)
  end function defined_scafun_dll

  pure subroutine destroy_scafun_dll (this)
    type(scafun_dll), intent(inout) :: this
    type(scafun_dll) :: default
    !if (this%so_handle /= 0) call dll_close (this%so_handle)
    if (associated(this%p)) deallocate(this%p)
    this = default  ! assign default initialization values
  end subroutine destroy_scafun_dll

  function eval_scafun_dll(this, x) result(f)
    type(scafun_dll), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8) :: f
#if (defined(COMPAQ_COMPILER) || defined(PATHSCALE_COMPILER))
    f = call_dll_scafun (%VAL(this%f_addr), x, this%p)
#else
    f = call_dll_scafun (this%f_addr, x, this%p)
#endif
  end function eval_scafun_dll

  pure subroutine scafun_dll_copy (dest, src)
    type(scafun_dll), intent(out) :: dest
    type(scafun_dll), intent(in)  :: src
    dest%so_handle = src%so_handle
    dest%f_addr = src%f_addr
    if (associated(src%p)) then
      allocate(dest%p(size(src%p)))
      dest%p = src%p
    end if
  end subroutine scafun_dll_copy

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SPECIFIC PROCEDURES FOR POLYNOMIAL FUNCTIONS
 !!
 !! This particular implementation avoids use of the exponentiation operator
 !! by using Horner's method for polynomial evaluation.  As a consequence it
 !! treats all powers from the minimum (or 0) to the maximum (or 0) as present
 !! with a zero-valued coefficient if necessary.  Thus it isn't terribly well
 !! suited to a polynomial like x**100.
 !!
 !! This type has a pointer component, and thus requires a defined assignment.
 !!

  pure subroutine scafun_poly_copy (dest, src)

    type(scafun_poly), intent(out) :: dest
    type(scafun_poly), intent(in)  :: src

    dest%emin = src%emin
    dest%emax = src%emax
    dest%x0   = src%x0
    if (associated(src%c)) then
      allocate(dest%c(lbound(src%c,1):ubound(src%c,1)))
      dest%c = src%c
    end if

  end subroutine scafun_poly_copy

  subroutine define_scafun_poly (this, c, e, x0)

    type(scafun_poly), intent(out) :: this
    real(r8), intent(in) :: c(:)
    integer,  intent(in) :: e(:)
    real(r8), intent(in), optional :: x0

    integer :: j

    INSIST( size(c) > 0 .and. size(c) == size(e) )

    this%emin = min(0, minval(e))
    this%emax = max(0, maxval(e))

    allocate(this%c(this%emin:this%emax))

    this%c = 0.0_r8
    do j = 1, size(e)
      this%c(e(j)) = this%c(e(j)) + c(j)
    end do

    if (present(x0)) this%x0 = x0

  end subroutine define_scafun_poly

  pure logical function defined_scafun_poly (this)
    type(scafun_poly), intent(in) :: this
    defined_scafun_poly = .false.
    if (.not.associated(this%c)) return
    if (lbound(this%c,1) /= this%emin) return
    if (ubound(this%c,1) /= this%emax) return
    if (this%emin > 0) return
    if (this%emax < 0) return
    defined_scafun_poly = .true.
  end function defined_scafun_poly

  pure subroutine destroy_scafun_poly (this)
    type(scafun_poly), intent(inout) :: this
    type(scafun_poly) :: default
    if (associated(this%c)) deallocate(this%c)
    this = default  ! assign default initialization values
  end subroutine destroy_scafun_poly

  function eval_scafun_poly(this, x) result(f)

    type(scafun_poly), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8) :: f

    integer  :: j
    real(r8) :: z, w

    !! Polynomial terms with non-negative powers.
    f = this%c(this%emax)
    if (this%emax > 0) then
      z = x(1) - this%x0
      do j = this%emax, 1, -1
        f = this%c(j-1) + z*f
      end do
    end if

    !! Polynomial terms with negative powers.
    if (this%emin < 0) then
      w = this%c(this%emin)
      z = 1.0_r8/(x(1) - this%x0)
      do j = this%emin, -2
        w = this%c(j+1) + z*w
      end do
      f = f + z*w
    end if

  end function eval_scafun_poly

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SPECIFIC PROCEDURES FOR MULTIVARIATE POLYNOMIAL FUNCTIONS
 !!
 !! This is the naive implementation that merely forms a sum of the specified
 !! polynomial terms, making heavy use of the exponentiation operator.
 !!
 !! N.B. - it is universal convention when dealing with polynomial expressions
 !! that x^0 is defined to be 1 regardless of the value of x.  (We are dealing
 !! exclusively with integer exponents in this context.)  Fortran, however,
 !! has the convention that 0.0**0 is undefined, thus we must explicitly check
 !! for 0 exponent values and avoid exponentation in those cases.
 !!
 !! N.B. - We should really be using a multivariate form of Horner's scheme
 !! to evaluate the polynomial.  It should be significantly more numerically
 !! stable and efficient.  I (NNC) have a test implementation, but the results
 !! on these counts are unexpected and it requires further investigation.
 !!
 !! This type has a pointer component, and thus requires a defined assignment.
 !!

  pure subroutine scafun_mpoly_copy (dest, src)

    type(scafun_mpoly), intent(out) :: dest
    type(scafun_mpoly), intent(in)  :: src

    if (associated(src%x0)) then
      allocate(dest%x0(size(src%x0)))
      dest%x0 = src%x0
    end if
    if (associated(src%expon)) then
      allocate(dest%expon(size(src%expon,1),size(src%expon,2)))
      dest%expon = src%expon
    end if
    if (associated(src%coef)) then
      allocate(dest%coef(size(src%coef)))
      dest%coef = src%coef
    end if

  end subroutine scafun_mpoly_copy

  subroutine define_scafun_mpoly (this, coef, expon, x0)

    type(scafun_mpoly), intent(out) :: this
    real(r8), intent(in) :: coef(:)
    integer,  intent(in) :: expon(:,:)
    real(r8), intent(in), optional :: x0(:)

    INSIST( size(coef) > 0 )
    INSIST( size(coef) == size(expon,dim=2) )
    INSIST( size(expon,dim=1) > 1 )

    allocate(this%x0(size(expon,dim=1)))
    if (present(x0)) then
      INSIST( size(x0) == size(this%x0) )
      this%x0 = x0
    else
      this%x0 = 0.0_r8
    end if

    allocate(this%expon(size(expon,1),size(expon,2)))
    this%expon = expon

    allocate(this%coef(size(coef)))
    this%coef = coef

  end subroutine define_scafun_mpoly

  pure logical function defined_scafun_mpoly (this)
    type(scafun_mpoly), intent(in) :: this
    defined_scafun_mpoly = .false.
    if (.not.associated(this%coef)) return
    if (.not.associated(this%expon)) return
    if (size(this%expon,2) /= size(this%coef)) return
    if (.not.associated(this%x0)) return
    if (size(this%expon,1) /= size(this%x0)) return
    if (size(this%x0) < 2) return
    defined_scafun_mpoly = .true.
  end function defined_scafun_mpoly

  pure subroutine destroy_scafun_mpoly (this)
    type(scafun_mpoly), intent(inout) :: this
    type(scafun_mpoly) :: default
    if (associated(this%coef)) deallocate(this%coef)
    if (associated(this%expon)) deallocate(this%expon)
    if (associated(this%x0)) deallocate(this%x0)
    this = default  ! assign default initialization values
  end subroutine destroy_scafun_mpoly

  function eval_scafun_mpoly(this, x) result(f)

    type(scafun_mpoly), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8) :: f

    integer  :: i, j
    real(r8) :: t

    ASSERT( size(x) >= size(this%x0) )

    f = 0.0_r8
    do j = 1, size(this%coef)
      t = this%coef(j)
      do i = 1, size(this%x0)
        if (this%expon(i,j) /= 0) then
          t = t * (x(i)-this%x0(i))**this%expon(i,j)
        end if
      end do
      f = f + t
    end do

  end function eval_scafun_mpoly

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SPECIFIC PROCEDURES FOR SMOOTH STEP FUNCTIONS
 !!

  subroutine define_scafun_sstep (this, x0, y0, x1, y1)

    type(scafun_sstep), intent(out) :: this
    real(r8), intent(in) :: x0, y0, x1, y1

    INSIST( x0 < x1 )
    
    this%x0 = x0
    this%y0 = y0
    this%x1 = x1
    this%y1 = y1

  end subroutine define_scafun_sstep

  pure logical function defined_scafun_sstep (this)
    type(scafun_sstep), intent(in) :: this
    defined_scafun_sstep = .true.
  end function defined_scafun_sstep

  pure subroutine destroy_scafun_sstep (this)
    type(scafun_sstep), intent(inout) :: this
    type(scafun_sstep) :: default
    this = default  ! assign default initialization values
  end subroutine destroy_scafun_sstep

  function eval_scafun_sstep(this, x) result(f)
    type(scafun_sstep), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8) :: f, s
    ASSERT( size(x) >= 1 )
    if (x(1) <= this%x0) then
      f = this%y0
    else if (x(1) >= this%x1) then
      f = this%y1
    else
      s = (x(1) - this%x0) / (this%x1 - this%x0)
      f = this%y0 + (this%y1 - this%y0)*s*s*(3 - 2*s)
    end if
  end function eval_scafun_sstep

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SPECIFIC PROCEDURES FOR THERMO-MECHANICS DENSITY FUNCTIONS
 !!

  pure subroutine scafun_tm_density_copy (dest, src)
    type(scafun_tm_density), intent(out) :: dest
    type(scafun_tm_density), intent(in)  :: src
    dest%d0 = src%d0
    call scafun_copy_aux (dest%int_cte, src%int_cte)
  end subroutine scafun_tm_density_copy

  subroutine define_scafun_tm_density (this, d0, T0, cte, stat, errmsg)
    type(scafun_tm_density), intent(out) :: this
    real(r8), intent(in) :: d0  ! reference density
    real(r8), intent(in) :: T0  ! reference temperature
    type(scafun), intent(in) :: cte ! linear coefficient of thermal expansion
    integer, intent(out) :: stat
    character(*), intent(out) :: errmsg
    this%d0 = d0
    call create_scafun_antideriv (cte, T0, 0.0_r8, this%int_cte, stat, errmsg)
  end subroutine define_scafun_tm_density

  pure logical function defined_scafun_tm_density (this)
    type(scafun_tm_density), intent(in) :: this
    defined_scafun_tm_density = defined_scafun_aux(this%int_cte)
  end function defined_scafun_tm_density

  pure subroutine destroy_scafun_tm_density (this)
    type(scafun_tm_density), intent(inout) :: this
    type(scafun_tm_density) :: default
    call destroy_scafun_aux (this%int_cte)
    this = default  ! assign default initialization values
  end subroutine destroy_scafun_tm_density

  function eval_scafun_tm_density (this, x) result(f)
    type(scafun_tm_density), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8) :: f
    f =  this%d0 * exp(-3*eval(this%int_cte,x))
  end function eval_scafun_tm_density


  elemental logical function is_const_scafun (this)
    type(scafun), intent(in) :: this
    is_const_scafun = (this%dyn_type == TYPE_IS_CONST)
  end function is_const_scafun

  subroutine create_scafun_antideriv (f, x0, g0, g, stat, errmsg)

    type(scafun), intent(in) :: f
    real(r8), intent(in) :: x0, g0
    type(scafun), intent(out) :: g
    integer, intent(out) :: stat
    character(len=*), intent(out) :: errmsg

    integer :: j

    select case (f%dyn_type)
    case (TYPE_IS_CONST)

      call create_scafun_poly (g, c=(/g0,f%const%c/), e=(/0,1/), x0=x0)

    case (TYPE_IS_POLY)

      if (f%poly%emin < 0) then
        if (f%poly%c(-1) /= 0.0_r8) then
          stat = -1
          errmsg = 'unable to create antiderivative for poly with 1/x term'
          return
        end if
      end if

      g%dyn_type = TYPE_IS_POLY
      allocate(g%poly)

      g%poly%emin = min(0, f%poly%emin + 1)
      if (f%poly%emax == 0 .and. f%poly%c(0) == 0.0_r8) then
        g%poly%emax = 0
      else
       g%poly%emax = f%poly%emax + 1
      end if

      allocate(g%poly%c(g%poly%emin:g%poly%emax))

      do j = g%poly%emin, g%poly%emax
        if (j == 0) cycle
        g%poly%c(j) = f%poly%c(j-1)/j
      end do

      g%poly%x0 = f%poly%x0
      g%poly%c(0) = 0.0_r8
      g%poly%c(0) = g0 - eval_scafun_poly(g%poly,(/x0/))

    case default

      stat = -1
      errmsg = 'cannot create antiderivative for this type of function'
      return

    end select

    stat = 0
    errmsg = ''

  end subroutine create_scafun_antideriv

  subroutine create_scafun_product (f, g, fg, stat, errmsg)

    type(scafun), intent(in) :: f, g
    type(scafun), intent(out) :: fg
    integer, intent(out) :: stat
    character(len=*), intent(out) :: errmsg

    integer :: j

    if (is_const_scafun(f)) then
    
      select case (g%dyn_type)
      case (TYPE_IS_CONST)
        fg = g
        fg%const%c = f%const%c * fg%const%c
      case (TYPE_IS_POLY)
        fg = g
        fg%poly%c = f%const%c * fg%poly%c
      case (TYPE_IS_MPOLY)
        fg = g
        fg%mpoly%coef = f%const%c * fg%mpoly%coef
      case default
        stat = -1
        errmsg = 'cannot create product for this type of second function argument'
        return
      end select
      
    else
    
      stat = -1
      errmsg = 'first function argument must be a constant'
      return
      
    end if

    stat = 0
    errmsg = ''

  end subroutine create_scafun_product

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! Procedures for manipulating scafun list variables.
 !!

  subroutine append_to_list0 (head, f)

    type(scafun_list), pointer :: head
    type(scafun), intent(in) :: f

    type(scafun_list), pointer :: l

    if (associated(head)) then
      l => head
      do while (associated(l%next))
        l => l%next
      end do
      allocate(l%next)
      l%next%f = f
    else
      allocate(head)
      head%f = f
    end if

  end subroutine append_to_list0

  subroutine append_to_list1 (head, f)

    type(scafun_vector_list), pointer :: head
    type(scafun), pointer :: f(:)

    type(scafun_vector_list), pointer :: l

    if (associated(head)) then
      l => head
      do while (associated(l%next))
        l => l%next
      end do
      allocate(l%next)
      l%next%f => f
    else
      allocate(head)
      head%f => f
    end if
    f => null()

  end subroutine append_to_list1

  subroutine destroy_scafun_list (head)

    type(scafun_list), pointer :: head

    type(scafun_list), pointer :: first

    do while (associated(head))
      first => head
      head => head%next
      call destroy (first%f)
      deallocate(first)
    end do

  end subroutine destroy_scafun_list

  subroutine destroy_scafun_vector_list (head)

    type(scafun_vector_list), pointer :: head

    type(scafun_vector_list), pointer :: first

    do while (associated(head))
      first => head
      head => head%next
      call destroy (first%f)
      deallocate(first%f)
      deallocate(first)
    end do

  end subroutine destroy_scafun_vector_list

  subroutine convert_list_to_array0 (head, array)

    type(scafun_list), pointer :: head
    type(scafun), intent(out) :: array(:)

    integer :: j
    type(scafun_list), pointer :: first

    do j = 1, size(array)
      INSIST(associated(head))
      array(j) = head%f
      call destroy (head%f)
      first => head
      head => head%next
      deallocate(first)
    end do
    INSIST(.not.associated(head))

  end subroutine convert_list_to_array0

  subroutine convert_list_to_array1 (head, array)

    type(scafun_vector_list), pointer :: head
    type(scafun), intent(out) :: array(:,:)

    integer :: j
    type(scafun_vector_list), pointer :: first

    do j = 1, size(array,dim=2)
      INSIST(associated(head))
      INSIST(associated(head%f))
      INSIST(size(head%f) == size(array,dim=1))
      array(:,j) = head%f(:)
      call destroy (head%f)
      deallocate(head%f)
      first => head
      head => head%next
      deallocate(first)
    end do
    INSIST(.not.associated(head))

  end subroutine convert_list_to_array1

end module scalar_functions
