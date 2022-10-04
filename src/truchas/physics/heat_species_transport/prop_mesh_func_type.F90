!!
!! PROP_MESH_FUNC
!!
!! The derived type prop_mesh_func is an opaque structure that encapsulates the
!! information required for mesh-wide evaluation of a scalar property on
!! mesh cells as a function of the cell-based state variables.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module prop_mesh_func_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use matl_mesh_func_type
  use avg_matl_prop_type
  implicit none
  private
  
!  public :: pmf_create, pmf_eval, pmf_eval_deriv, destroy, defined
!  public :: pmf_set_harmonic_average
  
  type, public :: prop_mesh_func
    private
    type(matl_mesh_func), pointer :: mmf => null()  ! reference only -- not owned
    type(region), pointer :: reg(:) => null()
    logical :: harmonic_average = .false.
    !! Data facilitating random access evaluation.
    type(cell_prop), pointer :: cprop(:) => null()
  contains
    procedure :: init
    generic   :: compute_value => pmf_eval, pmf_eval_one
    procedure, private :: pmf_eval, pmf_eval_one
    procedure :: compute_deriv => pmf_eval_deriv
    procedure :: defined => defined_prop_mesh_func
    procedure :: set_harmonic_average => pmf_set_harmonic_average
    !FIXME need finalizer or replace pointer components with allocatables
  end type prop_mesh_func
  
  type :: region
    integer :: mfirst = 1
    integer :: nmatl ! number of non-void materials
    !type(mat_prop), pointer :: mp(:) => null()
    type(avg_matl_prop), allocatable :: prop  !TODO: must this be allocatable?
  end type
  
  type :: cell_prop
    !type(mat_prop), pointer :: mprop(:) => null()
    type(avg_matl_prop), pointer :: mprop => null()
    real(r8), pointer :: vfrac(:) => null()
  end type cell_prop
  
!  interface pmf_eval
!    module procedure pmf_eval, pmf_eval_one
!  end interface
!  
!  interface destroy
!    module procedure destroy_prop_mesh_func
!  end interface
!  
!  interface defined
!    module procedure defined_prop_mesh_func
!  end interface
  
contains

  subroutine init(this, mmf, name, stat, errmsg)
  
    use material_model_driver, only: matl_model

    class(prop_mesh_func), intent(out) :: this
    type(matl_mesh_func), intent(in), target :: mmf
    character(*), intent(in) :: name
    integer, intent(out) :: stat
    character(*), intent(out) :: errmsg !FIXME: allocatable
    
    integer :: i
    integer, pointer :: material_id(:)
    character(:), allocatable :: errmsg2
    
    this%mmf => mmf
    
    allocate(this%reg(mmf%num_reg()))
    
    do i = 1, mmf%num_reg()
      material_id => mmf%reg_matl(i)
      !allocate(this%reg(i)%mp(size(material_id)))
      if (material_id(1) == 0) this%reg(i)%mfirst = 2  ! initial void material -- skip it
      associate (matl_ids => material_id(this%reg(i)%mfirst:))
        if (size(matl_ids) > 0) then
          !call matl_model%alloc_avg_matl_prop(name, matl_ids, this%reg(i)%prop, errmsg2)
          allocate(this%reg(i)%prop)
          call this%reg(i)%prop%init(name, matl_ids, matl_model, stat, errmsg2)
          if (stat /= 0) then
            stat = 1
            errmsg = 'unable to initialize region: ' // errmsg2
            return
          end if
        end if
      end associate
      ! this%reg(i)%prop not allocated => single material region of void
      
      !do m = this%reg(i)%mfirst, size(material_id)
      !  call mp_create (this%reg(i)%mp(m), material_id(m), property_id, stat, errmsg)
      !  if (stat /= 0) then
      !    errmsg = 'unable to initialize region: ' // trim(errmsg)
      !    !call destroy (this)
      !    return
      !  end if
      !end do
    end do
    
    call init_cprop (this)
    
    stat = 0
    errmsg = ''
    
  end subroutine init
  
  subroutine init_cprop(this)
    use base_mesh_class
    type(prop_mesh_func), intent(inout) :: this
    integer :: i, j, n
    class(base_mesh), pointer :: mesh
    integer, pointer :: cell(:)
    real(r8), pointer :: vfrac(:,:)
    class(avg_matl_prop), pointer :: mprop
    mesh => this%mmf%mesh_ptr()
    allocate(this%cprop(mesh%ncell))
    do i = 1, this%mmf%num_reg()
      cell  => this%mmf%reg_cell(i)
      vfrac => this%mmf%reg_vol_frac(i)
      !! The non-void material properties for this region
      !mprop => this%reg(i)%mp(this%reg(i)%mfirst:)
      mprop => null()
      if (allocated(this%reg(i)%prop)) mprop => this%reg(i)%prop
      do j = 1, size(cell)
        n = cell(j)
        this%cprop(n)%mprop => mprop
        if (associated(vfrac)) this%cprop(n)%vfrac => vfrac(j,this%reg(i)%mfirst:)
      end do
    end do
  end subroutine init_cprop
  
  subroutine pmf_set_harmonic_average (this)
    class(prop_mesh_func), intent(inout) :: this
    this%harmonic_average = .true.
  end subroutine pmf_set_harmonic_average
  
  subroutine pmf_eval (this, state, value)
  
    class(prop_mesh_func), intent(in) :: this
    real(r8), intent(in)  :: state(:,:)
    real(r8), intent(out) :: value(:)
    
    integer :: i, j, n
    real(r8), pointer :: vfrac(:,:)
    integer,  pointer :: cell(:)
    !type(mat_prop), pointer :: mp
    
    value = 0.0_r8
    
    do i = 1, this%mmf%num_reg()
    
      vfrac => this%mmf%reg_vol_frac(i)
      cell  => this%mmf%reg_cell(i)
      
      if (size(cell) == 0) cycle  ! nothing to do here
      
      if (associated(vfrac)) then ! multi-material mesh region
      
        if (this%harmonic_average) then
        
          INSIST(.false.) ! not implemented
          
!          allocate(mask(size(cell)))
!          mask = .false.
!          do m = this%reg(i)%mfirst, size(this%reg(i)%mp)
!            mp => this%reg(i)%mp(m)
!            do j = 1, size(cell)
!              n = cell(j)
!              if (vfrac(j,m) > 0.0_r8) then
!                mask(j) = .true.
!                call mp_eval (mp, state(n,:), val)
!                value(n) = value(n) + vfrac(j,m) / val
!              end if
!            end do
!          end do
!          do j = 1, size(cell)
!            n = cell(j)
!            if (mask(j)) value(n) = 1.0_r8 / value(n)
!          end do
!          deallocate(mask)
        
        else  ! straight volume fraction averaging
      
          do j = 1, size(cell)
            n = cell(j)
            call this%reg(i)%prop%compute_value(vfrac(j,this%reg(i)%mfirst:), state(n,:), value(n))
          end do
          
!          do m = this%reg(i)%mfirst, size(this%reg(i)%mp)
!            mp => this%reg(i)%mp(m)
!            do j = 1, size(cell)
!              n = cell(j)
!              if (vfrac(j,m) > 0.0_r8) then
!                call mp_eval (mp, state(n,:), val)
!                value(n) = value(n) + vfrac(j,m) * val
!              end if
!            end do
!          end do
        
        end if
      
       ! do j = 1, size(cell)
       !   n = cell(j)
       !   tmp = 0.0_r8
       !   do m = this%reg(i)%mfirst, size(this%reg(i)%mp)
       !     mp => this%reg(i)%mp(m)
       !     if (vfrac(j,m) > 0.0_r8) then
       !       call mp_eval (mp, state(n,:), val)
       !       tmp = tmp + vfrac(j,m) * val
       !     end if
       !   end do
       !   value(n) = tmp
       ! end do
      
      else  ! single-material mesh region
      
        !FIXME: Need a single-material version of compute_value
        if (this%reg(i)%mfirst == 1) then ! non-void
          do j = 1, size(cell)
            n = cell(j)
            call this%reg(i)%prop%compute_value([1.0_r8], state(n,:), value(n))
          end do
        end if
          
!        if (this%reg(i)%mfirst == 1) then ! non-void
!          mp => this%reg(i)%mp(1)
!          do j = 1, size(cell)
!            n = cell(j)
!            call mp_eval(mp, state(n,:), value(n))
!          end do
!        end if
      
      end if
      
    end do

  end subroutine pmf_eval
  
  
  subroutine pmf_eval_deriv (this, state, index, value)
  
    class(prop_mesh_func), intent(in) :: this
    real(r8), intent(in)  :: state(:,:)
    integer,  intent(in)  :: index
    real(r8), intent(out) :: value(:)
    
    integer :: i, j, n
    real(r8), pointer :: vfrac(:,:)
    integer,  pointer :: cell(:)
    !type(mat_prop), pointer :: mp
    
    value = 0.0_r8
    
    do i = 1, this%mmf%num_reg()
    
      vfrac => this%mmf%reg_vol_frac(i)
      cell  => this%mmf%reg_cell(i)
      
      if (size(cell) == 0) cycle  ! nothing to do here
      
      if (associated(vfrac)) then ! multi-material mesh region
      
        if (this%harmonic_average) then
        
          !! Not yet implemented
          INSIST(.false.)
        
        else  ! straight volume fraction averaging
      
          do j = 1, size(cell)
            n = cell(j)
            call this%reg(i)%prop%compute_deriv(vfrac(j,this%reg(i)%mfirst:), state(n,:), index, value(n))
          end do

!          do m = this%reg(i)%mfirst, size(this%reg(i)%mp)
!            mp => this%reg(i)%mp(m)
!            do j = 1, size(cell)
!              n = cell(j)
!              if (vfrac(j,m) > 0.0_r8) then
!                call mp_eval_deriv (mp, state(n,:), index, val)
!                value(n) = value(n) + vfrac(j,m) * val
!              end if
!            end do
!          end do
        
        end if
      
      else  ! single-material mesh region
      
        !FIXME: Need a single-material version of compute_value
        if (this%reg(i)%mfirst == 1) then ! non-void
          do j = 1, size(cell)
            n = cell(j)
            call this%reg(i)%prop%compute_deriv([1.0_r8], state(n,:), index, value(n))
          end do
        end if
      
!        if (this%reg(i)%mfirst == 1) then ! non-void
!          mp => this%reg(i)%mp(1)
!          do j = 1, size(cell)
!            n = cell(j)
!            call mp_eval_deriv (mp, state(n,:), index, value(n))
!          end do
!        end if
      
      end if
      
    end do

  end subroutine pmf_eval_deriv

  subroutine pmf_eval_one (this, n, state, value)
  
    class(prop_mesh_func), intent(in) :: this
    integer, intent(in) :: n
    real(r8), intent(in)  :: state(:)
    real(r8), intent(out) :: value
    
    real(r8), pointer :: vfrac(:)
    !type(mat_prop), pointer :: mprop(:)
    class(avg_matl_prop), pointer :: mprop
    
    ASSERT(n >= 1 .and. n <= size(this%cprop))
    
    mprop => this%cprop(n)%mprop
    vfrac => this%cprop(n)%vfrac
    
    value = 0.0_r8
    if (associated(vfrac)) then
      if (this%harmonic_average) then
        !TODO -- implement harmonic averaging.
        INSIST(.false.)
      else
        call this%cprop(n)%mprop%compute_value(vfrac, state, value)
        !do m = 1, size(mprop)
        !  if (vfrac(m) > 0.0_r8) then
        !    call mp_eval (mprop(m), state, val)
        !    value = value + vfrac(m) * val
        !  end if
        !end do
      end if
    else  ! single-material cell, but possibly void
      if (associated(mprop)) call mprop%compute_value([1.0_r8], state, value)
      !if (size(mprop) > 0) call mp_eval (mprop(1), state, value)
    end if
  
  end subroutine pmf_eval_one

!  subroutine destroy_prop_mesh_func (this)
!    type(prop_mesh_func), intent(inout) :: this
!    integer :: j
!    type(prop_mesh_func) :: default
!    !! N.B.  The prop_mesh_func object does not own the target of the MATL_MESH_FUNC pointer
!    !! component and so this routine must not deallocate/destroy it.
!    if (associated(this%reg)) then
!      do j = 1, size(this%reg)
!        if (associated(this%reg(j)%mp)) deallocate(this%reg(j)%mp)
!      end do
!      deallocate(this%reg)
!    end if
!    if (associated(this%cprop)) deallocate(this%cprop)
!    this = default  ! assign default initialization values
!  end subroutine destroy_prop_mesh_func
  
  elemental logical function defined_prop_mesh_func (this)
    class(prop_mesh_func), intent(in) :: this
    defined_prop_mesh_func = associated(this%reg)
  end function defined_prop_mesh_func

end module prop_mesh_func_type
