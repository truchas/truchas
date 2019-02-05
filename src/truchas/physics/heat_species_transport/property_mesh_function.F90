!!
!! PROPERTY_MESH_FUNCTION
!!
!! The derived type PROP_MF is an opaque structure that encapsulates the
!! information required for mesh-wide evaluation of a scalar property on
!! mesh cells as a function of the cell-based state variables.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! No defined assignment is provided; do not use instances in an assignment
!! statement unless you really know what the default assignment is doing.
!!
!!  CALL PMF_CREATE (THIS, MMF, PROPERTY_ID, STAT, ERRMSG) configures the
!!    PROP_MF object THIS to describe property PROPERTY_ID on the mesh given
!!    the mesh-wide distribution of materials specified by the MATL_MESH_FUNC object
!!    MMF.  The integer STAT returns a nonzero value if an error condition
!!    occurred, and an explanatory message is returned in the character string
!!    ERRMSG.  The structure THIS maintains a pointer to the MMF argument and
!!    so if the actual argument is not a pointer it must be given the target
!!    attribute (IMPORTANT!)
!!
!!  CALL PMF_EVAL (THIS, STATE, VALUE) evaluates the property function 
!!    described by the PROP_MF object THIS on all mesh cells given the
!!    cell-based values of the state variables provided by the rank-2
!!    real array STATE.  The computed property values are returned in
!!    the cell-based real array VALUE.  STATE(j,:) holds the values of
!!    the state variables on cell j, and VALUE(j) returns the property
!!    value on cell j.
!!
!!  CALL DESTROY (THIS) deallocates any allocated storage associated with the
!!    PROP_MF object THIS, returning it to its default initialization state.
!!
!! TREATMENT OF VOID
!!
!! In our usage of material mesh functions, the special material ID 0 is used
!! to denote void, or the absence of material.  There is no actual material
!! associated with ID 0, and so no associated property functions to evaluate.
!! Material ID 0 is skipped during evaluation of a property mesh function, so
!! that it makes no contribution to the value; e.g. the returned value for a
!! completely void cell will be 0.  (Note that the material ID list acquired
!! from the material mesh function is sorted, so that void will be the first
!! element in the list if it occurs at all.)
!!

#include "f90_assert.fpp"

module property_mesh_function

  use kinds, only: r8
  use matl_mesh_func_type
  use material_property
  implicit none
  private
  
  public :: pmf_create, pmf_eval, pmf_eval_deriv, destroy, defined
  public :: pmf_set_harmonic_average
  
  type, public :: prop_mf
    private
    type(matl_mesh_func), pointer :: mmf => null()
    type(region), pointer :: reg(:) => null()
    logical :: harmonic_average = .false.
    !! Data facilitating random access evaluation.
    type(cell_prop), pointer :: cprop(:) => null()
  end type prop_mf
  
  type :: region
    integer :: mfirst = 1
    type(mat_prop), pointer :: mp(:) => null()
  end type
  
  type :: cell_prop
    type(mat_prop), pointer :: mprop(:) => null()
    real(r8), pointer :: vfrac(:) => null()
  end type cell_prop
  
  interface pmf_eval
    module procedure pmf_eval, pmf_eval_one
  end interface
  
  interface destroy
    module procedure destroy_prop_mf
  end interface
  
  interface defined
    module procedure defined_prop_mf
  end interface
  
contains

  subroutine pmf_create (this, mmf, property_id, stat, errmsg)
  
    type(prop_mf), intent(out) :: this
    type(matl_mesh_func), intent(in), target :: mmf
    integer, intent(in) :: property_id
    integer, intent(out) :: stat
    character(len=*), intent(out) :: errmsg
    
    integer :: i, m
    integer, pointer :: material_id(:)
    
    this%mmf => mmf
    
    allocate(this%reg(mmf%num_reg()))
    
    do i = 1, mmf%num_reg()
      material_id => mmf%reg_matl(i)
      allocate(this%reg(i)%mp(size(material_id)))
      if (material_id(1) == 0) this%reg(i)%mfirst = 2  ! initial void material -- skip it
      do m = this%reg(i)%mfirst, size(material_id)
        call mp_create (this%reg(i)%mp(m), material_id(m), property_id, stat, errmsg)
        if (stat /= 0) then
          errmsg = 'unable to initialize region: ' // trim(errmsg)
          call destroy (this)
          return
        end if
      end do
    end do
    
    call init_cprop (this)
    
    stat = 0
    errmsg = ''
    
  end subroutine pmf_create
  
  subroutine init_cprop (this)
    use base_mesh_class
    type(prop_mf), intent(inout) :: this
    integer :: i, j, n
    class(base_mesh), pointer :: mesh
    type(mat_prop), pointer :: mprop(:)
    integer, pointer :: cell(:)
    real(r8), pointer :: vfrac(:,:)
    mesh => this%mmf%mesh_ptr()
    allocate(this%cprop(mesh%ncell))
    do i = 1, this%mmf%num_reg()
      cell  => this%mmf%reg_cell(i)
      vfrac => this%mmf%reg_vol_frac(i)
      !! The non-void material properties for this region
      mprop => this%reg(i)%mp(this%reg(i)%mfirst:)
      do j = 1, size(cell)
        n = cell(j)
        this%cprop(n)%mprop => mprop
        if (associated(vfrac)) this%cprop(n)%vfrac => vfrac(j,this%reg(i)%mfirst:)
      end do
    end do
  end subroutine init_cprop
  
  subroutine pmf_set_harmonic_average (this)
    type(prop_mf), intent(inout) :: this
    this%harmonic_average = .true.
  end subroutine pmf_set_harmonic_average
  
  subroutine pmf_eval (this, state, value)
  
    type(prop_mf), intent(in) :: this
    real(r8), intent(in)  :: state(:,:)
    real(r8), intent(out) :: value(:)
    
    integer :: i, j, m, n
    real(r8) :: val
    real(r8), pointer :: vfrac(:,:)
    integer,  pointer :: cell(:)
    type(mat_prop), pointer :: mp
    logical, allocatable :: mask(:)
    
    value = 0.0_r8
    
    do i = 1, this%mmf%num_reg()
    
      vfrac => this%mmf%reg_vol_frac(i)
      cell  => this%mmf%reg_cell(i)
      
      if (size(cell) == 0) cycle  ! nothing to do here
      
      if (associated(vfrac)) then ! multi-material mesh region
      
        if (this%harmonic_average) then
        
          allocate(mask(size(cell)))
          mask = .false.
          do m = this%reg(i)%mfirst, size(this%reg(i)%mp)
            mp => this%reg(i)%mp(m)
            do j = 1, size(cell)
              n = cell(j)
              if (vfrac(j,m) > 0.0_r8) then
                mask(j) = .true.
                call mp_eval (mp, state(n,:), val)
                value(n) = value(n) + vfrac(j,m) / val
              end if
            end do
          end do
          do j = 1, size(cell)
            n = cell(j)
            if (mask(j)) value(n) = 1.0_r8 / value(n)
          end do
          deallocate(mask)
        
        else  ! straight volume fraction averaging
      
          do m = this%reg(i)%mfirst, size(this%reg(i)%mp)
            mp => this%reg(i)%mp(m)
            do j = 1, size(cell)
              n = cell(j)
              if (vfrac(j,m) > 0.0_r8) then
                call mp_eval (mp, state(n,:), val)
                value(n) = value(n) + vfrac(j,m) * val
              end if
            end do
          end do
        
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
      
        if (this%reg(i)%mfirst == 1) then ! non-void
          mp => this%reg(i)%mp(1)
          do j = 1, size(cell)
            n = cell(j)
            call mp_eval(mp, state(n,:), value(n))
          end do
        end if
      
      end if
      
    end do

  end subroutine pmf_eval
  
  
  subroutine pmf_eval_deriv (this, state, index, value)
  
    type(prop_mf), intent(in) :: this
    real(r8), intent(in)  :: state(:,:)
    integer,  intent(in)  :: index
    real(r8), intent(out) :: value(:)
    
    integer :: i, j, m, n
    real(r8) :: val
    real(r8), pointer :: vfrac(:,:)
    integer,  pointer :: cell(:)
    type(mat_prop), pointer :: mp
    logical, allocatable :: mask(:)
    
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
      
          do m = this%reg(i)%mfirst, size(this%reg(i)%mp)
            mp => this%reg(i)%mp(m)
            do j = 1, size(cell)
              n = cell(j)
              if (vfrac(j,m) > 0.0_r8) then
                call mp_eval_deriv (mp, state(n,:), index, val)
                value(n) = value(n) + vfrac(j,m) * val
              end if
            end do
          end do
        
        end if
      
      else  ! single-material mesh region
      
        if (this%reg(i)%mfirst == 1) then ! non-void
          mp => this%reg(i)%mp(1)
          do j = 1, size(cell)
            n = cell(j)
            call mp_eval_deriv (mp, state(n,:), index, value(n))
          end do
        end if
      
      end if
      
    end do

  end subroutine pmf_eval_deriv

  subroutine pmf_eval_one (this, n, state, value)
  
    type(prop_mf), intent(in) :: this
    integer, intent(in) :: n
    real(r8), intent(in)  :: state(:)
    real(r8), intent(out) :: value
    
    integer :: m
    real(r8) :: val
    real(r8), pointer :: vfrac(:)
    type(mat_prop), pointer :: mprop(:)
    
    ASSERT(n >= 1 .and. n <= size(this%cprop))
    
    mprop => this%cprop(n)%mprop
    vfrac => this%cprop(n)%vfrac
    
    value = 0.0_r8
    if (associated(vfrac)) then
      if (this%harmonic_average) then
        !TODO -- implement harmonic averaging.
        INSIST(.false.)
      else
        do m = 1, size(mprop)
          if (vfrac(m) > 0.0_r8) then
            call mp_eval (mprop(m), state, val)
            value = value + vfrac(m) * val
          end if
        end do
      end if
    else  ! single-material cell, but possibly void
      if (size(mprop) > 0) call mp_eval (mprop(1), state, value)
    end if
  
  end subroutine pmf_eval_one

  subroutine destroy_prop_mf (this)
    type(prop_mf), intent(inout) :: this
    integer :: j
    type(prop_mf) :: default
    !! N.B.  The PROP_MF object does not own the target of the MATL_MESH_FUNC pointer
    !! component and so this routine must not deallocate/destroy it.
    if (associated(this%reg)) then
      do j = 1, size(this%reg)
        if (associated(this%reg(j)%mp)) deallocate(this%reg(j)%mp)
      end do
      deallocate(this%reg)
    end if
    if (associated(this%cprop)) deallocate(this%cprop)
    this = default  ! assign default initialization values
  end subroutine destroy_prop_mf
  
  elemental logical function defined_prop_mf (this)
    type(prop_mf), intent(in) :: this
    defined_prop_mf = associated(this%reg)
  end function defined_prop_mf

end module property_mesh_function
