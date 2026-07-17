!!
!! PROP_CELL_FUNC
!!
!! The derived type prop_cell_func is an opaque structure that encapsulates the
!! information required for mesh-wide evaluation of a scalar property on
!! mesh cells as a function of the cell-based state variables.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module prop_cell_func_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use cell_prop_func_class
  use matl_mesh_func_type
  use avg_matl_prop_type
  implicit none
  private
  
  type, extends(cell_prop_func), public :: prop_cell_func
    private
    type(matl_mesh_func), pointer :: mmf => null()  ! reference only -- not owned
    type(region), allocatable :: reg(:)
    logical :: harmonic_average = .false.
    !! Data facilitating random access evaluation.
    type(cell_prop), allocatable :: cprop(:)
  contains
    procedure :: init
    procedure :: compute_value_t => pmf_eval_t
    procedure :: compute_value_c => pmf_eval
    procedure :: compute_value_tc => pmf_eval_tc
    procedure :: compute_value_cell => pmf_eval_one
    procedure :: compute_value_cell_t => pmf_eval_one_t
    procedure :: compute_value_cell_tc => pmf_eval_one_tc
    procedure, private :: pmf_eval, pmf_eval_one, pmf_eval_one_t, pmf_eval_one_tc, &
        pmf_eval_t, pmf_eval_tc
    procedure :: compute_deriv_t => pmf_eval_deriv_t
    procedure :: compute_deriv_c => pmf_eval_deriv
    procedure :: compute_deriv_tc => pmf_eval_deriv_tc
    procedure, private :: pmf_eval_deriv, pmf_eval_deriv_t, pmf_eval_deriv_tc
    procedure :: defined => defined_prop_cell_func
    procedure :: set_harmonic_average => pmf_set_harmonic_average
  end type prop_cell_func
  
  type :: region
    integer :: mfirst = 1
    type(avg_matl_prop), allocatable :: prop
  end type
  
  type :: cell_prop
    integer :: ireg = 0
    real(r8), pointer :: vfrac(:) => null() ! reference only -- not owned
  end type cell_prop
  
contains

  subroutine init(this, mmf, name, stat, errmsg)
  
    use material_model_driver, only: matl_model

    class(prop_cell_func), intent(out) :: this
    type(matl_mesh_func), intent(in), target :: mmf
    character(*), intent(in) :: name
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    
    integer :: i
    integer, pointer :: material_id(:)
    character(:), allocatable :: errmsg2
    
    this%mmf => mmf
    
    allocate(this%reg(mmf%num_reg()))
    
    do i = 1, mmf%num_reg()
      material_id => mmf%reg_matl(i)
      if (material_id(1) == 0) this%reg(i)%mfirst = 2  ! initial void material -- skip it
      associate (matl_ids => material_id(this%reg(i)%mfirst:))
        if (size(matl_ids) > 0) then
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
    end do
    
    call init_cprop (this)
    
    stat = 0
    errmsg = ''
    
  end subroutine init
  
  subroutine init_cprop(this)
    use base_mesh_class
    type(prop_cell_func), intent(inout) :: this
    integer :: i, j, n
    class(base_mesh), pointer :: mesh
    integer, pointer :: cell(:)
    real(r8), pointer :: vfrac(:,:)
    mesh => this%mmf%mesh_ptr()
    allocate(this%cprop(mesh%ncell))
    do i = 1, this%mmf%num_reg()
      cell  => this%mmf%reg_cell(i)
      vfrac => this%mmf%reg_vol_frac(i)
      do j = 1, size(cell)
        n = cell(j)
        this%cprop(n)%ireg = i
        if (associated(vfrac)) this%cprop(n)%vfrac => vfrac(j,this%reg(i)%mfirst:)
      end do
    end do
  end subroutine init_cprop
  
  subroutine pmf_set_harmonic_average (this)
    class(prop_cell_func), intent(inout) :: this
    this%harmonic_average = .true.
  end subroutine pmf_set_harmonic_average
  
  subroutine pmf_eval_t (this, temp, value)

    class(prop_cell_func), intent(in) :: this
    real(r8), intent(in)  :: temp(:)
    real(r8), intent(out) :: value(:)

    integer :: i, j, n
    real(r8) :: state(1)
    real(r8), pointer :: vfrac(:,:)
    integer,  pointer :: cell(:)

    ASSERT(size(temp) == size(value))

    value = 0.0_r8

    do i = 1, this%mmf%num_reg()

      vfrac => this%mmf%reg_vol_frac(i)
      cell  => this%mmf%reg_cell(i)

      if (size(cell) == 0) cycle  ! nothing to do here

      if (associated(vfrac)) then ! multi-material mesh region

        if (this%harmonic_average) then

          INSIST(.false.) ! not implemented

        else  ! straight volume fraction averaging

          do j = 1, size(cell)
            n = cell(j)
            if (n > size(value)) cycle
            state(1) = temp(n)
            call this%reg(i)%prop%compute_value(vfrac(j,this%reg(i)%mfirst:), state, value(n))
          end do

        end if

      else  ! single-material mesh region

        if (this%reg(i)%mfirst == 1) then ! non-void
          do j = 1, size(cell)
            n = cell(j)
            if (n > size(value)) cycle
            state(1) = temp(n)
            call this%reg(i)%prop%compute_value([1.0_r8], state, value(n))
          end do
        end if

      end if

    end do

  end subroutine pmf_eval_t

  subroutine pmf_eval_tc (this, temp, conc, value)

    class(prop_cell_func), intent(in) :: this
    real(r8), intent(in)  :: temp(:)
    real(r8), intent(in)  :: conc(:,:)
    real(r8), intent(out) :: value(:)

    integer :: i, j, n
    real(r8) :: state(1+size(conc,dim=2))
    real(r8), pointer :: vfrac(:,:)
    integer,  pointer :: cell(:)

    ASSERT(size(temp) == size(value))
    ASSERT(size(conc,dim=1) == size(value))

    value = 0.0_r8

    do i = 1, this%mmf%num_reg()

      vfrac => this%mmf%reg_vol_frac(i)
      cell  => this%mmf%reg_cell(i)

      if (size(cell) == 0) cycle  ! nothing to do here

      if (associated(vfrac)) then ! multi-material mesh region

        if (this%harmonic_average) then

          INSIST(.false.) ! not implemented

        else  ! straight volume fraction averaging

          do j = 1, size(cell)
            n = cell(j)
            if (n > size(value)) cycle
            state(1) = temp(n)
            state(2:) = conc(n,:)
            call this%reg(i)%prop%compute_value(vfrac(j,this%reg(i)%mfirst:), state, value(n))
          end do

        end if

      else  ! single-material mesh region

        if (this%reg(i)%mfirst == 1) then ! non-void
          do j = 1, size(cell)
            n = cell(j)
            if (n > size(value)) cycle
            state(1) = temp(n)
            state(2:) = conc(n,:)
            call this%reg(i)%prop%compute_value([1.0_r8], state, value(n))
          end do
        end if

      end if

    end do

  end subroutine pmf_eval_tc

  subroutine pmf_eval (this, state, value)
  
    class(prop_cell_func), intent(in) :: this
    real(r8), intent(in)  :: state(:,:)
    real(r8), intent(out) :: value(:)
    
    integer :: i, j, n
    real(r8), pointer :: vfrac(:,:)
    integer,  pointer :: cell(:)
    
    value = 0.0_r8
    
    do i = 1, this%mmf%num_reg()
    
      vfrac => this%mmf%reg_vol_frac(i)
      cell  => this%mmf%reg_cell(i)
      
      if (size(cell) == 0) cycle  ! nothing to do here
      
      if (associated(vfrac)) then ! multi-material mesh region
      
        if (this%harmonic_average) then
        
          INSIST(.false.) ! not implemented
        
        else  ! straight volume fraction averaging
      
          do j = 1, size(cell)
            n = cell(j)
            if (n > size(value)) cycle
            call this%reg(i)%prop%compute_value(vfrac(j,this%reg(i)%mfirst:), state(n,:), value(n))
          end do
        
        end if
      
      else  ! single-material mesh region
      
        if (this%reg(i)%mfirst == 1) then ! non-void
          do j = 1, size(cell)
            n = cell(j)
            if (n > size(value)) cycle
            call this%reg(i)%prop%compute_value([1.0_r8], state(n,:), value(n))
          end do
        end if
      
      end if
      
    end do

  end subroutine pmf_eval
  
  
  subroutine pmf_eval_deriv (this, state, index, value)
  
    class(prop_cell_func), intent(in) :: this
    real(r8), intent(in)  :: state(:,:)
    integer,  intent(in)  :: index
    real(r8), intent(out) :: value(:)
    
    integer :: i, j, n
    real(r8), pointer :: vfrac(:,:)
    integer,  pointer :: cell(:)
    
    value = 0.0_r8
    
    do i = 1, this%mmf%num_reg()
    
      vfrac => this%mmf%reg_vol_frac(i)
      cell  => this%mmf%reg_cell(i)
      
      if (size(cell) == 0) cycle  ! nothing to do here
      
      if (associated(vfrac)) then ! multi-material mesh region
      
        if (this%harmonic_average) then
        
          INSIST(.false.) ! not implemented
        
        else  ! straight volume fraction averaging
      
          do j = 1, size(cell)
            n = cell(j)
            if (n > size(value)) cycle
            call this%reg(i)%prop%compute_deriv(vfrac(j,this%reg(i)%mfirst:), state(n,:), index, value(n))
          end do
        
        end if
      
      else  ! single-material mesh region
      
        if (this%reg(i)%mfirst == 1) then ! non-void
          do j = 1, size(cell)
            n = cell(j)
            if (n > size(value)) cycle
            call this%reg(i)%prop%compute_deriv([1.0_r8], state(n,:), index, value(n))
          end do
        end if
      
      end if
      
    end do

  end subroutine pmf_eval_deriv

  subroutine pmf_eval_deriv_t (this, temp, index, value)

    class(prop_cell_func), intent(in) :: this
    real(r8), intent(in)  :: temp(:)
    integer,  intent(in)  :: index
    real(r8), intent(out) :: value(:)

    integer :: i, j, n
    real(r8) :: state(1)
    real(r8), pointer :: vfrac(:,:)
    integer,  pointer :: cell(:)

    ASSERT(size(temp) == size(value))

    value = 0.0_r8

    do i = 1, this%mmf%num_reg()

      vfrac => this%mmf%reg_vol_frac(i)
      cell  => this%mmf%reg_cell(i)

      if (size(cell) == 0) cycle  ! nothing to do here

      if (associated(vfrac)) then ! multi-material mesh region

        if (this%harmonic_average) then

          INSIST(.false.) ! not implemented

        else  ! straight volume fraction averaging

          do j = 1, size(cell)
            n = cell(j)
            if (n > size(value)) cycle
            state(1) = temp(n)
            call this%reg(i)%prop%compute_deriv(vfrac(j,this%reg(i)%mfirst:), state, index, value(n))
          end do

        end if

      else  ! single-material mesh region

        if (this%reg(i)%mfirst == 1) then ! non-void
          do j = 1, size(cell)
            n = cell(j)
            if (n > size(value)) cycle
            state(1) = temp(n)
            call this%reg(i)%prop%compute_deriv([1.0_r8], state, index, value(n))
          end do
        end if

      end if

    end do

  end subroutine pmf_eval_deriv_t

  subroutine pmf_eval_deriv_tc (this, temp, conc, index, value)

    class(prop_cell_func), intent(in) :: this
    real(r8), intent(in)  :: temp(:)
    real(r8), intent(in)  :: conc(:,:)
    integer,  intent(in)  :: index
    real(r8), intent(out) :: value(:)

    integer :: i, j, n
    real(r8) :: state(1+size(conc,dim=2))
    real(r8), pointer :: vfrac(:,:)
    integer,  pointer :: cell(:)

    ASSERT(size(temp) == size(value))
    ASSERT(size(conc,dim=1) == size(value))

    value = 0.0_r8

    do i = 1, this%mmf%num_reg()

      vfrac => this%mmf%reg_vol_frac(i)
      cell  => this%mmf%reg_cell(i)

      if (size(cell) == 0) cycle  ! nothing to do here

      if (associated(vfrac)) then ! multi-material mesh region

        if (this%harmonic_average) then

          INSIST(.false.) ! not implemented

        else  ! straight volume fraction averaging

          do j = 1, size(cell)
            n = cell(j)
            if (n > size(value)) cycle
            state(1) = temp(n)
            state(2:) = conc(n,:)
            call this%reg(i)%prop%compute_deriv(vfrac(j,this%reg(i)%mfirst:), state, index, value(n))
          end do

        end if

      else  ! single-material mesh region

        if (this%reg(i)%mfirst == 1) then ! non-void
          do j = 1, size(cell)
            n = cell(j)
            if (n > size(value)) cycle
            state(1) = temp(n)
            state(2:) = conc(n,:)
            call this%reg(i)%prop%compute_deriv([1.0_r8], state, index, value(n))
          end do
        end if

      end if

    end do

  end subroutine pmf_eval_deriv_tc

  subroutine pmf_eval_one (this, n, state, value)
  
    class(prop_cell_func), intent(in) :: this
    integer, intent(in) :: n
    real(r8), intent(in)  :: state(:)
    real(r8), intent(out) :: value
    
    real(r8), pointer :: vfrac(:)
    integer :: ireg
    
    ASSERT(n >= 1 .and. n <= size(this%cprop))
    
    ireg = this%cprop(n)%ireg
    vfrac => this%cprop(n)%vfrac
    
    value = 0.0_r8
    if (associated(vfrac)) then
      if (this%harmonic_average) then
        INSIST(.false.) ! not implemented
      else
        call this%reg(ireg)%prop%compute_value(vfrac, state, value)
      end if
    else  ! single-material cell, but possibly void
      if (allocated(this%reg(ireg)%prop)) call this%reg(ireg)%prop%compute_value([1.0_r8], state, value)
    end if
  
  end subroutine pmf_eval_one

  subroutine pmf_eval_one_t (this, n, temp, value)

    class(prop_cell_func), intent(in) :: this
    integer, intent(in) :: n
    real(r8), intent(in)  :: temp
    real(r8), intent(out) :: value

    real(r8) :: state(1)

    state(1) = temp
    call pmf_eval_one(this, n, state, value)

  end subroutine pmf_eval_one_t

  subroutine pmf_eval_one_tc (this, n, temp, conc, value)

    class(prop_cell_func), intent(in) :: this
    integer, intent(in) :: n
    real(r8), intent(in)  :: temp
    real(r8), intent(in)  :: conc(:)
    real(r8), intent(out) :: value

    real(r8) :: state(1+size(conc))

    state(1) = temp
    state(2:) = conc
    call pmf_eval_one(this, n, state, value)

  end subroutine pmf_eval_one_tc
  
  elemental logical function defined_prop_cell_func (this)
    class(prop_cell_func), intent(in) :: this
    defined_prop_cell_func = allocated(this%reg)
  end function defined_prop_cell_func

end module prop_cell_func_type
