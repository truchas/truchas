module region_factory

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use base_mesh_class
  use parameter_list_type
  use region_class

  public :: alloc_region

contains

  subroutine alloc_region(reg, mesh, params, stat, errmsg)
  
    class(region), allocatable, intent(out) :: reg
    class(base_mesh), intent(in) :: mesh
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    character(:), allocatable :: context, rtype

    context = 'processing ' // params%path() // ': '

    call params%get('type', rtype, stat, errmsg)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if
    
    select case (rtype)
    case ('cell-set')
      block
        use cell_set_region_type, only: alloc_cell_set_region
        integer, allocatable :: setids(:)
        integer :: bitmask
        logical :: c
        call params%get('cell-set-ids', setids, stat, errmsg)
        if (stat /= 0) then
          errmsg = context // errmsg
          return
        end if
        call mesh%get_cell_set_bitmask(setids, bitmask, stat, errmsg)
        if (stat /= 0) then
          errmsg = context // 'invalid "cell-set-ids" value: ' // errmsg
          return
        end if
        call params%get('complement', c, stat, errmsg, default=.false.)
        if (stat /= 0) then
          errmsg = context // errmsg
          return
        end if
        call alloc_cell_set_region(reg, bitmask, c)
      end block
      
    case ('half-plane')
      block
        use half_plane_region_type, only: alloc_half_plane_region
        real(r8), allocatable :: p(:), n(:)
        logical :: c
        call params%get('point', p, stat, errmsg)
        if (stat /= 0) then
          errmsg = context // errmsg
          return
        else if (size(p) /= 2) then
          stat = 1
          errmsg = context // 'a 2-vector for "point" is required'
          return
        end if
        call params%get('normal', n, stat, errmsg)
        if (stat /= 0) then
          errmsg = context // errmsg
          return
        else if (size(n) /= 2) then
          stat = 1
          errmsg = context // 'a 2-vector for "normal" is required'
          return
        else if (norm2(n) == 0) then
          stat = 1
          errmsg = context // '"normal" must be non-zero'
          return
        end if
        call params%get('complement', c, stat, errmsg, default=.false.)
        if (stat /= 0) then
          errmsg = context // errmsg
          return
        end if
        call alloc_half_plane_region(reg, p, n, c)
      end block
      
    case ('box')
      block
        use box_region_type, only: alloc_box_region
        real(r8), allocatable :: x(:), y(:)
        logical :: c
        call params%get('lower-corner', x, stat, errmsg)
        if (stat /= 0) then
          errmsg = context // errmsg
          return
        else if (size(x) /= 2) then
          stat = 1
          errmsg = context // 'a 2-vector for "lower-box-corner" is required'
          return
        end if
        call params%get('upper-corner', y, stat, errmsg)
        if (stat /= 0) then
          errmsg = context // errmsg
          return
        else if (size(y) /= 2) then
          stat = 1
          errmsg = context // 'a 2-vector for "upper-box-corner" is required'
          return
        end if
        call params%get('complement', c, stat, errmsg, default=.false.)
        if (stat /= 0) then
          errmsg = context // errmsg
          return
        end if
        call alloc_box_region(reg, x, y, c)
      end block

    case ('disk')
      block
        use disk_region_type, only: alloc_disk_region
        real(r8), allocatable :: x(:)
        real(r8) :: r
        logical :: c
        call params%get('center', x, stat, errmsg)
        if (stat /= 0) then
          errmsg = context // errmsg
          return
        else if (size(x) /= 2) then
          stat = 1
          errmsg = context // 'a 2-vector for "center" is required'
          return
        end if
        call params%get('radius', r, stat, errmsg)
        if (stat /= 0) then
          errmsg = context // errmsg
          return
        else if (r <= 0) then
          stat = 1
          errmsg = context // '"radius" must be > 0'
          return
        end if
        call params%get('complement', c, stat, errmsg, default=.false.)
        if (stat /= 0) then
          errmsg = context // errmsg
          return
        end if
        call alloc_disk_region(reg, x, r, c)
      end block
    
    case ('background')
      call alloc_background_region(reg)

    case default
      stat = 1
      errmsg = context // 'unknown "type" value: ' // rtype
      return
    end select
  
  end subroutine alloc_region

end module region_factory
