!!
!! INTEGER_MAP_TYPE
!!
!! This module defines an integer_map container that stores integer (key, value)
!! value pairs. It is implemented using a weak AVL (WAVL) binary search tree
!! with bottom up rebalancing.
!!
!! Haeupler, B., Sen, S., and Tarjan, S. E. 2015. Rank-Balanced Trees.
!! ACM Trans. Algorithms 11, 4. DOI:https://doi.org/10.1145/2689412
!!
!! Copyright 2022 Neil N. Carlson <neil.n.carlson@gmail.com>
!! Use subject to the MIT license: https://opensource.org/licenses/MIT
!!

module integer_map_type

  implicit none
  private

  type, public :: integer_map
    private
    type(rbt_node), pointer :: root => null()
  contains
    procedure :: set => map_set
    procedure :: val => map_val
    procedure :: remove => map_remove
    procedure :: contains => map_contains
    final :: integer_map_delete
    !! Primarily testing procedures
  end type

  type :: rbt_node
    integer :: key, val
    integer :: rank = 0
    type(rbt_node), pointer :: left => null(), right => null()
  contains
    procedure :: left_rank_diff
    procedure :: right_rank_diff
    procedure :: is_leaf
  end type

contains

!!!! integer_map PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Final subroutine for integer_map objects
  elemental subroutine integer_map_delete(this)
    type(integer_map), intent(inout) :: this
    call dealloc_rbt_node(this%root)
  contains
    pure recursive subroutine dealloc_rbt_node(node)
      type(rbt_node), pointer :: node
      if (associated(node)) then
        call dealloc_rbt_node(node%left)
        call dealloc_rbt_node(node%right)
        deallocate(node)
      end if
    end subroutine
  end subroutine

  !! Returns true if the map contains an element for the given key
  logical function map_contains(this, key)
    class(integer_map), intent(in) :: this
    integer, intent(in) :: key
    type(rbt_node), pointer :: r
    r => this%root
    map_contains = .true.
    do while (associated(r))
      if (key < r%key) then
        r => r%left
      else if (key > r%key) then
        r => r%right
      else
        return
      end if
    end do
    map_contains = .false.
  end function

  !! Set the value for the given KEY in the map to VAL.
#ifdef GNU_PR69563
  subroutine map_set(this, key, val)
#else
  elemental subroutine map_set(this, key, val)
#endif
    class(integer_map), intent(inout) :: this
    integer, intent(in) :: key, val
    call rbt_insert(this%root, key, val)
  end subroutine

  !! Insert the given (key, val) element into the subtree at ROOT. If an
  !! element with the given key already exists, overwrite its value with
  !! the given val. The ROOT pointer may be modified as a result of the
  !! insertion and rebalancing of the tree.

  pure recursive subroutine rbt_insert(root, key, val)
    type(rbt_node), pointer, intent(inout) :: root
    integer, intent(in) :: key, val
    if (.not.associated(root)) then
      allocate(root)
      root%key = key
      root%val = val
    else if (key < root%key) then
      call rbt_insert(root%left, key, val)  ! output root%left may be a 0-child
      if (root%left_rank_diff() == 0) call insert_left_rebalance(root)
    else if (key > root%key) then
      call rbt_insert(root%right, key, val) ! output root%right may be a 0-child
      if (root%right_rank_diff() == 0) call insert_right_rebalance(root)
    else  ! key == root%key
      root%val = val
    end if
  end subroutine
  
  integer function map_val(this, key)
    class(integer_map), intent(in) :: this
    integer, intent(in) :: key
    type(rbt_node), pointer :: r
    r => this%root
    do while (associated(r))
      if (key < r%key) then
        r => r%left
      else if (key > r%key) then
        r => r%right
      else  ! key == r%key
        map_val = r%val
        return
      end if
    end do
    block ! no element with given key
      use,intrinsic :: iso_fortran_env, only: error_unit
      write(error_unit,'(a)') 'integer_map: out-of-bounds'
      stop 1
    end block
  end function

  !! WAVL procedure for rebalancing the subtree at ROOT after an insertion
  !! in its left subtree resulted in a left 0-child (assumed).

  pure subroutine insert_left_rebalance(root)
    type(rbt_node), pointer, intent(inout) :: root
    if (root%right_rank_diff() == 1) then ! root is 0,1
      root%rank = root%rank + 1 ! makes root 1,2
    else ! root is 0,2
      if (root%left%right_rank_diff() == 2) then
        call rotate_right(root) ! new root is 1,1
        root%right%rank = root%right%rank - 1
      else
        call rotate_left(root%left)
        call rotate_right(root) ! new root is 1,1
        root%rank = root%rank + 1
        root%left%rank = root%left%rank - 1
        root%right%rank = root%right%rank - 1
      end if
    end if
  end subroutine

  !! WAVL procedure for rebalancing the subtree at ROOT after an insertion
  !! in its right subtree resulted in a right 0-child (assumed).

  pure subroutine insert_right_rebalance(root)
    type(rbt_node), pointer, intent(inout) :: root
    if (root%left_rank_diff() == 1) then ! root is 1,0
      root%rank = root%rank + 1 ! makes root 2,1
    else ! root is 2,0
      if (root%right%left_rank_diff() == 2) then
        call rotate_left(root)
        root%left%rank = root%left%rank - 1
      else
        call rotate_right(root%right)
        call rotate_left(root)
        root%rank = root%rank + 1
        root%left%rank = root%left%rank - 1
        root%right%rank = root%right%rank - 1
      end if
    end if
  end subroutine

  !! Remove element with given key from the map
  elemental subroutine map_remove(this, key)
    class(integer_map), intent(inout) :: this
    integer, intent(in) :: key
    call rbt_remove(this%root, key)
  end subroutine

  !! Remove element with the given KEY from the subtree at ROOT. The ROOT
  !! pointer may be modified as a result of the removal and rebalancing
  !! of the tree.

  pure recursive subroutine rbt_remove(root, key)
    type(rbt_node), pointer, intent(inout) :: root
    integer, intent(in) :: key
    type(rbt_node), pointer :: next, temp
    if (.not.associated(root)) return
    if (key < root%key) then
      call rbt_remove(root%left, key)
      call remove_left_rebalance(root)
    else if (key > root%key) then
      call rbt_remove(root%right, key)
      call remove_right_rebalance(root)
    else  ! key == root%key
      if (associated(root%left) .and. associated(root%right)) then
        !! Replace root element with next larger element and then
        !! delete that element from the right subtree.
        next => root%right
        do while (associated(next%left))
          next => next%left
        end do
        root%key = next%key
        root%val = next%val
        call rbt_remove(root%right, next%key)
        call remove_right_rebalance(root)
      else if (associated(root%left)) then
        !! Replace root node with its left child; new root is a 2 or 3-child
        temp => root
        root => root%left
        deallocate(temp)
      else if (associated(root%right)) then
        !! Replace root node with its right child; new root is a 2 or 3-child
        temp => root
        root => root%right
        deallocate(temp)
      else  ! root is a leaf; just delete it
        deallocate(root)
      end if
    end if
  end subroutine

  !! WAVL procedure for rebalancing the subtree at ROOT after a removal
  !! in its left subtree resulted in it becoming a 2,2 leaf or with a
  !! possible left 3-child.

  pure subroutine remove_left_rebalance(root)
    type(rbt_node), pointer, intent(inout) :: root
    if (root%is_leaf()) then
      root%rank = 0 ! makes root 1,1
    else if (root%left_rank_diff() == 3) then
      if (root%right_rank_diff() == 2) then ! root is 3,2
        root%rank = root%rank - 1 ! makes root 2,1 and a 2 or 3-child
      else if (root%right%right_rank_diff() == 1) then ! 3,1 root with 1,1 or 2,1 right child
        call rotate_left(root)  ! new root is 1,2 and a 1 or 2-child
        root%rank = root%rank + 1
        root%left%rank = root%left%rank - 1
        if (root%left%is_leaf()) root%left%rank = root%left%rank - 1
      else if (root%right%left_rank_diff() == 1) then ! 3,1 root with 1,2 right child
        call rotate_right(root%right)
        call rotate_left(root)  ! new root is 2,2 and a 1 or 2-child
        root%rank = root%rank + 2
        root%left%rank = root%left%rank - 2
        root%right%rank = root%right%rank - 1
      else  ! 3,1 root with a 2,2 right child
        root%right%rank = root%right%rank - 1  ! root%right is 1,1
        root%rank = root%rank - 1 ! root is 2,1 and a 2 or 3-child
      end if
    end if
  end subroutine

  !! WAVL procedure for rebalancing the subtree at ROOT after a removal
  !! in its right subtree resulted in it becoming a 2,2 leaf or with a
  !! possible right 3-child.

  pure subroutine remove_right_rebalance(root)
    type(rbt_node), pointer, intent(inout) :: root
    if (root%is_leaf()) then
      root%rank = 0 ! makes root 1,1
    else if (root%right_rank_diff() == 3) then
      if (root%left_rank_diff() == 2) then ! root is 2,3
        root%rank = root%rank - 1 ! root is 1,2 and a 2 or 3-child
      else if (root%left%left_rank_diff() == 1) then ! 1,3 root with 1,1 or 1,2 left child
        call rotate_right(root)  ! new root is 2,1 and a 1 or 2-child
        root%rank = root%rank + 1
        root%right%rank = root%right%rank - 1
        if (root%right%is_leaf()) root%right%rank = root%right%rank - 1
      else if (root%left%right_rank_diff() == 1) then ! 1,3 root with 2,1 left child
        call rotate_left(root%left)
        call rotate_right(root)  ! new root is 2,2 and a 1 or 2-child
        root%rank = root%rank + 2
        root%right%rank = root%right%rank - 2
        root%left%rank = root%left%rank - 1
      else
        root%left%rank = root%left%rank - 1  ! root%left is 1,1
        root%rank = root%rank - 1 ! root is 1,2 and a 2 or 3-child
      end if
    end if

  end subroutine

  !! Return the rank difference of the left child.
  !! The rank of a null child is -1 by convention.
  pure integer function left_rank_diff(this)
    class(rbt_node), intent(in) :: this
    if (associated(this%left)) then
      left_rank_diff = this%rank - this%left%rank
    else
      left_rank_diff = this%rank + 1
    end if
  end function

  !! Return the rank difference of the right child
  !! The rank of a null child is -1 by convention.
  pure integer function right_rank_diff(this)
    class(rbt_node), intent(in) :: this
    if (associated(this%right)) then
      right_rank_diff = this%rank - this%right%rank
    else
      right_rank_diff = this%rank + 1
    end if
  end function

  !! Return true if the node THIS is a leaf
  pure logical function is_leaf(this)
    class(rbt_node), intent(in) :: this
    is_leaf = .not.(associated(this%left) .or. associated(this%right))
  end function

  !! Right-rotate subtree at ROOT; left child becomes root.
  pure subroutine rotate_right(root)
    type(rbt_node), pointer, intent(inout) :: root
    type(rbt_node), pointer :: pivot
    pivot => root%left
    root%left => pivot%right
    pivot%right => root
    root => pivot
  end subroutine

  !! Left-rotate subtree at ROOT; right child becomes root.
  pure subroutine rotate_left(root)
    type(rbt_node), pointer, intent(inout) :: root
    type(rbt_node), pointer :: pivot
    pivot => root%right
    root%right => pivot%left
    pivot%left => root
    root => pivot
  end subroutine

!!!! PROCEDURES USEFUL FOR TESTING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Serialize the ranked binary tree data. The data is written sequentially
  !! into ARRAY using pre-order traversal of the binary tree, with the node
  !! value followed by its rank. The size of ARRAY must be at least twice
  !! the size of the tree.

  subroutine serialize(this, array)
    class(integer_map), intent(in) :: this
    integer, intent(out) :: array(:)
    integer :: pos
    !ASSERT(size(array) >= 2*this%size())
    pos = 1
    call serialize_rbt(this%root)
  contains
    recursive subroutine serialize_rbt(root)
      type(rbt_node), pointer, intent(in) :: root
      if (associated(root)) then
        array(pos)   = root%val
        array(pos+1) = root%rank
        pos = pos + 2
        call serialize_rbt(root%left)
        call serialize_rbt(root%right)
      end if
    end subroutine
  end subroutine

  !! Deserialize data that describes a ranked binary tree. ARRAY contains
  !! the data according to a pre-order traversal of the binary tree, with
  !! node value followed by its rank, such as generated by the corresponding
  !! serialize procedure.

  subroutine deserialize(this, array)
    class(integer_map), intent(out) :: this
    integer, intent(in) :: array(:)
    integer :: pos
    pos = 1
    call deserialize_rbt(this%root, huge(array))
  contains
    recursive subroutine deserialize_rbt(root, hi)
      type(rbt_node), pointer, intent(out) :: root
      integer, intent(in), optional :: hi
      root => null()
      if (pos > ubound(array,1)) return
      if (array(pos) > hi) return
      allocate(root)
      root%val  = array(pos)
      root%rank = array(pos+1)
      pos = pos + 2
      call deserialize_rbt(root%left, root%val)
      call deserialize_rbt(root%right, hi)
    end subroutine
  end subroutine

!!!! DEBUGGING PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Print the internal ranked binary tree to stdout
  subroutine print_rbt(this)
    class(integer_map), intent(in) :: this
    print *, '------------------'
    call print_rbt_(this%root, 0)
    print *, '------------------'
  contains
    recursive subroutine print_rbt_(root, level)
      type(rbt_node), pointer, intent(in) :: root
      integer, intent(in) :: level
      if (associated(root)) then
        call print_rbt_(root%left, level+1)
        print '(a,i0,"[",i0,":",i0,",",i0,"]")', repeat('  ', level), &
            root%val, root%rank, root%left_rank_diff(), root%right_rank_diff()
        call print_rbt_(root%right, level+1)
      end if
    end subroutine
  end subroutine

  !! Returns true if all the ranks satisfy the WAVL rank rules
  logical function check_ranks(this)
    class(integer_map), intent(in) :: this
    check_ranks = check_ranks_(this%root)
  contains
    recursive logical function check_ranks_(root) result(okay)
      type(rbt_node), pointer, intent(in) :: root
      integer :: d(2)
      if (associated(root)) then
        d = [root%left_rank_diff(), root%right_rank_diff()]
        if (minval(d) >= 1 .and. maxval(d) <= 2) then
          okay = check_ranks_(root%left) .and. check_ranks_(root%right)
        else
          okay = .false.
        end if
        if (.not.associated(root%left) .and. .not.associated(root%right)) &
            okay = okay .and. root%rank == 0
      else
        okay = .true.
      end if
    end function
  end function

end module integer_map_type
