!!
!! INTEGER_SET_TYPE_WAVL
!!
!! This module defines an INTEGER_SET container that stores unique integer
!! values. This type has the same interface as the same-named type from the
!! original INTEGER_SET_TYPE module but uses a weak AVL binary tree internally.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! December 2020
!!
!! This new implementation retains the original interface but replaces the
!! internal simple linear search on an ordered linked list data structure
!! with a weak AVL (WAVL) binary search tree and bottom up rebalancing to
!! store the values. As such it is suitable for handling very large sets and
!! provides enormous speedups over the original in some Truchas use cases.
!!
!! Haeupler, B., Sen, S., and Tarjan, S. E. 2015. Rank-Balanced Trees.
!! ACM Trans. Algorithms 11, 4. DOI:https://doi.org/10.1145/2689412
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module integer_set_type_wavl

  implicit none
  private

  type, public :: integer_set
    private
    type(rbt_node), pointer :: root => null()
  contains
    procedure :: is_empty => set_is_empty
    procedure :: size => set_size
    procedure, private :: set_add, set_add_set, set_add_array
    generic   :: add => set_add, set_add_set, set_add_array
    procedure :: remove => set_remove
    procedure :: copy_to_array
    procedure, private, pass(rhs) :: set_to_array
    generic :: assignment(=) => set_to_array
    final :: integer_set_delete
    !! Primarily testing procedures
    procedure :: init => deserialize
    procedure :: serialize
    !! Debugging procedures
    procedure :: print_rbt
    procedure :: check_ranks
  end type

  type :: rbt_node
    integer :: val
    integer :: rank = 0
    type(rbt_node), pointer :: left => null(), right => null()
  contains
    procedure :: left_rank_diff
    procedure :: right_rank_diff
    procedure :: is_leaf
  end type

  !! To access the values of a set one may copy them into an array using
  !! the copy_to_array method or by assignment to an allocatable array.
  !! An alternative that doesn't require copying them into a separate
  !! array is the following iterator type. Inorder traversal of the binary
  !! tree is very simple using recursion, except that we would require
  !! returning control back to the caller after each iteration. This would
  !! be possible with coroutines, but without them the iterator needs to
  !! explicitly manage the path from the root to the current node using
  !! a stack that would otherwise be automatically handled by recursion.

  type :: node_stack
    type(node_stack_item), pointer :: head => null()
  contains
    procedure :: push
    procedure :: pop
    procedure :: top
    procedure :: is_empty
    final :: delete_node_stack
  end type

  type :: node_stack_item
    type(rbt_node), pointer :: node => null()
    type(node_stack_item), pointer :: next => null()
  end type

  type, public :: integer_set_iterator
    type(node_stack) :: path
  contains
    procedure :: begin
    procedure :: value
    procedure :: next
    procedure :: at_end
  end type

contains

!!!! NODE_STACK TYPE BOUND PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Push a pointer to NODE onto the stack
  subroutine push(this, node)
    class(node_stack), intent(inout) :: this
    type(rbt_node), pointer, intent(in) :: node
    type(node_stack_item), pointer :: rest
    rest => this%head
    allocate(this%head)
    this%head%node => node
    this%head%next => rest
  end subroutine

  !! Pop the top node pointer off the stack
  pure subroutine pop(this)
    class(node_stack), intent(inout) :: this
    type(node_stack_item), pointer :: head
    if (associated(this%head)) then
      head => this%head
      this%head => this%head%next
      deallocate(head)
    end if
  end subroutine

  !! Return a pointer to the top node on the stack
  function top(this) result(node)
    class(node_stack), intent(in) :: this
    type(rbt_node), pointer :: node
    if (associated(this%head)) then
      node => this%head%node
    else
      node => null()
    end if
  end function

  !! Return true if the stack is empty
  pure logical function is_empty(this)
    class(node_stack), intent(in) :: this
    is_empty = .not.associated(this%head)
  end function

  !! Final subroutine for NODE_STACK objects
  elemental subroutine delete_node_stack(this)
    type(node_stack), intent(inout) :: this
    do while (.not.is_empty(this))
      call pop(this)
    end do
  end subroutine

!!!! INTEGER_SET_ITERATOR PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine begin(this, set)
    class(integer_set_iterator), intent(out) :: this
    type(integer_set), intent(in) :: set
    type(rbt_node), pointer :: node
    node => set%root
    do while (associated(node))
      call this%path%push(node)
      node => node%left
    end do
  end subroutine

  logical function at_end(this)
    class(integer_set_iterator), intent(in) :: this
    at_end = this%path%is_empty()
  end function

  integer function value(this)
    class(integer_set_iterator), intent(in) :: this
    type(rbt_node), pointer :: curr
    curr => this%path%top()
    value = curr%val
  end function

  subroutine next(this)
    class(integer_set_iterator), intent(inout) :: this
    type(rbt_node), pointer :: curr, prev
    if (at_end(this)) return  ! no-op
    curr => this%path%top()
    if (associated(curr%right)) then ! advance to next left-most node
      curr => curr%right
      do while (associated(curr))
        call this%path%push(curr)
        curr => curr%left
      end do
    else
      do
        prev => curr
        call this%path%pop
        if (this%path%is_empty()) return ! done
        curr => this%path%top()
        if (associated(prev, curr%left)) return
      end do
    end if
  end subroutine

!!!! INTEGER_SET PROCEDURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Final subroutine for integer_set objects
  elemental subroutine integer_set_delete(this)
    type(integer_set), intent(inout) :: this
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

!  !! Returns the number of elements in the set.
!  integer function set_size(this)
!    class(integer_set), intent(in) :: this
!    set_size = rbt_size(this%root)
!  contains
!    recursive integer function rbt_size(root) result(n)
!      type(rbt_node), pointer, intent(in) :: root
!      if (associated(root)) then
!        n = 1 + rbt_size(root%left) + rbt_size(root%right)
!      else
!        n = 0
!      end if
!    end function
!  end function

  !! Returns the number of elements in the set.
  elemental integer function set_size(this)
    class(integer_set), intent(in) :: this
    if (associated(this%root)) then
      set_size = rbt_size(this%root)
    else
      set_size = 0
    end if
  contains
    pure recursive integer function rbt_size(root) result(n)
      type(rbt_node), intent(in) :: root
      n = 1
      if (associated(root%left))  n = n + rbt_size(root%left)
      if (associated(root%right)) n = n + rbt_size(root%right)
    end function
  end function

  !! Returns true if the set is empty
  elemental logical function set_is_empty(this)
    class(integer_set), intent(in) :: this
    set_is_empty = .not.associated(this%root)
  end function

  !! Copy the ordered set values to ARRAY. Its size must be sufficiently
  !! large to store the values. Unused array elements are left as is.

  subroutine copy_to_array(this, array)
    class(integer_set), intent(in) :: this
    integer, intent(inout) :: array(:)
    type(integer_set_iterator) :: iter
    integer :: n
    n = 0
    call iter%begin(this)
    do while (.not.iter%at_end())
      n = n + 1
      array(n) = iter%value()
      call iter%next
    end do
  end subroutine

  !! Defined assignment of the set to an allocatable rank-1 array.
  !! The array is allocated/reallocated to the correct size, and the
  !! sorted values in the set are written to the array.  If the set
  !! is empty the array is allocated with 0 size.

  subroutine set_to_array(lhs, rhs)
    integer, allocatable, intent(inout) :: lhs(:)
    class(integer_set), intent(in) :: rhs
    type(integer_set_iterator) :: iter
    integer :: n
    n = rhs%size()
    if (allocated(lhs)) then
      if (size(lhs) /= n) deallocate(lhs)
    end if
    if (.not.allocated(lhs)) allocate(lhs(n))
    call iter%begin(rhs)
    n = 0
    do while (.not.iter%at_end())
      n = n + 1
      lhs(n) = iter%value()
      call iter%next
    end do
  end subroutine

  !! Add VALUE to the set
#ifdef GNU_PR69563
  subroutine set_add (this, value)
#else
  elemental subroutine set_add (this, value)
#endif
    class(integer_set), intent(inout) :: this
    integer, intent(in) :: value
    call rbt_insert(this%root, value)
  end subroutine

  !! Add values from SET
  subroutine set_add_set(this, set)
    class(integer_set), intent(inout) :: this
    class(integer_set), intent(in) :: set
    type(integer_set_iterator) :: iter
    call iter%begin(set)
    do while (.not.iter%at_end())
      call this%add(iter%value())
      call iter%next
    end do
  end subroutine

  !! Add values from ARRAY
  subroutine set_add_array(this, array)
    class(integer_set), intent(inout) :: this
    integer, intent(in) :: array(:)
    integer :: j
    do j = 1, size(array)
      call this%add(array(j))
    end do
  end subroutine

  !! Insert VAL into the subtree at ROOT. The ROOT pointer may be modified
  !! as a result of the insertion and rebalancing of the tree.

  pure recursive subroutine rbt_insert(root, val)
    type(rbt_node), pointer, intent(inout) :: root
    integer, intent(in) :: val
    if (.not.associated(root)) then
      allocate(root)
      root%val = val
    else if (val < root%val) then
      call rbt_insert(root%left, val)  ! output root%left may be a 0-child
      if (root%left_rank_diff() == 0) call insert_left_rebalance(root)
    else if (val > root%val) then
      call rbt_insert(root%right, val) ! output root%right may be a 0-child
      if (root%right_rank_diff() == 0) call insert_right_rebalance(root)
    end if
  end subroutine

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

  !! Remove VAL from the set
  elemental subroutine set_remove(this, val)
    class(integer_set), intent(inout) :: this
    integer, intent(in) :: val
    call rbt_remove(this%root, val)
  end subroutine

  !! Remove VAL from the subtree at ROOT. The ROOT pointer may be modified
  !! as a result of the removal and rebalancing of the tree.

  pure recursive subroutine rbt_remove(root, val)
    type(rbt_node), pointer, intent(inout) :: root
    integer, intent(in) :: val
    type(rbt_node), pointer :: next, temp
    if (.not.associated(root)) return
    if (val < root%val) then
      call rbt_remove(root%left, val)
      call remove_left_rebalance(root)
    else if (val > root%val) then
      call rbt_remove(root%right, val)
      call remove_right_rebalance(root)
    else  ! val == root%val
      if (associated(root%left) .and. associated(root%right)) then
        !! Replace root value with next larger value and then
        !! delete that value from the right subtree.
        next => root%right
        do while (associated(next%left))
          next => next%left
        end do
        root%val = next%val
        call rbt_remove(root%right, next%val)
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
    class(integer_set), intent(in) :: this
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
    class(integer_set), intent(out) :: this
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
    class(integer_set), intent(in) :: this
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
    class(integer_set), intent(in) :: this
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

end module integer_set_type_wavl
