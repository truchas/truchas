
    integer :: j, tmp

    ASSERT(size(array) == size(perm))

    do j = 1, size(perm)
      perm(j) = j
    end do

    !! Create the heap from the bottom up.
    do j = size(perm)/2, 1, -1 ! start with the last parent
      call sift_down (j, size(perm))
    end do

    do j = size(perm), 2, -1
      tmp = perm(j)
      perm(j) = perm(1)
      perm(1) = tmp
      if (j > 2) call sift_down (1, j-1)
    end do

  contains

    subroutine sift_down (first, last)
      integer, intent(in) :: first, last
      integer :: root, child, proot
      root = first
      proot = perm(root)
      child = 2*root
      do while (child <= last)
        if (child < last) then  ! there is a sibling, take the larger
          if (array(perm(child)) < array(perm(child+1))) child = child + 1
        end if
        if (array(proot) >= array(perm(child))) exit  ! tree at root is a heap
        perm(root) = perm(child)
        root = child
        child = 2*root
      end do
      perm(root) = proot
    end subroutine

