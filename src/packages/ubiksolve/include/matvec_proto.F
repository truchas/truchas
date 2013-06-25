! -*- Mode: f90 -*-

    ! <arguments>
interface
   ! perform matrix-vector multiplication, y = Ax
   subroutine MatVec (x, y, status)
     use ubik_types_module
     use ubik_vector_type_module
     type(vector_type), intent(inout) :: x
     real(real_type), intent(inout), dimension(:), target :: y
     integer(int_type), intent(out) :: status
   end subroutine MatVec
end interface
