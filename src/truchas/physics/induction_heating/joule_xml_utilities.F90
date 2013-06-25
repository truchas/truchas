!!
!! JOULE_XML_UTILITIES
!!
!! This module provides the generic procedure WRITE_VAR that writes a FILEVAR
!! tag to an xml file with the actual data going to a look-aside file.  This
!! is intended to put a sane face on the tbrook crap, and is for use by
!! EM_DATA_PROXY only!
!!
!! Neil N. Carlson <nnc@newmexico.com>
!! 26 Jan 2005
!!
!! USAGE
!!
!! CALL WRITE_VAR (BXML, BBIN, NAME, VAR, STATUS) writes a FILEVAR tag to
!!   the xml brook BXML.  The character string NAME is used for the value of
!!   the name attribute of the tag.  The data passed in VAR is written to the
!!   look-aside brook BBIN.  VAR may be a double precision scalar or rank-1
!!   array, or an integer scalar.  This may be extended to other kinds by
!!   adding the appropriate specific procedure.  If an I/O errror occurs, a
!!   non-zero value is returned in the integer STATUS.
!! 

module joule_xml_utilities

  use brook_module
  use tbrook_module
  use string_utilities, only: i_to_c

  public :: write_var

  interface write_var
    module procedure write_var_D0, write_var_D1
    module procedure write_var_I0
  end interface

contains

  subroutine write_var_D0 (bxml, bbin, name, var, status)
    type(brook), intent(inout), target :: bxml, bbin
    character(len=*), intent(in) :: name
    double precision, intent(in) :: var
    integer, intent(out) :: status
    call tbrook_writexmltag (bxml, XMLTag='FILEVAR', &
        XMLAttributes='Name="' // trim(name) // '" DataType="d" Rank="0" Shape="" Offset="' // &
            i_to_c(brook_offset(bbin)) // '"', &
        scope=TB_SCOPE_LOCAL, istatus= status)
    if (status /= 0) return
    call tbrook_write (b=bbin, variable=var, scope=TB_SCOPE_LOCAL, istatus=status)
    if (status /= 0) return
  end subroutine write_var_D0

  subroutine write_var_I0 (bxml, bbin, name, var, status)
    type(brook), intent(inout), target :: bxml, bbin
    character(len=*), intent(in) :: name
    integer, intent(in) :: var
    integer, intent(out) :: status
    call tbrook_writexmltag (bxml, XMLTag='FILEVAR', &
        XMLAttributes='Name="' // trim(name) // '" DataType="i" Rank="0" Shape="" Offset="' // &
            i_to_c(brook_offset(bbin)) // '"', &
        scope=TB_SCOPE_LOCAL, istatus= status)
    if (status /= 0) return
    call tbrook_write (b=bbin, variable=var, scope=TB_SCOPE_LOCAL, istatus=status)
    if (status /= 0) return
  end subroutine write_var_I0

  subroutine write_var_D1 (bxml, bbin, name, var, status)
    type(brook), intent(inout), target :: bxml, bbin
    character(len=*), intent(in) :: name
    double precision, intent(in) :: var(:)
    integer, intent(out) :: status
    call tbrook_writexmltag (bxml, XMLTag='FILEVAR', &
        XMLAttributes='Name="' // trim(name) // '" DataType="d" Rank="1" Shape="' // &
            i_to_c(size(var)) // '" Offset="' // i_to_c(brook_offset(bbin)) // '"', &
        scope=TB_SCOPE_LOCAL, istatus= status)
    if (status /= 0) return
    call tbrook_write (b=bbin, variable=var, scope=TB_SCOPE_LOCAL, istatus=status)
    if (status /= 0) return
  end subroutine write_var_D1

end module joule_xml_utilities
