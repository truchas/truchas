MODULE CODE_MODULE
  !=======================================================================
  ! Purpose:
  !
  !    Define variables that describe the code: it's name, version,
  !    etc, and global flag variables.
  !
  ! Author(s): Douglas B. Kothe (dbk@lanl.gov)
  !=======================================================================
  use parameter_module, only: string_len

  implicit none

  private

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  
  ! character strings defining build and run time code information
  Character (LEN=string_len), Public, Save :: code_name
  Character (LEN=string_len), Public, Save :: code_version
  Character (LEN=string_len), Public, Save :: libraries
  Character (LEN=string_len), Public, Save :: build_date
  Character (LEN=string_len), Public, Save :: build_host
  Character (LEN=string_len), Public, Save :: run_architecture
  Character (LEN=string_len), Public, Save :: run_host

  ! the next two sometimes get quite long and can overflow string_len
  Character (LEN=1024), Public, Save :: build_architecture
  Character (LEN=1024), Public, Save :: build_flags

END MODULE CODE_MODULE
