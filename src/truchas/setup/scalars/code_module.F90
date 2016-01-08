!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
  character(string_len), Public, Save :: code_name
  character(string_len), Public, Save :: code_version
  character(string_len), Public, Save :: libraries
  character(string_len), Public, Save :: build_date
  character(string_len), Public, Save :: build_host
  character(string_len), Public, Save :: run_architecture
  character(string_len), Public, Save :: run_host

  ! the next two sometimes get quite long and can overflow string_len
  character(1024), Public, Save :: build_architecture
  character(1024), Public, Save :: build_flags

END MODULE CODE_MODULE
