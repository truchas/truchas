module output_data_module
  ! holds TBROOK data
  use parameter_module, only: string_len, mops
  use tbrook_module,  only: brook, B_IFORM_BINARY
  use truchas_env, only: prefix, indir=>input_dir, outdir=>output_dir, input_file, title
  use output_control, only: retain_last_step
  implicit none
  public
  
  integer, save, public :: out_tbrook, &
                           err_tbrook, &
                           aux_lun, &
                           int_lun, &
                           msh_lun, &
                           tab_lun, &
                           tty_tbrook, &
                           tty_lun, &
                           out_lun, &
                           err_lun

  type(Brook), save, public, target :: odm_b_out, &
                                       odm_b_err, &
                                       odm_b_aux, &
                                       odm_b_int, &
                                       odm_b_tty

  logical, save, public :: cycle_tag_open = .false.

  character(len=32), save, public :: xml_data_format = 'binary'
  integer,           save, public :: xml_iformat = B_IFORM_BINARY ! Internal use only
  !integer, dimension(mops), save, PUBLIC :: XML_OUTPUT_DT_MULTIPLIER = 1
  !logical, save, public :: retain_last_step = .false.
  
  logical, save, public :: enable_tbrook_output = .false.

end module output_data_module
