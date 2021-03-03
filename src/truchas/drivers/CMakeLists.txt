target_sources(truchas PRIVATE
  runinfo.c
  drivers.F90
  physics_module.F90
  physical_constants.F90
  time_step_module.F90
  time_step_sync_type.F90
  simulation_event_queue.F90
  sim_event_queue_type.F90
  hijack_truchas.F90
)
