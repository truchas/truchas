target_sources(truchas PRIVATE
  data_mapper_class.F90
  kuprat_mapper_type.F90
  kuprat/gm_mesh_type.F90
  kuprat/grid_mapping_module.F90
  kuprat/grid_mapping_utils.F90
  kuprat/hpsort.F90
  kuprat/overlap_module.F90
)

if(USE_PORTAGE)
  target_sources(truchas PRIVATE portage_mapper_type.F90 portage/truchas_portage.cc)
endif()
