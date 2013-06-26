# TRUCHAS_ROOT is a relative path to the root of the T source tree
TRUCHAS_ROOT = ../../../..

-include ${TRUCHAS_ROOT}/configuration

.SUFFIXES :
.SUFFIXES : .c .f90 .o .so

all : asciiwriter.so ensight.so gridmap.so exodus.so

# --- shared objects -----------------------------------------------------------

asciiwriter.so : asciiwriter.o
	${CL} ${CSHARED_FLAGS} asciiwriter.o ${CCOMPILER_LIBRARIES} -o asciiwriter.so
	cp asciiwriter.so ${TRUCHAS_ROOT}/bin

ensight.so : ensight.o
	${CL} ${CSHARED_FLAGS} ensight.o ${CCOMPILER_LIBRARIES} -o ensight.so
	cp ensight.so ${TRUCHAS_ROOT}/bin

gridmap.so : gridmap_c.o gridmap_f.o hpsort.o overlap_module.o gm_mesh_type.o grid_mapping_utils.o grid_mapping_module.o
	${FL} ${FSHARED_FLAGS} gridmap_c.o gridmap_f.o hpsort.o overlap_module.o gm_mesh_type.o grid_mapping_utils.o grid_mapping_module.o ${FCOMPILER_LIBRARIES} -o gridmap.so
	cp gridmap.so ${TRUCHAS_ROOT}/bin

exodus.so : exodus_c.o exodus_f.o exodus.o exodus_mesh_type.o exodus_mesh_reader.o exodus_mesh_writer.o exodus_errors.o exodus_truchas_hack.o exodus_mesh_utilities.o string_utilities.o
	${FL} ${FSHARED_FLAGS} exodus_c.o exodus_f.o exodus.o exodus_mesh_type.o exodus_mesh_reader.o exodus_mesh_writer.o exodus_errors.o exodus_truchas_hack.o exodus_mesh_utilities.o string_utilities.o ${NETCDF_LIBRARIES} ${FCOMPILER_LIBRARIES} -o exodus.so
	cp exodus.so ${TRUCHAS_ROOT}/bin

# ------------------------------------------------------------------------------

clean :
	rm -rf *.o *.g *.g90 *.so *.mod ${TRUCHAS_ROOT}/bin/asciiwriter.so ${TRUCHAS_ROOT}/bin/ensight.so ${TRUCHAS_ROOT}/bin/gridmap.so ${TRUCHAS_ROOT}/bin/exodus.so

diagnostic:
	@echo ${BUILD_DIR}

# ------------------------------------------------------------------------------

%.o : %.c
	${CC} -c ${CFLAGS} $<

%.o : %.f90
	${FC} -c ${FFLAGS} $<

vpath %.f90 ${BUILD_DIR}

# asciiwriter dependencies
asciiwriter.o : asciiwriter.c

# ensight dependencies
ensight.o : ensight.c

# gridmap dependencies
gridmap_c.o : gridmap_c.c
gridmap_f.o : gridmap_f.f90 grid_mapping_module.o
grid_mapping_module.o : grid_mapping_module.f90 hpsort.o gm_mesh_type.o grid_mapping_utils.o overlap_module.o
grid_mapping_utils.o : grid_mapping_utils.f90 gm_mesh_type.o
gm_mesh_type.o : gm_mesh_type.f90
hpsort.o : hpsort.f90
overlap_module.o : overlap_module.f90

# exodus dependencies
exodus_c.o : exodus_c.c
exodus_f.o : exodus_f.f90 exodus.o
exodus.o : exodus.f90 exodus_mesh_type.o exodus_mesh_reader.o exodus_mesh_writer.o exodus_errors.o exodus_truchas_hack.o exodus_mesh_utilities.o
exodus_mesh_type.o : exodus_mesh_type.f90
exodus_mesh_reader.o : exodus_mesh_reader.f90 string_utilities.o exodus_mesh_type.o exodus_errors.o
exodus_mesh_writer.o : exodus_mesh_writer.f90 string_utilities.o exodus_mesh_type.o exodus_errors.o
exodus_errors.o : exodus_errors.f90
exodus_truchas_hack.o : exodus_truchas_hack.f90
exodus_mesh_utilities.o : exodus_mesh_utilities.f90 exodus_mesh_type.o
string_utilities.o : string_utilities.f90
