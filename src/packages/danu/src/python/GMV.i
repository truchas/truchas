// -*- c -*-
//
// Danu Python Interface
//
//

%module GMV
%{

//System Includes
#include <string.h>

// Python and NumPy
#include <Python.h>
#include <numpy/arrayobject.h>

#include "except.h"

#include "GMVIface.h"
#include "gmvwrite.h"

#define GMV_DFLT_TIME  -99999.0
#define GMV_DFLT_CYCLE -1

extern MeshFile * allocate_MeshFile(const char *, const char *);
extern void       deallocate_MeshFile(MeshFile *);

extern DataFile * allocate_DataFile(const char *, const MeshFile *);
extern void       deallocate_DataFile(DataFile *);

extern int is_gmvobject_valid(void *);

extern void open_gmvwrite_file(void *, const char *, int, int, int);
extern void close_gmvwrite_file(void *);

%}

/* Handling char ** correctly striaght from the interwebs
* This tells SWIG to treat char ** as a special case
*/
%typemap(in) char ** {
  /* Check if is a list */
  if (PyList_Check($input)) {
    int size = PyList_Size($input);
    int i = 0;
    $1 = (char **) malloc((size+1)*sizeof(char *));
    for (i = 0; i < size; i++) {
      PyObject *o = PyList_GetItem($input,i);
      if (PyString_Check(o))
        $1[i] = PyString_AsString(PyList_GetItem($input,i));
      else {
        PyErr_SetString(PyExc_TypeError,"list must contain strings");
        free($1);
        return NULL;
      }
    }
    $1[i] = 0;
  } else {
    PyErr_SetString(PyExc_TypeError,"not a list");
    return NULL;
  }
}

/* This cleans up the char ** array we malloc'd before the function call */
%typemap(freearg) char ** {
  free((char *) $1);
}


/* Basic error handling and HID types */
%include exception_handling.i

/* NumPY SWIG interface */
%include "numpy.i"
%init %{
import_array();
%}

/* include the numpy.i fragments */
%fragment("NumPy_Fragments");


/* Include GMV File object */
%include "GMVIface.h"

/* Extend GMVFile object */
%extend MeshFile {

  /* NumPy type maps */
  %apply(double* IN_ARRAY1, int DIM1){(double* x_in, int nx),
                                      (double* y_in, int ny),
                                      (double* z_in, int nz)};

  %apply(int* IN_ARRAY2, int DIM1, int DIM2){(int* condata, int ncells, int elem_order)};


  //Constructor
  MeshFile(const char *filename, const char * mesh_type=NULL, int isize=4, int rsize=8) {
    MeshFile *gfile=NULL;

    if ( is_gmvwrite_active() ) {
      throw_exception("Can not have more the one GMV file object active");
    }
    else {
      gfile = allocate_MeshFile(filename,mesh_type);
      gfile->nnodes=-1;
      gfile->ncells=-1;
      gfile->isize=isize;
      gfile->rsize=rsize;
    }

    return gfile;

  }
  ~MeshFile() {
    if ( is_gmvwrite_active() ) {
      close_gmvwrite_file($self);
    }
    deallocate_MeshFile($self);
  }
  int is_valid() {
    return is_gmvobject_valid($self);
  }
  int open(int binary_flag=0) {
    int status;
    if ( is_gmvwrite_active() ) {
      throw_exception("Can not have more than one GMV file open.");
      status = -1;
    }
    else {
      open_gmvwrite_file($self,$self->name,binary_flag,$self->isize,$self->rsize);
      status =  0;
    }
    return status;
  }
  int close() {
    int status;
    if ( is_gmvwrite_active() ) {
      close_gmvwrite_file($self);
      status = 0;
    }
    else {
      throw_exception("GMV write file is not active. Can not close.");
      status = -1;
    }
    return status;
  }
  void write(double* x_in, int nx, double* y_in, int ny, double* z_in, int nz, int* condata, int ncells, int elem_order) {
    int i,nc,list[8];
    int nv;
    int *ptr;

    if ( is_gmvobject_valid($self) ) {
      if ( (nx != ny) || (ny != nz) || (nx != nz) ) {
        throw_exception("Mismatched node vector lengths");
      }
      else {
        $self->nnodes=nx;
        $self->ncells=ncells;
      
        gmvwrite_node_data(&$self->nnodes,x_in,y_in,z_in);
        gmvwrite_cell_header(&$self->ncells);

        nv = mesh_type_elem_order($self->type);

        if ( nv != elem_order ) {
          throw_exception("Invalid connectivity array");
        }
        else {
          ptr=condata;
          for(nc=0; nc<ncells; nc++ ) {
            for(i=0;i<elem_order;i++) {
              list[i] = *ptr;
              ptr++;
            }
            if (list[0] == list[1]) { /* tet element */
              gmvwrite_cell_type("ptet4", 4, &(list[1]));
            } else if (list[4] == list[5]) { /* pyramid element */
              gmvwrite_cell_type("ppyrmd5", 5, list);
            } else if (list[5] == list[6]) { /* wedge element */
              i = list[1]; list[1] = list[3]; list[3] = i;  /* swap 1 and 3 */
              i = list[2]; list[2] = list[4]; list[4] = i;  /* swap 2 and 4 */
              gmvwrite_cell_type("pprism6", 6, list);
            } else { /* hex element */
              gmvwrite_cell_type("phex8", 8, list);
            }
          }
        }
  
      }

    }

  }

  /* Clear out the type maps for MeshFile extension */
  %clear (double* x_in, int nx),
         (double* y_in, int ny),
         (double* z_in, int nz);

  %clear (int* condata, int ncells, int elem_order);
};

%extend DataFile {

  %apply(double* IN_ARRAY1, int DIM1){(double* cell_data, int in_ncells),
                                      (double* node_data, int in_nnodes)}
  %apply(int* IN_ARRAY1, int DIM1){(int* cell_ids, int in_ncells),
                                   (int* node_ids, int in_nnodes)}
  %apply(double* IN_ARRAY2, int DIM1, int DIM2){(double* cell_vdata, int in_ndims, int in_ncells)}
  %apply(double* IN_ARRAY2, int DIM1, int DIM2){(double* node_vdata, int in_ndims, int in_nnodes)} 

  DataFile(const char * filename, 
           MeshFile *mesh,
           int cycle=GMV_DFLT_CYCLE,
           double t=GMV_DFLT_TIME,
           int isize=4, int rsize=8) {
  
    DataFile *file = NULL;

    if ( is_gmvwrite_active() ) {
      throw_exception("Can not have more than one GMV file object");
    }
    else {
      file = allocate_DataFile(filename,mesh);
      file->isize=isize;
      file->rsize=rsize;
      file->cycle=cycle;
      file->time=t;
      file->nnodes=mesh->nnodes;
      file->ncells=mesh->ncells;
    }

    return file;
  }
  ~DataFile() {
    if ( is_gmvwrite_active() ) {
      close_gmvwrite_file($self);
    }
    deallocate_DataFile($self);
  }
  int is_valid() {
    return is_gmvobject_valid($self);
  }
  int open(int binary_flag=0) {
    int status;

    if ( is_gmvwrite_active() ) {
      throw_exception("Can not have more than one GMV file open.");
      status = -1;
    }
    else {
      open_gmvwrite_file($self,$self->name,binary_flag,$self->isize,$self->rsize);
      gmvwrite_nodes_fromfile($self->mesh->name,$self->mesh->nnodes);
      gmvwrite_cells_fromfile($self->mesh->name,$self->mesh->ncells);
      if ( $self->cycle != GMV_DFLT_CYCLE ) {
        gmvwrite_cycleno($self->cycle);
      }
      if ( $self->time != GMV_DFLT_TIME ) {
        gmvwrite_probtime($self->time);
      }
      status =  0;
    }
    return status;
  }
  int close() {
    int status;
    if ( is_gmvwrite_active() ) {
      close_gmvwrite_file($self);
      status = 0;
    }
    else {
      throw_exception("GMV write file is not active. Can not close.");
      status = -1;
    }
    return status;
  }
  int begin_variable_block() {
     
    int ret = 0;
 
    if ( is_gmvobject_valid($self) ) {
      gmvwrite_variable_header();
    }
    else {
      throw_exception("This object is not allowed to write.");
    }

    return ret;

  }
  int end_variable_block() {
 
    int ret = 0;
 
    if ( is_gmvobject_valid($self) ) {
      gmvwrite_variable_endvars();
    }
    else {
      throw_exception("This object is not allowed to write.");
    }

    return ret;

  }
  int write_cell_data(const char *varname, double* cell_data, int in_ncells) {
   
    int ret = 0;
 
    if ( is_gmvobject_valid($self) ) {

      if ( in_ncells != $self->ncells ) {
        throw_exception("Mismatched node data length");
      }
      else {
        gmvwrite_variable_name_data(GMV_DATA_CELL,varname,cell_data);
      }
    }
    else {
      throw_exception("This object is not allowed to write.");
    }

    return ret;
  }
  int write_cell_data_vectors(const char *varname, char **comp_names, double* cell_vdata, int in_ndims, int in_ncells) {
    int ret = 0;
    int i;
    double *ptr;
    if (is_gmvobject_valid($self) ) {

      if ( in_ncells != $self->ncells ) {
        throw_exception("Mismatched number of cells. Expecting a ndim x ncell array.");
      }
      else {
        gmvwrite_vector_name(varname,GMV_DATA_CELL,in_ndims,1);
        i=0;
        while(comp_names[i]) {
          gmvwrite_vector_compname(comp_names[i]);
          i++;
        }
        gmvwrite_vector_data(GMV_DATA_CELL,in_ndims,cell_vdata);
      }
    }
    else {
      throw_exception("This object is not allowed to write");
    }

    return 0;
  }
  int write_node_data(const char * varname, double* node_data, int in_nnodes) {
    int ret = 0; 

    if ( is_gmvobject_valid($self) ) {
      if ( in_nnodes != $self->nnodes ) {
        throw_exception("Mismatched node data length");
      }
      else {
        gmvwrite_variable_name_data(GMV_DATA_NODE,varname,node_data);
      }
    }
    else {
      throw_exception("This object is not allowed to write");
    }

    return ret;

  }
  int write_node_data_vectors(const char *varname, char **comp_names, double* node_vdata, int in_ndims, int in_nnodes) {
    int ret = 0;
    int i;

    if (is_gmvobject_valid($self) ) {

      if ( in_nnodes != $self->nnodes ) {
        throw_exception("Mismatched number of nodes. Expecting ndim x nnodes array.");
      }
      else {
        gmvwrite_vector_name(varname,GMV_DATA_NODE,in_ndims,1);
        i=0;
        while(comp_names[i]) {
          gmvwrite_vector_compname(comp_names[i]);
          i++;
        }
        gmvwrite_vector_data(GMV_DATA_NODE,in_ndims,node_vdata);
      }
    }
    else {
      throw_exception("This object is not allowed to write");
    }

    return 0;
  }
  int begin_vector_block() {
    
    int ret = 0;

    if (is_gmvobject_valid($self) ) {
      gmvwrite_vector_header();
    }
    else {
      throw_exception("This object is not allowed to write");
    }

    return ret;
  }
  int end_vector_block() {
 
    int ret = 0;

    if (is_gmvobject_valid($self) ) {
      gmvwrite_vector_endvars();
    }
    else {
      throw_exception("This object is not allowed to write");
    }

    return ret;

  }
  int write_cellids(int* cell_ids, int in_ncells) {
    int ret = 0;

    if ( is_gmvobject_valid($self) ) {

      if ( in_ncells != $self->ncells ) {
        throw_exception("Mismatched number of cells");
      }
      else {
        gmvwrite_cellids(cell_ids);
      }
    }
    else {
      throw_exception("This object is not allowed to write");
    }

    return ret;
  }
  int write_nodeids(int* node_ids, int in_nnodes) {
    int ret = 0;

    if ( is_gmvobject_valid($self) ) {
 
      if ( in_nnodes != $self->nnodes ) {
        throw_exception("Mismatched number of nodes");
      }
      else {
        gmvwrite_nodeids(node_ids);
      }
    }
    else {
      throw_exception("This object is not allowed to write");
    }

    return ret;
  }
  int write_materials(int nmats, char **mat_names, int* cell_ids, int in_ncells) {
    int ret = 0;
    int i;

    if ( is_gmvobject_valid($self) ) {

      if ( in_ncells != $self->ncells ) {
        throw_exception("Mismatched number of cells");
      }
      else {
        gmvwrite_material_header(nmats,GMV_DATA_CELL);
        i=0;
        while(mat_names[i]) {
          gmvwrite_material_name(mat_names[i]);
          i++;
        }
        gmvwrite_material_ids(cell_ids,GMV_DATA_CELL);
      }
    }
    else {
      throw_exception("This object is not allowed to write");
    }

    return ret;
  }
  int begin_flag_block() {
    int ret = 0;
    
    if (is_gmvobject_valid($self) ) {
      gmvwrite_flag_header();
    }
    else {
      throw_exception("This object is not allowed to write");
    }

    return ret;
  }
  int end_flag_block() {
 
    int ret = 0;

    if (is_gmvobject_valid($self) ) {
      gmvwrite_flag_endflag();
    }
    else {
      throw_exception("This object is not allowed to write");
    }

    return ret;

  }
  int write_cell_flag(const char* flag_name, char **id_names, int num_ids, int* cell_ids, int in_ncells) {
    int ret = 0;
    int i;

    if ( is_gmvobject_valid($self) ) {

      if ( in_ncells != $self->ncells ) {
        throw_exception("Mismatched number of cells");
      }
      else {
        gmvwrite_flag_name(flag_name,num_ids,GMV_DATA_CELL);
        i=0;
        while(id_names[i]) {
          gmvwrite_flag_subname(id_names[i]);
          i++;
        }
        gmvwrite_flag_data(GMV_DATA_CELL,cell_ids);
      }

    }
    else {
      throw_exception("This object is not allowed to write");
    }

    return ret;
  }
  int write_node_flag(const char* flag_name, char **id_names, int num_ids, int* node_ids, int in_nnodes) {
    int ret = 0;
    int i;

    if ( is_gmvobject_valid($self) ) {

      if ( in_nnodes != $self->nnodes ) {
        throw_exception("Mismatched number of nodes");
      }
      else {
        gmvwrite_flag_name(flag_name,num_ids,GMV_DATA_NODE);
        i=0;
        while(id_names[i]) {
          gmvwrite_flag_subname(id_names[i]);
          i++;
        }
        gmvwrite_flag_data(GMV_DATA_NODE,node_ids);
      }
    }
    else {
      throw_exception("This object is not allowed to write");
    }

    return ret;
  }


  %clear(double* cell_data, int in_ncells),
        (double* node_data, int in_nnodes);
  %clear(int* cell_ids, int in_ncells),
        (int* node_ids, int in_nnodes);
  %clear(double* cell_vdata, int in_ncells, int in_ndims),
        (double* node_vdata, int in_nnodes, int in_ndims);                               

};
