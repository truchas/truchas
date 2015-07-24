// -*- c -*-
//
// Danu Python Interface
//
//

%module Danu
%{

//System Includes
#include <string.h>

//HDF5 Includes
#include <hdf5.h>

// Python and NumPy
#include <Python.h>
#include <numpy/arrayobject.h>

#include <danu.h>

// Danu Python Includes
#include <DanuPyUtils.h>
#include <DanuPyException.h>
#include <DanuNumPyUtils.h>
#include <DanuH5Object.h>
#include <DanuFile.h>
#include <DanuMesh.h>
#include <DanuSimulations.h>
#include <DanuOutput.h>


%}

/* Basic error handling and HID types */
%include exception_handling.i
%include hid_t.i
%include herr_t.i

/* NumPY SWIG interface */
%include "numpy.i"
%init %{
import_array();
%}

/* include the numpy.i fragments */
%fragment("NumPy_Fragments");


// HDF5 Object definition and extensions
%include H5Obj.i

// File Object definition and extensions
%include File.i

// Output File Object definition and extensions
%include Output.i

// Mesh Object definition and extensions
%include Mesh.i

//  Simulation Objects definition and extensions
%include Simulations.i

%apply (int* IN_ARRAY2, int DIM1, int DIM2){(int* elements, int n1, int n2)};
int danu_hex_save(struct DanuH5Handle *h, int* elements, int n1, int n2,
        const char *dataset_name);
%clear (int* elements, int n1, int n2);

struct DanuH5Handle * danu_h5_create_handle();
void danu_h5_free_handle(struct DanuH5Handle *h);
void danu_h5_open(struct DanuH5Handle *h, const char *name,
        const char *group_name);
%apply (int* IN_ARRAY1, int DIM1){(int* array1d, int n)};
void danu_h5_save_int_1d_array(struct DanuH5Handle *h,
        const char *dataset_name, int* array1d, int n);
%clear (int* array1d, int n);
void danu_h5_close(struct DanuH5Handle *h);
