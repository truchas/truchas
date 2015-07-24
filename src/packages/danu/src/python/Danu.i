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
int danu_hex_save(const char *name, int* elements, int n1, int n2, const char
    *group_name, const char *dataset_name);
%clear (int* elements, int n1, int n2);
