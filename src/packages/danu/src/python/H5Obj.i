/* -*- c -*- */

/* HDF5 Object interface */
%include exception_handling.i
%include hid_t.i


%{

#include "DanuH5Object.h"

%}

// HDF5 Object definition
%include "DanuH5Obj.h";

%extend H5Obj {

        // Constructor
        H5Obj( hid_t hid, const char *name) {
            H5Obj *h5 = allocate_h5_object(hid,name);
            if ( h5 == NULL ) {
              throw_exception("Failed to allocate H5 object");
            }
            return h5;
        }
        // Destructor
        ~H5Obj() {
            DANU_FREE($self->name);
            DANU_FREE($self);
        }
        // Methods attached to H5Objects
        void write_attribute(const char *name, PyObject *value) {
          herr_t err = writePythonAttribute($self,name,value);
          if ( H5_RETURN_FAIL(err) ) {
            throw_exception("Failed to write attribute");
          }
        }
        PyObject *read_attribute(const char *name) {
          return readPythonAttribute($self,name);
        }
        PyObject *attributes() {
          return readPythonAllAttributes($self);
        }
};
