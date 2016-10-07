/* -*- c -*- */

// Output Object definition and extensions wrapped by SWIG
%include "DanuOutputObject.h"
%extend Output {

        // Constructor
        Output( const char *filename, char mode='\0') {
             Output *fh = construct_output_object(filename,mode);
             if (fh == NULL) {
               throw_exception("Failed to create an Output object");
             }
             return fh;
        }
        // Destructor
        ~Output() {
           deconstruct_output_object($self);
        }
        // Methods attached to Output objects
        void close() {
           hid_t fid = GET_H5OBJECT_ID($self->h5);
           if ( H5_ISA_VALID_ID(fid) ) {
             danu_file_close(fid);
           }
        }
        void flush(int flag=0) {
          hid_t fid = GET_H5OBJECT_ID($self->h5);
          if ( flag ) {
            danu_file_flush_global(fid);
          }
          else {
            danu_file_flush(fid);
          }
        }
        PyObject * get_attribute(const char *name) {
          H5Obj *h5 = $self->h5;
          return readPythonAttribute(h5,name);
        }
        PyObject * set_attribute(const char *name,PyObject *value) {
          H5Obj *h5 = $self->h5;
          herr_t err = writePythonAttribute(h5,name,value);
          if ( H5_RETURN_OK(err) ) {
            return Py_BuildValue("");
          }
          else {
            return NULL;
          }
        }
        PyObject * attributes() {
          H5Obj *h5 = $self->h5;
          return readPythonAllAttributes(h5);
        }
        hid_t hid() {
          return GET_H5OBJECT_ID($self->h5);
        }
 };
