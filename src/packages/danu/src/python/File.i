/* -*- C -*- */

// File Object definition and extensions
%include "DanuFileObject.h";

%extend File {
        // Constructor
        File( const char *filename, char mode='\0') {
             File *fh = construct_file_object(filename,mode);
             if (fh == NULL) {
               throw_exception("Failed to create an File object");
             }
             return fh;
        }
        // Destructor
        ~File() {
           deconstruct_file_object($self);
        }
        // Methods attached to File objects
        void close() {
           hid_t fid = GET_H5OBJECT_ID(self->h5);
           if ( H5_ISA_VALID_ID(fid) ) {
             danu_file_close(fid);
           }
        }
        void flush(int flag=0) {
          hid_t fid = GET_H5OBJECT_ID(self->h5);
          if ( flag ) {
            danu_file_flush_global(fid);
          }
          else {
            danu_file_flush(fid);
          }
        }
        PyObject * get_attribute(const char *name) {
          return readPythonAttribute($self->h5,name);
        }
        PyObject * set_attribute(const char *name,PyObject *value) {
          herr_t err = writePythonAttribute($self->h5,name,value);
          if ( H5_RETURN_OK(err) ) {
            return Py_BuildValue("");
          }
          else {
            return NULL;
          }
        }
        PyObject * attributes() {
          return readPythonAllAttributes($self->h5);
        }
        void set_offset(const char *dataset_name,int offset) {
          hid_t fid = $self->h5->hid;
          hid_t id = danu_dataset_open(fid,dataset_name);
          herr_t err = danu_set_offset(id,offset);
          if ( DANU_RETURN_FAIL(err) ) {
            throw_exception("Failed to set offset");
          }
          danu_dataset_close(id);
        }
        int get_offset(const char *dataset_name) {
          int offset;
          hid_t fid = $self->h5->hid;
          hid_t id = danu_dataset_open(fid,dataset_name);
          herr_t err = danu_get_offset(id,&offset);
          danu_dataset_close(id);
          if ( DANU_RETURN_FAIL(err) ) {
            throw_exception("Failed to read offset");
          }
          return offset;
        }

};
