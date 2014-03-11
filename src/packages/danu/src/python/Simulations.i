/* -*- c -*- */
// SWIG wrapped structs
%include "DanuSimulationObjects.h";

%extend Simulation {

        // Constructor
        Simulation( Output *fh, const char *sim_name) {
             Simulation *sim=create_simulation_object(fh,sim_name,1);
             if (!sim) {
               throw_exception("Failed to create Simulation object");
             }
             return sim;
        }
        // Destructor
        ~Simulation() {
           destroy_simulation_object($self);
        }
        PyObject * get_attribute(const char *name) {
          H5Obj *h5 = (H5Obj*) GET_SIMOBJ_H5OBJECT($self);
          return readPythonAttribute(h5,name);
        }
        PyObject * set_attribute(const char *name,PyObject *value) {
          H5Obj *h5 = GET_SIMOBJ_H5OBJECT($self);
          herr_t err = writePythonAttribute(h5,name,value);
          if ( H5_RETURN_OK(err) ) {
            return Py_BuildValue("");
          }
          else {
            return NULL;
          }
        }
        PyObject * attributes() {
          H5Obj *h5 = GET_SIMOBJ_H5OBJECT($self);
          return readPythonAllAttributes(h5);
        }
        hid_t hid() {
          return GET_SIMOBJ_HID($self);
        }
        int data_exists(const char * data_name) {
          hid_t sid=GET_SIMOBJ_HID($self);
          int flag=0;
          if ( DANU_RETURN_FAIL(data_exists(sid,data_name,&flag)) ) {
            throw_exception("Failed to stat non-series data");
          }
          return flag;
        }
        int data_count() {
          hid_t sid=GET_SIMOBJ_HID($self);
          int cnt=-1;
          if ( DANU_RETURN_FAIL(data_count(sid,&cnt)) ) {
            throw_exception("Failed to stat number of non-series datasets");
          }
          return cnt;
        }
        PyObject * data_list() {
          hid_t sid=GET_SIMOBJ_HID($self);
          int num = 0;
          int num_found;
          char **names;
          PyObject *list=NULL;
          herr_t err = data_count(sid,&num);
          if ( H5_RETURN_OK(err) && num ) {
            names = DANU_MALLOC(char *, num);
            err = data_list(sid,num,names,&num_found);
            list = convertCharListToPyList((const char * const *)names,num);
          }
          else {
            if ( H5_RETURN_FAIL(err) ) {
              throw_exception("Failed to determine the number of datasets");
            }
            else {
              list = Py_BuildValue("[]");
            }
          }

          return list;
        }
        void data_write(const char *data_name, PyObject *obj) {
          herr_t stat;
          PyArrayObject *array;
          int i,ndims;
          int *int_size;
          int is_new_obj;
          int typecode;
          void *data;
          if(is_array(obj)) {
            typecode=array_type(obj);
            array=obj_to_array_contiguous_allow_conversion(obj,typecode,&is_new_obj);
            ndims=array_numdims(array);
            data=array_data(array);
            int_size=DANU_MALLOC(int,ndims);
            for(i=0;i<ndims;i++)
              int_size[i]=(int)array_size(array,i);
            stat=simulation_data_write_iface($self,data_name,ndims,int_size,typecode,data);
            if(DANU_RETURN_FAIL(stat)) {
              throw_exception("Failed to write non-series dataset");
            }
            DANU_FREE(int_size);
            if ( is_new_obj && array ) {
              Py_DECREF(array);
            }
          }
          else {
            throw_exception("Invalid input in data_write. Input is not a NumPy array.");
          }
        }
        PyObject * data_read(const char *data_name) {
          PyObject *array=NULL;
          hid_t sid=GET_SIMOBJ_HID($self);
          hid_t loc;
          herr_t stat;
          int dexists;
          int typecode,ndims,i;
          hsize_t *dimensions;
          int *dims;
          npy_intp *npy_dims;
          void *data;
          stat=data_exists(sid,data_name,&dexists);
          if ( DANU_RETURN_OK(stat) && dexists ) {
            loc=data_open_group(sid);
            ndims=danu_dataset_rank(loc,data_name);
            if(ndims>0) {
              dimensions=DANU_MALLOC(hsize_t,ndims);
              npy_dims=DANU_MALLOC(npy_intp,ndims);
              dims=DANU_MALLOC(int,ndims);
              danu_dataset_dimensions(loc,data_name,ndims,dimensions);
              typecode=numpy_dataset_typecode(loc,data_name);
              for(i=0;i<ndims;i++) {
                npy_dims[i]=(npy_intp)dimensions[i];
                dims[i]=(int)dimensions[i];
              }
              array=PyArray_SimpleNew(ndims,npy_dims,typecode);
              if(array) {
                data=array_data(array);
                if(DANU_RETURN_FAIL(simulation_data_read_iface($self,data_name,ndims,dims,typecode,data) ) ) {
                  throw_exception("Failed to read non-series dataset.");
                }
              }
              else {
                throw_exception("Failed to allocate NumPy array.");
              }
              DANU_FREE(dims);
              DANU_FREE(npy_dims);
              DANU_FREE(dimensions);
            }
            else {
              throw_exception("Failed to determine dataset number of dimensions");
            }
            danu_group_close(loc);
          }
          else {
            if ( DANU_RETURN_FAIL(stat) ) {
              throw_exception("Failed to stat non-series dataset.");
            }
            else {
              throw_exception("Non-series does not exist.");
            }
          }
          return array;
        }
        char * get_sequence_name(int id) {
          char * seqname=sequence_get_name(id);
          if ( seqname == NULL ) {
            throw_exception("Failed to define sequence name");
          }
          return seqname;
        }
        Sequence * get_sequence(const char *seqname) {
          hid_t sid = GET_SIMOBJ_HID($self);
          int exists = 0;
          hid_t gid;
          Sequence * seq=NULL;
          double time;
          int cycle;

          if ( DANU_RETURN_OK(sequence_exists(sid,seqname,&exists))) {
            if ( exists ) {

              seq=create_sequence_object($self,seqname,-1,0.0);
            }
            else { 
              throw_exception("Sequence group does not exist.");
            }
          }
          else {
            throw_exception("Failed to stat sequence group");
          }

          return seq;
        }
        Sequence * get_nextSequence(int cycle, double t) {
          Sequence * seq=NULL;
          seq=create_sequence_object($self,NULL,cycle,t);
          return seq;
        }
        int sequence_count() {
          hid_t sid=GET_SIMOBJ_HID($self);
          int cnt=-1;
          herr_t stat=sequence_count(sid,&cnt);
          if (DANU_RETURN_FAIL(stat)){
            throw_exception("Failed to return sequence count");
          }
          return cnt;
        }
        int sequence_exists(const char * seriesname) {
          hid_t sid=GET_SIMOBJ_HID($self);
          herr_t stat;
          int flag=0;

          if (!seriesname) {
            throw_exception("Invalid input. Must specify id or group name");
          }

          if ( seriesname ) {
            stat=sequence_exists(sid,seriesname,&flag);
            if ( DANU_RETURN_FAIL(stat) ) {
              throw_exception("Failed to determine existence of series group");
            }
          }

          return flag;

        }
        PyObject * sequence_list() {
          hid_t sid = GET_SIMOBJ_HID($self);
          int num = 0;
          int num_found;
          char **names;
          PyObject *list;
          herr_t err = sequence_count(sid,&num);
          if ( DANU_RETURN_OK(err) && num ) {
            names = DANU_MALLOC(char *, num);
            err = sequence_list(sid,num,names,&num_found);
            if ( DANU_RETURN_OK(err) ) {
              list = convertCharListToPyList((const char * const *)names,num);
            }
            else {
              throw_exception("Failed to read sequence group list");
              list=NULL;
            }
          }
          else {
            if ( H5_RETURN_FAIL(err) ) {
              throw_exception("Failed to determine the sequence list length");
              list=NULL;
            }
            else {
              list = Py_BuildValue("[]");
            }
          }

          return list;
        }

        Mesh * link_mesh(Output * fo, const char * mesh_name) {
          hid_t sid = GET_SIMOBJ_HID($self);
          hid_t fid = GET_OUTPUT_OBJ_HID(fo);
          hid_t mid;
          herr_t r1,r2;
          Mesh * mesh = NULL;
          tmesh_t m_type;
          telem_t e_type;

          mid = simulation_link_mesh(fid,sid,mesh_name);
          if ( H5_ISA_INVALID_ID(mid) ) {
            throw_exception("Failed to link mesh");
          }
          else {
             r1 = mesh_get_type(mid,&m_type);
             r2 = mesh_get_elementtype(mid,&e_type);
          }
            
          return mesh;
        }
        
        Mesh * open_mesh_link() {
          hid_t sid = GET_SIMOBJ_HID($self);
          hid_t mid;
          Mesh * mesh = NULL;
          herr_t r1,r2;

          mid = simulation_open_mesh_link(sid);
          if ( H5_ISA_VALID_ID(mid) ) {
            mesh = construct_mesh_object_by_id(mid);
          }
          else {
            throw_exception("Failed to opne link mesh");
          }

          return mesh;
        }
        int mesh_link_exists() {
          hid_t sid = GET_SIMOBJ_HID($self);
          int flag=0;
          
          if ( DANU_RETURN_FAIL(simulation_mesh_link_exists(sid,&flag) ) ) {
            throw_exception("Failed to stat mesh link");
          }
          return flag;
        }
};

/* Methods related to Simulations for the Output Objects */
%extend Output {

        Simulation * add_simulation(const char * sim_name) {
          Simulation *sim=create_simulation_object($self,sim_name,1);
          return sim;
        }
        Simulation * get_simulation(const char * sim_name) {
          Simulation *sim=create_simulation_object($self,sim_name,0);
          return sim;
        }
        int simulation_count() {
          hid_t hid = GET_OUTPUT_OBJ_HID($self);
          int num = -1;
          if(DANU_RETURN_FAIL(simulation_count(hid,&num))) {
            throw_exception("Failed to count number of simulations");
          }
          return num;
        }
        int simulation_exists(const char *sim_name) {
          hid_t hid = GET_OUTPUT_OBJ_HID($self);
          int flag = -1;
          if ( DANU_RETURN_FAIL(simulation_exists(hid,sim_name,&flag)) ) {
            throw_exception("Failed to stat simulation group");
          }
          return flag;
        }
        PyObject * simulation_list() {
          hid_t hid = GET_OUTPUT_OBJ_HID($self);
          int num = 0;
          int num_found;
          char **names;
          PyObject *list;
          herr_t err = simulation_count(hid,&num);
          if ( DANU_RETURN_OK(err) && num ) {
            names = DANU_MALLOC(char *, num);
            err = simulation_list(hid,num,names,&num_found);
            list = convertCharListToPyList((const char * const *)names,num);
          }
          else {
            if ( DANU_RETURN_FAIL(err) ) {
              throw_exception("Failed to determine the mesh list length");
              list=NULL;
            }
            else {
              list = Py_BuildValue("[]");
            }
          }

          return list;
        }
        
};

/* Sequence Object extension */
%extend Sequence {

        // Constructor
        Sequence(Simulation *sim, const char *seriesname ) {
          return create_sequence_object(sim,seriesname,-1,0.0); 
        }
        // Destructor
        ~Sequence() {
           destroy_sequence_object($self);
        }
        PyObject * get_attribute(const char *name) {
          H5Obj *h5 = (H5Obj*) GET_SEQOBJ_H5OBJECT($self);
          return readPythonAttribute(h5,name);
        }
        PyObject * set_attribute(const char *name,PyObject *value) {
          H5Obj *h5 = GET_SEQOBJ_H5OBJECT($self);
          herr_t err = writePythonAttribute(h5,name,value);
          if ( H5_RETURN_OK(err) ) {
            return Py_BuildValue("");
          }
          else {
            return NULL;
          }
        }
        PyObject * attributes() {
          H5Obj *h5 = GET_SEQOBJ_H5OBJECT($self);
          return readPythonAllAttributes(h5);
        }
        hid_t hid() {
          return GET_SEQOBJ_HID($self);
        }
        void data_write(const char *data_name, PyObject *obj) {
          herr_t stat;
          PyArrayObject *array;
          int i,ndims;
          int *int_size;
          int is_new_obj;
          int typecode;
          void *data;
          if(is_array(obj)) {
            typecode=array_type(obj);
            array=obj_to_array_contiguous_allow_conversion(obj,typecode,&is_new_obj);
            ndims=array_numdims(array);
            data=array_data(array);
            int_size=DANU_MALLOC(int,ndims);
            for(i=0;i<ndims;i++)
              int_size[i]=(int)array_size(array,i);
            stat=sequence_data_write_iface($self,data_name,ndims,int_size,typecode,data);
            if(DANU_RETURN_FAIL(stat)) {
              throw_exception("Failed to write series dataset");
            }
            DANU_FREE(int_size);
            if ( is_new_obj && array ) {
              Py_DECREF(array);
            }
          }
          else {
            throw_exception("Invalid input in data_write. Input is not a NumPy array.");
          }
        }
        PyObject * data_read(const char *data_name) {
          PyObject *array=NULL;
          hid_t nsid=GET_SEQOBJ_HID($self);
          herr_t stat;
          int dexists=0;
          int danu_typecode,npy_type, ndims,i;
          hsize_t *dimensions;
          int *dims;
          npy_intp *npy_dims;
          void *data;
          stat=simulation_data_exists(nsid,data_name,&dexists);
          if ( DANU_RETURN_OK(stat) && dexists ) {
            ndims=simulation_data_rank(nsid,data_name);
            if(ndims>0) {
              /* Find out about the size and type of the data */
              dimensions=DANU_MALLOC(hsize_t,ndims);
              simulation_data_dimensions(nsid,data_name,ndims,dimensions);
              simulation_data_type(nsid,data_name,&danu_typecode);
              npy_type=numpy_dataset_typecode(nsid,data_name);
              /* This is a mess, THREE different types for sizes! */
              npy_dims=DANU_MALLOC(npy_intp,ndims);
              dims=DANU_MALLOC(int,ndims);
              for(i=0;i<ndims;i++) {
                npy_dims[i]=(npy_intp)dimensions[i];
                dims[i]=(int)dimensions[i];
              }
              array=PyArray_SimpleNew(ndims,npy_dims,npy_type);
              if(array) {
                data=array_data(array);
                if(DANU_RETURN_FAIL(sequence_data_read_iface($self,data_name,ndims,dims,npy_type,data) ) ) {
                  throw_exception("Failed to read non-series dataset.");
                }
              }
              else {
                throw_exception("Failed to allocate NumPy array.");
              }
              DANU_FREE(dims);
              DANU_FREE(npy_dims);
              DANU_FREE(dimensions);
            }
            else {
              throw_exception("Failed to determine dataset number of dimensions");
            }
          }
          else {
            if ( DANU_RETURN_FAIL(stat) ) {
              throw_exception("Failed to stat non-series dataset.");
            }
            else {
              throw_exception("Non-series does not exist.");
            }
          }
          return array;
        }
        int data_exists(const char *data_name) {
          hid_t nsid = GET_SEQOBJ_HID($self);
          int flag = 0;
          if ( DANU_RETURN_FAIL(simulation_data_exists(nsid,data_name,&flag) ) ) {
            throw_exception("Failed to stat data set.");
          }
          return flag;
        }
        PyObject * data_list() {
          hid_t nsid = GET_SEQOBJ_HID($self);
          int num = 0;
          int num_found;
          char **names;
          PyObject *list=NULL;
          herr_t err = simulation_data_count(nsid,&num);
          if ( DANU_RETURN_OK(err) && num ) {
            names = DANU_MALLOC(char *, num);
            err = simulation_data_list(nsid,num,names,&num_found);
            if ( DANU_RETURN_OK(err) ) {
              list = convertCharListToPyList((const char * const *)names,num);
            }
            else {
              throw_exception("Failed to read sequence data list");
              list=NULL;
            }
          }
          else {
            if ( H5_RETURN_FAIL(err) ) {
              throw_exception("Failed to determine the sequence data list length");
              list=NULL;
            }
            else {
              list = Py_BuildValue("[]");
            }
          }

          return list;
        }
        PyObject * data_attributes(const char *data_name) {
          hid_t nsid = GET_SEQOBJ_HID($self);
          hid_t id = simulation_open_data(nsid,data_name);
          H5Obj *h5 = allocate_h5_object(id,data_name);
          PyObject *ret_value=readPythonAllAttributes(h5);
          deallocate_h5_object(h5);
          danu_dataset_close(id);
          return ret_value;
        }
        PyObject * get_data_attribute(const char *data_name, const char * attr_name) {
          hid_t nsid = GET_SEQOBJ_HID($self);
          hid_t id = simulation_open_data(nsid,data_name);
          H5Obj *h5 = allocate_h5_object(id,data_name);
          PyObject *ret_value= readPythonAttribute(h5,attr_name);
          deallocate_h5_object(h5);
          danu_dataset_close(id);
          return ret_value;
        }
        int get_data_ndims(const char *data_name) {
          hid_t nsid = GET_SEQOBJ_HID($self);
          int rank = simulation_data_rank(nsid,data_name);
          return rank;
        }
        PyObject * get_data_dimensions(const char *data_name) {
            hid_t nsid = GET_SEQOBJ_HID($self);
            herr_t stat;
            int i;
            int rank = simulation_data_rank(nsid,data_name);
            hsize_t *dims;
            int * data;
            PyObject *array=NULL;
            if ( rank > 0 ) {
                dims=DANU_MALLOC(hsize_t,rank);
                data=DANU_MALLOC(int,rank);
                stat = simulation_data_dimensions(nsid,data_name,rank,dims);
                for(i=0;i<rank;i++)
                    data[i]=(int)dims[i];
                if ( DANU_RETURN_FAIL(stat) ) {
                    throw_exception("Failed to determine data set dimensions");
                    return array;
                }
                array=convertIntListToPyList(data,rank);
                DANU_FREE(dims);
                DANU_FREE(data);
            }
            else {
                throw_exception("Failed to determine data set rank");
            }
            return array;
        }
        int data_attribute_exists(const char *data_name, const char * attr_name) {
          hid_t nsid = GET_SEQOBJ_HID($self);
          hid_t id = simulation_open_data(nsid,data_name);
          int flag = 0;
          if ( DANU_RETURN_FAIL(danu_attr_exists(id,attr_name,&flag) ) ) {
            throw_exception("Failed to stat attribute");
          }
          danu_dataset_close(id);
          return flag;
        }





};

/* Probe extensions */
%extend Probe {

        // Constructor
        Probe(Simulation *sim, const char *pname, PyObject *obj ) {
          Probe *probe=NULL;
          PyArrayObject *array=NULL;
          int ndims,len, num;
          int typecode, is_new_obj;
          void * data;

          if ( is_array(obj) ) {
            typecode=array_type(obj);
            array=obj_to_array_contiguous_allow_conversion(obj,typecode,&is_new_obj);
            ndims=array_numdims(array);
            if ( ndims == 2 ) {
              data=array_data(array); 
              len=(int)array_size(array,PROBE_LEN_IDX);
              num=(int)array_size(array,PROBE_NUM_IDX);
              probe=create_probe_object(sim,pname,len,num,typecode,data);
            }
            else {
              throw_exception("Incorrect size for Probe dataset");
            }
          }
          else {
            throw_exception("Invalid input. Require a NumPy array");
          }
          if (is_new_obj && array) {
            Py_DECREF(array);
          }
          return probe; 
        }
        // Destructor
        ~Probe() {
          destroy_probe_object($self);
        }
        PyObject * get_attribute(const char *name) {
          H5Obj *h5 = (H5Obj*) GET_PRBOBJ_H5OBJECT($self);
          return readPythonAttribute(h5,name);
        }
        PyObject * set_attribute(const char *name,PyObject *value) {
          H5Obj *h5 = GET_PRBOBJ_H5OBJECT($self);
          herr_t err = writePythonAttribute(h5,name,value);
          if ( H5_RETURN_OK(err) ) {
            return Py_BuildValue("");
          }
          else {
            return NULL;
          }
        }
        PyObject * attributes() {
          H5Obj *h5 = GET_PRBOBJ_H5OBJECT($self);
          return readPythonAllAttributes(h5);
        }
        hid_t hid() {
          return GET_PRBOBJ_HID($self);
        }
        herr_t write(PyObject *obj) {
          hid_t pid = GET_PRBOBJ_HID($self);
          PyArrayObject *array=NULL;
          herr_t stat=DANU_FAILURE;
          int ndims,len,num;
          int probe_len;
          int typecode, is_new_obj;
          void * data;

          if ( is_array(obj) ) {
            typecode=array_type(obj);
            array=obj_to_array_contiguous_allow_conversion(obj,typecode,&is_new_obj);
            ndims=array_numdims(array);
            len=(int)array_size(array,PROBE_LEN_IDX);
            num=(int)array_size(array,PROBE_NUM_IDX);
            stat = probe_data_length2(pid,&probe_len);
            if ( (ndims == 2) && (len == probe_len) && (DANU_RETURN_OK(stat)) ) {
              data=array_data(array); 
              stat=probe_data_write_iface($self,num,typecode,data);
            }
            else if ( (ndims != 2 ) || (len != probe_len) ){
              throw_exception("Incorrect size for Probe dataset");
            }
            else {
              throw_exception("Failed to find probe data length");
            }
          }
          else {
            throw_exception("Invalid input. Require a NumPy array");
          }
          if (is_new_obj && array) {
            Py_DECREF(array);
          }
          return stat; 
        }
        PyObject * read() {
          hid_t pid = GET_PRBOBJ_HID($self);
          PyObject *array=NULL;
          herr_t stat=DANU_FAILURE;
          int len,num;
          npy_intp npy_dims[2];
          int typecode,danu_type;
          void * data;
  
          stat = probe_data_length2(pid,&len);
          stat &= probe_data_num2(pid,&num);
          stat  &= probe_data_type2(pid,&danu_type);

          if ( DANU_RETURN_OK(stat) ) {
            typecode=numpy_convert_typecode(danu_type);
            npy_dims[PROBE_LEN_IDX]= (npy_intp)len;
            npy_dims[PROBE_NUM_IDX]= (npy_intp)num;
            array=PyArray_SimpleNew(2,npy_dims,typecode);
            if(array) {
              data=array_data(array);
              stat=probe_data_read_iface($self,typecode,data);
              if ( DANU_RETURN_FAIL(stat) ) {
                throw_exception("Failed to read probe data");
              }
            }
            else {
              throw_exception("Failed to allocate NumPy array");
            }
          }
          else {
            throw_exception("Failed to stat probe dataset size or type");
          }

          return array; 
        }
        int length() {
          hid_t pid = GET_PRBOBJ_HID($self);
          int len=-1;
          if ( DANU_RETURN_FAIL(probe_data_length2(pid,&len) ) ) {
            throw_exception("Failed to determine probe data length");
          }
          return len;
        }
        int number() {
          hid_t pid = GET_PRBOBJ_HID($self);
          int num=-1;
          if ( DANU_RETURN_FAIL(probe_data_num2(pid,&num) ) ) {
            throw_exception("Failed to determine number of probe data vectors");
          }
          return num;
        }



};

/* Probe related methods for Simulation objects */

%extend Simulation {

  int probe_exists(const char *pname) {
    hid_t sid = GET_SIMOBJ_HID($self);
    int flag=0;
    herr_t err = probe_exists(sid,pname,&flag);
    if ( DANU_RETURN_FAIL(err) ) {
      throw_exception("Failed to stat probe dataset");
    }
    return flag;
  }
  int probe_count() {
    hid_t sid = GET_SIMOBJ_HID($self);
    int cnt=0;
    herr_t err = probe_count(sid,&cnt);
    if ( DANU_RETURN_FAIL(err) ) {
      throw_exception("Failed to count number of probe datasets");
    }
    return cnt;
  }
  PyObject * probe_list() {
    hid_t sid = GET_SIMOBJ_HID($self);
    int num = 0;
    int num_found;
    char **names;
    PyObject *list;
    herr_t err = probe_count(sid,&num);
    if ( DANU_RETURN_OK(err) && num ) {
       names = DANU_MALLOC(char *, num);
       err = probe_list(sid,num,names,&num_found);
       if ( DANU_RETURN_OK(err) ) {
          list = convertCharListToPyList((const char * const *)names,num);
       }
       else {
          throw_exception("Failed to read sequence group list");
          list=NULL;
       }
    }
    else {
       if ( DANU_RETURN_FAIL(err) ) {
         throw_exception("Failed to determine the mesh list length");
         list=NULL;
        }
        else {
          list = Py_BuildValue("[]");
        }
     }

     return list;
  }
  Probe * probe_open(const char *pname) {
    hid_t sid = GET_SIMOBJ_HID($self);
    Probe *probe=NULL;
    herr_t pid = probe_open_data(sid,pname);

    if ( H5_ISA_VALID_ID(pid) ) {
      probe=allocate_probe_object(pid,pname);
    }
    else {
      throw_exception("Failed to open probe dataset");
    }

    return probe;
  }

};
