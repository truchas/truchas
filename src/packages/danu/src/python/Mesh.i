/* -*- c -*- */

/* Danu files to wrap */
%include "../danu_mesh_types.h";

%include "DanuMeshObject.h";

%extend Mesh {

        /* NumPy typemaps */
        %apply(double* IN_ARRAY1, int DIM1){(double* coord_buf, int n)};
        %apply(int* IN_ARRAY2, int DIM1, int DIM2){(int* int_buf, int n1, int n2)};
        %apply(double* IN_ARRAY1, int DIM1){(double* x_in, int nx),
                                            (double* y_in, int ny),
                                            (double* z_in, int nz)};

        %apply(double* ARGOUT_ARRAY1, int DIM1){(double* out_buf, int n)};
                                            
        // Constructor
        Mesh( Output *fh, const char *meshname, tmesh_t mesh_type=INVALID_MESH, telem_t elem_type=INVALID_ELEM) {
             Mesh *m = construct_mesh_object(fh,meshname,mesh_type,elem_type);
             if (m == NULL) {
               throw_exception("Failed to create an Output object");
             }
             return m;
        }
        // Destructor
        ~Mesh() {
           deconstruct_mesh_object($self);
        }
        PyObject * get_attribute(const char *name) {
          H5Obj *h5 = GET_MOBJECT_H5OBJECT($self);
          return readPythonAttribute(h5,name);
        }
        PyObject * set_attribute(const char *name,PyObject *value) {
          H5Obj *h5 = GET_MOBJECT_H5OBJECT($self);
          herr_t err = writePythonAttribute(h5,name,value);
          if ( H5_RETURN_OK(err) ) {
            return Py_BuildValue("");
          }
          else {
            return NULL;
          }
        }
        PyObject * attributes() {
          H5Obj *h5 = GET_MOBJECT_H5OBJECT($self);
          return readPythonAllAttributes(h5);
        }
        hid_t hid() {
          return GET_MOBJECT_HID($self);
        }
        void write_coordinates(double* x_in, int nx, double* y_in=NULL, int ny, double* z_in=NULL, int nz) {

          hid_t mid = GET_MOBJECT_HID($self);
          herr_t stat;
          int in_ndim;
          
          in_ndim=1;
          if ( y_in != NULL )
            in_ndim++;

          if (z_in != NULL )
            in_ndim++;

          if ( in_ndim != $self->dim ) {
            throw_exception("Incorrect number of coordinate arrays");
            return;
          }

          if ( ( $self->dim > 1 ) && ( nx != ny ) ) {
            throw_exception("Mismatched array sizes (x.size() != y.size())");
            return;
          }

          if ( ($self->dim == 3 ) && ( nx != nz ) ) {
            throw_exception("Mismatched array sizes (x.size() != z.size())");
            return;
          }

          stat = mesh_write_coordinates(mid,nx,x_in,y_in,z_in);
          if ( DANU_RETURN_FAIL(stat) ) {
            throw_exception("Failed to write mesh coordinates");
          }

        }
        void read_coordinates(int idx, double* out_buf, int n) {
          hid_t mid = GET_MOBJECT_HID($self);
          herr_t stat;
          int nnodes;
          stat = mesh_get_nnodes(mid,&nnodes);
          if ( DANU_RETURN_FAIL(stat) ) {
            throw_exception("Failed to define number of nodes");
            return;
          }
          if ( nnodes > n ) {
            throw_exception("NumPy array request too small");
            return;
          }
          if ( nnodes < n ) {
            printf("Requested NumPy array for coordinates is too large.");
          }
          stat = mesh_read_coordinates_byindex(mid,idx,out_buf);
          if ( DANU_RETURN_FAIL(stat) ) {
            throw_exception("Failed to read mesh coordinates");
          }
        }
        PyObject * coordinates() {
          hid_t mid = GET_MOBJECT_HID($self);
          int i, ndim, nnodes;
          npy_intp dims;
          herr_t stat;
          PyObject *list=NULL;
          PyObject *coord[3];
          double *read_coord;

          stat = mesh_get_dimension(mid,&ndim);
          if ( DANU_RETURN_FAIL(stat) ) {
            throw_exception("Failed to stat mesh dimensions");
            return list;
          }
          stat = mesh_get_nnodes(mid,&nnodes);
          if ( DANU_RETURN_FAIL(stat) ) {
            throw_exception("Failed to stat mesh number of nodes");
            return list;
          }

          /* Initialize the coord pointers */
          for(i=0;i<3;i++) 
            coord[i]=NULL;

          /* Loop through and get arrays */
          dims=(npy_intp)nnodes;
          for(i=0;i<ndim;i++) {
            coord[i]=PyArray_SimpleNew(1,&dims,NPY_DOUBLE);
            if(!coord[i]) {
              Py_DECREF(coord[i]);
              throw_exception("Failed to create coordinate NumPy array");
              return list;
           }
           read_coord=(double*)array_data(coord[i]);
           stat=mesh_read_coordinates_byindex(mid,i,read_coord);
           if ( DANU_RETURN_FAIL(stat) ) {
             throw_exception("Failed to read coordinate dataset");
             return list;
           }
          }
          list=buildPyList(coord,ndim);
          return list;
        }
        void write_connectivity(int *int_buf, int n1, int n2) {
          hid_t mid = GET_MOBJECT_HID($self);
          herr_t stat;
          stat=mesh_write_connectivity(mid,n1,int_buf);
          if ( DANU_RETURN_FAIL(stat) ) {
            throw_exception("Failed to write connectivity data");
          }
        }
        PyObject * read_connectivity() {
          hid_t mid = GET_MOBJECT_HID($self);
          herr_t stat;
          int *out_data;
          int nelem, elem_order;
          PyObject *array = NULL;
          npy_intp dims[2];
          

          stat=mesh_get_nelem(mid,&nelem);
          if ( DANU_RETURN_FAIL(stat) ) {
            throw_exception("Failed to find connectivity number elements");
            return array;
          }
          stat=mesh_get_elem_order(mid,&elem_order);
          if ( DANU_RETURN_FAIL(stat) ) {
            throw_exception("Failed to find connectivity element order");
            return array;
          }
          dims[0]=(npy_intp)nelem;
          dims[1]=(npy_intp)elem_order;
          array=PyArray_SimpleNew(2,dims,NPY_INT);
          if(!array) {
            throw_exception("Failed to create connectivity array");
            return NULL;
          }
          out_data=(int*)array_data(array);
          stat=mesh_read_connectivity(mid,out_data);
          if (DANU_RETURN_FAIL(stat)) {
            Py_DECREF(array);
            throw_exception("Failed to read connectivity data");
            return NULL;
          }
          return array;
        }
        int connectivity_offset() {
          hid_t mid = GET_MOBJECT_HID($self);
          herr_t stat;
          int offset;
          stat=mesh_connectivity_get_offset(mid,&offset);
          if ( DANU_RETURN_FAIL(stat) ) {
            throw_exception("Failed to read connectivity offset");
          }
          return offset;
        }
        int nnodes(){
          hid_t mid = GET_MOBJECT_HID($self);
          herr_t stat;
          int n = -1;
          stat=mesh_get_nnodes(mid,&n);
          if ( DANU_RETURN_FAIL(stat) ) {
            throw_exception("Failed to read number of nodes");
          }
          return n;
        }
        int nelem() {
          hid_t mid = GET_MOBJECT_HID($self);
          herr_t stat;
          int n = -1;
          stat=mesh_get_nelem(mid,&n);
          if ( DANU_RETURN_FAIL(stat) ) {
            throw_exception("Failed to read number of elements");
          }
          return n;
        }
        int elem_order() {
          hid_t mid = GET_MOBJECT_HID($self);
          herr_t stat;
          int n = -1;
          stat=mesh_get_elem_order(mid,&n);
          if ( DANU_RETURN_FAIL(stat) ) {
            throw_exception("Failed to element order");
          }
          return n;
        }
        /* Clear NumPy typemaps */
        %clear (double* coord_buf, int n);
        %clear (int* int_buf, int n1, int n2);
        %clear (double* x_in, int nx),
               (double* y_in, int ny),
               (double* z_in, int nz);

        %clear (double* out_buf, int n);

};

/* Mesh related methods for Output objects */
%extend Output {

        Mesh * get_mesh(const char *meshname) {
          Mesh * mesh=NULL;
          hid_t hid = GET_OUTPUT_OBJ_HID($self);
          hid_t mid = mesh_open(hid,meshname);
          tmesh_t mesh_type;
          telem_t elem_type;
          if ( H5_ISA_VALID_ID(mid) ) {
            mesh_get_type(mid,&mesh_type);
            mesh_get_elementtype(mid,&elem_type);
            mesh = construct_mesh_object($self,meshname,mesh_type,elem_type);
          }
          else {
            throw_exception("Failed to open mesh");
          }
          return mesh;
        }
        Mesh * add_unstruct_mesh(const char *meshname, telem_t elem_type=INVALID_ELEM ) {
          tmesh_t mesh_type = UNSTRUCTURED_MESH;
          Mesh *m = construct_mesh_object($self,meshname,mesh_type,elem_type);
          return m;
        }
        int mesh_count() {
          hid_t hid = GET_OUTPUT_OBJ_HID($self);
          int num;
          herr_t err = mesh_count(hid,&num);
          if ( H5_RETURN_FAIL(err) ) {
            throw_exception("Failed to retrieve the mesh count");
          }
          return num;
        }
        int mesh_exists(const char *mesh_name) {
          hid_t hid = GET_OUTPUT_OBJ_HID($self);
          int flag = -1;
          herr_t err  = mesh_exists(hid,mesh_name,&flag);
          if (H5_RETURN_FAIL(err) ) {
            throw_exception("Failed to status existence of mesh");
          }
          return flag;
        }
        PyObject * mesh_list() {
          hid_t hid = GET_OUTPUT_OBJ_HID($self);
          int num = 0;
          int num_found;
          char **names;
          PyObject *list;
          herr_t err = mesh_count(hid,&num);
          if ( H5_RETURN_OK(err) && num ) {
            names = DANU_MALLOC(char *, num);
            err = mesh_list(hid,num,names,&num_found);
            list = convertCharListToPyList((const char * const *)names,num);
          }
          else {
            if ( H5_RETURN_FAIL(err) ) {
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
