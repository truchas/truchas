#include <hdf5.h>
#include <stdlib.h>

#define FILE        "vlstra.h5"

int verbosity = 1;
int num_errs = 0;
#define FAIL -1



static void test_write_vl_string_attribute(void)
{
    hid_t file, root, dataspace, att;
    hid_t type;
    herr_t ret;
    //char *string_att[2] = {"This is array attribute",
    //                             "It has two strings"};
    char *string_att[1] = {"This is my test string"}; 
    hsize_t   dims[1] = {1};

    file = H5Fcreate(FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /* Create a datatype to refer to. */
    type = H5Tcopy (H5T_C_S1);
    ret = H5Tset_size (type, H5T_VARIABLE);
    root = H5Gopen(file, "/",H5P_DEFAULT);
    dataspace = H5Screate_simple(1, dims, NULL);
    att = H5Acreate(root, "test_string_array", type, dataspace, H5P_DEFAULT,H5P_DEFAULT);
    ret = H5Awrite(att, type, &string_att);
    ret = H5Aclose(att);
    ret = H5Gclose(root);
    ret = H5Tclose(type);
    ret = H5Sclose(dataspace);
    ret = H5Fclose(file);
    return;
}

static void test_read_vl_string_attribute(void)
{
    hid_t file, root, att;
    hid_t type, ftype;
    herr_t ret;
    char *string_attr[1];
    H5T_class_t type_class;
    size_t size;
    htri_t size_var;

    file = H5Fopen(FILE, H5F_ACC_RDONLY, H5P_DEFAULT);
#ifdef EIP 
We can use H5Tget_native_type to disciver memory type

    /* Create a datatype to refer to. */
    type = H5Tcopy (H5T_C_S1);

    ret = H5Tset_size (type, H5T_VARIABLE);
#endif 
    root = H5Gopen(file, "/",H5P_DEFAULT);

    att = H5Aopen_name(root, "test_string_array");

    ftype = H5Aget_type(att);

    type_class = H5Tget_class (ftype);   

    if (type_class == H5T_STRING) printf ("File datatype has class H5T_STRING\n");
    size = H5Tget_size(ftype);

    printf(" Size is of the file datatype returned by H5Tget_size %d \n This is a size of char pointer\n Use H5Tis_variable_str call instead \n", size);

    if((size_var = H5Tis_variable_str(ftype)) == 1)
        printf(" to find if string has variable size \n");

    type = H5Tget_native_type(ftype, H5T_DIR_ASCEND);
    ret = H5Aread(att, type, &string_attr);

    //fprintf(stderr, "Attribute read was \n '%s'\n '%s' \n", string_attr[0], string_attr[1]);
    fprintf(stderr, "Attribute read was \n '%s'\n ", string_attr[0]);
    free(string_attr[0]);
    //free(string_attr[1]);
    ret = H5Aclose(att);
    ret = H5Tclose(type);
    ret = H5Gclose(root);
    ret = H5Fclose(file);
    return;
}

int main (int argc, char *argv[])
{
      test_write_vl_string_attribute();

      test_read_vl_string_attribute();

      return 0;
}     


