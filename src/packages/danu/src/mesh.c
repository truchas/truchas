/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
 * mesh.c
 *
 *  DANU Mesh Datasets, Groups and Files
 *
 *
 *  Purpose:
 *           This source file defines functions that create mesh datasets, 
 *           groups and files.
 *
 *
 */
#include <stdarg.h>
#include <string.h> 

#include <hdf5.h>


#include <danu_error.h>
#include <danu_memory.h>
#include <danu_utils.h>
#include <danu_group.h>
#include <danu_dataset.h>
#include <danu_attribute.h>
#include <danu_offset.h>
#include <danu_mesh.h>

/* Private Defines */

/* Coordinate Dataset layout and size */
#define COORD_SIZE_DIM        2
#define COORD_SIZE_NNODES_IDX 0
#define COORD_SIZE_DIM_IDX    1

/* Specific coordinate hyperslab indices */
#define COORD_X_IDX 0
#define COORD_Y_IDX 1
#define COORD_Z_IDX 2

/* Connectivity dataset layout and size */
#define CONNECT_SIZE_DIM            2
#define CONNECT_SIZE_NELEM_IDX      0
#define CONNECT_SIZE_ELEMORDER_IDX  1

/* Specific connectivity hyperslab indices */
#define CONNECT_NELEM_IDX     0
#define CONNECT_ELEMORDER_IDX 1

/* Default connectivity index */
#define CONNECT_OFFSET_DFLT 0

/* Private function prototypes */ 
herr_t elem_order_to_name(int elem_order, int dim, char * elem_name);
herr_t elem_type_to_name(telem_t etype, char * type_name);
herr_t elem_string_to_type(const char *type_name, telem_t *etype);
herr_t elem_order_to_type(int dim, int elem_order, telem_t * etype);
herr_t elem_type_to_order(telem_t  etype, int *elem_order);
herr_t mesh_type_to_string(tmesh_t mtype, char * type_name);
herr_t mesh_string_to_type(const char *type_name, tmesh_t *type);

herr_t mesh_set_type(hid_t mid, tmesh_t type);
herr_t mesh_set_elementtype(hid_t mid, telem_t type);
herr_t mesh_set_dimension(hid_t mid, int dim);
herr_t mesh_set_nnodes(hid_t mid, int nnodes);
herr_t mesh_set_nelem(hid_t mid, int nelem);
herr_t mesh_set_elem_order(hid_t mid, int order);




 /*
 * Routine: herr_t elem_order_to_name(int elem_order, int dim, char *elem_name)
 * Purpose: Return the string name associated with an integer element order
 * Description: Calling routine passes in an integer that represents the number
 *              of nodes for a particular element type. That element type string
 *              name is returned in elem_name. Routine is used to set mesh
 *              attributes. The elem_name is initalized to NULL. Calling 
 *              routine should check the return value before accessing the
 *              pointer.
 *
 * Parameters:
 *           elem_order       IN              Number of nodes for element type
 *           dim              IN              Dimension of the element
 *           elem_name        OUT             The string name of the element type
 *                              
 * Returns: A negative value if an error occurs, otherwise returns zero 
 * Errors: Possible error conditions are non-positive element order or a an element
 *          type is not found. 
 */ 
herr_t elem_order_to_name(int elem_order, int dim, char * elem_name)
{
    herr_t status = DANU_FAILURE;

    /* Check input */
    if ( dim <= 0 || dim > 3 ) {
        DANU_ERROR_MESS("invalid dimension type");
        return status;
    }

    if ( elem_order <= 0 || elem_order > 8 ) {
        DANU_ERROR_MESS("Invalid element order value");
        return status;
    }

    if ( DANU_BAD_PTR(elem_name) ) {
        DANU_ERROR_MESS("Invalid string pointer");
        return status;
    }

    /* Initialize to NULL */
    elem_name[0] = '\0';
    switch(elem_order)  {
        case TRI_ELEM_ORDER:
            sprintf(elem_name,TRI_ELEM_NAME);
            break;
        case HEX_ELEM_ORDER:
            sprintf(elem_name,HEX_ELEM_NAME);
            break;
        case 4:
            if ( dim == 2 ) {
                sprintf(elem_name,QUAD_ELEM_NAME);
            }
            else if ( dim == 3 ) {
                sprintf(elem_name,TET_ELEM_NAME);
            }
            else {
                DANU_ERROR_MESS("Unknown dimension value for 4-order element");
            }
            break;
        default:
            DANU_ERROR_MESS("Mesh element type not found");
            break;
    }

    if ( strlen(elem_name) > 0 ) {
        status = 0;
    }

    return status;
}
 /*
 * Routine: herr_t elem_order_to_type(int dim, int elem_order, telem_t *etype)
 * Purpose: Return the string name associated with an integer element order
 * Description: Calling routine passes in the dimension and an integer that 
 *              represents the number of nodes for a particular element type. 
 *              It sets etype to the matching enumerated type element. The 
 *              return value indicates if an error occurred. 
 *
 * Parameters:
 *           dim              IN              Dimension of the element
 *           elem_order       IN              Number of nodes for element type
 *           etype            OUT             The element type
 *                              
 * Returns: A negative value if an error occurs, otherwise returns zero 
 * Errors: Possible error conditions are non-positive element order or an element
 *          type is not found. 
 */
herr_t elem_order_to_type(int dim, int elem_order, telem_t *etype)
{
    herr_t status = DANU_FAILURE;

    if ( DANU_BAD_PTR(etype) ) {
        DANU_ERROR_MESS("Invalid element type pointer");
        return status;
    }
    *etype = INVALID_ELEM;


    if ( DIM_IS_INVALID(dim) ) {
        DANU_ERROR_MESS("Invalid dimension value");
        return status;
    }

    if ( EORDER_IS_INVALID(elem_order) ) {
        DANU_ERROR_MESS("Invalid element order value");
        return status;
    }

    switch(dim) {
        case 1:
	    if ( elem_order == LINE_ELEM_ORDER ) {
                *etype = LINE_ELEM;
	    }
            break;
        case 2:
            if ( elem_order == TRI_ELEM_ORDER ) {
                *etype = TRI_ELEM;
            }
            else if ( elem_order == QUAD_ELEM_ORDER ) {
                *etype = QUAD_ELEM;
            }
            break;
        default:
            if ( elem_order == TET_ELEM_ORDER ) {
                *etype = TET_ELEM;
            }
            else if ( elem_order == HEX_ELEM_ORDER ) {
                *etype = HEX_ELEM;
            }
    }

    if ( *etype != INVALID_ELEM ) {
        status = 0;
    }
    else {
        DANU_ERROR_MESS("Failed to match order to type")
    }

    return status;

}    
/*
 * Routine: herr_t elem_type_to_order(telem_t etype, int *order)
 * Purpose: Return the element order associated with element type 
 * Description: Calling routine passes an element enumerated type and
 *              sets order to the appropriate value. If an error occurs
 *              the function will return  negative value. 
 *
 * Parameters:
 *           etype       IN              The element type
 s           order       OUT             Number of nodes for element type
 *                              
 * Returns: A negative value if an error occurs, otherwise returns zero 
 * Errors: Possible error conditions are invalid element type or an element
 *          order is not found. 
 */
herr_t elem_type_to_order(telem_t etype, int *order)
{
    herr_t status = DANU_FAILURE;

    if ( DANU_BAD_PTR(order) ) {
        DANU_ERROR_MESS("Invalid element order pointer");
        return status;
    }

    if ( ETYPE_IS_INVALID(etype) ) {
        DANU_ERROR_MESS("Invalid element type ");
        return status;
    }

    *order = -1;

    switch(etype) {
        case LINE_ELEM:
            *order = LINE_ELEM_ORDER;
            break;
        case TRI_ELEM:
            *order = TRI_ELEM_ORDER;
            break;
        case QUAD_ELEM:
            *order = QUAD_ELEM_ORDER;
            break;
        case TET_ELEM:
            *order = TET_ELEM_ORDER;
            break;
        default:
            *order = HEX_ELEM_ORDER;
    }

    if ( *order != INVALID_ELEM ) {
        status = 0;
    }
    else {
        DANU_ERROR_MESS("Failed to match element type to element order");
    }

    return status;
} 



 /*
 * Routine: herr_t mesh_string_to_type(const char * type_name, tmesh_t *type)
 * Purpose: Convert string name of mesh type to enumerated type temsh_t
 * Description: Calling routine passes in string that describes the mesh
 *              type. This is translated in this routine to an 
 *              enumerated type tmesh_t. A negative value is returned if an
 *              error occurs. The value of type is initialized to INVAILD_MESH
 *              before the string comparison.
 *
 * Parameters:
 *           type_name       IN             String describing mesh type 
 *           type            OUT            Enumerated tmesh_t that matches type_name
 *                              
 * Returns: A negative value if an error occurs, otherwise returns zero 
 * Errors: Possible error conditions are invalid mesh type string or invalid 
 *         pointer. 
 */ 
herr_t mesh_string_to_type(const char * type_name, tmesh_t *type)
{
    herr_t status = DANU_FAILURE;

    if ( DANU_BAD_STRING(type_name) ) {
        DANU_ERROR_MESS("Invaild mesh type name string");
        return status;
    }

    if ( DANU_BAD_PTR(type) ) {
        DANU_ERROR_MESS("Invalid mesh type pointer");
        return status;
    }
    *type = INVALID_MESH;

    if ( strcmp(UNSTRUCT_MESH_NAME,type_name) == 0 ) {
        *type = UNSTRUCTURED_MESH;
    }
    else if ( strcmp(STRUCTURED_MESH_NAME,type_name) == 0 ) {
        *type = STRUCTURED_MESH;
    }
    else {
        DANU_ERROR_MESS("Can not match mesh type name to known value");
    }

    if ( *type != INVALID_MESH) {
        status = 0;
    }

    return status;
}

 /*
 * Routine: herr_t mesh_type_to_name(tmesh_t mtype, char *type_name)
 * Purpose: Return the string name associated with a mesh type
 * Description: Calling routine passes in an mesh type value and this function
 *              returns the string value either STRUCTURED or UNSTRUCTURED
 *              in type_name. Routine is used to set mesh attributes. Calling
 *              routine should check the return value before accessing the pointer.
 *
 * Parameters:
 *           mtype       IN              Mesh type
 *           type_name   OUT             Mesh type string name (STRUCTURED or UNSTRUCTURED)
 *                              
 * Returns: A negative value if an error occurs, otherwise returns zero. 
 * Errors: Possible error conditions are invalid mesh type or invalid string pointer.
 *
 */
 herr_t mesh_type_to_name(tmesh_t mtype, char *type_name)
 {
     herr_t status = DANU_FAILURE;

     /* Check input */
     if ( DANU_BAD_PTR(type_name) ) {
         DANU_ERROR_MESS("Invalid string pointer");
         return status;
     }

     type_name[0] = '\0';
     switch(mtype) {
         case UNSTRUCTURED_MESH: 
              sprintf(type_name, UNSTRUCT_MESH_NAME);
              break;
         case STRUCTURED_MESH: 
              sprintf(type_name, STRUCTURED_MESH_NAME);
              break;
         default:     
              DANU_ERROR_MESS("Invalid mesh type value");
              break;
     }

     if (strlen(type_name) > 0 ) {
         status = 0;
     }

     return status;
 }

/*
 * Routine: herr_t elem_type_to_name(telem_t etype, char *elem_name)
 * Purpose: Return the string name associated with a elem type
 * Description: Calling routine passes in an elem type value and this function
 *              returns the string value associated with that element type.
 *              Routine is used to set mesh attributes. Calling
 *              routine should check the return value before accessing the pointer.
 *
 * Parameters:
 *           etype       IN              Element type
 *           elem_name   OUT             Element type string name 
 *                                           TRI,QUAD,TET,HEX
 *                              
 * Returns: A negative value if an error occurs, otherwise returns zero. 
 * Errors: Possible error conditions are invalid elem type or invalid string pointer.
 *
 */
 herr_t elem_type_to_name(telem_t etype, char *elem_name)
 {
     herr_t status = DANU_FAILURE;

     /* Check input */
     if ( DANU_BAD_PTR(elem_name) ) {
         DANU_ERROR_MESS("Invalid string pointer to elem_name");
         return status;
     }

     elem_name[0] = '\0';
     switch(etype) {
         case LINE_ELEM:
              sprintf(elem_name,LINE_ELEM_NAME);
              break;
         case TRI_ELEM:
              sprintf(elem_name,TRI_ELEM_NAME);
              break;
         case QUAD_ELEM:
              sprintf(elem_name,QUAD_ELEM_NAME);
              break;
         case TET_ELEM:
              sprintf(elem_name,TET_ELEM_NAME);
              break;
         case HEX_ELEM:
              sprintf(elem_name,HEX_ELEM_NAME);
              break;
         default:
              DANU_ERROR_MESS("Invalid element type...element name not found");
              break;
     }

     if ( strlen(elem_name) > 0 ) {
         status = 0;
     }

     return status;
 }
/*
 * Routine: herr_t elem_string_to_type(const char *type_name, telem_t *etype)
 * Purpose: Return the element enumerated type that matches type_name
 * Description: Given a element type string the functioin compares that string
 *              against known element types. If a match is found, the corresponding
 *              element enumerated type is set to the appropriate value. 
 *              If no match is found, then type is INVALID_ELEM. The type_name
 *              must be an exact match to the element name. See the tou_mesh.h 
 *              header file for string name definitions.
 *
 * Parameters:
 *           type_name       IN              String containing element type name
 *           etype           OUT             Element enumerated type that 
 *                                            matches string name 
 *                                           TRI,QUAD,TET,HEX
 *                              
 * Returns: A negative value if an error occurs, otherwise returns zero. 
 * Errors: Possible error conditions are invalid element name type or invalid string pointer.
 *
 */
herr_t elem_string_to_type(const char *type_name, telem_t *etype)
{
    herr_t status = DANU_FAILURE;

    if ( DANU_BAD_STRING(type_name) ) {
       DANU_ERROR_MESS("Invaild element type name string");
       return status;
    }

    if ( DANU_BAD_PTR(etype) ) {
        DANU_ERROR_MESS("Invalid element type pointer");
        return status;
    }
    *etype = INVALID_ELEM;

    if ( strcmp(LINE_ELEM_NAME,type_name) == 0 ) {
        *etype = TRI_ELEM;
    }
    else if ( strcmp(TRI_ELEM_NAME,type_name) == 0 ) {
        *etype = TRI_ELEM;
    }
    else if ( strcmp(QUAD_ELEM_NAME,type_name) == 0 ) {
        *etype = QUAD_ELEM;
    }
    else if ( strcmp(TET_ELEM_NAME,type_name) == 0 ) {
        *etype = TET_ELEM;
    }
    else if ( strcmp(HEX_ELEM_NAME,type_name) == 0 ) {
        *etype = HEX_ELEM;
    }
    else {
        DANU_ERROR_MESS("Can not match element name type name to known value");
    }

    if ( *etype != INVALID_ELEM ) {
        status = 0;
    }

    return status;
}





 /* PUBLIC FUNCTIONS */   

 /*
 * Routine: herr_t mesh_create_root_group(hid_t fid)
 * Purpose: Create the root mesh group under the file fid
 * Description: Create the root mesh group with name MESH_ROOT_GROUP_NAME
 *              under the root group in file fid. Group will be closed before 
 *              returning. Errors are raised if fid is an invalid HDF5 identifier
 *              or the routine fails to create the group.
 *
 * Parameters:
 *           fid                IN              HDF5 identifier for the root file
 *                              
 * Returns: A negative value if an error occurs, otherwise returns zero 
 *          an error occurs.
 * Errors: Possible error conditions are invalid fid input or if the group is not 
 *         created.
 */ 
herr_t mesh_create_root_group(hid_t fid)
{
    herr_t status = DANU_FAILURE;

    char mesh_root_group_name[] = MESH_ROOT_GROUP_NAME;
    hid_t gid = danu_group_create(fid,mesh_root_group_name);

    if ( H5_ISA_VALID_ID(gid) ) {
        danu_group_close(gid);
        status = 0;
    }
    else {
        DANU_ERROR_MESS("Failed to create the mesh root group");
    }

    return status;
}
 /*
 * Routine: hid_t mesh_open_root_group(hid_t fid)
 * Purpose: Open the root mesh group under teh file fid
 * Description: Open the root mesh group with name MESH_ROOT_GROUP_NAME
 *              under the root group in file fid and returns the HDF5 identifier.  
 *              Errors are raised if fid is an invalid HDF5 identifier
 *              or the routine fails to open the group.
 *
 * Parameters:
 *           fid                IN              HDF5 identifier for the root file
 *                              
 * Returns: An invalid HDF5 identifier if an error occurs. 
 *     
 * Errors: Possible error conditions are invalid fid input or if the group is not 
 *         opened.
 */ 
hid_t mesh_open_root_group(hid_t fid)
{
    char mesh_root_group_name[] = MESH_ROOT_GROUP_NAME;
    hid_t gid = danu_group_open(fid,mesh_root_group_name);

    if ( H5_ISA_INVALID_ID(gid) ) {
        DANU_ERROR_MESS("Failed to open the mesh root group");
    }

    return gid;
}
    
/*
* Routine: hid_t mesh_count(hid_t fid,int *nmeshes)
* Purpose: Returns the number of mesh subgroups under the root file fid.
* Description: Opens the root mesh group and counts the number of HDF5 links under
*              that group. Routine assumes that only mesh groups are under the 
*              root mesh group. Errors are raised if fid is an invalid HDF5 identifier
*              or the routine fails to open the root mesh group. nmeshes is
*              initially set to -1. Calling routines should check the return
*              flag for errors.
*
* Parameters:
*           fid                IN              HDF5 identifier for the root file
*          *nmeshes            OUT             Number of mesh subgroups found 
*                              
* Returns: A negative value if an error occurs, otherwise returns zero.
*     
* Errors: Possible error conditions are invalid fid input, if the pointer to
*         nmeshes or if the group is not opened. 
*/ 
herr_t mesh_count(hid_t fid, int *nmeshes)
{
    herr_t status = DANU_FAILURE;

    hid_t gid = mesh_open_root_group(fid);
    hsize_t nlinks;

    if ( DANU_BAD_PTR(nmeshes) ) {
        DANU_ERROR_MESS("Invalid pointer to nmeshes");
        return status;
    }
    *nmeshes = -1;

    if ( H5_ISA_VALID_ID(gid) ) {
        status = danu_group_get_nlinks(gid,&nlinks);
        *nmeshes = (int) nlinks;
    }
    else {
        DANU_ERROR_MESS("Failed to open mesh root group");
    }

    return status;
}
 /*
* Routine: herr_t mesh_list(hid_t fid,int num, char **meshnames)
* Purpose: Returns the names of all the mesh subgroups found under root file fid.
* Description: Opens the root mesh group and returns all the subgroup names found
*              under that group. Routine assumes that only mesh groups are under
*              the root group. Calling routines should use mesh_count to determine
*              the correct value for num. Errors can occur if the fid is not valid, the
*              root group is not opened or the input for meshnames is not valid.
*              The routine will allocate memory for each name. The calling routine
*              is responsible for deleting the pointers created.
*
* Parameters:
*           fid                IN    HDF5 identifier for the root file
*           num                IN    Number of pointers in meshnames
*           meshnames          OUT   Array of pointers to the mesh names
*                              
* Returns: A negative value if an error occurs, otherwise returns zero.
*     
* Errors: Possible error conditions are invalid fid input, fails to open the 
*         root mesh group, fails to traverse the root mesh group or invalid 
*         pointers to the meshnames.
*         
*/ 
herr_t mesh_list(hid_t fid, int num, char **meshnames, int *num_found)
{
    herr_t status = DANU_FAILURE;

    hid_t gid = mesh_open_root_group(fid);

    if ( H5_ISA_VALID_ID(gid) ) {

      status = danu_group_get_subgroups(gid,num,meshnames,num_found);
      danu_group_close(gid);

    }
    else {
       DANU_ERROR_MESS("Failed to open mesh root group");
    }


    return status;
}
 /*
* Routine: herr_t mesh_exists(hid_t fid, const char *meshname, int *exist)
* Purpose: Determines the existence of meshname under root file fid.
* Description: Opens the root mesh group and searches the subgroup meshname.
*              Sets exist to TRUE if found, FALSE otherwise. The return code
*              will be negative if an error occurs while searching the for the
*              subgroup. The default value of exist is FALSE
*
* Parameters:
*           fid                IN    HDF5 identifier for the root file
*           meshname           IN    Target mesh name
*           exist              OUT   Flag (TRUE or FALSE) indicating if mesh exists
*                              
* Returns: A negative value if an error occurs, otherwise returns zero.
*     
* Errors: Possible error conditions are invalid fid input, fails to open the 
*         root mesh group or fails to traverse the root mesh group.
*         
*/ 
herr_t mesh_exists(hid_t fid, const char *meshname, int *exist)
{
    herr_t status = DANU_FAILURE;

    hid_t gid = mesh_open_root_group(fid);
    hbool_t found;

    if ( DANU_BAD_PTR(exist) ) {
        DANU_ERROR_MESS("Invalid pointer to exist");
        return status;
    }
    *exist = FALSE;

    if ( H5_ISA_VALID_ID(gid) ) {

        found = danu_group_exists(gid,meshname);
        if ( found ) {
            *exist = TRUE;
        }
        danu_group_close(gid);
        status = 0;
    }
    else {
        DANU_ERROR_MESS("Failed to open root mesh group");
    }

    return status;
}

/*
* Routine: herr_t mesh_create(hid_t fid, 
                              const char *meshname,
                              int dim,
                              tmesh_t mtype,
                              telem_t etype,
			      hid_t *mid)
* Purpose: Create a mesh subgroup under root file fid
* Description: Create the mesh subgroup, with all the appropriate subgroups, of
*              mesh type mtype and element type etype.
*              Attempting to create an existing mesh meshname is considered an
*              error. The return code will indicate if the
*              routine create all sub groups and datasets correctly.
*
* Parameters:
*           fid                IN     HDF5 identifier for the root file
*           meshname           IN     New mesh name
*           mtype              IN     STRUCTURED or UNSTRUCTURED
*           etype              IN     Mesh type LINE/TRI/QUAD/TET/HEX
*           mid                OUT    Pointer to mesh HDF5 group identifier
*                              
* Returns: A negative value if an error occurs, otherwise returns zero.
*     
* Errors: Possible error conditions are invalid fid input, fails to open the 
*         root mesh group, fails to create the subgroups for coordinates or connectivity,
*         or the mesh/element type is not valid 
*/
herr_t mesh_create(hid_t fid, 
                   const char * meshname,
                   tmesh_t mtype, 
                   telem_t etype,
		   hid_t *mid)
{
    herr_t status = DANU_FAILURE;

    hid_t gid = mesh_open_root_group(fid);
    int dim;
    int order;

    if ( H5_ISA_INVALID_ID(gid) ) {
        DANU_ERROR_MESS("Failed to open the root mesh group");
        return status;
    }

    if ( DANU_RETURN_FAIL(elem_type_to_order(etype,&order) ) ) {
        DANU_ERROR_MESS("Invalid element type");
        return status;
    }

    *mid = danu_group_create(gid,meshname);

    if ( H5_ISA_VALID_ID(*mid) ) {

	/* Determine the dimension from the element type */
	switch(etype) {
	    case(LINE_ELEM):
		dim = 1;
		break;
	    case(TRI_ELEM):
		dim = 2;
		break;
            case(QUAD_ELEM):
		dim = 2;
		break;
            case(TET_ELEM):
		dim = 3;
		break;
            case(HEX_ELEM):
		dim = 3;
		break;
	    default:
		dim = -1;
		DANU_ERROR_MESS("Invalid element type");
		break;
	}

        if ( dim > 0 ) {
          mesh_set_type(*mid,mtype);
          mesh_set_elementtype(*mid,etype);
          mesh_set_dimension(*mid,dim);
          mesh_set_elem_order(*mid,order);
	  status = DANU_SUCCESS;
	}

    }
    else {
        DANU_ERROR_MESS("Failed to create a mesh subgroup");
    }

    return status;
}
/*
* Routine: herr_t mesh_open(hid_t fid, const char *meshname) 
* Purpose: Open  mesh subgroup under root file fid
* Description: Open the mesh root group under fid, then open the mesh
*              subgroup meshname under the root group. Return teh HDF5
*              identifier for meshname. Routine will check the existence
*              of mesh name before opening the group.
*
* Parameters:
*           fid                IN    HDF5 identifier for the root file
*           meshname           IN    Mesh name
*                              
* Returns: An invalid HDF5 identifer is returned if an error occurs.
*     
* Errors: Possible error conditions are invalid fid input, fails to open the 
*         root mesh group, fails to open the mesh subgroup.
*/
hid_t mesh_open(hid_t fid,const char *meshname)
{
    hid_t mesh = H5I_BADID;
    hid_t gid;
    int exists = FALSE;

    if ( H5_RETURN_OK(mesh_exists(fid,meshname,&exists) ) ) {

        if ( exists ) {
            gid = mesh_open_root_group(fid);
            mesh = danu_group_open(gid,meshname);
            danu_group_close(gid);
        }
        else {
            DANU_ERROR_MESS("Attempting to open mesh that does not exist");
        }

    }
    else {
        DANU_ERROR_MESS("Failed to stat the mesh group");
    }

    return mesh;
}
/*
* Routine: herr_t mesh_add_unstructured(hid_t fid, 
                                        const char *meshname,
                                        int elemorder,
                                        int dim,
			                hid_t *mid)
* Purpose: Create an unstructured mesh subgroup under root file fid
* Description: Create the mesh subgroup, with all the appropriate subgroups, for 
*              an unstructured mesh type with element order elemorder and dimension dim.
*              Attempting to create an existing
*              mesh meshname is considered an error. The return code will indicate if the
*              routine create all sub groups and datasets correctly.
*
* Parameters:
*           fid                IN     HDF5 identifier for the root file
*           meshname           IN     New mesh name
*           elemorder          IN     Element order 2,3,4,8
*           dim                IN     Dimension of the mesh
*           mid                OUT    Pointer to mesh HDF5 group identifier
*                              
* Returns: A negative value if an error occurs, otherwise returns zero.
*     
* Errors: Possible error conditions are invalid fid input, fails to open the 
*         root mesh group, fails to create the subgroups for coordinates or connectivity,
*         or the mesh/element type is not valid 
*/
herr_t mesh_add_unstructured(hid_t fid,
	                     const char *meshname,
			     int elemorder,
			     int dim,
			     hid_t *mid)
{
    herr_t status = DANU_FAILURE;
    telem_t elem;
    tmesh_t type = UNSTRUCTURED_MESH;
    hid_t id;
    herr_t err;

    if ( DIM_IS_INVALID(dim) ) {
	DANU_ERROR_MESS("Invalid dimension argument");
	*mid = H5I_INVALID_HID;
	return status;
    }

    if ( EORDER_IS_INVALID(elemorder) ) {
	DANU_ERROR_MESS("Invalid element order argument");
	*mid = H5I_INVALID_HID;
	return status;
    }

    err = elem_order_to_type(dim, elemorder, &elem);

    if ( DANU_RETURN_FAIL(err) ) {
	DANU_ERROR_MESS("Mismatched dimension and element order values");
	*mid = H5I_INVALID_HID;
	return status;
    }

    status = mesh_create(fid,meshname,type,elem,mid);

    return status;
}
/*
* Routine: herr_t mesh_create_coordinates(hid_t mid, int dim, int nnodes) 
* Purpose: Create the coordinate dataset under mesh mid
* Description: Create the 2D dataset of size dim x nnodes under the mesh group
*              mid. Both the dim and nnodes are checked for strictly positive
*              values before the dataset is created. The dataset is closed
*              before the routine returns. The return flag will be neagtive
*              if an error is detected.
*
* Parameters:
*           mid                IN    HDF5 identifier for the mesh group
*           dim                IN    Dimension of the mesh
*           nnodes             IN    Number of nodes (vertices) 
*                              
* Returns: An invalid HDF5 identifier is returned if an error occurs.
*     
* Errors: Possible error conditions are invalid mid input, either dim or
*         nnodes is not strictly positive or the routine fails to 
*         create the dataset.
*/
hid_t mesh_create_coordinates(hid_t mid, int dim, int nnodes)
{
    herr_t status = DANU_FAILURE;

    hsize_t size[COORD_SIZE_DIM];
    hid_t cid;
    char coord_name[]=MESH_NODAL_DATA_NAME;

    if ( dim <= 0 || dim > 3 ) {
        DANU_ERROR_MESS("Invalid dimension value");
        return status;
    }

    if ( nnodes <= 0 ) {
        DANU_ERROR_MESS("Invalid number of nodes");
        return status;
    }

    /* Define the size of the dataset */
    size[COORD_SIZE_NNODES_IDX] = (hsize_t) nnodes;
    size[COORD_SIZE_DIM_IDX]    = (hsize_t) dim;

    cid = danu_dataset_create_double(mid,coord_name,COORD_SIZE_DIM,size);

    return cid;
}
/*
* Routine: hid_t mesh_open_coordinates(hid_t mid) 
* Purpose: Open the coordinate dataset under mesh mid
* Description: Open the 2D dataset of size dim x nnodes under the mesh group
*              mid. The HDF5 identifier for this dataset is returned upon
*              successful completion of the routine. An invalid identifier
*              indicates an error occured.
*
* Parameters:
*           mid                IN    HDF5 identifier for the mesh group
*                              
* Returns: Returns the HDF5 identifier of the coordinate dataset. An invalid
*          HDF5 identifier is returned if an error occurs.
*     
* Errors: Possible error conditions are invalid mid input or failure
*         to open the dataset.
*/
hid_t mesh_open_coordinates(hid_t mid)
{
    hid_t coord = H5I_BADID;

    char coord_name[] = MESH_NODAL_DATA_NAME;

    coord = danu_dataset_open(mid,coord_name);

    if ( H5_ISA_INVALID_ID(coord) ) {
        DANU_ERROR_MESS("Failed to open the coordinate dataset");
    }

    return coord;
}
/*
* Routine: herr_t mesh_close_coordinates(hid_t coord_id) 
* Purpose: Close the coordinate dataset coord_id
* Description: Close the coordinate dataset. Will return negative value if
*              error occurs or the input is invalid.
*
* Parameters:
*           coord_id                IN    HDF5 identifier for the mesh group
*                              
* Returns: A negative value if error is detected, otherwise returns 0.
*     
* Errors: Possible error conditions are invalid mid input or failure
*         to close the dataset.
*/
herr_t mesh_close_coordinates(hid_t coord_id)
{
    herr_t status = DANU_FAILURE;

    if ( H5_ISA_VALID_ID(coord_id) ) {
      status = danu_dataset_close(coord_id);
    }
    else {
      DANU_ERROR_MESS("Invalid coordinate H5 id. Can not close dataset.");
    }

    return status;
}
/*
* Routine: herr_t mesh_write_coordinates(hid_t mid, int nnodes, double * x, double *y, double *z) 
* Purpose: Write the coordinate dataset to the mesh group mid
* Description: Write the 2D dataset of size nnodes x dim under the mesh group
*              mid. The routine accepts y-,z-coordinates as fixed but optional arguments.
*              Calling routines that are writing 1d and 2d meshes must pass in
*              NULLs for y and z coordinate vectors. The dimension of the mesh
*              must match the number of non-NULL vectors.
*              The number of nnodes is set as an attribute of the mesh
*              group mid. The routine closes the dataset once it is written.
*              If an error occurs the functions returns a negative value.
*
* Parameters:
*           mid                IN    HDF5 identifier for the mesh group
*           nnodes             IN    Number of nodes (length of the coordinate arrays)
*           x                  IN    x-coordinates
*           y                  IN    y-coordinates (OPTIONAL)
*           z                  IN    z-coordinates (OPTIONAL)
*                              
* Returns: A negative value if an error occurs, otherwise the return value is zero.
*     
* Errors: Possible error conditions are invalid mid input, invalid nnodes,
*         dimension does not match the number arguments or failure
*         to open and write the dataset.
*/
herr_t mesh_write_coordinates(hid_t mid, int nnodes, double * x, double *y, double *z)
{
    herr_t status = DANU_FAILURE;

    int i,dim, ndim;
    double *coordinate[3];
    int coordinate_index[3];

    hid_t data;
    dslab_t *slab;
    hsize_t *offset, *size, n_nodes;
    herr_t loop_stat;

    if ( DANU_RETURN_FAIL(mesh_get_dimension(mid,&dim)) ) {
        DANU_ERROR_MESS("Failed to determine the mesh dimension");
        return status;
    }

    if ( nnodes <= 0 ) {
        DANU_ERROR_MESS("Invalid number of nodes");
        return status;
    }

    if ( DANU_BAD_PTR(x) ) {
        DANU_ERROR_MESS("Invalid x coordinate pointer");
        return status;
    }
    coordinate[COORD_X_IDX] = x;
    coordinate_index[0] = COORD_X_IDX;

    /* Now check other coordinate vectors */
    if ( dim >= 2 && DANU_BAD_PTR(y) ) {
        DANU_ERROR_MESS("Invalid y coordinate pointer");
        return status;
    }
    coordinate[COORD_Y_IDX] = y;
    coordinate_index[1] = COORD_Y_IDX;

    if ( dim == 3 && DANU_BAD_PTR(z) ) {
        DANU_ERROR_MESS("Invalid z coordinate pointer");
        return status;
    }
    coordinate[COORD_Z_IDX] = z;
    coordinate_index[2] = COORD_Z_IDX;

    /* Compare the dim with the number of valid pointers */
    ndim = 0;
    for(i=0; i < 3; i++) {
      if ( coordinate[i] != NULL ) {
	ndim++;
      }
    }

    if ( ndim != dim ) {
      DANU_ERROR_MESS("Mesh dimension mismatch detected");
      return status;
    }

    /* Create the dataset and open the dataset */
    if ( DANU_RETURN_FAIL(mesh_create_coordinates(mid,dim,nnodes) ) ) {
        DANU_ERROR_MESS("Failed to create the coordinate dataset");
        return status;
    }
    data = mesh_open_coordinates(mid);

    if ( H5_ISA_VALID_ID(data) ) {

        n_nodes = (hsize_t) nnodes;

        /* Define the hyper-slab */ 
        slab   = danu_slab_alloc(data);
        offset = DANU_MALLOC(hsize_t, COORD_SIZE_DIM);
        size   = DANU_MALLOC(hsize_t, COORD_SIZE_DIM);

        size[COORD_SIZE_NNODES_IDX] = n_nodes;
        size[COORD_SIZE_DIM_IDX]    = 1;

        for(i=0;i<COORD_SIZE_DIM;i++)
            offset[i] = 0;


        /* Loop through the coordinates and write each coordinate vector */
        i=0;
        loop_stat = 0;
        while( (i<dim) && (loop_stat >= 0 ) ) {

            if ( coordinate[i] != NULL ) {
                offset[COORD_SIZE_DIM_IDX] = coordinate_index[i];

                danu_slab_contiguous(slab,offset,size);
                loop_stat = danu_slab_select(data,slab);
                if ( DANU_RETURN_OK(loop_stat) ) {
                    loop_stat = danu_dataset_write(data,slab,H5T_NATIVE_DOUBLE,1,&n_nodes,coordinate[i]);
                }
                offset[COORD_SIZE_DIM_IDX] = 0;
            }
            i++;
        }

        /* Set the nnodes attribute if ALL coordinates writes returned with OK stat */
        if ( i == dim && ( loop_stat >= 0) ) {
            status = loop_stat;
            mesh_set_nnodes(mid,nnodes);
        }
        else {
            DANU_ERROR_MESS("Failed to write coordinate data");
        }

        DANU_FREE(offset);
        DANU_FREE(size);
        danu_slab_free(slab);
	mesh_close_coordinates(data);
        
    }
    else {
        DANU_ERROR_MESS("Can not open coordinate dataset");
    }

    return status;

}
/*
* Routine: herr_t mesh_write_coordinates_1d(hid_t mid, int nnodes, double * x) 
* Purpose: Write the coordinate dataset to the 1D mesh group mid
* Description: A wrapper routine to mesh_write_coordinates that passes NULL
*              values for the y and z coordinates. See mesh_write_coordinates
*              for more information.
*
* Parameters:
*           mid                IN    HDF5 identifier for the mesh group
*           nnodes             IN    Number of nodes (length of the coordinate arrays)
*           x                  IN    x-coordinates
*                              
* Returns: A negative value if an error occurs, otherwise the return value is zero.
*     
* Errors: Possible error conditions are invalid mid input, invalid nnodes,
*         dimension does not match the number arguments or failure
*         to open and write the dataset.
*/
herr_t mesh_write_coordinates_1d(hid_t mid, int nnodes, double *x)
{
  int dim;
  herr_t status = mesh_get_dimension(mid,&dim);

  if ( DANU_RETURN_FAIL(status) ) {
    DANU_ERROR_MESS("Failed to define the mesh dimension");
    return status;
  }

  if ( dim != 1 ) {
    DANU_ERROR_MESS("Mesh dimension does not match coordinate write size");
    return DANU_FAILURE;
  }

  status = mesh_write_coordinates(mid,nnodes,x,NULL,NULL);

  return status;

}
/*
* Routine: herr_t mesh_write_coordinates_2d(hid_t mid, int nnodes, double * x, double *y) 
* Purpose: Write the coordinate dataset to the 2D mesh group mid
* Description: A wrapper routine to mesh_write_coordinates that passes NULL
*              value the z coordinates. See mesh_write_coordinates
*              for more information.
*
* Parameters:
*           mid                IN    HDF5 identifier for the mesh group
*           nnodes             IN    Number of nodes (length of the coordinate arrays)
*           x                  IN    x-coordinates
*           y                  IN    y-coordinates
*                              
* Returns: A negative value if an error occurs, otherwise the return value is zero.
*     
* Errors: Possible error conditions are invalid mid input, invalid nnodes,
*         dimension does not match the number arguments or failure
*         to open and write the dataset.
*/
herr_t mesh_write_coordinates_2d(hid_t mid, int nnodes, double *x,double *y)
{
  int dim;
  herr_t status = mesh_get_dimension(mid,&dim);

  if ( DANU_RETURN_FAIL(status) ) {
    DANU_ERROR_MESS("Failed to define the mesh dimension");
    return status;
  }

  if ( dim != 2 ) {
    DANU_ERROR_MESS("Mesh dimension does not match coordinate write size");
    return DANU_FAILURE;
  }

  status = mesh_write_coordinates(mid,nnodes,x,y,NULL);

  return status;

}

/*
* Routine: herr_t mesh_read_coordinates(hid_t mid, double * x, double *y, double *z) 
* Purpose: Read the coordinate dataset from the mesh group mid
* Description: Read the x, (possibly y and z) coordinates from the 2D dataset containing the
*              coordinate data under mesh group mid. The routine assumes that the 
*              calling routine has provided appropriate sized memory allocations
*              for the x,y, and z pointers. The dimension of the mesh is checked against
*              the number of pointers passed in and an error is raised if they do not
*              match. For meshes that are not 3D, the calling routine should pass
*              in a NULL pointer for y or z.
*
* Parameters:
*           mid                IN    HDF5 identifier for the mesh group
*           x                  OUT   x-coordinates
*           y                  OUT   y-coordinates (OPTIONAL, ignored if NULL)
*           z                  OUT   z-coordinates (OPTIONAL, ignored if NULL)
*                              
* Returns: A negative value if an error occurs, otherwise the return value is zero.
*     
* Errors: Possible error conditions are invalid mid input, invalid pointers,
*         dimension does not match the number arguments or failure
*         to open and read the dataset.
*/
herr_t mesh_read_coordinates(hid_t mid, double *x, double *y, double *z)
{
    herr_t status = DANU_FAILURE;

    hid_t data;

    double *coordinate[3];
    int     coordinate_index[3];
    
    int i,ndim,dim,nnodes;
    dslab_t *slab;
    hsize_t *offset, *size, n_nodes;
    herr_t loop_stat;

    /* Read the dimension ... routine will check the id */
    if (DANU_RETURN_FAIL(mesh_get_dimension(mid,&dim) ) ) {
        DANU_ERROR_MESS("Failed to read mesh dimension");
        return status;
    }

    if ( dim <=0 || dim > 3 ) {
        DANU_ERROR_MESS("Invalid mesh dimension .... invalid mesh");
        return status;
    }


    /* Read in the number of nodes...need this for the hyperslab */
    if ( DANU_RETURN_FAIL(mesh_get_nnodes(mid,&nnodes) ) ){
        DANU_ERROR_MESS("Failed to read the number of nodes");
        return status;
    }

    /* Check the pointers and define the coordinate placeholder arrays */
    if ( DANU_BAD_PTR(x) ) {
        DANU_ERROR_MESS("Invalid x-coordinate pointer");
        return status;
    }
    coordinate[COORD_X_IDX] = x;
    coordinate_index[0] = COORD_X_IDX;


    if ( dim >= 2 && DANU_BAD_PTR(y) ) {
        DANU_ERROR_MESS("Invalid y-coordinate pointer");
        return status;
    }
    coordinate[COORD_Y_IDX] = y;
    coordinate_index[1] = COORD_Y_IDX;

    if ( dim == 3 && DANU_BAD_PTR(z) ) {
        DANU_ERROR_MESS("invalid z-coordinate pointer");
        return status;
    }
    coordinate[COORD_Z_IDX] = z;
    coordinate_index[2] = COORD_Z_IDX;

    /* Compare mesh dimension attribute to number of valid pointers */
    ndim = 0;
    for( i = 0; i < 3; i++) {
      if ( coordinate[i] != NULL ) {
	ndim++;
      }
    }
    if ( ndim != dim ) {
      DANU_ERROR_MESS("Mesh read dimension mismatch");
      return status;
    }


    /* Now ready to open and read */
    data = mesh_open_coordinates(mid);
    if ( H5_ISA_VALID_ID(data) ) {

        n_nodes = (hsize_t ) nnodes;

        /* Define the hyperslab */
        slab   = danu_slab_alloc(data);
        offset = DANU_MALLOC(hsize_t, COORD_SIZE_DIM);
        size   = DANU_MALLOC(hsize_t, COORD_SIZE_DIM);

        /* Set the size and offset arrays */
        size[COORD_SIZE_NNODES_IDX] = nnodes;
        size[COORD_SIZE_DIM_IDX]    = 1;

        for(i=0;i<COORD_SIZE_DIM;i++)
            offset[i] = 0;

        /* Loop through each coordinate and read the data */ 
         i=0;
        loop_stat = 0;
        while( (i<dim) && (loop_stat >= 0 ) ) {

            if ( coordinate[i] != NULL ) {
                offset[COORD_SIZE_DIM_IDX] = coordinate_index[i];

                danu_slab_contiguous(slab,offset,size);
                loop_stat = danu_slab_select(data,slab);
                if ( DANU_RETURN_OK(loop_stat) ) {
                    loop_stat = danu_dataset_read(data,slab,H5T_NATIVE_DOUBLE,1,&n_nodes,coordinate[i]);
                }
                offset[COORD_SIZE_DIM_IDX] = 0;
            }
            i++;
        }

        /* Set the return status to OK if the loop exited without an error */
        if ( i == dim && ( loop_stat >= 0) ) {
            status = loop_stat;
        }
        else {
            DANU_ERROR_MESS("Failed to read coordinate data");
        }

        DANU_FREE(offset);
        DANU_FREE(size);
        danu_slab_free(slab);
        
    }
    else {
        DANU_ERROR_MESS("Can not open coordinate dataset");
    }

    return status;
}
/*
* Routine: herr_t mesh_read_coordinates_1d(hid_t mid, double * x) 
* Purpose: Read the coordinate dataset of a 1D mesh into buffer x
* Description: A wrapper routine to mesh_read_coordinates that passes NULL
*              values for the y and z coordinates. See mesh_write_coordinates
*              for more information.
*
* Parameters:
*           mid                IN    HDF5 identifier for the mesh group
*           x                  OUT   x-coordinates
*                              
* Returns: A negative value if an error occurs, otherwise the return value is zero.
*     
* Errors: Possible error conditions are invalid mid input, invalid nnodes,
*         dimension does not match the number arguments or failure
*         to open and read the dataset.
*/
herr_t mesh_read_coordinates_1d(hid_t mid, double *x)
{
  int dim;
  herr_t status = mesh_get_dimension(mid,&dim);

  if ( DANU_RETURN_FAIL(status) ) {
    DANU_ERROR_MESS("Failed to define the mesh dimension");
    return status;
  }

  if ( dim != 1 ) {
    DANU_ERROR_MESS("Mesh dimension does not match coordinate read size");
    return DANU_FAILURE;
  }

  status = mesh_read_coordinates(mid,x,NULL,NULL);

  return status;

}
/*
* Routine: herr_t mesh_read_coordinates_2d(hid_t mid, int nnodes, double * x, double *y) 
* Purpose: Read the coordinate dataset of a 2D mesh to buffers x and y 
* Description: A wrapper routine to mesh_read_coordinates that passes NULL
*              value the z coordinates. See mesh_write_coordinates
*              for more information.
*
* Parameters:
*           mid                IN    HDF5 identifier for the mesh group
*           x                  OUT   x-coordinates
*           y                  OUT   y-coordinates
*                              
* Returns: A negative value if an error occurs, otherwise the return value is zero.
*     
* Errors: Possible error conditions are invalid mid input, 
*         dimension does not match the number arguments or failure
*         to open and read the dataset.
*/
herr_t mesh_read_coordinates_2d(hid_t mid, double *x,double *y)
{
  int dim;
  herr_t status = mesh_get_dimension(mid,&dim);

  if ( DANU_RETURN_FAIL(status) ) {
    DANU_ERROR_MESS("Failed to define the mesh dimension");
    return status;
  }

  if ( dim != 2 ) {
    DANU_ERROR_MESS("Mesh dimension does not match coordinate read size");
    return DANU_FAILURE;
  }

  status = mesh_read_coordinates(mid,x,y,NULL);

  return status;

}

/*
* Routine: herr_t mesh_read_coordinates_by_index(hid_t mid, int index, int nnodes, double * buff) 
* Purpose: Read single index coordinate dataset of a mesh to buffers buff 
* Description: Reads a single index dataset from the coordinates. Index is zero-based
*              and can not exceed the dimensions of the mesh. 
*
* Parameters:
*           mid                IN    HDF5 identifier for the mesh group
*           index              IN    index of the coordinate dataset to read
*           buf                OUT   data buffer
*                              
* Returns: A negative value if an error occurs, otherwise the return value is zero.
*     
* Errors: Possible error conditions are invalid mid input, index exceeds the 
*         dimensions of the mesh, or fail to open the mesh coordinate group.
*/
herr_t mesh_read_coordinates_byindex(hid_t mid, int index, double * buf)
{
  int dim, nnodes, i;
  herr_t status = DANU_FAILURE;
  hid_t coord;
  dslab_t *slab;
  hsize_t *offset, *size, n_nodes;

  /* Check the input */
  if ( DANU_RETURN_FAIL(mesh_get_dimension(mid,&dim))) {
     DANU_ERROR_MESS("Failed to define mesh dimensions");
     return DANU_FAILURE;
  }

  if ( DANU_RETURN_FAIL(mesh_get_nnodes(mid,&nnodes))) {
     DANU_ERROR_MESS("Failed to define number of nodes");
     return DANU_FAILURE;
  }

  if ( index < 0 || index >= dim ) {
    DANU_ERROR_MESS("Invalid index");
    return DANU_FAILURE;
  }

  if ( DANU_BAD_PTR(buf) ) {
    DANU_ERROR_MESS("Invalid buffer pointer");
    return DANU_FAILURE;
  }

  /* Define the data slab */
  coord = mesh_open_coordinates(mid);
  if ( H5_ISA_VALID_ID(coord) ) {

    n_nodes = (hsize_t) nnodes;

    slab = danu_slab_alloc(coord);
    offset = DANU_MALLOC(hsize_t, COORD_SIZE_DIM);
    size   = DANU_MALLOC(hsize_t, COORD_SIZE_DIM);

    size[COORD_SIZE_NNODES_IDX] = n_nodes;
    size[COORD_SIZE_DIM_IDX]    = 1;

    for(i=0;i<COORD_SIZE_DIM;i++)
      offset[i] = 0;
    offset[COORD_SIZE_DIM_IDX] = index;

    danu_slab_contiguous(slab,offset,size);
    status = danu_dataset_read(coord,slab,H5T_NATIVE_DOUBLE,1,&n_nodes,buf);

    danu_slab_free(slab);
    DANU_FREE(offset);
    DANU_FREE(size);
    mesh_close_coordinates(coord);

  }
  else {
    DANU_ERROR_MESS("Failed to open mesh coordinate group");
  }

  return status;

}

/*
* Routine: hid_t mesh_create_connectivity(hid_t mid, int nelem) 
* Purpose: Create the connectivity dataset under mesh mid and return the connectivity id.
* Description: Create connectivity dataset under the mesh group
*              mid. Note this dataset is ONLY required for UNSTRUCTURED mesh types
*
* Parameters:
*           mid                IN    HDF5 identifier for the mesh group
*           nelem              IN    Number of elements
*                              
* Returns: Returns an HDF5 identifier mapped to the connectivity dataset. 
*     
* Errors: Possible error conditions are invalid mid input, or nelem
*         is not strictly positive or the routine fails to 
*         create the dataset.
*/
herr_t mesh_create_connectivity(hid_t mid, int nelem)
{

    hid_t cid = H5I_INVALID_HID;
    char  conn_name[] = MESH_CONNECT_DATA_NAME;
    hsize_t size[CONNECT_SIZE_DIM];
    int elem_order; 

    /* Check input */
    if ( nelem <= 0 ) {
      DANU_ERROR_MESS("Invalid number of mesh elements.")
       return cid;
    }

    if ( DANU_RETURN_FAIL(mesh_get_elem_order(mid,&elem_order) ) ) {
      DANU_ERROR_MESS("Failed to determine mesh element order.");
      return cid;
    }


    /* Define the dimensions of the dataset */
    size[CONNECT_SIZE_NELEM_IDX] = (hsize_t) nelem;
    size[CONNECT_SIZE_ELEMORDER_IDX] = (hsize_t) elem_order;

    cid = danu_dataset_create_int(mid,conn_name,CONNECT_SIZE_DIM,size);
    
    /* Default offset is 0, can be changed with another set call */
    if ( DANU_RETURN_FAIL(danu_set_offset(cid,0)) ) {
      DANU_ERROR_MESS("Failed to write the offset for connectivity dataset");
    }

    return cid;
}
/*
* Routine: hid_t mesh_open_connectivity(hid_t mid) 
* Purpose: Open the connectivity dataset under mesh mid
* Description: Open the connectivity dataset under the mesh group
*              mid. The HDf5 identifier for this dataset is returned upon
*              successful completion of the routine. An invalid identifier
*              indicates an error occured.
*
* Parameters:
*           mid                IN    HDF5 identifier for the mesh group
*                              
* Returns: Returns the HDF5 idenitifier for the connectivity group.
*          An invalid HDF5 identifier is returned if an error occurs.
*     
* Errors: Possible error conditions are invalid mid input or failure
*         to open the group.
*/
hid_t mesh_open_connectivity(hid_t mid)
{
    hid_t conn = H5I_BADID;

    char conn_name[] = MESH_CONNECT_DATA_NAME;

    conn = danu_dataset_open(mid,conn_name);

    if ( H5_ISA_INVALID_ID(conn) ) {
        DANU_ERROR_MESS("Failed to open the connectivity dataset");
    }

    return conn;
}
/*
* Routine: hid_t mesh_write_connectivity(hid_t mid, int nelem, const int *data) 
* Purpose: Write connectivity data (nelem x elem_order) under the connectivity dataset of mesh group mid
* Description: Open the connnectivity dataset under mesh group mid and write a 2D dataset sized
*              nelems x elem_order.  
*
* Parameters:
*           mid                IN    HDF5 identifier for the mesh group
*           nelem              IN    Number of elements
*           data               IN    Data (nelem x elem_order) to write into dataset
*                              
* Returns: Returns a negative value if an error occurs, otherwise returns a zero.
*     
* Errors: Possible error conditions are invalid HDF5 identifier, invalid pointers, the
*         dataset already exists, fail to open the mesh or connectivity group and fail
*         to write the data.
*/
herr_t mesh_write_connectivity(hid_t mid, int nelem, const int *data)
{
    
    herr_t status = DANU_FAILURE;
    hsize_t size[CONNECT_SIZE_DIM];
    hid_t id;
    int elem_order;

    /* Check the input */
    if ( DANU_BAD_PTR(data) ) {
      DANU_ERROR_MESS("Invalid pointer to connectivity data.");
      return status;
    }

    if ( nelem < 0 ) {
      DANU_ERROR_MESS("Invalid number of elements.");
      return status;
    }

    if ( DANU_RETURN_OK(mesh_get_elem_order(mid,&elem_order)) ) {
      size[CONNECT_SIZE_NELEM_IDX]=(hsize_t)nelem;
      size[CONNECT_SIZE_ELEMORDER_IDX]=(hsize_t)elem_order;
      status = danu_data_write_int(mid,MESH_CONNECT_DATA_NAME,CONNECT_SIZE_DIM,size,data);

      if ( status == DANU_SUCCESS ) {
	id = mesh_open_connectivity(mid);
	status &= danu_set_offset(id,0);
	danu_dataset_close(id);
      }

      status &= mesh_set_nelem(mid,nelem);

    }

    return status;

}

/*
* Routine: hid_t mesh_read_connectivity(hid_t mid, int *data) 
* Purpose: Read connectivity data under mesh group mid
* Description: Open the connnectivity dataset under mesh group mid and read
*              the data. Routine assumes the pointer data is sufficently 
*              sized to hold the connectivity data.
*
* Parameters:
*           mid                IN    HDF5 identifier for the mesh group
*           data               OUT   Pointer to data
*                              
* Returns: Returns a negative value if an error occurs, otherwise returns a zero.
*     
* Errors: Possible error conditions are invalid HDF5 identifier, invalid pointers, the
*         dataset does not exist, fail to open the mesh or connectivity group and fail
*         to read the data.
*/
herr_t mesh_read_connectivity(hid_t mid, int *data)
{
    
    herr_t status = DANU_FAILURE;
    hsize_t size[CONNECT_SIZE_DIM];
    int elem_order, nelem;

    /* Check the input */
    if ( DANU_BAD_PTR(data) ) {
        DANU_ERROR_MESS("Invalid pointer to connectivity data");
        return status;
    }

    status = mesh_get_elem_order(mid,&elem_order);
    status &= mesh_get_nelem(mid,&nelem);
    if ( DANU_RETURN_OK(status) ) {
      size[CONNECT_SIZE_NELEM_IDX]=(hsize_t)nelem;
      size[CONNECT_SIZE_ELEMORDER_IDX]=(hsize_t)elem_order;
      status = danu_data_read_int(mid,MESH_CONNECT_DATA_NAME,2,size,data);
    }
    else {
      DANU_ERROR_MESS("Failed to retrieve mesh connectivity attributes.");
    }

    return status;

}
/*
* Routine: herr_t mesh_connectivity_size(hid_t mid, int *nelem, int * elem_order)
* Purpose: Return the size of an existing connectivity dataset.
* Description: Open an existing connectivity data set and set the
*               the number of elements and element order. Error is raised
*               if the dataset does not exist, mid is an invalid id or
*               pointers are NULL. Routine is a wrapper to mesh_get_* calls.
*
* Parameters:
*           mid                IN    HDF5 identifier for the mesh group
*           *nelem             OUT   Number of elements
*           *elem_order        OUT   Element order (number of vertices)
*                              
* Returns: Returns a negative value if an error occurs, otherwise returns a zero.
*     
* Errors: Possible error conditions are invalid HDF5 identifier, invalid pointers, the
*         dataset does not exist, fail to open the mesh or connectivity dataset.
*/
herr_t mesh_connectivity_size(hid_t mid, int *nelem, int *elem_order)
{
  herr_t stat = DANU_FAILURE;

  stat = mesh_get_nelem(mid,nelem);
  stat &= mesh_get_elem_order(mid,elem_order);
  return stat;
}
/*
* Routine: herr_t mesh_connectivity_set_offset(hid_t mid, int offset)
* Purpose: Set the offset attribute in the connectivity dataset
* Description: Open an existing connectivity data set and set the
*              offset attribute. 
*
* Parameters:
*           mid               IN    HDF5 identifier for the mesh group
*           offset            IN    Offset value
*                              
* Returns: Returns a negative value if an error occurs, otherwise returns a zero.
*     
* Errors: Possible error conditions are invalid HDF5 identifier, invalid pointers, the
*         dataset does not exist, fail to open the mesh or connectivity dataset.
*/
herr_t mesh_connectivity_set_offset(hid_t mid, int offset)
{
  herr_t stat = DANU_FAILURE;
  hid_t id;

  id = mesh_open_connectivity(mid);
  if (H5_ISA_VALID_ID(id)) {
    stat = danu_set_offset(id,offset);
  }
  else {
    DANU_ERROR_MESS("Failed to open Connectivity dataset");
  }
 
  return stat;
}
/*
* Routine: herr_t mesh_connectivity_get_offset(hid_t mid, int *offset)
* Purpose: Fetch the offset attribute in the connectivity dataset
* Description: Open an existing connectivity data set and get the
*              offset attribute value. Calling routine should check 
*              the return before using the offset value.
*
* Parameters:
*           mid              IN    HDF5 identifier for the mesh group
*           *offset          OUT    Pointer to the offset value
*                              
* Returns: Returns a negative value if an error occurs, otherwise returns a zero.
*     
* Errors: Possible error conditions are invalid HDF5 identifier, invalid pointers, the
*         dataset does not exist, fail to open the mesh or connectivity dataset.
*/
herr_t mesh_connectivity_get_offset(hid_t mid, int *offset)
{
  herr_t stat = DANU_FAILURE;
  hid_t id;

  if (DANU_BAD_PTR(offset)) {
    DANU_ERROR_MESS("Invalid offset pointer");
    return stat;
  }

  /* Set to dummy value */
  *offset=0xFFFF;

  id = mesh_open_connectivity(mid);
  if (H5_ISA_VALID_ID(id)) {
    stat = danu_get_offset(id,offset);
  }
  else {
    DANU_ERROR_MESS("Failed to open Connectivity dataset");
  }
 
  return stat;
}


/*
* Routine: herr_t mesh_set_type(hid_t mid, tmesh_t type) 
* Purpose: Set the mesh type to STRUCTURED or UNSTRUCTURED
* Description: Write the mesh type as a string attribute to the mesh
*              subgroup mid. 
*
* Parameters:
*           mid                IN    HDF5 identifier for the mesh group
*           type               IN    TOUT Mesh type (STRUCTURED or UNSTRUCTURED)
*                              
* Returns: A neagtive value is returned if an error occurs, otherwise returns a zero
*     
* Errors: Possible error conditions are invalid mid input or invalid mesh type
*       
*/
herr_t mesh_set_type(hid_t mid, tmesh_t type)
{
    herr_t status = DANU_FAILURE;

    char attr_label[] = MESH_TYPE_ATTR_NAME;
    char buffer[64];

    if ( DANU_RETURN_OK(mesh_type_to_name(type,buffer)) ) {
        status = danu_attr_write_string(mid,attr_label,buffer);
    }
    else {
        DANU_ERROR_MESS("Failed to write mesh type attribute");
    }

    return status;
}
/*
* Routine: herr_t mesh_get_type(hid_t mid, tmesh_t * type) 
* Purpose: Get the mesh type of mesh mid.
* Description: Read the string attribute MESH_TYPE_ATTR_NAME of HDF5 group mid.
*              This string is then converted to the enumerated type tmesh_t.
*
* Parameters:
*           mid                IN     HDF5 identifier for the mesh group
*           type               OUT    TOUT Mesh type (STRUCTURED or UNSTRUCTURED)
*                              
* Returns: A neagtive value is returned if an error occurs, otherwise returns a zero
*     
* Errors: Possible error conditions are invalid mid input or fail to read the attribute.
*       
*/
herr_t mesh_get_type(hid_t mid, tmesh_t *type)
{
    herr_t status = DANU_FAILURE;

    char attr_label[] = MESH_TYPE_ATTR_NAME;
    char buffer[64];

    /* Check input */
    if ( DANU_BAD_PTR(type) ) {
        DANU_ERROR_MESS("Invalid type pointer");
        return status;
    }

    if ( H5_RETURN_FAIL(danu_attr_read_string(mid,attr_label,buffer,64) ) ) {
        DANU_ERROR_MESS("Failed to read mesh attribute for mesh type");
    }
    else {
        status = mesh_string_to_type(buffer,type);
    }

    return status;

}

/*
* Routine: herr_t mesh_set_elementtype(hid_t mid, telem_t type) 
* Purpose: Set the mesh element type to TRI, QUAD, TET, HEX 
* Description: Write the mesh element type as a string attribute to the mesh
*              subgroup mid.  
*
* Parameters:
*           mid                IN    HDF5 identifier for the mesh group
*           type               IN    TOUT Mesh type (enumerated type)
*                              
* Returns: A neagtive value is returned if an error occurs, otherwise returns a zero
*     
* Errors: Possible error conditions are invalid mid input or invalid element type
*       
*/
herr_t mesh_set_elementtype(hid_t mid, telem_t type)
{
    herr_t status = DANU_FAILURE;

    char attr_label[] = MESH_ETYPE_ATTR_NAME;
    char buffer[64];

    if ( DANU_RETURN_OK(elem_type_to_name(type,buffer)) ) {
        status = danu_attr_write_string(mid,attr_label,buffer);
    }
    else {
        DANU_ERROR_MESS("Failed to write mesh element type attribute");
    }

    return status;
}
/*
* Routine: herr_t mesh_get_elementtype(hid_t mid, telem_t * type) 
* Purpose: Get the mesh element type of mesh mid.
* Description: Read the string attribute MESH_ETYPE_ATTR_NAME of HDF5 group mid.
*              This string is then converted to the enumerated type telem_t.
*
* Parameters:
*           mid                IN     HDF5 identifier for the mesh group
*           type               OUT    TOUT Mesh type (TRI, QUAD, TET and HEX)
*                              
* Returns: A neagtive value is returned if an error occurs, otherwise returns a zero
*     
* Errors: Possible error conditions are invalid mid input or fail to read the attribute.
*       
*/
herr_t mesh_get_elementtype(hid_t mid, telem_t *type)
{
    herr_t status = DANU_FAILURE;

    char attr_label[] = MESH_ETYPE_ATTR_NAME;
    char buffer[64];

    /* Check input */
    if ( DANU_BAD_PTR(type) ) {
        DANU_ERROR_MESS("Invalid type pointer");
        return status;
    }
    *type = INVALID_ELEM;

    if ( H5_RETURN_FAIL(danu_attr_read_string(mid,attr_label,buffer,64)) ) {
        DANU_ERROR_MESS("Failed to read mesh attribute for mesh element type");
    }
    else {
        status = elem_string_to_type(buffer,type);
    }


    return status;

}
/*
* Routine: herr_t mesh_set_dimension(hid_t mid, int dim) 
* Purpose: Set the mesh dimension to dim
* Description: Write the mesh dimension as an integer  attribute equal
*              to dim.
*
* Parameters:
*           mid         IN    HDF5 identifier for the mesh group
*           dim         IN    Mesh dimension
*                              
* Returns: A neagtive value is returned if an error occurs, otherwise returns a zero
*     
* Errors: Possible error conditions are invalid mid input or invalid dimension value.
*       
*/
herr_t mesh_set_dimension(hid_t mid, int dim)
{
    herr_t status = DANU_FAILURE;

    char attr_label[] = MESH_DIM_ATTR_NAME;

    if ( dim < 1 || dim > 4 ) {
        DANU_ERROR_MESS("Invalid dimension value for a mesh");
        return status;
    }

    status = danu_attr_write_int(mid,attr_label,dim);

    return status;
}
/*
* Routine: herr_t mesh_get_dimension(hid_t mid, int * dim) 
* Purpose: Get the mesh dimension by reading an attribute of mid 
* Description: Read the integer attribute MESH_DIM_ATTR_NAME of HDF5 group mid.
* 
*
* Parameters:
*           mid                IN     HDF5 identifier for the mesh group
*           dim                OUT    Mesh attribute that defines the dimension 
*                              
* Returns: A neagtive value is returned if an error occurs, otherwise returns a zero
*     
* Errors: Possible error conditions are invalid mid input or fail to read the attribute.
*       
*/
herr_t mesh_get_dimension(hid_t mid, int * dim)
{
    herr_t status = DANU_FAILURE;

    char attr_label[] = MESH_DIM_ATTR_NAME;

    /* Check input */
    if ( DANU_BAD_PTR(dim) ) {
        DANU_ERROR_MESS("Invalid dimension pointer");
        return status;
    }

    status = danu_attr_read_int(mid,attr_label,dim);

    return status;

}
/*
* Routine: herr_t mesh_set_nnodes(hid_t mid, int nnodes) 
* Purpose: Set the mesh number of nodes attribute
* Description: Write the number of nodes for a mesh, mid, as an integer
*              attribute equal to nnodes.
*
* Parameters:
*           mid         IN    HDF5 identifier for the mesh group
*           nnodes      IN    Number of nodes for mesh 
*                              
* Returns: A neagtive value is returned if an error occurs, otherwise returns a zero
*     
* Errors: Possible error conditions are invalid mid input or invalid nnodes value.
*          
*/
herr_t mesh_set_nnodes(hid_t mid, int nnodes)
{
    herr_t status = DANU_FAILURE;

    char attr_label[] = MESH_NNODES_ATTR_NAME;

    if ( nnodes <= 0  ) {
        DANU_ERROR_MESS("Invalid number of nodes value for a mesh");
        return status;
    }

    status = danu_attr_write_int(mid,attr_label,nnodes);

    return status;
}
/*
* Routine: herr_t mesh_get_nnodes(hid_t mid, int * nnodes) 
* Purpose: Get the number of nodes by reading an attribute of mid 
* Description: Read the integer attribute MESH_NNODES_ATTR_NAME of HDF5 group mid.
* 
*
* Parameters:
*           mid      IN     HDF5 identifier for the mesh group
*           nnodes   OUT    Mesh attribute that defines the  number of nodes
*                              
* Returns: A neagtive value is returned if an error occurs, otherwise returns a zero
*     
* Errors: Possible error conditions are invalid mid input or fail to read the attribute.
*       
*/
herr_t mesh_get_nnodes(hid_t mid, int * nnodes)
{
    herr_t status = DANU_FAILURE;

    char attr_label[] = MESH_NNODES_ATTR_NAME;

    /* Check input */
    if ( DANU_BAD_PTR(nnodes) ) {
        DANU_ERROR_MESS("Invalid number of nodes pointer");
        return status;
    }

    status = danu_attr_read_int(mid,attr_label,nnodes);

    return status;

}
/*
* Routine: herr_t mesh_set_nelem(hid_t mid, int nelem) 
* Purpose: Set the mesh number of elements attribute
* Description: Write the number of elements for a mesh, mid, as an integer
*              attribute equal to nelem. This routine does NOT preform
*              a consistency check, that is is ELEM_ORDER*NELEM = NNODES.
*              It simply writes the attribute value passed in.
*
* Parameters:
*           mid         IN    HDF5 identifier for the mesh group
*           nelem       IN    Number of elements for mesh mid
*                              
* Returns: A neagtive value is returned if an error occurs, otherwise returns a zero
*     
* Errors: Possible error conditions are invalid mid input or non-positive nelem value.
*          
*/
herr_t mesh_set_nelem(hid_t mid, int nelem)
{
    herr_t status = DANU_FAILURE;

    char attr_label[] = MESH_NELEM_ATTR_NAME;

    if ( nelem <= 0  ) {
        DANU_ERROR_MESS("Invalid number of elements value for a mesh");
        return status;
    }

    status = danu_attr_write_int(mid,attr_label,nelem);

    return status;
}
/*
* Routine: herr_t mesh_get_nelem(hid_t mid, int * nelem) 
* Purpose: Get the number of elements by reading an attribute of mid 
* Description: Read the integer attribute MESH_NELEM_ATTR_NAME of HDF5 group mid.
* 
*
* Parameters:
*           mid      IN     HDF5 identifier for the mesh group
*           nelem    OUT    Mesh attribute that defines the  number of elements
*                              
* Returns: A neagtive value is returned if an error occurs, otherwise returns a zero
*     
* Errors: Possible error conditions are invalid mid input or fail to read the attribute.
*       
*/
herr_t mesh_get_nelem(hid_t mid, int * nelem)
{
    herr_t status = DANU_FAILURE;

    char attr_label[] = MESH_NELEM_ATTR_NAME;

    /* Check input */
    if ( DANU_BAD_PTR(nelem) ) {
        DANU_ERROR_MESS("Invalid number of nelem pointer");
        return status;
    }
    *nelem = -1;

    status = danu_attr_read_int(mid,attr_label,nelem);

    return status;

}
/*
* Routine: herr_t mesh_set_elem_order(hid_t mid, int order) 
* Purpose: Set the element order attribute for mesh mid.
* Description: Write the element order (number of nodes per element) for a mesh,
*              mid, as an integer attribute equal to order.
*              This routine does NOT preform
*              a consistency check, that is is ELEM_ORDER matches ELEMENT TYPE.
*              It simply writes the attribute value passed in.
*
* Parameters:
*           mid         IN    HDF5 identifier for the mesh group
*           order       IN    Element order, number of nodes per element
*                              
* Returns: A neagtive value is returned if an error occurs, otherwise returns a zero
*     
* Errors: Possible error conditions are invalid mid input or non-positive order value.
*          
*/
herr_t mesh_set_elem_order(hid_t mid, int order)
{
    herr_t status = DANU_FAILURE;

    char attr_label[] = MESH_EORDER_ATTR_NAME;

    if ( order <= 0  ) {
        DANU_ERROR_MESS("Invalid element order value for a mesh");
        return status;
    }

    status = danu_attr_write_int(mid,attr_label,order);

    return status;
}
/*
* Routine: herr_t mesh_get_elem_order(hid_t mid, int * order) 
* Purpose: Get the element order by reading an attribute of mid 
* Description: Read the integer attribute MESH_NELEM_ATTR_NAME of HDF5 group mid.
* 
*
* Parameters:
*           mid      IN     HDF5 identifier for the mesh group
*           order    OUT    Mesh attribute that defines the element order
*                              
* Returns: A neagtive value is returned if an error occurs, otherwise returns a zero
*     
* Errors: Possible error conditions are invalid mid input or fail to read the attribute.
*       
*/
herr_t mesh_get_elem_order(hid_t mid, int * order)
{
    herr_t status = DANU_FAILURE;

    char attr_label[] = MESH_EORDER_ATTR_NAME;

    /* Check input */
    if ( DANU_BAD_PTR(order) ) {
        DANU_ERROR_MESS("Invalid number of nelem pointer");
        return status;
    }
    *order = -1;

    status = danu_attr_read_int(mid,attr_label,order);

    return status;

}
/*
* Routine: hid_t mesh_create_hex_unstruct(fid,mesh_name,nnodes,x,y,z,nelem,conn) 
* Purpose: Create a mesh group with coordinates (x,y,z) and connectivity data conn.
* Description: A wrapper function to create an unstructured HEX mesh with one function 
*              call. BETA TESTING.
*
* Parameters:
*           fid       IN     HDF5 identifier for the output file
*           mesh_name IN     Mesh name
*           nnodes    IN     Number of nodes (vertices)
*           x,y,z     IN     Pointers to the x,y, and z coordinates of the mesh
*           nelem     IN     Number of elements
*           conn      IN     Pointer to the connectivity array
*                              
* Returns: A neagtive value is returned if an error occurs, otherwise returns a
*          valid HDf5 group identifier.
*     
* Errors: Possible error conditions are invalid input or fail to write the data.
*       
*/
hid_t mesh_create_hex_unstruct(hid_t fid, 
                               const char * mesh_name,
                               int nnodes,
                               double *x,
                               double *y,
                               double *z,
                               int nelem,
                               int *conn)
{
  hid_t mid = H5I_INVALID_HID;
  herr_t status;
  int exists;

  /* Check input */
  if ( H5_ISA_INVALID_ID(fid) ) {
    DANU_ERROR_MESS("Invalid file identifier");
    return mid;
  }

  if ( DANU_BAD_STRING(mesh_name) ) {
    DANU_ERROR_MESS("Invalid mesh name.");
    return mid;
  }

  mesh_exists(fid,mesh_name,&exists);
  if ( exists ) {
    DANU_ERROR_MESS("Mesh group exists. Will not over write.");
    return mid;
  }

  if ( nnodes <= 0 ) {
    DANU_ERROR_MESS("Invalid number of nodes");
    return mid;
  }

  if ( DANU_BAD_PTR(x) || DANU_BAD_PTR(y) || DANU_BAD_PTR(z) ) {
    DANU_ERROR_MESS("Invalid coordinate pointer");
    return mid;
  }

  if ( nelem <= 0 ) {
    DANU_ERROR_MESS("Invalid number of elements");
    return mid;
  }

  if ( DANU_BAD_PTR(conn) ) {
    DANU_ERROR_MESS("Invalid connectivity pointer");
    return mid;
  }


  /* Steps to create a hex mesh */
  status=mesh_create(fid,mesh_name,UNSTRUCTURED_MESH,HEX_ELEM,&mid);
  if ( DANU_RETURN_OK(status) ) {
    status &= mesh_write_coordinates(mid,nnodes,x,y,z);
    status &= mesh_write_connectivity(mid,nelem,conn);

    if ( DANU_RETURN_FAIL(status) ) {
      DANU_ERROR_MESS("Write coordinates or connectivity failed.");
      danu_group_close(mid);
      mid=H5I_INVALID_HID;
    }
  }
  else {
    DANU_ERROR_MESS("Failed to create HEX mesh");
  }

  return mid;
}




