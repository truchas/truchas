/* Fortran-compatible interface routines to gmvwrite.c.
   Function names have a single appended underscore that many (but not all)
   Fortran compilers will expect for external functions.  Functions that are
   passed character arrays include a hidden pass-by-value integer argument
   appended to the end of the argument list.  Most Fortran compilers pass the
   length of the actual character string in this argument; again, this is not
   necessarily true for all compilers, but I've yet to find one that doesn't.
*/
   
#include <string.h>
#include <stdlib.h>

#include <fgmvwrite.h>

char tmpname[9];  /* 8-character string + terminating null */

void make_tmpname (char *string, int length)
{
  int i;
  for (i=0; (i<length) && (i<8); i++)
    tmpname[i] = string[i];
  for (i=length; i<8; i++)
    tmpname[i] = ' ';
  tmpname[8] = '\0';
}

void gmvwrite_openfile_ir_f (char filenam[], int *isize, int *rsize, int len)
{
  char *tmpbuf;
  tmpbuf = (char *)malloc((1+len)*sizeof(char));
  strncpy(tmpbuf,filenam,len);
  tmpbuf[len] = '\0';
  gmvwrite_openfile_ir(tmpbuf,*isize,*rsize);
  free(tmpbuf);
}

void gmvwrite_openfile_ir_ascii_f (char filenam[], int *isize, int *rsize, int len)
{
  char *tmpbuf;
  tmpbuf = (char *)malloc((1+len)*sizeof(char));
  strncpy(tmpbuf,filenam,len);
  tmpbuf[len] = '\0';
  gmvwrite_openfile_ir_ascii(tmpbuf,*isize,*rsize);
  free(tmpbuf);
}

void gmvwrite_closefile_f (void)
{
  gmvwrite_closefile ();
}

void gmvwrite_probtime_f (double *ptime)
{
  gmvwrite_probtime (*ptime);
}

void gmvwrite_cycleno_f (int *cyclenum)
{
  gmvwrite_cycleno (*cyclenum);
}

void gmvwrite_node_data_f (void *nndes, void *x, void *y, void *z)
{
  gmvwrite_node_data (nndes, x, y, z);
}

void gmvwrite_nodeids_f (void *nodeids)
{
  gmvwrite_nodeids (nodeids);
}

void gmvwrite_cell_header_f (void *ncells)
{
  gmvwrite_cell_header(ncells);
}

void gmvwrite_cell_type_f (char cell_type[], int *nverts, void *nodes, int len)
{
  make_tmpname(cell_type, len);
  gmvwrite_cell_type (tmpname, *nverts, nodes);
}

void gmvwrite_cellids_f (void *cellids)
{
  gmvwrite_cellids (cellids);
}

void gmvwrite_material_header_f (int *nmats, int *data_type)
{
  gmvwrite_material_header (*nmats, *data_type);
}

void gmvwrite_material_name_f (char matname[], int len)
{
  make_tmpname(matname, len);
  gmvwrite_material_name(tmpname);
}

void gmvwrite_material_ids_f (int *matids, int *data_type)
{
  gmvwrite_material_ids(matids, *data_type);
}

void gmvwrite_variable_header_f (void)
{
  gmvwrite_variable_header();
}

void gmvwrite_variable_name_data_f (int *data_type, char varname[], void *vids, int len)
{
  make_tmpname(varname, len);
  gmvwrite_variable_name_data(*data_type, tmpname, vids);
}

void gmvwrite_variable_endvars_f (void)
{
  gmvwrite_variable_endvars();
}

void gmvwrite_flag_header_f (void)
{
  gmvwrite_flag_header();
}

void gmvwrite_flag_name_f (char flagname[], int *numtypes, int *data_type, int len)
{
  make_tmpname(flagname, len);
  gmvwrite_flag_name(tmpname, *numtypes, *data_type);
}

void gmvwrite_flag_subname_f (char subname[], int len)
{
  make_tmpname(subname, len);
  gmvwrite_flag_subname(tmpname);
}

void gmvwrite_flag_data_f (int *data_type, int flag_data[])
{
  gmvwrite_flag_data(*data_type, flag_data);
}

void gmvwrite_flag_endflag_f (void)
{
  gmvwrite_flag_endflag();
}

