/* ************************************************************************* *
*  Genernal Mesh Viewer                                                      *
*  ************************************************************************* */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "GMVIface.h"

#define GMV_MALLOC(type,num)    (type *) malloc(num*sizeof(type))  
#define GMV_FREE(ptr)           free(ptr) 

extern void gmvwrite_openfile_ir_ascii(const char *, int, int );
extern void gmvwrite_closefile();


/* Global static constant that prevents clashing in the gmvwrite library */
typedef unsigned long long pointer_t;
static int gmvwrite_active = 0;
static pointer_t gmvwrite_pointer = 0;

MeshFile * allocate_MeshFile(const char *filename, const char *mesh_type) {
   MeshFile *mesh=NULL;

   mesh=GMV_MALLOC(MeshFile,1);
   mesh->name=GMV_MALLOC(char,strlen(filename)+1);
   sprintf(mesh->name,"%s",filename);
   if ( mesh_type != NULL ) {
     mesh->type=GMV_MALLOC(char,strlen(mesh_type)+1);
     sprintf(mesh->type,"%s",mesh_type);
   }
   else {
     mesh->type=GMV_MALLOC(char,4);
     sprintf(mesh->type,"%s","hex");
   }

   return mesh;
}

void deallocate_MeshFile(MeshFile *file) {

  if ( !file->name) GMV_FREE(file->name);
  if ( !file->type ) GMV_FREE(file->type);

  GMV_FREE(file);

}

MeshFile * copy_MeshFile(const MeshFile *src) {
  
  MeshFile *dest;

  dest = allocate_MeshFile(src->name,src->type);
  dest->nnodes=src->nnodes;
  dest->ncells=src->ncells;
  dest->isize=src->isize;
  dest->rsize=src->rsize;

  return dest;

}
  
DataFile * allocate_DataFile(const char * filename, const MeshFile *mesh) {

  DataFile *file;

  file=GMV_MALLOC(DataFile,1);
  file->mesh=copy_MeshFile(mesh);

  file->name=GMV_MALLOC(char, strlen(filename)+1);
  sprintf(file->name,"%s",filename);


  return file;

}

void deallocate_DataFile(DataFile *file) {

  if ( !file->name ) GMV_FREE(file->name);
  if ( !file->mesh ) deallocate_MeshFile(file->mesh);

  GMV_FREE(file);

}

DataFile * copy_DataFile(DataFile *src) {
  DataFile *dest;

  dest=allocate_DataFile(src->name,src->mesh);
  dest->nnodes=src->nnodes;
  dest->ncells=src->ncells;
  dest->isize=src->isize;
  dest->rsize=src->rsize;

  return dest;

}
int mesh_type_elem_order(const char *type) {

   if ( strcmp(type,"hex") == 0 ) {
     return 8;
   }
   else if ( strcmp(type,"phex8") == 0 ) {
     return 8;
   }
   else if ( strcmp(type,"quad") == 0 ) {
     return 4;
   }
   else {
    printf("(%s: %d) Error: Invalid mesh type. Can not determine number of vertices",
	     __FILE__, __LINE__);
    return -1;
  }
}

int is_gmvwrite_active() {
  return gmvwrite_active;
}

int is_gmvobject_valid(void * ptr) {
  int status = 0;
  pointer_t cpy_ptr = (pointer_t) ptr;
  if ( is_gmvwrite_active() && (cpy_ptr == gmvwrite_pointer) ) {
    status = 1;
  }
  return status;
}

void open_gmvwrite_file(void *file_ptr, const char * filename, int binary, int isize, int rsize) {
  if (binary) {
    gmvwrite_openfile_ir(filename,isize,rsize);
  }
  else {
    gmvwrite_openfile_ir_ascii(filename,isize,rsize);
  }
  gmvwrite_active=1;
  gmvwrite_pointer=(pointer_t) file_ptr;
}

void close_gmvwrite_file(void *file_ptr) {
  pointer_t ptr = (pointer_t) file_ptr;

  if ( ptr == gmvwrite_pointer ) {
    gmvwrite_closefile();
    gmvwrite_pointer = 0;
    gmvwrite_active=0;
  }
  else {
    printf("(%s %d): Attempt to close file that is not active\n",
	   __FILE__,__LINE__);
  }
}


