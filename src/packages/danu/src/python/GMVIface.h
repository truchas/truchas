/* SWIG processed header file */
#ifndef _GMV_OBJECT_H_
#define _GMV_OBJECT_H_

#define GMV_DATA_CELL 0
#define GMV_DATA_NODE 1
#define GMV_DATA_FACE 2


int is_gmvwrite_active();


typedef struct {

  char * name;
  char * type;
  int nnodes;
  int ncells;
  int isize;
  int rsize;

} MeshFile;

typedef struct {

  char * name;
  MeshFile * mesh;
  int nnodes;
  int ncells;
  int cycle;
  double time;
  int isize;
  int rsize;

} DataFile;


#endif
