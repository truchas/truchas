// exodus.c

#include <Python.h>
#include <numpy/arrayobject.h>
#include <string.h>

static PyObject *GetSizes(PyObject *self, PyObject *args)
{
  char *path;
  int nnodes, ncells;
  
  if (!PyArg_ParseTuple(args, "s:GetSizes", &path)) return NULL;
  
  parser_read_exodus_mesh_size_(path, &nnodes, &ncells, strlen(path));
  
  return Py_BuildValue("ii", nnodes, ncells);
}


static PyObject *GetMesh(PyObject *self, PyObject *args)
{
  char *path;
  PyArrayObject *connect, *coord, *cblock;
  
  if (!PyArg_ParseTuple(args, "sOOO:GetMesh", &path, &connect, &coord, &cblock)) return NULL;
  
  parser_read_exodus_mesh_(path, connect->data, coord->data, cblock->data, strlen(path));
  
  Py_RETURN_NONE;
}

static PyObject *GetCellFields(PyObject *self, PyObject *args)
{
  int n_dimensions;
  int n_cells;
  int n_nodes;
  int n_vertices;
  PyArrayObject *connectivity;
  PyArrayObject *coordinates;
  PyArrayObject *cell_centroids;
  PyArrayObject *cell_n_neighbors;

  if (!PyArg_ParseTuple(args, "iiiiOOOO:GetCellFields",
			&n_cells,
			&n_nodes,
			&n_dimensions,
			&n_vertices,
			&connectivity,
			&coordinates,
			&cell_centroids,
			&cell_n_neighbors))
    return NULL;

  getcellfields_(&n_dimensions,
		 &n_cells,
		 &n_nodes,
		 &n_vertices,
		 connectivity->data,
		 coordinates->data,
		 cell_centroids->data,
		 cell_n_neighbors->data);

  Py_RETURN_NONE;
}

static PyMethodDef exodus_methods[] = {
  { "GetSizes",      GetSizes,      METH_VARARGS },
  { "GetMesh",       GetMesh,       METH_VARARGS },
  { "GetCellFields", GetCellFields, METH_VARARGS },
  { NULL, NULL }
};

void initexodus(void)
{
  Py_InitModule("exodus", exodus_methods);
}
