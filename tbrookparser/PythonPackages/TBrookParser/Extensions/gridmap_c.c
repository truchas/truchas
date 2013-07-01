// gridmap.c

#include <Python.h>
#include <numpy/arrayobject.h>
 
static PyObject *getfield(PyObject *self, PyObject *args)
{
  char *map_rule;
  int max_warnings;
  double default_value;
  int n_dimensions;
  int n_cells_a;
  int n_vertices_a;
  int n_nodes_a;
  PyArrayObject *data_a;
  PyArrayObject *connectivity_a;
  PyArrayObject *coordinates_a;
  PyArrayObject *cblocks_a;
  int n_cells_b;
  int n_vertices_b;
  int n_nodes_b;
  PyArrayObject *data_b;
  PyArrayObject *connectivity_b;
  PyArrayObject *coordinates_b;
  PyArrayObject *cblocks_b;
  double integral_a;
  double integral_b;
  int ierr;

  if (!PyArg_ParseTuple(args, "sidiiiiOOOOiiiOOOO:getfield",
			&map_rule,
			&max_warnings,
			&default_value,
			&n_dimensions,
			&n_cells_a,
			&n_vertices_a,
			&n_nodes_a,
			&data_a,
			&connectivity_a,
			&coordinates_a,
			&cblocks_a,
			&n_cells_b,
			&n_vertices_b,
			&n_nodes_b,
			&data_b,
			&connectivity_b,
			&coordinates_b,
			&cblocks_b))
    return NULL;

  mapfield_(&n_cells_a,
	    &n_vertices_a,
	    connectivity_a->data,
	    &n_nodes_a,
	    coordinates_a->data,
	    cblocks_a->data,
	    &n_cells_b,
	    &n_vertices_b,
	    connectivity_b->data,
	    &n_nodes_b,
	    coordinates_b->data,
	    cblocks_b->data,
	    &n_dimensions,
	    &max_warnings,
	    map_rule,
	    data_a->data,
	    data_b->data,
	    &integral_a,
	    &integral_b,
	    &default_value,
	    &ierr,
	    strlen(map_rule));
  
  return Py_BuildValue("ddi", integral_a, integral_b, ierr);
}

static PyMethodDef gridmap_methods[] = {
  { "getfield", getfield, METH_VARARGS },
  { NULL, NULL }
};


void initgridmap(void)
{
  Py_InitModule("gridmap", gridmap_methods);
}
