// ensight.c

#include <Python.h>
#include <numpy/arrayobject.h>

static PyObject *getconnectivity(PyObject *self, PyObject *args)
{
  PyArrayObject *pmeshconnectivity;
  PyArrayObject *pconnectivity;
  PyArrayObject *pcellids;
  PyArrayObject *pvtxids;
  PyArrayObject *pivtx;
  int ncellsmesh;
  int nnodesmesh;
  int ncellsb;
  int nnodesb;
  int nvc;
  int i;
  int j;
  int k;
  int *meshconnectivity;
  int *connectivity;
  int *cellids;
  int *vtxids;
  int *ivtx;

  if (!PyArg_ParseTuple(args, "OOOOOiiiii:getconnectivity",
			&pmeshconnectivity,
			&pcellids,
			&pvtxids,
			&pconnectivity,
			&pivtx,
			&ncellsmesh,
			&ncellsb,
			&nnodesmesh,
			&nnodesb,
			&nvc))
    return NULL;

  meshconnectivity = (int *)pmeshconnectivity->data;
  connectivity     = (int *)pconnectivity->data;
  cellids          = (int *)pcellids->data;
  vtxids           = (int *)pvtxids->data;
  ivtx             = (int *)pivtx->data;

  for ( i = 0; i < nnodesb; i++ )
    ivtx[vtxids[i]-1] = i + 1;

  for ( i = 0; i < ncellsb; i++ )
    for ( j = 0; j < nvc; j++) {
      k = meshconnectivity[(cellids[i]-1) * nvc + j];
      connectivity[i * nvc + j] = ivtx[k-1];
    }

  Py_RETURN_NONE;
}

static PyMethodDef ensight_methods[] = {
  { "getconnectivity", getconnectivity, METH_VARARGS },
  { NULL, NULL }
};

void initensight(void)
{
  Py_InitModule("ensight", ensight_methods);
}
