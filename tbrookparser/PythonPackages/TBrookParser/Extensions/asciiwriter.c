// asciiwriter.c

#include <Python.h>
#include <numpy/arrayobject.h>
#include <stdio.h>

static int printme(FILE *fp,
		   char type,
		   int irank,
		   int *shape,
		   int iCurrent,
		   void *data,
		   char *format,
		   char *hdr,
		   int nperline,
		   int debug)
{
  /**
   * Recursive array printing function.  Note that currently only the following
   * types are recognized:
   *   'i' = integer
   *   'f' = float
   *   'd' = double
   *
   * Returns the number of bytes written.
   * Returns -1 on error.
   **/

  int offset = 0;

  /**
   * Catch a type error immediately **/
  switch (type) {
  case 'i':
  case 'f':
  case 'd':
    break;
  default:
    return(-1);
  }
  
  if(debug) printf("In printme with:type=%c rank=%d, iC=%d, shape(ic)=%d, nperline=%d\n", type, irank, iCurrent, shape[iCurrent], nperline);

  if (irank == 0 ) {
    /**
     * Rank 0 gets printed.
     * We should be really never be here
     **/
    fprintf(fp,"%s",hdr);
    fprintf(fp,format,data);
  }
  else if (iCurrent == irank-1) {
    /**
     * Rank 1
     * This is where the bulk of the writing gets done
     **/
    int i, j;
    i = shape[iCurrent];
    fprintf(fp,"%s",hdr);
    switch (type) {
    case 'i':
      {
	int *pi;
	pi = (int *) data;
	for(j=0; j<i; j++) {
	  if(j && (j%nperline == 0)) fprintf(fp,"\n");
	  fprintf(fp,format,pi[j]);
	}
	if((j%nperline==0)) fprintf(fp,"\n");
	offset = offset + i*sizeof(int);
      }
      break;
    case 'f':
      {
	float *pf;
	pf = (float *) data;
	for(j=0; j<i; j++) {
	  if(j && (j%nperline == 0)) fprintf(fp,"\n");
	  fprintf(fp,format,pf[j]);
	}
        if((j%nperline == 0)) fprintf(fp,"\n");
	offset = offset + i*sizeof(float);
      }
      break;
    case 'd':
      {
	double *pd = (double *) data;
	for(j=0; j<i; j++) {
	  if((j && j%nperline == 0)) fprintf(fp,"\n");
	  fprintf(fp,format,pd[j]);
	}
	if(j%nperline == 0) fprintf(fp,"\n");
	offset = offset + i*sizeof(double);
      }
      break;
    default:
      printf("Unknown operator %c\n",type);
      return(-1);
    }
  } else {
    /**
     * Recursively call self till we reach rank 1 **/
    int i, j;
    char *d;
    i = shape[iCurrent];
    d = (char *) data;
    for(j=0; j<i; j++) {
      int ioff;
      d = (char *) data + offset;
      ioff = printme(fp, type, irank, shape, iCurrent+1, (void *)d, format, hdr, nperline, debug);
      if(ioff<0)
	return(-1);
      else
	offset = offset + ioff;
    }
  }
  return(offset);
}

static int printquery(FILE *fp,
		      int rank,
		      int *shape,
		      int face,
		      void *posns,
		      void *indices,
		      int indlen,
		      void *data,
		      char *cformat,
		      char *dformat)
{
  int h;
  int i;
  int j;
  int k;
  int kk;
  int offset;
  int *pind;
  double *pd;
  double *pps;
  double xpos;
  double ypos;
  double zpos;
  double thisdpt;

  pd      = (double *) data; 
  pind    = (int *) indices;
  pps     = (double *) posns;
  offset  = 3;

  if ( rank == 1 ) {
    /**
     * rank 1 arrays
     **/

    i = shape[0];
    for (j=0; j<=indlen-1; j++) {
      k    = pind[j]-1;
      xpos = *(pps + k*offset);
      ypos = *(pps + k*offset + 1);
      zpos = *(pps + k*offset + 2);
      fprintf(fp,cformat,k+1,xpos,ypos,zpos);
      fprintf(fp,dformat,pd[k]);
      fprintf(fp,"   \n");
    }
  } else if ( rank == 2) {
    i = shape[1];
    for (j=0; j<=indlen-1; j++) {
      kk   = pind[j]-1;
      xpos = *(pps + kk*offset);
      ypos = *(pps + kk*offset + 1);
      zpos = *(pps + kk*offset + 2);
      fprintf(fp,cformat,kk+1,xpos,ypos,zpos);
      for (k=0; k<i; k++) {
	thisdpt = *(pd + kk*i + k);
	fprintf(fp,dformat,thisdpt);
      }
      fprintf(fp,"   \n");
    }
  } else if ( rank == 3) {
    h = shape[1];
    i = shape[2];
    for (j=0; j<=indlen-1; j++) {
      kk   = pind[j]-1;
      xpos = *(pps + kk*offset);
      ypos = *(pps + kk*offset + 1);
      zpos = *(pps + kk*offset + 2);
      fprintf(fp,cformat,kk+1,xpos,ypos,zpos);
      for (k=0; k<i; k++) {
	thisdpt = *(pd + (kk*h + (face-1))*i + k);
	fprintf(fp,dformat,thisdpt);
      }
      fprintf(fp,"   \n");
    }
  }
  return(offset);
}

static int printreducedfield(FILE *fp,
			     int rank,
			     int *shape,
			     int face,
			     void *posns,
			     void *indices,
			     int indlen,
			     void *data,
			     char *cformat,
			     char *dformat)
{
  int h;
  int i;
  int j;
  int k;
  int kk;
  int offset;
  int *pind;
  double *pd;
  double *pps;
  double xpos;
  double ypos;
  double zpos;
  double thisdpt;

  pd      = (double *) data; 
  pind    = (int *) indices;
  pps     = (double *) posns;
  offset  = 3;

  if ( rank == 1 ) {
    /**
     * rank 1 arrays
     **/
    i = shape[0];
    for (j=0; j<=indlen-1; j++) {
      k    = pind[j]-1;
      xpos = *(pps + k*offset);
      ypos = *(pps + k*offset + 1);
      zpos = *(pps + k*offset + 2);
      fprintf(fp,cformat,k+1,xpos,ypos,zpos);
      fprintf(fp,dformat,pd[j]);
      fprintf(fp,"   \n");
    }
  } else if ( rank == 2) {
    i = shape[1];
    for (j=0; j<=indlen-1; j++) {
      kk   = pind[j]-1;
      xpos = *(pps + kk*offset);
      ypos = *(pps + kk*offset + 1);
      zpos = *(pps + kk*offset + 2);
      fprintf(fp,cformat,kk+1,xpos,ypos,zpos);
      for (k=0; k<i; k++) {
	thisdpt = *(pd + j*i + k);
	fprintf(fp,dformat,thisdpt);
      }
      fprintf(fp,"   \n");
    }
  } else if ( rank == 3) {
    h = shape[1];
    i = shape[2];
    for (j=0; j<=indlen-1; j++) {
      kk   = pind[j]-1;
      xpos = *(pps + kk*offset);
      ypos = *(pps + kk*offset + 1);
      zpos = *(pps + kk*offset + 2);
      fprintf(fp,cformat,kk+1,xpos,ypos,zpos);
      for (k=0; k<i; k++) {
	thisdpt = *(pd + (j*h + (face-1))*i + k);
	fprintf(fp,dformat,thisdpt);
      }
      fprintf(fp,"   \n");
    }
  }
  return(offset);
}

static int printprobe(FILE *fp,
		      int rank,
		      int *shape,
		      void *indices,
		      int indlen,
		      void *data,
		      char *iformat,
		      char *dformat)
{
  int i;
  int j;
  int k;
  int kk;
  int *pind;
  double *pd;
  double thisdpt;
  int result = 0;

  pd      = (double *) data; 
  pind    = (int *) indices;

  if ( rank == 2) {
    i = shape[1];
    for (j=0; j<=indlen-1; j++) {
      kk     = j;
      for (k=0; k<i; k++) {
	thisdpt = *(pd + kk*i + k);
	fprintf(fp,dformat,thisdpt);
      }
      fprintf(fp,iformat,pind[j]);
      fprintf(fp,"\n");
    }
  } else {
    printf("Rank != 2 arrays not valid in probe mode");
    result = -1;
  }
  return(result);
}

static PyObject *writeQuery(PyObject *self, PyObject *args)
{
  int         rank;
  int         indlen;
  int         face;
  char       *hdr;
  PyObject   *file;
  int         debug;
  char       *coordsfmt;
  char       *fieldfmt;
  PyArrayObject *pshape;
  PyArrayObject *pdata;
  PyArrayObject *pindices;
  PyArrayObject *pposns;

  if (!PyArg_ParseTuple(args,"iOiOOiOsOssi:writeQuery",
			&rank,
			&pshape,
			&face,
			&pdata,
			&pindices,
			&indlen,
			&pposns,
			&hdr,
			&file,
			&coordsfmt,
			&fieldfmt,
			&debug))
    return NULL;

  fprintf(PyFile_AsFile(file), "%s", hdr);  

  if (rank <= 1) {
    if (indlen == pshape->data[0]) {
      if (printreducedfield(PyFile_AsFile(file), rank, (int *)pshape->data, 
	face, pposns->data, pindices->data, indlen, pdata->data, coordsfmt, 
	fieldfmt) < 0)
	return NULL;
    } else {
      if (printquery(PyFile_AsFile(file), rank, (int *)pshape->data, face, 
	pposns->data, pindices->data, indlen, pdata->data, coordsfmt, fieldfmt)
	< 0)
	return NULL;
    }
  } else if (rank == 2) {
    if (indlen == pshape->data[0]) {
      if (printreducedfield(PyFile_AsFile(file), rank, (int *)pshape->data, 
	face, pposns->data, pindices->data, indlen, pdata->data, coordsfmt, 
	fieldfmt) < 0)
	return NULL;
    } else 
      if (printquery(PyFile_AsFile(file), rank, (int *)pshape->data, face, 
	pposns->data, pindices->data, indlen, pdata->data, coordsfmt, fieldfmt)
	 < 0)
	return NULL;
  } else if (rank == 3) {
    if (indlen == pshape->data[0]) {
      if (printreducedfield(PyFile_AsFile(file), rank, (int *)pshape->data, 
	face, pposns->data, pindices->data, indlen, pdata->data, coordsfmt, 
	fieldfmt) < 0)
	return NULL;
    } else {
      if (printquery(PyFile_AsFile(file), rank, (int *)pshape->data, face, 
	pposns->data, pindices->data, indlen, pdata->data, coordsfmt, fieldfmt) 
	< 0)
	return NULL;
    }
  } else {
    fprintf(stderr, "Rank > 3 arrays not valid in query mode");
    return NULL;
  }

  Py_RETURN_NONE;
}

static PyObject *writeProbe(PyObject *self, PyObject *args)
{
  int         rank;
  int         cyclen;
  char       *hdr;
  PyObject   *file;
  int         debug;
  char       *cyclefmt;
  char       *probefmt;
  PyArrayObject *pshape;
  PyArrayObject *pdata;
  PyArrayObject *pcycles;

  if (!PyArg_ParseTuple(args,"iOOOisOssi:writeProbe",
			&rank,
			&pshape,
			&pdata,
			&pcycles,
			&cyclen,
			&hdr,
			&file,
			&cyclefmt,
			&probefmt,
			&debug))
    return NULL;

  fprintf(PyFile_AsFile(file), "%s", hdr);  

  if (rank == 2) {
    if (printprobe(PyFile_AsFile(file), rank, (int *)pshape->data, pcycles->data, cyclen, pdata->data, cyclefmt, probefmt) < 0)
      return NULL;
  } else {
    fprintf(stderr, "Rank !=2 arrays not valid in probe mode: rank=%d", rank);
    return NULL;
  }

  Py_RETURN_NONE;
}

static PyObject *writeArrayASCII(PyObject *self, PyObject *args)
{
  char       *type;
  int         rank;
  char       *hdr;
  char       *format;
  int         nperline;
  PyObject   *file;
  int         debug;
  FILE       *fp;
  PyArrayObject *pshape;
  PyArrayObject *pdata;

  if (!PyArg_ParseTuple(args,"siOOssiOi:writeArrayASCII",
			&type,
			&rank,
			&pshape,
			&pdata,
			&hdr,
			&format,
			&nperline,
			&file,
			&debug))
    return NULL;

  if (printme(PyFile_AsFile(file), *type, rank, (int *)pshape->data, 0, pdata->data, format, hdr, nperline, debug) < 0)
    return NULL;

  Py_RETURN_NONE;
}

static PyMethodDef asciiwriter_methods[] = {
  { "writeQuery",      writeQuery,      METH_VARARGS },
  { "writeProbe",      writeProbe,      METH_VARARGS },
  { "writeArrayASCII", writeArrayASCII, METH_VARARGS },
  { NULL, NULL }
};

void initasciiwriter(void)
{
  Py_InitModule("asciiwriter", asciiwriter_methods);
}
