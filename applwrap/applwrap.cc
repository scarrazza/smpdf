/*******************************************
 * A Simple APPLgrid python wrapper
 * S. Carrazza - April 2015
 ********************************************/

#include <Python.h>
#include <LHAPDF/LHAPDF.h>
#include <appl_grid/appl_grid.h>
using std::vector;
using std::cout;
using std::endl;

extern "C" void evolvepdf_(const double&,const double&,double*);
extern "C" double alphaspdf_(const double&);

static PyObject* py_loadpdf(PyObject* self, PyObject* args)
{
  int nrep;  
  char* setname;
  
  PyArg_ParseTuple(args, "si", &setname, &nrep);
  LHAPDF::initPDFSet(setname, nrep);

  return Py_BuildValue("");
}


static PyObject* py_convolute(PyObject* self, PyObject* args)
{  
  int pto;
  char *file;
  PyArg_ParseTuple(args,"si", &file, &pto);
  
  appl::grid g(file);  
  vector<double> xsec = g.vconvolute(evolvepdf_,alphaspdf_,pto);

  PyObject *out = PyList_New(xsec.size());
  for (int i = 0; i < (int) xsec.size(); i++)
    PyList_SET_ITEM(out, i, PyFloat_FromDouble(xsec[i]));

  return out;
}

static PyMethodDef applwrap_methods[] = {
  {"loadpdf", py_loadpdf, METH_VARARGS},
  {"convolute", py_convolute, METH_VARARGS},
  {NULL, NULL}
};

extern "C" void initapplwrap()
{
  (void) Py_InitModule("applwrap", applwrap_methods);
}
