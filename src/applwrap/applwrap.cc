/*******************************************
 * A Simple APPLgrid python wrapper
 * S. Carrazza - April 2015
 ********************************************/

#include <Python.h>
#include <LHAPDF/LHAPDF.h>
#include <LHAPDF/Exceptions.h>
#include <appl_grid/appl_grid.h>
#include <appl_grid/appl_igrid.h>
using std::vector;
using std::cout;
using std::endl;
using std::exception;

// I hate singletons - sc
appl::grid *_g = nullptr;
vector<LHAPDF::PDF*> _pdfs;
int _imem = 0;

extern "C" void evolvepdf_(const double& x,const double& Q, double* pdf)
{  
  for (int i = 0; i < 13; i++) 
    {
      const int id = i-6;
      pdf[i] = _pdfs[_imem]->xfxQ(id, x, Q);
    }
}

extern "C" double alphaspdf_(const double& Q)
{
  return _pdfs[_imem]->alphasQ(Q);
}

static PyObject* py_initpdf(PyObject* self, PyObject* args)
{
  char* setname;
  PyArg_ParseTuple(args, "s", &setname);
  
  for (int i = 0; i < (int) _pdfs.size(); i++)
    if (_pdfs[i]) delete _pdfs[i];
  _pdfs.clear();

  try
  {
    _pdfs = LHAPDF::mkPDFs(setname);
  }
  catch (LHAPDF::Exception e)
  {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  }
  _imem = 0;

  return Py_BuildValue("");
}

static PyObject* py_pdfreplica(PyObject* self, PyObject* args)
{
  int nrep;  
  PyArg_ParseTuple(args, "i", &nrep);
  _imem = nrep;

  return Py_BuildValue("");
}

static PyObject* py_initobs(PyObject* self, PyObject* args)
{  
  char *file;
  PyArg_ParseTuple(args,"s", &file);
    
  if (_g) delete _g;
  try
  {
    _g = new appl::grid(file);
  }
  catch(appl::grid::exception e)
  {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  }

  return Py_BuildValue("");
}

static PyObject* py_convolute(PyObject* self, PyObject* args)
{  
  int pto;
  PyArg_ParseTuple(args,"i", &pto);
    
  if (!_g) exit(-1);  
  vector<double> xsec = _g->vconvolute(evolvepdf_,alphaspdf_,pto);

  PyObject *out = PyList_New(xsec.size());
  for (int i = 0; i < (int) xsec.size(); i++)
    PyList_SET_ITEM(out, i, PyFloat_FromDouble(xsec[i]));

  return out;
}

static PyObject* py_getobsq(PyObject* self, PyObject* args)
{  
  int pto, bin;
  PyArg_ParseTuple(args,"ii", &pto, &bin);

  vector<double> Q;

  int iorder = pto; 
  if (_g->calculation() == 1) // if aMCfast change iorder
    {
      if (pto == 0) iorder = 3;
      if (pto == 1) iorder = 0;
    }

  appl::igrid const *igrid = _g->weightgrid(iorder,bin);  
  for (int ix1 = 0; ix1 < igrid->Ny1(); ix1++)
    for (int ix2 = 0; ix2 < igrid->Ny2(); ix2++)
      for (int t = 0; t < igrid->Ntau(); t++)
	for (int ip = 0; ip < _g->subProcesses(0); ip++)
	  {

	    const bool zero_weight = (*(const SparseMatrix3d*) const_cast<appl::igrid*>(igrid)->weightgrid(ip))(t,ix1,ix2) == 0;
	    
	    if (!zero_weight)
	      Q.push_back(sqrt(igrid->fQ2(igrid->gettau(t))));	    
	  }
  
  double sum = 0;
  for (int i = 0; i < (int) Q.size(); i++) sum += Q[i];
  sum /= Q.size();
  
  return Py_BuildValue("d", sum);
}

static PyObject* py_getnbins(PyObject* self, PyObject* args)
{  
  int nbins;
  try
  {
    nbins = _g->Nobs();
  }
  catch(appl::grid::exception e)
  {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  }

  return Py_BuildValue("i", nbins);
}

static PyMethodDef applwrap_methods[] = {
  {"initpdf", py_initpdf, METH_VARARGS},
  {"pdfreplica", py_pdfreplica, METH_VARARGS},
  {"initobs", py_initobs, METH_VARARGS},
  {"convolute", py_convolute, METH_VARARGS},
  {"getobsq", py_getobsq, METH_VARARGS},
  {"getnbins",py_getnbins, METH_VARARGS},
  {NULL, NULL}
};

extern "C" void initapplwrap()
{
  (void) Py_InitModule("applwrap", applwrap_methods);
}
