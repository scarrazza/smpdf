/*******************************************
 * A Simple APPLgrid python wrapper
 * S. Carrazza - April 2015
 ********************************************/

#include "Python.h"
#include <LHAPDF/LHAPDF.h>
#include <LHAPDF/Exceptions.h>
#include <appl_grid/appl_grid.h>
#include <appl_grid/appl_igrid.h>
using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::exception;

// I hate singletons - sc
appl::grid *_g = nullptr;
vector<LHAPDF::PDF*> _pdfs;
int _imem = 0;
static PyObject *_custom_xfxq = NULL;

extern "C" void evolvepdf_(const double& x,const double& Q, double* pdf)
{
  for (int i = 0; i < 13; i++)
  {
      const int id = i-6;
      pdf[i] = _pdfs[_imem]->xfxQ(id, x, Q);
   }
}


extern "C" void evolvepdf_python_(const double& x,const double& Q, double* pdf)
{
  PyObject *py_result;
  PyObject *arglist;
  double result;

  for (int i = 0; i < 13; i++)
  {
    const int id = i-6;
    arglist = Py_BuildValue("(iidd)", _imem, id, x, Q);
    py_result = PyObject_CallObject(_custom_xfxq, arglist);
    if (!py_result){
        PyErr_Print();
        //TODO: Can we raise python exception somewhere instead of exiting?
        exit(1);
    }
    Py_DECREF(arglist);
    result =  PyFloat_AsDouble(py_result);
    Py_DECREF(py_result);
    pdf[i] = result;
  }
}
extern "C" double alphaspdf_(const double& Q)
{
  return _pdfs[0]->alphasQ(Q);
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

static PyObject* py_setverbosity(PyObject *self, PyObject *args)
{
  int ver;
  PyArg_ParseTuple(args, "i", &ver);
  LHAPDF::setVerbosity(ver);

  return Py_BuildValue("");
}

static PyObject* py_xfxQ(PyObject *self, PyObject *args)
{
  int rep, fl;
  double x, Q, res;
  PyArg_ParseTuple(args, "iidd", &rep, &fl, &x, &Q);

  if (_pdfs[rep] == 0)
    {
      PyErr_SetString(PyExc_ValueError, "PDF not allocated");
      return NULL;
    }

  try
  {
    res = _pdfs[rep]->xfxQ(fl, x, Q);
  }
  catch (LHAPDF::Exception e)
  {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  }

  return Py_BuildValue("d", res);
}

static PyObject* py_q2Min(PyObject *self, PyObject *args)
{
  double res;

  if (_pdfs[0] == 0)
    {
      PyErr_SetString(PyExc_ValueError, "PDF not allocated");
      return NULL;
    }

  try
  {
    res = _pdfs[0]->q2Min();
  }
  catch (LHAPDF::Exception e)
  {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  }

  return Py_BuildValue("d", res);
}

static PyObject* py_pdfreplica(PyObject* self, PyObject* args)
{
  int nrep;
  PyArg_ParseTuple(args, "i", &nrep);
  _imem = nrep;

  return Py_BuildValue("");
}



// https://docs.python.org/3.4/extending/extending.html#calling-python-functions-from-c

static PyObject *
py_set_custom_xfxQ(PyObject *dummy, PyObject *args)
{
    PyObject *result = NULL;
    PyObject *temp;

    if (PyArg_ParseTuple(args, "O:set_custom_xfxQ", &temp)) {
        if (temp == Py_None){
           temp = NULL;
	} else if (!PyCallable_Check(temp)) {
            PyErr_SetString(PyExc_TypeError, "parameter must be callable or None");
            return NULL;
        }
        Py_XINCREF(temp);         /* Add a reference to new callback */
        Py_XDECREF(_custom_xfxq);  /* Dispose of previous callback */
        _custom_xfxq = temp;       /* Remember new callback */
        /* Boilerplate to return "None" */
        Py_INCREF(Py_None);
        result = Py_None;
    }
    return result;
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
  vector<double> xsec;
  PyArg_ParseTuple(args,"i", &pto);

  if (!_g) exit(-1);
  if (_custom_xfxq){
     xsec = _g->vconvolute(evolvepdf_python_,alphaspdf_,pto);
  }else{
     xsec = _g->vconvolute(evolvepdf_,alphaspdf_,pto);
  }
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
  if (_g->calculation() == appl::grid::AMCATNLO) // if aMCfast change iorder
    iorder = (pto == 0) ? 3:0;

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

static PyObject* py_setlhapdfpath(PyObject* self, PyObject* args)
{
  char *file;
  PyArg_ParseTuple(args,"s", &file);
  LHAPDF::setPaths(file);

  return Py_BuildValue("");
}

static PyObject* py_getlhapdfpath(PyObject* self, PyObject* args)
{
  vector<string> res = LHAPDF::paths();

  PyObject *out = PyList_New(res.size());
  for (int i = 0; i < (int) res.size(); i++)
    PyList_SET_ITEM(out, i, PyUnicode_FromString(res[i].c_str()));

  return out;
}

static PyObject* py_lhapdf_version(PyObject* self,  PyObject* noargs){
   string version = LHAPDF::version();
   return Py_BuildValue("s#", version.c_str(), version.length());
}

static PyMethodDef applwrap_methods[] = {
  {"setlhapdfpath", py_setlhapdfpath, METH_VARARGS, "set lhapdf path"},
  {"getlhapdfpath", py_getlhapdfpath, METH_VARARGS, "get lhapdf path"},
  {"setverbosity",  py_setverbosity,  METH_VARARGS, "set verbosity"},
  {"initpdf", py_initpdf, METH_VARARGS, "init pdf"},
  {"xfxQ", py_xfxQ, METH_VARARGS, "get xfxQ"},
  {"set_custom_xfxQ", py_set_custom_xfxQ, METH_VARARGS, "Set custom xfxQ for PDF"},
  {"q2Min", py_q2Min, METH_VARARGS, "get q2min"},
  {"pdfreplica", py_pdfreplica, METH_VARARGS, "set pdf replica"},
  {"initobs", py_initobs, METH_VARARGS, "init obs"},
  {"convolute", py_convolute, METH_VARARGS, "convolute"},
  {"getobsq", py_getobsq, METH_VARARGS, "get observable q"},
  {"getnbins",py_getnbins, METH_VARARGS, "get number of bins"},
  {"lhapdf_version",py_lhapdf_version, METH_NOARGS, "get LHAPDF version"},
  {NULL, NULL, 0, NULL}
};

static struct PyModuleDef applwrap_module = {
  PyModuleDef_HEAD_INIT,
  "applwrap",
  "A appl wrapper",
  -1,
  applwrap_methods
};

PyMODINIT_FUNC
PyInit_applwrap(void)
{
  return PyModule_Create(&applwrap_module);
}
