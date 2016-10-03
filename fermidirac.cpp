#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include "Python.h"
#include "numpy/ndarraytypes.h"
#include "numpy/ufuncobject.h"
#include "numpy/halffloat.h"

#include <cmath>

#include "fermidirac.hpp"


static PyMethodDef FDMethods[] = {
        {NULL, NULL, 0, NULL}
};

// F_{1/2}(x)

static void float_fd_pos(char **args, npy_intp *dimensions,
                         npy_intp* steps, void* data) {
  FermiDirac<float> fd;
  npy_intp n = dimensions[0];

  char *in = args[0], *out = args[1];
  npy_intp in_step = steps[0], out_step = steps[1];

  #pragma omp parallel for
  for( npy_intp i=0; i<n; ++i ) {
    float tmp = *(float*)in;
    *((float*)out) = fd.fermi_half(tmp);
    in  += in_step;
    out += out_step;
  }
}

static void double_fd_pos(char **args, npy_intp *dimensions,
                          npy_intp* steps, void* data) {
  FermiDirac<double> fd;
  npy_intp n = dimensions[0];

  char *in = args[0], *out = args[1];
  npy_intp in_step = steps[0], out_step = steps[1];

  #pragma omp parallel for
  for( npy_intp i=0; i<n; ++i ) {
    double tmp = *(double*)in;
    *((double*)out) = fd.fermi_half(tmp);
    in  += in_step;
    out += out_step;
  }
}

static void longdouble_fd_pos(char **args, npy_intp *dimensions,
                              npy_intp* steps, void* data) {
  FermiDirac<long double> fd;
  npy_intp n = dimensions[0];

  char *in = args[0], *out = args[1];
  npy_intp in_step = steps[0], out_step = steps[1];

  #pragma omp parallel for
  for( npy_intp i=0; i<n; ++i ) {
    long double tmp = *(long double*)in;
    *((float*)out) = fd.fermi_half(tmp);
    in  += in_step;
    out += out_step;
  }
}


// F_{0}(x)

static void float_fd_zero(char **args, npy_intp *dimensions,
                          npy_intp* steps, void* data) {
  FermiDirac<float> fd;
  npy_intp n = dimensions[0];

  char *in = args[0], *out = args[1];
  npy_intp in_step = steps[0], out_step = steps[1];

  #pragma omp parallel for
  for( npy_intp i=0; i<n; ++i ) {
    float tmp = *(float*)in;
    *((float*)out) = fd.fermi_zero(tmp);
    in  += in_step;
    out += out_step;
  }
}

static void double_fd_zero(char **args, npy_intp *dimensions,
                           npy_intp* steps, void* data) {
  FermiDirac<double> fd;
  npy_intp n = dimensions[0];

  char *in = args[0], *out = args[1];
  npy_intp in_step = steps[0], out_step = steps[1];

  #pragma omp parallel for
  for( npy_intp i=0; i<n; ++i ) {
    double tmp = *(double*)in;
    *((double*)out) = fd.fermi_zero(tmp);
    in  += in_step;
    out += out_step;
  }
}

static void longdouble_fd_zero(char **args, npy_intp *dimensions,
                               npy_intp* steps, void* data) {
  FermiDirac<long double> fd;
  npy_intp n = dimensions[0];

  char *in = args[0], *out = args[1];
  npy_intp in_step = steps[0], out_step = steps[1];

  #pragma omp parallel for
  for( npy_intp i=0; i<n; ++i ) {
    long double tmp = *(long double*)in;
    *((float*)out) = fd.fermi_zero(tmp);
    in  += in_step;
    out += out_step;
  }
}

// F_{-1/2}(x)

static void float_fd_neg(char **args, npy_intp *dimensions,
                         npy_intp* steps, void* data) {
  FermiDirac<float> fd;
  npy_intp n = dimensions[0];

  char *in = args[0], *out = args[1];
  npy_intp in_step = steps[0], out_step = steps[1];

  #pragma omp parallel for
  for( npy_intp i=0; i<n; ++i ) {
    float tmp = *(float*)in;
    *((float*)out) = fd.fermi_minus_half(tmp);
    in  += in_step;
    out += out_step;
  }
}

static void double_fd_neg(char **args, npy_intp *dimensions,
                          npy_intp* steps, void* data) {
  FermiDirac<double> fd;
  npy_intp n = dimensions[0];

  char *in = args[0], *out = args[1];
  npy_intp in_step = steps[0], out_step = steps[1];

  #pragma omp parallel for
  for( npy_intp i=0; i<n; ++i ) {
    double tmp = *(double*)in;
    *((double*)out) = fd.fermi_minus_half(tmp);
    in  += in_step;
    out += out_step;
  }
}

static void longdouble_fd_neg(char **args, npy_intp *dimensions,
                              npy_intp* steps, void* data) {
  FermiDirac<long double> fd;
  npy_intp n = dimensions[0];

  char *in = args[0], *out = args[1];
  npy_intp in_step = steps[0], out_step = steps[1];

  #pragma omp parallel for
  for( npy_intp i=0; i<n; ++i ) {
    long double tmp = *(long double*)in;
    *((float*)out) = fd.fermi_minus_half(tmp);
    in  += in_step;
    out += out_step;
  }
}

PyUFuncGenericFunction pos_funcs[3]  = {&float_fd_pos,  &double_fd_pos,  &longdouble_fd_pos};
PyUFuncGenericFunction zero_funcs[3] = {&float_fd_zero, &double_fd_zero, &longdouble_fd_zero};
PyUFuncGenericFunction neg_funcs[3]  = {&float_fd_neg,  &double_fd_neg,  &longdouble_fd_neg};

static char types[6] = { NPY_FLOAT,      NPY_FLOAT,
                         NPY_DOUBLE,     NPY_DOUBLE,
                         NPY_LONGDOUBLE, NPY_LONGDOUBLE };

static void *data[1] = {NULL};

#if PY_VERSION_HEX >= 0x03000000
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "fermidirac",
    NULL,
    -1,
    FDMethods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC PyInit_fermidirac(void)
{
    PyObject *m, *half, *zero, *minus_half, *d;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }

    import_array();
    import_umath();

    half = PyUFunc_FromFuncAndData(pos_funcs, data, types, 3, 1, 1,
                                   PyUFunc_None, "half",
                                    "half_docstring", 0);

    zero = PyUFunc_FromFuncAndData(zero_funcs, data, types, 3, 1, 1,
                                   PyUFunc_None, "zero",
                                    "zero_docstring", 0);

    minus_half = PyUFunc_FromFuncAndData(neg_funcs, data, types, 3, 1, 1,
                                         PyUFunc_None, "minus_half",
                                         "half_docstring", 0);

    d = PyModule_GetDict(m);

    PyDict_SetItemString(d, "half", half);
    PyDict_SetItemString(d, "zero", zero);
    PyDict_SetItemString(d, "minus_half", minus_half);

    Py_DECREF(half);
    Py_DECREF(zero);
    Py_DECREF(minus_half);

    return m;
}
#else
PyMODINIT_FUNC initfermidirac(void)
{
    PyObject *m, *half, *zero, *minus_half, *d;
    m = Py_InitModule("fermidirac", FDMethods);
    if (m == NULL) {
        return;
    }

    import_array();
    import_umath();

    half = PyUFunc_FromFuncAndData(pos_funcs, data, types, 3, 1, 1,
                                   PyUFunc_None, "half",
                                    "half_docstring", 0);

    zero = PyUFunc_FromFuncAndData(zero_funcs, data, types, 3, 1, 1,
                                   PyUFunc_None, "zero",
                                    "zero_docstring", 0);

    minus_half = PyUFunc_FromFuncAndData(neg_funcs, data, types, 3, 1, 1,
                                         PyUFunc_None, "minus_half",
                                         "half_docstring", 0);

    axpy = PyUFunc_FromFuncAndData(funcs, data, types, 3, 1, 1,
                                   PyUFunc_None, "axpy",
                                   "axpy_docstring", 0);

    d = PyModule_GetDict(m);

    PyDict_SetItemString(d, "half", half);
    PyDict_SetItemString(d, "zero", zero);
    PyDict_SetItemString(d, "minus_half", minus_half);

    Py_DECREF(half);
    Py_DECREF(zero);
    Py_DECREF(minus_half);
}
#endif





