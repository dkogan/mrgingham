#define NPY_NO_DEPRECATED_API NPY_API_VERSION

#include <stdbool.h>
#include <Python.h>
#include <structmember.h>
#include <numpy/arrayobject.h>
#include <signal.h>

#include "ChESS.h"

// Python is silly. There's some nuance about signal handling where it sets a
// SIGINT (ctrl-c) handler to just set a flag, and the python layer then reads
// this flag and does the thing. Here I'm running C code, so SIGINT would set a
// flag, but not quit, so I can't interrupt the solver. Thus I reset the SIGINT
// handler to the default, and put it back to the python-specific version when
// I'm done
#define SET_SIGINT() struct sigaction sigaction_old;                    \
do {                                                                    \
    if( 0 != sigaction(SIGINT,                                          \
                       &(struct sigaction){ .sa_handler = SIG_DFL },    \
                       &sigaction_old) )                                \
    {                                                                   \
        PyErr_SetString(PyExc_RuntimeError, "sigaction() failed");      \
        goto done;                                                      \
    }                                                                   \
} while(0)
#define RESET_SIGINT() do {                                             \
    if( 0 != sigaction(SIGINT,                                          \
                       &sigaction_old, NULL ))                          \
        PyErr_SetString(PyExc_RuntimeError, "sigaction-restore failed"); \
} while(0)

#define PYMETHODDEF_ENTRY(function_prefix, name, args) {#name,          \
                                                        (PyCFunction)function_prefix ## name, \
                                                        args,           \
                                                        function_prefix ## name ## _docstring}


static PyObject* py_ChESS_response_5(PyObject* NPY_UNUSED(self),
                                     PyObject* args)
{
    PyObject* result = NULL;
    SET_SIGINT();

    PyArrayObject* image = NULL;
    if(!PyArg_ParseTuple( args, "O&", PyArray_Converter, &image ))
        goto done;

    npy_intp* dims    = PyArray_DIMS   (image);
    npy_intp* strides = PyArray_STRIDES(image);
    int       ndims   = PyArray_NDIM   (image);
    if( ndims < 2 )
    {
        PyErr_Format(PyExc_RuntimeError, "The input image array must have at least 2 dims (extra ones will be broadcasted); got %d",
                     PyArray_NDIM(image));
        goto done;
    }
    if( PyArray_TYPE(image) != NPY_UINT8 )
    {
        PyErr_SetString(PyExc_RuntimeError, "The input image array must contain 8-bit unsigned data");
        goto done;
    }
    if( strides[ndims-1] != 1 )
    {
        PyErr_SetString(PyExc_RuntimeError, "Image rows must live in contiguous memory");
        goto done;
    }

    // This doesn't actually need to be true, but assuming it makes the
    // broadcasting code below much simpler
    if( !PyArray_IS_C_CONTIGUOUS(image) )
    {
        PyErr_SetString(PyExc_RuntimeError, "The input image array must be C-contiguous");
        goto done;
    }

    PyArrayObject* response =
        (PyArrayObject*)PyArray_SimpleNew(ndims, dims, NPY_INT16);
    if(response == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError, "Couldn't allocate reponse");
        goto done;
    }

    int Nslices = 1;
    for(int i=0; i<ndims-2; i++)
        Nslices *= dims[i];

    int imagesize = dims[ndims-1] * dims[ndims-2];
    for(int i=0; i<Nslices; i++)
    {
        // Everything is C-contiguous, so this is easy
        mrgingham_ChESS_response_5( (int16_t*)&PyArray_BYTES(response)[i * imagesize*2],
                                    (uint8_t*)&PyArray_BYTES(image)   [i * imagesize],
                                    dims[ndims-1], dims[ndims-2],
                                    strides[ndims-2]);
    }

    result = Py_BuildValue("O", response);
    Py_DECREF(response);

 done:
    Py_XDECREF(image);
    RESET_SIGINT();
    return result;
}

PyMODINIT_FUNC init_mrgingham(void)
{
    static const char py_ChESS_response_5_docstring[] =
#include "ChESS_response_5.docstring.h"
        ;
    static PyMethodDef methods[] =
        {
         PYMETHODDEF_ENTRY(py_,ChESS_response_5, METH_VARARGS),
         {}
        };

    PyObject* module = Py_InitModule3("_mrgingham", methods,
                                      "Chessboard-detection routines");
    import_array();
}
