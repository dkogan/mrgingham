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

    PyArrayObject* response =
        (PyArrayObject*)PyArray_SimpleNew(ndims, dims, NPY_INT16);
    if(response == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError, "Couldn't allocate reponse");
        goto done;
    }

    // broadcast through all the slices
    { // need a block to pacify the compiler
        npy_intp islice[ndims];
        islice[ndims - 1] = 0;
        islice[ndims - 2] = 0;

        void loop_dim(int idim)
        {
            if(idim < 0)
            {
                int16_t* data_response = (int16_t*)PyArray_GetPtr(response, islice);
                uint8_t* data_image    = (uint8_t*)PyArray_GetPtr(image,    islice);

                mrgingham_ChESS_response_5( data_response, data_image,
                                            dims[ndims-1], dims[ndims-2],
                                            strides[ndims-2]);
                return;
            }

            for(islice[idim]=0; islice[idim] < dims[idim]; islice[idim]++)
                loop_dim(idim-1);
        }

        // The last 2 dimensions index each slice (x,y inside each image). The
        // dimensions before that are for broadcasting
        loop_dim(ndims-3);
    }

    result = Py_BuildValue("O", response);
    Py_DECREF(response);

 done:
    Py_XDECREF(image);
    RESET_SIGINT();
    return result;
}

PyMODINIT_FUNC initmrgingham(void)
{
    static const char py_ChESS_response_5_docstring[] =
#include "ChESS_response_5.docstring.h"
        ;
    static PyMethodDef methods[] =
        {
         PYMETHODDEF_ENTRY(py_,ChESS_response_5, METH_VARARGS),
         {}
        };

    PyObject* module = Py_InitModule3("mrgingham", methods,
                                      "Chessboard-detection routines");
    import_array();
}
