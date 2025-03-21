#define NPY_NO_DEPRECATED_API NPY_API_VERSION

#include <stdbool.h>
#include <Python.h>
#include <structmember.h>
#include <numpy/arrayobject.h>
#include <signal.h>

#include "ChESS.h"
#include "mrgingham_pywrap_cplusplus_bridge.h"

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

#define PYMETHODDEF_ENTRY(name, c_function_name, args) {#name,          \
                                                        (PyCFunction)c_function_name, \
                                                        args,           \
                                                        name ## _docstring}


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
        PyErr_SetString(PyExc_RuntimeError, "Couldn't allocate response");
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

    result = (PyObject*)response;

 done:
    Py_XDECREF(image);
    RESET_SIGINT();
    return result;
}

static
bool add_points__find_points(int* xy, int N, double scale, void* cookie)
{
    PyObject** result = (PyObject**)cookie;
    *result = PyArray_SimpleNew(2,
                                ((npy_intp[]){N, 2}),
                                NPY_DOUBLE);
    if(*result == NULL) return false;

    double* out_data = (double*)PyArray_BYTES((PyArrayObject*)*result);
    for(int i=0; i<2*N; i++)
        out_data[i] = scale * (double)xy[i];
    return true;
}
static PyObject* find_points(PyObject* NPY_UNUSED(self),
                                         PyObject* args,
                                         PyObject* kwargs)
{
    PyArrayObject* image               = NULL;
    PyObject*      result              = NULL;
    int            image_pyramid_level = 0;
    int            blobs               = 0;
    int            debug               = 0;

    SET_SIGINT();

    char* keywords[] = { "image", "image_pyramid_level", "blobs",
                         "debug",
                         NULL };

    if(!PyArg_ParseTupleAndKeywords( args, kwargs,
                                     "O&|ipp",
                                     keywords,
                                     PyArray_Converter, &image,
                                     &image_pyramid_level,
                                     &blobs,
                                     &debug,
                                     NULL))
        goto done;

    if(blobs && image_pyramid_level != 0)
    {
        PyErr_Format(PyExc_RuntimeError, "blob detector requires that image_pyramid_level == 0");
        goto done;
    }

    npy_intp* dims    = PyArray_DIMS   (image);
    npy_intp* strides = PyArray_STRIDES(image);
    int       ndims   = PyArray_NDIM   (image);
    if( ndims != 2 )
    {
        PyErr_Format(PyExc_RuntimeError, "The input image array must have exactly 2 dims (broadcasting not supported here); got %d",
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

    if(! find_chessboard_corners_from_image_array_C(PyArray_DIMS(image)[0],
                                                    PyArray_DIMS(image)[1],
                                                    PyArray_STRIDES(image)[0],
                                                    PyArray_BYTES(image),

                                                    image_pyramid_level,
                                                    blobs,
                                                    debug,
                                                    &add_points__find_points, &result) )
    {
        if(result == NULL)
        {
            // Detector didn't find anything. I don't flag an error,
            // but simply return a no-data array
            result = PyArray_SimpleNew(2,
                                       ((npy_intp[]){0, 2}),
                                       NPY_DOUBLE);
        }
        else
        {
            // an actual error occured. I complain
            Py_DECREF(result);
            result = NULL;
            PyErr_SetString(PyExc_RuntimeError, "find_chessboard_corners_from_image_array_C() failed");
            goto done;
        }
    }

 done:
    Py_XDECREF(image);
    RESET_SIGINT();
    return result;
}

static
bool add_points__find_board(double* xy, int N, void* cookie)
{
    PyObject** result = (PyObject**)cookie;
    *result = PyArray_SimpleNew(2,
                                ((npy_intp[]){N, 2}),
                                NPY_DOUBLE);
    if(*result == NULL) return false;

    double* out_data = (double*)PyArray_BYTES((PyArrayObject*)*result);
    memcpy(out_data, xy, 2*N*sizeof(double));
    return true;
}
static PyObject* find_board(PyObject* NPY_UNUSED(self),
                                 PyObject* args,
                                 PyObject* kwargs)
{
    PyArrayObject* image               = NULL;
    PyObject*      result              = NULL;
    int            image_pyramid_level = -1;
    int            gridn               = 10;
    int            blobs               = 0;
    int            debug               = 0;
    const char*    debug_sequence_string = NULL;
    int debug_sequence_x = -1;
    int debug_sequence_y = -1;

    SET_SIGINT();

    char* keywords[] = { "image", "image_pyramid_level", "gridn", "blobs",
                         "debug", "debug_sequence",
                         NULL };

    if(!PyArg_ParseTupleAndKeywords( args, kwargs,
                                     "O&|iipps",
                                     keywords,
                                     PyArray_Converter, &image,
                                     &image_pyramid_level, &gridn,
                                     &blobs,
                                     &debug,
                                     &debug_sequence_string,
                                     NULL))
        goto done;

    if(blobs && image_pyramid_level != 0)
    {
        PyErr_Format(PyExc_RuntimeError, "blob detector requires that image_pyramid_level == 0");
        goto done;
    }

    if(debug_sequence_string != NULL)
    {
        // parse string into INTEGER,INTEGER
        int nbytesread;
        int res = sscanf(debug_sequence_string,
                         "%d,%d%n",
                         &debug_sequence_x,
                         &debug_sequence_x,
                         &nbytesread);

        if(!(res == 2 && debug_sequence_string[nbytesread] == '\0'))
        {
            PyErr_Format(PyExc_RuntimeError, "Couldn't parse debug_sequence as an 'INTEGER,INTEGER' string");
            goto done;

        }
    }


    npy_intp* dims    = PyArray_DIMS   (image);
    npy_intp* strides = PyArray_STRIDES(image);
    int       ndims   = PyArray_NDIM   (image);
    if( ndims != 2 )
    {
        PyErr_Format(PyExc_RuntimeError, "The input image array must have exactly 2 dims (broadcasting not supported here); got %d",
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
    if(gridn < 2)
    {
        PyErr_SetString(PyExc_RuntimeError, "gridn value must be >= 2");
        goto done;
    }

    if(! find_chessboard_from_image_array_C(PyArray_DIMS(image)[0],
                                            PyArray_DIMS(image)[1],
                                            PyArray_STRIDES(image)[0],
                                            PyArray_BYTES(image),

                                            gridn,
                                            image_pyramid_level,
                                            blobs,
                                            debug,
                                            debug_sequence_x,
                                            debug_sequence_y,

                                            &add_points__find_board, &result) )
    {
        // This is allowed to fail. We possibly found no chessboard. This is
        // sloppy since it ignore other potential errors, but there shouldn't be
        // any in this path
        Py_XDECREF(result);
        if( result == NULL )
        {
            result = Py_None;
            Py_INCREF(result);
        }
    }

 done:
    Py_XDECREF(image);
    RESET_SIGINT();
    return result;
}

static const char ChESS_response_5_docstring[] =
#include "ChESS_response_5.docstring.h"
    ;
static const char find_points_docstring[] =
#include "find_points.docstring.h"
    ;
static const char find_board_docstring[] =
#include "find_board.docstring.h"
    ;

static const char find_chessboard_corners_docstring[] =
#include "find_chessboard_corners.docstring.h"
    ;
static const char find_chessboard_docstring[] =
#include "find_chessboard.docstring.h"
    ;


static PyMethodDef methods[] =
    {
     PYMETHODDEF_ENTRY(ChESS_response_5,        py_ChESS_response_5, METH_VARARGS),
     PYMETHODDEF_ENTRY(find_points,             find_points,         METH_VARARGS | METH_KEYWORDS),
     PYMETHODDEF_ENTRY(find_board,              find_board,          METH_VARARGS | METH_KEYWORDS),

     // These are compatibility functions. Same implementation as the above, but
     // different names. Meant to keep old code running
     PYMETHODDEF_ENTRY(find_chessboard_corners, find_points,         METH_VARARGS | METH_KEYWORDS),
     PYMETHODDEF_ENTRY(find_chessboard,         find_board,          METH_VARARGS | METH_KEYWORDS),
     {}
    };

#if PY_MAJOR_VERSION == 2

__attribute__((visibility("default")))
PyMODINIT_FUNC initmrgingham(void)
{
    Py_InitModule3("mrgingham", methods,
                   "Chessboard-detection routines");
    import_array();
}

#else

static struct PyModuleDef module_def =
    {
     PyModuleDef_HEAD_INIT,
     "mrgingham",
     "Chessboard-detection routines",
     -1,
     methods
    };

PyMODINIT_FUNC PyInit_mrgingham(void)
{
    import_array();
    return PyModule_Create(&module_def);
}

#endif
