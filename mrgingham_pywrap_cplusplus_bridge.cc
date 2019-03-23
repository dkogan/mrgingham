#define NPY_NO_DEPRECATED_API NPY_API_VERSION
#define PY_ARRAY_UNIQUE_SYMBOL mrgingham_
#define NO_IMPORT_ARRAY

#include <numpy/arrayobject.h>

#include "find_chessboard_corners.hh"
#include "mrgingham.hh"

#include "mrgingham_pywrap_cplusplus_bridge.h"
#include "mrgingham-internal.h"

// These are wrappers. And the reason this file exists at all is because
// mrgingham has a c++ api, while the python wrappers are in C

extern "C"
PyObject* find_chessboard_corners_from_image_array_C( // in
                                                     PyObject* image,

                                                     // set to 0 to just use the image
                                                     int image_pyramid_level )
{
    cv::Mat cvimage(PyArray_DIMS((PyArrayObject*)image)[0], PyArray_DIMS((PyArrayObject*)image)[1],
                    CV_8UC1,
                    PyArray_BYTES((PyArrayObject*)image),
                    PyArray_STRIDES((PyArrayObject*)image)[0]);

    std::vector<mrgingham::PointInt> out_points;

    bool result = find_chessboard_corners_from_image_array( &out_points, cvimage, image_pyramid_level );
    if( !result )
        return NULL;

    npy_intp dims[2] = {(npy_intp)out_points.size(), 2};
    PyObject* out = PyArray_SimpleNew(2, dims, NPY_DOUBLE);

    double* out_data = (double*)PyArray_BYTES((PyArrayObject*)out);
    for(int i=0; i<(int)out_points.size(); i++)
    {
        out_data[2*i + 0] = (double)out_points[i].x / (double)FIND_GRID_SCALE;
        out_data[2*i + 1] = (double)out_points[i].y / (double)FIND_GRID_SCALE;
    }

    return out;
}

extern "C"
PyObject* find_chessboard_from_image_array_C( // in
                                             PyObject* image,

                                             // set to 0 to just use the image
                                             int image_pyramid_level )
{
    cv::Mat cvimage(PyArray_DIMS((PyArrayObject*)image)[0], PyArray_DIMS((PyArrayObject*)image)[1],
                    CV_8UC1,
                    PyArray_BYTES((PyArrayObject*)image),
                    PyArray_STRIDES((PyArrayObject*)image)[0]);

    std::vector<mrgingham::PointDouble> out_points;

    bool result = mrgingham::find_chessboard_from_image_array( out_points, cvimage, image_pyramid_level );
    if( !result )
        return NULL;

    npy_intp dims[2] = {(npy_intp)out_points.size(), 2};
    PyObject* out = PyArray_SimpleNew(2, dims, NPY_DOUBLE);

    double* out_data = (double*)PyArray_BYTES((PyArrayObject*)out);
    for(int i=0; i<(int)out_points.size(); i++)
    {
        out_data[2*i + 0] = out_points[i].x;
        out_data[2*i + 1] = out_points[i].y;
    }

    return out;
}
