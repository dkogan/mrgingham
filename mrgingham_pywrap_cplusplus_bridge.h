#pragma once

#include <Python.h>

#ifdef __cplusplus
extern "C" {
#endif

// These are wrappers. And the reason this file exists at all is because
// mrgingham has a c++ api, while the python wrappers are in C

PyObject* find_chessboard_corners_from_image_array_C( // in
                                                     PyObject* image,

                                                     // set to 0 to just use the image
                                                     int image_pyramid_level );

PyObject* find_chessboard_from_image_array_C( // in
                                             PyObject* image,

                                             // set to 0 to just use the image
                                             int image_pyramid_level );

#ifdef __cplusplus
}
#endif
