#include "find_chessboard_corners.hh"
#include "mrgingham.hh"

#include "mrgingham_pywrap_cplusplus_bridge.h"
#include "mrgingham-internal.h"

// These are wrappers. And the reason this file exists at all is because
// mrgingham has a c++ api, while the python wrappers are in C. AND there's
// another complicating factor: on CentOS7.2 if I build a C++ thing that talks
// to numpy it barfs:
//
//   g++ -Wall -Wextra -std=gnu99 -Wno-cast-function-type -I/usr/include/opencv -fvisibility=hidden -Wno-unused-function -Wno-missing-field-initializers -Wno-unused-parameter -Wno-strict-aliasing -Wno-int-to-pointer-cast -Wno-unused-variable -DOLD_OPENCV=1 -MMD -MP -g -fno-omit-frame-pointer -DVERSION='"1.1"' -pthread -fno-strict-aliasing -g -pipe -Wall -Wp,-D_FORTIFY_SOURCE=2 -fexceptions -fstack-protector-strong --param=ssp-buffer-size=4 -grecord-gcc-switches -m64 -mtune=generic -D_GNU_SOURCE -fPIC -fwrapv -DNDEBUG -g -pipe -Wall -Wp,-D_FORTIFY_SOURCE=2 -fexceptions -fstack-protector-strong --param=ssp-buffer-size=4 -grecord-gcc-switches -m64 -mtune=generic -D_GNU_SOURCE -fPIC -fwrapv -fPIC -I/usr/include/python2.7 -O3 -c -o/dev/null mrgingham_pywrap_cplusplus_bridge.cc  |& head -n 50
//   cc1plus: warning: command line option '-std=gnu99' is valid for C/ObjC but not for C++ [enabled by default]
//   In file included from /usr/include/numpy/ndarraytypes.h:7:0,
//                    from /usr/include/numpy/ndarrayobject.h:17,
//                    from /usr/include/numpy/arrayobject.h:15,
//                    from mrgingham_pywrap_cplusplus_bridge.cc:5:
//   /usr/include/numpy/npy_common.h:188:2: error: #error Must use Python with unicode enabled.
//    #error Must use Python with unicode enabled.
//     ^
//
// As a result, this file doesn't touch Python, and does stuff with callbacks.
// Yuck. Callbacks will make stuff slower. With a newer distro, the python
// unicode thing maybe isn't a problem. You can revert the commit on 2019/3/25
// if you want to try it again

extern "C"
bool find_chessboard_corners_from_image_array_C( // in
                                                int Nrows, int Ncols,
                                                int stride,
                                                char* imagebuffer, // this is const

                                                // set to 0 to just use the image
                                                int image_pyramid_level,

                                                bool (*init)(int N),
                                                bool (*add_point)(int i, double x, double y) )
{
    cv::Mat cvimage(Nrows, Ncols, CV_8UC1,
                    imagebuffer, stride);

    std::vector<mrgingham::PointInt> out_points;

    bool result = find_chessboard_corners_from_image_array( &out_points, cvimage, image_pyramid_level );
    if( !result ) return false;

    if(!(*init)(out_points.size())) return false;

    for(int i=0; i<(int)out_points.size(); i++)
        if(!( (*add_point)( i,
                            (double)out_points[i].x / (double)FIND_GRID_SCALE,
                            (double)out_points[i].y / (double)FIND_GRID_SCALE)))
            return false;
    return true;
}

extern "C"
bool find_chessboard_from_image_array_C( // in
                                        int Nrows, int Ncols,
                                        int stride,
                                        char* imagebuffer, // this is const

                                        // set to 0 to just use the image. Set
                                        // to <0 to try automatically find a
                                        // good scaling level. Try this first
                                        int image_pyramid_level,

                                        bool (*init)(int N),
                                        bool (*add_point)(int i, double x, double y) )
{
    cv::Mat cvimage(Nrows, Ncols, CV_8UC1,
                    imagebuffer, stride);

    std::vector<mrgingham::PointDouble> out_points;

    bool result = find_chessboard_from_image_array( out_points, cvimage, image_pyramid_level );
    if( !result ) return false;

    if(!(*init)(out_points.size())) return false;

    for(int i=0; i<(int)out_points.size(); i++)
        if(!( (*add_point)( i,
                            out_points[i].x, out_points[i].y)))
            return false;
    return true;
}
