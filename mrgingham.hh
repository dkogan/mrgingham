#pragma once

#include <opencv2/core/core.hpp>
#include <vector>
#include "point.hh"


namespace mrgingham
{

    struct debug_sequence_t
    {
        bool     dodebug;
        PointInt pt;
        debug_sequence_t() :
            dodebug(false),
            pt()
        {}
    };

    // This is currently hard-coded to find black-on-white blobs (look at
    // blobColor in find_blobs.cc)
    bool find_circle_grid_from_image_array( std::vector<mrgingham::PointDouble>& points_out,
                                            const cv::Mat& image,
                                            const int gridn,
                                            bool debug = false,
                                            debug_sequence_t debug_sequence = debug_sequence_t());
    bool find_circle_grid_from_image_file( std::vector<mrgingham::PointDouble>& points_out,
                                           const char* filename,
                                           const int gridn,
                                           bool debug = false,
                                           debug_sequence_t debug_sequence = debug_sequence_t());

    // set image_pyramid_level=0 to just use the image as is.
    //
    // image_pyramid_level > 0 cut down the image by a factor of 2 that many
    // times. So for example, level==2 means each dimension is cut down by a
    // factor of 4.
    //
    // image_pyramid_level < 0 means we try several levels, taking the first one
    // that produces results
    //
    // If we want to refine each reported point, pass a pointer to a buffer into
    // refinement_level. I'll realloc() the buffer as needed, and I'll return
    // the pyramid level of each point on exit.
    //
    // *refinement_level is managed by realloc(). IT IS THE CALLER'S
    // *RESPONSIBILITY TO free() IT
    //
    // Returns the pyramid level where we found the grid, or <0 on failure
    int  find_chessboard_from_image_array( std::vector<mrgingham::PointDouble>& points_out,
                                           signed char**                        refinement_level,
                                           const int                            gridn,
                                           const cv::Mat&                       image,
                                           int                                  image_pyramid_level  = -1,
                                           bool                                 debug                = false,
                                           debug_sequence_t                     debug_sequence = debug_sequence_t(),
                                           const char*                          debug_image_filename = NULL);

    // set image_pyramid_level=0 to just use the image as is.
    //
    // image_pyramid_level > 0 cut down the image by a factor of 2 that many
    // times. So for example, level==2 means each dimension is cut down by a
    // factor of 4.
    //
    // image_pyramid_level < 0 means we try several levels, taking the first one
    // that produces results
    //
    // If we want to refine each reported point, pass a pointer to a buffer into
    // refinement_level. I'll realloc() the buffer as needed, and I'll return
    // the pyramid level of each point on exit.
    //
    // *refinement_level is managed by realloc(). IT IS THE CALLER'S
    // *RESPONSIBILITY TO free() IT
    //
    // Returns the pyramid level where we found the grid, or <0 on failure
    int find_chessboard_from_image_file( std::vector<mrgingham::PointDouble>& points_out,
                                         signed char**                        refinement_level,
                                         const int                            gridn,
                                         const char*                          filename,
                                         int                                  image_pyramid_level = -1,
                                         bool                                 debug               = false,
                                         debug_sequence_t                     debug_sequence = debug_sequence_t());

    bool find_grid_from_points( std::vector<mrgingham::PointDouble>& points_out,
                                const std::vector<mrgingham::PointInt>& points,
                                const int gridn,
                                bool debug = false,
                                const debug_sequence_t& debug_sequence = debug_sequence_t());
};
