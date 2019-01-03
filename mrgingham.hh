#pragma once

#include <opencv2/core/core.hpp>
#include <vector>
#include "point.hh"


// I look for white-on-black dots


namespace mrgingham
{
    bool find_circle_grid_from_image_array( std::vector<mrgingham::PointDouble>& points_out,
                                            const cv::Mat& image, bool debug = false );
    bool find_circle_grid_from_image_file( std::vector<mrgingham::PointDouble>& points_out,
                                           const char* filename, bool debug = false );

    // set image_pyramid_level=0 to just use the image as is.
    //
    // image_pyramid_level > 0 cut down the image by a factor of 2 that many
    // times. So for example, level==2 means each dimension is cut down by a
    // factor of 4.
    //
    // image_pyramid_level < 0 means we try several levels, taking the first one
    // that produces results
    bool find_chessboard_from_image_array( std::vector<mrgingham::PointDouble>& points_out,
                                           const cv::Mat& image,
                                           int image_pyramid_level = -1,
                                           bool do_refine = true,
                                           bool debug = false);

    // set image_pyramid_level=0 to just use the image as is.
    //
    // image_pyramid_level > 0 cut down the image by a factor of 2 that many
    // times. So for example, level==2 means each dimension is cut down by a
    // factor of 4.
    //
    // image_pyramid_level < 0 means we try several levels, taking the first one
    // that produces results
    bool find_chessboard_from_image_file( std::vector<mrgingham::PointDouble>& points_out,
                                          const char* filename,
                                          int image_pyramid_level = -1,
                                          bool do_refine = true,
                                          bool debug = false );

    bool find_grid_from_points( std::vector<mrgingham::PointDouble>& points_out,
                                const std::vector<mrgingham::PointInt>& points,
                                bool debug = false);
};
