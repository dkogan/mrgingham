#pragma once

#include <opencv2/core/core.hpp>
#include <vector>
#include "point.hh"

namespace mrgingham
{
    bool find_circle_grid_from_image_array( std::vector<mrgingham::PointDouble>& points_out,
                                            const cv::Mat& image );
    bool find_circle_grid_from_image_file( std::vector<mrgingham::PointDouble>& points_out,
                                           const char* filename );
    bool find_grid_from_points( std::vector<mrgingham::PointDouble>& points_out,
                                const std::vector<mrgingham::Point>& points );
};
