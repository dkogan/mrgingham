#pragma once

#include <vector>
#include <opencv2/core/core.hpp>
#include "point.hh"



// these all output the points scaled by FIND_GRID_SCALE in points[]. The debug
// dumping with dodump=true does NOT scale by FIND_GRID_SCALE
bool find_chessboard_corners_from_image_array( std::vector<mrgingham::Point>* points,
                                               const cv::Mat& image,

                                               // set to 0 to just use the
                                               // image. Values > 0 cut down the
                                               // image by a factor of 2 that
                                               // many times. So for example,
                                               // level==2 means each dimension
                                               // is cut down by a factor of 4
                                               int image_pyramid_level,
                                               bool dodump = false);

bool find_chessboard_corners_from_image_file( std::vector<mrgingham::Point>* points,
                                              const char* filename,

                                              // set to 0 to just use the
                                              // image. Values > 0 cut down the
                                              // image by a factor of 2 that
                                              // many times. So for example,
                                              // level==2 means each dimension
                                              // is cut down by a factor of 4
                                              int image_pyramid_level,
                                              bool dodump = false);
