#pragma once

#include <vector>
#include <opencv2/core/core.hpp>
#include "point.hh"



// these all output the points scaled by FIND_GRID_SCALE in points[].
bool find_chessboard_corners_from_image_array( // out

                                               // integers scaled up by
                                               // FIND_GRID_SCALE to get more
                                               // resolution
                                               std::vector<mrgingham::PointInt>* points,

                                               // in
                                               const cv::Mat& image,

                                               // set to 0 to just use the
                                               // image. Values > 0 cut down the
                                               // image by a factor of 2 that
                                               // many times. So for example,
                                               // level==2 means each dimension
                                               // is cut down by a factor of 4
                                               int image_pyramid_level,
                                               bool debug = false);

bool find_chessboard_corners_from_image_file( // out

                                              // integers scaled up by
                                              // FIND_GRID_SCALE to get more
                                              // resolution
                                              std::vector<mrgingham::PointInt>* points,

                                              // in
                                              const char* filename,

                                              // set to 0 to just use the
                                              // image. Values > 0 cut down the
                                              // image by a factor of 2 that
                                              // many times. So for example,
                                              // level==2 means each dimension
                                              // is cut down by a factor of 4
                                              int image_pyramid_level,
                                              bool debug = false);
