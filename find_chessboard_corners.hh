#pragma once

#include <vector>
#include <opencv2/core/core.hpp>
#include "point.hh"



// these all output the points scaled by FIND_GRID_SCALE in points[].
bool find_chessboard_corners_from_image_array( // out

                                               // integers scaled up by
                                               // FIND_GRID_SCALE to get more
                                               // resolution
                                               std::vector<mrgingham::PointInt>* points_scaled_out,

                                               // in
                                               const cv::Mat& image_input,

                                               // set to 0 to just use the
                                               // image. Values > 0 cut down the
                                               // image by a factor of 2 that
                                               // many times. So for example,
                                               // level==2 means each dimension
                                               // is cut down by a factor of 4
                                               int image_pyramid_level,
                                               bool debug = false,
                                               const char* debug_image_filename = NULL);

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

int refine_chessboard_corners_from_image_array( // out/int

                                                // initial coordinates on input,
                                                // refined coordinates on output
                                                std::vector<mrgingham::PointDouble>* points,

                                                // if(!point_refinable[ipoint])
                                                // then that point shouldn't be
                                                // refined. If we try and fail
                                                // to refine a point, we set
                                                // point_refinable[ipoint] to
                                                // false
                                                bool* point_refinable,

                                                // in
                                                const cv::Mat& image_input,

                                                int image_pyramid_level,
                                                bool debug = false,
                                                const char* debug_image_filename = NULL);
