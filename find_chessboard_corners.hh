#pragma once

#include <vector>
#include <opencv2/core/core.hpp>
#include "point.hh"



// these all output the points scaled by FIND_GRID_SCALE
bool find_chessboard_corners_from_image_array( std::vector<mrgingham::Point>* points,
                                               const cv::Mat& image,
                                               bool dodump = false);
bool find_chessboard_corners_from_image_file( std::vector<mrgingham::Point>* points,
                                              const char* filename,
                                              bool dodump = false);
