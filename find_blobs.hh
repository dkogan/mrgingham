#pragma once

#include <vector>
#include <opencv2/core/core.hpp>
#include "point.hh"


// I look for white-on-black dots


// these all output the points scaled by FIND_GRID_SCALE
bool find_blobs_from_image_array( std::vector<mrgingham::Point>* points,
                                  const cv::Mat& image,
                                  bool dodump = false);
bool find_blobs_from_image_file( std::vector<mrgingham::Point>* points,
                                 const char* filename,
                                 bool dodump = false);
