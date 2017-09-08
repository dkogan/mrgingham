#pragma once

#include <vector>
#include <opencv2/core/core.hpp>
#include "point.hh"

void find_blobs_from_image_array( std::vector<Point>* points,
                                  const cv::Mat& image,
                                  bool dodump = false);
bool find_blobs_from_image_file( std::vector<Point>* points,
                                 const char* filename,
                                 bool dodump = false);
