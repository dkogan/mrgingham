#pragma once

#include <vector>
#include <opencv2/core/core.hpp>
#include "point.hh"


// This is currently hard-coded to find black-on-white blobs (look at
// blobColor in find_blobs.cc)

namespace mrgingham
{

// these all output the points scaled by FIND_GRID_SCALE
bool find_blobs_from_image_array( std::vector<mrgingham::PointInt>* points,
                                  const cv::Mat& image,
                                  bool dodump = false);
bool find_blobs_from_image_file( std::vector<mrgingham::PointInt>* points,
                                 const char* filename,
                                 bool dodump = false);
}
