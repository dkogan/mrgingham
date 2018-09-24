#include <opencv2/features2d/features2d.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "point.hh"
#include "mrgingham-internal.h"

using namespace mrgingham;

__attribute__((visibility("default")))
bool find_blobs_from_image_array( std::vector<Point>* points,
                                  const cv::Mat& image,
                                  bool dodump )
{
    cv::SimpleBlobDetector::Params blobDetectorParams;
    blobDetectorParams.minArea             = 40;
    blobDetectorParams.maxArea             = 80000;
    blobDetectorParams.minDistBetweenBlobs = 15;
    blobDetectorParams.blobColor           = 0; // black-on-white dots


    cv::SimpleBlobDetector* blobDetector =
        cv::SimpleBlobDetector::create(blobDetectorParams);

    std::vector<cv::KeyPoint> keypoints;
    blobDetector->detect(image, keypoints);

    for(std::vector<cv::KeyPoint>::iterator it = keypoints.begin();
        it != keypoints.end();
        it++)
    {
        if( dodump )
        {
            printf("%f %f\n", it->pt.x, it->pt.y);
        }
        else
        {
            points->push_back( Point((int)(it->pt.x * FIND_GRID_SCALE + 0.5),
                                     (int)(it->pt.y * FIND_GRID_SCALE + 0.5)));
        }
    }

    return true;
}

__attribute__((visibility("default")))
bool find_blobs_from_image_file( std::vector<Point>* points,
                                 const char* filename,
                                 bool dodump )
{
    cv::Mat image = cv::imread(filename, CV_LOAD_IMAGE_GRAYSCALE);
    if( image.data == NULL )
    {
        fprintf(stderr, "%s:%d in %s(): Couldn't open image '%s'."
                " Sorry.\n", __FILE__, __LINE__, __func__, filename);
        return false;
    }

    return find_blobs_from_image_array( points, image, dodump );
}
