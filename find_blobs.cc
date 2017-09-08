#include <opencv2/features2d/features2d.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "point.hh"


void find_blobs_from_image_array( std::vector<Point>* points,
                                  const cv::Mat& image )
{
    cv::SimpleBlobDetector::Params blobDetectorParams;
    blobDetectorParams.maxArea *= 16;

    cv::SimpleBlobDetector blobDetector(blobDetectorParams);

    std::vector<cv::KeyPoint> keypoints;
    blobDetector.detect(image, keypoints);

    for(std::vector<cv::KeyPoint>::iterator it = keypoints.begin();
        it != keypoints.end();
        it++)
    {
        points->push_back( Point(it->pt.x, it->pt.y));
    }
}

bool find_blobs_from_image_file( std::vector<Point>* points,
                                 const char* filename )
{
    cv::Mat image = cv::imread(filename, CV_LOAD_IMAGE_COLOR);
    if( image.data == NULL )
        return false;

    find_blobs_from_image_array( points, image );
    return true;
}
