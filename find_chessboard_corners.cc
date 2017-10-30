#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "point.hh"
#include "mrgingham_internal.h"

extern "C"
{
#include "ChESS.h"
}


#define QUALITY_LEVEL(x) ((double)(((unsigned int)x) >> 6))


using namespace mrgingham;

bool find_chessboard_corners_from_image_array( std::vector<Point>* points,
                                               const cv::Mat& image,
                                               bool dodump )
{
    if( !image.isContinuous() )
    {
        fprintf(stderr, "%s:%d in %s(): I can only handle continuous arrays (stride == width) currently."
                " Sorry.\n", __FILE__, __LINE__, __func__);
        return false;
    }

    if( image.type() != CV_8U )
    {
        fprintf(stderr, "%s:%d in %s(): I can only handle CV_8U arrays currently."
                " Sorry.\n", __FILE__, __LINE__, __func__);
        return false;
    }

    int w = image.cols;
    int h = image.rows;

    uint8_t* imageData = image.data;

    cv::Mat response( h, w, CV_16S );
    int16_t* responseData = (int16_t*)response.data;


    ChESS_response_5( responseData, imageData, w, h );


    // Alrighty. I have responses. The following is heuristic. I do mostly what
    // the post-response-map cv::goodFeaturesToTrack() does
    int16_t maxVal = 0;
    for(int i=0; i<w*h; i++)
        if( responseData[i] > maxVal) maxVal = responseData[i];

    if( maxVal <= 0 )
    {
        fprintf(stderr, "%s:%d in %s(): All responses negative..."
                " Sorry.\n", __FILE__, __LINE__, __func__);
        return false;
    }

    cv::Mat responseDilated;
    cv::threshold( response, response, QUALITY_LEVEL(maxVal), 0, cv::THRESH_TOZERO );

#define DILATION_R 1 /* Mat() below implies this */
    cv::dilate( response, responseDilated, cv::Mat() );
    int16_t* responseDilatedData = (int16_t*)responseDilated.data;

    // The dilation results are maxima over the given window. Any points where
    // dilation == image are local maxima that I grab. I make sure to avoid the
    // edges where the dilation kernel didn't fit comletely
    for( int y = DILATION_R; y < h - DILATION_R; y++ )
        for( int x = DILATION_R; x < w - DILATION_R; x++ )
        {
            int16_t val         = responseData[ x + y*w ];
            int16_t val_dilated = responseDilatedData[ x + y*w ];
            if( val != 0 && val == val_dilated )
            {
                if( dodump )
                    printf("%d %d\n", x, y);
                else
                    points->push_back( Point(x * FIND_GRID_SCALE,
                                             y * FIND_GRID_SCALE));
            }
        }

    return true;
}

bool find_chessboard_corners_from_image_file( std::vector<Point>* points,
                                              const char* filename,
                                              bool dodump )
{
    cv::Mat image = cv::imread(filename, CV_LOAD_IMAGE_COLOR);
    if( image.data == NULL )
        return false;

    return find_chessboard_corners_from_image_array( points, image, dodump );
}
