#include "mrgingham.hh"
#include "find_blobs.hh"
#include "find_chessboard_corners.hh"

#include <opencv2/highgui/highgui.hpp>


namespace mrgingham
{
    __attribute__((visibility("default")))
    bool find_circle_grid_from_image_array( std::vector<PointDouble>& points_out,
                                            const cv::Mat& image,
                                            bool debug)
    {
        std::vector<Point> points;
        find_blobs_from_image_array(&points, image);
        return find_grid_from_points(points_out, points, debug);
    }

    __attribute__((visibility("default")))
    bool find_circle_grid_from_image_file( std::vector<PointDouble>& points_out,
                                           const char* filename,
                                           bool debug)
    {
        std::vector<Point> points;
        find_blobs_from_image_file(&points, filename);
        return find_grid_from_points(points_out, points, debug);
    }

    // same as below, but with a valid image_pyramid_level. The main function
    // will try several of these if image_pyramid_level<0
    static bool _find_chessboard_from_image_array( std::vector<PointDouble>& points_out,
                                                   const cv::Mat& image,
                                                   int image_pyramid_level,
                                                   bool do_refine,
                                                   bool debug)
    {
        std::vector<Point> points;
        find_chessboard_corners_from_image_array(&points, image, image_pyramid_level, debug);
        return find_grid_from_points(points_out, points, debug);

    }

    __attribute__((visibility("default")))
    bool find_chessboard_from_image_array( std::vector<PointDouble>& points_out,
                                           const cv::Mat& image,
                                           int image_pyramid_level,
                                           bool do_refine,
                                           bool debug)

    {
        if( image_pyramid_level >= 0)
            return _find_chessboard_from_image_array( points_out,
                                                      image,
                                                      image_pyramid_level,
                                                      do_refine,
                                                      debug);

        for( image_pyramid_level=2; image_pyramid_level>=0; image_pyramid_level--)
        {
            bool result =
                _find_chessboard_from_image_array( points_out,
                                                   image,
                                                   image_pyramid_level,
                                                   do_refine,
                                                   debug);
            if(result)
                return true;
        }
        return false;
    }

    __attribute__((visibility("default")))
    bool find_chessboard_from_image_file( std::vector<PointDouble>& points_out,
                                          const char* filename,
                                          int image_pyramid_level,
                                          bool do_refine,
                                          bool debug)
    {
        cv::Mat image = cv::imread(filename, CV_LOAD_IMAGE_GRAYSCALE);
        if( image.data == NULL )
        {
            fprintf(stderr, "%s:%d in %s(): Couldn't open image '%s'."
                    " Sorry.\n", __FILE__, __LINE__, __func__, filename);
            return false;
        }

        std::vector<Point> points;
        return find_chessboard_from_image_array(points_out, image, image_pyramid_level, do_refine, debug);
    }
};
