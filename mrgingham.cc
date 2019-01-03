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
        std::vector<PointInt> points;
        find_blobs_from_image_array(&points, image);
        return find_grid_from_points(points_out, points, debug);
    }

    __attribute__((visibility("default")))
    bool find_circle_grid_from_image_file( std::vector<PointDouble>& points_out,
                                           const char* filename,
                                           bool debug)
    {
        std::vector<PointInt> points;
        find_blobs_from_image_file(&points, filename);
        return find_grid_from_points(points_out, points, debug);
    }

    static bool _find_chessboard_from_image_array( std::vector<PointDouble>& points_out,
                                                   const cv::Mat& image,
                                                   int image_pyramid_level,
                                                   bool do_refine,
                                                   bool debug,
                                                   const char* debug_image_filename)
    {
        std::vector<PointInt> points;
        find_chessboard_corners_from_image_array(&points, image, image_pyramid_level, debug, debug_image_filename);
        if(!find_grid_from_points(points_out, points, debug))
            return false;

        // we found a grid! If we're not trying to refine the locations, or if
        // we can't refine them, we're done
        if(!do_refine || image_pyramid_level == 0)
            return true;

        // Alright, I need to refine each intersection. Big-picture logic:
        //
        //   for(points)
        //   {
        //       zoom = next_from_current;
        //       while(update corner coord using zoom)
        //           zoom = next_from_current;
        //   }
        //
        // It would be more efficient to loop through the zoom levels once, so I
        // move that to the outer loop
        //
        //   for(zoom)
        //   {
        //       for(points)
        //       {
        //           if(!this point is refinable) continue;
        //           refine();
        //       }
        //       if( no points remain refinable )
        //           break;
        //   }

        int N = points_out.size();
        bool refinable[N];
        for(int i=0; i<N; i++) refinable[i] = true;

        while(image_pyramid_level--)
        {
            int Nrefined =
                refine_chessboard_corners_from_image_array( &points_out,
                                                            refinable,
                                                            image, image_pyramid_level,
                                                            debug, debug_image_filename);
            if(debug)
                fprintf(stderr, "Refining to level %d... Nrefined=%d\n", image_pyramid_level, Nrefined);
            if(Nrefined <= 0)
                break;
        }
        return true;
    }

    __attribute__((visibility("default")))
    bool find_chessboard_from_image_array( std::vector<PointDouble>& points_out,
                                           const cv::Mat& image,
                                           int image_pyramid_level,
                                           bool do_refine,
                                           bool debug,
                                           const char* debug_image_filename)

    {
        if( image_pyramid_level >= 0)
            return _find_chessboard_from_image_array( points_out,
                                                      image,
                                                      image_pyramid_level,
                                                      do_refine,
                                                      debug, debug_image_filename);

        for( image_pyramid_level=2; image_pyramid_level>=0; image_pyramid_level--)
        {
            bool result =
                _find_chessboard_from_image_array( points_out,
                                                   image,
                                                   image_pyramid_level,
                                                   do_refine,
                                                   debug, debug_image_filename);
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

        std::vector<PointInt> points;
        return find_chessboard_from_image_array(points_out, image, image_pyramid_level, do_refine, debug, filename);
    }
};
