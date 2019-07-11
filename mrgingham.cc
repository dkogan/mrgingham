#include "mrgingham.hh"
#include "find_blobs.hh"
#include "find_chessboard_corners.hh"

#include <opencv2/highgui/highgui.hpp>


namespace mrgingham
{
    __attribute__((visibility("default")))
    bool find_circle_grid_from_image_array( std::vector<PointDouble>& points_out,
                                            const cv::Mat& image,
                                            bool     debug,
                                            debug_sequence_t debug_sequence)
    {
        std::vector<PointInt> points;
        find_blobs_from_image_array(&points, image);
        return find_grid_from_points(points_out, points,
                                     debug, debug_sequence);
    }

    __attribute__((visibility("default")))
    bool find_circle_grid_from_image_file( std::vector<PointDouble>& points_out,
                                           const char* filename,
                                           bool     debug,
                                           debug_sequence_t debug_sequence)
    {
        std::vector<PointInt> points;
        find_blobs_from_image_file(&points, filename);
        return find_grid_from_points(points_out, points,
                                     debug, debug_sequence);
    }

    // *refinement_level is managed by realloc(). IT IS THE CALLER'S
    // *RESPONSIBILITY TO free() IT
    static bool _find_chessboard_from_image_array( std::vector<PointDouble>& points_out,
                                                   signed char** refinement_level,
                                                   const cv::Mat& image,
                                                   int image_pyramid_level,
                                                   bool     debug,
                                                   debug_sequence_t debug_sequence,
                                                   const char* debug_image_filename)
    {
        const bool do_refine = (refinement_level != NULL);

        std::vector<PointInt> points;
        find_chessboard_corners_from_image_array(&points, image, image_pyramid_level, debug, debug_image_filename);
        if(!find_grid_from_points(points_out, points,
                                  debug, debug_sequence))
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
        //           if(this point is refinable)
        //             refine();
        //       if( no points remain refinable )
        //           break;
        //   }

        int N = points_out.size();
        *refinement_level = (signed char*)realloc((void*)*refinement_level, N*sizeof(**refinement_level));
        assert(*refinement_level);
        for(int i=0; i<N; i++)
            (*refinement_level)[i] = (signed char)image_pyramid_level;

        while(image_pyramid_level--)
        {
            int Nrefined =
                mrgingham::
                refine_chessboard_corners_from_image_array( &points_out,
                                                            *refinement_level,
                                                            image, image_pyramid_level,
                                                            debug, debug_image_filename);
            if(debug)
                fprintf(stderr, "Refining to level %d... Nrefined=%d\n", image_pyramid_level, Nrefined);
            if(Nrefined <= 0)
                break;
        }
        return true;
    }

    // *refinement_level is managed by realloc(). IT IS THE CALLER'S
    // *RESPONSIBILITY TO free() IT
    __attribute__((visibility("default")))
    int find_chessboard_from_image_array( std::vector<PointDouble>& points_out,
                                          signed char** refinement_level,
                                          const cv::Mat& image,
                                          int image_pyramid_level,
                                          bool debug,
                                          debug_sequence_t debug_sequence,
                                          const char* debug_image_filename)

    {
        if( image_pyramid_level >= 0)
            return
                _find_chessboard_from_image_array( points_out,
                                                   refinement_level,
                                                   image,
                                                   image_pyramid_level,
                                                   debug, debug_sequence,
                                                   debug_image_filename)
                ? image_pyramid_level : -1;

        for( image_pyramid_level=3; image_pyramid_level>=0; image_pyramid_level--)
        {
            int result = _find_chessboard_from_image_array( points_out,
                                                            refinement_level,
                                                            image,
                                                            image_pyramid_level,
                                                            debug, debug_sequence,
                                                            debug_image_filename)
                ? image_pyramid_level : -1;
            if(result >= 0) return result;
        }
        return -1;
    }

    // *refinement_level is managed by realloc(). IT IS THE CALLER'S
    // *RESPONSIBILITY TO free() IT
    __attribute__((visibility("default")))
    int find_chessboard_from_image_file( std::vector<PointDouble>& points_out,
                                         signed char** refinement_level,
                                         const char* filename,
                                         int image_pyramid_level,
                                         bool debug,
                                         debug_sequence_t debug_sequence)
    {
        cv::Mat image = cv::imread(filename, CV_LOAD_IMAGE_GRAYSCALE);
        if( image.data == NULL )
        {
            fprintf(stderr, "%s:%d in %s(): Couldn't open image '%s'."
                    " Sorry.\n", __FILE__, __LINE__, __func__, filename);
            return -1;
        }

        std::vector<PointInt> points;
        return find_chessboard_from_image_array(points_out,
                                                refinement_level,
                                                image, image_pyramid_level,
                                                debug, debug_sequence,
                                                filename);
    }
};
