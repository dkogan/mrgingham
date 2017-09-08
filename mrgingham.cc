#include "mrgingham.hh"
#include "find_blobs.hh"


namespace mrgingham
{
    bool find_circle_grid_from_image_array( std::vector<PointDouble>& points_out,
                                            const cv::Mat& image )
    {
        std::vector<Point> points;
        find_blobs_from_image_array(&points, image);
        return find_grid_from_points(points_out, points);
    }

    bool find_circle_grid_from_image_file( std::vector<PointDouble>& points_out,
                                           const char* filename )
    {
        std::vector<Point> points;
        find_blobs_from_image_file(&points, filename);
        return find_grid_from_points(points_out, points);
    }
};
