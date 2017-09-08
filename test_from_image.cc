#include "find_grid.hh"
#include "find_blobs.hh"
#include <stdio.h>

using namespace mrgingham;

static bool find_grid_from_image_array( std::vector<PointDouble>& points_out,
                                        const cv::Mat& image )
{
    std::vector<Point> points;
    find_blobs_from_image_array(&points, image);
    return find_grid_from_points(points_out, points);
}

static bool find_grid_from_image_file( std::vector<PointDouble>& points_out,
                                       const char* filename )
{
    std::vector<Point> points;
    find_blobs_from_image_file(&points, filename);
    return find_grid_from_points(points_out, points);
}

int main(int argc, char* argv[])
{
    if( argc != 2 )
    {
        fprintf(stderr, "missing arg: need image filename on the cmdline\n");
        return 1;
    }

    std::vector<PointDouble>& points_out;
    bool result = find_grid_from_image_file(points_out, argv[1]);

    if( result )
    {
        for(int i=0; i<points_out.size(); i++)
            printf("%f %f\n", points_out[i].x, points_out[i].y);
        return 0;
    }
    return 1;
}
