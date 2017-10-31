#include "find_chessboard_corners.hh"
#include <stdio.h>

using namespace mrgingham;

int main(int argc, char* argv[])
{
    const char* usage = "Usage: %s image_pyramid_level imagefilename\n";

    if( argc != 3 )
    {
        fprintf(stderr, usage, argv[0]);
        return 1;
    }

    int image_pyramid_level = atoi(argv[1]);
    if( image_pyramid_level < 0 || image_pyramid_level > 10)
    {
        fprintf(stderr, usage, argv[0]);
        return 1;
    }

    const char* filename = argv[2];

    std::vector<Point> points;
    find_chessboard_corners_from_image_file(&points, filename, image_pyramid_level, true);
    return 0;
}
