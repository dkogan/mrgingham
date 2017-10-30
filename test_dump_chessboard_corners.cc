#include "find_chessboard_corners.hh"
#include <stdio.h>

using namespace mrgingham;

int main(int argc, char* argv[])
{
    if( argc != 2 )
    {
        fprintf(stderr, "missing arg: need image filename on the cmdline\n");
        return 1;
    }

    std::vector<Point> points;
    find_chessboard_corners_from_image_file(&points, argv[1], true);
    return 0;
}
