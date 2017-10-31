#include "mrgingham.hh"
#include <stdio.h>

using namespace mrgingham;

int main(int argc, char* argv[])
{
    const char* usage = "Usage: %s blobs|chessboard imagefile\n";

    if( argc != 3 )
    {
        fprintf(stderr, usage, argv[0]);
        return 1;
    }

    bool doblobs;
    if(      strcmp(argv[1], "blobs") == 0 )
        doblobs = true;
    else if( strcmp(argv[1], "chessboard") == 0 )
        doblobs = false;
    else
    {
        fprintf(stderr, usage, argv[0]);
        return 1;
    }

    std::vector<PointDouble> points_out;
    bool result;
    if(doblobs)
        result = find_circle_grid_from_image_file(points_out, argv[2]);
    else
        result = find_chessboard_from_image_file (points_out, argv[2]);

    if( result )
    {
        for(int i=0; i<(int)points_out.size(); i++)
            printf("%f %f\n", points_out[i].x, points_out[i].y);
        return 0;
    }
    return 1;
}
