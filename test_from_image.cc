#include "mrgingham.hh"
#include <stdio.h>

using namespace mrgingham;

int main(int argc, char* argv[])
{
    if( argc != 2 )
    {
        fprintf(stderr, "missing arg: need image filename on the cmdline\n");
        return 1;
    }

    std::vector<PointDouble> points_out;
    bool result = find_grid_from_image_file(points_out, argv[1]);

    if( result )
    {
        for(int i=0; i<(int)points_out.size(); i++)
            printf("%f %f\n", points_out[i].x, points_out[i].y);
        return 0;
    }
    return 1;
}
