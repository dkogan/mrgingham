#include "find_grid.hh"
#include <stdio.h>

using namespace mrgingham;

static bool read_points( std::vector<Point>* points, const char* file )
{
    FILE* fp = fopen(file, "r");
    if( fp == NULL )
    {
        fprintf(stderr, "couldn't open '%s'\n", file);
        return false;
    }

    while(1)
    {
        double x,y;
        int Nread = fscanf(fp, "%lf %lf", &x, &y);
        if(Nread != 2)
            break;

        Point pt( (int)( x * FIND_GRID_SCALE + 0.5 ),
                  (int)( y * FIND_GRID_SCALE + 0.5 ) );
        points->push_back(pt);
    }
    fclose(fp);
    return true;
}


int main(int argc, char* argv[])
{
    if( argc != 2 )
    {
        fprintf(stderr, "missing arg: need filename that contains ascii points\n");
        return 1;
    }

    std::vector<Point> points;
    if( !read_points(&points, argv[1]) )
        return 1;

    std::vector<PointDouble> points_out;
    bool result = find_grid_from_points(points_out, points);

    if( result )
    {
        for(int i=0; i<(int)points_out.size(); i++)
            printf("%f %f\n", points_out[i].x, points_out[i].y);
        return 0;
    }
    return 1;
}
