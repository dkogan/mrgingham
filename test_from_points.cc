#include "find_grid.hh"
#include <stdio.h>

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

    return find_grid_from_points(points) ? 0 : 1;
}
