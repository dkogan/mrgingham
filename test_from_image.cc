#include "find_grid.hh"
#include "find_blobs.hh"
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

static bool find_grid_from_image_array( const cv::Mat& image )
{
    std::vector<Point> points;
    find_blobs_from_image_array(&points, image);
    return find_grid_from_points(points);
}

static bool find_grid_from_image_file( const char* filename )
{
    std::vector<Point> points;
    find_blobs_from_image_file(&points, filename);
    return find_grid_from_points(points);
}

int main(int argc, char* argv[])
{
    if( argc != 2 )
    {
        fprintf(stderr, "missing arg: need image filename on the cmdline\n");
        return 1;
    }

    return find_grid_from_image_file(argv[1]) ? 0 : 1;
}
