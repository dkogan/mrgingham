#include "mrgingham.hh"
#include <stdio.h>
#include <getopt.h>

using namespace mrgingham;

int main(int argc, char* argv[])
{
    const char* usage = "Usage: %s --blobs|--chessboard imagefile\n";

    struct option opts[] = {
        { "blobs",      no_argument, NULL, 'B' },
        { "chessboard", no_argument, NULL, 'C' },
        { "help",       no_argument, NULL, 'h' },
        {}
    };


    bool have_doblobs = false;
    bool doblobs      = false; // to pacify compiler

    int opt;
    do
    {
        // "h" means -h does something
        opt = getopt_long(argc, argv, "h", opts, NULL);
        switch(opt)
        {
        case -1:
            break;

        case 'h':
            printf(usage, argv[0]);
            return 0;

        case 'C':
            if( have_doblobs )
            {
                fprintf(stderr, usage, argv[0]);
                return 1;

            }
            have_doblobs = true;
            doblobs      = (opt == 'B');
            break;

        case '?':
            fprintf(stderr, usage, argv[0]);
            return 1;
        }
    }  while( opt != -1 );

    if( !have_doblobs || optind != argc-1)
    {
        fprintf(stderr, usage, argv[0]);
        return 1;
    }

    const char* filename = argv[argc-1];

    std::vector<PointDouble> points_out;
    bool result;
    if(doblobs)
        result = find_circle_grid_from_image_file(points_out, filename);
    else
        result = find_chessboard_from_image_file (points_out, filename);

    if( result )
    {
        for(int i=0; i<(int)points_out.size(); i++)
            printf("%f %f\n", points_out[i].x, points_out[i].y);
        return 0;
    }
    return 1;
}
