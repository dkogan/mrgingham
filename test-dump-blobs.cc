#include "find_blobs.hh"
#include <stdio.h>

using namespace mrgingham;

int main(int argc, char* argv[])
{
    if( argc != 2 )
    {
        fprintf(stderr, "missing arg: need image filename on the cmdline\n");
        return 1;
    }

    std::vector<PointInt> points;
    find_blobs_from_image_file(&points, argv[1], true);
    return 0;
}
