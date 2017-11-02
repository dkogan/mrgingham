#include "mrgingham.hh"
#include <stdio.h>
#include <getopt.h>

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace mrgingham;

int main(int argc, char* argv[])
{
    const char* usage =
        "Usage: %s [--clahe] [--blur radius] [--level l] --blobs|--chessboard imagefile\n"
        "\n"
        "  --blobs or --chessboard are required; these tell the tool what to do.\n"
        "\n"
        "  --clahe is optional: it will pre-process the image with an adaptive histogram\n"
        "  equalization step. This is useful if the calibration board has a lighting\n"
        "  gradient across it.\n"
        "\n"
        "  --blur radius   applies a blur (after --clahe, if given) to the image before\n"
        "  processing\n"
        "\n"

        "  --level l   applies a downsampling to the image before processing it (after\n"
        "  --clahe and --blur, if given) to the image before processing. Level 0 means\n"
        "  'use the original image'. Level > 0 means downsample by 2**level. Level < 0\n"
        "  means 'try several different levels until we find one that works'.\n";

    struct option opts[] = {
        { "blobs",      no_argument,       NULL, 'B' },
        { "chessboard", no_argument,       NULL, 'C' },
        { "blur",       required_argument, NULL, 'b' },
        { "clahe",      no_argument,       NULL, 'c' },
        { "level",      required_argument, NULL, 'l' },
        { "help",       no_argument,       NULL, 'h' },
        {}
    };


    bool have_doblobs        = false;
    bool doblobs             = false; // to pacify compiler
    bool doclahe             = false;
    int  blur_radius         = -1;
    int  image_pyramid_level = -1;

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

        case 'B':
        case 'C':
            if( have_doblobs )
            {
                fprintf(stderr, usage, argv[0]);
                return 1;

            }
            have_doblobs = true;
            doblobs      = (opt == 'B');
            break;

        case 'c':
            doclahe = true;
            break;

        case 'b':
            blur_radius = atoi(optarg);
            if(blur_radius <= 0)
            {
                fprintf(stderr, usage, argv[0]);
                return 1;

            }
            break;

        case 'l':
            image_pyramid_level = atoi(optarg);
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

    if( doblobs && image_pyramid_level >= 0)
    {
        fprintf(stderr, "warning: 'image_pyramid_level' only implemented for chessboards. Will be ignored\n");
    }

    const char* filename = argv[argc-1];


    cv::Mat image = cv::imread(filename, CV_LOAD_IMAGE_GRAYSCALE);
    if( image.data == NULL )
    {
        fprintf(stderr, "Couldn't open image '%s'\n", filename);
        return 1;
    }

    if( doclahe )
    {
        cv::Ptr<cv::CLAHE> clahe = cv::createCLAHE();
        clahe->setClipLimit(8);

        clahe->apply(image, image);
    }
    if( blur_radius > 0 )
    {
        cv::blur( image, image,
                  cv::Size(1 + 2*blur_radius,
                           1 + 2*blur_radius));
    }

    std::vector<PointDouble> points_out;
    bool result;
    if(doblobs)
        result = find_circle_grid_from_image_array(points_out, image);
    else
        result = find_chessboard_from_image_array (points_out, image, image_pyramid_level);

    if( result )
    {
        for(int i=0; i<(int)points_out.size(); i++)
            printf("%f %f\n", points_out[i].x, points_out[i].y);
        return 0;
    }
    return 1;
}
