#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <stdio.h>
#include <getopt.h>
#include "find_chessboard_corners.hh"

using namespace mrgingham;

int main(int argc, char* argv[])
{
    const char* usage =
        "Usage: %s [--clahe] [--blur radius] [--level l] image\n"
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
        "  means 'try several different levels until we find one that works. This is the.\n"
        "  default.\n"
        "\n";

    struct option opts[] = {
        { "blur",    required_argument, NULL, 'b' },
        { "clahe",   no_argument,       NULL, 'c' },
        { "level",   required_argument, NULL, 'l' },
        { "help",    no_argument,       NULL, 'h' },
        {}
    };


    bool        doclahe             = false;
    int         blur_radius         = -1;
    int         image_pyramid_level = -1;

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

        case 'c':
            doclahe = true;
            break;

        case 'b':
            blur_radius = atoi(optarg);
            if(blur_radius <= 0)
            {
                fprintf(stderr, "blur_radius < 0\n");
                fprintf(stderr, usage, argv[0]);
                return 1;

            }
            break;

        case 'l':
            image_pyramid_level = atoi(optarg);
            break;

        case '?':
            fprintf(stderr, "Unknown option\n");
            fprintf(stderr, usage, argv[0]);
            return 1;
        }
    } while( opt != -1 );

    if( optind != argc-1)
    {
        fprintf(stderr, "Need a single image on the cmdline\n");
        fprintf(stderr, usage, argv[0]);
        return 1;
    }


    cv::Ptr<cv::CLAHE> clahe;

    if(doclahe)
    {
        clahe = cv::createCLAHE();
        clahe->setClipLimit(8);
    }

    const char* filename = argv[argc-1];

    cv::Mat image = cv::imread(filename, CV_LOAD_IMAGE_GRAYSCALE);
    if( image.data == NULL )
    {
        fprintf(stderr, "Couldn't open image '%s'\n", filename);
        return 1;
    }

    if( doclahe )
        clahe->apply(image, image);
    if( blur_radius > 0 )
    {
        cv::blur( image, image,
                  cv::Size(1 + 2*blur_radius,
                           1 + 2*blur_radius));
    }

    std::vector<PointInt> points;
    if(!find_chessboard_corners_from_image_array (&points, image, image_pyramid_level, true, filename))
    {
        fprintf(stderr, "Error computing the corners!\n");
        return 1;
    }

    return 0;
}
