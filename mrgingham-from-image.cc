#include "mrgingham.hh"
#include <stdio.h>
#include <getopt.h>
#include <glob.h>
#include <pthread.h>

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace mrgingham;

struct mrgingham_thread_context_t
{
    const glob_t* _glob;
    int           Njobs;
    bool          doclahe;
    int           blur_radius;
    bool          doblobs;
    bool          do_refine;
    int           gridn;
    bool          debug;
    debug_sequence_t debug_sequence;
    int           image_pyramid_level;
} ctx;

static void* worker( void* _ijob )
{
    // Worker thread. Processes images from the glob. Writes point detections
    // back out on the other end.

    int ijob = *(int*)(&_ijob);

    cv::Ptr<cv::CLAHE> clahe;

    if(ctx.doclahe)
    {
        clahe = cv::createCLAHE();
        clahe->setClipLimit(8);
    }

    // The buffer. I'll realloc() this as I go. MUST free at the end
    signed char* refinement_level = NULL;

    for(int i_image=ijob; i_image<(int)ctx._glob->gl_pathc; i_image += ctx.Njobs)
    {
        const char* filename = ctx._glob->gl_pathv[i_image];

        cv::Mat image = cv::imread(filename,
                                   cv::IMREAD_IGNORE_ORIENTATION |
                                   cv::IMREAD_GRAYSCALE);
        if( image.data == NULL )
        {
            fprintf(stderr, "Couldn't open image '%s'\n", filename);
            flockfile(stdout);
            {
                printf("## Couldn't open image '%s'\n", filename);
                printf("%s - - -\n", filename);
            }
            funlockfile(stdout);
            break;
        }

        if( ctx.doclahe )
        {
            // CLAHE doesn't by itself use the full dynamic range all the time.
            // I explicitly apply histogram equalization and then CLAHE
            cv::equalizeHist(image, image);
            clahe->apply(image, image);
        }
        if( ctx.blur_radius > 0 )
        {
            cv::blur( image, image,
                      cv::Size(1 + 2*ctx.blur_radius,
                               1 + 2*ctx.blur_radius));
        }

        if( ctx.debug )
        {
            do
            {
                char basename[1024];

                const char* last_slash = strrchr(filename, '/');
                const char* basename_start = last_slash ? &last_slash[1] : filename;

                char* basename_start_end = stpncpy(basename, basename_start, sizeof(basename));
                if(&basename[sizeof(basename)] == basename_start_end)
                {
                    fprintf(stderr, "--debug file dump overran filename buffer! Not dumping files\n");
                    break;
                }

                // basename is now just the FILENAME with no directory. It still has
                // an extension
                char* last_dot = strrchr(basename, '.');
                if(last_dot)
                    *last_dot = '\0';

                char filename_out[1024];
                if(snprintf(filename_out, sizeof(filename_out),
                            "/tmp/%s_preprocessed.png",
                            basename) >= (int)sizeof(filename_out))
                {
                    fprintf(stderr, "--debug file dump overran filename buffer! Not dumping files\n");
                    break;
                }

                cv::imwrite(filename_out, image);
                fprintf(stderr, "Wrote preprocessed image to %s\n", filename_out);

            } while(0);
        }
        std::vector<PointDouble> points_out;
        bool result;
        int found_pyramid_level; // need this because ctx.image_pyramid_level could be -1

        if(ctx.doblobs)
        {
            result = find_circle_grid_from_image_array(points_out,
                                                       image, ctx.gridn,
                                                       ctx.debug, ctx.debug_sequence);
            // ctx.image_pyramid_level == 0 here. cmdline parser makes sure.
            found_pyramid_level = 0;
        }
        else
        {
            found_pyramid_level =
                find_chessboard_from_image_array (points_out,
                                                  ctx.do_refine ? &refinement_level : NULL,
                                                  ctx.gridn,
                                                  image,
                                                  ctx.image_pyramid_level,
                                                  ctx.debug, ctx.debug_sequence,
                                                  filename);
            result = (found_pyramid_level >= 0);
        }

        flockfile(stdout);
        {
            if( result )
            {
                for(int i=0; i<(int)points_out.size(); i++)
                    printf( "%s %f %f %d\n", filename,
                            points_out[i].x,
                            points_out[i].y,
                            (refinement_level == NULL) ? found_pyramid_level : (int)refinement_level[i]);
            }
            else
                printf("%s - - -\n", filename);
        }
        funlockfile(stdout);
    }

    free(refinement_level);

    return NULL;
}

int main(int argc, char* argv[])
{
    const char* usage =
        "Usage: %s [--debug] [--debug-sequence x,y]\n"
        "                   [--jobs N] [--noclahe] [--blur radius] [--gridn N]\n"
        "                   [--level l] [--blobs] imageglobs imageglobs ...\n"
        "\n"
        "  By default we look for a chessboard. By default we apply adaptive histogram\n"
        "  equalization, then blur with a radius of 1. We then use an adaptive level of\n"
        "  downsampling when looking for the chessboard.\n"
        "\n"
        "  We will look for chessboards, unless --blobs is given\n"
        "\n"
        "  --noclahe is optional: unless given, we will pre-process the image with an\n"
        "  adaptive histogram equalization step (using the CLAHE algorithm). This is\n"
        "  useful if the calibration board has a lighting gradient across it.\n"
        "\n"
        "  --blur radius   applies a blur (after CLAHE) to the image before processing.\n"
        "  By default we will blur with a radius of 1. To disable, set the radius to <= 0\n"
        "\n"
        "  --level l   applies a downsampling to the image before processing it (after\n"
        "  CLAHE and --blur, if given) to the image before processing. Level 0 means\n"
        "  'use the original image'. Level > 0 means downsample by 2**level. Level < 0\n"
        "  means 'try several different levels until we find one that works. This is the\n"
        "  default.\n"
        "\n"
        "  --no-refine  By default, the coordinates of reported corners are re-detected at\n"
        "  less-downsampled zoom levels to improve their accuracy. If we do not want to do\n"
        "  that, pass --no-refine\n"
        "\n"
        "  --jobs N  will parallelize the processing N-ways. -j is a synonym. This is like\n"
        "  GNU make, except you're required to explicitly specify a job count.\n"
        "\n"
        "  We detect an NxN grid of corners, where N defaults to 10. To select a different\n"
        "  value, pass --gridn N\n"
        "\n"
        "  The images are given as (multiple) globs. The output is a vnlog with columns\n"
        "  filename,x,y. All filenames matched in the glob will appear in the output.\n"
        "  Images for which no chessboard pattern was found appear as a single record\n"
        "  with null x and y.\n"
        "\n"
        "  For debugging, pass in --debug. This will dump the various intermediate results\n"
        "  into /tmp and it will report more stuff on the console.\n"
        "  For debugging sequence candidates, pass in --debug-sequence x,y where 'x,y' are the\n"
        "  approximate image coordinates of the start of a given sequence (corner on the edge of\n"
        "  a chessboard. This doesn't need to be exact; mrgingham will report on the nearest\n"
        "  corner\n";

    struct option opts[] = {
        { "blobs",             no_argument,       NULL, 'B' },
        { "blur",              required_argument, NULL, 'b' },
        { "noclahe",           no_argument,       NULL, 'C' },
        { "level",             required_argument, NULL, 'l' },
        { "no-refine",         no_argument,       NULL, 'R' },
        { "jobs",              required_argument, NULL, 'j' },
        { "gridn",             required_argument, NULL, 'N' },
        { "debug",             no_argument,       NULL, 'd' },
        { "debug-sequence",    required_argument, NULL, 'D' },
        { "help",              no_argument,       NULL, 'h' },
        {}
    };


    bool        doblobs             = false;
    bool        doclahe             = true;
    bool        do_refine           = true;
    bool        debug               = false;
    bool        debug_sequence      = false;
    PointInt    debug_sequence_pt;
    int         blur_radius         = 1;
    int         image_pyramid_level = -1;
    int         jobs                = 1;
    int         gridn               = 10;

    int opt;
    do
    {
        // "h" means -h does something
        opt = getopt_long(argc, argv, "hj:b:l:", opts, NULL);
        switch(opt)
        {
        case -1:
            break;

        case 'h':
            printf(usage, argv[0]);
            return 0;

        case 'B':
            doblobs = true;
            break;

        case 'C':
            doclahe = false;
            break;

        case 'R':
            do_refine = false;
            break;

        case 'd':
            debug = true;
            break;

        case 'D':
            debug_sequence = true;
            if( 2 != sscanf(optarg, "%d,%d",
                            &debug_sequence_pt.x,
                            &debug_sequence_pt.y) )
            {
                fprintf(stderr, "I could not parse 'x,y' from --debug-sequence '%s'. Giving up\n",
                        optarg);
                fprintf(stderr, usage, argv[0]);
                return -1;
            }
            break;

        case 'N':
            gridn = atoi(optarg);
            break;

        case 'b':
            blur_radius = atoi(optarg);
            break;

        case 'l':
            image_pyramid_level = atoi(optarg);
            break;

        case 'j':
            jobs = atoi(optarg);
            break;

        case '?':
            fprintf(stderr, "Unknown option\n");
            fprintf(stderr, usage, argv[0]);
            return 1;
        }
    } while( opt != -1 );

    if( optind > argc-1)
    {
        fprintf(stderr, "Not enough arguments: need image globs\n");
        fprintf(stderr, usage, argv[0]);
        return 1;
    }
    if( jobs <= 0 )
    {
        fprintf(stderr, "The job count must be a positive integer\n");
        fprintf(stderr, usage, argv[0]);
        return 1;
    }
    if( doblobs && image_pyramid_level >= 0)
    {
        fprintf(stderr, "ERROR: 'image_pyramid_level' only implemented for chessboards.\n");
        return 1;
    }
    if(gridn < 2)
    {
        fprintf(stderr, "--gridn value must be >= 2\n");
        return 1;
    }

    glob_t _glob;
    int doappend = 0;
    for( int iopt_glob = optind; iopt_glob<argc; iopt_glob++ )
    {
        const char* imageglob = argv[iopt_glob];
        int globresult =
            glob(imageglob,
                 doappend |
                 GLOB_ERR | GLOB_MARK | GLOB_NOSORT | GLOB_TILDE_CHECK,
                 NULL, &_glob);
        if(globresult == GLOB_NOMATCH)
        {
            fprintf(stderr, "'%s' matched no files!\n", imageglob);
            return 1;
        }
        if(globresult != 0)
        {
            fprintf(stderr, "globbing '%s' failed!\n", imageglob);
            return 1;
        }

        doappend = GLOB_APPEND;
    }

    if(debug && _glob.gl_pathc != 1)
    {
        fprintf(stderr, "When debugging, pass one image at a time. Got %d instead\n",
                (int)_glob.gl_pathc);
        return 1;
    }


    printf("## generated with");
    for(int i=0; i<argc; i++)
        printf(" %s", argv[i]);
    printf("\n");

    printf("# filename x y level\n");

    // I'm done with the preliminaries. I now spawn the child threads. Note that
    // in this implementation it is important that these are THREADS and not a
    // fork. I want to make sure that the image output is atomic. To do that I
    // use flockfile(), and each child thread writes directly to stdout.
    // flockfile() does not work in a fork, but does work in a thread
    ctx._glob               = &_glob;
    ctx.Njobs               = jobs;
    ctx.doclahe             = doclahe;
    ctx.blur_radius         = blur_radius;
    ctx.doblobs             = doblobs;
    ctx.do_refine           = do_refine;
    ctx.gridn               = gridn;
    ctx.debug               = debug;

    ctx.debug_sequence.dodebug = debug_sequence;
    ctx.debug_sequence.pt      = debug_sequence_pt;

    ctx.image_pyramid_level = image_pyramid_level;

    pthread_t thread[jobs];
    for(int i=0; i<jobs; i++)
        pthread_create(&thread[i], NULL, &worker, (void*)i);

    for(int i=0; i<jobs; i++)
        pthread_join(thread[i], NULL);

    globfree(&_glob);
    return 0;
}
