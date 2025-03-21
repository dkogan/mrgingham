#pragma once

#ifdef __cplusplus
extern "C" {
#endif

// These are wrappers. And the reason this file exists at all is because
// mrgingham has a c++ api, while the python wrappers are in C

bool find_chessboard_corners_from_image_array_C( // in
                                                int Nrows, int Ncols,
                                                int stride,
                                                char* imagebuffer, // this is const

                                                // set to 0 to just use the image. Set
                                                // to <0 to try automatically find a
                                                // good scaling level. Try this first
                                                int image_pyramid_level,
                                                bool doblobs,
                                                bool debug,

                                                bool (*add_points)(int* xy, int N, double scale, void* cookie),
                                                void* cookie );

bool find_chessboard_from_image_array_C( // in
                                        int Nrows, int Ncols,
                                        int stride,
                                        char* imagebuffer, // this is const

                                        const int gridn,

                                        // set to 0 to just use the image. Set
                                        // to <0 to try automatically find a
                                        // good scaling level. Try this first
                                        int image_pyramid_level,
                                        bool doblobs,
                                        bool debug,
                                        int debug_sequence_x,
                                        int debug_sequence_y,

                                        bool (*add_points)(double* xy, int N, void* cookie),
                                        void* cookie );

#ifdef __cplusplus
}
#endif
