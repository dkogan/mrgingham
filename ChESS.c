/*
  This is the reference implementation from this paper:

  https://arxiv.org/abs/1301.5491

  Dima made a few modifications.

  Obtained from here:

  http://www-sigproc.eng.cam.ac.uk/Main/SB476Chess

  There's a more full-featured GPL-licensed implementation on that page
*/

/**
 * The ChESS corner detection algorithm
 *
 * Copyright 2010-2012 Stuart Bennett
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#include <stdint.h>
#include <stdlib.h>

/**
 * Perform the ChESS corner detection algorithm with a 5 px sampling radius
 *
 * @param    response  output response image. Densely-packed
 *                     signed-16-bits-per-pixel image if size (w,h). Densely-
 *                     packed means the stride doesn't apply
 * @param    image     input image. Assumed 8 bits (1 byte) per pixel. Not
 *                     densely-packed: the stride applies
 * @param    w         image width
 * @param    h         image height
 * @param    stride    the length (in bytes) of each row in memory of the input
 *                     image. If stored densely, w == stride
 */
void ChESS_response_5(      int16_t* restrict response,
                      const uint8_t* restrict image,
                      int w, int h, int stride )
{
    int x, y;
    // funny bounds due to sampling ring radius (5) and border of previously applied blur (2)
    for (y = 7; y < h - 7; y++)
        for (x = 7; x < w - 7; x++) {
            const unsigned offset_input    = x + y * stride;
            const unsigned offset_response = x + y * w;
            uint8_t circular_sample[16];

            circular_sample[2] = image[offset_input - 2 - 5 * stride];
            circular_sample[1] = image[offset_input - 5 * stride];
            circular_sample[0] = image[offset_input + 2 - 5 * stride];
            circular_sample[8] = image[offset_input - 2 + 5 * stride];
            circular_sample[9] = image[offset_input + 5 * stride];
            circular_sample[10] = image[offset_input + 2 + 5 * stride];
            circular_sample[3] = image[offset_input - 4 - 4 * stride];
            circular_sample[15] = image[offset_input + 4 - 4 * stride];
            circular_sample[7] = image[offset_input - 4 + 4 * stride];
            circular_sample[11] = image[offset_input + 4 + 4 * stride];
            circular_sample[4] = image[offset_input - 5 - 2 * stride];
            circular_sample[14] = image[offset_input + 5 - 2 * stride];
            circular_sample[6] = image[offset_input - 5 + 2 * stride];
            circular_sample[12] = image[offset_input + 5 + 2 * stride];
            circular_sample[5] = image[offset_input - 5];
            circular_sample[13] = image[offset_input + 5];

            // purely horizontal local_mean samples
            uint16_t local_mean = (image[offset_input - 1] + image[offset_input] + image[offset_input + 1]) * 16 / 3;

            uint16_t sum_response = 0;
            uint16_t diff_response = 0;
            uint16_t mean = 0;

            int sub_idx;
            for (sub_idx = 0; sub_idx < 4; ++sub_idx) {
                uint8_t a = circular_sample[sub_idx];
                uint8_t b = circular_sample[sub_idx + 4];
                uint8_t c = circular_sample[sub_idx + 8];
                uint8_t d = circular_sample[sub_idx + 12];

                sum_response += abs(a - b + c - d);
                diff_response += abs(a - c) + abs(b - d);
                mean += a + b + c + d;
            }

            response[offset_response] = sum_response - diff_response - abs(mean - local_mean);
        }
}
