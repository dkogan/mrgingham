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

#include <stddef.h>	// size_t
#include <stdint.h>
#include <stdlib.h>	// abs

/**
 * Perform the ChESS corner detection algorithm with a 5 px sampling radius
 *
 * @param	w	image width
 * @param	h	image height
 * @param	image	input image
 * @param	response	output response image
 */
void corner_detect5(const size_t w, const size_t h, const uint8_t image[w * h], int16_t response[w * h])
{
	int x, y;
	// funny bounds due to sampling ring radius (5) and border of previously applied blur (2)
	for (y = 7; y < h - 7; y++)
		for (x = 7; x < w - 7; x++) {
			const unsigned offset = x + y * w;
			uint8_t circular_sample[16];

			circular_sample[2] = image[offset - 2 - 5 * w];
			circular_sample[1] = image[offset - 5 * w];
			circular_sample[0] = image[offset + 2 - 5 * w];
			circular_sample[8] = image[offset - 2 + 5 * w];
			circular_sample[9] = image[offset + 5 * w];
			circular_sample[10] = image[offset + 2 + 5 * w];
			circular_sample[3] = image[offset - 4 - 4 * w];
			circular_sample[15] = image[offset + 4 - 4 * w];
			circular_sample[7] = image[offset - 4 + 4 * w];
			circular_sample[11] = image[offset + 4 + 4 * w];
			circular_sample[4] = image[offset - 5 - 2 * w];
			circular_sample[14] = image[offset + 5 - 2 * w];
			circular_sample[6] = image[offset - 5 + 2 * w];
			circular_sample[12] = image[offset + 5 + 2 * w];
			circular_sample[5] = image[offset - 5];
			circular_sample[13] = image[offset + 5];

			// purely horizontal local_mean samples
			uint16_t local_mean = (image[offset - 1] + image[offset] + image[offset + 1]) * 16 / 3;

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

			response[offset] = sum_response - diff_response - abs(mean - local_mean);
		}
}
