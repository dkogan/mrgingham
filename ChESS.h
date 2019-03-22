#pragma once

/*
  This is the reference implementation from this paper:

  https://arxiv.org/abs/1301.5491

  Dima made a few modifications.

  Obtained from here:

  http://www-sigproc.eng.cam.ac.uk/Main/SB476Chess

  There's a more full-featured GPL-licensed implementation on that page
*/


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
void mrgingham_ChESS_response_5(      int16_t* response,
                                const uint8_t* image,
                                int w, int h,
                                int stride);
