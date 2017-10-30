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
 * @param    response  output response image
 * @param    image     input image
 * @param    w         image width
 * @param    h         image height
 */
void ChESS_response_5(      int16_t* restrict response,
                      const uint8_t* restrict image,
                      int w, int h );
