#!/usr/bin/python3

r'''Generate a .fig file with a chessboard with a given geometry

SYNOPSIS

  ./generate-chessboard-fig.py 10 14 > board-10-14.fig
  fig2dev -L pdf board-10-14.fig > board-10-14.pdf
'''



import sys
import argparse
import re
import os

def parse_args():

    parser = \
        argparse.ArgumentParser(description = __doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--blobs',
                        action = 'store_true',
                        help='''If given, we generate a grid of circles, not a
                        chessboard''')
    parser.add_argument('--circle-radius-cells',
                        type=float,
                        default = 0.2,
                        help='''Applies only if --blobs. Specifies the radius of
                        each circle, in the units of grid spacing: radius == 0.5
                        would create circles that are large enough such that
                        adjacent ones touch.''')
    parser.add_argument('Wcorners',
                        type=int,
                        help='''The number of chessboard corners in the
                        horizontal direction''')
    parser.add_argument('Hcorners',
                        type=int,
                        help='''The number of chessboard corners in the vertical
                        direction''')

    args = parser.parse_args()
    return args

args = parse_args()


import numpy as np
import numpysane as nps


# The fig format is documented in many places. For instance:
#   https://web.archive.org/web/20070920204655/http://epb.lbl.gov/xfig/fig-format.html
#
# Here I simply took the chessboards I drew in xfig, and extended those .fig
# files using this script. The header and chessboard-square definitions come
# directly from the .fig files I made in xfig

header = fr'''#FIG 3.2  Produced by the mrgingham chessboard generator: {sys.argv[0]}
Portrait
Center
Metric
Letter
100.00
Single
-2
1200 2'''

# Arbitrary units. fig2dev appears to compute the bounding box in integer units,
# so if this is too small, the bounding box will be off
square_size = 10000

def black_cell(x,y, *,
               double_width  = False,
               double_height = False):

    X0 = (x+0)*square_size
    X1 = (x+1)*square_size
    Y0 = (y+0)*square_size
    Y1 = (y+1)*square_size
    if double_width:
        X1 += square_size
    if double_height:
        Y1 += square_size

    print("2 2 0 0 0 0 50 -1 20 0.000 0 0 -1 0 0 5")
    print(f"  {X0} {Y0} {X1} {Y0} {X1} {Y1} {X0} {Y1} {X0} {Y0}")

def black_circle(x,y,R):

    '''Units of all x,y is "cells". Units of R is "fig units"'''

    X = x * square_size
    Y = y * square_size

    print(f"1 3 0 0 0 0 50 -1 20 0 1 0 {X} {Y} {R} {R} {X} {Y} {X+R} {Y}")




Wcells = args.Wcorners + 3
Hcells = args.Hcorners + 3

print(header)

if not args.blobs:

    if (args.Wcorners // 2) * 2 != args.Wcorners or \
       (args.Hcorners // 2) * 2 != args.Hcorners:
        print("Wcorners and Hcorners are assumed to be even",
              file = sys.stderr)
        sys.exit(1)

    # top/bottom edges
    for y in (0,Hcells-2):
        for x in range(2,Wcells-2,2):
            black_cell(x,y,
                       double_height = True)

    # left/right edges
    for x in (0,Wcells-2):
        for y in range(2,Hcells-2,2):
            black_cell(x,y,
                  double_width = True)

    # middle
    for y in range(2,Hcells-2,2):
        for x in range(3,Wcells-2,2):
            black_cell(x,y)
    for y in range(3,Hcells-2,2):
        for x in range(2,Wcells-2,2):
            black_cell(x,y)
else:

    R = int(round(square_size * args.circle_radius_cells))

    for y in range(0,args.Hcorners):
        for x in range(0,args.Wcorners):
            black_circle(x,y,R)
