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

# arbitrary units
square_size = 100

def black_cell(x,y, *,
               double_width  = False,
               double_height = False):
    print("2 2 0 0 0 0 50 -1 20 0.000 0 0 -1 0 0 5")

    X0 = (x+0)*square_size
    X1 = (x+1)*square_size
    Y0 = (y+0)*square_size
    Y1 = (y+1)*square_size
    if double_width:
        X1 += square_size
    if double_height:
        Y1 += square_size

    print(f"  {X0} {Y0} {X1} {Y0} {X1} {Y1} {X0} {Y1} {X0} {Y0}")

if (args.Wcorners // 2) * 2 != args.Wcorners or \
   (args.Hcorners // 2) * 2 != args.Hcorners:
    print("Wcorners and Hcorners are assumed to be even",
          file = sys.stderr)
    sys.exit(1)


Wcells = args.Wcorners + 3
Hcells = args.Hcorners + 3

print(header)

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
