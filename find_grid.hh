#pragma once

#include <vector>
#include "point.hh"

#define FIND_GRID_SCALE 1000 /* Voronoi diagram is integer-only, so I scale-up
                                to get more resolution */

namespace mrgingham
{
    bool find_grid_from_points( // out
                               std::vector<PointDouble>& points_out,

                               // in
                               const std::vector<Point>& points );
};
