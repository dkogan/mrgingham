#pragma once

namespace mrgingham
{
    struct PointInt
    {
        int x,y;
        PointInt(int _x=0, int _y=0) : x(_x), y(_y) {}
    };

    struct PointDouble
    {
        double x,y;
        PointDouble(double _x=0, double _y=0) : x(_x), y(_y) {}
    };
};
