#pragma once

namespace mrgingham
{
    struct Point
    {
        int x,y;
        Point(int _x, int _y) : x(_x), y(_y) {}
        Point() {}
    };

    struct PointDouble
    {
        double x,y;
        PointDouble(double _x, double _y) : x(_x), y(_y) {}
        PointDouble() {}
    };
};
