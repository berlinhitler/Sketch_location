#include "PointM.h"

PointM::PointM()
{
    x = 0;
    y = 0;
    angle = 0;
    length = 0;
}

PointM::PointM(double _x, double _y, double _angle, double _length)
{
    x = _x;
    y = _y;
    angle = _angle;
    length = _length;
}

void PointM::setMe(PointM p)
{
    x = p.x;
    y = p.y;
    angle = p.angle;
    length = p.length;
}

PointM::~PointM()
{
    //dtor
}
