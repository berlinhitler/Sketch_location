#include "Cell.h"

Cell::Cell()
{
    //ctor
    MaxX = 0; MinX = 0; MaxY = 0; MinY = 0; MaxTheta = 0; MinTheta = 0;
    isInteresting = false;
}

void Cell::setMe(Cell _cell)
{
    MaxX = _cell.MaxX;
    MinX = _cell.MinX;
    MaxY = _cell.MaxY;
    MinY = _cell.MinY;
    MaxL = _cell.MaxL;
    MinL = _cell.MinL;
    MaxTheta = _cell.MaxTheta;
    MinTheta = _cell.MinTheta;
    isInteresting = _cell.isInteresting;
}

Cell::~Cell()
{
    //dtor
}
