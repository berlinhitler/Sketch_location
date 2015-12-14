#include "SearchSpaceInSHX_SY.h"
#include <cmath>
#include <iostream>

SearchSpaceInSHX_SY::SearchSpaceInSHX_SY()
{
    //ctor
    isInteresting = false;
    maxSHX = 1;
    minSHX = 1;
    maxSY = 1;
    minSY = 1;
    ld = 0.2;//结束区间宽度
}

void SearchSpaceInSHX_SY::setMe(SearchSpaceInSHX_SY source)
{
    maxSHX = source.maxSHX;
    minSHX = source.minSHX;
    maxSY = source.maxSY;
    minSY = source.minSY;
    ld = source.ld;
    isInteresting = source.isInteresting;
}
void SearchSpaceInSHX_SY::Initial()
{
    maxSHX = 0.2;//默认的Shx范围
    minSHX = -0.2;
    maxSY = 1.5;
    minSY = 0.5;
    ld = 0.2;
}
void SearchSpaceInSHX_SY::Divid()
{
    children.clear();
    SearchSpaceInSHX_SY child0;
    child0.maxSHX = (maxSHX + minSHX) / 2;
    child0.minSHX = minSHX;
    child0.maxSY = (maxSY + minSY) / 2;
    child0.minSY = minSY;
    child0.ld = ld;
    children.push_back(child0);
    SearchSpaceInSHX_SY child1;
    child1.maxSHX = maxSHX;
    child1.minSHX = (maxSHX + minSHX) / 2;
    child1.maxSY = (maxSY + minSY) / 2;
    child1.minSY = minSY;
    child1.ld = ld;
    children.push_back(child1);
    SearchSpaceInSHX_SY child2;
    child2.maxSHX = maxSHX;
    child2.minSHX = (maxSHX + minSHX) / 2;
    child2.maxSY = maxSY;
    child2.minSY = (maxSY + minSY) / 2;
    child2.ld = ld;
    children.push_back(child2);
    SearchSpaceInSHX_SY child3;
    child3.maxSHX = (maxSHX + minSHX) / 2;
    child3.minSHX = minSHX;
    child3.maxSY = maxSY;
    child3.minSY = (maxSY + minSY) / 2;
    child3.ld = ld;
    children.push_back(child3);
}
bool SearchSpaceInSHX_SY::isSmallEnough()
{
    double minlSHX = ld;
    double minlSY = ld;
    if (std::abs(maxSHX - minSHX) < minlSHX && std::abs(maxSY - minSY) < minlSY)
    {
        return true;
    }
    return false;
}

SearchSpaceInSHX_SY::~SearchSpaceInSHX_SY()
{
    //dtor
}
