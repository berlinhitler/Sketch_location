#ifndef SEARCHSPACEINSHX_SY_H
#define SEARCHSPACEINSHX_SY_H
#include <vector>

class SearchSpaceInSHX_SY
{
    public:
        std::vector<SearchSpaceInSHX_SY> children;
        bool isInteresting;
        double maxSHX, minSHX, maxSY, minSY, ld;
        SearchSpaceInSHX_SY();
        void setMe(SearchSpaceInSHX_SY);
        void Initial();
        void Divid();
        bool isSmallEnough();
        virtual ~SearchSpaceInSHX_SY();
    protected:
    private:
};

#endif // SEARCHSPACEINSHX_SY_H
