#ifndef CELL_H
#define CELL_H


class Cell
{
    public:
        double MaxX, MinX, MaxY, MinY, MaxL, MinL, MaxTheta, MinTheta;
        bool isInteresting;
        Cell();
        void setMe(Cell);
        virtual ~Cell();
    protected:
    private:
};

#endif // CELL_H
