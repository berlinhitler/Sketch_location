#ifndef POINTM_H
#define POINTM_H


class PointM
{
    public:
        double x,y,angle,length;
        PointM();
        PointM(double, double, double, double);
        void setMe(PointM);
        virtual ~PointM();
    protected:
    private:
};

#endif // POINTM_H
