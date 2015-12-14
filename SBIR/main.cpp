#include <iostream>
#include <cvipl/cvipl.h>
#include <cvipl/lem/lem.h>
#include <cvipl/misc.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>
#include <algorithm>
#include <cmath>

#include <cvipl/lem/i3.h>
#include <cvipl/lem/linePt2_nodisp.h>
#include <cvipl/lem/selectPt_merit.h>
#include <cvipl/lem/tif_thin.h>
#include <cvipl/lem/IPSthin.h>
#include <cvipl/lem/tif_IO.h>
#include <cvipl/lem/nD_arrays.h>


#include "Cell.h"
#include "PointM.h"
#include "SearchSpaceInSHX_SY.h"

using namespace cvipl;
using namespace cvipl::lem;
using namespace cv;
using namespace std;
using namespace boost::filesystem;


#pragma region function define
typedef pair<string, double> PAIR;
bool cmp(const PAIR &x, const PAIR &y)
{
    return x.second > y.second;
}
vector<PAIR> ranking;

struct Result
{
    double Sx,MinSHx,MaxSHx,MinSy,MaxSy,r,Tx,Ty,beta,d,Rate;
};
bool sort_cmp(const Result &x, const Result &y)
{
    if (x.Rate != y.Rate)
        return x.Rate > y.Rate;
    return x.d < x.d;
}

void match(path, LEM);
void from_sketch(string fname, LEM *le);
void Search(double, double, double, int,LEM, LEM);
void Search(vector<PointM>, PointM, SearchSpaceInSHX_SY, double, int);
void toPointM(LEM , vector<PointM> *);
void RearrangeByDistance(vector<PointM> *);
void ImageDistanceTransform(vector<PointM>, LEM);
void Cal3DDis(vector<PointM> , vector<PointM> );
void Calculate4ParaMtoO(PointM, vector<PointM> *);
void Calculate4ParaOtoI(PointM, PointM);
void doTransform1(PointM, SearchSpaceInSHX_SY, Cell*);
void doTransform2(PointM, SearchSpaceInSHX_SY, Cell*);
void TransformR(double, Cell &);
void TransformSX(PointM &, double sx);
void TransformSS(PointM , SearchSpaceInSHX_SY, Cell &);
void TransformTx(double , Cell &);
void TransformTy(double , Cell &);
double argMax1(double, double, double);
double argMax2(double, double);
void drawSearchResult(int, double, double, LEM, LEM);
void doTransformDraw(PointM, SearchSpaceInSHX_SY &, Cell & );
int ii = 0;

bool canReturn;
vector<PointM> TDdis;
vector<Result> result;
vector<PointM> sk;
double lowdelta = 0.2, thresholdDis = 8.7, thresholdPer = 0.81;
double sTx, sTy, sRt, iTx, iTy, iRtb, iRt, iSx;
int beta = 1, NowRN = 0, MaxRN = 10;
#pragma endregion function define



int main(int argc, char *argv[]) {
    vector<path> ssfpms;
    LEM sketch;
    Mat s;
    // generate LEM of sketch.
    sketch.from_image_file(argv[1], 20);
    //from_sketch(argv[1],&sketch);
    //LEM sk;
    //sk.from_file(argv[1]);
    sketch.draw(s);
    sketch.to_file("sketch.ssfpm");//for test
    //sk.to_file("sketch2.ssfpm");
    // load all pre-generated ssfpm of queryset.

    try{
        for ( recursive_directory_iterator end, dir(argv[2]); dir != end; ++dir )
        {
            if (!is_directory(dir->path()))
            {
                string name = dir->path().string();
                if (name.find(".ssfpm") != string::npos)
                {
                    ssfpms.push_back(dir->path());
                }
            }
        }
    }
    catch(const boost::filesystem::filesystem_error& e)
    {
        string name = argv[2];
        ssfpms.push_back(name);
    }
    // to compare and match.
    int z = 0;
    for(auto &ss : ssfpms)
    {
        Mat p;
        match(ss,sketch);
        if (ranking.size() <= 20)
        {
            sort(ranking.begin(), ranking.end(), cmp);
        }
        else
        {
            sort(ranking.begin(), ranking.end(), cmp);
            ranking.pop_back();
        }
        z++;
    }
    imwrite("Sketch.jpg",s); //for test
    //sketch.to_file("sketch.ssfpm");
    //toPointM(sketch, new vector<PointM>);

	return 0;
}


void match(path picture, LEM sketch)
{
    LEM pic;
    double score = 0;
    string name = picture.stem().string();
    pic.from_image_file(picture.string(),60);
    //pic.from_ssfpM_file(picture.string());
    //from_sketch(picture.string(),&pic);
    // matching score calculation.


    vector<PointM> m;
    toPointM(pic,&m);
    cout<<"done1"<<endl;
    ImageDistanceTransform(m,pic);
    cout<<"done2"<<endl;
    thresholdDis = 30;
    thresholdPer = 0.9;
    lowdelta = 0.2;
    Search(0, 1, 1, -1, sketch, pic);
    cout<<"done3 "<<result.size()<<endl;
    for(Result r : result)
    {
        cout<<r.Sx<<"\t"<<r.MinSHx<<"\t"<<r.MaxSHx<<"\t"<<r.MinSy<<"\t"<<r.MaxSy<<"\t"<<r.r<<"\t"<<r.Tx<<"\t"<<r.Ty<<"\t"<<r.beta<<"\t"<<r.d<<"\t"<<r.Rate<<endl;
    }
    drawSearchResult(ii,0,1,pic,sketch);
    cout<<"done4"<<endl;
    ranking.push_back(PAIR(name,score));

    Mat p;
    pic.to_file("sketch2.ssfpm");
    pic.draw(p);
    imwrite("Pic.jpg",p);
}


void Search(double pl, double ph, double mpt, int betaIn, LEM sketch, LEM image)
{

    canReturn = false;
    vector<PointM> sls;
    vector<PointM> ils;
    toPointM(sketch, &sls);
    toPointM(image, &ils);
    SearchSpaceInSHX_SY ssss;
    ssss.Initial();
    ssss.ld = lowdelta;
    Calculate4ParaMtoO(sls[0],&sls);
    //cout<<ils.size() * pl<<endl;
    for (int i = int(ils.size() * pl); i < ils.size() * ph; i++)
    {
        //cout<<"test"<<endl;
        NowRN = 0;
        PointM slpm = sls[0];
        PointM ilpm = ils[i];
        //cout<<i<<" "<<ilpm.x<<" "<<ilpm.y<<" "<<ilpm.length<<" "<<ilpm.angle<<endl;
        //cout<<i<<" "<<slpm.x<<" "<<slpm.y<<" "<<slpm.length<<" "<<slpm.angle<<endl;
        Calculate4ParaOtoI(slpm, ilpm);
        //cout<<iTx<<" "<<iTy<<" "<<iRtb<<" "<<iRt<<" "<<iSx<<endl;
        ssss.maxSY = iSx * 1;
        ssss.minSY = iSx / 3;
        //cout<<ssss.maxSY<<" "<<ssss.minSY<<endl;
        {
            ssss.isInteresting = true;
            ssss.Divid();
            for (int j = 0; j < 4; j++)
            {
                Search(sls,ilpm,ssss.children[j],mpt,betaIn);
            }
        }
        if (canReturn) return;
    }
}

void Search(vector<PointM> sls, PointM ilpm, SearchSpaceInSHX_SY ssss, double mpt, int betaIn)
{
    if (NowRN >= MaxRN) return;
    if (canReturn) return;
    double iSyT = (ssss.minSY + ssss.maxSY) / 2;
    int iNumber = 0;
    PointM slpm = sls[0];
    Calculate4ParaOtoI(slpm,ilpm);
    //cout<<iSx<<endl;
    if (iSx >= 2.5 || iSx <= 0.3)
    return;
    double ahd = 0;
    //cout<<"sls size: "<<sls.size()<<endl;
    for (int j = 0; j < sls.size() * mpt; j++)
    {
        slpm = sls[j];
        //cout << slpm.x << " " << slpm.y << " " << slpm.angle << " " << slpm.length <<endl;
        Cell cell;
        if (betaIn == 1)
        {
            doTransform1(slpm,ssss,&cell);
        }
        else
        {
            doTransform2(slpm,ssss,&cell);
        }
        int minT = (int)(cell.MinTheta * 180 / M_PI);
        int maxT = (int)(cell.MaxTheta * 180 / M_PI) + 1;
        double hd = 10000;
        //cout << cell.MinX << " " << cell.MaxX << " " << cell.MinY << " " << cell.MaxY << " " << minT << " " << maxT<<endl;
        int cou = 0;
        for (PointM p : TDdis)
        {
            if (p.x >= cell.MinX && p.x <= cell.MaxX && p.y >= cell.MinY && p.y <= cell.MaxY && p.angle >= minT && p.angle <= maxT)
            {
                if (p.length < hd)
                {
                    hd = p.length;
                }
                cou++;
            }
        }
        //cout<<"Hitted "<<cou <<" hd " <<hd<<endl;
        if (hd <= thresholdDis)
        {
            ahd += hd;
            cell.isInteresting = true;
            iNumber++;
        }
    }
    double iRate = (double) iNumber/(sls.size() * mpt);
    //cout<<iNumber<<" "<<sls.size() * mpt<<" iRate: "<<iRate<<endl;
    ahd = ahd / sls.size();
    //cout<<"ahd: "<<ahd<<endl<<endl;
    if(iRate >= thresholdPer)
    {
        ssss.isInteresting = true;
        if (ssss.isSmallEnough())
        {
            double angleR = iRtb;
            //insert
            Result r;
            r.Sx = iSx;
            r.MinSHx = ssss.minSHX;
            r.MaxSHx = ssss.maxSHX;
            r.MinSy = ssss.minSY;
            r.MaxSy = ssss.maxSY;
            r.r = angleR;
            r.Tx = iTx;
            r.Ty = iTy;
            r.beta = betaIn;
            r.d = ahd;
            r.Rate = iRate;
            result.push_back(r);
            NowRN++;
            return;
        }
        else
        {
            ssss.Divid();
            for (int j = 0; j < 4; j++)
                Search(sls,ilpm,ssss.children[j],mpt,betaIn);
        }
    }
    else
        return;
}

void drawSearchResult(int ii, double pl, double ph, LEM image, LEM sketch)
{
    SearchSpaceInSHX_SY ssss;
    Cell cell;
    sort(result.begin(), result.end(), sort_cmp);
    if (result.size() == 0)
    {
        return;
    }

    ssss.minSHX = result[ii].MinSHx;
    ssss.maxSHX = result[ii].MaxSHx;
    ssss.minSY = result[ii].MinSy;
    ssss.maxSY = result[ii].MaxSy;
    iSx = result[ii].Sx;
    iRt = result[ii].r;
    iTx = result[ii].Tx;
    iTy = result[ii].Ty;
    beta = result[ii].beta;
    vector<PointM> sls;
    toPointM(sketch, &sls);
    double dx,dy;
    int x1 = 0, y1 = 0, x2 = 0, y2 = 0;
    Mat pic = cv::Mat(image.img_h,image.img_w,CV_8UC3,CV_RGB(255,255,255));
    for (int i = 0; i < sls.size(); i++)
    {
        Calculate4ParaMtoO(sls[0],&sls);
        cell.MinX = sls[i].x;
        cell.MaxX = sls[i].x;
        cell.MinY = sls[i].y;
        cell.MaxY = sls[i].y;
        cell.MinL = sls[i].length;
        cell.MaxL = sls[i].length;
        cell.MinTheta = sls[i].angle;
        cell.MaxTheta = sls[i].angle;

        doTransformDraw(sls[i],ssss,cell);

        PointM pm;
        pm.x = (cell.MaxX + cell.MinX) / 2;
        pm.y = (cell.MaxY + cell.MinY) / 2;
        pm.angle = (cell.MaxTheta + cell.MinTheta) / 2;
        pm.length = (cell.MaxL + cell.MinL) / 2;
        dx = abs(pm.length * sin(pm.angle) / 2);
        dy = abs(pm.length * cos(pm.angle) / 2);
        if (pm.angle <= M_PI / 2)
        {
            x1 = (int)(pm.x + dx);
            y1 = (int)(pm.y + dy);
            x2 = (int)(pm.x - dx);
            y2 = (int)(pm.y - dy);
        }
        else
        {
            x1 = (int)(pm.x - dx);
            y1 = (int)(pm.y + dy);
            x2 = (int)(pm.x + dx);
            y2 = (int)(pm.y - dy);
        }
        //cout<<x1<<" "<<y1<<" "<<x2<<" "<<y2<<endl;
        line(pic, cv::Point(x1,y1),cv::Point(x2,y2), cv::Scalar(0,0,0));

    }
    imwrite("result.jpg",pic);
}

void doTransformDraw(PointM slpm, SearchSpaceInSHX_SY &ssss, Cell &cell)
{
    sTx = 0;
    sTy = 0;
    TransformR(sRt * beta, cell);
    PointM slpmb;
    slpmb.setMe(slpm);
    slpmb.x = cell.MinX;
    slpmb.y = cell.MinY;
    slpmb.angle = cell.MinTheta;
    TransformSX(slpmb, iSx);
    cell.MinTheta = slpmb.angle;
    cell.MaxTheta = slpmb.angle;
    cell.MinL = slpmb.length;
    cell.MaxL = slpmb.length;
    cell.MaxX = slpmb.x;
    cell.MinX = slpmb.x;
    cell.MaxY = slpmb.y;
    cell.MinY = slpmb.y;
    TransformSS(slpmb, ssss, cell);
    TransformR(iRt, cell);
    TransformTx(iTx, cell);
    TransformTy(iTy, cell);
    if (cell.MinTheta > cell.MaxTheta)
    {
        double tb = cell.MinTheta;
        cell.MinTheta = cell.MaxTheta;
        cell.MaxTheta = tb;
    }
    cell.MinX = (int)cell.MinX;
    cell.MaxX = (int)cell.MaxX + 1;
    cell.MinL = (int)cell.MinL;
    cell.MaxL = (int)cell.MaxL + 1;
}

#pragma region Search_Transformation
void doTransform1(PointM slpm, SearchSpaceInSHX_SY ssss, Cell * cell)
{
    //cell赋初值
    cell->MinX = slpm.x;
    cell->MaxX = slpm.x;
    cell->MinY = slpm.y;
    cell->MaxY = slpm.y;
    cell->MinL = slpm.length;
    cell->MaxL = slpm.length;
    cell->MinTheta = slpm.angle;
    cell->MaxTheta = slpm.angle;

    TransformR(sRt, *cell);
    PointM slpmb;
    slpmb.setMe(slpm);
    slpmb.x = cell->MinX;
    slpmb.y = cell->MinY;
    slpmb.angle = cell->MinTheta;
    TransformSX(slpmb, iSx);
    TransformSS(slpmb, ssss, *cell);
    TransformR(iRt, *cell);
    TransformTx(iTx, *cell);
    TransformTy(iTy, *cell);

    //规整变换之后的cell
    if (cell->MinTheta < 0) cell->MinTheta += M_PI;
    if (cell->MaxTheta < 0) cell->MaxTheta += M_PI;
    if (cell->MinTheta > cell->MaxTheta)
    {
        double tb = cell->MinTheta;
        cell->MinTheta = cell->MaxTheta;
        cell->MaxTheta = tb;
    }
    if (cell->MinL > cell->MaxL)
    {
        double tb = cell->MinL;
        cell->MinL = cell->MaxL;
        cell->MaxL = tb;
    }
    if (cell->MinX > cell->MaxX)
    {
        double tb = cell->MinX;
        cell->MinX = cell->MaxX;
        cell->MaxX = tb;
    }
    if (cell->MinY > cell->MaxY)
    {
        double tb = cell->MinY;
        cell->MinY = cell->MaxY;
        cell->MaxY = tb;
    }
}

void doTransform2(PointM slpm, SearchSpaceInSHX_SY ssss, Cell * cell)
{
    //cell赋初值
    cell->MinX = slpm.x;
    cell->MaxX = slpm.x;
    cell->MinY = slpm.y;
    cell->MaxY = slpm.y;
    cell->MinL = slpm.length;
    cell->MaxL = slpm.length;
    cell->MinTheta = slpm.angle;
    cell->MaxTheta = slpm.angle;
    //cout<<cell->MinX<<" "<<cell->MaxX<<" "<<cell->MinY<<" "<<cell->MaxY<<" "<<cell->MinL<<" "<<cell->MaxL<<" "<<cell->MinTheta<<" "<<cell->MaxTheta<<endl;
    //cout<<-mRt<<endl;
    TransformR(-sRt, *cell);
    //cout<<cell->MinX<<" "<<cell->MaxX<<" "<<cell->MinY<<" "<<cell->MaxY<<" "<<cell->MinL<<" "<<cell->MaxL<<" "<<cell->MinTheta<<" "<<cell->MaxTheta<<endl;
    //cout<<-mRt<<endl;
    PointM slpmb;
    slpmb.setMe(slpm);
    slpmb.x = cell->MinX;
    slpmb.y = cell->MinY;
    slpmb.angle = cell->MinTheta;
    //cout<<slpmb.x<<" "<<slpmb.y<<" "<<slpmb.angle<<" "<<slpmb.length<<" "<<endl;
    //cout<<slpmb.length<<endl;
    TransformSX(slpmb, iSx);
    //cout<<slpmb.x<<" "<<slpmb.y<<" "<<slpmb.angle<<" "<<slpmb.length<<" "<<endl;
    //cout<<iSx<<endl;
    //cout<<cell->MinX<<" "<<cell->MaxX<<endl;
    TransformSS(slpmb, ssss, *cell);
    //cout<<cell->MinX<<" "<<cell->MaxX<<" "<<cell->MinY<<" "<<cell->MaxY<<" "<<cell->MinL<<" "<<cell->MaxL<<" "<<cell->MinTheta<<" "<<cell->MaxTheta<<endl;
    //cout<<cell->MinX<<" "<<cell->MaxX<<endl;
    TransformR(iRt, *cell);
    //cout<<cell->MinX<<" "<<cell->MaxX<<" "<<cell->MinY<<" "<<cell->MaxY<<" "<<cell->MinL<<" "<<cell->MaxL<<" "<<cell->MinTheta<<" "<<cell->MaxTheta<<endl;
    TransformTx(iTx, *cell);
    //cout<<cell->MinX<<" "<<cell->MaxX<<" "<<cell->MinY<<" "<<cell->MaxY<<" "<<cell->MinL<<" "<<cell->MaxL<<" "<<cell->MinTheta<<" "<<cell->MaxTheta<<endl;
    TransformTy(iTy, *cell);
    //cout<<cell->MinX<<" "<<cell->MaxX<<" "<<cell->MinY<<" "<<cell->MaxY<<" "<<cell->MinL<<" "<<cell->MaxL<<" "<<cell->MinTheta<<" "<<cell->MaxTheta<<endl;

    //cout<<"Transformed "<<cell->MinX<<" "<<cell->MaxX<<" "<<cell->MinY<<" "<<cell->MaxY<<" "<<cell->MinL<<" "<<cell->MaxL<<" "<<cell->MinTheta<<" "<<cell->MaxTheta<<endl;

    //规整变换之后的cell
    if (cell->MinTheta < 0) cell->MinTheta += M_PI;
    if (cell->MaxTheta < 0) cell->MaxTheta += M_PI;
    if (cell->MinTheta > cell->MaxTheta)
    {
        double tb = cell->MinTheta;
        cell->MinTheta = cell->MaxTheta;
        cell->MaxTheta = tb;
    }
    if (cell->MinL > cell->MaxL)
    {
        double tb = cell->MinL;
        cell->MinL = cell->MaxL;
        cell->MaxL = tb;
    }
    if (cell->MinX > cell->MaxX)
    {
        double tb = cell->MinX;
        cell->MinX = cell->MaxX;
        cell->MaxX = tb;
    }
    if (cell->MinY > cell->MaxY)
    {
        double tb = cell->MinY;
        cell->MinY = cell->MaxY;
        cell->MaxY = tb;
    }
    //cout<<cell->MinX<<" "<<cell->MaxX<<" "<<cell->MinY<<" "<<cell->MaxY<<" "<<cell->MinL<<" "<<cell->MaxL<<" "<<cell->MinTheta<<" "<<cell->MaxTheta<<endl;
}

void Calculate4ParaMtoO(PointM slpm, vector<PointM> * ls)
{
    sTx = -slpm.x;
    sTy = -slpm.y;
    sRt = slpm.angle + M_PI / 2;
    //cout<<slpm.angle<<" "<<sRt<<" "<<sTx<<" "<<sTy<<endl;
    for (int i = 0; i < ls->size(); i++)
    {
        (*ls)[i].x += sTx;
        (*ls)[i].y += sTy;
    }

}

void Calculate4ParaOtoI(PointM slpm, PointM ilpm)
{
    iTx = ilpm.x;
    iTy = ilpm.y;
    iRtb = ilpm.angle;
    iRt = iRtb;
    iSx = ilpm.length/slpm.length;
    //cout<<iTx<<" "<<iTy<<" "<<iRtb<<" "<<iRt<<" "<<iSx<<endl;
    //cout << ilpm.length<<" iSx "<<slpm.length<<endl;
}

void TransformR(double r, Cell &cell)
{
    Cell cellB;
    cellB.setMe(cell);
    do
    {
        if (r >= M_PI * 2)
        {
            r -= M_PI * 2;
        }
        else if (r < 0)
        {
            r += M_PI * 2;
        }
        else
        {
            break;
        }
    } while (true);
    //cout<<"r is: "<<r<<endl;
    if (r >= 0 && r < M_PI / 2)
    {
        cell.MinX = cellB.MinX * cos(r) - cellB.MaxY * sin(r);
        cell.MaxX = cellB.MaxX * cos(r) - cellB.MinY * sin(r);
        cell.MinY = cellB.MinX * sin(r) + cellB.MinY * cos(r);
        cell.MaxY = cellB.MaxX * sin(r) + cellB.MaxY * cos(r);
        cell.MinTheta = cell.MinTheta - r;
        cell.MaxTheta = cell.MaxTheta - r;
    }
    else if (r >= M_PI / 2 && r < M_PI)
    {
        cell.MinX = cellB.MaxX * cos(r) - cellB.MaxY * sin(r);
        cell.MaxX = cellB.MinX * cos(r) - cellB.MinY * sin(r);
        cell.MinY = cellB.MinX * sin(r) + cellB.MaxY * cos(r);
        cell.MaxY = cellB.MaxX * sin(r) + cellB.MinY * cos(r);
        cell.MinTheta = cell.MinTheta - r;
        cell.MaxTheta = cell.MaxTheta - r;
    }
    else if (r >= M_PI && r < M_PI / 2 * 3)
    {
        cell.MinX = cellB.MaxX * cos(r) - cellB.MinY * sin(r);
        cell.MaxX = cellB.MinX * cos(r) - cellB.MaxY * sin(r);
        cell.MinY = cellB.MaxX * sin(r) + cellB.MaxY * cos(r);
        cell.MaxY = cellB.MinX * sin(r) + cellB.MinY * cos(r);
        cell.MinTheta = cell.MinTheta - r + M_PI;
        cell.MaxTheta = cell.MaxTheta - r + M_PI;
    }
    else
    {
        cell.MinX = cellB.MinX * cos(r) - cellB.MinY * sin(r);
        cell.MaxX = cellB.MaxX * cos(r) - cellB.MaxY * sin(r);
        cell.MinY = cellB.MaxX * sin(r) + cellB.MinY * cos(r);
        cell.MaxY = cellB.MinX * sin(r) + cellB.MaxY * cos(r);
        cell.MinTheta = cell.MinTheta - r + M_PI;
        cell.MaxTheta = cell.MaxTheta - r + M_PI;
    }
    //cout<<cell.MinX <<" ------ " << cell.MaxX<<endl;
    if (cell.MaxTheta < 0)
    {
        while (cell.MinTheta < 0)
        {
            cell.MinTheta += M_PI;
        }
        while (cell.MaxTheta < 0)
        {
            cell.MaxTheta += M_PI;
        }
    }
    else if (cell.MinTheta < 0)
    {
        double tb = cell.MinTheta;
        cell.MinTheta = cell.MaxTheta;
        cell.MaxTheta = tb + M_PI;
        if (cell.MinTheta > cell.MaxTheta)
        {
            tb = cell.MinTheta;
            cell.MinTheta = cell.MaxTheta;
            cell.MaxTheta = tb;
        }
    }
}

void TransformSX(PointM &slpm, double sx)
{
    slpm.x = slpm.x * sx;
    double angle = slpm.angle;
    if (slpm.angle > M_PI / 2)
    {
        slpm.angle = slpm.angle - M_PI;
        slpm.angle = atan(tan(slpm.angle) * sx);
        if (slpm.angle < 0) slpm.angle = M_PI + slpm.angle;
    }
    else if (slpm.angle == M_PI / 2)
    {
        slpm.angle = M_PI / 2;
    }
    else
    {
        slpm.angle = atan(tan(slpm.angle) * sx);
        if (slpm.angle < 0) slpm.angle = M_PI + slpm.angle;
    }
    slpm.length = slpm.length * sqrt(sin(angle) * sin(angle) * sx * sx + cos(angle) * cos(angle));
    //cout<<slpm.length<<"  in function "<<sin(angle) * sin(angle) * sx * sx<<" "<< cos(angle) * cos(angle)<<endl;
}


void TransformSS(PointM slpm, SearchSpaceInSHX_SY ssss, Cell &cell)
{
    cell.MinX = slpm.x + slpm.y * min(ssss.minSHX, ssss.maxSHX);
    cell.MaxX = slpm.x + slpm.y * max(ssss.minSHX, ssss.maxSHX);
    //cout<<slpm.x<<" "<<slpm.y<<" "<<cell.MinX<<" "<<cell.MaxX<<endl;
    cell.MinY = slpm.y;
    cell.MaxY = slpm.y;
    cell.MinL = slpm.length;
    cell.MaxL = slpm.length;
    cell.MaxTheta = atan(tan(slpm.angle) + ssss.maxSHX);
    cell.MinTheta = atan(tan(slpm.angle) + ssss.minSHX);
    double cotTheta = 1 / tan(slpm.angle);
    if (slpm.angle <= M_PI / 2)
    {
        if (ssss.minSHX >= 0)
        {
            cell.MinL = slpm.length * sqrt(1 + 2 * sin(slpm.angle) * cos(slpm.angle) * ssss.minSHX + sin(slpm.angle) * sin(slpm.angle) * ssss.minSHX * ssss.minSHX);
            cell.MaxL = slpm.length * sqrt(1 + 2 * sin(slpm.angle) * cos(slpm.angle) * ssss.maxSHX + sin(slpm.angle) * sin(slpm.angle) * ssss.maxSHX * ssss.maxSHX);
        }
        else
        {
            if (ssss.minSHX < -cotTheta & -cotTheta < ssss.maxSHX)
            {
                cell.MinL = slpm.length * sin(slpm.angle);
                double mm = argMax1(ssss.maxSHX, ssss.minSHX, cotTheta);
                cell.MaxL = slpm.length * sqrt(1 + 2 * sin(slpm.angle) * cos(slpm.angle) * mm + sin(slpm.angle) * sin(slpm.angle) * mm * mm);//isisis
            }
            else if (ssss.maxSHX <= -cotTheta)
            {
                cell.MaxL = slpm.length * sqrt(1 + 2 * sin(slpm.angle) * cos(slpm.angle) * ssss.minSHX + sin(slpm.angle) * sin(slpm.angle) * ssss.minSHX * ssss.minSHX);
                cell.MinL = slpm.length * sqrt(1 + 2 * sin(slpm.angle) * cos(slpm.angle) * ssss.maxSHX + sin(slpm.angle) * sin(slpm.angle) * ssss.maxSHX * ssss.maxSHX);
            }
            else if (-cotTheta <= ssss.minSHX)
            {
                cell.MinL = slpm.length * sqrt(1 + 2 * sin(slpm.angle) * cos(slpm.angle) * ssss.minSHX + sin(slpm.angle) * sin(slpm.angle) * ssss.minSHX * ssss.minSHX);
                cell.MaxL = slpm.length * sqrt(1 + 2 * sin(slpm.angle) * cos(slpm.angle) * ssss.maxSHX + sin(slpm.angle) * sin(slpm.angle) * ssss.maxSHX * ssss.maxSHX);
            }
        }
    }
    else
    {
        if (ssss.minSHX >= 0)
        {
            if (ssss.minSHX < -cotTheta & -cotTheta < ssss.maxSHX)
            {
                cell.MinL = slpm.length * sin(slpm.angle);
                double mm = argMax1(ssss.maxSHX, ssss.minSHX, cotTheta);
                cell.MaxL = slpm.length * sqrt(1 + 2 * sin(slpm.angle) * cos(slpm.angle) * mm + sin(slpm.angle) * sin(slpm.angle) * mm * mm);//isisis
            }
            else if (ssss.maxSHX <= -cotTheta)
            {
                cell.MaxL = slpm.length * sqrt(1 + 2 * sin(slpm.angle) * cos(slpm.angle) * ssss.minSHX + sin(slpm.angle) * sin(slpm.angle) * ssss.minSHX * ssss.minSHX);
                cell.MinL = slpm.length * sqrt(1 + 2 * sin(slpm.angle) * cos(slpm.angle) * ssss.maxSHX + sin(slpm.angle) * sin(slpm.angle) * ssss.maxSHX * ssss.maxSHX);
            }
            else if (-cotTheta <= ssss.minSHX)
            {
                cell.MinL = slpm.length * sqrt(1 + 2 * sin(slpm.angle) * cos(slpm.angle) * ssss.minSHX + sin(slpm.angle) * sin(slpm.angle) * ssss.minSHX * ssss.minSHX);
                cell.MaxL = slpm.length * sqrt(1 + 2 * sin(slpm.angle) * cos(slpm.angle) * ssss.maxSHX + sin(slpm.angle) * sin(slpm.angle) * ssss.maxSHX * ssss.maxSHX);
            }
        }
        else
        {
            cell.MaxL = slpm.length * sqrt(1 + 2 * sin(slpm.angle) * cos(slpm.angle) * ssss.minSHX + sin(slpm.angle) * sin(slpm.angle) * ssss.minSHX * ssss.minSHX);
            cell.MinL = slpm.length * sqrt(1 + 2 * sin(slpm.angle) * cos(slpm.angle) * ssss.maxSHX + sin(slpm.angle) * sin(slpm.angle) * ssss.maxSHX * ssss.maxSHX);
        }
    }


    if (slpm.y >= 0)
    {
        cell.MinY = slpm.y * ssss.minSY;
        cell.MaxY = slpm.y * ssss.maxSY;
    }
    else
    {
        cell.MinY = slpm.y * ssss.maxSY;
        cell.MaxY = slpm.y * ssss.minSY;
    }
    if (cell.MaxTheta <= M_PI / 2)
    {
        cell.MinTheta = atan(tan(cell.MinTheta) / ssss.minSY);
        cell.MaxTheta = atan(tan(cell.MaxTheta) / ssss.maxSY);
        cell.MinL = cell.MinL * sqrt(cos(cell.MinTheta) * cos(cell.MinTheta) * ssss.minSY * ssss.minSY + sin(cell.MinTheta) * sin(cell.MinTheta));
        cell.MaxL = cell.MaxL * sqrt(cos(cell.MaxTheta) * cos(cell.MaxTheta) * ssss.maxSY * ssss.maxSY + sin(cell.MaxTheta) * sin(cell.MaxTheta));
    }
    else if (cell.MinTheta < M_PI / 2 & M_PI / 2 < cell.MaxTheta)
    {
        cell.MinTheta = atan(tan(cell.MinTheta) / ssss.minSY);
        cell.MaxTheta = atan(tan(cell.MaxTheta) / ssss.minSY);
        double mm = argMax2(cell.MaxTheta, cell.MinTheta);
        cell.MinL = cell.MinL * sqrt(cos(mm) * cos(mm) * ssss.minSY * ssss.minSY + sin(mm) * sin(mm));//isisis
        cell.MaxL = cell.MaxL * ssss.maxSY;
    }
    else if (M_PI / 2 <= cell.MinTheta)
    {
        cell.MinTheta = atan(tan(cell.MinTheta) / ssss.maxSY);
        cell.MaxTheta = atan(tan(cell.MaxTheta) / ssss.minSY);
        cell.MinL = cell.MinL * sqrt(cos(cell.MaxTheta) * cos(cell.MaxTheta) * ssss.minSY * ssss.minSY + sin(cell.MaxTheta) * sin(cell.MaxTheta));
        cell.MaxL = cell.MaxL * sqrt(cos(cell.MinTheta) * cos(cell.MinTheta) * ssss.maxSY * ssss.maxSY + sin(cell.MinTheta) * sin(cell.MinTheta));
    }
}

void TransformTx(double tx, Cell &cell)
{
    cell.MinX += tx;
    cell.MaxX += tx;
}

void TransformTy(double ty, Cell &cell)
{
    cell.MinY += ty;
    cell.MaxY += ty;
}

double argMax1(double ma, double mi, double cot)
{
    if (abs(ma - (-cot)) > abs(mi - (-cot)))
    {
        return ma;
    }
    else
    {
        return mi;
    }

}

double argMax2(double ma, double mi)
{
    if (ma > mi)
    {
        return ma;
    }
    else
    {
        return mi;
    }
}

#pragma endregion Search_Transformation


#pragma region transformations
void toPointM(LEM le, vector<PointM> *ps)
{
    PointM p;
    for (auto &l : le.lines)
    {
        for (int i = 0; i < l.pts.size()-1; i++)
        {
            //cout<<l.pts[i].x<<" "<<l.pts[i].y<<endl;
            double xMid = (double)(l.pts[i].y + l.pts[i + 1].y) / 2;
            double yMid = (double)(l.pts[i].x + l.pts[i + 1].x) / 2;
            //cout<<xMid<<" "<<yMid<<endl;
            double dx = (double)(l.pts[i].y - l.pts[i + 1].y);
            double dy = (double)(l.pts[i].x - l.pts[i + 1].x);
            double angle = 0;
            if (dy == 0)
                angle = 90;
            else
                angle = 180 * atan(dx / dy) / M_PI;
            if (angle == 180) angle = 0;
            if (angle < 0) angle = 180 + angle;
            double length = sqrt(dx * dx + dy * dy);
            p.x = xMid;
            p.y = yMid;
            p.angle = angle * M_PI / 180;
            p.length = length;
            ps->push_back(p);
            //cout<<l.pts[i].x<<" "<<l.pts[i].y<<endl;
            //cout<<p.x<<" "<<p.y<<" "<<p.angle<<" "<<p.length<<" "<<endl<<endl;
        }
    }
    //cout<<le.lines.size()<<" "<<ps->size()<<endl;
    RearrangeByDistance(ps);
}

void RearrangeByDistance(vector<PointM> *ps)
{
    /*
    for (int i = 0; i < ps->size(); i++)
    {
        PointM p = ps->operator[](i);
        cout<<p.x<<" "<<p.y<<" "<<p.angle<<" "<<p.length<<endl;
    }
    cout<<endl;
    */
    for (int i = 1; i < ps->size(); i++)
    {
        for (int j = i; j > 0; j--)
        {
            PointM pN = ps->operator[](j);
            PointM pP = ps->operator[](j-1);
            if (pN.length > pP.length)
            {
                ps->operator[](j) = pP;
                ps->operator[](j-1) = pN;
            }
            else
            {
                break;
            }
        }
    }
    /*
    for (int i = 0; i < ps->size(); i++)
    {
        PointM p = ps->operator[](i);
        cout<<p.x<<" "<<p.y<<" "<<p.angle<<" "<<p.length<<endl;
    }
    cout<<endl;
    */
}

void ImageDistanceTransform(vector<PointM> pd, LEM l)
{
    vector<PointM> ps;
    vector<PointM> ls;
    for (int i =0; i< l.img_w; i++)
    {
        for (int j = 0; j< l.img_h; j++)
        {
            PointM p;
            p.x = i;
            p.y = j;
            ps.push_back(p);
        }
    }
    for(auto &q : pd)
    {
        if (q.angle >= 180 - 0.0001)
            q.angle = 0;
        ls.push_back(q);
        //cout<<p.x<<" "<<p.y<<" "<<p.angle<<" "<<p.length<<endl;
    }
    //cout<<ps.size()<<" "<<ls.size()<<endl;
    //cout<<"done"<<endl;
    Cal3DDis(ps, ls);
}

//need to fix for speed
void Cal3DDis(vector<PointM> ps, vector<PointM> ls)
{
    cout<<180*ps.size()*ls.size()<<endl;
    TDdis.clear();
    int z = 0;
    double minL = 10000;
    for (int k = 0; k< 180; k++)
    {
        for (int i =0; i < ps.size(); i++)
        {
            PointM pd;
            PointM p = ps[i];
            pd.x = (int)p.x;
            pd.y = (int)p.y;
            double dMin = 0;
            for (int j = 0; j < ls.size(); j++)
            {
                PointM l = ls[j];
                double d1 = sqrt((p.x - l.x) * (p.x - l.x) + (p.y - l.y) * (p.y - l.y));
                //cout<<d1<<endl;
                double d2 = abs(l.angle * 180/ M_PI - k);
                //cout<<d2<<endl;
                if (d1 <= 3) d1 = d1 / 3;
                else if (d1 <= 6) d1 = d1 / 2;
                //cout<<d1<<endl;
                if (d2 <= 5) d2 = d2 / 5;
                else if (d2 <= 10) d2 = d2 / 3;
                //cout<<d2<<endl;
                double d = sqrt(d1 * d1 + pow(d2, 4) / 900);
                //cout<<d<<endl;
                if (j == 0 || d < dMin)
                {
                    dMin = d;
                }
                if (round(dMin * 10000.0)/10000.0 == 0) break;
            }
            pd.angle = k;
            pd.length = round(dMin * 10000.0)/10000.0;
            //cout<<"Length "<<pd.length<<endl;
            if (pd.length < minL) minL = pd.length;
            TDdis.push_back(pd);
        }
    }
    //cout<<TDdis.size()<<"  "<< minL<<" "<<z<<endl;
}
#pragma endregion transformations

#pragma region read_from_sketch
void from_sketch(string fname, LEM *le)
{
    Mat image = imread(fname, CV_LOAD_IMAGE_GRAYSCALE);
	if(!image.data) cvipl::error("image read failed");

	const Mat *cv_image = &image;
	if(image.type() == CV_8UC3) {
		cv_image = new Mat();
		cvtColor(image,*cv_image,CV_RGB2GRAY);
	}
	if(cv_image->type() != CV_8UC1) warn_once("image passed to LEM::from_opencv_mat failed to convert to CV_8UC1, results could be unexpected");

	int width = image.cols;
	int height = image.rows;
	le->img_w = width;
	le->img_h = height;
	/*
	printf("height: %d\n", height);
	printf("width: %d\n", width);
	printf("depth: %d\n", image.depth());
	printf("step: %d\n", image.step[0]);
	printf("step: %d\n", image.step[1]);
	printf("channels: %d\n", image.channels());
	*/

	// find the start of the short (not including path) filename
	/*
	int namelen = strlen(img_name);
	int short_name_start = 0;
	for(int i = 0; i < namelen; i++) {
		if(img_name[i] == '/') short_name_start = i+1;
	}

	char out_name[200];
	sprintf(out_name, "%s/%s\0", output_dir, img_name + short_name_start);
	namelen = strlen(out_name);
	*/
	unsigned char **pic = alloc_2D_array<unsigned char>(height, width);


	// copy opencv image to buffTiff
	unsigned char **bufTiff = alloc_2D_array<unsigned char>(height, width);
	for(int r = 0; r < height; r++) {
		memcpy(bufTiff[r], cv_image->data + (r * cv_image->step[0]), sizeof(unsigned char) * width);
	}

    // create inverse grayscale image for input sketch
    Mat img = image.clone();
	for(int r = 0; r < height; r++) {
		memcpy(img.data + (r * image.step[0]), bufTiff[r], sizeof(unsigned char) * width);
	}
    //imwrite("0.jpg",img);


    Mat p;
    bitwise_not(img,p);
    //imwrite("1.jpg",p);
    for(int r = 0; r < height; r++) {
		memcpy(pic[r], p.data + (r * p.step[0]), sizeof(unsigned char) * width);
	}



	//Ctif_thin thin;
	//thin.Generate_thin(bufTiff, width, height, 20, pic);	// 20


	//- display image here
	/*
	Mat img_thin = image.clone();
	for(int r = 0; r < height; r++) {
		memcpy(img_thin.data + (r * image.step[0]), pic[r], sizeof(unsigned char) * width);
	}
	imwrite("thin.jpg", img_thin);
	waitKey(0);
	*/
	//



	vector< LEM_Polyline_type > lin_pts;
	CTiffLin lin;
	lin.Generate_Lin(pic, width, height, &lin_pts);
	/*
	//- display image here
	img_thin = image.clone();
	for(int r = 0; r < height; r++) {
		memcpy(img_thin.data + (r * image.step[0]), pic[r], sizeof(unsigned char) * width);
	}
	// make all LIN pts grey
	img_thin = Scalar(0);
	for(int i = 0; i < lin_pts.size(); i++) {
		for(int j = 0; j < lin_pts[i].pts.size(); j++) {
			LEM_Point pt = lin_pts[i].pts[j];
			img_thin.at<unsigned char>(pt.y, pt.x) = 255;
			cout << pt << endl;
		}
		printf("LIN PTS %d\n", i);
		if(lin_pts[i].polygon) printf("POLYGON!\n");
		else printf("NOT POLY!\n");
		cout << endl << endl;
		//imshow("lin", img_thin);
		//waitKey(0);
		if(lin_pts[i].polygon) continue;
		for(int j = 0; j < lin_pts[i].pts.size(); j++) {
			LEM_Point pt = lin_pts[i].pts[j];
			img_thin.at<unsigned char>(pt.y, pt.x) = 100;
			if(j == 0 || j == lin_pts[i].pts.size()-1) img_thin.at<unsigned char>(pt.y, pt.x) = 255;
		}
	}
	imshow("lin", img_thin);
	waitKey(0);
	// */

	vector< LEM_Polyline > pts_lines;
	CI3::Generate_FPM(&lin_pts, &pts_lines);

	CSelectPtM::GenerateSSFPM(&pts_lines, le);

	free_2D_array(pic);
	free_2D_array(bufTiff);

	// remove last (repeated) point from LEM_Polylines
	for(int l = 0; (unsigned)l < le->lines.size(); l++) {
		le->lines[l].pts.pop_back();
		assert(le->lines[l].con.back() == -100);
		le->lines[l].con.pop_back();
	}
}
#pragma endregion read_from_sketch
