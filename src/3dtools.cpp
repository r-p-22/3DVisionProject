#include "3dtools.h"
//#include "CImg.h"
using namespace cimg_library;

void snake::assign(vector<vector<int> > ps, CImg<unsigned char> &img)
{
    size = ps.size();
    for(int i = 0; i<size; ++i)
    {
    	Eigen::Vector2f v(ps[i][0],ps[i][1]);
        points.push_back(v);
    }
    image = img;
    preProcess();
}

void snake::display(CImg<unsigned char> &img, CImgDisplay &disp, const unsigned char *const colorP, const unsigned char *const color)
{
    for(int i = 0; i<size-1; ++i)
    {
        img.draw_circle(points[i][0], points[i][1], 2, colorP);
        img.draw_line(points[i][0],points[i][1],points[i+1][0],points[i+1][1],color);
    }
    img.draw_circle(points[size-1][0], points[size-1][1], 2, colorP);
    img.draw_line(points[size-1][0],points[size-1][1],points[0][0],points[0][1],color);
}

void snake::setParams(float _gamma, float _alpha, float _beta)
{
    gamma=_gamma; alpha=_alpha; beta=_beta;
}

void snake::preProcess()
{
    CImgList<float> gradient = image.get_gradient();
    CImg<float> norm(image.width(),image.height());
    norm = gradient(0).sqr() + gradient(1).sqr();
    norm.sqrt();
    norm*=-1;
    extGrad = norm.get_gradient();
}

void snake::update()
{
    vector <Eigen::Vector2f> totGrads;
    float norm;
    for(int i=0; i<size; ++i)
    {
    	Eigen::Vector2f p2, p4, grad;
        int i1 = (i-2+size)%size;
        int i2 = (i-1+size)%size;
        int i3 = (i+1+size)%size;
        int i4 = (i+2+size)%size;
        p2 = points[i2] - 2*points[i] + points[i3];
        p4 = points[i1] - 4*points[i2] + 6*points[i] - 4*points[i3] + points[i4];
        grad = -alpha*p2 + beta*p4;
        grad[0]+=extGrad(0,points[i][0],points[i][1]);
        grad[1]+=extGrad(1,points[i][0],points[i][1]);
        totGrads.push_back(grad);
        norm+=grad[0]*grad[0]+grad[1]+grad[1];
    }
    for(int i=0; i<size; ++i)
        points[i] = points[i]-gamma*totGrads[i];

   // cout<<sqrt(norm)<<endl;
}

bool snake::stop() {return end;}




