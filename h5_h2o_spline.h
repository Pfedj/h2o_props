#ifndef _H5_H2O_SPLINE_H
#define _H5_H2O_SPLINE_H

#include <string>
#include <iostream>
#include <exception>
#include <matrix3_2/inc/cspline.h>
#include "h5_h2o_data.h"
#include "h5_h2o_alg.h"

namespace H5_H2O
{

enum{
    RHO_I=0,
    TMPR_I,
    NVARS
};

template <class Tval>
class Th2o_1d_spline: public spline3<Tval,spline_grid_v <Tval> >
{
public:
    std::string inp_tag;
    std::string out_tag;
    Th2o_1d_spline(const Tval* x, const Tval* y, int n, const char* ot, const char* it="Pressure"):
        spline3<Tval,spline_grid_v <Tval> >(x,y,n),
        inp_tag(it),
        out_tag(ot)
    {
    }
};

template<class Tval>
class Th2o_liquid_1d_splines_set
{
public:
    Th2o_1d_spline<Tval> Density;
    Th2o_1d_spline<Tval> Enthalpy;
    Th2o_1d_spline<Tval> Temperature;
    Th2o_liquid_1d_splines_set(Th2o_data<Tval>& d):
        Density(d.liquid.Saturated.Pressure.data,d.liquid.Saturated.Density.data,d.liquid.Saturated.Density.nx,"LiquidSaturatedDensity"),
        Enthalpy(d.liquid.Saturated.Pressure.data,d.liquid.Saturated.Enthalpy.data,d.liquid.Saturated.Enthalpy.nx,"LiquidSaturatedEnthalpy"),
        Temperature(d.liquid.Saturated.Pressure.data,d.liquid.Saturated.Temperature.data,d.liquid.Saturated.Temperature.nx,"H2OSaturatedTemperature")
    {
    }
};

template<class Tval>
class Th2o_gas_1d_splines_set
{
public:
    Th2o_1d_spline<Tval> Density;
    Th2o_1d_spline<Tval> Enthalpy;
    Th2o_gas_1d_splines_set(Th2o_data<Tval>& d):
        Density(d.liquid.Saturated.Pressure.data,
                d.gas.Saturated.Density.data,
                d.gas.Saturated.Density.nx,
                "GasSaturatedDensity",""),
        Enthalpy(d.liquid.Saturated.Pressure.data,
                 d.gas.Saturated.Enthalpy.data,
                 d.gas.Saturated.Enthalpy.nx,
                 "GasSaturatedEnthalpy","")
    {
    }
};

template<class Tval, int Nvars>
class Th2o_2d_spline
{
    Tval *data[Nvars];
    Tval *datax[Nvars];
    Tval *datay[Nvars];
    int nx,ny;
    Tval x1,x2,y1,y2;
    Tval dx,dy,idx,idy,idx2,idy2,idxy;
    void check_data(Th2o_2d_table<Tval> *d[Nvars])
    {
        int i,j;
        for(i=0;i<Nvars;i++)
        {
            for(j=0;j<Nvars;j++)
            {
                if(d[i]->nx != d[j]->nx || d[i]->ny != d[j]->ny)
                {
                    throw std::exception("Invalid data set input when creating vector 2d spline\n");
                }
            }
        }
    }
public:
    std::string out_tag[Nvars];
    std::string inp_tag_x;
    std::string inp_tag_y;
    Th2o_2d_spline()
    {
        int i;
        for(i=0;i<Nvars;i++)
        {
            data[i]=datax[i]=datay[i]=0;
        }
    }
    Th2o_2d_spline(Th2o_2d_table<Tval> *d[Nvars], std::string ot[Nvars], const char* itx="Phi", const char* ity="Eta")
    {
        init(d,ot,itx,ity);
    }
    void init(Th2o_2d_table<Tval> *d[Nvars], std::string ot[Nvars], const char* itx="Phi", const char* ity="Eta")
    {
        try
        {
            check_data(d);
        }
        catch(std::exception x)
        {
            std::cerr<<x.what();
            exit(-1);
        }
        x1=d[0]->x1;
        x2=d[0]->x2;
        y1=d[0]->y1;
        y2=d[0]->y2;
        nx=d[0]->nx;
        ny=d[0]->ny;
        dx=(x2-x1)/(nx-1);
        dy=(y2-y1)/(ny-1);
        idx=1./dx;
        idy=1./dy;
		idx2=idx*idx;
		idy2=idy*idy;
		idxy=idx*idy;
        inp_tag_x=itx;
        inp_tag_y=ity;
        int i;
        int siz=nx*ny;
        int siz1=siz*sizeof(Tval);
        for(i=0;i<Nvars;++i)
        {
            data[i]=new Tval[siz];
            datax[i]=new Tval[siz];
            datay[i]=new Tval[siz];
            memcpy(data[i],d[i]->data,siz1);
            memcpy(datax[i],d[i]->datax,siz1);
            memcpy(datay[i],d[i]->datay,siz1);
            out_tag[i]=ot[i];
        }
    }
    ~Th2o_2d_spline()
    {
        int i;
        for(i=0;i<Nvars;++i)
        {
            if(data [i]!=0)delete[]data [i];
            if(datax[i]!=0)delete[]datax[i];
            if(datay[i]!=0)delete[]datay[i];
        }
    }

    inline void operator()(Tval x, Tval y, Tf_xy<Tval> res[Nvars])const
    {
        Tval x_loc,y_loc;
        int i,j;
        if(x>=x2)
        {
            x=x2;
            i=nx-2;
            x_loc=1.;
        }
        else if(x<=x1)
        {
            x=x1;
            i=0;
            x_loc=0;
        }
        else
        {
            i=int((x-x1)*idx);
            x_loc=(x-(x1+dx*i))*idx;
        }
        if(y>=y2)
        {
            y=y2;
            j=ny-2;
            y_loc=1.;
        }
        else if(y<=y1)
        {
            y=y1;
            j=0;
            y_loc=0.;
        }
        else
        {
            j=int((y-y1)*idy);
            y_loc=(y-(y1+dy*j))*idy;
        }
        int k[4];           
        /*     2 3
               0 1    */
        k[0]=i+j*nx;
        k[1]=k[0]+1;
        k[2]=k[0]+nx;
        k[3]=k[2]+1;
        Tval o[40];
        /*      c01 c11
                c00 c10         */
        Tval c00,c01,c10,c11;
        Tval cx00,cx01,cx10,cx11;
        Tval cy00,cy01,cy10,cy11;
        int l;
        
o[0] = 3.*x_loc;
o[1] = x_loc*x_loc;
o[2] = -2.*o[1];
o[3] = o[0] + o[2];
o[4] = -1.*x_loc*o[3];
o[5] = 1. + o[4];
o[6] = 3.*y_loc;
o[7] = y_loc*y_loc;
o[8] = -2.*o[7];
o[9] = o[6] + o[8];
o[10] = -1.*y_loc*o[9];
o[11] = 1. + o[10];
o[12] = x_loc*o[1];
o[13] = x_loc + o[2] + o[12];
o[14] = -1.*o[1];
o[15] = o[12] + o[14];
o[16] = y_loc*o[7];
o[17] = y_loc + o[8] + o[16];
o[18] = -1.*o[7];
o[19] = o[16] + o[18];
o[20] = -3.*x_loc;
o[21] = -4.*x_loc;
o[22] = 3. + o[21];
o[23] = -1.*x_loc*o[22];
o[24] = 2.*o[1];
o[25] = o[20] + o[23] + o[24];
o[26] = 3.*o[1];
o[27] = 1. + o[21] + o[26];
o[28] = -2.*x_loc;
o[29] = o[26] + o[28];
o[30] = -3.*y_loc;
o[31] = -4.*y_loc;
o[32] = 3. + o[31];
o[33] = -1.*y_loc*o[32];
o[34] = 2.*o[7];
o[35] = o[30] + o[33] + o[34];
o[36] = 3.*o[7];
o[37] = 1. + o[31] + o[36];
o[38] = -2.*y_loc;
o[39] = o[36] + o[38];

        for(l=0;l<Nvars;l++)
        {
			res[l].x=x;
			res[l].y=y;

c00 =data [l][k[0]];
c10 =data [l][k[1]];
c01 =data [l][k[2]];
c11 =data [l][k[3]];

cx00=datax[l][k[0]]*dx;
cx10=datax[l][k[1]]*dx;
cx01=datax[l][k[2]]*dx;
cx11=datax[l][k[3]]*dx;

cy00=datay[l][k[0]]*dy;
cy10=datay[l][k[1]]*dy;
cy01=datay[l][k[2]]*dy;
cy11=datay[l][k[3]]*dy;

res[l].f = x_loc*y_loc*c11*o[3]*o[9] + y_loc*c01*o[5]*o[9] +\
 x_loc*c10*o[3]*o[11] + c00*o[5]*o[11] + y_loc*cx01*o[9]*o[\
13] + cx00*o[11]*o[13] + y_loc*cx11*o[9]*o[15] + cx10*o[11]\
*o[15] + x_loc*cy10*o[3]*o[17] + cy00*o[5]*o[17] + x_loc*cy\
11*o[3]*o[19] + cy01*o[5]*o[19];

res[l].dfdx = y_loc*c11*o[3]*o[9] + c10*o[3]*o[11] + cy10*o[\
3]*o[17] + cy11*o[3]*o[19] + x_loc*y_loc*c11*o[9]*o[22] + x\
_loc*c10*o[11]*o[22] + x_loc*cy10*o[17]*o[22] + x_loc*cy11*\
o[19]*o[22] + y_loc*c01*o[9]*o[25] + c00*o[11]*o[25] + cy00\
*o[17]*o[25] + cy01*o[19]*o[25] + y_loc*cx01*o[9]*o[27] + c\
x00*o[11]*o[27] + y_loc*cx11*o[9]*o[29] + cx10*o[11]*o[29];

res[l].dfdy = x_loc*c11*o[3]*o[9] + c01*o[5]*o[9] + cx01*o[9\
]*o[13] + cx11*o[9]*o[15] + x_loc*y_loc*c11*o[3]*o[32] + y_\
loc*c01*o[5]*o[32] + y_loc*cx01*o[13]*o[32] + y_loc*cx11*o[\
15]*o[32] + x_loc*c10*o[3]*o[35] + c00*o[5]*o[35] + cx00*o[\
13]*o[35] + cx10*o[15]*o[35] + x_loc*cy10*o[3]*o[37] + cy00\
*o[5]*o[37] + x_loc*cy11*o[3]*o[39] + cy01*o[5]*o[39];

res[l].dfdx*=idx;
res[l].dfdy*=idy;

        }   
    }
    inline void d2(Tval x, Tval y, Tf_xy<Tval> res[Nvars])const
    {
        Tval x_loc,y_loc;
        int i,j;
        if(x>=x2)
        {
            x=x2;
            i=nx-2;
            x_loc=1.;
        }
        else if(x<=x1)
        {
            x=x1;
            i=0;
            x_loc=0;
        }
        else
        {
            i=int((x-x1)*idx);
            x_loc=(x-(x1+dx*i))*idx;
        }
        if(y>=y2)
        {
            y=y2;
            j=ny-2;
            y_loc=1.;
        }
        else if(y<=y1)
        {
            y=y1;
            j=0;
            y_loc=0.;
        }
        else
        {
            j=int((y-y1)*idy);
            y_loc=(y-(y1+dy*j))*idy;
        }
        int k[4];           
        /*     2 3
               0 1    */
        k[0]=i+j*nx;
        k[1]=k[0]+1;
        k[2]=k[0]+nx;
        k[3]=k[2]+1;
        Tval o[40];
        /*      c01 c11
                c00 c10         */
        Tval c00,c01,c10,c11;
        Tval cx00,cx01,cx10,cx11;
        Tval cy00,cy01,cy10,cy11;
        int l;
        

o[0] = 2.*x_loc;
o[1] = -1. + o[0];
o[2] = -1. + y_loc;
o[3] = pow(o[2],2.);
o[4] = 2.*y_loc;
o[5] = 1. + o[4];
o[6] = -2.*y_loc;
o[7] = 3. + o[6];
o[8] = y_loc*y_loc;
o[9] = 4.*x_loc;
o[10] = -3. + o[9];
o[11] = -3. + o[4];
o[12] = 3.*x_loc;
o[13] = 6.*x_loc;
o[14] = -4.*x_loc;
o[15] = 3. + o[14];
o[16] = -1. + x_loc;
o[17] = pow(o[16],2.);
o[18] = 1. + o[0];
o[19] = -1. + o[4];
o[20] = 4.*y_loc;
o[21] = -3. + o[20];
o[22] = -2.*x_loc;
o[23] = 3. + o[22];
o[24] = x_loc*x_loc;
o[25] = -3. + o[0];
o[26] = -4.*y_loc;
o[27] = 3. + o[26];
o[28] = 3.*y_loc;
o[29] = 6.*y_loc;
o[30] = 3.*o[8];
o[31] = 1. + o[26] + o[30];
o[32] = 3.*o[24];
o[33] = 1. + o[14] + o[32];
        for(l=0;l<Nvars;l++)
        {
			res[l].x=x;
			res[l].y=y;

c00 =data [l][k[0]];
c10 =data [l][k[1]];
c01 =data [l][k[2]];
c11 =data [l][k[3]];

cx00=datax[l][k[0]];
cx10=datax[l][k[1]];
cx01=datax[l][k[2]];
cx11=datax[l][k[3]];

cy00=datay[l][k[0]];
cy10=datay[l][k[1]];
cy01=datay[l][k[2]];
cy11=datay[l][k[3]];


res[l].dfdx2 = -4.*dy*x_loc*y_loc*cy10*o[3] + 6.*dy*y_loc*cy\
00*o[1]*o[3] - 4.*x_loc*c10*o[3]*o[5] + 6.*c00*o[1]*o[3]*o[\
5] - 4.*dy*x_loc*cy11*o[2]*o[8] + 6.*dy*cy01*o[1]*o[2]*o[8]\
 + (-6. + 12.*x_loc)*c01*o[7]*o[8] - 2.*c10*o[3]*o[5]*o[10]\
 + 4.*x_loc*c11*o[8]*o[11] + 2.*c11*o[8]*o[10]*o[11] + 2.*d\
x*cx00*o[3]*o[5]*(-2. + o[12]) + 2.*dx*cx10*o[3]*o[5]*(-1. \
+ o[12]) + dx*cx01*o[7]*o[8]*(-4. + o[13]) + dx*cx11*o[7]*o\
[8]*(-2. + o[13]) + 2.*dy*y_loc*cy10*o[3]*o[15] + 2.*dy*cy1\
1*o[2]*o[8]*o[15];
res[l].dfdy2 = -4.*dx*x_loc*y_loc*cx01*o[17] - 4.*y_loc*c01*\
o[17]*o[18] + 6.*dx*x_loc*cx00*o[17]*o[19] + 6.*c00*o[17]*o\
[18]*o[19] - 2.*c01*o[17]*o[18]*o[21] - 4.*dx*y_loc*cx11*o[\
16]*o[24] + 6.*dx*cx10*o[16]*o[19]*o[24] + (-6. + 12.*y_loc\
)*c10*o[23]*o[24] + 4.*y_loc*c11*o[24]*o[25] + 2.*c11*o[21]\
*o[24]*o[25] + 2.*dx*x_loc*cx01*o[17]*o[27] + 2.*dx*cx11*o[\
16]*o[24]*o[27] + 2.*dy*cy00*o[17]*o[18]*(-2. + o[28]) + 2.\
*dy*cy01*o[17]*o[18]*(-1. + o[28]) + dy*cy10*o[23]*o[24]*(-\
4. + o[29]) + dy*cy11*o[23]*o[24]*(-2. + o[29]);
res[l].dfdxy = 6.*(dy*x_loc*cy00*o[16]*o[31] - 1.*dy*x_loc*c\
y10*o[16]*o[31] + y_loc*(x_loc*(6.*(-1. + x_loc + y_loc - 1\
.*x_loc*y_loc)*c01 - 6.*c10 + 6.*x_loc*c10 + 6.*y_loc*c10 -\
 6.*x_loc*y_loc*c10 + 6.*c11 - 6.*x_loc*c11 - 6.*y_loc*c11 \
+ 6.*x_loc*y_loc*c11 + 2.*dx*cx10 - 3.*dx*x_loc*cx10 - 2.*d\
x*y_loc*cx10 + 3.*dx*x_loc*y_loc*cx10 - 2.*dx*cx11 + 3.*dx*\
x_loc*cx11 + 2.*dx*y_loc*cx11 - 3.*dx*x_loc*y_loc*cx11 + 2.\
*dy*cy01 - 2.*dy*x_loc*cy01 - 3.*dy*y_loc*cy01 + 3.*dy*x_lo\
c*y_loc*cy01 - 2.*dy*cy11 + 2.*dy*x_loc*cy11 + 3.*dy*y_loc*\
cy11 - 3.*dy*x_loc*y_loc*cy11 + 6.*c00*o[2]*o[16]) + dx*cx0\
0*o[2]*o[33] - 1.*dx*cx01*o[2]*o[33]));

res[l].dfdxy*=idxy;
res[l].dfdx2*=idx2;
res[l].dfdy2*=idy2;

		}
    }
};

template<class Tval>
class Th2o_gas_2d_spline: public Th2o_2d_spline<Tval,2>
{
public:
    Th2o_gas_2d_spline(Th2o_gas_data<Tval>& d)
    {
        Th2o_2d_table<Tval> *t[]={&d.Density,&d.Temperature};
        std::string ot[]={std::string("Density"),std::string("Temperature")};
        init(t,ot);
    }
};

template<class Tval>
class Th2o_fluid_2d_spline: public Th2o_2d_spline<Tval,2>
{
public:
    Th2o_fluid_2d_spline(Th2o_fluid_data<Tval>& d)
    {
        Th2o_2d_table<Tval> *t[]={&d.Density,&d.Temperature};
        std::string ot[]={std::string("Density"),std::string("Temperature")};
        init(t,ot);
    }
};

template<class Tval>
class Th2o_liquid_2d_spline: public Th2o_2d_spline<Tval,2>
{
public:
    Th2o_liquid_2d_spline(Th2o_liquid_data<Tval>& d)
    {
        Th2o_2d_table<Tval> *t[]={&d.Density,&d.Temperature};
        std::string ot[]={std::string("Density"),std::string("Temperature")};
        init(t,ot);
    }
};

template<class Tval>
class Th2o_gas_splines: public Th2o_gas_2d_spline<Tval> 
{
public:
    Th2o_gas_1d_splines_set<Tval> Saturated;
    Th2o_gas_splines(Th2o_data<Tval>& d): 
        Th2o_gas_2d_spline<Tval> (d.gas),
        Saturated(d)
    {
    }
};

template<class Tval>
class Th2o_fluid_splines: public Th2o_fluid_2d_spline<Tval> 
{
public:
    Th2o_fluid_splines(Th2o_data<Tval>& d): 
        Th2o_fluid_2d_spline<Tval> (d.fluid)
    {
    }
};

template<class Tval>
class Th2o_liquid_splines: public Th2o_liquid_2d_spline<Tval> 
{
public:
    Th2o_liquid_1d_splines_set<Tval> Saturated;
    Th2o_liquid_splines(Th2o_data<Tval>& d): 
        Th2o_liquid_2d_spline<Tval> (d.liquid),
        Saturated(d)
    {
    }
};

template<class Tval>
class Th2o_splines
{
public:
    Th2o_gas_splines<Tval> gas;
    Th2o_fluid_splines<Tval> fluid;
    Th2o_liquid_splines<Tval> liquid;
    Th2o_data_parameters<Tval> attrs;
    Th2o_splines(Th2o_data<Tval>& d):
        gas(d),
        fluid(d),
        liquid(d),
        attrs(d.attrs)
    {
        hw700=HWsat(attrs.MinimalPressure);
    }
    inline Tval Tsat(Tval P)const
    {
        return attrs.CriticalTemperature*liquid.Saturated.Temperature(P);
    }
    inline Tval dTsat(Tval P)const
    {
        return attrs.CriticalTemperature*liquid.Saturated.Temperature.d(P);
    }
    inline Tval d2Tsat(Tval P)const
    {
        return attrs.CriticalTemperature*liquid.Saturated.Temperature.d2(P);
    }
    inline Tval HWsat(Tval P)const
    {
        return attrs.CriticalEnthalpyW*liquid.Saturated.Enthalpy(P);
    }
    inline Tval dHWsat(Tval P)const
    {
        return attrs.CriticalEnthalpyW*liquid.Saturated.Enthalpy.d(P);
    }
    inline Tval d2HWsat(Tval P)const
    {
        return attrs.CriticalEnthalpyW*liquid.Saturated.Enthalpy.d2(P);
    }
    inline Tval HSsat(Tval P)const
    {
        return attrs.CriticalEnthalpyS*gas.Saturated.Enthalpy(P);
    }
    inline Tval dHSsat(Tval P)const
    {
        return attrs.CriticalEnthalpyS*gas.Saturated.Enthalpy.d(P);
    }
    inline Tval d2HSsat(Tval P)const
    {
        return attrs.CriticalEnthalpyS*gas.Saturated.Enthalpy.d2(P);
    }
    inline Tval RoWsat(Tval P)const
    {
        return attrs.CriticalDensity*liquid.Saturated.Density(P);
    }
    inline Tval dRoWsat(Tval P)const
    {
        return attrs.CriticalDensity*liquid.Saturated.Density.d(P);
    }
    inline Tval d2RoWsat(Tval P)const
    {
        return attrs.CriticalDensity*liquid.Saturated.Density.d2(P);
    }
    inline Tval RoSsat(Tval P)const
    {
        return attrs.CriticalDensity*gas.Saturated.Density(P);
    }
    inline Tval dRoSsat(Tval P)const
    {
        return attrs.CriticalDensity*gas.Saturated.Density.d(P);
    }
    inline Tval d2RoSsat(Tval P)const
    {
        return attrs.CriticalDensity*gas.Saturated.Density.d2(P);
    }
    Tval hw700;
};

}

#endif