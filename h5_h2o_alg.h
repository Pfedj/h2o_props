#ifndef _H5_H2O_ALG_H
#define _H5_H2O_ALG_H

namespace H5_H2O

{

template<class Tval>
struct Tf_x
{
    Tval x,f,df,d2f;
    Tf_x():x(0),f(0),df(0),d2f(0){}
};

template<class val_t>
struct Tf_xy
{
    Tf_xy():x(0),y(0),f(0),dfdx(0),dfdy(0),dfdx2(0),dfdy2(0),dfdxy(0){}
    Tf_xy(val_t x_, val_t y_): x(x_),y(y_),f(0),dfdx(0),dfdy(0),dfdx2(0),dfdy2(0),dfdxy(0){}
    val_t x,y; //! ���������
    val_t f;   //! �������� �������
    val_t dfdx,dfdy;    //! �������� �����������
    val_t dfdx2,dfdy2,dfdxy; //! �������� ������ �����������
    /*! �������������� ������ ����������� ��� �������� � ����� ����������� */
    inline void convert_to_ph(const Tf_xy<val_t>& x_ph,const Tf_xy<val_t>& y_ph, val_t p, val_t h)
    {
        x=p;
        y=h;
        val_t dfdx_1 = dfdx*x_ph.dfdx + dfdy*y_ph.dfdx;
        val_t dfdy_1 = dfdx*x_ph.dfdy + dfdy*y_ph.dfdy;
        dfdx=dfdx_1;
        dfdy=dfdy_1;
        dfdx2=dfdy2=dfdxy=0.;
    }
    inline Tf_xy<val_t>& operator*=(const val_t& v)
    {
        f*=    v;
        dfdx*= v;
        dfdy*= v;
        dfdx2*=v;
        dfdy2*=v;
        dfdxy*=v;
        return *this;
    }
    inline Tf_xy<val_t>& operator/=(const val_t& x)
    {
        val_t v=1./x;
        f*=    v;
        dfdx*= v;
        dfdy*= v;
        dfdx2*=v;
        dfdy2*=v;
        dfdxy*=v;
        return *this;
    }
};

template <class val_t>
void Bezier_sq_lu_interpolation(const Tf_xy<val_t>& px,const Tf_xy<val_t>& py,Tf_xy<val_t> &res)
{
    val_t x=py.x;
    val_t x0=px.x;
    val_t y=px.y;
    val_t y0=py.y;
    res.x=x;
    res.y=y;
    val_t o[40];
    
/*-------------OPTIMIZED CODE GENERATED BY EXTERNAL CODEGENERATOR - !!!DO NOT EDIT!!! -------------------------------*/
    
o[0] = -1.*y0;
o[1] = y + o[0];
o[2] = pow(o[1],2.);
o[3] = o[1]*o[2];
o[4] = -1.*x0;
o[5] = -1.*y;
o[6] = x + y0 + o[4] + o[5];
o[7] = pow(o[6],2.);
o[8] = o[6]*o[7];
o[9] = 1/o[8];
o[10] = -2.*px.f;
o[11] = 2.*py.f;
o[12] = px.dfdx + py.dfdy;
o[13] = -1.*o[6]*o[12];
o[14] = o[10] + o[11] + o[13];
o[15] = 1/o[7];
o[16] = -3.*px.f;
o[17] = 3.*py.f;
o[18] = 2.*py.dfdy;
o[19] = px.dfdx + o[18];
o[20] = -1.*o[6]*o[19];
o[21] = o[16] + o[17] + o[20];
o[22] = -6.*px.f;
o[23] = 6.*py.f;
o[24] = pow(o[7],2.);
o[25] = 1/o[24];
o[26] = 2.*px.dfdx;
o[27] = 6.*px.f;
o[28] = -6.*py.f;
o[29] = px.dfdxy*x;
o[30] = -1.*px.dfdxy*x0;
o[31] = -1.*px.dfdxy*y;
o[32] = px.dfdxy*y0;
res.f = py.f + py.dfdy*o[1] - 1.*o[3]*o[9]*o[14] - 1.*o[2]*o[15]*o[21];
res.dfdx = py.dfdx + py.dfdxy*o[1] + o[2]*o[9]*(-1.*(px.dfdx + 3.*py.dfdx + 2.*(py.dfdy + py.dfdxy*(-1.*x + x0 + y + o[0\
])))*o[6] + o[22] + o[23]) + o[3]*o[25]*(o[22] + o[23] - 1.*o[6]*(2.*py.dfdx - 1.*py.dfdxy*x + py.dfdxy*x0 + py.dfdxy*y\
- 1.*py.dfdxy*y0 + o[18] + o[26]));
res.dfdy = py.dfdy - 3.*o[2]*o[9]*o[14] - 2.*o[1]*o[15]*o[21] + o[2]*o[9]*(o[27] + o[28] + o[6]*(px.dfdx + 3.*px.dfdy + \
o[18] + o[29] + o[30] + o[31] + o[32])) + o[3]*o[25]*(o[27] + o[28] + o[6]*(2.*px.dfdy + o[18] + o[26] + o[29] + o[30] \
+ o[31] + o[32]));

/*-------------END OF OPTIMIZED CODE GENERATED BY EXTERNAL CODEGENERATOR----------------------------------------------*/

}


template <class val_t>
void Bezier_sq_ld_interpolation(const Tf_xy<val_t>& px,const Tf_xy<val_t>& py,Tf_xy<val_t> &res)
{
    val_t x=py.x;
    val_t x0=px.x;
    val_t y=px.y;
    val_t y0=py.y;
    res.x=x;
    res.y=y;
    val_t o[43];

/*-------------OPTIMIZED CODE GENERATED BY EXTERNAL CODEGENERATOR - !!!DO NOT EDIT!!! -------------------------------*/

o[0] = -1.*y0;
o[1] = y + o[0];
o[2] = pow(o[1],2.);
o[3] = o[1]*o[2];
o[4] = -1.*x0;
o[5] = -1.*y;
o[6] = x + y0 + o[4] + o[5];
o[7] = pow(o[6],2.);
o[8] = o[6]*o[7];
o[9] = 1/o[8];
o[10] = -2.*px.f;
o[11] = 2.*py.f;
o[12] = px.dfdx + py.dfdy;
o[13] = -1.*o[6]*o[12];
o[14] = o[10] + o[11] + o[13];
o[15] = 1/o[7];
o[16] = -3.*px.f;
o[17] = 3.*py.f;
o[18] = 2.*py.dfdy;
o[19] = px.dfdx + o[18];
o[20] = -1.*o[6]*o[19];
o[21] = o[16] + o[17] + o[20];
o[22] = -6.*px.f;
o[23] = 6.*py.f;
o[24] = pow(o[7],2.);
o[25] = 1/o[24];
o[26] = 2.*px.dfdx;
o[27] = 6.*px.f;
o[28] = -6.*py.f;
o[29] = px.dfdxy*x;
o[30] = -1.*px.dfdxy*x0;
o[31] = -1.*px.dfdxy*y;
o[32] = px.dfdxy*y0;
res.f = py.f + py.dfdy*o[1] - 1.*o[3]*o[9]*o[14] - 1.*o[2]*o[15]*o[21];
res.dfdx = py.dfdx + py.dfdxy*o[1] + o[2]*o[9]*(-1.*(px.dfdx + 3.*py.dfdx + 2.*(py.dfdy + py.dfdxy*(-1.*x + x0 + y + o[0\
])))*o[6] + o[22] + o[23]) + o[3]*o[25]*(o[22] + o[23] - 1.*o[6]*(2.*py.dfdx - 1.*py.dfdxy*x + py.dfdxy*x0 + py.dfdxy*y\
- 1.*py.dfdxy*y0 + o[18] + o[26]));
res.dfdy = py.dfdy - 3.*o[2]*o[9]*o[14] - 2.*o[1]*o[15]*o[21] + o[2]*o[9]*(o[27] + o[28] + o[6]*(px.dfdx + 3.*px.dfdy + \
o[18] + o[29] + o[30] + o[31] + o[32])) + o[3]*o[25]*(o[27] + o[28] + o[6]*(2.*px.dfdy + o[18] + o[26] + o[29] + o[30] \
+ o[31] + o[32]));

/*-------------END OF OPTIMIZED CODE GENERATED BY EXTERNAL CODEGENERATOR----------------------------------------------*/

}

}

#endif