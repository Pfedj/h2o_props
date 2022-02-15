#ifndef _H5_H2O_PROPS_H
#define _H5_H2O_PROPS_H

#include "h5_h2o_alg.h"
#include "h5_h2o_spline.h"

namespace H5_H2O

{

template<class Tval>
struct Tsat_props
{
    typedef Tf_x<Tval> fx_t;
    fx_t T;
    fx_t RoS;
    fx_t RoW;
    fx_t HS;
    fx_t HW;
};


template<class Tval>
struct Th2o_common_props_data_base
{
    Tval P;
    Tval H;
    Tval T;
    Tval Ro;
    Tval dRo_dP;
    Tval dRo_dH;
    Tval dT_dP;
    Tval dT_dH;
    Tval* vals_p[8];
    Th2o_common_props_data_base()
    {
        vals_p[0]=&P;
        vals_p[1]=&H;
        vals_p[2]=&T;
        vals_p[3]=&Ro;
        vals_p[4]=&dRo_dP;
        vals_p[5]=&dRo_dH;
        vals_p[6]=&dT_dP;
        vals_p[7]=&dT_dH;
        size=8;
    }
    int size;
};

template<class Tval>
struct Th2o_common_props_data_ref: public Th2o_common_props_data_base<Tval*>
{
protected:
    bool cr;
public:
    Th2o_common_props_data_ref(bool x=false):Th2o_common_props_data_base<Tval*>()
    {
        if(x)
        {   
            P=new Tval[8];
            H=P+1;
            T=P+2;
            Ro=P+3;
            dRo_dP=P+4;
            dRo_dH=P+5;
            dT_dP=P+6;
            dT_dH=P+7;
        }
        cr=x;
    }
    ~Th2o_common_props_data_ref()
    {
        if(cr) 
        {
            delete P;
        }
    }
    inline Th2o_common_props_data_ref<Tval> operator[](int i)
    {
        Th2o_common_props_data_ref<Tval> r;
        for(int k=0;k<size;k++)
        {
            *(r.vals_p[k])=*(vals_p[k])+i;
        }
    }
    inline Th2o_common_props_data_ref<Tval> operator[](int i)const
    {
        Th2o_common_props_data_ref<Tval> r;
        for(int k=0;k<size;k++)
        {
           *(r.vals_p[k])=*(vals_p[k])+i;
        }
        return r;
    }
};

template <class Tval>
struct Th2o_common_props_data_arr: public Th2o_common_props_data_ref<Tval>
{
    int arr_size;
    Th2o_common_props_data_arr(int len):Th2o_common_props_data_ref<Tval>()
    {
        arr_size=len;
        for(int i=0;i<size;i++)
        {
            *vals_p[i]=new Tval[len];
        }
    }
    ~Th2o_common_props_data_arr()
    {
        for(int i=0;i<size;i++)
        {
            delete *vals_p[i];
        }
    }
};

template<class Tval>
Th2o_common_props_data_ref<Tval> create_h2o_data_ref(Th2o_common_props_data_base<Tval>& x)
{
     Th2o_common_props_data_ref<Tval> res;
     for(int i=0;i<x.size;i++)
     {
        *(res.vals_p[i])=x.vals_p[i];
        //*((Tval**)(&res)+i)=(Tval*)(&x)+i;
     }
     return res;
}

template<class val_t>
struct Tgeneral_props_list
{
    val_t P,H,rho,T,drho_dp,drho_dh,drho_dps,dT_dp,dT_dps,dT_dh,Cp,dh_dps;
    val_t dh_dps2,drho_dps2;
    val_t drho_dpp,drho_dhh,drho_dph;
    val_t dT_dpp,dT_dhh,dT_dph;
};

template<class Tval>
class Th2o_calculator
{
    typedef Th2o_splines<Tval> h2o_splines_t;
    const h2o_splines_t& S;
    Tval p700;
    Tval pcrit;
    Tval hw700;
    Tval hwcrit;
    Tval idh700;
    Tval idp700;
    Tval rhokr;
    Tval Tkrit;
    
    template<typename func>
    void get_FE(Tval P, Tval H, func Hsat, func dHsat, Tf_xy<Tval>& Phi, Tf_xy<Tval>& Eta, Tf_x<Tval>& Hs)const
    {
        if(P<p700)P=p700;
        Phi.x=P;
        Phi.y=H;
        Eta.x=P;
        Eta.y=H;
        if(P<pcrit)
        {
            Hs.x=P;
            Hs.f=(S.*Hsat)(P);
            Hs.df=(S.*dHsat)(P);
            Phi.f=(S.HWsat(P)-hw700)*idh700;
            Phi.dfdx=S.dHWsat(P)*idh700;
            Phi.dfdy=0.;
            Tval tmp=1./Hs.f;
            Eta.f=H*tmp;
            Eta.dfdx=-H*Hs.df*tmp*tmp;
            Eta.dfdy=tmp;
        }
        else
        {
            Hs.x=pcrit;
            Hs.f=(S.*Hsat)(pcrit);
            Hs.df=0.;
            Phi.f=idp700*(P-pcrit);
            Phi.dfdx=idp700;
            Phi.dfdy=0.;
            Eta.dfdx=0.;
            Eta.dfdy=1./Hs.f;
            Eta.f=H*Eta.dfdy;
        }
    }
    
    /*!Преобразование безразмерных термодинамических свойств в размерные, включая производные первых двух порядков*/
    template<int NV>
    inline void convert_to_original(const Tf_xy<Tval> F[NV],
        const Tval s[NV], const Tf_xy<Tval>& Phi, const Tf_xy<Tval>& Eta, Tf_xy<Tval> res[NV])const 
    {
    /*...*/
        Tval o[6];
o[0] = pow(Eta.dfdx,2);
o[1] = pow(Phi.dfdx,2);
o[2] = pow(Eta.dfdy,2);
o[3] = pow(Phi.dfdy,2);
        for(int i=0;i<NV;i++)
        {
res[i].f = F[i].f*s[i];
res[i].dfdx = (Eta.dfdx*F[i].dfdy + F[i].dfdx*Phi.dfdx)*s[i];
res[i].dfdy = (Eta.dfdy*F[i].dfdy + F[i].dfdx*Phi.dfdy)*s[i];
res[i].dfdx2 = s[i]*(Eta.dfdx2*F[i].dfdy + 2.*Eta.dfdx*F[i].dfdxy*Phi.dfdx + F[i].dfdx*Phi\
.dfdx2 + F[i].dfdy2*o[0] + F[i].dfdx2*o[1]);
res[i].dfdy2 = s[i]*(Eta.dfdy2*F[i].dfdy + 2.*Eta.dfdy*F[i].dfdxy*Phi.dfdy + F[i].dfdx*Phi\
.dfdy2 + F[i].dfdy2*o[2] + F[i].dfdx2*o[3]);
res[i].dfdxy = (Eta.dfdxy*F[i].dfdy + Eta.dfdy*(Eta.dfdx*F[i].dfdy2 + F[i].dfdxy*Phi.dfdx)\
 + F[i].dfdx*Phi.dfdxy + (Eta.dfdx*F[i].dfdxy + F[i].dfdx2*Phi.dfdx)*Phi.dfdy)*s[i];        }

//
    }
    inline void saturated_water_props_P(Tval P,Tgeneral_props_list<Tval>& wt,int& err)const 
    {
        err=0;
//        val_t Fi((P-p700)/(pkrit-p700));
        if(P>pcrit)
        {
            P=pcrit;
            //P=pkrit;
            err=1;
        }
        else if(P<p700)
        {
            P=p700;
            //P=p700;
            err=2;
        }
        wt.P=P;
        wt.rho=S.RoWsat(P);       
        wt.T=S.Tsat(P);           
        wt.drho_dps=S.dRoWsat(P); 
        wt.dT_dps=S.dTsat(P);     
        wt.H=S.HWsat(P);          
        wt.dh_dps=S.dHWsat(P);    
        wt.dh_dps2=S.d2HWsat(P);  
        Tf_xy<Tval> Fi,Eta;
        get_FE(P,wt.H,&h2o_splines_t::HWsat,&h2o_splines_t::dHWsat,Fi,Eta,Tf_x<Tval>());
        Tf_xy<Tval> w[2];
        Tf_xy<Tval> dervs[2];
        Tval scales[2];
        scales[0]=rhokr;
        scales[1]=Tkrit;
        S.liquid(Fi.f,Eta.f,w);
        S.liquid.d2(Fi.f,Eta.f,w);
        
        convert_to_original<2>(w,scales,Fi,Eta,dervs);
        
        wt.drho_dp=dervs[0].dfdx;
        wt.drho_dh=dervs[0].dfdy;
        wt.drho_dpp=dervs[0].dfdx2;
        wt.drho_dhh=dervs[0].dfdy2;
        wt.drho_dph=dervs[0].dfdxy;
        wt.dT_dp=dervs[1].dfdx;
        wt.dT_dh=dervs[1].dfdy;
        wt.dT_dpp=dervs[1].dfdx2;
        wt.dT_dhh=dervs[1].dfdy2;
        wt.dT_dph=dervs[1].dfdxy;
        wt.Cp=1./wt.dT_dh;
    }
    inline void saturated_steam_props_P(Tval P,Tgeneral_props_list<Tval>& wt,int& err)const 
    {
        err=0;
//        val_t Fi((P-p700)/(pkrit-p700));
        if(P>pcrit)
        {
            P=pcrit;
            //P=pkrit;
            err=1;
        }
        else if(P<p700)
        {
            P=p700;
            //P=p700;
            err=2;
        }
        wt.P=P;
        wt.rho=S.RoSsat(P);       
        wt.T=S.Tsat(P);           
        wt.drho_dps=S.dRoSsat(P); 
        wt.dT_dps=S.dTsat(P);     
        wt.H=S.HSsat(P);          
        wt.dh_dps=S.dHSsat(P);    
        wt.dh_dps2=S.d2HSsat(P);  
        Tf_xy<Tval> Fi,Eta;
        get_FE(P,wt.H,&h2o_splines_t::HSsat,&h2o_splines_t::dHSsat,Fi,Eta,Tf_x<Tval>());
        Tf_xy<Tval> w[2];
        Tf_xy<Tval> dervs[2];
        Tval scales[2];
        scales[0]=rhokr;
        scales[1]=Tkrit;
        S.gas(Fi.f,Eta.f,w);
        S.gas.d2(Fi.f,Eta.f,w);
        
        convert_to_original<2>(w,scales,Fi,Eta,dervs);
        
        wt.drho_dp=dervs[0].dfdx;
        wt.drho_dh=dervs[0].dfdy;
        wt.drho_dpp=dervs[0].dfdx2;
        wt.drho_dhh=dervs[0].dfdy2;
        wt.drho_dph=dervs[0].dfdxy;
        wt.dT_dp=dervs[1].dfdx;
        wt.dT_dh=dervs[1].dfdy;
        wt.dT_dpp=dervs[1].dfdx2;
        wt.dT_dhh=dervs[1].dfdy2;
        wt.dT_dph=dervs[1].dfdxy;
        wt.Cp=1./wt.dT_dh;
    }
    inline void critical_props_H(Tval H, Tgeneral_props_list<Tval>& wt,int& err)const
    {
        Tf_xy<Tval> Fi,Eta;
        auto f_hs=&h2o_splines_t::HSsat;
        auto d_hs=&h2o_splines_t::dHSsat;
        if(H<S.attrs.CriticalEnthalpyW)
        {
            f_hs=&h2o_splines_t::HWsat;
            d_hs=&h2o_splines_t::dHWsat;
        }
        get_FE(pcrit,H,&h2o_splines_t::HSsat,&h2o_splines_t::dHSsat,Fi,Eta,Tf_x<Tval>());
        Tf_xy<Tval> r[2],dervs[2];
        S.fluid(Fi.f,Eta.f,r);
        S.fluid.d2(Fi.f,Eta.f,r);
        Tval scales[2]={rhokr,Tkrit};
        convert_to_original(r,scales,Fi,Eta,dervs);
        wt.rho=dervs[0].f;
        wt.drho_dp=dervs[0].dfdx;
        wt.drho_dh=dervs[0].dfdy;
        wt.drho_dpp=dervs[0].dfdx2;
        wt.drho_dhh=dervs[0].dfdy2;
        wt.drho_dph=dervs[0].dfdxy;
        wt.T=dervs[1].f;
        wt.dT_dp=dervs[1].dfdx;
        wt.dT_dh=dervs[1].dfdy;
        wt.dT_dpp=dervs[1].dfdx2;
        wt.dT_dhh=dervs[1].dfdy2;
        wt.dT_dph=dervs[1].dfdxy;
        wt.Cp=1./wt.dT_dh;
    }
public:
    Th2o_calculator(const Th2o_splines<Tval>& s): S(s)
    {
        p700=s.attrs.MinimalPressure;
        hw700=s.hw700;
        pcrit=s.attrs.CriticalPressure;
        hwcrit=s.attrs.CriticalEnthalpyW;
        idh700=1./(hwcrit-hw700);
        idp700=1./(pcrit-p700);
        rhokr=s.attrs.CriticalDensity;
        Tkrit=s.attrs.CriticalTemperature;
    }
    inline void water_props_PH(Tval P, Tval H, Th2o_common_props_data_ref<Tval>& res)const
    {
        *res.P=P;
        *res.H=H;
        Tf_xy<Tval> Phi;
        Tf_xy<Tval> Eta;
        if(P>=p700)
        {
            if(P<pcrit)
            {
                Tgeneral_props_list<Tval> w;
                int ierr;
                saturated_water_props_P(P,w,ierr);
                get_FE(P,H,&h2o_splines_t::HWsat,&h2o_splines_t::dHWsat,Phi,Eta,Tf_x<Tval>());
                if(Eta.f<=1.)
                /*!Equilibrium liquid state*/
                {
                    Tf_xy<Tval> r[2];
                    S.liquid(Phi.f,Eta.f,r);
                    *res.Ro=rhokr*r[0].f;
                    *res.T=Tkrit*r[1].f;
                    *res.dRo_dP=rhokr*(r[0].dfdx*Phi.dfdx+r[0].dfdy*Eta.dfdx);
                    *res.dRo_dH=rhokr*(r[0].dfdx*Phi.dfdy+r[0].dfdy*Eta.dfdy);
                    *res.dT_dP =Tkrit*(r[1].dfdx*Phi.dfdx+r[1].dfdy*Eta.dfdx);
                    *res.dT_dH =Tkrit*(r[1].dfdx*Phi.dfdy+r[1].dfdy*Eta.dfdy);
                }
                else
                /*!UnEquilibrium (overheated) liquid state, the artifical extrapolation by the Bezier*/
                {
                    Tf_xy<Tval> wsat[NVARS],fluid[NVARS];
                    auto f_hs=&h2o_splines_t::HSsat;
                    auto d_hs=&h2o_splines_t::dHSsat;
                    if(H<S.attrs.CriticalEnthalpyW)
                    {
                        f_hs=&h2o_splines_t::HWsat;
                        d_hs=&h2o_splines_t::dHWsat;
                    }
                    get_FE(P,H,f_hs,d_hs,Phi,Eta,Tf_x<Tval>());
                    S.liquid(Phi.f,1.,wsat);
                    S.liquid.d2(Phi.f,1.,wsat);
                    S.fluid(1.,Eta.f,fluid);
                    S.fluid.d2(1.,Eta.f,fluid);
                    Tf_xy<Tval> px[NVARS],py[NVARS],r[NVARS];
                    Tval sc[]={rhokr,Tkrit};
                    for(int i=0;i<2;++i)
                    {
                        Bezier_sq_lu_interpolation(fluid[i],wsat[i],r[i]);
                        r[i].convert_to_ph(Phi,Eta,P,H);
                        r[i]*=sc[i];
                    }
                    *res.Ro=r[0].f;
                    *res.T=r[1].f;
                    *res.dRo_dP=r[0].dfdx;
                    *res.dRo_dH=r[0].dfdy;
                    *res.dT_dP =r[1].dfdx;
                    *res.dT_dH =r[1].dfdy;
                }
            }
            else
            /*!Supercrtitical state*/
            {
                get_FE(P,H,&h2o_splines_t::HWsat,&h2o_splines_t::dHWsat,Phi,Eta,Tf_x<Tval>());
                Tf_xy<Tval> r[2];
                S.fluid(Phi.f,Eta.f,r);
                *res.Ro=rhokr*r[0].f;
                *res.T=Tkrit*r[1].f;
                *res.dRo_dP=rhokr*(r[0].dfdx*Phi.dfdx+r[0].dfdy*Eta.dfdx);
                *res.dRo_dH=rhokr*(r[0].dfdx*Phi.dfdy+r[0].dfdy*Eta.dfdy);
                *res.dT_dP =Tkrit*(r[1].dfdx*Phi.dfdx+r[1].dfdy*Eta.dfdx);
                *res.dT_dH =Tkrit*(r[1].dfdx*Phi.dfdy+r[1].dfdy*Eta.dfdy);
            }
        }
    }
    inline void steam_props_PH(Tval P, Tval H, Th2o_common_props_data_ref<Tval>& res)const
    {
        *res.P=P;
        *res.H=H;
        Tf_xy<Tval> Phi;
        Tf_xy<Tval> Eta;
        if(P>=p700)
        {
            if(P<pcrit)
            {
                Tgeneral_props_list<Tval> w;
                int ierr;
                saturated_steam_props_P(P,w,ierr);
                get_FE(P,H,&h2o_splines_t::HSsat,&h2o_splines_t::dHSsat,Phi,Eta,Tf_x<Tval>());
                if(Eta.f>=1.)
                /*!Equilibrium vapor state*/
                {
                    Tf_xy<Tval> r[2];
                    S.gas(Phi.f,Eta.f,r);
                    *res.Ro=rhokr*r[0].f;
                    *res.T=Tkrit*r[1].f;
                    *res.dRo_dP=rhokr*(r[0].dfdx*Phi.dfdx+r[0].dfdy*Eta.dfdx);
                    *res.dRo_dH=rhokr*(r[0].dfdx*Phi.dfdy+r[0].dfdy*Eta.dfdy);
                    *res.dT_dP =Tkrit*(r[1].dfdx*Phi.dfdx+r[1].dfdy*Eta.dfdx);
                    *res.dT_dH =Tkrit*(r[1].dfdx*Phi.dfdy+r[1].dfdy*Eta.dfdy);
                }
                else
                /*!UnEquilibrium (underheated) vapor state, the artifical extrapolation by the Bezier*/
                {
                    Tf_xy<Tval> wsat[NVARS],fluid[NVARS];
                    auto f_hs=&h2o_splines_t::HSsat;
                    auto d_hs=&h2o_splines_t::dHSsat;
                    if(H<S.attrs.CriticalEnthalpyW)
                    {
                        f_hs=&h2o_splines_t::HWsat;
                        d_hs=&h2o_splines_t::dHWsat;
                    }
                    get_FE(P,H,f_hs,d_hs,Phi,Eta,Tf_x<Tval>());
                    S.gas(Phi.f,1.,wsat);
                    S.gas.d2(Phi.f,1.,wsat);
                    S.fluid(1.,Eta.f,fluid);
                    S.fluid.d2(1.,Eta.f,fluid);
                    Tf_xy<Tval> px[NVARS],py[NVARS],r[NVARS];
                    Tval sc[]={rhokr,Tkrit};
                    for(int i=0;i<2;++i)
                    {
                        Bezier_sq_ld_interpolation(fluid[i],wsat[i],r[i]);
                        r[i].convert_to_ph(Phi,Eta,P,H);
                        r[i]*=sc[i];
                    }
                    *res.Ro=r[0].f;
                    *res.T=r[1].f;
                    *res.dRo_dP=r[0].dfdx;
                    *res.dRo_dH=r[0].dfdy;
                    *res.dT_dP =r[1].dfdx;
                    *res.dT_dH =r[1].dfdy;
                }
            }
            /*!Supercrtitical state*/
            else
            {
                get_FE(P,H,&h2o_splines_t::HWsat,&h2o_splines_t::dHWsat,Phi,Eta,Tf_x<Tval>());
                Tf_xy<Tval> r[2];
                S.fluid(Phi.f,Eta.f,r);
                *res.Ro=rhokr*r[0].f;
                *res.T=Tkrit*r[1].f;
                *res.dRo_dP=rhokr*(r[0].dfdx*Phi.dfdx+r[0].dfdy*Eta.dfdx);
                *res.dRo_dH=rhokr*(r[0].dfdx*Phi.dfdy+r[0].dfdy*Eta.dfdy);
                *res.dT_dP =Tkrit*(r[1].dfdx*Phi.dfdx+r[1].dfdy*Eta.dfdx);
                *res.dT_dH =Tkrit*(r[1].dfdx*Phi.dfdy+r[1].dfdy*Eta.dfdy);
            }
        }
    }
};

}

#endif

