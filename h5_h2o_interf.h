#ifndef _H5_H2O_INTER_H
#define _H5_H2O_INTER_H

#ifdef _DEBUG
#undef _DEBUG
#define _DEBUG 1
#else
#define _DEBUG 0
#endif

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include "h5_h2o_props.h"


namespace H5_H2O

{

template<class Tval>
class Th5_h2o_single
{
protected:
    const Th2o_splines<Tval> *spl;
    const Th2o_calculator<Tval> *calc;
public:
    Th5_h2o_single():spl(nullptr),calc(nullptr){}
    Th5_h2o_single(const char* f):spl(nullptr),calc(nullptr)
    {
        ld(f);
    }
    inline void ld(const char* f)
    {
        if(spl==nullptr)
        {
            spl=new Th2o_splines<Tval>(Th2o_data<Tval>(f));
            calc=new Th2o_calculator<Tval>(*spl);
        }
    }
    ~Th5_h2o_single()
    {
        if(spl!=nullptr) delete spl;
        if(calc!=nullptr) delete calc;
    }
    virtual void water_props_PH(Tval P, Tval H, Th2o_common_props_data_ref<Tval>& r)const
    {
        calc->water_props_PH(P,H,r);
    }
    virtual void steam_props_PH(Tval P, Tval H, Th2o_common_props_data_ref<Tval>& r)const
    {
        calc->steam_props_PH(P,H,r);
    }
    template<class fun>
    inline  Tf_x<Tval> Prop_sat(Tval P, fun f, fun df, fun d2f)
    {
        Tf_x<Tval> res;
        res.x=P;
        res.f=(spl->*f)(P);
        res.df=(spl->*df)(P);
        res.d2f=(spl->*d2f)(P);
        return res;
    }
};

template<class Tval>
class Th5_h2o_parallel: public Th5_h2o_single<Tval>
{
protected:
    template<class Tel_fun>
    struct run_element
    {
        const Th5_h2o_single<Tval>* el;
        Tel_fun fun;
        const Tval* P_in;
        const Tval* H_in;
        Th2o_common_props_data_ref<Tval> out;
        void operator()(const tbb::blocked_range<int>& range)const
        {
            for(int i=range.begin();i!=range.end();++i)
            {
                (el->*fun)(P_in[i],H_in[i],out[i]);
            }
        }
    };
    template<class Tel_fun>
    inline run_element<Tel_fun> generate_run_element(const Th5_h2o_single<Tval>& _el, Tel_fun f)const
    {
        run_element<Tel_fun> res;
        res.el=&_el;
        res.fun=f;
        return res;
    }
public:
    Th5_h2o_parallel(): Th5_h2o_single<Tval>(){}
    Th5_h2o_parallel(const char* f): Th5_h2o_single<Tval>(f){}
    virtual void water_props_PH(const Tval *P, const Tval *H, Th2o_common_props_data_ref<Tval>& r, int n)const
    {
        auto runner=generate_run_element(*(Th5_h2o_single<Tval>*)this,&Th5_h2o_single<Tval>::water_props_PH);
        runner.P_in=P;
        runner.H_in=H;
        runner.out=r;
        tbb::parallel_for(tbb::blocked_range<int>(0,n),runner);
    }
    virtual void steam_props_PH(const Tval *P, const Tval *H, Th2o_common_props_data_ref<Tval>& r, int n)const
    {
        auto runner=generate_run_element(*(Th5_h2o_single<Tval>*)this,&Th5_h2o_single<Tval>::steam_props_PH);
        runner.P_in=P;
        runner.H_in=H;
        runner.out=r;
        tbb::parallel_for(tbb::blocked_range<int>(0,n),runner);
    }
};


}

#endif