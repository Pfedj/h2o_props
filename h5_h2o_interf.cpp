#include "h5_h2o_interf.h"
#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/numeric.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/python/list.hpp>
#include <boost/python/dict.hpp>
#include <vector>


using namespace H5_H2O;
namespace bp=boost::python;
namespace bn=boost::python::numeric;

template class Th5_h2o_single<double>;
template class Th5_h2o_parallel<double>;

boost::shared_ptr<Th5_h2o_parallel<double> > ptr_h2o_par(nullptr);

void load_par(const char* f)
{
    ptr_h2o_par.reset(new Th5_h2o_parallel<double>(f));
}

template< typename T >
std::vector< T > to_std_vector( const bp::object& iterable )
{
    return std::vector< T >( bp::stl_input_iterator< T >( iterable ),
                             bp::stl_input_iterator< T >( ) );
}

template<class T>
inline bp::object to_python_array(const std::vector<T>& v)
{
    return bn::array(v.begin(),v.end());
}

template<class T>
inline bp::list to_python_array(const T* v, int n)
{
    bp::list res;
    for(int i=0;i<n;i++)
    {
        res.append(v[i]);
    }
    return res;
}

typedef Th2o_common_props_data_base<double> h2o_data_t;
typedef Th2o_common_props_data_ref<double> h2o_data_ref_t;
typedef Th2o_common_props_data_arr<double> h2o_data_arr_t;
typedef Tf_x<double> fx_t;
typedef Tsat_props<double> h2o_data_sat_t;

struct empty{};

h2o_data_t h2o_water_s(double P, double H)
{
    h2o_data_t res;
    ((Th5_h2o_single<double>*)(ptr_h2o_par.get()))->water_props_PH(P,H,create_h2o_data_ref(res));
    return res;
}

void h2o_water_s_multi(double *P, double *H, h2o_data_t* res, int N)
{
    auto xx=((Th5_h2o_single<double>*)(ptr_h2o_par.get()));
    for(int i=0;i<N;i++)
    {
      xx->water_props_PH(P[i],H[i],create_h2o_data_ref(res[i]));
    }
}

void h2o_sat_multi(double *P, double* sat, int N)
{
      Tsat_props<double> res[N];
      typedef Th2o_splines<double> h2os_t;
      for(int i=0;i<N;i++)
      {
       sat[i].T=ptr_h2o_par->Prop_sat(P[i],H2O_PROP_SAT_FUNCS(Tsat));
       sat[i].HW=ptr_h2o_par->Prop_sat(P[i],H2O_PROP_SAT_FUNCS(HWsat));
       sat[i].HS=ptr_h2o_par->Prop_sat(P[i],H2O_PROP_SAT_FUNCS(HSsat));
       sat[i].RoW=ptr_h2o_par->Prop_sat(P[i],H2O_PROP_SAT_FUNCS(RoWsat));
       sat[i].RoS=ptr_h2o_par->Prop_sat(P[i],H2O_PROP_SAT_FUNCS(RoSsat));    
      }

}

template<class T>
void spline_proph(T& P, T& H,	T &V, T &Temp, T &XMAS,
	T & EPS1, T& EPS2, T & VW,
	T & VS, T & HW, T& HS,
	T & VD, T & HD,T & CP,
	T & TS)
{
      const h2o_data_t& res_ = h2o_water_s(P,H);
      const h2o_data_sat_t& sat_ = h2o_sat_s(P);
      V = 1./res_.Ro;
      Temp = res_.T;
      VW = 1./sat_.RoW;      
      VS = 1./sat_.RoS;
      HW = sat_.HW;
      HS = sat_.HS;
      TS = sat_.T;
	XMAS = (H- HW)/(HW - HS);
      // EPS1 - dV/dP = -1 drho_dp/rho^2 {M**3/(KG*Pa)} 
      // EPS2 - dV/dH = -1 drho_dh/rho^2{M**3/J}
      T rev_rho = 1./(res_.Ro*res_.Ro);
      EPS1 = -1 * res_.dRo_dP * rev_rho;
      EPS2 = -1 * res_.dRo_dH * rev_rho;
      // VD = VS-VW  {M**3/KG} 
      // HD = HW-HS  {J/KG}    
      VD = VS - VW;
      HD = HW - HS;
      // Cp = dH_dT = (dT_dH)^-1
      Cp = -1./res_.dT_dH;     
      
}

void spline_proph_multi(T* P, T* H,	T *V, T *Temp, T *XMAS,
	T * EPS1, T* EPS2, T * VW,
	T * VS, T * HW, T* HS,
	T * VD, T * HD,T * CP,
	T * TS, int N)
{   
      h2o_data_t res[N];
      h2o_data_sat_t sat[N];
      
      h2o_water_s_multi(P, H, res, N);
      h2o_sat_multi(P, H, sat, N);
      
      for(int i=0;i<N;i++)
      {
            V[i] = 1./res_.Ro[i];
            Temp[i] = res_.T[i];
            VW[i] = 1./sat_.RoW[i];      
            VS[i] = 1./sat_.RoS[i];
            HW[i] = sat_.HW[i];
            HS[i] = sat_.HS[i];
            TS[i] = sat_.T[i];
            XMAS[i] = (H[i]- HW[i])/(HW[i] - HS[i]);
            // EPS1 - dV/dP = -1 drho_dp/rho^2 {M**3/(KG*Pa)} 
            // EPS2 - dV/dH = -1 drho_dh/rho^2{M**3/J}
            T rev_rho[i] = 1./(res_.Ro[i]*res_.Ro[i]);
            EPS1[i] = -1 * res_.dRo_dP[i] * rev_rho[i];
            EPS2[i] = -1 * res_.dRo_dH[i] * rev_rho[i];
            // VD = VS-VW  {M**3/KG} 
            // HD = HW-HS  {J/KG}    
            VD[i] = VS[i] - VW[i];
            HD[i] = HW[i] - HS[i];
            // Cp = dH_dT = (dT_dH)^-1
            Cp[i] = -1./res_.dT_dH[i];
      }      
}
      

static bp::object co;

template<class fun>
bp::object h2o_common_p(bp::object P, bp::object H, fun f)
{
    auto vP=to_std_vector<double>(P);
    auto vH=to_std_vector<double>(H);
    int n=vP.size();
    h2o_data_arr_t res(n);
    (ptr_h2o_par.get()->*f)(&vP[0],&vH[0],res,n);
    bp::object o=co();
    o.attr("P")=to_python_array(res.P,n);
    o.attr("H")=to_python_array(res.H,n);
    o.attr("Ro")=to_python_array(res.Ro,n);
    o.attr("dRo_dP")=to_python_array(res.dRo_dP,n);
    o.attr("dRo_dH")=to_python_array(res.dRo_dH,n);
    o.attr("T")=to_python_array(res.T,n);
    o.attr("dT_dP")=to_python_array(res.dT_dP,n);
    o.attr("dT_dH")=to_python_array(res.dT_dH,n);
    return o;
}

bp::object h2o_water_p(bp::object P, bp::object H)
{
    return h2o_common_p(P,H,&Th5_h2o_parallel<double>::water_props_PH);
}

bp::object h2o_steam_p(bp::object P, bp::object H)
{
    return h2o_common_p(P,H,&Th5_h2o_parallel<double>::steam_props_PH);
}

#define H2O_PROP_SAT_FUNCS(fun) &h2os_t::fun,&h2os_t::d##fun,&h2os_t::d2##fun

h2o_data_sat_t h2o_sat_s(double P)
{
    Tsat_props<double> res;
    typedef Th2o_splines<double> h2os_t;
    res.T=ptr_h2o_par->Prop_sat(P,H2O_PROP_SAT_FUNCS(Tsat));
    res.HW=ptr_h2o_par->Prop_sat(P,H2O_PROP_SAT_FUNCS(HWsat));
    res.HS=ptr_h2o_par->Prop_sat(P,H2O_PROP_SAT_FUNCS(HSsat));
    res.RoW=ptr_h2o_par->Prop_sat(P,H2O_PROP_SAT_FUNCS(RoWsat));
    res.RoS=ptr_h2o_par->Prop_sat(P,H2O_PROP_SAT_FUNCS(RoSsat));
    return res;
}

BOOST_PYTHON_MODULE(h5_h2o)
{
    bn::array::set_module_and_type( "numpy", "ndarray");
    bp::class_<h2o_data_t>("h2o_props_base",bp::init<>())
                        .def_readwrite("P",&h2o_data_t::P)
                        .def_readwrite("H",&h2o_data_t::H)
                        .def_readwrite("Ro",&h2o_data_t::Ro)
                        .def_readwrite("dRo_dP",&h2o_data_t::dRo_dP)
                        .def_readwrite("dRo_dH",&h2o_data_t::dRo_dH)
                        .def_readwrite("T",&h2o_data_t::T)
                        .def_readwrite("dT_dP",&h2o_data_t::dT_dP)
                        .def_readwrite("dT_dH",&h2o_data_t::dT_dH)
                        ;
                        
    bp::class_<h2o_data_sat_t>("h2o_data_sat_t",bp::init<>())
                        .def_readwrite("T",&h2o_data_sat_t::T)
                        .def_readwrite("HW",&h2o_data_sat_t::HW)
                        .def_readwrite("HS",&h2o_data_sat_t::HS)
                        .def_readwrite("RoW",&h2o_data_sat_t::RoW)
                        .def_readwrite("RoS",&h2o_data_sat_t::RoS)
                        ;

    bp::class_<fx_t>("fx_t",bp::init<>())
                        .def_readwrite("f",&fx_t::f)
                        .def_readwrite("df",&fx_t::df)
                        .def_readwrite("d2f",&fx_t::d2f)
                        .def_readwrite("x",&fx_t::x)
                        ;
                        
    co=bp::class_<empty>("",bp::init<>());
                        
    bp::def("h2o_water_s",h2o_water_s);
    bp::def("h2o_water_p",h2o_water_p);
    bp::def("h2o_steam_p",h2o_steam_p);
    bp::def("load",load_par);
    bp::def("h2o_sat_s",h2o_sat_s);
}
//