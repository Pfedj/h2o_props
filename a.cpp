#include "h5_h2o_data.h"
#include "h5_h2o_spline.h"
#include "h5_h2o_props.h"

#include <boost/python.hpp>

using namespace H5_H2O;

Th2o_calculator<double>* h2o_c;

void water_props_PH(double P, double H, Th2o_common_props_data_base<double>& res)
{
    h2o_c->water_props_PH(P,H,create_h2o_data_ref(res));
}

namespace boost{
 namespace python {     
  double extract(double* x){return *x;}
}
}

/*BOOST_PYTHON_MODULE(h5_h2o)
{
    using namespace boost::python;
    static Th2o_data<double> h2o_data;
    h2o_data.ld("a.h5","H2O");
    static Th2o_splines<double> h2o_spls(h2o_data);
    static Th2o_calculator<double> h2o(h2o_spls);
    h2o_c=&h2o;
    typedef Th2o_common_props_data_base<double> Tprp;
    class_<Tprp>("h2o_prop_set",init<>())
        .def_readwrite("P",&Tprp::P)
        .def_readwrite("H",&Tprp::H)
        .def_readwrite("Ro",&Tprp::Ro)
        .def_readwrite("dRo_dP",&Tprp::dRo_dP)
        .def_readwrite("dRo_dH",&Tprp::dRo_dH)
        .def_readwrite("T",&Tprp::T)
        .def_readwrite("dT_dP",&Tprp::dT_dP)
        .def_readwrite("dT_dH",&Tprp::dT_dH)
        ;
    def("water_props_PH",&water_props_PH);
}*/

int main()
{
    using namespace H5_H2O;
    static Th2o_data<double> h2o_data;
    h2o_data.ld("F:\\PROJECTS\\persey\\enisey_proj\\code\\H2Oprops_v\\a.h5","H2O");
    static Th2o_splines<double> h2o_spls(h2o_data);
    static Th2o_calculator<double> h2o(h2o_spls);
    h2o_c=&h2o;
    double hw=h2o_spls.HWsat(1e5);
    double hs=h2o_spls.HSsat(1e5);
    Th2o_common_props_data_ref<double> p(true);
    h2o.water_props_PH(1e5,1e5,p);
    return 0;
}