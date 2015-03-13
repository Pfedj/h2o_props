/*! Читалка данных из файла .h5 с таблицами свойств воды */

#ifndef _H2O_H5_DATA_H
#define _H2O_H5_DATA_H

#include <vector>
#include <string>
#include <h5i.h>

namespace H5_H2O

{

template<class Tval>
class Th2o_1d_table
{
    std::vector<Tval> data_buf;
public:
    const Tval* data;
    int nx;
    Th2o_1d_table(){}
    inline void ld(H5::CommonFG& file, const char* path)
    {
        std::string p(path);
        H5::Ld(file,(p+"/data").c_str(),data_buf,&nx);
        data=&data_buf[0];
    }
};

template<class Tval>
class Th2o_2d_table
{
    std::vector<Tval> data_buf;
    std::vector<Tval> datax_buf;
    std::vector<Tval> datay_buf;
    int dims[2];
public:
    const Tval* data;
    const Tval* datax;
    const Tval* datay;
    int nx,ny;
    Tval x1,x2,y1,y2;
    Th2o_2d_table(){}
    inline void ld(H5::CommonFG& file, const char* path)
    {
        std::string p(path);
        H5::Ld(file,(p+"/data").c_str(),data_buf,dims);
        H5::Ld(file,(p+"/datax").c_str(),datax_buf);
        H5::Ld(file,(p+"/datay").c_str(),datay_buf);
        data=&data_buf[0];
        datax=&datax_buf[0];
        datay=&datay_buf[0];
        H5::LdA(file,path,"nx",nx);
        H5::LdA(file,path,"ny",ny);
        H5::LdA(file,path,"x1",x1);
        H5::LdA(file,path,"x2",x2);
        H5::LdA(file,path,"y1",y1);
        H5::LdA(file,path,"y2",y2);
    }
};

template<class Tval>
class Th2o_common_1d_tables
{
public:
    Th2o_1d_table<Tval> Density;
    Th2o_1d_table<Tval> Enthalpy;
    inline void ld(H5::CommonFG& f, const char* path)
    {
        std::string p(path);
        Density.ld(f,(p+"/Density").c_str());
        Enthalpy.ld(f,(p+"/Enthalpy").c_str());
    }
};

template<class Tval>
class Th2o_gas_1d_tables:public Th2o_common_1d_tables<Tval>
{
};

template<class Tval>
class Th2o_liquid_1d_tables:public Th2o_common_1d_tables<Tval>
{
public:
    Th2o_1d_table<Tval> Pressure;
    Th2o_1d_table<Tval> Temperature;
    inline void ld(H5::CommonFG& f, const char* path)
    {
        std::string p(path);
        Th2o_common_1d_tables<Tval>::ld(f,path);
        Pressure.ld(f,(p+"/Pressure").c_str());
        Temperature.ld(f,(p+"/Temperature").c_str());
    }
};

template<class Tval>
class Th2o_common_2d_tables
{
public:
    Th2o_2d_table<Tval> Density;
    Th2o_2d_table<Tval> Temperature;
    inline void ld(H5::CommonFG& file, const char* path)
    {
        std::string p(path);
        Density.ld(file,(p+"/Density").c_str());
        Temperature.ld(file,(p+"/Temperature").c_str());
    }
};

template<class Tval>
class Th2o_fluid_2d_tables: public Th2o_common_2d_tables<Tval>
{
};

template<class Tval>
class Th2o_gas_2d_tables: public Th2o_common_2d_tables<Tval>
{
};

template<class Tval>
class Th2o_liquid_2d_tables: public Th2o_common_2d_tables<Tval>
{
};

template<class Tval>
class Th2o_fluid_data: public Th2o_fluid_2d_tables<Tval>
{
public:
    inline void ld(H5::CommonFG& f, const char* path)
    {
        std::string p(path);
        Th2o_fluid_2d_tables<Tval>::ld(f,p.c_str());
    }
};

template<class Tval>
class Th2o_gas_data: public Th2o_gas_2d_tables<Tval>
{
public:
    Th2o_gas_1d_tables<Tval> Saturated;
    inline void ld(H5::CommonFG& f, const char* path)
    {
        std::string p(path);
        Th2o_gas_2d_tables<Tval>::ld(f,p.c_str());
        Saturated.ld(f,(p+"/Saturated").c_str());
    }
};

template<class Tval>
class Th2o_liquid_data: public Th2o_liquid_2d_tables<Tval>
{
public:
    Th2o_liquid_1d_tables<Tval> Saturated;
    inline void ld(H5::CommonFG& f, const char* path)
    {
        std::string p(path);
        Th2o_liquid_2d_tables<Tval>::ld(f,p.c_str());
        Saturated.ld(f,(p+"/Saturated").c_str());
    }
};

template<class Tval>
class Th2o_data_parameters
{
public:
    Tval CriticalDensity;
    Tval CriticalEnthalpy;
    Tval CriticalEnthalpyW;
    Tval CriticalEnthalpyS;
    Tval CriticalTemperature;
    Tval CriticalPressure;
    Tval MinimalPressure;
};

template<class Tval>
class Th2o_data
{
public:
    Th2o_gas_data<Tval> gas;
    Th2o_fluid_data<Tval> fluid;
    Th2o_liquid_data<Tval> liquid;
    Th2o_data_parameters<Tval> attrs;
    Th2o_data(){}
    Th2o_data(const char* f){ld(f,"H2O");}
    inline void ld(const char* fn, const char* path)
    {
        H5::TH5File f(fn);
        std:: string p(path);
        gas.ld(f,(p+"/gas").c_str());
        fluid.ld(f,(p+"/fluid").c_str());
        liquid.ld(f,(p+"/liquid").c_str());
        H5::LdA(f,path,"CriticalDensity",attrs.CriticalDensity);
        H5::LdA(f,path,"CriticalEnthalpy",attrs.CriticalEnthalpy);
        H5::LdA(f,path,"CriticalEnthalpyW",attrs.CriticalEnthalpyW);
        H5::LdA(f,path,"CriticalEnthalpyS",attrs.CriticalEnthalpyS);
        H5::LdA(f,path,"CriticalTemperature",attrs.CriticalTemperature);
        H5::LdA(f,path,"CriticalPressure",attrs.CriticalPressure);
        H5::LdA(f,path,"MinimalPressure",attrs.MinimalPressure);
    }
};

}

#endif
