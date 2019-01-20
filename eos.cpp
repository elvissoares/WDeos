#include "eos.hpp"
#include "cpp_constants_CGS.hpp"

#include <cmath>
#include <stdio.h>
#include <string>
#include <fstream>

#include "util/roots.hpp"

void EvaluateEoS(const double& rho_input, const double& T_input, const std::vector<double> &X_input);
{
    
    rho = rho_input;
    T = T_input;
    
    Yi = Ion_MolarAbundance(X_input);
    Ye = ElectronsPerBaryons(X_input);
    
    if (ions_gas){
        ions.Evaluate(rho,T,Yi);
        
        P += ions.p;
        dPdrho += ions.dpdrho;
        dPdT += ions.dpdT;
        
        u += ions.u;
        dudrho += ions.dudrho;
        dudT += ions.dudT;
        
        s += ions.s;
        dsdrho += ions.dsdrho;
        dsdT += ions.dsdT;
    }
    
    if (photons_gas){
        photons.Evaluate(rho,T);
        
        P += photons.p;
        dPdrho += photons.dpdrho;
        dPdT += photons.dpdT;
        
        u += photons.u;
        dudrho += photons.dudrho;
        dudT += photons.dudT;
        
        s += photons.s;
        dsdrho += photons.dsdrho;
        dsdT += photons.dsdT;
    }
    
    if (pairs_gas){
        pairs.Evaluate(rho,T,Ye);
        
        P += pairs.p;
        dPdrho += pairs.dpdrho;
        dPdT += pairs.dpdT;
        
        u += pairs.u;
        dudrho += pairs.dudrho;
        dudT += pairs.dudT;
        
        s += pairs.s;
        dsdrho += pairs.dsdrho;
        dsdT += pairs.dsdT;
    }
    
    gamma = AdiabaticIndex();
    cs = SoundSpeed();
    cv = HeatCapacity_Volume;
    
}

//Calculate the mean molecular weight of the gas
double EoS::Ion_MolarAbundance(const std::vector<double> &X)
{
    //For complete ionization
    const int Z12 = 6;
    const int A12 = 12;

    const int Z16 = 8;
    const int A16 = 16;

    return (X[0]/A12 + X[1]/A16 + X[2]/24.);  //Carbon + Oxygen + Ashes
}

double EoS::ElectronsPerBaryons(const std::vector<double> &X)
{
    //For complete ionization
    const int Z12 = 6;
    const int A12 = 12;

    const int Z16 = 8;
    const int A16 = 16;

    return (X[0]*Z12/A12 + X[1]*Z16/A16 + X[2]*0.5); //Carbon + Oxygen + Ashes
}


double EoS::AdiabaticIndex()
{
    return (rho/P)*dPdrho;
}

double EoS::SoundSpeed()
{
    return (dPdrho - dPdT*dsdrho/dsdT);
}

double EoS::HeatCapacity_Volume()
{
    return (dudT);
}