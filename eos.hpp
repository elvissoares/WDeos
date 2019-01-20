#ifndef _EOS_HPP_
#define _EOS_HPP_

#include <iostream>
#include <vector>
#include <stdlib.h>

#include "ions.hpp"
#include "photons.hpp"
#include "fermions.hpp"

struct EoS {

    bool ions_gas, photons_gas, fermions_gas;

    EoS(bool ig = true, bool pg = true, bool eg = true): ions_gas(ig), photons_gas(pg), pairs_gas(eg) {}

    Fermions pairs;
    Ions ions;
    Photons photons;
    
    void EvaluateEoS(const double& rho, const double& T, const std::vector<double> &X);

    double Ion_MolarAbundance(const std::vector<double> &X);
    
    double ElectronsPerBaryons(const std::vector<double> &X);
    
    double AdiabaticIndex();
    
    double SoundSpeed();
    
    double HeatCapacity_Volume();

};

#endif