#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cassert>
#include <vector>
#include <sstream>
#include <iomanip>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_num_threads() 0
#define omp_get_thread_num() 0
#endif

#include "cpp_constants_CGS.hpp"
#include "fermions.hpp"
#include "ions.hpp"
#include "photons.hpp"

#define tab "\t"

#define test1 1 // Calculate the electrons+ions+photons eos 
#define test2 0 // Make a table to the thermodynamic quantities for an electron-positron pairs eos as a function of log rho and log T
#define test3 0 // Make a table to the Helmholtz energy and its derivatives for an electron-positron pairs eos as a function of log rho and log T

using namespace std;

int main (void) {
    
    std::cerr << "The max number of threads is " << omp_get_max_threads() << std::endl;
   
   	time_t begin, end;

    begin = time(0);

#if test1 // Calcula as qtds termodinâmicas para diferentes temperaturas
{  
  double T[8] = {1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11};

  using PhysConstants::N_A;
  using PhysConstants::a;
  using PhysConstants::k_B;

  string OutName = "../output/pressure.dat";
  ofstream OutFile(OutName.c_str());  

  string OutName2 = "../output/energydensity.dat";
  ofstream OutFile2(OutName2.c_str());

  string OutName3 = "../output/entropydensity.dat";
  ofstream OutFile3(OutName3.c_str());

  string OutName4 = "../output/dPdrho.dat";
  ofstream OutFile4(OutName4.c_str());

  string OutName5 = "../output/dUdrho.dat";
  ofstream OutFile5(OutName5.c_str());

  string OutName6 = "../output/dPdT.dat";
  ofstream OutFile6(OutName6.c_str());
    
  string OutName7 = "../output/dUdT.dat";
  ofstream OutFile7(OutName7.c_str());

  string OutName8 = "../output/relation1.dat";
  ofstream OutFile8(OutName8.c_str());

  string OutName9 = "../output/relation2.dat";
  ofstream OutFile9(OutName9.c_str());
    
  string OutName10 = "../output/relation3.dat";
  ofstream OutFile10(OutName10.c_str());

  Fermions pairs;    
  Ions ions;
  Photons photons;
    
  OutFile << "#log(rho/mu)   P" << std::endl;
  OutFile4 << "#log(rho/mu)   dPdrho" << std::endl;
  OutFile5 << "#log(rho/mu)   dUdrho" << std::endl;
  OutFile6 << "#log(rho/mu)   dPdT" << std::endl;
  OutFile7 << "#log(rho/mu)   dUdT" << std::endl;
  OutFile8 << "#log(rho/mu)   P - rho^2 dEdrho - T dPdT" << std::endl;
  OutFile9 << "#log(rho/mu)   dUdT - T dSdT" << std::endl;
  OutFile10 << "#log(rho/mu)   dSdrho + rho^-2 PdT" << std::endl;

  double A = 12., Z = 6.;
    
  double rho, nion, P, U, S, dUdrho, dUdT, dPdT, dSdrho, dSdT;
        
  for(double logrho = -4; logrho < 11; logrho+=0.1){

      rho = pow(10,logrho);
      
      OutFile << pow(10,logrho) << " ";
      OutFile2 << pow(10,logrho) << " ";
      OutFile3 << pow(10,logrho) << " ";
      OutFile4 << pow(10,logrho) << " ";
      OutFile5 << pow(10,logrho) << " ";
      OutFile6 << pow(10,logrho) << " ";
      OutFile7 << pow(10,logrho) << " ";

      OutFile8 << pow(10,logrho) << " ";
      OutFile9 << pow(10,logrho) << " ";
      OutFile10 << pow(10,logrho) << " ";

      for(unsigned int i = 0; i < 8; i++){
        pairs.Evaluate(rho,T[i],Z/A);
        ions.Evaluate(rho,T[i],1/A);
        photons.Evaluate(rho,T[i]);

        nion = N_A*rho/A;
          
        S = pairs.s + photons.s + ions.s;

        P = pairs.p + photons.p + ions.p;
        U = pairs.u + photons.u + ions.u;
        dUdrho = pairs.dudrho + photons.dudrho + ions.dudrho;
        dPdT = pairs.dpdT + photons.dpdT + ions.dpdT;
        dUdT = pairs.dudT + photons.dudT + ions.dudT; 
        dSdrho = pairs.dsdrho + photons.dsdrho + ions.dsdrho;
          
        dSdT = pairs.dsdT + photons.dsdT + ions.dsdT;

        OutFile << P << " ";
        OutFile2 << U << " ";
        OutFile3 << S << " ";
        OutFile4 << pairs.dpdrho + ions.dpdrho << " ";
        OutFile5 << dUdrho << " ";
        OutFile6 << dPdT << " ";
        OutFile7 << dUdT << " ";
        OutFile8 << fabs(P - Q(rho)*dUdrho - T[i]*dPdT)/P  << " ";
        OutFile9 << fabs(dUdT - T[i]*dSdT)/dUdT  << " ";
        OutFile10 << fabs(dSdrho + dPdT/Q(rho))/fabs(dPdT/Q(rho))  << " ";
      }

      OutFile << endl;
      OutFile2 << endl;
      OutFile3 << endl;
      OutFile4 << endl;
      OutFile5 << endl;
      OutFile6 << endl;
      OutFile7 << endl;
      OutFile8 << endl;
      OutFile9 << endl;
      OutFile10 << endl;
  }
    
  OutFile.close(); OutFile2.close(); OutFile3.close(); OutFile4.close();
  OutFile5.close(); OutFile6.close(); OutFile7.close();OutFile8.close();OutFile9.close();
} 
#endif 

#if test2 // Make a table to the thermodynamic quantities for an electron-positron pairs eos as a function of log rho and log T
{
  double logrho_lo = -6, logrho_hi = 11;
  double Dlogrho = 0.1;
    
  double logT_lo = 4, logT_hi = 11;
  double DlogT= 0.1;
    
  double logrho = logrho_lo;     
  double logT = logT_lo;
    
  string OutName = "../output/eos-pairs.dat";
  ofstream OutFile(OutName.c_str());

  OutFile << "# EoS grid of quantities for by log(rho/mu) and log(T) "<< endl;
  OutFile << "# Number of x1 points =" << (logrho_hi-logrho_lo)/Dlogrho +1 << endl;
  OutFile << "# Number of x2 points =" << (logT_hi-logT_lo)/DlogT +1 << endl;
  OutFile << "# log(rho/mu_e)  logT  logP  logu  logs" << std::endl;

  Fermions pairs;  
  
  for(double logT = logT_lo; logT<logT_hi; logT+=DlogT){

    for(double logrho = logrho_lo; logrho<logrho_hi; logrho+=Dlogrho){

      pairs.Evaluate(pow(10,logrho),pow(10,logT),1.);
      
      OutFile << logrho << tab << logT << tab 
        << log10(pairs.p) << tab 
        << log10(pairs.u) << tab 
        << log10(pairs.s) << std::endl;
    }
  }  
  
} 
#endif 

#if test3 // Make a table to the Helmholtz energy and its derivatives for an electron-positron pairs eos as a function of log rho and log T
{
    double logrho_lo = -6, logrho_hi = 11;
    double drho = 0.2;
    
    double logT_lo = 4, logT_hi = 11;
    double dT = 0.2;
    
    string OutName = "../output/free_energy-pairs.dat";
    ofstream OutFile(OutName.c_str());
    
    OutFile << "# EoS grid of quantities for by Rho/mu_e and T "<< endl;
    OutFile << "# Number of x1 points =" << (logrho_hi-logrho_lo)/drho +1 << endl;
    OutFile << "# Number of x2 points =" << (logT_hi-logT_lo)/dT +1 << endl;
    OutFile << "# Rho/mu_e  T  F  dFdrho  dFdT  d2FdRho2  d2FdT2   d2FdRhodT" << std::endl;
    
    Fermions pairs;
    
    double Rho, T, F, dFdRho, dFdT, d2FdRho2, d2FdT2, d2FdRhodT, d3FdRho2dT, d3FdT2dRho, d4FdT2dRho2 ;
    
    OutFile << std::setprecision(12);

    for(double logT = logT_lo; logT<=logT_hi; logT+=dT){
        T = pow(10,logT);
        
        for(double logrho = logrho_lo; logrho<=logrho_hi; logrho+=drho){
            
            Rho = pow(10,logrho);
            
            pairs.Evaluate(Rho,T,1.);
            pairs.Evaluate_FreeEnergyDerivatives();
            
            F = pairs.u - T*pairs.s;
            dFdRho = pairs.p/Q(Rho);
            dFdT = -pairs.s;
            d2FdRho2 = pairs.dpdrho/pow(Rho,2) - 2*pairs.p/pow(Rho,3);
            d2FdT2 = -pairs.dsdT;
            d2FdRhodT = pairs.dpdT/Q(Rho);
            d3FdRho2dT = pairs.d3fdrho2dT;
            d3FdT2dRho = pairs.d3fdT2drho;
            d4FdT2dRho2 = pairs.d4fdrho2dT2;

            
            OutFile << logrho << tab << logT << tab
            << F << tab
            << dFdRho << tab
            << dFdT << tab
            << d2FdRho2 << tab
            << d2FdT2 << tab
            << d2FdRhodT << tab
            << d3FdRho2dT  << tab
            << d3FdT2dRho << tab
            << d4FdT2dRho2  << std::endl;
        }
        
        OutFile << std::endl;
    }  
    
} 
#endif
    
	end = time(0);

	std::cerr << "Time elapsed: " << difftime(end,begin) << " s"<< std::endl;
    
	return 0;
}



