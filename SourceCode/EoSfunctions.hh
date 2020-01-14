#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <string>
#include <sstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <fstream>
#include "CONSTANTSandCONVERSIONS.hh"
#include <vector>
using namespace std;
vector< vector<double> > NCHIDATA;
double SCANMIN;
double SCANMAX;
/*==============================================================================================================
==== Emperical Nuclear parameters - working in units of MeV and fm^-1 for energies and fm for distances ========
==============================================================================================================*/
//const double rho0 = 0.153;					// staturation density (fm^-3)
//double Kompress = 200*MeVtoinvFM;				// Compression modulus (fm^-1)
//double mass = 938*MeVtoinvFM;					// nucleon mass (fm^-1) 
//double massOmega = 783*MeVtoinvFM;				// omega mass (fm^-1) 
//double massSigma = 550*MeVtoinvFM;				// approximate sigma mass (fm^-1) 
//double mstar = 0.75*mass;					// effective mass
//const double asymm = 32.5*MeVtoinvFM;				// symmetry energy coefficient (fm^-1)
//const double BperA = -16.3*MeVtoinvFM;				// Binding energy (fm^-1)
double Pi = pi;
struct Rparams{
	double Msig; double Mw;	double Mr; double Mn; double Me; double gs; double gw; double gr;
	double bethe; double cethe; double nB; double Easymm; double Kcompress; double BperA;	double gamma; double asymmetryfactor;
	double Mlambda; double xHadronomega; double xHadronsigma; double xHadronrho; double Mmu; 
	double Mximinus; double Msigmaminus; double Msigmazero; double Msigmaplus; double Mxizero; 
	double DMDens; double GchibDM; double MpiDM; double MchiDM; double LambdaDM;
};
/*===============================================================================================================
================================= functions to be integrated ====================================================
===============================================================================================================*/
double I3(double k, void *param){			//relevant integral for solving set of eqns
	double Integrand = 0;
	double m = *(double *) param;
	
	double E = sqrt(k*k + m*m);			

	Integrand = pow(k,4) / pow(E,3);
	return Integrand;
}

double I2(double k, void *param){			//relevant integral for solving set of eqns
	double Integrand = 0;
	double m = *(double *) param;
	
	double k2 = k*k;				//squared cuz ima lazy
	double E = sqrt(k2 + m*m);			

	Integrand = k2* E;
	return Integrand;
}

double I1(double k, void *param){			//relevant integral for solving set of eqns
	double Integrand = 0;
	double m = *(double *) param;
	
	double k2 = k*k;				//squared cuz ima lazy
	double E = sqrt(k2 + m*m);			

	Integrand = k2* m / E;
	return Integrand;
}

double IP(double k, void *param){			//relevant integral for solving set of eqns
	double Integrand = 0;
	double m = *(double *) param;
	
	double k2 = k*k;				//squared cuz ima lazy
	double E = sqrt(k2 + m*m);			

	Integrand = k2* k2 / E;
	return Integrand;
}

/*=================================== DM CONSTRAINED PARAMETERS ===========================================*/
double Gchi(double sigchi2, double mchi, double mpi){
	double arg = sqrt(sigchi2) * mpi*mpi/( sqrt(20.)* mchi);
	return sqrt(arg);
}
double Gpi(double sigchin, double mn, double mchi, double mpi, double gchi){
	double arg = sqrt(sigchin)* (mn + mchi)*mpi*mpi/( sqrt(3.)*gchi*mn*mchi);
	return sqrt(arg);
}

double Gchib(double sigchin, double mn, double mchi, double mpi){
	double arg = pi*sigchin* pow((mn + mchi),2.)*pow(mpi,4.)/pow(mn*mchi,2.);
	return pow(arg, 0.25);
}
double Lambdachi(double sigchi2, double mchi){
	double arg = sqrt(sigchi2*mchi*mchi*64*M_PI);
	return arg;
}
double Mpib(double sigchin, double mn, double mchi){
	double arg = pow(mn*mchi,2.)/(pow(mn + mchi,2.)*pi*sigchin);
	return pow(arg, 0.25);
}


