//#include <gsl/gsl_integration.h>
//#include "LHAPDF/LHAPDF.h"
#include <string>
#include <sstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <vector>
#include "CONSTANTSandCONVERSIONS.hh"

/*====================== EMPERICAL NS PARAMETERS =======================*/
const double Msolar = 1.98892e30*convkgtoGeV;						// mass of Sol in GeV
const double TNS = 1e5*convKtoGeV;								// NS temperature in K->GeV
const double TDM = TNS;
const double vbar = 220e5*convcmtoinvGeV/convstoinvGeV;					// mean thermal DM speed in cm/s->fraction of c - MB distrib
const double v2 = vbar*vbar;									// the square pops up more frequently
const double Nmass = 0.938;								// Neutron mass in GeV
double rhoDM = 0.47259/pow(convcmtoinvGeV,3.);						// GeV/cm^3 -> GeV^4 - mean DM density at 1 kpc from Galactic center

//double pFermi = Planckbar*pow((3*pi*pi*CoreDens), 1./3.)/clight;//GeV
const double pFermi = 0.426;									// Fermi momentum, neutrons, in GeV
const double T1 = 3.5;									// containment time	
const double TMAX = 1e10;


double tcon(double mDM, double sigchin, double Mns, double Rns){
	double sigcrit = Nmass*Rns*Rns/Mns;
	return 3.*pi*mDM*pow(Rns, 3./2.)*sigcrit/(4.*Nmass*sqrt(2.*GNewton0*Mns)*sigchin*pow(convcmtoinvGeV,2.))*sqrt(mDM/Nmass)/convstoinvGeV/(60*60*24*365);
}
double radius(double mDM, double mred, double t, double sigchin, double CDens, double Rns){
	double eta = 8*pi*sqrt(2.) * pow(mred, 3.) * pow(CDens, 2.) * sigchin*pow(convcmtoinvGeV,2.) *Rns*Rns* GNewton0 / (3.*mDM * pFermi);
	double r = Rns/sqrt(1.+eta*(t-T1)*365*24*60*60*convstoinvGeV);
	return r;
}	

/*============= ORIGINAL EVOLUTION =============*/
double Nchi0(double Cc, double Cs, double time){
	return Cc/Cs * (exp(Cs * time) - 1.);
}

/* ============ GEOMETRIC ======================*/
double tgeo(double mDM, double mred, double sigchin, double sigchi2, double Cc, double Cs, double CDens, double Rns){ 
	double time = T1;
	double tGEO = TMAX*2;
	while(time <= 1e10){
		double rad = radius(mDM, mred, time, sigchin, CDens, Rns)/Rns;
		double rfid = sqrt( Nchi0(Cc, Cs, time)* sigchi2*pow(convcmtoinvGeV,2.)/(pi * Rns *Rns));
		if (rad <= 1.001*rfid && rad >= 0.999 *rfid ){
			tGEO = time;	
		}
		time *= 1.0001;
	}
	if(tGEO < TMAX){
		return tGEO;
	}
	else if (tGEO > TMAX){
		return TMAX*2;
	}
}
double Nchigeo(double mDM, double mred, double geotime, double time, double sigchin, double sigchi2, double Cc, double Cs, double CDens, double Rns){
	double eta = 8*pi*sqrt(2.) * pow(mred, 3.) * pow(CDens, 2.) * sigchin*pow(convcmtoinvGeV,2.) *Rns*Rns* GNewton0 / (3*mDM * pFermi);		//NATURAL
	double radgeo = radius(mDM, mred, geotime, sigchin, CDens, Rns);												//NATURAL
	double Ngeo = pi* pow(radgeo, 2)/(sigchi2*pow(convcmtoinvGeV,2.));
	return Cc*time + Cs/(60*60*24*365)*convinvstoGeV*Ngeo*Rns*Rns/(pow(radgeo, 2)*eta) * log(1 + eta*(time - T1)*365*24*60*60*convstoinvGeV);
}


/* ============ THERMALIZATION ===================*/
double ttherm(double mDM, double mn, double mred, double sigchin, double CDens){
	return pFermi/(6.*sqrt(2.)*TNS*CDens*sigchin*pow(convcmtoinvGeV,2.)) * pow( (mn/mred), 3)*pow(mDM/mn,2.)/convstoinvGeV/(60*60*24*365);
}
double Nchitherm(double mDM, double mred, double thermtime, double time, double sigchin, double sigchi2, double Cc, double Cs, double CDens){
	double radtherm = sqrt(9.* TNS/(4.*pi*GNewton0*CDens*Nmass*mDM));
	return (Cc + Cs*pi*pow(radtherm, 2)/(sigchi2*pow(convcmtoinvGeV,2.) ) )* time;
}
double PauliBlocking(double mred, double Mns, double Rns){
	double xi = 2.* sqrt(GNewton0*Mns/Rns)*mred;
	if (xi < 1.){
		return xi;}
	else{ return 1.;}
}


