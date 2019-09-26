#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>
#include <string>
#include <fstream>
#include "CONSTANTSandCONVERSIONS.hh"
#include <time.h>
vector< vector<double> > EOSDATA;
vector<double> PRESSUREDATA;
vector<double> ENERGYDATA;

using namespace std;
double Msolar = 1.98e30;
double Msolar0 = Msolar*convkgtoGeV*1e-42;
double GNewtonPrime = GNewton0*1e10/convGeVtoinvcm*1e30;
double EnergyAux = 0.0;
double dummyplug;

double interpolateEOS(double p){
	int rows = EOSDATA.size();
	int n = 0;
	n = rows -1;
	double PRESS[rows];
	double ENDENS[rows];

	for (int i = 1; i < rows; i++){
		PRESS[i] = EOSDATA[i][2];
		ENDENS[i] = EOSDATA[i][1];
	}

	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	gsl_spline *spline = gsl_spline_alloc(gsl_interp_akima, n);
	gsl_spline_init (spline, PRESS, ENDENS, n);
	
	double eop = gsl_spline_eval(spline, p, acc);
	EnergyAux = eop;
	gsl_spline_free (spline);
	gsl_interp_accel_free (acc);
	return eop;
}

typedef struct{ double m; double r;} tov_results;
int tov_rhs(const double r, const double y[], double f[], void * params){
	if( r ==0){
		f[0] = 0.;
		f[1] = 0.;
	}
	else{
	double p = y[0];
        double m = y[1];
	double eps = interpolateEOS(p);
	f[0] = -(GNewtonPrime/(r*r))*(eps + p)*(m*Msolar0 + 4*pi*gsl_pow_3(r)*p)/(1.-2*m*Msolar0*GNewtonPrime/r);
	f[1] = 4*pi*gsl_pow_2(r)*eps/Msolar0;
	}
	return GSL_SUCCESS;
}
int tovinv_rhs(const double p, const double y[], double f[], void * params){
	double r = y[0];
	double m = y[1];
	double eps = interpolateEOS(p);
	f[0] = r * (2.0 * GNewtonPrime*m*Msolar0 - r) /((eps + p)*GNewtonPrime * (m*Msolar0 + 4.0 * pi * p * gsl_pow_3(r)));
	f[1] = f[0] * 4.0 * pi * eps * gsl_pow_2(r)/Msolar0;
	return GSL_SUCCESS;
}

int tov_solver_solve(tov_results * results, double CentralDensity, string EOSMODEL){

	double r;                   /* Radius (independent variable) */
	double p;
	double m;
	int status;                 /* Status returned by GSL functions */

/*==============================================================================
============================ READ DENSITY ======================================
* =============================================================================*/
	int rows = EOSDATA.size()-1;
	double BDENS[rows-1];
	double PRESS[rows-1];

	for (int i = 1; i < rows; i++){
		PRESS[i] = EOSDATA[i][2];
		BDENS[i] = EOSDATA[i][0];
	}

	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, rows);
	gsl_spline_init (spline, BDENS, PRESS, rows);
	double PRESSURE = gsl_spline_eval(spline, CentralDensity, acc);

	gsl_spline_free (spline);
	gsl_interp_accel_free (acc);

/*==============================================================================
============================ SET UP INTEGRATOR ONE ==============================
* =============================================================================*/

	double y[2];                /* ODE state */
	double eps_abs = 1.0e-18;   /* Absolute precision goal */
	double eps_rel = 1.0e-18;   /* Relative precision goal */
	const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk8pd;
	double h = 1e-6;

	gsl_odeiv2_step * s = gsl_odeiv2_step_alloc (T, 2);
	gsl_odeiv2_control * c = gsl_odeiv2_control_y_new (eps_abs, eps_rel);
	gsl_odeiv2_evolve * e = gsl_odeiv2_evolve_alloc (2);

	gsl_odeiv2_system sys;
	sys.params = NULL;
	sys.dimension = 2;
	sys.jacobian = NULL;
 	sys.function = &tov_rhs;

	r = 0.0;
	double rMax = 1e30;
	double pmin = PRESS[1];
	y[0] = PRESSURE;
	y[1] = 0.0;

/*===============================================================
 * =========== SET UP LIMITS BY MODEL TABLES =====================
 * =============================================================*/
	double lowerlim = 0;
	if(EOSMODEL == "DDHd_Y"){
		lowerlim = 1e-3;
	}
	else if(EOSMODEL == "GM1_Y"){
		lowerlim = 1e-8;
	}
	else if(EOSMODEL == "SigOm"){
		lowerlim = 1e-12;
	}
	else{
		cout << "Model not recognized - limits could not be set" << '\n';
	}

	while (y[0] > lowerlim){

		status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &r, rMax, &h, y);

		if (status != GSL_SUCCESS){
			cout << "ERROR: return value = " << status << '\n';
			break;
		}
	}
	gsl_odeiv2_evolve_free (e);
	gsl_odeiv2_control_free (c);
	gsl_odeiv2_step_free (s);
	results->r = r;
	results->m = y[1];

	return GSL_SUCCESS;
}



int main(int argc, char* argv[]){
/*==============================================================================
============================ NAMING SCHEME AND VARIABLES =======================
* =============================================================================*/
	int i;
	double BaryonDensity;
	double t1, t2;
	t1 = clock();
	tov_results results;
	int fno;// = argv[1];
	int rows = 0;	
	int FILES = 0;
	double pr1 = 0;
	double pr2 = 0;
	double pr3 = 0;
	double pr4 = 0;
	double pr5 = 0;
	double pr6 = 0;
	double pr7 = 0;
	double pr1p = 0;
	double pr2p = 0;
	double pr3p = 0;

	double title1 = atof(argv[1]);
	double title2 = atof(argv[2]);
	double title3 = atof(argv[3]);
	double title4 = atof(argv[4]);
	string STAR = argv[5];
	string EoSMODEL = argv[6];
	stringstream ss;
	ss << title1 << title2 << title3 << title4;
	string fnos = ss.str();

/*==============================================================================
============================ READ IN AND FIT EoS FILES USING INTERPOLATOR ======
* =============================================================================*/
	string particle = "BOSONS";
	stringstream liltitle1;
	liltitle1 <<"./" << particle << "/EoSFiles/"<<title4<<"GeV/EoS_SigOmDM_"<<fnos<<"p.dat";

	string	BIGTITLE1 = liltitle1.str();
	ifstream dater;
	dater.open(BIGTITLE1);

	if (!dater){
		cout << "Corresponding EoS File Not Found: " << " " << BIGTITLE1 <<  '\n';
		return 1;
	}

	int columns = 0;

	while(!dater.eof()){
		dater >> pr1 >> pr2 >> pr3 >> pr4 >> pr5 >> pr6;
		if(pr2> 1.1*pr2p and pr3>1.1*pr3p){
			vector<double> row;
			row.push_back(pr1);
			row.push_back(pr2);
			row.push_back(pr3);
			EOSDATA.push_back(row);
			columns = row.size();
			pr2p = pr2;
			pr3p = pr3;
		//cout << pr1 << " " << pr2 << " " << pr3 << '\n';
		}
		else{continue;}
	}

	dater.close();

	rows = EOSDATA.size();
	double BaryonDensityMax = EOSDATA[rows-3][0];
	cout << "rows: " << rows << " columns: " << columns << '\n';

//	cout << '\n' << "DATA COPYING DONE..." << '\n';

/*==============================================================================
============================ PREPARE OUTPUT ====================================
* =============================================================================*/
	stringstream liltitle2;
	liltitle2 << "./" << particle << "/MRFiles/" << title4 << "GeV/MR_" << fnos << "p.dat";
	string BIGTITLE2 = liltitle2.str();
	ofstream myfile2;
	myfile2.open(BIGTITLE2);
	int jk = 0;

/*==============================================================================
======== INTEGRATE THE EoS FOR ALL CENTRAL DENSITIES TO GET A SEQUENCE =========
* =============================================================================*/
	BaryonDensity = 0.01;

	while (BaryonDensity <= BaryonDensityMax){
		tov_solver_solve(&results, BaryonDensity, EoSMODEL);

		if((results.r*1e-3 >=0.) and (results.r*1e-3<=20)){
			myfile2 << results.m << " " << results.r*1e-3 << " " << BaryonDensity << " " << pr4 << " " << pr5 << " " << pr6 << '\n';
		}
		
		if (results.m <= 0.0){
			myfile2 << 0 << " " << 0 << " " << BaryonDensity << " " << pr4 << " " << pr5 << " " << pr6 << '\n';
			myfile2.close();
			cout << '\n' << "Unphysical Masses, Excluded" << '\n';
			return 0;
		}
	
//		cout << "M: " << results.m << ", R: " << results.r*1e-3 << ", rho0_c: " << BaryonDensity << 
//				", sigchin: " << pr4 << ", sigchi2: " << pr5 << ", DMmass: " << pr6 << '\n';

		BaryonDensity *= 1.1;
	}

	myfile2.close();
	cout << '\a' << '\n';
	t2 = clock();
	double tdiff = (t2 - t1)/CLOCKS_PER_SEC;
	if (tdiff <60){cout << "TIME: " << tdiff << '\n';}
	else if(tdiff>= 60){ tdiff*=1./60; cout << "TIME: " << tdiff << '\n';}
	cout << fno << '\n';
	return 0;
}
