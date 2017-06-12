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
	double PRESS[rows] = {0.0};
	double ENDENS[rows] = {0.0};

	for (int i = 1; i < rows; i++){
		PRESS[i] = EOSDATA[i][2];
		ENDENS[i] = EOSDATA[i][1];
	}

	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, n);
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

int tov_solver_solve(tov_results * results, double CentralDensity){

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
//		cout << p[i] << " " << e[i] << '\n';
	}

	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, rows);
	gsl_spline_init (spline, BDENS, PRESS, rows);
	double PRESSURE = gsl_spline_eval(spline, CentralDensity, acc);
//	cout << CentralDensity << " " << PRESSURE << " " << rows << '\n';

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

	while (y[0] > 1e-12){

		status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &r, rMax, &h, y);

		if (status != GSL_SUCCESS){
			cout << "ERROR: return value = " << status << '\n';
			break;
		}
	//	cout << "R: " << r*1e-3 << ",  M: " << y[1] << ",  P: " << y[0] << ",  E: " << EnergyAux << '\n';
//		cin >> dummyplug;
	}
	gsl_odeiv2_evolve_free (e);
	gsl_odeiv2_control_free (c);
	gsl_odeiv2_step_free (s);
	results->r = r;
	results->m = y[1];

/*==============================================================================
============================ SET UP INTEGRATOR TWO ==============================
* =============================================================================*/
/*
	gsl_odeiv2_step * s2 = gsl_odeiv2_step_alloc (T, 2);
	gsl_odeiv2_control * c2 = gsl_odeiv2_control_y_new (eps_abs, eps_rel);
	gsl_odeiv2_evolve * e2 = gsl_odeiv2_evolve_alloc (2);
 	sys.function = &tovinv_rhs;

	h *= -1;
	p = y[0];
	//cout << p << '\n';
	y[0] = r;
	while (p > 0.0){
		results->r = y[0];
		results->m = y[1];
		status = gsl_odeiv2_evolve_apply(e2, c2, s2, &sys, &p, 0.0, &h, y);

		if (status != GSL_SUCCESS){
			cout << "ERROR: return value = " << status << '\n';
			break;
		}
//		cout << "gwhynuno: " << p << " " << PressAux << " " << y[1] << " " << y[0] << " " << '\n';
//		cin >> kappa;
	}
	gsl_odeiv2_evolve_free (e2);
	gsl_odeiv2_control_free (c2);
	gsl_odeiv2_step_free (s2);
*/

	return GSL_SUCCESS;
}



int main(int argc, char* argv[]){
	int i;
	double BaryonDensity;
	double t1, t2;
	t1 = clock();
	tov_results results;
	int fno;// = argv[1];
//	cout << '\n' << '\n' << "FILENUMBER: "; cin >> fno;
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

	int title1 = atoi(argv[1]);
	int title2 = atoi(argv[2]);
	int title3 = atoi(argv[3]);
	int title4 = atoi(argv[4]);
	stringstream ss;
	ss << title1 << title2 << title3 << title4;
	string fnos = ss.str();

	string particle = "FERMIONS";
	stringstream liltitle1;
	liltitle1 <<"./" << particle << "/EoSFiles/"<<title4<<"GeV/EoS_SigOmDM_"<<fnos<<"p.dat";
	string	BIGTITLE1 = liltitle1.str();
	ifstream dater;
	dater.open(BIGTITLE1);

	if (!dater){
		cout << "not open: " << " " << BIGTITLE1 <<  '\n';
		return 1;
	}

	int columns = 0;

	while(!dater.eof()){
		dater >> pr1 >> pr2 >> pr3 >> pr4 >> pr5 >> pr6;
		if(pr2> pr2p and pr3>pr3p){
			vector<double> row;
			row.push_back(pr1);
			row.push_back(pr2);
			row.push_back(pr3);
			EOSDATA.push_back(row);
			columns = row.size();
			pr2p = pr2;
			pr3p = pr3;
//			cout << pr1 << " " << pr2 << " " << pr3 << '\n';
		}
		else{continue;}
	}

	dater.close();

	rows = EOSDATA.size();
	cout << "rows: " << rows << " columns: " << columns << '\n';

	cout << '\n' << "DATA COPYING DONE..." << '\n';


/*==================================================================*/

	stringstream liltitle2;
	liltitle2 << "./" << particle << "/MRFiles/" << title4 << "GeV/MR_" << fnos << "p.dat";
	string BIGTITLE2 = liltitle2.str();
	ofstream myfile2;
	myfile2.open(BIGTITLE2);
	int jk = 0;

	BaryonDensity = 0.01;
//	cout << "DENS: "; cin >> BaryonDensity; cout << '\n';
	while (BaryonDensity <= 1.4){
		tov_solver_solve(&results, BaryonDensity);

		if((results.r*1e-3 >=10.) and (results.r*1e-3<=50)){
			myfile2 << results.m << " " << results.r*1e-3 << " " << BaryonDensity << '\n';// << " " << pr4 << " " << pr5 << " " << pr6 << '\n';
		}
		else if (results.r*1e-3 <= 9){
			cout << "Radius lower than relevant range" << '\n';
			myfile2.close();
			return 0;
		}

		cout << "M: " << results.m << ", R: " << results.r*1e-3 << ", rho0_c: " << BaryonDensity << 
			", sigchin: " << pr4 << ", sigchi2: " << pr5 << ", DMmass: " << pr6 << '\n';
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
