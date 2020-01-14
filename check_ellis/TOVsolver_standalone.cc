//solves TOV equations for an EOS with the following format:
// [1] nB [fm^-3] [2] energy density [MeV/fm^3] [3] pressure [MeV/fm^3]
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
vector< vector<double> > M_OF_R_NUC;
vector< vector<double> > EOSDATA;
vector<double> RADIUSDATA;
vector< vector<double> > EOSDATAbare;

using namespace std;
double Msolar = 1.98e30;
double Msolar0 = Msolar*convkgtoGeV*1e-42;
double GNewtonPrime = GNewton0*1e10/convGeVtoinvcm*1e30;
double EnergyAux = 0.0;
double DarkMatterSphereRadius=0.0;
double PMAXBARE = 0.0;
double dummyplug;

//A few splines are needed:
//1 used for interpolating epsilon(p) in tov_rhs [DM EOS]
gsl_interp_accel *acc_EofP = gsl_interp_accel_alloc();
gsl_spline *spline_EofP;
//2 used for interpolating p(BDens) in tov_solver_solve [DM EOS]
gsl_interp_accel *acc_PofnB = gsl_interp_accel_alloc();
gsl_spline *spline_PofnB;


typedef struct{ double m; double r;} tov_results;
int tov_rhs(const double r, const double y[], double f[], void * params){
	if( r ==0){
		f[0] = 0.;
		f[1] = 0.;
	}
	else{
	double p = y[0];
        double m = y[1];
	vector<double> m_and_r;
        m_and_r.push_back(m); m_and_r.push_back(r);
        M_OF_R_NUC.push_back(m_and_r);
	double eps=0.0;
//        cout << p << " " << r << '\n';
	eps = gsl_spline_eval(spline_EofP,p,acc_EofP);
	f[0] = -(GNewtonPrime/(r*r))*(eps + p)*(m*Msolar0 + 4*pi*gsl_pow_3(r)*p)/(1.-2*m*Msolar0*GNewtonPrime/r);
	f[1] = 4*pi*gsl_pow_2(r)*eps/Msolar0;
	}
	return GSL_SUCCESS;
}

int tov_solver_solve(tov_results * results, double CentralDensity, double lowerlim){

	double r;                   /* Radius (independent variable) */
	double p;
	double m;
	int status;                 /* Status returned by GSL functions */

/*==============================================================================
============================ READ DENSITY ======================================
* =============================================================================*/
	//spline 2
	double PRESSURE = gsl_spline_eval(spline_PofnB, CentralDensity, acc_PofnB);
/*==============================================================================
============================ SET UP INTEGRATOR ONE ==============================
* =============================================================================*/

	double y[2];                /* ODE state */
	double eps_abs = 1e-18;   /* Absolute precision goal */
	double eps_rel = 1e-18;   /* Relative precision goal */
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
	double rMax = 1e20;
	//start integrating from the center, which always has DM
	y[0] = PRESSURE;
	y[1] = 0.0;

	int dcount = 0;
        //while the pressure is greater than the table's lower lim, keep integrating
	while (y[0] >= lowerlim){
		double yprev = y[0];
                //arguments are evolution function, control function, step function, system of equations, 
		//independent variable, upper bound on rMax, tolerance, step size, result
		status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &r, rMax, &h, y);
		
		//cout << y[0] << '\n';

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
	double Density = 0.0;
	double t1, t2;
	t1 = clock();
	tov_results results;
	int fno;// = argv[1];
	int rows = 0;	
	int rowsbare = 0;	
	int FILES = 0;
        //EOS file should hold 6 columns   [1]:Density [fm^-3]   [2]:Pressure [MeV fm^-3]   [3]:Energy Density [MeV fm^-3]
	//                                 [4]:SIGCHIN [cm^2]    [5]:SIGCHI2 [cm^2]         [6]:Mchi [GeV]
	double pr1 = 0;	double pr2 = 0;	double pr3 = 0;	double pr4 = 0;	double pr5 = 0;	double pr6 = 0;
        // only the first three are relevant for the MR sequence, so initialize these to be < 0 for monotonicity
	double pr1p = -2.; 
	double pr2p = -2.;
	double pr3p = -2.;

	string EOSFILE = argv[1];   //name of EOSFILE 
//	double lowlim = atof(argv[2]);  //tell the solve the lower lim in density

/*==============================================================================
=========================== READ IN AND FIT EoSs FILES USING INTERPOLATOR ======
* =============================================================================*/
	stringstream liltitle1;
	liltitle1 << EOSFILE;

	string	BIGTITLE1 = liltitle1.str();
	ifstream dater;
	dater.open(BIGTITLE1);
//        dater.ignore(1, '\n');

	if (!dater){
		cout << "Corresponding EoS File Not Found: " << " " << BIGTITLE1 <<  '\n';
		return 1;
	}

	int columns = 0;
	while(!dater.eof()){
		dater >> pr1 >> pr2 >> pr3 ;
//                cout << pr1 << " " << pr1p << " " << pr2 << " " << pr2p << '\n';
		if(pr2> 1.*pr2p and pr3>1.*pr3p){
//                        cout << pr1 << " " << pr2 << " " << pr3 << " " << pr4 << " " << pr5 << " " << pr6 << '\n';
			vector<double> row;
			row.push_back(pr1);
			row.push_back(pr2);
			row.push_back(pr3);
			EOSDATA.push_back(row);
			columns = row.size();
                        pr1p = pr1;
			pr2p = pr2;
			pr3p = pr3;
		}
		else{continue;}
	}

	dater.close();

	rows = EOSDATA.size();
	double DensityMax = EOSDATA[rows-2][0];

	int nE = 0;
	nE = rows -1;
	double PRESS[rows];
	double ENDENS[rows];
	double BDENS[rows];
	for (int i = 1; i < rows; i++){
		PRESS[i] = EOSDATA[i][2];
		ENDENS[i] = EOSDATA[i][1];
		BDENS[i] = EOSDATA[i][0];
	}

	//define spline 1
	spline_EofP = gsl_spline_alloc(gsl_interp_linear, nE);
	gsl_spline_init (spline_EofP, PRESS, ENDENS, nE);
	//spline 2
	spline_PofnB = gsl_spline_alloc(gsl_interp_linear, nE);
	gsl_spline_init (spline_PofnB, BDENS, PRESS, nE);
/*==============================================================================
============================ PREPARE OUTPUT ====================================
* =============================================================================*/
	stringstream liltitle2;
	liltitle2 << "./DMMR_" << EOSFILE;
	string BIGTITLE2 = liltitle2.str();
	ofstream myfile2;
	myfile2.open(BIGTITLE2);
	int jk = 0;

/*==============================================================================
======== INTEGRATE THE EoS FOR ALL CENTRAL DENSITIES TO GET A SEQUENCE =========
* =============================================================================*/
	Density = 0.001;

//	while (Density <= DensityMax){
//                double lowlim= gsl_spline_eval(spline_PofnB, Density, acc_PofnB)*1e-10;
//		tov_solver_solve(&results, Density, lowlim);
//
//		if((results.r*1e-3 >=0.0) and (results.r*1e-3<=200)){
//			myfile2 << results.m << " " << results.r*1e-3 << " " << Density << '\n';
//                	cout << Density << " " << results.m << " " << results.r*1e-3 << '\n';
//		}
//		if (results.m < 0.0){
//			myfile2 << 0 << " " << 0 << " " << Density << " " << pr4 << " " << pr5 << " " << pr6 << '\n';
//			myfile2.close();
//			return 0;
//		}
//		
//
//		Density *= 1.05;
//	}
        Density = atof(argv[2]); DensityMax = atof(argv[3]); double steps_in_density = atof(argv[4]);
	double delta_density = (DensityMax - Density)/steps_in_density;
	if (delta_density <=0.0){ delta_density = 10.0;}
	while (Density <= DensityMax){
                double lowlim= gsl_spline_eval(spline_PofnB, Density, acc_PofnB)*1e-10; //FIXME: may want to set this to 1e-10*central_pressure instead
//        	cout << "calling solver " << lowlim << '\n';
		tov_solver_solve(&results, Density, lowlim);
//        	cout << "solver done \n";

		if((results.r*1e-3 >=0.0) and (results.r*1e-3<=200)){
			myfile2 << results.m << " " << results.r*1e-3 << " " << Density << '\n';
                	cout << Density << " " << results.m << " " << results.r*1e-3 << '\n';
		}
		if (results.m < 0.0){
			myfile2 << 0 << " " << 0 << " " << Density << " " << pr4 << " " << pr5 << " " << pr6 << '\n';
                	cout << Density << " " << results.m << " " << results.r*1e-3 << '\n';
			myfile2.close();
			return 0;
		}
		

		Density += delta_density;
	}
//        cout << "free spline memories \n";
	//free spline memories
	//1
	gsl_spline_free (spline_EofP);
	gsl_interp_accel_free (acc_EofP);
	//2
	gsl_spline_free (spline_PofnB);
	gsl_interp_accel_free (acc_PofnB);
	myfile2.close();
	return 0;
}

