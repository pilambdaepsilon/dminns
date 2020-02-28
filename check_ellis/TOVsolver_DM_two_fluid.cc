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
#include "EOSDM_util.cc"
#include <time.h>
#define VERBOSE 0

using namespace std;
double Msolar = 1.98e30;
double Msolar0 = Msolar*convkgtoGeV*1e-42;
double GNewtonPrime = GNewton0*1e10/convGeVtoinvcm*1e30;

//A few splines are needed:
//1: used for interpolating epsilon(p) in tov_rhs [NUC EOS]
gsl_interp_accel *acc_EofP = gsl_interp_accel_alloc();
gsl_spline *spline_EofP;
//2: used for interpolating p(nB) in tov_solver_solve [NUC EOS]
gsl_interp_accel *acc_PofnB = gsl_interp_accel_alloc();
gsl_spline *spline_PofnB;
//3: used for interpolating nB(p) in tov_rhs [NUC EOS]
gsl_interp_accel *acc_nBofP = gsl_interp_accel_alloc();
gsl_spline *spline_nBofP;
//4: used for interpolating epsilon(p) in tov_rhs [DM EOS]
gsl_interp_accel *acc_EofP_DM = gsl_interp_accel_alloc();
gsl_spline *spline_EofP_DM;
//5: used for interpolating p(nB) in tov_solver_solve [DM EOS]
gsl_interp_accel *acc_PofnC_DM = gsl_interp_accel_alloc();
gsl_spline *spline_PofnC_DM;
//6: used for interpolating nChi(P) in loop over central density [DM EOS]
gsl_interp_accel *acc_nCofP_DM = gsl_interp_accel_alloc();
gsl_spline *spline_nCofP_DM;

int check_which_EOS = 0;
//Make a struct to hold the results of m(r) and r from profile
typedef struct{ double m; double r; double m_DM;} tov_results;
/*=======================================================
        _         
   _ __| |__  ___ 
  | '__| '_ \/ __|
  | |  | | | \__ \
  |_|  |_| |_|___/
=======================================================*/
//Define the rhs for the TOV equations - This uses the two fluid approach, where DM and baryons couple only through a common gravitational potential
int tov_rhs(const double r, const double y[], double f[], void * params){
        //enforce central BC: all masses are zero at the origin to avoid singularities
	if(r ==0){
		f[0] = 0.;
		f[1] = 0.;
		f[2] = 0.;
		f[3] = 0.;
	}
	//otherwise, set the pressure to the central one, specified in function call in main(), and integrate
	else{
	double p    = y[0];									//the pressure is one of the integration variables
        double m    = y[1];									//the mass is the other integration variable
	double p_DM = y[2];									//the pressure is one of the integration variables
	double m_DM = y[3];									//the pressure is one of the integration variables
        //use only physical pressures, floor it at zero
	if(p < 0.0){
		p = 0.0;
	}
	if(p_DM < 0.0){
		p_DM = 0.0;
	}
	//get the central number & energy densities for each component
	double nB = gsl_spline_eval(spline_nBofP,p,acc_nBofP);					//get the baryon mass density using a spline
	double nChi = gsl_spline_eval(spline_nCofP_DM,p_DM,acc_nCofP_DM);			//get the baryon mass density using a spline
	double eps    = 0.0;									//initialize the energy density (baryons)
	double eps_DM = 0.0;									//initialize the energy density (DM)
	eps = gsl_spline_eval(spline_EofP,p,acc_EofP);
	eps_DM = gsl_spline_eval(spline_EofP_DM,p_DM,acc_EofP_DM);
	//make the rhs of the TOV equations, coupled through the second equation (common gravitational potential due to total mass)
	//in two fluid approach, the ode for total pressure decouples to two separate equations (first and third below) corresponding to
	//pressure for each sector. But the total mass equation remains the same, with contribution from both sectors
	f[0] = -(GNewtonPrime/(r*r))*(eps + p)*(m*Msolar0 + 4*pi*gsl_pow_3(r)*(p+p_DM))/(1.-2*m*Msolar0*GNewtonPrime/r);	//baryon pressure, through total mass
	f[1] = 4*pi*gsl_pow_2(r)*(eps + eps_DM)/Msolar0;									//total mass, including both sectors
	f[2] = -(GNewtonPrime/(r*r))*(eps_DM + p_DM)*(m*Msolar0 + 4*pi*gsl_pow_3(r)*(p_DM+p))/(1.-2*m*Msolar0*GNewtonPrime/r);	//DM pressure, through total mass
	f[3] = 4*pi*gsl_pow_2(r)*(eps_DM )/Msolar0;										//DM mass
	}
	return GSL_SUCCESS;
}
/*=======================================================
   _       _                       _             
  (_)_ __ | |_ ___  __ _ _ __ __ _| |_ ___  _ __ 
  | | '_ \| __/ _ \/ _` | '__/ _` | __/ _ \| '__|
  | | | | | ||  __/ (_| | | | (_| | || (_) | |   
  |_|_| |_|\__\___|\__, |_|  \__,_|\__\___/|_|   
                   |___/                   
=======================================================*/
//Make a callable function which integrates the above rhs for a given EOS for a given central number density
int tov_solver_solve(tov_results * results, double CentralDensity, double CentralDMDensity, double lowerlim){
	
	if(VERBOSE ==1){
		cout << "STARTING INTEGRATOR \n";
	}

	//this uses the (decoupled) two fluid approach. We integrate the TOV equations separately for DM and SM components
	// with a commom gravitational potential
	double r;          								        /* Radius (independent variable) */
	double p;								               	/* Pressure (integration variable) */
	double m;                   								/* Mass (integration variable) */
	double m_DM;                   								/* Mass (integration variable) */
	int status;                								/* Status returned by GSL functions */

	/*==============================================================================
	============================ GET CENTRAL PRESSURES =============================
	* =============================================================================*/
	double PRESSURE = 0.0; double DMPRESSURE=0.0;						//initialize the pressures
	PRESSURE = gsl_spline_eval(spline_PofnB, CentralDensity, acc_PofnB);		//spline #2
	DMPRESSURE = gsl_spline_eval(spline_PofnC_DM, CentralDMDensity, acc_PofnC_DM);	//spline #5
//	cout << "PRESSURE: " << PRESSURE << "--- DMPRESSURE: " << DMPRESSURE << '\n';
	/*==============================================================================
	============================ SET UP THE INTEGRATOR =============================
	* =============================================================================*/
	double y[4];            								/* ODE variables (p and m) are held here: p,m,p_DM,m_DM*/
	double eps_abs = 1e-18;  								/* Absolute precision goal */
	double eps_rel = 1e-18;									/* Relative precision goal */
	const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk8pd;					//use RK8 which is probably overkill
	double h = 1e-6;									//step size in ode solver

	gsl_odeiv2_step * s = gsl_odeiv2_step_alloc (T, 4);					//allocate memory to the step
	gsl_odeiv2_control * c = gsl_odeiv2_control_y_new (eps_abs, eps_rel);			//set the preferred tolerances
	gsl_odeiv2_evolve * e = gsl_odeiv2_evolve_alloc (4);					//allocate memory to the evolution

	//configure the ode system
	if(VERBOSE ==1){
		cout << "CONFIGURING ODES \n";
	}
	gsl_odeiv2_system sys;									
	sys.params = NULL;								
	sys.dimension = 4;									//two dimensional system
	sys.jacobian = NULL;									//no coordinate transformations necessary
 	sys.function = &tov_rhs;								//the rhs specified above
	//set up boundary conditions
	r = 0.0;										//start integrating at r = 0
	double rMax = 1e20;									//don't go beyond this radius
	y[0] = PRESSURE;									//set central p (y[0]) and m (y[1])
	y[1] = 0.0;
	y[2] = DMPRESSURE;
	y[3] = 0.0;

        //while the (baryon) pressure is greater than the lower limit (specified in loop over CentralDensity), keep integrating
	if(VERBOSE ==1){
		cout << "INTEGRATING \n";
	}
	while (y[0] >= lowerlim){
                //arguments are evolution function, control function, step function, system of equations, 
		//independent variable, upper bound on rMax, tolerance, step size, result
		status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &r, rMax, &h, y);
		
		//throw the error in case of failure
		if (status != GSL_SUCCESS){
			cout << "ERROR: return value = " << status << '\n';
			break;
		}
	}
	//free the integration system, ready for next round of integration
	gsl_odeiv2_evolve_free (e);
	gsl_odeiv2_control_free (c);
	gsl_odeiv2_step_free (s);
	//it puts the lotion in the basket
	results->r = r;
	results->m = y[1];
	results->m_DM = y[3];
	if(VERBOSE ==1){
		cout << "DONE - GATHERING RESULTS \n";
	}

	return GSL_SUCCESS;
}



/*=======================================================
 _ __ ___   __ _(_)_ __  
| '_ ` _ \ / _` | | '_ \
| | | | | | (_| | | | | |
|_| |_| |_|\__,_|_|_| |_|
=======================================================*/
int main(int argc, char* argv[]){
	if(VERBOSE ==1){
		cout << "STARTING \n";
	}
	/*==============================================================================
	============================ NAMING SCHEME AND VARIABLES =======================
	* =============================================================================*/
	double Density = 0.0;									//Controls the central boundary condition -baryons
	double DMDensity = 0.0;									//Controls the central boundary condition -DM
	tov_results results;									//intialize struct
	int rows = 0;										//counts rows 
	int columns = 0;									//counts columns
        //EOS file should hold 6 columns   [1]:Density [fm^-3]   [2]:Pressure [MeV fm^-3]   [3]:Energy Density [MeV fm^-3]
	//                                 [4]:SIGCHIN [cm^2]    [5]:SIGCHI2 [cm^2]         [6]:Mchi [GeV]
	//name of EOSFILE
	string EOSFILEnuc = argv[1];  
	string EOSFILEdm = argv[6];  

	/*==============================================================================
	=========================== READ IN AND FIT EoS FILES USING INTERPOLATOR =======
	* =============================================================================*/
	if(VERBOSE ==1){ cout << "READING EOS \n";}
	bool PRINTEOS = 0;
	rows = readEOS(EOSFILEnuc, PRINTEOS);
	//set the max density somewhere near the table's upper limit
	double DensityMax = EOSDATA[rows-2][0];
	
	//dummy: number of indices used in the interpolation
	double PRESS[rows]; double ENDENS[rows]; double BDENS[rows];
	//copy the EOS into arrays, GSL doesn't like vectors
	for (int i = 0; i < rows; i++){
		PRESS[i] = EOSDATA[i][2]; ENDENS[i] = EOSDATA[i][1]; BDENS[i] = EOSDATA[i][0];
	}

	if(VERBOSE ==1){ cout << "FITTING NUC EOS \n";}
	int nE = rows-1;
	//define spline 1 - eps(p) for nuclear EOS
	spline_EofP = gsl_spline_alloc(gsl_interp_linear, nE);
	gsl_spline_init (spline_EofP, PRESS, ENDENS, nE);
	//spline 2 - p(nB) for nuclear EOS
	spline_PofnB = gsl_spline_alloc(gsl_interp_linear, nE);
	gsl_spline_init (spline_PofnB, BDENS, PRESS, nE);
	//spline 3 - nB(p) for nuclear EOS
	spline_nBofP = gsl_spline_alloc(gsl_interp_linear, nE);
	gsl_spline_init (spline_nBofP, PRESS, BDENS, nE);
	double nBtrial = gsl_spline_eval(spline_nBofP,PRESS[0],acc_nBofP);					//get the baryon mass density using a spline
	cout << "ROWS NUC: " << rows << '\n';

	/*==============================================================================
	=========================== READ IN AND FIT DM EoS =============================
	* =============================================================================*/

	EOSDATA.clear();
	if(VERBOSE ==1){ cout << "READING DM EOS \n";}
	bool PRINTDMEOS = 0;
	rows = readEOS(EOSFILEdm, PRINTDMEOS);

	//set the max density somewhere near the table's upper limit
	double DMDensityMax = EOSDATA[rows-2][0];
	
	//dummy: number of indices used in the interpolation
	double DMPRESS[rows]; double DMENDENS[rows]; double DMDENS[rows];
	//copy the EOS into arrays, GSL doesn't like vectors
	for (int i = 0; i < rows; i++){
		DMPRESS[i] = EOSDATA[i][2]; DMENDENS[i] = EOSDATA[i][1]; DMDENS[i] = EOSDATA[i][0];
	}


	if(VERBOSE ==1){ cout << "FITTING DM EOS \n";}
	nE = rows-1;

	//define spline 4 - eps(p) for DM EOS
	spline_EofP_DM = gsl_spline_alloc(gsl_interp_linear, nE);
	gsl_spline_init (spline_EofP_DM, DMPRESS, DMENDENS, nE);
	//spline 5 - p(nC) for DM EOS
	spline_PofnC_DM = gsl_spline_alloc(gsl_interp_linear, nE);
	gsl_spline_init (spline_PofnC_DM, DMDENS, DMPRESS, nE);
	//spline 6 - nC(p) for nuclear EOS
	spline_nCofP_DM = gsl_spline_alloc(gsl_interp_linear, nE);
	gsl_spline_init (spline_nCofP_DM, DMPRESS, DMDENS, nE);

	EOSDATA.clear();
	rows = readEOS(EOSFILEdm, PRINTEOS);
	cout << "ROWS DM: " << rows << '\n';

	//==============================================================================
	//============================ PREPARE OUTPUT ====================================
	// =============================================================================
	stringstream liltitle2;
	liltitle2 << "./DMMR_" << EOSFILEnuc;
	string BIGTITLE2 = liltitle2.str();
	ofstream myfile2;
	myfile2.open(BIGTITLE2);

	//==============================================================================
	//======== INTEGRATE THE EoS FOR ALL CENTRAL DENSITIES TO GET A SEQUENCE =========
	// =============================================================================
        Density = atof(argv[2]); DensityMax = atof(argv[3]); double steps_in_density = atof(argv[4]);
//	DMDensity = atof(argv[5]);   						//Central Baryon Density
	double DMFraction = atof(argv[5]);   					//Central DM Density (fraction of baryon density)
//	double DMmass = atof(argv[6]);   							//DMmass
//	double eta = atof(argv[7]);  								//Coupling to mediator mass ratio (squared)
	double delta_density = (DensityMax - Density)/steps_in_density;
	if (delta_density <=0.0){ delta_density = 10.0;}
	if(VERBOSE ==1){
		cout << "STARTING LOOP \n";
	}
	while (Density <= DensityMax){
		DMDensity = Density*DMFraction;
                double lowlim= 0.0;
		lowlim = gsl_spline_eval(spline_PofnB, Density, acc_PofnB)*1e-10;					//the lower limit in the pressure is p_atm = 10^-10 p_c
															//(10 orders of magnitude lower than the central)
		//==========================================================
		//======== INTEGRATE THE USING TWO-FLUID APPROACH  =========
		//==========================================================
		//M_OF_R_NUC is a global vector which holds the profiles
//		cout << "Density: " << Density << " --- DMDensity: "<< DMDensity << '\n';
		tov_solver_solve(&results, Density, DMDensity, lowlim);
		//keep only relevant models
		if((results.r*1e-3 >=0.0) and (results.r*1e-3<=200)){
			myfile2 << results.m << " " << results.r*1e-3 << " " << Density << '\n';
//                	cout << "nB M R M_DM:  " << Density << " " << results.m << " " << results.r*1e-3 << " " << results.m_DM << '\n';
                	cout << Density << " " << results.m << " " << results.r*1e-3 << " " << results.m_DM << '\n';
		}
		//floor unphysical models to zero mass and radius (used to debug - if MR relation goes to origin, then unphysical models built)
		if (results.m < 0.0){
			myfile2 << 0 << " " << 0 << " " << Density << '\n';
			myfile2.close();
			return 0;
		}

		Density += delta_density;
	}
	//free spline memories
	//1
	gsl_spline_free (spline_EofP);
	gsl_interp_accel_free (acc_EofP);
	//2
	gsl_spline_free (spline_PofnB);
	gsl_interp_accel_free (acc_PofnB);
	//3
	gsl_spline_free (spline_nBofP);
	gsl_interp_accel_free (acc_nBofP);
	//4
	gsl_spline_free (spline_EofP_DM);
	gsl_interp_accel_free (acc_EofP_DM);
	//5
	gsl_spline_free (spline_PofnC_DM);
	gsl_interp_accel_free (acc_PofnC_DM);
	//6
	gsl_spline_free (spline_nCofP_DM);
	gsl_interp_accel_free (acc_nCofP_DM);

	myfile2.close();
	return 0;
}

