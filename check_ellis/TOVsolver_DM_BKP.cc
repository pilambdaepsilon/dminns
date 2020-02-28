//solves TOV equations for an EOS with the following format:
// [1] nB [fm^-3] [2] energy density [MeV/fm^3] [3] pressure [MeV/fm^3]
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
//#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_roots.h>
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
#define VERBOSE 0

vector< vector<double> > M_OF_R_NUC;
vector< vector<double> > EOSDATA;

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
//4: used for interpolating nChi(nB) in loop over central density [DM EOS]
gsl_interp_accel *acc_nChiofnB = gsl_interp_accel_alloc();
gsl_spline *spline_nChiofnB;
//5: used for interpolating epsilon(p) in tov_rhs [DM EOS]
gsl_interp_accel *acc_EofP_DM = gsl_interp_accel_alloc();
gsl_spline *spline_EofP_DM;
//6: used for interpolating p(nB) in tov_solver_solve [DM EOS]
gsl_interp_accel *acc_PofnB_DM = gsl_interp_accel_alloc();
gsl_spline *spline_PofnB_DM;

int check_which_EOS = 0;
//          _..._
//        .'     '.
//       / \     / \
//      (  |     |  )
//      (`"`  "  `"`)
//       \         /
//        \  ___  /
//         '.___.'
struct Rparams{
	double dm_mass;
	double dm_coupling_eta;
	double central_dm_density;
	double schwarzschild_component;
};
double DM_DENSITY_PROFILE(double NX, void *params){
	//SET PARAMETERS
	double dm_mass = ((struct Rparams *) params)->dm_mass;
	double dm_coupling_eta = ((struct Rparams *) params)->dm_coupling_eta;
	double central_dm_density = ((struct Rparams *) params)->central_dm_density;
	double schwarzschild_component = ((struct Rparams *) params)->schwarzschild_component;
	double central_dm_chempot = sqrt(dm_mass*dm_mass + pow(3*pi*pi * central_dm_density, 2./3));
//	cout << "PARAMS: " << dm_mass << " " << dm_coupling_eta << " " << central_dm_density << " " << schwarzschild_component << '\n';
	//SET VARIABLES TO BE SOLVED FOR
	//SET FUNCTIONS TO BE SOLVED
	double f0 = (sqrt(dm_mass*dm_mass + pow(3*pi*pi* NX,2./3)) + dm_coupling_eta * NX - central_dm_chempot*sqrt(schwarzschild_component));
	
	return f0;
}
//Make a struct to hold the results of m(r) and r from profile
typedef struct{ double m; double r;} tov_results;
/*=======================================================
        _         
   _ __| |__  ___ 
  | '__| '_ \/ __|
  | |  | | | \__ \
  |_|  |_| |_|___/
=======================================================*/
//Define the rhs for the TOV equations
int tov_rhs(const double r, const double y[], double f[], void * params){
        int DMFLAG = *(int *)params;								//a flag that decides which EOS to use in integration
	if( r ==0){
		f[0] = 0.;
		f[1] = 0.;
		f[2] = 0.;
	}
	else{
	double p    = y[0];									//the pressure is one of the integration variables
        double m    = y[1];									//the mass is the other integration variable
	double p_DM = y[2];									//the pressure is one of the integration variables
//	cout << "prhs: " << p << " " << p_DM << '\n';
	if(p < 0.0){
		p = 0.0;
	}
	if(p_DM < 0.0){
		p_DM = 0.0;
	}

	double nB = gsl_spline_eval(spline_nBofP,p,acc_nBofP);					//get the baryon mass density using a spline
	vector<double> m_and_r;									//make a vector in which to temporarily store m(r) and r
        m_and_r.push_back(r); m_and_r.push_back(m); m_and_r.push_back(nB);			//push them back into the vector to extract profile
	m_and_r.push_back(p); 
        M_OF_R_NUC.push_back(m_and_r);								//push back the entire profile set into a global vector of vectors
	double eps=0.0;										//initialize the energy density
	//if DM flag is on, then use the DM (+nuclear) EOS to specify the energy density
        if(DMFLAG==1){
//		cout << eps << " " << p << '\n';
		if(VERBOSE ==1 and check_which_EOS == 0){
			cout << "DM EOS \n";check_which_EOS ++;
		}
		eps = gsl_spline_eval(spline_EofP,p,acc_EofP);
//		eps = gsl_spline_eval(spline_EofP_DM,p,acc_EofP_DM);
	}
	//otherwise use the nuclear EOS
	else{
		if(VERBOSE ==1 and check_which_EOS == 0){
			cout << "NUC EOS \n";check_which_EOS ++;
		}
		eps = gsl_spline_eval(spline_EofP,p,acc_EofP);
	}
	//make the rhs
	f[0] = -(GNewtonPrime/(r*r))*(eps + p)*(m*Msolar0 + 4*pi*gsl_pow_3(r)*p)/(1.-2*m*Msolar0*GNewtonPrime/r);
	f[1] = 4*pi*gsl_pow_2(r)*eps/Msolar0;
//	f[1] = 4*pi*gsl_pow_2(r)*(eps_B + eps_DM)/Msolar0;
//	f[2] = -(GNewtonPrime/(r*r))*(epsi_DM + p_DM)*(m*Msolar0 + 4*pi*gsl_pow_3(r)*p_DM)/(1.-2*m*Msolar0*GNewtonPrime/r);
	f[2] = 0.0;
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
int tov_solver_solve(tov_results * results, double CentralDensity, double lowerlim, int use_DMEOS){
	
	if(VERBOSE ==1){
		cout << "STARTING INTEGRATOR \n";
	}

	double r;          								        /* Radius (independent variable) */
	double p;								               	/* Pressure (integration variable) */
	double m;                   								/* Mass (integration variable) */
	int status;                								/* Status returned by GSL functions */

	/*==============================================================================
	============================ GET CENTRAL PRESSURE ==============================
	* =============================================================================*/
	double PRESSURE = 0.0;									//initialize the pressure
	if(use_DMEOS ==1){
		PRESSURE = gsl_spline_eval(spline_PofnB_DM, CentralDensity, acc_PofnB_DM);		//spline #2
		if(VERBOSE ==1){
			cout << "GETTING PRESSURE DM \n";
		}
	}
	else{
		PRESSURE = gsl_spline_eval(spline_PofnB, CentralDensity, acc_PofnB);		//spline #2
		if(VERBOSE ==1){
			cout << "GETTING PRESSURE NUC\n";
		}
	}
	/*==============================================================================
	============================ SET UP THE INTEGRATOR =============================
	* =============================================================================*/
	double y[3];            								/* ODE variables (p and m) are held here */
	double eps_abs = 1e-18;  								/* Absolute precision goal */
	double eps_rel = 1e-18;									/* Relative precision goal */
	const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk8pd;					//use RK8 which is probably overkill
	double h = 1e-8;									//step size in ode solver

	gsl_odeiv2_step * s = gsl_odeiv2_step_alloc (T, 3);					//allocate memory to the step
	gsl_odeiv2_control * c = gsl_odeiv2_control_y_new (eps_abs, eps_rel);			//set the preferred tolerances
	gsl_odeiv2_evolve * e = gsl_odeiv2_evolve_alloc (3);					//allocate memory to the evolution

	//configure the ode system
	if(VERBOSE ==1){
		cout << "CONFIGURING ODES \n";
	}
	gsl_odeiv2_system sys;									
	sys.params = &use_DMEOS;								//the only parameter is the DM EOS flag
	sys.dimension = 3;									//two dimensional system
	sys.jacobian = NULL;									//no coordinate transformations necessary
 	sys.function = &tov_rhs;								//the rhs specified above
	//set up boundary conditions
	r = 0.0;										//start integrating at r = 0
	double rMax = 1e20;									//don't go beyond this radius
	y[0] = PRESSURE;									//set central p (y[0]) and m (y[1])
	y[1] = 0.0;
	y[2] = 0.0;

        //while the pressure is greater than the lower limit (specified in loop over CentralDensity), keep integrating
	if(VERBOSE ==1){
		cout << "INTEGRATING \n";
	}
	while (y[0] >= lowerlim){
	//	cout << "pint: " << y[0] << " ";
		double yprev = y[0];
                //arguments are evolution function, control function, step function, system of equations, 
		//independent variable, upper bound on rMax, tolerance, step size, result
		status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &r, rMax, &h, y);
		
	//	cout << "pint2: " << y[0] << "\n";
		//throw the error
		if (status != GSL_SUCCESS){
			cout << "ERROR: return value = " << status << '\n';
			break;
		}
	}
	//free the integration system, ready for next round of integration
	gsl_odeiv2_evolve_free (e);
	gsl_odeiv2_control_free (c);
	gsl_odeiv2_step_free (s);
	//it puts the (lotion)results in the (basket)struct
	results->r = r;
	results->m = y[1];
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
	double Density = 0.0;									//Controls the central boundary condition
	tov_results results;									//intialize struct
	int rows = 0;										//counts rows 
	int columns = 0;									//counts columns
        //EOS file should hold 6 columns   [1]:Density [fm^-3]   [2]:Pressure [MeV fm^-3]   [3]:Energy Density [MeV fm^-3]
	//                                 [4]:SIGCHIN [cm^2]    [5]:SIGCHI2 [cm^2]         [6]:Mchi [GeV]
	double pr1 = 0.0; double pr2 = 0.0; double pr3 = 0.0; double pr4 = 0.0; double pr5 = 0.0; double pr6 = 0.0;
        // only the first three are relevant for the MR sequence, so initialize these to be < 0 for monotonicity
	double pr1p = -2.; double pr2p = -2.; double pr3p = -2.;

	string EOSFILE = argv[1];   								//name of EOSFILE 
	double CentralDMDensity = atof(argv[2]);   						//CentralDMDensity
	double DMmass = atof(argv[3]);   							//DMmass
	double eta = atof(argv[4]);  								//Coupling to mediator mass ratio (squared)

	/*==============================================================================
	=========================== READ IN AND FIT EoSs FILES USING INTERPOLATOR ======
	* =============================================================================*/
	if(VERBOSE ==1){
		cout << "READING EOS \n";
	}
	stringstream file_name_eos;
	file_name_eos << EOSFILE;
	string	MAINTITLE = file_name_eos.str();
	ifstream dater;										//input stream called dater
	dater.open(MAINTITLE);									//open the input for reading

	if (!dater){
		cout << "Corresponding EoS File Not Found: " << " " << MAINTITLE <<  '\n';
		return 1;
	}
	//while the input stream is open, read three parameters per line
	vector<double> rowvector; rowvector.push_back(pr1); rowvector.push_back(pr2); rowvector.push_back(pr3);
	EOSDATA.push_back(rowvector);
	rowvector.clear();
	pr1p = pr1; pr2p = pr2; pr3p = pr3;
	double trial_shift = 0.0;
	while(!dater.eof()){
		dater >> pr1 >> pr2 >> pr3 ;
		//ensure monotonicity
		if(pr2> 1.1*pr2p and pr3>1.1*pr3p){
//			cout << pr1 << " " << pr2 << " " << pr3 << '\n';
			//make a vector for the row
			//push each value into the vector row
			rowvector.push_back(pr1);	
			rowvector.push_back(pr2 + trial_shift);
			rowvector.push_back(pr3 +trial_shift);
			//push the entire vector into the vector of vectors for the EOS
			EOSDATA.push_back(rowvector);
			//count the columns
			columns = rowvector.size();
			//update the monotonicity dummy variables
                        pr1p = pr1;
			pr2p = pr2;
			pr3p = pr3;
			//clear the row vector
			rowvector.clear();
		}
		else{continue;}
	}
	//close the input stream
	dater.close();
	//count the rows
	rows = EOSDATA.size();
	//set the max density somewhere near the table's upper limit
	double DensityMax = EOSDATA[rows-2][0];
	
	//dummy: number of indices used in the interpolation
	int nE = rows-1;
	double PRESS[rows];
	double ENDENS[rows];
	double BDENS[rows];
	//copy the EOS into arrays, GSL doesn't like vectors
	for (int i = 0; i < rows; i++){
		PRESS[i] = EOSDATA[i][2];
		ENDENS[i] = EOSDATA[i][1];
		BDENS[i] = EOSDATA[i][0];
//		cout << i << " " << BDENS[i] << " " << ENDENS[i] << " " << PRESS[i] << '\n';
	}

	if(VERBOSE ==1){
		cout << "FITTING EOS \n";
	}
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
//	cout << "NBTRIAL: " << nBtrial << " " << BDENS[0] <<  " " << Density << " " << DensityMax << '\n';
	//==============================================================================
	//============================ PREPARE OUTPUT ====================================
	// =============================================================================
	stringstream liltitle2;
	liltitle2 << "./DMMR_" << EOSFILE;
	string BIGTITLE2 = liltitle2.str();
	ofstream myfile2;
	myfile2.open(BIGTITLE2);

	//==============================================================================
	//======== INTEGRATE THE EoS FOR ALL CENTRAL DENSITIES TO GET A SEQUENCE =========
	// =============================================================================
//	Density = 0.001; 
//	DensityMax = 1.19;

        Density = 0.55; DensityMax = Density + 0.001;
	if(VERBOSE ==1){
		cout << "STARTING LOOP \n";
	}
	while (Density <= DensityMax){
                double lowlim= 0.0;
                double kFermiDM = pow((3*pi*pi *CentralDMDensity*1e39*pow(convinvcmtoGeV,3.0)), (1./3.));	//the Fermi momentum depends on the (specified) central density in GeV
		lowlim = gsl_spline_eval(spline_PofnB, Density, acc_PofnB)*1e-10;					//the lower limit in the pressure is p_atm = 10^-10 p_c
//		cout << "P: " << PRESS[0] << '\n';
        //        double lowlim= PRESS[0];					//the lower limit in the pressure is p_atm = 10^-10 p_c
		//===========================================
		//======== INTEGRATE THE NUCLEAR EoS  =========
		//==========================================
		//M_OF_R_NUC is a global vector which holds the profiles
                int use_DMEOS = 0;
		tov_solver_solve(&results, Density, lowlim, use_DMEOS);
		//print
		if((results.r*1e-3 >=0.0) and (results.r*1e-3<=200)){
			myfile2 << results.m << " " << results.r*1e-3 << " " << Density << '\n';
	               	cout << "NUC: " << Density << " " << results.m << " " << results.r*1e-3 << '\n';
		}

		int array_size = M_OF_R_NUC.size();									//get the size of M_OF_R_NUC
		vector <double> NCHI_OF_R; vector <double> NB_OF_R; vector <double> R;					//Make new vectors to extract the profiles
		double nchidumb = 1000.0; double nbdumb = 1000.0;							//make dummy variables to ensure monotonicity
		//===========================================
		//==== GET RADIAL PROFILES OF NB, NCHI  =======
		// ==========================================
                for(int i = 0; i <  array_size; i++){
			double gtt = 1.0 -  (2.0*GNewton/(clight*clight)) * M_OF_R_NUC[i][1]*Msolar/M_OF_R_NUC[i][0];    // 00 component of metric in mks (actually unitless)
                        double DMchempot = sqrt(kFermiDM*kFermiDM + DMmass*DMmass);					//the DM chemical potential in GeV
			double nb_aux = M_OF_R_NUC[i][2];
                        double nchiofr = DMchempot * (1./sqrt(gtt)) * (1./eta) * pow(convGeVtoinvcm,3.0)*1e-39;
//======================= SET NCHIOFR USING NEWTON-RAPHSON =============================
			const gsl_root_fsolver_type *T;		//DECLARE VARIABLE FOR TYPE OF SOLVER
			gsl_root_fsolver *S;		//DECLARE NAME FOR SOLVER WORKSPACE
			int status;				//USED TO UPDATE ON STATUS (IS ACTUALLY BOOLEAN)
			size_t iter =0;
			const size_t n = 2;			//DIMENSIONALITY OF PROBLEM
			struct Rparams p = {DMmass, eta, CentralDMDensity*1e39*pow(convinvcmtoGeV,3.0), gtt};		//INITIALIZE PARAMETERS
			gsl_function func; func.function = &DM_DENSITY_PROFILE; func.params=&p;	//MAKE A FUNCTION FOR SOLVER OUT OF ROSENBROCK, OF DIMENSION n, WITH PARAMS p

			double x_lo = 0.0, x_hi = 5.0;
	
			T = gsl_root_fsolver_brent;	//MAKE THE TYPE OF SOLVER A HYBRID S SOLVER
			S = gsl_root_fsolver_alloc(T);	//MAKE THE WORKSPACE FOR SOLVER TYPE T OF DIMENSION 2
			gsl_root_fsolver_set (S, &func, x_lo, x_hi);

			double rootie_tootie = 0;
			while(status = GSL_CONTINUE && iter < 1000){
				iter ++;
//				print_state(iter, S);				//PRINT THE CURRENT STATUS
			        status = gsl_root_fsolver_iterate(S);
				rootie_tootie = gsl_root_fsolver_root(S);
				x_lo = gsl_root_fsolver_x_lower (S);
				x_hi = gsl_root_fsolver_x_upper (S);
				status = gsl_root_test_interval (x_lo, x_hi, 0, 0.001);
//				printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n", iter, x_lo, x_hi, r, r - CentralDMDensity, x_hi - x_lo);
	
			}

//			cout << "status = " << gsl_strerror(status) << " " << rootie_tootie << '\n';
			nchiofr = rootie_tootie/(1e39*pow(convinvcmtoGeV,3.0));
			gsl_root_fsolver_free(S);			//CLEAR MEMORY
//======================= SET NCHIOFR USING NEWTON-RAPHSON =============================
			//ensure monotonicity in the baryon density
			//if(1.1*M_OF_R_NUC[i][2] < nbdumb && nb_aux > 0.0){
			if(1.1*M_OF_R_NUC[i][2] < nbdumb){
				//push the profiles into the vectors
				NB_OF_R.push_back(nb_aux); NCHI_OF_R.push_back(nchiofr); R.push_back(M_OF_R_NUC[i][0]*1e-3);
				//print: Radius                          Mass                     nB                         nX                pressure                    gtt 
				cout << M_OF_R_NUC[i][0]*1e-3 << " " << M_OF_R_NUC[i][1] << " " << M_OF_R_NUC[i][2] << " " << nchiofr << " " << M_OF_R_NUC[i][3] << " " << 1.0/sqrt(1.0 - (2.0*GNewton/(clight*clight)) * M_OF_R_NUC[i][1]*Msolar/M_OF_R_NUC[i][0]) << '\n';
				nbdumb = nb_aux;
			}
		}

		//================================================
		//===== INTERPOLATE NCHI AS A FUNCTION OF NB =======
		//===============================================
		//copy the vectors NB_OF_R and NCHI_OF_R into arrays - GSL
		int size_of_nB_array = NB_OF_R.size();									//not the size of the original EOS table
		vector <double> NCHIMONO; vector <double> NBMONO;
		double nchiprev =-1.0;
		for(int i = size_of_nB_array-1; i >= 0; i--){
			if(NCHI_OF_R[i] > nchiprev){
				NCHIMONO.push_back(NCHI_OF_R[i]);
				NBMONO.push_back(NB_OF_R[i]);
				nchiprev = NCHI_OF_R[i];
//				cout << NB_OF_R[i] << " " << NCHI_OF_R[i] << '\n';
			}
		}
		int size_of_monotonized_nB_array = NBMONO.size();
		double NBARR[size_of_monotonized_nB_array]; double NCHIARR[size_of_monotonized_nB_array]; double RARR[size_of_monotonized_nB_array];	//TODO: that's probably not the safest
		//index as increasing baryon density instead
		for(int i = 0; i < size_of_monotonized_nB_array; i++){
			NBARR[i] = NBMONO[i];
			NCHIARR[i] = NCHIMONO[i];
//			RARR[(size_of_monotonized_nB_array-1)-i] = R[i];
		}
//		NBARR[0] = BDENS[0]; NCHIARR[0] = BDENS[0];
		for(int i=0; i < size_of_monotonized_nB_array; i++){
//			cout << i << " " << NBARR[i] << " " << NCHIARR[i] << '\n';
		}

		//spline 4 - nchi(nB) which couples DM density to baryon density
		spline_nChiofnB = gsl_spline_alloc(gsl_interp_linear, size_of_monotonized_nB_array);
		gsl_spline_init (spline_nChiofnB, NBARR, NCHIARR, size_of_monotonized_nB_array);
                M_OF_R_NUC.clear();											//clear the big profile vector - not needed anymore
		//=====================================================
		//===== CREATE A NEW EOS - MODIFIED BY DM TERMS =======
		// ====================================================
//
//		vector <double> PDM; vector <double> EDM; vector <double> NBDM;						//make vectors in which to hold the new EOS
//                double DensityAux = NBARR[1];										//a dummy variable for baryon density
//		double del_DensityAux = pow((Density/DensityAux), 1./10.0);
////		cout << "DAUX: " << DensityAux << " " << del_DensityAux << " " << BDENS[0] << '\n';
//		double P_nuc = 0.0; double E_nuc = 0.0; double Ptot = 0.0; double Etot = 0.0; double DMDensity = 0.0;
//		double P_kinDM = 0.0; double E_kinDM = 0.0;
//                PDM.push_back(0.0); EDM.push_back(0.0); NBDM.push_back(0.0);	
////		cout << DensityAux << " " << E_nuc << " " << Etot << '\n';
//                while(DensityAux < Density){
//			//get the DM density from the baryon density in fm^-3
//                        DMDensity = gsl_spline_eval(spline_nChiofnB, NBARR[0], acc_nChiofnB) * 1e39*pow(convinvcmtoGeV,3.0);
//			//get the nuclear and DM EOSs
//                        P_nuc = gsl_spline_eval(spline_PofnB, DensityAux, acc_PofnB);
//			E_nuc = gsl_spline_eval(spline_EofP, P_nuc, acc_EofP);
//			P_kinDM = 1./(3*Pi*Pi)*Int_eyeP(0, kFermiDM, DMmass, 1000);	//using integral definition of P_kin
//			E_kinDM = 1./(Pi*Pi)*eye2(kFermiDM, DMmass);			//using analytic result for E_kin
//			Ptot = P_nuc + 0.5*eta*DMDensity*DMDensity * (1e3*pow(convGeVtoinvcm,3.0)*1e-39) + P_kinDM;	//MeV/fm^3
//			Etot = E_nuc + 0.5*eta*DMDensity*DMDensity * (1e3*pow(convGeVtoinvcm,3.0)*1e-39) + E_kinDM;	//MeV/fm^3
////			cout << DensityAux << " " << E_nuc << " " << Etot << '\n';
//			//push back the new EOS into vectors
//                        PDM.push_back(Ptot); EDM.push_back(Etot); NBDM.push_back(DensityAux);	
//			DensityAux *= 1.001;
//		}
//                DMDensity = gsl_spline_eval(spline_nChiofnB, NBARR[0], acc_nChiofnB) * 1e39*pow(convinvcmtoGeV,3.0);
//		P_kinDM = 1./(3*Pi*Pi)*Int_eyeP(0, kFermiDM, DMmass, 1000);	//using integral definition of P_kin
//		E_kinDM = 1./(Pi*Pi)*eye2(kFermiDM, DMmass);			//using analytic result for E_kin
//                P_nuc = gsl_spline_eval(spline_PofnB, Density, acc_PofnB); E_nuc = gsl_spline_eval(spline_EofP, P_nuc, acc_EofP);
//		Ptot = P_nuc + 0.5*eta*DMDensity*DMDensity * (1e3*pow(convGeVtoinvcm,3.0)*1e-39) + P_kinDM;	//MeV/fm^3
//		Etot = E_nuc + 0.5*eta*DMDensity*DMDensity * (1e3*pow(convGeVtoinvcm,3.0)*1e-39) + E_kinDM;	//MeV/fm^3
////		cout << Density << " " << E_nuc << " " << Etot << '\n';
//		//push back the new EOS into vectors
//                PDM.push_back(Ptot); EDM.push_back(Etot); NBDM.push_back(Density);	
//		//copy the new EOS into arrays for interpolation
//		double BDENS_DM[PDM.size()]; double ENDENS_DM[PDM.size()]; double PRESS_DM[PDM.size()];
//		for(int i = 0; i < PDM.size(); i++){
//			BDENS_DM[i] = NBDM[i]; ENDENS_DM[i] = EDM[i]; PRESS_DM[i] = PDM[i];
////			cout << i << " " << BDENS_DM[i] << " " << ENDENS_DM[i] << " " << PRESS_DM[i] << '\n';
//		}
//		//interpolate the new EOS in the same way as the nuclear EOS
//		spline_EofP_DM = gsl_spline_alloc(gsl_interp_linear, PDM.size());
//		gsl_spline_init (spline_EofP_DM, PRESS_DM, ENDENS_DM, PDM.size());
//		spline_PofnB_DM = gsl_spline_alloc(gsl_interp_linear, PDM.size());
//		gsl_spline_init (spline_PofnB_DM, BDENS_DM, PRESS_DM, PDM.size());
////		cout << "PCHDM: " << gsl_spline_eval(spline_PofnB_DM, BDENS_DM[100], acc_PofnB_DM) << " " << BDENS_DM[PDM.size() - 1] << '\n';
//		Density = BDENS_DM[PDM.size() - 1];
//		NBDM.clear(); PDM.clear(); EDM.clear();
//		//====================================================
//		//===== SOLVE THE TOV EQUATIONS WITH THE NEW EOS =======
//		//===================================================
//                use_DMEOS = 1;
//		check_which_EOS = 0;
//		tov_solver_solve(&results, Density, lowlim, use_DMEOS);
//
		if((results.r*1e-3 >=0.0) and (results.r*1e-3<=200)){
			myfile2 << results.m << " " << results.r*1e-3 << " " << Density << '\n';
//                	cout << "DM:  " << Density << " " << results.m << " " << results.r*1e-3 << '\n';
		}
		if (results.m < 0.0){
			myfile2 << 0 << " " << 0 << " " << Density << " " << pr4 << " " << pr5 << " " << pr6 << '\n';
			myfile2.close();
			return 0;
		}

		Density *= 1.05;
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
	gsl_spline_free (spline_nChiofnB);
	gsl_interp_accel_free (acc_nChiofnB);
	//5
	gsl_spline_free (spline_EofP_DM);
	gsl_interp_accel_free (acc_EofP_DM);
	//6
	gsl_spline_free (spline_PofnB_DM);
	gsl_interp_accel_free (acc_PofnB_DM);

	myfile2.close();
	return 0;
}
