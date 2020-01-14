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

int main(int argc, char* argv[]){
	string EOSfile = argv[1];
	bool PRINTEOS = 0;
	int rows = readEOS(EOSfile, PRINTEOS);
	cout << "NOROWS: " << rows << '\n';
	//copy the EOS into arrays, GSL doesn't like vectors
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
	cout << "NBTRIAL: " << nBtrial << " " << BDENS[0] <<  " " << 1.5 << " " << DensityMax << '\n';
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

	EOSDATA.clear();
	return 0;
}
