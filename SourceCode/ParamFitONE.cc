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

int main(int argc, char* argv[]){
	vector< vector<double> > DATA95;
	double pr1 = 0; double pr2=0; double pr3=0; double pr1p=0; double pr2p=0; double pr3p=0; 

	int ParamStart = atoi(argv[1]);
	int ParamEnd = atoi(argv[2]);
	string star = argv[3];
	string model = argv[4];
	string particle = "BOSONS";

	stringstream liltitle1;
	stringstream liltitle2;

	liltitle2 <<"./ParameterSpace" << "/ParamSpace_" << particle << "_" << ParamStart << ParamEnd << "_" << star << "_" << model << "_" << "95.dat";

	string	BIGTITLE2 = liltitle2.str();

	ifstream dater2;

	dater2.open(BIGTITLE2);

	if (!dater2){
		cout << "2 Sigma EoS File Not Found: " << " " << BIGTITLE2 <<  '\n';
		return 1;
	}

	int columns2 = 0;


	pr1 = 0; pr2 = 0; pr3 = 0; pr1p = 0; pr2p = 0; pr3p = 0;
	while(!dater2.eof()){
		dater2 >> pr1 >> pr2 >> pr3;
		if(pr1> pr1p and pr3>pr3p){
			vector<double> row;
			row.push_back(pr1);
			row.push_back(pr3);
			DATA95.push_back(row);
			columns2 = row.size();
			pr1p = pr1;
			pr3p = pr3;
		cout << pr1 << " " << pr3 << '\n';
		}
		else{continue;}
	}

	dater2.close();

	int rows95 = DATA95.size();
	int n95 = 0;
	n95 = rows95;
	double MASS95[rows95] = {0.0};
	double SIGCHIN95[rows95] = {0.0};

	for (int i = 1; i < rows95; i++){
		MASS95[i] = DATA95[i][1];
		SIGCHIN95[i] = DATA95[i][0];
	}


//	cout << "WORKS" << '\n';
	double mass = 1;
	stringstream FitFileName2;

	FitFileName2 << "ParameterSpace/ParamFit_" << particle << "_" << ParamStart << ParamEnd << "_" << star << "_" << model << "_95.dat";

	string FitFileNameStr2 = FitFileName2.str();

	ofstream FitFile2;
	FitFile2.open(FitFileNameStr2);

	gsl_interp_accel *acc2 = gsl_interp_accel_alloc();
	gsl_spline *spline2 = gsl_spline_alloc(gsl_interp_steffen, n95);
	gsl_spline_init (spline2, MASS95, SIGCHIN95, n95);
	while(mass <= float(ParamEnd)){
		double XsecOfMass = gsl_spline_eval(spline2, mass, acc2);
//		cout << XsecOfMass << " " << mass << '\n';
		FitFile2 << XsecOfMass << " " << mass << '\n';
		mass += 1;
	}
		

	gsl_spline_free (spline2);
	gsl_interp_accel_free (acc2);
	
	FitFile2.close();

	return 0;
}

