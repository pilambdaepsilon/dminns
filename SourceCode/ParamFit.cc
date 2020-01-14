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
	vector< vector<double> > DATA68;
	vector< vector<double> > DATA95;
	double pr1 = 0; double pr2=0; double pr3=0; double pr1p=0; double pr2p=0; double pr3p=0; 

	int ParamStart = atoi(argv[1]);
	int ParamEnd = atoi(argv[2]);
	string star = argv[3];
	string model = argv[4];
	string particle = "BOSONS";

	stringstream liltitle1;
	stringstream liltitle2;

	liltitle1 <<"./ParameterSpace" << "/ParamSpace_" << particle << "_" << ParamStart << ParamEnd << "_" << star << "_" << model << "_" << "68.dat";
	liltitle2 <<"./ParameterSpace" << "/ParamSpace_" << particle << "_" << ParamStart << ParamEnd << "_" << star << "_" << model << "_" << "95.dat";

	string	BIGTITLE1 = liltitle1.str();
	string	BIGTITLE2 = liltitle2.str();

	ifstream dater1;
	ifstream dater2;

	dater1.open(BIGTITLE1);
	dater2.open(BIGTITLE2);

	if (!dater1){
		cout << "1 Sigma File Not Found: " << " " << BIGTITLE1 <<  '\n';
		return 1;
	}
	if (!dater2){
		cout << "2 Sigma EoS File Not Found: " << " " << BIGTITLE2 <<  '\n';
		return 1;
	}

	int columns1 = 0;
	int columns2 = 0;

	while(!dater1.eof()){
		dater1 >> pr1 >> pr2 >> pr3;
		if(pr3>pr3p){
			vector<double> row;
			row.push_back(pr1);
			row.push_back(pr3);
			DATA68.push_back(row);
			columns1 = row.size();
			pr1p = pr1;
			pr3p = pr3;
			cout << pr1 << " " << pr3 << '\n';
		}
		else{continue;}
	}

	pr1 = 0; pr2 = 0; pr3 = 0; pr1p = 0; pr2p = 0; pr3p = 0;
	while(!dater2.eof()){
		dater2 >> pr1 >> pr2 >> pr3;
		if(pr3>pr3p){
			vector<double> row;
			row.push_back(pr1);
			row.push_back(pr3);
			DATA95.push_back(row);
			columns1 = row.size();
			pr1p = pr1;
			pr3p = pr3;
		cout << pr1 << " " << pr3 << '\n';
		}
		else{continue;}
	}

	dater1.close();
	dater2.close();

	int rows68 = DATA68.size();
	int rows95 = DATA95.size();
	int n68 = 0;
	int n95 = 0;
	n68 = rows68;
	n95 = rows95;
	cout << rows68 << " " << rows95 << '\n';
	double MASS68[rows68] = {0.0};
	double SIGCHIN68[rows68] = {0.0};

	double MASS95[rows95] = {0.0};
	double SIGCHIN95[rows95] = {0.0};

	for (int i = 1; i < rows68; i++){
		MASS68[i] = DATA68[i][1];
		SIGCHIN68[i] = DATA68[i][0];
	}
	for (int i = 1; i < rows95; i++){
		MASS95[i] = DATA95[i][1];
		SIGCHIN95[i] = DATA95[i][0];
	}

	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	gsl_spline *spline = gsl_spline_alloc(gsl_interp_steffen, n68);
	gsl_spline_init (spline, MASS68, SIGCHIN68, n68);

	double mass = 1;
	stringstream FitFileName1;
	stringstream FitFileName2;

	FitFileName1 << "ParameterSpace/ParamFit_" << particle << "_" << ParamStart << ParamEnd << "_" << star << "_" << model << "_68.dat";
	FitFileName2 << "ParameterSpace/ParamFit_" << particle << "_" << ParamStart << ParamEnd << "_" << star << "_" << model << "_95.dat";

	string FitFileNameStr1 = FitFileName1.str();
	string FitFileNameStr2 = FitFileName2.str();

	ofstream FitFile1;
	ofstream FitFile2;
	FitFile1.open(FitFileNameStr1);
	FitFile2.open(FitFileNameStr2);

//	cout << "WORKS: " << DATA68[rows68-1][1] << '\n';
	while(mass <= DATA68[rows68-1][1]){
		double XsecOfMass = gsl_spline_eval(spline, mass, acc);
//		cout << XsecOfMass << " " << mass << '\n';
		FitFile1 << XsecOfMass << " " << mass << '\n';
		mass += 1;
	}

	mass = 1;
	gsl_spline_free (spline);
	gsl_interp_accel_free (acc);

	gsl_interp_accel *acc2 = gsl_interp_accel_alloc();
	gsl_spline *spline2 = gsl_spline_alloc(gsl_interp_steffen, n95);
	gsl_spline_init (spline2, MASS95, SIGCHIN95, n95);
	
	while(mass <= DATA95[rows95-1][1]){
		double XsecOfMass = gsl_spline_eval(spline2, mass, acc2);
//		cout << XsecOfMass << " " << mass << '\n';
		FitFile2 << XsecOfMass << " " << mass << '\n';
		mass += 1;
	}
		

	gsl_spline_free (spline2);
	gsl_interp_accel_free (acc2);
	
	FitFile1.close();
	FitFile2.close();

	return 0;
}
