#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>
#include <fstream>
#include "CONSTANTSandCONVERSIONS.hh"
#include <time.h>
/*======================================================================================================
 * ====================== MAKE GLOBAL VECTOR TO STORE DATA BEING COMPARED ==============================
 * ====================================================================================================*/

int main(int argc, char* argv[]){
string particle = "BOSONS";
/*======================================================================================================
 * ====================== DECLARE VARIABLES THAT STORE RELEVANT PARAMETERS =============================
 * ====================================================================================================*/
	double MASSscan = 0;
	double RADIUSscan = 0;
	double SIGCHIN = 0; double SIGCHI2 = 0; double DMMASS = 0;
	double pr1 = 0; double pr2 = 0; double pr3 = 0;

	int coeff = atoi(argv[1]);
	int sigchin = atoi(argv[2]);
	int sigchi2 = atoi(argv[3]);
	int mass = atoi(argv[4]);
	string star = argv[5];
	string model = argv[6];

/*======================================================================================================
 * ====================== STORE CONFIDENCE INTERVALS (F.OZEL et. al.) ==================================
 * ====================================================================================================*/
	stringstream liltitleMRprob68;
	stringstream liltitleMRprob95;

	liltitleMRprob68 << "../MRprob/MRprob_" << star << "_68.dat";	
	liltitleMRprob95 << "../MRprob/MRprob_" << star << "_95.dat";	
	string BIGtitleMRprob68 = liltitleMRprob68.str();
	string BIGtitleMRprob95 = liltitleMRprob95.str();
	ifstream dater68;
	ifstream dater95;
	dater68.open(BIGtitleMRprob68);
	dater95.open(BIGtitleMRprob95);
	if (!dater68 or !dater95){
		cout << "Corresponding MR Confidence Interval Files Not Found... EXITING" << '\n';
		return 1;

	}

	vector< vector<double> > MRBENCH68;
	vector< vector<double> > MRBENCH95;

	double pr168 = 0; double pr268 = 0; double pr368 = 0;
	double pr195 = 0; double pr295 = 0; double pr395 = 0;
	int col68 = 0; int col95 = 0;

	while(!dater68.eof()){
		dater68 >> pr168 >> pr268 >> pr368;
		vector<double> col;
		col.push_back(pr168);
		col.push_back(pr268);
		MRBENCH68.push_back(col);
		col68 = col.size();
	}

	while(!dater95.eof()){
		dater95 >> pr195 >> pr295 >> pr395;
		vector<double> col;
		col.push_back(pr195);
		col.push_back(pr295);
		MRBENCH95.push_back(col);
		col95 = col.size();
	}
	int rows68 = MRBENCH68.size();
	int rows95 = MRBENCH95.size();

/*======================================================================================================
 * ======================    RUN THROUGH ALL CALCULATED PARAMETERS AND COMPARE 	  ======================
 *  ===================== EACH MASS AND RADIUS TO THE ONES IN CONFIDENCE INTERVAL ======================
 * ====================================================================================================*/
	vector< vector<double> > MRDATA;
	stringstream liltitle;
	liltitle << "./"<< particle <<"/MRFiles/"<< mass <<"GeV/MR_"<< coeff << sigchin << sigchi2 << mass << "p.dat";
	string BIGTITLE = liltitle.str();
//	cout << BIGTITLE << '\n';
	ifstream dater;
	dater.open(BIGTITLE);
	if (!dater){
		MRDATA.clear();
		cout << "MR FILE NOT FOUND... EXITING..." << '\n';
		cout << BIGTITLE << '\n';
		return 1;
	}
	int columns = 0;
	while(!dater.eof()){
		dater >> pr1 >> pr2 >> pr3;
		vector<double> row;
		row.push_back(pr1);
		row.push_back(pr2);

		MRDATA.push_back(row);
		columns = row.size();
	}

	int rows = MRDATA.size();
	dater.close();

	int count68 = 0;
	int count95 = 0;
	for(int sc = 0; sc < rows; sc++){
		MASSscan = MRDATA[sc][0];
		RADIUSscan = MRDATA[sc][1];
		SIGCHIN = coeff*pow(10., -1*sigchin);
		SIGCHI2 = pow(10., -1*sigchi2);
		DMMASS = mass;
		for (int bsc = 0; bsc < rows68; bsc++){ 
			double relerrMass = fabs( MRDATA[sc][0] - MRBENCH68[bsc][0] )/max( MRDATA[sc][0], MRBENCH68[bsc][0] );
			double relerrRadius = fabs( MRDATA[sc][1] - MRBENCH68[bsc][1] )/max( MRDATA[sc][1], MRBENCH68[bsc][1] );
			/*======================================================================================================
			 * ======== COUNT THE NUMBER OF TIMES THERE IS A MATCH WITH A PAIR OF M-R IN CONFIDENCE INTERVAL =======
			 * ====================================================================================================*/
			double error_toleranceM = 5e-3;
			double error_toleranceR = 5e-3;
			if (relerrMass <= error_toleranceM and relerrRadius <= error_toleranceR){
				cout << "M: " << MRDATA[sc][0] << ", R: " << MRDATA[sc][1] << '\n';
				cout << "Mcomp: " << MRBENCH68[bsc][0] << ", Rcomp: " << MRBENCH68[bsc][1] << '\n';
				count68 ++;
			}
			else{continue;}
		}
		for (int bsc = 0; bsc < rows95; bsc++){ 
			double relerrMass = fabs( MRDATA[sc][0] - MRBENCH95[bsc][0] )/max( MRDATA[sc][0], MRBENCH95[bsc][0] );
			double relerrRadius = fabs( MRDATA[sc][1] - MRBENCH95[bsc][1] )/max( MRDATA[sc][1], MRBENCH95[bsc][1] );
			/*======================================================================================================
			 * ======== COUNT THE NUMBER OF TIMES THERE IS A MATCH WITH A PAIR OF M-R IN CONFIDENCE INTERVAL =======
			 * ====================================================================================================*/
			double error_tolerance = 5e-3;
			if (relerrMass <= error_tolerance and relerrRadius <= error_tolerance){
				//cout << "M: " << MRDATA[sc][0] << ", R: " << MRDATA[sc][0] << '\n';
				count95 ++;
			}
		}

 	}
	/*======================================================================================================
	 * ==== IF THE COUNT IS ZERO, THEN THESE PARAMETERS ARE EXCLUDED (M-R OUTSIDE OF ACCEPTABLE RANGE) =====
	 * ====================================================================================================*/
	if (count68 == 0){
		cout << "::Excluded 68%:: " << "sigchin: " << SIGCHIN << ", sigchi2: " << SIGCHI2 << ", mass: " << DMMASS << '\n';
	}
	if (count95 == 0){
		cout << "::Excluded 95%:: " << "sigchin: " << SIGCHIN << ", sigchi2: " << SIGCHI2 << ", mass: " << DMMASS << '\n';
	}
	MRDATA.clear();

	return 0;
}


