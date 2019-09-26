#include <string>
#include <algorithm>
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
string particle = "FERMIONS";
/*======================================================================================================
 * ====================== DECLARE VARIABLES THAT STORE RELEVANT PARAMETERS =============================
 * ====================================================================================================*/
	double MASSscan = 0;
	double RADIUSscan = 0;
	double SIGCHIN = 0; double SIGCHI2 = 0; double DMMASS = 0;
	double pr1 = 0; double pr2 = 0; double pr3 = 0;

	int startsig = atoi(argv[1]);
	int endsig = atoi(argv[2]);
	int sigchi2 = atoi(argv[3]);
	int startm = atoi(argv[4]);
	int endm = atoi(argv[5]);
	string star = argv[6];
	string model = argv[7];
	int massinc = atoi(argv[8]);

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

	stringstream datfs68;
	stringstream datfs95;
	datfs68 << "./ParameterSpace/ParamSpace_" << particle << "_" << startm << endm -1 << "_" << star << "_" << model << "_68.dat";
	datfs95 << "./ParameterSpace/ParamSpace_" << particle << "_" << startm << endm -1 << "_" << star << "_" << model << "_95.dat";
	string filetit68 = datfs68.str();
	string filetit95 = datfs95.str();
	ofstream myfile68;
	myfile68.open(filetit68);
	ofstream myfile95;
	myfile95.open(filetit95);

/*======================================================================================================
 * ======================    RUN THROUGH ALL CALCULATED PARAMETERS AND COMPARE 	  ======================
 *  ===================== EACH MASS AND RADIUS TO THE ONES IN CONFIDENCE INTERVAL ======================
 * ====================================================================================================*/
cout << "Increment: " << massinc << '\n';
for(int k = startm; k<endm; k+=massinc){
//if(k==startm){
//	massinc = 100;
//}
//else{massinc = 100;}
double SIGCHINPREV68 = 1e-45;
double SIGCHINPREV95 = 1e-45;
for(int j = endsig-1; j >= startsig; j--){
//cout << j << '\n';
double i = 0.9;
vector<double> accepted68;
vector<double> accepted95;
while(i < 9.8){
	i += 0.1;
	vector< vector<double> > MRDATA;
	stringstream liltitle;
	if (floor(i) == i){
		liltitle << "./"<< particle <<"/MRFiles/"<< k <<"GeV/MR_"<< int(i) << j << sigchi2 << k << "p.dat";
	}
	else{
		liltitle << "./"<< particle <<"/MRFiles/"<< k <<"GeV/MR_"<< i << j << sigchi2 << k << "p.dat";
	}
	
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
		SIGCHIN = i*pow(10., -1*j);
		SIGCHI2 = pow(10., -1*sigchi2);
		DMMASS = k;
		for (int bsc = 0; bsc < rows68; bsc++){ 
			double relerrMass = fabs( MRDATA[sc][0] - MRBENCH68[bsc][0] )/max( MRDATA[sc][0], MRBENCH68[bsc][0] );
			double relerrRadius = fabs( MRDATA[sc][1] - MRBENCH68[bsc][1] )/max( MRDATA[sc][1], MRBENCH68[bsc][1] );
			/*======================================================================================================
			 * ======== COUNT THE NUMBER OF TIMES THERE IS A MATCH WITH A PAIR OF M-R IN CONFIDENCE INTERVAL =======
			 * ====================================================================================================*/
			double error_tolerance = 1e-2;
			if (relerrMass <= error_tolerance and relerrRadius <= error_tolerance){
				//cout << "M: " << MRDATA[sc][0] << ", R: " << MRDATA[sc][0] << '\n';
				count68 ++;
				accepted68.push_back(SIGCHIN);
			}
		}
		for (int bsc = 0; bsc < rows95; bsc++){ 
			double relerrMass = fabs( MRDATA[sc][0] - MRBENCH95[bsc][0] )/max( MRDATA[sc][0], MRBENCH95[bsc][0] );
			double relerrRadius = fabs( MRDATA[sc][1] - MRBENCH95[bsc][1] )/max( MRDATA[sc][1], MRBENCH95[bsc][1] );
			/*======================================================================================================
			 * ======== COUNT THE NUMBER OF TIMES THERE IS A MATCH WITH A PAIR OF M-R IN CONFIDENCE INTERVAL =======
			 * ====================================================================================================*/
			double error_tolerance = 1e-2;
			if (relerrMass <= error_tolerance and relerrRadius <= error_tolerance){
				//cout << "M: " << MRDATA[sc][0] << ", R: " << MRDATA[sc][0] << '\n';
				count95 ++;
				accepted95.push_back(SIGCHIN);
			}
		}

 	}
	if (count68 != 0){
		SIGCHINPREV68 = *max_element(accepted68.begin(), accepted68.end());
//		cout << "68: " << SIGCHINPREV68 << " " << DMMASS << '\n';
	}
	if (count95 != 0){
		SIGCHINPREV95 = *max_element(accepted95.begin(), accepted95.end());
//		cout << "95: " << SIGCHINPREV95 << " " << DMMASS << '\n';
	}
	/*======================================================================================================
	 * ==== IF THE COUNT IS ZERO, THEN THESE PARAMETERS ARE EXCLUDED (M-R OUTSIDE OF ACCEPTABLE RANGE) =====
	 * ====================================================================================================*/
//	if (count68 == 0 and SIGCHIN < SIGCHINPREV68){
//		cout << "::Excluded 68%:: " << "sigchin: " << SIGCHIN << ", sigchi2: " << SIGCHI2 << ", mass: " << DMMASS << '\r' << flush;
//		SIGCHINPREV68 = SIGCHIN;
//	}
//	if (count95 == 0 and SIGCHIN < SIGCHINPREV95){
//		cout << "::Excluded 95%:: " << "sigchin: " << SIGCHIN << ", sigchi2: " << SIGCHI2 << ", mass: " << DMMASS << '\r' << flush;
//		SIGCHINPREV95 = SIGCHIN;
//	}

//	cout << accepted68[0] << '\n';
	MRDATA.clear();
	accepted68.clear();
	accepted95.clear();
}
}

		myfile68 << SIGCHINPREV68 << " " << SIGCHI2 << " " << DMMASS << '\n';
		myfile95 << SIGCHINPREV95 << " " << SIGCHI2 << " " << DMMASS << '\n';
		cout << "68: " << SIGCHINPREV68 << " " << SIGCHI2 << " " << DMMASS << '\n';
		cout << "95: " << SIGCHINPREV95 << " " << SIGCHI2 << " " << DMMASS << '\n';
}
	myfile68.close();
	myfile95.close();

	return 0;
}

