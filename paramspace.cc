#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>
#include <fstream>
#include "CONSTANTSandCONVERSIONS.hh"
#include <time.h>

vector< vector<double> > MRDATA;

int main(int argc, char* argv[]){

	double pr1 = 0;
	double pr2 = 0;
	double pr3 = 0;
	double pr4 = 0;
	double pr5 = 0;
	double pr6 = 0;
	double MASSscan = 0;
	double RADIUSscan = 0;
	double SIGCHINscan = 0; double SIGCHI2scan = 0; double DMMASSscan = 0;


	int endsig = atoi(argv[3]);
	int endm = atoi(argv[6]);
	int title3 = atoi(argv[4]);
	int title4 = atoi(argv[5]);
	stringstream datfs;
	datfs << "params" << title4 << endm-1 << ".dat";
	string filetit = datfs.str();
	ofstream myfile;
	myfile.open(filetit);
for(int k = title4; k<endm; k++){
	int title2 = atoi(argv[2]);
for(int j = title2; j < endsig; j++){
	int title1 = atoi(argv[1]);
for(int i = 0; i < 9; i++){
	stringstream ss;
	ss << title1 << title2 << title3 << title4;
	string fno = ss.str();

	stringstream liltitle;
	liltitle << "./MRFiles/"<<title4<<"GeV/MR_"<<fno<<"p.dat";
	string BIGTITLE = liltitle.str();
	ifstream dater;
	dater.open(BIGTITLE);
	if (!dater){
		MRDATA.clear();
		cout << "FILE NOT FOUND... EXITING..." << '\n';
		return 1;
	}
	int columns = 0;
	while(!dater.eof()){
		dater >> pr1 >> pr2 >> pr3 >> pr4 >> pr5 >> pr6;
		vector<double> row;
		row.push_back(pr1);
		row.push_back(pr2);
		row.push_back(pr4);
		row.push_back(pr5);
		row.push_back(pr6);

		MRDATA.push_back(row);
		columns = row.size();
	}

	int rows = MRDATA.size();
	dater.close();

	cout << "rows: " << rows << " columns: " << columns << '\n';

	cout << '\n' << "DATA COPYING DONE..." << '\n';


	int count = 0;
	for(int i = 0; i < rows; i++){
		MASSscan = MRDATA[i][0];
		RADIUSscan = MRDATA[i][1];
		SIGCHINscan = MRDATA[i][2];
		SIGCHI2scan = MRDATA[i][3];
		DMMASSscan = MRDATA[i][4];
		double relerrmass = fabs( MASSscan - 2.01 )/max(MASSscan, 2.01);
		if (RADIUSscan >= 11. and RADIUSscan <=15. and relerrmass <1e-1){
			cout << "M: " << MASSscan << ", R: " << RADIUSscan << '\n';
			count ++;
		}
 	}
	if (count == 0){
		cout << "::Excluded:: " << "sigchin: " << SIGCHINscan << ", sigchi2: " << SIGCHI2scan << ", mass: " << DMMASSscan << '\n';
		myfile << SIGCHINscan << " " << SIGCHI2scan << " " << DMMASSscan << '\n';
	}
	MRDATA.clear();
	title1 ++;
}
title2++;
}
title4++;
}
	myfile.close();

	return 0;
}
