vector< vector<double> > EOSDATA;
int readEOS(string EOSFILE, bool print_EOS){
	int rows = 0;										//counts rows 
	int columns = 0;									//counts columns
        //EOS file should hold 6 columns   [1]:Density [fm^-3]   [2]:Pressure [MeV fm^-3]   [3]:Energy Density [MeV fm^-3]
	//                                 [4]:SIGCHIN [cm^2]    [5]:SIGCHI2 [cm^2]         [6]:Mchi [GeV]
	double pr1 = 0.0; double pr2 = 0.0; double pr3 = 0.0; double pr4 = 0.0; double pr5 = 0.0; double pr6 = 0.0;
        // only the first three are relevant for the MR sequence, so initialize these to be < 0 for monotonicity
	double pr1p = -2.; double pr2p = -2.; double pr3p = -2.;

	stringstream file_name_eos;
	file_name_eos << EOSFILE;
	string	MAINTITLE = file_name_eos.str();
	ifstream dater;										//input stream called dater
	dater.open(MAINTITLE);									//open the input for reading

	if (!dater){
		cout << "Corresponding EoS File Not Found: " << " " << MAINTITLE <<  '\n';
		return 1;
	}
	//while the input stream is open, read three parameters per line, but first push back a zero
	vector<double> rowvector; rowvector.push_back(pr1); rowvector.push_back(pr2); rowvector.push_back(pr3);
	EOSDATA.push_back(rowvector);
	rowvector.clear();
	pr1p = pr1; pr2p = pr2; pr3p = pr3;
	double trial_shift = 0.0;

	while(!dater.eof()){
		dater >> pr1 >> pr2 >> pr3 ;
		//ensure monotonicity
		if(pr2> 1.*pr2p and pr3>1.*pr3p){
			if(print_EOS){
				cout << pr1 << " " << pr2 << " " << pr3 << '\n';
			}
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

	return rows;
}
