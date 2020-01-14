#include "EoSrootfinder.hh"

int main(int argc, char* argv[]){
	double title1 = atof(argv[1]);
	double title2 = atof(argv[2]);
	double title3 = atof(argv[3]);
	double title4 = atof(argv[4]);
	string STAR = argv[5];
	string EoSMODEL = argv[6];
        string particle = argv[7];

/*====================== TAKE IN CONFIGURATION FILE TO GET STAR AND MODEL DETAILS ==========================*/
	string starname; double dNS_Earth; double dNS_GC; double NSage; double DMAmbientDensity;

	string storevalues[5];
	int incr = 0;
	stringstream starconf;
	starconf << "./ConfigurationFiles/" << STAR << ".txt";
	string starconfstr = starconf.str();
	ifstream is_file;
	is_file.open(starconfstr);
	if (!is_file){
		cout << "Stellar Configuration File Not Found" << '\n';
	}
	string line;
	
	while(getline(is_file, line) ){
		istringstream is_line(line);
		string key;
		if(getline(is_line, key, '=') ){
			string value;
			if(getline(is_line, value) ){
				storevalues[incr] = value;
				incr++;
			}
		}
	}

	starname = storevalues[0]; dNS_Earth = stof(storevalues[1]); dNS_GC = stof(storevalues[2]); NSage = stof(storevalues[3]); 
	DMAmbientDensity = stof(storevalues[4]);
	cout << "Star: " << starname << '\n' <<
	       	"Distance to Earth: " << dNS_Earth << '\n' <<
		"Distance to Galactic Center: " << dNS_GC << '\n' << 
		"Age: " << NSage << '\n' << 
		"Ambient DM Density [NFW]: " << DMAmbientDensity << '\n';


/*====================== SET UP OUTPUT DIRECTORY AND PARTICLE TYPE ==========================*/
	stringstream ss;
	ss << title1 << title2 << title3 << title4;
	string fno = ss.str();

	cout << '\n' << "fno: " << fno << '\n';
	
	int calls = 1000;

	stringstream liltitle1b;
	liltitle1b << "./" << particle << "/EoSFiles/" << title4 << "GeV/EoS_SigOmDM_" << fno << "p.dat";
	string BIGTITLEb = liltitle1b.str();
	ofstream DATERclean;
	DATERclean.open(BIGTITLEb);

/*====================== SET UP INITIAL CONDITIONS ==========================*/
	double NMass = 0.938;

/*====================== READ AND FIT EXTERNAL EOS MODEL WITH CUBIC SPLINE ==========================*/
	vector< vector<double> > EOSDATA;

	stringstream liltitleEoS;
	liltitleEoS <<"./" << "EoSFORTIN/EOS_"<<EoSMODEL<<".dat";
	string	BIGTITLE1 = liltitleEoS.str();
	ifstream EoSMODELdater;
	EoSMODELdater.open(BIGTITLE1);

	if (!EoSMODELdater){
		cout << '\n' << "EoS Model File not Found: " << " " << BIGTITLE1 << "... EXITING" << '\n';
		return 1;
	}

	int columnsEoS = 0;
	double pr1EoS =0; double pr2EoS =0; double pr3EoS = 0; double pr2pEoS = 0; double pr3pEoS = 0;
	while(!EoSMODELdater.eof()){
		EoSMODELdater >> pr1EoS >> pr2EoS >> pr3EoS;
		if(pr2EoS> pr2pEoS and pr3EoS>pr3pEoS){
			vector<double> row;
			row.push_back(pr1EoS);
			row.push_back(pr2EoS);
			row.push_back(pr3EoS);
			EOSDATA.push_back(row);
			columnsEoS = row.size();
			pr2pEoS = pr2EoS;
			pr3pEoS = pr3EoS;
//			cout << pr1EoS << " " << pr2EoS << " " << pr3EoS << '\n';
		}
		else{continue;}
	}

	EoSMODELdater.close();

	int rowsEoS = EOSDATA.size() - 1;
	int nEoS = rowsEoS - 1;
	BaryonDensity = EOSDATA[0][0];
	double BaryonDensMax = EOSDATA[nEoS-1][0];
	double PRESS[rowsEoS];
	double BARYONDENSITY[rowsEoS];
	double ENDENS[rowsEoS];

	for (int i = 1; i < rowsEoS; i++){
		PRESS[i] = EOSDATA[i][2];
		ENDENS[i] = EOSDATA[i][1];
		BARYONDENSITY[i] = EOSDATA[i][0];
	}

	gsl_interp_accel *accE = gsl_interp_accel_alloc();
	gsl_spline *splineE = gsl_spline_alloc(gsl_interp_akima, nEoS);
	gsl_spline_init (splineE, BARYONDENSITY, ENDENS, nEoS);
	
	gsl_interp_accel *accP = gsl_interp_accel_alloc();
	gsl_spline *splineP = gsl_spline_alloc(gsl_interp_akima, nEoS);
	gsl_spline_init (splineP, BARYONDENSITY, PRESS, nEoS);

/*================= READ IN DENSITY INFORMATION FROM STELLAR CAPTURE CODE ========*/
	stringstream liltitle2;
	liltitle2 << "./" << particle << "/DmNoFiles/" << title4 << "GeV/nchi" << fno << "pSKINNY.dat";
	string BIGTITLE2 = liltitle2.str();
	ifstream dater;
	dater.open(BIGTITLE2);
	if (!dater){
		cout << "Corresponding Number Density File Not Found.. EXITING" << '\n';
		return 1;

	}
	
	double SIGCHIN = 0.0;
	double SIGCHI2 = 0.;
	double DMmass =0.;
	vector< vector<double> > NCHIDATA;
	int columns = 0;
	double pr1 =0; double pr2 =0; double pr3 =0; double pr4 =0; double pr5 =0; double pr6 =0; double pr7 =0;

	while(!dater.eof()){
		dater >> pr1 >> pr2 >> pr3 >> pr4 >> pr5 >> pr6 >> pr7;
		vector<double> col;
		col.push_back(pr1);
		col.push_back(pr3);
		col.push_back(pr4);
		SIGCHIN = pr5*1e26; 
		SIGCHI2 = pr6*1e26;
		DMmass = pr7*1e3*MeVtoinvFM;
		NCHIDATA.push_back(col);
		columns = col.size();
	}

	cout << "EoS Model: " << EoSMODEL << '\n';
	int rows = NCHIDATA.size();
	dater.close();

	double maxtimeforNS = NCHIDATA[rows-1][0];


	double DMDensity = 0;//PARAMSET[DMTimeSlice][2];			//fm^-3
	double dummyscannB = 0;
	double dummyscanT =0;
	
/*================= SOLVE THE EOS + DM EFFECT AT ALL BARYON DENSITIES USING FITTED EOS ========*/
	while (BaryonDensity < BaryonDensMax){
		AlphaAsymmetry = 0.;
		for(int scan = 0; scan < rows; scan++){
			dummyscannB = NCHIDATA[scan][2];
			double relerrDENS = fabs(dummyscannB - BaryonDensity)/max(dummyscannB, BaryonDensity);
			if(relerrDENS <= 5e-3){
				DMDensity = NCHIDATA[scan][1];
			}
		}
		
		double GPI = 0; double GCHI = 0; double GCHIb = 0; double LAMBDACHIb = 0;


		GCHI = 1;//Gchi(SIGCHI2, DMmass, MPI);
		const double MPI = Mpib(SIGCHIN, NMass, DMmass);
		GPI = Gpi(SIGCHIN, NMass, DMmass, MPI, GCHI);
		GCHIb = Gchib(SIGCHIN, NMass, DMmass, MPI);
		LAMBDACHIb = Lambdachi(SIGCHI2, DMmass);



		
/*============ USE FIT TO FIND UNMODIFIED EOS ==============================================================================*/
		double EnergyDensity = gsl_spline_eval(splineE, BaryonDensity, accE);
		double Pressure = gsl_spline_eval(splineP, BaryonDensity, accP);

/*============ DARK MATTER MODIFICATION OF EQUATION OF STATE ==============================================================================*/
		double EnergyDensityDM = 0;
		double PressureDM = 0;

		double kFermiDM = pow((3*pi*pi *DMDensity), (1./3.));
		double PDM = 0.0;
		double KineticEnergyDM = 0.0;

		if(particle == "FERMIONS"){
			PDM = 1./(3*Pi*Pi)*Int_eyeP(0, kFermiDM, DMmass, calls);
			KineticEnergyDM = 1./(Pi*Pi)*eye2(kFermiDM, DMmass);

			EnergyDensityDM = 0.5 * pow((GPI/MPI), 2.) * pow(BaryonDensity, 2)
			+ 0.5 * pow((GCHI/MPI), 2.) * pow(DMDensity, 2)
			+ 0.5 * (GPI/MPI)*(GCHI/MPI)* BaryonDensity*DMDensity
			+ KineticEnergyDM; 

			PressureDM = 0.5 * pow((GPI/MPI), 2.) * pow(BaryonDensity, 2)
			+ 0.5 * pow((GCHI/MPI), 2.) * pow(DMDensity, 2)
			+ 0.5 * (GPI/MPI)*(GCHI/MPI)* BaryonDensity*DMDensity
			+ PDM;
		//	cout << "Pressure: " << PressureDM << "   EnergyDensity: " << EnergyDensityDM << '\n';
		//	cout << "PDM: " << PDM << "   EkDM: " << KineticEnergyDM*2./3. << '\n';
		}

		else if (particle == "BOSONS"){
			double DMterm1a = 0; double DMterm2a = 0; double DMterm3a = 0; double DMterm4a = 0;
			double DMterm1b = 0; double DMterm2b = 0; double DMterm3b = 0; double DMterm4b = 0;

//			DMDensity = 0.0;
//			double alpha = GCHIb*GCHIb;
//			double beta = 2*DMmass*DMmass;
			//double gamma = LAMBDACHIb*DMDensity/(12.*GCHIb) + (GCHIb*DMDensity + GCHIb*BaryonDensity)*pow(DMmass/MPI,2.);
			double alpha = LAMBDACHIb/6.;
			double beta = DMmass*DMmass;
			double gamma = 2*GCHIb*GCHIb*DMDensity*(BaryonDensity + DMDensity)/pow(MPI,2.);

			double factor1 = 9*pow(alpha,2.)*gamma + sqrt(3.)*sqrt(27*pow(alpha,4.)*pow(gamma,2.) + 4*pow(alpha,3.)*pow(beta,3.));
			double Pi0 = GCHIb*(BaryonDensity + DMDensity)/pow(MPI,2.);
			double Chi0 = sqrt( (sqrt(4*alpha*gamma + beta*beta) - beta)/(2*alpha) );

			DMterm1a = 0.5*MPI*MPI*GCHIb*GCHIb*pow(BaryonDensity,2.)/pow((MPI*MPI + 0.5*GCHIb*GCHIb*DMDensity/DMmass), 2.);
			DMterm1b = 0.5*MPI*MPI*Pi0*Pi0;

			DMterm2a = 0.25*DMmass*DMDensity;
			DMterm2b = 0.5*DMmass*DMmass*Chi0*Chi0;

			DMterm3a = 0.5*pow(GCHIb, 4.)*pow(BaryonDensity,2.)*DMDensity/(DMmass*pow((MPI*MPI + 0.5*GCHIb*GCHIb*DMDensity/DMmass), 2.));
			DMterm3b = GCHIb*GCHIb*Pi0*Pi0*Chi0*Chi0;

			DMterm4a = 0.25*LAMBDACHIb/24. * pow((DMDensity/DMmass), 2.);
			DMterm4b = LAMBDACHIb/24.*pow(Chi0,4.);

			EnergyDensityDM = DMterm1b - DMterm2b - DMterm3b + DMterm4b + DMmass*DMDensity;
			PressureDM = DMterm1b + DMterm2b + DMterm3b - DMterm4b;

			cout << "====================================================" << '\n';
			cout << "nChi: " << DMDensity << " ... nB: " << BaryonDensity << '\n';
			cout << "Chi0: " << Chi0 << " ... Pi0: " << Pi0 << '\n';
			cout << "alpha: " << alpha << " ... beta: " << beta << " ... gamma: " << gamma << '\n';
			cout << "zero-point: " << DMmass*DMDensity << " ... beta: " << beta << " ... gamma: " << gamma << '\n';
			cout << "Lambda: " << LAMBDACHIb << " ... Gchi: " << GCHIb << " ... Mpi: " << MPI << '\n';
			cout << "1[prev/new]: " << DMterm1a << " " << DMterm1b << " " << '\n';
			cout << "2: " << DMterm2a << " " << DMterm2b << " " << '\n';
			cout << "3: " << DMterm3a << " " << DMterm3b << " " << '\n';
			cout << "4: " << DMterm4a << " " << DMterm4b << " " << '\n';
			cout << "E: " << EnergyDensityDM << " ... P: " << PressureDM << '\n';
		}

/*================= ADD THEEM UP ========*/
		EnergyDensity += EnergyDensityDM/MeVtoinvFM;
		Pressure += PressureDM/MeVtoinvFM;

/*============================== READ OUT AND OUTPUT TO CHECK =============================================*/
		if (isfinite(EnergyDensity) and isfinite(Pressure) and Pressure >= 0.0){
			DATERclean << BaryonDensity << " " << EnergyDensity << " " << Pressure << " " 
				<< SIGCHIN*1e-26 << " " << SIGCHI2*1e-26 << " " << DMmass*1e-3/MeVtoinvFM << '\n';

//			cout << "Energy Density: " << EnergyDensity << "   :: Baryon Density: " << BaryonDensity << "   :: Pressure: " << Pressure << '\n';
		}


		BaryonDensity += 0.001;
	}

	DATERclean.close();

	gsl_spline_free(splineE);
	gsl_interp_accel_free(accE);
	gsl_spline_free(splineP);
	gsl_interp_accel_free(accP);
	return 0;
	
}

