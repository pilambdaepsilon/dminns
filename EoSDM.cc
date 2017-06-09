#include "EoSrootfinder.hh"

int main(int argc, char* argv[]){
	int title1 = atoi(argv[1]);
	int title2 = atoi(argv[2]);
	int title3 = atoi(argv[3]);
	int title4 = atoi(argv[4]);

	stringstream ss;
	ss << title1 << title2 << title3 << title4;
	string fno = ss.str();
	
//	cout << '\n' << '\n' << "WHAT SHALL WE CALL THIS FILE?: "; cin >> fno;
//	cout << '\n' << "WHAT IS THE FILENAME OF THE NUMBER DENSITY FILE? - Naming convention is [55241] for SIGCHIN~-55, SIGCHI2~-24, MDM~ 1: "; cin >> fno2;
	int calls = 1000;
	string particle = "BOSONS";

	stringstream liltitle1a;
	stringstream liltitle1b;
	liltitle1a << "./" << particle << "/EoSFiles/" << title4 << "GeV/SigOm" << fno << "pDM.dat";
	liltitle1b << "./" << particle << "/EoSFiles/" << title4 << "GeV/EoS_SigOmDM_" << fno << "p.dat";
	string BIGTITLEa = liltitle1a.str();
	string BIGTITLEb = liltitle1b.str();
	ofstream DATER;
	ofstream DATERclean;
	//DATER.open(BIGTITLEa);
	DATERclean.open(BIGTITLEb);

	//SOLVE THE PROBLEM AT SATURATION FOR CONSTANTS OF THE THEORY *******************************************!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*************!*!*!*!*!**!
	SolveSaturation();

	double ChemicalPotentialNeutron = 1.; double ChemicalPotentialProton = 1.; double ChemicalPotentialElectron = 2.; double ChemicalPotentialMuon = 2.;
	ScalarField = 0.1; VectorField = 1.; IsoVectorField = 0.5; BaryonDensity = 0.001; double Asigma = SaturationVariables[3]; double Aomega = SaturationVariables[4];
	double Arho = SaturationVariables[5]; B = SaturationVariables[6]; C = SaturationVariables[7]; HadronCouplingOmega = SaturationVariables[8]; double kFermiNeutron = 1.;
	double kFermiProton = 1.; double kFermiElectron = 1.; double kFermiMuon = 1.; double kFermiLambda = 1.; double kFermiXiMinus = 1.; double kFermiSigmaMinus = 1.;
	double kFermiSigmaZero = 1.; double kFermiSigmaPlus = 1.; double kFermiXiZero = 1.;

	HadronCouplingRho = HadronCouplingSigma;
/*================= SOLVE FOR FIELD AMPLITUDES WITH KNOWN CONSTANTS FROM SATURATION ========*/
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

	int rows = NCHIDATA.size();
	dater.close();

	cout << "rows: " << rows << " columns: " << columns << '\n';


	double maxtimeforNS = NCHIDATA[rows-1][0];


	double NSage = NCHIDATA[0][0];//PARAMSET[DMTimeSlice][0];			//years
	double DMDensity = 0;//PARAMSET[DMTimeSlice][2];			//fm^-3
	double dummyscannB = 0;
	double dummyscanT =0;
	cout << NSage << " years" << '\n' << '\n';
	
	while (BaryonDensity <= 1.6){
		AlphaAsymmetry = 0.;
		for(int scan = 0; scan < rows; scan++){
			dummyscannB = NCHIDATA[scan][2];
			double relerrDENS = fabs(dummyscannB - BaryonDensity)/max(dummyscannB, BaryonDensity);
			if(relerrDENS <= 5e-3){
				DMDensity = NCHIDATA[scan][1];
				cout << BaryonDensity << " " << dummyscannB << " " << DMDensity << '\n';
			}
		}
		
		const gsl_multiroot_fsolver_type *T2;		//DECLARE VARIABLE FOR TYPE OF SOLVER
		gsl_multiroot_fsolver *S2;		//DECLARE NAME FOR SOLVER WORKSPACE

		int status2;				//USED TO UPDATE ON STATUS (IS ACTUALLY BOOLEAN)
		size_t i2, iter2 =0;
	
		const size_t n2 = 15;			//DIMENSIONALITY OF PROBLEM
				//  Msig,  Mw,   Mr,  Mn,   Me,  gs, gw, gr,   b,    c,    I3,  Jn,  Jp,  Je,  ban, bap, Qe,  Qp, nB             Easym, Kcompress   Eo 

		double GPI = 0; double GCHI = 0; double GCHIb = 0; double LAMBDACHIb = 0;

		struct Rparams p2 = {ScalarMass, VectorMass, IsoVectorMass, NMass, ElectronMass, Asigma, Aomega, Arho, B, C, BaryonDensity, AsymmetryEnergy,
			Compressibility, BindingPerNucleon0, MstarPerMn0, AlphaAsymmetry, LambdaMass, HadronCouplingOmega, HadronCouplingSigma, HadronCouplingRho,
			MuonMass, XiMinusMass, SigmaMinusMass, SigmaZeroMass, SigmaPlusMass, XiZeroMass};

		gsl_multiroot_function func2 = {&GM2_f, n2, &p2};	//MAKE A FUNCTION FOR SOLVER OUT OF ROSENBROCK, OF DIMENSION n, WITH PARAMS p

		double y_initial[n2] = {ScalarField, VectorField, IsoVectorField, ChemicalPotentialNeutron, ChemicalPotentialElectron,  
			kFermiNeutron, kFermiProton, kFermiElectron, kFermiMuon, kFermiLambda, kFermiSigmaMinus, kFermiSigmaZero, kFermiSigmaPlus,
		       	kFermiXiMinus, kFermiXiZero};

		gsl_vector *y = gsl_vector_alloc(n2);	//MAKE A VECTOR FOR SOLVER TO STORE VARIABLES
	
		gsl_vector_set(y, 0, y_initial[0]);
		gsl_vector_set(y, 1, y_initial[1]);	//SET VARIABLES TO BE SOLVED FOR TO INITIAL VALUES
		gsl_vector_set(y, 2, y_initial[2]);
		gsl_vector_set(y, 3, y_initial[3]);
		gsl_vector_set(y, 4, y_initial[4]);
		gsl_vector_set(y, 5, y_initial[5]);
		gsl_vector_set(y, 6, y_initial[6]);
		gsl_vector_set(y, 7, y_initial[7]);
		gsl_vector_set(y, 8, y_initial[8]);
		gsl_vector_set(y, 9, y_initial[9]);
		gsl_vector_set(y, 10, y_initial[10]);
		gsl_vector_set(y, 11, y_initial[11]);
		gsl_vector_set(y, 12, y_initial[12]);
		gsl_vector_set(y, 13, y_initial[13]);
		gsl_vector_set(y, 14, y_initial[14]);
	
		T2 = gsl_multiroot_fsolver_hybrids;	//MAKE THE TYPE OF SOLVER A HYBRID S SOLVER
		S2 = gsl_multiroot_fsolver_alloc(T2, n2);	//MAKE THE WORKSPACE FOR SOLVER TYPE T OF DIMENSION n2
		gsl_multiroot_fsolver_set(S2, &func2, y);	//USE THE SOLVER TO SOLVE FUNCTIONS func AND FOR VARIABLES x (both have to be gsl_multiroots
	
		while(status2 = GSL_CONTINUE && iter2 < 1000){
			iter2 ++;
			status2 = gsl_multiroot_fsolver_iterate (S2);	//ITERATE THE WORKSPACE UNTIL THE THINGS IS SOLVED

			if(status2) break;				//CHECK IF IT'S WORKING. IF NOT, STOP
			status2 = gsl_multiroot_test_residual(S2->f, 1e-7);	//SET TOLERANCE FOR CONVERGENCE
		}

		ScalarField = gsl_vector_get(S2->x,0);
		VectorField = gsl_vector_get(S2->x,1);
		IsoVectorField = gsl_vector_get(S2->x,2);

		ChemicalPotentialNeutron = gsl_vector_get(S2->x,3);
		ChemicalPotentialElectron = gsl_vector_get(S2->x,4);
		ChemicalPotentialProton = ChemicalPotentialNeutron - ChemicalPotentialElectron;

		kFermiNeutron = abs(gsl_vector_get(S2->x, 5));
		kFermiProton = abs(gsl_vector_get(S2->x, 6));
		kFermiElectron = abs(gsl_vector_get(S2->x, 7));
		kFermiMuon = abs(gsl_vector_get(S2->x, 8));
		kFermiLambda = abs(gsl_vector_get(S2->x, 9));
		kFermiSigmaMinus = abs(gsl_vector_get(S2->x, 10));
		kFermiSigmaZero = abs(gsl_vector_get(S2->x, 11));
		kFermiSigmaPlus = abs(gsl_vector_get(S2->x, 12));
		kFermiXiMinus = abs(gsl_vector_get(S2->x, 13));
		kFermiXiZero = abs(gsl_vector_get(S2->x, 14));
		double kFermiDM = pow((3*pi*pi *DMDensity), (1./3.));

 		double NeutronDensity = 1./(3*Pi*Pi)*pow(kFermiNeutron, 3.);
 		double ProtonDensity = 1./(3*Pi*Pi)*pow(kFermiProton, 3.);
 		double ElectronDensity = 1./(3*Pi*Pi)*pow(kFermiElectron, 3.);
 		double MuonDensity = 1./(3*Pi*Pi)*pow(kFermiMuon, 3.);
 		double LambdaDensity = 1./(3*Pi*Pi)*pow(kFermiLambda, 3.);
		double SigmaMinusDensity = 1./(3*Pi*Pi)*pow(kFermiSigmaMinus, 3.);
		double SigmaZeroDensity = 1./(3*Pi*Pi)*pow(kFermiSigmaZero, 3.);
		double SigmaPlusDensity = 1./(3*Pi*Pi)*pow(kFermiSigmaPlus,3.);
		double XiMinusDensity = 1./(3*Pi*Pi)*pow(kFermiXiMinus, 3.);
		double XiZeroDensity = 1./(3*Pi*Pi)*pow(kFermiXiZero,3.);
		double TotalDensity = NeutronDensity + ProtonDensity + ElectronDensity + MuonDensity + LambdaDensity + SigmaMinusDensity + SigmaZeroDensity 
			+ SigmaPlusDensity + XiMinusDensity + XiZeroDensity;


		AlphaAsymmetry = ((NeutronDensity - ProtonDensity) + HadronCouplingRho*(2*SigmaMinusDensity - 2*SigmaPlusDensity + XiMinusDensity - XiZeroDensity))/BaryonDensity;

		double xScalarField = ScalarField*HadronCouplingSigma;
		double kFermi = pow((3*Pi*Pi*BaryonDensity*0.5),1./3); 
		double Mstar = NMass - ScalarField;
		double MstarLambda = LambdaMass - xScalarField;
		double MstarSigmaMinus = SigmaMinusMass - xScalarField;
		double MstarSigmaZero = SigmaZeroMass - xScalarField;
		double MstarSigmaPlus = SigmaPlusMass - xScalarField;
		double MstarXiMinus = XiMinusMass - xScalarField;
		double MstarXiZero = XiZeroMass - xScalarField;
		GCHI = Gchi(SIGCHI2, DMmass, MPI);
		GPI = Gpi(SIGCHIN, Mstar, DMmass, MPI, GCHI);
		GCHIb = Gchib(SIGCHIN, Mstar, DMmass, MPI);
		LAMBDACHIb = Lambdachi(SIGCHI2, DMmass);

		double ScalarDensity = Int_eye1(0.,kFermiNeutron, Mstar,calls) + Int_eye1(0.,kFermiProton, Mstar,calls) + HadronCouplingSigma*(eye1(kFermiLambda, MstarLambda) + 
				eye1(kFermiXiMinus, MstarXiMinus) + eye1(kFermiSigmaMinus, MstarSigmaMinus) + eye1(kFermiSigmaZero, MstarSigmaZero) + 
				eye1(kFermiSigmaPlus, MstarSigmaPlus) + eye1(kFermiXiZero, MstarXiZero));
		ScalarDensity *= 1./(Pi*Pi);
		double ScalarDensity2 = 1./(Pi*Pi)*Int_eye1(0.,kFermi, Mstar,calls);

		/*GET THE ENERGY DENSITY CONTRIBUTION FOR EACH PARTICLE*/
		double EOmega = 0.5*VectorField*VectorField/Aomega;
		double EOmega2 = 0.5*Aomega*pow((NeutronDensity + ProtonDensity) + HadronCouplingOmega*(LambdaDensity + XiMinusDensity + SigmaMinusDensity + SigmaZeroDensity 
					+ SigmaPlusDensity + XiZeroDensity), 2.);

		double ERho = 0.5*IsoVectorField*IsoVectorField/Arho;
		double ERho2 = 1./8*Arho*pow(AlphaAsymmetry*BaryonDensity, 2.);
	
		double ESigma = 0.5*ScalarField*ScalarField/Asigma;
		double ESigma2 = 0.5*Asigma*pow(ScalarDensity, 2.);
		double USigma = (1./3.)*B*NMass*pow(ScalarField, 3.) + (0.25) * C*pow(ScalarField, 4.) ; 
 
		double EkinProton = 1./(Pi*Pi) * Int_eye2(0.,kFermiProton, Mstar,calls);
		double EkinNeutron = 1./(Pi*Pi) *Int_eye2(0.,kFermiNeutron, Mstar,calls);
		double EkinElectron = 1./(Pi*Pi)*Int_eye2(0.,kFermiElectron, ElectronMass,calls);
		double EkinMuon = 1./(Pi*Pi)*Int_eye2(0.,kFermiMuon, MuonMass,calls);
		double EkinLambda = 1./(Pi*Pi)*eye2(kFermiLambda, MstarLambda);
		double EkinXiMinus = 1./(Pi*Pi)*eye2(kFermiXiMinus, MstarXiMinus);
		double EkinSigmaMinus = 1./(Pi*Pi)*eye2(kFermiSigmaMinus, MstarSigmaMinus);
		double EkinSigmaZero = 1./(Pi*Pi)*eye2(kFermiSigmaZero, MstarSigmaZero);
		double EkinSigmaPlus = 1./(Pi*Pi)*eye2(kFermiSigmaPlus, MstarSigmaPlus);
		double EkinXiZero = 1./(Pi*Pi)*eye2(kFermiXiZero, MstarXiZero);
		double KineticEnergyDM = 1./(Pi*Pi)*eye2(kFermiDM, DMmass);

		double Ekin = (EkinNeutron + EkinProton + EkinElectron + EkinMuon + EkinLambda + EkinXiMinus + EkinSigmaMinus + EkinSigmaZero + EkinSigmaPlus + EkinXiZero);
		double Ekin2 = 2./(Pi*Pi)*eye2(kFermi, Mstar);

		/*GET THE PRESSURE CONTRIBUTION FOR EACH PARTICLE*/
		double PNeutron = 1./(3*Pi*Pi)*Int_eyeP(0, kFermiNeutron, Mstar, calls);
		double PProton = 1./(3*Pi*Pi)*Int_eyeP(0, kFermiProton, Mstar, calls);
		double PElectron = 1./(3*Pi*Pi)*Int_eyeP(0, kFermiElectron, ElectronMass, calls);
		double PMuon = 1./(3*Pi*Pi)*Int_eyeP(0, kFermiMuon, MuonMass, calls);

		double PLambda = 1./(3*Pi*Pi)*Int_eyeP(0, kFermiLambda, MstarLambda, calls);
		double PXiMinus = 1./(3*Pi*Pi)*Int_eyeP(0, kFermiXiMinus, MstarXiMinus, calls);
		double PSigmaMinus = 1./(3*Pi*Pi)*Int_eyeP(0, kFermiSigmaMinus, MstarSigmaMinus, calls);
		double PSigmaZero = 1./(3*Pi*Pi)*Int_eyeP(0, kFermiSigmaZero, MstarSigmaZero, calls);
		double PSigmaPlus = 1./(3*Pi*Pi)*Int_eyeP(0, kFermiSigmaPlus, MstarSigmaPlus, calls);
		double PXiZero = 1./(3*Pi*Pi)*Int_eyeP(0, kFermiXiZero, MstarXiZero, calls);
		double PElectron2 = 1./(3*Pi*Pi)*ChemicalPotentialElectron*pow((ChemicalPotentialElectron*ChemicalPotentialElectron - ElectronMass*ElectronMass),3./2)
				- EkinElectron;
		double PMuon2 = 1./(3*Pi*Pi)*ChemicalPotentialMuon*pow((ChemicalPotentialMuon*ChemicalPotentialMuon - MuonMass*MuonMass),3./2)
				- EkinMuon;

		double PDM = 1./(3*Pi*Pi)*Int_eyeP(0, kFermiDM, DMmass, calls);

		double PFermions = (PNeutron + PProton + PElectron + PMuon + PLambda + PXiMinus + PSigmaMinus + PSigmaZero + PSigmaPlus + PXiZero);
		
/*======================================================================================================*/
		double Elambsigzero = ChemicalPotentialNeutron - HadronCouplingOmega*VectorField;
		double Esigmin = ChemicalPotentialNeutron + ChemicalPotentialElectron + HadronCouplingRho*IsoVectorField - HadronCouplingOmega*VectorField;
		double Esigplus = ChemicalPotentialProton - HadronCouplingRho*IsoVectorField - HadronCouplingOmega*VectorField;
		double Eximin = ChemicalPotentialProton + 2*ChemicalPotentialElectron + 0.5*HadronCouplingRho*IsoVectorField - HadronCouplingOmega*VectorField;
		double Exi0 = ChemicalPotentialNeutron - 0.5*HadronCouplingRho*IsoVectorField - HadronCouplingOmega*VectorField;		

/*============ ADD THEM UP ==============================================================================*/
		double EnergyDensity = Ekin + EOmega + ERho + ESigma + USigma;
		double Pressure = EOmega + ERho - ESigma - USigma + PFermions;  


/*============ DARK MATTER MODIFICATION OF EQUATION OF STATE ==============================================================================*/
		double EnergyDensityDM = 0;
		double PressureDM = 0;

		if(particle == "FERMIONS"){
			EnergyDensityDM = 0.5 * pow((GPI/MPI), 2.) * pow(BaryonDensity, 2)
			+ 0.5 * pow((GCHI/MPI), 2.) * pow(DMDensity, 2)
			+ 0.5 * (GPI/MPI)*(GCHI/MPI)* BaryonDensity*DMDensity
			+ KineticEnergyDM; 

			PressureDM = 0.5 * pow((GPI/MPI), 2.) * pow(BaryonDensity, 2)
			+ 0.5 * pow((GCHI/MPI), 2.) * pow(DMDensity, 2)
			+ 0.5 * (GPI/MPI)*(GCHI/MPI)* BaryonDensity*DMDensity
			+ PDM;
		}

		else if (particle == "BOSONS"){
			EnergyDensityDM = 0.5*MPI*MPI*GCHIb*GCHIb*pow(BaryonDensity,2.)/pow((MPI*MPI + GCHIb*GCHIb*DMDensity/DMmass), 2.)
			+ 0.5*DMmass*DMDensity
			+ pow(GCHIb, 4.)*pow(BaryonDensity,2.)*DMDensity/(DMmass*pow((MPI*MPI + GCHIb*GCHIb*DMDensity/DMmass), 2.))
			+ LAMBDACHIb/24. * DMDensity*DMDensity/pow(DMmass,2.);

			PressureDM = 0.5*MPI*MPI*GCHIb*GCHIb*pow(BaryonDensity,2.)/pow((MPI*MPI + GCHIb*GCHIb*DMDensity), 2.)
			+ 0.5*DMmass*DMmass*DMDensity
			+ pow(GCHIb, 4.)*pow(BaryonDensity,2.)*DMDensity/pow((MPI*MPI + GCHIb*GCHIb*DMDensity), 2.)
			+ LAMBDACHIb/24. * DMDensity*DMDensity;
		}

		EnergyDensity += EnergyDensityDM;
		Pressure += PressureDM;
		double BindingPerNucleon = (EnergyDensity/BaryonDensity - NMass);



/*============================== READ OUT AND OUTPUT TO CHECK =============================================*/
		DATER << BaryonDensity << " " << (BindingPerNucleon)/MeVtoinvFM << " " << Pressure/MeVtoinvFM << " " << EnergyDensity/MeVtoinvFM << " " << AlphaAsymmetry << " " << 
			Mstar/MeVtoinvFM << " " << VectorField/MeVtoinvFM << " " << -IsoVectorField/MeVtoinvFM << " " << ScalarField/MeVtoinvFM << " " << 
			ChemicalPotentialElectron/MeVtoinvFM << " " <<(ChemicalPotentialNeutron-NMass)/MeVtoinvFM<< " "<<(ChemicalPotentialProton-NMass)/MeVtoinvFM << " " << 
			NeutronDensity/BaryonDensity << " " << ProtonDensity/BaryonDensity << " " << ElectronDensity/BaryonDensity << " " << 
			LambdaDensity/BaryonDensity << " " << MuonDensity/BaryonDensity << " " << XiMinusDensity/BaryonDensity << " " << 
			SigmaMinusDensity/BaryonDensity << " " << SigmaZeroDensity/BaryonDensity << " " << SigmaPlusDensity/BaryonDensity << " " <<

			(Elambsigzero)/MeVtoinvFM << " " << MstarLambda/MeVtoinvFM
			<< " " << MstarSigmaZero/MeVtoinvFM << " " << 
			(Esigplus)/MeVtoinvFM << " " << 
			MstarSigmaPlus/MeVtoinvFM << " " << 
			(Esigmin)/MeVtoinvFM <<
			" " << MstarSigmaMinus/MeVtoinvFM << " " << MstarXiMinus/MeVtoinvFM <<
			" " << (Eximin)/MeVtoinvFM <<
			" " << MstarXiZero/MeVtoinvFM << " " << XiZeroDensity << " " << 
			(Exi0)/MeVtoinvFM << '\n'; 

		DATERclean << BaryonDensity << " " << EnergyDensity/MeVtoinvFM << " " << Pressure/MeVtoinvFM << " " << SIGCHIN*1e-26 << " " << SIGCHI2*1e-26 << " " 
			<< DMmass*1e-3/MeVtoinvFM << '\n';

		if(BaryonDensity >= 0.999*0.153 && BaryonDensity <= 1.0001*0.153){
			cout << BaryonDensity << " " << EnergyDensity << " " << BindingPerNucleon/MeVtoinvFM << " " << '\n';
			print_state2(iter2, S2);
			cout << "status2 = " << gsl_strerror(status2) << '\n';
			cout << "EKn: " << EkinNeutron << ", EKp: " << EkinProton << ", Eke: " << EkinElectron << 
				", Ekmu: " << EkinMuon << '\n' << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << '\n';
			cout << '\n' << '\n' << EOmega << " :OMEGA: " << EOmega2 << '\n' <<
				ERho << " :RHO: " << ERho2 << '\n' <<
				ESigma << " :SIGMA: " << ESigma2 << " U: " << USigma << '\n' <<
				ScalarDensity << " :nS: " << ScalarDensity2 << '\n' << 
				Ekin << " :EK: " << Ekin2 << '\n' <<  
				PElectron << " :Pe: " << PElectron2 << '\n' << "==============================" << '\n' << 
				"ne: " << ElectronDensity << '\n' << "nmu: " << MuonDensity << '\n' << "nN: " << NeutronDensity << '\n' << "nP: " << ProtonDensity << '\n' <<
				"total: " << TotalDensity/BaryonDensity << '\n' <<
				"EDensity: " << EnergyDensity << '\n' << 
				"NSage: " << NSage << '\n';

				
		}


		gsl_multiroot_fsolver_free(S2);			//CLEAR MEMORY
		gsl_vector_free(y);				//ON EVERYTHING



		BaryonDensity += 0.001;
	}

	DATER.close();
	DATERclean.close();
	return 0;
	
}

