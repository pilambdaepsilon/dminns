#include "EoSrootfinder.hh"

int main(){
	cout << "MeV->1/fm: " << MeVtoinvFM << '\n';
	int calls = 1000;
	ofstream DATER;
	ofstream DATERclean;
	DATER.open("./EoSFiles/EoSFull.dat");		//OPEN DATA FILES TO WRITE TO
	DATERclean.open("./EoSFiles/SigOm.dat");

	//SOLVE THE PROBLEM AT SATURATION FOR CONSTANTS OF THE THEORY
	SolveSaturation();

	double ChemicalPotentialNeutron = 1.; double ChemicalPotentialProton = 1.; double ChemicalPotentialElectron = 2.; double ChemicalPotentialMuon = 2.;
	ScalarField = 0.1; VectorField = 1.; IsoVectorField = 0.5; BaryonDensity = 0.001; double Asigma = SaturationVariables[3]; double Aomega = SaturationVariables[4];
	double Arho = SaturationVariables[5]; B = SaturationVariables[6]; C = SaturationVariables[7]; HadronCouplingOmega = SaturationVariables[8]; double kFermiNeutron = 1.;
	double kFermiProton = 1.; double kFermiElectron = 1.; double kFermiMuon = 1.; double kFermiLambda = 1.; double kFermiXiMinus = 1.; double kFermiSigmaMinus = 1.;
	double kFermiSigmaZero = 1.; double kFermiSigmaPlus = 1.; double kFermiXiZero = 1.;

	HadronCouplingRho = HadronCouplingSigma;
//	HadronCouplingRho = HadronCouplingOmega;
/*================= SOLVE FOR FIELD AMPLITUDES WITH KNOWN CONSTANTS FROM SATURATION ========*/
		double ych = 0.;double xch = 0.;
		double pch = 0.;
		double peach = 0.; double each = 1.;
	while (BaryonDensity <= 1.6){
		AlphaAsymmetry = 0.;
		
		const gsl_multiroot_fsolver_type *T2;		//DECLARE VARIABLE FOR TYPE OF SOLVER
		gsl_multiroot_fsolver *S2;		//DECLARE NAME FOR SOLVER WORKSPACE

		int status2;				//USED TO UPDATE ON STATUS (IS ACTUALLY BOOLEAN)
		size_t i2, iter2 =0;
	
		const size_t n2 = 15;			//DIMENSIONALITY OF PROBLEM
				//  Msig,  Mw,   Mr,  Mn,   Me,  gs, gw, gr,   b,    c,    I3,  Jn,  Jp,  Je,  ban, bap, Qe,  Qp, nB             Easym, Kcompress   Eo 
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

		double PFermions = (PNeutron + PProton + PElectron + PMuon + PLambda + PXiMinus + PSigmaMinus + PSigmaZero + PSigmaPlus + PXiZero);
		
/*======================================================================================================*/
		double Elambsigzero = ChemicalPotentialNeutron - HadronCouplingOmega*VectorField;
		double Esigmin = ChemicalPotentialNeutron + ChemicalPotentialElectron + HadronCouplingRho*IsoVectorField - HadronCouplingOmega*VectorField;
		double Esigplus = ChemicalPotentialProton - HadronCouplingRho*IsoVectorField - HadronCouplingOmega*VectorField;
		double Eximin = ChemicalPotentialProton + 2*ChemicalPotentialElectron + 0.5*HadronCouplingRho*IsoVectorField - HadronCouplingOmega*VectorField;
		double Exi0 = ChemicalPotentialNeutron - 0.5*HadronCouplingRho*IsoVectorField - HadronCouplingOmega*VectorField;		

/*============ ADD THEM UP ==============================================================================*/
		double EnergyDensity = Ekin + EOmega + ERho + ESigma + USigma;
//		cout << Ekin << " " << EOmega << " " << ERho << " " << ESigma << " " << USigma << '\n';

		double BindingPerNucleon = (EnergyDensity/BaryonDensity - NMass);

		double Pressure = EOmega + ERho - ESigma - USigma + PFermions;  
		ych = 100*xch + 120*xch*xch - 0*pow(xch, 3.);
		pch = 0.1*xch + 220*xch*xch - 0*pow(xch, 3.);
		double kappa = 20.428/pow(6.242e5,2./3);
		double gammer = 5./3;
		peach = kappa*pow(each, gammer); 
		xch += 0.01;
		each *= 2.01;

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
			(Exi0)/MeVtoinvFM << " " << xch << " " << ych << " " << pch << " " << each << " " << peach <<  
 			'\n'; 

		DATERclean << BaryonDensity << " " << EnergyDensity/MeVtoinvFM << " " << Pressure/MeVtoinvFM << '\n';

		if(BaryonDensity >= 0.999*0.7 && BaryonDensity <= 1.0001*0.7){
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
				"PLambda: " << PLambda << '\n';
				
		}


		gsl_multiroot_fsolver_free(S2);			//CLEAR MEMORY
		gsl_vector_free(y);				//ON EVERYTHING



		BaryonDensity += 0.001;
	}

	DATER.close();
	DATERclean.close();
	return 0;
	
}

