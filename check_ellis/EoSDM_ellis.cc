#include "EoSrootfinder.hh"

int main(int argc, char* argv[]){
	double GCHI_MPI2 = atof(argv[1]);    //eta
//	double SIGCHIN = atof(argv[2]);    //DM nucleon cross-section
//	double SIGCHI2 = atof(argv[1]);    //DM self interaction cross-section [cm^2]
	double DMmass    = atof(argv[2]);  //DM mass
//	string EoSMODEL  = argv[5];        //Nuclear EOS model
  //      string PARTICLE  = argv[6];        //DM particle species
        double SIGCHI2 = pow(GCHI_MPI2*DMmass*convinvGeVtocm,2.0);
        cout << "Self-interaction cross-section: " << SIGCHI2 << '\n';

	stringstream ss;
	cout << "DM Properties: M[GeV]: " << DMmass << "   (self-coupling/mediator mass)^2 [GeV^-2]: " << GCHI_MPI2 << "   Cross-section[cm^2]: " << SIGCHI2 << '\n';
	ss << GCHI_MPI2 << "_" << DMmass << "_" << SIGCHI2;
	if(SIGCHI2/DMmass >= 2e-24){
		cout << "Cross-section/mass combination violates Bullet Cluster bound [2e-24 cm^2/GeV] " << SIGCHI2/DMmass << '\n';
		return 0;
	}
	string fno = ss.str();
/*====================== SET UP INITIAL CONDITIONS ==========================*/
	double NMass = 0.938;

/*====================== BUILD THE DARK MATTER EOS ==========================*/
        double EKinetic_DM = 0.0;	//DM kinetic + rest mass energy
	double PKinetic_DM = 0.0;	//DM kinetic + rest mass pressure
	double EnergyDensity = 0.0;	//total DM energy density
	double Pressure = 0.0;	//total DM pressure
	double EnergyDensityDM = 0.0;	//total DM energy density
	double PressureDM = 0.0;	//total DM pressure
	double DMDensity = 0.001;		//DM Density [fm^-3]
	double DMDensMax = 100;		//Max Density to build sequence
	int calls = 1000;               //Calls to integrator
//	GCHI_MPI2 = sqrt(SIGCHI2*pow(convcmtoinvGeV,2.0)) / DMmass;
        //open file to write DM EOS in
	stringstream liltitle1b;
	liltitle1b << "./EOSDM_" << fno << ".dat";
	string BIGTITLEb = liltitle1b.str();
	ofstream DATERclean;
	DATERclean.open(BIGTITLEb);
	while (DMDensity < DMDensMax){
		
/*============ USE FIT TO FIND UNMODIFIED EOS ==============================================================================*/
//		double EnergyDensity = gsl_spline_eval(splineE, BaryonDensity, accE);
//		double Pressure = gsl_spline_eval(splineP, BaryonDensity, accP);

/*============ DARK MATTER MODIFICATION OF EQUATION OF STATE ==============================================================================*/
		double kFermiDM = pow((3*pi*pi *DMDensity*1e39*pow(convinvcmtoGeV,3.0)), (1./3.));			//DM Fermi momementum
		PKinetic_DM = 1./(3*Pi*Pi)*Int_eyeP(0, kFermiDM, DMmass, calls);	//using integral definition of P_kin
		EKinetic_DM = 1./(Pi*Pi)*eye2(kFermiDM, DMmass);			//using analytic result for E_kin

		EnergyDensityDM = 0.5 * GCHI_MPI2 * pow(DMDensity*1e39*pow(convinvcmtoGeV,3.0), 2.0);
//		+ 0.5 * pow((GCHI/MPI), 2.) * pow(DMDensity, 2)
//		+ 0.5 * (GPI/MPI)*(GCHI/MPI)* BaryonDensity*DMDensity
//		+ EKinetic_DM; 

		PressureDM = 0.5 * GCHI_MPI2 * pow(DMDensity*1e39*pow(convinvcmtoGeV,3.0), 2.0);
//		+ 0.5 * pow((GCHI/MPI), 2.) * pow(DMDensity, 2)
//		+ 0.5 * (GPI/MPI)*(GCHI/MPI)* BaryonDensity*DMDensity
//		+ PKinetic_DM;
//		cout << DMDensity << " " << PKinetic_DM << " " << EKinetic_DM << " " << EnergyDensityDM << " " << PressureDM << " " << GCHI_MPI2 <<  '\n';

/*================= ADD THEM UP ========*/
		EnergyDensity = (EnergyDensityDM + EKinetic_DM)*1e3*pow(convGeVtoinvcm,3.0)*1e-39;
		Pressure = (PressureDM + PKinetic_DM)*1e3*pow(convGeVtoinvcm,3.0)*1e-39;

		DATERclean << DMDensity << " " << EnergyDensity << " " << Pressure << '\n'; 
//			<< SIGCHIN << " " << SIGCHI2 << " " << DMmass*1e-3/MeVtoinvFM << '\n';
/*============================== READ OUT AND OUTPUT TO CHECK =============================================*/


		DMDensity *= 1.01;
	}

	return 0;
	
}

