#include "dmnofunc.hh"
#include <vector>

int main(int argc, char* argv[]){
/*================================================================================================================================================
  _____  
 |  __ | 
 | |__)|
 |  ___/ 
 | |     
 |_| ull in cross-sections and other relevant parameters from user input and from OPE code
==================================================================================================================================================*/
	double MNS = 0;
	double RNS = 0;
	double coefficient = atof(argv[1]);
	double powsigchin = atof(argv[2]);
	double powsigchi2 = atof(argv[3]);
	double DMmass = atof(argv[4]);
	string STAR = argv[5];
	string EoSMODEL = argv[6];

	double title1 = atof(argv[1]);
	double title2 = atof(argv[2]);
	double title3 = atof(argv[3]);
	double title4 = atof(argv[4]);
	int fno;
	string fnos;
	string particle = "BOSONS";
/*====================== TAKE IN CONFIGURATION FILE TO GET STAR AND MODEL DETAILS ==========================*/
	string starname; double dNS_Earth; double dNS_GC; double NSage;

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
	rhoDM = stof(storevalues[4])/pow(convcmtoinvGeV,3.);
	cout << "Star: " << starname << '\n' <<
	       	"Distance to Earth: " << dNS_Earth << '\n' <<
		"Distance to Galactic Center: " << dNS_GC << '\n' << 
		"Age: " << NSage << '\n' << 
		"Ambient DM Density [NFW]: " << rhoDM*pow(convcmtoinvGeV,3.) << '\n';

/*================================= TAKE IN THE MR SEQUENCE WITHOUT ANY DM ===========================*/
	vector< vector<double> > MRDATA;
	double pr1, pr2, pr3, pr4;
	ifstream dater;
	stringstream liltitle1;
	liltitle1 << "./" << "/MRSequences/MR_" << EoSMODEL << ".dat";

	string BIGTITLE1 = liltitle1.str();
	dater.open(BIGTITLE1);
	if (!dater){
		cout << '\n' << "MR Sequence File Not Found... EXITING" << '\n';
		return 1;
	}
	int columns = 0;
	while(!dater.eof()){
		dater >> pr1 >> pr2 >> pr3;
		vector<double> col;
		col.push_back(pr1*Msolar);
		col.push_back(pr2*1e5*convcmtoinvGeV);
		col.push_back(pr3*1e39*pow(convinvcmtoGeV,3.));
		MRDATA.push_back(col);
		columns = col.size();
		col.clear();
	}

	int rows = MRDATA.size();
	dater.close();

	cout << "EoS Model: " << EoSMODEL << '\n';
	
	cout << '\n';
	cout << rows << '\n';
	double SIGCHIN = coefficient*pow(10., -1.*powsigchin);
//2	double SIGCHI2 = 2*DMmass*pow(10., -24);
	double SIGCHI2 = pow(10.,-1.*powsigchi2);
	
	if(SIGCHI2 > 2.*DMmass*1e-24){
		cout << "Mass/Self-interaction Combination breaks Bullet Cluster Limit, Try again \n";
		return 1;
	}

	cout << SIGCHIN << " " << SIGCHI2 << " " << DMmass << '\n';

	double ReducedMass = DMmass * Nmass/(DMmass + Nmass);

	double dummyscan = 0.;
	stringstream ss;
	ss << title1 << title2 << title3 << title4;
	fnos = ss.str();
	cout << fnos << '\n';

	stringstream liltitle2;
	liltitle2 << "./" << particle << "/DmNoFiles/" << title4 << "GeV/nchi" << fnos << "pSKINNY.dat";
	string BIGTITLE2 = liltitle2.str();
	ofstream myfile2;						//my file
	myfile2.open(string(BIGTITLE2));
	double TIME = 10;
	double DMNO = 0;
	double DMNOplus = 0;
	
	double CoreDens = 0.01*1e39*pow(convinvcmtoGeV, 3.);
	double FinalCoreDens = 1.6*1e39*pow(convinvcmtoGeV, 3.);
	double AvDens = 0;
	double CaptureRate = 0.;
	double SelfCapture = 0.;
	double DensWeightFactor = 0.557095;		//obtained as weight of central density in average of upper lim on central dens and vacuum...
while (CoreDens <= FinalCoreDens){
	AvDens = CoreDens*DensWeightFactor;	
	for(int scan = 0; scan < rows; scan++){
		dummyscan = MRDATA[scan][2];
		double relerrDENS = fabs(dummyscan - CoreDens)/max(dummyscan, CoreDens);
		if(relerrDENS <= 1e-3){
			MNS = MRDATA[scan][0];
			RNS = MRDATA[scan][1];
		}
	}
	MNS = MNS/Msolar;
	RNS = RNS/(1e5*convcmtoinvGeV);
	rhoDM *= pow(convcmtoinvGeV,3.);
	pFermi = pow((3.*pi*pi*AvDens),(1./3));
	double PauliNeutrons = PauliBlocking(ReducedMass, pFermi, MNS*Msolar, RNS*1e5*convcmtoinvGeV);
	CaptureRate = (7.43e31)*rhoDM/DMmass * MNS*RNS/(1 - 2.94*MNS/RNS)*min(3.78e-9*(SIGCHIN/1e-55)*MNS/(RNS*RNS),1.)*PauliNeutrons;
	SelfCapture = (4.65e-4)*rhoDM/DMmass*(SIGCHI2/1e-24)*(MNS/RNS)/(1-2.94*MNS/RNS);

	MNS = MNS*Msolar;
	RNS = RNS*(1e5*convcmtoinvGeV);
	rhoDM /= pow(convcmtoinvGeV,3.);
	TIME = 10;
	DMNO = 0;
	DMNOplus = 0;
/*================================================================================================================================================
   _____ 
  / ____|
 | |  __ 
 | | |_ |
 | |__| |
  \_____| et the relevant timescales for DM capture in the NS
==================================================================================================================================================*/

	double TCONTAIN = tcon(DMmass, SIGCHIN, MNS, RNS);
	if (TCONTAIN > 3.5){
		TCONTAIN = 3.5;
	}
	double TGEO = tgeo(DMmass, ReducedMass, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture, AvDens, RNS);
	double TTHERM = ttherm(DMmass, Nmass, ReducedMass, SIGCHIN, AvDens);
	double times[4] = {TGEO, TTHERM, TCONTAIN, TMAX};
	double t1 = 0;double t2=0;double t3 = 0;double t4 = 0;
	int counter = 0;
	vector<double> UNSAFExsecschin;
	vector<double> UNSAFEmasses;
	vector<double> UNSAFExsecschi2;

	vector<double> timevec (times, times+4);
	double RTH = radius(DMmass, ReducedMass, TTHERM, SIGCHIN, AvDens, RNS)/convcmtoinvGeV;
	double RTHEST = sqrt(9.* TNS/(4.*pi*GNewton0*AvDens*Nmass*DMmass))/convcmtoinvGeV;
/*================= Sort the time scales for proper evolution =====================================*/
	stable_sort(timevec.begin(), timevec.end());
	t1 = timevec[0];
	t2 = timevec[1];
	t3 = timevec[2];
	t4 = timevec[3];
	double tprev = 0;
	int count = 0;
//	cout << timevec[0] << " " << timevec[1] << " " << timevec[2] << " " << timevec[3] << '\n';
//	cout << "TCONTAIN: " << TCONTAIN << ", TGEO: " << TGEO << ", TTHERM: " << TTHERM << ", TMAX: " << TMAX << '\n';
/*================================================================================================================================================
 _____ 
/ ____|
| |     
| |     
| |____ 
\______| alculate the number of DM particles in the star by regime of capture
==================================================================================================================================================*/
	int chandcountb = 0;
	int chandcountf = 0;
	while(TIME <= TMAX){

		if (TIME < t1){
			DMNOplus = 0;
		}
 
		else if (TIME >= t1 && TIME < t2){
			if(t1 == TGEO){
				DMNOplus = Nchigeo(DMmass, ReducedMass, TGEO, TIME, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture, AvDens, RNS) - Nchigeo(DMmass, ReducedMass, TGEO, tprev, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture, AvDens, RNS);
				if (SIGCHI2 <= 1e-45){
					DMNOplus = CaptureRate*(TIME - tprev);
				}
				if(count == 0){
//				cout << "T1: Geometric" << '\n';
				count = 1;}
			}
			else if(t1 == TTHERM){
				DMNOplus = Nchitherm(DMmass, ReducedMass, TTHERM, TIME, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture, AvDens) - Nchitherm(DMmass, ReducedMass, TTHERM, tprev, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture, AvDens);
				if (SIGCHI2 <= 1e-45){
					DMNOplus = CaptureRate*(TIME - tprev);
				}
				if(count == 0){
//				cout << "T1: Thermalization" << '\n';
				count = 1;}
			}
			else if(t1 == TMAX){
				DMNOplus = 0;
//				cout << "DONE, t1 = 10^10 years." << '\n';
			}
			else{
				DMNOplus = Nchi0(CaptureRate, SelfCapture, TIME) - Nchi0(CaptureRate, SelfCapture, tprev);
				if (SIGCHI2 <= 1e-45){
					DMNOplus = CaptureRate*(TIME - tprev);
				}
				if(count == 0){
//				cout << "T1: Regular Evolution" << '\n';
				count = 1;}
			}
		}

		else if (TIME >= t2 && TIME < t3){
			if(t2 == TGEO){
				DMNOplus = Nchigeo(DMmass, ReducedMass, TGEO, TIME, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture, AvDens, RNS) - Nchigeo(DMmass, ReducedMass, TGEO, tprev, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture, AvDens, RNS);
				if (SIGCHI2 <= 1e-45){
					DMNOplus = CaptureRate*(TIME - tprev);
				}
				if(count == 1){
//				cout << "T2: Geometric" << '\n';
				count = 2;}
			}
			else if(t2 == TTHERM){
				DMNOplus = Nchitherm(DMmass, ReducedMass, TTHERM, TIME, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture, AvDens) - Nchitherm(DMmass, ReducedMass, TTHERM, tprev, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture, AvDens);
				if (SIGCHI2 <= 1e-45){
					DMNOplus = CaptureRate*(TIME - tprev);
				}
				if(count == 1){
//				cout << "T2: Thermalization" << '\n';
				count = 2;}
			}
			else if(t2 == TMAX){
				DMNOplus = 0;
//				cout << "DONE, t2 = 10^10 years" << '\n';
			}
			else{
				DMNOplus = Nchi0(CaptureRate, SelfCapture, TIME) - Nchi0(CaptureRate, SelfCapture, tprev);
				if (SIGCHI2 <= 1e-45){
					DMNOplus = CaptureRate*(TIME - tprev);
				}
				if(count == 1){
//				cout << "T2: Regular Evolution" << '\n';
				count = 2;}
			}
		}

		else if (TIME >= t3 && TIME <t4){
			if(t3 == TGEO){
				DMNOplus = Nchigeo(DMmass, ReducedMass, TGEO, TIME, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture, AvDens, RNS) - Nchigeo(DMmass, ReducedMass, TGEO, tprev, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture, AvDens, RNS);
				if (SIGCHI2 <= 1e-45){
					DMNOplus = CaptureRate*(TIME - tprev);
				}
				if(count == 2){
//				cout << "T3: Geometric" << '\n';
				count = 3;}
			}
			else if(t3 == TTHERM){
				DMNOplus = Nchitherm(DMmass, ReducedMass, TTHERM, TIME, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture, AvDens) - Nchitherm(DMmass, ReducedMass, TTHERM, tprev, SIGCHIN, SIGCHI2, CaptureRate, SelfCapture, AvDens);
				if (SIGCHI2 <= 1e-45){
					DMNOplus = CaptureRate*(TIME - tprev);
				}
				if(count == 2){
//				cout << "T3: Thermalization" << '\n';
				count = 3;}
			}
			else if(t3 == TMAX){
				DMNOplus = 0;
//				cout << "DONE, t3 = 10^10 years" << '\n';
			}
			else{
				DMNOplus = Nchi0(CaptureRate, SelfCapture, TIME) - Nchi0(CaptureRate, SelfCapture, tprev);
				if (SIGCHI2 <= 1e-45){
					DMNOplus = CaptureRate*(TIME - tprev);
				}
				if(count == 2){
//				cout << "T3: Regular Evolution" << '\n';
				count = 3;}
			}
		}

		else if (TIME >= t4){
			DMNOplus = 0;
			if(counter == 0){
//				cout << "DONE, t4 = 10^10 years" << '\n';
				counter++;
			}
		}

		else{
			DMNOplus = 0;
//			cout << "NO EVOLUTION" << '\n';
		}


		//if (DMNOplus != DMNOplus){ DMNOplus = 0.0;}	
		DMNO += DMNOplus;
		double Rchi = radius(DMmass, ReducedMass, TIME, SIGCHIN, AvDens, RNS)*convinvGeVtocm*1e13;		//DM sphere radius in fm
		double Volumechi = 4./3. *pi * pow(Rchi, 3.);
		double relerrTIME = fabs(TIME - NSage)/max(TIME,NSage);
		double DMDENS = DMNO/(Volumechi);			//DM number density in inverse fm^3
		if(particle == "FERMIONS"){
			double pFDM = pow((3*pi*pi*DMDENS),(1./3));
			double PauliDM = PauliBlocking(DMmass/2., pFDM, MNS, RNS);
			SelfCapture = (4.65e-4)*rhoDM/DMmass*(SIGCHI2/1e-24)*(MNS/RNS)/(1-2.94*MNS/RNS)*PauliDM;
		}
		if (isfinite(DMNO) and isfinite(DMDENS)){
			if (relerrTIME <= 5e-2){
				myfile2 << TIME << " " << DMNO << " " << DMDENS*1e-39*pow(convGeVtoinvcm,3.) << " " << CoreDens*1e-39*pow(convGeVtoinvcm,3.) << " " 
				<< SIGCHIN << " " << SIGCHI2 << " " << DMmass << '\n';
//				cout << TIME << " " << DMDENS*1e-39*pow(convGeVtoinvcm,3.) << " " << CoreDens*1e-39*pow(convGeVtoinvcm,3.) << 
//					" " << MNS/Msolar << " " << RNS/(1e5*convcmtoinvGeV) << '\n';
			}
		}
		tprev = TIME;
		TIME *= 1.1;
	}
//	cout << CoreDens << '\n';
	CoreDens *= 1.01;
	timevec.clear();
	}
	myfile2.close();
	MRDATA.clear();
	return 0;
}
