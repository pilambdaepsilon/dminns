#include "EoSintegrals.hh"
#include <complex>
double SaturationVariables[9] = {0.0};

double BaryonDensity = 0.153;
double ScalarMass = 550*MeVtoinvFM;
double VectorMass = 783*MeVtoinvFM;
double IsoVectorMass = 800*MeVtoinvFM;
double NMass = 938*MeVtoinvFM;
double LambdaMass = 1115*MeVtoinvFM;
double XiMinusMass = 1315*MeVtoinvFM;
double SigmaMinusMass = 1190*MeVtoinvFM;
double SigmaZeroMass = 1190*MeVtoinvFM;
double SigmaPlusMass = 1190*MeVtoinvFM;
double XiZeroMass = 1315*MeVtoinvFM;
double ElectronMass = 0.511*MeVtoinvFM;
double MuonMass = 105.7*MeVtoinvFM;
double HadronCouplingSigma = 0.6;
double HadronCouplingRho = 0.;
double HadronCouplingOmega = 0.;
double ScalarField = 0.0;
double VectorField = 0.0;
double IsoVectorField = 0.0;
double B = 0.0;
double C = 0.0;
double AsymmetryEnergy = 32.5*MeVtoinvFM;
double Compressibility = 240*MeVtoinvFM;
double BindingPerNucleon0 = -16.3*MeVtoinvFM;
double MstarPerMn0 = 0.78;
double AlphaAsymmetry = 0.;			// Aplha = 0 is for symmetric matter 

//DEFINE FUNCTION FOR SATURATION CONSTANTS
int GM1_f(const gsl_vector *x, void *params, gsl_vector *f){
	//SET PARAMETERS
	double Mn = ((struct Rparams *) params)->Mn;
	double Me = ((struct Rparams *) params)->Me;
	double nB = ((struct Rparams *) params)->nB;
	double Easymm = ((struct Rparams *) params)->Easymm;
	double Kcompress = ((struct Rparams *) params)->Kcompress;
	double BperA = ((struct Rparams *) params)->BperA;
	double gamma = ((struct Rparams *) params)->gamma;
	double asymmetryfactor = 0;
	double xHadronsigma = ((struct Rparams *) params)->xHadronsigma;
	
	//SET VARIABLES TO BE SOLVED FOR
	double gsigma = gsl_vector_get(x,0);
	double gomega = gsl_vector_get(x,1);
	double grho = gsl_vector_get(x,2);
	double asigma = gsl_vector_get(x,3);
	double aomega = gsl_vector_get(x,4);
	double arho = gsl_vector_get(x,5);
	double b = gsl_vector_get(x,6);
	double c = gsl_vector_get(x,7);

	double xOmega = gsl_vector_get(x,8);

	//SET FUNCTIONS TO BE SOLVED
	double mstar = Mn*gamma;
	double kF = pow((3*Pi*Pi*nB*0.5),1./3.); 
	double kFactor = sqrt(kF*kF + mstar*mstar);	
	int calls = 1000;
	double Integral1 = eye1(kF, mstar); 
	double Integral2 = eye2(kF, mstar); 
	double Integral3 = eye3(kF, mstar); 
	double Pref7 = 1./asigma + 2*b*Mn*gsigma + 3*c*pow(gsigma,2.) + 2./(Pi*Pi)*Integral3;
	
	double f0 = gomega - aomega*nB;												//*
	double f1 = grho + 0.5*arho*asymmetryfactor*nB;										//*
	double f2 = gsigma + b* Mn*asigma*pow(gsigma, 2.) + c*asigma*pow(gsigma,3.) - 2./(Pi*Pi)*asigma*Integral1;		//*

	double f3 = nB*aomega - Mn - BperA + kFactor;
	double f4 = gsigma - (1 - gamma)*Mn;											//*
 	double f5 = nB*(Mn +BperA) - 2./(Pi*Pi)*Integral2 - (1./3.)*b*Mn*pow(gsigma,3.) - 0.25*c*pow(gsigma,4.) - 0.5*gsigma*gsigma/asigma - 0.5*aomega*nB*nB;	//*

	double f6 = Easymm - arho*nB/(8.) - kF*kF/(6*kFactor);													//*
 	double f7 = Kcompress - 9*aomega*nB - 3*kF*kF/kFactor + 6*pow(kF,3.)/(Pi*Pi)*mstar*mstar/pow(kFactor, 2.)/Pref7; 			//

	double f8 = -28.*MeVtoinvFM - xOmega*gomega + xHadronsigma *gsigma;

	gsl_vector_set(f, 0, f0);
	gsl_vector_set(f, 1, f1);
	gsl_vector_set(f, 2, f2);
	gsl_vector_set(f, 3, f3);
	gsl_vector_set(f, 4, f4);
	gsl_vector_set(f, 5, f5);
	gsl_vector_set(f, 6, f6);
	gsl_vector_set(f, 7, f7);
	gsl_vector_set(f, 8, f8);
	
	return GSL_SUCCESS;
}

//DEFINE FUNCTION FOR EVOLUTION with nB
int GM2_f(const gsl_vector *x, void *params, gsl_vector *f){
	//SET PARAMETERS
	double Mn = ((struct Rparams *) params)->Mn;
	double Me = ((struct Rparams *) params)->Me;
	double Mmu = ((struct Rparams *) params)->Mmu;
	double as = ((struct Rparams *) params)->gs;
	double aw = ((struct Rparams *) params)->gw;
	double ar = ((struct Rparams *) params)->gr;
	double bethe = ((struct Rparams *) params)->bethe;
	double cethe = ((struct Rparams *) params)->cethe;
	double nB = ((struct Rparams *) params)->nB;
	double asymmetryfactor = ((struct Rparams *) params)->asymmetryfactor;
	double Mlambda = ((struct Rparams*) params)->Mlambda;
	double xHadronsigma = ((struct Rparams *) params)->xHadronsigma;
	double xHadronomega = ((struct Rparams *) params)->xHadronomega;
	double xHadronrho = ((struct Rparams *) params)->xHadronrho;
	double Mximinus = ((struct Rparams *) params)->Mximinus;
	double Msigmaminus = ((struct Rparams *) params) ->Msigmaminus;
	double Msigmazero = ((struct Rparams *) params) ->Msigmazero;
	double Msigmaplus = ((struct Rparams *) params)->Msigmaplus;
	double Mxizero = ((struct Rparams *) params)->Mxizero;

	//SET VARIABLES TO BE SOLVED FOR
	double gsigma = gsl_vector_get(x,0);
	double gomega = gsl_vector_get(x,1);
	double grho = gsl_vector_get(x,2);

	double mun = gsl_vector_get(x,3);
	double mue = gsl_vector_get(x,4);

	double kFn = gsl_vector_get(x,5);
	double kFp= gsl_vector_get(x,6);
	double kFe = gsl_vector_get(x,7);
	double kFmu = gsl_vector_get(x,8);
	double kFlambda = gsl_vector_get(x, 9);
	double kFsigmaminus = gsl_vector_get(x,10);
	double kFsigmazero = gsl_vector_get(x,11);
	double kFsigmaplus = gsl_vector_get(x,12);
	double kFximinus = gsl_vector_get(x,13);
	double kFxizero = gsl_vector_get(x,14); 

	//SET FUNCTIONS TO BE SOLVED
	double xgsigma = xHadronsigma*gsigma;
	double xgomega = xHadronomega*gomega;
	double xgrho = xHadronrho*grho;
	double Mnstar = Mn - gsigma;
	double Mlambdastar = Mlambda - xgsigma;
	double Msigmaminusstar = Msigmaminus - xgsigma;
	double Msigmazerostar = Msigmazero - xgsigma;
	double Msigmaplusstar = Msigmaplus - xgsigma;
	double Mximinusstar = Mximinus - xgsigma;
	double Mxizerostar = Mxizero - xgsigma;

	double nN = 1./(3*Pi*Pi)*pow(kFn,3.);
	double nP = pow(kFp,3.)*1./(3*Pi*Pi);
	double nE = pow(kFe,3.)*1./(3*Pi*Pi);
	double nMu = pow(kFmu,3.)*1./(3*Pi*Pi);
	double nLambda = pow(kFlambda,3.)*1./(3*Pi*Pi);
	double nXiminus = pow(kFximinus, 3.)*1./(3*Pi*Pi);
	double nSigmaminus = pow(kFsigmaminus, 3.)*1./(3*Pi*Pi);
	double nSigmazero = pow(kFsigmazero, 3.)*1./(3*Pi*Pi);
	double nSigmaplus = pow(kFsigmaplus, 3.)*1./(3*Pi*Pi);
	double nXizero = pow(kFxizero, 3.)*1./(3*Pi*Pi);

	double mup = mun - mue;
	double mulambda = mun;
	double musigmaminus = mun + mue;
	double musigmazero = mun;
	double musigmaplus = mun - mue;
	double muximinus = mun + mue;
	double muxizero = mun;

	double EFelectron = mue;
	double EFmuon = mue;
	double EFlambda  = mulambda - xgomega;
	double EFximinus = muximinus - xgomega + 0.5*xgrho;
	double EFsigmaminus = musigmaminus - xgomega + xgrho;
	double EFsigmazero = musigmazero - xgomega;
	double EFsigmaplus = musigmaplus - xgomega - xgrho;
	double EFxizero = muxizero - xgomega - 0.5*xgrho;
	double nCharge = nP - nE - nMu - nSigmaminus + nSigmaplus;
//	if(nCharge >= 1e-3){
//		cout << '\n' << '\n' << "CHECK!" << '\n' << "Baryons: " << nB << ", Protons: " << nP << ", Electrons: " << nE << 
//			", Muons: " << nMu << ", Sigmas(-) " << nSigmaminus << ", Total: " << nCharge << '\n';}

	double Integral1 = eye1(kFn, Mnstar) + eye1(kFp, Mnstar) + xHadronsigma*(eye1(kFlambda, Mlambdastar) + eye1(kFximinus, Mximinusstar) 
	+ eye1(kFsigmaminus, Msigmaminusstar) + eye1(kFsigmazero, Msigmazerostar) + eye1(kFsigmaplus, Msigmaplusstar) + eye1(kFxizero, Mxizerostar)); 

	double f0 = gomega - aw*(nP + nN) - aw*xHadronomega*(nLambda + nXiminus + nXizero + nSigmaminus + nSigmazero + nSigmaplus);
	double f1 = grho + 0.5*ar*(nN - nP) + 0.5*ar*xHadronrho*(nXiminus + 2*nSigmaminus - 2*nSigmaplus - nXizero);
	double f2 = gsigma + bethe* Mn*as*pow(gsigma, 2.) + cethe*as*pow(gsigma,3.) - 1./(Pi*Pi)*as*Integral1;

	double f3 = mun - sqrt(kFn*kFn + Mnstar*Mnstar) - gomega + 0.5*grho;
	double f4 = mup - sqrt(kFp*kFp + Mnstar*Mnstar) - gomega - 0.5*grho;
	double f5 = mue - sqrt(kFe*kFe + Me*Me);

	double f6 = nB - nN - nP - nLambda - nXiminus - nSigmaminus - nSigmazero - nSigmaplus - nXizero;		//BARYON DENSITY
	double f7 = nP - nE - nMu - nSigmaminus + nSigmaplus;				//CHARGE NEUTRALITY

	double f8 = kFmu;
	double f9 = kFlambda;										//ZERO MOMENTUM  UNTIL POPULATED
	double f10 = kFsigmaminus;
	double f11 = kFsigmazero;
	double f12 = kFsigmaplus;
	double f13 = kFximinus;
	double f14 = kFxizero;

	//Solve for particle species parameters after it has been populated
	if (EFmuon > Mmu){ f8 = EFmuon - sqrt(kFmu*kFmu + Mmu*Mmu);}
	if (EFlambda > Mlambdastar){ f9 = EFlambda - sqrt(kFlambda*kFlambda + Mlambdastar*Mlambdastar);}
	if (EFsigmaminus > Msigmaminusstar){ f10 = EFsigmaminus - sqrt(kFsigmaminus*kFsigmaminus + Msigmaminusstar*Msigmaminusstar);}
	if (EFsigmazero > Msigmazerostar){ f11 = EFsigmazero - sqrt(kFsigmazero*kFsigmazero + Msigmazerostar*Msigmazerostar);}
	if (EFsigmaplus > Msigmaplusstar){ f12 = EFsigmaplus - sqrt(kFsigmaplus*kFsigmaplus + Msigmaplusstar*Msigmaplusstar);}
	if (EFximinus > Mximinusstar){ f13 = EFximinus - sqrt(kFximinus*kFximinus + Mximinusstar*Mximinusstar);}
	if (EFxizero > Mxizerostar){ f14 = EFxizero - sqrt(kFxizero*kFxizero + Mxizerostar*Mxizerostar);}

	gsl_vector_set(f, 0, f0);
	gsl_vector_set(f, 1, f1);
	gsl_vector_set(f, 2, f2);
	gsl_vector_set(f, 3, f3);
	gsl_vector_set(f, 4, f4);
	gsl_vector_set(f, 5, f5);
	gsl_vector_set(f, 6, f6);
	gsl_vector_set(f, 7, f7);
	gsl_vector_set(f, 8, f8);
	gsl_vector_set(f, 9, f9);
	gsl_vector_set(f, 10, f10);
	gsl_vector_set(f, 11, f11);	
	gsl_vector_set(f, 12, f12);
	gsl_vector_set(f, 13, f13);
	gsl_vector_set(f, 14, f14);
	return GSL_SUCCESS;
}


int print_state(size_t iter, gsl_multiroot_fsolver *S){
cout << "iter = " << iter << '\n' << 
	" x = " << "gsigma: " << gsl_vector_get(S->x, 0)/MeVtoinvFM << ", gomega: " << gsl_vector_get(S->x, 1)/MeVtoinvFM << ", grho: " << gsl_vector_get(S->x, 2)/MeVtoinvFM 
	<< ", asigma: " << gsl_vector_get(S->x, 3) << ", aomega: " <<
	gsl_vector_get(S->x, 4) << ", arho: " << gsl_vector_get(S->x, 5) << ", b: " << gsl_vector_get(S->x, 6) << ", c: " << gsl_vector_get(S->x, 7) <<
       ", xomega: " << gsl_vector_get(S->x, 8) << '\n' << 
	" F(x) = " << gsl_vector_get(S->f, 0) << ", " << gsl_vector_get(S->f, 1) << ", " << gsl_vector_get(S->f, 2) << ", " << gsl_vector_get(S->f, 3) << ", " <<
	gsl_vector_get(S->f, 4) << ", " << gsl_vector_get(S->f, 5) << ", " << gsl_vector_get(S->f, 6) << ", " << gsl_vector_get(S->f, 7) << ", " << gsl_vector_get(S->f, 8) 
	<< '\n' <<
	 "============================================================" << '\n'; 
}
int print_state2(size_t iter2, gsl_multiroot_fsolver *S2){
cout << "iter = " << iter2 << '\n' << 
	" x = " << "gsigma: " << gsl_vector_get(S2->x, 0)/MeVtoinvFM << ", gomega: " << gsl_vector_get(S2->x, 1)/MeVtoinvFM << ", grho: " <<
       	gsl_vector_get(S2->x, 2)/MeVtoinvFM << ", MuN: " << gsl_vector_get(S2->x, 3)/MeVtoinvFM << ", MuE: " << gsl_vector_get(S2->x, 4)/MeVtoinvFM <<
       	", kFn: " << gsl_vector_get(S2->x, 5) << ", kFp: " << gsl_vector_get(S2->x, 6) << ", kFe: " << gsl_vector_get(S2->x, 7) << 
	", kFmu: " << gsl_vector_get(S2->x,8) << ", kFlambda: " << gsl_vector_get(S2->x, 9) << ", kFSigma(-): " << gsl_vector_get(S2->x,10) << 
	", kFSigma(0): " << gsl_vector_get(S2->x,11) << 
	", kFSigma(+): " << gsl_vector_get(S2->x, 12) << ", kFXi(-): " << gsl_vector_get(S2->x, 13) << ", kFXi(0): " << gsl_vector_get(S2->x,14) 
	<< '\n' <<
	" F(x) = " << gsl_vector_get(S2->f, 0) << ", " << gsl_vector_get(S2->f, 1) << ", " << gsl_vector_get(S2->f, 2) << ", "  << gsl_vector_get(S2->f,3) << ", " <<
	gsl_vector_get(S2->f,4) << ", " << gsl_vector_get(S2->f, 5) << ", " << gsl_vector_get(S2->f,6) << ", " << gsl_vector_get(S2->f,7) << ", " << 
	gsl_vector_get(S2->f, 8) << ", " << gsl_vector_get(S2->f, 9) << ", " << gsl_vector_get(S2->f, 11) << ", " << gsl_vector_get(S2->f, 12) << ", " << 
	gsl_vector_get(S2->f, 13) << ", " << gsl_vector_get(S2->f, 14)
	<< '\n' <<
	 "============================================================" << '\n'; 
}

void SolveSaturation(){
	struct Rparams p = {ScalarMass, VectorMass,  IsoVectorMass, NMass, ElectronMass, ScalarField, VectorField, IsoVectorField, B, C, BaryonDensity,
	       	AsymmetryEnergy, Compressibility, BindingPerNucleon0, MstarPerMn0, AlphaAsymmetry, LambdaMass, HadronCouplingOmega, HadronCouplingSigma, 
		HadronCouplingRho, MuonMass, XiMinusMass, SigmaMinusMass, SigmaZeroMass, SigmaPlusMass, XiZeroMass};
	const gsl_multiroot_fsolver_type *T;	//DECLARE VARIABLE FOR TYPE OF SOLVER
	gsl_multiroot_fsolver *S;		//DECLARE NAME FOR SOLVER WORKSPACE

	int status;				//USED TO UPDATE ON STATUS (IS ACTUALLY BOOLEAN)
	size_t i, iter =0;

	const size_t n = 9;			//DIMENSIONALITY OF PROBLEM
				//  Msig,  Mw,   Mr,  Mn,   Me,  gs, gw, gr,   b,    c,    I3,  Jn,  Jp,  Je,  ban, bap, Qe,  Qp, nB             Easym, Kcompress   Eo 

	gsl_multiroot_function func = {&GM1_f, n, &p};	//MAKE A FUNCTION FOR SOLVER OUT OF ROSENBROCK, OF DIMENSION n, WITH PARAMS p
	double x_initial[9] = {1, 1, 1, 1, 1, 1, 1e-5, -1e-5, 0.1};	//MAKE INITIAL GUESS FOR VARIABLES TO BE SOLVED FOR
	gsl_vector *x = gsl_vector_alloc(n);	//MAKE A VECTOR FOR SOLVER TO STORE VARIABLES
	
	gsl_vector_set(x, 0, x_initial[0]);
	gsl_vector_set(x, 1, x_initial[1]);	//SET VARIABLES TO BE SOLVED FOR TO INITIAL VALUES
	gsl_vector_set(x, 2, x_initial[2]);
	gsl_vector_set(x, 3, x_initial[3]);
	gsl_vector_set(x, 4, x_initial[4]);
	gsl_vector_set(x, 5, x_initial[5]);
	gsl_vector_set(x, 6, x_initial[6]);
	gsl_vector_set(x, 7, x_initial[7]);
	gsl_vector_set(x, 8, x_initial[8]);
	
	T = gsl_multiroot_fsolver_hybrids;	//MAKE THE TYPE OF SOLVER A HYBRID S SOLVER
	S = gsl_multiroot_fsolver_alloc(T, n);	//MAKE THE WORKSPACE FOR SOLVER TYPE T OF DIMENSION 2
	gsl_multiroot_fsolver_set(S, &func, x);	//USE THE SOLVER TO SOLVE FUNCTIONS func AND FOR VARIABLES x (both have to be gsl_multiroots
	
		//print_state(iter, S);			//PRINT THE ORIGINAL STATE
	while(status = GSL_CONTINUE && iter < 1000){
		iter ++;
		status = gsl_multiroot_fsolver_iterate (S);	//ITERATE THE WORKSPACE UNTIL THE THINGS IS SOLVED

		if(status) break;				//CHECK IF IT'S WORKING. IF NOT, STOP
		status = gsl_multiroot_test_residual(S->f, 1e-7);	//SET TOLERANCE FOR CONVERGENCE
	}

/*========================== PULL CONSTANTS FROM SOLUTION =================================*/
	for(int vari = 0; vari < n; vari++){
		SaturationVariables[vari] = gsl_vector_get(S->x, vari);
	}
	print_state(iter, S);				//PRINT THE CURRENT STATUS
	cout << "status1 = " << gsl_strerror(status) << '\n';
	gsl_multiroot_fsolver_free(S);			//CLEAR MEMORY
	gsl_vector_free(x);				//ON EVERYTHING

}
