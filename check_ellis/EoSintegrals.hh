#include "EoSfunctions.hh"

double Int_eye1(double kmin, double kmax, double Mst, int calls){
	gsl_function EYE1;

	EYE1.function = &I1;
	EYE1.params = &Mst;

	gsl_integration_workspace * WEYE = gsl_integration_workspace_alloc(calls);
	double relerr0 = 1e-7;

	double Ione = 0;
	double Ioneerr = 0;

	gsl_integration_qags(&EYE1, kmin, kmax, 0., relerr0, calls, WEYE, &Ione, &Ioneerr);
	gsl_integration_workspace_free(WEYE);
	
	return Ione;
}

double Int_eye2(double kmin, double kmax, double Mst, int calls){
	gsl_function EYE2;

	EYE2.function = &I2;
	EYE2.params = &Mst;

	gsl_integration_workspace * WEYE = gsl_integration_workspace_alloc(calls);
	double relerr0 = 1e-7;

	double Itwo = 0;
	double Itwoerr = 0;

	gsl_integration_qags(&EYE2, kmin, kmax, 0., relerr0, calls, WEYE, &Itwo, &Itwoerr);
	gsl_integration_workspace_free(WEYE);
	
	return Itwo;
}

double Int_eye3(double kmin, double kmax, double Mst, int calls){
	gsl_function EYE3;

	EYE3.function = &I3;
	EYE3.params = &Mst;

	gsl_integration_workspace * WEYE = gsl_integration_workspace_alloc(calls);
	double relerr0 = 1e-7;

	double Ithree = 0;
	double Ithreeerr = 0;

	gsl_integration_qags(&EYE3, kmin, kmax, 0., relerr0, calls, WEYE, &Ithree, &Ithreeerr);
	gsl_integration_workspace_free(WEYE);
	
	return Ithree;
}

double Int_eyeP(double kmin, double kmax, double Mst, int calls){
	gsl_function EYEP;

	EYEP.function = &IP;
	EYEP.params = &Mst;

	gsl_integration_workspace * WEYE = gsl_integration_workspace_alloc(calls);
	double relerr0 = 1e-7;

	double Ipee = 0;
	double Ipeeerr = 0;

	gsl_integration_qags(&EYEP, kmin, kmax, 0., relerr0, calls, WEYE, &Ipee, &Ipeeerr);
	gsl_integration_workspace_free(WEYE);
	
	return Ipee;
}

/* =================== ANALYTIC APPROXIMATIONS FOR NUMERICAL INTEGRALS =========================*/
double eye3(double kF, double m){
	double x = kF/m;
	double t = sqrt(1 + x*x);
	return m*m*(0.5*x*t + x/t - (1.5)*log(x+t));
}
	 
double eye1(double kF, double m){
	double x = kF/m;
	double t = sqrt(1 + x*x);
	return 0.5*pow(m, 3)*(x*t - log(x+t));
}
double eye2(double kF, double m){
	double x = kF/m;
	double t = sqrt(1 + x*x);
//        cout << "KF: " << x << " " << t << '\n';
	return 1./8* pow(m, 4)*(x*t*(1.0 + 2.0*pow(x, 2.0)) - log(x+t));
}

double eyeP(double kF, double m){
	double x = kF/m;
	double t = sqrt(kF*kF + m*m);
	return 1./8* (kF *t*(2*kF*kF - 3*m*m) + 3*pow(m,4.)*log(t + kF));
}
