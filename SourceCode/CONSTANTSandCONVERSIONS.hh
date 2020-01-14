using namespace std;

/*============== UNIT CONVERSION natural units where c = hBAR = kB  = 1 ===========*/
double pi = M_PI;
const double convGeVtokg = 1.78e-27;
const double convGeVtoinvcm = 5.06e13;
const double convGeVtoinvs = 1.52e24;
const double convGeVtoK =1.16e13;

const double convcmtoinvGeV = 5.06e13;
const double convstoinvGeV = 1.52e24;

const double convinvGeVtocm = 1.98e-14;
const double convinvGeVtos = 6.58e-25;

const double convkgtoGeV = 5.62e26;
const double convinvcmtoGeV = 1.98e-14;
const double convinvstoGeV = 6.58e-25;
const double convKtoGeV =8.62e-14;
const double convGeVtoinvFM = 0.0008065*(2*pi)*1e-3;			

/*===================== CONSTANTS ==========================*/
const double GNewton = 6.67408e-11;					// m^3/(kg*s^2)
const double kBoltzmann = 1.3806e-23;					// m^2*kg/(s^2*K)
const double clight = 3e8;						// m/s
const double Planckbar = (6.626e-34)/(2*pi);				// m^2*kg/s

const double GNewton0 = (6.67408e-11)*1e6*pow(convcmtoinvGeV,3.)/(convkgtoGeV*pow(convstoinvGeV,2.));			// natural, 1/M_planck^2
const double kBoltzmann0 = 1.3806e-23*1e4*pow(convcmtoinvGeV,2.)*convkgtoGeV/(pow(convstoinvGeV,2.)*convKtoGeV);	// natural, kB = 1
const double clight0 = 3e8*1e2*convcmtoinvGeV/convstoinvGeV;								// natural, c = 1
const double Planckbar0 = Planckbar*1e4*pow(convcmtoinvGeV,2)*convkgtoGeV/convstoinvGeV;				// natural, hBar = 1
const double MeVtoinvFM = 0.0008065*(2*pi);			// approximate conversion factor

