#include "linear_interpolation.h"

/**
 * Liquid water density as a function of temperature, data taken from NIST:
 * http://webbook.nist.gov/cgi/fluid.cgi?P=1&TLow=273.15&THigh=373.15&TInc=5&Applet=on&ID=C67561&Action=Load&Type=IsoBar&TUnit=K&PUnit=atm&DUnit=kg%2Fm3&HUnit=kJ%2Fmol&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm&RefState=DEF
 */
double rho_ch3oh( const double T )
{
	int n = 14;
	double data_T[] = { 273.15, 278.15, 283.15, 288.15, 293.15, 298.15, 303.15,
	                    308.15, 313.15, 318.15, 323.15, 328.15, 333.15,
	                    337.63 };
	double data_D[] = { 809.73, 805.05, 800.37, 795.69, 791.01, 786.33, 781.63,
	                    776.91, 772.17, 767.39, 762.58, 757.72, 752.81,
	                    748.36 };
	return linear_interpolation( T, data_T, data_D, n );
}
