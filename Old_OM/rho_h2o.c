#include "linear_interpolation.h"

/**
 * Liquid water density as a function of temperature, data taken from NIST:
 * http://webbook.nist.gov/cgi/fluid.cgi?P=1&TLow=273.15&THigh=373.15&TInc=5&Applet=on&ID=C7732185&Action=Load&Type=IsoBar&TUnit=K&PUnit=atm&DUnit=kg%2Fm3&HUnit=kJ%2Fmol&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm&RefState=DEF
 */
double rho_h2o( const double T )
{
	int n = 21;
	double data_T[] = { 273.16, 278.16, 283.16, 288.16, 293.16, 298.16, 303.16,
	                    308.16, 313.16, 318.16, 323.16, 328.16, 333.16, 338.16,
	                    343.16, 348.16, 353.16, 358.16, 363.16, 368.16,
	                    373.12 };
	double data_D[] = { 999.84, 999.97, 999.70, 999.10, 998.21, 997.05, 995.65,
	                    994.03, 992.21, 990.21, 988.03, 985.69, 983.19, 980.55,
	                    977.76, 974.84, 971.78, 968.60, 965.30, 961.88,
	                    958.37 };
	return linear_interpolation( T, data_T, data_D, n );
}
