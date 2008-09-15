/**
 * This package contains functions to calculate enthalpies, specific heats,
 * densities and partial pressures of components of interest for our
 * modelling activities.
 * Our components are oxygen, nitrogen, carbon dioxide, water and methanol, of
 * which only the latter two can occur in liquid form.
 */
encapsulated package Thermo
/**
 * Importing some useful definitions from the Modelica standard library.
 */
import Modelica.SIunits.Temperature;
import Modelica.SIunits.Pressure;
import Modelica.SIunits.MolarInternalEnergy;
import Modelica.SIunits.MolarHeatCapacity;
import Modelica.SIunits.MolarEntropy;
import Modelica.SIunits.PartialPressure;
import Modelica.SIunits.Density;
import Modelica.SIunits.MolarMass;
import Modelica.SIunits.MoleFraction;
// Definition lacking from Modelica library.
public type MolarEnthalpy = MolarInternalEnergy;

protected

/**
 * Molar masses of components.
 */
constant MolarMass mw_o2 = 31.9988e-3;
constant MolarMass mw_n2 = 28.01348e-3;
constant MolarMass mw_co2 = 44.0095e-3;
constant MolarMass mw_h2o = 18.0153e-3;
constant MolarMass mw_ch3oh = 32.0419e-3;


/**
* An abstract function containing the parameters for the Shomate equation to
* obtain the properties of oxygen.
* Source: http://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447&Units=SI&Mask=1#Thermo-Gas
*/
partial function OxygenShomateParameters
	input Temperature T(min=298.0, max=6000.0);
	protected
		constant Real A = 29.65900;
		constant Real B = 6.137261;
		constant Real C = -1.186521;
		constant Real D = 0.095780;
		constant Real E = -0.219663;
		constant Real F = -9.861391;
		constant Real G = 237.9480;
		constant Real H = 0.000000;
end OxygenShomateParameters;

/**
* An abstract function containing the parameters for the Shomate equation to
* obtain the properties of nitrogen.
* Source: http://webbook.nist.gov/cgi/cbook.cgi?ID=C7727379&Units=SI&Mask=1#Thermo-Gas
*/
partial function NitrogenShomateParameters
	input Temperature T(min=298.0, max=6000.0);
	protected
		constant Real A = 26.09200;
		constant Real B = 8.218801;
		constant Real C = -1.976141;
		constant Real D = 0.159274;
		constant Real E = 0.044434;
		constant Real F = -7.989230;
		constant Real G = 221.0200;
		constant Real H = 0.000000;
end NitrogenShomateParameters;

/**
* An abstract function containing the parameters for the Shomate equation to
* obtain the properties of carbon dioxide.
* Source: http://webbook.nist.gov/cgi/cbook.cgi?ID=C124389&Units=SI&Mask=1#Thermo-Gas
*/
partial function CarbonDioxideShomateParameters
	input Temperature T(min=298.0, max=1200.0);
	protected
		constant Real A = 24.99735;
		constant Real B = 55.18696;
		constant Real C = -33.69137;
		constant Real D = 7.948387;
		constant Real E = -0.136638;
		constant Real F = -403.6075;
		constant Real G = 228.2431;
		constant Real H = -393.5224;
end CarbonDioxideShomateParameters;

/**
* An abstract function containing the parameters for the Shomate equation to
* obtain the properties of liquid water.
* Source: http://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Units=SI&Mask=2#Thermo-Condensed
*/
partial function LiquidWaterShomateParameters
	input Temperature T(min=298.0, max=500.0);
	protected
		constant Real A = -203.6060;
		constant Real B = 1523.290;
		constant Real C = -3196.413;
		constant Real D = 2474.455;
		constant Real E = 3.855326;
		constant Real F = -256.5478;
		constant Real G = -488.7163;
		constant Real H = -285.8304;
end LiquidWaterShomateParameters;

/**
* This abstract function uses the Shomate parameters to calculate the enthalpy.
* Note that the Shomate parameters are still unspecified.
*/
partial function ShomateEnthalpy
	output MolarEnthalpy h;
	protected
		Real t;
	algorithm
		t := T / 1000.0;
		h := ( A*t + (B*t^2)/2 + (C*t^3)/3 + (D*t^4)/4 - E/t + F - H ) * 1000;
end ShomateEnthalpy;

/**
* This abstract function uses the Shomate parameters to calculate the specific
* heat. Note that the Shomate parameters are still unspecified.
*/
partial function ShomateSpecificHeat
	output MolarHeatCapacity cp;
	protected
		Real t;
	algorithm
		t := T / 1000.0;
		cp := A + B*t + C*t^2 + D*t^3 + E/t^2;
end ShomateSpecificHeat;

/**
* This abstract function uses the Shomate parameters to calculate the standard
* entropy. Note that the Shomate parameters are still unspecified.
*/
partial function ShomateEntropy
	output MolarEntropy s;
	protected
		Real t;
	algorithm
		t := T / 1000.0;
		s := A*log(t) + B*t + (C*t^2)/2 + (D*t^3)/3 - E/(2*t^2) + G;
end ShomateEntropy;


// The enthalpy of oxygen.
function h_o2
	extends ShomateEnthalpy;
	extends OxygenShomateParameters;
end h_o2;

// The specific heat of oxygen.
function cp_o2
	extends ShomateSpecificHeat;
	extends OxygenShomateParameters;
end cp_o2;


// The enthalpy of nitrogen.
function h_n2
	extends ShomateEnthalpy;
	extends NitrogenShomateParameters;
end h_n2;

// The standard enthalpy of formation of carbon dioxide.
constant MolarEnthalpy dhf_co2 = -393510.0;

// The specific heat of nitrogen.
function cp_n2
	extends ShomateSpecificHeat;
	extends NitrogenShomateParameters;
end cp_n2;


// The enthalpy of carbon dioxide.
function h_co2
	extends ShomateEnthalpy;
	extends CarbonDioxideShomateParameters;
end h_co2;

// The specific heat of carbon dioxide.
function cp_co2
	extends ShomateSpecificHeat;
	extends CarbonDioxideShomateParameters;
end cp_co2;

// The standard enthalpy of formation of liquid water.
constant MolarEnthalpy dhf_h2o_liq = -285830.0;

// The enthalpy of liquid water.
function h_h2o_liq
	extends ShomateEnthalpy;
	extends LiquidWaterShomateParameters;
end h_h2o_liq;

// The specific heat of liquid water.
function cp_h2o_liq
	extends ShomateSpecificHeat;
	extends LiquidWaterShomateParameters;
end cp_h2o_liq;

// The standard enthalpy of formation of water vapour.
constant MolarEnthalpy dhf_h2o_gas = -241826.0;

/**
* These parameters are necessary to calculate the enthalpy or specific heat of
* water vapour. Results are reliable between 300 and 2500 kelvin.
* Source: Perry's Chemical Engineers' Handbook, 7th edition, page 2-163.
*
* Note: parameter D has been added to those from Perry, so that enthalpy would
* be 0 at 25°C.
*/
partial function WaterVapourParameters
	input Temperature T(min=300.0, max=2500.0);
	protected
		constant Real A = 8.22;
		constant Real B = 0.00015;
		constant Real C = 0.00000134;
		constant Real D = 10326.2823;
end WaterVapourParameters;

// The enthalpy of gaseous water.
function h_h2o_gas
	extends WaterVapourParameters;
	output MolarEnthalpy h;
	algorithm
		// NOTE conversion from calories.
		h := ( A*T + (B*T^2)/2 + (C*T^3)/3 ) * 4.184 - D;
end h_h2o_gas;

// The specific heat of gaseous water.
function cp_h2o_gas
	extends WaterVapourParameters;
	output MolarHeatCapacity cp;
	algorithm
		// NOTE conversion from calories.
		cp := ( A + B*T + C*T^2 ) * 4.184;
end cp_h2o_gas;


// The standard enthalpy of formation of liquid methanol.
constant MolarEnthalpy dhf_ch3oh_liq = -238400.0;

/**
* These parameters are necessary to calculate the enthalpy or specific heat of
* liquid methanol. Results are reliable between 175.47 and 400 kelvin.
* Source: Perry's Chemical Engineers' Handbook, 7th edition, page 2-171.
*
* Note: parameter D has been added to those from Perry, so that enthalpy would
* be 0 at 25°C.
*/
partial function LiquidMethanolParameters
	input Temperature T(min=175.47, max=400.0);
	protected
		constant Real A = 1.058E5;
		constant Real B = -362.23;
		constant Real C = 0.9379;
		constant Real D = 23730.2384;
end LiquidMethanolParameters;

// The enthalpy of liquid methanol.
function h_ch3oh_liq
	extends LiquidMethanolParameters;
	output MolarEnthalpy h;
	algorithm
		h := ( A*T + (B*T^2)/2 + (C*T^3)/3 ) / 1000 - D;
end h_ch3oh_liq;

// The specific heat of liquid methanol.
function cp_ch3oh_liq
	extends LiquidMethanolParameters;
	output MolarHeatCapacity cp;
	algorithm
		cp := ( A + B*T + C*T^2 ) / 1000;
end cp_ch3oh_liq;

// The standard enthalpy of formation of gaseous methanol.
constant MolarEnthalpy dhf_ch3oh_gas = -205000.0;

/**
* These parameters are necessary to calculate the enthalpy or specific heat of
* gaseous methanol. Results are reliable between 200 and 1500 kelvin.
* Source: Perry's Chemical Engineers' Handbook, 7th edition, page 2-179.
*/
partial function GaseousMethanolParameters
	input Temperature T(min=200.0, max=1500.0);
	protected
		constant Real C1 = 0.3925E5;
		constant Real C2 = 0.879E5;
		constant Real C3 = 1.9165E3;
		constant Real C4 = 0.5365E5;
		constant Real C5 = 896.7;
end GaseousMethanolParameters;

/**
* This function returns the specific heat of gaseous methanol at a given
* temperature.
* The formula used in this function is the one reported in Perry's Chemical
* Engineers' Handbook, 7th edition, page 2-182.
*/
function cp_ch3oh_gas
	extends GaseousMethanolParameters;
	output MolarHeatCapacity cp;
	algorithm
		cp :=( C1 + C2*( C3/T/( sinh(C3/T) ) )^2
		          + C4*( C5/T/cosh(C5/T) )^2 ) / 1000;
end cp_ch3oh_gas;

/**
* This function returns the enthalpy of gaseous methanol at a given
* temperature.
* The general integral used in this function is the analytic integral of the
* formula reported in Perry's Chemical Engineers' Handbook, 7th edition, page
* 2-182, and has been obtained with Maple.
*/
function h_ch3oh_gas
	input Temperature T(min=200.0, max=1500.0);
	output MolarEnthalpy h;
	protected
		/**
		* This is the general integral of the function for specific heat.
		*/
		function generalIntegral
			extends GaseousMethanolParameters;
			output MolarEnthalpy y;
			algorithm
				y := ( C1*T + C2*C3*cosh(C3/T)/sinh(C3/T)
				            - C4*C5*sinh(C5/T)/cosh(C5/T) ) / 1000;
		end generalIntegral;
	algorithm
		h := generalIntegral( T ) - generalIntegral( 298.15 );
end h_ch3oh_gas;


/**
* This abstract function implements Antoine's law for vapour pressure.
*/
partial function AntoineLaw
	input Temperature T;
	output PartialPressure p;
	protected
		parameter Real A;
		parameter Real B;
		parameter Real C;
	algorithm
		// This is the Antoine Law itself. Note the conversion bar->Pa.
		p := 10.0^( A - B/(T+C) ) * 1e5;
end AntoineLaw;

/**
* This function returns the vapour pressure of water given the temperature.
* Source: http://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Mask=4#Thermo-Phase
* (data taken from Stull, 1947).
*/
function p_h2o
	extends AntoineLaw( T.min=255.8, T.max=373.0,
	                    A = 4.65430, B = 1435.264, C = -64.848 );
end p_h2o;

/**
* This function returns the vapour pressure of methanol given the temperature.
* Source: http://webbook.nist.gov/cgi/cbook.cgi?ID=C67561&Mask=4#Thermo-Phase
* (data taken from Ambrose and Sprake, 1970).
*/
function p_ch3oh
	extends AntoineLaw( T.min=288.0, T.max=356.83,
	                    A = 5.20409, B = 1581.341, C = -33.50 );
end p_ch3oh;

/**
 * This function returns the bubble pressure of a water-methanol solution at a
 * given temperature.
 */
function bubblePressure
	input MoleFraction x_ch3oh;
	input Temperature T;
	output Pressure p;
	protected
		MoleFraction x_h2o = 1.0 - x_ch3oh;
	algorithm
		p := p_h2o(T) * x_h2o + p_ch3oh(T) * x_ch3oh;
end bubblePressure;

/**
 * This function returns the dew pressure of a water-methanol mixture at a
 * given temperature. Note that there could be other non-condensing species
 * such as nitrogen in the gaseous mixture, so the the water and methanol
 * fractions do not necessarily sum up to one.
 */
function dewPressure
	input MoleFraction y_ch3oh;
	input MoleFraction y_h2o;
	input Temperature T;
	output Pressure p;
	algorithm
		p := 1.0 / ( y_h2o/p_h2o(T) + y_ch3oh/p_ch3oh(T) );
end dewPressure;

/**
 * TODO
 */
function LinearInterpolation
	input Real x;
	input Real[:] dataX;
	input Real[:] dataY;
	output Real y;
	protected
		Integer pos;
		Real f;
	algorithm
		if x < dataX[1] then
			y := dataY[1];
			return;
		elseif x >= dataX[end] then
			y := dataY[end];
			return;
		else
			pos := 1;
			while x > dataX[pos+1] loop
				pos := pos + 1;
			end while;
			f := (x - dataX[pos]) / (dataX[pos+1] - dataX[pos]);
			y := dataY[pos] + f * ( dataY[pos+1] - dataY[pos] );
		end if;

end LinearInterpolation;

/**
 * This function returns the density of water in liquid phase at one standard
 * atmosphere of pressure.
 */
function rho_h2o
	input Temperature T;
	output Density rho;
	external "C" annotation(Include="#include \"/data/portsys/projects/modelling/Modelica/rho_h2o.c\"");
end rho_h2o;

/**
 * This function returns the density of methanol in liquid phase at one standard
 * atmosphere of pressure.
 */
function rho_ch3oh
	input Temperature T;
	output Density rho;
	external "C" annotation(Include="#include \"/data/portsys/projects/modelling/Modelica/rho_ch3oh.c\"");
end rho_ch3oh;

/**
 * Dispatcher functions: depending on inputs, they call the previous functions.
 * First index for components:
 *   1) Methanol
 *   2) Water
 *   3) Oxygen
 *   4) Carbon dioxide
 *   5) Nitrogen
 * Second index for phase:
 *   1) Gaseous
 *   2) Liquid
 * FIXME it should be done with enumerations, but they do not work yet.
 * When they do, change the code.
 */

public

/*
type C     = enumeration ( ch3oh, h2o, o2, co2 );
type Phase = enumeration ( gas, liq );
*/

function mw
	input Integer n "Component";
	output MolarMass m;
	algorithm
		if n == 1 then
			m := mw_ch3oh;
		elseif n == 2 then
			m := mw_h2o;
		elseif n == 3 then
			m := mw_o2;
		elseif n == 4 then
			m := mw_co2;
		elseif n == 5 then
			m := mw_n2;
		else
			m := 0.0;	// FIXME assert and terminate don't work.
		end if;
end mw;

function dhf
	input Integer n "Component";
	input Integer p = 1 "Phase, gaseous by default";
	output MolarEnthalpy f;
	algorithm
		if n == 1 and p == 1 then
			f := dhf_ch3oh_gas;
		elseif n == 1 and p == 2 then
			f := dhf_ch3oh_liq;
		elseif n == 2 and p == 1 then
			f := dhf_h2o_gas;
		elseif n == 2 and p == 2 then
			f := dhf_h2o_liq;
		elseif n == 3 then
			f := dhf_o2;
		elseif n == 4 then
			f := dhf_co2;
		elseif n == 5 then
			f := dhf_n2;
		else
			f := 0.0;	// FIXME assert and terminate don't work.
		end if;
end dhf;

function h
	input Temperature T;
	input Integer n "Component";
	input Integer p = 1 "Phase, gaseous by default";
	output MolarEnthalpy H;
	algorithm
		if n == 1 and p == 1 then
			H := h_ch3oh_gas(T);
		elseif n == 1 and p == 2 then
			H := h_ch3oh_liq(T);
		elseif n == 2 and p == 1 then
			H := h_h2o_gas(T);
		elseif n == 2 and p == 2 then
			H := h_h2o_liq(T);
		elseif n == 3 and p == 1 then
			H := h_o2(T);
		elseif n == 4 and p == 1 then
			H := h_co2(T);
		elseif n == 5 and p == 1 then
			H := h_n2(T);
		else
			H := 1e20;	// FIXME assert and terminate don't work.
		end if;
end h;

function cp
	input Temperature T;
	input Integer n "Component";
	input Integer p = 1 "Phase, gaseous by default";
	output MolarHeatCapacity CP;
	algorithm
		if n == 1 and p == 1 then
			CP := cp_ch3ocp_gas(T);
		elseif n == 1 and p == 2 then
			CP := cp_ch3ocp_liq(T);
		elseif n == 2 and p == 1 then
			CP := cp_h2o_gas(T);
		elseif n == 2 and p == 2 then
			CP := cp_h2o_liq(T);
		elseif n == 3 and p == 1 then
			CP := cp_o2(T);
		elseif n == 4 and p == 1 then
			CP := cp_co2(T);
		elseif n == 5 and p == 1 then
			CP := cp_n2(T);
		else
			CP := 1e20;	// FIXME assert and terminate don't work.
		end if;
end cp;

function p_vap
	input Temperature T;
	input Integer n "Component";
	output PartialPressure p;
	algorithm
		if n == 1 then
			p := p_ch3oh(T);
		elseif n == 2 then
			p := p_h2o(T);
		else
			p := 1e20;	// FIXME assert and terminate don't work.
		end if;
end p_vap;

function rho
	input Temperature T;
	input Integer n "Component";
	input Integer p = 1 "Phase, gaseous by default";
	output Density RHO;
	import Modelica.Constants.R;
	protected constant Pressure p_env = 101325.0;
	algorithm
		if n == 1 and p == 2 then
			RHO := rho_ch3oh(T);
		elseif n == 2 and p == 2 then
			RHO := rho_h2o(T);
		elseif n >= 1 and n <= 5 and p == 1 then
			RHO := mw(n) * p_env / R / T;	// NOTE Assuming ideal gas.
		else
			RHO := 0.0;	// FIXME assert and terminate don't work.
		end if;
end rho;

end Thermo;
