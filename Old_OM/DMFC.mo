type MolarFlowRate = Real(final quantity="MolarFlowRate", final unit="mol/s");
type HeatTransferCoefficient = Real( final quantity="HeatTransferCoefficient",
                                     final unit="W/K" );


/**
 * A connector for the various units; ensures continuity of enthalpy and
 * molar flows. The composition array is arranged so:
 * 1) Methanol
 * 2) Water
 * 3) Oxygen
 * 4) Carbon dioxide
 * 5) Nitrogen
 */
connector Flow

	import Modelica.SIunits.EnthalpyFlowRate;
	import Modelica.SIunits.MoleFraction;
	/* H is the amount of enthalpy necessary to bring the components from
	 * standard conditions (298.15, liquid water and methanol) to their actual
	 * one at constant, environmental pressure. */
	flow EnthalpyFlowRate H;
	// This is the total molar flow rate.
	flow MolarFlowRate F;
	// This is the overall mole fraction.
	MoleFraction[5] z;

end Flow;


/**
 * This class defines the general properties of all stirred tanks, that is a
 * single representative temperature, an array of amounts of substance with five
 * elements (for methanol, water, oxygen, carbon dioxide and nitrogen), and the
 * flows with which the tank is connected to other ones.
 */
model StirredTank

	import Modelica.SIunits.AmountOfSubstance;
	import Modelica.SIunits.Pressure;
	import Modelica.SIunits.Temperature;
	import Modelica.SIunits.Heat;
	import Modelica.SIunits.HeatCapacity;
	import Modelica.SIunits.HeatFlowRate;
	import Modelica.SIunits.MoleFraction;
	import Modelica.SIunits.Volume;
	import Modelica.Constants.R;

	import Thermo.p_vap;
	import Thermo.rho;
	import Thermo.cp;
	import Thermo.dhf;

	outer Temperature T_env "Environment temperature.";
	outer Pressure p_env "Environment pressure.";

	parameter String name = "Nameless" "Tank identifier.";
	parameter HeatCapacity Cp "Heat capacity of tank (glass, lid, ...).";
	parameter HeatTransferCoefficient h_k = 0.0 "Default: perfect insulation.";
	parameter Integer m "Number of flows.";
	parameter Volume V "Total volume of the tank.";

	Temperature T(start=T_env, fixed=true) "Representative tank temperature.";
	Heat Q "Heat exchanged with environment (usually negative).";
	Flow[m] flows "Connections with other objects.";

	/* This declaration allows the initialisation algorithm to adjust the
	 * composition as it pleases. As we need only one component to be set at
	 * initialisation time, child classes may decide to redeclare this
	 * variable's "fixed" attribute to enforce particular initialisation
	 * policies, e.g. "fill with air", "fill with water" and so on. */
	AmountOfSubstance[5] n "Moles of each species.";
	AmountOfSubstance n_tot "Total moles.";
	MoleFraction[5] z "Overall mole fractions";

	Volume V_l "Amount of liquid volume.";
	AmountOfSubstance n_l "Total moles of liquid.";
	MoleFraction[5] x "Molar fraction in liquid phase.";

	Volume V_g "Amount of gaseous volume.";
	AmountOfSubstance n_g "Total moles of gas.";
	MoleFraction[5] y "Molar fraction in gas phase.";

	protected
		HeatFlowRate enthalpyInflow, sensibleHeat, latentHeat, dmHeat;

	equation
		/* NOTE: In this class I do not know yet where the fittings are, where
		 * they are connected and whether they draw from a specific phase; the
		 * relation between their and the tank's composition *cannot* be put in
		 * this class. */

		// Convenience variable for the total number of moles.
		n_tot = sum(n);
		// Overall molar fraction.
		z = n / n_tot;

		// Mass balance, piece of cake.
		der(n) = {sum( flows[j].F * flows[j].z[i] for j in 1:m ) for i in 1:5};

		/* The next section deals with the energy balance. It is a bit intricate
		 * because we have an unknown number of inflows, plus the tank has in
		 * general two phases. */
		// Complete energy balance
		Q + enthalpyInflow = Cp * der( T ) + sensibleHeat + latentHeat ;
		// Exchange with environment.
		Q = h_k * ( T_env - T );
		// The sum of enthalpy inflows from all flows.
		enthalpyInflow = sum( flows[j].H for j in 1:m );
		// The specific heat (J/K) of the components in the tank, liquid and gas
		sensibleHeat = der( T ) * ( sum( cp(T,i,1) * y[i] * n_g for i in 1:5 )
		               + sum( cp(T,i,2) * x[i] * n_l for i in 1:2 ) )
		               + dmHeat;
		/* The heat spent on warming incoming matter to tank temperature T.
		 * Note that methanol and water are a bit more complex as they change
		 * phase according to temperature. */
		dmHeat = sum( ( h(T,i,1)-h(298.15,i,1) ) * der( n[i] ) for i in 3:5 ) +
		         sum( ( h(T,i,1)-h(298.15,i,1) ) * der( y[i]*n_g ) +
		              ( h(T,i,2)-h(298.15,i,2) ) * der( x[i]*n_l )
		              for i in 1:2 );
		// The heat to bring water and methanol into gas phase
		latentHeat = sum( ( dhf(i,1) - dhf(i,2) ) * der( y[i]*n_g )
		                  for i in 1:2 );


		/* The next equations deal with the compositions of the gaseous and
		 * liquid phases, and the liquid / gaseous volumes. */
		// Consistency: sum of gas and liquid volumes is constant.
		V = V_l + V_g;
		// The relationship between liquid volume, moles and compositions.
		// NOTE: assumed no mixture volume.
		V_l = n_l * sum( x[i]*mw(i)/rho(T,i,2) );
		// The relationship between gaseous volume, moles and compositions.
		// NOTE: assumed no mixture volume.
		V_g = n_g * sum( y[i]*mw(i)/rho(T,i,1) );
		// Molar fractions sum to 1.
		sum( x ) = 1.0;
		sum( y ) = 1.0;
		// No gases in liquid phase.
		for i in 3:5 loop
			x[i] = 0.0;
		end for;
		// Raoult's law for methanol and water.
		for i in 1:2 loop
			p_env * y[i] = p_vap(T, i) * x[i];
		end for;
		// Relation between fractions in phases and total moles.
		for i in 1:5 loop
			n[i] = y[i] * n_g + x[i] * n_l;
		end for;


		/* Sanity checks: component cannot be less than zero, nor can molar
		 * fractions exit the [0, 1] interval. */
		for i in 1:5 loop
			assert( n[i] >= 0.0,
			        name + ": component "+i+"has a negative quantity." );
			assert( x[i] >= 0.0 and x[i] <= 1.0, name +
			        ": liquid molar fraction"+i+" is outside boundaries." );
			assert( y[i] >= 0.0 and y[i] <= 1.0, name +
			        ": gaseous molar fraction "+i+"is outside boundaries." );
		end for;
		// Termination criterion: liquid moles go below zero.
		assert( n_l >= 0.0, name +
		       ": liquid moles have become less than zero." );
		// Termination criterion: gaseous moles go below zero.
		assert( n_g >= 0.0, name + ": gas moles have become less than zero." );
		// Termination criterion: liquid phase goes below zero.
		assert( V_l >= 0.0, name +
		        ": liquid phase has become less than zero." );
		// Termination criterion: gaseous phase goes below zero.
		assert( V_g >= 0.0, name + ": gas phase has become less than zero." );
		/* This checks that all compositions in each mass flow sum to one, and
		 * that each element is between 0 and 1.
		 * Note that these checks cannot be put in the Flow class, because
		 * connector classes have no equations. */
		for i in 1:m loop
			assert( sum( { flows[i].z for i in 1:m } ) == 1.0, name +
			             ": sum of components in flow "+i+" is not unitary." );
			for j in 1:5 loop
				assert( flows[i].z[j] <= 1.0 and flows[i].z[j] >= 0.0, name +
				": component "+j+" of flow "+i+" is "+flows[i].z[j]+"." );
			end for;
		end for;

end StirredTank;


/**
 * This class extends the generic stirred tank in specifying a geometry for the
 * tank, a vertical cylinder with a given horizontal section A.
 */
model VerticalCylindricalTank

	extends StirredTank;
	import Modelica.SIunits.Area;
	import Modelica.SIunits.Length;

	parameter Area A;
	Length level;

	equation
		level = V_l / A;

end VerticalCylindricalTank;


/**
 * This class extends the generic stirred tank in specifying a geometry for the
 * tank, a horizontal cylinder with a given vertical section A.
 */
model HorizontalCylindricalTank

	extends StirredTank;
	import Modelica.SIunits.Area;
	import Modelica.SIunits.Length;
	import Modelica.Constants.pi;

	parameter Area A;
	Length level;
	Length L = 2.0*sqrt(A/pi) "The inner diameter, or also height of the tank";

	equation
		pi * V_l / V = 2.0 * (2.0*level/L-1.0)*sqrt(level/L-(level/L)^2) +
		               asin(2.0*level/L-1.0) - asin(-1.0);

end HorizontalCylindricalTank;


/**
 * A pipe segment has uniform properties and only two flows, one entering and
 * one leaving. The one leaving has the overall composition of the segment, the
 * one entering is left unspecified and depends on whatever the segment is
 * connected to.
 * Note that in special cases both flows may be positive or negative (e.g. a
 * gas content being chilled or heated with no forced flow). In case a flow is
 * exactly zero, a given value for the composition vector is given to make sure
 * that no redundant equations are generated.
 */
model PipeSegment

	extends StirredTank(final m=2);

	equation
		for i in 1:m loop
			if flows[i].F < 0.0 then	// If flow is leaving...
				flows[i].z = z;
			/* Special case: no flow at all. Set all compositions to same value.
			 * This way, equations will be generated from both ends of the
			 * connection, but will be the same. */
			elseif flows[i].F == 0.0 then
				for j in 1:5 loop
					flows[i].z[j] = 0.2;
				end for;
			end if;
		end for;

end PipeSegment;


/**
 * This class gathers a number n of pipe segments into one single pipe. The only
 * job done by the class is joining the segments by connecting their flows in
 * sequence.
 * The pipe's connectors are segments[1].flow[1] and segments[end].flow[2].
 */
model Pipe

	parameter Integer n "Number of pipe segments";
	/* The segments are replaceable in case one defined a child class of
	 * PipeSegment that one wanted to use. */
	replaceable PipeSegment[n] segments;

	equation
		for i in 1:(n-1) loop
			connect( segment[i].flows[2], segment[i+1].flows[1] );
		end for;

end Pipe;


/**
 * A flow controller (such as a pump or a gas MFC) sets a certain molar (or
 * mass) flow between its two Flow connectors.
 */
model FlowController

	import Modelica.SIunits.MassFlowRate;
	import Thermo.mw;

	MolarFlowRate F;
	MassFlowRate m;

	Flow[2] flows;

	equation
		m = sum( {flows[1].z[i] * F * mw(i) for i in 1:5} );
		flows[1].F = F;
		connect( flows[1], flows[2] );

end FlowController;

/**
 * A subclass of VerticalCylindricalTank describing our laboratory fuel tank,
 * with volume of 1 liter, inner diameter of 8 cm, and two connections, one to
 * process and one to atmosphere.
 */
model FuelTank

	import Modelica.Constants.pi;
	extends VerticalCylindricalTank( name="Fuel tank",
	                                 final m = 2, V = 1E-3, A = 16E-4*pi,
	                                 Cp = 0.0, h_k = 0.0 );

	equation
		// This flow connects with the atmosphere, dry air for starters.
		if flows[1].F < 0.0 then	// If flow is leaving...
			flows[1].z = y;
		else
			flows[1].z = z_air;
		end if;
		// This is the flow for liquid, usually an outlet.
		if flows[2].F < 0.0 then	// If flow is leaving...
			flows[2].z = x;
		/* Dummy equation for the case in which there is no flow. This equation
		 * must be the same in all models. */
		elseif flows[2].F == 0.0 then
			flows[2].z = { 0.2, 0.2, 0.2, 0.2, 0.2 };
		end if;

	initial equation
		n[1] / V_l = 1E-3; // 1 M liquid solution
		level = V / A / 2.0; // Half-full
		n[4] = 0.0;	// No CO2
		y[5] / 0.79 = y[3] / 0.21; // N2/O2 in air proportions

end FuelTank;

/**
 * This model is our trusty mixer. It has four flows in liquid phase and one
 * connecting the gas phase to the atmosphere. Its internal diameter is TODO
 */
model Mixer

	extends VerticalCylindricalTank( name="Mixer",
	                                 final m = 5, V = 5E-4, A = /* TODO*/pi,
	                                 Cp = 0.0, h_k = 0.0 );
	equation
		for i in 1:4 loop
			if flows[i].F < 0.0 then	// If flow is leaving...
				flows[i].z = x;
			elseif flows[i].F == 0.0 then
				flows[i].z = { 0.2, 0.2, 0.2, 0.2, 0.2 };
			end if;
		end for;

		if flows[5].F < 0.0 then	// If flow is leaving...
			flows[5].z = y;
		elseif flows[5].F == 0.0 then
			flows[5].z = { 0.2, 0.2, 0.2, 0.2, 0.2 };
		end if;

	initial equation
		n[1] = 0.0;	// No methanol in initial solution
		n[4] = 0.0;	// No CO2
		level = V / A / 2.0; // Half-full
		y[5] / 0.79 = y[3] / 0.21; // N2/O2 in air proportions

end Mixer;

/**
 * TODO
 */
model Separator

	extends HorizontalCylindricalTank( name="Separator",
	                                   final m = 3, V = 2.5E-4, A = /*TODO*/pi,
	                                   Cp = 0.0, h_k = 0.0 );

	equation
		/* The first flow is the feed, which enters at middle height. In case
		 * the flow were leaving instead of entering, the composition will be
		 * that of vapour if the level is low, or that of liquid if the level is
		 * high. */
		// Flow leaving, level more than half:
		if flows[1].F < 0.0 and level > 0.5 * L then
			flows[1].z = x;
		// Flow leaving, level less than half:
		elseif flows[1].F < 0.0 then
			flows[1].z = y;
		elseif flows[1].F == 0.0 then
			flows[1].z = { 0.2, 0.2, 0.2, 0.2, 0.2 };
		end if;

		/* The second flow is on the top of the tank. If the gas volume is
		 * larger than zero, it will output vapour composition, otherwise the
		 * liquid one.*/
		 // Flow leaving, gas is present:
		if flows[2].F < 0.0 and V_g > 0.0 then
			flows[2].z = y;
		// Flow leaving, gas is not present:
		elseif flows[2].F < 0.0 then
			flows[2].z = x;
		elseif flows[2].F == 0.0 then
			flows[2].z = { 0.2, 0.2, 0.2, 0.2, 0.2 };
		end if;

		/* The third flow is on the bottom of the tank. If the liquid volume is
		 * larger than zero, it will output liquid composition, otherwise the
		 * vapour one. */
		 // Flow leaving, gas is present:
		if flows[3].F < 0.0 and V_l > 0.0 then
			flows[3].z = x;
		// Flow leaving, gas is not present:
		elseif flows[3].F < 0.0 then
			flows[3].z = y;
		elseif flows[3].F == 0.0 then
			flows[3].z = { 0.2, 0.2, 0.2, 0.2, 0.2 };
		end if;

	initial equation
		n[1] = 0.0;	// No methanol in initial solution
		n[4] = 0.0;	// No CO2
		level = V / A / 2.0; // Half-full
		y[5] / 0.79 = y[3] / 0.21; // N2/O2 in air proportions

end Separator;


/////////////// TESTS ////////////////////

class StupidTank

	import Modelica.SIunits.AmountOfSubstance;
	import Modelica.SIunits.MoleFraction;

	Flow[2] flows;
	Real[5] n(start={1,0,0,1,0}) "Moles of each species.";
	Real n_tot "Total moles.";
	Real[5] z "Overall mole fractions";
	protected
		Real temp[5];

equation
		n_tot = sum(n);
		z = n / n_tot;
//		der(n) = flows[1].F * flows[1].z + flows[2].F * flows[2].z;
		for i in 1:5 loop
			temp[i] = sum( { (flows[j].F) * (flows[j].z[i]) for j in 1:2 } );
		end for;
		der(n) = temp;


		flows[1].z = z;
		flows[2].z = z;
		flows[1].H = 0;
		flows[2].H = 0;

end StupidTank;


class Test

	import Modelica.SIunits.Pressure;
	import Modelica.SIunits.Temperature;
	StupidTank tank;

	equation
		tank.flows[1].F = -2;
		tank.flows[2].F = 1;

end Test;
