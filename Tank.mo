

package Tank 
type MolarFlowRate = Real(final quantity="MolarFlowRate", final unit="mol/s") 
    annotation (Documentation(info="<html>
<p>Just a definition lacking from the standard library.</p>
</html>"));
  
type HeatTransferCoefficient = Real (final quantity="HeatTransferCoefficient",
                                     final unit = "W/(K.m2)") 
                                                         annotation (
      Documentation(info="<html>
<p>Just a definition lacking from the standard library.</p>
</html>"));
  
  connector CheckPoint "What passes through a control surface" 
    import Modelica.SIunits.MoleFraction;
    import Modelica.SIunits.EnthalpyFlowRate;
    import Thermo.MolarEnthalpy;
    import Thermo.AllSpecies;
    
    MolarFlowRate F;
    MoleFraction[size(AllSpecies,1)] z(each min=0, each max=1);
    EnthalpyFlowRate H;
    MoleFraction[size(AllSpecies,1)] z_local(each min=0, each max=1);
    MolarEnthalpy h_local;
    
    annotation (Documentation(info="<html>
<p>This is a connector for the various tank units; it ensures continuity of enthalpy and
molar flows. It consists of two flow variables, the <em>enthalpy flow</em> and the total
<em>molar flow</em>, and of the <em>composition array</em>. In addition, there are two
other variables, the composition and the molar enthalpy of the tank to which the CheckPoint
is connected.</p>
<p>The enthalpy flow is defined as the enthalpy necessary to bring the
components in the molar-flow array from standard conditions, which is defined as 298.15 K
with water and methanol in liquid phase and other components in gas phase, to their actual
conditions of temperature and phase, at environmental pressure.</p>
<p>The Checkpoint connectors should <em>not</em> be connected between tanks; a FlowConnector
should instead be between these, to whose sides the CheckPoints are connected. This is
necessary because at some point it must be decided what composition and enthalpy content
to use, the one up- or downstream of the connection, and this decision cannot be made in
either tank object (because it lacks the values for the other one).</p>
</html>"), Icon(Ellipse(extent=[-80,80; 80,-80], style(
            pattern=0,
            thickness=2,
            gradient=3,
            fillColor=1,
            rgbfillColor={255,0,0}))));
  end CheckPoint;
  
  model Plug "A class that blocks a flow connection" 
    import Thermo.AllSpecies 
    annotation (Icon(Polygon(points=[-60,40; -60,-40; 60,40; 60,-40; -60,40], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4,
            fillColor=7,
            rgbfillColor={255,255,255}))));
    CheckPoint c annotation (extent=[-72,-10; -52,10]);
  protected 
    constant Integer k = size(AllSpecies, 1);
  equation 
    c.F = 0;
    c.z = zeros(k);
    c.H = 0;
    
  end Plug;
  
  model FlowConnector "A class to determine flow direction" 
  public 
    CheckPoint port1 annotation (extent=[-58,-10; -38,10]);
    CheckPoint port2 annotation (extent=[18,-10; 38,10]);
  equation 
    // Connection part.
    port1.F + port2.F = 0;
    port1.H + port2.H = 0;
    port1.z = port2.z;
    
    // Determining which side is taken for the values of composition and molar enthalpy.
    if port1.F >= 0 then  // If flow is entering into the tank on side one, or is zero:
      port2.H = port2.F * port2.h_local;
      port2.z = port2.z_local;
    else // If flow is entering from the second connection
      port1.H = port1.F * port1.h_local;
      port1.z = port1.z_local;
    end if;
    
    annotation (Documentation(info="<html>
<p><tt>FlowConnector</tt>s are objects that operate similarly to Modelica connectors,
but contain additional equations. They have two <tt>CheckPoint</tt> connectors which
lead to two tanks. Once given the direction of the overall molar flow, its task is to
decide which composition <tt>z_local</tt> and specific molar enthalpy <tt>h_local</tt>
should be used to set the flow's composition <tt>z</tt> and enthalpic flow <tt>H</tt>.</p>
<p>Each subclass of <tt>StirredTank</tt> should set the values of <tt>z_local</tt> and
<tt>h_local</tt> for each <tt>CheckPoint</tt> according to which phase would be drawn
from the connection (gas, liquid, mixed or depending on some variables, such as
liquid level).</p>
</html>"), Icon(
        Polygon(points=[-40,0; -20,-40; -20,40; -40,0],   style(
            pattern=0,
            fillColor=0,
            rgbfillColor={0,0,0})),
        Polygon(points=[20,0; 0,-40; 0,40; 20,0],     style(
            pattern=0,
            fillColor=0,
            rgbfillColor={0,0,0}))),
      Diagram);
  end FlowConnector;
  
  annotation (uses(Modelica(version="2.2.1")), Documentation(info="<html>
<p>This package contains various models related to fluid flow in stirred tanks.</p>
</html>"));
  
  model MethanolSolution "A model to use as a generic source" 
    import Thermo.mw;
    import Thermo.rho;
    import Thermo.h;
    import Modelica.SIunits.Temperature;
    import Modelica.SIunits.Concentration;
    import Modelica.SIunits.MoleFraction;
    
    outer Temperature T_env "Environment temperature.";
    
    parameter Concentration C = 1000 
      "Concentration of methanol in water, mol/m�.";
    parameter Temperature T = T_env "Source temperature.";
    
    MoleFraction x_ch3oh "Molar fraction of methanol.";
    MoleFraction x_h2o "Molar fraction of water.";
    
    CheckPoint p "Connection point of the source" 
      annotation (extent=[-10,-10; 10,10]);
  equation 
    assert( C >= 0, "Negative concentration given in MethanolSolution object.");
    assert( C <= rho(T,1,2)/mw(1), "Methanol concentration over limit (" + String(mw(1)/rho(T,1,2)) + " mol/m�).");
    
    C = x_ch3oh / ( x_ch3oh*mw(1)/rho(T,1,2) + x_h2o*mw(2)/rho(T,2,2));
    x_ch3oh + x_h2o = 1.0;
    
    p.z_local = {x_ch3oh, x_h2o, 0, 0, 0};
    p.h_local = x_ch3oh*h(T,1,2) + x_h2o*h(T,2,2);
    
    annotation (Icon(Ellipse(extent=[-40,40; 40,-40], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4,
            fillColor=7,
            rgbfillColor={255,255,255}))));
  end MethanolSolution;
  
  model EnvironmentPort "A flow connection to environment conditions." 
    import Thermo.dhf;
    import Thermo.h;
    import Thermo.p_vap;
    import Thermo.MolarEnthalpy;
    import Modelica.SIunits.Temperature;
    import Modelica.SIunits.Pressure;
    import Modelica.SIunits.MoleFraction;
    
    outer Real RH_env "The relative humidity in the laboratory, in percent.";
    outer Temperature T_env "The temperature in the laboratory.";
    outer Pressure p_env "The atmospheric pressure.";
    
    MoleFraction z_h2o "Molar fraction of water in environment air";
    MoleFraction z_o2 "Molar fraction of oxygen in environment air";
    MoleFraction z_n2 "Molar fraction of nitrogen in environment air";
    
    MolarEnthalpy h_n2 "The molar enthalpy of nitrogen.";
    MolarEnthalpy h_o2 "The molar enthalpy of oxygen.";
    MolarEnthalpy h_h2o "The molar enthalpy of water vapour.";
    MolarEnthalpy h_air "The molar enthalpy of air in the current conditions.";
    
    CheckPoint c annotation (extent=[-100,-60; -80,-40]);
  equation 
    z_o2 / 0.21 = z_n2 / 0.79;
    z_h2o + z_o2 + z_n2 = 1.0;
    z_h2o = RH_env/100 * p_vap(T_env, 2)/p_env;
    
    h_h2o = h(T_env, 2, 1) + dhf(2, 1) - dhf(2,2);
    h_o2  = h(T_env, 3, 1);
    h_n2  = h(T_env, 5, 1);
    
    h_air = h_h2o*z_h2o + h_o2*z_o2 + h_n2*z_n2;
    
    c.z_local = {0.0, z_h2o, z_o2, 0.0, z_n2};
    c.h_local = h_air;
    
    annotation (Documentation(info="<html>
<p>This dummy object is used whenever a system flow goes to the environment, e.g.
the gas outlet of separators.</p>
</html>"), Icon(
        Polygon(points=[20,-40; 46,42; 74,-40; 20,-40], style(
            pattern=0,
            gradient=1,
            fillColor=58,
            rgbfillColor={0,127,0},
            fillPattern=8)),
        Polygon(points=[-80,-20; -10,-20; -44,82; -80,-20], style(
            pattern=0,
            gradient=1,
            fillColor=58,
            rgbfillColor={0,127,0},
            fillPattern=7)),
        Rectangle(extent=[-54,-20; -36,-74], style(
            pattern=0,
            gradient=1,
            fillColor=46,
            rgbfillColor={127,127,0})),
        Rectangle(extent=[40,-40; 52,-72], style(
            pattern=0,
            gradient=1,
            fillColor=46,
            rgbfillColor={127,127,0})),
        Ellipse(extent=[2,94; 48,50], style(
            pattern=0,
            gradient=3,
            fillColor=51,
            rgbfillColor={255,255,85}))),
      Diagram);
  end EnvironmentPort;
  
  partial model StirredTank "A generic stirred tank with an undefined shape." 
    
    import Modelica.SIunits.AmountOfSubstance;
    import Modelica.SIunits.Pressure;
    import Modelica.SIunits.Temperature;
    import Modelica.SIunits.Heat;
    import Modelica.SIunits.HeatCapacity;
    import Modelica.SIunits.HeatFlowRate;
    import Modelica.SIunits.MoleFraction;
    import Modelica.SIunits.CoefficientOfHeatTransfer;
    import Modelica.SIunits.Area;
    import Thermo.MolarEnthalpy;
    import Modelica.SIunits.Volume;
    import Modelica.Constants.R;
    
    import Thermo.p_vap;
    import Thermo.rho;
    import Thermo.cp;
    import Thermo.h;
    import Thermo.dhf;
    import Thermo.mw;
    import Thermo.AllSpecies;
    import Thermo.LiquidSpecies;
    import Thermo.GasSpecies;
    import Thermo.LiquidPhase;
    import Thermo.GasPhase;
    
    outer Temperature T_env "Environment temperature.";
    outer Pressure p_env "Environment pressure.";
    
    parameter String name = "Nameless" "Tank identifier.";
    parameter HeatCapacity Cp = 100 "Heat capacity of tank (glass, lid, ...).";
    parameter CoefficientOfHeatTransfer k_h = 0.0 
      "Default: perfect insulation.";
    parameter Area A_sur = 1e-3 
      "Surface area of contact between tank and environment";
    parameter Integer m = 1 "Number of flows, default must be != 0 for Dymola.";
    parameter Volume V = 1E-3 "Total volume of the tank.";
    
    Temperature T(start=T_env, fixed=true) "Representative tank temperature.";
    Heat Q "Heat exchanged with environment (usually negative).";
    CheckPoint[m] flows "Connections with other objects.";
    
    AmountOfSubstance[k] n(each min=0) "Moles of each species.";
    AmountOfSubstance n_tot(min=0) "Total moles.";
    MolarEnthalpy h_tot "Overall molar enthalpy in the tank.";
    MoleFraction[k] z(each min=0, each max=1) "Overall mole fractions";
    
    Volume V_l(min=0) "Amount of liquid volume.";
    AmountOfSubstance n_l(min=0) "Total moles of liquid.";
    MolarEnthalpy h_l "The molar enthalpy of the solution in liquid phase.";
    MoleFraction[k] x(each min=0, each max=1) "Molar fraction in liquid phase.";
    
    Volume V_g(min=0) "Amount of gaseous volume.";
    AmountOfSubstance n_g(min=0) "Total moles of gas.";
    MolarEnthalpy h_g "The molar enthalpy of the mixture in gas phase.";
    MoleFraction[k] y(each min=0, each max=1) "Molar fraction in gas phase.";
    
  protected 
    HeatFlowRate enthalpyInflow;
    HeatFlowRate sensibleHeat;
    HeatFlowRate latentHeat;
    HeatFlowRate newMassHeat;
    constant Real eps = 1e-8 "Constraint-violation tolerance (numerical noise)";
    constant Integer k = size(AllSpecies, 1);
    
  equation 
    /****************************
   *** MASS BALANCE SECTION ***
   ****************************/
    
    // Convenience variable for the total number of moles.
    n_tot = sum(n);
    
    // Overall molar fraction.
    z = n / n_tot;
    
    // Mass balance, piece of cake.
    der(n) = {sum(flows[j].F * flows[j].z[i] for j in 1:m) for i in AllSpecies};
    
    /******************************
   *** ENERGY BALANCE SECTION ***
   ******************************/
    
    // Complete energy balance.
    Q + enthalpyInflow = sensibleHeat + latentHeat + newMassHeat;
    
    // Exchange with environment.
    Q = A_sur * k_h * (T_env - T);
    
    // The sum of enthalpy inflows from all flows.
    enthalpyInflow = sum(flows[i].H for i in 1:m);
    
    // The specific heat (J/K) of the components in the tank, liquid and gas
    sensibleHeat = der(T) * (sum(cp(T,i,GasPhase) * y[i] * n_g for i in AllSpecies)+
                             sum(cp(T,i,LiquidPhase) * x[i] * n_l for i in LiquidSpecies)+
                             Cp);
    
    // The heat to bring water and methanol into gas phase
    latentHeat = sum((dhf(i,GasPhase) - dhf(i,LiquidPhase)) * der(y[i]*n_g) for i in LiquidSpecies);
    
    /* The heat spent on warming incoming matter to tank temperature T.
   * Note that methanol and water are a bit more complex as they change
   * phase according to temperature. */
    newMassHeat = sum(h(T,i,GasPhase) * der(n[i]) for i in GasSpecies)  +
                  sum(h(T,i,GasPhase) * der(y[i]*n_g) + h(T,i,LiquidPhase) * der(x[i]*n_l) for i in LiquidSpecies);
    
    /******************************
   *** MOLAR ENTHALPY SECTION ***
   ******************************/
    
    /* The molar enthalpy present in gas phase; includes heat of vaporisation
   * water and methanol. */
    h_g = (sum(h(T,i,GasPhase)*y[i] for i in AllSpecies) +
           sum((dhf(i,GasPhase)-dhf(i,LiquidPhase))*y[i] for i in LiquidSpecies));
    
    // The molar enthalpy present in liquid phase; includes only water and methanol.
    h_l = (sum(h(T,i,LiquidPhase)*x[i] for i in LiquidSpecies));
    
    // The overall molar enthalpy present in the tank.
    h_tot = (h_g*n_g + h_l*n_l)/n_tot;
    
    /************************************
   *** CHEMICAL EQUILIBRIUM SECTION ***
   ************************************/
    
    // Consistency: sum of gas and liquid volumes is constant.
    V = V_l + V_g;
    
    // The relationship between liquid volume, moles and compositions.
    // NOTE: assumed no mixture volume.
    V_l = n_l * sum(x[i]*mw(i)/rho(T,i,LiquidPhase) for i in LiquidSpecies);
    
    // The relationship between gaseous volume, moles and compositions.
    // NOTE: assumed no mixture volume.
    V_g = n_g * sum(y[i]*mw(i)/rho(T,i,GasPhase) for i in AllSpecies);
    
    // Molar fractions sum to 1.
    sum(x) = 1.0;
    sum(y) = 1.0;
    
    // Raoult's law for methanol and water.
    for i in LiquidSpecies loop
      p_env * y[i] = p_vap(T, i) * x[i];
    end for;
    // No gases in liquid phase.
    for i in GasSpecies loop
      x[i] = 0.0;
    end for;
    
    // Relation between fractions in phases and total moles.
    for i in AllSpecies loop
      n[i] = y[i] * n_g + x[i] * n_l;
    end for;
    
    /***************************
  *** SANITY CHECK SECTION ***
  ****************************/
    
    for i in AllSpecies loop
      assert( n[i] > -eps, "Negative amount of component "+String(i)+" in "+name+".");
    end for;
    
    assert( V_g > -eps, "Negative gas volume in "+name+".");
    assert( V_l > -eps, "Negative liquid volume in "+name+".");
    
    annotation (Documentation(info="<html>
<p>This class defines the general properties of all stirred tanks, that is a
single representative temperature, an array of amounts of substance with five
elements, the flows with which the tank is connected to other ones.</p>
<p>This class models:
<ul>
<li>Mass balance (set of differential equations);</li>
<li>Heat balance (differential equation);</li>
<li>Gas-liquid equilibrium (set of algebraic equations);</li>
<li>Sanity checks (various conditional statements).</li>
</ul>
</p>
<p>The number of flows <em>m</em> cannot be set to zero; if you need such a
isolated tank, set it to 1 and set <tt>flow[1].H = 0</tt>, <tt>flow[1].F = 
zeros(5)</tt>.</p>
</html>"), Icon);
  end StirredTank;
  
  partial model VerticalCylindricalTank 
    "A cylindrical tank in horizontal position." 
    extends StirredTank(name="Vertical cylindrical tank");
    import Modelica.SIunits.Area;
    import Modelica.SIunits.Length;
    
    parameter Area A "The cross-section of the tank.";
    Length level(min=0) "The liquid level, measured from the bottom.";
    
  equation 
      level = V_l / A;
    annotation (Documentation(info="<html>
<p>This class extends the generic stirred tank in specifying a geometry for the
tank, a vertical cylinder with a given horizontal section <em>A</em>.</p>
<p>This class contains a <em>level</em> variable, which is related to the degree
of filling of the tank by a linear expression.</p>
</html>"), Icon(           Rectangle(extent=[-40,0; 40,-80],  style(
            color=0,
            rgbcolor={0,0,0},
            thickness=2,
            fillColor=4,
            rgbfillColor={0,255,255},
            fillPattern=1)),
                           Rectangle(extent=[-40,80; 40,0], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=2,
            fillColor=7,
            rgbfillColor={255,255,255}))));
  end VerticalCylindricalTank;
  
  partial model HorizontalCylindricalTank 
    "A cylindrical tank in horizontal position." 
    extends StirredTank(name="Horizontal cylindrical tank");
    
    import Modelica.SIunits.Area;
    import Modelica.SIunits.Length;
    import Modelica.Constants.pi;
    import Modelica.Math.asin;
    
    parameter Area A "The cross-section of the tank";
    parameter Length d = 2.0*sqrt(A/pi) 
      "The inner diameter, or also height of the tank";
    Length level(min=0,max=d) "The liquid level, measured from the bottom.";
    
    /* This function provides the fraction of area of a circle that is
   * on a side of a given chord; the chord is identified with the
   * diameter-fraction argument x. */
  protected 
    function CircleAreaFraction 
        input Length x "Diameter Fraction";
        output Area y "Area fraction";
    algorithm 
        y := 1/pi * (asin(2*x-1)+2*(2*x-1)*sqrt(x-x^2) + pi/2);
    end CircleAreaFraction;
    
  equation 
    V_l = CircleAreaFraction(level/d) * V;
    
    annotation (Documentation(info="<html>
<p>This class extends the generic stirred tank in specifying a geometry for the
tank, a horizontal cylinder with a given horizontal section <em>A</em>.</p>
<p>This class contains a <em>level</em> variable, which is related to the degree
of filling of the tank by a fairly complex expression, which however does not
deviate significantly from linearity (maximum deviation is about 6%).</p>
</html>"), Icon(           Rectangle(extent=[-80,0; 80,-30],  style(
            color=0,
            rgbcolor={0,0,0},
            thickness=2,
            fillColor=67,
            rgbfillColor={85,255,255},
            fillPattern=1)),
                           Rectangle(extent=[-80,30; 80,0], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=2,
            fillColor=7,
            rgbfillColor={255,255,255}))));
  end HorizontalCylindricalTank;
  
  model FuelTank "The methanol fuel tank" 
    import Modelica.Constants.pi;
    extends VerticalCylindricalTank(name="Fuel tank", final m=2, V = 1E-3, A = (42E-3)^2*pi);
    
    CheckPoint bottomFlow = flows[1] annotation (extent=[-10,-90; 10,-70]);
    CheckPoint topFlow = flows[2] annotation (extent=[-10,70; 10,90]);
  equation 
  // TODO add if-equation to compensate in case of empty/full tank.
    flows[1].z_local = x;
    flows[1].h_local = h_l;
    flows[2].z_local = y;
    flows[2].h_local = h_g;
    
  initial equation 
    //x[1] = 0.2;  // FIXME why does the line below not work? Start bug hunt!
    n[1] / V_l = 1E3;  // 1 M liquid solution
    level =  0.5 * V / A;   // Half-full
    z[4] = 0.0;  // No CO2
    y[5] / 79 = y[3] / 21;  // N2/O2 in air proportions
    annotation (Icon);
  end FuelTank;
  
  model Mixer "Our laboratory mixer" 
    import Modelica.Constants.pi;
    extends VerticalCylindricalTank(name="Mixer", final m = 5, V = 5E-4, A = (37E-3)^2*pi);
    annotation (Documentation(info="<html>
<p>This is the mixer, which has four process streams (plus one connecting to the
atmosphere).</p>
<p>A set of initial equations is defined to determine initial conditions. These are:</p>
<ul>
<li>There is no methanol at the beginning;</li>
<li>Oxygen-to-Nitrogen ratio in gas phase as in dry air (0.21/0.79), but this does not
exclude presence of water in gas phase;</li>
<li>No CO<sub>2</sub> present;</li>
<li>The liquid level is half of the tank.</li>
</ul>
<p>Our laboratory's mixer has an internal diameter of 74 millimetres and a volume of
0.5 litres.</p>
</html>"),   Icon);
    
    CheckPoint flow1 = flows[1] "Connection at the bottom of the tank." 
                     annotation (extent=[-40,-90; -20,-70]);
    CheckPoint flow2 = flows[2] "Connection at the bottom of the tank." 
                     annotation (extent=[-20,-90; 0,-70]);
    CheckPoint flow3 = flows[3] "Connection at the bottom of the tank." 
                     annotation (extent=[0,-90; 20,-70]);
    CheckPoint flow4 = flows[4] "Connection at the bottom of the tank." 
                     annotation (extent=[20,-90; 40,-70]);
    CheckPoint topFlow = flows[5] 
      "Connection with the environment, at the top of the tank." 
      annotation (extent=[-10,70; 10,90]);
  equation 
    // The first four flows are connected with the liquid phase.
    for i in 1:4 loop
      flows[i].z_local = x;
      flows[i].h_local = h_l;
    end for;
    
    // Last connection is to atmosphere.
    flows[5].z_local = y;
    flows[5].h_local = h_g;
    
  initial equation 
    z[1] = 0.0; // (Almost) no methanol in initial solution
    level = 0.5 * V / A; // Half-full
    z[4] = 0.0; // No CO2 in gas
    y[5] / 0.79 = y[3] / 0.21; // N2/O2 in air proportions
    
  end Mixer;
  
  model Separator "A gas-liquid separator" 
    import Modelica.Constants.pi;
    extends HorizontalCylindricalTank( name="Separator", final m = 3, V = 0.25E-3, A = (23E-3)^2*pi);
    
    CheckPoint feed = flows[1] "Connection at mid-height in the tank." 
      annotation (extent=[-90,-10; -70,10]);
    CheckPoint gasOutlet = flows[2] "Connection at the top of the tank." 
      annotation (extent=[60,20; 80,40]);
    CheckPoint liquidOutlet = flows[3] "Connection on the bottom of the tank." 
      annotation (extent=[60,-40; 80,-20]);
  equation 
    // Set inlet composition in case of backflow: use average composition.
    flows[1].z_local = z;
    flows[1].h_local = h_tot;
    
    // Set gas outlet: always gas.
    flows[2].z_local = y;
    flows[2].h_local = h_g;
    
    // Set liquid outlet: always liquid.
    flows[3].z_local = x;
    flows[3].h_local = h_l;
    
  initial equation 
    z[1] = 0.0;  // No methanol in initial solution
    z[2] = 0.95;
    z[4] = 0.0;  // No CO2
    y[5] / 0.79 = y[3] / 0.21;  // N2/O2 in air proportions
    
    annotation (Documentation(info="<html>
<p>The separator is a stirred tank with an inlet and two outlets, separating liquid phase
from gas phase.</p>
<p>In case of return flow from the separator inlet, the gaseous composition and molar enthalpy
will be taken if the level is larger than zero; otherwise, the liquid ones will be used. If
this causes numerical problems, an approach might be attempted to make a mixture of these for
a brief range of level, or to use the overall composition <tt>z</tt>.</p>
<p>The gas outlet will output gas unless there is no more gas in the separator, in which case
it will overflow liquid; specularly the liquid outlet.</p>
<p>Our laboratory's separators have an internal diameter of 46 millimetres and a volume of
0.25 litres.</p>
</html>"), Icon);
  end Separator;
  
  model PipeSegment "A pipe segment, modelled as a tank." 
    extends StirredTank(m=2, name="Pipe segment");
  equation 
    for i in 1:m loop
      flows[i].z_local = z;
      flows[i].h_local = h_tot;
    end for;
    
    annotation (Documentation(info="<html>
<p>A pipe segment is modelled as a stirred tank, and the only additional piece of
information compared to the <tt>StirredTank</tt> class is the specification of the
tank compositions given to the <tt>CheckPoint</tt> connectors: these are the
overall composition <tt>z</tt> and enthalpy content <tt>h_tot</tt>, since it is
assumed that in a pipe the two phases travel together.</p>
<p>A more sophisticated implementation could make use of more complex formulations
to consider the drift of one phase with respect to the other.</p>
<p>A child class could declare some <tt>start</tt> values for composition, instead
of letting the initialisation algorithm do as it pleases; in will be useful in that
case to use the <tt>final</tt> keyword.</p>
</html>"));
  end PipeSegment;
  
  model Pipe "A series of PipeSegments or similar objects" 
    parameter Integer n = 10 "Number of segments";
    replaceable PipeSegment[n] segments(each V=1e-4);
    CheckPoint side1 = segments[1].flows[1] 
      annotation (extent=[-98,-10; -78,10]);
    CheckPoint side2 = segments[end].flows[2] 
                                            annotation (extent=[78,-10; 98,10]);
  protected 
    FlowConnector[n-1] connections;
  equation 
    for i in 1:(n-1) loop
      connect( segments[i].flows[2], connections[i].port1);
      connect( connections[i].port2, segments[i+1].flows[1]);
    end for;
    annotation (Documentation(info="<html>
<p>A pipe is a sequence of <tt>PipeSegment</tt> objects, intercalated by <tt>FlowConnector</tt>s.
The array of pipe segments is <tt>redeclarable</tt> so that the class may be used with different
implementations or children classes of <tt>PipeSegment</tt>.</p>
<p>To access the extreme connectors of the pipe, the connectors <tt>side1</tt> and <tt>side2</tt>
are available.</p>
</html>"), Icon(
        Rectangle(extent=[-80,20; 80,-20], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4,
            fillColor=7,
            rgbfillColor={255,255,255})),
        Line(points=[-60,20; -60,-20], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4)),
        Line(points=[-20,20; -20,-20], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4)),
        Line(points=[-40,20; -40,-20], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4)),
        Line(points=[20,20; 20,-20], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4)),
        Line(points=[40,20; 40,-20], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4)),
        Line(points=[60,20; 60,-20], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4)),
        Line(points=[0,20; 0,-20], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4))));
  end Pipe;
  
  model HeatExchangerPipeSegment 
    "A section of a pipe passing through a heat exchanger." 
    extends PipeSegment(m=3, name="Pipe section in heat exchanger");
    
  equation 
    /* These variables are unused, but must be initialised to some dummy
   * value, or the count of equations and variables will not square.*/
    flows[3].F = 0;
    flows[3].z = zeros(size(Thermo.AllSpecies,1));
    
    annotation (Documentation(info="<html>
<p>This kind of segments have a third connector through which no material flow
passes, rather only enthalpy flow, which can then be set to a function of
temperature differences.</p>
</html>"));
  end HeatExchangerPipeSegment;
  
  model HeatExchangerPipe "A pipe passing through an heat exchanger." 
    extends Pipe(redeclare HeatExchangerPipeSegment segments[n]) 
    annotation (Documentation(info="<html>
<p>This particular pipe consists of pipe elements that have a third connector,
which is used to transfer heat with other elements.</p>
</html>"));
  end HeatExchangerPipe;
  
  model HeatExchanger "A heat exchanger with two sides" 
    
  annotation (Icon(Ellipse(extent=[-60,60; 60,-60], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=2,
            fillColor=7,
            rgbfillColor={255,255,255})),
                           Line(points=[-60,0; -40,0; 0,40; 0,-40; 40,0; 60,0],
                style(
          color=0,
          rgbcolor={0,0,0},
          thickness=4,
          fillColor=67,
          rgbfillColor={85,255,255},
          fillPattern=1)),
        Text(
          extent=[10,70; 14,62],
          string="1",
          style(color=0, rgbcolor={0,0,0})),
        Text(
          extent=[-66,16; -62,8],
          string="1",
          style(color=0, rgbcolor={0,0,0})),
        Text(
          extent=[62,-10; 66,-18],
          style(color=0, rgbcolor={0,0,0}),
          string="2"),
        Text(
          extent=[-14,-60; -10,-68],
          style(color=0, rgbcolor={0,0,0}),
          string="2")),     Documentation(info="<html>
<p>This is a simple heat exchanger with two sides (left-right and top-bottom
in the figure) and no phase separation.</p>
</html>"));
    
    import Modelica.SIunits.Area;
    import Modelica.SIunits.Volume;
    import Modelica.SIunits.HeatFlowRate;
    import Modelica.SIunits.CoefficientOfHeatTransfer;
    
    CheckPoint f11 "Fluid 1, side1" 
                                  annotation (extent=[-70,-10; -50,10]);
    CheckPoint f21 "Fluid 2, side 1" 
                                   annotation (extent=[-10,50; 10,70]);
    CheckPoint f22 "Fluid 2, side 2" 
                                   annotation (extent=[-10,-70; 10,-50]);
    CheckPoint f12 "Fluid 1, side 2" 
                                   annotation (extent=[50,-10; 70,10]);
    
    parameter String name="Heat exchanger" "Name identifying the unit.";
    parameter Area A=1e-2 "The total heat transfer area of the exchanger.";
    parameter CoefficientOfHeatTransfer U = 10 "The heat-transfer coefficient.";
    parameter Integer steps = 5 
      "The number of subvolumes in which the two sides are divided.";
    parameter Volume V_1 = 1e-4 "The total volume of the first side.";
    parameter Volume V_2 = 1e-4 "The total volume of the second side.";
    
    HeatFlowRate Q "The heat moving from the first to the second side.";
    
    HeatExchangerPipe side1(n=steps, segments(each V=V_1/steps));
    HeatExchangerPipe side2(n=steps, segments(each V=V_2/steps));
    
  equation 
    // Connect pipe connectors to the heat exchanger's.
    connect( f11, side1.segments[1].flows[1]);
    connect( f12, side1.segments[end].flows[2]);
    connect( f21, side2.segments[1].flows[1]);
    connect( f22, side2.segments[end].flows[2]);
    
    /* Connect the pipe segments two by two, using the enthalpy element
   * in the third connectors to transport heat.*/
    for i in 1:steps loop
      side1.segments[i].flows[3].H + side2.segments[i].flows[3].H = 0;
      side1.segments[i].flows[3].H = (A/steps) * U * (side2.segments[i].T - side1.segments[i].T);
    end for;
    
    Q = sum({side2.segments[i].flows[3].H for i in 1:steps});
    
  end HeatExchanger;
  
  model FlowController "A unit modelling a pump or a MFC" 
    
    import Thermo.mw;
    import Thermo.AllSpecies;
    
    MolarFlowRate F;
    Modelica.SIunits.MassFlowRate m;
    CheckPoint inlet "Unit inlet" annotation (extent=[-12,-10; 8,10]);
    CheckPoint outlet "Unit outlet" annotation (extent=[-10,50; 10,70]);
  equation 
    m = sum({inlet.z[i] * F * mw(i) for i in AllSpecies});
    F = inlet.F;
    connect( inlet, outlet);
    
    annotation (Icon(
        Ellipse(extent=[-60,60; 60,-60], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=2,
            fillColor=7,
            rgbfillColor={255,255,255}))),
                              Documentation(info="<html>
<p>This class allows to set a certain overall molar or mass flow. It is not
immediately possible to set a <em>volume</em> flow, because this would entail
calculating the phase equilibrium, which we are not doing here (though child
classes could specialize).</p>
</html>"));
  end FlowController;
  
  model GasFlowController "A flow controller with only gas phase" 
    extends FlowController;
    annotation (Icon,     Documentation(info="<html>
<p>This class implements a mass flow controller with field volumetric units. Since
there are two different standards (the actual \"standard\" at 0� Celsius and the Norm
at 70� Fahrenheit), it is necessary to adjust the reference temperature; the default
assumes zero Celsius (\"standard\" value).</p>
<p>The flow assumes that all components are in gas phase and takes their density from
the Thermo library, where the ideal gas law is (usually) assumed.</p>
</html>"));
    import Thermo.rho;
    import Thermo.mw;
    import Thermo.AllSpecies;
    import Thermo.GasPhase;
    import Modelica.SIunits.VolumeFlowRate;
    import Modelica.SIunits.Temperature;
    import Modelica.SIunits.Conversions.from_degC;
    import Modelica.SIunits.Conversions.from_degF;
    
    VolumeFlowRate V "Volumetric flow rate.";
    parameter Temperature T = 273.15 "Reference temperature for standard flow.";
  equation 
    V = sum({inlet.z[i] * F * mw(i) / rho(T, i, GasPhase) for i in AllSpecies});
    
  end GasFlowController;
  
  model Pump "A liquid pump" 
    extends FlowController;
    annotation (Icon,        Documentation(info="<html>
<p>This class implements a liquid pump with field volumetric units. It will be
necessary to somehow specify its operating temperature to make it work: for
example, depending on flow direction, it might be the temperature of the element 
upstream or the one downstream.</p>
<p>The flow assumes that only water and methanol are present and are completely
in liquid phase; it takes their density from the Thermo library.</p>
</html>"));
    import Thermo.rho;
    import Thermo.mw;
    import Thermo.LiquidPhase;
    import Thermo.LiquidSpecies;
    import Modelica.SIunits.VolumeFlowRate;
    import Modelica.SIunits.Temperature;
    import Modelica.SIunits.Conversions.from_degC;
    import Modelica.SIunits.Conversions.from_degF;
    
    outer parameter Temperature T_env;
    parameter Temperature T = T_env 
      "Temperature (FIXME should not be a reference!).";
    VolumeFlowRate V "Standard volume flow rate.";
    
  equation 
      V = sum({inlet.z[i] * F * mw(i) / rho(T, i, LiquidPhase) for i in LiquidSpecies});
    
  end Pump;
  
  package Tests "A series of tests for units defined in the Tank library" 
    
    model TestFuelTank "A simple test for the FuelTank class" 
      
      FlowConnector flowConnector annotation (extent=[-32,-24; -2,8]);
      EnvironmentPort environmentPort annotation (extent=[28,-10; 50,16]);
      annotation (Diagram);
      FuelTank fuelTank(T(
                        start = 350),k_h=1000) annotation (extent=[-80,6; -40,44]);
      FlowConnector flowConnector1 annotation (extent=[-26,36; -2,62]);
      
      inner parameter Real RH_env = 50;
      inner parameter Modelica.SIunits.Pressure p_env = 101325;
      inner parameter Modelica.SIunits.Temperature T_env = 298.15;
      
      EnvironmentPort environmentPort1 annotation (extent=[18,48; 48,78]);
    equation 
      
      fuelTank.bottomFlow.F = -1;
      
      connect(flowConnector1.port1, fuelTank.topFlow) annotation (points=[-19.76,
            49; -60,49; -60,40.2], style(pattern=0, thickness=2));
      connect(fuelTank.bottomFlow, flowConnector.port1) annotation (points=[-60,9.8;
            -60,-8; -24.2,-8],      style(pattern=0, thickness=2));
      connect(flowConnector.port2, environmentPort.c) annotation (points=[-12.8,-8;
            18,-8; 18,-3.5; 29.1,-3.5],
                                      style(pattern=0, thickness=2));
      connect(environmentPort1.c, flowConnector1.port2) annotation (points=[19.5,
            55.5; 6,55.5; 6,49; -10.64,49],style(pattern=0, thickness=2));
    end TestFuelTank;
    
    model TestMixer "A Vodka mixer" 
      
      Mixer mixer annotation (extent=[-50,0; -10,42]);
      annotation (Diagram);
      EnvironmentPort environmentPort1 annotation (extent=[-2,48; 16,64]);
      FlowConnector flowConnector annotation (extent=[-66,-20; -46,0]);
      FlowConnector flowConnector1 annotation (extent=[-28,42; -8,62]);
      EnvironmentPort environmentPort2 annotation (extent=[4,-56; 24,-36]);
      Plug plug annotation (extent=[-16,-10; 4,10]);
      Plug plug1 annotation (extent=[-14,-30; 6,-10]);
      FlowConnector flowConnector3 annotation (extent=[-26,-64; -4,-38]);
      inner parameter Real RH_env = 30;
      inner parameter Modelica.SIunits.Pressure p_env = 101325;
      inner parameter Modelica.SIunits.Temperature T_env = 298.15;
      
      MethanolSolution source annotation (extent=[-86,8; -66,28]);
    equation 
      connect(flowConnector.port2, mixer.flow1) annotation (points=[-53.2,-10;
            -40,-10; -40,4.2; -36,4.2], style(pattern=0, thickness=2));
      connect(flowConnector1.port1, mixer.topFlow) annotation (points=[-22.8,52;
            -30,52; -30,37.8], style(pattern=0, thickness=2));
      connect(flowConnector1.port2, environmentPort1.c) annotation (points=[-15.2,52;
            -1.1,52],                 style(pattern=0, thickness=2));
      connect(plug.c, mixer.flow4) annotation (points=[-12.2,6.10623e-16; -24,
            6.10623e-16; -24,4.2],
          style(
          pattern=0,
          thickness=2,
          fillColor=46,
          rgbfillColor={127,127,0},
          fillPattern=7));
      connect(plug1.c, mixer.flow3) annotation (points=[-10.2,-20; -28,-20; -28,
            4.2],
          style(
          pattern=0,
          thickness=2,
          fillColor=46,
          rgbfillColor={127,127,0},
          fillPattern=7));
      connect(flowConnector3.port1, mixer.flow2) annotation (points=[-20.28,-51;
            -32,-51; -32,4.2], style(
          pattern=0,
          thickness=2,
          fillColor=46,
          rgbfillColor={127,127,0},
          fillPattern=7));
      connect(flowConnector3.port2, environmentPort2.c) annotation (points=[-11.92,
            -51; 5,-51], style(
          pattern=0,
          thickness=2,
          fillColor=46,
          rgbfillColor={127,127,0},
          fillPattern=7));
      flowConnector.port2.F = 1;
      flowConnector3.port2.F = 0.5;
      
      connect(source.p, flowConnector.port1) annotation (points=[-76,18; -76,
            -10; -60.8,-10], style(pattern=0, thickness=2));
    end TestMixer;
    
    model TestSeparator "Test suite for the separator model." 
      
      inner parameter Real RH_env = 30;
      inner parameter Modelica.SIunits.Pressure p_env = 101325;
      inner parameter Modelica.SIunits.Temperature T_env = 298.15;
      Separator separator annotation (extent=[-36,-24; 16,16]);
      EnvironmentPort environmentPort annotation (extent=[48,20; 72,44]);
      EnvironmentPort environmentPort1 annotation (extent=[46,-40; 68,-16]);
      FlowConnector BottomSeparatorConnector 
                                  annotation (extent=[18,-44; 38,-24]);
      FlowConnector TopSeparatorConnector 
                                   annotation (extent=[20,16; 40,36]);
      FlowConnector FuelTankToSeparatorConnector 
                                   annotation (extent=[-58,-14; -38,6]);
      MethanolSolution source(T=350) annotation (extent=[-72,4; -52,24]);
    equation 
      
      annotation (Diagram);
      connect(BottomSeparatorConnector.port2, environmentPort1.c) 
                                                       annotation (points=[30.8,-34;
            47.1,-34],                    style(
          pattern=0,
          thickness=2,
          fillPattern=1));
      connect(BottomSeparatorConnector.port1, separator.liquidOutlet) 
                                                           annotation (points=[23.2,-34;
            8.2,-34; 8.2,-10],               style(
          pattern=0,
          thickness=2,
          fillPattern=1));
      connect(TopSeparatorConnector.port2, environmentPort.c) 
                                                       annotation (points=[32.8,26;
            49.2,26],                     style(
          pattern=0,
          thickness=2,
          fillPattern=1));
      connect(TopSeparatorConnector.port1, separator.gasOutlet) 
                                                         annotation (points=[25.2,26;
            8,26; 8,2; 8.2,2],              style(
          pattern=0,
          thickness=2,
          fillPattern=1));
      connect(FuelTankToSeparatorConnector.port2, separator.feed) 
                                                    annotation (points=[-45.2,-4;
            -41.6,-4; -41.6,-4; -38,-4; -38,-4; -30.8,-4],
                                               style(
          pattern=0,
          thickness=2,
          fillPattern=1));
      
      source.p.F = -1;
      der(separator.level) = 0;
      connect(source.p, FuelTankToSeparatorConnector.port1) annotation (points=[-62,14;
            -62,-4; -52.8,-4],     style(pattern=0, thickness=2));
    end TestSeparator;
    
    model TestPipe 
      inner parameter Real RH_env = 30;
      inner parameter Modelica.SIunits.Pressure p_env = 101325;
      inner parameter Modelica.SIunits.Temperature T_env = 298.15;
      Pipe pipe(n=10);
      MethanolSolution source(T=320);
      FlowConnector inlet;
      EnvironmentPort nature;
      FlowConnector outlet;
    equation 
      connect(source.p, inlet.port1);
      connect(inlet.port2, pipe.segments[1].flows[1]);
      connect(pipe.segments[end].flows[end], outlet.port1);
      connect(outlet.port2, nature.c);
      
      source.p.F = -1;
    initial equation 
      for i in 1:pipe.n loop
        pipe.segments[i].z[1]=0.01;
        pipe.segments[i].z[2]=0.9;
        pipe.segments[i].z[4]=0;
        pipe.segments[i].z[3]/21 = pipe.segments[i].z[5]/79;
      end for;
      
    end TestPipe;
    
    model TestHeatExchangerPipe 
      inner parameter Real RH_env = 30;
      inner parameter Modelica.SIunits.Pressure p_env = 101325;
      inner parameter Modelica.SIunits.Temperature T_env = 298.15;
      HeatExchangerPipe pipe(n=10);
      MethanolSolution source(T=320);
      FlowConnector inlet;
      EnvironmentPort nature;
      FlowConnector outlet;
    equation 
      connect(source.p, inlet.port1);
      connect(inlet.port2, pipe.segments[1].flows[1]);
      connect(pipe.segments[end].flows[2], outlet.port1);
      connect(outlet.port2, nature.c);
      
      source.p.F = -1;
      for i in 1:pipe.n loop
        pipe.segments[i].flows[3].H = 20;
      end for;
      
    initial equation 
      for i in 1:pipe.n loop
        pipe.segments[i].z[1]=0;
        pipe.segments[i].z[2]=0.95;
        pipe.segments[i].z[4]=0;
        pipe.segments[i].z[3]/21 = pipe.segments[i].z[5]/79;
      end for;
      
    end TestHeatExchangerPipe;
    
    model TestHeatExchangerPipe_Gas 
      
      inner parameter Real RH_env = 30;
      inner parameter Modelica.SIunits.Pressure p_env = 101325;
      inner parameter Modelica.SIunits.Temperature T_env = 298.15;
      HeatExchangerPipe pipe(n=10,segments(T(each start = 290)));
      EnvironmentPort source;
      FlowConnector inlet;
      EnvironmentPort nature;
      FlowConnector outlet;
    equation 
      connect(source.c, inlet.port1);
      connect(inlet.port2, pipe.segments[1].flows[1]);
      connect(pipe.segments[end].flows[2], outlet.port1);
      connect(outlet.port2, nature.c);
      
      source.c.F = -0.1;
      for i in 1:pipe.n loop
        pipe.segments[i].flows[3].H = 20;
      end for;
      
    initial equation 
      for i in 1:pipe.n loop
        pipe.segments[i].z[1]=0;
        pipe.segments[i].z[2]=0.99;
        pipe.segments[i].n[4]=0.0;
        pipe.segments[i].z[3]/21 = pipe.segments[i].z[5]/79;
      end for;
      
    end TestHeatExchangerPipe_Gas;
    
    model TestHeatExchanger 
      
      inner parameter Real RH_env = 30;
      inner parameter Modelica.SIunits.Pressure p_env = 101325;
      inner parameter Modelica.SIunits.Temperature T_env = 298.15;
      
      HeatExchanger hx(steps=3, side1(segments(T(each start=320))),
        U=5000)                   annotation (extent=[-16,-4; 12,24]);
      MethanolSolution methanolSolution(T=340) 
                                        annotation (extent=[-64,0; -44,20]);
      annotation (Diagram);
      EnvironmentPort solutionOutlet  annotation (extent=[36,6; 52,22]);
      FlowConnector to_f11        annotation (extent=[-42,0; -22,20]);
      FlowConnector to_f12         annotation (extent=[12,0; 32,20]);
      FlowConnector to_f21    annotation (extent=[0,36; 20,56]);
      EnvironmentPort coolingInlet     annotation (extent=[26,-26; 44,-10]);
      FlowConnector to_f22   annotation (extent=[2,-32; 22,-12]);
      EnvironmentPort coolingOutlet   annotation (extent=[30,42; 46,58]);
    equation 
      connect(to_f11.port2, hx.f11)                   annotation (points=[-29.2,10;
            -19.8,10; -19.8,10; -10.4,10],   style(pattern=0, thickness=2));
      connect(hx.f12, to_f12.port1)                    annotation (points=[6.4,10;
            11.8,10; 11.8,10; 17.2,10],     style(pattern=0, thickness=2));
      connect(to_f12.port2, solutionOutlet.c) 
        annotation (points=[24.8,10; 36.8,10], style(pattern=0, thickness=2));
      connect(to_f22.port2, coolingInlet.c)     annotation (points=[14.8,-22;
            26.9,-22], style(pattern=0, thickness=2));
      connect(methanolSolution.p, to_f11.port1) 
        annotation (points=[-54,10; -36.8,10], style(pattern=0, thickness=2));
      methanolSolution.p.F = -1;
      coolingInlet.c.F = -10;
      
    initial equation 
      for i in 1:hx.steps loop
        hx.side1.segments[i].z[1]=0.001;
        hx.side1.segments[i].z[3]=0.001;
        hx.side1.segments[i].z[4]=0.001;
        hx.side1.segments[i].z[5]=0.001;
        
        hx.side2.segments[i].z[1]=0.001;
        hx.side2.segments[i].z[2]=0.001;
        hx.side2.segments[i].z[3]/21 = hx.side2.segments[i].z[5]/79;
        hx.side2.segments[i].z[4]=0.001;
      end for;
      
    equation 
      connect(to_f22.port1, hx.f22) annotation (points=[7.2,-22; -2,-22; -2,1.6],
          style(pattern=0, thickness=2));
      connect(to_f21.port2, coolingOutlet.c) 
        annotation (points=[12.8,46; 30.8,46], style(pattern=0, thickness=2));
      connect(to_f21.port1, hx.f21) annotation (points=[5.2,46; -2,46; -2,18.4],
          style(pattern=0, thickness=2));
    end TestHeatExchanger;
    
    model TestEnvironment 
      
      annotation (Diagram);
      inner parameter Real RH_env = 60;
      inner parameter Modelica.SIunits.Pressure p_env = 101325;
      inner parameter Modelica.SIunits.Temperature T_env = 298.15;
      Pipe pipe(n=2) annotation (extent=[-30,-10; 10,26]);
      EnvironmentPort inlet annotation (extent=[-54,24; -34,44]);
      FlowConnector enteringConnector annotation (extent=[-56,-2; -36,18]);
      FlowConnector exitingConnector annotation (extent=[26,-2; 46,18]);
      EnvironmentPort outlet annotation (extent=[52,2; 74,26]);
    equation 
      connect(exitingConnector.port1, pipe.side2) annotation (points=[31.2,8;
            19.4,8; 19.4,8; 7.6,8], style(pattern=0, thickness=2));
      connect(enteringConnector.port2, pipe.side1) annotation (points=[-43.2,8;
            -35.4,8; -35.4,8; -27.6,8], style(pattern=0, thickness=2));
      connect(inlet.c, enteringConnector.port1) annotation (points=[-53,29;
            -64.5,29; -64.5,8; -50.8,8], style(pattern=0, thickness=2));
      connect(outlet.c, exitingConnector.port2) 
        annotation (points=[53.1,8; 38.8,8], style(pattern=0, thickness=2));
      inlet.c.F = -1;
    initial equation 
      for i in 1:pipe.n loop
        pipe.segments[i].z[1]=0.1;
        pipe.segments[i].z[2]=0.1;
        pipe.segments[i].z[4]=0.1;
        pipe.segments[i].z[3]/21 = pipe.segments[i].z[5]/79;
      end for;
      
    end TestEnvironment;
    
    model SystemWithoutFC 
      Mixer mixer annotation (extent=[-70,-86; -30,-62]);
      FuelTank fuelTank annotation (extent=[32,-88; 52,-68]);
      EnvironmentPort environmentPort annotation (extent=[70,-70; 88,-52]);
      FlowConnector flowConnector annotation (extent=[48,-76; 68,-56]);
      annotation (Diagram);
      FlowConnector flowConnector1 annotation (extent=[-24,-96; -4,-76]);
      Pump pump annotation (extent=[-12,-102; 8,-82]);
      FlowConnector flowConnector2 annotation (extent=[12,-102; 32,-82]);
      EnvironmentPort environmentPort1 annotation (extent=[-36,-52; -20,-32]);
      Pump pump1 annotation (extent=[-84,-64; -64,-44]);
      Pipe pipe annotation (extent=[-50,-16; -18,20]);
      HeatExchanger heatExchanger annotation (extent=[-4,-8; 16,12]);
      FlowConnector flowConnector3 annotation (extent=[-74,-102; -54,-82]);
      FlowConnector flowConnector4 annotation (extent=[-70,-8; -50,12]);
      FlowConnector flowConnector5 annotation (extent=[-18,-8; 2,12]);
      FlowConnector flowConnector6 annotation (extent=[-52,-56; -34,-38]);
      EnvironmentPort environmentPort2 annotation (extent=[20,-32; 38,-14]);
      GasFlowController gasFlowController annotation (extent=[-18,-38; 2,-18]);
      FlowConnector flowConnector7 annotation (extent=[-12,-24; 8,-4]);
      FlowConnector flowConnector8 annotation (extent=[-2,-38; 18,-18]);
      FlowConnector flowConnector9 annotation (extent=[6,8; 26,28]);
      EnvironmentPort environmentPort3 annotation (extent=[24,14; 40,30]);
      Separator separator annotation (extent=[32,-12; 62,16]);
      FlowConnector flowConnector10 annotation (extent=[16,-8; 36,12]);
      FlowConnector flowConnector11 annotation (extent=[4,-90; 24,-70]);
      FlowConnector flowConnector12 annotation (extent=[56,4; 76,24]);
      EnvironmentPort environmentPort4 annotation (extent=[72,10; 88,26]);
    equation 
      connect(flowConnector.port2, environmentPort.c) annotation (points=[60.8,
            -66; 70.9,-66; 70.9,-65.5], style(pattern=0, thickness=2));
      connect(flowConnector.port1, fuelTank.topFlow) annotation (points=[53.2,
            -66; 42,-66; 42,-70], style(pattern=0, thickness=2));
      connect(flowConnector1.port1, mixer.flow2) annotation (points=[-18.8,-86;
            -26,-86; -26,-98; -52,-98; -52,-83.6], style(pattern=0, thickness=2));
      connect(pump.outlet, flowConnector1.port2) 
        annotation (points=[-2,-86; -11.2,-86], style(pattern=0, thickness=2));
      connect(pump.inlet, flowConnector2.port1) annotation (points=[-2.2,-92;
            17.2,-92], style(pattern=0, thickness=2));
      connect(flowConnector2.port2, fuelTank.bottomFlow) annotation (points=[
            24.8,-92; 42,-92; 42,-86], style(pattern=0, thickness=2));
      connect(pump1.inlet, flowConnector3.port1) annotation (points=[-74.2,-54;
            -74,-54; -74,-92; -68.8,-92], style(pattern=0, thickness=2));
      connect(flowConnector3.port2, mixer.flow1) annotation (points=[-61.2,-92;
            -56,-92; -56,-83.6], style(pattern=0, thickness=2));
      connect(flowConnector4.port1, pump1.outlet) annotation (points=[-64.8,2;
            -74,2; -74,-48], style(pattern=0, thickness=2));
      connect(flowConnector4.port2, pipe.side1) annotation (points=[-57.2,2;
            -54.92,2; -54.92,2; -52.64,2; -52.64,2; -48.08,2], style(pattern=0,
            thickness=2));
      connect(pipe.side2, flowConnector5.port1) annotation (points=[-19.92,2;
            -18.14,2; -18.14,2; -16.36,2; -16.36,2; -12.8,2], style(pattern=0,
            thickness=2));
      connect(flowConnector5.port2, heatExchanger.f11) annotation (points=[-5.2,2;
            -3.9,2; -3.9,2; -2.6,2; -2.6,2; 5.55112e-16,2],    style(pattern=0,
            thickness=2));
      connect(flowConnector6.port1, mixer.topFlow) annotation (points=[-47.32,
            -47; -50,-47; -50,-64.4], style(pattern=0, thickness=2));
      connect(flowConnector7.port2, heatExchanger.f22) annotation (points=[0.8,-14;
            6,-14; 6,-4],      style(pattern=0, thickness=2));
      connect(flowConnector7.port1, gasFlowController.outlet) annotation (
          points=[-6.8,-14; -8,-14; -8,-22], style(pattern=0, thickness=2));
      connect(flowConnector8.port1, gasFlowController.inlet) 
        annotation (points=[3.2,-28; -8.2,-28], style(pattern=0, thickness=2));
      connect(flowConnector8.port2, environmentPort2.c) annotation (points=[
            10.8,-28; 12,-28; 12,-27.5; 20.9,-27.5], style(pattern=0, thickness=
             2));
      connect(flowConnector9.port1, heatExchanger.f21) annotation (points=[11.2,18;
            6,18; 6,8],     style(pattern=0, thickness=2));
      connect(environmentPort3.c, flowConnector9.port2) 
        annotation (points=[24.8,18; 18.8,18], style(pattern=0, thickness=2));
      connect(flowConnector10.port1, heatExchanger.f12) annotation (points=[21.2,2;
            18.9,2; 18.9,2; 16.6,2; 16.6,2; 12,2],         style(pattern=0,
            thickness=2));
      connect(flowConnector10.port2, separator.feed) annotation (points=[28.8,2;
            30.35,2; 30.35,2; 31.9,2; 31.9,2; 35,2],
                                         style(pattern=0, thickness=2));
      connect(mixer.flow3, flowConnector11.port1) annotation (points=[-48,-83.6;
            -48,-94; -30,-94; -30,-80; 9.2,-80], style(pattern=0, thickness=2));
      connect(flowConnector11.port2, separator.liquidOutlet) annotation (points=[16.8,-80;
            32,-80; 32,-38; 57.5,-38; 57.5,-2.2],           style(pattern=0,
            thickness=2));
      connect(flowConnector6.port2, environmentPort1.c) annotation (points=[
            -40.48,-47; -38,-47; -38,-47; -35.2,-47], style(pattern=0,
            thickness=2));
      connect(environmentPort4.c, flowConnector12.port2) 
        annotation (points=[72.8,14; 68.8,14], style(pattern=0, thickness=2));
      connect(flowConnector12.port1, separator.gasOutlet) annotation (points=[61.2,14;
            58,14; 58,6.2; 57.5,6.2],          style(pattern=0, thickness=2));
    end SystemWithoutFC;
  end Tests;
end Tank;
