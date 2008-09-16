

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
    
    MolarFlowRate F;
    MoleFraction z[5];
    EnthalpyFlowRate H;
    MoleFraction z_tank[5];
    MolarEnthalpy h_tank;
    
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
    
    annotation (Icon(Polygon(points=[-60,40; -60,-40; 40,40; 40,-40; -60,40],
            style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4))));
    CheckPoint c annotation (extent=[-72,-10; -52,10]);
  equation 
    c.F = 0;
    c.z = zeros(5);
    c.H = 0;
    
  end Plug;
  
  model FlowConnector "A class to determine flow direction" 
  public 
    CheckPoint port1 annotation (extent=[-92,-10; -72,10]);
    CheckPoint port2 annotation (extent=[72,-10; 92,10]);
  equation 
    // Connection part.
    port1.F + port2.F = 0;
    port1.H + port2.H = 0;
    port1.z = port2.z;
    
    // Determining which side is taken for the values of composition and molar enthalpy.
    if port1.F >= 0 then  // If flow is entering into the tank on side one, or is zero:
      port2.H = port2.F * port2.h_tank;
      port2.z = port2.z_tank;
    else // If flow is entering from the second connection
      port1.H = port1.F * port1.h_tank;
      port1.z = port1.z_tank;
    end if;
    
    annotation (Documentation(info="<html>
<p><tt>FlowConnector</tt>s are objects that operate similarly to Modelica connectors,
but contain additional equations. They have two <tt>CheckPoint</tt> connectors which
lead to two tanks. Once given the direction of the overall molar flow, its task is to
decide which composition <tt>z_tank</tt> and specific molar enthalpy <tt>h_tank</tt>
should be used to set the flow's composition <tt>z</tt> and enthalpic flow <tt>H</tt>.</p>
<p>Each subclass of <tt>StirredTank</tt> should set the values of <tt>z_tank</tt> and
<tt>h_tank</tt> for each <tt>CheckPoint</tt> according to which phase would be drawn
from the connection (gas, liquid, mixed or depending on some variables, such as
liquid level).</p>
</html>"), Icon(
        Ellipse(extent=[40,40; 80,-40], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=2,
            arrow=3,
            fillColor=7,
            rgbfillColor={255,255,255},
            fillPattern=1)),
        Rectangle(extent=[-60,40; 60,-40], style(
            pattern=0,
            thickness=2,
            arrow=3,
            fillColor=7,
            rgbfillColor={255,255,255},
            fillPattern=1)),
        Line(points=[-60,40; 60,40], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=2,
            fillColor=0,
            rgbfillColor={0,0,0},
            fillPattern=1)),
        Line(points=[-60,-40; 60,-40], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=2,
            fillColor=0,
            rgbfillColor={0,0,0},
            fillPattern=1)),
        Ellipse(extent=[-80,40; -40,-40], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=2,
            arrow=3,
            fillColor=7,
            rgbfillColor={255,255,255})),
        Polygon(points=[-30,-2; -10,-32; -10,28; -30,-2], style(
            pattern=0,
            fillColor=0,
            rgbfillColor={0,0,0})),
        Polygon(points=[46,-2; 26,-32; 26,28; 46,-2], style(
            pattern=0,
            fillColor=0,
            rgbfillColor={0,0,0}))),
      Diagram);
  end FlowConnector;
  
  annotation (uses(Modelica(version="2.2.1")), Documentation(info="<html>
<p>This package contains various models related to fluid flow in stirred tanks.</p>
</html>"));
  
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
    
    MolarEnthalpy h_air "The molar enthalpy of air in the current conditions.";
    
    CheckPoint c annotation (extent=[-100,-100; -80,-80]);
  equation 
    z_o2 / 0.21 = z_n2 / 0.79;
    z_h2o + z_o2 + z_n2 = 1.0;
    z_h2o = RH_env/100 * p_vap(T_env, 1)/p_env;
    
    h_air = ( dhf(2, 1) - dhf(2,2) + h(T_env, 2, 1) - h(298.15, 2, 1))  * z_h2o
            + ( h(T_env, 3, 1) - h(298.15, 3, 1))  * z_o2
            + ( h(T_env, 5, 1) - h(298.15, 5, 1))  * z_n2;
    
    c.z_tank = {0.0, z_h2o, z_o2, 0.0, z_n2};
    c.h_tank = h_air;
    
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
    
    outer Temperature T_env "Environment temperature.";
    outer Pressure p_env "Environment pressure.";
    
    parameter String name = "Nameless" "Tank identifier.";
    parameter HeatCapacity Cp = 1000 "Heat capacity of tank (glass, lid, ...).";
    parameter CoefficientOfHeatTransfer h_k = 0.0 
      "Default: perfect insulation.";
    parameter Area A_sur = 1e-3 
      "Surface area of contact between tank and environment";
    parameter Integer m = 1 "Number of flows, default must be != 0 for Dymola.";
    parameter Volume V = 1E-3 "Total volume of the tank.";
    
    Temperature T(start=T_env, fixed=true) "Representative tank temperature.";
    Heat Q "Heat exchanged with environment (usually negative).";
    CheckPoint[m] flows "Connections with other objects.";
    
    AmountOfSubstance[5] n "Moles of each species.";
    AmountOfSubstance n_tot "Total moles.";
    MolarEnthalpy h_tot "Overall molar enthalpy in the tank.";
    MoleFraction[5] z "Overall mole fractions";
    
    Volume V_l "Amount of liquid volume.";
    AmountOfSubstance n_l "Total moles of liquid.";
    MolarEnthalpy h_l "The molar enthalpy of the solution in liquid phase.";
    MoleFraction[5] x "Molar fraction in liquid phase.";
    
    Volume V_g "Amount of gaseous volume.";
    AmountOfSubstance n_g "Total moles of gas.";
    MolarEnthalpy h_g "The molar enthalpy of the mixture in gas phase.";
    MoleFraction[5] y "Molar fraction in gas phase.";
    
  protected 
    HeatFlowRate enthalpyInflow;
    HeatFlowRate sensibleHeat;
    HeatFlowRate latentHeat;
    HeatFlowRate dmHeat;
    
  equation 
    /****************************
   *** MASS BALANCE SECTION ***
   ****************************/
    
    // Convenience variable for the total number of moles.
    n_tot = sum(n);
    
    // Overall molar fraction.
    z = n / n_tot;
    
    // Mass balance, piece of cake.
    der(n) = {sum(flows[j].F * flows[j].z[i] for j in 1:m) for i in 1:5};
    
    /******************************
   *** ENERGY BALANCE SECTION ***
   ******************************/
    
    // Complete energy balance.
    Q + enthalpyInflow = Cp * der(T) + sensibleHeat + latentHeat;
    
    // Exchange with environment.
    Q = A_sur * h_k * (T_env - T);
    
    // The sum of enthalpy inflows from all flows.
    enthalpyInflow = sum(flows[i].H for i in 1:m);
    
    // The specific heat (J/K) of the components in the tank, liquid and gas
    sensibleHeat = der(T) * (sum(cp(T,i,1) * y[i] * n_g for i in 1:5)
                   + sum(cp(T,i,2) * x[i] * n_l for i in 1:2)) + dmHeat;
    
    /* The heat spent on warming incoming matter to tank temperature T.
   * Note that methanol and water are a bit more complex as they change
   * phase according to temperature. */
    dmHeat = sum( (h(T,i,1)-h(298.15,i,1)) * der(n[i]) for i in 3:5)  +
             sum( (h(T,i,1)-h(298.15,i,1)) * der(y[i]*n_g) +
                  (h(T,i,2)-h(298.15,i,2)) * der(x[i]*n_l) for i in 1:2);
    
    // The heat to bring water and methanol into gas phase
    latentHeat = sum((dhf(i,1) - dhf(i,2)) * der(y[i]*n_g) for i in 1:2);
    
    /******************************
   *** MOLAR ENTHALPY SECTION ***
   ******************************/
    
    /* The molar enthalpy present in gas phase; includes heat of vaporisation
   * water and methanol. */
    h_g = (sum(h(T,i,1) * y[i] for i in 1:5) + sum((dhf(i,1)-dhf(i,2))*y[i] for i in 1:2));
    
    // The molar enthalpy present in liquid phase; includes only water and methanol.
    h_l = (sum(h(T,i,2) * x[i] for i in 1:2));
    
    // The overall molar enthalpy present in the tank.
    h_tot = ( h_g*n_g + h_l*n_l)  / n_tot;
    
    /************************************
   *** CHEMICAL EQUILIBRIUM SECTION ***
   ************************************/
    
    // Consistency: sum of gas and liquid volumes is constant.
    V = V_l + V_g;
    
    // The relationship between liquid volume, moles and compositions.
    // NOTE: assumed no mixture volume.
    V_l = n_l * sum(x[i]*mw(i)/rho(T,i,2) for i in 1:2);
    
    // The relationship between gaseous volume, moles and compositions.
    // NOTE: assumed no mixture volume.
    V_g = n_g * sum(y[i]*mw(i)/rho(T,i,1) for i in 1:5);
    
    // Molar fractions sum to 1.
    sum(x) = 1.0;
    sum(y) = 1.0;
    
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
    
    /****************************
   *** SANITY CHECK SECTION ***
   ****************************/
    
    for i in 1:5 loop
      assert( n[i] > -1e-9, name + ": component "+String(i)+" has a negative quantity.");
    end for;
    
    // Termination criterion: liquid phase goes to zero.
    assert( V_l > -1e-9, "The liquid phase in "+name+" has become negative.");
    assert( V_g > -1e-9, "The gas phase in "+name+" has become negative.");
    
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
</html>"), Icon(Rectangle(extent=[-40,80; 40,-80],  style(
            color=0,
            rgbcolor={0,0,0},
            thickness=2)), Rectangle(extent=[-40,0; 40,-80],  style(
            color=0,
            rgbcolor={0,0,0},
            thickness=2,
            fillColor=4,
            rgbfillColor={0,255,255},
            fillPattern=1))));
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
</html>"), Icon(Rectangle(extent=[-80,40; 80,-20], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=2)), Rectangle(extent=[-80,12; 80,-20], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=2,
            fillColor=67,
            rgbfillColor={85,255,255},
            fillPattern=1))));
  end HorizontalCylindricalTank;
  
  model FuelTank "The methanol fuel tank" 
    import Modelica.Constants.pi;
    extends VerticalCylindricalTank(name="Fuel tank", final m=2, V = 1E-3, A = (42E-3)^2*pi);
    
    CheckPoint bottomFlow = flows[1] annotation (extent=[-10,-90; 10,-70]);
    CheckPoint topFlow = flows[2] annotation (extent=[-10,70; 10,90]);
  equation 
  // TODO add if-equation to compensate in case of empty/full tank.
    flows[1].z_tank = x;
    flows[1].h_tank = h_l;
    flows[2].z_tank = y;
    flows[2].h_tank = h_g;
    
  initial equation 
    //x[1] = 0.2;  // FIXME why does the line below not work? Start bug hunt!
    n[1] / V_l = 1E3;  // 1 M liquid solution
    level =  0.5 * V / A;   // Half-full
    n[4] = 0.0;  // No CO2
    y[5] / 0.79 = y[3] / 0.21;  // N2/O2 in air proportions
    annotation (Icon);
  end FuelTank;
  
  model TestFuelTank "A simple test for the FuelTank class" 
    
    FlowConnector flowConnector annotation (extent=[-32,-24; -2,8]);
    EnvironmentPort environmentPort annotation (extent=[28,-10; 50,16]);
    annotation (Diagram);
    FuelTank fuelTank(T(
                      start = 350),h_k=1000) annotation (extent=[-80,6; -40,44]);
    FlowConnector flowConnector1 annotation (extent=[-26,36; -2,62]);
    
    inner parameter Real RH_env = 50;
    inner parameter Modelica.SIunits.Pressure p_env = 101325;
    inner parameter Modelica.SIunits.Temperature T_env = 298.15;
    
    EnvironmentPort environmentPort1 annotation (extent=[18,48; 48,78]);
  equation 
    
    fuelTank.bottomFlow.F = -1;
    
    connect(flowConnector1.port1, fuelTank.topFlow) annotation (points=[-23.84,
          49; -60,49; -60,40.2], style(pattern=0, thickness=2));
    connect(fuelTank.bottomFlow, flowConnector.port1) annotation (points=[-60,9.8;
          -60,-8; -29.3,-8],      style(pattern=0, thickness=2));
    connect(flowConnector.port2, environmentPort.c) annotation (points=[-4.7,-8;
          18,-8; 18,-8.7; 29.1,-8.7],
                                    style(pattern=0, thickness=2));
    connect(environmentPort1.c, flowConnector1.port2) annotation (points=[19.5,
          49.5; 6,49.5; 6,49; -4.16,49], style(pattern=0, thickness=2));
  end TestFuelTank;
  
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
      flows[i].z_tank = x;
      flows[i].h_tank = h_l;
    end for;
    
    // Last connection is to atmosphere.
    flows[5].z_tank = y;
    flows[5].h_tank = h_g;
    
  initial equation 
    n[1] = 0.1; // (Almost) no methanol in initial solution
    level = 0.5 * V / A; // Half-full
    y[4] = 1e-3; // (Almost) no CO2 in gas
    y[5] / 0.79 = y[3] / 0.21; // N2/O2 in air proportions
    
  end Mixer;
  
  model TestMixer "A Vodka mixer" 
    
    Mixer mixer annotation (extent=[-50,0; -10,42]);
    FuelTank fuelTank annotation (extent=[-88,12; -56,44]);
    annotation (Diagram);
    EnvironmentPort environmentPort annotation (extent=[-50,68; -30,88]);
    EnvironmentPort environmentPort1 annotation (extent=[-6,62; 14,82]);
    FlowConnector flowConnector annotation (extent=[-66,-20; -46,0]);
    FlowConnector flowConnector1 annotation (extent=[-28,42; -8,62]);
    FlowConnector flowConnector2 annotation (extent=[-70,44; -50,64]);
    EnvironmentPort environmentPort2 annotation (extent=[4,-52; 24,-32]);
    Plug plug annotation (extent=[-16,-10; 4,10]);
    Plug plug1 annotation (extent=[-14,-30; 6,-10]);
    FlowConnector flowConnector3 annotation (extent=[-26,-64; -4,-38]);
    inner parameter Real RH_env = 30;
    inner parameter Modelica.SIunits.Pressure p_env = 101325;
    inner parameter Modelica.SIunits.Temperature T_env = 298.15;
    
  equation 
    connect(flowConnector.port2, mixer.flow1) annotation (points=[-47.8,-10; 
          -40,-10; -40,4.2; -36,4.2], style(pattern=0, thickness=2));
    connect(flowConnector.port1, fuelTank.bottomFlow) annotation (points=[-64.2,
          -10; -72,-10; -72,15.2], style(pattern=0, thickness=2));
    connect(flowConnector2.port1, fuelTank.topFlow) annotation (points=[-68.2,
          54; -70,54; -70,40.8; -72,40.8], style(pattern=0, thickness=2));
    connect(flowConnector2.port2, environmentPort.c) annotation (points=[-51.8,
          54; -50,54; -50,69; -49,69], style(pattern=0, thickness=2));
    connect(flowConnector1.port1, mixer.topFlow) annotation (points=[-26.2,52;
          -30,52; -30,37.8], style(pattern=0, thickness=2));
    connect(flowConnector1.port2, environmentPort1.c) annotation (points=[-9.8,52; 
          -6,52; -6,63; -5,63],     style(pattern=0, thickness=2));
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
    connect(flowConnector3.port1, mixer.flow2) annotation (points=[-24.02,-51; 
          -32,-51; -32,4.2], style(
        pattern=0,
        thickness=2,
        fillColor=46,
        rgbfillColor={127,127,0},
        fillPattern=7));
    connect(flowConnector3.port2, environmentPort2.c) annotation (points=[-5.98,
          -51; 5,-51], style(
        pattern=0,
        thickness=2,
        fillColor=46,
        rgbfillColor={127,127,0},
        fillPattern=7));
    flowConnector.port2.F = 1;
    flowConnector3.port2.F = 0.5;
    
  end TestMixer;
  
  model Separator "A gas-liquid separator" 
    import Modelica.Constants.pi;
    extends HorizontalCylindricalTank( name="Separator", final m = 3, V = 0.25E-3, A = (23E-3)^2*pi);
    
    CheckPoint feed = flows[1] "Connection at mid-height in the tank." 
      annotation (extent=[-90,0; -70,20]);
    CheckPoint gasOutlet = flows[2] "Connection at the top of the tank." 
      annotation (extent=[60,30; 80,50]);
    CheckPoint liquidOutlet = flows[3] "Connection on the bottom of the tank." 
      annotation (extent=[60,-30; 80,-10]);
  equation 
    // Set inlet composition in case of backflow: use average composition.
    flows[1].z_tank = z;
    flows[1].h_tank = h_tot;
    
    // Set gas outlet: always gas.
    flows[2].z_tank = y;
    flows[2].h_tank = h_g;
    
    // Set liquid outlet: always liquid.
    flows[3].z_tank = x;
    flows[3].h_tank = h_l;
    
  initial equation 
    n[1] = 0.0;  // No methanol in initial solution
    n[2] = 3;
    n[4] = 0.0;  // No CO2
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
  
  model TestSeparator "Test suite for the separator model." 
    
    inner parameter Real RH_env = 30;
    inner parameter Modelica.SIunits.Pressure p_env = 101325;
    inner parameter Modelica.SIunits.Temperature T_env = 298.15;
    Separator separator annotation (extent=[-16,-30; 34,12]);
    EnvironmentPort environmentPort annotation (extent=[70,22; 90,42]);
    EnvironmentPort environmentPort1 annotation (extent=[54,-38; 74,-18]);
    FlowConnector BottomSeparatorConnector 
                                annotation (extent=[26,-46; 46,-26]);
    FlowConnector TopSeparatorConnector 
                                 annotation (extent=[40,14; 60,34]);
    FlowConnector FuelTankToSeparatorConnector 
                                 annotation (extent=[-38,-16; -18,4]);
    FuelTank fuelTank annotation (extent=[-68,4; -48,24]);
    EnvironmentPort environmentPort2 annotation (extent=[-30,56; -10,76]);
    FlowConnector topFuelTankConnector 
                                 annotation (extent=[-54,46; -34,66]);
  equation 
    
    annotation (Diagram);
    connect(BottomSeparatorConnector.port2, environmentPort1.c) 
                                                     annotation (points=[44.2,
          -36; 50,-36; 50,-37; 55,-37], style(
        pattern=0,
        thickness=2,
        fillPattern=1));
    connect(BottomSeparatorConnector.port1, separator.liquidOutlet) 
                                                         annotation (points=[
          27.8,-36; 26.5,-36; 26.5,-13.2], style(
        pattern=0,
        thickness=2,
        fillPattern=1));
    connect(TopSeparatorConnector.port2, environmentPort.c) 
                                                     annotation (points=[58.2,
          24; 64.6,24; 64.6,23; 71,23], style(
        pattern=0,
        thickness=2,
        fillPattern=1));
    connect(TopSeparatorConnector.port1, separator.gasOutlet) 
                                                       annotation (points=[41.8,24; 
          34,24; 34,-0.6; 26.5,-0.6],     style(
        pattern=0,
        thickness=2,
        fillPattern=1));
    connect(FuelTankToSeparatorConnector.port2, separator.feed) 
                                                  annotation (points=[-19.8,-6; 
          -15.85,-6; -15.85,-6.9; -11,-6.9], style(
        pattern=0,
        thickness=2,
        fillPattern=1));
    connect(fuelTank.bottomFlow, FuelTankToSeparatorConnector.port1) 
                                                       annotation (points=[-58,6; 
          -58,-6; -36.2,-6],    style(
        pattern=0,
        thickness=2,
        fillPattern=1));
    connect(topFuelTankConnector.port2, environmentPort2.c) 
                                                      annotation (points=[-35.8,56; 
          -32,56; -32,57; -29,57],     style(
        pattern=0,
        thickness=2,
        fillPattern=1));
    connect(topFuelTankConnector.port1, fuelTank.topFlow) 
                                                    annotation (points=[-52.2,56; 
          -56,56; -56,22; -58,22],     style(
        pattern=0,
        thickness=2,
        fillPattern=1));
    fuelTank.flows[1].F = -0.1;
    
    der(separator.level) = 0;
  end TestSeparator;
  
  model PipeSegment "A pipe segment, modelled as a tank." 
    extends StirredTank(m=2, name="Pipe segment");
  equation 
    for i in 1:m loop
      flows[i].z_tank = z;
      flows[i].h_tank = h_tot;
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
  
  partial model Pipe "A series of PipeSegments or similar objects" 
    parameter Integer n = 1 "Number of segments";
    replaceable PipeSegment segments[n];
  protected 
    FlowConnector[n-1] connections;
  equation 
    for i in 1:(n-1) loop
      connect( segments[i].flows[2], connections[i].f[1]);
      connect( connections[i].f[2], segments[i+1].flows[1]);
    end for;
    annotation (Documentation(info="<html>
<p>A pipe is a sequence of <tt>PipeSegment</tt> objects, intercalated by <tt>FlowConnector</tt>s.
The array of pipe segments is <tt>redeclarable</tt> so that the class may be used with different
implementations or children classes of <tt>PipeSegment</tt>.</p>
<p>To access the extreme connectors of the pipe, use <tt>pipe.segments[1].flows[1]</tt> and 
<tt>pipe.segments[end].flows[end]</tt>.</p>
</html>"));
  end Pipe;
  
  model HeatExchangerPipeSegment 
    "A section of a pipe passing through a heat exchanger." 
    extends PipeSegment(m=3, name="Pipe section in heat exchanger");
    
  equation 
    // Flow and composition are not important for the third flow.
    flows[3].F = 0;
    flows[3].z = z;
    
    annotation (Documentation(info="<html>
<p>This kind of segments have a third connector through which no material flow
passes, rather only enthalpy flow, which can then be set to a function of
temperature differences.</p>
</html>"));
  end HeatExchangerPipeSegment;
  
  model HeatExchangerPipe "A pipe passing through an heat exchanger." 
    extends Pipe(redeclare HeatExchangerPipeSegment segments[n]);
    
    annotation (Documentation(info="<html>
<p>This particular pipe consists of pipe elements that have a third connector,
which is used to transfer heat with other elements.</p>
</html>"));
  end HeatExchangerPipe;
  
  partial model HeatExchanger "A heat exchanger with two sides" 
    
  annotation (Icon(Ellipse(extent=[-60,60; 60,-60], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=67,
          rgbfillColor={85,255,255},
          fillPattern=1)), Line(points=[-60,0; -20,0; -10,40; 10,-40; 20,0; 60,
            0], style(
          color=0,
          rgbcolor={0,0,0},
          thickness=4,
          fillColor=67,
          rgbfillColor={85,255,255},
          fillPattern=1))), Documentation(info="<html>
<p>This is a simple heat exchanger with two sides (left-right and top-bottom
in the figure) and no phase separation.</p>
</html>"));
    
    import Modelica.SIunits.Area;
    import Modelica.SIunits.CoefficientOfHeatTransfer;
    
    CheckPoint f11 "Fluid 1, side1" 
                                  annotation (extent=[-70,-10; -50,10]);
    CheckPoint f21 "Fluid 2, side 1" 
                                   annotation (extent=[-10,50; 10,70]);
    CheckPoint f22 "Fluid 2, side 2" 
                                   annotation (extent=[-10,-70; 10,-50]);
    CheckPoint f12 "Fluid 1, side 2" 
                                   annotation (extent=[50,-10; 70,10]);
    
    parameter String name="Heat exchanger";
    parameter Area A=1e-2;
    parameter CoefficientOfHeatTransfer U = 10;
    parameter Integer steps = 1;
    
    HeatExchangerPipe side1(n=steps);
    HeatExchangerPipe side2(n=steps);
    
  equation 
    // Connect pipe connectors to the heat exchanger's.
    connect( f11, side1.segments[1].flows[1]);
    connect( f12, side1.segments[end].flows[end]);
    connect( f21, side2.segments[1].flows[1]);
    connect( f22, side2.segments[end].flows[end]);
    
    /* Connect the pipe segments two by two, using the enthalpy element
   * in the third connectors to transport heat.*/
    for i in 1:steps loop
      side1.segments[i].flows[3].H + side2.segments[i].flows[3].H = 0;
      side1.segments[i].flows[3].H = (A/steps) * U * (side1.segments[i].T - side2.segments[i].T);
    end for;
    
  end HeatExchanger;
  
  model FlowController "A unit modelling a pump or a MFC" 
    
  import Thermo.mw;
    
    MolarFlowRate F;
    Modelica.SIunits.MassFlowRate m;
    CheckPoint inlet "Unit inlet" annotation (extent=[-70,-10; -50,10]);
    CheckPoint outlet "Unit outlet" annotation (extent=[50,-10; 70,10]);
  equation 
    m = sum({inlet.z[i] * F * mw(i) for i in 1:5});
    F = inlet.F;
    connect( inlet, outlet);
    
    annotation (Icon(
        Polygon(points=[-60,-80; -44,-40; 44,-40; 60,-80; -60,-80], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=2,
            fillColor=4,
            rgbfillColor={0,255,255},
            fillPattern=1)),
        Ellipse(extent=[-60,60; 60,-60], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=2,
            fillColor=4,
            rgbfillColor={0,255,255},
            fillPattern=1)),
        Line(points=[0,10; 0,-10], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=2,
            fillColor=4,
            rgbfillColor={0,255,255},
            fillPattern=1)),
        Line(points=[-10,0; 10,0], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=2,
            fillColor=4,
            rgbfillColor={0,255,255},
            fillPattern=1))), Documentation(info="<html>
<p>This class allows to set a certain overall molar or mass flow. It is not
immediately possible to set a <em>volume</em> flow, because this would entail
calculating the phase equilibrium, which we are not doing here (though child
classes could specialize).</p>
</html>"));
  end FlowController;
  
  model GasFlowController "A flow controller with only gas phase" 
    extends FlowController;
    annotation (Icon(Text(
          extent=[-10,68; 8,60],
          style(
            color=0,
            rgbcolor={0,0,0},
            thickness=2,
            fillColor=4,
            rgbfillColor={0,255,255},
            fillPattern=1),
          string="Gas")), Documentation(info="<html>
<p>This class implements a mass flow controller with field volumetric units. Since
there are two different standards (the actual \"standard\" at 0� Celsius and the Norm
at 70� Fahrenheit), it is necessary to adjust the reference temperature; the default
assumes zero Celsius (\"standard\" value).</p>
<p>The flow assumes that all components are in gas phase and takes their density from
the Thermo library, where the ideal gas law is (usually) assumed.</p>
</html>"));
    import Thermo.rho;
    import Thermo.mw;
    import Modelica.SIunits.VolumeFlowRate;
    import Modelica.SIunits.Temperature;
    import Modelica.SIunits.Conversions.from_degC;
    import Modelica.SIunits.Conversions.from_degF;
    
    VolumeFlowRate V "Volumetric flow rate.";
    parameter Temperature T = 273.15 "Reference temperature for standard flow.";
  equation 
    V = sum({inlet.z[i] * F * mw(i) / rho(T, i, 1) for i in 1:5});
    
  end GasFlowController;
  
  model Pump "A liquid pump" 
    extends FlowController;
    annotation (Icon(Text(
          extent=[-14,70; 14,60],
          style(
            color=0,
            rgbcolor={0,0,0},
            thickness=2,
            fillColor=4,
            rgbfillColor={0,255,255},
            fillPattern=1),
          string="Liquid")), Documentation(info="<html>
<p>This class implements a liquid pump with field volumetric units. It will be
necessary to somehow specify its operating temperature to make it work: for
example, depending on flow direction, it might be the temperature of the element 
upstream or the one downstream.</p>
<p>The flow assumes that only water and methanol are present and are completely
in liquid phase; it takes their density from the Thermo library.</p>
</html>"));
    import Thermo.rho;
    import Thermo.mw;
    import Modelica.SIunits.VolumeFlowRate;
    import Modelica.SIunits.Temperature;
    import Modelica.SIunits.Conversions.from_degC;
    import Modelica.SIunits.Conversions.from_degF;
    
    VolumeFlowRate V "Standard volume flow rate.";
    Temperature T "Temperature (not a reference!)."; // TODO couple somehow to entering composition/enthalpy?
    
  equation 
      V = sum({inlet.z[i] * F * mw(i) / rho(T, i, 2) for i in 1:2});
    
  end Pump;
end Tank;
