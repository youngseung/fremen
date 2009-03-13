

package Tank 
  
type MolarFlowRate = Real(final quantity="MolarFlowRate", final unit="mol/s") 
    annotation (Documentation(info="<html>
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
</html>"), Icon(Ellipse(extent=[-100,100; 100,-100],
                                                 style(
            pattern=0,
            thickness=2,
            gradient=3,
            fillColor=1,
            rgbfillColor={255,0,0}))));
  end CheckPoint;
  
  model Plug "A class that blocks a flow connection" 
    import Thermo.AllSpecies;
    annotation (Icon(Polygon(points=[-80,60; -80,-60; 100,60; 100,-60; -80,60],
                                                                              style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4,
            fillColor=7,
            rgbfillColor={255,255,255}))), Documentation(info="<html>
<p>This simple item sets the connected flow to zero.</p>
</html>"));
    CheckPoint c annotation (extent=[-100,-10; -80,10]);
  protected 
    constant Integer k = size(AllSpecies, 1);
  equation 
    c.F = 0;
    c.z = zeros(k);
    c.H = 0;
    
  end Plug;
  
  model FlowConnector "A class to determine flow direction" 
    
  public 
    CheckPoint port1 annotation (extent=[-82,-10; -62,10]);
    CheckPoint port2 annotation (extent=[62,-10; 82,10]);
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
        Polygon(points=[-64,0; -20,-100; -20,100; -64,0], style(
            pattern=0,
            fillColor=0,
            rgbfillColor={0,0,0})),
        Polygon(points=[64,0; 20,-100; 20,100; 64,0], style(
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
    import Thermo.Methanol;
    import Thermo.Water;
    import Thermo.LiquidPhase;
    import Modelica.SIunits.Temperature;
    import Modelica.SIunits.Concentration;
    import Modelica.SIunits.MoleFraction;
    
    outer Temperature T_env "Environment temperature.";
    
    parameter Concentration C = 1000 
      "Concentration of methanol in water, mol/m.";
    parameter Temperature T = T_env "Source temperature.";
    
    MoleFraction x_ch3oh "Molar fraction of methanol.";
    MoleFraction x_h2o "Molar fraction of water.";
    
    CheckPoint p "Connection point of the source" 
      annotation (extent=[-20,-20; 20,20]);
  equation 
    assert( C >= 0, "==> Negative concentration given in MethanolSolution object.");
    assert( C <= rho(T,Methanol,LiquidPhase)/mw(Methanol), "==> Methanol concentration over limit (" + String(mw(Methanol)/rho(T,Methanol,LiquidPhase)) + " mol/m).");
    
    C = x_ch3oh / ( x_ch3oh*mw(Methanol)/rho(T,Methanol,LiquidPhase) + x_h2o*mw(Water)/rho(T,Water,LiquidPhase));
    x_ch3oh + x_h2o = 1.0;
    
    p.z_local = {x_ch3oh, x_h2o, 0, 0, 0};
    p.h_local = x_ch3oh*h(T,Methanol,LiquidPhase) + x_h2o*h(T,Water,LiquidPhase);
    
    annotation (Icon(Ellipse(extent=[-100,100; 100,-100],
                                                      style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4,
            fillColor=7,
            rgbfillColor={255,255,255}))), Documentation(info="<html>
<p>This item is a source for methanol-water solutions. Parameter <tt>C</tt>
allows to set the concentration in moler per <em>cubic metre</em>; note that
this is 1000 times the normal scale (1M = 1000 mol/m).</p>
</html>"));
  end MethanolSolution;
  
  model EnvironmentPort "A flow connection to environment conditions." 
    import Thermo.Water;
    import Thermo.Oxygen;
    import Thermo.Nitrogen;
    import Thermo.GasPhase;
    import Thermo.h;
    import Thermo.K;
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
    z_h2o = RH_env/100 * K(T_env, Water);
    
    h_h2o = h(T_env, Water, GasPhase);
    h_o2  = h(T_env, Oxygen, GasPhase);
    h_n2  = h(T_env, Nitrogen, GasPhase);
    
    h_air = h_h2o*z_h2o + h_o2*z_o2 + h_n2*z_n2;
    
    c.z_local = {0.0, z_h2o, z_o2, 0.0, z_n2};
    c.h_local = h_air;
    
    annotation (Documentation(info="<html>
<p>This dummy object is used whenever a system flow goes to the environment, e.g.
the gas outlet of separators.</p>
<p>If fluid is drawn from this object, it will have the composition of air with
the relative humidity set in the outer variable <tt>RH_env</tt>, which can have
values from 0 (dry air) to 100 (saturated air).</p>
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
  
  partial model Equilibrium "a class describing liquid-vapour equilibria." 
    
    import Modelica.SIunits.AmountOfSubstance;
    import Modelica.SIunits.Temperature;
    import Modelica.SIunits.MoleFraction;
    
    import Thermo.K;
    import Thermo.LiquidSpecies;
    import Thermo.GasSpecies;
    import Thermo.AllSpecies;
    import Thermo.Methanol;
    import Thermo.Water;
    import Thermo.rr;
    import Thermo.speciesName;
    
    Temperature T "Representative temperature.";
    MoleFraction beta;
    
    MoleFraction[size(AllSpecies, 1)] z(each min=0, each max=1) 
      "Overall molar fraction.";
    MoleFraction[size(AllSpecies, 1)] x(each min=0, each max=1) 
      "Liquid molar fraction.";
    MoleFraction[size(AllSpecies, 1)] y(each min=0, each max=1) 
      "Gaseous molar fraction.";
    
  equation 
    // Guard: if K_water == 1, Thermo.rr will return a singularity.
    if K(T,Water) >= 1 or rr(z[Methanol],z[Water],T) > 1 then
      beta = 1;
    else
      beta = Thermo.rr(z[Methanol],z[Water],T);
    end if;
    
    // Component material balance.
    z = beta*y + (1-beta)*x;
    
    /* Notes on mole-fraction consistency:
   * We assume that z is given and sums to 1. Then, it follows that:
   * 1) sum(x) = 1 is linearly dependent with sum(z) = 1, the Rachford-Rice relation and material balance.
   * 2) sum(y) = 1 is also linearly dependent, since we are using the Rachford-Rice relation.
   * Therefore, there is no particular need to write the consistency explicitly. */
      y[LiquidSpecies] = {K(T, i) * x[i] for i in LiquidSpecies};
      x[GasSpecies]    = zeros(size(GasSpecies,1));
    
    annotation (Documentation(info="<html>
<p>This class allows to calculate the phase equilibrium of components in
the <tt>Thermo</tt> library, and has n<sub>c</sub> unsaturated degrees of 
freedom, for example temperature and one complete composition (since 
compositions sum to one, only n<sub>c</sub>-1 variables are needed to set
a composition).</p>
<p>In order to calculate the physically correct values for the compositions,
one has to solve the Rachford-Rice equation as indicated by Whitson &
Michelsen (1989), namely finding the root of h(&beta;) in the interval
(beta<sub>min</sub>, beta<sub>max</sub>). If the resulting &beta; is too large
(>1) it will have to be capped; if the conditions are outside the validity
range of Rachford-Rice (typically, temperature is above water's boiling point)
&beta; must be set forcefully.</p>
</html>"), Icon,
      DymolaStoredErrors);
  end Equilibrium;
  
  partial model ExtensiveBalances 
    "Equations for the balance of extensive properties" 
    
    import Modelica.SIunits.AmountOfSubstance;
    import Modelica.SIunits.InternalEnergy;
    import Modelica.SIunits.Heat;
    import Thermo.AllSpecies;
    
    parameter Integer m(min=1) "Number of flows.";
    CheckPoint[m] flows "Connections with other objects.";
    
    AmountOfSubstance[size(AllSpecies,1)] n(each min=0) 
      "Accumulated moles of each species.";
    InternalEnergy U "The accumulated internal energy.";
    Heat Q "Heat exchanged with environment.";
    
  equation 
    der(n) = {sum(flows[j].F * flows[j].z[i] for j in 1:m) for i in AllSpecies};
    der(U) = sum(flows[j].H for j in 1:m) + Q;
    
  end ExtensiveBalances;
  
  partial model StirredTank "A generic stirred tank with an undefined shape." 
    extends Equilibrium(T(start=T_env,fixed=true));
    extends ExtensiveBalances(m=1);
    
    import Modelica.SIunits.AmountOfSubstance;
    import Modelica.SIunits.InternalEnergy;
    import Modelica.SIunits.Pressure;
    import Modelica.SIunits.Temperature;
    import Modelica.SIunits.HeatCapacity;
    import Modelica.SIunits.CoefficientOfHeatTransfer;
    import Modelica.SIunits.Area;
    import Modelica.SIunits.Volume;
    
    import Thermo.MolarEnthalpy;
    import Thermo.rho;
    import Thermo.h;
    import Thermo.mw;
    import Thermo.AllSpecies;
    import Thermo.LiquidSpecies;
    import Thermo.GasSpecies;
    import Thermo.LiquidPhase;
    import Thermo.GasPhase;
    
    outer parameter Temperature T_env = 298.15 "Environment temperature";
    
    parameter String name = "Unnamed stirred tank" "Tank identifier.";
    parameter HeatCapacity Cp = 100 "Heat capacity of tank (glass, lid, ...).";
    parameter CoefficientOfHeatTransfer k_h = 0.0 
      "Default: perfect insulation.";
    parameter Area A_sur = 1e-3 
      "Surface area of contact between tank and environment";
    parameter Volume V = 1E-3 "Total volume of the tank.";
    
    AmountOfSubstance n_tot(min=0) = sum(n) "Total moles.";
    
    AmountOfSubstance n_g_tot(min=0) "Total moles of gas.";
    AmountOfSubstance n_l_tot(min=0) "Total moles of liquid.";
    
    AmountOfSubstance[size(AllSpecies,1)] n_l = x*n_l_tot;
    AmountOfSubstance[size(AllSpecies,1)] n_g = y*n_g_tot;
    
    MolarEnthalpy h_tot = (h_g*n_g_tot + h_l*n_l_tot)/n_tot 
      "Overall molar enthalpy in the tank.";
    MolarEnthalpy h_l = sum(h(T,i,LiquidPhase)*x[i] for i in LiquidSpecies) 
      "The molar enthalpy of the solution in liquid phase.";
    MolarEnthalpy h_g = sum(h(T,i,GasPhase)*y[i] for i in AllSpecies) 
      "The molar enthalpy of the mixture in gas phase.";
    
    InternalEnergy tankSensibleHeat = Cp*(T-T_env) 
      "The heat accumulated in the tank's solid parts, such as glass, lid, etc.";
    InternalEnergy gasSpeciesEnergy = sum(h(T,i,GasPhase)*n_g_tot*y[i] for i in AllSpecies) 
      "The internal energy of the gaseous species.";
    InternalEnergy liquidSpeciesEnergy = sum(h(T,i,LiquidPhase)*n_l_tot*x[i] for i in LiquidSpecies) 
      "The internal energy of the liquid species, including latent heat.";
    
    Volume V_l(min=0) = sum(n_l[i]*mw(i)/rho(T,i,LiquidPhase) for i in LiquidSpecies) 
      "Amount of liquid volume.";
    Volume V_g(min=0) = sum(n_g[i]*mw(i)/rho(T,i,GasPhase) for i in AllSpecies) 
      "Amount of gaseous volume.";
    
  equation 
    if beta >= 1 or Thermo.K(T, Thermo.Water) > 1 then
      n_g_tot = n_tot;
      n_l_tot = 0;
    elseif beta <= 0 then
      n_g_tot = 0;
      n_l_tot = n_tot;
    else
      n_g_tot = beta*n_tot;
      n_l_tot = (1-beta)*n_tot;
    end if;
    
    // Exchange of heat with the environment.
    Q = A_sur * k_h * (T_env - T);
    
    // Allocation of internal energy. See each term for a description.
    U = tankSensibleHeat + gasSpeciesEnergy + liquidSpeciesEnergy;
    
    // Sum of gas and liquid volumes is constant and equal to the total volume.
    V = V_l + V_g;
    
    n = z * sum(n);
    
    annotation (Documentation(info="<html>
<p>This class inherits from the <tt>Equilibrium</tt> and <tt>ExtensiveBalances</tt>
classes, and combines them adding the constraint of a fixed volume.</p>
<p>The number of flows <em>m</em> cannot be set to zero; if you need such a
isolated tank, set it to 1 and set <tt>flow[1].H = 0</tt>, <tt>flow[1].F = 0</tt>.</p>
</html>"), Icon,
      DymolaStoredErrors);
  end StirredTank;
  
  partial model VerticalCylindricalTank 
    "A cylindrical tank in horizontal position." 
    extends StirredTank(name="Vertical cylindrical tank",V_l(start=V_l_0,fixed=true));
    import Modelica.SIunits.Area;
    import Modelica.SIunits.Length;
    import Modelica.SIunits.Volume;
    
    parameter Area A "The cross-section of the tank.";
    parameter Volume V_l_0 = 5e-4 "Initial solution volume";
    Length level(min=0) "The liquid level, measured from the bottom.";
    
  equation 
    level = V_l / A;
    assert( V_l_0 >= 0 and V_l_0 <= V, "==> Bad starting solution volume in "+name+": V_0 = "+String(V_l_0)+".");
    assert( V_l >= 0, "==> "+name+" ran out of liquid.");
    assert( V_l <= V, "==> "+name+" overflowed.");
    
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
    
    import Modelica.SIunits.Volume;
    import Modelica.SIunits.Concentration;
    import Thermo.Methanol;
    import Thermo.Water;
    import Thermo.Oxygen;
    import Thermo.CarbonDioxide;
    import Thermo.Nitrogen;
    
    CheckPoint bottomFlow = flows[1] annotation (extent=[-10,-90; 10,-70]);
    CheckPoint topFlow = flows[2] annotation (extent=[-10,70; 10,90]);
  equation 
    flows[1].z_local = x;
    flows[1].h_local = h_l;
    flows[2].z_local = y;
    flows[2].h_local = h_g;
    
  initial equation 
    x[Methanol] = 1;
    n[CarbonDioxide] = 0.0;  // No CO2, anywhere
    n[Nitrogen] / 79 = n[Oxygen] / 21;  // N2/O2 in air proportions
    annotation (Icon);
  end FuelTank;
  
  model Mixer "Our laboratory mixer" 
    import Modelica.Constants.pi;
    extends VerticalCylindricalTank(name="Mixer", final m=5, V=5E-4, A=(37E-3)^2*pi, V_l_0=2.5E-4);
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
</html>"),   Icon,
      Diagram);
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
    
    import Modelica.SIunits.Concentration;
    import Modelica.SIunits.Length;
    import Thermo.Methanol;
    import Thermo.Water;
    import Thermo.Oxygen;
    import Thermo.CarbonDioxide;
    import Thermo.Nitrogen;
    
    parameter Concentration c_0 = 1000 "Initial concentration";
    
    Concentration c(start=c_0,fixed=true) = n[Methanol]/V_l 
      "Initial concentration";
    
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
    y[CarbonDioxide] = 0.0; // No CO2
    y[Nitrogen] / 0.79 = y[Oxygen] / 0.21; // N2/O2 in air proportions
    
  end Mixer;
  
  model Separator "A gas-liquid separator" 
    import Modelica.Constants.pi;
    extends HorizontalCylindricalTank( name="Separator", final m = 3, V = 0.25E-3, A = (23E-3)^2*pi);
    
    CheckPoint feed = flows[1] "Connection at mid-height in the tank." 
      annotation (extent=[-100,-10; -80,10]);
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
    CheckPoint side1 
      annotation (extent=[-100,-10; -80,10]);
    CheckPoint side2                        annotation (extent=[80,-10; 100,10]);
  protected 
    FlowConnector[n-1] connections;
  equation 
    connect(side1, segments[1].flows[1]);
    connect(side2, segments[end].flows[2]);
    for i in 1:(n-1) loop
      connect(segments[i].flows[2], connections[i].port1);
      connect(connections[i].port2, segments[i+1].flows[1]);
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
    
  annotation (Icon(Ellipse(extent=[-80,80; 80,-80], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=2,
            fillColor=7,
            rgbfillColor={255,255,255})),
                           Line(points=[-80,0; -40,0; 0,58; 0,-60; 40,0; 80,0],
                style(
          color=0,
          rgbcolor={0,0,0},
          thickness=4,
          fillColor=67,
          rgbfillColor={85,255,255},
          fillPattern=1)),
        Text(
          extent=[10,100; 14,92],
          string="1",
          style(color=0, rgbcolor={0,0,0})),
        Text(
          extent=[-96,16; -92,8],
          string="1",
          style(color=0, rgbcolor={0,0,0})),
        Text(
          extent=[92,-10; 96,-18],
          style(color=0, rgbcolor={0,0,0}),
          string="2"),
        Text(
          extent=[-14,-90; -10,-98],
          style(color=0, rgbcolor={0,0,0}),
          string="2")),     Documentation(info="<html>
<p>This is a simple heat exchanger with two sides (left-right and top-bottom
in the figure) and no phase separation.</p>
<p>Both sides are assumed to have either liquid phase or two-phase. For gas
phase, see class <tt>Cooler</tt>.</p>
</html>"));
    
    import Modelica.SIunits.Area;
    import Modelica.SIunits.Volume;
    import Modelica.SIunits.HeatFlowRate;
    import Modelica.SIunits.CoefficientOfHeatTransfer;
    
    CheckPoint f11 "Fluid 1, side1" 
                                  annotation (extent=[-100,-10; -80,10]);
    CheckPoint f21 "Fluid 2, side 1" 
                                   annotation (extent=[-10,80; 10,100]);
    CheckPoint f22 "Fluid 2, side 2" 
                                   annotation (extent=[-10,-100; 10,-80]);
    CheckPoint f12 "Fluid 1, side 2" 
                                   annotation (extent=[80,-10; 100,10]);
    
    parameter String name="Heat exchanger" "Name identifying the unit.";
    parameter Area A=1e-2 "The total heat transfer area of the exchanger.";
    parameter CoefficientOfHeatTransfer U = 10 "The heat-transfer coefficient.";
    parameter Integer steps = 5 
      "The number of subvolumes in which the two sides are divided.";
    parameter Volume V_1 = 1e-4 "The total volume of the first side.";
    parameter Volume V_2 = 1e-4 "The total volume of the second side.";
    
    HeatFlowRate Q "The heat moving from the first to the second side.";
    
    replaceable HeatExchangerPipe side1(n=steps, segments(each V=V_1/steps));
    replaceable HeatExchangerPipe side2(n=steps, segments(each V=V_2/steps));
    
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
  
  model Cooler "A heat exchanger with only gas phase on side 2" 
    
    annotation (Icon(
        Text(
          extent=[-14,56; -10,48],
          string="1",
          style(color=0, rgbcolor={0,0,0})),
        Text(
          extent=[-92,16; -88,8],
          string="1",
          style(color=0, rgbcolor={0,0,0})),
        Text(
          extent=[88,16; 92,8],
          style(color=0, rgbcolor={0,0,0}),
          string="2"),
        Text(
          extent=[-14,-48; -10,-56],
          style(color=0, rgbcolor={0,0,0}),
          string="2"),
        Rectangle(extent=[-80,40; 80,-40], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4)),
        Ellipse(extent=[-60,20; -20,-20], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4,
            fillColor=7,
            rgbfillColor={255,255,255})),
        Ellipse(extent=[20,20; 60,-20], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4,
            fillColor=7,
            rgbfillColor={255,255,255})),
        Polygon(points=[-28,16; -28,-16; 28,16; 28,-16; -28,16], style(
            pattern=0,
            fillColor=7,
            rgbfillColor={255,255,255})),
        Line(points=[-28,16; 28,-16], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4,
            fillColor=7,
            rgbfillColor={255,255,255},
            fillPattern=8)),
        Line(points=[28,16; -28,-16], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4,
            fillColor=7,
            rgbfillColor={255,255,255},
            fillPattern=8))),
                            Documentation(info="<html>
<p>This is a simple heat exchanger with two sides (left-right and top-bottom
in the figure) and no phase separation.</p>
<p>Both sides are assumed to have either liquid phase or two-phase. For gas
phase, see class <tt>Cooler</tt>.</p>
</html>"));
    import Modelica.SIunits.Area;
    import Modelica.SIunits.Volume;
    import Modelica.SIunits.HeatFlowRate;
    import Modelica.SIunits.CoefficientOfHeatTransfer;
    
    CheckPoint f11 "Fluid 1, side1" 
                                  annotation (extent=[-100,-10; -80,10]);
    CheckPoint f21 "Fluid 2, side 1" 
                                   annotation (extent=[-10,40; 10,60]);
    CheckPoint f22 "Fluid 2, side 2" 
                                   annotation (extent=[-10,-60; 10,-40]);
    CheckPoint f12 "Fluid 1, side 2" 
                                   annotation (extent=[80,-10; 100,10]);
    
    parameter String name="Heat exchanger" "Name identifying the unit.";
    parameter Area A=1e-2 "The total heat transfer area of the cooler.";
    parameter CoefficientOfHeatTransfer U = 1000 
      "The heat-transfer coefficient.";
    parameter Integer steps = 5 
      "The number of subvolumes in which the two sides are divided.";
    parameter Volume V_1 = 1e-4 "The total volume of the first side.";
    parameter Volume V_2 = 1e-4 "The total volume of the second side.";
    
    HeatFlowRate Q "The heat moving from the first to the second side.";
    
    replaceable HeatExchangerPipe side1(n=steps, segments(each V=V_1/steps));
    replaceable HeatExchangerGasPipe side2(n=steps, segments(each V=V_2/steps));
    
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
    
  end Cooler;
  
  model FlowController "A unit modelling a pump or a MFC" 
    
    import Thermo.mw;
    import Thermo.AllSpecies;
    
    MolarFlowRate F;
    Modelica.SIunits.MassFlowRate m;
    CheckPoint inlet "Unit inlet" annotation (extent=[-12,-10; 8,10]);
    CheckPoint outlet "Unit outlet" annotation (extent=[-10,90; 10,110]);
  equation 
    m = sum({inlet.z[i] * F * mw(i) for i in AllSpecies});
    F = inlet.F;
    
    // Connect the two CheckPoints.
    inlet.F + outlet.F = 0;
    inlet.H + outlet.H = 0;
    inlet.z = outlet.z;
    inlet.z_local = outlet.z_local;
    inlet.h_local = outlet.h_local;
    
    annotation (Icon(
        Ellipse(extent=[-100,100; 100,-100],
                                         style(
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
</html>"),
      Diagram);
  end FlowController;
  
  model GasFlowController "A flow controller with only gas phase" 
    extends FlowController;
    annotation (Icon,     Documentation(info="<html>
<p>This class implements a mass flow controller with field volumetric units. Since
there are two different standards (the actual \"standard\" at 0 Celsius and the Norm
at 70 Fahrenheit), it is necessary to adjust the reference temperature; the default
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
  
  model Pump "A pump, with only liquid phase." 
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
    import Thermo.h;
    import Thermo.LiquidPhase;
    import Thermo.LiquidSpecies;
    import Modelica.SIunits.VolumeFlowRate;
    import Modelica.SIunits.Temperature;
    import Modelica.SIunits.Conversions.from_degC;
    import Modelica.SIunits.Conversions.from_degF;
    
    outer parameter Temperature T_env;
    Temperature T "Temperature of the passing flow.";
    VolumeFlowRate V "Standard volume flow rate.";
    
  equation 
    inlet.H = inlet.F * (sum(inlet.z[i]*h(T, i, LiquidPhase) for i in LiquidSpecies));
    V = sum({inlet.z[i] * F * mw(i) / rho(T, i, LiquidPhase) for i in LiquidSpecies});
    
  end Pump;
  
  package Tests "A series of tests for units defined in the Tank library" 
    
    model TestEquilibrium "Test for the Equilibrium class" 
      
      parameter Modelica.SIunits.Temperature T_start = 273.15;
      parameter Modelica.SIunits.MoleFraction ch3oh = 0.2;
      parameter Modelica.SIunits.MoleFraction h2o = 0.4;
      
      model MyEquilibrium 
        extends Equilibrium;
      end MyEquilibrium;
      
      MyEquilibrium eq;
      
    equation 
      eq.T = T_start+time;
      
      eq.z={ch3oh, h2o,0,0,0};
      
      annotation (experiment(StopTime=105), experimentSetupOutput(
            doublePrecision=true, derivatives=false));
    end TestEquilibrium;
    
    model TestExtensiveBalances 
      extends ExtensiveBalances(m=1);
      
    equation 
      for i in 1:m loop
        flows[i].h_local = 0;
        flows[i].z_local = 0*flows[i].z_local;
        flows[i].z = ones(size(flows[i].z,1))/size(flows[i].z,1);
        flows[i].H = 1;
        flows[i].F = 1;
      end for;
      Q = -0.5;
      
    end TestExtensiveBalances;
    
    model TestStirredtank "Test for the generic stirred tank" 
      
      inner parameter Real RH_env = 60;
      inner parameter Modelica.SIunits.Pressure p_env = 101325;
      inner parameter Modelica.SIunits.Temperature T_env = 298.15;
      
      import Thermo.Methanol;
      import Thermo.Water;
      import Thermo.Oxygen;
      import Thermo.CarbonDioxide;
      import Thermo.Nitrogen;
      
      model MyStirredTank 
        extends StirredTank(m=2);
      end MyStirredTank;
      
      parameter Modelica.SIunits.MoleFraction z_water_0 = 0.5 
        "Initial water fraction";
      
      MyStirredTank tank;
    equation 
      for i in 1:tank.m loop
        tank.flows[i].z_local = tank.z;
        tank.flows[i].h_local = tank.h_tot;
      end for;
      
      tank.flows[1].F=1;
      tank.flows[1].z={0.1,0.2,0,0,0.7};
      tank.flows[1].H=0;
      
      tank.flows[2].z=tank.z;
      tank.flows[2].H=tank.flows[2].F*tank.h_tot;
      
    initial equation 
      tank.z[Methanol] = 0;
      tank.z[Water] = z_water_0;
      tank.z[CarbonDioxide] = 0;
      tank.z[Oxygen] / 0.21 = tank.z[Nitrogen] / 0.79;
      annotation (experiment(Algorithm="Dassl"), experimentSetupOutput(
            doublePrecision=true));
    end TestStirredtank;
    
    model TestFuelTank "A simple test for the FuelTank class" 
      
      FlowConnector bottomConnector 
                                  annotation (extent=[-28,-24; -12,-16]);
      EnvironmentPort environmentPort annotation (extent=[18,-28; 48,4]);
      annotation (Diagram,
        experiment(Algorithm="Dassl"),
        experimentSetupOutput(doublePrecision=true, derivatives=false));
      FuelTank fuelTank                        annotation (extent=[-80,6; -40,44]);
      FlowConnector topConnector   annotation (extent=[-26,56; -12,64]);
      
      inner parameter Real RH_env = 50;
      inner parameter Modelica.SIunits.Pressure p_env = 101325;
      inner parameter Modelica.SIunits.Temperature T_env = 298.15;
      
      EnvironmentPort environmentPort1 annotation (extent=[18,52; 48,84]);
    equation 
      
      fuelTank.bottomFlow.F = -1;
      
      connect(topConnector.port1, fuelTank.topFlow)   annotation (points=[-24.04,
            60; -60,60; -60,40.2], style(pattern=0, thickness=2));
      connect(fuelTank.bottomFlow, bottomConnector.port1) 
                                                        annotation (points=[-60,9.8;
            -60,-20; -25.76,-20],   style(pattern=0, thickness=2));
      connect(bottomConnector.port2, environmentPort.c) 
                                                      annotation (points=[-14.24,
            -20; 19.5,-20],           style(pattern=0, thickness=2));
      connect(environmentPort1.c, topConnector.port2)   annotation (points=[19.5,60;
            -13.96,60],                    style(pattern=0, thickness=2));
    end TestFuelTank;
    
    model TestMixer "A Vodka mixer" 
      
      Mixer mixer annotation (extent=[-50,0; -10,42]);
      annotation (Diagram,
        experiment,
        experimentSetupOutput(doublePrecision=true));
      EnvironmentPort environmentPort1 annotation (extent=[-2,48; 16,64]);
      FlowConnector flowConnector annotation (extent=[-60,-12; -46,0]);
      FlowConnector flowConnector1 annotation (extent=[-26,46; -12,58]);
      EnvironmentPort environmentPort2 annotation (extent=[4,-56; 24,-36]);
      Plug plug annotation (extent=[-16,-10; 4,10]);
      Plug plug1 annotation (extent=[-14,-30; 6,-10]);
      FlowConnector flowConnector3 annotation (extent=[-22,-58; -8,-44]);
      inner parameter Real RH_env = 30;
      inner parameter Modelica.SIunits.Pressure p_env = 101325;
      inner parameter Modelica.SIunits.Temperature T_env = 298.15;
      
      MethanolSolution source annotation (extent=[-86,8; -66,28]);
    equation 
      mixer.flow1.F = 1;
      mixer.flow2.F = 0.5;
      
      connect(flowConnector.port2, mixer.flow1) annotation (points=[-47.96,-6;
            -40,-6; -40,4.2; -36,4.2],  style(pattern=0, thickness=2));
      connect(flowConnector1.port1, mixer.topFlow) annotation (points=[-24.04,
            52; -30,52; -30,37.8],
                               style(pattern=0, thickness=2));
      connect(flowConnector1.port2, environmentPort1.c) annotation (points=[-13.96,
            52; -1.1,52],             style(pattern=0, thickness=2));
      connect(plug.c, mixer.flow4) annotation (points=[-15,6.10623e-16; -24,
            6.10623e-16; -24,4.2],
          style(
          pattern=0,
          thickness=2,
          fillColor=46,
          rgbfillColor={127,127,0},
          fillPattern=7));
      connect(plug1.c, mixer.flow3) annotation (points=[-13,-20; -28,-20; -28,
            4.2],
          style(
          pattern=0,
          thickness=2,
          fillColor=46,
          rgbfillColor={127,127,0},
          fillPattern=7));
      connect(flowConnector3.port1, mixer.flow2) annotation (points=[-20.04,-51;
            -32,-51; -32,4.2], style(
          pattern=0,
          thickness=2,
          fillColor=46,
          rgbfillColor={127,127,0},
          fillPattern=7));
      connect(flowConnector3.port2, environmentPort2.c) annotation (points=[-9.96,
            -51; 5,-51], style(
          pattern=0,
          thickness=2,
          fillColor=46,
          rgbfillColor={127,127,0},
          fillPattern=7));
      connect(source.p, flowConnector.port1) annotation (points=[-76,18; -76,-6;
            -58.04,-6],      style(pattern=0, thickness=2));
    end TestMixer;
    
    model TestSeparator "Test suite for the separator model." 
      
      inner parameter Real RH_env = 30;
      inner parameter Modelica.SIunits.Pressure p_env = 101325;
      inner parameter Modelica.SIunits.Temperature T_env = 298.15;
      Separator separator annotation (extent=[-36,-24; 16,16]);
      EnvironmentPort environmentPort annotation (extent=[48,20; 72,44]);
      EnvironmentPort environmentPort1 annotation (extent=[46,-40; 68,-16]);
      FlowConnector BottomSeparatorConnector 
                                  annotation (extent=[24,-38; 38,-30]);
      FlowConnector TopSeparatorConnector 
                                   annotation (extent=[24,22; 38,30]);
      FlowConnector FuelTankToSeparatorConnector 
                                   annotation (extent=[-56,-8; -42,0]);
      MethanolSolution source(T=350) annotation (extent=[-72,4; -62,14]);
      annotation (Diagram);
      
    equation 
      source.p.F = -1;
      der(separator.V_l) = 0;
      
      connect(BottomSeparatorConnector.port2, environmentPort1.c) 
                                                       annotation (points=[36.04,
            -34; 47.1,-34],               style(
          pattern=0,
          thickness=2,
          fillPattern=1));
      connect(BottomSeparatorConnector.port1, separator.liquidOutlet) 
                                                           annotation (points=[25.96,
            -34; 8.2,-34; 8.2,-10],          style(
          pattern=0,
          thickness=2,
          fillPattern=1));
      connect(TopSeparatorConnector.port2, environmentPort.c) 
                                                       annotation (points=[36.04,26;
            49.2,26],                     style(
          pattern=0,
          thickness=2,
          fillPattern=1));
      connect(TopSeparatorConnector.port1, separator.gasOutlet) 
                                                         annotation (points=[25.96,26;
            8,26; 8,2; 8.2,2],              style(
          pattern=0,
          thickness=2,
          fillPattern=1));
      connect(FuelTankToSeparatorConnector.port2, separator.feed) 
                                                    annotation (points=[-43.96,
            -4; -38.68,-4; -38.68,-4; -33.4,-4],
                                               style(
          pattern=0,
          thickness=2,
          fillPattern=1));
      connect(source.p, FuelTankToSeparatorConnector.port1) annotation (points=[-67,9;
            -67,-4; -54.04,-4],    style(pattern=0, thickness=2));
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
      connect(inlet.port2, pipe.side1);
      connect(pipe.side2, outlet.port1);
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
    
    model TestGasPipe 
      inner parameter Real RH_env = 30;
      inner parameter Modelica.SIunits.Pressure p_env = 101325;
      inner parameter Modelica.SIunits.Temperature T_env = 298.15;
      GasPipe pipe(n=10, segments(T(each start = 340)));
      EnvironmentPort nature_in;
      FlowConnector inlet;
      EnvironmentPort nature_out;
      FlowConnector outlet;
    equation 
      connect(nature_in.c, inlet.port1);
      connect(inlet.port2, pipe.segments[1].flows[1]);
      connect(pipe.segments[end].flows[end], outlet.port1);
      connect(outlet.port2, nature_out.c);
      
      nature_in.c.F = -1;
    initial equation 
      for i in 1:pipe.n loop
        pipe.segments[i].z[1]=0;
        pipe.segments[i].z[2]=0;
        pipe.segments[i].z[4]=0;
        pipe.segments[i].z[3]/21 = pipe.segments[i].z[5]/79;
      end for;
      
    end TestGasPipe;
    
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
      connect(to_f11.port2, hx.f11)                   annotation (points=[-24.8,10;
            -19.8,10; -19.8,10; -14.6,10],   style(pattern=0, thickness=2));
      connect(hx.f12, to_f12.port1)                    annotation (points=[10.6,10;
            11.8,10; 11.8,10; 14.8,10],     style(pattern=0, thickness=2));
      connect(to_f12.port2, solutionOutlet.c) 
        annotation (points=[29.2,10; 36.8,10], style(pattern=0, thickness=2));
      connect(to_f22.port2, coolingInlet.c)     annotation (points=[19.2,-22;
            26.9,-22], style(pattern=0, thickness=2));
      connect(methanolSolution.p, to_f11.port1) 
        annotation (points=[-54,10; -39.2,10], style(pattern=0, thickness=2));
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
      connect(to_f22.port1, hx.f22) annotation (points=[4.8,-22; -2,-22; -2,
            -2.6],
          style(pattern=0, thickness=2));
      connect(to_f21.port2, coolingOutlet.c) 
        annotation (points=[17.2,46; 30.8,46], style(pattern=0, thickness=2));
      connect(to_f21.port1, hx.f21) annotation (points=[2.8,46; -2,46; -2,22.6],
          style(pattern=0, thickness=2));
    end TestHeatExchanger;
    
    model TestCooler 
      
      inner parameter Real RH_env = 30;
      inner parameter Modelica.SIunits.Pressure p_env = 101325;
      inner parameter Modelica.SIunits.Temperature T_env = 298.15;
      
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
      Cooler hx annotation (extent=[-20,-6; 10,26]);
    equation 
      connect(to_f12.port2, solutionOutlet.c) 
        annotation (points=[29.2,10; 36.8,10], style(pattern=0, thickness=2));
      connect(to_f22.port2, coolingInlet.c)     annotation (points=[19.2,-22;
            26.9,-22], style(pattern=0, thickness=2));
      connect(methanolSolution.p, to_f11.port1) 
        annotation (points=[-54,10; -39.2,10], style(pattern=0, thickness=2));
      connect(to_f21.port2, coolingOutlet.c) 
        annotation (points=[17.2,46; 30.8,46], style(pattern=0, thickness=2));
      connect(hx.f11, to_f11.port2) annotation (points=[-18.5,10; -24.8,10],
          style(pattern=0, thickness=2));
      connect(hx.f12, to_f12.port1) 
        annotation (points=[8.5,10; 14.8,10], style(pattern=0, thickness=2));
      connect(hx.f22, to_f22.port1) annotation (points=[-5,2; -5,-22; 4.8,-22],
          style(pattern=0, thickness=2));
      connect(to_f21.port1, hx.f21) annotation (points=[2.8,46; -5,46; -5,18],
          style(pattern=0, thickness=2));
      
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
      
    end TestCooler;
    
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
      
      Real h2o_liquid = pipe.segments[1].x[2]*pipe.segments[1].n_l;
      Real h2o_gas = pipe.segments[1].y[2]*pipe.segments[1].n_g;
    equation 
      connect(exitingConnector.port1, pipe.side2) annotation (points=[28.8,8;
            19.4,8; 19.4,8; 8,8],   style(pattern=0, thickness=2));
      connect(enteringConnector.port2, pipe.side1) annotation (points=[-38.8,8;
            -35.4,8; -35.4,8; -28,8],   style(pattern=0, thickness=2));
      connect(inlet.c, enteringConnector.port1) annotation (points=[-53,29;
            -64.5,29; -64.5,8; -53.2,8], style(pattern=0, thickness=2));
      connect(outlet.c, exitingConnector.port2) 
        annotation (points=[53.1,8; 43.2,8], style(pattern=0, thickness=2));
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
      
      import Modelica.SIunits.VolumeFlowRate;
      
      Mixer mixer annotation (extent=[-70,-86; -30,-62]);
      FuelTank fuelTank annotation (extent=[32,-88; 52,-68]);
      EnvironmentPort environmentPort annotation (extent=[70,-70; 86,-54]);
      FlowConnector flowConnector annotation (extent=[52,-70; 62,-62]);
      annotation (Diagram);
      FlowConnector flowConnector1 annotation (extent=[-24,-84; -12,-76]);
      Pump fuelPump 
                annotation (extent=[-10,-96; 2,-84]);
      FlowConnector flowConnector2 annotation (extent=[16,-94; 26,-86]);
      EnvironmentPort environmentPort1 annotation (extent=[-36,-52; -20,-36]);
      Pump solutionPump 
                 annotation (extent=[-80,-56; -68,-44]);
      FlowConnector flowConnector3 annotation (extent=[-70,-98; -60,-90]);
      FlowConnector flowConnector5 annotation (extent=[-14,-2; -6,6]);
      FlowConnector flowConnector6 annotation (extent=[-52,-52; -40,-44]);
      EnvironmentPort environmentPort2 annotation (extent=[20,-32; 36,-16]);
      GasFlowController gasFlowController annotation (extent=[-14,-34; -2,-22]);
      FlowConnector flowConnector7 annotation (extent=[-6,-14; 4,-6]);
      FlowConnector flowConnector8 annotation (extent=[4,-32; 14,-24]);
      FlowConnector flowConnector9 annotation (extent=[8,12; 18,20]);
      EnvironmentPort environmentPort3 annotation (extent=[24,12; 40,28]);
      Separator separator annotation (extent=[32,-12; 62,16]);
      FlowConnector flowConnector10 annotation (extent=[18,-2; 30,6]);
      FlowConnector flowConnector11 annotation (extent=[-24,-64; -12,-56]);
      FlowConnector flowConnector12 annotation (extent=[60,10; 70,18]);
      EnvironmentPort environmentPort4 annotation (extent=[72,10; 88,26]);
      Plug plug annotation (extent=[-24,-76; -12,-64]);
      Cooler cooler annotation (extent=[-4,-4; 16,8]);
      inner parameter Real RH_env = 60;
      inner parameter Modelica.SIunits.Pressure p_env = 101325;
      inner parameter Modelica.SIunits.Temperature T_env = 298.15;
      parameter VolumeFlowRate airFlow = 50e-3/60 
        "Cooler air flow, default 50 sL/min";
      parameter VolumeFlowRate solution = 30e-6/60 
        "Methanol solution flow, default 30 scc/min";
      parameter VolumeFlowRate methanol = 2e-6/60 
        "Pure methanol flow, default 2 scc/min";
    equation 
      connect(flowConnector.port2, environmentPort.c) annotation (points=[60.6,-66;
            70.8,-66],                  style(pattern=0, thickness=2));
      connect(flowConnector.port1, fuelTank.topFlow) annotation (points=[53.4,-66;
            42,-66; 42,-70],      style(pattern=0, thickness=2));
      connect(flowConnector1.port1, mixer.flow2) annotation (points=[-22.32,-80;
            -26,-80; -26,-96; -52,-96; -52,-83.6], style(pattern=0, thickness=2));
      connect(fuelPump.outlet, flowConnector1.port2) 
        annotation (points=[-4,-84; -4,-80; -13.68,-80],
                                                style(pattern=0, thickness=2));
      connect(fuelPump.inlet, flowConnector2.port1) 
                                                annotation (points=[-4.12,-90;
            17.4,-90], style(pattern=0, thickness=2));
      connect(flowConnector2.port2, fuelTank.bottomFlow) annotation (points=[24.6,-90;
            42,-90; 42,-86],           style(pattern=0, thickness=2));
      connect(solutionPump.inlet, flowConnector3.port1) 
                                                 annotation (points=[-74.12,-50;
            -74,-50; -74,-94; -68.6,-94], style(pattern=0, thickness=2));
      connect(flowConnector3.port2, mixer.flow1) annotation (points=[-61.4,-94;
            -56,-94; -56,-83.6], style(pattern=0, thickness=2));
      connect(flowConnector6.port1, mixer.topFlow) annotation (points=[-50.32,
            -48; -50,-48; -50,-64.4], style(pattern=0, thickness=2));
      connect(flowConnector7.port1, gasFlowController.outlet) annotation (
          points=[-4.6,-10; -8,-10; -8,-22], style(pattern=0, thickness=2));
      connect(flowConnector8.port1, gasFlowController.inlet) 
        annotation (points=[5.4,-28; -8.12,-28],style(pattern=0, thickness=2));
      connect(environmentPort3.c, flowConnector9.port2) 
        annotation (points=[24.8,16; 16.6,16], style(pattern=0, thickness=2));
      connect(flowConnector10.port2, separator.feed) annotation (points=[28.32,2;
            30.91,2; 30.91,2; 33.5,2],   style(pattern=0, thickness=2));
      connect(flowConnector11.port2, separator.liquidOutlet) annotation (points=[-13.68,
            -60; 0,-60; 0,-38; 57.5,-38; 57.5,-2.2],        style(pattern=0,
            thickness=2));
      connect(flowConnector6.port2, environmentPort1.c) annotation (points=[-41.68,
            -48; -35.2,-48],                          style(pattern=0,
            thickness=2));
      connect(environmentPort4.c, flowConnector12.port2) 
        annotation (points=[72.8,14; 68.6,14], style(pattern=0, thickness=2));
      connect(flowConnector12.port1, separator.gasOutlet) annotation (points=[61.4,14;
            58,14; 58,6.2; 57.5,6.2],          style(pattern=0, thickness=2));
      connect(mixer.flow4, flowConnector11.port1) annotation (points=[-44,-83.6;
            -44,-88; -36,-88; -36,-60; -22.32,-60],
                                               style(pattern=0, thickness=2));
      connect(plug.c, mixer.flow3) annotation (points=[-23.4,-70; -32,-70; -32,
            -92; -48,-92; -48,-83.6], style(pattern=0, thickness=2));
      connect(environmentPort2.c, flowConnector8.port2) annotation (points=[20.8,-28;
            12.6,-28],                                     style(pattern=0,
            thickness=2));
      connect(flowConnector5.port2, cooler.f11) annotation (points=[-7.12,2; -3,
            2],                      style(pattern=0, thickness=2));
      connect(cooler.f12, flowConnector10.port1) annotation (points=[15,2;
            19.68,2],                   style(pattern=0, thickness=2));
      connect(flowConnector7.port2, cooler.f22) annotation (points=[2.6,-10; 6,
            -10; 6,-1], style(pattern=0, thickness=2));
      connect(cooler.f21, flowConnector9.port1) 
        annotation (points=[6,5; 6,16; 9.4,16], style(pattern=0, thickness=2));
      gasFlowController.V = airFlow;
      solutionPump.V = solution;
      fuelPump.V = methanol;
      separator.level =0.5*separator.d;
      
      connect(solutionPump.outlet, flowConnector5.port1) annotation (points=[
            -74,-44; -74,2; -12.88,2], style(pattern=0, thickness=2));
    end SystemWithoutFC;
    
  end Tests;
  
end Tank;
