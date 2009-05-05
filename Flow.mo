import " Units.mo";


package Flow 
  
  connector TemperatureInput =  Modelica.Blocks.Interfaces.RealInput(redeclare 
        type SignalType = Modelica.SIunits.Temperature, start=298.15) annotation(defaultComponentName="T");
  connector VolumeFlowRateInput=Modelica.Blocks.Interfaces.RealInput(redeclare 
        type SignalType = Modelica.SIunits.VolumeFlowRate) 
                                                   annotation(defaultComponentName="V");
  connector ConcentrationOutput=Modelica.Blocks.Interfaces.RealOutput(redeclare 
        type SignalType = Modelica.SIunits.Concentration) annotation(defaultComponentName="c");
  connector TemperatureOutput = Modelica.Blocks.Interfaces.RealOutput(redeclare 
        type SignalType = Modelica.SIunits.Temperature, start=298.15) annotation(defaultComponentName="T");
  
  connector VolumeOutput =      Modelica.Blocks.Interfaces.RealOutput(redeclare 
        type SignalType = Modelica.SIunits.Volume) annotation(defaultComponentName="V");
  connector FlowPort "What passes through a control surface" 
    
    import Units.MolarFlow;
    import Modelica.SIunits.EnthalpyFlowRate;
    
    flow MolarFlow[size(Thermo.Molecules.All,1)] n;
    flow EnthalpyFlowRate H;
    
    annotation (Documentation(info="<html>
<p>This is a connector for the various tank units; it ensures continuity of enthalpy and
molar flows. It consists of two flow variables, the <em>enthalpy flow</em> and the array of 
<em>molar flows</em>.</p>
<p>The enthalpy flow is defined as the enthalpy necessary to bring the components in the 
molar-flow array from standard conditions, which is defined as starting from fundamental 
elements in their native state, to the actual ones. Using this definition, it is fairly 
easy to model chemical reactions.</p>
</html>"), Icon(Rectangle(extent=[-100,100; 100,-100],style(
            pattern=0,
            fillColor=62,
            rgbfillColor={0,127,127})), Text(
          extent=[-100,160; 100,100],
          string="%name",
          style(color=62, rgbcolor={0,127,127}))));
  end FlowPort;
  
  annotation (uses(Modelica(version="2.2.1")), Documentation(info="<html>
<p>This package contains various models related to fluid flow in stirred tanks.</p>
</html>"));
  
  model MethanolSolution "Source of methanol solution with given concentration" 
    import Thermo.mw;
    import Thermo.rho;
    import Thermo.h;
    import Thermo.Molecules.Incondensable;
    import Thermo.Molecules.Water;
    import Thermo.Molecules.Methanol;
    import Thermo.Phases.Gas;
    import Thermo.Phases.Liquid;
    import Modelica.SIunits.Temperature;
    import Modelica.SIunits.MoleFraction;
    import Modelica.SIunits.Concentration;
    
    FlowPort outlet "Methanol solution" 
      annotation (extent=[-20,-20; 20,20]);
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
</html>"),
      Diagram);
    
    outer parameter Temperature T_env = 298.15 "Environment temperature";
    
    parameter Concentration C = 1000 "Methanol concentration";
    parameter Temperature T = T_env "Temperature";
    
    MoleFraction x_ch3oh "Methanol molar fraction";
    MoleFraction x_h2o(start=1) "Water molar fraction";
    
  equation 
    assert(C >= 0, "==> Negative concentration given in MethanolSolution object.");
    assert(C <= rho(T,Methanol,Liquid)/mw(Methanol), "==> Methanol concentration over limit (" + String(mw(Methanol)/rho(T,Methanol,Liquid)) + " mol/m³).");
    
    C = x_ch3oh / (x_ch3oh*mw(Methanol)/rho(T,Methanol,Liquid) + x_h2o*mw(Water)/rho(T,Water,Liquid));
    x_ch3oh + x_h2o = 1.0;
    
    outlet.n[Incondensable] = zeros(size(Incondensable, 1));
    outlet.n[Methanol] * x_h2o = outlet.n[Water] * x_ch3oh;
    outlet.H = outlet.n[Methanol]*h(T,Methanol,Liquid) + outlet.n[Water]*h(T,Water,Liquid);
    
  end MethanolSolution;
  
  model PureMethanolSource "Source of pure methanol" 
    import Thermo.h;
    import Thermo.Molecules.Incondensable;
    import Thermo.Molecules.Water;
    import Thermo.Molecules.Methanol;
    import Thermo.Phases.Liquid;
    import Thermo.mw;
    import Thermo.rho;
    import Modelica.SIunits.Temperature;
    import Modelica.SIunits.AmountOfSubstance;
    import Modelica.SIunits.Mass;
    import Modelica.SIunits.Volume;
    
    FlowPort outlet "Pure methanol" 
      annotation (extent=[-20,-20; 20,20]);
    annotation (Icon(Ellipse(extent=[-100,100; 100,-100], style(
            pattern=0,
            thickness=4,
            fillColor=1,
            rgbfillColor={255,0,0}))),     Documentation(info="<html>
<p>This item is a source for a pure methanol stream.</p>
<p>This model also keeps track of how much methanol has been
released in terms of moles, mass and volume.</p>
</html>"),
      Diagram);
    
    outer parameter Temperature T_env = 298.15 "Enviroment temperature";
    
    AmountOfSubstance n(start=0,fixed=true) 
      "Released methanol moles from simulation start";
    Mass m = n*mw(Methanol) "Released methanol mass";
    Volume V = m / rho(T_env, Methanol, Liquid) "Released methanol volume";
    
  equation 
    der(n) = -outlet.n[Methanol];
    
    outlet.n[Incondensable] = zeros(size(Incondensable, 1));
    outlet.n[Water] = 0;
    outlet.H = outlet.n[Methanol]*h(T_env,Methanol,Liquid);
    
  end PureMethanolSource;
  
  model EnvironmentPort "A flow connection to environment conditions." 
    
    import Thermo.h;
    import Thermo.K;
    import Thermo.Molecules.Methanol;
    import Thermo.Molecules.Water;
    import Thermo.Molecules.Oxygen;
    import Thermo.Molecules.CarbonDioxide;
    import Thermo.Molecules.Nitrogen;
    import Thermo.Phases.Gas;
    import Thermo.Phases.Liquid;
    import Units.MolarEnthalpy;
    import Modelica.SIunits.Temperature;
    import Modelica.SIunits.Pressure;
    import Modelica.SIunits.MoleFraction;
    
    FlowPort outlet "Environment-air port" 
                 annotation (extent=[80,-60; 100,-40]);
    annotation (Documentation(info="<html>
<p>This object generates a gas flow corresponding to ambient air, accounting
also for humidity.</p>
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
            rgbfillColor={255,255,85})),
        Text(
          extent=[-100,160; 100,100],
          style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=43,
            rgbfillColor={255,85,85},
            fillPattern=1),
          string="%name")),
      Diagram);
    
    outer Real RH_env "Environment relative humidity";
    outer Temperature T_env "Environment temperature";
    outer Pressure p_env "Environment pressure";
    
    MoleFraction y_h2o "Water molar fraction";
    MoleFraction y_o2 "Oxygen molar fraction";
    MoleFraction y_n2 "Nitrogen molar fraction";
    
    MolarEnthalpy h_air "Air molar enthalpy";
    
  protected 
    MolarEnthalpy h_n2 "Nitrogen molar enthalpy";
    MolarEnthalpy h_o2 "Oxygen molar enthalpy";
    MolarEnthalpy h_h2o "Water-vapour molar enthalpy";
    
  equation 
    y_o2 / 0.21 = y_n2 / 0.79; // The O2/N2 ratio.
    y_h2o + y_o2 + y_n2 = 1.0; // Fractions sum to 1.
    y_h2o = RH_env/100 * K(T_env, Water); // Humidity of air
    
    h_h2o = h(T_env, Water, Gas);
    h_o2  = h(T_env, Oxygen, Gas);
    h_n2  = h(T_env, Nitrogen, Gas);
    h_air = h_h2o*y_h2o + h_o2*y_o2 + h_n2*y_n2;
    
    // Write with * instead of / to avoid division by zero.
    outlet.n[Water] * y_o2 = outlet.n[Oxygen] * y_h2o;
    outlet.n[Nitrogen] / y_n2 = outlet.n[Oxygen] / y_o2;
    outlet.n[Methanol] = 0;
    outlet.n[CarbonDioxide] = 0;
    outlet.H = h_air*sum(outlet.n);
    
  end EnvironmentPort;
  
  model SinkPort 
    
    annotation (Diagram, Icon(
        Rectangle(extent=[-100,100; 100,-100], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4,
            fillColor=7,
            rgbfillColor={255,255,255},
            fillPattern=1)),
        Line(points=[-40,60; 60,60; 60,-60; -40,-60], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4,
            fillColor=7,
            rgbfillColor={255,255,255},
            fillPattern=1)),
        Rectangle(extent=[-60,0; 60,0], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4,
            fillColor=7,
            rgbfillColor={255,255,255},
            fillPattern=1))),
      Documentation(info="<html>
<p>This very simple model is a terminal for flows leaving the system and
about which we do not care much.</p>
</html>"));
    FlowPort inlet "inlet for flow to discard" 
                      annotation (extent=[-120,-30; -60,30]);
  end SinkPort;
  
  model FlowController "A unit modelling a pump or a MFC" 
    
    import Thermo.mw;
    import Thermo.Molecules.All;
    import Modelica.SIunits.MassFlowRate;
    import Units.MolarFlow;
    
    annotation (Icon(
        Ellipse(extent=[-100,100; 100,-100],
                                         style(
            color=0,
            rgbcolor={0,0,0},
            thickness=2,
            fillColor=7,
            rgbfillColor={255,255,255})),Text(
          extent=[-100,160; 100,100],
          style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=43,
            rgbfillColor={255,85,85},
            fillPattern=1),
          string="%name")),   Documentation(info="<html>
<p>This class allows to set or read a certain overall molar or mass flow. It is not
immediately possible to set a <em>volume</em> flow, because this would entail
calculating the phase equilibrium, which we are not doing here since it is
computationally onerous; see the child classes <tt>Pump</tt> and 
<tt>GasFlowController</tt> if you need volume.</p>
</html>"),
      Diagram);
    FlowPort inlet "Unit inlet"   annotation (extent=[-10,-10; 10,10]);
    FlowPort outlet "Unit outlet"   annotation (extent=[-10,90; 10,110]);
    
    MolarFlow F "Molar flow rate";
    MassFlowRate m "Mass flow rate";
    
  equation 
    m = sum({inlet.n[i] * mw(i) for i in All});
    F = sum(inlet.n);
    
    connect(inlet, outlet) annotation (points=[5.55112e-16,5.55112e-16; 
          5.55112e-16,25; 5.55112e-16,25; 5.55112e-16,50; 5.55112e-16,100; 
          5.55112e-16,100], style(color=62, rgbcolor={0,127,127}));
  end FlowController;
  
  model GasFlowController "A flow controller with only gas phase" 
    extends FlowController;
    
    import Thermo.rho;
    import Thermo.mw;
    import Thermo.Molecules.All;
    import Thermo.Phases.Gas;
    import Modelica.SIunits.VolumeFlowRate;
    import Modelica.SIunits.Temperature;
    
    annotation (Icon,     Documentation(info="<html>
<p>This class implements a mass flow controller with volumetric units (\"field units\").
Since there are at least <em>two</em> different standards (the one named \"standard\" at 
0 Celsius and the one named \"norm\" at 70 Fahrenheit or 21.111... Celsius), it is 
necessary to provide the reference temperature; the default assumes zero Celsius 
(the so-called \"standard\" value).</p>

<p>The flow assumes that all components are in gas phase and takes their density from
the Thermo library, where the ideal gas law is (usually) assumed.</p>
</html>"), 
      Diagram);
    
    parameter Temperature T_ref = 273.15 "Reference temperature";
    
    VolumeFlowRateInput V "Volumetric flow rate"
      annotation (extent=[110,-10; 90,10]);
    
  equation 
    V = sum({inlet.n[i] * mw(i) / rho(T_ref, i, Gas) for i in All});
    
  end GasFlowController;
  
  model Pump "A pump, with only liquid phase." 
    extends FlowController;
    
    import Thermo.rho;
    import Thermo.mw;
    import Thermo.h;
    import Thermo.Phases.Liquid;
    import Thermo.Molecules.Condensable;
    import Modelica.SIunits.VolumeFlowRate;
    import Modelica.SIunits.Temperature;
    
    annotation (Icon,        Documentation(info="<html>
<p>This class implements a liquid pump with field volumetric units.</p>
<p>The flow assumes that only water and methanol are present and are completely
in liquid phase; it takes their density from the Thermo library.</p>
</html>"), 
      Diagram);
    
    Temperature T(start=298.15) "Temperature of the flow";
    VolumeFlowRateInput V "Volumetric flow rate" 
      annotation (extent=[110,-10; 90,10]);
    
  equation 
    V = sum({inlet.n[i] * mw(i) / rho(T, i, Liquid) for i in Condensable});
    // This is to find the temperature.
    inlet.H = (sum(inlet.n[i]*h(T, i, Liquid) for i in Condensable));
    
  end Pump;
  
    model FlowTemperature "Calculates a flow's temperature" 
      extends Modelica.Icons.RotationalSensor;
    
      import Modelica.SIunits.MoleFraction;
      import Modelica.SIunits.Temperature;
      import Thermo.Molecules.All;
      import Thermo.Molecules.Incondensable;
      import Thermo.Molecules.Condensable;
      import Thermo.Phases.Gas;
      import Thermo.Phases.Liquid;
      import Thermo.Molecules.Water;
      import Thermo.Molecules.Methanol;
      import Thermo.h;
      import Thermo.rr;
      import Thermo.K;
      import Units.MolarFlow;
    
      annotation (Diagram, Icon(
          Text(extent=[-100,100; 100,140],string="%name",
            style(color=0, rgbcolor={0,0,0}))),
      Documentation(info="<html>
<p>This basic unit takes a flow and returns it unchanged, while actually
performing an equilibrium calculation and figuring out the temperature of
the flow given its associated enthalpic flow.</p>
<p>This unit can be used in different situations, both to set or to measure
the temperature of a flow. The unit can also work with reverse flows, but
it might give divide-by-zero's during an inversion.</p>
<p><strong>Note</strong>: remember to use this unit as sparingly as possible,
since every measurement entails a full calculation of the multicomponent
equilibrium.</p>
</html>"));
      FlowPort inlet "Entering flow" 
                     annotation (extent=[-90,-10; -70,10]);
      FlowPort outlet "Exiting flow" 
                     annotation (extent=[70,-10; 90,10]);
    
      MoleFraction beta "Vapour fraction";
    
      MolarFlow[size(All,1)] vapour "Gas-phase flows";
      MolarFlow[size(All,1)] liquid "Liquid-phase flows";
    
  protected 
      MoleFraction z_m(start=0) "Overall methanol molar fraction";
      MoleFraction z_w(start=0.1) "Overall water molar fraction";
    
  public 
      TemperatureOutput T "The flow temperature"
      annotation (extent=[46,-74; 30,-54]);
    equation 
      // Liquid and vapour sum to overall flow
      liquid + vapour       = inlet.n;
      // Condensable species are present in both phases
      vapour[Condensable]   = {inlet.n[i]*beta*K(T,i) / (1 + beta*(K(T,i) -1)) for i in Condensable};
      // Incondensable species are only in gas phase
      liquid[Incondensable] = zeros(size(Incondensable,1));
    
      /* Enforce definition of molar fractions. Avoid using divisions
   * in order to avoid division-by-zero errors. */
      inlet.n[Methanol] = sum(inlet.n) * z_m;
      inlet.n[Water]    = sum(inlet.n) * z_w;
    
      // Combine enthalpic flow with temperature, thermodynamic data and material flows.
      inlet.H = sum( vapour[i] * h(T, i, Gas)    for i in All)
              + sum( liquid[i] * h(T, i, Liquid) for i in Condensable);
    
      /* Use the beta returned by Thermo.rr only if:
   * 1) there can be a phase equilibrium at all at this temperature;
   * 2) the returned beta is less than 1. */
      if K(T,Water) >= 1 or rr(z_m, z_w, T) >= 1 then
        beta = 1;
      else
        beta = rr(z_m, z_w, T);
      end if;
    
      connect(inlet, outlet) 
                          annotation (points=[-80,5.55112e-16; 0,-4.87687e-22; 
          0,5.55112e-16; 80,5.55112e-16],  style(color=62, rgbcolor={0,127,127}));
    end FlowTemperature;
  
    model FlowConcentration "Calculates a liquid flow's methanol concentration" 
      extends FlowTemperature;
    
      import Modelica.SIunits.Concentration;
      import Thermo.Molecules.Condensable;
      import Thermo.Phases.Liquid;
      import Thermo.Molecules.Methanol;
      import Thermo.rho;
      import Thermo.mw;
    
      annotation (Diagram, Icon,
      Documentation(info="<html>
<p>Adds the ability to read the methanol concentration <em>in the liquid 
phase</em> to the temperature measurement of <tt>FlowTemperature</tt>.</p>
<p>If there is no liquid flow, then the reported value is zero.</p>
</html>"));
    
    ConcentrationOutput c annotation (extent=[-46,-74; -30,-54]);
    equation 
      // Methanol flow is concentration times volumetric flow
      liquid[Methanol] = c * sum(liquid[i]*mw(i)/rho(T,i,Liquid) for i in Condensable);
    
    end FlowConcentration;
  
  model Mixer "A unit mixing four molar flows." 
    
    import Modelica.SIunits.Temperature;
    import Modelica.SIunits.AmountOfSubstance;
    import Modelica.SIunits.Concentration;
    import Modelica.SIunits.Volume;
    import Thermo.Molecules.Methanol;
    import Thermo.Molecules.Incondensable;
    import Thermo.Molecules.Condensable;
    import Thermo.Phases.Liquid;
    import Thermo.h;
    import Thermo.mw;
    import Thermo.rho;
    
    FlowPort outlet "The mixer's outlet" 
                          annotation (extent=[-90,-10; -70,10]);
    FlowPort fuelInlet "The methanol-feed inlet" 
                           annotation (extent=[-10,-90; 10,-70]);
    FlowPort loopInlet "The anode loop's inlet" 
                           annotation (extent=[-10,70; 10,90]);
    FlowPort waterInlet "The water-recovery inlet" 
                           annotation (extent=[70,-10; 90,10]);
    annotation (Diagram, Icon(
        Ellipse(extent=[-80,80; 80,-80], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4,
            fillColor=7,
            rgbfillColor={255,255,255})),
        Rectangle(extent=[-80,0; 80,80], style(
            pattern=0,
            thickness=4,
            fillColor=7,
            rgbfillColor={255,255,255},
            fillPattern=1)),
        Line(points=[-80,0; -80,80; 80,80; 80,0], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4,
            fillColor=7,
            rgbfillColor={255,255,255},
            fillPattern=1)),
        Line(points=[0,6; 0,-54], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4,
            fillColor=7,
            rgbfillColor={255,255,255},
            fillPattern=1)),
        Line(points=[0,6; -52,36], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4,
            fillColor=7,
            rgbfillColor={255,255,255},
            fillPattern=1)),
        Line(points=[0,6; 52,36], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4,
            fillColor=7,
            rgbfillColor={255,255,255},
            fillPattern=1)),
        Text(
          extent=[-100,160; 100,100],
          style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=43,
            rgbfillColor={255,85,85},
            fillPattern=1),
          string="%name")),
      Documentation(info="<html>
<p>The mixer features three input flows, one for the anodic loop, one for recovered
water and another for methanol inlet, and an output for the anodic loop. Its main
role is thought to be the loop's main mass and energy balance.</p>
<p>The outlet compositions are set to be the same as the mass balance's molar fractions, 
implying a perfect mixing; the enthalpy flow is also proportional to the internal-energy
holdup.</p>
<p>It is possible to set the initial methanol concentration to some specific value, 
by default it is 1 M.</p>
</html>"));
    
    AmountOfSubstance n[size(Thermo.Molecules.All,1)] "Molar holdup";
    Modelica.SIunits.InternalEnergy U "Internal energy";
    Concentration c(start=1000) "Methanol concentration";
    Temperature T(start=298.15) "Mixer temperature";
    VolumeOutput V(start=5E-6) "Solution volume"
      annotation (extent=[-52,-82; -36,-62]);
    
  equation 
    der(U) = fuelInlet.H + loopInlet.H + waterInlet.H + outlet.H;
    der(n) = fuelInlet.n + loopInlet.n + waterInlet.n + outlet.n;
    
    // Bind outlet's n to composition in holdup
    outlet.n[1:end-1] / sum(outlet.n) = n[1:end-1] / sum(n);
    // Bind outlet's H to specific internal energy and outlet flow
    outlet.H / sum(outlet.n) = U / sum(n);
    
    U = sum(n[i]*h(T,i,Liquid) for i in Condensable);
    V = sum( n[i]*mw(i)/rho(T, i, Liquid) for i in Condensable);
    c = n[Methanol] / V;
    
  initial equation 
    n[Incondensable] = zeros(size(Incondensable,1));
    
  end Mixer;
  
  model Separator 
    import Modelica.SIunits.Temperature;
    import Thermo.Molecules.Condensable;
    import Thermo.Phases.Liquid;
    import Thermo.h;
    
    FlowPort inlet "Two-phase inlet" 
                   annotation (extent=[-110,-10; -90,10]);
    FlowPort gasOutlet "Single-phase gas outlet" 
                       annotation (extent=[60,30; 80,50]);
    FlowPort liquidOutlet "Single-phase liquid outlet" 
                          annotation (extent=[60,-50; 80,-30]);
    annotation (Icon(
        Ellipse(extent=[-100,40; -60,-40], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4,
            fillColor=7,
            rgbfillColor={255,255,255},
            fillPattern=1)),
        Ellipse(extent=[60,40; 100,-40], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4,
            fillColor=7,
            rgbfillColor={255,255,255},
            fillPattern=1)),
        Rectangle(extent=[-80,40; 80,-40], style(
            color=0,
            rgbcolor={255,255,255},
            thickness=0,
            fillColor=7,
            rgbfillColor={255,255,255},
            fillPattern=1)),
        Line(points=[-80,40; 80,40], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4)),
        Line(points=[-80,-40; 80,-40], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4)),
        Text(
          extent=[-100,100; 100,40],
          style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=43,
            rgbfillColor={255,85,85},
            fillPattern=1),
          string="%name")), Diagram,
      Documentation(info="<html>
<p>The separator unit simply splits a flow in its gaseous and liquid components. 
The separation criterion is straightforwardly the liquid-vapor equilibrium.</p>
</html>"));
    TemperatureOutput T "Separator temperature" annotation (extent=[100,-10; 120,10]);
  protected 
    FlowTemperature ft annotation (extent=[-10,-10; 10,10]);
    
  equation 
    liquidOutlet.n = -ft.liquid;
    liquidOutlet.H = sum(h(T, i, Liquid)* liquidOutlet.n[i] for i in Condensable);
    
    connect(ft.inlet, inlet) annotation (points=[-8,6.10623e-16; -54,
          6.10623e-16; -54,5.55112e-16; -100,5.55112e-16], style(color=62,
          rgbcolor={0,127,127}));
    connect(ft.outlet, gasOutlet) annotation (points=[8,6.10623e-16; 40,
          6.10623e-16; 40,40; 70,40], style(color=62, rgbcolor={0,127,127}));
    connect(ft.outlet, liquidOutlet) annotation (points=[8,6.10623e-16; 40,
          6.10623e-16; 40,-40; 70,-40], style(color=62, rgbcolor={0,127,127}));
    connect(T, ft.T) annotation (points=[110,5.55112e-16; 78,0; 46,0; 46,-6.4; 3.8,
          -6.4], style(color=3, rgbcolor={0,0,255}));
  end Separator;
  
  partial model AbstractHeatExchanger "An abstract heat exchanger" 
    import Modelica.SIunits.HeatFlowRate;
    import Modelica.SIunits.Area;
    import Modelica.SIunits.Temperature;
    import Units.HeatTransferCoefficient;
    
    FlowPort hot_1 "Port for hot flow on side 1" 
                   annotation (extent=[-100,20; -80,40]);
    FlowPort hot_2 "Port for hot flow on side 2" 
                    annotation (extent=[-100,-40; -80,-20]);
    annotation (Diagram, Icon(
        Rectangle(extent=[-100,40; 100,-40],
                                           style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4,
            fillColor=7,
            rgbfillColor={255,255,255})),
        Text(
          extent=[-100,100; 100,40],
          style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=43,
            rgbfillColor={255,85,85},
            fillPattern=1),
          string="%name"),
        Line(points=[-60,40; -60,-40], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4)),
        Line(points=[20,40; 20,-40], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4)),
        Line(points=[-20,40; -20,-40], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4)),
        Line(points=[60,40; 60,-40], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4)),
        Line(points=[-40,40; -40,-40], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4)),
        Line(points=[40,40; 40,-40], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4)),
        Line(points=[0,40; 0,-40], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4)),
        Line(points=[80,40; 80,-40], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4)),
        Line(points=[-80,40; -80,-40], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4))),
      Documentation(info="<html>
<p>This is the interface for a generic heat exchanger. It constraints the two sides to
be able to exchange a heat duty <tt>Q</tt>, but no mass.</p>
<p>All subclasses may be counter-current or co-current exchangers depending on the
connections of the flows to the two sides.</p>
</html>"));
    FlowPort cold_1 "Port for cold flow on side 2" 
                   annotation (extent=[80,20; 100,40]);
    FlowPort cold_2 "Port for cold flow on side 2" 
                   annotation (extent=[80,-40; 100,-20]);
    TemperatureOutput T_hot_1 "Temperature of hot stream on side 1"
      annotation (extent=[-80,40; -100,60]);
    TemperatureOutput T_hot_2 "Temperature of hot stream on side 2"
      annotation (extent=[-80,-60; -100,-40]);
    TemperatureOutput T_cold_1 "Temperature of cold stream on side 1"
      annotation (extent=[80,40; 100,60]);
    TemperatureOutput T_cold_2 "Temperature of cold stream on side 2"
      annotation (extent=[80,-60; 100,-40]);
    
    parameter Area A = 2.3E-2 "Effective heat-transfer area";
    parameter HeatTransferCoefficient U = 188 "Heat-transfer coefficient";
    
    HeatFlowRate Q "Heating duty from hot side to cold side";
    
  equation 
    hot_1.H  + hot_2.H  - Q = 0;
    cold_1.H + cold_2.H + Q = 0;
    hot_1.n  + hot_2.n      = 0*hot_1.n;
    cold_1.n + cold_2.n     = 0*hot_1.n;
    
  end AbstractHeatExchanger;
  
  model LMTDHeatExchanger "A heat exchanger based on the LMTD" 
    extends AbstractHeatExchanger;
    
    import Modelica.SIunits.Temperature;
    import Modelica.Math.log;
    
    Temperature LMTD = (DT1-DT2)/log(DT1/DT2) "Log mean temperature difference";
    
  protected 
    Temperature DT1 = T_hot_1 - T_cold_1 "Temperature difference on side 1";
    Temperature DT2 = T_hot_2 - T_cold_2 "Temperature difference on side 2";
    
  protected 
    FlowTemperature hot_1_T "Temperature for hot flow on side 1" 
      annotation (extent=[-60,20; -40,40]);
    FlowTemperature cold_1_T "Temperature for cold flow on side 1" 
      annotation (extent=[40,20; 60,40]);
    FlowTemperature hot_2_T "Temperature for hot flow on side 2" 
      annotation (extent=[-60,-40; -40,-20]);
    FlowTemperature cold_2_T "Temperature for cold flow on side 2" 
      annotation (extent=[40,-40; 60,-20]);
    annotation (Diagram, DymolaStoredErrors,
      Documentation(info="<html>
<p>This heat exchanger assumes that LMTD theory is valid and uses it to
calculate the heating or cooling duty.</p>
<p>The main limitation is that the specific heat capacity of the two fluids
must be constant (not necessarily the same): this exchanger is ok for air,
methanol solution (there is usually too little evaporation to substantially
change the c<sub>p</sub>), but <em>not</em> for cathode outlet, since water
condensation is very important there.</p>
</html>"));
    SinkPort sinkPort annotation (extent=[20,-80; 40,-60]);
    
  equation 
    Q = U * A * LMTD;
    
    connect(hot_2_T.outlet, sinkPort.inlet)    annotation (points=[-42,-30; 0,
          -30; 0,-70; 21,-70], style(color=62, rgbcolor={0,127,127}));
    connect(cold_2_T.inlet, sinkPort.inlet)   annotation (points=[42,-30; 0,-30; 
          0,-70; 21,-70], style(color=62, rgbcolor={0,127,127}));
    connect(cold_1_T.inlet, sinkPort.inlet)  annotation (points=[42,30; 0,30; 0,
          -70; 21,-70], style(color=62, rgbcolor={0,127,127}));
    connect(hot_1_T.outlet, sinkPort.inlet)   annotation (points=[-42,30; 0,30; 
          0,-70; 21,-70], style(color=62, rgbcolor={0,127,127}));
    connect(cold_1, cold_1_T.outlet) 
      annotation (points=[90,30; 58,30], style(color=62, rgbcolor={0,127,127}));
    connect(cold_2_T.outlet, cold_2) annotation (points=[58,-30; 90,-30], style(
          color=62, rgbcolor={0,127,127}));
    connect(hot_2_T.inlet, hot_2) annotation (points=[-58,-30; -90,-30], style(
          color=62, rgbcolor={0,127,127}));
    connect(hot_1_T.inlet, hot_1) annotation (points=[-58,30; -90,30], style(
          color=62, rgbcolor={0,127,127}));
    connect(hot_1_T.T, T_hot_1) annotation (points=[-46.2,23.6; -46,23.6; -46,
          20; -40,20; -40,50; -90,50], style(color=3, rgbcolor={0,0,255}));
    connect(cold_1_T.T, T_cold_1) annotation (points=[53.8,23.6; 54,23.6; 54,20; 
          40,20; 40,50; 90,50], style(color=3, rgbcolor={0,0,255}));
    connect(cold_2_T.T, T_cold_2) annotation (points=[53.8,-36.4; 54,-36.4; 54,
          -50; 90,-50], style(color=3, rgbcolor={0,0,255}));
    connect(hot_2_T.T, T_hot_2) annotation (points=[-46.2,-36.4; -46,-36.4; -46,
          -50; -90,-50], style(color=3, rgbcolor={0,0,255}));
  end LMTDHeatExchanger;
  
  model DiscretisedHeatExchangerStep 
    "A step in a larger, discretised heat exchanger" 
    extends AbstractHeatExchanger;
    
    import Modelica.SIunits.Temperature;
    
    Temperature T_hot =  (hot_1_T.T  + hot_2_T.T)/2 
      "Average hot-side temperature";
    Temperature T_cold = (cold_1_T.T + cold_2_T.T)/2 
      "Average cold-side temperature";
    
  protected 
    FlowTemperature hot_1_T "Temperature of the hot flow on side 1" 
      annotation (extent=[-60,20; -40,40]);
    FlowTemperature cold_1_T "Temperature of the cold flow on side 1" 
      annotation (extent=[40,20; 60,40]);
    FlowTemperature hot_2_T "Temperature of the hot flow on side 2" 
      annotation (extent=[-60,-40; -40,-20]);
    FlowTemperature cold_2_T "Temperature of the cold flow on side 2" 
      annotation (extent=[40,-40; 60,-20]);
    annotation (Diagram, Documentation(info="<html>
<p>To implement a discretised heat exchanger, a single step is implemented 
here as a simplistic heat exchanger. It could be based on LMTD, but using
the average temperature is numerically more robust and allows Dymola to
perform algebraic manipulation.</p>
</html>"));
    SinkPort sink "Sink element" 
                      annotation (extent=[20,-80; 40,-60]);
    
  equation 
    Q = U*A*(T_hot-T_cold);
    
    connect(hot_2_T.inlet, hot_2) annotation (points=[-58,-30; -90,-30], style(
          color=62, rgbcolor={0,127,127}));
    connect(hot_1_T.inlet, hot_1) annotation (points=[-58,30; -90,30], style(
          color=62, rgbcolor={0,127,127}));
    connect(cold_1_T.outlet, cold_1) annotation (points=[58,30; 90,30], style(
          color=62, rgbcolor={0,127,127}));
    connect(cold_2_T.outlet, cold_2) annotation (points=[58,-30; 90,-30], style(
          color=62, rgbcolor={0,127,127}));
    connect(sink.inlet, hot_2_T.outlet) annotation (points=[21,-70; 0,-70; 0,-30;
          -42,-30],      style(color=62, rgbcolor={0,127,127}));
    connect(sink.inlet, cold_2_T.inlet) annotation (points=[21,-70; 0,-70; 0,-30;
          42,-30],      style(color=62, rgbcolor={0,127,127}));
    connect(sink.inlet, hot_1_T.outlet) annotation (points=[21,-70; 0,-70; 0,30;
          -42,30], style(color=62, rgbcolor={0,127,127}));
    connect(sink.inlet, cold_1_T.inlet) annotation (points=[21,-70; 0,-70; 0,30;
          42,30], style(color=62, rgbcolor={0,127,127}));
    connect(cold_1_T.T, T_cold_1) annotation (points=[53.8,23.6; 54,23.6; 54,20; 
          40,20; 40,50; 90,50],   style(color=3, rgbcolor={0,0,255}));
    connect(cold_2_T.T, T_cold_2) annotation (points=[53.8,-36.4; 54,-36.4; 54,
          -50; 90,-50],                    style(color=3, rgbcolor={0,0,255}));
    connect(hot_2_T.T, T_hot_2) annotation (points=[-46.2,-36.4; -46,-36.4; -46,
          -50; -90,-50],                      style(color=3, rgbcolor={0,0,255}));
    connect(hot_1_T.T, T_hot_1) annotation (points=[-46.2,23.6; -46,23.6; -46,20; -40,20; 
          -40,50; -90,50], style(color=3, rgbcolor={0,0,255}));
  end DiscretisedHeatExchangerStep;
  
  model DiscretisedHeatExchanger 
    "A heat exchanger made up of many discrete sections" 
    extends AbstractHeatExchanger;
    
    parameter Integer n(min=3) = 5 "Number of discretisation units";
    
  protected 
    DiscretisedHeatExchangerStep[n] steps(each A=A/n, each U=U) 
      "The steps the exchanger is divided in";
    
  protected 
    SinkPort sink annotation (extent=[20,-80; 40,-60]);
    annotation (Diagram, Documentation(info="<html>
<p>This heat exchanger is much more complex than a LMTD exchanger, and should be
used only for condensing flows where LMTD theory is not valid.</p>
</html>"), 
      Icon);
  equation 
    Q = sum(steps[i].Q for i in 1:n);
    
  // Connecting => temperatures <=
    connect(T_hot_1, steps[1].T_hot_1);
    connect(T_hot_2, steps[end].T_hot_2);
    connect(T_cold_1, steps[1].T_cold_1);
    connect(T_cold_2, steps[end].T_cold_2);
    
  // Connecting => flows <=
    // Connect first step to exchanger's side 1
    connect(steps[1].hot_1,  hot_1);
    connect(steps[1].cold_1, cold_1);
    
    // Connect the elements among each other
    for i in 1:(n-1) loop
      connect(steps[i].hot_2,  steps[i+1].hot_1);
      connect(steps[i].cold_2, steps[i+1].cold_1);
    end for;
    
    // Connect last step to exchanger's side 2
    connect(steps[n].hot_2,  hot_2);
    connect(steps[n].cold_2, cold_2);
    
  // Connecting => sink <=
    // Add some degrees of freedom; it's a "dirty hack" but it works fine.
    connect(sink.inlet, cold_2) annotation (points=[21,-70; 0,-70; 0,-30; 90,
          -30], style(color=62, rgbcolor={0,127,127}));
    connect(sink.inlet, hot_2) annotation (points=[21,-70; 0,-70; 0,-30; -90,
          -30], style(color=62, rgbcolor={0,127,127}));
    
  end DiscretisedHeatExchanger;
  
  partial model AbstractCooler "An abstract heat exchanger" 
    
    import Modelica.SIunits.Temperature;
    import Modelica.SIunits.VolumeFlowRate;
    
    VolumeFlowRate V_air = mfc.V "Coolant flow rate";
    Temperature T_coolant_in =  exchanger.T_cold_2 "Coolant inlet temperature";
    Temperature T_coolant_out = exchanger.T_cold_1 "Coolant outlet temperature";
    
    FlowPort inlet "Inlet to the cooler" 
                   annotation (extent=[-100,-6; -88,6]);
    FlowPort outlet "Outlet from the cooler" 
                    annotation (extent=[88,-6; 100,6]);
    annotation (Diagram, Icon(
        Rectangle(extent=[-90,20; 90,-20], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4,
            fillColor=7,
            rgbfillColor={255,255,255})),
        Ellipse(extent=[-60,16; -28,-16], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4,
            fillColor=7,
            rgbfillColor={255,255,255})),
        Ellipse(extent=[28,16; 60,-16], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4,
            fillColor=7,
            rgbfillColor={255,255,255})),
        Rectangle(extent=[-44,18; -26,-18], style(
            pattern=0,
            thickness=4,
            fillColor=7,
            rgbfillColor={255,255,255},
            fillPattern=1)),
        Rectangle(extent=[26,18; 44,-18], style(
            pattern=0,
            thickness=4,
            fillColor=7,
            rgbfillColor={255,255,255},
            fillPattern=1)),
        Line(points=[-44,16; 44,-16], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4,
            fillColor=0,
            rgbfillColor={0,0,0},
            fillPattern=1)),
        Line(points=[-44,-16; 44,16], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4,
            fillColor=0,
            rgbfillColor={0,0,0},
            fillPattern=1)),
        Text(
          extent=[-102,82; 100,20],
          style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=43,
            rgbfillColor={255,85,85},
            fillPattern=1),
          string="%name")),
      Documentation(info="<html>
<p>This is the interface for a generic air cooler. It only specifies that the outlet
species is equal to the inlet ones, and that the difference in enthalpic flow is the
cooling or heating duty.</p>
</html>"));
    TemperatureOutput T_process_in "Process inlet temperature" annotation (extent=[-88,8; 
          -100,20]);
    TemperatureOutput T_process_oun "Process outlet temperature" annotation (extent=[88,8; 
          100,20]);
    replaceable AbstractHeatExchanger exchanger 
      "The heat exchanger implementing the cooler" 
      annotation (extent=[-64,28; 16,108]);
  protected 
    EnvironmentPort env "Environmental air source" 
      annotation (extent=[0,0; 20,20]);
    SinkPort airSink "Air outlet sink" annotation (extent=[50,70; 70,90]);
    
  public 
    GasFlowController mfc "Mass flow controller for cooling air" 
      annotation (extent=[20,20; 40,40]);
    TemperatureInput T_ref "Reference temperature for the process outlet"
      annotation (extent=[100,-20; 88,-8]);
  equation 
    connect(exchanger.hot_2, outlet) annotation (points=[-60,56; -60,
          -2.22045e-16; 94,-2.22045e-16], style(color=62, rgbcolor={0,127,127}));
    connect(inlet, exchanger.hot_1) annotation (points=[-94,-2.22045e-16; -80,
          -2.22045e-16; -80,80; -60,80], style(color=62, rgbcolor={0,127,127}));
    connect(airSink.inlet, exchanger.cold_1) annotation (points=[51,80; 12,80],
        style(color=62, rgbcolor={0,127,127}));
    connect(env.outlet, mfc.inlet) annotation (points=[19,5; 30,5; 30,30],
        style(color=62, rgbcolor={0,127,127}));
    connect(mfc.outlet, exchanger.cold_2) annotation (points=[30,40; 30,56; 12,
          56],
        style(color=62, rgbcolor={0,127,127}));
    connect(exchanger.T_hot_1, T_process_in) annotation (points=[-60,88; -72,88; 
          -72,14; -94,14], style(color=3, rgbcolor={0,0,255}));
    connect(exchanger.T_cold_1, T_process_oun) annotation (points=[12,88; 40,88; 
          40,94; 80,94; 80,14; 94,14], style(color=3, rgbcolor={0,0,255}));
  end AbstractCooler;
  
  model LMTDCooler "A heat exchanger based on the LMTD" 
    extends AbstractCooler(redeclare LMTDHeatExchanger exchanger);
    annotation (Diagram, Documentation(info="<html>
<p>A cooler that uses a LMTD implementation. It can for instance be applied as the
anode-loop cooler or as a recuperating heat exchanger on the anode side.</p>
</html>"),
      Icon);
  end LMTDCooler;
  
  model DiscretisedCooler 
    extends AbstractCooler(redeclare DiscretisedHeatExchanger exchanger);
    annotation (Documentation(info="<html>
<p>A cooler using a discretised implementation. It can for instance be applied as the
cathode-loop cooler (condenser).</p>
</html>"));
  end DiscretisedCooler;
  
  partial model FuelCell "A generic DMFC" 
    
    import Modelica.SIunits.Area;
    import Modelica.SIunits.Current;
    import Modelica.SIunits.CurrentDensity;
    import Modelica.SIunits.Concentration;
    import Modelica.SIunits.DiffusionCoefficient;
    import Modelica.SIunits.FaradayConstant;
    import Modelica.SIunits.HeatCapacity;
    import Modelica.SIunits.Length;
    import Modelica.SIunits.MoleFraction;
    import Modelica.SIunits.PartialPressure;
    import Modelica.SIunits.Pressure;
    import Modelica.SIunits.Temperature;
    import Modelica.SIunits.Voltage;
    import Modelica.Constants.eps;
    import Modelica.Electrical.Analog.Interfaces.PositivePin;
    import Modelica.Electrical.Analog.Interfaces.NegativePin;
    
    import Thermo.mw;
    import Thermo.rho;
    import Thermo.Molecules.All;
    import Thermo.Molecules.Incondensable;
    import Thermo.Molecules.Methanol;
    import Thermo.Molecules.Water;
    import Thermo.Molecules.Oxygen;
    import Thermo.moleculeName;
    import Thermo.Phases.Liquid;
    
    import Units.MassTransportCoefficient;
    import Units.MolarFlow;
    import Units.MolarFlux;
    
    annotation (Icon(Rectangle(extent=[-100,60; 100,0], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4,
            fillColor=47,
            rgbfillColor={255,170,85})), Rectangle(extent=[-100,2; 100,-60],
            style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4,
            fillColor=71,
            rgbfillColor={85,170,255}))), Diagram,
      DymolaStoredErrors,
      Documentation(info="<html>
<p>This class implements a DMFC fuel cell from the point of view of reactant flows. A
modelling of the voltage is <em>not</em> included, and must be implemented by child classes.
A temperature port and two ports for methanol concentration on the anode side (inlet and
outlet) are featured.</p>
 
<p>The two inlets and the two outlets are connected to the \"nexus\", an internal protected
(i.e. invisible to the user) object, that accounts for components lost in reactions and
energy that leaves as electric power (=I*V).</p>
 
<p>There are fundamentally three ways by which components can appear or disappear in 
streams:</p>
<ul>
<li>Reaction: anode loses methanol and water, cathode loses oxygen and obtains water; proportional
to current.</li>
<li>Water drag: water leaves the anode and arrives at the cathode; proportional to the reaction, and
therefore to current.</li>
<li>Methanol crossover: anode loses methanol, cathode loses oxygen and obtains water; proportional
to the crossover current.</li>
</ul>
 
<p>The drag factor is taken from Schaffer et al., the temperature dependence is inferred to be
0.025 K<sup>-1</sup> from Ren and Gottesfeld. Although we have a N115 membrane, Ren and Gottesfeld
did not investigate these, so this is a bit of a leap of faith.</p>
 
<p>The methanol diffusion through the membrane is assumed to be constant, as Kallio et al. did not
report a particularly strong dependence; though they did not have many points, they did indeed
cover most of our range of interest.</p>
 
<p>The crossover current is calculated as proportional to the catalyst-layer concentration of
methanol on the anode, the methanol diffusion coefficient in the membrane and the inverse of the 
membrane's thickness. In turn, the difference between bulk and catalyst-layer concentration of 
methanol is proportional to the sum of crossover and reaction current densities.</p>
 
<p>The class calculates some quantities of interest, such as the anodic methanol 
concentration, and the cathodic partial pressures of oxygen and water. Note that all these are
based on the <em>exiting</em> flow.</p>
 
<p>In order to calculate the catalyst-layer methanol anodic concentration, the methanol diffusion 
coefficient in water is required. It is assumed to vary linearly with temperature.</p>
 
<p>Parameters have been taken from Krewer et al., unless differently stated.</p>
 
<h3>References</h3>
<ul>
<li>Ulrike Krewer, Hae-Kwon Yoon, and Hee-Tak Kim: Basic model for membrane electrode assembly 
design for direct methanol fuel cells, Journal of Power Sources, 760-772, 2008.</li>
<li>Thomas Schaffer, Thomas Tschinder, Viktor Hacker, and J&uuml;rgen O. Besenhard: Determination of 
methanol diffusion and electroosmotic drag coefficients in proton-exchange-membranes for DMFC, 
Journal of Power Sources 153(2), 210-216, feb 2006.</li>
<li>Xiaoming Ren, and Shimshon Gottesfeld: Electro-osmotic Drag of Water in Poly(perfluorosulfonic 
acid) Membranes, Journal of the Electrochemical Society 148(1), A87-A93, January 2001.</li>
<li>T. Kallio, K. Kisko, K. Kontturi, R. Serimaa, F. Sundholm, and G. Sundholm: Relationship Between 
Methanol Permeability and Structure of Different Radiation-Grafted Membranes, Fuel Cells: From 
Fundamentals to Systems 4(4), 328-336, December 2004.</li>
</ul>
 
</html>"));
    FlowPort cathode_inlet "The cathode flow's inlet" 
      annotation (extent=[-110,20; -90,40]);
    FlowPort cathode_outlet "The cathode flow's outlet" 
      annotation (extent=[90,20; 110,40]);
    FlowPort anode_inlet "The anode flow's inlet" 
      annotation (extent=[-110,-40; -90,-20]);
    FlowPort anode_outlet "The anode flow's outlet" 
      annotation (extent=[90,-40; 110,-20]);
    PositivePin plus "Pole connected to the cathode" 
                                      annotation (extent=[-70,50; -50,70]);
    NegativePin minus "Pole connected to the anode" 
                                    annotation (extent=[50,50; 70,70]);
  protected 
    FlowTemperature cathodeT "Temperature measurement on the cathode" 
      annotation (extent=[60,20; 80,40]);
    FlowConcentration anodeOutletTC annotation (extent=[60,-40; 80,-20]);
    FlowConcentration anodeInletTC annotation (extent=[-74,-40; -54,-20]);
    SinkPort nexus "Connection of all flows" 
                      annotation (extent=[-32,-10; -12,10]);
    
  public 
    outer parameter Pressure p_env = 101325 "Environment pressure";
    outer parameter Temperature T_env = 298.15 "Enviroment temperature";
    
    parameter Length d_M = 142E-6 "Membrane thickness";
    parameter Area A = 26E-4 "Membrane active area";
    parameter HeatCapacity Cp = 100 "Overall heat capacity of the stack";
    parameter DiffusionCoefficient D_M = 6E-10 
      "Methanol diffusion coefficient in the membrane";
    parameter Boolean enableSanityChecks = true 
      "Whether to activate checks for some non-negative quantities";
    
    // Parameters for N115 membrane.  
    Real k_drag = 4 + 0.025*(T-303.15) "Drag factor for N115";
    MassTransportCoefficient k_ad = 15.6E-6*T/333 "Mass transport coefficient";
    
    Current I = -plus.i "Cell current (generator convention)";
    Voltage V = plus.v - minus.v "Cell voltage";
    CurrentDensity i = I/A "Cell current density";
    
    MolarFlow n_H = I/F "Proton flow through the membrane";
    MolarFlow n_x "Crossover methanol flow";
    MolarFlux N_H = n_H / A "Proton flux";
    MolarFlux N_x = n_x / A "Crossover methanol flux";
    MolarFlow n_drag_h2o = n_H * k_drag "Drag water flow";
    MolarFlux N_drag_h2o = n_drag_h2o / A "Drag water flux";
    
    // KEEP THE INITIAL VALUE, or initialisation will crash on assertion.
    Concentration c_a(start=1000) = anodeOutletTC.c 
      "Methanol concentration, outlet is representative";
    Concentration c_ac(start=100) "Catalyst-layer methanol concentration";
    MoleFraction x_ac "Catalyst-layer methanol molar fraction";
    
    PartialPressure p_o2 "Oxygen partial pressure, outlet is representative";
    PartialPressure p_h2o 
      "Cathodic water partial pressure, outlet is representative";
    TemperatureOutput T "Representative stack temperature" annotation (extent=[100,-8; 120,12]);
    
  protected 
    constant FaradayConstant F = 96485.3415;
    /* This group of vectors represents the coefficients by which
   * proton (*_H) and crossover-methanol (*_M) flows must be 
   * multiplied to  find the associated production terms for all
   * species on cathode and anode; consumption terms are obviously
   * negative. */
    // NOTE remember that methanol reacts to 2 H2O, -3/2 O2, CO2!
    Real[:] cathode_H = {0, 1/2+k_drag, -1/4, 0, 0};
    Real[:] anode_H = {-1/6, -1/6-k_drag, 0, 1/6, 0};
    constant Real[:] cathode_M = {0, 2, -3/2, 1, 0};
    constant Real[:] anode_M = {-1, 0, 0, 0, 0};
    
  equation 
    // Anode-side mass balance, accounting for reaction, drag and crossover
    anode_inlet.n + anode_outlet.n + anode_H*n_H + anode_M*n_x = zeros(size(All,1));
    
    // Cathode-side mass balance, accounting for reaction, drag and crossover
    cathode_inlet.n + cathode_outlet.n + cathode_H*n_H + cathode_M*n_x = zeros(size(All,1));
    
    // The energy "lost" from the heat balance is the electrical power.
    nexus.inlet.H = I*V + der(T)*Cp;
    
    // Definition of oxygen partial pressure. On the denominator, the sum of vapours (methanol and water) and gases (all others).
    p_o2 = p_env * cathodeT.inlet.n[Oxygen] / (sum(cathodeT.vapour) + sum(cathodeT.inlet.n[Incondensable]));
    
    // Definition of water partial pressure (on the cathode side). On the denominator, the sum of vapours (methanol and water) and gases (all others).
    p_h2o = p_env * cathodeT.vapour[Water]  / (sum(cathodeT.vapour) + sum(cathodeT.inlet.n[Incondensable]));
    
    // Relate catalyst-layer methanol molar fraction (in liquid phase) and concentration
    x_ac = c_ac * (x_ac*mw(Methanol)/rho(T,Methanol,Liquid) + (1-x_ac)*mw(Water)/rho(T,Water,Liquid));
    
    // Methanol transport: binds c_a, c_ac and i (N_x is a function of c_ac, N_drag_ch3oh of i).
    k_ad * (c_a-c_ac) = N_x + N_H/6;
    
    // Crossover methanol flux.
    N_x = D_M/d_M * c_ac;
    
    // Set cathode outlet temperature equal to anode outlet temperature
    cathodeT.T = anodeOutletTC.T;
    
    // Charge balance
    plus.i + minus.i = 0;
    
    if enableSanityChecks then
      // Sanity check: crash simulation if conditions are unphysical
      assert( c_ac >= 0, "==> Methanol catalyst-layer concentration is negative ("+String(c_ac)+" mol/m^3) at temperature "+String(T)+" K, bulk concentration "+String(c_a)+" mol/m^3, inlet concentration "+String(anodeInletTC.c)+".");
      
      for i in All loop
        assert( cathode_outlet.n[i] < eps, "==> "+moleculeName(i)+" is entering from the cathode outlet.");
        assert( anode_outlet.n[i] < eps, "==> "+moleculeName(i)+" is entering from the anode outlet.");
        assert( cathode_inlet.n[i] > -eps, "==> "+moleculeName(i)+" is exiting from the cathode inlet.");
        assert( anode_inlet.n[i] > -eps, "==> "+moleculeName(i)+" is exiting from the anode inlet.");
      end for;
    end if;
    
    connect(cathodeT.outlet, cathode_outlet) 
      annotation (points=[78,30; 100,30], style(color=62, rgbcolor={0,127,127}));
    connect(cathode_inlet, nexus.inlet) annotation (points=[-100,30; -46,30; 
          -46,0; -31,0; -31,4.44089e-16],
                                      style(color=62, rgbcolor={0,127,127}));
    connect(cathodeT.inlet, nexus.inlet)       annotation (points=[62,30; -40,
          30; -40,4.44089e-16; -31,4.44089e-16],
                                             style(color=62, rgbcolor={0,127,127}));
    connect(anodeOutletTC.outlet, anode_outlet) annotation (points=[78,-30; 100,
          -30], style(color=62, rgbcolor={0,127,127}));
    connect(anodeOutletTC.inlet, nexus.inlet)    annotation (points=[62,-30; 
          -40,-30; -40,4.44089e-16; -31,4.44089e-16],
                                                  style(color=62, rgbcolor={0,127,
            127}));
    connect(anodeInletTC.inlet, anode_inlet) annotation (points=[-72,-30; -100,
          -30], style(color=62, rgbcolor={0,127,127}));
    connect(anodeInletTC.outlet, nexus.inlet)    annotation (points=[-56,-30; 
          -46,-30; -46,4.44089e-16; -31,4.44089e-16],
                                                  style(color=62, rgbcolor={0,127,
            127}));
    connect(T, cathodeT.T) annotation (points=[110,2; 80,2; 80,23.6; 73.8,23.6], 
        style(color=3, rgbcolor={0,0,255}));
  end FuelCell;
  
  model ConstantVoltageFuelCell "A simplified DMFC with constant voltage" 
    extends FuelCell;
    import Modelica.SIunits.Voltage;
    
    annotation (Documentation(info="<html>
<p>This trivial class inherits from the <tt>FuelCell</tt> class and allows to set a 
constant voltage for the cell.</p>
</html>"));
    
    parameter Voltage V0 = 0.5 "Cell voltage";
    
  equation 
    V = V0;
    
  end ConstantVoltageFuelCell;
  
  model TheveninFuelCell "A DMFC with Thevenin-like voltage" 
    extends FuelCell;
    import Modelica.SIunits.Voltage;
    import Modelica.SIunits.Resistance;
    
    parameter Voltage V0 = 0.7 "Open-circuit voltage";
    parameter Resistance R = 0.005 "Internal resistance";
    
    annotation (Documentation(info="<html>
<p>This class implements a voltage model that emulates a Thevenin equivalent circuit. It is possible
to set the open-circuit voltage and the specific resistance of the cell, from which resistivity and
overall resistance will be calculated (it is trivial to modify the class so that resistivity or
resistance can be set instead).</p>
<p>The default value of the area-specific resistance, 13.3×10<sup>-6</sup>, is valid for a Nafion 
N155 membrane.</p>
<p>Note that the open-circuit value to set is not the one measured on the actual cell, but the one 
that would result by extrapolating the characteristic of the ohmic region to the value of no 
current.</p>
</html>"));
  equation 
    V = V0 - R*I;
  end TheveninFuelCell;
  
  package Test "Package of test cases" 
    model FlowTemperatureTest "A test case for the temperature sensor" 
      
      replaceable FlowTemperature measurement 
                                      annotation (extent=[-20,0; 0,20]);
      annotation (Diagram);
      EnvironmentPort env "Atmospheric air" 
                                      annotation (extent=[-92,20; -72,40]);
      SinkPort sink "Dumpster" 
                        annotation (extent=[40,6; 48,14]);
      MethanolSolution solution "Source of methanol solution" 
                                        annotation (extent=[-80,-10; -68,2]);
      inner parameter Modelica.SIunits.Pressure p_env = 101325;
      inner parameter Modelica.SIunits.Temperature T_env = 298.15;
      inner Real RH_env = time;
      
    equation 
      /* Running from time 0 to 1 will test negative flows and
   * crossing of zero-flow condition. */
      sum(env.outlet.n) = -10*(time-0.5);
      sum(solution.outlet.n) = -(time-0.5);
      
      connect(sink.inlet, measurement.outlet) 
        annotation (points=[40.4,10; 5.55112e-16,10], style(color=62, rgbcolor=
              {0,127,127}));
      connect(env.outlet, measurement.inlet) 
                                           annotation (points=[-73,25; -60,25;
            -60,10; -20,10],   style(color=62, rgbcolor={0,127,127}));
      connect(solution.outlet, measurement.inlet) annotation (points=[-74,-4;
            -60,-4; -60,10; -20,10], style(color=62, rgbcolor={0,127,127}));
    end FlowTemperatureTest;
    
    model FlowConcentrationTest 
      "Test case for the more detailed concentration sensor" 
      extends FlowTemperatureTest(redeclare FlowConcentration measurement);
    end FlowConcentrationTest;
    
    model MixerTest "Test for the mixer unit" 
      inner parameter Modelica.SIunits.Pressure p_env = 101325;
      inner parameter Modelica.SIunits.Temperature T_env = 298.15;
      inner parameter Real RH_env = 60;
      
      Mixer mixer(
        c(fixed=true),
        V(fixed=true),
        T(fixed=true)) 
                  annotation (extent=[-20,0; 0,20]);
      MethanolSolution anodicLoop(T=330) "Solution coming from the anodic loop"
        annotation (extent=[-20,40; 0,60]);
      PureMethanolSource fuelTank "Methanol from the fuel tank" 
        annotation (extent=[0,-40; 20,-20]);
      MethanolSolution condenser(C=0, T=310) 
        "Water recovered from the cathode outlet" 
        annotation (extent=[20,0; 40,20]);
      FlowTemperature flowTemperature annotation (extent=[-64,-26; -44,-6]);
      SinkPort sinkPort annotation (extent=[-28,-20; -20,-12]);
    equation 
      
      annotation (Diagram);
      sum(fuelTank.outlet.n) = -0.1;
      sum(anodicLoop.outlet.n) = -1;
      sum(condenser.outlet.n) = -0.4;
      sum(mixer.outlet.n) = -1.5;
      
      connect(flowTemperature.outlet, sinkPort.inlet) 
        annotation (points=[-46,-16; -27.6,-16], style(color=62, rgbcolor={0,
              127,127}));
      connect(mixer.outlet, flowTemperature.inlet) annotation (points=[-18,10; 
            -72,10; -72,-16; -62,-16], style(color=62, rgbcolor={0,127,127}));
      connect(condenser.outlet, mixer.waterInlet) 
        annotation (points=[30,10; -2,10], style(color=62, rgbcolor={0,127,127}));
      connect(anodicLoop.outlet, mixer.loopInlet) 
        annotation (points=[-10,50; -10,18], style(color=62, rgbcolor={0,127,
              127}));
      connect(mixer.fuelInlet, fuelTank.outlet) 
        annotation (points=[-10,2; -10,-30; 10,-30], style(color=62, rgbcolor={
              0,127,127}));
    end MixerTest;
    
    model SeparatorTest "Test case for the separator unit" 
      
      Separator separator annotation (extent=[-22,-12; 26,38]);
    protected 
      SinkPort liquidSink 
                        annotation (extent=[46,-16; 54,-8]);
      annotation (Diagram);
      SinkPort gasSink   annotation (extent=[46,30; 54,38]);
      EnvironmentPort env             annotation (extent=[-54,28; -34,48]);
      MethanolSolution solution         annotation (extent=[-78,8; -68,18]);
    public 
      inner parameter Modelica.SIunits.Pressure p_env = 101325;
      inner parameter Modelica.SIunits.Temperature T_env = 298.15;
      inner parameter Real RH_env = 60;
      
    equation 
      sum(env.outlet.n) = -1;
      sum(solution.outlet.n) = -2;
      
      connect(env.outlet, separator.inlet)        annotation (points=[-35,33;
            -35,32.5; -22,32.5; -22,13], style(color=62, rgbcolor={0,127,127}));
      connect(solution.outlet, separator.inlet) 
        annotation (points=[-73,13; -22,13], style(color=62, rgbcolor={0,127,
              127}));
      connect(separator.liquidOutlet, liquidSink.inlet) 
        annotation (points=[18.8,3; 18,3; 18,-12; 46.4,-12], style(color=62,
            rgbcolor={0,127,127}));
      connect(separator.gasOutlet, gasSink.inlet)    annotation (points=[18.8,23; 
            18.4,23; 18.4,34; 46.4,34], style(color=62, rgbcolor={0,127,127}));
    end SeparatorTest;
    
    partial model AbstractHeatExchangerTest "Generic test for heat exchangers" 
      
      import Units.MolarFlow;
      
      inner parameter Modelica.SIunits.Pressure p_env = 101325;
      inner parameter Modelica.SIunits.Temperature T_env = 298.15;
      inner parameter Real RH_env = 60;
      
    protected 
      EnvironmentPort env             annotation (extent=[0,-58; 22,-34]);
      SinkPort coldSink "Sink for the cold flow" 
                        annotation (extent=[62,28; 70,36]);
      annotation (Diagram,
        experiment(StopTime=80),
        experimentSetupOutput);
    public 
                     replaceable AbstractHeatExchanger exchanger 
        "The heat exchanger" annotation (extent=[-40,-20; 40,60]);
    protected 
      SinkPort hotSink "Sink for the hot flow" 
                        annotation (extent=[-28,-52; -20,-44]);
      MethanolSolution methanolSolution(T=330) 
        annotation (extent=[-80,22; -60,42]);
      
    public 
      parameter MolarFlow n_air = 1;
      parameter MolarFlow n_sol = 1;
      
    equation 
      sum(env.outlet.n) = -n_air;
      sum(methanolSolution.outlet.n) = -n_sol;
      
      connect(methanolSolution.outlet, exchanger.hot_1) annotation (points=[-70,32; 
            -53,32; -53,32; -36,32],     style(color=62, rgbcolor={0,127,127}));
      connect(hotSink.inlet, exchanger.hot_2) annotation (points=[-27.6,-48; 
            -36,-48; -36,8], style(color=62, rgbcolor={0,127,127}));
      connect(env.outlet, exchanger.cold_2) annotation (points=[20.9,-52; 36,
            -52; 36,8], style(color=62, rgbcolor={0,127,127}));
      connect(exchanger.cold_1, coldSink.inlet) annotation (points=[36,32; 49.2,
            32; 49.2,32; 62.4,32], style(color=62, rgbcolor={0,127,127}));
    end AbstractHeatExchangerTest;
    
    model LMTDHeatExchangerTest "Test for the LMTD-based heat exchanger" 
      extends AbstractHeatExchangerTest(redeclare LMTDHeatExchanger exchanger);
    end LMTDHeatExchangerTest;
    
    model DiscretisedHeatExchangerStepTest 
      "Test for the single-step heat exchanger" 
      extends AbstractHeatExchangerTest(redeclare DiscretisedHeatExchangerStep 
          exchanger);
    end DiscretisedHeatExchangerStepTest;
    
    model DiscretisedHeatExchangerTest 
      "Test for the discretised heat exchanger" 
      extends AbstractHeatExchangerTest(redeclare DiscretisedHeatExchanger 
          exchanger(n=5));
    end DiscretisedHeatExchangerTest;
    
    partial model AbstractCoolerTest "Generic test for air coolers" 
      
      import Modelica.SIunits.VolumeFlowRate;
      
      inner parameter Modelica.SIunits.Pressure p_env = 101325;
      inner parameter Modelica.SIunits.Temperature T_env = 298.15;
      inner parameter Real RH_env = 60;
      
    protected 
      SinkPort sink     annotation (extent=[60,-4; 68,4]);
      annotation (Diagram,
        experiment(StopTime=80),
        experimentSetupOutput);
      MethanolSolution sol(T=330) annotation (extent=[-100,-20; -80,0]);
    public 
      replaceable AbstractCooler cooler annotation (extent=[-20,-20; 20,20]);
      Pump pump annotation (extent=[-60,-20; -40,0]);
      parameter VolumeFlowRate solution = 2E-6/60; // 2 ml/min
      parameter VolumeFlowRate coolingStart = 1E-3/60; // 1 L/min
      parameter VolumeFlowRate coolingStop =  10E-3/60; // 10 L/min
      
    equation 
      pump.V = solution;
      cooler.V_air = coolingStart + time*(coolingStop-coolingStart);
      
      connect(cooler.outlet, sink.inlet) annotation (points=[18.8,1.06581e-15; 
            39.4,1.06581e-15; 39.4,3.88578e-17; 60.4,3.88578e-17], style(color=
              62, rgbcolor={0,127,127}));
      connect(pump.outlet, cooler.inlet) annotation (points=[-50,5.55112e-16; 
            -30,5.55112e-16; -30,1.06581e-15; -18.8,1.06581e-15], style(color=
              62, rgbcolor={0,127,127}));
      connect(sol.outlet, pump.inlet) annotation (points=[-90,-10; -50,-10],
          style(color=62, rgbcolor={0,127,127}));
    end AbstractCoolerTest;
    
    model LMTDCoolerTest "Test for the LMTD-based air cooler" 
      extends AbstractCoolerTest(redeclare LMTDCooler cooler);
      
    end LMTDCoolerTest;
    
    model DiscretisedCoolerTest "Test for the discretised air cooler" 
      extends AbstractCoolerTest(redeclare DiscretisedCooler cooler(exchanger(n=10)));
    end DiscretisedCoolerTest;
    
    partial model CellTest "Generic test suite for fuel-cell models" 
      
      import Modelica.Electrical.Analog.Sources.ConstantCurrent;
      
      inner parameter Modelica.SIunits.Pressure p_env = 101325 
        "Environment pressure";
      inner parameter Modelica.SIunits.Temperature T_env = 298.15 
        "Enviroment temperature";
      inner parameter Real RH_env = 60 "Environment relative humidity";
      
      parameter Modelica.SIunits.VolumeFlowRate anodeFlow = 30E-6/60 
        "Anodic volumetric flow rate";
      parameter Modelica.SIunits.VolumeFlowRate cathodeFlow = 30E-5/60 
        "Cathodic volumetric flow rate";
      parameter Modelica.SIunits.Temperature anodeInletTemperature = 330 
        "Anodic inlet temperature";
      
      replaceable FuelCell fuelCell 
                        annotation (extent=[6,0; 42,34]);
      MethanolSolution methanolSolution(T=anodeInletTemperature) 
                                        annotation (extent=[-66,-30; -54,-18]);
      annotation (Diagram);
      Pump pump "Pump for the anode flow" annotation (extent=[-42,-30; -30,-18]);
      EnvironmentPort air annotation (extent=[-68,18; -48,38]);
      GasFlowController blower annotation (extent=[-42,18; -34,26]);
      SinkPort anodeSink annotation (extent=[62,10; 68,16]);
      SinkPort cathodeSink annotation (extent=[62,18; 68,24]);
      ConstantCurrent I_cell(I=5) annotation (extent=[12,48; 34,72]);
      Modelica.Electrical.Analog.Basic.Ground ground 
        "Negative pole to zero voltage" annotation (extent=[38,40; 58,60]);
    equation 
      pump.V = anodeFlow;
      blower.V = cathodeFlow;
      
      connect(methanolSolution.outlet, pump.inlet) 
                                              annotation (points=[-60,-24; -36,
            -24],        style(color=62, rgbcolor={0,127,127}));
      connect(blower.outlet, fuelCell.cathode_inlet) annotation (points=[-38,26; 
            -18,26; -18,22.1; 6,22.1], style(color=62, rgbcolor={0,127,127}));
      connect(air.outlet, blower.inlet) 
                                   annotation (points=[-49,23; -50,23; -50,22;
            -38,22],        style(color=62, rgbcolor={0,127,127}));
      connect(cathodeSink.inlet, fuelCell.cathode_outlet)    annotation (points=[62.3,21;
            52.15,21; 52.15,22.1; 42,22.1], style(color=62, rgbcolor={0,127,127}));
      connect(anodeSink.inlet, fuelCell.anode_outlet)    annotation (points=[62.3,13; 
            52.15,13; 52.15,11.9; 42,11.9], style(color=62, rgbcolor={0,127,127}));
      connect(ground.p, I_cell.n) 
        annotation (points=[48,60; 34,60], style(color=3, rgbcolor={0,0,255}));
      connect(I_cell.p, fuelCell.plus) annotation (points=[12,60; 12,27.2; 13.2,
            27.2],
          style(color=3, rgbcolor={0,0,255}));
      connect(I_cell.n, fuelCell.minus) annotation (points=[34,60; 34.8,60;
            34.8,27.2], style(color=3, rgbcolor={0,0,255}));
      connect(pump.outlet, fuelCell.anode_inlet) annotation (points=[-36,-18; 
            -16,-18; -16,11.9; 6,11.9], style(color=62, rgbcolor={0,127,127}));
    end CellTest;
    
    model ConstantVoltageCellTest "Test for the constant-voltage model" 
      extends CellTest(redeclare ConstantVoltageFuelCell fuelCell);
      annotation (Diagram);
    end ConstantVoltageCellTest;
    
    model TheveninCellTest "Test for the Thevenin-circuit model" 
      extends CellTest(redeclare TheveninFuelCell fuelCell);
    end TheveninCellTest;
    
    annotation (Documentation(info="<html>
<p>This is the package of test cases for the <tt>Flow</tt> library. Every class
defined in Flow and meant for general use should have a test case. Each test case
should either be a complete or an abstract model, so that checking Flow.Test with F8
is a quick check that the library is Ok; obviously this is not a final check, but
can help catch regressions induced in other classes by some change.</p>
</html>"));
  end Test;
  
end Flow;

