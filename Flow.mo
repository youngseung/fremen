

package Flow 
  
type MolarFlowRate = Real(final quantity="MolarFlowRate", final unit="mol/s") 
    annotation (Documentation(info="<html>
<p>Just a definition lacking from the standard library.</p>
</html>"));
  
  connector FlowPort "What passes through a control surface" 
    
    flow MolarFlowRate[size(Thermo.AllSpecies,1)] n;
    flow Modelica.SIunits.EnthalpyFlowRate H;
    
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
</html>"), Icon(Rectangle(extent=[-100,100; 100,-100], style(
            pattern=0,
            fillColor=1,
            rgbfillColor={255,0,0}))));
  end FlowPort;
  
  connector TemperaturePort "A connector providing temperature values." 
    
    Modelica.SIunits.Temperature T "The temperature, in kelvin.";
    
    annotation (Diagram, Icon(
                        Ellipse(extent=[-100,100; 100,-100],style(
            pattern=0,
            gradient=3,
            fillColor=1,
            rgbfillColor={255,0,0}))));
  end TemperaturePort;
  
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
    FlowPort c   annotation (extent=[-90,-10; -70,10]);
  equation 
    c.n = zeros(size(AllSpecies, 1));
    c.H = 0;
    
  end Plug;
  
  annotation (uses(Modelica(version="2.2.1")), Documentation(info="<html>
<p>This package contains various models related to fluid flow in stirred tanks.</p>
</html>"));
  
  model MethanolSolution "A model to use as a generic source" 
    import Thermo.mw;
    import Thermo.rho;
    import Thermo.h;
    import Thermo.GasSpecies;
    import Thermo.Water;
    import Thermo.Methanol;
    import Thermo.GasPhase;
    import Thermo.LiquidPhase;
    import Modelica.SIunits.Temperature;
    import Modelica.SIunits.MoleFraction;
    
    outer Temperature T_env;
    
    parameter Modelica.SIunits.Concentration C = 1000 
      "Concentration of methanol in water, mol/m�.";
    parameter Temperature T = T_env "Source temperature.";
    
    MoleFraction x_ch3oh "Molar fraction of methanol.";
    MoleFraction x_h2o "Molar fraction of water.";
    
    FlowPort c "Connection point of the source" 
      annotation (extent=[-20,-20; 20,20]);
  equation 
    assert(C >= 0, "Negative concentration given in MethanolSolution object.");
    assert(C <= rho(T,Methanol,LiquidPhase)/mw(Methanol), "Methanol concentration over limit (" + String(mw(Methanol)/rho(T,Methanol,LiquidPhase)) + " mol/m�).");
    
    C = x_ch3oh / (x_ch3oh*mw(Methanol)/rho(T,Methanol,LiquidPhase) + x_h2o*mw(Water)/rho(T,Water,LiquidPhase));
    x_ch3oh + x_h2o = 1.0;
    
    c.n[GasSpecies] = zeros(size(GasSpecies, 1));
    c.n[Methanol] / x_ch3oh = c.n[Water] / x_h2o;
    c.H = c.n[Methanol]*h(T,Methanol,LiquidPhase) + c.n[Water]*h(T,Water,LiquidPhase);
    
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
  
  model PureMethanolSource "A model to use as a generic source" 
    import Thermo.h;
    import Thermo.GasSpecies;
    import Thermo.Water;
    import Thermo.Methanol;
    import Thermo.LiquidPhase;
    
    outer Modelica.SIunits.Temperature T_env;
    
    FlowPort c "Connection point of the source" 
      annotation (extent=[-20,-20; 20,20]);
  equation 
    c.n[GasSpecies] = zeros(size(GasSpecies, 1));
    c.n[Water] = 0;
    c.H = c.n[Methanol]*h(T_env,Methanol,LiquidPhase);
    
    annotation (Icon(Ellipse(extent=[-100,100; 100,-100], style(
            color=67,
            rgbcolor={85,255,255},
            thickness=4,
            fillColor=7,
            rgbfillColor={255,255,255}))), Documentation(info="<html>
<p>This item is a source for methanol-water solutions. Parameter <tt>C</tt>
allows to set the concentration in moler per <em>cubic metre</em>; note that
this is 1000 times the normal scale (1M = 1000 mol/m).</p>
</html>"));
  end PureMethanolSource;
  
  model EnvironmentPort "A flow connection to environment conditions." 
    
    import Thermo.h;
    import Thermo.p_vap;
    import Thermo.MolarEnthalpy;
    import Thermo.Methanol;
    import Thermo.Water;
    import Thermo.Oxygen;
    import Thermo.CarbonDioxide;
    import Thermo.Nitrogen;
    import Thermo.GasPhase;
    import Thermo.LiquidPhase;
    import Modelica.SIunits.Temperature;
    import Modelica.SIunits.Pressure;
    import Modelica.SIunits.MoleFraction;
    
    outer Real RH_env "The relative humidity in the laboratory, in percent.";
    outer Temperature T_env "The temperature in the laboratory.";
    outer Pressure p_env "The atmospheric pressure.";
    
    MoleFraction y_h2o "Molar fraction of water in environment air";
    MoleFraction y_o2 "Molar fraction of oxygen in environment air";
    MoleFraction y_n2 "Molar fraction of nitrogen in environment air";
    
    MolarEnthalpy h_air "The molar enthalpy of air in the current conditions.";
    
    FlowPort c   annotation (extent=[-100,-60; -80,-40]);
  protected 
    MolarEnthalpy h_n2 "The molar enthalpy of nitrogen.";
    MolarEnthalpy h_o2 "The molar enthalpy of oxygen.";
    MolarEnthalpy h_h2o "The molar enthalpy of water vapour.";
    
  equation 
    y_o2 / 0.21 = y_n2 / 0.79; // The O2/N2 ratio.
    y_h2o + y_o2 + y_n2 = 1.0; // Fractions sum to 1.
    y_h2o = RH_env/100 * p_vap(T_env, 2)/p_env; // Humidity of air
    
    h_h2o = h(T_env, Water, GasPhase);
    h_o2  = h(T_env, Oxygen, GasPhase);
    h_n2  = h(T_env, Nitrogen, GasPhase);
    h_air = h_h2o*y_h2o + h_o2*y_o2 + h_n2*y_n2;
    
    c.n[Water] / y_h2o = c.n[Oxygen] / y_o2;
    c.n[Nitrogen] / y_n2 = c.n[Oxygen] / y_o2;
    c.n[Methanol] = 0;
    c.n[CarbonDioxide] = 0;
    c.H = h_air*sum(c.n);
    
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
            fillPattern=1))));
    FlowPort flowPort annotation (extent=[-120,-30; -60,30]);
  end SinkPort;
  
  model FlowController "A unit modelling a pump or a MFC" 
    
    import Thermo.mw;
    import Thermo.AllSpecies;
    
    MolarFlowRate F;
    Modelica.SIunits.MassFlowRate m;
    FlowPort inlet "Unit inlet"   annotation (extent=[-12,-10; 8,10]);
    FlowPort outlet "Unit outlet"   annotation (extent=[-10,90; 10,110]);
  equation 
    m = sum({inlet.n[i] * mw(i) for i in AllSpecies});
    F = sum(inlet.n);
    
    connect(inlet, outlet) 
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
    
    VolumeFlowRate V "Volumetric flow rate.";
    parameter Temperature T = 273.15 "Reference temperature for standard flow.";
  equation 
    V = sum({inlet.n[i] * mw(i) / rho(T, i, GasPhase) for i in AllSpecies});
    
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
    
    Temperature T "Temperature of the passing flow.";
    VolumeFlowRate V "Standard volume flow rate.";
    
  equation 
    V = sum({inlet.n[i] * mw(i) / rho(T, i, LiquidPhase) for i in LiquidSpecies});
    // This is to find the temperature.
    inlet.H = (sum(inlet.n[i]*h(T, i, LiquidPhase) for i in LiquidSpecies));
    
  end Pump;
  
    model FlowTemperature "A unit that calculates the temperature of a flow." 
    
      annotation (Diagram, Icon(
          Ellipse(extent=[-100,100; 100,-100], style(
              color=0,
              rgbcolor={0,0,0},
              thickness=4,
              fillColor=7,
              rgbfillColor={255,255,255},
              fillPattern=1)),
          Rectangle(extent=[-38,46; 42,36], style(
              color=0,
              rgbcolor={0,0,0},
              thickness=4,
              fillColor=7,
              rgbfillColor={0,0,0},
              fillPattern=1)),
          Rectangle(extent=[-4,36; 6,-48], style(
              color=0,
              rgbcolor={0,0,0},
              thickness=4,
              fillColor=7,
              rgbfillColor={0,0,0},
              fillPattern=1))));
    
      import Thermo.h;
      import Thermo.p_vap;
      import Thermo.AllSpecies;
      import Thermo.GasSpecies;
      import Thermo.LiquidSpecies;
      import Thermo.GasPhase;
      import Thermo.LiquidPhase;
    
      FlowPort inlet annotation (extent=[-110,-10; -90,10]);
      FlowPort outlet 
                     annotation (extent=[90,-10; 110,10]);
      TemperaturePort Tm                annotation (extent=[-10,90; 10,110]);
    
      Modelica.SIunits.Temperature T(start=298.15) = Tm.T;
    
      outer Modelica.SIunits.Pressure p_env;
    
      MolarFlowRate[size(LiquidSpecies,1)] vapour(each min=0);
      MolarFlowRate[size(LiquidSpecies,1)] condensate(each min=0);
    
      Real beta "The vapour fraction";
    
      Real z_m = inlet.n[1]/sum(inlet.n);
      Real z_w = inlet.n[2]/sum(inlet.n);
    
    equation 
      connect(inlet, outlet) 
                          annotation (points=[-100,5.55112e-16; 5,5.55112e-16;
          5,5.55112e-16; 100,5.55112e-16],
                                         style(pattern=0));
    
      beta = Thermo.rachfordRice(z_m, z_w, T);
      vapour + condensate = inlet.n[LiquidSpecies];
      vapour = { inlet.n[i]*beta*p_vap(T,i)/p_env / (1+beta*(p_vap(T,i)/p_env -1)) for i in LiquidSpecies};
    
      // Note that this works only because LiquidSpecies is a vector starting with 1, otherwise the second term would not work.
      inlet.H = sum( inlet.n[i]*h(T, i, GasPhase) for i in GasSpecies) +
                sum( vapour[i]*h(T, i, GasPhase) + condensate[i]*h(T, i, LiquidPhase) for i in LiquidSpecies);
    
    end FlowTemperature;
  
  model Separator 
    import Thermo.AllSpecies;
    import Thermo.GasSpecies;
    import Thermo.LiquidSpecies;
    import Thermo.LiquidPhase;
    import Thermo.h;
    
    FlowPort inlet annotation (extent=[-110,-10; -90,10]);
    FlowPort gasOutlet annotation (extent=[60,30; 80,50]);
    FlowPort liquidOutlet annotation (extent=[60,-50; 80,-30]);
    TemperaturePort Tm annotation (extent=[-10,30; 10,50]);
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
            thickness=4))), Diagram);
    FlowTemperature ft annotation (extent=[-10,-10; 10,10]);
  equation 
    connect(ft.Tm, Tm) annotation (points=[6.10623e-16,10; 5.55112e-16,10;
          5.55112e-16,40], style(color=1, rgbcolor={255,0,0}));
    connect(ft.inlet, inlet) annotation (points=[-10,6.10623e-16; -54,
          6.10623e-16; -54,5.55112e-16; -100,5.55112e-16], style(pattern=0));
    connect(ft.outlet, gasOutlet) annotation (points=[10,6.10623e-16; 40,
          6.10623e-16; 40,40; 70,40], style(pattern=0));
    connect(ft.outlet, liquidOutlet) annotation (points=[10,6.10623e-16; 40,
          6.10623e-16; 40,-40; 70,-40], style(pattern=0));
    liquidOutlet.n[GasSpecies] = zeros(size(GasSpecies,1));
    
    liquidOutlet.n[LiquidSpecies] = -ft.condensate;
    liquidOutlet.H = sum(h(ft.T, i, LiquidPhase) * liquidOutlet.n[i] for i in LiquidSpecies);
    
  end Separator;
  
  model Cooler "A simplified heat exchanger" 
    
    FlowPort inlet annotation (extent=[-100,-6; -88,6]);
    FlowPort outlet annotation (extent=[88,-6; 100,6]);
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
            fillPattern=1))));
    TemperaturePort inletTemperature "Temperature of the entering flow" 
      annotation (extent=[-92,20; -80,32]);
    TemperaturePort outletTemperature "Temperature of the exiting flow" 
      annotation (extent=[80,20; 92,32]);
    Modelica.SIunits.Heat coolingDuty "The heat removed by the cooler.";
  protected 
    FlowTemperature flowTemperature_inlet annotation (extent=[-68,-10; -48,10]);
    FlowTemperature flowTemperature_outlet annotation (extent=[50,-10; 70,10]);
    SinkPort sinkPort annotation (extent=[6,-52; 26,-32]);
  equation 
    connect(flowTemperature_inlet.inlet, inlet) annotation (points=[-68,
          6.10623e-16; -76,6.10623e-16; -76,-2.22045e-16; -94,-2.22045e-16],
        style(
        pattern=0,
        thickness=4,
        fillColor=7,
        rgbfillColor={255,255,255},
        fillPattern=1));
    connect(inletTemperature, flowTemperature_inlet.Tm) annotation (points=[-86,26;
          -58,26; -58,10],     style(
        pattern=0,
        thickness=4,
        fillColor=7,
        rgbfillColor={255,255,255},
        fillPattern=1));
    connect(flowTemperature_outlet.Tm, outletTemperature) annotation (points=[60,10;
          60,26; 86,26],        style(
        pattern=0,
        thickness=4,
        fillColor=7,
        rgbfillColor={255,255,255},
        fillPattern=1));
    connect(flowTemperature_outlet.outlet, outlet) annotation (points=[70,
          6.10623e-16; 80,6.10623e-16; 80,-2.22045e-16; 94,-2.22045e-16], style(
        pattern=0,
        thickness=4,
        fillColor=7,
        rgbfillColor={255,255,255},
        fillPattern=1));
    
    flowTemperature_inlet.outlet.n + flowTemperature_outlet.inlet.n = 0*flowTemperature_inlet.outlet.n;
    flowTemperature_inlet.outlet.H = flowTemperature_outlet.inlet.H + coolingDuty;
    
    connect(flowTemperature_inlet.outlet, flowTemperature_outlet.inlet) 
      annotation (points=[-48,6.10623e-16; 1,6.10623e-16; 1,6.10623e-16; 50,
          6.10623e-16], style(
        pattern=0,
        thickness=2,
        fillColor=7,
        rgbfillColor={255,255,255},
        fillPattern=1));
    connect(sinkPort.flowPort, flowTemperature_inlet.outlet) annotation (points=[7,-42; 0,
          -42; 0,6.10623e-16; -48,6.10623e-16],          style(
        pattern=0,
        thickness=2,
        fillColor=7,
        rgbfillColor={255,255,255},
        fillPattern=1));
  end Cooler;
  
  model Mixer "A unit mixing four molar flows." 
    
    import Modelica.SIunits.Temperature;
    import Modelica.SIunits.AmountOfSubstance;
    import Thermo.LiquidPhase;
    import Thermo.LiquidSpecies;
    import Thermo.h;
    
    outer Temperature T_env "The environment temperature.";
    
    parameter Temperature T_0 = T_env "The initial temperature.";
    parameter AmountOfSubstance n_MeOH_0 = 0.0 "Initial moles of methanol.";
    parameter AmountOfSubstance n_H2O_0 = 10.0 "Initial moles of water.";
    
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
            fillPattern=1))));
    
    AmountOfSubstance n[size(Thermo.AllSpecies,1)](start={n_MeOH_0,n_H2O_0,0,0,0}) 
      "The molar holdup";
    Modelica.SIunits.InternalEnergy U(start=0) "The mixer's internal energy";
    
  equation 
    der(U) = fuelInlet.H + loopInlet.H + waterInlet.H + outlet.H;
    der(n) = fuelInlet.n + loopInlet.n + waterInlet.n + outlet.n;
    
    // Bind outlet's n to composition in holdup
    outlet.n[1:end-1] / sum(outlet.n) = n[1:end-1] / sum(n);
    // Bind outlet's H to specific internal energy and outlet flow
    outlet.H / sum(outlet.n) = U / sum(n);
    
  initial equation 
    U = sum(n[i]*h(T_0,i,LiquidPhase) for i in LiquidSpecies);
    
  end Mixer;
  
  partial model FuelCell "A generic DMFC" 
    
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
      DymolaStoredErrors);
    import Modelica.SIunits.Length;
    import Modelica.SIunits.DiffusionCoefficient;
    import Modelica.SIunits.Velocity;
    import Modelica.SIunits.Area;
    import Modelica.SIunits.Current;
    import Modelica.SIunits.CurrentDensity;
    import Modelica.SIunits.Concentration;
    import Modelica.SIunits.Temperature;
    import Modelica.SIunits.Voltage;
    import Modelica.SIunits.Pressure;
    import Modelica.SIunits.PartialPressure;
    import Modelica.SIunits.Resistance;
    import Modelica.SIunits.Resistivity;
    
    import Thermo.mw;
    import Thermo.LiquidPhase;
    import Thermo.AllSpecies;
    import Thermo.LiquidSpecies;
    import Thermo.GasSpecies;
    import Thermo.Methanol;
    import Thermo.Water;
    import Thermo.Oxygen;
    
    outer Pressure p_env "The environment pressure";
    
    parameter Length d_M = 40E-6 "Membrane thickness";
    parameter DiffusionCoefficient D_M = 1E-10 
      "Methanol diffusion coefficient in the membrane.";
    parameter Velocity k_ad = 20E-6 "Mass transport coefficient";
    parameter Real k_drag = 2 "Drag factor coefficient";
    parameter Area A = 26E-4 "Membrane area";
    
    FlowPort cathode_inlet "The cathode flow's inlet" 
      annotation (extent=[-110,20; -90,40]);
    FlowPort cathode_outlet "The cathode flow's outlet" 
      annotation (extent=[90,20; 110,40]);
    FlowPort anode_inlet "The anode flow's inlet" 
      annotation (extent=[-110,-40; -90,-20]);
    FlowPort anode_outlet "The cathode flow's outlet" 
      annotation (extent=[90,-40; 110,-20]);
    TemperaturePort Tm "Cell temperature connector" 
                                          annotation (extent=[-10,-10; 10,10]);
    CurrentDensity i "Current density due to reaction";
    Current I = i * A "The overall cell current";
    CurrentDensity i_c "Equivalent current density due to methanol crossover";
    Current I_c = i_c * A "The cell's crossover current";
    Voltage V "Cell voltage";
    
    Concentration c_a "Anodic methanol concentration";
    Concentration c_ac "Catalyst-layer anodic methanol concentration";
    PartialPressure p_o2 "Oxygen partial pressure";
    PartialPressure p_h2o "Water vapour partial pressure in the cathode";
    
  protected 
    FlowTemperature cathodeT "Temperature measurement on the cathode" 
      annotation (extent=[60,20; 80,40]);
    FlowTemperature anodeT "temperature measurement on the anode" 
      annotation (extent=[60,-40; 80,-20]);
    SinkPort nexus "Connection of all flows" 
                      annotation (extent=[-40,-10; -20,10]);
    
  protected 
    constant Modelica.SIunits.FaradayConstant F = 96485.3415;
    /* This group of constant vectors represents the coefficients by which
   * current and crossover current must be multiplied to find the flows
   * associated to reaction, water drag and crossover.
   * Note: the *_drag terms are parameters, not constants, since k_drag
   * can be modified. */
    constant Real[:] cathode_current = {0, 1/2/F, -1/4/F, 0, 0};
    constant Real[:] anode_current = {-1/6/F, -1/6/F, 0, 1/6/F, 0};
    parameter Real[:] cathode_drag = {0, k_drag/F, 0, 0, 0};
    parameter Real[:] anode_drag = {0, -k_drag/F, 0, 0, 0};
    constant Real[:] cathode_crossover = {0, 1/3/F, -1/4/F, 1/6/F, 0};
    constant Real[:] anode_crossover = {-1/6/F, 0, 0, 0, 0};
    
    Temperature T = Tm.T "Cell temperature, a helper variable";
    
  equation 
    connect(anodeT.outlet, anode_outlet) 
      annotation (points=[80,-30; 100,-30], style(color=3, rgbcolor={0,0,255}));
    connect(cathodeT.outlet, cathode_outlet) 
      annotation (points=[80,30; 100,30], style(color=3, rgbcolor={0,0,255}));
    connect(anodeT.Tm, Tm) annotation (points=[70,-20; 70,-16; 46,-16; 46,
          5.55112e-16; 5.55112e-16,5.55112e-16],
                style(color=1, rgbcolor={255,0,0}));
    connect(nexus.flowPort, anode_inlet) 
                                        annotation (points=[-39,4.44089e-16; 
          -50,0; -60,0; -60,-30; -100,-30],
                                     style(pattern=0));
    connect(cathode_inlet, nexus.flowPort) 
                                        annotation (points=[-100,30; -60,30; 
          -60,0; -39,0; -39,4.44089e-16],
                            style(pattern=0));
    connect(cathodeT.inlet, nexus.flowPort)    annotation (points=[60,30; -50,
          30; -50,4.44089e-16; -39,4.44089e-16],
                         style(pattern=0));
    connect(anodeT.inlet, nexus.flowPort)    annotation (points=[60,-30; -50,
          -30; -50,4.44089e-16; -39,4.44089e-16],
                         style(pattern=0));
    connect(cathodeT.Tm, anodeT.Tm) annotation (points=[70,40; 70,44; 46,44; 46,
          -16; 70,-16; 70,-20],
                        style(color=1, rgbcolor={255,0,0}));
    anode_inlet.n + anode_outlet.n + (anode_current+anode_drag)*I + anode_crossover*I_c = zeros(size(AllSpecies,1));
    
    // Anode-side mass balance, accounting for reaction, drag and crossover
    // Cathode-side mass balance, accounting for reaction, drag and crossover
    cathode_inlet.n + cathode_outlet.n + (cathode_current+cathode_drag)*I + cathode_crossover*I_c = zeros(size(AllSpecies,1));
    
    // The energy "lost" from the heat balance is the electrical power.
    nexus.flowPort.H = I*V;
    
    // Definition of anodic methanol concentration: using _outlet_ values
    c_a = anodeT.condensate[Methanol] / sum(anodeT.condensate[i]*mw(i)/Thermo.rho(T,i,LiquidPhase) for i in LiquidSpecies);
    // Definition of oxygen partial pressure: using _outlet_ values
    p_o2 = p_env * cathodeT.inlet.n[Oxygen] / (sum(cathodeT.vapour) + sum(cathodeT.inlet.n[GasSpecies]));
    // Definition of water partial pressure (on the cathode side): using _outlet_ values
    p_h2o = p_env * cathodeT.vapour[Water]  / (sum(cathodeT.vapour) + sum(cathodeT.inlet.n[GasSpecies]));
    
    // Methanol transport: binds c_a, c_ac and i (i_c is a function of c_ac).
    k_ad * (c_a-c_ac) = (i_c + i)/(6*F);
    
    // Equivalent crossover current density in A/m^2.
    i_c = 6*F*(D_M/d_M * c_ac);
    
  end FuelCell;
  
  model ConstantVoltageFuelCell "A simplified DMFC with constant voltage" 
    extends FuelCell;
    
    parameter Modelica.SIunits.Voltage V_cell = 0.5 "Cell voltage";
  equation 
    V = V_cell;
    
  end ConstantVoltageFuelCell;
  
  package Test "Package of test cases" 
    model FlowTemperatureTest "A test case for the temperature sensor" 
      
      FlowTemperature flowTemp        annotation (extent=[-20,0; 0,20]);
      
      annotation (Diagram);
      inner parameter Modelica.SIunits.Pressure p_env = 101325;
      inner parameter Modelica.SIunits.Temperature T_env = 298.15;
      inner Real RH_env = time;
      
      EnvironmentPort environmentPort annotation (extent=[-50,26; -30,46]);
      SinkPort sinkPort annotation (extent=[40,6; 48,14]);
    equation 
      sum(flowTemp.inlet.n) = 1;
      
      connect(sinkPort.flowPort, flowTemp.outlet) 
        annotation (points=[40.4,10; 5.55112e-16,10],
                                                    style(pattern=0));
      connect(environmentPort.c, flowTemp.inlet) 
                                           annotation (points=[-49,31; -63.5,31;
            -63.5,10; -20,10], style(pattern=0));
    end FlowTemperatureTest;
    
    model SeparatorTest 
      Separator separator annotation (extent=[-22,-12; 26,38]);
      SinkPort liquidSink 
                        annotation (extent=[46,-16; 54,-8]);
      annotation (Diagram);
      SinkPort gasSink   annotation (extent=[46,30; 54,38]);
      EnvironmentPort environmentPort annotation (extent=[-52,46; -32,66]);
      
      inner parameter Modelica.SIunits.Pressure p_env = 101325;
      inner parameter Modelica.SIunits.Temperature T_env = 298.15;
      inner parameter Real RH_env = 60;
      MethanolSolution methanolSolution annotation (extent=[-78,8; -68,18]);
    equation 
      
      sum(environmentPort.c.n) = -1;
      sum(methanolSolution.c.n) = -2;
      
      connect(environmentPort.c, separator.inlet) annotation (points=[-51,51;
            -51,32.5; -22,32.5; -22,13], style(pattern=0));
      connect(methanolSolution.c, separator.inlet) 
        annotation (points=[-73,13; -22,13], style(pattern=0));
      connect(separator.liquidOutlet, liquidSink.flowPort) 
        annotation (points=[18.8,3; 18,3; 18,-12; 46.4,-12], style(pattern=0));
      connect(separator.gasOutlet, gasSink.flowPort) annotation (points=[18.8,
            23; 18.4,23; 18.4,34; 46.4,34], style(pattern=0));
    end SeparatorTest;
    
    model CoolerTest 
      
      inner parameter Modelica.SIunits.Pressure p_env = 101325;
      inner parameter Modelica.SIunits.Temperature T_env = 298.15;
      inner parameter Real RH_env = 60;
      
      EnvironmentPort environmentPort annotation (extent=[-60,20; -40,40]);
      SinkPort sinkPort annotation (extent=[60,-4; 68,4]);
      Cooler cooler annotation (extent=[-24,-20; 20,20]);
      annotation (Diagram);
    equation 
      connect(sinkPort.flowPort, cooler.outlet) annotation (points=[60.4,
            3.88578e-17; 38,3.88578e-17; 38,1.06581e-15; 18.68,1.06581e-15],
          style(
          thickness=2,
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
      connect(cooler.inlet, environmentPort.c) annotation (points=[-22.68,
            1.06581e-15; -76.05,1.06581e-15; -76.05,25; -59,25], style(
          thickness=2,
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
      
      sum(environmentPort.c.n) = -1;
      cooler.outletTemperature.T = cooler.inletTemperature.T - time;
      
    end CoolerTest;
    
    model MixerTest 
      inner parameter Modelica.SIunits.Pressure p_env = 101325;
      inner parameter Modelica.SIunits.Temperature T_env = 298.15;
      inner parameter Real RH_env = 60;
      
      Mixer mixer annotation (extent=[-20,0; 0,20]);
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
      sum(fuelTank.c.n) = -0.1;
      
      sum(anodicLoop.c.n) = -1;
      sum(condenser.c.n) = -0.4;
      sum(mixer.outlet.n) = -1.5;
      connect(flowTemperature.outlet, sinkPort.flowPort) 
        annotation (points=[-44,-16; -27.6,-16], style(thickness=2));
      connect(mixer.outlet, flowTemperature.inlet) annotation (points=[-18,10;
            -72,10; -72,-16; -64,-16], style(pattern=0, thickness=2));
      connect(condenser.c, mixer.waterInlet) 
        annotation (points=[30,10; -2,10], style(pattern=0));
      connect(anodicLoop.c, mixer.loopInlet) 
        annotation (points=[-10,50; -10,18], style(pattern=0));
      connect(mixer.fuelInlet, fuelTank.c) 
        annotation (points=[-10,2; -10,-30; 10,-30], style(pattern=0));
    end MixerTest;
    
    model CVCellTest 
      inner parameter Modelica.SIunits.Pressure p_env = 101325;
      inner parameter Modelica.SIunits.Temperature T_env = 298.15;
      inner parameter Real RH_env = 60;
      
      parameter Modelica.SIunits.VolumeFlowRate anodeFlow = 30E-6/60 
        "Anodic volumetric flow rate";
      parameter Modelica.SIunits.VolumeFlowRate cathodeFlow = 30E-5/60 
        "Cathodic volumetric flow rate";
      parameter Modelica.SIunits.Temperature anodeInletTemperature = 330 
        "Anodic inlet temperature";
      parameter Modelica.SIunits.Current I_cell = 0 "Cell current";
      
      ConstantVoltageFuelCell fuelCell 
                        annotation (extent=[6,0; 42,34]);
      MethanolSolution methanolSolution annotation (extent=[-66,-30; -54,-18]);
      annotation (Diagram);
      Pump pump "Pump for the anode flow" annotation (extent=[-42,-30; -30,-18]);
      Cooler heater annotation (extent=[-28,2; -8,22]);
      EnvironmentPort air annotation (extent=[-66,26; -46,46]);
      GasFlowController blower annotation (extent=[-42,18; -34,26]);
      SinkPort anodeSink annotation (extent=[62,10; 68,16]);
      SinkPort cathodeSink annotation (extent=[62,18; 68,24]);
    equation 
      connect(methanolSolution.c, pump.inlet) annotation (points=[-60,-24;
            -36.12,-24], style(color=1, rgbcolor={255,0,0}));
      connect(heater.outlet, fuelCell.anode_inlet) annotation (points=[-8.6,12; 
            -1.3,12; -1.3,11.9; 6,11.9], style(color=1, rgbcolor={255,0,0}));
      connect(heater.inlet, pump.outlet) annotation (points=[-27.4,12; -36,12; 
            -36,-18], style(color=1, rgbcolor={255,0,0}));
      connect(blower.outlet, fuelCell.cathode_inlet) annotation (points=[-38,26; 
            -18,26; -18,22.1; 6,22.1], style(color=1, rgbcolor={255,0,0}));
      connect(air.c, blower.inlet) annotation (points=[-65,31; -70.5,31; -70.5,
            22; -38.08,22], style(color=1, rgbcolor={255,0,0}));
      connect(cathodeSink.flowPort, fuelCell.cathode_outlet) annotation (points=[62.3,21;
            52.15,21; 52.15,22.1; 42,22.1],          style(color=1, rgbcolor={
              255,0,0}));
      connect(anodeSink.flowPort, fuelCell.anode_outlet) annotation (points=[62.3,13; 
            52.15,13; 52.15,11.9; 42,11.9],            style(color=1, rgbcolor=
              {255,0,0}));
      
      fuelCell.I = I_cell;
      pump.V = anodeFlow;
      heater.outletTemperature.T = anodeInletTemperature;
      blower.V = cathodeFlow;
      
    end CVCellTest;
    
  end Test;
  
end Flow;

