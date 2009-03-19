

package Flow 
  
type MolarFlow = Real(final quantity="Molar flow rate", final unit="mol/s") 
    annotation (Documentation(info="<html>
<p>Just a definition lacking from the standard library.</p>
</html>"));
  
type MolarFlux = Real(final quantity="Molar flux", final unit="mol/(m2.s)") 
    annotation (Documentation(info="<html>
<p>Just a definition lacking from the standard library.</p>
</html>"));
  
type ArealResistance = Real (final quantity="Areal resistance", final unit="Ohm.m2") 
    annotation (Documentation(info="<html>
<p>The resistance per unit area. Notice that the overall resistance is obtained <em>dividing</em>, not
multiplying this unit by area. That is because of the rule of sum of resistances in parallel.</p>
</html>"));
  
type MassTransportCoefficient = Real (final quantity="Mass transport coefficient", final unit
        =                                                                                     "m/s") 
    annotation (Documentation(info="<html>
<p>The coefficient for mass transport, to multiply by concentration difference to obtain the 
molar flux.</p>
</html>"));
  
  connector FlowPort "What passes through a control surface" 
    
    flow MolarFlow[size(Thermo.AllSpecies,1)] n;
    flow Modelica.SIunits.EnthalpyFlowRate H;
    
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
  
  connector TemperaturePort "A connector providing temperature values." 
    
    Modelica.SIunits.Temperature T "Temperature";
    
    annotation (Diagram, Icon(
                        Ellipse(extent=[-100,100; 100,-100],style(
            pattern=0,
            gradient=3,
            fillColor=1,
            rgbfillColor={255,0,0})), Text(
          extent=[-100,158; 100,100],
          string="%name",
          style(color=1, rgbcolor={255,0,0}))),
      Documentation(info="<html>
<p>This is a simple port to make the temperature of an object available for graphic
connections.</p>
</html>"));
  end TemperaturePort;
  
  connector ConcentrationPort "A connector providing concentration values." 
    
    Modelica.SIunits.Concentration c "Methanol concentration";
    
    annotation (Diagram, Icon(
                        Ellipse(extent=[-100,100; 100,-100], style(
            pattern=0,
            gradient=3,
            fillColor=71,
            rgbfillColor={85,170,255})), Text(
          extent=[-100,160; 100,100],
          string="%name",
          style(color=69, rgbcolor={0,128,255}))),
      Documentation(info="<html>
<p>This is a simple port to make the methanol concentration of a flow available
for graphic connections.</p>
</html>"));
  end ConcentrationPort;
  
  model Plug "A class that blocks a flow connection" 
    import Thermo.AllSpecies;
    annotation (Icon(Polygon(points=[-80,60; -80,-60; 100,60; 100,-60; -80,60],
                                                                              style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4,
            fillColor=7,
            rgbfillColor={255,255,255}))), Documentation(info="<html>
<p>This simple item sets the connected flows to zero, effectively \"plugging\"
connectors on other items.</p>
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
    
    parameter Modelica.SIunits.Concentration C = 1000 "Methanol concentration";
    parameter Temperature T = T_env "Temperature";
    
    MoleFraction x_ch3oh "Methanol molar fraction";
    MoleFraction x_h2o(start=1) "Water molar fraction";
    
    FlowPort c "Connection point of the source" 
      annotation (extent=[-20,-20; 20,20]);
  equation 
    assert(C >= 0, "==> Negative concentration given in MethanolSolution object.");
    assert(C <= rho(T,Methanol,LiquidPhase)/mw(Methanol), "==> Methanol concentration over limit (" + String(mw(Methanol)/rho(T,Methanol,LiquidPhase)) + " mol/m�).");
    
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
            pattern=0,
            thickness=4,
            fillColor=1,
            rgbfillColor={255,0,0}))),     Documentation(info="<html>
<p>This item is a source for a pure methanol stream.</p>
</html>"));
  end PureMethanolSource;
  
  model EnvironmentPort "A flow connection to environment conditions." 
    
    import Thermo.h;
    import Thermo.K;
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
    
    outer Real RH_env "Environment relative humidity";
    outer Temperature T_env "Environment temperature";
    outer Pressure p_env "Environment pressure";
    
    MoleFraction y_h2o "Water molar fraction";
    MoleFraction y_o2 "Oxygen molar fraction";
    MoleFraction y_n2 "Nitrogen molar fraction";
    
    MolarEnthalpy h_air "Air molar enthalpy";
    
    FlowPort c   annotation (extent=[-100,-60; -80,-40]);
  protected 
    MolarEnthalpy h_n2 "Nitrogen molar enthalpy";
    MolarEnthalpy h_o2 "Oxygen molar enthalpy";
    MolarEnthalpy h_h2o "Water-vapour molar enthalpy";
    
  equation 
    y_o2 / 0.21 = y_n2 / 0.79; // The O2/N2 ratio.
    y_h2o + y_o2 + y_n2 = 1.0; // Fractions sum to 1.
    y_h2o = RH_env/100 * K(T_env, Water); // Humidity of air
    
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
<p>This object generates a gas flow corresponding to ambient air, including the
effect of humidity.</p>
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
<p>This very simple object is a terminal for flows leaving the system and
about which we do not care much.</p>
</html>"));
    FlowPort flowPort annotation (extent=[-120,-30; -60,30]);
  end SinkPort;
  
  model FlowController "A unit modelling a pump or a MFC" 
    
    import Thermo.mw;
    import Thermo.AllSpecies;
    
    MolarFlow F "Molar flow rate";
    Modelica.SIunits.MassFlowRate m "Mass flow rate";
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
            rgbfillColor={255,255,255})),Text(
          extent=[-100,160; 100,100],
          style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=43,
            rgbfillColor={255,85,85},
            fillPattern=1),
          string="%name")),   Documentation(info="<html>
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
<p>This class implements a mass flow controller with volumetric units (\"field units\").
Since there are two different standards (the actual \"standard\" at 0 Celsius and the \"norm\"
at 70 Fahrenheit), it is necessary to adjust the reference temperature; the default
assumes zero Celsius (so-called \"standard\" value).</p>
<p>The flow assumes that all components are in gas phase and takes their density from
the Thermo library, where the ideal gas law is (usually) assumed.</p>
</html>"));
    import Thermo.rho;
    import Thermo.mw;
    import Thermo.AllSpecies;
    import Thermo.GasPhase;
    import Modelica.SIunits.VolumeFlowRate;
    import Modelica.SIunits.Temperature;
    
    VolumeFlowRate V "Volumetric flow rate";
    parameter Temperature T = 273.15 "Reference temperature";
  equation 
    V = sum({inlet.n[i] * mw(i) / rho(T, i, GasPhase) for i in AllSpecies});
    
  end GasFlowController;
  
  model Pump "A pump, with only liquid phase." 
    extends FlowController;
    annotation (Icon,        Documentation(info="<html>
<p>This class implements a liquid pump with field volumetric units. It will be
necessary to somehow specify its operating temperature to make it work, since
volume changes with temperature: for example, depending on flow direction, it 
might be the temperature of the element upstream or the one downstream.</p>
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
    
    outer Temperature T_env "Environment temperature";
    
    Temperature T(start=T_env) "Temperature of the passing flow.";
    VolumeFlowRate V "Standard volume flow rate.";
    
  equation 
    V = sum({inlet.n[i] * mw(i) / rho(T, i, LiquidPhase) for i in LiquidSpecies});
    // This is to find the temperature.
    inlet.H = (sum(inlet.n[i]*h(T, i, LiquidPhase) for i in LiquidSpecies));
    
  end Pump;
  
    model FlowTemperature "A unit that calculates the temperature of a flow." 
    
      import Modelica.SIunits.MoleFraction;
    
      import Thermo.h;
      import Thermo.AllSpecies;
      import Thermo.GasSpecies;
      import Thermo.LiquidSpecies;
      import Thermo.GasPhase;
      import Thermo.LiquidPhase;
      import Thermo.Water;
      import Thermo.Methanol;
      import Thermo.rr;
      import Thermo.K;
    
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
              fillPattern=1))),
      Documentation(info="<html>
<p>This basic unit takes a flow and returns it unchanged, while actually
performing an equilibrium calculation and figuring out the temperature of
the flow given its associated enthalpic flow.</p>
<p>This unit can be used in different situations to extract the temperature
of a flow, which is also presented with a temperature port in order to
be used in the graphical editor. It can also be used to <em>set</em> the
temperature, provided that the enthalpic flow can be modified in some other
unit.</p>
</html>"));
      FlowPort inlet annotation (extent=[-110,-10; -90,10]);
      FlowPort outlet 
                     annotation (extent=[90,-10; 110,10]);
      TemperaturePort Tm                annotation (extent=[-10,90; 10,110]);
    
      outer parameter Modelica.SIunits.Pressure p_env;
      outer parameter Modelica.SIunits.Temperature T_env;
    
      Modelica.SIunits.Temperature T(start=T_env) = Tm.T "Flow temperature";
      MoleFraction beta "Vapour fraction";
    
      MolarFlow[size(LiquidSpecies,1)] vapour(each min=0) "Vapour flows";
      MolarFlow[size(LiquidSpecies,1)] condensate(each min=0) "Liquid flows";
    
      MoleFraction z_m = inlet.n[Methanol]/sum(inlet.n) 
      "Methanol molar fraction";
      MoleFraction z_w = inlet.n[Water]/sum(inlet.n) "Water molar fraction";
    
    equation 
      if K(T,Water) >= 1 or rr(z_m, z_w, T) >= 1 then
        beta = 1;
      else
        beta = rr(z_m, z_w, T);
      end if;
    
      vapour + condensate = inlet.n[LiquidSpecies];
      vapour = { inlet.n[i]*beta*K(T,i) / (1+beta*(K(T,i) -1)) for i in LiquidSpecies};
    
      // Note that this works only because LiquidSpecies is a vector starting with 1, otherwise the second term would not work.
      inlet.H = sum( inlet.n[i]*h(T, i, GasPhase) for i in GasSpecies) +
                sum( vapour[i]*h(T, i, GasPhase) + condensate[i]*h(T, i, LiquidPhase) for i in LiquidSpecies);
    
      connect(inlet, outlet) 
                          annotation (points=[-100,5.55112e-16; 5,5.55112e-16;
          5,5.55112e-16; 100,5.55112e-16],
                                         style(pattern=0));
    
    end FlowTemperature;
  
    model FlowConcentration 
    "A unit that calculates the methanol concentration of a (liquid) flow." 
    
      import Modelica.SIunits.Concentration;
      import Thermo.AllSpecies;
      import Thermo.LiquidSpecies;
      import Thermo.LiquidPhase;
      import Thermo.Methanol;
      import Thermo.rho;
      import Thermo.mw;
    
      annotation (Diagram, Icon(
          Ellipse(extent=[-100,100; 100,-100], style(
              color=0,
              rgbcolor={0,0,0},
              thickness=4,
              fillColor=7,
              rgbfillColor={255,255,255},
              fillPattern=1)),
        Ellipse(extent=[-36,38; 60,-28], style(
            color=0,
            rgbcolor={0,0,0},
            thickness=4)),
        Rectangle(extent=[-2,56; 82,-40], style(
            pattern=0,
            thickness=2,
            fillColor=7,
            rgbfillColor={255,255,255},
            fillPattern=1))),
      Documentation(info="<html>
<p>This basic unit takes a flow and returns it unchanged, while measuring
temperature and concentration of methanol in the liquid phase of the flow.</p>
<p>If there is no liquid flow, then the reported value is zero.</p>
</html>"));
      FlowPort inlet annotation (extent=[-110,-10; -90,10]);
      FlowPort outlet 
                     annotation (extent=[90,-10; 110,10]);
      FlowTemperature ft annotation (extent=[-10,-10; 10,10]);
      ConcentrationPort cm "Measurement on the methanol concentration" 
      annotation (extent=[-10,-110; 10,-90]);
      TemperaturePort Tm "Measurement on the flow temperature" 
      annotation (extent=[-10,90; 10,110]);
    equation 
      connect(ft.outlet, outlet) annotation (points=[10,6.10623e-16; 68,
          6.10623e-16; 68,5.55112e-16; 100,5.55112e-16], style(
        pattern=0,
        thickness=2,
        fillColor=7,
        rgbfillColor={255,255,255},
        fillPattern=1));
      connect(ft.inlet, inlet) annotation (points=[-10,6.10623e-16; -42,
          6.10623e-16; -42,5.55112e-16; -100,5.55112e-16], style(
        pattern=0,
        thickness=2,
        fillColor=7,
        rgbfillColor={255,255,255},
        fillPattern=1));
      if sum(ft.condensate) > 0.0 or sum(ft.condensate) < 0.0 then
        cm.c = ft.condensate[Methanol] / sum(ft.condensate[i]*mw(i)/rho(ft.Tm.T,i,LiquidPhase) for i in LiquidSpecies);
      else
        cm.c = 0;
      end if;
    
      connect(ft.Tm, Tm) annotation (points=[6.10623e-16,10; 5.55112e-16,10;
          5.55112e-16,100], style(
        pattern=0,
        thickness=2,
        fillColor=7,
        rgbfillColor={255,255,255},
        fillPattern=1));
    end FlowConcentration;
  
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
<p>The separator unit splits a flow in its gaseous and liquid components</p>
<p>It uses a <tt>FlowTemperature</tt> object internally for phase calculations,
and has an additional temperature port connected to the one on the
<tt>FlowTemperature</tt> object so that it is not necessary to place an
additional object down- or upstream to measure the temperature.</p>
</html>"));
    FlowTemperature ft annotation (extent=[-10,-10; 10,10]);
  equation 
    connect(ft.Tm, Tm) annotation (points=[6.10623e-16,10; 5.55112e-16,10;
          5.55112e-16,40], style(color=1, rgbcolor={255,0,0}));
    connect(ft.inlet, inlet) annotation (points=[-10,6.10623e-16; -54,
          6.10623e-16; -54,5.55112e-16; -100,5.55112e-16], style(color=62,
          rgbcolor={0,127,127}));
    connect(ft.outlet, gasOutlet) annotation (points=[10,6.10623e-16; 40,
          6.10623e-16; 40,40; 70,40], style(color=62, rgbcolor={0,127,127}));
    connect(ft.outlet, liquidOutlet) annotation (points=[10,6.10623e-16; 40,
          6.10623e-16; 40,-40; 70,-40], style(color=62, rgbcolor={0,127,127}));
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
<p>This is the most basic form for a heat exchanger. It simply relates the loss in enthalpic
flow to the temperatures of inlet and outlet flows. The enthalpy loss is routed to a protected
(i.e. invisible to the user) sink object.</p>
</html>"));
    TemperaturePort inT "Temperature of the entering flow" 
      annotation (extent=[-92,20; -80,32]);
    TemperaturePort outT "Temperature of the exiting flow" 
      annotation (extent=[80,20; 92,32]);
    Modelica.SIunits.Heat coolingDuty "Removed heat";
    outer Modelica.SIunits.Temperature T_env "Environment temperature";
  protected 
    FlowTemperature flowTemperature_inlet annotation (extent=[-68,-10; -48,10]);
    FlowTemperature flowTemperature_outlet annotation (extent=[50,-10; 70,10]);
    SinkPort sinkPort annotation (extent=[8,-50; 28,-30]);
  equation 
    connect(flowTemperature_inlet.inlet, inlet) annotation (points=[-68,
          6.10623e-16; -76,6.10623e-16; -76,-2.22045e-16; -94,-2.22045e-16],
        style(
        color=62,
        rgbcolor={0,127,127},
        fillColor=7,
        rgbfillColor={255,255,255},
        fillPattern=1));
    connect(inT, flowTemperature_inlet.Tm)              annotation (points=[-86,26;
          -58,26; -58,10], style(
        color=1,
        rgbcolor={255,0,0},
        thickness=2,
        fillColor=7,
        rgbfillColor={255,255,255},
        fillPattern=1));
    connect(flowTemperature_outlet.Tm, outT)              annotation (points=[60,10;
          60,26; 86,26], style(
        color=1,
        rgbcolor={255,0,0},
        thickness=2,
        fillColor=7,
        rgbfillColor={255,255,255},
        fillPattern=1));
    connect(flowTemperature_outlet.outlet, outlet) annotation (points=[70,
          6.10623e-16; 80,6.10623e-16; 80,-2.22045e-16; 94,-2.22045e-16], style(
        color=62,
        rgbcolor={0,127,127},
        fillColor=7,
        rgbfillColor={255,255,255},
        fillPattern=1));
    connect(flowTemperature_inlet.outlet, flowTemperature_outlet.inlet) 
      annotation (points=[-48,6.10623e-16; 2,-3.36456e-22; 2,6.10623e-16; 50,
          6.10623e-16],                                     style(
        color=62,
        rgbcolor={0,127,127},
        fillColor=7,
        rgbfillColor={255,255,255},
        fillPattern=1));
    connect(sinkPort.flowPort, flowTemperature_inlet.outlet) annotation (points=[9,-40; 0,
          -40; 0,6.10623e-16; -48,6.10623e-16], style(
        color=62,
        rgbcolor={0,127,127},
        fillColor=7,
        rgbfillColor={255,255,255},
        fillPattern=1));
    inlet.n + outlet.n = 0*inlet.n;
    
    inlet.H = outlet.H + coolingDuty;
    
  end Cooler;
  
  model Mixer "A unit mixing four molar flows." 
    
    import Modelica.SIunits.Temperature;
    import Modelica.SIunits.AmountOfSubstance;
    import Modelica.SIunits.Concentration;
    import Modelica.SIunits.Volume;
    import Thermo.Methanol;
    import Thermo.GasSpecies;
    import Thermo.LiquidSpecies;
    import Thermo.LiquidPhase;
    import Thermo.h;
    import Thermo.mw;
    import Thermo.rho;
    
    outer parameter Temperature T_env = 298.15 "Environment temperature";
    
    parameter Temperature T_0 = T_env "Initial temperature";
    parameter Volume V_0 = 5E-6 "Initial volume";
    parameter Concentration c_0 = 1000.0 "Initial methanol concentration";
    
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
water and another for methanol inlet, and an output for the anodic loop. It is the 
only unit featuring a full mass and energy balance.</p>
<p>The outlet compositions are the same as the mass balance's molar fractions, 
implying a perfect mixing; the enthalpy flow is also proportional to the internal-energy
holdup.</p>
<p>It is possible to set the initial methanol concentration to some specific value, 
by default it is 1 M.</p>
</html>"));
    
    AmountOfSubstance n[size(Thermo.AllSpecies,1)] "Molar holdup";
    Modelica.SIunits.InternalEnergy U "Internal energy";
    Concentration c(start=c_0, fixed=true) "Methanol concentration";
    Temperature T(start=T_0, fixed=true) "Mixer temperature";
    Volume V(start=V_0, fixed=true) "Solution volume";
    
  equation 
    der(U) = fuelInlet.H + loopInlet.H + waterInlet.H + outlet.H;
    der(n) = fuelInlet.n + loopInlet.n + waterInlet.n + outlet.n;
    
    // Bind outlet's n to composition in holdup
    outlet.n[1:end-1] / sum(outlet.n) = n[1:end-1] / sum(n);
    // Bind outlet's H to specific internal energy and outlet flow
    outlet.H / sum(outlet.n) = U / sum(n);
    
    U = sum(n[i]*h(T,i,LiquidPhase) for i in LiquidSpecies);
    V = sum( n[i]*mw(i)/rho(T, i, LiquidPhase) for i in LiquidSpecies);
    c = n[Methanol] / V;
    
  initial equation 
    n[GasSpecies] = zeros(size(GasSpecies,1));
    
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
    replaceable Electrochemistry.ReactionModelling reaction 
      "Electrochemical model of the cell" annotation (extent=[-10,50; 10,70]);
    import Modelica.SIunits.Length;
    import Modelica.SIunits.DiffusionCoefficient;
    import Modelica.SIunits.Area;
    import Modelica.SIunits.Current;
    import Modelica.SIunits.CurrentDensity;
    import Modelica.SIunits.Concentration;
    import Modelica.SIunits.Temperature;
    import Modelica.SIunits.Voltage;
    import Modelica.SIunits.PartialPressure;
    import Modelica.SIunits.Pressure;
    import Modelica.Constants.eps;
    import Modelica.Electrical.Analog.Interfaces.PositivePin;
    import Modelica.Electrical.Analog.Interfaces.NegativePin;
    
    import Thermo.mw;
    import Thermo.LiquidPhase;
    import Thermo.AllSpecies;
    import Thermo.LiquidSpecies;
    import Thermo.GasSpecies;
    import Thermo.Methanol;
    import Thermo.Water;
    import Thermo.Oxygen;
    
    outer Pressure p_env "Environment pressure";
    
    parameter Length d_M = 142E-6 "Membrane thickness";
    parameter Area A = 26E-4 "Membrane active area";
    parameter DiffusionCoefficient D_M = 6E-10 
      "Methanol diffusion coefficient in the membrane";
    
    // Parameters for N115 membrane.  
    Real k_drag = 4 + 0.025*(T-303.15) "Drag factor for N115";
    MassTransportCoefficient k_ad = 15.6E-6*T/333 "Mass transport coefficient";
    
    Current I = -plus.i "Cell current";
    Voltage V = plus.v - minus.v "Cell voltage";
    CurrentDensity i = I/A "Cell current density";
    
    MolarFlow n_H = I/F "Proton flow through the membrane";
    MolarFlux N_x "Crossover methanol flux";
    MolarFlow n_x = N_x * A "Crossover methanol flow";
    
    Concentration c_a = anodeOutletTC.cm.c 
      "Methanol concentration, outlet is representative";
    Concentration c_ac(start=500) "Catalyst-layer methanol concentration";
    PartialPressure p_o2 "Oxygen partial pressure, outlet is representative";
    PartialPressure p_h2o 
      "Cathodic water partial pressure, outlet is representative";
    Temperature T = cathodeT.Tm.T "Representative stack temperature";
    
  protected 
    constant Modelica.SIunits.FaradayConstant F = 96485.3415;
    /* This group of vectors represents the coefficients by which
   * proton and crossover-methanol flows must be multiplied to 
   * find the flows associated to drag and crossover respectively. */
    Real[:] cathode_D = {0, +k_drag, 0, 0, 0};
    Real[:] anode_D = {0, -k_drag, 0, 0, 0};
    constant Real[:] cathode_M = {0, 2, -3/2, 1, 0};
    constant Real[:] anode_M = {-1, 0, 0, 0, 0};
  protected 
    FlowTemperature cathodeT "Temperature measurement on the cathode" 
      annotation (extent=[60,20; 80,40]);
    FlowConcentration anodeOutletTC annotation (extent=[60,-40; 80,-20]);
    FlowConcentration anodeInletTC annotation (extent=[-74,-40; -54,-20]);
    SinkPort nexus "Connection of all flows" 
                      annotation (extent=[-32,-10; -12,10]);
    
  equation 
    // Anode-side mass balance, accounting for reaction, drag and crossover
    anode_inlet.n + anode_outlet.n + reaction.anode_r + anode_D*n_H + anode_M*n_x = zeros(size(AllSpecies,1));
    
    // Cathode-side mass balance, accounting for reaction, drag and crossover
    cathode_inlet.n + cathode_outlet.n + reaction.cathode_r + cathode_D*n_H + cathode_M*n_x = zeros(size(AllSpecies,1));
    
    // The energy "lost" from the heat balance is the electrical power.
    nexus.flowPort.H = I*V;
    
    // Definition of oxygen partial pressure. On the denominator, the sum of vapours (methanol and water) and gases (all others).
    p_o2 = p_env * cathodeT.inlet.n[Oxygen] / (sum(cathodeT.vapour) + sum(cathodeT.inlet.n[GasSpecies]));
    
    // Definition of water partial pressure (on the cathode side). On the denominator, the sum of vapours (methanol and water) and gases (all others).
    p_h2o = p_env * cathodeT.vapour[Water]  / (sum(cathodeT.vapour) + sum(cathodeT.inlet.n[GasSpecies]));
    
    // Methanol transport: binds c_a, c_ac and i (N_x is a function of c_ac).
    k_ad * (c_a-c_ac) = N_x + i/6/F;
    
    // Crossover methanol flux.
    N_x = D_M/d_M * c_ac;
    
    // Sanity check: crash simulation if conditions are unphysical
    assert( c_ac >= 0, "==> Methanol catalyst-layer concentration is negative ("+String(c_ac)+").");
    assert( max(cathode_outlet.n) < eps, "==> Some components are entering from the cathode outlet.");
    assert( max(anode_outlet.n) < eps, "==> Some components are entering from the anode outlet.");
    
    connect(cathodeT.outlet, cathode_outlet) 
      annotation (points=[80,30; 100,30], style(color=62, rgbcolor={0,127,127}));
    connect(cathode_inlet, nexus.flowPort) 
                                        annotation (points=[-100,30; -46,30; 
          -46,0; -31,0; -31,4.44089e-16],
                                      style(color=62, rgbcolor={0,127,127}));
    connect(cathodeT.inlet, nexus.flowPort)    annotation (points=[60,30; -40,
          30; -40,4.44089e-16; -31,4.44089e-16],
                                             style(color=62, rgbcolor={0,127,127}));
    connect(anodeOutletTC.outlet, anode_outlet) annotation (points=[80,-30; 100,
          -30], style(color=62, rgbcolor={0,127,127}));
    connect(anodeOutletTC.inlet, nexus.flowPort) annotation (points=[60,-30; 
          -40,-30; -40,4.44089e-16; -31,4.44089e-16],
                                                  style(color=62, rgbcolor={0,127,
            127}));
    connect(cathodeT.Tm, anodeOutletTC.Tm) annotation (points=[70,40; 54,40; 54,
          -20; 70,-20],        style(color=1, rgbcolor={255,0,0}));
    connect(anodeInletTC.inlet, anode_inlet) annotation (points=[-74,-30; -100,
          -30], style(color=62, rgbcolor={0,127,127}));
    connect(anodeInletTC.outlet, nexus.flowPort) annotation (points=[-54,-30; 
          -46,-30; -46,4.44089e-16; -31,4.44089e-16],
                                                  style(color=62, rgbcolor={0,127,
            127}));
    connect(reaction.n, minus) 
      annotation (points=[10,60; 60,60], style(color=3, rgbcolor={0,0,255}));
    connect(reaction.p, plus) 
      annotation (points=[-10,60; -60,60], style(color=3, rgbcolor={0,0,255}));
  end FuelCell;
  
  model ConstantVoltageFuelCell "A simplified DMFC with constant voltage" 
    extends FuelCell(redeclare Electrochemistry.ConstantVoltage reaction);
    
    annotation (Documentation(info="<html>
<p>This trivial class inherits from the <tt>FuelCell</tt> class and allows to set a 
constant voltage for the cell.</p>
</html>"));
  end ConstantVoltageFuelCell;
  
  model TheveninFuelCell "A DMFC with Thevenin-like voltage" 
    extends FuelCell(redeclare Electrochemistry.Thevenin reaction);
    
    annotation (Documentation(info="<html>
<p>This class implements a voltage model that emulates a Thevenin equivalent circuit. It is possible
to set the open-circuit voltage and the specific resistance of the cell, from which resistivity and
overall resistance will be calculated (it is trivial to modify the class so that resistivity or
resistance can be set instead).</p>
<p>The default value of the area-specific resistance, 13.3�10<sup>-6</sup>, is valid for a Nafion 
N155 membrane.</p>
<p>Note that the open-circuit value to set is not the one measured on the actual cell, but the one 
that would result by extrapolating the characteristic of the ohmic region to the value of no 
current.</p>
</html>"));
  end TheveninFuelCell;
  
  
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
      // Running from time 0 to 1 will test negative flows too, i.e..
      sum(flowTemp.inlet.n) = time-0.5;
      
      connect(sinkPort.flowPort, flowTemp.outlet) 
        annotation (points=[40.4,10; 5.55112e-16,10], style(color=62, rgbcolor=
              {0,127,127}));
      connect(environmentPort.c, flowTemp.inlet) 
                                           annotation (points=[-49,31; -63.5,31;
            -63.5,10; -20,10], style(color=62, rgbcolor={0,127,127}));
    end FlowTemperatureTest;
    
    model FlowConcentrationTest 
      annotation (Diagram);
      
      SinkPort gasSink   annotation (extent=[46,16; 54,24]);
      EnvironmentPort environmentPort annotation (extent=[-52,46; -32,66]);
      MethanolSolution methanolSolution annotation (extent=[-78,8; -68,18]);
      FlowConcentration flowConcentration annotation (extent=[-6,10; 14,30]);
      
      inner parameter Modelica.SIunits.Pressure p_env = 101325;
      inner parameter Modelica.SIunits.Temperature T_env = 298.15;
      inner parameter Real RH_env = 60;
      
      parameter MolarFlow solutionFlow = 1;
      parameter MolarFlow airFlow = 1;
      
    equation 
      connect(flowConcentration.inlet, methanolSolution.c) annotation (points=[
            -6,20; -40,20; -40,13; -73,13], style(
          color=62,
          rgbcolor={0,127,127},
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
      connect(flowConcentration.inlet, environmentPort.c) annotation (points=[
            -6,20; -56,20; -56,51; -51,51], style(
          color=62,
          rgbcolor={0,127,127},
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
      connect(flowConcentration.outlet, gasSink.flowPort) annotation (points=[14,20;
            46.4,20], style(
          color=62,
          rgbcolor={0,127,127},
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
      
      sum(environmentPort.c.n) = -airFlow;
      sum(methanolSolution.c.n) = -solutionFlow;
      
    end FlowConcentrationTest;
    
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
            -51,32.5; -22,32.5; -22,13], style(color=62, rgbcolor={0,127,127}));
      connect(methanolSolution.c, separator.inlet) 
        annotation (points=[-73,13; -22,13], style(color=62, rgbcolor={0,127,
              127}));
      connect(separator.liquidOutlet, liquidSink.flowPort) 
        annotation (points=[18.8,3; 18,3; 18,-12; 46.4,-12], style(color=62,
            rgbcolor={0,127,127}));
      connect(separator.gasOutlet, gasSink.flowPort) annotation (points=[18.8,23;
            18.4,23; 18.4,34; 46.4,34], style(color=62, rgbcolor={0,127,127}));
    end SeparatorTest;
    
    model CoolerTest 
      
      inner parameter Modelica.SIunits.Pressure p_env = 101325;
      inner parameter Modelica.SIunits.Temperature T_env = 298.15;
      inner parameter Real RH_env = 60;
      
      parameter Modelica.SIunits.Temperature T_set = 320;
      
      EnvironmentPort environmentPort annotation (extent=[-60,20; -40,40]);
      SinkPort sinkPort annotation (extent=[60,-4; 68,4]);
      Cooler cooler annotation (extent=[-24,-20; 20,20]);
      annotation (Diagram);
    equation 
      connect(sinkPort.flowPort, cooler.outlet) annotation (points=[60.4,
            3.88578e-17; 38,3.88578e-17; 38,1.06581e-15; 18.68,1.06581e-15],
          style(
          color=62,
          rgbcolor={0,127,127},
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
      connect(cooler.inlet, environmentPort.c) annotation (points=[-22.68,
            1.06581e-15; -76.05,1.06581e-15; -76.05,25; -59,25], style(
          color=62,
          rgbcolor={0,127,127},
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
      
      sum(environmentPort.c.n) = -1;
      cooler.outT.T = cooler.inT.T - time;
      
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
        annotation (points=[-44,-16; -27.6,-16], style(color=62, rgbcolor={0,
              127,127}));
      connect(mixer.outlet, flowTemperature.inlet) annotation (points=[-18,10;
            -72,10; -72,-16; -64,-16], style(color=62, rgbcolor={0,127,127}));
      connect(condenser.c, mixer.waterInlet) 
        annotation (points=[30,10; -2,10], style(color=62, rgbcolor={0,127,127}));
      connect(anodicLoop.c, mixer.loopInlet) 
        annotation (points=[-10,50; -10,18], style(color=62, rgbcolor={0,127,
              127}));
      connect(mixer.fuelInlet, fuelTank.c) 
        annotation (points=[-10,2; -10,-30; 10,-30], style(color=62, rgbcolor={
              0,127,127}));
    end MixerTest;
    
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
      MethanolSolution methanolSolution annotation (extent=[-66,-30; -54,-18]);
      annotation (Diagram);
      Pump pump "Pump for the anode flow" annotation (extent=[-42,-30; -30,-18]);
      Cooler heater annotation (extent=[-28,2; -8,22]);
      EnvironmentPort air annotation (extent=[-66,26; -46,46]);
      GasFlowController blower annotation (extent=[-42,18; -34,26]);
      SinkPort anodeSink annotation (extent=[62,10; 68,16]);
      SinkPort cathodeSink annotation (extent=[62,18; 68,24]);
      ConstantCurrent I_cell(I=5) annotation (extent=[12,48; 34,72]);
      Modelica.Electrical.Analog.Basic.Ground ground 
        "Negative pole to zero voltage" annotation (extent=[38,40; 58,60]);
    equation 
      connect(methanolSolution.c, pump.inlet) annotation (points=[-60,-24;
            -36.12,-24], style(color=62, rgbcolor={0,127,127}));
      connect(heater.outlet, fuelCell.anode_inlet) annotation (points=[-8.6,12; 
            -1.3,12; -1.3,11.9; 6,11.9], style(color=62, rgbcolor={0,127,127}));
      connect(heater.inlet, pump.outlet) annotation (points=[-27.4,12; -36,12; 
            -36,-18], style(color=62, rgbcolor={0,127,127}));
      connect(blower.outlet, fuelCell.cathode_inlet) annotation (points=[-38,26; 
            -18,26; -18,22.1; 6,22.1], style(color=62, rgbcolor={0,127,127}));
      connect(air.c, blower.inlet) annotation (points=[-65,31; -70.5,31; -70.5,
            22; -38.08,22], style(color=62, rgbcolor={0,127,127}));
      connect(cathodeSink.flowPort, fuelCell.cathode_outlet) annotation (points=[62.3,21;
            52.15,21; 52.15,22.1; 42,22.1], style(color=62, rgbcolor={0,127,127}));
      connect(anodeSink.flowPort, fuelCell.anode_outlet) annotation (points=[62.3,13;
            52.15,13; 52.15,11.9; 42,11.9], style(color=62, rgbcolor={0,127,127}));
      pump.V = anodeFlow;
      heater.outT.T = anodeInletTemperature;
      blower.V = cathodeFlow;
      connect(ground.p, I_cell.n) 
        annotation (points=[48,60; 34,60], style(color=3, rgbcolor={0,0,255}));
      connect(I_cell.p, fuelCell.plus) annotation (points=[12,60; 12,27.2; 13.2,
            27.2],
          style(color=3, rgbcolor={0,0,255}));
      connect(I_cell.n, fuelCell.minus) annotation (points=[34,60; 34.8,60;
            34.8,27.2], style(color=3, rgbcolor={0,0,255}));
    end CellTest;
    
    model ConstantVoltageCellTest "Test for the constant-voltage model" 
      extends CellTest(redeclare ConstantVoltageFuelCell fuelCell, I_cell(I=5));
    end ConstantVoltageCellTest;
    
    model TheveninCellTest "Test for the Thevenin-circuit model" 
      extends CellTest(redeclare TheveninFuelCell fuelCell);
    end TheveninCellTest;
    
  end Test;
  
end Flow;

