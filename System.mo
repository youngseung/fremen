

package System 
  
type MolarFlowRate = Real(final quantity="MolarFlowRate", final unit="mol/s") 
    annotation (Documentation(info="<html>
<p>Just a definition lacking from the standard library.</p>
</html>"));
  
  connector FlowPort "What passes through a control surface" 
    
    flow MolarFlowRate[size(Thermo.AllSpecies,1)] n(each min=0);
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
    import Modelica.SIunits.Temperature;
    import Modelica.SIunits.MoleFraction;
    
    outer parameter Temperature T_env;
    
    parameter Modelica.SIunits.Concentration C = 1000 
      "Concentration of methanol in water, mol/m.";
    parameter Temperature T = T_env "Source temperature.";
    
    MoleFraction x_ch3oh "Molar fraction of methanol.";
    MoleFraction x_h2o "Molar fraction of water.";
    
    FlowPort c "Connection point of the source" 
      annotation (extent=[-20,-20; 20,20]);
  equation 
    assert(C >= 0, "Negative concentration given in MethanolSolution object.");
    assert(C <= rho(T,1,2)/mw(1), "Methanol concentration over limit (" + String(mw(1)/rho(T,1,2)) + " mol/m).");
    
    C = x_ch3oh / (x_ch3oh*mw(1)/rho(T,1,2) + x_h2o*mw(2)/rho(T,2,2));
    x_ch3oh + x_h2o = 1.0;
    
    c.n[GasSpecies] = zeros(size(GasSpecies, 1));
    c.n[Methanol] / x_ch3oh = c.n[Water] / x_h2o;
    c.H = c.n[Methanol]*h(T,1,2) + c.n[Water]*h(T,2,2);
    
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
    import Thermo.dhf;
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
    
    MolarEnthalpy h_n2 "The molar enthalpy of nitrogen.";
    MolarEnthalpy h_o2 "The molar enthalpy of oxygen.";
    MolarEnthalpy h_h2o "The molar enthalpy of water vapour.";
    MolarEnthalpy h_air "The molar enthalpy of air in the current conditions.";
    
    FlowPort c   annotation (extent=[-100,-60; -80,-40]);
  equation 
    y_o2 / 0.21 = y_n2 / 0.79; // The O2/N2 ratio.
    y_h2o + y_o2 + y_n2 = 1.0; // Fractions sum to 1.
    y_h2o = RH_env/100 * p_vap(T_env, 2)/p_env; // Humidity of air
    
    h_h2o = h(T_env, Water, GasPhase) + dhf(2, GasPhase) - dhf(2,LiquidPhase);
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
    annotation (Diagram);
    FlowPort flowPort annotation (extent=[-100,-10; -80,10]);
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
      import Thermo.dhf;
      import Thermo.p_vap;
      import Thermo.AllSpecies;
      import Thermo.GasSpecies;
      import Thermo.LiquidSpecies;
      import Thermo.GasPhase;
      import Thermo.LiquidPhase;
    
      FlowPort inlet annotation (extent=[-110,-10; -90,10]);
      FlowPort outlet 
                     annotation (extent=[90,-10; 110,10]);
      Modelica.SIunits.Temperature T(start=298.15) = Tm.T;
    
      outer Modelica.SIunits.Pressure p_env;
    
      MolarFlowRate[k] vapour(each min=0);
      MolarFlowRate[k] condensate(each min=0);
      TemperaturePort Tm                annotation (extent=[-10,90; 10,110]);
  protected 
      constant Integer k = size(LiquidSpecies,1);
    
    equation 
      vapour + condensate = inlet.n[LiquidSpecies];
    
      if sum(inlet.n[GasSpecies]) <= 0.0 and sum( p_vap(T, i)*inlet.n[i] for i in LiquidSpecies)  < p_env*sum(inlet.n) then
        /* There are no gas species and the liquids' vapour pressure cannot sustain
     * a gas phase (their sum is less than the environment pressure).
     * Therefore, there is only liquid and no gas phase.
     */
        vapour = zeros(k);
      elseif max(inlet.n[i]/p_vap(T,i) for i in LiquidSpecies)  < sum(inlet.n)/p_env then
        /* No species that can be present in liquid phase has a gas fraction sufficient for
     * condensation. Therefore, there is only gas phase.
     */
        condensate = zeros(k);
      else
        // Gas-liquid equilibrium (Raoult's law): y * p_env = x * p_vap
        // FIXME not enough! I must find a way to solve this correctly. I must select the right solution out of n-1.
        vapour / (sum(inlet.n[GasSpecies]) + sum(vapour)) * p_env = { condensate[i]/sum(condensate)*p_vap(T,i) for i in LiquidSpecies};
      end if;
    
      // Note that this works only because LiquidSpecies is a vector starting with 1, otherwise the second term would not work.
      inlet.H = sum( inlet.n[i] * h(T, i, GasPhase) for i in GasSpecies) +
                sum( vapour[i]*h(T, i, GasPhase) + condensate[i]*h(T, i, LiquidPhase) for i in LiquidSpecies);
    
      connect(inlet, outlet) 
                          annotation (points=[-100,5.55112e-16; 5,5.55112e-16; 
          5,5.55112e-16; 100,5.55112e-16],
                                         style(pattern=0));
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
          5.55112e-16,40], style(pattern=0));
    connect(ft.inlet, inlet) annotation (points=[-10,6.10623e-16; -54,
          6.10623e-16; -54,5.55112e-16; -100,5.55112e-16], style(pattern=0));
    connect(ft.outlet, gasOutlet) annotation (points=[10,6.10623e-16; 40,
          6.10623e-16; 40,40; 70,40], style(pattern=0));
    connect(ft.outlet, liquidOutlet) annotation (points=[10,6.10623e-16; 40,
          6.10623e-16; 40,-40; 70,-40], style(pattern=0));
    
    liquidOutlet.n[GasSpecies] = zeros(size(GasSpecies,1));
    
    liquidOutlet.n[LiquidSpecies] = -ft.condensate;
    liquidOutlet.H = sum(h(ft.T, i, LiquidPhase) * ft.condensate[i] for i in LiquidSpecies);
    
  end Separator;
  
  model Mixer "A unit mixing four molar flows." 
    FlowPort port1        annotation (extent=[-90,0; -70,20]);
    FlowPort port2         annotation (extent=[-10,-90; 10,-70]);
    FlowPort port3         annotation (extent=[-10,70; 10,90]);
    FlowPort port4         annotation (extent=[70,0; 90,20]);
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
  equation 
    connect(port1, port3) annotation (points=[-80,10; 0,10; 0,80; 5.55112e-16,
          80], style(pattern=0));
    connect(port1, port2) annotation (points=[-80,10; 0,10; 0,-80; 5.55112e-16,
          -80], style(pattern=0));
    connect(port1, port4) annotation (points=[-80,10; 80,10], style(pattern=0));
  end Mixer;
  
  package Test "Package of test cases" 
    model FlowTemperatureTest "A test case for the temperature sensor" 
      
      FlowTemperature fc              annotation (extent=[-20,0; 0,20]);
      
      annotation (Diagram);
      inner parameter Modelica.SIunits.Pressure p_env = 101325;
      inner parameter Modelica.SIunits.Temperature T_env = 298.15;
      inner parameter Real RH_env = 60;
      
      EnvironmentPort environmentPort annotation (extent=[-50,26; -30,46]);
      SinkPort sinkPort annotation (extent=[40,0; 60,20]);
    equation 
      sum(fc.inlet.n)  = 1;
      
      connect(sinkPort.flowPort, fc.outlet) 
        annotation (points=[41,10; 5.55112e-16,10], style(pattern=0));
      connect(environmentPort.c, fc.inlet) annotation (points=[-49,31; -63.5,31;
            -63.5,10; -20,10], style(pattern=0));
    end FlowTemperatureTest;
    
    model SeparatorTest 
      Separator separator annotation (extent=[-26,-10; 22,40]);
      SinkPort liquidSink 
                        annotation (extent=[62,-22; 82,-2]);
      annotation (Diagram);
      SinkPort gasSink   annotation (extent=[62,34; 82,54]);
      EnvironmentPort environmentPort annotation (extent=[-52,46; -32,66]);
      
      inner parameter Modelica.SIunits.Pressure p_env = 101325;
      inner parameter Modelica.SIunits.Temperature T_env = 298.15;
      inner parameter Real RH_env = 60;
    equation 
      connect(liquidSink.flowPort, separator.liquidOutlet) 
        annotation (points=[63,-12; 42,-12; 42,5; 14.8,5], style(pattern=0));
      connect(gasSink.flowPort, separator.gasOutlet) 
        annotation (points=[63,44; 40,44; 40,25; 14.8,25], style(pattern=0));
      connect(environmentPort.c, separator.inlet) annotation (points=[-51,51;
            -51,14.5; -26,14.5; -26,15], style(pattern=0));
      
      sum( environmentPort.c.n)  = -1;
    end SeparatorTest;
  end Test;
  
end System;
