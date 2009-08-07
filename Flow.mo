                                                                  /**
 * © Federico Zenith, 2008-2009.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */


package Flow 
  
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
  
  package IO "Input and output variables" 
    
  connector ConcentrationInput =Modelica.Blocks.Interfaces.RealInput(redeclare 
          type SignalType=Modelica.SIunits.Concentration) annotation(defaultComponentName="c");
  connector CurrentInput =      Modelica.Blocks.Interfaces.RealInput(redeclare 
          type SignalType=Modelica.SIunits.Current)                   annotation(defaultComponentName="I");
  connector PressureInput =     Modelica.Blocks.Interfaces.RealInput(redeclare 
          type SignalType=Modelica.SIunits.Pressure) 
                                                   annotation(defaultComponentName="p");
  connector TemperatureInput =  Modelica.Blocks.Interfaces.RealInput(redeclare 
          type SignalType=Units.Temperature) annotation(defaultComponentName="T");
  connector VolumeFlowRateInput=Modelica.Blocks.Interfaces.RealInput(redeclare 
          type SignalType=Modelica.SIunits.VolumeFlowRate) 
                                                   annotation(defaultComponentName="V");
  connector ConcentrationOutput=Modelica.Blocks.Interfaces.RealOutput(redeclare 
          type SignalType=Modelica.SIunits.Concentration) annotation(defaultComponentName="c");
  connector PressureOutput =    Modelica.Blocks.Interfaces.RealOutput(redeclare 
          type SignalType=Modelica.SIunits.Pressure) 
                                                   annotation(defaultComponentName="p");
  connector TemperatureOutput = Modelica.Blocks.Interfaces.RealOutput(redeclare 
          type SignalType=Units.Temperature) annotation(defaultComponentName="T");
    
  connector VolumeFlowRateOutput = 
                                Modelica.Blocks.Interfaces.RealOutput(redeclare 
          type SignalType=Modelica.SIunits.VolumeFlowRate) annotation(defaultComponentName="V");
  connector VolumeOutput =      Modelica.Blocks.Interfaces.RealOutput(redeclare 
          type SignalType=Modelica.SIunits.Volume) annotation(defaultComponentName="V");
  end IO;
  
  model Sink "A general-purpose flow sink" 
    
    annotation (Diagram, Icon(
        Rectangle(extent=[-100,100; 100,-100], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=7,
            rgbfillColor={255,255,255},
            fillPattern=1)),
        Line(points=[-40,60; 60,60; 60,-60; -40,-60], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=7,
            rgbfillColor={255,255,255},
            fillPattern=1)),
        Rectangle(extent=[-60,0; 60,0], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=7,
            rgbfillColor={255,255,255},
            fillPattern=1))),
      Documentation(info="<html>
<p>This very simple model is a terminal for flows leaving the system and
about which we do not care much.</p>
</html>"));
    FlowPort inlet "inlet for flow to discard" 
                      annotation (extent=[-120,-30; -60,30]);
  end Sink;
  
  annotation (uses(Modelica(version="2.2.1")), Documentation(info="<html>
<p>This package contains various models related to fluid flow in stirred tanks.</p>
</html>"));
  
  package Sources "Sources of flows of various compositions" 
    
  model Solution "Source of methanol solution with given concentration" 
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
    annotation (Icon(Ellipse(extent=[-100,100; 100,-100], style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=7,
              rgbfillColor={255,255,255}))),
                                           Documentation(info="<html>
<p>This item is a source for methanol-water solutions. Parameter <tt>C</tt>
allows to set the concentration in moler per <em>cubic metre</em>; note that
this is 1000 times the normal scale (1M = 1000 mol/m).</p>
</html>"),
      Diagram);
      
    parameter Concentration C = 1000 "Methanol concentration";
    parameter Temperature T = 298.15 "Temperature";
      
    MoleFraction x_ch3oh "Methanol molar fraction";
    MoleFraction x_h2o "Water molar fraction";
      
  equation 
    assert(C >= 0, "==> Negative concentration given in MethanolSolution object.");
    assert(C <= rho(T,Methanol,Liquid)/mw(Methanol), "==> Methanol concentration over limit (" + String(mw(Methanol)/rho(T,Methanol,Liquid)) + " mol/m3).");
      
    C = x_ch3oh / (x_ch3oh*mw(Methanol)/rho(T,Methanol,Liquid) + x_h2o*mw(Water)/rho(T,Water,Liquid));
    x_ch3oh + x_h2o = 1.0;
      
    outlet.n[Incondensable] = zeros(size(Incondensable, 1));
    outlet.n[Methanol] / x_ch3oh = outlet.n[Water] / x_h2o;
    outlet.H = outlet.n[Methanol]*h(T,Methanol,Liquid) + outlet.n[Water]*h(T,Water,Liquid);
      
  end Solution;
    
  model Methanol "Source of pure methanol" 
      import Thermo.h;
      import Thermo.Molecules.Incondensable;
      import Thermo.Molecules.Water;
      import Thermo.Molecules.Methanol;
      import Thermo.Phases.Liquid;
      import Thermo.mw;
      import Thermo.rho;
      import Units.Temperature;
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
<p>This item is a source for a pure methanol stream. This model also keeps track 
of how much methanol has been released in terms of moles, mass and volume.</p>
</html>"),
      Diagram);
      
    outer Temperature T_env = 298.15 "Enviroment temperature";
      
    AmountOfSubstance n(start=0,fixed=true) 
        "Released methanol moles from simulation start";
    Mass m = n*mw(Methanol) "Released methanol mass";
    Volume V = m / rho(T_env, Methanol, Liquid) "Released methanol volume";
      
  equation 
    der(n) = -outlet.n[Methanol];
      
    outlet.n[Incondensable] = zeros(size(Incondensable, 1));
    outlet.n[Water] = 0;
    outlet.H = outlet.n[Methanol]*h(T_env,Methanol,Liquid);
      
  end Methanol;
    
  model Environment "A flow connection to environment conditions." 
      
    Units.MolarFlow F "Total exchanged molar flow";
    Thermo.Air air "Data about environment air";
      
    FlowPort outlet "Environment-air port" 
                 annotation (extent=[80,-10; 100,10]);
    annotation (defaultComponentName="env", Documentation(info="<html>
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
      
  equation 
    outlet.n = F * air.y;
    outlet.H = F * air.H;
      
  end Environment;
  end Sources;
  
  package Measurements "Measurements on a flow" 
    
  partial model FlowController "A unit modelling a pump or a MFC" 
      
    import Thermo.mw;
    import Thermo.Molecules.All;
    import Modelica.SIunits.MassFlowRate;
    import Units.MolarFlow;
      
    annotation (Icon(
        Ellipse(extent=[-100,100; 100,-100], style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=7,
              rgbfillColor={255,255,255})),
                                         Text(
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
      
    IO.VolumeFlowRateInput V "Volumetric flow rate" 
      annotation (extent=[110,-10; 90,10]);
  equation 
    m = sum({inlet.n[i] * mw(i) for i in All});
    F = sum(inlet.n);
      
    connect(outlet, inlet) annotation (points=[5.55112e-16,100; 5.55112e-16,75;
            5.55112e-16,75; 5.55112e-16,50; 5.55112e-16,5.55112e-16;
            5.55112e-16,5.55112e-16], style(color=62, rgbcolor={0,127,127}));
  end FlowController;
    
  model GasFlowController "A flow controller with only gas phase" 
    extends FlowController;
      
    import Thermo.rho;
    import Thermo.mw;
    import Thermo.Molecules.All;
    import Thermo.Phases.Gas;
    import Modelica.SIunits.Temperature;
      
    parameter Temperature T_ref = 273.15 "Reference temperature";
      
    annotation (defaultComponentName="mfc", Icon,     Documentation(info="<html>
<p>This class implements a mass flow controller with volumetric units (\"field units\").
Since there are at least <em>two</em> different standards (the one named \"standard\" at 
0 Celsius and the one named \"norm\" at 70 Fahrenheit or 21.111... Celsius), it is 
necessary to provide the reference temperature; the default assumes zero Celsius 
(the so-called \"standard\" value).</p>
 
<p>The flow assumes that all components are in gas phase and takes their density from
the Thermo library, where the ideal gas law is (usually) assumed.</p>
</html>"),
      Diagram);
    IO.VolumeFlowRateInput V "Volumetric flow rate" 
      annotation (extent=[110,-10; 90,10]);
      
  equation 
    V = sum({inlet.n[i] * mw(i) / rho(T_ref, i, Gas) for i in All});
      
      connect(outlet, inlet) annotation (points=[5.55112e-16,100; 5.55112e-16,
            75; 5.55112e-16,75; 5.55112e-16,50; 5.55112e-16,5.55112e-16;
            5.55112e-16,5.55112e-16], style(color=62, rgbcolor={0,127,127}));
  end GasFlowController;
    
  model LiquidPump "A pump, only for liquid phase" 
    extends FlowController;
      
      import Thermo.rho;
      import Thermo.mw;
      import Thermo.h;
      import Thermo.Phases.Liquid;
      import Thermo.Phases.Gas;
      import Thermo.Molecules.Condensable;
      import Thermo.Molecules.Incondensable;
      import Modelica.SIunits.VolumeFlowRate;
      import Units.Temperature;
      
    annotation (Icon,        Documentation(info="<html>
<p>This class implements a liquid pump.</p>
<p>The pump takes density values from the Thermo library, and assumes only water and methanol
are present and in liquid phase.</p>
<p>With the appropriate degrees of freedom removed, it can work as a flow measurement as well.</p>
</html>"),
      Diagram);
    IO.VolumeFlowRateInput V "Volumetric flow rate" 
      annotation (extent=[110,-10; 90,10]);
      
    Temperature T;
      
  equation 
    V = sum(inlet.n[i] * mw(i) / rho(T, i, Liquid) for i in Condensable);
    inlet.H = sum(inlet.n[i] * h(T, i, Liquid) for i in Condensable);
      
  end LiquidPump;
    
  model PeristalticPump "A pump for two-phase flow" 
    extends FlowController;
      
    import Thermo.rho;
    import Thermo.mw;
    import Thermo.h;
    import Thermo.Phases.Liquid;
    import Thermo.Phases.Gas;
    import Thermo.Molecules.Condensable;
    import Thermo.Molecules.All;
      
    annotation (Icon,        Documentation(info="<html>
<p>This class implements a peristaltic pump, which can handle two phases.</p>
<p>The pump takes density values from the Thermo library.</p>
<p>With the appropriate degrees of freedom removed, it can work as a flow measurement as well.</p>
</html>"),
      Diagram(Text(
            extent=[-20,60; 20,40],
            string="Deactivated in code",
            style(color=0, rgbcolor={0,0,0}))));
    FlowTemperature T annotation (extent=[-50,32; -30,52], rotation=90);
      
  equation 
    // Deactivate old connection
    // NOTE that inlet is an _outside_ connector, so it needs a minus.
    -inlet.n + T.inlet.n = zeros(size(inlet.n,1));
    -inlet.H + T.inlet.H = 0;
      
    V = sum(T.liquid[i]*mw(i)/rho(T.T,i,Liquid) for i in Condensable) +
        sum(T.vapour[i]*mw(i)/rho(T.T,i,Gas) for i in All);
      
    connect(inlet, T.inlet) annotation (points=[5.55112e-16,5.55112e-16; -20,0;
            -40,0; -40,34], style(color=62, rgbcolor={0,127,127}));
    connect(T.outlet, outlet) annotation (points=[-40,50; -40,100; 5.55112e-16,
            100], style(color=62, rgbcolor={0,127,127}));
  end PeristalticPump;
    
    model FlowTemperature "Calculates a flow's temperature" 
      extends Modelica.Icons.RotationalSensor;
      
      import Modelica.SIunits.MoleFraction;
      import Thermo.Molecules.Incondensable;
      import Thermo.Molecules.Condensable;
      import Thermo.Molecules.All;
      import Thermo.Phases.Gas;
      import Thermo.Phases.Liquid;
      import Thermo.Molecules.Water;
      import Thermo.Molecules.Methanol;
      import Thermo.h;
      import Thermo.rr;
      import Thermo.K;
      import Units.MolarFlow;
      
      annotation (defaultComponentName="T",Diagram, Icon(
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
      IO.TemperatureOutput T "The flow temperature" 
      annotation (extent=[10,-90; -10,-70], rotation=270);
    equation 
      // Liquid and vapour sum to overall flow
      liquid + vapour       = inlet.n;
      // Condensable species are present in both phases
      for i in Condensable loop
        vapour[i] * (1 + beta*(K(T,i) -1)) = inlet.n[i]*beta*K(T,i);
      end for;
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
            0,5.55112e-16; 80,5.55112e-16],style(color=62, rgbcolor={0,127,127}));
    end FlowTemperature;
    
    model FlowConcentration "Calculates a liquid flow's methanol concentration" 
      extends FlowTemperature;
      
      import Thermo.Molecules.Condensable;
      import Thermo.Phases.Liquid;
      import Thermo.Molecules.Methanol;
      import Thermo.rho;
      import Thermo.mw;
      
      annotation (defaultComponentName="TC",Diagram, Icon,
      Documentation(defaultComponentName="c", info="<html>
<p>Adds the ability to read the methanol concentration <em>in the liquid 
phase</em> to the temperature measurement of <tt>FlowTemperature</tt>.</p>
<p>If there is no liquid flow, then the reported value is zero.</p>
</html>"));
      IO.ConcentrationOutput c 
                          annotation (extent=[-10,70; 10,90], rotation=90);
      
    equation 
      // Methanol flow is concentration times volumetric flow
      liquid[Methanol] = c * sum(liquid[i]*mw(i)/rho(T,i,Liquid) for i in Condensable);
      
    end FlowConcentration;
    
    package Test 
      model LiquidPumpTest 
        
        Sources.Solution source 
          annotation (extent=[-100,-10; -80,10]);
        Sink sink annotation (extent=[60,0; 80,20]);
        Flow.Measurements.LiquidPump pump annotation (extent=[-10,-10; 10,10]);
        Measurements.FlowTemperature T annotation (extent=[-40,-10; -20,10]);
        annotation (Diagram);
        
        inner parameter Modelica.SIunits.Pressure p_env = 101325;
        inner parameter Units.Temperature T_env = 298.15;
        inner parameter Units.RelativeHumidity RH_env = 60;
        
        parameter Modelica.SIunits.VolumeFlowRate V = 1;
        
      equation 
        pump.V = V;
        
        connect(pump.outlet, sink.inlet) annotation (points=[6.10623e-16,10;
              30.5,10; 30.5,10; 61,10],
                                   style(color=62, rgbcolor={0,127,127}));
        connect(T.outlet, pump.inlet) annotation (points=[-22,6.10623e-16; -16,
              -3.36456e-22; -16,6.10623e-16; 6.10623e-16,6.10623e-16], style(
              color=62, rgbcolor={0,127,127}));
        connect(T.inlet, source.outlet) annotation (points=[-38,6.10623e-16;
              -70,6.10623e-16; -70,6.66134e-16; -90,6.66134e-16],
                                                              style(color=62,
              rgbcolor={0,127,127}));
      end LiquidPumpTest;
      
      model PeristalticPumpTest 
        
      protected 
        Sources.Solution source(T=320) 
          annotation (extent=[-80,-30; -60,-10]);
        Sink sink annotation (extent=[46,6; 54,14]);
      public 
        Flow.Measurements.PeristalticPump pump 
                                          annotation (extent=[10,-10; 30,10]);
        Measurements.FlowTemperature T annotation (extent=[-20,-10; 0,10]);
        annotation (Diagram);
        
        inner parameter Modelica.SIunits.Pressure p_env = 101325;
        inner parameter Units.Temperature T_env = 298.15;
        inner parameter Units.RelativeHumidity RH_env = 60;
        
        parameter Modelica.SIunits.VolumeFlowRate V = 0.01;
        parameter Real l_to_g_molratio = 0.5;
        
      protected 
        Sources.Environment env annotation (extent=[-80,10; -60,30]);
      equation 
        pump.V = V;
        sum(source.outlet.n) / sum(env.outlet.n) = time;
        
        connect(pump.outlet, sink.inlet) annotation (points=[20,10; 33.2,10; 33.2,
              10; 46.4,10],        style(color=62, rgbcolor={0,127,127}));
        connect(T.outlet, pump.inlet) annotation (points=[-2,6.10623e-16; 4,
              -3.36456e-22; 4,6.10623e-16; 20,6.10623e-16], style(color=62,
              rgbcolor={0,127,127}));
        connect(T.inlet, source.outlet) annotation (points=[-18,6.10623e-16; -40,
              6.10623e-16; -40,-20; -70,-20], style(color=62, rgbcolor={0,127,127}));
        connect(env.outlet, T.inlet) annotation (points=[-61,20; -40,20; -40,
              6.10623e-16; -18,6.10623e-16], style(color=62, rgbcolor={0,127,127}));
      end PeristalticPumpTest;
      
      model FlowTemperatureTest "A test case for the temperature sensor" 
        
        replaceable Measurements.FlowTemperature measurement 
                                        annotation (extent=[-20,0; 0,20]);
        annotation (Diagram);
      protected 
        Flow.Sources.Environment env "Atmospheric air" 
                                        annotation (extent=[-92,20; -72,40]);
        Flow.Sink sink "Dumpster" 
                          annotation (extent=[40,6; 48,14]);
      public 
        Flow.Sources.Solution solution "Source of methanol solution" 
                                          annotation (extent=[-80,-10; -68,2]);
        inner parameter Modelica.SIunits.Pressure p_env = 101325;
        inner parameter Units.Temperature T_env = 298.15;
        inner Units.RelativeHumidity RH_env = time;
        
      equation 
        /* Running from time 0 to 1 will test negative flows and
   * crossing of zero-flow condition. */
        sum(env.outlet.n) = -10*(time-0.5);
        sum(solution.outlet.n) = -(time-0.5);
        
        connect(sink.inlet, measurement.outlet) 
          annotation (points=[40.4,10; -2,10],          style(color=62, rgbcolor=
                {0,127,127}));
        connect(env.outlet, measurement.inlet) 
                                             annotation (points=[-73,30; -60,30;
              -60,10; -18,10],   style(color=62, rgbcolor={0,127,127}));
        connect(solution.outlet, measurement.inlet) annotation (points=[-74,-4;
              -60,-4; -60,10; -18,10], style(color=62, rgbcolor={0,127,127}));
      end FlowTemperatureTest;
      
      model FlowConcentrationTest 
        "Test case for the more detailed concentration sensor" 
        extends Flow.Measurements.Test.FlowTemperatureTest(
                                    redeclare Measurements.FlowConcentration 
            measurement);
      end FlowConcentrationTest;
    end Test;
  end Measurements;
  
  package UnitOperations "Unit operations" 
    model Separator "Splits a flow in two parts" 
      
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
              fillColor=7,
              rgbfillColor={255,255,255},
              fillPattern=1)),
          Ellipse(extent=[60,40; 100,-40], style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=7,
              rgbfillColor={255,255,255},
              fillPattern=1)),
          Rectangle(extent=[-80,40; 80,-40], style(
              color=7,
              rgbcolor={255,255,255},
              fillColor=7,
              rgbfillColor={255,255,255},
              fillPattern=1)),
          Line(points=[-80,40; 80,40], style(color=0, rgbcolor={0,0,0})),
          Line(points=[-80,-40; 80,-40], style(color=0, rgbcolor={0,0,0})),
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
      IO.TemperatureOutput T "Separator temperature" 
                                                  annotation (extent=[100,-10; 120,10]);
    protected 
      Measurements.FlowTemperature ft 
                         annotation (extent=[-10,-10; 10,10]);
      
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
      connect(T, ft.T) annotation (points=[110,5.55112e-16; 78,0; 60,0; 60,-20;
            0,-20; 0,-8; 4.996e-16,-8],
                   style(color=3, rgbcolor={0,0,255}));
    end Separator;
    
    model Burner "An adiabatic combustor" 
      
      import Units.MolarFlow;
      import Thermo.Molecules.Methanol;
      import Thermo.Molecules.Water;
      import Thermo.Molecules.Oxygen;
      import Thermo.Molecules.CarbonDioxide;
      import Thermo.Molecules.Nitrogen;
      
      FlowPort inlet "Burner inlet" annotation (extent=[-108,-10; -88,10]);
      FlowPort outlet "Burner outlet" annotation (extent=[92,-10; 112,10]);
      annotation (
        Diagram,
        Icon(
          Ellipse(extent=[-100,40; -20,-40], style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=7,
              rgbfillColor={255,255,255})),
          Ellipse(extent=[20,40; 100,-40], style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=7,
              rgbfillColor={255,255,255})),
          Rectangle(extent=[-60,40; 60,-40], style(
              color=7,
              rgbcolor={255,255,255},
              thickness=2,
              fillColor=7,
              rgbfillColor={255,255,255},
              fillPattern=1)),
          Rectangle(extent=[-40,40; 40,-40], style(
              pattern=0,
              thickness=2,
              fillColor=45,
              rgbfillColor={255,128,0},
              fillPattern=10)),
          Line(points=[-60,40; 60,40; 60,40], style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=7,
              rgbfillColor={255,255,255},
              fillPattern=1)),
          Line(points=[-60,-40; 60,-40; 60,-40], style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=7,
              rgbfillColor={255,255,255},
              fillPattern=1))),
        Documentation(info="<html>
<p>This unit is a simple reactor that converts all available methanol and oxygen to water
and carbon dioxide. It is smart enough to figure out which of the reactants is limiting,
but not much else.</p>
<p>In particular, no technology is assumed: it could be an actual flame burner, or a
catalytic-bed converter.</p>
</html>"));
      
    protected 
      MolarFlow reaction "Reaction rate for CH3OH+3/2O2 -> 2H2O+CO2";
      Measurements.FlowTemperature T "Flow temperature measurement" 
        annotation (extent=[40,-10; 60,10]);
    public 
      IO.TemperatureOutput T_out(start=373.15) "Temperature after combustion" 
        annotation (extent=[40,-60; 60,-40], rotation=270);
    protected 
      Flow.Sink sink annotation (extent=[0,-50; 20,-30]);
    equation 
      sink.inlet.H = 0;
      sink.inlet.n[Methanol] = reaction;
      sink.inlet.n[Water] = -2*reaction;
      sink.inlet.n[Oxygen] = 1.5*reaction;
      sink.inlet.n[CarbonDioxide] = -reaction;
      sink.inlet.n[Nitrogen] = 0;
      
      reaction = min(inlet.n[Methanol], inlet.n[Oxygen]/1.5);
      
      connect(T_out, T.T) annotation (points=[50,-50; 50,-8], style(color=3,
            rgbcolor={0,0,255}));
      connect(T.outlet, outlet) annotation (points=[58,6.10623e-16; 61,
            6.10623e-16; 61,5.55112e-16; 102,5.55112e-16], style(color=62,
            rgbcolor={0,127,127}));
      connect(T.inlet, inlet) annotation (points=[42,6.10623e-16; -48,
            6.10623e-16; -48,5.55112e-16; -98,5.55112e-16], style(color=62,
            rgbcolor={0,127,127}));
      connect(sink.inlet, inlet) annotation (points=[1,-40; -20,-40; -20,
            5.55112e-16; -98,5.55112e-16], style(color=62, rgbcolor={0,127,127}));
    end Burner;
    
    model Mixer "A unit mixing four molar flows." 
      
      import Modelica.SIunits.Area;
      import Units.Temperature;
      import Modelica.SIunits.AmountOfSubstance;
      import Modelica.SIunits.Concentration;
      import Modelica.SIunits.InternalEnergy;
      import Modelica.SIunits.Volume;
      import Thermo.Molecules.Methanol;
      import Thermo.Molecules.Water;
      import Thermo.Molecules.Incondensable;
      import Thermo.Molecules.Condensable;
      import Thermo.Phases.Liquid;
      import Thermo.h;
      import Thermo.mw;
      import Thermo.rho;
      
      import g = Modelica.Constants.g_n;
      
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
              fillColor=7,
              rgbfillColor={255,255,255},
              fillPattern=1)),
          Line(points=[0,6; 0,-54], style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=7,
              rgbfillColor={255,255,255},
              fillPattern=1)),
          Line(points=[0,6; -52,36], style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=7,
              rgbfillColor={255,255,255},
              fillPattern=1)),
          Line(points=[0,6; 52,36], style(
              color=0,
              rgbcolor={0,0,0},
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
      IO.PressureOutput p 
        "Hydrostatic pressure measured at the bottom of the mixer" 
                       annotation (extent=[-54,-78; -40,-60], rotation=0);
      
      parameter Area A = 50E-4 "Mixer cross-sectional area";
      
      AmountOfSubstance n[size(Thermo.Molecules.All,1)] "Molar holdup";
      InternalEnergy U "Energy holdup";
      Concentration c(start=1000) "Methanol concentration";
      Temperature T "Mixer temperature";
      Volume V(start=5E-6, fixed=true) "Solution volume";
      
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
      
      p = rho(T,Water,Liquid)*g*(V/A);
      
    initial equation 
      n[Incondensable] = zeros(size(Incondensable,1));
      
    end Mixer;
    
    
    package Membrane 
      partial model Equilibrium "A unit mixing four molar flows." 
        
        import Units.MolarEnthalpy;
        import Modelica.SIunits.MoleFraction;
        import Thermo.Molecules.Methanol;
        import Thermo.Molecules.Water;
        import Thermo.Molecules.Oxygen;
        import Thermo.Molecules.CarbonDioxide;
        import Thermo.Molecules.Nitrogen;
        import Thermo.Molecules.Incondensable;
        import Thermo.Molecules.Condensable;
        import Thermo.Molecules.All;
        import Thermo.Phases.Liquid;
        import Thermo.Phases.Gas;
        import Thermo.h;
        import Thermo.K;
        
        Units.Temperature T "Unit temperature";
        
        MolarEnthalpy h_tot = beta*h_g + (1-beta)*h_l "Average molar enthalpy";
        MolarEnthalpy h_g = sum(h(T,i,Gas)*y[i] for i in All) 
          "Molar enthalpy in gas";
        MolarEnthalpy h_l = sum(h(T,i,Liquid)*x[i] for i in Condensable) 
          "Molar enthalpy in liquid";
        
        MoleFraction beta "Fraction of moles in gas phase";
        
        MoleFraction x[size(All,1)] "Liquid molar fraction";
        MoleFraction y[size(All,1)] "Gaseous molar fraction";
        MoleFraction z[size(All,1)] "Overall molar fraction";
        
      equation 
        x[Condensable]   = {z[i]/(1+beta*(K(T, i)-1)) for i in Condensable};
        x[Incondensable] = zeros(size(Incondensable,1));
        
        y[Condensable]   = {K(T,i)*x[i] for i in Condensable};
        y[Incondensable] * beta = z[Incondensable];
        
        if K(T,Water) >= 1 or Thermo.rr(z[Methanol], z[Water], T) >= 1 then
          beta = 1;
        else
          beta = Thermo.rr(z[Methanol], z[Water], T);
        end if;
        
      end Equilibrium;
      
      partial model Balances "A unit mixing four molar flows." 
        
        import Modelica.SIunits.AmountOfSubstance;
        import Thermo.Molecules.All;
        import Thermo.Molecules.Incondensable;
        import Units.MolarFlow;
        
        Modelica.SIunits.InternalEnergy U "Energy holdup";
        AmountOfSubstance n[size(All,1)] "Molar holdup";
        AmountOfSubstance n_tot = sum(n) "Total number of moles";
        MolarFlow F_env = sum(envPort.n) 
          "Mole exchange through the environment port";
        
        FlowPort outlet "The mixer's outlet" 
                              annotation (extent=[-90,-10; -70,10]);
        FlowPort fuelInlet "The methanol-feed inlet" 
                               annotation (extent=[-10,-90; 10,-70]);
        FlowPort inlet "The anodic-loop inlet" 
                               annotation (extent=[70,-10; 90,10]);
        FlowPort envPort "Connection with the environment" 
                              annotation (extent=[-10,70; 10,90]);
        annotation (Diagram, Icon(
            Text(
              extent=[-100,160; 100,100],
              style(
                color=0,
                rgbcolor={0,0,0},
                fillColor=43,
                rgbfillColor={255,85,85},
                fillPattern=1),
              string="%name")));
      protected 
        Sink accumulator "Accumulation sink" 
                                     annotation (extent=[0,-10; 20,10]);
      protected 
        MolarFlow dryIn = sum(inlet.n[Incondensable]) "Inflow of dry gases";
        
      equation 
        der(U) = accumulator.inlet.H;
        der(n) = accumulator.inlet.n;
        
        connect(accumulator.inlet, fuelInlet) 
                                      annotation (points=[1,4.44089e-16; 0,
              4.44089e-16; 0,-80; 5.55112e-16,-80],
                                style(color=62, rgbcolor={0,127,127}));
        connect(accumulator.inlet, outlet) 
                                   annotation (points=[1,4.44089e-16; 0,
              4.44089e-16; 0,5.55112e-16; -80,5.55112e-16],
                                             style(color=62, rgbcolor={0,127,127}));
        connect(accumulator.inlet, inlet) 
                                  annotation (points=[1,4.44089e-16; 0,
              4.44089e-16; 0,0; 80,0; 80,5.55112e-16],
                               style(color=62, rgbcolor={0,127,127}));
        connect(accumulator.inlet, envPort) 
                                    annotation (points=[1,4.44089e-16; 0,
              4.44089e-16; 0,80; 5.55112e-16,80],
                                     style(color=62, rgbcolor={0,127,127}));
      end Balances;
      
      partial model VolumeSum 
        import Modelica.SIunits.Volume;
        
        parameter Volume V = 5E-6 "Total physical volume";
        
        Volume V_g "Gas-phase volume";
        Volume V_l "Liquid-phase volume";
        
      equation 
        V = V_g + V_l;
        
      end VolumeSum;
      
      model Mixer "A unit mixing four molar flows." 
        extends Equilibrium;
        extends Balances;
        extends VolumeSum;
        
        import Units.MolarFlow;
        import Modelica.Constants.eps;
        import Modelica.SIunits.AmountOfSubstance;
        import Modelica.SIunits.Concentration;
        import Modelica.SIunits.EnthalpyFlowRate;
        import Modelica.SIunits.VolumeFraction;
        import Thermo.Molecules.Methanol;
        import Thermo.Molecules.Water;
        import Thermo.Molecules.Oxygen;
        import Thermo.Molecules.CarbonDioxide;
        import Thermo.Molecules.Nitrogen;
        import Thermo.Molecules.Incondensable;
        import Thermo.Molecules.Condensable;
        import Thermo.Molecules.All;
        import Thermo.Phases.Liquid;
        import Thermo.Phases.Gas;
        
        annotation (Diagram, Icon(
            Ellipse(extent=[-80,80; 80,-80], style(
                color=0,
                rgbcolor={0,0,0},
                fillColor=7,
                rgbfillColor={255,255,255},
                fillPattern=10)),
            Rectangle(extent=[-80,0; 80,80], style(
                pattern=0,
                fillColor=7,
                rgbfillColor={255,255,255})),
            Rectangle(extent=[-80,0; 80,80], style(
                color=0,
                rgbcolor={0,0,0},
                pattern=0,
                fillColor=7,
                rgbfillColor={255,255,255},
                fillPattern=10)),
            Line(points=[-80,0; -80,80; 80,80; 80,0], style(
                color=0,
                rgbcolor={0,0,0},
                fillColor=7,
                rgbfillColor={255,255,255},
                fillPattern=1)),
            Line(points=[0,6; 0,-54], style(
                color=0,
                rgbcolor={0,0,0},
                thickness=2,
                fillColor=7,
                rgbfillColor={255,255,255},
                fillPattern=1)),
            Line(points=[0,6; -52,36], style(
                color=0,
                rgbcolor={0,0,0},
                thickness=2,
                fillColor=7,
                rgbfillColor={255,255,255},
                fillPattern=1)),
            Line(points=[0,6; 52,36], style(
                color=0,
                rgbcolor={0,0,0},
                thickness=2,
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
              string="%name")));
        
        parameter Concentration c_0 = 1000 "Initial concentration";
        parameter VolumeFraction liquidFraction = 0.5 
          "Fraction of liquid volume";
        
        AmountOfSubstance n_l_tot = (1-beta)*n_tot 
          "Total moles in liquid phase";
        AmountOfSubstance[:] n_l = x*n_l_tot "Moles in liquid phase";
        
        AmountOfSubstance n_g_tot = beta*n_tot "Total moles in gas phase";
        AmountOfSubstance[:] n_g = y*n_g_tot "Moles in gas phase";
        
        Concentration c(start=c_0,fixed=true) "Methanol concentration";
        
      protected 
        Thermo.Air air "Air composition and enthalpy";
        
        MolarFlow L "Solution flow exchanged with the environment";
        MolarFlow wetOut = -dryIn/(1-sum(y[Condensable])) 
          "Outflow of incoming dry gases, after humidification";
        EnthalpyFlowRate H_wet = sum(Thermo.h(T, i, Gas)*(-inlet.n[i]) for i in Incondensable)
                               + sum(Thermo.h(T, i, Gas)*y[i]*wetOut   for i in Condensable) 
          "Enthalpy of the humidified dry inflow at unit temperature";
        
      equation 
        z = if noEvent( sum(n) > eps) then  n / sum(n) else 0*n;
        c = if noEvent(V_l > eps) then n_l[Methanol] / V_l else 0;
        
        V_g = sum(n_g[i]*Thermo.mw(i)/Thermo.rho(T,i,Gas)    for i in All);
        V_l = sum(n_l[i]*Thermo.mw(i)/Thermo.rho(T,i,Liquid) for i in Condensable);
        
        // Setting main-outlet composition: set only n-1 components, the last is dependent
        outlet.n[2:end] = sum(outlet.n) * z[2:end];
        outlet.H = sum(outlet.n) * h_tot;
        
        if V_g > 0 or L > 0 then
          envPort.n = semiLinear(F_env, air.y, y);
          envPort.H = semiLinear(F_env, air.H, h_g);
          der(L) = 0;
        else
          envPort.n[Incondensable] = -inlet.n[Incondensable];
          envPort.n[Condensable] = y[Condensable]*wetOut + x[Condensable]*L;
          envPort.H = H_wet + h_l*L;
          der(V_g) = 0;
        end if;
        
      initial equation 
        V_l = liquidFraction*V;
        y[Incondensable] = air.y[Incondensable];
        
      end Mixer;
      
      package Test 
        model EquilibriumTest 
          import Modelica.SIunits.MoleFraction;
          import Modelica.SIunits.Temperature;
          
          parameter MoleFraction z_methanol = 0.02;
          parameter MoleFraction z_water = 0.1;
          parameter Temperature T_start = 275;
          parameter Temperature T_stop = 330;
          
          model MyEquilibrium 
            extends Equilibrium;
          end MyEquilibrium;
          
          MyEquilibrium eq;
          
        equation 
          eq.z[1:2] = {z_methanol, z_water};
          der(eq.z[3:4]) = {0,0};
          sum(eq.z) = 1;
          
          eq.T = T_start + time*(T_stop-T_start);
          
        end EquilibriumTest;
        
        model BalancesTest 
          extends Balances;
        end BalancesTest;
        
        model VolumeSumTest 
          import Modelica.SIunits.MoleFraction;
          import Modelica.SIunits.Temperature;
          
          model MyVolumeBalance 
            extends VolumeSum;
          end MyVolumeBalance;
          
          MyVolumeBalance eq;
          
        equation 
          eq.V_g = 0 + time*eq.V;
          
        end VolumeSumTest;

        model MixerTest "Test for the mixer unit" 
          
          import Modelica.SIunits.VolumeFlowRate;
          
          inner parameter Modelica.SIunits.Temperature T_env = 298.15;
          inner parameter Units.RelativeHumidity RH_env = 60;
          
          parameter VolumeFlowRate airFlow = 1E-6;
          parameter VolumeFlowRate solutionIn = 1E-6;
          parameter VolumeFlowRate solutionOut = 1E-6;
          parameter VolumeFlowRate fuel = 0.1E-6;
          
          Mixer mixer annotation (extent=[-20,0; 0,20]);
          Sources.Solution anodicLoop(C=700, T=330) 
            "Solution coming from the anodic loop" 
            annotation (extent=[40,24; 50,34]);
          Sources.Methanol fuelTank "Methanol from the fuel tank" 
            annotation (extent=[4,-36; 16,-24]);
          Sink sinkPort     annotation (extent=[-64,12; -56,20], rotation=180);
          Sources.Environment env annotation (extent=[80,-20; 60,0]);
          Measurements.GasFlowController mfc 
            annotation (extent=[22,-16; 34,-4], rotation=0);
          Measurements.LiquidPump pump_out 
                                 annotation (extent=[-40,4; -28,16]);
          annotation (Diagram, experiment(StopTime=6));
          Measurements.LiquidPump pump_in 
                                 annotation (extent=[24,24; 34,34], rotation=180);
          Measurements.LiquidPump fuel_pump 
                                 annotation (extent=[-16,-36; -4,-24], rotation=0);
          Sink Overflow     annotation (extent=[-14,28; -6,36],  rotation=90);
        equation 
          mfc.V = airFlow;
          pump_in.V = solutionIn;
          pump_out.V = solutionOut;
          fuel_pump.V = fuel;
          
          connect(env.outlet, mfc.inlet) annotation (points=[61,-10; 28,-10], style(
                color=62, rgbcolor={0,127,127}));
          connect(mfc.outlet, mixer.inlet) annotation (points=[28,-4; 14,-4; 14,
                10; -2,10],
                        style(color=62, rgbcolor={0,127,127}));
          connect(pump_out.inlet, mixer.outlet) 
                                            annotation (points=[-34,10; -18,10],
              style(color=62, rgbcolor={0,127,127}));
          connect(sinkPort.inlet, pump_out.outlet) 
                                               annotation (points=[-56.4,16; -34,16],
              style(color=62, rgbcolor={0,127,127}));
          connect(pump_in.inlet, anodicLoop.outlet) annotation (points=[29,29; 45,
                29], style(color=62, rgbcolor={0,127,127}));
          connect(pump_in.outlet, mixer.inlet) annotation (points=[29,24; 14,24; 
                14,10; -2,10],
                            style(color=62, rgbcolor={0,127,127}));
          connect(fuelTank.outlet, fuel_pump.inlet) annotation (points=[10,-30; -10,
                -30], style(color=62, rgbcolor={0,127,127}));
          connect(fuel_pump.outlet, mixer.fuelInlet) annotation (points=[-10,-24; 
                -10,2], style(color=62, rgbcolor={0,127,127}));
          connect(mixer.envPort, Overflow.inlet) annotation (points=[-10,18; 
                -10,28.4], style(color=62, rgbcolor={0,127,127}));
        end MixerTest;
      end Test;
    end Membrane;
    
    package HeatExchangers "Various types of heat exchangers" 
      
    partial model Abstract "An abstract heat exchanger" 
        
      FlowPort hot_1 "Port for hot flow on side 1" 
                     annotation (extent=[-100,20; -80,40]);
      FlowPort hot_2 "Port for hot flow on side 2" 
                      annotation (extent=[-100,-40; -80,-20]);
      annotation (defaultComponentName="exchanger", Diagram, Icon(
          Rectangle(extent=[-100,40; 100,-40], style(
                color=0,
                rgbcolor={0,0,0},
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
          Line(points=[-60,40; -60,-40], style(color=0, rgbcolor={0,0,0})),
          Line(points=[20,40; 20,-40], style(color=0, rgbcolor={0,0,0})),
          Line(points=[-20,40; -20,-40], style(color=0, rgbcolor={0,0,0})),
          Line(points=[60,40; 60,-40], style(color=0, rgbcolor={0,0,0})),
          Line(points=[-40,40; -40,-40], style(color=0, rgbcolor={0,0,0})),
          Line(points=[40,40; 40,-40], style(color=0, rgbcolor={0,0,0})),
          Line(points=[0,40; 0,-40], style(color=0, rgbcolor={0,0,0})),
          Line(points=[80,40; 80,-40], style(color=0, rgbcolor={0,0,0})),
          Line(points=[-80,40; -80,-40], style(color=0, rgbcolor={0,0,0}))),
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
      IO.TemperatureOutput T_hot_1 "Temperature of hot stream on side 1" 
        annotation (extent=[-80,40; -100,60]);
      IO.TemperatureOutput T_hot_2 "Temperature of hot stream on side 2" 
        annotation (extent=[-80,-60; -100,-40]);
      IO.TemperatureOutput T_cold_1 "Temperature of cold stream on side 1" 
        annotation (extent=[80,40; 100,60]);
      IO.TemperatureOutput T_cold_2 "Temperature of cold stream on side 2" 
        annotation (extent=[80,-60; 100,-40]);
        
      parameter Modelica.SIunits.Area A = 2.3E-2 "Effective heat-transfer area";
      parameter Units.HeatTransferCoefficient U = 188 
          "Heat-transfer coefficient";
        
      Modelica.SIunits.HeatFlowRate Q "Heating duty from hot side to cold side";
        
    equation 
      hot_1.H  + hot_2.H  - Q = 0;
      cold_1.H + cold_2.H + Q = 0;
      hot_1.n  + hot_2.n      = 0*hot_1.n;
      cold_1.n + cold_2.n     = 0*hot_1.n;
        
    end Abstract;
      
    model LMTD "A heat exchanger based on the LMTD" 
      extends Flow.UnitOperations.HeatExchangers.Abstract;
        
        import Units.Temperature;
        import Modelica.Constants.eps;
        import Modelica.Math.log;
        import Flow.Measurements.FlowTemperature;
        
      IO.TemperatureOutput LMTD "Log mean temperature difference";
        
      protected 
      Temperature DT1 = T_hot_1 - T_cold_1 "Temperature difference on side 1";
      Temperature DT2 = T_hot_2 - T_cold_2 "Temperature difference on side 2";
        
      protected 
      FlowTemperature hot_1_T "Temperature for hot flow on side 1" 
        annotation (extent=[-60,40; -40,20]);
      FlowTemperature cold_1_T "Temperature for cold flow on side 1" 
        annotation (extent=[40,40; 60,20], rotation=0);
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
      Flow.Sink sinkPort 
                        annotation (extent=[20,-80; 40,-60]);
        
    equation 
      Q = U * A * LMTD;
        
      if noEvent(abs(DT1) > eps and abs(DT2) > eps) then
        LMTD = (DT1-DT2)/log(DT1/DT2);
      else
        LMTD = 0;
      end if;
        
      connect(hot_2_T.outlet, sinkPort.inlet)    annotation (points=[-42,-30; 0,
            -30; 0,-70; 21,-70], style(color=62, rgbcolor={0,127,127}));
      connect(cold_2_T.inlet, sinkPort.inlet)   annotation (points=[42,-30; 0,-30;
            0,-70; 21,-70], style(color=62, rgbcolor={0,127,127}));
      connect(cold_1_T.inlet, sinkPort.inlet)  annotation (points=[42,30; 0,30; 0,
            -70; 21,-70], style(color=62, rgbcolor={0,127,127}));
      connect(hot_1_T.outlet, sinkPort.inlet)   annotation (points=[-42,30; 0,30;
            0,-70; 21,-70], style(color=62, rgbcolor={0,127,127}));
      connect(cold_1, cold_1_T.outlet) 
        annotation (points=[90,30; 70,30; 70,30; 58,30],
                                           style(color=62, rgbcolor={0,127,127}));
      connect(cold_2_T.outlet, cold_2) annotation (points=[58,-30; 90,-30], style(
            color=62, rgbcolor={0,127,127}));
      connect(hot_2_T.inlet, hot_2) annotation (points=[-58,-30; -90,-30], style(
            color=62, rgbcolor={0,127,127}));
      connect(hot_1_T.inlet, hot_1) annotation (points=[-58,30; -90,30], style(
            color=62, rgbcolor={0,127,127}));
      connect(hot_1_T.T, T_hot_1) annotation (points=[-50,38; -50,50; -90,50],
                                         style(color=3, rgbcolor={0,0,255}));
      connect(cold_1_T.T, T_cold_1) annotation (points=[50,38; 50,50; 90,50],
                                  style(color=3, rgbcolor={0,0,255}));
      connect(cold_2_T.T, T_cold_2) annotation (points=[50,-38; 50,-50; 90,-50],
                          style(color=3, rgbcolor={0,0,255}));
      connect(hot_2_T.T, T_hot_2) annotation (points=[-50,-38; -50,-50; -90,-50],
                           style(color=3, rgbcolor={0,0,255}));
    end LMTD;
      
    model DiscretisedStep "A step in a larger, discretised heat exchanger" 
      extends Flow.UnitOperations.HeatExchangers.Abstract;
        
      import Units.Temperature;
      import Flow.Measurements.FlowTemperature;
        
      Temperature T_hot =  (hot_1_T.T  + hot_2_T.T)/2 
          "Average hot-side temperature";
      Temperature T_cold = (cold_1_T.T + cold_2_T.T)/2 
          "Average cold-side temperature";
        
      protected 
      FlowTemperature hot_1_T "Temperature of the hot flow on side 1" 
        annotation (extent=[-60,40; -40,20]);
      FlowTemperature cold_1_T "Temperature of the cold flow on side 1" 
        annotation (extent=[40,40; 60,20]);
      FlowTemperature hot_2_T "Temperature of the hot flow on side 2" 
        annotation (extent=[-60,-40; -40,-20]);
      FlowTemperature cold_2_T "Temperature of the cold flow on side 2" 
        annotation (extent=[40,-40; 60,-20]);
      annotation (defaultComponentName="step", Diagram, Documentation(info="<html>
<p>To implement a discretised heat exchanger, a single step is implemented 
here as a simplistic heat exchanger. It could be based on LMTD, but using
the average temperature is numerically more robust and allows Dymola to
perform algebraic manipulation.</p>
</html>"));
      Flow.Sink sink "Sink element" 
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
      connect(cold_1_T.T, T_cold_1) annotation (points=[50,38; 50,50; 90,50],
                                    style(color=3, rgbcolor={0,0,255}));
      connect(cold_2_T.T, T_cold_2) annotation (points=[50,-38; 50,-50; 90,-50],
                                             style(color=3, rgbcolor={0,0,255}));
      connect(hot_2_T.T, T_hot_2) annotation (points=[-50,-38; -50,-50; -90,-50],
                                                style(color=3, rgbcolor={0,0,255}));
      connect(hot_1_T.T, T_hot_1) annotation (points=[-50,38; -50,50; -90,50],
                             style(color=3, rgbcolor={0,0,255}));
    end DiscretisedStep;
      
    model Discretised "A heat exchanger made up of many discrete sections" 
      extends Flow.UnitOperations.HeatExchangers.Abstract;
        
      parameter Integer n(min=3) = 10 "Number of discretisation units";
        
      protected 
      Flow.UnitOperations.HeatExchangers.DiscretisedStep[n] steps(each A=A/n,
            each U=U) "The steps the exchanger is divided in";
        
      protected 
      Flow.Sink sink 
                    annotation (extent=[20,-80; 40,-60]);
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
        
    end Discretised;
      
      package Test 
        partial model AbstractHeatExchangerTest 
          "Generic test for heat exchangers" 
          
          import Modelica.SIunits.VolumeFlowRate;
          
          inner parameter Modelica.SIunits.Pressure p_env = 101325;
          inner parameter Units.Temperature T_env = 298.15;
          inner parameter Units.RelativeHumidity RH_env = 60;
          
        protected 
          Sources.Environment env         annotation (extent=[100,-60; 80,-40]);
          Sink coldSink "Sink for the cold flow" 
                            annotation (extent=[62,28; 70,36]);
          annotation (Diagram,
            experiment(StopTime=80),
            experimentSetupOutput);
        public 
          replaceable Abstract exchanger "The heat exchanger" 
                                 annotation (extent=[-40,-20; 40,60]);
        protected 
          Sink hotSink "Sink for the hot flow" 
                            annotation (extent=[-28,-52; -20,-44]);
          Sources.Solution methanolSolution(   T=330) 
            annotation (extent=[-80,-40; -60,-20]);
          Measurements.LiquidPump pump 
                    annotation (extent=[-80,0; -60,20]);
        public 
          Measurements.GasFlowController mfc 
                                annotation (extent=[52,-42; 72,-22]);
          
        public 
          parameter VolumeFlowRate air = 100E-3/60 "Full scale of cooling air";
          parameter VolumeFlowRate sol = 10E-6/60 "Loop solution";
          
        equation 
          mfc.V = 0.01*air+0.99*air*time;
          pump.V = sol;
          
          connect(hotSink.inlet, exchanger.hot_2) annotation (points=[-27.6,-48;
                -36,-48; -36,8], style(color=62, rgbcolor={0,127,127}));
          connect(exchanger.cold_1, coldSink.inlet) annotation (points=[36,32;
                49.2,32; 49.2,32; 62.4,32],
                                       style(color=62, rgbcolor={0,127,127}));
          connect(env.outlet, mfc.inlet) annotation (points=[81,-50; 62,-50; 62,-32],
              style(color=62, rgbcolor={0,127,127}));
          connect(mfc.outlet, exchanger.cold_2) annotation (points=[62,-22; 62,
                8; 36,8],
                       style(color=62, rgbcolor={0,127,127}));
          connect(pump.outlet, exchanger.hot_1) annotation (points=[-70,20; -70,
                32; -36,32],
                         style(color=62, rgbcolor={0,127,127}));
          connect(methanolSolution.outlet, pump.inlet) annotation (points=[-70,-30;
                -70,10], style(color=62, rgbcolor={0,127,127}));
        end AbstractHeatExchangerTest;
        
        model LMTDHeatExchangerTest "Test for the LMTD-based heat exchanger" 
          extends AbstractHeatExchangerTest(redeclare LMTD exchanger);
        end LMTDHeatExchangerTest;
        
        model DiscretisedHeatExchangerStepTest 
          "Test for the single-step heat exchanger" 
          extends AbstractHeatExchangerTest(redeclare DiscretisedStep exchanger);
        end DiscretisedHeatExchangerStepTest;
        
        model DiscretisedHeatExchangerTest 
          "Test for the discretised heat exchanger" 
          extends AbstractHeatExchangerTest(redeclare Discretised exchanger);
        end DiscretisedHeatExchangerTest;
      end Test;
    end HeatExchangers;
    
    package Coolers "Various types of coolers" 
      
    partial model Abstract "An abstract cooler for a process stream" 
        
      FlowPort inlet "Inlet to the cooler" 
                     annotation (extent=[-100,-6; -88,6]);
      FlowPort outlet "Outlet from the cooler" 
                      annotation (extent=[88,-6; 100,6]);
      annotation (defaultComponentName="cooler", Diagram, Icon(
          Rectangle(extent=[-90,20; 90,-20], style(
                color=0,
                rgbcolor={0,0,0},
                fillColor=7,
                rgbfillColor={255,255,255})),
          Ellipse(extent=[-60,16; -28,-16], style(
                color=0,
                rgbcolor={0,0,0},
                fillColor=7,
                rgbfillColor={255,255,255})),
          Ellipse(extent=[28,16; 60,-16], style(
                color=0,
                rgbcolor={0,0,0},
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
                fillColor=7,
                rgbfillColor={255,255,255},
                fillPattern=1)),
          Line(points=[-44,16; 44,-16], style(
                color=0,
                rgbcolor={0,0,0},
                fillColor=0,
                rgbfillColor={0,0,0},
                fillPattern=1)),
          Line(points=[-44,-16; 44,16], style(
                color=0,
                rgbcolor={0,0,0},
                fillColor=0,
                rgbfillColor={0,0,0},
                fillPattern=1)),
          Text(
            extent=[-102,82; 100,20],
            string="%name",
              style(
                color=0,
                rgbcolor={0,0,0},
                fillColor=43,
                rgbfillColor={255,85,85},
                fillPattern=1))),
        Documentation(info="<html>
<p>This is the most generic interface for an air cooler, whose heat-exchange model 
is left completely unspecified.</p><p>The cooler provides two flow ports and two temperature
outputs, which will output the temperature of the associate flow; the temperature
input is going to be used by some (yet unspecified) internal system to set the outlet
temperature.</p>
</html>"));
      IO.TemperatureOutput T_process_in "Process inlet temperature" 
                                                                 annotation (extent=[-88,-20;
            -100,-8]);
      IO.TemperatureOutput T_process_out "Process outlet temperature" 
                                                                   annotation (extent=[88,-20;
            100,-8]);
      IO.TemperatureInput T_ref "Reference temperature for the process outlet" 
        annotation (extent=[10,-40; -10,-20],rotation=90);
    end Abstract;
      
    model Simple "A simple cooler implementation" 
      extends Flow.UnitOperations.Coolers.Abstract;
        
      import Flow.Measurements.FlowTemperature;
        
      protected 
      Flow.Sink sink "Makes up for lost heat" 
                                             annotation (extent=[0,80; 20,100]);
      annotation (Diagram, Icon,
        Documentation(info="<html>
<p>This is a straightforward implementation of the abstract cooler that sets the
cooler's outlet temperature according to a given reference. There is a build-in,
customisable lag after which the outlet temperature will reach the reference value.</p>
<p>The outlet temperature can reach only values between inlet and environment
temperature; if the reference is outside these limits, the boolean output 
<tt>isSaturated</tt> will allow an external control algorithm to notice this and avoid
a wind-up situation.</p>
</html>"));
      FlowTemperature T_out "Outlet temperature measurement" 
        annotation (extent=[40,60; 60,80]);
      FlowTemperature T_in "Inlet temperature measurement" 
        annotation (extent=[-70,60; -50,80]);
      Modelica.Blocks.Continuous.FirstOrder lag(       initType=Modelica.Blocks.
            Types.Init.SteadyState, T=120) 
          "A lag representing the inner control algorithm setting the outlet temperature"
        annotation (extent=[20,0; 40,20],   rotation=0);
      Modelica.Blocks.Nonlinear.VariableLimiter limiter 
          "Keeps requested temperature within reason" 
        annotation (extent=[-20,0; 0,20]);
      constant Real eps = 0.01;
      outer Modelica.SIunits.Temperature T_env;
        
      public 
      Modelica.Blocks.Interfaces.BooleanOutput isSaturated 
          "Whether the current set point is saturated" 
        annotation (extent=[20,-40; 40,-20], rotation=270);
      Modelica.SIunits.HeatFlowRate Q = sink.inlet.H "Cooling duty";
        
    equation 
      limiter.limit2 = T_env;
        
      // No material loss
      sink.inlet.n = 0*sink.inlet.n;
        
      // Set the two variables to be equal
      // NOTE must do it here, these are both output variables.
      T_process_out = lag.y;
      connect(T_in.inlet, inlet) annotation (points=[-68,70; -80,70; -80,
              -2.22045e-16; -94,-2.22045e-16],
                                             style(color=62, rgbcolor={0,127,127}));
      connect(T_in.T, T_process_in) annotation (points=[-60,62; -60,-14; -94,-14],
          style(color=3, rgbcolor={0,0,255}));
      connect(T_out.T, T_process_out) annotation (points=[50,62; 50,-14; 94,-14],
          style(color=3, rgbcolor={0,0,255}));
      connect(T_out.outlet, outlet) annotation (points=[58,70; 76,70; 76,
              -2.22045e-16; 94,-2.22045e-16],
                                            style(color=62, rgbcolor={0,127,127}));
      connect(sink.inlet, T_in.outlet) 
        annotation (points=[1,90; -20,90; -20,70; -52,70], style(color=62,
            rgbcolor={0,127,127}));
      connect(T_out.inlet, T_in.outlet) 
        annotation (points=[42,70; -52,70], style(color=62, rgbcolor={0,127,127}));
      connect(limiter.limit1, T_in.T) annotation (points=[-22,18; -40,18; -40,40;
            -60,40; -60,62],
                         style(color=74, rgbcolor={0,0,127}));
      connect(limiter.u, T_ref) annotation (points=[-22,10; -56,10; -56,-16; 0,
              -16; 0,-30; -5.55112e-16,-30],
                                           style(color=74, rgbcolor={0,0,127}));
      connect(lag.u, limiter.y) 
        annotation (points=[18,10; 1,10], style(color=74, rgbcolor={0,0,127}));
    algorithm 
        
      when T_ref < T_env-eps or T_ref > T_process_in+eps then
        isSaturated := true;
      elsewhen T_ref > T_env+eps and T_ref < T_process_in-eps then
        isSaturated := false;
      end when;
        
    end Simple;
      
    partial model withAbstractExchanger 
        "An abstract exchanger-based cooler for a process stream" 
      extends Flow.UnitOperations.Coolers.Abstract;
        
      import Units.Temperature;
        
      Modelica.SIunits.VolumeFlowRate V_air = mfc.V "Coolant flow rate";
      Temperature T_coolant_in =  exchanger.T_cold_2 
          "Coolant inlet temperature";
      Temperature T_coolant_out = exchanger.T_cold_1 
          "Coolant outlet temperature";
        
        annotation (defaultComponentName="cooler", Documentation(info="<html>
<p>This is the interface for an air cooler based on a detailed but not specified heat-exchanger
model. It includes an internal PI controller that sets the process' outlet temperature by 
manipulating the coolant flow. The coolant itself enters at environment temperature.</p>
</html>"));
        
      replaceable Flow.UnitOperations.HeatExchangers.Abstract exchanger 
          "The heat exchanger implementing the cooler" 
        annotation (extent=[-64,28; 16,108]);
      Measurements.GasFlowController mfc "Mass flow controller for cooling air"
        annotation (extent=[60,46; 40,66], rotation=270);
      protected 
      Flow.Sources.Environment env "Environmental air source" 
        annotation (extent=[100,46; 80,66], rotation=0);
      Flow.Sink airSink "Air outlet sink" 
                                         annotation (extent=[50,70; 70,90]);
      public 
      Control.CoolerControl K annotation (extent=[12,4; 28,20]);
    equation 
      connect(exchanger.hot_2, outlet) annotation (points=[-60,56; -60,
              -2.22045e-16; 94,-2.22045e-16],
                                            style(color=62, rgbcolor={0,127,127}));
      connect(inlet, exchanger.hot_1) annotation (points=[-94,-2.22045e-16; -80,
              -2.22045e-16; -80,80; -60,80],
                                           style(color=62, rgbcolor={0,127,127}));
      connect(airSink.inlet, exchanger.cold_1) annotation (points=[51,80; 12,80],
          style(color=62, rgbcolor={0,127,127}));
      connect(env.outlet, mfc.inlet) annotation (points=[81,56; 50,56],
          style(color=62, rgbcolor={0,127,127}));
      connect(mfc.outlet, exchanger.cold_2) annotation (points=[40,56; 12,56],
          style(color=62, rgbcolor={0,127,127}));
      connect(exchanger.T_hot_1, T_process_in) annotation (points=[-60,88; -74,88;
            -74,-14; -94,-14],
                             style(color=3, rgbcolor={0,0,255}));
      connect(exchanger.T_hot_2, T_process_out) annotation (points=[-60,48; -70,
            48; -70,-14; 94,-14], style(color=3, rgbcolor={0,0,255}));
      connect(K.V, mfc.V) annotation (points=[29.6,12; 50,12; 50,46], style(color=
             3, rgbcolor={0,0,255}));
      connect(K.T_r, T_ref) annotation (points=[10.4,12; 0,12; 0,-30;
              -5.55112e-16,-30],
                               style(
          color=3,
          rgbcolor={0,0,255},
          pattern=3));
      connect(exchanger.T_hot_2, K.T_m) annotation (points=[-60,48; -70,48; -70,
            -14; 20,-14; 20,2.4], style(color=3, rgbcolor={0,0,255}));
    end withAbstractExchanger;
      
    model LMTD "A cooler implemented with a LMTD heat exchanger" 
      extends Flow.UnitOperations.Coolers.withAbstractExchanger(redeclare 
            HeatExchangers.LMTD exchanger);
      annotation (Diagram, Documentation(info="<html>
<p>A cooler using a LMTD implementation. It can for instance be applied as the
anode-loop cooler or as a recuperating heat exchanger on the anode side.</p>
</html>"),
        Icon);
    end LMTD;
      
    model Discretised "A cooler implemented with a discretised heat exchanger" 
      extends Flow.UnitOperations.Coolers.withAbstractExchanger(redeclare 
            HeatExchangers.Discretised exchanger);
      annotation (Documentation(info="<html>
<p>A cooler using a discretised implementation. It can for instance be applied as the
cathode-loop cooler (condenser).</p>
</html>"));
    end Discretised;
      
      package Test 
        partial model AbstractCoolerTest "Generic test for air coolers" 
          
          import Modelica.SIunits.VolumeFlowRate;
          import Units.Temperature;
          
          inner parameter Modelica.SIunits.Pressure p_env = 101325;
          inner parameter Modelica.SIunits.Temperature T_env = 298.15;
          inner parameter Units.RelativeHumidity RH_env = 60;
          
        protected 
          Sink sink         annotation (extent=[60,-4; 68,4]);
          annotation (Diagram,
            experiment(StopTime=80),
            experimentSetupOutput);
          Sources.Solution sol(   T=330) 
                                      annotation (extent=[-100,-20; -80,0]);
        public 
          replaceable Abstract cooler       annotation (extent=[-20,-20; 20,20]);
          Measurements.LiquidPump pump 
                    annotation (extent=[-60,-20; -40,0]);
          parameter VolumeFlowRate solution = 10E-6/60; // 10 ml/min
          parameter Temperature target = 315;
          
        equation 
          pump.V = solution;
          cooler.T_ref = target;
          
          connect(cooler.outlet, sink.inlet) annotation (points=[18.8,
                1.06581e-15; 39.4,1.06581e-15; 39.4,3.88578e-17; 60.4,
                3.88578e-17],                                          style(color=
                  62, rgbcolor={0,127,127}));
          connect(pump.outlet, cooler.inlet) annotation (points=[-50,
                5.55112e-16; -30,5.55112e-16; -30,1.06581e-15; -18.8,
                1.06581e-15],                                         style(color=
                  62, rgbcolor={0,127,127}));
          connect(sol.outlet, pump.inlet) annotation (points=[-90,-10; -50,-10],
              style(color=62, rgbcolor={0,127,127}));
        end AbstractCoolerTest;
        
        model SimpleCoolerTest 
          extends AbstractCoolerTest(redeclare Simple cooler);
          annotation (experiment(StopTime=3600));
          
        end SimpleCoolerTest;
        
        model LMTDCoolerTest "Test for the LMTD-based air cooler" 
          extends AbstractCoolerTest(redeclare LMTD cooler);
          
          annotation (experiment(StopTime=5000), experimentSetupOutput);
        end LMTDCoolerTest;
        
        model DiscretisedCoolerTest "Test for the discretised air cooler" 
          extends AbstractCoolerTest(redeclare Discretised cooler);
        end DiscretisedCoolerTest;
      end Test;
    end Coolers;
    
    package Stack "Models of fuel cells" 
      
    partial model Abstract "A generic DMFC stack" 
        
      import Modelica.SIunits.Area;
      import Modelica.SIunits.Current;
      import Modelica.SIunits.CurrentDensity;
      import Modelica.SIunits.Concentration;
      import Modelica.SIunits.DiffusionCoefficient;
      import Modelica.SIunits.Efficiency;
      import Modelica.SIunits.HeatCapacity;
      import Modelica.SIunits.Length;
      import Modelica.SIunits.MoleFraction;
      import Modelica.SIunits.PartialPressure;
      import Modelica.SIunits.Pressure;
      import Modelica.SIunits.StoichiometricNumber;
      import Modelica.SIunits.Temperature;
      import Modelica.SIunits.Voltage;
      import Modelica.Constants.eps;
      import Modelica.Constants.R;
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
      import Units.F;
        
      annotation (defaultComponentName="cell", Icon(Rectangle(extent=[-100,60; 100,0], style(
                color=0,
                rgbcolor={0,0,0},
                fillColor=47,
                rgbfillColor={255,170,85})),
                                           Rectangle(extent=[-100,2; 100,-60], style(
                color=0,
                rgbcolor={0,0,0},
                fillColor=71,
                rgbfillColor={85,170,255}))),
                                            Diagram,
        DymolaStoredErrors,
        Documentation(info="<html>
<p>This class implements a DMFC stack from the point of view of reactant flows. A
modelling of the voltage is <em>not</em> included, and must be implemented by child classes.
A temperature port and two ports for methanol concentration on the anode side (inlet and
outlet) are featured.</p>
 
<h3>Modelled Phenomena</h3>
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
coefficient in water is required. It is assumed to vary exponentially with temperature.</p>
 
<p>Parameters have been taken from Krewer et al., unless differently stated.</p>

<h3>Implementation details</h3>
<p>The two inlets and the two outlets are connected to the \"nexus\", an internal protected
(i.e. invisible to the user) object, that accounts for components lost in reactions and
energy that leaves as electric power (=I*V).</p>
 
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
      Flow.FlowPort cathode_inlet "The cathode flow's inlet" 
        annotation (extent=[-110,20; -90,40]);
      Flow.FlowPort cathode_outlet "The cathode flow's outlet" 
        annotation (extent=[90,20; 110,40]);
      Flow.FlowPort anode_inlet "The anode flow's inlet" 
        annotation (extent=[-110,-40; -90,-20]);
      Flow.FlowPort anode_outlet "The anode flow's outlet" 
        annotation (extent=[90,-40; 110,-20]);
      PositivePin plus "Pole connected to the cathode" 
                                        annotation (extent=[-70,50; -50,70]);
      NegativePin minus "Pole connected to the anode" 
                                      annotation (extent=[50,50; 70,70]);
      protected 
      Flow.Measurements.FlowTemperature cathodeT 
          "Cathode temperature measurement" 
        annotation (extent=[60,20; 80,40]);
      Flow.Measurements.FlowConcentration anodeTC 
          "Anode temperature and concentration measurement" 
                                      annotation (extent=[60,-40; 80,-20]);
      Flow.Sink nexus "Connection of all flows" 
                        annotation (extent=[-32,-10; -12,10]);
        
      public 
      outer Pressure p_env "Environment pressure";
      outer Temperature T_env "Enviroment temperature";
        
      parameter Integer cells = 1 "Number of cells";
      parameter Length d = 142E-6 "Membrane thickness";
      parameter Area A = 26E-4 "Membrane active area";
      parameter HeatCapacity Cp = 24.7 "Overall heat capacity of the stack";
      parameter DiffusionCoefficient D = 6E-10 
          "Methanol diffusion coefficient in the membrane";
      parameter Boolean enableSanityChecks = true 
          "Whether to activate checks for some non-negative quantities";
        
      // Parameters for N115 membrane.  
      Real k_d = 4.2 + (T-303.15)/40 "Drag factor for N115";
      MassTransportCoefficient k_m = 15.6E-6*exp(2436*(1/T-1/333)) 
          "Mass transport coefficient";
        
      Real a = k_m/(1+k_m*d/D) 
          "Partial derivative of crossover flux wrt. concentration";
      Real aAn = a*A*cells 
          "Partial derivative of crossover flow wrt. concentration";
      Real b = 1/(1+k_m*d/D) 
          "Opposite of partial derivative of crossover flux wrt. anodic reaction rate";
        
      Voltage V_rev "Reversible voltage";
      Efficiency eta_thermo = V/V_rev "Electrochemical efficiency";
      Efficiency eta_use = i / (i+6*F*N_x) 
          "Fraction of lost methanol reacting on anode";
      Efficiency eta_total = eta_thermo*eta_use "Overall cell efficiency";
        
      Current I = -plus.i "Cell current (generator convention)";
      Voltage V = plus.v - minus.v "Cell voltage";
      CurrentDensity i = I/A "Cell current density";
        
      MolarFlow n_H = cells*I/F "Proton flow through the membrane";
      MolarFlow n_x "Crossover methanol flow";
      MolarFlux N_H = n_H / A / cells "Proton flux";
      MolarFlux N_x = n_x / A / cells "Crossover methanol flux";
      MolarFlow n_drag_h2o = n_H * k_d "Drag water flow";
      MolarFlux N_drag_h2o = n_drag_h2o / A / cells "Drag water flux";
        
      // KEEP THE INITIAL VALUE, or initialisation will crash on assertion.
      Concentration c(start=1000) = anodeTC.c 
          "Methanol concentration, outlet is representative";
      Concentration c_cl(start=100) "Catalyst-layer methanol concentration";
        
      PartialPressure p_o2 "Oxygen partial pressure, outlet is representative";
      PartialPressure p_h2o 
          "Cathodic water partial pressure, outlet is representative";
      Flow.IO.TemperatureOutput T "Representative stack temperature" 
                                                             annotation (extent=[100,-8; 120,12]);
        
      /* This group of vectors represents the coefficients by which
   * proton (*_nu) and crossover-methanol (*_xi) flows must be 
   * multiplied to  find the associated production terms for all
   * species on cathode and anode; consumption terms are obviously
   * negative. */
      protected 
      StoichiometricNumber[:] cathode_nu = {0, 1/2+k_d, -1/4, 0, 0};
      StoichiometricNumber[:] anode_nu = {-1/6, -1/6-k_d, 0, 1/6, 0};
      constant StoichiometricNumber[:] cathode_xi = {0, 2, -3/2, 1, 0};
      constant StoichiometricNumber[:] anode_xi = {-1, 0, 0, 0, 0};
        
    equation 
      // Anode-side mass balance, accounting for reaction, drag and crossover
      anode_inlet.n + anode_outlet.n + anode_nu*n_H + anode_xi*n_x = zeros(size(All,1));
        
      // Cathode-side mass balance, accounting for reaction, drag and crossover
      cathode_inlet.n + cathode_outlet.n + cathode_nu*n_H + cathode_xi*n_x = zeros(size(All,1));
        
      // The energy "lost" from the heat balance is the electrical power.
      nexus.inlet.H = I*V + der(T)*Cp;
        
      // Definition of oxygen partial pressure. On the denominator, the sum of vapours (methanol and water) and gases (all others).
      p_o2 = p_env * cathodeT.inlet.n[Oxygen] / (sum(cathodeT.vapour) + sum(cathodeT.inlet.n[Incondensable]));
        
      // Definition of water partial pressure (on the cathode side). On the denominator, the sum of vapours (methanol and water) and gases (all others).
      p_h2o = p_env * cathodeT.vapour[Water]  / (sum(cathodeT.vapour) + sum(cathodeT.inlet.n[Incondensable]));
        
      // Methanol transport: binds c, c_cl and i (N_x is a function of c_ac, N_H of i).
      k_m * (c-c_cl) = N_x + N_H/6;
        
      // Crossover methanol flux.
      N_x = D/d * c_cl;
        
      /* The reversible voltage; the terms are:
   * - Standard reaction enthalpy, minus
   * - Temperature times Standard reaction entropy, minus
   * - Correction factor for oxygen activity = 0.2, 2.4142 = log(1/0.2^1.5)
   * - All multiplied by the number of cells
   * Other activities are assumed unitary. */
      V_rev = -(-726770 - T_env*(-81.105) + R*T_env*2.4142)/6/F * cells;
        
      // Setting the temperatures of cathode and anode to be equal
      /* NOTE: a connect() would prettier, but both these variables are outputs
   * and can therefore not be connected to each other.
   * It is neither a good idea to connect them both to FuelCell.T, since
   * this would give two sources to one output.
   * Therefore, we use an equation. */
      anodeTC.T = cathodeT.T;
        
      // Charge balance
      plus.i + minus.i = 0;
        
      if enableSanityChecks then
        // Sanity check: crash simulation if conditions are unphysical
        assert( c_cl >= 0, "==> Methanol catalyst-layer concentration is negative ("+String(c_cl)+" mol/m^3) at temperature "+String(T)+" K, bulk concentration "+String(c)+" mol/m^3.");
          
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
      connect(anodeTC.outlet, anode_outlet)       annotation (points=[78,-30; 100,
            -30], style(color=62, rgbcolor={0,127,127}));
      connect(anodeTC.inlet, nexus.inlet)          annotation (points=[62,-30;
              -40,-30; -40,4.44089e-16; -31,4.44089e-16],
                                                    style(color=62, rgbcolor={0,127,
              127}));
      connect(T, cathodeT.T) annotation (points=[110,2; 70,2; 70,22],
          style(color=3, rgbcolor={0,0,255}));
      connect(anode_inlet, nexus.inlet) annotation (points=[-100,-30; -46,-30;
              -46,4.44089e-16; -31,4.44089e-16],
                                               style(color=62, rgbcolor={0,127,
              127}));
    end Abstract;
      
    model ConstantVoltage "A simplified DMFC stack with constant voltage" 
      extends Abstract;
      import Modelica.SIunits.Voltage;
        
      annotation (Documentation(info="<html>
<p>This trivial class inherits from the <tt>Stack.Abstract</tt> class and allows to set a 
constant voltage for the stack.</p>
</html>"));
        
      parameter Voltage V0 = 0.5 "Cell voltage";
        
    equation 
      V = V0;
        
    end ConstantVoltage;
      
    model Thevenin "A DMFC stack with Thevenin-like voltage" 
      extends Abstract;
      import Modelica.SIunits.Voltage;
      import Modelica.SIunits.Resistance;
        
      parameter Voltage V0 = 0.7 "Open-circuit voltage";
      parameter Resistance R = 0.005 "Internal resistance";
        
      annotation (Documentation(info="<html>
<p>This class implements a voltage model that emulates a Thevenin equivalent circuit. It is possible
to set the open-circuit voltage and the specific resistance of the stack.</p>
<p>Note that the open-circuit value to set is not the one measured on the actual cell, but the one 
that would result by extrapolating the characteristic of the ohmic region to the value of no 
current.</p>
</html>"));
    equation 
      V = V0 - R*I;
    end Thevenin;
      
      package Test 
        partial model AbstractStackTest 
          "Generic test suite for fuel-cell models" 
          
          import Modelica.Electrical.Analog.Sources.ConstantCurrent;
          
          inner parameter Modelica.SIunits.Pressure p_env = 101325 
            "Environment pressure";
          inner parameter Modelica.SIunits.Temperature T_env = 298.15 
            "Enviroment temperature";
          inner parameter Units.RelativeHumidity RH_env = 60 
            "Environment relative humidity";
          
          parameter Modelica.SIunits.VolumeFlowRate anodeFlow = 30E-6/60 
            "Anodic volumetric flow rate";
          parameter Modelica.SIunits.VolumeFlowRate cathodeFlow = 30E-5/60 
            "Cathodic volumetric flow rate";
          parameter Units.Temperature anodeInletTemperature = 330 
            "Anodic inlet temperature";
          
          replaceable Abstract stack 
                            annotation (extent=[6,0; 42,34]);
          Sources.Solution methanolSolution(   T=anodeInletTemperature) 
                                            annotation (extent=[-66,-30; -54,-18]);
          annotation (Diagram);
          Measurements.LiquidPump pump "Pump for the anode flow" 
                                              annotation (extent=[-42,-30; -30,-18]);
          Sources.Environment air 
                              annotation (extent=[-70,28; -50,48]);
          Measurements.GasFlowController blower 
                                   annotation (extent=[-40,34; -32,42]);
          Sink anodeSink     annotation (extent=[62,10; 68,16]);
          Sink cathodeSink     annotation (extent=[62,18; 68,24]);
          ConstantCurrent I_cell(I=5) annotation (extent=[12,48; 34,72]);
          Modelica.Electrical.Analog.Basic.Ground ground 
            "Negative pole to zero voltage" annotation (extent=[38,40; 58,60]);
        equation 
          pump.V = anodeFlow;
          blower.V = cathodeFlow;
          
          connect(methanolSolution.outlet, pump.inlet) 
                                                  annotation (points=[-60,-24; -36,
                -24],        style(color=62, rgbcolor={0,127,127}));
          connect(blower.outlet, stack.cathode_inlet)    annotation (points=[-36,42;
                -18,42; -18,22.1; 6,22.1], style(color=62, rgbcolor={0,127,127}));
          connect(air.outlet, blower.inlet) 
                                       annotation (points=[-51,38; -36,38],
                                style(color=62, rgbcolor={0,127,127}));
          connect(cathodeSink.inlet, stack.cathode_outlet)       annotation (points=[62.3,21;
                52.15,21; 52.15,22.1; 42,22.1], style(color=62, rgbcolor={0,127,127}));
          connect(anodeSink.inlet, stack.anode_outlet)       annotation (points=[62.3,13;
                52.15,13; 52.15,11.9; 42,11.9], style(color=62, rgbcolor={0,127,127}));
          connect(ground.p, I_cell.n) 
            annotation (points=[48,60; 34,60], style(color=3, rgbcolor={0,0,255}));
          connect(I_cell.p, stack.plus)    annotation (points=[12,60; 12,27.2;
                13.2,27.2],
              style(color=3, rgbcolor={0,0,255}));
          connect(I_cell.n, stack.minus)    annotation (points=[34,60; 34.8,60;
                34.8,27.2], style(color=3, rgbcolor={0,0,255}));
          connect(pump.outlet, stack.anode_inlet)    annotation (points=[-36,-18;
                -16,-18; -16,11.9; 6,11.9], style(color=62, rgbcolor={0,127,127}));
        end AbstractStackTest;
        
        model ConstantVoltageStackTest "Test for the constant-voltage model" 
          extends AbstractStackTest(redeclare ConstantVoltage stack,
              anodeFlow = stack.cells*30E-6/60,
              cathodeFlow = stack.cells*30E-5/60);
          annotation (Diagram);
        end ConstantVoltageStackTest;
        
        model TheveninStackTest "Test for the Thevenin-circuit model" 
          extends AbstractStackTest(redeclare Thevenin stack(cells=3,V0=2.1,R=0.015),
              anodeFlow = stack.cells*30E-6/60,
              cathodeFlow = stack.cells*30E-5/60);
        end TheveninStackTest;
      end Test;
    end Stack;
    
    package Test 
      model SeparatorTest "Test case for the separator unit" 
        
        Separator separator annotation (extent=[-22,-12; 26,38]);
      protected 
        Flow.Sink liquidSink 
                          annotation (extent=[46,-16; 54,-8]);
        annotation (Diagram);
        Sink gasSink       annotation (extent=[46,30; 54,38]);
        Sources.Environment env         annotation (extent=[-54,22; -34,42]);
        Sources.Solution solution         annotation (extent=[-78,8; -68,18]);
      public 
        inner parameter Modelica.SIunits.Pressure p_env = 101325;
        inner parameter Units.Temperature T_env = 298.15;
        inner parameter Units.RelativeHumidity RH_env = 60;
        
      equation 
        sum(env.outlet.n) = -1;
        sum(solution.outlet.n) = -2;
        
        connect(env.outlet, separator.inlet)        annotation (points=[-35,32;
              -22,32; -22,13],             style(color=62, rgbcolor={0,127,127}));
        connect(solution.outlet, separator.inlet) 
          annotation (points=[-73,13; -22,13], style(color=62, rgbcolor={0,127,
                127}));
        connect(separator.liquidOutlet, liquidSink.inlet) 
          annotation (points=[18.8,3; 18,3; 18,-12; 46.4,-12], style(color=62,
              rgbcolor={0,127,127}));
        connect(separator.gasOutlet, gasSink.inlet)    annotation (points=[18.8,23;
              18.4,23; 18.4,34; 46.4,34], style(color=62, rgbcolor={0,127,127}));
      end SeparatorTest;
      
      model BurnerTest 
        
        import Modelica.SIunits.VolumeFlowRate;
        
        inner parameter Modelica.SIunits.Pressure p_env = 101325 
          "Environment pressure";
        inner parameter Units.Temperature T_env = 298.15 
          "Enviroment temperature";
        inner parameter Units.RelativeHumidity RH_env = 60 
          "Environment relative humidity";
        
        parameter VolumeFlowRate solution = 10E-6/60; // 10 ml/min
        parameter VolumeFlowRate air =      1E-3/60; // 1  l/min
        
        Burner burner annotation (extent=[-2,-20; 38,20]);
        Sources.Solution methanolSolution annotation (extent=[-100,20; -80,40]);
        Measurements.LiquidPump pump 
                  annotation (extent=[-60,40; -40,20]);
        annotation (Diagram, Documentation(info="<html>
</html>"));
        Sources.Environment env 
                            annotation (extent=[-100,-20; -80,0]);
        Measurements.GasFlowController mfc 
                              annotation (extent=[-60,-20; -40,0]);
        Sink sinkPort     annotation (extent=[80,-10; 100,10]);
        Measurements.FlowTemperature T_in 
                             annotation (extent=[-30,-10; -10,10]);
      equation 
        pump.V = solution;
        mfc.V = air;
        connect(pump.inlet, methanolSolution.outlet) annotation (points=[-50,30;
              -90,30], style(
            color=62,
            rgbcolor={0,127,127},
            fillColor=45,
            rgbfillColor={255,128,0},
            fillPattern=10));
        connect(mfc.inlet, env.outlet) annotation (points=[-50,-10; -81,-10],
            style(
            color=62,
            rgbcolor={0,127,127},
            fillColor=45,
            rgbfillColor={255,128,0},
            fillPattern=10));
        connect(T_in.outlet, burner.inlet) annotation (points=[-12,6.10623e-16;
              -7,6.10623e-16; -7,1.22125e-15; -1.6,1.22125e-15], style(
            color=62,
            rgbcolor={0,127,127},
            fillColor=45,
            rgbfillColor={255,128,0},
            fillPattern=10));
        connect(T_in.inlet, pump.outlet) annotation (points=[-28,6.10623e-16;
              -40,6.10623e-16; -40,20; -50,20],
                                            style(
            color=62,
            rgbcolor={0,127,127},
            fillColor=45,
            rgbfillColor={255,128,0},
            fillPattern=10));
        connect(T_in.inlet, mfc.outlet) annotation (points=[-28,6.10623e-16;
              -40,6.10623e-16; -40,5.55112e-16; -50,5.55112e-16],
                                                              style(
            color=62,
            rgbcolor={0,127,127},
            fillColor=45,
            rgbfillColor={255,128,0},
            fillPattern=10));
        connect(burner.outlet, sinkPort.inlet) annotation (points=[38.4,
              1.22125e-15; 60.2,1.22125e-15; 60.2,4.44089e-16; 81,4.44089e-16],
            style(color=62, rgbcolor={0,127,127}));
      end BurnerTest;
      
      model MixerTest "Test for the mixer unit" 
        inner parameter Modelica.SIunits.Pressure p_env = 101325;
        inner parameter Units.Temperature T_env = 298.15;
        inner parameter Units.RelativeHumidity RH_env = 60;
        
        Mixer mixer(
          c(fixed=true),
          T(fixed=true),
          V(fixed=true, start=500E-6)) 
                    annotation (extent=[-20,0; 0,20]);
        Sources.Solution anodicLoop(   T=330) 
          "Solution coming from the anodic loop" 
          annotation (extent=[-20,40; 0,60]);
        Sources.Methanol fuelTank "Methanol from the fuel tank" 
          annotation (extent=[0,-40; 20,-20]);
        Sources.Solution condenser(   C=0, T=310) 
          "Water recovered from the cathode outlet" 
          annotation (extent=[20,0; 40,20]);
        Measurements.FlowTemperature flowTemperature 
                                        annotation (extent=[-64,-26; -44,-6]);
        Sink sinkPort     annotation (extent=[-28,-20; -20,-12]);
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
      
    end Test;
  end UnitOperations;
  
end Flow;

