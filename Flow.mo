within ;
          /**
 * (c) Federico Zenith, 2008-2010.
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

    flow Units.MolarFlow[Thermo.Species] n;
    flow Modelica.SIunits.EnthalpyFlowRate H;

    annotation (Documentation(info="<html>
<p>This is a connector for the various tank units; it ensures continuity of enthalpy and
molar flows. It consists of two flow variables, the <em>enthalpy flow</em> and the array of 
<em>molar flows</em>.</p>
<p>The enthalpy flow is defined as the enthalpy necessary to bring the components in the 
molar-flow array from standard conditions, which is defined as starting from fundamental 
elements in their native state, to the actual ones. Using this definition, it is fairly 
easy to model chemical reactions.</p>
</html>"), Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
              {100,100}}), graphics={Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,255},
            pattern=LinePattern.None,
            fillColor={0,127,127},
            fillPattern=FillPattern.Solid), Text(
            extent={{-100,160},{100,100}},
            lineColor={0,127,127},
            textString="%name")}));
  end FlowPort;

  connector PressurePort "Add pressure value to flows"

    Modelica.SIunits.Pressure p "Pressure";
    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}), graphics={Rectangle(
            extent={{-100,100},{100,-100}},
            fillColor={85,255,255},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None), Text(
            extent={{-100,160},{100,100}},
            lineColor={85,255,255},
            textString="%name")}),       Documentation(info="<html>
<p>Adds the value of pressure to the connection of molar and enthalpic flows.</p>
</html>"));
  end PressurePort;

  package IO "Input and output variables"

  connector ConcentrationInput = input Modelica.SIunits.Concentration 
   annotation (defaultComponentName="c",
    Icon(graphics={Polygon(
            points={{-100,100},{100,0},{-100,-100},{-100,100}},
            lineColor={0,0,127},
            fillColor={0,0,127},
            fillPattern=FillPattern.Solid)},
         coordinateSystem(extent={{-100,-100},{100,100}}, preserveAspectRatio=true, initialScale=0.2)),
    Diagram(coordinateSystem(
          preserveAspectRatio=true, initialScale=0.2,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={Polygon(
            points={{0,50},{100,0},{0,-50},{0,50}},
            lineColor={0,0,127},
            fillColor={0,0,127},
            fillPattern=FillPattern.Solid), Text(
            extent={{-10,85},{-10,60}},
            lineColor={0,0,127},
            textString="%name")}));
  connector CurrentInput = input Modelica.SIunits.Current 
   annotation (defaultComponentName="I",
    Icon(graphics={Polygon(
            points={{-100,100},{100,0},{-100,-100},{-100,100}},
            lineColor={0,0,127},
            fillColor={0,0,127},
            fillPattern=FillPattern.Solid)},
         coordinateSystem(extent={{-100,-100},{100,100}}, preserveAspectRatio=true, initialScale=0.2)),
    Diagram(coordinateSystem(
          preserveAspectRatio=true, initialScale=0.2,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={Polygon(
            points={{0,50},{100,0},{0,-50},{0,50}},
            lineColor={0,0,127},
            fillColor={0,0,127},
            fillPattern=FillPattern.Solid), Text(
            extent={{-10,85},{-10,60}},
            lineColor={0,0,127},
            textString="%name")}));
  connector PressureInput = input Modelica.SIunits.Pressure 
  annotation (defaultComponentName="I",
    Icon(graphics={Polygon(
            points={{-100,100},{100,0},{-100,-100},{-100,100}},
            lineColor={0,0,127},
            fillColor={0,0,127},
            fillPattern=FillPattern.Solid)},
         coordinateSystem(extent={{-100,-100},{100,100}}, preserveAspectRatio=true, initialScale=0.2)),
    Diagram(coordinateSystem(
          preserveAspectRatio=true, initialScale=0.2,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={Polygon(
            points={{0,50},{100,0},{0,-50},{0,50}},
            lineColor={0,0,127},
            fillColor={0,0,127},
            fillPattern=FillPattern.Solid), Text(
            extent={{-10,85},{-10,60}},
            lineColor={0,0,127},
            textString="%name")}));
  connector TemperatureInput = input Units.Temperature 
  annotation (defaultComponentName="T",
    Icon(graphics={Polygon(
            points={{-100,100},{100,0},{-100,-100},{-100,100}},
            lineColor={0,0,127},
            fillColor={0,0,127},
            fillPattern=FillPattern.Solid)},
         coordinateSystem(extent={{-100,-100},{100,100}}, preserveAspectRatio=true, initialScale=0.2)),
    Diagram(coordinateSystem(
          preserveAspectRatio=true, initialScale=0.2,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={Polygon(
            points={{0,50},{100,0},{0,-50},{0,50}},
            lineColor={0,0,127},
            fillColor={0,0,127},
            fillPattern=FillPattern.Solid), Text(
            extent={{-10,85},{-10,60}},
            lineColor={0,0,127},
            textString="%name")}));
  connector VolumeFlowRateInput = input Modelica.SIunits.VolumeFlowRate 
  annotation (defaultComponentName="V",
    Icon(graphics={Polygon(
            points={{-100,100},{100,0},{-100,-100},{-100,100}},
            lineColor={0,0,127},
            fillColor={0,0,127},
            fillPattern=FillPattern.Solid)},
         coordinateSystem(extent={{-100,-100},{100,100}}, preserveAspectRatio=true, initialScale=0.2)),
    Diagram(coordinateSystem(
          preserveAspectRatio=true, initialScale=0.2,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={Polygon(
            points={{0,50},{100,0},{0,-50},{0,50}},
            lineColor={0,0,127},
            fillColor={0,0,127},
            fillPattern=FillPattern.Solid), Text(
            extent={{-10,85},{-10,60}},
            lineColor={0,0,127},
            textString="%name")}));
  connector ConcentrationOutput = output Modelica.SIunits.Concentration 
  annotation (defaultComponentName="c",
    Icon(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={Polygon(
            points={{-100,100},{100,0},{-100,-100},{-100,100}},
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}),
    Diagram(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={Polygon(
            points={{-100,50},{0,0},{-100,-50},{-100,50}},
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid), Text(
            extent={{30,110},{30,60}},
            lineColor={0,0,127},
            textString="%name")}));
  connector PressureOutput = output Modelica.SIunits.Pressure 
  annotation (defaultComponentName="p",
    Icon(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={Polygon(
            points={{-100,100},{100,0},{-100,-100},{-100,100}},
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}),
    Diagram(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={Polygon(
            points={{-100,50},{0,0},{-100,-50},{-100,50}},
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid), Text(
            extent={{30,110},{30,60}},
            lineColor={0,0,127},
            textString="%name")}));
  connector TemperatureOutput = output Units.Temperature 
  annotation (defaultComponentName="T",
    Icon(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={Polygon(
            points={{-100,100},{100,0},{-100,-100},{-100,100}},
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}),
    Diagram(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={Polygon(
            points={{-100,50},{0,0},{-100,-50},{-100,50}},
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid), Text(
            extent={{30,110},{30,60}},
            lineColor={0,0,127},
            textString="%name")}));

  connector VolumeFlowRateOutput = output Modelica.SIunits.VolumeFlowRate 
  annotation (defaultComponentName="V",
    Icon(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={Polygon(
            points={{-100,100},{100,0},{-100,-100},{-100,100}},
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}),
    Diagram(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={Polygon(
            points={{-100,50},{0,0},{-100,-50},{-100,50}},
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid), Text(
            extent={{30,110},{30,60}},
            lineColor={0,0,127},
            textString="%name")}));
  connector VolumeOutput = output Modelica.SIunits.Volume 
  annotation (defaultComponentName="V",
    Icon(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={Polygon(
            points={{-100,100},{100,0},{-100,-100},{-100,100}},
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}),
    Diagram(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={Polygon(
            points={{-100,50},{0,0},{-100,-50},{-100,50}},
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid), Text(
            extent={{30,110},{30,60}},
            lineColor={0,0,127},
            textString="%name")}));
  end IO;

  model Sink "A general-purpose flow sink"

    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{
              -100,-100},{100,100}}),
                        graphics),
                         Icon(coordinateSystem(preserveAspectRatio=false,
            extent={{-100,-100},{100,100}}), graphics={
          Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Line(points={{-40,60},{60,60},{60,-60},{-40,-60}}, color={0,0,0}),
          Rectangle(
            extent={{-60,0},{60,0}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}),
      Documentation(info="<html>
<p>This very simple model is a terminal for flows leaving the system and
about which we do not care much.</p>
</html>"));
    FlowPort inlet "inlet for flow to discard" 
                      annotation (Placement(transformation(extent={{-120,-30},{
              -60,30}}, rotation=0)));
  end Sink;

  annotation (uses(Modelica(version="3.1"), Units(version="1"),
      Control(version="1"),
      Thermo(version="1")),                    Documentation(info="<html>
<p>This package contains various models related to fluid flow in stirred tanks.</p>
</html>"),
    version="1",
    conversion(noneFromVersion=""));

  package Sources "Sources of flows of various compositions"

  model Solution "Source of methanol solution with given concentration"
      import Thermo.mw;
      import Thermo.rho;
      import Thermo.h;
      import Thermo.Species;
      import Thermo.Phases;
      import Thermo.Incondensables;
      import Modelica.SIunits.Temperature;
      import Modelica.SIunits.MoleFraction;
      import Modelica.SIunits.Concentration;

    FlowPort outlet "Methanol solution" 
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, rotation=
               0)));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={Ellipse(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}),
                                           Documentation(info="<html>
<p>This item is a source for methanol-water solutions. Parameter <tt>C</tt>
allows to set the concentration in moler per <em>cubic metre</em>; note that
this is 1000 times the normal scale (1M = 1000 mol/m).</p>
</html>"),
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                100,100}}),
              graphics));

    parameter Concentration C = 1000 "Methanol concentration";
    parameter Temperature T = 298.15 "Temperature";

    MoleFraction x_ch3oh "Methanol molar fraction";
    MoleFraction x_h2o(start=0.9) "Water molar fraction";

  equation
    assert(C >= 0, "==> Negative concentration given in MethanolSolution object.");
    assert(C <= rho(T,Species.Methanol,Phases.Liquid)/mw(Species.Methanol),
           "==> Methanol concentration over limit (" + String(mw(Species.Methanol)/rho(T,Species.Methanol,Phases.Liquid)) + " mol/m3).");

    C = x_ch3oh / (x_ch3oh*mw(Species.Methanol)/rho(T,Species.Methanol,Phases.Liquid)
      + x_h2o*mw(Species.Water)/rho(T,Species.Water,Phases.Liquid));
    x_ch3oh + x_h2o = 1.0;

    outlet.n[Incondensables] = zeros(size(Incondensables, 1));
    outlet.n[Species.Methanol] / x_ch3oh = outlet.n[Species.Water] / x_h2o;
    outlet.H = outlet.n[Species.Methanol]*h(T,Species.Methanol,Phases.Liquid)
             + outlet.n[Species.Water]*h(T,Species.Water,Phases.Liquid);

  end Solution;

  model Methanol "Source of pure methanol"
      import Thermo.h;
      import Thermo.Species;
      import Thermo.Phases;
      import Thermo.mw;
      import Thermo.Incondensables;
      import Thermo.rho;
      import Units.Temperature;
      import Modelica.SIunits.AmountOfSubstance;
      import Modelica.SIunits.Mass;
      import Modelica.SIunits.Volume;

    FlowPort outlet "Pure methanol" 
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, rotation=
               0)));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={Ellipse(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,255},
              pattern=LinePattern.None,
              lineThickness=1,
              fillColor={255,0,0},
              fillPattern=FillPattern.Solid)}),
                                           Documentation(info="<html>
<p>This item is a source for a pure methanol stream. This model also keeps track 
of how much methanol has been released in terms of moles, mass and volume.</p>
</html>"),
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                100,100}}),
              graphics));

    outer Temperature T_env "Enviroment temperature";

    AmountOfSubstance n(start=0,fixed=true)
        "Released methanol moles from simulation start";
    Mass m = n*mw(Species.Methanol) "Released methanol mass";
    Volume V = m / rho(T_env, Species.Methanol, Phases.Liquid)
        "Released methanol volume";

  equation
    der(n) = -outlet.n[Species.Methanol];

    outlet.n[Incondensables] = zeros(size(Incondensables, 1));
    outlet.n[Species.Water] = 0;
    outlet.H = outlet.n[Species.Methanol]*h(T_env,Species.Methanol,Phases.Liquid);

  end Methanol;

  model Environment "A flow connection to environment conditions."

    Units.MolarFlow F "Total exchanged molar flow";
    Thermo.Air air "Data about environment air";

    FlowPort outlet "Environment-air port" 
                 annotation (Placement(transformation(extent={{80,-10},{100,10}},
              rotation=0)));
    annotation (defaultComponentName="env", Documentation(info="<html>
<p>This object generates a gas flow corresponding to ambient air, accounting
also for humidity.</p>
</html>"), Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics={
            Polygon(
              points={{20,-40},{46,42},{74,-40},{20,-40}},
              lineColor={0,0,0},
              pattern=LinePattern.None,
              fillPattern=FillPattern.VerticalCylinder,
              fillColor={0,127,0}),
            Polygon(
              points={{-80,-20},{-10,-20},{-44,82},{-80,-20}},
              lineColor={0,0,0},
              pattern=LinePattern.None,
              fillPattern=FillPattern.VerticalCylinder,
              fillColor={0,127,0}),
            Rectangle(
              extent={{-54,-20},{-36,-74}},
              lineColor={0,0,0},
              pattern=LinePattern.None,
              fillPattern=FillPattern.VerticalCylinder,
              fillColor={127,127,0}),
            Rectangle(
              extent={{40,-40},{52,-72}},
              lineColor={0,0,0},
              pattern=LinePattern.None,
              fillPattern=FillPattern.VerticalCylinder,
              fillColor={127,127,0}),
            Ellipse(
              extent={{2,94},{48,50}},
              lineColor={0,0,0},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Sphere,
              fillColor={255,255,85}),
            Text(
              extent={{-100,160},{100,100}},
              lineColor={0,0,0},
              fillColor={255,85,85},
              fillPattern=FillPattern.Solid,
              textString="%name")}),
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                100,100}}),
              graphics));

  equation
    outlet.n = F * air.y;
    outlet.H = F * air.H;

  end Environment;
  end Sources;

  package Measurements "Measurements on a flow"

  partial model FlowController "A unit modelling a pump or a MFC"

    import Thermo.mw;
    import Modelica.SIunits.MassFlowRate;
    import Units.MolarFlow;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={Ellipse(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid), Text(
              extent={{-100,160},{100,100}},
              lineColor={0,0,0},
              fillColor={255,85,85},
              fillPattern=FillPattern.Solid,
              textString="%name")}),
                              Documentation(info="<html>
<p>This class allows to set or read a certain overall molar or mass flow. It is not
immediately possible to set a <em>volume</em> flow, because this would entail
calculating the phase equilibrium, which we are not doing here since it is
computationally onerous; see the child classes <tt>Pump</tt> and 
<tt>GasFlowController</tt> if you need volume.</p>
</html>"),
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                100,100}}),
              graphics));
    FlowPort inlet "Unit inlet"   annotation (Placement(transformation(extent={
                {-10,-10},{10,10}}, rotation=0)));
    FlowPort outlet "Unit outlet"   annotation (Placement(transformation(extent=
               {{-10,90},{10,110}}, rotation=0)));

    MolarFlow F "Molar flow rate";
    MassFlowRate m "Mass flow rate";

    IO.VolumeFlowRateInput V "Volumetric flow rate" 
      annotation (Placement(transformation(extent={{110,-10},{90,10}}, rotation=
               0)));
  equation
    m = sum({inlet.n[i] * mw(i) for i in Thermo.Species});
    F = sum(inlet.n);

    connect(outlet, inlet) annotation (Line(points={{5.55112e-16,100},{
              5.55112e-16,75},{5.55112e-16,75},{5.55112e-16,50},{5.55112e-16,
              5.55112e-16},{5.55112e-16,5.55112e-16}}, color={0,127,127}));
  end FlowController;

  model GasFlowController "A flow controller with only gas phase"
    extends FlowController;

    import Thermo.rho;
    import Thermo.mw;
    import Thermo.Species;
    import Thermo.Phases;
    import Modelica.SIunits.Temperature;

    parameter Temperature T_ref = 273.15 "Reference temperature";

    annotation (defaultComponentName="mfc", Icon(coordinateSystem(
              preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                                                 graphics),
                                                      Documentation(info="<html>
<p>This class implements a mass flow controller with volumetric units (\"field units\").
Since there are at least <em>two</em> different standards (the one named \"standard\" at 
0 Celsius and the one named \"norm\" at 70 Fahrenheit or 21.111... Celsius), it is 
necessary to provide the reference temperature; the default assumes zero Celsius 
(the so-called \"standard\" value).</p>
 
<p>The flow assumes that all components are in gas phase and takes their density from
the Thermo library, where the ideal gas law is (usually) assumed.</p>
</html>"),
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                100,100}}),
              graphics));
    IO.VolumeFlowRateInput V "Volumetric flow rate" 
      annotation (Placement(transformation(extent={{110,-10},{90,10}}, rotation=
               0)));

  equation
    V = sum({inlet.n[i] * mw(i) / rho(T_ref, i, Phases.Gas) for i in Species});

      connect(outlet, inlet) annotation (Line(points={{5.55112e-16,100},{
              5.55112e-16,75},{5.55112e-16,75},{5.55112e-16,50},{5.55112e-16,
              5.55112e-16},{5.55112e-16,5.55112e-16}}, color={0,127,127}));
  end GasFlowController;

  model LiquidPump "A pump, only for liquid phase"
    extends FlowController;

      import Thermo.rho;
      import Thermo.mw;
      import Thermo.h;
      import Thermo.Phases;
      import Thermo.Condensables;
      import Modelica.SIunits.VolumeFlowRate;
      import Units.Temperature;

    annotation (Icon(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
                -100},{100,100}}), graphics={Polygon(
              points={{-100,-100},{-80,-40},{80,-40},{100,-100},{-100,-100}},
              smooth=Smooth.None,
              pattern=LinePattern.None,
              lineColor={0,0,0},
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid), Ellipse(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}),
                             Documentation(info="<html>
<p>This class implements a liquid pump.</p>
<p>The pump takes density values from the Thermo library, and assumes only water and methanol
are present and in liquid phase.</p>
<p>With the appropriate degrees of freedom removed, it can work as a flow measurement as well.</p>
</html>"),
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                100,100}}),
              graphics));
    IO.VolumeFlowRateInput V "Volumetric flow rate" 
      annotation (Placement(transformation(extent={{110,-10},{90,10}}, rotation=
               0)));

    Temperature T;

  equation
    V = sum(inlet.n[i] * mw(i) / rho(T, i, Phases.Liquid) for i in Condensables);
    inlet.H = sum(inlet.n[i] * h(T, i, Phases.Liquid) for i in Condensables);

  end LiquidPump;

  model PeristalticPump "A pump for two-phase flow"
    extends FlowController;

    import Thermo.rho;
    import Thermo.mw;
    import Thermo.h;
    import Thermo.Phases;
    import Thermo.Condensables;
    import Thermo.Species;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={Polygon(
              points={{-100,-100},{-80,-40},{80,-40},{100,-100},{-100,-100}},
              smooth=Smooth.None,
              pattern=LinePattern.None,
              lineColor={0,0,0},
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid), Ellipse(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}),
                             Documentation(info="<html>
<p>This class implements a peristaltic pump, which can handle two phases.</p>
<p>The pump takes density values from the Thermo library.</p>
<p>With the appropriate degrees of freedom removed, it can work as a flow measurement as well.</p>
</html>"),
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                100,100}}), graphics={Text(
              extent={{-20,60},{20,40}},
              lineColor={0,0,0},
              textString="Deactivated in code")}));
    FlowTemperature T annotation (Placement(transformation(
            origin={-40,42},
            extent={{-10,-10},{10,10}},
            rotation=90)));

  equation
    // Deactivate old connection
    // NOTE that inlet is an _outside_ connector, so it needs a minus.
    -inlet.n + T.inlet.n = zeros(size(inlet.n,1));
    -inlet.H + T.inlet.H = 0;

    V = sum(T.liquid[i]*mw(i)/rho(T.T,i,Phases.Liquid) for i in Condensables) +
        sum(T.vapour[i]*mw(i)/rho(T.T,i,Phases.Gas)    for i in Species);

    connect(inlet, T.inlet) annotation (Line(points={{5.55112e-16,5.55112e-16},
              {-20,0},{-40,0},{-40,34}}, color={0,127,127}));
    connect(T.outlet, outlet) annotation (Line(points={{-40,50},{-40,100},{
              5.55112e-16,100}}, color={0,127,127}));
  end PeristalticPump;

    model FlowTemperature "Calculates a flow's temperature"
      extends Modelica.Icons.RotationalSensor;

      import Modelica.SIunits.MoleFraction;
      import Thermo.Incondensables;
      import Thermo.Condensables;
      import Thermo.Species;
      import Thermo.Phases;
      import Thermo.h;
      import Thermo.rr;
      import Thermo.K;
      import Units.MolarFlow;

      annotation (defaultComponentName="T",Diagram(graphics),
                                                    Icon(graphics={Text(
              extent={{-100,100},{100,140}},
              lineColor={0,0,0},
              textString="%name")}),
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
                     annotation (Placement(transformation(extent={{-90,-10},{
                -70,10}}, rotation=0)));
      FlowPort outlet "Exiting flow" 
                     annotation (Placement(transformation(extent={{70,-10},{90,
                10}}, rotation=0)));

      MoleFraction beta "Vapour fraction";

      MolarFlow[Species] vapour "Gas-phase flows";
      MolarFlow[Species] liquid "Liquid-phase flows";

    protected
      MoleFraction z_m(start=0) "Overall methanol molar fraction";
      MoleFraction z_w(start=0.1) "Overall water molar fraction";

    public
      IO.TemperatureOutput T "The flow temperature" 
      annotation (Placement(transformation(
            origin={0,-80},
            extent={{-10,10},{10,-10}},
            rotation=270)));
    equation
      // Liquid and vapour sum to overall flow
      liquid + vapour = inlet.n;
      // Condensable species are present in both phases
      for i in Condensables loop
        vapour[i] * (1 + beta*(K(T,i) -1)) = inlet.n[i]*beta*K(T,i);
      end for;
      // Incondensable species are only in gas phase
      liquid[Incondensables] = zeros(size(Incondensables,1));

      /* Enforce definition of molar fractions. Avoid using divisions
   * in order to avoid division-by-zero errors. */
      inlet.n[Species.Methanol] = sum(inlet.n) * z_m;
      inlet.n[Species.Water]    = sum(inlet.n) * z_w;

      // Combine enthalpic flow with temperature, thermodynamic data and material flows.
      inlet.H = sum( vapour[i] * h(T, i, Phases.Gas)    for i in Species)
              + sum( liquid[i] * h(T, i, Phases.Liquid) for i in Condensables);

      /* Use the beta returned by Thermo.rr only if:
   * 1) there can be a phase equilibrium at all at this temperature;
   * 2) the returned beta is less than 1. */
      if K(T,Species.Water) >= 1 or rr(z_m, z_w, T) >= 1 then
        beta = 1;
      else
        beta = rr(z_m, z_w, T);
      end if;

      connect(inlet, outlet) 
                          annotation (Line(points={{-80,5.55112e-16},{0,
              -4.87687e-22},{0,5.55112e-16},{80,5.55112e-16}}, color={0,127,127}));
    end FlowTemperature;

    model FlowConcentration "Calculates a liquid flow's methanol concentration"
      extends FlowTemperature;

      import Thermo.Condensables;
      import Thermo.Phases;
      import Thermo.Species;
      import Thermo.rho;
      import Thermo.mw;

      annotation (defaultComponentName="TC",Diagram(graphics),
                                                     Icon(graphics),
      Documentation(defaultComponentName="c", info="<html>
<p>Adds the ability to read the methanol concentration <em>in the liquid 
phase</em> to the temperature measurement of <tt>FlowTemperature</tt>.</p>
<p>If there is no liquid flow, then the reported value is zero.</p>
</html>"));
      IO.ConcentrationOutput c 
                          annotation (Placement(transformation(
            origin={0,80},
            extent={{-10,-10},{10,10}},
            rotation=90)));

    equation
      // Methanol flow is concentration times volumetric flow
      liquid[Species.Methanol] = c * sum(liquid[i]*mw(i)/rho(T,i,Phases.Liquid) for i in Condensables);

    end FlowConcentration;

    model MethanolInAir "A measurement of methanol in a gas stream"
      extends Modelica.Icons.TranslationalSensor;

      import Modelica.SIunits.Density;
      import Modelica.Constants.R;
      import Thermo.Species;
      import Thermo.mw;

      outer Modelica.SIunits.Pressure p_env "Environment pressure";

      Density exposure "Methanol mass in gas volume";

      Real chronic_toxicity = exposure / chronic "Chronic toxicity limit is 1";
      Real acute_toxicity =   exposure / acute "Acute toxicity limit is 1";
      Real Mak_toxicity =   exposure / MAK "MAK-Wert, limit is 1";

      FlowPort inlet annotation (Placement(transformation(extent={{-80,-30},{-60,-10}}),
            iconTransformation(extent={{-80,-30},{-60,-10}})));
      FlowPort outlet annotation (Placement(transformation(extent={{60,-30},{80,-10}}),
            iconTransformation(extent={{60,-30},{80,-10}})));
      annotation (Icon(graphics={Text(
              extent={{-100,80},{100,20}},
              lineColor={0,0,0},
              textString="%name")}),
                                  Diagram(graphics),
        Documentation(info="<html>
<p>This unit measures the concentration of methanol in air of a flow passing through,
and helps to draw conclusions on its toxicity.</p>
<p>Three values are provided, <tt>chronic_toxicity</tt>, <tt>acute_toxicity</tt>, and
<tt>Mak_toxicity</tt> which at less than unitary values indicate a concentration that is 
safe respectively permanently, for a period of one hour, or for German workplaces.</p>
<p>The numbers can also be read as the dilution required to bring the toxicity to
acceptable levels. E.g. a value of 10 indicates that the gas volume has to be diluted
10 times to become safe for breathing.</p>
<h3>Effects of Poisoning</h3>
<p>Acute methanol poisoning on humans is known and affects mainly the nervous system.
Chronic poisoning is extrapolated by tests on rats, and is assumed to induce 
teratogenicity.</p>
<h3>Source</h3>
<p>The given concentration limits are 4&nbsp;mg/m<sup>3</sup> for chronic and 
28&nbsp;mg/m<sup>3</sup> for acute exposure; data were taken from 
<a href=\"http://www.oehha.org/air/allrels.html\">the Office of Environmental Health Hazard 
Assessment</a> of the state of California.
The MAK value was taken from <a href=\"http://www.enius.de/schadstoffe/methanol.html\">this
link</a>; note that the link says 260&nbsp;mg/m<sup>3</sup>, but there are other numbers
reported. An official one should be sought.</p>

</html>"));
    protected
      constant Density chronic =  4E-6 "Chronic reference exposure level";
      constant Density acute = 28E-6 "Acute (1 h) reference exposure level";
      constant Density MAK = 260E-6 "Maximale Arbeitsplatzkonzentration";

      FlowTemperature T 
        annotation (Placement(transformation(extent={{-12,-32},{12,-8}})));
    equation

      exposure = T.vapour[Species.Methanol] * mw(Species.Methanol) / (sum(T.vapour) * R * T.T / p_env);

      connect(inlet, T.inlet) annotation (Line(
          points={{-70,-20},{-9.6,-20}},
          color={0,127,127},
          smooth=Smooth.None));
      connect(T.outlet, outlet) annotation (Line(
          points={{9.6,-20},{70,-20}},
          color={0,127,127},
          smooth=Smooth.None));
    end MethanolInAir;

    package Test
      model LiquidPumpTest

        Sources.Solution source 
          annotation (Placement(transformation(extent={{-100,-10},{-80,10}},
                rotation=0)));
        Sink sink annotation (Placement(transformation(extent={{60,0},{80,20}},
                rotation=0)));
        Flow.Measurements.LiquidPump pump annotation (Placement(transformation(
                extent={{-10,-10},{10,10}}, rotation=0)));
        Measurements.FlowTemperature T annotation (Placement(transformation(
                extent={{-40,-10},{-20,10}}, rotation=0)));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                  -100},{100,100}}),
                            graphics));

        inner parameter Modelica.SIunits.Pressure p_env = 101325;
        inner parameter Units.Temperature T_env = 298.15;
        inner parameter Units.RelativeHumidity RH_env = 60;

        parameter Modelica.SIunits.VolumeFlowRate V = 1;

      equation
        pump.V = V;

        connect(pump.outlet, sink.inlet) annotation (Line(points={{6.10623e-16,
                10},{30.5,10},{30.5,10},{61,10}}, color={0,127,127}));
        connect(T.outlet, pump.inlet) annotation (Line(points={{-22,6.10623e-16},
                {-16,-3.36456e-22},{-16,6.10623e-16},{6.10623e-16,6.10623e-16}},
              color={0,127,127}));
        connect(T.inlet, source.outlet) annotation (Line(points={{-38,
                6.10623e-16},{-70,6.10623e-16},{-70,6.66134e-16},{-90,
                6.66134e-16}}, color={0,127,127}));
      end LiquidPumpTest;

      model PeristalticPumpTest

      protected
        Sources.Solution source(T=320) 
          annotation (Placement(transformation(extent={{-80,-30},{-60,-10}},
                rotation=0)));
        Sink sink annotation (Placement(transformation(extent={{46,6},{54,14}},
                rotation=0)));
      public
        Flow.Measurements.PeristalticPump pump 
                                          annotation (Placement(transformation(
                extent={{10,-10},{30,10}}, rotation=0)));
        Measurements.FlowTemperature T annotation (Placement(transformation(
                extent={{-20,-10},{0,10}}, rotation=0)));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                  -100},{100,100}}),
                            graphics));

        inner parameter Modelica.SIunits.Pressure p_env = 101325;
        inner parameter Units.Temperature T_env = 298.15;
        inner parameter Units.RelativeHumidity RH_env = 60;

        parameter Modelica.SIunits.VolumeFlowRate V = 0.01;
        parameter Real l_to_g_molratio = 0.5;

      protected
        Sources.Environment env annotation (Placement(transformation(extent={{
                  -80,10},{-60,30}}, rotation=0)));
      equation
        pump.V = V;
        sum(source.outlet.n) / sum(env.outlet.n) = time;

        connect(pump.outlet, sink.inlet) annotation (Line(points={{20,10},{33.2,
                10},{46.4,10}},           color={0,127,127}));
        connect(T.outlet, pump.inlet) annotation (Line(points={{-2,6.10623e-16},
                {4,-3.36456e-22},{4,6.10623e-16},{20,6.10623e-16}}, color={0,
                127,127}));
        connect(T.inlet, source.outlet) annotation (Line(points={{-18,
                6.10623e-16},{-40,6.10623e-16},{-40,-20},{-70,-20}}, color={0,
                127,127}));
        connect(env.outlet, T.inlet) annotation (Line(points={{-61,20},{-40,20},
                {-40,6.10623e-16},{-18,6.10623e-16}}, color={0,127,127}));
      end PeristalticPumpTest;

      model FlowTemperatureTest "A test case for the temperature sensor"

        replaceable Measurements.FlowTemperature measurement 
                                        annotation (Placement(transformation(
                extent={{-20,0},{0,20}}, rotation=0)));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent=
                  {{-100,-100},{100,100}}),
                            graphics));
      protected
        Flow.Sources.Environment env "Atmospheric air" 
                                        annotation (Placement(transformation(
                extent={{-92,20},{-72,40}}, rotation=0)));
        Flow.Sink sink "Dumpster" 
                          annotation (Placement(transformation(extent={{40,6},{
                  48,14}}, rotation=0)));
      public
        Flow.Sources.Solution solution "Source of methanol solution" 
                                          annotation (Placement(transformation(
                extent={{-80,-10},{-68,2}}, rotation=0)));
        inner parameter Modelica.SIunits.Pressure p_env = 101325;
        inner parameter Units.Temperature T_env = 298.15;
        inner Units.RelativeHumidity RH_env = time;

      protected
        constant Modelica.SIunits.Time second = 1 "To get adimensional time";
      equation
        /* Running from time 0 to 1 will test negative flows and
   * crossing of zero-flow condition. */
        sum(env.outlet.n) = -10*(time/second-0.5);
        sum(solution.outlet.n) = -(time/second-0.5);

        connect(sink.inlet, measurement.outlet) 
          annotation (Line(points={{40.4,10},{-2,10}}, color={0,127,127}));
        connect(env.outlet, measurement.inlet) 
                                             annotation (Line(points={{-73,30},
                {-60,30},{-60,10},{-18,10}}, color={0,127,127}));
        connect(solution.outlet, measurement.inlet) annotation (Line(points={{
                -74,-4},{-60,-4},{-60,10},{-18,10}}, color={0,127,127}));
      end FlowTemperatureTest;

      model FlowConcentrationTest
        "Test case for the more detailed concentration sensor"
        extends Flow.Measurements.Test.FlowTemperatureTest(redeclare
            Measurements.FlowConcentration measurement);
      end FlowConcentrationTest;

      model GasMethanolTest "A test case for the methanol in air sensor"

        annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}),
                            graphics),
          experiment,
          experimentSetupOutput);
      protected
        Flow.Sources.Environment env "Atmospheric air" 
                                        annotation (Placement(transformation(
                extent={{-92,20},{-72,40}}, rotation=0)));
        Flow.Sink sink "Dumpster" 
                          annotation (Placement(transformation(extent={{40,6},{48,14}},
                           rotation=0)));
      public
        Flow.Sources.Solution solution "Source of methanol solution" 
                                          annotation (Placement(transformation(
                extent={{-80,-10},{-68,2}}, rotation=0)));
        inner parameter Modelica.SIunits.Pressure p_env = 101325;
        inner Units.Temperature T_env;
        inner Units.RelativeHumidity RH_env = 60;

      protected
        constant Modelica.SIunits.Time second = 1 "To get adimensional time";
      public
        MethanolInAir methanolInAir 
          annotation (Placement(transformation(extent={{-20,2},{0,22}})));
      equation
        /* Running from time 0 to 1 will test negative flows and
   * crossing of zero-flow condition. */
        sum(env.outlet.n) = -1;
        sum(solution.outlet.n) = -1;
        T_env = 300 + 80*time;

        connect(methanolInAir.inlet, solution.outlet) annotation (Line(
            points={{-17,10},{-44,10},{-44,-4},{-74,-4}},
            color={0,127,127},
            smooth=Smooth.None));
        connect(methanolInAir.inlet, env.outlet) annotation (Line(
            points={{-17,10},{-44,10},{-44,30},{-73,30}},
            color={0,127,127},
            smooth=Smooth.None));
        connect(methanolInAir.outlet, sink.inlet) annotation (Line(
            points={{-3,10},{40.4,10}},
            color={0,127,127},
            smooth=Smooth.None));
      end GasMethanolTest;
    end Test;
  end Measurements;

  package UnitOperations "Unit operations"
    model Mixer "A unit mixing four molar flows."

      import Modelica.SIunits.AmountOfSubstance;
      import Modelica.SIunits.Concentration;
      import Modelica.SIunits.InternalEnergy;
      import Modelica.SIunits.Volume;
      import Thermo.Species;
      import Thermo.Incondensables;
      import Thermo.Condensables;
      import Thermo.Phases;
      import Thermo.h;
      import Thermo.mw;
      import Thermo.rho;

      FlowPort outlet "The mixer's outlet" 
                            annotation (Placement(transformation(extent={{-90,
                -10},{-70,10}}, rotation=0)));
      FlowPort fuelInlet "The methanol-feed inlet" 
                             annotation (Placement(transformation(extent={{-10,-10},
                {10,10}},       rotation=180,
            origin={0,-80})));
      FlowPort loopInlet "The anode loop's inlet" 
                             annotation (Placement(transformation(extent={{-10,
                70},{10,90}}, rotation=0)));
      FlowPort waterInlet "The water-recovery inlet" 
                             annotation (Placement(transformation(extent={{70,
                -10},{90,10}}, rotation=0)));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
                -100},{100,100}}),
                          graphics),
                           Icon(coordinateSystem(preserveAspectRatio=false,
              extent={{-100,-100},{100,100}}), graphics={
            Ellipse(
              extent={{-80,80},{80,-80}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Rectangle(
              extent={{-80,0},{80,80}},
              lineColor={0,0,255},
              pattern=LinePattern.None,
              lineThickness=1,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Line(points={{-80,0},{-80,80},{80,80},{80,0}}, color={0,0,0}),
            Line(points={{0,6},{0,-54}}, color={0,0,0}),
            Line(points={{0,6},{-52,36}}, color={0,0,0}),
            Line(points={{0,6},{52,36}}, color={0,0,0}),
            Text(
              extent={{-100,160},{100,100}},
              lineColor={0,0,0},
              fillColor={255,85,85},
              fillPattern=FillPattern.Solid,
              textString="%name")}),
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

      AmountOfSubstance n[Species] "Molar holdup";
      InternalEnergy U "Energy holdup";
      Concentration c(start=1000) "Methanol concentration";
      Volume V(start=5E-6, fixed=true) "Solution volume";

      IO.TemperatureOutput T "Mixer temperature" 
        annotation (Placement(transformation(extent={{-80,60},{-100,80}},
              rotation=0)));
    equation
      der(U) = fuelInlet.H + loopInlet.H + waterInlet.H + outlet.H;
      der(n) = fuelInlet.n + loopInlet.n + waterInlet.n + outlet.n;

      // Bind outlet's n to composition in holdup
      outlet.n[1:end-1] / sum(outlet.n) = n[1:end-1] / sum(n);
      // Bind outlet's H to specific internal energy and outlet flow
      outlet.H / sum(outlet.n) = U / sum(n);

      U = sum(n[i]*h(T, i, Phases.Liquid) for i in Condensables);
      V = sum(n[i]*mw(i)/rho(T, i, Phases.Liquid) for i in Condensables);
      c = n[Species.Methanol] / V;

      assert(V > sqrt(Modelica.Constants.eps), "==> Mixer ran out of solution");

    initial equation
      n[Incondensables] = zeros(size(Incondensables,1));

    end Mixer;

    model HydrostaticMixer
      "A mixer producing a hydrostatic pressure for measurement"
      extends Mixer;

      import Thermo.mw;
      import Thermo.Condensables;
      import g = Modelica.Constants.g_n;

      Modelica.SIunits.Mass m;

      IO.PressureOutput p
        "Hydrostatic pressure measured at the bottom of the mixer" 
                       annotation (Placement(transformation(extent={{-54,-78},{
                -40,-60}}, rotation=0)));

      parameter Modelica.SIunits.Area A = 50E-4 "Mixer cross-sectional area";

    equation
      m = sum( n[i]*mw(i) for i in Condensables);
      p = m * g / A;

      annotation (Documentation(info="<html>
<p>This extension to the standard-issue Mixer class adds the possibility of inferring
the amount of solution by measuring the hydrostatic pressure in the tank.</p>
</html>"));
    end HydrostaticMixer;

    model ElasticMixer "Mixer with capability for volume expansion"
      extends Mixer;

      import Modelica.SIunits.Length;

      parameter Modelica.SIunits.Volume V0 = 5E-6
        "Maximum nominal volume at no pressure";
      parameter Modelica.SIunits.Stress E = 1.9E6 "Young modulus";
      parameter Length d_in = 3.2E-3 "Tubing internal diameter";
      parameter Length d_out = 4.8E-3 "Tubing external diameter";

      PressurePort pressure "The mixer's internal pressure" annotation (Placement(
            transformation(extent={{68,68},{88,88}}), iconTransformation(extent={{68,
                68},{88,88}})));
    equation
      pressure.p = max(0, (d_out/d_in - 1)*E*((V/V0)^(1/3) - 1));

      annotation (Diagram(graphics), Icon(coordinateSystem(preserveAspectRatio=true,
              extent={{-100,-100},{100,100}}), graphics),
        Documentation(info="<html>
<p>This mixer is a sort of \"balloon\" that can increase its volume when pressure is
increased.</p>
<p>It is assumed through a stretch of imagination that the tubing acts as a CSTR more than
a PFR. Take concentration values with a grain of salt.</p>
<p>The default numbers are related to the values for Tygon 3350 tubing, where
an overpressure of 10 kPa results in a 3% volume increase.</p>
</html>"));
    end ElasticMixer;

    model AbstractSeparator "Generic separator"

      import Thermo.Condensables;
      import Thermo.Phases;
      import Thermo.h;

      FlowPort inlet "Two-phase inlet" 
                     annotation (Placement(transformation(extent={{-110,-10},{
                -90,10}}, rotation=0)));
      FlowPort gasOutlet "Single-phase gas outlet" 
                         annotation (Placement(transformation(extent={{60,30},{
                80,50}}, rotation=0)));
      FlowPort liquidOutlet "Single-phase liquid outlet" 
                            annotation (Placement(transformation(extent={{-10,-10},{
                10,10}},        rotation=180,
            origin={70,-40})));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={
            Ellipse(
              extent={{-100,40},{-60,-40}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Ellipse(
              extent={{60,40},{100,-40}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Rectangle(
              extent={{-80,40},{80,-40}},
              lineColor={255,255,255},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Line(points={{-80,40},{80,40}}, color={0,0,0}),
            Line(points={{-80,-40},{80,-40}}, color={0,0,0}),
            Text(
              extent={{-100,100},{100,40}},
              lineColor={0,0,0},
              fillColor={255,85,85},
              fillPattern=FillPattern.Solid,
              textString="%name")}),
                              Diagram(coordinateSystem(preserveAspectRatio=true,
                       extent={{-100,-100},{100,100}}),
                                      graphics),
        Documentation(info="<html>
<p>An abstract representation of the separation of a flow.</p>
</html>"));
      IO.TemperatureOutput T "Separator temperature" 
                                                  annotation (Placement(
            transformation(extent={{100,-10},{120,10}}, rotation=0)));
    protected
      Measurements.FlowTemperature ft 
                         annotation (Placement(transformation(extent={{-10,-10},
                {10,10}}, rotation=0)));

    equation
      liquidOutlet.H = sum(h(T, i, Phases.Liquid) * liquidOutlet.n[i] for i in Condensables);
      connect(ft.inlet, inlet) annotation (Line(points={{-8,6.10623e-16},{-54,
              6.10623e-16},{-54,5.55112e-16},{-100,5.55112e-16}}, color={0,127,
              127}));
      connect(ft.outlet, gasOutlet) annotation (Line(points={{8,6.10623e-16},{
              40,6.10623e-16},{40,40},{70,40}}, color={0,127,127}));
      connect(ft.outlet, liquidOutlet) annotation (Line(points={{8,6.10623e-16},
              {40,6.10623e-16},{40,-40},{70,-40}}, color={0,127,127}));
      connect(T, ft.T) annotation (Line(points={{110,5.55112e-16},{78,0},{60,0},
              {60,-20},{0,-20},{0,-8},{6.10623e-16,-8}},
                                                       color={0,0,255}));
    end AbstractSeparator;

    model Separator "Splits a flow in two parts"
      extends AbstractSeparator;

    equation
      liquidOutlet.n = -ft.liquid;

      annotation (Documentation(info="<html>
<p>The separator unit simply splits a flow in its gaseous and liquid components. 
The separation criterion is straightforwardly the liquid-vapor equilibrium.</p>
</html>"));
    end Separator;

    model CapillarySeparator "Separates a flow depending on backpressure"
      extends AbstractSeparator;

      import Modelica.SIunits.Pressure;
      import Thermo.rho;
      import Thermo.mw;
      import Thermo.Phases;
      import Thermo.Species;
      import Thermo.Condensables;

      parameter Modelica.SIunits.Length dh_flow = 200E-6
        "Hydraulic diameter of hydrophobic channels";
      parameter Modelica.SIunits.Length dh_sep = 10E-6
        "Hydraulic diameter of hydrophilic channels";
      parameter Modelica.SIunits.Angle theta_flow = 2.2
        "Contact angle of hydrophobic channels";
      parameter Modelica.SIunits.Angle theta_sep = 0.35
        "Contact angle of hydrophilic channels";

      Pressure pc_flow "Capillary pressure in hydrophobic channels";
      Pressure pc_sep "Capillary pressure in hydrophilic channels";
      Modelica.SIunits.SurfaceTension sigma "Surface tension of water with air";

      Real R_gl "Gas-liquid volumetric ratio";
      Real recovery "Fraction of lost liquid water";

      PressurePort backPressure "Pressure from the liquid outlet" annotation (
          Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={40,-40}), iconTransformation(extent={{30,-50},{50,-30}})));
    protected
      Real fuzzifier "Continuous change between false (0) and true (1)";
      parameter Pressure p_eps = 5 "Small value for pressure";

    equation
      sigma = 0.076 - 0.00017*(T-273.15); // From Microfluidics, it is in Celsius there!
      pc_flow = - 4 * sigma / dh_flow * cos(theta_flow);
      pc_sep  = - 4 * sigma / dh_sep  * cos(theta_sep);

      R_gl = sum( ft.vapour[i] * mw(i) / rho(ft.T, i, Phases.Gas) for i in Species)  /
             sum( ft.liquid[i] * mw(i) / rho(ft.T, i, Phases.Liquid) for i in Condensables);

      recovery = -liquidOutlet.n[Species.Water] / ft.liquid[Species.Water];

      fuzzifier = min(1, max(0, (pc_flow - backPressure.p + p_eps)/(2*p_eps)));

      liquidOutlet.n = -ft.liquid * fuzzifier;

      connect(ft.inlet, inlet) annotation (Line(points={{-8,6.10623e-16},{-54,
              6.10623e-16},{-54,5.55112e-16},{-100,5.55112e-16}}, color={0,127,
              127}));
      connect(ft.outlet, gasOutlet) annotation (Line(points={{8,6.10623e-16},{
              40,6.10623e-16},{40,40},{70,40}}, color={0,127,127}));
      connect(ft.outlet, liquidOutlet) annotation (Line(points={{8,6.10623e-16},
              {40,6.10623e-16},{40,-40},{70,-40}}, color={0,127,127}));
      connect(T, ft.T) annotation (Line(points={{110,5.55112e-16},{78,0},{60,0},
              {60,-20},{0,-20},{0,-8},{6.10623e-16,-8}},
                                                       color={0,0,255}));
      annotation (Diagram(graphics), Documentation(info="<html>
<p>This separator will produce a liquid phase only if the capillary pressure
difference is higher than the backpressure.</p>
</html>"));
    end CapillarySeparator;

    model Burner "An adiabatic combustor"

      import Units.MolarFlow;
      import Thermo.Species;

      FlowPort inlet "Burner inlet" annotation (Placement(transformation(extent=
               {{-108,-10},{-88,10}}, rotation=0)));
      FlowPort outlet "Burner outlet" annotation (Placement(transformation(
              extent={{92,-10},{112,10}}, rotation=0)));
      annotation (
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}),
                graphics),
        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                100,100}}), graphics={
            Ellipse(
              extent={{-100,40},{-20,-40}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Ellipse(
              extent={{20,40},{100,-40}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Rectangle(
              extent={{-60,40},{60,-40}},
              lineColor={255,255,255},
              lineThickness=0.5,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Rectangle(
              extent={{-40,40},{40,-40}},
              lineColor={0,0,255},
              pattern=LinePattern.None,
              lineThickness=0.5,
              fillColor={255,128,0},
              fillPattern=FillPattern.CrossDiag),
            Line(points={{-60,40},{60,40},{60,40}}, color={0,0,0}),
            Line(points={{-60,-40},{60,-40},{60,-40}}, color={0,0,0})}),
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
        annotation (Placement(transformation(extent={{40,-10},{60,10}},
              rotation=0)));
    public
      IO.TemperatureOutput T_out(start=373.15) "Temperature after combustion" 
        annotation (Placement(transformation(
            origin={50,-50},
            extent={{-10,-10},{10,10}},
            rotation=270)));
    protected
      Flow.Sink sink annotation (Placement(transformation(extent={{0,-50},{20,
                -30}}, rotation=0)));
    equation
      sink.inlet.H = 0;
      sink.inlet.n[Species.Methanol] = reaction;
      sink.inlet.n[Species.Water] = -2*reaction;
      sink.inlet.n[Species.Oxygen] = 1.5*reaction;
      sink.inlet.n[Species.CarbonDioxide] = -reaction;
      sink.inlet.n[Species.Nitrogen] = 0;

      reaction = min(inlet.n[Species.Methanol], inlet.n[Species.Oxygen]/1.5);

      connect(T_out, T.T) annotation (Line(points={{50,-50},{50,-8}}, color={0,
              0,255}));
      connect(T.outlet, outlet) annotation (Line(points={{58,6.10623e-16},{61,
              6.10623e-16},{61,5.55112e-16},{102,5.55112e-16}}, color={0,127,
              127}));
      connect(T.inlet, inlet) annotation (Line(points={{42,6.10623e-16},{-48,
              6.10623e-16},{-48,5.55112e-16},{-98,5.55112e-16}}, color={0,127,
              127}));
      connect(sink.inlet, inlet) annotation (Line(points={{1,-40},{-20,-40},{
              -20,5.55112e-16},{-98,5.55112e-16}}, color={0,127,127}));
    end Burner;

    package HeatExchangers "Various types of heat exchangers"

    partial model Abstract "An abstract heat exchanger"

      FlowPort hot_1 "Port for hot flow on side 1" 
                     annotation (Placement(transformation(extent={{-100,20},{
                  -80,40}}, rotation=0)));
      FlowPort hot_2 "Port for hot flow on side 2" 
                      annotation (Placement(transformation(extent={{-100,-40},{
                  -80,-20}}, rotation=0)));
      annotation (defaultComponentName="exchanger", Diagram(coordinateSystem(
                preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                                                            graphics),
                                                             Icon(
              coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                  100,100}}), graphics={
              Rectangle(
                extent={{-100,40},{100,-40}},
                lineColor={0,0,0},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Text(
                extent={{-100,100},{100,40}},
                lineColor={0,0,0},
                fillColor={255,85,85},
                fillPattern=FillPattern.Solid,
                textString="%name"),
              Line(points={{-60,40},{-60,-40}}, color={0,0,0}),
              Line(points={{20,40},{20,-40}}, color={0,0,0}),
              Line(points={{-20,40},{-20,-40}}, color={0,0,0}),
              Line(points={{60,40},{60,-40}}, color={0,0,0}),
              Line(points={{-40,40},{-40,-40}}, color={0,0,0}),
              Line(points={{40,40},{40,-40}}, color={0,0,0}),
              Line(points={{0,40},{0,-40}}, color={0,0,0}),
              Line(points={{80,40},{80,-40}}, color={0,0,0}),
              Line(points={{-80,40},{-80,-40}}, color={0,0,0})}),
        Documentation(info="<html>
<p>This is the interface for a generic heat exchanger. It constraints the two sides to
be able to exchange a heat duty <tt>Q</tt>, but no mass.</p>
<p>All subclasses may be counter-current or co-current exchangers depending on the
connections of the flows to the two sides.</p>
</html>"));
      FlowPort cold_1 "Port for cold flow on side 2" 
                     annotation (Placement(transformation(extent={{80,20},{100,
                  40}}, rotation=0)));
      FlowPort cold_2 "Port for cold flow on side 2" 
                     annotation (Placement(transformation(extent={{80,-40},{100,
                  -20}}, rotation=0)));
      IO.TemperatureOutput T_hot_1 "Temperature of hot stream on side 1" 
        annotation (Placement(transformation(extent={{-80,40},{-100,60}},
                rotation=0)));
      IO.TemperatureOutput T_hot_2 "Temperature of hot stream on side 2" 
        annotation (Placement(transformation(extent={{-80,-60},{-100,-40}},
                rotation=0)));
      IO.TemperatureOutput T_cold_1 "Temperature of cold stream on side 1" 
        annotation (Placement(transformation(extent={{80,40},{100,60}},
                rotation=0)));
      IO.TemperatureOutput T_cold_2 "Temperature of cold stream on side 2" 
        annotation (Placement(transformation(extent={{80,-60},{100,-40}},
                rotation=0)));

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
        annotation (Placement(transformation(extent={{-60,40},{-40,20}},
                rotation=0)));
      FlowTemperature cold_1_T "Temperature for cold flow on side 1" 
        annotation (Placement(transformation(extent={{40,40},{60,20}}, rotation=
                 0)));
      FlowTemperature hot_2_T "Temperature for hot flow on side 2" 
        annotation (Placement(transformation(extent={{-60,-40},{-40,-20}},
                rotation=0)));
      FlowTemperature cold_2_T "Temperature for cold flow on side 2" 
        annotation (Placement(transformation(extent={{40,-40},{60,-20}},
                rotation=0)));
      annotation (Diagram(graphics),
                           DymolaStoredErrors,
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
                        annotation (Placement(transformation(extent={{20,-80},{
                  40,-60}}, rotation=0)));

    equation
      Q = U * A * LMTD;

      if noEvent(abs(DT1) > eps and abs(DT2) > eps) then
        LMTD = (DT1-DT2)/log(DT1/DT2);
      else
        LMTD = 0;
      end if;

      connect(hot_2_T.outlet, sinkPort.inlet)    annotation (Line(points={{-42,
                -30},{0,-30},{0,-70},{21,-70}}, color={0,127,127}));
      connect(cold_2_T.inlet, sinkPort.inlet)   annotation (Line(points={{42,
                -30},{0,-30},{0,-70},{21,-70}}, color={0,127,127}));
      connect(cold_1_T.inlet, sinkPort.inlet)  annotation (Line(points={{42,30},
                {0,30},{0,-70},{21,-70}}, color={0,127,127}));
      connect(hot_1_T.outlet, sinkPort.inlet)   annotation (Line(points={{-42,
                30},{0,30},{0,-70},{21,-70}}, color={0,127,127}));
      connect(cold_1, cold_1_T.outlet) 
        annotation (Line(points={{90,30},{70,30},{70,30},{58,30}}, color={0,127,
                127}));
      connect(cold_2_T.outlet, cold_2) annotation (Line(points={{58,-30},{90,
                -30}}, color={0,127,127}));
      connect(hot_2_T.inlet, hot_2) annotation (Line(points={{-58,-30},{-90,-30}},
              color={0,127,127}));
      connect(hot_1_T.inlet, hot_1) annotation (Line(points={{-58,30},{-90,30}},
              color={0,127,127}));
      connect(hot_1_T.T, T_hot_1) annotation (Line(points={{-50,38},{-50,50},{
                -90,50}}, color={0,0,255}));
      connect(cold_1_T.T, T_cold_1) annotation (Line(points={{50,38},{50,50},{
                90,50}}, color={0,0,255}));
      connect(cold_2_T.T, T_cold_2) annotation (Line(points={{50,-38},{50,-50},
                {90,-50}}, color={0,0,255}));
      connect(hot_2_T.T, T_hot_2) annotation (Line(points={{-50,-38},{-50,-50},
                {-90,-50}}, color={0,0,255}));
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
        annotation (Placement(transformation(extent={{-60,40},{-40,20}},
                rotation=0)));
      FlowTemperature cold_1_T "Temperature of the cold flow on side 1" 
        annotation (Placement(transformation(extent={{40,40},{60,20}}, rotation=
                 0)));
      FlowTemperature hot_2_T "Temperature of the hot flow on side 2" 
        annotation (Placement(transformation(extent={{-60,-40},{-40,-20}},
                rotation=0)));
      FlowTemperature cold_2_T "Temperature of the cold flow on side 2" 
        annotation (Placement(transformation(extent={{40,-40},{60,-20}},
                rotation=0)));
      annotation (defaultComponentName="step", Diagram(graphics),
                                                        Documentation(info="<html>
<p>To implement a discretised heat exchanger, a single step is implemented 
here as a simplistic heat exchanger. It could be based on LMTD, but using
the average temperature is numerically more robust and allows Dymola to
perform algebraic manipulation.</p>
</html>"));
      Flow.Sink sink "Sink element" 
                        annotation (Placement(transformation(extent={{20,-80},{
                  40,-60}}, rotation=0)));

    equation
      Q = U*A*(T_hot-T_cold);

      connect(hot_2_T.inlet, hot_2) annotation (Line(points={{-58,-30},{-90,-30}},
              color={0,127,127}));
      connect(hot_1_T.inlet, hot_1) annotation (Line(points={{-58,30},{-90,30}},
              color={0,127,127}));
      connect(cold_1_T.outlet, cold_1) annotation (Line(points={{58,30},{90,30}},
              color={0,127,127}));
      connect(cold_2_T.outlet, cold_2) annotation (Line(points={{58,-30},{90,
                -30}}, color={0,127,127}));
      connect(sink.inlet, hot_2_T.outlet) annotation (Line(points={{21,-70},{0,
                -70},{0,-30},{-42,-30}}, color={0,127,127}));
      connect(sink.inlet, cold_2_T.inlet) annotation (Line(points={{21,-70},{0,
                -70},{0,-30},{42,-30}}, color={0,127,127}));
      connect(sink.inlet, hot_1_T.outlet) annotation (Line(points={{21,-70},{0,
                -70},{0,30},{-42,30}}, color={0,127,127}));
      connect(sink.inlet, cold_1_T.inlet) annotation (Line(points={{21,-70},{0,
                -70},{0,30},{42,30}}, color={0,127,127}));
      connect(cold_1_T.T, T_cold_1) annotation (Line(points={{50,38},{50,50},{
                90,50}}, color={0,0,255}));
      connect(cold_2_T.T, T_cold_2) annotation (Line(points={{50,-38},{50,-50},
                {90,-50}}, color={0,0,255}));
      connect(hot_2_T.T, T_hot_2) annotation (Line(points={{-50,-38},{-50,-50},
                {-90,-50}}, color={0,0,255}));
      connect(hot_1_T.T, T_hot_1) annotation (Line(points={{-50,38},{-50,50},{
                -90,50}}, color={0,0,255}));
    end DiscretisedStep;

    model Discretised "A heat exchanger made up of many discrete sections"
      extends Flow.UnitOperations.HeatExchangers.Abstract;

      parameter Integer n(min=3) = 10 "Number of discretisation units";

      protected
      Flow.UnitOperations.HeatExchangers.DiscretisedStep[n] steps(each A=A/n,
            each U=U) "The steps the exchanger is divided in";

      protected
      Flow.Sink sink 
                    annotation (Placement(transformation(extent={{20,-80},{40,
                  -60}}, rotation=0)));
      annotation (Diagram(graphics),
                           Documentation(info="<html>
<p>This heat exchanger is much more complex than a LMTD exchanger, and should be
used only for condensing flows where LMTD theory is not valid.</p>
</html>"),
        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                  100,100}}),
             graphics));
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
      connect(sink.inlet, cold_2) annotation (Line(points={{21,-70},{0,-70},{0,
                -30},{90,-30}}, color={0,127,127}));
      connect(sink.inlet, hot_2) annotation (Line(points={{21,-70},{0,-70},{0,
                -30},{-90,-30}}, color={0,127,127}));

    end Discretised;

      package Test
        partial model AbstractHeatExchangerTest
          "Generic test for heat exchangers"

          import Modelica.SIunits.VolumeFlowRate;

          inner parameter Modelica.SIunits.Pressure p_env = 101325;
          inner parameter Units.Temperature T_env = 298.15;
          inner parameter Units.RelativeHumidity RH_env = 60;

        protected
          Sources.Environment env         annotation (Placement(transformation(
                  extent={{100,-60},{80,-40}}, rotation=0)));
          Sink coldSink "Sink for the cold flow" 
                            annotation (Placement(transformation(extent={{62,28},
                    {70,36}}, rotation=0)));
          annotation (Diagram(coordinateSystem(preserveAspectRatio=false,
                  extent={{-100,-100},{100,100}}),
                              graphics),
            experiment(StopTime=80),
            experimentSetupOutput);
        public
          replaceable Abstract exchanger "The heat exchanger" 
                                 annotation (Placement(transformation(extent={{
                    -40,-20},{40,60}}, rotation=0)));
        protected
          Sink hotSink "Sink for the hot flow" 
                            annotation (Placement(transformation(extent={{-28,
                    -52},{-20,-44}}, rotation=0)));
          Sources.Solution methanolSolution(   T=330) 
            annotation (Placement(transformation(extent={{-80,-40},{-60,-20}},
                  rotation=0)));
          Measurements.LiquidPump pump 
                    annotation (Placement(transformation(extent={{-80,0},{-60,
                    20}}, rotation=0)));
        public
          Measurements.GasFlowController mfc 
                                annotation (Placement(transformation(extent={{
                    52,-42},{72,-22}}, rotation=0)));

        public
          parameter VolumeFlowRate air = 100E-3/60 "Full scale of cooling air";
          parameter VolumeFlowRate sol = 10E-6/60 "Loop solution";

        equation
          mfc.V = 0.01*air+0.99*air*time;
          pump.V = sol;

          connect(hotSink.inlet, exchanger.hot_2) annotation (Line(points={{-27.6,
                  -48},{-36,-48},{-36,8}},       color={0,127,127}));
          connect(exchanger.cold_1, coldSink.inlet) annotation (Line(points={{36,32},
                  {49.2,32},{49.2,32},{62.4,32}},        color={0,127,127}));
          connect(env.outlet, mfc.inlet) annotation (Line(points={{81,-50},{62,
                  -50},{62,-32}}, color={0,127,127}));
          connect(mfc.outlet, exchanger.cold_2) annotation (Line(points={{62,-22},
                  {62,8},{36,8}},      color={0,127,127}));
          connect(pump.outlet, exchanger.hot_1) annotation (Line(points={{-70,20},
                  {-70,32},{-36,32}},     color={0,127,127}));
          connect(methanolSolution.outlet, pump.inlet) annotation (Line(points=
                  {{-70,-30},{-70,10}}, color={0,127,127}));
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
                     annotation (Placement(transformation(extent={{-100,-6},{
                  -88,6}}, rotation=0)));
      FlowPort outlet "Outlet from the cooler" 
                      annotation (Placement(transformation(extent={{88,-6},{100,
                  6}}, rotation=0)));
      annotation (defaultComponentName="cooler", Diagram(coordinateSystem(
                preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                                                         graphics),
                                                          Icon(coordinateSystem(
                preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
              graphics={
              Rectangle(
                extent={{-90,20},{90,-20}},
                lineColor={0,0,0},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{-60,16},{-28,-16}},
                lineColor={0,0,0},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{28,16},{60,-16}},
                lineColor={0,0,0},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-44,18},{-26,-18}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                lineThickness=1,
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{26,18},{44,-18}},
                lineColor={0,0,255},
                pattern=LinePattern.None,
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Line(points={{-44,16},{44,-16}}, color={0,0,0}),
              Line(points={{-44,-16},{44,16}}, color={0,0,0}),
              Text(
                extent={{-102,82},{100,20}},
                lineColor={0,0,0},
                fillColor={255,85,85},
                fillPattern=FillPattern.Solid,
                textString="%name")}),
        Documentation(info="<html>
<p>This is the most generic interface for an air cooler, whose heat-exchange model 
is left completely unspecified.</p><p>The cooler provides two flow ports and two temperature
outputs, which will output the temperature of the associate flow; the temperature
input is going to be used by some (yet unspecified) internal system to set the outlet
temperature.</p>
</html>"));
      IO.TemperatureOutput T_process_in "Process inlet temperature" 
                                                                 annotation (Placement(
              transformation(extent={{-88,-20},{-100,-8}}, rotation=0)));
      IO.TemperatureOutput T_process_out "Process outlet temperature" 
                                                                   annotation (Placement(
              transformation(extent={{88,-20},{100,-8}}, rotation=0)));
      IO.TemperatureInput T_ref "Reference temperature for the process outlet" 
        annotation (Placement(transformation(
              origin={0,-30},
              extent={{-10,10},{10,-10}},
              rotation=90)));
    end Abstract;

    model Simple "A simple cooler implementation"
      extends Flow.UnitOperations.Coolers.Abstract;

      import Flow.Measurements.FlowTemperature;

      protected
      Flow.Sink sink "Makes up for lost heat" 
                                             annotation (Placement(
              transformation(extent={{0,80},{20,100}}, rotation=0)));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
                  -100},{100,100}}),
                          graphics),
                           Icon(coordinateSystem(preserveAspectRatio=false,
                extent={{-100,-100},{100,100}}),
                                graphics),
        Documentation(info="<html>
<p>This is a straightforward implementation of the abstract cooler that sets the
cooler's outlet temperature according to a given reference. There is a build-in,
customisable lag after which the outlet temperature will reach the reference value.</p>
</html>"));
      FlowTemperature T_out "Outlet temperature measurement" 
        annotation (Placement(transformation(extent={{40,60},{60,80}}, rotation=
                 0)));
      FlowTemperature T_in "Inlet temperature measurement" 
        annotation (Placement(transformation(extent={{-70,60},{-50,80}},
                rotation=0)));
      Modelica.Blocks.Continuous.FirstOrder lag(       initType=Modelica.Blocks.
            Types.Init.SteadyState, T=120)
          "A lag representing the inner control algorithm setting the outlet temperature"
        annotation (Placement(transformation(extent={{20,0},{40,20}}, rotation=
                  0)));
      Modelica.Blocks.Nonlinear.VariableLimiter limiter
          "Keeps requested temperature within reason" 
        annotation (Placement(transformation(extent={{-20,0},{0,20}}, rotation=
                  0)));
      outer Modelica.SIunits.Temperature T_env;
      constant Modelica.SIunits.Temperature eps = 0.01;

      public
      Modelica.SIunits.HeatFlowRate Q = sink.inlet.H "Cooling duty";
      parameter Boolean perfectControl = false "Neglect internal dynamics";

    equation
      limiter.limit2 = T_env;

      // No material loss
      sink.inlet.n = 0*sink.inlet.n;

      // Set the two variables to be equal
      // NOTE must do it here, these are both output variables.
      T_process_out = if perfectControl then T_ref else lag.y;
      connect(T_in.inlet, inlet) annotation (Line(points={{-68,70},{-80,70},{
                -80,-2.22045e-16},{-94,-2.22045e-16}}, color={0,127,127}));
      connect(T_in.T, T_process_in) annotation (Line(points={{-60,62},{-60,-14},
                {-94,-14}}, color={0,0,255}));
      connect(T_out.T, T_process_out) annotation (Line(points={{50,62},{50,-14},
                {94,-14}}, color={0,0,255}));
      connect(T_out.outlet, outlet) annotation (Line(points={{58,70},{76,70},{
                76,-2.22045e-16},{94,-2.22045e-16}}, color={0,127,127}));
      connect(sink.inlet, T_in.outlet) 
        annotation (Line(points={{1,90},{-20,90},{-20,70},{-52,70}}, color={0,
                127,127}));
      connect(T_out.inlet, T_in.outlet) 
        annotation (Line(points={{42,70},{-52,70}}, color={0,127,127}));
      connect(limiter.limit1, T_in.T) annotation (Line(points={{-22,18},{-40,18},
                {-40,40},{-60,40},{-60,62}}, color={0,0,127}));
      connect(limiter.u, T_ref) annotation (Line(points={{-22,10},{-56,10},{-56,
                -16},{0,-16},{0,-30},{5.55112e-16,-30}},  color={0,0,127}));
      connect(lag.u, limiter.y) 
        annotation (Line(points={{18,10},{1,10}}, color={0,0,127}));

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
        annotation (Placement(transformation(extent={{-64,28},{16,108}},
                rotation=0)));
      Measurements.GasFlowController mfc "Mass flow controller for cooling air"
        annotation (Placement(transformation(
              origin={50,56},
              extent={{-10,10},{10,-10}},
              rotation=270)));
      protected
      Flow.Sources.Environment env "Environmental air source" 
        annotation (Placement(transformation(extent={{100,46},{80,66}},
                rotation=0)));
      Flow.Sink airSink "Air outlet sink" 
                                         annotation (Placement(transformation(
                extent={{50,70},{70,90}}, rotation=0)));
      public
      Control.CoolerControl K annotation (Placement(transformation(extent={{12,
                  4},{28,20}}, rotation=0)));
    equation
      connect(exchanger.hot_2, outlet) annotation (Line(points={{-60,56},{-60,
                -2.22045e-16},{94,-2.22045e-16}}, color={0,127,127}));
      connect(inlet, exchanger.hot_1) annotation (Line(points={{-94,
                -2.22045e-16},{-80,-2.22045e-16},{-80,80},{-60,80}}, color={0,
                127,127}));
      connect(airSink.inlet, exchanger.cold_1) annotation (Line(points={{51,80},
                {12,80}}, color={0,127,127}));
      connect(env.outlet, mfc.inlet) annotation (Line(points={{81,56},{50,56}},
              color={0,127,127}));
      connect(mfc.outlet, exchanger.cold_2) annotation (Line(points={{40,56},{
                12,56}}, color={0,127,127}));
      connect(exchanger.T_hot_1, T_process_in) annotation (Line(points={{-60,88},
                {-74,88},{-74,-14},{-94,-14}}, color={0,0,255}));
      connect(exchanger.T_hot_2, T_process_out) annotation (Line(points={{-60,
                48},{-70,48},{-70,-14},{94,-14}}, color={0,0,255}));
      connect(K.V, mfc.V) annotation (Line(points={{29.6,12},{50,12},{50,46}},
              color={0,0,255}));
      connect(K.T_r, T_ref) annotation (Line(
            points={{10.4,12},{0,12},{0,-30},{-5.55112e-16,-30}},
            color={0,0,255},
            pattern=LinePattern.Dot));
      connect(exchanger.T_hot_2, K.T_m) annotation (Line(points={{-60,48},{-70,
                48},{-70,-14},{20,-14},{20,2.4}}, color={0,0,255}));
    end withAbstractExchanger;

    model LMTD "A cooler implemented with a LMTD heat exchanger"
      extends Flow.UnitOperations.Coolers.withAbstractExchanger(redeclare
            HeatExchangers.LMTD exchanger);
      annotation (Diagram(graphics),
                           Documentation(info="<html>
<p>A cooler using a LMTD implementation. It can for instance be applied as the
anode-loop cooler or as a recuperating heat exchanger on the anode side.</p>
</html>"),
        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                  100,100}}),
             graphics));
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
          Sink sink         annotation (Placement(transformation(extent={{60,-4},
                    {68,4}}, rotation=0)));
          annotation (Diagram(coordinateSystem(preserveAspectRatio=false,
                  extent={{-100,-100},{100,100}}),
                              graphics),
            experiment(StopTime=80),
            experimentSetupOutput);
          Sources.Solution sol(   T=330) 
                                      annotation (Placement(transformation(
                  extent={{-100,-20},{-80,0}}, rotation=0)));
        public
          replaceable Abstract cooler       annotation (Placement(
                transformation(extent={{-20,-20},{20,20}}, rotation=0)));
          Measurements.LiquidPump pump 
                    annotation (Placement(transformation(extent={{-60,-20},{-40,
                    0}}, rotation=0)));
          parameter VolumeFlowRate solution = 10E-6/60; // 10 ml/min
          parameter Temperature target = 315;

        equation
          pump.V = solution;
          cooler.T_ref = target;

          connect(cooler.outlet, sink.inlet) annotation (Line(points={{18.8,
                  1.06581e-15},{39.4,1.06581e-15},{39.4,3.88578e-17},{60.4,
                  3.88578e-17}}, color={0,127,127}));
          connect(pump.outlet, cooler.inlet) annotation (Line(points={{-50,
                  5.55112e-16},{-30,5.55112e-16},{-30,1.06581e-15},{-18.8,
                  1.06581e-15}}, color={0,127,127}));
          connect(sol.outlet, pump.inlet) annotation (Line(points={{-90,-10},{
                  -50,-10}}, color={0,127,127}));
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
      import Modelica.Constants.R;
      import Modelica.Electrical.Analog.Interfaces.PositivePin;
      import Modelica.Electrical.Analog.Interfaces.NegativePin;

      import Thermo.mw;
      import Thermo.rho;
      import Thermo.Species;
      import Thermo.speciesName;
      import Thermo.Condensables;
      import Thermo.Incondensables;

      import Units.MassTransportCoefficient;
      import Units.MolarFlow;
      import Units.MolarFlux;
      import Units.F;

      annotation (defaultComponentName="cell", Icon(coordinateSystem(
              preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
              graphics={Rectangle(
                extent={{-100,60},{100,0}},
                lineColor={0,0,0},
                fillColor={255,170,85},
                fillPattern=FillPattern.Solid), Rectangle(
                extent={{-100,2},{100,-60}},
                lineColor={0,0,0},
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid)}),
                                            Diagram(coordinateSystem(
              preserveAspectRatio=true,  extent={{-100,-100},{100,100}}),
                                                    graphics),
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
 
<p>The methanol diffusion through membrane and diffusion layer are taken from our results for
Johnson-Matthey MEAs.</p>
 
<p>The crossover current is calculated as proportional to the catalyst-layer concentration of
methanol on the anode. In turn, the difference between bulk and catalyst-layer concentration of 
methanol is proportional to the sum of crossover and reaction current densities.</p>
 
<p>The class calculates some quantities of interest, such as the anodic methanol 
concentration, and the cathodic partial pressures of oxygen and water. Note that all these are
based on the <em>exiting</em> flow.</p>
 
<p>Parameters have been taken from our publication about parameter regression, unless differently stated.</p>

<p>Stack heat capacity is given as 25 J/K, that corresponds to the heat capacity of two graphite monopolar
plates as used in our laboratory, weighing 53 g together, with a heat capacity of 466.25 J/kgK.
Remember to adjust this parameter if you change the number of cells, or temperature transients will become
very fast for no particular reason.</p>


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
        annotation (Placement(transformation(extent={{-110,20},{-90,40}},
                rotation=0)));
      Flow.FlowPort cathode_outlet "The cathode flow's outlet" 
        annotation (Placement(transformation(extent={{90,20},{110,40}},
                rotation=0)));
      Flow.FlowPort anode_inlet "The anode flow's inlet" 
        annotation (Placement(transformation(extent={{-110,-40},{-90,-20}},
                rotation=0)));
      Flow.FlowPort anode_outlet "The anode flow's outlet" 
        annotation (Placement(transformation(extent={{90,-40},{110,-20}},
                rotation=0)));
      PositivePin plus "Pole connected to the cathode" 
                                        annotation (Placement(transformation(
                extent={{-70,50},{-50,70}}, rotation=0)));
      NegativePin minus "Pole connected to the anode" 
                                      annotation (Placement(transformation(
                extent={{50,50},{70,70}}, rotation=0)));
      protected
      Flow.Measurements.FlowTemperature cathodeT
          "Cathode temperature measurement" 
        annotation (Placement(transformation(extent={{60,20},{80,40}}, rotation=
                 0)));
      Flow.Measurements.FlowConcentration anodeTC
          "Anode temperature and concentration measurement" 
                                      annotation (Placement(transformation(
                extent={{60,-40},{80,-20}}, rotation=0)));
      Flow.Sink nexus "Connection of all flows" 
                        annotation (Placement(transformation(extent={{-32,-10},
                  {-12,10}}, rotation=0)));

      public
      outer Pressure p_env "Environment pressure";
      outer Temperature T_env "Enviroment temperature";

      parameter Integer cells = 1 "Number of cells";
      parameter Area A = 26E-4 "Membrane active area";
      parameter HeatCapacity Cp = 25 "Overall heat capacity of the stack";
      parameter MassTransportCoefficient k_x_333 = 2.27E-6
          "Mass transport coefficient across the membrane at 333 K";
      parameter MassTransportCoefficient k_m_333 = 5.86E-6
          "Mass transport coefficient across the bulk at 333 K";
      parameter Boolean enableSanityChecks = true
          "Whether to activate checks for some non-negative quantities";
      parameter Real k_d_303 = 4.2 "Drag factor at 303 K";

      // Parameters for N115 membrane.
      Real k_d = k_d_303 + (T-303.15)/40 "Drag factor for N115";
      MassTransportCoefficient k_m = k_m_333*exp(1395*(1/333-1/T))
          "Mass transport coefficient";
      MassTransportCoefficient k_x = k_x_333*exp(1395*(1/333-1/T))
          "Cross-over Mass transport coefficient";

      Real a = k_m*b "Partial derivative of crossover flux wrt. concentration";
      Real aAn = a*A*cells
          "Partial derivative of crossover flow wrt. concentration";
      Real b = 1/(1+k_m/k_x)
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
                                                             annotation (Placement(
              transformation(extent={{100,-8},{120,12}}, rotation=0)));

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

      constant Real eps = 1E-10 "Maximum numerical noise for molar flows";

    equation
      // Anode-side mass balance, accounting for reaction, drag and crossover
      anode_inlet.n + anode_outlet.n + anode_nu*n_H + anode_xi*n_x = 0*anode_inlet.n;

      // Cathode-side mass balance, accounting for reaction, drag and crossover
      cathode_inlet.n + cathode_outlet.n + cathode_nu*n_H + cathode_xi*n_x = 0*cathode_inlet.n;

      // The energy "lost" from the heat balance is the electrical power.
      nexus.inlet.H = I*V + der(T)*Cp;

      // Definition of oxygen partial pressure. On the denominator, the sum of vapours (methanol and water) and gases (all others).
      p_o2 = p_env * cathodeT.inlet.n[Species.Oxygen] / (sum(cathodeT.vapour) + sum(cathodeT.inlet.n[Incondensables]));

      // Definition of water partial pressure (on the cathode side). On the denominator, the sum of vapours (methanol and water) and gases (all others).
      p_h2o = p_env * cathodeT.vapour[Species.Water]  / (sum(cathodeT.vapour) + sum(cathodeT.inlet.n[Incondensables]));

      // Methanol transport: binds c, c_cl and i (N_x is a function of c_ac, N_H of i).
      k_m * (c-c_cl) = N_x + N_H/6;

      // Crossover methanol flux.
      N_x = k_x * c_cl;

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

        for i in Species loop
          assert( cathode_outlet.n[i] < eps, "==> "+speciesName(i)+" is entering from the cathode outlet.");
          assert( anode_outlet.n[i] < eps, "==> "+speciesName(i)+" is entering from the anode outlet.");
          assert( cathode_inlet.n[i] > -eps, "==> "+speciesName(i)+" is exiting from the cathode inlet.");
          assert( anode_inlet.n[i] > -eps, "==> "+speciesName(i)+" is exiting from the anode inlet.");
        end for;
      end if;

      connect(cathodeT.outlet, cathode_outlet) 
        annotation (Line(points={{78,30},{100,30}}, color={0,127,127}));
      connect(cathode_inlet, nexus.inlet) annotation (Line(points={{-100,30},{
                -46,30},{-46,0},{-31,0},{-31,4.44089e-16}}, color={0,127,127}));
      connect(cathodeT.inlet, nexus.inlet)       annotation (Line(points={{62,30},
                {-40,30},{-40,4.44089e-16},{-31,4.44089e-16}},     color={0,127,
                127}));
      connect(anodeTC.outlet, anode_outlet)       annotation (Line(points={{78,
                -30},{100,-30}}, color={0,127,127}));
      connect(anodeTC.inlet, nexus.inlet)          annotation (Line(points={{62,-30},
                {-40,-30},{-40,4.44089e-16},{-31,4.44089e-16}},      color={0,
                127,127}));
      connect(T, cathodeT.T) annotation (Line(points={{110,2},{70,2},{70,22}},
              color={0,0,255}));
      connect(anode_inlet, nexus.inlet) annotation (Line(points={{-100,-30},{
                -46,-30},{-46,4.44089e-16},{-31,4.44089e-16}}, color={0,127,127}));
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
      parameter Resistance R = 0.05 "Internal resistance";

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
                            annotation (Placement(transformation(extent={{6,0},
                    {42,34}}, rotation=0)));
          Sources.Solution methanolSolution(   T=anodeInletTemperature) 
                                            annotation (Placement(
                transformation(extent={{-66,-30},{-54,-18}}, rotation=0)));
          annotation (Diagram(coordinateSystem(preserveAspectRatio=false,
                  extent={{-100,-100},{100,100}}),
                              graphics));
          Measurements.LiquidPump pump(T(start=anodeInletTemperature))
            "Pump for the anode flow"         annotation (Placement(
                transformation(extent={{-42,-30},{-30,-18}}, rotation=0)));
          Sources.Environment air 
                              annotation (Placement(transformation(extent={{-70,
                    28},{-50,48}}, rotation=0)));
          Measurements.GasFlowController blower 
                                   annotation (Placement(transformation(extent=
                    {{-40,34},{-32,42}}, rotation=0)));
          Sink anodeSink     annotation (Placement(transformation(extent={{62,
                    10},{68,16}}, rotation=0)));
          Sink cathodeSink     annotation (Placement(transformation(extent={{62,
                    18},{68,24}}, rotation=0)));
          ConstantCurrent I_cell(I=5) annotation (Placement(transformation(
                  extent={{12,48},{34,72}}, rotation=0)));
          Modelica.Electrical.Analog.Basic.Ground ground
            "Negative pole to zero voltage" annotation (Placement(
                transformation(extent={{38,40},{58,60}}, rotation=0)));
        equation
          pump.V = anodeFlow;
          blower.V = cathodeFlow;

          connect(methanolSolution.outlet, pump.inlet) 
                                                  annotation (Line(points={{-60,
                  -24},{-36,-24}}, color={0,127,127}));
          connect(blower.outlet, stack.cathode_inlet)    annotation (Line(
                points={{-36,42},{-18,42},{-18,22.1},{6,22.1}}, color={0,127,
                  127}));
          connect(air.outlet, blower.inlet) 
                                       annotation (Line(points={{-51,38},{-36,
                  38}}, color={0,127,127}));
          connect(cathodeSink.inlet, stack.cathode_outlet)       annotation (Line(
                points={{62.3,21},{52.15,21},{52.15,22.1},{42,22.1}}, color={0,
                  127,127}));
          connect(anodeSink.inlet, stack.anode_outlet)       annotation (Line(
                points={{62.3,13},{52.15,13},{52.15,11.9},{42,11.9}}, color={0,
                  127,127}));
          connect(ground.p, I_cell.n) 
            annotation (Line(points={{48,60},{34,60}}, color={0,0,255}));
          connect(I_cell.p, stack.plus)    annotation (Line(points={{12,60},{12,
                  27.2},{13.2,27.2}}, color={0,0,255}));
          connect(I_cell.n, stack.minus)    annotation (Line(points={{34,60},{
                  34.8,60},{34.8,27.2}}, color={0,0,255}));
          connect(pump.outlet, stack.anode_inlet)    annotation (Line(points={{-36,-18},
                  {-16,-18},{-16,11.9},{6,11.9}},          color={0,127,127}));
        end AbstractStackTest;

        model ConstantVoltageStackTest "Test for the constant-voltage model"
          extends AbstractStackTest(redeclare ConstantVoltage stack,
              anodeFlow = stack.cells*30E-6/60,
              cathodeFlow = stack.cells*30E-5/60);
          annotation (Diagram(coordinateSystem(preserveAspectRatio=false,
                  extent={{-100,-100},{100,100}}),
                              graphics));
        end ConstantVoltageStackTest;

        model TheveninStackTest "Test for the Thevenin-circuit model"
          extends AbstractStackTest(redeclare Thevenin stack(cells=3,V0=2.1,R=0.015),
              anodeFlow = stack.cells*30E-6/60,
              cathodeFlow = stack.cells*30E-5/60);
          annotation (experiment(StopTime=40, Algorithm="Dassl"));
        end TheveninStackTest;
      end Test;
    end Stack;

    package Test
      partial model AbstractSeparatorTest "Test case for the separator unit"

        replaceable AbstractSeparator separator 
                                    annotation (Placement(transformation(extent={{-22,
                  -12},{26,38}}, rotation=0)));
      protected
        Flow.Sink liquidSink 
                          annotation (Placement(transformation(extent={{76,-14},
                  {84,-6}}, rotation=0)));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
                  -100},{100,100}}),
                            graphics));
        Sink gasSink       annotation (Placement(transformation(extent={{76,26},
                  {84,34}}, rotation=0)));
        Sources.Environment env         annotation (Placement(transformation(
                extent={{-54,22},{-34,42}}, rotation=0)));
        Sources.Solution solution         annotation (Placement(transformation(
                extent={{-78,8},{-68,18}}, rotation=0)));
      public
        inner parameter Modelica.SIunits.Pressure p_env = 101325;
        inner parameter Units.Temperature T_env = 298.15;
        inner parameter Units.RelativeHumidity RH_env = 60;

        Measurements.FlowTemperature T_liquid 
          annotation (Placement(transformation(extent={{40,-20},{60,0}})));
        Measurements.FlowTemperature T_gas 
          annotation (Placement(transformation(extent={{40,20},{60,40}})));
      equation
        sum(env.outlet.n) = -1;
        sum(solution.outlet.n) = -2;

        connect(env.outlet, separator.inlet)        annotation (Line(points={{
                -35,32},{-22,32},{-22,13}}, color={0,127,127}));
        connect(solution.outlet, separator.inlet) 
          annotation (Line(points={{-73,13},{-22,13}}, color={0,127,127}));
        connect(separator.gasOutlet, T_gas.inlet) annotation (Line(
            points={{18.8,23},{20,23},{20,30},{42,30}},
            color={0,127,127},
            smooth=Smooth.None));
        connect(T_gas.outlet, gasSink.inlet) annotation (Line(
            points={{58,30},{76.4,30}},
            color={0,127,127},
            smooth=Smooth.None));
        connect(separator.liquidOutlet, T_liquid.inlet) annotation (Line(
            points={{18.8,3},{20,3},{20,-10},{42,-10}},
            color={0,127,127},
            smooth=Smooth.None));
        connect(T_liquid.outlet, liquidSink.inlet) annotation (Line(
            points={{58,-10},{76.4,-10}},
            color={0,127,127},
            smooth=Smooth.None));
      end AbstractSeparatorTest;

      model SeparatorTest
        extends AbstractSeparatorTest(redeclare Separator separator);
      end SeparatorTest;

      model CapillarySeparatorTest
        extends AbstractSeparatorTest(redeclare CapillarySeparator separator);

      equation
        separator.backPressure.p = separator.delta_pc - exp(-time*10) * 1000 * sin(time*100);
      end CapillarySeparatorTest;

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

        Burner burner annotation (Placement(transformation(extent={{-2,-20},{38,
                  20}}, rotation=0)));
        Sources.Solution methanolSolution annotation (Placement(transformation(
                extent={{-100,20},{-80,40}}, rotation=0)));
        Measurements.LiquidPump pump 
                  annotation (Placement(transformation(extent={{-60,40},{-40,20}},
                rotation=0)));
        annotation (Diagram(graphics),
                             Documentation(info="<html>
</html>"));
        Sources.Environment env 
                            annotation (Placement(transformation(extent={{-100,
                  -20},{-80,0}}, rotation=0)));
        Measurements.GasFlowController mfc 
                              annotation (Placement(transformation(extent={{-60,
                  -20},{-40,0}}, rotation=0)));
        Sink sinkPort     annotation (Placement(transformation(extent={{80,-10},
                  {100,10}}, rotation=0)));
        Measurements.FlowTemperature T_in 
                             annotation (Placement(transformation(extent={{-30,
                  -10},{-10,10}}, rotation=0)));
      equation
        pump.V = solution;
        mfc.V = air;
        connect(pump.inlet, methanolSolution.outlet) annotation (Line(points={{
                -50,30},{-90,30}}, color={0,127,127}));
        connect(mfc.inlet, env.outlet) annotation (Line(points={{-50,-10},{-81,
                -10}}, color={0,127,127}));
        connect(T_in.outlet, burner.inlet) annotation (Line(points={{-12,
                6.10623e-16},{-7,6.10623e-16},{-7,1.22125e-15},{-1.6,
                1.22125e-15}}, color={0,127,127}));
        connect(T_in.inlet, pump.outlet) annotation (Line(points={{-28,
                6.10623e-16},{-40,6.10623e-16},{-40,20},{-50,20}}, color={0,127,
                127}));
        connect(T_in.inlet, mfc.outlet) annotation (Line(points={{-28,
                6.10623e-16},{-40,6.10623e-16},{-40,5.55112e-16},{-50,
                5.55112e-16}}, color={0,127,127}));
        connect(burner.outlet, sinkPort.inlet) annotation (Line(points={{38.4,
                1.22125e-15},{60.2,1.22125e-15},{60.2,4.44089e-16},{81,
                4.44089e-16}}, color={0,127,127}));
      end BurnerTest;

      model MixerTest "Test for the mixer unit"
        inner parameter Modelica.SIunits.Pressure p_env = 101325;
        inner parameter Units.Temperature T_env = 298.15;
        inner parameter Units.RelativeHumidity RH_env = 60;

        replaceable Mixer mixer(
          c(fixed=true),
          T(fixed=true),
          V(fixed=true, start=500E-6)) 
                    annotation (Placement(transformation(extent={{-20,0},{0,20}},
                rotation=0)));
        Sources.Solution anodicLoop(   T=330)
          "Solution coming from the anodic loop" 
          annotation (Placement(transformation(extent={{-20,40},{0,60}},
                rotation=0)));
        Sources.Methanol fuelTank "Methanol from the fuel tank" 
          annotation (Placement(transformation(extent={{0,-40},{20,-20}},
                rotation=0)));
        Sources.Solution condenser(   C=0, T=310)
          "Water recovered from the cathode outlet" 
          annotation (Placement(transformation(extent={{20,0},{40,20}},
                rotation=0)));
        Measurements.FlowTemperature flowTemperature 
                                        annotation (Placement(transformation(
                extent={{-64,-26},{-44,-6}}, rotation=0)));
        Sink sinkPort     annotation (Placement(transformation(extent={{-28,-20},
                  {-20,-12}}, rotation=0)));
      equation

        annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                  -100},{100,100}}),
                            graphics));
        sum(fuelTank.outlet.n) = -0.1;
        sum(anodicLoop.outlet.n) = -1;
        sum(condenser.outlet.n) = -0.4;
        sum(mixer.outlet.n) = -1.5;

        connect(flowTemperature.outlet, sinkPort.inlet) 
          annotation (Line(points={{-46,-16},{-27.6,-16}}, color={0,127,127}));
        connect(mixer.outlet, flowTemperature.inlet) annotation (Line(points={{
                -18,10},{-72,10},{-72,-16},{-62,-16}}, color={0,127,127}));
        connect(condenser.outlet, mixer.waterInlet) 
          annotation (Line(points={{30,10},{-2,10}}, color={0,127,127}));
        connect(anodicLoop.outlet, mixer.loopInlet) 
          annotation (Line(points={{-10,50},{-10,18}}, color={0,127,127}));
        connect(mixer.fuelInlet, fuelTank.outlet) 
          annotation (Line(points={{-10,2},{-10,-30},{10,-30}}, color={0,127,
                127}));
      end MixerTest;

      model HydrostaticMixerTest
        extends MixerTest(redeclare HydrostaticMixer mixer);
      end HydrostaticMixerTest;

      model ElasticMixerTest
        extends MixerTest(redeclare ElasticMixer mixer(V0=5E-4));
      end ElasticMixerTest;
    end Test;
  end UnitOperations;

  package IntegratedOperations
    partial model Equilibrium "Modelling of intensive properties"

      import Units.MolarEnthalpy;
      import Modelica.SIunits.MoleFraction;
      import Thermo.Species;
      import Thermo.Incondensables;
      import Thermo.Condensables;
      import Thermo.Phases;
      import Thermo.h;
      import Thermo.K;

      Units.Temperature T "Unit temperature";

      MolarEnthalpy h_tot = beta*h_g + (1-beta)*h_l "Average molar enthalpy";
      MolarEnthalpy h_g = sum(h(T, i, Phases.Gas)*y[i] for i in Species)
        "Molar enthalpy in gas";
      MolarEnthalpy h_l = sum(h(T, i, Phases.Liquid)*x[i] for i in Condensables)
        "Molar enthalpy in liquid";

      MoleFraction beta "Fraction of moles in gas phase";

      MoleFraction x[Species] "Liquid molar fraction";
      MoleFraction y[Species] "Gaseous molar fraction";
      MoleFraction z[Species] "Overall molar fraction";

    equation
      x[Condensables]   = {z[i]/(1+beta*(K(T, i)-1)) for i in Condensables};
      x[Incondensables] = 0*x[Incondensables];

      y[Condensables]   = {K(T,i)*x[i] for i in Condensables};
      y[Incondensables] = if noEvent(beta > 0) then z[Incondensables]/beta else 0*y[Incondensables];

      if K(T, Species.Water) >= 1 or Thermo.rr(z[Species.Methanol], z[Species.Water], T) >= 1 then
        beta = 1;
      else
        beta = Thermo.rr(z[Species.Methanol], z[Species.Water], T);
      end if;

      annotation (Documentation(info="<html>
<p>This class models the equilibrium relationships between liquid and gas
compositions in a two-phase mixture.</p>
<p>Equilibrium is calculated with the <tt>Thermo</tt> library's analytical
Rachford-Rice solution, which gives a given &beta; value once the 
<tt><b>z</b></tt> values are known.</p>
<p>The model has 6 degrees of freedom, corresponding to the five compositions 
and temperature.</p>
</html>"));
    end Equilibrium;

    partial model Balances "Modelling of mass and energy balances"

      import Modelica.SIunits.AmountOfSubstance;
      import Modelica.Constants.eps;
      import Thermo.Species;
      import Thermo.Incondensables;
      import Units.MolarFlow;

      Modelica.SIunits.InternalEnergy U "Energy holdup";
      AmountOfSubstance n[Species](each start=sqrt(eps)) "Molar holdup";
      AmountOfSubstance n_tot = sum(n) "Total number of moles";
      MolarFlow F_env(start=0) = sum(envPort.n)
        "Mole exchange through the environment port";
      MolarFlow L(start=0) "Liquid exchanged with the environment";

      FlowPort outlet "The mixer's outlet" 
                            annotation (Placement(transformation(extent={{-90,
                -10},{-70,10}}, rotation=0)));
      FlowPort fuelInlet "The methanol-feed inlet" 
                             annotation (Placement(transformation(extent={{-10,
                -90},{10,-70}}, rotation=0)));
      FlowPort inlet "The anodic-loop inlet" 
                             annotation (Placement(transformation(extent={{70,
                -10},{90,10}}, rotation=0)));
      FlowPort envPort "Connection with the environment" 
                            annotation (Placement(transformation(extent={{-10,
                70},{10,90}}, rotation=0)));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}),
                          graphics),
                           Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={Text(
              extent={{-100,160},{100,100}},
              lineColor={0,0,0},
              fillColor={255,85,85},
              fillPattern=FillPattern.Solid,
              textString="%name")}),
        Documentation(info="<html>
<p>This class models extensive balances of a mixer. Given four ports,
it calculates and keeps track of the accumulation of substance and
internal energy.</p>
<p>It also defines a few handy quantities, such as molar exchange through
the environment port, <tt>F_env</tt>, molar exchange with the environment 
<tt>L</tt> (which cannot be defined at this point, since equilibrium 
calculations have not yet been introduced), and the total molar inflow of
dry gases.</p>
</html>"));
    protected
      Sink accumulator "Accumulation sink" 
                                   annotation (Placement(transformation(extent=
                {{0,-10},{20,10}}, rotation=0)));
    protected
      MolarFlow dryIn = sum(inlet.n[Incondensables]) "Inflow of dry gases";

    equation
      der(U) = accumulator.inlet.H;
      der(n) = accumulator.inlet.n;

      connect(accumulator.inlet, fuelInlet) 
                                    annotation (Line(points={{1,4.44089e-16},{0,
              4.44089e-16},{0,-80},{5.55112e-16,-80}}, color={0,127,127}));
      connect(accumulator.inlet, outlet) 
                                 annotation (Line(points={{1,4.44089e-16},{0,
              4.44089e-16},{0,5.55112e-16},{-80,5.55112e-16}}, color={0,127,127}));
      connect(accumulator.inlet, inlet) 
                                annotation (Line(points={{1,4.44089e-16},{0,
              4.44089e-16},{0,0},{80,0},{80,5.55112e-16}}, color={0,127,127}));
      connect(accumulator.inlet, envPort) 
                                  annotation (Line(points={{1,4.44089e-16},{0,
              4.44089e-16},{0,80},{5.55112e-16,80}}, color={0,127,127}));
    end Balances;

    partial model VolumeSum "Modelling of volume balance"
      import Modelica.SIunits.Volume;

      parameter Volume V = 5E-6 "Total physical volume";

      Volume V_g "Gas-phase volume";
      Volume V_l "Liquid-phase volume";

    equation
      V = V_g + V_l;

      annotation (Documentation(info="<html>
<p>This simple class includes the variables used to keep track
of the liquid and gaseous volumes, simply by declaring them and
their relationship, i.e. that their sum is equal to a given
parameter.</p>
</html>"));
    end VolumeSum;

    partial model EquilibriumAndBalances
      "Joining equilibrium and mass and energy balances"
      extends Equilibrium;
      extends Balances;

      import Units.MolarFlow;
      import Modelica.Constants.eps;
      import Modelica.SIunits.AmountOfSubstance;
      import Modelica.SIunits.EnthalpyFlowRate;
      import Thermo.Species;
      import Thermo.Incondensables;
      import Thermo.Condensables;
      import Thermo.Phases;

      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}),
                          graphics),
                           Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={Text(
              extent={{-100,160},{100,100}},
              lineColor={0,0,0},
              fillColor={255,85,85},
              fillPattern=FillPattern.Solid,
              textString="%name")}),
        Documentation(info="<html>
<p>This partial class is a first step towards the completed mixer. It does
not yet include the volume balance, but joins the mass and energy balances
with the vapour-liquid equilibrium relations.</p>

<p>Joining the mass balances with the equilibrium is straightforward, as
<tt>z</tt> is trivially related to the molar holdup <tt>n</tt>.

<p>The class defines a new boolean variable, <tt>overflow</tt>, which indicates when
the sole gas flow is insufficient to satisfy the molar balance, and liquid
must therefore be expelled from the mixer as well. The rules to define when
<tt>overflow</tt> is true or false are left to be defined in a child class.</p>

<p>Depending on the state of overflow, the composition of flow through the
enviroment port is set. With no overflow, the composition is either that of
air (when entering) or that of the gas phase (when exiting). In case of overflow,
instead, the composition is a bit trickier: there is one gaseous component,
consisting of all the incondensable gases entering the system along with
associated vapours from the liquid phase, and then the liquid phase itself.
There is no case for entering composition when in overflow, since that would
immediately mean that <tt>overflow</tt> must be false, falling back to the
two-phase model.</p>

<p>The class also defines the exiting composition, usually the average 
composition in the mixer (<tt>z</tt>), but in case of overflow the liquid
composition: the reason for this is that the model will be numerically more
stable, with a negligibly small error being introduced.</p>
</html>"));

      parameter Boolean mixedOutlet = true
        "Whether to use z instead of x in outlet";

      Boolean overflow "Whether tank is full and overflowing";

      AmountOfSubstance n_l_tot = (1-beta)*n_tot "Total moles in liquid phase";
      AmountOfSubstance[:] n_l = x*n_l_tot "Moles in liquid phase";

      AmountOfSubstance n_g_tot = beta*n_tot "Total moles in gas phase";
      AmountOfSubstance[:] n_g = y*n_g_tot "Moles in gas phase";

    protected
      Thermo.Air air "Environmental conditions";
      MolarFlow wetOut = - dryIn/(1-sum(y[Condensables]))
        "Entering dry gases, exiting after humidification";
      EnthalpyFlowRate H_wet "Enthalpy flow associated to wetOut";

    equation
      z = n / sum(n);

      H_wet = - sum(Thermo.h(T, i, Phases.Gas)*inlet.n[i] for i in Incondensables)
              + sum(Thermo.h(T, c, Phases.Gas)*y[c]*wetOut for c in Condensables);

      if mixedOutlet then
        outlet.n[2:end] = sum(outlet.n) * (if overflow then x[2:end] else z[2:end]);
        outlet.H        = sum(outlet.n) * (if overflow then h_l else h_tot);
      else
        outlet.n[2:end] = sum(outlet.n) * x[2:end];
        outlet.H        = sum(outlet.n) * h_l;
      end if;

      if overflow then
        envPort.n[Incondensables] = -inlet.n[Incondensables];
        envPort.n[Condensables]   = wetOut*y[Condensables] + L*x[Condensables];
        envPort.H                 = H_wet + L*h_l;
      else
        // NOTE: L*x is supposed to be 0 here, but it helps with the initialisation.
        envPort.n[Condensables]   = semiLinear(F_env, air.y[Condensables],   y[Condensables]) + L*x[Condensables];
        envPort.n[Incondensables] = semiLinear(F_env, air.y[Incondensables], y[Incondensables]);
        envPort.H                 = semiLinear(F_env, air.H,                 h_g);
      end if;

    end EquilibriumAndBalances;

    model Mixer "A unit mixing four molar flows."
      extends EquilibriumAndBalances;
      extends VolumeSum;

      import Modelica.Constants.eps;
      import Modelica.SIunits.Concentration;
      import Thermo.Species;
      import Thermo.Condensables;
      import Thermo.Phases;
      import Thermo.rho;
      import Thermo.mw;

      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}),
                          graphics),
                           Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={
            Ellipse(
              extent={{-80,80},{80,-80}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.CrossDiag),
            Rectangle(
              extent={{-80,0},{80,80}},
              lineColor={0,0,255},
              pattern=LinePattern.None,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Rectangle(
              extent={{-80,0},{80,80}},
              lineColor={0,0,0},
              pattern=LinePattern.None,
              fillColor={255,255,255},
              fillPattern=FillPattern.CrossDiag),
            Line(points={{-80,0},{-80,80},{80,80},{80,0}}, color={0,0,0}),
            Line(
              points={{0,6},{0,-54}},
              color={0,0,0},
              thickness=0.5),
            Line(
              points={{0,6},{-52,36}},
              color={0,0,0},
              thickness=0.5),
            Line(
              points={{0,6},{52,36}},
              color={0,0,0},
              thickness=0.5),
            Text(
              extent={{-100,160},{100,100}},
              lineColor={0,0,0},
              fillColor={255,85,85},
              fillPattern=FillPattern.Solid,
              textString="%name")}),
        Documentation(info="<html>
<p>The model of an integrated separator and mixer, in which gas is preferentially
dumped to the environment, but liquid can be removed as well when all gas has
been exhausted.</p>

<p>The class integrates the volume balance with the other balances (mass and energy)
and the equilibrium relations through the definitions of volumes <tt>V_g</tt> and 
<tt>V_l</tt>, which depend on gas and liquid composition on one side, and on temperature
on the other.</p>

<p>The class also defines the relation between internal energy and temperature through
the overall specific enthalpy.</p>

<p>The most important relation is the definition of the <tt>overflow</tt> variable: it
is defined to be true when the gaseous volume is very small (not exactly zero, as this
would cause divide-by-zero errors e.g. in the calculations of gas molar fractions), 
<em>and</em> liquid is flowing out through the environment port. It is not enough that
there is no more gas left in the system to state that we are in an overflom condition:
if we do not specify that liquid can go only out, we will not be able to return back to 
the non-overflow state, because no gas will remain in the mixer in the overflow 
condition by construction (see how this was implemented in model <tt>EquilibriumAndBalances</tt>).
As a result, we would see liquid at the same condition as in the mixer <em>entering</em>
the mixer from the environment, which is a clear nonsense.</p>

<p>The model also defines a convenience variable for methanol concentration.</p>

<p>The model is initialised so than the initial gas composition is equal to that in the
atmosphere (in its dry part), and that the water moles are about half of the mixer
volume. Unfortunately, direct initialisation by setting a liquid or gaseous volume is
numerically unstable. Finally, the concentration start value is strictly enforced 
(<tt>fixed</tt>).</p>

</html>"));

      outer Modelica.SIunits.Pressure p_env "Environment pressure";
      outer Modelica.SIunits.Temperature T_env "Environment temperature";

      Concentration c(start=1000,fixed=true) "Methanol concentration";

    equation
      c = n_l[Species.Methanol]/V_l;

      V_g = sum(n_g[i]*mw(i)/rho(T, i, Phases.Gas) for i in Species);
      V_l = sum(n_l[i]*mw(i)/rho(T, i, Phases.Liquid) for i in Condensables);

      U + p_env*V = h_tot*sum(n);

      overflow = V_g < V/1000 and L < 0;

    initial equation
      y[Species.Oxygen] / 0.21 = y[Species.Nitrogen] / 0.79;
      y[Species.CarbonDioxide] = 385E-6;

      n[Species.Water] = 0.5*V*rho(T, Species.Water, Phases.Liquid)/mw(Species.Water);

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

      protected
        constant Modelica.SIunits.Time second = 1 "To get adimensional time";
      equation
        eq.z[1:2] = {z_methanol, z_water};
        der(eq.z[3:4]) = {0,0};
        sum(eq.z) = 1;

        eq.T = T_start + time/second*(T_stop-T_start);

      end EquilibriumTest;

      model BalancesTest

        model MyBalances
          extends Balances;
        end MyBalances;

        MyBalances bal annotation (Placement(transformation(extent={{-20,-20},{
                  20,20}}, rotation=0)));
        Sources.Solution solution annotation (Placement(transformation(extent={
                  {60,-10},{80,10}}, rotation=0)));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent=
                  {{-100,-100},{100,100}}),
                            graphics));
      equation

        bal.L = 0;

        sum(solution.outlet.n) = -1;

        connect(solution.outlet, bal.inlet) annotation (Line(points={{70,
                6.66134e-16},{44,6.66134e-16},{44,1.22125e-15},{16,1.22125e-15}},
              color={0,127,127}));

      end BalancesTest;

      model VolumeSumTest
        import Modelica.SIunits.MoleFraction;
        import Modelica.SIunits.Temperature;

        model MyVolumeBalance
          extends VolumeSum;
        end MyVolumeBalance;

        MyVolumeBalance eq;

      protected
        constant Modelica.SIunits.Time second = 1 "To get adimensional time";
      equation
        eq.V_g = 0 + time/second*eq.V;

      end VolumeSumTest;

      model EquilibriumAndBalancesTest "Test for the mixer unit"

        import Units.MolarFlow;
        import Modelica.SIunits.VolumeFlowRate;
        import Thermo.Molecules.All;

        inner parameter Modelica.SIunits.Temperature T_env = 298.15;
        inner parameter Units.RelativeHumidity RH_env = 60;
        inner parameter Modelica.SIunits.Pressure p_env = 101325;

        parameter VolumeFlowRate airFlow = 1E-6;
        parameter VolumeFlowRate solutionOut = 2E-6;
        parameter VolumeFlowRate solutionIn = 1E-6;
        parameter VolumeFlowRate fuel = 0.1E-6;

        model MyEquilibriumAndBalances
          extends EquilibriumAndBalances;
        end MyEquilibriumAndBalances;

        Sources.Solution anodicLoop(C=700, T=330)
          "Solution coming from the anodic loop" 
          annotation (Placement(transformation(extent={{40,24},{50,34}},
                rotation=0)));
        Sources.Methanol fuelTank "Methanol from the fuel tank" 
          annotation (Placement(transformation(extent={{4,-36},{16,-24}},
                rotation=0)));
        Sources.Environment env annotation (Placement(transformation(extent={{
                  80,-20},{60,0}}, rotation=0)));
        Measurements.GasFlowController mfc 
          annotation (Placement(transformation(extent={{22,-16},{34,-4}},
                rotation=0)));
        annotation (Diagram(graphics),
                             experiment);
        Measurements.LiquidPump pump_in 
                               annotation (Placement(transformation(
              origin={29,29},
              extent={{-5,-5},{5,5}},
              rotation=180)));
        Measurements.LiquidPump fuel_pump 
                               annotation (Placement(transformation(extent={{
                  -16,-36},{-4,-24}}, rotation=0)));
        MyEquilibriumAndBalances mixer annotation (Placement(transformation(
                extent={{-20,0},{0,20}}, rotation=0)));
        Sink sinkPort     annotation (Placement(transformation(
              origin={-76,16},
              extent={{-4,-4},{4,4}},
              rotation=180)));
        Sink sink annotation (Placement(transformation(
              origin={-10,40},
              extent={{-4,-4},{4,4}},
              rotation=90)));
        Measurements.PeristalticPump pump_out 
          annotation (Placement(transformation(extent={{-52,4},{-40,16}},
                rotation=0)));
      equation

        mixer.overflow = false;
        mixer.F_env = -1E-5;
        mixer.U = mixer.h_tot * sum(mixer.n);

        mfc.V = airFlow;
        pump_in.V = solutionIn;
        pump_out.V = solutionOut;
        fuel_pump.V = fuel;

        connect(env.outlet, mfc.inlet) annotation (Line(points={{61,-10},{28,
                -10}}, color={0,127,127}));
        connect(pump_in.inlet, anodicLoop.outlet) annotation (Line(points={{29,
                29},{45,29}}, color={0,127,127}));
        connect(fuelTank.outlet, fuel_pump.inlet) annotation (Line(points={{10,
                -30},{-10,-30}}, color={0,127,127}));
        connect(mixer.inlet, pump_in.outlet) annotation (Line(points={{-2,10},{
                14,10},{14,24},{29,24}}, color={0,127,127}));
        connect(mfc.outlet, mixer.inlet) annotation (Line(points={{28,-4},{14,
                -4},{14,10},{-2,10}}, color={0,127,127}));
        connect(fuel_pump.outlet, mixer.fuelInlet) annotation (Line(points={{
                -10,-24},{-10,2}}, color={0,127,127}));
        connect(sink.inlet, mixer.envPort) annotation (Line(points={{-10,36.4},
                {-10,18}}, color={0,127,127}));

      initial equation
        mixer.T = T_env;

      equation
        connect(pump_out.inlet, mixer.outlet) annotation (Line(points={{-46,10},
                {-18,10}}, color={0,127,127}));
        connect(pump_out.outlet, sinkPort.inlet) annotation (Line(points={{-46,
                16},{-72.4,16}}, color={0,127,127}));
      end EquilibriumAndBalancesTest;

      model MixerTest "Test for the mixer unit"

        import Units.MolarFlow;

        inner parameter Modelica.SIunits.Temperature T_env = 298.15;
        inner parameter Units.RelativeHumidity RH_env = 60;
        inner parameter Modelica.SIunits.Pressure p_env = 101325;

        parameter MolarFlow airFlow = 1E-6;
        parameter MolarFlow solutionOut = 1E-3;
        parameter MolarFlow fuel = 0.1E-3;

        Mixer mixer annotation (Placement(transformation(extent={{-20,0},{0,20}},
                rotation=0)));
        Sources.Solution anodicLoop(C=700, T=330)
          "Solution coming from the anodic loop" 
          annotation (Placement(transformation(extent={{40,24},{50,34}},
                rotation=0)));
        Sources.Methanol fuelTank "Methanol from the fuel tank" 
          annotation (Placement(transformation(extent={{4,-36},{16,-24}},
                rotation=0)));
        Sink sinkPort     annotation (Placement(transformation(
              origin={-60,16},
              extent={{-4,-4},{4,4}},
              rotation=180)));
        Sources.Environment env annotation (Placement(transformation(extent={{
                  80,-20},{60,0}}, rotation=0)));
        Measurements.GasFlowController mfc 
          annotation (Placement(transformation(extent={{22,-16},{34,-4}},
                rotation=0)));
        Measurements.PeristalticPump pump_out 
                               annotation (Placement(transformation(extent={{
                  -40,4},{-28,16}}, rotation=0)));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent=
                  {{-100,-100},{100,100}}),
                            graphics),
                             experiment(StopTime=500));
        Measurements.LiquidPump pump_in 
                               annotation (Placement(transformation(
              origin={29,29},
              extent={{-5,-5},{5,5}},
              rotation=180)));
        Measurements.LiquidPump fuel_pump 
                               annotation (Placement(transformation(extent={{
                  -16,-36},{-4,-24}}, rotation=0)));
        Sink Overflow     annotation (Placement(transformation(
              origin={-10,32},
              extent={{-4,-4},{4,4}},
              rotation=90)));
        Modelica.Blocks.Sources.Sine sine(
          freqHz=0.01,
          amplitude=0.005,
          offset=1E-3) annotation (Placement(transformation(extent={{-16,54},{
                  -4,66}}, rotation=0)));
      equation
        mfc.F = airFlow;
        pump_in.F = sine.y;
        pump_out.F = solutionOut;
        fuel_pump.F = fuel;

        connect(env.outlet, mfc.inlet) annotation (Line(points={{61,-10},{28,
                -10}}, color={0,127,127}));
        connect(mfc.outlet, mixer.inlet) annotation (Line(points={{28,-4},{14,
                -4},{14,10},{-2,10}}, color={0,127,127}));
        connect(pump_out.inlet, mixer.outlet) 
                                          annotation (Line(points={{-34,10},{
                -18,10}}, color={0,127,127}));
        connect(sinkPort.inlet, pump_out.outlet) 
                                             annotation (Line(points={{-56.4,16},
                {-34,16}}, color={0,127,127}));
        connect(pump_in.inlet, anodicLoop.outlet) annotation (Line(points={{29,
                29},{45,29}}, color={0,127,127}));
        connect(pump_in.outlet, mixer.inlet) annotation (Line(points={{29,24},{
                14,24},{14,10},{-2,10}}, color={0,127,127}));
        connect(fuelTank.outlet, fuel_pump.inlet) annotation (Line(points={{10,
                -30},{-10,-30}}, color={0,127,127}));
        connect(fuel_pump.outlet, mixer.fuelInlet) annotation (Line(points={{-10,-24},
                {-10,2}},          color={0,127,127}));
        connect(mixer.envPort, Overflow.inlet) annotation (Line(points={{-10,18},
                {-10,28.4}}, color={0,127,127}));

      initial equation
        mixer.T = T_env;

      end MixerTest;
    end Test;
  end IntegratedOperations;
end Flow;

