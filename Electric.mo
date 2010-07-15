within ;
              /**
 * © Federico Zenith, 2009.
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


package Electric "Components for electric interaction"
  model CellEmulator "A Thevenin emulator of a fuel cell"
    extends Modelica.Electrical.Analog.Interfaces.OnePort;

    parameter Modelica.SIunits.Voltage V = 5 "Open-circuit voltage";
    parameter Modelica.SIunits.Resistance R = 0.3 "Internal resistance";

    annotation (
      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
              100}}), graphics={Rectangle(
            extent={{-100,42},{100,-40}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            lineThickness=1), Text(
            extent={{-100,20},{100,-20}},
            lineColor={0,0,0},
            lineThickness=1,
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            textString="DMFC")}),
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}),
              graphics),
      Documentation(info="<html>
<p>A trivial fuel-cell emulator based on a Th&eacute;venin equivalent circuit.</p>
<p>The parameters (5&nbsp;V, 0.3&nbsp;&Omega;) are taken for the design values of the
cell stack ordered from Baltic Fuel Cells.</p>
</html>"));
  equation

    v = V + R*i;

  end CellEmulator;

  model PowerLoad "A typical laptop power load"
    extends Modelica.Blocks.Sources.CombiTimeTable(
    table = transpose({ {0, 10, 10, 20, 20, 30, 30, 36, 36, 40, 40, 49, 49, 60, 60, 80, 80, 110, 111, 120, 120, 130, 130, 140, 140, 160, 160, 180},
             {24, 24, 26, 26, 24, 24, 35, 35, 25, 25, 22, 22, 36, 36, 25, 25, 36,  36,  24,  24,  40,  40,  32,  32,  30,  30,  40,  40}}),
    extrapolation = Modelica.Blocks.Types.Extrapolation.Periodic);
    annotation (Documentation(info="<html>
<p>This data was produced by Martin Behrendt and is unsourced; probably just
a good guess.</p>
</html>"));
  end PowerLoad;
  annotation (uses(Modelica(version="3.1"), Units(version="1")),
    version="1",
    conversion(noneFromVersion=""));
  partial model AbstractBattery "A generic battery"

    annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
              -100},{100,100}}), graphics={Rectangle(
            extent={{-100,60},{100,-60}},
            pattern=LinePattern.None,
            lineColor={0,0,0},
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid), Polygon(
            points={{-18,-30},{-10,10},{0,6},{4,26},{16,26},{4,-6},{-6,-4},{-18,
                -30}},
            lineColor={0,0,0},
            smooth=Smooth.None,
            fillColor={255,255,0},
            fillPattern=FillPattern.Solid)}), Diagram(coordinateSystem(
            preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics),
      Documentation(info="<html>
<p>A generic battery with a state of charge, a voltage and a capacity, 
whose completion is left to child classes.</p>
<p>To instantiate a new class, it is necessary to specify:</p>
<ul>
<li>How voltage V is defined, typically as a function of SoC and current;</li>
<li>How the SoC is influenced by other variables, typically current;</li>
</ul>
<p>The definition of the relationship of the SoC with other variables could
also account for self-discharge effects.</p>
</html>"));

    Modelica.Electrical.Analog.Interfaces.PositivePin p "Positive pole" 
      annotation (Placement(transformation(extent={{-60,48},{-40,68}}),
          iconTransformation(extent={{-60,48},{-40,68}})));
    Modelica.Electrical.Analog.Interfaces.NegativePin n "Negative pole" 
      annotation (Placement(transformation(extent={{40,48},{60,68}}),
          iconTransformation(extent={{40,48},{60,68}})));

    parameter Units.Capacity C = 3600 "Capacity";

    Units.StateOfCharge SoC(start=0.5, min=0, max=1) "State of charge";
    Modelica.SIunits.Voltage V = p.v - n.v "Voltage";
    Modelica.SIunits.Current I = n.i "Current";

  equation
    p.i + n.i = 0; // Charge balance

  end AbstractBattery;

  partial model NoSelfDischargeBattery
    "A battery model without any self-discharge"
    extends AbstractBattery;

  equation
    C * der(SoC) = I;

    annotation (Documentation(info="<html>
<p>This battery implements a trivial model of the state of charge,
which assumes that no charge is being lost and no self-discharge is 
present.</p>
<p>To get a working battery, one has to specify the relationship 
defining the voltage.</p>
</html>"));
  end NoSelfDischargeBattery;

  model Li_ionBattery
    extends NoSelfDischargeBattery;

    import Modelica.Constants.R;
    import Modelica.Math.exp;
    import Modelica.Math.tanh;
    import Modelica.Math.asinh;
    import Modelica.SIunits.Area;
    import Modelica.SIunits.ElectricPotential;
    import Modelica.SIunits.PotentialDifference;
    import Modelica.SIunits.StoichiometricNumber;
    import Modelica.SIunits.CurrentDensity;
    import Units.F;

    parameter Modelica.SIunits.Temperature T = 298.15;

    parameter Integer cells = 1 "Number of cells in series";

    parameter Area A = 2.8E-3 "Electrode area";
    parameter CurrentDensity i0_a = 1.1 "Exchange current at anode";
    parameter CurrentDensity i0_c = 0.8 "Exchange current at cathode";

    ElectricPotential E0_a "Open-circuit potential at anode";
    ElectricPotential E0_c "Open-circuit potential at cathode";

    StoichiometricNumber x "Composition at the anode";
    StoichiometricNumber y "Composition at the cathode";

    PotentialDifference eta_a "Anodic overpotential";
    PotentialDifference eta_c "Cathodic overpotential";

  protected
    parameter Real k = 0.5 * F / R / T "Handy constant";

  equation
    V = cells * (E0_c - E0_a + eta_c - eta_a);
  /* FIXME why is this commented out?
  E0_a = -0.132 + 1.41*exp(-3.52*x);
  E0_c = 4.06279 + 0.0677504*tanh(-21.8502*y+12.8268)
       - 0.105734*(1/(1.00167-y)^0.379571 - 1.576)
       - 0.045*exp(-71.69*y^8) + 0.01*exp(-200*(y-0.19));
*/
    E0_a = -0.16 + 1.32*exp(-3*x) + 10*exp(-2000*x);
    E0_c = 4.19829 + 0.0565661*tanh(-14.5546*y+8.60942)
         - 0.0275479*(1/(0.998432-y)^0.492465 - 1.90111)
         - 0.157123*exp(-0.04738*y^8) + 0.810239*exp(-40*(y-0.133875));

    x = 0.0045 + 0.5655*SoC;
    y = 0.80   - 0.63  *SoC;

    eta_a = asinh(I/(A*i0_a)/2) / k;
    eta_c = asinh(I/(A*i0_c)/2) / k;

    annotation (Documentation(info="<html>
<p>This model is a lithium-ion battery as that modelled by Martin Behrendt.
The model for the overall open-circuit voltage is exceedingly parametrised,
but it works.</p>
<h4>References</h4>
<p>Doyle, M.; Newman, J.; Gozdz, A. N.; Schmutz, C. N. and Tarascon, J.: 
<em>Comparison of Modeling Predictions with Experimental Data from Plastic 
Lithium Ion Cells</em>, Journal of the Electrochemical Society, 1996, 143, 
1890&ndash;1903.</p>
</html>"));
  end Li_ionBattery;

  model PulseWidthModulator
    "Given a duty ratio, produces a pulse-width-modulation signal"

    Modelica.Blocks.Interfaces.RealInput D "Fraction of time switch is closed" 
      annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
    Modelica.Blocks.Interfaces.BooleanOutput switch
      "Whether to close the switch" 
      annotation (Placement(transformation(extent={{100,-20},{140,20}})));
    annotation (defaultComponentName="pwm", Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,
              -100},{100,100}}),
                           graphics), Icon(coordinateSystem(preserveAspectRatio=true,
            extent={{-100,-100},{100,100}}), graphics={Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,0},
            lineThickness=1,
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid), Text(
            extent={{-100,40},{100,-40}},
            lineColor={0,0,0},
            lineThickness=1,
            textString="PWM")}));

    parameter Modelica.SIunits.Time T = 1E-6 "Duty cycle";

  equation
    assert( T > 0, "Duty cycle must be a positive time");
    switch = rem(time, T) < D*T;

  end PulseWidthModulator;

  partial model Converter "A generic DC/DC converter"
    Modelica.SIunits.Voltage v1 "Voltage drop over the left port";
    Modelica.SIunits.Voltage v2 "Voltage drop over the right port";
    Modelica.SIunits.Current i1
      "Current flowing from pos. to neg. pin of the left port";
    Modelica.SIunits.Current i2
      "Current flowing from pos. to neg. pin of the right port";

    Modelica.Electrical.Analog.Interfaces.PositivePin p1
      "Positive pin of the left port (potential p1.v > n1.v for positive voltage drop v1)"
                                                                                                        annotation (Placement(
          transformation(extent={{-110,40},{-90,60}}, rotation=0)));
    Modelica.Electrical.Analog.Interfaces.NegativePin n1
      "Negative pin of the left port"              annotation (Placement(
          transformation(extent={{-90,-60},{-110,-40}}, rotation=0)));
    Modelica.Electrical.Analog.Interfaces.PositivePin p2
      "Positive pin of the right port (potential p2.v > n2.v for positive voltage drop v2)"
                                                                                                         annotation (Placement(
          transformation(extent={{110,40},{90,60}}, rotation=0)));
    Modelica.Electrical.Analog.Interfaces.NegativePin n2
      "Negative pin of the right port"              annotation (Placement(
          transformation(extent={{90,-60},{110,-40}}, rotation=0)));
    Modelica.Blocks.Interfaces.RealInput D "Duty ratio" annotation (Placement(
          transformation(
          extent={{-20,-20},{20,20}},
          rotation=270,
          origin={0,120}), iconTransformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={0,110})));
    annotation (
      Diagram(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={
          Polygon(
            points={{-120,53},{-110,50},{-120,47},{-120,53}},
            lineColor={160,160,164},
            fillColor={160,160,164},
            fillPattern=FillPattern.Solid),
          Line(points={{-136,50},{-111,50}}, color={160,160,164}),
          Polygon(
            points={{127,-47},{137,-50},{127,-53},{127,-47}},
            lineColor={160,160,164},
            fillColor={160,160,164},
            fillPattern=FillPattern.Solid),
          Line(points={{111,-50},{136,-50}}, color={160,160,164}),
          Text(
            extent={{112,-44},{128,-29}},
            lineColor={160,160,164},
            textString="i2"),
          Text(
            extent={{118,52},{135,67}},
            lineColor={0,0,0},
            textString="i2"),
          Polygon(
            points={{120,53},{110,50},{120,47},{120,53}},
            lineColor={0,0,0},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={160,160,164}),
          Line(points={{111,50},{136,50}}, color={0,0,0}),
          Line(points={{-136,-49},{-111,-49}}, color={160,160,164}),
          Polygon(
            points={{-126,-46},{-136,-49},{-126,-52},{-126,-46}},
            lineColor={160,160,164},
            fillColor={160,160,164},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-127,-46},{-110,-31}},
            lineColor={160,160,164},
            textString="i1"),
          Text(
            extent={{-136,53},{-119,68}},
            lineColor={160,160,164},
            textString="i1")}),
      Documentation(revisions="<html>
<ul>
<li><i> 1998   </i>
       by Christoph Clauss<br> initially implemented<br>
       </li>
</ul>
</html>", info="<html>

</html>"),
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}}), graphics={
          Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,0},
            lineThickness=1,
            fillPattern=FillPattern.Solid,
            fillColor={255,255,255}),
          Text(
            extent={{-100,100},{0,0}},
            lineColor={0,0,0},
            lineThickness=1,
            textString="="),
          Text(
            extent={{0,0},{100,-100}},
            lineColor={0,0,0},
            lineThickness=1,
            textString="="),
          Line(
            points={{-100,-100},{100,100}},
            color={0,0,0},
            thickness=1,
            smooth=Smooth.None)}));
  equation
    v1 = p1.v - n1.v;
    v2 = p2.v - n2.v;
    i1 = p1.i;
    i2 = p2.i;

  end Converter;

  model BuckBoost "Buck-boost converter"
    extends Converter;
    annotation (
      Diagram(graphics),
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}}), graphics),
      Documentation(info="<html>
<p>This class implements a canonical buck-boost converter, which can be controlled
by setting a duty ratio <i>D</i>, converted appropriately by an internal modulator
in opening/closing signals for the internal switch.</p>
</html>"));
    Modelica.Electrical.Analog.Basic.Inductor inductance(L=10.6E-3) 
                                                               annotation (
        Placement(transformation(
          extent={{-16,-16},{16,16}},
          rotation=270,
          origin={0,24})));
    Modelica.Electrical.Analog.Ideal.IdealClosingSwitch switch 
      annotation (Placement(transformation(extent={{-70,30},{-30,70}})));
    PulseWidthModulator pwm(T=1E-5) 
                            annotation (Placement(transformation(
          extent={{-6,-6},{6,6}},
          rotation=270,
          origin={-50,80})));
    Modelica.Electrical.Analog.Ideal.IdealDiode diode annotation (Placement(
          transformation(
          extent={{-20,-20},{20,20}},
          rotation=180,
          origin={50,50})));
    Modelica.Electrical.Analog.Basic.Resistor resistance(R=0.12) annotation (
        Placement(transformation(
          extent={{-16,-16},{16,16}},
          rotation=270,
          origin={0,-24})));
  equation
    connect(inductance.p, switch.n) 
                                  annotation (Line(
        points={{3.36094e-15,40},{3.36094e-15,50},{-30,50}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(pwm.switch, switch.control) annotation (Line(
        points={{-50,72.8},{-50,64}},
        color={255,0,255},
        smooth=Smooth.None));
    connect(pwm.D, D) annotation (Line(
        points={{-50,87.2},{-50,94},{0,94},{0,120},{1.11022e-15,120}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(p2, diode.p) annotation (Line(
        points={{100,50},{70,50}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(diode.n, inductance.p) 
                                 annotation (Line(
        points={{30,50},{3.36094e-15,50},{3.36094e-15,40}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(p1, switch.p) annotation (Line(
        points={{-100,50},{-70,50}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(resistance.p, inductance.n) annotation (Line(
        points={{3.36094e-15,-8},{-2.51717e-15,-8},{-2.51717e-15,8}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(resistance.n, n1) annotation (Line(
        points={{-2.51717e-15,-40},{0,-40},{0,-50},{-100,-50}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(resistance.n, n2) annotation (Line(
        points={{-2.51717e-15,-40},{0,-40},{0,-50},{100,-50}},
        color={0,0,255},
        smooth=Smooth.None));
  end BuckBoost;

  model IdealCuk "An ideal, lossless C'uk converter"
    extends Converter;
    PulseWidthModulator pwm annotation (Placement(transformation(
          extent={{-6,-6},{6,6}},
          rotation=180,
          origin={0,0})));
    Modelica.Electrical.Analog.Basic.Inductor inductor_in(L=10E-6) 
      annotation (Placement(transformation(extent={{-80,30},{-40,70}})));
    Modelica.Electrical.Analog.Basic.Capacitor capacitor(C=1E-3) 
      annotation (Placement(transformation(extent={{-20,30},{20,70}})));
    Modelica.Electrical.Analog.Ideal.IdealClosingSwitch switch 
      annotation (Placement(transformation(extent={{-20,-20},{20,20}},
          rotation=270,
          origin={-30,0})));
    annotation (Diagram(graphics));
    Modelica.Electrical.Analog.Ideal.IdealDiode diode annotation (Placement(
          transformation(
          extent={{-20,-20},{20,20}},
          rotation=270,
          origin={30,0})));
    Modelica.Electrical.Analog.Basic.Inductor inductor_out(L=10E-6) 
      annotation (Placement(transformation(extent={{40,30},{80,70}})));
  equation
    connect(inductor_in.p, p1) annotation (Line(
        points={{-80,50},{-100,50}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(switch.p, inductor_in.n) annotation (Line(
        points={{-30,20},{-30,50},{-40,50}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(pwm.D, D) annotation (Line(
        points={{7.2,-1.17037e-15},{16,-1.17037e-15},{16,80},{1.11022e-15,80},{
            1.11022e-15,120}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(switch.control, pwm.switch) annotation (Line(
        points={{-16,-1.68349e-15},{-7.2,-1.68349e-15},{-7.2,5.93059e-16}},
        color={255,0,255},
        smooth=Smooth.None));
    connect(capacitor.p, switch.p) annotation (Line(
        points={{-20,50},{-30,50},{-30,20}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(n1, switch.n) annotation (Line(
        points={{-100,-50},{-30,-50},{-30,-20}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(switch.n, diode.n) annotation (Line(
        points={{-30,-20},{-30,-50},{30,-50},{30,-20}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(n2, diode.n) annotation (Line(
        points={{100,-50},{30,-50},{30,-20}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(capacitor.n, inductor_out.p) annotation (Line(
        points={{20,50},{40,50}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(inductor_out.n, p2) annotation (Line(
        points={{80,50},{100,50}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(diode.p, capacitor.n) annotation (Line(
        points={{30,20},{30,50},{20,50}},
        color={0,0,255},
        smooth=Smooth.None));
  end IdealCuk;

  package Test "Test package"

    model NoSelfDischargeBatteryTest

      model MyNoSelfDischargeBattery
        extends NoSelfDischargeBattery;
      equation
        V = 1.1 + 0.2*SoC;
      end MyNoSelfDischargeBattery;

      MyNoSelfDischargeBattery myNoSelfDischargeBattery(C=1000) 
        annotation (Placement(transformation(extent={{-40,-80},{40,0}})));
      Modelica.Electrical.Analog.Basic.Ground ground 
        annotation (Placement(transformation(extent={{56,-20},{76,0}})));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                -100},{100,100}}),      graphics));
      Modelica.Electrical.Analog.Sources.SineCurrent sineCurrent(
        freqHz=0.001,
        I=2,
        phase=1.5707963267949) 
        annotation (Placement(transformation(extent={{-20,20},{20,60}})));
    equation
      connect(myNoSelfDischargeBattery.n, ground.p) annotation (Line(
          points={{20,-16.8},{20,5.55112e-16},{66,5.55112e-16}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(sineCurrent.n, myNoSelfDischargeBattery.n) annotation (Line(
          points={{20,40},{20,25.8},{20,25.8},{20,11.6},{20,-16.8},{20,-16.8}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(sineCurrent.p, myNoSelfDischargeBattery.p) annotation (Line(
          points={{-20,40},{-20,-16.8}},
          color={0,0,255},
          smooth=Smooth.None));
    end NoSelfDischargeBatteryTest;

    model Li_ionBatteryTest

      Modelica.Electrical.Analog.Basic.Ground ground 
        annotation (Placement(transformation(extent={{56,-20},{76,0}})));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                -100},{100,100}}),      graphics), experiment(StopTime=1000));
      Modelica.Electrical.Analog.Sources.SineCurrent sineCurrent(
        freqHz=0.001,
        I=2,
        phase=1.5707963267949) 
        annotation (Placement(transformation(extent={{-20,20},{20,60}})));
      Li_ionBattery li_ionBattery 
        annotation (Placement(transformation(extent={{-40,-80},{40,0}})));
    equation
      connect(li_ionBattery.n, ground.p) annotation (Line(
          points={{20,-16.8},{20,0},{66,0},{66,5.55112e-16}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(sineCurrent.n, li_ionBattery.n) annotation (Line(
          points={{20,40},{20,25.8},{20,25.8},{20,11.6},{20,-16.8},{20,-16.8}},
          color={0,0,255},
          smooth=Smooth.None));

      connect(sineCurrent.p, li_ionBattery.p) annotation (Line(
          points={{-20,40},{-20,-16.8}},
          color={0,0,255},
          smooth=Smooth.None));
    end Li_ionBatteryTest;

    model TestPulseWidthModulator "A test for the pulse-width modulator class"

      PulseWidthModulator pwm(T=0.01) 
        annotation (Placement(transformation(extent={{-40,-40},{40,40}})));

    equation
      pwm.D = time-0.1;

      annotation (experiment(StopTime=1.2, NumberOfIntervals=2000));
    end TestPulseWidthModulator;

    model TestConverter "A test case for a generic converter model"

      CellEmulator cellEmulator annotation (Placement(transformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={-80,0})));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                -100},{100,100}}),
                             graphics), experiment(StopTime=0.005,
            NumberOfIntervals=1000000));
      Modelica.Electrical.Analog.Basic.Resistor resistor(R=10) 
                                                              annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={62,0})));
      Modelica.Electrical.Analog.Basic.Ground ground 
        annotation (Placement(transformation(extent={{-100,-60},{-80,-40}})));
      Modelica.Electrical.Analog.Basic.Capacitor inletCapacitor(C=0.5e-3) 
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-50,0})));
      Modelica.Electrical.Analog.Basic.Capacitor outletCapacitor(C=0.5e-3) 
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={84,0})));
      Modelica.Blocks.Sources.Step D(
        height=0.2,
        offset=0.4,
        startTime=0.002) 
        annotation (Placement(transformation(extent={{-40,40},{-20,60}})));
      replaceable Converter converter 
        annotation (Placement(transformation(extent={{-18,-18},{18,18}})));
    equation
      connect(ground.p, cellEmulator.n) annotation (Line(
          points={{-90,-40},{-80,-40},{-80,-20}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(outletCapacitor.p, resistor.p) annotation (Line(
          points={{84,-10},{84,-28},{62,-28},{62,-10}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(outletCapacitor.n, resistor.n) annotation (Line(
          points={{84,10},{84,28},{62,28},{62,10}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(converter.p2, resistor.n) annotation (Line(
          points={{18,9},{40,9},{40,28},{62,28},{62,10}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(resistor.p, converter.n2) annotation (Line(
          points={{62,-10},{62,-28},{40,-28},{40,-9},{18,-9}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(cellEmulator.p, inletCapacitor.p) annotation (Line(
          points={{-80,20},{-80,28},{-50,28},{-50,10}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(cellEmulator.n, inletCapacitor.n) annotation (Line(
          points={{-80,-20},{-80,-28},{-50,-28},{-50,-10}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(converter.n1, inletCapacitor.n) annotation (Line(
          points={{-18,-9},{-40,-9},{-40,-28},{-50,-28},{-50,-10}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(converter.p1, inletCapacitor.p) annotation (Line(
          points={{-18,9},{-40,9},{-40,28},{-50,28},{-50,10}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(converter.D, D.y) annotation (Line(
          points={{-5.66214e-16,19.8},{-5.66214e-16,50},{-19,50}},
          color={0,0,127},
          smooth=Smooth.None));
    end TestConverter;

    model TestBuckBoost "Test for a buck-boost converter"
      extends TestConverter(redeclare BuckBoost converter(pwm(T=1E-4)),
                                                           D(startTime=0.5),
        inletCapacitor(C=10E-3),
        outletCapacitor(C=10E-3));
      annotation (experiment(Interval=1e-05));
    end TestBuckBoost;

    model TestIdealCuk "Test case for an ideal Cuk converter"
      extends TestConverter(redeclare IdealCuk converter(pwm(T=1E-5)), D(
            startTime=0.005));
      annotation (experiment(StopTime=0.005, NumberOfIntervals=100000));
    end TestIdealCuk;
  end Test;
end Electric;
