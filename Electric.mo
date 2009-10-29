within ;
package Electric "Components for electric interaction"
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
  /*
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
  end Test;
end Electric;
