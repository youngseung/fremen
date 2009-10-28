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
            extent={{-80,40},{80,-40}},
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None), Polygon(
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
</html>"));

    Modelica.Electrical.Analog.Interfaces.PositivePin p "Positive pole" 
      annotation (Placement(transformation(extent={{-60,30},{-40,50}}),
          iconTransformation(extent={{-60,30},{-40,50}})));
    Modelica.Electrical.Analog.Interfaces.NegativePin n "Negative pole" 
      annotation (Placement(transformation(extent={{40,30},{60,50}}),
          iconTransformation(extent={{40,30},{60,50}})));

    parameter Units.Capacity C "Capacity";

    Real SoC(start=0.5) "State of charge";
    Modelica.SIunits.Voltage V = p.v - n.v "Voltage";

  equation
    p.i + n.i = 0; // Charge balance

  end AbstractBattery;

  model Li_ionBattery
    extends AbstractBattery;
  end Li_ionBattery;
end Electric;
