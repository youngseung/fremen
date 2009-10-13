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
  annotation (uses(Modelica(version="3.1")),
    version="1",
    conversion(noneFromVersion=""));
end Electric;
