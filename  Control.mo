package Control "Controllers for the DMFC system" 
  block CoolerControl "PI controller with saturation" 
    extends Modelica.Blocks.Interfaces.SVcontrol;
    
    import Modelica.SIunits.Time;
    import Modelica.SIunits.VolumeFlowRate;
    
    annotation (defaultComponentName="K",Diagram,
      Icon(
        Line(points=[-80,76; -80,-92],   style(color=8)),
        Polygon(points=[-80,88; -88,66; -72,66; -80,88],     style(color=8,
               fillColor=8)),
        Line(points=[-90,-82; 82,-82],   style(color=8)),
        Polygon(points=[90,-82; 68,-74; 68,-90; 90,-82],     style(color=8,
               fillColor=8)),
        Line(points=[-80,-82; -80,-22; -80,-22; 30,58; 80,58]),
        Text(
          extent=[-20,-22; 80,-62],
          style(color=8),
          string="PI+sat")),
      Documentation(info="<html>
<p>This controller is used in the <tt>AbstractCooler</tt> class to set
the process outlet temperature by manipulating the coolant air flow.</p>
<p>The controller includes a standard Modelica PI controller and a 
saturation limit between zero and a configurable maximum.</p>

<h3>Tuning Procedure</h3>
<p>This controller has been tuned with the Skogestad rules for PI controllers of 
first-order processes.</p>
<p>Running some test runs with the child classes of <tt>AbstractHeatExchangerTest</tt>,
configured with realistic parameters, it was found that in the range of interest of
most variables the process outlet temperature was reduced by 1 K roughly every 
additional 1 L/min of coolant flow. Interpreting this coincidence as a sign of the
heavens, and knowing that feedback will even out any modelling inaccuracies, -1 K/(L/min), 
in SI units 60000 K/(m^3/s), was takes as the value of <em>k</em>.</p>
<p>The value of &tau; (the process time constant) was estimated by considering the
weight and heat capacity of the heat exchangers as provided by IMM, since their steel will 
hold a much larger part of the internal energy compared to the flows, and the heat
capacity of the coolant flow in typical conditions; a value of 10 minutes (600 s)
was found to be an acceptable approximation for &tau;.</p>

<p>The desired response time &tau;<sub>c</sub> was set to 3 minutes, not too aggressive
to avoid variable saturation during transients.</p>

<h3>References</h3>
Skogestad, Sigurd: <em>Simple analytic rules for model reduction and PID controller 
tuning</em>, Journal of Process Control, 13 (2003) 291-309.</p>
</html>"));
  protected 
    Modelica.Blocks.Math.Feedback difference "Calculates the offset" 
      annotation (extent=[-10,-10; 10,10]);
    Modelica.Blocks.Continuous.PI PI(k=Kc, T=tau) "The actual PI controller" 
                                 annotation (extent=[30,-10; 50,10]);
    Modelica.Blocks.Nonlinear.Limiter limiter(uMax=maxFlow, uMin=0) 
      "Limitations of the MFC" annotation (extent=[70,-10; 90,10]);
    
  public 
    parameter Real Kc = -6.6667E-5 "Proportionality constant";
    parameter Time tau = 600 "Integral time";
    parameter VolumeFlowRate maxFlow = 100E-3/60 "Maximum flow rate";
    
  equation 
    connect(u_s, difference.u1) annotation (points=[-120,1.11022e-15; -76,
          1.11022e-15; -76,6.66134e-16; -8,6.66134e-16], style(color=74, rgbcolor=
           {0,0,127}));
    connect(u_m, difference.u2) annotation (points=[-1.11022e-15,-120; 
          -1.11022e-15,-65; 6.66134e-16,-65; 6.66134e-16,-8], style(color=74,
          rgbcolor={0,0,127}));
    connect(difference.y, PI.u) annotation (points=[9,6.10623e-16; 16,
          6.10623e-16; 16,0; 22,0; 22,6.66134e-16; 28,6.66134e-16],
                                                       style(color=74, rgbcolor={
            0,0,127}));
    connect(PI.y, limiter.u) annotation (points=[51,6.10623e-16; 59.5,
          6.10623e-16; 59.5,6.66134e-16; 68,6.66134e-16], style(color=74,
          rgbcolor={0,0,127}));
    connect(limiter.y, y) annotation (points=[91,6.10623e-16; 97.5,6.10623e-16; 
          97.5,5.55112e-16; 110,5.55112e-16], style(color=74, rgbcolor={0,0,127}));
  end CoolerControl;
  annotation (uses(Modelica(version="2.2.1")));
end Control;
