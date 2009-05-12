package Control "Controllers for the DMFC system" 
  block CoolerControl "PI controller with saturation" 
    extends Modelica.Blocks.Interfaces.BlockIcon;
    
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
      annotation (extent=[-50,-10; -30,10]);
    Modelica.Blocks.Continuous.PI PI(k=Kc, T=tau) "The actual PI controller" 
                                 annotation (extent=[-10,-10; 10,10]);
    Modelica.Blocks.Nonlinear.Limiter limiter(uMax=maxFlow, uMin=0) 
      "Limitations of the MFC" annotation (extent=[40,-10; 60,10]);
    
  public 
    parameter Real Kc = -6.6667E-5 "Proportionality constant";
    parameter Time tau = 600 "Integral time";
    parameter VolumeFlowRate maxFlow = 100E-3/60 "Maximum flow rate";
    
    Flow.TemperatureInput T_r "Reference"
      annotation (extent=[-140,-20; -100,20]);
    Flow.TemperatureInput T_m "Measurement"
      annotation (extent=[-20,-140; 20,-100], rotation=90);
    Flow.VolumeFlowRateOutput V "Manipulable variable"
      annotation (extent=[100,-20; 140,20]);
  equation 
    connect(difference.y, PI.u) annotation (points=[-31,6.10623e-16; -14,
          6.10623e-16; -14,0; -8,0; -8,6.66134e-16; -12,6.66134e-16],
                                                       style(color=74, rgbcolor={
            0,0,127}));
    connect(PI.y, limiter.u) annotation (points=[11,6.10623e-16; 29.5,
          6.10623e-16; 29.5,6.66134e-16; 38,6.66134e-16], style(color=74,
          rgbcolor={0,0,127}));
    connect(difference.u2, T_m) annotation (points=[-40,-8; -40,-40; 0,-40; 0,
          -120; 1.11022e-15,-120], style(color=74, rgbcolor={0,0,127}));
    connect(difference.u1, T_r) annotation (points=[-48,6.66134e-16; -62,
          6.66134e-16; -62,1.11022e-15; -120,1.11022e-15], style(color=74, 
          rgbcolor={0,0,127}));
    connect(limiter.y, V) annotation (points=[61,6.10623e-16; 83.5,6.10623e-16; 
          83.5,1.11022e-15; 120,1.11022e-15], style(color=74, rgbcolor={0,0,127}));
  end CoolerControl;
  annotation (uses(Modelica(version="2.2.1")));
  block CathodeLambdaControl "Feedforward controller for the cathodic flow" 
    extends Modelica.Blocks.Interfaces.BlockIcon;
    
    import Modelica.Constants.R;
    import Modelica.SIunits.MoleFraction;
    import Modelica.SIunits.Temperature;
    import Modelica.SIunits.Pressure;
    
    parameter Real lambda = 2 "Reactant excess ratio";
    
    parameter Real aA "Partial derivative of I_x wrt. c";
    parameter Real b "Partial derivative of I_x wrt. I";
    
    parameter MoleFraction x_O2_env = 0.2 
      "Oxygen molar fraction in environment";
    parameter Pressure p_env = 101325 "Environment pressure";
    parameter Tempearture T_env = 298.15 "Environment temperature";
    
    annotation (Diagram, Icon(Text(
          extent=[100,100; -100,-100], 
          style(color=3, rgbcolor={0,0,255}), 
          string="V = f(I,c)")));
    Flow.ConcentrationInput c "Estimate of concentration"
      annotation (extent=[-140,-70; -100,-30]);
    Flow.CurrentInput I "Measurement of current"
      annotation (extent=[-140,30; -100,70]);
    Flow.VolumeFlowRateOutput V annotation (extent=[100,-20; 140,20]);
  protected 
    constant Modelica.SIunits.FaradayConstant F = 96485.3415;
    
    Modelica.Blocks.Math.Add add(k1=(1 - b), k2=aA) 
      "Sums the measurement with the appropriate weights"
      annotation (extent=[-20,-10; 0,10]);
    Modelica.Blocks.Math.Gain gain(k=lambda/(4*F)/x_O2_env*(R*T_env/p_env))
      annotation (extent=[40,-10; 60,10]);
  equation 
    connect(add.y, gain.u) annotation (points=[1,6.10623e-16; 6.5,6.10623e-16; 
          6.5,6.66134e-16; 38,6.66134e-16], style(color=74, rgbcolor={0,0,127}));
    connect(gain.y, V) annotation (points=[61,6.10623e-16; 73.5,6.10623e-16; 
          73.5,1.11022e-15; 120,1.11022e-15], style(color=74, rgbcolor={0,0,127}));
    connect(add.u1, I) annotation (points=[-22,6; -60,6; -60,50; -120,50], 
        style(color=74, rgbcolor={0,0,127}));
    connect(add.u2, c) annotation (points=[-22,-6; -60,-6; -60,-50; -120,-50], 
        style(color=74, rgbcolor={0,0,127}));
  end CathodeLambdaControl;
end Control;
