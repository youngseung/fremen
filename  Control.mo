package Control "Controllers for the DMFC system" 
  block CoolerControl "PI controller with saturation" 
    extends Modelica.Blocks.Interfaces.SVcontrol;
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
          string="PI+sat")));
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
