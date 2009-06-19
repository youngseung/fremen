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
saturation limit between configurable minimum and maximum. The minimum must
be set to a little higher than zero to avoid initialisation problems.</p>

<h3>Tuning Procedure</h3>
<p>This controller has been tuned with the Skogestad rules for PI controllers of 
first-order processes.</p>
<p>Running some test runs with the child classes of <tt>AbstractHeatExchangerTest</tt>,
configured with realistic parameters, it was found that in the range of interest of
most variables the process outlet temperature was reduced by 1 K roughly every 
additional 1 L/min of coolant flow. Interpreting this coincidence as a sign of the
heavens, and knowing that feedback will even out any modelling inaccuracies, -1 K/(L/min), 
in SI units 60000 K/(m^3/s), was taken as the value of <em>k</em>.</p>
<p>The value of &tau; (the process time constant) was estimated by considering the
weight and heat capacity of the heat exchangers as provided by IMM, since their steel will 
hold a much larger part of the internal energy compared to the flows, and the heat
capacity of the coolant flow in typical conditions; a value of 10 minutes (600 s)
was found to be an acceptable approximation for &tau;.</p>

<p>The desired response time &tau;<sub>c</sub> was set to 3 minutes, not too aggressive
to avoid variable saturation during transients.</p>

<h3>References</h3>
<p>Skogestad, Sigurd: <em>Simple analytic rules for model reduction and PID controller 
tuning</em>, Journal of Process Control, 13 (2003) 291-309.</p>
</html>"));
  protected 
    Modelica.Blocks.Math.Feedback difference "Calculates the offset" 
      annotation (extent=[-60,-10; -40,10]);
    Modelica.Blocks.Continuous.PI PI(T=tau, k=Kc) "The actual PI controller" 
                                 annotation (extent=[20,-10; 40,10]);
    Modelica.Blocks.Nonlinear.Limiter limiter(uMax=maxFlow, uMin=minFlow) 
      "Limitations of the MFC" annotation (extent=[60,-10; 80,10]);
    
  public 
    parameter Real Kc = 6.6667E-5 "Proportionality constant";
    parameter Time tau = 600 "Integral time";
    parameter VolumeFlowRate minFlow = 0.1E-3/60 "Minimum flow rate";
    parameter VolumeFlowRate maxFlow = 100E-3/60 "Maximum flow rate";
    
    Flow.TemperatureInput T_r "Reference" 
      annotation (extent=[-140,-20; -100,20]);
    Flow.TemperatureInput T_m "Measurement" 
      annotation (extent=[-20,-140; 20,-100], rotation=90);
    Flow.VolumeFlowRateOutput V "Manipulable variable" 
      annotation (extent=[100,-20; 140,20]);
  protected 
    Modelica.Blocks.Math.Gain signChange(k=-1) 
      "Change the sign to get the error signal" 
      annotation (extent=[-20,-10; 0,10]);
  equation 
    connect(difference.u2, T_m) annotation (points=[-50,-8; -50,-40; 0,-40; 0,
          -120; 1.11022e-15,-120], style(color=74, rgbcolor={0,0,127}));
    connect(difference.u1, T_r) annotation (points=[-58,6.66134e-16; -62,
          6.66134e-16; -62,1.11022e-15; -120,1.11022e-15], style(color=74,
          rgbcolor={0,0,127}));
    connect(limiter.y, V) annotation (points=[81,6.10623e-16; 83.5,6.10623e-16; 
          83.5,1.11022e-15; 120,1.11022e-15], style(color=74, rgbcolor={0,0,127}));
    connect(PI.u, signChange.y) annotation (points=[18,6.66134e-16; 12,
          6.66134e-16; 12,6.10623e-16; 1,6.10623e-16], style(color=74, rgbcolor=
           {0,0,127}));
    connect(signChange.u, difference.y) annotation (points=[-22,6.66134e-16; 
          -24,6.66134e-16; -24,6.10623e-16; -41,6.10623e-16], style(color=74,
          rgbcolor={0,0,127}));
    connect(limiter.u, PI.y) annotation (points=[58,6.66134e-16; 50,6.66134e-16; 
          50,6.10623e-16; 41,6.10623e-16], style(color=74, rgbcolor={0,0,127}));
  end CoolerControl;
  annotation (uses(Modelica(version="2.2.1")), Documentation(info="<html>
<p>A collection of controllers for system-wide control and for
some particular units, such as coolers.</p>
</html>"));
  block CathodeLambdaControl "Feedforward controller for the cathodic flow" 
    extends Modelica.Blocks.Interfaces.BlockIcon;
    
    import Modelica.Constants.R;
    import Modelica.SIunits.MoleFraction;
    import Units.Temperature;
    import Modelica.SIunits.Pressure;
    import Units.F;
    
    parameter Real lambda = 2 "Reactant excess ratio";
    
    parameter Real aA = 5E-3 "Partial derivative of I_x wrt. c";
    parameter Real b = 0.21 "Partial derivative of I_x wrt. I";
    
    outer parameter Pressure p_env = 101325 "Environment pressure";
    outer parameter Temperature T_env = 298.15 "Environment temperature";
    
    annotation (defaultComponentName="K", Diagram, Icon(Text(
          extent=[100,100; -100,-100],
          style(color=3, rgbcolor={0,0,255}),
          string="V = f(I,c)")), 
      Documentation(info="<html>
<p>This feedforward controller uses the current and the anodic 
methanol concentration (either a measurement or a value) to set
the volumetric cathodic inflow to a fuel cell.</p>
<p>The provided concentration value is supposed to be the one
in the anodic bulk, usually assumed to be the outlet one.</p>
<p>The controller allows to set an excess ratio &lambda;; other
parameters are the estimates for <tt>aA</tt> and <tt>b</tt>, necessary
to estimate the extent of cross-over current in the cell to compensate
for.</p>
</html>"));
    Flow.ConcentrationInput c "Concentration target" 
      annotation (extent=[-140,-80; -100,-40]);
    Flow.CurrentInput I "Current measurement" 
      annotation (extent=[-140,40; -100,80]);
    Flow.VolumeFlowRateOutput V "Air flow" 
                                annotation (extent=[100,-20; 140,20]);
  protected 
    constant MoleFraction x_O2_env = 0.2 "Oxygen molar fraction in air";
    
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
    connect(add.u1, I) annotation (points=[-22,6; -60,6; -60,60; -120,60],
        style(color=74, rgbcolor={0,0,127}));
    connect(add.u2, c) annotation (points=[-22,-6; -60,-6; -60,-60; -120,-60],
        style(color=74, rgbcolor={0,0,127}));
  end CathodeLambdaControl;
  
  block AnodeLambdaControl "Feedforward controller for the anodic flow" 
    extends Modelica.Blocks.Interfaces.BlockIcon;
    
    import Units.F;
    
    parameter Real lambda = 2 "Reactant excess ratio";
    
    parameter Real aA = 5E-3 "Partial derivative of I_x wrt. c";
    parameter Real b = 0.21 "Partial derivative of I_x wrt. I";
    
    annotation (defaultComponentName="K", Diagram, Icon(Text(
          extent=[100,100; -100,-100],
          style(color=3, rgbcolor={0,0,255}),
          string="V = f(I,c)")), 
      Documentation(info="<html>
<p>This feedforward controller uses the current and the anodic 
methanol concentration (either a measurement or a value) to set
the volumetric anodic inflow to a fuel cell.</p>
<p>The provided concentration value is supposed to be the one
in the anodic bulk, usually assumed to be the outlet one.</p>
<p>The controller allows to set an excess ratio &lambda;; other
parameters are the estimates for <tt>aA</tt> and <tt>b</tt>, necessary
to estimate the extent of cross-over current in the cell to compensate
for.</p>
</html>"));
    Flow.ConcentrationInput c "Concentration estimate" 
      annotation (extent=[-140,-80; -100,-40]);
    Flow.CurrentInput I "Current measurement" 
      annotation (extent=[-140,40; -100,80]);
    Flow.VolumeFlowRateOutput V "Solution flow" 
                                annotation (extent=[100,-20; 140,20]);
  protected 
    Modelica.Blocks.Math.Add add(k1=(1 - b), k2=aA) 
      "Sums the measurement with the appropriate weights" 
      annotation (extent=[-40,-4; -20,16]);
    Modelica.Blocks.Math.Gain gain(k=lambda/(6*F)) "Multiplies by excess ratio"
      annotation (extent=[0,-4; 20,16]);
    Modelica.Blocks.Math.Division divide 
      "Divides methanol molar flow with concentration" 
      annotation (extent=[50,-10; 70,10]);
  equation 
    connect(add.y, gain.u) annotation (points=[-19,6; -14.75,6; -14.75,6; -10.5,
          6; -10.5,6; -2,6],                style(color=74, rgbcolor={0,0,127}));
    connect(add.u1, I) annotation (points=[-42,12; -60,12; -60,60; -120,60],
        style(color=74, rgbcolor={0,0,127}));
    connect(add.u2, c) annotation (points=[-42,5.55112e-16; -60,5.55112e-16; 
          -60,-60; -120,-60],
        style(color=74, rgbcolor={0,0,127}));
    connect(divide.y, V) annotation (points=[71,6.10623e-16; 87.5,6.10623e-16; 
          87.5,1.11022e-15; 120,1.11022e-15], style(color=74, rgbcolor={0,0,127}));
    connect(gain.y, divide.u1) annotation (points=[21,6; 27.75,6; 27.75,6; 34.5,
          6; 34.5,6; 48,6],
                         style(color=74, rgbcolor={0,0,127}));
    connect(divide.u2, c) annotation (points=[48,-6; 36,-6; 36,-60; -120,-60],
        style(color=74, rgbcolor={0,0,127}));
  end AnodeLambdaControl;
  
  block FuelControl "Feedforward controller for the fuel flow" 
    extends Modelica.Blocks.Interfaces.BlockIcon;
    
    import Units.Temperature;
    import Thermo.Molecules.Methanol;
    import Thermo.Molecules.Water;
    import Thermo.Phases.Liquid;
    import Thermo.K;
    import Thermo.mw;
    import Thermo.rho;
    import Units.F;
    
    outer parameter Temperature T_env = 298.15 "Environment temperature";
    
    parameter Real aA = 5E-3 "Partial derivative of I_x wrt. c";
    parameter Real b = 0.21 "Partial derivative of I_x wrt. I";
    
    annotation (defaultComponentName="K", Diagram, Icon(Text(
          extent=[100,100; -100,-100],
          style(color=3, rgbcolor={0,0,255}),
          string="V = f(I,c,T)")), 
      Documentation(info="<html>
<p>This feedforward controller takes current, concentration and degasser
temperature and returns the appropriate volumetric flow rate of methanol from
the fuel tank.</p>
<p>The provided concentration value is supposed to be the one
in the anodic bulk, usually assumed to be the outlet one.</p>
<p>Two parameters are the estimates for <tt>aA</tt> and <tt>b</tt>, necessary
to estimate the extent of cross-over current in the cell to compensate for.</p>
</html>"));
    Flow.ConcentrationInput c "Concentration estimate" 
      annotation (extent=[-140,-20; -100,20]);
    Flow.CurrentInput I "Current measurement" 
      annotation (extent=[-140,40; -100,80]);
    Flow.VolumeFlowRateOutput V "Solution flow" 
                                annotation (extent=[100,-20; 140,20]);
    Flow.TemperatureInput T_deg "Degasser temperature" 
      annotation (extent=[-140,-80; -100,-40]);
    
  protected 
    Real fM "Degasser loss factor";
    Real n_to_V "Conversion factor from methanol moles to volume";
    
  equation 
    fM = K(T_deg, Methanol)*mw(Water)/rho(T_deg, Water, Liquid)/(6*F)/(1-K(T_deg,Water));
    n_to_V = mw(Methanol)/rho(T_env, Methanol, Liquid);
    
    V = n_to_V * (((1-b)*I + aA*c)/(6*F) + fM*I*c);
    
  end FuelControl;
  
  block WaterControl "Feedback controller for solution level" 
    extends Modelica.Blocks.Interfaces.BlockIcon;
    
    import Modelica.Constants.R;
    import g = Modelica.Constants.g_n;
    import Modelica.SIunits.Area;
    import Modelica.SIunits.Pressure;
    import Units.Temperature;
    import Modelica.SIunits.Time;
    import Thermo.dp_h2o_dt;
    import Thermo.mw;
    import Thermo.Molecules.Water;
    
    outer parameter Temperature T_env = 298.15 "Environment temperature";
    
    parameter Area A = 50E-4 "Mixer cross-sectional area";
    parameter Temperature T_0 = 325 "Nominal condenser temperature";
    parameter Time tau = 600 "Desired response time";
    
    annotation (defaultComponentName="K", Diagram, Icon(
        Line(points=[-80,76; -80,-92],   style(color=8)),
        Polygon(points=[-80,88; -88,66; -72,66; -80,88],     style(color=8,
               fillColor=8)),
        Line(points=[-90,-82; 82,-82],   style(color=8)),
        Polygon(points=[90,-82; 68,-74; 68,-90; 90,-82],     style(color=8,
               fillColor=8)),
        Line(points=[-80,-82; -80,-82; -80,-82; 80,58; 80,58]),
        Text(
          extent=[-20,-22; 80,-62],
          style(color=8),
          string="P")), 
      Documentation(info="<html>
<p>This controller takes the current hydrostatic pressure in the mixer,
the current volumetric flow in the cathode and the condenser temperature
to return a new reference for the condenser temperature, which is then
cascaded to the condenser controller.</p>
<p>The feedback part is in the measurement of the hydrostatic pressure,
which is (roughly) proportional to mixer solution content. At simulation
start, the controller memorises the initial value of this pressure and uses 
this initial value as a set-point.</p>
<p>The cathodic flow is not a measurement, but rather a value passed by the
cathodic flow &lambda; controller.</p>
<p>The condenser temperature is a relatively simple measurement, and is used
to calculate the vapour pressure of water in the condenser, which may change
significantly in the temperature range of interest.</p>
<p>An available parameter is &tau;, which sets the desired response time: this
should be set to at least 10 minutes (600 s), to allow the internal, cascaded
cooler controller time to settle. It is also possible to set the nominal 
target temperature and the mixer cross-sectional area, to convert hydrostatic
pressure in a volumetric estimate.</p>
<p>The controller is a P controller whose proportionality constant is obtained
from Skogestad's tuning rules.</p>

<h3>References</h3>
<p>Zenith, Federico, and Krewer, Ulrike: <em>Dynamics and Control of a DMFC 
System</em>, 7<sup>th</sup> Fuel Cell Science, Engineering and Technology 
Conference, June 2009, Newport Beach, USA.</p>
<p>Skogestad, Sigurd: <em>Simple analytic rules for model reduction and PID controller 
tuning</em>, Journal of Process Control, 13 (2003) 291-309.</p>
</html>"));
  public 
    Flow.TemperatureInput T_cond "Condenser temperature" 
      annotation (extent=[-140,-80; -100,-40]);
    Flow.TemperatureOutput T_ref "Target condenser temperature" 
      annotation (extent=[100,-20; 140,20]);
    Flow.PressureInput p_mix "Hydrostatic level measurement" 
      annotation (extent=[-140,40; -100,80]);
    Flow.VolumeFlowRateInput V_cath "Cathodic volumetric flow" 
      annotation (extent=[-140,-20; -100,20]);
  protected 
    discrete Pressure p_0 "Initial hydrostatic pressure";
    Pressure delta_p = p_mix - p_0 "Pressure variation";
    
  equation 
    T_ref = T_0 + 1/(tau * V_cath/R/T_env * dp_h2o_dt(T_cond,1)) * A/mw(Water)/g * delta_p;
    
    when initial() then
      p_0 = p_mix;
    end when;
    
  end WaterControl;
  
  block TemperatureControl "PID for temperature control" 
    extends Modelica.Blocks.Interfaces.BlockIcon;
    
    import Modelica.SIunits.Time;
    import Modelica.SIunits.VolumeFlowRate;
    import Units.Temperature;
    
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
<p>This PID controller tries to make the fuel-cell temperature converge to a given
set-point by manipulating the degasser reference temperature around a given nominal
value.</p>
<p>The effect of degasser temperature on fuel-cell temperature is assumed to be a unit-gain
process with two lags. The unit gain is justified since most of the flow exiting the anodic
side is the main heat input into the cell, so one degree more in the degasser will eventually
translate to one degree more in the fuel cell (in reality it will be a bit less).
The first lag is due to the solution in the mixer, which when
assumed to be 5 ml of (mostly) water results in a lag of 60 seconds. The second lag
is due to the material of the fuel-cell graphite plates, resulting in about 300 seconds.</p>
<p>From these two lags, the PID parameters are calculated with the Skogestad rules.</p>

<h3>References</h3>
<p>Skogestad, Sigurd: <em>Simple analytic rules for model reduction and PID controller 
tuning</em>, Journal of Process Control, 13 (2003) 291-309.</p>
</html>"));
    
  public 
    parameter Real Kc = 0.5 "Proportionality constant";
    parameter Time tau_I = 600 "Integral time";
    parameter Time tau_D = 30 "Derivative time";
    
    parameter Temperature T_0 = 315 "Nominal degasser temperature";
    parameter Temperature T_FC_ref = 333 
      "Set-point for the fuel cell's temperature";
    
    Flow.TemperatureInput T_m "Measurement" 
      annotation (extent=[-140,-20; -100,20]);
    Flow.TemperatureOutput T_deg_ref "Manipulable variable" 
      annotation (extent=[100,-20; 140,20]);
  protected 
    Modelica.Blocks.Sources.RealExpression nominalTemperature(y=T_0) 
      "The nominal temperature of the degasser" 
      annotation (extent=[20,0; 40,20]);
    Modelica.Blocks.Math.Add add annotation (extent=[60,-10; 80,10]);
    Modelica.Blocks.Continuous.PID PID(
      k=Kc,
      Ti=tau_I,
      Td=tau_D) annotation (extent=[20,-20; 40,0]);
  protected 
    Modelica.Blocks.Sources.RealExpression setPoint(y=T_FC_ref) 
      "The set-point for the fuel-cell temperature" 
      annotation (extent=[-80,10; -60,30]);
    Modelica.Blocks.Math.Feedback feedback annotation (extent=[-50,10; -30,30]);
  equation 
    connect(add.y, T_deg_ref) annotation (points=[81,6.10623e-16; 82,
          6.10623e-16; 82,1.11022e-15; 120,1.11022e-15], style(color=74,
          rgbcolor={0,0,127}));
    connect(nominalTemperature.y, add.u1) annotation (points=[41,10; 50,10; 50,
          6; 58,6], style(color=74, rgbcolor={0,0,127}));
    connect(PID.y, add.u2) annotation (points=[41,-10; 50,-10; 50,-6; 58,-6],
        style(color=74, rgbcolor={0,0,127}));
    connect(feedback.u1, setPoint.y) annotation (points=[-48,20; -59,20], style(
          color=74, rgbcolor={0,0,127}));
    connect(feedback.u2, T_m) annotation (points=[-40,12; -40,1.11022e-15; -120,
          1.11022e-15], style(color=74, rgbcolor={0,0,127}));
    connect(feedback.y, PID.u) annotation (points=[-31,20; -20,20; -20,-10; 18,
          -10], style(color=74, rgbcolor={0,0,127}));
  end TemperatureControl;
end Control;
