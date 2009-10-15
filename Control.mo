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


package Control "Controllers for the DMFC system"
  block CoolerControl "PI controller with saturation"
    extends Modelica.Blocks.Interfaces.BlockIcon;

    import Modelica.SIunits.Time;
    import Modelica.SIunits.VolumeFlowRate;

    annotation (defaultComponentName="K",Diagram(coordinateSystem(
            preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
          graphics),
      Icon(graphics={
          Line(points={{-80,76},{-80,-92}}, color={192,192,192}),
          Polygon(
            points={{-80,88},{-88,66},{-72,66},{-80,88}},
            lineColor={192,192,192},
            fillColor={192,192,192},
            fillPattern=FillPattern.Solid),
          Line(points={{-90,-82},{82,-82}}, color={192,192,192}),
          Polygon(
            points={{90,-82},{68,-74},{68,-90},{90,-82}},
            lineColor={192,192,192},
            fillColor={192,192,192},
            fillPattern=FillPattern.Solid),
          Line(points={{-80,-82},{-80,-22},{-80,-22},{30,58},{80,58}}),
          Text(
            extent={{-20,-22},{80,-62}},
            lineColor={192,192,192},
            textString=
                 "PI+sat")}),
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
      annotation (Placement(transformation(extent={{-60,-10},{-40,10}},
            rotation=0)));
    Modelica.Blocks.Continuous.PI PI(T=tau, k=Kc) "The actual PI controller" 
                                 annotation (Placement(transformation(extent={{
              20,-10},{40,10}}, rotation=0)));
    Modelica.Blocks.Nonlinear.Limiter limiter(uMax=maxFlow, uMin=minFlow)
      "Limitations of the MFC" annotation (Placement(transformation(extent={{60,
              -10},{80,10}}, rotation=0)));

  public
    parameter Real Kc = 6.6667E-5 "Proportionality constant";
    parameter Time tau = 120 "Integral time";
    parameter VolumeFlowRate minFlow = 0.1E-3/60 "Minimum flow rate";
    parameter VolumeFlowRate maxFlow = 100E-3/60 "Maximum flow rate";

    Flow.IO.TemperatureInput T_r "Reference" 
      annotation (Placement(transformation(extent={{-140,-20},{-100,20}},
            rotation=0)));
    Flow.IO.TemperatureInput T_m "Measurement" 
      annotation (Placement(transformation(
          origin={0,-120},
          extent={{-20,-20},{20,20}},
          rotation=90)));
    Flow.IO.VolumeFlowRateOutput V "Manipulable variable" 
      annotation (Placement(transformation(extent={{100,-20},{140,20}},
            rotation=0)));
  protected
    Modelica.Blocks.Math.Gain signChange(k=-1)
      "Change the sign to get the error signal" 
      annotation (Placement(transformation(extent={{-20,-10},{0,10}}, rotation=
              0)));
  equation
    connect(difference.u2, T_m) annotation (Line(points={{-50,-8},{-50,-40},{0,
            -40},{0,-120},{1.11022e-15,-120}}, color={0,0,127}));
    connect(difference.u1, T_r) annotation (Line(points={{-58,6.66134e-16},{-62,
            6.66134e-16},{-62,1.11022e-15},{-120,1.11022e-15}}, color={0,0,127}));
    connect(limiter.y, V) annotation (Line(points={{81,6.10623e-16},{83.5,
            6.10623e-16},{83.5,1.11022e-15},{120,1.11022e-15}}, color={0,0,127}));
    connect(PI.u, signChange.y) annotation (Line(points={{18,6.66134e-16},{12,
            6.66134e-16},{12,6.10623e-16},{1,6.10623e-16}}, color={0,0,127}));
    connect(signChange.u, difference.y) annotation (Line(points={{-22,
            6.66134e-16},{-24,6.66134e-16},{-24,6.10623e-16},{-41,6.10623e-16}}, 
          color={0,0,127}));
    connect(limiter.u, PI.y) annotation (Line(points={{58,6.66134e-16},{50,
            6.66134e-16},{50,6.10623e-16},{41,6.10623e-16}}, color={0,0,127}));
  end CoolerControl;
  annotation (uses(Modelica(version="3.1"), Flow(version="1")),
                                               Documentation(info="<html>
<p>A collection of controllers for system-wide control and for
some particular units, such as coolers.</p>
</html>"),
    version="1",
    conversion(noneFromVersion=""));
  block CathodeLambdaControl "Feedforward controller for the cathodic flow"
    extends Modelica.Blocks.Interfaces.BlockIcon;

    import Modelica.Constants.R;
    import Modelica.SIunits.MoleFraction;
    import Modelica.SIunits.Concentration;
    import Modelica.SIunits.Pressure;
    import Modelica.SIunits.Temperature;
    import Units.F;

    parameter Integer cells = 1 "Number of cells in the stack";
    parameter Real lambda = 2 "Reactant excess ratio";
    parameter Concentration c_est = 1500
      "Worst-case (highest) anodic concentration";

    parameter Real aA = 8.5E-9 "Partial derivative of n_x wrt. c";
    parameter Real b = 0.21 "Partial derivative of n_x wrt. n_H";

    outer Pressure p_env "Environment pressure";
    outer Temperature T_env "Environment temperature";

    constant MoleFraction x_O2_env = 0.2 "Oxygen molar fraction in air";

    annotation (defaultComponentName="K", Diagram(coordinateSystem(
            preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
          graphics),                               Icon(graphics={Text(
            extent={{100,100},{-100,-100}},
            lineColor={0,0,255},
            textString=
                 "V = f(I)")}),
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
    Flow.IO.CurrentInput I "Current measurement" 
      annotation (Placement(transformation(extent={{-140,-20},{-100,20}},
            rotation=0)));
    Flow.IO.VolumeFlowRateOutput V "Air flow" 
                                annotation (Placement(transformation(extent={{
              100,-20},{140,20}}, rotation=0)));
  protected
    Modelica.Blocks.Math.Add getO2Consumption(k1=(1 - b)/(4*F), k2=(3/2)*aA)
      "Sums the measurement with the appropriate weights" 
      annotation (Placement(transformation(extent={{-20,-20},{20,20}}, rotation
            =0)));
    Modelica.Blocks.Math.Gain convertToFlow(k=cells*lambda/x_O2_env*(R*T_env/
          p_env)) "Obtains air volumetric flow to feed" 
      annotation (Placement(transformation(extent={{50,-10},{70,10}}, rotation=
              0)));
    Modelica.Blocks.Sources.RealExpression c(y=c_est) "Concentration estimate" 
      annotation (Placement(transformation(extent={{-70,-24},{-50,0}}, rotation
            =0)));
  equation
    connect(getO2Consumption.y, convertToFlow.u) 
                           annotation (Line(points={{22,1.22125e-15},{6.5,
            1.22125e-15},{6.5,6.66134e-16},{48,6.66134e-16}}, color={0,0,127}));
    connect(convertToFlow.y, V) 
                       annotation (Line(points={{71,6.10623e-16},{73.5,
            6.10623e-16},{73.5,1.11022e-15},{120,1.11022e-15}}, color={0,0,127}));
    connect(getO2Consumption.u1, I) 
                       annotation (Line(points={{-24,12},{-80,12},{-80,
            1.11022e-15},{-120,1.11022e-15}}, color={0,0,127}));
    connect(c.y, getO2Consumption.u2) annotation (Line(points={{-49,-12},{-24,
            -12}}, color={0,0,127}));
  end CathodeLambdaControl;

  block AnodeLambdaControl "Feedforward controller for the anodic flow"
    extends Modelica.Blocks.Interfaces.BlockIcon;

    import Modelica.SIunits.Concentration;
    import Units.F;

    parameter Integer cells = 1 "Number of cells in the stack";
    parameter Real lambda = 5 "Reactant excess ratio";
    parameter Concentration c_est_an = 1500
      "Worst-case (highest) anodic concentration";
    parameter Concentration c_est_mix = 750
      "Worst-case (lowest) mixer concentration";

    parameter Real aA = 8.5E-9 "Partial derivative of n_x wrt. c";
    parameter Real b = 0.21 "Partial derivative of n_x wrt. n_H";

    annotation (defaultComponentName="K", Diagram(coordinateSystem(
            preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
          graphics),                               Icon(graphics={Text(
            extent={{100,100},{-100,-100}},
            lineColor={0,0,255},
            textString=
                 "V = f(I)")}),
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
    Flow.IO.CurrentInput I "Current measurement" 
      annotation (Placement(transformation(extent={{-140,-20},{-100,20}},
            rotation=0)));
    Flow.IO.VolumeFlowRateOutput V "Solution flow" 
                                annotation (Placement(transformation(extent={{
              100,-20},{140,20}}, rotation=0)));
  protected
    Modelica.Blocks.Math.Add add(            k2=aA, k1=(1 - b)/(6*F))
      "Sums the measurement with the appropriate weights" 
      annotation (Placement(transformation(extent={{-40,-14},{0,26}}, rotation=
              0)));
    Modelica.Blocks.Math.Gain excessRatio(k=cells*lambda)
      "Multiplies by excess ratio" 
      annotation (Placement(transformation(extent={{20,-4},{40,16}}, rotation=0)));
    Modelica.Blocks.Math.Division divide
      "Divides methanol molar flow with concentration" 
      annotation (Placement(transformation(extent={{60,-10},{80,10}}, rotation=
              0)));
    Modelica.Blocks.Sources.RealExpression c_mix(y=c_est_mix)
      "Concentration estimate in the mixer" 
      annotation (Placement(transformation(extent={{20,-10},{40,-30}}, rotation
            =0)));
    Modelica.Blocks.Sources.RealExpression c_an(y=c_est_an)
      "Concentration estimate in anode" 
      annotation (Placement(transformation(extent={{-80,4},{-60,-16}}, rotation
            =0)));
  equation
    connect(add.u1, I) annotation (Line(points={{-44,18},{-90,18},{-90,0},{-106,
            0},{-106,1.11022e-15},{-120,1.11022e-15}}, color={0,0,127}));
    connect(divide.y, V) annotation (Line(points={{81,6.10623e-16},{87.5,
            6.10623e-16},{87.5,1.11022e-15},{120,1.11022e-15}}, color={0,0,127}));
    connect(excessRatio.y, divide.u1) 
                               annotation (Line(points={{41,6},{45.25,6},{45.25,
            6},{49.5,6},{49.5,6},{58,6}}, color={0,0,127}));
    connect(c_mix.y, divide.u2) 
                            annotation (Line(points={{41,-20},{50,-20},{50,-6},
            {58,-6}}, color={0,0,127}));
    connect(add.y, excessRatio.u) annotation (Line(points={{2,6},{6,6},{6,6},{
            10,6},{10,6},{18,6}}, color={0,0,127}));
    connect(c_an.y, add.u2) annotation (Line(points={{-59,-6},{-51.5,-6},{-51.5,
            -6},{-44,-6}}, color={0,0,127}));
  end AnodeLambdaControl;

  block ReferenceFuelControl "Feedforward controller for the fuel flow"
    extends Modelica.Blocks.Interfaces.BlockIcon;

    import Modelica.SIunits.Concentration;
    import Modelica.SIunits.Temperature;
    import Thermo.Species;
    import Thermo.Phases;
    import Thermo.K;
    import Thermo.mw;
    import Thermo.rho;
    import Units.F;

    outer Temperature T_env "Environment temperature";

    parameter Integer cells = 1 "Number of cells in the stack";
    parameter Concentration c_ref = 1000 "Concentration set point";

    parameter Real aA = 8.5E-9 "Partial derivative of I_x wrt. c";
    parameter Real b = 0.21 "Partial derivative of I_x wrt. I";

    annotation (defaultComponentName="K", Diagram(coordinateSystem(
            preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics),
                                                   Icon(graphics={Text(
            extent={{100,100},{-100,-100}},
            lineColor={0,0,255},
            textString="V = f(I,T)")}),
      Documentation(info="<html>
<p>This feedforward controller takes current and degasser temperature and 
returns the appropriate volumetric flow rate of methanol from
the fuel tank.</p>
<p>Two parameters are the estimates for <tt>aA</tt> and <tt>b</tt>, necessary
to estimate the extent of cross-over current in the cell to compensate for.</p>
</html>"));
    Flow.IO.CurrentInput I "Current measurement" 
      annotation (Placement(transformation(extent={{-140,40},{-100,80}}, rotation=
             0)));
    Flow.IO.VolumeFlowRateOutput V "Solution flow" 
                                annotation (Placement(transformation(extent={{100,
              -20},{140,20}}, rotation=0)));
    Flow.IO.TemperatureInput T_deg "Degasser temperature" 
      annotation (Placement(transformation(extent={{-140,-80},{-100,-40}},
            rotation=0)));

  protected
    Real fM "Degasser loss factor";
    Real n_to_V "Conversion factor from methanol moles to volume";

  equation
    fM = K(T_deg, Species.Methanol)*mw(Species.Water)/rho(T_deg, Species.Water, Phases.Liquid)/(6*F)/(1-K(T_deg,Species.Water));
    n_to_V = mw(Species.Methanol)/rho(T_env, Species.Methanol, Phases.Liquid);

    V = n_to_V * ((1-b)/6/F*I + aA*c_ref + fM*I*c_ref)*cells;

  end ReferenceFuelControl;

  block WaterControl "Feedback controller for solution level"
    extends Modelica.Blocks.Interfaces.BlockIcon;

    import Modelica.Constants.R;
    import g = Modelica.Constants.g_n;
    import Modelica.SIunits.Area;
    import Modelica.SIunits.Pressure;
    import Modelica.SIunits.Temperature;
    import Modelica.SIunits.Time;
    import Thermo.dp_h2o_dt;
    import Thermo.mw;
    import Thermo.Species;

    outer Temperature T_env "Environment temperature";

    parameter Area A = 50E-4 "Mixer cross-sectional area";
    parameter Temperature T_0 = 325 "Nominal condenser temperature";
    parameter Time tau = 600 "Desired response time";

    annotation (defaultComponentName="K", Diagram(coordinateSystem(
            preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics),
                                                   Icon(graphics={
          Line(points={{-80,76},{-80,-92}}, color={192,192,192}),
          Polygon(
            points={{-80,88},{-88,66},{-72,66},{-80,88}},
            lineColor={192,192,192},
            fillColor={192,192,192},
            fillPattern=FillPattern.Solid),
          Line(points={{-90,-82},{82,-82}}, color={192,192,192}),
          Polygon(
            points={{90,-82},{68,-74},{68,-90},{90,-82}},
            lineColor={192,192,192},
            fillColor={192,192,192},
            fillPattern=FillPattern.Solid),
          Line(points={{-80,-82},{-80,-82},{-80,-82},{80,58},{80,58}}),
          Text(
            extent={{-20,-22},{80,-62}},
            lineColor={192,192,192},
            textString="P")}),
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
    Flow.IO.TemperatureInput T_cond "Condenser temperature" 
      annotation (Placement(transformation(extent={{-140,-80},{-100,-40}},
            rotation=0)));
    Flow.IO.TemperatureOutput T_ref "Target condenser temperature" 
      annotation (Placement(transformation(extent={{100,-20},{140,20}}, rotation=0)));
    Flow.IO.PressureInput p_mix "Hydrostatic level measurement" 
      annotation (Placement(transformation(extent={{-140,40},{-100,80}}, rotation=
             0)));
    Flow.IO.VolumeFlowRateInput V_cath "Cathodic volumetric flow" 
      annotation (Placement(transformation(extent={{-140,-20},{-100,20}},
            rotation=0)));
  protected
    discrete Pressure p_0 "Initial hydrostatic pressure";
    Pressure delta_p = p_mix - p_0 "Pressure variation";

  equation
    T_ref = T_0 + 1/(tau * V_cath/R/T_env * dp_h2o_dt(T_cond,1)) * A/mw(Species.Water)/g * delta_p;

    when initial() then
      p_0 = p_mix;
    end when;

  end WaterControl;

  block TemperatureControl "PID for temperature control"
    extends Modelica.Blocks.Interfaces.BlockIcon;

    import Modelica.SIunits.Time;
    import Modelica.SIunits.VolumeFlowRate;
    import Units.Temperature;

    annotation (defaultComponentName="K",Diagram(coordinateSystem(
            preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
          graphics),
      Icon(graphics={
          Line(points={{-80,76},{-80,-92}}, color={192,192,192}),
          Polygon(
            points={{-80,88},{-88,66},{-72,66},{-80,88}},
            lineColor={192,192,192},
            fillColor={192,192,192},
            fillPattern=FillPattern.Solid),
          Line(points={{-90,-82},{82,-82}}, color={192,192,192}),
          Polygon(
            points={{90,-82},{68,-74},{68,-90},{90,-82}},
            lineColor={192,192,192},
            fillColor={192,192,192},
            fillPattern=FillPattern.Solid),
          Line(points={{-80,-82},{-80,32},{-80,-20},{30,58},{80,58}}),
          Text(
            extent={{-20,-22},{80,-62}},
            lineColor={192,192,192},
            textString=
                 "PID+sat")}),
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
<p>Note that this controller implements a simple anti-windup strategy; if the manipulated
variable is reported to be saturated (through <tt>isSaturated</tt>), the integrator is
frozen.</p>
<h3>References</h3>
<p>Skogestad, Sigurd: <em>Simple analytic rules for model reduction and PID controller 
tuning</em>, Journal of Process Control, 13 (2003) 291-309.</p>
</html>"));

  public
    parameter Real Kc = 0.625 "Proportionality constant";
    parameter Time tau_I = 75 "Integral time";
    parameter Time tau_D = 60 "Derivative time";

    parameter Temperature T_deg_0 = 315 "Nominal degasser temperature";
    parameter Temperature T_FC_ref = 333
      "Set-point for the fuel cell's temperature";

    Flow.IO.TemperatureInput T_m "Measurement" 
      annotation (Placement(transformation(extent={{-140,-20},{-100,20}},
            rotation=0)));
    Flow.IO.TemperatureOutput T_deg_ref "Manipulable variable" 
      annotation (Placement(transformation(extent={{100,0},{140,40}}, rotation=
              0)));
    Modelica.Blocks.Interfaces.BooleanInput isSaturated
      "Whether input is saturated" annotation (Placement(transformation(extent=
              {{140,0},{100,-40}}, rotation=0)));

  protected
    Real int "Integral of the error";
    Temperature e = T_FC_ref -T_m "Measured error";

  equation
    T_deg_ref = T_deg_0 + Kc * (e + int/tau_I + der(e)*tau_D);

    // Anti-windup in case of saturation
    der(int) = if isSaturated then 0 else e;

  end TemperatureControl;

  model MingledTemperatureControl
    "Stack-temperature control for case of mingled outlets"
    extends Modelica.Blocks.Interfaces.BlockIcon;

    annotation (defaultComponentName="K",Diagram(coordinateSystem(
            preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics),
                                                  Icon(graphics={
          Line(points={{-76,-14},{82,-14},{82,-14}}, color={255,0,0}),
          Line(points={{-76,-78},{-76,-14},{-46,-14},{34,62},{84,62}}),
          Polygon(
            points={{-76,92},{-84,70},{-68,70},{-76,92}},
            lineColor={192,192,192},
            fillColor={192,192,192},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-76,-16},{84,-58}},
            lineColor={192,192,192},
            textString="PI+sat+lambda"),
          Polygon(
            points={{94,-78},{72,-70},{72,-86},{94,-78}},
            lineColor={192,192,192},
            fillColor={192,192,192},
            fillPattern=FillPattern.Solid),
          Line(points={{-86,-78},{86,-78}}, color={192,192,192}),
          Line(points={{-76,80},{-76,-88}}, color={192,192,192})}),
      Documentation(info="<html>
<p>This controller sets the anodic flow rate to maintain a certain cell
temperature. The controller will however maintain a minimum &lambda;, to 
prioritise the reaction to cell heating.</p>
<p>Internally, the controller uses a Skogestad PI controller, whose 
integrator is disconnected when overridden by the &lambda; control to 
avoid wind-up.</p>
</html>"));

    import Modelica.SIunits.Concentration;
    import Modelica.SIunits.Temperature;
    import Modelica.SIunits.Time;
    import Modelica.SIunits.Volume;
    import Modelica.SIunits.VolumeFlowRate;
    import Modelica.Constants.eps;
    import Flow.IO.CurrentInput;
    import Flow.IO.TemperatureInput;
    import Flow.IO.VolumeFlowRateOutput;
    import Units.F;

    CurrentInput I "Current in cell" annotation (Placement(transformation(extent={
              {-140,-80},{-100,-40}}, rotation=0)));
    TemperatureInput T_mix "Mixer's temperature" 
      annotation (Placement(transformation(extent={{-140,-20},{-100,20}},
            rotation=0)));
    TemperatureInput T_stack "Stack's temperature" 
      annotation (Placement(transformation(extent={{-140,40},{-100,80}}, rotation=
             0)));
    VolumeFlowRateOutput V "Anodic-loop flow" 
      annotation (Placement(transformation(extent={{100,-20},{140,20}}, rotation=0)));

    parameter VolumeFlowRate V_max = 1.5E-6 "Maximum flow allowed by pump";
    parameter Integer n = 1 "Number of cells in stack";
    parameter Temperature T_r = 330 "Target stack temperature";
    parameter Real lambda = 2 "Minimum lambda";
    parameter Concentration c = 800 "Worst-case (lowest) concentration";
    parameter Time tau_c = 600 "Desired closed-loop response";
    parameter Volume V_cp = 5.95E-6
      "Volume of solution with same heat capacity as one cell";
    parameter Real aA = 8.5E-9 "Partial derivative of n_x wrt. c";
    parameter Real b = 0.21 "Partial derivative of n_x wrt. n_H";

  protected
    Real K_c = n * V_cp * tau_c / second^2 "PI-controller gain";
    Time tau_I = 4*tau_c "PI-controller integrating time";

    VolumeFlowRate V_PI "Flow given by the PI controller";
    VolumeFlowRate V_lambda "Minimum (lambda-limited) flow";

    Temperature DeltaT "Temperature difference between stack and mixer";

    Real x "Internal integrator";
    Real e "Measured error";

    constant Time second = 1 "To get dimensions right";

  equation
    V = min(V_max, max(V_lambda, V_PI));
    V_lambda = lambda * n * (aA*c + (1-b) * I / (6*F)) / c;
    V_PI * DeltaT = K_c * (e + x/second);

    // Avoid division-by-zero errors
    DeltaT = if noEvent(abs(T_stack-T_mix) < eps) then  eps else T_stack-T_mix;

    e = T_stack - T_r;
    // Anti-windup strategy: integrate only when PI active
    der(x) = if V_PI < V_lambda or V_PI > V_max then 0 else e;

    assert( V>=0, "==> Negative flow resulting from control");

  end MingledTemperatureControl;

  block MingledFuelControl "Feedforward controller for the fuel flow"
    extends Modelica.Blocks.Interfaces.BlockIcon;

    import Modelica.SIunits.Concentration;
    import Modelica.SIunits.MoleFraction;
    import Modelica.SIunits.Pressure;
    import Modelica.SIunits.Temperature;
    import Modelica.Constants.R;
    import Thermo.Species;
    import Thermo.Phases;
    import Thermo.K;
    import Thermo.mw;
    import Thermo.rho;
    import Units.MolarFlow;
    import Units.F;

    outer Temperature T_env "Environment temperature";
    outer Pressure p_env "Environment pressure";

    parameter Integer cells = 1 "Number of cells in the stack";
    parameter Concentration c_ref = 1000 "Concentration set point";
    parameter Real lambda = 2 "Cathodic air excess ratio";

    parameter Real aA = 8.5E-9 "Partial derivative of I_x wrt. c";
    parameter Real b = 0.21 "Partial derivative of I_x wrt. I";

    annotation (defaultComponentName="K", Diagram(coordinateSystem(
            preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics),
                                                   Icon(graphics={Text(
            extent={{100,100},{-100,-100}},
            lineColor={0,0,255},
            textString="V = f(V,I,T)")}),
      Documentation(info="<html>
<p>This feedforward controller takes cathodic inflow, current and separator
temperature, and returns the appropriate volumetric flow rate of methanol from
the fuel tank.</p>
<p>Two parameters are the estimates for <tt>aA</tt> and <tt>b</tt>, necessary
to estimate the extent of cross-over current in the cell to compensate for.</p>
</html>"));
    Flow.IO.VolumeFlowRateInput V_cath "Cathodic inflow" 
      annotation (Placement(transformation(extent={{-140,40},{-100,80}}, rotation=
             0)));
    Flow.IO.TemperatureInput T_sep "Separator temperature" 
      annotation (Placement(transformation(extent={{-140,-80},{-100,-40}},
            rotation=0)));
    Flow.IO.CurrentInput I annotation (Placement(transformation(extent={{-140,-20},
              {-100,20}}, rotation=0)));
    Flow.IO.VolumeFlowRateOutput V "Solution flow" 
                                annotation (Placement(transformation(extent={{100,
              -20},{140,20}}, rotation=0)));

  protected
    constant MoleFraction y_O2 = 0.2 "Molar fraction of oxygen in air";
    MolarFlow n_cath = V_cath * p_env / R / T_env "Cathodic molar inflow";
    MoleFraction x "Approximate molar fraction of methanol";
    Real fM "Separator loss factor";
    Real n_to_V "Conversion factor from methanol moles to volume";

  equation
    x = c_ref * mw(Species.Water)/rho(T_sep, Species.Water, Phases.Liquid);
    fM = (1-y_O2/3/lambda)*K(T_sep, Species.Methanol)/(1-K(T_sep, Species.Water));
    n_to_V = mw(Species.Methanol)/rho(T_env, Species.Methanol, Phases.Liquid);

    V = n_to_V * ((1-b)/6/F*I + aA*c_ref)*cells + n_to_V * fM*n_cath*x;

  end MingledFuelControl;
end Control;
