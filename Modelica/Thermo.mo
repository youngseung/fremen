/**
 * © Federico Zenith, 2008-2009.
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


encapsulated package Thermo "Thermodynamic library"
  import Modelica;

  import Units.Temperature;
  import Modelica.SIunits.Pressure;
  import Modelica.SIunits.MolarInternalEnergy;
  import Modelica.SIunits.MolarHeatCapacity;
  import Modelica.SIunits.MolarEntropy;
  import Modelica.SIunits.PartialPressure;
  import Modelica.SIunits.Density;
  import Modelica.SIunits.MolarMass;
  import Modelica.SIunits.MoleFraction;
  import Modelica.SIunits.DiffusionCoefficient;
  import Units.MolarEnthalpy;

public
  type Species = enumeration(
      Methanol,
      Water,
      Oxygen,
      CarbonDioxide,
      Nitrogen) 
    annotation (Documentation(info="<html>
<p>This enumeration contains the species considered in our model.</p>
</html>"));

public
  constant Species[:] Condensables = Species.Methanol : Species.Water
    "All species that can be in both phases";
  constant Species[:] Incondensables = Species.Oxygen : Species.Nitrogen
    "All gaseous, incondensable species";

public
  type Phases = enumeration(
      Gas,
      Liquid)                                           annotation (Documentation(info="<html>
<p>This enumeration contains the phases considered in the model.</p>
</html>"));

protected
  record WaterVapourParameters "Parameters for water vapour"
    constant Real A=8.22;
    constant Real B=0.00015;
    constant Real C=0.00000134;
    constant Real D=10326.2823;
    constant Real T_min=300.0;
    constant Real T_max=2500.0;
    annotation (Documentation(info="<html>
<p>These parameters are necessary to calculate the enthalpy or specific heat of
water vapour. Results are reliable between 300 and 2500 kelvin.<br/>
Source: Perry's Chemical Engineers' Handbook, 7th edition, page 2-163.</p>
<p>Note: parameter D has been added to those from Perry, so that enthalpy would be 0 at 25&deg;C.</p>
</html>"));
  end WaterVapourParameters;

protected
  record GaseousMethanolParameters
    constant Real C1=0.3925E5;
    constant Real C2=0.879E5;
    constant Real C3=1.9165E3;
    constant Real C4=0.5365E5;
    constant Real C5=896.7;
    constant Real T_min=200.0;
    constant Real T_max=1500.0;
    annotation (Documentation(info="<html>
<p>These parameters are necessary to calculate the enthalpy or specific heat of
gaseous methanol. Results are reliable between 200 and 1500 kelvin.<br/>
Source: Perry's Chemical Engineers' Handbook, 7th edition, page 2-179.</p>
</html>"));
  end GaseousMethanolParameters;

protected
  record LiquidMethanolParameters
    constant Real A=1.058E5;
    constant Real B=-362.23;
    constant Real C=0.9379;
    constant Real D=23730.2384;
    constant Real T_min=175.47;
    constant Real T_max=400.0;
    annotation (Documentation(info="<html>
<p>These parameters are necessary to calculate the enthalpy or specific heat of
liquid methanol. Results are reliable between 175.47 and 400 kelvin.<br/>
Source: Perry's Chemical Engineers' Handbook, 7th edition, page 2-171.</p>
<p>Note: parameter D has been added to those from Perry, so that enthalpy would be 0 at 25&deg;C.</p>
</html>"));
  end LiquidMethanolParameters;

protected
  package Shomate "Parameters for shomate equations"

    record Parameters "Set of parameters for Shomate equations"
      Real A;
      Real B;
      Real C;
      Real D;
      Real E;
      Real F;
      Real G;
      Real H;
      Real T_min "Minimum temperature for data validity";
      Real T_max "Maximum temperature for data validity";
      annotation (Documentation(info="<html>
</html>"));
    end Parameters;

    constant Parameters O2( A=29.65900, B=6.137261, C=-1.186521, D=0.095780, E=-0.219663,
                            F=-9.861391, G=237.9480, H=0.000000, T_min=298.0, T_max=6000.0)
      "Shomate parameters for oxygen";

    constant Parameters N2( A=26.09200, B=8.218801, C=-1.976141, D=0.159274, E=0.044434,
                            F=-7.989230, G=221.0200, H=0.000000, T_min=298.0, T_max=6000.0)
      "Shomate parameters for nitrogen";

    constant Parameters CO2( A=24.99735, B=55.18696, C=-33.69137, D=7.948387, E=-0.136638,
                             F=-403.6075, G=228.2431, H=-393.5224, T_min=298.0, T_max=1200.0)
      "Shomate parameters for carbon dioxide";

    constant Parameters H2O( A=-203.6060, B=1523.290, C=-3196.413, D=2474.455, E=3.855326,
                             F=-256.5478, G=-488.7163, H=-285.8304, T_min=298.0, T_max=500.0)
      "Shomate parameters for liquid water";
    annotation (Documentation(info="<html>
<p>This package contains the parameters necessary for Shomate equations, used to
calculate the molar enthalpy, the molar heat capacity and the molar entropy of many
components.</p>
<p>Shomate parameters are instances of the same record structure <tt>Shomate.Parameters</tt>
and are set to be <tt>constant</tt>.</p>
<p>Data are from NIST, in detail:</p>
<ul>
<li><a href=\"http://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Units=SI&Mask=2#Thermo-Condensed\">Liquid water</a>;</li>
<li><a href=\"http://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447&Units=SI&Mask=1#Thermo-Gas\">Oxygen</a>;</li>
<li><a href=\"http://webbook.nist.gov/cgi/cbook.cgi?ID=C124389&Units=SI&Mask=1#Thermo-Gas\">Carbon dioxide</a>;</li>
<li><a href=\"http://webbook.nist.gov/cgi/cbook.cgi?ID=C7727379&Units=SI&Mask=1#Thermo-Gas\">Nitrogen</a>.</li>
</ul>
</html>"));
  end Shomate;

public
  function speciesName
    input Species species "The species code.";
    output String name "The species name.";
  algorithm
    if species == Species.Methanol then
      name := "methanol";
    elseif species == Species.Water then
      name := "water";
    elseif species == Species.Oxygen then
      name := "oxygen";
    elseif species == Species.CarbonDioxide then
      name := "carbon dioxide";
    elseif species == Species.Nitrogen then
      name := "nitrogen";
    else
      name := "unknown molecule ("+String(species)+")";
    end if;

    annotation (Documentation(info="<html>
<p>This function returns a string with the name of the species whose
code was set into the input.</p>
</html>"));
  end speciesName;

public
  function phaseName

    import Units.Phase;

    input Phases phase "The phase code.";
    output String name "The phase name.";
  algorithm
    if phase == Phases.Gas then
      name := "gas";
    elseif phase == Phases.Liquid then
      name := "liquid";
    else
      name := "unknown phase ("+String(phase)+")";
    end if;

    annotation (Documentation(info="<html>
<p>This function returns a string with the name of the phase whose
code was set into the input.</p>
</html>"));
  end phaseName;

protected
  function ShomateEnthalpy
    input Temperature T;
    input Shomate.Parameters p;
    output MolarEnthalpy h;
  protected
    Real t;
  algorithm
    t := T/1000.0;
    h := (p.A*t + (p.B*t*t)/2 + (p.C*t*t*t)/3 + (p.D*t*t*t*t)/4 - p.E/t + p.F - p.H)*1000;
    annotation (Documentation(info="<html>
<p>This abstract function uses the Shomate parameters to calculate the enthalpy.
Note that the Shomate parameters are still unspecified.</p>
<p>The exponents are obtained by subsequent multiplications, since this function 
is called relatively often.</p>
</html>"));
  end ShomateEnthalpy;

protected
  function ShomateHeatCapacity
    input Temperature T;
    input Shomate.Parameters p;
    output MolarHeatCapacity cp;
  protected
    Real t;
  algorithm
    t := T/1000.0;
    cp := p.A + p.B*t + p.C*t*t + p.D*t*t*t + p.E/t/t;
    annotation (Documentation(info="<html>
<p>This abstract function uses the Shomate parameters to calculate the specific
heat. Note that the Shomate parameters are still unspecified.</p>
</html>"));
  end ShomateHeatCapacity;

protected
  function h_h2o_gas "The enthalpy of gaseous water."
    input Temperature T;
    output MolarEnthalpy h;
  protected
    constant WaterVapourParameters p;
  algorithm
    // NOTE conversion from calories.
    h := (p.A*T + (p.B*T*T)/2 + (p.C*T*T*T)/3)*4.184 - p.D;
  end h_h2o_gas;

protected
  function cp_h2o_gas "The specific heat of gaseous water."
    input Temperature T;
    output MolarHeatCapacity cp;
  protected
    constant WaterVapourParameters p;
  algorithm
    // NOTE conversion from calories.
    cp := (p.A + p.B*T + p.C*T*T)*4.184;
  end cp_h2o_gas;

protected
  function h_ch3oh_liq "The enthalpy of liquid methanol."
    input Temperature T;
    output MolarEnthalpy h;
  protected
    constant LiquidMethanolParameters p;
  algorithm
    h := (p.A*T + (p.B*T*T)/2 + (p.C*T*T*T)/3)/1000 - p.D;
  end h_ch3oh_liq;

protected
  function cp_ch3oh_liq "The specific heat of liquid methanol."
    input Temperature T;
    output MolarHeatCapacity cp;
  protected
    constant LiquidMethanolParameters p;
  algorithm
    cp := (p.A + p.B*T + p.C*T*T)/1000;
  end cp_ch3oh_liq;

protected
  function cp_ch3oh_gas
    input Temperature T;
    output MolarHeatCapacity cp;
  protected
    constant GaseousMethanolParameters p;
  algorithm
    cp := (p.C1 + p.C2*(p.C3/T/(sinh(p.C3/T)))^2 + p.C4*(p.C5/T/cosh(p.C5/T))^2)/1000;
    annotation (Documentation(info="<html>
<p>This function returns the specific heat of gaseous methanol at a given temperature.</p>
<p>The formula used in this function is the one reported in Perry's Chemical Engineers' 
Handbook, 7th edition, page 2-182.</p>
</html>"));
  end cp_ch3oh_gas;

protected
  function h_ch3oh_gas
    input Temperature T;
    output MolarEnthalpy h;
  protected
    function generalIntegral
      input Temperature T;
      output MolarEnthalpy y;
      constant GaseousMethanolParameters p;
    algorithm
      y := (p.C1*T + 4*p.C2*p.C3^2/(2*p.C3*exp(2*p.C3/T)-2*p.C3) +
                     4*p.C4*p.C5^2/(2*p.C5*exp(2*p.C5/T)+2*p.C5)) / 1000;
    end generalIntegral;
  algorithm
    h := generalIntegral(T) - generalIntegral(298.15);
    annotation (Documentation(info="<html>
<p>This function returns the enthalpy of gaseous methanol at a given temperature.</p>
<p>The general integral used in this function is the analytic integral of the formula
reported in Perry's Chemical Engineers' Handbook, 7th edition, page 2-182, and has 
been obtained with Maxima.</p>
</html>"));
  end h_ch3oh_gas;

protected
  partial function AntoineLaw
    input Temperature T;
    output PartialPressure p;
  protected
    parameter Real A;
    parameter Real B;
    parameter Real C;
  algorithm
    // This is the Antoine Law itself. Note the conversion bar->Pa.
    p := 10.0^(A - B/(T + C))*1e5;
    annotation(derivative=AntoineLawDerivative);
    annotation (Documentation(info="<html>
<p>This abstract function implements Antoine's law for vapour pressure.</p>
</html>"));
  end AntoineLaw;

protected
  partial function AntoineLawDerivative
    import Modelica.Math.log;
    input Temperature T;
    input Real der_T;
    output Real der_p;
  protected
    constant Real k = log(10);
    parameter Real A;
    parameter Real B;
    parameter Real C;
  algorithm
    // This is Antoine Law's derivative with respect to time. Note the conversion bar->Pa.
    der_p := k*B*10.0^(A - B/(T + C))/(T+C)/(T+C)*1e5 * der_T;

  end AntoineLawDerivative;

protected
  function p_h2o
    extends AntoineLaw(
      T(min=255.8, max=373.0),
      A=4.65430,
      B=1435.264,
      C=-64.848);
    annotation(derivative=dp_h2o_dt);
    annotation (Documentation(info="<html>
<p>This function returns the vapour pressure of water given the temperature. Data from
<a href=\"http://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Mask=4#Thermo-Phase\">NIST</a>
(data taken from Stull, 1947).</p>
</html>"));
  end p_h2o;

public
  function dp_h2o_dt
    extends AntoineLawDerivative(
      T(min=255.8, max=373.0),
      A=4.65430,
      B=1435.264,
      C=-64.848);
  end dp_h2o_dt;

protected
  function p_ch3oh
    extends AntoineLaw(
      T(min=288.0, max=356.83),
      A=5.20409,
      B=1581.341,
      C=-33.50);
    annotation(derivative=dp_ch3oh_dt);
    annotation (Documentation(info="<html>
<p>This function returns the vapour pressure of methanol given the temperature. Data from
<a href=\"http://webbook.nist.gov/cgi/cbook.cgi?ID=C67561&Mask=4#Thermo-Phase\">NIST</a>
(data taken from Ambrose and Sprake, 1970).</p>
</html>"));
  end p_ch3oh;

protected
  function dp_ch3oh_dt
    extends AntoineLawDerivative(
      T(min=288.0, max=356.83),
      A=5.20409,
      B=1581.341,
      C=-33.50);
  end dp_ch3oh_dt;

protected
  function rho_h2o
    input Temperature T;
    output Density rho;
  protected
    constant Real a = 252.37;
    constant Real b = 6.5968;
    constant Real c = -18.288E-3;
    constant Real d = 15.222E-6;
  algorithm
    rho := a + b*T + c*T*T + d*T*T*T;
    annotation (derivative=drho_h2o_dt, Documentation(info="<html>
<p>This function returns the density of water in liquid phase at one 
standard atmosphere of pressure; the data has been interpolated with 
a cubic function, and deviates from data provided by NIST by at most
0.25 kg/m&sup3;.</p>
</html>"));
  end rho_h2o;

protected
  function drho_h2o_dt
    input Temperature T;
    input Real der_T;
    output Real der_rho;
  protected
    constant Real b = 6.5968;
    constant Real c = -18.288E-3;
    constant Real d = 15.222E-6;
  algorithm
    der_rho := (b + 2*c*T + 3*d*T*T)*der_T;
  end drho_h2o_dt;

protected
  function rho_ch3oh
    input Temperature T;
    output Density rho;
  protected
    constant Real a = 1069.1;
    constant Real b = -0.9488;
  algorithm
    rho := a + b*T;
    annotation (derivative=drho_ch3oh_dt, Documentation(info="<html>
<p>This function returns the density of methanol in liquid phase at 
one standard atmosphere of pressure; the data has been interpolated 
with a linear function, and deviates from data provided by NIST by at 
most 0.5 kg/m&sup3;.</p>
</html>"));
  end rho_ch3oh;

protected
  function drho_ch3oh_dt
    input Temperature T;
    input Real der_T;
    output Real der_rho;
  protected
    constant Real b = -0.9488;
  algorithm
    der_rho := b*der_T;
  end drho_ch3oh_dt;

public
  function mw
    input Species s;
    output MolarMass m;
  algorithm
    if s == Species.Methanol then
      m := 32.0419e-3;
    elseif s == Species.Water then
      m := 18.0153e-3;
    elseif s == Species.Oxygen then
      m := 31.9988e-3;
    elseif s == Species.CarbonDioxide then
      m := 44.0095e-3;
    elseif s == Species.Nitrogen then
      m := 28.01348e-3;
    else
      assert(false, "==> Bad input data: requested species "+String(s)+" is not implemented.");
    end if;
    annotation (Documentation(info="<html>
<p>This function returns the molecular weight of the chemical species. Data is from NIST, for:
<ul>
<li><a href=\"http://webbook.nist.gov/cgi/cbook.cgi?ID=C67561&Units=SI\">Methanol</a>;</li>
<li><a href=\"http://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Units=SI\">Water</a>;</li>
<li><a href=\"http://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447&Units=SI\">Oxygen</a>;</li>
<li><a href=\"http://webbook.nist.gov/cgi/cbook.cgi?ID=C124389&Units=SI\">Carbon dioxide</a>;</li>
<li><a href=\"http://webbook.nist.gov/cgi/cbook.cgi?ID=C7727379&Units=SI\">Nitrogen</a>.</li>
</ul></p>
</html>"));
  end mw;

protected
  function dhf
    input Species s;
    input Phases p;
    output MolarEnthalpy f;
  algorithm
    if s == Species.Methanol and p == Phases.Gas then
      f := -205000;
    elseif s == Species.Methanol and p == Phases.Liquid then
      f := -238400;
    elseif s == Species.Water and p == Phases.Gas then
      f := -241826;
    elseif s == Species.Water and p == Phases.Liquid then
      f := -285830;
    elseif s == Species.Oxygen and p == Phases.Gas then
      f := 0.0;
    elseif s == Species.CarbonDioxide and p == Phases.Gas then
      f := -393510.0;
    elseif s == Species.Nitrogen and p == Phases.Gas then
      f := 0.0;
    else
      assert(false, "==> Bad input data: "+speciesName(s)+" in "+phaseName(p)+" phase.");
    end if;
    annotation (Documentation(info="<html>
<p>This function returns the standard enthalpy of formation of the chemical species. Data is from NIST, for:
<ul>
<li>Methanol, <a href=\"http://webbook.nist.gov/cgi/cbook.cgi?ID=C67561&Units=SI&Mask=2#Thermo-Condensed\">liquid</a>
and <a href=\"http://webbook.nist.gov/cgi/cbook.cgi?ID=C67561&Units=SI&Mask=1#Thermo-Gas\">gaseous</a>;</li>
<li>Water, <a href=\"http://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Units=SI&Mask=2#Thermo-Condensed\">liquid</a> 
and <a href=\"http://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Units=SI&Mask=1#Thermo-Gas\">gaseous</a>;</li>
<li><a href=\"http://webbook.nist.gov/cgi/cbook.cgi?ID=C124389&Units=SI&Mask=1#Thermo-Gas\">Carbon dioxide</a>
(only in gas phase).</li>
</ul>
Nitrogen and oxygen are obviously zero, and are defined only in gas phase.</p>
</html>"));
  end dhf;

public
  function h
    input Temperature T;
    input Species s;
    input Phases p;
    output MolarEnthalpy H;
  protected
    constant Temperature T_ref = 298.15;
  algorithm
    if s == Species.Methanol and p == Phases.Gas then
      H := h_ch3oh_gas(T) - h_ch3oh_gas(T_ref)  + dhf(s,p);
    elseif s == Species.Methanol and p == Phases.Liquid then
      H := h_ch3oh_liq(T) - h_ch3oh_liq(T_ref) + dhf(s,p);
    elseif s == Species.Water and p == Phases.Gas then
      H := h_h2o_gas(T) - h_h2o_gas(T_ref) + dhf(s,p);
    elseif s == Species.Water and p == Phases.Liquid then
      H := ShomateEnthalpy(T, Shomate.H2O) - ShomateEnthalpy(T_ref, Shomate.H2O) + dhf(s,p);
    elseif s == Species.Oxygen and p == Phases.Gas then
      H := ShomateEnthalpy(T, Shomate.O2) - ShomateEnthalpy(T_ref, Shomate.O2);
    elseif s == Species.CarbonDioxide and p == Phases.Gas then
      H := ShomateEnthalpy(T, Shomate.CO2) - ShomateEnthalpy(T_ref, Shomate.CO2) + dhf(s,p);
    elseif s == Species.Nitrogen and p == Phases.Gas then
      H := ShomateEnthalpy(T, Shomate.N2) - ShomateEnthalpy(T_ref, Shomate.N2);
    else
      assert(false, "==> Bad input data: "+speciesName(s)+" in "+phaseName(p)+" phase.");
    end if;
    annotation(derivative=dh_dt, Documentation(info="<html>
<p>Returns the specific enthalpy of the given component at the given temperature and 
in the given phase; the reference state is always 298.15 K.</p><p>For non-elementary species,
such as water or carbon dioxide, the specific enthalpy formation for the element in 
its phase will be added.</p>
</html>"));
  end h;

protected
  function dh_dt "Derivative of the molar enthalpy function"
    input Temperature T;
    input Species s;
    input Phases p;
    input Real der_T "The rate of variation of temperature.";
    output Real der_h;
  algorithm
    der_h := cp(T, s, p)*der_T;
  end dh_dt;

public
  function cp
    input Temperature T;
    input Species s;
    input Phases p;
    output MolarHeatCapacity CP;
  algorithm
    if s == Species.Methanol and p == Phases.Gas then
      CP := cp_ch3oh_gas(T);
    elseif s == Species.Methanol and p == Phases.Liquid then
      CP := cp_ch3oh_liq(T);
    elseif s == Species.Water and p == Phases.Gas then
      CP := cp_h2o_gas(T);
    elseif s == Species.Water and p == Phases.Liquid then
      CP := ShomateHeatCapacity(T, Shomate.H2O);
    elseif s == Species.Oxygen and p == Phases.Gas then
      CP := ShomateHeatCapacity(T, Shomate.O2);
    elseif s == Species.CarbonDioxide and p == Phases.Gas then
      CP := ShomateHeatCapacity(T, Shomate.CO2);
    elseif s == Species.Nitrogen and p == Phases.Gas then
      CP := ShomateHeatCapacity(T, Shomate.N2);
    else
      assert(false, "==> Bad input data: "+speciesName(s)+" in "+phaseName(p)+" phase.");
    end if;
    annotation (Documentation(info="<html>
<p>Returns the molar heat capacity of the given component at the given 
temperature and in the given phase.</p>
</html>"));
  end cp;

public
  function p_vap
    input Temperature T;
    input Species s;
    output PartialPressure p;
  algorithm
    if s == Species.Methanol then
      p := p_ch3oh(T);
    elseif s == Species.Water then
      p := p_h2o(T);
    else
      assert(false, "==> Vapour pressure of species "+speciesName(s)+" requested.");
    end if;
    annotation(derivative=dp_vap_dt);
    annotation (Documentation(info="<html>
<p>Returns the vapour pressure of the given component at the given temperature.</p>
</html>"));
  end p_vap;

protected
  function dp_vap_dt
    input Temperature T;
    input Species s;
    input Real der_T;
    output Real der_p;
  algorithm
    if s == Species.Methanol then
      der_p := dp_ch3oh_dt(T,der_T);
    elseif s == Species.Water then
      der_p := dp_h2o_dt(T,der_T);
    else
      assert(false, "==> Error: vapour Pressure of species "+speciesName(s)+" requested.");
    end if;

  end dp_vap_dt;

public
  function K
    input Temperature T;
    input Species s;
    output Real Kvalue
      "The ratio y/x between liquid and gaseous molar fraction.";
  protected
    constant Pressure p_env = 101325 "The environment pressure.";
    constant Real gamma_ch3oh = 2.15
      "Methanol activity coefficient (diluted acqueous solution)";
  algorithm
    if s == Species.Methanol then
      Kvalue := gamma_ch3oh * p_ch3oh(T)/p_env;
    elseif s == Species.Water then
      Kvalue := p_h2o(T)/p_env;
    else
      assert(false, "==> Equilibrium constant of species "+speciesName(s)+" requested.");
    end if;

    annotation(derivative=dK_dt, Documentation(info="<html>
<p>This function provides the chemical-equilibrium y/x ratio of liquid versus gaseous
molar fraction.</p>
<p>Note that this function assumes environmental pressure (101325 Pa).</p>
</html>"));
    annotation (Documentation(info="<html>
<p>This function provides the chemical-equilibrium x/y ratio of liquid versus gaseous
molar fraction.</p>
<p>Note that this is the inverse of the usual y/x ratio: the reason is that this way 
it is possible to include components that are present in gaseous phase but not in 
liquid phase, by setting their K-value to be zero.</p>
<p>Note also that this function assumes environmental pressure (101.325 Pa).</p>
</html>"));
  end K;

protected
  function dK_dt
    input Temperature T;
    input Species s;
    input Real der_T "The derivative of temperature";
    output Real der_K "The derivative of K-values with respect to temperature";
  protected
    constant Pressure p_env = 101325 "Environment pressure";
    constant Real gamma_ch3oh = 2.15;
  algorithm
    if s == Species.Methanol then
      der_K := gamma_ch3oh * dp_ch3oh_dt(T,der_T)/p_env;
    elseif s == Species.Water then
      der_K := dp_h2o_dt(T,der_T)/p_env;
    else
      der_K := 0.0;
    end if;
  end dK_dt;

public
  function rho
    input Temperature T;
    input Species s;
    input Phases p;
    output Density RHO;
    import Modelica.Constants.R;
  protected
    constant Pressure p_env=101325.0;
  algorithm
    if p == Phases.Liquid then
      if s == Species.Methanol then
        RHO := rho_ch3oh(T);
      elseif s == Species.Water then
        RHO := rho_h2o(T);
      else
        assert(false, "==> Requested liquid-phase density for non-implemented component "+speciesName(s)+".");
      end if;
    elseif p == Phases.Gas then
      // Note: Assuming ideal gas.
      RHO := mw(s)*p_env/R/T;
    else
      assert(false, "==> Requested density for unknown phase ("+phaseName(p)+").");
    end if;
    annotation (derivative=drho_dt, Documentation(info="<html>
<p>Returns the density of the given component at the given temperature and in the given phase.
For gases, the ideal gas lawo is assumed; for liquids, a linear interpolation on data is
carried out.</p>
</html>"));
  end rho;

  annotation (uses(Modelica(version="3.1"), Units(version="1")),
                                               Documentation(info="<html>
<p>This package contains a number of property functions, such as molar enthalpy,
molar heat capacity, vapour pressure, molecular weight and so on.</p>
</html>"),
    version="1",
    conversion(noneFromVersion=""));

protected
  function drho_dt "Derivative of the density function"
    import Modelica.Constants.R;

    input Temperature T;
    input Species s;
    input Phases p;
    input Real der_T "The rate of variation of temperature";
    output Real der_rho "The derivative drho/dt";

  protected
    constant Pressure p_env=101325.0;
    constant Temperature dT = 0.01;
  algorithm
    if p == Phases.Liquid then
      if s == Species.Methanol then
        der_rho := drho_ch3oh_dt(T, der_T);
      elseif s == Species.Water then
        der_rho := drho_h2o_dt(T, der_T);
      else
        assert(false, "==> Density in liquid phase of gaseous component "+speciesName(s)+" requested.");
      end if;
    elseif p == Phases.Gas then
      // NOTE Assuming ideal gas.
      der_rho := -mw(s)*p_env/R/T/T*der_T;
    else
      assert(false, "==> Requested density for unknown phase ("+phaseName(p)+").");
    end if;
  end drho_dt;

public
  function rr "The vapour fraction of a methanol-water-gas mixture."
    import Modelica.Constants.inf;

    input MoleFraction z_methanol "The methanol fraction.";
    input MoleFraction z_water "The water fraction.";
    input Temperature T
      "The temperature at which the equilibrium is calculated.";
    output MoleFraction beta
      "The vapour molar fraction, or vaporisation ratio.";
  protected
    constant Pressure p_env = 101325 "The environment pressure.";
    Real C_methanol;
    Real C_water;
    Real a;
    Real b;
    Real c;
    Real delta;
    annotation (derivative=drr_dt, Documentation(info="<html>
<p>This function solves the Rachford-Rice equation for systems with methanol, water and
gaseous components. In order to solve the system analytically as a second-degree polynomial,
some assumptions have been made:</p>
<ul>
<li>No other components than methanol and water may condense;</li>
<li>All other components are gases whose molar fraction is the complement to one of
the fractions of methanol and water.</li>
</ul>
 
<h3>Rearranging the Rachford-Rice as a polynomial</h3>
<p>Assimilating all gases into one is a nice trick to reduce a 5-component problem into a
3-component one. The latter has only two solutions to the Rachford-Rice, and they are
easily found with the formula for 2nd-degree polynomials.</p>
 
<p>The actual polynomial is:<br/>
C<sub>H<sub>2</sub>O</sub>C<sub>CH<sub>3</sub>OH</sub> &beta;&sup2; + 
[C<sub>CH<sub>3</sub>OH</sub>z<sub>CH<sub>3</sub>OH</sub> + C<sub>H<sub>2</sub>O</sub>z<sub>H<sub>2</sub>O</sub>
+ (C<sub>CH<sub>3</sub>OH</sub> + C<sub>H<sub>2</sub>O</sub>)z<sub>g</sub>] &beta; + z<sub>g</sub> = 0</p>
<p>For a simpler notation, C<sub>i</sub> = K<sub>i</sub> - 1 has been used.</p>
 
 
<h3>Choice of the root of h(&beta;)</h3>
 
<p>The criterion to obtain an analytical solution to the Rachford-Rice relation can change
according to which components we actually have, and the temperature range in which we are. In
particular, when passing the boiling temperature of methanol, one of the asymptote of the
Rachford-Rice equation h(&beta;) moves from the right-hand to the left-hand plane, meaning
we have to pick a different solution when we pass beyond this temperature.</p>
 
<p>There are two cases:</p>
<ul>
<li>If C<sub>H<sub>2</sub>O</sub> is smaller than zero, yet C<sub>CH<sub>3</sub>OH</sub> is larger, we need to select the <em>second</em> or larger solution.</li>
<li>Otherwise, we select the <em>first</em> or smaller solution.</li>
</ul>
<p>Notice that the two cases of higher and lower solution have <em>the same mathematical expression</em>, since the coefficient of &beta;&sup2; changes sign in correspondance of the change in root selection; therefore, we do not need to change the sign in front of the root of &Delta;.</p>
 
<h3>References</h3>
<p>Whitson, Curtis H., and Michelsen, Michael L.: <em>The negative flash</em>, Fluid Phase Equilibria 53, 51&ndash;71, 1989.</p>
<p>Warren, John H., and Adewumi, Michael A.: <em>Polynomial Objective Functions for Flash Calculations: Binary, Ternary, and 
Quaternary Systems</em>, Industrial and Engineering Chemistry Research 32(7), 1528&ndash;1530, July 1993.</p>
</html>"));
  algorithm
    C_methanol := K(T, Species.Methanol) - 1;
    C_water    := K(T, Species.Water) - 1;

    assert( C_methanol > C_water, "==> Water is more volatile than methanol at temperature T="+String(T)+" K: this should not be possible.");

    a := C_methanol*C_water;
    b := C_methanol*z_methanol + C_water*z_water + (C_methanol + C_water)*(1.0 - z_methanol - z_water);
    c := 1.0 - z_methanol - z_water;
    delta := b*b - 4*a*c;

    if C_methanol == 0 then // Parameter a would be zero and would cause a singularity
      beta := - c / b;  // Analytic limit for C_methanol -> 0
    else
      beta := (-b-sqrt(delta))/(2*a);
    end if;

  end rr;

protected
  function drr_dt "The vapour fraction of a methanol-water-gas mixture."
    input MoleFraction z_methanol "The methanol fraction.";
    input MoleFraction z_water "The water fraction.";
    input Temperature T
      "The temperature at which the equilibrium is calculated.";
    input Real der_z_methanol "The time derivative of the methanol fraction.";
    input Real der_z_water "The time derivative of the water fraction.";
    input Real der_T "The time derivative of temperature.";
    output Real der_beta "The derivative of the vapour molar fraction.";

  protected
    constant Pressure p_env = 101325 "The environment pressure.";
    Real C_methanol;
    Real C_water;

    Real a;
    Real b;
    Real c;
    Real delta;
    Real beta;

    Real dCm_dt;
    Real dCw_dt;
    Real da_dCm;
    Real da_dCw;
    Real da_dt;

    Real db_dCm;
    Real db_dCw;
    Real db_dzm;
    Real db_dzw;
    Real db_dt;

    Real dc_dzm;
    Real dc_dzw;
    Real dc_dt;

    Real dbeta_da;
    Real dbeta_db;
    Real dbeta_dc;

  algorithm
    C_methanol := K(T, Species.Methanol) - 1;
    C_water    := K(T, Species.Water) - 1;

    a := C_methanol*C_water;
    b := C_methanol*z_methanol+C_water*z_water+(C_methanol+C_water)*(1.0-z_methanol-z_water);
    c := 1.0-z_methanol-z_water;
    delta := b*b - 4*a*c;
    beta := (-b-sqrt(delta))/(2*a);

    dCm_dt := dK_dt(T, Species.Methanol, der_T);
    dCw_dt := dK_dt(T, Species.Water, der_T);

    da_dCm := C_water;
    da_dCw := C_methanol;
    da_dt  := da_dCm*dCm_dt + da_dCw*dCw_dt;

    db_dCm := 1-z_water;
    db_dCw := 1-z_methanol;
    db_dzm := -C_water;
    db_dzw := -C_methanol;
    db_dt  := db_dCm*dCm_dt + db_dCw*dCw_dt + db_dzm*der_z_methanol + db_dzw*der_z_water;

    dc_dzm := -1;
    dc_dzw := -1;
    dc_dt  := dc_dzm*der_z_methanol + dc_dzw*der_z_water;

    dbeta_da := c/a/sqrt(delta) + (b+sqrt(delta))/2/a/a;
    dbeta_db := -(1+b/sqrt(delta))/2/a;
    dbeta_dc := 1/sqrt(delta);

    der_beta := dbeta_da*da_dt + dbeta_db*db_dt + dbeta_dc*dc_dt;

  end drr_dt;

public
  function pDoverT175 "Non-varying part of diffusivity coefficients"
    input Species i "First species";
    input Species j "Second species";
    output Real pDT "Constant part of the diffusivity coefficient";

  algorithm
    pDT := 0.1013 * sqrt(1E-3/mw(i)+1E-3/mw(j)) / (diffVol(i)^(1/3)+diffVol(j)^(1/3))^2;

    annotation (Documentation(info="<html>
<p>This function provides the constant part of a gas-phase diffusivity coefficient
estimated with Fuller's rule. The full-fledged diffusivity coefficient is obtained 
by multiplying by temperature in kelvin raised to 1.75 and dividing by pressure 
in pascal.</p>
<p>The order in which the components are provided is irrelevant.</p>
<p>Data obtained by Perry's Chemical Engineers' Handbook, 7<sup>th</sup> edition,
page 2-370.</p>
</html>"));
  end pDoverT175;

protected
  function diffVol "Molecular diffusion volumes"
    input Species s;
    output Real v;

  algorithm
    if s == Species.Methanol then
      v := 29.901;
    elseif s == Species.Water then
      v := 12.7;
    elseif s == Species.Oxygen then
      v := 16.6;
    elseif s == Species.CarbonDioxide then
      v := 26.9;
    elseif s == Species.Nitrogen then
      v := 17.9;
    else
      assert( false, "==> Unhandled component '"+speciesName(s)+"' requested.");
    end if;

    annotation (Documentation(info="<html>
<p>This simple function returns the molecular diffusion volumes of some
species. Data is taken from page 2-370 of Perry's Chemical Engineers' 
Handbook, 7<sup>th</sup> edition, table 2-400.</p>
</html>"));
  end diffVol;

public
  function D "Non-varying part of diffusivity coefficients"
    input Temperature T "Temperature";
    input Pressure p "Enviroment pressure";
    input Species i "First species";
    input Species j "Second species";
    output DiffusionCoefficient Dij
      "Constant part of the diffusivity coefficient";

  algorithm
    Dij := pDoverT175(i,j) * T^1.75 / p;

    annotation (Documentation(info="<html>
<p>This function provides the a gas-phase diffusivity coefficient estimated with Fuller's 
rule, given two components between which to calculate it.</p>
<p>The order in which the components are provided is irrelevant.</p>
<p>Data obtained by Perry's Chemical Engineers' Handbook, 7<sup>th</sup> edition,
page 2-370.</p>
</html>"));
  end D;

  model Air "An object calculating and providing air composition"
    extends Modelica.Icons.Record;

    import Units.RelativeHumidity;
    import Units.Temperature;
    import Units.MolarEnthalpy;
    import Modelica.SIunits.MoleFraction;

    outer RelativeHumidity RH_env "Environment relative humidity";
    outer Temperature T_env "Environment temperature";

    MoleFraction[Species] y "Molar fractions";
    MolarEnthalpy H "Molar enthalpy";

  equation
    sum(y) = 1;

    y[Species.Oxygen] / 0.21 = y[Species.Nitrogen] / 0.79;

    y[Species.CarbonDioxide] = 385E-6*sum(y[Incondensables]);

    y[Species.Methanol] = 0;

    y[Species.Water] = RH_env/100 * K(T_env, Species.Water);

    H = sum( h(T_env, i, Phases.Gas) * y[i] for i in Species);

    annotation (Documentation(info="<html>
<p>This object contains a vector of molar fractions in environment air and
the corresponding molar enthalpy of air; the values are adapted to current
values of <tt>T_env</tt> and <tt>RH_env</tt>.</p>
</html>"));
  end Air;

public
  package Test

  model Test_h "test case for specific enthalpy"
    import Units.Temperature;
    import Units.MolarEnthalpy;
    import Modelica.SIunits.Time;
    import Thermo.h;

    parameter Time blinkOfAnEye = 1E-6;

    parameter Temperature T_0 = 273.15 "Initial temperature";
    parameter Temperature T_f = 330 "Final temperature";
    Temperature T = T_0+time/second*(T_f-T_0)
        "Varying temperature,linear with time";

    Real[7] error;

    protected
    MolarEnthalpy enthalpy[7];
    MolarEnthalpy old_enthalpy[7];
    constant Time second = 1 "To get adimensional time";

  equation
    for i in Species loop
      enthalpy[i] = h(T, i, Phases.Gas);
    end for;
    enthalpy[6] = h(T, Species.Methanol, Phases.Liquid);
    enthalpy[7] = h(T, Species.Water, Phases.Liquid);

    for i in 1:7 loop
      old_enthalpy[i] = delay(enthalpy[i], blinkOfAnEye);
      error[i] = abs((enthalpy[i]-old_enthalpy[i])/blinkOfAnEye - der(enthalpy[i]));
    end for;

      annotation (experiment, experimentSetupOutput);

  end Test_h;

  model Test_cp "test case for specific enthalpy"
    import Units.Temperature;
    import Units.MolarEnthalpy;
    import Modelica.SIunits.MolarHeatCapacity;
    import Modelica.SIunits.Time;
    import Thermo.cp;
    import Thermo.h;

    parameter Time blinkOfAnEye = 1E-6;

    parameter Temperature T_0 = 273.15 "Initial temperature";
    parameter Temperature T_f = 330 "Final temperature";
    Temperature T = T_0+time/second*(T_f-T_0)
        "Varying temperature,linear with time";

    MolarHeatCapacity[7] error;

    protected
    MolarHeatCapacity heatCapacity[7];
    MolarEnthalpy enthalpy[7];
    constant Time second = 1 "To get adimensional time";

  equation
    for i in Species loop
      enthalpy[i] = h(T, i, Phases.Gas);
      heatCapacity[i] = cp(T, i, Phases.Gas);
    end for;
    enthalpy[6] = h(T, Species.Methanol, Phases.Liquid);
    heatCapacity[6] = cp(T, Species.Methanol, Phases.Liquid);
    enthalpy[7] = h(T, Species.Water, Phases.Liquid);
    heatCapacity[7] = cp(T, Species.Water, Phases.Liquid);

    for i in 1:7 loop
      error[i] = abs(der(enthalpy[i]) - heatCapacity[i]*der(T));
    end for;

      annotation (experiment, experimentSetupOutput);

  end Test_cp;

  model Test_p_vap "test case for species equilibrum constant"
    import Units.Temperature;
    import Modelica.SIunits.PartialPressure;
    import Modelica.SIunits.Time;

    parameter Time blinkOfAnEye = 1E-6;

    parameter Temperature T_0 = 273.15 "Initial temperature";
    parameter Temperature T_f = 330 "Final temperature";
    Temperature T = T_0+time/second*(T_f-T_0)
        "Varying temperature,linear with time";

    Real error_meoh = abs((p_ch3oh - old_p_ch3oh)/blinkOfAnEye - der(p_ch3oh));
    Real error_h2o =  abs((p_h2o - old_p_h2o)/blinkOfAnEye - der(p_h2o));

    protected
    PartialPressure p_h2o =   p_vap(T, Species.Water);
    PartialPressure p_ch3oh = p_vap(T, Species.Methanol);
    PartialPressure old_p_h2o =   delay(p_h2o, blinkOfAnEye);
    PartialPressure old_p_ch3oh = delay(p_ch3oh, blinkOfAnEye);
    constant Time second = 1 "To get adimensional time";

  end Test_p_vap;

  model Test_K "test case for species equilibrum constant"
    import Units.Temperature;
    import Modelica.SIunits.Time;

    parameter Time blinkOfAnEye = 1E-6;

    parameter Temperature T_0 = 273.15 "Initial temperature";
    parameter Temperature T_f = 330 "Final temperature";
    Temperature T = T_0+time/second*(T_f-T_0)
        "Varying temperature,linear with time";

    Real error_meoh = abs((K_ch3oh - old_K_ch3oh)/blinkOfAnEye - der(K_ch3oh));
    Real error_h2o =  abs((K_h2o - old_K_h2o)/blinkOfAnEye - der(K_h2o));

    protected
    Real K_h2o =   K(T, Species.Water);
    Real K_ch3oh = K(T, Species.Methanol);
    Real old_K_h2o =   delay(K_h2o, blinkOfAnEye);
    Real old_K_ch3oh = delay(K_ch3oh, blinkOfAnEye);
    constant Time second = 1 "To get adimensional time";

  end Test_K;

  model Test_rho "test case for specific enthalpy"
    import Units.Temperature;
    import Modelica.SIunits.Time;
    import Modelica.SIunits.Density;
    import Thermo.rho;

    parameter Time blinkOfAnEye = 1E-6;

    parameter Temperature T_0 = 273.15 "Initial temperature";
    parameter Temperature T_f = 330 "Final temperature";
    Temperature T = T_0+time/second*(T_f-T_0)
        "Varying temperature,linear with time";

    Real[7] error;

    protected
    MolarEnthalpy density[7];
    MolarEnthalpy old_density[7];
    constant Time second = 1 "To get adimensional time";

  equation
    for i in Species loop
      density[i] = rho(T, i, Phases.Gas);
    end for;
    density[6] = rho(T, Species.Methanol, Phases.Liquid);
    density[7] = rho(T, Species.Water, Phases.Liquid);

    for i in 1:7 loop
      old_density[i] = delay(density[i], blinkOfAnEye);
      error[i] = abs((density[i]-old_density[i])/blinkOfAnEye - der(density[i]));
    end for;

      annotation (experiment, experimentSetupOutput);

  end Test_rho;

    model TestRR
      import Units.Temperature;
      import Modelica.SIunits.MoleFraction;
      import Modelica.SIunits.Time;

      parameter Time blinkOfAnEye = 1E-6;
      parameter Temperature T_0 = 273.15;
      parameter Temperature T_f = 330;
      parameter MoleFraction zm_0 = 0;
      parameter MoleFraction zm_f = 0.5;
      parameter MoleFraction zw_0 = 0.;
      parameter MoleFraction zw_f = 0.3;

      Temperature T = T_0 + (T_f - T_0)*time/second;
      MoleFraction zm = zm_0 + (zm_f - zm_0)*time/second;
      MoleFraction zw = zw_0 + (zw_f - zw_0)*time/second;

      MoleFraction beta = rr(zm, zw, T);

      Real error = abs( (beta-old_beta)/blinkOfAnEye - der(beta));

    protected
      MoleFraction old_beta = delay(beta, blinkOfAnEye);
      constant Time second = 1 "To get adimensional time";

      annotation (experiment(NumberOfIntervals=5000),
                              experimentSetupOutput);
    end TestRR;

    model TestAir

      import Units.RelativeHumidity;
      import Modelica.SIunits.Temperature;
      import Modelica.SIunits.Time;

      parameter RelativeHumidity RH_env_start = 30 "initial relative humidity";
      parameter RelativeHumidity RH_env_stop = 60 "final relative humidity";
      parameter Temperature T_env_start = 275 "initial temperature";
      parameter Temperature T_env_stop = 350 "final temperature";

      inner Temperature T_env "Enviroment temperature";
      inner RelativeHumidity RH_env "Environment relative humidity";

      Air air;

    protected
      constant Time second = 1 "To get adimensional time";

    equation
      T_env = T_env_start + (T_env_stop-T_env_start)*time/second;
      RH_env = RH_env_start + (RH_env_stop-RH_env_start)*time/second;

    end TestAir;
  end Test;
end Thermo;
