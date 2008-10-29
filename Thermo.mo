encapsulated package Thermo "Thermodynamic library" 
  
  import Modelica.SIunits.Temperature;
  import Modelica.SIunits.Pressure;
  import Modelica.SIunits.MolarInternalEnergy;
  import Modelica.SIunits.MolarHeatCapacity;
  import Modelica.SIunits.MolarEntropy;
  import Modelica.SIunits.PartialPressure;
  import Modelica.SIunits.Density;
  import Modelica.SIunits.MolarMass;
  import Modelica.SIunits.MoleFraction;
  
public 
  type MolarEnthalpy = MolarInternalEnergy annotation (Documentation(info="<html>
<p>Definition lacking from Modelica library.</p>
</html>"));
  
  // Species
public 
  constant Integer[:] AllSpecies =    1:5;
  constant Integer[:] GasSpecies =    3:5;
  constant Integer[:] LiquidSpecies = 1:2;
  constant Integer Methanol =      1;
  constant Integer Water =         2;
  constant Integer Oxygen =        3;
  constant Integer CarbonDioxide = 4;
  constant Integer Nitrogen =      5;
  
  // Phases
public 
  constant Integer[:] AllPhases = 1001:1002;
  constant Integer GasPhase =     1001;
  constant Integer LiquidPhase =  1002;
  
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
<p>Note: parameter D has been added to those from Perry, so that enthalpy would be 0 at 25°C.</p>
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
<p>Note: parameter D has been added to those from Perry, so that enthalpy would be 0 at 25°C.</p>
</html>"));
  end LiquidMethanolParameters;
  
protected 
  record ShomateParameters "Set of parameters for Shomate equations" 
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
<p>This structure houses the parameters necessary for Shomate equations, used to
calculate the molar enthalpy, the molar heat capacity and the molar entropy of many
components. Data are from NIST, in detail:
<ul>
<li><a href=\"http://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Units=SI&Mask=2#Thermo-Condensed\">Liquid water</a>;</li>
<li><a href=\"http://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447&Units=SI&Mask=1#Thermo-Gas\">Oxygen</a>;</li>
<li><a href=\"http://webbook.nist.gov/cgi/cbook.cgi?ID=C124389&Units=SI&Mask=1#Thermo-Gas\">Carbon dioxide</a>;</li>
<li><a href=\"http://webbook.nist.gov/cgi/cbook.cgi?ID=C7727379&Units=SI&Mask=1#Thermo-Gas\">Nitrogen</a>.</li>
</ul>
</p>
</html>"));
  end ShomateParameters;
  
protected 
  constant ShomateParameters ShomateO2( A=29.65900, B=6.137261, C=-1.186521,
                                        D=0.095780, E=-0.219663, F=-9.861391,
                                        G=237.9480, H=0.000000,
                                        T_min=298.0, T_max=6000.0);
  
protected 
  constant ShomateParameters ShomateN2( A=26.09200, B=8.218801, C=-1.976141,
                                        D=0.159274, E=0.044434, F=-7.989230,
                                        G=221.0200, H=0.000000,
                                        T_min=298.0, T_max=6000.0);
protected 
  constant ShomateParameters ShomateCO2( A=24.99735, B=55.18696, C=-33.69137,
                                         D=7.948387, E=-0.136638, F=-403.6075,
                                         G=228.2431, H=-393.5224,
                                         T_min=298.0, T_max=1200.0);
protected 
  constant ShomateParameters ShomateH2O( A=-203.6060, B=1523.290, C=-3196.413,
                                         D=2474.455, E=3.855326, F=-256.5478,
                                         G=-488.7163, H=-285.8304,
                                         T_min=298.0, T_max=500.0);
public 
  function speciesName 
    input Integer n "The species code.";
    output String st "The species name.";
  algorithm 
    if n == Methanol then
      st := "methanol";
    elseif n == Water then
      st := "water";
    elseif n == Oxygen then
      st := "oxygen";
    elseif n == CarbonDioxide then
      st := "carbon dioxide";
    elseif n == Nitrogen then
      st := "nitrogen";
    else
      st := "unknown ("+String(n)+")";
    end if;
    
    annotation (Documentation(info="<html>
<p>This function returns a string with the name of the species whose
code was set into the input.</p>
</html>"));
  end speciesName;
  
public 
  function phaseName 
    input Integer n "The phase code.";
    output String st "The phase name.";
  algorithm 
    if n == GasPhase then
      st := "gas";
    elseif n == LiquidPhase then
      st := "liquid";
    else
      st := "unknown ("+String(n)+")";
    end if;
    
    annotation (Documentation(info="<html>
<p>This function returns a string with the name of the phase whose
code was set into the input.</p>
</html>"));
  end phaseName;
  
protected 
  function ShomateEnthalpy 
    input Temperature T;
    input ShomateParameters p;
    output MolarEnthalpy h;
  protected 
    Real t;
  algorithm 
    t := T/1000.0;
    h := (p.A*t + (p.B*t^2)/2 + (p.C*t^3)/3 + (p.D*t^4)/4 - p.E/t + p.F - p.H)*
      1000;
    annotation (Documentation(info="<html>
<p>This abstract function uses the Shomate parameters to calculate the enthalpy.
Note that the Shomate parameters are still unspecified.</p>
</html>"));
  end ShomateEnthalpy;
  
protected 
  function ShomateHeatCapacity 
    input Temperature T;
    input ShomateParameters p;
    output MolarHeatCapacity cp;
  protected 
    Real t;
  algorithm 
    t := T/1000.0;
    cp := p.A + p.B*t + p.C*t^2 + p.D*t^3 + p.E/t^2;
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
    h := (p.A*T + (p.B*T^2)/2 + (p.C*T^3)/3)*4.184 - p.D;
  end h_h2o_gas;
  
protected 
  function cp_h2o_gas "The specific heat of gaseous water." 
    input Temperature T;
    output MolarHeatCapacity cp;
  protected 
    constant WaterVapourParameters p;
  algorithm 
    // NOTE conversion from calories.
    cp := (p.A + p.B*T + p.C*T^2)*4.184;
  end cp_h2o_gas;
  
protected 
  function h_ch3oh_liq "The enthalpy of liquid methanol." 
    input Temperature T;
    output MolarEnthalpy h;
  protected 
    constant LiquidMethanolParameters p;
  algorithm 
    h := (p.A*T + (p.B*T^2)/2 + (p.C*T^3)/3)/1000 - p.D;
  end h_ch3oh_liq;
  
protected 
  function cp_ch3oh_liq "The specific heat of liquid methanol." 
    input Temperature T;
    output MolarHeatCapacity cp;
  protected 
    constant LiquidMethanolParameters p;
  algorithm 
    cp := (p.A + p.B*T + p.C*T^2)/1000;
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
    annotation (Documentation(info="<html>
<p>This abstract function implements Antoine's law for vapour pressure.</p>
</html>"));
  end AntoineLaw;
  
protected 
  partial function AntoineLawDerivative 
    input Temperature T;
    output Real der_p;
  protected 
    parameter Real A;
    parameter Real B;
    parameter Real C;
  algorithm 
    // This is Antoine Law's derivative with respect to temperature. Note the conversion bar->Pa.
    der_p := (B*10.0^(A - B/(T + C)))/(T+C)^2*1e5;
    
  end AntoineLawDerivative;
  
protected 
  function p_h2o 
    extends AntoineLaw(
      T(min=255.8, max=373.0),
      A=4.65430,
      B=1435.264,
      C=-64.848);
    annotation (Documentation(info="<html>
<p>This function returns the vapour pressure of water given the temperature. Data from
<a href=\"http://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Mask=4#Thermo-Phase\">NIST</a>
(data taken from Stull, 1947).</p>
</html>"));
  end p_h2o;
  
protected 
  function dp_h2o_dT 
    extends AntoineLawDerivative(
      T(min=255.8, max=373.0),
      A=4.65430,
      B=1435.264,
      C=-64.848);
  end dp_h2o_dT;
  
protected 
  function p_ch3oh 
    extends AntoineLaw(
      T(min=288.0, max=356.83),
      A=5.20409,
      B=1581.341,
      C=-33.50);
    annotation (Documentation(info="<html>
<p>This function returns the vapour pressure of methanol given the temperature. Data from
<a href=\"http://webbook.nist.gov/cgi/cbook.cgi?ID=C67561&Mask=4#Thermo-Phase\">NIST</a>
(data taken from Ambrose and Sprake, 1970).</p>
</html>"));
  end p_ch3oh;
  
protected 
  function dp_ch3oh_dT 
    extends AntoineLawDerivative(
      T(min=288.0, max=356.83),
      A=5.20409,
      B=1581.341,
      C=-33.50);
  end dp_ch3oh_dT;
  
protected 
  function LinearInterpolation 
    input Real x;
    input Real[:] dataX;
    input Real[:] dataY;
    output Real y;
  protected 
    Integer pos;
    Real f;
  algorithm 
    if x < dataX[1] then
      y := dataY[1];
    elseif x >= dataX[end] then
      y := dataY[end];
    else
      pos := 1;
      while x > dataX[pos + 1] loop
        pos := pos + 1;
      end while;
      f := (x - dataX[pos])/(dataX[pos + 1] - dataX[pos]);
      y := dataY[pos] + f*(dataY[pos + 1] - dataY[pos]);
    end if;
    
    annotation (Documentation(info="<html>
<p>This function provides a linear interpolation of an indipendent variable <em>x</em>
to a dependent variable <em>y</em>, provided two data sets <em>X</em> and <em>Y</em>.
When extrapolating, the extreme value of the <em>Y</em> data set is taken.</p>
</html>"));
  end LinearInterpolation;
  
protected 
  function rho_h2o 
    input Temperature T;
    output Density rho;
  protected 
    Real dataT[:]={273.16,278.16,283.16,288.16,293.16,298.16,303.16,308.16,
        313.16,318.16,323.16,328.16,333.16,338.16,343.16,348.16,353.16,358.16,
        363.16,368.16,373.12};
    Real dataD[:]={999.84,999.97,999.70,999.10,998.21,997.05,995.65,994.03,
        992.21,990.21,988.03,985.69,983.19,980.55,977.76,974.84,971.78,968.60,
        965.30,961.88,958.37};
  algorithm 
    rho := LinearInterpolation(T,dataT,dataD);
    annotation (Documentation(info="<html>
<p>This function returns the density of water in liquid phase at one 
standard atmosphere of pressure.</p>
</html>"));
  end rho_h2o;
  
protected 
  function rho_ch3oh 
    input Temperature T;
    output Density rho;
  protected 
    Real dataT[:]={273.15,278.15,283.15,288.15,293.15,298.15,303.15,308.15,
        313.15,318.15,323.15,328.15,333.15,337.63};
    Real dataD[:]={809.73,805.05,800.37,795.69,791.01,786.33,781.63,776.91,
        772.17,767.39,762.58,757.72,752.81,748.36};
  algorithm 
    rho := LinearInterpolation(T,dataT,dataD);
    annotation (Documentation(info="<html>
<p>This function returns the density of methanol in liquid phase at 
one standard atmosphere of pressure.</p>
</html>"));
  end rho_ch3oh;
  
public 
  function mw 
    input Integer n "Component";
    output MolarMass m;
  algorithm 
    if n == Methanol then
      m := 32.0419e-3;
    elseif n == Water then
      m := 18.0153e-3;
    elseif n == Oxygen then
      m := 31.9988e-3;
    elseif n == CarbonDioxide then
      m := 44.0095e-3;
    elseif n == Nitrogen then
      m := 28.01348e-3;
    else
      assert(false, "Bad input data: requested species "+String(n)+" is not implemented.");
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
  
public 
  function dhf 
    input Integer n "Component";
    input Integer p "Phase";
    output MolarEnthalpy f;
  algorithm 
    if n == Methanol and p == GasPhase then
      f := -205000;
    elseif n == Methanol and p == LiquidPhase then
      f := -238400;
    elseif n == Water and p == GasPhase then
      f := -241826;
    elseif n == Water and p == LiquidPhase then
      f := -285830;
    elseif n == Oxygen and p == GasPhase then
      f := 0.0;
    elseif n == CarbonDioxide and p == GasPhase then
      f := -393510.0;
    elseif n == Nitrogen and p == GasPhase then
      f := 0.0;
    else
      assert(false, "Bad input data: "+speciesName(n)+" in "+phaseName(p)+" phase.");
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
    input Integer n "Component";
    input Integer p "Phase";
    output MolarEnthalpy H;
  protected 
    constant Temperature T_ref = 298.15;
  algorithm 
    if n == Methanol and p == GasPhase then
      H := h_ch3oh_gas(T) - h_ch3oh_gas(T_ref) + dhf(Methanol, GasPhase) - dhf(Methanol, LiquidPhase);
    elseif n == Methanol and p == LiquidPhase then
      H := h_ch3oh_liq(T) - h_ch3oh_liq(T_ref);
    elseif n == Water and p == GasPhase then
      H := h_h2o_gas(T) - h_h2o_gas(T_ref) + dhf(Water, GasPhase) - dhf(Water, LiquidPhase);
    elseif n == Water and p == LiquidPhase then
      H := ShomateEnthalpy(T, ShomateH2O) - ShomateEnthalpy(T_ref, ShomateH2O);
    elseif n == Oxygen and p == GasPhase then
      H := ShomateEnthalpy(T, ShomateO2) - ShomateEnthalpy(T_ref, ShomateO2);
    elseif n == CarbonDioxide and p == GasPhase then
      H := ShomateEnthalpy(T, ShomateCO2) - ShomateEnthalpy(T_ref, ShomateCO2);
    elseif n == Nitrogen and p == GasPhase then
      H := ShomateEnthalpy(T, ShomateN2) - ShomateEnthalpy(T_ref, ShomateN2);
    else
      assert(false, "Bad input data: "+speciesName(n)+" in "+phaseName(p)+" phase.");
    end if;
    annotation(derivative=dh_dT, Documentation(info="<html>
<p>Returns the enthalpy of the given component at the given temperature and in
the given phase; the reference state is always 298.15 K.</p>
<p>Water and methanol are assumed to be in liquid phase at the reference state;
when asking their enthalpy in gas phase, the evaporation enthalpy will be added
to the result.</p>
</html>"));
    annotation (Documentation(info="<html>
<p>Returns the enthalpy of the given component at the given temperature and in
the given phase; the reference state is always 298.15 K.</p>
</html>"));
  end h;
  
protected 
  function dh_dT "Derivative of the enthalpy function" 
    input Temperature T;
    input Integer n "Component";
    input Integer p "Phase";
    input Real der_T "The rate of variation of temperature.";
    output MolarHeatCapacity c;
  algorithm 
    c := cp(T,n,p);
  end dh_dT;
  
public 
  function cp 
    input Temperature T;
    input Integer n "Component";
    input Integer p "Phase";
    output MolarHeatCapacity CP;
  algorithm 
    if n == Methanol and p == GasPhase then
      CP := cp_ch3oh_gas(T);
    elseif n == Methanol and p == LiquidPhase then
      CP := cp_ch3oh_liq(T);
    elseif n == Water and p == GasPhase then
      CP := cp_h2o_gas(T);
    elseif n == Water and p == LiquidPhase then
      CP := ShomateHeatCapacity(T, ShomateH2O);
    elseif n == Oxygen and p == GasPhase then
      CP := ShomateHeatCapacity(T, ShomateO2);
    elseif n == CarbonDioxide and p == GasPhase then
      CP := ShomateHeatCapacity(T, ShomateCO2);
    elseif n == Nitrogen and p == GasPhase then
      CP := ShomateHeatCapacity(T, ShomateN2);
    else
      assert(false, "Bad input data: "+speciesName(n)+" in "+phaseName(p)+" phase.");
    end if;
    annotation (Documentation(info="<html>
<p>Returns the molar heat capacity of the given component at the given 
temperature and in the given phase.</p>
</html>"));
  end cp;
  
public 
  function p_vap 
    input Temperature T;
    input Integer n "Component";
    output PartialPressure p;
  algorithm 
    if n == Methanol then
      p := p_ch3oh(T);
    elseif n == Water then
      p := p_h2o(T);
    else
      assert(false, "Vapour pressure of species "+speciesName(n)+" requested.");
    end if;
    annotation(derivative=dp_vap_dT);
    annotation (Documentation(info="<html>
<p>Returns the vapour pressure of the given component at the given temperature.</p>
</html>"));
  end p_vap;
  
public 
  function K 
    input Temperature T "Temperature";
    input Integer n "Component";
    output Real Kvalue 
      "The ratio y/x between liquid and gaseous molar fraction.";
  protected 
    constant Pressure p_env = 101325 "The environment pressure.";
  algorithm 
    if n == Methanol then
      Kvalue := p_ch3oh(T)/p_env;
    elseif n == Water then
      Kvalue := p_h2o(T)/p_env;
    else
      Kvalue := 1e20; // TODO implement solutibilies
    end if;
    
    annotation(derivative=dK_dT);
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
  function dK_dT 
    input Temperature T "Temperature.";
    input Integer n "Component.";
    input Real der_T "The derivative of temperature.";
    output Real der_K "The derivative of K-values with respect to temperature.";
  protected 
    constant Pressure p_env = 101325 "Environment pressure";
  algorithm 
    if n == Methanol then
      der_K := -p_env/p_ch3oh(T)^2*dp_ch3oh_dT(T);
    elseif n == Water then
      der_K := p_env/p_h2o(T)^2*dp_h2o_dT(T);
    else
      der_K := 0.0;
    end if;
  end dK_dT;
  
protected 
  function dp_vap_dT 
    input Temperature T;
    input Integer n "Component";
    input Real der_T;
    output Real der_p;
  algorithm 
    if n == 1 then
      der_p := dp_ch3oh_dT(T);
    elseif n == 2 then
      der_p := dp_h2o_dT(T);
    else
      assert(false, "Error: vapour Pressure of component "+String(n)+" requested.");
    end if;
    
  end dp_vap_dT;
  
public 
  function rho 
    input Temperature T;
    input Integer n "Component";
    input Integer p "Phase";
    output Density RHO;
    import Modelica.Constants.R;
  protected 
    constant Pressure p_env=101325.0;
  algorithm 
    if p == LiquidPhase then
      if n == Methanol then
        RHO := 791; // rho_ch3oh(T);
      elseif n == Water then
        RHO := 997; // rho_h2o(T);
      else
        assert(false, "Density in liquid phase of gaseous component "+String(n)+" requested.");
      end if;
    else
      // NOTE Assuming ideal gas.
      RHO := mw(n)*p_env/R/T;
    end if;
    annotation(derivative=drho_dT);
    annotation (Documentation(info="<html>
<p>Returns the density of the given component at the given temperature and in the given phase.
For gases, the ideal gas lawo is assumed; for liquids, a linear interpolation on data is
carried out.</p>
</html>"));
  end rho;
  
  annotation (uses(Modelica(version="2.2.1")), Documentation(info="<html>
<p>This package contains a number of property functions, such as molar enthalpy,
molar heat capacity, vapour pressure, molecular weight and so on.</p>
</html>"));
  
protected 
  function drho_dT "Derivative of the density function" 
    import Modelica.Constants.R;
    
    input Temperature T;
    input Integer n "The component code";
    input Integer p "Phase";
    input Real der_T "The rate of variation of temperature";
    output Real der_rho "The derivative drho/dT";
    
  protected 
    constant Pressure p_env=101325.0;
    constant Temperature dT = 0.01;
  algorithm 
    if p == 2 then
      if n == 1 then
        der_rho := 0; // (rho_ch3oh(T+dT)-rho_ch3oh(T))/dT; // FIXME use spline?
      elseif n == 2 then
        der_rho := 0; // (rho_h2o(T+dT)-rho_h2o(T))/dT;  // FIXME use spline?
      else
        assert(false, "Density in liquid phase of gaseous component "+String(n)+" requested.");
      end if;
    else
      // NOTE Assuming ideal gas.
      der_rho := -mw(n)*p_env/R/T^2*der_T;
    end if;
  end drho_dT;
  
end Thermo;
