import " Units.mo";


package Electrochemistry "Package containing electrochemical models" 
  
  partial model DynamicModelling 
    
    import Modelica.Electrical.Analog.Basic.Resistor;
    import Modelica.Electrical.Analog.Sources.ConstantVoltage;
    import Modelica.Electrical.Analog.Interfaces.PositivePin;
    import Modelica.Electrical.Analog.Interfaces.NegativePin;
    import Modelica.SIunits.Current;
    import Modelica.SIunits.Voltage;
    
    annotation (Diagram, Icon(
        Rectangle(extent=[-60,20; 60,-20], style(
            color=3,
            rgbcolor={0,0,255},
            fillColor=7,
            rgbfillColor={255,255,255})),
        Line(points=[-100,0; -60,0], style(
            color=3,
            rgbcolor={0,0,255},
            fillColor=7,
            rgbfillColor={255,255,255},
            fillPattern=1)),
        Line(points=[60,0; 100,0], style(
            color=3,
            rgbcolor={0,0,255},
            fillColor=7,
            rgbfillColor={255,255,255},
            fillPattern=1))));
    replaceable Electrode anode "Anodic model" 
      annotation (extent=[-70,-10; -50,10]);
    Resistor R(R=0.005) "Ohmic resistance" 
      annotation (extent=[-30,-10; -10,10]);
    replaceable Electrode cathode "Cathodic model" 
      annotation (extent=[50,-10; 70,10]);
    ConstantVoltage E(V=1.213) "Reversible voltage" 
      annotation (extent=[10,10; 30,30]);
    PositivePin p "Positive pin" annotation (extent=[-110,-10; -90,10]);
    NegativePin n "Negative pin" annotation (extent=[90,-10; 110,10]);
    
    Current I = -p.i;
    Voltage V = p.v - n.v;
    
  equation 
    connect(R.p, anode.n) 
      annotation (points=[-30,6.10623e-16; -42,-3.36456e-22; -42,6.10623e-16; 
          -50,6.10623e-16],                style(color=3, rgbcolor={0,0,255}));
    connect(p, cathode.n) annotation (points=[-100,5.55112e-16; -96,5.55112e-16; 
          -96,0; -90,0; -90,-20; 80,-20; 80,6.10623e-16; 70,6.10623e-16],
                                                                      style(
          color=3, rgbcolor={0,0,255}));
    connect(cathode.p, R.n) annotation (points=[50,6.10623e-16; 35,6.10623e-16; 
          35,6.10623e-16; 20,6.10623e-16; 20,6.10623e-16; -10,6.10623e-16],
        style(color=3, rgbcolor={0,0,255}));
    connect(anode.p, E.p) annotation (points=[-70,6.10623e-16; -80,6.10623e-16; 
          -80,20; 10,20], style(color=3, rgbcolor={0,0,255}));
    connect(E.n, n) annotation (points=[30,20; 90,20; 90,5.55112e-16; 100,
          5.55112e-16], style(color=3, rgbcolor={0,0,255}));
  end DynamicModelling;
  
  partial model Electrode "Generic electrode" 
    extends Modelica.Electrical.Analog.Interfaces.OnePort;
    
    import Thermo.AllSpecies;
    import Units.MolarFlow;
    
    annotation (Icon(
        Rectangle(extent=[-60,20; 60,-20], style(
            color=3,
            rgbcolor={0,0,255},
            fillColor=7,
            rgbfillColor={255,255,255})),
        Line(points=[-100,0; -60,0], style(color=3, rgbcolor={0,0,255})),
        Line(points=[60,0; 100,0], style(color=3, rgbcolor={0,0,255})),
        Text(extent=[-100,60; 100,100],     string="%name")),
        Documentation(info="<html>
<p>The generic interface of an electrode (be it anode or cathode).
</html>"));
    
    parameter Modelica.SIunits.Area A = 26E-4 "Active electrode area";
    
    MolarFlow[size(AllSpecies,1)] r 
      "Production of each species due to reaction";
  protected 
    constant Modelica.SIunits.FaradayConstant F = 96485.3415 
      "Faraday's constant";
  end Electrode;
  
  model KrewerAnode "Full model of Ulrike's DMFC anode" 
    extends Electrode;
    
    import Modelica.SIunits.Velocity;
    import Modelica.SIunits.Temperature;
    import Modelica.SIunits.Voltage;
    import Modelica.SIunits.Concentration;
    import Modelica.SIunits.CurrentDensity;
    import Modelica.SIunits.Current;
    import Modelica.Constants.R;
    import Thermo.AllSpecies;
    import Thermo.Methanol;
    import Thermo.Water;
    import Thermo.Oxygen;
    import Thermo.CarbonDioxide;
    import Thermo.Nitrogen;
    import Units.ArealCapacitance;
    import Units.SurfaceConcentration;
    import Units.ArealReactionRate;
    import Units.CatalystCoverage;
    
    parameter Boolean AllowFrozenEta = false 
      "Allows to handle open-circuit and low-current conditions";
    parameter Boolean FreezeCoveragesToo = true 
      "Block coverages too when the overvoltage is blocked";
    
    parameter ArealCapacitance C = 3348 "The electrode's areal capacitance";
    parameter Real alpha = 0.5 
      "Charge transfer coefficient for water adsorption";
    parameter Real betaCO = 0.5 
      "Symmetry parameter for Frumkin-Temkin adsorption on Pt";
    parameter Real betaOH = 0.5 
      "Symmetry parameter for Frumkin-Temkin adsorption on Ru";
    parameter Real gCO = 11 
      "Inhomogeneity factor for Frumkin-Temkin adsorption on Pt";
    parameter Real gOH = 0.43 
      "Interaction factor for Frumkin-Temkin adsorption on Ru";
    parameter SurfaceConcentration cPt = 0.117 
      "Active surface concentration of Pt";
    parameter SurfaceConcentration cRu = 0.165 
      "Active surface concentration of Ru";
    
    parameter Velocity r10 = 1.6E-4 
      "Reaction rate constant for methanol adsorption on Pt";
    parameter ArealReactionRate r20f = 7.2E-4 
      "Reaction rate constant for water adsorption on Ru";
    parameter ArealReactionRate r20b = 9.91E4 
      "Reaction rate constant for water desorption from Ru";
    parameter ArealReactionRate r30 = 0.19 
      "Reaction rate constant for carbon-dioxide production";
    
    ArealReactionRate r1 "Rate of methanol adsorption on Pt (4 protons)";
    ArealReactionRate r2f "Rate of water adsorption on Ru";
    ArealReactionRate r2b "Rate of water desorption from Ru";
    ArealReactionRate r2(stateSelect=StateSelect.always) 
      "Net rate of water adsorption on Ru (1 proton)";
    ArealReactionRate r3 "Rate of carbon-dioxide production (1 proton)";
    
    CatalystCoverage thetaCO(start=0.9) "Catalyst coverage of CO on Pt";
    CatalystCoverage thetaOH(start=0.01) "Catalyst coverage of OH on Ru";
    
    Concentration c "Catalyst-layer methanol concentration";
    Temperature T "Electrode temperature";
    Voltage eta(start=0.2) = v "Overvoltage";
    CurrentDensity ir = 4*F*r1 + F*r2 + F*r3 
      "Proton current density (if eta not frozen)";
    Boolean EtaIsFrozen(start=false) = AllowFrozenEta and eta <= 0 and i/A <= ir 
      "Whether eta has to be kept at a fixed value";
    
  equation 
    r1 = r10 * exp(-betaCO * gCO * (thetaCO - 0.5)) * c * (1 - thetaCO);
    
    r2f = r20f * exp(alpha*F*eta/R/T) * exp(-betaOH*gOH*(thetaOH-0.5)) * (1-thetaOH);
    r2b = r20b * exp(-(1-alpha)*F*eta/R/T) * exp((1-betaOH)*gOH*(thetaOH-0.5)) * thetaOH;
    r2  = r2f - r2b;
    
    r3 = r30 * exp((1 - betaCO) * gCO * (thetaCO - 0.5)) * thetaCO * thetaOH;
    
    if EtaIsFrozen and FreezeCoveragesToo then
      der(thetaCO) = 0;
      der(thetaOH) = 0;
    else
      cPt * der(thetaCO) = r1  - r3;
      cRu * der(thetaOH) = r2  - r3;
    end if;
    
    if EtaIsFrozen then
      der(eta)         = 0;
      r[Methanol]      = -i/6/F;
      r[Water]         = -i/6/F;
      r[CarbonDioxide] = i/6/F;
    else
      C * der(eta)     = i/A - ir;
      r[Methanol]      = -r1*A;
      r[Water]         = -r2*A;
      r[CarbonDioxide] = r3*A;
    end if;
    r[Oxygen]   = 0;
    r[Nitrogen] = 0;
    
    assert( AllowFrozenEta or eta >= 0, "==> Overvoltage became negative. Set AllowFrozenEta if you want to continue anyway.");
    
  initial equation 
    der(thetaOH) = 0;
    der(thetaCO) = 0;
    annotation (Documentation(info="<html>
<p>This class implements a DMFC anode as described by Krewer et al. [1]. Equations and parameters
are taken from there.</p>
 
<p>Krewer et al.'s model cannot work at open circuit, because reactions 1 (CO adsorption) and 3 
(CO<sub>2</sub> desorption) cannot reach a steady state. It was necessary to use hybrid modelling:
when overvoltage &eta; falls below zero, the model switches to a short circuit, with &eta;=V=0; when
then the conditions are present for an increase in &eta;, the model switches back to its full model.</p>
 
<p>Note on initialisation: it is not possible to <em>reliably</em> initialise &eta; to its steady state
at simulation start, so it is set to some &eta;<sub>0</sub>; the reason it is not possible is that, for a 
small yet not easily calculated range of current just above 0, &eta; cannot converge to a steady state.</p>
 
<p>Note on state selection: it is critical to set r<sub>2</sub> as a state, otherwise Dymola will choose 
&Theta;<sub>CO</sub>, &Theta;<sub>OH</sub> and &eta;; this will result in an imprecise value for 
r<sub>2</sub>, since it will be the difference between two large numbers, and numerical noise will be
amplified.</p>
 
<h3>References</h3>
<p>Krewer, Ulrike, Yoon, Hae-Kwon, and Kim, Hee-Tak: Basic model for membrane electrode assembly 
design for direct methanol fuel cells, Journal of Power Sources, 760-772, 2008.</p>
</html>"));
  end KrewerAnode;
  annotation (uses(Modelica(version="2.2.1")));
  
  model KrewerCathode "Ulrike's model of a DMFC cathode" 
    extends Electrode;
    annotation (Documentation(info="<html>
<p>This class implements a DMFC cathode as described by Krewer et al. with some modifications.</p>
 
<p>A term accounting for oxygen partial pressure was added; Krewer et al. used air on the cathode, 
so we use the dry-air value for p<sub>O<sub>2</sub></sub> as a reference. Also, the sign of &eta;
has been inverted so that it is </p>
 
<p>This model could misbehave if there is no reaction (neither main nor cross-over), inducing
a positive potential. However, this is unlikely to happen as even a small cross-over will put
the system back into place. If necessary, an inverse-reaction term could be inserted into
r<sub>c</sub>.<p>
 
<h3>References</h3> 
<p>Krewer, Ulrike, Yoon, Hae-Kwon, and Kim, Hee-Tak: Basic model for membrane electrode assembly 
design for direct methanol fuel cells, Journal of Power Sources, 760-772, 2008.</p>
</html>"));
    import Modelica.SIunits.Velocity;
    import Modelica.SIunits.Temperature;
    import Modelica.SIunits.Voltage;
    import Modelica.SIunits.Concentration;
    import Modelica.SIunits.CurrentDensity;
    import Modelica.SIunits.Current;
    import Modelica.SIunits.MoleFraction;
    import Modelica.Constants.R;
    import Thermo.AllSpecies;
    import Thermo.Methanol;
    import Thermo.Water;
    import Thermo.Oxygen;
    import Thermo.CarbonDioxide;
    import Thermo.Nitrogen;
    import Units.ArealCapacitance;
    import Units.ArealReactionRate;
    import Units.MolarFlux;
    
    parameter ArealCapacitance C = 907 "The electrode's areal capacitance";
    parameter Real alpha = 0.18 "Charge transfer coefficient";
    parameter ArealReactionRate r0 = 4E-6 "Reaction-rate constant";
    
    ArealReactionRate rc "Rate of reaction in =>oxidation<= direction";
    MoleFraction p_O2 "Catalyst-layer oxygen molar fraction";
    MolarFlux N_x "The methanol cross-over flux";
    Temperature T "Electrode temperature";
    Voltage eta = -v "Overvoltage";
    
  protected 
    constant MoleFraction p_O20 = 0.21*101325 
      "Oxygen molar fraction in reference conditions";
    constant Real[:] cathode_k = {0, 3, -3/2, 0, 0};
    
  equation 
    rc = -r0 * (p_O2/p_O20)^1.5 * exp(-(1-alpha)*F*eta/R/T);
    
    C * der(eta) = - i/A - 6*F*rc - 6*F*N_x;
    
    r = cathode_k * (-rc) * A;
    
  end KrewerCathode;
  
  model KrewerModel "Ulrike's electrochemical model" 
    extends DynamicModelling(redeclare KrewerAnode anode(AllowFrozenEta=true),redeclare 
        KrewerCathode cathode);
    annotation (Diagram, Documentation(info="<html>
<p>This class implements the model by Krewer et al., with a voltage generator, a resistance,
and the two electrodes. It also sets the anode's temperature equal to the cathode's.</p>
 
<h3>References</h3> 
<p>Krewer, Ulrike, Yoon, Hae-Kwon, and Kim, Hee-Tak: Basic model for membrane electrode assembly 
design for direct methanol fuel cells, Journal of Power Sources, 760-772, 2008.</p>
</html>"));
  equation 
    anode.T = cathode.T;
  end KrewerModel;
  
  model Diffusion "A one-dimensional model for multicomponent diffusion" 
    import Modelica.SIunits.Density;
    import Modelica.SIunits.MolarMass;
    import Modelica.SIunits.SurfaceTension;
    import Modelica.SIunits.Angle;
    import Modelica.SIunits.Temperature;
    import Modelica.SIunits.PartialPressure;
    import Modelica.SIunits.MoleFraction;
    import Modelica.SIunits.Velocity;
    import Modelica.SIunits.Pressure;
    import Modelica.SIunits.Length;
    import Modelica.SIunits.CurrentDensity;
    import Modelica.Constants.R;
    
    import Thermo.pDoverT175;
    import Thermo.AllSpecies;
    
    import Units.DynamicViscosity;
    import Units.MolarFlux;
    import Units.Permeability;
    import Units.Porosity;
    import Units.Saturation;
    import Units.CondensationCoefficient;
    import Units.ReactionRate;
    
    parameter Integer n(min=3) = 10 "Discretisation steps";
    
    parameter Length t = 0.5E-3 "Layer thickness";
    parameter Porosity eps = 0.5 "Layer porosity";
    parameter Permeability K = 2.55E-13 "Absolute permeability";
    parameter DynamicViscosity mu = 405E-6 "Liquid viscosity";
    parameter Saturation s_im= 0.1 "Immobile saturation";
    parameter SurfaceTension sigma = 64.4E-3 
      "Surface tension between water and air";
    parameter Angle theta = 1.0472 
      "Contact angle between water, air and backing";
    parameter CondensationCoefficient gamma = 900 
      "Condensation coefficient of water";
    
    parameter MolarFlux N_H "Proton flow"; // TODO find values
    parameter MolarFlux N_x "Crossover flow";
    
    parameter MoleFraction[size(AllSpecies,1)] y_bulk = {0, 0.2, 0.1, 0.05, 0.65} 
      "Bulk concentrations";
    MolarFlux[size(AllSpecies,1)] N_0 "Molar fluxes at catalyst layer";
    
    constant MolarMass M = Thermo.mw(Thermo.Water) "Water molar mass";
    constant Pressure p = 101325 "Atmospheric pressure";
    
    Temperature T "Temperature";
    Density rho = Thermo.rho(T, Thermo.Water, Thermo.LiquidPhase) 
      "Liquid density";
    PartialPressure p_vap = Thermo.p_vap(T, Thermo.Water) 
      "Liquid partial pressure";
    
    Velocity[n] liquidVelocity "Average liquid velocity";
    Saturation[n] s "Liquid water saturation (volumetric)";
    Saturation[n] S "Reduced liquid saturation";
    Permeability[n] K_rel "Relative permeability of liquid water";
    
    MoleFraction[n,size(AllSpecies,1)] y "Gas molar fractions";
    MolarFlux[n+1,size(AllSpecies,1)] N "Gas molar flux";
    ReactionRate[n] e "Water evaporation rate";
    
  protected 
    parameter Real epsFactor = eps*((eps-0.11)/(1-0.11))^0.785 
      "Effective diffusivity-porosity factor";
    parameter Length deltaX = t/n "Thickness of each sub-layer";
    
  equation 
    // Boundary condition for flow at catalyst layer
    N_0 = N_H*{0, 0.5, -0.25, 0, 0} + N_x*{0, 2, -1.5, 1, 0};
    
    // Continuity equation FIXME only water has R
    (p/R/T) * y[1,:] = (N_0-N[1,:])/deltaX;// + R[1];
    for i in 2:n loop
      (p/R/T) * der(y[i,:]) = (N[i-1,:]-N[i,:])/deltaX;// + R[i];
    end for;
    
    // Stefan-Maxwell diffusion
    for i in 1:n-1 loop
      for j in 1:size(AllSpecies,1) loop
        (p/R/T) * (y[i+1,j]-y[i,j]) / deltaX =
          sum( ((y[i+1,j]+y[i,j])/2*N[i,k] - (y[i+1,k]+y[i,k])/2*N[i,j]) /
               (pDoverT175(j,k) / p * T^1.75 * epsFactor * (1-(s[i+1]+s[i])/2)^2) 
               for k in 1:size(AllSpecies,1));
      end for;
    end for;
    // Stefan-Maxwell diffusion at last discretisation step (interface to bulk)
    for j in 1:size(AllSpecies,1) loop
      (p/R/T) * (y_bulk[j]-y[n,j]) / deltaX =
        sum( ((y_bulk[j]+y[n,j])/2*N[n,k] - (y_bulk[k]+y[n,k])/2*N[n,j]) /
             (pDoverT175(j,k) / p * T^1.75 * epsFactor * (1-s[n])^2) 
             for k in 1:size(AllSpecies,1));
    end for;
    
    // TODO all equations for the liquid phase
    
  end Diffusion;

  package Test 
    
    model KrewerAnodeTest 
      
      KrewerAnode anode(AllowFrozenEta=true) 
                        annotation (extent=[-20,-30; 40,30]);
      Modelica.Electrical.Analog.Basic.Ground ground 
        annotation (extent=[50,-20; 70,0]);
      annotation (Diagram);
      parameter Modelica.SIunits.Temperature T = 343 "Electrode temperature";
      parameter Modelica.SIunits.Concentration c = 500 "Methanol concentration";
      
      Modelica.Electrical.Analog.Sources.PulseCurrent pulseCurrent(period=30,
          startTime=10) 
                    annotation (extent=[-92,-30; -32,30]);
    equation 
      anode.c = c;
      anode.T = T;
      
      connect(anode.n, ground.p) annotation (points=[40,-9.4369e-16; 31,
            -9.4369e-16; 31,5.55112e-16; 60,5.55112e-16],    style(color=3,
            rgbcolor={0,0,255}));
      connect(pulseCurrent.n, anode.p) annotation (points=[-32,-9.4369e-16; -29,
            -9.4369e-16; -29,-9.4369e-16; -26,-9.4369e-16; -26,-9.4369e-16; -20,
            -9.4369e-16], style(color=3, rgbcolor={0,0,255}));
      connect(pulseCurrent.p, anode.n) annotation (points=[-92,-9.4369e-16; -92,
            60; 40,60; 40,-9.4369e-16], style(color=3, rgbcolor={0,0,255}));
    end KrewerAnodeTest;
    
    model KrewerCathodeTest 
      
      import Modelica.SIunits.Temperature;
      import Modelica.SIunits.MoleFraction;
      import Units.MolarFlux;
      KrewerCathode cathode 
                        annotation (extent=[-20,-30; 40,30]);
      Modelica.Electrical.Analog.Basic.Ground ground 
        annotation (extent=[50,-20; 70,0]);
      annotation (Diagram);
      Modelica.Electrical.Analog.Sources.PulseCurrent pulseCurrent(period=30,
          startTime=10) 
                    annotation (extent=[-92,-30; -32,30]);
      parameter Temperature T = 343 "Electrode temperature";
      parameter MoleFraction xO2 = 0.1 "Oxygen catalyst-layer concentration";
      parameter MolarFlux Nx = 0.0025 "Crossover methanol flux";
      
    equation 
      cathode.N_x = Nx;
      cathode.p_O2/101325 = xO2;
      cathode.T = T;
      
      connect(cathode.n, ground.p) 
                                 annotation (points=[40,-9.4369e-16; 31,
            -9.4369e-16; 31,5.55112e-16; 60,5.55112e-16],    style(color=3,
            rgbcolor={0,0,255}));
      connect(pulseCurrent.n, cathode.p) 
                                       annotation (points=[-32,-9.4369e-16; -26,
            4.29069e-22; -26,-9.4369e-16; -20,-9.4369e-16],
                          style(color=3, rgbcolor={0,0,255}));
      connect(cathode.n, pulseCurrent.p) 
                                       annotation (points=[40,-9.4369e-16; 40,
            60; -92,60; -92,-9.4369e-16], style(color=3, rgbcolor={0,0,255}));
    end KrewerCathodeTest;
    
    model KrewerModelTest 
      
      import Modelica.SIunits.Temperature;
      import Modelica.SIunits.MoleFraction;
      import Modelica.SIunits.Concentration;
      import Units.MolarFlux;
      Modelica.Electrical.Analog.Basic.Ground ground 
        annotation (extent=[30,-40; 50,-20]);
      annotation (Diagram);
      KrewerModel km(anode(AllowFrozenEta=true)) 
        annotation (extent=[-40,-50; 20,10]);
      Modelica.Electrical.Analog.Sources.PulseCurrent pulseCurrent(period=30,
          startTime=10,
        I=5,
        offset=10)  annotation (extent=[-40,0; 20,60]);
      parameter Temperature T = 343 "Cell temperature";
      parameter MoleFraction xO2 = 0.1 
        "Cathodic oxygen catalyst-layer concentration";
      parameter MolarFlux Nx = 0.0025 "Crossover methanol flux";
      parameter Concentration c = 500 "Anodic methanol concentration";
      
    equation 
      km.cathode.N_x = Nx;
      km.cathode.p_O2/101325 = xO2;
      km.anode.T = T;
      km.anode.c = c;
      
      connect(km.n, ground.p)       annotation (points=[20,-20; 20,17.5; 20,
            17.5; 20,5; 20,-20; 40,-20],
                                   style(color=3, rgbcolor={0,0,255}));
      connect(pulseCurrent.p, km.p) annotation (points=[-40,30; -40,-20],
                                   style(color=3, rgbcolor={0,0,255}));
      connect(pulseCurrent.n, km.n) annotation (points=[20,30; 20,17.5; 20,17.5;
            20,5; 20,-20; 20,-20], style(color=3, rgbcolor={0,0,255}));
    end KrewerModelTest;
    annotation (Documentation(info="<html>
<p>Test cases for the electrochemical models.</p>
</html>"));
  end Test;
end Electrochemistry;
