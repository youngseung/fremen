package Electrochemistry "Package containing electrochemical models" 
  
type CatalystCoverage = Real(final quantity="Catalyst coverage", final unit="", min=0, max=1) 
    annotation (Documentation(info="<html>
<p>Fraction of catalytic sites occupied by a certain species.</p>
</html>"));
  
type SurfaceConcentration = Real (final quantity="Surface concentration", final unit
        =                                                                            "mol/m2") 
    annotation (Documentation(info="<html>
<p>A unit commonly used for active concentrations of catalysts.</p>
</html>"));
  
type ArealCapacitance = Real (final quantity="Areal capacitance", final unit="F/m2") 
    annotation (Documentation(info="<html>
<p>Unit typically used to indicate the capacitance of charge double layers in electrodes.</p>
</html>"));
  
type ArealReactionRate = Real(final quantity="Areal reaction rate", final unit="mol/(m2.s)") 
    annotation (Documentation(info="<html>
<p>The rate of a reaction on a mole-per-surface-area basis.</p>
</html>"));
  
  model VoltageSource "Interface for voltage sources" 
    extends Modelica.Electrical.Analog.Interfaces.TwoPin;
    
    parameter Modelica.SIunits.Voltage V "Voltage";
    Modelica.SIunits.Current I = -p.i "Current";
    annotation (
      Icon(
        Ellipse(extent=[-50, 50; 50, -50], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=7,
            rgbfillColor={255,255,255})),
        Text(extent=[-100, -120; 100, -80], string="%V V"),
        Text(extent=[-100,80; 100,120],     string="%name"),
        Line(points=[-90,0; 90,0], style(color=0, rgbcolor={0,0,0})),
        Text(
          extent=[-120,50; -20,0],
          style(color=3, rgbcolor={0,0,255}),
          string="+"),
        Text(
          extent=[20,50; 120,0],
          style(color=3, rgbcolor={0,0,255}),
          string="-")),
      Window(
        x=0.31,
        y=0.09,
        width=0.6,
        height=0.6),
      Documentation(info="<html>
<p>A voltage source, reimplemented in order to follow the generator convention,
i.e. current <em>exiting</em> the + pole is positive.
This standard is easier for us as it keeps most values positive.</p>
</html>"),
      Diagram);
  equation 
    v = V;
  end VoltageSource;
  
  partial model ReactionModelling 
    "Generic modelling of electrochemical reactions" 
    extends Modelica.Electrical.Analog.Interfaces.TwoPin;
    import Thermo.AllSpecies;
    annotation (Icon, Documentation(info="<html>
<p>This model is the abstract representation of the reactions at both membrane
sides. What is always necessary to provide are current, voltage, and species
production or consumption; the latter part is implemented with two species
vectors, one for each side.</p>
<p>Note that this unit follows the <em>generator</em> convention, so that 
current <em>exiting</em> the + pole is positive.</p>
</html>"));
    ArealReactionRate[size(AllSpecies,1)] anode_r 
      "Production rates of species on the anode";
    ArealReactionRate[size(AllSpecies,1)] cathode_r 
      "Production rates of species on the cathode";
    
    Modelica.SIunits.Current I = -p.i;
    
  equation 
    p.i + n.i = 0;
    
  end ReactionModelling;
  
  partial model SteadyModelling "A steady-state electrochemical model" 
    extends ReactionModelling;
    
  protected 
    constant Modelica.SIunits.FaradayConstant F = 96485.3415;
    constant Real[:] cathode_k = {0, 1/2, -1/4, 0, 0};
    constant Real[:] anode_k = {-1/6, -1/6, 0, 1/6, 0};
  equation 
    anode_r   =   anode_k * I / F;
    cathode_r = cathode_k * I / F;
    annotation (Documentation(info="<html>
<p>This is the generic model that couples the reaction rates
directly to current, ignoring any phenomena such as overvoltage
transients. Constant-voltage and Th&eacute;venin models are
inherited from here.</p>
<p>In this kind of models, consumption of each species is 
directly related to current.</p>
</html>"));
  end SteadyModelling;
  
  model ConstantVoltage "A constant-voltage circuit" 
    extends SteadyModelling;
    annotation (Diagram, Icon(Ellipse(extent=[-40,40; 40,-40],  style(
            color=3,
            rgbcolor={0,0,255},
            fillColor=7,
            rgbfillColor={255,255,255})), Line(points=[-100,0; 100,0], style(
              color=3, rgbcolor={0,0,255}))),
      Documentation(info="<html>
<p>This class implements the simplest electrochemical modelling, namely a 
constant voltage which can be set as a parameter.</p>
<p>Unless we are specifically interested in dynamic modelling of the voltage,
this model is actually good enough in most cases: voltage in DMFCs does not
change much in the operating range (0.7 to 0.4) compared to the heat-neutral
voltage (&Delta;h<sub>r</sub>/6F = 1.16 V), and a controller unable to cope 
with an uncertainty of that order would anyway be unfit for the real world.</p>
</html>"));
    VoltageSource V(V=0.5) "Constant voltage" 
      annotation (extent=[-30,-30; 30,30]);
  equation 
    connect(V.p, p) annotation (points=[-30,-9.4369e-16; -66,-9.4369e-16; -66,
          5.55112e-16; -100,5.55112e-16], style(color=3, rgbcolor={0,0,255}));
    connect(V.n, n) annotation (points=[30,-9.4369e-16; 66,-9.4369e-16; 66,
          5.55112e-16; 100,5.55112e-16], style(color=3, rgbcolor={0,0,255}));
  end ConstantVoltage;
  
  model Thevenin "A Thevenin-equivalent circuit" 
    extends SteadyModelling;
    annotation (Diagram, Icon(Ellipse(extent=[-60,20; -20,-20], style(
            color=3,
            rgbcolor={0,0,255},
            fillColor=7,
            rgbfillColor={255,255,255})), Line(points=[-100,0; 4,0; 10,20; 20,-20;
              30,20; 40,-20; 50,20; 60,-20; 66,0; 100,0],      style(color=3,
              rgbcolor={0,0,255}))),
      Documentation(info="<html>
<p>This class implements a simple Th&eacute;venin equivalent circuit, allowing
to set its open-circuit voltage and its resistance.</p>
</html>"));
    Modelica.Electrical.Analog.Basic.Resistor R(R=0.005) 
      "Simple linear resistor" annotation (extent=[20,-20; 60,20]);
    VoltageSource V0(V=0.7) "Open-circuit voltage" 
      annotation (extent=[-60,-20; -22,20]);
  equation 
    connect(R.n, n) annotation (points=[60,1.22125e-15; 80,1.22125e-15; 80,
          5.55112e-16; 100,5.55112e-16], style(
        color=3,
        rgbcolor={0,0,255},
        fillColor=7,
        rgbfillColor={255,255,255},
        fillPattern=1));
    connect(V0.n, R.p) annotation (points=[-22,1.22125e-15; -11.5,1.22125e-15;
          -11.5,1.22125e-15; -1,1.22125e-15; -1,1.22125e-15; 20,1.22125e-15],
        style(color=3, rgbcolor={0,0,255}));
    connect(V0.p, p) annotation (points=[-60,1.22125e-15; -80,1.22125e-15; -80,
          5.55112e-16; -100,5.55112e-16], style(color=3, rgbcolor={0,0,255}));
  end Thevenin;
  
  partial model DynamicModelling 
    extends ReactionModelling;
    
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
      annotation (extent=[-56,2; -36,22]);
    Modelica.Electrical.Analog.Basic.Resistor R(R=0.005) "Ohmic resistance" 
      annotation (extent=[-10,2; 10,22]);
    replaceable Electrode cathode "Cathodic model" 
      annotation (extent=[36,2; 56,22]);
    VoltageSource E(V=1.213) "Reversible voltage" 
      annotation (extent=[-10,-30; 10,-10]);
  equation 
    anode_r   = anode.r;
    cathode_r = cathode.r;
    
    connect(E.p, p) annotation (points=[-10,-20; -80,-20; -80,0; -100,0; -100,
          5.55112e-16], style(color=3, rgbcolor={0,0,255}));
    connect(cathode.p, R.n) 
      annotation (points=[36,12; 10,12], style(color=3, rgbcolor={0,0,255}));
    connect(R.p, anode.n) 
      annotation (points=[-10,12; -36,12], style(color=3, rgbcolor={0,0,255}));
    connect(anode.p, n) annotation (points=[-56,12; -80,12; -80,40; 80,40; 80,
          5.55112e-16; 100,5.55112e-16], style(color=3, rgbcolor={0,0,255}));
    connect(cathode.n, E.n) annotation (points=[56,12; 70,12; 70,-20; 10,-20],
        style(color=3, rgbcolor={0,0,255}));
  end DynamicModelling;
  
  partial model Electrode "Generic electrode" 
    extends Modelica.Electrical.Analog.Interfaces.OnePort;
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
    
    import Thermo.AllSpecies;
    import Flow.MolarFlow;
    
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
    import Flow.MolarFlow;
    import Thermo.AllSpecies;
    import Thermo.Methanol;
    import Thermo.Water;
    import Thermo.Oxygen;
    import Thermo.CarbonDioxide;
    import Thermo.Nitrogen;
    
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
  
  model KrewerAnodeSteadyOH 
    "Model of Ulrike's DMFC anode without thetaOH dynamics" 
    extends Electrode;
    
    import Modelica.SIunits.Velocity;
    import Modelica.SIunits.Temperature;
    import Modelica.SIunits.Voltage;
    import Modelica.SIunits.Concentration;
    import Modelica.SIunits.CurrentDensity;
    import Modelica.SIunits.Current;
    import Modelica.Constants.R;
    import Flow.MolarFlow;
    import Thermo.AllSpecies;
    import Thermo.Methanol;
    import Thermo.Water;
    import Thermo.Oxygen;
    import Thermo.CarbonDioxide;
    import Thermo.Nitrogen;
    
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
    Boolean EtaIsFrozen = AllowFrozenEta and eta <= 0 and i/A <= ir 
      "Whether eta has to be kept at a fixed value";
    
  equation 
    r1 = r10 * exp(-betaCO * gCO * (thetaCO - 0.5)) * c * (1 - thetaCO);
    
    r2f = r20f * exp(alpha*F*eta/R/T) * exp(-betaOH*gOH*(thetaOH-0.5)) * (1-thetaOH);
    r2b = r20b * exp(-(1-alpha)*F*eta/R/T) * exp((1-betaOH)*gOH*(thetaOH-0.5)) * thetaOH;
    r2  = r2f - r2b;
    
    r3 = r30 * exp((1 - betaCO) * gCO * (thetaCO - 0.5)) * thetaCO * thetaOH;
    
    r2 = r3;
    
    if EtaIsFrozen and FreezeCoveragesToo then
      der(thetaCO) = 0;
    else
      cPt * der(thetaCO) = r1  - r3;
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
    
    annotation (Documentation(info="<html>
<p>This model is just like <tt>KrewerAnode</tt>, only it assumes a pseudo-steady
state for &Theta;<sub>OH</sub>.</p>
</html>"));
  end KrewerAnodeSteadyOH;
  
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
    import Flow.MolarFlow;
    import Flow.MolarFlux;
    import Thermo.AllSpecies;
    import Thermo.Methanol;
    import Thermo.Water;
    import Thermo.Oxygen;
    import Thermo.CarbonDioxide;
    import Thermo.Nitrogen;
    
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
    extends DynamicModelling(redeclare KrewerAnode anode(AllowFrozenEta=true),
                                                          redeclare 
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
  
  package Test 
    model VoltageSourceTest 
      
      Modelica.Electrical.Analog.Basic.Ground ground 
        annotation (extent=[50,-20; 70,0]);
      annotation (Diagram);
      parameter Modelica.SIunits.Temperature T = 343 "Electrode temperature";
      parameter Modelica.SIunits.Concentration c = 500 "Methanol concentration";
      
      VoltageSource source(V=10) annotation (extent=[-40,-30; 20,30]);
    equation 
      connect(source.n, ground.p) annotation (points=[20,-9.4369e-16; 42,
            -9.4369e-16; 42,5.55112e-16; 60,5.55112e-16], style(color=3,
            rgbcolor={0,0,255}));
    end VoltageSourceTest;
    
    model ConstantVoltageTest 
      
      Modelica.Electrical.Analog.Basic.Ground ground 
        annotation (extent=[50,-40; 70,-20]);
      annotation (Diagram);
      parameter Modelica.SIunits.Temperature T = 343 "Electrode temperature";
      parameter Modelica.SIunits.Concentration c = 500 "Methanol concentration";
      
      ConstantVoltage constantVoltage annotation (extent=[-60,-60; 20,20]);
      Modelica.Electrical.Analog.Sources.PulseCurrent I(
        period=30,
        startTime=10,
        I=5,
        offset=10)  annotation (extent=[-60,0; 20,78]);
    equation 
      connect(ground.p, constantVoltage.n) annotation (points=[60,-20; 40,-20;
            40,-20; 20,-20], style(color=3, rgbcolor={0,0,255}));
      connect(I.n, constantVoltage.n) annotation (points=[20,39; 20,-20], style(
            color=3, rgbcolor={0,0,255}));
      connect(I.p, constantVoltage.p) annotation (points=[-60,39; -60,-20],
          style(color=3, rgbcolor={0,0,255}));
    end ConstantVoltageTest;
    
    model TheveninTest 
      
      Modelica.Electrical.Analog.Basic.Ground ground 
        annotation (extent=[50,-58; 70,-38]);
      annotation (Diagram);
      parameter Modelica.SIunits.Temperature T = 343 "Electrode temperature";
      parameter Modelica.SIunits.Concentration c = 500 "Methanol concentration";
      
      Thevenin thevenin annotation (extent=[-60,-78; 20,2]);
      Modelica.Electrical.Analog.Sources.PulseCurrent I(
        period=30,
        startTime=10,
        I=5,
        offset=10)  annotation (extent=[-60,0; 20,78]);
    equation 
      connect(ground.p, thevenin.n) annotation (points=[60,-38; 20,-38], style(
            color=3, rgbcolor={0,0,255}));
      connect(I.p, thevenin.p) annotation (points=[-60,39; -60,-38; -60,-38],
          style(color=3, rgbcolor={0,0,255}));
      connect(I.n, thevenin.n) annotation (points=[20,39; 20,19.25; 20,19.25;
            20,-0.5; 20,-0.5; 20,-36], style(color=3, rgbcolor={0,0,255}));
    end TheveninTest;
    
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
    
    model KrewerAnodeSteadyOHTest 
      
      KrewerAnodeSteadyOH anode(AllowFrozenEta=true) 
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
    end KrewerAnodeSteadyOHTest;
    
    model KrewerCathodeTest 
      
      KrewerCathode cathode 
                        annotation (extent=[-20,-30; 40,30]);
      Modelica.Electrical.Analog.Basic.Ground ground 
        annotation (extent=[50,-20; 70,0]);
      annotation (Diagram);
      parameter Modelica.SIunits.Temperature T = 343 "Electrode temperature";
      parameter Modelica.SIunits.MoleFraction xO2 = 0.1 
        "Oxygen catalyst-layer concentration";
      parameter Flow.MolarFlux Nx = 0.0025 "Crossover methanol flux";
      
      Modelica.Electrical.Analog.Sources.PulseCurrent pulseCurrent(period=30,
          startTime=10) 
                    annotation (extent=[-92,-30; -32,30]);
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
      
      Modelica.Electrical.Analog.Basic.Ground ground 
        annotation (extent=[30,-40; 50,-20]);
      annotation (Diagram);
      parameter Modelica.SIunits.Temperature T = 343 "Cell temperature";
      parameter Modelica.SIunits.MoleFraction xO2 = 0.1 
        "Cathodic oxygen catalyst-layer concentration";
      parameter Flow.MolarFlux Nx = 0.0025 "Crossover methanol flux";
      parameter Modelica.SIunits.Concentration c = 500 
        "Anodic methanol concentration";
      
      KrewerModel km(anode(AllowFrozenEta=true)) 
        annotation (extent=[-40,-50; 20,10]);
      Modelica.Electrical.Analog.Sources.PulseCurrent pulseCurrent(period=30,
          startTime=10,
        I=5,
        offset=10)  annotation (extent=[-40,0; 20,60]);
    equation 
      km.cathode.N_x = Nx;
      km.cathode.p_O2/101325 = xO2;
      km.anode.T = T;
      km.anode.c = c;
      
      connect(km.n, ground.p) annotation (points=[20,-20; 20,-16; 24,-20; 40,
            -20],                                         style(
          color=3,
          rgbcolor={0,0,255},
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
      connect(pulseCurrent.p, km.p) annotation (points=[-40,30; -40,-20], style(
            color=3, rgbcolor={0,0,255}));
      connect(pulseCurrent.n, km.n) annotation (points=[20,30; 20,17.5; 20,17.5; 
            20,5; 20,-20; 20,-20], style(color=3, rgbcolor={0,0,255}));
    end KrewerModelTest;
    annotation (Documentation(info="<html>
<p>Test cases for the electrochemical models.</p>
</html>"));
  end Test;
end Electrochemistry;
