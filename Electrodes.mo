package Electrodes "Package containing electrode models" 
  
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
  
  partial model ReactionModelling 
    "Generic modelling of electrochemical reactions" 
    extends Modelica.Electrical.Analog.Interfaces.TwoPin;
    import Thermo.AllSpecies;
    annotation (Icon, Documentation(info="<html>
<p>This model is the abstract representation of the reactions at both membrane
sides. What is always necessary to provide are current, voltage, and species
production or consumption; the latter part is implemented with two species
vectors, one for each side.</p>
</html>"));
    ArealReactionRate[size(AllSpecies,1)] anode_r 
      "Production rates of species on the anode";
    ArealReactionRate[size(AllSpecies,1)] cathode_r 
      "Production rates of species on the cathode";
    
    Modelica.SIunits.Current I = p.i;
    
  end ReactionModelling;
  
  partial model SteadyModelling "A steady-state electrochemical model" 
    extends ReactionModelling;
    
  protected 
    constant Modelica.SIunits.FaradayConstant F = 96485.3415;
    /* This group of vectors represents the coefficients by which
   * current must be multiplied to find the production rate of
   * each species. */
    constant Real[:] cathode_k = {0, 1/2/F, -1/4/F, 0, 0};
    constant Real[:] anode_k = {-1/6/F, -1/6/F, 0, 1/6/F, 0};
  equation 
    anode_r   =   anode_k * I;
    cathode_r = cathode_k * I;
    annotation (Documentation(info="<html>
<p>This is the generic model that couples the reaction rates
directly to current, ignoring any phenomena such as overvoltage
transients. Constant-voltage and Th&eacute;venin models are
inherited from here.</p>
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
    Modelica.Electrical.Analog.Sources.ConstantVoltage V(V=0.5) 
      "Constant voltage" annotation (extent=[-20,-20; 20,20]);
  equation 
    connect(V.n, n) annotation (points=[20,1.22125e-15; 60,1.22125e-15; 60,
          5.55112e-16; 100,5.55112e-16], style(
        color=3, 
        rgbcolor={0,0,255}, 
        fillColor=7, 
        rgbfillColor={255,255,255}, 
        fillPattern=1));
    connect(V.p, p) annotation (points=[-20,1.22125e-15; -60,1.22125e-15; -60,
          5.55112e-16; -100,5.55112e-16], style(
        color=3, 
        rgbcolor={0,0,255}, 
        fillColor=7, 
        rgbfillColor={255,255,255}, 
        fillPattern=1));
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
    Modelica.Electrical.Analog.Sources.ConstantVoltage V0(V=0.7) 
      "Open-circuit voltage" annotation (extent=[-60,-20; -20,20]);
    Modelica.Electrical.Analog.Basic.Resistor R(R=0.005) 
      "Simple linear resistor" annotation (extent=[20,-20; 60,20]);
  equation 
    connect(R.n, n) annotation (points=[60,1.22125e-15; 80,1.22125e-15; 80,
          5.55112e-16; 100,5.55112e-16], style(
        color=3, 
        rgbcolor={0,0,255}, 
        fillColor=7, 
        rgbfillColor={255,255,255}, 
        fillPattern=1));
    connect(R.p, V0.n) annotation (points=[20,1.22125e-15; 0,1.22125e-15; 0,
          1.22125e-15; -20,1.22125e-15], style(
        color=3, 
        rgbcolor={0,0,255}, 
        fillColor=7, 
        rgbfillColor={255,255,255}, 
        fillPattern=1));
    connect(V0.p, p) annotation (points=[-60,1.22125e-15; -80,1.22125e-15; -80,
          5.55112e-16; -100,5.55112e-16], style(
        color=3, 
        rgbcolor={0,0,255}, 
        fillColor=7, 
        rgbfillColor={255,255,255}, 
        fillPattern=1));
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
      annotation (extent=[-70,-10; -50,10]);
    Modelica.Electrical.Analog.Basic.Resistor R(R=-0.005) "Ohmic resistance"
      annotation (extent=[10,-10; 30,10]);
    Modelica.Electrical.Analog.Sources.ConstantVoltage E(V=1.213) 
      "Reversible voltage" annotation (extent=[-30,-10; -10,10]);
    replaceable Electrode cathode "Cathodic model"
      annotation (extent=[50,-10; 70,10]);
  equation 
    anode_r = anode.r;
    cathode_r = cathode.r;
    
    connect(E.p, anode.n) annotation (points=[-30,6.10623e-16; -40,6.10623e-16; 
          -40,6.10623e-16; -50,6.10623e-16], style(
        color=3, 
        rgbcolor={0,0,255}, 
        fillColor=7, 
        rgbfillColor={255,255,255}, 
        fillPattern=1));
    connect(anode.p, p) annotation (points=[-70,6.10623e-16; -86,6.10623e-16; 
          -86,5.55112e-16; -100,5.55112e-16], style(
        color=3, 
        rgbcolor={0,0,255}, 
        fillColor=7, 
        rgbfillColor={255,255,255}, 
        fillPattern=1));
    connect(cathode.n, n) annotation (points=[70,6.10623e-16; 86,6.10623e-16; 
          86,5.55112e-16; 100,5.55112e-16], style(
        color=3, 
        rgbcolor={0,0,255}, 
        fillColor=7, 
        rgbfillColor={255,255,255}, 
        fillPattern=1));
    connect(R.p, E.n) annotation (points=[10,6.10623e-16; 5,6.10623e-16; 5,
          6.10623e-16; 0,6.10623e-16; 0,6.10623e-16; -10,6.10623e-16], style(
        color=3, 
        rgbcolor={0,0,255}, 
        fillColor=7, 
        rgbfillColor={255,255,255}, 
        fillPattern=1));
    connect(R.n, cathode.p) annotation (points=[30,6.10623e-16; 35,6.10623e-16; 
          35,6.10623e-16; 40,6.10623e-16; 40,6.10623e-16; 50,6.10623e-16], 
        style(
        color=3, 
        rgbcolor={0,0,255}, 
        fillColor=7, 
        rgbfillColor={255,255,255}, 
        fillPattern=1));
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
        Line(points=[60,0; 100,0], style(color=3, rgbcolor={0,0,255}))),
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
    Voltage eta(start=0.2) = -v "Overvoltage";
    CurrentDensity ir = 4*F*r1 + F*r2 + F*r3 
      "Proton current density (if eta not frozen)";
    discrete Boolean EtaIsFrozen(start=false);
    
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
      r[Methanol]      = -i/A/6/F;
      r[Water]         = -i/A/6/F;
      r[CarbonDioxide] = i/A/6/F;
    else
      C * der(eta)     = i/A - ir;
      r[Methanol]      = -r1;
      r[Water]         = -r2;
      r[CarbonDioxide] = r3;
    end if;
    r[Oxygen]   = 0;
    r[Nitrogen] = 0;
    
  algorithm 
    when eta <= 0 and not AllowFrozenEta then
      // Condition not allowed by the user, crash.
      assert(false, "==> Overvoltage became negative. Set AllowFrozenEta if you want to continue anyway.");
    elsewhen eta <= 0 then
      // Toggle to short circuit when eta goes below zero.
      EtaIsFrozen := true;
      reinit(eta,0);
    elsewhen EtaIsFrozen and i/A > ir then
      // Toggle back to Tafel when eta starts accumulating again.  
      EtaIsFrozen := false;
    end when;
    
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
    Voltage eta = v "Overvoltage";
    
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
    extends DynamicModelling(redeclare KrewerAnode anode, redeclare 
        KrewerCathode cathode);
    annotation (Diagram);
  end KrewerModel;

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
      connect(anode.n, pulseCurrent.p) annotation (points=[40,-9.4369e-16; 40,
            60; -92,60; -92,-9.4369e-16], style(color=3, rgbcolor={0,0,255}));
    end KrewerAnodeTest;
    
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
        annotation (extent=[50,-20; 70,0]);
      annotation (Diagram);
      parameter Modelica.SIunits.Temperature T = 343 "Cell temperature";
      parameter Modelica.SIunits.MoleFraction xO2 = 0.1 
        "Cathodic oxygen catalyst-layer concentration";
      parameter Flow.MolarFlux Nx = 0.0025 "Crossover methanol flux";
      parameter Modelica.SIunits.Concentration c = 500 
        "Anodic methanol concentration";
      
      KrewerModel km(anode(AllowFrozenEta=true))
        annotation (extent=[-12,-20; 28,20]);
      Modelica.Electrical.Analog.Sources.PulseCurrent pulseCurrent(period=30,
          startTime=10, 
        I=5, 
        offset=10)  annotation (extent=[-88,-30; -28,30]);
    equation 
      km.cathode.N_x = Nx;
      km.cathode.p_O2/101325 = xO2;
      km.cathode.T = T;
      km.anode.T = T;
      km.anode.c = c;
      
      connect(km.n, ground.p) annotation (points=[28,1.22125e-15; 44,
            1.22125e-15; 44,5.55112e-16; 60,5.55112e-16], style(
          color=3, 
          rgbcolor={0,0,255}, 
          fillColor=7, 
          rgbfillColor={255,255,255}, 
          fillPattern=1));
      connect(pulseCurrent.n, km.p) annotation (points=[-28,-9.4369e-16; -20,
            -9.4369e-16; -20,1.22125e-15; -12,1.22125e-15], style(
          color=3, 
          rgbcolor={0,0,255}, 
          fillColor=7, 
          rgbfillColor={255,255,255}, 
          fillPattern=1));
      connect(pulseCurrent.p, km.n) annotation (points=[-88,-9.4369e-16; -88,60; 
            28,60; 28,1.22125e-15], style(
          color=3, 
          rgbcolor={0,0,255}, 
          fillColor=7, 
          rgbfillColor={255,255,255}, 
          fillPattern=1));
    end KrewerModelTest;
  end Test;
end Electrodes;
