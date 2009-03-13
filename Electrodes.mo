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
  
type ArealReactionRate = Real(final quantity="Areal reaction rate", final unit="mol/m2s") 
    annotation (Documentation(info="<html>
<p>The rate of a reaction on a mole-per-surface-area basis.</p>
</html>"));
  
  model Thevenin "A Thevenin-equivalent circuit" 
    extends Modelica.Electrical.Analog.Interfaces.OnePort;
    
    parameter Modelica.SIunits.Voltage V0 = 0.7 "The open-circuit voltage";
    parameter Modelica.SIunits.Resistance R = 0.1 "The equivalent resistance";
    annotation (Diagram, Icon(Ellipse(extent=[-60,20; -20,-20], style(
            color=3,
            rgbcolor={0,0,255},
            fillColor=7,
            rgbfillColor={255,255,255})), Line(points=[-100,0; 14,0; 20,20; 30,
              -20; 40,20; 50,-20; 60,20; 70,-20; 76,0; 100,0], style(color=3,
              rgbcolor={0,0,255}))), 
      Documentation(info="<html>
<p>This class implements a simple Th&eacute;venin equivalent circuit, allowing
to set its open-circuit voltage and its resistance.</p>
</html>"));
  equation 
    v = V0-R*i;
  end Thevenin;
  
  partial model Electrode "Generic electrode" 
    extends Modelica.Electrical.Analog.Interfaces.OnePort;
    
    parameter Modelica.SIunits.Area A = 26E-4 "Active electrode area";
    annotation (Icon(
        Rectangle(extent=[-60,20; 60,-20], style(
            color=3,
            rgbcolor={0,0,255},
            fillColor=7,
            rgbfillColor={255,255,255})),
        Line(points=[-100,0; -60,0], style(color=3, rgbcolor={0,0,255})),
        Line(points=[60,0; 100,0], style(color=3, rgbcolor={0,0,255}))));
  protected 
    constant Modelica.SIunits.FaradayConstant F = 96485.3415;
  end Electrode;
  
  model ConstantElectrode "Electrode with constant overvoltage" 
    extends Electrode;
    
    Modelica.Electrical.Analog.Sources.ConstantVoltage eta(V=0) 
      annotation (extent=[-40,-40; 40,40]);
  equation 
    connect(eta.n, n) annotation (points=[40,2.44249e-15; 70,2.44249e-15; 70,
          5.55112e-16; 100,5.55112e-16], style(color=3, rgbcolor={0,0,255}));
    connect(eta.p, p) annotation (points=[-40,2.44249e-15; -70,2.44249e-15; -70,
          5.55112e-16; -100,5.55112e-16], style(color=3, rgbcolor={0,0,255}));
  end ConstantElectrode;
  
  model KrewerAnode "Full model of Ulrike's DMFC anode" 
    extends Electrode;
    
    import Modelica.SIunits.Velocity;
    import Modelica.SIunits.Temperature;
    import Modelica.SIunits.Voltage;
    import Modelica.SIunits.Concentration;
    import Modelica.SIunits.CurrentDensity;
    import Modelica.SIunits.Current;
    import Modelica.Constants.R;
    import Flow.MolarFlowRate;
    import Thermo.AllSpecies;
    import Thermo.Methanol;
    import Thermo.Water;
    import Thermo.Oxygen;
    import Thermo.CarbonDioxide;
    import Thermo.Nitrogen;
    
    parameter Boolean AllowFrozenEta = false "Allows to handle open circuit";
    parameter Boolean FreezeCoveragesToo = true 
      "Block coverages when the overvoltage is";
    
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
    
    MolarFlowRate[size(AllSpecies,1)] r "Rate of production of species";
    
    CatalystCoverage thetaCO(start=0.9) "Catalyst coverage of CO on Pt";
    CatalystCoverage thetaOH(start=0.01) "Catalyst coverage of OH on Ru";
    
    Concentration c "Catalyst-layer methanol concentration";
    Temperature T "Electrode temperature";
    Voltage eta(start=0.2) = v "Anodic overvoltage";
    CurrentDensity ir = 4*F*r1 + F*r2 + F*r3 
      "Proton current density (if eta not frozen)";
    Current Ir = ir*A "Reaction current";
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
 
<p>{1}: Krewer, Ulrike, Yoon, Hae-Kwon, and Kim, Hee-Tak: Basic model for membrane electrode assembly 
design for direct methanol fuel cells, Journal of Power Sources, 760-772, 2008.</p>
</html>"));
  end KrewerAnode;
  annotation (uses(Modelica(version="2.2.1")));
  
  package Test 
    model KrewerAnodeTest 
      
      KrewerAnode anode annotation (extent=[-20,-30; 40,30]);
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
    
  end Test;
end Electrodes;
