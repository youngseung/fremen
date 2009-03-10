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
  
type ArealCapacitance = Real (final quantity="Areal capacitance", final unit="F.m2") 
    annotation (Documentation(info="<html>
<p>Unit typically used to indicate the capacitance of charge double layers in electrodes.</p>
</html>"));
  
type ArealReactionRate = Real(final quantity="Areal reaction rate", final unit="mol/m2s") 
    annotation (Documentation(info="<html>
<p>The rate of a reaction on a mole-per-surface-area basis.</p>
</html>"));
  
  
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
    import Modelica.Constants.R;
    
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
    ArealReactionRate r2 "Net rate of water adsorption on Ru (1 proton)";
    ArealReactionRate r3 "Rate of carbon-dioxide production (1 proton)";
    
    CatalystCoverage thetaCO(start=0.5) "Catalyst coverage of CO on Pt";
    CatalystCoverage thetaOH(start=0.05) "Catalyst coverage of OH on Ru";
    
    Concentration c "Catalyst-layer methanol concentration";
    Temperature T "Electrode temperature";
    Voltage eta(start=0.2) = v "Anodic overvoltage";
    CurrentDensity ir = 4*F*r1 + F*r2 + F*r3 "Proton current density";
    
  equation 
    r1 = r10 * exp(-betaCO * gCO * (thetaCO - 0.5)) * c * (1 - thetaCO);
    
    r2f = r20f * exp(alpha*F*eta/R/T) * exp(-betaOH*gOH*(thetaOH-0.5)) * (1-thetaOH);
    r2b = r20b * exp(-(1-alpha)*F*eta/R/T) * exp((1-betaOH)*gOH*(thetaOH-0.5)) * thetaOH;
    r2  = r2f - r2b;
    
    r3 = r30 * exp((1 - betaCO) * gCO * (thetaCO - 0.5)) * thetaCO * thetaOH;
    
    cPt * der(thetaCO) = r1 - r3;
    cRu * der(thetaOH) = r2 - r3;
    
    C * der(eta) = i/A - ir;
    
    /* TODO decouple continuous model and OnePort::v/i
   * TODO if current is lower than minimum threshold (i_0)
   * TODO switch model to short circuit, else connect
   * eta to v and ir*A to i. */
    
  initial equation 
    der(thetaCO) = 0;
    der(thetaOH) = 0;
    der(eta) = 0;
    
  end KrewerAnode;
  annotation (uses(Modelica(version="2.2.1")));
  model KrewerAnodeSteady "Full model of Ulrike's DMFC anode" 
    extends Electrode;
    
    import Modelica.SIunits.Velocity;
    import Modelica.SIunits.Temperature;
    import Modelica.SIunits.Voltage;
    import Modelica.SIunits.Concentration;
    import Modelica.SIunits.CurrentDensity;
    import Modelica.Constants.R;
    
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
    ArealReactionRate r2 "Net rate of water adsorption on Ru (1 proton)";
    ArealReactionRate r3 "Rate of carbon-dioxide production (1 proton)";
    
    CatalystCoverage thetaCO "Catalyst coverage of CO on Pt";
    CatalystCoverage thetaOH "Catalyst coverage of OH on Ru";
    
    Concentration c "Catalyst-layer methanol concentration";
    Temperature T "Electrode temperature";
    Voltage eta = v "Anodic overvoltage";
    
  equation 
    r1 = r10 * exp(-betaCO * gCO * (thetaCO - 0.5)) * c * (1 - thetaCO);
    
    r2f = r20f * exp(alpha*F*eta/R/T) * exp(-betaOH*gOH*(thetaOH-0.5)) * (1-thetaOH);
    r2b = r20b * exp(-(1-alpha)*F*eta/R/T) * exp((1-betaOH)*gOH*(thetaOH-0.5)) * thetaOH;
    r2  = r2f - r2b;
    
    r3 = r30 * exp((1 - betaCO) * gCO * (thetaCO - 0.5)) * thetaCO * thetaOH;
    
    r1 = r3;
    r2 = r3;
    
    i = 4*F*r1 + F*r2 + F*r3;
    
  end KrewerAnodeSteady;
  
  package Test 
    model KrewerAnodeTest 
      
      KrewerAnode anode annotation (extent=[-20,-30; 40,30]);
      Modelica.Electrical.Analog.Basic.Ground ground 
        annotation (extent=[50,-20; 70,0]);
      annotation (Diagram);
      parameter Modelica.SIunits.Temperature T = 343 "Electrode temperature";
      parameter Modelica.SIunits.Concentration c = 500 "Methanol concentration";
      
      Modelica.Electrical.Analog.Sources.PulseCurrent pulseCurrent(offset=1, 
          period=5) annotation (extent=[-92,-30; -32,30]);
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
    
    model KrewerAnodeSteadyTest 
      
      KrewerAnodeSteady anode 
                        annotation (extent=[-60,-30; 0,30]);
      Modelica.Electrical.Analog.Basic.Ground ground 
        annotation (extent=[50,-20; 70,0]);
      annotation (Diagram);
      parameter Modelica.SIunits.Temperature T = 343 "Electrode temperature";
      parameter Modelica.SIunits.Concentration c = 500 "Methanol concentration";
      
    equation 
      anode.c = c;
      anode.T = T;
      
      connect(anode.n, ground.p) annotation (points=[-1.11022e-15,-9.4369e-16; 
            31,-9.4369e-16; 31,5.55112e-16; 60,5.55112e-16], style(color=3,
            rgbcolor={0,0,255}));
    end KrewerAnodeSteadyTest;
  end Test;
end Electrodes;
