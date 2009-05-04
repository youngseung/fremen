package Units "Collection of additional units" 
  
  annotation (uses(Modelica(version="2.2.1")), Documentation(info="<html>
<p>This package is a list of unit types that are useful in various parts
of our modelling effort. They are mostly self-documenting.</p>
</html>"));
  
  type ArealCapacitance = Real (final quantity="Areal capacitance", final unit="F/m2");
  type ArealReactionRate = Real(final quantity="Areal reaction rate",
                                final unit="mol/(m2.s)");
  type ArealResistance = Real (final quantity="Areal resistance", final unit="Ohm.m2");
  type CatalystCoverage = Real(final quantity="Catalyst coverage", final unit="", min=0, max=1);
  type CondensationCoefficient = Real(final quantity="Condensation coefficient",
                                      final unit="1/s", min=0);
  type DynamicViscosity = Real (final quantity="Dynamic viscosity", final unit="kg/ms", min=0);
  type HeatTransferCoefficient = Real(final quantity="Heat transfer coefficient",
                                      final unit="W/(m2.K)", min=0);
  type MassTransportCoefficient = Real(final quantity="Mass transport coefficient",
                                       final unit="m/s", min=0);
  type MolarEnthalpy = Modelica.SIunits.MolarInternalEnergy;
  type MolarFlow = Real(final quantity="Molar flow rate", final unit="mol/s");
  type MolarFlux = Real(final quantity="Molar flux", final unit="mol/(m2.s)");
  type Permeability = Real (final quantity="Permeability", final unit="m2", min=0);
  type Porosity = Real (final quantity="Porosity", final unit="", min=0, max=1);
  type ReactionRate = Real(final quantity="Reaction rate", final unit="mol/(m3.s)");
  type Saturation = Real (final quantity="Saturation", final unit="", min=0, max=1);
  type SurfaceConcentration = Real (final quantity="Surface concentration",
                                    final unit="mol/m2");
  
  type Molecule = Integer (final quantity="Molecule identifier", final min=1, final max=5) 
    "Identifier for various types of molecules";
  type Phase = Integer (final quantity="Phase identifier", final min=1001, final max=1002) 
    "Identifier for phases";
  
end Units;
