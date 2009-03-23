package Units "Collection of additional units" 
  
  annotation (uses(Modelica(version="2.2.1")));
  
  type ArealCapacitance = Real (final quantity="Areal capacitance", final unit="F/m2");
  type ArealReactionRate = Real(final quantity="Areal reaction rate",
                                final unit="mol/(m2.s)");
  type ArealResistance = Real (final quantity="Areal resistance", final unit="Ohm.m2");
  type CatalystCoverage = Real(final quantity="Catalyst coverage", final unit="", min=0, max=1);
  type MassTransportCoefficient = Real (final quantity="Mass transport coefficient",
                                        final unit="m/s");
  type MolarEnthalpy = Modelica.SIunits.MolarInternalEnergy;
  type MolarFlow = Real(final quantity="Molar flow rate", final unit="mol/s");
  type MolarFlux = Real(final quantity="Molar flux", final unit="mol/(m2.s)");
  type SurfaceConcentration = Real (final quantity="Surface concentration",
                                    final unit="mol/m2");
  
end Units;
