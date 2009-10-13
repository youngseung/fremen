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


package Units "Collection of additional units"

  annotation (uses(Modelica(version="3.1")),   Documentation(info="<html>
<p>This package is a list of unit types that are useful in various parts
of our modelling effort. Some constants are also included.</p>
</html>"),
    version="1",
    conversion(noneFromVersion=""));

  type ArealCapacitance = Real (final quantity="Areal capacitance", final unit="F/m2")
    "Used in electrochemical transients";
  type ArealReactionRate = Real(final quantity="Areal reaction rate",
                                final unit="mol/(m2.s)")
    "Important for membrane reaction modelling";
  type ArealResistance = Real (final quantity="Areal resistance", final unit="Ohm.m2")
    "Important property of fuel-cell membranes";
  type CatalystCoverage = Real(final quantity="Catalyst coverage", final unit="", min=0, max=1)
    "Fraction of catalyst occupied by some species";
  type CondensationCoefficient = Real(final quantity="Condensation coefficient",
                                      final unit="1/s", min=0)
    "How fast a liquid-vapour equilibrium is attained";
  type DynamicViscosity = Real (final quantity="Dynamic viscosity", final unit="kg/ms", min=0)
    "Property of fluids";
  constant Modelica.SIunits.FaradayConstant F = 96485.3415
    "The Faraday constant";
  type HeatTransferCoefficient = Real(final quantity="Heat transfer coefficient",
                                      final unit="W/(m2.K)", min=0)
    "Simplification of heat transfer processes";
  type MassTransportCoefficient = Real(final quantity="Mass transport coefficient",
                                       final unit="m/s", min=0)
    "Simplification of diffusion processes";
  type MolarEnthalpy = Modelica.SIunits.MolarInternalEnergy
    "Enthalpy per amount of substance";
  type MolarFlow = Real(final quantity="Molar flow rate", final unit="mol/s")
    "Flow of substance";
  type MolarFlux = Real(final quantity="Molar flux", final unit="mol/(m2.s)")
    "Flow per unit of cross-sectional area";
  type Permeability = Real (final quantity="Permeability", final unit="m2", min=0)
    "Parameter for fluid flow through porous media";
  type Porosity = Real (final quantity="Porosity", final unit="", min=0, max=1)
    "Fraction of pores in a medium";
  type ReactionRate = Real(final quantity="Reaction rate", final unit="mol/(m3.s)")
    "Rate of production of a species in a control volume";
  type RelativeHumidity = Real(final quantity="Relative humidity", final unit="1", min=0, max=100)
    "Percentage of maximum water-vapour quantity in air";
  type Saturation = Real (final quantity="Saturation", final unit="", min=0, max=1)
    "Indicates the liquid saturation in a two-phase setting";
  type SurfaceConcentration = Real (final quantity="Surface concentration",
                                    final unit="mol/m2")
    "Used to express catalyst concentration";

  type Temperature = Modelica.SIunits.Temperature(nominal=298.15, start=298.15)
    "Temperature, with nominal and start values" annotation (Documentation(info=
         "<html>
<p>This type allows to save one the nuisance to specify every time that
the first guess should be environment temperature, not absolute zero.</p>
<p>However, do <em>not</em> use that for parameters, since they will not
be editable any more in the simulation tab.</p>
</html>"));

end Units;
