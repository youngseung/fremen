#!/usr/bin/python3
#
# © Federico Zenith, 2008-2013.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#

# from enum import Enum # From Python 3.4
from math import log, sinh, cosh, exp, sqrt

# Define (pseudo-)enumerations
PHASES         = frozenset( ['liq', 'gas'] )
SPECIES        = frozenset( ['CH3OH', 'H2O', 'O2', 'CO2', 'N2'] )
CONDENSABLES   = frozenset( ['CH3OH', 'H2O'] )
INCONDENSABLES = frozenset( ['O2', 'CO2', 'N2'] )

class _Shomate:
  '''Class that packs the parameters of the Shomate equation and returns
enthalpy [J/mol], heat capacity [J/mol K] and entropy [J/mol K] according
to the same equation.'''
  def __init__(self, A, B, C, D, E, F, G, H):
    self._A = A
    self._B = B
    self._C = C
    self._D = D
    self._E = E
    self._F = F
    self._G = G
    self._H = H

  def h( self, T ):
    '''Enthalpy [J/mol] according to Shomate equation, temperature in K;
reference is set to 298.15 K.'''
    t = T/1000
    h = self._A*t + self._B*t**2/2 + self._C*t**3/3 + self._D*t**4/4 - \
        self._E/t + self._F - self._H
    return h*1000

  def cp( self, T ):
    '''Heat capacity [J/mol K] according to Shomate equation,
temperature in K.'''
    t = T/1000
    cp = self._A + self._B*t + self._C*t**2 + self._D*t**3 + self._E/t**2
    return cp
  
  def s( self, T ):
    '''Entropy [J/mol K] according to Shomate equation, temperature in K.'''
    t = T/1000
    s = self._A*log(t) + self._B*t + self._C*t**2/2 + self._D*t**3/3 - \
        self._E/t**2/2 + self._G
    return s

_SHOMATE = dict()
# NB: LIQUID water!
# http://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Units=SI&Mask=2#Thermo-Condensed
_SHOMATE['H2O'] = _Shomate( -203.6060, 1523.290, -3196.413, 2474.455,
                            3.855326, -256.5478, -488.7163, -285.8304 )
#http://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447&Units=SI&Mask=1#Thermo-Gas
_SHOMATE['O2']  = _Shomate( 29.65900, 6.137261, -1.186521, 0.095780,
                            -0.219663, -9.861391, 237.9480, 0.000000 )
#http://webbook.nist.gov/cgi/cbook.cgi?ID=C124389&Units=SI&Mask=1#Thermo-Gas
_SHOMATE['CO2'] = _Shomate( 24.99735, 55.18696, -33.69137, 7.948387,
                            -0.136638, -403.6075, 228.2431, -393.5224 )
#http://webbook.nist.gov/cgi/cbook.cgi?ID=C7727379&Units=SI&Mask=1#Thermo-Gas
_SHOMATE['N2']  = _Shomate( 26.09200, 8.218801, -1.976141, 0.159274,
                            0.044434, -7.989230, 221.0200, 0.000000 )

# Parameters for water properties in gas phase
# Perry's Chemical Engineers' Handbook, 7th edition, page 2-163
# Note: parameter D has been added to those from Perry,
# so that enthalpy would be 0 at 298 K (consistently with library)
_WG = {'A': 8.22, 'B': 0.00015, 'C': 0.00000134, 'D': 10326.2823}

def _cp_H2O_gas( T ):
  '''Returns the specific heat [J/mol K] of water in gas phase,
temperature in K.'''
  return (_WG['A'] + _WG['B']*T + _WG['C']*T**2)*4.184

def _h_H2O_gas( T ):
  '''Returns the enthalpy [J/mol] of water in gas phase, temperature in K.'''
  return (_WG['A']*T + _WG['B']*T**2/2 + _WG['C']*T**3/3)*4.184 - _WG['D']

# Parameters for methanol properties in gas phase
# Source: Perry's Chemical Engineers' Handbook, 7th edition, page 2-179
_MG = (0.3925E5, 0.879E5, 1.9165E3, 0.5365E5, 896.7)

def _cp_CH3OH_gas( T ):
  '''Returns the specific heat [J/mol K] of methanol in gas phase,
temperature in K.'''
  return (_MG[0] + _MG[1]*(_MG[2]/T/(sinh(_MG[2]/T)))**2 \
                 + _MG[3]*(_MG[4]/T/cosh(_MG[4]/T))**2 ) / 1000;

def _h_CH3OH_gas( T ):
  '''Returns the enthalpy [J/mol] of methanol in gas phase, temperature in K.'''
  return (_MG[0]*T + 4*_MG[1]*_MG[2]**2/(2*_MG[2]*exp(2*_MG[2]/T)-2*_MG[2]) + \
                     4*_MG[3]*_MG[4]**2/(2*_MG[4]*exp(2*_MG[4]/T)+2*_MG[4])) \
          / 1000;

# Parameters for methanol properties in liquid phase
# Source: Perry's Chemical Engineers' Handbook, 7th edition, page 2-171
_ML = {'A': 1.058E5, 'B': -362.23, 'C': 0.9379, 'D': 23730.2384}

def _cp_CH3OH_liq( T ):
  '''Returns the specific heat [J/mol K] of methanol in liquid phase,
temperature in K.'''
  return (_ML['A'] + _ML['B']*T + _ML['C']*T**2)/1000

def _h_CH3OH_liq( T ):
  '''Returns the enthalpy [J/mol] of methanol in liquid phase,
temperature in K.'''
  return (_ML['A']*T + _ML['B']*T**2/2 + _ML['C']*T**3/3)/1000 - _ML['D']

def mw( species ):
  'Returns the molecular weight of the species in kg/mol.'
  MW = {'CH3OH': 32.0419e-3,
        'H2O': 18.0153e-3,
        'O2': 31.9988e-3,
        'CO2': 44.0095e-3,
        'N2': 28.01348e-3}
  return MW[species]

# Standard enthalpies of formations of all species and phases
# http://webbook.nist.gov/cgi/cbook.cgi?ID=C67561&Units=SI&Mask=2#Thermo-Condensed
# http://webbook.nist.gov/cgi/cbook.cgi?ID=C67561&Units=SI&Mask=1#Thermo-Gas
# http://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Units=SI&Mask=2#Thermo-Condensed
# http://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Units=SI&Mask=1#Thermo-Gas
# http://webbook.nist.gov/cgi/cbook.cgi?ID=C124389&Units=SI&Mask=1#Thermo-Gas
_DHF = {'CH3OH': {'gas': -205000, 'liq': -238400},
        'H2O':    {'gas': -241826, 'liq': -285830},
        'O2':     {'gas': 0},
        'CO2':    {'gas': -393510},
        'N2':     {'gas': 0}}

def cp( species, phase, T ):
  '''The heat capacity of a species in the given phase at temperature T.
Units for temperature are K, for heat capacity J/mol K.'''
  # Sanity check: no data for incondensables in liquid phase
  assert( not( species in INCONDENSABLES and phase == 'liq' ) )
  if species in INCONDENSABLES or (species == 'H2O' and phase == 'liq' ):
    return _SHOMATE[species].cp( T )
  elif species == 'H2O' and phase == 'gas':
    return _cp_H2O_gas( T )
  elif species == 'CH3OH':
    if phase == 'liq':
      return _cp_CH3OH_liq( T )
    elif phase == 'gas':
      return _cp_CH3OH_gas( T )
  raise KeyError('Bad key: {0}_{1}'.format(species, phase))

def h( species, phase, T ):
  '''The enthalpy of a species in the given phase at temperature T.
Units for temperature are K, for enthalpy J/mol.
Reference is set to standard state and elements in their native form.'''
  T0 = 298.15
  if species in INCONDENSABLES or (species == 'H2O' and phase == 'liq' ):
    f = _SHOMATE[species].h
  elif species == 'H2O' and phase == 'gas':
    f = _h_H2O_gas
  elif species == 'CH3OH':
    if phase == 'liq':
      f = _h_CH3OH_liq
    elif phase == 'gas':
      f = _h_CH3OH_gas
  return f( T ) - f( T0 ) + _DHF[species][phase]

def _rho_H2O( T ):
  '''Density of liquid water, T in K and ρ in kg/m³.
The data has been interpolated with a cubic function, and deviates from
data provided by NIST by at most 0.25 kg/m³.'''
  a = 252.37
  b = 6.5968
  c = -18.288E-3
  d = 15.222E-6
  return a + b*T + c*T**2 + d*T**3

def _rho_CH3OH( T ):
  '''Density of liquid methanol, T in K and ρ in kg/m³.
The data has been interpolated with a linear function, and deviates from
data provided by NIST by at most 0.5 kg/m³.'''
  a = 1069.1
  b = -0.9488
  return a + b*T

def rho( species, phase, T ):
  '''Density of a species in a phase at a temperature in K.
If gas phase, ideal gas law is used. Density is returned as kg/m³.'''
  if phase == 'gas':
    p_env = 101325
    R = 8.314
    return mw(species)*p_env/R/T
  if species == 'CH3OH':
    return _rho_CH3OH( T )
  if species == 'H2O':
    return _rho_H2O( T )
  raise KeyError('Bad key: {0}_{1}'.format(species, phase)) 

def _AntoineLaw( A, B, C ):
  '''Returns a function of T (in K) that produces the partial pressure of a
liquid (in Pa).'''
  return lambda T: 10.0**(A - B/(T + C))*1e5;

def p_vap( species, T ):
  '''Returns the partial pressure in Pa of a liquid species at temperature T
in K.
Data from:
  CH₃OH: http://webbook.nist.gov/cgi/cbook.cgi?ID=C67561&Mask=4#Thermo-Phase
   (Ambrose and Sprake, 1970)
  H₂O: http://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Mask=4#Thermo-Phase
    (Stull, 1947)'''
  f = {'CH3OH': _AntoineLaw(5.20409, 1581.341, -33.5),
       'H2O':   _AntoineLaw(4.65430, 1435.264, -64.848)}
  return f[species](T)

def K( species, T ):
  'Liquid-vapour equilibrium constant for a species at temperature T in K.'
  p_env = 101325
  K = p_vap( species, T )/p_env
  if species == 'CH3OH':
    K *= 2.15 # Methanol activity coefficient γ for diluted solutions
  return K

def rr( zm, zw, T ):
  '''Calculates the mole fraction of vapour β for a mixture of methanol,
water, and incondensables at temperature T in K.
zm and zw are the overall molar fractions of methanol and water respectively.
It is assumed throughout that methanol is diluted.

Whitson, Curtis H., and Michelsen, Michael L.:
The negative flash, Fluid Phase Equilibria 53, 51-71, 1989.'''
  Cm = K('CH3OH', T) - 1
  Cw = K('H2O',   T) - 1
  assert( Cm > Cw )
  a = Cm*Cw
  b = Cm*zm + Cw*zw + (Cm+Cw)*(1-zw-zm)
  c = 1 - zw - zm
  delta = b**2 - 4*a*c

  if Cm == 0:
    return - c/b
  return ( -b - sqrt(delta) ) / (2*a)

if __name__ == '__main__':
  import matplotlib.pyplot as plt

  x = list( range( 274, 330 ) )

  plt.figure()
  for i in SPECIES:
    plt.plot( x, list(map( lambda T: h(i, 'gas', T), x )),
              label='{0}, gas'.format(i) )
  for i in CONDENSABLES:
    plt.plot( x, list(map( lambda T: h(i, 'liq', T), x )),
              label='{0}, liquid'.format(i) )
  plt.title('Enthalpies')
  plt.grid()
  plt.legend()

  plt.figure()
  for i in SPECIES:
    plt.plot( x, list(map( lambda T: cp(i, 'gas', T), x )),
              label='{0}, gas'.format(i) )
  for i in CONDENSABLES:
    plt.plot( x, list(map( lambda T: cp(i, 'liq', T), x )),
              label='{0}, liquid'.format(i) )
  plt.title('Heat capacities')
  plt.grid()
  plt.legend()

  plt.figure()
  for i in SPECIES:
    plt.semilogy( x, list(map( lambda T: rho(i, 'gas', T), x )),
                  label='{0}, gas'.format(i) )
  for i in CONDENSABLES:
    plt.semilogy( x, list(map( lambda T: rho(i, 'liq', T), x )),
                  label='{0}, liquid'.format(i) )
  plt.title('Densities')
  plt.grid()
  plt.legend()

  plt.figure()
  for i in CONDENSABLES:
    plt.plot( x, list(map( lambda T: p_vap(i, T), x )),
              label='{0}'.format(i) )
  plt.title('Vapour pressures')
  plt.grid()
  plt.legend()

  plt.figure()
  for i in CONDENSABLES:
    plt.plot( x, list(map( lambda T: K(i, T), x )),
              label='{0}'.format(i) )
  plt.title('Equilibrium constants')
  plt.grid()
  plt.legend()

  plt.show()

  # TODO RR test: exercise for 3D animation?
