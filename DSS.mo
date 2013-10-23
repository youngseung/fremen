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

package DSS "Derivative approximation functions"
  
  partial model FirstDer "Discretisation of values and first derivatives" 
    parameter Real dx "Step size";
    parameter Integer m "Number of points";
    
    Real[m] v "Values";
    Real[m] d "Derivatives";
  equation 
    assert( dx > 0, "==> Bad value for dx: "+String(dx));
    annotation (Documentation(info="<html>
<p>This partial class contains a vector of values and a corresponding vector
of approximated derivatives. The precise way of relating values and derivatives
is first established in child classes.</p>
</html>"));
  end FirstDer;
  
  model FirstDer_2 "Second-order approximation" 
    extends FirstDer;
  equation 
    assert( m >= 3, "==> Not enough discretisation points");
    
    d[1]     = (-3*v[1]+4*v[2]-v[3]) / (2*dx);
    d[2:m-1] = (v[3:m] - v[1:m-2]) / (2*dx);
    d[m]     = (v[m-2]-4*v[m-1]+3*v[m]) / (2*dx);
    annotation (Documentation(info="<html>
<p>This class contains a vector of values and a corresponding vector
of approximated derivatives, estimated with a O(&Delta;x<sup>2</sup>) 
approximation. Discretisation is assumed uniform with step &Delta;x.</p>
<p>This approximation is superaccurate (i.e. actually exact) for
polynomials up to the second degree.</p>
<h3>Reference</h3>
<p>Schiesser, 'The numerical method of lines', 1991, p.100.</p>
</html>"));
  end FirstDer_2;
  
  model FirstDer_4 "Fourth-order approximation" 
    extends FirstDer;
    annotation (Documentation(info="<html>
<p>This class contains a vector of values and a corresponding vector
of approximated derivatives, estimated with a O(&Delta;x<sup>4</sup>) 
approximation. Discretisation is assumed uniform with step &Delta;x.</p>
<p>This approximation is superaccurate (i.e. actually exact) for
polynomials up to the fourth degree.</p>
<h3>Reference</h3>
<p>Schiesser, 'The numerical method of lines', 1991, p.108.</p>
</html>"));
  equation 
    assert( m >= 5, "==> Not enough discretisation points");
    
    d[1]     = (-50*v[1]+96*v[2]-72*v[3]+32*v[4]-6*v[5]) / (24*dx);
    d[2]     = (-6*v[1]-20*v[2]+36*v[3]-12*v[4]+2*v[5]) / (24*dx);
    d[3:m-2] = (2*v[1:m-4]-16*v[2:m-3]+16*v[4:m-1]-2*v[5:m]) / (24*dx);
    d[m-1]   = (-2*v[m-4]+12*v[m-3]-36*v[m-2]+20*v[m-1]+6*v[m]) / (24*dx);
    d[m]     = (6*v[m-4]-32*v[m-3]+72*v[m-2]-96*v[m-1]+50*v[m]) / (24*dx);
    
  end FirstDer_4;
  
  model SecondDer 
    
    parameter Real dx "Step size";
    parameter Integer m "Number of points";
    
    parameter Boolean leftIsNeumann = false 
      "Whether the left boundary condition is of Neumann type";
    parameter Real dl = 0 
      "The value of the left Neumann boundary condition, if applicable";
    parameter Boolean rightIsNeumann = false 
      "Whether the right boundary condition is of Neumann type";
    parameter Real dr = 0 
      "The value of the right Neumann boundary condition, if applicable";
    
    Real[m] v "Values";
    Real[m] dd "Second derivatives";
    annotation (Documentation(info="<html>
<p>This class contains a vector of values and a corresponding vector
of approximated second-order derivatives, estimated with a O(&Delta;x<sup>4</sup>) 
approximation. Discretisation is assumed uniform with step &Delta;x.</p>
<p>This approximation is superaccurate (i.e. actually exact) for
polynomials up to the fourth degree.</p>
<p>This object allows to specify one of the boundary conditions directly as a
Neumann boundary condition, in which case a value has to be provided. At least one 
of the conditions, however, must always be a Dirichlet condition.</p>
<h3>Reference</h3>
<p>Schiesser, 'The numerical method of lines', 1991, p.112-115.</p>
</html>"));
  equation 
    assert( dx > 0, "==> Bad value for dx: "+String(dx));
    assert( m >= 6, "==> Not enough discretisation points");
    assert( not (leftIsNeumann and rightIsNeumann), "==> Cannot have two Neumann conditions");
    
    if leftIsNeumann then
      dd[1] = (-415/3*v[1] + 192*v[2] - 72*v[3] +64/3*v[4] -3*v[5] -100*dl*dx) / (24*dx*dx);
    else
      dd[1] = (90*v[1] - 308*v[2] +428*v[3] -312*v[4] +122*v[5] -20*v[6]) / (24*dx*dx);
    end if;
    dd[2] = (20*v[1] -30*v[2] -8*v[3] +28*v[4] -12*v[5] +2*v[6]) / (24*dx*dx);
    dd[3:(m-3)] = (-2*v[1:(m-5)] +32*v[2:(m-4)] -60*v[3:(m-3)] +32*v[4:(m-2)] -2*v[5:(m-1)]) / (24*dx*dx);
    dd[m-2] = (-2*v[m-4] +32*v[m-3] -60*v[m-2] +32*v[m-1] -2*v[m]) / (24*dx*dx);
    dd[m-1] = (2*v[m-5] -12*v[m-4] +28*v[m-3] -8*v[m-2] -30*v[m-1] +20*v[m]) / (24*dx*dx);
    if rightIsNeumann then
      dd[m] = (-3*v[m-4] +64/3*v[m-3] -72*v[m-2] +192*v[m-1] -415/3*v[m] +100*dr*dx) / (24*dx*dx);
    else
      dd[m] = (-20*v[m-5] +122*v[m-4] -312*v[m-3] +428*v[m-2] -308*v[m-1] +90*v[m]) / (24*dx*dx);
    end if;
    
  end SecondDer;

  package Test 
    
  public 
    model Test_2 
      
      parameter Real[:] v = (1:10).^2;
      FirstDer_2 d(m=10,dx=1);
      Real error = max(2*(1:10) - d.d);
      
    equation 
      d.v = v;
      
    end Test_2;
    
  public 
    model Test_4 
      
      parameter Real[:] v = (1:10).^4;
      FirstDer_4 d(m=10,dx=1);
      Real error = max(4*(1:10).^3 - d.d);
      
    equation 
      d.v = v;
      
    end Test_4;
  public 

    model Test_Der2 
      
      parameter Real[:] v = (1:10).^4;
      SecondDer dd(m=10,dx=1) "Dirichlet - Dirichlet";
      SecondDer dn(m=10,dx=1,rightIsNeumann=true,dr=4000) "Dirichlet - Neumann";
      SecondDer nd(m=10,dx=1,leftIsNeumann=true,dl=4) "Neumann - Dirichlet";
      Real[:] analytic = 12*(1:10).^2;
      Real error_dd = max(abs(analytic - dd.dd));
      Real error_dn = max(abs(analytic - dn.dd));
      Real error_nd = max(abs(analytic - nd.dd));
      Real error = max({error_dd, error_dn, error_nd});
      
    equation 
      dd.v = v;
      dn.v = v;
      nd.v = v;
      
    end Test_Der2;
  end Test;
end DSS;
