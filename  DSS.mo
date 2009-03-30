package DSS "Derivative approximation functions" 
  
  partial model FirstDer "Discretisation of values and first derivatives" 
    parameter Real dx "Step size";
    parameter Integer m "Number of points";
    
    Real[m] v "values";
    Real[m] d "derivatives";
  equation 
    
    annotation (Documentation(info="<html>
<p>This partial class contains a vector of values and a corresponding vector
of approximated derivatives. The precise way of relating values and derivatives
is first established in child classes.</p>
</html>"));
  end FirstDer;
  
  model FirstDer_2 "Second-order approximation" 
    extends FirstDer;
  equation 
    assert( m >= 3, "Not enough discretisation points");
    
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
<p>Schiesser, `The numerical method of lines', 1991, p.100.</p>
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
<p>Schiesser, `The numerical method of lines', 1991, p.108.</p>
</html>"));
  equation 
    assert( m >= 5, "Not enough discretisation points");
    
    d[1]     = (-50*v[1]+96*v[2]-72*v[3]+32*v[4]-6*v[5]) / (24*dx);
    d[2]     = (-6*v[1]-20*v[2]+36*v[3]-12*v[4]+2*v[5]) / (24*dx);
    d[3:m-2] = (2*v[1:m-4]-16*v[2:m-3]+16*v[4:m-1]-2*v[5:m]) / (24*dx);
    d[m-1]   = (-2*v[m-4]+12*v[m-3]-36*v[m-2]+20*v[m-1]+6*v[m]) / (24*dx);
    d[m]     = (6*v[m-4]-32*v[m-3]+72*v[m-2]-96*v[m-1]+50*v[m]) / (24*dx);
    
  end FirstDer_4;
  
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
  end Test;
end DSS;
