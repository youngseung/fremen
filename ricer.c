#include <assert.h>
#include <stdlib.h>

int nonZeroComponents( const double* const z_start,
                       const double* const K_start,
                       const size_t n_start,
                       double** const z,
                       double** const K,
                       size_t* const n );

void findBrackets( const double* const K,
                   const size_t n,
                   double* const beta_min,
                   double* const beta_max );

void findInterval( const double beta_max,
                   const double beta_min,
                   double* const a,
                   double* const b,
                   const double* const z,
                   const double* const K,
                   const size_t n );

double bisection( double a,
                  double b,
                  const double* const z,
                  const double* const K,
                  const size_t n,
                  const double reltol );

double h( const double beta,
          const double* const z,
          const double* const K,
          const size_t n );

/* This returns beta = n_g/n_tot, given z_i and K_i.
 * Note that K_i is defined, against convention, as x_i = K_i*y_i.
 */
double ricer( double* z_start, size_t nz, double* K_start, size_t nK )
{
	assert( nz == nK );
	size_t n = 0;
	double* z = 0;
	double* K = 0;
	nonZeroComponents(z_start, K_start, nz, &z, &K, &n);

	double beta_max, beta_min;
	findBrackets( K, n, &beta_max, &beta_min );

	if( beta_max == beta_min )
		return beta_max;

	double a, b; // Interval extremes
	findInterval( beta_max, beta_min, &a, &b, z, K, n );
	if( a == b )
		return a;

	double reltol = 1E-9;
	const double beta = bisection( a, b, z, K, n, reltol );
	free(z);
	free(K);
	return beta;
}

/* This function finds out how many components of the original vector are to be
 * kept, and allocates sub-vectors z and K, together with a count n of the
 * valid elements.
 */
int nonZeroComponents( const double* const z_start,
                       const double* const K_start,
                       const size_t n_start,
                       double** const z,
                       double** const K,
                       size_t* const n )
{
	size_t nonZero = 0;
	size_t i;

	for( i = 0; i < n_start; i++ ) {
		if (z_start[i] > 0)
			nonZero++;
	}

	if( nonZero == 0 )
		return 1; // Error code: no composition > 0

	*z = (double*)malloc(nonZero*sizeof(double));
	assert( *z );
	*K = (double*)malloc(nonZero*sizeof(double));
	assert( *K );
	*n = 0;

	for( i = 0; i < n_start; i++ ) {
		if( z_start[i] > 0.0 ) {
			(*z)[*n] = z_start[i];
			(*K)[*n] = K_start[i];
			(*n)++;
		}
	}

	return 0;
}

/* Given the vector of K_i, this function sets the pointers to beta_max and
 * beta_min accordingly.
 */
void findBrackets( const double* const K,
                   const size_t n,
                   double* const beta_min,
                   double* const beta_max )
{
	assert( n > 0 );
	double K_min = K[0];
	double K_max = K[0];
	size_t i;
	for( i = 1; i < n; i++ ) {
		if (K[i] < K_min)
			K_min = K[i];
		if (K[i] > K_max)
			K_max = K[i];
	}

	*beta_min = K_min / (K_min -1);
	*beta_max = K_max / (K_max -1);
	assert( *beta_max >= *beta_min );
}

/* Given beta_max and beta_min (at which h(beta) is infinite, and cannot
 * therefore be evaluated), this function finds an internal interval [a, b]
 * so that beta_min < a < b < beta_max, and h(a)*h(b) < 0.
 */
void findInterval( const double beta_max,
                   const double beta_min,
                   double* const a,
                   double* const b,
                   const double* const z,
                   const double* const K,
                   const size_t n )
{
	*a = (beta_max+beta_min) / 2.0;
	const double beta_to_approach = h(*a, z, K, n) < 0 ? beta_min : beta_max;
	*b = (beta_to_approach+*a) / 2.0;
	while( h(*b, z, K, n) * h(*a, z, K, n) >= 0 )
		*b = (beta_to_approach+*b) / 2.0;

	if( *a > *b ) {
		double temp = *a;
		*a = *b;
		*b = temp;
	}
}

/* This function performs a bisection method on function h() starting from
 * extremes a and b (a<b). It proceeds to refine the interval until it is
 * smaller than the specified tolerance.
 */
double bisection( double a,
                  double b,
                  const double* const z,
                  const double* const K,
                  const size_t n,
                  const double reltol )
{
	assert( h(a, z, K, n)*h(b, z, K, n) <= 0.0 );
	assert( b > a );

	double h_a = h(a, z, K, n);
	if( h_a == 0.0 )
		return a;
	double h_b = h(b, z, K, n);
	if( h_b == 0.0 )
		return b;

	double middle, h_middle;
	const double abstol = reltol*(b-a);
	while( b-a > abstol ) {
		middle = (b+a)/2.0;
		h_middle = h(middle, z, K, n);
		if( h_middle == 0.0 )
			return middle;
		else if( h_middle*h_a > 0.0 ) {
			a = middle;
			h_a = h_middle;
		} else {
			b = middle;
			h_b = h_middle;
		}
	}

	return (b+a)/2.0;
}

/* This function is the one that has to be zero for the modified Rachford-Rice
 * relation to hold.
 */
double h( const double beta,
          const double* const z,
          const double* const K,
          const size_t n )
{
	double sum = 0.0;
	size_t i;
	for( i = 0; i < n; i++ )
		sum += z[i]*(K[i]-1) / (beta + K[i]*(1-beta));
	return sum;
}

#if DEBUG
/* The main function will test whether ricer.c can work.
 */
int main()
{
	const size_t n = 5;
	double z[5] = { 0.2, 0.2, 0.0, 0.05, 0.0 };
	double K[5] = { 5.3, 31.2, 0.0, 0.0, 0.0 };
	ricer( z, n, K, n );
	return 0;
}
#endif
