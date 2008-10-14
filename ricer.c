#include <assert.h>
#include <stdlib.h>

/* Set the DEBUG macro variable to activate the debug messages. If you also set
 * the VERBOSE macro variable, additional messages will be displayed.
 * The debug() and debug_v() macros are used pretty much as printf().
 * For debugging, compile with:
 *   gcc -DDEBUG -ggdb -o test ricer.c
 * or
 *   gcc -DDEBUG -DVERBOSE -ggdb -o test ricer.c
 * and run "./test".
 */
#if DEBUG
        #include <stdio.h>
        #define debug(...) fprintf (stderr, __VA_ARGS__)
        #if VERBOSE
                #define debug_v(...) debug(__VA_ARGS__)
        #else
                #define debug_v(...) ;
        #endif
#else
        #define debug(...) ;
        #define debug_v(...) ;
#endif

int nonZeroComponents( const double* const z_start,
                       const double* const K_start,
                       const size_t n_start,
                       double** const z,
                       double** const K,
                       size_t* const n );

int findBrackets( const double* const K,
                  const size_t n,
                  double* const beta_max,
                  double* const beta_min );

void findFiniteInterval( const double beta_max,
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
 * Note that K_i is defined, against convention, as x_i = K_i*y_i. */
double ricer( double* z_start, size_t nz, double* K_start, size_t nK )
{
        assert( nz == nK );
        size_t n = 0;
        double* z = 0;
        double* K = 0;

        debug("Starting research of non-zero components...\n");
        switch( nonZeroComponents(z_start, K_start, nz, &z, &K, &n) ) {
        case 1: /* No memory freed, since none has been allocated. */
                return 0.0;
        case 0: /* Everything is well, proceed. */
                break;
        default: /* Unhandled case: something went really wrong. */
                exit(1);
        }
        /* If there is only one component, the Rachford-Rice equation does not
         * have any solution. A default value for beta is then returned; if
         * the only component's K value is less than one, it will be in liquid
         * phase; if it is larger than one, in gas phase; if it is exactly one,
         * it will split evenly between gas and liquid phase (though this is
         * really a freakish limit case). */
        if( n == 1 ) {
                double default_beta;
                if( K[0] > 1.0 )
                        default_beta = 0.0;
                else if( K[0] < 1.0 )
                        default_beta = 1.0;
                else
                        default_beta = 0.5;

                free( z );
                free( K );
                return default_beta;
        }
        debug("\tcomplete, %d non-zero components found.\n", n);

        double beta_max, beta_min;
        debug("Looking for asymptotes delimiting the physical solution...\n");
        switch( findBrackets( K, n, &beta_max, &beta_min ) ) {
        case 0: /* Everything is well, proceed. */
                break;
        default: /* Unhandled case: something went really wrong. */
                exit(1);
        }
        if( beta_max == beta_min ) {
                free( z );
                free( K );
                return beta_max;
        }
        debug("\tfound. They are %f and %f.\n", beta_min, beta_max);

        double a, b; // Interval extremes
        debug("Looking for a finite interval for the bisection algorithm...\n");
        findFiniteInterval( beta_max, beta_min, &a, &b, z, K, n );
        if( a == b ) {
                free( z );
                free( K );
                return a;
        }
        debug("\tfound. The interval is between %f and %f\n", a, b);

        double reltol = 1E-9;
        debug("Starting bisection algorithm...\n");
        const double beta = bisection( a, b, z, K, n, reltol );
        debug("\tfinished. Solution is %f, h(beta) = %f.\n", beta,
                                                           h(beta, z, K, n) );

        free(z);
        free(K);

        return beta;
}

/* This function finds out how many components of the original vector are to be
 * kept, and allocates sub-vectors z and K, together with a count n of the
 * valid elements.
 * It returns error codes:
 * - 0: if calculation completed successfully;
 * - 1: if there is no component with non-zero molar fraction. */
int nonZeroComponents( const double* const z_start,
                       const double* const K_start,
                       const size_t n_start,
                       double** const z,
                       double** const K,
                       size_t* const n )
{
        size_t nonZero, i;
        for( i = 0, nonZero = 0; i < n_start; i++ ) {
                if (z_start[i] != 0.0) {
                        debug_v("\tNon-zero component: %d, with z=%f, K=%f\n",
                                                    i, z_start[i], K_start[i] );
                        nonZero++;
                }
        }
        if( nonZero == 0 ) // Error code 1: no components available.
                return 1;

        *z = (double*)malloc(nonZero*sizeof(double));
        assert( *z );
        *K = (double*)malloc(nonZero*sizeof(double));
        assert( *K );

        for( i = 0, *n = 0; i < n_start; i++ ) {
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
 * No error codes are returned. */
int findBrackets( const double* const K,
                  const size_t n,
                  double* const beta_max,
                  double* const beta_min )
{
        /* We need at least two components to find a solution to the
         * Rachford-Rice relation. */
        assert( n >= 2 );
        double K_min = K[0];
        double K_max = K[0];
        size_t i;
        for( i = 1; i < n; i++ ) {
                if (K[i] < K_min)
                        K_min = K[i];
                if (K[i] > K_max)
                        K_max = K[i];
        }

        *beta_min = K_min / (K_min - 1);
        *beta_max = K_max / (K_max - 1);

        /* There is a chance that the betas will be the wrong way around. This
         * happens either when K_min > 1 or K_max < 1. In both of these cases
         * there can be only one phase, as sum(x_i) = 1 and sum(K_i x_i) = 1
         * cannot be satisfied at the same time by any combination of x_i: in
         * other words, the resulting beta will undoubtedly be either < 0 or
         * > 1 (no phase equilibrium). */
        if( *beta_min > *beta_max ) {
                double tmp = *beta_min;
                *beta_min = *beta_max;
                *beta_max = tmp;
        }

        return 0;
}

/* Given beta_max and beta_min (at which h(beta) is infinite, and cannot
 * therefore be evaluated), this function finds an internal interval [a, b]
 * so that beta_min < a < b < beta_max, and h(a)*h(b) < 0. */
void findFiniteInterval( const double beta_max,
                         const double beta_min,
                         double* const a,
                         double* const b,
                         const double* const z,
                         const double* const K,
                         const size_t n )
{
        assert( beta_max > beta_min );
        *a = (beta_max+beta_min) / 2.0;
        const double target_beta = h(*a, z, K, n) < 0 ? beta_max : beta_min;
        *b = (target_beta + *a) / 2.0;

        debug_v("findFiniteInterval: initial values %f and %f\n", *a, *b);
        while( h(*b, z, K, n) * h(*a, z, K, n) >= 0 ) {
                debug_v("\tCycling: h(%f) = %f, h(%f) = %f.\n",
                                       *a, h(*a, z, K, n), *b, h(*b, z, K, n) );
                debug_v("\tChanging b value from %f to %f.\n", *b,
                                                     (target_beta + *b) / 2.0 );
                *b = (target_beta + *b) / 2.0;
        }

        if( *a > *b ) {
                debug_v("findFiniteInterval: a (%f) and b (%f) are about to be "
                        "swapped.\n", *a, *b );
                double temp = *a;
                *a = *b;
                *b = temp;
        }
}

/* This function performs a bisection method on function h() starting from
 * extremes a and b (a<b). It proceeds to refine the interval until it is
 * smaller than the specified tolerance. */
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
 * relation to hold. */
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
/* The main function will test whether ricer.c can work. */
int main()
{
        const size_t n = 5;
        double z[5] = { 0.2, 0.2, 0.0, 0.05, 0.0 };
        double K[5] = { 5.3, 31.2, 0.0, 0.0, 0.0 };
        ricer( z, n, K, n );
        return 0;
}
#endif
