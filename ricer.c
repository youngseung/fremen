#include <assert.h>
#include <float.h>
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

#define max(a, b)       ((a) > (b) ? (a) : (b))
#define min(a, b)       ((a) < (b) ? (a) : (b))

int countComponents( const double* const z,
                     const size_t nz,
                     const double* const K,
                     const double z_g,
                     double* const default_beta );

void findAsymptotes( const double* const z,
                     const double* const K,
                     const size_t n,
                     const double z_g,
                     double* const beta_max,
                     double* const beta_min );

void findFiniteInterval( const double beta_min,
                         const double beta_max,
                         double* const a,
                         double* const b,
                         const double* const z,
                         const double* const K,
                         const size_t n,
                         const double z_g );

double bisection( double a,
                  double b,
                  const double* const z,
                  const double* const K,
                  const size_t n,
                  const double z_g,
                  const double reltol );

double h( const double beta,
          const double* const z,
          const double* const K,
          const size_t n,
          const double z_g );

/* This returns beta = n_g/n_tot, given z_i and K_i.
 * Note that K_i is defined, against convention, as x_i = K_i*y_i.
 * The provided compositions can be zero (for numerical noise, they could also
 * be slighly negative).
 */
double ricer( double* z, size_t nz, double* K, size_t nK, double z_g )
{
        assert( nz == nK );

        double single_beta;
        switch( countComponents(z, nz, K, z_g, &single_beta) ) {
        case 0: /* Fed bad data. Cannot find any sensible solution. Crash. */
                abort();
        case 1: /* Single component, return simple guess. */
                return single_beta;
        default: /* Normal situation, proceed. */
                break;
        }


        double beta_max, beta_min;
        debug("Looking for asymptotes delimiting the physical solution...\n");
        findAsymptotes( z, K, nz, z_g, &beta_max, &beta_min );
        /* This covers the case of only one component being present. */
        if( beta_min == beta_max )
                return beta_min;
        debug("\tfound. They are %f and %f.\n", beta_min, beta_max);

        double a, b; /* Interval extremes */
        debug("Looking for a finite interval for the bisection algorithm...\n");
        findFiniteInterval( beta_min, beta_max, &a, &b, z, K, nz, z_g );
        if( a == b )
                return a;
        debug("\tfound. The interval is between %f and %f\n", a, b);

        double reltol = 1E-9;
        debug("Starting bisection algorithm...\n");
        const double beta = bisection( a, b, z, K, nz, z_g, reltol );
        debug("\tfinished. Solution is %f, h(beta) = %e.\n", beta,
                                                       h(beta, z, K, nz, z_g) );

        return beta;
}

/* This function makes a preliminary count of the components provided to the
 * Rachford-Rice solution algorithm. If there is only one component set, the
 * variable pointed by single_beta will be set to zero or one, depending on the
 * K-value (smaller or larger than 1).
 * If the K-value is exactly 1, 0.5 will be returned, but that should strictly
 * speaking be an undetermined case; furthermore, in our model we do not expect
 * to get single-components mixtures of water or methanol at boiling
 * temperature.
 */
int countComponents( const double* const z,
                     const size_t n,
                     const double* const K,
                     const double z_g,
                     double* const single_beta )
{
        int nonZero = 0;
        size_t i;
        for( i = 0; i < n; i++ ) {
                if( z[i] > 0.0 ) {
                        nonZero++;
                        if( K[i] < 1.0 )
                                *single_beta = 0.0;
                        else if ( K[i] == 1.0 )
                                *single_beta = 0.5;
                        else
                                *single_beta = 1.0;
                }
        }
        if( z_g > 0.0 ) {
                nonZero++;
                *single_beta = 1.0;
        }
        return nonZero;
}

/* Given the vector of z_i and K_i, and the number of gas moles, this function
 * sets the pointers to beta_max and beta_min accordingly.
 * Components with zero composition are not considered, as they do not influence
 * the shape of the Rachford-Rice function.
 * If gas components are present, they are considered later on, since their
 * K-value would be infinite, while their corresponding value of beta is 1.0,
 * and is therefore manageable.
 *
 * No error codes are returned. */
void findAsymptotes( const double* const z,
                     const double* const K,
                     const size_t n,
                     const double z_g,
                     double* const beta_max,
                     double* const beta_min )
{// FIXME why do I get equal betas with water & gas?
        double K_min = DBL_MAX;
        double K_max = 0.0;

        size_t i;
        for( i = 0; i < n; i++ ) {
                if( z[i] > 0.0 ) {
                        K_min = min( K_min, K[i] );
                        K_max = max( K_max, K[i] );
                }
        }

        *beta_min = 1.0 / (1.0 - K_max);
        *beta_max = 1.0 / (1.0 - K_min);

        if( z_g > 0.0 )
                *beta_min = max( *beta_min, 0.0 );

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
}

/* Given beta_max and beta_min (at which h(beta) is infinite, and cannot
 * therefore be evaluated), this function finds an internal interval [a, b]
 * so that beta_min < a < b < beta_max, and h(a)*h(b) < 0.
 * First, the middle point between the two betas is taken. Then, it is checked
 * whether the value of h(beta) is positive or negative at that point, and a
 * target asymptote (i.e. the asymptote on the other side with respect to the
 * solution) is set. The algorithm then approaches the asymptote gradually by
 * bisection, until a point with different sign and finite value of h(beta) is
 * found.
 * Finally, if necessary, a and b are swapped (the algorithm gives have no
 * guarantee of them being in any order).
 *
 * The function returns no error codes. */
void findFiniteInterval( const double beta_min,
                         const double beta_max,
                         double* const a,
                         double* const b,
                         const double* const z,
                         const double* const K,
                         const size_t n,
                         const double z_g )
{
        assert( beta_max > beta_min );
        *a = (beta_min+beta_max) / 2.0;
        const double other_side_asymptote = h(*a, z, K, n, z_g) < 0 ? beta_min :
                                                                      beta_max;
        *b = (other_side_asymptote + *a) / 2.0;

        debug_v("findFiniteInterval: initial values %f and %f\n", *a, *b);
        while( h(*b, z, K, n, z_g) * h(*a, z, K, n, z_g) > 0 ) {
                debug_v("\tCycling: h(%f) = %f, h(%f) = %f.\n",
                             *a, h(*a, z, K, n, z_g), *b, h(*b, z, K, n, z_g) );
                debug_v("\tChanging b value from %f to %f.\n", *b,
                                            (other_side_asymptote + *b) / 2.0 );
                *b = (other_side_asymptote + *b) / 2.0;
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
                  const double z_g,
                  const double reltol )
{
        assert( h(a, z, K, n, z_g)*h(b, z, K, n, z_g) <= 0.0 );
        assert( b > a );

        double h_a = h(a, z, K, n, z_g);
        if( h_a == 0.0 )
                return a;
        double h_b = h(b, z, K, n, z_g);
        if( h_b == 0.0 )
                return b;

        double middle, h_middle;
        const double abstol = reltol*(b-a);
        while( b-a > abstol ) {
                middle = (b+a)/2.0;
                h_middle = h(middle, z, K, n, z_g);
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

/* This function is the one that has to be zero for the Rachford-Rice relation
 * to hold. Components with zero fraction are disregarded. */
double h( const double beta,
          const double* const z,
          const double* const K,
          const size_t n,
          const double z_g )
{
        double sum = 0.0;
        size_t i;
        for( i = 0; i < n; i++ ) {
                if( z[i] > 0.0 )
                        sum += z[i]*(K[i] - 1.0) / (1.0 + beta*(K[i]-1.0));
        }
        if( z_g > 0.0 )
                sum += z_g / beta;
        return sum;
}

#if DEBUG
/* The main function will test whether ricer.c can work. */
int main( int argc, char* argv[] )
{
        if( argc != 3 ) {
                printf("Usage: %s z_h2o z_ch3oh\n", argv[0]);
                return 1;
        }

        const size_t n = 2;
        double z[2] = { atof(argv[1]), atof(argv[2]) };
        double K[2] = { 0.1672, 0.031378 };
        double z_g = 1.0 - atof(argv[1]) - atof(argv[2]);
        debug("Beta: %f\n", ricer( z, n, K, n, z_g ) );
        return 0;
}
#endif
