#ifndef LINEAR_INTERPOLATION_H
#define LINEAR_INTERPOLATION_H

#include <assert.h>

/**
 * Interpolates data series data_x and data_y for a given point x.
 */
double linear_interpolation( const double x, const double* data_x,
                             const double* data_y, const int size )
{
	assert( size >= 2 );
	int i;
	for( i = 0; i < size-1; ++i )
		assert( data_x[i] < data_x[i+1] );

	if( x < data_x[0] )
		return data_y[0];
	if( x >= data_x[size-1] )
		return data_y[size-1];

	// Find the index of the first value inferior to x in the array
	int pos = 0;
	while( x > data_x[pos+1])
		pos++;
	assert( pos < size - 1 );

	// Linearly interpolate the value of y
	const double f = ( x - data_x[pos] ) / ( data_x[pos+1] - data_x[pos] );
	return ( data_y[pos] + f * ( data_y[pos+1] - data_y[pos] ) );
}

#endif
