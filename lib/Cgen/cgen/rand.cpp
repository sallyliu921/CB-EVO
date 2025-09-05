/*--------------------------------------------------------------------------*
 * Copyright 2002 by Paul D. Kundarewich, Michael Hutton, Jonathan Rose     *
 * and the University of Toronto. 											*
 * Use is permitted, provided that this attribution is retained  			*
 * and no part of the code is re-distributed or included in any commercial	*
 * product except by written agreement with the above parties.              *
 *                                                                          *
 * For more information, contact us directly:                               *
 *	  Paul D. Kundarewich (paul.kundarewich@utoronto.ca)					*
 *    Jonathan Rose  (jayar@eecg.toronto.edu)                               *
 *    Mike Hutton  (mdhutton@cs.toronto.edu, mdhutton@eecg.toronto.edu)     *
 *    Department of Electrical and Computer Engineering                     *
 *    University of Toronto, 10 King's College Rd.,                         *
 *    Toronto, Ontario, CANADA M5S 1A4                                      *
 *    Phone: (416) 978-6992  Fax: (416) 971-2286                            *
 *--------------------------------------------------------------------------*/

#ifdef VISUAL_C
#pragma warning(disable:4786)
#endif
 
#include "rand.h"
#include <algorithm>
#include "util.h"
#include "gen.h"
#include <map>
#include <math.h>

#ifdef VISUAL_C
static long g_current_seed = 0L;      // for MSVC substitute
const double M_PI = 3.14159265358979323846;
#endif

RANDOM_NUMBER::RANDOM_NUMBER()
{
	//  Initial state vector taken from random() man page.
	
	long state1[32] = { 3,
    0x9a319039, 0x32d9c024, 0x9b663182, 0x5da1f342,
    0x7449e56b, 0xbeb1dbb0, 0xab5c5918, 0x946554fd,
    0x8c2e680f, 0xeb3d799f, 0xb11ee0b7, 0x2d436b86,
    0xda672e2a, 0x1588ca88, 0xe369735d, 0x904f35f7,
    0xd7158fd6, 0x6fa6f051, 0x616e6b96, 0xac94efdc,
    0xde3b81e0, 0xdf0a6fb5, 0xf103bc02, 0x48f340fb,
    0x36413f93, 0xc622c298, 0xf5a42ab8, 0x8a88d77b,
    0xf5ad9d0e, 0x8999220b, 0x27fb47b9
	};

	int i;
	for (i = 0; i < 32; i++)
	{
		m_state1[i] = state1[i];
	}

	randomize(0);
};

RANDOM_NUMBER::RANDOM_NUMBER(long seed) 
{
	// Initial state vector taken from random() man page.
	
	long state1[32] = { 3,
    0x9a319039, 0x32d9c024, 0x9b663182, 0x5da1f342,
    0x7449e56b, 0xbeb1dbb0, 0xab5c5918, 0x946554fd,
    0x8c2e680f, 0xeb3d799f, 0xb11ee0b7, 0x2d436b86,
    0xda672e2a, 0x1588ca88, 0xe369735d, 0x904f35f7,
    0xd7158fd6, 0x6fa6f051, 0x616e6b96, 0xac94efdc,
    0xde3b81e0, 0xdf0a6fb5, 0xf103bc02, 0x48f340fb,
    0x36413f93, 0xc622c298, 0xf5a42ab8, 0x8a88d77b,
    0xf5ad9d0e, 0x8999220b, 0x27fb47b9
	};

	int i;
	for (i = 0; i < 32; i++)
	{
		m_state1[i] = state1[i];
	}

	randomize(seed);
}

RANDOM_NUMBER::~RANDOM_NUMBER()
{


}
//
//	Seed the random number generator.  Leave seed a parameter in case
//	we want to generate the same graph again.  Note that the seed should
//	always be kept as documentation for a circuit;  seed,n,k and this
//	algorithm totally determine a circuit.
//
//	PRE: seed is the random number seed. 
//	     this seed if 0 if we want to use clock ticks as the random seed
//	POST: the random function has been randomized
//	      m_seed contains the random seed
//	RETURNS: the random seed
//
long RANDOM_NUMBER::randomize(long seed)
{
#ifndef VISUAL_C
    //extern char *initstate();
    //extern char *setstate();
#endif
	m_seed = seed;
    long ticks;

    if (m_seed == 0) 
	{
		ticks = util_ticks();		 // Seed from clock (seconds since Jan 1, 1970)
        debug("Seeding randomize with clock-ticks (" << ticks << ")");
        m_seed = ticks; 
    } 
	else 
	{
        debug("Seeding randomize with supplied seed (" << m_seed << ")");
    }
#ifdef VISUAL_C
    //Warning("Need to fix poor-quality random number generation in MSVC.\n");
    //current_seed = seed;
#else
    initstate(m_seed, (char *) m_state1, 128);
    setstate((char *) m_state1);
#endif

    return m_seed;
}


#ifdef VISUAL_C
//
// a random number function to replace the one in visual c++ that 
// i've been told is bad.
//
extern 
long random(void)
{
    g_current_seed = g_current_seed * 214013L + 2531011L;
    return ((int) ((g_current_seed >>16) & 0x7fff));
}
#endif




// 
//  Return a random permutation on [n].
//  - Initialize 0..n-1 array; pick randomly on [i], i = n-1..0, and
//    swap with the end of the list. 
// 
// PRE: n >= 0
// RETURNS: a random permutation on [n]

NUM_ELEMENTS_VECTOR RANDOM_NUMBER::get_permutation
(
	const NUM_ELEMENTS& n
) const
{
	assert(n >= 0);
    NUM_ELEMENTS_VECTOR array(n);
    int i;
  
    for (i = 0; i < n; i++)
	{
        array[i] = i;
	}

	random_shuffle(array.begin(), array.end());

    return array;
}

// 
//  Choose uniformly from a weighted array.  i.e., with array [3,3,4]
//  return 0 with prob .3, 1 with prob .3 and 2 with prob .4.
//
//  If doing multiple queries, is better to use the cmf version, which 
//  takes only log(n) time; this one is linear.
//
//	RETURNS: a choice from the prob. mass function
NUM_ELEMENTS RANDOM_NUMBER::discrete_pmf(const NUM_ELEMENTS_VECTOR & pmf)
{
    NUM_ELEMENTS n = static_cast<NUM_ELEMENTS>(pmf.size());

	return discrete_pmf(n, pmf);
}

//  Choose uniformly from a weighted array.  i.e., with array [3,3,4]
//  return 0 with prob .3, 1 with prob .3 and 2 with prob .4.
//
//  If doing multiple queries, is better to use the cmf version, which 
//  takes only log(n) time; this one is linear.
//
//	RETURNS: a choice from the prob. mass function
NUM_ELEMENTS RANDOM_NUMBER::discrete_pmf(const NUM_ELEMENTS & n, const NUM_ELEMENTS_VECTOR & pmf)
{
    NUM_ELEMENTS i;
    NUM_ELEMENTS r, tot;
	NUM_ELEMENTS number = 0;
	m_cmf.clear();
	map<NUM_ELEMENTS, NUM_ELEMENTS>::iterator upper_bound_iter;

	tot = 0;
    for (i = 0; i < n; i++)
	{
		number = pmf[i];
		assert(number >= 0);
		tot += number;

		if (number > 0)
		{	
			m_cmf[tot]  = i;
		}
	}
	assert(tot > 0);

    r = static_cast<NUM_ELEMENTS>(random() % tot);		// choose our value 	



	//Finds the first element whose key greater than r. 
	// the index is stored in the second element of the cmf (total, index) pair
	upper_bound_iter = m_cmf.upper_bound(r);

	if (upper_bound_iter == m_cmf.end())
	{
		assert(r == tot);
		i = m_cmf[r];
	}
	else
	{
		i = upper_bound_iter->second;
	}

	assert(pmf[i] > 0);
    return i;				// return that range index 
}

// return a choice using the cummulative mass function
//
//  Choose uniformly from a cummalative weighted array.  i.e., a pmf [3,3,4]
//  would give a cmf [3,6,10]
//
// PRE: a is the number of choices
//      total is the sum over all weights
//      cmf is the cummulative mass function
//
//	RETURNS: a choice from the cmf 
NUM_ELEMENTS RANDOM_NUMBER::discrete_cmf
(
 	const NUM_ELEMENTS& n, 
	const double& total, 
	CUMULATIVE_MASS_FUNCTION & cmf
) const
{
	assert(total > 0);
    double r = 0.0;
	NUM_ELEMENTS index = 0;	
	CUMULATIVE_MASS_FUNCTION::const_iterator upper_bound_iter;

	r = random_double_number(1)* total;		// choose our value

	upper_bound_iter = cmf.upper_bound(r);

	if (upper_bound_iter == cmf.end())
	{
		// boundary condition
		assert(r == total);
		index = cmf[r];
	}
	else
	{
		index = upper_bound_iter->second;
	}

    return index;
}


//  Randomly choose a row and col from the matrix with the value at each 
//  coordinate being the weight that it will be chosen
// 
//  POST: row_choose, col_chosen are the randomly chosen coordinates
void RANDOM_NUMBER::discrete_pmf
(
	MATRIX& matrix_pmf, 
	NUM_ELEMENTS & row_chosen, 
	NUM_ELEMENTS & col_chosen
)
{
	INDEX_SIZE nCol = matrix_pmf.get_nColumns();
	NUM_ELEMENTS array_size;
	NUM_ELEMENTS array_position;
	const NUM_ELEMENTS_VECTOR & pmf = matrix_pmf.get_data();

	array_size = pmf.size();
    
	array_position = discrete_pmf(array_size, pmf);

	row_chosen = array_position / nCol;
	col_chosen = array_position - row_chosen*nCol;

	assert(row_chosen*nCol + col_chosen == array_position);
}
//  Randomly choose a x,y,z coordinate from the matrix with the value at each 
//  coordinate being the weight that it will be chosen
// 
//  POST: x_chosen,y_chosen,z_chosen are the coordinates chosen
void RANDOM_NUMBER::discrete_pmf
(
	MATRIX_3D& matrix_pmf, 
	NUM_ELEMENTS & x_chosen, 
	NUM_ELEMENTS & y_chosen,
	NUM_ELEMENTS & z_chosen
) 
{
    INDEX_SIZE nY = matrix_pmf.get_nY();
    INDEX_SIZE nZ = matrix_pmf.get_nZ();
	NUM_ELEMENTS array_position = 0;
	NUM_ELEMENTS temp_array_position = 0;

	const NUM_ELEMENTS_VECTOR & pmf = matrix_pmf.get_data();

	NUM_ELEMENTS array_size = pmf.size();
    
	array_position = discrete_pmf(array_size, pmf);

	x_chosen = array_position / (nY * nZ);
	temp_array_position = array_position - x_chosen*(nY * nZ);
	y_chosen = temp_array_position / nZ;
	z_chosen = temp_array_position - y_chosen * nZ;

	assert(x_chosen*(nY * nZ) + y_chosen* nZ + z_chosen == array_position);
	assert(matrix_pmf(x_chosen, y_chosen, z_chosen) > 0);
}
//  Randomly choose a row and col from the matrix.
//  Positive entries in the matrix have a prob. of being selected based on the value
//  of the entry.
//  Negative entries in the matrix have no prob. of being selected
// 
//  POST: row_choose, col_chosen are the randomly chosen coordinates
void RANDOM_NUMBER::discrete_pmf_with_negative_entries
(
	const MATRIX& matrix_pmf, 
	NUM_ELEMENTS & row_chosen, 
	NUM_ELEMENTS & col_chosen
) 
{
	INDEX_SIZE nRows = matrix_pmf.get_nRows();
	INDEX_SIZE nCol = matrix_pmf.get_nColumns();
	INDEX_SIZE row, col;
	NUM_ELEMENTS array_size;
	NUM_ELEMENTS array_position;
	NUM_ELEMENTS matrix_value;
	m_pmf.clear();

    for (row = 0; row < nRows; row++)
	for (col = 0; col < nCol ; col++)
	{
		matrix_value = matrix_pmf.get_value(row, col);
		matrix_value = MAX(0, matrix_value);
		m_pmf.push_back(matrix_value);
	}

	array_size = m_pmf.size();
    
	array_position = discrete_pmf(array_size, m_pmf);

	row_chosen = array_position / nCol;
	col_chosen = array_position - row_chosen*nCol;

	assert(row_chosen*nCol + col_chosen == array_position);
}

//  Randomly choose a row and col from the matrix.
//  Negative entries in the matrix have the positive prob. of being selected 
//  based on the absolute value of the entry
//  Positive entries in the matrix have no prob. of being selected
// 
//  POST: row_choose, col_chosen are the randomly chosen coordinates
void RANDOM_NUMBER::discrete_inverse_pmf_with_negative_entries
(
	const MATRIX& matrix_pmf, 
	NUM_ELEMENTS & row_chosen, 
	NUM_ELEMENTS & col_chosen
) 
{
	INDEX_SIZE nRows = matrix_pmf.get_nRows();
	INDEX_SIZE nCol = matrix_pmf.get_nColumns();
	INDEX_SIZE row, col;
	NUM_ELEMENTS array_size;
	NUM_ELEMENTS array_position;
	NUM_ELEMENTS matrix_value;
	m_pmf.clear();

    for (row = 0; row < nRows; row++)
	for (col = 0; col < nCol ; col++)
	{
		matrix_value = -1 * matrix_pmf.get_value(row, col);
		matrix_value = MAX(0, matrix_value);
		m_pmf.push_back(matrix_value);
	}

	array_size = m_pmf.size();
    
	array_position = discrete_pmf(array_size, m_pmf);

	row_chosen = array_position / nCol;
	col_chosen = array_position - row_chosen*nCol;

	assert(row_chosen*nCol + col_chosen == array_position);
}




//
// RETURNS: a random number between 0 and max_number
//
NUM_ELEMENTS RANDOM_NUMBER::random_number(const NUM_ELEMENTS& max_number) const
{
	NUM_ELEMENTS number = static_cast<NUM_ELEMENTS>(random() % (max_number+1));

	return number;
}
//
// RETURNS: a random double number between 0 and max_number
//
double RANDOM_NUMBER::random_double_number(const NUM_ELEMENTS& max_number) const
{
	assert(new_gen);
	// is this enough precision?
	double number = static_cast<double>(random() % 10000*max_number+1)/10000;

	return number;
}

//  Choose uniformly from a weighted array.  i.e., with array [3,3,4]
//  return 0 with prob .3, 1 with prob .3 and 2 with prob .4.
//
//  If doing multiple queries, is better to use the cmf version, which 
//  takes only log(n) time; this one is linear.
//
NUM_ELEMENTS RANDOM_NUMBER::discrete_pmf(const DOUBLE_VECTOR & pmf)
{
    NUM_ELEMENTS n = static_cast<NUM_ELEMENTS>(pmf.size());

	return discrete_pmf(n, pmf);
}

//  Choose uniformly from a weighted array.  i.e., with array [3,3,4]
//  return 0 with prob .3, 1 with prob .3 and 2 with prob .4.
//
//  If doing multiple queries, is better to use the cmf version, which 
//  takes only log(n) time; this one is linear.
//
//	RETURNS: a choice from the prob. mass function
NUM_ELEMENTS RANDOM_NUMBER::discrete_pmf(const NUM_ELEMENTS & n, const DOUBLE_VECTOR & pmf)
{
    NUM_ELEMENTS i;
    double r, tot;
	double number;
	m_double_cmf.clear();
	map<double, NUM_ELEMENTS>::iterator upper_bound_iter;

	tot = 0.0;
    for (i = 0; i < n; i++)
	{
		number = pmf[i];
		assert(number >= 0);

		// to avoid round off errors making the map non unique
		if (tot + number > tot)
		{	
			tot += number;
			m_double_cmf[tot]  = i;
		}
	}
	assert(tot > 0.0);

    r = random_double_number(1)* tot;		// choose our value

	upper_bound_iter = m_double_cmf.upper_bound(r);

	if (upper_bound_iter == m_double_cmf.end())
	{
		// boundary condition
		assert(r == tot);
		i = m_double_cmf[r];
	}
	else
	{
		i = upper_bound_iter->second;
	}

	assert(pmf[i] > 0);

    return i;				// return that range index 
}


// RETURN: the probability of the sample value being selected in the gaussian distribution
//         
// note: M_PI is pi defined in math.h
//
double RANDOM_NUMBER::get_gaussian_probability
(
	const double& sample, 
	const double& mean, 
	const double& standard_dev
) const
{
	double probability = 0.0;
	double exp_value = 0.0;

	// don't deal with this special case. outlaw it
	assert(standard_dev > 0.0);

	exp_value = (sample - mean) * (sample - mean)/ (2 * standard_dev * standard_dev);

	probability = 1/(standard_dev * sqrt(2.0 * M_PI)) * exp(- exp_value);


	return probability;
}
