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
#include "circ.h"
#include <map>
#include <math.h>

#ifdef VISUAL_C
static long g_current_seed = 0L;      // for MSVC substitute
const double M_PI = 3.14159265358979323846;
#endif

RANDOM_NUMBER::RANDOM_NUMBER()
{
	//  Initial state vector taken from random() man page.
	
	unsigned long state1[32] = { 3,
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
	
	unsigned long state1[32] = { 3,
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
	assert(Dlook_at);
	// is this enough precision?
	double number = static_cast<double>(random() % 10000*max_number+1)/10000;

	return number;
}
