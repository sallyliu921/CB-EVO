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



#ifndef random_number_H
#define random_number_H

#ifdef VISUAL_C
#pragma warning(disable:4786)	// stupid vc++ stl warnings
#endif 

#include "types.h"
#include "matrix.h"
#include <map>
//
// Class_name RANDOM_NUMBER
//
// Description
//		
/* 
 *  Routines for random objects:    M. Hutton.  Feb 95.
 *
 *  seed=randomize(seed)	-- randomize with seed, use clock if seed=0
 *                                 (automatic if not explicitly called)
 */

class RANDOM_NUMBER
{
public:
	RANDOM_NUMBER();
	RANDOM_NUMBER(long seed);
	~RANDOM_NUMBER();

	long randomize(long seed);

	long get_random_seed() const { return m_seed; }

	NUM_ELEMENTS 	random_number(const NUM_ELEMENTS& max_number) const;
	double 			random_double_number(const NUM_ELEMENTS& max_number) const;

private:
	long m_state1[32];		 //  Initial state vector taken from random() man page.
	long m_seed;
};


#ifdef VISUAL_C
long random(void);
#endif

#endif
