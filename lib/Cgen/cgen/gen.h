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



#ifndef GEN_H
#define GEN_H

#ifdef VISUAL_C
// to deal with stupid vc++ stl warnings.
#pragma warning(disable:4786)
#endif 


#include <fstream>
#include <assert.h>

const bool DEBUG    = true;
const bool DSPEC    = false;
const bool DPRE_DEGREE = false;
const bool DSEQ_LEVEL = false;
const bool new_gen	= true;				// all the bits of the code that need a new look
const bool DSANITY_CHECK = false;
const bool DDegree = false;
const bool DNODE_SPLIT = false;
const bool DSIMPLE_EDGE = false;
const bool DPOST_PROCESS = false;
const bool DMEDIC		= false;


#include <vector>
#include <string>
using namespace std;

#include "types.h"
#include "output.h"
#include "rand.h"
#include "gen_options.h"
#include "graph.h"

#define ABS(A)   ((A) < 0 ? -(A) : (A))
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#define MIN(A,B) ((A) > (B) ? (B) : (A))


extern GEN_OPTIONS * g_options;
extern RANDOM_NUMBER * g_rand_number_gen;

#endif
