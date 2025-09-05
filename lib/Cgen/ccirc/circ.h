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


#ifndef circ_H
#define circ_H

#ifdef VISUAL_C
// to deal with stupid vc++ stl warnings.
#pragma warning(disable:4786)
#endif 

const bool DCODE 	= false;
const bool DEBUG 	= false;
const bool DBLIF 	= false;		// for debuggin the parser
const bool DCONST	= false;		// for debuggin the graph constructor
const bool DMEDIC	= false;		// for debuggin the graph medic
const bool DOPTIONS = true;			// for debuggin the options
const bool DSYMBOL_TABLE = false;
const bool DCYCLE 	= false;		// for debuggin the cycler breaker
const bool DLEVEL 	= false;		// for debuggin the leveler
const bool DCONE_GEN= false;		// for debuggin the cone generator
const bool DCLUSTER = false;		// for debuggin the clustering of nodes
const bool DSUB_CLUSTER = false;	// for debuggin the sub clustering of nodes
const bool DCOLOUR 	= false;		// for debuggin colouring of the graph
const bool DDELAY_LEVELER = true;
const bool DCONNECTION = false;		// for debugging inter cluster connection matrix
const bool Dlook_at	= true;			// when everything is done and gone look at these things
const bool Dsingle_seq_level = true;


/* utility macros */

#define ABS(A)   ((A) < 0 ? -(A) : (A))
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#define MIN(A,B) ((A) > (B) ? (B) : (A))

/* GLOBAL INCLUDE */

#include <assert.h>
#include "edges_and_nodes.h"
#include "options.h"
#include "types.h"
#include <vector>
#include <string>
#include <iostream>
using namespace std;

extern OPTIONS * g_options;

#include "output.h"

#endif
