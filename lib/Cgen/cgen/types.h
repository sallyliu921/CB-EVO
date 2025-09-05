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


#ifndef TYPES_H
#define TYPES_H

/*
 *   Miscellaneous typedefs for circ/gen
 */


#include <vector>
#include <string>
#include <map>
using namespace std;


typedef long NUM_ELEMENTS;
typedef short DELAY_TYPE;
typedef short DEGREE_TYPE;
typedef short LEVEL_TYPE;
typedef short LUT_TYPE;
typedef short LENGTH_TYPE;
typedef short DEPTH_TYPE;
typedef long CLUSTER_NUMBER_TYPE;
typedef long ID_TYPE;
typedef long LOCALITY;
typedef double COST_TYPE;

typedef vector<NUM_ELEMENTS> 	NUM_ELEMENTS_VECTOR;
typedef NUM_ELEMENTS_VECTOR		DISTRIBUTION;
typedef NUM_ELEMENTS_VECTOR		SHAPE;
typedef vector<double>    DOUBLE_VECTOR;
typedef	map<double, NUM_ELEMENTS> CUMULATIVE_MASS_FUNCTION;

#endif
