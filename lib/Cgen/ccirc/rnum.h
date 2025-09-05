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



#ifndef RNUM_H
#define RNUM_H

using namespace std;
#include "circ.h"
#include "circuit.h"
#include <math.h>
#include "edges_and_nodes.h"

/* EXPORTS */
void rnum(CIRCUIT * circuit,double *m_R0,double *m_R0max,double *m_R0min); // iterate through each PI and returned the weighted rnum value (quick count method )

/* INTERNALS */
static void _rnum_of_node(CIRCUIT * circuit, NODE * node,double *n0, int *d0);
static void _reset_marks(CIRCUIT * circuit);
static void _mark_outcone(CIRCUIT * circuit, NODE * node);
static void _quick_count(CIRCUIT * circuit, double *n0, int *d0);
static int _count_marked_fanin(NODE *node);


#endif