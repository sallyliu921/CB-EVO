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



#ifndef cycle_breaker_H
#define cycle_breaker_H

#include "circuit.h"
#include "circ.h"

//
// Class_name CYCLE_BREAKER
//
// Description
// 		Breaks combinational cycles in the graph
//


class CYCLE_BREAKER
{
public:
	CYCLE_BREAKER();
	CYCLE_BREAKER(const CYCLE_BREAKER & another_cycle_breaker);
	CYCLE_BREAKER & operator=(const CYCLE_BREAKER & another_cycle_breaker);
	~CYCLE_BREAKER();
	void break_cycles(CIRCUIT * circuit);
private:
	CIRCUIT *		m_circuit;
	NUM_ELEMENTS	m_number_new_dff;
	NUM_ELEMENTS	m_number_moved_edges;

	bool			m_shown_breaking_cycles_warning;

	bool 	find_and_break_cycles(PORT * output_port, NODE_PTR_LIST & next_nodes, DEPTH_TYPE depth);
	bool	find_cycle_in_path_from_seq_nodes(NODE_PTR_LIST & next_nodes, DEPTH_TYPE depth);
	NODE *  find_latch_in_fanout(PORT * output_port);
	void 	break_cycle_by_adding_dff(PORT * source_port, EDGE * edge_to_break);

	void show_added_dff_warning();
	void show_breaking_cycles_warning();
	void show_number_dff_created_edges_moved_warning();
};


#endif
