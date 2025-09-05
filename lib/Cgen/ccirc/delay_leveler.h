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



#ifndef delay_leveler_H
#define delay_leveler_H

#include "circ.h"
#include "circuit.h"

//
// Class_name DELAY_LEVELER
//
// Description
//
// 		Calculates the combinational delay levels for all nodes in the circuit
//
//


class DELAY_LEVELER
{
public:
	DELAY_LEVELER();
	DELAY_LEVELER(const DELAY_LEVELER & another_delay_leveler);
	DELAY_LEVELER & operator=(const DELAY_LEVELER & another_delay_leveler);
	~DELAY_LEVELER();

	void calculate_and_label_combinational_delay_levels(CIRCUIT * circuit);

private:
	CIRCUIT * 	m_circuit;
	DELAY_TYPE 	m_max_combinational_delay;

	void calculate_and_label_combinational_delay_for_sequential_level(SEQUENTIAL_LEVEL * sequential_level);
	void get_unmarked_fanout(NODE * node, NODE_PTR_DEQUE & nodes_to_find_delay_for);
												
	void label_all_dff_with_0_delay();
	void label_edges_with_length();
	void find_comb_delay_of_node(NODE * node, NODE_PTR_DEQUE & nodes_to_find_delay_for);

	void add_nodes_to_sequential_level(SEQUENTIAL_LEVEL* sequential_level);
	void add_primary_inputs_to_0th_sequential_level(SEQUENTIAL_LEVEL * seq_level_0);
	LEVEL_TYPE get_max_comb_delay_level_of_fanin(const NODE * node) const;
	void check_sanity() const;
};


#endif
