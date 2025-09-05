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



#ifndef post_process_H
#define post_process_H

#include "gen.h"

//
// Class_name POST_PROCESS
//
// Description
//
//
//	Try to clean up some of the mistakes.
//


class POST_PROCESS
{
public:
	POST_PROCESS();
	POST_PROCESS(const POST_PROCESS & another_post_process);
	POST_PROCESS & operator=(const POST_PROCESS & another_post_process);
	~POST_PROCESS();

	void post_process(CIRCUIT_GRAPH * circuit);
private:
	CIRCUIT_GRAPH * m_circuit;
	CLUSTERS 	m_clusters;

	NUM_ELEMENTS m_nExcess_input_nodes_fixed,
				 m_nDouble_input_edge_nodes_fixed,
				 m_nNo_output_nodes_fixed,
				 m_nNo_input_nodes_fixed;

	void fix_number_of_primary_outputs();
	void try_and_find_output_for_node(NODE * node);
	bool try_and_find_input_for_node(CLUSTER* cluster, const DELAY_TYPE& delay_level, NODE * node);
	bool find_an_open_input_and_create_an_edge(NODE * source_node, 
									const CLUSTER_NUMBER_TYPE& sink_cluster_number,
									const LENGTH_TYPE& edge_length, const DELAY_TYPE& source_delay_level);
	EDGE * find_a_node_with_an_open_output(NODE * source_node, const CLUSTER_NUMBER_TYPE& sink_cluster_number,
									const LENGTH_TYPE& edge_length, const DELAY_TYPE& source_delay_level);

	void fix_nodes_with_double_inputs();
	void fix_nodes_with_too_many_inputs();
	void fix_nodes_with_no_inputs_and_outputs();
	void fix_nodes_with_no_outputs(NODE_LIST& nodes_with_no_outputs);
	void fix_nodes_with_no_inputs(NODE_LIST& nodes_with_no_inputs);
	void fix_flip_flops_with_connections_to_themselves();

	void print_warnings_as_to_what_was_fixed() const;
	void set_max_fanout_if_necessary(NODE * source_node);
};


#endif
