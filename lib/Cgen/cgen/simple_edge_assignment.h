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



#ifndef simple_edge_assignment_H
#define simple_edge_assignment_H


//
// Class_name SIMPLE_EDGE_ASSIGNER
//
// Description
//
//		Assigns the simple edges to the simple nodes.
//
//		As this stage the simple nodes exist at the various delay levels.
//		The simple nodes have a fanout degree.
//		The simple nodes that are PI, latched, DFF have been determined.
//		The super edges between level nodes have weights assigned to them
//		indicating the number of simple edges that will connect the simple nodes 
//		that belong to them.
//
//		All that is left is the assign the simple edges including some degree of locality.
//
#include "gen.h"

class SIMPLE_EDGE_ASSIGNER
{
public:
	SIMPLE_EDGE_ASSIGNER(CIRCUIT_GRAPH * circuit);
	~SIMPLE_EDGE_ASSIGNER();

	void assign_simple_edges(CLUSTER *cluster);
private:
	NUM_ELEMENTS 	m_sum_out,
					m_sum_in,
					m_nSrc,
					m_nDst;
	NODES 			m_src,
					m_dst;
	NUM_ELEMENTS_VECTOR m_choice_list;
	CLUSTER * 		m_cluster;
	CIRCUIT_GRAPH * m_circuit;



	void initialize_variables(const NUM_ELEMENTS& nNodes);
	void create_destination_list(const DELAY_TYPE& delay);
	void add_length1_intra_cluster_edges_to_source_list(const DELAY_TYPE& delay);
	void add_long_and_inter_cluster_edges_to_source_list(const DELAY_TYPE& delay);
	void add_length1_inter_cluster_edges_to_source_list(const DELAY_TYPE& delay);
	void add_source_node(NODE * node, const NUM_ELEMENTS& edges_to_assign);
	void do_high_fanouts();
	void for_the_destination_list_find_a_start_position_and_range(const LOCALITY& high_fanout_srcindex, 
														    LOCALITY& range, int& start);
	void do_one_high_fanout(const LOCALITY& srcindex, LOCALITY& range);
	void remove_dst_nodes_with_full_inputs(const int&  start);
	void establish_delay_level(const DELAY_TYPE& delay);
	void ensure_no_buffer_nodes();
	void no_unwanted_outputs(const DELAY_TYPE& delay);
	void assign_remaining_edges();
	bool choose_and_make_a_single_connection(const int& iDst);
	bool choose_and_make_two_connections(const int& iDst);
	void do_connect(const int& iSrc, const int& iDst);
	void show_progress();
	string get_status_string();
	NODES select_simple_nodes_from_source_cluster(const EDGE * super_edge) const;
	
	void add_latched_nodes_to_source_list();
	void add_dff_nodes_to_dst_list();
	void assign_dff_simple_edges();
	void remove_unused_latched_nodes();
	void output_debug_msg() const;
	void output_warning_msgs_to_user(const LEVEL_NODE * level_node) const;
};

#endif
