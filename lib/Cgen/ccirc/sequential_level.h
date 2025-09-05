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





#ifndef sequential_level_H
#define sequential_level_H

#include "circ.h"

class SEQUENTIAL_LEVEL;

typedef NODES 	DELAY_LEVEL;
typedef vector<DELAY_LEVEL> 	DELAY_LEVELS;
typedef vector<NUM_ELEMENTS>	SHAPE;
typedef vector<NUM_ELEMENTS>	DISTRIBUTION;

typedef vector<SEQUENTIAL_LEVEL *> SEQUENTIAL_LEVELS;

#include "matrix.h"

//
// Class_name SEQUENTIAL_LEVEL
//
// Description
//
//		Contains information about the sequential levels of the circuit.


class SEQUENTIAL_LEVEL
{
public:
	SEQUENTIAL_LEVEL();
	SEQUENTIAL_LEVEL(const LEVEL_TYPE & seq_level_number);
	~SEQUENTIAL_LEVEL();

	void set_max_combinational_delay(const DELAY_TYPE & max_delay);
	void add_node(NODE * node);
	void add_primary_input(PORT * primary_input);

	void set_is_clustered(const bool & is_clustered) { m_is_clustered = is_clustered;}

	DELAY_LEVELS &		get_delay_levels() { return m_delay_levels;}
	LEVEL_TYPE			get_sequential_level_number() const { return m_sequential_level_number;}
	NUM_ELEMENTS		get_nNodes() const { return m_number_of_nodes;}
	NUM_ELEMENTS		get_nPI() const { return m_PI.size(); }
	NUM_ELEMENTS		get_nPO() const { return m_PO.size(); }
	NUM_ELEMENTS		get_nInter_cluster_input_edges() const { return m_inter_cluster_input_edges.size(); }
	NUM_ELEMENTS		get_nInter_cluster_output_edges() const { return m_inter_cluster_output_edges.size(); }
	NUM_ELEMENTS		get_nLatched() const { return m_output_to_dff_edges.size() + 
													m_inter_cluster_output_to_dff_edges.size(); }
	DELAY_TYPE 			get_maximum_combinational_delay() const;
	NUM_ELEMENTS		get_max_width() const;

	SHAPE 				get_node_shape() const;
	SHAPE 				get_input_shape() const;
	SHAPE 				get_output_shape() const;
	SHAPE				get_PO_shape() const;
	SHAPE				get_latched_shape() const; 	// returns the position where nodes output to dff

	DISTRIBUTION		get_intra_cluster_edge_length_distribution() const;
	DISTRIBUTION		get_inter_cluster_input_edge_length_distribution() const;
	DISTRIBUTION		get_inter_cluster_output_edge_length_distribution() const;
	DISTRIBUTION		get_fanout_distribution(const NUM_ELEMENTS & max_fanout) const;


	EDGES&				get_internal_edges() { return m_internal_edges; }

	bool				is_clustered() const { return m_is_clustered; }

	void final_sanity_check();

	DISTRIBUTION		get_intra_cluster_input_edge_length_distribution(const DELAY_TYPE& delay_level);
	DISTRIBUTION		get_intra_cluster_output_edge_length_distribution(const DELAY_TYPE& delay_level);
	DISTRIBUTION		get_inter_cluster_input_edge_length_distribution(const DELAY_TYPE& delay_level);
	DISTRIBUTION		get_inter_cluster_output_edge_length_distribution(const DELAY_TYPE& delay_level);
private:
	SEQUENTIAL_LEVEL(const SEQUENTIAL_LEVEL & another_sequential_level);
	SEQUENTIAL_LEVEL & operator=(const SEQUENTIAL_LEVEL & another_sequential_level);

	DELAY_LEVELS		m_delay_levels;
	PORTS				m_PI;
	PORTS				m_dff_output_ports;
	PORTS				m_PO;
	PORT * 				m_clock_port;

	DELAY_TYPE			m_max_delay;
	LEVEL_TYPE			m_sequential_level_number;
	NUM_ELEMENTS		m_number_of_nodes;	
	bool				m_is_clustered;


	EDGES		m_internal_edges;	// edges that are internal to this sequential level and cluster
	EDGES		m_output_to_dff_edges;// edges that output to a dff
	EDGES		m_input_from_dff_edges;

	EDGES		m_inter_cluster_input_edges;
	EDGES		m_inter_cluster_output_edges;
	EDGES		m_inter_cluster_output_to_dff_edges;// edges that output to a dff 
	EDGES		m_inter_cluster_input_from_dff_edges;

	void add_node_to_delay_level(NODE * node);
	void add_edges_of_node(NODE * node);
	void add_input_edges(NODE * node);
	void add_input_edge(EDGE * edge, EDGES & edges_within_sequential_level,
						EDGES & inputs_to_dff_edges);
	void add_primary_input_edge(EDGE * edge, EDGES & edges_within_sequential_level);
	void add_output_edges(PORT * output_port);
	void add_output_to_dff_edge(EDGE * edge, EDGES & output_to_dff_edges);

	void add_latched_nodes_to_shape(SHAPE & latched_shape, const EDGES & output_to_dff_edges) const;

	void add_fanout_distribution_for_delay_level(DISTRIBUTION & fanout_distribution, 
												const DELAY_LEVEL & delay_level) const;
	void add_fanout_distribution_for_primary_inputs(DISTRIBUTION & fanout_distribution) const;

	DISTRIBUTION find_edge_length_distribution(const EDGES& edges) const;


	NUM_ELEMENTS get_number_of_non_clock_inputs(const NODES & nodes_at_a_delay_level) const;
	NUM_ELEMENTS get_number_of_outputs(const NODES & nodes_at_a_delay_level) const;
	NUM_ELEMENTS get_number_of_outputs(const PORTS & ports_at_a_delay_level) const;
};


#endif
