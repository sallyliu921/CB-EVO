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



#ifndef level_node_H
#define level_node_H

#include "types.h"
#include <vector>

class LEVEL_NODE;
typedef vector<LEVEL_NODE *> LEVEL_NODES;
#include "edges_and_nodes.h"
#include "spec.h"
#include "graph.h"

//
// Class_name LEVEL_NODE
//
// Description

//  A level structure is an intermediate step in the creation of a circuit.
//  This holds the nodes, and information about them, for one combinational
//  level of a cluster.
//



class LEVEL_NODE 
{
public:
	enum EDGE_TYPE {INTRA_CLUSTER_INPUT, INTRA_CLUSTER_OUTPUT, INTER_CLUSTER_INPUT, INTER_CLUSTER_OUPUT};
	LEVEL_NODE(const DELAY_TYPE & delay_level, CLUSTER_SPEC * comb_spec, 
			const NUM_ELEMENTS& nClusters,
			CLUSTER * my_cluster);
	~LEVEL_NODE();
	LEVEL_NODE(const LEVEL_NODE  & another_node);
	LEVEL_NODE & operator=(const LEVEL_NODE  & another_node);

	string			get_name() const { return m_name; }
	NODE*			get_node(const NUM_ELEMENTS& node_index) const;
	EDGES			get_input_edges() const { return m_input_edges; }
	EDGES			get_input_to_dff_edges() const { return m_input_to_dff_edges; }
	EDGES			get_output_edges() { return m_output_edges; }
	NUM_ELEMENTS	get_nNodes() const { return m_nNodes; }
	NUM_ELEMENTS	get_nCombinational_nodes() const { return (m_nNodes - m_nDFF); }
	NUM_ELEMENTS	get_nDFF() const { return m_nDFF; }
	NUM_ELEMENTS	get_nLatched() const { return m_nlatched; }
	NUM_ELEMENTS	get_nPO() const { return m_nPO; }
	CLUSTER_NUMBER_TYPE get_cluster_number() const { return m_cluster_number; }
	DELAY_TYPE		get_delay_level() const { return m_delay_level; }
	NUM_ELEMENTS	get_available_in() const { return m_available_in; }
	NUM_ELEMENTS	get_available_out() const { return m_available_out; }
	NUM_ELEMENTS	get_nIntra_cluster_edges() const 
					{ return m_nIntra_cluster_inputs+ m_nIntra_cluster_outputs; }
	NUM_ELEMENTS	get_nIntra_cluster_input() const { return m_nIntra_cluster_inputs; }
	NUM_ELEMENTS	get_nIntra_cluster_output() const { return m_nIntra_cluster_outputs; }
	NUM_ELEMENTS	get_nInter_cluster_input() const { return m_nInter_cluster_inputs; }
	NUM_ELEMENTS	get_nInter_cluster_output() const { return m_nInter_cluster_outputs; }
	NUM_ELEMENTS&	OutDegrees(const NUM_ELEMENTS& degree);
	NUM_ELEMENTS&	EdgeLengths(const LENGTH_TYPE& edge_length);
	DISTRIBUTION	get_fanout_distribution() const { return m_out_degrees; }
	NUM_ELEMENTS 	get_max_fanout() const { return m_out_degrees.size() -1; }
	NODES			get_nodes() const { return m_nodes; }
	NODES			get_PI() const;
	NODES			get_PO() const;
	NODES			get_DFF() const;
	NODES			get_latched_nodes_not_taken(const NUM_ELEMENTS& number_to_get) const;
	NODES			get_latched_nodes() const;
	CLUSTER *	get_my_cluster() const { return m_cluster; }
	
	void set_nNodes(const NUM_ELEMENTS &number_of_nodes);
	void set_available_in(const NUM_ELEMENTS& available_in)  { m_available_in = available_in; }
	void set_available_out(const NUM_ELEMENTS& available_out)  { m_available_out = available_out; }
	void set_new_max_fanout(const NUM_ELEMENTS& degree);
	void add_input_edge(EDGE * edge); 
	void add_output_edge(EDGE * edge); 
	void input_edge_has_changed(const EDGE * edge, const NUM_ELEMENTS& weight_addition);
	void output_edge_has_changed(const EDGE * edge, const NUM_ELEMENTS& weight_addition);
	void add_dff_input_edge(EDGE * edge);
	void add_output_to_dff_edge(EDGE * edge);

	EDGE * find_edge_with_source(const LEVEL_NODE * source_node) const;
	EDGE * find_edge_with_source(const CLUSTER_NUMBER_TYPE& source_cluster, 
								 const DELAY_TYPE& source_delay) const;

	bool	is_input_degree_satisfied() const;
	bool	is_output_degree_satisfied() const;
	bool	is_delay_satisfied() const;
	bool	is_an_input_edge(EDGE * edge) const;

	NUM_ELEMENTS	get_assigned_in() const { return m_assigned_in; }
	NUM_ELEMENTS	get_assigned_out() const { return m_assigned_out; }
	NUM_ELEMENTS	get_input_edge_weight() const { return m_nIntra_cluster_inputs+m_nInter_cluster_inputs;}
	NUM_ELEMENTS	get_output_edge_weight() const { return m_nIntra_cluster_outputs+m_nInter_cluster_outputs;}
	NUM_ELEMENTS	get_nSimple_input_edges() const;
	NUM_ELEMENTS	get_nSimple_comb_output_edges() const;
	NUM_ELEMENTS	get_input_to_spare() const;
	NUM_ELEMENTS	get_output_to_spare() const;
	NUM_ELEMENTS	get_input_difference() const;
	NUM_ELEMENTS	get_output_difference() const;
	NUM_ELEMENTS	get_delay_difference() const;
	NUM_ELEMENTS	get_delay_to_spare() const;
	NUM_ELEMENTS 	get_excess_number_of_inputs() const;
	NUM_ELEMENTS 	get_nOpen_inputs() const;
	NUM_ELEMENTS 	get_nOutputs_still_needed() const;
	NUM_ELEMENTS 	get_nInputs_still_needed() const;
	NUM_ELEMENTS 	get_nInputs_still_wanted() const;
	NUM_ELEMENTS	get_max_degree_possible() const;
	NUM_ELEMENTS	get_max_degree_assigned() const;
	NODE *			get_indiv_node_with_open_input();
	NUM_ELEMENTS    get_nExcess_inputs_on_indiv_nodes() const;
	NUM_ELEMENTS	get_nDouble_inputs_on_indiv_nodes() const;
	NUM_ELEMENTS	get_nDFF_loops() const;

	DISTRIBUTION	get_intra_cluster_input_edge_lengths() const { return m_intra_cluster_input_edge_lengths; }
	DISTRIBUTION	get_intra_cluster_output_edge_lengths() const { return m_intra_cluster_output_edge_lengths;}
	DISTRIBUTION	get_inter_cluster_input_edge_lengths() const { return m_inter_cluster_input_edge_lengths;}
	DISTRIBUTION	get_inter_cluster_output_edge_lengths() const { return m_inter_cluster_output_edge_lengths;}

	COST_TYPE		get_cost() const;
	COST_TYPE		get_output_cost() const;
	COST_TYPE		get_input_cost() const;
	COST_TYPE 		get_too_many_inputs_cost() const;
	COST_TYPE 		get_too_few_inputs_cost() const;
	COST_TYPE 		get_too_few_outputs_cost() const;
	COST_TYPE		get_wirelength_cost() const;
	COST_TYPE		get_nIndiv_nodes_too_few_inputs_cost() const;
	DISTRIBUTION 	get_number_of_unassigned_long_and_inter_cluster_edges_each_indiv_node_has() const;

	bool			only_dff_super_edge_is_this_level_node() const;


	void report_edge_lengths() const;
	void sanity_check() const;
	void set_lost(int new_lost) { lost = new_lost; }
	void add_PO(const NUM_ELEMENTS& nPO_to_add);

	void assign_degrees_to_nodes();
	void designate_PO();
	void increase_too_many_inputs_factor();	
	void increase_too_few_outputs_factor();	
	void increase_too_few_inputs_factor();	
	void decrease_too_many_inputs_factor();	
	void decrease_too_few_outputs_factor();	
	void decrease_too_few_inputs_factor();	

	bool reduce_weight_on_input_edges_to_eliminate_double_edge_connections();
	void reset_congestion_costs();
	void print_output_connections() const;
	void print_indiv_edge_input_connections() const;

	void remove_node(NODE * node_to_remove);
private:
	NUM_ELEMENTS 		m_nNodes;
	NUM_ELEMENTS 		m_nDFF;
	NODES				m_nodes;
	CLUSTER_NUMBER_TYPE m_cluster_number;
	DELAY_TYPE			m_delay_level;

	// level node to level node edges
	// the m_input_edges are ordered so that we may index into the array
	EDGES				m_input_edges;		
	EDGES				m_output_edges;

	EDGES				m_input_to_dff_edges;
	EDGES				m_output_to_dff_edges;

	NUM_ELEMENTS	m_nIntra_cluster_inputs, m_nIntra_cluster_outputs;
	NUM_ELEMENTS	m_nInter_cluster_inputs, m_nInter_cluster_outputs;

	CLUSTER *		m_cluster;

	DISTRIBUTION 	m_intra_cluster_input_edge_lengths; // Lengths of the assigned_out in-degrees.
	DISTRIBUTION 	m_intra_cluster_output_edge_lengths; 		// Lengths of the assigned_out out-degrees.
	DISTRIBUTION 	m_inter_cluster_input_edge_lengths; // Lengths of the assigned_out out-degrees.
	DISTRIBUTION 	m_inter_cluster_output_edge_lengths;// Lengths of the assigned_out out-degrees.
	DISTRIBUTION	m_out_degrees;						// Individual out-degrees assigned to level

	string			m_name;

	
	NUM_ELEMENTS	m_nPO, m_nlatched;
	
	NUM_ELEMENTS	m_min_in,  m_max_in;
	NUM_ELEMENTS	m_min_out, m_max_out;
	NUM_ELEMENTS	m_assigned_in, m_assigned_out;			// what the spec has assigned
	NUM_ELEMENTS	m_target_fanout;
	NUM_ELEMENTS	m_available_in, m_available_out;		// what we have left to assign
	NUM_ELEMENTS	m_kin;
	NUM_ELEMENTS	m_max_fanout;
	double			m_too_many_inputs_factor;
	double			m_too_few_outputs_factor;
	double			m_too_few_inputs_factor;
	double			m_mult_too_many_inputs_factor;
	double			m_mult_too_few_outputs_factor;
	double			m_mult_too_few_inputs_factor;

	NUM_ELEMENTS	get_total_edge_weight(const EDGES & edges) const;
	NUM_ELEMENTS 	get_edges_ok_for_reduction(DISTRIBUTION & edges_ok_for_reduction_pmf);
	NUM_ELEMENTS 	get_edges_for_reduction(DISTRIBUTION & edges_for_reduction_pmf) const;
	
  int lost;
};
#endif
