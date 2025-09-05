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



#ifndef final_improvement_H
#define final_improvement_H


//
// Class_name SIMPLE_EDGE_MOVE
//
// Description
//
// A final edge assignment move
//
// 
#include "gen.h"

class SIMPLE_EDGE_MOVE
{
public:
	enum SIMPLE_MOVE_TYPE {NONE, EDGE_ROTATION, EDGE_SWAP_BETWEEN_SOURCES,
							CHANGE_DFF_EDGE_ASSIGNMENT};
	SIMPLE_EDGE_MOVE();
	~SIMPLE_EDGE_MOVE();
	void set_circuit(CIRCUIT_GRAPH * circuit) { m_circuit = circuit; }

	void set_type(const SIMPLE_MOVE_TYPE& move_type) { m_type = move_type; }

	void set_edge_a(EDGE * edge_a) { m_edge_a = edge_a; }
	void set_destination_sink_node(NODE * node) { m_dst_sink_node = node; }
	void set_edge_b(EDGE * edge_b) { m_edge_b = edge_b; }

	SIMPLE_MOVE_TYPE get_type() const { return m_type; }
	EDGE * get_edge_a() const { return m_edge_a; }
	NODE * get_destination_sink_node() const { return m_dst_sink_node; }
	EDGE * get_edge_b() const { return m_edge_b; } 

	bool is_valid() const;
	void clear_move();

	void make_and_commit_move(NODES& bad_node_list, COST_TYPE& wirelength_cost, COST_TYPE& bad_node_cost);
	void make_move();
	void make_temporary_move();
	void unmake_temporary_move();
	void update_bad_node_list(NODES& bad_node_list);
	void update_costs(COST_TYPE& wirelength_cost, COST_TYPE& bad_node_cost);
	void add_remove_node_from_node_violation_list_if_necessary(NODE * sink_node, NODES& bad_node_list);
	void change_super_edges_invoved();

	COST_TYPE get_changed_bad_node_cost(const COST_TYPE & nEdges);
	COST_TYPE get_bad_node_cost_change_of_sink_nodes_involved();
	COST_TYPE get_bad_node_cost_of_sink_nodes_involved();

	COST_TYPE get_changed_wirelength_cost(const COST_TYPE& wirelength_cost, const COST_TYPE& wirelength_wanted,
										const COST_TYPE& max_index, const COST_TYPE& nEdges);
	COST_TYPE get_changed_wirelength_cost2(const COST_TYPE& wirelength_cost, const COST_TYPE& wirelength_wanted,
										const COST_TYPE& max_index, const COST_TYPE& nEdges);
	COST_TYPE get_wirelength_change_of_sink_nodes_involved(); 
	COST_TYPE get_wirelength_cost_of_sink_nodes_involved();
private:
	SIMPLE_MOVE_TYPE m_type;

	EDGE * 	m_edge_a;		// for both moves
	EDGE *  m_edge_b;		// for edge swaps
	NODE *  m_dst_sink_node;// for translations
	CIRCUIT_GRAPH * m_circuit;

	COST_TYPE m_dff_loop_penalty_multipler;

	void get_sink_nodes_involved(NODE *& sink_node_1, NODE *& sink_node_2);
	void get_source_nodes_involved(NODE *& source_node_1, NODE *& source_node_2);
	void get_sink_clusters(CLUSTER *& sink_1_cluster, CLUSTER *& sink_2_cluster);
	void get_wirelength_change_of_sink_nodes_involved(COST_TYPE& wirelength_change_sink_1,
		 											COST_TYPE& wirelength_changed_sink_2);
};


//
// Class_name SIMPLE_EDGE_ITERATOR
//
// Description
//
// A class that iteratively improve the final edge assignment
//
// 
class SIMPLE_EDGE_ITERATOR
{
public:
	SIMPLE_EDGE_ITERATOR();
	~SIMPLE_EDGE_ITERATOR();

	void iterate(CIRCUIT_GRAPH * circuit);
private:
	CIRCUIT_GRAPH *		m_circuit;
	CLUSTERS			m_clusters;
	NODES				m_nodes_with_violations;
	SIMPLE_EDGE_MOVE	m_move;
	bool				m_display_costs;
	EDGES				m_edges;
	COST_TYPE			m_max_index;
	COST_TYPE			m_nEdges;
	COST_TYPE			m_bad_node_cost;
	COST_TYPE			m_wirelength_cost;
	COST_TYPE			m_wirelength_wanted;

	COST_TYPE 			m_eta;
	COST_TYPE			m_dff_loop_penalty_multipler;
	//COST_TYPE			m_too_many_inputs
		

	COST_TYPE			m_important_bad_node_cost;


	COST_TYPE get_cost(); // the cost function
	COST_TYPE get_cost_using_wirelength_in_clusters(); // unused

	COST_TYPE get_changed_cost();
	COST_TYPE get_changed_cost2();


	// move generator functions
	void generate_move();
	void bad_node_generate_move();
	void generate_move_to_solve_dff_loop(NODE * sink_node);
	void generate_move_that_swaps_edge_between_dff(EDGE * edge_that_causes_loop);
	NODE* get_a_latched_node_in_same_cluster(const CLUSTER * cluster) const;


	void create_edges_to_rectify_delay_structure_with_individual_edge_assignment();
	void add_lost_edges_back_to_nodes();
	void try_to_fix_level_node(LEVEL_NODE * level_node);
	NODES get_lost_source_nodes(EDGE * edge);
	void make_list_of_node_that_have_violations();
	
	void update_temperature(const double& success_rate);
	void print_node_positions();
	void print_wirelength_stats() const;
	void update_cluster_to_reflect_which_simple_edges_they_contain();
};

#endif
