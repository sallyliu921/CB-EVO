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

#include "types.h"

class SEQUENTIAL_LEVEL;

typedef vector<SEQUENTIAL_LEVEL *> SEQUENTIAL_LEVELS;

#include "spec.h"
#include "graph.h"
#include "level_node.h"
#include "move.h"
#include "matrix3D.h"

//
// Class_name SEQUENTIAL_LEVEL
//
// Description

//  A class to create the delay structure.
//
//  It is an uber class.
//
//  It is big and messy and needs to be broken down into smaller chunks.
//
//

typedef vector<LEVEL_NODE *> DELAY_LEVEL;
typedef vector<DELAY_LEVEL *> DELAY_LEVEL_VECTOR;

class SEQUENTIAL_LEVEL 
{
public:
	SEQUENTIAL_LEVEL(CIRCUIT_GRAPH* circuit);
	~SEQUENTIAL_LEVEL();

	void create_delay_structure();

	void add_level_nodes(const LEVEL_NODES & level_nodes, const CLUSTER_NUMBER_TYPE& cluster_number);

private:
	CIRCUIT_GRAPH *		m_circuit;
	DELAY_LEVEL_VECTOR 	m_delay_levels;
	DELAY_TYPE			m_max_comb_delay;
	EDGES				m_edges;
	EDGES				m_inter_cluster_edges;
	NUM_ELEMENTS		m_max_edge_weight;
	CLUSTERS			m_clusters;
	NUM_ELEMENTS		m_nClusters;
	MATRIX 				m_Comb;
	MATRIX *			m_Comb_spec;
	COST_TYPE			m_Comb_spec_sum;
	MATRIX				m_violation_matrix;
	DISTRIBUTION		m_edge_lengths;
	NUM_ELEMENTS		m_number_loops_without_reduction;
	ID_TYPE				m_node_id;

	MATRIX_3D			m_edges_between_cluster;
	MATRIX_3D			m_possible_edges_between_cluster;
	MATRIX_3D			m_source_edge_length_matrix;
	MATRIX_3D			m_dest_edge_length_matrix;
	bool				m_display_too_many_inputs_msg;
	bool				m_sanity_check_fatal;
	DISTRIBUTION 		m_edge_weight_state;

	COST_TYPE 			m_lamda;


	COST_TYPE			m_alpha,m_beta,m_gamma;

	enum DELETION_SITUATION {LEVEL_NODES_HAVE_SPARE_WEIGHT, LEVEL_NODES_HAVE_WEIGHT, NO_EDGE_CAN_BE_DELETED};

	SEQUENTIAL_LEVEL(const SEQUENTIAL_LEVEL  & another_sequential_level_generator);
	SEQUENTIAL_LEVEL & operator=(const SEQUENTIAL_LEVEL  & another_sequential_level_generator);

	// main functions
	void create_super_edges_between_level_nodes();
	void finish_creating_initial_solution_to_combinational_delay_structure();
	void find_solution();
	bool try_to_create_delay_structure();
	bool is_valid_solution();
	void simulated_anneal();
	bool try_to_fix_graph();
	void report_solution_quality();
	void make_latched_to_dff_connections();


	//
	// Create the initial solution to the combinational delay structure
	//
	void create_super_edges_between_levels(DELAY_LEVEL* source_delay_level, DELAY_LEVEL* sink_delay_level);
	EDGE* create_combinational_edge (LEVEL_NODE * source_level_node, LEVEL_NODE * sink_level_node);
	EDGE* create_dff_edge (LEVEL_NODE * source_level_node, LEVEL_NODE * sink_level_node);
	void create_combinational_edge_with_max_weight(LEVEL_NODE * source_level_node, LEVEL_NODE * sink_level_node);
	void create_dff_edges(LEVEL_NODE* source_level_node, DISTRIBUTION & destination_cluster_numbers);
	void create_violation_matrix();

	void remove_inter_cluster_output_edges(CLUSTER * cluster);
	void remove_intra_cluster_edges(CLUSTER * cluster);

	DELETION_SITUATION get_deletion_situation() const;
	void update_deletion_situation(SEQUENTIAL_LEVEL::DELETION_SITUATION & deletion_situation) const;

	bool reduce_inter_cluster_edge_weight(CLUSTER* cluster_a, CLUSTER* cluster_b,
									const DELETION_SITUATION & deletion_situation);
	void reduce_intra_cluster_edge_weight(CLUSTER * cluster);

	NUM_ELEMENTS get_number_edges_that_can_be_deleted(const LEVEL_NODE * source_level_node, 
													const LEVEL_NODE * sink_level_node, 
													const EDGE * edge, 
													const NUM_ELEMENTS& source_violations,
													const NUM_ELEMENTS& sink_violations,
													const bool& only_look_for_spare) const;


	DISTRIBUTION get_intra_cluster_edge_length_violations(CLUSTER * cluster);
	DISTRIBUTION get_inter_cluster_input_edge_length_violations(CLUSTER * cluster);
	DISTRIBUTION get_inter_cluster_output_edge_length_violations(CLUSTER * cluster);


	void add_weight_to_edges_in_cluster(CLUSTER * cluster, const LENGTH_TYPE& edge_length, 
											NUM_ELEMENTS & number_of_violations);

	LENGTH_TYPE choose_edge_length_with_violation(const NUM_ELEMENTS_VECTOR & edge_length_input_violations,
									const NUM_ELEMENTS_VECTOR & edge_length_output_violations) const;

	EDGE * get_edge_we_can_delete(CLUSTER * cluster, const LENGTH_TYPE& edge_length);
	EDGE * find_length_1_edge_to_translate(CLUSTER * cluster) const;

	void reduce_edges_to_meet_total_number_of_edges();
	void find_and_remove_edge_lengths(NUM_ELEMENTS & number_to_remove, const LENGTH_TYPE& edge_length);

	bool is_Comb_satisfied();

	//
	// Iterative algorithm to improve the combinational delay structure solution
	//

	void create_edge_length_matrix();
	void update_violation_matrix(MOVE & move);
	void update_edge_length_data_structures(MOVE & move);
	void update_cluster_for_edge_length_data_structures(const CLUSTER_NUMBER_TYPE & source_cluster_number,
														const CLUSTER_NUMBER_TYPE & sink_cluster_number,
														const LENGTH_TYPE& edge_length);
	
	void calculate_Comb();


	// Iterative algorithm improvement. cost functions
	COST_TYPE get_cost() const;
	COST_TYPE get_normalized_cost() const;
	COST_TYPE get_unnormalized_cost() const;
	COST_TYPE get_changed_cost(MOVE& move, const COST_TYPE& current_cost);
	COST_TYPE get_cost_of_clusters_involved(const MOVE& move) const;
	void increase_congestion_cost_for_level_nodes();


	// Iterative algorithm improvement. move generation
	void generate_move(MOVE& move);
	void generate_intra_cluster_move(MOVE& move, CLUSTER * cluster);
	void generate_inter_cluster_move(MOVE& move, CLUSTER * source_cluster, CLUSTER * sink_cluster);
	void generate_edge_length_move(MOVE& move);
	void generate_source_sink_inter_cluster_pair(
			CLUSTER_NUMBER_TYPE& source_cluster_number, 
			CLUSTER_NUMBER_TYPE& sink_cluster_number) const;
	void generate_source_sink_pair_for_edge_that_will_accept_weight(const MOVE& move,
			CLUSTER_NUMBER_TYPE& source_cluster_id,
			CLUSTER_NUMBER_TYPE& sink_cluster_id) const;
	EDGE * get_edge_that_will_give_up_weight(const CLUSTER * source_cluster, const CLUSTER * sink_cluster) const;
	EDGE * get_edge_that_will_accept_weight(const MOVE& move, const LENGTH_TYPE& edge_length) const;
	void make_swap_move(EDGE * inter_cluster_edge, CLUSTER * sink_cluster);


	bool are_the_inter_cluster_edge_length_distributions_satisfied();
	void get_number_edges_between_two_clusters(const CLUSTER_NUMBER_TYPE & source_cluster_number,
								const CLUSTER_NUMBER_TYPE &sink_cluster_number,
								DISTRIBUTION & number_edges_between_two_clusters,
								DISTRIBUTION & number_poss_edges_between_two_clusters);
	void reset_graph();
	void reset_congestion_costs();
	void add_edges_to_nodes_that_need_outputs();



	void create_edge_length_state();
	void assign_lowest_cost_edge_weights();
	void set_graph_into_lowest_cost_state();

	// Sanity checking
	void sanity_check_all() const;
	void sanity_check_edge_length_matricies();
	void sanity_check_dff_connections() const;
	void sanity_check(CLUSTER * cluster) const;
	
	//
	// Post processing of the graph to eliminate node violations
	//
	bool eliminate_excess_number_of_inputs();
	EDGE * choose_edge_to_add_output_weight_to_level_node(LEVEL_NODE * level_node,
														EDGES& lost_edges_for_cluster);
	EDGE* find_an_edge_to_increase_weight(const CLUSTER_NUMBER_TYPE& source_cluster_number,
										const CLUSTER_NUMBER_TYPE& sink_cluster_number,
										const LENGTH_TYPE& edge_length,
										const DELAY_TYPE& source_delay_level) const;
	EDGES get_lost_edges_for_cluster(const EDGES& lost_edges,const CLUSTER_NUMBER_TYPE& cluster_number) const;

	// functions that can be used for debuging
	void output_inter_cluster_matricies_by_edge_length() const;
	void output_all_violations();
	void output_difference_between_current_and_spec(CLUSTER * cluster);
	void output_edge_debug_message(LEVEL_NODE * source_level_node, LEVEL_NODE * sink_level_node);
	void output_spare(CLUSTER * cluster);
	// alternative move generators
	// void new_generate_move(MOVE& move);
	// void new_generate_moveII(MOVE& move);

};
#endif
