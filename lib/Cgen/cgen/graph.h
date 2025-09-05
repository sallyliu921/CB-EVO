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


// 
//   Graph structure in cgen.  
//   Graph structures include in order of hierarchy CIRCUIT_GRAPH and CLUSTER.
//
//

#ifndef GRAPH_H
#define GRAPH_H

class CLUSTER;
class CIRCUIT_GRAPH;

#include "spec.h"
#include "level_node.h"
#include "types.h"


typedef vector<CLUSTER *> CLUSTERS;
typedef vector<CLUSTER *> CLUSTERS;

#include "sequential_level.h"

//
//  Cluster graph structure.
//
class CLUSTER
{
public:
	CLUSTER(CLUSTER_SPEC * new_comb_spec,
				const CLUSTER_NUMBER_TYPE & cluster_number, const NUM_ELEMENTS& nClusters, 
				CIRCUIT_GRAPH * circuit);
	~CLUSTER();

	COST_TYPE get_level_node_normalization_cost() const { return m_level_node_normalization_cost; }


	void add_edge(EDGE * edge) { m_edges.push_back(edge);}
	void add_intra_cluster_simple_edge(EDGE * edge) { m_intra_cluster_simple_edges.push_back(edge); }
	void add_inter_cluster_input_simple_edge(EDGE * edge) { m_inter_cluster_input_simple_edges.push_back(edge);}
	void add_inter_cluster_output_simple_edge(EDGE * edge){m_inter_cluster_output_simple_edges.push_back(edge);}
	void add_node(NODE * node) { m_nodes.push_back(node); }
	void inter_cluster_input_edge_has_changed(const EDGE * edge, const NUM_ELEMENTS& weight_addition);
	void inter_cluster_output_edge_has_changed(const EDGE * edge, const NUM_ELEMENTS& weight_addition);
	void intra_cluster_edge_has_changed(const EDGE * edge, const NUM_ELEMENTS& weight_addition);


	CLUSTER_NUMBER_TYPE get_cluster_number() const { return m_cluster_number; }
	NODES			get_nodes() const { return m_nodes; }
	NUM_ELEMENTS	get_nNodes() const { return m_nNodes; }
	NUM_ELEMENTS	get_nDFF() const { return m_nDFF; }
	NUM_ELEMENTS	get_nLatched() const { return m_nLatched; }
	NUM_ELEMENTS	get_nPI() const { return m_nPI; }
	NUM_ELEMENTS	get_nPO() const { return m_nPO; }
	NUM_ELEMENTS	get_nIntra_cluster_edges() const;
	NUM_ELEMENTS	get_nInter_cluster_input_edges() const;
	NUM_ELEMENTS	get_nInter_cluster_output_edges() const;
	NUM_ELEMENTS	get_nInter_cluster_input_edges(const LENGTH_TYPE& length) const;
	NUM_ELEMENTS	get_nIntra_cluster_edges(const LENGTH_TYPE& length) const;
	NUM_ELEMENTS	get_nInter_cluster_output_edges(const LENGTH_TYPE& length) const;
	NUM_ELEMENTS	get_kin() const { return m_kin; }
	DELAY_TYPE		get_delay() const { return m_delay; }
	CLUSTER_SPEC *		get_comb_spec() const { return m_comb_spec; }
	LENGTH_TYPE		get_max_edge_length() const;
	NUM_ELEMENTS	get_max_fanout() const { return m_max_fanout; }
	NUM_ELEMENTS	get_width() const;
	CIRCUIT_GRAPH*	get_my_circuit() const { return m_my_circuit; }

	
	LEVEL_NODES 		Level;		// Levels 0..delay

	DISTRIBUTION get_intra_cluster_edge_lengths() const;// { return m_intra_cluster_edge_lengths; }
	DISTRIBUTION get_inter_cluster_input_edge_lengths() const { return m_inter_cluster_input_edge_lengths; }
	DISTRIBUTION get_inter_cluster_output_edge_lengths() const { return m_inter_cluster_output_edge_lengths; }
	DISTRIBUTION get_input_shape_difference_from_spec() const;
	DISTRIBUTION get_output_shape_difference_from_spec() const;
	DISTRIBUTION get_delay_difference_by_delay_level() const;
	DISTRIBUTION get_intra_cluster_simple_edge_lengths() const;
	DISTRIBUTION get_inter_cluster_simple_input_edge_lengths() const;
	DISTRIBUTION get_inter_cluster_simple_output_edge_lengths() const;
	DISTRIBUTION get_degree_distribution();

	SHAPE		 get_latched_shape() const;

	COST_TYPE get_wirelength() const { return m_wirelength; }
	COST_TYPE get_cost() const;
	COST_TYPE get_wirelength_cost() const;
	COST_TYPE get_normalized_cost() const;
	COST_TYPE get_changed_wirelength_cost(const COST_TYPE& wirelength_change_of_sink_node) const;
	COST_TYPE get_level_node_cost() const;
	COST_TYPE get_shape_cost() const;
	COST_TYPE get_edge_length_cost() const;
	COST_TYPE get_normalized_edge_length_cost() const;
	COST_TYPE get_final_improvement_cost(COST_TYPE& total_bad_node_cost, const COST_TYPE& nEdges,
				const COST_TYPE& beta, const COST_TYPE& max_index);
	bool is_valid_solution() const;
	bool is_inter_cluster_edge_length_distributions_satisfied() const;
	bool is_intra_cluster_edge_length_distribution_satisfied() const;
	bool is_inter_cluster_input_edge_length_distribution() const;
	bool is_inter_cluster_output_edge_length_distribution() const;
	bool are_level_node_properties_satisfied() const;
	bool is_delay_satisfied() const;

	EDGE * get_edge(const DELAY_TYPE& source_delay_level, const LENGTH_TYPE & edge_length) const;
	
	EDGE * randomly_choose_an_intra_cluster_edge() const;
	void set_new_max_fanout(const NUM_ELEMENTS& degree);
	void add_PO(const NUM_ELEMENTS nPO) { m_nPO += nPO; }

	NUM_ELEMENTS get_degree_difference();
	NUM_ELEMENTS get_input_output_difference() const;
	NUM_ELEMENTS get_edge_length_difference() const;
	NUM_ELEMENTS get_po_difference() const;

	void output_nodes_to_blif(ofstream& blif_file) const;
	void output_comb_nodes_to_blif(const LEVEL_NODE * level_node, ofstream& blif_file) const;
	void output_comb_nodes_to_blif(ofstream& blif_file) const;
	void output_dff_nodes_to_blif(ofstream& blif_file) const;

	void output_nodes_to_vhdl(ofstream& vhd_file) const;
	void output_comb_nodes_to_vhdl(ofstream& vhdl_file) const;
	void output_dff_nodes_to_vhdl(ofstream& vhdl_file) const;
	void output_comb_nodes_to_vhdl(const LEVEL_NODE * level_node, ofstream& blif_file) const;


	void add_level_nodes_to_sequential_level(SEQUENTIAL_LEVEL * sequential_level);

	void remove_edge(EDGE * edge);
	void clear_simple_edges();

	void sanity_check() const;
private:
	string				m_name;			// name of the graph/circuit
	CLUSTER_NUMBER_TYPE m_cluster_number;
	DELAY_TYPE 			m_delay;		// combinational delay   
	NUM_ELEMENTS		m_max_fanout;
	NUM_ELEMENTS		m_kin; 			// maximum fanin to a node
	CLUSTER_SPEC *		m_comb_spec;	// Pointer to source spec for the graph 

	EDGES			m_edges;			// level node to level node edges
	EDGES			m_intra_cluster_simple_edges;
	EDGES			m_inter_cluster_input_simple_edges;
	EDGES			m_inter_cluster_output_simple_edges;
	NODES 			m_nodes;
	NUM_ELEMENTS	m_nLatched;			// Number of nodes attached to a dff
	
	NUM_ELEMENTS	m_nNodes;			// Number of nodes in the graph 
	NUM_ELEMENTS	m_nPI, m_nPO;       // Number of PI/PO nodes 
	NUM_ELEMENTS	m_nDFF;
	NUM_ELEMENTS	m_nIntra_cluster_edges;
	NUM_ELEMENTS	m_nInter_cluster_input_edges;
	NUM_ELEMENTS	m_nInter_cluster_output_edges;

	DISTRIBUTION	m_intra_cluster_edge_lengths;
	DISTRIBUTION	m_inter_cluster_input_edge_lengths;
	DISTRIBUTION	m_inter_cluster_output_edge_lengths;

	COST_TYPE		m_level_node_normalization_cost;
	COST_TYPE		m_wirelength;
	COST_TYPE		m_wirelength_wanted;

	CIRCUIT_GRAPH*	m_my_circuit;

	// recalculation from level nodes of the edge lengths
	void 		get_edge_lengths(DISTRIBUTION& intra_cluster, DISTRIBUTION& inter_cluster_input, 
								DISTRIBUTION& inter_cluster_output) const;
	
};


//
// Class_name CIRCUIT_GRAPH
//
// Description
// 		Global data structure for the generated circuit

class CIRCUIT_GRAPH
{
public:
	CIRCUIT_GRAPH(CIRCUIT_SPEC * circuit_spec);
	CIRCUIT_GRAPH(const CIRCUIT_GRAPH & another_circuit);
	CIRCUIT_GRAPH & operator=(const CIRCUIT_GRAPH & another_circuit);
	~CIRCUIT_GRAPH();

	void add_cluster(CLUSTER * cluster) { m_clusters.push_back(cluster); }

	NUM_ELEMENTS		get_nClusters() const { return m_clusters.size(); }
	DELAY_TYPE			get_delay() const { return m_max_combinational_delay; }
	NUM_ELEMENTS		get_nNodes() const { return m_number_seq_nodes + m_number_comb_nodes; }
	NUM_ELEMENTS		get_nComb() const { return m_number_comb_nodes;}
	NUM_ELEMENTS		get_nDFF() const { return m_number_seq_nodes;}
	NUM_ELEMENTS		get_nPI() const { return m_nPI; }
	NUM_ELEMENTS		get_nPO() const { return m_nPO; }
	NUM_ELEMENTS		get_nEdges() const { return m_number_edges; }
	NUM_ELEMENTS		get_width() const;
	NUM_ELEMENTS		get_max_index() const;
	EDGES&				get_simple_edges() { return m_simple_edges; }
	NUM_ELEMENTS		get_nDFF_loops() const;


	CLUSTER * 		get_cluster(const NUM_ELEMENTS & index) { return m_clusters[index]; }
	CLUSTERS		get_clusters() { return m_clusters; }
	NODES&			get_nodes() { return m_nodes; }
	NODES			get_PI() const;
	NODES			get_PO() const;

	ID_TYPE 		get_new_id();
	string			get_name() const { return m_name; }
	//NODE * 		create_node();
	EDGE * 			create_edge(LEVEL_NODE * source_level_node, LEVEL_NODE * sink_level_node);
	EDGE * 			create_edge(NODE * source_node, NODE * sink_node);
	EDGE * 			create_dff_edge(LEVEL_NODE * source_level_node, LEVEL_NODE * sink_level_node);
	void			add_node(NODE * node) { m_nodes.push_back(node); }

	CIRCUIT_SPEC*	get_circuit_spec() const { return m_circuit_spec; }
	SEQUENTIAL_LEVEL* get_sequential_level() { return m_sequential_level; }
	void add_lost_weight_edge(EDGE * edge) { m_lost_weight_edges.push_back(edge); }
	EDGES get_lost_weight_edges() { return m_lost_weight_edges; }
	void print_lost_edge_stats();
	COST_TYPE get_double_input_cost() const;
	void add_lost_source_node(NODE * node) { m_lost_source_nodes.push_back(node); }
	void remove_edge(EDGE * edge);

	bool is_combinational() const { return (m_number_seq_nodes == 0); }
	bool is_sequential() const { return (m_number_seq_nodes > 0); }

	bool is_the_inter_cluster_connection_matrix_satisfied();
	void create_inter_cluster_connection_matrix();


	void delete_node(NODE * node);
	void delete_unconnected_nodes();

	// the main steps in the algorithm
	
	// Step 1
	// creates the delay structure - the large scale connectivity in the circuit
	void create_delay_structure();

	// Step 2
	// partition the fanout amongst the level nodes
	void partition_degrees();

	// Step 3
	// split the level nodes into individual nodes and 
	// prepare the graph for final edge assignment
	void split_nodes();

	// Step 4
	// assign the simple edges between the individual nodes
	void assign_simple_edges();


	void output_final_circuit_quality();

	void output_blif();
	void output_vhdl();

	void final_sanity_check();
	void sanity_check();
private:

	string				m_name;
	NODES			 	m_nodes;
	//PORT_PTR_VECTOR		m_PI;
	//PORT_PTR_VECTOR		m_PO;

	EDGES				m_simple_edges;

	/*
	EDGE_PTR_VECTOR		m_intra_cluster_edges;
	EDGE_PTR_VECTOR		m_inter_cluster_edges;
	EDGE_PTR_VECTOR		m_clock_edges;
	*/
	EDGES				m_lost_weight_edges;
	NODES 				m_lost_source_nodes;


	NUM_ELEMENTS		m_number_seq_nodes;
	NUM_ELEMENTS		m_number_comb_nodes;
	NUM_ELEMENTS		m_number_edges;
	NUM_ELEMENTS		m_nPI;
	NUM_ELEMENTS		m_nPO;

	CLUSTERS			m_clusters;
	SEQUENTIAL_LEVEL *	m_sequential_level;
	//MATRIX				m_inter_cluster_connections;
	//MATRIX				m_inter_cluster_connection_matrix_for_dff;

	ID_TYPE				m_id;			// node id
	ID_TYPE				m_edge_id;		// simple edge id
	DELAY_TYPE			m_max_combinational_delay;
	LEVEL_TYPE			m_max_sequential_level;
	NUM_ELEMENTS		m_kin;


	void add_level_nodes(SEQUENTIAL_LEVEL * sequential_level);

	CIRCUIT_SPEC*		m_circuit_spec;


	bool m_sanity_check_display;
	NUM_ELEMENTS get_inter_cluster_connectivity_cost() const;
	void get_inter_cluster_connectivity_costs(NUM_ELEMENTS& cost_Comb, NUM_ELEMENTS& cost_Latched ) const;
	void create_sequential_levels();

	void output_primary_inputs_to_blif(ofstream& blif_file) const;
	void output_primary_outputs_to_blif(ofstream& blif_file) const;

	void output_top_most_entity(ofstream& vhdl_file) const;
};



#endif  /* GRAPH_H */
