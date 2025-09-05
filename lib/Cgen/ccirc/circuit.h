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



#ifndef circuit_H
#define circuit_H

class CIRCUIT;

#include "circ.h"
#include "sequential_level.h"
#include "cluster.h"
#include "degree_info.h"
#include "matrix.h"


//
// Class_name CIRCUIT
//
// Description
//
//		The main data class for circ.
//		It stores the basic data types such as nodes, edges, PI/PO
//		as well as the that that information in hierarchical form such as 
//		sequential levels, and clusters plus 
//		it stores some stat based classes
//
class CIRCUIT
{
public:
	CIRCUIT();
	CIRCUIT(const CIRCUIT & another_circuit);
	CIRCUIT & operator=(const CIRCUIT & another_circuit);
	~CIRCUIT();

	// circuit creation methods
	void	set_name(const string & new_name) { m_name = new_name;}
	PORT *	add_external_port(PORT * port, const PORT::EXTERNAL_TYPE & external_type);
	PORT *	create_and_add_external_port(const string & port_name, const PORT::EXTERNAL_TYPE & external_type);
	void 	remove_external_port(PORT * port_to_remove);
	void 	set_external_ports(PORTS & new_ports, PORT::IO_DIRECTION io_direction);

	NODE *	create_dff(const string & new_dff_name, const string & input_port_name);
	NODE * 	create_node(const string & node_name, const NODE::NODE_TYPE & node_type);
	EDGE * 	create_edge(PORT * source_port,PORT * sink_port,const LENGTH_TYPE & length);
	void 	remove_and_delete_edge(PORT * source_port, PORT * sink_port,
									EDGE * edge_to_delete);
	void	remove_from_node_list(NODE * node);
	void	set_global_clock(PORT * new_global_clock);
	void	erase_global_clock() { m_global_clock = 0; }


	// cluster methods
	void 	construct_clusters(const NUM_ELEMENTS & number_of_partitions);
	void 	generate_sub_clusters();
	void 	set_scaled_cost(const double& scaled_cost) { m_scaled_cost = scaled_cost; }

	// calculating the fanin/fanout degree
	void 	calculate_degree_information();
	void 	create_inter_cluster_matrix();

	// access methods for data in circuit
	PORTS	get_PI() {return m_PI;}
	PORTS	get_PI_with_clock();
	PORTS	get_PO() {return m_PO;}
	PORT *	get_global_clock() {return m_global_clock;}
	NODES&	get_nodes() { return m_nodes;}
	EDGES&	get_edges() { return m_edges;}
	NODES	get_dffs();
	
	NODES 	 			get_top_nodes();
	SEQUENTIAL_LEVEL*	get_sequential_level() { return m_sequential_level;}
	CLUSTERS  &			get_clusters() { return m_clusters;}
	DEGREE_INFO * 		get_degree_info() { return m_degree_info;}
	DELAY_TYPE 			get_maximum_combinational_delay();
	MATRIX &			get_inter_cluster_adjacency_matrix() { return m_inter_cluster_connections; }
	MATRIX &			get_inter_cluster_adjacency_matrix_for_dff() 
						{ return m_inter_cluster_connections_for_dff; }
	MATRIX 				get_inter_cluster_matrix(const LENGTH_TYPE& edge_length);

	string 			get_name() const 	{ return m_name;}
	NUM_ELEMENTS	get_nPI() const 	{ return m_PI.size(); }
	NUM_ELEMENTS	get_nPO() const 	{ return m_PO.size(); }
	NUM_ELEMENTS	get_nEdges() const;
	NUM_ELEMENTS	get_nEdges_without_clock_edges() const;
	NUM_ELEMENTS	get_nComb() const 	{ return m_number_comb_nodes;}
	NUM_ELEMENTS	get_nDFF() const 	{ return m_number_seq_nodes;}
	NUM_ELEMENTS	get_nNodes() const 	{ return m_number_comb_nodes+m_number_seq_nodes;}
	NUM_ELEMENTS	get_size() const   	{ return get_nNodes() + get_nPI(); }
	NUM_ELEMENTS	get_nClusters() const { return static_cast<NUM_ELEMENTS>(m_clusters.size()); }
	NUM_ELEMENTS	get_max_width() const;
	double			get_scaled_cost() const { return m_scaled_cost; }
	void get_cluster_stats(DISTRIBUTION& size_dist, DISTRIBUTION& nPI_dist,
							DISTRIBUTION& nDFF_dist, DISTRIBUTION& nIntra_cluster_edges_dist, 
							DISTRIBUTION& nInter_cluster_edges_dist, DISTRIBUTION& wirelength_dist);

	void	colour_graph(const NODE::COLOUR_TYPE & colour);
	void 	print_out_node_information() const;
	void 	draw_circuit_in_dot_format();
		
	void 	set_wirelength_approx(const COST_TYPE& wirelength_approx) { m_wirelength_approx = wirelength_approx;}
	COST_TYPE 	get_wirelength_approx() const { return m_wirelength_approx; }

	bool	is_sequential() const { return (m_number_seq_nodes > 0);}
	
	void	final_sanity_check();
private:
	string 				m_name;			// the name of the circuit	
	PORTS				m_PI;			// the collection of primary inputs
	PORTS				m_PO;			// the collection of primary ouputs
	PORT *				m_global_clock;	// in addition to being in the PO	
											// the clock is listed here
	NODES				m_nodes;
	EDGES				m_edges;

	NUM_ELEMENTS		m_max_fan_in;	// the max-fanin or lut-size of the nodes
	SEQUENTIAL_LEVEL* 	m_sequential_level;
	CLUSTERS			m_clusters;
	CLUSTER_NUMBER_TYPE m_max_clustering_depth;
	MATRIX				m_inter_cluster_connections;
	MATRIX	 			m_inter_cluster_connections_for_dff;

	NUM_ELEMENTS		m_number_seq_nodes;
	NUM_ELEMENTS		m_number_comb_nodes;
	NUM_ELEMENTS		m_number_edges;

	DEGREE_INFO	*		m_degree_info;
	double				m_scaled_cost;
	COST_TYPE			m_wirelength_approx;

	void	colour_nodes(const NODE::COLOUR_TYPE & colour);
	void 	colour_up_from_node(NODE * node, const NODE::COLOUR_TYPE & colour, NUM_ELEMENTS & number_nodes);
	void	decrement_node_count(NODE * node);
	void	increment_node_count(NODE * node);

	void divide_nodes_into_clusters();
	void divide_primary_inputs_into_clusters();
	void print_cluster_stats();
};

#endif
