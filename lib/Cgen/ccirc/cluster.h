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





#ifndef cluster_H
#define cluster_H

#include <vector>
using namespace std;

class CLUSTER;
typedef vector<CLUSTER *> CLUSTERS;	

#include "circ.h"
#include "sequential_level.h"
#include "degree_info.h"

//
// Class_name CLUSTER
//
// Description
//		Clusters are collection of nodes and primary inputs.
//		They are created after partitioning.
//		They are in essence mini-circuits.
//		I decided to create another class because I don't
//		want the circuit class to get to big.



class CLUSTER
{
public:
	enum CLUSTER_TYPE {REGULAR, SUB};
	CLUSTER();
	CLUSTER( const CLUSTER_NUMBER_TYPE & cluster_number,
			const DELAY_TYPE & max_comb_delay, const CLUSTER_TYPE & cluster_type);
	CLUSTER( const CLUSTER_NUMBER_TYPE & cluster_number, 
			const CLUSTER_NUMBERS & upper_sub_cluster_numbers,
			const CLUSTER_NUMBER_TYPE & sub_cluster_number,
			const DELAY_TYPE & max_comb_delay, const CLUSTER_TYPE & cluster_type);
	CLUSTER(const CLUSTER_NUMBER_TYPE & cluster_number, const DELAY_TYPE & max_comb_delay);
	~CLUSTER();

	void add_node(NODE * node);
	void add_primary_input(PORT * primary_input);
	void add_intra_cluster_edge(EDGE * edge);
	void add_inter_cluster_input_edge(EDGE * edge);
	void add_inter_cluster_output_edge(EDGE * edge);
	void set_depth(const CLUSTER_NUMBER_TYPE & depth) { m_depth = depth; }

	void calculate_degree_info();

	string				get_info() const;
	NUM_ELEMENTS		get_nNodes() const { return m_number_seq_nodes + m_number_comb_nodes; }
	NUM_ELEMENTS		get_nComb() const { return m_number_comb_nodes;}
	NUM_ELEMENTS		get_nDFF() const { return m_number_seq_nodes;}
	NUM_ELEMENTS		get_nPI() const { return m_PI.size(); }
	NUM_ELEMENTS		get_nPO() const { return m_PO.size(); }
	NUM_ELEMENTS		get_nIntra_cluster_edges() const { return m_intra_cluster_edges.size(); }
	NUM_ELEMENTS		get_nInter_cluster_edges() const { return m_inter_cluster_input_edges.size() +
																  m_inter_cluster_output_edges.size(); }
	NUM_ELEMENTS		get_nInter_cluster_input_edges() const { return m_inter_cluster_input_edges.size(); }
	NUM_ELEMENTS		get_nInter_cluster_output_edges() const { return m_inter_cluster_output_edges.size(); }
	NUM_ELEMENTS		get_size() const { return get_nPI() + get_nNodes(); }
	NUM_ELEMENTS		get_nLatched() const { return m_sequential_level->get_nLatched(); }

	EDGES&				get_intra_cluster_edges() { return m_intra_cluster_edges; }
	NODES &				get_nodes() 	{ return m_nodes; }
	PORTS 				get_PI() const	{ return m_PI; }
	CLUSTER_NUMBER_TYPE	get_cluster_number() const { return m_cluster_number; }
	CLUSTER_NUMBER_TYPE get_depth() const { return m_depth; }

	SHAPE 				get_node_shape() const;
	DEGREE_INFO *		get_degree_info() { return m_degree_info; }
	SEQUENTIAL_LEVEL *  get_sequential_level() { return m_sequential_level; }
	DISTRIBUTION		get_size_distribution(const CLUSTER_NUMBER_TYPE& depth);
	DISTRIBUTION		get_cluster_distance_distribution() const;
	DISTRIBUTION		get_PI_fanout_values() const;
	DISTRIBUTION		get_DFF_fanout_values() const;
	CLUSTER_NUMBER_TYPE get_max_clustering_depth() const { return m_max_clustering_depth; }
	double				get_wirelength_approx() const;
	void add_size_distribution(const CLUSTER_NUMBER_TYPE& depth, DISTRIBUTION & distribution);


	
	// splitting the cluster. maybe a good thing to get rid of all this?
	void	split_cluster(NUM_ELEMENTS number_times_to_split);
	void 	construct_sub_clusters(const NUM_ELEMENTS & number_of_partitions);
	CLUSTERS			get_sub_clusters() { return m_sub_clusters; }
	void add_sub_cluster_node(NODE * node);
	CLUSTER_NUMBERS		get_sub_cluster_numbers() const { return m_sub_cluster_numbers; }

	void final_sanity_check();
private:
	CLUSTER_TYPE			m_type;
	CLUSTER_NUMBER_TYPE 	m_cluster_number;
	CLUSTER_NUMBER_TYPE		m_depth;
	NODES					m_nodes;
	PORTS					m_PI;
	PORTS					m_PO;
	DELAY_TYPE				m_max_comb_delay;

	EDGES				m_intra_cluster_edges;
	EDGES				m_inter_cluster_input_edges;
	EDGES				m_inter_cluster_output_edges;

	NUM_ELEMENTS		m_number_seq_nodes;
	NUM_ELEMENTS		m_number_comb_nodes;

	SEQUENTIAL_LEVEL*	m_sequential_level;
	DEGREE_INFO * 		m_degree_info;


	void increment_node_count(const NODE * node);
	void add_edges(NODE * node);
	void add_input_edges(const NODE * node);
	void add_output_edges(const PORT * port);
	void add_non_clock_input(EDGE * edge);
	void print_stats();

	CLUSTER(const CLUSTER & another_cluster);
	CLUSTER & operator=(const CLUSTER & another_cluster);
	
	//
	// sub-cluster stuff. should remove it?
	//
	CLUSTER_NUMBERS 	m_sub_cluster_numbers;
	CLUSTERS			m_sub_clusters;
	CLUSTER_NUMBER_TYPE	m_max_clustering_depth; // for regular clusters
	void divide_nodes_into_clusters();
	void divide_intra_cluster_edges_into_clusters();
	void divide_inter_cluster_edges_into_clusters();
	void divide_primary_inputs_into_clusters();
};


#endif
