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





#include "cluster.h"
#include "node_partitioner.h"
#include "util.h"
#include <algorithm>

CLUSTER::CLUSTER()
{
	m_type				= CLUSTER::REGULAR;
	m_cluster_number	= -1;
	m_depth				= 0;
	m_number_seq_nodes 	= 0;
	m_number_comb_nodes	= 0;
	m_sequential_level	= 0;
	m_degree_info		= 0;
	m_max_comb_delay	= 0;
	m_max_clustering_depth = 0;

	m_degree_info = 0;
}
CLUSTER::CLUSTER
(
	const CLUSTER_NUMBER_TYPE & cluster_number,
	const DELAY_TYPE & max_comb_delay,
	const CLUSTER_TYPE & cluster_type
)
{
	assert(cluster_type == CLUSTER::REGULAR);

	m_type				= cluster_type;
	m_cluster_number	= cluster_number;
	m_depth				= 0;
	m_max_comb_delay 	= max_comb_delay;
	m_number_seq_nodes 	= 0;
	m_number_comb_nodes	= 0;
	m_degree_info 		= 0;
	m_max_clustering_depth = 0;

	m_sequential_level = new SEQUENTIAL_LEVEL(0);
	m_sequential_level->set_max_combinational_delay(max_comb_delay);
	m_sequential_level->set_is_clustered(true);
}
CLUSTER::CLUSTER
(
	const CLUSTER_NUMBER_TYPE & cluster_number,
	const CLUSTER_NUMBERS & upper_sub_cluster_numbers,
	const CLUSTER_NUMBER_TYPE & sub_cluster_number,
	const DELAY_TYPE & max_comb_delay,
	const CLUSTER_TYPE & cluster_type
)
{
	assert(cluster_type == CLUSTER::SUB);

	m_type				= cluster_type;
	m_cluster_number	= cluster_number;
	m_sub_cluster_numbers= upper_sub_cluster_numbers;
	m_sub_cluster_numbers.push_back(sub_cluster_number);
	m_depth				= static_cast<CLUSTER_NUMBER_TYPE>(m_sub_cluster_numbers.size());
	m_max_comb_delay 	= max_comb_delay;
	m_number_seq_nodes 	= 0;
	m_number_comb_nodes	= 0;
	m_degree_info 		= 0;
	m_sequential_level 	= 0;
	m_max_clustering_depth = 0;
}

CLUSTER::CLUSTER(const CLUSTER & another_cluster)
{
	m_type					= another_cluster.m_type;
	m_cluster_number		= another_cluster.m_cluster_number;
	m_sub_cluster_numbers	= another_cluster.m_sub_cluster_numbers;
	m_depth					= another_cluster.m_depth;
	m_nodes					= another_cluster.m_nodes;
	m_PI					= another_cluster.m_PI;
	m_PO					= another_cluster.m_PO;
	m_intra_cluster_edges	= another_cluster.m_intra_cluster_edges;
	m_inter_cluster_input_edges	= another_cluster.m_inter_cluster_input_edges;
	m_inter_cluster_output_edges= another_cluster.m_inter_cluster_output_edges;
	m_number_seq_nodes		= another_cluster.m_number_seq_nodes;
	m_number_comb_nodes 	= another_cluster.m_number_comb_nodes;
	m_sequential_level		= another_cluster.m_sequential_level;
	m_degree_info 			= another_cluster.m_degree_info;
	m_sub_clusters			= another_cluster.m_sub_clusters;
	m_max_clustering_depth 	= another_cluster.m_max_clustering_depth;
}

CLUSTER & CLUSTER::operator=(const CLUSTER & another_cluster)
{
	m_type					= another_cluster.m_type;
	m_cluster_number		= another_cluster.m_cluster_number;
	m_sub_cluster_numbers	= another_cluster.m_sub_cluster_numbers;
	m_depth					= another_cluster.m_depth;
	m_nodes					= another_cluster.m_nodes;
	m_PI					= another_cluster.m_PI;
	m_PO					= another_cluster.m_PO;
	m_intra_cluster_edges	= another_cluster.m_intra_cluster_edges;
	m_inter_cluster_input_edges	= another_cluster.m_inter_cluster_input_edges;
	m_inter_cluster_output_edges= another_cluster.m_inter_cluster_output_edges;
	m_number_seq_nodes		= another_cluster.m_number_seq_nodes;
	m_number_comb_nodes 	= another_cluster.m_number_comb_nodes;
	m_sequential_level		= another_cluster.m_sequential_level;
	m_degree_info 			= another_cluster.m_degree_info;
	m_sub_clusters			= another_cluster.m_sub_clusters;
	m_max_clustering_depth 	= another_cluster.m_max_clustering_depth;

	return (*this);
}

CLUSTER::~CLUSTER()
{
	if (m_degree_info)
	{
		delete m_degree_info;
		m_degree_info = 0;
	}
	
	if (m_sequential_level)
	{
		delete m_sequential_level;
		m_sequential_level = 0;
	}

	CLUSTERS::iterator cluster_iter;
	CLUSTER * cluster = 0;

	for (cluster_iter = m_sub_clusters.begin(); cluster_iter != m_sub_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);
		delete cluster;
	}
}

// 
// Adds a node to the cluster
//
// PRE: node is valid
// POST: the node, its edges, and ports are added to the cluster
//
void CLUSTER::add_node
(
	NODE * node
)
{
	assert(node);
	assert(m_type == CLUSTER::REGULAR);
	assert(node->get_cluster_number() == m_cluster_number);

	PORT * output_port = node->get_output_port();
	assert(output_port);

	debugif(DCLUSTER, "Adding node " << node->get_info() << " to cluster " << m_cluster_number);

	m_nodes.push_back(node);
	increment_node_count(node);

	m_sequential_level->add_node(node);
	add_edges(node);

	if (output_port->get_type() == PORT::EXTERNAL)
	{
		assert(output_port->get_external_type() == PORT::PO);
		debugif(DCLUSTER, "Adding PO port " << output_port->get_name() << " for node " << node->get_name() <<
				" to cluster " << m_cluster_number);
		m_PO.push_back(output_port);
	}

}


// 
// Adds the edges of this node to the cluster
//
// PRE: We have a regular cluster and the node is valid
// POST: the input and output edges of the nodes are added to the cluster
//
void CLUSTER::add_edges
(
	NODE * node
)
{
	assert(m_type == CLUSTER::REGULAR);
	const PORT * output_port = node->get_output_port();
	add_input_edges(node);
	add_output_edges(output_port);
}

// 
// Adds the input edges of this node to the cluster
//
// PRE: the node is valid
// POST: the input edges of the nodes are added to the cluster
//
void CLUSTER::add_input_edges
(
	const NODE * node
)
{
	assert(node);

	PORT * port = 0;
	EDGE * edge = 0;
	PORTS::iterator port_iter;

	PORTS input_ports = node->get_input_ports();

	for (port_iter = input_ports.begin(); port_iter != input_ports.end(); port_iter++)
	{
		port = *port_iter;
		assert(port);
		edge = port->get_edge();
		assert(edge);

		if (port->get_io_direction() != PORT::CLOCK)
		{
			add_non_clock_input(edge);
		}
	}
}


// 
// Adds the input edge that is not connected to a clock to the cluster
//
// PRE: the edge is valid
// POST: the input edge is added to the cluster
//
void CLUSTER::add_non_clock_input
(
	EDGE * edge
)
{
	assert(edge);
	NODE * sink_node = edge->get_sink_node();
	assert(sink_node);

	if (sink_node->get_type() != NODE::SEQ)
	{
		if (edge->is_intra_cluster())
		{
			add_intra_cluster_edge(edge);
		}
		else
		{
			assert(edge->is_inter_cluster());
			add_inter_cluster_input_edge(edge);
		}
	}
}

// 
// Adds the inter-cluster output edge to the cluster
//
// PRE: the output port is valid
// POST: if the edge is inter-cluster the edge connected to the output
//       port is added to the cluster
//
void CLUSTER::add_output_edges
(
	const PORT * output_port
)
{
	assert(output_port);

	EDGES::iterator edge_iter;
	EDGE * edge = 0;
	EDGES output_edges = output_port->get_edges();
	NODE * sink_node = 0;

	// look for ghost outputs
	for (edge_iter = output_edges.begin(); edge_iter != output_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		sink_node = edge->get_sink_node();
		assert(sink_node);

		if (sink_node->get_type() != NODE::SEQ)
		{
			if (edge->is_inter_cluster())
			{
				add_inter_cluster_output_edge(edge);
			}
			else
			{
				// intra cluster edges are alread taken care of with 
				// adding the input edges
			}
		}
	}
}

// 
// Adds the intra-cluster edge to the cluster
//
// PRE: the edge is valid and is intra-cluster
// POST: the edge is added to the cluster
//
void CLUSTER::add_intra_cluster_edge(EDGE * edge) 
{ 
	assert(edge && edge->is_intra_cluster());
	m_intra_cluster_edges.push_back(edge); 
	//debug("Adding intra cluster edge " << edge->get_name() << " to cluster " << m_cluster_number);
}

// 
// Adds the inter-cluster input edge to the cluster
//
// PRE: the edge is valid and inter-cluster
// POST: the edge is added to the cluster
//
void CLUSTER::add_inter_cluster_input_edge(EDGE * edge) 
{ 
	assert(edge);
	assert(edge->is_inter_cluster() || edge->is_inter_sub_cluster());
	assert(m_cluster_number == edge->get_sink_cluster_number());

	if (m_type == CLUSTER::SUB)
	{
		assert(edge->get_sink_sub_cluster_numbers() == m_sub_cluster_numbers);
	}

	//debug("Adding inter cluster input edge " << edge->get_name() << " to cluster " << m_cluster_number);

	m_inter_cluster_input_edges.push_back(edge); 
}
// 
// Adds the inter-cluster output edge to the cluster
//
// PRE: the edge is valid and inter-cluster
// POST: the edge is added to the cluster
//
void CLUSTER::add_inter_cluster_output_edge(EDGE * edge)
{ 
	assert(edge);
	assert(edge->is_inter_cluster() || edge->is_inter_sub_cluster());
	assert(m_cluster_number == edge->get_source_cluster_number());

	if (m_type == CLUSTER::SUB)
	{
		assert(edge->get_source_sub_cluster_numbers() == m_sub_cluster_numbers);
	}
	//debug("Adding inter cluster output edge " << edge->get_name() << " to cluster " << m_cluster_number);
	m_inter_cluster_output_edges.push_back(edge); 
}

// 
// Adds a primary input to the cluster
//
// PRE: the primary_input is valid and is a primary input
// POST: the port is added to the list of PIs
//       the output edges of the PI have been added to the cluster
//       
//
void CLUSTER::add_primary_input
(
	PORT * primary_input
)
{

	assert(primary_input && (primary_input->get_external_type() == PORT::PI || 
							 primary_input->get_external_type() == PORT::GI));
	assert(primary_input->get_sub_cluster_numbers() == m_sub_cluster_numbers);

	if (m_type == CLUSTER::REGULAR)
	{
		assert(m_sequential_level);

		//debugif(DCLUSTER, "Adding primary input " << primary_input->get_name() 
		//<< " to cluster " << m_cluster_number);
		m_sequential_level->add_primary_input(primary_input);
		add_output_edges(primary_input);
	}
	else
	{
		assert(m_type == CLUSTER::SUB);
	}

	m_PI.push_back(primary_input);
}

//
// Calculate the degree information of the cluster
//
// PRE: nothing
// POST: m_degree_info contains the fanin/fanout statistics
//
void CLUSTER::calculate_degree_info()
{
	m_degree_info = new DEGREE_INFO;
	assert(m_degree_info);

	m_degree_info->calculate_degree_information_for_cluster(this);
}

//
// Increments node counts
//
// PRE: node is valid
// POST: m_number_comb_nodes or m_number_seq_nodes has been incremented
//
void CLUSTER::increment_node_count
(
	const NODE * node
)
{
	assert(node);

	if (node->get_type() == NODE::COMB)
	{
		m_number_comb_nodes++;
	}
	else
	{
		assert(node->get_type() == NODE::SEQ);
		m_number_seq_nodes++;
	}
}

//
// Performs a sanity check on the cluster
// 
// PRE: nothing
// POST: The cluster number of all nodes is checked
// 		 If we have sub-clusters the number of nodes in the sub-clusters is checked
// 		 against the number of nodes in the cluster
// 		 The sub-cluster numbers are also checked
//
void CLUSTER::final_sanity_check()
{
	if (m_type == CLUSTER::REGULAR)
	{
		assert(m_sequential_level);
		m_sequential_level->final_sanity_check();
	}

	NODES::iterator node_iter;
	NODE * node = 0;
	CLUSTER_NUMBERS sub_cluster_numbers;
	NUM_ELEMENTS size = 0;
	NUM_ELEMENTS nPI = 0;
	CLUSTERS::iterator cluster_iter;
	CLUSTER * cluster = 0;
	
	for (node_iter = m_nodes.begin(); node_iter != m_nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		assert(node->get_cluster_number() == m_cluster_number);

		if (m_type == CLUSTER::SUB)
		{
			sub_cluster_numbers = node->get_sub_cluster_numbers();

			assert(m_depth >= 0 && static_cast<unsigned>(m_depth) < sub_cluster_numbers.size());
			assert(sub_cluster_numbers[m_depth] == m_sub_cluster_numbers.back());
		}
	}

	for (cluster_iter = m_sub_clusters.begin(); cluster_iter != m_sub_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);
		size += cluster->get_size();
		nPI += cluster->get_nPI();
	}

	if (! m_sub_clusters.empty())
	{
		assert(nPI == get_nPI());
		assert(size == get_size());
	}

}

//
// Prints the distribution of intra- and inter- cluster edges
//
void CLUSTER::print_stats()
{
	CLUSTERS::iterator cluster_iter;
	CLUSTER * cluster = 0;
	NUM_ELEMENTS nIntra_cluster, nInter_cluster;
	NUM_ELEMENTS size;
	double percentage = 0.0;

	for (cluster_iter = m_sub_clusters.begin(); cluster_iter != m_sub_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);

		size = cluster->get_size();
		nIntra_cluster = cluster->get_nIntra_cluster_edges();
		nInter_cluster = cluster->get_nInter_cluster_edges();
		percentage = 100 * static_cast<double>(nIntra_cluster)
							/static_cast<double>(nIntra_cluster + nInter_cluster);
		
		debug("Cluster " << get_info());
		debug("size\t\tnIntra_cluster\t\tnInter_cluster\t\t%intra_cluster");
		debug(size << "\t\t\t" << nIntra_cluster << "       \t\t\t" << nInter_cluster << 
				"       \t\t\t" << percentage << endl);
	}
}

// Get the size (number of nodes and primary inputs) at the specified partitioning depth
// 
// PRE: depth > 0
// RETURNS: the size of the sub-clusters at the specified depth
//
DISTRIBUTION CLUSTER::get_size_distribution
(
	const CLUSTER_NUMBER_TYPE& depth
)
{
	assert(m_type == CLUSTER::REGULAR);

	DISTRIBUTION size_distribution;
	CLUSTERS::iterator cluster_iter;
	CLUSTER * cluster = 0;
	NUM_ELEMENTS size = 0;

	if (depth == 0)
	{
		size = get_size();
		size_distribution.push_back(size);
	}
	else
	{
		for (cluster_iter = m_sub_clusters.begin(); cluster_iter != m_sub_clusters.end(); cluster_iter++)
		{
			cluster = *cluster_iter;
			assert(cluster);

			cluster->add_size_distribution(depth, size_distribution);
		}
	}

	return size_distribution;
}

//
// Adds the size of the clusters at the specified depth to the distribution
//
// PRE: depth > 0
// POST: if this cluster is at the specified depth its size has been added to the distribution
//       otherwise all sub-clusters at the specified depth below this cluster have been added
//       to the distribution
//
void CLUSTER::add_size_distribution
(
	const CLUSTER_NUMBER_TYPE& depth,
	DISTRIBUTION & distribution
)
{
	assert(depth >= m_depth);
	NUM_ELEMENTS size = 0;
	CLUSTERS::iterator cluster_iter;
	CLUSTER * cluster = 0;

	if (depth == m_depth)
	{
		size = get_size();

		distribution.push_back(size);
	}
	else
	{
		for (cluster_iter = m_sub_clusters.begin(); cluster_iter != m_sub_clusters.end(); cluster_iter++)
		{
			cluster = *cluster_iter;
			assert(cluster);

			cluster->add_size_distribution(depth, distribution);
		}
	}
}

// gets the node_shape of the cluster
//
// PRE: nothing
// RETURNS: the node_shape of the cluster
//
SHAPE CLUSTER::get_node_shape() const
{
	SHAPE shape;
	NODES::const_iterator node_iter;
	NODE * node = 0;
	DELAY_TYPE delay = 0;

	if (m_type == REGULAR)
	{
		assert(m_sequential_level);
		shape =  m_sequential_level->get_node_shape();
	}
	else
	{
		shape.resize(m_max_comb_delay+1, 0);
	
		for (node_iter = m_nodes.begin(); node_iter != m_nodes.end(); node_iter++)
		{
			node = *node_iter;
			assert(node);

			delay = node->get_max_comb_delay_level();

			assert(delay >= 0 && delay <= m_max_comb_delay);

			shape[delay]++;
		}
		shape[0] += m_PI.size();
	}

	return shape;
}

// returns a string that identifies the cluster
//
// PRE: nothing
// RETURNS: a string to identifies the cluster
string CLUSTER::get_info() const
{
	string info = util_long_to_string(m_cluster_number) + "  " + util_long_to_string(m_depth);
	CLUSTER_NUMBERS::const_iterator cluster_number_iter;
	
	info += "   " + util_long_to_string(m_cluster_number);

	if (m_type == CLUSTER::SUB)
	{

		for (cluster_number_iter = m_sub_cluster_numbers.begin();
			 cluster_number_iter != m_sub_cluster_numbers.end(); cluster_number_iter++)
		{
			info += " " + util_long_to_string(*cluster_number_iter);
		}
	}

	return info;
}

//
// for a cluster that has been repeatedly divided into sub-clusters
// this function returns a distribution that measures 
// the number of sub-clusters the edges cross
//
// ie. for a circuit that has been divided 3 times
//     we look at the sub-clusters of depth 3.
//     These clusters have sub-cluster number such as
//		   A) 0-2-3, 
//		   B) 0-2-4, 
//		   C) 0-3-3
//		   D) 1-2-3
//
//  From these numbers we can specified the number of clusters an edge crosses.
//
//  An edge that is internal to sub-cluster A has a distance of 0.
//  An edge that goes from sub-cluster A to sub-cluster B has distance of 1 since it 
//     crosses sub-clusters at the LSB.
//  An edge that goes from sub-cluster A to sub-cluster C has distance of 2 since it 
//     crosses sub-clusters one up from the LSB.
//  An edge that goes from sub-cluster A to sub-cluster D has a distance of 2 since it 
//     crosses sub-clusters two up from the LSB.
// 
// 
// PRE: nothing
// RETURNS: inter-cluster edge distance_distribution
//
DISTRIBUTION CLUSTER::get_cluster_distance_distribution() const
{
	DISTRIBUTION cluster_distances(m_max_clustering_depth+2, 0);
	EDGES::const_iterator edge_iter;
	EDGE * edge = 0;
	LENGTH_TYPE inter_cluster_distance = 0;

	// first, take care of the inter cluster edges
	cluster_distances[m_max_clustering_depth+1] = get_nInter_cluster_edges();

	// divide up the intra cluster edges
	for (edge_iter = m_intra_cluster_edges.begin(); edge_iter != m_intra_cluster_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		inter_cluster_distance = edge->get_inter_cluster_distance();
		assert(inter_cluster_distance >= 0 && 
				static_cast<unsigned>(inter_cluster_distance) < cluster_distances.size()-1);
		cluster_distances[inter_cluster_distance]++;
	}

	return cluster_distances;
}

// returns the measurement of Wirelength-approx
//
// PRE: the nodes have their horizontal positions set
// RETURNS: Wirelength-approx of nodes in the cluster
double CLUSTER::get_wirelength_approx() const
{
	NODES::const_iterator node_iter;
	NODE * node = 0;
	NUM_ELEMENTS wirelength_approx = 0;
	double total_wirelength_approx = 0.0;
	
	for (node_iter = m_nodes.begin(); node_iter != m_nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		wirelength_approx = node->get_wirelength_approx_cost();

		total_wirelength_approx += static_cast<double>(wirelength_approx);
	}

	return total_wirelength_approx;
}

//
// returns the fanout values of the primary inputs
//
// PRE: nothing
// RETURNS: pi fanout distribution
//
DISTRIBUTION CLUSTER::get_PI_fanout_values() const
{
	PORTS::const_iterator port_iter;
	PORT * pi = 0;
	DISTRIBUTION pi_distribution;
	NUM_ELEMENTS fanout = 0;

	for (port_iter = m_PI.begin(); port_iter != m_PI.end(); port_iter++) 
	{
		pi = *port_iter;
		assert(pi);

		fanout = pi->get_fanout_degree_to_combinational_nodes();
		assert(fanout >= 0);
		
		pi_distribution.push_back(fanout);
	}

	return pi_distribution;
}

//
// returns the fanout values of the flip-flops
//
// PRE: nothing
// RETURNS: DFF fanout distribution
//
DISTRIBUTION CLUSTER::get_DFF_fanout_values() const
{
	NODES::const_iterator node_iter;
	NODE * node = 0;
	DISTRIBUTION node_distribution;
	NUM_ELEMENTS fanout = 0;

	for (node_iter = m_nodes.begin(); node_iter != m_nodes.end(); node_iter++) 
	{
		node = *node_iter;
		assert(node);

		if (node->get_type() == NODE::SEQ)
		{

			fanout = node->get_fanout_degree_to_combinational_nodes();
			assert(fanout >= 0);
			
			node_distribution.push_back(fanout);
		}
	}

	return node_distribution;
}



// Get rid of all these sub-clusters stuff?

//
// 
// Adds a node to this sub-cluster
//
// PRE: Our cluster is a sub-cluster and the node is valid
// POST: the node, its edges, and ports are added to this cluster
//
void CLUSTER::add_sub_cluster_node
(
	NODE * node
)
{
	assert(node);
	assert(m_type == CLUSTER::SUB);
	assert(node->get_cluster_number() == m_cluster_number);
	assert(node->get_sub_cluster_numbers() == m_sub_cluster_numbers);

	PORT * output_port = node->get_output_port();
	assert(output_port);

	debugif(DCLUSTER, "Adding node " << node->get_info() << " to sub cluster " << get_info());

	m_nodes.push_back(node);
	increment_node_count(node);

	if (output_port->get_type() == PORT::EXTERNAL)
	{
		assert(output_port->get_external_type() == PORT::PO);
		debugif(DCLUSTER, "Adding PO port " << output_port->get_name() << " for node " << node->get_name() <<
				" to sub cluster " << get_info());
		m_PO.push_back(output_port);
	}
}
//
//
// Split the clusters into sub-clusters
//
// PRE: number_times_to_split is positive
// POST: The cluster is partitioned and re-partitioned number_times_to_split times
//
void CLUSTER::split_cluster(NUM_ELEMENTS number_times_to_split)
{
	assert(number_times_to_split >= 0);
	CLUSTERS::iterator cluster_iter;
	CLUSTER * cluster = 0;

	if (number_times_to_split == 0 || get_nNodes() == 0)
	{
		return;
	}
	debug("Subdividing cluster " << get_info());

	if (m_type == CLUSTER::REGULAR)
	{
		m_max_clustering_depth = number_times_to_split;
	}

	NODE_PARTITIONER node_partitioner;

	node_partitioner.partition_cluster(this);

	assert(! m_sub_clusters.empty());
	for (cluster_iter = m_sub_clusters.begin(); cluster_iter != m_sub_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);
		cluster->split_cluster(number_times_to_split-1);
	}
}
//
// Creates the sub-clusters
//
// PRE: number_of_partitions is the number of sub-clusters to be created
// POST: the sub-clusters have been created and the nodes, primary inputs, and edges have
//       been assigned to the sub-clusters
//
void CLUSTER::construct_sub_clusters
(
	const NUM_ELEMENTS & number_of_partitions
)
{
	NUM_ELEMENTS partition_number = 0;
	CLUSTER * cluster = 0;

	for (partition_number = 0; partition_number < number_of_partitions; partition_number++)
	{
		cluster = new CLUSTER(m_cluster_number, m_sub_cluster_numbers, 
								partition_number, m_max_comb_delay, CLUSTER::SUB);
		assert(cluster);
		m_sub_clusters.push_back(cluster);
	}

	divide_nodes_into_clusters();
	divide_primary_inputs_into_clusters();
	divide_intra_cluster_edges_into_clusters();
	divide_inter_cluster_edges_into_clusters();
	print_stats();
}

//
// The nodes in the cluster are sub-divided into sub-clusters
//
// PRE: the sub-clusters have been created
// POST: the cluster's nodes have been assigned to the sub-clusters
//
void CLUSTER::divide_nodes_into_clusters()
{
	NODES::iterator node_iter;
	NODE * node = 0;
	CLUSTER_NUMBER_TYPE cluster_number = 0;

	for (node_iter = m_nodes.begin(); node_iter != m_nodes.end(); node_iter++)
	{
		node = *node_iter;	
		assert(node);

		cluster_number = node->get_last_sub_cluster_number();
		assert(cluster_number >= 0 && cluster_number < static_cast<signed>(m_sub_clusters.size()));

		m_sub_clusters[cluster_number]->add_sub_cluster_node(node);
	}
}
//
// The intra-cluster edges in the cluster are assigned to the sub-clusters
//
// PRE: the sub-clusters have been created
// POST: the cluster's intra-cluster edges have been assigned to the sub-clusters as
//       to intra-sub-cluster and inter-sub-cluster edges 
//
void CLUSTER::divide_intra_cluster_edges_into_clusters()
{
	EDGES::iterator edge_iter;
	EDGE * edge = 0;
	CLUSTER_NUMBER_TYPE source_cluster, sink_cluster, sub_cluster;
	CLUSTER_NUMBERS source_sub_clusters, sink_sub_clusters;

	// divide up the intra cluster edges
	for (edge_iter = m_intra_cluster_edges.begin(); edge_iter != m_intra_cluster_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		source_sub_clusters = edge->get_source_sub_cluster_numbers();
		sink_sub_clusters = edge->get_sink_sub_cluster_numbers();

		assert(edge->is_intra_cluster());

		if (edge->is_intra_cluster_and_intra_sub_cluster())
		{
			assert(! source_sub_clusters.empty() && ! sink_sub_clusters.empty() );
			sub_cluster = sink_sub_clusters.back();
			assert(sub_cluster >= 0 && static_cast<unsigned>(sub_cluster) < m_sub_clusters.size());

			m_sub_clusters[sub_cluster]->add_intra_cluster_edge(edge);
		}
		else
		{
			// else we have an intra_cluster/inter sub cluster edge
			// add the edge to each of the subclusters
			
			assert(! source_sub_clusters.empty() && ! sink_sub_clusters.empty() );

			source_cluster = source_sub_clusters.back();
			sink_cluster = sink_sub_clusters.back();
			assert(source_sub_clusters.size() == sink_sub_clusters.size());
			assert(source_cluster != sink_cluster);

			assert(sink_cluster >= 0 && static_cast<unsigned>(sink_cluster) < m_sub_clusters.size());
			assert(source_cluster >= 0 && static_cast<unsigned>(source_cluster) < m_sub_clusters.size());

			m_sub_clusters[sink_cluster]->add_inter_cluster_input_edge(edge);
			m_sub_clusters[source_cluster]->add_inter_cluster_output_edge(edge);
		}
	}
}

//
// The inter-cluster edges in the cluster are assigned to the sub-clusters
//
// PRE: the sub-clusters have been created
// POST: the cluster's inter-cluster edges have been assigned to the sub-clusters
//
void CLUSTER::divide_inter_cluster_edges_into_clusters()
{
	EDGES::iterator edge_iter;
	EDGE * edge = 0;

	CLUSTER_NUMBERS source_sub_cluster_numbers, sink_sub_cluster_numbers;
	CLUSTER_NUMBER_TYPE source_sub_cluster, sink_sub_cluster;

	// divide up the intra cluster edges
	for (edge_iter = m_inter_cluster_input_edges.begin(); edge_iter != m_inter_cluster_input_edges.end(); 
			edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);
		assert(edge->is_inter_cluster() || edge->is_inter_sub_cluster());
		sink_sub_cluster_numbers = edge->get_sink_sub_cluster_numbers();
		assert(! sink_sub_cluster_numbers.empty());
		sink_sub_cluster = sink_sub_cluster_numbers.back();
		assert(sink_sub_cluster >= 0 && static_cast<unsigned>(sink_sub_cluster) < m_sub_clusters.size());
		m_sub_clusters[sink_sub_cluster]->add_inter_cluster_input_edge(edge);
	}

	for (edge_iter = m_inter_cluster_output_edges.begin(); edge_iter != m_inter_cluster_output_edges.end(); 
			edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);
		assert(edge->is_inter_cluster() || edge->is_inter_sub_cluster());
		
		source_sub_cluster_numbers = edge->get_source_sub_cluster_numbers();
		assert(! source_sub_cluster_numbers.empty());
		source_sub_cluster = source_sub_cluster_numbers.back();
		assert(source_sub_cluster >= 0 && static_cast<unsigned>(source_sub_cluster) < m_sub_clusters.size());
		m_sub_clusters[source_sub_cluster]->add_inter_cluster_output_edge(edge);
	}
}
//
//	assign the primary inputs to clusters.
//	a primary input goes to the clusters that it fansout to the most
//
//	PRE: the sub-clusters have been created
//	POST: the cluster's primary inputs are divided umongst the sub-clusters
//
void CLUSTER::divide_primary_inputs_into_clusters()
{

	CLUSTER_NUMBER_TYPE sub_cluster_number = 0;
	PORTS::iterator pi_iter;
	PORT * pi_port;

	for (pi_iter = m_PI.begin(); pi_iter != m_PI.end(); pi_iter++)
	{
		pi_port = *pi_iter;
		assert(pi_port);

		sub_cluster_number = pi_port->get_last_sub_cluster_number();
		assert(sub_cluster_number >= 0 && static_cast<unsigned>(sub_cluster_number) < m_sub_clusters.size());

		m_sub_clusters[sub_cluster_number]->add_primary_input(pi_port);
	}
}
