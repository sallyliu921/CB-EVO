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





const bool DCP = false;

#include "node_partitioner.h"

// extern "C" int HMETIS_PartKway(int number_of_nodes,int number_of_edges, int * node_weights, 
// 							int * start_of_edges_array, int * hyper_edges, int * edge_weights, 
// 								int number_parts, int ubfactor, int * options, int * results, int* edge_cut);

// extern "C" int HMETIS_PartRecursive(int number_of_nodes,int number_of_edges, int * node_weights, 
// 								int * start_of_edges_array, int * hyper_edges, int * edge_weights, 
// 									int number_parts, int ubfactor, int * options, int * results, int* edge_cut);


NODE_PARTITIONER::NODE_PARTITIONER()
{
	m_circuit = 0;
	m_cluster = 0;
}
NODE_PARTITIONER::NODE_PARTITIONER(const NODE_PARTITIONER & another_node_partitioner)
{
	m_circuit = another_node_partitioner.m_circuit;
	m_cluster = another_node_partitioner.m_cluster;
}

NODE_PARTITIONER & NODE_PARTITIONER::operator=(const NODE_PARTITIONER & another_node_partitioner)
{
	m_circuit = another_node_partitioner.m_circuit;
	m_cluster = another_node_partitioner.m_cluster;

	return (*this);
}


NODE_PARTITIONER::~NODE_PARTITIONER()
{
}

// 
// Partition the cluster
// PRE:  cluster is valid
// POST: the cluster has been partioned and sub-clusters have been created
//
COST_TYPE NODE_PARTITIONER::partition_cluster(CLUSTER * cluster)
{
	assert(cluster);
	m_cluster = cluster;
	m_partition_type = NODE_PARTITIONER::CLUSTER_TYPE;

	COST_TYPE cost = partition_graph();

	m_cluster = 0;

	return cost;

}
// 
// Partition the circuit
// PRE:  circuit is valid
// POST: the circuit has been partioned and clusters have been created
//
COST_TYPE NODE_PARTITIONER::partition_circuit(CIRCUIT * circuit)
{
	assert(circuit);
	m_circuit = circuit;
	m_partition_type = NODE_PARTITIONER::CIRCUIT_TYPE;

	COST_TYPE cost = partition_graph();

	m_circuit = 0;

	return cost;
}

//
// Returns: the size (number of nodes + number of primary inputs)
// 
int NODE_PARTITIONER::get_size()
{
	int size = 0;

	if (m_partition_type == NODE_PARTITIONER::CIRCUIT_TYPE)
	{
		size = m_circuit->get_size();
	}
	else
	{
		assert(m_partition_type == NODE_PARTITIONER::CLUSTER_TYPE);

		size = m_cluster->get_size();
	}

	return size;
}

//
// Returns: the number of edges (without the edges connecting the clock to the flip-flops)
//
int NODE_PARTITIONER::get_nEdges()
{
	int nEdges = 0;

	if (m_partition_type == NODE_PARTITIONER::CIRCUIT_TYPE)
	{
		nEdges = m_circuit->get_nEdges_without_clock_edges();
	}
	else
	{
		assert(m_partition_type == NODE_PARTITIONER::CLUSTER_TYPE);

		// might want to remove the clock edges from here as well
		assert(Dlook_at);
		nEdges = m_cluster->get_nIntra_cluster_edges();
	}

	return nEdges;
}

//
// Construct the clusters that resulted from the partitioning
// PRE: each node and primary input has its cluster labelled
// POST: the clusters have been created
//
void NODE_PARTITIONER::construct_clusters
(
	const NUM_ELEMENTS & number_partitions
)
{
	if (m_partition_type == NODE_PARTITIONER::CIRCUIT_TYPE)
	{
		m_circuit->construct_clusters(number_partitions);
	}
	else
	{
		assert(m_partition_type == NODE_PARTITIONER::CLUSTER_TYPE);
		m_cluster->construct_sub_clusters(number_partitions);
	}
}

// 
// Partition the graph
// PRE: circuit or cluster are valid
// POST: the circuit/cluster has been partioned and clusters/sub-clusters have been created
//
COST_TYPE NODE_PARTITIONER::partition_graph()
{
	int number_of_nodes = 0;
	int number_of_edges = 0;
	int * node_weights = 0;
	int * start_of_edges_array = 0;
	int * edges_info = 0;

	int * edge_weights = 0;	
	int ubfactor = g_options->get_ub_factor();
	// ubfactor: balancing factor. defined differently for k-way vs. bi-partitioning
		
	int options[9] = {0,0,0,0,0,0,0,0,0};
	int * results;			// where the results of the partitioning are returned
	int k_way_edge_cut = 0; // returns the number of hyperedges cut
	int rec_edge_cut = 0; // returns the number of hyperedges cut
	COST_TYPE k_way_scaled_cost = 0.0;
	COST_TYPE recursive_scaled_cost = 0.0;
	COST_TYPE cost = 0.0;


	debug("ubfactor is " << ubfactor);

	int number_partitions = g_options->get_nPartitions();

	// if the user wants a single cluster then 
	// just create create the cluster and return
	if (number_partitions == 1 || number_partitions == 0)
	{
		construct_single_cluster();
		return 0;
	}

	number_of_nodes = get_size();
	number_of_edges = get_nEdges();


	// Maps each output port to a number
	// This will allow us to use the hMetis API
	create_port_ptr_to_index_map();

	// create the data structures for hMetis
	start_of_edges_array = new int[number_of_edges+1];
	edges_info = new int[number_of_edges*2];
	results = new int[number_of_nodes];

	edge_weights = 0; 

	// fill in the data structures we pass to hMetis
	fill_in_edge_info(start_of_edges_array, edges_info);

	// Call the correct partitioner
	if (g_options->get_type_of_partitioning() == OPTIONS::KWAY)
	{
		debug("The total scaled cost \tk_way_scaled_cost");

		// disabled for 64bit usage 
		// HMETIS_PartKway(number_of_nodes, number_of_edges, node_weights, 
		// 	start_of_edges_array, edges_info, edge_weights, 
		// 	number_partitions, ubfactor, options, results, &k_way_edge_cut);

		k_way_scaled_cost = get_scaled_cost(number_of_edges, number_partitions, 
												number_of_nodes, edges_info, results);
		cost = k_way_scaled_cost;

		debug(number_partitions << "\t" << k_way_scaled_cost);
	}
	else
	{
		assert(g_options->get_type_of_partitioning() == OPTIONS::RECURSIVE_BI);

		debug("The total scaled cost \trecursive_scaled_cost");

		// HMETIS_PartRecursive(number_of_nodes, number_of_edges, node_weights, 
		// 	start_of_edges_array, edges_info, edge_weights, 
		// 	number_partitions, ubfactor, options, results, &rec_edge_cut);

		recursive_scaled_cost = get_scaled_cost(number_of_edges, number_partitions, 
												number_of_nodes, edges_info, results);
		cost = k_way_scaled_cost;

		debug(number_partitions << "\t" << recursive_scaled_cost);
	}

	// from the results we obtained from hMetis label each node with the partition they belong to
	label_nodes_with_the_partition_they_belong_to(results, number_of_nodes, number_partitions); 

	// construct the clusters 
	construct_clusters(number_partitions);

	m_circuit->set_scaled_cost(cost);

	delete [] results;
	results = 0;
	delete [] start_of_edges_array;
	start_of_edges_array = 0;
	delete [] edges_info;
	edges_info = 0;

	return cost;

}


// 
// Calculate the scaled cost for the partition we are given
// PRE: the results array has results
// RETURNS: the scaled cost  = sum over all clusters of (number of edges cut)/(number of nodes)
//
COST_TYPE NODE_PARTITIONER::get_scaled_cost
(
	const int & number_of_edges,
	const int & number_partitions,
	const int & number_of_nodes,
	int * edges_info,
	int * results
)
{
	COST_TYPE total_cost = 0.0;
	int edge_index = 0;
	int partition_index = 0;
	int edges_cut_for_partition = 0;
	int vertex_A, vertex_B;

	PARTITION_WEIGHTS_TYPE_ALSO partition_weights = get_partition_weights(results, number_of_nodes, 
																	number_partitions);
	
	// sum the cost over all the partitions
	for (partition_index = 0; partition_index < number_partitions; partition_index++)
	{
		edges_cut_for_partition = 0;

		debugif(DCP, "Looking at partition " << partition_index);
		// check all the edges to see if any of them cut the partition
		for (edge_index = 0; edge_index < number_of_edges; edge_index++)
		{
			vertex_A = edges_info[edge_index*2];
			vertex_B = edges_info[edge_index*2+1];

			if (results[vertex_A] == partition_index || results[vertex_B] == partition_index)
			{
				// check to see if the edge is cuts the partition
				if ((results[vertex_A] != partition_index || results[vertex_B] != partition_index))
				{
					debugif(DCP, "\tFound a edge " << vertex_A << " to " << vertex_B);
					edges_cut_for_partition++;
				}
			}
		}	
		// note some partitions have zero weight. this will cause the cost to have infinite cost
		total_cost += (static_cast<COST_TYPE>(edges_cut_for_partition))/
					 (static_cast<COST_TYPE>(partition_weights[partition_index]));
	}	

	total_cost = total_cost/static_cast<COST_TYPE>(number_partitions-1);

	return total_cost;
}		


		


	
//
// fill in a data structures we pass to hMetis that details what connects to what
//  If we are partitioning a circuit all edges (sans edges connected to clocks) are dealt with
//  If we are partitioning a cluster only intra-cluster edges are dealt with
//
//  PRE:   m_port_index_map has been filled in
//  POST:  edges_info and start_of_edges_array has been filled in
//
void NODE_PARTITIONER::fill_in_edge_info
(
	int * start_of_edges_array,
	int * edges_info
)
{
	assert(start_of_edges_array && edges_info);
	EDGES edges;
	EDGES::iterator edge_iter;

	// if we are partioning a circuit get all the edges, otherwise
	// get the edges internal to the cluster
	if (m_partition_type == NODE_PARTITIONER::CIRCUIT_TYPE)
	{
		edges = m_circuit->get_edges();
	}
	else
	{
		assert(m_partition_type == NODE_PARTITIONER::CLUSTER_TYPE);
		edges = m_cluster->get_intra_cluster_edges();
	}

	PORT * source_output_port = 0;
	NODE * sink_node = 0;
	PORT * sink_output_port = 0;
	int edge_index =  0;
	EDGE * edge = 0;
	NUM_ELEMENTS number_of_edges = get_nEdges();

	// for each edge write to the hMetis datastructure what each connects to
	debugif(DCP, "About to fill in edge information");
	for (edge_iter = edges.begin(); edge_iter != edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		source_output_port = edge->get_source();
		sink_node = edge->get_sink_node();
		assert(sink_node);
		sink_output_port = sink_node->get_output_port();
		assert(source_output_port && sink_output_port);

		// if both are ports are in the port map, add the edge info to the edges_info array
		if (m_port_index_map.find(source_output_port) != m_port_index_map.end() &&
			m_port_index_map.find(sink_output_port) != m_port_index_map.end())
		{

			edges_info[edge_index*2] = m_port_index_map[source_output_port];
			edges_info[edge_index*2+1] = m_port_index_map[sink_output_port];

			// set the edge ptr to point to the location of the vertice listing as per the hmetis doc
			start_of_edges_array[edge_index] = edge_index*2;
			edge_index++;
		}
		else
		{
			// we get here if the edge was a clock edge or if the edge was inter-cluster
			// an inter-cluster edge will not have one of its port in the port map
			if (source_output_port->get_io_direction() != PORT::CLOCK)
			{
				debug("Didn't find one or both of these ports " << source_output_port->get_name() << " and " 
						<< sink_output_port->get_name() << " in the list of external ports");
			}
		}
	}	

	assert(edge_index == number_of_edges);

	// close off the start_of_edges_array with the one past the last location of the edges_info
	start_of_edges_array[number_of_edges] = number_of_edges*2;

	debugif(DCP, "Finished filling in edge information");
}




//
// Calculates the weight of each cluster
// PRE: results has the results of partitioning
// RETURNS: an array that contains the weight of each cluster in terms of the number of nodes 
//
PARTITION_WEIGHTS_TYPE_ALSO NODE_PARTITIONER::get_partition_weights
(
	int * results,
	const int & number_of_nodes,
	const int & number_of_partitions
)
{
	PARTITION_WEIGHTS_TYPE_ALSO partition_weights(number_of_partitions,0);
	int partition_number = 0;
	int node_number = 0;
	int total_partition_weight = 0;

	for (node_number = 0; node_number < number_of_nodes; node_number++)
	{
		partition_number = results[node_number];
		assert(partition_number >= 0 && partition_number < number_of_partitions);

		partition_weights[partition_number]++; // each node has weight of 1

	}

	for (partition_number = 0; partition_number < number_of_partitions; partition_number++)
	{
		total_partition_weight += partition_weights[partition_number];
	}
	
	assert(total_partition_weight == number_of_nodes);

	return partition_weights;
}

//
// Map each output port in the circuit or cluster to a number
//
// POST: m_port_index_map contains a map of all relevant output ports to an index
//
void NODE_PARTITIONER::create_port_ptr_to_index_map()
{
	assert(m_circuit || m_cluster);
	NODES nodes;
	PORTS pi_ports;
	NODES::iterator node_iter;
	PORTS::iterator port_iter;
	NODE * node = 0;
	PORT * port = 0;
	int index = 0;


	// depending on whether we are partioning the whole circuit are just a cluster
	// get the number of nodes and PI we are partitioning
	if (m_partition_type == NODE_PARTITIONER::CIRCUIT_TYPE)
	{
		nodes = m_circuit->get_nodes();
		pi_ports = m_circuit->get_PI();
	}
	else
	{
		assert(m_partition_type == NODE_PARTITIONER::CLUSTER_TYPE);
		nodes = m_cluster->get_nodes();
		pi_ports = m_cluster->get_PI();
	}


	// map the nodes
	for (node_iter = nodes.begin(); node_iter != nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);
		port = node->get_output_port();
		assert(port);

		m_port_index_map[port] = index;
		index++;
	}

	// map the PIs
	for (port_iter = pi_ports.begin(); port_iter != pi_ports.end(); port_iter++)
	{
		port = *port_iter;
		assert(port);
		m_port_index_map[port] = index;
		index++;
	}
	assert(index = get_size());
}

//
// label each node with the partition they belong to 
// PRE: hMetis has been called. results array contains the results
// POST: each node and primary output has its partition labelled
void NODE_PARTITIONER::label_nodes_with_the_partition_they_belong_to
(
	int * results,
	const int & number_of_nodes,
	const int & number_of_partitions
)
{
	assert(m_circuit || m_cluster);
	NODES nodes;
	PORTS pi_ports;
	NODES::iterator node_iter;
	PORTS::iterator port_iter;
	NODE * node = 0;
	PORT * port = 0;
	CLUSTER_NUMBER_TYPE partition_number = 0;
	int node_number = 0;

	if (m_partition_type == NODE_PARTITIONER::CIRCUIT_TYPE)
	{
		nodes = m_circuit->get_nodes();
		pi_ports = m_circuit->get_PI();
	}
	else
	{
		assert(m_partition_type == NODE_PARTITIONER::CLUSTER_TYPE);
		nodes = m_cluster->get_nodes();
		pi_ports = m_cluster->get_PI();
	}


	for (node_iter = nodes.begin(); node_iter != nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);
		port = node->get_output_port();
		assert(port);

		node_number = m_port_index_map[port];
		partition_number = results[node_number];

		assert(partition_number >= 0 && partition_number < number_of_partitions);

		if (m_partition_type == NODE_PARTITIONER::CIRCUIT_TYPE)
		{
			node->set_cluster_number(partition_number);
		}
		else
		{
			assert(m_partition_type == NODE_PARTITIONER::CLUSTER_TYPE);
			node->add_sub_cluster_number(partition_number);
		}
	}

	for (port_iter = pi_ports.begin(); port_iter != pi_ports.end(); port_iter++)
	{
		port = *port_iter;
		assert(port);

		node_number = m_port_index_map[port];
		partition_number = results[node_number];

		assert(partition_number >= 0 && partition_number < number_of_partitions);

		if (m_partition_type == NODE_PARTITIONER::CIRCUIT_TYPE)
		{
			port->set_cluster_number(partition_number);
		}
		else
		{
			assert(m_partition_type == NODE_PARTITIONER::CLUSTER_TYPE);
			port->add_sub_cluster_number(partition_number);
		}
	}
}

// 
// Creates a single cluster
//
// PRE:  circuit is valid
// POST: each node and primary input is labeled with the cluster number 0
// 		 a single cluster has been created in circuit
//
void NODE_PARTITIONER::construct_single_cluster()
{
	assert(m_circuit);
	NODES nodes = m_circuit->get_nodes();
	PORTS pi_ports =  m_circuit->get_PI();
	NODES::iterator node_iter;
	PORTS::iterator port_iter;
	NODE * node = 0;
	PORT * port = 0;

	for (node_iter = nodes.begin(); node_iter != nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		node->set_cluster_number(0);
	}

	for (port_iter = pi_ports.begin(); port_iter != pi_ports.end(); port_iter++)
	{
		port = *port_iter;
		assert(port);

		port->set_cluster_number(0);
	}

	m_circuit->construct_clusters(1);
}
