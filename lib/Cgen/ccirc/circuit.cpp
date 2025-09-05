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





#include "circuit.h"
#include "graph_medic.h"
#include "drawer.h"
#include <algorithm>
using namespace std;

CIRCUIT::CIRCUIT()
{
	m_name				= "";

	m_global_clock		= 0;

	m_max_fan_in		= 0;

	m_number_seq_nodes 	= 0;
	m_number_comb_nodes = 0;
	m_number_edges		= 0;


	m_degree_info		= 0;

	m_sequential_level	=  new SEQUENTIAL_LEVEL(0);
	assert(m_sequential_level);
}

CIRCUIT::CIRCUIT(const CIRCUIT & another_circuit)
{
	assert(false);
	m_name				= another_circuit.m_name;
	m_max_fan_in		= another_circuit.m_max_fan_in;
	m_PI				= another_circuit.m_PI;
	m_PO				= another_circuit.m_PO;
	m_global_clock		= another_circuit.m_global_clock;
	m_nodes				= another_circuit.m_nodes;
	m_edges				= another_circuit.m_edges;
	m_sequential_level	= another_circuit.m_sequential_level;
	m_inter_cluster_connections = another_circuit.m_inter_cluster_connections;

	m_number_seq_nodes 	= another_circuit.m_number_seq_nodes;
	m_number_comb_nodes = another_circuit.m_number_comb_nodes;
	m_number_edges		= another_circuit.m_number_edges;

	m_degree_info		= another_circuit.m_degree_info;
}
CIRCUIT & CIRCUIT::operator=(const CIRCUIT & another_circuit)
{
	assert(false);
	m_name				= another_circuit.m_name;
	m_max_fan_in		= another_circuit.m_max_fan_in;
	m_PI				= another_circuit.m_PI;
	m_PO				= another_circuit.m_PO;
	m_global_clock		= another_circuit.m_global_clock;
	m_nodes				= another_circuit.m_nodes;
	m_edges				= another_circuit.m_edges;
	m_sequential_level	= another_circuit.m_sequential_level;
	m_inter_cluster_connections = another_circuit.m_inter_cluster_connections;

	m_number_seq_nodes 	= another_circuit.m_number_seq_nodes;
	m_number_comb_nodes = another_circuit.m_number_comb_nodes;
	m_number_edges		= another_circuit.m_number_edges;

	m_degree_info		= another_circuit.m_degree_info;

	return (*this);
}
CIRCUIT::~CIRCUIT()
{
	// must delete the graphs
	NODES::iterator node_iter;
	NODE * node;
	CLUSTERS::iterator cluster_iter;
	CLUSTER * cluster;

	for (node_iter = m_nodes.begin(); node_iter != m_nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);
		delete node;
	}

	delete m_degree_info;
	m_degree_info = 0;
	
	delete m_sequential_level;
	m_sequential_level = 0;

	for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);

		delete cluster;
	}
}

//
// RETURNS: all the PI with the clock PI
//
PORTS CIRCUIT::get_PI_with_clock()
{
	// terribly inefficient but it works
	PORTS pi = m_PI;

	if (m_global_clock)
	{
		pi.push_back(m_global_clock);
	}

	return pi;
}

//
// creates and adds a external port with name 
// port_name of type external_type
//
// PRE: port_name contains the name of the port
//      external_type is the type of port
// POST: an external port has been created
// RETURNS: the external port created
//
PORT *	CIRCUIT::create_and_add_external_port
(
	const string & port_name,
	const PORT::EXTERNAL_TYPE & external_type
)
{
	assert(external_type != PORT::NONE);

	// PI are outputs because they fanout
	// PO are outputs because they fanout and also because
	// making them outputs allows them to glue easily to the 
	// PI of other graphs
	
	PORT * port = new PORT(port_name, PORT::EXTERNAL, PORT::OUTPUT, external_type);
	assert(port);

	if (external_type == PORT::PI || external_type == PORT::GI)
	{
		m_PI.push_back(port);
	}
	else 
	{
		m_PO.push_back(port);
	}
		
	return port;
}

//
// adds the external port 
//
// PRE: port is valid
//      external_type is the type of port
// POST: the external port has been added
// RETURNS: the external port created
//
PORT *	CIRCUIT::add_external_port
(
	PORT * port,
	const PORT::EXTERNAL_TYPE & external_type
)
{
	assert(external_type != PORT::NONE);
	assert(port);

	// PI are outputs because they fanout
	// PO are outputs because they fanout and also because
	// making them outputs allows them to glue easily to the 
	// PI of other graphs
	
	// I haven't fully testing for external GI and GO.
	// note: I think for GI and GO you need at least a delay level they are attached to
	assert(external_type == PORT::PI || external_type == PORT::PO);

	if (external_type == PORT::PI || external_type == PORT::GI)
	{
		m_PI.push_back(port);
	}
	else 
	{
		m_PO.push_back(port);
	}
		
	return port;
}

//
// remove but do not delete the external port
// 
// PRE: port_to_remove is valid and external
// POST: the port has been removed from the list of external ports
//
void  CIRCUIT::remove_external_port
(
	PORT * port_to_remove
)
{
	assert(port_to_remove->get_type() == PORT::EXTERNAL);

	PORTS::iterator port_to_remove_iter;
	
	if (port_to_remove->get_external_type() == PORT::PI)
	{
		port_to_remove_iter = find(m_PI.begin(), m_PI.end(), port_to_remove);
		assert((*port_to_remove_iter) == port_to_remove );

		m_PI.erase(port_to_remove_iter);
	}
	else
	{
		assert(port_to_remove->get_external_type() == PORT::PO);

		// the clock port should not be still in the list of external ports
		// because it should have been recognized by set_global_clock
		debug("Removing external port " << port_to_remove->get_name());
		port_to_remove_iter = find(m_PO.begin(), m_PO.end(), port_to_remove);
		assert((*port_to_remove_iter) == port_to_remove );

		m_PO.erase(port_to_remove_iter);

	}	

	if (port_to_remove == m_global_clock)
	{
		m_global_clock = 0;
	}
}
//
// PRE: node is valid
// POST: the node has been removed from the list of nodes
//       the appropiate node count has been decremented
//
void CIRCUIT::remove_from_node_list
(
	NODE * node
)
{
	m_nodes.erase(remove(m_nodes.begin(), m_nodes.end(), node), m_nodes.end());
	decrement_node_count(node);
}

//
// PRE: new_global_clock is the global clock
// POST: new_global_clock has been designated as a clock
// 		 the global_clock has been removed from the list of primary inputs
// 		 m_global_clock points to the global clock
//
void CIRCUIT::set_global_clock
(
	PORT * new_global_clock
)
{
	assert(m_global_clock == 0);
	m_global_clock = new_global_clock;
	PORTS::iterator port_iter;

	// set the direction as clock. the direction is probably INPUT
	m_global_clock->set_direction(PORT::CLOCK);

	// find the clock in the primary inputs and delete it 
	// because the clock port is considered a special "primary input"
	port_iter = find(m_PI.begin(), m_PI.end(), m_global_clock);
	assert(*port_iter == m_global_clock);

	m_PI.erase(port_iter);
}

// sets the list of primary inputs or primary outputs to be 
// new ports
//
// PRE: new_ports contains a list of ports
//      io_direction is their direction
// POST: m_PI or m_PO is now new_ports
//
void CIRCUIT::set_external_ports
(	
	PORTS & new_ports,
	PORT::IO_DIRECTION io_direction
)
{
	if (io_direction == PORT::INPUT)
	{
		m_PI = new_ports;
	}
	else
	{
		m_PO = new_ports;
	}
}

//
// Do a DFS from the primary outputs, marking all nodes with the colour.  
// 
// PRE: all paths can be reach from the POs
// POST: the graph has been coloured
//
void CIRCUIT::colour_graph
(
	const NODE::COLOUR_TYPE & colour
)
{

	debugif(DCOLOUR, "Reseting the graph to colour " << static_cast<short>(colour));

	// first colour the nodes to NONE to allow all the nodes to be coloured
	// in the second part
	colour_nodes(NODE::NONE);
	
	// Now colour the nodes the colour the user wants
	colour_nodes(colour);
}


//
// Do a DFS from the primary outputs, marking all nodes with the colour.  
// 
// PRE: all paths can be reach from the POs
// POST: the graph has been coloured
//
void CIRCUIT::colour_nodes
(
	const NODE::COLOUR_TYPE & colour
)
{
	PORTS::iterator port_iter;
	PORT * output_port;
	NODE * node;
	NUM_ELEMENTS number_nodes = get_nNodes();
	
	for (port_iter = m_PO.begin(); port_iter != m_PO.end(); port_iter++)
	{
		output_port = *port_iter;
		assert(output_port);

		node = output_port->get_my_node();
	
		if (node->get_colour() != colour)
		{
			colour_up_from_node(node, colour, number_nodes);
		}
	}
	assert(number_nodes == 0);
}

// colour up from this node
//
// PRE: node is valid,
//      number_nodes contains the number of nodes still to colour
//      all paths can be reach from the POs
// POST: all nodes up from this node have been coloured
//
void CIRCUIT::colour_up_from_node
(
	NODE * node,
	const NODE::COLOUR_TYPE & colour,
	NUM_ELEMENTS & number_nodes
)
{
	PORTS	input_ports;
	NODE *			output_node;
	PORTS::iterator port_iter;

	assert(node);
	node->set_colour(colour);
	number_nodes--;

	debugif(DCOLOUR, node->get_name() << " = marked " << static_cast<short>(colour) );

	input_ports = node->get_input_ports();

	for (port_iter = input_ports.begin(); port_iter != input_ports.end(); port_iter++) 
	{
		assert(*port_iter);
		output_node = (*port_iter)->get_node_that_fanout_to_me();

		if (output_node && output_node->get_colour() != colour)
		{
			colour_up_from_node(output_node, colour, number_nodes);
		}
	}
}

// create a dff 
// only to be used once the graph is complete
//
// PRE: new_dff_name contains the new dff name
//      input_port_name contains the name of node that feeds the dff
// POST: the dff has been created and added to the circuit
// RETURNS: the new dff
// 
NODE * CIRCUIT::create_dff
(
	const string & new_dff_name, 
	const string & input_port_name
)
{
	NODE * node = 0;
	PORT * input_port = 0;
	PORT * clock_port = 0;

	debug("Creating a sequential node with node name : " << new_dff_name << " and input name " << input_port_name);

	if (! m_global_clock)
	{
		Fail("Trying to create dff.  No global clock exists.  Please create one.");
	}

	node = create_node(new_dff_name, NODE::SEQ);

	node->create_and_add_port(new_dff_name, PORT::INTERNAL, PORT::OUTPUT, PORT::NONE);

	input_port = node->create_and_add_port(new_dff_name, PORT::INTERNAL, PORT::INPUT, PORT::NONE);

	clock_port = node->create_and_add_port(m_global_clock->get_name(), 
											PORT::INTERNAL, PORT::CLOCK, PORT::NONE);

	// connect the global clock to the clock port
	create_edge(m_global_clock, clock_port, 1);


	return node;
}

//
// creates an edge
//
// PRE: source_port and sink_port are the ports
//      and are valid
//      length is the length of the edge
// POST: the edge has been created
// RETURNS: the edge
//
EDGE * CIRCUIT::create_edge
(
	PORT * source_port,
	PORT * sink_port,
	const LENGTH_TYPE & length
)
{
	assert(source_port && sink_port);

	EDGE * new_edge = new EDGE(source_port, sink_port, length);
	assert(new_edge);

	source_port->add_edge(new_edge);
	sink_port->set_edge(new_edge);

	m_edges.push_back(new_edge);
	m_number_edges++;

	return new_edge;
}

//
// deletes an edge
//
// PRE: source_port and sink_port are the edge's ports
//      edge_to_delete is the edge to delete 
// POST: the edge has been erased from the list of edges and deleted
//
void CIRCUIT::remove_and_delete_edge
(
	PORT * source_port,
	PORT * sink_port,
	EDGE * edge_to_delete
)
{
	assert(source_port);
	assert(sink_port);
	assert(sink_port->get_edge() == edge_to_delete);

	source_port->remove_edge(edge_to_delete);
	sink_port->remove_edge(edge_to_delete);
	
	debug("CIRCUIT::remove_and_delete_edge == DEBUG");
	m_edges.erase(remove(m_edges.begin(), m_edges.end(), edge_to_delete), m_edges.end());
	debug("CIRCUIT::remove_and_delete_edge == DEBUG");
	delete edge_to_delete;
	debug("CIRCUIT::remove_and_delete_edge == DEBUG");
	edge_to_delete = 0;
	debug("CIRCUIT::remove_and_delete_edge == DEBUG");

	m_number_edges--;
	debug("CIRCUIT::remove_and_delete_edge == DEBUG");
}

//
// creates a node
//
// PRE: node_name is the name of the node
//      node_type is the type of the node
// POST: the node without input and output ports has been created
//       a node counter have been incremented
// RETURNS: the node
//
NODE * CIRCUIT::create_node
(
	const string & node_name,
	const NODE::NODE_TYPE & node_type
)
{
	NODE * node;
	debugif(DCODE, "Creating a node with name = " << node_name);

	node = new NODE(node_name, node_type);
	assert(node);

	increment_node_count(node);
	m_nodes.push_back(node);

	return node;
}

// PRE: node is valid
// POST: m_number_comb_nodes or m_number_seq_nodes 
//       has been decremented
//
void CIRCUIT::decrement_node_count
(
	NODE * node
)
{
    assert(node);
	assert(node->get_type() == NODE::COMB || node->get_type() == NODE::SEQ);

    if (node->get_type() == NODE::COMB)
    {
        m_number_comb_nodes--;
    }
    else
    {
        m_number_seq_nodes--;
	}
}

// PRE: node is valid
// POST: m_number_comb_nodes or m_number_seq_nodes 
//       has been incremented
//
void CIRCUIT::increment_node_count
(
	NODE * node
)
{
    assert(node);
	assert(node->get_type() == NODE::COMB || node->get_type() == NODE::SEQ);

    if (node->get_type() == NODE::COMB)
    {
        m_number_comb_nodes++;
    }
    else
    {
        m_number_seq_nodes++;
	}
}

//
//	top nodes are nodes connected to PI or DFF.
// 
//  RETURNS: all the nodes connected to either a PI or DFF
//
NODES CIRCUIT::get_top_nodes()
{
	NODES top_nodes;
	NODES::iterator node_iter;
	NODES::iterator source_node_iter;
	PORTS::iterator pi_iter;
	EDGES edges;
	EDGES::iterator edge_iter;
	PORT * pi_port;
	NODE * node = 0;
	NODE * source_node = 0;
	EDGE * edge;
	NODES dffs = get_dffs();

	for (pi_iter = m_PI.begin(); pi_iter != m_PI.end(); pi_iter++)
	{
		pi_port = *pi_iter;
		assert(pi_port);

		edges = pi_port->get_edges();

		for (edge_iter = edges.begin(); edge_iter != edges.end(); edge_iter++)
		{
			edge = *edge_iter;
			assert(edge);
			node = edge->get_sink_node();
			assert(node);

			node_iter = find(top_nodes.begin(), top_nodes.end(), node);

			if (node_iter == top_nodes.end())
			{
				// we didn't find the node so add it to the list
				debugif(DCODE, "Nodes connected to PI: Adding " << node->get_name());
				top_nodes.push_back(node);
			}
		}
	}

	for (source_node_iter = dffs.begin(); source_node_iter != dffs.end(); source_node_iter++)
	{
		source_node = *source_node_iter;
		assert(source_node);

		edges = source_node->get_output_edges();

		for (edge_iter = edges.begin(); edge_iter != edges.end(); edge_iter++)
		{
			edge = *edge_iter;
			assert(edge);
			node = edge->get_sink_node();
			assert(node);

			node_iter = find(top_nodes.begin(), top_nodes.end(), node);

			if (node_iter == top_nodes.end())
			{
				// we didn't find the node so add it to the list
				debugif(DCODE, "Nodes connected to dff: Adding " << node->get_name());
				top_nodes.push_back(node);
			}
		}
	}

	return top_nodes;
}



//
// POST: information on the nodes has been printed out
//
void CIRCUIT::print_out_node_information() const
{
	NODES::const_iterator node_iter;
	NODE * node;

	for (node_iter = m_nodes.begin(); node_iter != m_nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);
		node->print_out_information();
	}
}



//
// RETURNS: the maximum combinational delay of the circuit
//
DELAY_TYPE CIRCUIT::get_maximum_combinational_delay()
{
	assert(m_sequential_level);

	DELAY_TYPE maximum_comb_delay = m_sequential_level->get_maximum_combinational_delay();

	return maximum_comb_delay;
}


// create the clusters and divides the nodes and primary inputs amongst the clusters
//
// PRE: the nodes have a cluster number and a delay level
// POST: the clusters have been created and assigned PI and nodes
// 
void CIRCUIT::construct_clusters
(
	const NUM_ELEMENTS & number_of_partitions
)
{
	NUM_ELEMENTS partition_number = 0;
	CLUSTER * cluster = 0;
	//NUM_ELEMENTS number_of_sequential_levels = get_nSeq_levels();
	//DELAYS max_delay_by_seq_level = get_max_delay_by_seq_level();
	DELAY_TYPE max_comb_delay = get_maximum_combinational_delay();

	debug("Status: Constructing " << number_of_partitions << " clusters");
	for (partition_number = 0; partition_number < number_of_partitions; partition_number++)
	{
		cluster = new CLUSTER(partition_number, max_comb_delay, CLUSTER::REGULAR);
		assert(cluster);
		m_clusters.push_back(cluster);
	}

	divide_nodes_into_clusters();
	divide_primary_inputs_into_clusters();
	print_cluster_stats();

	create_inter_cluster_matrix();
}

// divides the nodes amongst the clusters
//
// PRE: the nodes have a cluster number and a delay level
//      the clusters have been created
// POST: the nodes have been assigned to a cluster
// 
void CIRCUIT::divide_nodes_into_clusters()
{
	NODES::iterator node_iter;
	NODE * node = 0;
	CLUSTER_NUMBER_TYPE cluster_number = 0;

	for (node_iter = m_nodes.begin(); node_iter != m_nodes.end(); node_iter++)
	{
		node = *node_iter;	
		assert(node);

		cluster_number = node->get_cluster_number();
		assert(cluster_number >= 0 && cluster_number < static_cast<signed>(m_clusters.size()));

		m_clusters[cluster_number]->add_node(node);
	}
}


//
//	assign the primary inputs to clusters.
//
// PRE: the primary inputs have a cluster number and a delay level
//      the clusters have been created
// POST: the primary inputs have been assigned to a cluster
// 
void CIRCUIT::divide_primary_inputs_into_clusters()
{

	CLUSTER_NUMBER_TYPE cluster_number = 0;
	PORTS::iterator pi_iter;
	PORT * pi_port;

	for (pi_iter = m_PI.begin(); pi_iter != m_PI.end(); pi_iter++)
	{
		pi_port = *pi_iter;
		assert(pi_port);

		cluster_number = pi_port->get_cluster_number();
		assert(cluster_number >= 0 && static_cast<unsigned>(cluster_number) < m_clusters.size());
		m_clusters[cluster_number]->add_primary_input(pi_port);

	}
}

//
// POST: statistics on the clusters 
//	     (size/number of intra cluster edges/number of inter-cluster edges)
//	     have been printed out
//
void CIRCUIT::print_cluster_stats()
{
	CLUSTERS::iterator cluster_iter;
	CLUSTER * cluster = 0;
	NUM_ELEMENTS nIntra_cluster, nInter_cluster;
	NUM_ELEMENTS size;

	for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);

		size = cluster->get_size();
		nIntra_cluster = cluster->get_nInter_cluster_edges();
		nInter_cluster = cluster->get_nInter_cluster_edges();
		
		debug("Cluster " << cluster->get_cluster_number());
		debug("size\tnIntra_cluster\tnInter_cluster");
		debug(size << "   \t" << nIntra_cluster << "     \t\t" << nInter_cluster << endl);
	}
}

//
// Partitions the clusters into sub-clusters
//
// PRE: clusters have been created
// POST: the clusters have been split number_times_to_split
//
void CIRCUIT::generate_sub_clusters()
{
	debugSep;

	CLUSTERS::iterator cluster_iter;
	CLUSTER * cluster = 0;
	NUM_ELEMENTS number_times_to_split = 1;

	for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);

		cluster->split_cluster(number_times_to_split);
	}
}
// calculates the degree information for the circuit as a whole and
// calculates the degree information for the clusters themselves
//
// PRE: nothing
// POST: m_degree_info contains the fanin/fanout statistics
//       the clusters have calculated their degree_info
//
void CIRCUIT::calculate_degree_information()
{
	CLUSTERS::iterator cluster_iter;
	CLUSTER * cluster = 0;

	m_degree_info = new DEGREE_INFO;
	assert(m_degree_info);

	m_degree_info->calculate_degree_information_for_circuit(this);

	for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);

		cluster->calculate_degree_info();
	}
}


// creates the inter-cluster adjacentcy matrices
//
// PRE: the circuit has been created
//      clusters have been created
// POST: m_inter_cluster_connections contains the matrix that counts 
//       the number of connections to combinational nodes
//       m_inter_cluster_connections_for_dff contains the matrix that counts
//       the number of connections to flip-flops
void CIRCUIT::create_inter_cluster_matrix()
{
    CLUSTER_NUMBER_TYPE source_cluster_number, sink_cluster_number;
	EDGES::iterator edge_iter;
	EDGE * edge = 0;
	NODE * sink_node = 0;
	NUM_ELEMENTS nClusters = get_nClusters();
	m_inter_cluster_connections.resize(nClusters, nClusters);
	m_inter_cluster_connections_for_dff.resize(nClusters, nClusters);

    for (edge_iter = m_edges.begin(); edge_iter != m_edges.end(); edge_iter++)
    {
        edge = *edge_iter;
        assert(edge);
		sink_node = edge->get_sink_node();
		assert(sink_node);

		if (edge->is_inter_cluster() && ! edge->is_clock_edge())
		{
			sink_cluster_number = edge->get_sink_cluster_number();
			source_cluster_number = edge->get_source_cluster_number();
			assert(sink_cluster_number >= 0 && sink_cluster_number < nClusters);
			assert(source_cluster_number >= 0 && source_cluster_number < nClusters);

			if (sink_node->get_type() != NODE::SEQ)
			{
				m_inter_cluster_connections(source_cluster_number, sink_cluster_number) += 1;
			}
			else
			{
				m_inter_cluster_connections_for_dff(source_cluster_number, sink_cluster_number) += 1;
			}
		}
		else if (edge->is_intra_cluster() && ! edge->is_clock_edge() && sink_node->get_type() == NODE::SEQ)
		{
			sink_cluster_number = edge->get_sink_cluster_number();
			source_cluster_number = edge->get_source_cluster_number();
			assert(sink_cluster_number >= 0 && sink_cluster_number < nClusters);
			assert(source_cluster_number >= 0 && source_cluster_number < nClusters);

			m_inter_cluster_connections_for_dff(source_cluster_number, sink_cluster_number) += 1;
		}
    }
}



// Perform a final sanity check on the circuit
// 
// POST: a sanity check has been performed on all clusters
//
void CIRCUIT::final_sanity_check()
{

	// the sum of the edge lengths in a cluster should be equal to the number of inter cluster edges
	// the sum of the input shape should be equal to the sum of the output shape which should be equal 
	// to the number of inter cluster edges

	CLUSTERS::iterator cluster_iter;
	CLUSTER * cluster = 0;

	for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);

		cluster->final_sanity_check();
	}

}

//
// RETURNS: the dff in the circuit
//
NODES CIRCUIT::get_dffs()
{
	// must delete the graphs
	NODES::iterator node_iter;
	NODE * node = 0;
	NODES dffs;

	for (node_iter = m_nodes.begin(); node_iter != m_nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		if (node->get_type() == NODE::SEQ)
		{
			dffs.push_back(node);
		}
	}

	return dffs;
}

//
// RETURNS: the number of non-clock edges
//
NUM_ELEMENTS CIRCUIT::get_nEdges_without_clock_edges() const 	
{ 
	NUM_ELEMENTS nEdges = 0;

	if (m_global_clock)
	{
		assert(is_sequential());

		nEdges = m_number_edges - m_global_clock->get_fanout_degree();
	}
	else
	{
		nEdges = m_number_edges;
	}
	assert(nEdges >= 0);

	return nEdges; 
}

//
// RETURNS: the number of edges edges plus the edges that connect to the clock
//
NUM_ELEMENTS CIRCUIT::get_nEdges() const 
{
	return m_number_edges;
}

// gets the maximum number of nodes at a delay level 
// across all clusters
//
// PRE: clusters have been created 
// RETURNS: the maximum width across all clusters
NUM_ELEMENTS CIRCUIT::get_max_width() const
{
	CLUSTERS::const_iterator cluster_iter;
	CLUSTER * cluster = 0;
	SEQUENTIAL_LEVEL * seq_level = 0;
	NUM_ELEMENTS width = 0,
				 max_width = 0;

	for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);
		
		seq_level = cluster->get_sequential_level();
		assert(seq_level);

		width = seq_level->get_max_width();

		max_width = MAX(width, max_width);
	}
	debug("max width is " << max_width);


	return max_width;
}

// 
// POST: stats on the cluster have been calculated
//       and returned via:
//       size_dist, nPI_dist, nDFF_dist, nIntra_cluster_edges_dist,
//       nInter_cluster_edges_dist, wirelength_approx_dist
//
void CIRCUIT::get_cluster_stats
(
	DISTRIBUTION& size_dist, 
	DISTRIBUTION& nPI_dist,
	DISTRIBUTION& nDFF_dist, 
	DISTRIBUTION& nIntra_cluster_edges_dist, 
	DISTRIBUTION& nInter_cluster_edges_dist, 
	DISTRIBUTION& wirelength_approx_dist
)
{
	CLUSTERS::const_iterator cluster_iter;
	CLUSTER * cluster = 0;
	NUM_ELEMENTS size = 0, nPI = 0, nDFF = 0, 
				 nIntra_cluster_edges = 0, nInter_cluster_edges = 0, wirelength_approx = 0;

	for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);

		size = cluster->get_size();
		nPI = cluster->get_nPI(); 
		nDFF = cluster->get_nDFF();
		nIntra_cluster_edges = cluster->get_nIntra_cluster_edges();
		nInter_cluster_edges = cluster->get_nInter_cluster_edges(); 

		assert(Dlook_at);
		wirelength_approx = static_cast<NUM_ELEMENTS>(cluster->get_wirelength_approx());

		size_dist.push_back(size); 
		nPI_dist.push_back(nPI);
		nDFF_dist.push_back(nDFF);
		nIntra_cluster_edges_dist.push_back(nIntra_cluster_edges);
		nInter_cluster_edges_dist.push_back(nInter_cluster_edges);
		wirelength_approx_dist.push_back(wirelength_approx);
	}
}

//
//	RETURNS: the inter-cluster matrix of edges at edge_length
//
MATRIX CIRCUIT::get_inter_cluster_matrix(const LENGTH_TYPE& edge_length)
{
    CLUSTER_NUMBER_TYPE source_cluster_number, sink_cluster_number;
	EDGES::iterator edge_iter;
	EDGE * edge = 0;
	NODE * sink_node = 0;
	NUM_ELEMENTS nClusters = get_nClusters();
	MATRIX inter_cluster_connections(nClusters, nClusters);

    for (edge_iter = m_edges.begin(); edge_iter != m_edges.end(); edge_iter++)
    {
        edge = *edge_iter;
        assert(edge);
		sink_node = edge->get_sink_node();
		assert(sink_node);

		if (edge->is_inter_cluster() && ! edge->is_clock_edge() && edge->get_length() == edge_length)
		{
			sink_cluster_number = edge->get_sink_cluster_number();
			source_cluster_number = edge->get_source_cluster_number();
			assert(sink_cluster_number >= 0 && sink_cluster_number < nClusters);
			assert(source_cluster_number >= 0 && source_cluster_number < nClusters);

			inter_cluster_connections(source_cluster_number, sink_cluster_number) += 1;
		}
    }

	return inter_cluster_connections;
}

