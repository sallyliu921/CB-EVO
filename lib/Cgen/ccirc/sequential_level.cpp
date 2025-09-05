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





#include "sequential_level.h"
#include <algorithm>
#include <numeric>

const bool DSEQUENTIAL_LEVEL = false;

SEQUENTIAL_LEVEL::SEQUENTIAL_LEVEL()
{
	m_sequential_level_number	= -1;
	m_max_delay					= -1;
	m_number_of_nodes			= 0;
	m_is_clustered				= false;
}

SEQUENTIAL_LEVEL::SEQUENTIAL_LEVEL(const LEVEL_TYPE & seq_level_number)
{
	m_sequential_level_number	= seq_level_number;
	m_max_delay					= -1;
	m_number_of_nodes			= 0;
	m_is_clustered          	= false;
}

SEQUENTIAL_LEVEL::SEQUENTIAL_LEVEL(const SEQUENTIAL_LEVEL & another_sequential_level)
{
	m_delay_levels				= another_sequential_level.m_delay_levels;
	m_PI						= another_sequential_level.m_PI;
	m_PO						= another_sequential_level.m_PO;
	m_max_delay					= another_sequential_level.m_max_delay;
	m_sequential_level_number	= another_sequential_level.m_sequential_level_number;
	m_number_of_nodes			= another_sequential_level.m_number_of_nodes;
	m_is_clustered				= another_sequential_level.m_is_clustered;

    m_internal_edges			= another_sequential_level.m_internal_edges;			                       
    m_output_to_dff_edges		= another_sequential_level.m_output_to_dff_edges;
    m_input_from_dff_edges		= another_sequential_level.m_input_from_dff_edges;
    m_inter_cluster_input_edges	= another_sequential_level.m_inter_cluster_input_edges;
    m_inter_cluster_output_edges= another_sequential_level.m_inter_cluster_output_edges;
	m_inter_cluster_output_to_dff_edges = another_sequential_level.m_inter_cluster_output_to_dff_edges;
	m_inter_cluster_input_from_dff_edges = another_sequential_level.m_inter_cluster_input_from_dff_edges;
}   

SEQUENTIAL_LEVEL & SEQUENTIAL_LEVEL::operator=(const SEQUENTIAL_LEVEL & another_sequential_level)
{
	m_delay_levels				= another_sequential_level.m_delay_levels;
	m_PI						= another_sequential_level.m_PI;
	m_PO						= another_sequential_level.m_PO;
	m_max_delay					= another_sequential_level.m_max_delay;
	m_sequential_level_number	= another_sequential_level.m_sequential_level_number;
	m_number_of_nodes			= another_sequential_level.m_number_of_nodes;
	m_is_clustered				= another_sequential_level.m_is_clustered;

    m_internal_edges			= another_sequential_level.m_internal_edges;			                       
    m_output_to_dff_edges		= another_sequential_level.m_output_to_dff_edges;
    m_input_from_dff_edges		= another_sequential_level.m_input_from_dff_edges;
    m_inter_cluster_input_edges	= another_sequential_level.m_inter_cluster_input_edges;
    m_inter_cluster_output_edges= another_sequential_level.m_inter_cluster_output_edges;
	m_inter_cluster_output_to_dff_edges = another_sequential_level.m_inter_cluster_output_to_dff_edges;
	m_inter_cluster_input_from_dff_edges = another_sequential_level.m_inter_cluster_input_from_dff_edges;

	return (*this);
}

SEQUENTIAL_LEVEL::~SEQUENTIAL_LEVEL()
{
}


// adds the node to the sequential level
//
// PRE: the node is valid
// POST: the node has been added to its delay level
//       the edges of the node have been added
//       the node counter has been incremented
//       if the node has a primary output it has been added
//
void SEQUENTIAL_LEVEL::add_node
(
	NODE * node
)
{
	PORT * output_port = node->get_output_port();

	debugif(DSEQUENTIAL_LEVEL, "Adding node " << node->get_name() 
			<< " to seq level " << m_sequential_level_number);

	add_node_to_delay_level(node);
	add_edges_of_node(node);
	m_number_of_nodes++;

	if (output_port->get_type() == PORT::EXTERNAL) m_PO.push_back(output_port);
	
}

// adds a primary input
//
// PRE: the port is valid and is a primary input
// POST: if the port is not the clock it is added m_PI
//       otherwise it is assigned to m_clock_port
//
void  SEQUENTIAL_LEVEL::add_primary_input
(
	PORT * port
)
{
	assert(port && port->get_type() == PORT::EXTERNAL);

	if (port->get_io_direction() != PORT::CLOCK)
	{
		add_output_edges(port);
		m_PI.push_back(port);
	}
	else
	{
		assert(m_clock_port == 0);
		m_clock_port = port;
	}

}

// adds all the edges of the node 
//
// PRE: the node is valid
// POST: the input and output edges have been added as either 
//       intra- or inter-cluster input, or inter-cluster output edges
//
void SEQUENTIAL_LEVEL::add_edges_of_node
(
	NODE * node
)
{
	assert(node);
	PORT * output_port = node->get_output_port();

	add_input_edges(node);
	add_output_edges(output_port);
}

// adds all the input edges of the node
//
// PRE: the node is valid
// POST: all intra-cluster edges of the node have been added to m_internal_edges
//       all inter-cluster input edges of the node have been added to m_inter_cluster_input_edges
//
void SEQUENTIAL_LEVEL::add_input_edges
(
	NODE * node
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

		if (m_is_clustered && edge->is_inter_cluster())
		{
			add_input_edge(edge, m_inter_cluster_input_edges, m_inter_cluster_input_from_dff_edges);
		}
		else
		{
			add_input_edge(edge, m_internal_edges, m_input_from_dff_edges);
		}
	}
}
// adds the edge connected to the node to a list of edges
//
// PRE: the edge is valid
// POST: if the edge is not a clock edge it has been assigned 
//       to either comb_graph_edges or inputs_from_dff_edges
//
void SEQUENTIAL_LEVEL::add_input_edge
(
	EDGE * edge,
	EDGES & comb_graph_edges,
	EDGES & inputs_from_dff_edges
)
{
	NODE * source_node =  edge->get_source_node();
	NODE * sink_node = edge->get_sink_node();
	PORT * input_port = edge->get_sink();
	assert(sink_node && input_port);
	

	if (sink_node->get_type() == NODE::SEQ && input_port->get_io_direction() != PORT::CLOCK)
	{
		if (m_is_clustered && edge->is_inter_cluster())
		{
			debugif(DCLUSTER, "We have an inter cluster edge (" << edge->get_source_cluster_number() << ","
					<< edge->get_sink_cluster_number() << ")" 
					<< " ( " << edge->get_name() << " ) that outputs to a dff");
		}
	}

	if (source_node)
	{	
		if (source_node->get_type() == NODE::SEQ)
		{
			inputs_from_dff_edges.push_back(edge);
		}

		if (sink_node->get_type() == NODE::COMB)
		{
			comb_graph_edges.push_back(edge);
		}
	}
	else 
	{
		add_primary_input_edge(edge, comb_graph_edges);
	}
}

//
// if the edge is not connected to a clock or a flip-flop add it to the list of 
// intra-cluster combinational edges
//
// PRE: edge is valid
// POST: if the edge is not a clock edge it is added to 
//       comb_graph_edges
void SEQUENTIAL_LEVEL::add_primary_input_edge
(
	EDGE * edge,
	EDGES & comb_graph_edges
)
{
	assert(edge);
	NODE * sink_node = edge->get_sink_node();
	PORT * pi_port = edge->get_source();
	assert(sink_node && pi_port);

	if (pi_port->get_io_direction() != PORT::CLOCK && sink_node->get_type() == NODE::COMB)
	{
		comb_graph_edges.push_back(edge);
	}
}
					

// adds the edges of the output port to a list of edges
//
// PRE: the output port is valid
// POST: all inter-cluster combinational edges from the port have been added to m_inter_cluster_output_edges
//       all inter-cluster edges to dff from the port have been added to m_inter_cluster_output_to_dff_edges
//       all intra-cluster edges to dff from the port have been added to m_output_to_dff_edges
//
void SEQUENTIAL_LEVEL::add_output_edges
(
	PORT * output_port
)
{
	assert(output_port);

	EDGES::iterator edge_iter;
	EDGE * edge = 0;
	EDGES output_edges = output_port->get_edges();
	NODE * sink_node = 0;

	for (edge_iter = output_edges.begin(); edge_iter != output_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);
		sink_node = edge->get_sink_node();
		assert(sink_node);


		if (m_is_clustered && edge->is_inter_cluster())
		{
			debugif(DCLUSTER, "Edge " << edge->get_name() << " is inter cluster ");
			debugif(DCLUSTER, "is is between cluster " << edge->get_source_cluster_number() << " and "
					<< edge->get_sink_cluster_number());

			if (sink_node->get_type() == NODE::COMB)
			{
				m_inter_cluster_output_edges.push_back(edge);
			}
			else
			{
				add_output_to_dff_edge(edge, m_inter_cluster_output_to_dff_edges);
			}
		}
		else
		{
			// if the edge is internal it should be found by add_input_edges
			add_output_to_dff_edge(edge, m_output_to_dff_edges);
		}
	}
}


// adds the edge if the edge is an intra-cluster edge that connects to a flip-flop
//
// PRE: the edge is valid
// POST: if the edge was intra-cluster connecting to a flip flop it was added to output_to_dff_edges
//
void SEQUENTIAL_LEVEL::add_output_to_dff_edge
(
	EDGE * edge,
	EDGES & output_to_dff_edges
)
{
	assert(edge);

	NODE * sink_node = edge->get_sink_node();
	assert(sink_node);

	// internal edges are added by the add_input_edge
	if (sink_node->get_type() == NODE::SEQ)
	{
		output_to_dff_edges.push_back(edge);
	}
}

// adds the node to its delay level
// 
// PRE: the node is valid
//      the node has its delay level defined
// POST: the node has been added to its delay level
//
void SEQUENTIAL_LEVEL::add_node_to_delay_level
(
	NODE * node
)
{

	DELAY_TYPE max_comb_delay = node->get_max_comb_delay_level();
	assert(max_comb_delay >= 0 && max_comb_delay <= m_max_delay);
	assert(max_comb_delay < static_cast<signed>(m_delay_levels.size()));

	m_delay_levels[max_comb_delay].push_back(node);
}

// sets the maximum combinational delay in the circuit
// 
// PRE: max_delay contains the delay 
// POST: m_max_delay has been set
//       m_delay_levels has been resized 
//
void SEQUENTIAL_LEVEL::set_max_combinational_delay
(
	const DELAY_TYPE & max_delay
)
{ 
	assert(max_delay >= 0);
	m_max_delay = max_delay;
	m_delay_levels.resize(m_max_delay+1);

}

//
// RETURNS: the maximum combinational delay
//
DELAY_TYPE SEQUENTIAL_LEVEL::get_maximum_combinational_delay() const
{
 	assert(m_delay_levels.size() - 1 == static_cast<unsigned>(m_max_delay));
 
 	return m_max_delay;
}

//
// gets the number of nodes at each delay level
//
// PRE: the delay levels have been created
// RETURNS: the node shape
//
SHAPE SEQUENTIAL_LEVEL::get_node_shape() const
{
	SHAPE node_shape;

	DELAY_TYPE delay_level_index = 0;
	DELAY_TYPE max_delay_level_index = static_cast<signed>(m_delay_levels.size());

	for (delay_level_index = 0; delay_level_index < max_delay_level_index; delay_level_index++)
	{
		node_shape.push_back(m_delay_levels[delay_level_index].size());
	}

	if (m_sequential_level_number == 0 && max_delay_level_index > 0)
	{
		// we need to include the PI
		node_shape[0] += m_PI.size();
	}


	return node_shape;
}

//
// the number of inputs at each delay level
//
// PRE: the delay levels have been created
// RETURNS: the input shape
//
SHAPE SEQUENTIAL_LEVEL::get_input_shape() const
{

	SHAPE input_shape;

	DELAY_TYPE delay_level_index = 0;
	DELAY_TYPE max_delay_level_index = static_cast<signed>(m_delay_levels.size());
	NUM_ELEMENTS number_of_non_clock_inputs = 0;


	//debug("input_shape");
	for (delay_level_index = 0; delay_level_index < max_delay_level_index; delay_level_index++)
	{
		number_of_non_clock_inputs = get_number_of_non_clock_inputs(m_delay_levels[delay_level_index]);

		//debug("delay: " << delay_level_index << " nInputs: " << number_of_non_clock_inputs);
		input_shape.push_back(number_of_non_clock_inputs);
	}

	return input_shape;
}


// gets the number of non-clock inputs connected to the nodes in 
// nodes_at_a_delay_level
//
//
// PRE: nothing
// RETURNS: the number of non-clock inputs connected to the nodes in 
//          nodes_at_a_delay_level
//
NUM_ELEMENTS SEQUENTIAL_LEVEL::get_number_of_non_clock_inputs
(
	const NODES & nodes_at_a_delay_level
) const
{
	NUM_ELEMENTS number_of_non_clock_inputs = 0;
	NODES::const_iterator node_iter;
	NODE * node = 0;

	for (node_iter = nodes_at_a_delay_level.begin(); node_iter != nodes_at_a_delay_level.end(); 
			node_iter++)
	{
		node = *node_iter;
		assert(node);

		if (node->get_type() == NODE::SEQ)
		{
			// dff for our purposes have no inputs
			assert(Dlook_at);
			//number_of_non_clock_inputs += 1;
		}
		else
		{
			assert(node->get_type() == NODE::COMB);

			number_of_non_clock_inputs += node->get_fanin_degree();

		}
	}


	return number_of_non_clock_inputs;
}

// the number of outputs at each delay level
//
// PRE: the delay levels have been created
// RETURNS: the output shape
//
SHAPE SEQUENTIAL_LEVEL::get_output_shape() const
{


	DELAY_TYPE delay_level_index = 0;
	DELAY_TYPE max_delay_level_index = static_cast<signed>(m_delay_levels.size());
	SHAPE output_shape(max_delay_level_index, 0);
	NUM_ELEMENTS number_of_outputs = 0;

	//debug("output_shape");
	// we do not have an output shape for the last delay level
	for (delay_level_index = 0; delay_level_index < max_delay_level_index-1; delay_level_index++)
	{
		number_of_outputs = get_number_of_outputs(m_delay_levels[delay_level_index]);
		//debug("delay: " << delay_level_index << " nOutputs: " << number_of_outputs);
		output_shape[delay_level_index] = number_of_outputs;
	}

	if (m_sequential_level_number == 0 && max_delay_level_index > 0)
	{
		// we need to include the PI
		//
		//debug("pi: " << get_number_of_outputs(m_PI));
		output_shape[0] += get_number_of_outputs(m_PI);
	}


	return output_shape;
}

//
// gets the number of non-clock inputs connected to the ports in 
// ports_at_a_delay_level
//
// PRE: nothing
// RETURNS: the number of outputs connected to the ports in 
//          ports_at_a_delay_level
//
NUM_ELEMENTS SEQUENTIAL_LEVEL::get_number_of_outputs
(
	const PORTS & ports_at_a_delay_level
) const
{
	NUM_ELEMENTS number_of_outputs = 0;
	PORTS::const_iterator port_iter;
	PORT * port = 0;

	for (port_iter = ports_at_a_delay_level.begin(); port_iter != ports_at_a_delay_level.end(); 
			port_iter++)
	{
		port = *port_iter;
		assert(port);

		number_of_outputs += port->get_fanout_degree_to_combinational_nodes();
	}

	return number_of_outputs;
}

// the number of outputs at a delay level
//
// PRE: the delay levels have been created
// RETURNS: the number of output edges at a delay level
//
NUM_ELEMENTS SEQUENTIAL_LEVEL::get_number_of_outputs
(
	const NODES & nodes_at_a_delay_level
) const
{
	NUM_ELEMENTS number_of_outputs = 0;
	NODES::const_iterator node_iter;
	NODE * node = 0;

	for (node_iter = nodes_at_a_delay_level.begin(); node_iter != nodes_at_a_delay_level.end(); 
			node_iter++)
	{
		node = *node_iter;
		assert(node);

		number_of_outputs += node->get_fanout_degree_to_combinational_nodes();
	}

	return number_of_outputs;
}

//
// gets the number of primary outputs at each delay level
//
// PRE: the delay levels have been created
// RETURNS: the PO shape
//
SHAPE SEQUENTIAL_LEVEL::get_PO_shape() const
{
	SHAPE shape;

	if (m_delay_levels.empty())
	{
		return shape;
	}
	shape.resize(m_delay_levels.size(), 0);

	PORTS::const_iterator port_iter;
	PORT * port = 0;
	NODE * node = 0;
	DELAY_TYPE delay_level = 0;

	for (port_iter = m_PO.begin(); port_iter != m_PO.end(); port_iter++)
	{
		port = *port_iter;
		node = port->get_my_node();
		assert(node);

		delay_level = node->get_max_comb_delay_level();
		assert(delay_level >= 0 && delay_level < static_cast<signed>(m_delay_levels.size()));

		shape[delay_level]++;
	}

	return shape;
}


//
// gets the number of latched nodes at each delay level
//
// PRE: the delay levels have been created
// RETURNS: the latched shape
//
SHAPE SEQUENTIAL_LEVEL::get_latched_shape() const
{
	SHAPE latched_shape;

	if (m_delay_levels.empty())
	{
		return latched_shape;
	}
	latched_shape.resize(m_delay_levels.size(), 0);

	if 	(m_output_to_dff_edges.size() == 0 && m_inter_cluster_output_to_dff_edges.size() == 0)
	{
		return latched_shape;
	}
	add_latched_nodes_to_shape(latched_shape, m_output_to_dff_edges);

	if (m_is_clustered)
	{
		add_latched_nodes_to_shape(latched_shape, m_inter_cluster_output_to_dff_edges);
	}

	return latched_shape;
}

//
// for the edges in output_to_dff_edges adds the location of the source
// (the node that connects to the input of the DFF)
// to the latched shape
//
// PRE: the delay levels have been created
//      output_to_dff_edges contains edges 
//      that connect to the inputs of dff
// RETURNS: latched_shape
//
void SEQUENTIAL_LEVEL::add_latched_nodes_to_shape
(
	SHAPE & latched_shape,
	const EDGES & output_to_dff_edges
) const
{
	EDGES::const_iterator edge_iter;
	EDGE * edge = 0;
	NODE * source_node = 0;
	DELAY_TYPE delay_level = 0;


	if (output_to_dff_edges.empty())
	{
		return;
	}

	for (edge_iter = output_to_dff_edges.begin(); edge_iter != output_to_dff_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		source_node = edge->get_source_node();
		if (source_node)
		{
			delay_level = source_node->get_max_comb_delay_level();
			latched_shape[delay_level]++;
			assert(delay_level >= 0 && delay_level < static_cast<signed>(m_delay_levels.size()));
		}
		else
		{
			latched_shape[0]++;
		}
	}
}


//
// RETURNS: the number of edges at each edge length that are intra-cluster edges
//
DISTRIBUTION SEQUENTIAL_LEVEL::get_intra_cluster_edge_length_distribution() const
{
	//debug("intra cluster edge length");
	return find_edge_length_distribution(m_internal_edges);
}

//
// RETURNS: the number of edges at each edge length that are inter-cluster input edges
//
DISTRIBUTION SEQUENTIAL_LEVEL::get_inter_cluster_input_edge_length_distribution() const
{
	//debug("inter cluster input edge length");
	return find_edge_length_distribution(m_inter_cluster_input_edges);
}

//
// RETURNS: the number of edges at each edge length that are inter-cluster output edges
//
DISTRIBUTION SEQUENTIAL_LEVEL::get_inter_cluster_output_edge_length_distribution() const
{
	//debug("inter cluster output edge length");
	return find_edge_length_distribution(m_inter_cluster_output_edges);
}

// 
// PRE: the delay levels have been created
// RETURNS: the edge length distribution for the edges in 
//          edges
DISTRIBUTION SEQUENTIAL_LEVEL::find_edge_length_distribution
(
	const EDGES & edges
) const
{
	DISTRIBUTION edge_length_distribution(m_delay_levels.size(), 0);

	EDGES::const_iterator edge_iter;
	EDGE * edge = 0;
	DELAY_TYPE edge_length = 0;

	for (edge_iter = edges.begin(); edge_iter != edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		edge_length = edge->get_length();
		edge_length_distribution[edge_length]++;
		assert(edge_length >= 0 && static_cast<unsigned>(edge_length) < edge_length_distribution.size());
	}

	return edge_length_distribution;
}


//
// RETURNS: the number of edges at each edge length that are intra-cluster edges
//          that input into the delay level specified
//
DISTRIBUTION SEQUENTIAL_LEVEL::get_intra_cluster_input_edge_length_distribution(const DELAY_TYPE& delay_level)
{
	// inefficient but hey it works.
	EDGES input_edges_to_a_delay_level;
	// because we are not using 0 length edges we can leave it empty
	EDGES dff_edges; assert(Dlook_at);	
	EDGES::iterator edge_iter;
	EDGE * edge = 0;
	NODE * sink_node = 0;

	for (edge_iter = m_internal_edges.begin(); edge_iter != m_internal_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);
		sink_node = edge->get_sink_node();
		assert(sink_node);

		if (sink_node->get_max_comb_delay_level() == delay_level)
		{
			input_edges_to_a_delay_level.push_back(edge);
		}
	}

	return find_edge_length_distribution(input_edges_to_a_delay_level);
}

//
// RETURNS: the number of edges at each edge length that are intra-cluster edges
//          that ouput out of the delay level specified
//
DISTRIBUTION SEQUENTIAL_LEVEL::get_intra_cluster_output_edge_length_distribution(const DELAY_TYPE& delay_level)
{ 
	// inefficient but hey it works.
	EDGES output_edges_to_a_delay_level;
	// because we are not using 0 length edges we can leave it empty
	EDGES dff_edges; assert(Dlook_at);	
	EDGES::iterator edge_iter;
	EDGE * edge = 0;
	NODE * source_node = 0;

	for (edge_iter = m_internal_edges.begin(); edge_iter != m_internal_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);
		source_node = edge->get_source_node();

		if (source_node)
		{
			if (source_node->get_max_comb_delay_level() == delay_level)
			{
				output_edges_to_a_delay_level.push_back(edge);
			}
		}
		else
		{
			// we have a primary input if we are looking for delay level 0 nodes 
			// add the edge the output edges

			if (delay_level == 0)
			{
				output_edges_to_a_delay_level.push_back(edge);
			}
		}

	}

	return find_edge_length_distribution(output_edges_to_a_delay_level);
}
//
// RETURNS: the number of edges at each edge length that are inter-cluster edges
//          that input into the delay level specified
//
DISTRIBUTION SEQUENTIAL_LEVEL::get_inter_cluster_input_edge_length_distribution(const DELAY_TYPE& delay_level)
{ 
	// inefficient but hey it works.
	EDGES input_edges_to_a_delay_level;
	// because we are not using 0 length edges we can leave it empty
	EDGES dff_edges; assert(Dlook_at);	
	EDGES::iterator edge_iter;
	EDGE * edge = 0;
	NODE * sink_node = 0;

	for (edge_iter = m_inter_cluster_input_edges.begin(); edge_iter != m_inter_cluster_input_edges.end();
			edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);
		sink_node = edge->get_sink_node();
		assert(sink_node);

		if (sink_node->get_max_comb_delay_level() == delay_level)
		{
			input_edges_to_a_delay_level.push_back(edge);
		}
	}

	return find_edge_length_distribution(input_edges_to_a_delay_level);
}
//
// RETURNS: the number of edges at each edge length that are inter-cluster edges
//          that ouput out of the delay level specified
//
DISTRIBUTION SEQUENTIAL_LEVEL::get_inter_cluster_output_edge_length_distribution(const DELAY_TYPE& delay_level)

{
	// inefficient but hey it works.
	EDGES output_edges_to_a_delay_level;
	// because we are not using 0 length edges we can leave it empty
	EDGES dff_edges; assert(Dlook_at);	
	EDGES::iterator edge_iter;
	EDGE * edge = 0;
	NODE * source_node = 0;

	for (edge_iter = m_inter_cluster_output_edges.begin(); edge_iter != m_inter_cluster_output_edges.end();
			edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);
		source_node = edge->get_source_node();

		if (source_node)
		{
			if (source_node->get_max_comb_delay_level() == delay_level)
			{
				output_edges_to_a_delay_level.push_back(edge);
			}
		}
		else
		{
			// we have a primary input if we are looking for delay level 0 nodes 
			// add the edge the output edges

			if (delay_level == 0)
			{
				output_edges_to_a_delay_level.push_back(edge);
			}
		}

	}

	return find_edge_length_distribution(output_edges_to_a_delay_level);
}
//
// gets the fanout distribution, which is the number of nodes at each fanout 
// starting at 0
//
// PRE: max_fanout is the maximum fanout for the circuit
// RETURNS: fanout distribution
//
DISTRIBUTION SEQUENTIAL_LEVEL::get_fanout_distribution
(
	const NUM_ELEMENTS & max_fanout
) const
{

	DISTRIBUTION fanout_distribution(max_fanout+1, 0);
	DELAY_TYPE delay_level_index = 0;

	assert(m_max_delay + 1 == static_cast<signed>(m_delay_levels.size()));

	for (delay_level_index = 0; delay_level_index <= m_max_delay; delay_level_index++)
	{
		add_fanout_distribution_for_delay_level(fanout_distribution, m_delay_levels[delay_level_index]);
	}

	if (m_sequential_level_number == 0)
	{
		// we need to include the fanout of the PI
		add_fanout_distribution_for_primary_inputs(fanout_distribution);
	}


	return fanout_distribution;
}
//
// adds to fanout_distribution the 
// fanout distribution at the specified delay level
//
// PRE: fanout_distribution's size can handle the largest fanout 
// RETURNS: the delay level's fanouts have been added to
//          fanout_distribution
//
void SEQUENTIAL_LEVEL::add_fanout_distribution_for_delay_level
(
	DISTRIBUTION & fanout_distribution,
	const DELAY_LEVEL & delay_level
) const
{
	DELAY_LEVEL::const_iterator node_iter;
	NODE * node = 0;
	NUM_ELEMENTS fanout_number = 0;

	for (node_iter = delay_level.begin(); node_iter != delay_level.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		assert(Dlook_at);
		fanout_number = node->get_fanout_degree_to_combinational_nodes();

		assert(fanout_number >= 0 && fanout_number < static_cast<signed>(fanout_distribution.size()));

		fanout_distribution[fanout_number]++;
	}
}


//
// adds to fanout_distribution of the primary inputs
// to the fanout distribution
//
// PRE: fanout_distribution's size can handle the largest fanout 
// RETURNS: the primary inputs fanouts have been added to
//          fanout_distribution
//
void SEQUENTIAL_LEVEL::add_fanout_distribution_for_primary_inputs
(
	DISTRIBUTION & fanout_distribution
) const
{
	PORT * port = 0;
	PORTS::const_iterator port_iter;

	NUM_ELEMENTS fanout_number = 0;

	for (port_iter = m_PI.begin(); port_iter != m_PI.end(); port_iter++)
	{
		port = *port_iter;
		assert(port);
		assert(port->get_io_direction() != PORT::CLOCK);

		fanout_number = port->get_fanout_degree_to_combinational_nodes();

		assert(fanout_number >= 0 && fanout_number < static_cast<signed>(fanout_distribution.size()));

		fanout_distribution[fanout_number]++;
	}
}

//
// gets the maximum width of the circuit which is the 
// size of the largest delay level
//
// PRE: delay levels have been created
// RETURNS: maximum width 
//
NUM_ELEMENTS SEQUENTIAL_LEVEL::get_max_width() const
{
	SHAPE::const_iterator shape_iter;

	SHAPE node_shape = get_node_shape();

	shape_iter = max_element(node_shape.begin(), node_shape.end());

	NUM_ELEMENTS max_width = *shape_iter;

	return max_width;
}
//
// peforms a final sanity check
//
// PRE: nothing
// POST: the number of edges in the graph has been compared to 
//       the number of inputs and outputs in the graph
//
void SEQUENTIAL_LEVEL::final_sanity_check()
{
	
	// make sure that the sum of the input and output shape is equal to the sum of 
	//
	// the inter cluster input and output edges and 2 * the intra cluster edges
	//
	// and make sure that is equal to the sum of the fanout distribution
	SHAPE input_shape = get_input_shape();
	SHAPE output_shape = get_output_shape();

	DISTRIBUTION edge_lengths = get_intra_cluster_edge_length_distribution();
    DISTRIBUTION inter_cluster_input_edge_lengths = get_inter_cluster_input_edge_length_distribution();
    DISTRIBUTION inter_cluster_output_edge_lengths = get_inter_cluster_output_edge_length_distribution();

	//DISTRIBUTION fanout_dist = get_fanout_distribution();

	NUM_ELEMENTS shape_sum = 0;
	NUM_ELEMENTS edge_sum = 0;

	shape_sum = accumulate(input_shape.begin(), input_shape.end(), shape_sum);
	shape_sum = accumulate(output_shape.begin(), output_shape.end(), shape_sum);

	edge_sum = 2 * accumulate(edge_lengths.begin(), edge_lengths.end(), edge_sum);
	edge_sum = accumulate(inter_cluster_input_edge_lengths.begin(), inter_cluster_input_edge_lengths.end(), edge_sum);
	edge_sum = accumulate(inter_cluster_output_edge_lengths.begin(), inter_cluster_output_edge_lengths.end(), edge_sum);
	
	//debug("edge sum / shape sum " << edge_sum << " / " << shape_sum);
	assert(edge_sum == shape_sum);


	//if the sequential level is unclustered each delay level should have 
	//some nodes in the shape
}
