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





#include "edges_and_nodes.h"
#include <algorithm>
#include <limits.h>
#include "circ.h"

const string EDGE_CONNECTION_TEXT = "_TO_";	
const string EDGE_SUFFIX			= "_EDGE";	

EDGE::EDGE()
{
	m_source 		= 0;
	m_destination 	= 0;
	m_length 		= 0;
}

EDGE::EDGE (const EDGE & another_edge)
{
	m_source 		= another_edge.m_source;
	m_destination 	= another_edge.m_destination;
	m_length 		= another_edge.m_length;
	m_name			= another_edge.m_name;
}

EDGE & EDGE::operator=(const EDGE & another_edge)
{
	m_source 		= another_edge.m_source;
	m_destination 	= another_edge.m_destination;
	m_length 		= another_edge.m_length;
	m_name			= another_edge.m_name;

	return (*this);
}
EDGE::EDGE
(
	PORT * source_port, 
	PORT * sink_port, 
	const string & edge_name, 
	const LENGTH_TYPE& length
)
{
	m_source 		= source_port;
	m_destination 	= sink_port;
	m_name 			= edge_name;
	m_length		= length;
}

EDGE::EDGE
(
	PORT * source_port, 
	PORT * sink_port, 
	const LENGTH_TYPE&  length
)
{
	NODE * sink_node;

	m_source 		= source_port;
	m_destination 	= sink_port;
	m_length		= length;

	assert(m_destination);
	sink_node = m_destination->get_my_node();
	assert(sink_node);
	
	m_name 	= get_edge_name(sink_node->get_name(), source_port->get_name());
}

EDGE::~EDGE()
{
	// do not delete the ports
	// the nodes will delete their own ports
	m_source  = 0;
	m_destination = 0;
}

NODE * EDGE::get_sink_node() const
{
	assert(m_destination);
	return (m_destination->get_my_node());
}

NODE * EDGE::get_source_node() const
{
	assert(m_source);
	return (m_source->get_my_node());
}

//
// RETURNS: a name for the edge
//
string EDGE::get_edge_name
(
	const string & input_node_name,
	const string & output_port_name
) const
{
	assert(! input_node_name.empty());
	assert(! output_port_name.empty());

	return output_port_name + EDGE_CONNECTION_TEXT + input_node_name + EDGE_SUFFIX;
}


bool EDGE::is_inter_cluster() const
{
	NODE * sink_node = get_sink_node();
	NODE * source_node = get_source_node();
	CLUSTER_NUMBER_TYPE source_cluster_number, sink_cluster_number;
	bool is_inter_cluster_edge = false;

	assert(sink_node);
	sink_cluster_number = sink_node->get_cluster_number();

	if (source_node)
	{
		source_cluster_number = source_node->get_cluster_number();
	}
	else
	{
		// we have a primary input
		assert(m_source);
		source_cluster_number = m_source->get_cluster_number();
	}

	is_inter_cluster_edge = (source_cluster_number != sink_cluster_number);

	return is_inter_cluster_edge;
}

bool EDGE::is_intra_cluster() const
{
	return (! is_inter_cluster());
}


//
// PRE: sub-clusters exist
// RETURNS: whether the edge is internal to the cluster and all sub-clusters
//
bool EDGE::is_intra_cluster_and_intra_sub_cluster() const
{
	bool same_cluster_and_sub_cluster = 
		is_intra_cluster() && ! is_inter_sub_cluster();

	return same_cluster_and_sub_cluster;
}

//
// Does the edge cross sub-clusters?
//
// PRE: sub-clusters exist
// RETURNS: whether the edges crosses sub-clusters
//
bool EDGE::is_inter_sub_cluster() const
{
	NODE * sink_node = get_sink_node();
	assert(sink_node && m_source);
	bool inter_sub_cluster = false;

	CLUSTER_NUMBERS source_cluster_numbers = get_source_sub_cluster_numbers();
	CLUSTER_NUMBERS sink_cluster_numbers = sink_node->get_sub_cluster_numbers();

	// if the size of the sub cluster numbers is not the same they can't
	// possibly be in the same cluster
	inter_sub_cluster = source_cluster_numbers.size() != sink_cluster_numbers.size();


	// check the values
	if (! inter_sub_cluster)
	{
		inter_sub_cluster = (source_cluster_numbers != sink_cluster_numbers);
	}

	return inter_sub_cluster;
}

//
// Is the destination of the edge a flip-flop?
//
bool EDGE::is_sink_a_dff() const
{
	NODE * sink_node = get_sink_node();
	assert(sink_node);

	return (sink_node->get_type() == NODE::SEQ);
}

//
// Is the edge an edge that connects to a clock port of a node?
//
bool EDGE::is_clock_edge() const
{
	bool clock_edge = false;

	clock_edge = (m_source->get_io_direction() == PORT::CLOCK);

	return clock_edge;
}

//
// Get the cluster number of the sink node
//
CLUSTER_NUMBER_TYPE EDGE::get_sink_cluster_number() const
{
	const NODE * sink_node = get_sink_node();
	assert(sink_node);

	return (sink_node->get_cluster_number());
}

//
// Get the cluster number of the source node
//
CLUSTER_NUMBER_TYPE EDGE::get_source_cluster_number() const
{
	const NODE * source_node = get_source_node();

	if (source_node)
	{
		return (source_node->get_cluster_number());
	}
	else
	{
		assert(m_source->get_type() == PORT::EXTERNAL);

		return (m_source->get_cluster_number());
	}
}
//
// Get the sub-cluster numbers of the sink node
//
CLUSTER_NUMBERS EDGE::get_sink_sub_cluster_numbers() const
{
	const NODE * sink_node = get_sink_node();
	assert(sink_node);

	return (sink_node->get_sub_cluster_numbers());
}
//
// Get the sub-cluster numbers of the source node
//
CLUSTER_NUMBERS EDGE::get_source_sub_cluster_numbers() const
{
	NODE * source_node = get_source_node();

	if (source_node)
	{
		return (source_node->get_sub_cluster_numbers());
	}
	else
	{
		assert(m_source->get_type() == PORT::EXTERNAL);

		return (m_source->get_sub_cluster_numbers());
	}
}

// returns the inter-cluster distance of an edge
//
// ie. for a circuit that has been divided 3 times
//     we look at the sub-clusters of the source and sink 
//     clusters at a depth of 3
//
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
// RETURNS: the inter-cluster distance of the edge
//
LENGTH_TYPE EDGE::get_inter_cluster_distance() const
{
	CLUSTER_NUMBER_TYPE source_cluster_number = get_source_cluster_number();
	CLUSTER_NUMBER_TYPE sink_cluster_number = get_sink_cluster_number();
	CLUSTER_NUMBERS source_sub_cluster_numbers = get_source_sub_cluster_numbers();
	CLUSTER_NUMBERS sink_sub_cluster_numbers = get_sink_sub_cluster_numbers();
	assert(source_sub_cluster_numbers.size() == sink_sub_cluster_numbers.size());
	LENGTH_TYPE index = 0;

	// sub_cluster_size is the partitioning depth of our clusters
	LENGTH_TYPE sub_cluster_size = static_cast<LENGTH_TYPE>(source_sub_cluster_numbers.size());
	LENGTH_TYPE distance = sub_cluster_size + 1;

	if (source_cluster_number == sink_cluster_number)
	{
		// we are in the same cluster
		distance--;

		// check the sub cluster numbers 
		while (index < sub_cluster_size && source_sub_cluster_numbers[index] == sink_sub_cluster_numbers[index])
		{
			index++;
		}

		if (index != sub_cluster_size)
		{
			distance = distance - index;
		}
		else
		{
			distance = 0;
		}
	}
	assert(distance >= 0);

	return distance;
}


NODE::NODE()
{
	m_type			=	NODE::COMB;
	m_name			= 	"";
	m_output_port	=	0;

	m_delay_level	=	-1;
	m_colour_mark	=	NODE::UNMARKED;

	m_lut			=	0;

	m_horizontal_position = 0;
}


NODE::NODE(const string & node_name)
{
	m_name			= 	node_name;
	m_type			=	NODE::COMB;
	m_output_port	=	0;
	m_delay_level	=	-1;
	m_lut			=	0;
	m_colour_mark	=	NODE::UNMARKED;
	
	m_horizontal_position = 0;
}

NODE::NODE(const string & node_name, const NODE::NODE_TYPE & node_type)
{
	m_name			= 	node_name;
	m_type			=	node_type;
	m_output_port	=	0;
	m_delay_level	=	-1;
	m_lut			=	0;
	m_colour_mark	=	NODE::UNMARKED;

	m_horizontal_position = 0;
}

NODE::NODE(const NODE  & another_node)
{
	m_type				= another_node.m_type;
	m_name				= another_node.m_name;
	m_input_ports		= another_node.m_input_ports;
	m_output_port		= another_node.m_output_port;
	m_delay_level		= another_node.m_delay_level;
	m_colour_mark		= another_node.m_colour_mark;
	m_lut				= another_node.m_lut;
	m_cluster_number	= another_node.m_cluster_number;
	m_sub_cluster_numbers	= another_node.m_sub_cluster_numbers;
	m_horizontal_position = another_node.m_horizontal_position;
}

NODE& NODE::operator=(const NODE  & another_node)
{
	m_type				= another_node.m_type;
	m_name				= another_node.m_name;
	m_input_ports		= another_node.m_input_ports;
	m_output_port		= another_node.m_output_port;
	m_delay_level		= another_node.m_delay_level;
	m_colour_mark		= another_node.m_colour_mark;
	m_lut				= another_node.m_lut;
	m_cluster_number	= another_node.m_cluster_number;
	m_sub_cluster_numbers	= another_node.m_sub_cluster_numbers;
	m_horizontal_position = another_node.m_horizontal_position;

	return (*this);
}

NODE::~NODE()
{
	// delete all the input ports
	PORTS::iterator port_iter;
	for (port_iter = m_input_ports.begin(); port_iter != m_input_ports.end(); port_iter++)
	{
		// if we have an input port we should be the one to delete it
		assert(*port_iter);
		delete (*port_iter);
	}

	m_input_ports.clear();

	assert(m_output_port);
	delete m_output_port;
	m_output_port = 0;

	// if we have a lut delete it
	if (m_lut)
	{
		delete m_lut;
		m_lut = 0;
	}
}

// 
// Create an add a port to the node
//
// PRE: nothing
// POST: an input or output port has been created and added to the node
PORT * NODE::create_and_add_port
(
	const string & port_name,
	const PORT::PORT_TYPE & port_type,
	const PORT::IO_DIRECTION & io_direction,
	const PORT::EXTERNAL_TYPE & external_type
)
{
	PORT * port = new PORT(port_name, port_type, io_direction, external_type, this);
	assert(port);

	if (io_direction == PORT::INPUT || io_direction == PORT::CLOCK)
	{
		m_input_ports.push_back(port);
	}
	else 
	{
		// there should not be a previous output port
		assert(m_output_port == 0);

		m_output_port = port;
	}

	port->set_my_node(this);

	return port;
}


// 
// adds a port to the node
//
// PRE: nothing
// POST: an input or output port has been added to the node
void NODE::add_port
(
	PORT * port
)
{
	assert(port);

	if (port->get_io_direction() == PORT::INPUT)
	{
		m_input_ports.push_back(port);
	}
	else if (port->get_io_direction() == PORT::OUTPUT)
	{
		assert(! m_output_port);
		m_output_port = port;
	}
	else if (port->get_io_direction() == PORT::CLOCK)
	{
		m_input_ports.push_back(port);
	}
	else
	{
		Fail("Trying to add a port with unknown io direction");
	}	
}

// RETURNS: the D input port of the flip-flop
PORT * NODE::get_D_port() const
{
	assert(m_type == NODE::SEQ);
	assert(m_input_ports.size() == 2);
	
	return m_input_ports[1];

}

// RETURNS: the clock port of the flip-flop
PORT * NODE::get_clock_port() const
{
	assert(m_type == NODE::SEQ);
	assert(m_input_ports.size() == 2);
	assert(m_input_ports[0]->get_io_direction() == PORT::CLOCK);

	return m_input_ports[0];
}

// Removes the input port from the node
//
// PRE: port_to_remove is valid
// POST: port_to_remove has been removed from the node
//
void NODE::remove_input_port
(
	PORT * port_to_remove
)
{
	PORTS::iterator input_port_iter;

	input_port_iter = find(m_input_ports.begin(), m_input_ports.end(), port_to_remove);
	assert((*input_port_iter) == port_to_remove );

	debugif(DCODE, "Remove an input port " << port_to_remove->get_name()
			<< " in node " << m_name);

	m_input_ports.erase(input_port_iter);
}




//
// sets the lut for the node
//
// PRE: new_lut is valid
// POST: the lut has been set
//
void NODE::set_lut
(
 	LUT * new_lut
)
{
	// don't create mem leaks
	assert(m_lut);
	m_lut = new_lut;
}

// 
// RETURNS: a string with information about the node
//
string NODE::get_info() const
{
	string info_string = m_name;

	assert(m_type == NODE::COMB || m_type == NODE::SEQ);
	if (m_type == NODE::COMB)
	{
		info_string += " (comb)";
	}
	else 
	{
		info_string += " (seq)";
	}


	return info_string;
}

// 
// POST: a sub-cluster number has been added to the node
//
void NODE::add_sub_cluster_number(const CLUSTER_NUMBER_TYPE & cluster_number)
{
	debugif(DSUB_CLUSTER, "  Adding sub cluster number " << cluster_number << " to node " << get_info());
	m_sub_cluster_numbers.push_back(cluster_number);
}

// 
// POST: information about the node has been outputed
//
void NODE::print_out_information() const
{
	PORTS::const_iterator port_iter;
	PORT * port;

	Log("name " << m_name);

	for (port_iter = m_input_ports.begin(); 
		port_iter != m_input_ports.end(); port_iter++)
	{
		port = *port_iter;
		assert(port);

		Log("\tinput_port = " << port->get_name());
	}

	Log("Colour of Node = " << m_colour_mark);

}

EDGES NODE::get_output_edges() const 
{
	assert(m_output_port);
	return m_output_port->get_edges();
}

EDGES NODE::get_input_edges() const
{
	EDGES input_edges;
	PORTS::const_iterator port_iter;
	PORT * port  = 0;
	EDGE * edge  = 0;

	for (port_iter = m_input_ports.begin(); port_iter != m_input_ports.end(); port_iter++)
	{
		port = *port_iter;
		assert(port);
		edge = port->get_edge();
		assert(edge);

		input_edges.push_back(edge);
	}

	return input_edges;
}
EDGES NODE::get_input_edges_without_clock_edges() const
{
	EDGES input_edges;
	PORTS::const_iterator port_iter;
	PORT * port  = 0;
	EDGE * edge  = 0;

	for (port_iter = m_input_ports.begin(); port_iter != m_input_ports.end(); port_iter++)
	{
		port = *port_iter;
		assert(port);

		if (! port->is_clock_port())
		{
			edge = port->get_edge();
			assert(edge);

			input_edges.push_back(edge);
		}
	}

	return input_edges;
}

		
// 
// RETURNS: the Wirelength-approx of the node
//
NUM_ELEMENTS NODE::get_wirelength_approx_cost() const
{
	NUM_ELEMENTS total_wirelength_approx = 0,
	          	 wirelength_approx = 0;

	EDGES::const_iterator edge_iter;
	EDGE * edge = 0;
	NODE * source_node = 0;
	EDGES input_edges = get_input_edges_without_clock_edges();
	PORT * source_port = 0;
	
	for (edge_iter = input_edges.begin(); edge_iter != input_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		source_node = edge->get_source_node();

		if (source_node)
		{

			wirelength_approx = ABS(source_node->get_horizontal_position() - m_horizontal_position);
		}
		else
		{
			// we have a primary input

			source_port = edge->get_source();
			wirelength_approx = ABS(source_port->get_horizontal_position() - m_horizontal_position);
		}

		total_wirelength_approx += wirelength_approx;
	}

	return total_wirelength_approx;
}



//
// RETURNS: the number of nodes the node fanouts to
//
NUM_ELEMENTS NODE::get_fanout_degree() const
{
	NUM_ELEMENTS fanout_degree = 0;
	assert(m_output_port);

	fanout_degree = m_output_port->get_fanout_degree();

	return fanout_degree;
}
//
// RETURNS: the number of combinational nodes the node fanouts to
//
NUM_ELEMENTS NODE::get_fanout_degree_to_combinational_nodes() const
{
	NUM_ELEMENTS fanout_degree = 0;
	assert(m_output_port);

	fanout_degree = m_output_port->get_fanout_degree_to_combinational_nodes();

	return fanout_degree;
}

//
// RETURNS: the number of nodes that fanin to this node
//
NUM_ELEMENTS NODE::get_fanin_degree() const
{
	return m_input_ports.size();
}

PORT::PORT()
{
	m_port_type		= PORT::INTERNAL;
	m_external_type = PORT::NONE;
	m_io_direction	= UNKNOWN;
	m_my_node		= 0;
	m_cluster_number= -1;
}
PORT::PORT
(
	const string & port_name, 
	const PORT::PORT_TYPE & port_type, 
	const PORT::IO_DIRECTION & io_direction,
	const PORT::EXTERNAL_TYPE & external_type
)
{
	m_name			= port_name;
	m_port_type		= port_type;
	m_io_direction	= io_direction;
	m_external_type = external_type;
	m_my_node		= 0;
	m_cluster_number= -1;
}
PORT::PORT
(
	const string & port_name, 
	const PORT::PORT_TYPE & port_type, 
	const PORT::IO_DIRECTION & io_direction,
	const PORT::EXTERNAL_TYPE & external_type,
	NODE *	node_connected_to
)
{
	m_name			= port_name;
	m_port_type		= port_type;
	m_io_direction	= io_direction;
	m_external_type = external_type;
	m_my_node		= node_connected_to;
	m_cluster_number= -1;
	m_horizontal_position = 0;
}
PORT::PORT(const PORT & another_port)
{
	m_name			= another_port.m_name;
	m_port_type		= another_port.m_port_type;
	m_io_direction	= another_port.m_io_direction;
	m_external_type = another_port.m_external_type;
	m_my_node		= another_port.m_my_node;
	m_edges			= another_port.m_edges;
	m_cluster_number= another_port.m_cluster_number;
	m_horizontal_position = another_port.m_horizontal_position;
}
PORT & PORT::operator=(const PORT & another_port)
{
	m_name			= another_port.m_name;
	m_port_type		= another_port.m_port_type;
	m_io_direction	= another_port.m_io_direction;
	m_external_type = another_port.m_external_type;
	m_my_node		= another_port.m_my_node;
	m_edges			= another_port.m_edges;
	m_cluster_number= another_port.m_cluster_number;
	m_horizontal_position = another_port.m_horizontal_position;

	return (*this);
}
PORT::~PORT()
{
	
	// allow just the input ports to delete the edge
	if (m_io_direction == PORT::INPUT)
	{
		// nodes that are a constant 1 or 0 function have no edges
		if (! (m_edges.size() == 0))
		{
			assert(m_edges[0]);
			delete m_edges[0];
		}

	}
	
	m_edges.clear();

	m_my_node = 0;
}

//
// PRE: the port is an input edge (output edges have multiple edges)
// RETURNS: an edge
//
EDGE * PORT::get_edge() const
{
	return m_edges.front();
}

// 
// PRE: the port is an output edge
// RETURNS: the collection of edges attached to the port
// 
EDGES PORT::get_edges() const
{
	return m_edges;
}

//
// PRE: the port is an input edge (output edges have multiple edges)
// POST: the input edge is set
//
void PORT::set_edge
(
	EDGE * new_edge
)
{
	// this function should only be used with 
	// ports that can have a single edge

	assert(m_edges.size() == 0);
	assert(new_edge);
	m_edges.push_back(new_edge);
}

//
// PRE: the port is an output edge 
// POST: the edge has been added to the port
//
void PORT::add_edge
(
	EDGE * edge_to_add
)
{
	assert(edge_to_add);
	m_edges.push_back(edge_to_add);
}

//
// POST: the edge is removed from the port
void PORT::remove_edge
(
	EDGE * edge_to_remove
)
{
	EDGES::iterator ports_edge_reference;

	ports_edge_reference = find(m_edges.begin(), m_edges.end(), edge_to_remove);
	assert((*ports_edge_reference) == edge_to_remove );

	debugif(DCODE, "Remove a reference to edge " << edge_to_remove->get_name()
			<< " in port " << m_name);

	m_edges.erase(ports_edge_reference);
}
//
// RETURNS: the node that connects to this input port.
//          if there is no node (ie. the input port is connected to a PI) then return NULL
//
NODE * PORT::get_node_that_fanout_to_me() const
{
	EDGE * edge;
	PORT * output_port;
	NODE * output_node;

	assert(m_port_type == PORT::INTERNAL);
	assert(m_io_direction == PORT::INPUT || m_io_direction == PORT::CLOCK);

	edge = get_edge();
	assert(edge);

	output_port = edge->get_source();
	assert(output_port);
	output_node = output_port->get_my_node();

	// note: the node might be 0 if the output port is external
	return output_node;
}

//
// RETURNS: the output port that connects to this input port.
//
PORT *	PORT::get_output_port_that_fanout_to_me() const
{
	EDGE * edge;
	PORT * output_port;

	assert(m_port_type == PORT::INTERNAL);
	assert(m_io_direction == PORT::INPUT || m_io_direction == PORT::CLOCK);

	edge = get_edge();
	assert(edge);

	output_port = edge->get_source();
	assert(output_port);

	return output_port;
}

//
// RETURNS: the number of nodes the node fanouts to
//
NUM_ELEMENTS PORT::get_fanout_degree() const
{
	assert(m_io_direction == PORT::OUTPUT || m_port_type == PORT::EXTERNAL);

	return m_edges.size();
}


//
// RETURNS: the number of combinational nodes the node fanouts to
//
NUM_ELEMENTS PORT::get_fanout_degree_to_combinational_nodes() const
{
	assert(m_io_direction == PORT::OUTPUT);

	NUM_ELEMENTS fanout = 0;
	EDGES::const_iterator edge_iter;	
	EDGE * edge = 0;
	NODE * sink_node = 0;

	for (edge_iter = m_edges.begin(); edge_iter != m_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		sink_node = edge->get_sink_node();

		if (sink_node->get_type() != NODE::SEQ)
		{
			fanout++;
		}
	}

	return fanout;
}

//
// RETURNS: a dff in the fanout if one exists, else NULL
//          
NODE *  PORT::find_dff_in_fanout() const
{
	NODE * dff_found = 0;
	EDGES::const_iterator edge_iter;
	NODE * node_in_fanout;

	edge_iter = m_edges.begin();
	while ( (! dff_found) && (edge_iter != m_edges.end()) )
	{
		assert(*edge_iter);	
		node_in_fanout = (*edge_iter)->get_sink_node();

		if (node_in_fanout->get_type() == NODE::SEQ)
		{
			dff_found = node_in_fanout;
		}
		edge_iter++;
	}

	return dff_found;
}


bool PORT::is_connected_to_PI() const
{
	NODE * source_node = get_node_that_fanout_to_me();
	
	if (source_node)
	{
		// primary inputs are not connected to source nodes
		return false;
	}
	else
	{
		return true;
	}
}

