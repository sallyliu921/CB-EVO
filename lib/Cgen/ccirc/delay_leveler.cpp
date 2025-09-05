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





#include "delay_leveler.h"

DELAY_LEVELER::DELAY_LEVELER()
{
	m_max_combinational_delay = -1;

	m_circuit = 0;
}
DELAY_LEVELER::DELAY_LEVELER(const DELAY_LEVELER & another_delay_leveler)
{
	assert(false);
	m_max_combinational_delay	= another_delay_leveler.m_max_combinational_delay;
	m_circuit					= another_delay_leveler.m_circuit;
}

DELAY_LEVELER & DELAY_LEVELER::operator=(const DELAY_LEVELER & another_delay_leveler)
{
	assert(false);

	m_max_combinational_delay	= another_delay_leveler.m_max_combinational_delay;
	m_circuit					= another_delay_leveler.m_circuit;

	return (*this);
}

DELAY_LEVELER::~DELAY_LEVELER()
{
}

//
// Calculate and label the nodes in the graph with their delay level
//
// PRE:  circuit is valid
// POST: Each node has its delay level labelled
//       Edge edge has its length labelled
//       The maximum combinational delay of the circuit has been found
//       The sequential_level data structure has organized the nodes by delay level
//
void DELAY_LEVELER::calculate_and_label_combinational_delay_levels
(
	CIRCUIT * circuit
)
{
	assert(circuit);
	m_circuit = circuit;

	SEQUENTIAL_LEVEL * sequential_level = m_circuit->get_sequential_level();
	assert(sequential_level);

	debug("Status: Beginning combinational delay analysis.");

	m_circuit->colour_graph(NODE::UNMARKED);

	calculate_and_label_combinational_delay_for_sequential_level(sequential_level);

	// add the primary inputs and nodes to the sequential level datastructure
	add_primary_inputs_to_0th_sequential_level(sequential_level);
	add_nodes_to_sequential_level(sequential_level);

	label_edges_with_length();

	check_sanity();

	debug("Status: Combinational delay analysis is complete.");
}



//
// 	
// 	Calculate and label each node in the graph with its delay level.
//
//	PRE: sequential_level is valid
//	POST: all nodes have their delay level labelled
//	      the max_combinational_delay of the sequential level has been found and set

void DELAY_LEVELER::calculate_and_label_combinational_delay_for_sequential_level
(
	SEQUENTIAL_LEVEL * sequential_level
)
{
	NODES::iterator node_iter;
	NODE * top_node = 0;
	NODE_PTR_DEQUE nodes_to_find_delay_for;
	NODE * node = 0;
	NODES top_nodes = m_circuit->get_top_nodes();
	m_max_combinational_delay = 0;

	assert(sequential_level);

	label_all_dff_with_0_delay();

	for (node_iter = top_nodes.begin(); node_iter != top_nodes.end(); node_iter++)
	{
		top_node = *node_iter;
		assert(top_node);

		if (top_node->get_type() == NODE::SEQ)
		{
			// if the node is connected to a SEQ 
			// it should already have its delay and colour set
			assert(	top_node->get_colour() == NODE::MARKED &&
					top_node->get_max_comb_delay_level() == 0);
		}
		else
		{
			// the node is combinational
	
			// put the node on the queue of nodes we need to find delay for
			nodes_to_find_delay_for.push_back(top_node);

			// continue to find the delay level of nodes until the queue of nodes is empty
			while (! nodes_to_find_delay_for.empty())
			{
				// pop a node off the stack
				node = nodes_to_find_delay_for.front();
				nodes_to_find_delay_for.pop_front();

				// if it is unmarked look for its delay level
				if (node->get_colour() == NODE::UNMARKED)
				{
					find_comb_delay_of_node(node, nodes_to_find_delay_for);
				}
			}
		}
	}

	// set the maximum combinational delay we saw during labelling
	sequential_level->set_max_combinational_delay(m_max_combinational_delay);
}

//
// Find and set the delay level of the node
// 
// PRE: node is valid
// POST: the delay level of this nodes has been set and the delay level of all nodes in the fanin
//       of this node has been set
//
void DELAY_LEVELER::find_comb_delay_of_node
(
	NODE * node,
	NODE_PTR_DEQUE & nodes_to_find_delay_for
)
{
	assert(node && node->get_colour() == NODE::UNMARKED);

	PORTS input_ports = node->get_input_ports();
	PORTS::iterator port_iter;
	PORT * port = 0;
	DELAY_TYPE max_comb_delay_of_fanin = 0;
	NODE * fanin_node = 0;

	// push the unmarked fanout onto the queue
	get_unmarked_fanout(node, nodes_to_find_delay_for);

	// look for the maximum combinational delay of the fanin 
	for (port_iter = input_ports.begin(); port_iter != input_ports.end(); port_iter++)
	{
		port = *port_iter;

		fanin_node = port->get_node_that_fanout_to_me();

		// note: if the fanin node is NULL it is a PI 
		// We don't care about PI because their combinational delay is 0 
		// and we are looking for the maximum combinational delay
		if (fanin_node)
		{
			// if the node is combinational then it might have combinal delay 
			// else the node is a dff and does not have combinational delay
			if (fanin_node->get_type() == NODE::COMB)
			{

				// if this node has not had its delay level calculated. calculate it.
				if (fanin_node->get_colour() == NODE::UNMARKED)
				{
					find_comb_delay_of_node(fanin_node, nodes_to_find_delay_for);
				}
			
				// the fanin node should now have its max combinational delay determined
				assert(fanin_node->get_colour() == NODE::MARKED && fanin_node->get_max_comb_delay_level() >= 0);
				
				max_comb_delay_of_fanin = MAX(max_comb_delay_of_fanin, fanin_node->get_max_comb_delay_level());
			}
				
		}
	}

	node->set_max_comb_delay_level(max_comb_delay_of_fanin+1);
	node->set_colour(NODE::MARKED);

	m_max_combinational_delay = MAX(m_max_combinational_delay, max_comb_delay_of_fanin+1);
}

//
// gets all unmarked nodes in the fanout of this node and 
//
// PRE: node is valid
// POST: all unmarked combinational nodes in the fanin of this node are 
//       enqueued in nodes_to_find_delay_for
//
void DELAY_LEVELER::get_unmarked_fanout
(
	NODE * node, 
	NODE_PTR_DEQUE & nodes_to_find_delay_for
)
{

	assert(node);

	EDGES output_edges = node->get_output_edges();
	EDGES::iterator edge_iter;
	EDGE * edge = 0;
	NODE * fanout_node = 0;

	for (edge_iter = output_edges.begin(); edge_iter != output_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		fanout_node = edge->get_sink_node();
		assert(fanout_node);

		if (fanout_node->get_colour() == NODE::UNMARKED  && fanout_node->get_type() == NODE::COMB)
		{
			nodes_to_find_delay_for.push_back(fanout_node);
		}
		else
		{
			// the node should either be a dff or MARKED
			assert( fanout_node->get_type() == NODE::SEQ || fanout_node->get_colour() == NODE::UNMARKED);
		}
	}
}

//
// sanity check on the circuit to make sure everything is ok
//
// PRE: circuit is valid
// POST: all nodes have their delay level checked
//
void DELAY_LEVELER::check_sanity() const
{
	assert(m_circuit);

	NODES & nodes = m_circuit->get_nodes();
	NODES::const_iterator node_iter;
	NODE * node = 0;
	DELAY_TYPE node_comb_delay = -1;

	for (node_iter = nodes.begin(); node_iter != nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		node_comb_delay = node->get_max_comb_delay_level();


		if (node->get_type() == NODE::SEQ)
		{
			assert(node_comb_delay == 0);
		}
		else
		{
			assert(node->get_colour() == NODE::MARKED);
			assert(node_comb_delay == 1 + get_max_comb_delay_level_of_fanin(node));
		}
	}


    debug("Sanity checks on delay levels done (successful) ...");
}	

//
// Find the max combinational delay of the fanin to this node
//
// PRE: node is valid
// Returns: the maximum combinational delay of the fanin
//
LEVEL_TYPE DELAY_LEVELER::get_max_comb_delay_level_of_fanin
(
	const NODE * node
) const
{

	assert(node);

	PORTS input_ports = node->get_input_ports();
	PORTS::const_iterator port_iter;
	PORT * port = 0;
	DELAY_TYPE max_comb_delay = 0;
	NODE * fanin_node = 0;


	// look at the combinational delay of the fanin
	for (port_iter = input_ports.begin(); port_iter != input_ports.end(); port_iter++)
	{
		port = *port_iter;

		fanin_node = port->get_node_that_fanout_to_me();

		if (fanin_node)
		{
			if (fanin_node->get_type() == NODE::COMB)
			{
				max_comb_delay = MAX(max_comb_delay, fanin_node->get_max_comb_delay_level());
			}
		}
	}

	return max_comb_delay;
}



//
// Label each edge in the graph with the number of delay levels it crosses
// For edges that connect to the input of a flip-flop label that edge of length 0.
//
// PRE: circuit is valid
//      nodes have their delay level labelled
// POST: all edges in the circuit have their length labelled
//
void DELAY_LEVELER::label_edges_with_length()
{
	assert(m_circuit);
	EDGES & edges =  m_circuit->get_edges();
	EDGES::const_iterator edge_iter;
	EDGE * edge = 0;
	NODE * source_node = 0;
	NODE * sink_node = 0;
	LENGTH_TYPE edge_length = 0;
	DELAY_TYPE source_delay = 0;
	DELAY_TYPE sink_delay = 0;


	for (edge_iter = edges.begin(); edge_iter != edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		source_node = edge->get_source_node();
		sink_node = edge->get_sink_node();

		if (source_node)
		{
			if (sink_node->get_type() == NODE::COMB)
			{
				source_delay = source_node->get_max_comb_delay_level();
				sink_delay = sink_node->get_max_comb_delay_level();

				edge_length = sink_delay - source_delay;
			}
			else
			{
				// the edge is an edge to a dff.
				// the length of this edge is 0
				edge_length = 0;
			}
		}
		else
		{
			// we have a pi

			sink_delay = sink_node->get_max_comb_delay_level();

			edge_length = sink_delay;
		}

		assert(edge_length >= 0);

		edge->set_length(edge_length);
	}
}

//
// Label all the dff with a combinational delay of 0
//
// PRE: circuit is valid
// POST: all dff have a combinational delay of 0 and have been coloured as 
//       marked
void DELAY_LEVELER::label_all_dff_with_0_delay()
{
	assert(m_circuit);
	NODES dff_list = m_circuit->get_dffs();
	NODES::iterator node_iter;
	NODE * node = 0;

	for (node_iter = dff_list.begin(); node_iter != dff_list.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);
		assert(node->get_type() == NODE::SEQ);
		node->set_max_comb_delay_level(0);
		node->set_colour(NODE::MARKED);
	}
}

//
// add the nodes to the sequential level datastructure
// which organizes the nodes by delay level
//
// PRE: each node has its delay level defined
// POST: sequential_level has all the nodes in the circuit organized by delay level
// 
void DELAY_LEVELER::add_nodes_to_sequential_level
(
	SEQUENTIAL_LEVEL*  sequential_level
)
{
	assert(m_circuit);

	NODES & nodes = m_circuit->get_nodes();
	NODES::iterator node_iter;
	NODE * node = 0;

	for (node_iter = nodes.begin(); node_iter != nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		sequential_level->add_node(node);
	}

}

// 
// add the primary inputs to the sequential level datastructure
//
// PRE: seq_level_0 is valid
// POST: PIs have been added to seq_level_0
//
void DELAY_LEVELER::add_primary_inputs_to_0th_sequential_level
(
	SEQUENTIAL_LEVEL * seq_level_0
)
{

	PORTS primary_inputs = m_circuit->get_PI();
	PORTS::iterator port_iter;
	PORT * port = 0;

	for (port_iter = primary_inputs.begin(); port_iter != primary_inputs.end(); port_iter++)
	{
		port = *port_iter;
		assert(port);

		seq_level_0->add_primary_input(port);
	}
}
