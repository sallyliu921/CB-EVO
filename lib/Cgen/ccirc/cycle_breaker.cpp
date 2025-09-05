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





#include "cycle_breaker.h"



CYCLE_BREAKER::CYCLE_BREAKER()
{
	m_circuit			=	0;
	m_number_new_dff	=	0;
	m_number_moved_edges=	0;
	m_shown_breaking_cycles_warning = false;
}
CYCLE_BREAKER::CYCLE_BREAKER(const CYCLE_BREAKER & another_cycle_breaker)
{
	m_circuit			= another_cycle_breaker.m_circuit;
	m_number_new_dff	= another_cycle_breaker.m_number_new_dff;
	m_number_moved_edges= another_cycle_breaker.m_number_moved_edges;
	m_shown_breaking_cycles_warning = another_cycle_breaker.m_shown_breaking_cycles_warning;
}
CYCLE_BREAKER & CYCLE_BREAKER::operator=(const CYCLE_BREAKER & another_cycle_breaker)
{
	m_circuit			= another_cycle_breaker.m_circuit;
	m_number_new_dff	= another_cycle_breaker.m_number_new_dff;
	m_number_moved_edges= another_cycle_breaker.m_number_moved_edges;
	m_shown_breaking_cycles_warning = another_cycle_breaker.m_shown_breaking_cycles_warning;

	return (*this);
}
CYCLE_BREAKER::~CYCLE_BREAKER()
{
}

//
// Detects and breaks any combinational cycles that may exist in the graph
// PRE: circuit is valid
// POST: all combinational cycles have been broken
//
void CYCLE_BREAKER::break_cycles
(
	CIRCUIT * circuit
)
{
	assert(circuit);

	bool found_cycle = true;
	PORTS::iterator port_iter;
	NODE_PTR_LIST next_nodes;
	PORTS PI;

	m_number_new_dff = 0;
	m_number_moved_edges = 0;
	m_circuit = circuit;

	PI = m_circuit->get_PI();
	port_iter = PI.begin();

	// for each PI break cycles connected to it
	while (port_iter != PI.end() )
	{
		next_nodes.clear();
		debugsepif(DCYCLE);
		debugif(DCYCLE,"Looking for a cycle in the graph");
		
		// colour all the nodes as unexplored
		m_circuit->colour_graph(NODE::WHITE);

		// find and break cycles in the graph
		found_cycle = find_and_break_cycles(*port_iter, next_nodes, 0);

		// if we didn't find a cycle look to the next PI, 
		if (! found_cycle)
		{
			// we should not have any unexplored flip-flop nodes
			assert(next_nodes.empty());
			port_iter++;
		}
		else
		{
			// we didn't find a cycle. don't look to the next PI.
			// see if we have further cycles
			show_breaking_cycles_warning();
		}
	}

	show_number_dff_created_edges_moved_warning();
	debugSep;

}

/*
 **  Search the graph, looking for a cycle, from specified node.
 **
 **  Variation on the white/grey/black algorithm (see CLR book).  Do
 **  a DFS, with unexplored nodes white, active nodes grey and explored
 **  black.  A combinational cycle exists when a grey node uncovers a 
 **  grey fanout.  Variation is that DFFs are queued for later processing,
 **  rather than explored DFS; because comb cycles can't go through them.
 **/

bool CYCLE_BREAKER::find_and_break_cycles
(
	PORT * output_port,
	NODE_PTR_LIST & next_nodes,
	DEPTH_TYPE depth
)
{
	EDGES output_edges;
	EDGES::iterator edge_iter;
	EDGE * output_edge;
	NODE * node_in_fanout;
	PORT * node_in_fanout_output_port;
	NODE::COLOUR_TYPE colour;
	bool found_cycle = false;

	// get the edges that fanout from the output port
	assert(output_port);
	output_edges = output_port->get_edges();

	debugif(DCYCLE,"Looking at fanout from port " << output_port->get_name()); 

	// explore all these edges 
	for (edge_iter = output_edges.begin(); edge_iter != output_edges.end(); edge_iter++)
	{
		output_edge = (*edge_iter);
		assert(output_edge);

		node_in_fanout = output_edge->get_sink_node();
		assert(node_in_fanout);

		colour = node_in_fanout->get_colour();

		// what colour is the node that we fanout too?
	
		if (colour == NODE::BLACK)
		{
			// OK, not on a combinational cycle
			debugif(DCYCLE,"Found node " << node_in_fanout->get_name() << 
					" marked as black");
		}
		else if (colour == NODE::WHITE)
		{
			// we have not seen this node before
			if (node_in_fanout->get_type() == NODE::SEQ)
			{
				// if the node is a dff then this edge is not on a combinational cycle.
				// store the dff and we will come back to it later
				debugif(DCYCLE,"Queuing " << node_in_fanout->get_name() << " as node " <<
						next_nodes.size() << " for later");
				next_nodes.push_back(node_in_fanout);
			}
			else
			{
				debugif(DCYCLE,"...marking " << node_in_fanout->get_name() << " as grey");
				node_in_fanout->set_colour(NODE::GREY);
				node_in_fanout_output_port = node_in_fanout->get_output_port();

				// look for cycles that fanout from this node
				found_cycle = find_and_break_cycles(node_in_fanout_output_port, next_nodes,depth+1);
				if (found_cycle)
				{
					debugif(DCYCLE,"returning...");
					return true;
				}
				debugif(DCYCLE,"Marking " << node_in_fanout->get_name() << " as black");
				node_in_fanout->set_colour(NODE::BLACK);
			}
		}
		else if (colour == NODE::GREY)
		{
			// we found a cycle. break it.

			Verbose("Breaking combinational cycle by inserting DFF on edge "
					<< output_edge->get_name());
        	break_cycle_by_adding_dff(output_port, output_edge);
			debugif(DCYCLE,"returning...");
			return true;
		}
		else
		{
			assert(false);
		}
	}

	// we only get here if a cycle wasn't found above
	assert(! found_cycle);

	found_cycle = find_cycle_in_path_from_seq_nodes(next_nodes, depth);

	return found_cycle;
}

//
// look for cycles that exists in the fanout from the flip-flops
bool CYCLE_BREAKER::find_cycle_in_path_from_seq_nodes
(
	NODE_PTR_LIST & next_nodes,
	DEPTH_TYPE depth
)
{
	PORT * seq_node_output_port;
	NODE * seq_node;
	bool found_cycle = false;

	while (depth == 0 && ! next_nodes.empty())
	{
		seq_node = next_nodes.front();
		next_nodes.pop_front();
		
		// seq->comb->seq cycles are ok
		debugsepif(DCYCLE);
		debugif(DCYCLE,"Working from sequential node " << seq_node->get_name());
		debugif(DCYCLE,"...marking " << seq_node->get_name() << " as black");
		seq_node->set_colour(NODE::BLACK);

		seq_node_output_port = seq_node->get_output_port();
		found_cycle = find_and_break_cycles(seq_node_output_port, next_nodes, 0);

		if (found_cycle)
		{
			debugif(DCYCLE,"returning...");
			return true;
		}
	}

	return false;
}

// breaks a combinational cycle by adding inserting a DFF in the cycle.
//
//
// source_port				        source_port
//	    |                                |
//      | <-- edge_to_break              |
//      |                                |
//    destination_input_port     dff that gets inserted
//                  					 |
//					                     |
//				                 destination_input_port   
//
// PRE: a combinational cycle exists 
// POST: the combinational cycle has been broken by inserting a dff

void CYCLE_BREAKER::break_cycle_by_adding_dff
(
	PORT * source_port, 
	EDGE * edge_to_break
)
{
	NODE * dff_to_insert = 0;
	string new_latch_name;
	PORT * destination_input_port;

	assert(source_port);
	assert(edge_to_break);

	destination_input_port = edge_to_break->get_sink();
	assert(destination_input_port);

	// remove from the source and dest. ports the edge that we want to break
	m_circuit->remove_and_delete_edge(source_port,destination_input_port,edge_to_break);

	debug("Adding new flip-flop to the graph"); 

	// first look for an existing dff off the source port
	dff_to_insert = source_port->find_dff_in_fanout();

	if (dff_to_insert)
	{
		// found one!
		debug("Using the existing latch " << dff_to_insert->get_name());
		m_circuit->create_edge(dff_to_insert->get_output_port(),
								destination_input_port,1);
		m_number_moved_edges++;
	}
	else
	{
		// didn't find one.
		new_latch_name = source_port->get_name() + "_L";
		debug("Creating a new latch " << new_latch_name);
		dff_to_insert = m_circuit->create_dff(new_latch_name,source_port->get_name());
		m_circuit->create_edge(source_port, dff_to_insert->get_D_port(), 1);
		m_circuit->create_edge(dff_to_insert->get_output_port(), destination_input_port,1);

		m_number_new_dff++;
	}
	// rename the destination input port with the name of the output node
	destination_input_port->set_name(dff_to_insert->get_name());

	show_added_dff_warning();
}

void CYCLE_BREAKER::show_added_dff_warning()
{
	debug("Adding new latch to graph with size=" << m_circuit->get_size() <<
		  "\tn= " << m_circuit->get_nNodes() <<
		  "\te= " << m_circuit->get_nEdges() << 
		  "\tnumber dff = " << m_circuit->get_nDFF());
}

void CYCLE_BREAKER::show_breaking_cycles_warning()
{
	if ((m_number_new_dff + m_number_moved_edges == 1) && ! m_shown_breaking_cycles_warning)
	{
		Warning("Breaking combinational cycles (set 'verbose' to watch).");
		m_shown_breaking_cycles_warning = true;
	}
}

void CYCLE_BREAKER::show_number_dff_created_edges_moved_warning()
{
	if (m_number_new_dff + m_number_moved_edges > 0) 
	{
		Warning("Added " << m_number_new_dff << " dff, moved " << m_number_moved_edges
				<< " edges to make graph cycle-free.");
	}
}
