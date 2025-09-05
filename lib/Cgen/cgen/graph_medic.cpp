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




#include "graph_medic.h"
#include <algorithm>

GRAPH_MEDIC::GRAPH_MEDIC(CIRCUIT_GRAPH * circuit)
{
	m_graph					= circuit;
	m_problem_nodes			= GRAPH_MEDIC::NODES_NOT_CONNECTED_TO_PO;
	m_number_comb_deleted 	= 0;
	m_number_seq_deleted	= 0;
	m_number_unreachable	= 0;
	m_number_isolated		= 0;

	m_shown_eliminate_node_warning = false;
}


GRAPH_MEDIC::GRAPH_MEDIC(const GRAPH_MEDIC & another_graph_medic)
{
	m_graph				= another_graph_medic.m_graph; 
	m_nodes_to_delete	= another_graph_medic.m_nodes_to_delete;

	m_number_comb_deleted 	= another_graph_medic.m_number_comb_deleted;
	m_number_seq_deleted	= another_graph_medic.m_number_seq_deleted;
	m_number_unreachable	= another_graph_medic.m_number_unreachable;
	m_number_isolated		= another_graph_medic.m_number_isolated;

	m_problem_nodes		= another_graph_medic.m_problem_nodes;
	m_shown_eliminate_node_warning = another_graph_medic.m_shown_eliminate_node_warning;
}
	
GRAPH_MEDIC & GRAPH_MEDIC::operator=(const GRAPH_MEDIC & another_graph_medic)
{
	m_graph				= another_graph_medic.m_graph; 
	m_nodes_to_delete	= another_graph_medic.m_nodes_to_delete;

	m_number_comb_deleted 	= another_graph_medic.m_number_comb_deleted;
	m_number_seq_deleted	= another_graph_medic.m_number_seq_deleted;
	m_number_unreachable	= another_graph_medic.m_number_unreachable;
	m_number_isolated		= another_graph_medic.m_number_isolated;

	m_problem_nodes		= another_graph_medic.m_problem_nodes;
	m_shown_eliminate_node_warning = another_graph_medic.m_shown_eliminate_node_warning;

	return (*this);
}

GRAPH_MEDIC::~GRAPH_MEDIC()
{
}

// PRE: graph exists
// POST: all nodes are reachable from the PI or PO
//       no unconnected nodes or nodes exists
//
void GRAPH_MEDIC::delete_unusable_nodes()
{
	// delete the unconnected external nodes
	m_graph->delete_unconnected_nodes();

	// look for nodes that are unreachable from the PO
	delete_unreachable_nodes_from_outputs();

	// look for nodes that are unreachable from the PI
	delete_unreachable_nodes_from_inputs();

	// assert: all unreachable nodes are in m_nodes_to_delete

	delete_queued_nodes();

	// delete the unconnected external nodes
	m_graph->delete_unconnected_nodes();

	print_node_deletion_stats();
}

//
// Delete any nodes enqueued on m_nodes_to_delete
// POST: all nodes in m_nodes_to_delete have been deleted
//       deleted node counters have been updated
//
void GRAPH_MEDIC::delete_queued_nodes()
{
	NODE * node_to_delete;
	NODES::iterator node_iter;

	while(! m_nodes_to_delete.empty() )
	{
		// pop the first node off the list
		node_to_delete = m_nodes_to_delete.back();
		m_nodes_to_delete.pop_back();

		assert(node_to_delete);


		// remove all other references from the queque
		node_iter = find(m_nodes_to_delete.begin(), m_nodes_to_delete.end(), node_to_delete);
		while (node_iter != m_nodes_to_delete.end())
		{
			m_nodes_to_delete.erase(node_iter);	
			node_iter = find(m_nodes_to_delete.begin(), m_nodes_to_delete.end(), node_to_delete);
		}

		if (node_to_delete->is_combinational())
		{
			m_number_comb_deleted++;
		}
		else
		{
			m_number_seq_deleted++;
		}

		delete_node(node_to_delete);
	}
}

// this function looks at all nodes to see if they are buffer/inverter nodes that can be deleted.
// A buffer/inverter node can be delete if its source node has a fanout of 1 and is internal.
//
// Note: these nodes should have been removed during tech mapping but for whatever 
// reason they sometimes still exist.
//
// This function might change the name of the output nodes if the buffer node feeds a primary output.
// This can be avoided if we rename the source node to the name of the primary output but
// this complicates matters plus I am not sure if it is the right thing to do.
// If a dff fed out to a buffer and the buffer was a PO the dff would be renamed.
//
// Which is better. To rename the outputs or to rename the comb/dff?
//
//
//          0		- output_node_of_source_node: internal and nOutputs = 1
//          |
//			O		- node: buffer node, needs to be removed
//         /|
//		  O	O 		- nodes in fanout:
//
//
//
//			O		- output_node_of_source_node
//         /|
//        O	O		- nodes in fanout 
//
// PRE: buffer/inverter nodes may exist
// POST: all deletable buffer/inverter nodes have been deleted
//
void GRAPH_MEDIC::delete_buffer_and_inverter_nodes()
{

	assert(new_gen);
	//unimplemented in cgen
}





//
// Delete any unreachable node from the primary outputs
//
// PRE: nodes might exist that are unreachable from the POs
// POST: nodes that are unreachable from the POs have been 
//       queue for deletion
//
void GRAPH_MEDIC::delete_unreachable_nodes_from_outputs()
{

	m_problem_nodes = GRAPH_MEDIC::NODES_NOT_CONNECTED_TO_PO;

	debugif(DMEDIC, "Marking nodes up from the PO");
	mark_up(NODE::MARKED);

	debugif(DMEDIC, "Queing for deletion any unmarked node down from the PI");
	eliminate_down();
	
	debugif(DMEDIC, "Unmarking the nodes up from the PO");
	mark_up(NODE::UNMARKED);
}
//
// Colours the nodes up from the primary outputs
//
// PRE: the graph has been created
// POST: all nodes reachable from the POs are coloured
//
void GRAPH_MEDIC::mark_up
(
	const NODE::COLOUR_TYPE & colour
)
{
	NODES::iterator node_iter;
	NODE * po_node;
	NODES PO = m_graph->get_PO();
	
	// mark the nodes up from the primary outputs
	for (node_iter = PO.begin(); node_iter != PO.end(); node_iter++)
	{
		po_node = *node_iter;
		assert(po_node);

		if (po_node->get_colour() != colour)
		{
			mark_up_from_node(po_node, colour);
		}
	}

}

//
// Colours the nodes up from this node
//
// PRE: node is valid
// POST: all nodes up from this node that are reachable are coloured
//
void GRAPH_MEDIC::mark_up_from_node
(
	NODE * node,
	const NODE::COLOUR_TYPE & colour
)
{
	EDGES	input_edges;
	NODE *	source_node = 0;
	EDGES::iterator edge_iter;

	assert(node);
	node->set_colour(colour);

	debugif(DMEDIC, node->get_name() << " = marked " << static_cast<short>(colour) );

	input_edges = node->get_input_edges();

	for (edge_iter = input_edges.begin(); edge_iter != input_edges.end(); edge_iter++) 
	{
		assert(*edge_iter);
		source_node = (*edge_iter)->get_source_node();
		assert(source_node);

		if (source_node->get_colour() != colour)
		{
			mark_up_from_node(source_node, colour);
		}
	}
}
//
// Eliminates any UNMARKED node in the graph starting 
// by added that node to a list of nodes queue for deletion
//
// Starts from the PI
//
// PRE: nodes in the graph may be marked UNMARKED 
// POST: all nodes reachable from the PI marked UNMARKED
//       have been queue for deletion
//
void GRAPH_MEDIC::eliminate_down()
{
	NODES::iterator node_iter;
	NODES PI = m_graph->get_PI();
	NODE * pi_node;

	// delete any unmarked node down from the PIs
	for (node_iter = PI.begin(); node_iter != PI.end(); node_iter++) 
	{
		pi_node = *node_iter;
		assert(pi_node);

		if_unmarked_queue_for_deletion(pi_node);
		eliminate_down_from_node(pi_node);
	}
}


//
// Eliminates any UNMARKED node in the graph starting 
// down from this output node
//
// PRE: nodes in the graph may be marked UNMARKED 
// POST: all nodes reachable from this output node down 
//       and marked UNMARKED have been queue for deletion
void GRAPH_MEDIC::eliminate_down_from_node
(
	NODE * node
)
{
	EDGES edges;	
	EDGES::iterator edge_iter;	
	EDGE * edge;
	NODE * sink_node = 0;

	assert(node);
	edges = node->get_output_edges();
	
	for(edge_iter = edges.begin(); edge_iter != edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		sink_node = edge->get_sink_node();
		assert(sink_node);

		// if we haven't seen the node, process it
		if (sink_node->get_colour() == NODE::MARKED || sink_node->get_colour() == NODE::UNMARKED)
		{	
			debugif(DMEDIC, "Eliminate down:  searching at " <<  node->get_name());
			if_unmarked_queue_for_deletion(sink_node);
			eliminate_down_from_node(sink_node);
		}
	}
}


//
// Deletes any node unreachable from the primary inputs
// 
// PRE: nodes might exist that are unreachable from the PIs
// POST: nodes that are unreachable from the PIs have been 
//       queue for deletion
void GRAPH_MEDIC::delete_unreachable_nodes_from_inputs()
{

	m_problem_nodes = GRAPH_MEDIC::NODES_NOT_CONNECTED_TO_PI;

	debugif(DMEDIC, "Marking nodes down from the PI");
	mark_down(NODE::MARKED);

	debugif(DMEDIC, "Queing for deletion up from the PO");
	eliminate_up();

	// reset the graph
	debugif(DMEDIC, "Unmarking nodes down from the PI");
	mark_down(NODE::UNMARKED);
}

//
// Colours node down from the PIs
//
// PRE: node is valid
// POST: all nodes reachable from the PIs are coloured
void GRAPH_MEDIC::mark_down
(
	const NODE::COLOUR_TYPE & colour
)
{ 
	NODES PI = m_graph->get_PI();
	NODES::iterator node_iter;
	NODE * pi_node;

	// mark the nodes down from the PI
	for (node_iter = PI.begin(); node_iter != PI.end(); node_iter++) 
	{
		pi_node = *node_iter;
		assert(pi_node);

		mark_down_from_node(pi_node, colour);
	}
}


//
// Colour down from this output node
//
// PRE: node is valid
// POST: all nodes down from this node that are reachable are coloured
void GRAPH_MEDIC::mark_down_from_node
(
	NODE * node,
	const NODE::COLOUR_TYPE & colour
)

{
	NODE * sink_node;
	EDGES edges;
	EDGES::iterator edge_iter;

	node->set_colour(colour);

	assert(node);
	edges = node->get_output_edges();

	// over all the edges look for more nodes to mark
	for (edge_iter = edges.begin(); edge_iter != edges.end(); edge_iter++) 
	{
		assert(*edge_iter);
		
		sink_node = (*edge_iter)->get_sink_node();
		assert(sink_node);

		// if we haven't visted this node before colour it
		if (sink_node->get_colour() != colour)
		{
			debugif(DMEDIC, node->get_name() << " = marked " << static_cast<short>(colour) ); 
			node->set_colour(colour);

			mark_down_from_node(sink_node, colour);
		}
	}
}

//
// Queues for deletion any UNMARKED node
//
// PRE: nodes in the graph may be marked UNMARKED 
// POST: all nodes reachable from the PO marked UNMARKED
//       have been queue for deletion
void GRAPH_MEDIC::eliminate_up()
{
	NODES::iterator node_iter;
	NODES PO = m_graph->get_PO();
	NODE * po_node;

	// delete up from the PO
	for (node_iter = PO.begin(); node_iter != PO.end(); node_iter++) 
	{
		po_node = *node_iter;
		assert(po_node);

		if (po_node->get_colour() != NODE::MARKED_VISITED &&
		   po_node->get_colour() != NODE::UNMARKED_VISITED)
		{
			if_unmarked_queue_for_deletion(po_node);
			eliminate_up_from_node(po_node);
		}
	}
}

// PRE: nodes in the graph may be marked UNMARKED 
// POST: all nodes reachable from this output node up  
//       and marked UNMARKED have been queue for deletion
void GRAPH_MEDIC::eliminate_up_from_node
(
	NODE * node
)
{
	EDGES::iterator edge_iter;
	EDGE *	edge;
	NODE *	source_node;

	EDGES input_edges = node->get_input_edges();

	for (edge_iter = input_edges.begin(); edge_iter != input_edges.end(); edge_iter++) 
	{
		// input nodes should be connected to an edge and the edge should
		// be connected to an output node

		edge = *edge_iter;
		assert(edge);
		
		source_node	 = edge->get_source_node();
		assert(source_node);

		// if we haven't seen it before process it
		if ((source_node->get_colour() == NODE::MARKED) ||
			(source_node->get_colour() == NODE::UNMARKED))
		{
			if_unmarked_queue_for_deletion(source_node);
			eliminate_up_from_node(source_node);
		}
	}
}


//
// If the node is unmarked queue the node for deletion
//
// Mark node as visited in both cases
// 
// PRE: node is valud and marked MARKED or UNMARKED
// POST: If the node was unmarked it is queued for deletion
//       node is marked visited
void GRAPH_MEDIC::if_unmarked_queue_for_deletion
(
	NODE * node
)
{
	assert(node);

	assert(	node->get_colour() == NODE::MARKED || 
			node->get_colour() == NODE::UNMARKED || node->get_colour() == NODE::NONE);

	if (node->get_colour() == NODE::MARKED) 
	{
		node->set_colour(NODE::MARKED_VISITED);
		
	}
	else
	{
		node->set_colour(NODE::UNMARKED_VISITED);
		show_node_deletion_warning(node);
		Verbose("Adding " << node->get_name() << " to the deletion list");

		m_nodes_to_delete.push_back(node);
		m_number_unreachable++;
	}
}

//
// Delete the node from all data structures
//
// PRE: node is valid
// POST: node has been deleted
void GRAPH_MEDIC::delete_node
(
	NODE * node
)
{

	assert(node);

	// detach the node from both its inputs and outputs
	detach_node_from_references_in_fanin(node);
	detach_node_from_references_in_fanout(node);

	m_graph->delete_node(node);
	node = 0;

}

//
// Detach node from fanin
//
// PRE: node may have fanin
// POST: the node has nothing than fans into it
// 
// to detach the node 01 we must remove its reference from the output nodes
// of nodes 02 and 03

/*

O2     O3
| \	  /	
|  \ /
|	O1
PO

*/
void GRAPH_MEDIC::detach_node_from_references_in_fanin
(
	NODE * node
)
{
	assert(node);

	EDGES::iterator input_edge_iter;
	NODE * source_node;
	EDGE * edge;
	EDGES input_edges = node->get_input_edges();

	// Over all the input nodes find the output node that fans out to 
	// input node and remove a reference to this node
	for (input_edge_iter = input_edges.begin(); 
		input_edge_iter != input_edges.end(); input_edge_iter++) 
	{
		edge = *input_edge_iter;
		assert(edge);

		source_node = edge->get_source_node();

		// remove and delete the edge between the nodes
		m_graph->remove_edge(edge);
		node->remove_input_edge(edge);
		source_node->remove_output_edge(edge);

		delete edge;
	}
}


// Detach node from fanout
//
// PRE: node may have fanout
// POST: the node has nothing that it fans out to
//
// to detach the node 01 we must remove its reference from nodes 02 and O3
// by deleting the input nodes that connect to it.
// this function also decrements the count of the number of edges

/*
   	       PI	
    O1    / 
    /\	 /	
   /  \ /
  O3   O2

*/
void GRAPH_MEDIC::detach_node_from_references_in_fanout
(
	NODE * node
)
{
	assert(node);

	NODE * sink_node;
	EDGES edges = node->get_output_edges();
	EDGE * edge;
	EDGES::iterator edge_iter;


	// Over all the edges that fanout from output node find and delete the input node
	// that the edges connect to. Remove the reference to the input node in its node
	for (edge_iter = edges.begin(); edge_iter != edges.end(); edge_iter++) 
	{
		edge = *edge_iter;
		assert(edge);

		sink_node = edge->get_sink_node();

		// remove and delete the edge between the nodes
		m_graph->remove_edge(edge);
		node->remove_output_edge(edge);
		sink_node->remove_input_edge(edge);

		delete edge;
	}
}




//
// Display warning messages to the user
//
void GRAPH_MEDIC::show_node_deletion_warning
(
	NODE * node
)
{
	assert(node);

	string warning_msg = "Deleting node " + node->get_name();

	if (m_problem_nodes == NODES_NOT_CONNECTED_TO_PO)
	{
		warning_msg += " because it is not reachable from a primary output";
	}
	else if (m_problem_nodes == NODES_NOT_CONNECTED_TO_PI)
	{
		warning_msg += " because it is not reachable from a primary input";
	}
	else
	{
		warning_msg += " because it is isolated in the graph";
	}


	if (! m_shown_eliminate_node_warning)
	{
		Warning(warning_msg);
		m_shown_eliminate_node_warning = true;
	}
	else
	{
		Verbose(warning_msg);
	}
}

//
// Prints stats on the number of nodes deleted
//
void GRAPH_MEDIC::print_node_deletion_stats()
{
	if ( (! g_options->is_no_warn()) && (! g_options->is_quiet()))
	{
		if (m_number_unreachable)
		{
			Warning("Circuit had " << m_number_unreachable << " unreachable nodes");
		}
		if (m_number_comb_deleted)
		{
			Warning("Deleted " << m_number_comb_deleted << 
					" combinational nodes from the graph");
		}
		if (m_number_seq_deleted)
		{
			Warning("Deleted " << m_number_seq_deleted << 
					" sequential nodes from the graph");
		}
	}
}

