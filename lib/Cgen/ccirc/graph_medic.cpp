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

GRAPH_MEDIC::GRAPH_MEDIC(CIRCUIT * circuit, SYMBOL_TABLE * symbol_table)
{
	m_graph					= circuit;
	m_symbol_table			= symbol_table;
	m_problem_nodes			= GRAPH_MEDIC::NODES_NOT_CONNECTED_TO_PO;
	m_number_comb_deleted 	= 0;
	m_number_seq_deleted	= 0;
	m_number_unreachable	= 0;
	m_number_isolated		= 0;

	m_shown_eliminate_node_warning = false;
}

GRAPH_MEDIC::GRAPH_MEDIC(CIRCUIT * circuit)
{
	m_graph					= circuit;
	m_symbol_table			= 0;
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
	m_symbol_table 		= another_graph_medic.m_symbol_table;
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
	m_symbol_table 		= another_graph_medic.m_symbol_table;
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
//       no unconnected nodes or ports exists
//
void GRAPH_MEDIC::delete_unusable_nodes()
{
	// look for unreachable nodes, add them to m_nodes_to_delete
	delete_unreachable_nodes_from_outputs();
	delete_unreachable_nodes_from_inputs();

	// assert: all unreachable nodes are in m_nodes_to_delete
	delete_queued_nodes();

	// delete the unconnected external ports
	delete_unconnected_external_ports();


	// delete any nodes or ports that are not connected to anything in the graph
	delete_isolated_nodes();
	delete_isolated_ports();

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

	while(! m_nodes_to_delete.empty() )
	{
		// pop the first node off the list
		node_to_delete = m_nodes_to_delete.back();
		m_nodes_to_delete.pop_back();

		assert(node_to_delete);

		if (node_to_delete->get_type() == NODE::COMB)
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
// This function might change the name of the output ports if the buffer node feeds a primary output.
// This can be avoided if we rename the source node to the name of the primary output but
// this complicates matters plus I am not sure if it is the right thing to do.
// If a dff fed out to a buffer and the buffer was a PO the dff would be renamed.
//
// Which is better. To rename the outputs or to rename the comb/dff?
//
//
//          0		- output_port_of_source_port: internal and nOutputs = 1
//          |
//			O		- node: buffer node, needs to be removed
//         /|
//		  O	O 		- nodes in fanout:
//
//
//
//			O		- output_port_of_source_port
//         /|
//        O	O		- nodes in fanout 
//
// PRE: buffer/inverter nodes may exist
// POST: all deletable buffer/inverter nodes have been deleted
//
void GRAPH_MEDIC::delete_buffer_and_inverter_nodes()
{

    NODE* node = 0;
	PORT * output_port = 0;
	PORT * output_port_of_node_above = 0;
	NODES & nodes	= m_graph->get_nodes();
	NODES::iterator node_iter;
	PORTS input_ports;

	debugif(DMEDIC, "Looking to see if there are any buffer or inverter nodes.");
	// look at all the nodes to see if they harbor a buffer node
	for (node_iter = nodes.begin(); node_iter != nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);
		input_ports = node->get_input_ports();

		if (have_buffer_node(node,input_ports))
		{
			output_port = node->get_output_port();
			assert(output_port);

			assert(input_ports[0]);
			output_port_of_node_above = input_ports[0]->get_output_port_that_fanout_to_me();
			assert(output_port_of_node_above);

			if (should_delete_buffer_node(output_port, output_port_of_node_above))
			{
				transfer_output_port_connections_to_node_above(output_port, output_port_of_node_above);
				debug("Found inverter-buffer node " << node->get_name() << ".  Marked the node for deletion.");
				m_nodes_to_delete.push_back(node);
			}
			else
			{
				debug("Found inverter-buffer node " << node->get_name() << " that cannot be deleted because the source node output port is not internal or else has fanout");
			}
		}
	}

	if (! m_nodes_to_delete.empty())
	{
		Warning("Buffer/inverter nodes found.  " << m_nodes_to_delete.size() << " nodes will be deleted.");
		delete_queued_nodes();
	}
}

//
//	Takes all of the connections from the output port and transfers them to the source nodes
//	output port
//  
//  PRE: the output port may have connections
//  POST: these connections have been transferred to output_port_of_node above
//
void GRAPH_MEDIC::transfer_output_port_connections_to_node_above
(
	PORT * output_port,
	PORT * output_port_of_node_above
)
{
	assert(output_port && output_port_of_node_above);
	assert(! (output_port_of_node_above->get_external_type() == PORT::PI &&
		   	  output_port->get_external_type() == PORT::PO));

	EDGES::iterator edge_iter;
	EDGE * edge;
	NODE * node_in_fanout;
	string source_node_name = output_port_of_node_above->get_name();
	PORT * node_in_fanout_input_port;

	EDGES output_edges = output_port->get_edges();


	// if the buffer node we want to delete has external connections transfer them to the node above
	if (output_port->get_type() == PORT::EXTERNAL)
	{
		transfer_external_types_to_node_above(output_port, output_port_of_node_above);
	}

	// now transfer all the buffer node's connections to the node above
	for (edge_iter = output_edges.begin(); edge_iter != output_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		node_in_fanout = edge->get_sink_node();
		assert(node_in_fanout);

		// add and create a new input port
		node_in_fanout_input_port = 
			node_in_fanout->create_and_add_port(source_node_name, PORT::INTERNAL,PORT::INPUT, PORT::NONE);

		m_graph->create_edge(output_port_of_node_above, node_in_fanout_input_port, 1);
	}
}

//
// transfers the external type to the source node output port
//
//  PRE: the output port is a PO
//       the output_port_of_node_above is not a PI
//       
//  POST: these external types have been transferred to output_port_of_node above
//
void GRAPH_MEDIC::transfer_external_types_to_node_above
(
	PORT * output_port,
	PORT * output_port_of_node_above
)
{
	assert(output_port && output_port_of_node_above);
	assert(output_port->get_external_type() == PORT::PO);
	assert(output_port_of_node_above->get_external_type() != PORT::PI);

	PORT::EXTERNAL_TYPE external_type;
	
	if (output_port_of_node_above->get_external_type() == PORT::PO)
	{
		assert(output_port_of_node_above->get_external_type() == PORT::PO);
		m_graph->remove_external_port(output_port_of_node_above);
	}

	debugif(DMEDIC,"Transfering external types from " << output_port->get_name() << 
			" to " << output_port_of_node_above->get_name());

	output_port_of_node_above->set_type(PORT::EXTERNAL);
	external_type = output_port->get_external_type();
	output_port_of_node_above->set_external_type(external_type);

	//we don't need to remove the external port because the node will be quequed
	//and the external port will be removed there.
	m_graph->add_external_port(output_port_of_node_above, external_type);
}


// 
// remove any external reference the node may have 
// (external references are PI or PO)
// 
// PRE: the node is valid
//      the node's output port may exist in a list of PO in the graph
// POST: the output port does not exist in a list of PO in the graph
//
void GRAPH_MEDIC::remove_any_external_reference
(
	NODE * node
)
{
	PORT * output_port;
	assert(node);

	output_port = node->get_output_port();
	assert(output_port);

	if (output_port->get_type() == PORT::EXTERNAL)
	{
		// we need to remove it from the PO list
		m_graph->remove_external_port(output_port);
	}
}

// 
// PRE: unconnected primary_inputs, primary_outputs, and global_clocks may exist in the graph
// POST: unconnected primary_inputs, primary_outputs, and global_clocks do not exist in the graph
void GRAPH_MEDIC::delete_unconnected_external_ports()
{
	delete_unconnected_primary_inputs();
	delete_unconnected_primary_outputs();
	delete_unconnected_global_clocks();
}

//
// delete unconnected primary inputs 
//
// PRE: unconnected primary_inputs may exist in the graph
// POST: unconnected primary_inputs do not exist in the graph
void GRAPH_MEDIC::delete_unconnected_primary_inputs()
{
	PORTS PI = m_graph->get_PI();
	PORTS::iterator port_iter;
	EDGES edges;
	PORT * output_port;

	// look for unconnected PI. Unconnected PI will have no edges
	port_iter = PI.begin();
	while (port_iter != PI.end())
	{
		output_port = *port_iter;
		assert(output_port);
		edges = output_port->get_edges();

		debugif(DMEDIC, "Looking at port name " << output_port->get_name());
		if (edges.empty())
		{
			Warning("The PI port " << output_port->get_name() 
					<< " has no connections.  Deleting");
			port_iter = PI.erase(port_iter);

			m_symbol_table->remove_port(output_port->get_name());
			delete output_port;
		}
		else
		{
			m_symbol_table->remove_port(output_port->get_name());
			port_iter++;
		}
	}

	m_graph->set_external_ports(PI, PORT::INPUT);

}


//
// check to make sure that if have a clock port that it is connected to a dff.
// the dff it was connected to might have disappeared
//
// PRE: unconnected global clocks may exist in the graph
// POST: unconnected global clocks do not exist in the graph
//
void GRAPH_MEDIC::delete_unconnected_global_clocks()
{
	PORT * output_port = m_graph->get_global_clock();
	EDGES edges; 

	if (output_port)
	{
		edges = output_port->get_edges();

		debugif(DMEDIC, "Looking at port name " << output_port->get_name());
		if (edges.empty())
		{
			Warning("Clock port " << output_port->get_name() << " has no connections.  Deleting");
			m_graph->erase_global_clock();
			delete output_port;
		}
		// now just remove it from the table whether or not we erased the clock
		m_symbol_table->remove_port(output_port->get_name());
	}
}


//
// Look for and deletes unconnected PO. 
// Unconnected PO will have no node and no edges.
//
// PRE: unconnected primary outputs may exist in the graph
// POST: unconnected primary outputs do not exist in the graph
//
void GRAPH_MEDIC::delete_unconnected_primary_outputs()
{
	PORTS PO = m_graph->get_PO();
	PORTS::iterator port_iter;
	PORT * output_port;
	NODE * node;
	
	port_iter = PO.begin();
	while (port_iter != PO.end())
	{
		output_port = *port_iter;
		assert(output_port);
		node = output_port->get_my_node();

		if (! node)
		{
			Warning("The PO port " << output_port->get_name() 
					<< " has no driving node.  Deleting");

			detach_output_port_from_references_in_fanout(output_port);		
			if (m_symbol_table->query_for_port(output_port->get_name()))
			{
				m_symbol_table->remove_port(output_port->get_name());
			}
			port_iter = PO.erase(port_iter);
			delete (output_port);
		}
		else
		{
			m_symbol_table->remove_port(output_port->get_name());
			port_iter++;
		}
	}
	m_graph->set_external_ports(PO, PORT::OUTPUT);
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
	PORTS::iterator port_iter;
	PORT * output_port;
	NODE * node;
	PORTS PO = m_graph->get_PO();
	
	// mark the nodes up from the primary outputs
	for (port_iter = PO.begin(); port_iter != PO.end(); port_iter++)
	{
		output_port = *port_iter;
		assert(output_port);

		node = output_port->get_my_node();

		if (node && node->get_colour() != colour)
		{
			mark_up_from_node(node, colour);
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
	PORTS	input_ports;
	NODE *			output_node;
	PORTS::iterator port_iter;

	assert(node);
	node->set_colour(colour);

	debugif(DMEDIC, node->get_name() << " = marked " << static_cast<short>(colour) );

	input_ports = node->get_input_ports();

	for (port_iter = input_ports.begin(); port_iter != input_ports.end(); port_iter++) 
	{
		assert(*port_iter);
		output_node = (*port_iter)->get_node_that_fanout_to_me();

		if (output_node && output_node->get_colour() != colour)
		{
			mark_up_from_node(output_node, colour);
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
	PORTS::iterator port_iter;
	PORTS PI = m_graph->get_PI_with_clock();
	PORT * output_port;

	// delete any unmarked node down from the PIs
	for (port_iter = PI.begin(); port_iter != PI.end(); port_iter++) 
	{
		output_port = *port_iter;
		assert(output_port);

		eliminate_down_from_output_port(output_port);
	}
}


//
// Eliminates any UNMARKED node in the graph starting 
// down from this output port
//
// PRE: nodes in the graph may be marked UNMARKED 
// POST: all nodes reachable from this output port down 
//       and marked UNMARKED have been queue for deletion
void GRAPH_MEDIC::eliminate_down_from_output_port
(
	PORT * output_port
)
{
	EDGES edges;	
	EDGES::iterator edge_iter;	
	EDGE * edge;
	PORT * input_port;
	NODE * node;

	assert(output_port);
	edges = output_port->get_edges();
	
	for(edge_iter = edges.begin(); edge_iter != edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		input_port = edge->get_sink();
		assert(input_port);

		node = input_port->get_my_node();
		assert(node);
	
		// if we haven't seen the node, process it
		if (node->get_colour() == NODE::MARKED || 
			node->get_colour() == NODE::UNMARKED)
		{	
			debugif(DMEDIC, "Eliminate down:  searching at " <<  node->get_name());
			if_unmarked_queue_for_deletion(node);

			output_port = node->get_output_port();
			assert(output_port);
			eliminate_down_from_output_port(output_port);
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
	PORTS PI = m_graph->get_PI_with_clock();
	PORTS::iterator port_iter;
	PORT * output_port;

	// mark the nodes down from the PI
	for (port_iter = PI.begin(); port_iter != PI.end(); port_iter++) 
	{
		output_port = *port_iter;

		assert(output_port);
		mark_down_from_output_port(output_port, colour);
	}
}


//
// Colour down from this output port
//
// PRE: node is valid
// POST: all nodes down from this node that are reachable are coloured
void GRAPH_MEDIC::mark_down_from_output_port
(
	PORT * output_port,
	const NODE::COLOUR_TYPE & colour
)

{
	PORT * input_port;
	NODE * node;
	EDGES edges;
	EDGES::iterator edge_iter;


	assert(output_port);
	edges = output_port->get_edges();

	// over all the edges look for more nodes to mark
	for (edge_iter = edges.begin(); edge_iter != edges.end(); edge_iter++) 
	{
		assert(*edge_iter);
		
		input_port = (*edge_iter)->get_sink();
		assert(input_port);

		node = input_port->get_my_node();
		assert(node);	

		// if we haven't visted this node before colour it
		if (node->get_colour() != colour)
		{
			debugif(DMEDIC, node->get_name() << " = marked " 
					<< static_cast<short>(colour) ); 
			node->set_colour(colour);

			output_port = node->get_output_port();
			assert(output_port);
			mark_down_from_output_port(output_port, colour);
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
	PORTS::iterator port_iter;
	PORTS PO = m_graph->get_PO();
	PORT * output_port;
	NODE * node;

	// delete up from the PO
	for (port_iter = PO.begin(); port_iter != PO.end(); port_iter++) 
	{
		output_port = *port_iter;
		assert(output_port);
		node = output_port->get_my_node();
		if (node)
		{
			if_unmarked_queue_for_deletion(node);
			eliminate_up_from_node(node);
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
	PORTS	input_ports;
	PORTS::iterator port_iter;
	EDGE *	edge;
	PORT *	output_port;
	NODE *	output_node;

	for (port_iter = input_ports.begin(); port_iter != input_ports.end(); port_iter++) 
	{
		// input ports should be connected to an edge and the edge should
		// be connected to an output port

		assert(*port_iter);
		edge = (*port_iter)->get_edge();
		assert(edge);
		
		output_port = edge->get_source();
		assert(output_port);

		// if the port is internal then the port should have a node
		if (output_port->get_type() == PORT::INTERNAL)
		{
			output_node = output_port->get_my_node();
			assert(output_node);

			// if we haven't seen it before process it
			if ((output_node->get_colour() == NODE::MARKED) ||
				(output_node->get_colour() == NODE::UNMARKED))
			{
				if_unmarked_queue_for_deletion(node);
				eliminate_up_from_node(output_node);
			}
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
			node->get_colour() == NODE::UNMARKED);

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
	// Visited nodes are removed from the symbol table so
	// we can spot isolated nodes later

	remove_from_symbol_table(node);
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

	PORT * output_port;
	assert(node);
	output_port = node->get_output_port();

	remove_any_external_reference(node);

	m_graph->remove_from_node_list(node);

	// detach the node from both its inputs and outputs
	detach_node_from_references_in_fanin(node);
	detach_output_port_from_references_in_fanout(output_port);


	delete node;
	node = 0;

}

//
// Detach node from fanin
//
// PRE: node may have fanin
// POST: the node has nothing than fans into it
// 
// to detach the node 01 we must remove its reference from the output ports
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

	PORTS::iterator input_port_iter;
	PORT * input_port;
	PORT * output_port;
	EDGE * edge;
	PORTS input_ports = node->get_input_ports();

	// Over all the input ports find the output port that fans out to 
	// input port and remove a reference to this node
	for (input_port_iter = input_ports.begin(); 
		input_port_iter != input_ports.end(); input_port_iter++) 
	{
		input_port = *input_port_iter;
		assert(input_port);

		edge = input_port->get_edge();
		assert(edge);

		output_port = edge->get_source();
		assert(output_port);	

		// remove and delete the edge between the ports
		m_graph->remove_and_delete_edge(output_port, input_port, edge);
	}
}


// Detach node from fanout
//
// PRE: node may have fanout
// POST: the node has nothing that it fans out to
//
// to detach the node 01 we must remove its reference from nodes 02 and O3
// by deleting the input ports that connect to it.
// this function also decrements the count of the number of edges

/*
   	       PI	
    O1    / 
    /\	 /	
   /  \ /
  O3   O2

*/
void GRAPH_MEDIC::detach_output_port_from_references_in_fanout
(
	PORT * output_port
)
{
	assert(output_port);

	PORT * input_port;
	NODE * fanout_node;
	EDGES edges = output_port->get_edges();
	EDGE * edge;
	EDGES::iterator edge_iter;


	// Over all the edges that fanout from output port find and delete the input port
	// that the edges connect to. Remove the reference to the input port in its node
	for (edge_iter = edges.begin(); edge_iter != edges.end(); edge_iter++) 
	{
		edge = *edge_iter;
		assert(edge);

		input_port = edge->get_sink();
		assert(input_port);

		// remove and delete the edge between the ports
		m_graph->remove_and_delete_edge(output_port,input_port,edge);

		fanout_node = input_port->get_my_node();
		assert(fanout_node);	

		fanout_node->remove_input_port(input_port);

		delete input_port;
	}
}

//
// Remove the node and its output port from the symbol table so we can
// spot isolated nodes and output ports.
// Note: We do not need to remove the input ports because
// they were not added to the symbol table and we do not 
// need to remove the edges because they should not be stranded
//
// PRE: the output port is in the symbol table.
// POST: the output port has been removed from the symbol table
//
void GRAPH_MEDIC::remove_from_symbol_table
(
	NODE * node
)
{
	PORT * output_port;

	assert(node);

	m_symbol_table->remove_node(node->get_name());

	output_port = node->get_output_port();
	assert(output_port);

	m_symbol_table->remove_port(output_port->get_name());
		
}

//
// Delete any isolated nodes that might exist in the symbol table
//
// PRE: isolated_nodes may exist
// POST: now they don't.
//
void GRAPH_MEDIC::delete_isolated_nodes()
{
	NODE * node;
	PORT * output_port;

	debugif(DMEDIC, "Looking for isolated nodes");

	m_problem_nodes = GRAPH_MEDIC::ISOLATED_NODES;

	while (! m_symbol_table->empty_of_nodes() )
	{
		node = m_symbol_table->front_node();
		m_symbol_table->pop_front_node();
		output_port = node->get_output_port();
		assert(output_port);
		m_symbol_table->remove_port(output_port->get_name());

		show_node_deletion_warning(node);
		delete_node(node);
		m_number_isolated++;
	}
}

// 
// Delete any isolated ports in the symbol table.
// Ports are deemed isolated because noone has claimed them from the table earlier
//
// PRE: isolated ports may exist
// POST: now they don't.
//
void GRAPH_MEDIC::delete_isolated_ports()
{
	PORT * port;
	debugif(DMEDIC, "Looking for isolated ports");

	// label the problem so that show_port_deletion_warning displays correct warning
	m_problem_nodes = GRAPH_MEDIC::ISOLATED_PORTS;

	while (! m_symbol_table->empty_of_ports() )
	{
		port = m_symbol_table->front_port();
		m_symbol_table->pop_front_port();
		show_port_deletion_warning(port);
		delete port;
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
		warning_msg += " because it is not connected to a primary output";
	}
	else if (m_problem_nodes == NODES_NOT_CONNECTED_TO_PI)
	{
		warning_msg += " because it is not connected to a primary input";
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
// Show warning messages to the user
//
void GRAPH_MEDIC::show_port_deletion_warning
(
	PORT * port
)
{
	assert(port);

	if (! m_shown_eliminate_node_warning)
	{
		Warning("Deleting isolated port " << port->get_name());
		m_shown_eliminate_node_warning = true;
	}
	else
	{
		Verbose("Deleting isolated port " << port->get_name());
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
			Warning("Circuit has " << m_number_unreachable << " unreachable nodes");
		}
		if (m_number_isolated)
		{
			Warning("Circuit has " << m_number_isolated << " isolated nodes");
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

//
//  Do sanity checks, once the graph has been defined, and everything is 
//  supposed to be OK.  
//  PRE: graph has been constructed, all unusable nodes have been deleted
//  POST: We are sure that:
//        a) If we have DFF there is a global clock
//        b) the number of inputs == the number of outputs
//        c) every DFF has just 2 inputs 
//        d) every DFF has a clock input
//        e) every combinational node has at least 1 input
//        f) no input to a combinational node is a clock
//        g) if a node has no fanout it must have a primary output
//
void GRAPH_MEDIC::check_sanity()
{
	PORTS PO = m_graph->get_PO();
	PORTS::iterator port_iter;
	PORT * output_port;
	NODE * node;
	NUM_ELEMENTS number_outputs = 0;
	NUM_ELEMENTS number_inputs  = 0;

	assert(m_graph->get_nNodes() == static_cast<signed>(m_graph->get_nodes().size()));
	assert(m_graph->get_nEdges() == static_cast<signed>(m_graph->get_edges().size()));
	m_graph->colour_graph(NODE::UNMARKED);

	if (m_graph->get_nDFF())
	{
		// if we have DFF we have a global clock
		assert(m_graph->get_global_clock());
		debugif(DMEDIC, "The number of dff = " << m_graph->get_nDFF());
		debugif(DMEDIC, "The number of dff is also = " << m_graph->get_nDFF());
	}

	// get the fan out from the primary inputs
	number_outputs += get_fanout_number_from_PI();

	for (port_iter = PO.begin(); port_iter != PO.end(); port_iter++)
	{
		output_port = *port_iter;
		assert(output_port);

		debugif(DMEDIC, "Checking the sanity off of the PO port " << output_port->get_name());

		node = (*port_iter)->get_my_node();
		assert(node);

		if (node->get_colour() == NODE::UNMARKED)
		{
			check_sanity_up_from_node(node, number_inputs, number_outputs);
		}
	}

	debugif(DMEDIC, "The number of total inputs is " << number_inputs);
	debugif(DMEDIC, "The number of total outputs is " << number_outputs);

	assert(number_inputs == number_outputs);

	// reset the graph
	mark_up(NODE::UNMARKED);
}

//
// Do a sanity check on this node and all nodes up from this node
//
//
//  PRE: node is valid
//  POST: Number_inputs and number_outputs have been updated to reflect the nodes visited.
//
//        We are sure that up from this node that:
//        a) If we have DFF there is a global clock
//        b) every DFF has just 2 inputs 
//        c) every DFF has a clock input
//        d) every combinational node has at least 1 input
//        e) no input to a combinational node is a clock
//        f) if a node has no fanout it must have a primary output
//
void GRAPH_MEDIC::check_sanity_up_from_node
(
	NODE * node,
    NUM_ELEMENTS & number_inputs,
	NUM_ELEMENTS & number_outputs
)
{
	PORT *			output_port;
	PORTS	input_ports;
	PORT *			input_port;
	PORT *			clk_port;
	NODE *			output_node;
	PORTS::iterator port_iter;
	EDGE *			clk_edge;

	assert(node);
	input_ports = node->get_input_ports();
	output_port = node->get_output_port();
	assert(output_port);

	debugif(DMEDIC, "Graph Sanity: name = " << node->get_name() << 
			"\tfanin size = " << input_ports.size() << 
			"\tfanout size = " << output_port->get_edges().size());

	debugif(DMEDIC, "Graph Sanity: name of output port is " << output_port->get_name());

	node->set_colour(NODE::MARKED);

	if (node->get_type() == NODE::SEQ)
	{
		// dff have 2 inputs, the clock and the D-input
		assert(input_ports.size() == 2);

		// the clock input is labeled clock
		clk_port = node->get_clock_port();
		assert(clk_port);
		assert(clk_port->get_io_direction() == PORT::CLOCK);

		// the clock input port is connected to a real global clock
		clk_edge = clk_port->get_edge();
		assert(clk_edge);
		assert(clk_edge->get_source());
		assert(clk_edge->get_source() == m_graph->get_global_clock());
	}
	else
	{
		// The COMB node must have a fanin
		assert(input_ports.size() > 0);


		// although buffers/inverters are ok for VLSI type circuits for FPGA research
		// we will not allow them
		//assert(input_ports.size() > 1);

		// if the fanin is 1 the only fanin should not be the clock
		// because such circuits do not have a defined sequential level
		// and i do not know how to deal with them

		if (input_ports.size() == 1)
		{
			assert(input_ports[0]);
			assert(input_ports[0]->get_io_direction() != PORT::CLOCK);
		}
	}

	// the nodes must have fanout if it is not a primary output
	
	if (output_port->get_type() == PORT::INTERNAL)
	{
		assert(output_port->get_edges().size() > 0);
	}

	number_outputs += output_port->get_edges().size();
	number_inputs  += input_ports.size();

	for (port_iter = input_ports.begin(); port_iter != input_ports.end(); port_iter++) 
	{
		input_port = *port_iter;
		assert(input_port);
		assert(input_port->get_edges().size() == 1);

		if (node->get_type() == NODE::COMB)
		{
			check_input_port_not_connected_to_global_clock(input_port);
		}

		debugif(DMEDIC, "Graph Sanity: Node " << node->get_name() << " has input with name " << input_port->get_name());

		output_node = input_port->get_node_that_fanout_to_me();

		if (output_node && output_node->get_colour() == NODE::UNMARKED)
		{
			check_sanity_up_from_node(output_node, number_inputs, number_outputs);
		}
	}
}

//
// Find the sum of the fanout degree over all of the pi 
//
// RETURNS: the total fanout from all PIs
//
NUM_ELEMENTS GRAPH_MEDIC::get_fanout_number_from_PI()
{
	PORTS PI = m_graph->get_PI_with_clock();
	PORTS::iterator port_iter;
	NUM_ELEMENTS fanout = 0;

	for (port_iter = PI.begin(); port_iter != PI.end(); port_iter++)
	{
		fanout += (*port_iter)->get_edges().size();
	}

	return fanout;
}

//
// We currently do not handle combinational nodes connected to clocks
// 
// PRE: input_port is valid
// POST: nothing if this input port is not connected to the clock
//       if it is, then we have fail and exit the program
//
void GRAPH_MEDIC::check_input_port_not_connected_to_global_clock
(
	PORT * input_port
)
{
	PORT * output_port;
	EDGE * edge;

	assert(input_port);
	edge = input_port->get_edge();
	assert(edge);
	output_port = edge->get_source();
	assert(output_port);

	if (output_port == m_graph->get_global_clock())
	{
		Fail("Combinational node has input port connected to the global clock");
	}
}

//
// Returns: whether or not we have a buffer node
//			 true: if the node is combinational with only 1 input 
//			 false: otherwise
//
bool GRAPH_MEDIC::have_buffer_node
(
	const NODE * node,
	const PORTS & input_ports
)
{
	return (node->get_type() == NODE::COMB && input_ports.size() == 1);
}

//
// Returns: whether or not to delete the buffer node.
//	         false: if the buffer separate a PI from a PO
//	         true : otherwise
// 
bool GRAPH_MEDIC::should_delete_buffer_node
(
	const PORT * output_port,
	const PORT * output_port_of_node_above
)
{
	bool delete_buffer_node = ! (output_port->get_external_type() == PORT::PO && 
								 output_port_of_node_above->get_external_type() == PORT::PI);

	return delete_buffer_node;
}
