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





#include "drawer.h"
#include <algorithm>
#include "util.h"

const string IO_STYLE  		= " [shape=diamond,style=filled,colour=darkslategray4]";
const string FAKE_IO_STYLE 	= " [shape=diamond,style=dotted]";
const string DFF_STYLE 		= " [shape=box,style=filled,colour=blue]";
const string TO_DFF_EDGE	= " [style=dotted]";

/*
#define valid(node)	(!node->invis && \
			    (node->type == PI || node->type == LOG || \
			     node->type == PO || node->type == DFF )) 
*/

DRAWER::DRAWER()
{
	m_display_illegal_char_msg = true;
}
DRAWER::DRAWER(const DRAWER & another_drawer)
{
	assert(false);
}

DRAWER & DRAWER::operator=(const DRAWER & another_drawer)
{
	assert(false);
	return (*this);
}

DRAWER::~DRAWER()
{
}



// 
// add the preamble to the dot output file
//
void DRAWER::preamble()
{
    m_output_file << "digraph G {";
    m_output_file << "size=\"10,7\";"; 			// the size of the image is 10,7 inches
    m_output_file << "orientation=landscape;";
    m_output_file << "concentrators=true;";
    m_output_file << "center=on;";
    m_output_file << "ratio=fill;";
	m_output_file << endl; 
}

// 
// add the ending statement to the dot file
//
void DRAWER::cleanup()
{
    m_output_file << "}" << endl;
}

//
// Draws a primary input
//
// PRE: draw_input is valid. m_output_file is open
// POST: a primary input statement has been outputed
//
void DRAWER::draw_input
(
	PORT * input_port
)
{
	NODE * my_node;
	assert(input_port);
	my_node = input_port->get_my_node();
	assert(my_node);

	m_output_file << "\t" << print_name(input_port->get_name()) << "_IN" << " -> " << 
		print_name(my_node->get_name()) << endl;

	// for now just output the IO_STYLE
	m_output_file << "\t" << print_name(input_port->get_name()) << "_IN" << IO_STYLE;
	m_output_file << endl;
}


// 
// draws the edges that fanin to the node
//
// PRE: node is valid. m_output_file is open
// POST: edge statements of the fanin to the node have been outputed
//
void DRAWER::draw_node_fanin_edges
(
	NODE * sink_node
)
{
	assert(sink_node);

	PORTS	input_ports;
	PORTS::iterator port_iter;
	NODE *	source_node = 0;
	PORT * input_port = 0;

	draw_node(sink_node);

	input_ports = sink_node->get_input_ports();

	for (port_iter = input_ports.begin(); port_iter != input_ports.end(); port_iter++) 
	{
		input_port = *port_iter;
		assert(input_port);


		source_node = input_port->get_node_that_fanout_to_me();

		if (source_node)
		{
			draw_edge(source_node, sink_node);
		}
		else if (input_port->get_output_port_that_fanout_to_me()->get_io_direction() != PORT::CLOCK)
		{
			// we have a primary input
			draw_input(input_port);
		}
	}
}

// 
// draws an node
//
// PRE: node is valid. m_output_file is open
// POST: an node statement has been output to the file
//
void DRAWER::draw_node
(
	NODE * node
)
{
	assert(node);
	assert(node->get_type() == NODE::COMB || node->get_type() == NODE::SEQ);
	
	// if the node is a dff defined the graph shape as a box
	if (node->get_type() == NODE::SEQ)
	{
		m_output_file << "\t" << print_name(node->get_name()) << DFF_STYLE;
	}
}
		
	
// 
// draws an edge
//
// PRE: source_node and sink node are valid. m_output_file is open
// POST: an edge statement has been output to the file
//
void DRAWER::draw_edge
(
	NODE * source_node,
	NODE * sink_node
)
{
	assert(source_node);
	assert(sink_node);

	m_output_file << "\t" << print_name(source_node->get_name()) << " -> " 
		<< print_name(sink_node->get_name());

	if (sink_node->get_type() == NODE::SEQ)
	{
		m_output_file << TO_DFF_EDGE;
	}

	m_output_file << ";" << endl;
}

//
// Draws the circuit
//
// PRE: Circuit is valid.
//      Delay levels have been defined
// POST: circuit_name.dot file has been created that contains a drawing
//       in dot format
//
void DRAWER::draw_full_graph
(
	CIRCUIT * circuit
)
{
	assert(circuit);

	string file_name = circuit->get_name() + ".dot";
	NODES& nodes = circuit->get_nodes();
	NODES::iterator node_iter;
	NODE * node = 0;

	m_output_file.open(file_name.c_str(), ios::out);

	if (! m_output_file.is_open())
	{
		Warning("Could not open output file " << file_name << ". Returning");
		return;
	}

	preamble();

	constrain_io_ranks(circuit);

	for (node_iter = nodes.begin(); node_iter != nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);
		draw_node_fanin_edges(node);
	}

	
	//draw_inputs();

	cleanup();

	debugif(DCLUSTER, "Finished drawing this graph");
	m_output_file.close();
	
}

//
// Print a name of a node with "[" and "]" removed because they cause 
// problems with dot
//
// PRE: node_name is the name of the node
// RETURN: node_name with illegal characters stripped out
//
//
string DRAWER::print_name
(
	string node_name
)
{
	unsigned int illegal_character_pos;
	string replacement_char = "_";
	

	illegal_character_pos = node_name.find('[',0);
	while (illegal_character_pos < node_name.size())
	{
		if (m_display_illegal_char_msg)
		{
			debug("Found a [ in the nodes names while outputing to dot format." <<
				  " Removing [ in all node nodes.");
			m_display_illegal_char_msg = false;
		}

		node_name.replace(illegal_character_pos,1,replacement_char);
		illegal_character_pos = node_name.find('[',0);
	}

	illegal_character_pos = node_name.find(']',0);
	while (illegal_character_pos < node_name.size())
	{
		if (m_display_illegal_char_msg)
		{
			debug("Found a ] in the nodes names while outputing to dot format." <<
				  " Removing ] in all node nodes.");
			m_display_illegal_char_msg = false;
		}
		node_name.replace(illegal_character_pos,1,replacement_char);
		illegal_character_pos = node_name.find(']',0);
	}

	return node_name;
}

//
// Make our picture array the nodes by delay level from top (0th delay level) to 
// bottom (dmax delay level)
//
// PRE: Circuit is valid
// POST: dot will draw our picture arrayed by delay level
//
void DRAWER::constrain_io_ranks
(
	CIRCUIT * circuit
)
{
	assert(circuit);
	SEQUENTIAL_LEVEL * seq_level = circuit->get_sequential_level(); 
	assert(seq_level);
	DELAY_LEVELS& delay_levels = seq_level->get_delay_levels();
	DELAY_TYPE delay,
			   max_delay = circuit->get_maximum_combinational_delay();
	PORTS PIs = circuit->get_PI();
	PORTS::const_iterator port_iter;
	PORT * port = 0;
	DELAY_LEVEL::const_iterator node_iter;

	NODE * node = 0;

	// 1st the dff.
	assert(max_delay > 0);
	m_output_file << "{rank=min;  /* delay: " << 0 << " */" << endl;
	for (node_iter = delay_levels[0].begin(); node_iter != delay_levels[0].end(); node_iter++)
	{
		node = *node_iter;
		assert(node);
		m_output_file << print_name(node->get_name()) << ";\n";
	}

	// 2nd the pi
	for (port_iter = PIs.begin(); port_iter != PIs.end(); port_iter++)
	{
		port = *port_iter;
		assert(port);
		m_output_file << print_name(port->get_name()) << "_IN" <<  ";\n";
	}
	m_output_file << "}\n";


	// 3rd output the rest of the names
	for (delay = 1; delay < max_delay; delay++)
	{

		m_output_file << "{rank=same;  /* delay: " << delay << " */" << endl;
		for (node_iter = delay_levels[delay].begin(); node_iter != delay_levels[delay].end(); node_iter++)
		{
			node = *node_iter;
			assert(node);
			m_output_file << print_name(node->get_name()) << ";\n";
		}
		m_output_file << "}\n";
	}


	m_output_file << "{rank=max;  /* delay: " << max_delay << " */" << endl;
	for (node_iter = delay_levels[max_delay].begin(); node_iter != delay_levels[max_delay].end(); node_iter++)
	{
		node = *node_iter;
		assert(node);
		m_output_file << print_name(node->get_name()) << ";\n";
	}
	m_output_file << "}\n";
}
