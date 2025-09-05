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



#include "gen_drawer.h"
#include <algorithm>
using namespace std;


const string IO_STYLE  		= " [shape=diamond,style=filled,colour=darkslategray4]";
const string FAKE_IO_STYLE 	= " [shape=diamond,style=dotted]";
const string DFF_STYLE 		= " [shape=box,style=filled,colour=blue]";
const string TO_DFF_EDGE	= " [style=dotted]";

GEN_DRAWER::GEN_DRAWER()
{
	m_circuit = 0;
}
GEN_DRAWER::GEN_DRAWER(const GEN_DRAWER & another_drawer)
{
}

GEN_DRAWER & GEN_DRAWER::operator=(const GEN_DRAWER & another_drawer)
{

	return (*this);
}

GEN_DRAWER::~GEN_DRAWER()
{
}


// 
// add the preamble to the dot output file
//
void GEN_DRAWER::preamble()
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
void GEN_DRAWER::cleanup()
{
    m_output_file << "}" << endl;
}

//
// Draws a primary input
//
// PRE: edge is valid. m_output_file is open
// POST: a primary input statement has been outputed
//
void GEN_DRAWER::draw_input
(
	EDGE * edge
)
{
	assert(edge);

	NODE * sink_node = edge->get_sink_node();
	NODE * source_node = edge->get_source_node();

	assert(sink_node && source_node);

	m_output_file << "\t" << print_name(source_node->get_name()) << "_IN" << " -> " << 
		print_name(sink_node->get_name()) << endl;

	// for now just output the IO_STYLE
	m_output_file << "\t" << print_name(source_node->get_name()) << "_IN" << IO_STYLE;
	m_output_file << endl;
}


// 
// draws the edges that fanin to the node
//
// PRE: sink_node is valid. m_output_file is open
// POST: edge statements of the fanin to the node have been outputed
//
void GEN_DRAWER::draw_node_fanin_edges
(
	NODE * sink_node
)
{
	assert(sink_node);

	EDGES input_edges;
	EDGES::iterator edge_iter;
	NODE *	source_node = 0;
	EDGE * input_edge = 0;

	draw_node(sink_node);

	input_edges = sink_node->get_input_edges();

	for (edge_iter = input_edges.begin(); edge_iter != input_edges.end(); edge_iter++) 
	{
		input_edge = *edge_iter;
		assert(input_edge);

		source_node = input_edge->get_source_node();

		if (source_node->is_PI())
		{
			draw_input(input_edge);
		}
		else
		{
			draw_edge(source_node, sink_node);
		}
	}
}

// 
// draws an node
//
// PRE: node is valid. m_output_file is open
// POST: an node statement has been output to the file
//
void GEN_DRAWER::draw_node
(
	NODE * node
)
{
	assert(node);
	
	// if the node is a dff defined the graph shape as a box
	if (node->is_DFF())
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
void GEN_DRAWER::draw_edge
(
	NODE * source_node,
	NODE * sink_node
)
{
	assert(source_node);
	assert(sink_node);

	m_output_file << "\t" << print_name(source_node->get_name()) << " -> " 
		<< print_name(sink_node->get_name());

	if (sink_node->is_DFF())
	{
		m_output_file << TO_DFF_EDGE;
	}

	m_output_file << ";" << endl;
}


//
// Draws the circuit
//
// PRE: circuit is valid.
// POST: circuit_name.dot file has been created that contains a drawing
//       in dot format
//
void GEN_DRAWER::draw_full_graph
(
	CIRCUIT_GRAPH * circuit
)
{
	assert(circuit);
	m_circuit = circuit;

	string file_name = circuit->get_name() + ".dot";
	NODES& nodes = circuit->get_nodes();
	NODES::iterator node_iter;
	NODE * node = 0;

	m_output_file.open(file_name.c_str(), ios::out);

	preamble();

	constrain_io_ranks();

	for (node_iter = nodes.begin(); node_iter != nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);
		draw_node_fanin_edges(node);
	}

	
	//draw_inputs();

	cleanup();

	debug("Finished drawing this graph");
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
string GEN_DRAWER::print_name
(
	string node_name
)
{
	unsigned int illegal_character_pos;
	string replacement_char = "_";


	illegal_character_pos = node_name.find('[',0);
	while (illegal_character_pos < node_name.size())
	{
		debug("Found a [ in the nodes names while outputing to dot format.  Removing for the dot file only.");

		node_name.replace(illegal_character_pos,1,replacement_char);
		illegal_character_pos = node_name.find('[',0);
	}

	illegal_character_pos = node_name.find(']',0);
	while (illegal_character_pos < node_name.size())
	{
		debug("Found a ] in the nodes names while outputing to dot format.  Removing for the dot file only.");
		node_name.replace(illegal_character_pos,1,replacement_char);
		illegal_character_pos = node_name.find(']',0);
	}

	return node_name;
}


//
// Make our picture array the nodes by delay level from top (0th delay level) to 
// bottom (dmax delay level)
//
// PRE: m_circuit is valid
// POST: dot will draw our picture arrayed by delay level
//
void GEN_DRAWER::constrain_io_ranks()
{
	DELAY_TYPE delay,
			   max_delay = m_circuit->get_delay();
	NODES& nodes = m_circuit->get_nodes();
	NODES::const_iterator node_iter;

	NODE * node = 0;

	// 1st the dff and po
	assert(max_delay > 0);
	m_output_file << "{rank=min;  /* delay: " << 0 << " */" << endl;
	for (node_iter = nodes.begin(); node_iter != nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		if (node->get_delay_level() == 0)
		{
			if (node->is_PI())
			{
				m_output_file << print_name(node->get_name()) << "_IN" <<  ";\n";
			}
			else 
			{
				assert(node->is_DFF());
				m_output_file << print_name(node->get_name()) << ";\n";
			}
		}
	}
	m_output_file << "}\n";

	// 2nd output the rest of the names
	for (delay = 1; delay < max_delay; delay++)
	{
		m_output_file << "{rank=same;  /* delay: " << delay << " */" << endl;
		for (node_iter = nodes.begin(); node_iter != nodes.end(); node_iter++)
		{
			node = *node_iter;
			assert(node);

			if (node->get_delay_level() == delay)
			{
				m_output_file << print_name(node->get_name()) << ";\n";
			}
		}
		m_output_file << "}\n";
	}


	m_output_file << "{rank=max;  /* delay: " << max_delay << " */" << endl;
	for (node_iter = nodes.begin(); node_iter != nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);
		if (node->get_delay_level() == max_delay)
		{
			m_output_file << print_name(node->get_name()) << ";\n";
		}
	}
	m_output_file << "}\n";
}
