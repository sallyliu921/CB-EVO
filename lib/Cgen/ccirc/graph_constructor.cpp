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




#include "graph_constructor.h"

// from the parser parse_blif.y
extern 		int yyparse();
extern 		int yyerror(char * error_text);
extern		long g_linenum;

const string EDGE_CONNECTION_TEXT = "_TO_";	
const string EDGE_SUFFIX			= "_EDGE";	

/*
 * What the graph constructor does:
 * 
 * 1. Add external PI as output ports
 * 2. Add external PO as output ports
 * 3. Add combinational or latch nodes nodes.
 * 		3.1 Add the node
 * 		3.2 Add the output port. Create it if it doesn't exist
 * 		3.3 Add the input port. Connect to the output port. Create output
 * 		    port if it doesn't exist
 * 3. Finds any global clock
 * 4. Uses the graph_medic to delete unused nodes
 */

GRAPH_CONSTRUCTOR::GRAPH_CONSTRUCTOR(const string & circuit_name, OPTIONS * options)
{
	m_options 				= options;
	m_symbol_table 			= new SYMBOL_TABLE;
	m_graph					= new CIRCUIT;	

	assert(m_symbol_table);
	assert(m_graph);
	assert(m_options);

	m_graph->set_name(circuit_name);
}


GRAPH_CONSTRUCTOR::GRAPH_CONSTRUCTOR(OPTIONS * options)
{
	m_options 				= options;
	m_symbol_table 			= new SYMBOL_TABLE;
	m_graph					= new CIRCUIT;	

	assert(m_symbol_table);
	assert(m_graph);
	assert(m_options);

	m_graph->set_name(m_options->get_circuit_name());
	
}

GRAPH_CONSTRUCTOR::GRAPH_CONSTRUCTOR(const GRAPH_CONSTRUCTOR & another_graph_constructor)
{
	m_options				= another_graph_constructor.m_options;
	m_graph					= another_graph_constructor.m_graph; 
	m_symbol_table 			= another_graph_constructor.m_symbol_table;
}
	
GRAPH_CONSTRUCTOR & GRAPH_CONSTRUCTOR::operator=(const GRAPH_CONSTRUCTOR & another_graph_constructor)
{
	m_options				= another_graph_constructor.m_options;
	m_graph					= another_graph_constructor.m_graph; 
	m_symbol_table 			= another_graph_constructor.m_symbol_table;

	return (*this);
}

GRAPH_CONSTRUCTOR::~GRAPH_CONSTRUCTOR()
{
	// do not delete the graph constructed, just delete the symbol table

	debugif(DCONST, "About to delete the graph constructor");
	delete m_symbol_table;
}


//
// The function adds the external ports to the graph.
// External ports are primary input/output ports and 
// the ghost input/output ports.
// 
// PRE: port_name is the name of the port
//      external_type is its type
// POST: the port has been added to the graph
//       the port has been added to the symbol table
// 
//
//
void GRAPH_CONSTRUCTOR::new_external_port
(
	const string & port_name,
	const PORT::EXTERNAL_TYPE & external_type
)
{
    PORT * port;

    port = m_symbol_table->query_for_port(port_name);
    if (!port) 
	{
        port = m_graph->create_and_add_external_port(port_name, external_type);
		m_symbol_table->insert_port(port_name, port);
    } 
	else 
	{
		Warning("Graph Constructor:  Duplicate external port '" << port_name 
				<< "' ignored");
    }
}

//
// New combination logic (lut) block
// 
// PRE: variable_name_stack is the names of the inputs and outputs
//      current_lut contains the contents of the lut
// POST: If the lut is not a constant function a combinational node and its input and output 
//       ports have been created in the graph
//       If we are storing luts, the lut has been stored
//       The node and output port have been added to the symbol table
//
void GRAPH_CONSTRUCTOR::new_combination_block
(
	VARIABLE_STACK_TYPE * variable_name_stack,
	LUT * current_lut
)
{
    string node_name; 
	NODE * node;
	PORTS input_ports;

    debugif(DCONST, "---------------------------------------------------");
    debugif(DCONST, "Creating a new combination block");

    assert(! variable_name_stack->empty());

    /*  If only one name is on the stack the node has no inputs:
		warn and return.
		It is a legal blif format but is unsupported.
	 	(The single name means that the node is either a constant 1 or 0.
		(The truth table below it should either be a single 1 or have no 
		entries to signify a constant 0 function))
     */

    if (variable_name_stack->size() == 1) 
	{
		node_name = variable_name_stack->back();
		variable_name_stack->pop_back();
		Warning("Graph Constructor: disconnected node " << node_name << " ignored");
		delete current_lut;
		current_lut = 0;
		return;
    }


    //  Last name on the stack is the output name and also the node name
	node_name = variable_name_stack->back();
	variable_name_stack->pop_back();

	node = m_graph->create_node(node_name, NODE::COMB);	
	m_symbol_table->insert_node(node_name, node);

	add_output_port(node, node_name);
	add_input_ports_and_connect_to_graph(node, variable_name_stack);

    debugif(DCONST, "Done creating new combinational block");
    debugif(DCONST, "---------------------------------------------------");

    if (m_options->is_store_luts())
	{
		node->set_lut(current_lut);
    }
	else
	{
    	delete current_lut;
    	current_lut = 0;
	}
}

//
// Add an output port for the node to the graph
// If the port already exists then the port should either be a primary output 
// or it should have been prev. created by an input port that it feeds 
//
// PRE: node is valid
//      node_name contains the name of the node
// POST: if the output port of the node exists it was added to the node
//       else an output port was created and added to the node and the symbol table
//
PORT * GRAPH_CONSTRUCTOR::add_output_port
(
	NODE *	node,
	const string & node_name
)
{
	PORT * output_port = 0;

	output_port = m_symbol_table->query_for_port(node_name);

	// check if the port exists if not create a new port
	if (! output_port)
	{
		debugif(DCONST, "Creating new output port for '" << node_name << "'");

		output_port = node->create_and_add_port(node_name, PORT::INTERNAL, 
												PORT::OUTPUT, PORT::NONE);
		m_symbol_table->insert_port(node_name, output_port);
	}
	else
	{
		debugif(DCONST,"Found existing output port for '" << node_name << "'");
		node->add_port(output_port);
		output_port->set_my_node(node);
	}

	return output_port;
}


// 	Creates the input ports and connects them by edges to 
// 	the output ports that feed them.
// 	Creates the output ports if necessary.
//
// PRE: node is valid
// 		The names on the stack are the input node names
// POST: the variable_name_stack is empty
//       the input ports have been created and added to the node and graph
//       edges have been created to connect the input port to output ports
//       the output ports were either found in the symbol table or created
// 	

void GRAPH_CONSTRUCTOR::add_input_ports_and_connect_to_graph
(
	NODE * node,
	VARIABLE_STACK_TYPE	* variable_name_stack
)
{
	string node_name = "";
	string input_port_name = "";
	string output_node_name = "";

	assert(node);
	node_name = node->get_name();
	
    //  Remaining stacked elements are the inputs to the node.   First
    //  make sure there are the right number of them.  Too few is ok.
    if (m_options->get_k() != 0) {
	if (static_cast<K_TYPE>(variable_name_stack->size()) > m_options->get_k()) 
	{
	    Error("Too many arguments " << variable_name_stack->size() 
				<< " for k= " << g_options->get_k());
	    yyerror("");
	}
    }

	// while input variable names exist find or create input ports
	while(! variable_name_stack->empty() )
	{
		output_node_name = variable_name_stack->front();
		variable_name_stack->pop_front();

		add_an_input_port_and_connect_to_graph(node, node_name, output_node_name);
	}	

}

//
// add an input port to the node and try to connected the inport to the graph
// if its desired output port does not exist, create it.
// 
// PRE: node is valid
// 		node_name is the name of the node
// 		output_node_name the name of the node that feeds into node_name
//
// POST: the input port have been created and added to the node and graph
//       an edge has been created to connect the input port to an output port
//       the output port was either found in the symbol table or created
PORT * GRAPH_CONSTRUCTOR::add_an_input_port_and_connect_to_graph
(
 	NODE * node,
	const string & node_name,
	const string & output_node_name
)
{

	PORT * input_port = 0;

    debugif(DCONST, "Adding an input port with name : " << output_node_name);

	input_port = node->create_and_add_port(output_node_name, PORT::INTERNAL,
											PORT::INPUT, PORT::NONE);
	connect_input_port_to_graph(input_port, node_name, output_node_name);

	return input_port;
}

// PRE: node is valid
//      node_name is the name of the node
//      clock_name is the name of the clock
// POST: a input port has been added to the node 
//       the input port has been connected by an edge to the global clock 
//       the global clock was or is now recognized as the global clock
//
//       note: if a clock is internal, the program fails and exits
void GRAPH_CONSTRUCTOR::add_clock_port_and_connect_to_clock
(
 	NODE * node,
	const string & node_name,
	const string & clock_name
)
{

	PORT * clock_port = 0;
	PORT * global_clock_port = 0;
	string edge_name;

    debugif(DCONST, "Adding an clock port with name : " << clock_name);

	clock_port = node->create_and_add_port(clock_name, 
											PORT::INTERNAL,
											PORT::CLOCK, PORT::NONE);


	debugif(DCONST,"Looking for global clock with name  " << clock_name);
	global_clock_port  = m_symbol_table->query_for_port(clock_name);

	if (global_clock_port == 0 || 
			global_clock_port->get_type() == PORT::INTERNAL)
	{
		// The clock is internal but it should have been created at 
		// the model statement line. It should have been listed at the .input
		// or maybe even .ginput line. Internal clocks are not supported.
		
		Fail("Did not find clock port or found internal clock port");
	}

	debugif(DCONST,"Found the global clock port");

	if (global_clock_port->get_io_direction() != PORT::CLOCK)
	{
		// the clock port is listed as a primary input.  
		// tell the circuit that we found the global clock port
		m_graph->set_global_clock(global_clock_port);
	}
	assert(m_graph->get_global_clock() == global_clock_port);

	add_edge_between_ports(node_name, clock_port, global_clock_port);
}
//
//	This function will hook the input port up to the output port that feeds it
//	with an edge.
//	If an output port does not exist, it will create one.
//
//	PRE: input port is valid
//	     input_node_name is the name of the node the input port is connected to
//	     output_port_name is the name of the node that fanin into the input_port
//	PORT: an edge is created between the input port and an output port
//	      if the output_port did not exist in the symbol table the output port
//	      is created and added to the symbol table.
//
void GRAPH_CONSTRUCTOR::connect_input_port_to_graph
(
	PORT * input_port,
	const string & input_node_name,
	const string & output_port_name
)
{	
	assert(input_port);

	string input_port_name = input_port->get_name();

	PORT * output_port 		= m_symbol_table->query_for_port(output_port_name);

	debugif(DCONST, "Looking for an output port '" << output_port_name << "'" <<
			" that feeds the input port");
	
	// Create an internal output port if none exists 
	// Note : 	If the output ports is external then it should already 
	// 			have been created

	if (!output_port)
	{
		debugif(DCONST, "No match has been found for output port. Creating one");
		output_port = new PORT(output_port_name, PORT::INTERNAL, PORT::OUTPUT, PORT::NONE);	
		assert(output_port);
		m_symbol_table->insert_port(output_port_name, output_port);
	}
	else
	{
		debugif(DCONST, "An output port match has been found for the input port");
		debugif((output_port->get_type() != PORT::INTERNAL && DCONST), "An external port match was found");
		assert(output_port->get_io_direction() == PORT::OUTPUT);
	}

	add_edge_between_ports(input_node_name, input_port, output_port);

}

//
//	Adds an edge between the ports
//
//	PRE: input_port and output_port are valid
//	PORT: an edge is created between the input port and an output port
//
void GRAPH_CONSTRUCTOR::add_edge_between_ports
(
	const string & input_node_name,
	PORT * input_port,
	PORT * output_port
)
{
	string edge_name = "";
	EDGE * edge	= 0;

	assert(input_port);
	assert(output_port);

	string output_port_name = output_port->get_name();

	// check to see if edge is duplicate
	edge_name = get_edge_name(input_node_name, output_port_name);
	edge = m_symbol_table->query_for_edge(edge_name);

	assert(Dlook_at); // i am not adding the edge to the symbol table. why?
					
	if (! edge)
	{
		debugif(DCONST, "Adding an edge from " << output_port_name << " to " 
				<< input_node_name << " with name " << edge_name);	
		edge = m_graph->create_edge(output_port, input_port, 1);
		assert(edge);
	}
	else
	{
		Warning("Duplicate edge " << edge_name <<  " ignored near/above line " 
				<<  g_linenum << " of input.");
	}

}

//
// 
// PRE: cube is an input line of the function
//      output_value is the output value of the line (is the function SOP or POS)
//      number_input_variables is the number of inputs
//      current_lut is valid if we have previously created the lut, otherwise it is NULL
// POST: if we are not storing luts, nothing
//       of we are storing luts then lut was created if it didn't exists and
//       this line of the function was added to the lut
//
//       This line of the function was expanded if that option was selected.
//       Errors are generated if there were too many inputs
//
void GRAPH_CONSTRUCTOR::new_truth_table_entry
(
	string & cube, 
	VALUE_TYPE output_value, 
	short number_input_variables,
	LUT * current_lut
)
{

	// don't do anything with the output_value right now

    debugif(DCONST, "Graph Const: new truth table entry: called with cube: " 
			<< cube << " value: " << output_value 
			<< " and number vars: " << number_input_variables);

    //  Basic error-checking: cube-size correct?
    assert(output_value == 0 || output_value == 1);

    if (static_cast<short>(cube.size()) != number_input_variables) 
	{
		cout << "cube size " << cube.size() << " #inputs " << number_input_variables << endl;
		yyerror("Wrong number of bits in this table entry");
    }

    //  If we're not storing the LUTs, we're done.  Free the space and return.

    if (m_options->is_store_luts()) 
	{
		debugif(DCONST, "Proceeding to store this cube");
    } 
	else 
	{
		debugif(DCONST, "LUTs not stored, returning");
		cube = "";
		return;
    }

    if (number_input_variables > m_options->get_k()) 
	{
		yyerror("Too many bits in this table entry");
    } 
	else if (m_options->is_expand_luts() && 
			number_input_variables < m_options->get_k()) 
	{
		//  Expand the lut with don't care's
		Verbose("Converting undersized cube '" << cube); 
		cube.resize(m_options->get_k()+1, '-');
		Verbose("to full-size '" << cube <<"'");
		number_input_variables = m_options->get_k();
    }

	if (! current_lut)
	{
		current_lut = new LUT;
		debugif(DCONST,"Start of a new lut");
	}

	assert(current_lut);
	current_lut->add_cube(cube, output_value);
}

//  
//  Set _value from the string coming in as the value of a cube.
//
//  PRE: the value of the function has not been previously set.
//  POST: the function is specified as either a sum-of-products (value=1) or 
//  products of sum (value=0)
//
//  If we have previously specified the output value an error is generated or if 
//  the value text is something other than 0 or 1
//
VALUE_TYPE GRAPH_CONSTRUCTOR::new_value
(
	const string & value_text,
	LUT * current_lut
)
{
	VALUE_TYPE value = -1;

    debugif(DCONST, "new_value: passed " << value_text);

    if (value_text == "1")
	{
    	// if already have a cube then all cubes should be in sum of products form
        if ( current_lut && current_lut->is_sum_of_products()) 
		{
        	yyerror("Both max/minterm specified for this function.");
        } 
		else
		{
			value = 1;
		}
    } 
	else if (value_text ==  "0")
	{
    	// if already have a cube then all cubes should be in product of sums form
        if (current_lut && (! current_lut->is_sum_of_products()) )
		{
        	yyerror("Both max/minterm specified for this function.");
        } 
		else
		{
			value = 0;
		}
    } 
	else 
	{
        yyerror("Illegal cover_value, must be 0 or 1");
    }
    
	debugif(DCONST, "new_value: done, returning " << value);

    return value;
}
//
//  Create a new flip_flop, much the way that we did the logic.  
//  The clock is ignored, and instead all flip_flops are connected to the
//  global clock.
//
//  PRE: input_port_name contains the name of the D-input
//       output_port_name contains the name of the Q-output/node
//       clk_name contains the name of the clock
//  POST: the flip-flop and its input and output ports have been created and added to the graph
//        the edges that connected the D-input to its output have been created and added to the graph
//        the node and output port have been added to the symbol table
//
void GRAPH_CONSTRUCTOR::new_flip_flop
(
	string & input_port_name, 
	string & output_port_name, 
	string & clk_name
)
{
	NODE * node;
	const string node_name = output_port_name;
	const string output_to_node_name = input_port_name;

    debugif(DCONST, "---------------------------------------------------\n");
	debugif(DCONST, "Creating a sequential node with name : " << node_name);

	node = m_graph->create_node(node_name, NODE::SEQ);	
	m_symbol_table->insert_node(node_name, node);
		
	add_output_port(node, node_name);

	// create a clock and input port and connect them to the global clock
	add_clock_port_and_connect_to_clock(node, node_name, clk_name);
	add_an_input_port_and_connect_to_graph(node, node_name, output_to_node_name);
}

// 
//
// Delete any nodes that are unreachable for the primary inputs and outputs,
// that are unconnected or osolated, and that are buffer/inverter nodes.
// 
// PRE: nodes might exist in the graph that are unusable.
// POST: these nodes no longer exist
//
void GRAPH_CONSTRUCTOR::delete_unusable_nodes()
{
	GRAPH_MEDIC graph_medic(m_graph, m_symbol_table);
	
	graph_medic.delete_unusable_nodes();
	graph_medic.delete_buffer_and_inverter_nodes();
	graph_medic.check_sanity();
}


// RETURNS: a name for the edge
string GRAPH_CONSTRUCTOR::get_edge_name
(
	const string & input_node_name,
	const string & output_port_name
)
{
	assert(! input_node_name.empty());
	assert(! output_port_name.empty());

	return output_port_name + EDGE_CONNECTION_TEXT + input_node_name + EDGE_SUFFIX;
}
