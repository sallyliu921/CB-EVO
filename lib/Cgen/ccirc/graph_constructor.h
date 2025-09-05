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



#ifndef graph_constructor_H
#define graph_constructor_H

#include "circ.h"
#include "circuit.h"
#include "symbol_table.h"
#include <deque>
#include "lut.h"
#include "graph_medic.h"

//
// Class_name GRAPH_CONSTRUCTOR
//
// Description
// 		Builds the graph
// 		Works with the parser parse_blif.y
//

typedef deque<string> VARIABLE_STACK_TYPE;

class GRAPH_CONSTRUCTOR
{
public:
	GRAPH_CONSTRUCTOR(OPTIONS * options);
	GRAPH_CONSTRUCTOR(const string & circuit_name, OPTIONS * options);
	GRAPH_CONSTRUCTOR(const GRAPH_CONSTRUCTOR & another_graph_constructor);
	GRAPH_CONSTRUCTOR & operator=(const GRAPH_CONSTRUCTOR & another_graph_constructor);
	~GRAPH_CONSTRUCTOR();

	void new_external_port(const string & port_name,
							const PORT::EXTERNAL_TYPE & external_type);
	void new_combination_block(VARIABLE_STACK_TYPE * variable_name_stack,
								LUT * current_lut);
	void new_flip_flop(	string & input_port_name, string & output_port_name, 
					string & clk_name);
	void new_truth_table_entry(string & cube, VALUE_TYPE output_value, 
							   short number_input_variables, LUT * current_lut);
	VALUE_TYPE new_value(const string & value_text, LUT * current_lut);

	CIRCUIT *	get_constructed_graph() { return m_graph;}

	void 		delete_unusable_nodes();
private:
	CIRCUIT	*				m_graph; 
	SYMBOL_TABLE *			m_symbol_table;
	OPTIONS	*				m_options;


	/* functions to create the graph */
	PORT * 	add_output_port(NODE * node, const string & node_name);
	void 	add_input_ports_and_connect_to_graph(NODE * node,
										VARIABLE_STACK_TYPE * variable_name_stack);
	PORT *  add_an_input_port_and_connect_to_graph(NODE * node, 
										const string & node_name,
				    					const string & output_node_name);
	void 	add_clock_port_and_connect_to_clock(NODE * node, 
										const string & node_name,
										const string & global_clock_node_name);

	void 	connect_input_port_to_graph(PORT * input_port, 
										const string & input_node_name,
										const string & output_port_name);

	void 	add_edge_between_ports(const string & input_node_name,
							PORT * input_port,
							PORT * output_port);

	bool	open_circuit_input_file();

	string	get_edge_name(const string & input_node_name,
							const string & output_port_name);

};


#endif
