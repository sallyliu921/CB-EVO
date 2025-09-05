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



#ifndef drawer_H
#define drawer_H

#include "circ.h"
#include "circuit.h"
#include <fstream>

//
// Class_name DRAWER
//
// Description
//
//		Draws pictures of the circuit.

class DRAWER
{
public:
	DRAWER();
	DRAWER(const DRAWER & another_drawer);
	DRAWER & operator=(const DRAWER & another_drawer);
	~DRAWER();

	void draw_full_graph(CIRCUIT * circuit);	
private:
	fstream 	m_output_file;
	bool		m_display_illegal_char_msg;

	void preamble();
	void draw_input(PORT * input_port);
	void draw_node(NODE * node);
	void draw_node_fanin_edges(NODE * sink_node);
	void draw_edge(NODE * source_node, NODE * sink_node);
	void constrain_io_ranks(CIRCUIT * circuit);
	void cleanup();

	string print_name(string node_name);
	void open_output_file();

};


#endif
