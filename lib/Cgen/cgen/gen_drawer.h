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



#ifndef gen_drawer_H
#define gen_drawer_H

#include "gen.h"

//
// Class_name GEN_DRAWER
//
// Description
//
//		Draws pictures of the circuit.

class GEN_DRAWER
{
public:
	GEN_DRAWER();
	GEN_DRAWER(const GEN_DRAWER & another_drawer);
	GEN_DRAWER & operator=(const GEN_DRAWER & another_drawer);
	~GEN_DRAWER();

	void draw_full_graph(CIRCUIT_GRAPH * circuit);	
private:
	fstream 		m_output_file;
	CIRCUIT_GRAPH* 	m_circuit;

	void preamble();
	void draw_input(EDGE * edge);
	void draw_node(NODE * node);
	void draw_node_fanin_edges(NODE * sink_node);
	void draw_edge(NODE * source_node, NODE * sink_node);
	void constrain_io_ranks();
	void cleanup();

	string print_name(string node_name);
	void open_output_file();
};


#endif
