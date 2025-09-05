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




#include "symbol_table.h"

SYMBOL_TABLE::SYMBOL_TABLE()
{
}

SYMBOL_TABLE::SYMBOL_TABLE(const SYMBOL_TABLE & another_symbol_table)
{
	assert(false);
}

SYMBOL_TABLE & SYMBOL_TABLE::operator=(const SYMBOL_TABLE & another_symbol_table)
{
	assert(false);
	return(*this);
}

SYMBOL_TABLE::~SYMBOL_TABLE()
{
}

void SYMBOL_TABLE::insert_port
(
	const string & port_name,
	PORT * port
)
{	
	debugif(DSYMBOL_TABLE,"Symbol Table: Inserting port name = '" << port_name << "'");
	m_port_symbol_table[port_name] = port;
}

void SYMBOL_TABLE::insert_edge
(
	const string & edge_name,
	EDGE * edge
)
{
	debugif(DSYMBOL_TABLE,"Symbol Table: Inserting edge name = '" << edge_name << "'");
	m_edge_symbol_table[edge_name] = edge;
}

void SYMBOL_TABLE::insert_node
(
	const string & node_name,
	NODE * node
)
{
	debugif(DSYMBOL_TABLE,"Symbol Table: Inserting node name = '" << node_name << "'");
	m_node_symbol_table[node_name] = node;
}

PORT *	 SYMBOL_TABLE::query_for_port
(
	const string & port_name
)
{
	debugif(DSYMBOL_TABLE,"Symbol Table: Query for port name = '" << port_name << "'");
	return m_port_symbol_table[port_name];
}

EDGE *	 SYMBOL_TABLE::query_for_edge
(
	const string & edge_name
)
{
	debugif(DSYMBOL_TABLE,"Symbol Table: Query for edge name = '" << edge_name << "'");
	return m_edge_symbol_table[edge_name];
}

NODE *	 SYMBOL_TABLE::query_for_node
(
	const string & node_name
)
{
	debugif(DSYMBOL_TABLE,"Symbol Table: Query for node name = '" << node_name << "'");
	return m_node_symbol_table[node_name];
}

void SYMBOL_TABLE::remove_port
(
	const string & port_name
)
{	
	PORT_HASH_TABLE_ITER port_to_remove;
	
	debugif(DSYMBOL_TABLE,"Symbol Table::Removing port name = '" << port_name << "'");
	debugif(DSYMBOL_TABLE,"SYMBOL_TABLE::remove_port == DEBUG");
	port_to_remove = m_port_symbol_table.find(port_name);
	if (port_to_remove == m_port_symbol_table.end()){
		debugif(DSYMBOL_TABLE,"SYMBOL_TABLE::port not found");
	}else{
		debugif(DSYMBOL_TABLE,"SYMBOL_TABLE::remove_port Table size == "<< m_port_symbol_table.size());
		if (m_port_symbol_table.size()) m_port_symbol_table.erase(port_to_remove); //risky, 2021/5/6, don't understand why segmentation error
		debugif(DSYMBOL_TABLE,"SYMBOL_TABLE::remove_port == DEBUG");
	}
	
}
 
void SYMBOL_TABLE::remove_edge
(
	const string & edge_name
)
{
	EDGE_HASH_TABLE_ITER edge_to_remove;

	debugif(DSYMBOL_TABLE,"Symbol Table::Removing edge name = '" << edge_name << "'");
	debugif(DSYMBOL_TABLE,"SYMBOL_TABLE::remove_edge == DEBUG");
	edge_to_remove = m_edge_symbol_table.find(edge_name);
	debugif(DSYMBOL_TABLE,"SYMBOL_TABLE::remove_edge Table size == "<< m_edge_symbol_table.size());
	if (m_edge_symbol_table.size()) m_edge_symbol_table.erase(edge_to_remove);
	debugif(DSYMBOL_TABLE,"SYMBOL_TABLE::remove_edge == DEBUG");
}

void SYMBOL_TABLE::remove_node
(
	const string & node_name
)
{
	NODE_HASH_TABLE_ITER node_to_remove;

	debugif(DSYMBOL_TABLE,"Symbol Table::Removing node name = '" << node_name << "'");
	debugif(DSYMBOL_TABLE,"SYMBOL_TABLE::remove_node == DEBUG");
	node_to_remove = m_node_symbol_table.find(node_name);
	debugif(DSYMBOL_TABLE,"SYMBOL_TABLE::remove_node Table size == "<< m_node_symbol_table.size());
	if (m_node_symbol_table.size()) m_node_symbol_table.erase(node_to_remove);
	debugif(DSYMBOL_TABLE,"SYMBOL_TABLE::remove_node == DEBUG");
}

NODE *	SYMBOL_TABLE::front_node()
{
	return (*m_node_symbol_table.begin()).second;
}

PORT * 	SYMBOL_TABLE::front_port()
{

	return (*m_port_symbol_table.begin()).second;
}

EDGE * 	SYMBOL_TABLE::front_edge()
{
	return (*m_edge_symbol_table.begin()).second;
}

void	SYMBOL_TABLE::pop_front_node()
{
	m_node_symbol_table.erase(m_node_symbol_table.begin());
}

void	SYMBOL_TABLE::pop_front_port()
{
	m_port_symbol_table.erase(m_port_symbol_table.begin());
}

void	SYMBOL_TABLE::pop_front_edge()
{
	m_edge_symbol_table.erase(m_edge_symbol_table.begin());
}

bool	SYMBOL_TABLE::empty_of_nodes()
{
	return m_node_symbol_table.empty();
}
bool	SYMBOL_TABLE::empty_of_ports()
{
	return m_port_symbol_table.empty();
}

bool	SYMBOL_TABLE::empty_of_edges()
{
	return m_edge_symbol_table.empty();
}
