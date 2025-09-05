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





#ifndef symbol_H
#define symbol_H

#include "circ.h"
//#include <hash_map>
#include <functional>
#include <unordered_map>

//
// Class_name SYMBOL_TABLE
//
// Description
//
//	Contains information during graph construction.
//	Based on hash_map from the stl library.
//
//	If you are trying to compile this for VC++ and do not have 
//	access to the hash_map class please download and use STLport.
//	


//#ifndef VISUAL_C
# if 0
template <>
class std::hash<string>
{
	public:
		size_t operator()(const string &str) const
		{
			hash<char const *> h;

			return (h(str.c_str()));
		}
};
#endif

struct eqstr
{
	bool operator()(const string& s1, const string& s2) const
	{
		return (s1 == s2);
	}
};

typedef unordered_map<string, PORT *, hash<string>, eqstr> PORT_HASH_TABLE;
typedef unordered_map<string, EDGE *, hash<string>, eqstr> EDGE_HASH_TABLE;
typedef unordered_map<string, NODE *, hash<string>, eqstr> NODE_HASH_TABLE;

typedef PORT_HASH_TABLE::iterator PORT_HASH_TABLE_ITER;
typedef NODE_HASH_TABLE::iterator NODE_HASH_TABLE_ITER;
typedef EDGE_HASH_TABLE::iterator EDGE_HASH_TABLE_ITER;

class SYMBOL_TABLE
{
public:
	SYMBOL_TABLE();
	~SYMBOL_TABLE();

	void	insert_port(const string & port_name, PORT * port);
	void	insert_edge(const string & edge_name, EDGE * edge);
	void	insert_node(const string & node_name, NODE * node);

	PORT *	query_for_port(const string & port_name);
	EDGE *	query_for_edge(const string & edge_name);
	NODE *	query_for_node(const string & node_name);

	void	remove_port(const string & port_name);
	void	remove_edge(const string & edge_name);
	void	remove_node(const string & node_name);

	// Returns the first element in the symbol table
	NODE *	front_node();
	PORT * 	front_port();
	EDGE * 	front_edge();

	// Removes the first element in the symbol table
	void	pop_front_node();
	void	pop_front_port();
	void	pop_front_edge();

	bool	empty_of_nodes();
	bool	empty_of_ports();
	bool	empty_of_edges();
private:
	PORT_HASH_TABLE		m_port_symbol_table;
	EDGE_HASH_TABLE		m_edge_symbol_table;
	NODE_HASH_TABLE		m_node_symbol_table;

	SYMBOL_TABLE(const SYMBOL_TABLE & another_symbol);
	SYMBOL_TABLE & operator=(const SYMBOL_TABLE & another_symbol);
};
 

#endif
