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



#ifndef graph_pieces_H
#define graph_pieces_H

#include "lut.h"
#include <set>
#include <deque>
#include <vector>
#include <string>
using namespace std;
#include "types.h"

class EDGE;
class PORT;
class NODE;

typedef vector<DELAY_TYPE> DELAYS;

typedef vector<EDGE> EDGE_VECTOR;
typedef vector<PORT> PORT_VECTOR;
typedef vector<NODE> NODE_VECTOR;

typedef vector<EDGE *> EDGES;
typedef vector<PORT *> PORTS;
typedef vector<NODE *> NODES;

typedef list<NODE *> NODE_PTR_LIST;
typedef list<PORT *> PORT_PTR_LIST;
typedef list<EDGE *> EDGE_PTR_LIST;

typedef set<CLUSTER_NUMBER_TYPE> CLUSTER_NUMBER_SET;
typedef vector<CLUSTER_NUMBER_TYPE> CLUSTER_NUMBERS;
typedef deque<NODE *> NODE_PTR_DEQUE;
//
// Class_name EDGE
//
// Description
// 		An edge connects input and output ports together
// 		with the combinational delay between its source port 
// 		and its sink port.
//		

class EDGE
{
public:
	EDGE();
	EDGE(const EDGE & another_edge);
	EDGE(PORT * source_port, PORT * sink_port, const string& edge_name, const LENGTH_TYPE& length);
	EDGE(PORT * source_port, PORT * sink_port, const LENGTH_TYPE& length);
	EDGE & operator=(const EDGE & another_edge);
	~EDGE();
	PORT *		get_source() const { return m_source;}
	PORT *		get_sink() const { return m_destination;}
	string		get_name() const { return m_name;}
	NODE *		get_sink_node() const;
	NODE *		get_source_node() const;
	string		get_edge_name(const string & input_node_name, const string & output_port_name) const;
	LENGTH_TYPE get_length() const { return m_length;}

	CLUSTER_NUMBER_TYPE get_sink_cluster_number() const;
	CLUSTER_NUMBER_TYPE get_source_cluster_number() const;
	CLUSTER_NUMBERS get_sink_sub_cluster_numbers() const;
	CLUSTER_NUMBERS get_source_sub_cluster_numbers() const;
	LENGTH_TYPE get_inter_cluster_distance() const;

	void 		set_length(const LENGTH_TYPE & new_length) { m_length = new_length;}
	bool		is_intra_cluster() const;
	bool		is_intra_cluster_and_intra_sub_cluster() const;
	bool		is_inter_cluster() const;
	bool		is_inter_sub_cluster() const;
	bool		is_sink_a_dff() const;
	bool		is_clock_edge() const;

private:
	PORT * 		m_source;
	PORT * 		m_destination;
	LENGTH_TYPE m_length;
	string		m_name;

};

//
// Class_name PORT
//
// Description
//		A port can be an input or output port or clock port.
//		PI and PO are output ports because they can have multiple edges
//		and it will make glueing subcircuits easier
//
//		The port can also be used in the future to add a pin capacitance
//

class PORT
{
public:
	enum PORT_TYPE {EXTERNAL, INTERNAL};
	enum IO_DIRECTION {INPUT, OUTPUT, CLOCK, UNKNOWN};
	enum EXTERNAL_TYPE {PO, PI, GI, GO, NONE};
	PORT();
	PORT(const string & port_name, const PORT_TYPE & port_type, 
		const IO_DIRECTION & io_direction, const EXTERNAL_TYPE & external_type);
	PORT(const string & port_name, const PORT_TYPE & port_type, 
		const IO_DIRECTION & io_direction, const EXTERNAL_TYPE & external_type,
		NODE * node_connected_to);
	PORT(const PORT & another_port);
	PORT & operator=(const PORT & another_port);
	~PORT();

	string				get_name() const {return m_name;}
	PORT_TYPE			get_type() const {return m_port_type;}
	IO_DIRECTION		get_io_direction() const {return m_io_direction;}
	NODE *				get_my_node() const {return m_my_node;}
	NODE *				get_node_that_fanout_to_me() const;        // used for input ports
	PORT *				get_output_port_that_fanout_to_me() const; //used for inputs ports
	EXTERNAL_TYPE		get_external_type() const {return m_external_type;}
	CLUSTER_NUMBER_TYPE get_cluster_number() const { return m_cluster_number;}
	CLUSTER_NUMBER_TYPE	get_last_sub_cluster_number() const { return m_sub_cluster_numbers.back(); }
	CLUSTER_NUMBERS	get_sub_cluster_numbers() const { return m_sub_cluster_numbers; }


	void set_my_node(NODE * new_node) { m_my_node = new_node;}
	void set_name(const string & new_name) { m_name = new_name;}
	void set_direction(const IO_DIRECTION & new_direction) {m_io_direction = new_direction;}
	void set_type(const PORT_TYPE new_port_type) { m_port_type = new_port_type;}
	void set_external_type(const EXTERNAL_TYPE & new_external_type) {m_external_type = new_external_type;}
	void set_cluster_number(const CLUSTER_NUMBER_TYPE & new_cluster_number) 
						{ m_cluster_number = new_cluster_number; }
	void add_sub_cluster_number(const CLUSTER_NUMBER_TYPE & new_cluster_number)
						{ m_sub_cluster_numbers.push_back(new_cluster_number); }
	void set_horizontal_position(const NUM_ELEMENTS& horizontal_pos) { m_horizontal_position = horizontal_pos; }

	// used for output ports as they can have multiple edges
	void			add_edge(EDGE * edge_to_add);
	EDGES get_edges() const;
	NUM_ELEMENTS	get_fanout_degree() const;
	NUM_ELEMENTS	get_fanout_degree_to_combinational_nodes() const;
	NUM_ELEMENTS 	get_horizontal_position() const { return m_horizontal_position; }


	// used for input ports as they have single edges
	void 	set_edge(EDGE * new_edge);
	void	remove_edge(EDGE * edge_to_remove);
	EDGE *	get_edge() const;

	NODE *	find_dff_in_fanout() const;

	bool	is_clock_port() const { return (m_io_direction == PORT::CLOCK); }
	bool	is_connected_to_PI() const;
private:
	string				m_name;
	PORT_TYPE			m_port_type;
	EXTERNAL_TYPE		m_external_type;
	IO_DIRECTION		m_io_direction;
	NODE *				m_my_node;
	EDGES				m_edges;
	CLUSTER_NUMBER_TYPE m_cluster_number;	// ports need cluster numbers because the pi
											// are sometimes seen as nodes for the purpose 
											// of shape
	CLUSTER_NUMBERS	m_sub_cluster_numbers;	
	NUM_ELEMENTS		m_horizontal_position;

	// obsolete
	// NUM_ELEMENTS	get_nGO() const;
};

//
// Class_name NODE
//
// Description
//		A node is a combination or sequential element.
//
class NODE 
{
public:
	enum NODE_TYPE {COMB, SEQ};
	enum COLOUR_TYPE {	NONE, 
						UNMARKED, MARKED, MARKED_VISITED, UNMARKED_VISITED,
					   	WHITE, GREY, BLACK,
						IN_PROGRESS, CLUSTER_DONE,MARKED_OUTCONE};
	NODE();
	NODE(const string & node_name);
	NODE(const string & node_name, const NODE_TYPE & node_type);
	NODE(const NODE & another_node);
	NODE& operator=(const NODE  & another_node);
	~NODE();

	PORT * 	create_and_add_port(const string & port_name, const PORT::PORT_TYPE & port_type,
								const PORT::IO_DIRECTION & io_direction, 
								const PORT::EXTERNAL_TYPE & external_type);
	void 	add_port(PORT * port_to_add);
	void	add_sub_cluster_number(const CLUSTER_NUMBER_TYPE & cluster_number);

	void set_lut(LUT *	new_lut);
	void set_type(const NODE_TYPE & new_type) {m_type = new_type;}
	void set_colour(const COLOUR_TYPE & new_colour_mark) {m_colour_mark = new_colour_mark;}
	void set_max_comb_delay_level(const DELAY_TYPE new_delay_level) { m_delay_level = new_delay_level;}
	void set_cluster_number(const CLUSTER_NUMBER_TYPE & new_cluster_number) 
								{ m_cluster_number = new_cluster_number; }
	void remove_input_port(PORT * port_to_remove); 	// doesn't delete the port
	void set_horizontal_position(const NUM_ELEMENTS& horizontal_pos) { m_horizontal_position = horizontal_pos; }


	string 			get_name() const { return m_name;}
	NODE_TYPE 		get_type() const { return m_type;}
	COLOUR_TYPE		get_colour() const { return m_colour_mark;}
	DELAY_TYPE		get_max_comb_delay_level() const { return m_delay_level;}

	CLUSTER_NUMBER_TYPE 	get_cluster_number() const { return m_cluster_number; }
	CLUSTER_NUMBER_TYPE		get_last_sub_cluster_number() const { return m_sub_cluster_numbers.back(); }
	CLUSTER_NUMBERS	get_sub_cluster_numbers() const { return m_sub_cluster_numbers; }
	EDGES			get_output_edges() const; 
	EDGES 			get_input_edges() const;
	EDGES 			get_input_edges_without_clock_edges() const;
	NUM_ELEMENTS    get_fanout_degree() const;
	NUM_ELEMENTS	get_fanout_degree_to_combinational_nodes() const;
	NUM_ELEMENTS	get_fanin_degree() const;
	NUM_ELEMENTS	get_horizontal_position() const { return m_horizontal_position; }
	NUM_ELEMENTS 	get_wirelength_approx_cost() const;

	string			get_info() const ;
	PORT *			get_output_port() const { return m_output_port;}
	PORTS			get_input_ports() const { return m_input_ports;}
	
	// for sequential nodes
	PORT *			get_D_port() const;
	PORT *			get_clock_port() const;

	void print_out_information() const;
private:
	NODE_TYPE		m_type;
	string			m_name;

	PORTS 			m_input_ports;
	PORT * 			m_output_port;

	DELAY_TYPE 		m_delay_level;
	COLOUR_TYPE		m_colour_mark;

	LUT * 			m_lut;

	CLUSTER_NUMBER_TYPE	m_cluster_number;
	CLUSTER_NUMBERS		m_sub_cluster_numbers;	

	NUM_ELEMENTS		m_horizontal_position;
};

#endif
