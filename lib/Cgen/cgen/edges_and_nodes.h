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



#ifndef EDGES_AND_NODES_H
#define EDGES_AND_NODES_H

#include "types.h"
#include <list>

class EDGE;
class NODE;
struct LUT_FUNC;

typedef list<EDGE *> EDGE_LIST;
typedef list<NODE *> NODE_LIST;

typedef vector<EDGE *> EDGES;
typedef vector<NODE *> NODES;

#include "level_node.h"

//  Hook for multiple LUT types later.
typedef enum {NODE_LUT, NODE_DFF} CTYPE;
struct LUT_FUNC
{
  int type;    	
};

//
// Class_name NODE
//
// Description
//
//  The indiv. nodes in the graph
//	In gen nodes are also primary inputs (but not primary outputs)
//  

class NODE
{
public:
	enum COLOUR_TYPE {	NONE, UNMARKED, MARKED, MARKED_VISITED, UNMARKED_VISITED};
	NODE(const CLUSTER_NUMBER_TYPE& cluster_number, const DELAY_TYPE& delay_level, 
		 const ID_TYPE & id, const NUM_ELEMENTS& kin, LEVEL_NODE * my_level_node);
	~NODE();

	CLUSTER_NUMBER_TYPE	get_cluster_number() const { return m_cluster_number; }
	ID_TYPE				get_id() const { return m_id; }
	NUM_ELEMENTS		get_fanin_degree() const { return m_fanin_degree; }	
	NUM_ELEMENTS		get_fanout_degree() const { return m_fanout_degree; }	
	DELAY_TYPE			get_delay_level() const { return m_delay_level; }
	NUM_ELEMENTS		get_edges_to_assign() const { return m_edges_to_assign; }
	NUM_ELEMENTS		get_edges_currently_assigning() const { return m_edges_currently_assigning; }
	NUM_ELEMENTS		get_unassigned_long_and_inter_cluster_edges() const 
							{ return m_unassigned_long_and_inter_cluster_edges; }
	NUM_ELEMENTS		get_nLost_simple_edges() const { return m_lost_edges; }
	LOCALITY 			get_index() const { return m_index; }
	LOCALITY			get_span() const { return m_span; }
	LOCALITY			get_rspan() const { return m_rspan; }
	LOCALITY			get_lspan() const { return m_lspan; }
	string				get_name() const;
	EDGES				get_input_edges() const { return m_input_edges; }
	EDGES				get_output_edges() const { return m_output_edges; }
	EDGE*				get_output_edge_that_connects_to_DFF() const;
	EDGE*				get_dff_input_edge() const;
	string				get_input_node_names() const;
	NUM_ELEMENTS		get_nInputs() const { return m_input_edges.size(); }
	NUM_ELEMENTS		get_nOutputs() const { return m_output_edges.size(); }
	NUM_ELEMENTS		get_nComb_Outputs() const;
	NUM_ELEMENTS		get_nDouble_inputs() const;
	NUM_ELEMENTS		get_double_input_cost() const;
	EDGES				get_double_input_edges() const;
	NUM_ELEMENTS		get_nExcess_inputs() const;
	COST_TYPE			get_wirelength_cost() const;
	NUM_ELEMENTS		get_nInputs_still_needed() const;
	LEVEL_NODE * 		get_my_level_node() const { return m_my_level_node; }
	COLOUR_TYPE			get_colour() const { return m_colour; }


	void set_fanout_degree(const NUM_ELEMENTS& fanout_degree);
	void set_edges_to_assign(const NUM_ELEMENTS& edges_to_assign) { m_edges_to_assign = edges_to_assign; }
	void set_edges_currently_assigning(const NUM_ELEMENTS& edges_to_assign) { m_edges_currently_assigning = 0; }
	void add_edges_to_assign(const NUM_ELEMENTS& edges_to_assign);
	void add_edges_currently_assigning(const NUM_ELEMENTS& edges_to_assign);
	void set_is_latched(const bool& is_latched) { m_is_latched = is_latched; }
	void set_is_PI(const bool& is_PI) { m_is_PI = is_PI; }
	void set_is_PO(const bool& is_PO) { m_is_PO = is_PO; }
	void set_is_dff(const bool& is_dff) { m_is_dff = is_dff; }
	void set_is_latch_taken(const bool& is_latch_taken) { m_is_latch_taken = is_latch_taken; }
	void set_colour(const NODE::COLOUR_TYPE& new_colour) { m_colour = new_colour; }
	void add_unassigned_long_and_inter_cluster_edges(const NUM_ELEMENTS& nEdges);
	void add_fanin_degree(const NUM_ELEMENTS& fanin) { m_fanin_degree += fanin; }
	void add_lost_simple_edges(const NUM_ELEMENTS& edges_lost) { m_lost_edges += edges_lost; }
	
	void set_index(const LOCALITY& new_index) { m_index = new_index; }
	void set_span(const LOCALITY& new_span) { m_span = new_span; }
	void set_rspan(const LOCALITY& new_rspan) { m_rspan = new_rspan; }
	void set_lspan(const LOCALITY& new_lspan) { m_lspan = new_lspan; }

	bool 	is_combinational() const { return ! m_is_dff; }
	bool 	is_DFF() const { return m_is_dff; }
	bool	is_PI() const { return m_is_PI; }
	bool	is_PO() const { return m_is_PO; }
	bool	is_latched() const { return m_is_latched;}
	bool	is_latch_taken() const { return m_is_latch_taken; }
	bool	is_marked() const { return m_mark; }
	bool	is_buffer() const { return (get_nInputs() == 1 && is_combinational()); }
	bool	is_effectively_buffer() const;
	bool	has_double_inputs() const;
	bool	has_excess_inputs() const;
	bool	needs_inputs() const;
	bool	is_dff_have_loop_back_to_itself() const;

	void	add_input_edge(EDGE * input_edge);
	void	add_output_edge(EDGE * output_edge);


	void	remove_input_edge(EDGE * input_edge);
	void	remove_output_edge(EDGE * output_edge);

	bool	are_nodes_connected(const NODE * possible_source_node) const;

	NODE*  	get_source_node_of_DFF() const; // if the node is a dff get its source node
private:
	string 				m_name;			// name of nodes smallest sub-circuit 

	CLUSTER_NUMBER_TYPE m_cluster_number;
	DELAY_TYPE  		m_delay_level;		// combinational delay level
	NUM_ELEMENTS		m_fanin_degree;	
	NUM_ELEMENTS		m_fanout_degree;
	NUM_ELEMENTS 		m_edges_to_assign;	// unconnected fanouts    

	
	LUT_FUNC 			m_function;		// see lut.c, not impemented yet     
	ID_TYPE   			m_id;			// unique id                

	EDGES				m_input_edges;
	EDGES				m_output_edges;

	bool  m_is_PI, m_is_PO;		// T if node is such a thing    
	bool  m_is_latched;
	bool  m_mark;				// generic mark
	bool  m_is_dff;
	bool  m_is_latch_taken;
	

	LOCALITY m_index;			// Used for locality
	LOCALITY m_span, m_lspan, m_rspan;	// Used for locality
	NUM_ELEMENTS 	m_unassigned_long_and_inter_cluster_edges;
	NUM_ELEMENTS	m_edges_currently_assigning;

	COLOUR_TYPE 	m_colour;

	NUM_ELEMENTS	m_kin;
	NUM_ELEMENTS	m_lost_edges;	
	EDGES			m_bad_edges;		// edges that cause double connections
										// or edges that are one too many inputs
	LEVEL_NODE *	m_my_level_node;

	NODE(const NODE  & another_node);
	NODE & operator=(const NODE  & another_node);
};	

//
// Class_name EDGE
//
// Description
//
//  Edge in a graph.
//	Edges can be simple edges (node to node) or super edges (level node to level node) 
//	Super edges represent an aggregate of simple edges.
//

class EDGE
{
public:
	enum EDGE_TYPE {NODE_TO_NODE, LEVEL_NODE_TO_LEVEL_NODE};
	EDGE();
	EDGE(LEVEL_NODE* source_node, LEVEL_NODE* sink_node, const ID_TYPE& id);
	EDGE(NODE* source_node, NODE* sink_node);
		 
	~EDGE();

	LENGTH_TYPE get_length() const { return m_length; }
	NUM_ELEMENTS get_weight() const { return m_weight; }
	NUM_ELEMENTS get_max_weight() const { return m_max_weight; }
	LEVEL_NODE * get_source_level_node() const { return m_source_level_node; }
	LEVEL_NODE * get_sink_level_node() const { return m_sink_level_node; }
	NODE * get_source_node() const { return m_source_node; }
	NODE * get_sink_node() const { return m_sink_node; }
	ID_TYPE	get_id() const { return m_id; }
	string get_info() const;
	NUM_ELEMENTS get_total_simple_edge_weight() const;
	EDGE * get_my_super_edge() const { return m_my_super_edge; }

	void set_weight(const NUM_ELEMENTS& weight);
	void unsafe_set_weight(const NUM_ELEMENTS& weight); // hack
	void set_weight_and_update_level_nodes(const NUM_ELEMENTS& weight); 
	void set_max_weight(const NUM_ELEMENTS& weight) { m_max_weight = weight; }
	void add_weight(const NUM_ELEMENTS& weight_change); // sets the weight and adjusts the level nodes 
	void add_simple_edge(EDGE * simple_edge);			// for super edges
	void add_weight_to_dff_super_edge(const NUM_ELEMENTS& weight_change);
	void set_length(const LENGTH_TYPE& length) { m_length = length; }

	bool is_intra_cluster() const;
	bool is_inter_cluster() const;
	bool is_inter_cluster_input(const CLUSTER_NUMBER_TYPE & input_cluster_number) const;
	bool is_inter_cluster_output(const CLUSTER_NUMBER_TYPE & output_cluster_number) const;
	bool is_sink_a_dff() const;

	bool has_spare_weight() const; // does the edge and level node have weight to spare?
	bool can_reduce_weight_on_edge() const; // does the edge have weight and delay to spare?
	bool is_valid_edge_to_move() const;

	NUM_ELEMENTS get_nUnneeded_edges() const; // how much weight does the edge have that it can delete
	NUM_ELEMENTS get_nEdge_to_sacrifice() const; //how much weight can we remove without sacraficing the delay 
												 // into the node

	CLUSTER_NUMBER_TYPE get_sink_cluster_number() const;
	CLUSTER_NUMBER_TYPE get_source_cluster_number() const;

	// get the contribution of the edge towards the excessive cost at its level nodes
	COST_TYPE get_selection_weight() const;

	void	set_sink_node(NODE * node);
	void 	set_my_super_edge(EDGE * super_edge);

	void 	check_sanity() const;

	string get_vhdl_signal_name() const;

private:
	EDGE_TYPE		m_type;
	LENGTH_TYPE		m_length;
	NUM_ELEMENTS	m_weight;
	NUM_ELEMENTS	m_max_weight;
	NODE 	*		m_source_node;
	NODE 	* 		m_sink_node;
	LEVEL_NODE *	m_source_level_node;
	LEVEL_NODE *	m_sink_level_node;
	EDGE * 			m_my_super_edge;
	ID_TYPE			m_id;



	EDGES			m_simple_edges;		// for super edges



	EDGE(const EDGE  & another_edge);
	EDGE & operator=(const EDGE  & another_edge);
};



#endif
