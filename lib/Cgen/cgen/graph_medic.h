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



#ifndef graph_medic_H
#define graph_medic_H

#include "gen.h"
#include "graph.h"

//
// Class_name GRAPH_MEDIC
//
// Description:
// 	This class cleans up the graph after being constructed.
// 	It cleans it up by lopping off any node that is unreachable
// 	from a PI or PO.
//
// 	It also can delete any buffer or inverter nodes
//
//

class GRAPH_MEDIC
{
public:
	GRAPH_MEDIC(CIRCUIT_GRAPH * circuit);
	GRAPH_MEDIC(const GRAPH_MEDIC & another_graph_medic);
	GRAPH_MEDIC & operator=(const GRAPH_MEDIC & another_graph_medic);
	~GRAPH_MEDIC();

	void		delete_buffer_and_inverter_nodes();
	void		delete_unusable_nodes();
	void		check_sanity();				// checks the sanity of the graph
private:
	CIRCUIT_GRAPH *	m_graph;
	NODES			m_nodes_to_delete;

	NUM_ELEMENTS	m_number_comb_deleted;
	NUM_ELEMENTS	m_number_seq_deleted;
	NUM_ELEMENTS	m_number_unreachable;
	NUM_ELEMENTS	m_number_isolated;

	enum PROBLEM_NODES	{ 	NODES_NOT_CONNECTED_TO_PO, 
							NODES_NOT_CONNECTED_TO_PI,
							ISOLATED_NODES};

	// used to display the right error message to the user
	PROBLEM_NODES			m_problem_nodes;

	// used to allow a single warning message about node deletion
	bool	m_shown_eliminate_node_warning;	

	// functions to remove unnecessary nodes 
	void	delete_unreachable_nodes_from_inputs();
	void	delete_unreachable_nodes_from_outputs();
	void	delete_unconnected_external_ports();
	void	delete_isolated_nodes();
	void	delete_isolated_ports();

	void	mark_up(const NODE::COLOUR_TYPE & colour);
	void	mark_up_from_node(NODE * node, const NODE::COLOUR_TYPE & colour);
	void	eliminate_down();
	void	eliminate_down_from_node(NODE * node);
	void 	mark_down(const NODE::COLOUR_TYPE & colour);
	void 	mark_down_from_node(NODE * node, const NODE::COLOUR_TYPE & colour);
	void 	eliminate_up();
	void	eliminate_up_from_node(NODE * node);

	void	if_unmarked_queue_for_deletion(NODE * node);

	void	delete_queued_nodes();
	void	delete_node(NODE * node);
	void 	detach_node_from_references_in_fanin(NODE * node);
	void	detach_node_from_references_in_fanout(NODE * node);

	// warnings to the user
	void	show_node_deletion_warning(NODE * node);
	void 	print_node_deletion_stats();

};


#endif
