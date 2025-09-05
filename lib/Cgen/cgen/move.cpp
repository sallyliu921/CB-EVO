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



#include "move.h"
#include "gen.h"

MOVE::MOVE()
{
	m_source_edge	= 0;
	m_sink_edge		= 0;
}
MOVE::MOVE(const MOVE & another_move)
{
	m_source_edge	= another_move.m_source_edge;
	m_sink_edge		= another_move.m_sink_edge;
}

MOVE & MOVE::operator=(const MOVE & another_move)
{
	m_source_edge	= another_move.m_source_edge;
	m_sink_edge		= another_move.m_sink_edge;

	return (*this);
}

MOVE::~MOVE()
{
}

// POST: the move has been made
void MOVE::make_move()
{

	assert(m_source_edge && m_sink_edge);
	assert(m_source_edge->can_reduce_weight_on_edge());

	m_source_edge->add_weight(-1);

	m_sink_edge->add_weight(1);

	m_old_source_edge = m_source_edge;
	m_old_sink_edge = m_sink_edge;

	m_source_edge = 0;
	m_sink_edge = 0;
}


// POST: the move has temporarily made
void MOVE::make_temporary_move
(
	MATRIX& violation_matrix
)
{
	m_source_edge->add_weight(-1);
	m_sink_edge->add_weight(1);

	update_violation_matrix(violation_matrix, m_source_edge, -1);
	update_violation_matrix(violation_matrix, m_sink_edge, +1);
}

// POST: the temporarily move has been unmade
void MOVE::unmake_temporary_move
(
	MATRIX& violation_matrix
)
{
	m_source_edge->add_weight(1);
	m_sink_edge->add_weight(-1);

	update_violation_matrix(violation_matrix, m_source_edge, +1);
	update_violation_matrix(violation_matrix, m_sink_edge, -1);
}

// POST: the violation matrix has been updated
void MOVE::update_violation_matrix
(
	MATRIX& violation_matrix,
	EDGE * edge,
	NUM_ELEMENTS change
)
{
	CLUSTER_NUMBER_TYPE source_cluster = edge->get_source_cluster_number();
	CLUSTER_NUMBER_TYPE sink_cluster = edge->get_sink_cluster_number();

	violation_matrix(source_cluster, sink_cluster) += change;
}


// a move is valid if it has a source edge and sink edge
// its source edge has weight to give
// its sink edge has weight to take
// and the source edge has delay to spare if the edge length is 1
//
// RETURN: true if value, false otherwise
bool MOVE::is_valid() const
{ 
	bool valid_move = 	(m_source_edge != 0 && m_sink_edge != 0) &&
						( m_source_edge != m_sink_edge) &&
						(m_source_edge->get_weight() > 0) &&
						(m_sink_edge->get_weight() < m_sink_edge->get_max_weight()) &&
						(m_source_edge->get_sink_level_node()->get_delay_difference() > 0 
						 || m_source_edge->get_length() != 1);

	return valid_move;
}


// RETURN: true if the source or sink changed clusters
bool MOVE::edge_changed_clusters() const
{
	bool cluster_change =  
		((m_source_edge->get_source_cluster_number() != m_sink_edge->get_source_cluster_number()) ||
		 (m_source_edge->get_sink_cluster_number() != m_sink_edge->get_sink_cluster_number()));
	
	return cluster_change;
}

// POST: the move has been cleared
void MOVE::clear()
{
	m_source_edge = 0;
	m_sink_edge = 0;
}
