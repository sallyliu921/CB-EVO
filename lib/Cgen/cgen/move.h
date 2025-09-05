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



#ifndef move_H
#define move_H

class MOVE;

#include "types.h"
#include "edges_and_nodes.h"
#include "matrix.h"

//
// Class_name MOVE
//
// Description
//
//	Moves for iteratively improving the combinational delay structure solution.
//	Holds information about potential moves, the cost of a potential move.
//	Also has the ability to execute the move.
//
//
//
//
//	ie. Weight is subtracted from the source edge and 
//		weight is added to the sink edge
//


class MOVE
{
public:
	MOVE();
	~MOVE();

	void make_move();
	void make_temporary_move(MATRIX& violation_matrix);
	void unmake_temporary_move(MATRIX& violation_matrix);
	COST_TYPE get_changed_cost(MATRIX& violation_matrix) const; // returns the cost difference

	EDGE * get_source_edge() const { return m_source_edge; }
	EDGE * get_sink_edge() const { return m_sink_edge; }

	void set_source_edge(EDGE * source_edge) { m_source_edge = source_edge; }
	void set_sink_edge(EDGE * sink_edge) { m_sink_edge = sink_edge; }

	bool is_valid() const;
	bool edge_changed_clusters() const;

	void clear();
private:
	MOVE(const MOVE & another_move);
	MOVE & operator=(const MOVE & another_move);

	EDGE * m_source_edge;
	EDGE * m_sink_edge;
	EDGE * m_old_source_edge;
	EDGE * m_old_sink_edge;

	void update_violation_matrix(MATRIX& violation_matrix, EDGE * edge, NUM_ELEMENTS change);
};


#endif
