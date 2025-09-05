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



#ifndef wirelength_character_H
#define wirelength_character_H


//
// Class_name WIRELENGTH_CHARACTER
//
// Description
//


#include "circ.h"
#include "circuit.h"
#include "rand.h"

class WIRELENGTH_MOVE
{
public:
	WIRELENGTH_MOVE();
	~WIRELENGTH_MOVE();

	void set_node_a(NODE * node_a) { m_node_a = node_a; }
	void set_port_a(PORT * port_a) { m_port_a = port_a; }
	void set_node_b(NODE * node_b) { m_node_b = node_b; }
	void set_port_b(PORT * port_b) { m_port_b = port_b; }

	NODE * get_node_a() const { return m_node_a; }
	PORT * get_port_a() const { return m_port_a; }
	NODE * get_node_b() const { return m_node_b; } 
	PORT * get_port_b() const { return m_port_b; } 

	bool is_valid() const;
	void clear();

	void make_and_commit_move();
	void make_move();
	void make_temporary_move();
	void unmake_temporary_move();

	COST_TYPE get_changed_wirelength_cost(const COST_TYPE& wirelength_cost,
										const COST_TYPE& max_index, const COST_TYPE& nEdges);
	COST_TYPE get_wirelength_change_of_nodes_involved(); 
	COST_TYPE get_wirelength_cost_of_nodes_involved();
	COST_TYPE get_wirelength_cost_of_node_and_fanout(NODE * node);
	COST_TYPE get_wirelength_cost_of_fanout_of_PI(const PORT * port) const;
private:
	NODE * 	m_node_a;
	PORT * 	m_port_a;
	NODE *  m_node_b;		
	PORT *  m_port_b;
};
class WIRELENGTH_CHARACTER
{
public:
	WIRELENGTH_CHARACTER();
	WIRELENGTH_CHARACTER(const WIRELENGTH_CHARACTER & another_wirelength_character);
	WIRELENGTH_CHARACTER & operator=(const WIRELENGTH_CHARACTER & another_wirelength_character);
	~WIRELENGTH_CHARACTER();

	double get_circuit_wirelength_approx(CIRCUIT * circuit);
private:
	CIRCUIT * m_circuit;
	NUM_ELEMENTS m_width;
	RANDOM_NUMBER * m_rand_number_gen;
	WIRELENGTH_MOVE m_move;

	NODES m_nodes;
	PORTS m_PI;
	CLUSTERS m_clusters;

	void set_horizontal_position_of_nodes();
	double get_wirelength_measurement() const;
	void iterate();
	void generate_move();
	COST_TYPE get_cost() const;
	COST_TYPE get_changed_cost(const COST_TYPE& old_cost);
	void update_temperature(double& temperature, const double& success_rate, const NUM_ELEMENTS& loops);
	void print_wirelength_results() const;
};


#endif
