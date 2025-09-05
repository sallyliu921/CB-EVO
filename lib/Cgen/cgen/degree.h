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


#ifndef degree_partitioner_H
#define degree_partitioner_H


#include "gen.h"

//
// Class_name DEGREE_MOVE
//
// Description
//
//	Holds information about potential moves, the cost of a potential move.
//	Also has the ability to execute the move.
//
//	A move is switch of fanout degrees between the source level node
//	(that gives up high degree fanouts) and the sink level node 
//	(that gives up low degree fanouts).
//


class DEGREE_MOVE
{
public:
	DEGREE_MOVE();
	~DEGREE_MOVE();

	DELAY_TYPE 		get_source_level() const { return m_source_level; }
	DELAY_TYPE 		get_dest_level() const { return m_dest_level; }
	NUM_ELEMENTS 	get_source_degree() const { return m_source_degree; }
	NUM_ELEMENTS 	get_dest_degree() const { return m_dest_degree; }

	void set_source_level(const DELAY_TYPE& source_level)  { m_source_level = source_level; }
	void set_dest_level(const DELAY_TYPE& dest_level)  { m_dest_level = dest_level; }
	void set_source_degree(const NUM_ELEMENTS& source_degree) {  m_source_degree = source_degree; }
	void set_dest_degree(const NUM_ELEMENTS& dest_degree)  { m_dest_degree = dest_degree; }

	bool is_valid() const;

	void set_invalid() { m_valid = false; }

	void reset();
private:
	DEGREE_MOVE(const DEGREE_MOVE & another_move);
	DEGREE_MOVE & operator=(const DEGREE_MOVE & another_move);


	DELAY_TYPE 		m_source_level;
	DELAY_TYPE 		m_dest_level;
	NUM_ELEMENTS 	m_source_degree;
	NUM_ELEMENTS 	m_dest_degree;

	bool			m_valid;

};

//
// Class_name DEGREE_PARTITIONER
//
// Description
//
//  Partition the fanout distribution for the cluster so that each level node gets 
//  an fanout for each individual node in it, and so that these fanouts
//  are feasible given the number of edges that are exiting that level.
//


class DEGREE_PARTITIONER
{
public:
	DEGREE_PARTITIONER();
	~DEGREE_PARTITIONER();
	NUM_ELEMENTS degree_partition(CLUSTER *cluster);

	COST_TYPE get_changed_cost();

private:
	DELAY_TYPE 		m_max_comb_delay;	// shortcut to cluster->delay               
	SHAPE 			m_target;	 		// target total out-degree for level i 
	SHAPE			m_unassigned_nodes; // number of unassigned nodes at level i            
	SHAPE			m_sum_assigned;		// sum of assigned out-deg for level i   
	SHAPE			m_maxdeg;	 		// the maximum  out-degree for level i 
	SHAPE			m_mindeg; 			// the minimum  out-degree for level i     
	DOUBLE_VECTOR 	m_avgdeg;			// the average  out-degree for level i     
	DISTRIBUTION 	m_fanout_distribution;	 		// the available out-degrees from spec     
	NUM_ELEMENTS 	m_max_fanout; 					// the maximum fanout left to be assigned or 
	                                                // the maximum fanout in the cluster
	                                                   
	NUM_ELEMENTS	m_num_out_degrees_to_assign;	// number of fanouts to assign
	DOUBLE_VECTOR	m_pmf; 				// prob. mass function for ranomly selecting certain things
	SHAPE			m_current_maxdeg;	 			// the maximum  out-degree for level i 

	double			m_delta;			

	DEGREE_MOVE		m_move;
	bool			m_iterating;
	NUM_ELEMENTS	m_max_latched_degree;

	CLUSTER *	m_cluster;


	bool			m_sanity_check_fatal;

	// Create the initial solution
	void create_initial_solution();
	void initialize_variables();
	void compute_targets();
	void make_initial_assignments();
	void make_latch_and_zero_degree_assignments();
	void initial_assignments(const DELAY_TYPE& lev);
	void show_vectors() const;
	void assign_one_out();
	void make_assignment(const NUM_ELEMENTS& deg, const DELAY_TYPE& lev);
	void assign_to_level(const NUM_ELEMENTS& degree, 
						const DELAY_TYPE& lev, const NUM_ELEMENTS& num);
	void recalculate_for_level(const DELAY_TYPE& lev);
	void recalculate_max_fanout();
	void recalculate_max_degree();


	// Iterative improvement
	void simulated_anneal();
	COST_TYPE get_cost() const;
	COST_TYPE get_cost_fanout_penalty() const;
	COST_TYPE get_cost_degree_fanout_edge_misassignment() const;
	COST_TYPE get_changed_cost(const COST_TYPE& cost);
	void generate_move();
	void make_move();
	NUM_ELEMENTS choose_source_degree( const DELAY_TYPE& source_level) const;
	NUM_ELEMENTS choose_dest_degree( const DELAY_TYPE& sink_level, const NUM_ELEMENTS& source_degree);

	// Post-processing
	void fix_degree_assignments();
	void fix_max_degree_assignments();
	void prevent_double_connections();

	void find_out_how_many_additional_long_edges_required(const LEVEL_NODE * level_node,
			 		NUM_ELEMENTS& nAdditional_long_edges_required,
				    DISTRIBUTION& fanout_degrees_that_need_long_edges) const;

	// Sanity check
	void sanity_check() const;
	void final_sanity_check();
	void check_degree_nNodes_match() const;


	NUM_ELEMENTS report_solution_quality() const;
	void cleanup();

};

#endif
