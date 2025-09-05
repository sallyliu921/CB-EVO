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





#ifndef statistic_reporter_H
#define statistic_reporter_H

#include "circ.h"
#include "circuit.h"
#include "degree_info.h"
#include <fstream>

//
// Class_name STATISTIC_REPORTER
//
// Description
//
//	Reports characterization statistics
//

class STATISTIC_REPORTER
{
public:
	STATISTIC_REPORTER();
	STATISTIC_REPORTER(const STATISTIC_REPORTER & another_statistic_reporter);
	STATISTIC_REPORTER & operator=(const STATISTIC_REPORTER & another_statistic_reporter);
	~STATISTIC_REPORTER();

	void report_stats(CIRCUIT * circuit);
private:
	CIRCUIT * 		m_circuit;
	DEGREE_INFO * 	m_degree_info;

	fstream m_output_file;

	void report_global_stats();
	void report_by_cluster_statistics();
	void report_degree_information(DEGREE_INFO * degree_info);
	void report_shape_information(SEQUENTIAL_LEVELS & sequential_levels, DEGREE_INFO * degree_info);
	void report_level_shape(SEQUENTIAL_LEVEL * seq_level, DEGREE_INFO * degree_info);
	void report_cluster_stastistics();
	void report_cluster(CLUSTER * cluster);
	void report_cluster_size_distribution(CLUSTER * cluster);

	void output_shape(const SHAPE & shape);
	void output_distribution(const DISTRIBUTION & distribution);

	void report_inter_cluster_adjacency_matrix();

	void report_delay_defining_stats(CLUSTER * cluster);


	void update_stats_regarding_edges_that_define_delay(NODE * node, NUM_ELEMENTS & number_of_ties,
	    		NUM_ELEMENTS & number_of_inter_cluster_edges, NUM_ELEMENTS & number_of_internal_edges);

	void report_edge_length_by_delay_level(SEQUENTIAL_LEVEL * seq_level);

	void report_reconvergence(CIRCUIT * circuit);
};


#endif
