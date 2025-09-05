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



#ifndef degree_info_H
#define degree_info_H

class DEGREE_INFO;

#include "circ.h"
#include "circuit.h"
#include "cluster.h"

//
// Class_name DEGREE_INFO
//
// Description
//
//		Calculates the fanin/fanout degree for circuits and clusters.
//		Also serves as the repository of such information
//		


class DEGREE_INFO
{
public:
	DEGREE_INFO();
	DEGREE_INFO(const DEGREE_INFO & another_degree_info);
	DEGREE_INFO & operator=(const DEGREE_INFO & another_degree_info);
	~DEGREE_INFO();


	void calculate_degree_information_for_circuit(CIRCUIT * circuit);
	void calculate_degree_information_for_cluster(CLUSTER * cluster);

	NUM_ELEMENTS get_maximum_fanout_degree() const { return m_maximum_fanout; }

	double get_avg_fanin_for_comb() 	const { return m_avg_fanin_for_comb; }
	double get_avg_fanout() 			const { return m_avg_fanout; }
	double get_avg_fanout_for_comb()	const { return m_avg_fanout_for_comb; }
	double get_avg_fanout_for_dff() 	const { return m_avg_fanout_for_dff; }
	double get_avg_fanout_for_pi() 	const { return m_avg_fanout_for_pi; }

	NUM_ELEMENTS get_high_degree_comb() 	const { return m_high_degree_comb; }
	NUM_ELEMENTS get_high_degree_pi() 		const { return m_high_degree_pi; }
	NUM_ELEMENTS get_high_degree_dff() 		const { return m_high_degree_dff; }

	NUM_ELEMENTS get_10plus_fanout_degree_comb()const { return m_10plus_degree_comb; }
	NUM_ELEMENTS get_10plus_fanout_degree_pi() 	const { return m_10plus_degree_pi; }
	NUM_ELEMENTS get_10plus_fanout_degree_dff()	const { return m_10plus_degree_dff; }
	double get_std_dev_comb_fanin() const { return m_std_dev_comb_fanin; }
	double get_std_dev_pi_fanout() const { return m_std_dev_pi_fanout; }
	double get_std_dev_comb_fanout() const { return m_std_dev_comb_fanout; }
	double get_std_dev_dff_fanout() const { return m_std_dev_dff_fanout; }
	double get_std_dev_fanout() const { return m_std_dev_fanout; }

	DISTRIBUTION get_fanout_distribution() 	const { assert(false); return m_fanout_distribution; }
	DISTRIBUTION get_fanin_distribution() 	const { assert(false); return m_fanin_distribution; }

private:
	CIRCUIT * 		m_circuit;

	NUM_ELEMENTS 	m_maximum_fanout;

	double			m_avg_fanin_for_comb;
	double			m_avg_fanout;
	double			m_avg_fanout_for_comb;
	double			m_avg_fanout_for_dff;
	double			m_avg_fanout_for_pi;
	NUM_ELEMENTS	m_high_degree_comb;
	NUM_ELEMENTS	m_high_degree_pi;
	NUM_ELEMENTS	m_high_degree_dff;
	NUM_ELEMENTS	m_10plus_degree_comb;
	NUM_ELEMENTS	m_10plus_degree_pi;
	NUM_ELEMENTS	m_10plus_degree_dff;

	double			m_std_dev_comb_fanin;
	double			m_std_dev_pi_fanout;
	double			m_std_dev_comb_fanout;
	double			m_std_dev_dff_fanout;
	double			m_std_dev_fanout;

    double			m_number_of_nodes;
    double			m_number_of_pi;
    double			m_number_of_comb;
    double			m_number_of_dff; 

	bool   			m_dff_exist;

	DISTRIBUTION	m_fanout_distribution;
	DISTRIBUTION	m_fanin_distribution;


	void calculate_degree_information(NODES & list_of_nodes, PORTS & pi);

	void sum_node_totals(const NODES & list_of_nodes, NUM_ELEMENTS & total_comb_fanout, 
						NUM_ELEMENTS & total_dff_fanout, NUM_ELEMENTS & total_comb_fanin);
	void sum_pi_totals(const PORTS & pi, NUM_ELEMENTS & total_pi_fanout);
	void calculate_averages(const NUM_ELEMENTS & total_comb_fanin,
							const NUM_ELEMENTS & total_comb_fanout,
							const NUM_ELEMENTS & total_dff_fanout,
							const NUM_ELEMENTS & total_pi_fanout);
	void calculate_std_deviations(const NODES & list_of_nodes, const PORTS & pi);

	void find_high_degree_fanout_nodes(const NODES & list_of_nodes);
	void find_high_degree_fanout_pi(const PORTS & pi);

	
	double square(const double & arg);
};


#endif
