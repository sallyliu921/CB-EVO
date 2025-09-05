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





#include "degree_info.h"
#include <math.h>

DEGREE_INFO::DEGREE_INFO()
{
	m_maximum_fanout		= 0;
	m_avg_fanin_for_comb	= 0;
	m_avg_fanout			= 0;
	m_avg_fanout_for_comb	= 0;
	m_avg_fanout_for_dff	= 0;
	m_avg_fanout_for_pi		= 0;
	m_std_dev_pi_fanout		= 0.0;
	m_std_dev_comb_fanout	= 0.0;
	m_std_dev_dff_fanout	= 0.0;
	m_std_dev_fanout		= 0.0;
	m_std_dev_comb_fanin	= 0.0;

	m_high_degree_comb		= 0;
	m_high_degree_pi		= 0;
	m_high_degree_dff		= 0;
	m_10plus_degree_comb	= 0;
	m_10plus_degree_pi		= 0;
	m_10plus_degree_dff		= 0;

    m_number_of_nodes		= 0.0;
    m_number_of_pi			= 0.0;
    m_number_of_comb		= 0.0;
    m_number_of_dff			= 0.0; 

	m_dff_exist				= false;
}
DEGREE_INFO::DEGREE_INFO(const DEGREE_INFO & another_degree_info)
{
	m_maximum_fanout		= another_degree_info.m_maximum_fanout;
	m_avg_fanin_for_comb	= another_degree_info.m_avg_fanin_for_comb;
	m_avg_fanout			= another_degree_info.m_avg_fanout;
	m_avg_fanout_for_comb	= another_degree_info.m_avg_fanout_for_comb;
	m_avg_fanout_for_dff	= another_degree_info.m_avg_fanout_for_dff;
	m_avg_fanout_for_pi		= another_degree_info.m_avg_fanout_for_pi;
	m_std_dev_pi_fanout		= another_degree_info.m_std_dev_pi_fanout;
	m_std_dev_comb_fanout	= another_degree_info.m_std_dev_comb_fanout;
	m_std_dev_dff_fanout	= another_degree_info.m_std_dev_dff_fanout;
	m_std_dev_fanout		= another_degree_info.m_std_dev_fanout;
	m_std_dev_comb_fanin	= another_degree_info.m_std_dev_comb_fanin;
	m_high_degree_comb		= another_degree_info.m_high_degree_comb;
	m_high_degree_pi		= another_degree_info.m_high_degree_pi;
	m_high_degree_dff		= another_degree_info.m_high_degree_dff;
	m_10plus_degree_comb	= another_degree_info.m_10plus_degree_comb;
	m_10plus_degree_pi		= another_degree_info.m_10plus_degree_pi;
	m_10plus_degree_dff		= another_degree_info.m_10plus_degree_dff;
	m_std_dev_comb_fanin	= another_degree_info.m_std_dev_comb_fanin;

	m_std_dev_pi_fanout		= another_degree_info.m_std_dev_pi_fanout;
	m_std_dev_comb_fanout   = another_degree_info.m_std_dev_comb_fanout;
	m_std_dev_dff_fanout    = another_degree_info.m_std_dev_dff_fanout;
	m_std_dev_fanout		= another_degree_info.m_std_dev_fanout;
}

DEGREE_INFO & DEGREE_INFO::operator=(const DEGREE_INFO & another_degree_info)
{
	m_maximum_fanout		= another_degree_info.m_maximum_fanout;
	m_avg_fanin_for_comb	= another_degree_info.m_avg_fanin_for_comb;
	m_avg_fanout			= another_degree_info.m_avg_fanout;
	m_avg_fanout_for_comb	= another_degree_info.m_avg_fanout_for_comb;
	m_avg_fanout_for_dff	= another_degree_info.m_avg_fanout_for_dff;
	m_avg_fanout_for_pi		= another_degree_info.m_avg_fanout_for_pi;
	m_high_degree_comb		= another_degree_info.m_high_degree_comb;
	m_high_degree_pi		= another_degree_info.m_high_degree_pi;
	m_high_degree_dff		= another_degree_info.m_high_degree_dff;
	m_std_dev_pi_fanout		= another_degree_info.m_std_dev_pi_fanout;
	m_std_dev_comb_fanout	= another_degree_info.m_std_dev_comb_fanout;
	m_std_dev_dff_fanout	= another_degree_info.m_std_dev_dff_fanout;
	m_std_dev_fanout		= another_degree_info.m_std_dev_fanout;
	m_std_dev_comb_fanin	= another_degree_info.m_std_dev_comb_fanin;
	m_10plus_degree_comb	= another_degree_info.m_10plus_degree_comb;
	m_10plus_degree_pi		= another_degree_info.m_10plus_degree_pi;
	m_10plus_degree_dff		= another_degree_info.m_10plus_degree_dff;
	m_std_dev_pi_fanout		= another_degree_info.m_std_dev_pi_fanout;
	m_std_dev_comb_fanout   = another_degree_info.m_std_dev_comb_fanout;
	m_std_dev_dff_fanout    = another_degree_info.m_std_dev_dff_fanout;
	m_std_dev_fanout		= another_degree_info.m_std_dev_fanout;

	return (*this);
}

DEGREE_INFO::~DEGREE_INFO()
{
}

//
// Calcuate the degree information for the total circuit 
// that will be put into a .stats file 
//
// PRE: circuit is valid
// POST: info has been calculated 
//
void DEGREE_INFO::calculate_degree_information_for_circuit
(
	CIRCUIT * circuit
)
{
	assert(circuit);

	NODES & nodes = circuit->get_nodes();
	PORTS pi = circuit->get_PI();

	m_number_of_pi		= static_cast<double>(pi.size());
	m_number_of_nodes	= static_cast<double>(circuit->get_nNodes());
	m_number_of_comb	= static_cast<double>(circuit->get_nComb());
	m_number_of_dff		= static_cast<double>(circuit->get_nDFF());

	m_dff_exist = (circuit->get_nDFF() > 0);

	calculate_degree_information(nodes, pi);
}

//
// Calculates the degree information for the cluster
// that will be put into a .stats file
//
// PRE: cluster is valid
// POST: info has been calculated
//  
void DEGREE_INFO::calculate_degree_information_for_cluster
(
	CLUSTER * cluster
)
{
	assert(cluster);

	NODES & nodes 	= cluster->get_nodes();
	PORTS pi		  		= cluster->get_PI();

	m_number_of_pi		= static_cast<double>(cluster->get_nPI());
	m_number_of_nodes	= static_cast<double>(cluster->get_nNodes());
	m_number_of_comb	= static_cast<double>(cluster->get_nComb());
	m_number_of_dff		= static_cast<double>(cluster->get_nDFF());

	m_dff_exist			= (cluster->get_nDFF() > 0);

	calculate_degree_information(nodes, pi);
}
//
// Calculates the degree information for nodes and pi 
//
// PRE: nodes and pi parameters contain the things we want to calculate for
// POST: inform has been calculated
//  
void DEGREE_INFO::calculate_degree_information
(
	NODES & nodes,
	PORTS & pi
)
{
	NUM_ELEMENTS total_comb_fanin	= 0;
	NUM_ELEMENTS total_comb_fanout 	= 0;
	NUM_ELEMENTS total_dff_fanout	= 0;
	NUM_ELEMENTS total_pi_fanout	= 0;

	sum_node_totals(nodes, total_comb_fanout, total_dff_fanout, total_comb_fanin);
	sum_pi_totals(pi, total_pi_fanout); 

	calculate_averages(total_comb_fanout, total_dff_fanout, total_pi_fanout, total_comb_fanin);

	calculate_std_deviations(nodes, pi);
	find_high_degree_fanout_nodes(nodes);
	find_high_degree_fanout_pi(pi);
}


// Finds fanin and fanout sums
// 
// PRE: nodes contains the nodes that we want to calculate for.
// POST: total_comb_fanout,total_dff_fanout,total_comb_fanin have been calculated
//
void DEGREE_INFO::sum_node_totals
(
	const NODES & nodes,
	NUM_ELEMENTS & total_comb_fanout,
	NUM_ELEMENTS & total_dff_fanout,
	NUM_ELEMENTS & total_comb_fanin
)
{
	NODES::const_iterator node_iter;
	NODE * node = 0;


	NUM_ELEMENTS fanout_degree = 0;

	for (node_iter = nodes.begin(); node_iter != nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		fanout_degree = node->get_fanout_degree();

		m_maximum_fanout = MAX(m_maximum_fanout, fanout_degree);

		if (node->get_type() == NODE::COMB)
		{
			total_comb_fanout += fanout_degree;
			total_comb_fanin  += node->get_fanin_degree();
		}
		else
		{
			assert(node->get_type() == NODE::SEQ);

			total_dff_fanout += fanout_degree;
		}
	}
}

// Finds the total pi fanout
// 
// PRE: pi contains the PIs that we want to calculate for.
// POST: total_pi_fanout has been calculated
//
void DEGREE_INFO::sum_pi_totals
(
	const PORTS & pi,
	NUM_ELEMENTS & total_pi_fanout
)
{
	PORTS::const_iterator port_iter;
	PORT * port = 0;
	NUM_ELEMENTS pi_fanout = 0;

	for (port_iter = pi.begin(); port_iter != pi.end(); port_iter++)
	{
		port = *port_iter;
		assert(port);

		if (port->get_io_direction() != PORT::CLOCK)
		{
			pi_fanout = port->get_fanout_degree();
			total_pi_fanout += pi_fanout;

			m_maximum_fanout = MAX(m_maximum_fanout, pi_fanout);
		}
	}
}


//
// Calculates the fanin/fanout averages
//
// PRE: total_comb_fanout, total_dff_fanout, 
//      total_pi_fanout, total_comb_fanin have been calculated
// POST: m_avg_fanin_for_comb,m_avg_fanout_for_comb,
//       m_avg_fanout_for_pi,m_avg_fanout 
//       have been calculated
void DEGREE_INFO::calculate_averages
(
	const NUM_ELEMENTS & total_comb_fanout,
	const NUM_ELEMENTS & total_dff_fanout,
	const NUM_ELEMENTS & total_pi_fanout,
	const NUM_ELEMENTS & total_comb_fanin
)
{
	if (m_dff_exist)
	{
		m_avg_fanout_for_dff	= static_cast<double>(total_dff_fanout)/m_number_of_dff;
	}
	
	m_avg_fanin_for_comb	= static_cast<double>(total_comb_fanin)/m_number_of_comb;
	m_avg_fanout_for_comb	= static_cast<double>(total_comb_fanout)/m_number_of_comb;
	m_avg_fanout_for_pi		= static_cast<double>(total_pi_fanout)/m_number_of_pi;

	m_avg_fanout			= static_cast<double>(total_comb_fanout + total_dff_fanout + total_pi_fanout)/
												(m_number_of_pi+m_number_of_nodes);

}

//
// Calculates the std. deviations of the fanin/fanout
//
// PRE: nodes and pi parameters contains the things we want to calculate for
// POST: std. dev. have been calculated
//
void DEGREE_INFO::calculate_std_deviations
(
	const NODES & nodes,
	const PORTS & pi
)
{

	NODES::const_iterator node_iter;
	NODE * node = 0;
	PORTS::const_iterator port_iter;
	PORT * port = 0;

	// decided to make the totals a double because i think it might be better to have 
	// avg - data_point to be in double instead of long
	double fanout_degree = 0.0;
	double total_sq_comb_fanout = 0.0;
	double total_sq_dff_fanout = 0.0;
	double total_sq_pi_fanout = 0.0;
	double total_sq_fanout = 0.0;
	double total_sq_fanin = 0.0;

	
	for (node_iter = nodes.begin(); node_iter != nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		fanout_degree = node->get_fanout_degree();

		if (node->get_type() == NODE::COMB)
		{
			total_sq_comb_fanout += square(m_avg_fanout_for_comb - fanout_degree);
			total_sq_fanin  += square(m_avg_fanin_for_comb - node->get_fanin_degree());
		}
		else
		{
			assert(node->get_type() == NODE::SEQ);
			total_sq_dff_fanout += square(m_avg_fanout_for_dff - fanout_degree);
		}
		
		total_sq_fanout += square(m_avg_fanout -  fanout_degree);

	}

	for (port_iter = pi.begin(); port_iter != pi.end(); port_iter++)
	{
		port = *port_iter;
		assert(port);

		if (port->get_io_direction() != PORT::CLOCK)
		{
			fanout_degree = port->get_fanout_degree();
			total_sq_pi_fanout += square(m_avg_fanout_for_pi -  fanout_degree);
			total_sq_fanout += square(m_avg_fanout - fanout_degree);
		}
	}

	if (m_dff_exist)
	{
		m_std_dev_dff_fanout	= sqrt(static_cast<double>(total_sq_dff_fanout)/m_number_of_dff);
	}


	m_std_dev_pi_fanout 	= sqrt(static_cast<double>(total_sq_pi_fanout)/m_number_of_pi);
	m_std_dev_comb_fanout 	= sqrt(static_cast<double>(total_sq_comb_fanout)/m_number_of_comb);

	m_std_dev_fanout		= sqrt(static_cast<double>(total_sq_fanout)/(m_number_of_pi+m_number_of_nodes));

	m_std_dev_comb_fanin	= sqrt(static_cast<double>(total_sq_fanin)/m_number_of_comb);
}

// 
// Calculate the number of nodes with fanout that is:
// a) Above 10
// b) One std. deviation above the average
//
// PRE: nodes contains the nodes we want to calculate for
// POST: m_10plus_degree_comb and m_high_degree_comb have been calculated
//
void DEGREE_INFO::find_high_degree_fanout_nodes
(
	const NODES & nodes
)
{
	NODES::const_iterator node_iter;
	NODE * node = 0;

	NUM_ELEMENTS fanout_degree = 0;

	for (node_iter = nodes.begin(); node_iter != nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		fanout_degree = node->get_fanout_degree();

		if (node->get_type() == NODE::COMB)
		{
			if (fanout_degree >= 10)
			{
				m_10plus_degree_comb++;
			}
			if (fanout_degree >= m_avg_fanout + m_std_dev_fanout)
			{
				m_high_degree_comb++;
			}

		}
		else
		{
			assert(node->get_type() == NODE::SEQ);

			if (fanout_degree >= 10)
			{
				m_10plus_degree_dff++;
			}
			if (fanout_degree >= m_avg_fanout + m_std_dev_fanout)
			{
				m_high_degree_dff++;
			}
		}

	}
}

// 
// Calculate the number of primary inputs with fanout that is:
// a) Above 10
// b) One std. deviation above the average
//
// PRE: pi contains the primary inputs we want to calculate for
// POST: m_10plus_degree_pi and m_high_degree_pi have been calculated
//
void DEGREE_INFO::find_high_degree_fanout_pi
(
	const PORTS & pi
)
{
	PORTS::const_iterator port_iter;
	PORT * port = 0;
	NUM_ELEMENTS pi_fanout = 0;

	for (port_iter = pi.begin(); port_iter != pi.end(); port_iter++)
	{
		port = *port_iter;
		assert(port);

		if (port->get_io_direction() != PORT::CLOCK)
		{
			pi_fanout = port->get_fanout_degree();

			if (pi_fanout >= 10)
			{
				m_10plus_degree_pi++;
			}
			if (pi_fanout >= m_avg_fanout + m_std_dev_fanout)
			{
				m_high_degree_pi++;
			}
		}
	}
}



// returns arg*arg
double DEGREE_INFO::square
(
	const double & arg
)
{
	return arg*arg;
}

