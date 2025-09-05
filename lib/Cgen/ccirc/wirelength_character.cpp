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





#include "wirelength_character.h"
#include <math.h>
#include "rand.h"


WIRELENGTH_CHARACTER::WIRELENGTH_CHARACTER()
{
	long seed = 0;
	m_rand_number_gen = new RANDOM_NUMBER(seed);
}
WIRELENGTH_CHARACTER::WIRELENGTH_CHARACTER(const WIRELENGTH_CHARACTER & another_wirelength_character)
{
	m_circuit	= another_wirelength_character.m_circuit;
}

WIRELENGTH_CHARACTER & WIRELENGTH_CHARACTER::operator=(const WIRELENGTH_CHARACTER & another_wirelength_character)
{
	m_circuit	= another_wirelength_character.m_circuit;

	return (*this);
}

WIRELENGTH_CHARACTER::~WIRELENGTH_CHARACTER()
{
	delete m_rand_number_gen;
}

// Calculate the Wirelength-approx for the circuit
//
// PRE: Circuit is valid. Clusters have been created.
// POST: The Wirelength-approx has be calculated and set for the circuit
// RETURN: Wirelength-approx of the circuit
//
COST_TYPE  WIRELENGTH_CHARACTER::get_circuit_wirelength_approx
(
	CIRCUIT * circuit
)
{
	assert(circuit);
	m_circuit = circuit;
	double wirelength = 0;


	set_horizontal_position_of_nodes();

	wirelength = get_wirelength_measurement();

	iterate();

	wirelength = get_wirelength_measurement();

	m_circuit->set_wirelength_approx(wirelength);

	//print_wirelength_results();

	return wirelength;
}

// Sets the initial position of the nodes
//
// PRE: Circuit is valid. Clusters have been created.
// POST: The horizontal position of all nodes in the circuit have been set
//
void WIRELENGTH_CHARACTER::set_horizontal_position_of_nodes()
{
	m_width = m_circuit->get_max_width();
	SEQUENTIAL_LEVEL* cluster_seq_level = 0;

	CLUSTERS& clusters =  m_circuit->get_clusters();
	CLUSTERS::iterator cluster_iter;
	CLUSTER * cluster = 0;
	DELAY_LEVELS delay_levels;
	DELAY_LEVEL delay_level;
	NODE * node = 0;
	NODES::const_iterator node_iter;
	NUM_ELEMENTS index = 0,
				 nNodes = 0,
				 horizontal_position = 0;
	PORTS		 PI; 
	PORTS::iterator port_iter;
	PORT * pi = 0;
	DELAY_TYPE delay = 0,
			   max_delay = m_circuit->get_maximum_combinational_delay();
	NUM_ELEMENTS avg = 0,
				 central_offset = 0;


	// set the horizontal positions cluster by cluster
	for (cluster_iter = clusters.begin(); cluster_iter != clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);

		cluster_seq_level = cluster->get_sequential_level();
		assert(cluster_seq_level);

		delay_levels = cluster_seq_level->get_delay_levels();

		assert(static_cast<unsigned>(max_delay) == delay_levels.size()-1);


		// set the horizontal position at each delay level from 1 to the max_delay
		for (delay=1; delay <= max_delay; delay++)
		{
			delay_level = delay_levels[delay];


			nNodes = static_cast<NUM_ELEMENTS>(delay_level.size());

			// we want to spread the nodes out across
			//
			if (nNodes > 0)
			{
				avg = 16/2 * m_width * (nNodes - 1)/nNodes;
				central_offset = 16/2 * m_width - avg;
				assert(central_offset >= 0);
			}

			for (index = 0; index < nNodes; index++)
			{
				node = delay_level[index];
				assert(node);

				horizontal_position = (16 * m_width * index) / nNodes + central_offset;
				assert(0 <= horizontal_position && horizontal_position <= 16 * m_width);

				node->set_horizontal_position(horizontal_position);
			}
		}

		//
		// do the horizontal placement of the 0th delay level separately
		//

		// PI. should they be in the scheme of horizontal placement with nodes
		// or should they have a separate placement?
		assert(Dlook_at);

		PI = cluster->get_PI();
		delay_level = delay_levels[0];

		index = 0;
		nNodes = static_cast<NUM_ELEMENTS>(PI.size() + delay_level.size());

		if (nNodes > 0)
		{
			avg = 16/2 * m_width * (nNodes - 1)/nNodes;
			central_offset = 16/2 * m_width - avg;
			assert(central_offset >= 0);
		}

		// set the horizontal position of all primary inputs
		for (port_iter = PI.begin(); port_iter != PI.end(); port_iter++, index++)
		{
			pi = *port_iter;
			assert(pi);

			horizontal_position = (16 * m_width * index) / nNodes + central_offset;
			assert(0 <= horizontal_position && horizontal_position < 16 * m_width);

			pi->set_horizontal_position(horizontal_position);
		}
		// set the horizontal position of all flip-flops at the 0th delay level
		for (node_iter = delay_level.begin(); node_iter != delay_level.end(); node_iter++, index++)
		{
			node = *node_iter;
			assert(node);

			horizontal_position = (16 * m_width * index) / nNodes + central_offset;
			assert(0 <= horizontal_position && horizontal_position < 16 * m_width);

			node->set_horizontal_position(horizontal_position);
		}
	}
}


// gets the unnormalized Wirelength-approx
//
// PRE: nodes and primary inputs have their horizontal position defined
// RETURNS: unnormalized Wirelengthapprox
//
COST_TYPE WIRELENGTH_CHARACTER::get_wirelength_measurement() const
{
	CLUSTERS& clusters =  m_circuit->get_clusters();
	CLUSTERS::const_iterator cluster_iter;
	CLUSTER * cluster = 0;

	NUM_ELEMENTS cluster_number = 0;
	double  wirelength = 0.0,
			total_wirelength = 0.0;

	for (cluster_iter = clusters.begin(); cluster_iter != clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);

		wirelength = cluster->get_wirelength_approx();

		//debug("Cluster " << cluster_number << " has wirelength " << wirelength);

		total_wirelength += wirelength;

		cluster_number++;
	}
	//debug("Normalized by width wirelength " << total_wirelength/static_cast<double>(m_width));
	//debug("Normalized by width and edges wirelength " << 
	//total_wirelength/static_cast<double>(16*m_width * 
	//static_cast<double>(m_circuit->get_nEdges_without_clock_edges())));

	return total_wirelength;
}
// performs a simulated anneal to try and lower the wirelengthapprox measured
//
// PRE: nodes and primary inputs have their horizontal position defined
// POST: the wirelengthapprox is hopefully lower
//
void WIRELENGTH_CHARACTER::iterate()
{
	m_nodes = m_circuit->get_nodes();
	m_PI = m_circuit->get_PI();
	m_clusters = m_circuit->get_clusters();

	//print_wirelength_results();

	double cost = get_cost(),
				 lowest_cost = cost,
				 changed_cost = 0;

	NUM_ELEMENTS loops_without_display_of_lowest_cost = 0,
				 inner_loops = 1,
				 loops = 1;
	const NUM_ELEMENTS final_iteratation_limit = static_cast<NUM_ELEMENTS>(m_nodes.size())*1000;
	const NUM_ELEMENTS inner_loop_limit = static_cast<NUM_ELEMENTS>(m_nodes.size());

	double random_number = 0.0;
	double temperature = 0.01;
	double success_rate = 0.0;
	NUM_ELEMENTS number_success = 0;


	debug("Initial cost is " << cost << endl);

	while (cost > 0 && loops < final_iteratation_limit)
	{
		while (inner_loops < inner_loop_limit)
		{
			generate_move();

			if (m_move.is_valid())
			{
				changed_cost = get_changed_cost(cost);

				//debug("Changed cost is " << changed_cost);

				random_number = m_rand_number_gen->random_double_number(1);

				if (changed_cost <= 0.0 || random_number < exp(- changed_cost/temperature))
				{
					//debug("Making the move. Changed cost " << changed_cost);
					number_success++;

					cost += changed_cost;

					m_move.make_and_commit_move();

					//debug("New cost is " << cost);

					if (cost < lowest_cost)
					{

						lowest_cost = cost;
						
						loops_without_display_of_lowest_cost++;
						if (loops_without_display_of_lowest_cost == 500)
						{
							cost = get_cost();

							assert(m_move.get_port_a() == 0 && m_move.get_port_b() == 0 && 
									m_move.get_node_a() == 0 && m_move.get_node_b() == 0);
							lowest_cost = cost;
							loops_without_display_of_lowest_cost = 0;
							debug("************************************* Lowest cost " << lowest_cost << " **********");
							debug("*************************************************** Wirelengthapprox " << get_wirelength_measurement() << " **********");
						}
					}
				}
			}
				
			inner_loops++;
			loops++;
		}


		success_rate = static_cast<double>(number_success)/static_cast<double>(inner_loop_limit);	

		if (temperature > 1e-32)
		{
			update_temperature(temperature, success_rate, loops);
		}
		number_success = 0;
		inner_loops = 0;
	}

	debug("********************** Final lowest cost " << lowest_cost << " **********");
	debug("********************** Final Wirelengthapprox " << get_wirelength_measurement() << " **********");
}


// generate a move 
// 
// PRE: nothing
// POST: m_move contains a valid move
//
void WIRELENGTH_CHARACTER::generate_move()
{
	m_move.clear();

	NUM_ELEMENTS choice = 0,
				 size = 0;
	PORT * pi_a = 0,
		 * pi_b = 0;
	NODE * node_a = 0,
		 * node_b = 0;
	NUM_ELEMENTS nNodes = static_cast<NUM_ELEMENTS>(m_nodes.size());
	NUM_ELEMENTS nPI = static_cast<NUM_ELEMENTS>(m_PI.size());
	PORTS PIs;

	CLUSTER_NUMBER_TYPE cluster_number = 0;
	DELAY_TYPE delay = 0;
	CLUSTER * cluster = 0;
	SEQUENTIAL_LEVEL * seq_level = 0;
	DELAY_LEVELS delay_levels;
	DELAY_LEVEL delay_level;

	size = nNodes + nPI;

	choice = m_rand_number_gen->random_number(size-1);

	assert(m_move.get_port_a() == 0 && m_move.get_port_b() == 0 && m_move.get_node_a() == 0 && 
			m_move.get_node_b() == 0);

	if (choice >= nNodes)
	{
		// we have a primary input

		choice -= nNodes;
		assert(choice >= 0 && static_cast<unsigned>(choice) < m_PI.size());
		pi_a = m_PI[choice];
		assert(pi_a);
		cluster_number = pi_a->get_cluster_number();
		delay = 0;

		m_move.set_port_a(pi_a);
	}
	else
	{
		// we have only nodes at our delay level
		assert(choice < nNodes);

		node_a = m_nodes[choice];
		assert(node_a);

		cluster_number = node_a->get_cluster_number();
		delay = node_a->get_max_comb_delay_level();

		m_move.set_node_a(node_a);
	}



	assert(cluster_number >= 0 && static_cast<unsigned>(cluster_number) < m_clusters.size());

	cluster = m_clusters[cluster_number];
	assert(cluster);

	seq_level = cluster->get_sequential_level();
	assert(seq_level);

	delay_levels = seq_level->get_delay_levels();
	assert(delay >= 0 && static_cast<unsigned>(delay) < delay_levels.size());
	delay_level = delay_levels[delay];
	nNodes = static_cast<unsigned>(delay_level.size());


	if (delay > 0)
	{
		choice = m_rand_number_gen->random_number(nNodes-1);
		node_b = delay_level[choice];
		assert(node_b);

		m_move.set_node_b(node_b);
	}
	else
	{
		assert(delay == 0);

		size = nNodes + cluster->get_nPI();
		choice = m_rand_number_gen->random_number(size-1);

			
		// a choice greater than the number of nodes indicates a PI
		if (choice >= nNodes)
		{
			// we have a primary input
				
			choice -= nNodes;
			
			PIs = cluster->get_PI();


			assert(choice >= 0 && static_cast<unsigned>(choice) < PIs.size());
			pi_b = PIs[choice];
			assert(pi_b);

			m_move.set_port_b(pi_b);
		}
		else
		{
			// we have a node
	
			node_b = delay_level[choice];
			assert(node_b);

			m_move.set_node_b(node_b);
		}
	}
}

// 
// get the cost for the iterative improvement algorithm
//
// PRE: horizontal position of nodes has been defined
// RETURNS: the normalized cost
//
COST_TYPE WIRELENGTH_CHARACTER::get_cost() const
{
	NUM_ELEMENTS nEdges = m_circuit->get_nEdges_without_clock_edges();
	COST_TYPE wirelength = get_wirelength_measurement()/(static_cast<double>(m_width) * nEdges * 16);
	return wirelength;
}

// get the changed cost for the move being considered
//
// PRE: m_move contains a valid move
// RETURNS: the normalized changed cost of the move
//
COST_TYPE WIRELENGTH_CHARACTER::get_changed_cost(const double& old_cost)
{
	NUM_ELEMENTS nEdges = m_circuit->get_nEdges_without_clock_edges();

	COST_TYPE changed_cost = m_move.get_wirelength_change_of_nodes_involved();
	changed_cost = changed_cost/(static_cast<double>(m_width) * nEdges * 16);
	
	return changed_cost;

	/* if you want to get the exact cost
	m_move.make_temporary_move();
	COST_TYPE changed_cost = get_cost();
	m_move.unmake_temporary_move();
	changed_cost -= old_cost;
	*/
}
WIRELENGTH_MOVE::WIRELENGTH_MOVE()
{
	m_node_a 	= 0;
	m_node_b 	= 0;
	m_port_a 	= 0;
	m_port_b 	= 0;
}

WIRELENGTH_MOVE::~WIRELENGTH_MOVE()
{
}

//
// RETURNS: whether or not we have a valid move
//
bool WIRELENGTH_MOVE::is_valid() const
{
	bool valid = (m_node_a != m_node_b || m_port_a != m_port_b);

	return valid;
}

void WIRELENGTH_MOVE::make_temporary_move()
{
	make_move();
}
void WIRELENGTH_MOVE::unmake_temporary_move()
{
	make_move();
}

void WIRELENGTH_MOVE::make_and_commit_move()
{
	make_move();

	clear();
}

//
// clear the move
//
void WIRELENGTH_MOVE::clear()
{
	m_node_a 	= 0;
	m_node_b 	= 0;
	m_port_a 	= 0;
	m_port_b 	= 0;
}


//
// Swap the horizontal position of items selected for movement
//
void WIRELENGTH_MOVE::make_move()
{
	NUM_ELEMENTS horizontal_position_a = 0,
				 horizontal_position_b = 0;

	if (m_node_a && m_node_b)
	{
		horizontal_position_a = m_node_a->get_horizontal_position();
		horizontal_position_b = m_node_b->get_horizontal_position();

		m_node_a->set_horizontal_position(horizontal_position_b);
		m_node_b->set_horizontal_position(horizontal_position_a);
	}
	else if (m_node_a && m_port_b)
	{
		horizontal_position_a = m_node_a->get_horizontal_position();
		horizontal_position_b = m_port_b->get_horizontal_position();

		m_node_a->set_horizontal_position(horizontal_position_b);
		m_port_b->set_horizontal_position(horizontal_position_a);
	}
	else if (m_port_a && m_node_b)
	{
		horizontal_position_a = m_port_a->get_horizontal_position();
		horizontal_position_b = m_node_b->get_horizontal_position();

		m_port_a->set_horizontal_position(horizontal_position_b);
		m_node_b->set_horizontal_position(horizontal_position_a);
	}
	else if (m_port_a && m_port_b)
	{
		horizontal_position_a = m_port_a->get_horizontal_position();
		horizontal_position_b = m_port_b->get_horizontal_position();

		m_port_a->set_horizontal_position(horizontal_position_b);
		m_port_b->set_horizontal_position(horizontal_position_a);
	}
	else
	{
		Fail("Something is weird with a wirelength move");
	}

	assert(horizontal_position_a != horizontal_position_b);
}


//
// PRE: the move is valid
// RETURN: the unnormalized cost all the nodes and/or primary inputs affected by the move
//         
COST_TYPE WIRELENGTH_MOVE::get_wirelength_cost_of_nodes_involved()
{
	COST_TYPE wirelength_cost = 0;

	if (m_node_a && m_node_b)
	{
		wirelength_cost = get_wirelength_cost_of_node_and_fanout(m_node_a) + 
						get_wirelength_cost_of_node_and_fanout(m_node_b);
	}
	else if (m_node_a && m_port_b)
	{
		wirelength_cost = get_wirelength_cost_of_node_and_fanout(m_node_a) + 
						get_wirelength_cost_of_fanout_of_PI(m_port_b);
	}
	else if (m_port_a && m_node_b)
	{
		wirelength_cost = get_wirelength_cost_of_node_and_fanout(m_node_b) + 
						get_wirelength_cost_of_fanout_of_PI(m_port_a);
	}
	else if (m_port_a && m_port_b)
	{
		wirelength_cost = get_wirelength_cost_of_fanout_of_PI(m_port_a) + 
						get_wirelength_cost_of_fanout_of_PI(m_port_b);
	}
	else
	{
		Fail("Something is weird with a wirelength move");
	}

	return wirelength_cost;
}
//
// PRE: the node is valid
// RETURNS: the wirelength cost of the node and its fanout
//
COST_TYPE WIRELENGTH_MOVE::get_wirelength_cost_of_node_and_fanout
(
	NODE * node
) 
{
	assert(node);

	PORT * port = 0;
	NODE * sink_node = 0;
	EDGES output_edges;
	EDGES::const_iterator edge_iter;
	EDGE * edge = 0;

	// 1st get the wirelength cost of the node
	COST_TYPE wirelength = node->get_wirelength_approx_cost();

	// 2nd get the wirelength cost of the fanout of the node 
	// because those sink node's will also have had their wirelength change
	// by the movement of this node
	port = node->get_output_port();
	output_edges = port->get_edges();

	for (edge_iter = output_edges.begin(); edge_iter != output_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		sink_node = edge->get_sink_node();
		assert(node);

		wirelength += sink_node->get_wirelength_approx_cost();
	}

	return wirelength;
}

//
// PRE: the port is valid
// RETURNS: the wirelength cost of the primary output and its fanout
//
COST_TYPE WIRELENGTH_MOVE::get_wirelength_cost_of_fanout_of_PI
(
	const PORT * port
) const
{
	assert(port);

	COST_TYPE wirelength = 0.0;

	NODE * sink_node = 0;
	EDGES output_edges;
	EDGES::const_iterator edge_iter;
	EDGE * edge = 0;

	output_edges = port->get_edges();

	for (edge_iter = output_edges.begin(); edge_iter != output_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		sink_node = edge->get_sink_node();
		assert(sink_node);

		wirelength += sink_node->get_wirelength_approx_cost();
	}

	return wirelength;
}

//
//  Find the change in wirelength cost for the nodes involved in
//  the move
//
//  PRE: a wirelength move has been defined
//  RETURNS: the change in wirelength cost
//
COST_TYPE WIRELENGTH_MOVE::get_wirelength_change_of_nodes_involved()
{
	make_temporary_move();

	COST_TYPE changed_cost = get_wirelength_cost_of_nodes_involved();

	unmake_temporary_move();

	changed_cost -= get_wirelength_cost_of_nodes_involved();

	return changed_cost;
}




//
//  Find the change in wirelength cost normalized by the 
//  number of edges and the max_index
//
//  PRE: a wirelength move has been defined
//  RETURNS: the normalized changed in wirelength cost
//
COST_TYPE WIRELENGTH_MOVE::get_changed_wirelength_cost
(
	const double& wirelength_cost, 
	const double& max_index,
	const double& nEdges
)
{

	double changed_wirelength = 0.0,
			  changed_wirelength_of_sink_nodes = 0.0,
			  changed_cost = 0.0;

	changed_wirelength_of_sink_nodes = get_wirelength_change_of_nodes_involved();
		
	changed_cost = changed_wirelength/(max_index*nEdges);

	//debug("wirelength Changed cost " << changed_cost);
	
	return changed_cost;
}

//
// updates the temperature
// 
void WIRELENGTH_CHARACTER::update_temperature
(
	double & temperature,
	const double& success_rate,
	const NUM_ELEMENTS& loops
)
{
	if (success_rate > 0.96) 
	{
	   temperature = temperature * 0.5; 
	}
	else if (success_rate > 0.8) 
	{
	   temperature = temperature * 0.9;
	}
	else if (success_rate > 0.15)
	{
	   temperature = temperature * 0.95;
	}
	else if (success_rate > 0.2)
	{
	   //temperature = temperature * 0.8; 
	   temperature = temperature * 0.95; 
	}
	else
	{
	   temperature = temperature * 0.8; 
	}


	if (loops % 100)
	{
		debug("success_rate " << success_rate << "\tNew temperature " << temperature);
	}
}

//
// Prints out the list of horizontal positions of all nodes and primary inputs
//
void WIRELENGTH_CHARACTER::print_wirelength_results() const
{
	NODES::const_iterator node_iter;
	NODE * node = 0;
	PORTS::const_iterator port_iter;
	PORT * port = 0;

	for (port_iter = m_PI.begin(); port_iter != m_PI.end(); port_iter++)
	{
		port = *port_iter;
		assert(port);

		debug(port->get_name() << "\t\t" << " Cluster " << port->get_cluster_number() << "\t\t" << port->get_horizontal_position());
	}

	for (node_iter = m_nodes.begin(); node_iter != m_nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		debug(node->get_name() << "\t\t" << " Cluster " << node->get_cluster_number() << "\t\t" << node->get_horizontal_position());
	}
}
