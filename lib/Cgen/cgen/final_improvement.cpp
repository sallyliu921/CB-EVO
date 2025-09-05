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



#include "final_improvement.h"
#include "util.h"
#include <algorithm>
#include <numeric>
#include <math.h>

SIMPLE_EDGE_ITERATOR::SIMPLE_EDGE_ITERATOR()
{
	assert(g_options);

	m_circuit = 0;
	m_display_costs = false;
	m_max_index = 0.0;
	m_eta = g_options->get_final_edge_assignment_eta();
	m_dff_loop_penalty_multipler = g_options->get_final_edge_assignment_dff_loop_penalty_multiplier();
	m_bad_node_cost = 0.0;
	m_wirelength_wanted = 0.0;

}
SIMPLE_EDGE_ITERATOR::~SIMPLE_EDGE_ITERATOR()
{
}

//  After creating the initial solution, we employ an iterative 
//  algorithm that selects certain edges as candidates for relocation 
//  and accepts or rejects proposed changes ("moves") based on a cost 
//  function. The cost function of the algorithm is as follows: 
//
//  cost = Eta*(Wirelength-approx - desired Wirelength-approx) +
//         Sum over all nodes n of (Number of Violations(n))
//
//  Edge Assign Here, desired wirelength is the desired Wirelength-approx
//  and Number of Violations is a function that returns the number of nodes 
//  that have no inputs, too many inputs, two or more connections from the 
//  same source node, or if the node is a flip-flop a connection to itself. 
//  The wirelength costs are normalized to the maximum horizontal position 
//  multiplied by the number of edges while the Number of Violations cost is 
//  normalized to the number of edges. We use the factor Eta to balance the 
//  goal of achieving the desired wirelength against the goal of having no 
//  node violations. The value is set to 0.01 because the wirelength 
//  cost is often much larger than the no node violation cost and while 
//  achieving the desired wirelength is important it is more important 
//  that we have no node violations because they create sizeable difficulties. 
//
//  We generate moves 95% of the time purely randomly while 5% of the time we 
//  target nodes with violations. When we generate a purely random move, we 
//  start by randomly selecting an edge in the graph. With this edge we attempt 
//  one of two possible move types with equal probability. Each move preserves 
//  the combinational delay structure and the number of edges that output from each node.
//
//	We employ a hill climbing iterator, in which all valid changes that 
//	improve the cost function are accepted, and valid bad moves are accepted 
//	with a probability that decreases exponentially with the change in cost. 
//	A move is valid if it will not create any flip-flops with loops back to 
//	themselves. The algorithm continues until the cost is zero or until the number 
//	of moves attempted is a hundred times the number of edges in the graph. 
//  
//  PRE: circuit is valid
//  POST: we have attempted to improve the cost
//
void SIMPLE_EDGE_ITERATOR::iterate
(
	CIRCUIT_GRAPH * circuit
)
{
	assert(circuit);
	m_circuit = circuit;
	m_move.set_circuit(m_circuit);
	m_clusters = m_circuit->get_clusters();
	m_max_index = static_cast<COST_TYPE>(m_circuit->get_max_index());
	m_edges = m_circuit->get_simple_edges();
	m_nEdges = static_cast<COST_TYPE>(m_edges.size());

	debug("Status: Trying to Improve Final Edge Assignment\n");

	create_edges_to_rectify_delay_structure_with_individual_edge_assignment();
	make_list_of_node_that_have_violations();

	m_circuit->sanity_check();
	dbg("\n");

	// choose the Wirelength-approx to target
	if (g_options->is_best_wirelength_settings() && m_circuit->is_combinational())
	{
		debug("The circuit is combinational.\n");
		debug("The best Wirelength-approx setting is to accept the initial wirelength.");
			
		// We will calculate the initial wirelength in a bit
	}
	else if (g_options->is_best_wirelength_settings() && m_circuit->is_sequential())
	{
		debug("The circuit is sequential.");
		debug("The best Wirelength-approx setting is to go for 0 wirelength");
		m_wirelength_wanted = 0;
	}
	else if (g_options->is_accept_initial_wirelength())
	{
		debug("Accepting initial wirelength.");
			
		// We will calculate the initial wirelength in a bit
	}
	else if (g_options->is_wirelength_from_stat_file())
	{
		m_wirelength_wanted = m_circuit->get_circuit_spec()->get_global_wirelength();
		debug("Wirelength from spec is " << m_wirelength_wanted);
	}
	else if (g_options->is_user_wirelength())
	{
		m_wirelength_wanted = g_options->get_wirelength_wanted();
		debug("Using the wirelength argumenent " << m_wirelength_wanted);
	}
	else
	{
		Fail("Unknown wirelength setting");
	}


	COST_TYPE cost = get_cost(),
				 lowest_cost = cost,
				 changed_cost = 0;

	NUM_ELEMENTS loops_without_display_of_lowest_cost = 0,
				 loops = 1;
	const NUM_ELEMENTS final_iteratation_limit = static_cast<NUM_ELEMENTS>(m_nEdges)*100;

	double random_number = 0.0;
	double temperature = g_options->get_final_edge_assignment_init_temperature();

	// set m_wirelength_wanted if we wanted the initial wirelength
	if (g_options->is_accept_initial_wirelength() || 
		(g_options->is_best_wirelength_settings() && m_circuit->is_combinational()))
	{
		m_wirelength_wanted = m_wirelength_cost;
		debug("The initial Wirelength-approx is " << m_wirelength_wanted);
	}
	dbg("\n\n");

	debug("The initial Wirelength-approx is " << m_wirelength_cost);
	debug("Initial node violation cost is " << m_bad_node_cost << endl);

	debug("Initial cost is " << cost);
	debug("Initial temperature " << temperature << endl);

	while (cost > 0 && loops < final_iteratation_limit)
	{
		generate_move();

		if (m_move.is_valid())
		{
			changed_cost = get_changed_cost();

			random_number = g_rand_number_gen->random_double_number(1);

			if (changed_cost <= 0.0 || random_number < exp(- changed_cost/temperature))
			{
				//debug("Making the move");

				cost += changed_cost;

				m_move.make_and_commit_move(m_nodes_with_violations, m_wirelength_cost, m_bad_node_cost);

				//cost += changed_cost;
				//debug("New cost is " << cost);

				if (cost < lowest_cost)
				{
					cost = get_cost();
					lowest_cost = cost;

					loops_without_display_of_lowest_cost++;
					if (loops_without_display_of_lowest_cost == 1000)
					{
						loops_without_display_of_lowest_cost = 0;
						debug("********************** Lowest cost " << lowest_cost << " **********");
						debug("********************** wirelength Lowest cost " << m_wirelength_cost << " **********");
					}
				}
			}
		}

		assert(new_gen);
		// need to experiment with the best value of beta
		// need to experiment with annealing schedule.
		//if (loops % iteratation_limit == 0)
		//{
		//	dbg("********* Decreasing beta from " << m_eta);
		//	m_eta = m_eta * 0.8;
		//	cost = get_cost();
		//}
		loops++;
	}

	debug("********************** Final lowest cost " << lowest_cost << " **********");
	debug("\n********************** Final node violation cost " << m_bad_node_cost << " **********");
	debug("********************** Final wirelength cost " << m_wirelength_cost << " **********");
	debug("********************** Final normalized wirelength cost " << m_wirelength_cost/(m_max_index*m_nEdges) << " **********");

	get_cost();

	
	debug("\nStatus: Finished trying to improve the final edge assignment\n");

	//m_circuit->sanity_check();

	m_display_costs = true;

	//print_node_positions();
	//print_wirelength_stats();
	

	update_cluster_to_reflect_which_simple_edges_they_contain();
}




// PRE: we possibly have nodes in the circuit that could not
//      make as many connections to sink nodes as they wanted to 
//      and therefore have "lost edges"
// POST: the nodes that wanted to make more connections have 
//       a second chance. for these nodes 
//       their edges_currently_assigning has been set to the number
//       of connections they were missing.
//       
void SIMPLE_EDGE_ITERATOR::add_lost_edges_back_to_nodes()
{
	NUM_ELEMENTS lost_edges = 0,
				 unused_long_and_inter_cluster_edges = 0;

	NODES& nodes = m_circuit->get_nodes();
	NODES::const_iterator node_iter;
	NODE * node = 0;

	for (node_iter = nodes.begin(); node_iter != nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		assert(node->get_edges_currently_assigning() == 0);
		assert(node->get_edges_to_assign() == 0);

		lost_edges = node->get_nLost_simple_edges();
		unused_long_and_inter_cluster_edges = node->get_unassigned_long_and_inter_cluster_edges();

		if (lost_edges > 0)
		{
			node->add_edges_currently_assigning(lost_edges);
		}
		if (unused_long_and_inter_cluster_edges > 0)
		{
			node->add_edges_currently_assigning(unused_long_and_inter_cluster_edges);
		}
	}
}


// generate a move that will try to fix a node violation
//
// For the 5% of moves that explicitly attempt to eliminate node violations, 
// we randomly select a node that is in violation. If the violation type 
// is either too many inputs or two or more connections from same source node, 
// we select one of the problem edges and attempt an edge rotation. If the 
// node violation is a flip-flop that connects to itself we attempt a double 
// edge swap where we choose the first edge that connects the flip-flop 
// to itself and the second edge from the list of all edges in the 
// graph that connect to flip-flops and whose choice will preserve the Latched specification.
//
// PRE: We have nodes with violations, m_nodes_with_violations.size() > 0
// POST: m_move has been set if we could make a move
//
void SIMPLE_EDGE_ITERATOR::bad_node_generate_move()
{

	NODE * sink_node = 0,
		 * source_node = 0,
		 * source_node2 = 0,
		 * node = 0;
	
	NUM_ELEMENTS nBad_nodes = static_cast<NUM_ELEMENTS>(m_nodes_with_violations.size()),
				 choice = 0,
				 edge_choice = 0,
			     nInputs = 0,
			     nOutputs = 0,
				 nNodes = 0;

	EDGES input_edges,
		  output_edges;

	EDGE * edge = 0;
	LEVEL_NODE * source_level_node = 0,
			   * sink_level_node = 0;

	assert(nBad_nodes != 0);

	// 
	// Step 1: Let's select a node with a node violation
	//
	
	// 1st choose the sink node with problems
	choice = g_rand_number_gen->random_number(nBad_nodes-1);
	assert(nBad_nodes > 0);
	assert(choice >= 0 && choice < nBad_nodes);

	node = m_nodes_with_violations[choice];
	assert(node);


	// our node does have a violation, right?
	if (! (node->has_double_inputs() || node->has_excess_inputs() || 
				node->needs_inputs() || node->is_dff_have_loop_back_to_itself()))
	{
		// no it doesn't.
		//
		// A node can be on the node violation list even though it does not have a
		// violation because of one special case where...
		// when we have a dff node loop violation 
		// and a different source node moves an edge over to the buffer 
		// that in effects resolves the dff loop
		// (although in truth dff->gate->same dff is not reall a good thing)
		// we don't remove the dff from the node violation list

		// we should resolve this problem more elegantly
		assert(new_gen);

		// remove the node from the list and return
		NODES::iterator node_iter = find(m_nodes_with_violations.begin(), m_nodes_with_violations.end(), node);
		m_nodes_with_violations.erase(node_iter);

		return;
	}

	// yes, we have a node with a violation
	assert(node->has_double_inputs() || node->has_excess_inputs() || node->needs_inputs()
				|| node->is_dff_have_loop_back_to_itself());


	// Is it a dff with a loop?
	if (node->is_dff_have_loop_back_to_itself())
	{
		// yes, generate a move for it and return
		generate_move_to_solve_dff_loop(node);
		return;
	}


	// Our node has problems with the edges that input into it
	// and is therefore a sink node
	sink_node = node;



	// Step 2: select an edge that we want to move
	if (sink_node->has_double_inputs())
	{
		assert(sink_node->is_combinational());

		// get the double input edges
		input_edges = sink_node->get_double_input_edges();
	}
	else if (sink_node->has_excess_inputs())
	{
		assert(sink_node->is_combinational());
		
		// get all the edges because none of them are double edges
		input_edges = sink_node->get_input_edges();
	}
	else
	{
		// don't make a move for this yet, just return
		assert(new_gen);
		assert(sink_node->needs_inputs());

		return;
	}
	nInputs = static_cast<NUM_ELEMENTS>(input_edges.size());

	assert(nInputs > 0);

	// make a choice of an edge
	edge_choice = g_rand_number_gen->random_number(nInputs-1);
	assert(edge_choice >= 0 && edge_choice < nInputs);

	edge = input_edges[edge_choice];
	assert(edge);
	
	m_move.set_edge_a(edge);


	// Step 3:
	// we now have our edge we want to move somewhere.
	// we have two choice for moves.
	// randomly pick a move type

	if (g_rand_number_gen->random_double_number(1) < 0.5 && sink_node->is_combinational())
	{
		// make a edge rotation
		// for the edge the source stays the same. the destination changes

		m_move.set_type(SIMPLE_EDGE_MOVE::EDGE_ROTATION);

		// pick a destination node in the sink level node to move this edge to
		sink_level_node = sink_node->get_my_level_node();
		assert(sink_level_node);

		nNodes = sink_level_node->get_nNodes();

		choice = g_rand_number_gen->random_number(nNodes-1);
		assert(choice >= 0 && choice < nNodes);

		sink_node = sink_level_node->get_node(choice);
		assert(sink_node);

		m_move.set_destination_sink_node(sink_node);
	}
	else
	{
		// make a edge swap.
		// two source nodes swap edges. 

		m_move.set_type(SIMPLE_EDGE_MOVE::EDGE_SWAP_BETWEEN_SOURCES);

		// source_node is the source of the sink node with a problem.
		source_node = edge->get_source_node();
		assert(source_node);

		source_level_node = source_node->get_my_level_node();
		assert(source_level_node);

		nNodes = source_level_node->get_nNodes();

		choice = g_rand_number_gen->random_number(nNodes-1);
		assert(choice >= 0 && choice < nNodes);

		source_node2 = source_level_node->get_node(choice);
		assert(source_node2);

		output_edges = source_node2->get_output_edges();
		nOutputs = static_cast<NUM_ELEMENTS>(output_edges.size());

		// while we have not found a node with at least one output edge
		// and there should be at least one because source_node has at least 1.
		// ie. we might choose a duplicate edge but at least it will be an edge
		while (nOutputs < 1)
		{
			choice = g_rand_number_gen->random_number(nNodes-1);
			assert(choice >= 0 && choice < nNodes);

			source_node2 = source_level_node->get_node(choice);
			assert(source_node2);

			output_edges = source_node2->get_output_edges();
			nOutputs = static_cast<NUM_ELEMENTS>(output_edges.size());
		}

		choice = g_rand_number_gen->random_number(nOutputs-1);
		assert(choice >= 0 && choice < nOutputs);

		edge = output_edges[choice];
		assert(edge);

		m_move.set_edge_b(edge);
	}
}

//  We generate moves 95% of the time purely randomly while 5% of the time we 
//  target nodes with violations. When we generate a purely random move, we 
//  start by randomly selecting an edge in the graph. With this edge we attempt 
//  one of two possible move types with equal probability. Each move preserves 
//  the combinational delay structure and the number of edges that output from each node.
//
//  In the first move type, defined as an edge rotation, we move the end 
//  point of the edge to a new sink node. The new sink node is randomly 
//  chosen from within the same level node as the old sink node. In 
//  the second move type, defined as a double edge swap, we select a 
//  second edge and move the end point of each edge to the other edge's 
//  sink node. The second edge is randomly chosen from the edges that output 
//  from the same level node that the first edge outputs from. 
//
//  For the 5% of moves that explicitly attempt to eliminate node violations, we 
//  randomly select a node that is in violation. 
// 
//  PRE: nothing
//  POST: m_move contains a move that may or may not be valid
void SIMPLE_EDGE_ITERATOR::generate_move()
{
	NODE * sink_node = 0,
		 * source_node = 0,
		 * source_node2 = 0;
	NUM_ELEMENTS nEdges = static_cast<NUM_ELEMENTS>(m_edges.size()),
				 choice = 0,
			     nOutputs = 0,
				 nNodes = 0;

	EDGES input_edges,
		  output_edges;

	EDGE * edge = 0;
	LEVEL_NODE * source_level_node = 0,
			   * sink_level_node = 0;
	double random_number = 0.0;

	SHAPE latched_shape;
	NODES latched_nodes;

	m_move.clear_move();


	if (m_bad_node_cost > 0 && g_rand_number_gen->random_double_number(1) < 0.05)
	{
		// we are going to try and fix a node with violations
		bad_node_generate_move();
		return;
	}


	// 1st. Choose the edge 
	
	// select an edge
	choice = g_rand_number_gen->random_number(nEdges-1);
	assert(nEdges > 0);
	assert(choice >= 0 && choice < nEdges);

	edge = m_edges[choice];
	assert(edge);


	// edge a is the edge we are going to move some way
	m_move.set_edge_a(edge);

	// find the sink node of edge a
	sink_node = edge->get_sink_node();
	assert(sink_node);


	// 2nd. Choose the move type
	//
	// we have the edge we want to move
	// randomly pick a type of move
	
	random_number = g_rand_number_gen->random_double_number(1);

	if (random_number < 0.5 || sink_node->is_DFF())
	{
	
		// make a edge swap.
		// two source nodes swap edges. 
		//debug("Making a edge swap");
		m_move.set_type(SIMPLE_EDGE_MOVE::EDGE_SWAP_BETWEEN_SOURCES);

		source_node = edge->get_source_node();
		assert(source_node);

		source_level_node = source_node->get_my_level_node();
		assert(source_level_node);

		nNodes = source_level_node->get_nNodes();

		choice = g_rand_number_gen->random_number(nNodes-1);
		assert(choice >= 0 && choice < nNodes);

		source_node2 = source_level_node->get_node(choice);
		assert(source_node2);

		output_edges = source_node2->get_output_edges();
		nOutputs = static_cast<NUM_ELEMENTS>(output_edges.size());

		while (nOutputs < 1)
		{
			choice = g_rand_number_gen->random_number(nNodes-1);
			assert(choice >= 0 && choice < nNodes);

			source_node2 = source_level_node->get_node(choice);
			assert(source_node2);

			output_edges = source_node2->get_output_edges();
			nOutputs = static_cast<NUM_ELEMENTS>(output_edges.size());
		}

		choice = g_rand_number_gen->random_number(nOutputs-1);
		assert(choice >= 0 && choice < nOutputs);

		edge = output_edges[choice];
		assert(edge);

		m_move.set_edge_b(edge);
	}
	else
	{
		// make a edge rotation
		// for the edge the source stays the same. the destination changes

		m_move.set_type(SIMPLE_EDGE_MOVE::EDGE_ROTATION);

		// pick a destination node in the sink level node to move this edge to
		sink_level_node = sink_node->get_my_level_node();
		assert(sink_level_node);

		nNodes = sink_level_node->get_nNodes();

		choice = g_rand_number_gen->random_number(nNodes-1);
		assert(choice >= 0 && choice < nNodes);

		sink_node = sink_level_node->get_node(choice);
		assert(sink_node);

		m_move.set_destination_sink_node(sink_node);
	}
}




// the cost function that gets uses a global version of wirelength
//
//
//  RETURN:
//
//  cost = Eta*(Wirelength-approx - desired Wirelength-approx) +
//         Sum over all nodes n of (Number of Violations(n))
COST_TYPE SIMPLE_EDGE_ITERATOR::get_cost()
{
	// should had a cost for having at least one unit delay edge
	assert(new_gen); 

	COST_TYPE too_many_inputs_cost = 0,
			  double_inputs_cost = 0,
			  no_inputs_cost = 0,
			  dff_loops_cost = 0,
			  bad_node_cost = 0,
			  cost = 0,
			  wirelength_cost = 0;

	CLUSTERS::const_iterator cluster_iter;
	CLUSTER * cluster = 0;
	DELAY_TYPE delay = 0,
				max_delay = m_circuit->get_delay();
	LEVEL_NODE * level_node = 0;

				 
	for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);

		level_node = cluster->Level[0];

		assert(level_node);

		dff_loops_cost = m_dff_loop_penalty_multipler*level_node->get_nDFF_loops();

		bad_node_cost += dff_loops_cost;

		for (delay = 0; delay <= max_delay; delay++)
		{
			level_node = cluster->Level[delay];
			assert(level_node);
	
			too_many_inputs_cost = level_node->get_nExcess_inputs_on_indiv_nodes();
			double_inputs_cost = level_node->get_nDouble_inputs_on_indiv_nodes();
			no_inputs_cost = level_node->get_nIndiv_nodes_too_few_inputs_cost();

			bad_node_cost += too_many_inputs_cost + double_inputs_cost + no_inputs_cost;

			wirelength_cost += level_node->get_wirelength_cost();
		}
	}

	m_bad_node_cost = bad_node_cost;
	m_wirelength_cost = wirelength_cost;

	bad_node_cost = bad_node_cost/m_nEdges;
	wirelength_cost = ABS(wirelength_cost-m_wirelength_wanted)/(m_max_index*m_nEdges);

	//debug("norm bad node cost " << bad_node_cost << "\tnorm wirelength cost: " << wirelength_cost);

	cost = m_eta*wirelength_cost+ (1-m_eta)*bad_node_cost;

	return cost;
}
// not used at present.
// the cost function where we target the wirelength-approx in each cluster
//
COST_TYPE SIMPLE_EDGE_ITERATOR::get_cost_using_wirelength_in_clusters()
{
	CLUSTERS::iterator cluster_iter;
	CLUSTER * cluster = 0;
	COST_TYPE cost = 0;

	for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);

		cost += cluster->get_final_improvement_cost(m_bad_node_cost, m_nEdges, m_eta, m_max_index);
	}

	return cost;
}
//
// PRE: add_lost_edges_back_to_nodes() has been called that added the lost edges
//      back to the node
// POST: there is no longer a difference between the weight assigned to any super edge
//       in the delay structure and the number of individual edges that have been created
//       between the two level nodes that the super edge connects to.
//
//       we prob. now have node violations
//
void SIMPLE_EDGE_ITERATOR::create_edges_to_rectify_delay_structure_with_individual_edge_assignment()
{
	CLUSTERS::iterator cluster_iter;
	CLUSTER * cluster = 0;
	DELAY_TYPE delay = 0,
				max_delay = m_circuit->get_delay();
	LEVEL_NODE * level_node = 0;
	NUM_ELEMENTS input_difference = 0;
		
	// add the lost edges back to nodes that should have had edges
	add_lost_edges_back_to_nodes();

	for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);

		// we should not have to correct the delay level 0
		for (delay = 1; delay <= max_delay; delay++)
		{
			level_node = cluster->Level[delay];
			assert(level_node);

			input_difference = level_node->get_nSimple_input_edges() - level_node->get_input_edge_weight();

			if (input_difference != 0)
			{
				assert(input_difference < 0);

				try_to_fix_level_node(level_node);
			}
		}
	}
}

// creates a list of nodes that have node violations
//
// POST: m_nodes_with_violations contains the node in the circuit that have
//       node violations
//
void SIMPLE_EDGE_ITERATOR::make_list_of_node_that_have_violations()
{
	NODES& nodes = m_circuit->get_nodes();
	NODES::const_iterator node_iter;
	NODE * node = 0;


	for (node_iter = nodes.begin(); node_iter != nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		if (node->has_double_inputs() || node->has_excess_inputs() || node->needs_inputs()
				|| node->is_dff_have_loop_back_to_itself())
		{
			if (node->is_dff_have_loop_back_to_itself())
			{
				debugif(DSIMPLE_EDGE,"Adding " << node->get_name() << " to the bad node list");
			}

			debugif(DSIMPLE_EDGE,"Adding " << node->get_name() << " to the bad node list");
			m_nodes_with_violations.push_back(node);
		}
	}
}



// PRE: level_node is valid
//      level_node has a difference between the input weight it should have 
//      vs. the input weight is currently has in terms of the delay structure
//
//      add_lost_edges_back_to_nodes() has been called that added the lost edges
//      back to the nodes
//
// POST: there is no longer a difference between the weight assigned to 
//       super edges that input into this level node and the number of individual
//       edges that connect to this level node
//
//       we prob. now have node violations for the indiv. nodes of this level node
//       
//
void SIMPLE_EDGE_ITERATOR::try_to_fix_level_node
(
	LEVEL_NODE * level_node
)
{
	assert(level_node);

	EDGES input_edges = level_node->get_input_edges();
	NUM_ELEMENTS dst_size = level_node->get_nNodes();
	EDGES::iterator edge_iter;
	EDGE * edge = 0, 
	     * indiv_edge = 0;
	NODE * source_node = 0,
		 * dst_node = 0;
	NODES nodes;
	NODES::const_iterator node_iter;
	NUM_ELEMENTS total_simple_edge_weight = 0,
				 super_edge_weight = 0,
				 weight_difference = 0,
				 nLost_nodes = 0,
				 lost_edge_weight = 0,
				 source_choice = 0,
				 dst_choice = 0;

	for (edge_iter = input_edges.begin(); edge_iter != input_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		super_edge_weight = edge->get_weight();
		total_simple_edge_weight = edge->get_total_simple_edge_weight();

		weight_difference = super_edge_weight - total_simple_edge_weight;

		assert(weight_difference >= 0);

		if (weight_difference > 0)
		{

			// get the nodes that should have had edges
			nodes = get_lost_source_nodes(edge);

			nLost_nodes = static_cast<NUM_ELEMENTS>(nodes.size());

			lost_edge_weight = 0;
			for (node_iter = nodes.begin(); node_iter != nodes.end(); node_iter++)
			{
				source_node = *node_iter;
				assert(source_node);
				lost_edge_weight += source_node->get_edges_currently_assigning();
			}
			//debug("The lost source nodes have a lost weight of " << lost_edge_weight);
			//debug("Need " << weight_difference);

			assert(lost_edge_weight >= weight_difference);

			// while we still have not assigned enough edges to satisfy what
			// the super edges wants randomly make connections even though they 
			// might create node violations
			while (weight_difference > 0)
			{
				// choose a source node
				source_choice = g_rand_number_gen->random_number(nLost_nodes-1);
				assert(source_choice >= 0 && source_choice < nLost_nodes);

				source_node = nodes[source_choice];

				// if we have not removed this node from the list
				// of source nodes make a connection from this source node 
				// to a destination node
				if (source_node != 0)
				{
					dst_choice = g_rand_number_gen->random_number(dst_size-1);
					dst_node = level_node->get_node(dst_choice);
					assert(dst_node);

					// force the connection to the destination node in the hope
					// that during iterative the node violation can be resolved
					weight_difference--;
					indiv_edge = m_circuit->create_edge(source_node, dst_node);
					indiv_edge->set_my_super_edge(edge);
					m_edges.push_back(indiv_edge);
					m_nEdges++;

					debugif(DSIMPLE_EDGE,"Connecting " << source_node->get_name() << " to " 
							<< dst_node->get_name());

					source_node->add_edges_currently_assigning(-1);

					// if the source node has no more edges to assign 
					// remove it from consideration by zeroing its pointer in
					// the list
					if (source_node->get_edges_currently_assigning() == 0)
					{
						nodes[source_choice] = 0;
					}
					assert(dst_node->is_combinational());
				}
			}
		}
	}
}


// PRE: add_lost_edges_back_to_nodes() has been called that added the lost edges
//      back to the nodes in the circuit
// RETURN: the nodes that have edges they want to connect to something
NODES SIMPLE_EDGE_ITERATOR::get_lost_source_nodes
(
	EDGE * edge
)
{
	assert(edge);

	NODE * node = 0;
	LEVEL_NODE  * source_level_node = edge->get_source_level_node(),
				* sink_level_node = edge->get_sink_level_node();
	assert(sink_level_node);
	NODES::const_iterator node_iter;
	NODES source_nodes = source_level_node->get_nodes();
	NODES lost_source_nodes;
				 

	for (node_iter = source_nodes.begin(); node_iter != source_nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		if (node->get_edges_currently_assigning() > 0)
		{
			lost_source_nodes.push_back(node);
		}

	}

	return lost_source_nodes;
}

// PRE: move is valid
// POST: nothing
// RETURN: changed_cost
COST_TYPE SIMPLE_EDGE_ITERATOR::get_changed_cost()
{
	COST_TYPE bad_node_changed_cost = m_move.get_changed_bad_node_cost(m_nEdges);
	COST_TYPE wirelength_changed_cost = m_move.get_changed_wirelength_cost(m_wirelength_cost, m_wirelength_wanted,
																		m_max_index, m_nEdges);
	COST_TYPE changed_cost = 0.0;

	changed_cost = m_eta*wirelength_changed_cost + (1 - m_eta)*bad_node_changed_cost;

	//debug("Changed cost " << changed_cost);

	return changed_cost;
}

SIMPLE_EDGE_MOVE::SIMPLE_EDGE_MOVE()
{
	m_type 		= NONE;
	m_edge_a 	= 0;
	m_edge_b 	= 0;
	m_dst_sink_node = 0;

	assert(g_options);
	m_dff_loop_penalty_multipler = g_options->get_final_edge_assignment_dff_loop_penalty_multiplier();
}


SIMPLE_EDGE_MOVE::~SIMPLE_EDGE_MOVE()
{
}

// 
// valid move if:
//  a) we have an edge type
//  b) the sink nodes are different
//  c) the edges are different
//	d) we have a move that will create a dff buffer
//	   i.e.

//			  dff (sink node a)
//			  |
//			  buffer (source node b)
//			  |
//			  edge wanting to connect to sink node a
//
//	RETURN: validity of move
bool SIMPLE_EDGE_MOVE::is_valid() const
{

	bool valid = (m_type != NONE && 
				  m_edge_a != m_edge_b && 
				  m_edge_a->get_sink_node() != m_dst_sink_node);

	if (m_type == SIMPLE_EDGE_MOVE::EDGE_SWAP_BETWEEN_SOURCES)
	{
		// we want to swap edge_a and edge_b
		// take a look at their sources and sink to see if it would be ok to swap them

		NODE * sink_node_a = m_edge_a->get_sink_node();
		NODE * sink_node_b = m_edge_b->get_sink_node();
		assert(sink_node_a && sink_node_b);
		
		valid = valid && ((sink_node_a->is_DFF() && sink_node_b->is_DFF()) || 
							(! sink_node_a->is_DFF() && ! sink_node_b->is_DFF()));

		NODE * source_node_a = m_edge_a->get_source_node();
		NODE * source_node_b = m_edge_b->get_source_node();

		// check to make sure we are not connecting dff to themselves
		if (sink_node_a->is_DFF() || sink_node_b->is_DFF())
		{
			assert(source_node_a && source_node_b);

			valid = valid && (source_node_a != sink_node_a) &&
							 (source_node_b != sink_node_a);

			valid = valid && ! ( sink_node_a->is_DFF() && source_node_b->is_buffer() && 
								source_node_b->are_nodes_connected(sink_node_a));

			valid = valid && ! ( sink_node_b->is_DFF() && source_node_a->is_buffer() && 
								source_node_a->are_nodes_connected(sink_node_b));
		}
	}

	
	return valid;
}
// POST: the move has been cleared
void SIMPLE_EDGE_MOVE::clear_move()
{
	m_type 		= NONE;
	m_edge_a 	= 0;
	m_edge_b 	= 0;
	m_dst_sink_node = 0;
}
void SIMPLE_EDGE_MOVE::make_temporary_move()
{
	make_move();
}
void SIMPLE_EDGE_MOVE::unmake_temporary_move()
{
	make_move();
}


// POST: the move has been made.
//       the super edges involved have been updated
//       the node violation list has been updated
//       the length of the edges have been updated
void SIMPLE_EDGE_MOVE::make_and_commit_move
(
	NODES& bad_node_list,
	COST_TYPE& wirelength_cost,
	COST_TYPE& bad_node_cost
)
{
	update_costs(wirelength_cost, bad_node_cost);

	if (m_type == SIMPLE_EDGE_MOVE::CHANGE_DFF_EDGE_ASSIGNMENT) 
	{
		change_super_edges_invoved();
	}

	make_move();
	update_bad_node_list(bad_node_list);
	m_dst_sink_node = 0;
	m_edge_a = 0;
	m_edge_b = 0;

	m_type = NONE;
}
	
void SIMPLE_EDGE_MOVE::change_super_edges_invoved()
{
	EDGE * super_edge_a = m_edge_a->get_my_super_edge();
	EDGE * super_edge_b = m_edge_b->get_my_super_edge();

	if (super_edge_a != super_edge_b)
	{
		// add the weight of the new edges
		// if the edge doesn't exist we might have to create it
		LEVEL_NODE * source_node_a = super_edge_a->get_source_level_node();
		LEVEL_NODE * source_node_b = super_edge_b->get_source_level_node();
		LEVEL_NODE * dff_node_a = super_edge_a->get_sink_level_node();
		LEVEL_NODE * dff_node_b = super_edge_b->get_sink_level_node();

		EDGE * new_super_edge_1 = dff_node_a->find_edge_with_source(source_node_b);
		EDGE * new_super_edge_2 = dff_node_b->find_edge_with_source(source_node_a);
	
		// remove the weight on the old super edge
		super_edge_a->add_weight_to_dff_super_edge(-1);
		super_edge_b->add_weight_to_dff_super_edge(-1);

		if (! new_super_edge_1)
		{
			new_super_edge_1 = m_circuit->create_edge(source_node_b, dff_node_a);
		}
		if (! new_super_edge_2)
		{
			new_super_edge_2 = m_circuit->create_edge(source_node_a, dff_node_b);
		}
		new_super_edge_1->add_weight_to_dff_super_edge(1);
		new_super_edge_2->add_weight_to_dff_super_edge(1);

	
	}
}



//
// POST: we have made the move
//
//note:
// sigh. this will screw up the data structures in cluster
//       in regards to what edges each cluster contains.
//       after our algorithm is over we should update the edges 
//       assigned to the clusters
//       i don't want to update them here because that would slow things down too much
void SIMPLE_EDGE_MOVE::make_move()
{
	NODE * sink_node_1 = 0,
		 * sink_node_2 = 0;

	if (m_type == SIMPLE_EDGE_MOVE::EDGE_ROTATION)
	{
		sink_node_1 = m_edge_a->get_sink_node();
		
		assert(sink_node_1 && m_dst_sink_node);

		sink_node_1->remove_input_edge(m_edge_a);

		m_edge_a->set_sink_node(m_dst_sink_node);

		m_dst_sink_node->add_input_edge(m_edge_a);

		// below make it possible to reverse this 
		// move by calling make temporary move again
		m_dst_sink_node = sink_node_1;
	}
	else if (m_type == SIMPLE_EDGE_MOVE::CHANGE_DFF_EDGE_ASSIGNMENT)
	{
		assert(m_edge_a && m_edge_b);

		EDGE * super_edge_a = m_edge_a->get_my_super_edge();
		EDGE * super_edge_b = m_edge_b->get_my_super_edge();

		get_sink_nodes_involved(sink_node_1, sink_node_2);
		assert(sink_node_1 && sink_node_2);

		sink_node_1->remove_input_edge(m_edge_a);
		sink_node_2->remove_input_edge(m_edge_b);

		m_edge_a->set_sink_node(sink_node_2);
		m_edge_b->set_sink_node(sink_node_1);

		sink_node_2->add_input_edge(m_edge_a);
		sink_node_1->add_input_edge(m_edge_b);

		m_edge_a->set_my_super_edge(super_edge_b);
		m_edge_b->set_my_super_edge(super_edge_a);
	}
	else
	{
		assert(m_type == SIMPLE_EDGE_MOVE::EDGE_SWAP_BETWEEN_SOURCES);
		assert(m_edge_a && m_edge_b);

		LENGTH_TYPE length_a = m_edge_a->get_length(),
					length_b = m_edge_b->get_length();
		EDGE * super_edge_a = m_edge_a->get_my_super_edge();
		EDGE * super_edge_b = m_edge_b->get_my_super_edge();

		get_sink_nodes_involved(sink_node_1, sink_node_2);
		assert(sink_node_1 && sink_node_2);

		sink_node_1->remove_input_edge(m_edge_a);
		sink_node_2->remove_input_edge(m_edge_b);

		m_edge_a->set_sink_node(sink_node_2);
		m_edge_b->set_sink_node(sink_node_1);

		sink_node_2->add_input_edge(m_edge_a);
		sink_node_1->add_input_edge(m_edge_b);

		assert((sink_node_1->is_DFF() && sink_node_2->is_DFF()) || 
				(! sink_node_2->is_DFF() && ! sink_node_2->is_DFF()));

		m_edge_a->set_length(length_b);
		m_edge_b->set_length(length_a);

		m_edge_a->set_my_super_edge(super_edge_b);
		m_edge_b->set_my_super_edge(super_edge_a);
	}
}

void SIMPLE_EDGE_MOVE::update_bad_node_list
(
	NODES& bad_node_list
)
{
	NODE * sink_node_1 = 0,
		 * sink_node_2 = 0;
	NODE * source_node_1 = 0,
		 * source_node_2 = 0;

	get_sink_nodes_involved(sink_node_1, sink_node_2);
	assert(sink_node_1 && sink_node_2);

	add_remove_node_from_node_violation_list_if_necessary(sink_node_1, bad_node_list);
	add_remove_node_from_node_violation_list_if_necessary(sink_node_2, bad_node_list);

	// make sure if the source nodes are dff that we haven't create dff loops
	get_source_nodes_involved(source_node_1, source_node_2);
	assert(source_node_1 && source_node_2);

	add_remove_node_from_node_violation_list_if_necessary(source_node_1, bad_node_list);
	add_remove_node_from_node_violation_list_if_necessary(source_node_2, bad_node_list);
}


void SIMPLE_EDGE_MOVE::update_costs
(
	COST_TYPE& wirelength_cost,
	COST_TYPE& bad_node_cost
)
{
	wirelength_cost += get_wirelength_change_of_sink_nodes_involved();
	bad_node_cost += get_bad_node_cost_change_of_sink_nodes_involved();
}


// PRE: node is valid
// POST: if the sink node has a node violation it is on the node violation list
//       else it is not on the node violation list
//
void SIMPLE_EDGE_MOVE::add_remove_node_from_node_violation_list_if_necessary
(
	NODE * sink_node,
	NODES& nodes_with_violations
)
{

	NODES::iterator node_iter;

	if (sink_node->has_double_inputs() || sink_node->has_excess_inputs() || sink_node->needs_inputs() 
			|| sink_node->is_dff_have_loop_back_to_itself())
	{
		if (find(nodes_with_violations.begin(), nodes_with_violations.end(), sink_node) 
				== nodes_with_violations.end())
		{
			nodes_with_violations.push_back(sink_node);
		}
		else
		{
			// we have seen it before. don't add it again
		}
	}
	else
	{
		// this node either a) was not in violation for b) is no longer in violation and needs to be 
		// removed from the node violation list
		// if it is there remove it

		node_iter = find(nodes_with_violations.begin(), nodes_with_violations.end(), sink_node);

		if (node_iter != nodes_with_violations.end())
		{
			nodes_with_violations.erase(node_iter);
		}
	}
}
COST_TYPE SIMPLE_EDGE_MOVE::get_bad_node_cost_of_sink_nodes_involved()
{
	NODE * sink_node_1 = 0,
		 * sink_node_2 = 0;
	COST_TYPE bad_node_cost = 0.0;

	get_sink_nodes_involved(sink_node_1, sink_node_2);

	bad_node_cost = sink_node_1->get_double_input_cost() + sink_node_2->get_double_input_cost() +
		sink_node_1->get_nExcess_inputs() + sink_node_2->get_nExcess_inputs() +
		sink_node_1->get_nInputs_still_needed() + sink_node_2->get_nInputs_still_needed();

	bad_node_cost += m_dff_loop_penalty_multipler*m_circuit->get_nDFF_loops();

	return bad_node_cost;
}


void SIMPLE_EDGE_MOVE::get_sink_nodes_involved
(
	NODE *& sink_node_1,
	NODE *& sink_node_2
)
{
	if (m_type == SIMPLE_EDGE_MOVE::EDGE_ROTATION)
	{
		sink_node_1 = m_dst_sink_node;
		sink_node_2 = m_edge_a->get_sink_node();
	}
	else
	{
		sink_node_1 = m_edge_a->get_sink_node();
		sink_node_2 = m_edge_b->get_sink_node();
	}
}
void SIMPLE_EDGE_MOVE::get_source_nodes_involved
(
	NODE *& source_node_1,
	NODE *& source_node_2
)
{
	source_node_1 = m_edge_a->get_source_node();
	if (m_type == SIMPLE_EDGE_MOVE::EDGE_ROTATION)
	{
		source_node_2 = m_edge_a->get_source_node();
	}
	else
	{
		source_node_2 = m_edge_b->get_source_node();
	}
}
void SIMPLE_EDGE_MOVE::get_sink_clusters
(
	CLUSTER *& sink_1_cluster, 
	CLUSTER *& sink_2_cluster
)
{
	NODE * sink_node_1 = 0,
		 * sink_node_2 = 0;

	get_sink_nodes_involved(sink_node_1, sink_node_2);
	assert(sink_node_1 && sink_node_2);

	LEVEL_NODE * level_node_1 = sink_node_1->get_my_level_node(),
			   * level_node_2 = sink_node_2->get_my_level_node();
	
	assert(level_node_1 && level_node_2);

	sink_1_cluster = level_node_1->get_my_cluster();
	sink_2_cluster = level_node_2->get_my_cluster();

	assert(sink_1_cluster && sink_2_cluster);
}

COST_TYPE SIMPLE_EDGE_MOVE::get_wirelength_cost_of_sink_nodes_involved()
{
	NODE * sink_node_1 = 0,
		 * sink_node_2 = 0;
	COST_TYPE wirelength_cost = 0.0;

	get_sink_nodes_involved(sink_node_1, sink_node_2);

	wirelength_cost = sink_node_1->get_wirelength_cost() + sink_node_2->get_wirelength_cost();

	return wirelength_cost;
}

COST_TYPE SIMPLE_EDGE_MOVE::get_wirelength_change_of_sink_nodes_involved()
{
	make_temporary_move();

	COST_TYPE changed_cost = get_wirelength_cost_of_sink_nodes_involved();

	unmake_temporary_move();

	changed_cost -= get_wirelength_cost_of_sink_nodes_involved();

	return changed_cost;
}


COST_TYPE SIMPLE_EDGE_MOVE::get_bad_node_cost_change_of_sink_nodes_involved()
{
	make_temporary_move();

	COST_TYPE changed_cost = get_bad_node_cost_of_sink_nodes_involved();

	unmake_temporary_move();

	changed_cost -= get_bad_node_cost_of_sink_nodes_involved();

	return changed_cost;
}




COST_TYPE SIMPLE_EDGE_MOVE::get_changed_bad_node_cost(const COST_TYPE& nEdges)
{
	COST_TYPE changed_cost = get_bad_node_cost_change_of_sink_nodes_involved();

	//debug("Bad node changed cost " << changed_cost);

	COST_TYPE bad_node_cost = changed_cost/nEdges;

	//debug("Bad node norm changed cost " << bad_node_cost);
	
	return bad_node_cost;
}

// get the changed cost using the global wirelength model.
COST_TYPE SIMPLE_EDGE_MOVE::get_changed_wirelength_cost
(
	const COST_TYPE& wirelength_cost, 
	const COST_TYPE& wirelength_wanted,
	const COST_TYPE& max_index,
	const COST_TYPE& nEdges
)
{

	COST_TYPE changed_wirelength = 0.0,
			  changed_wirelength_of_sink_nodes = 0.0,
			  new_wirelength = 0.0,
			  changed_cost = 0.0;


	changed_wirelength_of_sink_nodes = get_wirelength_change_of_sink_nodes_involved();
		
	new_wirelength = 	(wirelength_cost + changed_wirelength_of_sink_nodes);

	changed_wirelength = ABS(new_wirelength - wirelength_wanted) - ABS(wirelength_cost - wirelength_wanted);
		
	changed_cost = changed_wirelength/(max_index*nEdges);

	//debug("wirelength Changed cost " << changed_cost);
	
	return changed_cost;
}
// get the changed cost using the clustered wirelength model
COST_TYPE SIMPLE_EDGE_MOVE::get_changed_wirelength_cost2
(
	const COST_TYPE& wirelength_cost, 
	const COST_TYPE& wirelength_wanted,
	const COST_TYPE& max_index,
	const COST_TYPE& nEdges
)
{
	NODE * sink_node_1 = 0,
		 * sink_node_2 = 0;
	COST_TYPE changed_wirelength = 0;
	COST_TYPE wirelength_change_sink_1 = 0,
			  wirelength_change_sink_2 = 0,
			  wirelength_change = 0;

	get_sink_nodes_involved(sink_node_1, sink_node_2);

	CLUSTER * sink_1_cluster = 0,
			   * sink_2_cluster = 0;

	get_sink_clusters(sink_1_cluster, sink_2_cluster);

	if (sink_1_cluster != sink_2_cluster)
	{
		// the sink nodes are in different comb graphs
		// get the change of wirelength in each comb graph
		// and then find the how that changes each comb graphs wirelength cost

		get_wirelength_change_of_sink_nodes_involved(wirelength_change_sink_1, wirelength_change_sink_2);
	
		changed_wirelength = sink_1_cluster->get_changed_wirelength_cost(wirelength_change_sink_1);
		changed_wirelength += sink_2_cluster->get_changed_wirelength_cost(wirelength_change_sink_2);
	}
	else
	{
		// both sink nodes are in the same comb graph
		// get the change of wirelength of the comb graph
		// and then find the how that changes the comb graphs wirelength cost


		wirelength_change = get_wirelength_change_of_sink_nodes_involved();
		changed_wirelength = sink_1_cluster->get_changed_wirelength_cost(wirelength_change);
	}	

	changed_wirelength = changed_wirelength/(max_index*nEdges);

	//debug("wirelength Changed cost " << changed_cost);
	
	return changed_wirelength;
}


void SIMPLE_EDGE_MOVE::get_wirelength_change_of_sink_nodes_involved
(
	COST_TYPE& wirelength_change_sink_1,
	COST_TYPE& wirelength_changed_sink_2
)
{
	NODE * sink_node_1 = 0,
		 * sink_node_2 = 0;

	get_sink_nodes_involved(sink_node_1, sink_node_2);
	assert(sink_node_1 && sink_node_2);

	make_temporary_move();

	wirelength_change_sink_1 = sink_node_1->get_wirelength_cost();
	wirelength_changed_sink_2 = sink_node_2->get_wirelength_cost();

	unmake_temporary_move();

	wirelength_change_sink_1 -= sink_node_1->get_wirelength_cost();
	wirelength_changed_sink_2 -= sink_node_2->get_wirelength_cost();
}


// PRE: success_rate is the number of moves accepted
// POST: m_temperature has been updated
void SIMPLE_EDGE_ITERATOR::update_temperature
(
	const double& success_rate
)
{
	double m_temperature = 0.0;
	bool nodes_with_bad_things = false;

	if (success_rate > 0.96) 
	{
	   m_temperature = m_temperature * 0.5; 
	}
	else if (success_rate > 0.8) 
	{
	   m_temperature = m_temperature * 0.9;
	}
	else if (success_rate > 0.15 || nodes_with_bad_things)
	{
	   m_temperature = m_temperature * 0.95;
	}
	else 
	{
	   m_temperature = m_temperature * 0.8; 
	}
}

// prints the horizontal positions of nodes in the circuit
//
void SIMPLE_EDGE_ITERATOR::print_node_positions()
{
	NODES & nodes = m_circuit->get_nodes();
	NODES::const_iterator node_iter;
	NODE * node = 0;

	for (node_iter = nodes.begin(); node_iter != nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		debug("Node " << node->get_name() << " position " << node->get_index());	
	}
}

// POST: the wirelength-approx used by each cluster has been printed as well as the 
//       the wirelength-approx desired for each cluster
void SIMPLE_EDGE_ITERATOR::print_wirelength_stats() const
{
	CLUSTERS::const_iterator cluster_iter;
	CLUSTER * cluster = 0;
	CLUSTER_SPEC * comb_spec = 0;

	debug("wirelength wanted\t\twirelength found\n");

	for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);

		comb_spec = cluster->get_comb_spec();
		assert(comb_spec);

		debug(comb_spec->get_wirelength() << "\t\t" << cluster->get_wirelength());
	}
}


// generate a move to solve dff loops
//
// PRE: sink_node is a valid node that 
//      is a flip-flop that loops back to itself
// POST: m_move has been set
//       the move may or may not be valid
//
void SIMPLE_EDGE_ITERATOR::generate_move_to_solve_dff_loop
(
	NODE * sink_node
)
{
	LEVEL_NODE * sink_level_node = 0;
	EDGES input_edges;
	EDGE * edge_that_causes_loop = 0;

	assert(sink_node->is_dff_have_loop_back_to_itself());

	sink_level_node = sink_node->get_my_level_node();
	assert(sink_level_node);


	edge_that_causes_loop = sink_node->get_dff_input_edge();
	assert(edge_that_causes_loop);

	// see if we can swap an edge with another latched node/dff combo
	generate_move_that_swaps_edge_between_dff(edge_that_causes_loop);

	assert(new_gen);
	// What can we do about this dff loop?
	// Here are our choices as I see them:
	//
	// 1. if there is more than one latched node is the source cluster of
	//    this edge maybe we can swap this edge for its edge
	//    this won't change spec.
	//
	// 2. if the DFF was latched during node splitting because it has zero 
	//    fanout maybe we swap fanout with a PI node and make the PI node a latched node
	//    this won't change spec
	//
	// 3. we can change the DFF superstructure
	//    this will change spec
	//
	//
	// Right now we are oding only 1. Maybe investigate #2,3 ?
}


// Latch measure the connectivity between clusters.
// it does not measure the connectivity between level nodes in the various clusters
// 
// therefore we can change the latched node that connects to our flip-flop
// without changing the specification
// 
// try to find another latched node to switch edges with
//
// PRE: edge_that_causes_loop is an edge that connects to the input
//      of a flip-flop that causes a dff loop
// POST: m_move has been set 
//       The two edges in the move might be the same if 
//       we did not pick a different latched node
//
void SIMPLE_EDGE_ITERATOR::generate_move_that_swaps_edge_between_dff
(
	EDGE * edge_that_causes_loop
)
{
	assert(edge_that_causes_loop);


	NODE * latched_node_a = edge_that_causes_loop->get_source_node();
	assert(latched_node_a);

	LEVEL_NODE * level_node_a = latched_node_a->get_my_level_node();
	assert(level_node_a);
	CLUSTER * cluster = level_node_a->get_my_cluster();
	assert(cluster);

	NODE * another_latched_node = get_a_latched_node_in_same_cluster(cluster);

	EDGE * edge_b = another_latched_node->get_output_edge_that_connects_to_DFF();

	m_move.set_edge_a(edge_that_causes_loop);
	m_move.set_edge_b(edge_b);	
	m_move.set_type(SIMPLE_EDGE_MOVE::CHANGE_DFF_EDGE_ASSIGNMENT);
}

// PRE: cluster is valid
// RETURNS: a latched node in the cluster
//
NODE* SIMPLE_EDGE_ITERATOR::get_a_latched_node_in_same_cluster
(
	const CLUSTER * cluster
) const
{
	SHAPE latched_shape = cluster->get_latched_shape();
	NUM_ELEMENTS random_level = g_rand_number_gen->discrete_pmf(latched_shape);
	assert(random_level >= 0 && random_level <= cluster->get_delay());
	
	LEVEL_NODE * level_node_b = cluster->Level[random_level];

	NUM_ELEMENTS nLatched_nodes = level_node_b->get_nLatched();
	NUM_ELEMENTS random_node_num = g_rand_number_gen->random_number(nLatched_nodes-1);
	assert(random_node_num >= 0 && random_node_num < nLatched_nodes);

	NODES latched_nodes = level_node_b->get_latched_nodes();

	NODE * latched_node = latched_nodes[random_node_num];
	assert(latched_node);

	return latched_node;
}

// PRE: our clusters prob. don't have the simple edges they think they do.
// POST: the clusters have the edges they think they do
void SIMPLE_EDGE_ITERATOR::update_cluster_to_reflect_which_simple_edges_they_contain()
{
	EDGES::iterator edge_iter;
	EDGE * edge = 0;

	CLUSTERS::iterator cluster_iter;
	CLUSTER_NUMBER_TYPE source_cluster_number = 0, sink_cluster_number = 0;
	CLUSTER * cluster = 0;
	m_edges = m_circuit->get_simple_edges();


	// clear the simple edges the clusters have
	for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		cluster->clear_simple_edges();
	}


	// add them back
	for (edge_iter = m_edges.begin(); edge_iter != m_edges.end(); edge_iter++)
	{
		edge = *edge_iter;

		source_cluster_number = edge->get_source_cluster_number();
		sink_cluster_number = edge->get_sink_cluster_number();

		assert(source_cluster_number >= 0 && source_cluster_number < m_circuit->get_nClusters());
		assert(sink_cluster_number >= 0 && sink_cluster_number < m_circuit->get_nClusters());

		if (source_cluster_number == sink_cluster_number)
		{
			m_clusters[source_cluster_number]->add_intra_cluster_simple_edge(edge);
		}
		else
		{
			m_clusters[source_cluster_number]->add_inter_cluster_output_simple_edge(edge);
			m_clusters[sink_cluster_number]->add_inter_cluster_input_simple_edge(edge);
		}
	}
}

