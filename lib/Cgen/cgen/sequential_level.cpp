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




#include "gen.h"
#include "sequential_level.h"
#include "matrix3D.h"
#include <numeric>
#include <math.h>
#include <algorithm>

#ifndef VISUAL_C
#include <values.h>
#else
const long MAXLONG = 2147483647; // kludge magic number 2^31-1
#endif

SEQUENTIAL_LEVEL::SEQUENTIAL_LEVEL
(
	CIRCUIT_GRAPH * circuit
)
{
	assert(circuit);
	m_circuit = circuit;

	m_clusters 	= m_circuit->get_clusters();
	m_nClusters = m_circuit->get_nClusters();
    m_Comb.resize(m_nClusters, m_nClusters);

	CIRCUIT_SPEC * circuit_spec = m_circuit->get_circuit_spec();
	assert(circuit_spec);

	m_max_comb_delay = circuit_spec->get_max_comb_delay();
	m_Comb_spec = circuit_spec->get_Comb();
	m_Comb_spec_sum = 0.0;

	// create the delay levels which are a collection of the level nodes
	// in various clusters with the minimum possible delay
	DELAY_LEVEL* delay_level = 0;
	DELAY_TYPE delay;

	for (delay = 0; delay <= m_max_comb_delay; delay++)
	{
		delay_level = new DELAY_LEVEL;
		m_delay_levels.push_back(delay_level);
	}

	m_node_id = 0;

	m_display_too_many_inputs_msg = false;
	m_sanity_check_fatal		  = false;


	m_alpha = g_options->get_delay_structure_alpha();
	m_beta =  g_options->get_delay_structure_beta();
	m_gamma = g_options->get_delay_structure_gamma();
}
SEQUENTIAL_LEVEL::~SEQUENTIAL_LEVEL()
{
	EDGES::iterator edges_iterator;
	EDGE * edge = 0;
	DELAY_TYPE delay = 0;

	for (edges_iterator = m_edges.begin(); edges_iterator != m_edges.end(); edges_iterator++)
	{
		edge = *edges_iterator;
		assert(edge);

		delete edge;
	}
	for (delay = 0; delay <= m_max_comb_delay; delay++)
	{
		delete m_delay_levels[delay];
	}
	m_delay_levels.clear();

	m_edges.clear();
}

SEQUENTIAL_LEVEL::SEQUENTIAL_LEVEL(const SEQUENTIAL_LEVEL  & another_sequential_level_generator)
{

}
SEQUENTIAL_LEVEL& SEQUENTIAL_LEVEL::operator=(const SEQUENTIAL_LEVEL  & another_sequential_level_generator)
{
	assert(false);
	return (*this);
}

// Creates the delay structure - the large scale connectivity in the circuit 
//
//
// from Section 4.1 of Thesis
//
// The input to this step is the delay structure
// characterization which consists of the Comb and Latched matrices, dmax, and for each cluster the node
// shape, the input shape, the output shape, the latched shape, the intra-cluster edge length distribution,
// the inter-cluster input edge length distribution, and the inter-cluster output edge length distribution.
// 
// The delay structure of a circuit can be thought of as the large-scale connectivity between delay levels in
// the clusters of the circuit. Individual nodes in the delay levels of a cluster are aggregated into level
// nodes. Edges between the level nodes are aggregated into super edges where a super edge is an edge between
// two levels nodes with a weight equal to the number of edges intended to be formed between the two level
// nodes. With the aggregations, we define the delay structure of a circuit as the graph 
// delay structure =(Level Nodes, Super Edges), where Level Nodes represents the set of all level nodes 
// and Super Edges represents the set of all possible super edges in a circuit. 
//
// The desired output is the delay structure with the weights of the super edges assigned to ensure that 
// each node will be able to have its delay level correctly set, the delay structure will not force 
// node violations to be made during final edge assignment, and deviations from the given characterization 
// are minimized. 
// 
// 
// PRE: g_rand_number_gen is valid 
//      g_options is valid 
//
// POST: we have created the large scale connectivity in the circuit.
//
//
void SEQUENTIAL_LEVEL::create_delay_structure() 
{
	create_super_edges_between_level_nodes();
	find_solution();
}

/*

    This function creates the super edges between level nodes and assigns the maximum
	weight to each edge based on the specification.
   
	The algorithm to create the initial solution proceeds in four steps. 
	In Step 1, we create the level nodes and super edges between the level nodes. 
	(the level nodes were created previously in the cluster's construction)
	In Step 2, we assign the maximum possible weight to each super edge based on the input specifications. 
	If the super edge is intra-cluster, we limit the maximum weight by the cluster's 
	intra-cluster edge length distribution, input shape, output shape, and by the possible 
	number of unique connections between nodes in the source level node and nodes in the sink 
	level node. If the super edge is inter-cluster, we limit the maximum weight by the source cluster's
	inter-cluster output edge length distribution and output shape, the sink cluster's inter-cluster 
	input edge length distribution and input shape, and by the possible number of unique connections 
	between nodes in the source level node and nodes in the sink level node. 
	The maximum weight assigned to a super edge is used as a starting point in the creation of the 
	initial solution and as an upper bound on super edge weight when we later try to improve the 
	initial solution in Steps 3 and Step 4.

   
   e.g. a picture a super edges and level nodes

   o = level nodes
   

   o    o   o 
   | \/ |\ /| \
   |/ \ | /\|  |
   o    o   o  |
   | \  |  /|  |
   |  \ | / | /
   o    o   o


PRE:  level nodes have been created
POST: we have super edges between level nodes.
      each super edge has been assigned the maximum weight permitted by the specification

*/

void SEQUENTIAL_LEVEL::create_super_edges_between_level_nodes()
{
	DELAY_TYPE source_level = 0;
	DELAY_TYPE sink_level = 0;


	debugif(DSEQ_LEVEL, "Creating super edges between level nodes");

	assert(static_cast<unsigned>(m_max_comb_delay +1) == m_delay_levels.size());
	for (source_level = 0; source_level < m_max_comb_delay; source_level++)
	for (sink_level = source_level+1; sink_level <= m_max_comb_delay; sink_level++)
	{

		//debug("Source delay " << source_level);
		//debug("Sink delay " << sink_level);
		create_super_edges_between_levels(m_delay_levels[source_level], m_delay_levels[sink_level]);
	}

	sanity_check_all();
}


// create all combinational edges between the source and sink delay levels.
//
// PRE:  source_delay_level and sink_delay_level are valid 
//       sink delay level > source delay level
// POST: we have super edges between the level nodes with each super edge having the
//       maximum weight possible under the circuit specification
//
void SEQUENTIAL_LEVEL::create_super_edges_between_levels
(
	DELAY_LEVEL * source_delay_level,
	DELAY_LEVEL * sink_delay_level
)
{
	assert(source_delay_level && sink_delay_level);
	DELAY_LEVEL::iterator source_level_node_iter;
	DELAY_LEVEL::iterator sink_level_node_iter;
	LEVEL_NODE * source_level_node = 0;
	LEVEL_NODE * sink_level_node = 0;

	assert(source_delay_level->size() == m_clusters.size());
	assert(sink_delay_level->size() == m_clusters.size());

	for (source_level_node_iter = source_delay_level->begin(); source_level_node_iter != 
									source_delay_level->end(); source_level_node_iter++)
	{
		source_level_node = *source_level_node_iter;
		assert(source_level_node);

		for (sink_level_node_iter = sink_delay_level->begin(); sink_level_node_iter != 
										sink_delay_level->end(); sink_level_node_iter++)
		{
			sink_level_node = *sink_level_node_iter;
			assert(sink_level_node);

			create_combinational_edge_with_max_weight(source_level_node, sink_level_node);	
		}
	}

	// just check out results to make sure that each level node 
	// has enough delay defining edges

	// if we were assigning length 1 edges make sure that we assigned enough 
	// length 1 edges to ensure that the nodes have their delay level defined
	// (each individual node needs at least one length 1 edge)
	if (sink_level_node->get_delay_level() - source_level_node->get_delay_level() == 1)
	{
		for (sink_level_node_iter = sink_delay_level->begin(); sink_level_node_iter != 
										sink_delay_level->end(); sink_level_node_iter++)
		{
			sink_level_node = *sink_level_node_iter;
			assert(sink_level_node->get_delay_difference() >= 0);
		}
	}
}


//
//	creates the edge and set the edge weight as the maximum allowable
//	based on the specification and constraints
// 
//  PRE: source_level_node and sink_level_node are valid
//  POST: a super edge has been created between the level node
//
void SEQUENTIAL_LEVEL::create_combinational_edge_with_max_weight
(
	LEVEL_NODE * source_level_node,
	LEVEL_NODE * sink_level_node
)
{
	assert(source_level_node && sink_level_node);

	CLUSTER_NUMBER_TYPE sink_cluster_number = sink_level_node->get_cluster_number();
	CLUSTER_NUMBER_TYPE source_cluster_number = source_level_node->get_cluster_number();
	assert(source_cluster_number >= 0 && source_cluster_number < static_cast<signed>(m_clusters.size()));
	assert(sink_cluster_number >= 0 && sink_cluster_number < static_cast<signed>(m_clusters.size()));

	NUM_ELEMENTS sink_edge_weight, source_edge_weight, edge_weight;
	NUM_ELEMENTS inter_cluster_connections = 0;
	LENGTH_TYPE edge_length = 0;
	CLUSTER_SPEC * source_comb_spec = 0;
	CLUSTER_SPEC * sink_comb_spec = 0;

	CIRCUIT_SPEC * circ_spec = m_circuit->get_circuit_spec();
	assert(circ_spec);

	CLUSTER * source_cluster = m_clusters[source_cluster_number];
	assert(source_cluster->get_cluster_number() == source_cluster_number);

	EDGE * edge = create_combinational_edge(source_level_node, sink_level_node);
	assert(edge);

	edge_length = edge->get_length();

	if (edge->is_intra_cluster())
	{
		sink_comb_spec = m_clusters[sink_cluster_number]->get_comb_spec();
		edge_weight = sink_comb_spec->get_intra_cluster_edge_lengths()[edge_length];
	}
	else
	{
		assert(edge->is_inter_cluster());

		source_comb_spec = m_clusters[source_cluster_number]->get_comb_spec();
		sink_comb_spec = m_clusters[sink_cluster_number]->get_comb_spec();

		sink_edge_weight = sink_comb_spec->get_inter_cluster_input_edge_lengths()[edge_length];
		source_edge_weight = source_comb_spec->get_inter_cluster_output_edge_lengths()[edge_length];

		edge_weight = MIN(source_edge_weight, sink_edge_weight);

		// edge weight is also the minimum of the number of connections between the 2 clusters
		inter_cluster_connections = (*m_Comb_spec)(source_cluster_number, sink_cluster_number);
		edge_weight = MIN(edge_weight, inter_cluster_connections);
	}

	
	// unnecessary constraint based on sink_level_node->get_assigned_in() constraint
	// edge_weight = MIN(edge_weight, sink_level_node->get_nNodes() * kin);
	//

	edge_weight = MIN(edge_weight, source_level_node->get_assigned_out());
	edge_weight = MIN(edge_weight, sink_level_node->get_assigned_in());

	// make sure we get no double connections between node nodes
	edge_weight = MIN(edge_weight, source_level_node->get_nNodes() * sink_level_node->get_nNodes());

	assert(edge_weight >= 0);

	//debug("Setting weight to sink node " << sink_level_node->get_name() << " edge " << edge->get_info() << " to " << edge_weight);
	edge->set_max_weight(edge_weight);
	edge->set_weight_and_update_level_nodes(edge_weight);
}
	

// Creates a super edge between level ndoes
//
// PRE: source and sink level nodes are valid
// POST: a super edge has been created 
//       the edge has been added to m_edges
//       the edge has 0 weight
//
EDGE* SEQUENTIAL_LEVEL::create_combinational_edge
(
	LEVEL_NODE * source_level_node,
	LEVEL_NODE * sink_level_node
)
{
	assert(source_level_node && sink_level_node);

	//output_edge_debug_message(source_level_node, sink_level_node);

	EDGE * edge = m_circuit->create_edge(source_level_node, sink_level_node);
	assert(edge);
	edge->set_weight_and_update_level_nodes(0);

	// this id is important for quick edge_weight_state updates
	assert(static_cast<signed>(m_edges.size()) == edge->get_id());

	m_edges.push_back(edge);

	if (edge->is_inter_cluster())
	{
		m_inter_cluster_edges.push_back(edge);
	}

	return edge;
}

// 
// PRE: source_level_node and sink_level_node are valid
// POST: we have created a dff super edge been the two level nodes
//       the weight of the edge is zero.
//
EDGE* SEQUENTIAL_LEVEL::create_dff_edge
(
	LEVEL_NODE * source_level_node,
	LEVEL_NODE * sink_level_node
)
{
	assert(source_level_node && sink_level_node);
	assert(sink_level_node->get_delay_level() == 0);

	//output_edge_debug_message(source_level_node, sink_level_node);

	EDGE * edge = m_circuit->create_edge(source_level_node, sink_level_node);
	assert(edge);
	edge->set_weight(0);

	// this id is important for quick edge_weight_state updates
	assert(static_cast<signed>(m_edges.size()) == edge->get_id());

	m_edges.push_back(edge);

	if (edge->is_inter_cluster())
	{
		m_inter_cluster_edges.push_back(edge);
	}

	return edge;
}

// output a debug message about the creation of an edge
void SEQUENTIAL_LEVEL::output_edge_debug_message
(
	LEVEL_NODE * source_level_node,
	LEVEL_NODE * sink_level_node
)
{
	assert(source_level_node && sink_level_node);
	debug("Adding an edge from cluster/delay level " << 
			source_level_node->get_cluster_number() << " " << 
			source_level_node->get_delay_level() << " to " << 
			sink_level_node->get_cluster_number() << " " << 
			sink_level_node->get_delay_level());
}


// adds a cluster's level nodes to m_delay_levels
//
// PRE: the level nodes have been created
//      we have previously added the level nodes 
//      of clusters 0..cluster_number-1 to m_delay_levels
// POST: the level nodes of cluster cluster_number have been a
//       added to m_delay_levels
//      
void SEQUENTIAL_LEVEL::add_level_nodes
(
	const LEVEL_NODES& level_nodes,
	const CLUSTER_NUMBER_TYPE& cluster_number
)
{
	DELAY_TYPE delay = 0;	
	LEVEL_NODE * level_node = 0;
	DELAY_LEVEL * delay_level = 0;

	// we have the correct number of level nodes
	assert(static_cast<signed>(level_nodes.size()) == m_max_comb_delay + 1);

	for (delay = 0; delay <= m_max_comb_delay; delay++)
	{
		level_node = level_nodes[delay];
		assert(level_node);

		delay_level = m_delay_levels[delay];
		assert(delay_level);
		
		// make sure we are adding level nodes in order of cluster 
		assert(static_cast<signed>(delay_level->size()) == cluster_number);	

		//debug("In cluster " << cluster_number << ".  Adding level node at delay level " << delay);

		delay_level->push_back(level_node);
	}
}


// Have we satisfied all specifications?
// 
// PRE: we have a solution
// RETURN: true if we have, false otherwise
//        
bool SEQUENTIAL_LEVEL::is_valid_solution()
{
	CLUSTER* cluster = 0;
	bool valid_solution =  true;
	CLUSTERS::iterator cluster_iter = m_clusters.begin(); 

	while (valid_solution && (cluster_iter != m_clusters.end()))
	{
		cluster = *cluster_iter;
		assert(cluster);

		valid_solution = cluster->is_valid_solution();

		cluster_iter++;
	}

	valid_solution = valid_solution && is_Comb_satisfied();

	return valid_solution;	
}


// Have we satisfied all specifications?
// 
// PRE: we have a solution
// RETURN: true if we have, false otherwise
//        
bool SEQUENTIAL_LEVEL::is_Comb_satisfied()
{
	bool inter_cluster_connection_matrix_satisfied = false;
	calculate_Comb();

	inter_cluster_connection_matrix_satisfied = (m_Comb == *m_Comb_spec);
	assert(inter_cluster_connection_matrix_satisfied);

	return inter_cluster_connection_matrix_satisfied;
}


// Calculate the current Comb matrix
//
// PRE: m_inter_cluster_edges contains all the inter-cluster edges in the circuit
// POST: m_Comb has been calculated
//
void SEQUENTIAL_LEVEL::calculate_Comb()
{
	CLUSTER_NUMBER_TYPE source_cluster_number, sink_cluster_number;
	EDGES::iterator edge_iter;
	EDGE * edge = 0;
	NUM_ELEMENTS weight = 0;

	m_Comb.clear();

	for (edge_iter = m_inter_cluster_edges.begin(); edge_iter != m_inter_cluster_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		source_cluster_number = edge->get_source_cluster_number();
		sink_cluster_number = edge->get_sink_cluster_number();

		assert(source_cluster_number != sink_cluster_number);

		weight = edge->get_weight();

		m_Comb(source_cluster_number, sink_cluster_number) += weight;
	}
}


// Creates an initial solution
// Iterates to improve the solution
// Solution is accepted if we have no node violations are could
// resolve the violations
//
// PRE: we have level nodes and super edges
//      the super edges been assigned the maximum possible weight
//
// POST: we have either found an acceptable solution or else
//       we have failed after trying a number of times
//
//
void SEQUENTIAL_LEVEL::find_solution()
{
	NUM_ELEMENTS attempts = 0;
	NUM_ELEMENTS max_attempts = 15;

	// try to find a solution
	bool success = try_to_create_delay_structure();

	while (! success && attempts < max_attempts)
	{
		attempts++;
		reset_graph();
		success = try_to_create_delay_structure();
	}

	if (attempts == max_attempts)
	{
		Fail("We tried the maximum number of times (" << max_attempts << 
			") to get a solution for " << m_circuit->get_name() << " and failed.");
	}
}


// try to solve the level graph structure problem.
//
// 1. Create an initial solution
// 2. Try to improve the solution by employing an iterative algorithm that selects 
//    certain edges as candidates for relocation and accepts or rejects proposed changes 
//    based on a cost function.
// 3. Post-process the solution to try and fix the any node violations
// 4. If we have no remaining node violations, create the edges that connect nodes to flip-flops
// 5. Call it a day. check sanity. pass the solution to the next stage.
//
//
bool SEQUENTIAL_LEVEL::try_to_create_delay_structure()
{
	bool success = true;
	m_display_too_many_inputs_msg = false;
	m_sanity_check_fatal = false;


	debug("===============================================\n");
	debug("\nStep 1: Creating the Delay Structure");
	debug(Sep);

	debug("Status: Creating the combinational delay structure first\n");

	if (! is_valid_solution())
	{
		finish_creating_initial_solution_to_combinational_delay_structure();

		simulated_anneal();

		m_display_too_many_inputs_msg = true;
		sanity_check_all();

		success = try_to_fix_graph();

		if (! success)
		{
			Warning("Reseting the delay structure and trying again.");
			return false;
		}
	}

	debugsep;
	is_valid_solution();
	//output_all_violations();

	sanity_check_all();

	debug("\nStatus: Finished creating the combinational delay structure\n");
	debug("Solution Quality for Delay Structure Creation");
	debug("---------------------------------------------");
	report_solution_quality();

	make_latched_to_dff_connections();

	m_sanity_check_fatal = true;
	sanity_check_all();

	debug("Status: Finished creating the delay structure\n");
	debug(Sep);

	return true;
}


// Finish creating the initial solution to the combinational delay structure.
//
//
// Part of creating the initial combinational delay structure.
// (see 4.1.1 of thesis)
// 
// Previously, in our four part series "Creating the Initial Combinational
// Delay Structure" we have done Parts 1 and 2.
// i.e.
// We have level nodes and super edges and the super edges have 
// been assigned the maximum possible weight.
// 
// We now present Steps 3 and 4....
//
// In Step 3, we randomly reduce weight on the inter-cluster edges to meet Comb, and in each cluster
// the inter-cluster input edge length distribution and inter-cluster output edge length 
// distribution while making sure that each level node has enough unit edges to define the 
// delay level of its nodes. The delay level of a node is defined if it has a unit edge from 
// a node in the delay level above. In Step 4, we randomly reduce weight on the intra-cluster 
// edges to meet the intra-cluster edge length distribution in each cluster while making sure 
// that each level node has enough unit edges to define the delay level of its nodes.
//
//
// PRE: inter-cluster super edges exist.
//      the total number of inter-cluster super edges may be greater than it should be.
//      intra-cluster super edges exist.
//      the total number of intra-cluster super edges may be greater than it should be.
//
// POST: we have an initial combinational delay structure solution
//
void SEQUENTIAL_LEVEL::finish_creating_initial_solution_to_combinational_delay_structure()
{
	DELETION_SITUATION deletion_situation = SEQUENTIAL_LEVEL::LEVEL_NODES_HAVE_SPARE_WEIGHT;
	CLUSTERS::iterator cluster_iter;
	CLUSTER * cluster = 0;
	CLUSTER_NUMBER_TYPE source_cluster_number, sink_cluster_number;
	m_number_loops_without_reduction = 0;
	NUM_ELEMENTS loop_limit = 30000,
	             spare_weight_limit = 1000;
	bool reduced_weight;


	// these functions need a rewrite
	// i don't like the algorithms to create the initial solution
	assert(new_gen);

	// Step 3: Reduce weight on inter-cluster edges
	//

	create_violation_matrix();
	if (! m_violation_matrix.is_positive())
	{
		// if we have modified the characterization (edge length for instance) without paying attention to the 
		// input and output shapes we can sometimes have less inter-cluster edges than we should.
		// This will sometimes cause the violation matrix to have non-zero entries
		
		// temporarily zero those negative entries so that we 
		// can remove the edges that should not be in the graph
		//
	
		// later on should we maybe not delete a few edges to make up for the negative entries?
		assert(new_gen);

		m_violation_matrix.zero_negative_entries();
	}


	// update whether or not we have edges that we can delete without 
	// violating a constraint or specification
	update_deletion_situation(deletion_situation);

	while (! m_violation_matrix.is_zero() && m_number_loops_without_reduction < loop_limit)
	{
		// get a source and sink cluster number with violations
		g_rand_number_gen->discrete_pmf(m_violation_matrix, source_cluster_number, sink_cluster_number);

		assert(m_violation_matrix(source_cluster_number, sink_cluster_number) > 0);

		assert(m_clusters[source_cluster_number]->get_cluster_number() == source_cluster_number);
		assert(m_clusters[sink_cluster_number]->get_cluster_number() == sink_cluster_number);

		// try to reduce weight on edges between the clusters to remove the violation
		reduced_weight = reduce_inter_cluster_edge_weight(m_clusters[source_cluster_number], 
													m_clusters[sink_cluster_number], 
													deletion_situation);
		// were we successful?
		if (reduced_weight)
		{
			// yes we were, just reduce the weight by 1	
			m_violation_matrix(source_cluster_number, sink_cluster_number) -= 1;

			// update our delation situation
			update_deletion_situation(deletion_situation);

			m_number_loops_without_reduction = 0;
		}
		else
		{
			m_number_loops_without_reduction++;
		}
	
		// we are not reducing weight on edges
		// tell the functions that they can reduce weight on edges that do not have spare weight
		if (m_number_loops_without_reduction > spare_weight_limit && 
			deletion_situation == LEVEL_NODES_HAVE_SPARE_WEIGHT)
		{
			deletion_situation = LEVEL_NODES_HAVE_WEIGHT;
		}
	}

	if (! m_violation_matrix.is_zero())
	{
		reduce_edges_to_meet_total_number_of_edges();
	}

	// create the final violation matrix
	create_violation_matrix();

	debug("Comb\n" << m_Comb);
	debug("Spec of Comb\n" << *m_Comb_spec);
	debug("Violation matrix (Comb- Spec of Comb)\n" << m_violation_matrix);

	assert(m_violation_matrix == m_Comb - *m_Comb_spec);



	// Step 4: Reduce weight on intra-cluster edges
	//
	// now sculpt with the cluster to satisfy the properties of each cluster
	for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);

		reduce_intra_cluster_edge_weight(cluster);
	}	

	create_edge_length_matrix();
}

//
// Create the violation matrix = current Comb - spec of Comb
//
// PRE:  combinational super edges have been created
// POST: m_violation_matrix has been created
// 		 m_Comb_spec_sum has been calculated
//
void SEQUENTIAL_LEVEL::create_violation_matrix()
{
	calculate_Comb();

	m_Comb_spec_sum = static_cast<COST_TYPE>(m_Comb_spec->get_absolute_sum());

	m_violation_matrix = m_Comb - *m_Comb_spec;

	//debug("Inter cluster connection matrix\n" << m_Comb);
	//debug("Spec cluster connext matrix\n" << *m_Comb_spec);
	//debug("Violation matrix\n" << m_violation_matrix);
}


// Reduce weight on the inter-cluster edges between the souce cluster
// and the sink cluster
//
// a good type of edge to reduce is one that both clusters
// have an excess number of in terms of their edge-length distributions
// try to find and delete such an edge
//
// PRE: source_cluster and sink_cluster are valid
// 		there are too many connections between the source and sink clusters
//      There are edges we can delete in the circuit
//
// POST: if successful we have reduced the weight on an edge between the source cluster
//       and the sink cluster
//
// RETURNS: success on reducing weight
//
bool SEQUENTIAL_LEVEL::reduce_inter_cluster_edge_weight
(
	CLUSTER* source_cluster,
	CLUSTER* sink_cluster,
	const DELETION_SITUATION & deletion_situation
)
{
	bool reduced_weight = false;
	DELAY_TYPE delay = 0;
	LEVEL_NODE * source_level_node = 0;
	LEVEL_NODE * sink_level_node = 0;
	EDGE * edge = 0;
	NUM_ELEMENTS number_edges_to_delete = 0;
	LENGTH_TYPE edge_length = 0;


	// get the difference between the current edge length violations and their spec
	DISTRIBUTION inter_cluster_output_edge_length_violations = 
		get_inter_cluster_output_edge_length_violations(source_cluster);
	DISTRIBUTION inter_cluster_input_edge_length_violations = 
		get_inter_cluster_input_edge_length_violations(sink_cluster);

	// they should be the same size
	assert(inter_cluster_input_edge_length_violations.size() == 
		   inter_cluster_output_edge_length_violations.size());


	// do any of our super edges in the circuit have weight to spare?
	bool spare_edges = (deletion_situation == SEQUENTIAL_LEVEL::LEVEL_NODES_HAVE_SPARE_WEIGHT);
	assert(deletion_situation != SEQUENTIAL_LEVEL::NO_EDGE_CAN_BE_DELETED);

	// try to find an edge_length that both clusters have an excess of.
	edge_length = choose_edge_length_with_violation(inter_cluster_input_edge_length_violations, 
									inter_cluster_output_edge_length_violations);

	if (edge_length < 0)
	{
		// we were not able to find such an edge
		reduced_weight = false;
		return reduced_weight;
	}
			
			
	// we have an excess number of these edges
	assert(inter_cluster_input_edge_length_violations[edge_length] > 0 && 
		   inter_cluster_output_edge_length_violations[edge_length] > 0);

	// randomly pick a delay level
	delay = (short) g_rand_number_gen->random_number(m_max_comb_delay-edge_length);
	assert((delay+edge_length) >= 0 && (delay+edge_length) <= m_max_comb_delay);
	

	// find the level nodes that could harbour such an edge
	source_level_node = source_cluster->Level[delay];
	sink_level_node = sink_cluster->Level[delay+edge_length];
	edge = sink_level_node->find_edge_with_source(source_level_node);
	assert(edge);

	// can we find an edge to delete between the level nodes with the edge length?
	number_edges_to_delete = get_number_edges_that_can_be_deleted(source_level_node, sink_level_node, 
											edge, inter_cluster_output_edge_length_violations[edge_length],
											inter_cluster_input_edge_length_violations[edge_length], spare_edges);

	if ( number_edges_to_delete > 0)
	{
		// yes 
	
		assert(source_level_node->get_inter_cluster_output_edge_lengths()[edge_length] > 0);
		assert(sink_level_node->get_inter_cluster_input_edge_lengths()[edge_length] > 0);
		assert(edge->get_weight() > 0);
		if (edge->get_length() == 1)
		{
			assert(sink_level_node->get_delay_difference() > 0);
		}

		// just delete one edge
		number_edges_to_delete = 1;
		edge->add_weight(-number_edges_to_delete);

		sink_level_node->sanity_check();

		reduced_weight = true;
		inter_cluster_input_edge_length_violations[edge_length] -= number_edges_to_delete;
		inter_cluster_output_edge_length_violations[edge_length] -= number_edges_to_delete;
	}

	return reduced_weight;
}



// we want to reduce weight between two clusters.
// a good type of edge to reduce is one that both clusters
// have an excess number of in terms of their edge-length distributions
// try to find such an edge
//
//
// PRE:  Inter_cluster_input_edge_length_violations contains 
//         the difference between current and spec for sink cluster
//       Inter_cluster_output_edge_length_violations contains 
//         the difference between current and spec for source cluster
// RETURN: the edge length where we could find an excessive number of edges
//          in both distribution
//          If we couldn't find such an edge length then we return -1 
//
LENGTH_TYPE SEQUENTIAL_LEVEL::choose_edge_length_with_violation
(
	const NUM_ELEMENTS_VECTOR & inter_cluster_input_edge_length_violations,
	const NUM_ELEMENTS_VECTOR & inter_cluster_output_edge_length_violations
) const
{
	LENGTH_TYPE edge_length = 0;
	NUM_ELEMENTS edge_length_size = inter_cluster_input_edge_length_violations.size();
	NUM_ELEMENTS_VECTOR edge_length_violation_pmf(edge_length_size, 0);

	assert(inter_cluster_input_edge_length_violations.size() == 
		   inter_cluster_output_edge_length_violations.size());

	// for each edge length find if we can find 
	// an edge length that has excessive number of inter-cluster input edges
	// and excessive number of inter-cluster output edges
	for (edge_length = 0; edge_length < edge_length_size; edge_length++)
	{
		if (inter_cluster_input_edge_length_violations[edge_length] > 0 && inter_cluster_output_edge_length_violations[edge_length] > 0)
		{
			// we can eliminate an edge at this edge length
			edge_length_violation_pmf[edge_length] 
				= MIN(inter_cluster_input_edge_length_violations[edge_length],
					  inter_cluster_output_edge_length_violations[edge_length]);
		}
		else
		{
			edge_length_violation_pmf[edge_length] = 0;
		}
	}


	// if we have no matches at all just return
	if (accumulate(edge_length_violation_pmf.begin(), edge_length_violation_pmf.end(), 0) == 0)
	{
		// we did not have an edge length violation match
		// return -1 to signify that we couldn't find an edge length
		return -1;
	}

	// assert: we have some edges lengths that have excess number of inter-cluster input
	// and inter-cluster output edges


	// select one of the edge lengths
	edge_length = (short) g_rand_number_gen->discrete_pmf(edge_length_size, edge_length_violation_pmf);
	assert(edge_length >= 0 && edge_length <= edge_length_size - 1);



	// we are able to return a valid length

	return edge_length;
}

//
// PRE: source_level_node, sink_level_node, and edge are valid
//      source_violations is the number of excessive edges 
//		  of the right length in the source cluster
//		sink_violations is the number of excessive edges 
//		  of the right lenght in the sink cluster
//		only_look_for_spare is true if we want the function 
//		  to only look for edges that the level nodes could spare
//		  and false if we want the function to be more ruthless
//      
// RETURN: the number of edges that can be deleted from the edge
NUM_ELEMENTS SEQUENTIAL_LEVEL::get_number_edges_that_can_be_deleted
(
	const LEVEL_NODE * source_level_node,
	const LEVEL_NODE * sink_level_node,
	const EDGE * edge,
	const NUM_ELEMENTS& source_violations,
	const NUM_ELEMENTS& sink_violations,
	const bool& only_look_for_spare
) const
{
	NUM_ELEMENTS spare_edges = 0;
	assert(edge && source_level_node && sink_level_node);
	LENGTH_TYPE edge_length = edge->get_length();

	spare_edges = edge->get_weight();

	if (only_look_for_spare)
	{
		spare_edges = MIN(spare_edges, source_level_node->get_output_to_spare());
		spare_edges = MIN(spare_edges, sink_level_node->get_input_to_spare());
	}
	spare_edges = MIN(spare_edges, source_violations);
	spare_edges = MIN(spare_edges, sink_violations);

	if (edge_length == 1)
	{
		// need to check the delay
		spare_edges = MIN(spare_edges, sink_level_node->get_delay_to_spare());
	}
	
	// check to make sure everything is ok
	assert(spare_edges >= 0);


	return spare_edges;
}


// 
// PRE: cluster is valid
// RETURNS: current inter-cluster input edge length distribution - 
//          spec    inter-cluster input edge length distribution
//
DISTRIBUTION SEQUENTIAL_LEVEL::get_inter_cluster_input_edge_length_violations
(
	CLUSTER * cluster
)
{
	assert(cluster);
	CLUSTER_SPEC * comb_spec = cluster->get_comb_spec();
	assert(comb_spec);
	LENGTH_TYPE edge_length;
	DISTRIBUTION edge_length_violations = cluster->get_inter_cluster_input_edge_lengths();
	DISTRIBUTION spec_edge_length = comb_spec->get_inter_cluster_input_edge_lengths();

    assert(edge_length_violations.size() == spec_edge_length.size());

	if (m_display_too_many_inputs_msg)
	{
		cout << "Existing inter_cluster_input_edge_length\n";
		copy(edge_length_violations.begin(),edge_length_violations.end(),
				ostream_iterator<NUM_ELEMENTS>(cout, " "));
		cout << "\nSpec inter_cluster_input_edge_length\n";
		copy(spec_edge_length.begin(), spec_edge_length.end(), 
				ostream_iterator<NUM_ELEMENTS>(cout, " ")); 
	}


	for (edge_length = 0; edge_length < static_cast<signed>(edge_length_violations.size()); edge_length++)
	{
		// sometimes edge_length_violations[edge_length] < spec_edge_length[edge_length]
		// contrary to expectation.
		//
		// the above sometimes happens when we change the specification because
		// the new input and output shapes cannot accomodate the edges that want to be inserted into
		// the graph
		//
		// therefore because of the above we cannot use the assert below
		//assert(edge_length_violations[edge_length] >= spec_edge_length[edge_length]);

		edge_length_violations[edge_length] -= spec_edge_length[edge_length];
	}


	if (m_display_too_many_inputs_msg)
	{
		cout << "\ninput edge length violations for cluster " << cluster->get_cluster_number() << "\n";
		copy(edge_length_violations.begin(), edge_length_violations.end(), 
				ostream_iterator<NUM_ELEMENTS>(cout, " ")); cout << endl;
	}

	return edge_length_violations;
}
// 
// PRE: cluster is valid
// RETURNS: current inter-cluster output edge length distribution - 
//          spec    inter-cluster output edge length distribution
//
DISTRIBUTION SEQUENTIAL_LEVEL::get_inter_cluster_output_edge_length_violations
(
	CLUSTER * cluster
)
{
	assert(cluster);
	CLUSTER_SPEC * comb_spec = cluster->get_comb_spec();
	assert(comb_spec);
	LENGTH_TYPE edge_length;
	DISTRIBUTION edge_length_violations = cluster->get_inter_cluster_output_edge_lengths();
	DISTRIBUTION spec_edge_length = comb_spec->get_inter_cluster_output_edge_lengths();

    assert(edge_length_violations.size() == spec_edge_length.size());

	if (m_display_too_many_inputs_msg)
	{
		cout << "Existing inter_cluster_output_edge_length\n";
		copy(edge_length_violations.begin(),edge_length_violations.end(),
				ostream_iterator<NUM_ELEMENTS>(cout, " "));

		cout << "\nSpec inter_cluster_output_edge_length\n";
		copy(spec_edge_length.begin(), spec_edge_length.end(), 
				ostream_iterator<NUM_ELEMENTS>(cout, " ")); 

		cout << endl;
	}


	for (edge_length = 0; edge_length < static_cast<signed>(edge_length_violations.size()); edge_length++)
	{
		// sometimes edge_length_violations[edge_length] < spec_edge_length[edge_length]
		// contrary to expectation.
		//
		// the above sometimes happens when we change the specification because
		// the new input and output shapes cannot accomodate the edges that want to be inserted into
		// the graph
		//
		// therefore because of the above we cannot use the assert below
		//assert(edge_length_violations[edge_length] >= spec_edge_length[edge_length]);

		edge_length_violations[edge_length] -= spec_edge_length[edge_length];
	}


	if (m_display_too_many_inputs_msg)
	{
		cout << "output edge length violations for cluster " << cluster->get_cluster_number() << "\n";
		copy(edge_length_violations.begin(), edge_length_violations.end(), 
				ostream_iterator<NUM_ELEMENTS>(cout, " ")); cout << endl;
	}

	return edge_length_violations;
}
// 
// PRE: cluster is valid
// RETURNS: current intra-cluster edge length distribution - 
//          spec    intra-cluster edge length distribution
//
DISTRIBUTION SEQUENTIAL_LEVEL::get_intra_cluster_edge_length_violations
(
	CLUSTER * cluster
)
{
	assert(cluster);
	CLUSTER_SPEC * comb_spec = cluster->get_comb_spec();
	assert(comb_spec);
	LENGTH_TYPE edge_length;
	DISTRIBUTION edge_length_violations = cluster->get_intra_cluster_edge_lengths();
	DISTRIBUTION spec_edge_length = comb_spec->get_intra_cluster_edge_lengths();

    assert(edge_length_violations.size() == spec_edge_length.size());

	if (m_display_too_many_inputs_msg)
	{
		cout << "Existing intra_cluster_edge_length\n";
		copy(edge_length_violations.begin(), edge_length_violations.end(),
				ostream_iterator<NUM_ELEMENTS>(cout, " "));
		cout << "\nSpec intra_cluster_edge_length\n";
		copy(spec_edge_length.begin(), spec_edge_length.end(), 
				ostream_iterator<NUM_ELEMENTS>(cout, " ")); 
	}

	for (edge_length = 0; edge_length < static_cast<signed>(edge_length_violations.size()); edge_length++)
	{
		// sometimes edge_length_violations[edge_length] < spec_edge_length[edge_length]
		// contrary to expectation.
		//
		// the above sometimes happens when we change the specification because
		// the new input and output shapes cannot accomodate the edges that want to be inserted into
		// the graph
		//
		// therefore because of the above we cannot use the assert below
		//assert(edge_length_violations[edge_length] >= spec_edge_length[edge_length]);
		//assert(edge_length_violations[edge_length] >= spec_edge_length[edge_length]);

		edge_length_violations[edge_length] -= spec_edge_length[edge_length];
	}


	if (m_display_too_many_inputs_msg)
	{
		cout << endl;
		cout << "edge length violations\n";
		copy(edge_length_violations.begin(), edge_length_violations.end(), 
				ostream_iterator<NUM_ELEMENTS>(cout, " ")); cout << endl;
	}

	return edge_length_violations;
}

//
// PRE: cluster is valid
// POST: we either have no intra-cluster edge violations OR 
//       we have removed as many edges as we could 
//       from this cluster
//
//		 We have issued a warning if we have less intra-cluster edges than we should.
//		 This can occur due to input-output shape/edge length conflict that 
//		 can prevent us from  adding enough intra-cluster edges during 
//       create_super_edges_between_level_nodes()
//		 
void SEQUENTIAL_LEVEL::reduce_intra_cluster_edge_weight
(
	CLUSTER * cluster
)
{
	assert(cluster);
	NUM_ELEMENTS number_of_violations = 0;
	const NUM_ELEMENTS number_edges_to_delete = 1;
	EDGE * edge = 0;
	LENGTH_TYPE edge_length = 0;

	DISTRIBUTION edge_length_violations = get_intra_cluster_edge_length_violations(cluster);

	// look over all the edge lengths 
	for (edge_length = static_cast<signed>(edge_length_violations.size()-1); edge_length > 0; edge_length--)
	{
		number_of_violations = edge_length_violations[edge_length];

		// if because of input-output shape/edge length conflict prevented us
		// from adding enough edges issue a warning for now
		if (number_of_violations < 0)
		{
			// make the number of violations positive
			number_of_violations *= -1; 

			Warning("We have " << number_of_violations << " less edges of length " 
					<< edge_length << " in cluster " << cluster->get_cluster_number() 
					<< " than we should have due to shape/edge length conflict.");
			number_of_violations = 0;

			assert(new_gen);
			
			// we have less edges that we should have.
			// this is because of input-output shape/edge length distribution mismatch
			//
			// we have to decide whether we want to violate the input shape or the edge length 
			// distribution
			// for right now do nothing, ie. violate the edge length distributions
			/*
			 * unused code to force edge addition in the cluster
			add_weight_to_edges_in_cluster(cluster, edge_length, number_of_violations);

			if (number_of_violations > 0)
			{
				Warning("After trying to force edge addition we still have " 
						<< number_of_violations << " less edges of length " 
						<< edge_length << " in cluster " << cluster->get_cluster_number());
			}
			number_of_violations = 0;
			*/
		}


		// while we have too many edges at this edge length
		while (number_of_violations > 0)
		{
			// try to get an edge to delete
			edge = get_edge_we_can_delete(cluster, edge_length);

			if (edge)
			{
				//we have such an edge

				// if the edge is length 1 make sure
				// we are not going to violate the number of nodes 
				// needed to define delay
				if (edge->get_length() == 1) assert(edge->get_sink_level_node()->get_delay_difference() > 0);

				// reduce weight on the edge
				edge->add_weight(-number_edges_to_delete);

				// is everything still ok with the sink node?
				edge->get_sink_level_node()->sanity_check();

				number_of_violations -= number_edges_to_delete;
				edge_length_violations[edge_length] -= number_edges_to_delete;
			}
			else
			{
				// we couldn't not find an edge to delete.  need to continue
				// with an extra edge
				number_of_violations = 0;
			}

		}

		if (number_of_violations != 0)
		{
			//debug("We still have " << number_of_violations << " violations ");
			output_spare(cluster);
		}
		assert(number_of_violations == 0);

	}
}
//
// PRE: cluster is valid
// POST: we have outputed for this cluster
//       how much each delay level can spare
//
void SEQUENTIAL_LEVEL::output_spare
(
	CLUSTER * cluster
)
{
	assert(cluster);

	debug("Looking at cluster " << cluster->get_cluster_number());
	LEVEL_NODE * level_node = 0;
	DELAY_TYPE delay = 0;

	for (delay = 0; delay <= m_max_comb_delay; delay++)
	{
		debug("\tLooking at delay " << delay);

		level_node = cluster->Level[delay];

		debug("\t\tinputs difference " << level_node->get_input_difference());
		debug("\t\toutput difference " << level_node->get_output_difference());
		debug("\t\tdelay difference " << level_node->get_delay_difference());
	
		level_node->report_edge_lengths();	
		
	}
}


// check the sanity fo the cluster
//
// PRE: cluster is valid
// POST: If we have level nodes node_violations i.e. level nodes with:
//         - too many inputs
//		   - not enough inputs
//		   - not enough outputs
//		 We have warned the user if 
//		 m_display_too_many_inputs_msg is set to true
//
//		 If we found node violations and m_sanity_check_fatal is 
//		 true we have Failed and exited the program
//
//       We have also sanity checked the clusters and the level nodes within
//       the clusters.
//
void SEQUENTIAL_LEVEL::sanity_check
(
	CLUSTER * cluster
) const
{
	debugif(DSANITY_CHECK, "Sanity checking ");
	assert(cluster);
	CLUSTER_SPEC * comb_spec = cluster->get_comb_spec();
	assert(comb_spec);
	LEVEL_NODE * level_node = 0;
	DELAY_TYPE delay = 0;
	NUM_ELEMENTS too_many_inputs = 0,
				 level_nodes_with_too_many_inputs = 0,
				 level_nodes_that_need_inputs = 0,
				 level_nodes_that_need_outputs = 0;
	bool warning_found = false;

	for (delay = 0; delay <= m_max_comb_delay; delay++)
	{
		//debug("\tLooking at delay " << delay);

		level_node = cluster->Level[delay];
		assert(level_node);

		level_node->sanity_check();

		if (delay > 0)
		{
			too_many_inputs = level_node->get_excess_number_of_inputs();
			if (too_many_inputs > 0 && m_display_too_many_inputs_msg)
			{
				Verbose("We have " << too_many_inputs << " too many inputs inputs to cluster "
						<< cluster->get_cluster_number() << " level " << delay);
				level_nodes_with_too_many_inputs++;
				warning_found = true;
			}
		}
		if (level_node->get_nOutputs_still_needed() > 0)
		{
			Verbose("We have " << level_node->get_nOutputs_still_needed() << " too few outputs to cluster "
					<< cluster->get_cluster_number() << " level " << delay);
			level_nodes_that_need_outputs++;
			warning_found = true;
		}
		if (level_node->get_nInputs_still_needed()  > 0)
		{
			Verbose("We have " << level_node->get_nInputs_still_needed() << " too few inputs to cluster "
					<< cluster->get_cluster_number() << " level " << delay);
			level_nodes_that_need_inputs++;
			warning_found = true;
		}
	}

	if (warning_found && m_display_too_many_inputs_msg)
	{
		debug("Have to post-process the combinational delay structure in cluster "
				<< cluster->get_cluster_number());
		
		if (level_nodes_with_too_many_inputs > 0)
		{
			Verbose("  Have " << level_nodes_with_too_many_inputs << " level nodes with too many inputs");
		}
		if (level_nodes_that_need_inputs > 0)
		{
			Verbose("  Have " << level_nodes_that_need_inputs << " level nodes with too few inputs");
		}
		if (level_nodes_that_need_outputs > 0)
		{
			Verbose("  Have " << level_nodes_that_need_outputs << " level nodes with too few outputs");
		}
	}

	cluster->sanity_check();

	if (warning_found && m_sanity_check_fatal)
	{
		Fail("Couldn't fix errors.");
	}

}


// 1. randomly try to find a delay level
// 2. iterate sequentially through delay levels looking for an delay level with spare
// 3. iterator sequentially through delay levels looking for a delay level that has 
//    a edge we can delete.
// we should have found an edge by this time
// 
// PRE: cluster is valid
// POST: edge length is the edge_length we want to 
//       delete it the cluster
// RETURN: if sucessful, the edge we can delete, otherwise 0
//
EDGE * SEQUENTIAL_LEVEL::get_edge_we_can_delete
(
	CLUSTER * cluster,
	const LENGTH_TYPE& edge_length
)
{
	assert(edge_length >= 0 && edge_length <= m_max_comb_delay);

	DELAY_TYPE source_delay = (short) g_rand_number_gen->random_number(m_max_comb_delay-edge_length),
	           sink_delay = source_delay + edge_length,
	           old_source_delay = 0;
	LEVEL_NODE * source_level_node = 0,
	           * sink_level_node = 0;
	EDGE * edge = 0;
	bool have_gone_through_all_edges = false;
	
	assert(	source_delay >= 0 && source_delay <= m_max_comb_delay &&
			sink_delay >= 0 && sink_delay <= m_max_comb_delay);

	//debug("randomly trying source delay level " << source_delay);
	edge = cluster->get_edge(source_delay, edge_length);

	if (! edge->has_spare_weight())
	{
		// the edge and the level nodes did not had some weight to spare
		// try to increase delay until we find an edge and level node that does
		// if no edge can be found, just try to find an edge that we can remove some weight from
		
		old_source_delay = source_delay;
		have_gone_through_all_edges = false;

		while (! edge->has_spare_weight() && ! have_gone_through_all_edges)
		{
			//debug("delay level " << source_delay << " didn't have any edges to spare.");
			source_delay = (source_delay +1) % (m_max_comb_delay-edge_length+1);

			edge = cluster->get_edge(source_delay, edge_length);
			assert(edge);

			have_gone_through_all_edges = (source_delay == old_source_delay);
		}

		if (have_gone_through_all_edges)
		{
			// we didn't find a delay level with spare
			// now just try to edge we can remove some weight from 
		
			have_gone_through_all_edges = false;

			while (! edge->can_reduce_weight_on_edge() && ! have_gone_through_all_edges)
			{
				//debug("delay level " << source_delay << " didn't have any weight at all.");
				source_delay = (source_delay +1) % (m_max_comb_delay-edge_length+1);
				edge = cluster->get_edge(source_delay, edge_length);
				assert(edge);

				have_gone_through_all_edges = (source_delay == old_source_delay);
			}

			// we should have found an edge by 
			// now. 
			// we get here very occasionally without having found an
			// edge with edge_length 1
			assert(new_gen);
			if (have_gone_through_all_edges)
			{
				Warning("Tried and fail to delete from cluster " << cluster->get_cluster_number() 
						<< " a length " << edge_length << "...continuing with an extra edge");
				if (g_options->is_verbose())
				{
					output_all_violations();
				}

				return 0;
			}	
		}
		else
		{
			//debug("delay level " << source_delay << " had edges to spare.");
		}
	}
	else
	{
		//debug("We have spare at delay level " << source_delay);
	}


	// just checking to make sure everything is right
	sink_delay = source_delay+edge_length;

	assert(	source_delay >= 0 && source_delay <= m_max_comb_delay &&
			sink_delay >= 0 && sink_delay <= m_max_comb_delay);

	source_level_node = edge->get_source_level_node();
	sink_level_node = edge->get_sink_level_node();
	assert(source_level_node && sink_level_node);


	if (edge_length == 1)
	{
		assert(sink_level_node->get_delay_difference() > 0);
	}
	assert(source_level_node->get_intra_cluster_output_edge_lengths()[edge_length] > 0);
	assert(sink_level_node->get_intra_cluster_input_edge_lengths()[edge_length] > 0);
	assert(edge->get_weight() > 0);

	return edge;
}


// a first attempt at trying to solve
// the input-output shape/edge length conflict by 
// violating the input-output shape to add more weight to edges
//
// haven't tested this function
//
// PRE: cluster is valid
//      edge_length is the length of edge we want to add to the cluster
//      number_of_violations is the number of edges that need to be added
// POST: we have tried to add edges to the graph in order to 
//       to eliminate the node violations
//
void SEQUENTIAL_LEVEL::add_weight_to_edges_in_cluster
(
	CLUSTER * cluster, 
	const LENGTH_TYPE& edge_length, 
	NUM_ELEMENTS & number_of_violations
)
{
	EDGE * edge = 0;
	DELAY_TYPE sink_delay_level = 0,
			   source_delay = 0,
			   max_delay = 0;
	LEVEL_NODE * source_level_node = 0,
			   * sink_level_node = 0;

	NUM_ELEMENTS max_edge_weight = 0,
				 old_max_edge_weight = 0,
				 edges_to_add = 0;

					
	for (source_delay = 0; source_delay < max_delay; source_delay++)
	{
		edge = cluster->get_edge(source_delay, edge_length);
		assert(edge);

		sink_delay_level = source_delay + edge_length;
		assert(sink_delay_level >= 0 && sink_delay_level <= max_delay);

		sink_level_node = cluster->Level[sink_delay_level];
		sink_level_node = cluster->Level[source_delay];
	
		old_max_edge_weight = edge->get_max_weight();
		
		// make sure we get no double connections between individual nodes
		max_edge_weight = source_level_node->get_nNodes() * sink_level_node->get_nNodes();
		

		// set all edges to the maximum so that we don't
		// constrain moves in this cluster by just making one 
		// edge have the capicity to accept edges
		if (max_edge_weight > old_max_edge_weight)
		{
			edge->set_max_weight(max_edge_weight);
		}

		// now see if we can add further edges
		edges_to_add = MIN(max_edge_weight - old_max_edge_weight, number_of_violations);
		edge->set_weight_and_update_level_nodes(edges_to_add);

		number_of_violations -= edges_to_add;
		assert(number_of_violations >= 0);
	}
}

// RETURNS: the unnormalized cost
//
COST_TYPE SEQUENTIAL_LEVEL::get_unnormalized_cost() const
{
	// get_cost is unnormalized for now
	return get_cost();
}

// RETURNS: for now it returns the unnormalized cost
//
COST_TYPE SEQUENTIAL_LEVEL::get_cost() const
{
	CLUSTERS::const_iterator cluster_iter;
	COST_TYPE cost_Level_Node = 0.0,
			  cost_Edge_Length = 0.0,
			  cost_Comb = 0.0,
			  cost = 0.0;
	CLUSTER * cluster = 0;

	// find costs
	for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);

		// cost_Level_Node = cost_Level_Shape + cost_Problem_Node
		cost_Level_Node += cluster->get_level_node_cost();

		cost_Edge_Length += cluster->get_edge_length_cost();
	}

	cost_Comb = m_violation_matrix.get_absolute_sum();
	
	cost = m_alpha*cost_Comb + m_beta*cost_Edge_Length + m_gamma*cost_Level_Node;

	return cost;
}

// RETURNS: the cost just considering the clusters involved in the move 
//          + the violation matrix cost
//
COST_TYPE SEQUENTIAL_LEVEL::get_cost_of_clusters_involved
(
	const MOVE& move 
 ) const
{
	COST_TYPE cost_Level_Node = 0.0,
			  cost_Edge_Length = 0.0,
			  cost_Comb = 0.0,
			  cost = 0.0;

	EDGE * source_edge = move.get_source_edge();
	EDGE * sink_edge = move.get_sink_edge();

	assert(source_edge && sink_edge);

	CLUSTER_NUMBER_TYPE 
		source_edge_source_cluster = source_edge->get_source_cluster_number(),
		source_edge_sink_cluster   = source_edge->get_sink_cluster_number(),
		sink_edge_source_cluster = sink_edge->get_source_cluster_number(),
		sink_edge_sink_cluster   = sink_edge->get_sink_cluster_number();
							
	// cost_Level_Node = cost_Level_Shape + cost_Problem_Node

	cost_Level_Node += m_clusters[source_edge_source_cluster]->get_level_node_cost(); 
	cost_Edge_Length += m_clusters[source_edge_source_cluster]->get_edge_length_cost();


	cost_Level_Node += m_clusters[source_edge_sink_cluster]->get_level_node_cost(); 
	cost_Edge_Length += m_clusters[source_edge_sink_cluster]->get_edge_length_cost();

	if (source_edge_source_cluster != sink_edge_source_cluster)
	{
		cost_Level_Node += m_clusters[sink_edge_source_cluster]->get_level_node_cost(); 
		cost_Edge_Length += m_clusters[sink_edge_source_cluster]->get_edge_length_cost();
	}

	if (source_edge_sink_cluster != sink_edge_sink_cluster)
	{
		cost_Level_Node += m_clusters[sink_edge_sink_cluster]->get_level_node_cost(); 
		cost_Edge_Length += m_clusters[sink_edge_sink_cluster]->get_edge_length_cost();
	}

	cost_Comb = m_violation_matrix.get_absolute_sum();
	
	cost = m_alpha*cost_Comb + m_beta*cost_Edge_Length + m_gamma*cost_Level_Node;

	return cost;
}
//
// 
// RETURNS: return the normalized cost
//          
//          NOTE: I don't know how to normalize cost_Problem node 
//          	  with the congestion factors multiplying 
//                the cost_Problem_node to very large numbers
//
COST_TYPE SEQUENTIAL_LEVEL::get_normalized_cost() const
{
	CLUSTERS::const_iterator cluster_iter;
	COST_TYPE cost_Level_Node = 0.0,
			  cost_Edge_Length = 0.0,
			  cost_Comb = 0.0,
			  cost_cluster = 0.0,
			  cost = 0.0;
	CLUSTER * cluster = 0;

	assert(g_options);
	double phi = g_options->get_delay_structure_phi();
	double lamda = g_options->get_delay_structure_lamda();
	double nEdges = static_cast<double>(m_circuit->get_nEdges());

	for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);

		// cost_Level_Node = cost_Level_Shape + cost_Problem_Node
		cost_Level_Node += cluster->get_level_node_cost();

		cost_Edge_Length += cluster->get_edge_length_cost();
	}


	cost_cluster = phi* cost_Level_Node/nEdges + (1-phi)* cost_Edge_Length/nEdges;

	cost_Comb = m_violation_matrix.get_absolute_sum();
	cost_Comb = cost_Comb/nEdges;

	cost = (1-lamda)*cost_Comb + lamda*cost_cluster;

	return cost;
}


// get the changed cost of the move
//
// PRE: move is valid
// POST: nothing
// RETURN: changed_cost
//
COST_TYPE SEQUENTIAL_LEVEL::get_changed_cost
(
	MOVE & move,
	const COST_TYPE& current_cost
) 
{
	move.make_temporary_move(m_violation_matrix);

	COST_TYPE changed_cost = get_cost_of_clusters_involved(move);

	move.unmake_temporary_move(m_violation_matrix);

	//changed_cost -= current_cost;
	changed_cost -= get_cost_of_clusters_involved(move);

	//debug("Changed cost " << changed_cost);
	return changed_cost	;
}

//	After creating the initial solution, we employ an iterative 
//	algorithm that selects certain edges as candidates for relocation 
//	and accepts or rejects proposed changes based on a cost function. 
//	In the following paragraphs, we will describe the cost function, 
//	how edges are selected and modified ("moves"), the overall structure 
//	of the algorithm, and the performance of the algorithm. The cost function 
//	of the iterative algorithm is as follows: 
//  
//     cost = cost_Comb + cost_Edge_Length + cost_Level_Shape cost_Problem_Node 
//
//  Here, COST_Comb measures the absolute difference between the current Comb 
//  and its specification, cost_Edge_Length measures the absolute difference 
//  between the current edge length distributions and their specifications, 
//  cost_Level_Shape measures the absolute difference between the current input 
//  and output shapes and their specifications with congestion factors multiplying 
//  this cost at each level node to penalize level nodes that have too many inputs 
//  or outputs. Finally, cost_Problem_Node measures the number of node violations that 
//  would be forced to be made during final edge assignment because of the number of 
//  edges input or output of a level node. 
//
//  The algorithm selects super edges to change in the graph as follows: First, we randomly 
//  select a super edge to move with edges that contribute more to the COST_Level_Shape having 
//  a higher probability of being selected.  We weight the probability of an edge being 
//  selected where we give higher selection weights to super edges between source level 
//  nodes that have too many outputs and too many inputs and give lower selection weights 
//  to super edges that have source level nodes with not enough outputs or sink level nodes 
//  with not enough inputs. We also add a small bias to the selection weight of each 
//  super edge to ensure that each edge has a small probability of being chosen.
//
//	Second, after selecting an edge to move, we randomly choose a new destination for 
//	the edge from the possible super-edges that can accept the edge with its given length. 
//	A super edge can accept the individual edge to move if it has the same length and if 
//	the super edge is not already at its maximum possible weight (the weight we calculated 
//	above during the creation of the initial solution).
//
//	To improve the solution, we employ a simulated annealing-like algorithm. In the algorithm, 
//	all valid moves that improve the cost function are accepted, and valid bad moves 
//	are accepted with a probability that decreases exponentially with an increasing change 
//	of cost. A move is valid if, after the move, the sink level nodes involved will still 
//	have enough unit edges to correctly define the delay level of their individual nodes. 
//	The algorithm terminates when the cost is zero or after the number of moves attempted 
//	is fifty times the number of edges in the graph. During iteration, to help resolve 
//	contentions for level nodes among the edges every |E(G)|/2 moves we increase the 
//	congestion costs of level nodes that have too many inputs or outputs and decrease 
//	the congestion costs of level nodes that no longer have too many inputs or outputs. 
//	Also during iteration, after attempting moves equal to forty times the number of edges, 
//	we decrease the temperature to perform more of a greedy optimization.
//
void SEQUENTIAL_LEVEL::simulated_anneal()
{
	assert(g_options);


	MOVE move;
	COST_TYPE cost = get_cost(),
			  changed_cost = 0.0,
			  lowest_cost = cost;
	NUM_ELEMENTS loops = 0;
	NUM_ELEMENTS loops_without_display_of_lowest_cost = 0;
	const NUM_ELEMENTS nEdges = static_cast<NUM_ELEMENTS>(m_edges.size());
	const NUM_ELEMENTS iteratation_limit = 100*nEdges;

	double random_number = 0.0;
	double temperature = g_options->get_delay_structure_init_temperature();

	debug("Initial cost with node violations is " << cost);
	debug("Initial temperature " << temperature << endl);

	debug("Initial Solution Quality for Combinational Delay Structure Creation");
	debug("----------------------------------------------------------------------");

	report_solution_quality();

	create_edge_length_state();

	while (cost > 0 && loops < iteratation_limit)
	{
		generate_move(move);

		if (move.is_valid())
		{
			changed_cost = get_changed_cost(move, cost);

			random_number = g_rand_number_gen->random_double_number(1);

			if (changed_cost <= 0.0 || random_number < exp(- changed_cost/temperature)) 
			{
				if (move.edge_changed_clusters())
				{
					update_violation_matrix(move);
					update_edge_length_data_structures(move);
				}
				move.make_move();

				cost += changed_cost;

				if (cost < lowest_cost)
				{
					cost = get_cost();
					lowest_cost = cost;
			
					if (! g_options->is_fast_gen())
					{
						create_edge_length_state();
					}
					
					loops_without_display_of_lowest_cost++;
					if (loops_without_display_of_lowest_cost == 10)
					{
						loops_without_display_of_lowest_cost = 0;
						debug("********************** Lowest cost " << lowest_cost << " **********");
					}
				}
			}
		}

		if (loops > 0 && loops % (nEdges/2) == 0)
		{
			if (loops % (40*nEdges) == 0)
			{
				reset_congestion_costs();

				temperature /= 2;
				temperature = MAX(temperature, 1);
				if (temperature > 1)
				{
					debug("New temperature: " << temperature);
				}
			}
			else
			{
				increase_congestion_cost_for_level_nodes();
			}
			

			// if there are still nodes with too many inputs/too few outputs 
			// let's look to increase the penalty factors to 
			// force edges off the inputs or onto the outputs of the 
			// congested or decongested level nodes

			if (! g_options->is_fast_gen())
			{
				set_graph_into_lowest_cost_state();
			}

			debug("********************** Lowest cost " << get_cost() << " **********");

			cost = get_cost();
			lowest_cost = cost;
		}

		loops++;

	}

	reset_congestion_costs();
	if (! g_options->is_fast_gen())
	{
		set_graph_into_lowest_cost_state();
	}

	lowest_cost = get_cost();

	debug("********************** Level Final Lowest cost " << lowest_cost << " **********");
	debug("********************** Level Final Unnormalized cost " << get_unnormalized_cost() << " **********");
	dbg("\n");
}


// generate a move for the algorithm
//
//
// In this move generator we spend part of the time 
// making inter-cluster moves and part of the time making intra-cluster moves.
//
// For inter-cluster moves, if our violation matrix is not zero or if 
// our edge length distributions are not satisfied we spend part of time
// trying to make moves to specifically satisfy them
//
// It is a little convoluted but it seems to work better than a random experiment.
// At least for the experiments I did.
//
//
// POST: move contains a move that might be valid or invalid
void SEQUENTIAL_LEVEL::generate_move
(
	MOVE& move
)
{
	move.clear();

	CLUSTER_NUMBER_TYPE source_cluster = 0,
						sink_cluster = 0;
	double percentage_time_for_inter_cluster = 0.25;
	if (m_nClusters == 1)
	{
		percentage_time_for_inter_cluster = 0;
	}

	double random_number = g_rand_number_gen->random_double_number(1);

	if (random_number < percentage_time_for_inter_cluster)
	{
		// make an inter cluster move
		random_number = g_rand_number_gen->random_double_number(1);

		if (random_number < 0.5)
		{
			if (! m_violation_matrix.is_zero())
			{
				// try to fix the inter cluster connectivity.
				g_rand_number_gen->discrete_pmf_with_negative_entries(m_violation_matrix, 
																	source_cluster, sink_cluster);
				generate_inter_cluster_move(move, m_clusters[source_cluster], m_clusters[sink_cluster]);
			}
			else if (! are_the_inter_cluster_edge_length_distributions_satisfied())
			{
				// try to fix the edge length distribution for each cluster
				generate_edge_length_move(move);
			}
			else
			{
				// make an inter-cluster translation
				generate_source_sink_inter_cluster_pair(source_cluster, sink_cluster);
				generate_inter_cluster_move(move, m_clusters[source_cluster], m_clusters[sink_cluster]);
			}
		}
		else
		{
			// make an inter-cluster translation
			generate_source_sink_inter_cluster_pair(source_cluster, sink_cluster);
			generate_inter_cluster_move(move, m_clusters[source_cluster], m_clusters[sink_cluster]);
		}
	}
	else 
	{
		// make a translation for an intra-cluster edge
		source_cluster = g_rand_number_gen->random_number(m_nClusters-1);
		generate_intra_cluster_move(move, m_clusters[source_cluster]);
	}
}

//
// Generate an intra-cluster move in the cluster
//
// An intra-cluster move consists of "moving" an individual edge in the cluster
// while keeping the length of the edge constant
// In the move source edge is the edge that will give up weight.
// In the move sink edge is the edge that will accept weight.
//
// PRE: cluster is valid
// POST: a move has been generated
//
void SEQUENTIAL_LEVEL::generate_intra_cluster_move
(
	MOVE& move,
	CLUSTER * cluster
)
{
	if (cluster->get_nIntra_cluster_edges() == 0)
	{
		return;
	}

	LEVEL_NODE * sink_level_node = 0;
	LENGTH_TYPE edge_length = 0;
	EDGE * edge = 0;
	bool is_valid_source = false;
	CLUSTER_NUMBER_TYPE cluster_number = cluster->get_cluster_number();

	//debug("Generative move");

	edge = get_edge_that_will_give_up_weight(cluster, cluster);
	
	if (edge == 0)
	{
		return;
	}

	assert(edge);
	edge_length = edge->get_length();

	// we need delay to spare to move an edge of length 1
	is_valid_source = edge->can_reduce_weight_on_edge();

	if (! is_valid_source)
	{
		return;
	}

	edge_length = edge->get_length();
	move.set_source_edge(edge);
	

	DELAY_TYPE destination_edge_source_delay = (short) g_rand_number_gen->random_number(m_max_comb_delay-edge_length);
	assert(destination_edge_source_delay+edge_length >= 0 && 
			destination_edge_source_delay+edge_length <= m_max_comb_delay);

	sink_level_node = cluster->Level[destination_edge_source_delay+edge_length];
	assert(sink_level_node);
	edge = sink_level_node->find_edge_with_source(cluster_number, destination_edge_source_delay);
	assert(edge);
	
	move.set_sink_edge(edge);
}


//
// Generate an inter-cluster move in the cluster
//
// An intra-cluster move consists of "moving" an individual edge in to a
// new location while keeping the length of the edge constant
// and keeping the source and sink clusters of the edge the same.
// ie. The move is more of translation up or down the delay levels.
//
// In the move source edge is the super edge that will give up weight.
// In the move sink edge is the super edge that will accept weight.
//
// PRE: source_cluster and sink_cluster are valid
// POST: a move has been generated
//
void SEQUENTIAL_LEVEL::generate_inter_cluster_move
(
	MOVE& move,
	CLUSTER * source_cluster,
	CLUSTER * sink_cluster
)
{
	assert(source_cluster && sink_cluster);

	if (source_cluster->get_nInter_cluster_output_edges() == 0)
	{
		return;
	}

	LENGTH_TYPE edge_length = 0;
	EDGE * edge = 0;
	bool is_valid_source = false;

	//debug("Generating inter cluster move");

	// get an edge that will give up weight
	edge = get_edge_that_will_give_up_weight(source_cluster, sink_cluster);
	if (edge == 0)
	{
		return;
	}
	assert(edge);
	edge_length = edge->get_length();

	// let's try out the modification that we need delay to spare to move an edge of length 1
	is_valid_source = edge->can_reduce_weight_on_edge();

	if (! is_valid_source)
	{
		// lets exchange the length 1 inter cluster edge with an intra cluster edge
		//debug("Going to try to swap inter cluster edge here!!!!");

		// this procedure might not be necessary and i
		// is it still workwhile because it enchances randomness?
	
		assert(new_gen);

		make_swap_move(edge, sink_cluster);
		is_valid_source = edge->can_reduce_weight_on_edge();
	}


	if (! is_valid_source)
	{
		return;
	} 

	edge_length = edge->get_length();

	move.set_source_edge(edge);

	edge = get_edge_that_will_accept_weight(move, edge_length);
	
	move.set_sink_edge(edge);
}

//
// Get a super edge that will give up weight in the move.
//
// High cost super edges have a higher probability of having an individual edge taken from it.
// High cost edges are edges whose source level node has too many outputs and whose
// sink level node have too many inputs.
//
//
// PRE: source_cluster and sink_cluster are valid
// RETURNS: an edge to move if one could be selected
//          otherwise 0
//
//
EDGE * SEQUENTIAL_LEVEL::get_edge_that_will_give_up_weight
(
	const CLUSTER * source_cluster,
	const CLUSTER * sink_cluster
) const
{
	DELAY_TYPE source_level = 0;
	DELAY_TYPE sink_level = 0;
	LEVEL_NODE * sink_level_node = 0;
	vector<EDGE *> edge_ptr_map;
	NUM_ELEMENTS edge_number = 0;
	EDGE * edge = 0;
	CLUSTER_NUMBER_TYPE source_cluster_number = source_cluster->get_cluster_number();

	COST_TYPE selection_weight = 0,
	          min_selection_weight = static_cast<COST_TYPE>(MAXLONG),
			  max_selection_weight = 0,
			  total_selection_weight = 0;
	CUMULATIVE_MASS_FUNCTION edge_selection_cmf;

	// find the minimum and maximimum selection weight among the edges 
	for (source_level = 0; source_level < m_max_comb_delay; source_level++)
	for (sink_level = source_level+1; sink_level <= m_max_comb_delay; sink_level++)
	{
		sink_level_node = sink_cluster->Level[sink_level];
		assert(sink_level_node);
		edge = sink_level_node->find_edge_with_source(source_cluster_number, source_level);
		assert(edge);

		if (edge->get_weight() > 0)
		{
			selection_weight = edge->get_selection_weight();
			if (selection_weight < min_selection_weight)
			{
				min_selection_weight = selection_weight;
			}
			if (selection_weight > max_selection_weight)
			{
				max_selection_weight = selection_weight;
			}
		}
	}

	// if our min selection_weight is less than 0, add a bias to all the edge costs.
	if (min_selection_weight < 0.0)
	{
		min_selection_weight = -(min_selection_weight);
	}
	else
	{
		min_selection_weight = 0;
	}

	// add a small prob. to make sure all edges have a chance at being selected
	min_selection_weight += (min_selection_weight + max_selection_weight)*0.05;
	
	for (source_level = 0; source_level < m_max_comb_delay; source_level++)
	for (sink_level = source_level+1; sink_level <= m_max_comb_delay; sink_level++)
	{
		sink_level_node = sink_cluster->Level[sink_level];
		assert(sink_level_node);

		edge = sink_level_node->find_edge_with_source(source_cluster_number, source_level);
		assert(edge);

		// if the edge has weight to give
		if (edge->get_weight() > 0)
		{
			selection_weight = edge->get_selection_weight();
			selection_weight = selection_weight + min_selection_weight;

			total_selection_weight += selection_weight;

			edge_selection_cmf[total_selection_weight] = edge_number;

			edge_number++;

			edge_ptr_map.push_back(edge);
		}
	}


	// if no edge had any weight to to give or any selection weight
	// then return that no edge could be chosen
	if (total_selection_weight == 0)
	{
		return 0;
	}

	// randomly pick an edge
	edge_number = g_rand_number_gen->discrete_cmf(edge_number, total_selection_weight, edge_selection_cmf);
	edge = edge_ptr_map[edge_number];

	assert(edge && edge->get_weight() >= 0);

	return edge;
}

// Get a super edge that will accept weight in the move.
//
// PRE: edge_length is the length we want our super edge to have
// RETURNS: an edge to move if one could be selected
//          otherwise 0
//
//
EDGE * SEQUENTIAL_LEVEL::get_edge_that_will_accept_weight
(
	const MOVE& move,
	const LENGTH_TYPE& edge_length
) const
{
	DELAY_TYPE source_delay;
	CLUSTER_NUMBER_TYPE source_cluster_id = 0,
	                    sink_cluster_id = 0;
	EDGE * edge = 0;
	CLUSTER * sink_cluster = 0;
	LEVEL_NODE * sink_level_node = 0;

	generate_source_sink_pair_for_edge_that_will_accept_weight(move, source_cluster_id, sink_cluster_id);

	source_delay = (short) g_rand_number_gen->random_number(m_max_comb_delay-edge_length);
	assert(source_delay + edge_length >= 0 && source_delay+edge_length <= m_max_comb_delay);

	sink_cluster = m_clusters[sink_cluster_id];
	assert(sink_cluster);

	sink_level_node = sink_cluster->Level[source_delay+edge_length];
	assert(sink_level_node);

	edge = sink_level_node->find_edge_with_source(source_cluster_id, source_delay);
	assert(edge && edge->get_length() == edge_length);

	return edge;
}

// Get the source and sink clusters for the super edge that we will
// accept weight in the move
//
// Needs a rewrite
//
//
// PRE: nothing
// POST: source_cluster_id and sink_cluster_id have been calculated
//
void SEQUENTIAL_LEVEL::generate_source_sink_pair_for_edge_that_will_accept_weight
(
 	const MOVE& move,
	CLUSTER_NUMBER_TYPE& source_cluster_id,
	CLUSTER_NUMBER_TYPE& sink_cluster_id
 ) const
{
	double random_number = g_rand_number_gen->random_double_number(1);
	double percentage_time_for_cluster_to_cluster_moves = 0.20;
	EDGE * source_edge = 0;

	// if our inter-cluster connectivity is not quite right 
	if (m_violation_matrix.is_positive())
	{
		// move the weight of an invidual edge to a different source-sink cluster pair
		if ( random_number < percentage_time_for_cluster_to_cluster_moves)
		{
			generate_source_sink_inter_cluster_pair(source_cluster_id, sink_cluster_id);
		}
		else
		{
			// just make a inter cluster translation within the level nodes
			source_edge = move.get_source_edge();
			assert(source_edge);
			source_cluster_id = source_edge->get_source_cluster_number();
			sink_cluster_id = source_edge->get_sink_cluster_number();
		}
	}
	else
	{
		// try to make a move to satisfy the inter cluster connection matrix
		g_rand_number_gen->discrete_inverse_pmf_with_negative_entries(m_violation_matrix, 
																	source_cluster_id, sink_cluster_id);
	}
}



// Length 1 edges are special because they can define the delay level of a node.
// Because they are special I do not allow a move that will leave less
// length 1 edges than nodes.
//
// (I should maybe add delay to the cost function.
// It might give me much better solutions).
//
//	We can't move an inter-cluster edges because the sink node it points to 
//	doesn't have enough length 1 edges.
//
//	Take an length 1 edge from somewhere in the cluster and 
//	add it to the inter cluster edge's sink level node so that we can move the 
//	inter-cluster edge.
//
// PRE: We want to move a length 1 inter-cluster edge from a level node.
//      But we can't because the level node has just the right number of length 1 edges. 
//      One length 1 edge less would mean we would not have enough edges to 
//      define the delay level of the individual edges.
//
// POST: IF we could find a length 1 intra-cluster edge to move to the inter-cluster edge's sink node
//       we have moved it. We can now move the inter-cluster edge
//
//       If we couldn't. Then do nothing.
//
//
//
void SEQUENTIAL_LEVEL::make_swap_move
(
	EDGE * inter_cluster_edge,
	CLUSTER * sink_cluster
)
{
	assert(inter_cluster_edge && sink_cluster);

	MOVE move;
	LEVEL_NODE * inter_cluster_source_node = inter_cluster_edge->get_source_level_node();
	LEVEL_NODE * inter_cluster_sink_node = inter_cluster_edge->get_sink_level_node();
	LEVEL_NODE * intra_cluster_sink_node = 0;
	EDGE * intra_cluster_edge = 0;
	EDGE * edge = 0;
	CLUSTER_NUMBER_TYPE sink_cluster_id = sink_cluster->get_cluster_number();
	DELAY_TYPE source_delay = inter_cluster_sink_node->get_delay_level() - 1;

	assert(inter_cluster_sink_node && inter_cluster_source_node);
	
	// meant for inter cluster length 1 edges where the sink has no delay to spare
	assert(inter_cluster_sink_node->get_delay_difference() == 0 && inter_cluster_edge->get_length() == 1);

	
	// look for an intra cluster length 1 edge to add weight to the sink node 
	// so that we can move the inter cluster edge

	intra_cluster_edge = find_length_1_edge_to_translate(sink_cluster);

	if (intra_cluster_edge == 0)
	{
		// could not find an intra cluster edge to move.  
		// just return
		return;
	}
	intra_cluster_sink_node = intra_cluster_edge->get_sink_level_node();
	assert(intra_cluster_sink_node);

	// we have an acceptable edge
	assert(intra_cluster_edge->is_intra_cluster() && 
			intra_cluster_edge->get_length() == 1 && 
			intra_cluster_sink_node->get_delay_difference() > 0);


	edge = inter_cluster_sink_node->find_edge_with_source(sink_cluster_id, source_delay);
	assert(edge);


	move.set_source_edge(intra_cluster_edge);
	move.set_sink_edge(edge);

	assert(new_gen);
	// adding this requirement, which might just defeat
	// the above code because we are getting 
	// weight added to level nodes with no simple nodes
	// look at it again later
	//
	if (move.is_valid())
	{

		//COST_TYPE swap_costs = move.get_changed_cost();
		//wanted_move.set_swap_costs(swap_costs);

		//debug("Making a swap move.");
		move.make_move();
	}
}


// look for an intra cluster length 1 edge that can be moved
// within the cluster
//
// PRE: cluster is valid
// RETURNS: a length 1 intra-cluster edge that can be moved if one could
//          be found, 0 otherwise.
EDGE * SEQUENTIAL_LEVEL::find_length_1_edge_to_translate
(
	CLUSTER * cluster
) const
{
	const LENGTH_TYPE edge_length = 1;
	assert(edge_length <= m_max_comb_delay);

	DELAY_TYPE source_delay = (short) g_rand_number_gen->random_number(m_max_comb_delay-edge_length);
	DELAY_TYPE sink_delay = source_delay + edge_length;
	DELAY_TYPE old_source_delay = 0;
	LEVEL_NODE * source_level_node = 0;
	LEVEL_NODE * sink_level_node = 0;
	EDGE * edge = 0;
	bool have_gone_through_all_edges;
	
	assert(	source_delay >= 0 && source_delay <= m_max_comb_delay &&
			sink_delay >= 0 && sink_delay <= m_max_comb_delay);


	edge = cluster->get_edge(source_delay, edge_length);

	if (! edge->has_spare_weight())
	{
		// the edge and the level nodes did not had some weight to spare
		// try to increase delay until we find an edge and level node that does
		// if no edge can be found, just try to find an edge that we can remove some weight from
		
		old_source_delay = source_delay;
		have_gone_through_all_edges = false;

		while (! edge->has_spare_weight() && ! have_gone_through_all_edges)
		{
			source_delay = (source_delay +1) % (m_max_comb_delay-edge_length+1);

			edge = cluster->get_edge(source_delay, edge_length);
			assert(edge);

			have_gone_through_all_edges = (source_delay == old_source_delay);
		}

		if (have_gone_through_all_edges)
		{
		
			have_gone_through_all_edges = false;

			while (! edge->can_reduce_weight_on_edge() && ! have_gone_through_all_edges)
			{
				source_delay = (source_delay +1) % (m_max_comb_delay-edge_length+1);
				edge = cluster->get_edge(source_delay, edge_length);
				assert(edge);

				have_gone_through_all_edges = (source_delay == old_source_delay);
			}

			if (have_gone_through_all_edges)
			{
				// out of luck
				return 0;
			}
		}
		else
		{
			//debug("delay level " << source_delay << " had edges to spare.");
		}

	}
	else
	{
		//debug("We have spare at delay level " << source_delay);
	}

	// just checking to make sure everything is right
	sink_delay = source_delay+edge_length;

	assert(	source_delay >= 0 && source_delay <= m_max_comb_delay &&
			sink_delay >= 0 && sink_delay <= m_max_comb_delay);

	source_level_node = edge->get_source_level_node();
	sink_level_node = edge->get_sink_level_node();
	assert(source_level_node && sink_level_node);


	if (edge_length == 1)
	{
		assert(sink_level_node->get_delay_difference() > 0);
	}
	assert(source_level_node->get_intra_cluster_output_edge_lengths()[edge_length] > 0);
	assert(sink_level_node->get_intra_cluster_input_edge_lengths()[edge_length] > 0);
	assert(edge->get_weight() > 0);

	return edge;
}




//
// examine all the edges to see if one them has weight that can 
// be reduced without violating the lower bounds on the level nodes 
// AND/OR an edge that can be reduced without violating the delay
// 
// PRE:  nothing
// POST: deletion_situation has been updated 
void SEQUENTIAL_LEVEL::update_deletion_situation
(
	SEQUENTIAL_LEVEL::DELETION_SITUATION & deletion_situation
) const
{

	if (deletion_situation == SEQUENTIAL_LEVEL::NO_EDGE_CAN_BE_DELETED)
	{
		return;
	}

	// start with the assumption that no edge can be deleted.
	deletion_situation = SEQUENTIAL_LEVEL::NO_EDGE_CAN_BE_DELETED;

	EDGES::const_iterator edges_iter = m_edges.begin();
	bool edge_has_spare_weight = false; 
	bool can_reduce_weight_on_edge = false;
	EDGE * edge = 0;
	
	while (edges_iter != m_edges.end() && deletion_situation != SEQUENTIAL_LEVEL::LEVEL_NODES_HAVE_SPARE_WEIGHT)
	{
		edge = *edges_iter;
		assert(edge);
		edge_has_spare_weight = edge->has_spare_weight();
		can_reduce_weight_on_edge = edge->can_reduce_weight_on_edge();

		if (edge_has_spare_weight)
		{
			deletion_situation = SEQUENTIAL_LEVEL::LEVEL_NODES_HAVE_SPARE_WEIGHT;
		}
		else if (can_reduce_weight_on_edge)
		{
			deletion_situation = SEQUENTIAL_LEVEL::LEVEL_NODES_HAVE_WEIGHT;
		}

		edges_iter++;
	}
}


//
// PRE: there are at least 2 clusters
// RETURNS: two source and sink cluster numbers that are not equal
//
void SEQUENTIAL_LEVEL::generate_source_sink_inter_cluster_pair
(
	CLUSTER_NUMBER_TYPE& source_cluster_number, 
	CLUSTER_NUMBER_TYPE& sink_cluster_number
) const
{
	source_cluster_number = g_rand_number_gen->random_number(m_nClusters-1);
	sink_cluster_number = g_rand_number_gen->random_number(m_nClusters-1);

	while (sink_cluster_number == source_cluster_number)
	{
		sink_cluster_number = g_rand_number_gen->random_number(m_nClusters-1);
	}
}

//
// POST: For each cluster:
// 			outputs if it has too many or too little edges of a certain length
//          outputs what is the difference between the current shapes and their specs
//
//		 We are sure that sure that the total number of inter-cluster input edge length violations 
//		 in the circuit is equal to the total number of inter-cluster output edge length violations
//
void SEQUENTIAL_LEVEL::output_all_violations()
{
	CLUSTERS::iterator cluster_iter;
	CLUSTER * cluster = 0;
	NUM_ELEMENTS cluster_number = 0;
	DISTRIBUTION my_input_violations; 
	DISTRIBUTION my_output_violations;
	NUM_ELEMENTS edge_length_size = m_max_comb_delay + 1;
	DISTRIBUTION input_violations_by_edge_length(edge_length_size, 0);
	DISTRIBUTION output_violations_by_edge_length(edge_length_size, 0);
	NUM_ELEMENTS edge_length;

	for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end(); cluster_iter++)
	{
		debug("********* Cluster " << cluster_number << " **********\n");
		cluster = *cluster_iter;
		assert(cluster);

		output_difference_between_current_and_spec(cluster);

		my_input_violations = get_inter_cluster_input_edge_length_violations(cluster);
		my_output_violations = get_inter_cluster_output_edge_length_violations(cluster);

		for (edge_length = 0; edge_length < edge_length_size; edge_length++)
		{
			input_violations_by_edge_length[edge_length] += my_input_violations[edge_length];
			output_violations_by_edge_length[edge_length] += my_output_violations[edge_length];
		}

		cluster_number++;
	}	

	// check to make sure the total number of inter-cluster input edge length violations in
	// the circuit is equal to the total number of inter-cluster output edge length violations
	for (edge_length = 0; edge_length < edge_length_size; edge_length++)
	{
		assert(input_violations_by_edge_length[edge_length] == output_violations_by_edge_length[edge_length]);
	}
}
//
// PRE:  cluster is valid
// POST: For the cluster:
// 			outputs if it has too many or too little edges of a certain length
//          outputs what is the difference between the current shapes and their specs
//
void SEQUENTIAL_LEVEL::output_difference_between_current_and_spec
(
	CLUSTER * cluster
)
{
	DISTRIBUTION input_violations = get_inter_cluster_input_edge_length_violations(cluster);
	DISTRIBUTION output_violations = get_inter_cluster_output_edge_length_violations(cluster);
	DISTRIBUTION input_difference = cluster->get_input_shape_difference_from_spec();
	DISTRIBUTION output_difference = cluster->get_output_shape_difference_from_spec();
	DISTRIBUTION delay_difference = cluster->get_delay_difference_by_delay_level();

	cout << "Difference between current and spec for: \n\n";

	cout << "  Inter-cluster input edge lengths\n  ";
	copy(input_violations.begin(),input_violations.end(), 
			ostream_iterator<NUM_ELEMENTS>(cout, " ")); cout << endl;

	cout << "Inter-cluster output edge lengths\n  ";
	copy(output_violations.begin(), output_violations.end(), 
			ostream_iterator<NUM_ELEMENTS>(cout, " ")); cout << endl;

	cout << "Input Shape\n  ";
	copy(input_difference.begin(), input_difference.end(), 
			ostream_iterator<NUM_ELEMENTS>(cout, " ")); cout << endl;

	cout << "Output Shape\n  ";
	copy(output_difference.begin(), output_difference.end(), 
			ostream_iterator<NUM_ELEMENTS>(cout, " ")); cout << endl;

	cout << "\nDifference between the number of delay defining edges and the number needed\n  ";
	copy(delay_difference.begin(), delay_difference.end(), 
			ostream_iterator<NUM_ELEMENTS>(cout, " ")); cout << endl;

	cout << endl;
}




//
// POST: we have outputed the inter-cluster connectivity matrices by edge length
//
void SEQUENTIAL_LEVEL::output_inter_cluster_matricies_by_edge_length() const
{
	assert(m_circuit);

	LENGTH_TYPE edge_length_size = m_circuit->get_delay() + 1,
	            edge_length = 0;
	NUM_ELEMENTS weight = 0;
	
	vector<MATRIX> connectivity_by_edge_length(edge_length_size);
	EDGES::const_iterator edge_iter;
	EDGE * edge = 0;
	CLUSTER_NUMBER_TYPE source_cluster_number = 0,
	                    sink_cluster_number = 0;

	for (edge_length = 0; edge_length < edge_length_size; edge_length++)
	{
		connectivity_by_edge_length[edge_length].resize(m_nClusters, m_nClusters);
	}

	for (edge_iter = m_inter_cluster_edges.begin(); edge_iter != m_inter_cluster_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge && edge->is_inter_cluster());

		edge_length = edge->get_length();
		assert(edge_length >= 0 && edge_length < edge_length_size);

		weight = edge->get_weight();
		assert(weight >= 0);

		source_cluster_number = edge->get_source_cluster_number();
		sink_cluster_number = edge->get_sink_cluster_number();
		assert(source_cluster_number >= 0 && source_cluster_number < m_nClusters);
		assert(sink_cluster_number >= 0 && sink_cluster_number < m_nClusters);


		connectivity_by_edge_length[edge_length](source_cluster_number, sink_cluster_number) +=
			weight;
	}

	for (edge_length = 1; edge_length < edge_length_size; edge_length++)
	{
		debug("Connectivity matrix. edge length = " << edge_length << endl 
				<< connectivity_by_edge_length[edge_length] << endl);
	}

}


//
// PRE: move is valid
// POST: we have updated the violation matrix
//
void SEQUENTIAL_LEVEL::update_violation_matrix
(
	MOVE & move
)
{
	//debug("Violation matrix before move \n" << m_violation_matrix);
	EDGE * source_edge = move.get_source_edge();
	EDGE * sink_edge  = move.get_sink_edge();
	assert(source_edge && sink_edge);
	CLUSTER_NUMBER_TYPE source_cluster_number 	= source_edge->get_source_cluster_number();
	CLUSTER_NUMBER_TYPE sink_cluster_number 	= source_edge->get_sink_cluster_number();

	m_violation_matrix(source_cluster_number, sink_cluster_number) -= 1;

	source_cluster_number = sink_edge->get_source_cluster_number();
	sink_cluster_number = sink_edge->get_sink_cluster_number();

	m_violation_matrix(source_cluster_number, sink_cluster_number) += 1;

	//debug("Violation matrix after move \n" << m_violation_matrix);
}
//
// PRE: move is valid
// POST: we have updated the edge length violation structures
//
void SEQUENTIAL_LEVEL::update_edge_length_data_structures
(
	MOVE & move
)
{
	//debug("Violation matrix before move \n" << m_violation_matrix);
	EDGE * source_edge = move.get_source_edge();
	EDGE * sink_edge  = move.get_sink_edge();
	assert(source_edge && sink_edge);
	assert(new_gen);
	assert(source_edge->is_inter_cluster() && sink_edge->is_inter_cluster());
	LENGTH_TYPE edge_length = source_edge->get_length();
	assert(sink_edge->get_length() == edge_length);

	CLUSTER_NUMBER_TYPE source_source_cluster_number= source_edge->get_source_cluster_number();
	CLUSTER_NUMBER_TYPE source_sink_cluster_number 	= source_edge->get_sink_cluster_number();
	CLUSTER_NUMBER_TYPE dest_source_cluster_number 	= sink_edge->get_source_cluster_number();
	CLUSTER_NUMBER_TYPE dest_sink_cluster_number 	= sink_edge->get_sink_cluster_number();

	m_edges_between_cluster(source_source_cluster_number, source_sink_cluster_number, edge_length) -= 1;
	m_edges_between_cluster(dest_source_cluster_number, dest_sink_cluster_number, edge_length) += 1;
	//debug("Violation matrix after move \n" << m_violation_matrix);

	update_cluster_for_edge_length_data_structures(source_source_cluster_number, source_sink_cluster_number,
			edge_length);
	update_cluster_for_edge_length_data_structures(dest_source_cluster_number, dest_sink_cluster_number,
			edge_length);
}


//
// PRE: m_violation_matrix is not all zero which may imply that 
//      we have too many inter-cluster edges in the graph
// POST: we have removed the excess number of edges from the 
//		 graph for all edge length except maybe for edges of length 1
//
//       If we have length 1 edges we might have 
//       removed less if we needed some of the edges to define delay
void SEQUENTIAL_LEVEL::reduce_edges_to_meet_total_number_of_edges()
{
	DISTRIBUTION total_number_of_edges(m_max_comb_delay+1, 0);
	DISTRIBUTION inter_cluster_spec(m_max_comb_delay+1, 0);
	DISTRIBUTION inter_cluster_input, inter_cluster_output;
	EDGES::iterator edge_iter;
	EDGE * edge = 0;
	LENGTH_TYPE edge_length;
	NUM_ELEMENTS edge_weight;
	CLUSTERS::iterator cluster_iter;
	CLUSTER * cluster = 0;
	CLUSTER_SPEC * comb_spec = 0;
	NUM_ELEMENTS number_to_remove = 0;

	// find the total number of edges of each edge length
	for (edge_iter = m_inter_cluster_edges.begin(); edge_iter != m_inter_cluster_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		edge_length = edge->get_length();
		edge_weight = edge->get_weight();
		assert(edge_weight >= 0 && edge_length >= 0 && edge_length <= m_max_comb_delay);

		total_number_of_edges[edge_length] += edge_weight;
	}


	// find the total edges of each edge length that should be there according to spec
	for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end();
			cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);
		comb_spec = cluster->get_comb_spec();
		assert(comb_spec);

		inter_cluster_input = comb_spec->get_inter_cluster_input_edge_lengths();
		inter_cluster_output = comb_spec->get_inter_cluster_output_edge_lengths();

		for (edge_length = 1; edge_length <= m_max_comb_delay; edge_length++)
		{
			inter_cluster_spec[edge_length] += inter_cluster_input[edge_length];
			inter_cluster_spec[edge_length] += inter_cluster_output[edge_length];
		}
	}


	// because each edge produces an inter-cluster input edge length and an
	// inter-cluster output edge length divided the number of edges at each edge length by 2
	for (edge_length = 1; edge_length <= m_max_comb_delay; edge_length++)
	{
		assert(inter_cluster_spec[edge_length] % 2 == 0);
		inter_cluster_spec[edge_length] = inter_cluster_spec[edge_length] / 2;
	}

	// remove edges to come down to spec
	for (edge_length = m_max_comb_delay; edge_length >= 1; edge_length--)
	{
		number_to_remove = total_number_of_edges[edge_length] - inter_cluster_spec[edge_length];

		if (number_to_remove > 0)
		{
			find_and_remove_edge_lengths(number_to_remove, edge_length);
		}
		else if (number_to_remove < 0)
		{
			assert(new_gen);
			// should add edges in this case?
			// input-ouput shape/edge length conflict
		}
		else
		{
			// things are good
		}
	}
}


// PRE: number_to_remove is the number of edges to remove
//      edge_length is the length of edge to target
// POST: we have removed number_to_remove weight from the 
//       edges 
//
//       If we have length 1 edges we might have 
//       removed less if we needed some of the edges to define delay
//
void SEQUENTIAL_LEVEL::find_and_remove_edge_lengths
(
	NUM_ELEMENTS & number_to_remove,
	const LENGTH_TYPE& edge_length
)
{
	assert(number_to_remove >= 0);
	assert(edge_length >= 0 && edge_length < m_max_comb_delay +1);
	EDGES::iterator edge_iter = m_edges.begin();
	NUM_ELEMENTS edge_weight = 0;
	NUM_ELEMENTS weight_to_remove = 0;
	CLUSTER_NUMBER_TYPE source_cluster = 0;
	CLUSTER_NUMBER_TYPE sink_cluster = 0;
	EDGE * edge = 0;


	while (number_to_remove > 0 && edge_iter != m_edges.end())
	{
		edge = *edge_iter;
		assert(edge);

		if (edge->get_length() == edge_length && edge->is_inter_cluster())
		{
			edge_weight = edge->get_weight();

			weight_to_remove  = MIN(edge_weight, number_to_remove);

			if (edge_length == 1)
			{
				weight_to_remove  = MIN(weight_to_remove, edge->get_sink_level_node()->get_delay_to_spare());
			}

			if (weight_to_remove > 0)
			{
				edge->add_weight(-weight_to_remove);
				source_cluster = edge->get_source_cluster_number();
				sink_cluster = edge->get_sink_cluster_number();
				m_violation_matrix(source_cluster, sink_cluster) -= weight_to_remove;

				number_to_remove -= weight_to_remove;
			}
		}

		edge_iter++;
	}

	assert(number_to_remove == 0);
}

// RETURNS: are the inter-cluster edge length distributions satisfied?
//
bool SEQUENTIAL_LEVEL::are_the_inter_cluster_edge_length_distributions_satisfied()
{
	CLUSTER* cluster = 0;
	CLUSTERS::iterator cluster_iter = m_clusters.begin(); 
	bool edge_lengths_satisfied = true;

	while (edge_lengths_satisfied && cluster_iter != m_clusters.end())
	{
		cluster = *cluster_iter;
		assert(cluster);

		edge_lengths_satisfied = cluster->is_inter_cluster_edge_length_distributions_satisfied();

		cluster_iter++;
	}

	return edge_lengths_satisfied;
}

void SEQUENTIAL_LEVEL::create_edge_length_matrix()
{
	CLUSTER* source_cluster = 0;
	CLUSTER* sink_cluster = 0;
	CLUSTER_SPEC * source_comb_spec = 0;
	CLUSTER_SPEC * sink_comb_spec = 0;
	DISTRIBUTION violations_input_edge_length;
	DISTRIBUTION violations_output_edge_length;
	NUM_ELEMENTS edge_length_size = m_max_comb_delay + 1;
	NUM_ELEMENTS source_cluster_number, sink_cluster_number;
	NUM_ELEMENTS edge_length = 0;
	m_edges_between_cluster.clear();
	m_possible_edges_between_cluster.clear();
	m_source_edge_length_matrix.clear();
	m_dest_edge_length_matrix.clear();
	m_edges_between_cluster.resize(m_nClusters, m_nClusters, edge_length_size);
	m_possible_edges_between_cluster.resize(m_nClusters, m_nClusters, edge_length_size);
	m_source_edge_length_matrix.resize(m_nClusters, m_nClusters, edge_length_size);
	m_dest_edge_length_matrix.resize(m_nClusters, m_nClusters, edge_length_size);
	NUM_ELEMENTS source_number = 0;
	NUM_ELEMENTS destination_number = 0;
	DISTRIBUTION number_edges_between_two_clusters;
	DISTRIBUTION number_poss_edges_between_two_clusters;

	//output_all_violations();

	// create a matrix: Cluster x Edge length violations for that cluster
	for (source_cluster_number = 0; source_cluster_number < m_nClusters; source_cluster_number++)
	for (sink_cluster_number = 0; sink_cluster_number < m_nClusters; sink_cluster_number++)
	{
		if (source_cluster_number != sink_cluster_number)
		{
			number_edges_between_two_clusters.clear();
			number_poss_edges_between_two_clusters.clear();
			number_edges_between_two_clusters.resize(edge_length_size, 0);
			number_poss_edges_between_two_clusters.resize(edge_length_size, 0);

			//debug("source/sink cluster " << source_cluster_number << " / " << sink_cluster_number);
			source_cluster = m_clusters[source_cluster_number];
			sink_cluster = m_clusters[sink_cluster_number];
			assert(source_cluster && sink_cluster);

			source_comb_spec = source_cluster->get_comb_spec();
			sink_comb_spec = sink_cluster->get_comb_spec();
			assert(source_comb_spec && sink_comb_spec);
			violations_input_edge_length = get_inter_cluster_input_edge_length_violations(sink_cluster);
			violations_output_edge_length = get_inter_cluster_output_edge_length_violations(source_cluster);
			get_number_edges_between_two_clusters(source_cluster_number, sink_cluster_number,
												number_edges_between_two_clusters,
												number_poss_edges_between_two_clusters);



			for (edge_length = 1; edge_length < edge_length_size; edge_length++)
			{
				m_edges_between_cluster(source_cluster_number, sink_cluster_number, edge_length) = 
									number_edges_between_two_clusters[edge_length];

				m_possible_edges_between_cluster(source_cluster_number, sink_cluster_number, edge_length) = 
									number_poss_edges_between_two_clusters[edge_length];

				if (number_edges_between_two_clusters[edge_length] < 1 ||
					violations_output_edge_length[edge_length] < 0 ||
					violations_input_edge_length[edge_length] < 0)

				{
					if (number_edges_between_two_clusters[edge_length] > 1)
					{
						source_number = 1;
					}
					else
					{
						//debug("The number between clusters at edge length " << edge_length << " is " << number_edges_between_two_clusters[edge_length]);
						source_number = 0;
					}
				}
				else
				{
					source_number = 4 * MAX(violations_output_edge_length[edge_length],
										violations_input_edge_length[edge_length]);

					if (source_number < 4)
					{
						source_number = 4;
					}
				}
				//debug("final source number = " << source_number);


				m_source_edge_length_matrix(source_cluster_number, sink_cluster_number, edge_length) = 
					source_number;
			
				if (number_poss_edges_between_two_clusters[edge_length] < 1 ||
					violations_output_edge_length[edge_length] > 0 ||
					violations_input_edge_length[edge_length] > 0)

				{
					if (number_poss_edges_between_two_clusters[edge_length] > 1)
					{
						destination_number = 1;
					}
					else
					{
						destination_number = 0;
					}
				}
				else
				{
					destination_number = 4 * MAX(-1 * violations_output_edge_length[edge_length],
										-1 * violations_input_edge_length[edge_length]);

					if (destination_number < 4)
					{
						destination_number = 4;
					}
				}
				assert(destination_number >= 0);

				m_dest_edge_length_matrix(source_cluster_number, sink_cluster_number, edge_length) = 
					destination_number;
			}
		}
	}

	sanity_check_edge_length_matricies();
}

// Is this a good sanity check?
//
// POST: We have no negative entries in:
//			 m_source_edge_length_matrix
//			 m_dest_edge_length_matrix
//			 m_edges_between_cluster
//			 m_possible_edges_between_cluster
//		 The diagonals of the above matrices are zero.
//
void SEQUENTIAL_LEVEL::sanity_check_edge_length_matricies()
{
	NUM_ELEMENTS edge_length_size = m_max_comb_delay + 1;
	NUM_ELEMENTS source_cluster_number; 
	NUM_ELEMENTS edge_length = 0;

	assert(m_source_edge_length_matrix.get_absolute_sum() >= 0);
	assert(m_dest_edge_length_matrix.get_absolute_sum() >= 0);
	assert(m_edges_between_cluster.get_absolute_sum() >= 0);
	assert(m_possible_edges_between_cluster.get_absolute_sum() >= 0);

	for (source_cluster_number = 0; source_cluster_number < m_nClusters; source_cluster_number++)
	{
		for (edge_length = 1; edge_length < edge_length_size; edge_length++)
		{
			assert(m_source_edge_length_matrix(source_cluster_number, source_cluster_number, edge_length)
					== 0);
			assert(m_dest_edge_length_matrix(source_cluster_number, source_cluster_number, edge_length)
					== 0);
			assert(m_edges_between_cluster(source_cluster_number, source_cluster_number, edge_length)
					== 0);
			assert(m_possible_edges_between_cluster(source_cluster_number, source_cluster_number, edge_length)
					== 0);
		}
	}
}

// PRE: source_cluster_number,  sink_cluster_number, edge_length are valid
// POST: m_source_edge_length_matrix and m_dest_edge_length_matrix have been updated
//
void SEQUENTIAL_LEVEL::update_cluster_for_edge_length_data_structures
(
	const CLUSTER_NUMBER_TYPE & source_cluster_number,
	const CLUSTER_NUMBER_TYPE & sink_cluster_number,
	const LENGTH_TYPE& edge_length
)
{
	assert(edge_length >= 0 && edge_length < m_max_comb_delay+1);
	CLUSTER * source_cluster = m_clusters[source_cluster_number];
	CLUSTER * sink_cluster = m_clusters[sink_cluster_number];
	assert(source_cluster && sink_cluster);
	DISTRIBUTION violations_input_edge_length = get_inter_cluster_input_edge_length_violations(sink_cluster);
	DISTRIBUTION violations_output_edge_length = get_inter_cluster_output_edge_length_violations(source_cluster);
	NUM_ELEMENTS number_edges_between_two_clusters = 0;
	NUM_ELEMENTS number_poss_edges_between_two_clusters = 0;
	NUM_ELEMENTS source_number = 0;
	NUM_ELEMENTS destination_number = 0;

	number_edges_between_two_clusters = 
			m_edges_between_cluster(source_cluster_number, sink_cluster_number, edge_length);

	if (number_edges_between_two_clusters < 1 || violations_output_edge_length[edge_length] < 0 || 
													violations_input_edge_length[edge_length] < 0)
	{
		if (number_edges_between_two_clusters > 1)
		{
			source_number = 1;
		}
		else
		{
			source_number = 0;
		}
	}
	else
	{
		source_number = 4 * MAX(violations_output_edge_length[edge_length],
							violations_input_edge_length[edge_length]);

		if (source_number < 4)
		{
			source_number = 4;
		}
	}
	//debug("final source number = " << source_number);


	m_source_edge_length_matrix(source_cluster_number, sink_cluster_number, edge_length) = source_number;
	
	number_poss_edges_between_two_clusters = 
		m_possible_edges_between_cluster(source_cluster_number, sink_cluster_number, edge_length);

	if (number_poss_edges_between_two_clusters < 1 ||
		violations_output_edge_length[edge_length] > 0 ||
		violations_input_edge_length[edge_length] > 0)

	{
		if (number_poss_edges_between_two_clusters > 1)
		{
			destination_number = 1;
		}
		else
		{
			destination_number = 0;
		}
	}
	else
	{
		destination_number = 4 * MAX(-1 * violations_output_edge_length[edge_length],
							-1 * violations_input_edge_length[edge_length]);

		if (destination_number < 4)
		{
			destination_number = 4;
		}
	}
	assert(destination_number >= 0);

	m_dest_edge_length_matrix(source_cluster_number, sink_cluster_number, edge_length) = 
		destination_number;
}

// generate a move that will attempt to resolve differences between
// the current edge length distributions and their spec
//
// PRE: nothing
// POST: move has an edge length move
//
void SEQUENTIAL_LEVEL::generate_edge_length_move
(
	MOVE& move
)
{
	NUM_ELEMENTS source_cluster_number, sink_cluster_number;
	NUM_ELEMENTS dest_source_cluster_number, dest_sink_cluster_number;
	NUM_ELEMENTS source_edge_length = 0,
				 sink_edge_length = 0;
	LENGTH_TYPE edge_length = 0,
				edge_length_size = m_max_comb_delay + 1;

	LEVEL_NODE * sink_level_node = 0;
	EDGE * edge = 0;
	DELAY_TYPE source_delay = 0;

	// pick a source and sink cluster with an edge length that needs moving
	g_rand_number_gen->discrete_pmf(m_source_edge_length_matrix,
									source_cluster_number, sink_cluster_number, source_edge_length);
	g_rand_number_gen->discrete_pmf(m_dest_edge_length_matrix,
								dest_source_cluster_number, dest_sink_cluster_number, sink_edge_length);
	while (source_edge_length != sink_edge_length)
	{
		g_rand_number_gen->discrete_pmf(m_source_edge_length_matrix,
										source_cluster_number, sink_cluster_number, source_edge_length);
		g_rand_number_gen->discrete_pmf(m_dest_edge_length_matrix,
									dest_source_cluster_number, dest_sink_cluster_number, sink_edge_length);
	}
	assert(source_cluster_number != sink_cluster_number);
	assert(dest_source_cluster_number != dest_sink_cluster_number);
	assert(source_edge_length == sink_edge_length);

	edge_length = static_cast<short>(source_edge_length);
	assert(edge_length >= 0 && edge_length < edge_length_size);

	assert(m_clusters[source_cluster_number]->get_nInter_cluster_output_edges(edge_length) > 0);
	assert(m_clusters[sink_cluster_number]->get_nInter_cluster_input_edges(edge_length) > 0);

	//
	//	We now have a source and sink cluster and an edge length to move 
	//	now we just need to exact edge to move
	//

	// find an actual edge that can fit to the source cluster with that edge length
	source_delay = (short) g_rand_number_gen->random_number(m_max_comb_delay-edge_length);
	assert(source_delay+edge_length >= 0 && source_delay+edge_length <= m_max_comb_delay);
	sink_level_node = m_clusters[sink_cluster_number]->Level[edge_length+source_delay];
	assert(sink_level_node);
	edge = sink_level_node->find_edge_with_source(source_cluster_number, source_delay);
	assert(edge);

	NUM_ELEMENTS attempts = 0;
	while (! edge->can_reduce_weight_on_edge() && attempts < m_max_comb_delay)
	{
		source_delay = (short) g_rand_number_gen->random_number(m_max_comb_delay-edge_length);

		sink_level_node = m_clusters[sink_cluster_number]->Level[edge_length+source_delay];
		assert(sink_level_node);
		edge = sink_level_node->find_edge_with_source(source_cluster_number, source_delay);
		assert(edge);
		attempts++;
	}

	if (! edge->can_reduce_weight_on_edge())
	{
		return;
	}
	assert(edge->is_inter_cluster() && edge->can_reduce_weight_on_edge());
	//debug("Source edge:  " << edge->get_source_cluster_number() << " to " << edge->get_sink_cluster_number());

	move.set_source_edge(edge);

	// choose a destination sink cluster number and delay level
	source_delay = (short) g_rand_number_gen->random_number(m_max_comb_delay-edge_length);
	assert(source_delay+edge_length >= 0 && source_delay+edge_length <= m_max_comb_delay);
	sink_level_node = m_clusters[dest_sink_cluster_number]->Level[edge_length+source_delay];
	assert(sink_level_node);
	edge = sink_level_node->find_edge_with_source(dest_source_cluster_number, source_delay);
	assert(edge && edge->is_inter_cluster());

	move.set_sink_edge(edge);
}


// get the number of edges between two clusters
//
// PRE: source_cluster_number and sink_cluster_number are valid
// POST: number_edges_between_two_clusters is the number of edges between the two clusters
//       for each edge length
//       number_poss_edges_between_two_clusters is the number of edges possible between the two 
//       clusters at each edge length
//
void SEQUENTIAL_LEVEL::get_number_edges_between_two_clusters
(
	const CLUSTER_NUMBER_TYPE& source_cluster_number,
	const CLUSTER_NUMBER_TYPE& sink_cluster_number,
	DISTRIBUTION & number_edges_between_two_clusters,
	DISTRIBUTION & number_poss_edges_between_two_clusters
)
{
	LENGTH_TYPE edge_length_size = m_max_comb_delay + 1;
	LENGTH_TYPE edge_length = 0;
	CLUSTER * sink_cluster = m_clusters[sink_cluster_number];
	LEVEL_NODE * sink_level_node = 0;
	EDGE * edge = 0;
	DELAY_TYPE source_delay = 0;

	//debug("source/sink " << source_cluster_number << " / " << sink_cluster_number);
	for (edge_length = 1; edge_length < edge_length_size; edge_length++)
	{
		//debug("edge length = " << edge_length);
	for (source_delay = 0; source_delay < m_max_comb_delay - edge_length; source_delay++)
	{
		sink_level_node = sink_cluster->Level[source_delay+edge_length];
		assert(sink_level_node);
		edge = sink_level_node->find_edge_with_source(source_cluster_number, source_delay);
		assert(edge);

		number_edges_between_two_clusters[edge_length] += edge->get_weight();
		number_poss_edges_between_two_clusters[edge_length] += edge->get_max_weight();
	}
		//debug("number/poss " << number_edges_between_two_clusters[edge_length] << " / " <<
				//number_poss_edges_between_two_clusters[edge_length]);
	}
}

// make the connections between the level nodes with latched indiv. nodes
// and level nodes with flip-flops
// 
//
// From the thesis...
//
// With the global combinational delay structure complete, 
// the delay structure is finished by forming the connections between 
// the level nodes with latched nodes and the level nodes with flip-flops. 
// The input into this phase of the algorithm is the matrix Latched and the 
// latched shapes in each cluster. The output is the finished delay structure graph. 
//
// Each connection in the graph is formed in three steps. First, we randomly select a source and 
// sink cluster from the Latched matrix. Next, we randomly select the source level node from 
// the set of all level nodes with latched nodes as dictated by the latched shape in the source cluster. 
// Finally, a connection is made between the source level node and the 0 th delay level in the sink cluster. 
//
// In actuality, we do something equivalent but slightly difference from above.
//
// We iterative through the level nodes with latched nodes.
// For such a level node we randomly pick the destination clusters
// they connect to and then make the latched level node to flip-flop 
// level node connections.
//
// It's just easier this way.
// 
//
// POST: we have completed the delay structure

void SEQUENTIAL_LEVEL::make_latched_to_dff_connections()
{
	CLUSTER* cluster = 0;
	CIRCUIT_SPEC * circuit_spec = m_circuit->get_circuit_spec();
	LEVEL_NODE * source_level_node = 0;
	DELAY_TYPE delay_level = 0;
	DISTRIBUTION destination_cluster_numbers;
	MATRIX * Latched = 0;
	MATRIX_TYPE row_matrix;
	NUM_ELEMENTS nLatched = 0;
	NUM_ELEMENTS source_cluster = 0,
				destination_cluster = 0;

	assert(circuit_spec);
	Latched = circuit_spec->get_Latched();
	assert(Latched && Latched->is_positive());

	debug("Status: Creating the large scale connections to flip-flops\n");

	// for all clusters assigned all 
	for (source_cluster = 0; source_cluster < m_nClusters; source_cluster++)
	{
		cluster = m_clusters[source_cluster];
		assert(cluster);

		// the row matrix contains all of the destination dff locations for the source cluster
		// we will sample this distribution randomly to select which level nodes 
		// will connect to what level nodes with dffs

		row_matrix = Latched->get_row(source_cluster);

		// iterate down the delay levels looking for latched nodes
		for (delay_level = 0; delay_level <= m_max_comb_delay; delay_level++)
		{
			source_level_node = cluster->Level[delay_level];
			assert(source_level_node);

			nLatched = source_level_node->get_nLatched();

			while (nLatched > 0)
			{
				// the level node has latched simple nodes.
				// pick a destination cluster at random and add it to the list of 
				// destination clusters
				destination_cluster = g_rand_number_gen->discrete_pmf(m_nClusters, row_matrix);
				assert(destination_cluster >= 0 && destination_cluster < m_nClusters);
				row_matrix[destination_cluster] -= 1;
				assert(row_matrix[destination_cluster] >= 0);

				destination_cluster_numbers.push_back(destination_cluster);

				nLatched--;
			}

			// now add connections between the source level node and the destination 
			// dffs
			create_dff_edges(source_level_node, destination_cluster_numbers);
			assert(destination_cluster_numbers.empty());
		}

		assert(accumulate(row_matrix.begin(), row_matrix.end(), 0) == 0);
	}

	sanity_check_dff_connections();
}

//
// create the super edges for the level node with latches to 
// level nodes with flip-flops
//
//
// PRE: source_level_node is the level node that has latched nodes
//      destination_cluster_numbers is the list of destination clusters that we 
//       want to make connections to.
//
// POST: we have made connections between the source level node and the level nodes 
//       in all clusters listed in destination_cluster_numbers
//
void SEQUENTIAL_LEVEL::create_dff_edges
(
	LEVEL_NODE * source_level_node,
	DISTRIBUTION & destination_cluster_numbers
)
{
	assert(m_circuit && source_level_node);
	EDGE * edge = 0;
	CLUSTER_NUMBER_TYPE destination_cluster = 0;
	NUM_ELEMENTS edge_weight = 0;
	LEVEL_NODE * dff_sink_node = 0;

	// sort the destination cluster numbers so that we can efficiently 
	// make edge connections
	sort(destination_cluster_numbers.begin(), destination_cluster_numbers.end());

	while (! destination_cluster_numbers.empty())
	{

		destination_cluster = destination_cluster_numbers.back();
		destination_cluster_numbers.pop_back();
		edge_weight = 1;

		// find out the weight of the edge that will connect the latched level node to the 
		// dff level node
		while (! destination_cluster_numbers.empty() && 
				destination_cluster == destination_cluster_numbers.back())
		{
			destination_cluster_numbers.pop_back();
			edge_weight++;
		}

		assert(destination_cluster >= 0 && static_cast<unsigned>(destination_cluster) < m_clusters.size());
		dff_sink_node = m_clusters[destination_cluster]->Level[0];
		assert(dff_sink_node);

		// create the edge
		debugif(DSEQ_LEVEL, "Creating a dff conection between clusters " 
				<< source_level_node->get_cluster_number() << " and " 
				<< destination_cluster << " from source delay level "  << 
				source_level_node->get_delay_level() << " with weight " << edge_weight);

		edge = create_dff_edge(source_level_node, dff_sink_node);
		edge->set_max_weight(edge_weight);
		edge->set_weight(edge_weight);
		dff_sink_node->input_edge_has_changed(edge, edge_weight);
		
		// don't update the source because in node splitting and in degrees we first
		// concern outselves with the combinational structure 1st then look towards dff.
	}
}

// POST: we have performed a sanity check on all clusters
void SEQUENTIAL_LEVEL::sanity_check_all() const
{
	// check to make sure everything is ok
	CLUSTER* cluster = 0;
	CLUSTERS::const_iterator cluster_iter;

	for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);

		sanity_check(cluster);
	}
}

// POST: we are sure that the input weight coming into the 0th delay levels is the 
//       number of flip-flops in the level node
void SEQUENTIAL_LEVEL::sanity_check_dff_connections() const
{
	// check to make sure everything is ok
	CLUSTER* cluster = 0;
	CLUSTERS::const_iterator cluster_iter;
	LEVEL_NODE * dff_level_node = 0;
	NUM_ELEMENTS input_weight = 0;

	for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);

		dff_level_node = cluster->Level[0];
		assert(dff_level_node);

		input_weight = dff_level_node->get_nIntra_cluster_input() + dff_level_node->get_nInter_cluster_input();

		assert(input_weight == dff_level_node->get_nDFF());
	}
}

//
// create a state of the simulated anneal
// the use of this function will slow things down
//
// POST: m_edge_weight_state contains the weight assigned to each edge

void SEQUENTIAL_LEVEL::create_edge_length_state()
{
	EDGES::iterator edge_iter;
	EDGE * edge = 0;
	NUM_ELEMENTS edge_weight = 0;
	m_edge_weight_state.clear();

	for (edge_iter = m_edges.begin(); edge_iter != m_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);
		edge_weight = edge->get_weight();
		assert(edge_weight >= 0);

		m_edge_weight_state.push_back(edge_weight);
	}
}


// recreate the lowest cost solution by reassigning the lowest cost state's
// edges weight to the graph
//
// PRE:  m_edge_weight_state is the lowest cost state
// POST: our combinational level graph structure is in the last
//       saved lowest cost state
void SEQUENTIAL_LEVEL::assign_lowest_cost_edge_weights()
{
	EDGES::iterator edge_iter;
	EDGE * edge = 0;
	NUM_ELEMENTS edge_weight = 0;
	NUM_ELEMENTS edge_number = 0;

	for (edge_iter = m_edges.begin(); edge_iter != m_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);
		assert(edge->get_id() == edge_number);

		edge_weight = m_edge_weight_state[edge_number];
		assert(edge_weight >= 0);

		edge->set_weight_and_update_level_nodes(edge_weight);
		edge_number++;
	}
}
//
// removes edges from level nodes that have too many input edges
// 
// examines each level node to see if it has too many input edges coming 
// into the node.  if it does, reduce weight on some of the super edges 
// input into the level node until either we have seen all the level nodes or
// until we fail
//
bool SEQUENTIAL_LEVEL::eliminate_excess_number_of_inputs()
{
	assert(new_gen);

	CLUSTERS::iterator cluster_iter;
	CLUSTER * cluster = 0;
	DELAY_TYPE delay = 0;
	LEVEL_NODE * level_node = 0;
	NUM_ELEMENTS excess_number_of_inputs = 0;
	bool success = true;

	for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);

		for (delay = 0; delay <= m_max_comb_delay; delay++)
		{
			level_node = cluster->Level[delay];
			assert(level_node);

			excess_number_of_inputs = level_node->get_excess_number_of_inputs();
			
			if (excess_number_of_inputs > 0)
			{
				success = level_node->reduce_weight_on_input_edges_to_eliminate_double_edge_connections();


				if (! success)
				{
					return success;
				} 

				debug("  Reduced weight on an input edge into level node in cluster " <<
						cluster->get_cluster_number() << " delay level " << level_node->get_delay_level() 
						<< " because " << endl << 
						"  the input edge had too much weight and would have caused node violations.");
			}
		}
	}

	return success;
}

// add edges to the level nodes that need more outputs 
// if an edge cannot be added to those nodes try to add PO
//
// 1st try to use edges that had weight reduced from them to solve
//     but for now just try the edges in the cluster
// 2nd just try to add weight to an edge out of the problem level node
//     
// 
// PRE: we may have level nodes that need more outputs
// POST: we have tried to add an edge to level nodes that need outputs
//       if we couldn't, we tried to add primary outputs to those level nodes
//       if we couldn't, we failed and exited the program
//
void SEQUENTIAL_LEVEL::add_edges_to_nodes_that_need_outputs()
{
	LEVEL_NODE * level_node = 0;
	DELAY_TYPE delay = 0;
	NUM_ELEMENTS outputs_needed = 0,
				 weight = 0;
	CLUSTERS::iterator cluster_iter;
	CLUSTER * cluster = 0;

	EDGE * edge = 0;
	CLUSTER_NUMBER_TYPE cluster_number = 0;

 	EDGES lost_edges = m_circuit->get_lost_weight_edges(),
 		  lost_edges_for_cluster;


	for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);

		cluster_number = cluster->get_cluster_number();
 		lost_edges_for_cluster = get_lost_edges_for_cluster(lost_edges, cluster_number);

		// try to assign outputs at the lower parts of the graph first
		// as they will be the toughest to find outputs for
		for (delay = m_max_comb_delay; delay >= 0; delay--)
		{
			level_node = cluster->Level[delay];
			assert(level_node);

			outputs_needed = level_node->get_nOutputs_still_needed();

			if (outputs_needed > 0)
			{
				debug("  Adding " << outputs_needed << " output edges to cluster " 
						<< cluster_number << " delay level " << delay << " to fix too few output problem"); 
			}
			while (outputs_needed > 0)
			{
				// try to find a lost edge that will serve our purpose.
				// otherwise just find an edge that will serve our purpose
				edge = choose_edge_to_add_output_weight_to_level_node(level_node, lost_edges_for_cluster);

				if (edge)
				{
					assert(edge->get_source_level_node() == level_node);
					Verbose("Adding weight to edge " << edge->get_info());

					weight = edge->get_weight() + 1;
					edge->unsafe_set_weight(weight);
					edge->get_source_level_node()->output_edge_has_changed(edge, 1);
					edge->get_sink_level_node()->input_edge_has_changed(edge, 1);
					outputs_needed--;
				}
				else
				{

					if (delay == 0)
					{
						Fail("Tried to add a primary output to a primary input");
					}
					// maybe we should go latch hunting instead?
					debug("  Needed to add " << outputs_needed << " PO to level node" << 
						  " in cluster " << cluster->get_cluster_number() << 
						  " delay level " << level_node->get_delay_level());
					level_node->add_PO(outputs_needed);
					outputs_needed = 0;
				}
			}
			assert(level_node->get_nOutputs_still_needed() == 0);
		}
	}	
}



// POST: the graph has been set into the last saved lowest cost state
//       the violation matrix has been recalculated
//       the edge length matrices have been recalculated
//
void SEQUENTIAL_LEVEL::set_graph_into_lowest_cost_state()
{
	assign_lowest_cost_edge_weights();
	create_violation_matrix();
	create_edge_length_matrix();
}


// POST: all congestion factors have been reset to their default value
void SEQUENTIAL_LEVEL::reset_congestion_costs()
{

	DELAY_TYPE delay_level = 0;
	DELAY_LEVEL::iterator delay_level_iter;
	LEVEL_NODE * level_node = 0;


	for (delay_level = 0; delay_level <= m_max_comb_delay; delay_level++)
	{
		for (delay_level_iter = m_delay_levels[delay_level]->begin(); 
			 delay_level_iter != m_delay_levels[delay_level]->end(); 
			 delay_level_iter++)
		{
			level_node = *delay_level_iter;
			assert(level_node);

			level_node->reset_congestion_costs();
		}

	}
}


// Increase the congestion cost for the level nodes
//
// PRE: 
// POST: each level nodes that has 
//       - too many inputs
//       - not enough inputs
//       - not enough outputs
//       have had their congestion factors for that cost increased
//       
//       otherwise they have had their congestion factors for that cost decreased
//
//
void SEQUENTIAL_LEVEL::increase_congestion_cost_for_level_nodes()
{

	DELAY_TYPE delay_level = 0;
	DELAY_LEVEL::iterator delay_level_iter;
	LEVEL_NODE * level_node = 0;
	NUM_ELEMENTS nLevel_nodes_with_excess_number_of_inputs = 0,
				 nLevel_nodes_with_Outputs_still_needed = 0,
				 nLevel_nodes_with_too_few_inputs = 0;
	NUM_ELEMENTS nInputs_in_excess = 0,
				 nOutputs_still_needed = 0,
				 nInputs_still_needed = 0;
	NUM_ELEMENTS Total_nInputs_in_excess = 0,
				 Total_nOutputs_still_needed = 0,
				 Total_nInputs_still_needed = 0;


	for (delay_level = 0; delay_level <= m_max_comb_delay; delay_level++)
	{
		for (delay_level_iter = m_delay_levels[delay_level]->begin(); 
			 delay_level_iter != m_delay_levels[delay_level]->end(); 
			 delay_level_iter++)
		{
			level_node = *delay_level_iter;
			assert(level_node);

			nInputs_in_excess = level_node->get_excess_number_of_inputs();
			nOutputs_still_needed = level_node->get_nOutputs_still_needed();
			nInputs_still_needed = level_node->get_nInputs_still_wanted();

			if (nInputs_in_excess > 0)
			{
				nLevel_nodes_with_excess_number_of_inputs ++;
				level_node->increase_too_many_inputs_factor();
			}
			else
			{
				level_node->decrease_too_many_inputs_factor();
			}
			if (nOutputs_still_needed > 0)
			{
				nLevel_nodes_with_Outputs_still_needed++;
				level_node->increase_too_few_outputs_factor();
			}
			else
			{
				level_node->decrease_too_few_outputs_factor();
			}
			if (nInputs_still_needed > 0)
			{
				nLevel_nodes_with_too_few_inputs++;
				level_node->increase_too_few_inputs_factor();
			}
			else
			{
				level_node->decrease_too_few_inputs_factor();
			}
			if (nInputs_still_needed > 0)
			{
				nLevel_nodes_with_too_few_inputs++;
				level_node->increase_too_few_inputs_factor();
			}
			else
			{
				level_node->decrease_too_few_inputs_factor();
			}

			Total_nInputs_in_excess += nInputs_in_excess;
			Total_nOutputs_still_needed += nOutputs_still_needed;
			Total_nInputs_still_needed += nInputs_still_needed;


		}
	}

	if (nLevel_nodes_with_excess_number_of_inputs > 0)
	{
		debug("There are " << nLevel_nodes_with_excess_number_of_inputs << 
				" level nodes with " << Total_nInputs_in_excess << 
				" too many inputs. Increasing too many inputs penalty factor.");
	}
	if (nLevel_nodes_with_Outputs_still_needed > 0)
	{
		debug("There are " << nLevel_nodes_with_Outputs_still_needed << 
				" level nodes with " << Total_nOutputs_still_needed << 
				" too few outputs. Increasing too few outputs penalty factor.");
	}
	if (nLevel_nodes_with_too_few_inputs > 0)
	{
		debug("There are " << nLevel_nodes_with_too_few_inputs << 
				" level nodes with " << Total_nInputs_still_needed << 
				" too few inputs. Increasing too few inputs penalty factor.");
	}
}


// Reset all the edges to have their maximum possible weight
//
//
// POST: all edges have their maximum possible weight
// 
void SEQUENTIAL_LEVEL::reset_graph()
{
	EDGES::iterator edges_iterator;
	EDGE * edge = 0;
	NUM_ELEMENTS max_weight = 0;

	for (edges_iterator = m_edges.begin(); edges_iterator != m_edges.end(); edges_iterator++)
	{
		edge = *edges_iterator;
		assert(edge);

		max_weight = edge->get_max_weight();

		edge->set_weight_and_update_level_nodes(max_weight);
	}
}


// Choose an edge to add weight to so that we can 
// eliminate node violations due to not having enough output edges
//
// 1st try the lost_edges_for_cluster
// 2nd try any old edge
// 
// PRE:level_node is valid
//     edges_for_cluster is a list of edges in the cluster
// RETURN: an edge we can add output weight to, else 0
//       
EDGE * SEQUENTIAL_LEVEL::choose_edge_to_add_output_weight_to_level_node
(
	LEVEL_NODE * level_node,
	EDGES& lost_edges_for_cluster
)
{
	assert(level_node);
	EDGES::iterator edges_iterator;
	EDGE * edge = 0;
	EDGE * edge_found = 0;
	bool success = false;
	LENGTH_TYPE edge_length;
	DELAY_TYPE delay_level = level_node->get_delay_level();
	CLUSTER_NUMBER_TYPE source_cluster_number = level_node->get_cluster_number(),
						sink_cluster_number = 0;

	// try all the lost edges - the edges that lost weight when we 
	// tried to eliminate node violations
	for (edges_iterator = lost_edges_for_cluster.begin(); 
		! success && edges_iterator != lost_edges_for_cluster.end(); edges_iterator++)
	{
		edge = *edges_iterator;
		assert(edge);

		edge_length = edge->get_length();

		if (edge_length + delay_level <= m_max_comb_delay)
		{
			// this edge will fit into this level node

			if (edge->is_intra_cluster())
			{
				// see if a level node internal to this cluster can 
				// handle an edge of this length from this level node
				
				edge_found = 
					find_an_edge_to_increase_weight(source_cluster_number, source_cluster_number, 
													edge_length, delay_level);

				success = (edge_found != 0);
			}
			else
			{
				assert(edge->is_inter_cluster());

				// things get more tricky.
				// we have an edge length and multiple clusters
				// right now major hack
				// later we can refine
				assert(new_gen);

				for (sink_cluster_number = 0; ! success && sink_cluster_number < m_nClusters; 
						sink_cluster_number++)
				{
					edge_found = find_an_edge_to_increase_weight(source_cluster_number, sink_cluster_number, 
														edge_length, delay_level);
					success = (edge_found != 0);
				}
			}
		}
	}

	if (success)
	{
		assert(edge_found->get_source_level_node() == level_node);
		return edge_found;
	}
	else
	{
		// ok. look over all edge length and clusters for an edge to add 
		// this is bad.
		Verbose("No lost edge could be found to add to output. Looking over all edge lengths and clusters");

		Verbose("First look to intra-cluster edges");

		for (edge_length = 1; ! success && edge_length <= m_max_comb_delay-delay_level; edge_length++)
		{
			edge_found = find_an_edge_to_increase_weight(source_cluster_number, source_cluster_number, 
														edge_length, delay_level);
			success = (edge_found != 0);
		}

		if (success)
		{
			assert(edge_found->get_source_level_node() == level_node);
			return edge_found;
		}

		//debug("Next look to inter-cluster edges");
		for (edge_length = 1; ! success && edge_length <= m_max_comb_delay-delay_level; edge_length++)
		for (sink_cluster_number = 0; ! success && sink_cluster_number < m_nClusters; sink_cluster_number++)
		{
			edge_found = find_an_edge_to_increase_weight(source_cluster_number, sink_cluster_number, 
														edge_length, delay_level);
			success = (edge_found != 0);
		}

		Verbose("No edge could be added to. We are screwed.");
	}

	return 0;
}


//
// PRE: cluster_number is valid
// RETURNS: the edges in lost_edges that have a source of cluster_number
//
EDGES SEQUENTIAL_LEVEL::get_lost_edges_for_cluster
(
	const EDGES& lost_edges,
	const CLUSTER_NUMBER_TYPE& cluster_number
) const
{
	EDGES lost_edges_for_cluster;
	EDGES::const_iterator edges_iterator;
	EDGE * edge = 0;

	for (edges_iterator = lost_edges.begin(); edges_iterator != lost_edges.end(); edges_iterator++)
	{
		edge = *edges_iterator;
		assert(edge);

		if (edge->get_source_cluster_number() == cluster_number)
		{
			lost_edges_for_cluster.push_back(edge);
		}
	}
	return lost_edges_for_cluster;
}

// PRE: source_cluster_number, source_cluster_number, 
//      edge_length, source_delay_level are all valid
// RETURN: an edge that can have its weight increased if one could be found, 
//         0 otherwise
//
EDGE * SEQUENTIAL_LEVEL::find_an_edge_to_increase_weight
(
	const CLUSTER_NUMBER_TYPE& source_cluster_number,
	const CLUSTER_NUMBER_TYPE& sink_cluster_number,
	const LENGTH_TYPE& edge_length,
	const DELAY_TYPE& source_delay_level
) const
{
	assert(sink_cluster_number >= 0 && sink_cluster_number < m_nClusters);

	DELAY_TYPE sink_delay_level = source_delay_level + edge_length;
	assert(sink_delay_level >= 0 && sink_delay_level <= m_max_comb_delay);
	LEVEL_NODE * sink_level_node = m_clusters[sink_cluster_number]->Level[sink_delay_level];
	NUM_ELEMENTS number_open_inputs = 0;
	EDGE * edge_found = 0;

	//debug("We are looking at an edge from source cluster " << source_cluster_number << 
	//" to sink cluster " << sink_cluster_number << " with edge length " << edge_length << 
	//" at source delay " << source_delay_level << " and sink delay " << sink_delay_level);

	number_open_inputs = sink_level_node->get_nOpen_inputs();

	//debug("The sink level node has " << number_open_inputs << " open inputs" << 
	//" as it has " << sink_level_node->get_nNodes() << " nodes and " 
	//<< sink_level_node->get_input_edge_weight() << " simple inputs");

	if (number_open_inputs > 0)
	{
		edge_found = sink_level_node->find_edge_with_source(source_cluster_number, source_delay_level);
		assert(edge_found);
	}

	return edge_found;
}

// try to fix the node violations in the graph
//
// POST: if sucessful we have no level node 
//       with too many inputs or not enough outputs
//     
// RETURNS: true if sucessful
bool SEQUENTIAL_LEVEL::try_to_fix_graph()
{

	bool success = eliminate_excess_number_of_inputs();

	if (! success)
	{
		return false;
	}
	//m_circuit->print_lost_edge_stats();

	add_edges_to_nodes_that_need_outputs(); 

	return success;
}

//
// PRE: nothing
// POST: we have reported on the quality of the solution
//       m_violation_matrix has been updated
//
void SEQUENTIAL_LEVEL::report_solution_quality()
{
	CLUSTERS::iterator cluster_iter;
	CLUSTER * cluster = 0;

	NUM_ELEMENTS input_output_shape_cost = 0,
				edge_length_cost = 0,
				cost_Comb = 0,
				po_cost = 0,
				total_cost = 0;
	double norm_total_cost = 0.0;

	create_violation_matrix();

	cost_Comb = m_violation_matrix.get_absolute_sum();


	for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);

		po_cost += cluster->get_po_difference();
		edge_length_cost += static_cast<NUM_ELEMENTS>(cluster->get_edge_length_cost());
		input_output_shape_cost += static_cast<NUM_ELEMENTS>(cluster->get_shape_cost());
	}

	total_cost = po_cost + edge_length_cost + input_output_shape_cost + cost_Comb;

	// should (maybe) really differentiate between the edges to combinational nodes and the number of edges
	assert(new_gen);
	norm_total_cost = static_cast<double>(total_cost)/static_cast<double>(m_circuit->get_nEdges());


	debug("cost_Level_Node: " << input_output_shape_cost);
	debug("cost_Edge_length: " << edge_length_cost);
	debug("cost_PO: " << po_cost);
	debug("cost_Comb: " << cost_Comb);

	debug("total_cost: " << total_cost);
	debug("normalized_total_cost: " << norm_total_cost);
	debug(sep);
	

	/*
	 * unused debug statements
	debug(Sep);
	debug("Inter cluster connection matrix\n" << m_Comb);
	debug("Spec cluster connext matrix\n" << *m_Comb_spec);
	debug("Violation matrix\n" << m_violation_matrix);

		NUM_ELEMENTS cluster_number = 0;
		for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end(); cluster_iter++)
		{
			debug("********* Cluster " << cluster_number << " **********\n");
			cluster = *cluster_iter;
			assert(cluster);
		
			debug("Comb graph cost " << cluster->get_cost());
			debug("\tComb graph edge length cost " << cluster->get_edge_length_cost());
			debug("\tComb graph level node cost " << cluster->get_level_node_cost() << endl);
			output_difference_between_current_and_spec(cluster);
			cluster_number++;
		}	
	*/
}

/*
// this move generator is completely random
// it is a little less fast than the directed move generator
void SEQUENTIAL_LEVEL::new_generate_moveIII
(
	MOVE& move
)
{
	DISTRIBUTION edge_selection_weight_pmf;

	NUM_ELEMENTS nEdges = static_cast<NUM_ELEMENTS>(m_edges.size());
	NUM_ELEMENTS edge_number = 0;
	EDGE * source_edge = 0,
		 * sink_edge = 0;
	CLUSTER_NUMBER_TYPE source_cluster = 0;
	LENGTH_TYPE source_edge_length = 0;

	edge_number = g_rand_number_gen->random_number(nEdges-1);

	source_edge = m_edges[edge_number];
	assert(source_edge);
	move.set_source_edge(source_edge);


	if (source_edge->is_inter_cluster())
	{
		// just randomly select the translocation location
		source_edge_length = edge->get_length();
		nEdges = m_inter_cluster_edges_by_length[source_edge_length].get_size();

		edge_number = g_rand_number_gen->random_number(nEdges-1);
		sink_edge = m_inter_cluster_edges_by_length[source_edge_length][edge_number];
		assert(sink_edge);
		while (! sink_edge_is_ok_to_move_to(source_edge,sink_edge))
		{
			edge_number = g_rand_number_gen->random_number(nEdges-1);
			sink_edge = m_inter_cluster_edges_by_length[source_edge_length][edge_number];
			assert(sink_edge);
		}
	}
	else
	{
		source_cluster = static_cast<NUM_ELEMENTS>(source_edge->get_sink_level_node()->get_cluster_number());
		nEdges = static_cast<NUM_ELEMENTS>(m_intra_cluster_edges_by_cluster[source_cluster].get_size());
		edge_number = g_rand_number_gen->random_number(nEdges-1);
		sink_edge = m_intra_cluster_edges_by_cluster[source_cluster][edge_number];
		assert(sink_edge);
		while (! sink_edge_is_ok_to_move_to(source_edge,sink_edge))
		{
			edge_number = g_rand_number_gen->random_number(nEdges-1);
			sink_edge = m_intra_cluster_edges_by_cluster[source_cluster][edge_number];
			assert(sink_edge);
		}
	}

	move.set_sink_edge(sink_edge);
}

// this new move generator looks at all the edges and compute selection weights
// to choose which edge will give up a simple edge
//
// it is not as fast as the directed move generator
void SEQUENTIAL_LEVEL::new_generate_moveII
(
	MOVE& move
)
{
	DOUBLE_VECTOR edge_selection_weight_pmf;

	EDGES::iterator edge_iter;
	NUM_ELEMENTS edge_number = 0;
	NUM_ELEMENTS nEdges = static_cast<NUM_ELEMENTS>(m_edges.size());
	COST_TYPE min_selection_weight = static_cast<COST_TYPE>(MAXLONG),
			max_selection_weight = 0,
			selection_weight = 0;
	EDGE * source_edge = 0,
		 * sink_edge = 0;
	// find the min and max selection weights
	for (edge_iter = m_edges.begin(); edge_iter != m_edges.end(); edge_iter++)
	{
		if ((*edge_iter)->get_weight() > 0) 
		{
			selection_weight = (*edge_iter)->get_selection_weight();
			if (selection_weight < min_selection_weight)
			{
				min_selection_weight = selection_weight;
			}
			if (selection_weight > max_selection_weight)
			{
				max_selection_weight = selection_weight;
			}
		}
	}

	// if our min selection_weight is less than 0, find the bias to all the edge selection_weights.
	if (min_selection_weight < 0.0)
	{
		min_selection_weight = -(min_selection_weight);
	}
	else
	{
		min_selection_weight = 0;
	}
	min_selection_weight += (min_selection_weight + max_selection_weight)*0.05;

	// create the pmf distribution
	for (edge_iter = m_edges.begin(); edge_iter != m_edges.end(); edge_iter++)
	{
		// check if we have any edges to remove
		if ((*edge_iter)->get_weight() == 0)
		{
			selection_weight = 0;
		}
		else
		{
			//selection_weight = (*edge_iter)->get_selection_weight();
			selection_weight = 1;
			selection_weight = selection_weight + min_selection_weight;
		}

		assert(selection_weight >= 0);
		edge_selection_weight_pmf.push_back(selection_weight);
	}
	edge_number = g_rand_number_gen->discrete_pmf(nEdges, edge_selection_weight_pmf);

	source_edge = m_edges[edge_number];
	assert(source_edge);
	move.set_source_edge(source_edge);

	// just randomly select the translocation location
	edge_number = g_rand_number_gen->random_number(nEdges-1);
	sink_edge = m_edges[edge_number];
	assert(sink_edge);
	while (! sink_edge_is_ok_to_move_to(source_edge,sink_edge))
	{
		edge_number = g_rand_number_gen->random_number(nEdges-1);
		sink_edge = m_edges[edge_number];
		assert(sink_edge);
	}
	move.set_sink_edge(sink_edge);
}

bool SEQUENTIAL_LEVEL::sink_edge_is_ok_to_move_to
(
	const EDGE * source_edge,
	const EDGE * sink_edge
) const
{
	bool sink_edge_is_ok = (source_edge->get_length() == sink_edge->get_length());

	if (source_edge->is_inter_cluster())
	{
		sink_edge_is_ok = sink_edge_is_ok && sink_edge->is_inter_cluster();
	}
	else
	{
		assert(source_edge->is_intra_cluster());
		sink_edge_is_ok = sink_edge_is_ok && sink_edge->is_intra_cluster() &&
						  (source_edge->get_source_cluster_number() == sink_edge->get_source_cluster_number());
	}

	return sink_edge_is_ok;
}
*/

