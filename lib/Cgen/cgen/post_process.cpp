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



#include "post_process.h"
#include "graph_medic.h"


POST_PROCESS::POST_PROCESS()
{
	m_circuit	= 0;
	m_nExcess_input_nodes_fixed = 0;
	m_nDouble_input_edge_nodes_fixed = 0;
	m_nNo_output_nodes_fixed = 0;
	m_nNo_input_nodes_fixed = 0;
}
POST_PROCESS::POST_PROCESS(const POST_PROCESS & another_post_process)
{
	m_circuit	= another_post_process.m_circuit;
}

POST_PROCESS & POST_PROCESS::operator=(const POST_PROCESS & another_post_process)
{
	m_circuit	= another_post_process.m_circuit;

	return (*this);
}

POST_PROCESS::~POST_PROCESS()
{
}



// PRE: circuit is valid
//      the circuit may still have node violations
// POST: No node that has multiple connections from the same source node
//       All nodes have number of inputs <= kin
// 		 All nodes have an output or else are POs
//       All nodes have inputs or else or PIs
void POST_PROCESS::post_process
(
	CIRCUIT_GRAPH * circuit
)
{

	debug("Status: Seeing if we need to post-process the circuit");

	assert(circuit);
	m_circuit = circuit;

	GRAPH_MEDIC medic(m_circuit);

	fix_nodes_with_double_inputs();

	fix_nodes_with_too_many_inputs();

	fix_nodes_with_no_inputs_and_outputs();

	assert(new_gen); //fix_number_of_primary_outputs();

	fix_flip_flops_with_connections_to_themselves();

	medic.delete_unusable_nodes();
	
	print_warnings_as_to_what_was_fixed();

	debug("\nStatus: Post processing done\n");
}


// PRE: we have nodes with no outputs
// POST: all nodes have an output edge, else they are PO
void POST_PROCESS::fix_nodes_with_no_outputs
(
	NODE_LIST& nodes_with_no_outputs
)
{
	NODE_LIST::iterator node_iter;
	NODE * node = 0;

	for (node_iter = nodes_with_no_outputs.begin(); node_iter != nodes_with_no_outputs.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		assert(node->get_nOutputs() == 0 && ! node->is_PO());

		try_and_find_output_for_node(node);
	}
}

// PRE: we have nodes with no inptus
// POST: all nodes have at least one input, if they are not PIs
void POST_PROCESS::fix_nodes_with_no_inputs
(
	NODE_LIST& nodes_with_no_inputs
)
{
	CLUSTER * cluster = 0;
	DELAY_TYPE  lev = 0;
	LEVEL_NODE * level_node = 0;
	NODE_LIST::iterator node_iter;
	NODE * node = 0;
	NUM_ELEMENTS nClusters = m_circuit->get_nClusters(),
				 choice = 0;
	bool success = false;

	for (node_iter = nodes_with_no_inputs.begin(); node_iter != nodes_with_no_inputs.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		assert(node->get_nInputs() == 0 && ! node->is_PI());

		lev = node->get_delay_level();
		level_node = node->get_my_level_node();
		cluster = level_node->get_my_cluster();

		// try intra-cluster
		success = try_and_find_input_for_node(cluster, lev, node);

		while (! success)
		{
			choice = g_rand_number_gen->random_number(nClusters-1);
			assert(choice >= 0 && choice < nClusters);
			cluster = m_clusters[choice];
			
			success = try_and_find_input_for_node(cluster, lev, node);
		}

	}
}

// PRE: we may have nodes that have no inputs that should have inputs
//      we may have nodes that have no outputs that should have outputs
// POST: all nodes have an output or else are POs
//       all nodes have inputs or else or PIs
void POST_PROCESS::fix_nodes_with_no_inputs_and_outputs()
{
	NODES& indiv_nodes = m_circuit->get_nodes();
	NODES::iterator node_iter;
	NODE_LIST no_input_nodes, no_output_nodes;
	NODE_LIST::iterator sink_node_iter, source_node_iter;
	NODE * node = 0,
		 * source_node = 0,
		 * sink_node = 0;
	bool match_found = false;

	// find all nodes with no inputs or outputs
	for (node_iter = indiv_nodes.begin(); node_iter != indiv_nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);
		
		if (node->get_nInputs() == 0 && ! node->is_PI())
		{
			// DFF shouldn't produce these types of errors
			assert(node->is_combinational());
			no_input_nodes.push_back(node);
		}
		else if (node->get_nOutputs() == 0 && ! node->is_PO())
		{
			no_output_nodes.push_back(node);
		}
	}

	// attach the nodes with no inputs to the nodes with no outputs
	// is this right? who knows but it will hopefully prevent 
	// more node violations
	if (! no_input_nodes.empty() && ! no_output_nodes.empty())
	{
		sink_node_iter = no_input_nodes.begin();
		sink_node = *sink_node_iter;
		assert(sink_node);

		while (sink_node_iter != no_input_nodes.end())
		{

			source_node_iter = no_output_nodes.begin();
			match_found = false;

			// look to see if we can make any connections from a node with no output
			// to a node with no input
			while (! match_found && source_node_iter != no_output_nodes.end())
			{
				source_node = *source_node_iter;

				if (source_node->get_delay_level() < sink_node->get_delay_level())
				{
					// these nodes have no inputs or outputs they should not be connected
					assert(! sink_node->are_nodes_connected(source_node));

					match_found = true;
					m_circuit->create_edge(source_node, sink_node);

					// should make sure the edge has its super edge set.
					//
					assert(new_gen);

					set_max_fanout_if_necessary(source_node);

					// remove the source and sink nodes from the lists
					sink_node_iter = no_input_nodes.erase(sink_node_iter);
					no_output_nodes.erase(source_node_iter);

					m_nNo_input_nodes_fixed++;
					m_nNo_output_nodes_fixed++;
				}
				else
				{
					source_node_iter++;
				}
			}

			sink_node_iter++;
		}
	}

	// fix each of these nodes separately.
	if (! no_output_nodes.empty())
	{
		fix_nodes_with_no_outputs(no_output_nodes);
	}
	if (! no_input_nodes.empty())
	{
		fix_nodes_with_no_inputs(no_input_nodes);
	}
}

// this node does not have an input.
// find it an input
//
// PRE: node_to_find_input_for is a node that needs an input
//      cluster is the source cluster that we want to find a node that will
//      feed our node
//      delay_level_of_sink_node is the delay level of the sink node
// POST: we have found and added an input to this node
//
bool POST_PROCESS::try_and_find_input_for_node
(
	CLUSTER* cluster, 
	const DELAY_TYPE& delay_level_of_sink_node, 
	NODE * node_to_find_input_for
)
{
	assert(cluster && node_to_find_input_for);

	// should do this smarter.
	// but i don't have time. quick and dirty random method
	assert(new_gen);

	assert(node_to_find_input_for->get_nInputs() == 0 && 
			! node_to_find_input_for->is_PI());


	LEVEL_NODE * level_node = 0;
	NUM_ELEMENTS nNodes = 0,
				 choice = 0;
	DELAY_TYPE lev;
	NODE * source_node = 0;
	NODES indiv_nodes;
	bool success = false;


	while (! success)
	{

		// select a delay level for the source node
		lev = (short) g_rand_number_gen->random_number(delay_level_of_sink_node-1);
		assert(lev >= 0 && lev < delay_level_of_sink_node);

		level_node = cluster->Level[lev];
		assert(level_node);

		nNodes = level_node->get_nNodes();

		if (nNodes > 0)
		{

			// select a node from the level node
			choice = g_rand_number_gen->random_number(nNodes-1);

			indiv_nodes = level_node->get_nodes();

			source_node = indiv_nodes[choice];
			assert(source_node);

			// if the nodes are not connected make a make connection
			if (! node_to_find_input_for->are_nodes_connected(source_node))
			{
				success = true;

				debugif(DPOST_PROCESS,"Fixing " << node_to_find_input_for->get_name() << 
						" with no input by assigning it a random edge from " << source_node->get_name());
				m_nNo_input_nodes_fixed++;
				m_circuit->create_edge(source_node, node_to_find_input_for);
				set_max_fanout_if_necessary(source_node);
				m_nNo_input_nodes_fixed++;
			}
		}
	}

	return success;
}
//
// PRE: node does not have any outputs 
// POST: the node has an output, else we have made the node a primary output
//       
// note: the making a node a PO is hopefully rare
void POST_PROCESS::try_and_find_output_for_node
(
	NODE * node
)
{
	assert(node);

	bool success = false;
	LENGTH_TYPE edge_length = 0;

	CLUSTER_NUMBER_TYPE sink_cluster_number = 0;
	NUM_ELEMENTS nClusters = m_circuit->get_nClusters();
	DELAY_TYPE max_comb_delay = m_circuit->get_delay();
	CLUSTER * cluster = 0;

	LEVEL_NODE * level_node = level_node = node->get_my_level_node();
	DELAY_TYPE delay_level = node->get_delay_level();

	sink_cluster_number = node->get_cluster_number();
	debugif(DPOST_PROCESS,"First look to intra-cluster edges");
	for (edge_length = 1; ! success && edge_length <= max_comb_delay-delay_level; edge_length++)
	{
		success = find_an_open_input_and_create_an_edge(node, sink_cluster_number, edge_length, delay_level);
	}

	debugif(DPOST_PROCESS,"Next look to inter-cluster edges");
	for (edge_length = 1; ! success && edge_length <= max_comb_delay-delay_level; edge_length++)
	for (sink_cluster_number = 0; ! success && sink_cluster_number < nClusters; sink_cluster_number++)
	{
		success = find_an_open_input_and_create_an_edge(node, sink_cluster_number, edge_length, delay_level);
	}

	if (! success)
	{
		cluster = level_node->get_my_cluster();

		Warning("No node with an input could be found.");
		Warning("Setting it to be a primary output.");
		node->set_is_PO(true);
		assert(new_gen); // maybe we should add it to the level node as well
		cluster->add_PO(1);
	}
	m_nNo_output_nodes_fixed++;
}

// PRE: source node is the node we want to find an output for
//      and source_delay level is its delay level
//      sink_cluster_number is the sink cluster we want to find an input in
//      edge_length is the length of edge we want 
// POST: if we could find an open input we have created an edge
// RETURN: success
bool POST_PROCESS::find_an_open_input_and_create_an_edge
(
	NODE * source_node,
	const CLUSTER_NUMBER_TYPE& sink_cluster_number,
	const LENGTH_TYPE& edge_length,
	const DELAY_TYPE& source_delay_level
)
{
	NODE * sink_node = 0;
	EDGE * edge = 0;
	bool edge_found = false;

	assert(sink_cluster_number >= 0 && sink_cluster_number < m_circuit->get_nClusters());

	DELAY_TYPE sink_delay_level = source_delay_level + edge_length;
	LEVEL_NODE * sink_level_node = m_clusters[sink_cluster_number]->Level[sink_delay_level];

	sink_node = sink_level_node->get_indiv_node_with_open_input();

	if (sink_node && ! sink_node->are_nodes_connected(source_node))
	{
		debugif(DPOST_PROCESS,"Found a simple node with an input at level node " << sink_level_node->get_name());

		edge = m_circuit->create_edge(source_node, sink_node);
		set_max_fanout_if_necessary(source_node);

		edge_found = (edge != 0);
	}

	return edge_found;
}


// we added primary outputs.  see if we can remove some of the primary outputs 
// to reobtain our original number
//
// 1st look at nodes at the same delay level for a primary input
// 2nd look at the cluster
// 3rd look at other clusters.
// 
void POST_PROCESS::fix_number_of_primary_outputs()
{
	// haven't code this yet
	assert(false);
	assert(new_gen);


}

// PRE: We may have nodes that have multiple inputs from the same source node
// POST: No node that has multiple connections from the same source node
void POST_PROCESS::fix_nodes_with_double_inputs()
{
	CLUSTERS::iterator cluster_iter;
	CLUSTER * cluster = 0;
	LEVEL_NODE * level_node = 0;
	DELAY_TYPE	lev = 0;
	NODES indiv_nodes;
	NODES::iterator node_iter;
	NODE * node = 0,
		 * source_node = 0;
	NUM_ELEMENTS choice = 0;
	EDGE * edge = 0;

	m_clusters = m_circuit->get_clusters();
	NUM_ELEMENTS nDouble_input_edges = 0;
	EDGES double_input_edges;

	// examine the nodes in all clusters
	for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);

		for (lev = 0; lev <= cluster->get_delay(); lev++)
		{
			level_node = cluster->Level[lev];
			assert(level_node);

			// get the indiv. nodes in the clusters
			indiv_nodes = level_node->get_nodes();

			for (node_iter = indiv_nodes.begin(); node_iter != indiv_nodes.end(); node_iter++)
			{
				node = *node_iter;
				assert(node);
				
				if (node->get_nDouble_inputs() > 0)
				{
					m_nDouble_input_edge_nodes_fixed++;
				}

				// while the node has multiple connections
				// from the same source nodes remove one of those edges
				while (node->get_nDouble_inputs() > 0)
				{
					double_input_edges = node->get_double_input_edges();
					nDouble_input_edges = static_cast<NUM_ELEMENTS>(double_input_edges.size());

					debugif(DPOST_PROCESS,"Node " << node->get_name() << " has " << nDouble_input_edges 
							<< " double input edges.  Trying to remove one");

					choice = g_rand_number_gen->random_number(nDouble_input_edges-1);
					assert(choice >= 0 && choice < nDouble_input_edges);

					edge = double_input_edges[choice];
					assert(edge);

					source_node = edge->get_source_node();
					assert(source_node);

					m_circuit->remove_edge(edge);
					node->remove_input_edge(edge);
					source_node->remove_output_edge(edge);

					delete edge;
				}

				assert(node->get_nDouble_inputs() == 0);
			}
		}
	}
}

//
// PRE: We may have nodes that have too many inputs for the kin value
// POST: All nodes have number of inputs <= kin
void POST_PROCESS::fix_nodes_with_too_many_inputs()
{
	CLUSTERS::iterator cluster_iter;
	CLUSTER * cluster = 0;
	LEVEL_NODE * level_node = 0;
	DELAY_TYPE	lev = 0;
	NODES indiv_nodes;
	NODES::iterator node_iter;
	EDGES::iterator edge_iter;
	NODE * node = 0,
		 * source_node = 0;
	NUM_ELEMENTS choice = 0;
	EDGE * edge = 0;
	NUM_ELEMENTS_VECTOR source_node_output_weight;

	m_clusters = m_circuit->get_clusters();
	NUM_ELEMENTS source_node_weight = 0,
				 nInputs = 0;
	EDGES input_edges;

	for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);

		for (lev = 0; lev <= cluster->get_delay(); lev++)
		{
			level_node = cluster->Level[lev];
			assert(level_node);

			indiv_nodes = level_node->get_nodes();


			for (node_iter = indiv_nodes.begin(); node_iter != indiv_nodes.end(); node_iter++)
			{
				node = *node_iter;
				assert(node);

				// while we have nodes have too many inputs
				// remove one
				while (node->get_nExcess_inputs() > 0)
				{
					input_edges = node->get_input_edges();
					nInputs = node->get_nInputs();
					source_node_output_weight.clear();
	
					
					debugif(DPOST_PROCESS,"Node " << node->get_name() << " has " << nInputs << 
							" inputs.  Trying to fix");
					m_nExcess_input_nodes_fixed++;

					// we want to pick edges with source nodes that have 
					// lots of edges to give
					// so 1st find the number of outputs the source nodes have 
					// and then use this weight as a pmf 
					// to randonmly select from
					for (edge_iter = input_edges.begin(); edge_iter != input_edges.end(); edge_iter++)
					{
						edge = *edge_iter;
						assert(edge);

						source_node = edge->get_source_node();
						assert(source_node);

						source_node_weight = source_node->get_nOutputs();

						source_node_output_weight.push_back(source_node_weight);
					}


					choice = g_rand_number_gen->discrete_pmf(source_node_output_weight);

					edge = input_edges[choice];
					assert(edge);

					source_node = edge->get_source_node();
					assert(source_node);

					node->remove_input_edge(edge);
					source_node->remove_output_edge(edge);
				}

				assert(node->get_nExcess_inputs() == 0);
			}
		}
	}
}

// POST: warnings about what needed to be fixed have been outputed
void POST_PROCESS::print_warnings_as_to_what_was_fixed() const
{
	if (m_nExcess_input_nodes_fixed > 0 || m_nDouble_input_edge_nodes_fixed > 0 ||
		m_nNo_output_nodes_fixed > 0 || m_nNo_input_nodes_fixed > 0)
	{
		dbg("\n");
		Warning("Postprocessing was necessary.\n" <<
				"-------------------------------------");

		if (m_nExcess_input_nodes_fixed > 0)
		{
			Warning("Needed to fix " << m_nExcess_input_nodes_fixed << " nodes with excessive number of inputs");
		}
		if (m_nDouble_input_edge_nodes_fixed > 0)
		{
			Warning("Needed to fix " << m_nDouble_input_edge_nodes_fixed << " nodes with more than one connection to the same source node.");
		}
		if (m_nNo_output_nodes_fixed > 0)
		{
			Warning("Needed to fix " << m_nNo_output_nodes_fixed << " nodes with no output edges that were not PO");
		}
		if (m_nNo_input_nodes_fixed > 0)
		{
			Warning("Needed to fix " << m_nNo_input_nodes_fixed << " nodes that had no inputs that were not PI");
		}

	}
}

// sets the max fanout of the cluster if our modified source node
// has changed it
//
// PRE: source_node is a node that has had an edge added to it
// POST: if source_node's fanout is now the biggest
//       the source cluster maximum fanout has been updated
//
void POST_PROCESS::set_max_fanout_if_necessary
(
	NODE * source_node
)
{
	assert(source_node);

	CLUSTER_NUMBER_TYPE source_cluster_number = source_node->get_cluster_number();
	CLUSTER * source_cluster = m_clusters[source_cluster_number];
	NUM_ELEMENTS source_node_fanout = source_node->get_nOutputs();

	assert(source_cluster);
		
	if (source_node_fanout > source_cluster->get_max_fanout())
	{
		debugif(DPOST_PROCESS,"The fanout of node " << source_node->get_name() << " is now greater than the "
				<< "fanout of the cluster because of post-processing.");
		source_cluster->set_new_max_fanout(source_node_fanout);
	}
}


// PRE: we might have nodes that loop back to themselves
// POST: no such node exists in the circuit
void POST_PROCESS::fix_flip_flops_with_connections_to_themselves()
{
	// for now just delete the edge
	assert(new_gen);

	NODES& indiv_nodes = m_circuit->get_nodes();
	NODES::iterator node_iter;
	NODE * node = 0,
		 * source_node = 0;
	EDGE * edge = 0;

	// find all the flip-flops with loops back to themselves.
	// delete the flip-flop and the possible buffer node it connects to
	for (node_iter = indiv_nodes.begin(); node_iter != indiv_nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);
		
		if (node->is_dff_have_loop_back_to_itself())
		{
			assert(node->is_DFF());
			edge = node->get_dff_input_edge();
			assert(edge);
			source_node = edge->get_source_node();

			m_circuit->remove_edge(edge);
			node->remove_input_edge(edge);
			source_node->remove_output_edge(edge);

			// the node should be deleted by graph medic if it is unconnected to anything

			delete edge;
		}
	}
}

