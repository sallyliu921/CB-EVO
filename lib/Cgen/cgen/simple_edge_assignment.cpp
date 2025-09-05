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


#define sanity assert

#define HIGHDEG 0
#define LOC 0
#define FULL 0

#include "simple_edge_assignment.h"
#include "util.h"
#include <limits.h>
#include <algorithm>
#include <numeric>
#include "rand.h"
#include <math.h>
#include "fp.h"

const NUM_ELEMENTS CHOICELIST_SIZE = 128;

SIMPLE_EDGE_ASSIGNER::SIMPLE_EDGE_ASSIGNER(CIRCUIT_GRAPH * circuit)
{
	assert(circuit);
	m_circuit = circuit;
	m_sum_out 	= 0;
	m_sum_in	= 0;
	m_nSrc		= 0;
	m_nDst		= 0;
	m_cluster= 0;
}
SIMPLE_EDGE_ASSIGNER::~SIMPLE_EDGE_ASSIGNER()
{
	m_sum_out = 0;
	m_cluster = 0;
}
//
//  Mainline for edge-assignment part of node-splitting.
//
//  Proceed level-by-level, and assign edges from preceding levels.
//  
//  The last stage in the algorithm is final edge assignment. 
//  The input into the algorithm is the pre-edge assignment structure 
//  from node_splitter.cpp.  
//  The output from this step is the graph where the individual edges have been 
//  assigned to the indiv. nodes are where node violations may exist.
//  The initial solution is created in two parts. 
//  First, we create the connections to the combinational (individual) nodes and 
//  secondly we create the connections between the latched nodes and flip-flops.
//
//  We make the connections to combinational nodes in the graph by visiting each 
//  level node in turn and forming connections to the individual nodes it contains in 
//  four steps. Our method to construct the initial solution is based on and evolves 
//  from Hutton s work . 
//
//  In Step 1, we create a list of source and  destination nodes for the connections 
//  we want to form. The destination node list we create from the level node's individual 
//  nodes. The source node list we create by sampling individual nodes with unassigned fanout 
//  from the source level nodes, which are all the level nodes with a super edge input 
//  into the level node. The number of nodes sampled from each source level node is the weight 
//  of its corresponding super edge. When making connections, we remove nodes in the 
//  destination list once they can no longer accept an input and nodes in the source list 
//  once the numbers of connections that have been made from it equal the number of times 
//  it was sampled. 
//
//  In Step 2, we ensure that each destination node has its delay level well defined by 
//  assigning an edge between it and a source node that is a unit edge length apart. 
//  The source node is chosen by randomly sampling the source node list for nodes 
//  that are a unit distant apart and choosing the source node that is closest in 
//  horizontal position to the destination node. This process of sampling the source 
//  nodes for a node that we can make a connection to is called the edge connection process.
// 
//  In Step 3, we ensure that the destination nodes are not single-input buffers by assigning 
//  them a second edge with the edge connection process. The only restriction on this node 
//  selection process is that the destination node cannot make a second connection to the 
//  same source node. 
//
//  In Step 4, we form the remainder of the connections. We randomly select destination nodes 
//  and make connections using the edge connection process until either we can no longer make 
//  further connections without creating node violations or we run out of source nodes. If 
//  any source nodes remain, we randomly select destination nodes and make connections in 
//  spite of node violations in the hope that they will be resolved later during iterative 
//  improvement.
//
//  Differences from thesis: 
//   
//   1. We don't add all the source nodes to the list all at once.
//      We add them when they are needed
//   2.	Before the delay defining step if the source nodes have a high fanout 
//		in relation to the number of destination nodes 
//		we do edge assignment for high fanout nodes first
//		before making sure that they have their delay level defined because 
//		high fanout nodes will be hard to find edges for.
//
void SIMPLE_EDGE_ASSIGNER::assign_simple_edges(CLUSTER * cluster)
{
	assert(cluster);
	m_cluster = cluster;
	NUM_ELEMENTS nNodes = m_cluster->get_nNodes();
	NUM_ELEMENTS nStranded_nodes = 0;
	LEVEL_NODE * level_node = 0;

    DELAY_TYPE delay, lev;

    debugif(DSIMPLE_EDGE, Sep << "Assigning actual edges to the graph");

	initialize_variables(nNodes);

    delay = m_cluster->get_delay();
    for (lev = 1; lev <= delay; lev++) 
	{
		level_node = m_cluster->Level[lev];
		assert(level_node);

		// if the input weight into this level node will cause a "too many inputs" node violations
		// then fail because this should have been dealt with before this
		if (level_node->get_input_edge_weight() > level_node->get_nNodes() * m_cluster->get_kin())
		{
			Fail("Have too many inputs for the number of nodes at the level. " <<
				"This should have been fixed before getting here");
		}

		// if we have destination nodes to assign edges to 
		if (level_node->get_nNodes() > 0)
		{
			// Step 1: create the source and destination list 
			//         for the source list only add source nodes that 
			//         can be reached by a length 1 intra-cluster edge
			add_length1_intra_cluster_edges_to_source_list(lev-1);
			create_destination_list(lev);

			output_debug_msg();

			//  Before the delay defining step if the source nodes have a high fanout 
			//  in relation to the number of destination nodes 
			//	we do edge assignment for high fanout nodes first
			//	before making sure that they have their delay level defined because 
			//	high fanout nodes will be hard to find edges for.
			//  Otherwise just proceed to delay-defining step .

			if (m_sum_out > 2 * m_nDst) 
			{
				do_high_fanouts();
			}

			// Add the rest of the length 1 edges to the source list
			add_length1_inter_cluster_edges_to_source_list(lev);

			// Step 2: establish the delay level of the destination nodes 
			establish_delay_level(lev);

			// Add the long and inter-cluster edges to the source list in anticipation of Step 3
			add_long_and_inter_cluster_edges_to_source_list(lev);

			output_debug_msg();

			m_sum_in = MIN(m_sum_in, m_sum_out);

			// Step 3. Ensure that the destination nodes are not single-input 
			// buffers by assigning them a second edge
			ensure_no_buffer_nodes();

			//  Step 4. Form the remainder of the connections. 
			assign_remaining_edges();

			// warn the user about standed source nodes and edges not being assigned
			output_warning_msgs_to_user(level_node);
			nStranded_nodes += m_nSrc;

			no_unwanted_outputs(lev);

			debugif(DSIMPLE_EDGE,"\nWe assigned Level node " << level_node->get_name() << "  " <<
					level_node->get_nSimple_input_edges() << " simple input edges.");
		
			debugif(DSIMPLE_EDGE,"Phase A spec was " << level_node->get_input_edge_weight());
			debugif(DSIMPLE_EDGE,"The orig spec was " << level_node->get_assigned_in());


			m_src.clear();
			m_dst.clear();
			m_choice_list.clear();
			m_src.resize(nNodes, 0);
			m_dst.resize(nNodes, 0);
			m_choice_list.resize(CHOICELIST_SIZE, 0);
			m_nSrc = 0;
			m_nDst = 0;
			m_sum_out = 0;
			m_sum_in = 0;
		}
		else
		{
			debugif(DSIMPLE_EDGE,"Level " << lev << " has zero nodes.");
		}
	}
	// Summarize for the user the number of stranded source nodes we
	// had assigning edges in this cluster
	if (nStranded_nodes > 0)
	{
		Verbose("Stranded " << nStranded_nodes << " nodes with output connections.");
		nStranded_nodes = 0;
	}

    debugif(DSIMPLE_EDGE,"Checking if any dff exist for this cluster");

	if (m_cluster->get_nDFF() > 0)
	{
		debugif(DSIMPLE_EDGE,Sep);
    	debugif(DSIMPLE_EDGE,"We have " << m_cluster->get_nDFF() << " in this circuit");
		debugif(DSIMPLE_EDGE,"Input weight into dff is " << m_cluster->Level[0]->get_nIntra_cluster_input()
				+ m_cluster->Level[0]->get_nInter_cluster_input());
		add_dff_nodes_to_dst_list();
		add_latched_nodes_to_source_list();
		debugif(DSIMPLE_EDGE,"We have " << m_nSrc << " latched nodes and " << m_nDst << " dff to connect to");
		assert(m_nSrc == m_nDst);
		assign_dff_simple_edges();
		remove_unused_latched_nodes();
		assert(m_nSrc == 0 && m_nDst == 0);
		m_nSrc = 0;
		m_nDst = 0;
	}

    m_src.clear();
    m_dst.clear();
    debugif(DSIMPLE_EDGE, Sep);
}

// POST: variables have been initialized
//
void SIMPLE_EDGE_ASSIGNER::initialize_variables
(
	const NUM_ELEMENTS& nNodes
)
{
	debugif(DSIMPLE_EDGE, Sep);
	debugif(DSIMPLE_EDGE,"Working on cluster " << m_cluster->get_cluster_number());

    m_sum_out = m_sum_in = 0;

	m_nSrc	= 0; 
	m_nDst	= 0;
    m_src.resize(nNodes, 0);
    m_dst.resize(nNodes, 0);
    m_choice_list.resize(CHOICELIST_SIZE, 0);
}

//
//  Add the nodes from specified level to the candidate-list for edge
//  sources.   
//
//  POST: source nodes reachable by a lenght 1 intra-cluster edge have 
//        been added to the source list
void SIMPLE_EDGE_ASSIGNER::add_length1_intra_cluster_edges_to_source_list(const DELAY_TYPE& delay_level)
{
	assert(delay_level >= 0);

    LEVEL_NODE* level_node = 0;
    NODE* node	= 0;
	NUM_ELEMENTS node_index = 0;
	NUM_ELEMENTS nNodes = 0, edges_to_assign = 0, nLong_edges = 0;
   
    assert(m_nSrc == 0 || m_sum_out == 0);

    level_node = m_cluster->Level[delay_level];
	assert(level_node);
	nNodes = level_node->get_nNodes();

    debugif(DSIMPLE_EDGE,"Adding to source list intra-cluster edges from level " 
			<< delay_level << " with nodes " << nNodes);

    for (node_index = 0; node_index < nNodes; node_index++)
	{
		node = level_node->get_node(node_index);
		assert(node);

		edges_to_assign = node->get_edges_to_assign();
		nLong_edges = node->get_unassigned_long_and_inter_cluster_edges();

		debugif(DSIMPLE_EDGE, "..node " << node->get_id() << " has to_assign=" << edges_to_assign <<
				", long=" << nLong_edges << ", nout=" << node->get_fanout_degree());
		
		assert(edges_to_assign + nLong_edges <= node->get_fanout_degree());

		if (edges_to_assign > 0) 
		{ 
			add_source_node(node, edges_to_assign);
			node->add_edges_to_assign(-edges_to_assign);
			assert(node->get_edges_to_assign() == 0);
		} 
		else 
		{
			debugif(DSIMPLE_EDGE, "Ignoring node with out=0");
		}
    }

	dbgif(DSIMPLE_EDGE, "\n");
}

// add the source node to the source list
//
// PRE: node is valid, edges_to_assign > 0
// POST: if the node was not already on the list it has been added
//       the number of edges_currently_assigning has been increased
void SIMPLE_EDGE_ASSIGNER::add_source_node
(
	NODE * node,
	const NUM_ELEMENTS& edges_to_assign
)
{
	assert(node);
	assert(edges_to_assign > 0);

	if (node->get_edges_currently_assigning() == 0)
	{
		// the  node is not on the source list. add it.
		if (static_cast<unsigned>(m_nSrc) < m_src.size())
		{
			m_src[m_nSrc++] = node;
		}
		else
		{
			assert(static_cast<unsigned>(m_nSrc) == m_src.size());

			m_src.push_back(node);
			m_nSrc++;
		}
		debugif(DSIMPLE_EDGE,"Adding node " << node->get_id() << 
				" with name " << node->get_name() << " to source list at index " << m_nSrc-1);
	}
	else
	{
		debugif(DSIMPLE_EDGE,"Increasing the number of out edges for node " << node->get_name() << " on source list");
	}

	node->add_edges_currently_assigning(edges_to_assign);
	m_sum_out += edges_to_assign;
}

//
//  Add long intra-cluster and inter-cluster edges to the list
// 
// POST: source nodes that are long and are inter-cluster edges have
//       been added to the source list 
//
void SIMPLE_EDGE_ASSIGNER::add_long_and_inter_cluster_edges_to_source_list(const DELAY_TYPE& delay)
{
	assert(delay >= 0);
	LEVEL_NODE * level_node = m_cluster->Level[delay];
	assert(level_node);
	EDGES input_edges = level_node->get_input_edges();
	EDGES::iterator edge_iter;
	EDGE * edge = 0;
	NODE * node = 0;
	NODES nodes;

	debugif(DSIMPLE_EDGE,"Adding long and inter-cluster edges to source list");

	for (edge_iter = input_edges.begin(); edge_iter != input_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		if (edge->get_length() != 1 && edge->get_weight() > 0)
		{
			nodes = select_simple_nodes_from_source_cluster(edge);

			while (! nodes.empty())
			{
				node = nodes.back();
				nodes.pop_back();
				assert(node);
				
				node->add_unassigned_long_and_inter_cluster_edges(-1);

				add_source_node(node,1);
			}
		}
	}

    debugif(DSIMPLE_EDGE,"Done with allocating long and inter-cluster edges to source list," << get_status_string() << endl);
}


// 
// POST: source nodes that are inter-cluster edges a length 1 distant have
//       been added to the source list 
//       
void SIMPLE_EDGE_ASSIGNER::add_length1_inter_cluster_edges_to_source_list(const DELAY_TYPE& delay)
{
	assert(delay >= 0);
	LEVEL_NODE * level_node = m_cluster->Level[delay];
	assert(level_node);
	EDGES input_edges = level_node->get_input_edges();
	EDGES::const_iterator edge_iter;
	EDGE * edge = 0;
	NODE * node = 0;
	NODES nodes;
	NUM_ELEMENTS length_1_added = 0;

	debugif(DSIMPLE_EDGE,"Adding length 1 inter cluster edges to source list"); 

	for (edge_iter = input_edges.begin(); edge_iter != input_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		if (edge->is_inter_cluster() && edge->get_length() == 1 &&
			edge->get_weight() > 0)

		{
			nodes = select_simple_nodes_from_source_cluster(edge);

			while (! nodes.empty())
			{
				node = nodes.back();
				nodes.pop_back();
				assert(node);
				
				// node is not on the list
				// add it to the list
				node->add_unassigned_long_and_inter_cluster_edges(-1);

				add_source_node(node,1);

				length_1_added++;
			}
		}
	}

	debugif(DSIMPLE_EDGE,"Added " << length_1_added << " inter cluster length 1 edges to the source list");
	debugif(DSIMPLE_EDGE,"Have " << m_nSrc << " nodes with " << m_sum_out << " edges out; " 
					<< m_nDst << " nodes with " << m_sum_in << " possible edges in.\n");
}


// return indiv. nodes from the source level node of the super edge.
// the number of simple nodes to return is the weight of the super edge.
//
// if not enough simple nodes with long edges exist then reduce the number of 
// simple nodes to return
//
// PRE: super_edge is valid
// RETURN: simple nodes that can connect the indiv. nodes in the sink level node
//
NODES SIMPLE_EDGE_ASSIGNER::select_simple_nodes_from_source_cluster
(
	const EDGE * super_edge 
) const
{
	LEVEL_NODE * source_level_node = super_edge->get_source_level_node();
	NUM_ELEMENTS nSimple_edges_to_assign = super_edge->get_weight();
	DISTRIBUTION nLongDistribution = 
		source_level_node->get_number_of_unassigned_long_and_inter_cluster_edges_each_indiv_node_has();
	NUM_ELEMENTS nNodes = source_level_node->get_nNodes();
	NUM_ELEMENTS choice = 0;
	NODE * node = 0;
	NODES nodes;
	NUM_ELEMENTS total_available = 0;

	// to prevent double connections limit the nLong to the 
	// number of nodes in the destination level node
	NUM_ELEMENTS max_nLong = super_edge->get_sink_level_node()->get_nNodes();

	for (choice = 0; static_cast<unsigned>(choice) < nLongDistribution.size(); choice++)
	{
		nLongDistribution[choice] = MIN(nLongDistribution[choice], max_nLong);
	}

	// we will just randomly select one of the simple nodes that have 
	// long edges assigned it as the source simple node 

	total_available = accumulate(nLongDistribution.begin(), nLongDistribution.end(), 0);
	assert(total_available >= 0);
	
	debugif(DSIMPLE_EDGE,"Need " << nSimple_edges_to_assign << " simple edges from level node "
			<< source_level_node->get_name());
	debugif(DSIMPLE_EDGE,"Have " << total_available << " long edges available");

	if (nSimple_edges_to_assign > total_available)
	{
		Verbose("Super edge wanted " << (nSimple_edges_to_assign - total_available) << 
				" more simple nodes than available.");
	}
	nSimple_edges_to_assign = MIN(nSimple_edges_to_assign, total_available);

	assert(new_gen);	
	//assert(accumulate(nLongDistribution.begin(), nLongDistribution.end(), 0) >= nSimple_edges_to_assign);
	//below is a problem that shouldn't happen
	while (nSimple_edges_to_assign > 0)
	{
		choice = g_rand_number_gen->discrete_pmf(nNodes, nLongDistribution);
		assert(choice >= 0 && choice <= nNodes);
		assert(nLongDistribution[choice] > 0);
		nLongDistribution[choice]--;
		
		node = source_level_node->get_node(choice);

		nodes.push_back(node);

		nSimple_edges_to_assign--;

	}

	return nodes;
}

//
//  Add the nodes at this level to the list of nodes requiring fanout.
//  
//  POST: the destination list has been created 
//
void SIMPLE_EDGE_ASSIGNER::create_destination_list(const DELAY_TYPE& delay_level) 
{
	assert(delay_level > 0);
	assert(m_nDst == 0 && m_sum_in == 0);
	LEVEL_NODE* level_node = m_cluster->Level[delay_level];
	NODE *node = 0;
	NUM_ELEMENTS node_index = 0;
	NUM_ELEMENTS nNodes = 0;


	nNodes = level_node->get_nNodes();

	debugif(DSIMPLE_EDGE, "Adding " << nNodes << " nodes to the destination list\n");
	
	for (node_index = 0; node_index < nNodes; node_index++)
	{
		node = level_node->get_node(node_index);
		assert(node);
		m_dst[m_nDst++] = node;
	}
	m_sum_in += nNodes * m_cluster->get_kin();
}


//
//  Work first with the src list, to allocate edges from high fanout 
//  nodes evenly across the dst list.
//
//
//  PRE: m_src and m_dst are the source and destination lists and 
//       and have been assigned nodes
//
//  POST: we have assigned intra-cluster edges between the high fanout source
//        nodes and the destination nodes if we could
//
void SIMPLE_EDGE_ASSIGNER::do_high_fanouts()
{
    NODE* source_node = 0;
	LOCALITY i;
	LOCALITY range = 0;

    debugif(DSIMPLE_EDGE,"Allocating edges for high-fanout nodes in source list");

    for (i = 0; i < m_nSrc; i++) 
	{
		source_node = m_src[i];

		debugif(DSIMPLE_EDGE,"Src-node " << i << " has to assign fanout of " <<  
				source_node->get_edges_currently_assigning());

		// if we have a high fanout degree
		if (m_nDst > m_nSrc && source_node->get_edges_currently_assigning() > m_nDst/4) 
		{
			debugif(DSIMPLE_EDGE,"..Is high-degree---assigning now");

			range = 5*source_node->get_fanout_degree()/4;

			do_one_high_fanout(i, range);

			if (source_node->get_edges_currently_assigning() > m_nDst/8) 
			{
				debugif(DSIMPLE_EDGE,"..Second attempt at high-degree on this node");

				range = 3*source_node->get_fanout_degree()/2;
				do_one_high_fanout(i, range);
			}
		}
    }
    debugif(DSIMPLE_EDGE,"Done with allocating for high fanout," << get_status_string() << endl);

}


// 
//  Allocate edges from a single high-degree node source_node on the src list
//
//  Choose a band of nodes of from the destination list and try to assign edges
//  between the source node and the nodes in the destination list
// 
//  POST: we have assigned edges between a high fanout source node and the 
//        nodes in the destination list if we could
//
void SIMPLE_EDGE_ASSIGNER::do_one_high_fanout
(
	const LOCALITY& high_fanout_srcindex, 
	LOCALITY& range
)
{
    NODE *source_node = 0, 
		*dst_node = 0;
    int start, i, edges_currently_assigning;
    NUM_ELEMENTS_VECTOR perm;
	NUM_ELEMENTS kin = m_cluster->get_kin();

    debugif(HIGHDEG, "..Doing assignment for " << high_fanout_srcindex << 
			", with range " << range << " of " << m_nDst);
		
    source_node = m_src[high_fanout_srcindex];
    edges_currently_assigning = source_node->get_edges_currently_assigning();

	// find out where we should look for a node to connect to in the destination list
	for_the_destination_list_find_a_start_position_and_range(high_fanout_srcindex, range, start);
    assert(edges_currently_assigning <= range);

	// can a random permutation of numbers that we will use when we walk
	// randomly through the destination list
	perm = g_rand_number_gen->get_permutation(range);

	// make connections between the source node and destination nodes
    for (i = 0; i < edges_currently_assigning; i++) 
	{
		dst_node = m_dst[start + perm[i]];
		assert(dst_node->get_fanin_degree() < kin);

		if (dst_node->are_nodes_connected(source_node))
		{
			debugif(HIGHDEG, "......Already connected to " << start + perm[i] << ", continuing"); 
		} 
		else 
		{
			debugif(HIGHDEG, "......Connecting src[" << high_fanout_srcindex << "] " << 
					"to dst[" << start+perm[i] << "] port " << 1 + dst_node->get_fanin_degree());

			assert(dst_node->get_fanin_degree() < kin);
			do_connect(high_fanout_srcindex, start + perm[i]);
		}
    }

	perm.clear();

    remove_dst_nodes_with_full_inputs(start);
}


// PRE: high_fanout_srcindex is the index of the source node
//      range has been assigned a value
// POST: start and range have values consistent with the size of m_dst
void SIMPLE_EDGE_ASSIGNER::for_the_destination_list_find_a_start_position_and_range
(
	const LOCALITY& high_fanout_srcindex, 
	LOCALITY& range,
	int& start
)
{
	int mid = 0,
		end = 0;
	//
	// Find the correct range and mid coordinates
	// 
    if (range > m_nDst) 
	{
		range = m_nDst;
		start = 0; 
		mid = m_nDst / 2;
		end = m_nDst;
    } 
	else 
	{
#if NEW
		// Find closest index to the source index
		mid = 0;
		best = abs(m_dst[0]->get_index() - source_node->get_index());
		for (i = 0; i < m_nDst; i++) 
		{
		   curr = abs(m_dst[i]->get_index() - source_node->get_index());
		   if (curr < best) 
		   {
				best = curr;
				mid = i;
		   }
		}
#else
        mid = (int) nint((double) m_nDst * high_fanout_srcindex / m_nSrc);
#endif
        debugif(HIGHDEG, "....Node for " << high_fanout_srcindex << " of " << m_nSrc << 
				" in src[] is " << mid << " of " << m_nDst << " in dst[]");

		start = mid - (range % 2 + range/2);
		end   = mid + range/2;

		debugif(HIGHDEG, "......initial start=" << start << ", end=" << end);

		if (start < 0) 
		{
			end += (-start);
			start = 0;
		}
		if (end > m_nDst) 
		{
			start -= (end - m_nDst);
			end = m_nDst;
		}
		start = MAX(start,0);
		end = MIN(end,m_nDst);
		range = end - start;
		debugif(HIGHDEG, "......final start=" << start << ", end=" << end << " diff=" << end-start <<
				" for range " << range);
		sanity(end <= m_nDst && start + range <= m_nDst && start >= 0);
    }  
}

//
//  Remove destination nodes that can no longer accept an input
//
//  POST: no node in m_dst has inputs >= kin
//
void SIMPLE_EDGE_ASSIGNER::remove_dst_nodes_with_full_inputs(const int & start)
{
    int i, removed;
	NUM_ELEMENTS kin = m_cluster->get_kin();

    debugif(DSIMPLE_EDGE,"..Removing full nodes from dstlist.  " << m_sum_in << 
			" inputs/" << m_nDst << " nodes at start");

    removed = 0;
    i = start; 

    while (i < m_nDst) 
	{
		debugif(HIGHDEG && FULL, "....i=" << i << ", removed=" << removed << ", nDst=" << m_nDst);

		m_dst[i] = m_dst[i+removed];
		assert(m_dst[i]->get_fanin_degree() <= kin);

		// if our node has its inputs full remove it from the list
		if ((m_dst[i]->get_fanin_degree()  == kin) || 
		   (m_dst[i]->is_DFF() && m_dst[i]->get_fanin_degree() == 1))
		{
			debugif(DSIMPLE_EDGE,"....removing " << i + removed);
			removed += 1;
			m_nDst -= 1;
		}
		else 
		{
			i += 1;
		}
    }
    debugif(DSIMPLE_EDGE,"..Done removing full nodes.  removed " << removed << " leaving " 
			<< m_sum_in << " inputs/" << m_nDst << " nodes" << endl);
}


//
//  Take a linear pass through the destination list to make sure that each node
//  gets assigned at least one node from the previous comb delay level.
//
//  PRE:  
//  POST: each node should be connected to at least one length 1 edge
//        if this couldn't be done, then we have warned the user
//
void SIMPLE_EDGE_ASSIGNER::establish_delay_level(const DELAY_TYPE& delay_level)
{
	assert(delay_level > 0);
    NODE *dst_node = 0;
    bool ok = false;
	NUM_ELEMENTS i;

	assert(new_gen);

    debugif(DSIMPLE_EDGE,"Establishing combinational delay level for dst nodes");

    for (i = 0; i < m_nDst; i++) 
	{
		dst_node = m_dst[i];

		// if this nodes does not have an edge, try to give it one
		if (dst_node->get_fanin_degree() == 0)
		{
			debugif(DSIMPLE_EDGE,"..Choosing fanins with span " << dst_node->get_span() << 
					" (nin=" << dst_node->get_fanin_degree() << ")");

			ok = choose_and_make_two_connections(i);
			if (!ok) 
			{
				debugif(DSIMPLE_EDGE,"....Didn't take -- trying again for just one");
				ok = choose_and_make_a_single_connection(i);
			}
			if (!ok) 
			{
				Log("Unimportant sanity violation at establish delay in simple_edge_assignment.C");
			}
		}
		else 
		{
			assert(dst_node->get_fanin_degree() > 0);
			debugif(DSIMPLE_EDGE,"..Node dst[" << i << "]==" << dst_node->get_id() 
					<< ", lev " << delay_level << " is OK");
		} 
    }

    debugif(DSIMPLE_EDGE,"Done with allocating to establish delay level, " << get_status_string() << endl);
		
	// remove any destination nodes that now have their inputs full
	remove_dst_nodes_with_full_inputs(0);
}

//
//  Take another pass through the dstlist to make sure that each node
//  gets assigned at least two in-edges.  
//
//  This is almost exactly the same step as ensure_delay(), except
//  for counting setting the tolerance to 2.
//  One other thing that is important is that it is now possible for
//  the initial guess to already be conneted to the node.  So we look
//  one node to the left or right if this has happened.
//
//
//  POST: each node should be connected to at least two edges
//        if this couldn't be done, then we have warned the user
//
void SIMPLE_EDGE_ASSIGNER::ensure_no_buffer_nodes()
{
    NODE *dst_node = 0;
    int i, prev_nDst;
    bool ok;

    debugif(DSIMPLE_EDGE,"Allocating edges to endure no buffers," << get_status_string());

    for (i = 0; m_nSrc>0 && m_sum_out>0 && m_nDst>0 && m_sum_in>0 && i<m_nDst; i++) 
	{
		prev_nDst = m_nDst;
		dst_node = m_dst[i];
		debugif(DSIMPLE_EDGE,"working on node " << i << ", " << get_status_string());

		if (dst_node->get_fanin_degree() >= 2) 
		{
			debugif(DSIMPLE_EDGE,"..Node " << i << " is OK");
		} 
		else if (m_nSrc == 1 && dst_node->are_nodes_connected(m_src[0]))
		{
			Warning("We have a buffer node");
		} 
		else 
		{
			debugif(DSIMPLE_EDGE,"Choosing a fanin for node " << i);
			ok = choose_and_make_a_single_connection(i);
			sanity(ok);
			if (m_nDst < prev_nDst) 
			{
				debugif(DSIMPLE_EDGE,"Compensating for decrease in nDst by decrementing i");
				i -= 1;
			}
		}
    }

    if (m_nDst == 0 || m_nSrc == 0 || m_sum_out == 0) 
	{
		// Happens because of the "cheap hack" above
		while (m_nSrc) 
		{
			// cheap hack needs to be removed
			assert(false);
			--m_nSrc;
			m_src[m_nSrc]->set_edges_to_assign(0);
			m_src[m_nSrc]->set_edges_currently_assigning(0);
		}
		m_nDst = m_nSrc = m_sum_in = m_sum_out = 0;
    }
    debugif(DSIMPLE_EDGE,"Done removing possibility of buffers, " << get_status_string() << endl);
}


//
//  Now continue to assign edges until we run out:  Choose a random dst 
//  node, then choose a random fanin node for it.  Continue until we 
//  run out of src nodes.
//  
//  POST: we can no longer assign an edge between the ndoes in the source list and 
//        nodes in the destination list
//
void SIMPLE_EDGE_ASSIGNER::assign_remaining_edges()
{
    int iDst;
    bool ok;
    debugif(DSIMPLE_EDGE,"Allocating edges for remaining fanin nodes in dst list");

    while (m_sum_out > 0 && m_nDst > 0 && m_nSrc > 0)
	{
		iDst = random() % m_nDst;
		debugif(DSIMPLE_EDGE, get_status_string());
		debugif(DSIMPLE_EDGE, "Looking for fanin for dst[" << iDst << "] " <<
				m_dst[iDst]->get_name() <<" (has " << m_dst[iDst]->get_fanin_degree() << " inputs now)");

		sanity(m_dst[iDst]->get_fanin_degree() < m_cluster->get_kin());
		ok = choose_and_make_a_single_connection(iDst);
		if (!ok) 
		{
			debugif(DSIMPLE_EDGE,"Lost node " << iDst << " from dst list because of fanin problems");
			m_dst[iDst] = m_dst[--m_nDst];
		}
		if (DEBUG && m_sum_out % 10 == 0) 
		{
			show_progress();
		}
    }

    debugif(DSIMPLE_EDGE,"Done with assigning remaining edges, " << get_status_string() << endl);
}



//
//  It is possible, at this point, that there are unresolved edges
//  because (a) the nodes on dst and src are already connected, or
//  (b) because we ran out of fanin ports before we ran out of out
//  edges.  
//  
//  PRE:  the only nodes left of the source list are the nodes that 
//        could not connect to anything in the destination list
//  POST: nodes in the source list have had edges_currently_assigning reset to zero
//        the source nodes have been added to the lost node list
//        the source and destination list have been cleared
//		  
void SIMPLE_EDGE_ASSIGNER::no_unwanted_outputs(const DELAY_TYPE& delay)
{
    NODE *node = 0;
    int i, 
        num_lost = 0;
	NUM_ELEMENTS edges_currently_assinging = 0;

    debugif(DSIMPLE_EDGE,"Checking for unwanted outputs at level " << delay << 
			" of " << m_cluster->get_delay());

	// if we have stranded source nodes warn the user and 
	if (m_nSrc > 0 && g_options->is_verbose())
	{
		debug("\n****Stranded " << m_nSrc << " source nodes and " << m_sum_out <<
				" out-edges at end of edge-assignment");
	}

    for (i = 0; i < m_nSrc; i++) 
	{
		node = m_src[i];

		debugif(DSIMPLE_EDGE,node->get_name() << " with edges curr assigning " << 
				node->get_edges_currently_assigning());

		edges_currently_assinging = node->get_edges_currently_assigning();
		assert(edges_currently_assinging > 0);

		node->add_lost_simple_edges(edges_currently_assinging);
		num_lost += edges_currently_assinging;

		m_circuit->add_lost_source_node(node);

		node->add_edges_currently_assigning(- edges_currently_assinging);
		assert(node->get_edges_currently_assigning() == 0);

		debugif(DSIMPLE_EDGE,node->get_name() << " with edges curr assigning " << 
				edges_currently_assinging << " to_assigns, " << 
				node->get_unassigned_long_and_inter_cluster_edges() << " long of " << node->get_fanout_degree());

		if (node->get_nOutputs() > 0)
		{
			debugif(DSIMPLE_EDGE,"....is OK, has at least one output");
		} 
		else if (node->get_unassigned_long_and_inter_cluster_edges() > 0) 
		{
			debugif(DSIMPLE_EDGE,"....could get lucky, it still has long edges");
		} 
		else 
		{
			debugif(DSIMPLE_EDGE,"...This node will not be connected to anything.");
		}
	
		// need to fix below
		assert(new_gen);
		//if (delay - node->get_delay_level() > 1)
		//{
			// we have a source node that was added to this level as a long node.
			// give it back it's long edges
		//}

		assert(node->get_unassigned_long_and_inter_cluster_edges() + 
			   node->get_edges_to_assign() <= node->get_fanout_degree());
    }

	assert(new_gen);
    if (num_lost > 0) 
	{
        Verbose("Lost " << num_lost << " edges at delay level " << delay);
    }

    m_nSrc = m_nDst = 0;
    m_sum_out = m_sum_in = 0;

    debugif(DSIMPLE_EDGE,"Done with checking for unwanted outputs. src[], dst[] zeroed");
}


//  Choose a single fanin for the specified node.
//
//  PRE: iDst is the position in the destination list of the destination node 
//  POST: we have made a connections if we could
//  RETURN: true if we could make a connection, false otherwise
//
bool SIMPLE_EDGE_ASSIGNER::choose_and_make_a_single_connection(const int& iDst) 
{
    int j, best, bestcost, choice, cost, i;
    NODE *source_node = 0, *dst_node = 0;
    int useloc;
	DISTRIBUTION permutation;

    dst_node = m_dst[iDst];
	assert(dst_node);

    assert(dst_node->get_fanin_degree() < m_cluster->get_kin());

    if (m_nSrc == 0) 
	{
		Verbose("Ran out of sources in edge allocation");
		return false;
    }

	i = 0;
	permutation = g_rand_number_gen->get_permutation(m_nSrc);
	best = permutation[i++];

    debugif(DSIMPLE_EDGE,"..Initial choice for fanin is " << best << " of " <<  m_nSrc);
    source_node = m_src[best];

	while (dst_node->are_nodes_connected(source_node) && i < m_nSrc)
	{
		best = permutation[i];
    	source_node = m_src[best];
        debugif(DSIMPLE_EDGE,"....already connected, taking neighbour " << best << " instead");
		i++;
	}

    if (dst_node->are_nodes_connected(source_node))
	{
		debugif(DSIMPLE_EDGE,"Couldn't find a starting point for node " << iDst << " " << 
				dst_node->get_name() << ", giving up");
		return false;
    }

	// we found a valid connection.
	// get the cost of the connection
    bestcost = abs(source_node->get_index() - dst_node->get_index());

    //  Some special cases that arise are when the number of sources
    //  we are choosing is much greater than the number of dst, or 
    //  if locality (calculated from n) was underestimated because the
    //  circuit is heavy on one level.  To get better behaviour, we bump 
    //  up the locality in either case.  Note that this is not the same
    //  as just globally increasing locality, as we don't do it when 
    //  the levels are balanced.


    useloc = m_cluster->get_comb_spec()->get_mike_hutton_locality();
    if (m_nSrc > m_nDst) 
	{
		useloc = useloc * m_nSrc / m_nDst;
    }
    if (m_nSrc > m_cluster->get_nNodes()/ 4) 
	{
		useloc = useloc * 2;
    }


	// now see if we can find a better cost

    //  New locality strategy:  Take the *first* valid connection between
    //  lspan and rspan.  If we never get such a node, take the best 
    //  after L tries (as before).

    debugif(DSIMPLE_EDGE,"..Need fanin for " << iDst << "; guess " << best << ", cost " << bestcost);
    for (j = 0; bestcost > 0 && j < useloc; j++) 
	{
        choice = random() % m_nSrc;
        cost = abs(m_src[choice]->get_index() - dst_node->get_index());
		if (dst_node->are_nodes_connected(m_src[choice]))
		{
			cost = INT_MAX;
		} 
		else if (cost < bestcost) 
		{
			best = choice;
			source_node = m_src[best];
			bestcost = cost;

			// the the src node is within our span. the cost of the node is 0
			if (m_src[choice]->get_index() > dst_node->get_lspan() && 
				m_src[choice]->get_index() < dst_node->get_rspan() ) 
			{
				cost = 0;
			}
        }
    }

    debugif(DSIMPLE_EDGE,"Will connect from src " << best << " with cost " << bestcost);
    do_connect(best, iDst);
    return true;
}


//
//  Choose two fanins at distance span apart for node.
//
//  PRE: iDst is the position in the destination list of the destination node 
//  POST: we have made to connections if we could, otherwise no connections have been made
//  RETURN: true if two edges could be assigned, false otherwise
//  
bool SIMPLE_EDGE_ASSIGNER::choose_and_make_two_connections(const int & iDst) 
{
    int i, j, iTmp, start, choice, nchoices;
    NODE *dst_node = 0;
	NUM_ELEMENTS_VECTOR::iterator new_end;
	NUM_ELEMENTS_VECTOR pmf;

    dst_node = m_dst[iDst];
	assert(dst_node);

    if (m_cluster->get_kin() - (dst_node->get_fanin_degree()) < 2) 
	{
		// we can only make one connection to the destination node.
		// we cannot make two connections
		debugif(DSIMPLE_EDGE,"....only one node to assign.  returning");
		return false;
    }


    //  This is the number of times we will sample from the source list
    nchoices = m_cluster->get_comb_spec()->get_mike_hutton_locality();
    nchoices = MAX(nchoices, m_nSrc/10);
    nchoices = MAX(nchoices, 10);
    nchoices = MIN(nchoices, CHOICELIST_SIZE-1);
    nchoices = MIN(nchoices, m_nSrc);
	assert(nchoices >= 0 && nchoices < CHOICELIST_SIZE);

    debugif(LOC,"\nWill make " << nchoices << " fanin choices of " << m_nSrc << 
			" for node at index " << dst_node->get_index() << ", span " << dst_node->get_span())

    //  Build a list of choices.

    debugif(LOC, "Choosing...");


	// build a choice list from the source list
	pmf.resize(m_nSrc, 1);
	i = 0;
	while (i < nchoices)
	{
		choice = random() % m_nSrc;

		if (pmf[choice] == 1)
		{
			debugif(LOC, "..chose " << choice << " at index " << m_src[choice]->get_index() << " "
					<< m_src[choice]->get_name());
			m_choice_list[i] = choice;
			pmf[choice] = 0;
			i++;
		}
    }

    //  Prune the list of things that are no good.
	//  ie. whose span is greater than the difference of the indici
    debugif(LOC, "Pruning...");
    for (i = 0; i < nchoices; i++) 
	{
		if ( abs(m_src[m_choice_list[i]]->get_index() - dst_node->get_index()) > dst_node->get_span()) 
		{
			debugif(LOC,"..deleting choice " << m_choice_list[i] << " with index " << 
					m_src[m_choice_list[i]]->get_index() << " -- " << m_src[m_choice_list[i]]->get_name());
			m_choice_list[i] = m_choice_list[nchoices - 1];
			nchoices -= 1;
			i -= 1;
		}
    }

    if (nchoices <= 1) 
	{
		debugif(DSIMPLE_EDGE,"....only one choice available.  returning");
		return false;
    }


	// note: i almost never get here because the span it too small of the width
	// maybe something to reconsider
	assert(new_gen);

    debugif(LOC, "Sorting and Uniquing the node list...");
	sort(m_choice_list.begin(), m_choice_list.begin() + nchoices);
	new_end = unique(m_choice_list.begin(), m_choice_list.begin() + nchoices);
	nchoices = distance(m_choice_list.begin(), new_end);

    //  Sort the list based on index.  Simple bubblesort.
    debugif(LOC, "Sorting...");
    for (i = 0; i < nchoices - 1; i++) 
	{
        for (j = 0; j < nchoices - 1; j++) 
		{
			if (m_src[m_choice_list[j]]->get_index() > m_src[m_choice_list[j+1]]->get_index()) 
			{
				iTmp = m_choice_list[j];
				m_choice_list[j] = m_choice_list[j+1];
				m_choice_list[j+1] = iTmp;
			}
		}
    }

#if LOC 
    debugif(DSIMPLE_EDGE,"LOC list for node at " << dst_node->get_index() << " span " << dst_node->get_span() << ": ");
	dbg("Indici: ");
    for (i = 0; i < nchoices; i++) 
	{
		dbg(m_src[m_choice_list[i]]->get_index() << " "); 
    }
	dbg("\nnames: ");
    for (i = 0; i < nchoices; i++) 
	{
		dbg(m_src[m_choice_list[i]]->get_name() << " "); 
    }
    dbg("\n");
#endif

	// narrow down our choices until our start and end indicies are within span.
	assert(new_gen);
    start = 0;
    while ((m_src[m_choice_list[nchoices-1]]->get_index() - m_src[m_choice_list[start]]->get_index() 
		> dst_node->get_span()) && (nchoices - start > 2 )) 
	{
		choice = random() % 2;
		if (choice == 0) 
		{
			start += 1;
		} 
		else 
		{
			nchoices -= 1;
		}
    }

    if (start == nchoices - 1) 
	{
		return false;
    }


	// set the new span for the destination node
    dst_node->set_lspan(m_src[m_choice_list[start]]->get_index());
    dst_node->set_rspan(m_src[m_choice_list[nchoices - 1]]->get_index());
    debugif(DSIMPLE_EDGE,"Loc-connecting indici " << m_src[m_choice_list[start]]->get_index() << 
			" and " << m_src[m_choice_list[nchoices - 1]]->get_index() <<
			" to " << dst_node->get_index() << " with span " << dst_node->get_span());

	assert(m_choice_list[nchoices - 1] != m_choice_list[start]);


	// connect the destination node with the two end point of the choices

    // MUST do the larger number first, because the list changes
    do_connect(m_choice_list[nchoices - 1], iDst);
    do_connect(m_choice_list[start], iDst);
    dst_node->set_span(dst_node->get_rspan() - dst_node->get_lspan() + 1);

    debugif(DSIMPLE_EDGE, sep);
	
    return true;
}


//
//  Make an connection and update counters and lists.
//
//  PRE: iSrc in the index in the source list
//       iDst is the index in the destination list
//  POST: an simple edge has been created between two nodes
//  	  counters and list have been updated
//
void SIMPLE_EDGE_ASSIGNER::do_connect(const int& iSrc, const int& iDst)
{
	NODE * dst_node = m_dst[iDst];
	NODE * source_node = m_src[iSrc];

	assert(dst_node && source_node);

	LEVEL_NODE * source_level_node = source_node->get_my_level_node();
	LEVEL_NODE * sink_level_node = dst_node->get_my_level_node();
	EDGE * super_edge = sink_level_node->find_edge_with_source(source_level_node);
	EDGE * edge = 0;

    sanity(dst_node->get_fanin_degree() < m_cluster->get_kin());

    debugif(DSIMPLE_EDGE,"Connecting source " << source_node->get_name() << " at index " 
			<< iSrc << " with dst " << dst_node->get_name() << " at index " << iDst);

    edge = m_circuit->create_edge(source_node, dst_node);
	edge->set_my_super_edge(super_edge);

    source_node->add_edges_currently_assigning(-1);

	if (dst_node->is_DFF())
	{
    	debugif(DSIMPLE_EDGE,"designating latch " << dst_node->get_name() << " as latched");

		source_node->set_is_latch_taken(true);
	}
    debugif(DSIMPLE_EDGE,"....source node " << source_node->get_name() << " now has "
			<< source_node->get_edges_currently_assigning() << " edges currently assigning");

    m_sum_out -= 1;
    if (source_node->get_edges_currently_assigning() == 0) 
	{
		m_nSrc -= 1;
        m_src[iSrc] = m_src[m_nSrc];
        debugif(DSIMPLE_EDGE,"..Removing exhausted fanout " << source_node->get_name() << 
				" from source list (" << m_nSrc << " nodes left)");
    }

    dst_node->add_fanin_degree(1);
    m_sum_in -= 1;
    if (dst_node->get_fanin_degree() == m_cluster->get_kin() || dst_node->is_DFF())
	{
		m_dst[iDst] = m_dst[--m_nDst];
		debugif(DSIMPLE_EDGE,"..Removing full fanin " << iDst << " from dst[] (" << m_nDst << " left)");
    }
    debugif(FULL, "Done with do_connect for " << iSrc << " and " << iDst);
}



//
//  Debug routine to show progress on a level.
//
void SIMPLE_EDGE_ASSIGNER::show_progress()
{
#if DEBUG
    int i, check;

    dbg("\n");
    debugsepif(DSIMPLE_EDGE);
    debugif(DSIMPLE_EDGE,"Current status:");
    dbg("src (" << m_sum_out << "/" << m_nSrc << "):  ");
    check = 0;
    for (i = 0; i < m_nSrc; i++) 
	{
		dbg(m_src[i]->to_assign << " ");
		check += m_src[i]->to_assign;
    }
    dbg("\n");
    sanity(check == m_sum_out);
    dbg("dst (" << sum_in << "/" << m_nDst << "):  ");
    check = 0;
    for (i = 0; i < m_nDst; i++) 
	{
		dbg(m_cluster->kin - (dst[i]->nin + dst[i]->nGI) << " ");
		check += m_cluster->kin - (dst[i]->nin + dst[i]->nGI);
    }
    dbg("\n");
	sanity(check == sum_in);
    debugsepif(DSIMPLE_EDGE);
    dbg("\n");
#endif
}



//
// RETURN: returns a status string for the edge assignment
string SIMPLE_EDGE_ASSIGNER::get_status_string()
{
	string status = "src=" + util_long_to_string(m_sum_out) + " out/" + util_long_to_string(m_nSrc) + 
		" nodes, dst=" + util_long_to_string(m_sum_in) + " in/" + util_long_to_string(m_nDst) + " nodes";
	return status;
}

//
// POST: the latched nodes have been added to the source list
//
void SIMPLE_EDGE_ASSIGNER::add_latched_nodes_to_source_list()
{
	LEVEL_NODE 	* sink_level_node = 0,
				* source_level_node = 0;
	EDGES input_edges;
	NODES source_nodes;
	EDGES::iterator edge_iter;
	EDGE * edge = 0;
	NODE * source_node = 0;
	NUM_ELEMENTS edge_weight = 0;

	sink_level_node = m_cluster->Level[0];
	assert(sink_level_node);

	input_edges = sink_level_node->get_input_to_dff_edges();

	for (edge_iter = input_edges.begin(); edge_iter != input_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		edge_weight = edge->get_weight();

		source_level_node = edge->get_source_level_node();
		assert(source_level_node);
		
		assert(new_gen);
		// we might want to make to get more latched nodes than
		// the number of dff so that we get good locality assignment.
		// but what we have to watch for is that we only take
		// the number of latched nodes that we are alloted from any source 
		// level node.
		//


		// we need to get the right latch nodes
		source_nodes = source_level_node->get_latched_nodes_not_taken(edge_weight);

		assert(source_nodes.size() >= static_cast<unsigned>(edge_weight));

		while (! source_nodes.empty())
		{
			source_node = source_nodes.back();
			source_nodes.pop_back();
			assert(source_node);
			
			assert(source_node->get_edges_currently_assigning() == 0);
			assert(source_node->is_latched() && ! source_node->is_latch_taken());
			
			// node is not on the list
			// add it to the list
			
			add_source_node(source_node,1);
		}
	}
}

// POST: flip-flops have been added to the destination list
void SIMPLE_EDGE_ASSIGNER::add_dff_nodes_to_dst_list()
{
	assert(m_nDst == 0 && m_sum_in == 0);
	LEVEL_NODE* level_node = m_cluster->Level[0];
	
	m_dst = level_node->get_DFF();
	m_nDst = level_node->get_nDFF();

	m_sum_in += level_node->get_nDFF();
}

// POST: we have made all the latched node to flip-flop connections 
void SIMPLE_EDGE_ASSIGNER::assign_dff_simple_edges()
{
    int iDst;
    debugif(DSIMPLE_EDGE,"Allocating edges for dff connections\n");
	bool success = true;

	// we should do the intra cluster dff edges first

	while (m_nDst > 0)
	{
		// assign the last element on the list
		iDst = m_nDst-1;
		// only one connection for a dff
		assert(m_dst[iDst]->get_fanin_degree() < 1);
		success = choose_and_make_a_single_connection(iDst);
		assert(success);
    }

    debugif(DSIMPLE_EDGE,"Done with assigning remaining edges, " << get_status_string());
}

// there will be some latched not that are not used because the source nodes will
// always be greater than or equal to the number of dst nodes for all clusters except the last
//
// POST: we have removed the latched nodes from the source list 
void SIMPLE_EDGE_ASSIGNER::remove_unused_latched_nodes()
{
	NODE * node = 0;
	NUM_ELEMENTS i;


	for (i = 0; i < m_nSrc; i++)
	{
		node = m_src[i];
		assert(node);

		node->add_edges_currently_assigning(-1);
		m_sum_out--;

		m_src.pop_back();

		Fail("Didn't use node " << node->get_name());
	}
}

// output a debug message
//
void SIMPLE_EDGE_ASSIGNER::output_debug_msg() const
{
	debugif(DSIMPLE_EDGE,"Have " << m_nSrc << " nodes with " << m_sum_out << " edges out; " 
			<< m_nDst << " nodes with " << m_sum_in << " possible edges in.\n");
}

// warn the user about stranded nodes or differences between the 
// delay structure with super edges and the delay structure with individual edges
void SIMPLE_EDGE_ASSIGNER::output_warning_msgs_to_user
(
	const LEVEL_NODE * level_node
) const
{
	assert(level_node);

	// if we have not assigned all the edges we should have, warn the user
	NUM_ELEMENTS input_difference = level_node->get_nSimple_input_edges() - 
									level_node->get_input_edge_weight();
	if (input_difference != 0 && g_options->is_verbose())
	{
		debug("Level node " << level_node->get_name() << " had a simple edge input difference from phase A of " << input_difference);
	}
}
