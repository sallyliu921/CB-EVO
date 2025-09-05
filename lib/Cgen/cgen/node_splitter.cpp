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



#include "node_splitter.h"
#include <algorithm>
#include <numeric>
#include "fp.h"


NODE_SPLITTER::NODE_SPLITTER()
{
	m_cluster = 0;
}

NODE_SPLITTER::~NODE_SPLITTER()
{
	m_cluster = 0;
}

//  Mainline for the first part of node-splitting
//
//  The next step in the generation process is to split the 
//  level nodes into individual nodes (which will ultimately become 
//  4-input LUTs, primary inputs, and flip-flops in the final generated circuit). 
//  The graph is also prepared for final edge assignment. The input into this phase 
//  of the algorithm is the fanout-assigned delay graph structure 
//  with the latched shape, the primary output shape, and the number of 
//  primary inputs and flip-flops. The output from this is a graph where the 
//  individual nodes have been created at each level node and assigned 
//  a logic type from one of flip-flop, primary input, or LUT; 
//  where each node has an assigned fanout; where each
//	node has a horizontal position and where all latched nodes and 
//	nodes that are primary outputs have been designated. 
//	This structure is called the pre-edge assignment structure
//
// PRE: We have assigned the fanout distribution amongst the level nodes
// POST: We are reading for final edge assignment
//
NUM_ELEMENTS NODE_SPLITTER::split_nodes(CLUSTER * cluster)
{
	assert(cluster);
	m_cluster = cluster;
	DELAY_TYPE lev = 0;
	NUM_ELEMENTS quality = 0;

    debugif(DNODE_SPLIT, "Splitting level_nodes into simple nodes");
    debugif(DNODE_SPLIT,"Allocating all nodes now");
		
    debugif(DNODE_SPLIT,"About to assign fanout and long and inter cluster edges...");
    for (lev = 0; lev < m_cluster->get_delay(); lev++) 
	{
		if (m_cluster->Level[lev]->get_nNodes() > 0)
		{
			// Step 1: assign nodes a fanout
			assign_fanouts(lev);

			// Step 1.5: designate which nodes have long/inter-cluster edges
			designate_long_and_inter_cluster_edges(lev);

			// Step 2: designate which nodes are latched
			designate_latched_nodes(lev);
		}
    }
	// Step 2: designate which nodes are latched
	designate_latched_nodes(m_cluster->get_delay());

	// Step 3: designate logic types for nodes
	//         all nodes at levels > 1 are combinational by default
	designate_dff_and_pi_nodes();

	// Step 4: designate which nodes have primary ouputs
	designate_PO();

	// Step 5: assign each node a horizontal position
	assign_nodes_a_horizontal_position(); 

    debugif(DNODE_SPLIT, Sep);
    display();

	final_sanity_check();

	quality = report_solution_quality();

    debugif(DNODE_SPLIT,"All done with node splitting and fanout allocation to nodes" << Sep);

	m_cluster = 0;

	return quality;
}


//
//  Assign each node an index.
//  The index is the horizontal position
//
//  The locality index is a balanced division into the numbers 0..16*width
//  where width is the width of the maximum delay-level in the graph.
//
//  PRE: m_cluster is valid
//  POST: each node has an index, and left and right span, and a span
//
void NODE_SPLITTER::assign_nodes_a_horizontal_position()
{
    DELAY_TYPE delay_level = 0;
    int width;
    LEVEL_NODE *level_node = 0;
    NODE * node = 0;
	NUM_ELEMENTS nNodes = 0;
	NUM_ELEMENTS j, index;
	NUM_ELEMENTS central_offset = 0,
				 avg = 0;
	CIRCUIT_GRAPH * circuit = m_cluster->get_my_circuit();
	assert(circuit);

	assert(new_gen); // should width be comb graph specific or circuit specific
					// in the past it was comb graph specific

	width = circuit->get_width();

    for (delay_level = 0; delay_level<= m_cluster->get_delay(); delay_level++) 
	{
		level_node = m_cluster->Level[delay_level];
		assert(level_node);

		nNodes = level_node->get_nNodes();

		if (nNodes > 0)
		{
			// the 16 in the width is just a factor to spread the nodes out a little.
			// it was used in hutton's work
			//
	
			// now as to how nodes are assigned horizontal positions:
			//
			// nodes in the delay level are uniformly spread out over the level delay
			//
			//
			// eg.
			// width = 20
			// nNodes at a level 7
			//
			// 0   20/7   2*20/7   3*20/7 ...  6*20/7
			//
			// therefore to correct the bias we find move the centre to the middle of the width ie.
			//
			// take 0 + 6*20/7 /2 to find the middle 
			//
			avg = 16/2 * width * (nNodes - 1)/nNodes;
			central_offset = 16/2 * width - avg;
			assert(central_offset >= 0);
		}

		for (j = 0; j < nNodes; j++) {

			node = level_node->get_node(j);
			assert(node);

			index = (16 * width * j) / nNodes + central_offset;
			assert(0 <= index && index < 16 * width);

			node->set_index(index);
			node->set_lspan(index);
			node->set_rspan(index);
			node->set_span(1);
			debugif(DNODE_SPLIT, "....node " << node->get_id() << " lev=" << delay_level << " index=" << index);
		}
    }
}


//
//  Assign the individual fanouts to the level nodes
//
//  We assign horizontal positions according to Hutton's method, in which 
//  high fanout nodes are assigned positions that balance them across each 
//  level node so as to not skew wirelengthApprox. The only change we make 
//  to the algorithm is that we assign a slightly different set of horizontal 
//  positions. The input into the algorithm is the level nodes where the individual nodes 
//  within have been assigned a fanout value. The output from the algorithm is the level nodes 
//  where the individual nodes within have been assigned a horizontal position. 
//  We assign horizontal positions to the nodes in each level node in three steps. 
//
//  In Step 1, we create a list of horizontal positions. In Step 2, we create 
//  a binary tree probability mass function (pmf) that assigns a probability 
//  to each horizontal position. To calculate this distribution we first take 
//  a balanced binary tree, of size equal to the number of individual nodes, and 
//  label its interior nodes with the number of leaves in the interior node's 
//  subtree and label its leaves with the number one. Next, we traverse this 
//  tree in-order and record the labels into a vector. An example of a vector 
//  of length 15 is [1,2,1,4,1,2,1,8,1,2,1,4,1,2,1]. Last in this step, we 
//  create the pmf from this vector by normalizing each entry in the vector 
//  by the total sum over all entries. At this point, each entry in the pmf 
//  corresponds to an entry in the list of horizontal positions. In Step 3, 
//  we assign a horizontal position to each node by sampling the pmf for an 
//  index into the list of horizontal positions in order of nodes of highest 
//  fanout to nodes of lowest fanout. After each sample the entry in the pmf 
//  is zeroed and the vector renormalized in order to keep the sum of the 
//  probabilities in the pmf equal to one. As Hutton noted, this method will 
//  most likely assign the highest fanout to the centre, with the next two 
//  fanouts to the quartiles.
//
//  POST: each individual node has a fanout
//
void NODE_SPLITTER::assign_fanouts
(
	const DELAY_TYPE& lev
)
{
    debugif(DNODE_SPLIT,"Assigning out-degree values to the nodes at level " << lev);

    debugif(DNODE_SPLIT, "Assigning high-fanout nodes to level " << lev);
    assign_high_fanouts(lev);

    debugif(DNODE_SPLIT, "Assigning low-fanout nodes to level " << lev);
    assign_low_fanouts(lev);

	display_fanout_results(lev);
}


//
//  Create a binary tree pmf.  This is a pmf which models an array-stored
//  binary tree with weights multiplying by k as we go up the tree.
//    e.g.   1  2  1  4  1  2  1  2   for length 8
//    e.g.   1  2  1  4  1  2  1  8  1  2  1  4  1  2  1  2   for length 16
//  are binary-tree distn's of weight 2.  For weight k we mult by k
//  rather than 2 as the height of the tree increases.  
//
//  Problem that this pmf will need numbers of size log2(n)*log2(mult) bits.  We
//  need to, thus, deflate this value to < 32 bits, hence the deflation 
//  calculation.  Also want to balance the array, so that the top of the
//  tree is at n/2, so calculate to next highest pwr of 2, and shift.
// 
//  PRE: n is the length of the array wanted
//       mult is the multiplying factor
//
DISTRIBUTION NODE_SPLITTER::tree_pmf(const int& n, int mult)
{
    int i, j, inc, sum;
    int deflate, len, shift;
    double deflate_val;
    double dval;
	unsigned int bits;
	double max_value = 0.0;
	double total = 0.0;
	bool should_output = false;

	DISTRIBUTION pmf;
	DOUBLE_VECTOR dpmf;

    debugif(should_output,"Calculating a tree-pmf for n=" << n << " mult=" << mult);

    deflate = 1;
    deflate_val = 1.0;
    bits = static_cast<int>(log2(static_cast<double>(n)) * log2(static_cast<double>(mult))) + 4;
    debugif(should_output,"Estimate I need " << bits << " bits for the eventual cdf");

    if (bits >= 8*sizeof(int)) 
	{
		deflate = bits - 8*sizeof(int);
		deflate_val = exp2(-deflate);
		debugif(should_output,"That isn't enough.  Will deflate all numbers by " << deflate_val << " (" << deflate << " bits)");
		assert(deflate_val > 0);
    } 
	else 
	{
		debugif(should_output,"OK -- numbers are small enough to not deflate");
    }

    len = (int) ceil(exp2(ceil(log2((double) n))));
    debugif(should_output,"Rounding n=" << n << " to next highest power of 2 len=" << len << " for balance");
	assert(len >= n);

    dpmf.resize(len, 0);
    pmf.resize(len, 0);

    shift = (len-n)/2;
    debugif(should_output,"Will shift the resulting array by " << shift);

    for (i = 0; i < len; i++) 
	{
		dpmf[i] = 1.0;
	}

    for (inc = 1; inc <= len/2 ; inc = inc * 2 ) 
	{
        for (i = inc; i < len; i+= inc) 
		{
			dpmf[i] = dpmf[i] * (double) mult;
			assert(dpmf[i] >= 0);
		}
    }
 
	for (i = 0; i < len; i++)
	{
		total += deflate_val*dpmf[i];
		max_value = MAX(max_value, dpmf[i]);
	}
	assert(total > 0);
	if (total > exp2(31.0))
	{
		deflate_val = deflate_val * exp2(31.0)/total;
		Verbose("Max value=" << max_value << ".  Tot=" << total << ".  This is too large!");
		Verbose("New deflate value= " << deflate_val);
	}

    assert(n + shift <= len);
    assert(shift >= 0);

    sum = 0;
    debugif(should_output,"pmf: ");
    for (j = 0, i = shift; i < n + shift; i++, j++) 
	{
		dval = static_cast<int>(dpmf[i] * deflate_val);
		assert(dval >= 0);

		dbgif(should_output, dval << "/");
		// usually sizeof(int) is 4 bytes. 
		// means the assert is
		// log2 of the real value is less than 31 bits
		assert(log2(dval) < (double) (8*sizeof(int)-1));

		pmf[j] = (int) nint(dval);
		sum += pmf[j];
		dbgif(should_output, pmf[j] << " ");
    }
    assert(sum > 0);
    dbgif(should_output, "\n");

	dpmf.clear();

    return pmf;
}


//
//  Assign high fanouts probabilistically to a binary-tree distribution. 
// 
//  PRE: lev is the delay level we want to assign fanouts to
//  POST: we have assigned all fanouts equal to or above 4 to indiv. nodes
//
void NODE_SPLITTER::assign_high_fanouts
(
	const DELAY_TYPE& lev
)
{

	// a problem we have is that sometimes we have pmf that is too large.
	// it forces some of the values to 0 which causes
	// us to not have any weight in the pmf past a certain point

	int i;
    NUM_ELEMENTS deg, choice;
    NUM_ELEMENTS num_assigned = 0;
    NODE * node = 0;
    LEVEL_NODE *level_node = m_cluster->Level[lev];
	assert(level_node);
	NUM_ELEMENTS nNodes = level_node->get_nNodes();
	DISTRIBUTION fanout_distribution = level_node->get_fanout_distribution();
	NUM_ELEMENTS total = 0;

    debugif(DNODE_SPLIT,"Allocating a binary-tree pmf for degree allocations");
    DISTRIBUTION pmf = tree_pmf(level_node->get_nNodes(), 8);
    DISTRIBUTION alternative_pmf(level_node->get_nNodes(), 0);

    debugif(DNODE_SPLIT,"Assigning out-degrees one-by-one");

    debugif(DNODE_SPLIT,"pmf:  ");

	assert(static_cast<unsigned>(nNodes) <= pmf.size());
    for (i = 0; i < nNodes; i++) 
	{
		dbgif(DNODE_SPLIT,pmf[i] << " ");
    }   dbgif(DNODE_SPLIT,"\n");


	// our altenative pmf has a value of 1 for 
	// where our tree_pmf has values of zero
	for (i=0; i < nNodes; i++)
	{
		if (pmf[i] == 0)
		{
			alternative_pmf[i] =1;
		}
	}


    deg = fanout_distribution.size() - 1;

	total = accumulate(pmf.begin(), pmf.end(), 0);

	// find out what is the highest existing fanout 
	// for this level node
    while (fanout_distribution[deg] == 0 && deg > 0)
        --deg;

	// assign all fanouts above 4
    while (deg >= 4 && num_assigned < nNodes)
	{
        debugif(DNODE_SPLIT, "..working on fanout " << deg);

		choice = g_rand_number_gen->discrete_pmf(nNodes, pmf);

        debugif(DNODE_SPLIT, "....got back choice of " << choice << " (from " << nNodes << ") val=" << pmf[choice]);
        debugif(DNODE_SPLIT, "..assigning fanout of " << deg << " to node " << choice);
		assert(pmf[choice] > 0);
		total -= pmf[choice];
		pmf[choice] = 0;

		// if we have no more weight with our tree 
		// use the alternative pmf
		if (total == 0)
		{ 
			pmf = alternative_pmf;
		}
		for (i = 0; i < nNodes; i++) 
		{
			dbgif(DNODE_SPLIT,pmf[i] << " ");
		}   dbgif(DNODE_SPLIT,"\n");
        assert(0 <= choice && choice < nNodes);

        node = level_node->get_node(choice);
		assert(node);
		
		num_assigned++;
        node->set_fanout_degree(deg);
		node->set_edges_to_assign(deg);

        level_node->OutDegrees(deg) -= 1;

        while (level_node->OutDegrees(deg) == 0 && deg > 0)
            --deg;
    }
}



//
//  PRE: lev is the delay level we want to assign fanouts to
//  POST: All fanout nodes below 4 have been assigned purely randomly, not based on binary tree.
//
void NODE_SPLITTER::assign_low_fanouts(const DELAY_TYPE& lev)
{
    int sum, sum1, deg;
    DISTRIBUTION pmf;
	NUM_ELEMENTS max_low_fanout, fanout;
    NODE *node = 0;
    LEVEL_NODE *level_node = 0;
	NUM_ELEMENTS node_number = 0;

    level_node = m_cluster->Level[lev];
	assert(level_node);

	max_low_fanout = MIN(level_node->get_max_fanout(), 3);
	assert(max_low_fanout >= 0 && max_low_fanout <= 3);
	

    dbgif(DNODE_SPLIT,"Assigning remaining low out-degrees randomly:  ");
	for (fanout = 0; fanout <= max_low_fanout; fanout++)
	{
		dbgif(DNODE_SPLIT,level_node->OutDegrees(fanout) << " ");
	}
	dbgif(DNODE_SPLIT,"\n");

    sum = 0;
    pmf = level_node->get_fanout_distribution();

    for (node_number = 0; node_number < level_node->get_nNodes(); node_number++)
	{
		node = level_node->get_node(node_number);
		assert(node);

		if (node->get_fanout_degree() == 0) 
		{
			// we haven't assigned fanout yet.

			sum1 = accumulate(pmf.begin(), pmf.end(), 0);
			if (sum1 > 0) 
			{
				deg = g_rand_number_gen->discrete_pmf(max_low_fanout+1, pmf);
				assert(deg >= 0 && deg <= max_low_fanout);
				assert(level_node->OutDegrees(deg) > 0);
				level_node->OutDegrees(deg) -= 1;
			} 
			else 
			{
				Fail("We ran out of out-degrees, that should not have happened.");
				// in the past mh assigned unit degrees
			}
			debugif(DNODE_SPLIT, "..assigning fanout " << deg << " to node " << node->get_id());
			node->set_fanout_degree(deg);
			node->set_edges_to_assign(deg);
			pmf[deg] -= 1;
		} 
		else 
		{
			debugif(DNODE_SPLIT, "..fanout " << node->get_fanout_degree() << " already assigned to node " << node->get_id());
		}
		sum += node->get_fanout_degree();
    }
    debugif(DNODE_SPLIT,"Have " << sum << " actual out-edges assigned of expected " << level_node->get_available_out());
    assert(sum == level_node->get_available_out());
}


//
//  PRE: lev is the delay level we want to display for 
//  POST: we have displayed  the fanout assignments for this lev
//
void NODE_SPLITTER::display_fanout_results
(
	const DELAY_TYPE& lev
) const
{
	LEVEL_NODE * level_node = 0;
	NODE * node = 0;
	NUM_ELEMENTS i;

	level_node = m_cluster->Level[lev];
	assert(level_node);

    debugif(DNODE_SPLIT,"Results for level " << lev << ": ");

    for (i = 0; i < level_node->get_nNodes(); i++) 
	{
		node = level_node->get_node(i);
		assert(node);
		dbgif(DNODE_SPLIT,node->get_fanout_degree() << " ");
    }
    dbgif(DNODE_SPLIT,"\n");
}

//
//  Move the appropriate number of to_assign edges to long_edges.
// 
//  PRE: delay level is valid
//  POST: the indiv. nodes that will have long/inter-cluster edges have been designated
//
void NODE_SPLITTER::designate_long_and_inter_cluster_edges(const DELAY_TYPE& delay_level)
{
    LEVEL_NODE *level_node = 0;
    NODE *node = 0;
    NUM_ELEMENTS nLong, longs_to_assign;
    NUM_ELEMENTS node_number = 0;
	NUM_ELEMENTS nNodes = 0;
	NUM_ELEMENTS choice = 0;
    NUM_ELEMENTS_VECTOR edges_left_to_assign_dist;
    NUM_ELEMENTS edges_left_to_assign;
	NODES nodes;

    level_node = m_cluster->Level[delay_level];
	assert(level_node);

    nLong = level_node->get_available_out() - level_node->get_intra_cluster_output_edge_lengths()[1];
	// note: we don't need to subtract 0 length latches because they are never added to the 
	// level node's available out

	nNodes = level_node->get_nNodes();

	assert(nLong >= 0);

    if (nLong <= 0) 
	{
		debugif(DNODE_SPLIT,"Level " << delay_level << " has no long edges to assign, returning");
		return;
    }
    debugif(DNODE_SPLIT, sep << "Need to assign " << nLong << " long edges from " << nNodes << " nodes, " <<
			level_node->get_available_out() << " edges");

    // First look to the nodes that will be forced to have longs
	// because of their high fanout 
    for (node_number =  0; node_number < nNodes; node_number++)
	{
		node = level_node->get_node(node_number);
		assert(node);

		// if the fanout degree of the node is greater than the number of 
		// nodes at the delay level below it will need long edges
		
		// note: we need to modify this to reflect the fact we now have inter cluster
		// edges which may be different than long intra-cluster edges
		assert(new_gen);

		longs_to_assign = node->get_fanout_degree() - m_cluster->Level[delay_level+1]->get_nNodes();
		if (longs_to_assign > 0)
		{
			// this node needs inter cluster or long edges
			debugif(DNODE_SPLIT,"MUST assign " << longs_to_assign << " long/inter-cluster to node with fanout " 
					<< node->get_fanout_degree());
			
			node->add_unassigned_long_and_inter_cluster_edges(longs_to_assign);
			node->add_edges_to_assign(-longs_to_assign);
			nLong -= longs_to_assign;
		}
    }
	assert(nLong >= 0);

    if (nLong == 0) 
	{
		debugif(DNODE_SPLIT,"All long edges assigned, returning");
		return;
    }

    // Randomly assign the rest of the long edges to the nodes
	// with nodes that have more edges to assign having a high probability 
	// of being assigned a long edge
    debugif(DNODE_SPLIT,"Choosing random " << nLong << " subset of " 
			<< level_node->get_available_out() << " edges");
    debugif(DNODE_SPLIT,"Allocating " << nLong << " remaining long edges from k-sample");

	// form a prob. mass function. weight = edges left to assign
    for (node_number =  0; node_number < nNodes; node_number++)
	{
		node = level_node->get_node(node_number);
		assert(node);

		edges_left_to_assign = node->get_edges_to_assign();
		edges_left_to_assign_dist.push_back(edges_left_to_assign);
    }


	// randomly pick the nodes that will have long edges 
	while (nLong > 0)
	{
		choice = g_rand_number_gen->discrete_pmf(nNodes, edges_left_to_assign_dist);
		assert(choice >= 0 && choice < nNodes);
		assert(edges_left_to_assign_dist[choice] > 0);

		node = level_node->get_node(choice);
		assert(node);
		assert(node->get_edges_to_assign() > 0);
		node->add_edges_to_assign(-1);
		node->add_unassigned_long_and_inter_cluster_edges(1);
		nLong--;

		edges_left_to_assign_dist[choice]--;
	}

    debugif(DNODE_SPLIT, sep);
}

//
//  Designate some number of the degree-0 nodes at level delay as latched.
//  latched means they output to the input of a flip-flop
//
//  PRE: delay_level is valid
//  POST: all latched nodes at this delay level have been designated
//
void NODE_SPLITTER::designate_latched_nodes
(
	const DELAY_TYPE& delay_level
)
{
    LEVEL_NODE *level_node = 0;
    NUM_ELEMENTS nNodes, nLatched, nZero_degee_latched_nodes;
    NODE *node = 0;
	NUM_ELEMENTS node_number = 0;
	NUM_ELEMENTS_VECTOR permutation;
	NUM_ELEMENTS max_latched_degree = 0;

    level_node = m_cluster->Level[delay_level];
	assert(level_node);

    nLatched = level_node->get_nLatched();
	nNodes = level_node->get_nNodes();

    if (nLatched == 0)
	return;


    debugif(DNODE_SPLIT,"Designating " << nLatched << " nodes at level " << delay_level << " to be latched");

   
	// first find out how many zero degree nodes exist
	nZero_degee_latched_nodes = 0;
    for (node_number = 0; node_number < nNodes; node_number++)
	{
		node = level_node->get_node(node_number);
		assert(node);

		if (node->get_fanout_degree() == 0) 
		{
			nZero_degee_latched_nodes++;
		}
    }
	
	// if we have more latched nodes that zero degree nodes warn the user
    if (nLatched - nZero_degee_latched_nodes > 0)
	{
		// should be warn or verbose
		assert(new_gen);
		Verbose("Warning: There are  " << nLatched - nZero_degee_latched_nodes << 
				" latched nodes of " << nLatched << " that have fanout greater than 0 at level " 
				<< delay_level);
	}


	// assign them randomly
	permutation = g_rand_number_gen->get_permutation(nNodes);

	max_latched_degree = 0;
	while (nLatched > 0)
	{
		for (node_number = 0; nLatched > 0 && node_number < nNodes; node_number++)
		{
			assert(permutation[node_number] >= 0 && permutation[node_number] < nNodes);

			node = level_node->get_node(permutation[node_number]);
			assert(node);

			if (node->get_fanout_degree() == max_latched_degree)
			{
				assert(! node->is_latched());
				node->set_is_latched(true);
				nLatched--;
			}
		}

		if (nLatched > 0)
		{
			max_latched_degree++;
			assert(max_latched_degree <= level_node->get_max_fanout());
		}
	}
 
    debugif(DNODE_SPLIT,"Done designating latches.  Max latched degree is " << max_latched_degree);
}


// Post-processing verbose display of the results.
//
void NODE_SPLITTER::display()
{
    LEVEL_NODE *level_node = 0;
    NODE *node = 0;
	DELAY_TYPE lev;
	NUM_ELEMENTS node_number = 0;

    debugif(DNODE_SPLIT, "\nResults of final degree assignment to nodes:\n");
    for (lev = 0; lev <= m_cluster->get_delay(); lev++) 
	{
		dbgif(DNODE_SPLIT,"= (" << lev << ") ");
		level_node = m_cluster->Level[lev];
		assert(level_node);

		for (node_number = 0; node_number < level_node->get_nNodes(); node_number++)
		{
			node = level_node->get_node(node_number);
			assert(node);

			dbgif(DNODE_SPLIT,node->get_fanout_degree() << " ");
			if (node_number+1 % 35 == 0 && node_number < level_node->get_nNodes())
			{
				dbgif(DNODE_SPLIT,"\n      ");
			}
		}
		dbgif(DNODE_SPLIT,"\n");
    }
    dbgif(DNODE_SPLIT,"\n");
}

// At the 0th delay level designate which nodes are primary inputs
// and which nodes are flip-flops
//
// Because the fanout characteristics of pi and dff are very different
// (pi usually have higher fanout with a wider standard deviation)
// we use the fanout statistics to help designate them correctly.
//
// PRE: fanouts have been assigned to the nodes
// POST: all nodes at the 0th delay level for the cluster in question
//       have been designated as either PI or flip-flop
//
void NODE_SPLITTER::designate_dff_and_pi_nodes()

{
	// we use a gaussian fanout distribution to model 
	// the pi and fanout statistics which is not 
	// quite correct most of the time
	//
	// work needs to be done to investigate other distributions
	//
	assert(new_gen);

	LEVEL_NODE * level_node = m_cluster->Level[0];
	assert(level_node);

	NODE * node = 0;
	NUM_ELEMENTS i;
	NUM_ELEMENTS nDFF = level_node->get_nDFF();
	NUM_ELEMENTS nNodes = level_node->get_nNodes();
	NUM_ELEMENTS nPI = m_cluster->get_nPI();
	NUM_ELEMENTS fanout_degree = 0;
	double pi_mean, dff_mean, pi_std_dev, dff_std_dev;
	double pi_probability, dff_probability, total_probability, random_number;
	double high_degree_boundry, low_degree_boundry;

	get_pi_and_dff_statistics(pi_mean, dff_mean, 
							pi_std_dev, dff_std_dev, 
							high_degree_boundry, low_degree_boundry);


	assert(nPI == level_node->get_nCombinational_nodes());
	assert(nPI + nDFF == level_node->get_nNodes());


	// we cannot have an unlatched pi (combinational) node at delay level 0 with 0 degree 
	// this causes a primary inputs to be hooked up to a primary output
	// which does not make sense.

	// 1. make 0 degree unlatched nodes dff
	for (i = 0; i < nNodes; i++)
	{
		node = level_node->get_node(i);	
		assert(node);

		if (node->get_fanout_degree() == 0 && ! node->is_latched())
		{
			assert(! node->is_PI());
			designate_node_as_dff(node);
			nDFF--;
		}
		else
		{
			assert(node->get_fanout_degree() != 0 || node->is_latched());
		}

	}

	if (nDFF < 0)
	{
		Fail("We had to assign too many dff at 0 delay level because we had too many 0 degree nodes.");
	}

	// 2. Assign the highest and lowest degrees next.
	// 	  
	//   
	for (i = 0; i < nNodes; i++)
	{
		node = level_node->get_node(i);	
		assert(node);

		if (! node->is_DFF() && ! node->is_PI())
		{

			fanout_degree = node->get_fanout_degree();

			if (fanout_degree >= high_degree_boundry || fanout_degree <= low_degree_boundry)
			{
				if (nDFF > 0 && nPI > 0)
				{
					// select what the node will be randomly based on their fanout 
					// distributions. our distribution will have valid and real mean and std deviation
					// if we get to this point

					assert(new_gen); // should really me sampling without replacement.

					pi_probability = g_rand_number_gen->get_gaussian_probability(fanout_degree, 
																				pi_mean, pi_std_dev);
					dff_probability = g_rand_number_gen->get_gaussian_probability(fanout_degree, 
																				dff_mean, dff_std_dev);

					total_probability = pi_probability + dff_probability;

					if (total_probability > 0)
					{
						pi_probability = pi_probability/total_probability;
						dff_probability = dff_probability/total_probability;

						random_number = g_rand_number_gen->random_double_number(1);

						if (random_number < pi_probability)
						{
							designate_node_as_pi(node);
							nPI--;
						}
						else
						{
							designate_node_as_dff(node);
							nDFF--;
						}
					}
					else
					{
						// we are far away from either distribution.
						deal_with_weird_special_case(fanout_degree, pi_mean, pi_std_dev, 
													dff_mean, dff_std_dev, high_degree_boundry,
													low_degree_boundry, nDFF, nPI, node);
					}
				}
				else if (nDFF > 0)
				{
					designate_node_as_dff(node);
					nDFF--;
				}
				else
				{
					assert(nPI > 0);
					designate_node_as_pi(node);
					nPI--;
				}
			}
		}
	}

	// 3. - iterate over the nodes that have degrees greater than 0.
	//    - figure out the probabilities that the node is a pi or dff.
	//    - randomly choose based on those weights
	for (i = 0; i < nNodes; i++)
	{
		node = level_node->get_node(i);	
		assert(node);

		if (! node->is_DFF() && ! node->is_PI())
		{

			if (nDFF > 0 && nPI > 0)
			{
				// select what the node will be randomly based on their fanout 
				// distributions
				//

				assert(new_gen); // should really me sampling without replacement.

				pi_probability = g_rand_number_gen->get_gaussian_probability(fanout_degree, 
																			pi_mean, pi_std_dev);
				dff_probability = g_rand_number_gen->get_gaussian_probability(fanout_degree, 
																			dff_mean, dff_std_dev);

				total_probability = pi_probability + dff_probability;

				assert(total_probability > 0);

				pi_probability = pi_probability/total_probability;
				dff_probability = dff_probability/total_probability;

				random_number = g_rand_number_gen->random_double_number(1);

				if (random_number < pi_probability)
				{
					designate_node_as_pi(node);
					nPI--;
				}
				else
				{
					designate_node_as_dff(node);
					nDFF--;
				}
			}
			else if (nDFF > 0)
			{
				designate_node_as_dff(node);
				nDFF--;
			}
			else
			{
				assert(nPI > 0);
				designate_node_as_pi(node);
				nPI--;
			}
		}

	}

	assert(nDFF == 0 && nPI == 0);
	assert(static_cast<unsigned>(level_node->get_nDFF()) == level_node->get_DFF().size());
	assert(static_cast<unsigned>(m_cluster->get_nPI()) == level_node->get_PI().size());
	assert(level_node->get_nDFF() + m_cluster->get_nPI() == level_node->get_nNodes());

	//debug("number dff assigned " << level_node->get_DFF().size());

}


// At the 0th delay level designate which nodes are primary inputs
// and which nodes are flip-flops
//
// Do so randomly instead of like designate_dff_and_pi_nodes() where we use
// the fanout statistics.
//
// This function produces about a < 5% worse wirelength than designate_dff_and_pi_nodes()
// in the experiments that i ran.
//
// PRE: fanouts have been assigned to the nodes
// POST: all nodes at the 0th delay level for the cluster in question
//       have been designated as either PI or flip-flop
//
void NODE_SPLITTER::randomly_designate_dff_and_pi_nodes()
{
	assert(new_gen);

	LEVEL_NODE * level_node = m_cluster->Level[0];
	assert(level_node);

	NODE * node = 0;
	NUM_ELEMENTS i,
				 random_number = 0;
	
	NUM_ELEMENTS nDFF = level_node->get_nDFF(),
				 nNodes = level_node->get_nNodes(),
				 nPI = m_cluster->get_nPI();

	// we cannot have an unlatched pi (combinational) node at delay level 0 with 0 degree 
	// this causes a primary inputs to be hooked up to a primary output
	// which does not make sense.

	// 1. make 0 degree nodes dff
	for (i = 0; i < nNodes; i++)
	{
		node = level_node->get_node(i);	
		assert(node);

		if (node->get_fanout_degree() == 0 && ! node->is_latched())
		{
			assert(! node->is_PI());
			designate_node_as_dff(node);
			nDFF--;
		}
		else
		{
			assert(node->get_fanout_degree() != 0 || node->is_latched());
		}
	}

	if (nDFF < 0)
	{
		Fail("We had to assign too many dff at 0 delay level because we had too many 0 degree nodes.");
	}

	for (i = 0; i < nNodes; i++)
	{
		node = level_node->get_node(i);	
		assert(node);

		if (! node->is_DFF() && ! node->is_PI())
		{
			random_number = g_rand_number_gen->random_number(1);

			if (nPI == 0 || (random_number = 1 && nDFF > 0))
			{
				designate_node_as_dff(node);
				nDFF--;
			}
			else
			{
				assert(nPI > 0);
				designate_node_as_pi(node);
				nPI--;
			}
		}
	}


	assert(nDFF == 0 && nPI == 0);
	assert(static_cast<unsigned>(level_node->get_nDFF()) == level_node->get_DFF().size());
	assert(static_cast<unsigned>(m_cluster->get_nPI()) == level_node->get_PI().size());
	assert(level_node->get_nDFF() + m_cluster->get_nPI() == level_node->get_nNodes());

}

// a node has a fanout that is far from the average
// of the pi and dff fanout statistics.
//
// designated it as either a PI and DFF based on which distribution is 
// closer to reaching it
//
// PRE: node is valid
//      node has a fanout that is far from avg. of the 
//      PI fanout statistics and DFF fanout statistics
// POST: node has been designated as either a PI or DFF
void NODE_SPLITTER::deal_with_weird_special_case
(
	const double& fanout_degree, 
	const double& pi_mean, 
	const double& pi_std_dev, 
	const double& dff_mean, 
	const double& dff_std_dev, 
	const double& high_degree_boundry,
	const double& low_degree_boundry,
	NUM_ELEMENTS& nDFF, 
	NUM_ELEMENTS& nPI, 
	NODE * node
)
{

	assert(new_gen);
	// weird case do occur when we have 0 total probility
	//
	// they occur when we are far from either mean

	if (fanout_degree >= high_degree_boundry)
	{
		// choose the higher of pi or dff
	
		if ((pi_mean + pi_std_dev) > (dff_mean + dff_std_dev))
		{
			designate_node_as_pi(node);
			nPI--;
		}
		else
		{
			designate_node_as_dff(node);
			nDFF--;
		}
	}
	else
	{
		assert(fanout_degree <= low_degree_boundry);

		// choose the lower of pi or dff 
		if (pi_mean - pi_std_dev < dff_mean - dff_std_dev)
		{
			designate_node_as_pi(node);
			nPI--;
		}
		else
		{
			designate_node_as_dff(node);
			nDFF--;
		}
	}
}

// 
// PRE: node is not a PO or flip-flop
// POST: node is a primary input
void NODE_SPLITTER::designate_node_as_pi
(
	NODE * node
)
{
	assert(! node->is_PO() || node->is_DFF());
	node->set_is_PI(true);
}
// 
// PRE: node is not a PI or flip-flop
// POST: node is a flip-flop
void NODE_SPLITTER::designate_node_as_dff
(
	NODE * node
)
{
	assert(! node->is_DFF() && ! node->is_PI());
	node->set_is_dff(true);
}



// PRE: 
// POST: all PO have been designated in the cluster
void NODE_SPLITTER::designate_PO()
{
	LEVEL_NODE * level_node = 0;
	DELAY_TYPE lev;

    for (lev = 0; lev <= m_cluster->get_delay(); lev++) 
	{
		level_node = m_cluster->Level[lev];
		assert(level_node);

		level_node->designate_PO();
    }
}

//
//
// POST: We are sure that:
//       - each node at the 0th delay level either either a PI of a flip-flop
//       - no node at the 0th delay level is both a PI and a POnode is a flip-flop
//       We have warned the user if we will be forced to connect a flip-flop to itself
//       during final edge assignment.
//
void NODE_SPLITTER::final_sanity_check() const
{
	NODE * node = 0;
	LEVEL_NODE * level_node = 0;
	DELAY_TYPE delay_level = 0;
	NUM_ELEMENTS i, nNodes, edges_to_assign;

	// check the 0th delay level
	level_node = m_cluster->Level[0];
	assert(level_node);
	nNodes = level_node->get_nNodes();

	for (i = 0; i < nNodes; i++)
	{
		node = level_node->get_node(i);	
		assert(node->get_delay_level() == 0);

		if (node->is_PI() && node->is_PO())
		{
			Warning("Node " << node->get_name() << " is a primary output and primary output");
			assert(false);
		}

		if (node->is_combinational() && ! node->is_PI())
		{
			Warning("Node " << node->get_name() << " is a combinational node at the 0th delay level "
					<< "and is not a PI");
			assert(false);
		}
		
		// check if we are going to force DFF to connect to themselves
		if (node->is_DFF() && node->is_latched() &&
			level_node->only_dff_super_edge_is_this_level_node())
		{
			Warning("Loop ahead");
		}
	}


	for (delay_level = 0; delay_level < m_cluster->get_delay(); delay_level++)
	{
		level_node = m_cluster->Level[delay_level];
		assert(level_node);
		nNodes = level_node->get_nNodes();

		for (i = 0; i < nNodes; i++)
		{
			node = level_node->get_node(i);	

			if (node->get_fanout_degree() == 0 && ! (node->is_PO() || node->is_latched()))
			{
				Warning("Node " << node->get_name() << 
						" has 0 fanout and is not latched or a primary output.");
				assert(false);
			}

			// edges to assign are length 1 edges (the long edges are not included)
			// if we have more length 1 edges for each node that the number of nodes at the next
			// delay level we will be forced to make double connections
			edges_to_assign	= node->get_edges_to_assign();
			assert(edges_to_assign <= m_cluster->Level[delay_level+1]->get_nNodes());
		}
	}
}

// get the fanout statistics for PIs and flip flops
//
// POST: pi_mean, dff_mean, pi_std_dev, dff_std_dev, high_degree_boundry, low_degree_boundry
//       have all been obtained from comb_spec
//
void NODE_SPLITTER::get_pi_and_dff_statistics
(
	double& pi_mean, 
	double& dff_mean, 
	double& pi_std_dev, 
	double& dff_std_dev,
	double& high_degree_boundry, 
	double& low_degree_boundry
)
{
	CLUSTER_SPEC * comb_spec = m_cluster->get_comb_spec();
	assert(comb_spec);

	pi_mean = comb_spec->get_mean_fanout_pi();
	dff_mean = comb_spec->get_mean_fanout_dff();

	pi_std_dev = comb_spec->get_std_dev_fanout_pi();
	dff_std_dev = comb_spec->get_std_dev_fanout_dff();

	// we want the lowest number of the mean + std_deviation for both distributions.
	// at this point all numbers above the highest degree boundry should be assigned to the 
	// other distribution
	high_degree_boundry = MIN(pi_mean + pi_std_dev, dff_mean + dff_std_dev);


	// we want the highest number of the mean - std_deviation for both distributions.
	// at this point all numbers below the lowest degree boundry should be assigned to the 
	// other distribution

	low_degree_boundry = MIN(pi_mean - pi_std_dev, dff_mean - dff_std_dev);
	low_degree_boundry = MAX(low_degree_boundry, 0);
}


// We have reported on the cost of the solution
//
// cost = the difference from the fanout distribution assigned to the 
// individual nodes and the spec
//
NUM_ELEMENTS NODE_SPLITTER::report_solution_quality() const
{
	assert(m_cluster);
	CLUSTER_SPEC * comb_spec = m_cluster->get_comb_spec();
	assert(comb_spec);
	
	DELAY_TYPE lev = 0,
				max_comb_delay = m_cluster->get_delay();
	LEVEL_NODE * level_node = 0;
	NUM_ELEMENTS fanout_degree = 0,
				 max_fanout = m_cluster->get_max_fanout();

	DISTRIBUTION degree_distribution(max_fanout + 1, 0);
	DISTRIBUTION spec_degree_distribution = comb_spec->get_fanout_distribution();
				 spec_degree_distribution.resize(max_fanout +1, 0);
	NUM_ELEMENTS degree_difference = 0;
	NUM_ELEMENTS spec_max_fanout_degree = max_fanout,
					max_fanout_degree = max_fanout;
	NODES::const_iterator node_iter;
	NODES nodes;
	NODE * node = 0;

	assert(spec_degree_distribution.size() == degree_distribution.size());


	// find out the degree distribution
    for (lev = 0; lev < max_comb_delay; lev++) 
	{
		level_node = m_cluster->Level[lev];
		assert(level_node);

		nodes = level_node->get_nodes();

		for (node_iter = nodes.begin(); node_iter != nodes.end(); node_iter++)
		{
			node = *node_iter;
			assert(node);

			fanout_degree = node->get_fanout_degree();
			degree_distribution[fanout_degree] += 1;
		}
	}

	// compute the cost
	while (max_fanout_degree >= 0)
	{
		// find the highest degree
		while (max_fanout_degree >= 0 && degree_distribution[max_fanout_degree] == 0)
		{
			max_fanout_degree--;
		}

		if (max_fanout_degree >= 0)
		{
			// we have a degree
			assert(degree_distribution[max_fanout_degree] > 0);

			// find the highest spec degree
			while (spec_max_fanout_degree >= 0 && spec_degree_distribution[spec_max_fanout_degree] == 0)
			{
				spec_max_fanout_degree--;
			}

			// because the number of nodes did not change if we had a real degree 
			// we should have a spec degree
			assert(spec_max_fanout_degree >= 0);
			assert(spec_degree_distribution[spec_max_fanout_degree] > 0);

			// if the degrees differ take the distance between them as an indication
			// of cost
			if (max_fanout_degree != spec_max_fanout_degree)
			{
				degree_difference += ABS(max_fanout_degree - spec_max_fanout_degree);
			}

			degree_distribution[max_fanout_degree]--;
			spec_degree_distribution[spec_max_fanout_degree]--;
		}
	}

	debug("Difference between Fanout Distribution and spec: " << degree_difference);

	return degree_difference;
}
