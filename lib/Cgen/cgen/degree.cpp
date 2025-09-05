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



#include "degree.h"
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <math.h>

DEGREE_PARTITIONER::DEGREE_PARTITIONER()
{

	m_max_comb_delay 	= 0;
	m_max_fanout		= 0; 
	m_num_out_degrees_to_assign = 0;
	m_cluster		= 0;
	m_iterating			= false;
	m_sanity_check_fatal= false;
	m_max_latched_degree = 0;

	assert(g_options);
	m_delta = g_options->get_degree_partitioning_delta();
}



//  Mainline for degree partition.
//
//  After creating the delay structure, the fanout distribution of each cluster 
//  (which is a set of fanout values that ultimately will be assigned to each individual node) 
//  is partitioned among the level nodes in the cluster. The input into this degree partitioning 
//  step is the delay structure (generated above in sequential_level->create_delay_structure()) 
//  and the fanout distribution. 
//  The output is the delay structure graph with the fanout distribution partitioned among the 
//  level nodes such that the sum of the fanout degrees assigned to each level node matching 
//  the number of edges that leave the level node.
//
//  The degree partitioning occurs in three steps. 
//
//  In the first step, an initial assignment is made. 
//  In the second step, the degree partition is iteratively improved. 
//  In the third step, post processing is done to enforce the equality described above.
//
//  PRE: cluster is valid
//  POST: fanouts have been assigned the all level nodes in the cluster
//
NUM_ELEMENTS DEGREE_PARTITIONER::degree_partition
(
	CLUSTER *cluster
)
{
    debugif(DDegree, Sep << "Inside degree partition");

	assert(cluster);
	m_cluster = cluster;
	NUM_ELEMENTS quality = 0;

	initialize_variables();
    
    debugif(DDegree,"Initial fanouts to assign " << m_num_out_degrees_to_assign);
    if (m_num_out_degrees_to_assign == 0) 
	{
		Verbose("Cluster has no nodes, leaving degree-allocation early (OK)");
		return 0;
    }

	create_initial_solution();

    debugif(DDegree,sep << "Iterating towards adherance");

	recalculate_max_degree();


	simulated_anneal();

	// Post-process if necessary 
	if (get_cost() > 0)
	{
		debugif(DDegree, sep << "Fixing degree assignments for cluster " << m_cluster->get_cluster_number());
		fix_degree_assignments();
		fix_max_degree_assignments();
		sanity_check();
	}

	prevent_double_connections();

	final_sanity_check();

	quality = report_solution_quality();

	cleanup();

	return quality;
}


// Step 1 - Create the initial solution for degree partitioning
//
// The initial fanout assignment in each cluster is made based on the fanout distribution 
// and the node shape and output shape of the delay structure. 
// It occurs in five steps with Steps 1, 4, and 5 being based on Hutton's work. 
// In Step 1, the level node at delay level dmax is assigned a fanout of 0 for 
// each node because nodes with maximum combinational delay are either latched or 
// are a primary output. In Step 2, level nodes with latches are assigned the 
// lowest unassigned fanout degrees. In Step 3, if any zero degree fanouts remain, 
// they are assigned to level nodes with POs. 
// In Step 4, pre-assignment of low fanout degrees are made to level nodes with a 
// small ratio of output degree to number of nodes. In Step 5, the remainder of the 
// fanouts are assigned beginning with the largest fanouts based on the fanout capacity of the level nodes.
// 
// PRE: 
// POST: we have created the initial solution
//
void DEGREE_PARTITIONER::create_initial_solution()
{
	DELAY_TYPE lev = 0;


    debugif(DDegree,"Computing the target out-degree for each level");
    compute_targets();

	// Step 1,2,3 and 4
	make_initial_assignments();


	// Step 5: the remainder of the 
	// fanouts are assigned beginning with the largest fanouts based on 
	// the fanout capacity of the level nodes.
	//
    debugif(DDegree,"Beginning the main assignment loop");
    while (m_num_out_degrees_to_assign > 0 && m_max_fanout > 1) 
	{
        debugif(DDegree, sep << "Have " << m_num_out_degrees_to_assign << " out-degrees remaining to assign");
        //show_vectors();
        assign_one_out();
    }
	assert(accumulate(m_fanout_distribution.begin(), m_fanout_distribution.end(), 0) == 
			m_num_out_degrees_to_assign);

    debugif(DDegree,sep << "Cleaning up with unit out-degrees");

    for (lev = 0; lev < m_max_comb_delay; lev++) 
	{
		assert(new_gen); // added fanout dist term
		while (m_unassigned_nodes[lev] > 0 && m_fanout_distribution[1] > 0)
		{
			debugif(DDegree,"..assigning unit out-degree to level " << lev);
			make_assignment(1, lev);
		}
		assert(m_unassigned_nodes[lev] == 0);
    }

	assert(m_fanout_distribution[1] == 0);

	check_degree_nNodes_match();
}

//
//  Compute the target degree assignments for each level.  We will try to 
//  come as close as possible to these targets.
//
//  PRE:  nothing
//  POST: m_targets contains the output degree of the 
//        individual level nodes
void DEGREE_PARTITIONER::compute_targets()
{
	DELAY_TYPE lev = 0;

    debugif(DDegree,"Computing target out-degree");
    m_target.resize(m_max_comb_delay+1, 0);

	for (lev = 0; lev < m_max_comb_delay; lev++) 
	{
	    m_target[lev] = m_cluster->Level[lev]->get_available_out();
	}

}

// 
// PRE: m_cluster is valid
// POST: variables have been assigned values and
//       assigned memory space
void DEGREE_PARTITIONER::initialize_variables()
{
    m_max_comb_delay = m_cluster->get_delay();
    m_num_out_degrees_to_assign = m_cluster->get_nNodes();

	m_max_latched_degree = 0;
    m_max_fanout = m_cluster->get_max_fanout();
    m_fanout_distribution = m_cluster->get_comb_spec()->get_fanout_distribution();
    m_unassigned_nodes.resize(m_max_comb_delay + 1, 0);
    m_sum_assigned.resize(m_max_comb_delay+1, 0);

	// there should only be one fanout for each node
	assert(accumulate(m_fanout_distribution.begin(), m_fanout_distribution.end(), 0) == 
			m_num_out_degrees_to_assign);

    m_maxdeg.resize(m_max_comb_delay+1, 0);
    m_mindeg.resize(m_max_comb_delay+1, 0);
    m_avgdeg.resize(m_max_comb_delay+1, 0);
    m_pmf.resize(m_max_comb_delay+1, 0);
	m_current_maxdeg.resize(m_max_comb_delay+1, 0);

}



// 
// Steps 1 through 3 in the creation of the initial solution
//
// POST: the level node at delay level dmax is assigned a fanout of 0 for 
//       each node
//       level nodes with latches are assigned the lowest unassigned fanout degrees
//       if any zero degree fanouts remain, they are assigned to level nodes with POs. 
void
DEGREE_PARTITIONER::make_latch_and_zero_degree_assignments()
{
    int lev;
    LEVEL_NODE *level_node = 0;
	NUM_ELEMENTS check = 0;
	NUM_ELEMENTS nLatched = 0;
	NUM_ELEMENTS nPO = 0;
	NUM_ELEMENTS to_assign = 0;


	// Step 1 - the level node at delay level dmax is assigned a fanout of 0 for 
	// each node because nodes with maximum combinational delay are either latched or 
	// are a primary output.
	level_node = m_cluster->Level[m_max_comb_delay];
	assert(level_node);
	level_node->OutDegrees(0) = level_node->get_nNodes();
	m_num_out_degrees_to_assign -= level_node->get_nNodes();
	m_unassigned_nodes[m_max_comb_delay] = 0;
	m_fanout_distribution[0] -= level_node->get_nNodes();
	assert(m_fanout_distribution[0] >= 0);

    for (lev = m_max_comb_delay-1; lev >= 0; lev--) 
	{
		assert(lev >= 0 && lev < m_max_comb_delay);

		level_node = m_cluster->Level[lev];
		assert(level_node);

		m_unassigned_nodes[lev] = level_node->get_nNodes();

		nLatched = level_node->get_nLatched();

		debugif(DDegree, "Assigning " << nLatched << " latched outputs at level " << lev);

		// Step 2 - level nodes with latches are assigned the lowest unassigned fanout degrees
		while (nLatched > 0)
		{
			to_assign = MIN(nLatched, m_fanout_distribution[m_max_latched_degree]);

			assign_to_level(m_max_latched_degree, lev, to_assign);
			assert(m_fanout_distribution[m_max_latched_degree] >= 0);

			nLatched -= to_assign;
	
			if (m_fanout_distribution[m_max_latched_degree] == 0 && nLatched > 0)
			{
				m_max_latched_degree++;
			}
		}

		check += m_unassigned_nodes[lev];
    }
    assert(check == m_num_out_degrees_to_assign);

	// Step 3 - if any zero degree fanouts remain, 
	// they are assigned to level nodes with POs. 
	// the only nodes that have zero fanouts are PO.
	lev = m_max_comb_delay-1;
	while (m_fanout_distribution[0] > 0)
	{
		assert(lev >= 0 && lev < m_max_comb_delay);

		level_node = m_cluster->Level[lev];
		assert(level_node);

		nPO	= level_node->get_nPO();
		nPO = MIN(nPO, m_fanout_distribution[0]);
		nPO = MIN(nPO, m_unassigned_nodes[lev]);

		debugif(DDegree,"Assigning " << nPO << " PO outputs at level " << lev);

		assign_to_level(0, lev, nPO);
		check -= nPO;

		lev--;
    }

	assert(accumulate(m_fanout_distribution.begin(), m_fanout_distribution.end(), 0) == 
			m_num_out_degrees_to_assign);
	
    assert(check == m_num_out_degrees_to_assign);
}


// Steps 1 through 4 in the creation of the initial solution
//
// In Step 1, the level node at delay level dmax is assigned a fanout of 0 for 
// each node because nodes with maximum combinational delay are either latched or 
// are a primary output. In Step 2, level nodes with latches are assigned the 
// lowest unassigned fanout degrees. In Step 3, if any zero degree fanouts remain, 
// they are assigned to level nodes with POs. 
// In Step 4, pre-assignment of low fanout degrees are made to level nodes with a 
// small ratio of output degree to number of nodes. 
//
// PRE: nothing
// POST: Steps 1 through 4 done.

void DEGREE_PARTITIONER::make_initial_assignments()
{
	DELAY_TYPE lev = 0;
	NUM_ELEMENTS check = 0;

    debugif(DDegree,"Making initial assignments");

	// Step 1 and 2 and 3
    make_latch_and_zero_degree_assignments();

    show_vectors();


	// Step 4 pre-assignment of low fanout degrees are made to level nodes with a 
	// small ratio of output degree to number of nodes
    for (lev = 0; lev < m_max_comb_delay; lev++) 
	{
        initial_assignments(lev);
    }

	show_vectors();

    check = 0;
    for (lev = 0; lev < m_max_comb_delay; lev++) 
	{
		debugif(DDegree, "Computing max/min/avg for level " << lev);

		recalculate_for_level(lev);
		check += m_unassigned_nodes[lev];
    }
    assert(check == m_num_out_degrees_to_assign);

	recalculate_max_fanout();
}




//
//  Look for levels which have the ratio of out_edges/out_nodes <= 2, and
//  assign 0, 1, 2 out-degrees to these to get rid of them.  This is to make
//  sure we don't accidentally (even though the prob is pretty low) give them
//  a high-degree node and then get stuck making a bunch of zeros later.
// 
//  PRE: lev is the delay level we want to assign values to
//  POST:  we have made preassignment of 0,1,2 degrees to the level node
// 
void DEGREE_PARTITIONER::initial_assignments
(
	const DELAY_TYPE& lev
)
{
	assert(lev < m_max_comb_delay);

    NUM_ELEMENTS num_zeros, num_ones, num_twos, assigned;
    double lowtarget;
    LEVEL_NODE *level_node = m_cluster->Level[lev];

	assert(level_node);

    assigned = level_node->get_available_out();

	if (level_node->get_nNodes() == 0)
	{
    	debugif(DDegree, "Pre-assign: no nodes at level " << lev);
		assert(assigned == 0);
		return;
	}


    debugif(DDegree, "Pre-assign: lev " << lev << " with " << m_unassigned_nodes[lev] << " nodes has target " << m_target[lev]);

	lowtarget = static_cast<double>(m_target[lev]);

	// low fanout levels are levels with fanout less than 2*nodes
    if ((lowtarget / static_cast<double>(m_unassigned_nodes[lev])) > 2.0) 
	{
		debugif(DDegree,"..no pre-assign for this level -- lots of edges");
        return;
    } 

    num_zeros = MAX(m_unassigned_nodes[lev] - m_target[lev], 0);
	num_zeros = MIN(num_zeros, m_fanout_distribution[0]);
    assert(num_zeros >= 0);

	assert(new_gen);

	// not sure if this is a good thing in newgen
	// the "m_unassigned_nodes[lev] > m_cluster->pSpec->cNodes/m_cluster->delay" term
	// is to prevent the wide levels in a thin combinational circuit from being
	// pinned down to low fanout.  We need high fanout for this level because we 
	// need things to connect
    if (num_zeros == 0 && m_unassigned_nodes[lev] > m_cluster->get_nNodes()/m_max_comb_delay) 
	{
		debugif(DDegree,"..no pre-assignment for this level -- lots of nodes ");
        return;
    }


	assign_to_level(0, lev, num_zeros);

	if (m_fanout_distribution.size() - 1 < 2)
	{
		debugif(DDegree,"..max fanout less than 2. no further pre-assignments");
		return;
	}

    if (m_target[lev] == 2 * m_unassigned_nodes[lev]) 
	{
        num_twos = m_unassigned_nodes[lev];
    } 
	else 
	{
        num_twos = MIN(m_target[lev] - m_unassigned_nodes[lev], m_unassigned_nodes[lev]);
    }
    num_twos = MAX(num_twos, 0);
	num_twos = MIN(num_twos, m_fanout_distribution[2]);
    num_ones = m_unassigned_nodes[lev] - num_twos;
	num_ones = MIN(num_ones, m_fanout_distribution[1]);

    debugif(DDegree, "..Assigning " << num_zeros << " zeros, " << num_ones << " ones and " << num_twos <<
			" twos (totalling " << m_unassigned_nodes[lev] << ")");
   
    if (num_ones + 2*num_twos > m_target[lev]) 
	{
		if (DDegree)
		{
			Verbose("Noticing an over-assignment of fanout for level " << lev);
		}
    }

    assert(num_ones >= 0);
    assert(num_twos >= 0);

	assign_to_level(1, lev, num_ones);
	assign_to_level(2, lev, num_twos);
}


//
//  Debug progress of the assignment.
//
void DEGREE_PARTITIONER::show_vectors() const
{

    int i;
	LEVEL_NODE * level_node = 0;

    debugif(DDegree,"=================== Degree Assignments =====================");
    debugif(DDegree,"  Orig   Orig Orig Assd Assd Trgt | Curr  To  | Max  Avg   Min");
    debugif(DDegree,"  Nodes nLat  sum  sum zero  sum | num  Assn | deg  deg   deg");
    debugif(DDegree," ---- ---- ---- ---- ---- ---- | ---- ---- | ---- ---- ----");
    for (i = 0; i <= m_max_comb_delay; i++) 
	{
		level_node = m_cluster->Level[i];
		assert(level_node);

		debugif(DDegree," " << setprecision(3) << setw(4) << 
			level_node->get_nNodes() << " " << setw(4) << 
			level_node->get_nLatched() << " " << setw(4) <<
			level_node->get_available_out() << " " << setw(4) <<
			m_sum_assigned[i] << " " << setw(4) << 
			level_node->OutDegrees(0) << " " << setw(4) <<
			m_target[i] << " | " << setw(4) << 

			m_unassigned_nodes[i] << " " << setw(4) << 
			m_target[i] - m_sum_assigned[i] << " | " << setw(4) <<
			m_maxdeg[i] << " " << setw(4) << m_avgdeg[i] << " " << m_mindeg[i]);
    }
    dbgif(DDegree,Sep << "Remaining out-deg: ");
    for (i = 0; i <= m_max_fanout; i++) 
	{ 
		dbgif(DDegree,m_fanout_distribution[i] << " ");
    }
    dbgif(DDegree,"\n");
}


// 
//  Assign one out-degree to some level.  
//
//  This is the "standard" allocation of a single out-degree in the main loop.
//
//  First we look for needy levels:  if there is a level which has a min-deg
//  requirement of what we have to assign, just give it that degree right away.
//  Otherwise, build a pmf of desire for the degree (based on avg degree 
//  needed) and choose probabilistically who gets it.
//
//  PRE: m_max_fanout is a degree that needs to be assigned
//  POST: we have assigned m_max_fanout to a level node
//        we have assigned 0,1,2 degree assignments to the level node
//        if it needed it after our assignment 
///
void DEGREE_PARTITIONER::assign_one_out()
{

    int deg, i;
	DELAY_TYPE choice = 0;
	COST_TYPE tot_pmf = 0.0;

    debugif(DDegree,sep);
    deg = m_max_fanout;
    //debugif(DDegree,"Trying to assign out-degree " << deg << 
	//" (are " << m_fanout_distribution[deg] << " of them).");

	assert(deg > 1);
	assert(m_fanout_distribution[deg] > 0);

	// we start with high degree values.
	// as such if highest_deg < avdeg we just give that level the avgdeg
	// and let it as such
	//
	// the same thing applies with m_mindeg
	//
	// otherwise create the pmf

	
    for (i = 0; i < m_max_comb_delay; i++)  
	{
		if (deg <= m_avgdeg[i]) 
		{
			make_assignment(deg, i);

			return;
		} 

		if (deg < m_maxdeg[i])
		{
			m_pmf[i] = m_avgdeg[i]*m_avgdeg[i];
			tot_pmf += m_pmf[i];
		}
		else 
		{
			m_pmf[i] = 0;
		}

    }

	
	NUM_ELEMENTS nChoices = static_cast<NUM_ELEMENTS>(m_max_comb_delay);

	if (tot_pmf > 0.0)
	{
		choice = (short) g_rand_number_gen->discrete_pmf(nChoices, m_pmf);
	}
	else
	{
		 choice = (short) g_rand_number_gen->discrete_pmf(nChoices, m_unassigned_nodes);
	}

    debugif(DDegree,"..chose to assign " << deg << " to level " << choice);

    make_assignment(deg, choice);
}


//
//  Make the assignment decided upon in assign_one_out.
//
//  Note that this is only done for single degrees.  Assignments of multiple
//  zero/one/two in the precalculation, or in recalculate are done internally 
//  in those places.
// 
//  PRE: deg is a degree to assign
//       lev is the delay level to assign it to
//  POST: we have assigned the degree
//        we have recalculated:
//			- the max. fanout to assign 
//			- the various min, max degree statistics for the level nodes
//        we have assigned 0,1,2 degree assignments to the level
//        if it need them after our assignment
void
DEGREE_PARTITIONER::make_assignment
(
	const NUM_ELEMENTS& deg, 
	const DELAY_TYPE& lev
)
{
	assert(m_cluster);
    assert(lev >= 0 && lev <= m_max_comb_delay);
	assert(deg >= 0 && deg <= m_max_fanout);
	assert(m_fanout_distribution[deg] > 0);
	assert(deg == m_max_fanout);

    debugif(DDegree,"Making the degree asssignment of " << deg << " to level " << lev);

	assign_to_level(deg, lev, 1);
 
	recalculate_max_fanout();

    recalculate_for_level(lev);

	recalculate_max_fanout();
}


//
//  Re-calculate the status vectors after making a degree assignment 
//  to a level.
// 
//
//  PRE: lev is the level node to recalculate for
//  POST: we have re-calculate the status vectors
//        we have assigned 0,1,2 degrees to this level if 
//        it need it after our assignment
void
DEGREE_PARTITIONER::recalculate_for_level
(
	const DELAY_TYPE& lev
)
{
	assert(m_cluster);

    LEVEL_NODE *level_node = 0;
    int tmp, num_ones;
    int oldmaxdeg, oldmindeg;
    double oldavgdeg;
	NUM_ELEMENTS maxdegree_without_double_inputs = 0;

    level_node = m_cluster->Level[lev];
	assert(level_node);

	// just for display purposes. not used in any calculation below
    oldmaxdeg = m_maxdeg[lev];
    oldmindeg = m_mindeg[lev];
    oldavgdeg = m_avgdeg[lev];
   
    debugif(DDegree,"..Level " << lev << " has sum=" << m_sum_assigned[lev]);

    //  If the average out-degree requirement falls below 1.5, assign 
    //  some ones to get it back up to 2
	//  

    if (m_unassigned_nodes[lev] > 0 && m_target[lev] - m_sum_assigned[lev] < 1.5 * m_unassigned_nodes[lev] )
	{
		num_ones = 2 * m_unassigned_nodes[lev] - (m_target[lev] - m_sum_assigned[lev]);
		num_ones = MAX(num_ones, 0);
		num_ones = MIN(num_ones, m_unassigned_nodes[lev]);
		num_ones = MIN(num_ones, m_fanout_distribution[1]);
		assert(num_ones >= 0);
		debugif(DDegree,"** Assigning " << num_ones << " ones to level " << lev << " with " 
				<< m_unassigned_nodes[lev] << " nodes, target " << m_target[lev]);
		debugif(DDegree,"fanout distribution of 1 " << m_fanout_distribution[1]);
		assert(m_unassigned_nodes[lev] >= num_ones && num_ones >= 0);

		assign_to_level(1, lev, num_ones);
    }

    //  If we have assigned all out-degrees for this level, assign everything
    //  to zero, so we don't have any residual effects from rounding or slop
    //  assignments.

    if (m_unassigned_nodes[lev] == 0) 
	{
		debugif(DDegree,"..out of nodes, updating everything to 0");
		m_maxdeg[lev] = 0;
		m_avgdeg[lev] = 0.0;
		m_mindeg[lev] = 0;
		return;
    }

    // Remainder of routine re-calculates min/max/avg out
  
	// max degree:
	// m_maxdeg is also limited by max_out bounds subtract what we have assigned 
	// subtracted m_unassigned_nodes which should at least a fanout of degree 1 (minus the max fanout
	// node we are considering)

	maxdegree_without_double_inputs = level_node->get_max_degree_possible();
	
    m_maxdeg[lev] = m_target[lev] - m_sum_assigned[lev] - (m_unassigned_nodes[lev] - 1);
    m_maxdeg[lev] = MIN(m_maxdeg[lev], m_cluster->get_comb_spec()->get_max_fanout());
	m_maxdeg[lev] = MIN(m_maxdeg[lev], maxdegree_without_double_inputs);
	m_maxdeg[lev] = MAX(m_maxdeg[lev], 0);

	// average degree:
    m_avgdeg[lev] = 	static_cast<double>(m_target[lev] - m_sum_assigned[lev]) / 
					static_cast<double>(m_unassigned_nodes[lev]);

    if (m_avgdeg[lev] < 0.0) 
	{
		debugif(DDegree,"..Forcing avgdeg, maxdeg, mindeg all to 0");
		m_avgdeg[lev] = 0.0;
		m_maxdeg[lev] = 0;
		m_mindeg[lev] = 0;
   	} 
	
	// min degree:
    if (m_unassigned_nodes[lev] == 1) 
	{
		m_mindeg[lev] = level_node->get_available_out() - m_sum_assigned[lev];
		debugif(DDegree,"..initial mindeg is " << m_mindeg[lev]);
    } 
	else 
	{
		// if all the rest of the nodes have max_fanout except one
		// what is the minimum degree possible
        tmp = (m_unassigned_nodes[lev] - 1) * m_max_fanout;
        m_mindeg[lev] = level_node->get_available_out() - tmp;
		debugif(DDegree,"..inital mindeg is " << m_mindeg[lev]);
    }
    m_mindeg[lev] = MAX(m_mindeg[lev], 0);

    m_mindeg[lev] = MIN(m_mindeg[lev], m_maxdeg[lev]);
    m_mindeg[lev] = MIN(m_mindeg[lev], m_max_fanout);
	debugif(DDegree,"....mindeg now " << m_mindeg[lev]);

    m_avgdeg[lev] = MAX(m_avgdeg[lev], m_mindeg[lev]);

	assert(m_avgdeg[lev] >= 0);
	assert(m_mindeg[lev] >= 0);
	assert(m_maxdeg[lev] >= 0);

    debugif(DDegree,"..updating maxdeg from " << oldmaxdeg << " to " << m_maxdeg[lev]);
    debugif(DDegree,"..updating avgdeg from " << oldavgdeg << " to " << m_avgdeg[lev]);
    debugif(DDegree,"..updating mindeg from " << oldmindeg << " to " << m_mindeg[lev]);
}



// 
// Sanity check the degree assignment
//
// POST: we have checked the sanity or 
//       failed and exited if something was wrong
//
void DEGREE_PARTITIONER::final_sanity_check() 
{
	m_sanity_check_fatal = true;
	sanity_check();
}
// Sanity check the degree assignment
//
// POST: We have warned the user if we saw something funny
//       in terms of data structures being consistent or 
//       in terms of node violations.
// 		
// 		 If we didn't have to warn then we are sure that:
// 		 - check_degree_nNodes_matched returned no errors
//		 - we have a latched node or PO for every 0 degree node
//		 - we have not assigned a fanout that will force node violations
//		 - m_sum_assigned is consistent with what has been assigned to the level node
//		 - m_fanout_distribution is consistent with m_num_out_degrees_to_assign
//		 - every fanout is between 0 and max_fanout
//
//       If did had degree assignments that will cause node violations
//       and if m_sanity_check_fatal then
//       we failed and exited the program
//
void DEGREE_PARTITIONER::sanity_check() const
{
    int j, sum, tot_out, tot_in, lost_edges;
	NUM_ELEMENTS max_fanout = m_cluster->get_max_fanout();
    LEVEL_NODE *level_node = 0;
	bool warning_found = false;
	DELAY_TYPE lev = 0;
	NODES nodes;
	NODES::iterator node_iter;
	NODE * node = 0;
	NUM_ELEMENTS fanout_degree = 0;

    debugif(DDegree,Sep << "Final results:");
    show_vectors();

	assert(accumulate(m_fanout_distribution.begin(), m_fanout_distribution.end(), 0) == 
			m_num_out_degrees_to_assign);

	check_degree_nNodes_match();


    dbgif(DDegree, Sep << "Results of degree assignments:");
   
    lost_edges = 0;
    tot_out = tot_in = 0;
    for (lev = 0; lev <= m_max_comb_delay; lev++) 
	{
		// we have no left over degrees. 
		assert(m_unassigned_nodes[lev] == 0);	

		level_node = m_cluster->Level[lev];
		assert(level_node);

		// make sure we either have a latch or a po for each 0 degree node
		assert(level_node->get_nLatched() + level_node->get_nPO() >= level_node->OutDegrees(0));

        sum = 0;
		dbgif(DDegree,"-- " << lev << " (" << level_node->get_nNodes() << "/" << m_sum_assigned[lev] << "):  ");
		for (j = 0; j <= max_fanout; j++) 
		{
			dbgif(DDegree,level_node->OutDegrees(j) << " ");
			sum += level_node->OutDegrees(j) * j;
		}
		dbgif(DDegree,"\n");
		tot_out += sum;
		debugif(DDegree,"....need " << sum << " == " << m_sum_assigned[lev]); 
		assert(sum==m_sum_assigned[lev]);

		if (level_node->get_max_degree_assigned() > level_node->get_max_degree_possible())
		{
			Warning("We have a max poss degree of " << m_maxdeg[lev] << 
				" but we have a assigned this level a degree of " << level_node->get_max_degree_assigned());
			warning_found = true;
		}

		// check to make sure out max degree from comb graph is ok
		nodes = level_node->get_nodes();
		for(node_iter = nodes.begin(); node_iter != nodes.end(); node_iter++)
		{
			node = *node_iter;
			assert(node);

			fanout_degree = node->get_fanout_degree();
			assert(fanout_degree >= 0 && fanout_degree <= m_max_fanout);
		}
		
    }
    debugif(DDegree, sep << "All done with degree sanity" << Sep);

	if (warning_found && m_sanity_check_fatal)
	{
		Fail("Couldn't fix errors.");
	}
}

// PRE: degree is teh degre to assign
//      lev is the delay level
//      num is the number to assign
// POST: we have assigned the degrees
//

void DEGREE_PARTITIONER::assign_to_level
(
	const NUM_ELEMENTS& degree,
	const DELAY_TYPE& lev, 
	const NUM_ELEMENTS& num
)
{
	if (num == 0)
	{
		return;
	}
	assert(m_cluster);
    assert(lev >= 0 && lev < m_max_comb_delay);

	if (! m_iterating)
	{
		assert(degree >= 0 && degree <= m_max_fanout);
		assert(m_fanout_distribution[degree] >= 0);
		assert(m_num_out_degrees_to_assign >= 0);
	}


    LEVEL_NODE *level_node = m_cluster->Level[lev];
	assert(level_node);

	m_unassigned_nodes[lev] -= num;

	level_node->OutDegrees(degree) += num;
	m_num_out_degrees_to_assign -= num;
	m_sum_assigned[lev] += degree*num;
	m_fanout_distribution[degree] -= num;

	assert(level_node->OutDegrees(degree) >= 0);

	if (! m_iterating)
	{
		assert(m_fanout_distribution[degree] >= 0);
		assert(m_num_out_degrees_to_assign >= 0);
		assert(m_unassigned_nodes[lev] >= 0);
	}
}

// cleanup
//
// POST: we are ready for another call to degree_partition
//
void DEGREE_PARTITIONER::cleanup()
{
	m_max_comb_delay 	= 0;
	m_max_fanout		= 0; 
	m_num_out_degrees_to_assign = 0;

    m_fanout_distribution.clear();
    m_unassigned_nodes.clear();
    m_sum_assigned.clear();

    m_maxdeg.clear();
    m_mindeg.clear();
    m_avgdeg.clear();
    m_pmf.clear();
	m_current_maxdeg.clear();

	m_cluster = 0;
}


// POST: We are sure that:
// 		 - m_sum_assigned is consistent with the degrees assigned to the level node
// 		 - m_num_out_degrees_to_assign is consistent with the degrees we have still
// 		   to assign
void DEGREE_PARTITIONER::check_degree_nNodes_match() const
{
	int i;
	NUM_ELEMENTS degree = 0;
	NUM_ELEMENTS check = 0;
	LEVEL_NODE * level_node = 0;
	NUM_ELEMENTS max_degree = m_cluster->get_max_fanout();
	NUM_ELEMENTS level_check = 0;
	NUM_ELEMENTS level_sum = 0;

    for (i = 0; i <= m_max_comb_delay; i++) 
	{
		level_node = m_cluster->Level[i];
		assert(level_node);
	

		level_check = 0;
		level_sum = 0;
		for (degree = 0; degree <= max_degree; degree++)
		{
			check += level_node->OutDegrees(degree);
			level_check += level_node->OutDegrees(degree);
			level_sum += level_node->OutDegrees(degree) * degree;
		}
		assert(level_check == level_node->get_nNodes());
		assert(level_sum == m_sum_assigned[i]);
	}

	assert(check == m_cluster->get_nNodes());

	assert(accumulate(m_fanout_distribution.begin(), m_fanout_distribution.end(), 0) == 
				m_num_out_degrees_to_assign);
}


// a quenched simulated anneal.
// After the initial assignment, the degree assignment is 
// improved by swapping fanout degrees between level nodes in an 
// attempt to improve the solution quality defined as: 
//
// cost =  cost_Degree_Fanout_Edge_Misassignment + cost_Fanout_Penalty
//
// Here cost_Fanout_Edge_Missassignment measures for each level node the 
// absolute difference between the sum of the fanout degrees assigned and 
// the number of output edges that were assigned in Section 4.1 and 
// cost_Fanout_Penality is a cost that penalizes level nodes with fanouts 
// that will force node violations in the final edge assignment 
// (described below in Section 4.4). The degrees to be swapped are 
// chosen by first, randomly choosing a non-zero fanout degree from the 
// level nodes with higher total fanout than number of output edges 
// assigned. Next, we randomly choose a non-zero lower fanout degree 
// from the level nodes with lower total fanout than the number of 
// output edges assigned. Zero fanout degrees are never chosen 
// because they can only be assigned to either a latched node or a node 
// with a primary output and hence their movement is more complicated. 
//
// 
// We employ a hill climbing iterator, in which all changes that improve 
// the cost function are accepted, and bad moves are accepted with a probability 
// that decreases exponentially with the change in cost. Moves are continually 
// generated in the algorithm until there is no change in the lowest cost for 
// 10000 iterations or until the total cost is zero.

void DEGREE_PARTITIONER::simulated_anneal()
{
	COST_TYPE cost = get_cost();
	COST_TYPE lowest_cost = cost;
	COST_TYPE changed_cost = 0.0;
	NUM_ELEMENTS loops_without_lowest_cost_change = 0;
	DELAY_TYPE source_level;
	NUM_ELEMENTS source_degree;
	NUM_ELEMENTS iteration_limit = 10000;

	m_iterating = true;

	// reset the max fanout to the maximum fanout in the cluster
	m_max_fanout = m_cluster->get_max_fanout();

	double random_number = 0.0;
	double temperature = g_options->get_degree_init_temperature();

	debugif(DDegree,"Initial cost is " << cost);


	debug("********************** Degree cost_Fanout_Penalty " << get_cost_fanout_penalty() << " **********");
	debug("********************** Degree cost_Degree_Fanout_Edge_Misassignment " 
			<< get_cost_degree_fanout_edge_misassignment() << " **********");

	while (cost >= 1.0 && loops_without_lowest_cost_change < iteration_limit)
	{
		generate_move();

		source_level = m_move.get_source_level();
		source_degree = m_move.get_source_degree();

		if (m_move.is_valid())
		{
			changed_cost = get_changed_cost();


			assert(new_gen);
			// the loops_without_lowest_cost_change is kind of weird.
			// should investigate making this algorithm a real 
			// sim. annealing algorithm

			// we are nearing the end of the algorithm
			// just accept moves that lower the cost
			if (loops_without_lowest_cost_change > 8000)
			{
				random_number = 1;
			}
			else
			{
				random_number = g_rand_number_gen->random_double_number(1);
			}

			if (changed_cost <= 0.0 || random_number < exp(- changed_cost/temperature))
			{
				assert(m_cluster->Level[source_level]->OutDegrees(source_degree) > 0);
				make_move();
				check_degree_nNodes_match();
			
				//debugif(DDegree,"cost " <<  cost);

				cost += changed_cost;
				//show_vectors();

				if (cost < lowest_cost)
				{
					lowest_cost = cost;
					//debugif(DDegree,"********************** Lowest cost " << lowest_cost << " **********");
					
					loops_without_lowest_cost_change = 0;
				}
				else
				{
					loops_without_lowest_cost_change++;
				}
			}
			else
			{
				loops_without_lowest_cost_change++;
			}
		}
		else
		{
			loops_without_lowest_cost_change++;
		}

	}

	debug("********************** Degree Final cost " << cost << " *******************");
	debugif(DDegree,"********************** Degree Final Lowest cost " << lowest_cost << " **********");
	debug("********************** Degree cost_Fanout_Penalty " << get_cost_fanout_penalty() << " **********");
	debug("********************** Degree cost_Degree_Fanout_Edge_Misassignment " 
			<< get_cost_degree_fanout_edge_misassignment() << " **********");
}

//
//
// cost =  cost_Degree_Fanout_Edge_Misassignment + cost_Fanout_Penalty
//
// RETURN: cost of the solution
//
COST_TYPE DEGREE_PARTITIONER::get_cost() const
{
	COST_TYPE cost = 0.0,
			  cost_Fanout_Penalty = 0.0,
			  cost_Degree_Fanout_Edge_Misassignment = 0.0;

	cost_Degree_Fanout_Edge_Misassignment = get_cost_degree_fanout_edge_misassignment();

	cost_Fanout_Penalty = get_cost_fanout_penalty();

	cost =  m_delta*cost_Degree_Fanout_Edge_Misassignment + (1-m_delta)*cost_Fanout_Penalty;

	return cost;
}


//
// RETURN: the change of cost the m_move will result in
//
COST_TYPE DEGREE_PARTITIONER::get_changed_cost()
{

	//debugif(DDegree,"source level " << source_level << " dest_level " << dest_level << " source degree " << source_degree << " sink degree " << dest_degree );

	COST_TYPE changed_cost = 0.0;
	NUM_ELEMENTS changed_cost_Degree_Fanout_Edge_Misassignment = 0,
			     changed_cost_Fanout_Penalty = 0;
	NUM_ELEMENTS too_many_input_cost = 0;
	LEVEL_NODE * level_node = 0;

	DELAY_TYPE source_level = m_move.get_source_level(),
				 dest_level = m_move.get_dest_level();
	NUM_ELEMENTS source_degree = m_move.get_source_degree(),
				 dest_degree = m_move.get_dest_degree(),
				 degree_difference = source_degree - dest_degree;

	assert(m_cluster && source_level != dest_level);

	changed_cost_Degree_Fanout_Edge_Misassignment = 	
				ABS(m_target[source_level] - m_sum_assigned[source_level] + degree_difference) +
				ABS(m_target[dest_level] - m_sum_assigned[dest_level] - degree_difference) - 
				ABS(m_target[source_level] - m_sum_assigned[source_level]) -
				ABS(m_target[dest_level] - m_sum_assigned[dest_level]);

	//debugif(DDegree,"changed cost from sum/target mismatch=" << int_changed_cost);


	// if our source degree is more than what it should be at the source level
	// calculate the fanout penalty cost
	if (source_degree > m_maxdeg[source_level] || dest_degree > m_maxdeg[source_level])
	{
		level_node = m_cluster->Level[source_level];
		assert(level_node);
	
		// temporarily make the move
		level_node->OutDegrees(source_degree) -= 1;	
		level_node->OutDegrees(dest_degree) += 1;	
		
		// cost after the temp. move
		too_many_input_cost = level_node->get_max_degree_assigned() - m_maxdeg[source_level];
		too_many_input_cost = MAX(too_many_input_cost, 0);
		changed_cost_Fanout_Penalty += too_many_input_cost;
		
		// unmake the temporary move
		level_node->OutDegrees(source_degree) += 1;	
		level_node->OutDegrees(dest_degree) -= 1;	
		assert(m_current_maxdeg[source_level] == level_node->get_max_degree_assigned());

		// subtract the cost before move
		too_many_input_cost = (m_current_maxdeg[source_level] - m_maxdeg[source_level]);
		too_many_input_cost = MAX(too_many_input_cost, 0);

		changed_cost_Fanout_Penalty -= too_many_input_cost;

	}
	

	// if our source degree is more than what it would be at the destination level
	// calculate the fanout penalty cost
	if (source_degree > m_maxdeg[dest_level] || dest_degree > m_maxdeg[dest_level])
	{
		level_node = m_cluster->Level[dest_level];
		assert(level_node);

		// temporarily make the move
		level_node->OutDegrees(dest_degree) -= 1;	
		level_node->OutDegrees(source_degree) += 1;	
		
		// cost after the move
		too_many_input_cost = level_node->get_max_degree_assigned() - m_maxdeg[dest_level];
		too_many_input_cost = MAX(too_many_input_cost, 0);
		changed_cost_Fanout_Penalty += too_many_input_cost;
		
		// unmake the temporary move
		level_node->OutDegrees(dest_degree) += 1;	
		level_node->OutDegrees(source_degree) -= 1;	
		assert(m_current_maxdeg[dest_level] == level_node->get_max_degree_assigned());
	
		// subtract the cost before move
		too_many_input_cost = (m_current_maxdeg[dest_level] - m_maxdeg[dest_level]);
		too_many_input_cost = MAX(too_many_input_cost, 0);
		changed_cost_Fanout_Penalty -= too_many_input_cost;
	}

	//debugif(DDegree,"total changed cost=" << int_changed_cost);
	changed_cost = m_delta*static_cast<double>(changed_cost_Degree_Fanout_Edge_Misassignment) + 
			      (1-m_delta)*static_cast<double>(changed_cost_Fanout_Penalty); 
	return changed_cost;
}


//
// RETURNS: cost_Fanout_Penalty
//
COST_TYPE DEGREE_PARTITIONER::get_cost_fanout_penalty() const
{

	DELAY_TYPE i;
	LEVEL_NODE * level_node = 0;
	NUM_ELEMENTS number_zero_degrees_possible = 0;
	COST_TYPE cost_Fanout_Penalty = 0.0;

	// sum up the cost over all delay levels
    for (i = 0; i < m_max_comb_delay; i++) 
	{
		level_node = m_cluster->Level[i];
		assert(level_node);

		number_zero_degrees_possible = level_node->get_nLatched() + level_node->get_nPO();

		// add a penalty if we have too many zero degree nodes
		if (level_node->OutDegrees(0) > number_zero_degrees_possible)
		{
			cost_Fanout_Penalty = static_cast<double>(level_node->OutDegrees(0) - number_zero_degrees_possible);
		}
	
		// add a penalty if the maximum degree at this level is greater that 
		// the maxdegree possible for this level 
		if (m_current_maxdeg[i] > m_maxdeg[i])
		{
			assert(m_current_maxdeg[i] == level_node->get_max_degree_assigned());
			//Warning("Our max degree is greater than possible for this level");
			cost_Fanout_Penalty  += static_cast<double>(m_current_maxdeg[i] - m_maxdeg[i]);
		}
	}
	
	return cost_Fanout_Penalty;
}


COST_TYPE DEGREE_PARTITIONER::get_cost_degree_fanout_edge_misassignment() const
{
	DELAY_TYPE i;
	COST_TYPE cost_Degree_Fanout_Edge_Misassignment = 0.0;

	// sum up the cost over all delay levels
    for (i = 0; i < m_max_comb_delay; i++) 
	{

		// cost_Degree_Fanout_Edge_Misassignment
		// first part of the cost is the difference between output degree of the level node
		// and the sum of the degrees assigned to the level node
		cost_Degree_Fanout_Edge_Misassignment +=	static_cast<double>(ABS(m_target[i] - m_sum_assigned[i]));
	}

	return cost_Degree_Fanout_Edge_Misassignment;
}

// generate a move in the iterative algorithm
// 
// The degrees to be swapped are chosen by first, 
// randomly choosing a non-zero fanout degree from the level nodes 
// with higher total fanout than number of output edges assigned. 
// Next, we randomly choose a non-zero lower fanout degree from the 
// level nodes with lower total fanout than the number of output edges assigned. 
// Zero fanout degrees are never chosen because they can only be assigned to either 
// a latched node or a node with a primary output and hence their movement is 
// more complicated. 
//
void DEGREE_PARTITIONER::generate_move()
{
	m_move.reset();

	DELAY_TYPE source_level = 0; 
	DELAY_TYPE dest_level = 0;
	NUM_ELEMENTS source_degree = 0; 
	NUM_ELEMENTS dest_degree = 0;
	NUM_ELEMENTS_VECTOR source_level_pmf(m_max_comb_delay,0);
	NUM_ELEMENTS_VECTOR dest_level_pmf(m_max_comb_delay,0);
	DELAY_TYPE delay = 0;
	NUM_ELEMENTS difference = 0;

	for (delay = 0; delay < m_max_comb_delay; delay++)
	{
		difference = m_target[delay] - m_sum_assigned[delay];

		// pos differences means the level node is underassigned so make them sinks
		// neg differences means the level node is overassigned so make them sources
		if (difference > 0)
		{
			dest_level_pmf[delay] = difference;
		}
		else
		{
			source_level_pmf[delay] = -difference;
		}
	}

	// hack.  i think with available instead of assigned in 
	// we sometimes get degrees that don't sum to zero
	
	if (accumulate(source_level_pmf.begin(), source_level_pmf.end(), 0) == 0 ||
		accumulate(dest_level_pmf.begin(), dest_level_pmf.end(), 0) == 0)
	{
		m_move.set_invalid();
		return;
	}
		
	source_level = (short) g_rand_number_gen->discrete_pmf(m_max_comb_delay, source_level_pmf);
	dest_level = (short) g_rand_number_gen->discrete_pmf(m_max_comb_delay, dest_level_pmf);

	//debugif(DDegree,"source level " << source_level << " sink level " << dest_level);

	// impossible for both levels to be the same
	assert(source_level != dest_level);

	m_move.set_source_level(source_level);
	m_move.set_dest_level(dest_level);

	source_degree = choose_source_degree(source_level);
	dest_degree = choose_dest_degree(dest_level, source_degree);

	m_move.set_source_degree(source_degree);
	m_move.set_dest_degree(dest_degree);


	assert(m_cluster->Level[source_level]->OutDegrees(source_degree) > 0);
}

DEGREE_MOVE::DEGREE_MOVE()
{
	m_source_level = 0;
	m_dest_level = 0;
	m_source_degree = 0;
	m_dest_degree = 0;

	m_valid = true;
}

DEGREE_MOVE::~DEGREE_MOVE()
{
}

	
// 0 degree moves are more complicated because they need
// to be in level nodes with PO and Latches.
//
// For now, bar their movement.
//
// RETURNS: if the move is acceptable
//
bool DEGREE_MOVE::is_valid() const
{
	// just don't let 0 degree nodes move for right now.
	// think about how we want to let latch and PO with 0 degree
	// nodes move later.
	bool valid_move = m_valid && (m_source_degree != 0 && m_dest_degree != 0);

	return valid_move;
}

// make the move
//
// PRE: move is valid
// POST: the move has been made
//
void DEGREE_PARTITIONER::make_move()
{
	DELAY_TYPE  source_level = m_move.get_source_level();
	DELAY_TYPE  dest_level = m_move.get_dest_level();
	NUM_ELEMENTS source_degree = m_move.get_source_degree();
	NUM_ELEMENTS dest_degree = m_move.get_dest_degree();

	assign_to_level(source_degree, source_level, -1);
	assign_to_level(dest_degree, source_level, +1);
	assign_to_level(dest_degree, dest_level, -1);
	assign_to_level(source_degree, dest_level,+1);

	m_current_maxdeg[source_level] = m_cluster->Level[source_level]->get_max_degree_assigned();
	m_current_maxdeg[dest_level] = m_cluster->Level[dest_level]->get_max_degree_assigned();
	assert(m_current_maxdeg[source_level] == m_cluster->Level[source_level]->get_max_degree_assigned());
	assert(m_current_maxdeg[dest_level] == m_cluster->Level[dest_level]->get_max_degree_assigned());
	
	m_move.reset();
}


// choose a source degree for the move from the source delay level
// a source degree can be anything that is greater than 0
//  
// RETURN: a source degree
//
NUM_ELEMENTS DEGREE_PARTITIONER::choose_source_degree
(
	const DELAY_TYPE& source_level
) const
{
	LEVEL_NODE * level_node = m_cluster->Level[source_level];
	assert(level_node);

	DISTRIBUTION fanout_distribution = level_node->get_fanout_distribution();
	NUM_ELEMENTS dist_size = fanout_distribution.size();
	NUM_ELEMENTS i;
	NUM_ELEMENTS source_degree = 0;

	for (i = 0; i < dist_size; i++)
	{
		assert(fanout_distribution[i] >= 0);

		if (fanout_distribution[i] > 0)
		{
			fanout_distribution[i] = 1;
		}
	}

	assert(accumulate(fanout_distribution.begin(), fanout_distribution.end(), 0) > 0);

	source_degree = g_rand_number_gen->discrete_pmf(dist_size, fanout_distribution);

	return source_degree;
}

// choose a destination degree for the move from the destination delay level
//  
// PRE: dest_level is a delay level with fanouts
//      m_target[delay] - m_sum_assigned[delay] > 0
// RETURN: a destination degree
NUM_ELEMENTS DEGREE_PARTITIONER::choose_dest_degree
(
	const DELAY_TYPE& dest_level,
	const NUM_ELEMENTS& source_degree
)
{
	LEVEL_NODE * level_node = m_cluster->Level[dest_level];
	assert(level_node);

	NUM_ELEMENTS_VECTOR fanout_distribution = level_node->get_fanout_distribution();
	NUM_ELEMENTS dist_size = fanout_distribution.size();
	NUM_ELEMENTS i;
	NUM_ELEMENTS dest_degree = 0;
	NUM_ELEMENTS tot_under_source_degree = 0;
	NUM_ELEMENTS tot = 0;
	

	// give every out degree less than source degree (except 0) an equal chance of being selected
	for (i = 0; i < dist_size; i++)
	{
		assert(fanout_distribution[i] >= 0);

		if (fanout_distribution[i] > 0)
		{
			fanout_distribution[i] = 1;
		}

		if (i < source_degree)
		{
			tot_under_source_degree += fanout_distribution[i];
		}
		tot += fanout_distribution[i];
	}

	// if we have fanout values under the value of the source degree
	// choose among them
	//
	if (tot_under_source_degree >  0)
	{
		dest_degree = g_rand_number_gen->discrete_pmf(source_degree, fanout_distribution);
	}
	else if (tot > 0)
	{
		// we don't have any fanouts less than the source degree but we do
		// have fanouts. choose one
		//
		dest_degree = g_rand_number_gen->discrete_pmf(dist_size, fanout_distribution);
	}
	else
	{
		// no fanout above zero has been assigned to this level node.
		//
		// we should at least have zero degree fanouts or else 
		// we should not have selected the destination delay level
		assert(level_node->OutDegrees(0) > 0);
		dest_degree = 0;
	}

	return dest_degree;
}


// Our level nodes now have super edges and a degree assignement.
// We need to check to see if the sum of the output degrees is equal to 
// the weight of the out edges.
//
// If not....
//
// we have a choice to make.  We can either 
// move edges around in the comb graph or change degree assignments in order
// to fix the degree assignments.
//
// I've choosen to change the degree assignments.
// More experimentation needed.
//
// if a level node has an overassignment of total output degree to edges
// 	we will preferentially select higher fanouts first 
//
// if a level node has an underassignment of total output degree to edges
// 	we will preferentially select lower fanouts first
// 
//
// POST: for each level node in the cluster the sum of the fanout 
//       degrees is equal to the number of output edges
//
//       if we had over assigned level nodes we may have had to,
//       in desparate circumstances, had to add a primary output
//       or a length 1 intra-cluster edge
//
void DEGREE_PARTITIONER::fix_degree_assignments()
{
	DELAY_TYPE lev;
	assert(get_cost() > 0);
	LEVEL_NODE * level_node = 0;
	DISTRIBUTION fanout_distribution;
	NUM_ELEMENTS double_connection_weight = 0,
				 fanout_size = 0,
				 degree	= 0,
				 choice = 0,
				 total = 0,
				 max_weight = 0,
			     edge_to_fanout_difference = 0,
				 nOver_assigned = 0,
				 nUnder_assigned = 0,
				 nPossible = 0;
	EDGE * edge = 0;
	bool ok = true,
		have_open_inputs,
		have_double_connections;


    for (lev = 0; lev < m_max_comb_delay; lev++) 
	{
		level_node =  m_cluster->Level[lev];
		assert(level_node);

		edge_to_fanout_difference = m_target[lev] - m_sum_assigned[lev];
		

		if (edge_to_fanout_difference > 0)
		{
			// we have a level nodes with more edges that fanouts
			// increase the fanouts
			nUnder_assigned = edge_to_fanout_difference;

			debugif(DDegree,"We are under assigned by " << nUnder_assigned << " fanouts at level " << lev);

			fanout_distribution = level_node->get_fanout_distribution();
			fanout_size = static_cast<NUM_ELEMENTS>(fanout_distribution.size());
		
			// try not to touch the zero degree nodes
			// because they are PO and latched nodes which do 
			// not in generally fanout to any combinational nodes
			fanout_distribution[0] = 0;

			if (DDegree)
			{
				debug("Original fanout dist");
				copy(fanout_distribution.begin(), fanout_distribution.end(), 
						ostream_iterator<NUM_ELEMENTS>(cout, " ")); cout << endl;
				debugsep;
			}
		
			// while we still have edge-fanout mismatched
			while (nUnder_assigned > 0)
			{
				// choose a fanout and increase its degree by 1
				if (accumulate(fanout_distribution.begin(), fanout_distribution.end(), 0) == 0)
				{
					// we need to increase a zero degree assignment
					assert(level_node->OutDegrees(0) > 0);
					choice = 0;
				}
				else
				{
					choice = g_rand_number_gen->discrete_pmf(fanout_size, fanout_distribution);
					assert(choice >= 1);
				}
				debugif(DDegree,"Decreasing degree " << choice << " increasing degree " << choice +1);
				assert(level_node->OutDegrees(choice) > 0);

				// resize the fnaout distributions if need be
				if (choice >= fanout_size-1)
				{
					m_cluster->set_new_max_fanout(choice+1);
					fanout_size = choice+2;
					m_fanout_distribution.resize(fanout_size, 0);	
					fanout_distribution.resize(fanout_size, 0);
				}

				// change the degree assignments
				assign_to_level(choice, lev, -1);
				assign_to_level(choice+1, lev, +1);
				assert(level_node->OutDegrees(choice) >= 0);
				assert(level_node->OutDegrees(choice+1) > 0);

				if (choice != 0)
				{
					fanout_distribution[choice] = fanout_distribution[choice]-1;
				}
				fanout_distribution[choice+1] += 1;
				assert(fanout_distribution[choice] >= 0);
				assert(fanout_distribution[choice+1] > 0);

				nUnder_assigned--;
				
				if(DDegree)
				{
					debug("after correction");
					copy(fanout_distribution.begin(), fanout_distribution.end(), 
							ostream_iterator<NUM_ELEMENTS>(cout, " ")); cout << endl;
					debug("Number under assigned now " << nUnder_assigned);
				}
			}
			fanout_distribution = level_node->get_fanout_distribution();
			if(DDegree)
			{
				debug("level now correct");
				debug("fixed fanout dist:");
				copy(fanout_distribution.begin(), fanout_distribution.end(), 
						ostream_iterator<NUM_ELEMENTS>(cout, " ")); cout << endl;
				debugsep;
			}

		}
		else if (edge_to_fanout_difference < 0)
		{
			// we have a level nodes with fanouts than edges 
			// degrees the fanouts
			
			nOver_assigned = -edge_to_fanout_difference;
			debugif(DDegree,"We are over assigned by " << nOver_assigned << " fanouts at level " << lev);

			fanout_distribution = level_node->get_fanout_distribution();
			fanout_size = static_cast<NUM_ELEMENTS>(fanout_distribution.size());
			if(DDegree)
			{
				debug("Original fanout dist");
				copy(fanout_distribution.begin(), fanout_distribution.end(), 
						ostream_iterator<NUM_ELEMENTS>(cout, " ")); cout << endl;
				debugsep;
			}


			// don't touch the zero or 1 degree nodes
			// a zero degree node means a latched node or a PO
			// and we don't want to add more of these
			fanout_distribution[0] = 0;
			fanout_distribution[1] = 0;
			total += 0;

			// spike the distribution towards high fanout
			for (degree = 1; degree < fanout_size; degree++)
			{
				fanout_distribution[degree] = fanout_distribution[degree] * degree;
			}

			ok = true;
			// while we have more fanout degree than edges
			while (nOver_assigned > 0)
			{
				// if we have no other fanouts but zero and 1 degrees
				if (accumulate(fanout_distribution.begin(), fanout_distribution.end(), 0) == 0)
				{
					// we have no degrees above 1
					// we should still have out degrees 1, right? 
					assert(level_node->OutDegrees(1) > 0);

					// do we still have latches or PO we can assign to have zero fanout degree
					if (level_node->get_nLatched() + level_node->get_nPO() > level_node->OutDegrees(0) ||
						m_cluster->Level[lev+1]->get_nNodes() == 0)
					{
						choice = 1;
					}
					else
					{
						// we are screwed.  we need to steal a PO from somewhere else OR we need to add
						// an edge. lets add an edge. increment length 1 intra cluster edge by nOver_assigned
						// hack

						edge = m_cluster->Level[lev+1]->find_edge_with_source(level_node);

						double_connection_weight 
							= m_cluster->Level[lev]->get_nNodes() * 
							m_cluster->Level[lev+1]->get_nNodes();

						have_open_inputs = (m_cluster->Level[lev+1]->get_nOpen_inputs() > 1);
						have_double_connections = (edge->get_weight() + 1 > double_connection_weight);

						if (! have_double_connections && have_open_inputs)
						{
							nPossible =  MIN(nOver_assigned, m_cluster->Level[lev+1]->get_nOpen_inputs());
							nPossible = MIN(nOver_assigned, double_connection_weight);
							max_weight = MAX(edge->get_max_weight(), edge->get_weight() + nPossible);

							edge->set_max_weight(max_weight);

							Warning("We need to add " << nPossible << " edges to cluster " << 
									m_cluster->get_cluster_number() << " at source delay level " << lev);
						
							edge->add_weight(nPossible);
							nOver_assigned -= nPossible;
							ok = false;
						}
						else
						{
							// we are really screwed.
							Warning("Look at section. Sigh. In degree assignment.  Had to add " << nOver_assigned << " POs");
							choice = 1;
							m_cluster->Level[lev]->add_PO(1);
						}
					}
				}
				else
				{
					choice = g_rand_number_gen->discrete_pmf(fanout_size, fanout_distribution);
					assert(choice > 1);
				}

				if (ok)
				{
					debugif(DDegree,"Decreasing degree " << choice << " increasing degree " << choice -1);
					assert(level_node->OutDegrees(choice) > 0);


					assign_to_level(choice, lev, -1);
					assign_to_level(choice-1, lev, +1);
					assert(level_node->OutDegrees(choice) >= 0);
					assert(level_node->OutDegrees(choice-1) > 0);

					if (choice > 1)
					{
						// increase by the new weight

						if (choice > 2)
						{
							fanout_distribution[choice-1] += choice-1;
							assert(fanout_distribution[choice-1] > 0);
						}
						fanout_distribution[choice] = fanout_distribution[choice]-choice;
						assert(fanout_distribution[choice] >= 0);
					}
					nOver_assigned--;
					
					if (DDegree)
					{
						debug("after correction");
						copy(fanout_distribution.begin(), fanout_distribution.end(), 
								ostream_iterator<NUM_ELEMENTS>(cout, " ")); cout << endl;
						debugif(DDegree,"Number overassigned now " << nOver_assigned);
					}
				}
			}
			fanout_distribution = level_node->get_fanout_distribution();
			if(DDegree)
			{
				debug("level now correct");
				debug("fixed fanout dist");
				copy(fanout_distribution.begin(), fanout_distribution.end(), 
						ostream_iterator<NUM_ELEMENTS>(cout, " ")); cout << endl;
				debugsep;
			}
		}
		else
		{
			// do nothing. we are ok
		}
		assert(m_target[lev] == m_sum_assigned[lev]);
		m_current_maxdeg[lev] = level_node->get_max_degree_assigned();
	}
}
	

//
// Sometimes we will have assigned a degree or two that is 
//
// PRE: m_maxdeg is the maximum degree possible without node violations
// POST: all nodes in the cluster have a fanout that is less than m_maxdeg for
//       their delay level OR else we have failed and exited the program
//
void DEGREE_PARTITIONER::fix_max_degree_assignments()
{
	DELAY_TYPE lev;
	LEVEL_NODE * level_node = 0;
	DISTRIBUTION fanout_distribution;
	NUM_ELEMENTS fanout_size,
					degree	= 0,
					choice = 0,
					max_degree_assigned = 0;

    for (lev = 0; lev < m_max_comb_delay; lev++) 
	{
		level_node =  m_cluster->Level[lev];
		assert(level_node);

		if (level_node->get_max_degree_assigned() > m_maxdeg[lev])
		{
			debugif(DDegree,"We have a max poss degree of " << m_maxdeg[lev] << " at level " << lev << 
				" but we have a assigned this level a degree of " << level_node->get_max_degree_assigned());
			debugif(DDegree,"Trying to fix this.");

			fanout_distribution = level_node->get_fanout_distribution();
			fanout_size = static_cast<NUM_ELEMENTS>(fanout_distribution.size());
			if(DDegree)
			{
				debug("Original fanout dist");
				copy(fanout_distribution.begin(), fanout_distribution.end(), 
						ostream_iterator<NUM_ELEMENTS>(cout, " ")); cout << endl;
				debugsep;
			}


			// spike the sampling distribution towards high fanout. 
			// eliminate any degree more than maxdeg
			for (degree = 1; degree < fanout_size; degree++)
			{
				fanout_distribution[degree] = fanout_distribution[degree] * degree;

				if (degree >= m_maxdeg[lev])
				{
					fanout_distribution[degree] = 0;
				}
			}


			max_degree_assigned = level_node->get_max_degree_assigned();

			while (max_degree_assigned > m_maxdeg[lev])
			{
				debugif(DDegree,"new fanout distribution");

				if (DDegree)
				{
					copy(fanout_distribution.begin(), fanout_distribution.end(), 
						ostream_iterator<NUM_ELEMENTS>(cout, " ")); cout << endl;
				}

				if (accumulate(fanout_distribution.begin(), 
							fanout_distribution.end(), 0) == 0)
				{
					// we don't have any degrees under maxdeg
					// split the max_degree_assign
				
					level_node->print_output_connections();
					Fail("We have a double input where we shouldn't");
				}
				// pick a degree less than m_maxdeg[lev] to increase
				choice = g_rand_number_gen->discrete_pmf(fanout_size, fanout_distribution);
				assert(choice >= 0 && choice < m_maxdeg[lev]);

				debugif(DDegree,"Decreasing degree " << max_degree_assigned << " increasing degree " << choice);
				assert(level_node->OutDegrees(choice) > 0);

				// decrease the max_degree assigned and increase the choice 
				assign_to_level(max_degree_assigned, lev, -1);
				assign_to_level(max_degree_assigned-1, lev, +1);
				assign_to_level(choice, lev, -1);
				assign_to_level(choice+1, lev, +1);

				// max sure everything is ok with our level node
				assert(level_node->OutDegrees(max_degree_assigned) >= 0);
				assert(level_node->OutDegrees(choice) >= 0);
				assert(level_node->OutDegrees(max_degree_assigned-1) > 0);
				assert(level_node->OutDegrees(choice+1) > 0);

				// fix sampling distribution
				if (choice > 0)
				{
					fanout_distribution[choice] -= choice;
					assert(fanout_distribution[choice] >= 0);
					
				}
				else
				{
					fanout_distribution[0] -= 1;
					assert(fanout_distribution[0] >= 0);
				}

				// if we can still increase the degree at (choice+1) put it back into the
				// fanout distribution for sampling
				if (m_maxdeg[lev] > choice+1)
				{
					fanout_distribution[choice+1] += choice+1;
				}

				// update max_degree assigned
				max_degree_assigned = level_node->get_max_degree_assigned();
			}
			debugif(DDegree,"level now correct");
		}
		else
		{
			// do nothing. we are ok
		}
		assert(m_target[lev] == m_sum_assigned[lev]);
		m_current_maxdeg[lev] = max_degree_assigned;
	}
}


//
// decrease m_max_fanout until we have a degree to assign
// 
// POST: the maximum fanout left to be assigned is m_max_fanout
void DEGREE_PARTITIONER::recalculate_max_fanout()
{
	while (m_fanout_distribution[m_max_fanout] == 0 && m_max_fanout > 1)
	{
		m_max_fanout -= 1;
		debugif(DDegree,"..Decreasing m_max_fanout to " << m_max_fanout);
	}
}
	

DEGREE_PARTITIONER::~DEGREE_PARTITIONER()
{
}


void DEGREE_PARTITIONER::prevent_double_connections()
{
	DELAY_TYPE lev;
	LEVEL_NODE * level_node = 0;
	NUM_ELEMENTS nAdditional_long_edges_required = 0;
	NUM_ELEMENTS degree;
	NUM_ELEMENTS replacement_degree = 0;
	NUM_ELEMENTS nNodes_at_next_level = 0;
	NUM_ELEMENTS choice = 0;
	DISTRIBUTION fanout_degrees_that_need_long_edges;
	DISTRIBUTION fanout_distribution;
	NUM_ELEMENTS dist_size;

    for (lev = 0; lev < m_max_comb_delay; lev++) 
	{
		level_node =  m_cluster->Level[lev];
		assert(level_node);

		if (level_node->get_available_out() > 0)
		{
			find_out_how_many_additional_long_edges_required(level_node, 
					nAdditional_long_edges_required, fanout_degrees_that_need_long_edges);

			// while we still need long edges reduce the degree on 
			// the high fanout nodes

			assert(new_gen); // right thing to do?

			while (nAdditional_long_edges_required > 0)
			{
				// we do not have enough long edges because some nodes have too high a fanout for the number
				// of nodes at the next level
				// we can either covert some of the length 1 intra cluster edges to longs OR
				// we can just reduce the degree on the nodes
				// I am going to reduce the degrees on the nodes

				debugif(DDegree,"We are using " << nAdditional_long_edges_required 
						<< " too many long or inter cluster edges.");
				nNodes_at_next_level = m_cluster->Level[lev+1]->get_nNodes();
				assert(nNodes_at_next_level > 0);

				fanout_distribution = level_node->get_fanout_distribution();
				dist_size = fanout_distribution.size();

				choice = g_rand_number_gen->random_number(fanout_degrees_that_need_long_edges.size() - 1);

				degree = fanout_degrees_that_need_long_edges[choice];
				assert(degree > 1 && degree < dist_size);
				assert(fanout_distribution[degree] > 0);
				assert(degree - nNodes_at_next_level > 0);

				// look for a degree to increase below that of the degree choosen for reduction
				replacement_degree = g_rand_number_gen->discrete_pmf(degree-1, fanout_distribution);

				assert(new_gen);	// should we try and not increase the zero degree nodes?
				assert(replacement_degree >= 0 && replacement_degree < degree);
				assert(fanout_distribution[replacement_degree] > 0);

				assign_to_level(degree, lev, -1);
				assign_to_level(degree-1, lev, +1);
				assign_to_level(replacement_degree, lev, -1);
				assign_to_level(replacement_degree+1, lev, +1);

				find_out_how_many_additional_long_edges_required(level_node, 
						nAdditional_long_edges_required, fanout_degrees_that_need_long_edges);
			}
		}
	}
}

// 
// PRE: level node is valid
// 		the delay level of the level node is < delay of the circuit
//
// POST: nAdditional_long_edges_required is the difference between the
//       number of long edges that the level nodes needs and the number of long
//       edges that have been assigned
//
//       fanout_degrees_that_need_long_edges contains the fanout degrees that need 
//       long edges
//
void DEGREE_PARTITIONER::find_out_how_many_additional_long_edges_required
(
	const LEVEL_NODE * level_node,
	NUM_ELEMENTS& nAdditional_long_edges_required,
	DISTRIBUTION& fanout_degrees_that_need_long_edges
) const
{
	assert(level_node);

	fanout_degrees_that_need_long_edges.clear();

	NUM_ELEMENTS degree = 0;
	DISTRIBUTION fanout_distribution = level_node->get_fanout_distribution();
	NUM_ELEMENTS dist_size = fanout_distribution.size();
	DELAY_TYPE   delay = level_node->get_delay_level();
	NUM_ELEMENTS number_of_long_edges_needed = 0,
	             nLong_edges_have = 0,
				 total_number_of_long_edges_needed = 0;

	assert(delay < m_max_comb_delay);
	NUM_ELEMENTS nNodes_at_next_level = m_cluster->Level[delay+1]->get_nNodes();

	// find out how many long edges will be needed by the individual nodes
	// in this level node
	// a long edge is an intra cluster edge of length > 1 or an inter-cluster edge
	//
	// long edges are needed if a degree assigned to the level node is more than the 
	// number of individual nodes at the next delay level 
	for (degree = 1; degree < dist_size; degree++)
	{
		if (fanout_distribution[degree] > 0)
		{
			// the number of long edges needed for this fanout =
			// (the number of nodes with this fanout * the number of longs edges that will be 
			// needed by this fanout)
			//
			// the number of longs edges that will be needed by this fanout = 
			// (degree - the number of nodes at the next level)
			//
			number_of_long_edges_needed = fanout_distribution[degree]* (degree - nNodes_at_next_level);

			if (number_of_long_edges_needed > 0)
			{
				total_number_of_long_edges_needed += number_of_long_edges_needed;

				fanout_degrees_that_need_long_edges.push_back(degree);
			}
		}
	}

	level_node->sanity_check();

	assert(level_node->get_intra_cluster_output_edge_lengths().size() >= 2);
	assert(level_node->get_intra_cluster_output_edge_lengths()[1] >= 0);

	// find out how many long edges to we have
	nLong_edges_have = level_node->get_available_out() - 
						level_node->get_intra_cluster_output_edge_lengths()[1];


	nAdditional_long_edges_required = MAX(0, total_number_of_long_edges_needed - nLong_edges_have);
}

// POST: m_maxdeg has been reset for all delay levels 
//         to the maximum degree possible without node violations
//		 m_current_maxdeg has been reset for all delay levels
//		   to the maximum degree assigned
//
void DEGREE_PARTITIONER::recalculate_max_degree()
{
	DELAY_TYPE lev;
	NUM_ELEMENTS maxdegree_without_double_inputs = 0;
	LEVEL_NODE * level_node = 0;

    for (lev = 0; lev < m_max_comb_delay; lev++) 
	{
		level_node = m_cluster->Level[lev];
		assert(level_node);

		maxdegree_without_double_inputs = level_node->get_max_degree_possible();
	
		m_maxdeg[lev] = maxdegree_without_double_inputs;
		m_current_maxdeg[lev] = level_node->get_max_degree_assigned();
	}
}			

// POST: we have reported to the user the quality of the solution
//
NUM_ELEMENTS DEGREE_PARTITIONER::report_solution_quality() const
{

	// only currently reporting deviation from specification
	// might want to add the cost_Degree_Fanout_Edge_Misassignment
	// even though it is not quite related to spec
	assert(new_gen);



	assert(m_cluster);
	CLUSTER_SPEC * comb_spec = m_cluster->get_comb_spec();
	assert(comb_spec);
	
	DELAY_TYPE lev;
	LEVEL_NODE * level_node = 0;
	NUM_ELEMENTS deg = 0;
	NUM_ELEMENTS max_fanout = m_cluster->get_max_fanout();

	DISTRIBUTION degree_distribution(max_fanout + 1, 0);
	DISTRIBUTION spec_degree_distribution = comb_spec->get_fanout_distribution();
				 spec_degree_distribution.resize(max_fanout +1, 0);
	NUM_ELEMENTS degree_difference = 0;
	NUM_ELEMENTS spec_max_fanout_degree = max_fanout,
					max_fanout_degree = max_fanout;

	assert(spec_degree_distribution.size() == degree_distribution.size());


	// find out the degree distribution
    for (lev = 0; lev < m_max_comb_delay; lev++) 
	{
		level_node = m_cluster->Level[lev];
		assert(level_node);
		
		for (deg = 0; deg <= max_fanout; deg++)
		{
			degree_distribution[deg] += level_node->OutDegrees(deg);
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

	debug("\nDifference between Fanout Distribution and spec: " << degree_difference);

	return degree_difference;
}
//
// POST: we have reset the move
void DEGREE_MOVE::reset()
{
	m_source_level = 0;
	m_dest_level = 0;
	m_source_degree = 0;
	m_dest_degree = 0;
	m_valid = true;
}
