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

//
// Class_name NODE_SPLITTER
//
// Description
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

class NODE_SPLITTER
{
public:
	NODE_SPLITTER();
	~NODE_SPLITTER();
	NUM_ELEMENTS split_nodes(CLUSTER * cluster);

private:
	CLUSTER * m_cluster;

	void assign_nodes_a_horizontal_position();
	void assign_fanouts(const DELAY_TYPE& lev);
	void assign_high_fanouts(const DELAY_TYPE& lev);
	DISTRIBUTION tree_pmf(const int& n, int mult);
	void assign_low_fanouts(const DELAY_TYPE& lev);
	void designate_long_and_inter_cluster_edges(const DELAY_TYPE& delay_level);
	void designate_latched_nodes(const DELAY_TYPE& delay_level);
	void designate_dff_and_pi_nodes();
	void randomly_designate_dff_and_pi_nodes();
	void designate_PO();
	void designate_node_as_dff(NODE * node);
	void designate_node_as_pi(NODE * node);
	void display();
	void display_fanout_results(const DELAY_TYPE& lev) const;
	void get_pi_and_dff_statistics(double& pi_mean, double& dff_mean, double& pi_std_dev, double& dff_std_dev,
									double& high_degree_boundry, double& low_degree_boundry);
	void deal_with_weird_special_case(const double& fanout_degree, const double& pi_mean, 
										const double& pi_std_dev, const double& dff_mean, 
										const double& dff_std_dev, const double& high_degree_boundry,
										const double& low_degree_boundry, NUM_ELEMENTS& nDFF, 
										NUM_ELEMENTS& nPI, NODE * node);
	void final_sanity_check() const;
	NUM_ELEMENTS report_solution_quality() const;
};

/*  old header:
 *
 *  First part of the nodes splitting step.  Allocate the nodes for each
 *  level, and assign specific out-degrees to them from level_node->OutDegrees,
 *  which is the integer partition of level_node->assigned_out.
 *  
 *  For this step, the levels are treated independently.  We need to
 *  split each level-node $N_i$ into $n_i$ individual nodes, and assign
 *  each of these a fanout from the list of available fanouts for $N_i$.
 *  This would be trivial, were it not for the necessity to introduce
 *  {\em locality} into the resulting circuit.
 *   
 *  Because of the way that real circuits are designed, an inherent
 *  local structure develops in the circuit.  Nodes tend to exhibit a
 *  clustered behaviour, whereby nodes in a cluster tend to accept fanin
 *  from approximately the same set of nodes as other nodes in their
 *  cluster.  This local clustering is described by Rent's Rule
 *  \cite{LandmanRusso} and theoretical models to describe it have been
 *  proposed by Donath \cite{Donath} and others.  Without some method
 *  of modeling local behaviour, our circuits would be ``too random" and
 *  hence too difficult for CAD tools.
 *   
 *  Our approach to introducing locality into the generation algorithm is
 *  to define a geometric proximity for nodes so that, when we choose
 *  edge-connections between nodes in the next step, we can use this
 *  metric to choose between otherwise identical nodes.  This can also
 *  be viewed as trying to generate graphs which will look good when
 *  displayed as pictures such as Figure~\ref{complete}.  The geometric
 *  proximity we use is the sorted order within the linear list of nodes
 *  within each level,  The measure of ``goodness" of an
 *  edge is then measured as the distance between the source and
 *  destination nodes in their levels node-lists, relative to
 *  competitors.  As a result, the order in which fanouts are assigned
 *  within the node list becomes important, because placing high-fanout
 *  nodes in an unbalanced way into the node list will skew the effects
 *  of locality measurement in Step 5.
 *   
 *  The locality index assigned to each of the $n_i$ nodes in the nodelist
 *  for level $i$ is a scaled proportion in the maximum sized level.  Thus
 *  if the maximum level contains 100 nodes, and the current level 10, then
 *  its nodes will have locality indices 5, 15, 25, ..., 95.
 *  
 *  Our goal in assigning fanouts to nodes in the list is to distribute the
 *  high fanout nodes well for maximum locality generation.  To do this,
 *  we sample a ``binary tree distribution" to allocate fanouts, in order
 *  from the highest to lowest fanout.  For x[] the inorder traversal of
 *  a complete binary tree, define the value of x[j], j=1..n, as the number
 *  of nodes beneath $x[j]$ in the tree.  For example, the binary tree pmf of
 *  length 16 is [1,2,1,4,1,2,1,8,1,2,1,4,1,2,1,2].   In the most likely case,
 *  then, the highest fanout node would be assigned in the middle, the next
 *  two highest fanouts at the quartiles, and so on.  The probabilistic
 *  sampling means we don't get exactly the same result each time.
 *   
 *  Once the fanouts have been assigned, we allocate the $p_i$ long
 *  (non-unit) edges between nodes among the $f_i$ total fanout edges
 *  represented by the node fanouts.  This is done probabilistically
 *  by calculating a $p$-subset of $f_i$, and assigning linearly through
 *  the fanouts in the node list.
 *   
 *  At the conclusion of Step 4, each node in level $i$ has an assigned
 *  fanout $f_j$, and a distribution of the number of each length which
 *  are represented by $f_i$, but no actual edges have been assigned
 *  between nodes at different levels in the graph.
 *
 *  Historical note:  JP had a more complicated way of doing this, which 
 *  I re-wrote in July 96.  M. Hutton.
 */
