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




#include "graph.h"
#include <algorithm>
#include <numeric>
#include "gen.h"
#include "degree.h"
#include "node_splitter.h"
#include "simple_edge_assignment.h"

CLUSTER::CLUSTER
(
	CLUSTER_SPEC * comb_spec,
	const CLUSTER_NUMBER_TYPE & cluster_number,
	const NUM_ELEMENTS& nClusters,
	CIRCUIT_GRAPH * circuit
)
{
	assert(circuit);
	m_my_circuit = circuit;

	DELAY_TYPE	level = 0;
	LEVEL_NODE*	level_node = 0;

	m_comb_spec = comb_spec;
	m_cluster_number = cluster_number;
	m_nNodes 		= m_comb_spec->get_nNodes();
	m_nDFF			= m_comb_spec->get_nDFF();
	m_nPI			= m_comb_spec->get_nPI();
	m_nPO			= m_comb_spec->get_nPO();
	m_nLatched		= m_comb_spec->get_nLatched();
	m_delay			= m_comb_spec->get_delay();
	m_kin			= m_comb_spec->get_kin();
	m_max_fanout	= m_comb_spec->get_max_fanout();

	for (level= 0; level <= m_delay; level++)
	{
		level_node = new LEVEL_NODE(level, m_comb_spec, nClusters, this);
		Level.push_back(level_node);
	}
	m_intra_cluster_edge_lengths.resize(m_delay+1, 0);
	m_inter_cluster_input_edge_lengths.resize(m_delay+1, 0);
	m_inter_cluster_output_edge_lengths.resize(m_delay+1, 0);

	m_nIntra_cluster_edges = 0;
	m_nInter_cluster_input_edges = 0;
	m_nInter_cluster_output_edges = 0;

	SHAPE input_shape = m_comb_spec->get_level_in();
	SHAPE output_shape = m_comb_spec->get_level_out();
	NUM_ELEMENTS level_node_normalization_cost;
	level_node_normalization_cost = accumulate(input_shape.begin(), input_shape.end(), 0);
	level_node_normalization_cost += accumulate(output_shape.begin(), output_shape.end(), 0);
	m_level_node_normalization_cost = static_cast<COST_TYPE>(level_node_normalization_cost);
	m_wirelength_wanted				= m_comb_spec->get_wirelength();
}

CLUSTER::~CLUSTER()
{
	LEVEL_NODES::iterator level_node_iter;
	LEVEL_NODE * level_node = 0;

	for (level_node_iter = Level.begin(); level_node_iter != Level.end(); level_node_iter++)
	{
		level_node = *level_node_iter;
	
		delete level_node;
	}
}

// RETURN: true if the all specification are met for the combinational delay structure
//         false otherwise
bool CLUSTER::is_valid_solution() const
{
	bool valid_solution =
		is_intra_cluster_edge_length_distribution_satisfied() &&
		is_inter_cluster_input_edge_length_distribution() &&
		is_inter_cluster_output_edge_length_distribution() &&
		are_level_node_properties_satisfied();
	
	return valid_solution;
}
// RETURN: true if the the spec for the inter-cluster input and output edge length 
//         distributions are satisfied, false otherwise
bool CLUSTER::is_inter_cluster_edge_length_distributions_satisfied() const
{
	bool valid_solution =
		is_inter_cluster_input_edge_length_distribution() &&
		is_inter_cluster_output_edge_length_distribution();
	
	return valid_solution;
}
// RETURN: true if the the spec for the intra-cluster edge length 
//         distributions is satisfied, false otherwise
bool CLUSTER::is_intra_cluster_edge_length_distribution_satisfied() const
{
	bool edge_length_distribution_satisfied = false;

	const DISTRIBUTION edge_length_distribution = get_intra_cluster_edge_lengths();
	const DISTRIBUTION distribution_to_match = m_comb_spec->get_intra_cluster_edge_lengths();

	assert(edge_length_distribution.size() == distribution_to_match.size());

	/*
	debug("The edge length distribution spec.");
	copy(distribution_to_match.begin(), distribution_to_match.end(), 
				ostream_iterator<NUM_ELEMENTS>(cout, " ")); cout << endl;

	debug("The edge length distribution found.");
	copy(edge_length_distribution.begin(), edge_length_distribution.end(), 
				ostream_iterator<NUM_ELEMENTS>(cout, " ")); cout << endl;
	*/

	edge_length_distribution_satisfied = (edge_length_distribution == distribution_to_match);

	return edge_length_distribution_satisfied;
}

// an inter-cluster edge that inputs into this cluster has changed.
// update data structures
//
// POST: date structures have been updated
//
void CLUSTER::inter_cluster_input_edge_has_changed
(
	const EDGE * edge,
	const NUM_ELEMENTS& weight_addition
)
{
	assert(edge);

	LENGTH_TYPE edge_length = edge->get_length();
	assert(edge_length >= 0 && edge_length <= m_delay);

	assert(edge->is_inter_cluster());
	m_inter_cluster_input_edge_lengths[edge_length] += weight_addition;
	m_nInter_cluster_input_edges += weight_addition;

	assert(m_inter_cluster_input_edge_lengths[edge_length] >= 0);
	assert(m_nInter_cluster_input_edges >= 0);
}
// an inter-cluster edge that outputs out of this cluster has changed.
// update data structures
//
// POST: date structures have been updated
//
void CLUSTER::inter_cluster_output_edge_has_changed
(
	const EDGE * edge,
	const NUM_ELEMENTS& weight_addition
)
{
	assert(edge);

	LENGTH_TYPE edge_length = edge->get_length();
	assert(edge_length >= 0 && edge_length <= m_delay);

	assert(edge->is_inter_cluster());
	m_inter_cluster_output_edge_lengths[edge_length] += weight_addition;

	m_nInter_cluster_output_edges += weight_addition;

	assert(m_nInter_cluster_output_edges >= 0);
	assert(m_inter_cluster_output_edge_lengths[edge_length] >= 0);
}
// an intra-cluster edge that internal to this cluster has changed.
// update data structures
//
// POST: date structures have been updated
//
void CLUSTER::intra_cluster_edge_has_changed
(
	const EDGE * edge,
	const NUM_ELEMENTS& weight_addition
)
{
	assert(edge);

	LENGTH_TYPE edge_length = edge->get_length();
	assert(edge_length >= 0 && edge_length <= m_delay);

	assert(edge->is_intra_cluster());
	m_intra_cluster_edge_lengths[edge_length] += weight_addition;

	m_nIntra_cluster_edges += weight_addition;

	assert(m_intra_cluster_edge_lengths[edge_length] >= 0);
	assert(m_nIntra_cluster_edges >= 0);
}


// get the edge length distributions
//
// POST: intra_cluster, inter_cluster_input, inter_cluster_output
//       contain the respective edge length distributions
//
void CLUSTER::get_edge_lengths
(
	DISTRIBUTION& intra_cluster, 
	DISTRIBUTION& inter_cluster_input, 
	DISTRIBUTION& inter_cluster_output
) const
{
	intra_cluster.resize(m_delay+1, 0);
	inter_cluster_input.resize(m_delay+1, 0);
	inter_cluster_output.resize(m_delay+1, 0);
	LENGTH_TYPE edge_length = 0;
	EDGES::const_iterator edge_iter;
	EDGE * edge = 0;

	for (edge_iter = m_edges.begin(); edge_iter != m_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);
		edge_length = edge->get_length();
		assert(edge_length > 0 && edge_length <= m_delay && (edge->get_weight() >= 0));

		if (edge->is_intra_cluster())
		{
			intra_cluster[edge_length] += edge->get_weight();
		}
		else if (edge->is_inter_cluster_input(m_cluster_number))
		{
			inter_cluster_input[edge_length] += edge->get_weight();
		}
		else 
		{
			assert(edge->is_inter_cluster_output(m_cluster_number));
			inter_cluster_output[edge_length] += edge->get_weight();
		}
	}

	assert(intra_cluster == m_intra_cluster_edge_lengths);
}

// RETURN: true if the inter-cluster input edge length distribution is equal to its spec
bool CLUSTER::is_inter_cluster_input_edge_length_distribution() const
{
	//debug("Checking inter cluster input edge length");
	bool edge_length_distribution_satisfied = false;

	const DISTRIBUTION distribution_to_match = m_comb_spec->get_inter_cluster_input_edge_lengths();

	assert(m_inter_cluster_input_edge_lengths.size() == distribution_to_match.size());

	
	/*
	debug("The edge length distribution spec.");
	copy(distribution_to_match.begin(), distribution_to_match.end(), 
				ostream_iterator<NUM_ELEMENTS>(cout, " ")); cout << endl;

	debug("The edge length distribution found.");
	copy(edge_length_distribution.begin(), edge_length_distribution.end(), 
				ostream_iterator<NUM_ELEMENTS>(cout, " ")); cout << endl;
				*/


	edge_length_distribution_satisfied = (m_inter_cluster_input_edge_lengths == distribution_to_match);

	return edge_length_distribution_satisfied;

}


// RETURN: true if the inter-cluster output edge length distribution is equal to its spec
bool CLUSTER::is_inter_cluster_output_edge_length_distribution() const
{
	bool edge_length_distribution_satisfied = false;

	//debug("checking inter cluster output edge length");

	const DISTRIBUTION edge_length_distribution = get_inter_cluster_output_edge_lengths();
	const DISTRIBUTION distribution_to_match = m_comb_spec->get_inter_cluster_output_edge_lengths();
	
	assert(edge_length_distribution.size() == distribution_to_match.size());

	/*
	debug("The edge length distribution spec.");
	copy(distribution_to_match.begin(), distribution_to_match.end(), 
				ostream_iterator<NUM_ELEMENTS>(cout, " ")); cout << endl;

	debug("The edge length distribution found.");
	copy(edge_length_distribution.begin(), edge_length_distribution.end(), 
				ostream_iterator<NUM_ELEMENTS>(cout, " ")); cout << endl;
				*/


	edge_length_distribution_satisfied = (edge_length_distribution == distribution_to_match);


	return edge_length_distribution_satisfied;
}

// RETURN: true if the input, output shapes are equal to their spec 
//         and all level nodes have enough length 1 edges to define the delay levels
//         of the indiv. nodes
//       
bool CLUSTER::are_level_node_properties_satisfied() const
{
	bool valid_solution = true;
	LEVEL_NODES::const_iterator level_nodes_iter = Level.begin();
	LEVEL_NODE * level_node = 0;

	//debug("checking are level nodes props satisfied");

	while (valid_solution && level_nodes_iter != Level.end())
	{
		level_node = *level_nodes_iter;
		assert(level_node);

		valid_solution = 	level_node->is_input_degree_satisfied() &&
		        			level_node->is_output_degree_satisfied() &&
			        		level_node->is_delay_satisfied();

		level_nodes_iter++;
	}

	return valid_solution;
}
// RETURN: true if all level nodes have enough length 1 edges to define the delay levels
//         of the indiv. nodes
bool CLUSTER::is_delay_satisfied() const
{
	bool valid_solution = true;
	LEVEL_NODES::const_iterator level_nodes_iter = Level.begin();
	LEVEL_NODE * level_node = 0;

	//debug("checking are level nodes props satisfied");

	while (valid_solution && level_nodes_iter != Level.end())
	{
		level_node = *level_nodes_iter;
		assert(level_node);

		valid_solution = level_node->is_delay_satisfied();

		level_nodes_iter++;
	}

	return valid_solution;
}
NUM_ELEMENTS CLUSTER::get_nIntra_cluster_edges() const
{
	//NUM_ELEMENTS nEdges = static_cast<NUM_ELEMENTS>(
	//accumulate(m_intra_cluster_edge_lengths.begin(), m_intra_cluster_edge_lengths.end(), 0));

	return m_nIntra_cluster_edges;
}
NUM_ELEMENTS CLUSTER::get_nInter_cluster_input_edges() const
{
	//NUM_ELEMENTS nEdges = static_cast<NUM_ELEMENTS>(
	//accumulate(m_inter_cluster_input_edge_lengths.begin(), m_inter_cluster_input_edge_lengths.end(), 0));

	return m_nInter_cluster_input_edges;
}
NUM_ELEMENTS CLUSTER::get_nInter_cluster_output_edges() const
{
	//NUM_ELEMENTS nEdges = static_cast<NUM_ELEMENTS>(
	//accumulate(m_inter_cluster_output_edge_lengths.begin(), m_inter_cluster_output_edge_lengths.end(), 0));

	return m_nInter_cluster_output_edges;
}

// RETURN: the number of intra-cluster edges that are internal to this cluster with 
//         edge length length 
NUM_ELEMENTS CLUSTER::get_nIntra_cluster_edges(const LENGTH_TYPE& length) const
{
	assert(length >= 0 && static_cast<unsigned>(length) < m_intra_cluster_edge_lengths.size());
	return m_intra_cluster_edge_lengths[length];
}

// RETURN: the number of inter-cluster input edges that enter into to this cluster with 
//         edge length length 
NUM_ELEMENTS CLUSTER::get_nInter_cluster_input_edges(const LENGTH_TYPE& length) const
{
	assert(length >= 0 && static_cast<unsigned>(length) < m_inter_cluster_input_edge_lengths.size());
	return m_inter_cluster_input_edge_lengths[length];
}
// RETURN: the number of inter-cluster output edges that exit out of this cluster with 
//         edge length length 
NUM_ELEMENTS CLUSTER::get_nInter_cluster_output_edges(const LENGTH_TYPE& length) const
{
	assert(length >= 0 && static_cast<unsigned>(length) < m_inter_cluster_output_edge_lengths.size());
	return m_inter_cluster_output_edge_lengths[length];
}


// RETURN: width = max of the node_shape for this cluster
NUM_ELEMENTS CLUSTER::get_width() const
{
	NUM_ELEMENTS width = 0;
    DELAY_TYPE delay_level = 0;

	for (delay_level = 0; delay_level<= m_delay; delay_level++) 
	{
		width = MAX(width, Level[delay_level]->get_nNodes());
	}

	return width;
}

// RETURNS: cost_Level_Shape + cost_Problem_Node 
//          for the cluster
COST_TYPE CLUSTER::get_level_node_cost() const
{
	COST_TYPE cost = 0.0;
	LEVEL_NODES::const_iterator level_node_iter;
	LEVEL_NODE * level_node = 0;

	for (level_node_iter = Level.begin(); level_node_iter != Level.end(); level_node_iter++)
	{
		level_node = *level_node_iter;
		assert(level_node);

		cost += level_node->get_cost();
	}
	//debug("Level node cost = " << cost);

	return cost;
}

// get the difference of the shape function from spec
// don't add in any congestion factors or node violations costs
//
// RETURN: the absolute difference between the input and output shapes and their spec
COST_TYPE CLUSTER::get_shape_cost() const
{
	COST_TYPE cost = 0.0;
	LEVEL_NODES::const_iterator level_node_iter;
	LEVEL_NODE * level_node = 0;

	for (level_node_iter = Level.begin(); level_node_iter != Level.end(); level_node_iter++)
	{
		level_node = *level_node_iter;
		assert(level_node);

		cost += level_node->get_input_cost() + level_node->get_output_cost();
	}
	//debug("Level node cost = " << cost);

	return cost;
}

// RETURN: the ABS difference between the intra-cluster, inter-cluster input, inter-cluster output 
//         edge length distributions and their specification 
//         normalized by the sum of their specs
COST_TYPE CLUSTER::get_normalized_edge_length_cost() const
{
	COST_TYPE cost = 0.0;
	COST_TYPE intra_edge_cost = 0.0;
	COST_TYPE inter_input_cost = 0.0;
	COST_TYPE inter_output_cost = 0.0;
	LENGTH_TYPE edge_length = 0;

	DISTRIBUTION spec_intra_cluster 		= m_comb_spec->get_intra_cluster_edge_lengths();
	DISTRIBUTION spec_inter_cluster_input	= m_comb_spec->get_inter_cluster_input_edge_lengths();
	DISTRIBUTION spec_inter_cluster_output	= m_comb_spec->get_inter_cluster_output_edge_lengths();
	double spec_nIntra_cluster_edges = static_cast<double>(m_comb_spec->get_nIntra_cluster_edges());
	double spec_nInter_cluster_inputs =  static_cast<double>(m_comb_spec->get_nInter_cluster_input_edges());
	double spec_nInter_cluster_outputs = static_cast<double>(m_comb_spec->get_nInter_cluster_output_edges());


	for (edge_length = 0; edge_length <= m_delay; edge_length++)
	{
		intra_edge_cost += ABS(m_intra_cluster_edge_lengths[edge_length] - spec_intra_cluster[edge_length]);
		inter_input_cost += ABS(m_inter_cluster_input_edge_lengths[edge_length] - 
								spec_inter_cluster_input[edge_length]);
		inter_output_cost += ABS(m_inter_cluster_output_edge_lengths[edge_length] -  
								spec_inter_cluster_output[edge_length]);
	}
	//debug("edge length cost = " << intra_edge_cost + inter_input_cost + inter_output_cost);
	cost = (intra_edge_cost + inter_input_cost + inter_output_cost)/
			(spec_nIntra_cluster_edges + spec_nInter_cluster_inputs + spec_nInter_cluster_outputs);

	//debug("Normalized edge length cost = " << cost);

	return cost;
}

// RETURN: the ABS difference between the intra-cluster, inter-cluster input, inter-cluster output 
//         edge length distributions and their specification 
COST_TYPE CLUSTER::get_edge_length_cost() const
{
	COST_TYPE cost = 0.0;
	LENGTH_TYPE edge_length = 0;

	DISTRIBUTION spec_intra_cluster 		= m_comb_spec->get_intra_cluster_edge_lengths();
	DISTRIBUTION spec_inter_cluster_input	= m_comb_spec->get_inter_cluster_input_edge_lengths();
	DISTRIBUTION spec_inter_cluster_output	= m_comb_spec->get_inter_cluster_output_edge_lengths();


	for (edge_length = 0; edge_length <= m_delay; edge_length++)
	{
		cost += ABS(m_intra_cluster_edge_lengths[edge_length] -  spec_intra_cluster[edge_length]);
		cost += ABS(m_inter_cluster_input_edge_lengths[edge_length] -  spec_inter_cluster_input[edge_length]);
		cost += ABS(m_inter_cluster_output_edge_lengths[edge_length] -  spec_inter_cluster_output[edge_length]);
	}

	//debug("Edge length cost = " << cost);

	return cost;
}

// return an edge in the cluster with the give source delay level
// and given edge length
//
// RETURN: an edge connecting from source_node to this level node if one could be found
//         otherwise 0
//
EDGE * CLUSTER::get_edge
(
	const DELAY_TYPE& source_delay_level, 
	const LENGTH_TYPE & edge_length
) const
{
	DELAY_TYPE sink_delay = source_delay_level + edge_length;
	assert(	source_delay_level >= 0 && source_delay_level <= m_delay &&
			sink_delay >= 0 && sink_delay <= m_delay);

	LEVEL_NODE * sink_level_node = Level[sink_delay];
	assert(sink_level_node);

	EDGE * edge = sink_level_node->find_edge_with_source(m_cluster_number, source_delay_level);
	assert(edge);
	LEVEL_NODE * source_level_node = Level[source_delay_level];
	assert(source_level_node == edge->get_source_level_node());
	assert(edge_length == edge->get_length());

	return edge;
}

// RETURN: the difference between the input shape and its spec
DISTRIBUTION CLUSTER::get_input_shape_difference_from_spec() const
{
	LEVEL_NODES::const_iterator level_node_iter;
	LEVEL_NODE * level_node = 0;
	DISTRIBUTION difference_by_delay_level;

	for (level_node_iter = Level.begin(); level_node_iter != Level.end(); level_node_iter++)
	{
		level_node = *level_node_iter;
		assert(level_node);

		difference_by_delay_level.push_back(level_node->get_input_difference());
	}

	return difference_by_delay_level;

}
// RETURN: the difference between the output shape and its spec
DISTRIBUTION CLUSTER::get_output_shape_difference_from_spec() const
{
	LEVEL_NODES::const_iterator level_node_iter;
	LEVEL_NODE * level_node = 0;
	DISTRIBUTION difference_by_delay_level;

	for (level_node_iter = Level.begin(); level_node_iter != Level.end(); level_node_iter++)
	{
		level_node = *level_node_iter;
		assert(level_node);

		difference_by_delay_level.push_back(level_node->get_output_difference());
	}

	return difference_by_delay_level;

}
// RETURN: the difference between the number of length 1 edges that input to this level
//         node and the number of indiv. nodes
DISTRIBUTION CLUSTER::get_delay_difference_by_delay_level() const
{
	LEVEL_NODES::const_iterator level_node_iter;
	LEVEL_NODE * level_node = 0;
	DISTRIBUTION difference_by_delay_level;

	for (level_node_iter = Level.begin(); level_node_iter != Level.end(); level_node_iter++)
	{
		level_node = *level_node_iter;
		assert(level_node);

		difference_by_delay_level.push_back(level_node->get_delay_difference());
	}

	return difference_by_delay_level;
}

// RETURN: the intra-cluster edge length distribution
DISTRIBUTION CLUSTER::get_intra_cluster_edge_lengths() const
{ 
	DISTRIBUTION intra_cluster; 
	DISTRIBUTION inter_cluster_input; 
	DISTRIBUTION inter_cluster_output;
	get_edge_lengths( intra_cluster, inter_cluster_input, inter_cluster_output);

	return m_intra_cluster_edge_lengths; 
}

//
// sets the maximum fanout
//
// POST: m_max_fanout has been set
//       the level node's max fanout has been set
//        
void CLUSTER::set_new_max_fanout
(
	const NUM_ELEMENTS& degree
)
{
	assert(degree >= 0);
	LEVEL_TYPE	level = 0;
	LEVEL_NODE*	level_node = 0;

	m_max_fanout = degree;

	for (level= 0; level <= m_delay; level++)
	{
		level_node = Level[level];
		assert(level_node);

		level_node->set_new_max_fanout(m_max_fanout);
	}
}

// POST: the nodes in the circuit have been output to a blif file
void CLUSTER::output_nodes_to_blif
(
	ofstream& blif_file
) const
{
	DELAY_TYPE	lev = 0;
	LEVEL_NODE * level_node = 0;

	blif_file << "#" << endl;

	output_dff_nodes_to_blif(blif_file);

	for (lev = 1; lev <= m_delay; lev++)
	{
		level_node = Level[lev];
		assert(level_node);
		output_comb_nodes_to_blif(level_node, blif_file);
	}

}

// POST: the combinational nodes in the circuit have been output to a blif file
void CLUSTER::output_comb_nodes_to_blif
(
	const LEVEL_NODE * level_node,
	ofstream& blif_file
) const
{
	assert(level_node);
	NODE * node = 0;
	NODES nodes = level_node->get_nodes();
	NODES::const_iterator node_iter;
	string input_node_names;
	NUM_ELEMENTS nInputs;

	for (node_iter = nodes.begin(); node_iter != nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		input_node_names  = node->get_input_node_names();
		nInputs = node->get_input_edges().size();
	

		blif_file << ".names " << input_node_names << " " << node->get_name() << "\n";

		while (nInputs > 0)
		{
			blif_file << "1";
			nInputs--;
		}
		blif_file << " 1" << endl;
	}

	blif_file << "#" << endl;
}

// POST: the flip-flop nodes in the circuit have been output to a blif file
void CLUSTER::output_dff_nodes_to_blif
(
	ofstream& blif_file
) const
{
	LEVEL_NODE * level_node = Level[0];
	assert(level_node);
	NODE * node = 0;
	NODES DFF_nodes = level_node->get_DFF();
	NODES::const_iterator node_iter;
	string input_node_names;

	for (node_iter = DFF_nodes.begin(); node_iter != DFF_nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);
		assert(node->is_DFF());
		input_node_names  = node->get_input_node_names();

		blif_file << ".latch " << input_node_names << " " << node->get_name() << " re clock 2\n";
	}

	blif_file << "#" << endl;
}

// PRE: m_edges contains valid super edges
// RETURN: a randomly chosen intra-cluster edge 
EDGE * CLUSTER::randomly_choose_an_intra_cluster_edge() const
{
	NUM_ELEMENTS random_edge_number = g_rand_number_gen->random_number(m_nIntra_cluster_edges-1);

	return m_edges[random_edge_number];
}

// PRE: we have assigned the fanouts to the indiv. nodes
// POST: m_max_fanout has been increased if need be
// RETURN: fanout distribution for the indiv. nodes
DISTRIBUTION CLUSTER::get_degree_distribution()
{
	NUM_ELEMENTS fanout_degree = 0;
	NODE * node = 0;
	NODES::const_iterator node_iter;
	DISTRIBUTION degree_distribution(m_max_fanout + 1, 0);

	for (node_iter = m_nodes.begin(); node_iter != m_nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		fanout_degree = node->get_nComb_Outputs();
		assert(fanout_degree >= 0);
		if (fanout_degree > m_max_fanout)
		{
			Warning("The fanout of node " << node->get_name() << " is greater than the "
					<< "fanout of the cluster.  This should not happen.");

			assert(new_gen);// find out if it does happen

			degree_distribution.resize(fanout_degree+1, 0);
			m_max_fanout = fanout_degree;
		}

		degree_distribution[fanout_degree] += 1;
	}

	return degree_distribution;
}

// RETURN: the difference between the fanout distribution assigned to the indiv.
//         nodes and the specification
NUM_ELEMENTS CLUSTER::get_degree_difference()
{
	NUM_ELEMENTS degree_difference = 0;
	DISTRIBUTION degree_distribution = get_degree_distribution();
	DISTRIBUTION spec_degree_distribution = m_comb_spec->get_fanout_distribution();
	spec_degree_distribution.resize(m_max_fanout +1, 0);

	NUM_ELEMENTS spec_max_fanout_degree = m_max_fanout,
					max_fanout_degree = m_max_fanout;

	assert(spec_degree_distribution.size() == degree_distribution.size());

	while (max_fanout_degree >= 0)
	{
		// find the highest degree assigned to in this comb graph
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

			degree_distribution[max_fanout_degree] -= 1;
			spec_degree_distribution[spec_max_fanout_degree] -= 1;
		}
	}

	return degree_difference;
}

// RETURN: the total of the sum of the abs difference between the input and output
//         shapes and their spec
NUM_ELEMENTS CLUSTER::get_input_output_difference() const
{
	NUM_ELEMENTS input_output_difference = 0,
				input_difference = 0,
				output_difference = 0,
				total_input_difference = 0,
				total_output_difference = 0;
	DELAY_TYPE	lev = 0;
	LEVEL_NODE*	level_node = 0;
	SHAPE level_in = m_comb_spec->get_level_in();
	SHAPE level_out = m_comb_spec->get_level_out();

	if (g_options->is_verbose()) debug("Looking at cluster " << m_cluster_number);

	for (lev= 0; lev <= m_delay; lev++)
	{
		level_node = Level[lev];
		assert(level_node);

		if (lev > 0)
		{
			input_difference = level_node->get_nSimple_input_edges() - level_in[lev];
			total_input_difference += ABS(input_difference);
		}
		output_difference = level_node->get_nSimple_comb_output_edges() - level_out[lev];
		total_output_difference += ABS(output_difference);

		if (g_options->is_verbose())
		{
			debug("lev: " << lev << "\tin diff: " << input_difference 
				<< "\tout diff: " << output_difference);
		}
	}

	if (g_options->is_verbose())
	{
		debug("total input difference: " << total_input_difference);
		debug("total output difference: " << total_output_difference);
	}
	input_output_difference = total_input_difference + total_output_difference;

	return input_output_difference;
}

// RETURN: the intra-cluster edge length distribution for the indiv. edges
DISTRIBUTION CLUSTER::get_intra_cluster_simple_edge_lengths() const
{
	EDGES::const_iterator edge_iter;
	EDGE * edge = 0;
	LENGTH_TYPE length = 0;

	DISTRIBUTION intra_cluster_edge_lengths(m_delay+1, 0);

	for (edge_iter = m_intra_cluster_simple_edges.begin(); 
		edge_iter != m_intra_cluster_simple_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		length = edge->get_length();
		assert(length >= 0 && length <= m_delay);

		intra_cluster_edge_lengths[length] += 1;
	}

	return intra_cluster_edge_lengths;
}
// RETURN: the inter-cluster input edge length distribution for the indiv. edges
DISTRIBUTION CLUSTER::get_inter_cluster_simple_input_edge_lengths() const
{
	EDGES::const_iterator edge_iter;
	EDGE * edge = 0;
	LENGTH_TYPE length = 0;

	DISTRIBUTION inter_cluster_input_edge_lengths(m_delay+1, 0);

	for (edge_iter = m_inter_cluster_input_simple_edges.begin(); 
		edge_iter != m_inter_cluster_input_simple_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		length = edge->get_length();
		assert(length >= 0 && length <= m_delay);

		inter_cluster_input_edge_lengths[length] += 1;
	}

	return inter_cluster_input_edge_lengths;
}
// RETURN: the inter-cluster output edge length distribution for the indiv. edges
DISTRIBUTION CLUSTER::get_inter_cluster_simple_output_edge_lengths() const
{
	EDGES::const_iterator edge_iter;
	EDGE * edge = 0;
	LENGTH_TYPE length = 0;

	DISTRIBUTION inter_cluster_output_edge_lengths(m_delay+1, 0);

	for (edge_iter = m_inter_cluster_output_simple_edges.begin(); 
		edge_iter != m_inter_cluster_output_simple_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		length = edge->get_length();
		assert(length >= 0 && length <= m_delay);

		inter_cluster_output_edge_lengths[length] += 1;
	}

	return inter_cluster_output_edge_lengths;
}
NUM_ELEMENTS CLUSTER::get_edge_length_difference() const
{
	NUM_ELEMENTS difference = 0;
	LENGTH_TYPE edge_length = 0;

	DISTRIBUTION intra_cluster_edge_lengths = get_intra_cluster_simple_edge_lengths();
	DISTRIBUTION inter_cluster_input_edge_lengths = get_inter_cluster_simple_input_edge_lengths();
	DISTRIBUTION inter_cluster_output_edge_lengths = get_inter_cluster_simple_output_edge_lengths();
	DISTRIBUTION spec_intra_cluster 		= m_comb_spec->get_intra_cluster_edge_lengths();
	DISTRIBUTION spec_inter_cluster_input	= m_comb_spec->get_inter_cluster_input_edge_lengths();
	DISTRIBUTION spec_inter_cluster_output	= m_comb_spec->get_inter_cluster_output_edge_lengths();

	for (edge_length = 1; edge_length <= m_delay; edge_length++)
	{
		difference += ABS(intra_cluster_edge_lengths[edge_length] -  spec_intra_cluster[edge_length]);
		difference += ABS(inter_cluster_input_edge_lengths[edge_length] -  spec_inter_cluster_input[edge_length]);
		difference += ABS(inter_cluster_output_edge_lengths[edge_length] -  spec_inter_cluster_output[edge_length]);
	}

	return difference;
}

// RETURN: the difference between the number of POs and the number of the spec says
NUM_ELEMENTS CLUSTER::get_po_difference() const 
{ 
	return ABS(m_nPO - m_comb_spec->get_nPO()); 
}



// POST: We are sure that our intra-cluster, 
//       inter-cluster input and inter-cluster output
//       edge length distributions are consistent with 
//       with the edge length distribution data structures 
//       in our level nodes.
//
void CLUSTER::sanity_check() const
{
	// check the distributions

	DISTRIBUTION intra_cluster_edge_lengths(m_delay + 1, 0);
	DISTRIBUTION inter_cluster_input_edge_lengths(m_delay + 1, 0);
	DISTRIBUTION inter_cluster_output_edge_lengths(m_delay + 1, 0);

	DISTRIBUTION	level_nodes_intra_cluster;
	DISTRIBUTION	level_nodes_inter_cluster_input;
	DISTRIBUTION	level_nodes_inter_cluster_output;

	LEVEL_NODE * level_node = 0;
	LEVEL_NODES::const_iterator level_node_iter;
	LENGTH_TYPE length = 0;

	for (level_node_iter = Level.begin(); level_node_iter != Level.end(); level_node_iter++)
	{
		level_node = *level_node_iter;
		assert(level_node);

		level_nodes_intra_cluster = level_node->get_intra_cluster_output_edge_lengths();
		level_nodes_inter_cluster_input = level_node->get_inter_cluster_input_edge_lengths();
		level_nodes_inter_cluster_output= level_node->get_inter_cluster_output_edge_lengths();

		assert(level_nodes_intra_cluster.size() == level_nodes_inter_cluster_input.size() &&
				level_nodes_intra_cluster.size() == level_nodes_inter_cluster_output.size());
		assert(static_cast<unsigned>(m_delay+1) == level_nodes_inter_cluster_output.size());

		for (length = 0; length <= m_delay; length++)
		{
			intra_cluster_edge_lengths[length] += level_nodes_intra_cluster[length];
			inter_cluster_input_edge_lengths[length] += level_nodes_inter_cluster_input[length];
			inter_cluster_output_edge_lengths[length] += level_nodes_inter_cluster_output[length];
		}
	}


	// make sure that they are equal
	for (length = 0; length <= m_delay; length++)
	{
		assert(intra_cluster_edge_lengths[length] == m_intra_cluster_edge_lengths[length]);
		assert(inter_cluster_input_edge_lengths[length] == m_inter_cluster_input_edge_lengths[length]);
		assert(inter_cluster_output_edge_lengths[length] == m_inter_cluster_output_edge_lengths[length]);
	}
}

// RETURN: the latched shape
SHAPE CLUSTER::get_latched_shape() const
{
	DELAY_TYPE lev = 0;
	LEVEL_NODE * level_node = 0;
	SHAPE latched_shape;
	NUM_ELEMENTS nLatches = 0;

	for (lev = 0; lev <= m_delay; lev++)
	{
		level_node = Level[lev];
		assert(level_node);

		nLatches = level_node->get_nLatched();
		latched_shape.push_back(nLatches);
	}

	return latched_shape;
}


// POST: the level nodes in this cluster have been added to the sequential level
void CLUSTER::add_level_nodes_to_sequential_level
(
	SEQUENTIAL_LEVEL * sequential_level
)
{	
	assert(sequential_level);

	sequential_level->add_level_nodes(Level, m_cluster_number);
}

// RETURN: the final normalized cost in final edge assignment
COST_TYPE CLUSTER::get_final_improvement_cost
(
	COST_TYPE & total_bad_node_cost,
	const COST_TYPE& nEdges,
	const COST_TYPE& beta,
	const COST_TYPE& max_index
) 
{
	LEVEL_NODES::const_iterator level_iter;
	LEVEL_NODE * level_node = 0;
	COST_TYPE too_many_inputs_cost = 0,
			  double_inputs_cost = 0,
			  no_inputs_cost = 0,
			  bad_node_cost = 0,
			  wirelength_cost = 0,
			  cost = 0;

	m_wirelength = 0;
	total_bad_node_cost = 0;

	for (level_iter = Level.begin(); level_iter != Level.end(); level_iter++)
	{
		level_node = *level_iter;
		assert(level_node);

		too_many_inputs_cost = level_node->get_nExcess_inputs_on_indiv_nodes();

		no_inputs_cost = level_node->get_nIndiv_nodes_too_few_inputs_cost();
		double_inputs_cost = level_node->get_nDouble_inputs_on_indiv_nodes();

		m_wirelength += level_node->get_wirelength_cost();

		bad_node_cost += too_many_inputs_cost + double_inputs_cost + no_inputs_cost;
		total_bad_node_cost += too_many_inputs_cost + double_inputs_cost + no_inputs_cost;
	}

	bad_node_cost = bad_node_cost/nEdges;
	wirelength_cost = ABS(m_wirelength-m_wirelength_wanted)/(max_index*nEdges);

	cost = beta*wirelength_cost + (1-beta)*bad_node_cost;

	return cost;
}

// RETURN: in final edge assignment the changed cost in this cluster
//         due to a change in wirelength-approx
COST_TYPE CLUSTER::get_changed_wirelength_cost
(
	const COST_TYPE& wirelength_change_of_sink_node
) const
{
	COST_TYPE changed_cost = (m_wirelength + wirelength_change_of_sink_node - m_wirelength_wanted);

	return changed_cost;
}



CIRCUIT_GRAPH::CIRCUIT_GRAPH(CIRCUIT_SPEC* circuit_spec)
{
	NUM_ELEMENTS i;
	CLUSTER_SPEC * cluster_spec = 0;
	CLUSTER * cluster = 0;
	assert(circuit_spec);
	m_circuit_spec = circuit_spec;

	m_name					= m_circuit_spec->get_name();
	m_max_combinational_delay= m_circuit_spec->get_max_comb_delay();
	m_number_seq_nodes			= m_circuit_spec->get_nDFF();
	m_number_comb_nodes			= m_circuit_spec->get_nComb();
	m_number_edges				= m_circuit_spec->get_nEdges();
	m_nPI						= m_circuit_spec->get_nPI();
	m_nPO						= m_circuit_spec->get_nPO();
	m_id						= 0;
	m_edge_id					= 0;
	
	NUM_ELEMENTS nClusters = m_circuit_spec->get_nClusters();
	for (i=0; i < nClusters; i++)
	{
		cluster_spec = m_circuit_spec->get_cluster_spec(i);
		assert(cluster_spec);
		cluster = new CLUSTER(cluster_spec, i, nClusters, this);
		m_clusters.push_back(cluster);
	}

	m_sequential_level 	= new SEQUENTIAL_LEVEL(this);
	add_level_nodes(m_sequential_level);
	m_kin				= m_circuit_spec->get_kin();

	m_sanity_check_display = false;
}

CIRCUIT_GRAPH::CIRCUIT_GRAPH(const CIRCUIT_GRAPH & another_circuit)
{
	assert(false);
	m_id						= 0;

}
CIRCUIT_GRAPH & CIRCUIT_GRAPH::operator=(const CIRCUIT_GRAPH & another_circuit)
{

	assert(false);

	return (*this);

}
CIRCUIT_GRAPH::~CIRCUIT_GRAPH()
{
	CLUSTERS::iterator cluster_iter;

	for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end(); cluster_iter++)
	{
		delete *cluster_iter;
		*cluster_iter = 0;
	}
	
	m_clusters.clear();

	delete m_sequential_level;

	// need to delete simple edges
	EDGES::iterator edge_iter;
	for (edge_iter = m_simple_edges.begin(); edge_iter != m_simple_edges.end(); edge_iter++)
	{
		delete *edge_iter;
		*edge_iter = 0;
	}


	// delete the indiv. nodes
	NODES::iterator node_iter;
	NODE * node = 0;

	for (node_iter = m_nodes.begin(); node_iter != m_nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		delete node;
	}
}


// post: m_id contains the next id
// return: a new id for a node
ID_TYPE	CIRCUIT_GRAPH::get_new_id()
{
	ID_TYPE id_to_return = m_id;

	m_id++;

	return id_to_return;
}


// POST: the level nodes for all clusters have been added to the sequential level
void CIRCUIT_GRAPH::add_level_nodes
(
	SEQUENTIAL_LEVEL * sequential_level
)
{
	CLUSTERS::iterator cluster_iter;
	CLUSTER * cluster;

	for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);

		cluster->add_level_nodes_to_sequential_level(sequential_level);
	}
}

// RETURN: the PI nodes
NODES CIRCUIT_GRAPH::get_PI() const
{
	NODES PI, PI_from_cluster;
	NODES::const_iterator node_iter;
	CLUSTERS::const_iterator cluster_iter;
	CLUSTER * cluster = 0;
	LEVEL_NODE * level_node = 0;
	NODE * node = 0;

	for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);
	
		level_node = cluster->Level[0];
		assert(level_node);
		PI_from_cluster = level_node->get_PI();

		for (node_iter = PI_from_cluster.begin(); node_iter != PI_from_cluster.end(); node_iter++)
		{
			node = *node_iter;
			assert(node);

			PI.push_back(node);
		}
	}

	return PI;
}


// POST: we have cleared our simple edges
void CLUSTER::clear_simple_edges()
{
	m_intra_cluster_simple_edges.clear();
	m_inter_cluster_input_simple_edges.clear();
	m_inter_cluster_output_simple_edges.clear();
}




// PRE: the edge is a simple edge
// POST: we have removed the edge from the list of simple edges
void CLUSTER::remove_edge
(
	EDGE * edge
)
{
	assert(edge);
	EDGES::iterator edge_iter;

	if (edge->is_intra_cluster())
	{
		edge_iter = find(m_intra_cluster_simple_edges.begin(), m_intra_cluster_simple_edges.end(), edge);
		assert(edge_iter != m_intra_cluster_simple_edges.end());

		m_intra_cluster_simple_edges.erase(edge_iter);
	}
	else if (edge->is_inter_cluster_input(m_cluster_number))
	{
		edge_iter = find(m_inter_cluster_input_simple_edges.begin(), 
						 m_inter_cluster_input_simple_edges.end(), edge);
		
		assert(edge_iter != m_inter_cluster_input_simple_edges.end());
		m_inter_cluster_input_simple_edges.erase(edge_iter);
	}
	else
	{
		assert(edge->is_inter_cluster_output(m_cluster_number));

		edge_iter = find(m_inter_cluster_output_simple_edges.begin(), 
						 m_inter_cluster_output_simple_edges.end(), edge);
		
		assert(edge_iter != m_inter_cluster_output_simple_edges.end());
		m_inter_cluster_output_simple_edges.erase(edge_iter);
	}
}
// RETURN: the PO nodes
NODES CIRCUIT_GRAPH::get_PO() const
{
	NODES::const_iterator node_iter;
	NODE * node = 0;
	NODES PO;

	for (node_iter = m_nodes.begin(); node_iter != m_nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		if (node->is_PO())
		{
			PO.push_back(node);
		}
	}

	return PO;
}
// Create the delay structure for the circuit
void CIRCUIT_GRAPH::create_delay_structure()
{
	assert(m_sequential_level);
	m_sequential_level->create_delay_structure();
}

// Partition the fanout distribution in each cluster among the level nodes
void CIRCUIT_GRAPH::partition_degrees()
{
	CLUSTERS::iterator cluster_iter;
	CLUSTER * cluster = 0;
	NUM_ELEMENTS quality = 0;

	DEGREE_PARTITIONER degree_partitioner;

	debug("Step 2: Partitioning the Fanout Distributions");
	debug("=============================================\n");

	for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);

		dbg(sep);
		debug("Status: Partitioning degrees in cluster " << cluster->get_cluster_number());
		quality += degree_partitioner.degree_partition(cluster);
	}
	dbg(sep);
	debug("Status: Finished Partitioning the Fanout Distributions.\n");
	debug("Total difference between fanout distributions and specs: " << quality);
	debug(Sep);
}
// Split the level nodes in each cluster into indiv. nodes and prepare the graph for
// final edge assignment
void CIRCUIT_GRAPH::split_nodes()
{
	CLUSTERS::iterator cluster_iter;
	CLUSTER * cluster = 0;
	NUM_ELEMENTS quality = 0;

	NODE_SPLITTER node_splitter;


	debug("\nStep 3: Splitting the Level Nodes into Individual Nodes");
	debug("=======================================================\n");
	for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);

		dbg(sep);
		debug("Status: Splitting nodes in cluster " << cluster->get_cluster_number() << endl);
		quality += node_splitter.split_nodes(cluster);
	}
	dbg(sep);
	debug("\nTotal difference between fanout distributions and specs: " << quality);
	debugSep;
}

// RETURN: width = across all clusters the maximum of the node_shape
NUM_ELEMENTS CIRCUIT_GRAPH::get_width() const
{
	CLUSTERS::const_iterator cluster_iter;
	CLUSTER * cluster = 0;
	NUM_ELEMENTS width = 0;

	for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);

		width = MAX(width, cluster->get_width());
	}

	return width;
}
// RETURN: max_width * 16
// 
NUM_ELEMENTS CIRCUIT_GRAPH::get_max_index() const
{
	// take a look at node splitter to see where we get the 16 from
	NUM_ELEMENTS width = get_width();
	NUM_ELEMENTS max_index = 16*width;

	return max_index;
}

// POST: final edge assignment has been done
void CIRCUIT_GRAPH::assign_simple_edges()
{
	CLUSTERS::iterator cluster_iter;
	CLUSTER * cluster = 0;

	SIMPLE_EDGE_ASSIGNER simple_edge_assigner(this);


	debug("\n\nStep 4: Final Edge Assignment");
	debug("=============================\n");

	for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);

		dbg(sep);
		debug("Status: Performing initial final edges assignment in cluster " << 
				cluster->get_cluster_number() << endl);
		simple_edge_assigner.assign_simple_edges(cluster);
	}

	debug(sep << "Status: Finished Performing Initial Final Edge Assignment");

	debug(Sep);
}
// POST: a blif file has been outputed to a file
void CIRCUIT_GRAPH::output_blif()
{
	CLUSTERS::iterator cluster_iter;
	CLUSTER * cluster = 0;
	ofstream blif_file;
	string file_name = m_name + "_clone.blif";

	blif_file.open(file_name.c_str(), ios::out);

	if (! blif_file.is_open())
	{
		Warning("Could not open output file " << file_name << 
				". Therefore, could not output a blif file.");
		return;
	}
	debug("\nStatus: Opening file " << file_name << " for output.");

	blif_file << "# Gen clone\n";
	blif_file << ".model " << m_name << endl;

	output_primary_inputs_to_blif(blif_file);
	output_primary_outputs_to_blif(blif_file);

	for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);
		
		cluster->output_nodes_to_blif(blif_file);
	}


	blif_file << ".end " << endl;

	blif_file.close();

	debug("Status: Done outputing circuit to blif file.\n");
}

// POST: the primary inputs have been outputed to the blif file
void CIRCUIT_GRAPH::output_primary_inputs_to_blif
(
	ofstream& blif_file
) const
{
	NODES PI_nodes;
	NODES::const_iterator node_iter;
	NODE * node = 0;
	NUM_ELEMENTS input_counter = 0;

	blif_file << ".inputs ";

	PI_nodes = get_PI();

	for (node_iter = PI_nodes.begin(); node_iter != PI_nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		blif_file << node->get_name() << " ";

		// vpack doesn't like long input lines
		if (input_counter == 7)
		{
			 blif_file << "\\" << endl;
			 input_counter = 0;
		}
		else
		{
			input_counter++;
		}
	}

	if (m_number_seq_nodes > 0)
	{
		blif_file << "clock ";
	}

	blif_file << endl;
}

// POST: the primary outputs have been outputed to the file
void CIRCUIT_GRAPH::output_primary_outputs_to_blif
(
	ofstream& blif_file
) const
{
	CLUSTERS::const_iterator cluster_iter;
	CLUSTER * cluster = 0;
	LEVEL_NODE * level_node = 0;
	DELAY_TYPE	lev = 0;
	NODES PO_nodes;
	NODES::const_iterator node_iter;
	NODE * node = 0;
	NUM_ELEMENTS output_counter = 0;

	blif_file << ".outputs ";

	for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);

		for (lev = 0; lev <= cluster->get_delay(); lev++)
		{
	
			level_node = cluster->Level[lev];
			assert(level_node);
			PO_nodes = level_node->get_PO();

			for (node_iter = PO_nodes.begin(); node_iter != PO_nodes.end(); node_iter++)
			{
				node = *node_iter;
				assert(node);

				blif_file << node->get_name() << " ";
				// vpack doesn't like long output lines
				if (output_counter == 7)
				{
					 blif_file << "\\" << endl;
					 output_counter = 0;
				}
				else
				{
					output_counter++;
				}
			}
		}
	}

	blif_file << endl;
}


// PRE: the source and sink level nodes are valid
// POST: a super edge has been created between the two level nodes
//       the super edge has been added to the level nodes and clusters
//       m_id has been incremented
// RETURN: the super edge
EDGE * CIRCUIT_GRAPH::create_edge
(
	LEVEL_NODE * source_level_node,
	LEVEL_NODE * sink_level_node
)
{
	assert(source_level_node && sink_level_node);

	CLUSTER_NUMBER_TYPE source_cluster_number = source_level_node->get_cluster_number(),
						 sink_cluster_number = sink_level_node->get_cluster_number();
	
	EDGE * edge  = new EDGE(source_level_node, sink_level_node, m_edge_id);
	m_edge_id++;

	if (sink_level_node->get_delay_level() > 0)
	{
		// we have a normal input node
		source_level_node->add_output_edge(edge);
		sink_level_node->add_input_edge(edge);
	}
	else
	{
		assert(sink_level_node->get_delay_level() == 0);
		// the input node is a dff
		sink_level_node->add_dff_input_edge(edge);
		source_level_node->add_output_to_dff_edge(edge);
	}


	if (edge->is_intra_cluster())
	{
		assert(source_cluster_number == sink_cluster_number);
		m_clusters[sink_cluster_number]->add_edge(edge);
	}
	else
	{
		assert(source_cluster_number != sink_cluster_number);
		m_clusters[sink_cluster_number]->add_edge(edge);
		m_clusters[source_cluster_number]->add_edge(edge);
	}

	return edge;
}
// PRE: the source and sink nodes are valid
// POST: a simple/indiv. edge has been created between the two nodes
//       the edge has been added to the level nodes and clusters and the super
//       edge it is a part of 
//       the edge has been added to m_simple_edges - the list of all simple edges
// RETURN: the edge
EDGE * CIRCUIT_GRAPH::create_edge
(
	NODE * source_node, 
	NODE * sink_node
)
{
	assert(source_node && sink_node);
	CLUSTER_NUMBER_TYPE source_cluster_number = source_node->get_cluster_number(),
						 sink_cluster_number = sink_node->get_cluster_number();

	assert(source_cluster_number >= 0 && sink_cluster_number >= 0);

	EDGE * edge  = new EDGE(source_node, sink_node);

	source_node->add_output_edge(edge);
	sink_node->add_input_edge(edge);

	if (edge->is_intra_cluster())
	{
		assert(source_cluster_number == sink_cluster_number);
		m_clusters[sink_cluster_number]->add_intra_cluster_simple_edge(edge);
	}
	else
	{
		assert(source_cluster_number != sink_cluster_number);
		m_clusters[sink_cluster_number]->add_inter_cluster_input_simple_edge(edge);
		m_clusters[source_cluster_number]->add_inter_cluster_output_simple_edge(edge);
	}


	m_simple_edges.push_back(edge);

	DELAY_TYPE source_delay_level = source_node->get_delay_level();
	DELAY_TYPE sink_delay = sink_node->get_delay_level();
	LEVEL_NODE * source_level_node 
		= m_clusters[source_cluster_number]->Level[source_delay_level];
	LEVEL_NODE * sink_level_node 
		= m_clusters[sink_cluster_number]->Level[sink_delay];
	EDGE * super_edge = 
		sink_level_node->find_edge_with_source(source_level_node);

	super_edge->add_simple_edge(edge);


	return edge;
}

// POST: the edges that lost weight have had their information outputed
void CIRCUIT_GRAPH::print_lost_edge_stats()
{
	EDGES lost_edges;
	EDGES::iterator edge_iter;
	EDGE * edge = 0;

	if (m_lost_weight_edges.size() > 0)
	{
		debug("The following edges lost weight");
	}

	for (edge_iter = m_lost_weight_edges.begin(); edge_iter != m_lost_weight_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		debug(edge->get_info());
	}
}

// POST: the final quality of the synthetic circuit has been outputed
void CIRCUIT_GRAPH::output_final_circuit_quality()
{
	NUM_ELEMENTS degree_cost = 0,
				input_output_shape_cost = 0,
				edge_length_cost = 0,
				po_cost = 0,
				cost_Comb = 0,
				cost_Latched = 0,
				total_cost = 0;
				// assert(new_gen);
				//po_dff_fanout_cost don't know what to do with this deviation from spec
				//latch shape cost
				//po shape cost

	double nEdges = m_circuit_spec->get_nEdges();
	double total_normalized_cost = 0.0;

	CLUSTERS::iterator cluster_iter;
	CLUSTER * cluster = 0;

	for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);

		degree_cost += cluster->get_degree_difference();
		input_output_shape_cost += cluster->get_input_output_difference();
		edge_length_cost += cluster->get_edge_length_difference();
		po_cost += cluster->get_po_difference();
	}
	get_inter_cluster_connectivity_costs(cost_Comb, cost_Latched);

	total_cost = degree_cost + input_output_shape_cost + edge_length_cost + cost_Comb + cost_Latched +
		 		po_cost;


	debug("\n\nFinal Solution Quality");
	debug("----------------------");

	debug("fanout_distribution_cost: " << degree_cost);
	debug("input_output_shape_cost: " << input_output_shape_cost);
	debug("edge_length_distribution_cost: " << edge_length_cost);
	debug("po_cost: " << po_cost);
	debug("Comb_cost: " << cost_Comb);
	debug("Latched_cost: " << cost_Latched);

	debug("total_cost: " << total_cost);

	assert(new_gen);	// need to see what the difference in the edges_with_clock_and_output_to_dff_edges	
	debug("\nnEdges: " << nEdges);
	total_normalized_cost = static_cast<double>(total_cost)/nEdges;
	debug("total_normalized_cost: " << total_normalized_cost);
	debug(sep);

	/*
	ofstream quality_file;
	string file_name = m_name + "_clone.blif";

	blif_file.open(file_name.c_str(), ios::out);

	blif_file << "# Gen clone\n";
	*/
}


// POST: cost_Comb contains the absolute sum of all entries in (Comb - Comb spec)
// POST: cost_Latched contains the absolute sum of all entries in (Latched - Latched spec)
void CIRCUIT_GRAPH::get_inter_cluster_connectivity_costs
(
	NUM_ELEMENTS& cost_Comb,
	NUM_ELEMENTS& cost_Latched
) const
{
	INDEX_SIZE nClusters = static_cast<INDEX_SIZE>(get_nClusters());
	CLUSTER_NUMBER_TYPE source_cluster = 0,
						sink_cluster = 0;
	EDGES::const_iterator edge_iter;
	EDGE * edge = 0;
	MATRIX Comb(nClusters, nClusters),
			Latched(nClusters, nClusters),
			violation_matrix;
	MATRIX	* Comb_spec = m_circuit_spec->get_Comb(),
			* Latched_spec = m_circuit_spec->get_Latched();

	for (edge_iter = m_simple_edges.begin(); edge_iter != m_simple_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		source_cluster = edge->get_source_cluster_number();
		sink_cluster = edge->get_sink_cluster_number();

		if (edge->is_sink_a_dff())
		{
			Latched(source_cluster, sink_cluster) += 1;
		}
		else if (edge->is_inter_cluster())
		{
			assert(! edge->is_sink_a_dff());

			Comb(source_cluster, sink_cluster) += 1;
		}
		else
		{
			// nothing
		}
	}

	//debug("\nConnectivity matrix:\n" << Comb);

	violation_matrix = Comb - *Comb_spec;

	cost_Comb = violation_matrix.get_absolute_sum();

	violation_matrix = Latched - *Latched_spec;

	cost_Latched = violation_matrix.get_absolute_sum();

	//debug("cost_Comb: " << comb_cost);
	//debug("cost_Latched: " << dff_cost);
}

// POST: a sanity check has been performed and displayed
void CIRCUIT_GRAPH::final_sanity_check()
{
	debug("\nFinal Sanity Checks and Reporting on Solution Quality")
	debug("=====================================================\n");
	debug("Final Sanity Check");
	debug("------------------");
	m_sanity_check_display = true;
	sanity_check();
}
// POST: We have outputed a warning if if:
// 	     - a node has inputs above kin
// 	     - a node has more than one connection to the same source node
// 	     - a node has less than 1 output and is not a PO
// 	     - a node is a PI and PO
// 	     - a node has less than 1 input and is not a PI
// 	     - a node has less than 2 inputs and is not a PI or DFF
// 	     - a node is a flip-flop with a loop back to itself
//
// 	    We have output detailed warnings if m_sanity_check_display is true.
void CIRCUIT_GRAPH::sanity_check()
{
	NODE * node = 0;
	NODES::const_iterator node_iter;
	NUM_ELEMENTS nInputs = 0,
				 nOutputs = 0;
	NUM_ELEMENTS nExcess_input_nodes = 0,
				 nDouble_input_edge_nodes = 0,
				 nNo_output_nodes = 0,
				 nNo_input_nodes = 0,
				 nBuffer_nodes = 0,
				 nDFF_loops = 0;
	EDGES::const_iterator edge_iter;
	EDGE * edge = 0;
	
	for (edge_iter = m_simple_edges.begin(); edge_iter != m_simple_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		edge->check_sanity();
	}


	for (node_iter = m_nodes.begin(); node_iter != m_nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		nInputs = node->get_nInputs();
		nOutputs = node->get_nOutputs();

		if (nInputs > m_kin)
		{
			debugif(m_sanity_check_display,"Sanity Check: node " << node->get_name() 
					<< " has " << nInputs - m_kin << " too many inputs");
			nExcess_input_nodes++;
		}
		if (node->get_nDouble_inputs() > 0)
		{
			debugif(m_sanity_check_display,"Sanity Check: node " << node->get_name() 
					<< " has " << node->get_nDouble_inputs() << 
					" more than one connections to the same source node.");
			nDouble_input_edge_nodes++;
		}

		if (nOutputs < 1 && ! node->is_PO())
		{
			debugif(m_sanity_check_display,"Sanity Check: node " << node->get_name() << " has " << nOutputs << " outputs and is not a primary outputs.");

			nNo_output_nodes++;
		}

		if (node->is_PI() && node->is_PO())
		{
			Fail("Sanity Check: node " << node->get_name() << 
					" is a primary output and primary output");
		}
		if (node->get_nInputs() < 1 && ! node->is_PI())
		{
			debugif(m_sanity_check_display,"Sanity Check: node " << node->get_name() << 
					" has no inputs and is not a PI");
			nNo_input_nodes++;
		}
		if (node->get_nInputs() < 2 && ! node->is_PI() && ! node->is_DFF())
		{
			if (g_options->is_verbose())
			{
				debugif(m_sanity_check_display,"Sanity Check: node " << node->get_name() << 
						" is a buffer node.");
			}
			nBuffer_nodes++;
		}
		if (node->is_DFF())
		{
			if (node->is_dff_have_loop_back_to_itself())
			{
				Warning("Sanity Check: DFF " << node->get_name() << " has a loop to itself ");
				nDFF_loops++;
			}
		}
	}


	if ((nExcess_input_nodes > 0) || (nDouble_input_edge_nodes > 0) ||
	    (nNo_output_nodes > 0) || (nNo_input_nodes > 0) ||
		(nBuffer_nodes > 0) || (nDFF_loops > 0))
	{
		dbgif(m_sanity_check_display,"\n");
		Warning("Sanity Checking Circuit\n" << "--------------------------------");
	}
	if (nExcess_input_nodes > 0)
	{
		Warning("There are " << nExcess_input_nodes << " nodes with too many inputs");
	}
	if (nDouble_input_edge_nodes > 0)
	{
		Warning("There are " << nDouble_input_edge_nodes << " nodes with more than one connection to the same source node.");
	}
	if (nNo_output_nodes > 0)
	{
		Warning("There are " << nNo_output_nodes << " nodes with no output edges that are not PO");
	}
	if (nNo_input_nodes > 0)
	{
		Warning("There are " << nNo_input_nodes << " nodes that have no inputs and are not PI");
	}
	if (nBuffer_nodes > 0)
	{
		Warning("There are " << nBuffer_nodes << " buffer nodes.");
	}
	if (nDFF_loops > 0)
	{
		Warning("There are " << nDFF_loops << " dff loops.");
	}
}

// RETURN: the number of flop-flops that loop back to themselves
NUM_ELEMENTS CIRCUIT_GRAPH::get_nDFF_loops() const
{
	CLUSTERS::const_iterator cluster_iter;
	LEVEL_NODE * level_node = 0;
	NUM_ELEMENTS nDFF_loops = 0;
	CLUSTER * cluster = 0;

	for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);

		level_node = cluster->Level[0];

		nDFF_loops += level_node->get_nDFF_loops();
	}

	return nDFF_loops;
}

// PRE: the edge is a simple edge
// POST: we have removed the edge from the list of simple edges
//       we have removed the edge from its clusters
void CIRCUIT_GRAPH::remove_edge
(
	EDGE * edge
)
{
	assert(edge);
	EDGES::iterator edge_iter;

	edge_iter = find(m_simple_edges.begin(), m_simple_edges.end(), edge);
	assert(edge_iter != m_simple_edges.end());
	m_simple_edges.erase(edge_iter);


	CLUSTER_NUMBER_TYPE source_cluster = edge->get_source_cluster_number(),
	 					sink_cluster = edge->get_sink_cluster_number();
	
	m_clusters[source_cluster]->remove_edge(edge);

	if (source_cluster != sink_cluster)
	{
		m_clusters[sink_cluster]->remove_edge(edge);
	}
}

// POST: any nodes that were not connected to anything have been deleted
void CIRCUIT_GRAPH::delete_unconnected_nodes()
{
	NODES::iterator node_iter = m_nodes.begin(); 
	NODE * node = 0;
	LEVEL_NODE * level_node = 0;

	while (node_iter != m_nodes.end())
	{
		node = *node_iter;
		if (node->get_nOutputs() == 0 && node->get_nInputs() == 0)
		{

			debug("Removing node " << node->get_name() << 
					" from the circuit because it has no inputs and no outputs");
					
			// node not connected to anything. remove 
			node_iter = m_nodes.erase(node_iter);

			// remove the node from the level node's list of nodes
			level_node = node->get_my_level_node();
			level_node->remove_node(node);

			if (node->is_combinational())
			{
				m_number_comb_nodes--;
			}
			else
			{
				m_number_seq_nodes--;
			}

			delete node;
		}
		else
		{
			node_iter++;
		}

	}
}

void CIRCUIT_GRAPH::delete_node
(
	NODE * node
)
{
	NODES::iterator node_iter;
	LEVEL_NODE * level_node = 0;

	node_iter = find(m_nodes.begin(), m_nodes.end(), node);
	assert(node_iter != m_nodes.end());

	// remove the node from the list of nodes 
	m_nodes.erase(node_iter);

	// remove the node from level nodes list of nodes
	level_node = node->get_my_level_node();
	level_node->remove_node(node);

	if (node->is_combinational())
	{
		m_number_comb_nodes--;
	}
	else
	{
		m_number_seq_nodes--;
	}

	delete node;
}

// POST: a vhdl file has been outputed to a file
void CIRCUIT_GRAPH::output_vhdl()
{

	CLUSTERS::iterator cluster_iter;
	CLUSTER * cluster = 0;
	ofstream vhdl_file;
	string file_name = m_name + "_clone.vhd";
	NODES::const_iterator node_iter;
	NODE * node = 0;

	vhdl_file.open(file_name.c_str(), ios::out);

	if (! vhdl_file.is_open())
	{
		Warning("Could not open output file " << file_name << 
				". Therefore, could not output a vhdl file.");
		return;
	}
	debug("\nStatus: Opening file " << file_name << " for output.");

	vhdl_file << "LIBRARY ieee;\n";
	vhdl_file << "USE ieee.std_logic_1164.all;\n\n";

	// define a flip-flop
	vhdl_file << "ENTITY flipflop IS\n";
	vhdl_file << "\tPORT (D,Clock: IN STD_LOGIC;\n";
	vhdl_file << "\t      Q		 : OUT STD_LOGIC);\n\n";
	vhdl_file << "END flipflop;\n";

	vhdl_file << "ARCHITECTURE Behavior OF flipflop IS\n";
	vhdl_file << "BEGIN\n";
	vhdl_file << "\tPROCESS(Clock)\n";
	vhdl_file << "\tBEGIN\n";
	vhdl_file << "\t\tIF Clock'EVENT AND Clock = '1' THEN\n";
	vhdl_file << "\t\t\t Q <= D;\n";
	vhdl_file << "\t\tEND IF;\n";
	vhdl_file << "\tEND PROCESS;\n";
	vhdl_file << "END Behavior;\n" << endl;


	output_top_most_entity(vhdl_file);

	vhdl_file << "architecture behavior of " << m_name << "_clone is \n";

	vhdl_file << "COMPONENT flipflop\n";
	vhdl_file << "\t" << "PORT (D,Clock: IN STD_LOGIC;\n";
	vhdl_file << "\t" << "      Q      : OUT STD_LOGIC);\n";
	vhdl_file << "END COMPONENT;\n";

	for (node_iter = m_nodes.begin(); node_iter != m_nodes.end(); node_iter++)
	{
		node = *node_iter;

		if (! node->is_PI())
		{
			if (node->is_PO())
			{
				vhdl_file << "\tSIGNAL tmp_" << node->get_name() << " : STD_LOGIC;\n";
			}
			else
			{
				vhdl_file << "\tSIGNAL " << node->get_name() << " : STD_LOGIC;\n";
			}
		}
	}



	vhdl_file << "BEGIN\n";

	for (cluster_iter = m_clusters.begin(); cluster_iter != m_clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);
		
		cluster->output_nodes_to_vhdl(vhdl_file);
	}

	NODES PO = get_PO();

	// output the output_signal <= tmp_signal_name
	for (node_iter = PO.begin(); node_iter != PO.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		vhdl_file << "\t" << node->get_name() << " <= tmp_" << node->get_name() << ";" << endl;
	}

	vhdl_file << "END behavior;\n";

	vhdl_file.close();

	debug("Status: Done outputing circuit to vhdl file.\n");
}
// POST: the primary inputs have been outputed to the vhdl file
void CIRCUIT_GRAPH::output_top_most_entity
(
	ofstream& vhdl_file
) const
{
	NODES::const_iterator node_iter;
	NODE * node = 0;
	NODES PI = get_PI(),
	      PO = get_PO();
	NUM_ELEMENTS po_number = 1,
	             nPO = static_cast<NUM_ELEMENTS>(PO.size());

	vhdl_file << "LIBRARY ieee;\n";
	vhdl_file << "USE ieee.std_logic_1164.all;\n\n";
	vhdl_file << "USE WORK.all;\n\n";

	// define a flip-flop
	vhdl_file << "ENTITY " << m_name <<  "_clone IS\n"; 
	vhdl_file << "\tPORT (" << endl;
	
	for (node_iter = PI.begin(); node_iter != PI.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		vhdl_file << "\t\t" << node->get_name() << ": IN \tSTD_LOGIC;" << endl;

	}
	if (is_sequential())
	{
		vhdl_file << "\t\t clock: IN \tSTD_LOGIC;";
	}
	vhdl_file << endl;

	assert(nPO >= 1);

	for (node_iter = PO.begin(); node_iter != PO.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		if (po_number != nPO)
		{
			vhdl_file << "\t\t" << node->get_name() << ": OUT \tSTD_LOGIC;" << endl;
		}
		else
		{
			vhdl_file << "\t\t" << node->get_name() << ": OUT \tSTD_LOGIC);" << endl;
		}

		po_number++;
	}
	vhdl_file << "END " << m_name << "_clone;" << endl << endl;
}
// POST: the nodes in the circuit have been output to a vhdl file
void CLUSTER::output_nodes_to_vhdl
(
	ofstream& vhdl_file
) const
{
	DELAY_TYPE	lev = 0;
	LEVEL_NODE * level_node = 0;

	output_dff_nodes_to_vhdl(vhdl_file);

	for (lev = 1; lev <= m_delay; lev++)
	{
		level_node = Level[lev];
		assert(level_node);
		output_comb_nodes_to_vhdl(level_node, vhdl_file);
	}

}
// POST: the flip-flop nodes in the circuit have been output to a vhdl file
void CLUSTER::output_dff_nodes_to_vhdl
(
	ofstream& vhdl_file
) const
{
	LEVEL_NODE * level_node = Level[0];
	assert(level_node);

	NODE * node = 0;
	NODES DFF_nodes = level_node->get_DFF();
	NODES::const_iterator node_iter;
	EDGE * edge = 0;

	for (node_iter = DFF_nodes.begin(); node_iter != DFF_nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);
		assert(node->is_DFF());

    	vhdl_file << "flip_flop_" << node->get_name() << " : flipflop port MAP (\n";

		edge  = node->get_dff_input_edge();

		vhdl_file <<  "\t" << edge->get_vhdl_signal_name() << ", clock, "
		 		  <<  node->get_name() << ");\n\n";
	}

	vhdl_file << endl;
}

// POST: the combinational nodes in the circuit have been output to a vhdl file
void CLUSTER::output_comb_nodes_to_vhdl
(
	const LEVEL_NODE * level_node,
	ofstream& vhdl_file
) const
{
	assert(level_node);
	NODE * node = 0;
	NODES nodes = level_node->get_nodes();
	NODES::const_iterator node_iter;
	NUM_ELEMENTS input_number = 0,
				 nInputs = 0;
	EDGES input_edges;
	EDGES::const_iterator edge_iter;
	EDGE * edge = 0;

	for (node_iter = nodes.begin(); node_iter != nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		if (node->is_PO())
		{
			vhdl_file << "\ttmp_" << node->get_name() << " <= NOT (";
		}
		else
		{
			vhdl_file << "\t" << node->get_name() << " <= NOT (";

		}

		input_edges  = node->get_input_edges();
		nInputs = static_cast<NUM_ELEMENTS>(input_edges.size());

		input_number = 0;
		for (edge_iter = input_edges.begin(); edge_iter != input_edges.end(); edge_iter++)
		{
			edge = *edge_iter;
			assert(edge);

			// if not the last input name
			if (input_number < nInputs - 1)
			{
				vhdl_file <<   edge->get_vhdl_signal_name() << " AND ";
			}
			else
			{
				vhdl_file <<   edge->get_vhdl_signal_name() << ");\n";
			}

			input_number++;
		}
	}

	vhdl_file << endl;
}
