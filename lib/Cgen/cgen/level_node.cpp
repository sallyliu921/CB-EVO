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




#include "level_node.h"
#include "gen.h"
#include <algorithm>
#include <numeric>
#include "util.h"

LEVEL_NODE::LEVEL_NODE
(
	const DELAY_TYPE & delay_level, 
	CLUSTER_SPEC * comb_spec,
	const NUM_ELEMENTS & nClusters,
	CLUSTER * my_cluster
)
{
	assert(comb_spec && my_cluster && g_options);

	NODE * node = 0;

	NUM_ELEMENTS max_fanout = comb_spec->get_max_fanout();
	DELAY_TYPE delay = comb_spec->get_delay();
	m_intra_cluster_output_edge_lengths.resize(delay + 1, 0); 
	m_intra_cluster_input_edge_lengths.resize(delay+1, 0);
	m_inter_cluster_input_edge_lengths.resize(delay + 1, 0);
	m_inter_cluster_output_edge_lengths.resize(delay + 1, 0);
	m_out_degrees.resize(max_fanout+1, 0);

	m_delay_level	= delay_level;
	m_cluster_number= comb_spec->get_cluster_number();

	m_nNodes = comb_spec->get_shape()[m_delay_level];

	if (m_delay_level == 0)
	{
		m_nDFF = comb_spec->get_nDFF();
	}
	else
	{
		m_nDFF = 0;
	}

	m_kin 			= comb_spec->get_kin();
	m_nPO			= comb_spec->get_PO_shape()[m_delay_level];
	m_nlatched		= comb_spec->get_latched_shape()[m_delay_level];
	m_min_in  		= 0;
	m_max_in		= 0;
	m_min_out 		= 0;
	m_max_out		= 0;
	m_assigned_in 	= comb_spec->get_level_in()[m_delay_level];
	m_assigned_out 	= comb_spec->get_level_out()[m_delay_level];
	m_target_fanout = 0;
	m_available_in 	= 0;
	m_available_out	= 0;

	m_nIntra_cluster_inputs 	= 0;
	m_nIntra_cluster_outputs	= 0;
	m_nInter_cluster_inputs 	= 0;
	m_nInter_cluster_outputs	= 0;

	m_cluster = my_cluster;
	CIRCUIT_GRAPH * circuit = m_cluster->get_my_circuit();
	assert(circuit);
	ID_TYPE node_id = 0;

	for (int i=0; i < m_nNodes; i++)
	{
		node_id = circuit->get_new_id();
		node = new NODE(m_cluster_number, m_delay_level, node_id, m_kin, this);
		circuit->add_node(node);
		m_cluster->add_node(node);
		m_nodes.push_back(node);
	}


	// create a name for the node based on cluster "C" and delay level "D"
	m_name = "C" + util_long_to_string(m_cluster_number) + "D" + util_long_to_string(m_delay_level);

	m_input_edges.resize(nClusters*m_delay_level,0);

	m_max_fanout = comb_spec->get_max_fanout();

	m_too_many_inputs_factor = g_options->get_delay_structure_init_too_many_inputs_factor();
	m_too_few_outputs_factor = g_options->get_delay_structure_init_too_few_outputs_factor();
	m_too_few_inputs_factor  = g_options->get_delay_structure_init_too_few_inputs_factor();

	m_mult_too_many_inputs_factor = g_options->get_delay_structure_multiplier_too_many_inputs_factor();
	m_mult_too_few_outputs_factor = g_options->get_delay_structure_multiplier_too_few_outputs_factor();
	m_mult_too_few_inputs_factor  = g_options->get_delay_structure_multiplier_too_few_inputs_factor();
}


LEVEL_NODE::~LEVEL_NODE()
{
}
LEVEL_NODE::LEVEL_NODE(const LEVEL_NODE  & another_node)
{
	m_nodes			= another_node.m_nodes;

	m_intra_cluster_output_edge_lengths	= another_node.m_intra_cluster_output_edge_lengths; 
	m_inter_cluster_input_edge_lengths 	= another_node.m_inter_cluster_input_edge_lengths;
	m_inter_cluster_output_edge_lengths	= another_node.m_inter_cluster_output_edge_lengths;

	m_out_degrees	= another_node.m_out_degrees;
	
	m_nNodes		= another_node.m_nNodes;
	m_nPO			= another_node.m_nPO;
	m_nlatched		= another_node.m_nlatched;
	m_min_in  		= another_node.m_min_in;
	m_max_in		= another_node.m_max_in;
	m_min_out 		= another_node.m_min_out;
	m_max_out		= another_node.m_max_out;
	m_assigned_in 	= another_node.m_assigned_in;
	m_assigned_out 	= another_node.m_assigned_out ;
	m_target_fanout = another_node.m_target_fanout;
	m_available_in	= another_node.m_available_in;
	m_available_out	= another_node.m_available_out;

	m_nInter_cluster_inputs 	= another_node.m_nInter_cluster_inputs;
	m_nInter_cluster_outputs	= another_node.m_nInter_cluster_outputs;

	m_kin						= another_node.m_kin;
}
LEVEL_NODE & LEVEL_NODE::operator=(const LEVEL_NODE  & another_node)
{
	m_nodes			= another_node.m_nodes;

	m_intra_cluster_output_edge_lengths = another_node.m_intra_cluster_output_edge_lengths; 
	m_inter_cluster_input_edge_lengths 	= another_node.m_inter_cluster_input_edge_lengths;
	m_inter_cluster_output_edge_lengths	= another_node.m_inter_cluster_output_edge_lengths;

	m_out_degrees	= another_node.m_out_degrees;
	
	m_nNodes		= another_node.m_nNodes;
	m_nPO			= another_node.m_nPO;
	m_nlatched		= another_node.m_nlatched;
	m_min_in  		= another_node.m_min_in;
	m_max_in		= another_node.m_max_in;
	m_min_out 		= another_node.m_min_out;
	m_max_out		= another_node.m_max_out;
	m_assigned_in 	= another_node.m_assigned_in;
	m_assigned_out 	= another_node.m_assigned_out ;
	m_target_fanout	= another_node.m_target_fanout;
	m_available_in	= another_node.m_available_in;
	m_available_out	= another_node.m_available_out;

	m_nInter_cluster_inputs 	= another_node.m_nInter_cluster_inputs;
	m_nInter_cluster_outputs	= another_node.m_nInter_cluster_outputs;

	m_kin						= another_node.m_kin;
	return (*this);
}

void LEVEL_NODE::set_nNodes(const NUM_ELEMENTS & number_of_nodes)
{
	m_nNodes = number_of_nodes;
	m_nodes.resize(m_nNodes, 0);
}

// RETURN: the total weight for all the edges
//
NUM_ELEMENTS LEVEL_NODE::get_total_edge_weight(const EDGES& edges) const
{
	NUM_ELEMENTS total_weight_of_edges = 0;
	EDGES::const_iterator edge_iter;
	EDGE * edge = 0;

	for (edge_iter = edges.begin(); edge_iter != edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge->get_weight() >= 0);
		total_weight_of_edges += edge->get_weight();
	}
	
	return total_weight_of_edges;
}

//
// RETURN: true if the sum of the output edges is equal to what the output shape
//         says it should be via m_assigned_out
//
bool LEVEL_NODE::is_output_degree_satisfied() const
{
	bool output_degree_satisfied;
	NUM_ELEMENTS total_weight_of_edges = 0;

	debugif(DSANITY_CHECK, "Checking to see if output degree satisfied");

	total_weight_of_edges = m_nIntra_cluster_outputs + m_nInter_cluster_outputs;

	output_degree_satisfied = (total_weight_of_edges == m_assigned_out);


	if (output_degree_satisfied)
	{
		debugif(DSANITY_CHECK, "Output degree satisfied");
	}
	else
	{
		debugif(DSANITY_CHECK, "Output degree not satisfied");
	}
	
	return output_degree_satisfied;
}

//
// RETURN: true if the sum of the input edges is equal to what the input shape
//         says it should be via m_assigned_in
//
bool LEVEL_NODE::is_input_degree_satisfied() const
{
	bool input_degree_satisfied = false;
	NUM_ELEMENTS total_weight_of_input_edges = 0;

	debugif(DSANITY_CHECK, "Checking to see if input degree satisfied");

	total_weight_of_input_edges = m_nIntra_cluster_inputs + m_nInter_cluster_inputs;

	input_degree_satisfied = (total_weight_of_input_edges == m_assigned_in);

	return input_degree_satisfied;
}

//
// RETURN: true if the the number of delay defining edes input into this node 
//         is equal to what it should be - the number of indiv. nodes
//
bool LEVEL_NODE::is_delay_satisfied() const
{
	bool delay_satisfied  = false;
	NUM_ELEMENTS total_weight_of_length_1_edges;

	if (m_delay_level == 0)
	{
		// pi or dff do not need their delay defined
		return true;
	}

	debugif(DSANITY_CHECK,"checking to see if delay requirements satisfied");

	total_weight_of_length_1_edges = m_intra_cluster_input_edge_lengths[1] +
									m_inter_cluster_input_edge_lengths[1];			
	delay_satisfied = (total_weight_of_length_1_edges >= m_nNodes);

	if (delay_satisfied)
	{
		debugif(DSANITY_CHECK,"Delay satisfied");
	}
	else
	{
		debug("Delay not satisfied");
		debug("Wanted at least " << m_nNodes << " got " << total_weight_of_length_1_edges);
	}

	return delay_satisfied;
}

// RETURN: the number of input edges above what the input shape specifies
// 
NUM_ELEMENTS LEVEL_NODE::get_input_to_spare() const
{
	NUM_ELEMENTS input_to_spare = 0;
	NUM_ELEMENTS input_difference = get_input_difference();

	input_to_spare = MAX(0, input_difference);

	//debug("\tInput to spare " << input_to_spare);
	
	return input_to_spare;
}

// RETURN: the difference between how many input edges have been assigned and the 
//         the input shape (m_assigned_in)
// 
NUM_ELEMENTS LEVEL_NODE::get_input_difference() const
{
	NUM_ELEMENTS input_difference = 0;
	NUM_ELEMENTS total_weight_of_input_edges = 0;

	total_weight_of_input_edges = m_nIntra_cluster_inputs + m_nInter_cluster_inputs;

	input_difference = (total_weight_of_input_edges - m_assigned_in);

	//debug("\tInput to difference " << input_to_difference);
	
	return input_difference;
}

// RETURN: the number of output edges above what the output shape specifies
// 
NUM_ELEMENTS LEVEL_NODE::get_output_to_spare() const
{
	NUM_ELEMENTS output_to_spare = 0;
	NUM_ELEMENTS output_difference = get_output_difference();

	output_to_spare = MAX(0, output_difference);

	//debug("\tOutput to spare " << output_to_spare);
	
	return output_to_spare;
}

// RETURN: the difference between how many output edges have been assigned and the 
//         the output shape (m_assigned_in)
// 
NUM_ELEMENTS LEVEL_NODE::get_output_difference() const
{
	NUM_ELEMENTS output_difference;
	NUM_ELEMENTS total_weight_of_edges = 0;

	total_weight_of_edges = m_nIntra_cluster_outputs + m_nInter_cluster_outputs;


	output_difference = (total_weight_of_edges - m_assigned_out);


	//debug("\tOutput to difference " << output_to_difference);
	
	return output_difference;
}

// RETURN: the number of extra length 1 edges that this level node
//         does not need
NUM_ELEMENTS LEVEL_NODE::get_delay_to_spare() const
{
	NUM_ELEMENTS delay_to_spare  = 0;
	NUM_ELEMENTS delay_difference  = get_delay_difference();


	delay_to_spare = MAX(0, delay_difference);

	//debug("\tDelay to spare " << delay_to_spare);

	return delay_to_spare;
}

// RETURN: the difference between the number of length 1 edge the level 
//         node needs to define the delay level of its indiv. nodes and the 
//         number of length 1 edges that it has being input into it
NUM_ELEMENTS LEVEL_NODE::get_delay_difference() const
{

	NUM_ELEMENTS delay_difference  = 0;
	NUM_ELEMENTS total_weight_of_length_1_edges = 0;

	if (m_delay_level == 0)
	{
		// pi or dff do not need their delay defined
		return 0;
	}

	total_weight_of_length_1_edges = m_intra_cluster_input_edge_lengths[1] +
									m_inter_cluster_input_edge_lengths[1];			
	delay_difference = (total_weight_of_length_1_edges - m_nNodes);

	//debug("\tDelay difference" << delay_difference);

	assert(delay_difference >= 0);

	return delay_difference;

}

// return an edge with the given source level node 
//
// RETURN: an edge connecting from source_node to this level node if one could be found
//         otherwise 0
//
EDGE * LEVEL_NODE::find_edge_with_source
(
	const LEVEL_NODE * source_node
) const
{
	// draw a picture to understand
	assert(source_node);
	EDGE * edge = 0;
	NUM_ELEMENTS edge_position = 0;
	bool edge_found = false;
	EDGES::const_iterator edge_iter;

	if (m_delay_level > 0)
	{
		CLUSTER_NUMBER_TYPE source_cluster = source_node->get_cluster_number();
		DELAY_TYPE source_delay_level = source_node->get_delay_level();

		edge_position = source_cluster * m_delay_level + source_delay_level;

		assert(edge_position >= 0 && static_cast<unsigned>(edge_position) < m_input_edges.size());
		edge = m_input_edges[edge_position];
		assert(edge->get_source_level_node() == source_node);
	}
	else
	{
		for (edge_iter = m_input_to_dff_edges.begin(); ! edge_found && edge_iter != m_input_to_dff_edges.end();
		 	edge_iter++)
		{
			edge = *edge_iter;
			assert(edge);

			edge_found = (edge->get_source_level_node() == source_node);
		}

		if (! edge_found)
		{
			edge = 0;
		}
	}
				

	return edge;
}
// return an edge with the given source level node 
//
// PRE: we have assigned input edges to this level node in a particular way
//      source_delay_level > m_delay_level
// RETURN: an edge connecting from source_node to this level node 
//
EDGE * LEVEL_NODE::find_edge_with_source
(
	const CLUSTER_NUMBER_TYPE& source_cluster, 
	const DELAY_TYPE& source_delay_level
) const
{
	// draw a picture to understand

	NUM_ELEMENTS edge_position = source_cluster * m_delay_level + source_delay_level;

	assert(edge_position >= 0 && static_cast<unsigned>(edge_position) < m_input_edges.size());
	EDGE * edge = m_input_edges[edge_position];
	//assert(edge->get_source_cluster_number() == source_cluster);
	//assert(edge->get_source_level_node()->get_delay_level() == source_delay_level);
	//assert(edge->get_sink_level_node() == this);

	return edge;
}


// POST: we have added this input edge to this level node 
void LEVEL_NODE::add_input_edge(EDGE * edge)
{
	assert(edge);

	LEVEL_NODE * source_level_node = edge->get_source_level_node();
	assert(source_level_node);
	CLUSTER_NUMBER_TYPE source_cluster = source_level_node->get_cluster_number();
	DELAY_TYPE source_delay_level = source_level_node->get_delay_level();

	NUM_ELEMENTS edge_position = source_cluster * m_delay_level + source_delay_level;

	assert(edge_position >= 0 && static_cast<unsigned>(edge_position) < m_input_edges.size());
	m_input_edges[edge_position] = edge;
}
void LEVEL_NODE::add_output_edge(EDGE * edge)
{ 
	assert(edge);
	m_output_edges.push_back(edge);
}

// an input edge has changed. update data members to reflect this
//
void LEVEL_NODE::input_edge_has_changed
(
	const EDGE * edge,
	const NUM_ELEMENTS& weight_addition
)
{
	assert(edge && m_cluster);

	LENGTH_TYPE edge_length = edge->get_length();
	assert(edge_length >= 0);

	if (edge->is_intra_cluster())
	{
		m_intra_cluster_input_edge_lengths[edge_length] += weight_addition;
		m_nIntra_cluster_inputs += weight_addition;

	}
	else
	{
		assert(edge->is_inter_cluster());

		m_inter_cluster_input_edge_lengths[edge_length] += weight_addition;
		m_nInter_cluster_inputs += weight_addition;

		m_cluster->inter_cluster_input_edge_has_changed(edge, weight_addition);
	}

	m_available_in += weight_addition;

	assert(m_nIntra_cluster_inputs >= 0 && m_nInter_cluster_inputs >= 0);


	/*

	Sanity checks

	NUM_ELEMENTS intra_edge_sum = 0;
	NUM_ELEMENTS inter_edge_sum = 0;
	
	intra_edge_sum = accumulate(m_intra_cluster_input_edge_lengths.begin(), 
							m_intra_cluster_input_edge_lengths.end(), 0);
	assert(intra_edge_sum == m_nIntra_cluster_inputs);

	inter_edge_sum = 0;
	inter_edge_sum = accumulate(m_inter_cluster_input_edge_lengths.begin(), 
							m_inter_cluster_input_edge_lengths.end(), 0);
	assert(inter_edge_sum == m_nInter_cluster_inputs);

	EDGES::const_iterator edge_iter;
	NUM_ELEMENTS input_sum = 0;

	if (m_delay_level > 0)
	{
		for (edge_iter = m_input_edges.begin(); edge_iter != m_input_edges.end(); edge_iter++)
		{
			if ((*edge_iter) != 0)
			{
				input_sum += (*edge_iter)->get_weight();
			}
		}
		assert(input_sum == inter_edge_sum + intra_edge_sum);
	}
	*/
}

// an output edge has changed. update data members to reflect this
//
void LEVEL_NODE::output_edge_has_changed
(
	const EDGE * edge,
	const NUM_ELEMENTS& weight_addition
)
{
	assert(edge && m_cluster);
	LENGTH_TYPE edge_length = edge->get_length();
	assert(edge_length >= 0);

	if (edge->is_intra_cluster())
	{
		m_intra_cluster_output_edge_lengths[edge_length] += weight_addition;
		m_nIntra_cluster_outputs += weight_addition;
		m_cluster->intra_cluster_edge_has_changed(edge, weight_addition);
	}
	else
	{
		assert(edge->is_inter_cluster());

		m_inter_cluster_output_edge_lengths[edge_length] += weight_addition;
		m_nInter_cluster_outputs += weight_addition;
		m_cluster->inter_cluster_output_edge_has_changed(edge, weight_addition);
	}

	m_available_out += weight_addition;
	assert(m_nIntra_cluster_outputs >= 0 && m_nInter_cluster_outputs >= 0);

	/*
	 *
	 * Sanity checks
	 *
	NUM_ELEMENTS edge_sum = 0;
	edge_sum = accumulate(m_intra_cluster_output_edge_lengths.begin(), 
							m_intra_cluster_output_edge_lengths.end(), 0);
	assert(edge_sum == m_nIntra_cluster_outputs);
	edge_sum = accumulate(m_inter_cluster_output_edge_lengths.begin(), 
							m_inter_cluster_output_edge_lengths.end(), edge_sum);
	NUM_ELEMENTS output_edge_weight = get_output_edge_weight();

	assert(output_edge_weight == edge_sum);
	*/
}

void LEVEL_NODE::report_edge_lengths() const
{
	cout << "\t\tIntra cluster input_edge_length\n\t\t";
	copy(m_intra_cluster_input_edge_lengths.begin(), m_intra_cluster_input_edge_lengths.end(), 
			ostream_iterator<NUM_ELEMENTS>(cout, " ")); cout << "\n";

	cout << "\t\tIntra cluster output_edge_length\n\t\t";
	copy(m_intra_cluster_output_edge_lengths.begin(), m_intra_cluster_output_edge_lengths.end(), 
			ostream_iterator<NUM_ELEMENTS>(cout, " ")); cout << "\n";

	cout << "\t\tInter cluster input_edge_length\n\t\t";
	copy(m_inter_cluster_input_edge_lengths.begin(), m_inter_cluster_input_edge_lengths.end(), 
			ostream_iterator<NUM_ELEMENTS>(cout, " ")); cout << "\n";

	cout << "\t\tInter cluster output_edge_length\n\t\t";
	copy(m_inter_cluster_output_edge_lengths.begin(), m_inter_cluster_output_edge_lengths.end(), 
			ostream_iterator<NUM_ELEMENTS>(cout, " ")); cout << endl;
}

// POST: We are sure that:
//        m_nIntra_cluster_inputs,m_nInter_cluster_inputs,
//        m_nIntra_cluster_outputs,m_nInter_cluster_outputs 
//        m_available_in, m_available_out
//        m_intra_cluster_input_edge_lengths,
//        m_intra_cluster_output_edge_lengths,
//        m_intra_cluster_output_edge_lengths,
//        m_inter_cluster_output_edge_lengths
//
//        are all consistent with the edges input and output out of this level node
//
//        Also:
//
//        m_nIntra_cluster_inputs + m_nInter_cluster_inputs = m_available_in
//        m_nIntra_cluster_outputs + m_nInter_cluster_outputs = m_available_out
//
//        delay is satisfied for our indiv. nodes
//
void LEVEL_NODE::sanity_check() const
{
	assert(m_cluster);
	CLUSTER_SPEC * comb_spec = m_cluster->get_comb_spec();
	assert(comb_spec);
	DELAY_TYPE delay = comb_spec->get_delay();

	NUM_ELEMENTS nIntra_cluster_inputs = 0, 
				nIntra_cluster_outputs = 0,
				nInter_cluster_inputs = 0,
				nInter_cluster_outputs = 0,
				available_in = 0, 
				available_out = 0;

	DISTRIBUTION intra_cluster_input_edge_lengths(delay +1, 0);
	DISTRIBUTION intra_cluster_output_edge_lengths(delay + 1, 0);
	DISTRIBUTION inter_cluster_input_edge_lengths(delay + 1, 0); 
	DISTRIBUTION inter_cluster_output_edge_lengths(delay + 1, 0);

	assert(m_delay_level <= delay);

	EDGES::const_iterator edge_iter;
	EDGE * edge = 0;
	LENGTH_TYPE length = 0;
	NUM_ELEMENTS weight = 0;

	// find the input distributions
	for (edge_iter = m_input_edges.begin(); edge_iter != m_input_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		length = edge->get_length();
		assert(length >= 0 && length <= m_delay_level);

		weight = edge->get_weight();
		assert(weight >= 0);

		available_in += weight;

		if (edge->is_intra_cluster())
		{
			nIntra_cluster_inputs += weight;
			intra_cluster_input_edge_lengths[length] += weight;
		}
		else
		{
			assert(edge->is_inter_cluster());

			nInter_cluster_inputs += weight;
			inter_cluster_input_edge_lengths[length] += weight;
		}
	}

	if (m_delay_level == 0)
	{
		// find the input distributions
		for (edge_iter = m_input_to_dff_edges.begin(); edge_iter != m_input_to_dff_edges.end(); edge_iter++)
		{
			edge = *edge_iter;
			assert(edge);

			length = edge->get_length();
			assert(length >= 0 && length <= m_delay_level);

			weight = edge->get_weight();
			assert(weight >= 0);

			available_in += weight;

			if (edge->is_intra_cluster())
			{
				nIntra_cluster_inputs += weight;
				intra_cluster_input_edge_lengths[length] += weight;
			}
			else
			{
				assert(edge->is_inter_cluster());

				nInter_cluster_inputs += weight;
				inter_cluster_input_edge_lengths[length] += weight;
			}
		}
	}

	// find the output distributions
	for (edge_iter = m_output_edges.begin(); edge_iter != m_output_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		length = edge->get_length();
		assert(length >= 0 && length <= delay-m_delay_level);

		weight = edge->get_weight();
		assert(weight >= 0);

		available_out += weight;

		if (edge->is_intra_cluster())
		{
			nIntra_cluster_outputs += weight;
			intra_cluster_output_edge_lengths[length] += weight;
		}
		else
		{
			assert(edge->is_inter_cluster());

			nInter_cluster_outputs += weight;
			inter_cluster_output_edge_lengths[length] += weight;
		}
	}
	// check our state variables

	assert(nIntra_cluster_inputs == m_nIntra_cluster_inputs);
	assert(nInter_cluster_inputs == m_nInter_cluster_inputs);
	assert(nIntra_cluster_outputs == m_nIntra_cluster_outputs);
	assert(nInter_cluster_outputs == m_nInter_cluster_outputs);
	assert(available_in == m_available_in);
	assert(available_out == m_available_out);
	assert(m_available_in == m_nIntra_cluster_inputs + m_nInter_cluster_inputs);
	assert(m_available_out == m_nIntra_cluster_outputs + m_nInter_cluster_outputs);
	assert(m_assigned_in == comb_spec->get_level_in()[m_delay_level]);
	assert(m_assigned_out == comb_spec->get_level_out()[m_delay_level]);

	for (length = 0; length <= delay; length++)
	{
		assert(intra_cluster_input_edge_lengths[length] == m_intra_cluster_input_edge_lengths[length]);
		assert(inter_cluster_input_edge_lengths[length] == m_inter_cluster_input_edge_lengths[length]);
		assert(intra_cluster_output_edge_lengths[length] == m_intra_cluster_output_edge_lengths[length]);
		assert(inter_cluster_output_edge_lengths[length] == m_inter_cluster_output_edge_lengths[length]);
	}

	if (m_delay_level != 0)
	{
		assert(m_assigned_in + get_input_difference() == get_input_edge_weight());
		assert(m_assigned_out + get_output_difference() == get_output_edge_weight());
	}

	// for now we don't add 0 length edges to our distributions
	assert(intra_cluster_output_edge_lengths[0] == 0);
	assert(inter_cluster_output_edge_lengths[0] == 0);


	if (m_nNodes == 0)
	{
		assert(m_available_in == 0);
		assert(m_available_out == 0);
	}


	assert(is_delay_satisfied());
}

		
// the cost for the combinational delay graph structure iterative algorithm
// note: we don't have to cost the double input cost because
// (the double input case is where we will force an invidual node to have 
//  two edges from the same node during final edge assignment)
// that is taken care of during edge creation when we set the maximum
// weight that an edge can attain.
//
// RETURNS: cost_Level_Shape + cost_Problem_Node 
//          for the level node
//
COST_TYPE LEVEL_NODE::get_cost() const
{

	// cost_Level_Shape = get_input_cost() + get_output_cost()
	//
	// cost_Problem_Node = 	get_too_many_inputs_cost() + 
	//					   	get_too_few_outputs_cost() + 
	//						get_too_few_inputs_cost();

	COST_TYPE cost = get_input_cost() + 
					 get_output_cost() + 
					 get_too_many_inputs_cost() + 
					 get_too_few_outputs_cost() + 
					 get_too_few_inputs_cost();
					 

	return cost;
}

// RETURN: the number of outputs still needed * congestion factor
COST_TYPE LEVEL_NODE::get_too_few_outputs_cost() const 
{
	return m_too_few_outputs_factor * get_nOutputs_still_needed();
}

// RETURN: difference between the output weight and what it should be
//
COST_TYPE LEVEL_NODE::get_output_cost() const
{
	NUM_ELEMENTS total_weight_of_output_edges = get_output_edge_weight();
	COST_TYPE cost = ABS(total_weight_of_output_edges - m_assigned_out);

	return cost;
}



// RETURN: difference between the input weight and what it should be
//
COST_TYPE LEVEL_NODE::get_input_cost() const
{	
	NUM_ELEMENTS total_weight_of_input_edges = get_input_edge_weight();
	COST_TYPE cost = 0;

	if (m_delay_level > 0)
	{
		cost = ABS(total_weight_of_input_edges - m_assigned_in);
	}
	else
	{
		cost = 0;
	}

	return cost;
}
// RETURN: the number of excessive inputs  * congestion factor
COST_TYPE LEVEL_NODE::get_too_many_inputs_cost() const
{
	COST_TYPE cost = 0;
	NUM_ELEMENTS excess_number_of_inputs = get_excess_number_of_inputs();

	if (excess_number_of_inputs > 0)
	{
		cost = m_too_many_inputs_factor * excess_number_of_inputs;
	}

	return cost;
}

// Increase/Decrease the congestion factors for node violations

void LEVEL_NODE::increase_too_many_inputs_factor()
{
	m_too_many_inputs_factor *= m_mult_too_many_inputs_factor;
}
void LEVEL_NODE::increase_too_few_outputs_factor()
{
	m_too_few_outputs_factor *= m_mult_too_few_outputs_factor;
}
void LEVEL_NODE::increase_too_few_inputs_factor()
{
	m_too_few_inputs_factor *= m_mult_too_few_inputs_factor;
}
void LEVEL_NODE::decrease_too_many_inputs_factor()
{
	m_too_many_inputs_factor /= 2;
	m_too_many_inputs_factor = MAX(1, m_too_many_inputs_factor);
}
void LEVEL_NODE::decrease_too_few_outputs_factor()
{
	m_too_few_outputs_factor /= 2;
	m_too_few_outputs_factor = MAX(1, m_too_few_outputs_factor);
}
void LEVEL_NODE::decrease_too_few_inputs_factor()
{
	m_too_few_inputs_factor /= 2;
	assert(new_gen);
	m_too_few_inputs_factor = MAX(1, m_too_few_inputs_factor);
}


// RETURN: the number of inputs needed * congestion factor
COST_TYPE LEVEL_NODE::get_too_few_inputs_cost() const
{

	COST_TYPE total_input_weight = get_input_edge_weight(),
			  not_enough_input = 0;

	// each node should have at least two inputs unless
	// we should have less than two inputs per node 
	
	if (m_delay_level != 0 && total_input_weight < m_assigned_in)
	{
		not_enough_input = 2*m_nNodes - total_input_weight;
		not_enough_input = MAX(not_enough_input, m_assigned_in-total_input_weight);
		not_enough_input = MAX(not_enough_input, 0);
		not_enough_input = m_too_few_inputs_factor*not_enough_input;
	}
	else
	{
		// right now do nothing.
	}

	return not_enough_input;
}

// POST: the congestion costs have been reset to their default values
void LEVEL_NODE::reset_congestion_costs()
{
	assert(g_options);
	m_too_many_inputs_factor = g_options->get_delay_structure_init_too_many_inputs_factor();
	m_too_few_outputs_factor = g_options->get_delay_structure_init_too_few_outputs_factor();
	m_too_few_inputs_factor  = g_options->get_delay_structure_init_too_few_inputs_factor();
}

// RETURN: the number of excessive number of inputs
//         which for all delay levels but the 0th is 
//         the difference between the input edge weight and
//         kin* the number of nodes
NUM_ELEMENTS LEVEL_NODE::get_excess_number_of_inputs() const
{
	NUM_ELEMENTS total_weight_of_input_edges = get_input_edge_weight();
	NUM_ELEMENTS excess_number_of_inputs = 0;

	if (m_delay_level != 0)
	{
		excess_number_of_inputs = total_weight_of_input_edges - m_kin*m_nNodes;
		excess_number_of_inputs = MAX(excess_number_of_inputs, 0);
	}
	else
	{
		excess_number_of_inputs = 0;
	}

	return excess_number_of_inputs;
}
// RETURN: get the number of inputs that will be available 
//         during final edge assignment based on the input edge weight
NUM_ELEMENTS LEVEL_NODE::get_nOpen_inputs() const
{
	NUM_ELEMENTS total_weight_of_input_edges = get_input_edge_weight();
	NUM_ELEMENTS open_inputs = 0;

	if (m_delay_level != 0)
	{
		open_inputs = m_kin*m_nNodes - total_weight_of_input_edges;
		open_inputs = MAX(open_inputs, 0);
	}
	else
	{
		open_inputs = 0;
	}

	return open_inputs;
}




// RETURN: the number of output edge weight that is needed to ensure that
//          each indiv. node will have enough inputs during final edge 
//          assignment
NUM_ELEMENTS LEVEL_NODE::get_nOutputs_still_needed() const
{
	NUM_ELEMENTS total_output_weight = get_output_edge_weight();
	NUM_ELEMENTS not_enough_output = 0;

	// each node should have at least one output unless it is a latch or a PO
	if (m_delay_level != m_cluster->get_delay())
	{
		not_enough_output = m_nNodes - m_nlatched - m_nPO - total_output_weight;
		not_enough_output = MAX(not_enough_output, 0);
	}
	else
	{
		not_enough_output = 0;
	}

	return not_enough_output;
}

// need at least 1 input for every comb node
// want 2 but need at least 1
//
// RETURN: the number of inputs we still want based on giving each 
//         combinational node at least 2 inputs so that 
//         there won't be buffers
NUM_ELEMENTS LEVEL_NODE::get_nInputs_still_wanted() const
{
	// for right now don't look at the 0thy delay level
	assert(new_gen);


	NUM_ELEMENTS total_input_weight = get_input_edge_weight();
	NUM_ELEMENTS not_enough_input = 0;

	// each node should have at least two inputs
	// unless the specification says it shouldn't.
	if (m_delay_level != 0 && total_input_weight < m_assigned_in)
	{
		not_enough_input = 2*m_nNodes - total_input_weight;
		not_enough_input = MAX(not_enough_input, (m_assigned_in-total_input_weight));
		not_enough_input = MAX(not_enough_input, 0);
	}
	else
	{
		// we should do something
		//not_enough_input = m_nDFF - total_input_weight;
	}

	return not_enough_input;
}
// need at least 1 input for every comb node
// want 2 but need at least 1
//
// RETURN: the number of inputs we still need to make sure
//         each node has at least one input
NUM_ELEMENTS LEVEL_NODE::get_nInputs_still_needed() const
{
	// for right now don't look at the 0thy delay level
	assert(new_gen);


	NUM_ELEMENTS total_input_weight = get_input_edge_weight();
	NUM_ELEMENTS not_enough_input = 0;

	// each node should have at least two inputs
	// but right now we will settle for 1
	if (m_delay_level != 0)
	{
		not_enough_input = 1*m_nNodes - total_input_weight;
		not_enough_input = MAX(not_enough_input, 0);
	}
	else
	{
		// we should do something
		//not_enough_input = m_nDFF - total_input_weight;
	}

	return not_enough_input;
}


// RETURN: the node violation cost related to having not enough
//         input edge weight so that during final edge assignment we
//         will have nodes with no inputs
//         
COST_TYPE LEVEL_NODE::get_nIndiv_nodes_too_few_inputs_cost() const
{
	NODES::const_iterator node_iter;
	NODE * node = 0;
	COST_TYPE inputs_needed = 0.0;

	for (node_iter = m_nodes.begin(); node_iter != m_nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		inputs_needed += node->get_nInputs_still_needed();
	}

	return inputs_needed;
}

// RETURN: the number of fanouts with value degree
NUM_ELEMENTS& LEVEL_NODE::OutDegrees
(
	const NUM_ELEMENTS& degree
)
{
	assert(degree >= 0 && static_cast<unsigned>(degree) < m_out_degrees.size());
	return m_out_degrees[degree];
}

// RETURN: the wirelength-approx cost contribution from this level node
COST_TYPE LEVEL_NODE::get_wirelength_cost() const
{
	COST_TYPE total_wirelength_cost = 0;
	NODES::const_iterator node_iter;
	NODE * node = 0;
	
	for (node_iter = m_nodes.begin(); node_iter != m_nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		total_wirelength_cost += node->get_wirelength_cost();
	}

	return total_wirelength_cost;



}
// PRE: degree > max fanout
// POST: m_out_degrees has been expanded
void LEVEL_NODE::set_new_max_fanout
(
	const NUM_ELEMENTS& degree
)
{
	assert(degree >= 0 && static_cast<unsigned>(degree) >= m_out_degrees.size());

	m_out_degrees.resize(degree+1, 0);
}


// RETURN: the node at position node_index
NODE* LEVEL_NODE::get_node
(
	const NUM_ELEMENTS& node_index
) const
{ 
	assert(node_index >= 0 && static_cast<unsigned>(node_index) < m_nodes.size());

	return m_nodes[node_index];
}

// POST: an edge that inputs to a flip-flop has been added 
void LEVEL_NODE::add_dff_input_edge(EDGE * edge)
{
	assert(edge);
	assert(m_delay_level == 0);

	m_input_to_dff_edges.push_back(edge);
}

// POST: an edge that outputs to a flip-flop has been added
void LEVEL_NODE::add_output_to_dff_edge(EDGE * edge)
{
	assert(edge);

	m_output_to_dff_edges.push_back(edge);
}

// PRE: we have not assigned any edges yet.
// RETURN: the number of long or inter-cluster edges that each indiv. node 
//         in this level node will connect to
DISTRIBUTION LEVEL_NODE::get_number_of_unassigned_long_and_inter_cluster_edges_each_indiv_node_has() const
{
	DISTRIBUTION long_distribution;
	NUM_ELEMENTS nLong;
	NODES::const_iterator node_iter;
	NODE * node = 0;

	//debug("Geting the long distribution");
	for (node_iter = m_nodes.begin(); node_iter != m_nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		nLong = node->get_unassigned_long_and_inter_cluster_edges();

		//debug("Node " << node->get_id() << " has nLong = " << nLong);

		long_distribution.push_back(nLong);
	}

	return long_distribution;
}

// RETURN: the PI nodes
NODES LEVEL_NODE::get_PI() const
{
	NODES::const_iterator node_iter;
	NODE * node = 0;
	NODES PI_nodes;

	for (node_iter = m_nodes.begin(); node_iter != m_nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		if (node->is_PI())
		{
			PI_nodes.push_back(node);
		}
	}

	return PI_nodes;
}

// RETURN: the nodes with POs
NODES LEVEL_NODE::get_PO() const
{
	NODES::const_iterator node_iter;
	NODE * node = 0;
	NODES PO_nodes;

	for (node_iter = m_nodes.begin(); node_iter != m_nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		if (node->is_PO())
		{
			PO_nodes.push_back(node);
		}
	}

	return PO_nodes;
}

// RETURN: the nodes with DFFs
NODES LEVEL_NODE::get_DFF() const
{
	NODES::const_iterator node_iter;
	NODE * node = 0;
	NODES DFF_nodes;

	for (node_iter = m_nodes.begin(); node_iter != m_nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		if (node->is_DFF())
		{
			DFF_nodes.push_back(node);
		}
	}

	assert(static_cast<unsigned>(m_nDFF) == DFF_nodes.size());

	return DFF_nodes;
}
// RETURN: the latched nodes that have not been assigned
NODES LEVEL_NODE::get_latched_nodes_not_taken(const NUM_ELEMENTS& number_to_get) const
{
	NODES::const_iterator node_iter;
	NODE * node = 0;
	NODES latched_nodes;
	NUM_ELEMENTS nLatch_nodes_to_obtain = number_to_get;

	for (node_iter = m_nodes.begin(); nLatch_nodes_to_obtain > 0 && node_iter != m_nodes.end(); 
		node_iter++)
	{
		node = *node_iter;
		assert(node);

		if (node->is_latched() && ! node->is_latch_taken())
		{
			latched_nodes.push_back(node);
			nLatch_nodes_to_obtain--;
		}
	}

	assert(nLatch_nodes_to_obtain == 0);

	return latched_nodes;
}
NODES LEVEL_NODE::get_latched_nodes() const
{
	NODES::const_iterator node_iter;
	NODE * node = 0;
	NODES latched_nodes;

	for (node_iter = m_nodes.begin(); node_iter != m_nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		if (node->is_latched())
		{
			latched_nodes.push_back(node);
		}
	}

	return latched_nodes;
}


//
//
// PRE: the latched nodes have been designated
// POST: the nodes with PO have been designated
void LEVEL_NODE::designate_PO()
{
	NODES::const_iterator node_iter;
	NODE * node = 0;
	NODES DFF_nodes;
	NUM_ELEMENTS nPO = m_nPO;
	NUM_ELEMENTS node_chosen;
	NUM_ELEMENTS po_count = 0;

	// we will treat the 0th delay level specially because
	// the combinational nodes at the 0th delay level cannot be PO
	// because it does not make sense
	
	// for delay_levels > 0
	// find all the combinational non-latched 0 degree nodes and make them a primary output
	// for delay level 0 find all of the dff with 0 fanout that are not latched
	// and make them po

	for (node_iter = m_nodes.begin(); node_iter != m_nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		if (node->get_fanout_degree() == 0 && ! node->is_latched())
		{
			po_count++;
			node->set_is_PO(true);
		}
	}

	nPO -= po_count;
	
	// if we still have PO to designate, designate them randomly
	while (nPO > 0)
	{
		node_chosen = g_rand_number_gen->random_number(m_nNodes-1);
		node = m_nodes[node_chosen];
		assert(node);

		// if we haven't assigned the node a PO and it isn't a combinational primary input
		if (! node->is_PO() && ! node->is_PI())
		{
			node->set_is_PO(true);
			po_count++;
			nPO--;
		}
	}
}

// print what connects to the the indiv. nodes in this level node
void LEVEL_NODE::print_indiv_edge_input_connections() const
{
	NODES::const_iterator node_iter;
	NODE * node = 0;

	for (node_iter = m_nodes.begin(); node_iter != m_nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		debug("dst node " << node->get_name());
		debug("inputs to node " << node->get_input_node_names() << endl);
	}
}

// print what the super output edges connects to
void LEVEL_NODE::print_output_connections() const
{
	EDGES::const_iterator edge_iter;
	LEVEL_NODE * sink_level_node = 0;

	debug("Node " << m_name << " has " << m_nNodes << " nodes");
	debug("Node " << m_name << " has outputs from: ");

	for (edge_iter = m_output_edges.begin(); edge_iter != m_output_edges.end(); edge_iter++)
	{
		if ((*edge_iter)->get_weight() > 0)
		{
			sink_level_node = (*edge_iter)->get_sink_level_node();
			assert(sink_level_node);

			if ((*edge_iter)->get_weight() > 0)
			{
				debug(sink_level_node->get_name() << " of weight " << (*edge_iter)->get_weight());
				debug("Sink node has " << sink_level_node->get_nNodes() << " nodes");
			}
		}
	}

	dbg("\n");
}

// PRE: 
// POST: edges_for_reduction_pmf has the weight
//       that each input edge can give with the contraints being:
//       1) that for each edge the source level node should have enough weight
//          to be able to give each of its indiv. nodes at least one output and 
//       2) that the level node should still be able to define the delay level of its nodes
// RETURN: the total weight that can be reduced from these edges
//
NUM_ELEMENTS LEVEL_NODE::get_edges_ok_for_reduction
(
	DISTRIBUTION & edges_ok_for_reduction_pmf
)
{

	EDGES::iterator edge_iter;
	EDGE * edge = 0;
	LEVEL_NODE * source_level_node = 0;
	NUM_ELEMENTS input_sum = 0,
				unneeded_edges = 0,
				total_unneeded_edges = 0,
				excess_number_of_inputs = get_excess_number_of_inputs();


	edges_ok_for_reduction_pmf.clear();

	Verbose("Node " << m_name << " has " << m_nNodes << " nodes and " << excess_number_of_inputs <<
			" too many inputs.");
	//debug("It has inputs from: ");


	for (edge_iter = m_input_edges.begin(); edge_iter != m_input_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);
		unneeded_edges = 0;

		if (edge->get_weight() > 0)
		{
			input_sum += edge->get_weight();
			source_level_node = edge->get_source_level_node();
			assert(source_level_node);

			//debug(source_level_node->get_name() << " of weight " << edge->get_weight());

			unneeded_edges = edge->get_nUnneeded_edges();
			assert(unneeded_edges >= 0);

			total_unneeded_edges += unneeded_edges;

				
		}
		edges_ok_for_reduction_pmf.push_back(unneeded_edges);
	}

	assert(m_input_edges.size() == edges_ok_for_reduction_pmf.size());
	assert(input_sum == get_input_edge_weight());

	return total_unneeded_edges;
}


//
// PRE: 
// POST: edges_for_reduction_pmf has the weight
//       that each input edge can give with the only contraint being 
//       that the level node should still be able to define the delay level of its nodes
// RETURN: the total weight that can be reduced from these edges
NUM_ELEMENTS LEVEL_NODE::get_edges_for_reduction
(
	DISTRIBUTION & edges_for_reduction_pmf
) const
{

	EDGES::const_iterator edge_iter;
	EDGE * edge = 0;
	LEVEL_NODE * source_level_node = 0;
	NUM_ELEMENTS input_sum = 0,
				unneeded_edges = 0,
				total_unneeded_edges = 0,
				excess_number_of_inputs = get_excess_number_of_inputs();


	edges_for_reduction_pmf.clear();

	debug("Node " << m_name << " has " << m_nNodes << " nodes and " << excess_number_of_inputs <<
			" too many inputs.");
	//debug("It has inputs from: ");

	for (edge_iter = m_input_edges.begin(); edge_iter != m_input_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);
		unneeded_edges = 0;

		if (edge->get_weight() > 0)
		{
			input_sum += edge->get_weight();
			source_level_node = edge->get_source_level_node();
			assert(source_level_node);

			//debug(source_level_node->get_name() << " of weight " << edge->get_weight());

			unneeded_edges = edge->get_nEdge_to_sacrifice();
			assert(unneeded_edges >= 0);

			total_unneeded_edges += unneeded_edges;
		}
		edges_for_reduction_pmf.push_back(unneeded_edges);
	}

	assert(m_input_edges.size() == edges_for_reduction_pmf.size());
	assert(input_sum == get_input_edge_weight());

	return total_unneeded_edges;
}

// try to reduce weight on the super edges that input into this level node so that 
// we won't have node violations during simple edge assignment
//
// RETURN: true if succesful in reducing weight by the amount needed to void 
//         double input node violations
bool LEVEL_NODE::reduce_weight_on_input_edges_to_eliminate_double_edge_connections()
{
	EDGE * edge = 0;
	DISTRIBUTION edges_ok_for_reduction_pmf;
	NUM_ELEMENTS edge_chosen = 0,
				total_unneeded_edges = 0,
				excess_number_of_inputs = get_excess_number_of_inputs();
	CIRCUIT_GRAPH * circuit = m_cluster->get_my_circuit();
	assert(circuit);

	assert(excess_number_of_inputs > 0);

	Verbose("Node " << m_name << " has " << m_nNodes << " nodes and " << excess_number_of_inputs <<
			" too many inputs.");
	//debug("It has inputs from: ");

	total_unneeded_edges = get_edges_ok_for_reduction(edges_ok_for_reduction_pmf);
	

	assert(new_gen); // clean up this code
	if (total_unneeded_edges < excess_number_of_inputs)
	{
		Verbose("Node " << m_name << " does not have enough edges to sacrifice nicely." 
				<< " Trying to remove all but the delay edges");
		total_unneeded_edges = get_edges_for_reduction(edges_ok_for_reduction_pmf);
		if (total_unneeded_edges < excess_number_of_inputs)
		{
			Verbose("Node " << m_name << " is screwed.");
			return false;
		}
	}
	assert(total_unneeded_edges >= excess_number_of_inputs);

	Verbose("\nNeed " << excess_number_of_inputs << 
			" fewer input edges. Have " << total_unneeded_edges << " to sacrifice");
	while (excess_number_of_inputs > 0)
	{
		edge_chosen = g_rand_number_gen->discrete_pmf(edges_ok_for_reduction_pmf);
		assert(edges_ok_for_reduction_pmf[edge_chosen] > 0);

		edge = m_input_edges[edge_chosen];
		assert(edge);

		Verbose("  Reducing weight on an edge");
		edge->add_weight(-1);
		edges_ok_for_reduction_pmf[edge_chosen]--;
		excess_number_of_inputs--;

		if (edge->get_length() == 1)
		{
			// recalculate the number of unneeded edges
			// because reduction of a length 1 edges may 
			// cause a reduction of the number of unneeded edges of 
			// other length 1 edges

			total_unneeded_edges = get_edges_ok_for_reduction(edges_ok_for_reduction_pmf);

			if (total_unneeded_edges < excess_number_of_inputs)
			{
				return false;
			}
		}

		// add the edges that lost weight to the list of edges that needed to be deleted
		// or have become "lost"
		circuit->add_lost_weight_edge(edge);
	}

	dbg("\n");

	return true;
}

// RETURN: the number of indiv. edges that input into 
//         the indiv. nodes of this level node
NUM_ELEMENTS LEVEL_NODE::get_nSimple_input_edges() const
{
	NODES::const_iterator node_iter;
	NODE * node = 0;
	NUM_ELEMENTS nInputs = 0;

	for (node_iter = m_nodes.begin(); node_iter != m_nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		if (node->is_PI())
		{
			assert(node->get_nInputs() == 0);
		}

		nInputs += node->get_nInputs();
	}

	return nInputs;
}

// RETURN: the number of simple edges that connect to combinational nodes
NUM_ELEMENTS LEVEL_NODE::get_nSimple_comb_output_edges() const
{
	NODES::const_iterator node_iter;
	NODE * node = 0;
	NUM_ELEMENTS nOutputs = 0;

	for (node_iter = m_nodes.begin(); node_iter != m_nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		nOutputs += node->get_nComb_Outputs();
	}

	return nOutputs;
}

// find the largest possible degree that can be assigned to a indiv node of this level 
// node.
//
// the largest possible degree is the number of nodes that 
// can be reached by the level node
//
//     o
//     | \  (can reach 7 nodes here)
//     |  o
//     o (can reach 5 nodes here)
//
//     max degree is 12
// 
//
// RETURN: the largest fanout possible without creating node violations
NUM_ELEMENTS LEVEL_NODE::get_max_degree_possible() const
{
	EDGES::const_iterator edge_iter;
	EDGE * edge = 0;
	NUM_ELEMENTS max_degree = 0;
	LEVEL_NODE * sink_level_node = 0;
	NUM_ELEMENTS edge_weight = 0;

	for (edge_iter = m_output_edges.begin(); edge_iter != m_output_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		edge_weight = edge->get_weight();

		if (edge_weight > 0)
		{
			sink_level_node = edge->get_sink_level_node();
			assert(sink_level_node);
			assert(sink_level_node->get_excess_number_of_inputs() == 0);

			max_degree += MIN(edge_weight, sink_level_node->get_nNodes());
		}
	}

	return max_degree;
}


// RETURN: the largest fanout degree assigned to this node.
NUM_ELEMENTS LEVEL_NODE::get_max_degree_assigned() const
{
	NUM_ELEMENTS max_degree = static_cast<NUM_ELEMENTS>(m_out_degrees.size() - 1);

	while (m_out_degrees[max_degree] == 0 && max_degree > 0)
	{
		max_degree--;
	}

	return max_degree;
}

// RETURN: the first indiv. node with an open input
NODE * LEVEL_NODE::get_indiv_node_with_open_input()
{
	NODES::iterator node_iter;
	NODE * node = 0;

	for (node_iter = m_nodes.begin(); node_iter != m_nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		if (node->get_nInputs() < m_kin)
		{
			return node;
		}
	}

	debugif(DPOST_PROCESS,"No indiv node with open input could be found within level node " << m_name);

	return 0;
}

// POST: nPO_to_add has been added to this level node and cluster
void LEVEL_NODE::add_PO
(
	const NUM_ELEMENTS& nPO_to_add
) 
{ 
	m_nPO += nPO_to_add; 
	assert(m_nPO >= 0);

	m_cluster->add_PO(nPO_to_add);
}

// RETURN: the number of excessive inputs on the indiv. nodes
NUM_ELEMENTS LEVEL_NODE::get_nExcess_inputs_on_indiv_nodes() const
{
	NUM_ELEMENTS total_excess_inputs = 0;
	NODES::const_iterator node_iter;
	NODE * node = 0;
	
	for (node_iter = m_nodes.begin(); node_iter != m_nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		total_excess_inputs += node->get_nExcess_inputs();
	}

	return total_excess_inputs;
}

// RETURN: the total number of edges input into the indiv. nodes within this 
//         level node that cause multiple connections from the same source nodes
NUM_ELEMENTS LEVEL_NODE::get_nDouble_inputs_on_indiv_nodes() const
{
	NUM_ELEMENTS nDouble_inputs = 0;
	NODES::const_iterator node_iter;
	NODE * node = 0;
	
	for (node_iter = m_nodes.begin(); node_iter != m_nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		nDouble_inputs += node->get_double_input_cost();
	}

	return nDouble_inputs;
}

// RETURN: the number of dff that loop back to themselves
//
NUM_ELEMENTS LEVEL_NODE::get_nDFF_loops() const
{
	assert(m_delay_level == 0);
	NODES::const_iterator node_iter;
	NODE * node = 0;
	NUM_ELEMENTS nDFF_loops = 0;

	if (m_nDFF == 0) 
	{
		return 0;
	}

	for (node_iter = m_nodes.begin(); node_iter != m_nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);

		if (node->is_DFF() && node->is_dff_have_loop_back_to_itself())
		{
			nDFF_loops++;
		}
	}

	return nDFF_loops;
}






// RETURN: true if we have only one dff super edge into this level node and
//         if the source of the super edge is ourselves
bool LEVEL_NODE::only_dff_super_edge_is_this_level_node() const
{
	assert(m_input_to_dff_edges.size() >= 1);

	// there is more than one super edge
	if (m_input_to_dff_edges.size() > 1) 
	{
		return false;
	}


	EDGES::const_iterator edge_iter = m_input_to_dff_edges.begin();

	EDGE * input_edge = *edge_iter;
	assert(input_edge);

	LEVEL_NODE * source_node = input_edge->get_source_level_node();
	assert(source_node);

	return (source_node == this);
}
void LEVEL_NODE::remove_node
(
	NODE * node_to_remove
)
{
	NODES::iterator node_iter;

	assert(node_to_remove);

	if (node_to_remove->is_DFF())
	{
		m_nDFF--;
	}
	if (node_to_remove->is_PO())
	{
		m_nPO--;
	}
	if (node_to_remove->is_latched())
	{
		m_nlatched--;
	}
	m_nNodes--;

	node_iter = find(m_nodes.begin(), m_nodes.end(), node_to_remove);

	assert(node_iter != m_nodes.end());

	m_nodes.erase(node_iter);
}
