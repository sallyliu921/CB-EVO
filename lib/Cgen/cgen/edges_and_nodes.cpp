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



#include "edges_and_nodes.h"
#include "gen.h"
#include "util.h"
#include <algorithm>

NODE::NODE
(
	const CLUSTER_NUMBER_TYPE& cluster_number,
	const DELAY_TYPE& delay_level,
	const ID_TYPE& id,
	const NUM_ELEMENTS& kin,
	LEVEL_NODE * my_level_node
)
{
	m_cluster_number = cluster_number;
	m_delay_level 	= delay_level;
	m_id 	= id;
	m_kin 	= kin;

	m_index 	= 0;
	m_span 	= 0;
	m_lspan 	= 0;
	m_rspan 	= 0;

	m_fanin_degree	= 0;
	m_fanout_degree	= 0;
	m_edges_to_assign=0;
	m_edges_currently_assigning = 0;
	m_lost_edges	 = 0;

	m_is_PI	= false;
	m_is_PO	= false;
	m_is_dff= false;
	m_is_latched = false;
	m_mark	= false;
	m_is_latch_taken = false;

	m_colour = NODE::NONE;

	m_unassigned_long_and_inter_cluster_edges = 0;

	m_my_level_node = my_level_node;
}

NODE::~NODE()
{

}

void NODE::set_fanout_degree
(
	const NUM_ELEMENTS& fanout_degree
) 
{ 
	// we haven't set the fanout before
	assert(m_fanout_degree == 0);
	m_fanout_degree = fanout_degree; 
}


// POST: adds the edges that this node will assign.
//       (both intra-cluster and inter-cluster)
//       
void NODE::add_edges_to_assign
(
	const NUM_ELEMENTS& edges_to_assign
) 
{ 
	m_edges_to_assign += edges_to_assign; 
	assert(m_edges_to_assign >= 0);
}

//
// RETURN: the number of fanouts to combinational nodes
//
NUM_ELEMENTS NODE::get_nComb_Outputs() const
{
	EDGES::const_iterator edge_iter;
	NUM_ELEMENTS nComb_Outputs = 0;
	NODE * sink_node = 0;
	EDGE * edge = 0;

	for (edge_iter = m_output_edges.begin(); edge_iter != m_output_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		sink_node = edge->get_sink_node();
		assert(sink_node);

		if (! sink_node->is_DFF())
		{
			nComb_Outputs++;
		}
	}

	return nComb_Outputs;
}


//
// RETURN: the total wirelength-approx of the connections to this node
//         (sum over all inputs to this node of the difference in horizontal positions)
//
COST_TYPE NODE::get_wirelength_cost() const
{
	COST_TYPE total_wirelength = 0,
	          wirelength = 0;

	EDGES::const_iterator edge_iter;
	EDGE * edge = 0;
	NODE * source_node = 0;
	
	for (edge_iter = m_input_edges.begin(); edge_iter != m_input_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		source_node = edge->get_source_node();
		assert(source_node);

		wirelength = ABS(source_node->get_index() - m_index);

		total_wirelength += wirelength;
	}

	return total_wirelength;
}

// the number of input edges still needed
// 
// Number of inputs needed:
// 	Combinational node at least 2 but we will settle for 1.
//  Flip-flops: 1
//  Primary inputs: 0
//
// RETURN: the number of additional indiv. edges that this node needs
NUM_ELEMENTS NODE::get_nInputs_still_needed() const
{
	NUM_ELEMENTS nInputs_needed = 0;
	NUM_ELEMENTS nInputs = get_nInputs();

	if (! is_DFF() && ! is_PI())
	{
		nInputs_needed = MAX(1-nInputs, 0); //+ 0.25*(MAX(2-nInputs, 0));
	}
	else if (! is_PI())
	{
		assert(is_DFF());
		nInputs_needed = MAX(1-nInputs, 0);
	}

	//debug("Node " << get_name() << " has " << nInputs << " inputs and needs " << nInputs_needed);

	return nInputs_needed;
}

// PRE: node is a flip-flop
// RETURN: the input edge of the flip-flop
EDGE* NODE::get_dff_input_edge() const 
{ 
	assert(is_DFF());
	assert(m_input_edges.size() == 1);
	return *(m_input_edges.begin()); 
}
// RETURN: get the source node of the DFF
NODE*  NODE::get_source_node_of_DFF() const
{
	EDGE * input_edge = get_dff_input_edge();

	assert(input_edge);

	NODE * source_node = input_edge->get_source_node();

	return source_node;
}

// add unassigned long and inter_cluster edges to this node
void NODE::add_unassigned_long_and_inter_cluster_edges(const NUM_ELEMENTS& nEdges)
{ 
	m_unassigned_long_and_inter_cluster_edges += nEdges; 
	assert(m_unassigned_long_and_inter_cluster_edges >= 0);
}

// add edges that we will assign shortly
void NODE::add_edges_currently_assigning(const NUM_ELEMENTS& edges_to_assign)
{
	m_edges_currently_assigning += edges_to_assign;
	assert(m_edges_currently_assigning >= 0);
}

// add a simple input edge to this node
void NODE::add_input_edge(EDGE * input_edge)
{
	assert(input_edge);

	if (m_is_dff)
	{
		assert(get_nInputs() == 0);
	}

	m_input_edges.push_back(input_edge);
}
// add a simple output edge to this node
void NODE::add_output_edge(EDGE * output_edge)
{
	assert(output_edge);

	m_output_edges.push_back(output_edge);
}

// removes a simple edge from this node
void NODE::remove_input_edge(EDGE * input_edge)
{
	assert(input_edge);
	//assert(is_combinational());

	EDGES::iterator edge_iter = find(m_input_edges.begin(), m_input_edges.end(), input_edge);

	assert(edge_iter != m_input_edges.end() && *edge_iter == input_edge);

	m_input_edges.erase(edge_iter);
}

// removes an simple output edge from this node
void NODE::remove_output_edge(EDGE * output_edge)
{
	assert(output_edge);

	EDGES::iterator edge_iter = find(m_output_edges.begin(), m_output_edges.end(), output_edge);

	assert(edge_iter != m_output_edges.end() && *edge_iter == output_edge);

	m_output_edges.erase(edge_iter);
}

// RETURN: a node for this name 
string NODE::get_name() const
{
	string name;
	NUM_ELEMENTS id = static_cast<NUM_ELEMENTS>(m_id);

	if (m_is_PI && m_is_dff)
	{
		name = "pi_ff";
	}
	else if (m_is_PI)
	{
		assert(! m_is_dff);
		name = "pi";
	}
	else if (m_is_dff)
	{
		assert(! m_is_PI);
		name = "ff";
	}
	else 
	{
		name = "le";
	}

	name += util_long_to_string(id) + "_"; 
	name += m_my_level_node->get_name();


	return name;
}
string NODE::get_input_node_names() const
{
	string input_node_names;
	string input_name;
	NODE * source_node = 0;
	EDGES::const_iterator edge_iter;
	EDGE * edge = 0;

	//debug(" for node " << get_name() << " there are " << m_input_edges.size() << " inputs. ");
	for (edge_iter = m_input_edges.begin(); edge_iter != m_input_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		source_node = edge->get_source_node();
		assert(source_node);

		input_name = source_node->get_name();

		input_node_names += input_name + " ";
	}

	return input_node_names;
}

// RETURN: true if the node has more than one connection from the same source node
bool NODE::has_double_inputs() const
{
	NODE * source_node = 0;
	EDGES::const_iterator edge_iter;
	EDGE * edge = 0;
	bool double_inputs = false;
	NODES source_nodes;

	for (edge_iter = m_input_edges.begin(); ! double_inputs && edge_iter != m_input_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		source_node = edge->get_source_node();
		assert(source_node);

		if (find(source_nodes.begin(), source_nodes.end(), source_node) != source_nodes.end())
		{
			double_inputs = true;
		}
		else
		{
			source_nodes.push_back(source_node);
		}
	}

	return double_inputs;
}

// RETURN: true if the node has more than kin inputs
bool NODE::has_excess_inputs() const
{
	NUM_ELEMENTS nInputs = get_nInputs();

	bool excess_inputs = (nInputs > m_kin);

	return excess_inputs;
}
// RETURN: true if the node is:
// 		   - a dff and has < 1 input
// 		   - a combinational node < 2 inputs
//
// 		   false otherwise
bool NODE::needs_inputs() const
{
	NUM_ELEMENTS nInputs = get_nInputs();

	bool needs_more_inputs = false;

	if (m_is_dff)
	{
		needs_more_inputs = (nInputs < 1);
	}
	else if (m_is_PI)
	{
		// needs no inputs
	}
	else
	{
		needs_more_inputs = (nInputs < 2);
	}


	return needs_more_inputs;
}

// the node is effectively a buffer if:
// 1. its a buffer
// 2. all the input edges come from the same source node
//
// RETURN: true if effectively a buffer node
bool NODE::is_effectively_buffer() const
{
	// if the node is a buffer... it's a buffer
	if (is_buffer())
	{
		return true;
	}

	//if the node is not a buffer if it is a dff or if it has no input edges
	if (! is_combinational() || m_input_edges.empty())
	{
		return false;
	}

	assert(m_input_edges.size() > 1);

	EDGES::const_iterator edge_iter = m_input_edges.begin();
	EDGE * first_edge = *edge_iter;
	EDGE * edge = 0;
	bool all_edges_are_the_same = true;

	assert(first_edge);
	edge_iter++;

	while (all_edges_are_the_same && edge_iter != m_input_edges.end())
	{
		edge = *edge_iter;
		assert(edge);

		all_edges_are_the_same = (edge == first_edge);
		edge_iter++;
	}
	assert(!all_edges_are_the_same || edge_iter != m_input_edges.end());

	return all_edges_are_the_same;
}


// if the node is a dff and connects to itself we have a loop
// if the node is a dff and connects to a buffer that connects back to itself we have a loop
//
// RETURN: true if the node is a dff with a loop
bool NODE::is_dff_have_loop_back_to_itself() const
{
	if (! m_is_dff)
	{
		return false;
	}

	NODE * source_node = get_source_node_of_DFF();
	assert(source_node);

	if (source_node == this)
	{
		// we have a loop
		return true;
	}

	// assert from this point forward we do not have a dff looping back to itself


	assert(new_gen);
	// test dffA->buffer->dffA
	// note: the buffer might be a buffer with two inputs from the same source node
	//
	// note2: we don't test for dff->dff->dffA
	//
	// ie.
	//  ____________ 
	//  |  _        |
	// dff/ \buffer-|
	//    \_/
	if (source_node->is_buffer() || source_node->is_effectively_buffer())
	{
		EDGES input_edge_of_buffer = source_node->get_input_edges();
		assert(! input_edge_of_buffer.empty());

		EDGE * input_edge = input_edge_of_buffer[0];
		assert(input_edge);

		NODE * source_node_of_buffer = input_edge->get_source_node();
		assert(source_node_of_buffer);

		if (source_node_of_buffer == this)
		{
			//Warning("Node " << get_name() << " has a dff loop with inputs " << get_input_node_names());
				
			// we have a loop
			return true;
		}
	} 

	return false;
}




// RETURN: the number of inputs above kin
NUM_ELEMENTS NODE::get_nExcess_inputs() const
{
	NUM_ELEMENTS nExcess_inputs = 0;

	nExcess_inputs = get_nInputs() - m_kin;
	nExcess_inputs = MAX(nExcess_inputs, 0);

	return nExcess_inputs;
}

// RETURNS: the number of times of a source node connects to this node above 
//          the one legal connections
//
// ie.
//
// a source node connection once to this node contributes 0
// a source node connection twice to this node contributes 1
// a source node connection 3 times to this node contributes 2
//
NUM_ELEMENTS NODE::get_nDouble_inputs() const
{
	NODE * source_node = 0;
	EDGES::const_iterator edge_iter;
	EDGE * edge = 0;
	NODES source_nodes;
	NUM_ELEMENTS nDouble_inputs = 0;


	// we count the number of double inputs funny.
	// not sure what is the right way to count double inputs
	// do we count the number of times a source node is repeated (present) or
	// do we count the number of double edges found?
	// these are two different methods
	assert(new_gen);

	for (edge_iter = m_input_edges.begin(); edge_iter != m_input_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		source_node = edge->get_source_node();
		assert(source_node);

		if (find(source_nodes.begin(), source_nodes.end(), source_node) != source_nodes.end())
		{
			// every time we come across this source node increase the nDouble_inputs
			// note: we don't have to add it to the list because we only need
			// to know that the source node already exists.
			nDouble_inputs++;
		}
		else
		{
			source_nodes.push_back(source_node);
		}
	}

	return nDouble_inputs;
}
NUM_ELEMENTS NODE::get_double_input_cost() const
{
	NUM_ELEMENTS cost = get_nDouble_inputs();
	
	if (cost > 1)
	{
		cost *= 10;
	}

	return cost;
}

// RETURN: the edges that contribute to this node having multiple connections from the same
//         source node(s)
EDGES NODE::get_double_input_edges() const
{
	NODE * source_node = 0;
	EDGES::const_iterator edge_iter;
	EDGE * edge = 0;
	NODES source_nodes,
		  double_source_nodes;
	EDGES double_input_edges;

	// 1st build the list of source nodes that feed the same input node more than once
	for (edge_iter = m_input_edges.begin(); edge_iter != m_input_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		source_node = edge->get_source_node();
		assert(source_node);

		if (find(source_nodes.begin(), source_nodes.end(), source_node) != source_nodes.end())
		{
			// we have a source node with multiple edges to the same input node
			// see if we have seen this source node before.
			// if not add it to the double source nodes list
			if (find(double_source_nodes.begin(), double_source_nodes.end(), source_node) 
					== double_source_nodes.end())
			{
				double_source_nodes.push_back(source_node);
			}
		}
		else
		{
			source_nodes.push_back(source_node);
		}
	}

	// 2nd add all edges that have a double input source node to the list of double input edges
	for (edge_iter = m_input_edges.begin(); edge_iter != m_input_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		source_node = edge->get_source_node();
		assert(source_node);

		if (find(double_source_nodes.begin(), double_source_nodes.end(), source_node) 
				!= double_source_nodes.end())
		{
			double_input_edges.push_back(edge);
		}
	}

	return double_input_edges;
}


// RETURN: the first edge found that ouputs to a dff
EDGE* NODE::get_output_edge_that_connects_to_DFF() const
{
	EDGES::const_iterator edge_iter;
	EDGE * output_edge = 0;
	EDGE * dff_edge = 0;
	NODE * sink_node = 0;
	bool edge_not_found = true;

	assert(m_is_latched);

	for (edge_iter = m_output_edges.begin(); edge_not_found && edge_iter != m_output_edges.end(); edge_iter++)
	{
		output_edge = *edge_iter;
		assert(output_edge);

		sink_node = output_edge->get_sink_node();
		assert(sink_node);

		if (sink_node->is_DFF())
		{
			edge_not_found = false;
			dff_edge = output_edge;
		}
	}

	return dff_edge;
}
EDGE::EDGE()
{
	m_type = NODE_TO_NODE;

	m_weight= 0;
	m_max_weight = 0;
	m_length= -1;

	m_source_node= 0;
	m_sink_node	= 0;
	m_source_level_node	= 0;
	m_sink_level_node	= 0;
	m_my_super_edge 	= 0;


}
EDGE::EDGE
(
	LEVEL_NODE* source_node,
	LEVEL_NODE* sink_node,
	const ID_TYPE& id
)
{
	assert(source_node && sink_node);

	m_type = EDGE::LEVEL_NODE_TO_LEVEL_NODE;
	m_weight= 0;
	m_max_weight = 0;

	m_source_node= 0;
	m_sink_node	= 0;
	m_source_level_node	= source_node;
	m_sink_level_node	= sink_node;
	m_my_super_edge 	= 0;
	m_id				= id;

	if (m_sink_level_node->get_delay_level() > 0)
	{
		m_length= m_sink_level_node->get_delay_level() - m_source_level_node->get_delay_level();
	}
	else
	{
		// we have a dff edge.  length is zero
		assert(m_sink_level_node->get_delay_level() == 0);
		m_length = 0;
	}
}
EDGE::EDGE
(
	NODE* source_node, 
	NODE* sink_node
)
{
	assert(source_node && sink_node);

	m_type = EDGE::NODE_TO_NODE;
	m_weight= 1;
	m_max_weight = 1;

	m_source_node= source_node;
	m_sink_node	= sink_node;
	m_source_level_node	= 0;
	m_sink_level_node	= 0;

	if (m_sink_node->get_delay_level() > 0)
	{
		m_length = m_sink_node->get_delay_level() - m_source_node->get_delay_level();
	}
	else
	{
		// we have a dff edge.  length is zero
		assert(m_sink_node->get_delay_level() == 0);
		m_length = 0;
	}
	m_my_super_edge 	= 0;
}

EDGE::~EDGE()
{
	
}


// RETURN: true if the edge is intra-cluster
bool EDGE::is_intra_cluster() const
{
	CLUSTER_NUMBER_TYPE sink_cluster_number = 0;
	CLUSTER_NUMBER_TYPE source_cluster_number = 0;

	if (m_type == EDGE::NODE_TO_NODE)
	{
		assert(m_source_node && m_sink_node);
		sink_cluster_number = m_sink_node->get_cluster_number();
		source_cluster_number = m_source_node->get_cluster_number();
	}
	else
	{
		assert(m_sink_level_node && m_source_level_node);
		sink_cluster_number = m_sink_level_node->get_cluster_number();
		source_cluster_number = m_source_level_node->get_cluster_number();
	}

	return (sink_cluster_number == source_cluster_number);
}

// RETURN: true if the edge is inter-cluster
bool EDGE::is_inter_cluster() const
{
	return (! is_intra_cluster());
}


// PRE: input_cluster_number is the sink cluster number
// RETURN: true if the edge inputs into the cluster in question
bool EDGE::is_inter_cluster_input
(
	const CLUSTER_NUMBER_TYPE & input_cluster_number
) const
{
	CLUSTER_NUMBER_TYPE sink_cluster_number;
	bool is_inter_cluster_input = false; 

	if (m_type == EDGE::NODE_TO_NODE)
	{
		assert(m_sink_node);
		sink_cluster_number = m_sink_node->get_cluster_number();
	}
	else
	{
		assert(m_sink_level_node);
		sink_cluster_number = m_sink_level_node->get_cluster_number();
	}

	is_inter_cluster_input = is_inter_cluster() && (sink_cluster_number == input_cluster_number);

	return is_inter_cluster_input;
}
// PRE: output_cluster_number is the source cluster number
// RETURN: true if the edge outputs out of the cluster in question
bool EDGE::is_inter_cluster_output
(
	const CLUSTER_NUMBER_TYPE & output_cluster_number
) const
{
	CLUSTER_NUMBER_TYPE source_cluster_number;
	bool is_inter_cluster_output = false; 

	if (m_type == EDGE::NODE_TO_NODE)
	{
		assert(m_source_node);
		source_cluster_number = m_source_node->get_cluster_number();
	}
	else
	{
		assert(m_source_level_node);
		source_cluster_number = m_source_level_node->get_cluster_number();
	}

	is_inter_cluster_output = is_inter_cluster() && (source_cluster_number == output_cluster_number);

	return is_inter_cluster_output;
}

// RETURN: true if the sink node of the edge is a flip-flop
bool EDGE::is_sink_a_dff() const
{
	bool is_dff = false; 

	if (m_type == EDGE::NODE_TO_NODE)
	{
		assert(m_sink_node);
		is_dff = m_sink_node->is_DFF();
	}
	else
	{
		assert(m_sink_level_node);
		is_dff = (m_sink_level_node->get_delay_level() == 0);
	}

	return is_dff;
}
// RETURN: the cluster number of the sink node of the edge 
CLUSTER_NUMBER_TYPE EDGE::get_sink_cluster_number() const
{
	if (m_type == NODE_TO_NODE)
	{
		return (m_sink_node->get_cluster_number());
	}
	else
	{
		return (m_sink_level_node->get_cluster_number());
	}

}
// RETURN: the cluster number of the source node of the edge 
CLUSTER_NUMBER_TYPE EDGE::get_source_cluster_number() const
{
	if (m_type == NODE_TO_NODE)
	{
		return (m_source_node->get_cluster_number());
	}
	else
	{
		return (m_source_level_node->get_cluster_number());
	}

}


// POST: the weight of the edge and its level nodes have been updated
void EDGE::set_weight_and_update_level_nodes(const NUM_ELEMENTS& weight)
{
	assert(m_type == EDGE::LEVEL_NODE_TO_LEVEL_NODE);
	assert(m_source_level_node && m_sink_level_node);
	NUM_ELEMENTS weight_addition = weight - m_weight;

	m_weight = weight; 

	assert(m_weight >= 0);
	assert(m_weight <= m_max_weight);

	if (weight_addition != 0)
	{
		m_source_level_node->output_edge_has_changed(this, weight_addition);
		m_sink_level_node->input_edge_has_changed(this, weight_addition);
	}
}
// POST: the weight of the edge has been set
void EDGE::set_weight(const NUM_ELEMENTS& weight)
{
	assert(m_type == EDGE::LEVEL_NODE_TO_LEVEL_NODE);
	assert(m_source_level_node && m_sink_level_node);

	m_weight = weight; 

	assert(m_weight >= 0);
	assert(m_weight <= m_max_weight);
}

// sets the weight 
// don't check to see if we are over the maximum weight
void EDGE::unsafe_set_weight(const NUM_ELEMENTS& weight)
{
	assert(m_type == EDGE::LEVEL_NODE_TO_LEVEL_NODE);
	assert(m_source_level_node && m_sink_level_node);

	m_weight = weight; 

	assert(m_weight >= 0);
}



// PRE: the weight change won't make the edge weight negative
// POST: the weight has been added to the edge
//        the level nodes have been updated
void EDGE::add_weight(const NUM_ELEMENTS& weight_change)
{
	assert(m_type == EDGE::LEVEL_NODE_TO_LEVEL_NODE);
	assert(m_source_level_node && m_sink_level_node);

	m_weight += weight_change;
	assert(m_weight >= 0);
	assert(m_weight <= m_max_weight);

	if (weight_change != 0)
	{
		m_source_level_node->output_edge_has_changed(this, weight_change);
		m_sink_level_node->input_edge_has_changed(this, weight_change);
		if (m_length == 1)
		{
			assert(m_sink_level_node->get_delay_difference() >= 0);
		}
	}
}


// the selection weight of an edge is g
//
// The selection weight of an edge is higher for edges that contribute more 
// to the costLevel_Shape.
// We give higher selection weights to super edges between source level nodes 
// that have too many outputs and too many inputs and give lower selection weights 
// to super edges that have source level nodes with not enough outputs or sink level 
// nodes with not enough inputs. We also add a small bias to the selection weight of 
// each super edge to ensure that each edge has a small probability of being chosen. 
//
// RETURN: the selection weight
COST_TYPE EDGE::get_selection_weight() const
{

	assert(m_source_level_node && m_sink_level_node);
	COST_TYPE cost_contribution = 0;

	cost_contribution = static_cast<double>(m_source_level_node->get_output_difference()) *
			static_cast<double>(m_weight)/static_cast<double>(m_source_level_node->get_output_edge_weight());

	cost_contribution += static_cast<double>(m_sink_level_node->get_input_difference()) *
				static_cast<double>(m_weight)/static_cast<double>(m_sink_level_node->get_input_edge_weight());

	return cost_contribution;
}

// RETURN: true if the super edge has weight that it can give up 
//         which is contrained by not going lower than the input and output shapes
//         and not violating the ability of the sink level node to define the 
//         delay level of its indiv. nodes
bool EDGE::has_spare_weight() const
{
	assert(m_source_level_node && m_sink_level_node);


	bool spare_weight = (m_weight > 0) && 
						(m_source_level_node->get_output_to_spare() > 0) &&
						(m_sink_level_node->get_input_to_spare() > 0);

	if (m_length == 1)
	{
		spare_weight = spare_weight && (m_sink_level_node->get_delay_to_spare() > 0);
	}

	return spare_weight;
}

// RETURN: true if the super edge has weight that can be reduced
//         which is contrained by the weight of the edge and by
//         not violating the ability of the sink level node to define the 
//         delay level of its indiv. nodes
bool EDGE::can_reduce_weight_on_edge() const
{
	assert(m_source_level_node && m_sink_level_node);

	bool spare_weight = (m_weight > 0);

	if (m_length == 1)
	{
		spare_weight = spare_weight && (m_sink_level_node->get_delay_difference() > 0);
	}

	return spare_weight;
}

// RETURN: the weight that the super edge can be reduced
//         (the number of simple edges that can removed)
//         which is contrained by:
//         - the weight of the edge
//         - making sure the source level node has enough outputs
//         - making sure the sink level node can define the 
//           delay level of its indiv. nodes
NUM_ELEMENTS EDGE::get_nUnneeded_edges() const
{
	assert(m_source_level_node && m_sink_level_node);
	assert(m_weight >= 0);

	NUM_ELEMENTS spare_weight = m_weight;
	NUM_ELEMENTS min_output_weight = m_source_level_node->get_output_edge_weight() - 
	(m_source_level_node->get_nNodes() - m_source_level_node->get_nLatched() - m_source_level_node->get_nPO());
	//assert(min_output_weight >= 0);
	min_output_weight = MAX(min_output_weight, 0);

	spare_weight = MIN(spare_weight, min_output_weight);

	if (m_length == 1)
	{
		spare_weight = MIN(spare_weight, m_sink_level_node->get_delay_to_spare());
	}

	assert(spare_weight >= 0);

	return spare_weight;
}
// RETURN: the edge weight that can be reduced on this super edge 
//         which is contrained by the weight of the edge and by
//         not violating the ability of the sink level node to define the 
//         delay level of its indiv. nodes
NUM_ELEMENTS EDGE::get_nEdge_to_sacrifice() const
{
	assert(m_source_level_node && m_sink_level_node);
	assert(m_weight >= 0);

	NUM_ELEMENTS spare_weight = m_weight;

	if (m_length == 1)
	{
		spare_weight = MIN(spare_weight, m_sink_level_node->get_delay_to_spare());
	}

	assert(spare_weight >= 0);

	return spare_weight;
}

// see if this node is connected to us.
// (linear search is ok. because max length is kin)
//
// RETURN: true if the possible_source_node is connected by a simple edge to this node
bool NODE::are_nodes_connected(const NODE * possible_source_node) const
{
	EDGES::const_iterator edge_iter;
	EDGE * edge = 0;
	NODE * source_node = 0;
	
	for (edge_iter = m_input_edges.begin(); edge_iter != m_input_edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);
	
		source_node = edge->get_source_node();	

		if (source_node == possible_source_node)
		{
			return true;
		}
	}

	return false;
}


// RETURN: a string to identify the edge
string EDGE::get_info() const
{
	string info;
	CLUSTER_NUMBER_TYPE source_cluster = get_source_cluster_number();
	CLUSTER_NUMBER_TYPE sink_cluster = get_sink_cluster_number();

	info += "From cluster " + util_long_to_string(source_cluster) + 
		" to cluster " + util_long_to_string(sink_cluster) + 
		" with length " + util_long_to_string(m_length);

	return info;
}
// PRE: the edge is a super edge
// POST: a simple edge has been added to this super edge
void EDGE::add_simple_edge
(
	EDGE * simple_edge
)
{
	assert(simple_edge);
	assert(m_type == EDGE::LEVEL_NODE_TO_LEVEL_NODE);
	simple_edge->set_my_super_edge(this);

	m_simple_edges.push_back(simple_edge);
}

// PRE: the edge is a super edge
// RETURN: the number of real simple edges aggregated under this super edge
NUM_ELEMENTS EDGE::get_total_simple_edge_weight() const
{
	NUM_ELEMENTS total_edge_weight = static_cast<NUM_ELEMENTS>(m_simple_edges.size());

	return total_edge_weight;
}
void EDGE::set_sink_node(NODE * node)
{
	assert(node);
	assert(m_type == EDGE::NODE_TO_NODE);

	m_sink_node = node;
}

// POST: we are sure that:
//  	 the edge length is correct
//
void EDGE::check_sanity() const
{
	// if we don't have a dff for a sink node
	if (m_sink_node->get_delay_level() > m_source_node->get_delay_level())
	{
		assert(m_length == (m_sink_node->get_delay_level() - m_source_node->get_delay_level()));
	}
	else
	{
		assert(m_length == 0);
	}
}
// PRE: we are a simple edge
// RETURN: we now belong to a super edge
void EDGE::set_my_super_edge(EDGE * super_edge)
{
	assert(super_edge);
	m_my_super_edge = super_edge;
}

// PRE: we are a super edge that connects to the sink of a level node with flip-flops
// POST: we have updated the weight of this edge and notified the sink level node
void EDGE::add_weight_to_dff_super_edge
(
	const NUM_ELEMENTS& weight_change
)
{
	assert(m_type == EDGE::LEVEL_NODE_TO_LEVEL_NODE);
	assert(m_source_level_node && m_sink_level_node);

	m_weight += weight_change;
	assert(m_weight >= 0);
	m_max_weight = MAX(0, m_weight);

	if (weight_change != 0)
	{
		m_sink_level_node->input_edge_has_changed(this, weight_change);
	}
}

string EDGE::get_vhdl_signal_name() const 
{
	assert(m_source_node);

	string sig_name;

	if (m_source_node->is_PO())
	{  
		sig_name = "tmp_" + m_source_node->get_name();

	}
	else
	{
		sig_name = m_source_node->get_name();
	}

	return sig_name;
}
