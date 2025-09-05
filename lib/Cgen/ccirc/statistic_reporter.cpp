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





#include "statistic_reporter.h"
#include <algorithm>
#include <numeric>
#include <iterator>
#include "degree_info.h"
#include "rnum.h"

STATISTIC_REPORTER::STATISTIC_REPORTER()
{
	m_circuit	= 0;
}

STATISTIC_REPORTER::STATISTIC_REPORTER(const STATISTIC_REPORTER & another_statistic_reporter)
{
	m_circuit		= another_statistic_reporter.m_circuit;
}

STATISTIC_REPORTER & STATISTIC_REPORTER::operator=(const STATISTIC_REPORTER & another_statistic_reporter)
{
	m_circuit	= another_statistic_reporter.m_circuit;

	return (*this);
}

STATISTIC_REPORTER::~STATISTIC_REPORTER()
{
}

void STATISTIC_REPORTER::report_stats
(
	CIRCUIT * circuit
)
{
	assert(circuit);
	m_circuit = circuit;
	DEGREE_INFO * degree_info = m_circuit->get_degree_info();
	assert(degree_info);
	SEQUENTIAL_LEVEL * sequential_level = m_circuit->get_sequential_level();
	assert(sequential_level);

	string circuit_name;
	string file_name;

	file_name = g_options->get_output_file_name();

	if (file_name.empty())
	{
		circuit_name = m_circuit->get_name();
		assert(! circuit_name.empty());
		file_name = circuit_name + ".stats";
	}

	Log("About to open the statistical results file: " << file_name);

	m_output_file.open(file_name.c_str(), ios::out);

	if (! m_output_file.is_open())
	{
		Warning("Could not open output file " << file_name << 
				". Therefore, could not output a .stats file\n");
		return;
	}

	m_output_file << "######################## BASIC ############################" << endl;
	m_output_file << "Circuit_Name:  " 	<< m_circuit->get_name()	<< endl;
    m_output_file << "Number_of_Nodes:  " 	<< m_circuit->get_size()	<< endl;
    m_output_file << "Number_of_Edges: " << m_circuit->get_nEdges_without_clock_edges() << endl;
    m_output_file << "Maximum_Delay: " 		<< m_circuit->get_maximum_combinational_delay() << endl;
    m_output_file << "Number_of_PI: " 	<< m_circuit->get_nPI()		<< endl;
    m_output_file << "Number_of_PO: " 	<< m_circuit->get_nPO()		<< endl;
    m_output_file << "Number_of_Combinational_Nodes: " 	<< m_circuit->get_nComb()	<< endl;
    m_output_file << "Number_of_DFF: " 	<< m_circuit->get_nDFF()	<< endl;
    m_output_file << "kin: " << g_options->get_k()	<< endl;
	
	// if we didn't calculate wirelength approx. don't print 0 but print not_calculated
	if (g_options->is_determine_wirelength_approx())
	{
		m_output_file << "Wirelength_approx: " << m_circuit->get_wirelength_approx()	<< endl;
	}

	if (m_circuit->get_global_clock()) 
	{ 
		m_output_file << "clock: " 	<< m_circuit->get_global_clock()->get_name() << endl;
	}
	
    //m_output_file << "Num_unusable: " 	<< m_circuit->num_unusable()	<< endl;
    //m_output_file << "Num_unreachable: " << m_circuit->num_unreachable()<< endl;


	//report_by_cluster_statistics();

	report_degree_information(degree_info);

	report_level_shape(sequential_level, degree_info);

	//report_cluster_stastistics();

	report_reconvergence(circuit);
}

void STATISTIC_REPORTER::report_reconvergence(CIRCUIT * circuit){

	m_output_file << "======================== SHAPE ============================" << endl;

	//report reconvergence value 
	double R0,R0min,R0max;
	rnum (circuit,&R0,&R0max,&R0min);
	m_output_file << "Reconvergence: " << R0 << endl;
	m_output_file << "Reconvergence_max: " << R0max	<< endl;
	m_output_file << "Reconvergence_min: " << R0min	<< endl;

}

void STATISTIC_REPORTER::report_by_cluster_statistics()
{
	DISTRIBUTION size, nPI, nDFF, nIntra_cluster_edges, nInter_cluster_edges, wirelength_approx;

	m_output_file << "======================== Cluster_Summary ==================\n";

	m_circuit->get_cluster_stats(size, nPI, nDFF, nIntra_cluster_edges, nInter_cluster_edges, wirelength_approx);

	m_output_file << "Number_of_nodes: "; output_distribution(size);
	m_output_file << "Number_of_pi: "; output_distribution(nPI);
	m_output_file << "Number_of_dff: "; output_distribution(nDFF);
	m_output_file << "Number_of_intra_cluster_edges: "; output_distribution(nIntra_cluster_edges);
	m_output_file << "Number_of_inter_cluster_edges: "; output_distribution(nInter_cluster_edges);
	if (g_options->is_determine_wirelength_approx())
	{
		m_output_file << "Wirelength_approx: "; output_distribution(wirelength_approx);
	}

	m_output_file << "Partitioned_scaled_cost: " << m_circuit->get_scaled_cost() << endl;
}


void STATISTIC_REPORTER::report_degree_information
(
	DEGREE_INFO * degree_info
)
{
	assert(degree_info);

	m_output_file << "======================== DEGREE ============================" << endl;

	m_output_file << "Avg_fanin_comb: " << degree_info->get_avg_fanin_for_comb() 
		<< " (" << degree_info->get_std_dev_comb_fanin() << ")" << endl;

	m_output_file << "Avg_fanout: " 	<< degree_info->get_avg_fanout()
		<< " (" << degree_info->get_std_dev_fanout() << ")" 	<< endl;

	m_output_file << "Avg_fanout_comb: "<< degree_info->get_avg_fanout_for_comb() 	
		<< " (" << degree_info->get_std_dev_comb_fanout() << ")" << endl;

	m_output_file << "Avg_fanout_pi: " 	<< degree_info->get_avg_fanout_for_pi() 
		<< " (" << degree_info->get_std_dev_pi_fanout() << ")" 	<< endl;

	m_output_file << "Avg_fanout_dff: " << degree_info->get_avg_fanout_for_dff() 
		<< " (" << degree_info->get_std_dev_dff_fanout() << ")"	<< endl;

	m_output_file << "Maximum_fanout: " << degree_info->get_maximum_fanout_degree() 	<< endl;

	m_output_file << "Number_of_high_degree_comb: " 	<< degree_info->get_high_degree_comb()<< endl;
	m_output_file << "Number_of_high_degree_pi: " 	<< degree_info->get_high_degree_pi() 	<< endl;
	m_output_file << "Number_of_high_degree_dff: " 	<< degree_info->get_high_degree_dff()	<< endl;

	m_output_file << "Number_of_10plus_degree_comb: " << degree_info->get_10plus_fanout_degree_comb() 	<< endl;
	m_output_file << "Number_of_10plus_degree_pi: " 	<< degree_info->get_10plus_fanout_degree_pi() 	<< endl;
	m_output_file << "Number_of_10plus_degree_dff: " 	<< degree_info->get_10plus_fanout_degree_dff() 	<< endl;

}

void STATISTIC_REPORTER::report_shape_information
(
	SEQUENTIAL_LEVELS & sequential_levels,
	DEGREE_INFO * degree_info
)
{
	SEQUENTIAL_LEVELS::iterator level_iter;
	SEQUENTIAL_LEVEL * seq_level;

	m_output_file << "======================== SHAPE ============================" << endl;

	for (level_iter = sequential_levels.begin(); level_iter != sequential_levels.end(); level_iter++)
	{
		seq_level = *level_iter;
		assert(seq_level);
		report_level_shape(seq_level, degree_info);
	}

}

void STATISTIC_REPORTER::report_level_shape
(
	SEQUENTIAL_LEVEL * seq_level,
	DEGREE_INFO * degree_info
)
{
	m_output_file << "======================== SHAPE ============================" << endl;

	assert(seq_level && degree_info);

	NUM_ELEMENTS max_fanout = degree_info->get_maximum_fanout_degree();
	SHAPE node_shape = seq_level->get_node_shape();
	SHAPE input_shape = seq_level->get_input_shape();
	SHAPE output_shape_by_level = seq_level->get_output_shape();
	SHAPE primary_output_shape = seq_level->get_PO_shape();
	SHAPE latched_shape	= seq_level->get_latched_shape();
	DISTRIBUTION intra_cluster_edge_length_distribution = 
					seq_level->get_intra_cluster_edge_length_distribution();
	DISTRIBUTION inter_cluster_input_edge_length_distribution = 
					seq_level->get_inter_cluster_input_edge_length_distribution();
	DISTRIBUTION inter_cluster_output_edge_length_distribution = 
					seq_level->get_inter_cluster_output_edge_length_distribution();
	DISTRIBUTION fanout_distribution = seq_level->get_fanout_distribution(max_fanout);

	m_output_file << "Node_shape: ";
	output_shape(node_shape);

	m_output_file << "Input_shape: ";
	output_shape(input_shape);

	m_output_file << "Output_shape: ";
	output_shape(output_shape_by_level);

	m_output_file << "Latched_shape: ";
	output_shape(latched_shape);

	m_output_file << "POshape: ";
	output_shape(primary_output_shape);



	if (seq_level->is_clustered())
	{
		m_output_file << "Intra_cluster_edge_length_distribution: ";
		output_distribution(intra_cluster_edge_length_distribution);

		m_output_file << "Inter_cluster_input_edge_length_distribution: ";
		output_distribution(inter_cluster_input_edge_length_distribution);

		m_output_file << "Inter_cluster_output_edge_length_distribution: ";
		output_distribution(inter_cluster_output_edge_length_distribution);
	}
	else
	{
		m_output_file << "Edge_length_distribution: ";
		output_distribution(intra_cluster_edge_length_distribution);
	}

	m_output_file << "Fanout_distribution: ";
	output_distribution(fanout_distribution);
}

void STATISTIC_REPORTER::report_cluster_stastistics()
{
	assert(m_circuit);

	CLUSTERS & clusters = m_circuit->get_clusters();
	CLUSTERS::iterator cluster_iter;
	CLUSTER * cluster = 0;

	m_output_file << "#################### Clusters ######################" << endl;


	m_output_file << "Number_of_Clusters: " << clusters.size()	<< endl;

	for (cluster_iter = clusters.begin(); cluster_iter != clusters.end(); cluster_iter++)
	{
		cluster = *cluster_iter;
		assert(cluster);

		report_cluster(cluster);
	}

	report_inter_cluster_adjacency_matrix();
}

void STATISTIC_REPORTER::report_cluster
(
	CLUSTER * cluster
)
{
	assert(cluster);
	SEQUENTIAL_LEVEL* sequential_level = cluster->get_sequential_level();
	CLUSTER_NUMBER_TYPE cluster_number = cluster->get_cluster_number();
	DEGREE_INFO * degree_info = cluster->get_degree_info();
	assert(degree_info && sequential_level);

	m_output_file << "#################### Cluster " <<  cluster_number << " ######################" << endl;

	m_output_file << "Number_of_Nodes: " << cluster->get_size()<< endl;
	m_output_file << "Number_of_Intra_cluster_edges: " << cluster->get_nIntra_cluster_edges()<< endl;
	m_output_file << "Number_of_Inter_cluster_edges: " << cluster->get_nInter_cluster_edges()<< endl;
	m_output_file << "Number_of_PI: " 	<< cluster->get_nPI()	<< endl;
	m_output_file << "Number_of_PO: " 	<< cluster->get_nPO()	<< endl;
	m_output_file << "Number_of_Comb: " 	<< cluster->get_nComb()	<< endl;
	m_output_file << "Number_of_DFF: " 	<< cluster->get_nDFF()	<< endl;
	m_output_file << "Number_of_Latched: " << cluster->get_nLatched() << endl;
	m_output_file << "Number_of_Inter_cluster_input_edges: " 
		<< cluster->get_nInter_cluster_input_edges() << endl;
	m_output_file << "Number_of_Inter_cluster_output_edges: " 
		<< cluster->get_nInter_cluster_output_edges() << endl;

	if (g_options->is_determine_wirelength_approx())
	{
		m_output_file << "Wirelength_approx: " 	<< cluster->get_wirelength_approx()	<< endl;
	}

	report_delay_defining_stats(cluster);

	report_degree_information(degree_info);

	if (g_options->is_display_pi_and_dff_distributions())
	{
		// display the values assigned to the pi and dffs

		DISTRIBUTION pi_fanout_distribution = cluster->get_PI_fanout_values();
		DISTRIBUTION dff_fanout_distribution = cluster->get_DFF_fanout_values();
		m_output_file << "pi_fanout: ";
		output_distribution(pi_fanout_distribution);
		m_output_file << "dff_fanout: ";
		output_distribution(dff_fanout_distribution);
	}

	assert(Dsingle_seq_level);

	report_level_shape(sequential_level, degree_info);
}

void STATISTIC_REPORTER::output_shape
(
	const SHAPE & shape
)
{

	m_output_file << "( ";
	copy(shape.begin(), shape.end(), ostream_iterator<NUM_ELEMENTS>(m_output_file, " "));
	m_output_file << ")" << endl;
}

void STATISTIC_REPORTER::output_distribution
(
	const DISTRIBUTION & distribution
)
{
	m_output_file << "( ";
	copy(distribution.begin(), distribution.end(), ostream_iterator<NUM_ELEMENTS>(m_output_file, " "));
	m_output_file << ")" << endl;
}


void STATISTIC_REPORTER::report_inter_cluster_adjacency_matrix()
{
	MATRIX &  inter_cluster_connections = m_circuit->get_inter_cluster_adjacency_matrix();
	MATRIX & inter_cluster_connections_for_dff = m_circuit->get_inter_cluster_adjacency_matrix_for_dff();
	INDEX_SIZE nRows, nCols;

	nRows = inter_cluster_connections.get_nRows();
	nCols = inter_cluster_connections.get_nColumns();

	m_output_file << "-------------------- Inter_cluster_adjacentcy_matrix_to_combinational_nodes --------------------\n\n";

	m_output_file << inter_cluster_connections;

	m_output_file << "\n-------------------- Inter_cluster_adjacentcy_matrix_to_dffs -----------------------------------\n\n";

	m_output_file << inter_cluster_connections_for_dff << endl;


	if (g_options->is_display_inter_cluster_matricies_at_each_edge_length())
	{
		m_output_file << "\n\nMatrices at each edge length\n";

		DELAY_TYPE max_delay = m_circuit->get_maximum_combinational_delay();
		LENGTH_TYPE edge_length = 0;

		for (edge_length = 0; edge_length <= max_delay; edge_length++)
		{

			m_output_file << "\nEdge_length: " << edge_length << endl;
			m_output_file << m_circuit->get_inter_cluster_matrix(edge_length);
		}
	}
}

void STATISTIC_REPORTER::report_delay_defining_stats
(
	CLUSTER * cluster
)
{
	assert(cluster);
	NODES & nodes = cluster->get_nodes();
	NODES::const_iterator node_iter;
	NODE * node = 0;
	NUM_ELEMENTS number_of_ties, number_of_inter_cluster_edges, number_of_internal_edges;
	number_of_ties = number_of_inter_cluster_edges = number_of_internal_edges = 0;

	for (node_iter = nodes.begin(); node_iter != nodes.end(); node_iter++)
	{
		node = *node_iter;
		assert(node);


		if (node->get_type() == NODE::COMB)
		{
			update_stats_regarding_edges_that_define_delay(node, number_of_ties, 
										number_of_inter_cluster_edges, number_of_internal_edges);
		}
	}

	assert(cluster->get_nComb() == number_of_ties + number_of_inter_cluster_edges + number_of_internal_edges);

	if (g_options->is_display_statistics_on_delay_defining_edges())
	{
		m_output_file << "---------Statistics for delay defining edges-----------" << endl;
		m_output_file << "Number_of internal edges " << number_of_internal_edges << endl;
		m_output_file << "Number of inter cluster edges " << number_of_inter_cluster_edges << endl;
		m_output_file << "Number of ties for delay defining edges " << number_of_ties << endl;
	}
}

void STATISTIC_REPORTER::update_stats_regarding_edges_that_define_delay
(
	NODE * node,
	NUM_ELEMENTS & number_of_ties,
	NUM_ELEMENTS & number_of_inter_cluster_edges,
	NUM_ELEMENTS & number_of_internal_edges
)
{
	assert(node);
	bool inter_cluster_edge_defines_delay = false;
	bool internal_edge_defines_delay = false;
	NODE * source_node = 0;
	DELAY_TYPE comb_delay = node->get_max_comb_delay_level();
	DELAY_TYPE comb_delay_of_source_node = 0;
	EDGES edges = node->get_input_edges();
	EDGES::iterator edge_iter;
	EDGE * edge = 0;

	for (edge_iter = edges.begin(); edge_iter != edges.end(); edge_iter++)
	{
		edge = *edge_iter;
		assert(edge);

		source_node = edge->get_source_node();

		if (source_node)
		{
			comb_delay_of_source_node = source_node->get_max_comb_delay_level();
		}
		else
		{
			comb_delay_of_source_node = 0;		// connects to a primary input
		}

		// check to see if this edge defines delay
		if (comb_delay_of_source_node+1 == comb_delay)
		{
			// we have a delay defining edge
			if (edge->is_inter_cluster())
			{
				inter_cluster_edge_defines_delay = true;
			}
			else
			{
				internal_edge_defines_delay = true;
			}
		}
	}

	// one of them must be true or else we are in trouble.
	assert(inter_cluster_edge_defines_delay || internal_edge_defines_delay);

	if (inter_cluster_edge_defines_delay && internal_edge_defines_delay)
	{
		number_of_ties++;
	}
	else if (inter_cluster_edge_defines_delay)
	{
		number_of_inter_cluster_edges++;
	}
	else
	{
		number_of_internal_edges++;
	}
}

void STATISTIC_REPORTER::report_edge_length_by_delay_level
(
	SEQUENTIAL_LEVEL * seq_level
)
{

	assert(seq_level);
	assert(m_circuit);
	DELAY_TYPE delay = 0;
	DISTRIBUTION edge_lengths;
	m_output_file << "edge_lengths_by_delay_level: " << "\n\n";

	for (delay = 0; delay <= m_circuit->get_maximum_combinational_delay(); delay++)
	{
		m_output_file << "Delay " << delay << "\n\n";
	
		m_output_file << "Intra_cluster_input_edge_lengths\n";
		edge_lengths = seq_level->get_intra_cluster_input_edge_length_distribution(delay);
		copy(edge_lengths.begin(), edge_lengths.end(), ostream_iterator<NUM_ELEMENTS>(m_output_file, " "));

		m_output_file << "\nIntra_cluster_output_edge_lengths\n";
		edge_lengths = seq_level->get_intra_cluster_output_edge_length_distribution(delay);
		copy(edge_lengths.begin(), edge_lengths.end(), ostream_iterator<NUM_ELEMENTS>(m_output_file, " "));

		m_output_file << "\nInter_cluster_input_edge_lengths\n";
		edge_lengths = seq_level->get_inter_cluster_input_edge_length_distribution(delay);	
		copy(edge_lengths.begin(), edge_lengths.end(), ostream_iterator<NUM_ELEMENTS>(m_output_file, " "));
		
		m_output_file << "\nInter_cluster_output_edge_lengths\n";
		edge_lengths = seq_level->get_inter_cluster_output_edge_length_distribution(delay);
		copy(edge_lengths.begin(), edge_lengths.end(), ostream_iterator<NUM_ELEMENTS>(m_output_file, " "));

		m_output_file << endl << endl;
	}
}
void STATISTIC_REPORTER::report_cluster_size_distribution
(
	CLUSTER * cluster
)
{
	assert(cluster);

	CLUSTER_NUMBER_TYPE max_cluster_depth = cluster->get_max_clustering_depth();
	CLUSTER_NUMBER_TYPE depth = 0;
	DISTRIBUTION size_distribution;
	SHAPE node_shape;

	for (depth = 0; depth <= max_cluster_depth; depth++)
	{
		m_output_file << "\ndepth " << depth << "\n";
		size_distribution = cluster->get_size_distribution(depth);
		copy(size_distribution.begin(), size_distribution.end(), ostream_iterator<NUM_ELEMENTS>(m_output_file, " "));
		m_output_file << endl;

		node_shape = cluster->get_node_shape();
	}

	debugsep;
}
