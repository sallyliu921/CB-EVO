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



#ifndef SPEC_H
#define SPEC_H

#include "types.h"
#include <string>
#include "matrix.h"
#include <fstream>

// Class_name CLUSTER_SPEC
//
// Description
//
//  Specification structure. Holds the specs for the cluster
//

class CLUSTER_SPEC;

typedef vector<CLUSTER_SPEC *> CLUSTER_SPECS;
typedef vector<DELAY_TYPE> DELAYS;

class CLUSTER_SPEC
{
public:
	CLUSTER_SPEC(const CLUSTER_NUMBER_TYPE& cluster_number, const NUM_ELEMENTS& kin,
				 const DELAY_TYPE& max_delay);
	~CLUSTER_SPEC();
	CLUSTER_SPEC(const CLUSTER_SPEC & another_comb_spec);
	CLUSTER_SPEC & operator=(const CLUSTER_SPEC & another_comb_spec);
	friend inline ostream& operator<<(ostream& str, const CLUSTER_SPEC& comb_spec);
	void read_specs(ifstream & spec_file);

	NUM_ELEMENTS	get_nEdges() const 
		{ return m_nIntra_cluster_edges + m_nInter_cluster_input_edges + m_nInter_cluster_output_edges; }

	string 				get_name() const { return m_name; }  
	bool 				is_exact() const { return m_is_exact; }
	NUM_ELEMENTS		get_nNodes() const { return m_nNodes; } 
	NUM_ELEMENTS		get_nIntra_cluster_edges() const { return m_nIntra_cluster_edges; }
	NUM_ELEMENTS		get_nInter_cluster_input_edges() const { return m_nInter_cluster_input_edges; } 
	NUM_ELEMENTS		get_nInter_cluster_output_edges() const { return m_nInter_cluster_output_edges; }
	NUM_ELEMENTS		get_nLatched() const { return m_nLatched; }
	NUM_ELEMENTS		get_nDFF() const { return m_nDFF; } 
	NUM_ELEMENTS		get_nPI() const { return m_nPI; } 
	NUM_ELEMENTS		get_nPO() const { return m_nPO; } 
	NUM_ELEMENTS		get_max_fanout() const { return m_max_fanout; } 
	NUM_ELEMENTS		get_kin() const { return m_kin; } 

	DELAY_TYPE			get_delay() const { return m_delay; } 
	CLUSTER_NUMBER_TYPE get_cluster_number() const { return m_cluster_number; }
	COST_TYPE			get_wirelength() const { return m_wirelength; }			// new iterative wirelength
	NUM_ELEMENTS		get_mike_hutton_locality() const { return m_mike_hutton_locality; }

	SHAPE				get_shape() const { return m_shape; }
	SHAPE				get_latched_shape() const { return m_latched_shape; }
	SHAPE				get_PO_shape() const { return m_PO_shape; }
	DISTRIBUTION		get_intra_cluster_edge_lengths() const 	
								{ return m_intra_cluster_edge_lengths; }
	DISTRIBUTION		get_inter_cluster_input_edge_lengths() const 
								{ return	m_inter_cluster_input_edge_lengths; }
	DISTRIBUTION		get_inter_cluster_output_edge_lengths() const 
								{ return	m_inter_cluster_output_edge_lengths; }
	DISTRIBUTION		get_fanout_distribution() const { return m_fanout_distribution; }
	DISTRIBUTION		get_pi_fanout_distribution() const { return m_pi_fanout_distribution; }
	DISTRIBUTION		get_dff_fanout_distribution() const { return m_dff_fanout_distribution; }

	double				get_mean_fanout_pi() const { return m_mean_fanout_pi; }
	double				get_mean_fanout_dff() const { return m_mean_fanout_dff; }
	double				get_std_dev_fanout_pi() const { return m_std_dev_fanout_pi; }
	double				get_std_dev_fanout_dff() const { return m_std_dev_fanout_dff; }

	SHAPE		get_Inter_cluster_input_shape() const { return m_inter_cluster_input_shape; }
	SHAPE		get_Inter_cluster_output_shape() const { return m_inter_cluster_output_shape; }
	////  Below only if exact is set
	SHAPE		get_level_in() const { return m_level_in; }
	SHAPE		get_level_out() const { return m_Level_out; }		


private:

	void read_number_elements(const char * text_field, long& variable, ifstream& spec_file);
	void read_global_comb_spec(ifstream & spec_file);
	void read_number_elements_vector(const char * text_field, NUM_ELEMENTS_VECTOR & distribution,
				     ifstream & spec_file);
	void read_number_elements_vector(const char * text_field, NUM_ELEMENTS_VECTOR & distribution,
				     ifstream & spec_file, const NUM_ELEMENTS& sum_to_match);
	void read_fanout_distribution(const char * text_field, DISTRIBUTION & distribution, ifstream & spec_file);
	void read_edge_lengths_by_delay_level(ifstream & spec_file);
	void read_degree_info(ifstream & spec_file);
	void read_mean_std_dev_pair(ifstream & spec_file, const string & text_to_find,
								double & mean, double & std_dev);
	void read_global_cluster_specs(ifstream & spec_file);
	void read_shape_specs(ifstream & spec_file);
	void determine_hutton_locality();

	string 			m_name;  
	bool 			m_is_exact;
	CLUSTER_NUMBER_TYPE m_cluster_number;
	NUM_ELEMENTS	m_nNodes, m_nCombinational, m_nDFF, m_nLatched, m_nPI, m_nPO;
	NUM_ELEMENTS	m_nIntra_cluster_edges,m_nInter_cluster_edges,
					m_nInter_cluster_input_edges, m_nInter_cluster_output_edges, 
					m_max_fanout, m_kin;
	NUM_ELEMENTS	m_mike_hutton_locality;
	DELAY_TYPE		m_delay;
	COST_TYPE		m_wirelength;

	double			m_mean_fanout_pi;
	double			m_mean_fanout_dff;
	double			m_std_dev_fanout_pi;
	double			m_std_dev_fanout_dff;

	SHAPE			m_shape;
	SHAPE			m_latched_shape;
	SHAPE			m_PO_shape;
	DISTRIBUTION	m_intra_cluster_edge_lengths;
	DISTRIBUTION	m_inter_cluster_input_edge_lengths;
	DISTRIBUTION	m_inter_cluster_output_edge_lengths;
	DISTRIBUTION	m_fanout_distribution;
	
	////  Below only if exact is set
	SHAPE			m_level_in;
	SHAPE			m_Level_out;		

	SHAPE			m_inter_cluster_input_shape;
	SHAPE			m_inter_cluster_output_shape;

	DISTRIBUTION	m_pi_fanout_distribution;
	DISTRIBUTION	m_dff_fanout_distribution;
};


//
// Class_name CIRCUIT_SPEC
//
// Description
//
//  Specification structures hold the user specifications of the circuit.
//


class CIRCUIT_SPEC
{
public:
	CIRCUIT_SPEC();
	~CIRCUIT_SPEC();
	CIRCUIT_SPEC(const CIRCUIT_SPEC & another_circuit_spec);
	CIRCUIT_SPEC & operator=(const CIRCUIT_SPEC & another_circuit_spec);
	friend inline ostream& operator<<(ostream& str, const CIRCUIT_SPEC& circuit_spec);

	void read_spec();


	bool is_exact() const { return m_is_exact; }

	string 				get_name() const { return m_name; }  
	NUM_ELEMENTS		get_nNodes() const  { return m_nNodes;}
	NUM_ELEMENTS		get_nEdges() const  { return m_nEdges;}
	NUM_ELEMENTS		get_nDFF() const  { return m_nDFF;}
	NUM_ELEMENTS		get_nComb() const  { return m_nCombinational;}
	NUM_ELEMENTS		get_nPI() const  { return m_nPI;} 
	NUM_ELEMENTS		get_nPO() const  { return m_nPO;} 
	NUM_ELEMENTS		get_max_fanout() const  { return m_max_fanout;}
	NUM_ELEMENTS		get_kin() const  { return m_kin;}
	NUM_ELEMENTS		get_nClusters() const  { return m_nClusters; }

	DELAY_TYPE			get_max_comb_delay() const  { return m_max_combinational_delay;}
	
	CLUSTER_SPEC*		get_cluster_spec(const NUM_ELEMENTS & index);

	MATRIX*				get_Comb() { return m_Comb; }
	MATRIX*				get_Latched() { return m_Latched;}
	double				get_mean_fanout_pi() const { return m_mean_fanout_pi; }
	double				get_mean_fanout_dff() const { return m_mean_fanout_dff; }
	double				get_std_dev_fanout_pi() const { return m_std_dev_fanout_pi; }
	double				get_std_dev_fanout_dff() const { return m_std_dev_fanout_dff; }
	double				get_global_wirelength() const { return m_global_wirelength; }

private:
	string 				m_name;
	bool				m_is_exact;
	NUM_ELEMENTS		m_nNodes;
	NUM_ELEMENTS		m_nEdges;
	NUM_ELEMENTS		m_nDFF;
	NUM_ELEMENTS		m_nCombinational;
	NUM_ELEMENTS		m_nPI; 
	NUM_ELEMENTS		m_nPO; 
	NUM_ELEMENTS		m_max_fanout;
	NUM_ELEMENTS		m_kin;
	double				m_global_wirelength;
	DELAY_TYPE			m_max_combinational_delay;
	bool				m_clock;
	double				m_mean_fanout_pi;
	double				m_mean_fanout_dff;
	double				m_std_dev_fanout_pi;
	double				m_std_dev_fanout_dff;

	NUM_ELEMENTS		m_nClusters;
	CLUSTER_SPECS		m_cluster_specs;
	MATRIX*				m_Comb;
	MATRIX*				m_Latched;

	void read_global_circuit_spec(ifstream& spec_file);
	void read_cluster_summary_spec(ifstream & spec_file);
	void read_shape_info(ifstream& spec_file);
	void read_degree_info(ifstream& spec_file);
	void read_sequential_info(ifstream& spec_file);
	void read_cluster_specs(ifstream& spec_file);
	void read_inter_cluster_connection_matrix(ifstream& spec_file);
	void read_matrix(ifstream& spec_file, MATRIX* matrix);
	void read_mean_std_dev_pair(ifstream & spec_file, const string & text_to_find,
								double & mean, double & std_dev);
	void read_number_elements(const char * text_field, long& variable, ifstream& spec_file);
	void read_number_elements(const char * text_field, short& variable, ifstream& spec_file);
};

/*
inline ostream& operator<<(ostream& stream, const CLUSTER_SPEC& comb_spec)
{
	stream 	<< comb_spec.name << " " << endl;
stream 	<< comb_spec.is_exact  << " " << comb_spec.nNodes  << " " 
			<< comb_spec.nIntra_cluster_edges << " " << comb_spec.nInter_cluster_edges << " " 
			<< comb_spec.nDFF << " " << comb_spec.nPI  << " " << comb_spec.nPO  << " "
			<< comb_spec.nGI  << " " << comb_spec.nGO  << " " << comb_spec.max_fanout << " " 
			<< comb_spec.kin << " " << comb_spec.delay  << " " 
			<< " " << comb_spec.wirelength << endl;

	copy(comb_spec.shape.begin(), comb_spec.shape.end(), ostream_iterator<NUM_ELEMENTS>(stream, " "));
	stream << endl;
	copy(comb_spec.DFF_shape.begin(), comb_spec.DFF_shape.end(), ostream_iterator<NUM_ELEMENTS>(stream, " "));
	stream << endl;
	copy(comb_spec.PO_shape.begin(), comb_spec.PO_shape.end(), ostream_iterator<NUM_ELEMENTS>(stream, " "));
	stream << endl;
	copy(comb_spec.GI_shape.begin(), comb_spec.GI_shape.end(), ostream_iterator<NUM_ELEMENTS>(stream, " "));
	stream << endl;
	copy(comb_spec.GO_shape.begin(), comb_spec.GO_shape.end(), ostream_iterator<NUM_ELEMENTS>(stream, " "));
	stream << endl;
	copy(comb_spec.edge_lengths.begin(), comb_spec.edge_lengths.end(), 
			ostream_iterator<NUM_ELEMENTS>(stream, " "));
	stream << endl;
	copy(comb_spec.inter_cluster_edge_lengths.begin(), comb_spec.inter_cluster_edge_lengths.end(), 
			ostream_iterator<NUM_ELEMENTS>(stream, " "));
	stream << endl;
	copy(comb_spec.fanout_distribution.begin(), comb_spec.fanout_distribution.end(), 
			ostream_iterator<NUM_ELEMENTS>(stream, " "));
	stream << endl;
	copy(comb_spec.Level_in.begin(), comb_spec.Level_in.end(), ostream_iterator<NUM_ELEMENTS>(stream, " "));
	stream << endl;
	copy(comb_spec.Level_out.begin(), comb_spec.Level_out.end(), ostream_iterator<NUM_ELEMENTS>(stream, " "));
	stream << endl;
	copy(comb_spec.inter_cluster_ghost_input_shape.begin(), comb_spec.inter_cluster_ghost_input_shape.end(), 
			ostream_iterator<NUM_ELEMENTS>(stream, " "));
	stream << endl;
	copy(comb_spec.inter_cluster_ghost_output_shape.begin(), comb_spec.inter_cluster_ghost_output_shape.end(), 
			ostream_iterator<NUM_ELEMENTS>(stream, " "));
	stream << endl;
	copy(comb_spec.inter_cluster_input_shape.begin(), comb_spec.inter_cluster_input_shape.end(), 
			ostream_iterator<NUM_ELEMENTS>(stream, " "));
	stream << endl;
	copy(comb_spec.inter_cluster_output_shape.begin(), comb_spec.inter_cluster_output_shape.end(), 
			ostream_iterator<NUM_ELEMENTS>(stream, " "));
	stream << endl;

	return stream;
}

inline ostream& operator<<(ostream& stream, const CLUSTER_SPEC& cluster_spec)
{
	int i;
	stream 	<< cluster_spec.m_name << " " << endl;
	stream 	<< cluster_spec.m_is_exact  << " " << cluster_spec.m_nNodes  << " " 
			<< cluster_spec.m_nIntra_cluster_edges << " " << cluster_spec.m_nInter_cluster_edges << " " 
			<< cluster_spec.m_nDFF << " " << cluster_spec.m_nPI  << " " << cluster_spec.m_nPO  << " "
			<< cluster_spec.m_max_fanout << " " << cluster_spec.m_kin << " " << cluster_spec.m_delay  << endl;

	for (i=0; i < cluster_spec.m_comb_specs.size(); i++)
	{
		stream << cluster_spec.m_comb_specs[i] << endl;
	}

	return stream;
	

}

inline ostream& operator<<(ostream& stream, const CIRCUIT_SPEC& circuit_spec)
{
	int i;
	stream 	<< circuit_spec.m_name << " " << endl;
	stream 	<< circuit_spec.m_is_exact  << " " << circuit_spec.m_nNodes  << " " 
			<< circuit_spec.m_nEdges << " " << circuit_spec.m_nDFF << " " 
			<< circuit_spec.m_nPI  << " " << circuit_spec.m_nPO  << " "
			<< circuit_spec.m_max_fanout << " " << circuit_spec.m_kin << " " << circuit_spec.m_delay  << endl;
	stream 	<< circuit_spec.m_nClusters;

	for (i=0; i < circuit_spec.m_cluster_specs.size(); i++)
	{
		stream << circuit_spec.m_cluster_specs[i] << endl;
	}

	assert(false);

	//stream << circuit_spec.m_inter_cluster_connection_matrix;
	return stream;

}
*/
#endif
