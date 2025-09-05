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




#include "spec.h"
#include "gen.h"
#include <algorithm>
#include <math.h>

long g_line_number;

void is_same_text(const string& string_have, const char * string_wanted)
{
	if (string_wanted != string_have)
	{
		debug("\nWanted string: " << string_wanted << ". Got string " << string_have << ".");
		debug("Look around line number " << g_line_number << endl);
		assert(false);
	}
}

void is_same_text(const string& string_have, const string&  string_wanted)
{
	if (string_wanted != string_have)
	{
		debug("\nWanted string: " << string_wanted << ". Got string " << string_have << ".");
		debug("Look around line number " << g_line_number << endl);
		assert(false);
	}
}
void is_same_text(const char& string_have, const char&  string_wanted)
{
	if (string_wanted != string_have)
	{
		debug("\nWanted string: " << string_wanted << ". Got string " << string_have << ".");
		debug("Look around line number " << g_line_number << endl);
		assert(false);
	}
}



CLUSTER_SPEC::CLUSTER_SPEC
(
	const CLUSTER_NUMBER_TYPE& cluster_number,
	const NUM_ELEMENTS& kin,
	const DELAY_TYPE& max_delay
)
{
	m_is_exact	= false;
	m_nNodes		= 0;
	m_nIntra_cluster_edges = 0;
	m_nInter_cluster_input_edges = 0;
	m_nInter_cluster_output_edges = 0;
	m_nDFF			= 0;
	m_nLatched		= 0;
	m_nPI			= 0;
	m_nPO			= 0;
	m_max_fanout	= 0;
	m_kin			= kin;
	m_delay			= max_delay;
	m_cluster_number= cluster_number;
	m_wirelength	= 0;
	m_mike_hutton_locality = 0;
}

CLUSTER_SPEC::CLUSTER_SPEC(const CLUSTER_SPEC & another_comb_spec)
{
	m_is_exact		= another_comb_spec.m_is_exact;
	m_nNodes		= another_comb_spec.m_nNodes;
	m_nIntra_cluster_edges			= another_comb_spec.m_nIntra_cluster_edges;
	m_nInter_cluster_input_edges	= another_comb_spec.m_nInter_cluster_input_edges;
	m_nInter_cluster_output_edges	= another_comb_spec.m_nInter_cluster_output_edges;
	m_nDFF			= another_comb_spec.m_nDFF;
	m_nLatched		= another_comb_spec.m_nLatched;
	m_nPI			= another_comb_spec.m_nPI;	   
	m_nPO			= another_comb_spec.m_nPO; 
	m_max_fanout	= another_comb_spec.m_max_fanout;	
	m_kin			= another_comb_spec.m_kin;		   
	m_delay			= another_comb_spec.m_delay;
	m_wirelength		= another_comb_spec.m_wirelength;   

	m_shape			= another_comb_spec.m_shape;
	m_latched_shape	= another_comb_spec.m_latched_shape;
	m_PO_shape		= another_comb_spec.m_PO_shape;
	m_fanout_distribution				= another_comb_spec.m_fanout_distribution;	
	m_level_in							= another_comb_spec.m_level_in;	
	m_Level_out							= another_comb_spec.m_Level_out;
	m_inter_cluster_input_shape 		= another_comb_spec.m_inter_cluster_input_shape;
	m_inter_cluster_output_shape 		= another_comb_spec.m_inter_cluster_output_shape;

	m_intra_cluster_edge_lengths 		= another_comb_spec.m_intra_cluster_edge_lengths; 
	m_inter_cluster_input_edge_lengths	= another_comb_spec.m_inter_cluster_input_edge_lengths;
	m_inter_cluster_output_edge_lengths = another_comb_spec.m_inter_cluster_output_edge_lengths;
}
CLUSTER_SPEC & CLUSTER_SPEC::operator=(const CLUSTER_SPEC & another_comb_spec)
{

	m_is_exact		= another_comb_spec.m_is_exact;
	m_nNodes		= another_comb_spec.m_nNodes;
	m_nIntra_cluster_edges	= another_comb_spec.m_nIntra_cluster_edges;
	m_nInter_cluster_input_edges	= another_comb_spec.m_nInter_cluster_input_edges;
	m_nInter_cluster_output_edges	= another_comb_spec.m_nInter_cluster_output_edges;
	m_nDFF			= another_comb_spec.m_nDFF;
	m_nLatched		= another_comb_spec.m_nLatched;
	m_nPI			= another_comb_spec.m_nPI;	   
	m_nPO			= another_comb_spec.m_nPO; 
	m_max_fanout	= another_comb_spec.m_max_fanout;	
	m_kin			= another_comb_spec.m_kin;		   
	m_delay			= another_comb_spec.m_delay;
	m_wirelength	= another_comb_spec.m_wirelength;   

	m_shape			= another_comb_spec.m_shape;
	m_latched_shape	= another_comb_spec.m_latched_shape;
	m_PO_shape		= another_comb_spec.m_PO_shape;
	m_fanout_distribution		= another_comb_spec.m_fanout_distribution;	
	m_level_in					= another_comb_spec.m_level_in;	
	m_Level_out					= another_comb_spec.m_Level_out;
	m_inter_cluster_input_shape = another_comb_spec.m_inter_cluster_input_shape;
	m_inter_cluster_output_shape= another_comb_spec.m_inter_cluster_output_shape;

	m_intra_cluster_edge_lengths 		= another_comb_spec.m_intra_cluster_edge_lengths; 
	m_inter_cluster_input_edge_lengths	= another_comb_spec.m_inter_cluster_input_edge_lengths;
	m_inter_cluster_output_edge_lengths = another_comb_spec.m_inter_cluster_output_edge_lengths;

	return (*this);
}

CLUSTER_SPEC::~CLUSTER_SPEC()
{





}

void CLUSTER_SPEC::read_specs
(
	ifstream & spec_file
)
{
	read_global_cluster_specs(spec_file);
	read_degree_info(spec_file);
	read_shape_specs(spec_file);
	determine_hutton_locality();
}

void CLUSTER_SPEC::read_shape_specs
(
	ifstream & spec_file
)
{
	string text;

	getline(spec_file, text, '\n');	// get the rest of the line
	is_same_text(text,"======================== SHAPE ============================");
	g_line_number++;

	read_number_elements_vector("Node_shape:", m_shape, spec_file, m_nNodes);
	read_number_elements_vector("Input_shape:", m_level_in, spec_file);
	read_number_elements_vector("Output_shape:", m_Level_out, spec_file);
	read_number_elements_vector("Latched_shape:", m_latched_shape, spec_file, m_nLatched);
	read_number_elements_vector("POshape:", m_PO_shape, spec_file, m_nPO);
	read_number_elements_vector("Intra_cluster_edge_length_distribution:", 
								m_intra_cluster_edge_lengths, spec_file, m_nIntra_cluster_edges);

	read_number_elements_vector("Inter_cluster_input_edge_length_distribution:", 
								m_inter_cluster_input_edge_lengths, spec_file, m_nInter_cluster_input_edges);

	read_number_elements_vector("Inter_cluster_output_edge_length_distribution:", 
								m_inter_cluster_output_edge_lengths, spec_file, m_nInter_cluster_output_edges);

	read_fanout_distribution("Fanout_distribution:", m_fanout_distribution, spec_file);

}



void CLUSTER_SPEC::determine_hutton_locality()
{
	// need to choose between the three localities. have no idea of the difference
	assert(new_gen);
	// comb.gen locality
    double loc1 = log(m_nNodes)/log(2);
    double loc2 = sqrt(m_nNodes)/3;
    //LOCALITY loc3 = 2 * exp(log(n)/log(10));// i have no idea what loc3 is used for

	assert(new_gen); // need look at this
    m_mike_hutton_locality = static_cast<NUM_ELEMENTS>(ceil(MAX(loc1, loc2)));

	// fsm.gen locality
    //m_locality = nint(6 + exp(log(n)/log(10))/exp(2));
}
void CLUSTER_SPEC::read_number_elements_vector
(
	const char * text_field, 
	NUM_ELEMENTS_VECTOR & num_element_vector,
    ifstream & spec_file
)
{
	int i = 0;
	string text;
	NUM_ELEMENTS number;

	debugif(DSPEC,"Reading " << text_field);

	spec_file >> text; is_same_text(text,text_field);
	spec_file >> text; is_same_text(text,"(");

	for (i = 0; i <= m_delay; i++)
	{
		spec_file >> number;
		num_element_vector.push_back(number);
	}

	spec_file >> text; is_same_text(text,")");


	g_line_number++;
}



void CLUSTER_SPEC::read_number_elements_vector
(
	const char * text_field, 
	NUM_ELEMENTS_VECTOR & num_element_vector,
    ifstream & spec_file,
	const NUM_ELEMENTS& sum_to_match
)
{
	int i = 0;
	string text;
	NUM_ELEMENTS number,
				 sum = 0;

	debugif(DSPEC,"Reading " << text_field);


	spec_file >> text; is_same_text(text,text_field);
	spec_file >> text; is_same_text(text,"(");

	for (i = 0; i <= m_delay; i++)
	{
		spec_file >> number;
		num_element_vector.push_back(number);
		sum += number;
	}

	spec_file >> text; is_same_text(text,")");

	if (sum != sum_to_match)
	{
		Fail("The sum of vector " << text_field << " around line number " << g_line_number <<
				" is not equal to its sum above.");
	}
			 


	g_line_number++;
}
	
void CLUSTER_SPEC::read_fanout_distribution
(
	const char * text_field, 
	DISTRIBUTION & distribution,
    ifstream & spec_file
)
{
	string text;
	NUM_ELEMENTS number;

	spec_file >> text; is_same_text(text,text_field);
	spec_file >> text; is_same_text(text, "(");

	debugif(DSPEC,"The size of max fanout is " << m_max_fanout);

	for (int i = 0; i <= m_max_fanout; i++)
	{
		spec_file >> number;
		distribution.push_back(number);
	}
	//copy_n(istream_iterator<NUM_ELEMENTS>(spec_file), max_fanout+1, back_inserter(distribution));  

	spec_file >> text; is_same_text(text,")");

	getline(spec_file, text, '\n');	// get the endl

	g_line_number++;
}

// just eat the edge lengths by delay level for now
void CLUSTER_SPEC::read_edge_lengths_by_delay_level
(
	ifstream & spec_file
)
{
	DELAY_TYPE delay_index = 0;	
	DELAY_TYPE delay_found = 0;	
	NUM_ELEMENTS number;
	string text;
	int i;

	spec_file >> text;
	is_same_text(text,"edge_lengths_by_delay_level:");

	for (delay_index = 0; delay_index <= m_delay; delay_index++)
	{
		spec_file >> text;
		is_same_text(text,"Delay");
		spec_file >> delay_found;
		assert(delay_found == delay_index);

		spec_file >> text;
		is_same_text(text,"Intra_cluster_input_edge_lengths");
		for (i = 0; i <= m_delay; i++) { spec_file >> number; }

		spec_file >> text;
		is_same_text(text,"Intra_cluster_output_edge_lengths");
		for (i = 0; i <= m_delay; i++) { spec_file >> number; }

		spec_file >> text;
		is_same_text(text,"Inter_cluster_input_edge_lengths");
		for (i = 0; i <= m_delay; i++) { spec_file >> number; }

		spec_file >> text;
		is_same_text(text,"Inter_cluster_output_edge_lengths");
		for (i = 0; i <= m_delay; i++) { spec_file >> number; }
	}
}




void CLUSTER_SPEC::read_global_cluster_specs
(
	ifstream & spec_file
)
{
	string text;

	spec_file >> text; is_same_text(text,"####################");
	spec_file >> text; is_same_text(text,"Cluster");
	spec_file >> m_cluster_number;
	debugif(DSPEC,"Cluster number is " << m_cluster_number);
	spec_file >> text;
	g_line_number++;

	read_number_elements("Number_of_Nodes:", m_nNodes, spec_file);
	read_number_elements("Number_of_Intra_cluster_edges:", m_nIntra_cluster_edges, spec_file);
	read_number_elements("Number_of_Inter_cluster_edges:", m_nInter_cluster_edges, spec_file);
	read_number_elements("Number_of_PI:", m_nPI, spec_file);
	read_number_elements("Number_of_PO:", m_nPO, spec_file);
	read_number_elements("Number_of_Comb:", m_nCombinational, spec_file);
	read_number_elements("Number_of_DFF:", m_nDFF, spec_file);
	read_number_elements("Number_of_Latched:", m_nLatched, spec_file);
	read_number_elements("Number_of_Inter_cluster_input_edges:", m_nInter_cluster_input_edges, spec_file);
	read_number_elements("Number_of_Inter_cluster_output_edges:", m_nInter_cluster_output_edges, spec_file);
	getline(spec_file, text, '\n'); // get the end of line

	char next_char = spec_file.peek();

 	if (next_char == 'W')
	{
 		spec_file >> text; is_same_text(text, "Wirelength_approx:");

		spec_file >> m_wirelength;
		getline(spec_file, text, '\n'); // get the rest of the line
		g_line_number++;
	}
}

void CLUSTER_SPEC::read_degree_info
(
	ifstream & spec_file
)
{
	string text;

	getline(spec_file, text, '\n'); 
	is_same_text(text,"======================== DEGREE ============================");
	g_line_number++;

	getline(spec_file, text, '\n'); // read the avg_fanin_comb
	getline(spec_file, text, '\n'); // read the avg_fanout
	getline(spec_file, text, '\n'); // read the avg_fanout_comb
	g_line_number += 3;

	read_mean_std_dev_pair(spec_file, "Avg_fanout_pi:", m_mean_fanout_pi, m_std_dev_fanout_pi);
	read_mean_std_dev_pair(spec_file, "Avg_fanout_dff:", m_mean_fanout_dff, m_std_dev_fanout_dff);

	read_number_elements("Maximum_fanout:", m_max_fanout, spec_file);
	debugif(DSPEC,"Found max_fanout:.  Max fanout is: " << m_max_fanout);
	getline(spec_file, text, '\n'); // get the rest of the max_fanout line

	getline(spec_file, text, '\n');
	getline(spec_file, text, '\n');
	getline(spec_file, text, '\n');
	getline(spec_file, text, '\n');
	getline(spec_file, text, '\n');
	getline(spec_file, text, '\n');
	g_line_number += 6;

	if (g_options->is_read_level0_fanout())
	{
		read_number_elements_vector("pi_fanout:", m_pi_fanout_distribution, spec_file);
		read_number_elements_vector("dff_fanout:", m_dff_fanout_distribution, spec_file);
	}
}	
void CLUSTER_SPEC::read_number_elements(const char * text_field, long& variable, ifstream& spec_file)
{
	string text;

	spec_file >> text;

	if (text != text_field)
	{
		debugif(DSPEC,"Should have found " << text_field << " but found " << text);
	}
	is_same_text(text,text_field);
	spec_file >> variable;

	g_line_number++;
}
void CLUSTER_SPEC::read_mean_std_dev_pair
(
	ifstream & spec_file,
	const string & text_to_find,
	double & mean,
	double & std_dev
)
{
	string text;

	spec_file >> text;
	is_same_text(text,text_to_find);

	getline(spec_file, text, '(');	

	if (text == " nan ")
	{
		// we should not have any pi or dff in this circuit
		mean = 0.0;
		std_dev = 0.0;

		// get the rest of the line
		getline(spec_file, text, '\n');	
	}
	else
	{
		// wasn't nan. 

		mean = atof(text.c_str());

		getline(spec_file, text, ')');	

		std_dev = atof(text.c_str());

		assert(new_gen);			// a good thing?
		std_dev = MAX(1, std_dev); 	// give the avg a little std deviation to deal with 
									// our imperfect world

		getline(spec_file, text, '\n');	
	}

	g_line_number++;
}

CIRCUIT_SPEC::CIRCUIT_SPEC()
{
	m_is_exact 	= false;
	m_nNodes	= 0;
	m_nEdges	= 0;
	m_nDFF		= 0;
	m_nCombinational		= 0;
	m_nPI		= 0; 
	m_nPO		= 0; 
	m_max_fanout= 0;
	m_kin		= 0;
	m_global_wirelength = 0;
	m_max_combinational_delay = 0;
	m_nClusters	= 0;
	m_clock		= false;
	m_Comb 		= 0;
	m_Latched	= 0;
}

CIRCUIT_SPEC::~CIRCUIT_SPEC()
{
	CLUSTER_SPECS::iterator cluster_spec_iter = m_cluster_specs.begin();

	for (cluster_spec_iter =  m_cluster_specs.begin(); cluster_spec_iter != m_cluster_specs.end();
			cluster_spec_iter++)
	{
		delete *cluster_spec_iter;
	}

	delete m_Comb;
	delete m_Latched;

	m_Comb = 0;
	m_Latched = 0;
}


// read the specs from the .stats file
//
// PRE: g_options contains the name of the input file
// POST: we have read the specs from the spec file else 
//       something was wrong with the .stats file
//
void CIRCUIT_SPEC::read_spec()
{
	ifstream spec_file;

	debugif(DSPEC,"Opening the spec file");

	assert(g_options);
	string file_name = g_options->get_input_file_name();
	assert(! file_name.empty());

	debug("File name of spec file is " << file_name);

	spec_file.open(file_name.c_str(), ios::in);

	if (! spec_file.good())
	{
		Fail("Something wrong with file " << file_name);
	}

	g_line_number = 1;
	read_global_circuit_spec(spec_file);
	read_cluster_summary_spec(spec_file);
	read_degree_info(spec_file);
	read_shape_info(spec_file);
	read_cluster_specs(spec_file);	
	read_inter_cluster_connection_matrix(spec_file);
	spec_file.close();
}

void CIRCUIT_SPEC::read_global_circuit_spec
(
	ifstream & spec_file
)
{
	string text;
	char next_char;

	getline(spec_file, text, '\n');
	g_line_number++;

	spec_file >> text;
	assert(text == "Circuit_Name:");
	spec_file >> m_name;
	debugif(DSPEC,"Circuit Name: " << m_name);
	g_line_number++;


	read_number_elements("Number_of_Nodes:", m_nNodes, spec_file);
	read_number_elements("Number_of_Edges:", m_nEdges, spec_file);
	read_number_elements("Maximum_Delay:", m_max_combinational_delay, spec_file);
	read_number_elements("Number_of_PI:", m_nPI, spec_file);
	read_number_elements("Number_of_PO:", m_nPO, spec_file);
	read_number_elements("Number_of_Combinational_Nodes:", m_nCombinational, spec_file);
	read_number_elements("Number_of_DFF:", m_nDFF, spec_file);
	read_number_elements("kin:", m_kin, spec_file);
	getline(spec_file, text, '\n'); // get the rest of the line

	next_char = spec_file.peek();
 	if (next_char == 'W')
	{
		spec_file >> text; is_same_text(text,"Wirelength_approx:");
		spec_file >> m_global_wirelength;
		getline(spec_file, text, '\n'); // get the newline
		g_line_number++;
	}

	next_char = spec_file.peek();
 	if (next_char == 'c')
	{
		spec_file >> text; is_same_text(text,"clock:");
		getline(spec_file, text, '\n'); // get the rest of the line
		g_line_number++;
	}
}
void CIRCUIT_SPEC::read_cluster_summary_spec
(
	ifstream & spec_file
)
{
	string text;

	
	getline(spec_file, text, '\n');	// get the rest of the line
	is_same_text(text,"======================== Cluster_Summary ==================");

	getline(spec_file, text, '\n'); // get the rest of the Number_of_Nodes line
	getline(spec_file, text, '\n'); // get the nPI line
	getline(spec_file, text, '\n'); // get the nDFF line
	getline(spec_file, text, '\n'); // get the nIntra_cluster_edges line
	getline(spec_file, text, '\n'); // get the nInter_cluster_edges line
	g_line_number += 6;

	spec_file >> text;
	if (text == "Wirelength_approx:")
	{
		spec_file >> text;
		getline(spec_file, text, '\n'); // get the rest of the wirelength_approx distribution
		g_line_number++;
	}

	getline(spec_file, text, '\n'); // get the partitioned scaled cost line
	g_line_number++;
}


void CIRCUIT_SPEC::read_degree_info
(
	ifstream & spec_file
)
{
	string text;

	getline(spec_file, text, '\n');	// get the rest of the line
	is_same_text(text,"======================== DEGREE ============================");

	getline(spec_file, text, '\n'); // read the avg_fanin_comb
	getline(spec_file, text, '\n'); // read the avg_fanout
	getline(spec_file, text, '\n'); // read the avg_fanout_comb
	g_line_number += 4;

	read_mean_std_dev_pair(spec_file, "Avg_fanout_pi:", m_mean_fanout_pi, m_std_dev_fanout_pi);
	read_mean_std_dev_pair(spec_file, "Avg_fanout_dff:", m_mean_fanout_dff, m_std_dev_fanout_dff);

	read_number_elements("Maximum_fanout:", m_max_fanout, spec_file);
	getline(spec_file, text, '\n'); // rest of max_fanout

	getline(spec_file, text, '\n'); // high comb
	getline(spec_file, text, '\n'); // high pi
	getline(spec_file, text, '\n'); // high dff
	getline(spec_file, text, '\n'); // 10plus comb
	getline(spec_file, text, '\n'); // 10plus pi
	getline(spec_file, text, '\n'); // 10plus dff
	g_line_number += 6;
}	

void CIRCUIT_SPEC::read_shape_info
(
	ifstream & spec_file
)
{
	string text;

	getline(spec_file, text, '\n');	// get the rest of the line
	is_same_text(text,"======================== SHAPE ============================");

	getline(spec_file, text, '\n'); // get the node_shape
	getline(spec_file, text, '\n'); // input shape
	getline(spec_file, text, '\n'); // output shape
	getline(spec_file, text, '\n'); // latched shape
	getline(spec_file, text, '\n'); // poshape
	getline(spec_file, text, '\n'); // edge length
	getline(spec_file, text, '\n'); // fanout 
	g_line_number += 8;
}



void CIRCUIT_SPEC::read_cluster_specs
(
	ifstream & spec_file
)
{
	string text;
	CLUSTER_SPEC * cluster_spec;

	getline(spec_file, text, '\n');
	is_same_text(text, "#################### Clusters ######################");
	g_line_number++;

	read_number_elements("Number_of_Clusters:", m_nClusters, spec_file);
	getline(spec_file, text, '\n'); // eat the rest of the Number of clusters line

	for (int i = 0; i < m_nClusters; i++)
	{
		cluster_spec = new CLUSTER_SPEC(i, m_kin, m_max_combinational_delay);
		assert(cluster_spec);
		cluster_spec->read_specs(spec_file);
		m_cluster_specs.push_back(cluster_spec);
	}

}

void CIRCUIT_SPEC::read_inter_cluster_connection_matrix
(
	ifstream& spec_file
)
{
	string text_file;	
	string text;

	getline(spec_file, text, '\n');
	is_same_text(text, 
			"-------------------- Inter_cluster_adjacentcy_matrix_to_combinational_nodes --------------------");
	g_line_number += 2;

	m_Comb = new MATRIX(m_nClusters, m_nClusters);
	read_matrix(spec_file, m_Comb);

	debugif(DSPEC, "inter cluster connection\n" << *m_Comb);

	getline(spec_file, text, '\n');
	getline(spec_file, text, '\n');
	getline(spec_file, text, '\n');
	getline(spec_file, text, '\n');
	is_same_text(text, 
			"-------------------- Inter_cluster_adjacentcy_matrix_to_dffs -----------------------------------");
	g_line_number += 2;

	m_Latched = new MATRIX(m_nClusters, m_nClusters);
	read_matrix(spec_file, m_Latched);

	debugif(DSPEC,"inter cluster connection for dff \n" << *m_Latched);
}

void CIRCUIT_SPEC::read_matrix
(
	ifstream& spec_file,
	MATRIX* matrix
)
{
	INDEX_SIZE row, column;

	for (row = 0; row < m_nClusters; row++)
	for (column= 0; column < m_nClusters; column++)
	{
		spec_file >> (*matrix)(row, column);
	}

	g_line_number += m_nClusters;
}
CLUSTER_SPEC* CIRCUIT_SPEC::get_cluster_spec(const NUM_ELEMENTS & index)
{ 
	assert(index >= 0 && index < static_cast<signed>(m_cluster_specs.size()));
	return m_cluster_specs[index];
}
void CIRCUIT_SPEC::read_mean_std_dev_pair
(
	ifstream & spec_file,
	const string & text_to_find,
	double & mean,
	double & std_dev
)
{
	string text;

	spec_file >> text; is_same_text(text, text_to_find);

	getline(spec_file, text, '(');	

	// check to see if the next string is number or not a number (nan)
	// nan usually indicates that we divided by zero
	if (text == " nan ")
	{
		// we should not have any pi or dff in this circuit
		mean = 0.0;
		std_dev = 0.0;

		// get the rest of the line
		getline(spec_file, text, '\n');	
	}
	else
	{
		// wasn't nan. 

		mean = atof(text.c_str());

		getline(spec_file, text, ')');	

		std_dev = atof(text.c_str());

		assert(new_gen);			// a good thing?
		std_dev = MAX(1, std_dev); 	// give the avg a little std deviation to deal with 
									// our imperfect world

		getline(spec_file, text, '\n');	
	}
	g_line_number++;
}
void CIRCUIT_SPEC::read_number_elements(const char * text_field, long& variable, ifstream& spec_file)
{
	string text;

	spec_file >> text;

	if (text != text_field)
	{
		debugif(DSPEC,"Should have found " << text_field << " but found " << text);
	}
	is_same_text(text,text_field);
	spec_file >> variable;

	g_line_number++;
}
void CIRCUIT_SPEC::read_number_elements(const char * text_field, short& variable, ifstream& spec_file)
{
	string text;

	spec_file >> text;

	if (text != text_field)
	{
		debugif(DSPEC,"Should have found " << text_field << " but found " << text);
	}
	is_same_text(text,text_field);
	spec_file >> variable;

	g_line_number++;
}
