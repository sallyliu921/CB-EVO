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



#ifndef options_H
#define options_H

#include <string>
using namespace std;

//
// Class_name OPTIONS
//
// Description
//		Reads in and holds the user's options
//		

typedef short K_TYPE;

class OPTIONS
{
public:
	enum TYPE_OF_PARTITIONING {RECURSIVE_BI, KWAY};

	OPTIONS();
	OPTIONS(const OPTIONS & another_options);
	OPTIONS & operator=(const OPTIONS & another_options);
	~OPTIONS();
	void process_options(int argc,char **argv);

	void print_options() const;

	string	get_circuit_name() const { return m_circuit_name;}
	string	get_output_file_name() const { return m_output_file_name; }
	string	get_input_file_name() const { return m_input_file_name;}

	bool 	is_verbose() const { return m_verbose; }
	bool	is_no_warn() const { return m_no_warn; }
	bool	is_quiet() const 	 { return m_quiet; }

	bool	is_draw_circuit() const { return m_draw; }
	bool    is_determine_wirelength_approx() const { return m_determine_wirelength_approx; }



	K_TYPE	get_k() const { return m_k;}
	int		get_nPartitions() const { return m_nPartitions; }
	int 	get_ub_factor() const { return m_ubfactor; }

	TYPE_OF_PARTITIONING 	get_type_of_partitioning() const { return m_partitioning_type;}
	



	/* unimportant options */

	bool 	is_store_luts() const { return m_store_luts;}
	bool	is_expand_luts() const { return m_expand_luts;}

	bool 	is_display_pi_and_dff_distributions() const { return m_display_pi_and_dff_distributions; }
	bool 	is_display_inter_cluster_matricies_at_each_edge_length() const 
				{ return m_display_inter_cluster_matricies_at_each_edge_length; } 
	bool 	is_display_statistics_on_delay_defining_edges() const 
				{ return m_display_statistics_on_delay_defining_edges; }
private:

	string					m_input_file_name;	
	string					m_output_file_name;	
	string 					m_circuit_name;

    K_TYPE  				m_k;					// define LUT-size for analysis
	TYPE_OF_PARTITIONING	m_partitioning_type;	// what kind of partitioning to do
	int						m_nPartitions;			// how many clusters to create
	int 					m_ubfactor;				// balancing factor. defined differently for 
													//   k-way vs. bi-partitioning

	bool					m_verbose;
	bool					m_no_warn;
	bool					m_quiet;



	void display_option_usage() const;
	void read_arguments(int argc, char ** argv);
	string get_circuit_name_from_filename(const string & file_name);


	// measure the wirelength approx
	bool m_determine_wirelength_approx;

	bool m_draw; 		// draw the circuit




	/* unimportant options */


	bool m_display_pi_and_dff_distributions;
	bool m_display_inter_cluster_matricies_at_each_edge_length;
	bool m_display_statistics_on_delay_defining_edges;
	
	// the lut contents are not used as present
	// these functions are for possible future work
	bool m_store_luts;				// store the luts are not. depends on k
	bool m_expand_luts;				// expand luts with don't cares

	bool additional_arguments(const int& argnum, const int& argc, const string& arg) const;

	// unimplemented options
	// lut_analysis
};


#endif
