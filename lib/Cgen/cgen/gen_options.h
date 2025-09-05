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
// Class_name GEN_OPTIONS
//
// Description
//		Reads in and holds the user's options
//		


class GEN_OPTIONS
{
public:
	GEN_OPTIONS();
	GEN_OPTIONS(const GEN_OPTIONS & another_options);
	GEN_OPTIONS & operator=(const GEN_OPTIONS & another_options);
	~GEN_OPTIONS();
	void process_options(int argc,char **argv);

	void print_options();

	string	get_input_file_name() const { return m_input_file_name;}
	string	get_result_file_name() const { return m_result_file_name; }


	long get_random_seed() const { return m_random_seed; }

	bool	is_output_vhdl() const { return m_output_vhdl; }
	bool 	is_verbose() const { return m_verbose; }
	bool	is_no_warn() const { return m_no_warnings; }
	bool	is_quiet() 	 const { return m_quiet; }
	bool	is_draw_graph()  const { return m_draw_graph; }
	bool 	is_fast_gen() const { return m_is_fast_gen; }

	double 	get_delay_structure_alpha() const { return m_alpha; }
	double 	get_delay_structure_beta() const { return m_beta; }
	double 	get_delay_structure_gamma() const { return m_gamma; }
	double 	get_delay_structure_lamda() const { return m_lamda; }
	double 	get_delay_structure_phi() const { return m_phi; }
	double	get_delay_structure_init_too_many_inputs_factor() const { return m_too_many_inputs_factor; }
	double 	get_delay_structure_init_too_few_outputs_factor() const { return m_too_few_outputs_factor; }
	double  get_delay_structure_init_too_few_inputs_factor() const { return m_too_few_inputs_factor; }
	double	get_delay_structure_multiplier_too_many_inputs_factor() const 
				{ return m_mult_too_many_inputs_factor; }
	double	get_delay_structure_multiplier_too_few_outputs_factor() const 
				{ return m_mult_too_few_outputs_factor; }
	double	get_delay_structure_multiplier_too_few_inputs_factor() const 
				{ return m_mult_too_few_inputs_factor ; }
	double get_delay_structure_init_temperature() const { return m_delay_structure_init_temperature; }

	double 	get_degree_partitioning_delta() const { return m_delta; }
	double  get_degree_init_temperature() const { return m_degree_init_temperature; }


	double 	get_final_edge_assignment_eta() const { return m_eta; }
	double  get_final_edge_assignment_init_temperature() const 
			   { return m_final_edge_assignment_init_temperature; }
	double  get_final_edge_assignment_dff_loop_penalty_multiplier() const 
			   { return m_dff_loop_penalty_multipler; }
	double	get_wirelength_wanted() const { return m_wirelength_wanted; }
	double  get_final_edge_assignment_mult_double_inputs_factor() const { return m_final_edge_assignment_mult_double_inputs_factor; }

	bool	is_best_wirelength_settings() const { return m_best_wirelength_settings; }
	bool	is_wirelength_from_stat_file() const { return m_wirelength_from_stat_file; }
	bool	is_accept_initial_wirelength() const { return m_accept_initial_wirelength; }
	bool	is_user_wirelength() const { return m_user_wirelength; }


	bool 	is_read_level0_fanout() const { return m_read_level0_fanout; }
private:

	long m_random_seed;			//  random seed  

	string m_input_file_name;	// the name of the input file
	string m_result_file_name;	// the name of the output file

	bool m_verbose;         	// Be more verbose in the messages output
	bool m_quiet;
	bool m_no_warnings;			// Supress warnings
	bool m_draw_graph;			// Draw a dot drawing of the clone
	bool m_output_vhdl;			// output a vhdl file
	bool m_is_fast_gen;			// faster but less accurate. not tested

	// Delay structure creation options
	double		m_alpha, m_beta,m_gamma;
	double		m_too_many_inputs_factor;
	double 		m_too_few_outputs_factor;
	double 		m_too_few_inputs_factor;
	double		m_mult_too_many_inputs_factor;
	double 		m_mult_too_few_outputs_factor;
	double 		m_mult_too_few_inputs_factor;
	double		m_delay_structure_init_temperature;

	// Degree Partitioning options
	double		m_delta;
	double 		m_degree_init_temperature;

	// Final Edge assignment options
	double		m_eta;
	double		m_dff_loop_penalty_multipler;
	double 		m_final_edge_assignment_init_temperature;
	double      m_final_edge_assignment_mult_double_inputs_factor;
  

	// Wirelength-approx options
	bool m_wirelength_from_stat_file;
	bool m_accept_initial_wirelength;
	bool m_best_wirelength_settings;// best for comb is to accept init. 
									// best for seq is to go for lowest wirelength
	bool m_user_wirelength;			// desired wirelength is the users wirelength
	double 	m_wirelength_wanted; // the users wirelength if supplied



	// these options are not used yet.
	bool m_should_generate_k_reg_graph;			//  Should we generate a k-regular graph?  
	int  m_kreg_n, m_kreg_i, m_kreg_o, m_kreg_e, m_kreg_d;  	//  k-regular parameters  
	bool m_choose_new_seed_for_every_circuit;        			//  Choose new seed for every circuit 
	bool m_sanity_checks_on;				//  Fail on sanity checks (off == just warn). 

	void display_option_usage();
	bool additional_arguments(const int& argnum, const int& argc, const string& arg) const;

	bool 	m_read_level0_fanout;
	double	m_lamda,m_phi;			 // unused
};

#endif
