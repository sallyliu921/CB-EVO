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



#include "gen_options.h"
#include "rand.h"
#include "util.h"
#include <assert.h>

extern RANDOM_NUMBER * g_rand_number_gen;

#ifdef VISUAL_C
// to deal with stupid vc++ stl warnings.
#pragma warning(disable:4786)
#endif 


GEN_OPTIONS::GEN_OPTIONS()
{
	m_choose_new_seed_for_every_circuit	= 0L;
	m_should_generate_k_reg_graph	= false;
	m_kreg_n						= 0;
	m_kreg_e						= 0;
	m_kreg_i						= 0;
	m_kreg_o						= 0;
	m_kreg_d						= 0;
	m_random_seed					= 0;
	m_verbose						= false;
	m_quiet							= false;
	m_sanity_checks_on				= true;
	m_no_warnings					= false;
	m_draw_graph 					= false;
	m_output_vhdl					= false;
	m_is_fast_gen					= false;

	// for delay structure creation
	m_alpha							= 1;	// cost_Comb
	m_beta							= 1;	// cost_Edge_Length
	m_gamma							= 1;	// cost_Level_Node
	m_too_many_inputs_factor		= 1;	// init. values
	m_too_few_outputs_factor		= 1;
	m_too_few_inputs_factor			= 1;
	m_mult_too_many_inputs_factor	= 1.5;
	m_mult_too_few_outputs_factor	= 1.5;
	m_mult_too_few_inputs_factor	= 1.5;
	m_delay_structure_init_temperature = 6;

	// for degree partitioning
	m_delta							= 0.5;
	m_degree_init_temperature 		= 1;

	// for final edge iteration
	m_eta							= 0.02;
	m_dff_loop_penalty_multipler 	= 10;
	m_final_edge_assignment_init_temperature = 0.0000001;
	m_final_edge_assignment_mult_double_inputs_factor = 1;

	// wirelength-approx 
	m_best_wirelength_settings		= true;
	m_wirelength_from_stat_file		= false;
	m_accept_initial_wirelength		= false;
	m_user_wirelength				= false;
	m_wirelength_wanted				= 0;

	m_result_file_name  			= "results.log";

	// unused
	m_read_level0_fanout			= false;
	m_lamda							= 0.95;  //unused
	m_phi							= 0.5;   //unused
}
GEN_OPTIONS::GEN_OPTIONS(const GEN_OPTIONS & another_options)
{
	assert(false);
	m_should_generate_k_reg_graph	= another_options.m_should_generate_k_reg_graph;
	m_kreg_n						= another_options.m_kreg_n;
	m_random_seed					= another_options.m_random_seed;
	m_choose_new_seed_for_every_circuit	= another_options.m_choose_new_seed_for_every_circuit;
	m_verbose						= another_options.m_verbose;
	m_sanity_checks_on				= another_options.m_sanity_checks_on;
	m_no_warnings					= another_options.m_no_warnings;
	m_draw_graph 					= another_options.m_draw_graph;

	m_input_file_name				= another_options.m_input_file_name;
	m_result_file_name				= another_options.m_result_file_name;

	m_best_wirelength_settings		= another_options.m_best_wirelength_settings;
	m_too_many_inputs_factor		= another_options.m_too_many_inputs_factor;
	m_too_few_outputs_factor		= another_options.m_too_few_outputs_factor;
	m_wirelength_wanted				= another_options.m_wirelength_wanted;
	m_wirelength_from_stat_file		= another_options.m_wirelength_from_stat_file;
	m_accept_initial_wirelength		= another_options.m_accept_initial_wirelength;

	m_read_level0_fanout			= false;
}

GEN_OPTIONS & GEN_OPTIONS::operator=(const GEN_OPTIONS & another_options)
{
	assert(false);
	m_should_generate_k_reg_graph	= another_options.m_should_generate_k_reg_graph;
	m_kreg_n						= another_options.m_kreg_n;
	m_random_seed					= another_options.m_random_seed;
	m_choose_new_seed_for_every_circuit	= another_options.m_choose_new_seed_for_every_circuit;
	m_verbose						= another_options.m_verbose;
	m_sanity_checks_on				= another_options.m_sanity_checks_on;
	m_no_warnings					= another_options.m_no_warnings;
	m_draw_graph 					= another_options.m_draw_graph;

	m_input_file_name				= another_options.m_input_file_name;
	m_result_file_name				= another_options.m_result_file_name;

	m_too_many_inputs_factor		= another_options.m_too_many_inputs_factor;
	m_too_few_outputs_factor		= another_options.m_too_few_outputs_factor;
	m_best_wirelength_settings		= another_options.m_best_wirelength_settings;
	m_wirelength_wanted				= another_options.m_wirelength_wanted;
	m_wirelength_from_stat_file		= another_options.m_wirelength_from_stat_file;
	m_accept_initial_wirelength		= another_options.m_accept_initial_wirelength;

	m_read_level0_fanout			= another_options.m_read_level0_fanout;

	return (*this);
}

GEN_OPTIONS::~GEN_OPTIONS()
{
}

// PRE: argc contains the number of command line arguments
//      argv contains the arguments
// POSTS: the arguments have been read in and processed
//
void  GEN_OPTIONS::process_options
(
	int argc,
	char **argv
)
{
	// randomize the random number generator.  the seed might be in which
	// case the clock is used as a seed.

	string arg,
	       next_arg;
    int argnum = 1;

	// we would have exited if we did not have at least 2 arguments
	assert(argc >= 2);

	arg = string(argv[argnum]);

	if (arg == "--help" || arg == "help" || arg == "-h" || arg == "-help")
	{
		display_option_usage();
		exit(0);
	} 
	else
	{
		m_input_file_name = string(argv[argnum]);
	}
	argnum++;

    while (argnum < argc) 
	{
		arg = string(argv[argnum]);
	
		if (arg =="--help" || arg =="help" || arg =="-h" || arg =="-help")
		{
			display_option_usage();
			exit(0);
		} 
		else if (arg =="--result") 
		{
			if (additional_arguments(argnum, argc, arg))
			{
				argnum++;
				m_result_file_name = string(argv[argnum]);
			}
		}
		else if (arg =="--verbose") 
		{
			m_verbose = true;
		} 
		else if (arg == "--quiet") 
		{
			m_quiet = true;
		} 
		else if (arg == "--nowarn") 
		{
			m_no_warnings = true;
			cout << "option: nowarn. Will not give warnings.\n";
		} 
		else if (arg == "--draw") 
		{
			m_draw_graph = true;
		} 
		else if (arg == "--fast") 
		{
			m_is_fast_gen = true;
		} 
		else if (arg == "--seed") 
		{
			if (additional_arguments(argnum, argc, arg))
			{
				argnum++;
				next_arg = string(argv[argnum]);
				m_random_seed = atol(next_arg.c_str());
				cout << "Option seed = " << m_random_seed << endl;
			}
		} 
		else if (arg == "--vhdl_output") 
		{
			m_output_vhdl = true;
		} 
		else if (arg == "--alpha") 
		{
			if (additional_arguments(argnum, argc, arg))
			{
				argnum++;
				next_arg = string(argv[argnum]);
				m_alpha = atof(next_arg.c_str());
				cout << "Option alpha = " << m_alpha << endl;
			}
		} 
		else if (arg == "--beta") 
		{
			if (additional_arguments(argnum, argc, arg))
			{
				argnum++;
				next_arg = string(argv[argnum]);
				m_beta = atof(next_arg.c_str());
				cout << "Option beta = " << m_beta << endl;
			}
		} 
		else if (arg == "--gamma") 
		{
			if (additional_arguments(argnum, argc, arg))
			{
				argnum++;
				next_arg = string(argv[argnum]);
				m_gamma = atof(next_arg.c_str());
				cout << "Option gamma = " << m_gamma << endl;
			}
		} 
		else if (arg == "--init_too_many_inputs_factor") 
		{
			if (additional_arguments(argnum, argc, arg))
			{
				argnum++;
				next_arg = string(argv[argnum]);
				m_too_many_inputs_factor = atof(next_arg.c_str());
				cout << "Option too_many_inputs = " << m_too_many_inputs_factor << endl;
			}
		} 
		else if (arg == "--init_too_few_outputs_factor") 
		{
			if (additional_arguments(argnum, argc, arg))
			{
				argnum++;
				next_arg = string(argv[argnum]);
				m_too_few_outputs_factor = atof(next_arg.c_str());
				cout << "Option too_few_outputs_factor = " << m_too_few_outputs_factor << endl;
			}
		} 
		else if (arg == "--init_too_few_inputs_factor") 
		{
			if (additional_arguments(argnum, argc, arg))
			{
				argnum++;
				next_arg = string(argv[argnum]);
				m_too_few_inputs_factor = atof(next_arg.c_str());
				cout << "Option too_few_inputs_factor = " << m_too_few_inputs_factor << endl;
			}
		} 
		else if (arg == "--mult_too_many_inputs_factor") 
		{
			if (additional_arguments(argnum, argc, arg))
			{
				argnum++;
				next_arg = string(argv[argnum]);
				m_mult_too_many_inputs_factor = atof(next_arg.c_str());
				cout << "Option too_many_inputs = " << m_mult_too_many_inputs_factor << endl;
			}
		} 
		else if (arg == "--mult_too_few_outputs_factor") 
		{
			if (additional_arguments(argnum, argc, arg))
			{
				argnum++;
				next_arg = string(argv[argnum]);
				m_mult_too_few_outputs_factor = atof(next_arg.c_str());
				cout << "Option too_few_outputs_factor = " << m_mult_too_few_outputs_factor << endl;
			}
		} 
		else if (arg == "--mult_too_few_inputs_factor") 
		{
			if (additional_arguments(argnum, argc, arg))
			{
				argnum++;
				next_arg = string(argv[argnum]);
				m_mult_too_few_inputs_factor = atof(next_arg.c_str());
				cout << "Option too_few_inputs_factor = " << m_mult_too_few_inputs_factor << endl;
			}
		} 
		else if (arg == "--delay_structure_init_temperature") 
		{
			if (additional_arguments(argnum, argc, arg))
			{
				argnum++;
				next_arg = string(argv[argnum]);
				m_delay_structure_init_temperature = atof(next_arg.c_str());
				cout << "Option delay_structure_init_temperature = " << m_gamma << endl;
			}
		} 
		else if (arg == "--delta") 
		{
			if (additional_arguments(argnum, argc, arg))
			{
				argnum++;
				next_arg = string(argv[argnum]);
				m_delta = atof(next_arg.c_str());
				cout << "Option delta = " << m_delta << endl;
			}
		} 
		else if (arg == "--degree_init_temperature") 
		{
			if (additional_arguments(argnum, argc, arg))
			{
				argnum++;
				next_arg = string(argv[argnum]);
				m_degree_init_temperature = atof(next_arg.c_str());
				cout << "Option delay_init_temperature = " << m_gamma << endl;
			}
		} 
		else if (arg == "--eta") 
		{
			if (additional_arguments(argnum, argc, arg))
			{
				argnum++;
				next_arg = string(argv[argnum]);
				m_eta = atof(next_arg.c_str());
				cout << "Option eta = " << m_eta << endl;
			}
		} 
		else if (arg == "--wirelength")
		{
			if (additional_arguments(argnum, argc, arg))
			{
				argnum++;
				next_arg = string(argv[argnum]);

				if (next_arg == "init")
				{
					cout << "Option accept intial solution wirelength " << endl;
					m_accept_initial_wirelength = true;
					m_best_wirelength_settings = false;
				}
				else if (next_arg == "stats") 
				{
					cout << "Option get the wirelength from the stats file" << endl;
					m_wirelength_from_stat_file = true;
					m_best_wirelength_settings = false;
				} 
				else
				{
					m_wirelength_wanted = static_cast<double>(atof(next_arg.c_str()));
					cout << "Option wirelength = " << m_wirelength_wanted << endl;
					m_user_wirelength			= true;
					m_best_wirelength_settings  = false;
				}
			}
		} 
		else if (arg == "--final_edge_assign_dff_loop_penalty_multiplier") 
		{
			if (additional_arguments(argnum, argc, arg))
			{
				argnum++;
				next_arg = string(argv[argnum]);
				m_dff_loop_penalty_multipler = atof(next_arg.c_str());
				cout << "Option final_edge_assign_dff_loop_penalty_multiplier = " 
						<< m_dff_loop_penalty_multipler << endl;
			}
		} 
		else if (arg == "--final_edge_assignment_init_temperature") 
		{
			if (additional_arguments(argnum, argc, arg))
			{
				argnum++;
				next_arg = string(argv[argnum]);
				m_final_edge_assignment_init_temperature = atof(next_arg.c_str());
				cout << "Option final_edge_assignment_init_temperature = " << m_final_edge_assignment_init_temperature << endl;
			}
		} 
		else
		{
	    	cerr << "Warning:  unknown option '" << arg << "' ignored." << endl;
		}
		
		argnum += 1;
	}
}
// POST: option usage has been displayed to the user
void GEN_OPTIONS::display_option_usage()
{
	cerr << "\nUsage:  cgen stats_file [Options ...]\n\n\n";
	cout << "See the external documentation for detailed" << endl;
	cout << "description of options." << endl << endl;

	cout << "General Options:\n";
	cout << "        [--help] \n";
	cout << "        [--nowarn]\n";
	cout << "        [--fast]\n";
	cout << endl;
	cout << "Delay Structure cost multipliers:\n";
	cout << "        [--alpha <float> ]     (alpha * cost_Comb)\n";
	cout << "        [--beta <float>  ]     (beta  * cost_Edge_Length)\n";
	cout << "        [--gamma <float> ]     (gamma * cost_Level_Node)\n";
	cout << "        [--init_too_many_inputs_factor <float>]\n";
	cout << "        [--init_double_input_factor <float> ]\n";
	cout << "        [--init_too_few_inputs_factor <float> ]\n" ;
	cout << "        [--init_too_few_outputs_factor <float> ]\n" ;
	cout << "        [--mult_too_many_inputs_factor <float> ]\n";
	cout << "        [--mult_double_input_factor <float> ]\n";
	cout << "        [--mult_too_few_inputs_factor <float> ]\n" ;
	cout << "        [--mult_too_few_outputs_factor <float> ]\n" ;
	cout << endl;
	cout << "Degree Partitioning cost trade-off:\n";
	cout << "        [--delta <float>]\n";
	cout << endl;
	cout << "Final Edge assignment cost trade-off:\n";
	cout << "        [--eta <float>]\n";
	cout << endl;
	cout << "Final Edge assignment Wirelength-approx:\n";
	cout <<	"        [--wirelength init | stats | <float> ]\n";
	cout << "        [--final_edge_assign_dff_loop_penalty_multiplier <float> ]\n";
	cout << endl;
	cout << "Output a dot drawning of the clone:\n";
	cout << "        [--draw]\n";
	cout << endl;
		

}


// POST: we have warned the user if there is no additional arguments 
// RETURN: true if there are additional arguments, false otherwise
bool GEN_OPTIONS::additional_arguments
(
	const int& argnum, 
	const int& argc,
	const string& arg
) const
{
	assert(argnum < argc);
	// if we don't have any additional arguments
	// print a warning
	bool additional_options = (argnum + 1 != argc);

	if (! additional_options)
	{
		cerr << "Warning:  no argument for --result. Ignoring." << endl;
	}	

	return additional_options;
}
