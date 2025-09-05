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




#include "options.h"
#include "util.h"
#include "circ.h"

#define Warning_for_options


OPTIONS::OPTIONS()
{
	m_input_file_name   = "";
	m_output_file_name	= "";
	m_circuit_name		= "";

    m_k					= 6;					
	m_partitioning_type = OPTIONS::KWAY;
	m_nPartitions		= 8;


	m_verbose 			= true;
	m_quiet				= false;
	m_no_warn			= false;
	

	m_store_luts		= false;


	m_determine_wirelength_approx = false;

    m_draw 				= false;

    m_expand_luts 		= false;	

	m_display_pi_and_dff_distributions = false;
	m_display_inter_cluster_matricies_at_each_edge_length = false;
	m_display_statistics_on_delay_defining_edges = false;

	m_ubfactor			= 20;		// balancing factor. defined differently for k-way vs. bi-partitioning
}

OPTIONS::OPTIONS(const OPTIONS & another_options)
{
	m_verbose 			= another_options.m_verbose;
	m_quiet				= another_options.m_quiet;
	m_no_warn			= another_options.m_no_warn;
	m_circuit_name		= another_options.m_circuit_name;

	m_output_file_name	= another_options.m_output_file_name;

    m_k					= another_options.m_k;
	m_store_luts		= false;
	m_partitioning_type = another_options.m_partitioning_type;
	m_nPartitions		= another_options.m_nPartitions;

	m_determine_wirelength_approx = another_options.m_determine_wirelength_approx;

    m_draw 					= another_options.m_draw;

    m_expand_luts 		= another_options.m_expand_luts;

	m_display_pi_and_dff_distributions = another_options.m_display_pi_and_dff_distributions;
	m_display_statistics_on_delay_defining_edges = another_options.m_display_statistics_on_delay_defining_edges;
	m_display_inter_cluster_matricies_at_each_edge_length = 
										another_options.m_display_inter_cluster_matricies_at_each_edge_length;


}

OPTIONS & OPTIONS::operator=(const OPTIONS & another_options)
{
	m_verbose 			= another_options.m_verbose;
	m_quiet				= another_options.m_quiet;
	m_no_warn			= another_options.m_no_warn;
	m_circuit_name		= another_options.m_circuit_name;

	m_output_file_name	= another_options.m_output_file_name;

    /* processing options and information*/
    m_k					= another_options.m_k;
	m_store_luts		= false;
	m_partitioning_type = another_options.m_partitioning_type;
	m_nPartitions		= another_options.m_nPartitions;

	m_determine_wirelength_approx = another_options.m_determine_wirelength_approx;

    m_draw 					= another_options.m_draw;

    m_expand_luts 		= another_options.m_expand_luts;

	m_display_pi_and_dff_distributions = another_options.m_display_pi_and_dff_distributions;
	m_display_statistics_on_delay_defining_edges = another_options.m_display_statistics_on_delay_defining_edges;
	m_display_inter_cluster_matricies_at_each_edge_length = 
										another_options.m_display_inter_cluster_matricies_at_each_edge_length;


	return (*this);
}
OPTIONS::~OPTIONS()

{
}

// PRE: argc contains the number of command line arguments
//      argv contains the arguments
// POSTS: the arguments have been read in and processed
//        the circuit name has been obtained fromt the filename
//
void OPTIONS::process_options
(
	int argc,
	char **argv
)
{
	read_arguments(argc, argv);

	m_circuit_name = get_circuit_name_from_filename(m_input_file_name);

}
// PRE: argc contains the number of command line arguments
//      argv contains the arguments
// POSTS: the arguments have been read in and processed
//
void OPTIONS::read_arguments
(
	int argc,
	char **argv
)
{
    short argnum;
    string arg, next_arg;
	string text;

    argnum = 1;

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
		else if (arg == "--partition_type") 
		{
			if (additional_arguments(argnum, argc, arg))
			{
				argnum++;
				next_arg = string(argv[argnum]);

				if (next_arg == "bi")
				{
					cout << "option: partition type: recursive bi-partitioning\n";
					m_partitioning_type = OPTIONS::RECURSIVE_BI;
				}
				else if (next_arg == "kway")
				{
					cout << "option: partition type: k-way partitioning\n";
					m_partitioning_type = OPTIONS::KWAY;
				}
				else
				{
					cerr << "Warning: unknown partitioning type found:'" << next_arg  
						<< "'.  Ignoring. "  << endl;
				}
			}
		} 
		else if (arg == "--out") 
		{
			if (additional_arguments(argnum, argc, arg))
			{
				argnum++;
				next_arg = string(argv[argnum]);
				m_output_file_name = string(argv[argnum]);
			}
			cout << "option: output file: " << m_output_file_name << endl;
		} 
		else if (arg == "--verbose") 
		{
			m_verbose = true;
		} 
		else if (arg == "--quiet") 
		{
			m_quiet = true;
		} 
		else if (arg =="--nowarn") 
		{
			m_no_warn = true;
			cout << "option: nowarn. Will not give warnings.\n";
		} 
		else if (arg == "--k") 
		{
			cerr << "Error: This function is not yet implemented. if it were k would equal: " << m_k << endl;
			assert(false);

			if (additional_arguments(argnum, argc, arg))
			{
				argnum++;
				next_arg = string(argv[argnum]);
				m_k = static_cast<K_TYPE>(atoi(next_arg.c_str()));
			}
        } 
		else if (arg == "--partitions") 
		{
			if (additional_arguments(argnum, argc, arg))
			{
				argnum++;
				next_arg = string(argv[argnum]);
				m_nPartitions = atoi(next_arg.c_str());
				cout << "option:  nPartitions: " << m_nPartitions << endl;
			}
        } 
		else if (arg == "--ubfactor") 
		{
			if (additional_arguments(argnum, argc, arg))
			{
				argnum++;
				next_arg = string(argv[argnum]);
				m_ubfactor = atoi(next_arg.c_str());
				cout << "option:  ubfactor: " << m_ubfactor << endl;
			}
        } 
		else if (arg == "--wirelength_approx")
		{
			m_determine_wirelength_approx = true;
			cout << "option:  Determine wirelength_approx" << endl;
        } 
		else if (arg == "--draw") 
		{
	    	m_draw = true;
		} 
		else if (arg == "--display_pi_and_dff_distributions") 
		{
	    	m_display_pi_and_dff_distributions = true;
        } 
		else if (arg == "--display_inter_cluster_matricies_at_each_edge_length") 
		{
	    	m_display_inter_cluster_matricies_at_each_edge_length = true;
        } 
		else if (arg == "--display_statistics_on_delay_defining_edges") 
		{
			m_display_statistics_on_delay_defining_edges = true;
        } 
		else if (arg == "--expand_luts") 
		{
			m_expand_luts = true;
		} 
		else 
		{
	    	cerr << "Warning:  unknown option '" << arg << "' ignored." << endl;
		}

		argnum += 1;
    }
}




// POST: option usage has been displayed to the user
void OPTIONS::display_option_usage() const
{
	cout << "\nUsage:  ccirc circuit.blif [Options...]\n\n\n";
	cout << "See the external documentation for detailed" << endl;
	cout << "description of options.\n\n" << endl;

	cout << "General Options:\n";
	cout << "        [--help] \n";
	cout << "        [--nowarn]\n";
	cout << "        [--out]\n";
	cout << endl;
	cout << "Partitioning Options:\n";
	cout << "        [--partition_type  bi | kway]\n";
	cout << "        [--partitions <int>]\n";
	cout << "        [--ubfactor <int>]";
	cout << endl;
	cout << "Calculate Wirelength-approx:\n";
	cout << "        [--wirelength_approx]\n";
	cout << endl;
	cout << "Output a dot drawning of the clone:\n";
	cout << "        [--draw]\n";
	cout << endl;
}

// PRE: file_name has the file name
// RETURNS: the circuit name obtained from the filename
string OPTIONS::get_circuit_name_from_filename
(
	const string & file_name
)
{
	// strip the directory name and then remove the suffix
	string circuit_name = util_strip_directory_name(file_name);
	circuit_name = util_strip_file_extension(circuit_name);

	cout << "option:  circuit name is '" << circuit_name << "'" << endl;

	assert(circuit_name != "");
	return circuit_name;
}	



// POST: we have warned the user if there is no additional arguments 
// RETURN: true if there are additional arguments, false otherwise
bool OPTIONS::additional_arguments
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
