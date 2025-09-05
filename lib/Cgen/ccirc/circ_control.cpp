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




#include "circ_control.h"
#include <fstream>
#include "graph_medic.h"
#include "cycle_breaker.h"
#include "delay_leveler.h"
#include "statistic_reporter.h"
#include "drawer.h"
#include "wirelength_character.h"
#include "node_partitioner.h"

extern		OPTIONS * g_options;
extern 		int yyparse();
extern 		FILE *yyin;

CIRCUIT *	g_parsed_graph;

CIRC_CONTROL::CIRC_CONTROL()
{
	m_circuit = 0;
	m_input_file = 0;
}

CIRC_CONTROL::CIRC_CONTROL(const CIRC_CONTROL & another_circ_control)
{
	m_circuit		=	another_circ_control.m_circuit; 
	m_input_file	= 	another_circ_control.m_input_file;
}
	
CIRC_CONTROL & CIRC_CONTROL::operator=(const CIRC_CONTROL & another_circ_control)
{
	m_circuit		=	another_circ_control.m_circuit; 
	m_input_file	= 	another_circ_control.m_input_file;

	return (*this);
}

CIRC_CONTROL::~CIRC_CONTROL()
{
	delete m_circuit;
}
void CIRC_CONTROL::help()
{
	cerr << "\nUsage:  ccirc circuit.blif  [Options ...]\n\n";
	cerr << "\n\nFor list of options type: ccirc --help\n" << endl;
}

// 
// Read in the circuits from the blif file
// 
// PRE: nothing
// POST: we have read in the graph or failed and exited
//
void CIRC_CONTROL::read_circuits()
{
	open_circuit_input_file();
	yyin = m_input_file;

	yyparse();

	close_circuit_input_file();

	// Check to see if we have a parsed the circuit and have a graph
	if (g_parsed_graph)
	{
		// yes we do
		m_circuit = g_parsed_graph;
		g_parsed_graph = 0;
	}
	else
	{
		Fail("Could not parse the graph");
	}


	Log("Finished reading in the circuits");
	Verbose("Circuit Stats: size=" << m_circuit->get_size() <<
			"\tn= " << m_circuit->get_nNodes() <<
			"\te= " << m_circuit->get_nEdges() << "\tnumber dff = " << 
			m_circuit->get_nDFF());
	debugSep;
}


// 
// This function controls the execution of ccirc
// 
// PRE: circuit is valid
// POST: ccirc is done
//
void CIRC_CONTROL::analyze_graphs()
{
	assert(m_circuit);

	NUM_ELEMENTS size;
	bool should_log = true;

	GRAPH_MEDIC medic(m_circuit);
	CYCLE_BREAKER cycle_breaker;
	DELAY_LEVELER delay_leveler;
	NODE_PARTITIONER node_partitioner;
	STATISTIC_REPORTER statistic_reporter;
	DRAWER drawer;
	WIRELENGTH_CHARACTER wirelength_characterizer;

	size = m_circuit->get_size();
	should_log = DEBUG || (size>1000);

	Logif(should_log,"Status: Looking to break combinational cycles");
	cycle_breaker.break_cycles(m_circuit);
	
	Logif(should_log,"Status: Calculating combinational delay");
	delay_leveler.calculate_and_label_combinational_delay_levels(m_circuit);

	Logif(should_log,"Status: Doing sanity checks on the graph");
	medic.check_sanity();

	// Logif(should_log,"Status: Partitioning");
	// node_partitioner.partition_circuit(m_circuit);

	Logif(should_log, "Status: Calculating degree (fanin/fanout) stats");
	m_circuit->calculate_degree_information();

	// Logif(should_log, "Status: Final sanity check");
	// m_circuit->final_sanity_check();

	if (g_options->is_determine_wirelength_approx())
	{
		Logif(should_log, "Status: Determining lowest wirelength_approx");
		wirelength_characterizer.get_circuit_wirelength_approx(m_circuit);
	}

	Logif(should_log, "Status: Analysis is complete");

	Logif(should_log,"Status: Reporting Statistics");
	statistic_reporter.report_stats(m_circuit);

	if (g_options->is_draw_circuit())
	{
		Logif(should_log,"Status: Drawing Circuits");
		drawer.draw_full_graph(m_circuit);
	}

	Logif(should_log,"Status: Done");
}



// 
// This function tries to open the the input file
// PRE: nothing
// POST: m_input_file is open else we have exited the program
//
void CIRC_CONTROL::open_circuit_input_file()
{
	string file_name = g_options->get_input_file_name();

	if (file_name.empty()) 
	{
		exit(1);
	}
	else
	{
		debug("Working with filename " << file_name <<  " and k = " << g_options->get_k());

		m_input_file = try_to_open_file(file_name);

		if (! m_input_file)
		{
			Fail("Cannot locate the specified input file.  Please try again.");
		}
	}
}

// 
// This function tries to open the the input file trying 3 locations
// 1. The current directory
// 2. The TEST directory in the current directory if it exists
// 3. The CIRCDIR directory where CIRCDIR is an env variable that points
//    to the directory where circuits are located
//  
//  PRE: file_name is not empty 
//  POST: we have opened a circuit file if we could
//  RETURN: a valid file pointer if we could open the file else NULL
//
FILE * CIRC_CONTROL::try_to_open_file(const string & file_name)
{
	FILE * input_file;

	// try the current directory
	if ( (input_file = try_to_open_file_in_a_directory("", file_name)))
	{
		return input_file;
	}

	// try the TEST directory
	if ( (input_file = try_to_open_file_in_a_directory("TEST", file_name)) )
	{
		return input_file;
	}

	// try the CIRCDIR directory.
	// the CIRCDIR directory is an environmental variable that points to the 
	// location of your circuit files
	const char * circuit_directory = getenv("CIRCDIR");
	
	if ( (input_file = try_to_open_file_in_a_directory(circuit_directory, file_name)) )
	{
		return input_file;
	}

	// the input file was not found
	return 0;

}

// This function tries to open the specified file in the specified directory
// it will also try to add .blif as an extension to see if the user forgot to 
// add that suffix
// PRE: file_name is not empty
// RETURN: a valid file pointer if we could open the file else NULL
//
FILE * CIRC_CONTROL::try_to_open_file_in_a_directory
(
	const char * directory,
	const string & file_name
)
{
	FILE * input_file;
	string full_input_file_name;
	string directory_name = "";

	// stlport and visual c++ don't like to initialize strings with null pointers
	if (directory) 
	{
		directory_name = string(directory); 
	}


	if (! directory_name.empty())
	{
		directory_name = directory_name + "/";
	}

	// try to open the file without any further extension
	full_input_file_name = directory_name + file_name;

	if ( (input_file = open_file(full_input_file_name)) )
	{
		return input_file;
	}

	// try to open the file with a .blif extension
	full_input_file_name = directory_name + file_name + string(".blif");

	if ( (input_file = open_file(full_input_file_name)) )
	{
		return input_file;
	}

	// the input file was not opened
	
	return 0;
}


// Try to actually open the file and returns a file pointer
//
// PRE: full_file_name contains the name of the file to open
// RETURNS: file pointer if the file opened or NULL if it could not be opened
//
FILE * CIRC_CONTROL::open_file
(
	const string & full_file_name
)
{
	debug("circ_control: trying to open file " << full_file_name);
	return ((FILE *) fopen(full_file_name.c_str(), "r"));
}

// PRE: m_input file is open
// POST: m_input_file is closed
//
void CIRC_CONTROL::close_circuit_input_file()
{
	fclose(m_input_file);
}
