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



#include <cstdio>
#include "gen_control.h"
#include "graph.h"
#include "spec.h"
#include "post_process.h"
#include "gen_drawer.h"
#include "final_improvement.h"

extern	GEN_OPTIONS * g_options;

GEN_CONTROL::GEN_CONTROL()
{
}

GEN_CONTROL::GEN_CONTROL(const GEN_CONTROL & another_gen_control)
{
}
	
GEN_CONTROL & GEN_CONTROL::operator=(const GEN_CONTROL & another_gen_control)
{

	return (*this);
}

GEN_CONTROL::~GEN_CONTROL()
{
}




// Read in the specs from the spec file and generate a clone circuit .
//
// The input to the generation process is the characterizations of the circuit which are:
// Comb, Latched, the number of clusters, and for each cluster the number of nodes, 
// the number of primary inputs, the number of flip-flops, the node shape, the primary output shape, 
// the latched shape, the input shape, the output shape, the fanout distribution, 
// the intra-cluster edge length distribution, and the inter-cluster input and output edge length 
// distributions. 
//
// The output from the generation process is a circuit in BLIF format ready to be placed and routed. 
// We set the function of each LUT to be a NAND gate to give each LUT a logic type. 
// The generation algorithm proceeds in four steps. 
//
// In Step 1, we create the delay structure of the circuit by assigning edges between the 
// level nodes in the circuit (recall that a level node contains all of the individual nodes 
// at each delay level). In Step 2, the fanout distribution in each cluster is partitioned among 
// the level nodes. In Step 3, the level nodes are split into individual nodes and the individual 
// nodes are prepared for final edge assignment. In Step 4, edges are assigned
// between individual nodes in the circuit based on the delay structure of the circuit and the fanouts
//assigned to the individual nodes. 
//
// PRE: nothing 
// POST: we have generated a clone circuit 

void GEN_CONTROL::generate_clone_circuit()
{
	CIRCUIT_SPEC * circuit_spec = new CIRCUIT_SPEC;
	POST_PROCESS post_processor;
	GEN_DRAWER drawer;
	SIMPLE_EDGE_ITERATOR simple_edge_iterator;

	circuit_spec->read_spec();

	CIRCUIT_GRAPH circuit(circuit_spec);

	// Step 1
	// creates the delay structure - the large scale connectivity in the circuit
	circuit.create_delay_structure();

	// Step 2
	// partition the fanout amongst the level nodes
	circuit.partition_degrees();

	// Step 3
	// split the level nodes into individual nodes and 
	// prepare the graph for final edge assignment
	circuit.split_nodes();

	// Step 4
	// assign the simple edges between the individual nodes
	circuit.assign_simple_edges();

	// iterate to improve the solution
	simple_edge_iterator.iterate(&circuit);

	// post-process the graph to remove node violations 
	// (nodes with too many inputs/not enough inputs/no outputs, etc)
	post_processor.post_process(&circuit);


	circuit.final_sanity_check();

	circuit.output_final_circuit_quality();

	if (g_options->is_output_vhdl())
	{
		circuit.output_vhdl();
	}
	else
	{
		circuit.output_blif();
	}

	if (g_options->is_draw_graph())
	{
		drawer.draw_full_graph(&circuit);
	}

	delete circuit_spec;
	circuit_spec = 0;
}


void GEN_CONTROL::help()
{
	cerr << "\nUsage:  cgen stats_file [Options ...]\n";
	cerr << "\n\nFor list of options type: cgen --help\n" << endl;
}
