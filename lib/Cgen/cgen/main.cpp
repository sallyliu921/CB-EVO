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



#include "gen.h"

#include "gen_control.h"
#include "gen_options.h"
#include "rand.h"
#include "util.h"
#include "gen_version.h"

GEN_OPTIONS * g_options;
RANDOM_NUMBER * g_rand_number_gen;

int main(int argc, char **argv)
{
	cout << "\n\n";
	cout << "CGEN Synthetic Circuit Generation Version " << gen_version();
	cout << endl;
	cout << "By Paul Kundarewich, Mike Hutton, and Jonathan Rose\n";
	cout << "This code is licensed only for non-commercial use.\n" << endl;

    int start_ticks = util_ticks();
	int end_ticks;
	long seed = 0;

	GEN_CONTROL gen_control;

    if (argc == 1) { gen_control.help(); exit(0); }

	g_options = new GEN_OPTIONS;
	g_options->process_options(argc, argv);

	seed = g_options->get_random_seed();
	g_rand_number_gen = new RANDOM_NUMBER(seed);


	gen_control.generate_clone_circuit();
 
	end_ticks = util_ticks();

	Log("\nElapsed time " << end_ticks - start_ticks << " seconds");


	delete g_rand_number_gen;
	delete g_options;

    return 0;
}

