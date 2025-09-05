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


#ifndef CURRENT_DATE
#define CURRENT_DATE "<unknown>"
#endif



#include "circ.h"
#include "options.h"
#include "circ_control.h"
#include "circ_version.h"
#include <iostream>

OPTIONS * g_options;

int main(int argc, char ** argv)
{
	cout << "\n\n";
	cout << "CCIRC Circuit Characterization Software Version " << circ_version();
	cout << endl;
	cout << "By Paul Kundarewich, Mike Hutton, and Jonathan Rose\n";
	cout << "This code is licensed only for non-commercial use.\n" << endl;

	CIRC_CONTROL circ_control;

    if (argc == 1) { circ_control.help(); exit(0); }

	g_options = new OPTIONS;
	assert(g_options);

	g_options->process_options(argc, argv);

	debug("Reading in the circuits");
	circ_control.read_circuits();

	debug("Analyzing the circuits"); debugSep;
	circ_control.analyze_graphs();

	return 0;
}
