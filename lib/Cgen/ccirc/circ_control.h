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



#ifndef circ_control_H
#define circ_control_H

#include "circ.h"
#include "circuit.h"
#include <cstdio>

//
// Class_name CIRC_CONTROl
//
// Description
//		Main control function for the circ program
//


class CIRC_CONTROL
{
public:
	CIRC_CONTROL();
	CIRC_CONTROL(const CIRC_CONTROL & another_circ_control);
	CIRC_CONTROL & operator=(const CIRC_CONTROL & another_circ_control);
	~CIRC_CONTROL();
	void read_circuits();
	void analyze_graphs();

	void print_report_on_circuits();
	void help();

private:
	CIRCUIT	*			m_circuit; 
	FILE * 				m_input_file;

	void open_circuit_input_file();
	FILE * try_to_open_file(const string & file_name);
	FILE * try_to_open_file_in_a_directory(const char * directory, 
											const string & file_name);
	FILE * open_file(const string & full_file_name);

	void close_circuit_input_file();
};


#endif
