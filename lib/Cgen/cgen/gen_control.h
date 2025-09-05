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



#ifndef gen_control_H
#define gen_control_H

//
// Class_name GEN_CONTROL
//
// Description
//		Main control function for the gen program
//

#include "gen.h"

class GEN_CONTROL
{
public:
	GEN_CONTROL();
	GEN_CONTROL(const GEN_CONTROL & another_gen_control);
	GEN_CONTROL & operator=(const GEN_CONTROL & another_gen_control);
	~GEN_CONTROL();
	void generate_clone_circuit();
	void help();
	void open_error_and_log_files();
private:
};


#endif
