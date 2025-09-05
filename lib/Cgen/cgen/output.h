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




#ifndef OUT_H
#define OUT_H

// Warnings
//

#define Warning(warning_output) \
		if (! g_options->is_no_warn() ) { \
			cerr << "Warning: " << warning_output << endl; \
			cout << "Warning: " << warning_output << endl; \
		} 

// Error and Fail statements 
//

#define Error(error_text) cerr << "Error: " << error_text << endl;
	
#define Fail(fail_text) { \
			cerr << endl << "Terminal Error: " << fail_text << endl; \
			exit(-1);	}

#define Verbose(verbose_output) \
		if (g_options->is_verbose() && (! g_options->is_quiet())) \
			cout << verbose_output << endl;
			

// Log 
//
#define Log(log_text) cout << log_text << endl;

#define Logif(should_output,log_text) \
		if (should_output) cout << log_text << endl;

// Debug statements
//

#define dbg(debug_output) \
		if (DEBUG) \
			cout << debug_output;

#define dbgif(should_output, debug_output) \
if (DEBUG && should_output) \
	cout << debug_output;


#define debug(debug_output)	\
		if (DEBUG) \
			cout << debug_output << endl;

#define debugif(should_output, debug_output) \
		if (should_output) \
			cout << debug_output << endl;

// Separators

const string sep = "-----------------------------------------------------------------------\n";
const string Sep = "=======================================================================\n";

#define debugsep \
		if (DEBUG) \
		cout << sep;

#define debugSep \
		if (DEBUG) \
		cout << Sep;


#endif
