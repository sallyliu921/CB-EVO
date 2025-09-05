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



#ifndef shape_analyzer_H
#define shape_analyzer_H

 
#include "circ.h"
#include "sequential_level.h"



//
// Class_name SHAPE_ANALYZER
//
// Description
//


class SHAPE_ANALYZER
{
public:
	enum SHAPE_TYPE {STRICT_DECREASING, CONICAL, CONICAL_WITH_BUMP, OTHER};
	enum DIRECTION_TYPE {NONE, INCREASING, DECREASING};
	SHAPE_ANALYZER();
	SHAPE_ANALYZER(const SHAPE_ANALYZER & another_shape_analyzer);
	SHAPE_ANALYZER & operator=(const SHAPE_ANALYZER & another_shape_analyzer);
	~SHAPE_ANALYZER();

	string get_shape_type_string(const DISTRIBUTION& shape);
	SHAPE_TYPE get_shape_type(const DISTRIBUTION& shape);
	DISTRIBUTION get_differential_shape(const DISTRIBUTION& shape);
	NUM_ELEMENTS get_number_peaks ( const DISTRIBUTION& shape);
	DISTRIBUTION get_low_pass_filter_of_shape(const DISTRIBUTION & shape);
private:
	
	//SHAPE_ANALYZER::DIRECTION get_direction( const NUM_ELEMENTS& present_value, const NUM_ELEMENTS& last_value);

	SHAPE_ANALYZER::SHAPE_TYPE determine_shape_type(const NUM_ELEMENTS& number_peaks,
					const SHAPE_ANALYZER::DIRECTION_TYPE& starting_direction);

	bool get_next_peak(NUM_ELEMENTS& position, NUM_ELEMENTS& lowest_value_left, NUM_ELEMENTS& lowest_value_right,
						const DISTRIBUTION& diff_shape);

	SHAPE_ANALYZER::DIRECTION_TYPE get_starting_direction(const DISTRIBUTION & shape);
};


#endif
