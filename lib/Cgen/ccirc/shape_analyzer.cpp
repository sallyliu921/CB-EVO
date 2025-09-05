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



#include "shape_analyzer.h"

 

SHAPE_ANALYZER::SHAPE_ANALYZER()
{
}
SHAPE_ANALYZER::SHAPE_ANALYZER(const SHAPE_ANALYZER & another_shape_analyzer)
{
}

SHAPE_ANALYZER & SHAPE_ANALYZER::operator=(const SHAPE_ANALYZER & another_shape_analyzer)
{

	return (*this);
}

SHAPE_ANALYZER::~SHAPE_ANALYZER()
{
}

string SHAPE_ANALYZER::get_shape_type_string(const DISTRIBUTION& shape)
{
	SHAPE_TYPE shape_type = get_shape_type(shape);
	string shape_string;
	enum SHAPE_TYPE {STRICT_DECREASING, CONICAL, CONICAL_WITH_BUMP, OTHER};

	if (shape_type == SHAPE_ANALYZER::STRICT_DECREASING)
	{
		shape_string = "STRICT_DECREASING";
	}
	else if (shape_type == SHAPE_ANALYZER::CONICAL)
	{
		shape_string = "CONICAL";
	}
	else if (shape_type == SHAPE_ANALYZER::CONICAL_WITH_BUMP)
	{
		shape_string = "CONICAL_WITH_BUMP";
	}
	else 
	{
		assert(shape_type == SHAPE_ANALYZER::OTHER);
		shape_string = "OTHER";
	}

	return shape_string;
}


DISTRIBUTION SHAPE_ANALYZER::get_differential_shape
(
	const DISTRIBUTION& shape
)
{
	assert(shape.size() > 2);
	
	DISTRIBUTION differential;

	NUM_ELEMENTS shape_size = static_cast<NUM_ELEMENTS>(shape.size());
	NUM_ELEMENTS shape_index = 0;
	NUM_ELEMENTS difference = 0;

	for (shape_index = 1; shape_index < shape_size; shape_index++)
	{
		difference = shape[shape_index] - shape[shape_index-1];
		differential.push_back(difference);
	}
	return differential;
}


bool SHAPE_ANALYZER::get_next_peak
(
	NUM_ELEMENTS& position,
	NUM_ELEMENTS& lowest_value_left,
	NUM_ELEMENTS& lowest_value_right,
	const DISTRIBUTION& diff_shape
)
{
	NUM_ELEMENTS diff_size = static_cast<NUM_ELEMENTS>(diff_shape.size());

	while (position < diff_size && diff_shape[position] <= 0)
	{
		position++;
	}
	lowest_value_left = position;

	// while we haven't hit the end of the array and haven't hit a peak value
	while (position < diff_size && diff_shape[position] >= 0)
	{
		position++;
	}

	if (position == diff_size)
	{
		// we didn't find a peak
		return false;
	}
	else
	{
		lowest_value_right = position;
		while (lowest_value_right < diff_size-1 && diff_shape[lowest_value_right] <= 0)
		{
			lowest_value_right++;
		}

		return true;
	}
}
		


NUM_ELEMENTS SHAPE_ANALYZER::get_number_peaks
(
	const DISTRIBUTION& shape
)
{
	assert(shape.size() > 2);
	
	DISTRIBUTION differential = get_differential_shape(shape);

	NUM_ELEMENTS diff_size = static_cast<NUM_ELEMENTS>(differential.size());
	//NUM_ELEMENTS difference = 0;
	NUM_ELEMENTS position = 0;
	NUM_ELEMENTS lowest_value_left = 0;
	NUM_ELEMENTS lowest_value_right = 0;
	NUM_ELEMENTS number_peaks = 0;
	bool peak_exists;


	while (position < diff_size)
	{
		peak_exists = get_next_peak(position, lowest_value_left, lowest_value_right, differential);

		if (peak_exists)
		{
			number_peaks++;
		}
		//shape[position] - AVG(shape[lowest_value_left] + shape[lowest_right])
	}


	return number_peaks;
}


SHAPE_ANALYZER::SHAPE_TYPE SHAPE_ANALYZER::get_shape_type
(
	const DISTRIBUTION& shape
)
{
	SHAPE_TYPE shape_type;

	NUM_ELEMENTS number_peaks = get_number_peaks(shape);
	DIRECTION_TYPE starting_direction = get_starting_direction(shape);

	shape_type = determine_shape_type(number_peaks, starting_direction);

	return shape_type;
}

SHAPE_ANALYZER::SHAPE_TYPE SHAPE_ANALYZER::determine_shape_type
(
	const NUM_ELEMENTS& number_peaks,
	const SHAPE_ANALYZER::DIRECTION_TYPE& starting_direction
)
{
	SHAPE_TYPE shape_type = SHAPE_ANALYZER::CONICAL;

	if (number_peaks == 0)
	{
		if (starting_direction == SHAPE_ANALYZER::DECREASING)
		{
			shape_type = SHAPE_ANALYZER::STRICT_DECREASING;
		}
		else
		{
			shape_type = SHAPE_ANALYZER::OTHER;
		}
	}
	else if (number_peaks == 1)
	{
		if (starting_direction == SHAPE_ANALYZER::INCREASING)
		{
			shape_type = SHAPE_ANALYZER::CONICAL;
		}
		else
		{
			shape_type = SHAPE_ANALYZER::OTHER;
		}
	}
	else if (number_peaks == 2)
	{
		if (starting_direction == SHAPE_ANALYZER::INCREASING)
		{
			shape_type = SHAPE_ANALYZER::CONICAL_WITH_BUMP;
		}
		else
		{
			shape_type = SHAPE_ANALYZER::OTHER;
		}
	}
	else
	{
		shape_type = SHAPE_ANALYZER::OTHER;
	}
	return shape_type;
}
		

SHAPE_ANALYZER::DIRECTION_TYPE SHAPE_ANALYZER::get_starting_direction
(
	const DISTRIBUTION & shape
)
{

	DISTRIBUTION differential = get_differential_shape(shape);
	NUM_ELEMENTS diff_size = static_cast<NUM_ELEMENTS>(differential.size());
	NUM_ELEMENTS diff_position = 0;

	while (diff_position < diff_size && differential[diff_position] == 0)
	{
		diff_position++;
	}

	if (diff_position == diff_size)
	{
		return SHAPE_ANALYZER::NONE;
	}
	else
	{
		if (differential[diff_position] > 0)
		{
			return SHAPE_ANALYZER::INCREASING;
		}
		else
		{
			return SHAPE_ANALYZER::DECREASING;
		}
	}
}

DISTRIBUTION SHAPE_ANALYZER::get_low_pass_filter_of_shape
(
	const DISTRIBUTION & shape
)
{
	DISTRIBUTION low_pass_filtered_shape = shape;
	assert(low_pass_filtered_shape.size() == shape.size());
	NUM_ELEMENTS max_size = static_cast<NUM_ELEMENTS>(shape.size() - 1);
	NUM_ELEMENTS shape_index = 1;

	for (shape_index = 1; shape_index <= max_size-1; shape_index++)
	{
		low_pass_filtered_shape[shape_index] = 
			(shape[shape_index-1] + 4*shape[shape_index] + shape[shape_index+1])/6;
	}

	low_pass_filtered_shape[0] = (4*shape[0] + shape[1])/5;
	low_pass_filtered_shape[max_size] = (shape[max_size-1] + 4*shape[max_size])/5;

	return low_pass_filtered_shape;
}

