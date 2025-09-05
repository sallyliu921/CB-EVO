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





#include "lut.h"
#include <assert.h>



LUT::LUT()
{
	m_output_value	=	-1;
}
LUT::LUT(const LUT & another_lut)
{
	m_cubes			= another_lut.m_cubes;
	m_output_value	= another_lut.m_output_value;
}
LUT & LUT::operator=(const LUT & another_lut)
{
	m_cubes			= another_lut.m_cubes;
	m_output_value	= another_lut.m_output_value;

	return (*this);
}
LUT::~LUT()
{
}

void LUT::add_cube
(
 	const string & new_cube,
	VALUE_TYPE output_value
)
{
	m_cubes.push_back(new_cube);
}

bool LUT::is_sum_of_products() const
{
	assert(m_output_value != -1);
	
	// the output value should be defined before this call

	return (m_output_value == 1);
}

