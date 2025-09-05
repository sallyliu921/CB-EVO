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



#ifndef lut_H
#define lut_H


#ifdef VISUAL_C
// to deal with stupid vc++ stl warnings.
#pragma warning(disable:4786)
#endif



#include <list>
#include <string>
#include <vector>
using namespace std;


//
// Class_name LUT
//
// Description
//
//	The contents of the LUT
//

typedef short VALUE_TYPE;
typedef vector<string> CUBE_VECTOR;

class LUT
{
public:
	LUT();
	LUT(const LUT & another_cube);
	LUT & operator=(const LUT & another_cube);
	~LUT();
	void add_cube(const string & new_cube, VALUE_TYPE output_value);
	bool is_sum_of_products() const;
private:
	// all the implicants that evaluate to the 	
	// output value 
	// i.e. if the output value is 1 all the cubes are in sum of product form
	// and  if the output value is 0 all the cubes are in product of sum form
	CUBE_VECTOR		m_cubes;
	VALUE_TYPE 		m_output_value;	

};


#endif
