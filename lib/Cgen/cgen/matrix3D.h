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



#ifndef matrix_3D_H
#define matrix_3D_H

#include <vector>
#include <iostream>
#include "types.h"

class MATRIX_3D;

//
// Class_name MATRIX_3D
//
// Description
//
//	A matrix class 
//

class MATRIX_3D 
{
public:
	MATRIX_3D();
	MATRIX_3D(const NUM_ELEMENTS & nX, const NUM_ELEMENTS & nY, const NUM_ELEMENTS & nZ);
	~MATRIX_3D();
	MATRIX_3D(const MATRIX_3D& another_matrix);
	MATRIX_3D& operator= (const MATRIX_3D& another_matrix);

	void clear();
	void resize(const NUM_ELEMENTS & nX, const NUM_ELEMENTS & nY, const NUM_ELEMENTS & nZ);
	
    // Access methods to get the (i,j) element:
	NUM_ELEMENTS & operator() (const NUM_ELEMENTS & x, const NUM_ELEMENTS & y, const NUM_ELEMENTS & z);
	bool operator==(const MATRIX_3D& another_matrix);
	MATRIX_3D operator-(const MATRIX_3D& another_matrix);

	NUM_ELEMENTS get_nX() const { return m_nX; }
	NUM_ELEMENTS get_nY() const { return m_nY; }
	NUM_ELEMENTS get_nZ() const { return m_nZ; }
	long		get_sum() const;
	long		get_absolute_sum() const;
	double 		get_squared_sum() const;

	void output();

	bool is_positive() const; // has all posible entries
	bool is_zero() const;	  // has all zero entries

	NUM_ELEMENTS_VECTOR &	get_data() { return m_data; }

private:
	NUM_ELEMENTS_VECTOR 	m_data;
	NUM_ELEMENTS 		m_nX;
	NUM_ELEMENTS 		m_nY;
	NUM_ELEMENTS 		m_nZ;
};

#endif
