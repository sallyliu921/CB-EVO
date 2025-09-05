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



#include "matrix3D.h"
#include <assert.h>

MATRIX_3D::MATRIX_3D()
{
	m_nX 	= 0;
	m_nY	= 0;
	m_nZ	= 0;
}
MATRIX_3D::MATRIX_3D(const NUM_ELEMENTS & nX, const NUM_ELEMENTS & nY, const NUM_ELEMENTS & nZ)
{
	assert(nX != 0 && nY != 0 && nZ != 0);
	m_nX 	= nX;
	m_nY	= nY;
	m_nZ	= nZ;

	m_data = NUM_ELEMENTS_VECTOR(m_nX*m_nY*m_nZ, 0);

}
// resize the matrix 
//
// PRE: nX nY, nZ > 0
// POST: m_data size is now m_nX*m_nY*m_n
//       m_nX, m_nY, m_nZ have been updated
//
void MATRIX_3D::resize(const NUM_ELEMENTS & nX, const NUM_ELEMENTS & nY, const NUM_ELEMENTS & nZ)
{
	assert(nX >= 0 && nY >= 0 && nZ >= 0);
	m_nX 	= nX;
	m_nY	= nY;
	m_nZ	= nZ;

	m_data.resize(m_nX*m_nY*m_nZ, 0);
}


MATRIX_3D::~MATRIX_3D()
{
}

MATRIX_3D::MATRIX_3D(const MATRIX_3D & another_matrix)
{
	m_data 	= another_matrix.m_data;
	m_nX  	= another_matrix.m_nX;
	m_nY 	= another_matrix.m_nY;
	m_nZ 	= another_matrix.m_nZ;
}


MATRIX_3D & MATRIX_3D::operator=(const MATRIX_3D & another_matrix)
{
	m_data	= another_matrix.m_data;
	m_nY 	= another_matrix.m_nY;
	m_nX  	= another_matrix.m_nX;
	m_nZ 	= another_matrix.m_nZ;

	return (*this);
}

// gets matrix[row,col]
//
// PRE: x, y, z are valid for this matrix
// RETURN: a reference to the value at x,y,z 
//
NUM_ELEMENTS & MATRIX_3D::operator() 
(
	const NUM_ELEMENTS & x, const NUM_ELEMENTS & y, const NUM_ELEMENTS & z
) 
{
	NUM_ELEMENTS coordinate = 0;
	assert(x < m_nX && y < m_nY && z < m_nZ);
	assert(x >= 0 && y >= 0 && z >= 0);
	coordinate = x*(m_nX*m_nZ) + y*m_nZ + z;
	assert(coordinate >= 0 && static_cast<unsigned>(coordinate) < m_data.size());
	return m_data[coordinate];
}

// 
//  overload the == operator to return whether the 
//  two matrices are equal
//  
//  PRE: both matricies are the same size
//  RETURNS: this matrix == another matrix
//
bool MATRIX_3D::operator==(const MATRIX_3D& another_matrix)
{
	NUM_ELEMENTS index = 0;
	NUM_ELEMENTS max_coord = m_nX*m_nY;

	for (index = 0; index < max_coord; index++)
	{
		if (m_data[index] != another_matrix.m_data[index])
		{
			return false;
		}
	}

	return true;
}


// zero the matrix
//
// PRE: nothing
// POST: the matrix has all zero entries
//
void MATRIX_3D::clear()
{
	NUM_ELEMENTS index = 0;
	NUM_ELEMENTS max_coord = m_nX*m_nY;

	for (index = 0; index < max_coord; index++)
	{
		m_data[index] = 0;
	}

}

long MATRIX_3D::get_sum() const
{
	NUM_ELEMENTS index = 0;
	NUM_ELEMENTS max_coord = m_nX*m_nY*m_nZ;
	long sum = 0;

	for (index = 0; index < max_coord; index++)
	{
		sum += m_data[index];
	}

	return sum;
}
long MATRIX_3D::get_absolute_sum() const
{
	NUM_ELEMENTS index = 0;
	NUM_ELEMENTS max_coord = m_nX*m_nY*m_nZ;
	long sum = 0;
	long value = 0;

	for (index = 0; index < max_coord; index++)
	{
		value = m_data[index];
		
		if (value < 0)
		{
			sum += -value;
		}
		else
		{
			sum += value;
		}
	}

	return sum;
}

// RETURN: sum over all coordinates of each data item squared
double MATRIX_3D::get_squared_sum() const
{

	NUM_ELEMENTS index = 0;
	NUM_ELEMENTS max_coord = m_nX*m_nY*m_nZ;
	double sum = 0;

	for (index = 0; index < max_coord; index++)
	{
		sum += m_data[index] * m_data[index];
	}

	return sum;
}

void MATRIX_3D::output()
{

	NUM_ELEMENTS index = 0;
	NUM_ELEMENTS max_coord = m_nX*m_nY*m_nZ;

	for (index = 0; index < max_coord; index++)
	{
		if (m_data[index] != 0)
		{
			if ((index % 40) == 0)
			{
				cout << m_data[index] << endl;
			}
			else
			{
				cout << m_data[index] << " ";
			}
		}
	}
}

//
// RETURNS: whether all entries in the matrix positive
//
bool MATRIX_3D::is_positive() const
{
	NUM_ELEMENTS index = 0;
	NUM_ELEMENTS max_coord = m_nX*m_nY*m_nZ;

	for (index = 0; index < max_coord; index++)
	{
		if (m_data[index] < 0)
		{
			return false;
		}
	}

	return true;
}
//
// RETURNS: whether all entries in the matrix are zero
//
bool MATRIX_3D::is_zero() const
{
	NUM_ELEMENTS index = 0;
	NUM_ELEMENTS max_coord = m_nX*m_nY*m_nZ;

	for (index = 0; index < max_coord; index++)
	{
		if (m_data[index] != 0)
		{
			return false;
		}
	}

	return true;
}
