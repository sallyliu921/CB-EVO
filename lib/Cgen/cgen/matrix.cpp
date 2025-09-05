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





#include "matrix.h"
#include <assert.h>

MATRIX::MATRIX()
{
	m_nColumns 	= 0;
	m_nRows		= 0;
}
MATRIX::MATRIX(const INDEX_SIZE & nRows, const INDEX_SIZE & nColumns)
{
	assert(nRows != 0 && nColumns != 0);
	m_nColumns 	= nColumns;
	m_nRows		= nRows;

	m_data = MATRIX_TYPE(m_nColumns*m_nRows, 0);

}

MATRIX::~MATRIX()
{
}

MATRIX::MATRIX(const MATRIX & another_matrix)
{
	m_data 		= another_matrix.m_data;
	m_nRows 	= another_matrix.m_nRows;
	m_nColumns  = another_matrix.m_nColumns;
}

// resize the matrix 
//
// PRE: nRows and nColumns > 0 
// POST: m_data size is now nRows*nColumns
//       m_nColumns and m_nRows have been updated 
void MATRIX::resize
(
	const INDEX_SIZE & nRows, 
	const INDEX_SIZE & nColumns
)
{
	assert(nRows > 0 && nColumns > 0);
	m_nColumns 	= nColumns;
	m_nRows		= nRows;

	m_data = MATRIX_TYPE(m_nColumns*m_nRows, 0);
}


MATRIX & MATRIX::operator=(const MATRIX & another_matrix)
{
	m_data 		= another_matrix.m_data;
	m_nRows 	= another_matrix.m_nRows;
	m_nColumns  = another_matrix.m_nColumns;

	return (*this);
}

//  overload the - operator to subtract one 
//  matrix from another
//
//  PRE: both matricies are the same size in terms of the 
//       number of rows and columns
//  RETURNS: this matrix - another matrix
//
MATRIX MATRIX::operator-(const MATRIX& another_matrix) const
{
	assert(m_data.size() == another_matrix.m_data.size());
	assert((m_nRows == another_matrix.m_nRows) && (m_nColumns == another_matrix.m_nColumns));

	MATRIX result(m_nRows, m_nColumns);

	INDEX_SIZE row, column;

	for (row = 0; row < m_nRows; row++)
	for (column = 0; column < m_nColumns; column++)
	{
		result(row, column) = m_data[row*m_nColumns + column] - another_matrix.m_data[row*m_nColumns + column];
	}

	return result;
}

// gets matrix[row,col]
//
// PRE: row and col are valid for this matrix
// RETURN: a reference to the value at row and col
//
INDEX_SIZE & MATRIX::operator() 
(
	const INDEX_SIZE & row, const INDEX_SIZE & col
) 
{
	assert(row < m_nRows && col < m_nColumns);
	assert(row >= 0 && col >= 0);
	return m_data[row*m_nColumns + col];
}

// gets matrix[row,col]
//
// PRE: row and col are valid for this matrix
// RETURN: the value at row and col
//
INDEX_SIZE MATRIX::get_value
(
	const INDEX_SIZE & row, const INDEX_SIZE & col
) const
{
	assert(row < m_nRows && col < m_nColumns);
	assert(row >= 0 && col >= 0);
	return m_data[row*m_nColumns + col];
}

// 
//  overload the == operator to return whether the 
//  two matrices are equal
//  
//  PRE: both matricies are the same size in terms of the 
//       number of rows and columns
//  RETURNS: this matrix == another matrix
//
bool MATRIX::operator==(const MATRIX& another_matrix) const
{
	assert(m_data.size() == another_matrix.m_data.size());
	assert((m_nRows == another_matrix.m_nRows) && (m_nColumns == another_matrix.m_nColumns));

	INDEX_SIZE index = 0;
	INDEX_SIZE max_coord = m_nColumns*m_nRows;

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
void MATRIX::clear()
{
	INDEX_SIZE index = 0;
	INDEX_SIZE max_coord = m_nColumns*m_nRows;

	for (index = 0; index < max_coord; index++)
	{
		m_data[index] = 0;
	}

}

//
// RETURNS: the sum over all entries in the matrix
//
long MATRIX::get_sum() const
{
	INDEX_SIZE index = 0;
	INDEX_SIZE max_coord = m_nColumns*m_nRows;
	long sum = 0;

	for (index = 0; index < max_coord; index++)
	{
		sum += m_data[index];
	}

	return sum;
}
//
// RETURNS: sum over all entries(absolute value(entry))
//
long MATRIX::get_absolute_sum() const
{
	INDEX_SIZE index = 0;
	INDEX_SIZE max_coord = m_nColumns*m_nRows;
	long sum = 0;
	long value = 0;

	for (index = 0; index < max_coord; index++)
	{
		value = m_data[index];
		
		if (value < 0)
		{
			value *= -1;
		}
		sum += value;
	}

	return sum;
}

// make all negative entries positive
//
// PRE: nothing
// POST: all entries in the matrix are positive
//
void MATRIX::zero_negative_entries()
{

	INDEX_SIZE index = 0;
	INDEX_SIZE max_coord = m_nColumns*m_nRows;

	for (index = 0; index < max_coord; index++)
	{
		if (m_data[index] < 0)
		{
			m_data[index] = 0;
		}
	}
}

//
// RETURNS: sum over all entries(entry*entry)
//
double MATRIX::get_squared_sum() const
{

	INDEX_SIZE index = 0;
	INDEX_SIZE max_coord = m_nColumns*m_nRows;
	double sum = 0;

	for (index = 0; index < max_coord; index++)
	{
		sum += m_data[index] * m_data[index];
	}

	return sum;
}

//
// RETURNS: whether all entries in the matrix positive
//
bool MATRIX::is_positive() const
{
	INDEX_SIZE index = 0;
	INDEX_SIZE max_coord = m_nColumns*m_nRows;

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
bool MATRIX::is_zero() const
{
	INDEX_SIZE index = 0;
	INDEX_SIZE max_coord = m_nColumns*m_nRows;

	for (index = 0; index < max_coord; index++)
	{
		if (m_data[index] != 0)
		{
			return false;
		}
	}

	return true;
}
//
// returns the specified row in the matrix
//
// RETURNS: matrix[row]
//
MATRIX_TYPE	MATRIX::get_row(const INDEX_SIZE & row) const
{
	assert(row < m_nRows);

	MATRIX_TYPE row_matrix;
	INDEX_SIZE data = 0;
	INDEX_SIZE col;

	for (col = 0; col < m_nColumns; col++)
	{
		data = m_data[row*m_nColumns + col];
		row_matrix.push_back(data);
	}

	assert(static_cast<unsigned>(m_nColumns) == row_matrix.size());

	return row_matrix;
}
