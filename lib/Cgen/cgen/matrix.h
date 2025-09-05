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



#ifndef matrix_H
#define matrix_H

#include <vector>
#include <iostream>
using namespace std;

typedef long INDEX_SIZE;

class MATRIX;

typedef vector<MATRIX *> MATRICES;
typedef vector<INDEX_SIZE> MATRIX_TYPE;

//
// Class_name MATRIX
//
// Description
//
//	A matrix class 
//

class MATRIX 
{
public:
	MATRIX();
	MATRIX(const INDEX_SIZE & nRows, const INDEX_SIZE & nColumns);
	~MATRIX();
	MATRIX(const MATRIX& another_matrix);
	MATRIX& operator= (const MATRIX& another_matrix);

	MATRIX_TYPE & get_data() { return m_data; }
	void resize(const INDEX_SIZE & nRows, const INDEX_SIZE & nColumns);
	void clear();
	
    // Access methods to get the (i,j) element:
	INDEX_SIZE & operator() (const INDEX_SIZE & row, const INDEX_SIZE & col);
	bool operator==(const MATRIX& another_matrix) const;
	MATRIX operator-(const MATRIX& another_matrix) const;
	friend inline ostream& operator<<(ostream& str, const MATRIX& matrix);
	INDEX_SIZE get_value(const INDEX_SIZE & row, const INDEX_SIZE & col) const;

	void zero_negative_entries();


	INDEX_SIZE get_nRows() const { return m_nRows; }
	INDEX_SIZE get_nColumns() const { return m_nColumns; }
	long		get_sum() const;
	long		get_absolute_sum() const;
	double 		get_squared_sum() const;


	MATRIX_TYPE	get_row(const INDEX_SIZE & row) const;

	bool is_positive() const; // has all posible entries
	bool is_zero() const;	  // has all zero entries
private:
	MATRIX_TYPE 	m_data;
	INDEX_SIZE 		m_nRows;
	INDEX_SIZE 		m_nColumns;
};

inline ostream& operator<<(ostream& stream, const MATRIX& matrix)
{
	for (INDEX_SIZE row = 0; row < matrix.get_nRows(); row++)
	{
		for (INDEX_SIZE col = 0; col < matrix.get_nColumns(); col++)
		{
			stream << matrix.get_value(row, col) << " ";
		}
		stream << "\n";
	}
	stream << endl;

	return stream;
	
}

#endif
