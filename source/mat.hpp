#ifndef MAT_HPP
#define MAT_HPP

#include <cstddef>

template<typename T> T ** new_matrix(int nrow, int ncol, T val);
template<typename T> T ** new_matrix(int nrow, int ncol);
template<typename T> void del_matrix(T ** mat, int nrow);

template<typename T> T ** new_matrix(int nrow, int ncol, T val)
{
	T ** mat = new T * [nrow];
	for (int i = 0; i < nrow; ++i) {
		mat[i] = new T[ncol];
		for (int j = 0; j < ncol; ++j)
			mat[i][j] = val;
	}
	return mat;
}

template<typename T> T ** new_matrix(int nrow, int ncol)
{
	T ** mat = new T * [nrow];
	for (int i = 0; i < nrow; ++i)
		mat[i] = new T[ncol];
	return mat;
}

template<typename T> void del_matrix(T ** mat, int nrow)
{
	for (int i = 0; i < nrow; ++i)
		delete[] mat[i];
	delete[] mat;
}

#endif