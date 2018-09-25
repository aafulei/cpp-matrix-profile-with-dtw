#ifndef PNT_HPP
#define PNT_HPP

#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>

// TO DO : allow INF print
template<typename T>
void aprint(T const* a, int n, int k = 0)
{
	int w = 1;
	if (n > 1)
		w = log10(n - 1) + 1;
	for (int i = 0; i < n; ++i)
		std::cout << "[" << std::setw(w) << (i + k) << "]" << " " << a[i] << (i % 5 == 4 ? '\n' : '\t');
	std::cout << std::endl;
}

template<typename T>
void mprint(T const*const* m, int nrow, int ncol)
{
	for (int i = 0; i < nrow; ++i)
		aprint(m[i], ncol, i * ncol);
}

void pprint(int i, int n, bool stay = false);

#endif