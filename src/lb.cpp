#include "lb.hpp"

double lower_bound_keogh(double const* q, double const* c, int m, int r, bool debug)
{
	if (m <= 0) {
		std::cerr << "Lower Bound Keogh Error : Bad Length m = " << m << std::endl;
		return BAD_LOWER_BOUND;
	}
	if (r < 0 || r >= m) {
		std::cerr << "Lower Bound Keogh Error : Bad Band r = " << r << " with Length m = " << m << std::endl;
		return BAD_LOWER_BOUND;
	}
	
	double * u = new double[m];
	double * l = new double[m];
	double sum = 0;
	double const* b;
	double const* e;

	for (int i = 0; i < m; ++i) {
		b = q + std::max(i - r, 0);
		e = q + std::min(i + r + 1, m);
		u[i] = *std::max_element(b, e);
		l[i] = *std::min_element(b, e);
		if (c[i] > u[i])
			sum += (c[i] - u[i]) * (c[i] - u[i]);
		else if (c[i] < l[i])
			sum += (l[i] - c[i]) * (l[i] - c[i]);
	}

	if (debug) {
		std::cout << "[Q]" << std::endl; aprint(q, m);
		std::cout << "[U]" << std::endl; aprint(u, m);
		std::cout << "[L]" << std::endl; aprint(l, m);
		std::cout << "[C]" << std::endl; aprint(c, m);
	}

	delete[] u;
	delete[] l;

	return sqrt(sum);
}

double lower_bound_front(double const* q, double const* c, int r, bool debug)
{
	double * uo = new double[r + 1];
	double * lo = new double[r + 1];
	double * un = new double[r + 1];
	double * ln = new double[r + 1];
	double so = 0;
	double sn = 0;
	double const* e = q + r + 1;

	for (int i = 0; i <= r; ++i, ++e) {
		uo[i] = *std::max_element(q, e);
		lo[i] = *std::min_element(q, e);
		if (c[i] > uo[i])
			so += (c[i] - uo[i]) * (c[i] - uo[i]);
		else if (c[i] < lo[i])
			so += (lo[i] - c[i]) * (lo[i] - c[i]);

		if (i == 0)
			un[i] = ln[i] = INFINITY;
		else {
			un[i] = *std::max_element(q + 1, e);
			ln[i] = *std::min_element(q + 1, e);
			if (c[i] > un[i]) 
				sn += (c[i] - un[i]) * (c[i] - un[i]);
			else if (c[i] < ln[i]) 
				sn += (ln[i] - c[i]) * (ln[i] - c[i]);
		}
	}

	if (debug) {
		std::cout << "[Q]  "; aprint(q,  r + 1);
		std::cout << "[UO] "; aprint(uo, r + 1);
		std::cout << "[LO] "; aprint(lo, r + 1);
		std::cout << "[C]  "; aprint(c,  r + 1);
		std::cout << "[UN] "; aprint(un, r + 1);
		std::cout << "[LN] "; aprint(ln, r + 1);
		std::cout << "[C]  "; aprint(c,  r + 1);
	}

	delete[] uo;
	delete[] lo;
	delete[] un;
	delete[] ln;

	return sn - so;
}

double lower_bound_back(double const* q, double const* c, int m, int r, bool debug)
{
	double * uo = new double[r + 1];
	double * lo = new double[r + 1];
	double * un = new double[r + 1];
	double * ln = new double[r + 1];
	double so = 0;
	double sn = 0;
	double const* b = q + m - 2 * r;
	double const*const e = q + m;

	for (int i = 0, k = m - r; i <= r; ++i, ++k, ++b) {
		if (i == r)
			uo[i] = lo[i] = INFINITY;
		else {
			uo[i] = *std::max_element(b, e);
			lo[i] = *std::min_element(b, e);
			if (c[k] > uo[i])
				so += (c[k] - uo[i]) * (c[k] - uo[i]);
			else if (c[k] < lo[i])
				so += (lo[i] - c[k]) * (lo[i] - c[k]);
		}
		
		un[i] = *std::max_element(b, e + 1);
		ln[i] = *std::min_element(b, e + 1);
		if (c[k] > un[i]) 
			sn += (c[k] - un[i]) * (c[k] - un[i]);
		else if (c[k] < ln[i]) 
			sn += (ln[i] - c[k]) * (ln[i] - c[k]);
	}

	if (debug) {
		std::cout << "[Q]  "; aprint(q + m - r, r + 1);
		std::cout << "[UO] "; aprint(uo, r + 1);
		std::cout << "[LO] "; aprint(lo, r + 1);
		std::cout << "[C]  "; aprint(c + m - r, r + 1);
		std::cout << "[UN] "; aprint(un, r + 1);
		std::cout << "[LN] "; aprint(ln, r + 1);
		std::cout << "[C]  "; aprint(c + m - r, r + 1);
	}

	delete[] uo;
	delete[] lo;
	delete[] un;
	delete[] ln;

	return sn - so;
}

// assert r <= 0.5 m
double lower_bound_update(double lbk, double const* q, double const* c, int m, int r, bool debug)
{
	if (m <= 0) {
		std::cerr << "Lower Bound Update Error : Bad Length m = " << m << std::endl;
		return BAD_LOWER_BOUND;
	}
	if (r < 0 || m < 2 * r) {
		std::cerr << "Lower Bound Update Error : Bad Band r = " << r << " with Length m = " << m << std::endl;
		return BAD_LOWER_BOUND;
	}

	double lbf = lower_bound_front(q, c, r, false);
	double lbb = lower_bound_back(q, c, m, r, false);
	double lbu = lbk * lbk + lbf + lbb;

	if (debug) {
		printf("LBK = %lf; LBF = %lf; LBB = %lf\n", lbk, lbf, lbb);
		if (lbu < 0)
			std::cerr << "Lower Bound Update Warning : " << lbu << " is Taken as 0" << std::endl;
	}

	if (lbu < 0)
		lbu = 0;

	return sqrt(lbu);
}

// #include <cstdio>
// #include <ctime>
// #include <random>

// #include "mat.hpp"
// #include "pnt.hpp"

// int main()
// {
// 	int n = 1000, m = n / 10, l = n - m + 1, r = 2;
// 	double * t = new double[n];
// 	std::default_random_engine e(0);
// 	std::uniform_real_distribution<double> u(-100, 100);
// 	for (int i = 0; i < n; ++i)
// 		t[i] = u(e);
// 	double ** lbk = new_matrix<double>(l, l);
// 	double ** lbu = new_matrix<double>(l, l);

// 	clock_t tb, te;
// 	printf("::::Compute Lower Bound Keogh Directly:::::::::\n");
// 	tb = clock();
// 	for (int i = 0; i < l; ++i) {
// 		for (int j = 0; j < l; ++j) {
// 			lbk[i][j] = lower_bound_keogh(t + i, t + j, m, r);
// 		}
// 		pprint(i, l);
// 	}
// 	te = clock();
// 	printf("Time = %.1lf s\n", (double)(te - tb) / CLOCKS_PER_SEC);
	
// 	tb = clock();
// 	for (int i = 0; i < l; ++i)
// 		lbu[i][i] = 0;
// 	for (int j = 1; j < l; ++j) {
// 		lbu[0][j] = lower_bound_keogh(t, t + j, m, r);
// 		int j0 = j;
// 		for (int i = 0; j < l - 1; ++i, ++j)
// 			lbu[i+1][j+1] = lower_bound_update(lbu[i][j], t + i, t + j, m, r);
// 		j = j0;
// 	}
// 	for (int i = 1; i < l; ++i) {
// 		lbu[i][0] = lower_bound_keogh(t + i, t, m, r);
// 		int i0 = i;
// 		for (int j = 0; i < l - 1; ++i, ++j) {
// 			lbu[i+1][j+1] = lower_bound_update(lbu[i][j], t + i, t + j, m, r);
// 		}
// 		i = i0;
// 	}
// 	te = clock();
// 	printf("::::Compute Lower Bound Keogh Incrementally::::\nTime = %.1lf s\n", (double)(te - tb) / CLOCKS_PER_SEC);

// 	printf("::::Consistency Check Begins:::::::::::::::::::\n");
// 	int safe = 10;
// 	for (int i = 0; i < l; ++i)
// 		for (int j = 0; j < l; ++j)
// 			if (std::abs(lbk[i][j] - lbu[i][j]) > 1e-4 && safe-- > 0)
// 				printf("Error : [%d,%d]\tLBK = %lf\tLBU = %lf\n", i, j, lbk[i][j], lbu[i][j]);
// 	printf("::::Consistency Check Ends:::::::::::::::::::::\n");

// 	del_matrix(lbk, l);
// 	del_matrix(lbu, l);

// 	return 0;
// }