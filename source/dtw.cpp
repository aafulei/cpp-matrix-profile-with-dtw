#include "dtw.hpp"

// TO DO : static variable to count steps to control print format
void hprint(int const*const* h, int i, int j)
{
	if (i < 0 || j < 0) {
		std::cerr << "History Print Error : Bad Index i = " << i << ", j = " << j << std::endl;
		return;
	}
	if (i == 0 && j == 0)
		printf("(0,0)");
	else {
		if (h[i][j] == FROM_LOWER_LEFT)
			hprint(h, i-1, j-1);
		else if (h[i][j] == FROM_LOWER)
			hprint(h, i-1, j);
		else if (h[i][j] == FROM_LEFT)
			hprint(h, i, j-1);
		printf("-->(%d, %d)", i, j);
	}
}

// by default m = 1
double squared_euclidean_distance(double const* a, double const* b, int m)
{
	if (m <= 0) {
		std::cerr << "Squared Euclidean Distance Error : Bad Length m = " << m << std::endl;
		return BAD_DISTANCE;
	}
	if (m == 1)
		return (a[0] - b[0]) * (a[0] - b[0]);
	double dif, sum = 0;
	for (int i = 0; i < m; ++i) {
		dif = a[i] - b[i];
		sum += dif * dif;
	}
	return sum;
}

// by default m = 1
double euclidean_distance(double const* a, double const* b, int m)
{
	if (m <= 0) {
		std::cerr << "Euclidean Distance Error : Bad Length m = " << m << std::endl;
		return BAD_DISTANCE;
	}
	if (m == 1)
		return std::abs(a[0] - b[0]);
	return sqrt(squared_euclidean_distance(a, b, m));
}

// TO DO : redesign DTW to allow very large data (rotate matrix by 45 degrees)
// when r == m - 1, the band is trivially no restriction
double dynamic_time_warping(double const* q, double const* c, int m, int r, bool ddebug, bool hdebug, bool well_defined)
{
	if (m <= 0) {
		std::cerr << "Dynamic Time Warping Error : Bad Length m = " << m << std::endl;
		return BAD_DISTANCE;
	}
	if (r < 0 || r >= m) {
		std::cerr << "Dynamic Time Warping Error : Bad Band r = " << r << " with Length m = " << m << std::endl;
		return BAD_DISTANCE;
	}
	if (r == 0)
		return euclidean_distance(q, c, m);
	double ** d = well_defined ? new_matrix<double>(m, m, INFINITY) : new_matrix<double>(m, m);
	int ** h =  well_defined ? new_matrix<int>(m, m, FROM_NOWHERE) : new_matrix<int>(m, m);
	d[0][0] = squared_euclidean_distance(q, c);
	for (int i = 1; i <= r; ++i) {
		d[i][0] = squared_euclidean_distance(q + i, c) + d[i-1][0];
		h[i][0] = FROM_LOWER;
	}
	if (r != m - 1)
		d[r+1][0] = INFINITY;
	for (int j = 1; j <= r; ++j) {
		d[0][j] = squared_euclidean_distance(q, c + j) + d[0][j-1];
		h[0][j] = FROM_LEFT;
	}
	if (r != m - 1)
		d[0][r+1] = INFINITY;
	for (int i = 1; i < m; ++i) {
		int b = std::max(i - r, 1);
		int e = std::min(i + r + 1, m);
		if (b != 1)
			d[i][b-1] = INFINITY;
		for (int j = b; j < e ; ++j) {
			d[i][j] = squared_euclidean_distance(q + i, c + j);
			std::array<double, 3> a = { d[i-1][j-1], d[i-1][j], d[i][j-1] };
			auto mit = std::min_element(a.begin(), a.end());
			d[i][j] += *mit;
			h[i][j] = mit - a.begin();
		}
		if (e != m)
			d[i][e] = INFINITY;
	}
	if (ddebug)
		mprint(d, m, m);
	if (hdebug) {
		hprint(h, m-1, m-1);
		std::cout << std::endl;
	}
	double dtw = sqrt(d[m-1][m-1]);
	del_matrix(d, m);
	del_matrix(h, m);
	return dtw;
}

double dynamic_time_warping(double const* q, double const* c, int m, bool ddebug, bool hdebug, bool well_defined)
{
	return dynamic_time_warping(q, c, m, m - 1, ddebug, hdebug, well_defined);
}

// #include <ctime>
// #include <random>

// int main()
// {
// 	int m = 1e4;
// 	double * q = new double[m];
// 	double * c = new double[m];
// 	std::default_random_engine e(0);
// 	std::uniform_real_distribution<double> u(-100, 100);
// 	for (int i = 0; i < m; ++i) {
// 		q[i] = u(e);
// 		c[i] = u(e);
// 	}
// 	std::cout <<"ED : "<< euclidean_distance(q, c, m) << std::endl;
// 	for (int r = 0; r < 100; ++r)
// 		std::cout <<"DTW (r = " << r << "): "<< dynamic_time_warping(q, c, m, r) << std::endl;
// 	return 0;
// }