#include "md.hpp"

int read(std::string fname, double * t, int n)
{
	if (n <= 0) {
		std::cerr << "Read Error : Length n = " << n << std::endl;
		return BAD_READ;
	}

	std::ifstream ifile(fname, std::ifstream::in);
	if (!ifile.good()) {
		std::cerr << "Read Error : File Name = " << fname << std::endl;
		return BAD_READ;
	}

	double point;
	int i = 0;
	while (ifile >> point && i < n) {
		t[i] = point;
		i += 1;
	}
	return i;
}

bool parse(int argc, char * argv[], double * t, int & n, int & m, int & l, int & f, int & g, int & r, std::default_random_engine & e)
{
	if (argc == 1) {
		n = DEF_LEN;
		t = (double *)realloc(t, n * sizeof(double));
		std::uniform_real_distribution<double> u(-100, 100);
		for (int i = 0; i < n; ++i)
			t[i] = u(e);
	}
	else if (argc <= 5) {
		int nread = read(argv[1], t, n);
		if (nread == BAD_READ)
			return false;
		if (nread != n) {
			n = nread;
			t = (double *)realloc(t, n * sizeof(double));
		}
		if (argc == 3) {
			r = atoi(argv[2]);
			if (r < 0)
				return false;
		}
		else if (argc == 4 || argc == 5) {
			if (atoi(argv[2]) < n) {
				n = atoi(argv[2]);
				if (n <= 0)
					return false;
				t = (double *)realloc(t, n * sizeof(double));
			}
			r = atoi(argv[3]);
			if (r < 0)
				return false;
			if (argc == 5) {
				double mm = atof(argv[4]);
				if (mm < 1.0)
					m = mm * n;
				else
					m = (int)mm;
				if (m <= 0)
					return false;
			}
		}
	}
	else
		return false;
	if (n < MIN_LEN)
		return false;
	if (m == NOT_SET)
		m = n / 10;
	else if (m > n / 2)
		return false;
	l = n - m + 1;
	if (r == NOT_SET)
		r = DEF_BAND;
	else if (r > m / 2)
		return false;
	f = m / 4 + 1;
	g = (l - f) * (l - f + 1);
	return true;
}

void usage()
{
	std::cerr << "Command Line Options" << std::endl;
	std::cerr << "# 1 : [command]" << std::endl;
	std::cerr << "# 2 : [command] [filename]" << std::endl;
	std::cerr << "# 3 : [command] [filename] [band width r]" << std::endl;
	std::cerr << "# 4 : [command] [filename] [length n] [band width r]" << std::endl;
	std::cerr << "# 5 : [command] [filename] [length n] [band width r] [subsequence length m]" << std::endl;
	std::cerr << "# 6 : [command] [filename] [length n] [band width r] [subsequence proportion mm]" << std::endl;
	std::cerr << "Require Time Series Data Length n > " << MIN_LEN << std::endl;
	std::cerr << "Require Sakoe-Chiba Band r in [0, m/2]" << std::endl;
	std::cerr << "Require Subsequence Length m in [1, n/2], or Subsequence Proportion mm in (0, 0.5]" << std::endl;
	std::cerr << "Max Time Series Data Length n = " << MAX_LEN  << std::endl;
	std::cerr << "Default Time Series Data Length n = " << DEF_LEN  << std::endl;
	std::cerr << "Default Sakoe-Chiba Band r = " << DEF_BAND << std::endl;
	std::cerr << "Default Subsequence Length m = n/10, or Subsequence Proportion mm = 0.1" << std::endl;
	std::cerr << "If No Data File (#1), Program will Generate a Random Time Series." << std::endl;
}

double ** mp_dtw_brutal_force(double const* t, int l, int m, int f, int r, double * mpv, int * mpi, bool make_return)
{
	double ** dtw = new_matrix<double>(l, l, INFINITY);
	clock_t tb, te;

	tb = clock();
	for (int i = 0; i < l; ++i) {
		for (int j = i + f; j < l; ++j) {
			dtw[i][j] = dynamic_time_warping(t + i, t + j, m, r);
			dtw[j][i] = dtw[i][j];
		}
		pprint(i, l);
	}
	for (int j = 0; j < l; ++j) {
		double min_dtw = INFINITY;
		int min_dtw_row = NO_ROW;
		for (int i = 0; i < l; ++i) {
			if (std::abs(i - j) < f)
				continue;
			else if (dtw[i][j] < min_dtw) {
				min_dtw = dtw[i][j];
				min_dtw_row = i;
			}
		}
		mpv[j] = min_dtw;
		mpi[j] = min_dtw_row;
	}
	te = clock();
	printf("DTW Brutal Force : Time = %.1lf s\n", (double)(te - tb) / CLOCKS_PER_SEC);

	if (make_return)
		return dtw;
	else
		del_matrix(dtw, l);
	return nullptr;
}

double ** lower_bound_direct(double const* t, int l, int m, int f, int r)
{
	double ** lb = new_matrix<double>(l, l, INFINITY);
	clock_t tb, te;
	tb = clock();
	for (int i = 0; i < l; ++i) {
		for (int j = i + f; j < l; ++j) {
			lb[i][j] = lower_bound_keogh(t + i, t + j, m, r);
			lb[j][i] = lower_bound_keogh(t + j, t + i, m, r);
			lb[i][j] = lb[j][i] = std::max(lb[i][j], lb[j][i]);
		}
		pprint(i, l);
	}	
	te = clock();
	printf("Lower Bound Direct : Time = %.1lf s\n", (double)(te - tb) / CLOCKS_PER_SEC);
	return lb;
}

double ** lower_bound_incremental(double const* t, int l, int m, int f, int r)
{
	double ** lb = new_matrix<double>(l, l, INFINITY);
	clock_t tb, te;
	tb = clock();
	for (int j = f; j < l; ++j) {
		lb[0][j] = lower_bound_keogh(t, t + j, m, r);
		int j0 = j;
		for (int i = 0; j < l - 1; ++i, ++j) {
			lb[i+1][j+1] = lower_bound_update(lb[i][j], t + i, t + j, m, r);
		}
		j = j0;
		pprint(j, 2 * l);
	}
	for (int i = f; i < l; ++i) {
		lb[i][0] = lower_bound_keogh(t + i, t, m, r);
		int i0 = i;
		for (int j = 0; i < l - 1; ++i, ++j) {
			lb[i+1][j+1] = lower_bound_update(lb[i][j], t + i, t + j, m, r);
		}
		i = i0;
		pprint(i + l, 2 * l);
	}
	for (int i = 0; i < l; ++i)
		for (int j = i + f; j < l; ++j)
			lb[i][j] = lb[j][i] = std::max(lb[i][j], lb[j][i]);
	te = clock();
	printf("Lower Bound Incremental : Time = %.1lf s\n", (double)(te - tb) / CLOCKS_PER_SEC);
	return lb;
}

void mp_dtw_lower_bound(double const* t, double const*const* lb, int l, int m, int f, int g, int r, double * mpv, int * mpi)
{
	clock_t tb, te;
	int save = 0;

	tb = clock();
	for (int j = 0; j < l; ++j) {
		double dtw, min_dtw = INFINITY;
		int min_dtw_row = NO_ROW;
		for (int i = 0; i < l; ++i) {
			if (std::abs(i - j) < f)
				continue;
			else if (lb[i][j] < min_dtw) {
				dtw = dynamic_time_warping(t + i, t + j, m, r);
				if (dtw < min_dtw) {
					min_dtw = dtw;
					min_dtw_row = i;
				}
			}
			else
				save += 1;
		}
		mpv[j] = min_dtw;
		mpi[j] = min_dtw_row;
		pprint(j, l);
	}
	te = clock();
	printf("DTW Lower Bound : Time = %.1lf s\n", (double)(te - tb) / CLOCKS_PER_SEC);
	printf("# Saved = %d (%.0lf%%)\n", save, 100.0 * save / g);
}

void mp_dtw_lower_bound_randomized(double const* t, double const*const* lb, int l, int m, int f, int g, int r, double * mpv, int * mpi, std::list<int> const& rows, std::default_random_engine & e)
{
	clock_t tb, te;
	int save = 0;
	double lower_bound, dtw, min_dtw = INFINITY;
	int min_dtw_row = NO_ROW;
	std::list<int>::iterator it, itb, ite;

	tb = clock();
	for (int j = 0; j < l; ++j) {		
		std::list<int> a(rows);
		it = a.begin();
		itb = std::next(it, std::max(j - f + 1, 0));
		ite = std::next(it, std::min(j + f, l));
		a.erase(itb, ite);
		min_dtw = INFINITY;
		min_dtw_row = NO_ROW;
		while (!a.empty()) {
			it = a.begin();
			std::uniform_int_distribution<unsigned> random_row(0, a.size() - 1);	// C++11 guarantees O(1) time for std::list.size()
			std::advance(it, random_row(e));
			lower_bound = lb[*it][j];
			if (lower_bound > min_dtw) {
				save += 1;
				it = a.erase(it);
				if (it == a.end())
					--it;
			}
			else {
				dtw = dynamic_time_warping(t + *it, t + j, m, r);
				if (dtw < min_dtw) {
					min_dtw = dtw;
					min_dtw_row = *it;
				}
				it = a.erase(it);
				if (it == a.end())
					--it;
			}
		}
		mpv[j] = min_dtw;
		mpi[j] = min_dtw_row;
		pprint(j, l);
	}
	te = clock();
	printf("DTW Randomized : Time  = %.1lf s\n", (double)(te - tb) / CLOCKS_PER_SEC);
	printf("# Saved = %d (%.0lf%%)\n", save, 100.0 * save / g);
}

void mp_dtw_lower_bound_simulated_annealing(double const* t, double const*const* lb, int l, int m, int f, int g, int r, double * mpv, int * mpi, std::list<int> const& rows, std::default_random_engine & e)
{
	clock_t tb, te;
	int save = 0;
	double lower_bound, prev_lower_bound, dtw, min_dtw = INFINITY;
	int min_dtw_row = NO_ROW;
	std::list<int>::iterator it, it0, itb, ite;
	std::bernoulli_distribution random_sign;
	bool sign;
	unsigned steps;
	bool accept;
	double prob;
	const double T = 100;
	
	tb = clock();
	for (int j = 0; j < l; ++j) {
		std::list<int> a(rows);
		it = a.begin();
		itb = std::next(it, std::max(j - f + 1, 0));
		ite = std::next(it, std::min(j + f, l));
		a.erase(itb, ite);
		it = a.begin();
		prev_lower_bound = lb[*it][j];
		min_dtw = dynamic_time_warping(t + *it, t + j, m, r);
		min_dtw_row = *it;
		while (!a.empty()) {
			it0 = it;
			sign = random_sign(e);
			std::uniform_int_distribution<unsigned> random_steps(1u, std::max((unsigned)a.size() / 4, 1u));
			steps = random_steps(e);
			while (steps-- != 0 && (sign ? it != prev(a.end()) : it != a.begin())) {
				if (sign)
					++it;
				else
					--it;
			}
			lower_bound = lb[*it][j];
			if (lower_bound > min_dtw) {
				save += 1;
				prev_lower_bound = lower_bound;
				it = a.erase(it);
				if (it == a.end())
					--it;
			}
			else {
				if (lower_bound < prev_lower_bound)
					accept = true;
				else {
					prob = exp((prev_lower_bound - lower_bound) / (a.size() * T));
					std::bernoulli_distribution random_accept(prob);
					accept = random_accept(e);
				}
				if (accept) {
					dtw = dynamic_time_warping(t + *it, t + j, m, r);
					if (dtw < min_dtw) {
						min_dtw = dtw;
						min_dtw_row = *it;
					}
					prev_lower_bound = lower_bound;
					it = a.erase(it);
					if (it == a.end())
						--it;
				}
				else
					it = it0;
			}
		}
		mpv[j] = min_dtw;
		mpi[j] = min_dtw_row;
		pprint(j, l);
	}
	te = clock();
	printf("DTW Simulated Annealing : Time  = %.1lf s\n", (double)(te - tb) / CLOCKS_PER_SEC);
	printf("# Saved = %d (%.0lf%%)\n", save, 100.0 * save / g);
}

void lower_bound_consistency_check(double ** lb1, double ** lb2, int l, double error, int safe)
{
	printf("::::Lower Bound Consistency Check Begins:::::::\n");
	for (int i = 0; i < l; ++i)
		for (int j = 0; j < l; ++j)
			if (safe-- > 0 && std::abs(lb1[i][j] - lb2[i][j]) > error)
				printf("Error : [%d,%d] LB(1) = %lf != LB(2) = %lf\n", i, j, lb1[i][j], lb2[i][j]);
	printf("::::Lower Bound Consistency Check Ends:::::::::\n");
}

void lower_boundedness_check(double ** lb, double ** dtw, int l, int f, double error, int safe)
{
	printf("::::Lower Boundedness Check Begins:::::::::::::\n");
	for (int i = 0; i < l; ++i)
		for (int j = i + f; j < l; ++j)
			if (safe-- > 0 && lb[i][j] >= (dtw[i][j] + error))
				printf("Error : [%d,%d] LB = %lf > DTW = %lf\n", i, j, lb[i][j], dtw[i][j]);
	printf("::::Lower Boundedness Check Ends:::::::::::::::\n");
}

void matrix_profile_consistency_check(double * mpv1, int * mpi1, double * mpv2, int * mpi2, int l, double error, int safe)
{
	printf("::::Matrix Profile Consistency Check Begins::::\n");
	for (int j = 0; j < l; ++j)
		if (safe-- > 0 && (std::abs(mpv1[j] - mpv2[j]) > error || mpi1[j] != mpi2[j]))
			printf("Error : Col %d : \tMP(1) = %.1lf with MPI(1) = %d differs from MP(2) = %.1lf with MPI(2) = %d\n", j, mpv1[j], mpi1[j], mpv2[j], mpi2[j]);
	printf("::::Matrix Profile Consistency Check Ends::::::\n");
}

int main(int argc, char * argv[])
{
	int n, m, l, f, g, r; 
	n = MAX_LEN;
	m = l = f = g = r = NOT_SET;
	double * t = (double *)calloc(n, sizeof(double));
	std::default_random_engine e(0);

	if (!parse(argc, argv, t, n, m, l, f, g, r, e)) {
		std::cerr << "Main Error : Program Terminates." << std::endl;
		usage();
		return EXIT_FAILURE;
	}
	
	printf("Length of Time Series Data n = %d\n", n);
	printf("Length of Subsequence m = %d\n", m);
	printf("Length of Matrix Profile l = %d\n", l);
	printf("Length of Forbidden Zone f = %d\n", f);
	printf("Entries to Compute g = %d\n", g);
	printf("Sakoe-Chiba Band r = %d\n", r);

	unsigned const error = 1e-4;
	double const safe = 10;
	int const seed = 0;
	std::list<int> rows;
	for (int i = 0; i < l; ++i)
		rows.push_back(i);

	double *  mpv_dbf = new double[l]; std::fill_n(mpv_dbf, l, INFINITY);
	int    *  mpi_dbf = new int   [l]; std::fill_n(mpi_dbf, l, NO_ROW);
	double ** dtw = mp_dtw_brutal_force(t, l, m, f, r, mpv_dbf, mpi_dbf, true);

	double ** lbd = lower_bound_direct(t, l, m, f, r);
	double ** lbi = lower_bound_incremental(t, l, m, f, r);	

	lower_bound_consistency_check(lbd, lbi, l, error, safe);
	lower_boundedness_check(lbi, dtw, l, f, error, safe);

	double *  mpv_dlb = new double[l]; std::fill_n(mpv_dlb, l, INFINITY);
	int    *  mpi_dlb = new int   [l]; std::fill_n(mpi_dlb, l, NO_ROW);
	mp_dtw_lower_bound(t, lbi, l, m, f, g, r, mpv_dlb, mpi_dlb);

	matrix_profile_consistency_check(mpv_dbf, mpi_dbf, mpv_dlb, mpi_dlb, l, error, safe);

	double *  mpv_dlr = new double[l]; std::fill_n(mpv_dlr, l, INFINITY);
	int    *  mpi_dlr = new int   [l]; std::fill_n(mpi_dlr, l, NO_ROW);
	mp_dtw_lower_bound_randomized(t, lbi, l, m, f, g, r, mpv_dlr, mpi_dlr, rows, e);

	double *  mpv_dls = new double[l]; std::fill_n(mpv_dls, l, INFINITY);
	int    *  mpi_dls = new int   [l]; std::fill_n(mpi_dls, l, NO_ROW);
	mp_dtw_lower_bound_simulated_annealing(t, lbi, l, m, f, g, r, mpv_dls, mpi_dls, rows, e);

	matrix_profile_consistency_check(mpv_dlr, mpi_dlr, mpv_dls, mpi_dls, l, error, safe);
	
	free(t);

	return EXIT_SUCCESS;
}