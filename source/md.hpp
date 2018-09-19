#ifndef MD_HPP
#define MD_HPP

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <list>
#include <random>
#include <string>

#include "dtw.hpp"
#include "lb.hpp"
#include "mat.hpp"
#include "pnt.hpp"

#define BAD_READ	-1
#define DEF_LEN		1000
#define DEF_BAND	1
#define MAX_LEN		10000
#define MIN_LEN		100
#define NO_ROW		-1
#define NOT_SET		-1


template<typename T>
int awrite(T const* a, int n, std::string fname = "out") {
    std::ofstream ofile(fname, std::ofstream::out);
    for (int i = 0; i < n; ++i)
    	ofile << a[i] << ' ';
    ofile.flush();
    ofile.close();
    return 0;
}

int read(std::string fname, double * t, int n);
bool parse(int argc, char * argv[], double * t, int & n, int & m, int & l, int & f, int & g, int & r, std::default_random_engine & e);
void usage();
double ** mp_dtw_brutal_force(double const* t, int l, int m, int f, int r, double * mpv, int * mpi, bool make_return);
double ** lower_bound_direct(double const* t, int l, int m, int f, int r);
double ** lower_bound_incremental(double const* t, int l, int m, int f, int r);
void mp_dtw_lower_bound(double const* t, double const*const* lb, int l, int m, int f, int g, int r, double * mpv, int * mpi);
void mp_dtw_lower_bound_randomized(double const* t, double const*const* lb, int l, int m, int f, int g, int r, double * mpv, int * mpi, std::list<int> const& rows, std::default_random_engine & e);
void mp_dtw_lower_bound_simulated_annealing(double const* t, double const*const* lb, int l, int m, int f, int g, int r, double * mpv, int * mpi, std::list<int> const& rows, std::default_random_engine & e);
void lower_bound_consistency_check(double ** lb1, double ** lb2, int l, double error, int safe);
void lower_boundedness_check(double ** lb, double ** dtw, int l, int f, double error, int safe);
void matrix_profile_consistency_check(double * mpv1, int * mpi1, double * mpv2, int * mpi2, int l, double error, int safe);


#endif