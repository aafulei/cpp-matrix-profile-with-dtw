#ifndef DTW_HPP
#define DTW_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <iostream>
#include <limits>

#include "mat.hpp"
#include "pnt.hpp"

#define BAD_DISTANCE	-1
#define FROM_NOWHERE	-1
#define FROM_LOWER_LEFT  0
#define FROM_LOWER 		 1
#define FROM_LEFT 		 2

void hprint(int const*const* h, int i, int j);
double squared_euclidean_distance(double const* a, double const* b, int m = 1);
double euclidean_distance(double const* a, double const* b, int m);
double dynamic_time_warping(double const* q, double const* c, int m, int r, bool ddebug = false, bool hdebug = false, bool well_defined = false);
double dynamic_time_warping(double const* q, double const* c, int m, bool ddebug = false, bool hdebug = false, bool well_defined = false);

#endif