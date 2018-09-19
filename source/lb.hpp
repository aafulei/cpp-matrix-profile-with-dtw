#ifndef LB_HPP
#define LB_HPP

#include <algorithm>
#include <cmath>
#include <iostream>

#include "pnt.hpp"

#define BAD_LOWER_BOUND -1

double lower_bound_keogh(double const* q, double const* c, int m, int r, bool debug = false);
double lower_bound_front(double const* q, double const* c, int r, bool debug = false);
double lower_bound_back(double const* q, double const* c, int m, int r, bool debug = false);
double lower_bound_update(double lbk, double const* q, double const* c, int m, int r, bool debug = false);

#endif