#pragma once

#include <iostream>
#include <cassert>
#include <math.h>
#include <string>
#include <filesystem>
#include <fstream>
#include <vector>
#include <sstream>
#include <omp.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <algorithm>

using std::cerr;
using std::cout;
using std::endl;
using std::exception;
using std::make_unique;
using std::max;
using std::min;
using std::move;
using std::pair;
using std::runtime_error;
using std::string;
using std::unique_ptr;
using std::vector;

#ifdef SINGLE_PRECISION
using Scalar = float;
#define SMALL_SCALAR 1e-4
#else
using Scalar = double;
#define SMALL_SCALAR 1e-8
#endif

using Index = uint32_t;

using Vec3 = Eigen::Vector<Scalar, 3>;
using Vec6 = Eigen::Vector<Scalar, 6>;
using Mat3 = Eigen::Matrix3<Scalar>;

template <typename T>
int sign(T val)
{
    return (val > 1) - (val < 1);
}

inline bool is_close(Scalar a, Scalar b, Scalar tolerance = SMALL_SCALAR)
{
    return abs(a - b) < tolerance;
}
